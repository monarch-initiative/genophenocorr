import logging
import typing

import hpotk

from phenopackets import GenomicInterpretation, Phenopacket, Disease as PPDisease

from genophenocorr.model import Patient, SampleLabels, Disease
from genophenocorr.model import VariantCoordinates, Variant, Genotype, Genotypes
from genophenocorr.model.genome import GenomeBuild, GenomicRegion, Strand
from ._api import VariantCoordinateFinder, FunctionalAnnotator
from ._audit import Notepad
from ._patient import PatientCreator
from ._phenotype import PhenotypeCreator


class PhenopacketVariantCoordinateFinder(VariantCoordinateFinder[GenomicInterpretation]):
    """
    `PhenopacketVariantCoordinateFinder` figures out :class:`genophenocorr.model.VariantCoordinates`
    and :class:`genophenocorr.model.Genotype` from `GenomicInterpretation` element of Phenopacket Schema.

    :param build: genome build to use in `VariantCoordinates
    :param hgvs_coordinate_finder: the coordinate finder to use for parsing HGVS expressions
    """

    def __init__(self, build: GenomeBuild,
                 hgvs_coordinate_finder: VariantCoordinateFinder[str]):
        self._logger = logging.getLogger(__name__)
        self._build = hpotk.util.validate_instance(build, GenomeBuild, 'build')
        self._hgvs_finder = hpotk.util.validate_instance(hgvs_coordinate_finder, VariantCoordinateFinder,
                                                         'hgvs_coordinate_finder')

    def find_coordinates(self, item: GenomicInterpretation) -> typing.Tuple[VariantCoordinates, Genotype]:
        """Creates a VariantCoordinates object from the data in a given Phenopacket

        Args:
            item (GenomicInterpretation): a genomic interpretation element from Phenopacket Schema
        Returns:
            tuple[VariantCall, Genotype]: A tuple of :class:`VariantCoordinates` and :class:`Genotype`
        """
        if not isinstance(item, GenomicInterpretation):
            raise ValueError(f"item must be a Phenopacket GenomicInterpretation but was type {type(item)}")

        variant_descriptor = item.variant_interpretation.variation_descriptor

        vc = None
        if self._vcf_is_available(variant_descriptor.vcf_record):
            # We have a VCF record.
            if not self._check_assembly(variant_descriptor.vcf_record.genome_assembly):
                raise ValueError(f"Variant id {variant_descriptor.id} for patient {item.subject_or_biosample_id} has a different Genome Assembly than what was given. " \
                                 + f"{variant_descriptor.vcf_record.genome_assembly} is not {self._build.identifier}.")
            contig = self._build.contig_by_name(variant_descriptor.vcf_record.chrom)
            start = int(variant_descriptor.vcf_record.pos) - 1
            ref = variant_descriptor.vcf_record.ref
            alt = variant_descriptor.vcf_record.alt
            end = start + len(ref)
            change_length = end - start

            region = GenomicRegion(contig, start, end, Strand.POSITIVE)
            vc = VariantCoordinates(region, ref, alt, change_length)
        elif self._cnv_is_available(variant_descriptor.variation):
            # We have a CNV.
            variation = variant_descriptor.variation
            seq_location = variation.copy_number.allele.sequence_location
            refseq_contig_name = seq_location.sequence_id.split(':')[1]
            contig = self._build.contig_by_name(refseq_contig_name)

            # Assuming SV coordinates are 1-based (VCF style),
            # so we subtract 1 to transform to 0-based coordinate system
            start = int(seq_location.sequence_interval.start_number.value) - 1
            end = int(seq_location.sequence_interval.end_number.value)
            ref = 'N'
            number = int(variation.copy_number.number.value)
            if number > 2:
                alt = '<DUP>'
            else:
                alt = '<DEL>'
            change_length = end - start

            region = GenomicRegion(contig, start, end, Strand.POSITIVE)
            vc = VariantCoordinates(region, ref, alt, change_length)
        elif len(variant_descriptor.expressions) > 0:
            # We have some expressions. Let's try to find the 1st expression with `hgvs.c` syntax.
            for expression in variant_descriptor.expressions:
                if expression.syntax == 'hgvs.c':
                    vc = self._hgvs_finder.find_coordinates(expression.value)
                    break
        else:
            raise ValueError("Unable to find variant coordinates.")

        # Last, parse genotype.
        genotype = variant_descriptor.allelic_state.label
        gt = self._map_geno_genotype_label(genotype)

        return vc, gt

    def _check_assembly(self, genome_assembly: str) -> bool:
        if '38' in genome_assembly and self._build.identifier == 'GRCh38.p13':
            return True
        elif ('37' in genome_assembly or '19' in genome_assembly) and self._build.identifier == 'GRCh37.p13':
            return True
        else:
            return False


    @staticmethod
    def _vcf_is_available(vcf_record) -> bool:
        """
        Check if we can parse data out of VCF record.
        """
        return (vcf_record.genome_assembly != ''
                and vcf_record.chrom != ''
                and vcf_record.pos >= 0
                and vcf_record.ref != ''
                and vcf_record.alt != '')

    @staticmethod
    def _cnv_is_available(variation):
        seq_location = variation.copy_number.allele.sequence_location
        return (seq_location.sequence_id != ''
                and seq_location.sequence_interval.start_number.value >= 0
                and seq_location.sequence_interval.end_number.value >= 0
                and variation.copy_number.number.value != '')

    @staticmethod
    def _map_geno_genotype_label(genotype: str) -> Genotype:
        """
        Mapping from labels of the relevant GENO terms that is valid as of Oct 2nd, 2023.
        """
        if genotype in ('heterozygous', 'compound heterozygous', 'simple heterozygous'):
            return Genotype.HETEROZYGOUS
        elif genotype == 'homozygous':
            return Genotype.HOMOZYGOUS_ALTERNATE
        elif genotype in ('hemizygous', 'hemizygous X-linked', 'hemizygous Y-linked'):
            return Genotype.HEMIZYGOUS
        else:
            raise ValueError(f'Unknown genotype {genotype}')


class PhenopacketPatientCreator(PatientCreator[Phenopacket]):
    """
    `PhenopacketPatientCreator` transforms `Phenopacket` into :class:`genophenocorr.model.Patient`.
    """

    def __init__(self, build: GenomeBuild,
                 phenotype_creator: PhenotypeCreator,
                 var_func_ann: FunctionalAnnotator,
                 hgvs_coordinate_finder: VariantCoordinateFinder[str]):
        self._logger = logging.getLogger(__name__)
        # Violates DI, but it is specific to this class, so I'll leave it "as is".
        self._coord_finder = PhenopacketVariantCoordinateFinder(build, hgvs_coordinate_finder)
        self._phenotype_creator = hpotk.util.validate_instance(phenotype_creator, PhenotypeCreator, 'phenotype_creator')
        self._func_ann = hpotk.util.validate_instance(var_func_ann, FunctionalAnnotator, 'var_func_ann')

    def process(self, inputs: Phenopacket, notepad: Notepad) -> Patient:
        """Creates a Patient from the data in a given Phenopacket

        Args:
            inputs (Phenopacket): A Phenopacket object
            notepad (Notepad): notepad to write down the issues
        Returns:
            Patient: A Patient object
        """
        sample_id = SampleLabels(label=inputs.subject.id, meta_label=inputs.id if len(inputs.id) > 0 else None)

        # Check phenotypes 
        pfs = notepad.add_subsection('phenotype-features')
        phenotypes = self._phenotype_creator.process(
            ((pf.type.id, not pf.excluded) for pf in inputs.phenotypic_features),
            pfs
        )
        
        # Check diseases
        diseases = self._add_diseases([dis for dis in inputs.diseases], pfs)

        # Check variants
        vs = notepad.add_subsection('variants')
        variants = self._add_variants(sample_id, inputs, vs)
    
        return Patient(sample_id, phenotypes=phenotypes, variants=variants, diseases=diseases)

    def _add_diseases(self, diseases: typing.Sequence[PPDisease], notepad: Notepad) -> typing.Sequence[Disease]:
        """Creates a list of Disease objects from the data in a given Phenopacket

        Args:
            diseases (Sequence[Disease]): A sequence of Phenopacket Disease objects
        Returns:
            Sequence[Dis]: A list of Disease objects
        """
        if len(diseases) == 0:
            notepad.add_warning(f"No diseases found.")
            return []
        final_diseases = []
        for i, dis in enumerate(diseases):
            if not dis.HasField("term"):
                raise ValueError('Could not find term in Disease.')
            term_id = hpotk.TermId.from_curie(dis.term.id)
            # Do not include excluded diseases if we decide to assume excluded if not included 
            final_diseases.append(Disease(term_id, dis.term.label, not dis.excluded))
        return final_diseases
        

    def _add_variants(self, sample_id: SampleLabels, pp: Phenopacket, notepad: Notepad) -> typing.Sequence[Variant]:
        """Creates a list of Variant objects from the data in a given Phenopacket

        Args:
            pp (Phenopacket): A Phenopacket object
        Returns:
            Sequence[Variant]: A list of Variant objects
        """
        variants = []
        for i, interp in enumerate(pp.interpretations):
            for genomic_interp in interp.diagnosis.genomic_interpretations:
                try:
                    vc, gt = self._coord_finder.find_coordinates(genomic_interp)
                except ValueError:
                    notepad.add_warning(f"Expected a VCF record, a VRS CNV, or an expression with `hgvs.c` but had an error retrieving any from patient {sample_id}",
                                        "Remove variant from testing")
                    continue

                try:
                    tx_annotations = self._func_ann.annotate(vc)
                except ValueError:
                    notepad.add_warning(f"Patient {pp.id} has an error with variant {vc.variant_key}",
                                        "Try again or remove variant form testing")
                    continue

                if len(tx_annotations) == 0:
                    notepad.add_warning(f"Patient {pp.id} has an error with variant {vc.variant_key}",
                                        "Remove variant from testing")
                    continue

                genotype = Genotypes.single(sample_id, gt)
                variants.append(Variant(vc, tx_annotations, genotype))

        if len(variants) == 0:
            notepad.add_warning(f"Patient {pp.id} has no variants to work with")

        return variants

