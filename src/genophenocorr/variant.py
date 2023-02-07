from curses.ascii import isdigit
import varcode as vc
import pyensembl
import re

class Variant:
    def __init__(self,ref, phenopack):
        self._phenopack = phenopack
        self._variant = self.__find_variant(reference = ref)
        
    def __find_variant(self, reference):
        if len(self._phenopack.interpretations) != 0:
            Interp = self._phenopack.interpretations[0]
            try:
                genInterp = Interp.diagnosis.genomic_interpretations[0]
                varInterp = genInterp.variant_interpretation.variation_descriptor.vcf_record
                contig = re.sub(r'[^0-9MXY]', '', varInterp.chrom)
                if len(contig) == 0 or (contig.isdigit() and (int(contig) == 0 or int(contig) >= 24)):
                    ## Chromosome can only be values 1-23, X, Y, or M
                    raise ValueError(f"Contig did not work: {varInterp.chrom}")
                start = varInterp.pos
                ref = varInterp.ref
                alt = varInterp.alt
                if reference.lower() == 'hg37' or reference.lower() == 'grch37' or reference.lower() == 'hg19':
                    ens = pyensembl.ensembl_grch37
                elif reference.lower() == 'hg38' or reference.lower() == 'grch38':
                    ens = pyensembl.ensembl_grch38
                else:
                    raise ValueError(f'Unknown Reference {reference}. Please use hg37 or hg38.')
                myVar = vc.Variant(contig, start, ref, alt, ensembl = ens)
                return myVar
            except AttributeError:
                raise AttributeError('Could not find Variants')
        else:
            raise ValueError('No Genomic Interpretations. Cannot find Variant.')

    @property
    def variant(self):
        return self._variant

    @property
    def variant_string(self):
        return self.variant.short_description

    @property
    def variant_type(self):
        all_types = []
        if not self.top_effect.short_description.endswith("*") and self.variant.is_snv:
            all_types.append('missense')
        if self.top_effect.short_description.endswith("*") and self.variant.is_snv:
            all_types.append('nonsense')
        if 'dup' in self.top_effect.short_description:
            all_types.append('duplication')
        if self.variant.is_deletion:
            all_types.append('deletion')
        if self.variant.is_insertion:
            all_types.append('insertion')
        if self.variant.is_transition:
            all_types.append('transition')
        if self.variant.is_transversion:
            all_types.append('transversion')
        if self.variant.is_indel:
            all_types.append('indel')
        return all_types


    @property
    def genomic_location(self):
        if self._variant is not None:
            return self._variant.start
        else:
            return None

    @property
    def effects(self):
        if self._variant is not None:
            return self.variant.effects()
        else:
            return None

    @property
    def top_effect(self):
        if self.effects is not None:
            return self.effects.top_priority_effect()
        else:
            return None

    @property
    def all_transcripts(self):
        if self._variant is not None:
            return self.variant.transcripts
        else:
            return None

    @property
    def all_transcript_ids(self):
        if self.variant.transcripts is not None:
            return self.variant.transcript_ids
        else:
            return None

    @property
    def top_effect_transcript(self):
        if self.top_effect is not None:
            return self.top_effect.transcript
        else:
            return None

    @property
    def top_effected_protein(self):
        if self.top_effect_transcript is not None:
            return self.top_effect_transcript.protein_id
        else:
            return None

    @property
    def protein_effect_location(self):
        # Currently only works with single amino acid substitutions
        loc = None
        if self.top_effected_protein is not None:
            pattern = re.compile(r'p\.[A-Z](\d+)[A-Z]')
            if pattern.match(self.top_effect.short_description):
                loc = int(re.sub(pattern, '\\g<1>', self.top_effect.short_description))
        return loc
            