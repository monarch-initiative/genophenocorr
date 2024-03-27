import os
import sys
import typing
import warnings

import hpotk
# pyright: reportGeneralTypeIssues=false
from google.protobuf.json_format import Parse
from phenopackets import Phenopacket
from tqdm import tqdm

from genophenocorr.model import Cohort
from genophenocorr.model.genome import GRCh37, GRCh38, GenomeBuild
from ._api import FunctionalAnnotator
from ._audit import NotepadTree
from ._patient import CohortCreator
from ._phenopacket import PhenopacketPatientCreator
from ._phenotype import PhenotypeCreator
from ._variant import VarCachingFunctionalAnnotator, VariantAnnotationCache
from ._vep import VepFunctionalAnnotator
from ._vv import VVHgvsVariantCoordinateFinder

VALIDATION_POLICIES = {'none', 'lenient', 'strict'}


def configure_caching_cohort_creator(
        hpo: hpotk.MinimalOntology,
        genome_build: str = 'GRCh38.p13',
        validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
        cache_dir: typing.Optional[str] = None,
        variant_fallback: str = 'VEP',
        timeout: int = 10,
        ) -> CohortCreator[Phenopacket]:
    """
    A convenience function for configuring a caching :class:`genophenocorr.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO. Other options are optional.

    :param hpo: a HPO instance.
    :param genome_build: name of the genome build to use, choose from `{'GRCh37.p13', 'GRCh38.p13'}`.
    :param validation_runner: an instance of the validation runner.
    :param cache_dir: path to the folder where we will cache the results fetched from the remote APIs or `None`
     if the data should be cached in `.cache` folder in the current working directory.
     In any case, the directory will be created if it does not exist (including non-existing parents).
    :param variant_fallback: the fallback variant annotator to use if we cannot find the annotation locally.
     Choose from ``{'VEP'}`` (just one fallback implementation is available at the moment).
    :param timeout: timeout in seconds for the VEP API
    """
    if cache_dir is None:
        cache_dir = os.path.join(os.getcwd(), '.genophenocorr_cache')
    os.makedirs(cache_dir, exist_ok=True)

    build = _configure_build(genome_build)
    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    functional_annotator = _configure_functional_annotator(cache_dir, variant_fallback, timeout)
    hgvs_annotator = VVHgvsVariantCoordinateFinder(build)
    pc = PhenopacketPatientCreator(build, phenotype_creator, functional_annotator, hgvs_annotator)

    return CohortCreator(pc)


def configure_caching_patient_creator(
        hpo: hpotk.MinimalOntology,
        genome_build: str = 'GRCh38.p13',
        validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
        cache_dir: typing.Optional[str] = None,
        variant_fallback: str = 'VEP',
        timeout: int = 10,
    ) -> PhenopacketPatientCreator:
    """
    A convenience function for configuring a caching :class:`genophenocorr.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO. Other options are optional.

    :param hpo: a HPO instance.
    :param genome_build: name of the genome build to use, choose from `{'GRCh37.p13', 'GRCh38.p13'}`.
    :param validation_runner: an instance of the validation runner.
    :param cache_dir: path to the folder where we will cache the results fetched from the remote APIs or `None`
     if the data should be cached in `.cache` folder in the current working directory.
     In any case, the directory will be created if it does not exist (including non-existing parents).
    :param variant_fallback: the fallback variant annotator to use if we cannot find the annotation locally.
     Choose from ``{'VEP'}`` (just one fallback implementation is available at the moment).
    :param timeout: timeout in seconds for the VEP API
    """
    warnings.warn('`configure_caching_patient_creator` was deprecated. '
                  'Use `configure_caching_cohort_creator` instead', DeprecationWarning, stacklevel=2)
    if cache_dir is None:
        cache_dir = os.path.join(os.getcwd(), '.genophenocorr_cache')

    build = _configure_build(genome_build)
    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    functional_annotator = _configure_functional_annotator(cache_dir, variant_fallback, timeout)
    hgvs_annotator = VVHgvsVariantCoordinateFinder(build)
    return PhenopacketPatientCreator(build, phenotype_creator, functional_annotator, hgvs_annotator)


def configure_cohort_creator(
        hpo: hpotk.MinimalOntology,
        genome_build: str = 'GRCh38.p13',
        validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
        variant_fallback: str = 'VEP',
        timeout: int = 10,
        ) -> CohortCreator[Phenopacket]:
    """
    A convenience function for configuring a non-caching :class:`genophenocorr.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO. Other options are optional

    :param hpo: a HPO instance.
    :param genome_build: name of the genome build to use, choose from `{'GRCh37.p13', 'GRCh38.p13'}`.
    :param validation_runner: an instance of the validation runner.
     if the data should be cached in `.cache` folder in the current working directory.
     In any case, the directory will be created if it does not exist (including non-existing parents).
    :param variant_fallback: the fallback variant annotator to use if we cannot find the annotation locally.
     Choose from ``{'VEP'}`` (just one fallback implementation is available at the moment).
    :param timeout: timeout in seconds for the VEP API
    """
    build = _configure_build(genome_build)

    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    functional_annotator = _configure_fallback_functional(variant_fallback, timeout)
    hgvs_annotator = VVHgvsVariantCoordinateFinder(build)
    pc = PhenopacketPatientCreator(build, phenotype_creator, functional_annotator, hgvs_annotator)

    return CohortCreator(pc)


def configure_patient_creator(hpo: hpotk.MinimalOntology,
                              genome_build: str = 'GRCh38.p13',
                              validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
                              variant_fallback: str = 'VEP',
                              validation: str = 'lenient') -> PhenopacketPatientCreator: #Rename to something more understandable by user
    """
                                ^^^ none, lenient, strict -
                                none = run unless unrunnable
                                lenient = fix what we can, abort unfixable
                                strict = abort at any issue
    A convenience function for configuring a non-caching :class:`genophenocorr.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO. Other options are optional

    :param hpo: a HPO instance.
    :param genome_build: name of the genome build to use, choose from `{'GRCh37.p13', 'GRCh38.p13'}`.
    :param validation_runner: an instance of the validation runner.
     if the data should be cached in `.cache` folder in the current working directory.
     In any case, the directory will be created if it does not exist (including non-existing parents).
    :param variant_fallback: the fallback variant annotator to use if we cannot find the annotation locally.
     Choose from ``{'VEP'}`` (just one fallback implementation is available at the moment).
    """
    warnings.warn('`configure_patient_creator` was deprecated. '
                  'Use `configure_cohort_creator` instead', DeprecationWarning, stacklevel=2)

    build = _configure_build(genome_build)

    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    functional_annotator = _configure_fallback_functional(variant_fallback)
    hgvs_annotator = VVHgvsVariantCoordinateFinder(build)
    return PhenopacketPatientCreator(build, phenotype_creator, functional_annotator, hgvs_annotator)


def _configure_build(genome_build: str) -> GenomeBuild:
    if genome_build == 'GRCh38.p13':
        return GRCh38
    elif genome_build == 'GRCh37.p13':
        return GRCh37
    else:
        raise ValueError(f'Unknown build {genome_build}. Choose from [\'GRCh37.p13\', \'GRCh38.p13\']')


def _setup_phenotype_creator(hpo: hpotk.MinimalOntology,
                             validator: typing.Optional[hpotk.validate.ValidationRunner]) -> PhenotypeCreator:
    if validator is None:
        # This will be the default validator
        validator = hpotk.validate.ValidationRunner((
            hpotk.validate.ObsoleteTermIdsValidator(hpo),
            hpotk.validate.AnnotationPropagationValidator(hpo),
            hpotk.validate.PhenotypicAbnormalityValidator(hpo)
        ))
    else:
        validator = hpotk.util.validate_instance(validator, hpotk.validate.ValidationRunner, 'validator')
    return PhenotypeCreator(hpo, validator)


def _configure_functional_annotator(
        cache_dir: str,
        variant_fallback: str,
        timeout: int,
        ) -> FunctionalAnnotator:

    # (2) FunctionalAnnotator
    # Setup fallback
    fallback = _configure_fallback_functional(variant_fallback, timeout)

    # Setup variant cache
    var_cache_dir = os.path.join(cache_dir, 'variant_cache')
    os.makedirs(var_cache_dir, exist_ok=True)
    var_cache = VariantAnnotationCache(var_cache_dir)

    # Assemble the final functional annotator
    return VarCachingFunctionalAnnotator(var_cache, fallback)


def _configure_fallback_functional(
        variant_fallback: str,
        timeout: int,
        ) -> FunctionalAnnotator:
    if variant_fallback == 'VEP':
        fallback = VepFunctionalAnnotator(timeout=timeout)
    else:
        raise ValueError(f'Unknown variant fallback annotator type {variant_fallback}')
    return fallback


def load_phenopacket_folder(pp_directory: str,
                            cohort_creator: CohortCreator[Phenopacket],
                            validation_policy: str = 'none') -> Cohort:
    """
    Creates a Patient object for each phenopacket formatted JSON file in the given directory `pp_directory`.

    :param pp_directory: path to a folder with phenopacket JSON files. An error is raised if the path does not point to
      a directory with at least one phenopacket.
    :param cohort_creator: cohort creator for turning a sequence of phenopacket into a :class:`genophenocorr.model.Cohort`.
    :param validation_policy: a `str` with the validation policy. The value must be one of `{'none', 'lenient', 'strict'}`
    :return: a cohort made of the phenopackets
    """
    # Check inputs before moving a finger.
    hpotk.util.validate_instance(cohort_creator, CohortCreator, 'cohort_creator')
    fpath_pp_abs = os.path.abspath(pp_directory)
    if not os.path.isdir(fpath_pp_abs):
        raise ValueError(f'`{fpath_pp_abs}` does not point to a directory')
    if validation_policy.lower() not in VALIDATION_POLICIES:
        raise ValueError(f'{validation_policy} must be one of {VALIDATION_POLICIES}')

    # Load phenopackets
    pps = _load_phenopacket_dir(pp_directory)
    if len(pps) < 1:
        raise ValueError(f"No phenopackets could be parsed from `{fpath_pp_abs}`")

    # Turn phenopackets into a cohort using the cohort creator.
    # Keep track of the progress by wrapping the list of phenopackets
    # with TQDM 😎
    cohort_iter = tqdm(pps, desc='Patients Created')
    notepad = cohort_creator.prepare_notepad(f'{len(pps)} phenopacket(s) found at `{pp_directory}`')
    cohort = cohort_creator.process(cohort_iter, notepad)


    validation_summary = _summarize_validation(validation_policy, notepad)
    print(os.linesep.join(validation_summary), file=sys.stderr)
    if validation_policy == 'none':
        # No validation
        return cohort
    elif validation_policy == 'lenient':
        if notepad.has_errors(include_subsections=True):
            raise ValueError('Cannot load cohort due to errors')
    elif validation_policy == 'strict':
        if notepad.has_errors_or_warnings(include_subsections=True):
            raise ValueError('Cannot load cohort due to warnings/errors')
    else:
        # Bug, please report to the developers.
        raise ValueError(f'Unexpected policy {validation_policy}')

    # create cohort from patients
    return cohort


def _summarize_validation(policy: str,
                          notepad: NotepadTree,
                          indent: int = 2):
    lines = [f'Validated under {policy} policy']
    n_errors = sum(node.error_count() for node in notepad.iterate_nodes())
    n_warnings = sum(node.warning_count() for node in notepad.iterate_nodes())
    if n_errors > 0 or n_warnings > 0:
        lines.append('Showing errors and warnings')
        for node in notepad.iterate_nodes():
            if node.has_errors_or_warnings(include_subsections=True):
                # We must report the node label even if there are no issues with the node.
                l_pad = ' ' * (node.level * indent)
                lines.append(l_pad + node.label)
                if node.has_errors():
                    lines.append(l_pad + ' errors:')
                    for error in node.errors():
                        lines.append(l_pad + ' ' + error.message
                                     + (f'. {error.solution}' if error.solution else ''))
                if node.has_warnings():
                    lines.append(l_pad + ' warnings:')
                    for warning in node.warnings():
                        lines.append(l_pad + ' ·' + warning.message
                                     + (f'. {warning.solution}' if warning.solution else ''))
    else:
        lines.append('No errors or warnings were found')
    return lines


def _load_phenopacket_dir(pp_dir: str) -> typing.Sequence[Phenopacket]:
    patients = []
    for patient_file in os.listdir(pp_dir):
        if patient_file.endswith('.json'):
            phenopacket_path = os.path.join(pp_dir, patient_file)
            pp = load_phenopacket(phenopacket_path)
            patients.append(pp)
    return patients


def load_phenopacket(phenopacket_path: str) -> Phenopacket:
    """
    Load phenopacket JSON file.

    :param phenopacket_path: a `str` pointing to phenopacket JSON file.
    """
    with open(phenopacket_path) as f:
        return Parse(f.read(), Phenopacket())
