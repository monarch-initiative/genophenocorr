import hpotk
import pytest


from genophenocorr.preprocessing import PhenopacketPatientCreator, PhenotypeCreator, CohortCreator
from genophenocorr.preprocessing import configure_patient_creator, configure_cohort_creator, load_phenopacket
from genophenocorr.preprocessing import load_phenopacket_folder
from genophenocorr.preprocessing import Level


class TestPhenopacketPatientCreator:

    @pytest.fixture
    def phenopacket_patient_creator(self, toy_hpo: hpotk.MinimalOntology) -> PhenopacketPatientCreator:
        return configure_patient_creator(toy_hpo)

    @pytest.mark.skip('Skipping online test')
    def test_load_phenopacket(self, phenopacket_patient_creator: PhenopacketPatientCreator):
        pp = load_phenopacket('../docs/data/simple_cohort/PMID_36446582_KBG12.json')
        patient = phenopacket_patient_creator.create_patient(pp)
        print(patient)


class TestPhenopacketCohortCreator:

    @pytest.fixture
    def phenopacket_cohort_creator(self, toy_hpo: hpotk.MinimalOntology) -> CohortCreator:
        return configure_cohort_creator(toy_hpo)

    @pytest.mark.skip('Skipping online test')
    def test_load_phenopacket(self, phenopacket_cohort_creator: CohortCreator):
        cohort = load_phenopacket_folder('docs/data/simple_cohort', phenopacket_cohort_creator)
        print(cohort)


class TestPhenotypeCreator:

    @pytest.fixture
    def phenotype_creator(self, toy_hpo: hpotk.MinimalOntology,
                          toy_validation_runner: hpotk.validate.ValidationRunner) -> PhenotypeCreator:
        return PhenotypeCreator(toy_hpo, toy_validation_runner)

    @pytest.mark.parametrize('curie, message, solution',
                             [
                                 (
                                         'WalterWhite',
                                         '#0 The CURIE WalterWhite has no colon `:` or underscore `_`',
                                         'Ensure the term ID consists of a prefix (e.g. `HP`) and id (e.g. `0001250`) joined by colon `:` or underscore `_`'
                                 ),
                                 (
                                         'NCIT_C12345',
                                         '#0 NCIT:C12345 is not an HPO term',
                                         'Remove non-HPO concepts from the analysis input'
                                 ),
                                 (
                                         'HP:999999',
                                         '#0 HP:999999 is not in HPO version `2022-10-05`',
                                         'Correct the HPO term or use the latest HPO for the analysis'
                                 ),
                             ])
    def test_input_curies(self, phenotype_creator: PhenotypeCreator,
                          curie: str, message: str, solution: str):
        inputs = ((curie, True),)

        notepad = phenotype_creator.prepare_notepad('top-level')
        phenotypes = phenotype_creator.process(inputs, notepad)

        assert len(phenotypes) == 0

        issues = list(notepad.issues)
        assert len(issues) == 2

        first = issues[0]
        assert first.message == message
        assert first.solution == solution


    def test_annotation_propagation(self, phenotype_creator: PhenotypeCreator):
        inputs = (
            ('HP:0001250', False), # Seizure
            ('HP:0002266', True),  # Focal clonic seizure
        )

        notepad = phenotype_creator.prepare_notepad('top-level')
        phenotypes = phenotype_creator.process(inputs, notepad)

        assert len(phenotypes) == 2

        issues = list(notepad.issues)
        assert len(issues) == 1

        first = issues[0]
        assert first.level == Level.ERROR
        assert first.message == 'Terms should not contain both present Focal clonic seizure [HP:0002266] and its present or excluded ancestor Seizure [HP:0001250]'
        assert first.solution is None
