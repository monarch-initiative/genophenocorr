import abc
import typing
from collections import namedtuple, defaultdict

import hpotk
import pandas as pd

from genophenocorr.model import VariantEffect, Patient, FeatureType
from .predicate import PatientCategory

PatientsByHPO = namedtuple('PatientsByHPO', field_names=['all_with_hpo', 'all_without_hpo'])


class GenotypePhenotypeAnalysisResult:
    """
    `GenotypePhenotypeAnalysisResult` summarizes results of genotype-phenotype correlation analysis of a cohort.
    """

    def __init__(self, n_usable: pd.Series,
                 all_counts: pd.DataFrame,
                 pvals: pd.Series,
                 corrected_pvals: typing.Optional[pd.Series],
                 phenotype_categories: typing.Iterable[PatientCategory],
                 question: str):
        self._n_usable = n_usable
        self._all_counts = all_counts
        self._pvals = pvals
        self._corrected_pvals = corrected_pvals
        self._phenotype_categories = tuple(phenotype_categories)
        self._question = question

    @property
    def n_usable(self) -> pd.Series:
        """
        Get a :class:`pandas.Series` with mapping from HPO term and number of patients where the term was either
         present or excluded, and are, thus, usable for genotype-phenotype correlation analysis.
        """
        return self._n_usable

    @property
    def all_counts(self) -> pd.DataFrame:
        """
        Get a :class:`pandas.DataFrame` with counts of patients in genotype and phenotype groups.

        An example for a genotype predicate that bins into two categories (`Yes` and `No`) based on presence
        of a missense variant in transcript, and phenotype predicate that checks presence/absence of a phenotype term::

                                   Has MISSENSE_VARIANT in NM_123456.7
                                   No       Yes
            HPO         Present
            HP:0001166  Yes        1        13
            HP:0001166  No         7        5
            HP:0003011  Yes        3        11
            HP:0003011  No         7        16
            HP:0012638  Yes        9        15
            HP:0012638  No         3        4

        """
        return self._all_counts

    @property
    def pvals(self) -> pd.Series:
        """
        Get a :class:`pandas.Series` with p values for each tested HPO term.
        """
        return self._pvals

    @property
    def corrected_pvals(self) -> typing.Optional[pd.Series]:
        """
        Get an optional :class:`pandas.Series` with p values for each tested HPO term after multiple testing correction.
        """
        return self._corrected_pvals

    @property
    def phenotype_categories(self) -> typing.Sequence[PatientCategory]:
        """
        Get a sequence of phenotype patient categories that can be investigated.
        """
        return self._phenotype_categories

    def summarize(self, hpo: hpotk.MinimalOntology,
                  phenotype_category: PatientCategory) -> pd.DataFrame:
        """
        Create a data frame with summary of the genotype phenotype analysis.

        The *rows* of the frame correspond to the analyzed HPO terms.

        The columns of the data frame have `Count` and `Percentage` per used genotype predicate.
        For instance, if we used :class:`genophenocorr.analysis.predicate.genotype.VariantEffectPredicate`
        which compares phenotype with and without variants with certain :class:`genophenocorr.model.VariantEffect`,
        we will have the following columns:

        =====  =======  =====  =======
        Yes             No
        --------------  --------------
        Count  Percent  Count  Percent
        =====  =======  =====  =======

        The final data frame may look something like this::

            MISSENSE_VARIANT on `NM_1234.5`                       No                  Yes
                                                                  Count    Percent    Count    Percent   p value   Corrected p value
            Arachnodactyly [HP:0001166]                           1         3.8461    13       50.0000   0.000781  0.020299
            Abnormality of the musculature [HP:0003011]           6        35.2941    11       64.7058   1.000000  1.000000
            Abnormal nervous system physiology [HP:0012638]       9        37.5000    15       62.5000   1.000000  1.000000
            ...                                                   ...      ...        ...      ...       ...       ...
        """
        if phenotype_category not in self._phenotype_categories:
            raise ValueError(f'Unknown phenotype category: {phenotype_category}. Use one of {self._phenotype_categories}')

        # Row index: a list of tested HPO terms
        pheno_idx = pd.Index(self._n_usable.index)
        # Column index: multiindex of counts and percentages for all genotype predicate groups
        geno_idx = pd.MultiIndex.from_product(
            iterables=(self._all_counts.columns, ('Count', 'Percent')),
            names=(self._question, None))

        # We'll fill this frame with data
        df = pd.DataFrame(index=pheno_idx, columns=geno_idx)

        # We must choose the counts with the `phenotype_category` of interest.
        # This depends on the used phenotype predicate.
        # Usually, we are interested patients who HAVE a feature,
        # but we can also investigate those who do NOT have the feature.
        # Hence, we need `phenotype_category`.
        counts = self._all_counts.loc[(slice(None), phenotype_category), :]
        counts.set_index(pheno_idx, inplace=True)

        # Fill the frame cells
        for col in geno_idx.levels[0]:
            cnt = counts[col]
            # Format `Count` as `N/M` string, where `N` is the sample count
            # and `M` is the number of usable samples - the samples where
            # both genotype and phenotype predicates were able to make a call.
            df[col, 'Count'] = cnt.map(str) + '/' + self._n_usable.map(str)
            df[col, 'Percent'] = cnt * 100 / self._n_usable

        # Add columns with p values and corrected p values (if present)
        df.insert(df.shape[1], ('', self._pvals.name), self._pvals)
        if self._corrected_pvals is not None:
            df.insert(df.shape[1], ('', self._corrected_pvals.name), self._corrected_pvals)

        # Format the index values: `HP:0001250` -> `Seizure [HP:0001250]`
        labeled_idx = df.index.map(lambda term_id: f'{hpo.get_term(term_id).name} [{term_id.value}]')

        # Last, sort by corrected p value or just p value
        df = df.set_index(labeled_idx)
        if self._corrected_pvals is not None:
            return df.sort_values(by=[('', self._corrected_pvals.name), ('', self._pvals.name)])
        else:
            return df.sort_values(by=('', self._pvals.name))


class CohortAnalysis(metaclass=abc.ABCMeta):
    """
    `CohortAnalysis` is a driver class for running genotype-phenotype correlation analyses.

    The class provides various methods to test genotype-phenotype correlations. All methods wrap results
    into :class:`GenotypePhenotypeAnalysisResult`.
    """

    @abc.abstractmethod
    def compare_by_variant_effect(self, effect: VariantEffect, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        """
        Compares variants with `effect` vs. variants with all other variant effects.

        :param effect: variant effect
        :param tx_id: the accession of the transcript of interest
        """
        pass

    @abc.abstractmethod
    def compare_by_variant_key(self, variant_key: str) -> GenotypePhenotypeAnalysisResult:
        """
        Compares variant with `variant_key` vs all the other variants.

        .. seealso::

          :attr:`genophenocorr.model.VariantCoordinates.variant_key`

        :param variant_key: a `str` with the variant key, e.g. ``X_12345_12345_C_G``
        """
        pass

    @abc.abstractmethod
    def compare_by_exon(self, exon_number: int, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        """
        Compares variants with `effect` vs. variants with all other variant effects.

        .. note::

          We use 1-based numbering to number the exons, not the usual 0-based numbering of the computer science.
          Therefore, the first exon of the transcript has ``exon_number==1``, the second exon is ``2``, and so on ...

        .. note::

          We do not check if the `exon_number` spans beyond the number of exons of the given `transcript_id`!
          Therefore, ``exon_number==10,000`` will effectively return :attr:`BooleanPredicate.FALSE` for *all* patients!!! 😱
          Well, at least the patients of the *Homo sapiens sapiens* taxon...

        :param exon_number: a positive `int` with exon number
        :param tx_id: the accession of the transcript of interest
        """
        pass

    @abc.abstractmethod
    def compare_by_variant_keys(self, a: str, b: str) -> GenotypePhenotypeAnalysisResult:
        """
        Compare genotype-phenotype correlation between groups with variant `a` vs. variant `b`.

        .. seealso::

          :attr:`genophenocorr.model.VariantCoordinates.variant_key`

        :param a: a `str` with variant key of variant `a`, e.g. ``X_12345_12345_C_G``
        :param b: a `str` with variant key of variant `b`
        """
        pass

    @abc.abstractmethod
    def compare_by_protein_feature_type(self, feature_type: FeatureType, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        """
        Compare genotype-phenotype correlation between variants that affect a given `feature_type` and the variants
        affecting the rest of the protein.

        Args:
            feature_type (FeatureType): the protein feature type of interest
            tx_id: the accession of the transcript of interest
        """
        pass

    @abc.abstractmethod
    def compare_by_protein_feature(self, feature: str, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        """
        Compare genotype-phenotype correlation between variants that affect a given `feature` and the variants
        affecting the rest of the protein.

        .. seealso::

          The protein feature names can be accessed at :attr:`genophenocorr.model.ProteinFeature.name`.

        Args:
            feature (string): feature identifier, e.g. ``DNA-binding``.
            tx_id (string): the accession of the transcript of interest
        """
        pass

    @staticmethod
    def _check_min_perc_patients_w_hpo(min_perc_patients_w_hpo: typing.Union[int, float],
                                       cohort_size: int) -> float:
        """
        Check if the input meets the requirements.
        """
        if isinstance(min_perc_patients_w_hpo, int):
            if min_perc_patients_w_hpo > 0:
                return min_perc_patients_w_hpo / cohort_size
            else:
                raise ValueError(f'`min_perc_patients_w_hpo` must be a positive `int` '
                                 f'but got {min_perc_patients_w_hpo}')
        elif isinstance(min_perc_patients_w_hpo, float):
            if 0 < min_perc_patients_w_hpo <= 1:
                return min_perc_patients_w_hpo
            else:
                raise ValueError(f'`min_perc_patients_w_hpo` must be a `float` in range (0, 1] '
                                 f'but got {min_perc_patients_w_hpo}')
        else:
            raise ValueError(f'`min_perc_patients_w_hpo` must be a positive `int` or a `float` in range (0, 1] '
                             f'but got {type(min_perc_patients_w_hpo)}')

    @staticmethod
    def _group_patients_by_hpo(phenotypic_features: typing.Iterable[hpotk.TermId],
                               patients: typing.Iterable[Patient],
                               hpo: hpotk.GraphAware,
                               missing_implies_excluded: bool) -> PatientsByHPO:
        all_with_hpo = defaultdict(list)
        all_without_hpo = defaultdict(list)
        for hpo_term in phenotypic_features:
            for patient in patients:
                found = False
                for pf in patient.present_phenotypes():
                    if hpo_term == pf.identifier or hpo.graph.is_ancestor_of(hpo_term, pf):
                        # Patient is annotated with `hpo_term` because `pf` is equal to `hpo_term`
                        # or it is a descendant of `hpo_term`.
                        all_with_hpo[hpo_term].append(patient)

                        # If one `pf` of the patient is found to be a descendant of `hpo`, we must break to prevent
                        # adding the patient to `present_hpo` more than once due to another descendant!
                        found = True
                        break
                if not found:
                    # The patient is not annotated by the `hpo_term`.

                    if missing_implies_excluded:
                        # The `hpo_term` annotation is missing, hence implicitly excluded.
                        all_without_hpo[hpo_term].append(patient)
                    else:
                        # The `hpo_term` must be explicitly excluded patient to be accounted for.
                        for ef in patient.excluded_phenotypes():
                            if hpo_term == ef.identifier or hpo.graph.is_descendant_of(hpo_term, ef):
                                all_with_hpo[hpo_term].append(patient)
                                break

        return PatientsByHPO(all_with_hpo, all_without_hpo)
