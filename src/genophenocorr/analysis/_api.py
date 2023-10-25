import abc
import typing
from collections import namedtuple, defaultdict

import hpotk
import pandas as pd

from genophenocorr.model import VariantEffect, Patient
from .predicate import PatientCategory

PatientsByHPO = namedtuple('PatientsByHPO', field_names=['all_with_hpo', 'all_without_hpo'])


class GenotypePhenotypeAnalysisResult:

    def __init__(self, n_usable: pd.Series,
                 all_counts: pd.DataFrame,
                 pvals: pd.Series,
                 corrected_pvals: typing.Optional[pd.Series],
                 question: str):
        self._n_usable = n_usable
        self._all_counts = all_counts
        self._pvals = pvals
        self._corrected_pvals = corrected_pvals
        self._question = question

    @property
    def n_usable(self) -> pd.Series:
        return self._n_usable

    @property
    def all_counts(self) -> pd.DataFrame:
        return self._all_counts

    @property
    def pvals(self) -> pd.Series:
        return self._pvals

    @property
    def corrected_pvals(self) -> typing.Optional[pd.Series]:
        return self._corrected_pvals

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
            df[col, 'Count'] = cnt
            df[col, 'Percent'] = cnt * 100 / self._n_usable

        # Add columns with p values and corrected p values
        df.insert(df.shape[1], ('', self._pvals.name), self._pvals)
        df.insert(df.shape[1], ('', self._corrected_pvals.name), self._corrected_pvals)

        # Format the index values: `HP:0001250` -> `Seizure [HP:0001250]`
        labeled_idx = df.index.map(lambda term_id: f'{hpo.get_term(term_id).name} [{term_id.value}]')

        # Last, sort by corrected p value
        return df.set_index(labeled_idx).sort_values(by=[('', self._corrected_pvals.name)])


class AbstractCohortAnalysis(metaclass=abc.ABCMeta):
    # TODO
    #  - factor out convenient methods from `CohortAnalysis`
    #  - finish splitting predicates to `genotype` and `phenotype`

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
