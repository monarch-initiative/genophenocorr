import abc
import typing

import pandas as pd

from statsmodels.stats import multitest

from genophenocorr.model import Patient

from .predicate import PolyPredicate, PatientCategory
from .predicate.phenotype import PhenotypePolyPredicate, P

from ._api import GenotypePhenotypeAnalysisResult
from ._stats import run_fisher_exact, run_recessive_fisher_exact


class GPAnalyzer(typing.Generic[P], metaclass=abc.ABCMeta):
    """
    `GPAnalyzer` calculates p values for genotype-phenotype correlation of phenotypic features of interest.
    """

    @abc.abstractmethod
    def analyze(
            self,
            patients: typing.Iterable[Patient],
            pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
            gt_predicate: PolyPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        """
        Test for association between `phenotypic_features` of interest groups or `patients` determined
        by the `predicate`.

        :return: :class:`GenotypePhenotypeAnalysisResult` object with the results.
        """
        pass

    @staticmethod
    def _count_patients(
            patients: typing.Iterable[Patient],
            pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
            gt_predicate: PolyPredicate,
    ) -> typing.Tuple[typing.Iterable[PatientCategory], pd.Series, typing.Mapping[P, pd.DataFrame]]:
        phenotypes = set()
        categories = set()
        for predicate in pheno_predicates:
            categories.update(predicate.get_categories())
            phenotypes.update(predicate.phenotypes)

        n_usable_patients = pd.Series(data=0, index=pd.Index(phenotypes))

        # Apply genotype and phenotype predicates
        counts = {}
        for ph_predicate in pheno_predicates:
            phenotypes = ph_predicate.phenotypes

            for phenotype in phenotypes:
                if phenotype not in counts:
                    # Make an empty frame for keeping track of the counts.
                    counts[phenotype] = pd.DataFrame(
                        data=0,
                        index=pd.Index(
                            data=ph_predicate.get_categories(),
                            name=ph_predicate.get_question(),
                        ),
                        columns=pd.Index(
                            data=gt_predicate.get_categories(),
                            name=gt_predicate.get_question(),
                        ),
                    )

            for patient in patients:
                pheno_cat = ph_predicate.test(patient)
                geno_cat = gt_predicate.test(patient)

                if pheno_cat is not None and geno_cat is not None:
                    counts[pheno_cat.phenotype].loc[pheno_cat.category, geno_cat.category] += 1
                    n_usable_patients[pheno_cat.phenotype] += 1

        return categories, n_usable_patients, counts


class FisherExactAnalyzer(typing.Generic[P], GPAnalyzer[P]):
    """
    `FisherExactAnalyzer` uses Fisher's exact test to calculate p value for phenotypic features of interest.

    Following the test, the code applies one of the multiple testing corrections provided by
    :func:`statsmodels.stats.multitest.multipletests` function at given `mtc_alpha`.

    :param p_val_correction: a `str` with name of the multiple testing correction method or `None` if no correction
      should be applied.
    :param mtc_alpha: a `float` in range :math:`(0, 1]` with the multiple testing correction alpha value.
    """

    def __init__(
            self,
            p_val_correction: typing.Optional[str] = None,
            mtc_alpha: float = .05,
    ):
        self._correction = p_val_correction
        self._mtc_alpha = mtc_alpha

    def analyze(
            self,
            patients: typing.Iterable[Patient],
            pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
            gt_predicate: PolyPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        # 1) Count the patients
        categories, n_usable, all_counts = GPAnalyzer._count_patients(
            patients, pheno_predicates, gt_predicate,
        )

        # 2) Statistical tests
        pheno_idx = pd.Index(n_usable.index, name='p_val')
        pvals = pd.Series(float('nan'), index=pheno_idx, name='p value')
        for phenotype in pheno_idx:
            counts = all_counts[phenotype]
            # TODO - this is where we must fail unless we have the contingency table of the right size!
            if counts.shape == (2, 2):
                pvals[phenotype] = run_fisher_exact(counts)
            elif counts.shape == (3, 2):
                pvals[phenotype] = run_recessive_fisher_exact(counts)
            else:
                raise ValueError(
                    "Invalid number of categories. "
                    f"A {counts.shape} table was created. Only (2, 2) and (3, 2) are valid sizes."
                )

        # 3) Multiple correction
        if self._correction is not None:
            _, pvals_corrected, _, _ = multitest.multipletests(pvals, alpha=self._mtc_alpha, method=self._correction)
            corrected_idx = pd.Index(n_usable.index, name='p_val_corrected')
            corrected_pvals_series = pd.Series(data=pvals_corrected, index=corrected_idx,
                                               name='Corrected p value')
        else:
            corrected_pvals_series = None

        # 4) Wrap up
        return GenotypePhenotypeAnalysisResult(
            n_usable=n_usable,
            all_counts=all_counts,
            pvals=pvals,
            corrected_pvals=corrected_pvals_series,
            phenotype_categories=categories,
            geno_predicate=gt_predicate,
        )
