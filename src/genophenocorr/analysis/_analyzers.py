import logging
import typing
import warnings

import hpotk

from statsmodels.stats import multitest
from pandas import DataFrame, MultiIndex
from collections import Counter

from genophenocorr.model import Cohort, FeatureType, VariantEffect, Patient

from .predicate import PolyPredicate
from .predicate.genotype import VariantEffectPredicate, VariantPredicate, ExonPredicate, ProtFeatureTypePredicate, ProtFeaturePredicate
# from .predicate.genotype import HOMOZYGOUS, HETEROZYGOUS, NO_VARIANT

from ._api import CohortAnalysis
from ._stats import run_recessive_fisher_exact, run_fisher_exact


class BaseCohortAnalysis(CohortAnalysis):
    """
    `CohortAnalysis` supports running a palette of genotype-phenotype correlation analyses for a :class:`Cohort`.

    The genotype-phenotype correlation is run only for the present HPO terms. The excluded phenotypic features are used
    to calculate the total number of the patients investigated for a term.

    .. note:

      `missing_implies_excluded` option assumes that a patients were investigated for all cohort phenotypic features,
      and unless the feature is not explicitly present, it is marked as excluded and used in calculating
      of the total number of the patients investigated for the feature.


    Attributes:
        analysis_cohort (Cohort): The Cohort object given for this analysis
        analysis_patients (list): A subset of patients with no structural variants if they are being removed from the analysis. Otherwise, all patients in the cohort.
        analysis_transcript (string): The given transcript ID that will be used in this analysis
        is_recessive (boolean): True if the variant is recessive. Default is False.
        missing_implies_excluded (boolean): True if we assume that a patient without a specific phenotype listed
          *does not* have the phenotype. Otherwise, the only excluded phenotypes are those that are excluded explicitly.
          Defaults to `False`.
        include_large_SV (boolean): True if we want to include large structural variants in the analysis. Defaults to True.
        min_perc_patients_w_hpo (int | float): Will only test phenotypes that are listed in at least this percent of patients.
          Can be a `float` in range :math:`(0,1]` or a positive `int` with the minimum count (included). Defaults to `.1` (10%).
    Methods:
        compare_by_variant_type(var_type1:VariantEffect, var_type2:Optional VariantEffect): Runs Fisher Exact analysis, finds any correlation between given variant effects and phenotypes
        compare_by_variant_type(variant1:string, variant2:Optional string): Runs Fisher Exact analysis, finds any correlation between given variants and phenotypes
        compare_by_variant_type(exon1:integer, exon2:Optional integer): Runs Fisher Exact analysis, finds any correlation between given affected exons and phenotypes
        compare_by_variant_type(feature1:FeatureType, feature2:Optional FeatureType): Runs Fisher Exact analysis, finds any correlation between given protein feature types and phenotypes
        compare_by_variant_type(feature1:string, feature2:Optional string): Runs Fisher Exact analysis, finds any correlation between given protein features and phenotypes
    """
    def __init__(self, cohort: Cohort, transcript: str, hpo: hpotk.MinimalOntology, recessive: bool = False,
                 missing_implies_excluded: bool = False, p_val_correction: str = 'bonferroni',
                include_large_SV: bool = True, min_perc_patients_w_hpo: typing.Union[float, int] = .1) -> None:
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort must be type Cohort but was type {type(cohort)}")
        if transcript not in cohort.all_transcripts:
            raise ValueError(f"Transcript {transcript} not found in Cohort")
        self._logger = logging.getLogger(__name__)
        self._cohort = cohort
        self._transcript = transcript
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._recessive = recessive
        self._missing_implies_excluded = missing_implies_excluded ##TODO: Figure out how to move this out of analysis
        self._include_SV = include_large_SV
        self._patient_list = list(cohort.all_patients) if include_large_SV else [pat for pat in cohort.all_patients if not all(var.variant_coordinates.is_structural() for var in pat.variants)]
        if len(self._patient_list) == 0:
            raise ValueError('No patients left for analysis!')
        self._percent_pats = self._check_min_perc_patients_w_hpo(min_perc_patients_w_hpo, len(self._patient_list))

        # As of now, `self._testing_hpo_terms` is a Sequence[TermId] with terms that were *present* in >=1 cohort member.
        self._testing_hpo_terms = self._remove_low_hpo_terms()
        self._patients_by_hpo = self._group_patients_by_hpo(self._testing_hpo_terms, self._patient_list,
                                                            self._hpo, self._missing_implies_excluded)
        self._correction = p_val_correction

    @property
    def analysis_cohort(self) -> Cohort:
        """
        Returns:
            Cohort: The Cohort object being analyzed
        """
        # TODO[0.3.0] - remove
        # This is because we may actually filter out some patients in the constructor and, in result,
        # not analyze the ENTIRE cohort. Exposing the cohort is, thus, superfluous and unless there are strong reasons
        # for keeping the cohort, I think it should be removed.
        warnings.warn('The `analysis_cohort` will be removed from the CohortAnalyzer API in 0.3.0',
                      DeprecationWarning, stacklevel=2)
        return self._cohort

    @property
    def analysis_patients(self) -> typing.Sequence[Patient]:
        """
        Returns:
            Sequence: A sequence of patients that are subject to the analysis.
        """
        return self._patient_list

    @property
    def analysis_transcript(self) -> str:
        """
        Returns:
            string: The accession ID of the transcript that will be used in this analysis.
        """
        return self._transcript

    @property
    def is_recessive(self) -> bool:
        """
        Returns:
            boolean: True if the analysis is assuming a disease with the Autosomal recessive mode of inheritance.
        """
        return self._recessive

    @property
    def include_unmeasured(self) -> bool:
        """
        Returns:
            boolean: True if we want to assume a patient without a specific phenotype listed does not have that phenotype.
            Otherwise, will only include those that explicitly say the phenotype was not observed.
        """
        # TODO[0.3.0] - remove
        warnings.warn('The `include_unmeasured` will be removed from the CohortAnalyzer API in `v0.3.0`.'
                      'Use `missing_implies_excluded` instead',
                      DeprecationWarning, stacklevel=2)
        return self.missing_implies_excluded

    @property
    def missing_implies_excluded(self) -> bool:
        """
        Returns:
            boolean: True if we assume that patient without a specific phenotype listed *does not* have the phenotype.
            Otherwise, the only excluded phenotypes are those that are excluded explicitly.
        """
        return self._missing_implies_excluded

    @property
    def include_large_SV(self):
        """
        Returns:
            boolean: True if we want to include large structural variants in the analysis.
        """
        return self._include_SV

    @property
    def min_perc_patients_w_hpo(self):
        """
        Returns:
            integer: Will only test phenotypes that are listed in at least this percent of patients.
        """
        # TODO[0.3.0] - remove
        warnings.warn('The `min_perc_patients_w_hpo` will be removed from the CohortAnalyzer API in 0.3.0',
                      DeprecationWarning, stacklevel=2)
        return round(self._percent_pats * 100)

    def _remove_low_hpo_terms(self) -> typing.Sequence[hpotk.TermId]:
        final_hpo_list = set()
        all_hpo_counts = Counter()
        all_no_hpo_counts = Counter()
        for pat in self.analysis_patients:
            all_hpo_counts.update(pat.present_phenotypes())
            if not self.include_unmeasured:
                all_no_hpo_counts.update(pat.excluded_phenotypes())
        if self.include_unmeasured:
            for hpo, count in all_hpo_counts.items():
                count = count if not None else 0
                if (count / len(self.analysis_patients)) >= (self.min_perc_patients_w_hpo / 100):
                    final_hpo_list.add(hpo)
        else:
            for hpo, count in all_hpo_counts.items():
                count = count if count is not None else 0
                no_count = all_no_hpo_counts.get(hpo) if all_no_hpo_counts.get(hpo) is not None else 0
                if (count / (count + no_count)) >= (self.min_perc_patients_w_hpo / 100):
                    final_hpo_list.add(hpo)
        if len(final_hpo_list) == 0:
            raise ValueError(f"No HPO terms found in over {self.min_perc_patients_w_hpo}% of patients.")

        return [p.identifier for p in final_hpo_list]

    def _run_dom_analysis(self, predicate: PolyPredicate, variable1, variable2):
        final_dict = dict()
        all_pvals = []
        # if not self.is_recessive:
        #     for hpo in self._testing_hpo_terms:
        #         with_hpo_var1_count = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) in (HOMOZYGOUS, HETEROZYGOUS)])
        #         not_hpo_var1_count = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) in (HOMOZYGOUS, HETEROZYGOUS)])
        #         if variable2 is None:
        #             with_hpo_var2_count = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) == NO_VARIANT])
        #             not_hpo_var2_count = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) == NO_VARIANT])
        #         else:
        #             with_hpo_var2_count = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable2) in (HOMOZYGOUS, HETEROZYGOUS)])
        #             not_hpo_var2_count = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable2) in (HOMOZYGOUS, HETEROZYGOUS)])
        #         if with_hpo_var1_count + not_hpo_var1_count == 0 or with_hpo_var2_count + not_hpo_var2_count == 0:
        #             self._logger.warning(f"Divide by 0 error with HPO {hpo.value}, not included in this analysis.")
        #             continue
        #         p_val = run_fisher_exact([[with_hpo_var1_count, not_hpo_var1_count], [with_hpo_var2_count, not_hpo_var2_count]])
        #         all_pvals.append(p_val)
        #         term_name = self._hpo.get_term(hpo).name
        #         final_dict[f"{hpo.value} ({term_name})"] = [with_hpo_var1_count, "{:0.2f}%".format((with_hpo_var1_count/(with_hpo_var1_count + not_hpo_var1_count)) * 100), with_hpo_var2_count, "{:0.2f}%".format((with_hpo_var2_count/(with_hpo_var2_count + not_hpo_var2_count)) * 100), p_val]
        # else:
        #     return ValueError(f"Run a recessive analysis")
        corrected_pvals = multitest.multipletests(all_pvals, method = self._correction)[1]
        return final_dict, corrected_pvals

    def _run_rec_analysis(self, predicate, variable1, variable2):
        final_dict = dict()
        all_pvals = []
        # if self.is_recessive:
            # for hpo in self._testing_hpo_terms:
                # if variable2 is None:
                #     with_hpo_var1_var1 = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) == HOMOZYGOUS])
                #     no_hpo_var1_var1 = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) == HOMOZYGOUS])
                #     with_hpo_var1_var2 = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) == HETEROZYGOUS])
                #     no_hpo_var1_var2 = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) == HETEROZYGOUS])
                #     with_hpo_var2_var2 = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) == NO_VARIANT])
                #     no_hpo_var2_var2 = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) == NO_VARIANT])
                # else:
                #     with_hpo_var1_var1 = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) == HOMOZYGOUS])
                #     no_hpo_var1_var1 = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) == HOMOZYGOUS])
                #     with_hpo_var1_var2 = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) == HETEROZYGOUS and predicate.test(pat, variable2) == HETEROZYGOUS])
                #     no_hpo_var1_var2 = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) == HETEROZYGOUS and predicate.test(pat, variable2) == HETEROZYGOUS])
                #     with_hpo_var2_var2 = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable2) == HOMOZYGOUS])
                #     no_hpo_var2_var2 = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable2) == HOMOZYGOUS])
                # if with_hpo_var1_var1 + no_hpo_var1_var1 == 0 or with_hpo_var1_var2 + no_hpo_var1_var2 == 0 or with_hpo_var2_var2 + no_hpo_var2_var2 == 0:
                #     self._logger.warning(f"Divide by 0 error with HPO {hpo.value}, not included in this analysis.")
                #     continue
                # p_val = run_recessive_fisher_exact([[with_hpo_var1_var1, no_hpo_var1_var1], [with_hpo_var1_var2, no_hpo_var1_var2], [with_hpo_var2_var2, no_hpo_var2_var2]] )
                # all_pvals.append(p_val)
                # term_name = self._hpo.get_term(hpo).name
                # final_dict[f"{hpo.value} ({term_name})"] = [with_hpo_var1_var1, "{:0.2f}%".format((with_hpo_var1_var1/(with_hpo_var1_var1 + no_hpo_var1_var1)) * 100), with_hpo_var1_var2, "{:0.2f}%".format((with_hpo_var1_var2/(no_hpo_var1_var2 + with_hpo_var1_var2)) * 100), with_hpo_var2_var2, "{:0.2f}%".format((with_hpo_var2_var2/(with_hpo_var2_var2 + no_hpo_var2_var2)) * 100), p_val]
        # else:
        #     return ValueError("Run a dominant analysis")
        corrected_pvals = multitest.multipletests(all_pvals, method = self._correction)[1]
        return final_dict, corrected_pvals

    def compare_by_variant_effect(self, effect: VariantEffect, tx_id: str):
        # We need to cache the previous transcript
        previous = self._transcript
        self._transcript = tx_id
        result = self.compare_by_variant_type(effect)
        self._transcript = previous
        return result

    def compare_by_variant_type(self, var_type1:VariantEffect, var_type2:typing.Optional[VariantEffect] = None):
        """Runs Fisher Exact analysis, finds any correlation between given variant effects across phenotypes.

        Args:
            var_type1 (VariantEffect):
            var_type2 (VariantEffect, Optional): If None, we compare between those with var_type1 and those without var_type1
        Returns:
            DataFrame: A pandas DataFrame showing the results of the analysis
        """
        var_effect_test = VariantEffectPredicate(self._transcript, var_type1)
        if self.is_recessive:
            final_dict, corrected_pvals = self._run_rec_analysis(var_effect_test, var_type1, var_type2)
            if var_type2 is None:
                col_name1 = f'Heterozygous {str(var_type1)}'
                col_name2 = f'No {str(var_type1)}'
            else:
                col_name1 = f'Heterozygous {str(var_type1)} and {str(var_type2)}'
                col_name2 = f'Homozygous {str(var_type2)}'
            index = MultiIndex.from_product([[f'Homozygous {str(var_type1)}', col_name1, col_name2], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(6, ('', "p-value")))
            final_df.insert(7, ('', "Corrected p-values"), corrected_pvals, True)
        else:
            final_dict, corrected_pvals = self._run_dom_analysis(var_effect_test, var_type1, var_type2)
            if var_type2 is None:
                col_name = f'Without {str(var_type1)}'
            else:
                col_name = f'With {str(var_type2)}'
            index = MultiIndex.from_product([[f'With {str(var_type1)}', col_name], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(4, ('', 'p-value')))
            final_df.insert(5, ('', "Corrected p-values"), corrected_pvals, True)
        return final_df.sort_values([('', 'Corrected p-values'), ('', 'p-value')])

    def compare_by_variant_key(self, variant_key: str):
        return self.compare_by_variant(variant_key)

    def compare_by_variant(self, variant1:str, variant2:str = None):
        """Runs Fisher Exact analysis, finds any correlation between given variants across phenotypes.

        Args:
            variant1 (string):
            variant2 (string, Optional): If None, we compare between those with variant1 and those without variant1
        Returns:
            DataFrame: A pandas DataFrame showing the results of the analysis
        """
        variant_test = VariantPredicate(self.analysis_transcript)
        if self.is_recessive:
            final_dict, corrected_pvals = self._run_rec_analysis(variant_test, variant1, variant2)
            if variant2 is None:
                col_name1 = f'Heterozygous {variant1}'
                col_name2 = f'No {variant1}'
            else:
                col_name1 = f'Heterozygous {variant1} and {variant2}'
                col_name2 = f'Homozygous {variant2}'
            index = MultiIndex.from_product([[f'Homozygous {variant1}', col_name1, col_name2], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(6, ('', 'p-value')))
            final_df.insert(7, ('', "Corrected p-values"), corrected_pvals, True)
        else:
            final_dict, corrected_pvals = self._run_dom_analysis(variant_test, variant1, variant2)
            if variant2 is None:
                col_name = f'Without {variant1}'
            else:
                col_name = f'With {variant2}'
            index = MultiIndex.from_product([[f'With {variant1}', col_name], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(4, ('', 'p-value')))
            final_df.insert(5, ('', "Corrected p-values"), corrected_pvals, True)
        return final_df.sort_values([('', 'Corrected p-values'), ('', 'p-value')])

    def compare_by_exon(self, exon1:int, exon2:typing.Optional[int] = None):
        """Runs Fisher Exact analysis, finds any correlation between given exons across phenotypes.

        Args:
            exon1 (integer):
            exon2 (integer, Optional): If None, we compare between those within exon1 and those outside exon1
        Returns:
            DataFrame: A pandas DataFrame showing the results of the analysis
        """
        exon_test = ExonPredicate(self.analysis_transcript, exon1)
        if self.is_recessive:
            final_dict, corrected_pvals = self._run_rec_analysis(exon_test, exon1, exon2)
            if exon2 is None:
                col_name1 = f'Heterozygous {exon1}'
                col_name2 = f'No {exon1}'
            else:
                col_name1 = f'Heterozygous {exon1} and {exon2}'
                col_name2 = f'Homozygous {exon2}'
            index = MultiIndex.from_product([[f'Homozygous {exon1}', col_name1, col_name2], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(6, ('', 'p-value')))
            final_df.insert(7, ('', "Corrected p-values"), corrected_pvals, True)
        else:
            final_dict, corrected_pvals = self._run_dom_analysis(exon_test, exon1, exon2)
            if exon2 is None:
                col_name = f'Outside Exon {exon1}'
            else:
                col_name = f'Inside Exon {exon2}'
            index = MultiIndex.from_product([[f'Inside Exon {exon1}', col_name], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(4, ('', 'p-value')))
            final_df.insert(5, ('', "Corrected p-values"), corrected_pvals, True)
        return final_df.sort_values([('', 'Corrected p-values'), ('', 'p-value')])

    def compare_by_protein_feature_type(self, feature1:FeatureType, feature2:typing.Optional[FeatureType] = None):
        """Runs Fisher Exact analysis, finds any correlation between given feature type across phenotypes.

        Args:
            feature1 (FeatureType):
            feature2 (FeatureType, Optional): If None, we compare between those with feature1 and those without feature1
        Returns:
            DataFrame: A pandas DataFrame showing the results of the analysis
        """
        feat_type_test = ProtFeatureTypePredicate(self.analysis_transcript, feature1)
        if self.is_recessive:
            final_dict, corrected_pvals = self._run_rec_analysis(feat_type_test, feature1, feature2)
            if feature2 is None:
                col_name1 = f'Heterozygous {feature1}'
                col_name2 = f'No {feature1}'
            else:
                col_name1 = f'Heterozygous {feature1} and {feature2}'
                col_name2 = f'Homozygous {feature2}'
            index = MultiIndex.from_product([[f'Homozygous {feature1}', col_name1, col_name2], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(6, ('', 'p-value')))
            final_df.insert(7, ('', "Corrected p-values"), corrected_pvals, True)
        else:
            final_dict, corrected_pvals = self._run_dom_analysis(feat_type_test, feature1, feature2)
            if feature2 is None:
                col_name = f'Outside {feature1.name}'
            else:
                col_name = f'Inside {feature2.name}'
            index = MultiIndex.from_product([[f'Inside {feature1.name}', col_name], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(4, ('', 'p-value')))
            final_df.insert(5, ('', "Corrected p-values"), corrected_pvals, True)
        return final_df.sort_values([('', 'Corrected p-values'), ('', 'p-value')])


    def compare_by_protein_feature(self, feature1:str, feature2:typing.Optional[str] = None):
        """Runs Fisher Exact analysis, finds any correlation between given feature and phenotypes.

        Args:
            feature1 (string):
            feature2 (string, Optional): If None, we compare between those within feature1 and those outside feature1
        Returns:
            DataFrame: A pandas DataFrame showing the results of the analysis
        """
        domain_test = ProtFeaturePredicate(self.analysis_transcript, feature1)
        if self.is_recessive:
            final_dict, corrected_pvals = self._run_rec_analysis(domain_test, feature1, feature2)
            if feature2 is None:
                col_name1 = f'Heterozygous {feature1}'
                col_name2 = f'No {feature1}'
            else:
                col_name1 = f'Heterozygous {feature1} and {feature2}'
                col_name2 = f'Homozygous {feature2}'
            index = MultiIndex.from_product([[f'Homozygous {feature1}', col_name1, col_name2], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(6, ('', 'p-value')))
            final_df.insert(7, ('', "Corrected p-values"), corrected_pvals, True)
        else:
            final_dict, corrected_pvals = self._run_dom_analysis(domain_test, feature1, feature2)
            if feature2 is None:
                col_name = f'Outside {feature1}'
            else:
                col_name = f'Inside {feature2}'
            index = MultiIndex.from_product([[f'Inside {feature1}', col_name], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(4, ('', 'p-value')))
            final_df.insert(5, ('', "Corrected p-values"), corrected_pvals, True)
        return final_df.sort_values([('', 'Corrected p-values'), ('', 'p-value')])


    def compare_by_variant_keys(self, a: str, b: str):
        raise NotImplementedError()
