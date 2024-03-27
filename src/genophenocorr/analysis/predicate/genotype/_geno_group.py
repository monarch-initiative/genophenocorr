import typing

import hpotk

from genophenocorr.model import Patient, FeatureType, VariantEffect
from genophenocorr.preprocessing import ProteinMetadataService

from .._api import PatientCategory, GroupingPredicate


# TODO - should we remove these three?
HETEROZYGOUS = PatientCategory(cat_id=0,
                               name='Heterozygous',
                               description="""
                               This sample has the tested attribute on one allele.
                               """)  #: :meta hide-value:

HOMOZYGOUS = PatientCategory(cat_id=1,
                             name='Homozygous',
                             description="""
                             This sample has the tested attribute on both alleles.
                             """)  #: :meta hide-value:

NO_VARIANT = PatientCategory(cat_id=2,
                             name='No Variant',
                             description="""
                             The sample does not have the tested attribute.
                             """)  #: :meta hide-value:



class VariantEffectsPredicate(GroupingPredicate):
    """
    `VariantEffectsPredicate` tests if the `patient` has at least one variant that is predicted to have
    one of the two functional `effects` on the transcript of interest.

    :param transcript_id: the accession of the transcript of interest.
    :param effect1: the first tested variant effect.
    :param effect2: the second tested variant effect.
    """

    def __init__(self, transcript_id: str,
                 effect1: VariantEffect, effect2: VariantEffect) -> None:
        self._tx_id = transcript_id
        self._effect1 = hpotk.util.validate_instance(effect1, VariantEffect, 'effect1')
        self._effect2 = hpotk.util.validate_instance(effect2, VariantEffect, 'effect2')

    def get_question(self) -> str:
        return f'{self._effect1.name} or {self._effect2.name} on {self._tx_id}'

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        results = [False, False]
        for variant in patient.variants:
            for ann in variant.tx_annotations:
                if ann.transcript_id == self._tx_id:
                    for var_eff in ann.variant_effects:
                        if var_eff == self._effect1:
                            results[0] = True
                        if var_eff == self._effect2:
                            results[1] = True
        if results == [True, False]:
            return GroupingPredicate.FIRST
        elif results == [False, True]:
            return GroupingPredicate.SECOND
        else:
            return None

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'VariantEffectPredicate(transcript_id={self._tx_id}, effect1={self._effect1}, effect2={self._effect2})'


class VariantsPredicate(GroupingPredicate):
    """
    `VariantsPredicate` tests if the `patient` has at least one allele of one of the variants described by
    the `variant_key1` and `variant_key2`.
    **If patient has both variant_keys, it is not included**

    .. note::

      You can get the variant key by calling :class:`genophenocorr.model.VariantCoordinates.variant_key`.

    :param variant_key1: a `str` with a variant key.
    :param variant_key2: a `str` with a variant key.
    """

    def __init__(self, variant_key1: str, variant_key2: str) -> None:
        self._variant_key1 = hpotk.util.validate_instance(variant_key1, str, 'variant_key1')
        self._variant_key2 = hpotk.util.validate_instance(variant_key2, str, 'variant_key2')

    def get_question(self) -> str:
        return f'>=1 allele of either variant {self._variant_key1} or variant {self._variant_key2}'

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        result = [False, False]
        for variant in patient.variants:
            if variant.variant_coordinates.variant_key == self._variant_key1:
                result[0] = True
            if variant.variant_coordinates.variant_key == self._variant_key2:
                result[1] = True
        if result == [True, False]:
            return GroupingPredicate.FIRST
        elif result == [False, True]:
            return GroupingPredicate.SECOND
        else:
            return None

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'VariantPredicate(variant_key1={self._variant_key1}, variant_key2={self._variant_key2})'


class ExonsPredicate(GroupingPredicate):
    """
    `ExonsPredicate` tests if the `patient` has a variant that affects one of the
    two given *n*-th exons of the transcript of interest.

    .. warning::

      We use 1-based numbering to number the exons, not the usual 0-based numbering of the computer science.
      Therefore, the first exon of the transcript has ``exon_number==1``, the second exon is ``2``, and so on ...

    .. warning::

      We do not check if the `exon_number` spans beyond the number of exons of the given `transcript_id`!
      Therefore, ``exon_number==10,000`` will effectively return `None` for *all* patients!!! 😱
      Well, at least the patients of the *Homo sapiens sapiens* taxon...

    :param transcript_id: the accession of the transcript of interest.
    :param exon1_number: a positive `int` of the first target exon.
    :param exon2_number: a positive `int` of the second target exon.
    """

    def __init__(self, transcript_id: str,
                 exon1_number: int, exon2_number: int) -> None:
        self._tx_id = transcript_id
        self._exon1_number = hpotk.util.validate_instance(exon1_number, int, 'exon1_number')
        self._exon2_number = hpotk.util.validate_instance(exon2_number, int, 'exon2_number')
        if self._exon1_number <= 0 or self._exon2_number <= 0:
            raise ValueError(f'`exon1_number` and `exon2_number` must both be a positive `int` but got {(self._exon1_number, self._exon2_number)}')

    def get_question(self) -> str:
        return f'Variant in exon {self._exon1_number} vs in exon {self._exon2_number} on {self._tx_id}'

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        result = [False, False]
        for variant in patient.variants:
            for ann in variant.tx_annotations:
                if ann.transcript_id == self._tx_id:
                    if ann.overlapping_exons is not None:
                        if self._exon1_number in ann.overlapping_exons:
                            result[0] = True
                        if self._exon2_number in ann.overlapping_exons:
                            result[1] = True

        if result == [True, False]:
            return GroupingPredicate.FIRST
        elif result == [False, True]:
            return GroupingPredicate.SECOND
        else:
            return None

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'ExonPredicate(tx_id={self._tx_id}, exon1_number={self._exon1_number}, exon2_number={self._exon2_number})'


class ProtFeatureTypesPredicate(GroupingPredicate):
    """
    `ProtFeatureTypesPredicate` tests if the `patient` has a variant that affects one of the two given
    :class:`FeatureType`s in the transcript of interest.

    :param transcript_id: the accession of the transcript of interest.
    :param feature_type1: an instance of the first target :class:`FeatureType`.
    :param feature_type2: an instance of the second target :class:`FeatureType`.
    """

    def __init__(self, transcript_id: str, feature_type1: FeatureType, feature_type2: FeatureType, protein_service: ProteinMetadataService) -> None:
        self._tx_id = transcript_id
        self._feature_type1 = hpotk.util.validate_instance(feature_type1, FeatureType, 'feature_type')
        self._feature_type2 = hpotk.util.validate_instance(feature_type2, FeatureType, 'feature_type2')
        if feature_type1 == feature_type2:
            raise ValueError(f"Cannot compare against same feature types. {feature_type1} equals {feature_type2}.")
        self._protein_service = hpotk.util.validate_instance(protein_service, ProteinMetadataService, 'protein_service')

    def get_question(self) -> str:
        return f'Variant that affects {self._feature_type1.name} protein feature type vs {self._feature_type2} on {self._tx_id}'

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        result = [False, False]
        for variant in patient.variants:
            for ann in variant.tx_annotations:
                if ann.transcript_id == self._tx_id:
                    prot_loc = ann.protein_effect_location
                    prot_id = ann.protein_id
                    if prot_id is not None and prot_loc is not None:
                        proteins = self._protein_service.annotate(prot_id)
                        for prot in proteins:
                            if prot.protein_id == prot_id:
                                for feat in prot.protein_features:
                                    if feat.feature_type == self._feature_type1:
                                        if prot_loc.overlaps_with(feat.info.region):
                                            result[0] = True
                                    if feat.feature_type == self._feature_type2:
                                        if prot_loc.overlaps_with(feat.info.region):
                                            result[1] = True
        if result == [True, False]:
            return GroupingPredicate.FIRST
        elif result == [False, True]:
            return GroupingPredicate.SECOND
        else:
            return None

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'ProtFeatureTypePredicate(tx_id={self._tx_id}, feature_type1={self._feature_type1}, feature_type2={self._feature_type2}, protein_service={self._protein_service})'


class ProtFeaturesPredicate(GroupingPredicate):
    """
    `ProtFeaturesPredicate` tests if the `patient` has a variant that overlaps with one of the two given
     protein features.

    The predicate needs the name of the protein feature.
    For instance, `EGF-like 2` for `FBN1 <https://www.uniprot.org/uniprotkb/P35555/entry#family_and_domains>`_

    :param transcript_id: the accession of the transcript of interest.
    :param protein_feature1_name: a `str` with the name of the first protein feature.
    :param protein_feature2_name: a `str` with the name of the second protein feature.
    """

    def __init__(self, transcript_id: str, protein_feature1_name: str, protein_feature2_name: str, protein_service: ProteinMetadataService) -> None:
        self._tx_id = transcript_id
        self._pf1_name = hpotk.util.validate_instance(protein_feature1_name, str, 'protein_feature1_name')
        self._pf2_name = hpotk.util.validate_instance(protein_feature2_name, str, 'protein_feature2_name')
        if self._pf1_name == self._pf2_name:
            raise ValueError(f"Cannot compare against the same features. {self._pf1_name} equals {self._pf2_name}.")
        self._protein_service = hpotk.util.validate_instance(protein_service, ProteinMetadataService, 'protein_service')

    def get_question(self) -> str:
        return f'Variant that affects {self._pf1_name} protein feature vs {self._pf2_name} on {self._tx_id}'

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        results = [False, False]
        for variant in patient.variants:
            for ann in variant.tx_annotations:
                if ann.transcript_id == self._tx_id:
                    prot_loc = ann.protein_effect_location
                    prot_id = ann.protein_id
                    if prot_id is not None and prot_loc is not None:
                        proteins = self._protein_service.annotate(prot_id)
                        for prot in proteins:
                            if prot.protein_id == prot_id:
                                for feat in prot.protein_features:
                                    if feat.info.name == self._pf1_name:
                                        if prot_loc.overlaps_with(feat.info.region):
                                            results[0] = True
                                    if feat.info.name == self._pf2_name:
                                        if prot_loc.overlaps_with(feat.info.region):
                                            results[1] = True

        if results == [True, False]:
            return GroupingPredicate.FIRST
        elif results == [False, True]:
            return GroupingPredicate.SECOND
        else:
            return None

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'ProtFeaturePredicate(tx_id={self._tx_id}, protein_feature1_name={self._pf1_name}, protein_feature2_name={self._pf2_name}, protein_service={self._protein_service})'
