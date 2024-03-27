import typing

import hpotk

from genophenocorr.model import Patient, FeatureType, VariantEffect
from genophenocorr.preprocessing import ProteinMetadataService

from .._api import Categorization, GenotypeBooleanPredicate


class VariantEffectPredicate(GenotypeBooleanPredicate):
    """
    `VariantEffectPredicate` tests if the `patient` has at least one variant that is predicted to have
    the functional `effect` on the transcript of interest.

    :param transcript_id: the accession of the transcript of interest.
    :param effect: the tested variant effect.
    """

    def __init__(self, transcript_id: str,
                 effect: VariantEffect) -> None:
        self._tx_id = transcript_id
        self._effect = hpotk.util.validate_instance(effect, VariantEffect, 'effect')

    def get_question(self) -> str:
        return f'{self._effect.name} on {self._tx_id}'

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        for variant in patient.variants:
            for ann in variant.tx_annotations:
                if ann.transcript_id == self._tx_id:
                    for var_eff in ann.variant_effects:
                        if var_eff == self._effect:
                            return GenotypeBooleanPredicate.YES

        return GenotypeBooleanPredicate.NO

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'VariantEffectPredicate(transcript_id={self._tx_id}, effect={self._effect})'


class VariantPredicate(GenotypeBooleanPredicate):
    """
    `VariantPredicate` tests if the `patient` has ar least one allele of the variant described by the `variant_key`.

    .. note::

      You can get the variant key by calling :class:`genophenocorr.model.VariantCoordinates.variant_key`.

    :param variant_key: a `str` with the variant key.
    """

    def __init__(self, variant_key: str) -> None:
        self._variant_key = hpotk.util.validate_instance(variant_key, str, 'variant_key')

    def get_question(self) -> str:
        return f'>=1 allele of the variant {self._variant_key}'

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        for variant in patient.variants:
            if variant.variant_coordinates.variant_key == self._variant_key:
                return GenotypeBooleanPredicate.YES

        return GenotypeBooleanPredicate.NO

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'VariantPredicate(variant_key={self._variant_key})'


class ExonPredicate(GenotypeBooleanPredicate):
    """
    `ExonPredicate` tests if the `patient` has a variant that affects *n*-th exon of the transcript of interest.

    .. warning::

      We use 1-based numbering to number the exons, not the usual 0-based numbering of the computer science.
      Therefore, the first exon of the transcript has ``exon_number==1``, the second exon is ``2``, and so on ...

    .. warning::

      We do not check if the `exon_number` spans beyond the number of exons of the given `transcript_id`!
      Therefore, ``exon_number==10,000`` will effectively return :attr:`GenotypeBooleanPredicate.FALSE`
      for *all* patients!!! 😱
      Well, at least the patients of the *Homo sapiens sapiens* taxon...

    :param transcript_id: the accession of the transcript of interest.
    :param exon_number: a positive `int` of the target exon.
    """

    def __init__(self, transcript_id: str,
                 exon_number: int) -> None:
        self._tx_id = transcript_id
        self._exon_number = hpotk.util.validate_instance(exon_number, int, 'exon_number')
        if self._exon_number <= 0:
            raise ValueError(f'`exon_number` must be a positive `int` but got {self._exon_number}')

    def get_question(self) -> str:
        return f'Variant in exon {self._exon_number} on {self._tx_id}'

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        for variant in patient.variants:
            for ann in variant.tx_annotations:
                if ann.transcript_id == self._tx_id:
                    if ann.overlapping_exons is not None:
                        if self._exon_number in ann.overlapping_exons:
                            return GenotypeBooleanPredicate.YES

        return GenotypeBooleanPredicate.NO

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'ExonPredicate(tx_id={self._tx_id}, exon_number={self._exon_number})'


class ProtFeatureTypePredicate(GenotypeBooleanPredicate):
    """
    `ProtFeatureTypePredicate` tests if the `patient` has a variant that affects a :class:`FeatureType`
    in the protein encoded by the transcript of interest.

    :param transcript_id: the accession of the transcript of interest.
    :param feature_type: an instance of the target :class:`FeatureType`.
    :param protein_service: an instance of a :class:`ProteinMetadataService`.
    """

    def __init__(self, transcript_id: str, feature_type: FeatureType, protein_service: ProteinMetadataService) -> None:
        self._tx_id = transcript_id
        self._feature_type = hpotk.util.validate_instance(feature_type, FeatureType, 'feature_type')
        self._protein_service = hpotk.util.validate_instance(protein_service, ProteinMetadataService, 'protein_service')

    def get_question(self) -> str:
        return f'Variant that affects {self._feature_type.name} feature type on protein encoded by transcript {self._tx_id}'

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

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
                                    if feat.feature_type == self._feature_type:
                                        if prot_loc.overlaps_with(feat.info.region):
                                            return GenotypeBooleanPredicate.YES

        return GenotypeBooleanPredicate.NO
        #TODO: Add a logger field, add a branch that handles the state where prot_id is set but prot_loc is not - gives warning

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'ProtFeatureTypePredicate(tx_id={self._tx_id}, feature_type={self._feature_type})'


class ProtFeaturePredicate(GenotypeBooleanPredicate):
    """
    `ProtFeaturePredicate` tests if the `patient` has a variant that overlaps with a protein feature.

    The predicate needs the name of the protein feature.
    For instance, `EGF-like 2` for `FBN1 <https://www.uniprot.org/uniprotkb/P35555/entry#family_and_domains>`_

    :param transcript_id: the accession of the transcript of interest.
    :param protein_feature_name: a `str` with the name of the protein feature.
    """

    def __init__(self, transcript_id: str, protein_feature_name: str, protein_service:ProteinMetadataService) -> None:
        self._tx_id = transcript_id
        self._pf_name = hpotk.util.validate_instance(protein_feature_name, str, 'protein_feature_name')
        self._protein_service = hpotk.util.validate_instance(protein_service, ProteinMetadataService, 'protein_service')

    def get_question(self) -> str:
        return f'Variant that affects {self._pf_name} feature on protein encoded by transcript {self._tx_id}'

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

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
                                    if feat.info.name == self._pf_name:
                                        if prot_loc.overlaps_with(feat.info.region):
                                            return GenotypeBooleanPredicate.YES

        return GenotypeBooleanPredicate.NO

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'ProtFeaturePredicate(tx_id={self._tx_id}, exon_number={self._pf_name})'
