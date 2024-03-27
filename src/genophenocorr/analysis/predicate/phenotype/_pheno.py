import abc
import typing

import hpotk

from genophenocorr.model import Patient

from .._api import PolyPredicate, PatientCategory, PatientCategories, C, Categorization


class Comparable(metaclass=abc.ABCMeta):
    """
    A protocol for annotating comparable types.
    """

    @abc.abstractmethod
    def __eq__(self, other: typing.Any) -> bool:
        pass

    @abc.abstractmethod
    def __lt__(self, other: typing.Any) -> bool:
        pass


P = typing.TypeVar('P', typing.Hashable, Comparable)
"""
Phenotype entity of interest, such as :class:`hpotk.model.TermId`, representing an HPO term or an OMIM/MONDO term.

However, phenotype entity can be anything as long as it is hashable and comparable.
"""


class PhenotypeCategorization(typing.Generic[P], Categorization):
    """"
    On top of the attributes of the `Categorization`, `PhenotypeCategorization` keeps track of the target phenotype `P`.
    """

    def __init__(
            self,
            category: PatientCategory,
            phenotype: P,
    ):
        super().__init__(category)
        self._phenotype = phenotype

    @property
    def phenotype(self) -> P:
        return self._phenotype

    def __repr__(self) -> str:
        return f"PhenotypeCategorization(" \
               f"category={self._category}, " \
               f"phenotype={self._phenotype})"

    def __str__(self) -> str:
        return repr(self)

    def __eq__(self, other) -> bool:
        return isinstance(other, PhenotypeCategorization) \
            and self._category == other._category \
            and self._phenotype == other._phenotype

    def __hash__(self) -> int:
        return hash((self._category, self._phenotype))


class PhenotypePolyPredicate(typing.Generic[P], PolyPredicate[PhenotypeCategorization[P]], metaclass=abc.ABCMeta):
    """
    `PhenotypePolyPredicate` investigates a patient in context of one or more phenotype categories `P`.

    Each phenotype category `P` can be a :class:`hpotk.model.TermId` representing an HPO term or an OMIM/MONDO term.

    Most of the time, only one category is investigated, and :attr:`phenotype` will return a sequence with
    one element only (e.g. *Arachnodactyly* `HP:0001166`). However, multiple categories can be sought as well,
    such as when the predicate bins the patient into one of discrete diseases
    (e.g. Geleophysic dysplasia, Marfan syndrome, ...). Then the predicate should return a sequence of disease
    identifiers.
    """

    @property
    @abc.abstractmethod
    def phenotypes(self) -> typing.Sequence[P]:
        """
        Get the phenotype entities of interest.
        """
        pass


class PropagatingPhenotypePredicate(PhenotypePolyPredicate[PhenotypeCategorization[hpotk.TermId]]):
    """
    `PropagatingPhenotypePredicate` tests if a patient is annotated with an HPO term.

    Note, `query` must be a term of the provided `hpo`!

    :param hpo: HPO object
    :param query: the HPO term to test
    :param missing_implies_excluded: `True` if lack of an explicit annotation implies term's absence`.
    """

    def __init__(self, hpo: hpotk.MinimalOntology,
                 query: hpotk.TermId,
                 missing_implies_excluded: bool = False):
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._query = hpotk.util.validate_instance(query, hpotk.TermId, 'phenotypic_feature')
        self._query_label = self._hpo.get_term_name(query)
        assert self._query_label is not None, f'Query {query} is in HPO'
        self._missing_implies_excluded = hpotk.util.validate_instance(missing_implies_excluded, bool,
                                                                      'missing_implies_excluded')
        self._present = PhenotypeCategorization(
            category=PatientCategories.YES,
            phenotype=query,
        )
        self._excluded = PhenotypeCategorization(
            category=PatientCategories.NO,
            phenotype=query,
        )

    def get_question(self) -> str:
        return f'Is {self._query_label} present in the patient?'

    @property
    def phenotypes(self) -> typing.Sequence[hpotk.TermId]:
        # We usually test just a single HPO term, so we return a tuple with a single member.
        return self._query,

    def get_categorizations(self) -> typing.Sequence[C]:
        return self._present, self._excluded

    def test(self, patient: Patient) -> typing.Optional[PhenotypeCategorization[P]]:
        self._check_patient(patient)

        if len(patient.phenotypes) == 0:
            return None

        for phenotype in patient.phenotypes:
            if phenotype.is_present:
                if any(self._query == anc for anc in self._hpo.graph.get_ancestors(phenotype, include_source=True)):
                    return self._present
            else:
                if self._missing_implies_excluded:
                    return self._excluded
                else:
                    if any(self._query == desc for desc in
                           self._hpo.graph.get_descendants(phenotype, include_source=True)):
                        return self._excluded

        return None

    def __repr__(self):
        return f'PropagatingPhenotypeBooleanPredicate(query={self._query})'


class DiseasePresencePredicate(PhenotypePolyPredicate[PhenotypeCategorization[hpotk.TermId]]):
    """
    `DiseasePresencePredicate` tests if the patient was diagnosed with a disease.

    The predicate tests if the patient's diseases include a disease ID formatted as a :class:`hpotk.model.TermId`.
    """
    # TODO: Lauren please implement and test.
    pass
