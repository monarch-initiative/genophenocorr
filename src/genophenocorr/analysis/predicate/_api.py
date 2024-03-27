import abc
import typing

import hpotk.util

from genophenocorr.model import Patient


class PatientCategory:
    """
    `PatientCategory` represents one of several exclusive discrete groups.

    Patient category has :attr:`cat_id`, a unique numeric identifier of the group,
    :attr:`name` with human-readable group name, and :attr:`description` with an optional verbose description.
    """

    def __init__(self, cat_id: int,
                 name: str,
                 description: typing.Optional[str] = None):
        self._cat_id = hpotk.util.validate_instance(cat_id, int, 'cat_id')
        self._name = hpotk.util.validate_instance(name, str, 'name')
        self._description = hpotk.util.validate_optional_instance(description, str, 'description')

    @property
    def cat_id(self) -> int:
        """
        Get an `int` with the unique numeric identifier of the group.
        """
        return self._cat_id

    @property
    def name(self) -> str:
        """
        Get a `str` with a human-readable name of the group.
        """
        return self._name

    @property
    def description(self) -> typing.Optional[str]:
        """
        Get a `str` with an optional detailed group description.
        """
        return self._description

    def __repr__(self) -> str:
        return f"PatientCategory(cat_id={self.cat_id}, " \
               f"name={self.name}, " \
               f"description={self.description})"

    def __str__(self) -> str:
        return self._name

    def __eq__(self, other) -> bool:
        return isinstance(other, PatientCategory) \
            and self.cat_id == other.cat_id \
            and self.name == other.name \
            and self.description == other.description

    def __hash__(self) -> int:
        return hash((self.cat_id, self.name, self.description))


class PatientCategories(metaclass=abc.ABCMeta):
    """
    A static utility class to serve common patient categories.
    """

    YES = PatientCategory(1, 'Yes', 'The patient belongs to the group.')
    """
    Category for a patient who *belongs* to the tested group.
    """

    NO = PatientCategory(0, 'No', 'The patient does not belong to the group.')
    """
    Category for a patient who does *not* belong to the tested group.
    """


class Categorization:
    """
    `Categorization` represents one of discrete group a :class:`genophenocorr.model.Patient` can be assigned into.
    """

    def __init__(
            self,
            category: PatientCategory,
    ):
        self._category = hpotk.util.validate_instance(category, PatientCategory, 'category')

    @property
    def category(self) -> PatientCategory:
        return self._category

    def __eq__(self, other):
        return isinstance(other, Categorization) and self._category == other._category

    def __hash__(self):
        return hash((self._category,))

    def __repr__(self):
        return f'Categorization(category={self._category})'

    def __str__(self):
        return repr(self)


C = typing.TypeVar('C', bound=Categorization)
"""
A generic bound for types that extend :class:`Categorization`.
"""


class PolyPredicate(typing.Generic[C], metaclass=abc.ABCMeta):
    """
    `PolyPredicate` bins a :class:`genophenocorr.model.Patient` into one of several discrete groups represented
    by :class:`Categorization`.

    The groups must be *exclusive* - the patient can be binned into one and only one group,
    and *exhaustive* - the groups must cover all possible scenarios.

    However, if the patient cannot be assigned into any meaningful category, `None` can be returned.

    Note, that `PolyPredicate` must be used on a happy path - testing a patient must inherently make sense.
    Predicate will *not* check if, for instance, the patient variants are compatible with a certain mode of inheritance.
    """

    @abc.abstractmethod
    def get_categorizations(self) -> typing.Sequence[C]:
        """
        Get a sequence of all categories which the `PolyPredicate` can produce.
        """
        pass

    def get_categories(self) -> typing.Sequence[PatientCategory]:
        """
        Get a sequence with :class:`PatientCategory` instances that the predicate can produce.
        """
        return tuple(c.category for c in self.get_categorizations())

    @abc.abstractmethod
    def get_question(self) -> str:
        """
        Prepare a `str` with the question the predicate can answer.
        """
        pass

    @abc.abstractmethod
    def test(self, patient: Patient) -> typing.Optional[C]:
        """
        Assign a `patient` into a categorization.

        Return `None` if the patient cannot be assigned into any meaningful category.
        """
        pass

    def n_categorizations(self) -> int:
        """
        Get the number of categorizations the predicate can produce.
        """
        return len(self.get_categorizations())

    @staticmethod
    def _check_patient(patient: Patient):
        """
        Check if the `patient` meets the predicate requirements.
        """
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")


class GenotypePolyPredicate(PolyPredicate[Categorization], metaclass=abc.ABCMeta):
    """
    `GenotypePolyPredicate` constrains `PolyPredicate` to investigate the genotype aspects
    of patients.
    """
    pass


class GenotypeBooleanPredicate(GenotypePolyPredicate, metaclass=abc.ABCMeta):
    """
    `GenotypeBooleanPredicate` tests if a :class:`genophenocorr.model.Patient` belongs to a genotype group
     and returns a boolean binning.
    """
    YES = Categorization(PatientCategories.YES)
    NO = Categorization(PatientCategories.NO)

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        """
        The predicate bins a patient into :class:`BooleanPredicate.NO` or :class:`BooleanPredicate.YES` category.
        """
        return GenotypeBooleanPredicate.YES, GenotypeBooleanPredicate.NO


class GroupingPredicate(GenotypePolyPredicate, metaclass=abc.ABCMeta):
    """
    `GroupingPredicate` tests if a :class:`genophenocorr.model.Patient` belongs to one of two groups and returns
    FIRST or SECOND based on which group Patient belongs in.
    """

    FIRST = Categorization(PatientCategory(0, 'First', 'The patient belongs in the first group.'))
    """
    Category for a patient who belongs in the first given tested group.
    """
    SECOND = Categorization(PatientCategory(1, 'Second', 'The patient belongs in the second group.'))
    """
    Category for a patient who belongs to the second given tested group.
    """

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        """
        The predicate bins a patient into :class:`GroupingPredicate.FIRST` or :class:`GroupingPredicate.SECOND` category.
        """
        return GroupingPredicate.FIRST, GroupingPredicate.SECOND
