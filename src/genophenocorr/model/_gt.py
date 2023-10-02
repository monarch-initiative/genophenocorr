import abc
import bisect
import enum
import typing
from typing import Iterator

import numpy as np


class Genotype(enum.Enum):
    """
    `Genotype` represents state of a variable locus in a diploid genome.
    """

    NO_CALL = ('.',)
    HOMOZYGOUS_REFERENCE = ('0/0',)
    HETEROZYGOUS = ('0/1',)
    HOMOZYGOUS_ALTERNATE = ('1/1',)
    HEMIZYGOUS = ('1',)

    def __init__(self, code: str):
        self._code = code

    @property
    def code(self) -> str:
        return self._code

    def __str__(self):
        return self.code


class Genotypes(typing.Sized, typing.Iterable):
    """
    `Genotypes` is a container for mapping between sample ID and its genotype.

    Use one of the static methods to create an instance:

    >>> gts = Genotypes.single('A', Genotype.HETEROZYGOUS)
    >>> gts = Genotypes.from_mapping({'A': Genotype.HETEROZYGOUS, 'B': Genotype.HOMOZYGOUS_ALTERNATE})

    There are 2 genotypes in the container:

    >>> len(gts)
    2

    You can get a genotype for a sample ID:

    >>> gts.for_sample('A')
    Genotype.HETEROZYGOUS

    You will get `None` if the sample is not present:

    >>> gts.for_sample('UNKNOWN')
    None

    You can iterate over sample-genotype pairs:

    >>> for sample_id, gt in gts:
    ...   print(sample_id, gt)
    A 0/1
    B 1/1
    """

    @staticmethod
    def empty():
        return EMPTY

    @staticmethod
    def single(sample_id: str, genotype: Genotype):
        """
        A shortcut for creating `Genotypes` for a single sample:

        >>> gts = Genotypes.single('A', Genotype.HOMOZYGOUS_ALTERNATE)

        >>> assert len(gts) == 1
        >>> assert gts.for_sample('A') == Genotype.HOMOZYGOUS_ALTERNATE
        """

        return Genotypes((sample_id,), (genotype,))

    @staticmethod
    def from_mapping(mapping: typing.Mapping[str, Genotype]):
        """
        Create `Genotypes` from mapping between sample IDs and genotypes.

        >>> gts = Genotypes.from_mapping({'A': Genotype.HETEROZYGOUS, 'B': Genotype.HOMOZYGOUS_ALTERNATE})

        >>> assert len(gts) == 2
        """

        return Genotypes(*Genotypes._preprocess_mapping(mapping))

    @staticmethod
    def _preprocess_mapping(genotypes: typing.Mapping[str, Genotype]) -> typing.Tuple[typing.Sequence[str], typing.Sequence[Genotype]]:
        samples = np.empty(shape=(len(genotypes),), dtype=object)
        gts = np.empty(shape=(len(genotypes),), dtype=object)

        for i, (sample, gt) in enumerate(genotypes.items()):
            samples[i] = sample
            gts[i] = gt

        indices = np.argsort(samples)

        return samples[indices], gts[indices]

    def __init__(self, samples: typing.Iterable[str], genotypes: typing.Iterable[Genotype]):
        self._samples = tuple(samples)
        self._gts = tuple(genotypes)
        if len(self._samples) != len(self._gts):
            raise ValueError(f'Mismatch between the sample and genotype count: {len(self._samples)} != {len(self._gts)}')

    def for_sample(self, sample_id: str) -> typing.Optional[Genotype]:
        """
        Get a genotype for a sample or `None` if the genotype is not present.

        :param sample_id: a `str` with sample's identifier.
        """
        idx = bisect.bisect_left(self._samples, sample_id)
        if idx != len(self._samples) and self._samples[idx] == sample_id:
            return self._gts[idx]
        return None

    def __len__(self) -> int:
        return len(self._samples)

    def __iter__(self) -> Iterator[typing.Tuple[str, Genotype]]:
        return zip(self._samples, self._gts)

    def __hash__(self):
        return hash((self._samples, self._gts))

    def __eq__(self, other):
        return (isinstance(other, Genotypes)
                and self._samples == other._samples
                and self._gts == other._gts)

    def __str__(self):
        return f'Genotypes({["{}={}".format(sid, gt) for sid, gt in zip(self._samples, self._gts)]})'

    def __repr__(self):
        return f'Genotypes(samples={self._samples}, genotypes={self._gts})'


EMPTY = Genotypes((), ())


class Genotyped(metaclass=abc.ABCMeta):
    """
    `Genotyped` entities
    """

    @property
    @abc.abstractmethod
    def genotypes(self) -> Genotypes:
        pass

    def genotype_for_sample(self, sample_id: str) -> typing.Optional[Genotype]:
        """
        Get a genotype for a sample or `None` if the genotype is not present.

        :param sample_id: a `str` with sample's identifier.
        """
        return self.genotypes.for_sample(sample_id)
