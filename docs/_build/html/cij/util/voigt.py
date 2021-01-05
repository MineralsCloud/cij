'''
Hashable representation of standard and voigt subscripts.
'''

from typing import NamedTuple
from enum import Enum, auto


VOIGT_TO_STANDARD = {
    1: (1, 1),
    2: (2, 2),
    3: (3, 3),

    4: (2, 3),
    5: (1, 3),
    6: (1, 2)
}


STANDARD_TO_VOIGT = dict((v, k) for k, v in VOIGT_TO_STANDARD.items())


class ElasticModulusCalculationType(Enum):
    LONGITUDINAL = auto()
    OFF_DIAGONAL = auto()
    SHEAR = auto()


class StrainRepresentation(NamedTuple):

    i: int
    j: int

    @classmethod
    def create(cls, i: int, j: int = None):
        if j is not None and i is not None:
            return cls.from_standard(i, j)
        elif j is None and type(i) == int:
            if i < 10:
                return cls.from_voigt(i)
            else:
                return cls.create(str(i))
        elif j is None and type(i) == str:
            return cls.create(*(int(k) for k in i))
        else:
            raise RuntimeError(f"Invalid indices ({i}{j})")

    @classmethod
    def from_voigt(cls, i):
        if i not in VOIGT_TO_STANDARD.keys():
            raise RuntimeError(f"Invalid voigt index {i}")
        return cls(*VOIGT_TO_STANDARD[i])

    @classmethod
    def from_standard(cls, i, j):
        i, j = sorted((i, j))
        if (i, j) not in VOIGT_TO_STANDARD.values():
            raise RuntimeError(f"Invalid standard index {i}{j}")
        return cls(i, j)

    @property
    def voigt(self):
        return STANDARD_TO_VOIGT[self]

    @property
    def v(self):
        return self.voigt

    @property
    def standard(self):
        return tuple(self)

    @property
    def s(self):
        return self.standard

    def __repr__(self):
        return "%d(%s)" % (self.v, "".join(str(i) for i in self.s))

    @classmethod
    def _(cls, *args):
        return cls.create(*args)
        

class ModulusRepresentation(NamedTuple):

    i: StrainRepresentation
    j: StrainRepresentation

    @classmethod
    def create(cls, *args):
        if len(args) == 4:
            return cls.from_standard(*args)
        elif len(args) == 2:
            return cls.from_voigt(*args)
        elif len(args) == 1 and type(args[0]) == str:
            return cls.create(*(int(c) for c in args[0]))
        elif len(args) == 1 and type(args[0]) == int:
            return cls.create(str(args[0]))
        else:
            raise RuntimeError(f"Invalid modulus representation {args}")

    @classmethod
    def from_voigt(cls, i, j):
        return cls(
            *sorted((
                StrainRepresentation.from_voigt(i),
                StrainRepresentation.from_voigt(j)
            ), key=lambda e: e.voigt)
        )

    @classmethod
    def from_standard(cls, i, j, k, l):
        return cls(
            *sorted((
                StrainRepresentation.from_standard(i, j),
                StrainRepresentation.from_standard(k, l),
            ), key=lambda e: e.voigt)
        )

    @property
    def voigt(self):
        return (self.i.voigt, self.j.voigt)

    @property
    def v(self):
        return self.voigt

    @property
    def standard(self):
        return (*self.i, *self.j)

    @property
    def s(self):
        return self.standard

    def __repr__(self):
        return "%s(%s)" % (
            "".join(str(i) for i in self.v),
            "".join(str(i) for i in self.s)
        )
    
    @classmethod
    def _(cls, *args):
        return cls.create(*args)

    @property
    def is_longitudinal(self) -> bool:
        return self.i == self.j and not self.is_shear

    @property
    def is_off_diagonal(self) -> bool:
        return not self.is_shear and not self.is_longitudinal

    @property
    def is_shear(self) -> bool:
        return self.i.voigt in {4, 5, 6} or self.j.voigt in {4, 5, 6}

    @property
    def multiplicity(self) -> int:
        return 1 \
            << (self.i != self.j) \
            << (self.i.i != self.i.j) \
            << (self.j.i != self.j.j)
    
    @property
    def calc_type(self) -> ElasticModulusCalculationType:
        if self.is_longitudinal:
            return ElasticModulusCalculationType.LONGITUDINAL
        if self.is_off_diagonal:
            return ElasticModulusCalculationType.OFF_DIAGONAL
        if self.is_shear:
            return ElasticModulusCalculationType.SHEAR

E_ = StrainRepresentation
C_ = ModulusRepresentation

