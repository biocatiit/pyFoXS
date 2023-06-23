"""
\file IMP/core/XYZ.h     \brief Simple XYZ decorator.

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

from IMP.algebra import Vector3D, Transformation3D
from IMP import core

class XYZ(core.Decorator):
    def __init__(self, particle: core.Particle):
        super().__init__(particle)

    @staticmethod
    def get_coordinate_key(i: int):
        assert 0 <= i < 3, "Out of range coordinate"
        return core.internal.xyzr_keys[i]

    @property
    def x(self) -> float:
        return self.get_particle().get_value(self.get_coordinate_key(0))

    @x.setter
    def x(self, value: float):
        self.get_particle().set_value(self.get_coordinate_key(0), value)

    @property
    def y(self) -> float:
        return self.get_particle().get_value(self.get_coordinate_key(1))

    @y.setter
    def y(self, value: float):
        self.get_particle().set_value(self.get_coordinate_key(1), value)

    @property
    def z(self) -> float:
        return self.get_particle().get_value(self.get_coordinate_key(2))

    @z.setter
    def z(self, value: float):
        self.get_particle().set_value(self.get_coordinate_key(2), value)

    def set_coordinate(self, i: int, value: float):
        assert 0 <= i < 3, "Out of range coordinate"
        self.get_particle().set_value(self.get_coordinate_key(i), value)

    def set_coordinates(self, v: Vector3D):
        self.x, self.y, self.z = v

    def get_coordinate(self, i: int) -> float:
        assert 0 <= i < 3, "Out of range coordinate"
        return self.get_particle().get_value(self.get_coordinate_key(i))

    def get_derivative(self, i: int) -> float:
        return self.get_particle().get_derivative(self.get_coordinate_key(i))

    def add_to_derivative(self, i: int, v: float, d):
        self.get_particle().add_to_derivative(self.get_coordinate_key(i), v, d)

    def add_to_derivatives(self, v: Vector3D, d):
        self.get_model().add_to_coordinate_derivatives(self.get_particle_index(), v, d)

    def get_coordinates(self) -> Vector3D:
        return Vector3D(self.x, self.y, self.z)

    def get_derivatives(self) -> Vector3D:
        return self.get_model().get_coordinate_derivatives(self.get_particle_index())

    def get_coordinates_are_optimized(self) -> bool:
        return all(
            self.get_particle().get_is_optimized(self.get_coordinate_key(i))
            for i in range(3)
        )

    def set_coordinates_are_optimized(self, tf: bool):
        for i in range(3):
            self.get_particle().set_is_optimized(self.get_coordinate_key(i), tf)

    def get_vector_to(self, b: 'XYZ') -> Vector3D:
        return b.get_coordinates() - self.get_coordinates()

    def show(self) -> str:
        return f"({', '.join(map(str, self.get_coordinates()))})"

def transform(a: XYZ, tr: Transformation3D):
    assert not isinstance(a, core.RigidBody), (
        "Particle is also a RigidBody, so this function would do the wrong "
        "thing; decorate your Particle as a RigidBody instead and call "
        "transform(RigidBody a, Transformation3D tr)"
    )
    a.set_coordinates(tr.get_transformed(a.get_coordinates()))

def get_dihedral(a: XYZ, b: XYZ, c: XYZ, d: XYZ) -> float:
    return core.internal.dihedral(
        a.get_coordinates(), b.get_coordinates(),
        c.get_coordinates(), d.get_coordinates()
    )

def set_vector_geometry(d: XYZ, v: Vector3D):
    d.set_coordinates(v)

def get_vector_geometry(d: XYZ) -> Vector3D:
    return d.get_coordinates()

def get_distance(a: XYZ, b: XYZ) -> float:
    return (a.get_coordinates() - b.get_coordinates()).get_magnitude()


class XYZR(XYZ):
    @staticmethod
    def do_setup_particle(m, pi, s):
        XYZ.setup_particle(m, pi, s.get_center())
        XYZR.do_setup_particle(m, pi, s.get_radius())

    @staticmethod
    def do_setup_particle(m, pi, r):
        m.add_attribute(XYZR.get_radius_key(), pi, r, False)

    @staticmethod
    def do_setup_particle(m, pi):
        if not XYZ.get_is_setup(m, pi):
            XYZ.setup_particle(m, pi)
        m.add_attribute(XYZR.get_radius_key(), pi, 0, False)

    def get_radius(self):
        return self.get_sphere().get_radius()

    def set_radius(self, r):
        self.get_model().get_sphere(self.get_particle_index())[3] = r

    def get_sphere(self):
        return self.get_model().get_sphere(self.get_particle_index())

    @staticmethod
    def get_radius_key():
        return core.internal.xyzr_keys[3]

    def add_to_radius_derivative(self, v, d):
        self.get_particle().add_to_derivative(XYZR.get_radius_key(), v, d)
