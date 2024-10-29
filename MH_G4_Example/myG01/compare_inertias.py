#!/usr/bin/python
import math
import sys

from dataclasses import dataclass

PREC = 0.01  # 1%

@dataclass
class Vector:
    x: float
    y: float
    z: float

    def __str__(self):
        return f'({self.x}, {self.y}, {self.z})'

    def mag(self):
        return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

    def __eq__(self, other):
        if isinstance(other, Vector):
            minx = min(abs(self.x), abs(other.x))
            if abs(self.x - other.x) > 0.01:
                return False
            miny = min(abs(self.y), abs(other.y))
            if abs(self.y - other.y) > 0.01:
                return False
            minz = min(abs(self.z), abs(other.z))
            if abs(self.z - other.z) > 0.01:
                return False
            return True
        else:
            return False


@dataclass
class InertiaMatrix:
    Ixx: float
    Iyy: float
    Izz: float

    def __str__(self):
        return f'({self.Ixx}, {self.Iyy}, {self.Izz})'

    def __eq__(self, other):
        if isinstance(other, InertiaMatrix):
            v1 = Vector(self.Ixx, self.Iyy, self.Izz)
            v2 = Vector(other.Ixx, other.Iyy, other.Izz)
            return v1 == v2
        else:
            return False


@dataclass
class Moments:
    Volume: float
    volUnit: str
    cm: Vector
    cmUnit: str
    II: InertiaMatrix
    inertiaUnit: str

    def __eq__(self, other):
        if isinstance(other, Moments):
            minvol = min(abs(self.Volume), abs(other.Volume))
            if abs(self.Volume - other.Volume) > PREC*minvol:
                return False
            if self.volUnit != other.volUnit:
                return False
            if self.cm != other.cm:
                return False
            if self.cmUnit != other.cmUnit:
                return False
            if self.II != other.II:
                return False
            if self.inertiaUnit != other.inertiaUnit:
                return False
            return True
        else:
            return False
            

def read_moments(filename) -> Moments | None:
        
    fd = open(filename, 'r')

    for line in fd.readlines():
        if "Total" in line and "Volume" in line:
            line = line.replace('(',' ')
            line = line.replace(',',' ')
            line = line.replace(')',' ')
            words = line.strip().split()
            volume = float(words[2])
            volUnit = words[3]
            cm = Vector(float(words[8]),
                        float(words[9]),
                        float(words[10]))
            cmUnit = words[12]
            II = InertiaMatrix(float(words[15]),
                               float(words[16]),
                               float(words[17]))
            inertiaUnit = words[18]
            if 'x10' in inertiaUnit:
                inertiaUnit += words[19]

            return Moments(volume, volUnit, cm, cmUnit, II, inertiaUnit)

    return None
                
            

    
file1 = sys.argv[1]
file2 = sys.argv[2]

moments1 = read_moments(file1)
moments2 = read_moments(file2)

if moments1 != moments2:
    print(f'{file1} {file2} have unequal moments:')
    print(f'file1: vol: {moments1.Volume} {moments1.volUnit}  cm: {moments1.cm} {moments1.cmUnit} Inertia: {moments1.II} {moments1.inertiaUnit}')
    print(f'file2: vol: {moments2.Volume} {moments2.volUnit}  cm: {moments2.cm} {moments2.cmUnit} Inertia: {moments2.II} {moments2.inertiaUnit}')



