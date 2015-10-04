#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
Slater Module
-------------

This module implement several classes in order to deal with atomic orbital,
electronic configuration and compute atomic orbital energies using the slater's
rule.

Important note: In the following, electronic sub-shell such as 2s, 3d ... are
called atomic orbital (AO) for simplicity. However, 2p, for example, is an
electronic sub-shell which contains 3 AO 2p_0, 2p_-1, 2p_1.

The implementation of this module is strongly inspired by the 
electronic_structure.core module of [pymatgen](http://pymatgen.org/) package.
"""

__author__ = "Germain Salvato Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"
__date__ = "Juin 2013"

import collections
import copy

class _AOTypeImpl(collections.namedtuple("_AOType", "name l")):
    """ implement an electronic shell type. Do not use directly. """

    __slots__ = ()

    def __int__(self):
        return self.l

    def __eq__(self, other):
        if other is None:
            return False
        else:
            return self.l == other.l

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

class AOType(object):
    """ Define electronic shell types """

    s = _AOTypeImpl("s", 0)
    p = _AOTypeImpl("p", 1)
    d = _AOTypeImpl("d", 2)
    f = _AOTypeImpl("f", 3)

    all_shells = (s, p, d, f)

    @staticmethod
    def from_string(st):
        """ 
        return an orbital from a string representation. Example "s", "p" ... 
        """
        for ao in AOType.all_shells:
            if str(ao) == st:
                return ao
        raise ValueError("Illegal electronic shell definition")

    @staticmethod
    def from_int(l):
        """ 
        return an orbital from an int. Example, 0, 1, 2 ...
        """
        if l > 3:
            raise ValueError("l must be lower than 3.")

        for ao in AOType.all_shells:
            if ao.l == l:
                return ao

class AO(object):
    """ 
    Define an AO (electronic subshell exactly) characterized by quantum numbers
    n and l.
    """

    def __init__(self, n, aoType, occ=0):
        try:
            self._n = int(n)
        except ValueError:
            raise ValueError("n must be an integer")

        if isinstance(aoType, str):
            self._AO = AOType.from_string(aoType)
        elif isinstance(aoType, _AOTypeImpl):
            self._AO = aoType
        elif isinstance(aoType, int):
            self._AO = AOType.from_int(aoType)
        else:
            raise TypeError("aoType must be either a string or an int or a AOType object.")

        if self._AO.l > self._n - 1:
            raise ValueError("n and l are unconsistant")

        try:
            occ = int(occ)
        except ValueError:
            raise ValueError("occ must be an integer")

        if occ > 4 * self._AO.l + 2:
            raise ValueError("occupation is too be in respect with this type of AO")
        self._occ = occ
        
    @staticmethod
    def from_string(st):
        """ 
        return an atomic orbtal object from a string such as '2s' or '2p^6'
        """

        st = st.strip()
        try:
            n = int(st[0])
        except ValueError:
            raise ValueError("syntax must be such as '2s' or '2p^6'")
        aotype = AOType.from_string(st[1])


        if len(st) == 4:
            try:
                occ = int(st[3])
            except ValueError:
                raise ValueError("Syntax must be such as '2s' or '2p^6'")
        else:
            occ = 0

        return AO(n, aotype, occ)

    def _get_occ(self):
        return self._occ

    def _set_occ(self, val):
        try:
            val = int(val)
        except ValueError:
            raise ValueError("Occuapation must be an integer")
        if val > 2 * self.degeneracy:
            raise ValueError("occupation is too be in respect with this type of AO")
        self._occ = val
    occ = property(_get_occ, _set_occ)
    """ occupation of the atomic orbital. """

    @property
    def n(self):
        """ main quantum number """
        return self._n

    @property
    def l(self):
        """ secondary quantum number """
        return self._AO.l

    @property
    def typ(self):
        """ AO type """
        return self._AO

    @property
    def name(self):
        """ AO name """
        return str(self._n) + self._AO.name

    @property
    def degeneracy(self):
        """ number of degenerated atomic orbital. """
        return 2 * self._AO.l + 1


    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def __eq__(self, other):
        return self.n == other.n and self.l == other.l
        
    def __hash__(self):
        return hash((self.n, self.l))

class Klechkowski(tuple):
    """ This class gives klechkowski rule """
    
    def __new__(cls):
        """ class constructor """
        k = (AO(1, AOType.s),
            AO(2, AOType.s), AO(2, AOType.p),
            AO(3, AOType.s), AO(3, AOType.p),
            AO(4, AOType.s), AO(3, AOType.d), AO(4, AOType.p),
            AO(5, AOType.s), AO(4, AOType.d), AO(5, AOType.p),
            AO(6, AOType.s), AO(4, AOType.f), AO(5, AOType.d), AO(6, AOType.p),
            AO(7, AOType.s), AO(5, AOType.f), AO(6, AOType.d), AO(7, AOType.p))
        return super(Klechkowski, cls).__new__(cls, k)

    def __str__(self):
        s = ""
        n = 1
        for ao in self:
            if ao.n > n:
                s += "\n"
                n = ao.n
            s += " " + str(ao)
        return s

class ElectronicConf(list):
    """ 
    This class describes an electronic configuration of an atom.
    
    WARNING: You must start with an atom and then build an ion.
    """

    klechkowski = Klechkowski()

    def __init__(self, *args, **kargs):
        """ initialize the object as a neutral atom """
        if len(args) != 0:
            raise TypeError("You must use a key=value syntax")
        if len(kargs) != 1:
            raise TypeError("__init__() takes exactly 1 key=value argument")
            
        if "nelec" in kargs:
            try:
                nelec = int(kargs["nelec"])
                if nelec > 118:
                    raise ValueError("nelec > 118")
            except ValueError:
                raise ValueError("nelec must be an integer")
            list.__init__(self, ElectronicConf.__buildConfig(nelec))
            self._z = nelec
        elif "aolist" in kargs:
            self._z = sum([ao.occ for ao in kargs["aolist"]])
            list.__init__(self, kargs["aolist"])
        else:
            raise ValueError("__init__() argument unknown")

    @property
    def nelec(self):
        """ number of electron """
        return sum([ao.occ for ao in self])

    @property
    def Z(self):
        """ atomic number """
        return self._z

    @property
    def q(self):
        """ atomic number """
        return self._z - self.nelec

    @property
    def maxn(self):
        """ maximum value of n quantum number """
        maxn = -1
        for ao in self:
            maxn = max(maxn, ao.n)
        return maxn

    @property
    def valence(self):
        """ valence electronic configuration """
        val = list()
        s = ""
        for ao in self:
            if ao.n == self.maxn or ao.occ != 2 * ao.degeneracy:
                val.append(ao)
                s += "%s^%d " % (ao.name, ao.occ)
        #return tuple([s[:-1], tuple(val)])
        return s.strip()
        
    @property
    def energy(self):
        """ 
        Compute the energy associated to the electronic configuration using the
        Slater's rule
        """
        data = self.computeEnergy(verb=False)
        e = 0
        for oa in self:
            e += oa.occ * data[oa][1]
        return e

    @staticmethod
    def __buildConfig(nelec):
        """ build electronic configuration """

        ne = 0
        ishell = 0
        while ne < nelec:
            ne += 2 * ElectronicConf.klechkowski[ishell].degeneracy
            ishell += 1

        config = list()
        reste = nelec
        for ao in ElectronicConf.klechkowski[:ishell]:
            config.append(copy.deepcopy(ao))
            config[-1].occ = min(2 * ao.degeneracy, reste)
            reste -= config[-1].occ

        return config

    @staticmethod
    def from_string(s):
        """ 
        Return an electronic configuration object form a string such as :

        1s^2 2s^2 2p^2
        """
        config = list()
        for el in s.split():
            ao, occ = el.split("^")
            n = int(ao[0])
            stype = ao[1]

            config.append(AO(n, AOType.from_string(stype), int(occ)))
        
        return ElectronicConf(aolist=config)

    def computeEnergy(self, verb=True):
        """ compute energy of all atomic orbitals """
        
        def eneOA(z, n):
            return -13.602 * (z / n)**2

        dat = dict()

        if verb:
            print(35 * "-")
            print("# AO   sigma     Z*    eps (eV)")
            print(35 * "-")
        for aoi in self:
            if aoi.name == "1s":
                sigma = 0.31 * (aoi.occ - 1)
            else:
                sigma = 0.35 * (aoi.occ - 1)

            for aoj in self:

                if aoj.n == aoi.n and aoj.l == aoi.l:
                    continue
                elif aoi.typ == AOType.s or aoi.typ == AOType.p:
                    if aoj.n < aoi.n - 1:
                        sigma += 1. * aoj.occ
                    elif aoj.n == aoi.n - 1:
                        sigma += .85 * aoj.occ
                    elif aoj.n == aoi.n and aoj.l < 2:
                        sigma += 0.35 * aoj.occ
                else:
                    if aoj.n <= aoi.n:
                        sigma += 1. * aoj.occ

            zstar = self._z - sigma
            e = eneOA(zstar, aoi.n)
            if verb:
                print("%4s%8.2f%8.2f%10.2f" % (aoi, sigma, zstar, e))

            dat[aoi] = (sigma, e)
        if verb:
            print(35 * "-")

        return dat

    def ionize(self, q):
        """ 
        Remove or add q electrons from the electronic configuration.

        If q > 0, abs(q) electrons are removed.
        If q < 0, abs(q) electrons are added.
        """

        def maxl(config, nval):
            """ return the max value of l for n=nval """
            lmax = -1
            for ao in config:
                if ao.n == nval:
                    lmax = max(lmax, ao.l)
            return lmax

        def adde(ao, reste):
            vac = 2 * ao.degeneracy - ao.occ
            if vac > 0:
                if reste <= vac:
                    ao.occ += reste
                    reste = 0
                else:
                    reste -= vac
                    ao.occ += vac
            return ao, reste

        ion = copy.deepcopy(self)

        if q > 0:
            # remove electrons
            reste = q
            while reste != 0:
                lmax = maxl(ion, ion.maxn)
                for ao in ion:
                    if ao.n == ion.maxn and ao.l == lmax:
                        if ao.occ > reste:
                            ao.occ -= reste
                            reste = 0
                        elif ao.occ == reste:
                            reste = 0
                            ion.remove(ao)
                        else:
                            reste -= ao.occ
                            ion.remove(ao)
        else:
            # add electrons following klechkowski rule
            reste = -q
            i = 0
            while reste != 0:
                ao = ElectronicConf.klechkowski[i]
                if ao in ion:
                    iao = ion.index(ao)
                    ion[iao], reste = adde(ion[iao], reste)
                else:
                    ion.insert(i, ao)
                    ion[i], reste = adde(ion[i], reste)

                i += 1
                
        return ion

    def toTex(self):
        """ return a LaTeX code of the electronic configuration """
        s = r"$"
        s += r"".join([r"\text{%s}^\text{%d}" % (ao.name, ao.occ) for ao in self])
        s += r"$"
        return s

    def __str__(self):
        return " ".join(["%s^%d" % (ao.name, ao.occ) for ao in self])

    def __repr__(self):
        return self.__str__()
        
class STO(AO):
    """
    This class define a Slater type orbital and rely on the atomic orbital class.
    """

    def __init__(self, n, aoType, occ=0):
        """
        Set up an atomic orbital from n and l value or OA type.
        """
        AO.__init__(self, n, aoType, occ)
    
    
    
    
    
