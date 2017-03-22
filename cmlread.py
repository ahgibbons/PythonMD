## Python 2.7

import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import copy

testTreePath = "/media/andrewsim/Samsung_T1/pythonMD/cyclohexane.cml"
testTree = ET.parse(testTreePath)
testRoot = testTree.getroot()

atomsArray= testRoot.getchildren()[0]

atoms = atomsArray.getchildren()

class Molecule:
    def __init__(self, cmlpath):
        cmlTree = ET.parse(cmlpath)
        molecule = cmlTree.getroot()
        [atomArray, bondArray] = molecule.getchildren()
        atoms = atomArray.getchildren()
        bonds = bondArray.getchildren()
        self.atoms = [Atom(ae) for ae in atoms]
        self.bonds = [Bond(be) for be in bonds]

    def total_mass(self):
        return sum([a.mass for a in self.atoms])
        
    def centre_of_mass(self):
        return (sum([a.pos*a.mass for a in self.atoms]) / self.total_mass())

    def radius_of_gyration_all_atoms(self):
        r_mean = sum([a.pos for a in self.atoms]) / len(self.atoms)
        r_g = sum([(a.pos - r_mean)**2 for a in self.atoms]) / len(self.atoms)
        return r_g

    def radius_of_gyration_carbon(self):
        r_mean = sum([a.pos for a in self.atoms if a.element_type=='C']) / len(self.atoms)
        r_g = sum([(a.pos - r_mean)**2 for a in self.atoms if a.element_type=='C']) / len(self.atoms)

    def monomer_bonds(self):
        c_bond = element_bonds['C']
        carbon_ids = [a.id for a in self.atoms if a.element_type=='C']
        candidate_links = [aid for aid in carbon_ids if self.atom_bonds(aid)!=4]
        if len(candidate_links)!=2:
            raise Exception("Invalid monomer configuration. Cannot determine linking points")
        return candidate_links
        
    def atom_bonds(self,atom_id):
        return sum([b.order for b in self.bonds if atom_id in b.atomRefs])

    def translate(self,dr):
        for atom in self.atoms:
            atom.translate(dr)

    def bond(self, molecule, sid, oid, order=1):
        curids = [int(aid[1:]) for aid in molecule.atoms]
        maxid = curids.sort(reverse=True)[0]
    
def genPolymer(monomer, number, sep=1):
    m 
    

class Atom:
    def __init__(self,atomElement):
        attribs = atomElement.attrib

        self.element_type = attribs['elementType']
        self.id = attribs['id']
        x = float(attribs['x3'])
        y = float(attribs['y3'])
        z = float(attribs['z3'])
        self.pos = np.array([x,y,z])
        self.mass = element_mass[self.element_type]

    def translate(self,dr):
        self.pos+=dr

class Bond:
    def __init__(self,bondElement):
        attribs = bondElement.attrib

        self.atomRefs = attribs['atomRefs2'].split()
        self.order = int(attribs['order'])

element_mass = {'H':1, 'C':12}
element_bonds = {'C':4}
