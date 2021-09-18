#!/usr/bin/python3


"""
MDH.compound
~~~~~~~~~~~~~~~~~~~~~~~~~

This module provides the :class:`~MDH.compound.Atom` class, the :class:`MDH.compound.Bond` class,
and the :class:`~MDH.compound.Compound` class to construct compound entity.

"""

import collections
from math import atan, sin, cos
import itertools
import numpy
import BASS_Cc_3 as BASS_1
import heapq


class Atom:

	""" Atom class describes the :class:`~MDH.compound.Atom` entity in the compound. """

    def __init__(self, x, y, z, atom_symbol, mass_difference, charge, atom_stereo_parity, hydrogen_count, stereo_care_box, 
    	valence, h0designator, atom_atom_mapping_number, inversion_retention_flag, atom_number, exact_change_flag, default_symbol, kat=""):
    	
    	"""Atom initializer.

    	:param str x: the x coordinate 
    	:param str y:
    	:param str z:
    	:param str atom_symbol:
    	:param str mass_difference:
    	:param str charge:
    	:param str atom_stereo_parity:
    	:param str hydrogen_count:
    	:param str stereo_care_box
    	:param str valence
    	:param str h0designator
    	:param str atom_atom_mapping_number
    	:param str inversion_retention_flag
    	:param str atom_number
    	:param str exact_change_flag
    	:param str default_symbol
    	:param str kat

    	"""

        self.x = float(x.strip())
        self.y = float(y.strip())
        self.z = float(z.strip())
        self.atom_symbol = atom_symbol.strip()
        self.mass_difference = mass_difference.strip()
        self.charge = int(charge.strip())
        self.atom_stereo_parity = atom_stereo_parity.strip()
        self.hydrogen_count = hydrogen_count.strip()
        self.stereo_care_box = stereo_care_box.strip()
        self.valence = valence.strip()
        self.h0designator = h0designator.strip()
        self.atom_atom_mapping_number = atom_atom_mapping_number.strip()
        self.inversion_retention_flag = inversion_retention_flag.strip()
        self.atom_number = atom_number
        self.exact_change_flag = exact_change_flag.strip()
        self.default_symbol = default_symbol.strip()
        self.neighbors = []
        self.color_0 = ""
        self.color = ""
        self.in_cycle = False
        self.bond_counts = 0
        self.group_id = 0
        self.double_bond_counts = 0
        self.color_tuple = tuple()
        self.color_layers = {}
        self.distance_to_R = 0
        self.kat = kat


    def update_symbol(self, symbol):
        
        self.atom_symbol = symbol


    def remove_neighbors(self, neighbors):
        
        for atom_index in neighbors:
            if atom_index in self.neighbors:
                self.neighbors.remove(atom_index)


    def add_neighbors(self, neighbors):
        
        for atom_index in neighbors:
            if atom_index not in self.neighbors:
                self.neighbors.append(atom_index)


    def update_stereochemistry(self, stereo):

        self.atom_stereo_parity = stereo
    

    def is_R(self):

        if "A" in self.default_symbol or "R" in self.default_symbol or "*" in self.default_symbol:
            if self.default_symbol not in notR:
                return True
        return False


class Bond:

	""" Bond class describes the :class:`~MDH.compound.Bond` entity in the compound. """

    def __init__(self, first_atom_number, second_atom_number, bond_type, bond_stereo, bond_topology, reacting_center_status):

    	"""Bond initializer.



    	"""

        self.first_atom_number = int(first_atom_number) - 1
        self.second_atom_number = int(second_atom_number) - 1
        self.bond_type = bond_type.strip()
        self.bond_stereo = bond_stereo.strip()
        self.bond_topology = bond_topology.strip()
        self.reacting_center_status = reacting_center_status.strip()

    
    def update_bond_type(self, bond_type):

        self.bond_type = bond_type

    def update_stereochemistry(self, stereo):

        self.bond_stereo = stereo


class Compound:

	""" Compound class describes the :class:`~MDH.compound.Compound` entity. """

    def __init__(self, entry, atoms, bonds):

    	"""Compound initializer.
			

    	"""

        self.atoms = atoms
        self.bonds = bonds
        self.entry = entry
        self.color_groups = collections.defaultdict(list)
        self.bond_lookup = {}
        self.has_cycle = False
        
        for bond in self.bonds:
            if bond.bond_type != "8":
                firstAtom, secondAtom = self.atoms[bond.first_atom_number], self.atoms[bond.second_atom_number]
                self.bond_lookup[(bond.first_atom_number, bond.second_atom_number)] = bond
                self.bond_lookup[(bond.second_atom_number, bond.first_atom_number)] = bond
                firstAtom.neighbors.append(bond.second_atom_number)
                firstAtom.bond_counts += int(bond.bond_type)
                secondAtom.neighbors.append(bond.first_atom_number)
                secondAtom.bond_counts += int(bond.bond_type)     
        
        self.find_cycles()
        self.calculate_distance_to_R_groups()
    
    @property
    def R_groups(self):
        
        rs = []
        for index, atom in enumerate(self.atoms):
            if atom.is_R():
                rs.append(index)
        return rs

    def calculate_distance_to_R_groups(self):
        
        """

        """
        distance_matrix = [len(self.heavy_atoms)] * len(self.heavy_atoms) 
        
        if self.R_groups:
            for r_index in self.R_groups:
                ques = [(0, r_index)]
                seen = {}
                while ques:
                    l, atom_index = heapq.heappop(ques)
                    if atom_index in seen and seen[atom_index] <= l:
                        continue
                    seen[atom_index] = l
                    distance_matrix[self.index_of_heavy_atoms[atom_index]] = min(distance_matrix[self.index_of_heavy_atoms[atom_index]], l)
                    for neighbor in self.atoms[atom_index].neighbors:
                        if self.atoms[neighbor].atom_symbol != "H":
                            if neighbor not in seen or seen[neighbor] > l + 1:
                                if distance_matrix[self.index_of_heavy_atoms[neighbor]] > l+1:
                                    heapq.heappush(ques, (l+1, neighbor))
        
            for i, l in enumerate(distance_matrix):
                self.heavy_atoms[i].distance_to_R = l


    @property
    def heavy_atoms(self):
       
        return [atom for atom in self.atoms if atom.atom_symbol != "H"]
    
    @property
    def index_of_heavy_atoms(self):

        return { self.heavy_atoms[i].atom_number: i for i in range(len(self.heavy_atoms)) } 

    @property
    def loose_matrix(self):

        loose_matrix = numpy.zeros((len(self.heavy_atoms), len(self.heavy_atoms)), dtype=numpy.uint8)
        for bond in self.bonds:
            atom_1, atom_2 = bond.first_atom_number, bond.second_atom_number
            if atom_1 in self.index_of_heavy_atoms and atom_2 in self.index_of_heavy_atoms:
                loose_matrix[self.index_of_heavy_atoms[atom_1]][self.index_of_heavy_atoms[atom_2]] = 
                loose_matrix[self.index_of_heavy_atoms[atom_2]][self.index_of_heavy_atoms[atom_1]] = int(bond.bond_type) if int(bond.bond_type) != 2 else 1
        return loose_matrix

    @property
    def strict_matrix(self):
        matrix = numpy.zeros((len(self.heavy_atoms), len(self.heavy_atoms)), dtype=numpy.uint8)
        for bond in self.bonds:
            atom_1, atom_2 = bond.first_atom_number, bond.second_atom_number
            if atom_1 in self.index_of_heavy_atoms and atom_2 in self.index_of_heavy_atoms:
                matrix[self.index_of_heavy_atoms[atom_1]][self.index_of_heavy_atoms[atom_2]] = 
                matrix[self.index_of_heavy_atoms[atom_2]][self.index_of_heavy_atoms[atom_1]] = int(bond.bond_type)
        return matrix


    @property
    def distance_matrix(self):
        if self.heavy_atoms:
            dis_matrix = numpy.ones((len(self.heavy_atoms), len(self.heavy_atoms)),dtype=numpy.uint16 ) * len(self.heavy_atoms)
            for bond in self.bonds:
                atom_1, atom_2 = bond.first_atom_number, bond.second_atom_number
                if atom_1 in self.index_of_heavy_atoms and atom_2 in self.index_of_heavy_atoms:
                    dis_matrix[self.index_of_heavy_atoms[atom_1]][self.index_of_heavy_atoms[atom_2]] = 
                    dis_matrix[self.index_of_heavy_atoms[atom_2]][self.index_of_heavy_atoms[atom_1]] = 1

            n = len(self.heavy_atoms)
            for k in range(n):
                for i in range(n):
                    for j in range(n):
                        if dis_matrix[i][j] > dis_matrix[i][k] + dis_matrix[k][j]:
                            dis_matrix[i][j] = dis_matrix[i][k] + dis_matrix[k][j]
            return dis_matrix
        return None

    def create_loose_color_tuple(self):

        for atom in self.heavy_atoms:
            components = {"C": 0, "N": 0, "O": 0, "S": 0}
            single_bond, aromatic_bond = 0, 0
            for neighbor_index in atom.neighbors:
                neighbor = self.atoms[neighbor_index]
                if neighbor.atom_symbol != "H":
                    if neighbor.atom_symbol in components:
                        components[neighbor.atom_symbol] += 1
                    bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                    if bond.bond_type == "4":
                        aromatic_bond += 1
                    if bond.bond_type in ["1", "2"]:
                        single_bond += 1
            atom.color_tuple = (components["C"], components["N"], components["S"], components["O"], single_bond, aromatic_bond)


    def colorCpd(self, R=True,stereo_resolved=True, atom_stereo=True, strict=True):

        self.firstRoundColor(R=R, stereo_resolved=stereo_resolved, atom_stereo=atom_stereo, strict=strict)
        self.curateSymmetry(stereo_resolved=stereo_resolved, strict=strict)
        self.colorMetals(stereo_resolved=stereo_resolved, strict=strict)
    
    def create_strict_color_tuple(self):

        for atom in self.heavy_atoms:
            components = {"C": 0, "N": 0, "O": 0, "S": 0}
            single_bond, double_bond, aromatic_bond = 0, 0, 0
            #print(atom.neighbors)
            for neighbor_index in atom.neighbors:
                neighbor = self.atoms[neighbor_index]
                if neighbor.atom_symbol != "H":
                    if neighbor.atom_symbol in components:
                        components[neighbor.atom_symbol] += 1
                    #print(atom.atom_number, neighbor_index)
                    bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                    if bond.bond_type == "4":
                        aromatic_bond += 1
                    if bond.bond_type == "1":
                        single_bond += 1
                    if bond.bond_type == "2":
                        double_bond += 1
            atom.color_tuple = (components["C"], components["N"], components["S"], components["O"], single_bond, double_bond, aromatic_bond)
    

    def update_mapped_atoms(self, another, mapped):
        
        for atom in self.atoms:
            atom.atom_atom_mapping_number = []
        for atom in another.atoms:
            atom.atom_atom_mapping_number = []

        if mapped:
            for atomI in mapped:
                for atomJ in mapped[atomI]:
                    self.atoms[atomI].atom_atom_mapping_number.append(atomJ)
                    another.atoms[atomJ].atom_atom_mapping_number.append(atomI)


    def find_strict_mappings(self, another, R=False, mapped=None):
        
        mappping_matrix = BASS_1.make_mapping_matrix_R if R else BASS_1.make_mapping_matrix 
        find_mappings = BASS_1.find_mappings
        
        self.create_strict_color_tuple()
        another.create_strict_color_tuple()
        substructures = []
        
        self.update_mapped_atoms(another, mapped)

        mmat = mappping_matrix(another, self, True, True)
        #print("strict mappings", sum(sum(row) for row in mmat))
        if mmat is not None:
            substructures = find_mappings(another.strict_matrix, another.distance_matrix, self.strict_matrix, self.distance_matrix, mmat)

        indexMatch = []

        for sub in substructures:
            thisMatch = {}
            for f, t in enumerate(sub):
                thisMatch[self.heavy_atoms[t].atom_number] = another.heavy_atoms[f].atom_number
            indexMatch.append(thisMatch)
        return indexMatch


    def find_loose_mappings(self, another, R=False):
        
        mappping_matrix = BASS_1.make_mapping_matrix_R if R else BASS_1.make_mapping_matrix
        find_mappings = BASS_1.find_mappings

        self.create_loose_color_tuple()
        another.create_loose_color_tuple()

        substructures = []
        mmat = mappping_matrix(another, self, True, True)
        #print("loose mappings", sum(sum(row) for row in mmat))
        if mmat is not None:
            substructures = find_mappings(another.loose_matrix, another.distance_matrix, self.loose_matrix, self.distance_matrix, mmat)
        
        indexMatch = []
        
        for sub in substructures:
            thisMatch = {}
            for f, t in enumerate(sub):
                thisMatch[self.heavy_atoms[t].atom_number] =  another.heavy_atoms[f].atom_number
            indexMatch.append(thisMatch)
        return indexMatch

    def match_loosely(self, another, R=False):

        indexMatch = self.find_loose_mappings(another, R=R)
        validMatch = []
        for thisMatch in indexMatch:
            flag = True
            reverseMatch = { thisMatch[key] : key for key in thisMatch }
            for f in sorted(thisMatch.keys()):
                t = thisMatch[f]
                if 1 in self.atoms[f].color_layers and 1 in another.atoms[t].color_layers and self.atoms[f].color_layers[1] == another.atoms[t].color_layers[1]:
                    continue
                
                three = {f}
                fNeighbor = self.find_double_bond_linked_atom(f)
                nNeighbor = another.find_double_bond_linked_atom(t)
                
                if fNeighbor == -1 and nNeighbor == -1:
                    flag = False
                    break

                elif fNeighbor == -1:
                    reverseFNeighbor = reverseMatch[nNeighbor]
                    three.add(reverseFNeighbor)
                    fNextNeighbor = self.find_double_bond_linked_atom(reverseFNeighbor)
                    if fNextNeighbor != -1:
                        three.add(fNextNeighbor)

                elif nNeighbor == -1 and fNeighbor in thisMatch:
                    three.add(fNeighbor)
                    reverseNNeighbor = thisMatch[fNeighbor]
                    nNextNeighbor = another.find_double_bond_linked_atom(reverseNNeighbor)
                    if nNextNeighbor != -1:
                        three.add(reverseMatch[nNextNeighbor])

                elif fNeighbor != -1 and nNeighbor != -1:
                    three.add(fNeighbor)
                    three.add(reverseMatch[nNeighbor])
                
                if len(three) < 3:
                    flag = False
                    break
                
                k = [self.atoms[i].atom_symbol for i in three].count("C")
                if k == 3:
                    flag = False
                    break
            if flag:
                validMatch.append(thisMatch)
        return []
    
    def find_double_bond_linked_atom(self, i):

        for neighbor_index in self.atoms[i].neighbors:
            if neighbor_index in self.index_of_heavy_atoms and self.bond_lookup[(i, neighbor_index)].bond_type == "2":
                return neighbor_index
        return -1

    def linear_circular_interchange(self, another):
        sAtoms = [atom.color_layers[1] if 1 in atom.color_layers else atom.atom_symbol for atom in self.heavy_atoms]
        aAtoms = [atom.color_layers[1] if 1 in atom.color_layers else atom.atom_symbol for atom in self.heavy_atoms]
        atomIntersect = list((collections.Counter(sAtoms) & collections.Counter(aAtoms)).elements())
        #print("atom name cnt: ", len(sAtoms), len(aAtoms), len(atomIntersect))
        if 2 <= len(sAtoms) - len(atomIntersect) <= 3:
            sinCycle = [ atom.in_cycle for atom in self.heavy_atoms].count(True)
            ainCycle = [ atom.in_cycle for atom in another.heavy_atoms].count(True)
            #print("incycle number: ", sinCycle, ainCycle)
            if sinCycle != ainCycle:
                return True
        return False

    def find_cycles(self, short_circuit=False, cutoff=40):
        atoms, not_cyclic, cyclic, all_cycles = self.atoms, [], [], []

        def prune():
            while 1:
                terminate = True
                for prune_atom in atoms:
                    if prune_atom.atom_number not in not_cyclic:
                        # for atom with cycle, it at least has 2 neighbors to form a cycle.
                        if len([neighbor for neighbor in prune_atom.neighbors if neighbor not in not_cyclic]) < 2:
                            not_cyclic.append(prune_atom.atom_number)
                            terminate = False
                if terminate:
                    break

        def search(start_index):
            # the paths is a list of list, which means each path is a list.
            paths, cycles = collections.deque([[start_index]]), []
            while paths:
                indices = paths.popleft()
                for neighbor in atoms[indices[-1]].neighbors:
                    if neighbor not in not_cyclic:
                        if neighbor not in indices and len(indices) < cutoff:
                            paths.append(list(indices) + [atoms[neighbor].atom_number])
                        elif neighbor == start_index and len(indices) != 2:
                            cycles.append(list(indices) + [atoms[neighbor].atom_number])
                            if short_circuit:
                                return cycles
            if cycles:
                return cycles
        prune()
        for atom in atoms:
            start_index = atom.atom_number
            if start_index not in not_cyclic:
                if not short_circuit or start_index not in cyclic:
                    search_result = search(start_index)
                    if search_result is None:
                        not_cyclic.append(start_index)
                        prune()
                    else:
                        all_cycles.append(search_result)
                        if short_circuit:
                            cyclic += [item for sublist in search_result for item in sublist]

        for cycle_list in all_cycles:
            for cycle in cycle_list:
                for index in cycle:
                    atoms[index].in_cycle = True
        
        if all_cycles:
            self.has_cycle = True

    def defineEZ(self):
        for bond in self.bonds:
            if bond.bond_type == "2":
                bond.bond_stereo = self.calculateEZ(bond)
    
    def calculateEZ(self, bond):

        vertical = False

        firstAtom = self.atoms[bond.first_atom_number]
        secondAtom = self.atoms[bond.second_atom_number]
        xf, yf = firstAtom.x, firstAtom.y
        xs, ys = secondAtom.x, secondAtom.y

        if xf != xs:
            slope = (yf-ys)/(xf-xs)
            b = yf - slope * xf 
        else:
            vertical = True

        firstIndex = [neighbor for neighbor in firstAtom.neighbors if neighbor != secondAtom.atom_number]
        secondIndex = [neighbor for neighbor in secondAtom.neighbors if neighbor != firstAtom.atom_number]
        
        if not len(firstIndex) or not len(secondIndex):
            return 0
        
        if len(firstIndex) == 1:
            firstH, firstL = self.atoms[firstIndex[0]], None
        else:
            firstH, firstL = self.getAtomRank(self.atoms[firstIndex[0]], self.atoms[firstIndex[1]], firstAtom)
        
        if len(secondIndex) == 1:
            secondH, secondL = self.atoms[secondIndex[0]], None
        else:
            secondH, secondL = self.getAtomRank(self.atoms[secondIndex[0]], self.atoms[secondIndex[1]], secondAtom)
        

        if (firstL and firstH.atom_number == firstL.atom_number) or (secondL and secondH.atom_number == secondL.atom_number):
            return 0

        if vertical:
            if (firstH.x - xs) * (secondH.x - xs) > 0:
                return 1
            else:
                return -1
        else:
            if (self.calculateY(slope, b, firstH) - firstH.y) * (self.calculateY(slope, b, secondH) - secondH.y) > 0:
                return 1
            else:
                return -1

    def calculateY(self, slope, b, atom):

        return atom.x * slope + b

    def getAtomRank(self, atomLeft, atomRight, pre):

        if atomicNumbers[atomLeft.default_symbol] > atomicNumbers[atomRight.default_symbol]:
            return [atomLeft, atomRight]
        elif atomicNumbers[atomLeft.default_symbol] < atomicNumbers[atomRight.default_symbol]:
            return [atomRight, atomLeft]
        else:
            leftNeighbors = [(atomLeft, pre.atom_number)]
            leftVisited = set([pre.atom_number])
            rightNeighbors = [(atomRight, pre.atom_number)]
            rightVisited = set([pre.atom_number])
            leftNeighborsAtomNumberList = self.getNeighborsAtomicNumber(leftNeighbors)
            rightNeighborsAtomNumberList = self.getNeighborsAtomicNumber(rightNeighbors)

            while leftNeighborsAtomNumberList == rightNeighborsAtomNumberList and leftNeighborsAtomNumberList:
                leftNeighborsNew = []
                for atom, preIndex in leftNeighbors:
                    leftVisited.add(atom.atom_number)
                    for neighbor in atom.neighbors:
                        if neighbor != preIndex and neighbor not in leftVisited:
                            leftNeighborsNew.append((self.atoms[neighbor], atom.atom_number))
                leftNeighbors = leftNeighborsNew
                leftNeighborsAtomNumberList = self.getNeighborsAtomicNumber(leftNeighbors)

                rightNeighborsNew = []
                for atom, preIndex in rightNeighbors:
                    rightVisited.add(atom.atom_number)
                    for neighbor in atom.neighbors:
                        if neighbor!= preIndex and neighbor not in rightVisited:
                            rightNeighborsNew.append((self.atoms[neighbor], atom.atom_number))
                rightNeighbors = rightNeighborsNew
                rightNeighborsAtomNumberList = self.getNeighborsAtomicNumber(rightNeighbors)


            if tuple(leftNeighborsAtomNumberList) > tuple(rightNeighborsAtomNumberList):
                return [atomLeft, atomRight]
            elif tuple(leftNeighborsAtomNumberList) < tuple(rightNeighborsAtomNumberList):
                return [atomRight, atomLeft]
            else:
                return [atomLeft, atomLeft]

    def getNeighborsAtomicNumber(self, leftNeighbors):


        leftList = []
        for atom, pre in leftNeighbors:
            leftList.extend([atomicNumbers[self.atoms[neighbor].default_symbol] for neighbor in atom.neighbors if neighbor != pre])
        leftList.sort(reverse=True)
        return leftList

    def isolatedCpd(self):

        for atom in self.atoms:
            if atom.bond_counts == 0:
                return True
        return False

    def connectedComponent(self):

        visited = set()
        group_id = 0

        def dfs(i, group_id):
            self.atoms[i].group_id = group_id
            visited.add(i)
            for neighbor in self.atoms[i].neighbors:
                dfs(neighbor, group_id)

        for index, atom in enumerate(self.atoms):
            if index not in visited:
                dfs(index, group_id)
                group_id += 1

        components = [[] for i in range(group_id)]
        for index, atom in enumerate(self.atoms):
            components[atom.group_id].append(index)
        return components

    def abnormalBond(self):

        abnormalAtoms = collections.defaultdict(list)
        for index, atom in enumerate(self.atoms):
            if atom.atom_symbol in atomBonds:
                bond_counts = atom.bond_counts
                bond_counts -= atom.charge
                if type(atomBonds[atom.atom_symbol]) is list:
                    if bond_counts not in atomBonds[atom.atom_symbol]:
                        abnormalAtoms[atom.atom_symbol].append(index)
                else:
                    if bond_counts > atomBonds[atom.atom_symbol]:
                        abnormalAtoms[atom.atom_symbol].append(index)
        return abnormalAtoms

    def curateInvalidN(self):

        abnormalAtoms = self.abnormalBond()
        if "N" in abnormalAtoms:
            for atom_index in abnormalAtoms["N"]:
                atom = self.atoms[atom_index]
                atom.charge += 1
                for neighbor in atom.neighbors:
                    if self.atoms[neighbor].atom_symbol == "O" and self.bond_lookup[(atom.atom_number, neighbor)].bond_type == "2":
                        self.atoms[neighbor].charge -= 1
                        self.bond_lookup[(atom.atom_number, neighbor)].bond_type = "1"

    def clearColor(self):

        for atom in self.atoms:
            atom.color_0 = ""
            atom.color = ""
            atom.color_layers = {}
        self.color_groups = collections.defaultdict(list)
    


    def firstRoundColor(self, isotope_resolved=False, stereo_resolved=False, charge=False, R=False, atom_stereo=False, strict=True, depth=5000):
        
        atomNeighbors = collections.defaultdict(list)
        color_groups = collections.defaultdict(list)
        excluded = self.isolatedMetalIndex + self.HIndex
        self.clearColor()

        for index, atom in enumerate(self.atoms):
            if index not in excluded:
                if not atom.is_R():
                    atom.color_0 += atom.atom_symbol
                else:
                    if R:
                        atom.color_0 += "R"
                    else:
                        continue
                if atom_stereo:
                    atom.color_0 += atom.atom_stereo_parity
                if charge:
                    atom.color_0 += str(atom.charge)
                if isotope_resolved:
                    atom.color_0 += atom.mass_difference
                atom.color = atom.color_0
                atomNeighbors[index].append(atom)
                color_groups[atom.color].append(atom)

        
        atomNeighborRepresentation = collections.defaultdict(str)
        numOfValidAtom = 0

        notComplete = []
        for index, atom in enumerate(self.atoms):
            if atom.color_0:
                color_components = []
                numOfValidAtom += 1
                for neighbor_index in atom.neighbors:
                    if self.atoms[neighbor_index].color_0:
                        connecting_bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                        bond_type = connecting_bond.bond_type
                        if not strict and bond_type == "2":
                            bond_type = "1"
                        if stereo_resolved:
                            color_components.append((self.atoms[neighbor_index].color_0, bond_type + connecting_bond.bond_stereo))
                        else:
                            color_components.append((self.atoms[neighbor_index].color_0, bond_type))
                notComplete.append(atom)
                atomNeighborRepresentation[index] = atom.color_0 + "".join(["(" + a + "," + b + ")" for (a, b) in sorted(color_components)])
        
        visited = collections.defaultdict(set)
        i = 0
        if depth != 5000:
            numOfValidAtom = depth

        while i < numOfValidAtom and notComplete:
        
            inner_color_groups = collections.defaultdict(list)
            for atom in notComplete:
                newNeighbors = []
                color_list = []
                for neighbor in atomNeighbors[atom.atom_number]:
                    color_list.append(atomNeighborRepresentation[neighbor.atom_number])
                    if neighbor.atom_number not in visited[atom.atom_number]:
                        visited[atom.atom_number].add(neighbor.atom_number)
                        newNeighbors.extend([self.atoms[newNeighbor_index] for newNeighbor_index in neighbor.neighbors if self.atoms[newNeighbor_index].color_0])
                if color_list:
                    color_list.sort()
                    atom.color += "(" + "".join(color_list) + ")"
                    atom.color_layers[i+1] = atom.color
                atomNeighbors[atom.atom_number] = newNeighbors
                inner_color_groups[atom.color].append(atom)
            if i > 3:
                newNotComplete = []
                for name in inner_color_groups.keys():
                    if len(inner_color_groups[name]) > 1:
                        newNotComplete.extend(inner_color_groups[name])
                notComplete = newNotComplete
            i += 1
        
        for atom in self.atoms:
            if atom.color_0:
                self.color_groups[atom.color].append(atom.atom_number)
    
    
    def checkSymmetry(self, stereo_resolved=False, strict=True):

        atomNeighborRepresentation = collections.defaultdict(str)
        for atom in self.atoms:
            if atom.color:
                color_components = []
                for neighbor_index in atom.neighbors:
                    if self.atoms[neighbor_index].color:
                        connecting_bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                        bond_type = connecting_bond.bond_type
                        if not strict and bond_type == "2":
                            bond_type = "1"
                        if stereo_resolved:
                            color_components.append((self.atoms[neighbor_index].color, bond_type + connecting_bond.bond_stereo))
                        else:
                            color_components.append((self.atoms[neighbor_index].color, bond_type))
                atomNeighborRepresentation[atom.atom_number] = atom.color + "".join(["(" + a + "," + b + ")" for (a, b) in sorted(color_components)])
        
        notValid = {}
        for name in self.color_groups.keys():
            if len(self.color_groups[name]) > 1:
                sameColors = self.color_groups[name]
                
                visited = collections.defaultdict(set)
                currentLayer = {atom_index: [atom_index] for atom_index in sameColors}
                count = 0
                flag = True
                while flag:
                    count += 1
                    targetatom_index = sameColors[0]
                    targetColorList = [atomNeighborRepresentation[atom_index] for atom_index in currentLayer[targetatom_index]]
                    targetColorList.sort()
                    targetColoring = ''.join(targetColorList)

                    for compareatom_index in currentLayer.keys():
                        if compareatom_index != targetatom_index:
                            compareColorList = [atomNeighborRepresentation[atom_index] for atom_index in currentLayer[compareatom_index]]
                            compareColorList.sort()
                            compareColoring = ''.join(compareColorList)

                            if targetColoring != compareColoring:
                                notValid[name] = (sameColors, count)
                                flag = False
                                break
                    if flag:
                        for atom_index in currentLayer.keys():
                            nextNeighbors = []
                            for neighbor_index in currentLayer[atom_index]:
                                if neighbor_index not in visited[atom_index]:
                                    visited[atom_index].add(neighbor_index)
                                    nextNeighbors.extend([newIndex for newIndex in self.atoms[neighbor_index].neighbors if self.atoms[newIndex].color])
                            currentLayer[atom_index] = nextNeighbors
                        targetLength = len(currentLayer[targetatom_index])
                        for atom_index in currentLayer.keys():
                            if len(currentLayer[atom_index]) != targetLength:
                                notValid[name] = (sameColors, count)
                                flag = False
                                break
                        if targetLength == 0:
                            flag = False
        return notValid

    def curateSymmetry(self, stereo_resolved=False, strict=True):

        notValid = self.checkSymmetry()
        if not notValid:
            return
        
        atomNeighborRepresentation = collections.defaultdict(str)

        for atom in self.atoms:
            if atom.color:
                color_components = []
                for neighbor_index in atom.neighbors:
                    if self.atoms[neighbor_index].color:
                        connecting_bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                        bond_type = connecting_bond.bond_type
                        if not strict and bond_type == "2":
                            bond_type = "1"
                        if stereo_resolved:
                            color_components.append((self.atoms[neighbor_index].color_0, bond_type + connecting_bond.bond_stereo))
                        else:
                            color_components.append((self.atoms[neighbor_index].color_0, bond_type))
                atomNeighborRepresentation[atom.atom_number] = atom.color + "".join(["(" + a + "," + b + ")" for (a, b) in sorted(color_components)])

        
        for name in notValid:
            del self.color_groups[name]

            notSameIndex, layer = notValid[name]
            
            for targetatom_index in notSameIndex:
                visited = set()
                currentLayer = [targetatom_index]
                targetColor = ""
                count = layer + 1
                while count  and currentLayer:
                    curColorList = [atomNeighborRepresentation[atom_index] for atom_index in currentLayer]
                    curColorList.sort()
                    targetColor += "".join(curColorList)
                    nextNeighbors = []
                    for neighbor_index in currentLayer:
                        if neighbor_index not in visited:
                            visited.add(neighbor_index)
                            nextNeighbors.extend([newIndex for newIndex in self.atoms[neighbor_index].neighbors if self.atoms[newIndex].color])
                    currentLayer = nextNeighbors
                    count -= 1
                self.color_groups[targetColor].append(targetatom_index)
                self.atoms[targetatom_index].color = targetColor

        if self.checkSymmetry():
            print(self.entry)
        #print(self.color_groups)

    def colorMetals(self, charge=False, isotope_resolved=False, stereo_resolved=False, strict=True):
        
        for atom_index in self.isolatedMetalIndex:
            atom = self.atoms[atom_index]
            colorComponents = collections.defaultdict(int)
            for neighbor_index in atom.neighbors:
                if self.atoms[neighbor_index].color:
                    connecting_bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                    bond_type = connecting_bond.bond_type
                    if not strict and bond_type == "2":
                        bond_type = "1"
                    if stereo_resolved:
                        component = "({0}, {1})".format(self.atoms[neighbor_index].color, bond_type+connecting_bond.bond_stereo)
                    else:
                        component = "({0}, {1})".format(self.atoms[neighbor_index].color, bond_type)
                    colorComponents[component] += 1
            atom.color_0 = atom.atom_symbol
            if charge:
                atom.color_0 += str(atom.charge)
            if isotope_resolved:
                atom.color_0 += atom.mass_difference
            atom.color = atom.color_0
            for name in sorted(colorComponents.keys()):
                atom.color += "({0})({1})".format(colorComponents[name], name)

    def colorH(self, stereo_resolved=False, charge=False, isotope_resolved=False):

        for index in self.HIndex:
            atom = self.atoms[index]
            colorComponents = []
            for neighbor_index in atom.neighbors:
                connecting_bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                if stereo_resolved:
                    colorComponents.append((self.atoms[neighbor_index].color, connecting_bond.bond_type + connecting_bond.bond_stereo))
                else:
                    colorComponents.append((self.atoms[neighbor_index].color, connecting_bond.bond_type))
            colorComponents.sort()
            atom.color_0 = atom.atom_symbol
            if charge:
                atom.color_0 += str(atom.charge)
            if isotope_resolved:
                atom.color_0 += atom.mass_difference
            atom.color = atom.color_0 + "".join(["(" + a + "," + b + ")" for (a, b) in sorted(colorComponents)])


    def changeAromaticBonds(self, cycles):

        atomInCyle = [atom for cycle in cycles for atom in cycle]
        for cycle in cycles:
            aromaticBonds = self.extractAromaticBond(cycle)
            for bond in aromaticBonds:
                bond.bond_type = "4"
        outBonds = self.extractOutdouble_bond_counts(atomInCyle)
        for bond in outBonds:
            bond.bond_type = "1"

    def extractOutdouble_bond_counts(self, atomInCycle):
        double_bond_countsOutCycle = []
        for atom_index in atomInCycle:
            if self.atoms[atom_index].atom_symbol == "C":
                for neighbor_index in self.atoms[atom_index].neighbors:
                    if neighbor_index not in atomInCycle:
                        if self.bond_lookup[(atom_index, neighbor_index)].bond_type == "2":
                            double_bond_countsOutCycle.append(self.bond_lookup[(atom_index, neighbor_index)])
        return double_bond_countsOutCycle

    def extractAromaticBond(self, cycle):
        aromaticBonds = []
        allPairs = list(itertools.combinations(cycle, 2))
        visited = set()
        for pair in allPairs:
            if (pair[0], pair[1]) not in visited and (pair[0], pair[1]) in self.bond_lookup:
                aromaticBonds.append(self.bond_lookup[(pair[0], pair[1])])
                visited.add((pair[0], pair[1]))
                visited.add((pair[1], pair[1]))
        return aromaticBonds

    def changedouble_bond_countsStereo(self, bonds):
        for bond in bonds:
            firstAtom = bond.first_atom_number
            secondAtom = bond.second_atom_number
            if self.bond_lookup[(firstAtom, secondAtom)].bond_type == "2":
                self.bond_lookup[(firstAtom, secondAtom)].bond_stereo = str(bond.bond_stereo)

    def metalString(self, noBond=True):
        if noBond:  
            stringCounter = collections.Counter([self.atoms[index].color_0 for index in self.isolatedMetalIndex])
        else:
            stringCounter = collections.Counter([self.atoms[index].color for index in self.isolatedMetalIndex])
        stringRepresent = ""
        #print(stringCounter)
        keys = sorted(stringCounter)
        for key in keys:
            stringRepresent += "({0})({1})".format(stringCounter[key], key)
        return stringRepresent
    
    def HString(self, details=False):
        if details:
            stringCounter = collections.Counter([self.atoms[index].color for index in self.HIndex])
        else:
            stringCounter = collections.Counter([self.atoms[index].color_0 for index in self.HIndex])
        stringRepresent = ""
        keys = sorted(stringCounter)
        for key in keys:
            stringRepresent += "({0})({1})".format(stringCounter[key], key)
        return stringRepresent

    @property
    def colorMetalNoBond(self):
        return self.backBoneString + self.metalString()

    @property
    def colorMetal(self):
        return self.backBoneString + self.metalString(noBond=False)
    
    @property
    def backBoneString(self):
        stringRepresent = ""
        keys = sorted(self.color_groups)
        for key in keys:
            stringRepresent += "({0})({1})".format(len(self.color_groups[key]), key)
        return stringRepresent

    def getMetal(self):

        metals = collections.defaultdict(list)
        for index, atom in enumerate(self.atoms):
            if atom.default_symbol in metalSymbols:
                metals[atom.default_symbol].append(index)
        return metals

    @property
    def isolatedMetalIndex(self):
        
        return [index for index, atom in enumerate(self.atoms) if atom.atom_symbol in metalSymbols]
    
    @property
    def HIndex(self):

        return [index for index, atom in enumerate(self.atoms) if atom.atom_symbol == "H"]


    @property
    def formula(self):
        
        res = ""
        counter = collections.defaultdict(int)
        for atom in self.atoms:
            if atom.default_symbol == "H":
                continue
            elif atom.is_R():
                counter["R"] += 1
            else:
                counter[atom.atom_symbol] += 1
        for char in sorted(counter.keys()):
            res += char + str(counter[char])
        return res

    def containsR(self):
        for atom in self.atoms:
            if atom.is_R():
                return True
        return False


