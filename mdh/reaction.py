#!/usr/bin/python3

"""
MDH.reaction
~~~~~~~~~~~~

This module provides the :class:`~MDH.reaction.Reaction` class entity.

"""

class Reaction:

    """ Reaction class describes the :class:`~MDH.reaction.Reaction` entity. """

    def __init__(self, reaction_name, one_side_compounds, the_other_side_compounds, ecs, atom_mappings, coefficients):
        """
        Reaction initializer.

        :param reaction_name: the reaction name.
        :type reaction_name: :py:class:`str`.
        :param one_side_compounds: the list of :class:`~MDH.compound.Compound` entities in one side of the reaction.
        :type one_side_compounds: :py:class:`list`.
        :param the_other_side_compounds: the list of :class:`~MDH.compound.Compound` entities in the other side of the reaction.
        :type the_other_side_compounds: :py:class:`list`.
        :param ecs: the list of Enzyme Commission numbers (EC numbers) of the reaction.
        :type ecs: :py:class:`list`.
        :param atom_mappings: the list of atom mappings between two sides of the reaction.
        :type atom_mappings: :py:class:`list`.
        :param coefficients: the dictionary of compound names and their corresponding coefficients in the reaction.
        :type coefficients: :py:class:`dict`.
        """
        self.reaction_name = reaction_name
        self.one_compound = one_side_compounds
        self.the_other_compound = the_other_side_compounds
        self.ecs = ecs
        self.atom_mappings = atom_mappings
        self.coefficients = coefficients
        
    @property
    def name(self):
        """
        To get the reaction name.

        :return: the reaction name.
        :rtype: :py:class:`str`.
        """
        return self.reaction_name







        