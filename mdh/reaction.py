#!/usr/bin/python3

"""
MDH.reaction
~~~~~~~~~~~~

This module provides the :class:`~mdh.reaction.Reaction` class entity.

"""


class Reaction:

    """ Reaction class describes the :class:`~mdh.reaction.Reaction` entity. """

    def __init__(self, reaction_name: str, one_side: list, the_other_side: list, ecs: dict, atom_mappings: list,
                 coefficients: dict) -> None:
        """
        Reaction initializer.

        :param reaction_name: the reaction name.
        :param one_side: the list of :class:`~mdh.compound.Compound` entities in one side of the reaction.
        :param the_other_side: the list of :class:`~mdh.compound.Compound` entities in the other side of the reaction.
        :param ecs: the dict of Enzyme Commission numbers (EC numbers) of the reaction.
        :param atom_mappings: the list of atom mappings between two sides of the reaction.
        :param coefficients: the dictionary of compound names and their corresponding coefficients in the reaction.
        """
        self.name = reaction_name
        self.one_side = one_side
        self.the_other_side = the_other_side
        self.ecs = ecs
        self.atom_mappings = atom_mappings
        self.coefficients = coefficients
        
    @property
    def name(self) -> str:
        """
        To get the reaction name.

        :return: the reaction name.
        """
        return self.name
