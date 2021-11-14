
Optical material properties
===========================

.. currentmodule:: arc.materials



Generic optical material
------------------------

This is abstract class that all specific optical materials inherit and implement
their own implementation of methods.

.. autoclass:: arc.materials.OpticalMaterial
    :members: __init__
    :exclude-members: __init__ 

    .. rubric:: Methods

    .. autosummary::
        :toctree: generated/

        ~OpticalMaterial.getN

    .. rubric:: Attributes

    .. autosummary::
        :toctree: generated/

        ~OpticalMaterial.name
        ~OpticalMaterial.sources
        ~OpticalMaterial.sourcesComment
        ~OpticalMaterial.sourcesRange


Solids
------

.. autoclass:: arc.materials.Sapphire
    :members: __init__
    :exclude-members: __init__

    .. rubric:: Methods

    .. autosummary::
        :toctree: generated/

        ~Sapphire.getN

    .. rubric:: Attributes

    .. autosummary::
        

        ~Sapphire.name
        ~Sapphire.sources
        ~Sapphire.sourcesComment
        ~Sapphire.sourcesRange


Gases
-----

.. autoclass:: arc.materials.Air
    :members: __init__
    :exclude-members: __init__

    .. rubric:: Methods

    .. autosummary::
        :toctree: generated/

        ~Air.getN

    .. rubric:: Attributes

    .. autosummary::
        :toctree: generated/

        ~Air.name
        ~Air.sources
        ~Air.sourcesComment
        ~Air.sourcesRange