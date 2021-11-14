Alkali atom functions
=====================

.. currentmodule:: arc.alkali_atom_functions


.. autoclass:: AlkaliAtom
    :members: __init__
    :exclude-members: __init__

    .. autosummary::
        :toctree: generated/


Elementary properties
---------------------

    .. rubric:: Methods

    .. autosummary::
        :toctree: generated/

        ~AlkaliAtom.getPressure
        ~AlkaliAtom.getNumberDensity
        ~AlkaliAtom.getAverageInteratomicSpacing
        ~AlkaliAtom.getAverageSpeed

    .. rubric:: Attributes

    .. autosummary::
        :toctree: generated/

        ~AlkaliAtom.I
        ~AlkaliAtom.Z
        ~AlkaliAtom.abundance
        ~AlkaliAtom.elementName
        ~AlkaliAtom.groundStateN
        ~AlkaliAtom.extraLevels
        
        ~AlkaliAtom.mass
        ~AlkaliAtom.meltingPoint
    

Internal structure of atom states
---------------------------------

    .. rubric:: Methods

    .. autosummary::
        :toctree: generated/

        ~AlkaliAtom.corePotential
        ~AlkaliAtom.effectiveCharge
        ~AlkaliAtom.potential
        ~AlkaliAtom.radialWavefunction

    .. rubric:: Attributes

    .. autosummary::
        :toctree: generated/

        ~AlkaliAtom.a1
        ~AlkaliAtom.a2
        ~AlkaliAtom.a3
        ~AlkaliAtom.a4
        ~AlkaliAtom.rc
        ~AlkaliAtom.alphaC
        ~AlkaliAtom.cpp_numerov


Energies of atom states
-----------------------

    .. rubric:: Methods

    .. autosummary::
        :toctree: generated/

        ~AlkaliAtom.getEnergy
        ~AlkaliAtom.getZeemanEnergyShift
        ~AlkaliAtom.getQuantumDefect
        ~AlkaliAtom.breitRabi

    .. rubric:: Attributes

    .. autosummary::
        :toctree: generated/
 
        ~AlkaliAtom.gI
        ~AlkaliAtom.gL       
        ~AlkaliAtom.hyperfineStructureData 

        ~AlkaliAtom.levelDataFromNIST
        ~AlkaliAtom.sEnergy
        ~AlkaliAtom.quantumDefect
        ~AlkaliAtom.minQuantumDefectN


Transitions between states
--------------------------


    .. rubric:: Methods

    .. autosummary::
        :toctree: generated/

        ~AlkaliAtom.getDipoleMatrixElement
        ~AlkaliAtom.getDipoleMatrixElementHFS
        ~AlkaliAtom.getTransitionWavelength
        ~AlkaliAtom.getTransitionFrequency
        ~AlkaliAtom.getRabiFrequency
        ~AlkaliAtom.getRabiFrequency2
        ~AlkaliAtom.getStateLifetime
        ~AlkaliAtom.getTransitionRate
        ~AlkaliAtom.getReducedMatrixElementJ_asymmetric
        ~AlkaliAtom.getReducedMatrixElementJ
        ~AlkaliAtom.getReducedMatrixElementL
        ~AlkaliAtom.getRadialMatrixElement
        ~AlkaliAtom.getQuadrupoleMatrixElement

        ~AlkaliAtom.getC6term
        ~AlkaliAtom.getC3term
        ~AlkaliAtom.getEnergyDefect
        ~AlkaliAtom.getEnergyDefect2
        ~AlkaliAtom.updateDipoleMatrixElementsFile
        ~AlkaliAtom.getRadialCoupling
        ~AlkaliAtom.getLiteratureDME

        ~AlkaliAtom.getSphericalMatrixElementHFStoFS
        ~AlkaliAtom.getDipoleMatrixElementHFStoFS
        ~AlkaliAtom.getMagneticDipoleMatrixElementHFS
        ~AlkaliAtom.getHFSCoefficients
        ~AlkaliAtom.getHFSEnergyShift
        ~AlkaliAtom.getBranchingRatio
        ~AlkaliAtom.getSaturationIntensity
        ~AlkaliAtom.getSaturationIntensityIsotropic
        ~AlkaliAtom.groundStateRamanTransition
        ~AlkaliAtom.twoPhotonRydbergExcitation
        ~AlkaliAtom.getLandegj
        ~AlkaliAtom.getLandegjExact
        ~AlkaliAtom.getLandegf
        ~AlkaliAtom.getLandegfExact
        ~AlkaliAtom.breitRabi

    .. rubric:: Attributes

    .. autosummary::
        :toctree: generated/

        ~AlkaliAtom.literatureDMEfilename
        ~AlkaliAtom.dipoleMatrixElementFile
        ~AlkaliAtom.quadrupoleMatrixElementFile
        

        
        

