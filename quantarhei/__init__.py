# -*- coding: utf-8 -*-
"""
    Quantarhei User Level Classes and Objects
    =========================================

    In Quantarhei, classes are losely grouped into three categories. First,
    there is agroup of classes, which represent basic concepts of quantum
    mechanics, provide access to important implementations of spectroscopic
    simulations and dynamics of open quantum systems, and classes which allow
    basic management of the simulation environment and numerical results.
    These classed are called **user level classes**, and they are all
    accessible in highest namespace level of the Quantarhei package.
    If you import Quantarhei like this:

    >>> import quantarhei as qr

    you can access user level classes through the qr. prefix, e.g.


    >>> manager = qr.Manager()
    >>> print(manager.version)
    0.0.35

    The list of user level classes is provided below. Tue latest and most
    uptodate information can be obtained by viewing the source code of the 
    root `__init__.py` file of the packages. All classes imported there are
    considered user level classes.
    
    
    Other Class Levels
    ------------------
    
    In this documentation we recognize two more groups (or levels) of classes.
    More specialized classes, which normal user does not need as often as the
    user level classes are called **advanced level classes**. These use the
    second level name space. For instance the class `SystemBathInteraction`
    is relatively rarely used directly. It is therefore *hidden* in the name
    space `qm` (as quantum mechanics) of the package. This class can be 
    instantiated e.g. like this
    
    >>> import quantarhei as qr
    >>> sbi = qr.qm.SystemBathInteraction()
    
    Advanced level classes are still intendend for relatively frequent use
    by the user. However, in order to reduced the *apparent* complexity of
    basic usage of Quantarhei, advanced level classes are documented in their
    respective sub-packages, one level deeper than user level classes. Complete
    documentation of advanced level classes is available in the Advanced Level
    Classes section of this documentation.
    
    Everything else in Quantarhei package goes under the banner of
    **expert level classes**. This includes all classes and objects used 
    internally in Quantarhei. We make every effort to document also this part
    of the package as completely as possible, but it is the last item on the
    list, so to say. The user is welcome to learn and use the expert level
    classes, but our aim is to structure Quantarhei in such a way, that this
    is not necessary. More on expert level classes in the section in 
    Quantarhei internals.

    User Level Objects and Convenience Functions
    ============================================
    
    Besides classes, Quantarhei also defines some user level objects and
    convenience functions. They are listed here under several categories
  
    Numeric types
    -------------

    .. toctree::
        :maxdepth: 2
        
        functions/numtypes        
    
    Convenience Functions
    ---------------------
    
    .. toctree::
        :maxdepth: 2
        
        functions/convenience
        
    
    Logging Functions and Loglevels
    -------------------------------

    .. toctree::
        :maxdepth: 2
        
        functions/logging

   
    .. 
        Builders
        --------
        
        Mode .......... represents a harmonic vibrational mode of a molecule
        Molecule ...... represents a molecule
        Aggregate ..... represents an aggregate of molecules
        PDBFile ....... reader and writter of structures from PDB format
        Disorder ...... class managing static disorder of molecular transition
                        energies
         
        Core classes
        ------------
        
        TimeAxis ......... linear axis of real values representing discrete time
        FrequencyAxis .... linear axis of real values representing discrete
                           frequency axis
        DFunction ........ discrete function
        
        
        Various managers
        ----------------
        
        Manager ............ the main behind-the-scenes manager of the package
        energy_units ....... energy units manager for use with the "with" construct
        frequency_units .... frequency units manager for use with 
                             the "with" construct
        eigenbasis_of ...... manager of the basis transformations to be used with 
                             the "with" construct
        set_current_units .. function to set current units globally
        
        ... to be continued


"""


###############################################################################
#
#
#            Imports of high level classes and functions 
#
#
###############################################################################

#
# Fix used numerical types
#
#import numpy
from .core.managers import Manager
m = Manager()

REAL = m.get_real_type() #numpy.float64
COMPLEX = m.get_complex_type() #numpy.complex128

LOG_URGENT = 0
LOG_REPORT = 3
LOG_INFO = 5
LOG_DETAIL = 7
LOG_QUICK = 9

#
# Builders
#
from .builders.modes import Mode
from .builders.molecules import Molecule
from .builders.molecule_test import TestMolecule
from .builders.aggregates import Aggregate
from .builders.aggregate_test import TestAggregate
from .builders.pdb import PDBFile
from .builders.disorder import Disorder

#
# Core classes
#
from .core.time import TimeAxis
from .core.frequency import FrequencyAxis
from .core.valueaxis import ValueAxis
from .core.dfunction import DFunction
#from .core.saveable import Saveable

from .core.saveable import Saveable
from .core.parcel import Parcel

#
# Various managers
#
from .core.managers import energy_units
from .core.managers import frequency_units
from .core.managers import length_units
from .core.managers import eigenbasis_of
from .core.managers import set_current_units

#
# Parallelization
#
from .core.parallel import distributed_configuration
from .core.parallel import start_parallel_region
from .core.parallel import close_parallel_region
from .core.parallel import parallel_function
from .core.parallel import block_distributed_range

###############################################################################
#                            SPECTROSCOPY
###############################################################################

#
# Linear absorption 
#
from .spectroscopy.abs2 import AbsSpectrum
from .spectroscopy.abscontainer import AbsSpectrumContainer
from .spectroscopy.abscalculator import AbsSpectrumCalculator
#
# Fluorescence
#
from .spectroscopy.fluorescence import FluorSpectrum
from .spectroscopy.fluorescence import FluorSpectrumContainer
from .spectroscopy.fluorescence import FluorSpectrumCalculator
#
# Linear dichroism
#
from .spectroscopy.linear_dichroism import LinDichSpectrum
from .spectroscopy.linear_dichroism import LinDichSpectrumContainer
from .spectroscopy.linear_dichroism import LinDichSpectrumCalculator
#
# Circular dichroism
#
from .spectroscopy.circular_dichroism import CircDichSpectrum
from .spectroscopy.circular_dichroism import CircDichSpectrumContainer
from .spectroscopy.circular_dichroism import CircDichSpectrumCalculator
#
# Fourier transform Two-Dimensional Spectra
#
from .spectroscopy.twod2 import TwoDSpectrum
from .spectroscopy.twodcontainer import TwoDSpectrumContainer
from .spectroscopy.twod2 import TwoDSpectrumCalculator
from .spectroscopy.twod2 import MockTwoDSpectrumCalculator


from .spectroscopy.pathwayanalyzer import LiouvillePathwayAnalyzer

from .spectroscopy.labsetup import LabSetup

###############################################################################
#                           QUANTUM MECHANICS
###############################################################################


#
# Operators
#
from .qm import StateVector
from .qm import DensityMatrix
from .qm import ReducedDensityMatrix
from .qm import BasisReferenceOperator
from .qm import Hamiltonian
from .qm import TransitionDipoleMoment

#
# Propagators
#
from .qm.propagators.poppropagator import PopulationPropagator
from .qm.propagators.svpropagator import StateVectorPropagator
from .qm import ReducedDensityMatrixPropagator

#
# Evolutions (time-dependent operators)
#
from .qm.propagators.statevectorevolution import StateVectorEvolution
from .qm import DensityMatrixEvolution
from .qm import ReducedDensityMatrixEvolution

#
# System-bath interaction
#
from .qm.corfunctions import CorrelationFunction
from .qm.corfunctions import SpectralDensity



###############################################################################
# Convenience functions
###############################################################################
#from .core.saveable import load
#from .core.saveable import read_info

from .core.parcel import save_parcel
from .core.parcel import load_parcel
from .core.parcel import check_parcel

from .core.units import convert
from .core.units import in_current_units

from .utils.vectors import normalize2
from .utils.vectors import norm 

from .utils.logging import printlog
from .utils.logging import loglevels2bool
from .utils.logging import log_urgent
from .utils.logging import log_report
from .utils.logging import log_info
from .utils.logging import log_detail
from .utils.logging import log_quick



