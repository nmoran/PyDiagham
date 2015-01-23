#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension
import os
import os.path
import platform
import sys

os.environ['CC'] = 'g++'; 
os.environ['CXX'] = 'g++';
os.environ['CPP'] = 'g++'; 
os.environ['LDSHARED'] = 'g++';

#diagham_src='/home/niall/local/Diagham/trunk/'
#diagham_build='/home/niall/local/Diagham/builds/debug_main/'

diagham_build='/home/niall/local/Diagham/builds/debug_main/'
diagham_src='/home/niall/local/Diagham/trunk/'

if platform.system() == 'Darwin' :
  shared_flag = '-bundle'
else:
  shared_flag = '-shared'
  
  
pythonlib = 'python%d.%d'%(sys.version_info[0],sys.version_info[1])
pythonlib = 'python2.6'

diagham_module = Extension('_diagham',
#                           sources=['diagham_wrap.cxx', '../../trunk/FQHE/src/HilbertSpace/BosonOnSphereTwoLandauLevels.cc'],
                           sources=['diagham_wrap.cxx','utils.cxx'],
			   			   include_dirs=[os.path.join(diagham_src,'FQHE/src/HilbertSpace/'),os.path.join(diagham_src,'src/Vector/'),os.path.join(diagham_src,'src/Matrix/'),os.path.join(diagham_src,'src/'),os.path.join(diagham_src,'FQHE/src/'),os.path.join(diagham_build,'src/')],
#			   			   library_dirs=['/Library/Frameworks/EPD64.framework/Versions/Current/lib'],
#			   			   libraries=['QHEHilbertSpace','QHEHamiltonian','Array','Matrix'],
						   libraries=[pythonlib,'lapack','blas'],
						   extra_compile_args=['-fPIC'],
#			   			   extra_link_args=['../debug_main/FQHE/src/HilbertSpace/libQHEHilbertSpace.a','../debug_main/FQHE/src/Hamiltonian/libQHEHamiltonian.a','../debug_main/src/Matrix/libArray.a','../debug_main/src/Matrix/libMatrix.a','../debug_main/src/Vector/libVector.a','../debug_main/src/GeneralTools/libGeneralTools.a',]
			   			   extra_link_args=		       [os.path.join(diagham_build,'Base/src/BitmapTools/BitmapPicture/libBitmapPicture.a'),
											os.path.join(diagham_build,'Base/src/BitmapTools/Color/libColor.a'),
											os.path.join(diagham_build,'FQHE/src/Architecture/ArchitectureOperation/libQHEArchitectureOperation.a'),
											os.path.join(diagham_build,'FQHE/src/FunctionBasis/libQHEFunctionBasis.a'),
											os.path.join(diagham_build,'FQHE/src/Hamiltonian/libQHEHamiltonian.a'),
											os.path.join(diagham_build,'FQHE/src/HilbertSpace/libQHEHilbertSpace.a'),
											os.path.join(diagham_build,'FQHE/src/MainTask/libQHEMainTask.a'),
											os.path.join(diagham_build,'FQHE/src/Operator/libQHEOperator.a'),
											os.path.join(diagham_build,'FQHE/src/QuantumNumber/libQHEQuantumNumber.a'),
											os.path.join(diagham_build,'FQHE/src/Tools/FQHEFiles/libQHEFiles.a'),
											os.path.join(diagham_build,'FQHE/src/Tools/FQHEMonteCarlo/libFQHEMonteCarlo.a'),
											os.path.join(diagham_build,'FQHE/src/Tools/FQHESpectrum/libQHESpectrum.a'),
											os.path.join(diagham_build,'FQHE/src/Tools/FQHEWaveFunction/libQHEWaveFunction.a'),
											os.path.join(diagham_build,'src/Architecture/ArchitectureOperation/libArchitectureOperation.a'),
											os.path.join(diagham_build,'src/Architecture/ClusterArchitecture/libClusterArchitecture.a'),
											os.path.join(diagham_build,'src/Architecture/libArchitecture.a'),
											os.path.join(diagham_build,'src/DMRGAlgorithm/libDMRGAlgorithm.a'),
											os.path.join(diagham_build,'src/FunctionBasis/libFunctionBasis.a'),
											os.path.join(diagham_build,'src/GeneralTools/libGeneralTools.a'),
											os.path.join(diagham_build,'src/Hamiltonian/DMRGHamiltonian/libDMRGHamiltonian.a'),
											os.path.join(diagham_build,'src/Hamiltonian/libHamiltonian.a'),
											os.path.join(diagham_build,'src/HilbertSpace/DMRGHilbertSpace/libDMRGHilbertSpace.a'),
											os.path.join(diagham_build,'src/HilbertSpace/libHilbertSpace.a'),
											os.path.join(diagham_build,'src/HilbertSpace/ManyBodyHilbertSpace/libManyBodyHilbertSpace.a'),
											os.path.join(diagham_build,'src/Interaction/InternalInteraction/libInternalInteraction.a'),
											os.path.join(diagham_build,'src/Interaction/libInteraction.a'),
											os.path.join(diagham_build,'src/LanczosAlgorithm/libLanczosAlgorithm.a'),
											os.path.join(diagham_build,'src/MainTask/libMainTask.a'),
											os.path.join(diagham_build,'src/MathTools/libMathTools.a'),
											os.path.join(diagham_build,'src/MathTools/NumericalAnalysis/libNumericalAnalysis.a'),
											os.path.join(diagham_build,'src/MathTools/RandomNumber/libRandomNumber.a'),
											os.path.join(diagham_build,'src/Matrix/libArray.a'),
											os.path.join(diagham_build,'src/Matrix/libMatrix.a'),
											os.path.join(diagham_build,'src/MCObservables/libMCObservables.a'),
											os.path.join(diagham_build,'src/Operator/libOperator.a'),
											os.path.join(diagham_build,'src/Options/libOptions.a'),
											os.path.join(diagham_build,'src/Output/libOutput.a'),
											os.path.join(diagham_build,'src/Polynomial/libPolynomial.a'),
											os.path.join(diagham_build,'src/QuantumNumber/libQuantumNumber.a'),
											os.path.join(diagham_build,'src/Tensor/libTensor.a'),
											os.path.join(diagham_build,'src/TensorProduct/libTensorProduct.a'),
											os.path.join(diagham_build,'src/Vector/libVector.a'),
											shared_flag,
											'-fPIC']
                           )

setup (name = 'diagham',
       version = '0.1',
       author      = "Niall Moran",
       description = """Simple swig example from docs""",
       ext_modules = [diagham_module],
       py_modules = ["diagham"],
       )

