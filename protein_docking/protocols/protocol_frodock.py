# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os

from protein_docking.constants import FRODOCKGRID
from pwem.protocols import EMProtocol
from pyworkflow.protocol import PointerParam, EnumParam

from protein_docking import Plugin
from pyworkflow.utils import yellowStr


class ProtFrodockProtein(EMProtocol):
    """
    Protocol to perform protein-protein docking using the FRODOCK docking tool
    """

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputPdbReceptor', PointerParam,
                      pointerClass='AtomStruct',
                      label="Receptor pdb", important=True,
                      help='The receptor pdb')
        form.addParam('inputPdbLigand', PointerParam,
                      pointerClass='AtomStruct',
                      label="ligand pdb", important=True,
                      help='The ligand pdb')
        form.addParam('interactionType', EnumParam,
                      choices=['Enzyme-Substrate', 'Antigen-Antibody', 'Unknown'],
                      default=2,
                      label="Type of interaction",
                      help='Type of interaction')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.mapGenerationStep)
        self._insertFunctionStep(self.dockingSearchStep)
        self._insertFunctionStep(self.clusteringStep)
        self._insertFunctionStep(self.createOutputStep)

    def mapGenerationStep(self):
        """
        All necessary potential maps must be pre-computed using FRODOCKGRID.
        Although vdw and electrostatics maps could be computed on the fly
        during the docking search, it is recommendable to create the maps
        beforehand in order to visualize them and check that they are consistent
        with the original structure. The precomputation of desolvation potential
        maps for receptor and ligand is always required
        """
        receptorPdbPath = os.path.abspath(self.inputPdbReceptor.get().getFileName())
        ligandPdbPath = os.path.abspath(self.inputPdbLigand.get().getFileName())
        interactionDict = ['E', 'A', 'None']
        interactionType = interactionDict[self.interactionType.get()]
        program = self._getProgram(FRODOCKGRID)

        # Creation of receptor vdw potential map
        print(yellowStr('Creation of receptor vdw potential map'), flush=True)
        program, args = self.getFrodockGridCommand(program=program,
                                                   pdbFile=receptorPdbPath,
                                                   outputSuffix='_W.ccp4')
        Plugin.runProgram(program, args)

        # Creation of the receptor electrostatic potential map
        print(yellowStr('Creation of the receptor electrostatic potential map'),
              flush=True)
        program, args = self.getFrodockGridCommand(program=program,
                                                   pdbFile=receptorPdbPath,
                                                   outputSuffix='_E.ccp4',
                                                   mValue=1,
                                                   tValue=interactionType)

        Plugin.runProgram(program, args)

        # Creation of the receptor desolvation potential map
        print(yellowStr('Creation of the receptor desolvation potential map'), flush=True)
        program, args = self.getFrodockGridCommand(program=program,
                                                   pdbFile=receptorPdbPath,
                                                   outputSuffix='_DS.ccp4',
                                                   mValue=3)

        Plugin.runProgram(program, args)

        # Creation of the ligand desolvation potential map
        print(yellowStr('Creation of the ligand desolvation potential map'),
              flush=True)
        program, args = self.getFrodockGridCommand(program=program,
                                                   pdbFile=ligandPdbPath,
                                                   outputSuffix='_DS.ccp4',
                                                   mValue=3)

        Plugin.runProgram(program, args)


    def dockingSearchStep(self):
        """Creation of the receptor electrostatic potential map"""
        pass

    def clusteringStep(self):
        pass

    def createOutputStep(self):
        pass

    # -----------------------Utils funcitons-------------------------------

    def _getProgram(self, programName):
        """ Return program binary. """
        return Plugin.getProgram(programName)

    def getFrodockGridCommand(self,  **kwargs):
        program = kwargs.get('program')
        pdbInputFile = kwargs.get('pdbFile')
        mValue = kwargs.get('mValue', None)
        tValue = kwargs.get('tValue')
        mValue = ' -m %s' % mValue if mValue is not None else ''
        tValue = ' -t %s' % tValue if tValue is not None else ''
        outputSuffix = kwargs.get('outputSuffix')

        outputFileName = os.path.basename(pdbInputFile).split('.')[0] + outputSuffix
        outputPdbFilePath = os.path.abspath(self._getExtraPath(outputFileName))

        params = '%s -o %s%s%s' % (pdbInputFile, outputPdbFilePath, mValue, tValue)
        return program,  params



