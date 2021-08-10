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
from pwem.protocols import EMProtocol
from pyworkflow.protocol import PointerParam, EnumParam


class ProtZrankProtein(EMProtocol):
    """
    Protocol to perform protein-protein docking using the ZRANK docking tool
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
        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
        pass