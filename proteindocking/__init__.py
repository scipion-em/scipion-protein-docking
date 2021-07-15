# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *          Erney Ramirez Aportela     (eramirez@cnb.csic.es)
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
from pyworkflow.utils import Environ, runJob, greenStr
import pwem
from .constants import PROTEIN_DOCKING_HOME

_logo = ""
_references = ['']
__version__ = '0.0.1'


class Plugin(pwem.Plugin):
    _url = "https://github.com/scipion-em/scipion-protein-docking"
    _homeVar = PROTEIN_DOCKING_HOME
    _pathVars = [PROTEIN_DOCKING_HOME]
    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(PROTEIN_DOCKING_HOME, 'frodock3-3.12')

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch frodock3. """
        environ = Environ(os.environ)

        environ.update({
            'PATH': Plugin.getHome(),
            'LD_LIBRARY_PATH': (str.join(cls.getHome(), 'lib') + ":" +
                                cls.getHome()),
        }, position=Environ.BEGIN)

        return environ

    @classmethod
    def getProgram(cls, program):
        """ Return the program binary that will be used. """
        path = cls.getVar(PROTEIN_DOCKING_HOME)
        if os.path.exists(path):
            binary = os.path.join(path, 'bin', program)

        return binary

    @classmethod
    def runProgram(cls, program, args):
        #runJob(None, program, args, env=cls.getEnviron())
        cmd = '%s %s' % (program, args)
        print("** Running command: %s" % greenStr(cmd), flush=True)
        os.system(cmd)


    @classmethod
    def defineBinaries(cls, env):
        # Add frodock3
        env.addPackage('frodock3', version='3.12',
                       tar='frodock3_linux64.tgz',
                       default=True)