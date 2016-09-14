import xspec
from SimOutput import *
import copy
import os
import shutil
import re
import numpy as np

class Sim:
## Handler class for simulating a fake spectrum and fitting it to a model
## Attributes
#   ModelParameters     : list of all parameters in the current model
#   Spectra             : list of all loaded spectra
#   BadChannels         : list(str) of channels to be ignored for each spectra
## Methods
#   Run() -> SimOutput  : run the simulation, return a (copy) SimOutput object
#   DefineSimModel()    : define the model used for simulating fake spectra
#                           must be overriden or will throw error
#   DefineFitModel()    : define the model used for fitting the fake spectra
#   _initParameters()    : get all fit parameters for the model and place the in ModelParameters
    def __init__(self, phaFiles, par, bad, model, exposure=-1, bakExposure=-1, useBak=True):
        #xspec = __import__("xspec")
        #### Attributes ####
        self.ModelParameters = []
        self.ParameterNames = []
        self.BadChannels = bad
        self.Exposure = float(exposure)
        self.BakExposure = float(bakExposure)
        self.UseBak = useBak
        #### Variables #####
        self._pha = phaFiles
        self.spec_model = model
        self._out = SimOutput()
        self._parmeter_values = par
        #### Constructor ###
        ## chatter settings
        xspec.Xset.logChatter = 0
        ## fit settings
        xspec.Fit.query = "no"
        xspec.Fit.statMethod = "pgstat"
        xspec.Fit.statTest = "chi"
        xspec.Fit.nIterations = 10000
        ## get model parameters for output
        self._out.ModelParameters = par

#### Methods ###########
    def DefineSimModel(self):
    # must be overridden by subclass or will throw error
        raise NotImplementedError("Simulation model not defined")
    
    def DefineFitModel(self):
        xspec.AllModels.clear()
        # define the model
        if self.spec_model != None:
            self._model = xspec.Model( self.spec_model )
        else:
            raise RuntimeError("No fit model selected")
        # set parameters
        self._initParameters()
        for (i,p) in enumerate(self.ModelParameters):
            if p!=None: p.values = self._parmeter_values[i]
        # check the number of parameters
        if len(self._parmeter_values) != len(self.ModelParameters):
            raise RuntimeError("Incorrect number of model parameters passed")
    
    def _initParameters(self):
        # must be overridden by subclass or will throw error
        raise NotImplementedError("Simulation model not defined")

    def Clear(self):
        xspec.AllData.clear()
        xspec.AllModels.clear()

    def Run(self, folder="", nr=0, silent=True, useStats=True, clear=True, plot=False, randomStart=False, error=False):
        if silent:
            xspec.Xset.chatter = 0
        # load spectras with background and response
        for pha in self._pha:
            s = xspec.Spectrum(pha)
            if not self.UseBak:
                s.background = None
        # fakeit
        self._doFake(folder, nr, useStats)
        # fit data
        self._doFit(randomStart, error)
        # computation complete
        if not silent:
            xspec.AllData.show()
            xspec.AllModels.show()
        if plot:
            xspec.Plot.device = "/xs"
            xspec.Plot.xAxis = "keV"
            xspec.Plot("eeufspec eemodel")
        if clear:
            # clear loaded spectra
            self.Clear()
        # restore chatter level
        xspec.Xset.chatter = 10
        # return the data (by value)
        return copy.deepcopy(self._out)

    def GetSteppar(self, par, pBounds, logStep=False, nrSteps=100, silent=True):
        arg = ""
        if silent: xspec.Xset.chatter = 0
        if logStep:
            arg += "log "
        if self._model.nParameters == 6:
            if par>2: par -= 1
            elif par == 2: return [[0]*nrSteps]
        arg += "%.0f %f %f %.0f" % (par+1, pBounds[0], pBounds[1], nrSteps)
        print("----    Running steppar: p={0}         ----".format(par))
        xspec.Fit.steppar( arg )
        xspec.Xset.chatter = 10
        return xspec.Fit.stepparResults("delstat")

    def LoadFak(self, faks, silent=True):
        if silent: xspec.Xset.chatter = 0
        print faks
        self.Clear()
        for fak in faks:
            s = xspec.Spectrum(fak)
        self._doFit()
        xspec.Xset.chatter = 10

    def _doFake(self, folder="", nr=0, useStats=False):
        # load the model used
        self.DefineSimModel()
        # no spectra are produced if no folder is provided
        if folder != "":
            nw = False
        else:
            nw = True
        fs = []
        i = 0
        abc = "abcdefghijklm"
        for x in self._pha:
            # set the full .fak path
            pt = "{0}/spec{1}{2}.fak".format(folder, nr, abc[i])
            fs.append( xspec.FakeitSettings(fileName=pt) )
            if self.Exposure != -1:
                fs[-1].exposure = str(self.Exposure)
            if self.BakExposure != -1:
                fs[-1].backExposure = str(self.BakExposure)
            i += 1
        xspec.AllData.fakeit(nSpectra=len(self._pha), settings=fs, applyStats=useStats, noWrite=nw)
    
    def _doFit(self, randomStart=False, error=False):
        # load model for fitting
        self.DefineFitModel()
        # ignore bad channels
        i = 1
        for bc in self.BadChannels:
            xspec.AllData.ignore( "{0}: {1}".format(str(i), bc) )
            i+=1
        # select inital start guess
        if randomStart:
            # select parameters to change
            map = self._getMap()
            # modify each free parameter
            for i in map:
                # find allowed spectrum for the paramter
                lim = [self.ModelParameters[i].values[2], self.ModelParameters[i].values[5]]
                # set parameter value
                rval = np.random.rand()*(lim[1]-lim[0]) + lim[0]
                self.ModelParameters[i].values = rval
        # fit spectrum to model
        xspec.Fit.perform()
        # do error
        if error: self._doError()
        # save data in a SimOutput object
        del self._out.FitParameters[:]
        del self._out.ParameterNames[:]
        for p in self.ParameterNames:
            self._out.ParameterNames.append(p)
        for (i,p) in enumerate(self.ModelParameters):
            if p!=None: self._out.FitParameters.append(p.values[0])
            else: self._out.FitParameters.append(self._parmeter_values[i])
        self._out.PGstat = xspec.Fit.statistic
        self._out.Chi2 = xspec.Fit.testStatistic
        self._out.Dof = xspec.Fit.dof
        
    def _doError(self):
        # select parameters to change
        map = self._getMap()
        oldq = xspec.Fit.query
        xspec.Fit.query = "yes"
        for i,m in enumerate(map):
            xspec.Fit.error("1.0 %s"%(i+1))
            # the error string should be "FFFFFFFFF"
        xspec.Fit.query = oldq

    def _getMap(self):
        # return parameters indices for all non-frozen parameters
        map = []
        for i,p in enumerate(self.ModelParameters):
            if p!=None and not p.frozen:
                map.append(i)
        return map

class GridSim(Sim):
## Simulation using the same model for simulating fake spectra and fitting
## Methods
#   DefineSimModel() : define the model used for producing fake spectra
    def __init__(self, pha_files, par, ign, model, exposure, bakExposure, useBak):
        Sim.__init__(self, pha_files, par, ign, model, exposure, bakExposure)
#### Methods ###########
    def DefineSimModel(self):
        self.DefineFitModel()
    
    def DefineFitModel(self):
        Sim.DefineFitModel(self)
        # freeze the norm
        self._model.Photosphere.norm.frozen = True

    def _initParameters(self):
        del self.ModelParameters[:]
        del self.ParameterNames[:]
        # get parameters
        p1 = self._model.Photosphere.Tau
        self.ModelParameters.append(p1)
        p2 = self._model.Photosphere.epsilon_e
        self.ModelParameters.append(p2)
        if self._model.nParameters==7: p3 = self._model.Photosphere.epsilon_B
        else: p3 = None
        self.ModelParameters.append(p3)
        p4 = self._model.Photosphere.LGRB
        self.ModelParameters.append(p4)
        p5 = self._model.Photosphere.Gamma
        self.ModelParameters.append(p5)
        p6 = self._model.Photosphere.z
        self.ModelParameters.append(p6)
        p7 = self._model.Photosphere.norm
        self.ModelParameters.append(p7)
        for p in self.ModelParameters:
            if p!=None: self.ParameterNames.append(p.name)
            else: self.ParameterNames.append("-")

class OffGridSim(GridSim):
## Simulation using different models for simulating fake spectra and fitting
## Methods
#   DefineSimModel() : define the model used for producing fake spectra
    def __init__(self, simMod, pha_files, par, ign, model, exposure, bakExposure, useBak):
        GridSim.__init__(self, pha_files, par, ign, model, exposure, bakExposure, useBak)
        self.simulationModel = simMod
    #### Methods ###########
    def DefineSimModel(self):
        xspec.AllModels.clear()
        p = re.search("(.*)lum(.*)", self.simulationModel.split("/")[-1])
        mod = "models/temp_lum" + p.group(2)
        while os.path.isfile(mod): mod = mod[:-5] + str(np.random.randint(9)) + mod[-5:]
        shutil.copy(self.simulationModel, mod)
        self._model = xspec.Model( "atable{ %s }"%mod )
        os.remove(mod)
        # set parameters
        self._initParameters()
        # set redshift
        self._model.Photosphere.z = self._parmeter_values[5]
        # set and freeze the norm
        self._model.Photosphere.norm = self._parmeter_values[-1]
        self._model.Photosphere.norm.frozen = True


class NativeSim(Sim):
    def __init__(self, pha_files, par, ign, model, exposure, bakExposure, useBak):
        if model==None:
            raise RuntimeError("A native XSPEC model must be specified")
        Sim.__init__(self, pha_files, par, ign, model, exposure, bakExposure)
#### Methods ###########
    def DefineSimModel(self):
        self.DefineFitModel()
    def _initParameters(self):
        del self.ModelParameters[:]
        del self.ParameterNames[:]
        for i in range(self._model.nParameters):
            p = self._model(i+1)
            self.ModelParameters.append(p)
            self.ParameterNames.append(p.name)
        xspec.AllModels.show()



