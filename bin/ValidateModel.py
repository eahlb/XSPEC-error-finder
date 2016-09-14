from Sim import *
from Singleton import *
import json
import matplotlib.pyplot as plt
import numpy as np
import pyfits
import os
import shutil
import random
import sys

class ValidateModel(Singleton):
## Manager object for xspec model validation and visualization of data
##      (This is a singleton object and only one instance is allowed,
##      any two objects of this type will be the same object)
## Attributes
#   Outputs                     : list of SimOutput produced by the sims
## Methods
#   RunSimulations()            : run all queued simulations
#   AddSimulation(nr, pha_files, par, ign=[], type="Grid", silent=False)
#                               : add a Sim with arguments given to list of sims
#                                   to be run nr of times
#   Clear()                     : removes all sims and outputs
#   SaveData( path )            : saves Outputs to path
#   LoadData( path )            : load outputs from path
#   PlotParameter( p, path )    : plot parameter p and save plot to path
#   PlotAllParameters( path )   : plot all paramters to path
    def __init__(self):
        #### Attributes ####
        self.Outputs = []
        self.ID = ""
        #### Variables #####
        self._sim = None
        self._nr = 0
        #### Verify file structure ###
        if not os.path.isdir("data"): raise IOError("No data files exist")
        if not os.path.isdir("data/temp"): os.makedirs("data/temp")
        if not os.path.isdir("saves"): os.makedirs("saves")
    
    #### Methods ###########
    def Clear(self):
        del self.Outputs[:]
        if self._sim != None:
            # remove temp bak file
            tbak = []
            tpha = []
            for p in self._sim._pha:
                if "data/temp" in p: tpha.append( p )
                b = self._getBakPath(p)
                if "data/temp/" in b: tbak.append( b )
            for t in tbak+tpha: os.remove(t)
            self._sim.Clear()
        self._sim = None
        self._nr = 0
    
    def GetSimList(self):
        # get list of all possible sims
        x = []
        for file in os.listdir("data"):
            if file.endswith(".ds"):
                x.append(file[:-3])
        return x
    
    def RunSimSeries(self, dataID, saveID, nrOfSims=1, exposure=-1, bakExposure=-1, useBak=True, noWrite=False, silent=True, prompt=True, useStats=True, clear=True, pro=None, SNR=None, box=None, eb=[], L=[], startindex=0, randomStart=False, error=False, model=None):
        i = 0
        while os.path.isfile( "data/%s%s.ds" % (dataID,i+1) ): i += 1
        modelSNR = []
        correctSNR = []
        if SNR != None:
            saveID = "%sSNR%s_" % (saveID,SNR)
        for j in range(i):
            (x,y) = self.LoadSimulation(ID="%s%s"%(dataID,j+1), nrOfSims=nrOfSims, exposure=exposure, bakExposure=bakExposure, useBak=useBak, SNR=SNR, model=model)
            modelSNR.append(x)
            correctSNR.append(y)
            if j+1<startindex: continue
            self.RunSimulations(ID="%s%s"%(saveID,j+1), noWrite=noWrite, silent=silent, prompt=prompt, useStats=useStats, clear=clear, randomStart=randomStart, error=error)
        # save snr
        with open("saves/%s_SNR"%saveID,'w') as f:
            f.write(str(modelSNR) + "#" + str(correctSNR))
            f.close()
        if pro != None:
            pro(saveID, nrOfSims=i, savePlot=True, plotSteppar=False, dsPrefix=dataID, stepAll=False, useFak=True)
        if box != None:
            for e in eb:
                for l in L:
                    box(saveID,e,l)
    
    def RunSimulations(self, ID="last", noWrite=True, silent=True, prompt=True, useStats=True, clear=True, randomStart=False, error=False):
        # run all queued simulations
        if ID == "":
            raise AttributeError("An ID must be given")
        # check if there is a sim loaded
        if self._sim == None:
            raise IOError("No simulation loaded")
        self.ID = ID
        if os.path.isdir( "saves/{0}".format(self.ID) ):
            # check if user want to overwrite data
            ans = ""
            while not ans in ["y", "n"] and prompt:
                ans = raw_input("Sim with ID={0} alredy exist. Overwrite (y/n)?: ".format(ID))
            if ans == "n":
                print("Simulation aborted")
                return
            # remove old folder if it exist
            shutil.rmtree( "saves/{0}".format(self.ID) )
        # create new folder
        os.makedirs( "saves/{0}".format(self.ID) )
        
        # clear list of outputs
        del self.Outputs[:]
        print("----    Simulations started: {0}         ----".format(ID))
        # run sim the right amount of times
        self.printProgress(0, self._nr)
        for j in range(0, self._nr):
            if noWrite:
                self.Outputs.append(self._sim.Run( nr=j+1, silent=silent, useStats=useStats, clear=clear, randomStart=randomStart, error=error ))
            else:
                self.Outputs.append(self._sim.Run( "saves/{0}".format( ID), nr=j+1, silent=silent, useStats=useStats, clear=clear, randomStart=randomStart, error=error ))
            self.printProgress(j+1, self._nr)
        print("")
        self.SaveData()
        if clear: self.Clear()
        print("----    All simulations complete: {0}    ----".format(ID))
    
    def LoadSimulation(self, ID="sim", nrOfSims=1, exposure=-1, bakExposure=-1, useBak=True, SNR=None, model=None, dataFolder="data"):
        self.Clear()
        # background exposure should be the same as source exposure if set to default
        if bakExposure == -1:
            bakExposure=exposure
        # get settings
        _nrSpec, _type, _param, _pha, _rsp, _bak, _ign = self._readDs(dataFolder, ID)
        
        ### calculate model SNR
        # create temporary pha files
        _tempPha = []
        for (i,p) in enumerate(_pha):
            dst = "data/temp/%sdet%sSNR%s.pha" % (ID,i,SNR)
            while os.path.isfile(dst): dst = dst[:-4] + str(np.random.randint(9)) + dst[-4:]
            shutil.copy(p,dst)
            _tempPha.append(dst)
        # change paths so that the correct files are used
        self._editPaths(_tempPha,_rsp,_bak)
        
        obsSNR = self.GetSNR(_tempPha,_bak)[0]
        if SNR != None:
            temp = [0.0]*len(obsSNR)
            temp[obsSNR.index(max(obsSNR))] = float(SNR)
            obsSNR = temp
        self._addSimulation(1, _tempPha, _param, ign=_ign, type=_type, exposure=-1, bakExposure=-1, useBak=True, model=model)
        self.RunSimulations(ID="temp_%s"%ID, noWrite=False, silent=True, prompt=False,useStats=False,clear=True)
        
        _fak = []
        _bfak = []
        abc = "abcdefghi"
        for (i,x) in enumerate(_pha):
            _fak.append( "saves/temp_%s/spec1%s.fak"%(ID,abc[i]) )
            _bfak.append( "saves/temp_%s/spec1%s_bkg.fak"%(ID,abc[i]) )
        (calcSNR, calcBak) = self.GetSNR(_fak,_bfak, noStats=True)
        shutil.rmtree("saves/temp_%s"%ID)
        index = obsSNR.index(max(obsSNR))
        tm = self._getExpTime(_pha[index])
        
        ### modify SNR to desired level
        if SNR == None: x=1.0
        else:           x = max( [(calcSNR[index]/SNR)**2, 0.01])
        corrSNR = calcSNR[index]/np.sqrt(x)
        if SNR!=None and corrSNR*1.01 < SNR: nrOfSims=10
        # x now hold correction factor for background
        
        ### create temporary files
        # bak files
        _newBak = []
        for (i,b) in enumerate(_bak):
            dst = "data/temp/%sdet%sSNR%s.bak" % (ID,i,SNR)
            while os.path.isfile(dst): dst = dst[:-4] + str(np.random.randint(9)) + dst[-4:]
            shutil.copy(b,dst)
            _newBak.append(dst)
        # reduce background
        self._setBackgroundCorrection(x, _newBak)
        # pha files
        _newPha = []
        for (i,p) in enumerate(_pha):
            dst = "data/temp/%sdet%sSNR%s.pha" % (ID,i,SNR)
            while os.path.isfile(dst): dst = dst[:-4] + str(np.random.randint(9)) + dst[-4:]
            shutil.copy(p,dst)
            _newPha.append(dst)
        # change paths so that the temporary files are used
        self._editPaths(_newPha,_rsp,_newBak)
        _pha = _newPha
        
        print "----    Model predicted SNR: %.2f    Corrected SNR: %.2f  (Norm->%.2f, t->%.2f, B->%.2f )  ----" % (calcSNR[index], (corrSNR if SNR!=None else calcSNR[index]), _param[-1], (tm if exposure==-1 else exposure),x)
        # make sim ready for running
        self._addSimulation(nrOfSims, _pha, _param, ign=_ign, type=_type, exposure=exposure, bakExposure=bakExposure, useBak=useBak, model=model)
        return (calcSNR[index], (corrSNR if SNR!=None else calcSNR[index]))

    def _readDs(self, folder, ID):
        _pha = []
        _rsp = []
        _bak = []
        _ign = []
        _param = []
        with open( "data/%s.ds"%(ID), 'r' ) as f:
            s = f.readline().rstrip() # nr of spectra
            _nrSpec = int(s)
            s = f.readline().rstrip() # type of sim
            _type = s
            s = f.readline().rstrip() # model parameters
            _param = json.loads(s)
            s = f.readline().rstrip() # header break (##)
            if s != "##":
                raise IOError("Attempt to load incorrect file format (header)")
            for i in range(0,_nrSpec):
                s = f.readline().rstrip() # pha
                _pha.append( "{0}/{1}".format(folder,s) )
                s = f.readline().rstrip() # rsp
                _rsp.append( "{0}/{1}".format(folder,s) )
                s = f.readline().rstrip() # bak
                _bak.append( "{0}/{1}".format(folder,s) )
                s = f.readline().rstrip() # channel selection
                _ign.append( s )
                s = f.readline().rstrip() # spectrum break (#)
                if s != "#":
                    raise IOError("Attempt to load incorrect file format (spectra nr:{0})".format(i))
        f.close()
        return _nrSpec, _type, _param, _pha, _rsp, _bak, _ign
    
    def GetSteppar(self, prefix, dsPrefix, simNr, parameter, pBounds=[], logStep=False, nrSteps=100, silent=True, useFak=False, fakIndex=None, useOldModel=False):
        # load sim
        self.Clear()
        ID = "%s%s" % (dsPrefix, simNr)
        _nrSpec, _type, _param, _pha, _rsp, _bak, _ign = self._readDs("data", ID)
        self._addSimulation(1 , _pha, _param, _ign, "Grid", -1, -1, True, useOldModel=useOldModel)
        if useFak and fakIndex!=None:
            folder = "saves/%s%s" % (prefix, simNr)
            abc = "abcdefg"
            faks = []
            for (i,x) in enumerate(_pha):
                faks.append( "%s/spec%s%s.fak" % (folder, fakIndex, abc[i]) )
            self._sim.LoadFak(faks)
        else:
            # run one sim without statistical noise
            tempID = "temp_%s_%s_%s" % (dsPrefix, simNr, parameter)
            # make sure tempID is unique
            while os.path.isdir( "saves/{0}".format(tempID) ):
                tempID += "_%s" % str( random.randint(0,9) )
            self.RunSimulations(ID=tempID, noWrite=True, silent=True ,prompt=False, useStats=False, clear=False)
            # delete tempID folder
            shutil.rmtree( "saves/{0}".format(tempID) )
        # get steppar
        res = self._sim.GetSteppar(parameter, pBounds, logStep, nrSteps,silent)
        self._sim.Clear()
        return res

    def _setBackgroundCorrection(self, x, bakFiles):
        for b in bakFiles:
            with pyfits.open(b,"update") as hd:
                for i in range( len(hd["SPECTRUM"].data) ):
                    hd["SPECTRUM"].data[i]["RATE"] = hd["SPECTRUM"].data[i]["RATE"]*x
                    hd["SPECTRUM"].data[i]["STAT_ERR"] = hd["SPECTRUM"].data[i]["STAT_ERR"]*x

    def _editPaths(self, phaFiles, rspFiles, bakFiles):
        i = 0
        for p in phaFiles:
            with pyfits.open(p,"update") as hd:
                hd[1].header["RESPFILE"] = rspFiles[i]
                hd[1].header["BACKFILE"] = bakFiles[i]
                i += 1
    
    def _getRspPath(self, phaPath):
        x = ""
        with pyfits.open(phaPath) as hd:
            x = hd[1].header["RESPFILE"]
        return x
    
    def _getBakPath(self, phaPath):
        x = ""
        with pyfits.open(phaPath) as hd:
            x = hd[1].header["BACKFILE"]
        return x
        
    def _getExpTime(self, phaPath):
        x = ""
        with pyfits.open(phaPath) as hd:
            x = hd[1].header["EXPOSURE"]
        return x

    def _addSimulation(self, nr, phaFiles, par, ign, type, exposure, bakExposure, useBak, model=None):
    # add sim to list of sims to be run
        if type == "Grid":
            if model==None: model = 'atable{models/TableModel_grid1296_r0free_res8_Nr8.fits}'
            self._sim = GridSim(phaFiles, par, ign, model=model, exposure=exposure, bakExposure=bakExposure, useBak=useBak)
        elif "." in type:
            # assume off-grid sim and that type is path to sim model
            if model==None: model = 'atable{models/TableModel_grid1296_r0free_res8_Nr8.fits}'
            self._sim = OffGridSim(type, phaFiles, par, ign, model, exposure, bakExposure, useBak)
        else:
            self._sim = NativeSim(phaFiles, par, ign, model=type, exposure=exposure, bakExposure=bakExposure, useBak=useBak)
        self._nr = int(nr)

    def GetSNR(self, files, bfiles, noStats=False):
        SNR = []
        BackRate = []
        for (j,file) in enumerate(files):
            signal = [] # cts/s
            background = [] #cst/s
            bvar = [] # cts/s
            t = 0.0 #s
            with pyfits.open(file) as hd:
                t = float( hd["SPECTRUM"].header["EXPOSURE"] )
                for i in range( 0, len(hd["SPECTRUM"].data) ):
                    if noStats: signal.append( hd["SPECTRUM"].data[i]["RATE"] )
                    else: signal.append( hd["SPECTRUM"].data[i]["COUNTS"]/t )
                hd.close()
            with pyfits.open(bfiles[j]) as hd:
                for i in range( 0, len(hd["SPECTRUM"].data) ):
                    background.append( hd["SPECTRUM"].data[i]["RATE"] )
                    bvar.append( hd["SPECTRUM"].data[i]["STAT_ERR"]**2 )
                hd.close()
            s = np.array(signal)
            b = np.array(background)
            S = (np.sum(s) - np.sum(b))*t
            B = np.sum(b)*t
            X = ( S )/( S+B )
            SNR.append( S / np.sqrt(B) )
            BackRate.append( [np.sum(b), np.sqrt(np.sum(bvar))] )
        #print "%s SNR: %.1f  X: %.3f" % (index, SNR, X)
        return (SNR, BackRate)
    
    def printProgress (self, iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 70):
        # Adopted from https://gist.github.com/vladignatyev/06860ec2040cb497f0f3
        """
        Call in a loop to create terminal progress bar
        @params:
        iterations  - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        """
        filledLength    = int(round(barLength * iteration / float(total)))
        percents        = round(100.00 * (iteration / float(total)), decimals)
        bar             = '#' * filledLength + '-' * (barLength - filledLength)
        sys.stdout.write('%s [%s] %s%s %s\r' % (prefix, bar, percents, '%', suffix)),
        sys.stdout.flush()
        if iteration == total:
            print("\n")

    def SaveData(self):
        path = "saves/{0}/output.do".format(self.ID)
        if not os.path.isdir( "saves/{0}".format(self.ID) ):
            os.makedirs( "saves/{0}".format(self.ID) )
        # open file
        with open(path,'w') as f:
        # save header
        # save Outputs
            for out in self.Outputs:
                f.write(str(out) + "\n")
        f.close()
        # add EBOUNDS to .fak
        print("----    All output objects saved: {0}    ----".format(self.ID))

    def LoadData(self, ID):
        # open file
        ##  Format
        #   [model parameters]
        #   [fit parameters]
        #   [parameter names]
        #   PG stat
        #   Chi^2
        #   nr. dof
        #   #
        self.ID = ID
        path = "saves/{0}/output.do".format(self.ID)
        with open(path, 'r') as f:
        # clear current list of Outputs
            del self.Outputs[:]
        # load header
        # load Outputs
            out = SimOutput()
            s = f.readline().rstrip()
            while True:
                # 1st - check if at eof
                # model parameters
                if s == "":
                    break
                out.ModelParameters = json.loads(s)
                # 2nd
                # fit parameters
                s = f.readline().rstrip()
                out.FitParameters = json.loads(s)
                # 3rd
                # parameter names
                s = f.readline().rstrip()
                out.ParameterNames = json.loads(s)
                # 4th
                # PG stat
                s = f.readline().rstrip()
                out.PGstat = json.loads(s)
                # 5th
                # Chi2
                s = f.readline().rstrip()
                out.Chi2 = json.loads(s)
                # 6th
                # nr dof
                s = f.readline().rstrip()
                out.Dof = json.loads(s)
                # end sign - #
                s = f.readline().rstrip()
                if s != "#":
                    raise IOError("Attempt to load incorrect file format")
                # read next line
                s = f.readline().rstrip()
                self.Outputs.append(out)
                out = SimOutput()
        f.close()
        print("----    Data loaded:{0}    ----".format(ID))
        return self.Outputs








