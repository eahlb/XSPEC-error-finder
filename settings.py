import os
import re
import pyfits
import shutil
import random
import sys
import string
import bisect
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib import rc
import matplotlib.pylab as lab
from pylab import *
from scipy.integrate import quad

def PlotCountSpectra(folder="forftp/res_8/", showPlot=False):
    path = "/Users/eahlb/coding/python lib/%s" % folder
    name = []
    spectra = []
    ebins = []
    ym = 0
    plt.figure(figsize=(18, 10))
    for file in os.listdir(path):
        if file.endswith(".txt"):
            #p = re.search("(.*)lum(.*)gamma(.*)tau(.*)epl(.*)ee(.*)eb(.*)ed(.*)", file)
            #tau = p.group(4); ee = p.group(6); eb = p.group(7); L = p.group(2); G = p.group(3)
            #if float(tau) != 15.0: continue
            name.append( file[:-4] )
            file = "%s%s" % (path,file)
            with open(file) as f:
                finals = f.read().split("Final spectrum")[1]
                ts = []
                te = []
                for line in finals.split("\n")[3:]:
                    # add to spectra
                    tl = line.split()
                    if len(tl)==0: break
                    te.append( float(tl[0]) )
                    ts.append( float(tl[1]) )
                nts = np.array(ts)
                nte = np.array(te)
                nts,nte = _getBins(nts,nte)
                spectra.append(nts)
                ebins.append(nte)
                ym = max(ym, np.amax(nts))
                plt.plot(nte,nts, label=name[-1][8:-16])
    plt.xscale('log')
    plt.xlabel(r"Energy (keV)")
    plt.xlim([1e-3, 1e8])
    plt.yscale('log')
    plt.ylabel(r"EF$_E$ (keV s$^{-1}$ m$^{-2}$)")
    plt.ylim([ym*1e-4, ym*1e1])
    plt.grid()
    plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., prop={'size':8})
    if showPlot: plt.show()
    else: plt.savefig("/Users/eahlb/coding/python lib/models/lastsim.png", bbox_inches='tight')
def _getBins(N_g, E_g, tdyn=1, z=1, gamma=100):
    c = 3e10 # cm s^-1
    ergtokev = 6.24*10**8 # Go from erg to KeV
    H0 = 69.9 #km/sMpc
    OmegaM = 0.286
    OmegaLambda = 0.714
    OmegaK = 0 # 1-OmegaM-OmegaLambda
    dH = divide(c,H0)*1e-5 #answer in Mpc
    Einv = lambda Z: divide(1,sqrt(OmegaM*(1+Z)**3+OmegaK*(1+Z)**2+OmegaLambda))
    dL = lambda z: multiply((1+z),multiply(dH,quad(Einv,0,z)[0])) # obtain value in Mpc
    Mpc_to_cm = 3.086e24
    d_L = multiply(dL(z),Mpc_to_cm)
    # Next, we assume that E_g and N_g are read from the output files from the simulations.
    # The constant "tdyn" can also b found in the simulation files.
    E_g2 = multiply(E_g,50) # Here we boost the energies with gamma
    energybins_lowerbound = E_g2
    energybins_upperbound = E_g2[1:].tolist()
    energybins_upperbound.append(multiply(E_g2[-1],divide(E_g2[1],E_g2[0])))
    energybins_upperbound = np.array(energybins_upperbound)
    d_E = np.subtract(energybins_upperbound,energybins_lowerbound)
    counts = N_g
    area = 4.0*np.pi*d_L**2
    counts = divide(counts,area*tdyn)
    shiftfactor = gamma/50
    newFlux = []
    j=len(energybins_lowerbound)-1
    newEnergy = energybins_lowerbound[j]/shiftfactor
    newFlux = np.zeros(len(energybins_lowerbound))
    while j>0 and newEnergy>energybins_lowerbound[0]:
        k=0
        while energybins_lowerbound[k]<newEnergy:
            k+=1
        k-=1
        flux1 = counts[k]
        flux2 = counts[k+1]
        energy1 = energybins_lowerbound[k]
        energy2 = energybins_lowerbound[k+1]
        flux = abs((flux2-flux1)/(energy2-energy1)*(newEnergy-energy1)+flux1)
        newFlux[j]=flux
        j-=1
        newEnergy = energybins_lowerbound[j]/shiftfactor
    energies = np.zeros(len(energybins_lowerbound))
    for i in range(0,len(energybins_lowerbound)):
        energies[i] = sqrt(multiply(energybins_upperbound[i],energybins_lowerbound[i]))
    x_plot = np.array([a*ergtokev/(1.+z) for a in energies])
    y_plot = np.array(newFlux)
    E_over_dE = energies/d_E
    EF_E = x_plot*y_plot*E_over_dE*gamma/(1.+z)
    return EF_E, x_plot

def CreateSims(grb, prefix, eb):
    #prefix = "eb1e-6_L1e2_"
    ## type of sim
    type  = "Grid"
    ## parameters to span
    '''tau   = [1.0, 5, 10, 35]
    ee    = [1e-2, 0.05, 0.1, 0.2, 0.3, 0.5]
    #eb    = [1e-6]
    L     = [1,10,100]
    gamma = [100, 150, 200, 250, 300, 400]'''
    tau   = [1.0, 5, ]
    ee    = [1e-2,  0.5]
    L     = [10, 100]
    gamma = [100, 250, 400]
    ## files to use
    '''pha   = ["grb100414097_bgo_01_int11.pha", "grb100414097_nai_07_int11.pha", "grb100414097_nai_09_int11.pha", "grb100414097_nai_11_int11.pha"]
    rsp   = ["grb100414097_bgo_01_int11.rsp", "grb100414097_nai_07_int11.rsp", "grb100414097_nai_09_int11.rsp", "grb100414097_nai_11_int11.rsp"]
    bak   = ["grb100414097_bgo_01_int11.bak", "grb100414097_nai_07_int11.bak", "grb100414097_nai_09_int11.bak", "grb100414097_nai_11_int11.bak"]'''
    ## ignore settings
    ign   = ["**-200.0 40000.0-**", "**-8.0 1000.0-**", "**-8.0 1000.0-**", "**-8.0 1000.0-**"]
    # create list of files from grb
    (pha, rsp, bak, ign) = _getGrbFiles(grb)
    nr    = len(pha)
    # check that the data folder exsits
    if not os.path.isdir( "data" ):
        raise IOError("no data folder")
    # check if the same prefix exists
    for file in os.listdir("data"):
        if file.endswith(".index"):
            if file == "{0}.index".format(prefix):
                # simulations exist, overwrite?
                ans = ""
                while not ans in ["y", "n"]:
                    ans = raw_input("Sims with prefix={0} alredy exist. Overwrite (y/n)?: ".format(prefix))
                if ans == "n":
                    raise IOError("Will not overwrite files".format(prefix))
    # create list of parameters
    params = []
    for _tau in tau:
        for _ee in ee:
            for _eb in eb:
                for _L in L:
                    for _gamma in gamma:
                        mp = [_tau, _ee, _eb, _L, _gamma, 1, 1]
                        params.append(mp)
    print(len(params))
    # create .ds files
    i = 0
    for mp in params:
        i += 1
        with open( "data/{0}{1}.ds".format(prefix, i), 'w' ) as f:
            f.write( "{0} \n".format(nr) )       # nr of spectra
            f.write( "{0} \n".format(type) )    # type of sim
            f.write( "{0} \n".format(mp) )      # model parameters
            f.write( "## \n" )                  # break
            x = 0
            for p in pha:
                f.write( "{0} \n".format(pha[x]) )
                f.write( "{0} \n".format(rsp[x]) )
                f.write( "{0} \n".format(bak[x]) )
                f.write( "{0} \n".format(ign[x]) )
                f.write( "# \n" )
                x += 1
        f.close()
    print(i)
    # create index file
    i = 0
    with open( "data/{0}.index".format(prefix) ,'w') as f:
        # nr of spectra
        f.write( "{0} simulations with prefix: {1} \n".format(len(params), prefix) )
        # type of simulation
        f.write( "Simulation type: {0} \n".format(type) )
        # break
        f.write( "## \n" )
        for mp in params:
            i += 1
            f.write( "{0} - {1} \n".format(i, mp) )
    f.close()

def _getGrbFiles(grb):
    ## ignore settings [b,n]
    ebounds   = ["**-200.0 40000.0-**", "**-8.0 1000.0-**"]
    pha = []; rsp = []; bak = []; ign = []
    for file in os.listdir("data/%s"%grb):
        if file.endswith(".pha"):
            stem = "%s/%s" % (grb, file[:-4])
            pha.append("%s.pha"%stem)
            if os.path.isfile("data/%s.bak"%stem):
                bak.append("%s.bak"%stem)
            else: raise  IOError("No background file exist")
            if os.path.isfile("data/%s.rsp"%stem):
                rsp.append("%s.rsp"%stem)
            elif os.path.isfile("data/%s.rsp2"%stem):
                rsp.append("%s.rsp2"%stem)
            else: raise  IOError("No response file exist")
            if file[0] == 'b':
                ign.append(ebounds[0])
            elif file[0] == 'n':
                ign.append(ebounds[1])
            else: raise  IOError("Unknown detector type")
    return pha, rsp, bak, ign

def CreateOffGrid(grb, prefix, grid): #SimID, grb):
    modFolder = "models/%s" % grid

    # check that the data folder exsits
    if not os.path.isdir( "data/%s"%grb ):
        raise IOError("no data folder")
    # check if the same prefix exists
    for file in os.listdir("data"):
        if file.endswith(".index"):
            if file == "{0}.index".format(prefix):
                # simulations exist, overwrite?
                ans = ""
                while not ans in ["y", "n"]:
                    ans = raw_input("Sims with prefix={0} alredy exist. Overwrite (y/n)?: ".format(prefix))
                if ans == "n":
                    raise IOError("Will not overwrite files".format(prefix))
    # create list of files from grb
    (pha, rsp, bak, ign) = _getGrbFiles(grb)
    # create list of models
    # create list of parameters
    models = []
    params = []
    for file in os.listdir( modFolder ):
        if file.endswith(".fits"):
            models.append( "%s/%s" % (modFolder,file) )
            p = re.search("(.*)lum(.*)gamma(.*)tau(.*)epl(.*)ee(.*)eb(.*)ed(.*)", file)
            tau = p.group(4)
            ee = p.group(6)
            eb = p.group(7)
            L = p.group(2)
            G = p.group(3)
            mp = [float(tau), float(ee), float(eb), float(L), float(G), 1, 1]
            params.append(mp)
    # create .ds files
    i = 0
    nr = len(pha)
    for mp in params:
        i += 1
        with open( "data/{0}{1}.ds".format(prefix, i), 'w' ) as f:
            f.write( "{0} \n".format(nr) )       # nr of spectra
            f.write( "{0} \n".format(models[i-1]) )    # type of sim
            f.write( "{0} \n".format(mp) )      # model parameters
            f.write( "## \n" )                  # break
            x = 0
            for p in pha:
                f.write( "{0} \n".format(pha[x]) )
                f.write( "{0} \n".format(rsp[x]) )
                f.write( "{0} \n".format(bak[x]) )
                f.write( "{0} \n".format(ign[x]) )
                f.write( "# \n" )
                x += 1
        f.close()
    print(i)
    # create index file
    i = 0
    with open( "data/{0}.index".format(prefix) ,'w') as f:
        # nr of spectra
        f.write( "{0} simulations with prefix: {1} \n".format(len(params), prefix) )
        # type of simulation
        f.write( "Simulation type: {0} \n".format(type) )
        # break
        f.write( "## \n" )
        for mp in params:
            i += 1
            f.write( "{0} - {1} \n".format(i, mp) )
    f.close()

def CreateModels():
    filenames =[]
    for file in os.listdir("/Users/eahlb/coding/python lib/forftp/res_8"): filenames.append(file)
    for i,file in enumerate(filenames):
        file = file[:-4]
        #file2 = file
        #while file==file2: file2=filenames[random.randint(0, len(filenames)-1)]
        p = re.search("(.*)lum(.*)gamma(.*)tau(.*)epl(.*)ee(.*)eb(.*)ed(.*)", file)
        tau = p.group(4)
        ee = p.group(6)
        eb = p.group(7)
        L = p.group(2)
        G = p.group(3)
        pt = p.group(1) + "lum%4%gamma%5%tau%1%epl" + p.group(5) + "ee%2%eb%3%ed" + p.group(8)
        with open("/Users/eahlb/coding/python lib/forftp/inputs/InputData%s"%i,'w') as f:
            f.write('''NSTEPS_1 = 1
NSTEPS_2 = 1
NSTEPS_3 = 1
NSTEPS_4 = 1
NSTEPS_5 = 1

NPAR = 5

NENERGIES = 785
INT_MET = 1
N_ADDSPEC = 0
REDSHIFT = 1

PARNAME_1 = "Tau"
PARNAME_2 = "epsilon_e"
PARNAME_3 = "epsilon_B"
PARNAME_4 = "LGRB"
PARNAME_5 = "Gamma"

PAR1_1 = %s

PAR_FILENAME1_1 = "%s"

PAR2_1 = %s

PAR_FILENAME2_1 = "%s"

PAR3_1 = %s

PAR_FILENAME3_1 = "%s"

PAR4_1 = %s

PAR_FILENAME4_1 = "%s"

PAR5_1 = %s

PAR_FILENAME5_1 = "%s"

PATH = "/home/eahlb/forftp/res_8/%s.txt"

OUTPUT_FITS = "TableModel_%s"''' % (float(tau),tau, float(ee),ee, float(eb),eb, float(L),L, float(G),G, pt,file))
        f.close()

def findstartstop(openfile):
    start = "Final spectrum"
    slut = "Electron distribution"
    slutindex = 0
    with open(openfile,'r') as fil:
        for num, line in enumerate(fil,1):
            if start in line:
                startindex = num+2 # Since "Final spectrum" is always three rows from the first row of data
            if slut in line:
                slutindex = num-1 # Since "Electron distribution" is always one row from the end of the photon distribution
            sista = float(num)
            if not slutindex==0:
                skipfoot = int(sista-slutindex)
            else:
                skipfoot = 0
    return startindex,skipfoot,slutindex

def ChangeEres():
    PATH = "/Users/eahlb/coding/python lib/forftp/"  #############
    os.chdir(PATH)
    fileending = '(?i)(DTF\d+[\.]*[\d]*)*'

    #Define new energy grid with eight times the original resolution
    res_factor=8
    j=np.arange(785)
    new_estep = 1.0415
    elow=1.385e-14
    energy_new = elow * new_estep**j
    #*********************Create file names********************************
    #Rootnames for file
    # rootfilename = 'lumLgammaGtauTeplPeeEebBedD.txt'

    #Loop over all free parameters
    n=0
    filelist = []

    for number,file in enumerate(os.listdir(PATH),1):
        if (re.match('(r0_(.*))*lum(.*)gamma(.*)tau(.*)epl(.*)ee(.*)eb(.*)ed(\d+[\.]*[\d]*)'+fileending+'\.txt',file,re.IGNORECASE)):
            bisect.insort(filelist,file)
    res8folder = PATH+'res_8'
    if os.path.isdir(res8folder):
        shutil.rmtree(res8folder)
    os.mkdir(res8folder)
    for f in filelist:
        infile = f
        startindex,skipfoot,slutindex = findstartstop(f)
        outfile= re.sub('.txt','_res'+str(res_factor)+'.txt',infile)
        #count input files
        n=n+1
        #Read original file
        energy, nphot = genfromtxt(f,skiprows=startindex,skip_footer=skipfoot-1,usecols=(0,1),unpack=True)
        # Get some information we need
        elen = len(energy)
        #Deal with the bin that has the annihilation line
        nphot_withline = nphot[54]
        nphot[54] = 0.5*(nphot[53]+nphot[55])
        line_photons = nphot_withline - nphot[54]
        
        #interpolate
        ff= interpolate.interp1d(energy,nphot,kind='slinear')
        nphot_new = ff(energy_new)
        
        #Add the line photons at the right energy
        linedist = np.abs(energy_new - 8.187e-7)
        line_index = np.argmin(linedist)
        if (energy_new[line_index]>8.187e-7):
            line_index = line_index - 1
        nphot_new[line_index] = nphot_new[line_index] + line_photons

        #Remove any spurious negative values resulting from the interpolation
        min_nphot = np.min(nphot)
        for i in range (len(energy_new)):
            if (nphot_new[i]<min_nphot):
                nphot_new[i]=0.

        #Normalise so that the total number of photons is conserved
        tot_nphot = np.sum(nphot)
        tot_nphot_new = np.sum(nphot_new)
        norm= tot_nphot/tot_nphot_new
        nphot_new = nphot_new * norm
        
        
        #Copy input file and remove the old low-res sepectrum AND the electron spectrum (this won't be included again)
        copy = 'cp ' + infile + ' ttt.txt'
        os.system(copy)
        removeold = 'sed \''+str(startindex-3)+',$d\'' ' ttt.txt >' + outfile
        os.system(removeold)
        os.system('rm ttt.txt')
        
        #Append the interpolated data to the output file
        writefile = open(outfile, 'a') # open file for appending text
        #pad BB part with 0s so that the lines will be correct for the table code
        for j in range(len(energy_new)-len(energy)):
            writefile.write('%10.3e %10.3e %10.3e \n' % (0., 0., 0.))
        #now write new spectrum
        writefile.write('\nFinal spectrum\nphotons: E_g     N_g     dNg_dt \n-------------------------------------\n')
        for i in range(len(energy_new)):
            writefile.write('%10.3e %10.3e %10.3e \n' % (energy_new[i], nphot_new[i], 0.))
        writefile.close()
        # Move the files to another directory. Make sure that this directory exists
        shutil.move(outfile,res8folder)

def GetAllSNR(ID, ending=".fak", maxind=None, noStats=False, showBkg=False):
    x = []
    l = []
    path = "saves/%s" % ID
    for file in os.listdir(path):
        if file.endswith(ending):
            if ending == ".fak":
                if "_bkg" in file: continue
                p = re.search("spec(.*)%s"%ending, file)
                if not ( int(p.group(1)[:-1]) in x ): x.append(int(p.group(1)[:-1]))
                if not ( p.group(1)[-1] in l ): l.append(p.group(1)[-1])
            elif ending == ".pha":
                x.append(file)
                l.append(file[:-4]+".bak")
            else: break
    if maxind==None: maxind=len(x)
    SNR = []
    for (i,index) in enumerate(sorted(x)[:maxind]):
        printProgress(i,maxind)
        if ending==".fak": SNR.append(GetSNR(index,l,path,ending,noStats=noStats,showBkg=showBkg)[0])
        elif ending==".pha": SNR.append(GetSNR(index,[l[i]],path,ending,noStats=noStats,showBkg=showBkg)[0])
    print "Mean SNR: %.2f (+- %.2f)                               " % (np.mean(SNR), np.std(SNR))
    printProgress(1,1)
    #print "All SNR: %s" % SNR
    #plt.plot(SNR)
    #plt.ylim([0,70])
    #plt.grid()
    #plt.show()

def GetSNR(index,sub,path,ending=".fak", bending=".fak", noStats=False, showBkg=False):
    SNR = []
    BackRate = []
    for si in sub:
        signal = [] # cts/s
        background = [] #cts/s
        bvar = [] # cts/s
        t = 0.0 #s
        if ending == ".fak":
            file = "%s/spec%s%s.fak" % (path,index,si)
            bfile = "%s/spec%s%s_bkg.fak" % (path,index,si)
        elif ending == ".pha":
            file = "%s/%s" % (path,index)
            bfile = "%s/%s" % (path,si)
        else: return 0
        with pyfits.open(file) as hd:
            t = float( hd["SPECTRUM"].header["EXPOSURE"] )
            for i in range( 0, len(hd["SPECTRUM"].data) ):
                if noStats: signal.append( hd["SPECTRUM"].data[i]["RATE"] )
                else: signal.append( hd["SPECTRUM"].data[i]["COUNTS"]/t )
            hd.close()
        with pyfits.open(bfile) as hd:
            for i in range( 0, len(hd["SPECTRUM"].data) ):
                background.append( hd["SPECTRUM"].data[i]["RATE"] )
                bvar.append( hd["SPECTRUM"].data[i]["STAT_ERR"]**2 )
            hd.close()
        s = np.array(signal)
        b = np.array(background)
        S = (np.sum(s) - np.sum(b))*t
        B = np.sum(b)*t
        if showBkg: print "%s Background rate: %.2f +- %.2f" % (index, np.sum(b), np.sqrt( np.sum(bvar) ))
        X = ( S )/( S+B )
        SNR.append( S / np.sqrt(B) )
        BackRate.append( [np.sum(b), np.sqrt(np.sum(bvar))] )
    print "%s SNR: %s" % (index, SNR)
    return (SNR, BackRate)

def printProgress (iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 70):
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

