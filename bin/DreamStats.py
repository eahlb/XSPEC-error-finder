from SimOutput import *
from bin.SimPlot import *
from bin.ValidateModel import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
import os.path
import json
import re
import copy
rc('font',size=10)
rc('font',family='serif')
rc('axes',labelsize=10)
notice = [0,1,2,3]
gridpoints = [[1.0, 5, 10, 35]
              ,[0.01, 0.05, 0.1, 0.2, 0.3, 0.5]
              ,[1e-6, 0.01, 0.05]
              ,[1.0, 10, 100]
              ,[100.0, 150, 200, 250, 300, 400]
              ,[]]

def GetTable(prefix, SNR, texname, ignore=[]):
    if len(SNR)==0: return
    print texname
    tex = r'''\begingroup
\def\arraystretch{1.1}%  1 is the default, change whatever you need
\begin{tabular}{cc|cc|ccc|c|c}
\hline \hline
& & \multicolumn{2}{c|}{Grid-points} & \multicolumn{3}{c|}{Triggered flags} & Deviation from & Confidence \\
SNR & & Used & On edge & No. 1 & No. 2 & No. 3 & true value & interval width \\
'''
    pc = [r"$\tau$", "$\epsilon_e$", "$\epsilon_B$", "L", "$\Gamma$"]
    com = ""
    for snr in SNR:
        fname = "%sSNR%s_"%(prefix,snr)
        if not os.path.isfile( "saves/%s1/output.do"%fname ): continue
        (modelSNR, corrSNR) = GetSNR(fname)
        tex += "\hline \n"
        sc = [""]*len(pc)
        used = [""]*len(pc)
        edge = [0]*len(pc)
        sc[0] = "\multirow{%.0f}{*}{%s}" % ( len(pc)-len(ignore), str(snr) )
        # flags
        fcounts = [[0]*len(pc), [0]*len(pc), [0]*len(pc)]
        ofvals = [list([]) for _ in xrange(len(pc))]
        confw = [list([]) for _ in xrange(len(pc))]
        nr = 0
        while os.path.isfile( "saves/%s%s/output.do" % (fname,nr+1) ): nr += 1
        for i in range(1, nr+1):
            if corrSNR!=None and corrSNR[i-1]-snr < -1e-7: continue
            tempPlot = SimPlot()
            # for each simulation
            data = _getData("{0}{1}".format(fname, i))
            tempPlot.AddDataParam(data)
            for j in range(len(pc)):
                if IsOnEdge(data, j): edge[j] += 1
                e, l, u = _getStats( tempPlot.data[j] )
                tempPlot.AddStats(j,e,l,u)
                # get flags
                flag = _getFlags(tempPlot.params[j], tempPlot.stats[j].est, tempPlot.stats[j].lower ,tempPlot.stats[j].upper,
                                tempPlot.data[j], tempPlot.gridpoints[j])
                for (k,f) in enumerate(flag):
                    if f == 'T': fcounts[k][j] += 1
                ofvals[j].append( abs(e-tempPlot.params[j]) )
                confw[j].append( u-l )
        com += "-  SNR%s=%s  " % ( snr,len(ofvals[0]) )
        used[0] = "\multirow{%.0f}{*}{%s}" % ( len(pc)-len(ignore), str(len(ofvals[0])) )
        f1c = ["{0:.0f}".format( float(i) ) for i in fcounts[0]]
        f2c = ["{0:.0f}".format( float(i) ) for i in fcounts[1]]
        f3c = ["{0:.0f}".format( float(i) ) for i in fcounts[2]]
        # off-value
        o1c = [GetTabFormat(v,i) for i,v in enumerate(ofvals)]
        # conf. int.
        c1c = [GetTabFormat(v,i) for i,v in enumerate(confw)]
        ep = 0
        for i,e in enumerate(edge):
            if not i in ignore: ep+=e
        if ep==0:
            edge = [""]*len(pc)
            edge[0] = "\multirow{%.0f}{*}{0}" % (len(pc)-len(ignore))
        for j in range(len(pc)):
            if not j in ignore:
                tex += '''%s & %s &%s&%s  &%s&%s&%s  &%s  &%s \\\\
                ''' % (sc[j],pc[j], used[j],edge[j], f1c[j],f2c[j],f3c[j], o1c[j], c1c[j])
    tex += r'''\hline \hline
\end{tabular}
\endgroup'''
    tex += "\n% tot={0} {1}".format(nr,com)
    with open("tables/%s.tex"%texname,'w') as f:
        f.write(tex)
        f.close()
    src = "/Users/eahlb/coding/python_lib/tables/%s.tex"%texname
    dst = "/Users/eahlb/Dropbox/KTH/SH204X/rapport/tables/%s.tex"%texname
    shutil.copy(src,dst)

def IsOnEdge(data, p):
    mp = data[0].ModelParameters
    if abs(mp[p]-gridpoints[p][0])<1e-8: return True
    if abs(mp[p]-gridpoints[p][-1])<1e-8: return True
    return False

def GetTabFormat(v,i):
    med = np.median(v)
    low = med - np.percentile(v,15.9)
    high = np.percentile(v,84.1) - med
    if i in [0,3,4]:
        r = "$ %s^{+%.1f}_{-%.1f} $"%("{:.1f}".format(med) if med>1e-1 else "<0.1", high, low)
    elif i==1:
        r = "$ %s^{+%.3f}_{-%.3f} $"%("{:.3f}".format(med) if med>1e-3 else "<0.001", high, low)
    else:
        r = "$ %s^{+%.4f}_{-%.4f} $"%("{:.4f}".format(med) if med>1e-4 else "<0.0001", high, low)
    return r

def GetGRBtable(grb):
    print grb
    v = ValidateModel()
    path = "data/%s" % grb
    detectors = []
    #get all detectors
    for file in os.listdir(path):
        if file.endswith(".pha"): detectors.append(file[:-4])
    # order the detector
    # read detector angles
    with open("%s/detector_angles.txt"%path, 'r') as f:
        angles = f.read()
        f.close()
    tex = r'''\begin{tabular}{cc|ccc}
\hline \hline
\multicolumn{2}{c|}{Detector} & \multirow{2}{*}{Viewing angle} & \multirow{2}{*}{SNR} & Background rate  \\
Type & No.  &  & & (counts$\,$s$^{-1}$) \\
\hline
'''
    for det in detectors:
        dname = "NaI" if det[0]=='n' else "BGO"
        dno = str(det[1]) if det[1].isdigit() else str(10+"abcdef".index(det[1]))
        view = re.findall("angle for %s = ([0-9.]+)"%det, angles)[0]
        (s, b) = v.GetSNR(["%s/%s.pha"%(path,det)], ["%s/%s.bak"%(path,det)], noStats=False)
        snr = "%.1f"%s[0] ; brate = "%.0f"%b[0][0]; berror = "%.0f"%b[0][1]
        t = "%.2f"%v._getExpTime("%s/%s.pha"%(path,det))
        tex += "%s & %s & $%s^{\circ}$	&%s & $%s \pm %s$   \\\\ \n" % (dname, dno, view, snr, brate, berror)
    # NaI & 7 & $45.1^{\circ}$	&55.1 & $1274.4 \pm 6.5 \,$ph$\,$s$^{-1}$ &1.88 s  \\
    tex += r'''\hline \hline
\end{tabular}'''
    with open("tables/%s.tex"%grb,'w') as f:
        f.write(tex)
        f.close()
    src = "/Users/eahlb/coding/python_lib/tables/%s.tex"%grb
    dst = "/Users/eahlb/Dropbox/KTH/SH204X/rapport/tables/%s.tex"%grb
    shutil.copy(src,dst)

def CompareInterpolation():
    series1 = "rapp_midgrid_LOG_eb1e-6_"
    series2 = "rapp_midgrid_LIN_eb1e-6_"
    params = [r"$\tau$", r"$\epsilon_e$", "", r"L", r"$\Gamma$"]
    #step = [np.linspace(1,35,35), np.linspace(0.00,0.5,11), None ,np.linspace(0,100,21), np.linspace(100,400,13)]
    snr = [150,100,40,20]
    grid = [ [2.5, 7.5, 22.5]
            ,[0.025, 0.15, 0.4]
            ,[1e-6]
            ,[5, 50]
            ,[275, 350]]
            #,[125, 175, 275, 350]]
    # get all snr
    (modelSNR, corrSNR) = GetSNR("%sSNR20_"%series1)
    bestSNR = []
    # find run sim with best snr for all points
    for m in modelSNR:
        best = 0
        for t in snr:
            if t < m*10.0001:
                best = t
                break
        bestSNR.append(best)
    # load all the data
    data1 = []
    data2 = []
    for i,b in enumerate(bestSNR):
        if b==0:
            data1.append(None)
            data2.append(None)
        else:
            data1.append( _getData("%sSNR%s_%s"%(series1,b,i+1)) )
            data2.append( _getData("%sSNR%s_%s"%(series2,b,i+1)) )
            #break
    # make one plot(s) for each parameter
    for i,pname in enumerate(params):
        if pname=="": continue
        #if pname!="L": continue
        nr = len(grid[i])
        f, ax = plt.subplots(1, nr, figsize=( min([(7.0/3.0)*nr,7.0]) , 2.5), sharey='row')
        f.tight_layout()
        ax[0].set_ylabel("Normalized frequency")
        f.subplots_adjust(wspace=0.1)
        # make one subplot for each parameter value
        for j,pval in enumerate(grid[i]):
            # find all point for this value
            v1 = []
            v2 = []
            for k in range(len(data1)):
                if data1[k]==None: continue
                if abs(pval-data1[k][0].ModelParameters[i]) < 1e-8:
                    # correct params
                    # find fit value
                    temp = _getParameters(data1[k])[i]
                    v1.append( _getStats(temp)[0] )
                    temp = _getParameters(data2[k])[i]
                    v2.append( _getStats(temp)[0] )
            ax[j].set_xlabel(params[i])
            # set plot limits
            ax[j].set_ylim(0, 1.1)
            dev = 1.1*max([abs(pval-max(v1+v2)), abs(pval-min(v1+v2))])
            ax[j].set_xlim(max([pval-dev,0]), min([pval+dev,gridpoints[i][-1]]))
            # plot the values
            x,y = _getHist(v2) # lin
            ax[j].plot(x,y, c="blue", ls="-", linewidth=1.5)
            x,y = _getHist(v1) # log
            ax[j].plot(x,y, c="red", ls="-", linewidth=1.5)
            _plotLine(ax[j], pval, c="black", linestyle="--", max=1.1)
        #ax[j].grid()
        plt.legend( ["Linear","Logarithmic"], loc = 'upper center', bbox_to_anchor = (0.0,0.1,1,1), bbox_transform = plt.gcf().transFigure, ncol=nr )
        f.tight_layout()
        #plt.show()
        # save the plot
        file = "intcomp_%s.eps" % params[i].strip('$').strip(r"\\")
        print file
        plt.savefig(file, format='eps', dpi=1000, bbox_inches='tight')
        plt.close()
        # save plot to dropbox
        src = "/Users/eahlb/coding/python_lib/%s" % file
        dst = "/Users/eahlb/Dropbox/KTH/SH204X/rapport/compplots/%s" % file
        shutil.move(src, dst)

def PlotDist(ID, fname, xlim=None):
    # get data
    tempPlot = SimPlot()
    # for each simulation
    data = _getData(ID)
    tempPlot.AddDataParam(data)
    tempPlot.title = "{0}".format(ID)
    for j in range(0,6):
        e, l, u = _getStats( tempPlot.data[j] )
        tempPlot.AddStats(j,e,l,u)
        # get flags
        flag = _getFlags(tempPlot.params[j], tempPlot.stats[j].est, tempPlot.stats[j].lower ,tempPlot.stats[j].upper,
                    tempPlot.data[j], tempPlot.gridpoints[j])
        tempPlot.AddFlag(j,flag)
        st = _getTitle(tempPlot.pnames[j], tempPlot.stats[j].est, tempPlot.stats[j].lower ,tempPlot.stats[j].upper,
                        tempPlot.params[j], tempPlot.flags[j])
        tempPlot.AddSubtitle(j,tempPlot.pnames[j])
    _setPlot(tempPlot,nrPlots=[0,1,3,4,5], xlim=xlim)
    #plt.show()
    # save the plot
    file = "%s.eps" % fname
    print file
    plt.savefig(file, format='eps', dpi=1000, bbox_inches='tight')
    plt.close()
    # save plot to dropbox
    src = "/Users/eahlb/coding/python_lib/%s" % file
    dst = "/Users/eahlb/Dropbox/KTH/SH204X/rapport/randplots/%s" % file
    shutil.move(src, dst)

def GetOffValuePlot(series, fname, conf=False, useOnlyBest=True, legend=None):
    SNR = [20,40,100,150]
    plots = [OffPlot(), OffPlot(), OffPlot(), OffPlot()]
    leg = [[] for _ in range(len(plots))]
    if legend!=None: leg[0] = legend
    for p in plots: p.snr = SNR
    # find value and error for each series
    for serie in series:
        for p in plots:
            p.res.append([])
            p.uerr.append([])
            p.lerr.append([])
        # find value for each snr
        for snr in SNR:
            prefix = "%sSNR%s_" % (serie, snr)
            # get the snr for the serie
            (modelSNR, corrSNR) = GetSNR(prefix)
            # load data with correct snr
            i = 1
            data = []
            while os.path.isfile( "saves/%s%s/output.do" % (prefix,i) ):
                # check snr
                if useOnlyBest: cond = modelSNR[i-1]-15.0 > -1e-5
                else: cond = True
                if cond and corrSNR[i-1]-snr > -1e-7:
                    data.append( _getData("%s%s"%(prefix,i)) )
                #if i==30: break
                i += 1
            # calculate res and errors
            values = [[],[],[],[]]
            for d in data:
                #raw_input(len(data))
                temp = _getParameters(d)
                # calculate for each parameter
                i = 0
                j = 0
                for p in temp[0:2]+temp[3:5]:
                    fp = _getStats(p)
                    mp = d[0].ModelParameters[j]
                    if not conf: v = abs(fp[0]-mp)   # off value
                    else: v = abs(fp[2]-fp[1])       # conf width
                    values[i].append(v)
                    i += 1
                    if j==1: j+=2
                    else: j+=1
            # add values to plot
            i = 0
            for v in values:
                r = np.percentile(v, 50)
                plots[i].res[-1].append( r )
                plots[i].lerr[-1].append( r-np.percentile(v,15.9) )
                plots[i].uerr[-1].append( np.percentile(v,84.1)-r )
                i += 1
    # plot data
    fig = plt.figure(figsize=(6.5,4))
    params = [r"$\tau$", r"$\epsilon_e$", r"L", r"$\Gamma$"]
    ylim = [[1,35],[0.01,0.5],[1,100],[100,400]]
    for i,p in enumerate(plots):
        ax = fig.add_subplot(2, 2, i+1)
        if i<2: low = False
        else:   low = True
        p.plot(ax, params[i], ylim[i], low, leg[i])
    plt.legend( leg[0], loc = 'upper center', bbox_to_anchor = (0.0,0.1,1,1), bbox_transform = plt.gcf().transFigure, ncol=len(series) )
    #plt.show()
    fig.tight_layout()
    # save the plot
    if conf: file = "confwidth_%s.eps" %fname
    else:    file = "offvalue_%s.eps" %fname
    print file
    plt.savefig(file, format='eps', dpi=1000, bbox_inches='tight')
    plt.close()
    # save plot to dropbox
    src = "/Users/eahlb/coding/python_lib/%s" % file
    dst = "/Users/eahlb/Dropbox/KTH/SH204X/rapport/plots/%s" % file
    shutil.move(src, dst)

def PlotValue(series, pi, pv, fname="plot_value", leg=None):
    # get all data with correct eb and L
    data = []
    if isinstance(series, list): lstPrefix = series
    else: lstPrefix = [series]
    for prefix in lstPrefix:
        data.append([])
        i = 1
        while os.path.isfile( "saves/%s%s/output.do" % (prefix,i) ):
            temp = _getData( "%s%s"%(prefix,i) )
            i += 1
            if abs(temp[0].ModelParameters[pi]-pv)<1e-10:
                data[-1].append( temp[0].FitParameters[pi] )
    f = plt.figure(figsize=(5, 3))
    plt.ylabel(r'$d(x)$')
    plt.xlabel(r'$x$')
    xvals = np.linspace(gridpoints[pi][0],gridpoints[pi][-1])
    if leg==None: leg=[""]*len(series)
    for i,v in enumerate(data):
        diff = []
        for x in xvals:
            d = 0
            for val in v: d += abs(x-val)/float(len(v))
            diff.append(d)
        plt.plot(xvals,diff,label=leg[i])
        _plotLine(plt, pv, c="black", max=max(diff))
    plt.axis('equal')
    plt.xlim([min(xvals),max(xvals)])
    #plt.ylim([min(xvals)*(1-p),max(xvals)*(1+p)])
    if leg[0]!='': plt.legend()
    plt.grid()
    plt.savefig("%s.eps"%fname, format='eps', dpi=1000, bbox_inches='tight')
    plt.show()

def PlotBox(series, eb, L, header=None, snr=0):
    # get all data with correct eb and L
    i = 1
    loads = []
    if isinstance(series, list): lstPrefix = series
    else: lstPrefix = [series]
    if not isinstance(L, list): L = [L]
    if not isinstance(eb, list): eb = [eb]
    for prefix in lstPrefix:
        (modelSNR, corrSNR) = GetSNR(prefix)
        while os.path.isfile( "saves/%s%s/output.do" % (prefix,i) ):
            temp = _getData( "%s%s"%(prefix,i) )
            # check snr
            if corrSNR!=None and not (corrSNR[i-1]-snr < -1e-7):
                loads.append(temp)
            i += 1
    # data now contains all relevant data
    # get the flags for all points
    for l in L:
        for e in eb:
            data = []
            for temp in loads:
                if abs(temp[0].ModelParameters[2]-e)<1e-10 and abs(temp[0].ModelParameters[3]-l)<1e-10:
                    data.append( temp )
            print "---  Getting flags       ---"
            flags = _getAllFlags(data)
            fig = plt.figure(figsize=(6.25,3.5))
            for (i,p) in enumerate([0,1,3,4]):
                ax = fig.add_subplot(2, 2, i, projection='3d')
                _setPlotBox(flags, data, e, l, ax, pFlag=p)
            if header==None: header=r"$\epsilon_b={0},$ $L={1}$  (order $\tau,\epsilon_e,L,\Gamma$)".format(e,l)
            fig.suptitle( header )
            fig.subplots_adjust(top=1, wspace=0.25)
            fig.tight_layout()
            file = "box_%s_lum%s_eb%s.eps" % (series, _getExpForm(l), _getExpForm(e))
            print file
            plt.savefig(file, format='eps', dpi=1000, bbox_inches='tight')
            plt.close()
            # save plot to dropbox
            src = "/Users/eahlb/coding/python_lib/%s" % file
            dst = "/Users/eahlb/Dropbox/KTH/SH204X/rapport/boxplots/%s" % file
            shutil.move(src, dst)

def _getExpForm(f):
    p = np.floor( np.log10(f) )
    r = f/(10**p)
    return "%.0fe%.0f" % (r,p)

def _setPlotBox(flags, data, eb, L, ax, pFlag=0):
    # get color of points
    g = [[],[],[]]
    y = [[],[],[]]
    r = [[],[],[]]
    w = [[],[],[]]
    na= [[],[],[]]
    for (i,f) in enumerate(flags):
        p = data[i][0].ModelParameters
        co = [p[4], p[1], p[0]]
        if _isRed(f,pFlag)     : t = r
        elif _isYellow(f,pFlag): t = y
        elif _isWhite(f,pFlag) : t = w
        elif _isGreen(f,pFlag) : t = g
        else                   : t = na
        for j in range(3):
            t[j].append( co[j] )
    ang = [25,60]
    # plot points
    ax.scatter(g[0],g[1],g[2],c='g',marker='o')
    ax.scatter(y[0],y[1],y[2],c='y',marker='o')
    ax.scatter(r[0],r[1],r[2],c='r',marker='o')
    ax.scatter(w[0],w[1],w[2],c='w',marker='o')
    ax.set_ylabel("\n     $\epsilon_e$     ", linespacing=2.4)
    ax.set_ylim(0,0.5)
    ax.set_xlabel("\n     $\Gamma$     ", linespacing=2.4)
    ax.set_xlim(400,100)
    ax.set_zlabel(r"     $\tau$     ", rotation=30)
    ax.set_zlim(0,35)
    ax.view_init(elev=ang[0], azim=ang[1])
    titles = [r"$(\tau)$",r"$(\epsilon_e)$",r"$(\epsilon_b)$",r"$($L$)$",r"$(\Gamma)$"]
    ax.set_title(titles[pFlag])
def _isRed(f,p):
    if f[p][0] == 'T':
        return True
    return False
def _isYellow(f,p):
    if f[p][1] == 'T':
        return True
    return False
def _isWhite(f,p):
    if f[p][2] == 'T':
        return True
    return False
def _isGreen(f,p):
    return True

def GetSnrEvolution(s=0):
    snr = [20,40,100,150]
    data = []
    if s==0:
        prefix = "rapp_eb5e-2_"
        index = 413
        param = 4
        ticks = [200,300,400]
        xlim = [180,420]
        lbl = [r"$\Gamma$"]*4
    elif s==1:
        prefix = "rapp_5par_"
        index = 36
        param = 2
        ticks = [0,0.05]
        xlim = [-0.01,0.06]
        lbl = [r"$\epsilon_b$"]*4
    elif s==2:
        prefix = "rapp_rand_eb1e-6_"
        index = 125
        param = 4
        ticks = [100,200,300,400]
        xlim = [80,420]
        lbl = [r"$\Gamma$"]*4
    elif s==3:
        prefix = "rapp_eb5e-2_"
        index = 409
        param = 0
        ticks = [20,25,30,35]
        xlim = [17,37]
        lbl = [r"$\tau$"]*4
    elif s==4:
        prefix = "rapp_eb5e-2_"
        index = 409
        param = 3
        ticks = [40,60,80,100]
        xlim = [15,105]
        lbl = [r"L"]*4
    elif s==5:
        prefix = "rapp_5par_"
        index = 197
        param = 1
        ticks = [0.0,0.1,0.2,0.3,0.5]
        xlim = [-0.05,0.55]
        lbl = [r"$\epsilon_e$"]*4
    elif s==6:
        prefix = "rapp_5par_"
        index = 575
        param = 2
        ticks = [0.01,0.05]
        xlim = [-0.01,0.06]
        lbl = [r"$\epsilon_b$"]*4
    elif s==7:
        prefix = "rapp_5par_"
        index = 1289
        param = 2
        ticks = [0.01,0.05]
        xlim = [-0.01,0.06]
        lbl = [r"$\epsilon_b$"]*4
    for t in snr:
        data.append(_getData( "%sSNR%s_%s"%(prefix,t,index) ))
    plotObjects = []
    for d in data:
        tempPlot = SimPlot()
        # for each simulation
        tempPlot.AddDataParam(d)
        for j in range(0,6):
            e, l, u = _getStats( tempPlot.data[j] )
            tempPlot.AddStats(j,e,l,u)
            # get flags
            flag = _getFlags(tempPlot.params[j], tempPlot.stats[j].est, tempPlot.stats[j].lower ,tempPlot.stats[j].upper,
                             tempPlot.data[j], tempPlot.gridpoints[j])
            tempPlot.AddFlag(j,flag)
            st = _getTitle(tempPlot.pnames[j], tempPlot.stats[j].est, tempPlot.stats[j].lower ,tempPlot.stats[j].upper, tempPlot.params[j], tempPlot.flags[j])
            tempPlot.AddSubtitle(j,st)
        plotObjects.append(tempPlot)
    # create subplot
    f, ax = plt.subplots(1, 4, figsize=(7, 2.5), sharey='row')
    f.tight_layout()
    ax[0].set_ylabel("Normalized frequency")
    f.subplots_adjust(wspace=0.1)
    # create plots
    #lbl = ["(a)","(b)","(c)","(d)"]
    for i,obj in enumerate(plotObjects):
        lw = 1.0
        # plot data
        x, y = _getHist(obj.data[param], 30)
        ax[i].plot(x,y,lw=lw)
        # plot conf
        _plotLine(ax[i], obj.stats[param].est, c="blue", linestyle="-", linewidth=lw)
        _plotLine(ax[i], obj.stats[param].lower, c="blue", linewidth=lw)
        _plotLine(ax[i], obj.stats[param].upper, c="blue", linewidth=lw)
        ax[i].set_xlabel( lbl[i] )
        ax[i].set_ylim( [0, 1.1] )
        ax[i].set_xlim(xlim)
        #ax[i].grid()
        ax[i].xaxis.set_ticks(ticks)
    #plt.show()
    # save the plot
    file = "snr_series%s.eps" %s
    print file
    plt.savefig(file, format='eps', dpi=1000, bbox_inches='tight')
    plt.close()
    # save plot to dropbox
    src = "/Users/eahlb/coding/python_lib/%s" % file
    dst = "/Users/eahlb/Dropbox/KTH/SH204X/rapport/plots/%s" % file
    shutil.move(src, dst)

def ProcessSeries(prefix, nrOfSims, savePlot=False, plotSteppar=False, dsPrefix=None, stepAll=False, useFak=False, useOldModel=False, nrBins=30, snr=None):
    # get data
    plotObjects = []
    (modelSNR, corrSNR) = GetSNR(prefix)
    for i in range(1, nrOfSims+1):
        if snr!=None and abs(corrSNR[i-1]-snr)>1e-7: continue
        tempPlot = SimPlot()
        # for each simulation
        data = _getData("{0}{1}".format(prefix, i))
        tempPlot.AddDataParam(data)
        tempPlot.title = "{0}{1}".format(prefix, i)
        if modelSNR != None and corrSNR != None:
            tempPlot.modelSNR = modelSNR[i-1]
            tempPlot.correctedSNR = corrSNR[i-1]
            tempPlot.title = r"%s (SNR: %.2f$\rightarrow$%.2f)" % (tempPlot.title, modelSNR[i-1], corrSNR[i-1])
        for j in range(0,6):
            e, l, u = _getStats( tempPlot.data[j] )
            tempPlot.AddStats(j,e,l,u)
            # get flags
            flag = _getFlags(tempPlot.params[j], tempPlot.stats[j].est, tempPlot.stats[j].lower ,tempPlot.stats[j].upper,
                             tempPlot.data[j], tempPlot.gridpoints[j])
            tempPlot.AddFlag(j,flag)
            st = _getTitle(tempPlot.pnames[j], tempPlot.stats[j].est, tempPlot.stats[j].lower ,tempPlot.stats[j].upper,
                             tempPlot.params[j], tempPlot.flags[j])
            tempPlot.AddSubtitle(j,st)
        plotObjects.append(tempPlot)

    # save series
    with open("saves/{}".format(prefix), 'w') as f:
        for i,p in enumerate(plotObjects):
            tempFlags = ["", "", "", "", ""]
            found = False
            for j in range(0,5):
                tempFlags[j] = p.flags[j]
                if 'T' in tempFlags[j]:
                    found = True
            if found:
                s = "%s # %s # %s \n" % ( json.dumps(i+1), json.dumps(p.params), json.dumps(tempFlags) )
                f.write(s)
    f.close()
    # save plots
    ## -- must map new indices to old if we want to plot steppar --
    if plotSteppar:
        stepp = _getStepparMatrix(prefix, dsPrefix, nrOfSims, stepAll, useFak=useFak, useOldModel=useOldModel)
        for (i,s) in enumerate(stepp): plotObjects[i].AddSteppars(s)
    if savePlot:
        tag = ""
        if stepAll:
            tag = "ALL"
        file = "plots/series_{0}{1}.pdf".format(prefix, tag)
        with PdfPages(file) as pdf:
            for i,p in enumerate(plotObjects):
                print( "Saving page {0}".format(i+1) )
                # save all plots in a pdf
                _setPlot(p)
                pdf.savefig()
                plt.close()
        # save plot to dropbox
        src = "/Users/eahlb/coding/python_lib/%s" % file
        dst = "/Users/eahlb/Dropbox/KTH/SH204X/python_lib/%s" % file
        shutil.copy(src, dst)

def SaveSteppar(prefix, dsPrefix, stepAll=False, nrOfSims=-1, useFak=False, fakIndex=None, firstSim=0, useOldModel=False):
    keeps = [] # [sim nr, param]
    if stepAll:
        if nrOfSims == -1:
            raise IOError("The number of sims must be specified")
        for sim in range(firstSim+1,nrOfSims+1):
            for i in range(0, 5):
                keeps.append( [sim, i] )
    else:
        # run steppar for possibly multipeaked spectra
        file = "saves/{}VERIFIED".format(prefix)
        if not os.path.isfile(file):
            print( "Verified flags for series %s does not exist. Using raw flags" % prefix )
            file = "saves/{}".format(prefix)
            if not os.path.isfile(file):
                print( "Flags for series %s does not exist. Run ProcessSeries" % prefix )
                return
        # find all [xTxx] flags
        with open(file) as sr:
            s = sr.readline().rstrip()
            while True:
                if s == "":
                    break
                t = s.split("#", 3)
                sim = int( t[0] )
                flag = json.loads( t[2].strip() )
                for i in range(0, len(flag)):
                    if sim > firstSim and sim < nrOfSims+1 and flag[i][1] == "T":
                        keeps.append( [sim, i] )
                s = sr.readline().rstrip()
        sr.close()
    # run steppar
    v = ValidateModel()
    res = []
    xs  = []
    for tuple in keeps:
        pb = [ gridpoints[tuple[1]][0], gridpoints[tuple[1]][-1] ]
        print "----     Stepping sim:%s     ----" % tuple[0]
        if useFak or fakIndex==None:
            fakIndex = _getBestFitIndex( "%s%s" % (prefix,tuple[0]), 5 )
        else: fakIndex=[None]
        s = []
        for f in fakIndex:
            s.append( v.GetSteppar(prefix, dsPrefix, tuple[0], tuple[1], pb, logStep=False, nrSteps=40,
                                   useFak=useFak, fakIndex=f, useOldModel=useOldModel) )
        st = _getParamNames()
        temp = []
        tempBounds = []
        for i in range(0,5):
            tempBounds.append( pb )
            if i == tuple[1]:
                temp.append( s )
            else:
                temp.append( np.zeros(len(s)) )
        x = np.linspace( pb[0], pb[1], len(s[0]) )
        res.append( s )
        xs.append( x.tolist() )
    # save steppars
    file = "saves/{}STEPPAR".format(prefix)
    if useFak:
        file += "_FAK_"
    if stepAll:
        file += "ALL"
    with open(file, 'w') as sw:
        sw.write( json.dumps(keeps) )
        sw.write( "#" )
        sw.write( json.dumps(res) )
        sw.write( "#" )
        sw.write( json.dumps(xs) )
    sw.close()

def GetSNR(ID):
    if os.path.isfile("saves/%s_SNR"%ID):
        with open("saves/%s_SNR"%ID, 'r') as f:
            temp = f.read().split('#')
            f.close()
        m = json.loads(temp[0])
        c = json.loads(temp[1])
        return (m,c)
    else:
        return (None, None)

def CheckParameterPeaks(ID, parameter, points, nrBins=30, plotSteppar=False, dsPrefix=None, useFak=False):
    groups, bestFit = _getGroups(ID, parameter, points)
    # group[i] now holds data object for that point
    file = "plots/peaks_%s_%s.pdf" % (ID, parameter)
    prefix = re.split('[0-9]+',ID, maxsplit=1)[0]
    simNr  = int( re.split('[a-z]+',ID, maxsplit=1)[1] )
    n = _getPlotNorm( _getParameters( _getData(ID) ) )
    with PdfPages(file) as pdf:
        for (i,g) in enumerate(groups):
            if plotSteppar:
                stepp = _getStepparMatrix( prefix, dsPrefix, simNr, stepAll=True, useFak=useFak, fakIndex=bestFit[i], firstSim=simNr-1 )
                _setParamPlot(ID, nrBins=nrBins, params=g, norm=n, stepp=stepp[0])
            else:
                _setParamPlot(ID, nrBins=nrBins, params=g, norm=n)
            pdf.savefig()
            plt.close()

def _getBestFitIndex(ID, nr=1):
    indices = []
    data = _getData(ID)
    modelParams = data[0].ModelParameters
    params = _getParameters(data)
    res = []
    for (i,x) in enumerate(params[0]):
        t = 0
        for j in range(5):
            t += abs( params[j][i] - modelParams[j] )/abs( gridpoints[j][0]-gridpoints[j][-1] )
        res.append(t)
    for i in range(nr):
        indices.append( res.index(min(res))+1 )
        res[res.index(min(res))] = max(res)+1
    return indices

def _getGroups(ID, parameter, points):
    # load the data
    data = _getParameters( _getData(ID) )
    # check which indices belong to which peak
    groups = [ [[] for i in range(len(data))] for j in range(len(points)+1) ]
    groups[0] = data
    for (i,d) in enumerate(data[parameter]):
        for (j,p) in enumerate(points):
            tol = max( [1e-3, 0.1*gridpoints[parameter][p]] )
            if abs( gridpoints[parameter][p]-d ) < tol:
                for k in range( 0, len(groups[j]) ):
                    groups[j+1][k].append( data[k][i] )
    # group[i] now holds data object for that point
    indices = []
    for d in groups:
        indices.append( data[-1].index(min(d[-1]))+1 )
    return groups, indices

def _getStepparMatrix(prefix, dsPrefix, nrOfSims, stepAll=False, useFak=False, fakIndex=None, firstSim=0, useOldModel=False):
    # check if saved steppars exist
    file = "saves/{}STEPPAR".format(prefix)
    if useFak:
        file += "_FAK_"
    if firstSim != 0:
        file += str( nrOfSims )
    if stepAll:
        if not os.path.isfile( "%sALL" % file ):
            print( "No saved steppars exist, saving steppars" )
            if dsPrefix != None:
                SaveSteppar( prefix, dsPrefix, stepAll, nrOfSims, useFak=useFak, fakIndex=fakIndex, firstSim=firstSim, useOldModel=useOldModel )
            else:
                raise IOError("No ds file is given and no saved (all) steppars exist")
    else:
        if not os.path.isfile( file ):
            print( "No saved steppars exist, saving steppars" )
            if dsPrefix != None:
                SaveSteppar( prefix, dsPrefix, stepAll, nrOfSims, useFak=useFak, fakIndex=fakIndex, firstSim=firstSim, useOldModel=useOldModel )
            else:
                raise IOError("No ds file is given and no saved steppars exist")
    # load steppars
    file = "saves/{}STEPPAR".format(prefix)
    if useFak:  file += "_FAK_"
    if stepAll: file += "ALL"
    with open(file) as sr:
        val = sr.readline().split("#")
        keeps = json.loads( val[0] )
        res = json.loads( val[1] )
        xs = json.loads( val[2] )
    sr.close()
    # build steppar data
    stepp = []
    for i in range(firstSim,nrOfSims):
        # check if steppar exist for this sim
        ind = []
        temp = []
        for j in range(0,len(keeps)):
            if keeps[j][0] == i+1:
                ind.append(j)
        # build steppar list
        if len(ind) != 0:
            temp = [ [],[],[],[],[] ]
            for j in ind:
                # keeps[j][1] parameter index
                temp[ keeps[j][1] ].append( xs[j] )
                temp[ keeps[j][1] ].append( res[j] )
        stepp.append(temp)
    return stepp

def _setParamPlot(ID="last", nrBins=10, title=None, data=None, stepp=[], stepOnly=False, mu=None, lowConf=None, upConf=None, params=None, norm=None):
    if title == None:
        title = "%s (%s)" % (ID,_getBestFitIndex(ID))
    if data == None:
        data = _getData(ID)
    if mu==None or lowConf==None or upConf==None:
        mu = []
        lowConf = []
        upConf = []
    if params == None:
        params = _getParameters(data)
    pNames = _getParamNames()
    pBounds= _getParamBounds()
    subtitles = []
    modelParams = [-1, -1, -1, -1, -1, -1]
    # get all parameters
    modelParams[:5] = data[0].ModelParameters[:5]
    # set PG stat bounds
    pBounds[5][0] = 0#.9*min( params[5] )
    pBounds[5][1] = 600#1.1*max( params[5] )
    # get subtitles
    for i in range(0,6):
        if len(mu) < 6:
            x, y, z = _getStats(params[i])
            mu.append(x)
            lowConf.append(y)
            upConf.append(z)
        flags = _getFlags( modelParams[i], mu[i] ,lowConf[i], upConf[i], params[i], gridpoints[i] )
        subtitles.append( _getTitle(pNames[i], mu[i], lowConf[i], upConf[i], modelParams[i], flags) )
    _setPlot(params, title, subtitles ,pBounds, nrBins, showLines=True, stepp=stepp, stepOnly=stepOnly, mu=mu, lowConf=lowConf, upConf=upConf, norm=norm)

def _displayFlags(flag, p):
    msg = ["Actual value not within confidence interval",
           "Confidence interval greater than 30% of estimate",
           "Parameter value possibly unbound",
           ""]
    for i in range(0,len(flag)):
        if flag[i] == 'T' and i in notice:
            nm = _getParamNames()[p]
            pad = ""
            for j in range(0, 12-len(nm)):
                pad += " "
            print "%s %s : %s" % (nm, pad , msg[i])

def _setPlot(obj, nrBins=30, nrPlots=None, showLines=True, norm=None, xlim=None):
    if nrPlots==None:
        nrPlots=range(6)
        fs = (23, 6)
    else: fs = (7.5,2.5)
    ticks = [[1,5,10,35],
            [0.0,0.1,0.2,0.3,0.5],
            [1e-6,0.01,0.05],
            [1,10,100],
            [100,200,300,400],
            [300,400,500,600,700]]
    # create subplot
    f, ax = plt.subplots(1, len(nrPlots), figsize=fs, sharey='row')
    f.tight_layout()
    ax[0].set_ylabel("Normalized frequency")
    if nrPlots==range(6): f.suptitle(obj.title)
    f.subplots_adjust(wspace=0.1)
    # create plots
    for i,j in enumerate(nrPlots):
        lw = 1.5
        # plot data
        if norm == None: x, y = _getHist(obj.data[j], nrBins)
        else: x, y = _getHist(obj.data[j], nrBins, norm=norm[j])
        ax[i].plot(x,y,lw=1.0)
        # plot steppar
        if not obj.steps[i].IsEmpty():
            tax = ax[i].twinx()
            nr = len(obj.steps[j].spectras)
            c = 0.0
            for tl in tax.get_yticklabels():
                tl.set_color('r')
            for (k,t) in enumerate(obj.steps[j].spectras):
                c = c + 0.6/float(nr)
                tax.plot(obj.steps[j].xs, obj.steps[j].spectras[k], color=(1.0,c,c), lw=lw)
            tax.yaxis.set_ticks(np.arange(-1, 1.75, 0.5))
            tax.set_ylim([ -1, 1.75])
        # plot conf
        _plotLine(ax[i], obj.stats[j].est, c="blue", linestyle="-", linewidth=lw)
        _plotLine(ax[i], obj.stats[j].lower, c="blue", linewidth=lw)
        _plotLine(ax[i], obj.stats[j].upper, c="blue", linewidth=lw)
        #ax[i].set_title( obj.subtitles[j] )
        ax[i].set_xlabel( obj.subtitles[j] )
        ax[i].set_ylim( [0, 1.1] )
        if obj.gridpoints[j] == []: pass #ax[i].set_xlim( 0, 700 )
        else:
            _plotLine(ax[i], obj.params[j], c="red", linestyle="-.", linewidth=lw)
            ax[i].set_xlim( 0.5*obj.gridpoints[j][0], 1.1*obj.gridpoints[j][-1] )
        ax[i].grid()
        ax[i].xaxis.major.locator.set_params(nbins=4)
        if xlim!=None and xlim[i]!=[]: ax[i].set_xlim( xlim[i] )
        i += 1

def _isMultipeaked(data, points):
    return _getPeaks(data, points) != 1

def _getPeaks(data, points, lowConf, upConf):
    # compare counts within 5% of gridpoints to 10%
    tol = 0.05
    abstol = 1e-3
    bin = np.zeros( len(points) )
    back = np.zeros( len(points) )
    for d in data:
        for i in range(0, len(points)):
            intol = abs(d-points[i]) < tol*points[i] or abs(d-points[i]) < abstol
            inconf = d < lowConf or d > upConf
            if intol and inconf:
                bin[i] += 1
    res = []
    for (i,p) in enumerate(points):
        if p < lowConf or p > upConf:
            res.append( bin[i] )
    peaks = 0
    # noramlize bin counts
    bin = bin / len(data)
    for b in bin:
        # in the counts in the bin is 5% of the maximum bin
        if b > 0.05:
            peaks += 1
    return peaks

def _getAllFlags(data):
    res = []
    for sim in data:
        temp = []
        p = _getParameters( sim )
        for i in range(5):
            value = sim[0].ModelParameters[i]
            m, l, u = _getStats( p[i] )
            f = _getFlags( value, m, l, u, p[i], gridpoints[i] )
            temp.append(f)
        res.append(temp)
    return res

def _getFlags(value, mu, lowConf, upConf, data, gridpoints):
    if gridpoints == []:
        return "xxxx"
    f = ""
    w = gridpoints[-1] - gridpoints[0]
    if 0 in notice:
        inConf = (value >= lowConf) and (value <= upConf)
        acc = abs( mu-value ) < 0.05*w
        if inConf or acc:
            f += "F"
        else:
            f += "T"
    else:
        f += "x"
    if 1 in notice:
        narrow = ( (upConf-lowConf) / w) > 0.4
        if narrow:
            f += "T"
        else:
            f += "F"
    else:
        f += "x"
    if 2 in notice:
        per = 0.01
        low = abs(lowConf-gridpoints[0]) < ( per*w )
        high = abs(upConf-gridpoints[-1]) < ( per*w )
        if low or high:
            f += "T"
        else:
            f+= "F"
    else:
        f += "x"
    if 3 in notice:
        f += str( _getPeaks(data,gridpoints,lowConf,upConf) )
    else:
        f += "x"
    return f

def _getConcentration(data, mu=None, lowConf=None, upConf=None, lvl=1):
    i = 0
    if mu==None or lowConf==None or upConf==None:
        mu, lowConf, upConf = _getStats(data)
    for x in data:
        if (x > lowConf) and (x < upConf):
            i += 1
    return float(i)/len(data)

def _getStats(data):
    # use the mean to estimate the value
    temp = []
    mu = np.median(data)
    lowConf = np.percentile(data, 15.9)
    upConf = np.percentile(data, 84.1)
    return mu, lowConf, upConf

def _getData(ID="last"):
    temp = ValidateModel()
    return temp.LoadData(ID)

def _getBins(data, nrBins):
    if nrBins == -1:
        b = _getBayesianBlocks(data)
        b = np.append(b, np.amax(data))
    else:
        if np.amin(data) != np.amax(data):
            b = np.linspace( np.amin(data), np.amax(data), nrBins+1)
            b = np.append(b, np.amax(b)+b[1]-b[0])
        else:
            b = [np.amin(data), np.amin(data)+1e-5]
    return b

def _getHist(data, nrBins=10, norm=None, normed=True):
    if data == []: return [0],[0]
    b = _getBins(data, nrBins)
    n = np.zeros( len(b)-1 )
    # bin the data
    for d in data:
        i = 0
        while d > b[i]:
            i += 1
        n[i] += 1
    # create x & y
    x = np.zeros( 2*len(n)+2 )
    y = np.zeros( 2*len(n)+2 )
    i = 1
    for d in n:
        y[2*i-1] = d
        x[2*i-1] = b[i-1]
        y[2*i]   = d
        x[2*i]   = b[i]
        i += 1
    # add endings to curves
    x[0] = x[1]
    x[-1] = x[-2]
    x = x-(x[2]-x[0])*0.5
    if norm == None: norm = np.amax(y)
    if normed: y = y/norm
    return x, y

def _getPlotNorm(data):
    n = []
    for d in data:
        x, y = _getHist(d, normed=False)
        n.append( np.amax(y) )
    return n

def _plotLine(loc, x, linewidth=1, linestyle="--", c="red", max=2e4):
    loc.plot([x, x], [0, max], linewidth=linewidth, linestyle=linestyle, c=c)

def _getTitle(name, mu, lowConf, upConf, modelParam=-1, flag=None):
    if modelParam != -1:
        mp = "=%s" % modelParam
    else:
        mp = ""
    t = r"%s$%s \,(\hat{\theta}=%.3g_{%.3g}^{%.3g}) $" % (name, mp, mu, lowConf, upConf)
    if flag != None and modelParam != -1:
        t += "[%s]" % flag
    return t

def _getParamBounds():
    return [[1,35], [0.01,0.5], [1e-6,0.05], [1,100], [100,400], [1,1000]]

def _getParamNames():
    return [r"$\tau$", r"$\epsilon_e$", r"$\epsilon_B$", r"$L$", r"$\Gamma$", "PG stat"]

def _getParameters(data):
    params = [[], [], [], [], [], []]
    for d in data:
        for i in range(0,5):
            params[i].append( d.FitParameters[i] )
        params[5].append( d.PGstat )
    return params

def _getBayesianBlocks(t):
    """Bayesian Blocks Implementation
        
        By Jake Vanderplas.  License: BSD
        Based on algorithm outlined in http://adsabs.harvard.edu/abs/2012arXiv1207.5578S
        
        Parameters
        ----------
        t : ndarray, length N
        data to be histogrammed
        
        Returns
        -------
        bins : ndarray
        array containing the (N+1) bin edges
        
        Notes
        -----
        This is an incomplete implementation: it may fail for some
        datasets.  Alternate fitness functions and prior forms can
        be found in the paper listed above.
        """
    # copy and sort the array
    t = np.sort(t)
    N = t.size
    # create length-(N + 1) array of cell edges
    edges = np.concatenate([t[:1], 0.5 * (t[1:] + t[:-1]), t[-1:]])
    block_length = t[-1] - edges
    # arrays needed for the iteration
    nn_vec = np.ones(N)
    best = np.zeros(N, dtype=float)
    last = np.zeros(N, dtype=int)
    #-----------------------------------------------------------------
    # Start with first data cell; add one cell at each iteration
    #-----------------------------------------------------------------
    for K in range(N):
        # Compute the width and count of the final bin for all possible
        # locations of the K^th changepoint
        width = block_length[:K + 1] - block_length[K + 1]
        count_vec = np.cumsum(nn_vec[:K + 1][::-1])[::-1]
        # evaluate fitness function for these possibilities
        fit_vec = count_vec * (np.log(count_vec) - np.log(width))
        fit_vec -= 4  # 4 comes from the prior on the number of changepoints
        fit_vec[1:] += best[:K]
        # find the max of the fitness: this is the K^th changepoint
        i_max = np.argmax(fit_vec)
        last[K] = i_max
        best[K] = fit_vec[i_max]
    #-----------------------------------------------------------------
    # Recover changepoints by iteratively peeling off the last block
    #-----------------------------------------------------------------
    change_points =  np.zeros(N, dtype=int)
    i_cp = N
    ind = N
    while True:
        i_cp -= 1
        change_points[i_cp] = ind
        if ind == 0:
            break
        ind = last[ind - 1]
    change_points = change_points[i_cp:]
    return edges[change_points]

class OffPlot:
    def __init__(self):
        self.res = []
        self.uerr = []
        self.lerr = []
        self.snr = []
    def plot(self, ax, p, ylim, low, leg):
        ax.set_ylabel(p)
        ax.xaxis.set_ticks(self.snr)
        ax.set_xlim([10,160])
        if low: ax.set_xlabel("SNR")
        ax.grid()
        ls = ['b-o','r--o','g-.o']
        if len(self.res) == 3: offset=[-4,0,4]
        else: offset=[-2,2]
        # print each series
        i = 0
        snr = np.array(self.snr)
        for (r,l,u,) in zip(self.res,self.lerr,self.uerr):
            ax.errorbar(snr+offset[i], r, yerr=[l, u], fmt=ls[i], linewidth=1.5)
            i += 1
        #if leg!=[]: ax.legend(leg, loc='upper center', bbox_to_anchor=(1.0, 1.35), ncol=len(leg))






