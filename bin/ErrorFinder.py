from ValidateModel import *
import re
import os.path
import numpy as np
import shutil
import copy

BASE_FOLDER = "/Users/eahlb/coding/"

TABLE_FOLDER = "/Users/eahlb/Dropbox/KTH/SH204X/rapport/tables/"
# this folder is a remnant from my thesis
# use at own risk

class ErrorFinder:
    def __init__(self, input_path, output_path, base_folder=None):
        if base_folder != None:
            global BASE_FOLDER
            BASE_FOLDER = base_folder
        self.__verify_structure()
        self.__first_line = ""
        self.pha_folder = None
        self.pha_names = None
        self.bak_folder = None
        self.bak_names = None
        self.rsp_folder = None
        self.rsp_names = None
        self.ignores = None
        self.model_path = None
        self.__indices = [] # list of points
        self.input_path = input_path
        self.output_path = output_path
        self.__parse_input(input_path)
        self.__initial_indices = copy.deepcopy(self.__indices)
        self.update_errors()
    
    def __parse_input(self, path):
        f = open(path)
        l = f.read().split("\n")
        f.close()
        self.pha_folder = l[0]
        self.pha_names = self.__parse_names(l[1])
        self.bak_folder = l[2]
        self.bak_names = self.__parse_names(l[3])
        self.rsp_folder = l[4]
        self.rsp_names = self.__parse_names(l[5])
        self.ignores = self.__parse_tuple(l[6])
        self.model_path = l[7]
        self.__parse_result(BASE_FOLDER+l[8])
    
    def __parse_names(self, names):
        final = []
        names = names.strip()
        names = names[1:-1]
        for tp in re.findall('\(.*?\)',names):
            final.append( self.__parse_tuple(tp) )
        return final
    
    def __parse_tuple(self, tp):
        temp = []
        for t in tp.strip()[1:-1].split(','):
            temp.append( t.strip() )
        return tuple(temp)
    
    def __parse_result(self, path):
        f = open(path)
        lines = f.read().split("\n")
        f.close()
        self.__first_line = lines[0]
        for l in lines[1:]:
            if l=="": continue
            p = _point(l)
            self.__indices.append(p)
    
    def __write_result(self, path):
        if self.__indices==[]: return
        s = self.__first_line
        for p in self.__indices:
            s += "\n" + str(p)
        # write to file
        f = open(path,'w')
        f.write(s)
        f.close()

    def __get_ds(self,p):
        beg = '''%s
Grid''' % len(self.pha_names[p])
        end = "##"
        for i in range( len(self.pha_names[p]) ):
            end += '''
%s%s
%s%s
%s%s
%s
#''' % (self.pha_folder,self.pha_names[p][i], self.rsp_folder,self.rsp_names[p][i], self.bak_folder,self.bak_names[p][i], self.ignores[i])
        return (beg,end)

    def __verify_structure(self):
        old_path = os.getcwd()
        os.chdir(BASE_FOLDER)
        if not os.path.isdir("saves"): os.makedirs("saves")
        if not os.path.isdir("data"): os.makedirs("data")
        if not os.path.isdir("data/temp"): os.makedirs("data/temp")
        os.chdir(old_path)

    def update_errors(self):
        if self.__indices==[]: return
        for p in self.__indices:
            p.update_errors( self.__get_ds(p.ind), self.model_path )
        self.__write_result(self.output_path)

    def write_comparison_table(self, tblname, indices):
        tex = r'''\begingroup
\def\arraystretch{1.1}%  1 is the default, change whatever you need
\begin{tabular}{cc|ccccc}
\hline \hline
Bin & & Value & Error &  MC  \\
'''
        get_err1 = lambda p, m: r"$^{+%.1f}_{-%.1f}$" % (p,m)
        get_err3 = lambda p, m: r"$^{+%.3f}_{-%.3f}$" % (p,m)
        for i in indices:
            erp = self.__indices[i]
            mcp = self.__initial_indices[i]
            b = "\multirow{4}{*}{%s}" % i
            # endline
            tex += r"\hline"
            tex += "\n"
            # tau
            tex += r"{:s} & $\tau$ & {:.1f} & {:s} & {:s} \\".format(b, erp.Tau[0], get_err1(erp.Tau[1],erp.Tau[2]), get_err1(mcp.Tau[1],mcp.Tau[2]))
            tex += "\n"
            # e_e
            tex += r" & $\epsilon_e$ & {:.3f} & {:s} & {:s} \\".format(erp.epsilon_e[0], get_err3(erp.epsilon_e[1],erp.epsilon_e[2]), get_err3(mcp.epsilon_e[1],mcp.epsilon_e[2]))
            tex += "\n"
            # L
            tex += r" & L & {:.1f} & {:s} & {:s} \\".format(erp.LRGB[0], get_err1(erp.LRGB[1],erp.LRGB[2]), get_err1(mcp.LRGB[1],mcp.LRGB[2]))
            tex += "\n"
            # gamma
            tex += r" & $\Gamma$ & {:.1f} & {:s} & {:s} \\".format(erp.Gamma[0], get_err1(erp.Gamma[1],erp.Gamma[2]), get_err1(mcp.Gamma[1],mcp.Gamma[2]))
            tex += "\n"
        
        tex += r'''\hline \hline
\end{tabular}
\endgroup'''
        with open("%s.tex"%tblname,'w') as f:
            f.write(tex)
            f.close()
        src = "%s.tex"%tblname
        dst = TABLE_FOLDER + src
        shutil.move(src,dst)

class _point:
    def __init__(self, row):
        self.ind = -1
        self.pgstat = -1
        self.dof = -1
        self.Tau = (-1,-1,-1)
        self.epsilon_e = (-1,-1,-1)
        self.LRGB = (-1,-1,-1)
        self.Gamma = (-1,-1,-1)
        self.z = -1
        self.norm = -1
        self.is_updated = False
        self.parse_row(row)
    
    def update_errors(self, ds, model_path):
        old_path = os.getcwd()
        os.chdir(BASE_FOLDER)
        # create temp sim file
        tempds = "data/temp.ds"
        while os.path.isfile(tempds): tempds = tempds[:-3] + str(np.random.randint(9)) + tempds[-3:]
        f = open(tempds,'w')
        f.write(ds[0]+'\n')
        f.write( str(self.get_parameters()) )
        f.write('\n'+ds[1])
        f.close()
        
        # run simulations
        v = ValidateModel()
        ID = tempds[5:-3]
        v.LoadSimulation(ID, nrOfSims=1e3, model=model_path, dataFolder=".")
        v.RunSimulations(ID, prompt=False)
        
        # get parameters
        out = v.LoadData(ID)
        params = self.extract_parameters(out)
        
        # update errors
        newpar = []
        for (i,par) in zip([0,1,3,4],[self.Tau, self.epsilon_e, self.LRGB, self.Gamma]):
            err = self.calculate_stats(i, params, par[0])
            newpar.append( (par[0], err[1], err[2]) )
        self.Tau = newpar[0]
        self.epsilon_e = newpar[1]
        self.LRGB = newpar[2]
        self.Gamma = newpar[3]
        
        # remove temp sim files and clear data
        v.Clear()
        shutil.rmtree("saves/%s"%ID)
        os.remove(tempds)
        os.chdir(old_path)
    
    def parse_row(self, row):
        s = []
        for e in row.split(): s.append( float(e) )
        self.ind = int(s[0])
        self.pgstat = s[1]
        self.dof = s[2]
        self.Tau = (s[3],s[4],s[5])
        self.epsilon_e = (s[6],s[7],s[8])
        self.LRGB = (s[9],s[10],s[11])
        self.Gamma = (s[12],s[13],s[14])
        self.z = s[15]
        self.norm = s[16]
    
    def get_parameters(self):
        return [self.Tau[0], self.epsilon_e[0], 1e-6, self.LRGB[0], self.Gamma[0], self.z, self.norm]
    
    def __str__(self):
        return "{:<11.0f}{:<11.1f}{:<11.0f}{:<11.2f}{:<11.2f}{:<11.2f}{:<11.2e}{:<11.2e}{:<11.2e}{:<11.2f}{:<11.2f}{:<11.2f}{:<11.2f}{:<11.2f}{:<11.2f}{:<11.2f}{:<11.2f}".format(self.ind,self.pgstat,self.dof, self.Tau[0],self.Tau[1],self.Tau[2], self.epsilon_e[0],self.epsilon_e[1],self.epsilon_e[2], self.LRGB[0],self.LRGB[1],self.LRGB[2], self.Gamma[0],self.Gamma[1],self.Gamma[2], self.z,self.norm)

    def extract_parameters(self, data):
        params = [[], [], [], [], [], []]
        for d in data:
            for i in range(0,5):
                params[i].append( d.FitParameters[i] )
            params[5].append( d.PGstat )
        return params

    def calculate_stats(self, i, params, p):
        data = params[i]
        mu = p
        low = mu-np.percentile(data, 15.9)
        up = np.percentile(data, 84.1)-mu
        if up<0.:
            up = max(data)-mu
            low = mu-np.percentile(data, 31.8)
        elif low<0.:
            low = mu-min(data)
            up = np.percentile(data, 68.2)-mu
        return mu, up, low

## testing
#e = ErrorFinder("input.txt", "output.txt")
#e.write_comparison_table("conf1",range(10))
#e.write_comparison_table("conf2",range(10))










