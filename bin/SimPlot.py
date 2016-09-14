class SimPlot:
    def __init__(self):
        self.title = ""
        self.subtitles = [""]*6
        self.modelSNR = None
        self.correctedSNR = None
        self.params = []                          # list of parameters(floats)
        self.data = [[], [], [], [], [], []]      # list of spectra   (arrays)
        self.stats = [ConfInt(0,0,0)]*6           # list of statistics(ConfInt)
        self.steps = [Steppars([])]*6             # list of statistics(ConfInt)
        self.xs = [[], [], [], [], [], []]        # list of x-plots   (arrays)
        self.ys = [[], [], [], [], [], []]        # list of y-plots   (arrays)
        self.flags = [""]*6                       # list of flags     (strings)
        self.pnames = [r"$\tau$", r"$\epsilon_e$", r"$\epsilon_B$", r"$L$", r"$\Gamma$", "PG stat"]
        self.gridpoints = [[1.0, 5, 10, 35]
              ,[0.01, 0.05, 0.1, 0.2, 0.3, 0.5]
              ,[1e-6, 0.01, 0.05]
              ,[1.0, 10, 100]
              ,[100.0, 150, 200, 250, 300, 400]
              ,[]]
    
    def AddSubtitle(self, index, str):
        self.subtitles[index] = str
    def AddDataParam(self, dataobj):
        self.data = [[], [], [], [], [], []]
        self.params = [0,0,0,0,0,0]
        for d in dataobj:
            for i in range(0,5):
                self.data[i].append( d.FitParameters[i] )
                self.params[i] = d.ModelParameters[i]
            self.data[5].append( d.PGstat )
            self.params[5] = -1
    def AddStats(self, index, est, lower, upper):
        self.stats[index] = ConfInt(est,lower,upper)
    def AddFlag(self, index, flag):
        self.flags[index] = flag
    def AddSteppars(self,steps):
        for (i,x) in enumerate(steps):
            if x != []:
                temp = Steppars(x[0])
                self.steps[i] = temp
                for s in x[1]: temp.AddSpectra(s)

class Steppars:
    def __init__(self, xs):
        self.spectras = []
        self.xs = xs
    def AddSpectra(self, spec):
        self.spectras.append(spec)
    def IsEmpty(self):
        return len(self.spectras) == 0

class ConfInt:
    def __init__(self, est, lower, upper):
        self.est = est
        self.lower = lower
        self.upper = upper
