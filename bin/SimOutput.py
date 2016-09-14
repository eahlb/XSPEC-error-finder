import json
class SimOutput:
## Data holder for simulation output. Also handle formating and save/load of data
## Attributes
#   ModelParameters     : parameters used to simulate the fake spectrum
#   FitParameters       : parameters found by fitting the model to fake spectrum
## Methods
#   __str__()           : output string, used for saving as ASCII
    def __init__(self):
        #### Attributes ####
        self.PGstat = -1.0
        self.Chi2 = -1.0
        self.Dof = -1
        self.ModelParameters = []
        self.FitParameters = []
        self.ParameterNames = []
        #### Variables #####
        #### Constructor ###
#### Methods ###########
    def __str__(self):
        s = json.dumps(self.ModelParameters) + "\n"
        s += json.dumps(self.FitParameters) + "\n"
        s += json.dumps(self.ParameterNames) + "\n"
        s += json.dumps(self.PGstat) + "\n"
        s += json.dumps(self.Chi2) + "\n"
        s += json.dumps(self.Dof) + "\n"
        s += "#"
        return s

#0 = self._model.Photosphere.Tau
#1 = self._model.Photosphere.epsilon_e
#2 = self._model.Photosphere.epsilon_B
#3 = self._model.Photosphere.LGRB
#4 = self._model.Photosphere.Gamma
#5 = self._model.Photosphere.z
#6 = self._model.Photosphere.norm