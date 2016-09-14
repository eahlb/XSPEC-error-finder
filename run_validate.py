# example use of validator
from bin.ValidateModel import *
import settings as s

# create .ds setting files
s.CreateSims("grb100414","smalltest",[1e-6])

# run the simulations
# in this example only one sim per point is done
v = ValidateModel()
v.RunSimSeries("smalltest", "smalltest", nrOfSims=1, silent=True, prompt=False, model='atable{BASE_FOLDER/TableModel.fits}')

# example of how one index may be loaded
# to get all data simply iterate over all indices
res = v.LoadData("smalltest1")
print res