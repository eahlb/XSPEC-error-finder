from bin.ErrorFinder import *

# run the error finder for input.txt
# the result is stored in output.txt
e = ErrorFinder("input.txt", "output.txt", base_folder="BASE_FOLDER/")
# the optinal base_folder argument can be set permanently in bin/ErrorFinder.py
# the number of sims to be used is modified in _point.update_errors (line 179 in bin/ErrorFinder.py)