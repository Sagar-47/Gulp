import os
import numpy as np
import subprocess
from scipy.optimize import minimize
from scipy.optimize import basinhopping



data4 = np.loadtxt('UD.dat', dtype='str').T #problem with gulp, sometimes it generates 2 neighbouring values not separated by space for bad intial guesses 
print(data4)
    
