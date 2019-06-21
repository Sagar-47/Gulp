import pandas as pd
import array as arr
import numpy as np
import matplotlib.pyplot as plt
raw_data=pd.read_table('1_ZnS.txt', sep='\s+',comment='%') 
numpy_array=raw_data.values 
print numpy_array
