# config.py
import platform
import os

def get_config():
   
    scoreLim = 10   # lowest number of nucleotieds to be considered a match 
    hpLEN = 7  # length of hairpin to be created
   # Define constants
    lowest_proabability = 0.4  # Parameter to control the percentage of lowest-probability hairpin-complement combinations combining with MMP9
    doPLOT=0 #1 for doing the very time consuming subplots of all diagonals with probabillity
    numTOP = 200 #this is the number of hairpin candidates we select for further analysis, started with 100. IF NaN then no restriction
    os_name = platform.system()
    
    # Directory setup
    if os_name == "Linux":
        dataDIR = '/home/michaelbh23/nupackCode/'
    else:
        dataDIR = 'E:/GENEdataNUPACK/JaumeCODEandDATA/'
        
    results_dir = os.path.join(dataDIR, 'results/')
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
        
    return scoreLim, hpLEN, dataDIR, results_dir, doPLOT, numTOP, lowest_proabability