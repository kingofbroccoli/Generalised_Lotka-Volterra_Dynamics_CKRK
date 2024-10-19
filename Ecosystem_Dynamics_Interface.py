import numpy as np
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path

def count_pattern_in_dir(directory, pattern='*'):
    count = 0
    for item in directory.iterdir():
        if (item.is_file() and item.match(pattern)):
            count += 1
    return count

plt.rcParams.update({
    "text.usetex": True,
    "font.family" : "serif",
    "font.serif" : ["Palatino", "New Century Schoolbook", "Bookman", "Computer Modern Roman"],
    #"figure.constrained_layout.use" : True,
    "figure.autolayout" : True,
    "figure.titlesize" : 20,
    "axes.labelsize" : 20,
    "legend.fontsize" : 15,
    "xtick.labelsize" : 15,
    "xtick.major.size" : 5,
    "ytick.labelsize" : 15,
    "ytick.major.size" : 5,
    })


S_list = ["40", "60", "80", "125", "175", "250", "350", "500", "700", "1000", "1400", "2000", "2800", "4000"] # ["1000", "4000"]
c = "2.00"
mu = "0.00" # ["2.00", "3.00"] # ["1.00", "1.50", "2.00", "2.50", "3.00"]
sigma_list = ["1.00"] # ["1.30", "1.50", "1.80", "3.00"]
tau = "-0.60"
N_ext = 300
N_previous_ext = 'This_is_the_kind_of_job_in_which_computers_are_way_way_better_than_me'
pattern = 'Log_Initial_Condition_Extraction_*_Measure_1.txt'
N_measures = 1

for sigma in sigma_list:
    for S in S_list:
        dir_path = Path.cwd() / f'./mu_{mu}_sigma_{sigma}/S_{S}_c_{c}/Initial_Conditions/'
        dir_path.mkdir(parents=True, exist_ok=True, mode = 0o777)
        N_previous_ext = count_pattern_in_dir(dir_path, pattern)
        print(f'mu_{mu}_sigma_{sigma}_tau_{tau} ---> {N_previous_ext}')
        output_file = Path.cwd() / f"Output.txt"
        with open(output_file, 'a') as of:
            #print(f"Stiamo per eseguire:\n./Ecosystem_Extract_and_Evolve {S} {sigma} {mu} {N_ext} {N_previous_ext} {N_measures}\n")     
            EaE_process = subprocess.run(
                ['./Ecosystem_Extract_and_Evolve', f'{S}', f'{sigma}', f'{mu}', f'{N_ext}', f'{N_previous_ext}', f'{N_measures}'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True # Altrimenti è in binario
            )
            print(f'mu_{mu}_sigma_{sigma}_tau_{tau}/S_{S}_c_{c} DONE\n')
            #print(EaE_process.returncode)
            of.write(f'\nOutput:\n {EaE_process.stdout}\n')
            #print(f'\nOutput:\n {EaE_process.stdout}\n')
            of.write(f'\nError: \n {EaE_process.stderr}\n')
            #print(f'\nError: \n {EaE_process.stderr}\n')