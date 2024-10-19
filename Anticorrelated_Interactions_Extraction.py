import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from pathlib import Path

def multiline_comment(comment):
    pass

def count_pattern_in_dir(directory, pattern='*'):
    count = 0
    for item in directory.iterdir():
        if (item.is_file() and item.match(pattern)):
            count += 1
    return count

def mean_standard_deviation(x):
    return np.std(x) / np.sqrt(x.size)

def correlated_sparse_random_matrix(S = 1000, c = 2, sigma = 1., tau=0, mean=[0, 0]):
    p_c = c/(S-1)
    cov = [[sigma, tau*sigma], [tau*sigma, sigma]]
    random_matrix = np.zeros((S,S))
    rng = np.random.default_rng()
    for i in range(S):
        # random_matrix[i][i] = d # Actually we are interested mostly on the non diagonal matrix
        for j in range(i+1, S, 1):
            if(rng.random() < p_c):
                random_matrix[i][j], random_matrix[j][i] = rng.multivariate_normal(mean, cov)
    return random_matrix

def save_sparse_matrix_txt(sparse_matrix, matrix_file_path):
    with open(matrix_file_path, 'w') as mf:
        for i in range(sparse_matrix.row.size):
            mf.write(f"{sparse_matrix.row[i]}\t{sparse_matrix.col[i]}\t{sparse_matrix.data[i]}\n")


S_list = [40, 60, 80, 125, 175, 250, 350, 500, 700, 1000, 1400, 2000, 2800, 4000] #[125, 250, 500, 1000] // {40, 60, 80, 125, 175, 250, 350, 500, 700, 1000, 1400, 2000, 2800, 4000};
c_str = "2.00"
c = float(c_str)
mu = "0.00"
sigma_list = ["1.00"]
tau = -3./5.
N_ext = 300 # 0 
N_prev = 'This_is_the_kind_of_job_in_which_computers_are_way_way_better_than_me'
pattern = 'Extraction_*_Interaction_Matrix_Sparse.txt'

for sigma in sigma_list:
    sigma_path = Path.cwd() / f"mu_{mu}_sigma_{sigma}" # f"sigma_{sigma}"
    #q_path.mkdir(parents=True, exist_ok=True, mode = 0o777)
    for S in S_list:
        s_path = sigma_path / f"S_{S}_c_{c_str}"
        #s_path.mkdir(parents=True, exist_ok=True, mode = 0o777)
        extractions_path = s_path / "Interactions_Extractions"
        extractions_path.mkdir(parents=True, exist_ok=True, mode = 0o777)
        N_prev = count_pattern_in_dir(extractions_path, pattern)
        print(f'mu_{mu}_sigma_{sigma}_tau_{tau} ---> {N_prev}')
        for ext in range(N_prev, N_prev+N_ext):
            csrm = correlated_sparse_random_matrix(S = S, c = c, sigma = float(sigma), tau=tau)
            sparse_matrix = coo_matrix(csrm)
            #print(sparse_matrix)
            matrix_file = extractions_path / f"Extraction_{ext+1}_Interaction_Matrix_Sparse.txt" 
            save_sparse_matrix_txt(sparse_matrix, matrix_file)
            print(f"S={S}, sigma={sigma}, ext={ext+1}")
