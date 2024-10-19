import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from pathlib import Path

def multiline_comment(comment):
    pass

def mean_standard_deviation(x):
    return np.std(x) / rng.sqrt(x.size)

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

def correlated_sparse_random_matrix_with_extras(S = 1000, c = 2, sigma = 1., tau=0, mean=[0, 0]):
    p_c = c/(S-1)
    cov = [[sigma, tau*sigma], [tau*sigma, sigma]]
    x = np.zeros(2*c*S)
    y = np.zeros(2*c*S)
    cnt = 0
    random_matrix = np.zeros((S,S))
    rng = np.random.default_rng()
    for i in range(S):
        # random_matrix[i][i] = d # Actually we are interested mostly on the non diagonal matrix
        for j in range(i+1, S, 1):
            if(random.random() < p_c):
                a, b = rng.multivariate_normal(mean, cov)
                random_matrix[i][j], random_matrix[j][i] = (a, b)
                x[cnt] = a
                y[cnt] = b
                cnt = cnt + 1
    x = np.resize(x, cnt)
    y = np.resize(y, cnt)
    return random_matrix, x, y

def save_sparse_matrix_txt(sparse_matrix, matrix_file_path):
    with open(matrix_file_path, 'w') as mf:
        for i in range(sparse_matrix.row.size):
            mf.write(f"{sparse_matrix.row[i]}\t{sparse_matrix.col[i]}\t{sparse_matrix.data[i]}\n")

def draw_axis(ax):
    xlim = ax.get_xlim()
    x_line = np.linspace(xlim[0], xlim[1], 10)
    ylim = ax.get_ylim()
    y_line = np.linspace(ylim[0], ylim[1], 10)
    zero_line = np.zeros(10)
    ax.plot(x_line, zero_line, label=None, linewidth=1.0, color="#000000", zorder=0)
    ax.plot(zero_line, y_line, label=None, linewidth=1.0, color="#000000", zorder=0)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

plt.rcParams.update({
    "text.usetex": True,
    "font.family" : "serif",
    "font.serif" : ["Palatino", "New Century Schoolbook", "Bookman", "Computer Modern Roman"],
    #"figure.constrained_layout.use" : True,
    "figure.autolayout" : True,
    "figure.titlesize" : 20,
    "axes.labelsize" : 25,
    "legend.fontsize" : 14,
    "xtick.labelsize" : 15,
    "xtick.major.size" : 5,
    "ytick.labelsize" : 15,
    "ytick.major.size" : 5,
    })

S_list = [40, 60, 80, 125, 175, 250, 350, 500, 700, 1000, 1400, 2000, 2800, 4000] #[125, 250, 500, 1000] // {40, 60, 80, 125, 175, 250, 350, 500, 700, 1000, 1400, 2000, 2800, 4000};
c = 2
sigma_list = ["0.90"]
tau = -3./5.
N_ext = 1 #200
N_prev = 'This_is_the_kind_of_job_in_which_computers_are_way_way_better_than_me'
pattern = 'Extraction_*_Interaction_Matrix_Sparse.txt'
graph_format = 'pdf'

x = np.array([])
y = np.array([])

for sigma in sigma_list:
    sigma_path = Path.cwd() / f"sigma_{sigma}" # f"sigma_{sigma}"
    #q_path.mkdir(parents=True, exist_ok=True, mode = 0o777)
    for S in S_list:
        s_path = sigma_path / f"S_{S}"
        #s_path.mkdir(parents=True, exist_ok=True, mode = 0o777)
        extractions_path = s_path / "Extractions"
        extractions_path.mkdir(parents=True, exist_ok=True, mode = 0o777)
        N_prev = count_pattern_in_dir(extractions_path, pattern)
        print(f'mu_{mu}_sigma_{sigma}_tau_{tau} ---> {N_prev}')
        eig_path = s_path / "Eigenvalues"
        eig_path.mkdir(parents=True, exist_ok=True, mode = 0o777)
        graphs_path = s_path / "Graphs"
        graphs_path.mkdir(parents=True, exist_ok=True, mode = 0o777)
        fig, ax = plt.subplots() # figsize=(6.4, 4.8), constrained_layout=True, dpi=250
        for ext in range(N_prev, N_prev+N_ext):
            csrm, tmp_x, tmp_y = correlated_sparse_random_matrix_with_extras(S = S, c = c, sigma = float(sigma), tau=tau)
            x = np.concatenate((x, tmp_x))
            y = np.concatenate((y, tmp_y))
            sparse_matrix = coo_matrix(csrm)
            matrix_file = extractions_path / f"Extraction_{ext+1}_Interaction_Matrix_Sparse.txt" 
            save_sparse_matrix_txt(sparse_matrix, matrix_file)
            print(f"S={S}, sigma={sigma}, ext={ext+1}")
            w = np.linalg.eigvals(csrm)
            ax.plot(np.real(w), np.imag(w), marker = ".", markersize = 2.0, linestyle = 'None', color='#bfff00') # Mattone '#ff6633' - Lime '#bfff00'
            eig_file = eig_path / f"Eigenvalues_Interaction_{ext+1}.npy"
            np.save(eig_file, w)        
            print(ext+1)
        draw_axis(ax)
        #ax.set_xlim(, )
        ax.set_xlabel("$\\Re(\\lambda)$") #fontsize=20
        ax.set_ylabel("$\\Im(\\lambda)$") #fontsize=20
        fig.suptitle(f"$S = {S} \\;,\\; \\sigma = {sigma} \\;,\\; \\tau = {tau}$")
        #graphs_path = graphs_path / f"Anticorrelated_Interaction_Spectra_S{S}.{graph_format}"
        fig.savefig(graphs_path / f"Anticorrelated_Interaction_Spectra_S{S}.{graph_format}")
        plt.close(fig)

        h, edgex, edgey = np.histogram2d(x, y, bins=[100, 100], density=True)
        edgex = edgex[:-1]
        edgey = edgey[:-1]
        xs, ys = np.meshgrid(edgex, edgey)
        # Filled contour
        fig, ax = plt.subplots() # figsize=(6.4, 4.8), constrained_layout=True, dpi=250
        filled_contour = ax.contourf(xs, ys, h, cmap='hot') # Per specificare colormap --> cmap = 'jet', 'viridis', 'hot', 'bone'
        ax.set_xlabel("$\\alpha_{ik}$") #fontsize=20
        ax.set_ylabel("$\\alpha_{ki}$") #fontsize=20
        fig.suptitle(f"$S = {S} \\;,\\; \\sigma = {sigma} \\;,\\; \\tau = {tau}$")
        fig.colorbar(filled_contour)
        fig.savefig(graphs_path / f"Anticorrelated_Distribution_S{S}.{graph_format}")
        plt.close(fig)
