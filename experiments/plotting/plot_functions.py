import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FixedLocator, LogFormatterMathtext
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.patches import Rectangle
from matplotlib.ticker import LogLocator, LogFormatterMathtext, NullLocator
from matplotlib.ticker import FixedLocator, FuncFormatter
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans"],
    "mathtext.fontset": "dejavusans",
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "xtick.minor.size": 2.5,
    "ytick.minor.size": 2.5,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    "legend.frameon": False,
    "legend.framealpha": 1.0,
    "legend.edgecolor": "black",
    "legend.fancybox": False,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
    "ps.fonttype": 42,
    "pdf.fonttype": 42,
    "text.usetex": False,
})

def read_data(name): # Get the step-size sequence from a file
    file = open(name, "r")
    n=len(file.readlines())
    file.close()
    T=np.zeros(n)
    file = open(name, "r")
    for line in file:
        x=line.split()
        T[int(x[-1])]=float(x[-2])
    file.close()
    deltaT=np.zeros(n-1)
    for i in range(n-1):
        deltaT[i] = T[i+1]-T[i]
    return(T,deltaT)

def read_prop(name, epses, y):  # Get the global error for different tolerance values
    error = []
    epses_number = []

    for j in range(len(epses)):
        eps = float(epses[j])
        if len(epses) == 150:
            eps0 = 1000 * eps
        else:
            eps0 = 100 * eps

        eps_string = "{:.6e}".format(eps)
        eps0_string = "{:.6e}".format(eps0)
        try:
            file = open(name.format(eps0_string, eps_string), "r")
    
            for line in file:
                x = line.split()
                y1 = []
    
                for k in range(len(x)):
                    y1.append(float(x[k]))
    
                y1 = np.array(y1)
    
            error.append(np.linalg.norm(y - y1))
            epses_number.append(eps)
            file.close()
        except:
            print(f"Missing file: {name}-{eps_string}")

    return(epses_number, error)

def read_article(name, epses, y):  # Get the global error
    error = []
    epses_number = []

    for j in range(len(epses)):
        eps = float(epses[j])
        eps_string = "{:.6e}".format(eps)

        try:
            file = open(name.format(eps_string), "r")
    
            for line in file:
                x = line.split()
                y1 = []
    
                for k in range(len(x)):
                    y1.append(float(x[k]))
    
                y1 = np.array(y1)
    
            error.append(np.linalg.norm(y-y1))
            epses_number.append(eps)
            file.close()
        except:
            print(f"Missing file: {name}-{eps_string}")

    return(epses_number, error)

def read_sequential(name, epses, y):  # Get the sequential global error
    error = []
    epses_number = []

    for j in range(len(epses)):
        eps = float(epses[j])
        eps_string = "{:.6e}".format(eps)

        try:
            file = open(name.format(eps_string), "r")
    
            for line in file:
                x = line.split()
                y1 = []
    
                for k in range(len(x)):
                    y1.append(float(x[k]))
    
                y1 = np.array(y1)
    
            error.append(np.linalg.norm(y-y1))
            epses_number.append(eps)
            file.close()
        except:
            print(f"Hiányzó fájl: {name}-{eps_string}")

    return(epses_number, error)

def read_prop_r(name, eps0s, y):  # Get the global error for different tolerance values
    error = []
    rs = []
    eps0s2 = []

    for j in range(len(eps0s)):
        eps0 = float(eps0s[j])
        eps0_string = "{:.6e}".format(eps0)
        try:
            file = open(name.format(eps0_string), "r")
    
            for line in file:
                x = line.split()
                y1 = []
    
                for k in range(len(x)):
                    y1.append(float(x[k]))
    
                y1 = np.array(y1)
    
            error.append(np.linalg.norm(y - y1))
            r = int((eps0/10**(-9))**(1/2))
            rs.append(r)
            eps0s2.append(eps0)
            file.close()
        except:
            print(f"Hiányzó fájl: {name.format(eps0_string)}")

    return(rs, eps0s2, error)

def load_data_hermite(result_dir="../fig14-16/results_heatlin-hermite_harmonic"):
    result_dir = Path(result_dir)
    files = sorted(result_dir.glob("r_*_p_*_test_*.dat"))
    records = []

    for filename in files:
        values = np.loadtxt(filename, ndmin=1)
        r = int(values[0])
        P = int(values[1])
        test = int(filename.name[-5])

        grid_time = float(values[2])
        prerefined_time = float(values[3])
        mgrit_time = float(values[4])
        total_time = float(values[5])
        iterations = float(values[6])

        records.append(
            {
                "r": r,
                "P": P,
                "test": test,
                "grid_time": grid_time,
                "prerefined_time": prerefined_time,
                "mgrit_time": mgrit_time,
                "total_time": total_time,
                "iterations": iterations,
                "filename": filename,
            }
        )

    p_values = np.array(sorted({record["P"] for record in records}),dtype=int)
    r_values = np.array(sorted({record["r"] for record in records}),dtype=int)

    p_to_index = {P: index for index, P in enumerate(p_values)}
    r_to_index = {r: index for index, r in enumerate(r_values)}
    
    shape = (len(p_values),len(r_values),)

    grouped_records = {}

    for record in records:
        key = (record["P"],record["r"])

        if key not in grouped_records:
            grouped_records[key] = {
                "grid_time": [],
                "prerefined_time": [],
                "mgrit_time": [],
                "total_time": [],
                "iterations": [],
            }

        grouped_records[key]["grid_time"].append(record["grid_time"])
        grouped_records[key]["prerefined_time"].append(record["prerefined_time"])
        grouped_records[key]["mgrit_time"].append(record["mgrit_time"])
        grouped_records[key]["total_time"].append(record["total_time"])
        grouped_records[key]["iterations"].append(record["iterations"])

    grid_time_array = np.full(shape,np.nan,dtype=float)
    prerefined_time_array = np.full(shape,np.nan,dtype=float)
    mgrit_time_array = np.full(shape,np.nan,dtype=float)
    total_time_array = np.full(shape,np.nan,dtype=float)
    iteration_array = np.full(shape,np.nan,dtype=float)

    for (P, r), measurements in grouped_records.items():
        row = p_to_index[P]
        column = r_to_index[r]

        grid_time_array[row, column] = np.median(measurements["grid_time"])
        prerefined_time_array[row, column] = np.median(measurements["prerefined_time"])
        mgrit_time_array[row, column] = np.median(measurements["mgrit_time"])
        total_time_array[row, column] = np.median(measurements["total_time"])
        iteration_array[row, column] = np.median(measurements["iterations"])

    return (grid_time_array,prerefined_time_array,mgrit_time_array,total_time_array,iteration_array,p_values,r_values)

def load_data_prev(result_dir="../fig14-16/results_heatlin-article_original"):
    result_dir = Path(result_dir)
    files = sorted(result_dir.glob("ntime_*_p_*_test_*.dat"))
    records = []

    for filename in files:
        values = np.loadtxt(filename, ndmin=1)
        P = int(values[1])
        test = int(filename.name[-5])

        total_time = float(values[3])
        iterations = float(values[4])

        records.append(
            {
                "P": P,
                "test": test,
                "total_time": total_time,
                "iterations": iterations,
                "filename": filename,
            }
        )

    p_values = np.array(sorted({record["P"] for record in records}),dtype=int)
    p_to_index = {P: index for index, P in enumerate(p_values)}
    
    shape = (len(p_values))

    grouped_records = {}

    for record in records:
        key = (record["P"])

        if key not in grouped_records:
            grouped_records[key] = {
                "total_time": [],
                "iterations": [],
            }
        grouped_records[key]["total_time"].append(record["total_time"])
        grouped_records[key]["iterations"].append(record["iterations"])

    total_time_array = np.full(shape,np.nan,dtype=float)
    iteration_array = np.full(shape,np.nan,dtype=float)

    for (P), measurements in grouped_records.items():
        row = p_to_index[P]

        total_time_array[row] = np.median(measurements["total_time"])
        iteration_array[row] = np.median(measurements["iterations"])

    return (total_time_array,iteration_array,p_values)

def plot_fig6():
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.6), layout="constrained")

    T, deltaT = read_data("../fig6/results/lin-hermite_harmonic_cut-1.000000e-05-1.000000e-07.txt")
    ax.plot(T[-31:-1], deltaT[-30:], color="C0", label="Truncated presolve step")

    T, deltaT = read_data("../fig6/results/lin-hermite_harmonic-1.000000e-05-1.000000e-07.txt")
    ax.plot(T[-31:-1], deltaT[-30:], color="C1", label="Full presolve step")

    ax.set_yscale("log")
    ax.yaxis.set_major_locator(FixedLocator([1e-2, 2e-2]))
    ax.yaxis.set_major_formatter(FixedFormatter([r"$10^{-2}$", r"$2\times 10^{-2}$"]))
    ax.yaxis.set_minor_formatter(NullFormatter())

    ax.legend(frameon=False, loc="lower left", handlelength=2.2, handletextpad=0.4, labelspacing=0.7)

    fig.savefig("figures/fig6.pdf", format="pdf", bbox_inches="tight")
    plt.show()

def plot_fig7(problem, method):
    fig, ax = plt.subplots(figsize=(3.5, 2.6), layout="constrained")
    if method == "IE":
        T, deltaT = read_data(f"../fig7/results/{problem}-original_{method}-1.000000e-07.txt")
        ax.plot(T[:-2], deltaT[:-1], color="C0", label=r"$r_n$")
        
        T, deltaT = read_data(f"../fig7/results/{problem}-article_{method}-1.000000e-07.txt")
        ax.plot(T[:-2], deltaT[:-1], color="C1", label=r"$r_{mod,n}$")
        
        T, deltaT = read_data(f"../fig7/results/{problem}-{method}-1.000000e-05-1.000000e-07.txt")
        ax.plot(T[:-2], deltaT[:-1], color="C2", label=r"$r$")
        
        T, deltaT = read_data(f"../fig7/results/{problem}-sequential_{method}-1.000000e-07.txt")
        ax.plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black", label=r"Sequential")
    elif method == "BS":
        T, deltaT = read_data(f"../fig7/results/{problem}-original_{method}-1.000000e-09.txt")
        ax.plot(T[:-2], deltaT[:-1], color="C0", label=r"$r_n$")
        
        T, deltaT = read_data(f"../fig7/results/{problem}-article_{method}-1.000000e-09.txt")
        ax.plot(T[:-2], deltaT[:-1], color="C1", label=r"$r_{mod,n}$")
        
        T, deltaT = read_data(f"../fig7/results/{problem}-{method}-1.000000e-06-1.000000e-09.txt")
        ax.plot(T[:-2], deltaT[:-1], color="C2", label=r"$r$")
        
        T, deltaT = read_data(f"../fig7/results/{problem}-sequential_{method}-1.000000e-09.txt")
        ax.plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black", label=r"Sequential")

    ax.set_yscale("log")
    
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
    ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))

    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.legend(frameon=False, loc="lower right", handlelength=2.2, handletextpad=0.4, labelspacing=0.5)
    if method == "IE":
        fig.savefig(f"figures/fig7_{problem}-{method}-1e-7.pdf", format="pdf", bbox_inches="tight")
    else:
        fig.savefig(f"figures/fig7_{problem}-{method}-1e-9.pdf", format="pdf", bbox_inches="tight")
    plt.show()

def plot_fig8(name1, method, epses):  # Plot the tolerance - global error function
    f = open("../fig8/exact_sols/{}-exact.txt".format(name1),"r")

    y = f.read().split()
    y = [float(numeric_string) for numeric_string in y]
    y = np.array(y)
    f.close()

    fig, ax = plt.subplots(figsize=(3.5, 2.6), layout="constrained")
    [tol, error] = read_article("../fig8/results/{}-article_{}-{}.txt".format(name1, method, "{}"),epses,y)
    ax.plot(tol, error, label=r"$r_n$")
    
    [tol, error] = read_article("../fig8/results/{}-original_{}-{}.txt".format(name1, method, "{}"),epses,y)
    ax.plot(tol, error, label=r"$r_{mod,n}$")

    [tol, error] = read_prop("../fig8/results/{}-{}-{}-{}.txt".format(name1, method, "{}", "{}"),epses,y)
    ax.plot(tol, error, label=r"r")

    [tol, error] = read_sequential("../fig8/results/{}-sequential_{}-{}.txt".format(name1, method, "{}"),epses,y)
    ax.plot(tol, error, "--", linewidth=1, color = "black", label="Sequential")

    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.minorticks_on()
    ax.legend(frameon=False, loc="lower right", handlelength=2.2, handletextpad=0.4, labelspacing=0.5)
    if method == "IE":
        fig.savefig(f"figures/fig8_{name1}-{method}.pdf", format="pdf", bbox_inches="tight")
    else:
        fig.savefig(f"figures/fig8_{name1}-{method}.pdf", format="pdf", bbox_inches="tight")
    plt.show()

def plot_fig9():
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 5.2), layout="constrained", sharex=True, sharey=True)
    T, deltaT = read_data(f"../fig9/results/vdp-hermite_harmonic_I-1.000000e-03-1.000000e-05.txt")
    axs[0,0].plot(T[:-2], deltaT[:-1], color="C0", label=r"TOL=$10^{-5}$")
    T, deltaT = read_data(f"../fig9/results/vdp-hermite_harmonic_I-1.000000e-05-1.000000e-07.txt")
    axs[0,0].plot(T[:-2], deltaT[:-1], color="C1", label=r"TOL=$10^{-7}$")
    T, deltaT = read_data(f"../fig9/results/vdp-sequential_I-1.000000e-05.txt")
    axs[0,0].plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black", label="Sequential")
    T, deltaT = read_data(f"../fig9/results/vdp-sequential_I-1.000000e-07.txt")
    axs[0,0].plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black")
    axs[0,0].set_yscale("log")
    axs[0,0].minorticks_on()
    axs[0,0].legend(frameon=False, loc="lower right", handlelength=2.2, handletextpad=0.4, labelspacing=0.5)
    
    T, deltaT = read_data(f"../fig9/results/vdp-hermite_harmonic_exp-1.000000e-03-1.000000e-05.txt")
    axs[0,1].plot(T[:-2], deltaT[:-1], color="C0", label=r"TOL=$10^{-5}$")
    T, deltaT = read_data(f"../fig9/results/vdp-hermite_harmonic_exp-1.000000e-05-1.000000e-07.txt")
    axs[0,1].plot(T[:-2], deltaT[:-1], color="C1", label=r"TOL=$10^{-7}$")
    T, deltaT = read_data(f"../fig9/results/vdp-sequential_exp-1.000000e-05.txt")
    axs[0,1].plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black", label="Sequential")
    T, deltaT = read_data(f"../fig9/results/vdp-sequential_exp-1.000000e-07.txt")
    axs[0,1].plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black")
    axs[0,1].set_yscale("log")
    axs[0,1].minorticks_on()
    axs[0,1].legend(frameon=False, loc="lower right", handlelength=2.2, handletextpad=0.4, labelspacing=0.5)

    T, deltaT = read_data(f"../fig9/results/vdp-hermite_harmonic-1.000000e-03-1.000000e-05.txt")
    axs[1,0].plot(T[:-2], deltaT[:-1], color="C0", label=r"TOL=$10^{-5}$")
    T, deltaT = read_data(f"../fig9/results/vdp-hermite_harmonic-1.000000e-05-1.000000e-07.txt")
    axs[1,0].plot(T[:-2], deltaT[:-1], color="C1", label=r"TOL=$10^{-7}$")
    T, deltaT = read_data(f"../fig9/results/vdp-sequential-1.000000e-05.txt")
    axs[1,0].plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black", label="Sequential")
    T, deltaT = read_data(f"../fig9/results/vdp-sequential-1.000000e-07.txt")
    axs[1,0].plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black")
    axs[1,0].set_yscale("log")
    axs[1,0].minorticks_on()
    axs[1,0].legend(frameon=False, loc="lower right", handlelength=2.2, handletextpad=0.4, labelspacing=0.5)

    T, deltaT = read_data(f"../fig9/results/vdp-hermite_harmonic_PID-1.000000e-03-1.000000e-05.txt")
    axs[1,1].plot(T[:-2], deltaT[:-1], color="C0", label=r"TOL=$10^{-5}$")
    T, deltaT = read_data(f"../fig9/results/vdp-hermite_harmonic_PID-1.000000e-05-1.000000e-07.txt")
    axs[1,1].plot(T[:-2], deltaT[:-1], color="C1", label=r"TOL=$10^{-7}$")
    T, deltaT = read_data(f"../fig9/results/vdp-sequential_PID-1.000000e-05.txt")
    axs[1,1].plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black", label="Sequential")
    T, deltaT = read_data(f"../fig9/results/vdp-sequential_PID-1.000000e-07.txt")
    axs[1,1].plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black")
    axs[1,1].set_yscale("log")
    axs[1,1].minorticks_on()
    axs[1,1].legend(frameon=False, loc="lower right", handlelength=2.2, handletextpad=0.4, labelspacing=0.5)
    
    
    fig.savefig("figures/fig9.pdf", format="pdf", bbox_inches="tight")
    plt.show()

def plot_fig10():
    epses = []
    for i in range(200):
        epses.append(10**(-9+i/50))
    f = open("../fig8/exact_sols/vdp-exact.txt","r")

    y = f.read().split()
    y = [float(numeric_string) for numeric_string in y]
    y = np.array(y)
    f.close()

    fig, ax = plt.subplots(figsize=(3.5, 2.6), layout="constrained")
    [tol, error] = read_prop("../fig10/results/vdp-hermite_harmonic_I-{}-{}.txt".format("{}", "{}"),epses,y)
    ax.plot(tol, error, label=r"I-controller")
    [tol, error] = read_sequential("../fig10/results/vdp-sequential_I-{}.txt".format("{}"),epses,y)
    ax.plot(tol, error, "--", linewidth=1, color = "black", label="Sequential-I")

    [tol, error] = read_prop("../fig10/results/vdp-hermite_harmonic_exp-{}-{}.txt".format("{}", "{}"),epses,y)
    ax.plot(tol, error, label=r"Exp-controller")
    [tol, error] = read_sequential("../fig10/results/vdp-sequential_exp-{}.txt".format("{}"),epses,y)
    ax.plot(tol, error, ":", linewidth=1, color = "black", label="Sequential-exp")

    [tol, error] = read_prop("../fig10/results/vdp-hermite_harmonic-{}-{}.txt".format("{}", "{}"),epses,y)
    ax.plot(tol, error, label=r"PI-controller")
    [tol, error] = read_sequential("../fig10/results/vdp-sequential-{}.txt".format("{}"),epses,y)
    ax.plot(tol, error, "-.", linewidth=1, color = "black", label="Sequential-PI")

    [tol, error] = read_prop("../fig10/results/vdp-hermite_harmonic_PID-{}-{}.txt".format("{}", "{}"),epses,y)
    ax.plot(tol, error, label=r"PID-controller")
    [tol, error] = read_sequential("../fig10/results/vdp-sequential_PID-{}.txt".format("{}"),epses,y)
    ax.plot(tol, error, linestyle=(0, (8, 2, 2, 2)), linewidth=1, color = "black", label="Sequential-PID")

    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.legend(frameon=False, loc="lower right", handlelength=2.2, handletextpad=0.4, labelspacing=0.5)
    fig.savefig("figures/fig10.pdf", format="pdf", bbox_inches="tight")
    plt.show()

def plot_fig11():
    fig, ax = plt.subplots(figsize=(3.5, 2.6), layout="constrained")
    T, deltaT = read_data(f"../fig11/results/lin-linear-1.000000e-03-1.000000e-05.txt")
    ax.plot(T[:-2], deltaT[:-1], color="C0", label=r"Uniform")
    T, deltaT = read_data(f"../fig11/results/lin-hermite_harmonic-1.000000e-03-1.000000e-05.txt")
    ax.plot(T[:-2], deltaT[:-1], color="C1", label=r"Hermite")
    T, deltaT = read_data(f"../fig11/results/lin-sequential-1.000000e-05.txt")
    ax.plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black", label=r"Sequential")
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
    ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.legend(frameon=False, loc="lower right", handlelength=2.2, handletextpad=0.4, labelspacing=0.7)
    fig.savefig("figures/fig11_lin.pdf", format="pdf", bbox_inches="tight")
    plt.show()

    fig, ax = plt.subplots(figsize=(3.5, 2.6), layout="constrained")
    T, deltaT = read_data(f"../fig11/results/LV-linear-1.000000e-03-1.000000e-05.txt")
    ax.plot(T[:-2], deltaT[:-1], color="C0", label=r"Uniform")
    T, deltaT = read_data(f"../fig11/results/LV-hermite_harmonic-1.000000e-03-1.000000e-05.txt")
    ax.plot(T[:-2], deltaT[:-1], color="C1", label=r"Hermite")
    T, deltaT = read_data(f"../fig11/results/LV-sequential-1.000000e-05.txt")
    ax.plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black", label=r"Sequential")
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
    ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_yticks([2e-2])
    ax.set_yticklabels([r"$2\cdot10^{-2}$"])
    ax.legend(frameon=False, loc="lower right", handlelength=2.2, handletextpad=0.4, labelspacing=0.7)
    fig.savefig("figures/fig11_LV.pdf", format="pdf", bbox_inches="tight")
    plt.show()

    fig, ax = plt.subplots(figsize=(3.5, 2.6), layout="constrained")
    T, deltaT = read_data(f"../fig11/results/vdp-linear-1.000000e-03-1.000000e-05.txt")
    ax.plot(T[:-2], deltaT[:-1], color="C0", label=r"Uniform")
    T, deltaT = read_data(f"../fig11/results/vdp-hermite_harmonic-1.000000e-03-1.000000e-05.txt")
    ax.plot(T[:-2], deltaT[:-1], color="C1", label=r"Hermite")
    T, deltaT = read_data(f"../fig11/results/vdp-sequential-1.000000e-05.txt")
    ax.plot(T[:-2], deltaT[:-1], "--", linewidth=1, color="black", label=r"Sequential")
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
    ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.legend(frameon=False, loc="lower right", handlelength=2.2, handletextpad=0.4, labelspacing=0.7)
    fig.savefig("figures/fig11_vdp.pdf", format="pdf", bbox_inches="tight")
    plt.show()

def plot_fig12(name,maxi):
    maxi_plot = maxi
    [T,deltaT] = read_data("{}".format(name))
    if maxi == None:
        maxi = 0
        for i in range(1,len(deltaT)-2):
            x = deltaT[i]/deltaT[i-1]
            y = deltaT[i+1]/deltaT[i]
            if x>maxi:
                maxi = x + 0.05
            if 1/x>maxi:
                maxi = 1/x + 0.05
            if y>maxi:
                maxi = y + 0.05
            if 1/y>maxi:
                maxi = 1/y + 0.05

    # Sűrűbb rács a simább ábra érdekében
    xs = np.exp(np.linspace(np.log(1/maxi), np.log(maxi), 1000))
    ys = np.exp(np.linspace(np.log(1/maxi), np.log(maxi), 1000))
    
    X, Y = np.meshgrid(xs, ys)
    
    h0 = 1.0
    h1 = X * h0
    h2 = Y * h1
    
    d0 = 2 * h0 * h1 / (h0 + h1)
    d1 = 2 * h1 * h2 / (h1 + h2)
    
    # A der1 feltétele
    condition_der1 = (
        (d0 + d1 < 2 * h1)
        | (2 * d0 + d1 < 3 * h1)
        | (d0 + 2 * d1 < 3 * h1)
        | (d0**2 + d0 * d1 + d1**2 - 6 * h1 * (d0 + d1) + 9 * h1**2 < 0)
    )
    
    der1 = np.where(condition_der1, 1, -1)
    
    # A der2 feltétele
    der2 = np.zeros_like(X, dtype=int)
    der2[(2 * d0 + d1 < 3 * h1) & (3 * h1 < d0 + 2 * d1)] = 1
    der2[(d0 + 2 * d1 < 3 * h1) & (3 * h1 < 2 * d0 + d1)] = -1
    
    region = np.empty_like(X, dtype=int)
    region[(der1 == 1) & (der2 == 1)] = 0
    region[(der1 == 1) & (der2 == -1)] = 1
    region[(der1 == 1) & (der2 == 0)] = 2
    region[der1 == -1] = 3
    
    # Jól elkülönülő, nyomtatásban is kellemes színek
    colors = ["#0072B2", "#D55E00", "#009E73", "#D1495B"]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, 4.5, 1), cmap.N)
    
    # Egyhasábos cikkábrához megfelelő méret
    fig, ax = plt.subplots(figsize=(3.8, 3.5), layout="constrained")
    
    ax.pcolormesh(
        X,
        Y,
        region,
        cmap=cmap,
        norm=norm,
        shading="auto",
        rasterized=True,
    )
    
    # A tartományok határvonalai
    ax.contour(
        X,
        Y,
        region,
        levels=[0.5, 1.5, 2.5],
        colors="black",
        linewidths=0.65,
    )
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    
    ax.set_xlabel(r"$x = H_{n}/H_{n-1}$")
    ax.set_ylabel(r"$y = H_{n+1}/H_n$")
    
    ax.set_xlim(xs.min(), xs.max())
    ax.set_ylim(ys.min(), ys.max())
    ax.set_aspect("equal", adjustable="box")
    
    # Segédvonalak
    ax.axvline(x=1, color="black", linestyle="--", linewidth=1.0)
    ax.axhline(y=1, color="black", linestyle="--", linewidth=1.0)
    ax.plot(xs, xs, color="black", linestyle=":", linewidth=1.2)
    
    # Beosztások mind a négy oldalon
    ax.tick_params(
        which="both",
        top=True,
        right=True,
        labeltop=False,
        labelright=False,
    )
    
    # Saját jelmagyarázat
    legend_elements = [
        Patch(facecolor=colors[0], edgecolor="black", linewidth=0.4,
              label=r"$\gamma'' > 0$"),
        Patch(facecolor=colors[1], edgecolor="black", linewidth=0.4,
              label=r"$\gamma'' < 0$"),
        Patch(facecolor=colors[2], edgecolor="black", linewidth=0.4,
              label=r"$\gamma'' change$"),
    ]
    
    
    [T,deltaT] = read_data("{}".format(name))
    for i in range(1,len(deltaT)-2):
        x = deltaT[i]/deltaT[i-1]
        y = deltaT[i+1]/deltaT[i]
        ax.plot(x,y,"*", color = "purple")
    
    ax.legend(
        handles=legend_elements,
        loc="lower right",
        borderpad=0.5,
        labelspacing=0.35,
        handlelength=1.4,
        handleheight=0.8,
        )
    ticks = [0.8, 0.9, 1.0, 1.1, 1.2]

    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(FixedLocator(ticks))
        axis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    
    fig.savefig(f"figures/fig12-{name[17:len(name)-4]}.pdf", format="pdf", bbox_inches="tight")
    plt.show()

def plot_fig13(names):
    colors = ["C0", "C1", "C2"]

    eps0s = []
    for i in range(3):
        for k in range(2,11):
            r = k*10**i
            eps0s.append(10**(-9)*r*r)

    fig, axs = plt.subplots(1, len(names), figsize=(3.5*len(names), 2.6), layout="constrained", sharex=True)
    axs = np.atleast_1d(axs)

    for i, name in enumerate(names):
        f = open("../fig8/exact_sols/{}-exact.txt".format(name), "r")
        y = f.read().split()
        y = [float(numeric_string) for numeric_string in y]
        y = np.array(y)
        f.close()

        [rs, eps0s2, error] = read_prop_r("../fig13/results/{}-hermite_harmonic-{}-1.000000e-09.txt".format(name, "{}"), eps0s, y)

        axs[i].plot(eps0s2, error, color=colors[i], marker="o", markersize=2.0)
        axs[i].set_xscale("log")
        axs[i].set_yscale("log")
        axs[i].minorticks_on()

    fig.savefig("figures/fig13.pdf", format="pdf", bbox_inches="tight")
    plt.show()

def plot_fig14(r, p_values_hermite, total_time_hermite, p_values_article, total_time_article, p_values_original, total_time_original):
    if r <= 10:
        index = r-2
    elif r<=100:
        index = int(r/10+7)
    else:
        index = int(r/100+16)
    fig, ax = plt.subplots(figsize=(3.5, 2.6), layout="constrained")
    ax.plot(p_values_original, total_time_original, marker="o", markersize=2.5, color = "C0", label = r"$r_n$")
    ax.plot(p_values_article, total_time_article, marker="o", markersize=2.5, color = "C1", label = r"$r_{mod,n}$")
    ax.plot(p_values_hermite, total_time_hermite[:,index], marker="o", markersize=2.5, color = "C2", label = r"$r$")
    ax.plot(p_values_hermite, 3.15967748*np.ones_like(p_values_hermite), "--", linewidth=1, color = "black", label = "Sequential")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xticks(p_values_hermite, labels=p_values_hermite)
    ax.legend(frameon=False, loc="lower left", handlelength=2.2, handletextpad=0.4, labelspacing=0.5)
    fig.savefig("figures/fig14.pdf", format="pdf", bbox_inches="tight")
    plt.show()

def plot_fig15(rs, p_values, total_time_array, mgrit_time_array, prerefined_time_array, grid_time_array):
    fig, axs = plt.subplots(1, len(rs), figsize=(3.5*len(rs), 2.6), layout="constrained", sharex=True, sharey=True)
    for i in range(len(rs)):
        r = rs[i]
        if r <= 10:
            index = r-2
        elif r<=100:
            index = int(r/10+7)
        else:
            index = int(r/100+16)
        axs[i].plot(p_values, grid_time_array[:,index], marker="o", markersize=2.5, color = "C0", label = "Presolve")
        axs[i].plot(p_values, prerefined_time_array[:,index], marker="o", markersize=2.5, color = "C1", label = "Prerefine")
        axs[i].plot(p_values, mgrit_time_array[:,index], marker="o", markersize=2.5, color = "C2", label = "MGRIT")
        axs[i].plot(p_values, total_time_array[:,index], marker="o", markersize=2.5, color = "C3", label = "Total")
        axs[i].plot(p_values, 3.15967748*np.ones_like(p_values), "--", linewidth=1, color = "black", label = "Sequential")
        axs[i].set_xscale("log", base=2)
        axs[i].set_yscale("log", base=2)
        axs[i].set_xticks(p_values, labels=p_values)
        axs[i].legend(frameon=False, loc="lower left", handlelength=2.2, handletextpad=0.4, labelspacing=0.3)
    fig.savefig("figures/fig15.pdf", format="pdf", bbox_inches="tight")
    plt.show()

def plot_fig16(Ps, r_values, total_time_array, mgrit_time_array, prerefined_time_array, grid_time_array, mgrit_iteration_array):
    fig, axs = plt.subplots(1, len(Ps), figsize=(3.5*len(Ps), 2.6), layout="constrained")
    axs = np.atleast_1d(axs)

    for i in range(len(Ps)):
        index = int(np.log2(Ps[i]))

        mgrit_times = mgrit_time_array[index,:]
        total_times = total_time_array[index,:]
        mgrit_iterations = mgrit_iteration_array[index,:]

        line1, = axs[i].plot(r_values, mgrit_times, marker="o", markersize=1.5, color="C2", label="MGRIT")
        line2, = axs[i].plot(r_values, total_times, marker="o", markersize=1.5, color="C3", label="Total")

        axs[i].set_xscale("log", base=10)
        axs[i].set_yscale("log", base=10)
        axs[i].xaxis.set_major_formatter(ScalarFormatter())

        ax_iter = axs[i].twinx()
        line3, = ax_iter.step(r_values, mgrit_iterations, where="mid", color="C0", linestyle="--", marker="o", markersize=2.0, label="MGRIT iterations")

        valid = np.isfinite(mgrit_times) & (mgrit_times > 0.0) & np.isfinite(mgrit_iterations)
        valid_times = mgrit_times[valid]
        valid_iterations = mgrit_iterations[valid]

        if valid_times.size > 0:
            min_time = np.min(valid_times)
            max_time = np.max(valid_times)
            min_iteration = np.min(valid_iterations)
            max_iteration = np.max(valid_iterations)
            time_axis_min, time_axis_max = axs[i].get_ylim()

            if max_time > min_time and max_iteration > min_iteration and time_axis_max > time_axis_min:
                log_axis_min = np.log10(time_axis_min)
                log_axis_max = np.log10(time_axis_max)

                min_position = (np.log10(min_time)-log_axis_min)/(log_axis_max-log_axis_min)
                max_position = (np.log10(max_time)-log_axis_min)/(log_axis_max-log_axis_min)

                iteration_axis_range = (max_iteration-min_iteration)/(max_position-min_position)
                iteration_axis_min = min_iteration-min_position*iteration_axis_range
                iteration_axis_max = iteration_axis_min+iteration_axis_range

                ax_iter.set_ylim(iteration_axis_min, iteration_axis_max)

            elif max_iteration == min_iteration:
                ax_iter.set_ylim(min_iteration-1, max_iteration+1)

        ax_iter.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax_iter.tick_params(axis="y", colors="C0")

        lines = [line1, line2, line3]
        labels = [line.get_label() for line in lines]
        axs[i].legend(lines, labels, frameon=False, loc="best", handlelength=2.2, handletextpad=0.4, labelspacing=0.7)

    fig.savefig("figures/fig16.pdf", format="pdf", bbox_inches="tight")
    plt.show()
