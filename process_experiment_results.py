# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker

def process_input(line, results):
    # split line on whitespaces
    line = line.split()
    a = ""
    
    

    # first element is path to file, split on '/'
    relevant_path = line[0][line[0].find('GAMBIT'):]
    path = relevant_path.split('/')
    gambit = path[0]
    runs_completed = int(line[1].split("/")[0])
    
    currDataSettings = f"{path[0]}__{path[3]}__{path[4]}__{path[5]}"
    
    if(path[3] == "F5"):
        a = path[4].split("_")[1]
        if(float(a) > 5.5) or (float(a) <= 4.125):
            # From results, we can see that both gambit_R and gambit_K can solve F5 with a>=5.5 for all selected problem sizes, so not interesting to plot anymore.
            # Similarly, below 4.125 neither can solve all selected problem sizes, so not interesting to plot anymore.
            a = ""

    runs_completed = int(line[1].split("/")[0])
    attempts = int(line[1].split("/")[1]) # TODO doublecheck this is correct (and following if statement)
    # lines that didn't complete at least 29 runs (but executed more than 1 run), just add empty tuple (if not already in dict)
    # check line length of 5 instead of 4, since if 0/0 happens, it does print avgEvals: , but no value (so len 5)
    if(len(line) <= 5):
        if(attempts > 0 and (currDataSettings not in results or (results[currDataSettings][0] is None and runs_completed > results[currDataSettings][3]))):
            results[currDataSettings] = (None, None, None, runs_completed)
        return (gambit, a)
    
    # lines that have less than 30 attempts should be ignored in results, but can be added as none to show they were attempted
    if(attempts < 30):
        if(currDataSettings not in results or (results[currDataSettings][0] is None and runs_completed > results[currDataSettings][3])):
            results[currDataSettings] = (None, None, None, runs_completed)
        return (gambit, a)

    currPopsize = int(line[3])
    currAvgEvals = float(line[5])
    popsize_mult = float(path[1].split("_")[1])
    allEvals = list(map(float, line[7:]))
    if((currDataSettings not in results) or (results[currDataSettings][0] is None) or (currAvgEvals < results[currDataSettings][1]) ): # (currPopsize < results[currDataSettings][0]) ):
        results[currDataSettings] = (currPopsize, currAvgEvals, popsize_mult, runs_completed, allEvals)
    return (gambit, a)
    
        


def get_ys_F1_F4(results, currF, x, gambit):
    fractions = ["0.25", "0.5", "0.75"]
    ys = [[-1]*len(x) for _ in range(len(fractions))]
    popsizes = [[-1]*len(x) for _ in range(len(fractions))]
    popsizeMults = [[-1]*len(x) for _ in range(len(fractions))]
    for key in results:
        settings = key.split("__")
        if(settings[1] == currF):
            if(results[key][0] is None or settings[0] != gambit):
                continue
            totalvars = sum(map(int, settings[3].split("_")))
            frac = settings[2]
            ys[fractions.index(frac)][x.index(totalvars)] = results[key][1]
            popsizes[fractions.index(frac)][x.index(totalvars)] = results[key][0]
            popsizeMults[fractions.index(frac)][x.index(totalvars)] = results[key][2]
    return (ys, popsizes, popsizeMults)

def get_ys_F5(results, x, gambit, a_values, currF="F5"):
    ys = [[0]*len(x) for _ in range(len(a_values))]
    popsizes = [[0]*len(x) for _ in range(len(a_values))]
    popsizeMults = [[0]*len(x) for _ in range(len(a_values))]
    for key in results:
        settings = key.split("__")
        if(settings[1] == currF):
            if(results[key][0] is None or settings[0] != gambit):
                continue
            totalvars = sum(map(int, settings[3].split("_")))
            a = settings[2].split("_")[1]
            if(a not in a_values):
                # skip results with a value not in a_values
                continue
            ys[a_values.index(a)][x.index(totalvars)] = results[key][1]
            popsizes[a_values.index(a)][x.index(totalvars)] = results[key][0]
            popsizeMults[a_values.index(a)][x.index(totalvars)] = results[key][2]
    return (ys, popsizes, popsizeMults)

def make_plot(results, currF, title, fn, gambit):
    fig, (ax1, ax2) = plt.subplots(2, figsize=(4, 7.5)) #plt.subplots(3, figsize=(6.4, 12.8))
    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    fig.suptitle(title)
    x = [40, 80, 120, 160]
    ys, popsizes, popsizeMults = get_ys_F1_F4(results, currF, x, gambit)
    ax2.loglog(x, ys[0], 'ro-', label='f=0.25', color='red', clip_on=False)
    ax2.loglog(x, ys[1], 'ro-', label='f=0.5', color='green', clip_on=False)
    ax2.loglog(x, ys[2], 'ro-', label='f=0.75', color='blue', clip_on=False)
    ax2.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.set_xticks(x)
    ax2.ticklabel_format(axis='x', style='plain')
    ax2.set_xlabel('Problem Size')
    ax2.set_ylabel('Number of Evaluations')
    ax2.set_xlim([40, 160])
    if currF in ["F1", "F2"]:
        ax2.set_ylim([10**3, 10**6])
    elif currF in ["F3", "F4"]:
        ax2.set_ylim([10**4, 10**6])
    ax2.legend()

    ax1.loglog(x, popsizes[0], 'ro-', label='f=0.25', color='red', linestyle='-', alpha=0.5, clip_on=False)
    ax1.loglog(x, popsizes[1], 'ro-', label='f=0.5', color='green', linestyle='--', alpha=0.5, clip_on=False)
    ax1.loglog(x, popsizes[2], 'ro-', label='f=0.75', color='blue', linestyle='-.', alpha=0.5, clip_on=False)
    ax1.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useMathText=False))  # Disable scientific notation for y-axis
    # ax1.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax1.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.set_xticks(x)
    ax1.ticklabel_format(axis='x', style='plain')
    # ax1.ticklabel_format(axis='y', style='plain')
    ax1.set_xlabel('Problem Size')
    ax1.set_ylabel('Population size')
    ax1.set_xlim([40, 160])
    if currF in ["F1", "F2"]:
        ax1.set_ylim([10, 100])
        ax1.set_yticks([10, 20, 30, 40, 60, 80, 100])
    elif currF in ["F3", "F4"]:
        ax1.set_ylim([10, 1000])
        ax1.set_yticks([10, 100, 1000])
        ax1.set_yticklabels(map(str, [10, 100, 1000]))


    ax1.legend()

    # ax0.loglog(x, popsizeMults[0], 'ro-', label='f=0.25', color='red', linestyle='-', alpha=0.5)
    # ax0.loglog(x, popsizeMults[1], 'ro-', label='f=0.5', color='green', linestyle='--', alpha=0.5)
    # ax0.loglog(x, popsizeMults[2], 'ro-', label='f=0.75', color='blue', linestyle='-.', alpha=0.5)
    # ax0.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax0.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    # ax0.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax0.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    # ax0.set_xticks(x)
    # ax0.ticklabel_format(style='plain')
    # ax0.set_xlabel('Problem size')
    # ax0.set_ylabel('Popsize multiplier')
    # ax0.legend()
    # plt.show()
    fig.savefig(f"plots/remake_for_thesis/{fn}.png", bbox_inches='tight', dpi=600)


def make_F5_plot(results, currF, title, fn, gambit, a_values):
    colors = ['red', 'green', 'blue', 'orange', 'purple', 'brown', 'cyan', 'olive', 'gray']
    fig, (ax0, ax1, ax2) = plt.subplots(3, figsize=(6.4, 12.8))
    fig.suptitle(title)
    x = [20, 40, 60]
    ys, popsizes, popsizeMults = get_ys_F5(results, x, gambit, a_values, currF)
    for i in range(len(a_values)):
        ax2.loglog(x, ys[i], 'ro-', label=f'a={a_values[i]}', color=colors[i%len(colors)])
    ax2.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.set_xticks(x)
    ax2.ticklabel_format(axis='x', style='plain')
    ax2.set_xlabel('Problem Size')
    ax2.set_ylabel('Number of Evaluations')
    ax2.legend()

    for i in range(len(a_values)):
        ax1.loglog(x, popsizes[i], 'ro-', label=f'a={a_values[i]}', linestyle=['-', "--", "-.", ":"][i%4], linewidth=8-(6*i/len(a_values)), alpha=0.5, color=colors[i%len(colors)])
    ax1.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.set_xticks(x)
    ax1.ticklabel_format(style='plain')
    ax1.set_xlabel('Problem Size')
    ax1.set_ylabel('Population size')
    ax1.legend()

    for i in range(len(a_values)):
        ax0.loglog(x, popsizeMults[i], 'ro-', label=f'a={a_values[i]}', linestyle=['-', "--", "-.", ":"][i%4], linewidth=8-(6*i/len(a_values)), alpha=0.5, color=colors[i%len(colors)])
    ax0.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax0.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax0.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax0.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax0.set_xticks(x)
    ax0.ticklabel_format(style='plain')
    ax0.set_xlabel('Problem size')
    ax0.set_ylabel('Popsize multiplier')
    ax0.legend()
    # plt.show()
    fig.savefig(f"plots/remake_for_thesis/{fn}.png", bbox_inches='tight', dpi=300)

def make_alt_F5_plot(results, currF, title, fn, gambit, a_values):
    fig, (ax1, ax2) = plt.subplots(2)
    fig.suptitle(title)
    x = [20, 40, 60]
    ys, popsizes, _ = get_ys_F5(results, x, gambit, a_values, currF)
    print(f"ys: {ys}")
    print(f"popsizes: {popsizes}")
    f5_colors = ['#66c2a5', '#fc8d62', '#8da0cb']
    xaxis_data = [[float(a_val) for idx, a_val in enumerate(a_values) if ys[idx][0] > 0],
                    [float(a_val) for idx, a_val in enumerate(a_values) if ys[idx][1] > 0],
                    [float(a_val) for idx, a_val in enumerate(a_values) if ys[idx][2] > 0]]
    ax1.plot(xaxis_data[0], [row[0] for row in ys if row[0] > 0], 'ro-', label='ps=20', color=f5_colors[0])
    ax1.plot(xaxis_data[1], [row[1] for row in ys if row[1] > 0], 'ro-', label='ps=40', color=f5_colors[1])
    ax1.plot(xaxis_data[2], [row[2] for row in ys if row[2] > 0], 'ro-', label='ps=60', color=f5_colors[2])
    # ax1.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax1.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    # ax1.set_xticks(a_values)
    # ax1.ticklabel_format(axis='x', style='plain')

    ax1.set_xlabel('a value')
    ax1.set_ylabel('Number of Evaluations')
    ax1.legend()

    xaxis_data = [[float(a_val) for idx, a_val in enumerate(a_values) if popsizes[idx][0] > 0],
                    [float(a_val) for idx, a_val in enumerate(a_values) if popsizes[idx][1] > 0],
                    [float(a_val) for idx, a_val in enumerate(a_values) if popsizes[idx][2] > 0]]
    ax2.plot(xaxis_data[0], [row[0] for row in popsizes if row[0] > 0], 'ro-', label='ps=20', color=f5_colors[0], linestyle='-', alpha=0.5)
    ax2.plot(xaxis_data[1], [row[1] for row in popsizes if row[1] > 0], 'ro-', label='ps=40', color=f5_colors[1], linestyle='--', alpha=0.5)
    ax2.plot(xaxis_data[2], [row[2] for row in popsizes if row[2] > 0], 'ro-', label='ps=60', color=f5_colors[2], linestyle='-.', alpha=0.5)
    # ax2.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax2.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    # ax2.set_xticks(a_values)
    # ax2.ticklabel_format(axis='x', style='plain')

    ax2.set_xlabel('a value')
    ax2.set_ylabel('Population size')
    ax2.legend()

    # plt.show()
    fig.savefig(f"plots/remake_for_thesis/{fn}.png", bbox_inches='tight', dpi=300)


def make_all_plots(results, gambit, a_values):
    # # Map gambit_K and gambit_R to names in line with what I will use in thesis (original and modified, respectively)
    # if(gambit == "GAMBIT_K"):
    #     title_gambit = "Original"
    # elif(gambit == "GAMBIT_R"):
    #     title_gambit = "Modified"

    # title_gambit = f"{gambit}: "
    title_gambit = ""

    # F1
    make_plot(results, "F1", f'{title_gambit}F1 OneMaxSphere', f"{gambit}_F1_OneMaxSphere", gambit)
    # F2
    make_plot(results, "F2", f'{title_gambit}F2 OneMaxEllipse', f"{gambit}_F2_OneMaxEllipse", gambit)
    # F3
    make_plot(results, "F3", f'{title_gambit}F3 TrapSphere', f"{gambit}_F3_TrapSphere", gambit)
    # F4
    make_plot(results, "F4", f'{title_gambit}F4 TrapEllipse', f"{gambit}_F4_TrapEllipse", gambit)

    # F5 alt
    make_alt_F5_plot(results, "F5", f'{gambit}: F5 Trap Block Ellipse', f"{gambit}_F5_TrapBlockEllipse_alt", gambit, a_values)



    # # F5
    # make_F5_plot(results, "F5", f'{gambit}: F5 Trap Block Ellipse', f"{gambit}_F5_TrapBlockEllipse", gambit, a_values)
    # plt.figure()
    # x = [20, 40, 60]
    # ys = [[20114, -1, -1],   # a = 1.15
    #     [21818, 53865, 100150],    # a = 2.75  
    #     [21530, 67267, 108328]]    # a = 3

    # plt.loglog(x, ys[0], 'ro-', label='a=1.15', color='red')
    # plt.loglog(x, ys[1], 'ro-', label='a=2.75', color='green')
    # plt.loglog(x, ys[2], 'ro-', label='a=3.00', color='blue')
    # plt.xticks(x)
    # plt.title('F5 Trap Block Ellipse')
    # plt.xlabel('Problem Size')
    # plt.ylabel('Number of Evaluations')
    # plt.legend()
    # plt.show()
    return

def main():
    # read in text file line by line
    results = {}
    a_values_R = []
    a_values_K = []
    with open('plots/experiment4/experiment4_results.txt', 'r') as f:
        lines = f.readlines()
        for line in lines:
            (whichgambit, a_val) = process_input(line, results)
            if(whichgambit == "GAMBIT_R" and a_val != "" and a_val not in a_values_R):
                a_values_R.append(a_val)
            if(whichgambit == "GAMBIT_K" and a_val != "" and a_val not in a_values_K):
                a_values_K.append(a_val)

    a_values_R = list(map(str, sorted(map(float, a_values_R))))
    a_values_K = list(map(str, sorted(map(float, a_values_K))))

    print(f"Starting to make plots...")
    make_all_plots(results, "GAMBIT_R", a_values_R)
    print(f"Done with GAMBIT_R, starting GAMBIT_K plots...")
    make_all_plots(results, "GAMBIT_K", a_values_K)

    results_list = []
    for key in results:
        if(key.split("__")[1] != "F5"):
            continue
        results_list.append(f"{key}: ({results[key][0]}, {results[key][1]}, {results[key][2]}, {results[key][3]})")
    print(*results_list, sep="\n")
    

if __name__ == "__main__":
    main()

