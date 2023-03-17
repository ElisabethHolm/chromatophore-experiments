# Code by Elisabeth Holm
# Analyzes chromatophore data to connect
# temporal and spatial chromatophore data
# Python 3.10.4
# March 2023 - Present

import math
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser  # argparse 1.4.0
from gooey import Gooey, GooeyParser  # Gooey 1.0.8.1
import numpy as np  # numpy 1.22.4
import pandas as pd  # pandas 1.5.3
import matplotlib.pyplot as plt  # matplotlib 3.5.2

filepath = "New_Vids/Sept_2022_Tests/TS-20220801155100769_367_12_383_456.xlsx"  # clear chrom pathway, ROI near electrode
chrom_IDs = []  # list of chrom IDs in order of how they appear on spreadsheet, populate in main()
args = []
window_size = 30
n_df = None  # df from neighborhood sheet (populated in load())
a_df = None  # df from activations sheet  (populated in load())


@Gooey  # The GUI Decorator goes here
def parse_args():
    global args

    ap = GooeyParser(formatter_class=ArgumentDefaultsHelpFormatter,
                     conflict_handler='resolve')

    ap.add_argument("-s", "--spreadsheet", default=filepath, widget="FileChooser",
                    help="File path of the spreadsheet (.xlsx file generated from pathways.py)")
    ap.add_argument("-w", "--window_size", default=window_size, type=int,
                    help="Size of the window to divide frames by")

    args = vars(ap.parse_args())


# load the relevant sheets from the file
def load(filename):
    # neighborhood spreadsheet
    n_df = pd.read_excel(f"{filename}", sheet_name="Neighbors")
    a_df = pd.read_excel(f"{filename}", sheet_name="Activations")

    return n_df, a_df


# set global nums based on args
def set_global_nums():
    global filepath
    global n_df
    global a_df
    global chrom_IDs
    global window_size

    filepath = args["spreadsheet"]
    window_size = args["window_size"]
    n_df, a_df = load(filepath)  # load the relevant sheets from the file
    chrom_IDs = a_df.loc[:,"Chromatophore ID"].values.tolist()  # get list of chrom IDs


# sort all activations chronologically
# return 2D list of activations in chronological order with items [chrom_ID, activation_frame]
# e.g [["C12", 23], ["C45", 28], ...]
def get_sorted_activations():
    time_activations = []
    activations = {}

    # iterate through each chromatophore, populate sorted_activations list
    for i in range(len(chrom_IDs)):
        row = a_df.iloc[i]
        chrom_ID = row[0]
        chrom_activations = []
        # iterate through all activation frames
        for a_i in range(1, len(row), 3):
            a_frame = row[a_i]
            # if this is an activation (not nan)
            if not row.isnull()[a_i]:
                time_activations.append([chrom_ID, a_frame])
                chrom_activations.append(a_frame)

        activations[chrom_ID] = chrom_activations

    # sort chronologically
    time_activations = sorted(time_activations, key=lambda l:l[1])

    return time_activations, activations


# sort activation frame num into windows
# return a dataframe with column names as chromatophore and rows as windows of frame nums
def make_act_histogram(activations):
    n_frames = len(pd.read_excel(f"{filepath}", sheet_name="Areas"))  # number of frames in vid
    n_windows = int(n_frames/window_size) + 1
    windows = np.linspace(0, n_frames, n_windows)  # number will be start of window (0-29, 30-59, etc)
    data = {"window start": windows}  # initialize data with windows as first column

    # iterate through chronological activations
    for chrom_ID, act_list in activations.items():
        chrom_bin_hist_data = [0]*len(windows)  # initialize a list of 0's

        # iterate through each activation, switch bit of whatever window it happens in to 1
        for act in act_list:
            window_i = int(act // window_size)  # index of window that activation falls under
            chrom_bin_hist_data[window_i] = 1  # indicate this chrom activates in this window

        data[chrom_ID] = chrom_bin_hist_data  # store column of data in data

    hist_a_df = pd.DataFrame(data)

    plot_act_histogram(hist_a_df, window_size)

    return hist_a_df


# plot histogram of activations
def plot_act_histogram(hist_a_df, window_length):
    num_windows = len(hist_a_df)
    windows = []
    num_window_activations = [0] * num_windows

    # iterate through all rows of hist_a_df
    for i in range(num_windows):
        row = hist_a_df.iloc[i]  # data for this window
        windows.append(str(int(row[0])) + "-" + str(int(row[0]) + window_length - 1))  # bar label
        win_data = row[1:]
        num_window_activations[i] = win_data.sum()  # number of chroms that activated in this window

    # plot graph
    plt.bar(windows, num_window_activations)
    plt.title('Number of chromatophore activations per frame window', fontsize=14)
    plt.xlabel('Range of Frames', fontsize=14)
    plt.ylabel('Number of Activations', fontsize=14)
    plt.grid(True)
    plt.show()


# P(A_i AND A_j)
# = # windows where A_i activates AND A_j activates / # windows
def get_p_A_i_and_A_j(df, ID_i, ID_j):
    rows_Ai_and_Aj = df.loc[(df[ID_i] == 1) & (df[ID_j] == 1)]
    count_Ai_and_Aj = len(rows_Ai_and_Aj)
    count_windows = len(df)

    p_A_i_and_A_j = count_Ai_and_Aj / count_windows
    return p_A_i_and_A_j


# P(A_i)
# = # windows where A_i activates / # windows
def get_p_A_i(df, ID_i):
    rows_Ai = df.loc[df[ID_i] == 1]
    count_Ai = len(rows_Ai)
    count_windows = len(df)

    p_A_i = count_Ai / count_windows
    return p_A_i


# return bool whether a given pair is dependent from one another
# in temporal activation data
# (if True, their activations happen somewhat at the same time usually)
def pair_time_dependent(hist_a_df, ID_i, ID_j):
    p_A_i_and_A_j = get_p_A_i_and_A_j(hist_a_df, ID_i, ID_j)  # P(A_i AND A_j)
    p_A_i = get_p_A_i(hist_a_df, ID_i)  # P(A_i)
    p_A_j = get_p_A_i(hist_a_df, ID_j)  # P(A_j)

    # if P(A_i AND A_j) = P(A_i)P(A_j), pair is independent
    difference = abs(p_A_i_and_A_j - p_A_i*p_A_j)
    # also don't count pair as dep if one of them never activated
    if difference < 0.01 or p_A_i*p_A_j == 0:
        dep = False  # they are independent
    else:
        dep = True  # they are dependent

    # print calculations for ID_i and ID_j being independent based on activation time(s)
    # print("Chrom_i:", ID_i, ", Chrom_j:", ID_j)
    # print("P(A_i):", p_A_i)
    # print("P(A_j):", p_A_j)
    # print("P(A_i AND A_j):", p_A_i_and_A_j)
    # print("P(A_i)*P(A_j):", p_A_i*p_A_j)
    # print("Difference:", difference)
    # print("Dependent:", dep)
    # print()

    return dep


# find all pairs that are time dependent on each other
# return dict with key chrom_ID : value {chrom IDs that this one is dependent with}
def get_time_dep_pairs(hist_a_df):
    time_dep = {}
    num_time_dep_pairs = 0

    # initialize time_dep dict with empty sets
    for chrom_ID in chrom_IDs:
        time_dep[chrom_ID] = set()

    # print("If A_i and A_j are independent, we expect P(A_i AND A_j) = P(A_i) * P(A_j)")
    # print()

    # check every pair
    for i in range(len(chrom_IDs)):
        ID_i = chrom_IDs[i]
        for j in (range(i + 1, len(chrom_IDs))):
            ID_j = chrom_IDs[j]

            # if pair is time dependent on each other
            if pair_time_dependent(hist_a_df, ID_i, ID_j):
                # add both chroms as dependent of each other
                time_dep[ID_i].add(ID_j)
                time_dep[ID_j].add(ID_i)
                num_time_dep_pairs += 1  # keep track of how many time_dep_pairs for P(T) later

    return time_dep, num_time_dep_pairs


# get neighborhood as a dict with key chrom_ID: value {neighbors}
def get_neighborhood():
    neighborhood = {}

    # initialize neighborhood dict with empty sets
    for chrom_ID in chrom_IDs:
        neighborhood[chrom_ID] = set()

    # iterate through all rows in Neighbors spreadsheet
    for i in range(len(chrom_IDs)):
        row = n_df.iloc[i]
        chrom_ID = row[0]
        # iterate through all neighbors
        for n_i in range(1, len(row)):
            # if the value isn't NaN
            if not row.isnull()[n_i]:
                neighborhood[chrom_ID].add(row[n_i])

    return neighborhood


# P(T AND N)
# # dependent neighbor pairs / # pairs
def get_p_T_and_N(time_dep, neighbors, num_pairs):
    num_T_and_N = 0

    # check every pair
    for i in range(len(chrom_IDs)):
        ID_i = chrom_IDs[i]
        for j in (range(i + 1, len(chrom_IDs))):
            ID_j = chrom_IDs[j]

            # if ID_j is both time-dependent and a neighbor of ID_i
            if ID_j in time_dep[ID_i] and ID_j in neighbors[ID_i]:
                num_T_and_N += 1

    p_T_and_N = num_T_and_N / num_pairs
    print("P(T AND N) = # pairs that are neighbors and time-dependent / # pairs")
    print("=", num_T_and_N, "/", num_pairs, "=", p_T_and_N)
    print()
    return p_T_and_N


# P(T)
# # time-dependent pairs / # pairs
def get_p_T(num_time_dep_pairs, num_pairs):
    p_T = num_time_dep_pairs / num_pairs
    print("P(T) = # time-dependent pairs / # pairs")
    print("=", num_time_dep_pairs, "/", num_pairs, "=", p_T)
    print()
    return p_T


# P(N)
# # neighboring pairs / # pairs
def get_p_N(neighborhood, num_pairs):
    unique_neighboring_pairs = set()

    # iterate through neighborhood, get number of unique pairs
    for ID_i, neighbors in neighborhood.items():
        for N_ID in neighbors:
            pair = (ID_i, N_ID)
            # if this pair isn't already in unique pairs (accounting for order not mattering)
            if (ID_i, N_ID) not in unique_neighboring_pairs \
                    and (N_ID, ID_i) not in unique_neighboring_pairs:
                # add the pair to unique pairs
                unique_neighboring_pairs.add(pair)

    num_neighboring_pairs = len(unique_neighboring_pairs)

    p_N = num_neighboring_pairs / num_pairs

    print("P(N) = # neighboring pairs / # pairs")
    print("=", num_neighboring_pairs, "/", num_pairs, "=", p_N)
    print()

    return p_N


# return bool whether T and N are dependent
# in temporal activation data
# (if True, T and N are dependent on one another)
def T_N_dependent(time_dep, neighborhood, num_time_dep_pairs):
    num_pairs = math.comb(len(chrom_IDs), 2)  # num chroms choose 2 = number of possible pairs

    print("T:", time_dep)
    print()
    print("N:", neighborhood)

    p_T_and_N = get_p_T_and_N(time_dep, neighborhood, num_pairs)  # P(T AND N)
    p_T = get_p_T(num_time_dep_pairs, num_pairs)  # P(T)
    p_N = get_p_N(neighborhood, num_pairs)  # P(N)

    # if P(T AND N) = P(T) * P(N), T and N are independent
    difference = abs(p_T_and_N - p_T*p_N)

    # print calculations for T and N being independent
    print("If T and N are independent, we expect P(T AND N) = P(T) * P(N)")
    print("P(T):", p_T)
    print("P(N):", p_N)
    print("P(T AND N):", p_T_and_N)
    print("P(T)*P(N):", p_T*p_N)
    print("Difference:", difference)
    print()

    if difference < 0.01:
        dep = False  # they are independent
    else:
        dep = True  # they are dependent

    return dep


def main():
    parse_args()  # parse input from user
    set_global_nums()  # set global nums based on user input

    # get activations from spreadsheet
    time_activations, activations = get_sorted_activations()

    # sort activation times into windows
    hist_a_df = make_act_histogram(activations)

    # find time-dependent pairs
    time_dep, num_time_dep_pairs = get_time_dep_pairs(hist_a_df)

    # get neighborhood from spreadsheet
    neighborhood = get_neighborhood()

    # calc if T and N are dependent
    T_N_dep = T_N_dependent(time_dep, neighborhood, num_time_dep_pairs)

    # display results
    if T_N_dep:
        print("T and N are dependent.")
    else:
        print("T and N are independent.")

if __name__ == "__main__":
    main()
