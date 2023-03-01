# Code by Elisabeth Holm
# Tracks the pathways of activation
# Python 3.10.4
# Feb 2023 - Present

import numpy as np  # numpy 1.22.4
import math
from collections import OrderedDict
import xlrd  # xlrd 2.0.1 for reading from xls files
from xlwt import Workbook  # xlwt 1.3.0, for writing to xls files
from xlutils.copy import copy  # xlwt 1.3.0
from xlutils.margins import number_of_good_cols  # xlwt 1.3.0
import statistics
from openpyxl.workbook import Workbook as openpyxlWorkbook

# 367 12 383 456  # ROI gilly chose

#file = "New_Vids/Sept_2022_Tests/TS-20220801155100769_391_61_363_354.xls"  # I chose this ROI, clear pathway, ROI near electrode
# later: use data from vid 6363 for all chroms going off at once
# 9896 - clear pathway
file = "New_Vids/Sept_2022_Tests/TS-20220801155100769_367_12_383_456.xls"

centroids = OrderedDict()  # dict with key: ID num, value: list [x,y] of its centroid coordinates
neighborhood = OrderedDict()  # dict with key: ID num, value: set of ID num(s) of its neighbors


# converts xls file to xlsx
def convert_xls_to_xlsx():
    # open file using xlrd (read-only)
    xls_wb = xlrd.open_workbook(file)
    # create a pyxl workbook that we'll copy the data over to
    pyxl_wb = openpyxlWorkbook()

    for i in range(xls_wb.nsheets):
        xlsSheet = xls_wb.sheet_by_index(i)
        sheet = pyxl_wb.active if i == 0 else pyxl_wb.create_sheet()
        sheet.title = xlsSheet.name

        for row in range(xlsSheet.nrows):
            for col in range(xlsSheet.ncols):
                sheet.cell(row=row + 1, column=col + 1).value = xlsSheet.cell_value(row, col)


# Get sheets in the workbook by index
areas_sheet = xls_wb.sheet_by_index(0)
centroids_sheet = xls_wb.sheet_by_index(1)


# get the centroids from the spreadsheet, populates the centroids dict so we can use that instead
def extract_centroids():
    global centroids
    global PIXELS_PER_MM

    # if using an older spreadsheet that doesn't have this data, use the most common val
    try:
        PIXELS_PER_MM = float(centroids_sheet.cell_value(6, 2))
    except:
        PIXELS_PER_MM = 58.3590729166667

    # Put the centroid of each chrom into the centroids dictionary
    for colNumber in range(1, centroids_sheet.ncols):
        cur_ID = centroids_sheet.cell_value(0, colNumber)  # get ID num of cur chrom
        #cur_centroid = centroids_sheet.col_values(colNumber)[1:3]  # get [x,y] centroid of cur chrom
        x = float(centroids_sheet.cell_value(1, colNumber))
        y = float(centroids_sheet.cell_value(2, colNumber))
        cur_centroid = [x, y]

        centroids[cur_ID] = cur_centroid  # populate centroid dict
        neighborhood[cur_ID] = set()  # prepare neighborhood dict to be filled later


extract_centroids()


# go through centroids sheet and determine which chroms are adjacent to each other
# based on a threshold Euclidean distance
def determine_neighbors():
    global neighborhood

    mm_nei_thresh = 3  # TODO play with thresh dist for neighbor

    ID_nums = list(centroids.keys())
    cents = list(centroids.values())

    # iterate through each chrom and find all its neighbors
    for i in range(len(ID_nums)):
        ID_num = ID_nums[i]
        centroid = cents[i]
        # check against unchecked pns
        for j in (range(i + 1, len(ID_nums))):  # pn = potential neighbor
            pn_ID = ID_nums[j]
            pn_centroid = cents[j]
            # if a pn is close enough to the current chrom, count it as a neighbor
            if ((math.dist(centroid, pn_centroid)/PIXELS_PER_MM) < mm_nei_thresh):
                # both chroms as neighbors of each other
                (neighborhood[ID_num]).add(pn_ID)
                neighborhood[pn_ID].add(ID_num)

    print(neighborhood)

determine_neighbors()


# compute the rolling variance time series for a single chrom
# if unable to compute variance (e.g not enough data points in window), say var = 0
# params: chromatophore ID num (e.g "C12"), the entire column of areas data (as a list) for chrom
def compute_var_data(chrom_ID, area_data):
    var_data = {}
    rolling_window_size = 20

    # generate variance data point for each
    for i in range(len(area_data)):
        # fill in first window-size frames as having a previous window variance of 0 (since there
        # aren't enough points to compute a variance
        if i < rolling_window_size - 1:
            window_var = 0
        else:
            # take window of data as the previous [window size] points before the ith point
            window_data = area_data[i - rolling_window_size: i - 1]

            # clear out empty strings from window data (that's how empty cells appear)
            while '' in window_data:
                window_data.remove('')

            if len(window_data) > 1:  # check that we have at least 2 data points, compute variance
                window_var = statistics.variance(window_data)  # rolling variance
                window_avg = statistics.mean(window_data)  # rolling avg
            else:  # if the window covers a gap in the data, put the variance as 0
                window_var = 0


        var_data[i] = window_var

    #print(var_data)
    #print()
    return var_data


# determine activation time(s) for a single chrom
# param: the entire column of variation data (as a dict of frame_num:var)) for that chromatophore
def find_activations(chrom_ID, var_data):
    # iterate through each variance data point
    for frame_num in range(len(var_data.keys())):
        var_pt = var_data[frame_num]

        # TODO
        # take avg slope between prev 3 data points -- testing frame 4, we'd avg slopes between
        # 1/2 2/3 3/4
        # if that avg slope exceeds 500 (or some threshold amount) then that counts as an event


# save the variance data to a new sheet in the spreadsheet
def save_variance_data(all_var_data):
    # copy the contents of xls file (since we can't really edit an existing spreadsheet,
    # we have to make a copy, edit that copy, then save that copy)
    new_workbook = copy(xls_wb)
    try:
        varSheet = new_workbook.add_sheet('Variances')
    except:
        varSheet = new_workbook.add_sheet('Variances')

    varSheet.write(0, 0, "frame")

    numFrames = len(list(all_var_data.values())[0])

    # add column in column 0 that says each frame number
    for frame_num in range(numFrames):
        varSheet.write(frame_num + 1, 0, frame_num)

    curCol = 1

    # iterate through each identified chromatophore
    for id_num, var_data in all_var_data.items():
        # add ID num on row 0
        varSheet.write(0, curCol, id_num)

        # add var time series in each column
        for frame_num, var in var_data.items():
            varSheet.write(frame_num + 1, curCol, var)

        curCol += 1

    # save the file (under the same name, thus replacing the unedited file)
    new_workbook.save(file)


# returns active frame #s (frames when chrom is active) and
# event frame #s (first time a chrom is activated)
def find_activations_min_max_method(chrom_ID, area_data):
    just_area_nums = area_data.copy()

    # clear out empty strings from window data (that's how empty cells appear)
    while '' in just_area_nums:
        just_area_nums.remove('')

    min_a = min(just_area_nums)
    max_a = max(just_area_nums)
    size_range = max_a - min_a
    percent_active_thresh = .05  # once chrom has reached this percent of its size range, count as active
    thresh_area = size_range * percent_active_thresh

    active_frames = []
    act_event_frames = []
    deact_event_frames = []

    for i in range(len(area_data)):
        # skip empty spreadsheet cells
        if area_data[i] == '':
            # if it loses the chrom, count that as deactivation # TODO check if this is what we want to do
            if i-1 in active_frames:
                deact_event_frames.append(i)
            continue

        # if chrom is active in this frame
        if (area_data[i] - min_a) >= thresh_area:
            active_frames.append(i)
            # if this is a frame when it switches from inactive to active
            if (i > 0 and (i-1 not in active_frames)) or i == 0:
                act_event_frames.append(i)
            # if it's still active
            if i == len(area_data) - 1:
                deact_event_frames.append("active at end of vid")

        # if this is a frame when it switches from active to inactive
        elif (i > 0 and (i-1 in active_frames)):
            deact_event_frames.append(i)

    return active_frames, act_event_frames, deact_event_frames


# for every chrom, determine its times of activation -- save in a binary version of the areas sheet
def make_activation_dataset():
    all_var_data = {}  # dict with key: chrom ID, value: that chrom's variance time series (a dict)

    for colNumber in range(1, areas_sheet.ncols):
        chrom_col = areas_sheet.col_values(colNumber)  # get column as a list
        chrom_ID = chrom_col[0]

        var_data = compute_var_data(chrom_ID, chrom_col[1:])  # create var time series for chrom
        all_var_data[chrom_ID] = var_data

        activation_frames = find_activations(chrom_ID, var_data)

        # doing the % size range method
        active_frames, act_event_frames, deact_event_frames = find_activations_min_max_method(chrom_ID, chrom_col[1:])
        #print(chrom_ID)
        #print(active_frames)
        #print(act_event_frames)
        #print(deact_event_frames)
        #print()

    #save_variance_data(all_var_data)

        # TODO make activation dataset from variance versus time dataset

make_activation_dataset()


def save_activation_data():
    # copy the contents of xls file (since we can't really edit an existing spreadsheet,
    # we have to make a copy, edit that copy, then save that copy)
    new_workbook = copy(xls_wb)

    activSheet = new_workbook.add_sheet('activation')

    # TODO populate the new spreadsheet with the activation data

    # save the file
    new_workbook.save(file)


#def identify_pathways():
    #print("wahoo")

#def predict_independence():
    #print("bat time")

'''
reader = csv.DictReader(open('bats.csv'))

sum_G1 = 0
sum_G2 = 0
sum_G3 = 0
sum_G4 = 0
sum_G5 = 0
sum_T = 0
num_bats = 0

sum_G1_given_T = 0
sum_G2_given_T = 0
sum_G3_given_T = 0
sum_G4_given_T = 0
sum_G5_given_T = 0

for bat in reader:
    num_bats += 1

    if bat['G1'].lower() == "true":
        sum_G1 += 1
    if bat['G2'].lower() == "true":
        sum_G2 += 1
    if bat['G3'].lower() == "true":
        sum_G3 += 1
    if bat['G4'].lower() == "true":
        sum_G4 += 1
    if bat['G5'].lower() == "true":
        sum_G5 += 1
    if bat['T'].lower() == "true":
        sum_T += 1
        if bat['G1'].lower() == "true":
            sum_G1_given_T += 1
        if bat['G2'].lower() == "true":
            sum_G2_given_T += 1
        if bat['G3'].lower() == "true":
            sum_G3_given_T += 1
        if bat['G4'].lower() == "true":
            sum_G4_given_T += 1
        if bat['G5'].lower() == "true":
            sum_G5_given_T += 1

prob_G1 = sum_G1/num_bats
prob_G2 = sum_G2/num_bats
prob_G3 = sum_G3/num_bats
prob_G4 = sum_G4/num_bats
prob_G5 = sum_G5/num_bats
prob_T = sum_T/num_bats

print("P(GX) -- Probability of each event independently happening (individual events/total num of bats)")
print("G1: " + str(prob_G1))
print("G2: " + str(prob_G2))
print("G3: " + str(prob_G3))
print("G4: " + str(prob_G4))
print("G5: " + str(prob_G5))
print("T:  " + str(prob_T))
print()

p_G1_given_T = sum_G1_given_T/num_bats
p_G2_given_T = sum_G2_given_T/num_bats
p_G3_given_T = sum_G3_given_T/num_bats
p_G4_given_T = sum_G4_given_T/num_bats
p_G5_given_T = sum_G5_given_T/num_bats

print("P(GX|T) -- Probability of GX appearing given T is true (sum_GX_given_T/num_bats)")
print("G1|T: " + str(p_G1_given_T))
print("G2|T: " + str(p_G2_given_T))
print("G3|T: " + str(p_G3_given_T))
print("G4|T: " + str(p_G4_given_T))
print("G5|T: " + str(p_G5_given_T))
print()

def predict_indep(p_GX, p_GX_given_T):
    difference = abs(p_GX*prob_T - p_GX_given_T)
    print("Difference:" + str(difference))
    if difference < 0.01 :
        print("Result: NOT independent of T")
    else:
        print("Result: independent of T")
    print()

print("Prediction of independence based on difference between P(GX)*P(T) and P(GX|T) (dependent if <0.01 apart from each other)")
print("G1 :")
predict_indep(prob_G1, p_G1_given_T)
print("G2 :")
predict_indep(prob_G2, p_G2_given_T)
print("G3 :")
predict_indep(prob_G3, p_G3_given_T)
print("G4 :")
predict_indep(prob_G4, p_G4_given_T)
print("G5 :")
predict_indep(prob_G5, p_G5_given_T)
'''
