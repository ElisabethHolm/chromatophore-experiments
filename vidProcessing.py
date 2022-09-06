# Code by Elisabeth Holm
# Calculates the changing areas of chromatophores from raw video data
# using computer vision (openCV)
# Python 3.10.4
# Spring 2022 - Present

# so I don't have to type it out every time:
# cd Downloads/Coding/Chromatophore_Research/dataProcessing
# python3 vidProcessing.py

# import the necessary packages
import imutils  # imutils 0.5.4
import cv2  # opencv-python 4.6.0.66
import numpy as np  # numpy 1.22.4
from scipy.spatial import distance as dist  # scipy 1.8.1
from collections import OrderedDict
import matplotlib.pyplot as plt  # matplotlib 3.5.2
import argparse  # argparse 1.4.0
from xlwt import Workbook  # xlwt 1.3.0

# these values will be set in setGlobalNums() based on cmd arguments
WINDOW_WIDTH = 1000
(ORIGINAL_WIDTH, ORIGINAL_HEIGHT) = (1920, 1080)
CHROM_PIXEL_DIAMETER = 103  # assumes zoom of 7 and WINDOW_WIDTH of 1000 pix
PIXELS_PER_MM = 1
MAX_DISAPPEARED = 30

nextObjectID = 0
disappeared = OrderedDict()  # objectID : num of frames it has disappeared for
matchedCentroids = OrderedDict()  # objectID : centroid
chromAreas = OrderedDict()  # objectID : area time series
initialCentroids = OrderedDict()  # object ID : centroid when it was first detected

# directory of raw video to process (relative to where this program file is)
#curVidPath = "15_41.811_12TEA+TTX.m4v"  # oldest vid, good w option B or C
curVidPath = "New_Vids/April_22/VID00062.AVI"  # main test -- new, no muscle, great w/ option C (outlines clear but hollow)
#curVidPath = "New_Vids/April_22/VID00062_cropped.mov"  # smaller portion of frame from vid in line above, used for testing program on small scale
#curVidPath = "New_Vids/April_22/45_4-20-22.AVI" # new, decent w/ option B (reduced noise but can't identify big one), or E with (7, 7)
#curVidPath = "New_Vids/April_22/VID00065.AVI"  # pulled from Puneeta on 6/23 (not from 6/23 experiment)
#curVidPath = "New_Vids/April_22/VID00056.AVI"  # pulled from Puneeta on 6/23 (not from 6/23 experiment), struggles w keeping ID nums
#curVidPath = "New_Vids/June_22/TS-20220630153941965.avi"  # from June lab visit

ap = argparse.ArgumentParser()
ap.add_argument("-z", "--zoom", default=7, help="zoom level of microscope (7, 10, 12.5, 16, 20, 25, 32, 40, 50, 63, or 90)")
ap.add_argument("-w", "--window_width", default=WINDOW_WIDTH, help="desired pixel width of frame")
ap.add_argument("-x", "--original_width", default=1920, help="width of original frame in pixels")
ap.add_argument("-y", "--original_height", default=1080, help="height of original frame in pixels")
ap.add_argument("-v", "--vid", default=curVidPath, help="file path of video")

args = vars(ap.parse_args())


# sets the global constants for
# window width and chromatophore pixel diameter based on
# arguments
def setGlobalNums():
	global WINDOW_WIDTH
	global CHROM_PIXEL_DIAMETER
	global PIXELS_PER_MM
	global ORIGINAL_WIDTH
	global ORIGINAL_HEIGHT
	global curVidPath
	global args

	WINDOW_WIDTH = int(args["window_width"])
	ORIGINAL_WIDTH = int(args["original_width"])
	ORIGINAL_HEIGHT = int(args["original_height"])
	zoom = float(args["zoom"])
	curVidPath = args["vid"]

	windowHeight = WINDOW_WIDTH * (ORIGINAL_HEIGHT/ORIGINAL_WIDTH)
	zoom7ChromHeightRatio = 5.5/30  # ratio of radius of expanded chrom over height of window

	# calculate ratio of chromatophore radius to height of window
	# if you want to use the conditional set above, take this line out
	chromToHeightRatio = zoom7ChromHeightRatio * (zoom / 7)

	CHROM_PIXEL_DIAMETER = int(windowHeight * chromToHeightRatio)

	# ensure CHROM_PIXEL_DIAMETER is odd (requirement for param it's passed into)
	if CHROM_PIXEL_DIAMETER % 2 == 0:
		CHROM_PIXEL_DIAMETER += 1

	# pixels per mm adjusting for zoom and window width
	# Note: at x50 zoom, number of pixels for 1 mm is 800.353 on 1920x1080 display
	PIXELS_PER_MM = (zoom / 50) * (800.353 / 1920) * WINDOW_WIDTH


# masks the current frame from the video to identify only the chromatophores
# returns the masked frame
def maskChromatophores(curFrame):
	# blue, green, red = cv2.split(curFrame)  # could use for doing things based on color

	# convert frame to grayscale
	gray_frame = cv2.cvtColor(curFrame, cv2.COLOR_BGR2GRAY)
	# apply Gaussian blur to reduce noise
	gray_frame = cv2.GaussianBlur(gray_frame, (5, 5), 0)

	# identifies darkest parts of frame in order to mask out chromatophores from background
	#threshold = cv2.adaptiveThreshold(gray_frame, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV, 5, 2)  # B. good, semi noisy/inaccurate blur (11, 11) !
	#threshold = cv2.adaptiveThreshold(gray_frame, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV, 5, 3)  # C. some noise but true to size, blur (3, 3) !!
	#threshold = cv2.adaptiveThreshold(gray_frame,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,cv2.THRESH_BINARY_INV, 3 ,1)  # D. higher threshhold, noisy but true to size, blur (13, 13)
	#threshold = cv2.adaptiveThreshold(gray_frame, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV, 3, 2)  # E. not super solid but good mask, blur (5, 5)
	threshold = cv2.adaptiveThreshold(gray_frame, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV, CHROM_PIXEL_DIAMETER, 10)
	# 2nd to last param should be roughly the num of pixels of the diameter of an avg chromatophore in the frame

	return gray_frame, threshold  # return gray/blurred and masked frame


# define and draw the contour lines based on the masked frame
# takes in the binary masked image, returns the list of chromatophores
def contour(threshold, origImg):
	# detect the contours on the binary image
	contours, hierarchy = cv2.findContours(image=threshold, mode=cv2.RETR_EXTERNAL, method=cv2.CHAIN_APPROX_NONE)

	# draw contours on the original image
	contourCopy = origImg.copy()
	cv2.drawContours(image=contourCopy, contours=contours, contourIdx=-1, color=(0, 255, 0),
					 thickness=cv2.FILLED, lineType=cv2.LINE_AA)

	return contourCopy, contours


# create rotated bounding rectangles for each contour in contours
# returns a list of bounding rects
def createRotatedBRects(contours):
	return [cv2.minAreaRect(c) for c in contours]


# draw rotated rectangle
def drawRotatedRect(rect, img):
	# calculate corners of box
	box = cv2.boxPoints(rect)
	box = np.int0(box)

	# draw contour of box on image
	cv2.drawContours(img, [box], 0, (0, 0, 255), 2)

	return img


# sort the given contours from left to right
# return the sorted contours and sorted bounding boxes
def sortLR(contours):
	# construct the list of bounding boxes
	boundingBoxes = createRotatedBRects(contours)

	# sort from left to right
	(sortedContours, sortedBoundingBoxes) = zip(*sorted(zip(contours, boundingBoxes),
		key=lambda b:b[1][0][0]))

	# return the list of sorted contours and bounding boxes
	return sortedContours, sortedBoundingBoxes


# sort the given contours from top to bottom
# return the sorted contours and bounding boxes
def sortTB(contours):
	# construct the list of bounding boxes
	boundingBoxes = createRotatedBRects(contours)

	# sort from top to bottom
	(sortedContours, sortedBoundingBoxes) = zip(*sorted(zip(contours, boundingBoxes),
		key=lambda b:b[1][0][1]))

	# return the list of sorted contours and list of sorted bounding boxes
	return sortedContours, sortedBoundingBoxes


# takes in the chromatophore contours found in the first frame (ff)
# of the video; builds and returns the list of chromatophores
# (setting up to later add area data for each chromatophore)
def createAreasDict(ff_contours):
	global chromAreas

	# create "new chromatophore" in chromAreas for every chromatophore in 1st frame
	for i in range(0, len(ff_contours)):
		# create new chromatophore in chromAreas
		chromAreas[i] = {}
		i += 1


# populate the matchedCentroids dictionary with (objectID: centroid)
# in the first frame of the video when there is
# no previous frame to compare to
def createMatchedCentroids(curCentroids):
	for centroid in curCentroids:
		register(centroid)


# put the centroids of each bounding box in a
# list with the center of each bounding box in the
# same index as the original bounding box
# takes in list of bounding boxes and contours
# note: index of bBoxes and contours must correspond
# returns centroids and a dict that maps each centroid to its contour
def getCentroids(bBoxes, contours):
	centroids = []
	centroidsToContours = {}

	for curIndex in range(0, len(bBoxes)):
		centroid = bBoxes[curIndex][0]
		centroids.append(centroid)
		centroidsToContours[centroid] = contours[curIndex]

	return centroids, centroidsToContours


# register a new chromatophore in
def register(centroid):
	global matchedCentroids
	global disappeared
	global nextObjectID
	global chromAreas

	# when registering an object we use the next available object
	# ID to store the centroid
	matchedCentroids[nextObjectID] = centroid
	disappeared[nextObjectID] = 0
	chromAreas[nextObjectID] = {}
	initialCentroids[nextObjectID] = centroid
	nextObjectID += 1


# deregister given object in matchedCentroids and disappeared dicts
def deregister(objectID):
	global matchedCentroids
	global disappeared

	# to deregister an object ID we delete the object ID from
	# our respective dictionaries
	del matchedCentroids[objectID]
	del disappeared[objectID]
	del chromAreas[objectID]


# compare current centroids with previous centroids to
# determine the closest one
# and return the calculated IDs of each current centroid.
# takes in the list of current centroids and dict of previous matched centroids
# returns dict of ID_num:current centroids matched
# code based on https://pyimagesearch.com/2018/07/23/simple-object-tracking-with-opencv/
def matchCentroids(curCentroids):
	global matchedCentroids
	global nextObjectID
	global disappeared

	# grab the set of object IDs and corresponding centroids
	objectIDs = list(matchedCentroids.keys())
	objectCentroids = list(matchedCentroids.values())

	# compute the distance between each pair of object
	# centroids and input centroids, respectively -- our
	# goal will be to match an input centroid to an existing
	# object centroid
	D = dist.cdist(np.array(objectCentroids), curCentroids)

	# in order to perform this matching we must (1) find the
	# smallest value in each row and then (2) sort the row
	# indexes based on their minimum values so that the row
	# with the smallest value is at the *front* of the index
	# list
	rows = D.min(axis=1).argsort()

	# next, we perform a similar process on the columns by
	# finding the smallest value in each column and then
	# sorting using the previously computed row index list
	cols = D.argmin(axis=1)[rows]

	# in order to determine if we need to update, register,
	# or deregister an object we need to keep track of which
	# of the rows and column indexes we have already examined
	usedRows = set()
	usedCols = set()
	# loop over the combination of the (row, column) index
	# tuples
	for (row, col) in zip(rows, cols):
		# if we have already examined either the row or
		# column value before, ignore it
		# val
		if row in usedRows or col in usedCols:
			continue
		# otherwise, grab the object ID for the current row,
		# set its new centroid, and reset the disappeared
		# counter
		objectID = objectIDs[row]
		matchedCentroids[objectID] = curCentroids[col]
		disappeared[objectID] = 0
		# indicate that we have examined each of the row and
		# column indexes, respectively
		usedRows.add(row)
		usedCols.add(col)

	# compute both the row and column index we have NOT yet
	# examined
	unusedRows = set(range(0, D.shape[0])).difference(usedRows)
	unusedCols = set(range(0, D.shape[1])).difference(usedCols)

	# loop over the unused row indexes
	for row in unusedRows:
		# grab the object ID for the corresponding row
		# index and increment the disappeared counter
		objectID = objectIDs[row]
		disappeared[objectID] += 1
		# check to see if the number of consecutive
		# frames the object has been marked "disappeared"
		# for warrants deregistering the object
		if disappeared[objectID] > MAX_DISAPPEARED:
			deregister(objectID)

	for col in unusedCols:
		register(curCentroids[col])


# takes in chromatophore contours, calculates area for each chromatophore
# adds area to appropriate area time series
# (the list of areas under the chrom w/ the correct ID num)
# returns updated chromAreas dict
def calcCurAreas(centroidsToContours, curFrameIndex):
	global chromAreas
	global matchedCentroids
	global disappeared
	global PIXELS_PER_MM

	# iterate through each centroid that's paired with an ID number
	# and add the area of each centroid's contour into chromAreas for this frame index
	for ID_num, centroid in matchedCentroids.items():
		# if the centroid is in the current frame, add its area to chromAreas
		# but if the centroid disappeared
		# then don't add a new area entry into chromAreas, just skip adding
		# the current frame as a data point for that chromatophore
		if disappeared[ID_num] < 1:
			curContour = centroidsToContours[centroid]  # get corresponding contour
			contourArea = cv2.contourArea(curContour)  # calculate area of corresponding contour

			# convert area from pixels^2 to mm^2
			contourArea = contourArea/PIXELS_PER_MM

			# add entry of area for this chrom's area in this frame
			chromAreas[ID_num][curFrameIndex] = contourArea


# draw ID numbers, bounding boxes, and centroids for each object in the frame
# returns frame with everything drawn in it
def drawIDNums(boundingRects, frame):
	global matchedCentroids

	# draw bounding boxes
	for box in boundingRects:
		frame = drawRotatedRect(box, frame)

	# draw centroids and ID nums
	for ID, centroid in matchedCentroids.items():
		# draw both the ID of the object and the centroid of the
		# object on the output frame
		text = str(ID)
		cv2.putText(frame, text, (int(centroid[0]) - 5, int(centroid[1]) - 5),
					cv2.FONT_HERSHEY_DUPLEX, 0.3, (0, 255, 0), 1)
		cv2.circle(frame, (int(centroid[0]), int(centroid[1])), 2, (0, 255, 0), -1)

	return frame


# filter out irrelevant data entries (ex: chrom that was only identified in
# a few frames)
def cleanUpData():
	global chromAreas

	# make a copy to iterate through to avoid "mutated during iteration" error
	origChromAreas = chromAreas.copy()

	for id_num, areas in origChromAreas.items():
		# if a certain chromatophore has less than 25 area data entries
		# (~1 sec since usually about 25fps), delete it from the dataset
		if len(areas) < 25:
			del chromAreas[id_num]
			del initialCentroids[id_num]


# format final dataset into a spreadsheet
# important note: the xlwt library currently only supports up to 256 columns
# so the amount of chromatophores in the selected ROI must be < 257
def formatData():
	global chromAreas

	# Workbook is created
	wb = Workbook()
	areasSheet = wb.add_sheet('Areas')
	centroidsSheet = wb.add_sheet('Centroids')

	# iterate through each identified chromatophore
	for id_num, areas in chromAreas.items():
		# add centroid to the centroids sheet
		(x, y) = initialCentroids[id_num]
		centroidsSheet.write(1, id_num, x)
		centroidsSheet.write(2, id_num, y)

		# add area time series to areas sheet
		for frame_num, area in areas.items():
			areasSheet.write(frame_num, id_num, area)


	directory = curVidPath[:-4] + ".xls"
	print(curVidPath)
	wb.save(directory)


# plots chrom areas over time (each chromatophore is
# one line) then shows and saves the graph to the computer
def plotChromAreas():
	global chromAreas

	for id_num, areas in chromAreas.items():
		x = areas.keys()
		y = areas.values()

		plt.plot(x, y, label=str(id_num))

	# label the graph
	plt.xlabel('Frame Index')  # naming the x-axis
	plt.ylabel('Area (in mm)')  # naming the y-axis
	plt.title('Areas of individual chromatophores for ' + curVidPath)  # title the graph
	plt.legend()  # show a legend on the plot

	# show the graph
	plt.show()

	# save the graph
	#curVid = curVidPath.replace("/", "_")
	#curVid = curVid.replace(".", "")
	#plt.savefig('areas_' + curVid + '.png')


# show images from each step of the process
def showImages(curFrameIndex, original, gray, threshold, contours, ID_labeled):
	cv2.imshow("Original" + str(curFrameIndex), original)
	cv2.imshow("Grayscale + Gaussian Blur" + str(curFrameIndex), gray)
	cv2.imshow("Threshold" + str(curFrameIndex), threshold)
	cv2.imshow("Contours" + str(curFrameIndex), contours)
	cv2.imshow("IDs + centroids + bounding boxes " + str(curFrameIndex), ID_labeled)

	cv2.waitKey(1)

	cv2.destroyWindow("Original" + str(curFrameIndex))
	cv2.destroyWindow("Grayscale + Gaussian Blur" + str(curFrameIndex))
	cv2.destroyWindow("Threshold" + str(curFrameIndex))
	cv2.destroyWindow("Contours" + str(curFrameIndex))
	cv2.destroyWindow("IDs + centroids + bounding boxes " + str(curFrameIndex))



# saves each generated image individually in the frame_cap folder
# showing each step of the process for the current frame
# to use, call at the end of the while True loop in processData()
def saveIndividualImages(curFrameIndex, original, gray, threshold, contours, ID_labeled):
	cv2.imwrite("frame_cap/original_" + str(curFrameIndex) + ".png", original)
	cv2.imwrite("frame_cap/gray_" + str(curFrameIndex) + ".png", gray)
	cv2.imwrite("frame_cap/thresh_" + str(curFrameIndex) + ".png", threshold)
	cv2.imwrite("frame_cap/contour_" + str(curFrameIndex) + ".png", contours)
	cv2.imwrite("frame_cap/ID_labeled_" + str(curFrameIndex) + ".png", ID_labeled)


# saves one image that shows each step of the process for the current frame
# saves to the frame_cap folder
# to use, call at the end of the while True loop in processData()
def saveImages(curFrameIndex, original, gray, threshold, contours, ID_labeled):
	# convert 2 channel images (grayscale) into 3 channel images (RGB)
	gray = cv2.cvtColor(gray, cv2.COLOR_GRAY2RGB)
	threshold = cv2.cvtColor(threshold, cv2.COLOR_GRAY2RGB)

	# stack images vertically
	processImg = np.concatenate((original, gray), axis=0)
	processImg = np.concatenate((processImg, threshold), axis=0)
	processImg = np.concatenate((processImg, contours), axis=0)
	processImg = np.concatenate((processImg, ID_labeled), axis=0)

	# save in frame_cap folder
	cv2.imwrite("frame_cap/process_" + str(curFrameIndex) + ".png", processImg)


# summary function of entire program
# takes in filepath to raw video and returns time series
# of each chromatophores' change in area
def processData(vidPath):
	global matchedCentroids
	global chromAreas

	# open the video
	vid = cv2.VideoCapture(vidPath)
	curFrameIndex = 0

	# run computer vision process for each frame in the video
	while True:
		# read current frame from the video
		ok, frame = vid.read()

		# stop loop if error reading frame or end of vid reached
		if not ok:
			break

		# resize frame to make it easier to compare different windows
		frame = imutils.resize(frame, width=WINDOW_WIDTH)

		# select region of interest
		if curFrameIndex == 0:
			# Get ROI, record its (x,y) and height/width
			# note: x and y coordinates are of top left corner
			[ROI_x, ROI_y, ROI_width, ROI_height] = cv2.selectROI("Select ROI", frame, fromCenter=False, showCrosshair=True)

		# resize the frame to only include the ROI
		frame = frame[ROI_y:ROI_y + ROI_height, ROI_x: ROI_x + ROI_width]

		# mask out the chromatophores and display the masked frame
		gray_frame, thresh_frame = maskChromatophores(frame)

		# find contours from threshold
		contour_frame, contours = contour(thresh_frame, frame)

		# sort contours from left to right
		sortedContours, sortedBoundingBoxes = sortLR(contours)

		# put current centroids of each bounding box (unpaired with ID num) in a list
		curCentroids, centroidsToContours = getCentroids(sortedBoundingBoxes, sortedContours)

		# if first frame, create dict of chromatophores and begin object tracking
		if curFrameIndex == 0:
			createMatchedCentroids(curCentroids)  # initialize matchedCentroids dict
		# otherwise match each current centroid with an ID num
		else:
			matchCentroids(curCentroids)

		# draw matched ID numbers and centroids onto original frame
		ID_frame = drawIDNums(sortedBoundingBoxes, frame.copy())

		# show images
		#showImages(curFrameIndex, frame, gray_frame, thresh_frame, contour_frame, ID_frame)
		cv2.imshow("ID frame for " + str(curFrameIndex), ID_frame)
		cv2.waitKey(1)
		cv2.destroyWindow("ID frame for " + str(curFrameIndex))

		# save generated images to the frame_cap folder
		saveImages(curFrameIndex, frame, gray_frame, thresh_frame, contour_frame, ID_frame)

		# add areas of chromatophores from current frame to chromAreas
		calcCurAreas(centroidsToContours, curFrameIndex)

		curFrameIndex += 1  # update frame index

	# filter out irrelevant data entries
	cleanUpData()

	# format final dataset and save as an .xls file in video's original directory
	formatData()

	# plot chromatophore areas
	#plotChromAreas()

	# close graph display
	#plt.close()

	#print(chromAreas)
	# close any open windows
	cv2.destroyAllWindows()


# use argument vals if any given, otherwise use default nums
setGlobalNums()

# run the program's summary function
processData(curVidPath)
