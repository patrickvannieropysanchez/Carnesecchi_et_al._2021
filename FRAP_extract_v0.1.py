from ij import IJ, ImagePlus
from ij.gui import ShapeRoi, Roi
from ij.io import Opener
from ij.measure import ResultsTable
from ij.plugin import Commands, WindowOrganizer, ZAxisProfiler, ZProjector
from ij.plugin.filter import ThresholdToSelection
from ij.plugin.frame import ThresholdAdjuster, RoiManager, ContrastAdjuster

import os, re, sys, time

show_steps = True
Commands().closeAll()

###################################
########### Definitions ###########
###################################
def Align(imp, roi, method = 5):

	imp2 = imp.duplicate()
	windowsizex = roi.getFloatHeight()
	windowsizey = roi.getFloatWidth()
	x0 = roi.getXBase()
	y0 = roi.getYBase()

	IJ.run(imp2, "Align slices in stack...", "method=5 windowsizex="+str(int(roi.getFloatHeight()))+
	" windowsizey="+str(int(roi.getFloatWidth()))+
	" x0="+str(int(roi.getXBase()))+
	" y0="+str(int(roi.getYBase()))+
	" swindow=10 subpixel=true itpmethod=1 ref.slice=1 show=true")
	
	IJ.run(imp2, "Select None", "")
	
	return imp2

def BleachT(imp):
	intensity_t = list(ZAxisProfiler().getPlot(imp, "time").getYValues())
	intensity_diff = [j-i for i, j in zip(intensity_t[:-1], intensity_t[1:])]
	t_min = intensity_diff.index(min(intensity_diff))+2 # python indexing + intensity_diff shift
	
	return t_min

def LargestRoi(roi):

	rois = ShapeRoi(roi).getRois()

	if len(rois) != 1:
		roi_areas = []
		for roi in rois:
			roi_areas.append(roi.getStatistics().area)
			
		index_max = roi_areas.index(max(roi_areas))
		roi = rois[index_max]

	return ShapeRoi(roi)
###################################
###################################
###################################

workingDir = "/Users/patrick/Documents/PhD Local/Microscopy/FRAP2020/"
for image in os.listdir(workingDir):
	start = time.time()
	rm = RoiManager().getInstance()
	rm.reset()
	table = ResultsTable()
	table.reset()
	
	if image.endswith(".nd2"):
		print "Processing", image
		imp = Opener().openUsingBioFormats(os.path.join(workingDir, image))
		imp.getProcessor().setAutoThreshold("Li", 1, 1)
		roi_ref = LargestRoi(ShapeRoi(ThresholdToSelection().convert(imp.getProcessor())))
		imp.setRoi(roi_ref)

		if show_steps == 1:
			imp.show()
		
		# Align stack
		imp_aligned = Align(imp, roi_ref)

		if show_steps == 1:
			imp_aligned.setTitle("Aligned_stack")
			imp_aligned.show()

		# Find bleaching frame index
		t_bleach = BleachT(imp)

		# Get pre-bleach ROI of whole nucleus based on median pre-bleach
		imp_median = ZProjector().run(imp_aligned.crop("1-"+str(t_bleach-1)), "median")

		imp_median.getProcessor().setAutoThreshold("Li dark b&w")

		roi_nuc = LargestRoi(ShapeRoi(ThresholdToSelection().convert(imp_median.getProcessor())))

		if show_steps == 1:
			imp_median.setTitle("Median")
			imp_median.show()

		## Import imp_stim and convert to roi_stim
		imp_stim = Opener().openUsingBioFormats(os.path.join(workingDir, image[0:4]+"_ROI.tif"))
		imp_stim.getProcessor().setThreshold(1, 2, 1)
		roi_stim = ShapeRoi(ThresholdToSelection().convert(imp_stim.getProcessor()))

		if show_steps == 1:
			imp_stim.setTitle("ROI Stimulation")
			imp_stim.show()
		
		# Define non-bleached and bleached ROIs
		roi_nuc_bleach = ShapeRoi(roi_nuc).and(ShapeRoi(roi_stim))

		roi_nuc_nonbleach = ShapeRoi(roi_nuc).not(ShapeRoi(roi_stim))

		imp_max = ZProjector().run(imp, "max")

		if show_steps == 1:
			imp_max.setTitle("Maximum")
			imp_max.show()
			
		imp_max.getProcessor().setAutoThreshold("Li light b&w")
		roi_bg = ThresholdToSelection().convert(imp_max.getProcessor())

		if show_steps == 1:
			rm = RoiManager().getInstance()
			rm.addRoi(roi_nuc)
			rm.addRoi(roi_stim)
			rm.addRoi(roi_nuc_nonbleach)
			rm.addRoi(roi_nuc_bleach)
			rm.addRoi(roi_bg)
		
		#### Main loop ####
		#### Measure intensities through frames of roi_nuc, roi_nuc_nonbleach, roi_nuc_bleach, roi_bg
		#### Main loop ####
		roi_list = {"Nucleus" : roi_nuc, "NonBleach" : roi_nuc_nonbleach, "Bleach" : roi_nuc_bleach, "BG": roi_bg}
		
		timestamps = open("/Users/patrick/Documents/PhD Local/Microscopy/FRAP2020/timestamps/timestamps.txt")
		
		time_list = []
		
		for t in timestamps:
			timeSplit = [int(s) for s in re.findall(r'\b\d+\b', t)]
			sec = timeSplit[0]*60 + timeSplit[1] + int(round(timeSplit[2]/10.0))/1000.0
			time_list.append(sec)
		
		for frame in range(1, imp.getNFrames()+1):
		
			timestamp = time_list[frame-1]
			table.setValue("Time", frame-1, timestamp)
			
			imp_aligned.setPosition(1, 1, frame)

			for roi in roi_list:
				imp_aligned.setRoi(roi_list[roi])
				table.setValue(roi+"Mean", frame-1, imp_aligned.getRoi().getStatistics().mean)
				table.setValue(roi+"Area", frame-1, imp_aligned.getRoi().getStatistics().area)

		if show_steps == 1:
			table.show("FRAP Results")
			WindowOrganizer().run("tile") # Tile all the images on the screen for visibility

		table.save("/Users/patrick/Documents/PhD Local/Microscopy/FRAP2020/Data_new/"+imp.getTitle()[0:4]+".txt")

		print "Task took ", time.time() - start, " sec to complete."
		
print("Finished")