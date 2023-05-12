//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				MOE Golgi Analyzer by Vincent Rigolot 																								//
//					version 1.0 28/06/2022 																											//
//	 																																				//
// 		for analysis of confocal fluorescence microscopy images of mammalian cells 																	//
//		requires ImageJ v1.4 with Bio-render plugin, 																								//
//		images should be as .nd2 format but it can easily be changed, simply search & replace all occurences of	".nd2" with your format				//
//																																					//
//		images should be organized with every replicate of a same test-condition in a unique folder													//
//		the macro will analyze the whole folder at once and will create a folder in it to save results												//
//																																					//
//		adapted to 3-channel confocal fluorescence images (200*200µm), with :																		//
//		 	C1 = nucleus labeling (e.g. DAPI, Hoechst, etc.)																						//	
//			C2 = signal of interest, the one you want to measure in whole cells & in the region of interest											//
//			C3 = region of interest (ROI) (e.g. an antibody directed against a particular organelle)												//
//																																					//
//		this macro will :						 																									//
//	 	count the cells according to C1 (user input of threshold values for C1) 																	//
//	 	create ROI(s) according to C3 (user input of threshold values, or manual setting of each image for C3) 										//
//	 	measure signal of C2 (mean min max grey values, integrated density, area) in whole cells (user input of threshold values for C2)			//
//		measure signal of C2 in ROI(s) 																												//
//	 	save results as a .csv file																													//
//																																					//
//		it will also create several .png images for each analyzed one : 																			//
//				C1+nucleusROI (to assess correct cell counting)																						//
//				C3+ROIC3 (to assess correct creation of ROI(s) from C3 signal)																		//
//				C2 (glow LUT) + ROIC3 (to assess correct thresholding of C2 signal)																	//
//				C2+ROIC3 																															//
//				merge C1+C2+C3																														//
//				 																																	//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//%%%%%%%%%%%%%%% INITIAL SETTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
run("Set Measurements...", "area mean min integrated display redirect=None decimal=3");

// Dialog box creation for user input
  ID = "NoID";
  threshMinC1=500;
  threshMaxC1=2000;
  threshMinC2=90;
  threshMaxC2=1800;
  threshMinC3=0;
  threshMaxC3=4096;
  Dialog.create("MOE Golgi Analyzer by Vincent");
  Dialog.addString("Unique ID for saving results :", ID);
  Dialog.addMessage("Threshold setting");
  Dialog.addNumber("lower threshold C1", threshMinC1);
  Dialog.addToSameRow();
  Dialog.addNumber("upper threshold C1", threshMaxC1);
  Dialog.addNumber("lower threshold C2", threshMinC2);
  Dialog.addToSameRow();
  Dialog.addNumber("upper threshold C2", threshMaxC2);
  Dialog.addCheckbox("Manual setting C3 ?", true);
  Dialog.addNumber("otherwise lower threshold C3", threshMinC3);
  Dialog.addToSameRow();
  Dialog.addNumber("upper threshold C3", threshMaxC3);
  Dialog.addMessage(" ");
  Dialog.addCheckbox("Create and save images ? (slow !)", false);
  Dialog.addCheckbox("Debug log", false);
  Dialog.show();
  ID = Dialog.getString();
  threshMinC1=Dialog.getNumber();
  threshMaxC1=Dialog.getNumber();
  threshMinC2=Dialog.getNumber();
  threshMaxC2=Dialog.getNumber();
  manualthreshC3 = Dialog.getCheckbox();
  threshMinC3= Dialog.getNumber();
  threshMaxC3 = Dialog.getNumber();
  makeImg = Dialog.getCheckbox();
  dbug = Dialog.getCheckbox();

rep=getDirectory("Select a folder - "+ID);
File.makeDirectory(rep+"results_"+ID);
filelist=getFileList(rep);

if (dbug == true) {
print("There are "+filelist.length+" files in rep");			
}

// Compte le nombre de fichiers nd2 dans rep
nbre_nd2=0;
for (k=0; k<filelist.length; k++) {
	if (indexOf(filelist[k], ".nd2")>0) {
		nbre_nd2++;
	}
}

if (dbug == true) {
print("There are"+nbre_nd2+".nd2 files in rep");
}


// Creation of arrays which size equals the number of .nd2 files in rep (one Array for each parameters of interest)
t_nombre_cellules=newArray(nbre_nd2);
t_imagename=newArray(nbre_nd2);
t_nameNoExt=newArray(nbre_nd2);
t_nombre_Golgi=newArray(nbre_nd2);
t_threshold_low_Golgi=newArray(nbre_nd2);
t_threshold_up_Golgi=newArray(nbre_nd2);
t_threshMinC2=newArray(nbre_nd2);
t_threshMaxC2=newArray(nbre_nd2);
t_RawIntDen_C2_OutOfGolgi=newArray(nbre_nd2);
t_RawIntDen_C2_Golgi=newArray(nbre_nd2);
t_IntDen_C2_Golgi=newArray(nbre_nd2);
t_IntDen_C2_OutOfGolgi=newArray(nbre_nd2);
t_mean_C2_Golgi=newArray(nbre_nd2);
t_min_C2_Golgi=newArray(nbre_nd2);
t_max_C2_Golgi=newArray(nbre_nd2);
t_mean_C2_OutOfGolgi=newArray(nbre_nd2);
t_min_C2_OutOfGolgi=newArray(nbre_nd2);
t_max_C2_OutOfGolgi=newArray(nbre_nd2);


//%%%%%%%%%%%%%%%% MAIN LOOP FOR FILE TREATMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (manualthreshC3 == false) {
	setBatchMode(true);
}

for (i=0; i<filelist.length; i++) {	
	run("Clear Results");
	
	if (dbug == true) {
	print("\n\nAnalysis of file number "+i+" on "+filelist.length);
	}
	
	if (indexOf(filelist[i], ".nd2")>0) {

		run("Bio-Formats Windowless Importer", "open="+rep+filelist[i]);
		
		if (dbug == true) {
		print("\n\nAnalysis of image number "+i+" on "+nbre_nd2);
		print("Image pathway : "+rep+filelist[i]);
		}
				
		// storage of the original filename
		t_imagename[i] = filelist[i];
		t_nameNoExt[i] = File.nameWithoutExtension;
		rename("image");
		
		if(roiManager("count")>0) {
			roiManager("Delete"); 
		}
		if(roiManager("count")>0) {
			roiManager("Delete"); 
		}

//%%%%%%%%%%%%%%%% TREATMENT OF C1 SIGNAL (CELL COUNTING)	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		run("Split Channels");
		selectWindow("C1-image");
		setThreshold(threshMinC1, threshMaxC1);
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Analyze Particles...", "size=20-Infinity show=Masks display exclude clear include add");
		
		if (makeImg == true) {
		roiManager("Save", rep+"results_"+ID+File.separator+t_nameNoExt[i]+"_ROIM_nucleus.zip");
		}

		// stores the number of cells in array t_nombre_cellules
		t_nombre_cellules[i]=roiManager("count");
		if (dbug == true) {
		print("After analyzing particles in C1, number of ROIs in ROI Manager = "+roiManager("count"));
		}

		
			if(roiManager("count")>0) {
				roiManager("Delete"); 
			}
			if(roiManager("count")>0) {
				roiManager("Delete"); 
			}
	
//%%%%%%%%%%%%%%%% TREATMENT OF C3 (CREATION OF ROI(s)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		selectWindow("C3-image");
		run("Median...", "radius=1");
		setAutoThreshold("Default dark");
		setOption("BlackBackground", true);
		

		// (manual) thresholding of C3, saves the values used
		selectWindow("C3-image");
		if (manualthreshC3 == true) {
		run("Threshold...");
		title = "Thresholding C3";
		msg = "Set your threshold as desired and click OK (NOT APPLY !).";
		waitForUser(title, msg);
		getThreshold(t_threshold_low_Golgi[i], t_threshold_up_Golgi[i]);
			if (t_threshold_up_Golgi[i] > 4096) {
			t_threshold_up_Golgi[i] = 4096;
			setThreshold(t_threshold_low_Golgi[i], 4096);
			}
		}
		else {
			t_threshold_low_Golgi[i] = threshMinC3;
			t_threshold_up_Golgi[i] = threshMaxC3;
				if (t_threshold_up_Golgi[i] > 4096) {
				t_threshold_up_Golgi[i] = 4096;
				}
			setThreshold(t_threshold_low_Golgi[i], t_threshold_up_Golgi[i]);
		}
		
		// conversion to mask
		run("Convert to Mask");
		// if desired add your binary treatment of C3 mask here
		run("Analyze Particles...", "size=1-Infinity display add");
		
		if (makeImg == true) {
		roiManager("Save", rep+"results_"+ID+File.separator+t_nameNoExt[i]+"_ROIM_C3.zip");
		}

		// stores the number of C3 ROI(s) in array t_nombre_Golgi
		t_nombre_Golgi[i]=roiManager("count");
		if (dbug == true) {
		print("After analyzing particles in C3, number of ROI(s) = "+roiManager("count"));
		}
		
		run("Clear Results");

		
//%%%%%%%%%%%%%%%% TREATMENT OF C2 AND MEASURE OF SIGNAL IN WHOLE CELLS AND ROI(s) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		selectWindow("C2-image");
		run("32-bit");
		setThreshold(threshMinC2, threshMaxC2);
		run("NaN Background");
		t_threshMinC2[i]=threshMinC2;
		t_threshMaxC2[i]=threshMaxC2;

		// measure whole cell
		selectWindow("C2-image");
		run("Measure");
		t_RawIntDen_C2_OutOfGolgi[i]=getResult("RawIntDen", 0);
		if (dbug == true) {
		print("RawIntDen_C2_OutofGolgi= "+t_RawIntDen_C2_OutOfGolgi[i]);
		}

		selectWindow("C2-image");
		run("Measure");
		t_IntDen_C2_OutOfGolgi[i]=getResult("IntDen", 0);
		if (dbug == true) {
		print("IntDen_C2_OutofGolgi= "+t_IntDen_C2_OutOfGolgi[i]);
		}
		
		t_mean_C2_OutOfGolgi[i]=getResult("Mean", 0);
		if (dbug == true) {
		print("Mean_C2_OutofGolgi= "+t_mean_C2_OutOfGolgi[i]);
		}

		t_min_C2_OutOfGolgi[i]=getResult("Min", 0);
		if (dbug == true) {
		print("Min_C2_OutofGolgi= "+t_min_C2_OutOfGolgi[i]);
		}

		t_max_C2_OutOfGolgi[i]=getResult("Max", 0);
		if (dbug == true) {
		print("Max_C2_OutofGolgi= "+t_max_C2_OutOfGolgi[i]);
		}
		
		run("Clear Results");

		// measure of C2 in ROIC3
		roiManager("Measure");
		
		// somme des RawIntDen (1 ligne par ROI) pour stocker RawIntDen totale des Golgi dans tableau t_RawIntDen_C2_Golgi
		// calculate the parameters sum of every ROI 
		
		for (j = 0; j < nResults; j++) {
			RawIntDen_j=getResult("RawIntDen", j);
			if (isNaN(RawIntDen_j) == true) {
				RawIntDen_j=0;
			}
			t_RawIntDen_C2_Golgi[i]+=RawIntDen_j;
			if (dbug == true) {
			print("RawIntDen_C2_Golgi after "+j+"+1 loop turn(s) (on "+nResults+") ="+t_RawIntDen_C2_Golgi[i]);			
			}

			IntDen_j=getResult("IntDen", j);
			if (isNaN(IntDen_j) == true) {
				IntDen_j=0;
			}
			t_IntDen_C2_Golgi[i]+=IntDen_j;
			if (dbug == true) {
			print("IntDen_C2_Golgi after "+j+"+1 loop turn (on "+nResults+") ="+t_IntDen_C2_Golgi[i]);			
			}

			mean_j=getResult("Mean", j);
			if (isNaN(mean_j) == true) {
				mean_j=0;
			}
			t_mean_C2_Golgi[i]+=mean_j;
			if (dbug == true) {
			print("mean_C2_Golgi after "+j+"+1 loop turn (on "+nResults+") ="+t_mean_C2_Golgi[i]);			
			}

			min_j=getResult("Min", j);
			if (isNaN(min_j) == true) { 
				min_j=0;
			}
			t_min_C2_Golgi[i]+=min_j;
			if (dbug == true) {
			print("min_C2_Golgi after "+j+"+1 loop turn (on "+nResults+") ="+t_min_C2_Golgi[i]);			
			}

			max_j=getResult("Max", j);
			if (isNaN(max_j) == true) { 
				max_j=0;
			}
			t_max_C2_Golgi[i]+=max_j;
			if (dbug == true) {
			print("max_C2_Golgi after "+j+"+1 loop turn (sur "+nResults+") ="+t_max_C2_Golgi[i]);			
			}
		}

		// meaning of integrated density and min mean and max grey values
		t_IntDen_C2_Golgi[i]=t_IntDen_C2_Golgi[i]/(j+1);
		t_mean_C2_Golgi[i]=t_mean_C2_Golgi[i]/(j+1);
		t_min_C2_Golgi[i]=t_min_C2_Golgi[i]/(j+1);
		t_max_C2_Golgi[i]=t_max_C2_Golgi[i]/(j+1);

		run("Clear Results");
		run("Close All");
		
	}
}

//%%%%%%%%%%%%%%	RESULTS SAVING	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Displays the results stored in the different arrays in the Results window
// then saves the results as a .csv file
for (i=0; i<nbre_nd2; i++) {
	setResult("Image name", i, t_nameNoExt[i]);
	setResult("Number of cells", i, t_nombre_cellules[i]);	
	setResult("Number of ROI(s)", i, t_nombre_Golgi[i]);
	setResult("threshold min C3", i, t_threshold_low_Golgi[i]);
	setResult("threshold max C3", i, t_threshold_up_Golgi[i]);
	setResult("threshold min C2", i, t_threshMinC2[i]);
	setResult("threshold max C2", i, t_threshMaxC2[i]);
	setResult("RawIntDen C2 in C3", i, t_RawIntDen_C2_Golgi[i]);
	setResult("RawIntDen C2 whole cells", i, t_RawIntDen_C2_OutOfGolgi[i]);
	setResult("IntDen C2 in C3", i, t_IntDen_C2_Golgi[i]);
	setResult("IntDen C2 whole cells", i, t_IntDen_C2_OutOfGolgi[i]);
	setResult("Mean Strepta in Golgi", i, t_mean_C2_Golgi[i]);
	setResult("Min Strepta in Golgi", i, t_min_C2_Golgi[i]);
	setResult("Max Strepta in Golgi", i, t_max_C2_Golgi[i]);
	setResult("Mean Strepta Out of Golgi", i, t_mean_C2_OutOfGolgi[i]);
	setResult("Min Strepta Out of Golgi", i, t_min_C2_OutOfGolgi[i]);
	setResult("Max Strepta Out of Golgi", i, t_max_C2_OutOfGolgi[i]);
	updateResults();
	// saving "Results" in "rep/sauvegarde_resultats/resultats_"ID".csv"
	saveAs("Results", rep+"results_"+ID+File.separator+"results_"+ID+".csv");
}

//%%%%%%%%%%%%%%	IMAGES CREATION	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (makeImg == true) {	

		title_warn = "MOE Golgi Analyzer by Vincent - Images Creation";
		msg_warn = "Measures saved. The macro will now create images. \n It can take a while ! Do not do anything until prompted by the final confirmation message";
		waitForUser(title_warn, msg_warn);
		
	setBatchMode(true);
	
	for (i=0; i<filelist.length; i++) {	
		if (indexOf(filelist[i], ".nd2")>0) { 

	// Thresholding initial image
	
			run("Bio-Formats Windowless Importer", "open="+rep+filelist[i]);
			rename("image");
			setSlice(1);
			getMinAndMax(mindapi, maxdapi);
			setMinAndMax(threshMinC1, maxdapi);
			setSlice(2);
			setMinAndMax(t_threshMinC2[i], t_threshMaxC2[i]);
			setSlice(3);
			setMinAndMax(t_threshold_low_Golgi[i], t_threshold_up_Golgi[i]);

	// Creation of different montages
		
			// merge C1 C2 C3 thresholded
			selectWindow("image");
			run("Duplicate...", "title=Montage111 duplicate");
			Stack.setDisplayMode("composite");
			Stack.setActiveChannels("111");
			run("Make Composite");
			saveAs("PNG", rep+"results_"+ID+File.separator+t_nameNoExt[i]+"_C1-C2-C3");
			close();		

			// merge C2 C3 thresholded
			selectWindow("image");
			run("Duplicate...", "title=Montage111 duplicate");
			Stack.setDisplayMode("composite");
			Stack.setActiveChannels("011");
			run("Make Composite");
			saveAs("PNG", rep+"results_"+ID+File.separator+t_nameNoExt[i]+"_C2-C3");
			close();	
			
			// C2 thresholded displayed with glow LUT
			if(roiManager("count")>0) { 
			roiManager("Delete"); 
			}
			if(roiManager("count")>0) {
			roiManager("Delete"); 
			}
			
			selectWindow("image");
			setSlice(2);
			run("Duplicate...", "title=Strepta-Glow");
			run("glow");
			roiManager("Open", rep+"results_"+ID+File.separator+t_nameNoExt[i]+"_ROIM_C3.zip");
			roiManager("Show All without labels");
			saveAs("PNG", rep+"results_"+ID+File.separator+t_nameNoExt[i]+"_C2glow_ROIC3");
			close();	

			// C1 + ROInucleus
			if(roiManager("count")>0) { 
			roiManager("Delete"); 
			}
			if(roiManager("count")>0) {
			roiManager("Delete"); 
			}

			roiManager("Open", rep+"results_"+ID+File.separator+t_nameNoExt[i]+"_ROIM_nucleus.zip");
			selectWindow("image");
			setSlice(1);
			run("Duplicate...", "title=DAPI-ROI");
			roiManager("Show All with labels");
			run("Flatten");
			saveAs("PNG", rep+"results_"+ID+File.separator+t_nameNoExt[i]+"_C1_ROInucleus");
			close();	

			// C2 thresholded + C3 ROI(s)
			if(roiManager("count")>0) { 
			roiManager("Delete"); 
			}
			if(roiManager("count")>0) {
			roiManager("Delete"); 
			}
			
			roiManager("Open", rep+"results_"+ID+File.separator+t_nameNoExt[i]+"_ROIM_C3.zip");
			selectWindow("image");
			setSlice(2);
			run("Duplicate...", "title=Strepta-ROI-Golgi");
			roiManager("Show All without labels");
			run("Flatten");
			saveAs("PNG", rep+"results_"+ID+File.separator+t_nameNoExt[i]+"_C2_ROIC3");
			close();		

			// C3 thresholded + C3 ROI(s)
			if(roiManager("count")>0) {
			roiManager("Delete"); 
			}
			if(roiManager("count")>0) {
			roiManager("Delete"); 
			}
			
			roiManager("Open", rep+"results_"+ID+File.separator+t_nameNoExt[i]+"_ROIM_C3.zip");
			selectWindow("image");
			setSlice(3);
			run("Duplicate...", "title = TGN-ROI-Golgi");
			run("Median...", "radius=1");
			roiManager("Show All without labels");
			run("Flatten");
			saveAs("PNG", rep+"results_"+ID+File.separator+t_nameNoExt[i]+"C3_ROIC3");
			close();
			
		}
		run("Close All");
	}
}

		title_success = "MOE Golgi Analyzer by Vincent - Success !";
		msg_success = "Operation completed with success !";
		waitForUser(title_success, msg_success);