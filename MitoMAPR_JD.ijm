//MitoMAPR.ijm from https://elifesciences.org/articles/49158/figures#scode1
//Slight modification by James Crichton for Jacob Day (University of Exeter) - 15/8/24

//Code segments and skeltonizes mitochondrial staining, then quantifies features of these skeletons
//Z-stacks are first projected (max intensity)

run("Fresh Start");//reset basic settings 

var SampleIdARRAY = newArray();
var ObjectsARRAY = newArray();
var JunctionsArray = newArray();
var JunctionsPerNetworkArray = newArray();
var NetworksARRAY = newArray();
var ObjectLengthARRAY = newArray()
var mitoAreaARRAY = newArray();
var MitoPartAreaArray = newArray();
var CoverageArray = newArray();
var CellAreaArray = newArray();


         
// MitoMAPR-------------------------------------------------------------------

    macro "MitoMAPR-1.2" {     

		mito_channel=3// channel with the mitochondial staining 

		//Open image
		img_path=File.openDialog("Open an image to analyse");
		open(img_path);
		imgName=getTitle();
		
		
		//Make directories for output
		//Images are .lif projects, so will make a directory for the project, and sub-directories for the images within that
		proj_dir=replace(img_path, ".lif.*$","");
		File.makeDirectory(proj_dir);
		img_dir=replace(imgName, "^.*lif - ","");
		dir3=proj_dir+"/"+img_dir;
		
		File.makeDirectory(dir3);
		
		//Make an 8-bit single image of the mitochondria for analysis
		
			Stack.getDimensions(width, height, channels, slices, frames);//Get dimensions of image opened
			
			//If z-stack then max intensity project
			if (slices>1){
				selectImage(imgName);
				run("Z Project...", "projection=[Max Intensity]");
				close(imgName);								
				projected_title=getTitle();
				close("original_img");
				selectWindow(projected_title);rename(imgName);//give the projection the name of the original image for furture continuity	
			}
			
		
		Stack.setDisplayMode("composite");
		setSlice(1);run("Grays");
		setSlice(2);run("Green");
		setSlice(3);run("Magenta");
		
		//Select region(s) of interest to analyse
		while (roiManager("count")==0){ //If this is -1, there is no selection, so repeat 
			setTool("rectangle");
			run("ROI Manager...");	
			roiManager("Show All");
			waitForUser ("Select regions of interest.\n  After each, press \"t\". \n When you're done press OK");
		}
		
		
		//Loop through each roi
		for (j = 0; j < roiManager("count"); j++) {
			
			roi_selection=j+1;
			print("processing roi "+roi_selection); 
			selectWindow(imgName);
			roiManager("select", j);
			roiManager("rename", "ROI_"+j);
			run("Duplicate...", "duplicate channels="+mito_channel);//copy the mictochnodrial channel at the ROI selected
//			run("Clear Outside");
			run("Select None");
			run("8-bit");
			
			dir4=dir3+"/ROI"+roi_selection;
			print(dir4);
			File.makeDirectory(dir4);
			
			rename("ROI.tif");
			selectWindow("ROI.tif");
			save(dir4+"/ROI.tif");
				
			    
			run("Duplicate...", "title=TheOne");
				selectWindow("TheOne");
				run("Grays");
				getDimensions(width, height, channels, slices, frames);
				getPixelSize(unit, pixelWidth, pixelHeight);
			
			
			//Image processing steps to try to 
			run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None*");//Modifies contrast levels
			run("Unsharp Mask...", "radius=2 mask=0.60");//Sharpens edges
			run("Median...", "radius=2");//
			
			//Create a binary mask
				run("Duplicate...", "title=Skeletor");
				selectWindow("Skeletor");
				run("32-bit");
				run("Make Binary");
				run("8-bit");
			
				saveAs(dir4+"/Mask.tif");//save this mask
				rename("Skeletor");
//			
			//Calculate roi measurements
			getStatistics(area, Mean, min, max); // These stats refer to the entire image frame, not a selection
				mitoArea = pow(pixelWidth, 2.0) * parseFloat(width) * parseFloat(height) * (Mean / parseFloat(max)) ;
				run("Restore Selection");
				getStatistics(area);
				TotalCellArea=area;
			    MitoCoverage = ((mitoArea/TotalCellArea)*100);
				run("Skeletonize");
				run("Red");
				selectWindow("TheOne");
				run("Add Image...", "image=Skeletor x=0 y=0 opacity=100 zero");
				size = toString(round(0.25*parseFloat(width)*parseFloat(pixelWidth)));
				run("Scale Bar...", "width=" + size + " height=4 font=14 color=White background=Black location=[Lower Right] bold overlay");
			
	
	
			
			selectWindow("Skeletor");
			run("Analyze Skeleton (2D/3D)", "prune=none show display");
	
			selectWindow("Tagged skeleton");
	        saveAs(dir4+"/Tagged skeletor.tif");
	        selectWindow("Skeletor-labeled-skeletons");
	        saveAs(dir4+"/labeled skeletor.tif");
			close("Tagged skeleton");
			close("Skeletor");
			close("Skeletor-labeled-skeletons");
	
			
	selectWindow("TheOne");
	saveAs(dir4+"/TheOne.tif");
	close();
	
			selectWindow("Results"); rows = nResults;
			ObjectCounts = newArray(rows);
			for (i=0; i<rows; i++) {
				ObjectCounts[i] = getResult("# Branches", i);
			}
			JunctionCounts = newArray(rows);
			for (i=0; i<rows; i++){
				JunctionCounts[i] = getResult("# Junctions", i);
			}
			selectWindow("Results"); run("Close");
			IJ.renameResults("Branch information", "Results");
			selectWindow("Results"); rows = nResults;
			ObjectLengths = newArray(rows);
			for (i=0; i<rows; i++) {
				ObjectLengths[i] = getResult("Branch length", i);
			}
	
			run("Close");
	
			
			
			Objects = (parseFloat(CountForMe(ObjectCounts)))*10;
			Networks = (parseFloat(countNetworksForMe(ObjectCounts)))*10;
			JunctionPoints = (AddForMe(JunctionCounts))*10;
			JunctionsPerNetwork = round(JunctionPoints/Networks);
	        MitoPartArea = mitoArea/Objects;
			ObjectLength = BeMeanForMe(ObjectLengths);
	
	
	Maths = newArray(  "Objects(OC)",
					   "ObjectLength(OL)",
					   "Networks(N)",
					   "JunctionPoints(JP)",
				       "JunctionsPerNetwork(JP/N)",
				       "MitochondrialFootprint(MF",
				       "ObjectArea(OA)",
				       "MitoCoverage(MC)",
				       "TotalCellArea(TCA)");
				     
	
	
	Numbers = newArray(Objects,
	                   ObjectLength,
					   Networks,
					   JunctionPoints,
					   JunctionsPerNetwork,
					   mitoArea,
					   MitoPartArea,
					   MitoCoverage,
					   TotalCellArea);
	
	Units = newArray("Counts",
	                  unit,
	                 "Counts",
	                 "Counts",
				     "Counts",
				      unit+" squared",
				      unit+" squared",
				      "Percent",
				      unit+" squared");
	
	
	Array.show("Data", Maths, Numbers, Units);
	
	selectWindow("Data");
	saveAs("results",dir4+"/Data.csv");
	
	close("Data.csv");       
	close("ROI.tif");

	 }
//	 	run("Close All");
	roiManager("save", dir3+"/Rois.zip");
	}
//}


//Do Maths for Me
function CountForMe(data) {
	entries = data.length;
	total = 0.0;
	for (i=0; i<entries; i++) 
		{
			total = total + 1;
		}
		return(total);
}

function BeMeanForMe(data) {
	entries = data.length;
	total = 0.0;
	for (i=0; i<entries; i++) {
		total = total + data[i];
	}
	ave = total/parseFloat(entries);
	return(ave);
}

function AddForMe(data) {
	entries = data.length;
	total = 0.0;
	for (i=0; i<entries; i++) {
		total = total + data[i];
	}
	sum = total;
	return(sum);
}

function countNetworksForMe(data) {
	entries = data.length;
	total = 0.0;
	for (i=0; i<entries; i++) {
		if (data[i] > 1) {
			total = total + 1;
		}
		else {
		}
	}
	return(total);
}

function CountNetworkBranchesForMe(data) {
	entries = data.length;
	total = 0.0;
	for (i=0; i<entries; i++) {
		if (data[i] > 1) {
			total = total + data[i];
		}
		else {
		}
	}
	return(total);
}

