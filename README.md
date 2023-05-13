# MOE-Golgi-Analyzer

DOI : 10.5281/zenodo.6778042

Description of the folder content :


1) The macro in .ijm format.
Suited for analysis of 3-channel confocal fluorescence microscopy images of mammalian cells (~200*200µm).        
Requires ImageJ v1.4 with Bio-render plugin.
Images should be as .nd2 format but it can easily be changed, simply search & replace all occurences of ".nd2" with your format in the macro code.
Images should be organized with every replicate of a same test-condition in a unique folder. The macro will analyze the whole folder at once and will create a folder in it to save results.


2) A folder named "example_data", it contains 3 representative images that can be used to test the macro. 
It also contains a results folder with representative data obtained via the analysis of these representative images with the macro (see Description of the macro for description of the results obtained)

____________________________

Description of the macro :

input : 3-channel image with 

C1 = nucleus labeling (e.g. DAPI, Hoechst, etc.)                                                                                
C2 = signal of interest, the one you want to measure in whole cells & in the region of interest                                            
C3 = region of interest (ROI) (e.g. an antibody directed against a particular organelle, in our case Golgi apparatus)                                                
                                                                                                                                                    
this macro will :                                                                                                                             
count the cells according to C1 (user input of threshold values for C1)                                                                     
create ROI(s) according to C3 (user input of threshold values, or manual setting of each image for C3)                                         
measure signal of C2 (mean min max grey values, integrated density, area) in whole cells (user input of threshold values for C2)   measure signal of C2 in ROI(s)                                                                                                                 
save results as a .csv file                                                                                                                    
                                                                                                                                                  
it will also create several .png images for each analyzed one :                                                                             
C1+nucleusROI (to assess correct cell counting)                                                                                        
C3+ROIC3 (to assess correct creation of ROI(s) from C3 signal)                                                                        
C2 (glow LUT) + ROIC3 (to assess correct thresholding of C2 signal)                                                                    
C2+ROIC3                                                                                                                             
merge C1+C2+C3
