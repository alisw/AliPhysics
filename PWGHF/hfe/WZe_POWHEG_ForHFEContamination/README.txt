
//======================================================================
//======================================================================
//======================================================================
ONLINE PART -> To generate .root files in alice grid
In Folder AliceGridAnalysis :

Make Changes in :

JDL                   : Change the path of alimoniter folder
base_powheg.input     : Change parameters such as PDG, Decay modes, PDF , atomic number etc. as required.
GenW_Pythia6_POWHEG.C : Change parameters such as y range for accepting electrons etc., as required.


//======================================================================
For running in ALICE grid :

1. Create the folder that is mentioned in the JDL file in the alimoniter page
2. Upload all contents of AliceGridAnalysis folder there
3. Click on JDL and "submit"
4. After the .root files are generated in as output, the offline analysis has to be done.


//======================================================================
For changing from W to Z make changes in files :

base_powheg.input     : The PDG to 23 from 24

In file simrun.sh : 
Line 139 use : ./runPowheg.sh -v Z $(expr $nevts \* 20) $CONFIG_SEED > powheg.log 2>&1

//======================================================================
//======================================================================
//======================================================================
OFFLINE PART  :
In Folder OfflineAnalysis

1. Analyze each output by : readPrimaries.C 
       Line 56 : PDG to be chnaged for Z (23) and W(24)
       Line 58 : Y range
2. Then output has to be merged separately for +ve and -ve
3. Plot_We.C : Normalisation  done               
       In the factor norm, the last term is N_tried/ N_sigma (mb).
       The values of N_tried and N_sigma are found at the end of the sim.log file in the output folder
4. Positive and negative are combined with the code Comb_We.C (The combination of Z can be included in the code if required)


