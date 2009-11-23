Calibration and Performace  train at GSI:

1. The calibration train will run on the regular basis.
2. The results of the calibration and the setup of the calibration 
   will be done in workspace.
3. Content of the workspace:
   3.a) ConfigOCDB.C macro - define the OCDB setup (Default and Specific 
        storage) + 
        the RecoParam used for refitting of the track (possibility to switch On/OFF different kind of corrections)
   3.b) CalibrateTPC.C macro - define the Tasks and components of the current calibration / performance train

4. run.list - list of runs used in the calibration
5. esdxxx.txt files - list of the esd for run xxx
6. Run directories 
   The calibration at GSI is done on the batch farm. The results are merged per run.



List of files:
Macros:
CalibrateTPC.C - standard in $ALICE_ROOT/TPC/macros   
ConfigOCDB.C   - example (to be modified according studies)  in $ALICE_ROOT/TPC/macros 
filterESD.C    - $ALICE_ROOT/TPC/macros/filterESD.C

Shell scripts:
balice.sh      - symbolic link to the AliRoot setup script
               - example (to be modified):
	         ln -sf ~/.balice64HEAD0108 balice.sh
alienSetup.sh  - symbolic link to the alien setup script
               - hint (don't modify unles you know what are you doing): 
	       - ln -sf $ALICE_ROOT/TPC/CalibMacros/alienSetupGSI.sh alienSetup.sh
submitCalib.sh - hint: copy of the (don't change it unless ...)
               - cp $ALICE_ROOT/TPC/CalibMacros/submitcalib.sh .

Run lists:
(No recepie for the moment how to create it)
run.list      - main source of the run numbers  (created from the logbook)      
esd.list      - the list of files -standard 
    

Run the train:
runTrainGSI.sh  - this is a recepie
                - $ALICE_ROOT/TPC/scripts/runTrainBatch.sh    