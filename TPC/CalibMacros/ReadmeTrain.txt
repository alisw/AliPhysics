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


   
     
    