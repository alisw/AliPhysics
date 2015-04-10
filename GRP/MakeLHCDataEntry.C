Float_t ProcessEnergy(TObjArray* const array, Double_t timeStart){

  //
  // Method to processo LHC Energy information
  // Only the first value is returned, provided that it is withing DAQ_time_start and DAQ_time_end
  //

  Int_t nCounts = array->GetEntries();
  Float_t energy = -1;
  Double_t timeEnergy = -1;
  Int_t indexEnergy = -1;
  Bool_t foundEnergy = kFALSE;

  Printf("Energy measurements = %d\n",nCounts);
  if (nCounts ==0){
    AliWarning("No Energy values found! Beam Energy remaining invalid!");
  }
  else{
    for (Int_t i = 0; i < nCounts; i++){
      AliDCSArray *dcs = (AliDCSArray*)array->At(i);
      if (dcs){
        if (dcs->GetTimeStamp()<=timeStart && dcs->GetTimeStamp()>=timeEnergy){// taking always the very last entry: of two measurements have the same timestamp, the last one is taken
          timeEnergy = dcs->GetTimeStamp();
          indexEnergy = i;
          foundEnergy = kTRUE;
        }
        else{
          break;
        }
      }
    }
    if (!foundEnergy){
      Printf("No value for the Energy found before start of run, the Energy will remain invalid");
    }
    else {
      AliDCSArray* dcs = (AliDCSArray*)array->At(indexEnergy);
      energy = (Float_t)(TMath::Nint(((Double_t)(dcs->GetInt(0)))*120/1000)); // sqrt(s)/2 energy in GeV
      Printf("Energy value found = %d (at %f), converting --> sqrt(s)/2 = %f (GeV)", dcs->GetInt(0),dcs->GetTimeStamp(),energy);
    }
  }

  return energy;
}

AliLHCClockPhase* ProcessLHCClockPhase(TObjArray *beam1phase,TObjArray *beam2phase, Double_t timeEnd)
{
  //
  // Method to process LHC-Clock Phase data
  // Only the values between DAQ_time_start and DAQ_time_end are kept
  //
  AliLHCClockPhase *phaseObj = new AliLHCClockPhase;

  Bool_t foundBeam1Phase = kFALSE, foundBeam2Phase = kFALSE;
  const Float_t threshold = 0.050; // we store the measurement only in case they differ with more 50ps from the previous one 

  //TString timeCreatedStr = GetRunParameter("time_created");
  //Double_t timeCreated = timeCreatedStr.Atof();
  Double_t timeCreated = 0.0;

  Int_t nCounts = beam1phase->GetEntries();
  Printf("Beam1 phase measurements = %d\n",nCounts);
  if (nCounts ==0){
    Printf("WARNING: No beam1 LHC clock phase values found!");
    delete phaseObj;
    return NULL;
  }
  else{
    Double_t prevPhase = 0;
    for (Int_t i = 0; i < nCounts; i++){
      AliDCSArray *dcs = (AliDCSArray*)beam1phase->At(i);
      if (dcs){
	      //if (dcs->GetTimeStamp()>=timeStart && dcs->GetTimeStamp()<=timeEnd) {
	      if (dcs->GetTimeStamp()>=timeCreated && dcs->GetTimeStamp()<=timeEnd) {
	  if ((i == 0) || (i == (nCounts-1)) ||
	      !foundBeam1Phase ||
	      (TMath::Abs(dcs->GetDouble(0)-prevPhase) > threshold)) {
	    prevPhase = dcs->GetDouble(0);
	    foundBeam1Phase = kTRUE;
	    Printf("B1 Clk Phase = %f at TS = %f", (Float_t)dcs->GetDouble(0),dcs->GetTimeStamp());  
	    phaseObj->AddPhaseB1DP((UInt_t)dcs->GetTimeStamp(),(Float_t)dcs->GetDouble(0));
	  }
	}
      }
    }
    if (!foundBeam1Phase){
      Printf("ERROR: No beam1 LHC clock phase values found within the run!");
      delete phaseObj;
      return NULL;
    }
  }

  nCounts = beam2phase->GetEntries();
  Printf("Beam2 phase measurements = %d\n",nCounts);
  if (nCounts ==0){
    Printf("WARNING: No beam2 LHC clock phase values found!");
    delete phaseObj;
    return NULL;
  }
  else{
    Double_t prevPhase = 0;
    for (Int_t i = 0; i < nCounts; i++){
      AliDCSArray *dcs = (AliDCSArray*)beam2phase->At(i);
      if (dcs){
	if (dcs->GetTimeStamp()>=timeCreated && dcs->GetTimeStamp()<=timeEnd) {
	  if ((i == 0) || (i == (nCounts-1)) ||
	      !foundBeam2Phase ||
	      (TMath::Abs(dcs->GetDouble(0)-prevPhase) > threshold)) {
	    prevPhase = dcs->GetDouble(0);
	    foundBeam2Phase = kTRUE;
	    Printf("B2 Clk Phase = %f at TS = %f", (Float_t)dcs->GetDouble(0),dcs->GetTimeStamp());  
	    phaseObj->AddPhaseB2DP((UInt_t)dcs->GetTimeStamp(),(Float_t)dcs->GetDouble(0));
	  }
	}
      }
    }
    if (!foundBeam2Phase){
      Printf("ERROR: No beam2 LHC clock phase values found within the run!");
      delete phaseObj;
      return NULL;
    }
  }

  return phaseObj;
}

void MakeLHCDataEntry(char* storageUri="local://$ALICE_ROOT/../AliRoot/OCDB", Int_t firstRun=0, Int_t lastRun=999999999)
{
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(storageUri);

  // Get time start from the simulated LHCData file
  Double_t timeStart = 0.0;
  Double_t timeEnd = 1.0e+10;

  TString fileName(gSystem->ExpandPathName("$ALICE_ROOT/../AliRoot/GRP/ShuttleInput/testShuttle_GRP_run_number_testShuttle_data.txt"));
  Printf("Getting the file %s", fileName.Data());

  const Int_t fgknLHCDP = 9;   // number of dcs dps from LHC data
  const char* fgkLHCDataPoints[fgknLHCDP] = {
    "LHC_Beam_Energy",
    "LHC_MachineMode",
    "LHC_BeamMode",
    "LHC_Beams_Particle_Type",
    "BPTX_Phase_Shift_B1",
    "BPTX_Phase_Shift_B2",
    "LHC_Particle_Type_B1",
    "LHC_Particle_Type_B2",
    "LHC_Data_Quality_Flag"
  };

  AliGRPObject *grpobj = new AliGRPObject();
	// grpobj->SetBeamEnergyIsSqrtSHalfGeV(); // new format
  //
  //Getting the LHC Data from DCS FXS
  //
  AliLHCReader lhcReader;

  // Processing data to be put in AliGRPObject

		// Energy
		Printf("*************Energy ");
		TObjArray* energyArray = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[0]);
		if (energyArray){			
			Float_t energy = ProcessEnergy(energyArray,timeStart);
			if (energy != -1.) {
				grpobj->SetBeamEnergy(energy);
				grpobj->SetBeamEnergyIsSqrtSHalfGeV(kTRUE);
			}
			delete energyArray;
		}
		else {
			AliError("Energy not found in LHC Data file!!!");
		}	

  Double_t timeBeamModeEnd = timeEnd;        // max validity for Beam Mode 
  Double_t timeMachineModeEnd = timeEnd;     // max validity for Machine Mode
  Double_t timeBeamEnd = timeEnd;            // max validity for Beam Type
  Double_t timeBeamTypeEnd[2] = {timeEnd, timeEnd}; // max validity for Beam Type1,2
  Double_t timeBeamModeStart = -1;    // min validity for Beam Mode
  Double_t timeMachineModeStart = -1; // min validity for Machine Mode
  Double_t timeBeamStart = -1;        // min validity for Beam Type
  Double_t timeBeamTypeStart[2] = {-1,-1};        // min validity for Beam Type1,2
  Int_t indexBeamMode = -1;                  // index of measurement used to set Beam Mode
  Int_t indexMachineMode = -1;               // index of measurement used to set Machine Mode
  Int_t indexBeam = -1;                      // index of measurement used to set Beam Type
  Int_t indexBeamType[2] = {-1, -1};                      // index of measurement used to set Beam Type1,2
  Bool_t foundBeamModeStart = kFALSE;        // flag to be set in case an entry for the Beam Mode is found before (or at) SOR
  Bool_t foundMachineModeStart = kFALSE;     // flag to be set in case an entry for the Machine Mode is found before (or at) SOR
  Bool_t foundBeamStart = kFALSE;            // flag to be set in case an entry for the Beam Type is found before (or at) SOR
  Bool_t foundBeamTypeStart[2] = {kFALSE, kFALSE};            // flag to be set in case an entry for the Beam Type1,2 is found before (or at) SOR
  Bool_t flagBeamMode = kFALSE;  //flag set true if a changed occurred in BeamMode
  Bool_t flagMachineMode = kFALSE;  //flag set true if a changed occurred in MachineMode
  Bool_t flagBeam = kFALSE;  //flag set true if a changed occurred in BeamType
  Bool_t flagBeamType[2] = {kFALSE, kFALSE};  //flag set true if a changed occurred in BeamType1,2

  Double_t arrayTimes[5]={2.E9, 2.E9, 2.E9, 2.E9, 2.E9}; // array to keep track of the times of the possible changes of the LHC DPs; each entry set to Wed May 18 2033, 03:33:20 GMT (ALICE should not be running anymore...)
  // arrayTimes elements order correspond to the one used in the array of the strings fgkLHCDataPoints, i.e.:
  // arrayTimes[0] --> MachineMode
  // arrayTimes[1] --> BeamMode
  // arrayTimes[2] --> BeamType (when written together)
  // arrayTimes[3] --> BeamType1 (when written separate)
  // arrayTimes[4] --> BeamType2 (when written separate)

  // BeamMode
  Printf("*************BeamMode (LHCState) ");
  TObjArray* beamModeArray = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[2]);
  Int_t nBeamMode = -1;
  if (beamModeArray){	
    nBeamMode = beamModeArray->GetEntries();	
    if (nBeamMode==0){
      Printf("Found zero entries for the Beam Mode, leaving it empty");
    }
    else{
      for (Int_t iBeamMode = 0; iBeamMode<nBeamMode; iBeamMode++){
        AliDCSArray* beamMode = (AliDCSArray*)beamModeArray->At(iBeamMode);
        if (beamMode){
          if (beamMode->GetTimeStamp()<=timeStart && beamMode->GetTimeStamp()>=timeBeamModeStart){// taking always the very last entry: of two measurements have the same timestamp, the last one is taken
            timeBeamModeStart = beamMode->GetTimeStamp();
            indexBeamMode = iBeamMode;
            foundBeamModeStart = kTRUE;
          }
          else {
            break;

          }
        }
      }
      if (!foundBeamModeStart){
        Printf("No value for the Beam Mode found before start of run, the Beam Mode will remain empty");
      }
      else {
        AliDCSArray* beamMode = (AliDCSArray*)beamModeArray->At(indexBeamMode);
        TObjString* beamModeString = beamMode->GetStringArray(0);
        Printf(Form("LHC State (corresponding to BeamMode) = %s (set at %f)",(beamModeString->String()).Data(),beamMode->GetTimeStamp()));
        grpobj->SetLHCState(beamModeString->String());
        if (indexBeamMode < nBeamMode-1){
          AliDCSArray* beamMode1 = (AliDCSArray*)beamModeArray->At(indexBeamMode+1);
          if (beamMode1){
            if (beamMode1->GetTimeStamp()<=timeStart){
              Printf("ERROR: you did not choose the correct value! there is still something before (or at) SOR, but later than this!");
            }
            else if (beamMode1->GetTimeStamp()>timeStart && beamMode1->GetTimeStamp()<=timeEnd){
              timeBeamModeEnd = beamMode1->GetTimeStamp();
              TObjString* beamModeString1 = beamMode1->GetStringArray(0);
              TString bmString0 = beamModeString->String();
              TString bmString1 = beamModeString1->String();
              if (bmString0.CompareTo(bmString1.Data(),TString::kIgnoreCase) == -1){
                Printf("WARNING: The beam mode changed from %s to %s during the run at timestamp %f! Setting it to %s and keeping track of the time of the change to set MaxTimeLHCValidity afterward",bmString0.Data(), bmString1.Data(), timeBeamModeEnd, bmString0.Data());
                flagBeamMode = kTRUE;
                arrayTimes[1]=timeBeamModeEnd;

              }
            }
          }
          else {
            Printf("Invalid pointer for the first entry for Beam Mode after the first valid one, not considering anything after what has already been found");
          }
        }
      }
    }
    delete beamModeArray;
  }
  else{
    Printf("ERROR: Beam mode array not found in LHC Data file!!!");
  }
		
  // MachineMode
  Printf("*************MachineMode ");
  TObjArray* machineModeArray = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[1]);
  Int_t nMachineMode = -1;
  if (machineModeArray){
    nMachineMode = machineModeArray->GetEntries();
    if (nMachineMode==0){
      Printf("No Machine Mode found, leaving it empty");
    }
    else{
      for (Int_t iMachineMode = 0; iMachineMode<nMachineMode; iMachineMode++){
        AliDCSArray* machineMode = (AliDCSArray*)machineModeArray->At(iMachineMode);
        if (machineMode){
          if (machineMode->GetTimeStamp()<=timeStart && machineMode->GetTimeStamp()>=timeMachineModeStart){// taking always the very last entry: of two measurements have the same timestamp, the last one is taken
            timeMachineModeStart = machineMode->GetTimeStamp();
            indexMachineMode = iMachineMode;
            foundMachineModeStart = kTRUE;
          }
          else{
            break;
          }
        }
      }
      if (!foundMachineModeStart){
        Printf("No value for the Machine Mode found before start of run, the Machine Mode will remain empty");
      }
      else {
        AliDCSArray* machineMode = (AliDCSArray*)machineModeArray->At(indexMachineMode);
        TObjString* machineModeString = machineMode->GetStringArray(0);
        Printf(Form("MachineMode = %s (set at %f)",(machineModeString->String()).Data(),machineMode->GetTimeStamp()));
        grpobj->SetMachineMode(machineModeString->String());
        if (indexMachineMode < nMachineMode-1){
          AliDCSArray* machineMode1 = (AliDCSArray*)machineModeArray->At(indexMachineMode+1);
          if (machineMode1){
            if (machineMode1->GetTimeStamp()>timeStart && machineMode1->GetTimeStamp()<=timeEnd){
              timeMachineModeEnd = machineMode1->GetTimeStamp();
              TObjString* machineModeString1 = machineMode1->GetStringArray(0);
              TString mmString0 = machineModeString->String();
              TString mmString1 = machineModeString1->String();
              if (mmString0.CompareTo(mmString1.Data(),TString::kIgnoreCase) == -1){
                Printf("WARNING: The machine mode changed from %s to %s during the run at timestamp %f! Setting it to %s and keeping track of the time of the change to set MaxTimeLHCValidity afterward",mmString0.Data(),mmString1.Data(),timeMachineModeEnd,mmString0.Data());
                flagMachineMode = kTRUE;
                arrayTimes[0]=timeMachineModeEnd;
              }
            }
          }
          else {
            Printf("Invalid pointer for the first entry for Machine Mode after the first valid one, not considering anything after what has already been found");
          }
        }
      }
    }
    delete machineModeArray;
  }
  else{
    Printf("ERROR: Machine mode array not found in LHC Data file!!!");
  }

  // BeamType1 and BeamType2 - both put in the same string
  Printf("*************BeamType ");
  TObjArray* beamArray = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[3]);
  if (beamArray){			
    Int_t nBeam = beamArray->GetEntries();
    if (nBeam==0){
      Printf("No Beam Type found, leaving it empty");
    }
    else{
      for (Int_t iBeam = 0; iBeam<nBeam; iBeam++){
        AliDCSArray* beam = (AliDCSArray*)beamArray->At(iBeam);
        if (beam){
          if (beam->GetTimeStamp()<=timeStart && beam->GetTimeStamp()>=timeBeamStart){// taking always the very last entry: of two measurements have the same timestamp, the last one is taken
            timeBeamStart = beam->GetTimeStamp();
            indexBeam = iBeam;
            foundBeamStart = kTRUE;
          }
          else{
            break;
          }
        }
      }
      if (!foundBeamStart){
        Printf("No value for the Beam Type found before start of run, the (common) Beam Type will remain empty");
      }
      else {
        AliDCSArray* beam = (AliDCSArray*)beamArray->At(indexBeam);
        TObjString* beamString = beam->GetStringArray(0);
        TString beamType = beamString->String();
        Printf(Form("Beam Type = %s",beamType.Data()));	
        if (beamType.CompareTo("PROTON",TString::kIgnoreCase) == 0){
          Printf("Setting beam type to p-p");
          grpobj->SetBeamType("p-p");
        }
        else { // if there is no PROTON beam, we suppose it is Pb, and we put A-A
          Printf("Setting beam type to A-A");
          grpobj->SetBeamType("A-A");
        }
        /*
           else if (beamType.CompareTo("LEAD82",TString::kIgnoreCase) == 0){
           Printf("Setting beam type to Pb-Pb");
           grpobj->SetBeamType("Pb-Pb");
           }
           else{
           Printf("ERROR: Beam Type not known, leaving it empty");
           }
           */
        if (indexBeam < nBeam-1){
          AliDCSArray* beam1 = (AliDCSArray*)beamArray->At(indexBeam+1);
          if (beam1){
            if (beam1->GetTimeStamp()>timeStart && beam1->GetTimeStamp()<=timeEnd){
              timeBeamEnd = beam1->GetTimeStamp();
              TObjString* beamString1 = beam1->GetStringArray(0);
              TString beamType1 = beamString1->String();
              if (beamType.CompareTo(beamType1.Data(),TString::kIgnoreCase) == -1){
                Printf("WARNING: The Beam Type changed from %s to %s during the run at timestamp %f! Setting it to %s and keeping track of the time of the change to set MaxTimeLHCValidity afterward",beamType.Data(),(beamString1->String()).Data(),timeBeamEnd,beamType.Data());
                flagBeam = kTRUE;
                arrayTimes[2] = timeBeamEnd;
              }
            }
          }
          else {
            Printf("Invalid pointer for the first entry for Beam Type after the first valid one, not considering anything after what has already been found");
          }
        }
      }
    }
    delete beamArray;
  }
  else{
    Printf("ERROR: Beam Type array not found in LHC Data file!!!");
  }		

  // BeamType1 and BeamType2 - in separete string
  Printf("*************BeamType, 1 and 2 ");
  Int_t indexBeamTypeString = 6;  // index of the string with the alias of BeanType1 in the array fgkLHCDataPoints
  TString combinedBeamType = "-";  // combined beam type, built from beam type 1 and beam type 2
  TString combinedBeamTypeFromLHC = "-";  // combined beam type, built from beam type 1 and beam type 2 AS SENT FROM LHC
  for (Int_t ibeamType = 0; ibeamType<2; ibeamType++){
    beamArray = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[indexBeamTypeString+ibeamType]);
    if (beamArray){			
      Int_t nBeam = beamArray->GetEntries();
      if (nBeam==0){
        Printf(Form("No Beam Type %s found, leaving it empty",fgkLHCDataPoints[indexBeamTypeString+ibeamType]));
      }
      else{
        for (Int_t iBeam = 0; iBeam<nBeam; iBeam++){
          AliDCSArray* beam = (AliDCSArray*)beamArray->At(iBeam);
          if (beam){
            if (beam->GetTimeStamp()<=timeStart && beam->GetTimeStamp()>=timeBeamTypeStart[ibeamType]){// taking always the very last entry: of two measurements have the same timestamp, the last one is taken
              timeBeamTypeStart[ibeamType] = beam->GetTimeStamp();
              indexBeamType[ibeamType] = iBeam;
              foundBeamTypeStart[ibeamType] = kTRUE;
            }
            else{
              break;
            }
          }
        }
        if (!foundBeamTypeStart[ibeamType]){
          Printf(Form("No value for the Beam Type %s found before start of run, the Beam Type %d will remain empty", fgkLHCDataPoints[indexBeamTypeString+ibeamType], ibeamType));
        }
        else {
          AliDCSArray* beam = (AliDCSArray*)beamArray->At(indexBeam);
          TObjString* beamString = beam->GetStringArray(0);
          TString beamType = beamString->String();
          Printf(Form("Beam Type (for %s) = %s", fgkLHCDataPoints[indexBeamTypeString+ibeamType], beamType.Data()));
          TString singleBeam = ParseBeamTypeString(beamType,ibeamType);
          Printf(Form("Single Beam Type for beam %d set to %s", ibeamType, singleBeam.Data()));
          grpobj->SetSingleBeamType(ibeamType, singleBeam);
          if (beamType.CompareTo("PROTON",TString::kIgnoreCase) == 0){
            Printf(Form("Setting beam %d for combined beam type to p", ibeamType));
            if (ibeamType == 0) combinedBeamType.Prepend("p"); 
            else combinedBeamType.Append("p");
          }
          else { // if there is no PROTON beam, we suppose it is Pb, and we put A-A
            Printf(Form("Setting beam %d for combined beam type to A",ibeamType));
            if (ibeamType == 0) combinedBeamType.Prepend("A");
            else combinedBeamType.Append("A");
          }
          if (ibeamType == 0) combinedBeamTypeFromLHC.Prepend(beamType); 
          else combinedBeamTypeFromLHC.Append(beamType);
          /*
             else if (beamType.CompareTo("LEAD82",TString::kIgnoreCase) == 0){
             Printf("Setting beam type to Pb-Pb");
             grpobj->SetSingleBeamType(ibeamType, "Pb-Pb");
             }
             else{
             Printf("ERROR: Beam Type not known, leaving it empty");
             }
             */
          if (indexBeamType[ibeamType] < nBeam-1){
            AliDCSArray* beam1 = (AliDCSArray*)beamArray->At(indexBeam+1);
            if (beam1){
              if (beam1->GetTimeStamp()>timeStart && beam1->GetTimeStamp()<=timeEnd){
                timeBeamTypeEnd[ibeamType] = beam1->GetTimeStamp();
                TObjString* beamString1 = beam1->GetStringArray(0);
                TString beamType1 = beamString1->String();
                if (beamType.CompareTo(beamType1.Data(),TString::kIgnoreCase) == -1){
                  Printf("WARNING: The Beam Type for %s changed from %s to %s during the run at timestamp %f! Setting it to %s and keeping track of the time of the change to set MaxTimeLHCValidity afterward",fgkLHCDataPoints[indexBeamTypeString+ibeamType],beamType.Data(),(beamString1->String()).Data(),timeBeamEnd,beamType.Data());
                  flagBeamType[ibeamType] = kTRUE;
                  arrayTimes[3+ibeamType] = timeBeamTypeEnd[ibeamType];
                }
              }
            }
            else {
              Printf(Form("Invalid pointer for the first entry for Beam Type %s after the first valid one, not considering anything after what has already been found",fgkLHCDataPoints[indexBeamTypeString+ibeamType]));
            }
          }
        }
      }
      delete beamArray;
    }
    else{
      AliError(Form("Beam Type %s array not found in LHC Data file!!!",fgkLHCDataPoints[indexBeamTypeString+ibeamType]));
    }		
  }
  Printf(Form("Setting combined beam type to %s",combinedBeamType.Data()));
  grpobj->SetBeamType(combinedBeamType);
  Printf(Form("Setting combined beam type form LHC to %s",combinedBeamTypeFromLHC.Data()));
  grpobj->SetBeamTypeFromLHC(combinedBeamTypeFromLHC);

  // Setting minTimeLHCValidity
  if (flagBeamMode == kTRUE || flagMachineMode == kTRUE || flagBeam == kTRUE || flagBeamType[0] == kTRUE || flagBeamType[1] == kTRUE){ 
    Double_t minTimeLHCValidity= TMath::MinElement(5,arrayTimes);
    Printf("WARNING: Setting MaxTimeLHCValidity to %f",minTimeLHCValidity);
    grpobj->SetMaxTimeLHCValidity(minTimeLHCValidity);
  }
  /* 
  // Old way to determine the Maximum Time during which the LHC info is valid
  if (timeBeamModeEnd!=0 || timeMachineModeEnd!=0 || timeBeamEnd !=0){
  Double_t minTimeLHCValidity;
  if (flagBeamMode == kFALSE && flagMachineMode == kFALSE && flagBeam == kTRUE){ // flagBeam only true --> it is the only one that changed
  minTimeLHCValidity = timeBeamEnd;
  }
  else if (flagBeamMode == kFALSE && flagMachineMode == kTRUE && flagBeam == kFALSE){ // flagMachineMode only true
  minTimeLHCValidity = timeMachineModeEnd;
  }
  else if (flagBeamMode == kTRUE && flagMachineMode == kFALSE && flagBeam == kFALSE){ // flagBeamMode only true
  minTimeLHCValidity = timeBeamModeEnd;
  }
  else if (flagBeamMode == kFALSE && flagMachineMode == kTRUE && flagBeam == kTRUE){ // flagBeam and flagMachineMode only true
  minTimeLHCValidity= TMath::Min(timeBeamEnd,timeMachineModeEnd);
  }
  else if (flagBeamMode == kTRUE && flagMachineMode == kFALSE && flagBeam == kTRUE){ // flagBeam and flagBeamMode only true
  minTimeLHCValidity= TMath::Min(timeBeamEnd,timeBeamModeEnd);
  }
  else if (flagBeamMode == kTRUE && flagMachineMode == kTRUE && flagBeam == kFALSE){ // flagMachineMode and flagBeamMode only true
  minTimeLHCValidity= TMath::Min(timeMachineModeEnd,timeBeamModeEnd);
  }
  else {
  Double_t arrayTimes[3] = {timeBeamModeEnd,timeMachineModeEnd,timeBeamEnd};// flagMachineMode and flagBeamMode and flagBeam 
  minTimeLHCValidity= TMath::MinElement(3,arrayTimes);
  }
  Printf("WARNING: Setting MaxTimeLHCValidity to %f",minTimeLHCValidity));
  grpobj->SetMaxTimeLHCValidity(minTimeLHCValidity);
  }
  */

  // Data Quality Flag --> storing start and end values of periods within the run during which the value was found to be FALSE
  Printf("*************Data Quality Flag ");
  TObjArray* dataQualityArray = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[8]);
  Int_t nDataQuality = -1;
  Double_t timeDataQualityStart = -1; // min validity for Data Quality Flag
  Int_t indexDataQuality = -1;               // index of first measurement used to set Data Quality Flag
  Bool_t foundDataQualityStart = kFALSE;     // flag to be set in case an entry for the Data Quality Flag is found before (or at) SOR

  if (dataQualityArray){
    nDataQuality = dataQualityArray->GetEntries();
    if (nDataQuality==0){
      Printf("No Data Quality Flag found, leaving it empty");
    }
    else{
      for (Int_t iDataQuality = 0; iDataQuality<nDataQuality; iDataQuality++){
        AliDCSArray* dataQuality = (AliDCSArray*)dataQualityArray->At(iDataQuality);
        if (dataQuality){
          if (dataQuality->GetTimeStamp()<=timeStart && dataQuality->GetTimeStamp()>=timeDataQualityStart){// taking always the very last entry: if two measurements have the same timestamp, the last one is taken
            timeDataQualityStart = dataQuality->GetTimeStamp();
            indexDataQuality = iDataQuality;
            foundDataQualityStart = kTRUE;
          }
          else{
            // we suppose here that if the first measurement is not before SOR, then none will be (they MUST be in chronological order!!!) 
            break;
          }
        }
      }
      if (!foundDataQualityStart){
        // The Data Quality Flag should be found and TRUE at the start of the run. For the time being, if it is not found, don't do anything, but it means there is a problem..
        Printf("No value for the Data Quality Flag found before start of run, the Data Quality Flag will remain empty");
      }
      else {
        // counting how many FALSE values there are
        Bool_t foundEndOfFalse = kFALSE;
        Int_t nFalse = 0;
        for (Int_t iDataQuality = indexDataQuality; iDataQuality < nDataQuality; iDataQuality ++){
          AliDCSArray* dataQuality = (AliDCSArray*)dataQualityArray->At(iDataQuality);
          Printf("dataQuality->GetTimeStamp() = %f, timeDataQualityStart = %f, timeEnd = %f", dataQuality->GetTimeStamp(), timeDataQualityStart, timeEnd );
          if (dataQuality->GetTimeStamp()>=timeDataQualityStart && dataQuality->GetTimeStamp()<=timeEnd){ // considering only values between the first valid and the end of the run
            Bool_t dataQualityFlag = dataQuality->GetBool(0);
            Printf("DataQuality = %d (set at %f)",(Int_t)dataQualityFlag,dataQuality->GetTimeStamp());
            if (dataQualityFlag != kTRUE){
              if (iDataQuality == indexDataQuality) {  // the first Data Quality value should be TRUE, but ignoring the problem now...
                Printf("ERROR: The first value for the Data Quality MUST be TRUE! Ignoring for now...");
              }
              nFalse++;
            }
          }
        }

        Printf(Form("Found %d FALSE values for the Data Quality Flag",nFalse));
        Double_t falses[nFalse*2];  // dimensioning this to the maximum possible, as if each false value was followed by a true one --> the false periods correspond to the number of falses

        Int_t iDataQuality = indexDataQuality;
        if (nFalse > 0){
          Int_t iFalse = 0;
          // filling the info about the periods when the flag was set to FALSE
          // starting, like for the other DPS, from the measurement closest to SOR (the index of which is iDataQuality)
          while (iDataQuality < nDataQuality){
            Printf("iDataQuality = %d",iDataQuality);
            AliDCSArray* dataQuality = (AliDCSArray*)dataQualityArray->At(iDataQuality);
            if (dataQuality->GetTimeStamp()>=timeDataQualityStart && dataQuality->GetTimeStamp()<=timeEnd){ // considering only values between the first valid and the end of the run
              Bool_t dataQualityFlag = dataQuality->GetBool(0);
              Printf("DataQuality = %d (set at %f)",(Int_t)dataQualityFlag,dataQuality->GetTimeStamp());
              if (dataQualityFlag == kTRUE){
                // found TRUE value, continuing
                iDataQuality++;
                continue;
              }
              else{
                /*
                // the check was already done before
                if (iDataQuality == indexDataQuality) {  // the first Data Quality value should be TRUE, but ignoring the problem now...
                Printf("ERROR: The first value for the Data Quality MUST be TRUE! Ignoring for now...");
                }
                */
                falses[iFalse*2] = dataQuality->GetTimeStamp();
                foundEndOfFalse = kFALSE;
                Int_t iDataQualityNext = iDataQuality+1;
                while (iDataQualityNext < nDataQuality){
                  AliDCSArray* dataQualityNext = (AliDCSArray*)dataQualityArray->At(iDataQualityNext);
                  if (dataQualityNext->GetTimeStamp()>timeDataQualityStart && dataQualityNext->GetTimeStamp()<=timeEnd && dataQualityNext->GetTimeStamp() > dataQuality->GetTimeStamp()){ // considering only values between the first valid and the end of the run, and subsequent to the current value
                    Bool_t dataQualityFlagNext = dataQualityNext->GetBool(0);
                    Printf("DataQualityNext = %d (set at %f)",(Int_t)dataQualityFlagNext,dataQualityNext->GetTimeStamp());
                    if (dataQualityFlagNext == kTRUE){
                      // found TRUE value, first FALSE period completed
                      foundEndOfFalse = kTRUE;
                      falses[iFalse*2+1] = dataQualityNext->GetTimeStamp();
                      iFalse++;
                      break;
                    }
                    iDataQualityNext++;
                  }
                }
                if (!foundEndOfFalse) {
                  Printf("Please, note that the last FALSE value lasted until the end of the run");
                  falses[iFalse*2+1] = timeEnd;
                  iFalse++;
                  break;
                }
                iDataQuality = iDataQualityNext+1;
              }
            }
          }
          grpobj->SetNFalseDataQualityFlag(iFalse);
          grpobj->SetFalseDataQualityFlagPeriods(falses);
        }
      }
    }
    delete dataQualityArray;
  }
  else{
    Printf("ERROR: Data Quality Flag array not found in LHC Data file!!!");
  }

  // Processing data to go to AliLHCData object
  AliLHCData* dt = new AliLHCData(fileName.Data(),timeStart,timeEnd);
  // storing AliLHCData in OCDB
  if (dt){
    Printf(Form("Filled %d records to AliLHCData object",dt->GetData().GetEntriesFast()));
    AliCDBMetaData md;
    md.SetResponsible("Ruben Shahoyan");
    md.SetComment("LHC data from the GRP preprocessor.");
    Bool_t result = kTRUE;
    AliCDBId id("GRP/GRP/LHCData", 0, AliCDBRunRange::Infinity());
    result = cdb->Put(dt, id, &md); 
    delete dt;
    if (!result){
      Printf("Problems in storing LHC Data - but not going into Error");
    }
  }

  // processing LHC Phase

  TObjArray *beam1phase = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[4]);
  TObjArray *beam2phase = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[5]);
  if (beam1phase == 0x0 || beam2phase == 0x0){
    Printf(Form("Problems in retrieving LHC Clock data from LHC file"));
    return 4;
  }			
  AliLHCClockPhase *phaseObj = ProcessLHCClockPhase(beam1phase,beam2phase,timeEnd);
  delete beam1phase;
  delete beam2phase;
  if (phaseObj){
    Printf(Form("LHC Phase found"));
    AliCDBMetaData mdPhase;
    mdPhase.SetResponsible("Cvetan Cheshkov");
    mdPhase.SetComment("LHC Clock Phase");
    Bool_t result = kTRUE;
    AliCDBId id("GRP/Calib/LHCClockPhase", 0, AliCDBRunRange::Infinity());
    result = cdb->Put(phaseObj, id, &mdPhase); 
    delete phaseObj;
    if (!result) return 3;
  }
  else return 4;

  return 0;
}
