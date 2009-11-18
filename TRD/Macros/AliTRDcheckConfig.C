//===================================================================================
// This is a macro to analyze TRD/Calib/DCS OCDB objects either
// from the grid for a given run number or from a local object.
// If you want to analyze data from the grid, please don't forget to
// have a valid alien token initialized and the file /tmp/gclient_env_$UID source'd.
//
// Arguments:
// The first argument is the runnumber (this is ignored in case of a local file),
// the second is a string that needs to contain either "grid" or "local". Further
// you can add either verbose or quiet to that string. If you don't, you'll be asked
// for all stuff individually wether you want to see it or not
// the thrid argument is the number of the ROC you (eventually) want to dump its data
// member of.
// The fourth one is the path and name of the local file you might want to look at.
//
// So the simplest way to use this macro is if you want to check the output of a given
// run from the OCDB:
// .x AliTRDcheckConfig.C(60111)
//
// An example for quickly checking a local file:
// .x AliTRDcheckConfig.C(0, "local quiet", 533, "$ALICE_ROOT/TRD/Calib/DCS/Run0_999999999_v0_s0.root")
//
// Please contact Frederick Kramer in case of problems
//===================================================================================

// This is the path one needs to change if the year is no longer 2009 
// and the runnumber cannot be found
TString alienOcdbPath("alien://folder=/alice/data/2009/OCDB/");

// Do not make changes below here unless you know what your doing

const Int_t nROC = 540;
const Int_t nROB = 8;
const Int_t nMCM = 18;
const Int_t cArraySize = 1000;

Bool_t errors = false;

Int_t AnalyzeArray(Int_t states[cArraySize], Int_t occur[cArraySize]) {
  long long srtIndx[cArraySize] = 0;
  
  TMath::Sort(cArraySize, occur, srtIndx);

  Int_t totalSum = 0, subSum = 0, iIndex = 0;
  for (Int_t i=0; i<cArraySize; i++) totalSum += occur[i];
  
  cout << "    The majority ("<< occur[srtIndx[0]] << " of " 
       << totalSum <<") is: " << states[srtIndx[0]] << endl;
  subSum = occur[srtIndx[0]];
  while (totalSum != subSum) {
    if (++iIndex > 999) {
      cout << "E : out of bounds." << endl;
      break;
    }
    Printf("    Next: %7d (%d)", states[srtIndx[iIndex]], occur[srtIndx[iIndex]]);
    subSum += occur[srtIndx[iIndex]];
  }
  return states[srtIndx[0]];
}



void FillItemInArray(Int_t states[cArraySize], Int_t occur[cArraySize], Int_t item, Bool_t allowNeg) {
  for (Int_t iArrPos=0; iArrPos<cArraySize; iArrPos++) {
    // if allowNeg is set then we change the number indicating that the item ws not set from -1 to -100
    // so that small negitive numbers can be sorted too
    if ((allowNeg && item == -100000) || (!allowNeg && item == -1)) break; // value not set
    if (states[iArrPos] == item) {
      occur[iArrPos]++;
      break;
    } else if (occur[iArrPos] == 0) {
      states[iArrPos] = item;
      occur[iArrPos]++;
      break;
    }
  }
}

void GetMajoritys(AliTRDCalDCS* calDCSObj) {
  
  Int_t gsmStates[cArraySize] = {0}, gsmOccur[cArraySize] = {0};
  Int_t nimStates[cArraySize] = {0}, nimOccur[cArraySize] = {0};
  Int_t nevStates[cArraySize] = {0}, nevOccur[cArraySize] = {0};
  Int_t nptStates[cArraySize] = {0}, nptOccur[cArraySize] = {0};
  
  for (Int_t i=0; i<cArraySize; i++) {
    gsmStates[i] = 0; gsmOccur[i]  = 0;
    nimStates[i] = 0; nimOccur[i]  = 0;
    nevStates[i] = 0; nevOccur[i]  = 0;
    nptStates[i] = 0; nptOccur[i]  = 0;  
  }
  
  for (Int_t i=0; i<nROC && i<calDCSObj->GetFEEArr()->GetSize(); i++) {
    AliTRDCalDCSFEE *idcsfee = calDCSObj->GetCalDCSFEEObj(i);
    if ((idcsfee == NULL) || (idcsfee->GetStatusBit() != 0)) continue;
    for (Int_t j=0; j<nROB; j++) {
      for (Int_t k=0; k<nMCM; k++) {
	Int_t igsm = idcsfee->GetMCMGlobalState(j,k);
	Int_t inim = idcsfee->GetMCMStateNI(j,k);
	Int_t inev = idcsfee->GetMCMEventCnt(j,k);
	Int_t inpt = idcsfee->GetMCMPtCnt(j,k);
	
	FillItemInArray(gsmStates, gsmOccur, igsm, false); 
	FillItemInArray(nimStates, nimOccur, inim, false); 
	FillItemInArray(nevStates, nevOccur, inev, false); 
	FillItemInArray(nptStates, nptOccur, inpt, false); 
      }
    }
  }
  
  cout << "I : Global MCM state statistics:" << endl;
  AnalyzeArray(gsmStates, gsmOccur);
  cout << "I : Network interface state statistics:" << endl;
  AnalyzeArray(nimStates, nimOccur);
  cout << "I : MCM Event counter reading statistics:" << endl;
  AnalyzeArray(nevStates, nevOccur);
  cout << "I : MCM PreTrigger counter reading statistics:" << endl;
  AnalyzeArray(nptStates, nptOccur);
  
  return;
}



void GetMajorityDifferences(AliTRDCalDCS* calDCSObj, AliTRDCalDCS* calDCSObj2) {
  
  Int_t gsmStates[cArraySize] = {0}, gsmOccur[cArraySize] = {0};
  Int_t nimStates[cArraySize] = {0}, nimOccur[cArraySize] = {0};
  Int_t nevStates[cArraySize] = {0}, nevOccur[cArraySize] = {0};
  Int_t nptStates[cArraySize] = {0}, nptOccur[cArraySize] = {0};
  
  for (Int_t i=0; i<cArraySize; i++) {
    gsmStates[i] = 0; gsmOccur[i]  = 0;
    nimStates[i] = 0; nimOccur[i]  = 0;
    nevStates[i] = 0; nevOccur[i]  = 0;
    nptStates[i] = 0; nptOccur[i]  = 0;  
  }
  
  for (Int_t i=0; i<nROC && i<calDCSObj->GetFEEArr()->GetSize() && i<calDCSObj2->GetFEEArr()->GetSize(); i++) {
    AliTRDCalDCSFEE *idcsfee = calDCSObj->GetCalDCSFEEObj(i);
    AliTRDCalDCSFEE *idcsfee2 = calDCSObj2->GetCalDCSFEEObj(i);
    if ((idcsfee == NULL) || (idcsfee2 == NULL) || 
      (idcsfee->GetStatusBit() != 0) /*|| (idcsfee2->GetStatusBit() != 0)*/) continue;
    for (Int_t j=0; j<nROB; j++) {
      for (Int_t k=0; k<nMCM; k++) {
	Int_t igsm = idcsfee->GetMCMGlobalState(j,k) - idcsfee2->GetMCMGlobalState(j,k);
	Int_t inim = idcsfee->GetMCMStateNI(j,k)     - idcsfee2->GetMCMStateNI(j,k);
	Int_t inev = idcsfee2->GetMCMEventCnt(j,k)   - idcsfee->GetMCMEventCnt(j,k);
	Int_t inpt = idcsfee2->GetMCMPtCnt(j,k)      - idcsfee->GetMCMPtCnt(j,k);
	
	// if they were set to -1, it means they were not actauuly set
	// change -1 to -100 to mean they werent set since the above 
	// can give negitives
	if (idcsfee->GetMCMGlobalState(j,k) == -1 && igsm == 0) igsm =-100000;
	if (idcsfee->GetMCMStateNI(j,k) == -1 && inim == 0)     inim =-100000;
	if (idcsfee->GetMCMEventCnt(j,k) == -1 && inev == 0)    inev =-100000;
	if (idcsfee->GetMCMPtCnt(j,k) == -1 && inpt == 0)       inpt =-100000;
	
	FillItemInArray(gsmStates, gsmOccur, igsm, true); 
	FillItemInArray(nimStates, nimOccur, inim, true); 
	FillItemInArray(nevStates, nevOccur, inev, true); 
	FillItemInArray(nptStates, nptOccur, inpt, true); 
      }
    }
  }
  
  cout << "I : Global MCM state difference statistics:" << endl;
  AnalyzeArray(gsmStates, gsmOccur);
  cout << "I : Network interface state difference statistics:" << endl;
  AnalyzeArray(nimStates, nimOccur);
  cout << "I : MCM Event counter difference statistics:" << endl;
  if (AnalyzeArray(nevStates, nevOccur) < 1) {
    cout << "E : There should have been some events recorded, but there weren't" << endl;
    errors = true;
  }
  cout << "I : MCM PreTrigger counter difference statistics:" << endl;
  if (AnalyzeArray(nptStates, nptOccur) < 1) {
    cout << "E : There should have been some events recorded, but there weren't" << endl;
    errors = true;
  }
  
  return;
}


void AliTRDcheckConfig(Int_t runNr=0, char *pathfile="nopathgiven"){

  AliCDBEntry *entry=0;
  TString pathfilets(pathfile);

  // get the source
  if(pathfilets.Contains("nopathgiven")) {
    cout << "I : Accessing grid storage for run number " << runNr << endl;
    cout << "I : Get CDBManager instance." << endl;
    AliCDBManager *man = AliCDBManager::Instance();
    cout << "I : SetDefaultStorage." << endl;
    man->SetDefaultStorage(alienOcdbPath);
    cout << "I : Get OCDB Entry." << endl;
    entry = man->Get("TRD/Calib/DCS", runNr);
    if (entry == NULL) {
      cout << endl << "ERROR: Unable to get the AliTRDCalDCS object from the OCDB for run number " << runNr << endl << endl;
      cout << "If the run number is correct, it could be that the year is no longer 2009 and" << endl;
      cout << "the path where the objects is stored has changed, check the top of this macro " << endl;
      cout << "to change the path." << endl;
      return;
    }
  } else {
    cout << "I : Accessing local storage" << endl;
    TFile *f = new TFile(pathfile);
    if(f != NULL) {
      entry = (AliCDBEntry*) f->Get("AliCDBEntry");
    }
    else {
      cout << "E : Cannot open file" << endl;
      return;
    }
  }
  
  TObject *objectCDB = (TObject*)entry->GetObject();
  if (objectCDB->IsA()->InheritsFrom("TObjArray")) {
    TObjArray *objArrayCDB = (TObjArray*)entry->GetObject();
  }
  
  // the CalDCS object
  AliTRDCalDCS *caldcs;
  AliTRDCalDCS *caldcs2;

  Bool_t sorandeor = true;

  caldcs = (AliTRDCalDCS*) objArrayCDB->At(0);
  caldcs2 = (AliTRDCalDCS*) objArrayCDB->At(1);

  if (caldcs == NULL && caldcs2 == NULL) {
    cout << "E : Niether the start or end of run files were in the root file.";
    return;
  } else if (caldcs != NULL && caldcs2 == NULL) {
    cout << "E : The EOR file was not in the root file.";
    errors = true;
    sorandeor = false;
  } else if (caldcs == NULL && caldcs2 != NULL) {
    cout << "E : The SOR file was not in the root file.";
    errors = true;
    sorandeor = false;
    caldcs = caldcs2;
  }

  cout << endl << "============ Non responding ROC Summary: ============" << endl;
  TString bitfivestr = " ROCs with status bit 5. These havn't responded to communication\nattempts over DIM. Most probably they just were off this is ok.\n    DCS IDs: ";
  Int_t lengthfive = bitfivestr.Length();
  TString bitfourstr = " ROCs with status bit 4! BAD! This might be due to a communication problem between fxsproxy and the feeserver(s) \n    DCS IDs: ";
  Int_t lengthfour = bitfourstr.Length();
  TString bitthreestr = " ROCs with status bit 3! BAD! data from fee server was old or corrupt.\n    DCS IDs: ";
  Int_t lengththree = bitthreestr.Length();
  TString bittwostr = " ROCs with status bit 2. These have been in states in which they cannot be read out, e.g. Standby.\n    DCS IDs: ";
  Int_t lengthtwo = bittwostr.Length();
  TString bitonestr = " ROCs with status bit 1! BAD! This means the chamber(s) didn't respont even though is should have been in a good state.\n    DCS IDs: ";
  Int_t lengthone = bitonestr.Length();

  Int_t nSB1=0, nSB2=0, nSB3=0, nSB4=0, nSB5=0, nTot=0;
  for (Int_t i=0; i<nROC && i<caldcs->GetFEEArr()->GetSize(); i++) {
    AliTRDCalDCSFEE *idcsfee;
    idcsfee = caldcs->GetCalDCSFEEObj(i);
    if (idcsfee != NULL) {
      Int_t sb = idcsfee->GetStatusBit();
      if (sb == 5) { bitfivestr  += i; bitfivestr  += "  "; nSB5++; }
      if (sb == 4) { bitfourstr  += i; bitfourstr  += "  "; nSB4++; errors = true; }
      if (sb == 3) { bitthreestr += i; bitthreestr += "  "; nSB3++; errors = true; }
      if (sb == 2) { bittwostr   += i; bittwostr   += "  "; nSB2++; }
      if (sb == 1) { bitonestr   += i; bitonestr   += "  "; nSB1++; errors = true; }
      nTot += 1;
    }
  }

  if (lengthfive < bitfivestr.Length()) cout << nSB5 << bitfivestr.Data() << endl << endl;
  if (lengthfour < bitfourstr.Length()) cout << nSB4 << bitfourstr.Data() << endl << endl;
  if (lengththree < bitthreestr.Length()) cout << nSB3 << bitthreestr.Data() << endl << endl;
  if (lengthtwo < bittwostr.Length()) cout << nSB2 << bittwostr.Data() << endl << endl;
  if (lengthone < bitonestr.Length()) cout << nSB1 << bitonestr.Data() << endl << endl;
  
  cout << "The remaining " << nTot-(nSB1+nSB2+nSB3+nSB4+nSB5) << " ROCs responded correctly in the start of run."<<endl;

  Int_t nChanged=0, nTot=0;
  for (Int_t i=0; i<nROC && i<caldcs->GetFEEArr()->GetSize(); i++) {
    AliTRDCalDCSFEE *idcsfee;
    idcsfee = caldcs->GetCalDCSFEEObj(i);
    idcsfee2 = caldcs2->GetCalDCSFEEObj(i);
    if (idcsfee != NULL && idcsfee2 != NULL) {
      Int_t sbd = idcsfee->GetStatusBit() - idcsfee2->GetStatusBit();
      if (sbd != 0) { 
	cout << "ROC " << i << " changed from state " << idcsfee->GetStatusBit() << " at start of the run to "  << idcsfee2->GetStatusBit() << " at the end of the run." << endl;
	nChanged++; 
      }
      nTot += 1;
    }
  }
  
  if (nChanged == 0) {
    cout << "No ROCs changed state between the start and end of the run" << endl;
  } else {
    cout << "E : " << nChanged << " out of " << nTot << " ROCs changed state during the run" << endl;
    errors = true; 
  }

  cout << endl << "============ Statistics from RSTATE: ============" << endl;
  cout<<"I : The majority entry is given as well as all other values," << endl;
  cout<<"    sorted according to their occurrence." << endl << endl;
  GetMajoritys(caldcs);
  if (sorandeor) GetMajorityDifferences(caldcs,caldcs2);

  cout << endl << "============ Global Configuraton: ============" << endl;
  cout<<"I : Anything not listed is not set, mixed numbers are indicated with a" << endl;
  cout<<"    value of -2 and strings are set to 'mixed' if they're mixed." << endl << endl;
  if (caldcs->GetGlobalNumberOfTimeBins() != -1)
    cout<<"Global number of time bins.........................: "<<caldcs->GetGlobalNumberOfTimeBins() << endl;
  if (caldcs->GetGlobalConfigTag() != -1)
    cout<<"Global configuration tag...........................: "<<caldcs->GetGlobalConfigTag() << endl;
  if (caldcs->GetGlobalSingleHitThres() != -1)
    cout<<"Global single hit threshold........................: "<<caldcs->GetGlobalSingleHitThres() << endl;
  if (caldcs->GetGlobalThreePadClustThres() != -1)
    cout<<"Global three pad cluster threshold.................: "<<caldcs->GetGlobalThreePadClustThres()<<endl;
  if (caldcs->GetGlobalSelectiveNoZS() != -1)
    cout<<"Global selective ZS (every i'th event).............: "<<caldcs->GetGlobalSelectiveNoZS() << endl;
  if (caldcs->GetGlobalTCFilterWeight() != -1)
    cout<<"Global tail cancellation filter weight.............: "<<caldcs->GetGlobalTCFilterWeight() << endl;
  if (caldcs->GetGlobalTCFilterShortDecPar() != -1)
    cout<<"Global tail cancellat. filter short decay parameter: "<<caldcs->GetGlobalTCFilterShortDecPar()<<endl;
  if (caldcs->GetGlobalTCFilterLongDecPar() != -1)
    cout<<"Global tail cancellation filt. long decay parameter: "<<caldcs->GetGlobalTCFilterLongDecPar()<<endl;
  if (caldcs->GetGlobalModeFastStatNoise() != -1)
    cout<<"Global fast statistics mode?.......................: "<<caldcs->GetGlobalModeFastStatNoise() << endl;
  if (caldcs->GetGlobalConfigVersion() != "")
    cout<<"Global configuration tag version...................: "<<caldcs->GetGlobalConfigVersion() << endl;
  if (caldcs->GetGlobalConfigName() != "")
    cout<<"Global configuration tag name......................: "<<caldcs->GetGlobalConfigName() << endl;
  if (caldcs->GetGlobalFilterType() != "")
    cout<<"Global filter type.................................: "<<caldcs->GetGlobalFilterType() << endl;
  if (caldcs->GetGlobalReadoutParam() != "")
    cout<<"Global readout parameter...........................: "<<caldcs->GetGlobalReadoutParam() << endl;
  if (caldcs->GetGlobalTestPattern() != "")
    cout<<"Global test pattern................................: "<<caldcs->GetGlobalTestPattern() << endl;
  if (caldcs->GetGlobalTrackletMode() != "")
    cout<<"Global tracklet mode...............................: "<<caldcs->GetGlobalTrackletMode() << endl;
  if (caldcs->GetGlobalTrackletDef() != "")
    cout<<"Global tracklet definition.........................: "<<caldcs->GetGlobalTrackletDef() << endl;
  if (caldcs->GetGlobalTriggerSetup() != "")
    cout<<"Global trigger setup...............................: "<<caldcs->GetGlobalTriggerSetup() << endl;
  if (caldcs->GetGlobalAddOptions() != "")
    cout<<"Global additional options..........................: "<<caldcs->GetGlobalAddOptions() << endl;
  
  cout << endl << "============ Error Summary: ============" << endl;
  if (errors) {
    cout<<"    I noticed some errors, please see above for the specifics." << endl;
  } else {
    cout<<"    I didn't notice any errors, but that doesn't mean there weren't any!" << endl;
  }
  

}
