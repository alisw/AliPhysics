//===================================================================================
// This is a macro to analyze TRD/Calib/DCS OCDB objects either
// from the grid for a given run number or from a local object.
// If you want to analyze data from the grid, please don't forget to
// have a valid alien token initialized
//
// Arguments:
// Either provide a run number as the first argument to access the 
// corresponding file on the grid
// or the path + filename as the second argument (and an arbitrary number as the first)
// to access a local file.
//
// Please note that leading zeros in the run number are not supported.
//
// Examples:
// .x AliTRDcheckConfig.C (60111)
// .x AliTRDcheckConfig.C (0, "$ALICE_ROOT/TRD/Calib/DCS/Run0_999999999_v0_s0.root")
//
// Please contact Frederick Kramer or Hans Beck in case of problems
//===================================================================================


const Int_t nROC = 540;
const Int_t nROB = 8;
const Int_t nMCM = 18;
const Int_t cArraySize = 1000;

Bool_t errors = false;
Int_t  calVer = 0;

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

void GetMajoritys(TObject* calDCSObj) {
  
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
  

  Int_t feeArrSiz = 0;
  if (calVer == 1) feeArrSiz = ((AliTRDCalDCS*)calDCSObj)->GetFEEArr()->GetSize();
  if (calVer == 2) feeArrSiz = ((AliTRDCalDCSv2*)calDCSObj)->GetFEEArr()->GetSize();

  for (Int_t i=0; i<nROC && i<feeArrSiz; i++) {
    TObject* idcsfee;
    if (calVer == 1) idcsfee = ((AliTRDCalDCS*)calDCSObj)->GetCalDCSFEEObj(i);
    if (calVer == 2) idcsfee = ((AliTRDCalDCSv2*)calDCSObj)->GetCalDCSFEEObj(i);

    if (idcsfee == NULL) continue;

    Int_t sbit = 0;
    if (calVer == 1) ((AliTRDCalDCSFEE*)idcsfee)->GetStatusBit();
    if (calVer == 2) ((AliTRDCalDCSFEEv2*)idcsfee)->GetStatusBit();

    if (sbit != 0) continue;

    for (Int_t j=0; j<nROB; j++) {
      for (Int_t k=0; k<nMCM; k++) {
	Int_t igsm = 0;
	Int_t inim = 0;
	Int_t inev = 0;
	Int_t inpt = 0;
	if (calVer == 1) {
	  igsm = ((AliTRDCalDCSFEE*)idcsfee)->GetMCMGlobalState(j,k);
	  inim = ((AliTRDCalDCSFEE*)idcsfee)->GetMCMStateNI(j,k);
	  inev = ((AliTRDCalDCSFEE*)idcsfee)->GetMCMEventCnt(j,k);
	  inpt = ((AliTRDCalDCSFEE*)idcsfee)->GetMCMPtCnt(j,k);
	}
	if (calVer == 2) {
	  igsm = ((AliTRDCalDCSFEEv2*)idcsfee)->GetMCMGlobalState(j,k);
	  inim = ((AliTRDCalDCSFEEv2*)idcsfee)->GetMCMStateNI(j,k);
	  inev = ((AliTRDCalDCSFEEv2*)idcsfee)->GetMCMEventCnt(j,k);
	  inpt = ((AliTRDCalDCSFEEv2*)idcsfee)->GetMCMPtCnt(j,k);
	}
	
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



void GetMajorityDifferences(TObject* calDCSObj, TObject* calDCSObj2) {
  
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

  Int_t feeArrSiz1 = 0;
  Int_t feeArrSiz2 = 0;
  if (calVer == 1) {
    feeArrSiz1 = ((AliTRDCalDCS*)calDCSObj)->GetFEEArr()->GetSize();
    feeArrSiz2 = ((AliTRDCalDCS*)calDCSObj2)->GetFEEArr()->GetSize();
  }
  if (calVer == 2) {
    feeArrSiz1 = ((AliTRDCalDCSv2*)calDCSObj)->GetFEEArr()->GetSize();
    feeArrSiz2 = ((AliTRDCalDCSv2*)calDCSObj2)->GetFEEArr()->GetSize();
  }

  for (Int_t i=0; i<nROC && i<feeArrSiz1 && i<feeArrSiz2; i++) {
    TObject* idcsfee;
    TObject* idcsfee2;

    if (calVer == 1) {
      idcsfee  = ((AliTRDCalDCS*)calDCSObj)->GetCalDCSFEEObj(i);
      idcsfee2 = ((AliTRDCalDCS*)calDCSObj2)->GetCalDCSFEEObj(i);
    }
    if (calVer == 2) {
      idcsfee  = ((AliTRDCalDCSv2*)calDCSObj)->GetCalDCSFEEObj(i);
      idcsfee2 = ((AliTRDCalDCSv2*)calDCSObj2)->GetCalDCSFEEObj(i);
    }
    if ((idcsfee == NULL) || (idcsfee2 == NULL)) continue;

    Int_t sbit = 0;
    if (calVer == 1) sbit = ((AliTRDCalDCSFEE*)idcsfee)->GetStatusBit();
    if (calVer == 2) sbit = ((AliTRDCalDCSFEEv2*)idcsfee)->GetStatusBit();
    if (sbit != 0) continue;

    for (Int_t j=0; j<nROB; j++) {
      for (Int_t k=0; k<nMCM; k++) {
	Int_t igsm, inim, inev, inpt, igsm1, inim1, inev1, inpt1, igsm2, inim2, inev2, inpt2;
	if (calVer == 1) {
	  igsm1 = ((AliTRDCalDCSFEE*)idcsfee)->GetMCMGlobalState(j,k);
	  inim1 = ((AliTRDCalDCSFEE*)idcsfee)->GetMCMStateNI(j,k);
	  inev1 = ((AliTRDCalDCSFEE*)idcsfee)->GetMCMEventCnt(j,k);
	  inpt1 = ((AliTRDCalDCSFEE*)idcsfee)->GetMCMPtCnt(j,k);
	  igsm2 = ((AliTRDCalDCSFEE*)idcsfee2)->GetMCMGlobalState(j,k);
	  inim2 = ((AliTRDCalDCSFEE*)idcsfee2)->GetMCMStateNI(j,k);
	  inev2 = ((AliTRDCalDCSFEE*)idcsfee2)->GetMCMEventCnt(j,k);
	  inpt2 = ((AliTRDCalDCSFEE*)idcsfee2)->GetMCMPtCnt(j,k);
	}
	if (calVer == 2) {
	  igsm1 = ((AliTRDCalDCSFEEv2*)idcsfee)->GetMCMGlobalState(j,k);
	  inim1 = ((AliTRDCalDCSFEEv2*)idcsfee)->GetMCMStateNI(j,k);
	  inev1 = ((AliTRDCalDCSFEEv2*)idcsfee)->GetMCMEventCnt(j,k);
	  inpt1 = ((AliTRDCalDCSFEEv2*)idcsfee)->GetMCMPtCnt(j,k);
	  igsm2 = ((AliTRDCalDCSFEEv2*)idcsfee2)->GetMCMGlobalState(j,k);
	  inim2 = ((AliTRDCalDCSFEEv2*)idcsfee2)->GetMCMStateNI(j,k);
	  inev2 = ((AliTRDCalDCSFEEv2*)idcsfee2)->GetMCMEventCnt(j,k);
	  inpt2 = ((AliTRDCalDCSFEEv2*)idcsfee2)->GetMCMPtCnt(j,k);
	}

	igsm = igsm1 - igsm2;
	inim = inim1 - inim2;
	inev = inev2 - inev1;
	inpt = inpt2 - inpt1;
	
	// if they were set to -1, it means they were not actauuly set
	// change -1 to -100 to mean they werent set since the above 
	// can give negitives
	if (igsm1 == -1 && igsm == 0) igsm =-100000;
	if (inim1 == -1 && inim == 0) inim =-100000;
	if (inev1 == -1 && inev == 0) inev =-100000;
	if (inpt1 == -1 && inpt == 0) inpt =-100000;
	
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
    man->SetDefaultStorageFromRun(runNr);

    cout << "I : Get OCDB Entry." << endl;
    entry = man->Get("TRD/Calib/DCS", runNr);
    if (entry == NULL) {
      cout << endl << "ERROR: Unable to get the AliTRDCalDCS object from the OCDB for run number " << runNr << "." << endl;
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
  
  Int_t iesor=0;
  for (iesor=0; iesor<3; iesor++) if(objArrayCDB->At(iesor)) break;
  if (iesor > 1) {
    cout << "E : Neither the start or end of run objects were in the root file.";
    return;
  }

  Bool_t hasSOR = (objArrayCDB->At(0));
  Bool_t hasEOR = (objArrayCDB->At(1));
  printf("SOR entry: %d, EOR entry: %d\n", hasSOR, hasEOR);

  if (!strcmp(objArrayCDB->At(iesor)->ClassName(),"AliTRDCalDCS"))   calVer = 1;
  else if (!strcmp(objArrayCDB->At(iesor)->ClassName(),"AliTRDCalDCSv2")) calVer = 2;
  else {  
    cout << "E : Object types undefined.";
    return;
  }

  Bool_t sorandeor = true;
  TObject *caldcs  = objArrayCDB->At(0);
  TObject *caldcs2 = objArrayCDB->At(1);

  if (caldcs == NULL && caldcs2 == NULL) {
    cout << "E : Neither the start or end of run objects were in the root file.";
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

  Int_t feeArrSiz = 0;
  if (calVer == 1) feeArrSiz = ((AliTRDCalDCS*)caldcs)->GetFEEArr()->GetSize();
  if (calVer == 2) feeArrSiz = ((AliTRDCalDCSv2*)caldcs)->GetFEEArr()->GetSize();

  Int_t nSB1=0, nSB2=0, nSB3=0, nSB4=0, nSB5=0, nTot=0;
  for (Int_t i=0; i<nROC && i<feeArrSiz; i++) {
    TObject* idcsfee;
    if (calVer == 1) idcsfee = ((AliTRDCalDCS*)caldcs)->GetCalDCSFEEObj(i);
    if (calVer == 2) idcsfee = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(i);
    if (idcsfee != NULL) {
      Int_t sb;
      if (calVer == 1) sb = ((AliTRDCalDCSFEE*)idcsfee)->GetStatusBit();
      if (calVer == 2) sb = ((AliTRDCalDCSFEEv2*)idcsfee)->GetStatusBit();
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
  for (Int_t i=0; i<nROC && i<feeArrSiz; i++) {
    TObject* idcsfee;
    TObject* idcsfee2;
    if (calVer == 1) {
      if (caldcs)  idcsfee  = ((AliTRDCalDCS*)caldcs)->GetCalDCSFEEObj(i);
      if (caldcs2) idcsfee2 = ((AliTRDCalDCS*)caldcs2)->GetCalDCSFEEObj(i);
    }
    if (calVer == 2) {
      if (caldcs)  idcsfee  = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(i);
      if (caldcs2) idcsfee2 = ((AliTRDCalDCSv2*)caldcs2)->GetCalDCSFEEObj(i);
    }
    if (idcsfee != NULL && idcsfee2 != NULL) {
      Int_t sbd1 = 0;
      Int_t sbd2 = 0;
      if (calVer == 1) {
	sbd1 = ((AliTRDCalDCSFEE*)idcsfee)->GetStatusBit();
	sbd2 = ((AliTRDCalDCSFEE*)idcsfee2)->GetStatusBit();
      }
      if (calVer == 2) {
	sbd1 = ((AliTRDCalDCSFEEv2*)idcsfee)->GetStatusBit();
	sbd2 = ((AliTRDCalDCSFEEv2*)idcsfee2)->GetStatusBit();
      }
      Int_t sbd = sbd1 - sbd2;
      if (sbd != 0) { 
	cout << "ROC " << i << " changed from state " << sbd1 << " at start of the run to "  << sbd2 << " at the end of the run." << endl;
	cout << "ROC " << i << " changed from state " << sbd1 << " at start of the run to "  << sbd2 << " at the end of the run." << endl;
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
 
  Int_t   gtb, gct, gsh, gtc, gsz, gfw, gfs, gfl, gsn;
  TString gcv, gcn, gft, grp, gtp, gtm, gtd, gts, gao;

  if (calVer == 1) {
    gtb = ((AliTRDCalDCS*)caldcs)->GetGlobalNumberOfTimeBins();
    gct = ((AliTRDCalDCS*)caldcs)->GetGlobalConfigTag();
    gsh = ((AliTRDCalDCS*)caldcs)->GetGlobalSingleHitThres();
    gtc = ((AliTRDCalDCS*)caldcs)->GetGlobalThreePadClustThres();
    gsz = ((AliTRDCalDCS*)caldcs)->GetGlobalSelectiveNoZS();
    gfw = ((AliTRDCalDCS*)caldcs)->GetGlobalTCFilterWeight();
    gfs = ((AliTRDCalDCS*)caldcs)->GetGlobalTCFilterShortDecPar();
    gfl = ((AliTRDCalDCS*)caldcs)->GetGlobalTCFilterLongDecPar();
    gsn = ((AliTRDCalDCS*)caldcs)->GetGlobalModeFastStatNoise();
    gcv = ((AliTRDCalDCS*)caldcs)->GetGlobalConfigVersion();
    gcn = ((AliTRDCalDCS*)caldcs)->GetGlobalConfigName();
    gft = ((AliTRDCalDCS*)caldcs)->GetGlobalFilterType();
    grp = ((AliTRDCalDCS*)caldcs)->GetGlobalReadoutParam();
    gtp = ((AliTRDCalDCS*)caldcs)->GetGlobalTestPattern();
    gtm = ((AliTRDCalDCS*)caldcs)->GetGlobalTrackletMode();
    gtd = ((AliTRDCalDCS*)caldcs)->GetGlobalTrackletDef();
    gts = ((AliTRDCalDCS*)caldcs)->GetGlobalTriggerSetup();
    gao = ((AliTRDCalDCS*)caldcs)->GetGlobalAddOptions();
  }
  if (calVer == 2) {
    gtb = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetNumberOfTimeBins();
    gct = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetConfigTag();
    gsh = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetSingleHitThres();
    gtc = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetThreePadClustThres();
    gsz = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetSelectiveNoZS();
    gfw = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetTCFilterWeight();
    gfs = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetTCFilterShortDecPar();
    gfl = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetTCFilterLongDecPar();
    gcv = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetConfigVersion();
    gcn = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetConfigName();
    gft = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetFilterType();
    grp = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetReadoutParam();
    gtp = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetTestPattern();
    gtm = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetTrackletMode();
    gtd = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetTrackletDef();
    gts = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetTriggerSetup();
    gao = ((AliTRDCalDCSv2*)caldcs)->GetCalDCSFEEObj(0)->GetAddOptions();
  }


  if (gtb != -1) cout<<"Global number of time bins.........................: "<< gtb << endl;
  if (gct != -1) cout<<"Global configuration tag...........................: "<< gct << endl;
  if (gsh != -1) cout<<"Global single hit threshold........................: "<< gsh << endl;
  if (gtc != -1) cout<<"Global three pad cluster threshold.................: "<< gtc << endl;
  if (gsz != -1) cout<<"Global selective ZS (every i'th event).............: "<< gsz << endl;
  if (gfw != -1) cout<<"Global tail cancellation filter weight.............: "<< gfs << endl;
  if (gfs != -1) cout<<"Global tail cancellat. filter short decay parameter: "<< gfs << endl;
  if (gfl != -1) cout<<"Global tail cancellation filt. long decay parameter: "<< gfl << endl;
  if (gsn != -1) cout<<"Global fast statistics mode?.......................: "<< gsn << endl;
  if (gcv != "") cout<<"Global configuration tag version...................: "<< gcv << endl;
  if (gcn != "") cout<<"Global configuration tag name......................: "<< gcn << endl;
  if (gft != "") cout<<"Global filter type.................................: "<< gft << endl;
  if (grp != "") cout<<"Global readout parameter...........................: "<< grp << endl;
  if (gtp != "") cout<<"Global test pattern................................: "<< gtp << endl;
  if (gtm != "") cout<<"Global tracklet mode...............................: "<< gtm << endl;
  if (gtd != "") cout<<"Global tracklet definition.........................: "<< gtd << endl;
  if (gts != "") cout<<"Global trigger setup...............................: "<< gts << endl;
  if (gao != "") cout<<"Global additional options..........................: "<< gao << endl;
  
  cout << endl << "============ Error Summary: ============" << endl;
  if (errors) {
    cout<<"    I noticed some errors, please see above for the specifics." << endl;
  } else {
    cout<<"    I didn't notice any errors, but that doesn't mean there weren't any!" << endl;
  }
  

}
