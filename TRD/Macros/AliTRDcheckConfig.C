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


const Int_t nROC = 540;
const Int_t nROB = 8;
const Int_t nMCM = 18;
const Int_t cArraySize = 1000;



void AnalyzeArray(Int_t states[cArraySize], Int_t occur[cArraySize]) {
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
}



void FillItemInArray(Int_t states[cArraySize], Int_t occur[cArraySize], Int_t item) {
  for (Int_t iArrPos=0; iArrPos<cArraySize; iArrPos++) {
    if (item == -1) break; // value not set
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
    gsmStates[i] = 0;
    gsmOccur[i]  = 0;
    nimStates[i] = 0;
    nimOccur[i]  = 0;
    nevStates[i] = 0;
    nevOccur[i]  = 0;
    nptStates[i] = 0;
    nptOccur[i]  = 0;  
  }
  
  for (Int_t i=0; i<nROC; i++) {
    AliTRDCalDCSFEE *idcsfee = calDCSObj->GetCalDCSFEEObj(i);
    if ((idcsfee == NULL) || (idcsfee->GetStatusBit() != 0)) continue;
    for (Int_t j=0; j<nROB; j++) {
      for (Int_t k=0; k<nMCM; k++) {
	Int_t igsm = idcsfee->GetMCMGlobalState(j,k);
        Int_t inim = idcsfee->GetMCMStateNI(j,k);
        Int_t inev = idcsfee->GetMCMEventCnt(j,k);
	Int_t inpt = idcsfee->GetMCMPtCnt(j,k);
     
        FillItemInArray(gsmStates, gsmOccur, igsm); 
        FillItemInArray(nimStates, nimOccur, inim); 
        FillItemInArray(nevStates, nevOccur, inev); 
        FillItemInArray(nptStates, nptOccur, inpt); 
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


void AliTRDcheckConfig(Int_t runNr=0, TString source="grid", Int_t dumpid=0, char *pathfile="nopathgiven"){

  AliCDBEntry *entry=0;

  // get the source
  if(source.Contains("grid")) {
    cout << "I : Accessing grid storage for run number " << runNr << endl;
    cout << "I : Get CDBManager instance." << endl;
    AliCDBManager *man = AliCDBManager::Instance();
    cout << "I : SetDefaultStorage." << endl;
    man->SetDefaultStorage("alien://folder=/alice/data/2008/LHC08d/OCDB/");
    cout << "I : Get OCDB Entry." << endl;
    entry = man->Get("TRD/Calib/DCSCONFIG", runNr);
  } else if(source.Contains("local")) {
    cout << "I : Accessing local storage" << endl;
    TFile *f = new TFile(pathfile);
    if(f != NULL) {
      entry = (AliCDBEntry*) f->Get("AliCDBEntry");
    }
    else {
      cout << "E : Cannot open file" << endl;
      return;
    }
  } else {
    cout << "E : Please specify a correct source (grid/local)" << endl;
    return;
  }


  Int_t dump[20], contSelect=-1;
  for(Int_t i=0; i<20; i++) dump[i]=0;
  if(!source.Contains("quiet") && !source.Contains("verbose")) {
    cout << "Q : Do you want to dump the AliCDBEntry? (1/0) ";
    cin >> dump[0];
    cout << "Q : Do you want to dump the AliTRDCalDCS Object? (1/0) ";
    cin >> dump[1];
    cout << "Q : Do you want to dump the AliTRDCalDCSFEE Object? (1/0) ";
    cin >> dump[2];
    if(dump[2] == 1) {
      cout << "Q : Which Detector ID? (0..540) ";
      cin >> dumpid;
    }
    cout << "Q : Do you want to get a status bit summary? (1/0) ";
    cin >> dump[3];
    cout << "Q : Do you want to get a status bit histogram? (1/0) ";
    cin >> dump[4];
    cout << "Q : Do you want to get a check on data integrity? (1/0) ";
    cin >> dump[5];
  }
  if(source.Contains("verbose")) for(Int_t i=0; i<20; i++) dump[i]=1;

  // check version
  //TObject *objectCDB = (TObject*)entry->GetObject();
  TObject *objectCDB = (TObject*)entry->GetObject();
  if (objectCDB->IsA()->InheritsFrom("TObjArray")) {
    cout << "I : It seems like you are checking a file containing both SOR and EOR informations." << endl;
    cout << "Q : Do you want to check SOR (0), EOR (1) or their differences (2)? ";
    cin >> contSelect;
    TObjArray *objArrayCDB = (TObjArray*)entry->GetObject();
  }



  // the AliCDBEntry
  if(dump[0] > 0) {
    cout << endl << "============ Dumping the AliCDBEntry ============" << endl;
    entry->Dump();
  }

  
  // the CalDCS object
  AliTRDCalDCS *caldcs;
  switch(contSelect) {
    case -1:
      caldcs = (AliTRDCalDCS*) entry->GetObject();
      break;    
    case 0:
      caldcs = (AliTRDCalDCS*) objArrayCDB->At(0);
      break;    
    case 1:
      caldcs = (AliTRDCalDCS*) objArrayCDB->At(1);
      break;    
    case 2:
      cout << "I : Sorry, not implemented yet. Stop." << endl;
      return;
      break;    
    default:
      cout << "E : Invalid key. Stop." << endl;
      return;
  }

    

  if(dump[1] > 0) {
    cout << endl << "============ Dumping the AliTRDCalDCS Object ============" << endl;
    caldcs->Dump();
  }

  // one CalDCSFEE object
  if(dump[2] > 0) {
    AliTRDCalDCSFEE *dcsfee;
    dcsfee = caldcs->GetCalDCSFEEObj(dumpid);
    cout <<endl<< "============ Dumping the AliTRDCalDCSFEE Object Nr. " << dumpid << " ============" << endl;
    if (dcsfee != NULL) {
      dcsfee->Dump();
    }
  }

  // fill histograms
  TH1F *hStatusBits = new TH1F("hStatusBits", "DCS FEE Status Bits", 10, -.5, 9.5);
  hStatusBits->SetFillColor(38);
  hStatusBits->GetXaxis()->SetTitle("Status Bit Value");
  hStatusBits->GetYaxis()->SetTitle("Occurrence");

  // get a status bit summary
  if(dump[3] > 0) {
    cout << endl << "============ Status Bit Summary: ============" << endl;
    TString bitfivestr = " ROCs with status bit 5. These have been not even responding to any DIM communication attempt. Most probably they just were off.\n    DCS IDs: ";
    Int_t lengthfive = bitfivestr.Length();
    TString bitfourstr = " ROCs with status bit 4! BAD! This might be due to a communication problem between fxsproxy and the feeserver(s) \n    DCS IDs: ";
    Int_t lengthfour = bitfourstr.Length();
    TString bitthreestr = " ROCs with status bit 3!!! VERY BAD! This might be due to a serious communication error between the fxsproxy program and the FEEservers.\n    DCS IDs: ";
    Int_t lengththree = bitthreestr.Length();
    TString bittwostr = " ROCs with status bit 2. These have been in states in which they cannot be read out, e.g. Standby.\n    DCS IDs: ";
    Int_t lengthtwo = bittwostr.Length();
    TString bitonestr = " ROCs with status bit 1! BAD! This means the chamber(s) didn't respont even though is should have been in a good state.\n    DCS IDs: ";
    Int_t lengthone = bitonestr.Length();

    Int_t nSB1=0, nSB2=0, nSB3=0, nSB4=0, nSB5=0;
    for (Int_t i=0; i<540; i++) {
      AliTRDCalDCSFEE *idcsfee;
      idcsfee = caldcs->GetCalDCSFEEObj(i);
      if (idcsfee != NULL) {
	Int_t sb = idcsfee->GetStatusBit();
	if (sb == 5) { bitfivestr  += i; bitfivestr  += "  "; nSB5++; }
	if (sb == 4) { bitfourstr  += i; bitfourstr  += "  "; nSB4++; }
	if (sb == 3) { bitthreestr += i; bitthreestr += "  "; nSB3++; }
	if (sb == 2) { bittwostr   += i; bittwostr   += "  "; nSB2++; }
	if (sb == 1) { bitonestr   += i; bitonestr   += "  "; nSB1++; }
	hStatusBits->Fill(sb);
      }
    }

    if (lengthfive < bitfivestr.Length()) cout << nSB5 << bitfivestr.Data() << endl << endl;
    else cout << "GOOD: No ROCs with status bit 5" << endl;
    if (lengthfour < bitfourstr.Length()) cout << nSB4 << bitfourstr.Data() << endl << endl;
    else cout << "GOOD: No ROCs with status bit 4" << endl;
    if (lengththree < bitthreestr.Length()) cout << nSB3 << bitthreestr.Data() << endl << endl;
    else cout << "GOOD: No ROCs with status bit 3" << endl;
    if (lengthtwo < bittwostr.Length()) cout << nSB2 << bittwostr.Data() << endl << endl;
    else cout << "GOOD: No ROCs with status bit 2" << endl;
    if (lengthone < bitonestr.Length()) cout << nSB1 << bitonestr.Data() << endl << endl;
    else cout << "GOOD: No ROCs with status bit 1" << endl;
  }


  // get a status bit histogram
  if(dump[4] > 0) {
    TCanvas *c1 = new TCanvas("c1");
    gPad->SetLogy();
    hStatusBits->Draw("HIST");
  }



  // get a check on data integrity
  if(dump[5] > 0) {
    cout << endl << "============ Data Integrity: ============" << endl;
    TH1F *tb = new TH1F("tbs", "", 100000, 0, 100000);
    TH1F *ct = new TH1F("cts", "", 100000, 0, 100000);
    const Int_t ngsm = nROB * nMCM;

    Int_t gsm_occ[nROC][ngsm][2];
    Int_t ni_occ[nROC][ngsm][2];
    Int_t ev_occ[nROC][ngsm][2];
    Int_t pt_occ[nROC][ngsm][2];

    for (Int_t i=0; i<nROC; i++) {
      for (Int_t j=0; j<ngsm; j++) {
	gsm_occ[i][j][0] = 0; // values
	gsm_occ[i][j][1] = 0; // counter
	ni_occ[i][j][0] = 0; // values
	ni_occ[i][j][1] = 0; // counter
	ev_occ[i][j][0] = 0; // values
	ev_occ[i][j][1] = 0; // counter
	pt_occ[i][j][0] = 0; // values
	pt_occ[i][j][1] = 0; // counter
      }
    }

    Int_t majGSM[2], majNI[2], majEV[2], majPT[2];
    majGSM[0] = 0; // value
    majGSM[1] = 0; // count
    majNI[0] = 0; // value
    majNI[1] = 0; // count
    majEV[0] = 0; // value
    majEV[1] = 0; // count
    majPT[0] = 0; // value
    majPT[1] = 0; // count


    // find the majority states/counters
    for (Int_t i=0; i<nROC; i++) {
      AliTRDCalDCSFEE *idcsfee;
      idcsfee = caldcs->GetCalDCSFEEObj(i);
      if (idcsfee != NULL) {
	if(idcsfee->GetStatusBit() == 0) {
	  tb->Fill(idcsfee->GetNumberOfTimeBins()-1);
	  ct->Fill(idcsfee->GetConfigTag()-1);
	  for (Int_t j=0; j<nROB; j++) {
	    for (Int_t k=0; k<nMCM; k++) {
	      Int_t igsm = idcsfee->GetMCMGlobalState(j, k);
	      Int_t ini  = idcsfee->GetMCMStateNI(j, k);
	      Int_t iev  = idcsfee->GetMCMEventCnt(j, k);
	      Int_t ipt  = idcsfee->GetMCMPtCnt(j, k);
	      // gsm
	      for (Int_t l=0; l<nROB*nMCM; l++) {
		if (gsm_occ[i][l][1] == 0){
		  gsm_occ[i][l][1] = 1;
		  gsm_occ[i][l][0] = igsm;
		  break;
		}
		else if (gsm_occ[i][l][0] == igsm) {
		  gsm_occ[i][l][1]++;
		  break;
		}
	      }
	      // ni
	      for (Int_t l=0; l<nROB*nMCM; l++) {
		if (ni_occ[i][l][1] == 0){
		  ni_occ[i][l][1] = 1;
		  ni_occ[i][l][0] = ini;
		  break;
		}
		else if (ni_occ[i][l][0] == ini) {
		  ni_occ[i][l][1]++;
		  break;
		}
	      }
	      // ev
	      for (Int_t l=0; l<nROB*nMCM; l++) {
		if (ev_occ[i][l][1] == 0){
		  ev_occ[i][l][1] = 1;
		  ev_occ[i][l][0] = iev;
		  break;
		}
		else if (ev_occ[i][l][0] == iev) {
		  ev_occ[i][l][1]++;
		  break;
		}
	      }
	      // pt
	      for (Int_t l=0; l<nROB*nMCM; l++) {
		if (pt_occ[i][l][1] == 0){
		  pt_occ[i][l][1] = 1;
		  pt_occ[i][l][0] = ipt;
		  break;
		}
		else if (pt_occ[i][l][0] == ipt) {
		  pt_occ[i][l][1]++;
		  break;
		}
	      }
	    }
	  }

	  for (Int_t j=0; j<ngsm; j++) {
	    // gsm
	    if (gsm_occ[i][j][1] > 0) {
	      if (gsm_occ[i][j][1] > majGSM[1]) {
		majGSM[0] = gsm_occ[i][j][0];
		majGSM[1] = gsm_occ[i][j][1];
	      }
	    }
	    // ni
	    if (ni_occ[i][j][1] > 0) {
	      if (ni_occ[i][j][1] > majNI[1]) {
		majNI[0] = ni_occ[i][j][0];
		majNI[1] = ni_occ[i][j][1];
	      }
	    }
	    // ev
	    if (ev_occ[i][j][1] > 0) {
	      if (ev_occ[i][j][1] > majEV[1]) {
		majEV[0] = ev_occ[i][j][0];
		majEV[1] = ev_occ[i][j][1];
	      }
	    }
	    // pt
	    if (pt_occ[i][j][1] > 0) {
	      if (pt_occ[i][j][1] > majPT[1]) {
		majPT[0] = pt_occ[i][j][0];
		majPT[1] = pt_occ[i][j][1];
	      }
	    }  
	  }
	}
      }
    }

    Int_t maxBinCT = ct->GetMaximumBin();
    //cout << "Majority number of configuration tags: " << maxBinCT << endl;
    Int_t maxBinTB = tb->GetMaximumBin();
    //cout << "Majority number of timebins: " << maxBinTB << endl;
    Int_t integrityProblem = 0;

    for (Int_t i=0; i<nROC; i++) {
      AliTRDCalDCSFEE *idcsfee;
      idcsfee = caldcs->GetCalDCSFEEObj(i);
      if (idcsfee != NULL) {
	if(idcsfee->GetStatusBit() == 0) {
	  Int_t ict = idcsfee->GetConfigTag();
	  if(ict != maxBinCT) {
	    cout << "ROB " << i << " has the config tag = " << ict 
	      << " what differs from the majority (=" << maxBinCT << ")!" << endl;
	    integrityProblem++;
	  }
	  Int_t itb = idcsfee->GetNumberOfTimeBins();
	  if(itb != maxBinTB) {
	    cout << "ROB " << i << " has number of time bins = " << itb 
	      << " what differs from the majority (=" << maxBinTB << ")!" << endl;
	    integrityProblem++;
	  }


	  for (Int_t j=0; j<ngsm; j++) {
	    // gsm
	    if ((gsm_occ[i][j][1] > 0) && (gsm_occ[i][j][0] != majGSM[0])) {
	      if (!((gsm_occ[i][j][0] == -1) && ((gsm_occ[i][j][1] == 6) || (gsm_occ[i][j][1] == 40)))) { 
		printf("ROC %3d %3d inconstistent global rstates found with value %2d\n",
		    i,gsm_occ[i][j][1],gsm_occ[i][j][0]);
		integrityProblem++;
	      }
	    }
	    // ni
	    if ((ni_occ[i][j][1] > 0) && (ni_occ[i][j][0] != majNI[0])) {
	      if (!((ni_occ[i][j][0] == -1) && ((ni_occ[i][j][1] == 6) || (ni_occ[i][j][1] == 40)))) { 
		printf("ROC %3d %3d inconstistent network interface states found with value %2d\n",
		    i,ni_occ[i][j][1],ni_occ[i][j][0]);
		integrityProblem++;
	      }
	    }
	    // ev
	    if ((ev_occ[i][j][1] > 0) && (ev_occ[i][j][0] != majEV[0])) {
	      if (!((ev_occ[i][j][0] == -1) && ((ev_occ[i][j][1] == 6) || (ev_occ[i][j][1] == 40)))) { 
		printf("ROC %3d %3d inconstistent event counters found with value %2d\n",
		    i,ev_occ[i][j][1],ev_occ[i][j][0]);
		integrityProblem++;
	      }
	    }
	    // pt
	    if ((pt_occ[i][j][1] > 0) && (pt_occ[i][j][0] != majPT[0])) {
	      if (!((pt_occ[i][j][0] == -1) && ((pt_occ[i][j][1] == 6) || (pt_occ[i][j][1] == 40)))) { 
		printf("ROC %3d %3d inconstistent pretrigger counters found with value %2d\n",
		    i,pt_occ[i][j][1],pt_occ[i][j][0]);
		integrityProblem++;
	      }
	    }

	  }
	}
      }
    }

    if(integrityProblem == 0) cout << "No problem with data integrity found." << endl;
    else cout << endl << "There were in total  " << integrityProblem << "  inconsistencies!" << endl;
  }


  cout << endl << "============ Statistics from RSTATE: ============" << endl;
  cout<<"I : The majority entry is given as well as all other values," << endl;
  cout<<"    sorted according to their occurrence." << endl << endl;
  GetMajoritys(caldcs);
 

  cout << endl << "============ Global Configuraton: ============" << endl;
  cout<<"I : Numbers are -1 if not set and -2 if mixed," << endl;
  cout<<"    strings are empty if not set and 'mixed' if so." << endl << endl;
  cout<<"Global number of time bins.........................: "<<caldcs->GetGlobalNumberOfTimeBins() << endl;
  cout<<"Global configuration tag...........................: "<<caldcs->GetGlobalConfigTag() << endl;
  cout<<"Global single hit threshold........................: "<<caldcs->GetGlobalSingleHitThres() << endl;
  cout<<"Global three pad cluster threshold.................: "<<caldcs->GetGlobalThreePadClustThres()<<endl;
  cout<<"Global selective ZS (every i'th event).............: "<<caldcs->GetGlobalSelectiveNoZS() << endl;
  cout<<"Global tail cancellation filter weight.............: "<<caldcs->GetGlobalTCFilterWeight() << endl;
  cout<<"Global tail cancellat. filter short decay parameter: "<<caldcs->GetGlobalTCFilterShortDecPar()<<endl;
  cout<<"Global tail cancellation filt. long decay parameter: "<<caldcs->GetGlobalTCFilterLongDecPar()<<endl;
  cout<<"Global fast statistics mode?.......................: "<<caldcs->GetGlobalModeFastStatNoise() << endl;
  cout<<"Global configuration tag version...................: "<<caldcs->GetGlobalConfigVersion() << endl;
  cout<<"Global configuration tag name......................: "<<caldcs->GetGlobalConfigName() << endl;
  cout<<"Global filter type.................................: "<<caldcs->GetGlobalFilterType() << endl;
  cout<<"Global readout parameter...........................: "<<caldcs->GetGlobalReadoutParam() << endl;
  cout<<"Global test pattern................................: "<<caldcs->GetGlobalTestPattern() << endl;
  cout<<"Global tracklet mode...............................: "<<caldcs->GetGlobalTrackletMode() << endl;
  cout<<"Global tracklet definition.........................: "<<caldcs->GetGlobalTrackletDef() << endl;
  cout<<"Global trigger setup...............................: "<<caldcs->GetGlobalTriggerSetup() << endl;
  cout<<"Global additional options..........................: "<<caldcs->GetGlobalAddOptions() << endl;

}
