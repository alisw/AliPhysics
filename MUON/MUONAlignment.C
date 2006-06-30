/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

// ---
// Macro for MUON alignment using physics tracks. The macro uses AliMUONAlignment
// class to calculate the alignment parameters.
// An array for the alignment parameters is created and can be filled with
// initial values that will be used as starting values by the alignment
// algorithm.
// Ny default the macro run over galice.root in the working directory. If a file list
// of galice.root is provided as third argument the macro will run over all files.
// The macro loop over the files, events and tracks. For each track
// AliMUONAlignment::ProcessTrack(AliMUONTrack * track) and then
// AliMUONAlignment::LocalFit(Int_t iTrack, Double_t *lTrackParam, Int_t
// lSingleFit) are called. After all tracks have been procesed and fitted
// AliMUONAlignment::GlobalFit(Double_t *parameters,Double_t *errors,Double_t *pulls)
// is called. The array parameters contains the obatained misalignement parameters.
// A realigned geometry is generated in a local CDB.
//
// Authors: B. Becker and J. Castillo
// ---

void MUONAlignment(Int_t nEvents = 100000, TString fileName = "galice.root", TString fileList = "")
{
 
  Double_t parameters[3*156];
  Double_t errors[3*156];
  Double_t pulls[3*156];
  for(Int_t k=0;k<3*156;k++) {
    parameters[k]=0.;
    errors[k]=0.;
    pulls[k]=0.;
  }

  Double_t trackParams[8] = {0.,0.,0.,0.,0.,0.,0.,0.};

  Int_t iPar = 0;
  // Set initial values here, good guess may help convergence
  // St 1 
  //   parameters[iPar++] =  0.010300 ;  parameters[iPar++] =  0.010600 ;  parameters[iPar++] =  0.000396 ;  

  bool bLoop = kFALSE;
  ifstream sFileList;
  if (fileList.Contains(".list")) {
    cout << "Reading file list: " << fileList.Data() << endl;
    bLoop = kTRUE;

    TString fullListName(".");
    fullListName +="/";
    fullListName +=fileList;
    sFileList.open(fileList.Data());
  }
  
  TH1F *fInvBenMom = new TH1F("fInvBenMom","fInvBenMom",200,-0.1,0.1); 
  TH1F *fBenMom = new TH1F("fBenMom","fBenMom",200,-40,40); 

  AliMUONAlignment* alig = new AliMUONAlignment();
  alig->InitGlobalParameters(parameters);

  AliMUONGeometryTransformer *transform = new AliMUONGeometryTransformer(true);
  transform->ReadGeometryData("volpath.dat", "transform.dat");
  alig->SetGeometryTransformer(transform);

  char cFileName[100];  
  AliMUONDataInterface amdi;

  Int_t lMaxFile = 1000;
  Int_t iFile = 0;
  bool bKeepLoop = kTRUE;
  while(bKeepLoop && iFile<lMaxFile){
    iFile++;
    if (bLoop) {
      sFileList.getline(cFileName,100);
      if (sFileList.eof()) bKeepLoop = kFALSE;
    }
    else {
      sprintf(cFileName,fileName.Data());
      bKeepLoop = kFALSE;
    }
    if (!strstr(cFileName,"galice.root")) continue;      
    cout << "Using file: " << cFileName << endl;
    amdi.SetFile(cFileName);
    Int_t nevents = amdi.NumberOfEvents();
    for(Int_t event = 0; event < nevents; event++) {
      amdi.GetEvent(event);
      Int_t ntracks = amdi.NumberOfRecTracks();
      cout << " there are " << ntracks << " tracks in event " << event << endl;
      Int_t iTrack=0;
      Int_t iTrackOk=0;
      AliMUONTrack* track = amdi.RecTrack(iTrack);
      while(track) {   
	Double_t invBenMom = track->GetTrackParamAtVertex()->GetInverseBendingMomentum();
	fInvBenMom->Fill(invBenMom);
	fBenMom->Fill(1./invBenMom);
	if (TMath::Abs(invBenMom)<=1.04) {
	  cout << "Track " << iTrack << endl;
	  alig->ProcessTrack(track);
	  cout << "Calling LocalFit" << endl;
	  alig->LocalFit(iTrackOk++,trackParams,0);
	}
	iTrack++;
	track = amdi.RecTrack(iTrack);
      }
    }
  }
  alig->GlobalFit(parameters,errors,pulls);

  cout << "Done with GlobalFit " << endl;

  // Store results
  Double_t DEid[156] = {0};
  Double_t MSDEx[156] = {0};
  Double_t MSDEy[156] = {0};
  Double_t MSDExt[156] = {0};
  Double_t MSDEyt[156] = {0};
  Double_t DEidErr[156] = {0};
  Double_t MSDExErr[156] = {0};
  Double_t MSDEyErr[156] = {0};
  Double_t MSDExtErr[156] = {0};
  Double_t MSDEytErr[156] = {0};
  Int_t lNDetElem = 4*2+4*2+18*2+26*2+26*2;
  Int_t lNDetElemCh[10] = {4,4,4,4,18,18,26,26,26,26};
  Int_t lSNDetElemCh[10] = {4,8,12,16,34,52,78,104,130,156};
  Int_t idOffset = 0; // 400
  Int_t lSDetElemCh = 0;
  for(Int_t iDE=0; iDE<lNDetElem; iDE++){
    DEidErr[iDE] = 0.;
    DEid[iDE] = idOffset+100;
    DEid[iDE] += iDE; 
    lSDetElemCh = 0;
    for(Int_t iCh=0; iCh<9; iCh++){
      lSDetElemCh += lNDetElemCh[iCh];
      if (iDE>=lSDetElemCh) {
	DEid[iDE] += 100;
	DEid[iDE] -= lNDetElemCh[iCh];
      }
    }
    MSDEx[iDE]=parameters[3*iDE+0];
    MSDEy[iDE]=parameters[3*iDE+1];
    MSDExt[iDE]=parameters[3*iDE+2];
    MSDEyt[iDE]=parameters[3*iDE+2];
    MSDExErr[iDE]=(Double_t)alig->GetParError(3*iDE+0);
    MSDEyErr[iDE]=(Double_t)alig->GetParError(3*iDE+1);
    MSDExtErr[iDE]=(Double_t)alig->GetParError(3*iDE+2);
    MSDEytErr[iDE]=(Double_t)alig->GetParError(3*iDE+2);
  }

  cout << "Let's create graphs ...  " << endl;

  TGraphErrors *gMSDEx = new TGraphErrors(lNDetElem,DEid,MSDEx,DEidErr,MSDExErr); 
  TGraphErrors *gMSDEy = new TGraphErrors(lNDetElem,DEid,MSDEy,DEidErr,MSDEyErr); 
  TGraphErrors *gMSDExt = new TGraphErrors(lNDetElem,DEid,MSDExt,DEidErr,MSDExtErr); 
  TGraphErrors *gMSDEyt = new TGraphErrors(lNDetElem,DEid,MSDEyt,DEidErr,MSDEytErr); 

  cout << "... graphs created, open file ...  " << endl;

  TFile *hFile = new TFile("measShifts.root","RECREATE");

  cout << "... file opened ...  " << endl;

  gMSDEx->Write("gMSDEx");
  gMSDEy->Write("gMSDEy");
  gMSDExt->Write("gMSDExt");
  gMSDEyt->Write("gMSDEyt");
  fInvBenMom->Write();
  fBenMom->Write();
  hFile->Close();
  
  cout << "... and closed!" << endl;
  // Re Align
  AliMUONGeometryTransformer *newTransform = alig->ReAlign(transform,parameters,true); 
  newTransform->WriteTransformations("transform2ReAlign.dat");
  
  // Generate realigned data in local cdb
  TClonesArray* array = newTransform->GetMisAlignmentData();
   
  // CDB manager
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("local://ReAlignCDB");
  
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON alignment objects with residual misalignment");
  AliCDBId id("MUON/Align/Data", 0, 0); 
  cdbManager->Put(array, id, cdbData);

} 
