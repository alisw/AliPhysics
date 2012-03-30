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

/// \ingroup macros
/// \file MUONAlignment.C
/// \brief Macro for MUON alignment using physics tracks. 
///
/// The macro uses the AliMUONAlignment class to calculate the alignment parameters.
/// An array for the alignment parameters is created and can be filled with
/// initial values that will be used as starting values by the alignment
/// algorithm.
///
/// By default the macro run over galice.root in the working directory. If a file list
/// of galice.root is provided as third argument the macro will run over all files.
/// The macro loop over the files, events and tracks. For each track
/// AliMUONAlignment::ProcessTrack(AliMUONTrack * track) and then
/// AliMUONAlignment::LocalFit(Int_t iTrack, Double_t *lTrackParam, Int_t
/// lSingleFit) are called. After all tracks have been procesed and fitted
/// AliMUONAlignment::GlobalFit(Double_t *parameters,Double_t *errors,Double_t *pulls)
/// is called. The array parameters contains the obatained misalignement parameters.
/// A realigned geometry is generated in a local CDB.
///
/// \author: B. Becker and J. Castillo

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUONAlignment.h"
#include "AliMUONTrack.h"
#include "AliMUONRecoParam.h"
#include "AliMUONTrackParam.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONESDInterface.h"
#include "AliMUONCDB.h"

#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliGeomManager.h"

#include <TString.h>
#include <TError.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <Riostream.h>

#include <fstream>

#endif

void MUONAlignment(Int_t nEvents = 100000, char* geoFilename = "geometry.root", TString esdFileName = "AliESDs.root", TString fileList = "")
{
 
  // Import TGeo geometry (needed by AliMUONTrackExtrap::ExtrapToVertex)
  if ( ! AliGeomManager::GetGeometry() ) {
    AliGeomManager::LoadGeometry(geoFilename);
    if (! AliGeomManager::GetGeometry() ) {
      Error("MUONAlignment", "getting geometry from file %s failed", geoFilename);
      return;
    }
  }
  
  // load necessary data from OCDB
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdbManager->SetRun(0);
  if (!AliMUONCDB::LoadField()) return;
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  
  // reset tracker for restoring initial track parameters at cluster
  AliMUONESDInterface::ResetTracker(recoParam);

  Double_t parameters[4*156];
  Double_t errors[4*156];
  Double_t pulls[4*156];
  for(Int_t k=0;k<4*156;k++) {
    parameters[k]=0.;
    errors[k]=0.;
    pulls[k]=0.;
  }

  Double_t trackParams[8] = {0.,0.,0.,0.,0.,0.,0.,0.};

  // Set initial values here, good guess may help convergence
  // St 1 
  //  Int_t iPar = 0;
  //  parameters[iPar++] =  0.010300 ;  parameters[iPar++] =  0.010600 ;  parameters[iPar++] =  0.000396 ;  

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

  AliMUONGeometryTransformer *transform = new AliMUONGeometryTransformer();
  transform->LoadGeometryData();
  alig->SetGeometryTransformer(transform);

  // Do alignment with magnetic field off
  alig->SetBFieldOn(kFALSE);
  
  // Set tracking station to use
  Bool_t bStOnOff[5] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
  Bool_t bChOnOff[10] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};

  // Set degrees of freedom to align (see AliMUONAlignment)
  alig->AllowVariations(bChOnOff);

  // Fix parameters or add constraints here
//   for (Int_t iSt=0; iSt<5; iSt++)
//     if (!bStOnOff[iSt]) alig->FixStation(iSt+1);
  for (Int_t iCh=0; iCh<10; iCh++)
    if (!bChOnOff[iCh]) alig->FixChamber(iCh+1);

  // Left and right sides of the detector are independent, one can choose to align 
  // only one side
  Bool_t bSpecLROnOff[2] = {kTRUE,kTRUE};
  alig->FixHalfSpectrometer(bChOnOff,bSpecLROnOff);

  alig->SetChOnOff(bChOnOff);
  alig->SetSpecLROnOff(bChOnOff);

  // Here we can fix some detection elements
  alig->FixDetElem(908);
  alig->FixDetElem(1020);

  // Set predifined global constrains: X, Y, P, XvsZ, YvsZ, PvsZ, XvsY, YvsY, PvsY
  Bool_t bVarXYT[9] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
  Bool_t bDetTLBR[4] = {kFALSE,kTRUE,kFALSE,kTRUE};
  //  alig->AddConstraints(bChOnOff,bVarXYT,bDetTLBR,bSpecLROnOff);


  char cFileName[100];  
    
  Int_t lMaxFile = 1000;
  Int_t iFile = 0;
  Int_t iEvent = 0;
  bool bKeepLoop = kTRUE;
  Int_t iTrackTot=0;
  Int_t iTrackOk=0;

  while(bKeepLoop && iFile<lMaxFile){
    iFile++;
    if (bLoop) {
      sFileList.getline(cFileName,100);
      if (sFileList.eof()) bKeepLoop = kFALSE;
    }
    else {
      sprintf(cFileName,esdFileName.Data());
      bKeepLoop = kFALSE;
    }
    if (!strstr(cFileName,"AliESDs")) continue;      
    cout << "Using file: " << cFileName << endl;
    
    // load ESD event
    TFile* esdFile = TFile::Open(cFileName); // open the file
    if (!esdFile || !esdFile->IsOpen()) {
      cout << "opening ESD file " << cFileName << "failed" << endl;
      continue;
    }
    TTree* esdTree = (TTree*) esdFile->Get("esdTree"); // get the tree
    if (!esdTree) {
      cout << "no ESD tree found" << endl;
      esdFile->Close();
      continue;
    }
    AliESDEvent* esdEvent = new AliESDEvent(); // link ESD event to the tree
    esdEvent->ReadFromTree(esdTree);

    Int_t nevents = esdTree->GetEntries();
    cout << "... with " << nevents << endl;
    for(Int_t event = 0; event < nevents; event++) {
      if (iEvent >= nEvents){
	bKeepLoop = kFALSE;
	break;
      }
      iEvent++;

      if (esdTree->GetEvent(event) <= 0) {
	cout << "fails to read ESD object for event " << event << endl;
	continue;
      }

      Int_t nTracks = Int_t(esdEvent->GetNumberOfMuonTracks());
      if (!event%100) cout << " there are " << nTracks << " tracks in event " << event << endl;
      for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
	AliESDMuonTrack* esdTrack = esdEvent->GetMuonTrack(iTrack);
	if (!esdTrack->ContainTrackerData()) continue;
	Double_t invBenMom = esdTrack->GetInverseBendingMomentum();
	fInvBenMom->Fill(invBenMom);
	fBenMom->Fill(1./invBenMom);
	if (TMath::Abs(invBenMom)<=1.04) {
	  AliMUONTrack track;
	  AliMUONESDInterface::ESDToMUON(*esdTrack, track);
	  alig->ProcessTrack(&track);
	  alig->LocalFit(iTrackOk++,trackParams,0);
	}
	iTrackTot++;
      }
    }
    delete esdEvent;
    esdFile->Close();
    cout << "Processed " << iTrackTot << " Tracks so far." << endl;
  }
  alig->GlobalFit(parameters,errors,pulls);

  cout << "Done with GlobalFit " << endl;

  // Store results
  Double_t DEid[156] = {0};
  Double_t MSDEx[156] = {0};
  Double_t MSDEy[156] = {0};
  Double_t MSDEz[156] = {0};
  Double_t MSDEp[156] = {0};
  Double_t DEidErr[156] = {0};
  Double_t MSDExErr[156] = {0};
  Double_t MSDEyErr[156] = {0};
  Double_t MSDEzErr[156] = {0};
  Double_t MSDEpErr[156] = {0};
  Int_t lNDetElem = 4*2+4*2+18*2+26*2+26*2;
  Int_t lNDetElemCh[10] = {4,4,4,4,18,18,26,26,26,26};
  // Int_t lSNDetElemCh[10] = {4,8,12,16,34,52,78,104,130,156};
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
    MSDEx[iDE]=parameters[4*iDE+0];
    MSDEy[iDE]=parameters[4*iDE+1];
    MSDEz[iDE]=parameters[4*iDE+3];
    MSDEp[iDE]=parameters[4*iDE+2];
    MSDExErr[iDE]=(Double_t)alig->GetParError(4*iDE+0);
    MSDEyErr[iDE]=(Double_t)alig->GetParError(4*iDE+1);
    MSDEzErr[iDE]=(Double_t)alig->GetParError(4*iDE+3);
    MSDEpErr[iDE]=(Double_t)alig->GetParError(4*iDE+2);
  }

  cout << "Let's create graphs ...  " << endl;

  TGraphErrors *gMSDEx = new TGraphErrors(lNDetElem,DEid,MSDEx,DEidErr,MSDExErr); 
  TGraphErrors *gMSDEy = new TGraphErrors(lNDetElem,DEid,MSDEy,DEidErr,MSDEyErr); 
  TGraphErrors *gMSDEz = new TGraphErrors(lNDetElem,DEid,MSDEz,DEidErr,MSDEzErr); 
  TGraphErrors *gMSDEp = new TGraphErrors(lNDetElem,DEid,MSDEp,DEidErr,MSDEpErr); 

  cout << "... graphs created, open file ...  " << endl;

  TFile *hFile = new TFile("measShifts.root","RECREATE");

  cout << "... file opened ...  " << endl;

  gMSDEx->Write("gMSDEx");
  gMSDEy->Write("gMSDEy");
  gMSDEz->Write("gMSDEz");
  gMSDEp->Write("gMSDEp");
  fInvBenMom->Write();
  fBenMom->Write();
  hFile->Close();
  
  cout << "... and closed!" << endl;
  // Re Align
  AliMUONGeometryTransformer *newTransform = alig->ReAlign(transform,parameters,true); 
  newTransform->WriteTransformations("transform2ReAlign.dat");
  
  // Generate realigned data in local cdb
  const TClonesArray* array = newTransform->GetMisAlignmentData();

  // 100 mum residual resolution for chamber misalignments?
  alig->SetAlignmentResolution(array,-1,0.01,0.01,0.004,0.003);

  cdbManager->SetSpecificStorage("MUON/Align/Data","local://ReAlignCDB");
  
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON alignment objects with residual misalignment");
  AliCDBId id("MUON/Align/Data", 0, AliCDBRunRange::Infinity()); 
  cdbManager->Put(const_cast<TClonesArray*>(array), id, cdbData);

}

