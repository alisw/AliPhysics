/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for the single muon and dimuon from HF analysis,
// using the classes AliMuonsHFHeader,
//                   AliMuonInfoStoreRD,
//                   AliDimuInfoStoreRD,
//                   AliMuonInfoStoreMC,
//                   AliDimuInfoStoreMC.
//
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include <TList.h>
#include <TClonesArray.h>
#include <TFile.h>

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDMuonTrack.h"
#include "AliAODTrack.h"
#include "AliMuonsHFHeader.h"
#include "AliMuonInfoStoreRD.h"
#include "AliMuonInfoStoreMC.h"
#include "AliDimuInfoStoreRD.h"
#include "AliDimuInfoStoreMC.h"
#include "AliAnalysisTaskSEMuonsHF.h"

class AliAnalysisTaskSE;

ClassImp(AliAnalysisTaskSEMuonsHF)

//_____________________________________________________________________________
AliAnalysisTaskSEMuonsHF::AliAnalysisTaskSEMuonsHF() :
AliAnalysisTaskSE(),
fAnaMode(0),
fIsOutputTree(kFALSE),
fIsUseMC(kFALSE),
fHeader(0),
fMuonClArr(0),
fDimuClArr(0),
fListHisHeader(0),
fListHisMuon(),
fListHisDimu(0)
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
AliAnalysisTaskSEMuonsHF::AliAnalysisTaskSEMuonsHF(const char *name) :
AliAnalysisTaskSE(name),
fAnaMode(0),
fIsOutputTree(kFALSE),
fIsUseMC(kFALSE),
fHeader(0),
fMuonClArr(0),
fDimuClArr(0),
fListHisHeader(0),
fListHisMuon(0),
fListHisDimu(0)
{
  //
  // Constructor
  //
}

//_____________________________________________________________________________
AliAnalysisTaskSEMuonsHF::~AliAnalysisTaskSEMuonsHF()
{
  //
  // Default destructor
  //
  if (fHeader)    { delete fHeader;     fHeader    =NULL; }
  if (fMuonClArr) { delete fMuonClArr;  fMuonClArr =NULL; }
  if (fDimuClArr) { delete fDimuClArr;  fDimuClArr =NULL; }

  if (fListHisHeader) { delete    fListHisHeader; fListHisHeader=NULL; }
  if (fListHisMuon)   { delete [] fListHisMuon;   fListHisMuon  =NULL; }
  if (fListHisDimu)   { delete [] fListHisDimu;   fListHisDimu  =NULL; }
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::Init()
{
  // Initialization
  // Setting and initializing the running mode and status

  AliMuonsHFHeader::SetAnaMode(fAnaMode);
  AliMuonsHFHeader::SetIsMC(fIsUseMC);
  if (!fHeader) {
    fHeader = new AliMuonsHFHeader();
    fHeader->SetName(AliMuonsHFHeader::StdBranchName());
  }

  if (!fMuonClArr) {
    if (fIsUseMC) { 
      fMuonClArr = new TClonesArray("AliMuonInfoStoreMC", 0);
      fMuonClArr->SetName(AliMuonInfoStoreMC::StdBranchName());
    } else {
      fMuonClArr = new TClonesArray("AliMuonInfoStoreRD", 0);
      fMuonClArr->SetName(AliMuonInfoStoreRD::StdBranchName());
    }
  }

  if (fAnaMode!=1 && !fDimuClArr) {
    if (fIsUseMC) {
      fDimuClArr = new TClonesArray("AliDimuInfoStoreMC", 0);
      fDimuClArr->SetName(AliDimuInfoStoreMC::StdBranchName());
    } else {
      fDimuClArr = new TClonesArray("AliDimuInfoStoreRD", 0);
      fDimuClArr->SetName(AliDimuInfoStoreRD::StdBranchName());
    }
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::UserCreateOutputObjects()
{
  // Create the output container

  if (!fListHisHeader) fListHisHeader = new TList();
  if (fAnaMode!=2 && !fListHisMuon) {
    if (fIsUseMC) fListHisMuon = new TList[AliMuonInfoStoreMC::NSources()];
    else fListHisMuon = new TList();
  }
  if (fAnaMode!=1 && !fListHisDimu) {
    if (fIsUseMC) fListHisDimu = new TList[AliDimuInfoStoreMC::NSources()];
    else fListHisDimu = new TList();
  }
  fHeader->CreateHistograms(fListHisHeader, fListHisMuon, fListHisDimu);

  if (fIsOutputTree) {
    AddAODBranch("AliMuonsHFHeader", &fHeader);
    AddAODBranch("TClonesArray", &fMuonClArr);
    if (fAnaMode!=1) AddAODBranch("TClonesArray", &fDimuClArr);
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::UserExec(Option_t *)
{
  // Execute analysis for current event:
  // muon event header & (di)muon info store

  AliVEvent *event = dynamic_cast<AliVEvent*>(InputEvent());
  //if (AODEvent() && IsStandardAOD()) event = dynamic_cast<AliVEvent*>(AODEvent());
  TString evName = event->IsA()->GetName();
  Bool_t isAOD = ((evName=="AliAODEvent") ? kTRUE : kFALSE);
  event = 0x0;

  Int_t ntrks = 0;
  AliAODEvent *aod = 0;
  AliESDEvent *esd = 0;
  TClonesArray *mcClArr = 0;
  AliMCEventHandler *mcH = 0;
  if (isAOD) {
    aod = dynamic_cast<AliAODEvent*>(InputEvent());
    //if (AODEvent() && IsStandardAOD()) aod = dynamic_cast<AliAODEvent*>(AODEvent());
    if (!aod) { AliError("AOD event not found. Nothing done!"); return; }
    if (fIsUseMC) {
      mcClArr = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      if (!mcClArr) { AliError("MC Array not found. Nothing done!"); return; }
    }
    fHeader->SetEvent(aod);
    ntrks = aod->GetNTracks();
  } else {
    esd = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!esd) { AliError("ESD event not found. Nothing done!"); return; }
    if (fIsUseMC) {
      mcH = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
      if (!mcH) { AliError("MC Handler not found. Nothing done!"); return; }
    }
    fHeader->SetEvent(esd);
    ntrks = esd->GetNumberOfMuonTracks();
  }

  fHeader->FillHistosEventH(fListHisHeader);

  fMuonClArr->Delete();
  TClonesArray &muonRef = *fMuonClArr;
  Int_t countN = fMuonClArr->GetEntriesFast();

  AliAODTrack        *trkAOD = 0;
  AliESDMuonTrack    *trkESD = 0;
  AliMuonInfoStoreRD *trkRD  = 0;
  AliMuonInfoStoreMC *trkMC  = 0;
  for (Int_t itrk=0; itrk<ntrks; itrk++) {  // loop over all tracks
    if (isAOD) {
      trkAOD = (AliAODTrack*)aod->GetTrack(itrk);
      if (!trkAOD->IsMuonTrack()) { trkAOD=0; trkRD=0; trkMC=0; continue; }
      if (fIsUseMC) trkMC = new AliMuonInfoStoreMC(trkAOD, mcClArr);
      else trkRD = new AliMuonInfoStoreRD(trkAOD);
      trkAOD = 0;
    } else {
      trkESD = (AliESDMuonTrack*)esd->GetMuonTrack(itrk);
      if (!trkESD->ContainTrackerData()) { trkESD=0; trkRD=0; trkMC=0; continue; }
      if (fIsUseMC) trkMC = new AliMuonInfoStoreMC(trkESD, esd, mcH);
      else trkRD = new AliMuonInfoStoreRD(trkESD);
      trkESD = 0;
    }

    if (trkRD) {
      if (fAnaMode!=2) fHeader->FillHistosMuonRD(fListHisMuon, trkRD);
      new(muonRef[countN++]) AliMuonInfoStoreRD(*trkRD);
    }
    if (trkMC) {
      if (fAnaMode!=2) fHeader->FillHistosMuonMC(fListHisMuon, trkMC);
      new(muonRef[countN++]) AliMuonInfoStoreMC(*trkMC);
    }

    if (trkRD) { delete trkRD; trkRD=0; }
    if (trkMC) { delete trkMC; trkMC=0; }
  }  // end loop of all tracks

  if (fAnaMode==1) return;
  fDimuClArr->Delete();
  countN = fDimuClArr->GetEntriesFast();
  TClonesArray &dimuRef = *fDimuClArr;

  AliDimuInfoStoreRD *dimuRD = 0;
  AliDimuInfoStoreMC *dimuMC = 0;
  ntrks = fMuonClArr->GetEntriesFast();
  for (Int_t itrk=0; itrk<ntrks-1; itrk++) {  // 1st loop over muon tracks
    for (Int_t jtrk=itrk+1; jtrk<ntrks; jtrk++) {  // 2nd loop ofver muon tracks
      if (fIsUseMC)
        dimuMC = new AliDimuInfoStoreMC((AliMuonInfoStoreMC*)fMuonClArr->At(itrk), (AliMuonInfoStoreMC*)fMuonClArr->At(jtrk));
      else
        dimuRD = new AliDimuInfoStoreRD((AliMuonInfoStoreRD*)fMuonClArr->At(itrk), (AliMuonInfoStoreRD*)fMuonClArr->At(jtrk));

      if (dimuRD) {
        fHeader->FillHistosDimuRD(fListHisDimu, dimuRD);
        if (fIsOutputTree) new(dimuRef[countN++]) AliDimuInfoStoreRD(*dimuRD);
      }
      if (dimuMC) {
        fHeader->FillHistosDimuMC(fListHisDimu, dimuMC);
        if (fIsOutputTree) new(dimuRef[countN++]) AliDimuInfoStoreMC(*dimuMC);
      }

      if (dimuRD) { delete dimuRD; dimuRD=0; }
      if (dimuMC) { delete dimuMC; dimuMC=0; }
    }  // end 2nd loop of muon tracks
  }  // end 1st loop of muon tracks

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::Terminate(Option_t *)
{
  // Terminate analysis

  TFile *fout = new TFile("muonsHF.root", "RECREATE");
  gDirectory->mkdir("header");
  gDirectory->cd("header");
  fListHisHeader->Write();
  gDirectory->cd("..");

  if (fAnaMode!=2) {
    gDirectory->mkdir("muon");
    gDirectory->cd("muon");
    if (fIsUseMC) {
      char *muonName[] = {"BottomMu", "CharmMu", "PrimaryMu", "SecondaryMu", "NotMu", "Undentified", "All"};
      for (Int_t i=AliMuonInfoStoreMC::NSources(); i--;) {
        gDirectory->mkdir(muonName[i]);
        gDirectory->cd(muonName[i]);
        fListHisMuon[i].Write();
        gDirectory->cd("..");
      }
    } else fListHisMuon->Write();
    gDirectory->cd("..");
  }

  if (fAnaMode!=1) {
    gDirectory->mkdir("dimu");
    gDirectory->cd("dimu");
    if (fIsUseMC) {
      char *dimuName[] = {"BBdiff", "Bchain", "DDdiff", "Dchain", "Resonance", "Background", "All"};
      for (Int_t i=AliDimuInfoStoreMC::NSources(); i--;) {
        gDirectory->mkdir(dimuName[i]);
        gDirectory->cd(dimuName[i]);
        fListHisDimu[i].Write();
        gDirectory->cd("..");
      }
    } else fListHisDimu->Write();
    gDirectory->cd("..");
  }

  fout->Write();
  fout->Close();
  return;
}
