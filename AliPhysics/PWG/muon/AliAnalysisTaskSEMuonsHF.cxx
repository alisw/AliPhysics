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

/* $Id$ */

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

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliESDHeader.h"
#include "AliAODHeader.h"
#include "AliESDMuonTrack.h"
#include "AliAODTrack.h"
#include "AliCentrality.h"
#include "AliMuonsHFHeader.h"
#include "AliMuonInfoStoreRD.h"
#include "AliMuonInfoStoreMC.h"
#include "AliDimuInfoStoreRD.h"
#include "AliDimuInfoStoreMC.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMuonTrackCuts.h"
#include "AliMuonPairCuts.h"
#include "AliAnalysisTaskSEMuonsHF.h"

ClassImp(AliAnalysisTaskSEMuonsHF)

//_____________________________________________________________________________
AliAnalysisTaskSEMuonsHF::AliAnalysisTaskSEMuonsHF() :
AliAnalysisTaskSE(),
fAnaMode(0),
fIsOutputTree(kFALSE),
fIsMC(kFALSE),
fIsFull(kFALSE),
fCutsMuon(0x0),
fCutsDimu(0x0),
fHeader(0),
fMuonClArr(0),
fDimuClArr(0),
fListOutput(0)
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
AliAnalysisTaskSEMuonsHF::AliAnalysisTaskSEMuonsHF(const char *name, const AliMuonTrackCuts& cutsMuon, const AliMuonPairCuts& cutsDimu) :
AliAnalysisTaskSE(name),
fAnaMode(0),
fIsOutputTree(kFALSE),
fIsMC(kFALSE),
fIsFull(kFALSE),
fCutsMuon(new AliMuonTrackCuts(cutsMuon)),
fCutsDimu(new AliMuonPairCuts(cutsDimu)),
fHeader(0),
fMuonClArr(0),
fDimuClArr(0),
fListOutput(0)
{
  //
  // Constructor
  //
  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskSEMuonsHF::~AliAnalysisTaskSEMuonsHF()
{
  //
  // Default destructor
  //
  if (fCutsMuon)   { delete fCutsMuon;   fCutsMuon  =NULL; }
  if (fCutsDimu)   { delete fCutsDimu;   fCutsDimu  =NULL; }
  if (fHeader)     { delete fHeader;     fHeader    =NULL; }
  if (fMuonClArr)  { delete fMuonClArr;  fMuonClArr =NULL; }
  if (fDimuClArr)  { delete fDimuClArr;  fDimuClArr =NULL; }
  if (fListOutput) { delete fListOutput; fListOutput=NULL; }
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::Init()
{
  // Initialization
  // Setting and initializing the running mode and status

  if (!fHeader) {
    fHeader = new AliMuonsHFHeader();
    fHeader->SetName(AliMuonsHFHeader::StdBranchName());
  } fHeader->SetAnaMode(fAnaMode);
    fHeader->SetIsMC(fIsMC);

  if (!fMuonClArr) {
    if (fIsMC) { 
      fMuonClArr = new TClonesArray("AliMuonInfoStoreMC", 0);
      fMuonClArr->SetName(AliMuonInfoStoreMC::StdBranchName());
    } else {
      fMuonClArr = new TClonesArray("AliMuonInfoStoreRD", 0);
      fMuonClArr->SetName(AliMuonInfoStoreRD::StdBranchName());
    }
  }

  if (fAnaMode!=1 && !fDimuClArr) {
    if (fIsMC) {
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

  if (!fListOutput) fListOutput = new TList();
  fHeader->CreateHistograms(fListOutput);

  if (fIsOutputTree) {
    AddAODBranch("AliMuonsHFHeader", &fHeader);
    AddAODBranch("TClonesArray", &fMuonClArr);
    if (fAnaMode!=1) AddAODBranch("TClonesArray", &fDimuClArr);
  }

  fCutsMuon->Print("mask");
  fCutsDimu->Print("mask");

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::UserExec(Option_t *)
{
  // Execute analysis for current event:
  // muon event header & (di)muon info store

  if (fIsMC) {
    if (fInputHandler->MCEvent()) {
      if (fInputHandler->MCEvent()->GetNumberOfTracks()<=0)
           { AliError("MC event not found. Nothing done!"); return; }
    } else { AliError("MC event not found. Nothing done!"); return; }
  } 

  if ( !fCutsMuon) { AliError("AliMuonTrackCuts should be loaded!"); return;}
  if ((!fCutsDimu) && (fAnaMode!=1)) { AliError("AliMuonPairCuts should be loaded!"); return; }

  Int_t ntrks = 0;
  AliAODEvent *aod = 0;
  AliESDEvent *esd = 0;
  if (((TString)InputEvent()->IsA()->GetName())=="AliAODEvent") {
    aod = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!aod) { AliError("AOD event not found. Nothing done!"); return; }
    if (!fIsMC && (aod->GetHeader()->GetEventType()!=7)) return;
    ntrks = aod->GetNumberOfTracks();
  } else {
    esd = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!esd) { AliError("ESD event not found. Nothing done!"); return; }
    if (!fIsMC && (esd->GetHeader()->GetEventType()!=7)) return;
    ntrks = esd->GetNumberOfMuonTracks();
  } if (fIsOutputTree) AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);

  fHeader->SetEventInfo(fInputHandler);
  fHeader->FillHistosEvnH(fListOutput);

  fMuonClArr->Delete();
  TClonesArray &muonRef = *fMuonClArr;
  Int_t countN = fMuonClArr->GetEntriesFast();

  Int_t trkID[ntrks];
  AliAODTrack        *trkAOD = 0;
  AliESDMuonTrack    *trkESD = 0;
  AliMuonInfoStoreRD *muonRD = 0;
  AliMuonInfoStoreMC *muonMC = 0;
  for (Int_t itrk=0; itrk<ntrks; itrk++) {  // loop over all tracks
    if (aod) {
      trkAOD = (AliAODTrack*)aod->GetTrack(itrk);
      if (!trkAOD->IsMuonTrack())        { trkAOD=0; continue; }
      if (fIsMC) muonMC = new AliMuonInfoStoreMC(trkAOD,fInputHandler->MCEvent(),fCutsMuon->GetSelectionMask(trkAOD),fIsFull);
      else muonRD = new AliMuonInfoStoreRD(trkAOD,fCutsMuon->GetSelectionMask(trkAOD));
      trkAOD = 0;
    } else {
      trkESD = (AliESDMuonTrack*)esd->GetMuonTrack(itrk);
      if (!trkESD->ContainTrackerData()) { trkESD=0; continue; }
      if (fIsMC) muonMC = new AliMuonInfoStoreMC(trkESD,fInputHandler->MCEvent(),fCutsMuon->GetSelectionMask(trkESD),fIsFull);
      else muonRD = new AliMuonInfoStoreRD(trkESD,fCutsMuon->GetSelectionMask(trkESD));
      trkESD = 0;
    } if (muonRD) {
      new(muonRef[countN++]) AliMuonInfoStoreRD(*muonRD);
      if (fAnaMode!=2) fHeader->FillHistosMuon(fListOutput, muonRD);
    } if (muonMC) {
      new(muonRef[countN++]) AliMuonInfoStoreMC(*muonMC);
      if (fAnaMode!=2) fHeader->FillHistosMuon(fListOutput, muonMC, muonMC->Source());
    } trkID[countN-1] = itrk;

    if (muonRD) { delete muonRD; muonRD=0; }
    if (muonMC) { delete muonMC; muonMC=0; }
  }  // end loop of all tracks

  if (fAnaMode==1) { aod=0; esd=0; PostData(1,fListOutput); return; }

  fDimuClArr->Delete();
  countN = fDimuClArr->GetEntriesFast();
  TClonesArray &dimuRef = *fDimuClArr;

  AliDimuInfoStoreRD *dimuRD = 0;
  AliDimuInfoStoreMC *dimuMC = 0;
  ntrks = fMuonClArr->GetEntriesFast();
  for (Int_t itrk=0; itrk<ntrks-1; itrk++) {  // 1st loop over muon tracks
    for (Int_t jtrk=itrk+1; jtrk<ntrks; jtrk++) {  // 2nd loop ofver muon tracks
      UInt_t selMask = 0;
      if (aod) selMask = fCutsDimu->GetSelectionMask(aod->GetTrack(    trkID[itrk]),aod->GetTrack(    trkID[jtrk]));
      else     selMask = fCutsDimu->GetSelectionMask(esd->GetMuonTrack(trkID[itrk]),esd->GetMuonTrack(trkID[jtrk]));

      if (fIsMC) dimuMC = new AliDimuInfoStoreMC((AliMuonInfoStoreMC*)fMuonClArr->At(itrk),(AliMuonInfoStoreMC*)fMuonClArr->At(jtrk),selMask);
      else       dimuRD = new AliDimuInfoStoreRD((AliMuonInfoStoreRD*)fMuonClArr->At(itrk),(AliMuonInfoStoreRD*)fMuonClArr->At(jtrk),selMask);

      if (dimuRD) {
        fHeader->FillHistosDimu(fListOutput, dimuRD);
        if (fIsOutputTree) new(dimuRef[countN++]) AliDimuInfoStoreRD(*dimuRD);
      }
      if (dimuMC) {
        fHeader->FillHistosDimu(fListOutput, dimuMC, dimuMC->Source());
        if (fIsOutputTree) new(dimuRef[countN++]) AliDimuInfoStoreMC(*dimuMC);
      }

      if (dimuRD) { delete dimuRD; dimuRD = 0; }
      if (dimuMC) { delete dimuMC; dimuMC = 0; }
    }  // end 2nd loop of muon tracks
  }  // end 1st loop of muon tracks

  aod = 0; esd = 0;
  PostData(1, fListOutput);

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::Terminate(Option_t *)
{
  // Terminate analysis

  // add the correction matrix

  return;
}

void AliAnalysisTaskSEMuonsHF::NotifyRun()
{
  // Notify of the input handler 

  if (fCutsMuon) fCutsMuon->SetRun(fInputHandler);
  if (fCutsDimu) fCutsDimu->SetRun(fInputHandler);

  return;
}

void AliAnalysisTaskSEMuonsHF::SetEvsHCuts(Double_t cuts[5])  const
{
  AliMuonsHFHeader::SetSelectionCuts(cuts);
}
