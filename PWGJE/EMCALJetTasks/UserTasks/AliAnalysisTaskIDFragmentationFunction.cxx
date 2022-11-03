// *************************************************************************
// *                                                                       *
// * Task for ID Fragmentation Function Analysis                           *
// *                                                                       *
// *************************************************************************


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

/* $Id: */

#include <vector>

#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TFile.h"
#include "TKey.h"
#include "TRandom3.h"
#include "TAxis.h"
#include "TVector3.h"
#include "TDatime.h"

#include "AliAODInputHandler.h" 
#include "AliAODHandler.h" 
#include "AliESDEvent.h"
#include "AliAODMCParticle.h"
#include "AliEmcalJet.h"
#include "AliAODJetEventBackground.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"

#include "AliAnalysisHelperJetTasks.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliVParticle.h"
#include "AliVEvent.h"

#include "FJ_includes.h"
#include "AliFJWrapper.h"

#include "AliAnalysisTaskMTFPID.h"
#include "AliPIDResponse.h"

#include "AliESDtrackCuts.h"

#include "AliPieceWisePoly.h"
#include "AliHelperClassFastSimulation.h"

#include "AliAnalysisTaskIDFragmentationFunction.h"

using std::cout;
using std::endl;
using std::cerr;

ClassImp(AliAnalysisTaskIDFragmentationFunction)

//____________________________________________________________________________
AliAnalysisTaskIDFragmentationFunction::AliAnalysisTaskIDFragmentationFunction()
   : AliAnalysisTaskEmcalJet("AliAnalysisTaskIDFragmentationFunction", kTRUE)
   ,fESD(0)
   ,fAOD(0)
   ,fAODJets(0)  
   ,fAODExtension(0)
   ,fNonStdFile("")
   ,fMinMultiplicity(-1)
   ,fMaxMultiplicity(-1)
   ,fNameTrackContainer("Tracks")
   ,fNameTrackContainerEfficiency("")
   ,fNameMCParticleContainer("")
   ,fNameJetContainer("")
   ,fNameMCParticleJetContainer("")
   ,fUseAODInputJets(kTRUE)
   ,fUsePhysicsSelection(kTRUE)
   ,fEvtSelectionMask(0)
   ,fEventClass(-1)
   ,fMaxVertexZ(10)
   ,fFFRadius(0)
   ,fFFMinLTrackPt(-1)
   ,fFFMaxTrackPt(-1)
   ,fFFMinnTracks(0)   
   ,fAvgTrials(0)
   ,fStudyTransversalJetStructure(kFALSE)
   ,fCommonHistList(0)
   ,fh1EvtSelection(0)
   ,fh1VtxSelection(0)
   ,fh1VertexNContributors(0)
   ,fh1VertexZ(0)
   ,fh1EvtMult(0)
   ,fh1EvtCent(0)
   ,fh1Xsec(0)
   ,fh1Trials(0)
   ,fh1PtHard(0)
   ,fh1PtHardTrials(0)
   ,fh1EvtsPtHardCut(0)
   ,fh1nRecJetsCuts(0)
   ,fh1nRecJetsCutsGroomed(0x0)
   ,fh1nRecJetsCuts2(0)
   ,fh1nRCinUnderground(0)
   ,fh1nGenJets(0)
   ,fh1TotJetEnergy(0)
   ,fhDCA_XY(0)
   ,fhDCA_Z(0)
   ,fhJetPtRefMultEta5(0)
   ,fhJetPtRefMultEta8(0)
   ,fhJetPtMultPercent(0)
   ,fh3trackDensity(0x0)
   ,fh2TrackDef(0x0)
   ,fRandom(0)
   ,fOnlyLeadingJets(kFALSE)
   ,fMCPtHardCut(-1.)
   ,fAnaUtils(0)
   // PID framework
   ,fNumInclusivePIDtasks(0)
   ,fNumJetPIDtasks(0)
   ,fNumJetUEPIDtasks(0)
   ,fNameInclusivePIDtask(0x0)
   ,fNameJetPIDtask(0x0)
   ,fNameJetUEPIDtask(0x0)
   ,fInclusivePIDtask(0x0)
   ,fJetPIDtask(0x0)
   ,fJetUEPIDtask(0x0)
   ,fUseInclusivePIDtask(kFALSE)
   ,fUseJetPIDtask(kFALSE)
   ,fUseJetUEPIDtask(kFALSE)
   ,fIsPP(kFALSE)
   ,fFillDCA(kFALSE)
   ,fDoGroomedJets(kFALSE)
   ,fBetaSoftDrop(0.0)
   ,fZSoftDrop(0.1)
   ,fUseFastSimulations(kFALSE)
   ,fastSimulationParameters("")
   ,fEffFunctions(0x0)
   ,fFastSimEffFactor(1.0)
   ,fFastSimRes(0.002)
   ,fFastSimResFactor(1.0)
   ,fFFChange(AliAnalysisTaskIDFragmentationFunction::kNoChange)
   ,fRCTrials(1)
   ,fUEMethods(0x0)
   ,fUseRealJetArea(kTRUE)
{
   // default constructor
  
  if (fFillDCA) {
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      fhDCA_XY_prim_MCID[i] = 0x0;
      fhDCA_Z_prim_MCID[i] = 0x0;
      
      fhDCA_XY_sec_MCID[i] = 0x0;
      fhDCA_Z_sec_MCID[i] = 0x0;
    }
  }
}

//_______________________________________________________________________________________________
AliAnalysisTaskIDFragmentationFunction::AliAnalysisTaskIDFragmentationFunction(const char *name) 
  : AliAnalysisTaskEmcalJet(name, kTRUE)
  ,fESD(0)
  ,fAOD(0)
  ,fAODJets(0)  
  ,fAODExtension(0)
  ,fNonStdFile("")
  ,fMinMultiplicity(-1)
  ,fMaxMultiplicity(-1)
  ,fNameTrackContainer("Tracks")
  ,fNameTrackContainerEfficiency("")
  ,fNameMCParticleContainer("")
  ,fNameJetContainer("")
  ,fNameMCParticleJetContainer("")  
  ,fUseAODInputJets(kTRUE)
  ,fUsePhysicsSelection(kTRUE)
  ,fEvtSelectionMask(0)
  ,fEventClass(-1)
  ,fMaxVertexZ(10)
  ,fFFRadius(0)
  ,fFFMinLTrackPt(-1)
  ,fFFMaxTrackPt(-1)
  ,fFFMinnTracks(0)  
  ,fAvgTrials(0)
  ,fStudyTransversalJetStructure(kFALSE)
  ,fCommonHistList(0)
  ,fh1EvtSelection(0)
  ,fh1VtxSelection(0)
  ,fh1VertexNContributors(0)
  ,fh1VertexZ(0)
  ,fh1EvtMult(0)
  ,fh1EvtCent(0)
  ,fh1Xsec(0)
  ,fh1Trials(0)
  ,fh1PtHard(0)
  ,fh1PtHardTrials(0)
  ,fh1EvtsPtHardCut(0)
  ,fh1nRecJetsCuts(0)
  ,fh1nRecJetsCutsGroomed(0x0)
  ,fh1nRecJetsCuts2(0)
  ,fh1nRCinUnderground(0)
  ,fh1nGenJets(0)
  ,fh1TotJetEnergy(0x0)
  ,fhDCA_XY(0)
  ,fhDCA_Z(0)
  ,fhJetPtRefMultEta5(0)
  ,fhJetPtRefMultEta8(0)
  ,fhJetPtMultPercent(0)
  ,fh3trackDensity(0x0)
  ,fh2TrackDef(0x0)
  ,fRandom(0)
  ,fOnlyLeadingJets(kFALSE)
  ,fMCPtHardCut(-1.)
  ,fAnaUtils(0)
  // PID framework
  ,fNumInclusivePIDtasks(0)
  ,fNumJetPIDtasks(0)
  ,fNumJetUEPIDtasks(0)
  ,fNameInclusivePIDtask(0x0)
  ,fNameJetPIDtask(0x0)
  ,fNameJetUEPIDtask(0x0)
  ,fInclusivePIDtask(0x0)
  ,fJetPIDtask(0x0)
  ,fJetUEPIDtask(0x0)
  ,fUseInclusivePIDtask(kFALSE)
  ,fUseJetPIDtask(kFALSE)
  ,fUseJetUEPIDtask(kFALSE)
  ,fIsPP(kFALSE)
  ,fFillDCA(kFALSE)
  ,fDoGroomedJets(kFALSE)
  ,fBetaSoftDrop(0.0)
  ,fZSoftDrop(0.1)
  ,fUseFastSimulations(kFALSE) 
  ,fastSimulationParameters("")
  ,fEffFunctions(0x0)
  ,fFastSimEffFactor(1.0)
  ,fFastSimRes(0.002)
  ,fFastSimResFactor(1.0)  
  ,fFFChange(AliAnalysisTaskIDFragmentationFunction::kNoChange)
  ,fRCTrials(1)
  ,fUEMethods(0x0)
  ,fUseRealJetArea(kTRUE)
{
  // constructor
  
  if (fFillDCA) {
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      fhDCA_XY_prim_MCID[i] = 0x0;
      fhDCA_Z_prim_MCID[i] = 0x0;
      
      fhDCA_XY_sec_MCID[i] = 0x0;
      fhDCA_Z_sec_MCID[i] = 0x0;
    }
  }
  
  DefineOutput(1,TList::Class());
}

//___________________________________________________________________________
AliAnalysisTaskIDFragmentationFunction::~AliAnalysisTaskIDFragmentationFunction()
{
  // destructor

  delete fRandom;
  
  delete [] fNameInclusivePIDtask;
  fNameInclusivePIDtask = 0x0;
    
  delete [] fInclusivePIDtask;
  fInclusivePIDtask = 0x0;
  
  delete [] fNameJetPIDtask;
  fNameJetPIDtask = 0x0;
    
  delete [] fJetPIDtask;
  fJetPIDtask = 0x0;
  
  delete [] fNameJetUEPIDtask;
  fNameJetUEPIDtask = 0x0;
    
  delete [] fJetUEPIDtask;
  fJetUEPIDtask = 0x0;  
  
  delete fAnaUtils;
  fAnaUtils = 0x0;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskIDFragmentationFunction::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // (taken from AliAnalysisTaskJetSpectrum2)
  // 
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials  = 1;

  fAvgTrials = 1;
  if(tree){
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      AliError("No current file");
      return kFALSE;
    }
    
    AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,ftrials);
    
    if (fUseInclusivePIDtask) {
      for (Int_t i = 0; i < fNumInclusivePIDtasks; i++)
        fInclusivePIDtask[i]->FillXsec(xsection);
    }
    
    if (fUseJetPIDtask) {
      for (Int_t i = 0; i < fNumJetPIDtasks; i++)
        fJetPIDtask[i]->FillXsec(xsection);
    }
    
    if (fUseJetUEPIDtask) {
      for (Int_t i = 0; i < fNumJetUEPIDtasks; i++)
        fJetUEPIDtask[i]->FillXsec(xsection);
    }    
  
    if(!fh1Xsec || !fh1Trials) {
      AliError("No Histogram fh1Xsec");
      return kFALSE;
    }
    
    fh1Xsec->Fill("<#sigma>",xsection);
    // construct a poor man average trials 
    Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
    if(ftrials>=nEntries && nEntries>0.)fAvgTrials = ftrials/nEntries;
  }

  // Set seed for backg study
  delete fRandom;
  fRandom = new TRandom3();
  TDatime* date = new TDatime();
  fRandom->SetSeed(date->Get());
  
  delete date;

  return kTRUE;
}

//__________________________________________________________________
void AliAnalysisTaskIDFragmentationFunction::UserCreateOutputObjects()
{
  // create output objects

  AliDebug(1, "Start creating user outputs");
  
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  
  //
  // Create histograms / output container
  //

  OpenFile(1);
  fCommonHistList = new TList();
  fCommonHistList->SetOwner(kTRUE);

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  
  // Histograms	
  fh1EvtSelection            = new TH1F("fh1EvtSelection", "Event Selection", 6, -0.5, 5.5);
  fh1EvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
  fh1EvtSelection->GetXaxis()->SetBinLabel(2,"event selection: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(3,"event class: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(4,"vertex Ncontr: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(5,"vertex z: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(6,"vertex type: rejected");
  
  fh1VtxSelection            = new TH1F("fh1VtxSelection", "Vertex Selection", 10, -1, 9);
  fh1VtxSelection->GetXaxis()->SetBinLabel(1,"Undef");
  fh1VtxSelection->GetXaxis()->SetBinLabel(2,"Primary");
  fh1VtxSelection->GetXaxis()->SetBinLabel(3,"Kink");
  fh1VtxSelection->GetXaxis()->SetBinLabel(4,"V0");
  fh1VtxSelection->GetXaxis()->SetBinLabel(5,"Cascade");
  fh1VtxSelection->GetXaxis()->SetBinLabel(6,"Multi");
  fh1VtxSelection->GetXaxis()->SetBinLabel(7,"SPD");
  fh1VtxSelection->GetXaxis()->SetBinLabel(8,"PileUpSPD");
  fh1VtxSelection->GetXaxis()->SetBinLabel(9,"PileUpTracks");
  fh1VtxSelection->GetXaxis()->SetBinLabel(10,"TPC");
  
  fh1VertexNContributors     = new TH1F("fh1VertexNContributors", "Vertex N contributors", 2500,-.5, 2499.5);
  fh1VertexZ                 = new TH1F("fh1VertexZ", "Vertex z distribution", 30, -15., 15.);
  fh1EvtMult 	             = new TH1F("fh1EvtMult","Event multiplicity, track pT cut > 150 MeV/c, |#eta| < 0.9",40,0.,200.);
  fh1EvtCent 	             = new TH1F("fh1EvtCent","centrality",100,0.,100.);
  
  
//   fh3trackDensity            = new TH3F("fh3TrackDensity", "Overview of Tracks in the TPC",161,85.0,246.0,90,TMath::Pi()/4.0,3.0*TMath::Pi()/4.0,360,-TMath::Pi(),TMath::Pi());
//   fh3trackDensity->GetXaxis()->SetTitle("Radial Distance");
//   fh3trackDensity->GetYaxis()->SetTitle("#theta");
//   fh3trackDensity->GetZaxis()->SetTitle("#phi");
  
  fh1Xsec                    = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fh1Trials                  = new TH1F("fh1Trials","trials from pyxsec.root",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fh1PtHard                  = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",350,-.5,349.5);
  fh1PtHardTrials            = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",350,-.5,349.5);
  
  fh1EvtsPtHardCut           = new TH1F("fh1EvtsPtHardCut", "#events before and after MC #it{p}_{T,hard} cut;;Events",2,0,2);
  fh1EvtsPtHardCut->GetXaxis()->SetBinLabel(1, "All");
  fh1EvtsPtHardCut->GetXaxis()->SetBinLabel(2, "#it{p}_{T,hard}");
  

  fh1nRecJetsCuts            = new TH1F("fh1nRecJetsCuts","reconstructed jets charged",240,-0.5,59.5);
  fh1nRecJetsCuts->GetXaxis()->SetTitle("p_{T}^{jet}");
  fh1nRecJetsCutsGroomed            = new TH1F("fh1nRecJetsCutsGroomed","reconstructed jets charged and groomed",240,-0.5,59.5);
  fh1nRecJetsCutsGroomed->GetXaxis()->SetTitle("p_{T}^{jet}");  
  fh1nRecJetsCuts2            = new TH1F("fh1nRecJetsCuts2","reconstructed jets full",40,-0.5,39.5);
  fh1nRecJetsCuts2->GetXaxis()->SetTitle("p_{T}^{jet}");
  fh1nRCinUnderground        = new TH1F("fh1nRCinUnderground", "random cones in the underground",240,-0.5,59.5);
  fh1nRCinUnderground->GetXaxis()->SetTitle("p_{T}^{jet}");
  fh1nGenJets                = new TH1F("fh1nGenJets","generated jets",40,-0.5,39.5);
  fh1nGenJets->GetXaxis()->SetTitle("p_{T}^{jet}");
  fh1TotJetEnergy            = new TH1F("fh1TotJetEnergy", "Total Jet Energy",1,0,1);
  fh2TrackDef                = new TH2F("fh2TrackDef","Comparison of Track definitions",50,0,400,50,0,2500);
  
 
  // ____________ define histograms ___________________

  fCommonHistList->Add(fh1EvtSelection);
  fCommonHistList->Add(fh1VtxSelection);
  fCommonHistList->Add(fh1EvtMult);
  fCommonHistList->Add(fh1EvtCent);
  fCommonHistList->Add(fh1VertexNContributors);
  fCommonHistList->Add(fh1VertexZ);    
  fCommonHistList->Add(fh1nRecJetsCuts);
  fCommonHistList->Add(fh1nRecJetsCutsGroomed);
  fCommonHistList->Add(fh1nRecJetsCuts2);
  fCommonHistList->Add(fh1nRCinUnderground);
  fCommonHistList->Add(fh1TotJetEnergy);
  fCommonHistList->Add(fh1Xsec);
  fCommonHistList->Add(fh1Trials);
  fCommonHistList->Add(fh1PtHard);
  fCommonHistList->Add(fh1PtHardTrials);
  fCommonHistList->Add(fh1EvtsPtHardCut);
//   fCommonHistList->Add(fh3trackDensity);
 
  if(GetJetContainer(GetNameMCParticleJetContainer())) fCommonHistList->Add(fh1nGenJets);
                                        
  fCommonHistList->Add(fh2TrackDef);
  
  // Default analysis utils
  fAnaUtils = new AliAnalysisUtils();
  
  // Not used yet, but to be save, forward vertex z cut to analysis utils object
  fAnaUtils->SetMaxVtxZ(fMaxVertexZ);

  // Load PID framework if desired
  AliDebug(1, "Loading PID framework");
  
  if (fUseJetPIDtask || fUseInclusivePIDtask || fUseJetUEPIDtask) {
    TObjArray* tasks = AliAnalysisManager::GetAnalysisManager()->GetTasks();
    if (!tasks) {
      AliError("ERROR loading PID tasks: Failed to retrieve tasks from analysis manager!");
      
      fUseInclusivePIDtask = kFALSE;
      fUseJetPIDtask = kFALSE;
      fUseJetUEPIDtask = kFALSE;
    }
    
    if (fUseInclusivePIDtask) {
      delete [] fInclusivePIDtask;
      fInclusivePIDtask = 0x0;
      
      if (fNumInclusivePIDtasks > 0) {
        fInclusivePIDtask = new AliAnalysisTaskMTFPID*[fNumInclusivePIDtasks];
        
        for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
          fInclusivePIDtask[i] = (AliAnalysisTaskMTFPID*)tasks->FindObject(fNameInclusivePIDtask[i].Data());
          
          if (!fInclusivePIDtask[i]) {
            AliErrorStream() << "ERROR Failed to load inclusive pid task" << std::endl;
            fUseInclusivePIDtask = kFALSE;
          }
        }
      }
      else {
        AliWarningStream() << "zero inclusive pid tasks!" << std::endl;
        fUseInclusivePIDtask = kFALSE;
      }
    }    
    
    if (fUseJetPIDtask) {
      delete [] fJetPIDtask;
      fJetPIDtask = 0x0;
      
      if (fNumJetPIDtasks > 0) {
        fJetPIDtask = new AliAnalysisTaskMTFPID*[fNumJetPIDtasks];
        
        for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
          fJetPIDtask[i] = (AliAnalysisTaskMTFPID*)tasks->FindObject(fNameJetPIDtask[i].Data());
          
          if (!fJetPIDtask[i]) {
            AliErrorStream() << "ERROR Failed to load jet pid task" << std::endl;
            fUseJetPIDtask = kFALSE;
          }
        }
      }
      else {
        AliWarningStream() << "zero jet pid tasks!" << std::endl;
        fUseJetPIDtask = kFALSE;
      }
    }
    
    if (fUseJetUEPIDtask) {
      delete [] fJetUEPIDtask;
      fJetUEPIDtask = 0x0;
      
      if (fNumJetUEPIDtasks > 0) {
        fJetUEPIDtask = new AliAnalysisTaskMTFPID*[fNumJetUEPIDtasks];
        
        for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
          fJetUEPIDtask[i] = (AliAnalysisTaskMTFPID*)tasks->FindObject(fNameJetUEPIDtask[i].Data());
          
          if (!fJetUEPIDtask[i]) {
            AliErrorStream() << "ERROR Failed to load jet underlying event pid task" << std::endl;
            fUseJetUEPIDtask = kFALSE;
          }
        }
      }
      else {
        AliWarningStream() << "zero jet underlying event pid tasks!" << std::endl;
        fUseJetUEPIDtask = kFALSE;
      }
    }    
   
  }
  
  const Int_t nRefMultBins = 100;
  const Double_t refMultUp = 100.;
  const Double_t refMultDown = 0.;
  
  const Int_t nJetPtBins = 20;
  const Double_t jetPtUp = 100.;
  const Double_t jetPtDown = 0.;
  
  const Int_t nCentBins = 12;
  const Double_t binsCentpp[nCentBins+1] =   { 0, 0.01, 0.1, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
  
  fhJetPtRefMultEta5 = new TH2F("fhJetPtRefMultEta5",
                                "Correlation between jet energy and event multiplicity (|#eta| < 0.5);Ref. mult. (|#eta| < 0.5);#it{p}_{T, jet}^{ch} (GeV/#it{c})",
                                nRefMultBins, refMultDown, refMultUp, nJetPtBins, jetPtDown, jetPtUp);
  fhJetPtRefMultEta8 = new TH2F("fhJetPtRefMultEta8",
                                "Correlation between jet energy and event multiplicity (|#eta| < 0.8);Ref. mult. (|#eta| < 0.8);#it{p}_{T, jet}^{ch} (GeV/#it{c})",
                                nRefMultBins, refMultDown, refMultUp, nJetPtBins, jetPtDown, jetPtUp);
  fhJetPtMultPercent  = new TH2F("fhJetPtMultPercent",
                                "Correlation between jet energy and event multiplicity percentile (V0M);Multiplicity Percentile (V0M);#it{p}_{T, jet}^{ch} (GeV/#it{c})",
                                nCentBins, binsCentpp, nJetPtBins, jetPtDown, jetPtUp);
  fCommonHistList->Add(fhJetPtRefMultEta5);
  fCommonHistList->Add(fhJetPtRefMultEta8);
  fCommonHistList->Add(fhJetPtMultPercent);
  
  if (fUseJetPIDtask && fFillDCA) {
    const Int_t nPtBins = 68;
    Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
           0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
           1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
           2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
           4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
           11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
           26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };
    
    const Int_t DCAbins = 320;
    const Double_t DCA_Z_max = 3.5;
    const Double_t DCA_XY_max = 2.5;
    
    fhDCA_XY = new TH2F("fhDCA_XY", "All rec. tracks;#it{p}_{T} (GeV/#it{c});DCA_{XY}", nPtBins, binsPt,  DCAbins, -DCA_XY_max, DCA_XY_max);
    fhDCA_Z = new TH2F("fhDCA_Z",   "All rec. tracks;#it{p}_{T} (GeV/#it{c});DCA_{Z}",  nPtBins, binsPt,  DCAbins, -DCA_Z_max,  DCA_Z_max);
    fCommonHistList->Add(fhDCA_XY);
    fCommonHistList->Add(fhDCA_Z);
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      fhDCA_XY_prim_MCID[i] = new TH2F(Form("fhDCA_XY_prim_MCID_%s", AliPID::ParticleShortName(i)),
                                       Form("Rec. %s (prim.);#it{p}_{T} (GeV/#it{c});DCA_{XY}", AliPID::ParticleLatexName(i)),
                                       nPtBins, binsPt,  DCAbins, -DCA_XY_max, DCA_XY_max);
      fhDCA_Z_prim_MCID[i]  = new TH2F(Form("fhDCA_Z_prim_MCID_%s", AliPID::ParticleShortName(i)),
                                       Form("Rec. %s (prim.);#it{p}_{T} (GeV/#it{c});DCA_{Z}", AliPID::ParticleLatexName(i)),
                                       nPtBins, binsPt,  DCAbins, -DCA_Z_max, DCA_Z_max);
      fCommonHistList->Add(fhDCA_XY_prim_MCID[i]);
      fCommonHistList->Add(fhDCA_Z_prim_MCID[i]);
      
      fhDCA_XY_sec_MCID[i] = new TH2F(Form("fhDCA_XY_sec_MCID_%s", AliPID::ParticleShortName(i)),
                                      Form("Rec. %s (sec.);#it{p}_{T} (GeV/#it{c});DCA_{XY}", AliPID::ParticleLatexName(i)),
                                      nPtBins, binsPt,  DCAbins, -DCA_XY_max, DCA_XY_max);
      fhDCA_Z_sec_MCID[i]  = new TH2F(Form("fhDCA_Z_sec_MCID_%s", AliPID::ParticleShortName(i)),
                                      Form("Rec. %s (sec.);#it{p}_{T} (GeV/#it{c});DCA_{Z}", AliPID::ParticleLatexName(i)),
                                      nPtBins, binsPt,  DCAbins, -DCA_Z_max, DCA_Z_max);
      fCommonHistList->Add(fhDCA_XY_sec_MCID[i]);
      fCommonHistList->Add(fhDCA_Z_sec_MCID[i]);
    }
  }
  
  
  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fCommonHistList->GetEntries(); ++i){
    TH1 *h1 = dynamic_cast<TH1*>(fCommonHistList->At(i));
    if (h1) h1->Sumw2();
    else{
      THnSparse *hnSparse = dynamic_cast<THnSparse*>(fCommonHistList->At(i));
      if(hnSparse) hnSparse->Sumw2();
    }
  }
  
  TH1::AddDirectory(oldStatus);

  AliDebugStream(1) << "Posting Output" << std::endl;
  
  PostData(1, fCommonHistList);
  
  // TODO should be moved to more appropriate interface
  ResetEffFunctions();
  InitialiseFastSimulationFunctions();
  
  AliDebugStream(1) << "Done" << std::endl;
}

//_______________________________________________
void AliAnalysisTaskIDFragmentationFunction::InitialiseFastSimulationFunctions()
{
  if (fUseFastSimulations && !fEffFunctions) {
    AliDebugStream(2) << "Fast Simulation parameters: " << fastSimulationParameters << std::endl;
    fEffFunctions = new TF1*[2*AliPID::kSPECIES];
    
    if (fastSimulationParameters != "") {
      AliPieceWisePoly::ReadFSParametersFromString(fastSimulationParameters, fEffFunctions);
    } else {
      //For electrons
      const Int_t parts_e = 4;
      Double_t cuts_e[parts_e-1] = {0.6,3.2,8.0};
      Int_t nparameters_e[parts_e] = {7,5,3,2};
      AliPieceWisePoly* pwp_e = new AliPieceWisePoly(parts_e,cuts_e,nparameters_e,0,50,0x0,2);
      TF1* func_e = new TF1("func_e",pwp_e,0,50,pwp_e->GetNOfParam());
      Double_t parameters_e[11] = {-1.22427e+00, 2.25003e+01, -9.00154e+01, 1.42536e+02, 1.98605e+00, -2.33708e+02, 1.74505e+02, -4.40750e-01, 1.43504e-01, -1.59226e-02, -6.14939e-04};
      func_e->SetParameters(parameters_e);
      fEffFunctions[2*AliPID::kElectron] = fEffFunctions[2*AliPID::kElectron+1] = func_e;
      
      //For protons
      const Int_t parts_p = 6;
      Double_t cuts_p[parts_p-1] = {0.4,1.6,2.5,8.0,12.0};
      Int_t nparameters_p[parts_p] = {5,4,4,2,5,2};
      AliPieceWisePoly* pwp_p = new AliPieceWisePoly(parts_p,cuts_p,nparameters_p,0,50,0x0,2);
      TF1* func_p = new TF1("func_p",pwp_p,0,50,pwp_p->GetNOfParam());
      Double_t parameters_p[12] = {-1.04124e+01, 1.72024e+02, -1.02722e+03, 2.61164e+03, -2.35742e+03, -6.20212e-01, 1.47330e-01, 8.64180e-02, -5.27106e-03, -1.91588e-01, 1.21507e-02, -2.85315e-04};
      func_p->SetParameters(parameters_p);  
      fEffFunctions[2*AliPID::kProton] = fEffFunctions[2*AliPID::kProton+1] = func_p;   

      //For kaons
      const Int_t parts_k = 5;
      Double_t cuts_k[parts_k-1] = {0.4,1.2,6,15};
      Int_t nparameters_k[parts_k] = {3,3,5,4,2};  
      AliPieceWisePoly* pwp_k = new AliPieceWisePoly(parts_k,cuts_k,nparameters_k,0,50,0x0,2);
      TF1* func_k = new TF1("func_k",pwp_k,0,50,pwp_k->GetNOfParam());
      Double_t parameters_k[9] = {-7.18856e-01, 5.10339e+00, -5.44263e+00, -4.80959e-01, 9.57122e-02, -1.80916e-02, 1.15958e-03, -6.17673e-03, 1.66119e-04};
      func_k->SetParameters(parameters_k);   
      fEffFunctions[2*AliPID::kKaon] = fEffFunctions[2*AliPID::kKaon+1] = func_k;

    //   For pions
      const Int_t parts_pi = 6;
      Double_t cuts_pi[parts_pi-1] = {0.8,1.6,3.0,10.0,12.0};
      Int_t nparameters_pi[parts_pi] = {9,4,4,3,3,2};   
      AliPieceWisePoly* pwp_pi = new AliPieceWisePoly(parts_pi,cuts_pi,nparameters_pi,0,50,0x0,2);
      TF1* func_pi = new TF1("func_pi",pwp_pi,0,50,pwp_pi->GetNOfParam());
      Double_t parameters_pi[15] = {-1.87482e-01, 8.99878e+00, -3.34776e+01, 5.73258e+01, -2.41936e+01, -4.02440e+01, 1.61212e+01, 5.09543e+01, -3.49975e+01, 7.70485e-02, -4.42996e-02, 3.16051e-01, -3.93133e-02, -5.79754e-04, 1.34446e-04};
      func_pi->SetParameters(parameters_pi);   
      fEffFunctions[2*AliPID::kPion] = fEffFunctions[2*AliPID::kPion+1] = func_pi;       
    }
  }
}

//_____________________________________________________________
Bool_t AliAnalysisTaskIDFragmentationFunction::FillHistograms() 
{	
  AliDebugStream(1) << "Start FillHistograms" << std::endl;
  
  AliDebugStream(1) << "Analyse Event #" << fEntry << std::endl;
  
  fMCEvent = MCEvent();
  
  if(!fMCEvent){
    AliDebugStream(3) << "MCEvent not found in the input" << std::endl;
  }
  
  // Extract pThard and nTrials in case of MC. 
  
  Double_t ptHard = 0.;
  Double_t nTrials = 1; // trials for MC trigger weight for real data

  if(fMCEvent && kFALSE) {
    AliGenEventHeader* genHeader = fMCEvent->GenEventHeader();
    if(genHeader){
      
      AliGenPythiaEventHeader*  pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
      AliGenHijingEventHeader*  hijingGenHeader = 0x0;
      
      if(pythiaGenHeader){
        AliDebugStream(3) << "pythiaGenHeader found" << std::endl;
        nTrials = pythiaGenHeader->Trials();
        ptHard  = pythiaGenHeader->GetPtHard();
      } else { // no pythia, hijing?
        AliDebugStream(3) << "no pythiaGenHeader found" << std::endl;
        
        hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);
        
        if(!hijingGenHeader){
          AliWarningStream() << "no pythiaGenHeader or hjingGenHeader found" << std::endl;
        } else {
            AliDebugStream(3) << "hijingGenHeader found" << std::endl;
        }
      }
    }
  }
  
  // Cut on pThard if fMCEvent and pThard >= 0 and fill histo with #evt before and after the cut
  if (fMCEvent) {
    // Before cut
    fh1EvtsPtHardCut->Fill(0.); 
    
    FillPIDTasksCutHisto(0.0, AliAnalysisTaskMTFPID::kMCPtHardCut);
    
    // Cut
    if (fMCPtHardCut >= 0. && ptHard >= fMCPtHardCut) {
      AliDebugStream(3) << "skipping event with pThard " << ptHard << " (>= " << fMCPtHardCut << ")" << std::endl;
      PostData(1, fCommonHistList);
      return kFALSE;
    }
    
    // After cut
    fh1EvtsPtHardCut->Fill(1.);
    
    FillPIDTasksCutHisto(1.0, AliAnalysisTaskMTFPID::kMCPtHardCut);    
  }
  
  // Trigger selection
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  
  if(!(inputHandler->IsEventSelected() & fEvtSelectionMask)){
    fh1EvtSelection->Fill(1.);
    AliDebugStream(1) << "Trigger Selection: event REJECTED .." << std::endl;
    PostData(1, fCommonHistList);
    return kFALSE;
  }
  
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fESD){
    AliDebugStream(3) << "ESDEvent not found in the input" << std::endl;
  }
  
  // get AOD event from input/ouput
  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if( handler && handler->InheritsFrom("AliAODInputHandler") ) {
    fAOD  =  ((AliAODInputHandler*)handler)->GetEvent();
    if(fUseAODInputJets) fAODJets = fAOD;
    AliDebugStream(1) << "AOD event from input" << std::endl;
  }
  else {
    handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( handler && handler->InheritsFrom("AliAODHandler") ) {
      fAOD = ((AliAODHandler*)handler)->GetAOD();
      fAODJets = fAOD;
      AliDebugStream(1) << "AOD event from output" << std::endl;
    }
  }
  
  if (!handler)
      handler = inputHandler;
  
  if(!fAODJets && !fUseAODInputJets){ // case we have AOD in input & output and want jets from output
    TObject* outHandler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( outHandler && outHandler->InheritsFrom("AliAODHandler") ) {
      fAODJets = ((AliAODHandler*)outHandler)->GetAOD();
      AliDebugStream(1) << "jets from output AOD" << std::endl;
    }
  }
  
  if(fNonStdFile.Length()!=0){
    // case we have an AOD extension - fetch the jets from the extended output
    
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
    fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);    
    if(!fAODExtension){
      AliDebugStream(1) << "AODExtension not found for " << fNonStdFile << std::endl;
    }
  }
  
  if(!fAOD){
    AliDebugStream(3) << "AODEvent not found" << std::endl;
  }
  if(!fAODJets){
    AliDebugStream(3) << "AODEvent with jet branch not found" << std::endl;
  }

  
  // event selection **************************************************
  // *** event class ***
  AliVEvent* evtForCentDetermination = handler->InheritsFrom("AliAODInputHandler") ? fAOD : InputEvent();
  
  Int_t multiplicity;
  if (strcmp(evtForCentDetermination->ClassName(), "AliESDEvent") == 0)
    multiplicity = ((AliESDEvent*)evtForCentDetermination)->GetNumberOfTPCTracks();
  else
    multiplicity = evtForCentDetermination->GetNumberOfESDTracks();
  
  if (GetMaxMultiplicity() != -1) {
    if (multiplicity < GetMinMultiplicity() || multiplicity > GetMaxMultiplicity())
      return kFALSE;
  }
  
  //TODO: Simplify getting centrality. Also should not depend on fIsPP, centrality can be always estimated
  
  Double_t centPercent = fCent;
  if (fEventClass > -1) {
    if (fCentBin != fEventClass) {
      // event not in selected event class, reject event
      AliDebugStream(1) << "event not in selected event class: event REJECTED ..." << std::endl;
      
      fh1EvtSelection->Fill(2.);
      PostData(1, fCommonHistList);
      return kFALSE;
    }
  }
  
  //Lines can be removed probably
  if (fCentEst.Contains("NoCentrality",TString::kIgnoreCase)) {
    centPercent = -1;
  }
  
  AliPIDResponse *pidResponse = inputHandler->GetPIDResponse();
  if (!pidResponse) {
    AliFatal("PIDResponse object was not created");
  }
  
  pidResponse->SetCurrentCentrality(centPercent);

  // Retrieve reference multiplicities in |eta|<0.8 and <0.5
  const Int_t refMult5 = 0; //((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb05();
  const Int_t refMult8 = 0; //((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08();
  const Double_t centPercentPP = fAnaUtils->GetMultiplicityPercentile(InputEvent(), "V0M");
  
  
  IncrementPIDTasksEventCounts(centPercent, AliAnalysisTaskMTFPID::kTriggerSel);

  // *** vertex cut ***
  const AliVVertex* primVtx = InputEvent()->GetPrimaryVertex();
	if (!primVtx) {
    AliErrorStream() << "Primary vertex not found " << std::endl;
		return kFALSE;
  }
	
  Int_t nTracksPrim = primVtx->GetNContributors();
  fh1VertexNContributors->Fill(nTracksPrim);
  
  AliDebugStream(1) << "primary vertex selection: " << nTracksPrim << std::endl;
  if(nTracksPrim <= 0) {
    AliDebugStream(1) << "primary vertex selection: event REJECTED..." << std::endl;
    fh1EvtSelection->Fill(3.);
    PostData(1, fCommonHistList);
    return kFALSE;
  }
  
  TString primVtxName(primVtx->GetName());
  if(primVtxName.CompareTo("TPCVertex",TString::kIgnoreCase) == 1){
    AliDebugStream(1) << "primary vertex selection: TPC vertex, event REJECTED..." << std::endl;
    fh1EvtSelection->Fill(5.);
    PostData(1, fCommonHistList);
    return kFALSE;
  }
  
  // Count events with trigger selection and vtx cut, note: Set centrality percentile fix to -1 for pp for PID framework
  IncrementPIDTasksEventCounts(centPercent, AliAnalysisTaskMTFPID::kTriggerSelAndVtxCut); 
  
  fh1VertexZ->Fill(primVtx->GetZ());
  
  if(TMath::Abs(primVtx->GetZ())>fMaxVertexZ) {
    AliDebugStream(1) << "primary vertex z = " << primVtx->GetZ() << ": event REJECTED..." << std::endl;
    fh1EvtSelection->Fill(4.);
    PostData(1, fCommonHistList);
    return kFALSE; 
  }
  
  // Count events with trigger selection, vtx cut and z vtx cut
  IncrementPIDTasksEventCounts(centPercent, AliAnalysisTaskMTFPID::kTriggerSelAndVtxCutAndZvtxCutNoPileUpRejection);
  
  // Store for each task, whether this task would tag this event as pile-up or not
  const Int_t arrSizeInclusive = TMath::Max(1, fNumInclusivePIDtasks);
  const Int_t arrSizeJet = TMath::Max(1, fNumJetPIDtasks);
  const Int_t arrSizeJetUE = TMath::Max(1, fNumJetUEPIDtasks);
  Bool_t isPileUpInclusivePIDtask[arrSizeInclusive];
  Bool_t isPileUpJetPIDtask[arrSizeJet];
  Bool_t isPileUpJetUEPIDtask[arrSizeJetUE];
  
  for (Int_t i = 0; i < arrSizeInclusive; i++) 
    isPileUpInclusivePIDtask[i] = kFALSE;  
  
  for (Int_t i = 0; i < arrSizeJet; i++) 
    isPileUpJetPIDtask[i] = kFALSE;
  
  for (Int_t i = 0; i < arrSizeJetUE; i++) 
    isPileUpJetUEPIDtask[i] = kFALSE;  
  
  // Check whether there is at least one task that does not reject the event (saves processing time in the following code)
  Bool_t isPileUpForAllInclusivePIDTasks = kTRUE;
  Bool_t isPileUpForAllJetPIDTasks = kTRUE;
  Bool_t isPileUpForAllJetUEPIDTasks = kTRUE;
  
  if (fUseInclusivePIDtask) {
    for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
      isPileUpInclusivePIDtask[i] = fInclusivePIDtask[i]->GetIsPileUp(InputEvent(), fInclusivePIDtask[i]->GetPileUpRejectionType());
      isPileUpForAllInclusivePIDTasks = isPileUpForAllInclusivePIDTasks && isPileUpInclusivePIDtask[i];
    }
  }
  
  if (fUseJetPIDtask) {
    for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
      isPileUpJetPIDtask[i] = fJetPIDtask[i]->GetIsPileUp(InputEvent(), fJetPIDtask[i]->GetPileUpRejectionType());
      isPileUpForAllJetPIDTasks = isPileUpForAllJetPIDTasks && isPileUpJetPIDtask[i];
    }
  }
  
  if (fUseJetUEPIDtask) {
    for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
      isPileUpJetUEPIDtask[i] = fJetUEPIDtask[i]->GetIsPileUp(InputEvent(), fJetUEPIDtask[i]->GetPileUpRejectionType());
      isPileUpForAllJetUEPIDTasks = isPileUpForAllJetUEPIDTasks && isPileUpJetUEPIDtask[i];
    }
  }
    
  // Count events with trigger selection, vtx cut, z vtx cut and after pile-up rejection (if enabled in that task)
  IncrementPIDTasksEventCounts(centPercent, AliAnalysisTaskMTFPID::kTriggerSelAndVtxCutAndZvtxCut, isPileUpInclusivePIDtask, isPileUpJetPIDtask, isPileUpJetUEPIDtask);
  
  AliDebugStream(1) << "event ACCEPTED ..." << std::endl;

  fh1EvtSelection->Fill(0.);
  
  const AliAODVertex* aodVertex = dynamic_cast<const AliAODVertex*>(primVtx);
  if (aodVertex) {
      fh1VtxSelection->Fill(aodVertex->GetType());
  }
  
  fh1EvtCent->Fill(centPercent);

  // Set centrality percentile fix to -1 for pp to be used for the PID framework
  if (fIsPP)
    centPercent = -1;
  
  // Call ConfigureTaskForCurrentEvent of PID tasks to ensure that everything is set up properly for the current event
  // (e.g. run/period dependence of max eta variation map)
  if (fUseInclusivePIDtask) {
    for (Int_t i = 0; i < fNumInclusivePIDtasks; i++)
      if (!isPileUpInclusivePIDtask[i])
        fInclusivePIDtask[i]->ConfigureTaskForCurrentEvent(InputEvent());
  }
  
  if (fUseJetPIDtask) {
    for (Int_t i = 0; i < fNumJetPIDtasks; i++)
      if (!isPileUpJetPIDtask[i])
        fJetPIDtask[i]->ConfigureTaskForCurrentEvent(InputEvent());
  }
  
  if (fUseJetUEPIDtask) {
    for (Int_t i = 0; i < fNumJetUEPIDtasks; i++)
      if (!isPileUpJetUEPIDtask[i])
        fJetUEPIDtask[i]->ConfigureTaskForCurrentEvent(InputEvent());
  }  
  
  
  //___ fill MC information __________________________________________________________________

  fh1Trials->Fill("#sum{ntrials}",fAvgTrials); 
  
  if (fUseInclusivePIDtask) {
    for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
      if (!isPileUpInclusivePIDtask[i])
        fInclusivePIDtask[i]->FillPythiaTrials(fAvgTrials);
    }
  }
  
  if (fUseJetPIDtask) {
    for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
      if (!isPileUpJetPIDtask[i])
        fJetPIDtask[i]->FillPythiaTrials(fAvgTrials);
    }
  }
  
  if (fUseJetUEPIDtask) {
    for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
      if (!isPileUpJetUEPIDtask[i])
        fJetUEPIDtask[i]->FillPythiaTrials(fAvgTrials);
    }
  }  

  if (fMCEvent) {
    fh1PtHard->Fill(ptHard);
    fh1PtHardTrials->Fill(ptHard,nTrials);
  }
  
  Bool_t tuneOnDataTPC = pidResponse->IsTunedOnData() &&
                        ((pidResponse->GetTunedOnDataMask() & AliPIDResponse::kDetTPC) == AliPIDResponse::kDetTPC);
  
  AliDebugStream(1) << "Starting processing..." << std::endl;
  
  //____ analysis, fill histos ___________________________________________________
  
  AliTrackContainer *trackContainer = GetTrackContainer(GetNameTrackContainer());
  AliMCParticleContainer *mcParticleContainer = GetMCParticleContainer(GetNameMCParticleContainer());
  Double_t trackEtaMin = -0.9;
  Double_t trackEtaMax = 0.9;
  if (trackContainer) {
    trackEtaMin = trackContainer->GetParticleEtaMin();
    trackEtaMax = trackContainer->GetParticleEtaMax();
  }
  
  //Prepare V0 Indexes
  
  AliAnalysisTaskPIDV0base* pidTaskWithV0Indexes = 0x0;
  
//   if (fUseInclusivePIDtask && fInclusivePIDtask[0]->GetDoTPCclusterStudies()) {
//     fInclusivePIDtask[0]->FillV0PIDlist(fESD);
//     pidTaskWithV0Indexes=fInclusivePIDtask[0];
//   }
//   
//   if (!pidTaskWithV0Indexes && fUseJetPIDtask && fJetPIDtask[0]->GetDoTPCclusterStudies()) {
//     fJetPIDtask[0]->FillV0PIDlist(fESD);
//     pidTaskWithV0Indexes=fJetPIDtask[0];    
//   }
//   
//   if (!pidTaskWithV0Indexes && fUseJetUEPIDtask && fJetUEPIDtask[0]->GetDoTPCclusterStudies()) {
//     fJetUEPIDtask[0]->FillV0PIDlist(fESD);
//     pidTaskWithV0Indexes=fJetUEPIDtask[0];    
//   } 
      
  // Fill efficiency for generated primaries and also fill histos for generated yields (primaries + all)
  // Efficiency, inclusive - particle level
  AliDebugStream(2) << "Starting Inclusive Efficiency..." << std::endl;
  
  if (fUseInclusivePIDtask && mcParticleContainer && !isPileUpForAllInclusivePIDTasks) {
    for (auto part : mcParticleContainer->accepted()) {

      if (!part)
        continue;
      
      // Define clean MC sample with corresponding particle level track cuts:
      // - MC-track must be in desired eta range
      // - MC-track must be physical primary
      // - Species must be one of those in question
        
      if (part->Eta() > trackEtaMax || part->Eta() < trackEtaMin)
        continue;
      
      Int_t mcID = AliAnalysisTaskMTFPID::PDGtoMCID(part->GetPdgCode());
      
      // Following lines are not needed - just keep other species (like casecades) - will end up in overflow bin
      // and only affect the efficiencies for all (i.e. not identified) what is desired!
      //if (mcID == AliPID::kUnknown)
      //  continue;
      
      if (!part->IsPhysicalPrimary())
        continue;

      Double_t pT = part->Pt();
      
      // AliAODMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
      Double_t chargeMC = part->Charge() / 3.;
      
      if (TMath::Abs(chargeMC) < 0.01)
        continue; // Reject neutral particles (only relevant, if mcID is not used)
      
      Double_t valuesGenYield[AliAnalysisTaskMTFPID::kGenYieldNumAxes] = { static_cast<Double_t>(mcID), pT, centPercent,
                                                                        -1, -1, -1, -1, -1, -1 };

      for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
        if (!isPileUpInclusivePIDtask[i] && fInclusivePIDtask[i]->IsInAcceptedEtaRange(TMath::Abs(part->Eta()))) {
          valuesGenYield[fInclusivePIDtask[i]->GetIndexOfChargeAxisGenYield()] = fInclusivePIDtask[i]->GetStoreCharge() ? chargeMC : -2;
          fInclusivePIDtask[i]->FillGeneratedYield(valuesGenYield);
        }
      }
      
      // Always store the charge for the efficiency (needed for geant-fluka correction)
      Double_t valuesEff[AliAnalysisTaskMTFPID::kEffNumAxes] = { static_cast<Double_t>(mcID), pT, part->Eta(), chargeMC,
                                                              centPercent, -1, -1, -1, -1, -1 };// no jet pT etc since inclusive spectrum 
      for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
        if (!isPileUpInclusivePIDtask[i])
          fInclusivePIDtask[i]->FillEfficiencyContainer(valuesEff, AliAnalysisTaskMTFPID::kStepGenWithGenCuts);
      }
      if (fUseFastSimulations) { 
        //Fast simulations for inclusive particles - only used to check the parametrization
        if (mcID == AliPID::kMuon || mcID >= AliPID::kSPECIES)
          continue;
        
        Bool_t posCharge = chargeMC > 0.0;
        
        Double_t xeff = fFastSimEffFactor * fEffFunctions[2*mcID+(Int_t)posCharge]->Eval(pT);
        Double_t x = fRandom->Rndm();
        if (x > xeff)
          continue;

        Double_t smearedPt = 1.0/(fRandom->Gaus(1.0/pT,fFastSimResFactor * fFastSimRes));

        Double_t valuesEffFastMeasObs[AliAnalysisTaskMTFPID::kEffNumAxes] = { static_cast<Double_t>(mcID), smearedPt, part->Eta(), chargeMC,
                                                                centPercent, -1, -1, -1, -1, -1 };// no jet pT etc since inclusive spectrum                                                                   
        for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
          if (!isPileUpInclusivePIDtask[i]) {
            fInclusivePIDtask[i]->FillEfficiencyContainer(valuesEffFastMeasObs, AliAnalysisTaskMTFPID::kStepRecWithRecCutsMeasuredObsPrimaries);
          }
        }
      }
    }
  }
  
  AliTrackContainer* trackContainerEfficiency = GetTrackContainer(GetNameTrackContainerEfficiency());
  if (fUseInclusivePIDtask && mcParticleContainer && trackContainerEfficiency && !isPileUpForAllInclusivePIDTasks) {
    //Efficiency, inclusive - detector level
    for(auto inclusiveaod : trackContainerEfficiency->accepted()) {
      // fill inclusive tracks XXX, they have the same track cuts!

      if(!inclusiveaod)
        continue;
      
      Double_t dEdxTPC = tuneOnDataTPC ? pidResponse->GetTPCsignalTunedOnData(inclusiveaod) 
                                        : inclusiveaod->GetTPCsignal();
      
      if (dEdxTPC <= 0)
        continue;      
      
      Bool_t survivedTPCCutMIGeo = AliAnalysisTaskMTFPID::TPCCutMIGeo(inclusiveaod, InputEvent());
      Bool_t survivedTPCnclCut = AliAnalysisTaskMTFPID::TPCnclCut(inclusiveaod);        //Included in the cut above
      
      Int_t label = TMath::Abs(inclusiveaod->GetLabel());
      
      // find MC track in our list, if available
      AliAODMCParticle* gentrack = mcParticleContainer->GetMCParticleWithLabel(label);

      // For efficiency: Reconstructed track has survived all cuts on the detector level (no cut on eta acceptance)
      // and has an associated MC track
      // -> Check whether associated MC track belongs to the clean MC sample defined above,
      //    i.e. survives the particle level track cuts
      if (gentrack && !fUseFastSimulations) {
        Int_t pdg = gentrack->GetPdgCode();
        Int_t mcID = AliAnalysisTaskMTFPID::PDGtoMCID(pdg);
        
        // Following lines are not needed - just keep other species (like casecades) - will end up in overflow bin
        // and only affect the efficiencies for all (i.e. not identified) what is desired!
        //if (mcID == AliPID::kUnknown)
        //  continue;
        
        // Fill efficiency for reconstructed primaries
        if (!gentrack->IsPhysicalPrimary())
          continue;
          
        if (gentrack->Eta() > trackEtaMax || gentrack->Eta() < trackEtaMin)
          continue;
          
        // AliAODMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
        // Always store the charge for the efficiency (needed for geant-fluka correction)
        Double_t value[AliAnalysisTaskMTFPID::kEffNumAxes] = { static_cast<Double_t>(mcID), gentrack->Pt(), gentrack->Eta(),
                                                            gentrack->Charge() / 3., centPercent, -1, -1, -1, -1, -1 };// no jet pT etc since inclusive spectrum 
        for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
          if (!isPileUpInclusivePIDtask[i] && (!fInclusivePIDtask[i]->GetUseTPCCutMIGeo() || survivedTPCCutMIGeo) && (!fInclusivePIDtask[i]->GetUseTPCnclCut() || survivedTPCnclCut))
            fInclusivePIDtask[i]->FillEfficiencyContainer(value, AliAnalysisTaskMTFPID::kStepRecWithGenCuts);
        }
            
        Double_t valueMeas[AliAnalysisTaskMTFPID::kEffNumAxes] = { static_cast<Double_t>(mcID), inclusiveaod->Pt(), inclusiveaod->Eta(),
                                                                static_cast<Double_t>(inclusiveaod->Charge()), centPercent,
                                                                -1, -1, -1, -1, -1 };// no jet pT etc since inclusive spectrum 
        for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
          if (!isPileUpInclusivePIDtask[i] && (!fInclusivePIDtask[i]->GetUseTPCCutMIGeo() || survivedTPCCutMIGeo) && (!fInclusivePIDtask[i]->GetUseTPCnclCut() || survivedTPCnclCut))
            fInclusivePIDtask[i]->FillEfficiencyContainer(valueMeas, AliAnalysisTaskMTFPID::kStepRecWithGenCutsMeasuredObs);
        }
      }
    }
  }  
  
  AliTrackContainer* v0TrackContainer = GetTrackContainer("V0TrackContainer"); 
  
  if (trackContainer) {
    Int_t accTracks = trackContainer->GetNAcceptedTracks();
    
    fh2TrackDef->Fill(accTracks, multiplicity);
    
    fh1EvtMult->Fill(accTracks);
  }
  
  if (fUseInclusivePIDtask && v0TrackContainer && !isPileUpForAllInclusivePIDTasks) {
    for (auto part : v0TrackContainer->accepted()) {
      if (!part)
        continue;
      Double_t dEdxTPC = tuneOnDataTPC ? pidResponse->GetTPCsignalTunedOnData(part)
                                        : part->GetTPCsignal();    
                                        
      if (dEdxTPC <= 0)
        continue;  

//       if (part->GetTPCsignalN() < 120)
//         continue;
      
      Int_t v0tag = (Int_t)pidTaskWithV0Indexes->GetV0tag(part->GetID());
      
      for (Int_t i=0;i<fNumInclusivePIDtasks;++i) {
        if (!isPileUpInclusivePIDtask[i]) {
          if (OverlapsWithAnyRecJet(part)) {
//             if (fJetPIDtask[i])
//               fJetPIDtask[i]->DoTPCclusterStudies(part, v0tag, multiplicity, dEdxTPC, 1, fNumInclusivePIDtasks > 1 ? trackContainer : 0x0);
          }
          else {
            fInclusivePIDtask[i]->DoTPCclusterStudies(part, v0tag, multiplicity, dEdxTPC, 1, 0x0);
          }
        }
      }      
    }
  }
  
  AliDebugStream(2) << "Process inclusive tracks..." << std::endl;
  
  if (fUseInclusivePIDtask && trackContainer && !isPileUpForAllInclusivePIDTasks) {
    for(auto part : trackContainer->accepted()) {
      if (!part) 
        continue;
      
      Int_t label = TMath::Abs(part->GetLabel());

//       find MC track in our list, if available
      AliAODMCParticle* gentrack = mcParticleContainer ? mcParticleContainer->GetMCParticleWithLabel(label) : 0x0;
      Int_t pdg = 0;
      
      if (gentrack)
        pdg = gentrack->GetPdgCode();      
      
      Double_t dEdxTPC = tuneOnDataTPC ? pidResponse->GetTPCsignalTunedOnData(part)
                                        : part->GetTPCsignal();
      
      if (dEdxTPC <= 0)
        continue;   
      
//       if (part->GetTPCsignalN() < 120)
//         continue;
      
//       for (Int_t i=0;i<fNumInclusivePIDtasks;++i) {
//         if (!isPileUpInclusivePIDtask[i]) {
//           if (fJetPIDtask[i]) {
//             if (fJetPIDtask[i])
//               fJetPIDtask[i]->DoTPCclusterStudies(part, 0, multiplicity, dEdxTPC, 0, fNumInclusivePIDtasks > 1 ? trackContainer : 0x0);
//           }
//           else {
//             fInclusivePIDtask[i]->DoTPCclusterStudies(part, 0, multiplicity, dEdxTPC, 0, 0x0);
//           }
//         }
//       }
      
      //ATTENTION: Only here for the TPCcluster studies, comment out or delete otherwise
      //continue; 
     
      Bool_t survivedTPCCutMIGeo = AliAnalysisTaskMTFPID::TPCCutMIGeo(part, InputEvent());
      Bool_t survivedTPCnclCut = AliAnalysisTaskMTFPID::TPCnclCut(part);   //Included above
      
      for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
        if (!isPileUpInclusivePIDtask[i] && (!fInclusivePIDtask[i]->GetUseTPCCutMIGeo() || survivedTPCCutMIGeo) && (!fInclusivePIDtask[i]->GetUseTPCnclCut() || survivedTPCnclCut)) {
            if (fInclusivePIDtask[i]->IsInAcceptedEtaRange(TMath::Abs(part->Eta()))) {
                fInclusivePIDtask[i]->ProcessTrack(part, pdg, centPercent, -1, kFALSE, kTRUE, -1, -1); // no jet pT etc since inclusive spectrum 
                fInclusivePIDtask[i]->FillDCA(part, dEdxTPC, primVtx, fMCEvent, -1.0);
            }
        }
      }
      
      if (gentrack && !fUseFastSimulations) {
        Int_t mcID = AliAnalysisTaskMTFPID::PDGtoMCID(pdg);
        // Always store the charge for the efficiency (needed for geant-fluka correction)
        Double_t valueRecAllCuts[AliAnalysisTaskMTFPID::kEffNumAxes] = { static_cast<Double_t>(mcID), part->Pt(), part->Eta(), 
                                                                      static_cast<Double_t>(part->Charge()), centPercent,
                                                                      -1, -1, -1, -1, -1 };
        for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
          if (!isPileUpInclusivePIDtask[i] && (!fInclusivePIDtask[i]->GetUseTPCCutMIGeo() || survivedTPCCutMIGeo) && (!fInclusivePIDtask[i]->GetUseTPCnclCut() || survivedTPCnclCut))
            if (!fUseFastSimulations) 
              fInclusivePIDtask[i]->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskMTFPID::kStepRecWithRecCutsMeasuredObs);
        }
        
        Double_t weight = IsSecondaryWithStrangeMotherMC(gentrack, mcParticleContainer)
                            ? GetMCStrangenessFactorCMS(gentrack, mcParticleContainer) : 1.0;
        for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
          if (!isPileUpInclusivePIDtask[i] && (!fInclusivePIDtask[i]->GetUseTPCCutMIGeo() || survivedTPCCutMIGeo) && (!fInclusivePIDtask[i]->GetUseTPCnclCut() || survivedTPCnclCut))
            fInclusivePIDtask[i]->FillEfficiencyContainer(valueRecAllCuts, 
                                                          AliAnalysisTaskMTFPID::kStepRecWithRecCutsMeasuredObsStrangenessScaled,
                                                          weight);
        }
        
        if (gentrack->IsPhysicalPrimary()) {
          // AliAODMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
          // Always store the charge for the efficiency (needed for geant-fluka correction)
          Double_t valueGenAllCuts[AliAnalysisTaskMTFPID::kEffNumAxes] = { static_cast<Double_t>(mcID), gentrack->Pt(), gentrack->Eta(), 
                                                                        gentrack->Charge() / 3., centPercent, -1, -1, -1, -1, -1 };
          
          Double_t valuePtResolution[AliAnalysisTaskMTFPID::kPtResNumAxes] = { -1, gentrack->Pt(), part->Pt(),
                                                                            gentrack->Charge() / 3., centPercent };
        
          for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
            if (!isPileUpInclusivePIDtask[i] && (!fInclusivePIDtask[i]->GetUseTPCCutMIGeo() || survivedTPCCutMIGeo) && (!fInclusivePIDtask[i]->GetUseTPCnclCut() || survivedTPCnclCut)) {
              fInclusivePIDtask[i]->FillEfficiencyContainer(valueRecAllCuts, 
                                                            AliAnalysisTaskMTFPID::kStepRecWithRecCutsMeasuredObsPrimaries);
              fInclusivePIDtask[i]->FillEfficiencyContainer(valueGenAllCuts, 
                                                            AliAnalysisTaskMTFPID::kStepRecWithRecCutsPrimaries);
              
              fInclusivePIDtask[i]->FillPtResolution(mcID, valuePtResolution);
            }
          }
        }
      }
    }    
  }
  
  AliDebugStream(2) << "Process jets..." << std::endl;
  
  AliJetContainer* mcJetContainer = GetJetContainer(GetNameMCParticleJetContainer()); 
  
  if (fUseJetPIDtask && mcJetContainer && !isPileUpForAllJetPIDTasks && !fUseFastSimulations) {
    AliDebugStream(2) << "Process Jets - efficiency, particle level.." << std::endl;
    
    if (fOnlyLeadingJets)
      mcJetContainer->SortArray();    
    
    for (Int_t i=0;i<fNumJetPIDtasks;i++) {
      
      if (!isPileUpJetPIDtask[i]) {

        for(auto jet : mcJetContainer->accepted()) {
            
          if(!jet) 
            continue;
          
          Float_t jetPt   = jet->Pt();
          
          if (mcJetContainer->GetRhoParameter())
            jetPt = jetPt - mcJetContainer->GetRhoVal() * jet->Area();  
          
          fh1nGenJets->Fill(jetPt);
          
          fJetPIDtask[i]->FillGenJets(centPercent, jetPt);
          
          for(Int_t it=0; it<jet->GetNumberOfTracks(); ++it) {
            AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(jet->Track(it));
            PerformJetMonteCarloAnalysisGeneratedYield(jet, trackVP, fJetPIDtask[i], centPercent);
          }
          
          if (fOnlyLeadingJets)
            break;
        }    
      }
    }
  }
  
  // Fast simulations 
  // TODO: Extend code so it can run all fast simulations in the same train (reducing output numbers and complexity)
  if (fUseJetPIDtask && mcJetContainer && !isPileUpForAllJetPIDTasks && fUseFastSimulations) {
    AliDebugStream(2) << "Do fast simulation" << std::endl;
    for (Int_t i=0;i<fNumJetPIDtasks;i++) {
      if (!isPileUpJetPIDtask[i]) {
        AliFJWrapper* wrapper = new AliFJWrapper("wrapper", "wrapper");        
        SetUpFastJetWrapperWithOriginalValues(wrapper);  
        Double_t jetMinPt = mcJetContainer->GetMinPt();
        Double_t jetMaxEta = mcJetContainer->GetMaxEta() - mcJetContainer->GetJetRadius();
        Double_t jetMinEta = mcJetContainer->GetMinEta() + mcJetContainer->GetJetRadius(); 
        AliHelperClassFastSimulation* fastSimulation = new AliHelperClassFastSimulation(fEffFunctions, fFastSimEffFactor, fFastSimResFactor, fFastSimRes, jetMinPt, jetMaxEta, jetMinEta, wrapper);

        Int_t nextFreeIndex;
        
        if (fFFChange == kLowPtEnhancement) {
          Int_t maxLabel = 0;
          for (auto part : mcParticleContainer->accepted()) {
            if (!part)
              continue;
                
            maxLabel = TMath::Max(maxLabel, part->GetLabel());
          }
            
          nextFreeIndex = maxLabel + 1;
        }

        AliDebugStream(4) << "Free index: " << nextFreeIndex << std::endl;
        
        for(auto jet : mcJetContainer->accepted()) {
          if(!jet) 
            continue;
          
          Double_t jetPt   = jet->Pt();
          
          if (mcJetContainer->GetRhoParameter())
            jetPt = jetPt - mcJetContainer->GetRhoVal() * jet->Area();  
          
          fh1nGenJets->Fill(jetPt);
          
          fJetPIDtask[i]->FillGenJets(centPercent, jetPt);
          
          AliAODMCParticle* leadingTrack = dynamic_cast<AliAODMCParticle*>(jet->GetLeadingTrack());
          
          Double_t leadingTrackPt = leadingTrack->Pt();
          Double_t smallestTrackPt = leadingTrackPt;
          
          const Int_t nOfJetTracks = jet->GetNumberOfTracks();
          
          for(Int_t it=0; it<nOfJetTracks; ++it) {
            AliAODMCParticle* track = dynamic_cast<AliAODMCParticle*>(jet->Track(it));
            if (track != leadingTrack && fFFChange == kLowPtDepletion) {
              if (fRandom->Rndm() > 0.75)
                continue;
            }
            FillEfficiencyContainerFromTrack(track, jet, centPercent, AliAnalysisTaskMTFPID::kStepGenWithGenCuts);
            
            fastSimulation->AddParticle(track);
            
            if (fFFChange == kLowPtEnhancement) {
              smallestTrackPt = TMath::Min(track->Pt(), smallestTrackPt);	
            }
          }
                            
          if (fFFChange == kLowPtEnhancement && nOfJetTracks > 1) {
            Double_t jetPhi = jet->Phi();
            Double_t jetTheta = jet->Theta();
            
            const Double_t survivalChanceSlope = 0.5/(leadingTrackPt - smallestTrackPt);
        
            for(Int_t it=0; it<nOfJetTracks; ++it) {
              AliAODMCParticle *track  = dynamic_cast<AliAODMCParticle*>(jet->Track(it));
                
              if (!track || jet->Track(it) == leadingTrack)
                continue;  
              
              Double_t pt = track->Pt();
              // Discard additional tracks randomly. Survival probability going linearly from 75% for the lowest to 25% for the highest momentum
              Double_t survivalChance = 0.75 - survivalChanceSlope * (pt - smallestTrackPt);
              if (fRandom->Rndm() > survivalChance)
                continue;
            
              Double_t phi = track->Phi();
              Double_t theta = track->Theta();
              
              // Rotate 180° around jet axis
              phi = 2.0 * jetPhi - phi;
              theta = 2.0 * jetTheta - theta;
              Double_t eta = -TMath::Log(TMath::Tan(0.5 * theta));
              
              Double_t px = pt * TMath::Cos(phi);
              Double_t py = pt * TMath::Sin(phi);
              Double_t pz = pt * TMath::SinH(eta);
              Double_t mass = track->M();
              Double_t energy = TMath::Sqrt(pt * pt + pz * pz + mass*mass);
              
              TParticle* tpart = new TParticle(track->GetPdgCode(), 0, 0, 0, 0, 0, px, py, pz, energy, 0.0, 0.0, 0.0, 0.0);
              AliMCParticle* mcParticle = new AliMCParticle(tpart);
              AliAODMCParticle* doubled_part = new AliAODMCParticle(mcParticle, nextFreeIndex);

              FillEfficiencyContainerFromTrack(doubled_part, jet, centPercent, AliAnalysisTaskMTFPID::kStepGenWithGenCuts);
              fastSimulation->AddParticle(doubled_part);
              
              nextFreeIndex++;	
            }
          }
        }
        
        AliDebugStream(4) << "Run the fast simulation" << std::endl;
        
        fastSimulation->Run();

        for (UInt_t j=0;j<fastSimulation->GetNJets();++j) {
          
          AliEmcalJet *jet = fastSimulation->GetJet(j);
        
          Double_t jetPt = jet->Pt();
        
          for (Int_t i=0;i<fNumJetPIDtasks;++i) {
            fJetPIDtask[i]->FillRecJets(centPercent, jetPt);
          }    

          for (UInt_t ic = 0;ic<fastSimulation->GetNParticlesOfJet(j);++ic) {
            AliAODMCParticle* part = fastSimulation->GetTrackOfJet(ic, j);
            
            if (!part)
              continue;
          
            FillEfficiencyContainerFromTrack(part, jet, centPercent, AliAnalysisTaskMTFPID::kStepRecWithRecCutsMeasuredObsPrimaries);
          }
        }
        
        delete fastSimulation;
        fastSimulation = 0x0;
      }
    }
  }
  
  if (trackContainer) {
    for (Int_t i=0;i<20;i++) {
      Double_t Eta = (trackContainer->GetMaxEta() - (GetFFRadius()/2.0)) * (2 * fRandom->Rndm() - 1.0);    //random eta value in range: [-dEtaMax+R/2,+dEtaConeMax-R/2]
      Double_t Phi = TMath::TwoPi() * fRandom->Rndm();        //random phi value in range: [0,2*Pi]
      
      Double_t jetPt = 0.0;
      for (auto track : trackContainer->accepted()) {
        if(!track)
          continue;
        
        Double_t dEta = track->Eta() - Eta;
        Double_t dPhi = track->Phi() - Phi;
        Double_t angulardistance = TMath::Sqrt(dEta * dEta + dPhi * dPhi); 
        if (angulardistance < TMath::Abs(GetFFRadius())) {
          jetPt += track->Pt();
        }
      }
      fh1nRCinUnderground->Fill(jetPt);  
    }
  }

  AliJetContainer* jetContainer = GetJetContainer(GetNameJetContainer());
  
  if (GetDoGroomedJets() && trackContainer) {
    AliFJWrapper* wrapper = new AliFJWrapper("wrapper","wrapper");
    SetUpFastJetWrapperWithOriginalValues(wrapper);    
    AliTrackIterableMomentumContainer itcont = trackContainer->accepted_momentum();
    for (AliTrackIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
      TLorentzVector pvec(it->first.Px(), it->first.Py(), it->first.Pz(), it->first.E());
      wrapper->AddInputVector(pvec.Px(), pvec.Py(), pvec.Pz(), pvec.E(), it.current_index());
    }
    wrapper->Run();
    wrapper->DoSoftDrop();
    std::vector<fastjet::PseudoJet> jets = wrapper->GetInclusiveJets();
    std::vector<fastjet::PseudoJet> groomedJets = wrapper->GetGroomedJets();
    for (UInt_t j=0;j<jets.size();++j) {
      if (jets[j].perp() < 5.0 || TMath::Abs(jets[j].eta()) > 0.5)
        continue;
      else {
        fh1nRecJetsCuts->Fill(jets[j].perp());
      }
    }
    for (UInt_t j=0;j<groomedJets.size();++j) {
      if (groomedJets[j].perp() < 5.0 || TMath::Abs(groomedJets[j].eta()) > 0.5)
        continue;
      else {
        fh1nRecJetsCutsGroomed->Fill(groomedJets[j].perp());
//         std::vector<fastjet::PseudoJet> constituents = groomedJets[j].constituents();
//         cout << "Groomed Jet: " << jets[j].perp() << " " << constituents.size() << endl;
/*        for (Int_t i=0;i<constituents.size();++i) {
          Int_t uid = constituents[i].user_index();
          if (uid == -1)    //Ghost particle
            continue;
          cout << uid << endl;
        }  */      
      }
    }    
  }
  
  
  if (jetContainer && !fUseFastSimulations) { 
    if (fOnlyLeadingJets)
      jetContainer->SortArray();
    for (auto jet : jetContainer->accepted()) {
      if (!jet)
        continue;
      
      Float_t jetPt   = jet->Pt();
      
      if (jetContainer->GetRhoParameter())
        jetPt = jetPt - jetContainer->GetRhoVal() * jet->Area();
      
//       fh1nRecJetsCuts->Fill(jetPt);
      
      fh1TotJetEnergy->Fill(0.5,jetPt);
      
      fhJetPtRefMultEta5->Fill(refMult5, jetPt);
      fhJetPtRefMultEta8->Fill(refMult8, jetPt);
      fhJetPtMultPercent->Fill(centPercentPP, jetPt);
      
      if (!(fUseJetPIDtask && !isPileUpForAllJetPIDTasks))
        continue;
      
      for (Int_t i=0;i<fNumJetPIDtasks;++i) {
        if (!isPileUpJetPIDtask[i])
          fJetPIDtask[i]->FillRecJets(centPercent, jetPt);
      }
      
      for(Int_t j=0; j<jet->GetNumberOfTracks(); ++j) {
        AliVTrack *track  = dynamic_cast<AliVTrack*>(jet->Track(j));
        if (!track)
          continue;        
        Bool_t trackRejectedByTask[arrSizeJet];
        Bool_t trackRejectedByAllTasks = kFALSE;
        if (track) {
          Bool_t survivedTPCCutMIGeo = AliAnalysisTaskMTFPID::TPCCutMIGeo(track, InputEvent());
          Bool_t survivedTPCnclCut = AliAnalysisTaskMTFPID::TPCnclCut(track);    //Included above
          Double_t dEdxTPC = tuneOnDataTPC ? pidResponse->GetTPCsignalTunedOnData(track) : track->GetTPCsignal();
          if (dEdxTPC < 0)
            continue;
          
          for (Int_t i=0;i<fNumJetPIDtasks;++i) {
            trackRejectedByTask[i] = isPileUpJetPIDtask[i] || !(!fJetPIDtask[i]->GetUseTPCCutMIGeo() || survivedTPCCutMIGeo) || !(!fJetPIDtask[i]->GetUseTPCnclCut() || survivedTPCnclCut) || !(fJetPIDtask[i]->IsInAcceptedEtaRange(TMath::Abs(track->Eta())));
            trackRejectedByAllTasks = trackRejectedByAllTasks && trackRejectedByTask[i];
          }

          if (!trackRejectedByAllTasks) {
            AnalyseJetTrack(track, jet, 0x0, trackRejectedByTask, centPercent, mcParticleContainer);
            FillDCA(track, mcParticleContainer);
//             for (Int_t i=0;i<fNumJetPIDtasks;++i) {
//               if (!trackRejectedByTask[i]) {
//                 fJetPIDtask[i]->FillDCA(track, dEdxTPC, primVtx, fMCEvent, jetPt);
//               }
//             }
          }
        }
      }
      
      if (fOnlyLeadingJets)
        break;         
    }
  }
  
  AliDebugStream(2) << "Process Underlying event..." << std::endl;
  if (fUseJetUEPIDtask && jetContainer && !isPileUpForAllJetUEPIDTasks && !fUseFastSimulations) {
    for (Int_t i=0;i<fNumJetUEPIDtasks;++i) {
      if (!isPileUpJetUEPIDtask[i]) {
        TString method = fUEMethods[i];
        AliAnalysisTaskMTFPID* task = fJetUEPIDtask[i];
        
        TList *jetUElist = 0x0;
        TList* mcJetUElist = 0x0;

        if (mcJetContainer) {
          if (fOnlyLeadingJets)
            mcJetContainer->SortArray();
          
          if (method.Contains("RC",TString::kIgnoreCase) || method.Contains("Random",TString::kIgnoreCase)) {
            mcJetUElist = GetUEJetsWithRandomConeMethod(mcJetContainer, TMath::Abs(GetFFRadius()), task->GetEtaAbsCutUp()); 
          }
          else if (method.Contains("PC",TString::kIgnoreCase) || method.Contains("Perpendicular",TString::kIgnoreCase)) {
            mcJetUElist = GetUEJetsWithPerpendicularConeMethod(mcJetContainer);
          }
        } 
            
        if (jetContainer) {
        //Get Underlying event jets
          if (fOnlyLeadingJets)
            jetContainer->SortArray();
          //Random Cones______________________________________________
          if (method.Contains("RC",TString::kIgnoreCase) || method.Contains("Random",TString::kIgnoreCase)) {
            jetUElist = GetUEJetsWithRandomConeMethod(jetContainer, TMath::Abs(GetFFRadius()), task->GetEtaAbsCutUp());  
          }
          //Perpendicular Cones______________________________________________
          else if (method.Contains("PC",TString::kIgnoreCase) || method.Contains("Perpendicular",TString::kIgnoreCase)) {
            jetUElist = GetUEJetsWithPerpendicularConeMethod(jetContainer);
          }
        }
        //Process Tracks
        if (jetUElist) {
          for (Int_t i=0;i<jetUElist->GetEntries();++i) { 
            AliEmcalJet* jet = (AliEmcalJet*)(jetUElist->At(i));
            Double_t UEPtDensity = 0.0;
            Float_t jetPt = jet->Pt();
            if (jetContainer->GetRhoParameter())
              jetPt = jetPt - jetContainer->GetRhoVal() * jet->Area();
            task->FillRecJets(centPercent, jetPt);
            task->FillJetArea(centPercent, jet->Area());
              
            TList *jetUEtracklist = GetTracksInCone(jet, trackContainer);
            if (!jetUEtracklist) {
              cout << "No track list for an underlying event jet. Check Code!" << endl;
              continue;
            }
            
            jetUEtracklist->SetOwner(kFALSE);

            for (Int_t j=0;j<jetUEtracklist->GetEntries();++j) {
              AliVTrack* UEtrack = (AliVTrack*)(jetUEtracklist->At(j));
              
              if (UEtrack) {        
                Bool_t survivedTPCCutMIGeo = AliAnalysisTaskMTFPID::TPCCutMIGeo(UEtrack, InputEvent());
                Bool_t survivedTPCnclCut = AliAnalysisTaskMTFPID::TPCnclCut(UEtrack);   //Included above
                Double_t dEdxTPC = tuneOnDataTPC ? pidResponse->GetTPCsignalTunedOnData(UEtrack) : UEtrack->GetTPCsignal();
                
                if ((!task->GetUseTPCCutMIGeo() || survivedTPCCutMIGeo) && (!task->GetUseTPCnclCut() || survivedTPCnclCut) && task->IsInAcceptedEtaRange(TMath::Abs(UEtrack->Eta())) && (dEdxTPC > 0.0)) {
                  AnalyseJetTrack(UEtrack, jet, task, 0x0, centPercent, mcParticleContainer);
//                   task->FillDCA(UEtrack, dEdxTPC, primVtx, fMCEvent, jetPt);
                }
                UEPtDensity += UEtrack->Pt();
              }
            }
            delete jetUEtracklist;
            jetUEtracklist = 0x0;
            UEPtDensity = UEPtDensity/(TMath::Abs(GetFFRadius()) * TMath::Abs(GetFFRadius()) * TMath::Pi());
            task->FillUEDensity(centPercent, UEPtDensity);
          }
          delete jetUElist;
          jetUElist = 0x0;
        }
        if (mcJetUElist) {
          for (Int_t i=0;i<mcJetUElist->GetEntries();++i) {
            AliEmcalJet* jet = (AliEmcalJet*)(mcJetUElist->At(i));
            
            task->FillGenJets(centPercent, jet->Pt());
            TList* mcJetUEtracklist = GetTracksInCone(jet, mcParticleContainer);
            if (!mcJetUEtracklist) {
              cout << "No track list for MC Tracks in Underlying event. Check Code!" << endl;
              continue;
            }
            mcJetUEtracklist->SetOwner(kFALSE);
            for (Int_t j=0;j<mcJetUEtracklist->GetEntries();++j) {
              AliVParticle* mcUEParticle = dynamic_cast<AliVParticle*>(mcJetUEtracklist->At(j));
              PerformJetMonteCarloAnalysisGeneratedYield(jet, mcUEParticle, task, centPercent);
            }
          }
        }
      }
    }
  }
  //End of underlying event        
  
  if (pidTaskWithV0Indexes)
    pidTaskWithV0Indexes->ClearV0PIDlist();
  
  //___________________
  
  AliDebugStream(2) << "Processing done, posting output data now..." << std::endl;
  
  //Post output data.
  PostData(1, fCommonHistList);
  
  if (fUseInclusivePIDtask) {
    for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
      fInclusivePIDtask[i]->PostOutputData();
    }
  }

  if (fUseJetPIDtask) {
    for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
      fJetPIDtask[i]->PostOutputData();
    }
  }

  if (fUseJetUEPIDtask) {
    for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
      fJetUEPIDtask[i]->PostOutputData();
    }
  }
  
  AliDebugStream(1) << "Event done!" << std::endl;
  
  return kTRUE;
}

//______________________________________________________________
void AliAnalysisTaskIDFragmentationFunction::Terminate(Option_t *) 
{
  AliDebugStream(1) << "Terminate!" << std::endl;
  // terminated
  if (fUseJetUEPIDtask) {
    for (Int_t i=0;i<fNumJetUEPIDtasks;++i) {
      fJetUEPIDtask[i]->NormalizeJetArea(TMath::Abs(GetFFRadius()));
    }
  }
}  

// _________________________________________________________________________________________________________
void AliAnalysisTaskIDFragmentationFunction::SetProperties(THnSparse* h, Int_t dim, const char** labels)
{
  // Set properties of THnSparse 

  for(Int_t i=0; i<dim; i++){
    h->GetAxis(i)->SetTitle(labels[i]);
    h->GetAxis(i)->SetTitleColor(1);
  }
}

// __________________________________________________________________________________________
void AliAnalysisTaskIDFragmentationFunction::SetProperties(TH1* h,const char* x, const char* y)
{
  //Set properties of histos (x and y title)

  h->SetXTitle(x);
  h->SetYTitle(y);
  h->GetXaxis()->SetTitleColor(1);
  h->GetYaxis()->SetTitleColor(1);
}

// _________________________________________________________________________________________________________
void AliAnalysisTaskIDFragmentationFunction::SetProperties(TH1* h,const char* x, const char* y, const char* z)
{
  //Set properties of histos (x,y and z title)

  h->SetXTitle(x);
  h->SetYTitle(y);
  h->SetZTitle(z);
  h->GetXaxis()->SetTitleColor(1);
  h->GetYaxis()->SetTitleColor(1);
  h->GetZaxis()->SetTitleColor(1);
}



//_______________________________________________________
Double_t AliAnalysisTaskIDFragmentationFunction::GetMCStrangenessFactor(Double_t pt) const
{
  // factor strangeness data/MC as function of pt from UE analysis (Sara Vallero)

  Double_t alpha = 1;

  if(0.150<pt && pt<0.200) alpha = 3.639;
  if(0.200<pt && pt<0.250) alpha = 2.097;
  if(0.250<pt && pt<0.300) alpha = 1.930;
  if(0.300<pt && pt<0.350) alpha = 1.932;
  if(0.350<pt && pt<0.400) alpha = 1.943;
  if(0.400<pt && pt<0.450) alpha = 1.993;
  if(0.450<pt && pt<0.500) alpha = 1.989;
  if(0.500<pt && pt<0.600) alpha = 1.963;
  if(0.600<pt && pt<0.700) alpha = 1.917;
  if(0.700<pt && pt<0.800) alpha = 1.861;
  if(0.800<pt && pt<0.900) alpha = 1.820;
  if(0.900<pt && pt<1.000) alpha = 1.741;
  if(1.000<pt && pt<1.500) alpha = 0.878;

  return alpha;
}

//__________________________________________________________________________________________________
Double_t AliAnalysisTaskIDFragmentationFunction::GetMCStrangenessFactorCMS(AliAODMCParticle* daughter,
                                                                           AliMCParticleContainer* mcParticleContainer) const
{
  // strangeness ratio MC/data as function of mother pt from CMS data in |eta|<2.0

  if(!mcParticleContainer) return 1;

  AliAODMCParticle* currentMother   = daughter;
  AliAODMCParticle* currentDaughter = daughter;


  // find first primary mother K0s, Lambda or Xi   
  while(1){

    Int_t daughterPDG   = currentDaughter->GetPdgCode();	

    Int_t motherLabel   = currentDaughter->GetMother();

    currentMother     = (AliAODMCParticle*) mcParticleContainer->GetMCParticleWithLabel(motherLabel);

    if(!currentMother){ 
      currentMother = currentDaughter; 
      break; 
    }

    Int_t motherPDG   = currentMother->GetPdgCode();	
 
    // phys. primary found ?  	
    if(currentMother->IsPhysicalPrimary()) break; 

    if(TMath::Abs(daughterPDG) == 321){ // K+/K- e.g. from phi (ref data not feeddown corrected)
      currentMother = currentDaughter; break; 
    }	 	
    if(TMath::Abs(motherPDG) == 310 ){ // K0s e.g. from phi (ref data not feeddown corrected)
      break; 
    } 	
    if(TMath::Abs(motherPDG) == 3212 && TMath::Abs(daughterPDG) == 3122){ // mother Sigma0, daughter Lambda (this case not included in feeddown corr.)
      currentMother = currentDaughter; break; 
    }

    currentDaughter = currentMother;
  }


  Int_t motherPDG   = currentMother->GetPdgCode();	
  Double_t motherGenPt = currentMother->Pt();	

  return AliAnalysisTaskMTFPID::GetMCStrangenessFactorCMS(motherPDG, motherGenPt);
}

//_______________________________________________________
Double_t AliAnalysisTaskIDFragmentationFunction::GetDistanceJetTrack(const AliEmcalJet* jet, const AliVParticle* track) const
{
  // Calculate the distance between jet and track in the eta-phi-plane.

  if (!fStudyTransversalJetStructure)
    return -1;  
  
  if (!jet || !track)
    return -1.;
  
  TVector3 v_track(track->Px(), track->Py(), track->Pz());
  TVector3 v_jet(jet->Px(), jet->Py(), jet->Pz());
  
  return v_track.DeltaR(v_jet);
}

//_______________________________________________________
Double_t AliAnalysisTaskIDFragmentationFunction::GetPerpendicularMomentumTrackJet(const AliEmcalJet* jet, const AliVParticle* track) const
{
  // Calculate the momentum of the track perpendicular to the jet axis.
  
  if (!fStudyTransversalJetStructure)
    return -1;
  
  if (!jet || !track)
    return -1.;
  
  TVector3 v_track(track->Px(), track->Py(), track->Pz());
  TVector3 v_jet(jet->Px(), jet->Py(), jet->Pz());
  
  // Important that the jet is the argument, otherwise one would get the transverse jet momentum w.r.t. the track direction
  return v_track.Perp(v_jet);
}


// _________________________________________________________________________________
Bool_t AliAnalysisTaskIDFragmentationFunction::IsSecondaryWithStrangeMotherMC(AliAODMCParticle* part, AliMCParticleContainer* mcParticleContainer)
{
  // Check whether particle is a secondary with strange mother, i.e. returns kTRUE if a strange mother is found
  // and the particle is NOT a physical primary. In all other cases kFALSE is returned
  
  if (!mcParticleContainer || !part)
    return kFALSE;
  
  if (part->IsPhysicalPrimary())
    return kFALSE;
  
  Int_t iMother = part->GetMother();
  if (iMother < 0)
    return kFALSE;
  
  
  AliAODMCParticle* partM = dynamic_cast<AliAODMCParticle*>(mcParticleContainer->GetMCParticleWithLabel(iMother));
  if (!partM) 
   return kFALSE;
  
  Int_t codeM = TMath::Abs(partM->GetPdgCode());
  Int_t mfl = Int_t(codeM / TMath::Power(10, Int_t(TMath::Log10(codeM))));
  if (mfl == 3 && codeM != 3) // codeM = 3 is for s quark
    return kTRUE;
  
  return kFALSE;
}

TList* AliAnalysisTaskIDFragmentationFunction::GetUEJetsWithRandomConeMethod(AliJetContainer* jetContainer, Double_t coneRadius, Double_t maxEtaTrack) {
  TList *jetUElist = new TList();
  jetUElist->SetOwner(kTRUE); 
  if (jetContainer) {
    AliEmcalJet* jetRC;
    for (auto jet : jetContainer->accepted()) {
      jetRC = GetRandomCone(jet, maxEtaTrack - coneRadius, 2.0 * coneRadius);
      if (jetRC) {
        jetUElist->Add(jetRC);
      }
      if (fOnlyLeadingJets)
        break;  
    }
  }
  return jetUElist;
}

TList* AliAnalysisTaskIDFragmentationFunction::GetUEJetsWithPerpendicularConeMethod(AliJetContainer* jetContainer) {
  TList *jetUElist = new TList();
  jetUElist->SetOwner(kTRUE); 
  if (jetContainer) {
    Double_t perpAngle = TMath::Pi()/2.0;
    for (auto jet : jetContainer->accepted()) {
      AliEmcalJet* jetPC = 0x0;
      jetPC = GetPerpendicularCone(jet,perpAngle);
      if (jetPC) {
        jetUElist->Add(jetPC);
      }
      jetPC = GetPerpendicularCone(jet,-perpAngle);
      if (jetPC) {
        jetUElist->Add(jetPC);
      }
      if (fOnlyLeadingJets)
        break;  
    }
  }
  return jetUElist;
}

AliEmcalJet* AliAnalysisTaskIDFragmentationFunction::GetRandomCone(AliEmcalJet* processedJet, Double_t dEtaConeMax, Double_t dDistance) const { 
  if (!processedJet)
    return 0x0;
  
  TLorentzVector vecRdCone;
  AliEmcalJet* jetRC = 0x0;                          //random cone candidate
  Double_t dEta, dPhi;                             //random eta and phi value for RC
  UInt_t i = 0;
  Double_t jetPt = processedJet->Pt();           //ATTENTION Do NOT subtract rho here - we want to have an exact copy of the processed jet. 
  Double_t jetM = processedJet->M();
  
  //Loop is running until the number of trials equals the maximum number of trials or the created random jet is not overlapping with any jet in the jetlist  
  do {
    delete jetRC;                                 //Delete previously created jet           
    dEta = dEtaConeMax * (2 * fRandom->Rndm() - 1.0);    //random eta value in range: [-dEtaConeMax,+dEtaConeMax]
    dPhi = TMath::TwoPi() * fRandom->Rndm();        //random phi value in range: [0,2*Pi]
    jetRC = new AliEmcalJet(jetPt, dEta, dPhi, jetM);             //new RC candidate
    ++i;          
  }
  while (i<fRCTrials && OverlapsWithAnyRecJet(jetRC, dDistance));   
  
  //If the last jet is also overlapping, delete it and set pointer to zero
  if(OverlapsWithAnyRecJet(jetRC, dDistance)) {
    delete jetRC;
    jetRC = 0x0;
  }
  
  if (GetUseRealJetArea() && jetRC)
    jetRC->SetArea(processedJet->Area());
  
  return jetRC;
}

AliEmcalJet* AliAnalysisTaskIDFragmentationFunction::GetPerpendicularCone(AliEmcalJet* processedJet, Double_t perpAngle) const {
  if (!processedJet)
    return 0x0;
  
  AliEmcalJet* perpJet = new AliEmcalJet(processedJet->Pt(), processedJet->Eta(), processedJet->Phi() + perpAngle, processedJet->M());
  
  if (GetUseRealJetArea())
    perpJet->SetArea(processedJet->Area());
  
  return perpJet;
}

Bool_t AliAnalysisTaskIDFragmentationFunction::OverlapsWithAnyRecJet(const AliVParticle* part, Double_t dDistance) const {
  
  if(!part) return kFALSE;
  
  AliJetContainer* jetContainer = GetJetContainer(GetNameJetContainer());
  
  if (dDistance < 0.0)
    dDistance = jetContainer->GetJetRadius();
  
  if (!jetContainer)
    return kFALSE;
  
  for(auto jet : jetContainer->accepted()){   //loop over all reconstructed jets in events      
    if(!jet){
      AliWarningStream() << "AliAnalysisTaskIDFragmentationFunction::OverlapsWithAnyRecJet jet pointer invalid!" << std::endl;
      continue;
    }
    if(IsParticleInCone(jet, part, dDistance) == kTRUE) return kTRUE;//RC and JC are overlapping
    
  }//end loop testing RC-JC overlap
  return kFALSE;//RC and JC are not overlapping -> good!
}

Bool_t AliAnalysisTaskIDFragmentationFunction::IsParticleInCone(const AliVParticle* part1, const AliVParticle* part2, Double_t dRMax) const {
// decides whether a particle is inside a jet cone, or has a minimum distance to a second track/axis
  if (!part1 || !part2)
    return kFALSE;

  TVector3 vecMom2(part2->Px(),part2->Py(),part2->Pz());
  TVector3 vecMom1(part1->Px(),part1->Py(),part1->Pz());
  Double_t dR = vecMom2.DeltaR(vecMom1); // = sqrt(dEta*dEta+dPhi*dPhi)
  return (dR < dRMax); // momentum vectors of part1 and part2 are closer than dRMax
}

TList* AliAnalysisTaskIDFragmentationFunction::GetTracksInCone(const AliEmcalJet* jet, AliParticleContainer* particleContainer) const {
  TList* jetTrackList = new TList();
  jetTrackList->SetOwner(kFALSE);
  if (!jet)
    return jetTrackList;
  
  if (!particleContainer)
    particleContainer = GetTrackContainer(GetNameTrackContainer());
  
  Double_t radius = TMath::Abs(GetFFRadius());
  
  if (GetUseRealJetArea() && jet->Area() > 0.0)
    radius = TMath::Sqrt(jet->Area()/TMath::Pi());
    
  
  for (auto track : particleContainer->accepted()) {

    if(!track)
      continue;

    if (IsParticleInCone(track,jet,radius)) {
      jetTrackList->Add(track);
    }
  }
  return jetTrackList;
}

void AliAnalysisTaskIDFragmentationFunction::PerformJetMonteCarloAnalysisGeneratedYield(AliEmcalJet* jet, AliVParticle* trackVP, AliAnalysisTaskMTFPID* task, Double_t centPercent, AliJetContainer* mcJetContainer) {
  if (!jet || !trackVP || !task)
    return;
  
  if (!mcJetContainer)
    mcJetContainer = GetJetContainer(GetNameMCParticleJetContainer());
  
  if (!mcJetContainer) {
    AliWarningStream() << "Monte-Carlo Jet Container not found." << std::endl;
  }
  
  Float_t jetPt   = jet->Pt();
  if (mcJetContainer->GetRhoParameter())
    jetPt = jetPt - mcJetContainer->GetRhoVal() * jet->Area();
  Float_t trackPt = trackVP->Pt();
  
  // Efficiency, jets - particle level
  AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(trackVP);
  if (!part) {
    AliErrorStream() << "expected ref track not found!" << std::endl;
    return;
  }
  
  Double_t trackEtaMin = -1.0 * (task->GetEtaAbsCutUp());
  Double_t trackEtaMax = task->GetEtaAbsCutUp();
  
  AliMCParticleContainer* mcParticleContainer = GetMCParticleContainer(GetNameTrackContainerEfficiency());
  if (mcParticleContainer) {
    trackEtaMin = mcParticleContainer->GetParticleEtaMin();
    trackEtaMax = mcParticleContainer->GetParticleEtaMax();
  }  
  
  // Fill efficiency for generated primaries and also fill histos for generated yields (primaries + all)
  if (part->Eta() > trackEtaMax || part->Eta() < trackEtaMin)
    return;
  
  if (!part->IsPhysicalPrimary())
    return;
  
  Int_t mcID = AliAnalysisTaskMTFPID::PDGtoMCID(part->GetPdgCode());

  // Following lines are not needed - just keep other species (like casecades) - will end up in overflow bin
  // and only affect the efficiencies for all (i.e. not identified) what is desired!
  //if (mcID == AliPID::kUnknown)
  //  continue;
  //
  //   Int_t iMother = part->GetMother();      
  //   if (iMother >= 0)
  //     continue; // Not a physical primary
  //
  
  Double_t z = -1., xi = -1.;
  AliAnalysisTaskMTFPID::GetJetTrackObservables(trackPt, jetPt, z, xi);
  Double_t distance = GetDistanceJetTrack(jet, trackVP);
  Double_t jT = GetPerpendicularMomentumTrackJet(jet, trackVP);
  
  // AliAODMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
  Double_t chargeMC = part->Charge() / 3.;
  
  if (TMath::Abs(chargeMC) < 0.01)
    return; // Reject neutral particles (only relevant, if mcID is not used)
  
  Double_t valuesGenYield[AliAnalysisTaskMTFPID::kGenYieldNumAxes] = { static_cast<Double_t>(mcID), trackPt, centPercent, jetPt, z, xi,
                                                                    chargeMC, distance, jT };

  if (task->IsInAcceptedEtaRange(TMath::Abs(part->Eta()))) {
    valuesGenYield[task->GetIndexOfChargeAxisGenYield()] = task->GetStoreCharge() ? chargeMC : -2;
    task->FillGeneratedYield(valuesGenYield);
  }
  
  // Always store the charge for the efficiency (needed for geant-fluka correction)
  Double_t valuesEff[AliAnalysisTaskMTFPID::kEffNumAxes] = { static_cast<Double_t>(mcID), trackPt, part->Eta(), chargeMC,
                                                          centPercent, jetPt, z, xi, distance, jT };                                                     
  task->FillEfficiencyContainer(valuesEff, AliAnalysisTaskMTFPID::kStepGenWithGenCuts);
}
  
void AliAnalysisTaskIDFragmentationFunction::AnalyseJetTrack(AliVTrack* track, AliEmcalJet* jet, AliAnalysisTaskMTFPID* task, const Bool_t* trackRejectedByTask, Double_t centPercent, AliMCParticleContainer* mcParticleContainer) {
  
  if (!track || (!task && !trackRejectedByTask)) {
    AliErrorStream() << "Cannot analyse track! Track: " << track << "; Task: " << task << "; RejectionArray: " << trackRejectedByTask << std::endl;
    return;
  }
    
  AliJetContainer* jetContainer = GetJetContainer(GetNameJetContainer());
    
  Float_t jetPt   = jet->Pt();
  if (jetContainer->GetRhoParameter())
    jetPt = jetPt - jetContainer->GetRhoVal() * jet->Area();    
    
  Double_t trackPt = track->Pt();
  
  Int_t label = TMath::Abs(track->GetLabel());
  
  Int_t pdg = 0;
  Int_t mcID = AliPID::kUnknown;
  
  Double_t z = -1., xi = -1.;
  AliAnalysisTaskMTFPID::GetJetTrackObservables(trackPt, jetPt, z, xi);
  Double_t distance = GetDistanceJetTrack(jet, track);
  Double_t jT = GetPerpendicularMomentumTrackJet(jet, track); 
  
  if (!mcParticleContainer)
    mcParticleContainer = GetMCParticleContainer(GetNameMCParticleContainer());

  AliAODMCParticle* gentrack = mcParticleContainer ? mcParticleContainer->GetMCParticleWithLabel(label) : 0x0;    
  
  //Get PDG Code
  if (gentrack) {
    pdg = gentrack->GetPdgCode();
  }
 
  if (task)
    task->ProcessTrack(track, pdg, centPercent, jetPt, kFALSE, kTRUE, distance, jT);
  else {
    for (Int_t i=0;i<fNumJetPIDtasks;++i) {
      if (!trackRejectedByTask[i])
        fJetPIDtask[i]->ProcessTrack(track, pdg, centPercent, jetPt, kFALSE, kTRUE, distance, jT);
    }
  }
  
  //Process 
  if (gentrack && !fUseFastSimulations) {
    // Secondaries, jets
    mcID = AliAnalysisTaskMTFPID::PDGtoMCID(pdg);
    
    Double_t valueRecAllCuts[AliAnalysisTaskMTFPID::kEffNumAxes] = { static_cast<Double_t>(mcID), trackPt, track->Eta(),
                                                                  static_cast<Double_t>(track->Charge()),
                                                                  centPercent, jetPt, z, xi, distance, jT };
                                                  
    Double_t weight = IsSecondaryWithStrangeMotherMC(gentrack, mcParticleContainer)
                        ? GetMCStrangenessFactorCMS(gentrack, mcParticleContainer) : 1.0;
                        
    if (task) {
      task->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskMTFPID::kStepRecWithRecCutsMeasuredObs);
      task->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskMTFPID::kStepRecWithRecCutsMeasuredObsStrangenessScaled,
                                            weight);
    }
    else {
      for (Int_t i=0;i<fNumJetPIDtasks;++i) {
        if (!trackRejectedByTask[i]) {
          fJetPIDtask[i]->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskMTFPID::kStepRecWithRecCutsMeasuredObs);
          fJetPIDtask[i]->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskMTFPID::kStepRecWithRecCutsMeasuredObsStrangenessScaled,
                                            weight);
        }
      }
    }                  
                        
    if (gentrack->IsPhysicalPrimary()) {
      Double_t genPt = gentrack->Pt();
      Double_t genZ = -1., genXi = -1.;
      AliAnalysisTaskMTFPID::GetJetTrackObservables(genPt, jetPt, genZ, genXi);
      Double_t genDistance = GetDistanceJetTrack(jet, gentrack);
      Double_t genJt = GetPerpendicularMomentumTrackJet(jet, gentrack);
      
      // Always store the charge for the efficiency (needed for geant-fluka correction)
      Double_t valueGenAllCuts[AliAnalysisTaskMTFPID::kEffNumAxes] = { static_cast<Double_t>(mcID), genPt, gentrack->Eta(), 
                                                                    gentrack->Charge() / 3., centPercent, jetPt, genZ, 
                                                                    genXi, genDistance, genJt };
      
      Double_t valuePtResolution[AliAnalysisTaskMTFPID::kPtResNumAxes] = { jetPt, genPt, trackPt, gentrack->Charge() / 3., centPercent }; 
      
      if (task) {
        task->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskMTFPID::kStepRecWithRecCutsMeasuredObsPrimaries);
        task->FillEfficiencyContainer(valueGenAllCuts, AliAnalysisTaskMTFPID::kStepRecWithRecCutsPrimaries);                                                     
        task->FillPtResolution(mcID, valuePtResolution); 
      }
      else {
        for (Int_t i=0;i<fNumJetPIDtasks;++i) {
          if (!trackRejectedByTask[i]) {
            fJetPIDtask[i]->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskMTFPID::kStepRecWithRecCutsMeasuredObsPrimaries);
            fJetPIDtask[i]->FillEfficiencyContainer(valueGenAllCuts, AliAnalysisTaskMTFPID::kStepRecWithRecCutsPrimaries);                                                     
            fJetPIDtask[i]->FillPtResolution(mcID, valuePtResolution); 
          }
        }
      }         
      
      //Efficiency, jets - detector level
      Double_t trackEtaMin, trackEtaMax;
      if (task) {
        trackEtaMin = -1.0 * (task->GetEtaAbsCutUp());
        trackEtaMax = task->GetEtaAbsCutUp();
      }
      else {
        trackEtaMin = -1.0 * (fJetPIDtask[0]->GetEtaAbsCutUp());
        trackEtaMax = fJetPIDtask[0]->GetEtaAbsCutUp();
      }        
      AliMCParticleContainer* trackContainerEfficiency = GetMCParticleContainer(GetNameTrackContainerEfficiency());
      if (trackContainerEfficiency) {
        trackEtaMin = trackContainerEfficiency->GetParticleEtaMin();
        trackEtaMax = trackContainerEfficiency->GetParticleEtaMax();
      } 
      
      if (gentrack->Eta() <= trackEtaMax || gentrack->Eta() >= trackEtaMin) {
        
        Double_t genZ = -1., genXi = -1.;
        Double_t genPt = gentrack->Pt();
        AliAnalysisTaskMTFPID::GetJetTrackObservables(genPt, jetPt, genZ, genXi);
        Double_t genDistance = GetDistanceJetTrack(jet, gentrack);
        Double_t genJt = GetPerpendicularMomentumTrackJet(jet, gentrack);
        
        Double_t measZ = -1., measXi = -1.;
        AliAnalysisTaskMTFPID::GetJetTrackObservables(trackPt, jetPt, measZ, measXi);
        Double_t measDistance = GetDistanceJetTrack(jet, track);
        Double_t measJt = GetPerpendicularMomentumTrackJet(jet, track);
        
        // AliAODMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
        // Always store the charge for the efficiency (needed for geant-fluka correction)
        Double_t value[AliAnalysisTaskMTFPID::kEffNumAxes] = { static_cast<Double_t>(mcID), genPt, gentrack->Eta(), gentrack->Charge() / 3., centPercent, jetPt, genZ, genXi, genDistance, genJt };                                                     
        
        Double_t valueMeas[AliAnalysisTaskMTFPID::kEffNumAxes] = { static_cast<Double_t>(mcID), trackPt, track->Eta(),
                                                                static_cast<Double_t>(track->Charge()),
                                                                centPercent, jetPt, measZ, measXi, measDistance, measJt }; 
        
        if (task) {
          task->FillEfficiencyContainer(value, AliAnalysisTaskMTFPID::kStepRecWithGenCuts);
          task->FillEfficiencyContainer(valueMeas, AliAnalysisTaskMTFPID::kStepRecWithGenCutsMeasuredObs);                                                   
        }
        else {
          for (Int_t i=0;i<fNumJetPIDtasks;++i) {
            if (!trackRejectedByTask[i]) {
              fJetPIDtask[i]->FillEfficiencyContainer(value, AliAnalysisTaskMTFPID::kStepRecWithGenCuts);
              fJetPIDtask[i]->FillEfficiencyContainer(valueMeas, AliAnalysisTaskMTFPID::kStepRecWithGenCutsMeasuredObs);       
            }
          }
        }             
      } 
    }                                                                        
  }    
  return;
}  

void AliAnalysisTaskIDFragmentationFunction::FillDCA(AliVTrack* track, AliMCParticleContainer* mcParticleContainer) {
  
  if (!fFillDCA)
    return;
  
  Double_t trackPt = track->Pt();
  
  Int_t label = TMath::Abs(track->GetLabel());

  if (!mcParticleContainer)
    mcParticleContainer = GetMCParticleContainer(GetNameMCParticleContainer());
  
  Int_t mcID = AliPID::kUnknown;

  AliAODMCParticle* gentrack = mcParticleContainer ? mcParticleContainer->GetMCParticleWithLabel(label) : 0x0;   
  
  if (gentrack) {
    Int_t pdg = gentrack->GetPdgCode();
    mcID = AliAnalysisTaskMTFPID::PDGtoMCID(pdg);
  }
  
  Double_t dca[2] = {0., 0.}; // 0: xy; 1: z
  
  Double_t v[3]   = {0, };
  Double_t pos[3] = {0, };
  const AliVVertex* primVtx = InputEvent()->GetPrimaryVertex();
  if (!primVtx) {
    AliWarningStream() << "Primary Vertex not found, do not fill DCA" << std::endl;
    return;
  }
  primVtx->GetXYZ(v);
  track->GetXYZ(pos);
  dca[0] = pos[0] - v[0];
  dca[1] = pos[1] - v[1];
  
  // "Unidentified" for data and MC
  fhDCA_XY->Fill(trackPt, dca[0]);
  fhDCA_Z->Fill(trackPt, dca[1]);
  
  // "Identified" for MC
  if (gentrack && mcID != AliPID::kUnknown) {
    // MC
    if (gentrack->IsPhysicalPrimary()) {
      fhDCA_XY_prim_MCID[mcID]->Fill(trackPt, dca[0]);
      fhDCA_Z_prim_MCID[mcID]->Fill(trackPt, dca[1]);
    }
    else {
      fhDCA_XY_sec_MCID[mcID]->Fill(trackPt, dca[0]);
      fhDCA_Z_sec_MCID[mcID]->Fill(trackPt, dca[1]);
    }
  } 
}

void AliAnalysisTaskIDFragmentationFunction::SetUpFastJetWrapperWithOriginalValues(AliFJWrapper* wrapper) {
  wrapper->SetAreaType(fastjet::active_area_explicit_ghosts);
  wrapper->SetGhostArea(0.005);
  wrapper->SetR(TMath::Abs(GetFFRadius()));
  //Currently not working, including AliEmcalJetTask fails
  //wrapper->SetAlgorithm(AliEmcalJetTask::ConvertToFJAlgo(jetContainer->GetJetAlgorithm()));
  //wrapper->SetRecombScheme(AliEmcalJetTask::ConvertToFJRecoScheme(jetContainer->GetRecombinationScheme()));
  wrapper->SetAlgorithm(fastjet::antikt_algorithm);
  wrapper->SetRecombScheme(fastjet::pt_scheme);
  wrapper->SetMaxRap(1);    
}

void AliAnalysisTaskIDFragmentationFunction::FillEfficiencyContainerFromTrack(AliAODMCParticle* part, AliEmcalJet* jet, Double_t centPercent, AliAnalysisTaskMTFPID::EffSteps step) {
	
	if (!part || !jet)
		return;
	
	Double_t jetPt = jet->Pt();
	Double_t trackPt = part->Pt();

	Double_t z = -1., xi = -1.;
	AliAnalysisTaskMTFPID::GetJetTrackObservables(trackPt, jetPt, z, xi);
	Double_t distance = GetDistanceJetTrack(jet, part);
	Double_t jT = GetPerpendicularMomentumTrackJet(jet, part);   
	
	Double_t values[AliAnalysisTaskMTFPID::kEffNumAxes] = { static_cast<Double_t>(AliAnalysisTaskMTFPID::PDGtoMCID(part->GetPdgCode())), trackPt, part->Eta(), part->Charge() / 3., centPercent, jetPt, z, xi, distance, jT}; 
	
	for (Int_t i=0;i<fNumJetPIDtasks;++i) {
		fJetPIDtask[i]->FillEfficiencyContainer(values, step);
	}	
}

void AliAnalysisTaskIDFragmentationFunction::FillPIDTasksCutHisto(Double_t value, AliAnalysisTaskMTFPID::CutHistoType histoType) {
  if (fUseInclusivePIDtask) {
    for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
      fInclusivePIDtask[i]->FillCutHisto(value, histoType);
    }
  }    
  
  if (fUseJetPIDtask) {
    for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
      fJetPIDtask[i]->FillCutHisto(value, histoType);
    }
  }
  
  if (fUseJetUEPIDtask) {
    for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
      fJetUEPIDtask[i]->FillCutHisto(value, histoType);
    }
  }
}

void AliAnalysisTaskIDFragmentationFunction::IncrementPIDTasksEventCounts(Double_t centPercent, AliAnalysisTaskMTFPID::EventCounterType eventCounterType, Bool_t* isPileUpInclusivePIDtask, Bool_t* isPileUpJetPIDtask, Bool_t* isPileUpJetUEPIDtask) {
  // Count events with trigger selection, note: Set centrality percentile fix to -1 for pp for PID framework
  if (fUseInclusivePIDtask) {
    for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
      if (isPileUpInclusivePIDtask && isPileUpInclusivePIDtask[i])
          continue;
      
      fInclusivePIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, eventCounterType);
    }
  }

  if (fUseJetPIDtask) {
    for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
      if (isPileUpJetPIDtask && isPileUpJetPIDtask[i])
        continue;
      
      fJetPIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, eventCounterType);
    }
  }
  
  if (fUseJetUEPIDtask) {
    for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
      if (isPileUpJetUEPIDtask && isPileUpJetUEPIDtask[i])
        continue;
      
      fJetUEPIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, eventCounterType);
    }
  }  
}
