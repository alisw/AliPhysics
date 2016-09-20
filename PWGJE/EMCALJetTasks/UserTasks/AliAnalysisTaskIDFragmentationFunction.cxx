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

#include "AliAnalysisTaskPID.h"
#include "AliPIDResponse.h"


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
   ,fCentralityEstimator("V0M")
   ,fNameTrackContainer("Tracks")
   ,fNameTrackContainerEfficiency("")
   ,fNameMCParticleContainer("")
   ,fNameJetContainer("")
   ,fNameMCParticleJetContainer("")
   ,fUseAODInputJets(kTRUE)
   ,fUsePhysicsSelection(kTRUE)
   ,fEvtSelectionMask(0)
   ,fEventClass(0)
   ,fMaxVertexZ(10)
   ,fFFRadius(0)
   ,fFFMinLTrackPt(-1)
   ,fFFMaxTrackPt(-1)
   ,fFFMinnTracks(0)   
   ,fAvgTrials(0)
   ,fStoreXi(0)
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
   ,fh1nGenJets(0)
   ,fhDCA_XY(0)
   ,fhDCA_Z(0)
   ,fhJetPtRefMultEta5(0)
   ,fhJetPtRefMultEta8(0)
   ,fhJetPtMultPercent(0)

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
   ,fRCTrials(1)
   ,fUEMethods(0x0)
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
  ,fCentralityEstimator("V0M")
  ,fNameTrackContainer("Tracks")
  ,fNameTrackContainerEfficiency("")
  ,fNameMCParticleContainer("")
  ,fNameJetContainer("")
  ,fNameMCParticleJetContainer("")  
  ,fUseAODInputJets(kTRUE)
  ,fUsePhysicsSelection(kTRUE)
  ,fEvtSelectionMask(0)
  ,fEventClass(0)
  ,fMaxVertexZ(10)
  ,fFFRadius(0)
  ,fFFMinLTrackPt(-1)
  ,fFFMaxTrackPt(-1)
  ,fFFMinnTracks(0)  
  ,fAvgTrials(0)
  ,fStoreXi(0)  
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
  ,fh1nGenJets(0)
  ,fhDCA_XY(0)
  ,fhDCA_Z(0)
  ,fhJetPtRefMultEta5(0)
  ,fhJetPtRefMultEta8(0)
  ,fhJetPtMultPercent(0)
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
  ,fRCTrials(1)
  ,fUEMethods(0x0)
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

  if(fRandom)               delete fRandom;
  
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
      Error("Notify","No current file");
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
  
    if(!fh1Xsec||!fh1Trials){
      Printf("%s:%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
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
  fRandom->SetSeed(0);

  return kTRUE;
}

//__________________________________________________________________
void AliAnalysisTaskIDFragmentationFunction::UserCreateOutputObjects()
{
  // create output objects

  if(fDebug > 1) Printf("AliAnalysisTaskIDFragmentationFunction::UserCreateOutputObjects()");

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
  fh1EvtMult 	             = new TH1F("fh1EvtMult","Event multiplicity, track pT cut > 150 MeV/c, |#eta| < 0.9",120,0.,12000.);
  fh1EvtCent 	             = new TH1F("fh1EvtCent","centrality",100,0.,100.);

  fh1Xsec                    = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fh1Trials                  = new TH1F("fh1Trials","trials from pyxsec.root",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fh1PtHard                  = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",350,-.5,349.5);
  fh1PtHardTrials            = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",350,-.5,349.5);
  
  fh1EvtsPtHardCut           = new TH1F("fh1EvtsPtHardCut", "#events before and after MC #it{p}_{T,hard} cut;;Events",2,0,2);
  fh1EvtsPtHardCut->GetXaxis()->SetBinLabel(1, "All");
  fh1EvtsPtHardCut->GetXaxis()->SetBinLabel(2, "#it{p}_{T,hard}");
  

  fh1nRecJetsCuts            = new TH1F("fh1nRecJetsCuts","reconstructed jets per event",10,-0.5,9.5);
  fh1nGenJets                = new TH1F("fh1nGenJets","generated jets per event",10,-0.5,9.5);
  
 
  // ____________ define histograms ___________________

  fCommonHistList->Add(fh1EvtSelection);
  fCommonHistList->Add(fh1VtxSelection);
  fCommonHistList->Add(fh1EvtMult);
  fCommonHistList->Add(fh1EvtCent);
  fCommonHistList->Add(fh1VertexNContributors);
  fCommonHistList->Add(fh1VertexZ);    
  fCommonHistList->Add(fh1nRecJetsCuts);
  fCommonHistList->Add(fh1Xsec);
  fCommonHistList->Add(fh1Trials);
  fCommonHistList->Add(fh1PtHard);
  fCommonHistList->Add(fh1PtHardTrials);
  fCommonHistList->Add(fh1EvtsPtHardCut);
 
  if(GetJetContainer(GetNameMCParticleJetContainer())) fCommonHistList->Add(fh1nGenJets);
  
  // Default analysis utils
  fAnaUtils = new AliAnalysisUtils();
  
  // Not used yet, but to be save, forward vertex z cut to analysis utils object
  fAnaUtils->SetMaxVtxZ(fMaxVertexZ);

  // Load PID framework if desired
  if(fDebug > 1) Printf("AliAnalysisTaskIDFragmentationFunction::UserCreateOutputObjects() -> Loading PID framework");
  
  if (fUseJetPIDtask || fUseInclusivePIDtask || fUseJetUEPIDtask) {
    TObjArray* tasks = AliAnalysisManager::GetAnalysisManager()->GetTasks();
    if (!tasks) {
      Printf("ERROR loading PID tasks: Failed to retrieve tasks from analysis manager!\n");
      
      fUseInclusivePIDtask = kFALSE;
      fUseJetPIDtask = kFALSE;
      fUseJetUEPIDtask = kFALSE;
    }
    
    if (fUseInclusivePIDtask) {
      delete [] fInclusivePIDtask;
      fInclusivePIDtask = 0x0;
      
      if (fNumInclusivePIDtasks > 0) {
        fInclusivePIDtask = new AliAnalysisTaskPID*[fNumInclusivePIDtasks];
        
        for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
          fInclusivePIDtask[i] = (AliAnalysisTaskPID*)tasks->FindObject(fNameInclusivePIDtask[i].Data());
          
          if (!fInclusivePIDtask[i]) {
            Printf("ERROR: Failed to load inclusive pid task!\n");
            fUseInclusivePIDtask = kFALSE;
          }
        }
      }
      else {
        Printf("WARNING: zero inclusive pid tasks!\n");
        fUseInclusivePIDtask = kFALSE;
      }
    }    
    
    if (fUseJetPIDtask) {
      delete [] fJetPIDtask;
      fJetPIDtask = 0x0;
      
      if (fNumJetPIDtasks > 0) {
        fJetPIDtask = new AliAnalysisTaskPID*[fNumJetPIDtasks];
        
        for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
          fJetPIDtask[i] = (AliAnalysisTaskPID*)tasks->FindObject(fNameJetPIDtask[i].Data());
          
          if (!fJetPIDtask[i]) {
            Printf("ERROR: Failed to load jet pid task!\n");
            fUseJetPIDtask = kFALSE;
          }
        }
      }
      else {
        Printf("WARNING: zero jet pid tasks!\n");
        fUseJetPIDtask = kFALSE;
      }
    }
    
    if (fUseJetUEPIDtask) {
      delete [] fJetUEPIDtask;
      fJetUEPIDtask = 0x0;
      
      if (fNumJetUEPIDtasks > 0) {
        fJetUEPIDtask = new AliAnalysisTaskPID*[fNumJetUEPIDtasks];
        
        for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
          fJetUEPIDtask[i] = (AliAnalysisTaskPID*)tasks->FindObject(fNameJetUEPIDtask[i].Data());
          
          if (!fJetUEPIDtask[i]) {
            Printf("ERROR: Failed to load jet underlying event pid task!\n");
            fUseJetUEPIDtask = kFALSE;
          }
        }
      }
      else {
        Printf("WARNING: zero jet underlying event pid tasks!\n");
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

  if(fDebug > 2) Printf("AliAnalysisTaskIDFragmentationFunction::UserCreateOutputObjects() -> Posting Output");
  
  PostData(1, fCommonHistList);
  
  if(fDebug > 2) Printf("AliAnalysisTaskIDFragmentationFunction::UserCreateOutputObjects() -> Done");
}

//_______________________________________________
void AliAnalysisTaskIDFragmentationFunction::Init()
{
  // Initialization
  if(fDebug > 1) Printf("AliAnalysisTaskIDFragmentationFunction::Init()");

}

//_____________________________________________________________
Bool_t AliAnalysisTaskIDFragmentationFunction::FillHistograms() 
{
  
  if(fDebug > 1) Printf("AliAnalysisTaskIDFragmentationFunction::FillHistograms()");
  
  
  if(fDebug > 1) Printf("Analysis event #%5d", (Int_t) fEntry);
  
  fMCEvent = MCEvent();
  if(!fMCEvent){
    if(fDebug>3) Printf("%s:%d MCEvent not found in the input", (char*)__FILE__,__LINE__);
  }
   
  
  // Extract pThard and nTrials in case of MC. 
  
  Double_t ptHard = 0.;
  Double_t nTrials = 1; // trials for MC trigger weight for real data
  Bool_t pythiaGenHeaderFound = kFALSE;

  if(fMCEvent){
    AliGenEventHeader* genHeader = fMCEvent->GenEventHeader();
    
    if(genHeader){
      AliGenPythiaEventHeader*  pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
      AliGenHijingEventHeader*  hijingGenHeader = 0x0;
      
      if(pythiaGenHeader){
        if(fDebug>3) Printf("%s:%d pythiaGenHeader found", (char*)__FILE__,__LINE__);
        pythiaGenHeaderFound = kTRUE;
        nTrials = pythiaGenHeader->Trials();
        ptHard  = pythiaGenHeader->GetPtHard();
      } else { // no pythia, hijing?
        
        if(fDebug>3) Printf("%s:%d no pythiaGenHeader found", (char*)__FILE__,__LINE__);
        
        hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);
        if(!hijingGenHeader){
          Printf("%s:%d no pythiaGenHeader or hjingGenHeader found", (char*)__FILE__,__LINE__);
        } else {
          if(fDebug>3) Printf("%s:%d hijingGenHeader found", (char*)__FILE__,__LINE__);
        }
      }
      
      //fh1Trials->Fill("#sum{ntrials}",fAvgTrials); 
    }
  }
  
  
  // Cut on pThard if fMCEvent and pThard >= 0 and fill histo with #evt before and after the cut
  if (fMCEvent) {
    // Before cut
    fh1EvtsPtHardCut->Fill(0.); 
    
    if (fUseInclusivePIDtask) {
      for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
        fInclusivePIDtask[i]->FillCutHisto(0., AliAnalysisTaskPID::kMCPtHardCut);
      }
    }    
    
    if (fUseJetPIDtask) {
      for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
        fJetPIDtask[i]->FillCutHisto(0., AliAnalysisTaskPID::kMCPtHardCut);
      }
    }
    
    if (fUseJetUEPIDtask) {
      for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
        fJetUEPIDtask[i]->FillCutHisto(0., AliAnalysisTaskPID::kMCPtHardCut);
      }
    }
    
    // Cut
    if (fMCPtHardCut >= 0. && ptHard >= fMCPtHardCut) {
      if (fDebug>3) Printf("%s:%d skipping event with pThard %f (>= %f)", (char*)__FILE__,__LINE__, ptHard, fMCPtHardCut);
      PostData(1, fCommonHistList);
      return kFALSE;
    }
    
    // After cut
    fh1EvtsPtHardCut->Fill(1.);

    if (fUseInclusivePIDtask) {
      for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
        fInclusivePIDtask[i]->FillCutHisto(1., AliAnalysisTaskPID::kMCPtHardCut);
      }
    }
    
    if (fUseJetPIDtask) {
      for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
        fJetPIDtask[i]->FillCutHisto(1., AliAnalysisTaskPID::kMCPtHardCut);
      }
    }
    
    if (fUseJetUEPIDtask) {
      for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
        fJetUEPIDtask[i]->FillCutHisto(1., AliAnalysisTaskPID::kMCPtHardCut);
      }
    }    
 
  }
  
  // Trigger selection
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  
  if(!(inputHandler->IsEventSelected() & fEvtSelectionMask)){
    fh1EvtSelection->Fill(1.);
    if (fDebug > 1 ) Printf(" Trigger Selection: event REJECTED ... ");
    PostData(1, fCommonHistList);
    return kFALSE;
  }
  
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fESD){
    if(fDebug>3) Printf("%s:%d ESDEvent not found in the input", (char*)__FILE__,__LINE__);
  }
  
  // get AOD event from input/ouput
  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if( handler && handler->InheritsFrom("AliAODInputHandler") ) {
    fAOD  =  ((AliAODInputHandler*)handler)->GetEvent();
    if(fUseAODInputJets) fAODJets = fAOD;
    if (fDebug > 1)  Printf("%s:%d AOD event from input", (char*)__FILE__,__LINE__);
  }
  else {
    handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( handler && handler->InheritsFrom("AliAODHandler") ) {
      fAOD = ((AliAODHandler*)handler)->GetAOD();
      fAODJets = fAOD;
      if (fDebug > 1)  Printf("%s:%d AOD event from output", (char*)__FILE__,__LINE__);
    }
  }
  
  if(!fAODJets && !fUseAODInputJets){ // case we have AOD in input & output and want jets from output
    TObject* outHandler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( outHandler && outHandler->InheritsFrom("AliAODHandler") ) {
      fAODJets = ((AliAODHandler*)outHandler)->GetAOD();
      if (fDebug > 1)  Printf("%s:%d jets from output AOD", (char*)__FILE__,__LINE__);
    }
  }
  
  if(fNonStdFile.Length()!=0){
    // case we have an AOD extension - fetch the jets from the extended output
    
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
    fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);    
    if(!fAODExtension){
      if(fDebug>1)Printf("AODExtension not found for %s",fNonStdFile.Data());
    }
  }
  
  if(!fAOD){
    Printf("%s:%d AODEvent not found", (char*)__FILE__,__LINE__);
    return kFALSE;
  }
  if(!fAODJets){
    Printf("%s:%d AODEvent with jet branch not found", (char*)__FILE__,__LINE__);
    return kFALSE;
  }

  
  // event selection **************************************************
  // *** event class ***
  AliVEvent* evtForCentDetermination = handler->InheritsFrom("AliAODInputHandler") ? fAOD : InputEvent();
  
  Double_t centPercent = -1;
  
  if(fEventClass>0){
    Int_t cl = 0;
    if(handler->InheritsFrom("AliAODInputHandler")){ 
      // since it is not supported by the helper task define own classes
      centPercent = ((AliAODHeader*)fAOD->GetHeader())->GetCentrality();
      cl = 1;
      if(centPercent>10) cl = 2;
      if(centPercent>30) cl = 3;
      if(centPercent>50) cl = 4;
    }
    else {
      cl = AliAnalysisHelperJetTasks::EventClass();
      if(fESD) centPercent = fESD->GetCentrality()->GetCentralityPercentile(fCentralityEstimator.Data()); 
    }
    
    if(cl!=fEventClass){
      // event not in selected event class, reject event
      if (fDebug > 1) Printf("%s:%d event not in selected event class: event REJECTED ...",(char*)__FILE__,__LINE__);
      fh1EvtSelection->Fill(2.);
      PostData(1, fCommonHistList);
      return kFALSE;
    }
  }
  
  if (fCentralityEstimator.Contains("NoCentrality",TString::kIgnoreCase)) {
    centPercent = -1;
  }
  else {
    if (!fIsPP) centPercent = evtForCentDetermination->GetCentrality()->GetCentralityPercentile(fCentralityEstimator.Data());
  }
  
  AliPIDResponse *fPIDResponse = inputHandler->GetPIDResponse();
  if (!fPIDResponse) {
    AliError("PIDResponse object was not created");
  }
  
  fPIDResponse->SetCurrentCentrality(centPercent);

  // Retrieve reference multiplicities in |eta|<0.8 and <0.5
  const Int_t refMult5 = ((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb05();
  const Int_t refMult8 = ((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08();
  const Double_t centPercentPP = fAnaUtils->GetMultiplicityPercentile(fAOD, "V0M");
  
  
  // Count events with trigger selection, note: Set centrality percentile fix to -1 for pp for PID framework
  if (fUseInclusivePIDtask) {
    for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
      fInclusivePIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, AliAnalysisTaskPID::kTriggerSel);
    }
  }

  if (fUseJetPIDtask) {
    for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
      fJetPIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, AliAnalysisTaskPID::kTriggerSel);
    }
  }
  
  if (fUseJetUEPIDtask) {
    for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
      fJetUEPIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, AliAnalysisTaskPID::kTriggerSel);
    }
  }  

  // *** vertex cut ***
  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
	if (!primVtx) {
		Printf("%s:%d Primary vertex not found", (char*)__FILE__,__LINE__);
		return kFALSE;
	}
	
  Int_t nTracksPrim = primVtx->GetNContributors();
  fh1VertexNContributors->Fill(nTracksPrim);
  
  
  if (fDebug > 1) Printf("%s:%d primary vertex selection: %d", (char*)__FILE__,__LINE__,nTracksPrim);
  if(nTracksPrim <= 0) {
    if (fDebug > 1) Printf("%s:%d primary vertex selection: event REJECTED...",(char*)__FILE__,__LINE__); 
    fh1EvtSelection->Fill(3.);
    PostData(1, fCommonHistList);
    return kFALSE;
  }
  
  TString primVtxName(primVtx->GetName());

  if(primVtxName.CompareTo("TPCVertex",TString::kIgnoreCase) == 1){
    if (fDebug > 1) Printf("%s:%d primary vertex selection: TPC vertex, event REJECTED...",(char*)__FILE__,__LINE__);
    fh1EvtSelection->Fill(5.);
    PostData(1, fCommonHistList);
    return kFALSE;
  }
  
  // Count events with trigger selection and vtx cut, note: Set centrality percentile fix to -1 for pp for PID framework
  if (fUseInclusivePIDtask) {
    for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
      fInclusivePIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, AliAnalysisTaskPID::kTriggerSelAndVtxCut);
    }
  }
  
  if (fUseJetPIDtask) {
    for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
      fJetPIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, AliAnalysisTaskPID::kTriggerSelAndVtxCut);
    }
  }
  
  if (fUseJetUEPIDtask) {
    for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
      fJetUEPIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, AliAnalysisTaskPID::kTriggerSelAndVtxCut);
    }
  }  

  
  fh1VertexZ->Fill(primVtx->GetZ());
  
  if(TMath::Abs(primVtx->GetZ())>fMaxVertexZ){
    if (fDebug > 1) Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVtx->GetZ()); 
    fh1EvtSelection->Fill(4.);
    PostData(1, fCommonHistList);
    return kFALSE; 
  }
  
  // Count events with trigger selection, vtx cut and z vtx cut, note: Set centrality percentile fix to -1 for pp for PID framework
  if (fUseInclusivePIDtask) {
    for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
      fInclusivePIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, AliAnalysisTaskPID::kTriggerSelAndVtxCutAndZvtxCutNoPileUpRejection);
    }
  }
  
  if (fUseJetPIDtask) {
    for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
      fJetPIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, AliAnalysisTaskPID::kTriggerSelAndVtxCutAndZvtxCutNoPileUpRejection);
    }
  }
  
  if (fUseJetUEPIDtask) {
    for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
      fJetUEPIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, AliAnalysisTaskPID::kTriggerSelAndVtxCutAndZvtxCutNoPileUpRejection);
    }
  }  
  
  
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
  
  // Count events with trigger selection, vtx cut, z vtx cut and after pile-up rejection (if enabled in that task)
  // Note: Set centrality percentile fix to -1 for pp for PID framework
  if (fUseInclusivePIDtask) {
    for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
      isPileUpInclusivePIDtask[i] = fInclusivePIDtask[i]->GetIsPileUp(fAOD, fInclusivePIDtask[i]->GetPileUpRejectionType());
      if (!isPileUpInclusivePIDtask[i])
        fInclusivePIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, AliAnalysisTaskPID::kTriggerSelAndVtxCutAndZvtxCut);
      
      isPileUpForAllInclusivePIDTasks = isPileUpForAllInclusivePIDTasks && isPileUpInclusivePIDtask[i];
    }
  }
  
  if (fUseJetPIDtask) {
    for (Int_t i = 0; i < fNumJetPIDtasks; i++) {
      isPileUpJetPIDtask[i] = fJetPIDtask[i]->GetIsPileUp(fAOD, fJetPIDtask[i]->GetPileUpRejectionType());
      if (!isPileUpJetPIDtask[i])
        fJetPIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, AliAnalysisTaskPID::kTriggerSelAndVtxCutAndZvtxCut);
      
      isPileUpForAllJetPIDTasks = isPileUpForAllJetPIDTasks && isPileUpJetPIDtask[i];
    }
  }
  
  if (fUseJetUEPIDtask) {
    for (Int_t i = 0; i < fNumJetUEPIDtasks; i++) {
      isPileUpJetUEPIDtask[i] = fJetUEPIDtask[i]->GetIsPileUp(fAOD, fJetUEPIDtask[i]->GetPileUpRejectionType());
      if (!isPileUpJetUEPIDtask[i])
        fJetUEPIDtask[i]->IncrementEventCounter(fIsPP ? -1. : centPercent, AliAnalysisTaskPID::kTriggerSelAndVtxCutAndZvtxCut);
      
      isPileUpForAllJetUEPIDTasks = isPileUpForAllJetUEPIDTasks && isPileUpJetUEPIDtask[i];
    }
  }
    
  
  
  if (fDebug > 1) Printf("%s:%d event ACCEPTED ...",(char*)__FILE__,__LINE__); 
  fh1EvtSelection->Fill(0.);
  fh1VtxSelection->Fill(primVtx->GetType());
  fh1EvtCent->Fill(centPercent);

  // Set centrality percentile fix to -1 for pp to be used for the PID framework
  if (fIsPP)
    centPercent = -1;
  
  // Call ConfigureTaskForCurrentEvent of PID tasks to ensure that everything is set up properly for the current event
  // (e.g. run/period dependence of max eta variation map)
  if (fUseInclusivePIDtask) {
    for (Int_t i = 0; i < fNumInclusivePIDtasks; i++)
      if (!isPileUpInclusivePIDtask[i])
        fInclusivePIDtask[i]->ConfigureTaskForCurrentEvent(fAOD);
  }
  
  if (fUseJetPIDtask) {
    for (Int_t i = 0; i < fNumJetPIDtasks; i++)
      if (!isPileUpJetPIDtask[i])
        fJetPIDtask[i]->ConfigureTaskForCurrentEvent(fAOD);
  }
  
  if (fUseJetUEPIDtask) {
    for (Int_t i = 0; i < fNumJetUEPIDtasks; i++)
      if (!isPileUpJetUEPIDtask[i])
        fJetUEPIDtask[i]->ConfigureTaskForCurrentEvent(fAOD);
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
  
  AliPIDResponse* pidResponse = 0x0;
  Bool_t tuneOnDataTPC = kFALSE;
  if (fUseJetPIDtask || fUseInclusivePIDtask || fUseJetUEPIDtask) {
    if (!inputHandler) {
      AliFatal("Input handler needed");
      return kFALSE;
    }
    else {
      // PID response object
      pidResponse = inputHandler->GetPIDResponse();
      if (!pidResponse) {
        AliFatal("PIDResponse object was not created");
        return kFALSE;
      }
      else {
        tuneOnDataTPC = pidResponse->IsTunedOnData() &&
                        ((pidResponse->GetTunedOnDataMask() & AliPIDResponse::kDetTPC) == AliPIDResponse::kDetTPC);
      }
    }
  }
  
  if(fDebug>2)Printf("%s:%d Starting processing...",(char*)__FILE__,__LINE__);
  
  //____ analysis, fill histos ___________________________________________________
  
  AliTrackContainer *trackContainer = GetTrackContainer(GetNameTrackContainer());
  AliMCParticleContainer *mcParticleContainer = GetMCParticleContainer(GetNameMCParticleContainer());
  Double_t trackEtaMin = -0.9;
  Double_t trackEtaMax = 0.9;
  if (trackContainer) {
    trackEtaMin = trackContainer->GetParticleEtaMin();
    trackEtaMax = trackContainer->GetParticleEtaMax();
  }
      
  // Fill efficiency for generated primaries and also fill histos for generated yields (primaries + all)
  // Efficiency, inclusive - particle level
  if(fDebug>2)Printf("%s:%d Starting Inclusive Efficiency...",(char*)__FILE__,__LINE__);
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
      
      Int_t mcID = AliAnalysisTaskPID::PDGtoMCID(part->GetPdgCode());
      
      // Following lines are not needed - just keep other species (like casecades) - will end up in overflow bin
      // and only affect the efficiencies for all (i.e. not identified) what is desired!
      //if (mcID == AliPID::kUnknown)
      //  continue;
      
      if (!part->IsPhysicalPrimary())
        continue;
      /*
      Int_t iMother = part->GetMother();
      if (iMother >= 0)
        continue; // Not a physical primary
      */

      Double_t pT = part->Pt();
      
      // AliAODMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
      Double_t chargeMC = part->Charge() / 3.;
      
      if (TMath::Abs(chargeMC) < 0.01)
        continue; // Reject neutral particles (only relevant, if mcID is not used)
      
      Double_t valuesGenYield[AliAnalysisTaskPID::kGenYieldNumAxes] = { static_cast<Double_t>(mcID), pT, centPercent,
                                                                        -1, -1, -1, -1, -1, -1 };

      for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
        if (!isPileUpInclusivePIDtask[i] && fInclusivePIDtask[i]->IsInAcceptedEtaRange(TMath::Abs(part->Eta()))) {
          valuesGenYield[fInclusivePIDtask[i]->GetIndexOfChargeAxisGenYield()] = fInclusivePIDtask[i]->GetStoreCharge() ? chargeMC : -2;
          fInclusivePIDtask[i]->FillGeneratedYield(valuesGenYield);
        }
      }
      
      // Always store the charge for the efficiency (needed for geant-fluka correction)
      Double_t valuesEff[AliAnalysisTaskPID::kEffNumAxes] = { static_cast<Double_t>(mcID), pT, part->Eta(), chargeMC,
                                                              centPercent, -1, -1, -1, -1, -1 };// no jet pT etc since inclusive spectrum 
      for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
        if (!isPileUpInclusivePIDtask[i])
          fInclusivePIDtask[i]->FillEfficiencyContainer(valuesEff, AliAnalysisTaskPID::kStepGenWithGenCuts);
      }
    }
  }
  
  AliTrackContainer* trackContainerEfficiency = GetTrackContainer(GetNameTrackContainerEfficiency());
  if (fUseInclusivePIDtask && mcParticleContainer && !isPileUpForAllInclusivePIDTasks) {
    //Efficiency, inclusive - detector level
    for(auto inclusiveaod : trackContainerEfficiency->accepted()) {
      // fill inclusive tracks XXX, they have the same track cuts!

      if(inclusiveaod){
        Double_t dEdxTPC = tuneOnDataTPC ? pidResponse->GetTPCsignalTunedOnData(inclusiveaod) 
                                          : inclusiveaod->GetTPCsignal();
        
        if (dEdxTPC <= 0)
          continue;
        
        Bool_t survivedTPCCutMIGeo = AliAnalysisTaskPID::TPCCutMIGeo(inclusiveaod, InputEvent());
        Bool_t survivedTPCnclCut = AliAnalysisTaskPID::TPCnclCut(inclusiveaod);
        
        Int_t label = TMath::Abs(inclusiveaod->GetLabel());
        
        // find MC track in our list, if available
        AliAODMCParticle* gentrack = mcParticleContainer->GetMCParticleWithLabel(label);
        Int_t pdg = 0;
        
        if (gentrack)
          pdg = gentrack->GetPdgCode();
        
        // For efficiency: Reconstructed track has survived all cuts on the detector level (no cut on eta acceptance)
        // and has an associated MC track
        // -> Check whether associated MC track belongs to the clean MC sample defined above,
        //    i.e. survives the particle level track cuts
        if (gentrack) {
          Int_t mcID = AliAnalysisTaskPID::PDGtoMCID(pdg);
          
          // Following lines are not needed - just keep other species (like casecades) - will end up in overflow bin
          // and only affect the efficiencies for all (i.e. not identified) what is desired!
          //if (mcID == AliPID::kUnknown)
          //  continue;
          
          // Fill efficiency for reconstructed primaries
          if (!gentrack->IsPhysicalPrimary())
            continue;
          /*
            *     Int_t iMother = gentrack->GetMother();
            *     if (iMother >= 0)
            *       continue; // Not a physical primary
            */
            
          if (gentrack->Eta() > trackEtaMax || gentrack->Eta() < trackEtaMin)
            continue;
            
          // AliAODMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
          // Always store the charge for the efficiency (needed for geant-fluka correction)
          Double_t value[AliAnalysisTaskPID::kEffNumAxes] = { static_cast<Double_t>(mcID), gentrack->Pt(), gentrack->Eta(),
                                                              gentrack->Charge() / 3., centPercent, -1, -1, -1, -1, -1 };// no jet pT etc since inclusive spectrum 
          for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
            if (!isPileUpInclusivePIDtask[i] &&
                ((!fInclusivePIDtask[i]->GetUseTPCCutMIGeo() && !fInclusivePIDtask[i]->GetUseTPCnclCut()) ||
                (survivedTPCCutMIGeo && fInclusivePIDtask[i]->GetUseTPCCutMIGeo()) ||
                (survivedTPCnclCut && fInclusivePIDtask[i]->GetUseTPCnclCut())))
              fInclusivePIDtask[i]->FillEfficiencyContainer(value, AliAnalysisTaskPID::kStepRecWithGenCuts);
          }
              
          Double_t valueMeas[AliAnalysisTaskPID::kEffNumAxes] = { static_cast<Double_t>(mcID), inclusiveaod->Pt(), inclusiveaod->Eta(),
                                                                  static_cast<Double_t>(inclusiveaod->Charge()), centPercent,
                                                                  -1, -1, -1, -1, -1 };// no jet pT etc since inclusive spectrum 
          for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
            if (!isPileUpInclusivePIDtask[i] &&
                ((!fInclusivePIDtask[i]->GetUseTPCCutMIGeo() && !fInclusivePIDtask[i]->GetUseTPCnclCut()) ||
                (survivedTPCCutMIGeo && fInclusivePIDtask[i]->GetUseTPCCutMIGeo()) ||
                (survivedTPCnclCut && fInclusivePIDtask[i]->GetUseTPCnclCut())))
              fInclusivePIDtask[i]->FillEfficiencyContainer(valueMeas, AliAnalysisTaskPID::kStepRecWithGenCutsMeasuredObs);
          }
        }
      }
    }
  }  
  
  if(fDebug>2)Printf("%s:%d Process inclusive tracks...",(char*)__FILE__,__LINE__);
  if (fUseInclusivePIDtask && trackContainer && !isPileUpForAllInclusivePIDTasks) {
    for(auto part : trackContainer->accepted()) {
      if (!part) 
        continue;

      Double_t dEdxTPC = tuneOnDataTPC ? pidResponse->GetTPCsignalTunedOnData(part)
                                        : part->GetTPCsignal();
      
      if (dEdxTPC <= 0)
        continue;
      
      Bool_t survivedTPCCutMIGeo = AliAnalysisTaskPID::TPCCutMIGeo(part, InputEvent());
      Bool_t survivedTPCnclCut = AliAnalysisTaskPID::TPCnclCut(part);
      
      Int_t label = TMath::Abs(part->GetLabel());

      // find MC track in our list, if available
      AliAODMCParticle* gentrack = mcParticleContainer ? mcParticleContainer->GetMCParticleWithLabel(label) : 0x0;
      Int_t pdg = 0;
      
      if (gentrack)
        pdg = gentrack->GetPdgCode();
      
      for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
        if (!isPileUpInclusivePIDtask[i] &&
            ((!fInclusivePIDtask[i]->GetUseTPCCutMIGeo() && !fInclusivePIDtask[i]->GetUseTPCnclCut()) ||
            (survivedTPCCutMIGeo && fInclusivePIDtask[i]->GetUseTPCCutMIGeo()) ||
            (survivedTPCnclCut && fInclusivePIDtask[i]->GetUseTPCnclCut()))) {
              if (fInclusivePIDtask[i]->IsInAcceptedEtaRange(TMath::Abs(part->Eta())))
                fInclusivePIDtask[i]->ProcessTrack(part, pdg, centPercent, -1, kFALSE, kTRUE, fStoreXi, -1, -1); // no jet pT etc since inclusive spectrum 
        }
      }
      
      if (gentrack) {
        Int_t mcID = AliAnalysisTaskPID::PDGtoMCID(pdg);
        // Always store the charge for the efficiency (needed for geant-fluka correction)
        Double_t valueRecAllCuts[AliAnalysisTaskPID::kEffNumAxes] = { static_cast<Double_t>(mcID), part->Pt(), part->Eta(), 
                                                                      static_cast<Double_t>(part->Charge()), centPercent,
                                                                      -1, -1, -1, -1, -1 };
        for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
          if (!isPileUpInclusivePIDtask[i] &&
              ((!fInclusivePIDtask[i]->GetUseTPCCutMIGeo() && !fInclusivePIDtask[i]->GetUseTPCnclCut()) ||
              (survivedTPCCutMIGeo && fInclusivePIDtask[i]->GetUseTPCCutMIGeo()) ||
              (survivedTPCnclCut && fInclusivePIDtask[i]->GetUseTPCnclCut())))
            fInclusivePIDtask[i]->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskPID::kStepRecWithRecCutsMeasuredObs);
        }
        
        Double_t weight = IsSecondaryWithStrangeMotherMC(gentrack, mcParticleContainer)
                            ? GetMCStrangenessFactorCMS(gentrack, mcParticleContainer) : 1.0;
        for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
          if (!isPileUpInclusivePIDtask[i] &&
              ((!fInclusivePIDtask[i]->GetUseTPCCutMIGeo() && !fInclusivePIDtask[i]->GetUseTPCnclCut()) ||
              (survivedTPCCutMIGeo && fInclusivePIDtask[i]->GetUseTPCCutMIGeo()) ||
              (survivedTPCnclCut && fInclusivePIDtask[i]->GetUseTPCnclCut())))
            fInclusivePIDtask[i]->FillEfficiencyContainer(valueRecAllCuts, 
                                                          AliAnalysisTaskPID::kStepRecWithRecCutsMeasuredObsStrangenessScaled,
                                                          weight);
        }
        
        if (gentrack->IsPhysicalPrimary()) {
          // AliAODMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
          // Always store the charge for the efficiency (needed for geant-fluka correction)
          Double_t valueGenAllCuts[AliAnalysisTaskPID::kEffNumAxes] = { static_cast<Double_t>(mcID), gentrack->Pt(), gentrack->Eta(), 
                                                                        gentrack->Charge() / 3., centPercent, -1, -1, -1, -1, -1 };
          
          Double_t valuePtResolution[AliAnalysisTaskPID::kPtResNumAxes] = { -1, gentrack->Pt(), part->Pt(),
                                                                            gentrack->Charge() / 3., centPercent };
        
          for (Int_t i = 0; i < fNumInclusivePIDtasks; i++) {
            if (!isPileUpInclusivePIDtask[i] &&
                ((!fInclusivePIDtask[i]->GetUseTPCCutMIGeo() && !fInclusivePIDtask[i]->GetUseTPCnclCut()) ||
                (survivedTPCCutMIGeo && fInclusivePIDtask[i]->GetUseTPCCutMIGeo()) ||
                (survivedTPCnclCut && fInclusivePIDtask[i]->GetUseTPCnclCut()))) {
              fInclusivePIDtask[i]->FillEfficiencyContainer(valueRecAllCuts, 
                                                            AliAnalysisTaskPID::kStepRecWithRecCutsMeasuredObsPrimaries);
              fInclusivePIDtask[i]->FillEfficiencyContainer(valueGenAllCuts, 
                                                            AliAnalysisTaskPID::kStepRecWithRecCutsPrimaries);
              
              fInclusivePIDtask[i]->FillPtResolution(mcID, valuePtResolution);
            }
          }
        }
      }
    }    
  }
  
  if(fDebug>2)Printf("%s:%d Process Jets - efficiency, particle level..",(char*)__FILE__,__LINE__);
  AliJetContainer* mcJetContainer = GetJetContainer(GetNameMCParticleJetContainer());
  if (fOnlyLeadingJets)
    mcJetContainer->SortArray();
  
  if (fUseJetPIDtask && mcJetContainer && !isPileUpForAllJetPIDTasks) {
    
    for (Int_t i=0;i<fNumJetPIDtasks;i++) {
      if (!isPileUpJetPIDtask[i]) {

        for(auto jet : mcJetContainer->accepted()) {
            
          if(!jet) 
            continue;
          
          Float_t jetPt   = jet->Pt();
          if (mcJetContainer->GetRhoParameter())
            jetPt = jetPt - mcJetContainer->GetRhoVal() * jet->Area();          
          
          fJetPIDtask[i]->FillGenJets(centPercent, jetPt);
        
          TClonesArray *arrayWithTracks = mcParticleContainer ? mcParticleContainer->GetArray() : 0x0;
  
          for(Int_t it=0; it<jet->GetNumberOfTracks(); ++it) {
        
            AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(jet->TrackAt(it,arrayWithTracks));
            PerformJetMonteCarloAnalysisGeneratedYield(jet, trackVP, fJetPIDtask[i], centPercent);
          }
          
          if (fOnlyLeadingJets)
            break;
        }
      }
    }
  }
  
  if(fDebug>2)Printf("%s:%d Process jets...",(char*)__FILE__,__LINE__);
  AliJetContainer* jetContainer = GetJetContainer(GetNameJetContainer());
  if (fOnlyLeadingJets)
    jetContainer->SortArray();
  
  if (fUseJetPIDtask && jetContainer && !isPileUpForAllJetPIDTasks) {
          
    for (auto jet : jetContainer->accepted()) {
      if (!jet)
        continue;

      Float_t jetPt   = jet->Pt();
      if (jetContainer->GetRhoParameter())
        jetPt = jetPt - jetContainer->GetRhoVal() * jet->Area();
    
      for (Int_t i=0;i<fNumJetPIDtasks;++i) {
        if (!isPileUpJetPIDtask[i])
          fJetPIDtask[i]->FillRecJets(centPercent, jetPt);
      }
      
      fhJetPtRefMultEta5->Fill(refMult5, jetPt);
      fhJetPtRefMultEta8->Fill(refMult8, jetPt);
      fhJetPtMultPercent->Fill(centPercentPP, jetPt);
      
      TClonesArray *arrayWithTracks = trackContainer->GetArray();

      for(Int_t j=0; j<jet->GetNumberOfTracks(); ++j) {

        AliVTrack *track  = dynamic_cast<AliVTrack*>(jet->TrackAt(j,arrayWithTracks));
        Bool_t trackRejectedByTask[arrSizeJet];
        Bool_t trackRejectedByAllTasks = kTRUE;
        if (track) {
          Bool_t survivedTPCCutMIGeo = AliAnalysisTaskPID::TPCCutMIGeo(track, InputEvent());
          Bool_t survivedTPCnclCut = AliAnalysisTaskPID::TPCnclCut(track); 
          Double_t dEdxTPC = tuneOnDataTPC ? pidResponse->GetTPCsignalTunedOnData(track) : track->GetTPCsignal();
          for (Int_t i=0;i<fNumJetPIDtasks;++i) {
            trackRejectedByTask[i] = isPileUpJetPIDtask[i] || !((!fJetPIDtask[i]->GetUseTPCCutMIGeo() && !fJetPIDtask[i]->GetUseTPCnclCut()) || (survivedTPCCutMIGeo && fJetPIDtask[i]->GetUseTPCCutMIGeo()) || (survivedTPCnclCut && fJetPIDtask[i]->GetUseTPCnclCut()) && fJetPIDtask[i]->IsInAcceptedEtaRange(TMath::Abs(track->Eta())) && (dEdxTPC > 0.0));
            trackRejectedByAllTasks = trackRejectedByAllTasks && trackRejectedByTask[i];
          }
        
          if (!trackRejectedByAllTasks) {
            AnalyseJetTrack(track, jet, 0x0, trackRejectedByTask, centPercent, mcParticleContainer);
            FillDCA(track, mcParticleContainer);
          }
        }
        
        if (fOnlyLeadingJets)
          break;   
      }
    }
  }
  
  if(fDebug>2)Printf("%s:%d Process Underlying event...",(char*)__FILE__,__LINE__);
  if (fUseJetUEPIDtask && jetContainer && !isPileUpForAllJetUEPIDTasks) {
    for (Int_t i=0;i<fNumJetUEPIDtasks;++i) {
      if (!isPileUpJetUEPIDtask[i]) {
        TString method = fUEMethods[i];
        AliAnalysisTaskPID* task = fJetUEPIDtask[i];
        
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
                Bool_t survivedTPCCutMIGeo = AliAnalysisTaskPID::TPCCutMIGeo(UEtrack, InputEvent());
                Bool_t survivedTPCnclCut = AliAnalysisTaskPID::TPCnclCut(UEtrack); 
                Double_t dEdxTPC = tuneOnDataTPC ? pidResponse->GetTPCsignalTunedOnData(UEtrack) : UEtrack->GetTPCsignal();
                
                if ((!task->GetUseTPCCutMIGeo() && !task->GetUseTPCnclCut()) || (survivedTPCCutMIGeo && task->GetUseTPCCutMIGeo()) || (survivedTPCnclCut && task->GetUseTPCnclCut()) && task->IsInAcceptedEtaRange(TMath::Abs(UEtrack->Eta())) && (dEdxTPC > 0.0))
                  AnalyseJetTrack(UEtrack, jet, task, 0x0, centPercent, mcParticleContainer);
                
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
  
  //___________________
  
  if(fDebug > 2) Printf("%s:%d Processing done, posting output data now...",(char*)__FILE__,__LINE__);
  
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
  
  if(fDebug > 2) Printf("%s:%d Done",(char*)__FILE__,__LINE__);
  
  return kTRUE;
}

//______________________________________________________________
void AliAnalysisTaskIDFragmentationFunction::Terminate(Option_t *) 
{
  // terminated
  if (fUseJetUEPIDtask) {
    for (Int_t i=0;i<fNumJetUEPIDtasks;++i) {
      fJetUEPIDtask[i]->NormalizeJetArea(TMath::Abs(GetFFRadius()));
    }
  }

  if(fDebug > 1) printf("AliAnalysisTaskIDFragmentationFunction::Terminate() \n");
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
Double_t AliAnalysisTaskIDFragmentationFunction::GetDistanceJetTrack(const AliEmcalJet* jet, const AliVParticle* track,
                                                                     Bool_t useInternalFlag) const
{
  // Calculate the distance between jet and track in the eta-phi-plane. If useInternalFlag is on and fStudyTransversalJetStructure is kFALSE,
  // the function just returns -1
  
  if (useInternalFlag && !fStudyTransversalJetStructure)
    return -1;
  
  if (!jet || !track)
    return -1.;
  
  TVector3 v_track(track->Px(), track->Py(), track->Pz());
  TVector3 v_jet(jet->Px(), jet->Py(), jet->Pz());
  
  return v_track.DeltaR(v_jet);
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

  return AliAnalysisTaskPID::GetMCStrangenessFactorCMS(motherPDG, motherGenPt);
}


//_______________________________________________________
Double_t AliAnalysisTaskIDFragmentationFunction::GetPerpendicularMomentumTrackJet(const AliEmcalJet* jet, const AliVParticle* track,
                                                                                  Bool_t useInternalFlag) const
{
  // Calculate the momentum of the track perpendicular to the jet axis. If useInternalFlag is on and fStudyTransversalJetStructure is kFALSE,
  // the function just returns -1
  
  if (useInternalFlag && !fStudyTransversalJetStructure)
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
  while (i<fRCTrials && IsRCJCOverlap(jetRC,dDistance));   
  
  //If the last jet is also overlapping, delete it and set pointer to zero
  if(IsRCJCOverlap(jetRC,dDistance)) {
    delete jetRC;
    jetRC = 0x0;
  }
  
  if (jetRC)
    jetRC->SetArea(processedJet->Area());
  
  return jetRC;
}

AliEmcalJet* AliAnalysisTaskIDFragmentationFunction::GetPerpendicularCone(AliEmcalJet* processedJet, Double_t perpAngle) const {
  if (!processedJet)
    return 0x0;
  
  Double_t dDistance = 2.0 * TMath::Abs(GetFFRadius());
  
  AliEmcalJet* perpJet = new AliEmcalJet(processedJet->Pt(), processedJet->Eta(), processedJet->Phi()+perpAngle, processedJet->M());
  
  perpJet->SetArea(processedJet->Area());
  
  if (!IsParticleInCone(processedJet,perpJet,dDistance))
    return perpJet;
  else
    return 0x0;
}

Bool_t AliAnalysisTaskIDFragmentationFunction::IsRCJCOverlap(const AliVParticle* part, Double_t dDistance) const {
  
  if(!part) return kFALSE;
  if(!dDistance) return kFALSE;
  
  AliJetContainer* jetContainer = GetJetContainer(GetNameJetContainer());
  
  for(auto jet : jetContainer->accepted()){   //loop over all reconstructed jets in events      
    if(!jet){
      if (fDebug>2) std::cout << "AliAnalysisTaskJetChem::IsRCJCOverlap jet pointer invalid!" << std::endl; 
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
  if(dR < dRMax) // momentum vectors of part1 and part2 are closer than dRMax
    return kTRUE;
  return kFALSE;
}

TList* AliAnalysisTaskIDFragmentationFunction::GetTracksInCone(const AliEmcalJet* jet, AliParticleContainer* particleContainer) const {
  TList* jetTrackList = new TList();
  jetTrackList->SetOwner(kFALSE);
  if (!jet)
    return jetTrackList;
  
  if (!particleContainer)
    particleContainer = GetTrackContainer(GetNameTrackContainer());
  
  for (auto track : particleContainer->accepted()) {

    if(!track)
      continue;

    if (IsParticleInCone(track,jet,TMath::Abs(GetFFRadius()))) {
      jetTrackList->Add(track);
    }
  }
  return jetTrackList;
}

//End of underlying event calculations

void AliAnalysisTaskIDFragmentationFunction::PerformJetMonteCarloAnalysisGeneratedYield(AliEmcalJet* jet, AliVParticle* trackVP, AliAnalysisTaskPID* task, Double_t centPercent, AliJetContainer* mcJetContainer) {
  if (!jet || !trackVP || !task)
    return;
  
  if (!mcJetContainer)
    mcJetContainer = GetJetContainer(GetNameMCParticleJetContainer());
  
  if (!mcJetContainer) {
    std::cout << "Monte-Carlo Jet Container not found." << std::endl;
  }
  
  Float_t jetPt   = jet->Pt();
  if (mcJetContainer->GetRhoParameter())
    jetPt = jetPt - mcJetContainer->GetRhoVal() * jet->Area();
  Float_t trackPt = trackVP->Pt();
  
  // Efficiency, jets - particle level
  AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(trackVP);
  if (!part) {
    AliError("expected ref track not found ");
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
  
  Int_t mcID = AliAnalysisTaskPID::PDGtoMCID(part->GetPdgCode());
  
  // Following lines are not needed - just keep other species (like casecades) - will end up in overflow bin
  // and only affect the efficiencies for all (i.e. not identified) what is desired!
  //if (mcID == AliPID::kUnknown)
  //  continue;
  
  if (!part->IsPhysicalPrimary())
    return;
  //
  //   Int_t iMother = part->GetMother();      
  //   if (iMother >= 0)
  //     continue; // Not a physical primary
  //
  
  Double_t z = -1., xi = -1.;
  AliAnalysisTaskPID::GetJetTrackObservables(trackPt, jetPt, z, xi, fStoreXi);
  Double_t distance = GetDistanceJetTrack(jet, trackVP);
  Double_t jT = GetPerpendicularMomentumTrackJet(jet, trackVP);
  
  // AliAODMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
  Double_t chargeMC = part->Charge() / 3.;
  
  if (TMath::Abs(chargeMC) < 0.01)
    return; // Reject neutral particles (only relevant, if mcID is not used)
  
  Double_t valuesGenYield[AliAnalysisTaskPID::kGenYieldNumAxes] = { static_cast<Double_t>(mcID), trackPt, centPercent, jetPt, z, xi,
                                                                    chargeMC, distance, jT };

  if (task->IsInAcceptedEtaRange(TMath::Abs(part->Eta()))) {
    valuesGenYield[task->GetIndexOfChargeAxisGenYield()] = task->GetStoreCharge() ? chargeMC : -2;
    task->FillGeneratedYield(valuesGenYield);
  }
  
  // Always store the charge for the efficiency (needed for geant-fluka correction)
  Double_t valuesEff[AliAnalysisTaskPID::kEffNumAxes] = { static_cast<Double_t>(mcID), trackPt, part->Eta(), chargeMC,
                                                          centPercent, jetPt, z, xi, distance, jT };
  task->FillEfficiencyContainer(valuesEff, AliAnalysisTaskPID::kStepGenWithGenCuts);
}
  
void AliAnalysisTaskIDFragmentationFunction::AnalyseJetTrack(AliVTrack* track, AliEmcalJet* jet, AliAnalysisTaskPID* task, const Bool_t* trackRejectedByTask, Double_t centPercent, AliMCParticleContainer* mcParticleContainer) {
  
  if (!track || (!task && !trackRejectedByTask)) {
    std::cout << "Cannot analyse track" << std::endl;
    std::cout << "Track: " << track << std::endl;
    std::cout << "Task: " << task << std::endl;
    std::cout << "Rejection Array: " << trackRejectedByTask << std::endl;
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
  AliAnalysisTaskPID::GetJetTrackObservables(trackPt, jetPt, z, xi, fStoreXi);
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
    task->ProcessTrack(track, pdg, centPercent, jetPt, kFALSE, kTRUE, fStoreXi, distance, jT);
  else {
    for (Int_t i=0;i<fNumJetPIDtasks;++i) {
      if (!trackRejectedByTask[i])
        fJetPIDtask[i]->ProcessTrack(track, pdg, centPercent, jetPt, kFALSE, kTRUE, fStoreXi, distance, jT);
    }
  }
  
  //Process 
  if (gentrack) {
    
    // Secondaries, jets
    mcID = AliAnalysisTaskPID::PDGtoMCID(pdg);
    
    Double_t valueRecAllCuts[AliAnalysisTaskPID::kEffNumAxes] = { static_cast<Double_t>(mcID), trackPt, track->Eta(),
                                                                  static_cast<Double_t>(track->Charge()),
                                                                  centPercent, jetPt, z, xi, distance, jT };
                                                  
    Double_t weight = IsSecondaryWithStrangeMotherMC(gentrack, mcParticleContainer)
                        ? GetMCStrangenessFactorCMS(gentrack, mcParticleContainer) : 1.0;
                        
    if (task) {
      task->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskPID::kStepRecWithRecCutsMeasuredObs);
      task->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskPID::kStepRecWithRecCutsMeasuredObsStrangenessScaled,
                                            weight);
    }
    else {
      for (Int_t i=0;i<fNumJetPIDtasks;++i) {
        if (!trackRejectedByTask[i])
          fJetPIDtask[i]->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskPID::kStepRecWithRecCutsMeasuredObs);
          fJetPIDtask[i]->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskPID::kStepRecWithRecCutsMeasuredObsStrangenessScaled,
                                            weight);
      }
    }                  
                        
    if (gentrack->IsPhysicalPrimary()) {
      Double_t genPt = gentrack->Pt();
      Double_t genZ = -1., genXi = -1.;
      AliAnalysisTaskPID::GetJetTrackObservables(genPt, jetPt, genZ, genXi, fStoreXi);
      Double_t genDistance = GetDistanceJetTrack(jet, gentrack);
      Double_t genJt = GetPerpendicularMomentumTrackJet(jet, gentrack);
      
      // Always store the charge for the efficiency (needed for geant-fluka correction)
      Double_t valueGenAllCuts[AliAnalysisTaskPID::kEffNumAxes] = { static_cast<Double_t>(mcID), genPt, gentrack->Eta(), 
                                                                    gentrack->Charge() / 3., centPercent, jetPt, genZ, 
                                                                    genXi, genDistance, genJt };
      
      Double_t valuePtResolution[AliAnalysisTaskPID::kPtResNumAxes] = { jetPt, genPt, trackPt, gentrack->Charge() / 3., centPercent }; 
      
      if (task) {
        task->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskPID::kStepRecWithRecCutsMeasuredObsPrimaries);
        task->FillEfficiencyContainer(valueGenAllCuts, AliAnalysisTaskPID::kStepRecWithRecCutsPrimaries);                                                     
        task->FillPtResolution(mcID, valuePtResolution); 
      }
      else {
        for (Int_t i=0;i<fNumJetPIDtasks;++i) {
          if (!trackRejectedByTask[i])
            fJetPIDtask[i]->FillEfficiencyContainer(valueRecAllCuts, AliAnalysisTaskPID::kStepRecWithRecCutsMeasuredObsPrimaries);
            fJetPIDtask[i]->FillEfficiencyContainer(valueGenAllCuts, AliAnalysisTaskPID::kStepRecWithRecCutsPrimaries);                                                     
            fJetPIDtask[i]->FillPtResolution(mcID, valuePtResolution); 
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
      AliMCParticleContainer* mcParticleContainer = GetMCParticleContainer(GetNameTrackContainerEfficiency());
      if (mcParticleContainer) {
        trackEtaMin = mcParticleContainer->GetParticleEtaMin();
        trackEtaMax = mcParticleContainer->GetParticleEtaMax();
      } 
      
      if (gentrack->Eta() <= trackEtaMax || gentrack->Eta() >= trackEtaMin) {
        
        Double_t genZ = -1., genXi = -1.;
        Double_t genPt = gentrack->Pt();
        AliAnalysisTaskPID::GetJetTrackObservables(genPt, jetPt, genZ, genXi, fStoreXi);
        Double_t genDistance = GetDistanceJetTrack(jet, gentrack);
        Double_t genJt = GetPerpendicularMomentumTrackJet(jet, gentrack);
        
        Double_t measZ = -1., measXi = -1.;
        AliAnalysisTaskPID::GetJetTrackObservables(trackPt, jetPt, measZ, measXi, fStoreXi);
        Double_t measDistance = GetDistanceJetTrack(jet, track);
        Double_t measJt = GetPerpendicularMomentumTrackJet(jet, track);
        
        // AliAODMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
        // Always store the charge for the efficiency (needed for geant-fluka correction)
        Double_t value[AliAnalysisTaskPID::kEffNumAxes] = { static_cast<Double_t>(mcID), genPt, gentrack->Eta(), gentrack->Charge() / 3., centPercent, jetPt, genZ, genXi, genDistance, genJt };                                                     
        
        Double_t valueMeas[AliAnalysisTaskPID::kEffNumAxes] = { static_cast<Double_t>(mcID), trackPt, track->Eta(),
                                                                static_cast<Double_t>(track->Charge()),
                                                                centPercent, jetPt, measZ, measXi, measDistance, measJt }; 
        
        if (task) {
          task->FillEfficiencyContainer(value, AliAnalysisTaskPID::kStepRecWithGenCuts);
          task->FillEfficiencyContainer(valueMeas, AliAnalysisTaskPID::kStepRecWithGenCutsMeasuredObs);                                                   
        }
        else {
          for (Int_t i=0;i<fNumJetPIDtasks;++i) {
            if (!trackRejectedByTask[i])
              fJetPIDtask[i]->FillEfficiencyContainer(value, AliAnalysisTaskPID::kStepRecWithGenCuts);
              fJetPIDtask[i]->FillEfficiencyContainer(valueMeas, AliAnalysisTaskPID::kStepRecWithGenCutsMeasuredObs);                                                      
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
    mcID = AliAnalysisTaskPID::PDGtoMCID(pdg);
  }
  
  Double_t dca[2] = {0., 0.}; // 0: xy; 1: z
  
  Double_t v[3]   = {0, };
  Double_t pos[3] = {0, };
  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
  if (!primVtx) {
    std::cout << "Primary Vertex not found, do not fill DCA" << std::endl;
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