// $Id$
/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////////////////
//
//  Task to match cluster from a jet to GA patch.....
//
//
//         andrew.john.castro@cern.ch
//////////////////////////////////////////////////////////////////////////////////



#include "AliAnalysisTaskEmcalJetTriggerMatcher.h"

// general ROOT includes                                                                                                                                                  
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVector3.h>
#include <TObjArray.h>

// AliROOT includes                                                                                                                         
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliEmcalJet.h"
#include "AliAODJet.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include <AliVEvent.h>
#include <AliVParticle.h>
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalParticle.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALTriggerPatchInfo.h"

// event handler (and pico's) includes                                                                                                      
#include <AliInputEventHandler.h>
#include <AliVEventHandler.h>
#include "AliESDInputHandler.h"
#include "AliPicoTrack.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalJetTriggerMatcher)

//________________________________________________________________________
AliAnalysisTaskEmcalJetTriggerMatcher::AliAnalysisTaskEmcalJetTriggerMatcher() : 
  AliAnalysisTaskEmcalJet("heavyF",kFALSE),
  event(0),
  fEventTrigEMCALL1Gamma1(0),
  fEventTrigEMCALL1Gamma2(0),
  fInputEvent(0x0),
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0),
  fJetHIpt(20.0),
  fTrackEta(0.9),
  fFillHists(0),
  fMatchDist(0),
  fTracksCont(0), fCaloClustersCont(0),
  fESD(0), fAOD(0),
  fMaxPatch(0),
  fhnTriggerInfo(0),  //trigger QA
  fMainPatchType(kManual),
  //fMainPatchType(kEmcalJet),
  //fMainTrigCat(kTriggerLevel1Jet),
  fMainTrigCat(kTriggerLevel1Gamma),
  //fMainTrigCat(kTriggerRecalcJet),  // Recalculated max trigger patch; does not need to be above trigger threshold
  //fMainTrigCat(kTriggerRecalcGamma),
  fMainTrigSimple(kFALSE),
  fHistTriggerBitInfo(0),
  fHistMaxTriggerBitInfo(0),
  fHistEventSelection(0),
  fHistRecalcGASize(0),
  fHistRecalcGAEnergy(0),
  fHistClusEvPatchE(0),
  fHistdEtaPatchvdPhiPatch(0),
  fHistRawJetPtvPatchE(0),
  fHistMatchedClusJet(0),
  fHistMatchdetadphi(0),
  fHistdEtaPatchvdPhiPatchCMtoGeo(0),
  fHistEventClusSpect(0),
  fhnJetTrigger(0x0)
{
  // Default constructor.
  SetMakeGeneralHistograms(kTRUE);

}

//________________________________________________________________________
AliAnalysisTaskEmcalJetTriggerMatcher::AliAnalysisTaskEmcalJetTriggerMatcher(const char *name) :
  AliAnalysisTaskEmcalJet(name,kTRUE),
  event(0),
  fEventTrigEMCALL1Gamma1(0),
  fEventTrigEMCALL1Gamma2(0),
  fInputEvent(0x0),
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0),
  fJetHIpt(20.0),
  fTrackEta(0.9),
  fFillHists(0),
  fMatchDist(0),
  fTracksCont(0), fCaloClustersCont(0),
  fESD(0), fAOD(0),
  fMaxPatch(0),
  fhnTriggerInfo(0),  //trigger QA
  fMainPatchType(kManual),
  //fMainPatchType(kEmcalJet),
  //fMainTrigCat(kTriggerLevel1Jet),
  fMainTrigCat(kTriggerLevel1Gamma),
  //fMainTrigCat(kTriggerRecalcJet),  // Recalculated max trigger patch; does not need to be above trigger threshold
  //fMainTrigCat(kTriggerRecalcGamma),
  fMainTrigSimple(kFALSE),
  fHistTriggerBitInfo(0),
  fHistMaxTriggerBitInfo(0),
  fHistEventSelection(0),
  fHistRecalcGASize(0),
  fHistRecalcGAEnergy(0),
  fHistClusEvPatchE(0),
  fHistdEtaPatchvdPhiPatch(0),
  fHistRawJetPtvPatchE(0),
  fHistMatchedClusJet(0),
  fHistMatchdetadphi(0),
  fHistdEtaPatchvdPhiPatchCMtoGeo(0),
  fHistEventClusSpect(0),
  fhnJetTrigger(0x0)
{ 

   SetMakeGeneralHistograms(kTRUE);
 
   DefineInput(0,TChain::Class());
   DefineOutput(1, TList::Class());
}

//_______________________________________________________________________
AliAnalysisTaskEmcalJetTriggerMatcher::~AliAnalysisTaskEmcalJetTriggerMatcher()
{
  // destructor
  //
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetTriggerMatcher::UserCreateOutputObjects()
{
  if (! fCreateHisto)
    return;
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

 //no jets, just analysis tracks and clusters
  fTracksCont       = GetParticleContainer(0);
  fCaloClustersCont = GetClusterContainer(0);
  fTracksCont->SetClassName("AliVTrack");
  fCaloClustersCont->SetClassName("AliVCluster");
  TString histname;
  TString name;
  TString title;

  if(fFillHists>0){
  fHistClusEvPatchE                 = new TH2F("ClusterEvPatchE","Cluster Energy v Patch E; Cluster Energy (GeV); Trigger Patch Energy (GeV)",100,0,100,100,0,100);
  fHistdEtaPatchvdPhiPatch          = new TH2F("dEtaPatchdPhiPatch; #eta; #phi","dEtaPatchdPhiPatch; #Delta#eta; #Delta#phi",100,-0.1,0.1,100,-0.1,0.1);
  fHistRawJetPtvPatchE              = new TH2F("RawJetPtvTriggerE","RawJetPtvTriggerE; Jet Pt (GeV); Trigger Patch Energy (GeV)",100,0,100,100,0,100);
  fHistTriggerBitInfo               = new TH1F("TriggerBitInfo","TriggerBitInfo",15,0.5,15.5);
  fHistMaxTriggerBitInfo            = new TH1F("MaxPatchTriggerInfo","MaxPatchTriggerInfo",15,0.5,15.5);
  fHistEventSelection               = new TH1F("EventSelectionQA","EventSelectionQA",18,0.5,18.5);
  fHistRecalcGASize                 = new TH2F("RecalcGAPatchSize","Patch_Size; #eta(Towers); #phi(Towers)",40,0,40,40,0,40);
  fHistRecalcGAEnergy               = new TH1F("ReclacGAEnergy","RecalcGAEnergy; Energy(GeV)",150,0,150);
  fHistdEtaPatchvdPhiPatchCMtoGeo   = new TH2F("dEtadPhi_GApatchCMtoGeo","dEtadPhi_GApatchCMtoGeo; #eta; #phi",100,-0.1,0.1,100,-0.1,0.1);
  fHistEventClusSpect               = new TH1F ("EventClusterSpectrum","ClusterEnergySpectrum",105,-5,100);
      
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(1,"IsLevel0");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(2,"IsJetlow");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(3,"IsJetHigh");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(4,"IsGammalow");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(5,"IsGammaHigh");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(6,"IsMainTrigger");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(7,"IsJetLowSimple");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(8,"IsJetHighSimple");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(9,"IsGammaLowSimple");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(10,"IsGammaHighSimple");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(11,"IsMainTriggerSimple");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(12,"IsOfflineSimple");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(13,"IsRecalcJet");
  fHistTriggerBitInfo->GetXaxis()->SetBinLabel(14,"IsRecalcGamma");

  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(1,"IsLevel0");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(2,"IsJetlow");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(3,"IsJetHigh");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(4,"IsGammalow");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(5,"IsGammaHigh");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(6,"IsMainTrigger");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(7,"IsJetLowSimple");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(8,"IsJetHighSimple");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(9,"IsGammaLowSimple");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(10,"IsGammaHighSimple");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(11,"IsMainTriggerSimple");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(12,"IsOfflineSimple");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(13,"IsRecalcJet");
  fHistMaxTriggerBitInfo->GetXaxis()->SetBinLabel(14,"IsRecalcGamma");
      
  fHistEventSelection->GetXaxis()->SetBinLabel(1,"AliVEvent::kMB");
  fHistEventSelection->GetXaxis()->SetBinLabel(2,"AliVEvent::kINT7");
  fHistEventSelection->GetXaxis()->SetBinLabel(3,"AliVEvent::kHighMult");
  fHistEventSelection->GetXaxis()->SetBinLabel(4,"AliVEvent::kEMC1");
  fHistEventSelection->GetXaxis()->SetBinLabel(5,"AliVEvent::kCINT5");
  fHistEventSelection->GetXaxis()->SetBinLabel(6,"AliVEvent::kEMC7");
  fHistEventSelection->GetXaxis()->SetBinLabel(7,"AliVEvent::kEMC8");
  fHistEventSelection->GetXaxis()->SetBinLabel(8,"AliVEvent::kEMCEJE");
  fHistEventSelection->GetXaxis()->SetBinLabel(9,"AliVEvent::kEMCEGA");
  fHistEventSelection->GetXaxis()->SetBinLabel(10,"AliVEvent::kCentral");
  fHistEventSelection->GetXaxis()->SetBinLabel(11,"AliVEvent::kSemiCentral");
  fHistEventSelection->GetXaxis()->SetBinLabel(12,"AliVEvent::kZED");
  fHistEventSelection->GetXaxis()->SetBinLabel(13,"AliVEvent::kINT8");
  fHistEventSelection->GetXaxis()->SetBinLabel(14,"AliVEvent::kFastOnly");
  fHistEventSelection->GetXaxis()->SetBinLabel(15,"AliVEvent::kAnyINT");
  fHistEventSelection->GetXaxis()->SetBinLabel(16,"AliVEvent::kAny");
  fHistEventSelection->GetXaxis()->SetBinLabel(17,"Event_Counter");
    
    
  // PT bins used to be (2000, -100, 300)
      
  Int_t fgkNCentBins = 21;
  Float_t kMinCent   = 0.;
  Float_t kMaxCent   = 105.;
      
  Int_t fgkNVZEROBins = 100;
  Float_t kMinVZERO   = 0.;
  Float_t kMaxVZERO   = 25000;
      
  const Int_t fgkNEPatch = 100;
  Float_t kMinEPatch = 0.;
  Float_t kMaxEPatch = 200.;
      
  const Int_t fgkNADC = 100;
  Float_t kMinADC = 0.;
  Float_t kMaxADC = 1500.;
      
  const Int_t fgkNEta = 10;
  const Int_t fgkNPhi = 10;
      
  const Int_t nDim = 8;//cent;V0mult;ptjet1;ptjet2;Epatch;ADCpatch;EtaPatch;PhiPatch
  const Int_t nBins[nDim] = {fgkNCentBins,fgkNVZEROBins,fgkNEta,fgkNPhi,fgkNEPatch,fgkNADC,fgkNEta,fgkNPhi};
  const Double_t xmin0[nDim]  = {kMinCent,kMinVZERO,-0.7,1.4,kMinEPatch,kMinADC,-0.7,1.4};
  const Double_t xmax0[nDim]  = {kMaxCent,kMaxVZERO,0.7,3.14,kMaxEPatch,kMaxADC, 0.7,3.14};
  fhnTriggerInfo = new THnSparseF("fhnTriggerInfo", "hnTriggerInfo;cent;V0mult;EtaCMPatch;PhiCMPatch;Epatch;ADCpatch;EtaGeoPatch;PhiGeoPatch",nDim,nBins,xmin0,xmax0);
      
  fOutput->Add(fhnTriggerInfo);
  fOutput->Add(fHistClusEvPatchE);
  fOutput->Add(fHistdEtaPatchvdPhiPatch);
  fOutput->Add(fHistRawJetPtvPatchE);
  fOutput->Add(fHistdEtaPatchvdPhiPatchCMtoGeo);
  fOutput->Add(fHistEventClusSpect);
  fOutput->Add(fHistTriggerBitInfo);
  fOutput->Add(fHistMaxTriggerBitInfo);
  fOutput->Add(fHistEventSelection);
  fOutput->Add(fHistRecalcGASize);
  fOutput->Add(fHistRecalcGAEnergy);
      
  }//Fill Histograms
    
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(!inputHandler) {
    AliFatal("Input handler needed");
    return;
  }
  
    
  UInt_t bitcoded1 = 0;
  bitcoded1 = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8 | 1<<9 | 1<<10 | 1<<11 ;
  fhnJetTrigger = NewTHnSparseDJetTrigger("fhnJetTrigger", bitcoded1);
  
  cout << "_______________Created Sparse__________________" << endl;
  fOutput->Add(fhnJetTrigger);

  PostData(1, fOutput);

}

//________________________________________________________
void AliAnalysisTaskEmcalJetTriggerMatcher::ExecOnce()
{
  //  Initialize the analysis
  AliAnalysisTaskEmcalJet::ExecOnce();
  
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;
} // end of ExecOnce

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetTriggerMatcher::Run()
{
  // check to see if we have any tracks
  if (!fTracks)  return kTRUE;
  if (!fJets)  return kTRUE;
    
  UInt_t trig = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(fFillHists>0){
      if(trig & AliVEvent::kMB)          fHistEventSelection->Fill(1);
      if(trig & AliVEvent::kINT7)        fHistEventSelection->Fill(2);
      if(trig & AliVEvent::kHighMult)    fHistEventSelection->Fill(3);
      if(trig & AliVEvent::kEMC1)        fHistEventSelection->Fill(4);
      if(trig & AliVEvent::kCINT5)       fHistEventSelection->Fill(5);
      if(trig & AliVEvent::kEMC7)        fHistEventSelection->Fill(6);
      if(trig & AliVEvent::kEMC8)        fHistEventSelection->Fill(7);
      if(trig & AliVEvent::kEMCEJE)      fHistEventSelection->Fill(8);
      if(trig & AliVEvent::kEMCEGA)      fHistEventSelection->Fill(9);
      if(trig & AliVEvent::kCentral)     fHistEventSelection->Fill(10);
      if(trig & AliVEvent::kSemiCentral) fHistEventSelection->Fill(11);
      if(trig & AliVEvent::kZED)         fHistEventSelection->Fill(12);
      if(trig & AliVEvent::kINT8)        fHistEventSelection->Fill(13);
      if(trig & AliVEvent::kFastOnly)    fHistEventSelection->Fill(14);
      if(trig & AliVEvent::kAnyINT)      fHistEventSelection->Fill(15);
      if(trig & AliVEvent::kAny)         fHistEventSelection->Fill(16);
      fHistEventSelection->Fill(17);
    }
  // what kind of event do we have: AOD or ESD?
  Bool_t useAOD;
  if (dynamic_cast<AliAODEvent*>(InputEvent())) useAOD = kTRUE;
  else useAOD = kFALSE;
  
  fEventTrigEMCALL1Gamma1 = kFALSE;
  fEventTrigEMCALL1Gamma2 = kFALSE;
  
  // if we have ESD event, set up ESD object
  if(!useAOD){
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) {
      AliError(Form("ERROR: fESD not available\n"));
      return kTRUE;
    }
  }

  // if we have AOD event, set up AOD object
  if(useAOD){
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) {
      AliError(Form("ERROR: fAOD not available\n"));
      return kTRUE;
    }
  }
    
  // get vertex information
  Double_t fvertex[3]={0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(fvertex);

  // create pointer to list of input event                                                                                                  
  TList *list = InputEvent()->GetList();
  if(!list) {
    AliError(Form("ERROR: list not attached\n"));
    return kTRUE;
  }

  // background density                                                                                                                                                                                                                               
  fRhoVal = fRho->GetVal();

  // initialize TClonesArray pointers to jets and tracks                                                                                    
  TClonesArray *jets = 0;
  // get Jets object                                                                                                                        
  jets = dynamic_cast<TClonesArray*>(list->FindObject(fJets));
  if(!jets){
    AliError(Form("Pointer to jets %s == 0", fJets->GetName()));
    return kTRUE;
  } // verify existence of jets

  // get number of jets and tracks
  const Int_t Njets = jets->GetEntries();
  if(Njets<1)     return kTRUE;
  
    event++;
    //cout<<"Event #: "<<event<<"  Number of Clusters: "<<fCaloClustersCont->GetNClusters()<<"  Number of Tracks: "<<fTracksCont->GetNParticles()<<"  Number Of Jets: "<< Njets <<endl;
    
  if (fCaloClustersCont) {
    fCaloClustersCont->ResetCurrentID();
    AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster();
    while(cluster) {
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);
      if(fFillHists>0) fHistEventClusSpect->Fill(cluster->E());
      cluster = fCaloClustersCont->GetNextAcceptCluster();
        
    }
  }
  
    //  Start Jet Analysis
    // initialize jet parameters
    Int_t ijethi=-1;
    Double_t highestjetpt=0.0;
    Double_t GAeta = 0.0;
    Double_t GAphi = 0.0;
    Double_t GAetaCM = -999.;
    Double_t GAphiCM = -999.;
    Double_t PatchEE = 0.0;
    Double_t VZEROAmp = (Double_t)(InputEvent()->GetVZEROData()->GetTriggerChargeA() + InputEvent()->GetVZEROData()->GetTriggerChargeC());
    Double_t PatchEtaMax, PatchEtaMin, PatchPhiMax, PatchPhiMin, PatchAreaE, PatchAreaP;
    // Get GA Trigger Info
    if(fTriggerPatchInfo) {
      if(fMainPatchType==kManual) ExtractMainPatch();
      else if(fMainPatchType==kEmcalJet){
        fMaxPatch = GetMainTriggerPatch(fMainTrigCat,fMainTrigSimple);
      }
      PatchEE = fMaxPatch->GetPatchE();
      GAphi = fMaxPatch->GetPhiGeo();
      GAeta = fMaxPatch->GetEtaGeo();
      GAetaCM = fMaxPatch->GetEtaCM();
      GAphiCM = fMaxPatch->GetPhiCM();
      PatchEtaMax = fMaxPatch->GetEtaMax();
      PatchEtaMin = fMaxPatch->GetEtaMin();
      PatchPhiMax = fMaxPatch->GetPhiMax();
      PatchPhiMin = fMaxPatch->GetPhiMin();
      PatchAreaE = PatchEtaMax - PatchEtaMin;
      PatchAreaP = PatchPhiMax - PatchPhiMin;
       
      if(fFillHists>0){
        if (fMaxPatch->IsLevel0() == 1 )            fHistMaxTriggerBitInfo->Fill(1);
        if (fMaxPatch->IsJetLow() == 1 )            fHistMaxTriggerBitInfo->Fill(2);
        if (fMaxPatch->IsJetHigh() == 1 )           fHistMaxTriggerBitInfo->Fill(3);
        if (fMaxPatch->IsGammaLow() == 1 )          fHistMaxTriggerBitInfo->Fill(4);
        if (fMaxPatch->IsGammaHigh() == 1 )         fHistMaxTriggerBitInfo->Fill(5);
        if (fMaxPatch->IsMainTrigger() == 1 )       fHistMaxTriggerBitInfo->Fill(6);
        if (fMaxPatch->IsJetLowSimple() == 1 )      fHistMaxTriggerBitInfo->Fill(7);
        if (fMaxPatch->IsJetHighSimple() == 1 )     fHistMaxTriggerBitInfo->Fill(8);
        if (fMaxPatch->IsGammaLowSimple() == 1 )    fHistMaxTriggerBitInfo->Fill(9);
        if (fMaxPatch->IsGammaHighSimple() == 1 )   fHistMaxTriggerBitInfo->Fill(10);
        if (fMaxPatch->IsMainTriggerSimple() == 1 ) fHistMaxTriggerBitInfo->Fill(11);
        if (fMaxPatch->IsOfflineSimple() == 1 )     fHistMaxTriggerBitInfo->Fill(12);
        if (fMaxPatch->IsRecalcJet() == 1 )         fHistMaxTriggerBitInfo->Fill(13);
        
        if (fMaxPatch->IsRecalcGamma() == 1 ){
          fHistMaxTriggerBitInfo->Fill(14);
          fHistRecalcGASize->Fill(PatchAreaE/0.014, PatchAreaP/0.014);
          fHistRecalcGAEnergy->Fill(PatchEE);
        }
      }
        
      Double_t var[8] = {fCent, VZEROAmp, GAetaCM, GAphiCM, fMaxPatch->GetPatchE(), (Double_t)fMaxPatch->GetADCAmp(), fMaxPatch->GetEtaGeo(), fMaxPatch->GetPhiGeo()
        };
      if(fFillHists>0) fhnTriggerInfo->Fill(var);
      if(fFillHists>0) fHistdEtaPatchvdPhiPatchCMtoGeo->Fill(GAetaCM - GAeta, GAphiCM - GAphi);
        
    }
    // **********************************************************************
    //                JET LOOP
    // **********************************************************************
    // loop over jets in the event and make appropriate cuts
    for (Int_t iJets = 0; iJets < Njets; ++iJets) {
      AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(iJets));
      if (!jet) continue;
      
      // phi of jet, constrained to 1.6 < Phi < 2.94
      float jetphi = jet->Phi();      // phi of jet
      // apply jet cuts
      if(!AcceptMyJet(jet)) continue;
      Int_t JetClusters = jet->GetNumberOfClusters();
      Int_t JetTracks = jet -> GetNumberOfTracks();
      // Initializations and Calculations
      Double_t jetptraw = jet->Pt();    				             // raw pT of jet
      Double_t jetPt = -500;                                     // initialize corr jet pt LOCAL
      Double_t jetarea = -500;					               // initialize jet area
      jetarea = jet->Area();		           		              // jet area
      jetPt = jet->Pt() - jetarea*fRhoVal;                    // semi-corrected pT of jet from GLOBAL rho value
        
      if(jet->Pt() > fJetHIpt) {


        Float_t jeteta = -999.;
        Float_t jetphi = -999.;
        jeteta = jet->Eta();
        jetphi = jet->Phi();
          
        Double_t deta = 999;
        Double_t dphi = 999;
          
     //******************************Cluster Matched
     //**************************************************************
     for(int iCluster = 0; iCluster <= fCaloClustersCont->GetNClusters(); iCluster++){
        //Get closest track to cluster to track matching!!!!!
        //AliVCluster *cluster = fCaloClustersCont->GetNextAcceptedCluster(iCluster);
        AliVCluster *clusMatch = fCaloClustersCont->GetCluster(iCluster);
        if(!clusMatch) continue;
        if(! IsJetCluster(jet, iCluster, kFALSE)) continue;
       
        Double_t mClusterE = clusMatch->E();
        Float_t pos_mc[3];
        clusMatch->GetPosition(pos_mc);  // Get cluster position
        TVector3 mcp(pos_mc);
        Double_t maxMatchedPatchEta = (fMatchDist * 0.014);
        Double_t maxMatchedPatchPhi = (fMatchDist * 0.014);
        Double_t mtchPosition = maxMatchedPatchPhi * maxMatchedPatchPhi + maxMatchedPatchEta * maxMatchedPatchEta;
        Double_t rGApatch = TMath::Sqrt(mtchPosition);
        Double_t rJetCluster = -999.;
            
        Double_t PatchClusterdiff = (mcp.Phi() - fMaxPatch->GetPhiGeo())*(mcp.Phi() - fMaxPatch->GetPhiGeo()) + (mcp.PseudoRapidity() - fMaxPatch->GetEtaGeo())*(mcp.PseudoRapidity() - fMaxPatch->GetEtaGeo());
        rJetCluster = TMath::Sqrt(PatchClusterdiff);

          //Matched Trigger Jet
        if ( rJetCluster <= rGApatch ){

          Double_t dEtaPatch = 0.0;
          Double_t dPhiPatch = 0.0;
          dEtaPatch = fMaxPatch->GetEtaGeo() - mcp.PseudoRapidity();
          dPhiPatch = fMaxPatch->GetPhiGeo() - mcp.Phi();
          if(fFillHists>0) fHistClusEvPatchE->Fill(mClusterE,fMaxPatch->GetPatchE());
          if(fFillHists>0) fHistdEtaPatchvdPhiPatch->Fill(dEtaPatch,dPhiPatch);
          if(fFillHists>0) fHistRawJetPtvPatchE->Fill(jetptraw,fMaxPatch->GetPatchE());
          Double_t JetTrig[18] = {jetptraw, jetPt, jeteta, jetphi, jet->E(), fMaxPatch->GetPatchE(), fMaxPatch->GetEtaGeo(), fMaxPatch->GetPhiGeo(), (Double_t)fMaxPatch->GetADCAmp(), mClusterE, mcp.PseudoRapidity(), mcp.Phi() };
          fhnJetTrigger->Fill(JetTrig);
        } // matching cluster of a jet to trigger patch
      } // cluster for
    } // highest pt jet cut
  } // LOOP over JETS in event

  return kTRUE;
  
}
//________________________________________________________________________
void AliAnalysisTaskEmcalJetTriggerMatcher::Terminate(Option_t *)
{
  cout<<"###########################"<<endl;
  cout<<"####   Task Finished   ####"<<endl;
  cout<<"###########################"<<endl;
  cout<<"###########################"<<endl;
} // end of terminate


//________________________________________________________________________
void AliAnalysisTaskEmcalJetTriggerMatcher::ExtractMainPatch() {
    
    //Find main trigger
    if(!fTriggerPatchInfo)
        return;
    
    //number of patches in event
    Int_t nPatch = fTriggerPatchInfo->GetEntriesFast();
    //extract main trigger patch
    Double_t emax = -1.;
    for (Int_t iPatch = 0; iPatch < nPatch; iPatch++) {
        AliEMCALTriggerPatchInfo *patch = (AliEMCALTriggerPatchInfo*)fTriggerPatchInfo->At( iPatch );
        if (!patch) continue;
        
        if(fFillHists>0){
          if (patch->IsLevel0() == 1 )            fHistTriggerBitInfo->Fill(1);
          if (patch->IsJetLow() == 1 )            fHistTriggerBitInfo->Fill(2);
          if (patch->IsJetHigh() == 1 )           fHistTriggerBitInfo->Fill(3);
          if (patch->IsGammaLow() == 1 )          fHistTriggerBitInfo->Fill(4);
          if (patch->IsGammaHigh() == 1 )         fHistTriggerBitInfo->Fill(5);
          if (patch->IsMainTrigger() == 1 )       fHistTriggerBitInfo->Fill(6);
          if (patch->IsJetLowSimple() == 1 )      fHistTriggerBitInfo->Fill(7);
          if (patch->IsJetHighSimple() == 1 )     fHistTriggerBitInfo->Fill(8);
          if (patch->IsGammaLowSimple() == 1 )    fHistTriggerBitInfo->Fill(9);
          if (patch->IsGammaHighSimple() == 1 )   fHistTriggerBitInfo->Fill(10);
          if (patch->IsMainTriggerSimple() == 1 ) fHistTriggerBitInfo->Fill(11);
          if (patch->IsOfflineSimple() == 1 )     fHistTriggerBitInfo->Fill(12);
          if (patch->IsRecalcJet() == 1 )         fHistTriggerBitInfo->Fill(13);
          if (patch->IsRecalcGamma() == 1 )       fHistTriggerBitInfo->Fill(14);
        }
        //Force ExtractMainPatch to return Recalc GA info
        if (patch->IsRecalcGamma() == 1 ){
            if(patch->GetPatchE()>emax) {
                fMaxPatch = patch;
                emax = patch->GetPatchE();
            }
        }
    }
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetTriggerMatcher::AcceptMyJet(AliEmcalJet *jet) {
  //applies all jet cuts except pt
  if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax)) return 0;
  if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax)) return 0;
  if (jet->Area()<fAreacut) return 0;
  // prevents 0 area jets from sneaking by when area cut == 0
  if (jet->Area()==0) return 0;
  //exclude jets with extremely high pt tracks which are likely misreconstructed
  if(jet->MaxTrackPt()>100) return 0;
  
  //passed all above cuts
  return 1;
}


THnSparse* AliAnalysisTaskEmcalJetTriggerMatcher::NewTHnSparseDJetTrigger(const char* name, UInt_t entries)
{
    // generate new THnSparseD JetQA, axes are defined in GetDimParamsJetQA()
    Int_t count = 0;
    UInt_t tmp = entries;
    while(tmp!=0){
        count++;
        tmp = tmp &~ -tmp;  // clear lowest bit
    }
    
    TString hnTitle(name);
    const Int_t dim = count;
    Int_t nbins[dim];
    Double_t xmin[dim];
    Double_t xmax[dim];
    
    Int_t i=0;
    Int_t c=0;
    while(c<dim && i<32){
        if(entries&(1<<i)){
            
            TString label("");
            GetDimParamsJetTrigger(i, label, nbins[c], xmin[c], xmax[c]);
            hnTitle += Form(";%s",label.Data());
            c++;
        }
        
        i++;
    }
    hnTitle += ";";
    
    return new THnSparseD(name, hnTitle.Data(), dim, nbins, xmin, xmax);
} // end of NewTHnSparseF JetTrigger



void AliAnalysisTaskEmcalJetTriggerMatcher::GetDimParamsJetTrigger(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
    // stores label and binning of axis for THnSparse
    const Double_t pi = TMath::Pi();
    
    switch(iEntry){
            
        case 0:
            label = "Jet Pt";
            nbins = 50;
            xmin = 0.;
            xmax = 250.;
            break;
            
        case 1:
            label = "Jet Corr Pt";
            nbins = 50;
            xmin = -50.;
            xmax = 200.;
            break;
            
        case 2:
            label = "Jet Eta";
            nbins = 24;
            xmin = -1.2;
            xmax =  1.2;
            break;
            
        case 3:
            label = "Jet Phi";
            nbins = 72;
            xmin = 0;
            xmax = 2*pi;
            break;
            
        case 4:
            label = "Jet E";
            nbins = 50;
            xmin = 0.;
            xmax = 250.;
            break;
            
        case 5:
            label = "Trigger Patch E";
            nbins = 100;
            xmin = 0.;
            xmax = 100.;
            break;
            
        case 6:
            label = "Trigger Patch Eta";
            nbins = 24;
            xmin = -1.2;
            xmax =  1.2;
            break;
            
        case 7:
            label = "Trigger Patch Phi";
            nbins = 72;
            xmin = 0;
            xmax = 2*pi;
            break;
            
        case 8:
            label = "Trigger Patch ADC";
            nbins = 150;
            xmin = 0;
            xmax = 1500;
            break;
            
        case 9:
            label = "Cluster E";
            nbins = 100;
            xmin = 0.;
            xmax = 100.;
            break;
            
        case 10:
            label = "Cluster Eta";
            nbins = 24;
            xmin = -1.2;
            xmax =  1.2;
            break;
            
        case 11:
            label = "Cluster Phi";
            nbins = 72;
            xmin = 0;
            xmax = 2*pi;
            break;
            
    } // end of switch
} // end of getting dim-params
