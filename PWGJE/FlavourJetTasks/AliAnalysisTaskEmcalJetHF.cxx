// $Id$
/****************************************************************************
*   Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
*                                                                          *
*   Author: The ALICE Off-line Project.                                    *
*   Contributors are mentioned in the code where appropriate.              *
*                                                                          *
*   Permission to use, copy, modify and distribute this software and its   *
*   documentation strictly for non-commercial purposes is hereby granted   *
*   without fee, provided that the above copyright notice appears in all   *
*   copies and that both the copyright notice and this permission notice   *
*   appear in the supporting documentation. The authors make no claims     *
*   about the suitability of this software for any purpose. It is          *
*   provided "as is" without express or implied warranty.                  *
****************************************************************************/
///////////////////////////////////////////////////////////////////////////////////////////
//  Task to tag HFE Jets using the EMCal
//          Andrew Castro (UTK)
//      andrew.john.castro@cern.ch
///////////////////////////////////////////////////////////////////////////////////////////
#include "AliAnalysisTaskEmcalJetHF.h"
// general ROOT includes                                                                                                                                                  
#include <TCanvas.h>
#include <TChain.h>
#include <TMath.h>
#include <TProfile.h>
#include <TAxis.h>
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
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalParticle.h"
#include "AliESDCaloCluster.h"
#include <AliESDtrackCuts.h>
#include "AliPID.h"
#include "AliTPCdEdxInfo.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliVVZERO.h"
#include "AliAnalysisTaskEmcal.h"
#include <AliInputEventHandler.h>
#include <AliVEventHandler.h>
#include "AliEmcalTriggerSetupInfo.h"
#include "AliAODHandler.h"
#include "AliESDInputHandler.h"
#include "AliPicoTrack.h"
#include "AliEventPoolManager.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
// PID includes                                                                                                                             
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliESDpid.h"
// magnetic field includes
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalJetHF)

//________________________________________________________________________
AliAnalysisTaskEmcalJetHF::AliAnalysisTaskEmcalJetHF() : 
  AliAnalysisTaskEmcalJet("heavyF",kFALSE), 
  event(0),
  fFillHists(0),
  fEventTrigEMCALL1Gamma1(0),
  fEventTrigEMCALL1Gamma2(0),
  fGlobalQA(0),
  fJetPID(0),
  fInputEvent(0x0),
  fCuts(0),
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0),
  fJetHIpt(20.0),
  fTrackPtCut(2.0),
  fTrackEta(0.9),
  fTrkQAcut(0),
  fM02max(0.35),
  fM02min(0.006),
  fesdTrackCuts(0),
  fPIDResponse(0x0), fTPCResponse(),
  fEsdtrackCutsITSTPC(),
  fEsdtrackCutsTPC(),
  fEsdtrackCutsITS(),
  fJetsCont(0), fTracksCont(0), fCaloClustersCont(0), fTracksJetCont(0), fCaloClustersJetCont(0),
  fESD(0), fAOD(0),
  fMaxPatch(0),
  fhnTriggerInfo(0),
  fMainPatchType(kManual),
  //fMainPatchType(kEmcalJet),
  fMainTrigCat(kTriggerLevel0),
  //fMainTrigCat(kTriggerLevel1Gamma),
  //fMainTrigCat(kTriggerRecalcJet),  // Recalculated max trigger patch; does not need to be above trigger threshold
  //fMainTrigCat(kTriggerRecalcGamma),
  fMainTrigSimple(kTRUE),
  fHistTriggerBitInfo(0),
  fHistMaxTriggerBitInfo(0),
  fHistEventSelection(0),
  fHistRecalcGASize(0),
  fHistRecalcGAEnergy(0),
  fHistCorrJetEvsPatchE(0),
  fHistClusEvPatchE(0),
  fHistdEtaPatchvdPhiPatch(0),
  fHistRawJetEvPatchE(0),
  fHistMatchedClusJet(0),
  fhnTrackClusterQA(0x0), fhnPIDHFTtoC(0x0), fhnJetTrigger(0x0)
{
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetHF::AliAnalysisTaskEmcalJetHF(const char *name) :
  AliAnalysisTaskEmcalJet(name,kTRUE),
  event(0),
  fFillHists(0),
  fEventTrigEMCALL1Gamma1(0),
  fEventTrigEMCALL1Gamma2(0),
  fGlobalQA(0),
  fJetPID(0),
  fInputEvent(0x0),
  fCuts(0),
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0),
  fJetHIpt(20.0),
  fTrackPtCut(2.0),
  fTrackEta(0.9),
  fTrkQAcut(0),
  fM02max(0.35),
  fM02min(0.006),
  fesdTrackCuts(0),
  fPIDResponse(0x0), fTPCResponse(),
  fEsdtrackCutsITSTPC(),
  fEsdtrackCutsTPC(),
  fEsdtrackCutsITS(),
  fJetsCont(0), fTracksCont(0), fCaloClustersCont(0), fTracksJetCont(0), fCaloClustersJetCont(0),
  fESD(0), fAOD(0),
  fMaxPatch(0),
  fhnTriggerInfo(0),
  fMainPatchType(kManual),
  //fMainPatchType(kEmcalJet),
  fMainTrigCat(kTriggerLevel0),
  //fMainTrigCat(kTriggerLevel1Gamma),
  //fMainTrigCat(kTriggerRecalcJet),  // Recalculated max trigger patch; does not need to be above trigger threshold
  //fMainTrigCat(kTriggerRecalcGamma),
  fMainTrigSimple(kTRUE),
  fHistTriggerBitInfo(0),
  fHistMaxTriggerBitInfo(0),
  fHistEventSelection(0),
  fHistRecalcGASize(0),
  fHistRecalcGAEnergy(0),
  fHistCorrJetEvsPatchE(0),
  fHistClusEvPatchE(0),
  fHistdEtaPatchvdPhiPatch(0),
  fHistRawJetEvPatchE(0),
  fHistMatchedClusJet(0),
  fhnTrackClusterQA(0x0), fhnPIDHFTtoC(0x0), fhnJetTrigger(0x0)
{ 
   SetMakeGeneralHistograms(kTRUE);
   DefineInput(0,TChain::Class());
   DefineOutput(1, TList::Class());
}
//_______________________________________________________________________
AliAnalysisTaskEmcalJetHF::~AliAnalysisTaskEmcalJetHF()
{
  // destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}
//________________________________________________________________________
void AliAnalysisTaskEmcalJetHF::UserCreateOutputObjects()
{
  if (! fCreateHisto) return;
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  //fJetsCont           = GetJetContainer(0);
  //if(fJetsCont) { //get particles and clusters connected to jets
  //  fTracksJetCont       = fJetsCont->GetParticleContainer();
  //  fCaloClustersJetCont = fJetsCont->GetClusterContainer();
  //}
 //else {        //no jets, just analysis tracks and clusters
  fTracksCont       = GetParticleContainer(0);
  fCaloClustersCont = GetClusterContainer(0);
//}
  fTracksCont->SetClassName("AliVTrack");
  fCaloClustersCont->SetClassName("AliVCluster");
  TString histname;
  TString name;
  TString title;
    
  if(fFillHists>0){
  fHistMatchedClusJet          = new TH1F("MatchedGAClusterEE","ClusterEE; Energy (GeV)",100,0,100);
  fHistCorrJetEvsPatchE        = new TH2F("Corr Jet E vs Trigger Patch E","CorrJetEvTriggerE; Corrected Jet Energy (GeV); Patch Energy (GeV)",100,0,100,100,0,100);
  fHistClusEvPatchE            = new TH2F("ClusterEvPatchE","Cluster Energy v Patch E; Cluster Energy (GeV); Patch Energy (GeV)",100,0,100,100,0,100);
  fHistdEtaPatchvdPhiPatch     = new TH2F("dEtaPatchdPhiPatch; #eta; #phi","dEtaPatchdPhiPatch; #Delta#eta; #Delta#phi ",100,-0.1,0.1,100,-0.1,0.1);
  fHistRawJetEvPatchE          = new TH2F("RaeJetEvTriggerE","RawJetEvTriggerE; Jet Energy (GeV); Patch Energy (GeV)",100,0,100,100,0,100);
  fHistTriggerBitInfo          = new TH1F("TriggerBitInfo","TriggerBitInfo; ; Counts",15,0.5,15.5);
  fHistMaxTriggerBitInfo       = new TH1F("MaxPatchTriggerInfo","MaxPatchTriggerInfo; ;Counts",15,0.5,15.5);
  fHistEventSelection          = new TH1F("EventSelectionQA","EventSelectionQA",17,0.5,17.5);
  fHistRecalcGASize            = new TH2F("RecalcGAPatchSize","Patch_Size; #eta(Towers); #phi(Towers)",40,0,40,40,0,40);
  fHistRecalcGAEnergy          = new TH1F("ReclacGAEnergy","RecalcGAEnergy; Energy(GeV); Counts",150,0,150);
      
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
    
  const Int_t fgkNCentBins = 21;
  Float_t kMinCent   = 0.;
  Float_t kMaxCent   = 105.;
  const Int_t fgkNVZEROBins = 100;
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
  //Trigger QA sparse in here for debugging purposes
  fhnTriggerInfo = new THnSparseF("fhnTriggerInfo", "hnTriggerInfo;cent;V0mult;EtaPatchCM;PhiPatchCM;Epatch;ADCpatch;EtaPatchGeo;PhiPatchGEo",nDim,nBins,xmin0,xmax0);
      
  fOutput->Add(fhnTriggerInfo);
  fOutput->Add(fHistMatchedClusJet);
  fOutput->Add(fHistCorrJetEvsPatchE);
  fOutput->Add(fHistClusEvPatchE);
  fOutput->Add(fHistdEtaPatchvdPhiPatch);
  fOutput->Add(fHistRawJetEvPatchE);
  fOutput->Add(fHistTriggerBitInfo);
  fOutput->Add(fHistMaxTriggerBitInfo);
  fOutput->Add(fHistEventSelection);
  fOutput->Add(fHistRecalcGASize);
  fOutput->Add(fHistRecalcGAEnergy);
      
  }//Fill Histograms

  // ****************************** PID *****************************************************                                               
  // set up PID handler                                                                                                                     
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(!inputHandler) {
    AliFatal("Input handler needed");
    return;
  }

  // PID response object                                                                                                                    
  //fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();                                                                         
  //  inputHandler->CreatePIDResponse(fIsMC);         // needed to create object, why though?                                                 
  fPIDResponse = inputHandler->GetPIDResponse();
  if (!fPIDResponse) {
    AliError("PIDResponse object was not created");
    return;
  }
  // ****************************************************************************************
  UInt_t bitcoded3 = 0;  // bit coded, see GetDimParamsPID() below
  bitcoded3 = 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<8 | 1<<9 | 1<<10 | 1<<15 | 1<<16 | 1<<17;
  fhnTrackClusterQA = NewTHnSparseDHF("fhnTrackClusterQA", bitcoded3);
  
  UInt_t bitcoded7 = 0;  // bit coded, see GetDimParamsPID() below
  bitcoded7 = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8 | 1<<9 | 1<<10 | 1<<11 | 1<<12 | 1<<13| 1<<14 | 1<<15 | 1<<16 | 1<<17;
  fhnPIDHFTtoC = NewTHnSparseDHF("fhnPIDHFTtoC", bitcoded7);
    
  UInt_t bitcoded8 = 0;
  bitcoded8 = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8 | 1<<9 | 1<<10 | 1<<11 | 1<<12 | 1<<13 | 1<<14 | 1<<15 | 1<<16 | 1<<17;
  fhnJetTrigger = NewTHnSparseDJetTrigger("fhnJetTrigger", bitcoded8);

  if(fJetPID>0) fOutput->Add(fhnPIDHFTtoC);
  if(fJetPID>0) fOutput->Add(fhnJetTrigger);
  if(fGlobalQA>0)  fOutput->Add(fhnTrackClusterQA);
  PostData(1, fOutput);
}
//________________________________________________________
void AliAnalysisTaskEmcalJetHF::ExecOnce()
{
  //  Initialize the analysis
  AliAnalysisTaskEmcalJet::ExecOnce();
  
  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;


} // end of ExecOnce

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetHF::Run()
{
  // check to see if we have any tracks
  if (!fTracks)  return kTRUE;
  if (!fJets)    return kTRUE;

  UInt_t trig = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(fFillHists>0)     if(trig & AliVEvent::kMB)          fHistEventSelection->Fill(1);
  if(fFillHists>0)     if(trig & AliVEvent::kINT7)        fHistEventSelection->Fill(2);
  if(fFillHists>0)     if(trig & AliVEvent::kHighMult)    fHistEventSelection->Fill(3);
  if(fFillHists>0)     if(trig & AliVEvent::kEMC1)        fHistEventSelection->Fill(4);
  if(fFillHists>0)     if(trig & AliVEvent::kCINT5)       fHistEventSelection->Fill(5);
  if(fFillHists>0)     if(trig & AliVEvent::kEMC7)        fHistEventSelection->Fill(6);
  if(fFillHists>0)     if(trig & AliVEvent::kEMC8)        fHistEventSelection->Fill(7);
  if(fFillHists>0)     if(trig & AliVEvent::kEMCEJE)      fHistEventSelection->Fill(8);
  if(fFillHists>0)     if(trig & AliVEvent::kEMCEGA)      fHistEventSelection->Fill(9);
  if(fFillHists>0)     if(trig & AliVEvent::kCentral)     fHistEventSelection->Fill(10);
  if(fFillHists>0)     if(trig & AliVEvent::kSemiCentral) fHistEventSelection->Fill(11);
  if(fFillHists>0)     if(trig & AliVEvent::kZED)         fHistEventSelection->Fill(12);
  if(fFillHists>0)     if(trig & AliVEvent::kINT8)        fHistEventSelection->Fill(13);
  if(fFillHists>0)     if(trig & AliVEvent::kFastOnly)    fHistEventSelection->Fill(14);
  if(fFillHists>0)     if(trig & AliVEvent::kAnyINT)      fHistEventSelection->Fill(15);
  if(fFillHists>0)     if(trig & AliVEvent::kAny)         fHistEventSelection->Fill(16);
    
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
    
  Double_t MagSign;
  Double_t MagF;
    
  if(!useAOD){
    // get magnetic field info for DCA
    MagF = fESD->GetMagneticField();
    MagSign = 1.0;
    if(MagF<0)MagSign = -1.0;
    // set magnetic field
    if (!TGeoGlobalMagField::Instance()->GetField()) {
        AliMagF* field = new AliMagF("Maps","Maps", MagSign, MagSign, AliMagF::k5kG);
        TGeoGlobalMagField::Instance()->SetField(field);
    }
  }
    
  if(useAOD){
    // get magnetic field info for DCA
    MagF = fAOD->GetMagneticField();
    MagSign = 1.0;
    if(MagF<0)MagSign = -1.0;
    // set magnetic field
    if (!TGeoGlobalMagField::Instance()->GetField()) {
        AliMagF* field = new AliMagF("Maps","Maps", MagSign, MagSign, AliMagF::k5kG);
        TGeoGlobalMagField::Instance()->SetField(field);
    }
  }

  // get centrality bin
  Int_t centbin = GetCentBin(fCent);
  //for pp analyses we will just use the first centrality bin
  if (centbin == -1)  centbin = 0;

  // get vertex information
  Double_t fvertex[3]={0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(fvertex);
  //Double_t zVtx=fvertex[2];

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
  //TClonesArray *tracks = 0;
  //TClonesArray *clusters = 0;
  //TClonesArray * clusterList = 0;
  
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
  if (fTracksCont) {
    fTracksCont->ResetCurrentID();
    AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    while(track) {
    track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    }
  }
  if (fCaloClustersCont) {
    fCaloClustersCont->ResetCurrentID();
    AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster();
    while(cluster) {
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);
      cluster = fCaloClustersCont->GetNextAcceptCluster();
    }
  }
  //  Start Jet Analysis
  // initialize jet parameters
  Double_t GAeta = 0.0;
  Double_t GAphi = 0.0;
  Double_t PatchEE = 0.0;
  Double_t VZEROAmp = (Double_t)(InputEvent()->GetVZEROData()->GetTriggerChargeA() + InputEvent()->GetVZEROData()->GetTriggerChargeC());
  Double_t PatchEtaMax, PatchEtaMin, PatchPhiMax, PatchPhiMin, PatchAreaE, PatchAreaP;
  // Get GA Trigger Info
  if(fTriggerPatchInfo) {
    if(fMainPatchType==kManual) ExtractMainPatch();
    else if(fMainPatchType==kEmcalJet) fMaxPatch = GetMainTriggerPatch(fMainTrigCat,fMainTrigSimple);

    PatchEE = fMaxPatch->GetPatchE();
    GAphi = fMaxPatch->GetPhiGeo();
    GAeta = fMaxPatch->GetEtaGeo();
    PatchEtaMax = fMaxPatch->GetEtaMax();
    PatchEtaMin = fMaxPatch->GetEtaMin();
    PatchPhiMax = fMaxPatch->GetPhiMax();
    PatchPhiMin = fMaxPatch->GetPhiMin();
    PatchAreaE = PatchEtaMax - PatchEtaMin;
    PatchAreaP = PatchPhiMax - PatchPhiMin;
    Double_t var[8] = {fCent, VZEROAmp, fMaxPatch->GetEtaCM(), fMaxPatch->GetPhiCM(), fMaxPatch->GetPatchE(), (Double_t)fMaxPatch->GetADCAmp(), fMaxPatch->GetEtaGeo(), fMaxPatch->GetPhiGeo() };
    
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
      if (fMaxPatch->IsRecalcGamma() == 1 )       fHistMaxTriggerBitInfo->Fill(14);
      fHistRecalcGASize->Fill(PatchAreaE/0.014, PatchAreaP/0.014);
      fhnTriggerInfo->Fill(var);
      fHistRecalcGAEnergy->Fill(PatchEE);
    }//Fill Patch Histograms and Trigger QA sparse
  }//Max Energy Patch info
    
  // **********************************************************************
  //                JET LOOP
  // **********************************************************************
  // loop over jets in the event and make appropriate cuts
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(iJets));
    if (!jet) continue;
      
    // phi of jet, constrained to 1.6 < Phi < 2.94
    Double_t jetphi = jet->Phi();                              // phi of jet
    // apply jet cuts
    if(!AcceptMyJet(jet)) continue;
      
    Int_t JetClusters = jet->GetNumberOfClusters();
    Int_t JetTracks = jet -> GetNumberOfTracks();
    Int_t tag = -999;
    // Initializations and Calculations
    Double_t jetptraw = jet->Pt();    				           // raw pT of jet
    Double_t jetPt = -500;                                     // initialize corr jet pt LOCAL
    Double_t jetarea = -500;					               // initialize jet area
    jetarea = jet->Area();		           		               // jet area
    jetPt = jet->Pt() - jetarea*fRhoVal;                       // semi-corrected pT of jet from GLOBAL rho value
     
    if(jet->Pt() > fJetHIpt) {
      if(!fTracksCont || !fCaloClustersCont) continue;
      Double_t dEdx = -99;
      Double_t EovP = -99;
      Double_t DCAxy = -999;
      Double_t DCAz = -999;
      Double_t jeteta = -999.;
      jeteta = jet->Eta();
      jetphi = jet->Phi();
      //******************************Cluster Matched To Closest Track
      //**************************************************************
      Int_t NumbTrackContainer = -999;
      NumbTrackContainer = fTracksCont->GetNParticles();
      for(int iTracks = 0; iTracks <= NumbTrackContainer; iTracks++){
        AliVTrack *AcceptedTrack =static_cast<AliVTrack*>(fTracksCont->GetParticle(iTracks));
        if(!AcceptedTrack){
          AliError(Form("Couldn't get AliVTrack Container %d\n", iTracks));
          continue;
        }
        if(!IsJetTrack(jet,iTracks,kFALSE))continue;
        //Get matched cluster
        Int_t emc1 = -999;
        emc1 = AcceptedTrack->GetEMCALcluster();//Get EMCal Cluster Matched Index
        if(emc1 < 0) continue;
        Int_t TPCNclus = -999, ITSNclus = -999;
        Double_t acceptTrackP = AcceptedTrack->P();
        Double_t acceptTrackPt = AcceptedTrack->Pt();
        Double_t acceptTrackEta = AcceptedTrack->Eta();
        Double_t acceptTrackPhi = AcceptedTrack->Phi();
        Double_t nSigmaElectron_TPC_at = fPIDResponse->NumberOfSigmasTPC(AcceptedTrack,AliPID::kElectron);
        Double_t nSigmaElectron_TOF_at = fPIDResponse->NumberOfSigmasTOF(AcceptedTrack,AliPID::kElectron);
        Double_t dEdxat = AcceptedTrack->GetTPCsignal();

        if(!useAOD){
          AliESDtrack *ESDacceptedTrack = static_cast<AliESDtrack*>(AcceptedTrack);
          if(!ESDacceptedTrack){
            AliError(Form("Couldn't get AliESDTrack %d\n", iTracks));
            continue;
          }
        //ESDacceptedTrack->GetImpactParameters(DCAxy_at, DCAz_at);
        TPCNclus = ESDacceptedTrack->GetTPCNcls();
        ITSNclus = ESDacceptedTrack->GetITSNcls();
        }
        if(useAOD){
          AliAODTrack *AODacceptedTrack = static_cast<AliAODTrack*>(AcceptedTrack);
          if(!AODacceptedTrack) continue;
          if(!fPIDResponse) continue;
          TPCNclus = AODacceptedTrack->GetTPCNcls();
          ITSNclus = AODacceptedTrack->GetITSNcls();
          // AODacceptedTrack->GetImpactParameters(DCAxy_at, DCAz_at);
        }
        if(fCaloClustersCont && emc1>=0) {
          AliVCluster *clusMatch = fCaloClustersCont->GetCluster(emc1);
          if(!clusMatch){
            AliError(Form("Couldn't get matched AliVCluster %d\n", emc1));
            continue;
          }
          Double_t m02 = -999.;
          m02 = clusMatch->GetM02();
          if(TPCNclus < 80) continue;
          if(ITSNclus < 3) continue;
          if(m02 <= fM02max && m02 >= fM02min) continue;
          Double_t mClusterE = clusMatch->E();
          Float_t pos_mc[3];
          clusMatch->GetPosition(pos_mc);  // Get cluster position
          TVector3 mcp(pos_mc);
          Double_t EovP_mc = -999;
          EovP_mc = mClusterE/acceptTrackP;            
        if(fJetPID>0){ 
          Double_t HF_tracks2[18] = {fCent, acceptTrackPt, acceptTrackP ,acceptTrackEta, acceptTrackPhi, EovP_mc, 0, 0, dEdxat,nSigmaElectron_TPC_at, nSigmaElectron_TOF_at,0 , jetPt, jet->Phi(), jet->Eta(),mClusterE,mcp.PseudoRapidity(),mcp.Phi()};
          fhnPIDHFTtoC->Fill(HF_tracks2);    // fill Sparse Histo with trigger entries
        }
        //Double_t minMatchedPatchEta = mcp.PseudoRapidity() - (6.0 * 0.014);
        Double_t maxMatchedPatchEta = (7.0 * 0.014);
        //Double_t minMatchedPatchPhi = mcp.Phi() - (6.0 * 0.014);
        Double_t maxMatchedPatchPhi = (7.0 * 0.014);
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
          Double_t EovPJetGA = mClusterE / acceptTrackP;
          Double_t nSigmaElectron_TPC_GA = fPIDResponse->NumberOfSigmasTPC(AcceptedTrack,AliPID::kElectron);
                
          if(fFillHists>0){
            fHistCorrJetEvsPatchE->Fill(jetPt,fMaxPatch->GetPatchE());
            fHistClusEvPatchE->Fill(mClusterE,fMaxPatch->GetPatchE());
            fHistdEtaPatchvdPhiPatch->Fill(dEtaPatch,dPhiPatch);
            fHistRawJetEvPatchE->Fill(jetptraw,fMaxPatch->GetPatchE());
          }//Fill Histos
          Double_t JetTrig[18] = {jetptraw, jetPt, jeteta, jetphi, jet->E(), fMaxPatch->GetPatchE(), fMaxPatch->GetEtaGeo(), fMaxPatch->GetPhiGeo(), (Double_t)fMaxPatch->GetADCAmp(), mClusterE, mcp.PseudoRapidity(), mcp.Phi(),acceptTrackP, acceptTrackPt, acceptTrackEta, acceptTrackPhi, EovPJetGA, nSigmaElectron_TPC_GA};
          if(fJetPID>0) fhnJetTrigger->Fill(JetTrig);
        }//Trigger Patch Matching 'if'
        
        //Start HFE Tagging for seprerate analysis tasks
        if(AcceptJetforTag(clusMatch, AcceptedTrack)){
            AliEmcalJet::EFlavourTag tag=AliEmcalJet::kSig1;
            jet->AddFlavourTag(tag);
        }
            
        }  //Cluster Matched to Track loop
        //AcceptedTrack = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
      } //loop over tracks from Jet
    } // highest pt jet cut
  } // LOOP over JETS in event
  if(fGlobalQA>0) CheckClusTrackMatchingQA();
  return kTRUE;
}//________________________________________________________________________End of Main Task
//_________________________________________________________________________________________
void AliAnalysisTaskEmcalJetHF::Terminate(Option_t *)
{
  cout<<"#######################"<<endl;
  cout<<"#### Task Finished ####"<<endl;
  cout<<"#######################"<<endl;
}//________________________________________________________________________End of Terminate
//_________________________________________________________________________________________
void AliAnalysisTaskEmcalJetHF::ExtractMainPatch() {
  //Find main trigger
  if(!fTriggerPatchInfo) return;
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
      if (patch->IsRecalcGamma() == 1 )       fHistTriggerBitInfo->Fill(13);
    }
    if (patch->IsRecalcGamma() == 1 ){
      if(patch->GetPatchE()>emax) {
        fMaxPatch = patch;
        emax = patch->GetPatchE();
      }
    }
  }
}//___________________________________________End of Extrating Highest Energy Trigger Patch
//_________________________________________________________________________________________
void AliAnalysisTaskEmcalJetHF::CheckClusTrackMatchingQA()
{
  if(!fTracksCont || !fCaloClustersCont) return;
  Int_t useAOD = 1;
  if (dynamic_cast<AliAODEvent*>(InputEvent())) useAOD = 1;
  else useAOD = 0;
    
  // if we have ESD event, set up ESD object
  if(useAOD == 0){
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) {
        AliError(Form("ERROR: fESD not available\n"));
    }
  }
  // if we have AOD event, set up AOD object
  if(useAOD == 1){
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) {
        AliError(Form("ERROR: fAOD not available\n"));
    }
  }
  Double_t pQA1 = -999.;
  Double_t nSigmaElectron_TPC_QA1 = -999.;
  Double_t nSigmaElectron_TOF_QA1 = -999.;
  Double_t dEdxQA1 = -999.;
  Double_t deta = 999;
  Double_t dphi = 999;
  //Get closest cluster to track
  fTracksCont->ResetCurrentID();
  AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
  while(track) {
    //if(!track) continue;
    if(!useAOD){
      AliESDtrack *ESDtrackQA1 = static_cast<AliESDtrack*>(track);
      if(!ESDtrackQA1) continue;
      if(!fPIDResponse) continue;
      pQA1 = track->P();
      nSigmaElectron_TPC_QA1 = fPIDResponse->NumberOfSigmasTPC(ESDtrackQA1,AliPID::kElectron);
      nSigmaElectron_TOF_QA1 = fPIDResponse->NumberOfSigmasTOF(ESDtrackQA1,AliPID::kElectron);
      dEdxQA1 = ESDtrackQA1->GetTPCsignal();
    }
    if(useAOD){
      AliAODTrack *AODtrackQA1 = static_cast<AliAODTrack*>(track);
      if(!AODtrackQA1) continue;
      if(!fPIDResponse) continue;
      pQA1 = track->P();
      nSigmaElectron_TPC_QA1 = fPIDResponse->NumberOfSigmasTPC(AODtrackQA1,AliPID::kElectron);
      nSigmaElectron_TOF_QA1 = fPIDResponse->NumberOfSigmasTOF(AODtrackQA1,AliPID::kElectron);
      dEdxQA1 = AODtrackQA1->GetTPCsignal();
    }
    //Get matched cluster
    Int_t emc1 = track->GetEMCALcluster();
    if(fCaloClustersCont && emc1>=0) {
      AliVCluster *clusMatch = fCaloClustersCont->GetCluster(emc1);
      if(!clusMatch) continue;
      if(clusMatch) {
        Double_t ClusterE_QA1 = clusMatch->E();
        Double_t EovPQA1 = ClusterE_QA1/pQA1;
        Float_t pos_mc1[3];
        clusMatch->GetPosition(pos_mc1);  // Get cluster position
        TVector3 mc1(pos_mc1);
        Double_t HF_tracks3[11] = {track->Pt(), track->P() , track->Eta(), track->Phi(), EovPQA1, dEdxQA1 ,nSigmaElectron_TPC_QA1, nSigmaElectron_TOF_QA1, clusMatch->E(), mc1.PseudoRapidity(),mc1.Phi()};
        fhnTrackClusterQA->Fill(HF_tracks3);
      }//clus matching
    }//matched cluster
  track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
  }//track loop
}//________________________________________________________________________End Event QA PID
//_________________________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHF::AcceptMyJet(AliEmcalJet *jet) {
  //applies all jet cuts except pt
  if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax)) return 0;
  if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax)) return 0;
  if (jet->Area()<fAreacut) return 0;
  // prevents 0 area jets from sneaking by when area cut == 0
  if (jet->Area()==0) return 0;
  //exclude jets with extremely high pt tracks which are likely misreconstructed
  if(jet->MaxTrackPt()>100) return 0;
  
  //passed all above cuts for jet acceptance
  return 1;
}//____________________________________________________________________End of Jet QA checks
//_________________________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHF::AcceptJetforTag(AliVCluster *clusMatch, AliVTrack *AcceptedTrack) {
  Double_t EovP_tag;
  if (AcceptedTrack->Pt() < 5) return 0;
  EovP_tag = clusMatch->E() / AcceptedTrack->P();
  if(EovP_tag < 0.8) return 0;
  if(EovP_tag > 1.2) return 0;
  // passed rough cuts for HFE stats
  return 1;
}//__________________________________________________________________End of Jet HFE tagging
//_________________________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHF::GetCentBin(Double_t cent) const
{  // Get centrality bin.

  Int_t centbin = -1;
  if (cent>=0 && cent<10)
    centbin = 0; 
  else if (cent>=10 && cent<20)
    centbin = 1;
  else if (cent>=20 && cent<30)
    centbin = 2;
  else if (cent>=30 && cent<40)
    centbin = 3;
  else if (cent>=40 && cent<50)
    centbin = 4;
  else if (cent>=50 && cent<90)
    centbin = 5;
  return centbin;
}//________________________________________End of Centrality Binning for Histos and pp Runs
//_________________________________________________________________________________________
THnSparse* AliAnalysisTaskEmcalJetHF::NewTHnSparseDHF(const char* name, UInt_t entries)
{
  // generate new THnSparseD PID, axes are defined in GetDimParams()                                                                                                     
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
      GetDimParamsHF(i, label, nbins[c], xmin[c], xmax[c]);
      hnTitle += Form(";%s",label.Data());
      c++;
    }
    i++;
  }
  hnTitle += ";";
  return new THnSparseD(name, hnTitle.Data(), dim, nbins, xmin, xmax);
}//________________________________________________________________end of NewTHnSparseF PID
//_________________________________________________________________________________________
THnSparse* AliAnalysisTaskEmcalJetHF::NewTHnSparseDJetTrigger(const char* name, UInt_t entries)
{
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
}//_________________________________________________________end of NewTHnSparseF JetTrigger
//_________________________________________________________________________________________
void AliAnalysisTaskEmcalJetHF::GetDimParamsHF(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
  // stores label and binning of axis for THnSparse                                                                                                                      
  const Double_t pi = TMath::Pi();
  switch(iEntry){

  case 0:
    label = "V0 centrality (%)";
    nbins = 10;
    xmin = 0.;
    xmax = 100.;
    break;

  case 1:
    label = "Track p_{T}";
    nbins = 300;
    xmin = 0.;
    xmax = 75.;
    break;

  case 2:
    label = "Track p";
    nbins = 300;
    xmin = 0.;
    xmax = 75.;
    break;

  case 3:
    label = "Track Eta";
    nbins = 48;
    xmin = -1.2;
    xmax = 1.2;
    break;

  case 4:
    label = "Track Phi";
    nbins = 72;
    xmin = 0;
    xmax = 2*pi;
    break;

  case 5:
    label = "E/p of track";
    nbins = 400;
    xmin = 0;
    xmax = 4.0;
    break;

 case 6:
    label = "DCA xy";
    nbins = 20;
    xmin = -10;
    xmax =  10;
    break;

  case 7:
    label = "DCA z";
    nbins = 20;
    xmin = -10;
    xmax = 10;
    break;

  case 8:                                                                                                                                               
    label = "dEdX of track - TPC";
    nbins = 300;
    xmin = 0;
    xmax = 300;
    break;

  case 9:                                                                                                                                                
    label = "nSigma electron TPC";
    nbins = 50;
    xmin = -5;
    xmax = 5;
    break;

   case 10:
    label = "nSigma electron TOF";
    nbins = 50;
    xmin = -5;
    xmax = 5;
    break;

   case 11:
    label = "nSigma electron Emcal";
    nbins = 50;
    xmin = -5;
    xmax = 5;
    break;
      
  case 12:
    label = "Jet pT";
    nbins = 40;
    xmin  = 0;
    xmax  = 200;
    break;
      
  case 13:
    label = "Jet Phi";
    nbins = 72;
    xmin = 0;
    xmax = 2*pi;
    break;
      
  case 14:
    label = "Jet Eta";
    nbins = 24;
    xmin = -1.2;
    xmax = 1.2;
    break;
      
  case 15:
    label = "Cluster Energy";
    nbins = 150;
    xmin = 0;
    xmax = 15;
    break;
      
  case 16:
    label = "Cluster Eta";
    nbins = 24;
    xmin = -1.2;
    xmax =  1.2;
    break;
      
  case 17:
    label = "Cluster Phi";
    nbins = 72;
    xmin = 0;
    xmax = 2*pi;
    break;

  } // end of switch
} //______________________________________________________________end of getting dim-params
//_________________________________________________________________________________________
void AliAnalysisTaskEmcalJetHF::GetDimParamsJetTrigger(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
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
    label = "Trigger E";
    nbins = 100;
    xmin = 0.;
    xmax = 100.;
    break;
            
  case 6:
    label = "Trigger Eta";
    nbins = 24;
    xmin = -1.2;
    xmax =  1.2;
    break;
            
  case 7:
    label = "Trigger Phi";
    nbins = 72;
    xmin = 0;
    xmax = 2*pi;
    break;
            
  case 8:
    label = "Trigger ADC";
    nbins = 100;
    xmin = 0;
    xmax = 1500;
    break;
            
  case 9:
    label = "Cluster E";
    nbins = 15;
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
            
  case 12:
    label = "Track P";
    nbins = 100;
    xmin = 0.;
    xmax = 100.;
  break;
            
  case 13:
    label= "Track Pt";
    nbins = 100;
    xmin = 0.;
    xmax = 100.;
    break;
            
  case 14:
    label="Track Eta";
    nbins = 24;
    xmin = -1.2;
    xmax =  1.2;
    break;
            
  case 15:
    label = "Track Phi";
    nbins = 72;
    xmin = 0;
    xmax = 2*pi;
    break;
            
  case 16:
    label = "EovP_Trigger";
    nbins = 160;
    xmin = 0;
    xmax = 1.6;
    break;
            
  case 17:
    label = "nSigma electron TPC";
    nbins = 50;
    xmin = -5;
    xmax = 5;
          
  } // end of switch
} //______________________________________________________________end of getting dim-params




