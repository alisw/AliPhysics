/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
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
// AliAnalysisTaskSE for D0 candidates (2Prongs)
// and hadrons correlations
//
// Author:
// Fabio Colamaria, fabio.colamaria@ba.infn.it (correlations)
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSED0Correlations.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"
#include "AliHFOfflineCorrelator.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSED0Correlations)


//________________________________________________________________________
AliAnalysisTaskSED0Correlations::AliAnalysisTaskSED0Correlations():
AliAnalysisTaskSE(),
  fNPtBinsCorr(0), 
  fBinLimsCorr(),
  fPtThreshLow(),
  fPtThreshUp(), 
  fLSBLowLim(), 
  fLSBUppLim(), 
  fRSBLowLim(), 
  fRSBUppLim(),
  fDaughTrackID(),
  fDaughTrigNum(),
  fSoftPiTrackID(),
  fSoftPiTrigNum(),
  fEvents(0),
  fAlreadyFilled(kFALSE),
  fNtrigD(0),
  fNsoftPi(0),
  fOutputMass(0),
  fOutputCorr(0),
  fOutputStudy(0),
  fNentries(0), 
  fCutsD0(0),
  fCutsTracks(0),
  fCorrelatorTr(0),
  fCorrelatorKc(0),
  fCorrelatorK0(0),
  fReadMC(0),
  fRecoTr(kTRUE),
  fRecoD0(kTRUE),
  fSelEvType(kFALSE),
  fMixing(kFALSE),
  fCounter(0),
  fNPtBins(1),
  fFillOnlyD0D0bar(0),
  fIsSelectedCandidate(0),
  fSys(0),
  fEtaForCorrel(0),
  fIsRejectSDDClusters(0),
  fFillGlobal(kFALSE),
  fMultEv(0.),
  fzVtx(0.),
  fSoftPiCut(kTRUE),
  fMEAxisThresh(kFALSE),
  fKaonCorr(kFALSE),
  fSignLeft_LowPt(0),
  fSignRight_LowPt(0),
  fSignLeft_HighPt(0),
  fSignRight_HighPt(0),
  fPoolNum(0),
  fSpeed(kTRUE),
  fMergePools(kFALSE),
  fUseDeff(kTRUE),
  fUseTrackeff(kTRUE),
  fPtAssocLimit(1.),
  fMinDPt(2.),
  fFillTrees(kNoTrees),
  fFractAccME(100),
  fAODProtection(1),
  fBranchD(),
  fBranchTr(),
  fBranchDCutVars(),
  fTreeD(0x0),
  fTreeTr(0x0),
  fTrackArray(0x0),
  fTrackArrayFilled(kFALSE)     
{
  // Default constructor

}

//________________________________________________________________________
AliAnalysisTaskSED0Correlations::AliAnalysisTaskSED0Correlations(const char *name,AliRDHFCutsD0toKpi* cutsD0, AliHFAssociatedTrackCuts* cutsTrk):
  AliAnalysisTaskSE(name),
  fNPtBinsCorr(0),  
  fBinLimsCorr(),
  fPtThreshLow(),
  fPtThreshUp(), 
  fLSBLowLim(), 
  fLSBUppLim(), 
  fRSBLowLim(), 
  fRSBUppLim(),
  fDaughTrackID(),
  fDaughTrigNum(),
  fSoftPiTrackID(),
  fSoftPiTrigNum(),
  fEvents(0),
  fAlreadyFilled(kFALSE),
  fNtrigD(0),
  fNsoftPi(0),
  fOutputMass(0),
  fOutputCorr(0),
  fOutputStudy(0),
  fNentries(0),
  fCutsD0(0),
  fCutsTracks(cutsTrk),
  fCorrelatorTr(0),
  fCorrelatorKc(0),
  fCorrelatorK0(0),
  fReadMC(0),
  fRecoTr(kTRUE),
  fRecoD0(kTRUE),
  fSelEvType(kFALSE),
  fMixing(kFALSE),
  fCounter(0),
  fNPtBins(1),
  fFillOnlyD0D0bar(0),
  fIsSelectedCandidate(0),
  fSys(0),
  fEtaForCorrel(0),
  fIsRejectSDDClusters(0),
  fFillGlobal(kFALSE),
  fMultEv(0.),
  fzVtx(0.),
  fSoftPiCut(kTRUE),
  fMEAxisThresh(kFALSE),
  fKaonCorr(kFALSE),
  fSignLeft_LowPt(0),
  fSignRight_LowPt(0),
  fSignLeft_HighPt(0),
  fSignRight_HighPt(0),
  fPoolNum(0),
  fSpeed(kTRUE),
  fMergePools(kFALSE),
  fUseDeff(kTRUE),
  fUseTrackeff(kTRUE),
  fPtAssocLimit(1.),
  fMinDPt(2.),
  fFillTrees(kNoTrees),
  fFractAccME(100),
  fAODProtection(1),
  fBranchD(),
  fBranchTr(),
  fBranchDCutVars(),
  fTreeD(0x0),
  fTreeTr(0x0),
  fTrackArray(0x0),
  fTrackArrayFilled(kFALSE)        
{
  // Default constructor

  fNPtBins=cutsD0->GetNPtBins();
    
  fCutsD0=cutsD0;

  // Output slot #1 writes into a TList container (mass with cuts)
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes into a TH1F container (number of events)
  DefineOutput(2,TH1F::Class());  //My private output
  // Output slot #3 writes into a AliRDHFD0toKpi container (cuts)
  DefineOutput(3,AliRDHFCutsD0toKpi::Class());  //My private output
  // Output slot #4 writes Normalization Counter 
  DefineOutput(4,AliNormalizationCounter::Class());
  // Output slot #5 writes into a TList container (correl output)
  DefineOutput(5,TList::Class());  //My private output
  // Output slot #6 writes into a TList container (correl advanced)
  DefineOutput(6,TList::Class());  //My private output
  // Output slot #7 writes into a AliHFAssociatedTrackCuts container (cuts)
  DefineOutput(7,AliHFAssociatedTrackCuts::Class());  //My private output
  // Output slot #8 writes into a TTree (D0)
  DefineOutput(8,TTree::Class());  //My private output
  // Output slot #9 writes into a TTree (Tracks)
  DefineOutput(9,TTree::Class());  //My private output
}

//________________________________________________________________________
AliAnalysisTaskSED0Correlations::AliAnalysisTaskSED0Correlations(const AliAnalysisTaskSED0Correlations &source):
  AliAnalysisTaskSE(source),
  fNPtBinsCorr(source.fNPtBinsCorr), 
  fBinLimsCorr(source.fBinLimsCorr),
  fPtThreshLow(source.fPtThreshLow),
  fPtThreshUp(source.fPtThreshUp), 
  fLSBLowLim(source.fLSBLowLim), 
  fLSBUppLim(source.fLSBUppLim), 
  fRSBLowLim(source.fRSBLowLim), 
  fRSBUppLim(source.fRSBUppLim),
  fDaughTrackID(source.fDaughTrackID),
  fDaughTrigNum(source.fDaughTrigNum),
  fSoftPiTrackID(source.fSoftPiTrackID),
  fSoftPiTrigNum(source.fSoftPiTrigNum),
  fEvents(source.fEvents),
  fAlreadyFilled(source.fAlreadyFilled),
  fNtrigD(source.fNtrigD),
  fNsoftPi(source.fNsoftPi),
  fOutputMass(source.fOutputMass),
  fOutputCorr(source.fOutputCorr),
  fOutputStudy(source.fOutputStudy),
  fNentries(source.fNentries), 
  fCutsD0(source.fCutsD0),
  fCutsTracks(source.fCutsTracks),
  fCorrelatorTr(source.fCorrelatorTr),
  fCorrelatorKc(source.fCorrelatorKc),
  fCorrelatorK0(source.fCorrelatorK0),
  fReadMC(source.fReadMC),
  fRecoTr(source.fRecoTr),
  fRecoD0(source.fRecoD0),
  fSelEvType(source.fSelEvType),
  fMixing(source.fMixing),
  fCounter(source.fCounter),
  fNPtBins(source.fNPtBins),
  fFillOnlyD0D0bar(source.fFillOnlyD0D0bar),
  fIsSelectedCandidate(source.fIsSelectedCandidate),
  fSys(source.fSys),
  fEtaForCorrel(source.fEtaForCorrel),
  fIsRejectSDDClusters(source.fIsRejectSDDClusters),
  fFillGlobal(source.fFillGlobal),
  fMultEv(source.fMultEv),
  fzVtx(source.fzVtx),
  fSoftPiCut(source.fSoftPiCut),
  fMEAxisThresh(source.fMEAxisThresh),
  fKaonCorr(source.fKaonCorr),
  fSignLeft_LowPt(source.fSignLeft_LowPt),
  fSignRight_LowPt(source.fSignRight_LowPt),
  fSignLeft_HighPt(source.fSignLeft_HighPt),
  fSignRight_HighPt(source.fSignRight_HighPt),
  fPoolNum(source.fPoolNum),
  fSpeed(source.fSpeed),
  fMergePools(source.fMergePools),
  fUseDeff(source.fUseDeff),
  fUseTrackeff(source.fUseTrackeff),
  fPtAssocLimit(source.fPtAssocLimit),
  fMinDPt(source.fMinDPt),
  fFillTrees(source.fFillTrees),
  fFractAccME(source.fFractAccME),
  fAODProtection(source.fAODProtection),
  fBranchD(source.fBranchD),
  fBranchTr(source.fBranchTr),
  fBranchDCutVars(source.fBranchDCutVars), 
  fTreeD(source.fTreeD),
  fTreeTr(source.fTreeTr),  
  fTrackArray(source.fTrackArray),
  fTrackArrayFilled(source.fTrackArrayFilled)   
{
  // Copy constructor
}

//________________________________________________________________________
AliAnalysisTaskSED0Correlations::~AliAnalysisTaskSED0Correlations()
{
  if (fOutputMass) {
    delete fOutputMass;
    fOutputMass = 0;
  }
  if (fOutputCorr) {
    delete fOutputCorr;
    fOutputCorr = 0;
  }
  if (fOutputStudy) {
    delete fOutputStudy;
    fOutputStudy = 0;
  }
  if (fCutsD0) {
    delete fCutsD0;
    fCutsD0 = 0;
  }
  if (fNentries){
    delete fNentries;
    fNentries = 0;
  }
  if (fCorrelatorTr) {
    delete fCorrelatorTr;
    fCorrelatorTr = 0;
  }
  if (fCorrelatorKc) {
    delete fCorrelatorKc;
    fCorrelatorKc = 0;
  }
  if (fCorrelatorK0) {
    delete fCorrelatorK0;
    fCorrelatorK0 = 0;
  }
  if (fCounter){
    delete fCounter;
    fCounter=0;
  }
}  

//______________________________________________________________________________
AliAnalysisTaskSED0Correlations& AliAnalysisTaskSED0Correlations::operator=(const AliAnalysisTaskSED0Correlations& orig)
{
// Assignment
  if (&orig == this) return *this; //if address is the same (same object), returns itself

  AliAnalysisTaskSE::operator=(orig); //Uses the AliAnalysisTaskSE operator to assign the inherited part of the class
  fNPtBinsCorr = orig.fNPtBinsCorr; 
  fBinLimsCorr = orig.fBinLimsCorr;
  fPtThreshLow = orig.fPtThreshLow;
  fPtThreshUp = orig.fPtThreshUp; 
  fLSBLowLim = orig.fLSBLowLim; 
  fLSBUppLim = orig.fLSBUppLim; 
  fRSBLowLim = orig.fRSBLowLim;  
  fRSBUppLim = orig.fRSBUppLim; 
  fDaughTrackID = orig.fDaughTrackID;
  fDaughTrigNum = orig.fDaughTrigNum;
  fSoftPiTrackID = orig.fSoftPiTrackID;
  fSoftPiTrigNum = orig.fSoftPiTrigNum;
  fEvents = orig.fEvents;
  fAlreadyFilled = orig.fAlreadyFilled;
  fNtrigD = orig.fNtrigD;
  fNsoftPi = orig.fNsoftPi;
  fOutputMass = orig.fOutputMass;
  fOutputCorr = orig.fOutputCorr;
  fOutputStudy = orig.fOutputStudy;
  fNentries = orig.fNentries; 
  fCutsD0 = orig.fCutsD0;
  fCutsTracks = orig.fCutsTracks;
  fCorrelatorTr = orig.fCorrelatorTr;
  fCorrelatorKc = orig.fCorrelatorKc;
  fCorrelatorK0 = orig.fCorrelatorK0;
  fReadMC = orig.fReadMC;
  fRecoTr = orig.fRecoTr;
  fRecoD0 = orig.fRecoD0;
  fSelEvType = orig.fSelEvType;
  fMixing = orig.fMixing;
  fCounter = orig.fCounter;
  fNPtBins = orig.fNPtBins;
  fFillOnlyD0D0bar = orig.fFillOnlyD0D0bar;
  fIsSelectedCandidate = orig.fIsSelectedCandidate;
  fSys = orig.fSys;
  fEtaForCorrel = orig.fEtaForCorrel;
  fIsRejectSDDClusters = orig.fIsRejectSDDClusters;
  fFillGlobal = orig.fFillGlobal;
  fMultEv = orig.fMultEv;
  fzVtx = orig.fzVtx;
  fSoftPiCut = orig.fSoftPiCut;
  fMEAxisThresh = orig.fMEAxisThresh;
  fKaonCorr = orig.fKaonCorr;
  fSignLeft_LowPt = orig.fSignLeft_LowPt;
  fSignRight_LowPt = orig.fSignRight_LowPt;
  fSignLeft_HighPt = orig.fSignLeft_HighPt;
  fSignRight_HighPt = orig.fSignRight_HighPt;
  fPoolNum = orig.fKaonCorr;
  fSpeed = orig.fKaonCorr;   
  fMergePools = orig.fMergePools;
  fUseDeff = orig.fUseDeff;
  fUseTrackeff = orig.fUseTrackeff;
  fPtAssocLimit = orig.fPtAssocLimit;
  fMinDPt = orig.fMinDPt;
  fFillTrees = orig.fFillTrees;
  fFractAccME = orig.fFractAccME; 
  fAODProtection = orig.fAODProtection;
  fBranchD = orig.fBranchD;
  fBranchTr = orig.fBranchTr;
  fBranchDCutVars = orig.fBranchDCutVars;
  fTreeD = orig.fTreeD;
  fTreeTr = orig.fTreeTr;
  fTrackArray = orig.fTrackArray;      
  fTrackArrayFilled = orig.fTrackArrayFilled;
  
  return *this; //returns pointer of the class
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSED0Correlations::Init() \n");
  
  //Copy of cuts objects
  AliRDHFCutsD0toKpi* copyfCutsD0 = new AliRDHFCutsD0toKpi(*fCutsD0);
  const char* nameoutput=GetOutputSlot(3)->GetContainer()->GetName();
  copyfCutsD0->SetName(nameoutput);

  fTrackArray = new TObjArray();
  
  //needer to clear completely the objects inside with Clear() method
  // Post the data
  PostData(3,copyfCutsD0);
  PostData(7,fCutsTracks);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::UserCreateOutputObjects()
{

  if(fFillTrees==kFillTrees) {

    fBranchD = new AliHFCorrelationBranchD();
    fBranchTr = new AliHFCorrelationBranchTr();

    fTreeD = new TTree("fTreeD","TTree for D0 mesons");
    fTreeD->Branch("branchD",&fBranchD);

    fTreeTr = new TTree("fTreeTr","TTree for Associated Tracks");
    fTreeTr->Branch("branchTr",&fBranchTr);

    PostData(8,fTreeD);
    PostData(9,fTreeTr);
  }
  
  if(fFillTrees==kFillCutOptTree) {

    fBranchDCutVars = new AliD0hCutOptim();

    fTreeD = new TTree("fTreeD","TTree for D0 mesons - Vars for Cut Optimization");
    fTreeD->Branch("branchD",&fBranchDCutVars);

    PostData(8,fTreeD);
  }  

  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Correlations::UserCreateOutputObjects() \n");

  //HFCorrelator creation and definition
  fCorrelatorTr = new AliHFCorrelator("CorrelatorTr",fCutsTracks,fSys,fCutsD0);//fSys=0 use multiplicity, =1 use centrality
  fCorrelatorKc = new AliHFCorrelator("CorrelatorKc",fCutsTracks,fSys,fCutsD0);
  fCorrelatorK0 = new AliHFCorrelator("CorrelatorK0",fCutsTracks,fSys,fCutsD0);
  fCorrelatorTr->SetDeltaPhiInterval(-TMath::Pi()/2,3*TMath::Pi()/2);// set the Delta Phi Interval you want (in this case -0.5Pi to 1.5 Pi)
  fCorrelatorKc->SetDeltaPhiInterval(-TMath::Pi()/2,3*TMath::Pi()/2);
  fCorrelatorK0->SetDeltaPhiInterval(-TMath::Pi()/2,3*TMath::Pi()/2);
  fCorrelatorTr->SetEventMixing(fMixing);// sets the analysis on a single event (kFALSE) or mixed events (kTRUE)
  fCorrelatorKc->SetEventMixing(fMixing);
  fCorrelatorK0->SetEventMixing(fMixing);
  fCorrelatorTr->SetAssociatedParticleType(1);// set 1 for correlations with hadrons, 2 with kaons, 3 with KZeros
  fCorrelatorKc->SetAssociatedParticleType(2);// set 1 for correlations with hadrons, 2 with kaons, 3 with KZeros
  fCorrelatorK0->SetAssociatedParticleType(3);// set 1 for correlations with hadrons, 2 with kaons, 3 with KZeros
  fCorrelatorTr->SetApplyDisplacementCut(2); //0: don't calculate d0; 1: return d0; 2: return d0/d0err
  fCorrelatorKc->SetApplyDisplacementCut(2);
  fCorrelatorK0->SetApplyDisplacementCut(0);
  fCorrelatorTr->SetUseMC(fReadMC);// sets Montecarlo flag
  fCorrelatorKc->SetUseMC(fReadMC);
  fCorrelatorK0->SetUseMC(fReadMC);
  fCorrelatorTr->SetUseReco(fRecoTr);// sets (if MC analysis) wheter to analyze Reco or Kinem tracks
  fCorrelatorKc->SetUseReco(fRecoTr);
  fCorrelatorK0->SetUseReco(fRecoTr);
  fCorrelatorKc->SetPIDmode(2); //switch for K+/- PID option
  if(fMixing && fSoftPiCut) {
    fCorrelatorTr->SetStoreInfoSoftPiME(kTRUE);
    fCorrelatorKc->SetStoreInfoSoftPiME(kTRUE);
  }
  Bool_t pooldefTr = fCorrelatorTr->DefineEventPool();// method that defines the properties ot the event mixing (zVtx and Multipl. bins)
  Bool_t pooldefKc = fCorrelatorKc->DefineEventPool();// method that defines the properties ot the event mixing (zVtx and Multipl. bins)
  Bool_t pooldefK0 = fCorrelatorK0->DefineEventPool();// method that defines the properties ot the event mixing (zVtx and Multipl. bins)
  if(!pooldefTr) AliInfo("Warning:: Event pool not defined properly");
  if(!pooldefKc) AliInfo("Warning:: Event pool not defined properly");
  if(!pooldefK0) AliInfo("Warning:: Event pool not defined properly");

  // Several histograms are more conveniently managed in a TList
  fOutputMass = new TList();
  fOutputMass->SetOwner();
  fOutputMass->SetName("listMass");

  fOutputCorr = new TList();
  fOutputCorr->SetOwner();
  fOutputCorr->SetName("correlationslist");

  fOutputStudy = new TList();
  fOutputStudy->SetOwner();
  fOutputStudy->SetName("controlplots");

  TString nameMass=" ",nameSgn=" ", nameBkg=" ", nameRfl=" ",nameMassWg=" ",nameSgnWg=" ", nameBkgWg=" ", nameRflWg=" ";

//for origin c case (or for data)
  for(Int_t i=0;i<fCutsD0->GetNPtBins();i++){

    nameMass="histMass_"; if(fReadMC) nameMass+="c_";
    nameMass+=i;
    nameMassWg="histMass_WeigD0Eff_"; if(fReadMC) nameMassWg+="c_";
    nameMassWg+=i;
    nameSgn="histSgn_"; if(fReadMC) nameSgn+="c_";
    nameSgn+=i;
    nameSgnWg="histSgn_WeigD0Eff_"; if(fReadMC) nameSgnWg+="c_";
    nameSgnWg+=i;
    nameBkg="histBkg_"; if(fReadMC) nameBkg+="c_";
    nameBkg+=i;
    nameBkgWg="histBkg_WeigD0Eff_"; if(fReadMC) nameBkgWg+="c_";
    nameBkgWg+=i;
    nameRfl="histRfl_"; if(fReadMC) nameRfl+="c_";
    nameRfl+=i;
    nameRflWg="histRfl_WeigD0Eff_"; if(fReadMC) nameRflWg+="c_";
    nameRflWg+=i;

    //histograms of invariant mass distributions

    //MC signal
    if(fReadMC){
      TH1F* tmpSt = new TH1F(nameSgn.Data(), "D^{0} invariant mass c - MC; M [GeV]; Entries",150,1.5648,2.1648);
      TH1F* tmpStWg = new TH1F(nameSgnWg.Data(), "D^{0} invariant mass c - MC; M [GeV] - weight 1/D0eff; Entries",150,1.5648,2.1648);
      tmpSt->Sumw2();
      tmpStWg->Sumw2();

      //Reflection: histo filled with D0Mass which pass the cut (also) as D0bar and with D0bar which pass (also) the cut as D0
      TH1F* tmpRt = new TH1F(nameRfl.Data(), "Reflected signal invariant mass c - MC; M [GeV]; Entries",150,1.5648,2.1648);
      TH1F* tmpRtWg = new TH1F(nameRflWg.Data(), "Reflected signal invariant mass c - MC - weight 1/D0eff; M [GeV]; Entries",150,1.5648,2.1648);
      TH1F* tmpBt = new TH1F(nameBkg.Data(), "Background invariant mass c - MC; M [GeV]; Entries",150,1.5648,2.1648);
      TH1F* tmpBtWg = new TH1F(nameBkgWg.Data(), "Background invariant mass c - MC - weight 1/D0eff; M [GeV]; Entries",150,1.5648,2.1648);
      tmpBt->Sumw2();
      tmpBtWg->Sumw2();
      tmpRt->Sumw2();
      tmpRtWg->Sumw2();
      fOutputMass->Add(tmpSt);
      fOutputMass->Add(tmpStWg);
      fOutputMass->Add(tmpRt);
      fOutputMass->Add(tmpRtWg);
      fOutputMass->Add(tmpBt);
      fOutputMass->Add(tmpBtWg);
    }

    //mass
    TH1F* tmpMt = new TH1F(nameMass.Data(),"D^{0} invariant mass c; M [GeV]; Entries",150,1.5648,2.1648);
    tmpMt->Sumw2();
    fOutputMass->Add(tmpMt);
    //mass weighted by 1/D0eff
    TH1F* tmpMtwg = new TH1F(nameMassWg.Data(),"D^{0} invariant mass c - weight 1/D0eff; M [GeV]; Entries",150,1.5648,2.1648);
    tmpMtwg->Sumw2();
    fOutputMass->Add(tmpMtwg);
    
    if(fFillTrees>0) { //multi-histo for mass, pT, centrality for offline code (use in place of TH1F to select centrality and change pT bins offline)
      nameMass="histMass2D_";  nameMass+=i;
      TH2F* hMass2D = new TH2F(nameMass.Data(),"Mass histogram vs centrality; Entries",150,1.5648,2.1648,100,0.,100.);
      hMass2D->Sumw2();
      fOutputMass->Add(hMass2D);

      nameMass="histMass2D_WeigD0Eff_";  nameMass+=i;      
      TH2F* hMass2DW = new TH2F(nameMass.Data(),"Mass histogram vs centrality - weight 1/D0eff; Entries",150,1.5648,2.1648,100,0.,100.);
      hMass2DW->Sumw2();
      fOutputMass->Add(hMass2DW);      
    }
    
  }

//for origin b case (no Bkg and Mass histos, here for weights you should use c+b efficiencies, while on data (on MC they're useless))
  for(Int_t i=0;i<fCutsD0->GetNPtBins();i++){

    nameSgn="histSgn_b_";
    nameSgn+=i;
    nameSgnWg="histSgn_WeigD0Eff_b_";
    nameSgnWg+=i;
    nameRfl="histRfl_b_";
    nameRfl+=i;
    nameRflWg="histRfl_WeigD0Eff_b_";
    nameRflWg+=i;

    //histograms of invariant mass distributions

    //MC signal
    if(fReadMC){
      TH1F* tmpSt = new TH1F(nameSgn.Data(), "D^{0} invariant mass b - MC; M [GeV]; Entries",150,1.5648,2.1648);
      TH1F* tmpStWg = new TH1F(nameSgnWg.Data(), "D^{0} invariant mass b - MC; M [GeV] - weight 1/D0eff; Entries",150,1.5648,2.1648);
      tmpSt->Sumw2();
      tmpStWg->Sumw2();

      //Reflection: histo filled with D0Mass which pass the cut (also) as D0bar and with D0bar which pass (also) the cut as D0
      TH1F* tmpRt = new TH1F(nameRfl.Data(), "Reflected signal invariant mass b - MC; M [GeV]; Entries",150,1.5648,2.1648);
      TH1F* tmpRtWg = new TH1F(nameRflWg.Data(), "Reflected signal invariant mass b - MC - weight 1/D0eff; M [GeV]; Entries",150,1.5648,2.1648);
      tmpRt->Sumw2();
      tmpRtWg->Sumw2();
      fOutputMass->Add(tmpSt);
      fOutputMass->Add(tmpStWg);
      fOutputMass->Add(tmpRt);
      fOutputMass->Add(tmpRtWg);
    }
  }

  const char* nameoutput=GetOutputSlot(2)->GetContainer()->GetName();

  fNentries=new TH1F(nameoutput, "Control plot", 20,-0.5,19.5);

  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nEventsSelected");
  fNentries->GetXaxis()->SetBinLabel(3,"nEventsGoodVtxPrim");
  fNentries->GetXaxis()->SetBinLabel(4,"mismatch AOD/dAOD");
  fNentries->GetXaxis()->SetBinLabel(5,"REJ: no prim vtx");
  fNentries->GetXaxis()->SetBinLabel(6,"REJ: Pile-up");  
  fNentries->GetXaxis()->SetBinLabel(7,"REJ: centrality");
  fNentries->GetXaxis()->SetBinLabel(8,"REJ: cent wrong ext");
  fNentries->GetXaxis()->SetBinLabel(9,"REJ: cent flatten");
  fNentries->GetXaxis()->SetBinLabel(10,"REJ: trigger class");
  fNentries->GetXaxis()->SetBinLabel(11,"REJ: trigger mask");
  fNentries->GetXaxis()->SetBinLabel(12,"REJ: zVtx>10cm");
  fNentries->GetXaxis()->SetBinLabel(13,"N. of 0SMH");    
  if(fIsRejectSDDClusters) fNentries->GetXaxis()->SetBinLabel(14,"SDD-Cls Rej");
  if(fReadMC) fNentries->GetXaxis()->SetBinLabel(15,"nEvsWithProdMech");
  fNentries->GetXaxis()->SetBinLabel(16,"D0 failed to be filled");
  fReadMC ? fNentries->GetXaxis()->SetBinLabel(17,"nTrueD0Selected(MC)") : fNentries->GetXaxis()->SetBinLabel(17,"Dstar<-D0");
  fNentries->GetXaxis()->SetBinLabel(18,"ptbin = -1");
  if(fSys==0) fNentries->GetXaxis()->SetBinLabel(19,"nCandSel(QualTr)");
  fNentries->GetXaxis()->SetBinLabel(20,"nCandSel(Cuts)");

  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(4)->GetContainer()->GetName()));
  fCounter->Init();

  CreateCorrelationsObjs(); //creates histos for correlations analysis

  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fNentries);
  PostData(4,fCounter);
  PostData(5,fOutputCorr);
  PostData(6,fOutputStudy);
  PostData(8,fTreeD);
  PostData(9,fTreeTr);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  //cout<<"I'm in UserExec"<<endl;


  //cuts order
  //     printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
  //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
  //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
  //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
  //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
  //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
  //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
  //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
  //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  fEvents++;

  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fNentries->Fill(3);
      return;
    }
  }

  TString bname="D0toKpi";

  TClonesArray *inputArray=0;

  fMultEv = 0.; //reset event multiplicity
  fzVtx = 0.; //reset event multiplicity
  fPoolNum = 0; //reset event pool

  fDaughTrackID.clear(); //removes daugher IDs from previous event
  fDaughTrigNum.clear(); //removes daugher trigger matchings from previous event
  fSoftPiTrackID.clear(); //removes soft pion IDs from previous event
  fSoftPiTrigNum.clear(); //removes soft pion trtigger matchings from previous event
  fTrackArray->Clear(); //removes associated tracks selected from previous event
  fTrackArrayFilled = kFALSE; //associated track array is now not filled
  
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent* aodFromExt = ext->GetAOD();
      inputArray=(TClonesArray*)aodFromExt->GetList()->FindObject(bname.Data());
    }
  } else if(aod) {
    inputArray=(TClonesArray*)aod->GetList()->FindObject(bname.Data());
  }

  if(!inputArray || !aod) {
    printf("AliAnalysisTaskSED0Correlations::UserExec: input branch not found!\n");
    return;
  }

  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;

  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if(fReadMC) {
    // load MC particles
    mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      printf("AliAnalysisTaskSED0Correlations::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSED0Correlations::UserExec: MC header branch not found!\n");
      return;
    }
  }

  //histogram filled with 1 for every AOD
  fNentries->Fill(0);
  fCounter->StoreEvent(aod,fCutsD0,fReadMC); 

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNentries->Fill(12);

  //Call IsEventSelected only for Reco! (and Data of course)
  if(fRecoD0 && !fCutsD0->IsEventSelected(aod)) {
    if(fCutsD0->GetWhyRejection()==0) fNentries->Fill(4); // no prim vertex
    if(fCutsD0->GetWhyRejection()==1) fNentries->Fill(5); // rejected for pileup
    if(fCutsD0->GetWhyRejection()==2) fNentries->Fill(6); // rejected for centrality
    if(fCutsD0->GetWhyRejection()==3) fNentries->Fill(7); // rejected for centrality wrong extim
    if(fCutsD0->GetWhyRejection()==4) fNentries->Fill(8); // rejected for centrality flattening
    if(fCutsD0->GetWhyRejection()==5) fNentries->Fill(9); // rejected for trigger class
    if(fCutsD0->GetWhyRejection()==6) fNentries->Fill(11); // rejected for vtx outside 10
    if(fCutsD0->GetWhyRejection()==7) fNentries->Fill(10); // rejected for trigger mask
    return;
  }

  //On Kine, instead of IsEventSelected just select on zVtx and trigger mask in pPb
  if(!fRecoD0) {
    Double_t zVtxMC = mcHeader->GetVtxZ();
    if(TMath::Abs(zVtxMC)>10) return;
    if(aod->GetTriggerMask()==0 && (aod->GetRunNumber()>=195344 && aod->GetRunNumber()<=195677)) return;
  }
  
  fNentries->Fill(1); //event selected after selection

  //Setting PIDResponse for associated tracks
  fCorrelatorTr->SetPidAssociated();
  fCorrelatorKc->SetPidAssociated();
  fCorrelatorK0->SetPidAssociated();

  //Selection on production type (MC)
  if(fReadMC && fSelEvType){ 

    Bool_t isMCeventgood = kFALSE;
            
    Int_t eventType = mcHeader->GetEventType();
    Int_t NMCevents = fCutsTracks->GetNofMCEventType();
               
    for(Int_t k=0; k<NMCevents; k++){
      Int_t * MCEventType = fCutsTracks->GetMCEventType();
          
      if(eventType == MCEventType[k]) isMCeventgood= kTRUE;
      ((TH1D*)fOutputStudy->FindObject("EventTypeMC"))->Fill(eventType);
    }
                
    if(NMCevents && !isMCeventgood){
      if(fDebug > 2) std::cout << "The MC event " << eventType << " not interesting for this analysis: skipping" << std::endl;
      return; 
    }
    fNentries->Fill(14); //event with particular production type                
  
  } //end of selection


  // Check the Nb of SDD clusters
  if (fIsRejectSDDClusters) { 
    Bool_t skipEvent = kFALSE;
    Int_t ntracks = 0;
    if (aod) ntracks = aod->GetNumberOfTracks();
    for(Int_t itrack=0; itrack<ntracks; itrack++) { // loop on tacks
      //    ... get the track
      AliAODTrack * track = (AliAODTrack*)aod->GetTrack(itrack);
      if(!track) {
	AliWarning("Error in casting to AOD track. Not a standard AOD?");
        continue;
      }
      if(TESTBIT(track->GetITSClusterMap(),2) || TESTBIT(track->GetITSClusterMap(),3) ){
	skipEvent=kTRUE;
	fNentries->Fill(13);
	break;
      }
    }
    if (skipEvent) return;
  }
  
  //HFCorrelators initialization (for this event)
  fCorrelatorTr->SetAODEvent(aod); // set the AOD event from which you are processing
  fCorrelatorKc->SetAODEvent(aod);
  fCorrelatorK0->SetAODEvent(aod);
  Bool_t correlatorONTr = fCorrelatorTr->Initialize(); // initialize the pool for event mixing
  Bool_t correlatorONKc = fCorrelatorKc->Initialize();
  Bool_t correlatorONK0 = fCorrelatorK0->Initialize();
  if(!correlatorONTr) {AliInfo("AliHFCorrelator (tracks) didn't initialize the pool correctly or processed a bad event"); return;}
  if(!correlatorONKc) {AliInfo("AliHFCorrelator (charged K) didn't initialize the pool correctly or processed a bad event"); return;}
  if(!correlatorONK0) {AliInfo("AliHFCorrelator (K0) didn't initialize the pool correctly or processed a bad event"); return;}

  if(fReadMC) {
    fCorrelatorTr->SetMCArray(mcArray); // set the TClonesArray *fmcArray for analysis on monte carlo
    fCorrelatorKc->SetMCArray(mcArray);
    fCorrelatorK0->SetMCArray(mcArray);
  }

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();

  //Pool definition
  Double_t MultipOrCent = fCorrelatorTr->GetCentrality();
  Double_t zVtxPosition = vtx1->GetZ();
  fzVtx = zVtxPosition;
  if(!fMergePools) fPoolNum = fCutsTracks->GetPoolBin(MultipOrCent, zVtxPosition);
  
  //vtx1->Print();
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) {
    fNentries->Fill(2);
  }

  //Reset flag for tracks distributions fill and counter of D0 triggers and of Soft pions
  fAlreadyFilled=kFALSE;
  fNtrigD=0;
  fNsoftPi=0;

  //***** Loop over D0 candidates *****
  Int_t nInD0toKpi = inputArray->GetEntriesFast();
  if(fDebug>2) printf("Number of D0->Kpi: %d\n",nInD0toKpi);

  if(fFillGlobal) { //loop on V0 and tracks for each event, to fill Pt distr. and InvMass distr.
/*
    TClonesArray *v0array = (TClonesArray*)aod->GetList()->FindObject("v0s");
    Int_t pdgCodes[2] = {211,211};
    Int_t idArrayV0[v0array->GetEntriesFast()][2];
    for(int iV0=0; iV0<v0array->GetEntriesFast(); iV0++) {
      for(int j=0; j<2; j++) {idArrayV0[iV0][j]=-2;}
      AliAODv0 *v0 = (AliAODv0*)v0array->UncheckedAt(iV0);
      if(SelectV0(v0,vtx1,2,idArrayV0)) { //option 2 = for mass inv plots only
        if(fReadMC && fRecoTr && (v0->MatchToMC(310,mcArray,2,pdgCodes)<0)) continue; //310 = K0s, 311 = K0 generico!!
        ((TH2F*)fOutputStudy->FindObject("hK0MassInv"))->Fill(v0->MassK0Short(),v0->Pt()); //invariant mass plot
        ((TH1F*)fOutputStudy->FindObject("hist_Pt_K0_AllEv"))->Fill(v0->Pt()); //pT distribution (in all events), K0 case
      }
    }

    for(Int_t itrack=0; itrack<aod->GetNTracks(); itrack++) { // loop on tacks
      AliAODTrack * track = aod->GetTrack(itrack);
      if(!track) {
	AliWarning("Error in casting to AOD track. Not a standard AOD?");
        continue;
      }
      //rejection of tracks
      if(track->GetID() < 0) continue; //discard negative ID tracks
      if(track->Pt() < fPtThreshLow.at(0) || track->Pt() > fPtThreshUp.at(0)) continue; //discard tracks outside pt range for hadrons/K
      if(!fCutsTracks->IsHadronSelected(track) || !fCutsTracks->CheckHadronKinematic(track->Pt(),0.1)) continue; //0.1 = dummy (d0 value, no cut on it for me)
      //pT distribution (in all events), charg and hadr cases
      ((TH1F*)fOutputStudy->FindObject("hist_Pt_Charg_AllEv"))->Fill(track->Pt()); 
      if(fCutsTracks->CheckKaonCompatibility(track,kFALSE,0,2)) ((TH1F*)fOutputStudy->FindObject("hist_Pt_Kcharg_AllEv"))->Fill(track->Pt());
    }
*/
  } //end of loops for global plot fill

  Int_t nSelectedloose=0,nSelectedtight=0;  

  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
  
  //Fill Event Multiplicity (needed only in Reco)
  fMultEv = (Double_t)(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.));

  //Fill control plots for event centrality and zVtx
  if(fCutsD0->GetUseCentrality()) ((TH1F*)fOutputStudy->FindObject("hCentralEvts"))->Fill(fCutsD0->GetCentrality(aod));
  ((TH1F*)fOutputStudy->FindObject("hZvtxEvts"))->Fill(vtx1->GetZ());

  //RecoD0 case ************************************************
  if(fRecoD0) {

    for (Int_t iD0toKpi = 0; iD0toKpi < nInD0toKpi; iD0toKpi++) {
      AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)inputArray->UncheckedAt(iD0toKpi);

      if(!(vHF->FillRecoCand(aod,d))) {//Fill the data members of the candidate only if they are empty.   
        fNentries->Fill(15); //monitor how often this fails 
        continue;
      }

      if(d->Pt() < fMinDPt) continue; //to save time and merging memory...

      if(d->GetSelectionMap()) if(!d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)){
  	fNentries->Fill(16);
  	continue; //skip the D0 from Dstar  
      }
    
      if(fCutsD0->IsInFiducialAcceptance(d->Pt(),d->Y(421))) {
        nSelectedloose++;
        nSelectedtight++;      
        if(fSys==0){
  	  if(fCutsD0->IsSelected(d,AliRDHFCuts::kTracks,aod)) fNentries->Fill(18);       
        }  
        Int_t ptbin=fCutsD0->PtBin(d->Pt());
        if(ptbin==-1) {fNentries->Fill(17); continue;} //out of bounds

        fIsSelectedCandidate=fCutsD0->IsSelected(d,AliRDHFCuts::kAll,aod); //D0 selected
        if(!fIsSelectedCandidate) continue;

        //D0 infos
        Double_t phiD0 = fCorrelatorTr->SetCorrectPhiRange(d->Phi());
                 phiD0 = fCorrelatorKc->SetCorrectPhiRange(d->Phi());  //bad usage, but returns a Double_t...
                 phiD0 = fCorrelatorK0->SetCorrectPhiRange(d->Phi());
        fCorrelatorTr->SetTriggerParticleProperties(d->Pt(),phiD0,d->Eta()); // sets the parameters of the trigger particles that are needed
        fCorrelatorKc->SetTriggerParticleProperties(d->Pt(),phiD0,d->Eta());
        fCorrelatorK0->SetTriggerParticleProperties(d->Pt(),phiD0,d->Eta());
        fCorrelatorTr->SetD0Properties(d,fIsSelectedCandidate); //sets special properties for D0
        fCorrelatorKc->SetD0Properties(d,fIsSelectedCandidate);
        fCorrelatorK0->SetD0Properties(d,fIsSelectedCandidate);

        if(!fReadMC) {
          if (TMath::Abs(d->Eta())<fEtaForCorrel) {
	    if(!fAlreadyFilled && !fFillTrees) ((TH1F*)fOutputStudy->FindObject(Form("hEvtsPerPool_%d",ptbin)))->Fill(fPoolNum+0.5);			
            if(!fMixing && !fAlreadyFilled) {
 	      ((TH1F*)fOutputStudy->FindObject("hZvtx"))->Fill(vtx1->GetZ());
	      ((TH1F*)fOutputStudy->FindObject(Form("hMultiplEvt_Bin%d",ptbin)))->Fill(fMultEv);
            }
	    if(fFillTrees==kNoTrees) CalculateCorrelations(d); //correlations on real data
	    if(fFillTrees==kFillTrees) FillTreeD0(d,aod); //only for offline correlations
	    if(fFillTrees==kFillCutOptTree) FillTreeD0ForCutOptim(d,aod); //only for offline correlations to optimize D0 cut variables
	  }
        } else { //correlations on MC -> association of selected D0 to MCinfo with MCtruth
          if (TMath::Abs(d->Eta())<fEtaForCorrel) {
            Int_t pdgDgD0toKpi[2]={321,211};
    	    Int_t labD0 = d->MatchToMC(421,mcArray,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not
            if (labD0>-1) {
  	      if(!fAlreadyFilled && !fFillTrees) ((TH1F*)fOutputStudy->FindObject(Form("hEvtsPerPool_%d",ptbin)))->Fill(fPoolNum+0.5);
              if(!fMixing && !fAlreadyFilled) {
		((TH1F*)fOutputStudy->FindObject("hZvtx"))->Fill(vtx1->GetZ());
                ((TH1F*)fOutputStudy->FindObject(Form("hMultiplEvt_Bin%d",ptbin)))->Fill(fMultEv); //Fill multiplicity histo
              }
	      CalculateCorrelations(d,labD0,mcArray);
	    }
          }
        }

        FillMassHists(d,mcArray,fCutsD0,fOutputMass,aod);
      }
    }
  }
  //End RecoD0 case ************************************************
  
  //MCKineD0 case ************************************************
  if(fReadMC && !fRecoD0) {

    for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) { //Loop over all the tracks of MCArray
      AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
      if (!mcPart) {
        AliWarning("Particle not found in tree, skipping"); 
        continue;
      } 
  
      if(TMath::Abs(mcPart->GetPdgCode()) == 421){  // THIS IS A D0
        if (fCutsD0->IsInFiducialAcceptance(mcPart->Pt(),mcPart->Y()) ) {
          nSelectedloose++;
          nSelectedtight++;      

          //Removal of cases in which D0 decay is not in Kpi!
	  if(mcPart->GetNDaughters()!=2) continue;
	  AliAODMCParticle* mcDau1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcPart->GetDaughter(0)));
	  AliAODMCParticle* mcDau2 = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcPart->GetDaughter(1)));
	  if(!mcDau1 || !mcDau2) continue;
	  Int_t pdg1 = TMath::Abs(mcDau1->GetPdgCode());
	  Int_t pdg2 = TMath::Abs(mcDau2->GetPdgCode());
          if(!((pdg1 == 211 && pdg2 == 321) || (pdg2 == 211 && pdg1 == 321))) continue;
          if(TMath::Abs(mcDau1->Eta())>0.8||TMath::Abs(mcDau2->Eta())>0.8) continue;
            //Check momentum conservation (to exclude 4-prong decays with tracks outside y=1.5)
            Double_t p1[3]  = {mcDau1->Px(),mcDau1->Py(),mcDau1->Pz()};
            Double_t p2[3]  = {mcDau2->Px(),mcDau2->Py(),mcDau2->Pz()};
            Double_t pD0[3] = {mcPart->Px(),mcPart->Py(),mcPart->Pz()};
            if(TMath::Abs( (p1[0]+p2[0]-pD0[0])*(p1[0]+p2[0]-pD0[0]) + (p1[1]+p2[1]-pD0[1])*(p1[1]+p2[1]-pD0[1]) + (p1[2]+p2[2]-pD0[2])*(p1[2]+p2[2]-pD0[2]) )>0.1) continue;

          if(fSys==0) fNentries->Fill(18);
          Int_t ptbin=fCutsD0->PtBin(mcPart->Pt());
          if(ptbin==-1) {fNentries->Fill(17); continue;} //out of bounds  
  
          //D0 infos
          Double_t phiD0 = fCorrelatorTr->SetCorrectPhiRange(mcPart->Phi());
                   phiD0 = fCorrelatorKc->SetCorrectPhiRange(mcPart->Phi());  //bad usage, but returns a Double_t...
                   phiD0 = fCorrelatorK0->SetCorrectPhiRange(mcPart->Phi());
          fCorrelatorTr->SetTriggerParticleProperties(mcPart->Pt(),phiD0,mcPart->Eta()); // sets the parameters of the trigger particles that are needed
          fCorrelatorKc->SetTriggerParticleProperties(mcPart->Pt(),phiD0,mcPart->Eta());
          fCorrelatorK0->SetTriggerParticleProperties(mcPart->Pt(),phiD0,mcPart->Eta());
          //fCorrelatorTr->SetD0Properties(mcPart,fIsSelectedCandidate); //needed for D* soft pions rejection, useless in MCKine
          //fCorrelatorKc->SetD0Properties(mcPart,fIsSelectedCandidate);
          //fCorrelatorK0->SetD0Properties(mcPart,fIsSelectedCandidate);
  
          if (TMath::Abs(mcPart->Eta())<fEtaForCorrel) {
  
            //Removal of D0 from D* feeddown! This solves also the problem of soft pions, now excluded
         /*   Int_t mother = mcPart->GetMother();
  	    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(mcArray->At(mother));
            if(!mcMoth) continue;
	    if(TMath::Abs(mcMoth->GetPdgCode())==413) continue;*/

            if (mcPart->GetPdgCode()==421) fIsSelectedCandidate = 1;
    	    else fIsSelectedCandidate = 2;

	    TString fillthis="histSgn_"; 
	    if(CheckD0Origin(mcArray,mcPart)==4) fillthis+="c_";
	    else if(CheckD0Origin(mcArray,mcPart)==5) fillthis+="b_";
            else continue;
            fillthis+=ptbin;
	    ((TH1F*)(fOutputMass->FindObject(fillthis)))->Fill(1.864);
	  
            CalculateCorrelationsMCKine(mcPart,mcArray);
            if(!fMixing) ((TH1F*)fOutputStudy->FindObject(Form("hMultiplEvt_Bin%d",ptbin)))->Fill(fMultEv); //Fill multiplicity histo
          }
        }
      }
    }
  }
  //End MCKineD0 case ************************************************

  if(fMixing && fFillTrees!=kFillTrees /* && fAlreadyFilled*/) { // update the pool for Event Mixing, if: enabled,  event is ok, at least a SelD0 found! (fAlreadyFilled's role!)
    Bool_t updatedTr = fCorrelatorTr->PoolUpdate();
    Bool_t updatedKc = fCorrelatorKc->PoolUpdate();
    Bool_t updatedK0 = fCorrelatorK0->PoolUpdate();
    if(!updatedTr || !updatedKc || !updatedK0) AliInfo("Pool was not updated");
  }
  if(fFillTrees==kFillTrees && fAlreadyFilled) FillTreeTracks(aod);
  
  fCounter->StoreCandidates(aod,nSelectedloose,kTRUE);  
  fCounter->StoreCandidates(aod,nSelectedtight,kFALSE);  
  delete vHF;

  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fNentries);
  PostData(4,fCounter);
  PostData(5,fOutputCorr);
  PostData(6,fOutputStudy);
  if(fFillTrees!=kNoTrees) PostData(8,fTreeD); //fill in case kFillTrees or kFillCutOptTree
  if(fFillTrees==kFillTrees) PostData(9,fTreeTr); //fill only in case kFillTrees
  
  return;
}

//____________________________________________________________________________
void AliAnalysisTaskSED0Correlations::FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi* cuts, TList *listout, AliAODEvent *aod) {
  //
  // function used in UserExec to fill mass histograms:
  //
  if (!fIsSelectedCandidate) return;

  if(fDebug>2)  cout<<"Candidate selected"<<endl;

  Double_t invmassD0 = part->InvMassD0(), invmassD0bar = part->InvMassD0bar();
  Int_t ptbin = cuts->PtBin(part->Pt());
  
  TString fillthis="";
  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t labD0=-1;
  if (fReadMC) labD0 = part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)

  //count candidates selected by cuts
  fNentries->Fill(19);
  //count true D0 selected by cuts
  if (fReadMC && labD0>=0) fNentries->Fill(16);

  if ((fIsSelectedCandidate==1 || fIsSelectedCandidate==3) && fFillOnlyD0D0bar<2) { //D0

    if(fReadMC){ //on MC
      if(labD0>=0 && CheckD0Origin(arrMC,(AliAODMCParticle*)arrMC->At(labD0))==4) {
  	AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();
	if (pdgD0==421){ //D0
	  fillthis="histSgn_c_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
          fillthis="histSgn_WeigD0Eff_c_";
          fillthis+=ptbin;
          Double_t effD0 = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
          if(!fUseDeff) effD0=1.;
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0,1./effD0);
  	} else{ //it was a D0bar
	  fillthis="histRfl_c_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
          fillthis="histRfl_WeigD0Eff_c_";
          fillthis+=ptbin;
          Double_t effD0 = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
          if(!fUseDeff) effD0=1.;          
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0,1./effD0);
  	}
      } else if(labD0>=0 && CheckD0Origin(arrMC,(AliAODMCParticle*)arrMC->At(labD0))==5) {
  	  AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	  Int_t pdgD0 = partD0->GetPdgCode();
	  if (pdgD0==421){ //D0
	    fillthis="histSgn_b_";
	    fillthis+=ptbin;
	    ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
            fillthis="histSgn_WeigD0Eff_b_";
            fillthis+=ptbin;
            Double_t effD0 = fCutsTracks->GetTrigWeightB(part->Pt(),fMultEv);
            if(!fUseDeff) effD0=1.;            
            ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0,1./effD0);
  	  } else{ //it was a D0bar
	    fillthis="histRfl_b_";
	    fillthis+=ptbin;
	    ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
            fillthis="histRfl_WeigD0Eff_b_";
            fillthis+=ptbin;
            Double_t effD0 = fCutsTracks->GetTrigWeightB(part->Pt(),fMultEv);
            if(!fUseDeff) effD0=1.;          
            ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0,1./effD0);
	  }
      } else {//background
  	fillthis="histBkg_c_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
        fillthis="histBkg_WeigD0Eff_c_";
        fillthis+=ptbin;
        Double_t effD0 = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
        if(!fUseDeff) effD0=1.; 
        ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0,1./effD0);
      }
    }else{ //on data
      fillthis="histMass_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      fillthis="histMass_WeigD0Eff_";
      fillthis+=ptbin;
      Double_t effD0 = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
      if(!fUseDeff) effD0=1.; 
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0,1./effD0);
      if(fFillTrees>0) {
	Double_t centFill = 0.;
	if(fCutsD0->GetUseCentrality()) centFill = fCutsD0->GetCentrality(aod);
        ((TH2F*)(listout->FindObject(Form("histMass2D_%d",ptbin))))->Fill(invmassD0,centFill);
        ((TH2F*)(listout->FindObject(Form("histMass2D_WeigD0Eff_%d",ptbin))))->Fill(invmassD0,centFill,1./effD0);
      }
      
    }
     
  }
  if (fIsSelectedCandidate>1 && (fFillOnlyD0D0bar==0 || fFillOnlyD0D0bar==2)) { //D0bar

    if(fReadMC){ //on MC
      if(labD0>=0 && CheckD0Origin(arrMC,(AliAODMCParticle*)arrMC->At(labD0))==4) {
  	AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();
	if (pdgD0==-421){ //D0
	  fillthis="histSgn_c_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
          fillthis="histSgn_WeigD0Eff_c_";
          fillthis+=ptbin;
          Double_t effD0 = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
          if(!fUseDeff) effD0=1.; 
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar,1./effD0);
  	} else{ //it was a D0bar
	  fillthis="histRfl_c_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
          fillthis="histRfl_WeigD0Eff_c_";
          fillthis+=ptbin;
          Double_t effD0 = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
          if(!fUseDeff) effD0=1.; 
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar,1./effD0);
  	}
      } else if(labD0>=0 && CheckD0Origin(arrMC,(AliAODMCParticle*)arrMC->At(labD0))==5) {
  	  AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	  Int_t pdgD0 = partD0->GetPdgCode();
	  if (pdgD0==-421){ //D0
	    fillthis="histSgn_b_";
	    fillthis+=ptbin;
	    ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
            fillthis="histSgn_WeigD0Eff_b_";
            fillthis+=ptbin;
            Double_t effD0 = fCutsTracks->GetTrigWeightB(part->Pt(),fMultEv);
            if(!fUseDeff) effD0=1.; 
            ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar,1./effD0);
  	  } else{ //it was a D0bar
	    fillthis="histRfl_b_";
	    fillthis+=ptbin;
	    ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
            fillthis="histRfl_WeigD0Eff_b_";
            fillthis+=ptbin;
            Double_t effD0 = fCutsTracks->GetTrigWeightB(part->Pt(),fMultEv);
            if(!fUseDeff) effD0=1.; 
            ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar,1./effD0);
	  }
      } else {//background
  	fillthis="histBkg_c_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
        fillthis="histBkg_WeigD0Eff_c_";
        fillthis+=ptbin;
        Double_t effD0 = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
        if(!fUseDeff) effD0=1.; 
        ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar,1./effD0);
      }
    }else{ //on data
      fillthis="histMass_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
      fillthis="histMass_WeigD0Eff_";
      fillthis+=ptbin;
      Double_t effD0 = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
      if(!fUseDeff) effD0=1.; 
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar,1./effD0);
      if(fFillTrees>0) {
	Double_t centFill = 0.;
	if(fCutsD0->GetUseCentrality()) centFill = fCutsD0->GetCentrality(aod);
        ((TH2F*)(listout->FindObject(Form("histMass2D_%d",ptbin))))->Fill(invmassD0bar,centFill);
        ((TH2F*)(listout->FindObject(Form("histMass2D_WeigD0Eff_%d",ptbin))))->Fill(invmassD0bar,centFill,1./effD0);
      }      
    }

  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Correlations: Terminate() \n");

  fOutputMass = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputMass) {     
    printf("ERROR: fOutputMass not available\n");
    return;
  }

  fNentries = dynamic_cast<TH1F*>(GetOutputData(2));
  
  if(!fNentries){
    printf("ERROR: fNEntries not available\n");
    return;
  }

  fCutsD0 = dynamic_cast<AliRDHFCutsD0toKpi*>(GetOutputData(3));
  if(!fCutsD0){
    printf("ERROR: fCuts not available\n");
    return;
  }

  fCounter = dynamic_cast<AliNormalizationCounter*>(GetOutputData(4));    
  if (!fCounter) {
    printf("ERROR: fCounter not available\n");
    return;
  }
  fOutputCorr = dynamic_cast<TList*> (GetOutputData(5));
  if (!fOutputCorr) {     
    printf("ERROR: fOutputCorr not available\n");
    return;
  }
  fOutputStudy = dynamic_cast<TList*> (GetOutputData(6));
  if (!fOutputStudy) {     
    printf("ERROR: fOutputStudy not available\n");
    return;
  }
  fCutsTracks = dynamic_cast<AliHFAssociatedTrackCuts*>(GetOutputData(7));
  if(!fCutsTracks){
    printf("ERROR: fCutsTracks not available\n");
    return;
  }

  return;
}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSED0Correlations::CheckD0Origin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //
  printf("AliAnalysisTaskSED0Correlations::CheckD0Origin() \n");
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;

  while (mother > 0){
    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcMoth){
      pdgGranma = mcMoth->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcMoth->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }
  
  if(isQuarkFound) {
    if(isFromB) return 5;
    else return 4;
  }
  else return 1;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::CreateCorrelationsObjs() {
//

  TString namePlot = "";

  //These for limits in THnSparse (one per bin, same limits). 
  //Vars: DeltaPhi, InvMass, PtTrack, Displacement, DeltaEta --> Last bin for pTassoc is to avoid the overflow!
  Int_t nBinsPhi[5] = {32,150,(int)(2*fPtAssocLimit+1),3,16}; 
  Double_t binMinPhi[5] = {-TMath::Pi()/2.,1.5848,0.,0.,-1.6};  //is the minimum for all the bins
  Double_t binMaxPhi[5] = {3.*TMath::Pi()/2.,2.1848,fPtAssocLimit+0.5,3.,1.6};  //is the maximum for all the bins

  //Vars: DeltaPhi, InvMass, DeltaEta
  Int_t nBinsMix[5] = {32,150,16,(int)(2*fPtAssocLimit+1),2};
  Double_t binMinMix[5] = {-TMath::Pi()/2.,1.5848,-1.6,0.,-0.5};  //is the minimum for all the bins
  Double_t binMaxMix[5] = {3.*TMath::Pi()/2.,2.1848,1.6,fPtAssocLimit+0.5,1.5};  //is the maximum for all the bins

  Int_t nPoolForHistos=1;
  if(!fMergePools) nPoolForHistos= fCutsTracks->GetNZvtxPoolBins()*fCutsTracks->GetNCentPoolBins(); //multeplicity of histos in case of correct pools treatment: sum(SE_i/ME_i)
 
  for(Int_t i=0;i<fNPtBinsCorr;i++) {

    //Modify n of bins with fast speed: in the "for" loop since bins can depend on pT (e.g. mass bin)
    //setting of mass bin is done at the end of the loop!
    if(fSpeed) { //these with fast speed
      if(i>=9) {nBinsPhi[0] = 32; nBinsPhi[1] = 67; nBinsPhi[3] = 1; nBinsPhi[4] = 16;}
      else {nBinsPhi[0] = 32; nBinsPhi[1] = 43; nBinsPhi[3] = 1; nBinsPhi[4] = 16;}
      binMinPhi[0] = -TMath::Pi()/2.; binMinPhi[1] = 1.5848; binMinPhi[3] = 0.; binMinPhi[4] = -1.6;
      binMaxPhi[0] = 3.*TMath::Pi()/2.; binMaxPhi[1] = 2.1848; binMaxPhi[3] = 3.; binMaxPhi[4] = 1.6;
    
      if(i>=9) {nBinsMix[0] = 32; nBinsMix[1] = 67; nBinsMix[2] = 16;}
      else {nBinsMix[0] = 32; nBinsMix[1] = 43; nBinsMix[2] = 16;} 
      binMinMix[0] = -TMath::Pi()/2.; binMinMix[1] = 1.5848; binMinMix[2] = -1.6;
      binMaxMix[0] = 3.*TMath::Pi()/2.; binMaxMix[1] = 2.1848; binMaxMix[2] = 1.6;
    }  	  

    if(!fMixing) {
     if(!fFillTrees) {
      for(Int_t k=0; k<nPoolForHistos; k++) {    	    
    	    
        //THnSparse plots: correlations for various invariant mass (MC and data)
        namePlot="hPhi_K0_Bin";
        namePlot+=i; namePlot+="_p"; namePlot+=k;

        THnSparseF *hPhiK = new THnSparseF(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiK->Sumw2();
        fOutputCorr->Add(hPhiK);

        namePlot="hPhi_Kcharg_Bin";
        namePlot+=i; namePlot+="_p"; namePlot+=k;

        THnSparseF *hPhiH = new THnSparseF(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiH->Sumw2();
        fOutputCorr->Add(hPhiH);

        namePlot="hPhi_Charg_Bin";
        namePlot+=i; namePlot+="_p"; namePlot+=k;

        THnSparseF *hPhiC = new THnSparseF(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiC->Sumw2();
        fOutputCorr->Add(hPhiC);
  
        //histos for c/b origin for D0 (MC only)
        if (fReadMC) {

          //generic origin for tracks
          namePlot="hPhi_K0_From_c_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiK_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiK_c->Sumw2();
          fOutputCorr->Add(hPhiK_c);

          namePlot="hPhi_Kcharg_From_c_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiH_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiH_c->Sumw2();
          fOutputCorr->Add(hPhiH_c);

          namePlot="hPhi_Charg_From_c_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiC_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiC_c->Sumw2();
          fOutputCorr->Add(hPhiC_c);
  
          namePlot="hPhi_K0_From_b_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiK_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiK_b->Sumw2();
          fOutputCorr->Add(hPhiK_b);

          namePlot="hPhi_Kcharg_From_b_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiH_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiH_b->Sumw2();
          fOutputCorr->Add(hPhiH_b);

          namePlot="hPhi_Charg_From_b_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiC_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiC_b->Sumw2();
          fOutputCorr->Add(hPhiC_b);

          //HF-only tracks (c for c->D0, b for b->D0)
          namePlot="hPhi_K0_HF_From_c_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;
  
          THnSparseF *hPhiK_HF_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiK_HF_c->Sumw2();
          fOutputCorr->Add(hPhiK_HF_c);

          namePlot="hPhi_Kcharg_HF_From_c_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiH_HF_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiH_HF_c->Sumw2();
          fOutputCorr->Add(hPhiH_HF_c);

          namePlot="hPhi_Charg_HF_From_c_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiC_HF_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiC_HF_c->Sumw2();
          fOutputCorr->Add(hPhiC_HF_c);

          namePlot="hPhi_K0_HF_From_b_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiK_HF_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiK_HF_b->Sumw2();
          fOutputCorr->Add(hPhiK_HF_b);

          namePlot="hPhi_Kcharg_HF_From_b_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiH_HF_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiH_HF_b->Sumw2();
          fOutputCorr->Add(hPhiH_HF_b);

          namePlot="hPhi_Charg_HF_From_b_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;
     
          THnSparseF *hPhiC_HF_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiC_HF_b->Sumw2();
          fOutputCorr->Add(hPhiC_HF_b);

          namePlot="hPhi_K0_NonHF_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiK_Non = new THnSparseF(namePlot.Data(), "Azimuthal correlation - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiK_Non->Sumw2();
          fOutputCorr->Add(hPhiK_Non);

          namePlot="hPhi_Kcharg_NonHF_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiH_Non = new THnSparseF(namePlot.Data(), "Azimuthal correlation - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiH_Non->Sumw2();
          fOutputCorr->Add(hPhiH_Non);

          namePlot="hPhi_Charg_NonHF_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiC_Non = new THnSparseF(namePlot.Data(), "Azimuthal correlation - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiC_Non->Sumw2();
          fOutputCorr->Add(hPhiC_Non);
        } //end of MC
  
        //modify here the mass axis of THnSparse! 
        if(fSpeed) {
      	  Int_t nBins; Double_t mBin;      
      	  if(i>=9) { //signal range is 1.7488 to 2.0008, plus 1 bin L and R for sidebands
      	    nBins = 67;
      	    mBin = 1.7488;
      	  }
          else { //signal range is 1.7968 to 1.9528, plus 1 bin L and R for sidebands
            nBins = 43;
      	    mBin = 1.7968;
          }
     	
          Double_t varBins[nBins+1];
          varBins[0] = 1.5848;
          varBins[1] = 1.6048;
      	  for(int j = 2; j<nBins-1; j++) {varBins[j]=mBin; mBin+=0.004;}
      	  varBins[nBins-1] = 2.1648;
      	  varBins[nBins] = 2.1848;
        
      	  ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_K0_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
      	  ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Kcharg_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
      	  ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
          if (fReadMC) {
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_K0_From_c_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Kcharg_From_c_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_From_c_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_K0_From_b_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Kcharg_From_b_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_From_b_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_K0_HF_From_c_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Kcharg_HF_From_c_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_HF_From_c_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_K0_HF_From_b_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Kcharg_HF_From_b_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_HF_From_b_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_K0_NonHF_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Kcharg_NonHF_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_NonHF_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
          }
        } //end of fSpeed

      } //end of pool multiplicity      
   
      //Resume the definition of histos
      if(!fSpeed) {
        //leading hadron correlations
        namePlot="hPhi_Lead_Bin";
        namePlot+=i;

        THnSparseF *hCorrLead = new THnSparseF(namePlot.Data(), "Leading particle correlations; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
        hCorrLead->Sumw2();
        fOutputCorr->Add(hCorrLead);

        if (fReadMC) {
          namePlot="hPhi_Lead_From_c_Bin";
          namePlot+=i;

          THnSparseF *hCorrLead_c = new THnSparseF(namePlot.Data(), "Leading particle correlations - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          hCorrLead_c->Sumw2();
          fOutputCorr->Add(hCorrLead_c);
  
          namePlot="hPhi_Lead_From_b_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrLead_b = new THnSparseF(namePlot.Data(), "Leading particle correlations - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          hCorrLead_b->Sumw2();
          fOutputCorr->Add(hCorrLead_b);
  
          namePlot="hPhi_Lead_HF_From_c_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrLead_HF_c = new THnSparseF(namePlot.Data(), "Leading particle correlations HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          hCorrLead_HF_c->Sumw2();
          fOutputCorr->Add(hCorrLead_HF_c);
  
          namePlot="hPhi_Lead_HF_From_b_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrLead_HF_b = new THnSparseF(namePlot.Data(), "Leading particle correlations HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          hCorrLead_HF_b->Sumw2();
          fOutputCorr->Add(hCorrLead_HF_b);

          namePlot="hPhi_Lead_NonHF_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrLead_Non = new THnSparseF(namePlot.Data(), "Leading particle correlations - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          hCorrLead_Non->Sumw2();
          fOutputCorr->Add(hCorrLead_Non);
        }
      
        //pT weighted correlations
        namePlot="hPhi_Weig_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrWeig = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted); #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
        fOutputCorr->Add(hCorrWeig);
  
        if (fReadMC) {
          namePlot="hPhi_Weig_From_c_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrWeig_c = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          fOutputCorr->Add(hCorrWeig_c);
  
          namePlot="hPhi_Weig_From_b_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrWeig_b = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          fOutputCorr->Add(hCorrWeig_b);
  
          namePlot="hPhi_Weig_HF_From_c_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrWeig_HF_c = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          fOutputCorr->Add(hCorrWeig_HF_c);
  
          namePlot="hPhi_Weig_HF_From_b_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrWeig_HF_b = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          fOutputCorr->Add(hCorrWeig_HF_b);

          namePlot="hPhi_Weig_NonHF_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrWeig_Non = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          fOutputCorr->Add(hCorrWeig_Non);
        }
      } //end of fSpeed
     } //end of !fFillTrees

     //pT distribution histos
     namePlot = "hist_Pt_Charg_Bin"; namePlot+=i;
     TH1F *hPtC = new TH1F(namePlot.Data(), "Charged track pT (in D0 evs); p_{T} (GeV/c)",240,0.,12.);
     hPtC->SetMinimum(0);
     fOutputStudy->Add(hPtC);

     namePlot = "hist_Pt_Kcharg_Bin"; namePlot+=i;
     TH1F *hPtH = new TH1F(namePlot.Data(), "Hadrons pT (in D0 evs); p_{T} (GeV/c)",240,0.,12.);
     hPtH->SetMinimum(0);
     fOutputStudy->Add(hPtH);

     namePlot = "hist_Pt_K0_Bin"; namePlot+=i;
     TH1F *hPtK = new TH1F(namePlot.Data(), "Kaons pT (in D0 evs); p_{T} (GeV/c)",240,0.,12.);
     hPtK->SetMinimum(0);
     fOutputStudy->Add(hPtK);

     //Events multiplicity
     namePlot = "hMultiplEvt_Bin"; namePlot+=i;
     TH1F *hMultEv = new TH1F(namePlot.Data(), "Event multiplicity",1500,0.,6000.);
     hMultEv->SetMinimum(0);
     fOutputStudy->Add(hMultEv);

    } //end of !fMixing

    if(fMixing && !fFillTrees) {
      for(Int_t k=0; k<nPoolForHistos; k++) {     	    
        //THnSparse plots for event mixing!
        namePlot="hPhi_K0_Bin";
        namePlot+=i; namePlot+="_p"; namePlot+=k; namePlot+="_EvMix";

        THnSparseF *hPhiK_EvMix = new THnSparseF(namePlot.Data(), "Az. corr. EvMix; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
        hPhiK_EvMix->Sumw2();
        fOutputCorr->Add(hPhiK_EvMix);

        namePlot="hPhi_Kcharg_Bin";
        namePlot+=i; namePlot+="_p"; namePlot+=k; namePlot+="_EvMix";
  
        THnSparseF *hPhiH_EvMix = new THnSparseF(namePlot.Data(), "Az. corr. EvMix; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
        hPhiH_EvMix->Sumw2();
        fOutputCorr->Add(hPhiH_EvMix);

        namePlot="hPhi_Charg_Bin";
        namePlot+=i; namePlot+="_p"; namePlot+=k; namePlot+="_EvMix";

        THnSparseF *hPhiC_EvMix = new THnSparseF(namePlot.Data(), "Az. corr. EvMix; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
        hPhiC_EvMix->Sumw2();
        fOutputCorr->Add(hPhiC_EvMix);  

        //modify here the mass axis of THnSparse! 
        if(fSpeed) {
      	  Int_t nBins; Double_t mBin;      
      	  if(i>=9) { //signal range is 1.7488 to 2.0008, plus 1 bin L and R for sidebands
      	    nBins = 67;
      	    mBin = 1.7488;
      	  }
          else { //signal range is 1.7968 to 1.9528, plus 1 bin L and R for sidebands
            nBins = 43;
      	    mBin = 1.7968;
          }
     	
          Double_t varBins[nBins+1];
          varBins[0] = 1.5848;
          varBins[1] = 1.6048;
      	  for(int j = 2; j<nBins-1; j++) {varBins[j]=mBin; mBin+=0.004;}
      	  varBins[nBins-1] = 2.1648;
      	  varBins[nBins] = 2.1848;
        
      	  ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_K0_Bin%d_p%d_EvMix",i,k)))->GetAxis(1)->Set(nBins, varBins);
      	  ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Kcharg_Bin%d_p%d_EvMix",i,k)))->GetAxis(1)->Set(nBins, varBins);
      	  ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d_EvMix",i,k)))->GetAxis(1)->Set(nBins, varBins);
      	  
        } //end of fSpeed   
      
      } //end of Mult pools
    } //end of Mix
 
    //both for SE and for ME
    //D* feeddown pions rejection histos
    namePlot = "hDstarPionsVsDmass_Bin"; namePlot+=i;
    TH2F *hDstarPions = new TH2F(namePlot.Data(), "Tracks rejected for D* inv.mass cut vs D inv mass; # Tracks",2,0.,2.,150,1.5848,2.1848);
    hDstarPions->GetXaxis()->SetBinLabel(1,"Not rejected");
    hDstarPions->GetXaxis()->SetBinLabel(2,"Rejected");
    hDstarPions->SetMinimum(0);
    fOutputStudy->Add(hDstarPions); 

    namePlot = "hDstarPionsVsdeltaPhi_Bin"; namePlot+=i;
    TH2F *hDstarPions2 = new TH2F(namePlot.Data(), "Tracks rejected for D* inv.mass cut vs deltaPhi; # Tracks",2,0.,2.,64,-TMath::Pi()/2.,3.*TMath::Pi()/2.);
    hDstarPions2->GetXaxis()->SetBinLabel(1,"Not rejected");
    hDstarPions2->GetXaxis()->SetBinLabel(2,"Rejected");
    hDstarPions2->SetMinimum(0);
    fOutputStudy->Add(hDstarPions2); 

    if(!fFillTrees) {
      //ME filling control plots
      namePlot="hEvtsPerPool_"; namePlot+=i;
      TH1F *hEvPerPool = new TH1F(namePlot.Data(), "Events With selD0 in ME pools",nPoolForHistos,0.,nPoolForHistos);
      hEvPerPool->SetMinimum(0);
      fOutputStudy->Add(hEvPerPool);
    }

  } //end of bin loop

  //out of bin loop
  TH1F *hCountC = new TH1F("hist_Count_Charg", "Charged track counter; # Tracks",6000,0.,6000.);
  hCountC->SetMinimum(0);
  fOutputStudy->Add(hCountC);

  TH1F *hCountH = new TH1F("hist_Count_Kcharg", "Hadrons counter; # Tracks",100,0.,100.);
  hCountH->SetMinimum(0);
  fOutputStudy->Add(hCountH);

  TH1F *hCountK = new TH1F("hist_Count_K0", "Kaons counter; # Tracks",100,0.,100.);
  hCountK->SetMinimum(0);
  fOutputStudy->Add(hCountK);

  TH1F *hZvtx = new TH1F("hZvtx", "z of Primary vtx (for events with selected D); # Events",48,-12.,12.);
  hZvtx->SetMinimum(0);
  fOutputStudy->Add(hZvtx);
  
  TH1F *hZvtxEvts = new TH1F("hZvtxEvts", "z of Primary vtx (for selected events); # Events",120,-30.,30.);
  hZvtxEvts->SetMinimum(0);
  fOutputStudy->Add(hZvtxEvts);
  
  TH1F *hCentralEvts = new TH1F("hCentralEvts","Centrality of events (for selected events); # Events",102,-1.,101.);
  hCentralEvts->SetMinimum(0);
  fOutputStudy->Add(hCentralEvts);  

  if (fReadMC) {
    TH1D *hEventTypeMC = new TH1D("EventTypeMC","EventTypeMC",100,-0.5,99.5);
    fOutputStudy->Add(hEventTypeMC); 
  }

  if (fFillGlobal) { //all-events plots
    //pt distributions
    TH1F *hPtCAll = new TH1F("hist_Pt_Charg_AllEv", "Charged track pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtCAll->SetMinimum(0);
    fOutputStudy->Add(hPtCAll);

    TH1F *hPtHAll = new TH1F("hist_Pt_Kcharg_AllEv", "Kaons pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtHAll->SetMinimum(0);
    fOutputStudy->Add(hPtHAll);

    TH1F *hPtKAll = new TH1F("hist_Pt_K0_AllEv", "K0 pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtKAll->SetMinimum(0);
    fOutputStudy->Add(hPtKAll);

    //K0 Invariant Mass plots
    TH2F *hK0MassInv = new TH2F("hK0MassInv", "K0 invariant mass; Invariant mass (MeV/c^{2}); pT (GeV/c)",200,0.4,0.6,100,0.,10.);
    hK0MassInv->SetMinimum(0);
    fOutputStudy->Add(hK0MassInv);
  }

  if(!fMixing) {
    //phi distributions
    TH1F *hPhiDistCAll = new TH1F("hist_PhiDistr_Charg", "Charged track phi distr. (All); p_{T} (GeV/c)",64,0,6.283);
    hPhiDistCAll->SetMinimum(0);
    fOutputStudy->Add(hPhiDistCAll);

    TH1F *hPhiDistHAll = new TH1F("hist_PhiDistr_Kcharg", "Kaons phi distr. (All); p_{T} (GeV/c)",64,0,6.283);
    hPhiDistHAll->SetMinimum(0);
    fOutputStudy->Add(hPhiDistHAll);

    TH1F *hPhiDistKAll = new TH1F("hist_PhiDistr_K0", "K0 phi distr. (All); p_{T} (GeV/c)",64,0,6.283);
    hPhiDistKAll->SetMinimum(0);
    fOutputStudy->Add(hPhiDistKAll);

    TH1F *hPhiDistDAll = new TH1F("hist_PhiDistr_D0", "D^{0} phi distr. (All); p_{T} (GeV/c)",64,0,6.283);
    hPhiDistDAll->SetMinimum(0);
    fOutputStudy->Add(hPhiDistDAll);

    //phi distributions
    TH1F *hEtaDistCAll = new TH1F("hist_EtaDistr_Charg", "Charged track eta distr. (All); p_{T} (GeV/c)",40,-1,1);
    hEtaDistCAll->SetMinimum(0);
    fOutputStudy->Add(hEtaDistCAll);

    TH1F *hEtaDistHAll = new TH1F("hist_EtaDistr_Kcharg", "Kaons eta distr. (All); p_{T} (GeV/c)",40,-1,1);
    hEtaDistHAll->SetMinimum(0);
    fOutputStudy->Add(hEtaDistHAll);

    TH1F *hEtaDistKAll = new TH1F("hist_EtaDistr_K0", "K0 eta distr. (All); p_{T} (GeV/c)",40,-1,1);
    hEtaDistKAll->SetMinimum(0);
    fOutputStudy->Add(hEtaDistKAll);

    TH1F *hEtaDistDAll = new TH1F("hist_EtaDistr_D0", "D^{0} eta distr. (All); p_{T} (GeV/c)",40,-1,1);
    hEtaDistDAll->SetMinimum(0);
    fOutputStudy->Add(hEtaDistDAll);
    
    //phivsEta
    TH2F *hPhiVsEtaDistCAll = new TH2F("hist_PhiVsEtaDistr_Charg", "Phi vs Eta distribution - Charged tracks",64,0,6.283,40,-1,1);
    hPhiVsEtaDistCAll->SetMinimum(0);
    fOutputStudy->Add(hPhiVsEtaDistCAll);     
    
    TH2F *hPhiVsEtaDistHAll = new TH2F("hist_PhiVsEtaDistr_Kcharg", "Phi vs Eta distribution - Charged Kaons",64,0,6.283,40,-1,1);
    hPhiVsEtaDistHAll->SetMinimum(0);
    fOutputStudy->Add(hPhiVsEtaDistHAll);  
    
    TH2F *hPhiVsEtaDistKAll = new TH2F("hist_PhiVsEtaDistr_K0", "Phi vs Eta distribution - K^{0}",64,0,6.283,40,-1,1);
    hPhiVsEtaDistKAll->SetMinimum(0);
    fOutputStudy->Add(hPhiVsEtaDistKAll);  
    
    TH2F *hPhiVsEtaDistDAll = new TH2F("hist_PhiVsEtaDistr_D0", "Phi vs Eta distribution - D^{0}",64,0,6.283,40,-1,1);
    hPhiVsEtaDistDAll->SetMinimum(0);
    fOutputStudy->Add(hPhiVsEtaDistDAll);  
    }

  //for MC analysis only
  for(Int_t i=0;i<fNPtBinsCorr;i++) {

    if (fReadMC && !fMixing) {

      //displacement histos
      namePlot="histDispl_K0_Bin"; namePlot+=i;
      TH1F *hDisplK = new TH1F(namePlot.Data(), "Kaons Displacement; DCA",150,0.,0.15);
      hDisplK->SetMinimum(0);
      fOutputStudy->Add(hDisplK);
  
      namePlot="histDispl_K0_HF_Bin";  namePlot+=i;
      TH1F *hDisplK_HF = new TH1F(namePlot.Data(), "Kaons Displacement (from HF decay only); DCA",150,0.,0.15);
      hDisplK_HF->SetMinimum(0);
      fOutputStudy->Add(hDisplK_HF);

      namePlot="histDispl_Kcharg_Bin"; namePlot+=i;
      TH1F *hDisplHadr = new TH1F(namePlot.Data(), "Hadrons Displacement; DCA",150,0.,0.15);
      hDisplHadr->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr);
  
      namePlot="histDispl_Kcharg_HF_Bin";  namePlot+=i;
      TH1F *hDisplHadr_HF = new TH1F(namePlot.Data(), "Hadrons Displacement (from HF decay only); DCA",150,0.,0.15);
      hDisplHadr_HF->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_HF);

      namePlot="histDispl_Charg_Bin"; namePlot+=i;
      TH1F *hDisplCharg = new TH1F(namePlot.Data(), "Charged tracks Displacement; DCA",150,0.,0.15);
      hDisplCharg->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg);
  
      namePlot="histDispl_Charg_HF_Bin";  namePlot+=i;
      TH1F *hDisplCharg_HF = new TH1F(namePlot.Data(), "Charged tracks Displacement (from HF decay only); DCA",150,0.,0.15);
      hDisplCharg_HF->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_HF);

      namePlot="histDispl_K0_From_c_Bin"; namePlot+=i;
      TH1F *hDisplK_c = new TH1F(namePlot.Data(), "Kaons Displacement - c origin; DCA",150,0.,0.15);
      hDisplK_c->SetMinimum(0);
      fOutputStudy->Add(hDisplK_c);
  
      namePlot="histDispl_K0_HF_From_c_Bin";  namePlot+=i;
      TH1F *hDisplK_HF_c = new TH1F(namePlot.Data(), "Kaons Displacement (from HF decay only) - c origin; DCA",150,0.,0.15);
      hDisplK_HF_c->SetMinimum(0);
      fOutputStudy->Add(hDisplK_HF_c);

      namePlot="histDispl_Kcharg_From_c_Bin"; namePlot+=i;
      TH1F *hDisplHadr_c = new TH1F(namePlot.Data(), "Hadrons Displacement - c origin; DCA",150,0.,0.15);
      hDisplHadr_c->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_c);
  
      namePlot="histDispl_Kcharg_HF_From_c_Bin";  namePlot+=i;
      TH1F *hDisplHadr_HF_c = new TH1F(namePlot.Data(), "Hadrons Displacement (from HF decay only) - c origin; DCA",150,0.,0.15);
      hDisplHadr_HF_c->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_HF_c);

      namePlot="histDispl_Charg_From_c_Bin"; namePlot+=i;
      TH1F *hDisplCharg_c = new TH1F(namePlot.Data(), "Charged tracks Displacement - c origin; DCA",150,0.,0.15);
      hDisplCharg_c->Sumw2();
      hDisplCharg_c->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_c);
  
      namePlot="histDispl_Charg_HF_From_c_Bin";  namePlot+=i;
      TH1F *hDisplCharg_HF_c = new TH1F(namePlot.Data(), "Charged tracks Displacement (from HF decay only) - c origin; DCA",150,0.,0.15);
      hDisplCharg_HF_c->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_HF_c);

      namePlot="histDispl_K0_From_b_Bin"; namePlot+=i;
      TH1F *hDisplK_b = new TH1F(namePlot.Data(), "Kaons Displacement - b origin; DCA",150,0.,0.15);
      hDisplK_b->SetMinimum(0);
      fOutputStudy->Add(hDisplK_b);
  
      namePlot="histDispl_K0_HF_From_b_Bin";  namePlot+=i;
      TH1F *hDisplK_HF_b = new TH1F(namePlot.Data(), "Kaons Displacement (from HF decay only) - b origin; DCA",150,0.,0.15);
      hDisplK_HF_b->SetMinimum(0);
      fOutputStudy->Add(hDisplK_HF_b);

      namePlot="histDispl_Kcharg_From_b_Bin"; namePlot+=i;
      TH1F *hDisplHadr_b = new TH1F(namePlot.Data(), "Hadrons Displacement - b origin; DCA",150,0.,0.15);
      hDisplHadr_b->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_b);

      namePlot="histDispl_Kcharg_HF_From_b_Bin";  namePlot+=i;
      TH1F *hDisplHadr_HF_b = new TH1F(namePlot.Data(), "Hadrons Displacement (from HF decay only) - b origin; DCA",150,0.,0.15);
      hDisplHadr_HF_b->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_HF_b);

      namePlot="histDispl_Charg_From_b_Bin"; namePlot+=i;
      TH1F *hDisplCharg_b = new TH1F(namePlot.Data(), "Charged tracks Displacement - b origin; DCA",150,0.,0.15);
      hDisplCharg_b->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_b);
  
      namePlot="histDispl_Charg_HF_From_b_Bin";  namePlot+=i;
      TH1F *hDisplCharg_HF_b = new TH1F(namePlot.Data(), "Charged tracks Displacement (from HF decay only) - b origin; DCA",150,0.,0.15);
      hDisplCharg_HF_b->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_HF_b);

      //origin of tracks histos
      namePlot="histOrig_Charg_Bin";  namePlot+=i;
      TH1F *hOrigin_Charm = new TH1F(namePlot.Data(), "Origin of charged tracks",9,0.,9.);
      hOrigin_Charm->SetMinimum(0);
      hOrigin_Charm->GetXaxis()->SetBinLabel(1,"Not HF");
      hOrigin_Charm->GetXaxis()->SetBinLabel(2,"D->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(3,"D->X->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(4,"c hadr.");
      hOrigin_Charm->GetXaxis()->SetBinLabel(5,"B->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(6,"B->X-># (X!=D)");
      hOrigin_Charm->GetXaxis()->SetBinLabel(7,"B->D->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(8,"B->D->X->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(9,"b hadr.");
      fOutputStudy->Add(hOrigin_Charm);

      namePlot="histOrig_Kcharg_Bin";  namePlot+=i;
      TH1F *hOrigin_Kcharg = new TH1F(namePlot.Data(), "Origin of hadrons",9,0.,9.);
      hOrigin_Kcharg->SetMinimum(0);
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(1,"Not HF");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(2,"D->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(3,"D->X->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(4,"c hadr.");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(5,"B->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(6,"B->X-># (X!=D)");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(7,"B->D->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(8,"B->D->X->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(9,"b hadr.");
      fOutputStudy->Add(hOrigin_Kcharg);

      namePlot="histOrig_K0_Bin";  namePlot+=i;
      TH1F *hOrigin_K = new TH1F(namePlot.Data(), "Origin of kaons",9,0.,9.);
      hOrigin_K->SetMinimum(0);
      hOrigin_K->GetXaxis()->SetBinLabel(1,"Not HF");
      hOrigin_K->GetXaxis()->SetBinLabel(2,"D->#");
      hOrigin_K->GetXaxis()->SetBinLabel(3,"D->X->#");
      hOrigin_K->GetXaxis()->SetBinLabel(4,"c hadr.");
      hOrigin_K->GetXaxis()->SetBinLabel(5,"B->#");
      hOrigin_K->GetXaxis()->SetBinLabel(6,"B->X-># (X!=D)");
      hOrigin_K->GetXaxis()->SetBinLabel(7,"B->D->#");
      hOrigin_K->GetXaxis()->SetBinLabel(8,"B->D->X->#");
      hOrigin_K->GetXaxis()->SetBinLabel(9,"b hadr.");
      fOutputStudy->Add(hOrigin_K);
    }

    if (fReadMC) {
      //origin of D0 histos
      namePlot="histOrig_D0_Bin";  namePlot+=i;
      TH1F *hOrigin_D0 = new TH1F(namePlot.Data(), "Origin of D0",2,0.,2.);
      hOrigin_D0->SetMinimum(0);
      hOrigin_D0->GetXaxis()->SetBinLabel(1,"From c");
      hOrigin_D0->GetXaxis()->SetBinLabel(2,"From b");
      fOutputStudy->Add(hOrigin_D0);

      //primary tracks (Kine & Reco)
      namePlot="hPhysPrim_Bin";  namePlot+=i;
      TH1F *hPhysPrim = new TH1F(namePlot.Data(), "Origin of hadrons",2,0.,2.);
      hPhysPrim->SetMinimum(0);
      hPhysPrim->GetXaxis()->SetBinLabel(1,"OK");
      hPhysPrim->GetXaxis()->SetBinLabel(2,"NO");
      fOutputStudy->Add(hPhysPrim);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::CalculateCorrelations(AliAODRecoDecayHF2Prong* d, Int_t labD0, TClonesArray* mcArray) {
//
// Method for correlations D0-hadrons study
//
  Int_t N_Charg = 0, N_KCharg = 0, N_Kaons = 0;
  Double_t mD0, mD0bar;
  Int_t origD0 = 0, PDGD0 = 0, ptbin = 0;
  d->InvMassD0(mD0,mD0bar);
  Double_t mInv[2] = {mD0, mD0bar};
  ptbin = PtBinCorr(d->Pt());

  if(ptbin < 0) return;

  //Fill of D0 phi distribution
  if (!fMixing) ((TH1F*)fOutputStudy->FindObject("hist_PhiDistr_D0"))->Fill(d->Phi());  
  if (!fMixing) ((TH1F*)fOutputStudy->FindObject("hist_EtaDistr_D0"))->Fill(d->Eta()); 
  if (!fMixing) ((TH2F*)fOutputStudy->FindObject("hist_PhiVsEtaDistr_D0"))->Fill(d->Phi(),d->Eta()); 

  //Origin of D0
  TString orig="";
  if(fReadMC) {
    origD0=CheckD0Origin(mcArray,(AliAODMCParticle*)mcArray->At(labD0));
    PDGD0 = ((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode();
    switch (CheckD0Origin(mcArray,(AliAODMCParticle*)mcArray->At(labD0))) {
      case 4:
        orig = "_From_c";
        ((TH1F*)fOutputStudy->FindObject(Form("histOrig_D0_Bin%d",ptbin)))->Fill(0.);
        break;
      case 5:
        orig = "_From_b";
        ((TH1F*)fOutputStudy->FindObject(Form("histOrig_D0_Bin%d",ptbin)))->Fill(1.);
        break;
      default:
        return;
    }
  }

  Double_t highPt = 0; Double_t lead[4] = {0,0,0,1};  //infos for leading particle (pt,deltaphi)

  //loop over the tracks in the pool 
  Bool_t execPoolTr = fCorrelatorTr->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
  Bool_t execPoolKc = fCorrelatorKc->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
  Bool_t execPoolK0 = fCorrelatorK0->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
		
  Int_t NofEventsinPool = 1;
  if(fMixing) {
    NofEventsinPool = fCorrelatorTr->GetNofEventsInPool(); 
    if(!execPoolTr) {
      AliInfo("Mixed event analysis: track pool is not ready");
      NofEventsinPool = 0;
    }
  }

  //Charged tracks
  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++) {// loop on events in the pool; if it is SE analysis, stops at one (index not needed there)
    Bool_t analyzetracksTr = fCorrelatorTr->ProcessAssociatedTracks(jMix);// process all the tracks in the aodEvent, by applying the selection cuts
    if(!analyzetracksTr) {
      AliInfo("AliHFCorrelator::Cannot process the track array");
      continue;
    }
	
    for(Int_t iTrack = 0; iTrack<fCorrelatorTr->GetNofTracks(); iTrack++){ // looping on track candidates

      Bool_t runcorrelation = fCorrelatorTr->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* track = fCorrelatorTr->GetAssociatedParticle();

      if(!fMixing) {
        Int_t idDaughs[2] = {((AliVTrack*)d->GetDaughter(0))->GetID(),((AliVTrack*)d->GetDaughter(1))->GetID()}; //IDs of daughters to be skipped
        if(track->GetID() == idDaughs[0] || track->GetID() == idDaughs[1]) continue; //discards daughters of candidate
      }
      if(track->Pt() < fPtThreshLow.at(ptbin) || track->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K

      if(fReadMC) {
        AliAODMCParticle* trkKine = (AliAODMCParticle*)mcArray->At(track->GetLabel());
        if (!trkKine) continue;
        if (!trkKine->IsPhysicalPrimary()) {
 	  ((TH1F*)fOutputStudy->FindObject(Form("hPhysPrim_Bin%d",ptbin)))->Fill(1.);  
  	  continue; //reject the Reco track if correspondent Kine track is not primary
        } else ((TH1F*)fOutputStudy->FindObject(Form("hPhysPrim_Bin%d",ptbin)))->Fill(0.);
      }

      Double_t effTr = track->GetWeight(); //extract track efficiency
      Double_t effD0 = 1.;
      if(fReadMC) {
        if(origD0==4) effD0 = fCutsTracks->GetTrigWeight(d->Pt(),fMultEv);
        if(origD0==5) effD0 = fCutsTracks->GetTrigWeightB(d->Pt(),fMultEv);
      } else effD0 = fCutsTracks->GetTrigWeight(d->Pt(),fMultEv);
      if(!fUseDeff) effD0=1.; 
      if(!fUseTrackeff) effTr=1.; 
      Double_t eff = effTr*effD0;

      if(!fMixing) {
        if(fSoftPiCut && !track->CheckSoftPi()) { //removal of soft pions
          if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mD0);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mD0bar);
          ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsdeltaPhi_Bin%d",ptbin)))->Fill(1.,fCorrelatorTr->GetDeltaPhi());
    	  continue; //in SE events, just reject the soft pion
        } else { //not a soft pion
          if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mD0);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mD0bar);
          ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsdeltaPhi_Bin%d",ptbin)))->Fill(0.,fCorrelatorTr->GetDeltaPhi());
        }
      }
      if(fMixing) { 
        if(fSoftPiCut && !fCutsTracks->InvMassDstarRejection(d,track,fIsSelectedCandidate)) { //removal of soft pions
          if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mD0);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mD0bar);
          ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsdeltaPhi_Bin%d",ptbin)))->Fill(1.,fCorrelatorTr->GetDeltaPhi());
          if(fMixing) FillSparsePlots(mcArray,mInv,origD0,PDGD0,track,ptbin,kTrack,1,1./eff); //in ME events, fill the THnSparse under the softpi hypothesis
    	  continue; 
        } else { //not a soft pion
          if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mD0);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mD0bar);
          ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsdeltaPhi_Bin%d",ptbin)))->Fill(0.,fCorrelatorTr->GetDeltaPhi());
        }
      }
 
      FillSparsePlots(mcArray,mInv,origD0,PDGD0,track,ptbin,kTrack,0,1./eff); //fills for charged tracks

      if(!fMixing) N_Charg++;

      //retrieving leading info...
      if(track->Pt() > highPt) {
        if(fReadMC && track->GetLabel()<1) continue;
        if(fReadMC && !(AliAODMCParticle*)mcArray->At(track->GetLabel())) continue;
        lead[0] = fCorrelatorTr->GetDeltaPhi();
        lead[1] = fCorrelatorTr->GetDeltaEta();
        lead[2] = fReadMC ? CheckTrackOrigin(mcArray,(AliAODMCParticle*)mcArray->At(track->GetLabel())) : 0;
        if(fReadMC) {
  	  if(origD0==4) lead[3] = 1./(track->GetWeight()*fCutsTracks->GetTrigWeight(d->Pt(),fMultEv)); //weight is 1./efficiency
	  if(origD0==5) lead[3] = 1./(track->GetWeight()*fCutsTracks->GetTrigWeightB(d->Pt(),fMultEv)); //weight is 1./efficiency
	} else lead[3] = 1./(track->GetWeight()*fCutsTracks->GetTrigWeight(d->Pt(),fMultEv));
        highPt = track->Pt();
      }

    } // end of tracks loop
  } //end of event loop for fCorrelatorTr

 if(fKaonCorr) { //loops for Kcharg and K0

  if(fMixing) {
    NofEventsinPool = fCorrelatorKc->GetNofEventsInPool(); 
    if(!execPoolKc) {
      AliInfo("Mixed event analysis: K+/- pool is not ready");
      NofEventsinPool = 0;
    }
  }

  //Charged Kaons loop
  for (Int_t jMix = 0; jMix < NofEventsinPool; jMix++) {// loop on events in the pool; if it is SE analysis, stops at one (index not needed there)
    Bool_t analyzetracksKc = fCorrelatorKc->ProcessAssociatedTracks(jMix);
    if(!analyzetracksKc) {
      AliInfo("AliHFCorrelator::Cannot process the K+/- array");
      continue;
    }  

    for(Int_t iTrack = 0; iTrack<fCorrelatorKc->GetNofTracks(); iTrack++){ // looping on charged kaons candidates

      Bool_t runcorrelation = fCorrelatorKc->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* kCharg = fCorrelatorKc->GetAssociatedParticle();

      if(!fMixing) {  
        Int_t idDaughs[2] = {((AliVTrack*)d->GetDaughter(0))->GetID(),((AliVTrack*)d->GetDaughter(1))->GetID()}; //IDs of daughters to be skipped
        if(kCharg->GetID() == idDaughs[0] || kCharg->GetID() == idDaughs[1]) continue; //discards daughters of candidate
      }
      if(kCharg->Pt() < fPtThreshLow.at(ptbin) || kCharg->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K
  
      if(!fMixing) {
        if(fSoftPiCut && !kCharg->CheckSoftPi()) { //removal of soft pions
          if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mD0);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mD0bar);
          ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsdeltaPhi_Bin%d",ptbin)))->Fill(1.,fCorrelatorKc->GetDeltaPhi());
    	  continue; //in SE events, just reject the soft pion
        } else { //not a soft pion
          if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mD0);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mD0bar);
          ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsdeltaPhi_Bin%d",ptbin)))->Fill(0.,fCorrelatorKc->GetDeltaPhi());
        }
      }
      if(fMixing) { 
        if(fSoftPiCut && !fCutsTracks->InvMassDstarRejection(d,kCharg,fIsSelectedCandidate)) { //removal of soft pions
          if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mD0);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mD0bar);
          ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsdeltaPhi_Bin%d",ptbin)))->Fill(1.,fCorrelatorKc->GetDeltaPhi());
          if(fMixing) FillSparsePlots(mcArray,mInv,origD0,PDGD0,kCharg,ptbin,kKCharg,1); //fills for charged tracks
    	  continue; 
        } else { //not a soft pion
          if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mD0);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mD0bar);
          ((TH2F*)fOutputStudy->FindObject(Form("hDstarPionsVsdeltaPhi_Bin%d",ptbin)))->Fill(0.,fCorrelatorKc->GetDeltaPhi());
        }
      }
       
      
      
      
      
      

      FillSparsePlots(mcArray,mInv,origD0,PDGD0,kCharg,ptbin,kKCharg,0); //fills for charged tracks

      if(!fMixing) N_KCharg++;

    } // end of charged kaons loop
  } //end of event loop for fCorrelatorKc

  if(fMixing) {
    NofEventsinPool = fCorrelatorK0->GetNofEventsInPool(); 
    if(!execPoolK0) {
      AliInfo("Mixed event analysis: K0 pool is not ready");
      NofEventsinPool = 0;
    }
  }

  //K0 loop
  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++) {// loop on events in the pool; if it is SE analysis, stops at one (index not needed there)
    Bool_t analyzetracksK0 = fCorrelatorK0->ProcessAssociatedTracks(jMix);
    if(!analyzetracksK0) {
      AliInfo("AliHFCorrelator::Cannot process the K0 array");
      continue;
    }  

    for(Int_t iTrack = 0; iTrack<fCorrelatorK0->GetNofTracks(); iTrack++){ // looping on k0 candidates

      Bool_t runcorrelation = fCorrelatorK0->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* k0 = fCorrelatorK0->GetAssociatedParticle();

      if(k0->Pt() < fPtThreshLow.at(ptbin) || k0->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K
  
      FillSparsePlots(mcArray,mInv,origD0,PDGD0,k0,ptbin,kK0,0); //fills for charged tracks

      if(!fMixing) N_Kaons++;

    } // end of charged kaons loop
  } //end of event loop for fCorrelatorK0

 } //end of 'if(fKaonCorr)'

  Double_t fillSpLeadD0[4] = {lead[0],mD0,lead[1],0.4}; //dummy value for threshold of leading!
  Double_t fillSpLeadD0bar[4] = {lead[0],mD0bar,lead[1],0.4};

  //leading track correlations fill
  if(!fMixing && !fSpeed) {
    if(fReadMC) {
      if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==421 && (fIsSelectedCandidate==1||fIsSelectedCandidate==3)) { //D0
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadD0,lead[3]); //c and b D0
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadD0,lead[3]); //c or b D0
        if(origD0==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadD0,lead[3]);  
        if(origD0==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadD0,lead[3]);  
        if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadD0,lead[3]);  //non HF  
      }
      if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==-421 && fIsSelectedCandidate>1) { //D0bar
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadD0bar,lead[3]);
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadD0bar,lead[3]); //c or b D0
        if(origD0==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadD0bar,lead[3]);  
        if(origD0==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadD0bar,lead[3]); 
        if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadD0bar,lead[3]);  //non HF  
      }
    } else {
        if(fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadD0,lead[3]); 
        if(fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadD0bar,lead[3]);
    }
  }
    //Fill of count histograms
  if (!fAlreadyFilled && !fMixing) { 
    ((TH1F*)fOutputStudy->FindObject("hist_Count_Charg"))->Fill(N_Charg);
    ((TH1F*)fOutputStudy->FindObject("hist_Count_Kcharg"))->Fill(N_KCharg);
    ((TH1F*)fOutputStudy->FindObject("hist_Count_K0"))->Fill(N_Kaons);
  }


  fAlreadyFilled=kTRUE; //at least a D0 analyzed in the event; distribution plots already filled

}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::CalculateCorrelationsMCKine(AliAODMCParticle* d, TClonesArray* mcArray) {
//
// Method for correlations D0-hadrons study
//
  Int_t N_Charg = 0, N_KCharg = 0, N_Kaons = 0;
  Double_t mD0 = 1.864, mD0bar = 1.864;
  Double_t mInv[2] = {mD0, mD0bar};
  Int_t origD0 = 0, PDGD0 = 0;
  Int_t ptbin = PtBinCorr(d->Pt());

  if(ptbin < 0) return;

  //Fill of D0 phi distribution
  if (!fMixing) ((TH1F*)fOutputStudy->FindObject("hist_PhiDistr_D0"))->Fill(d->Phi()); 
  if (!fMixing) ((TH1F*)fOutputStudy->FindObject("hist_EtaDistr_D0"))->Fill(d->Phi()); 
  if (!fMixing) ((TH2F*)fOutputStudy->FindObject("hist_PhiVsEtaDistr_D0"))->Fill(d->Phi(),d->Eta()); 
  
  //Origin of D0
  TString orig="";
  origD0=CheckD0Origin(mcArray,d);
  PDGD0 = d->GetPdgCode();
  switch (CheckD0Origin(mcArray,d)) {
    case 4:
      orig = "_From_c";
      ((TH1F*)fOutputStudy->FindObject(Form("histOrig_D0_Bin%d",ptbin)))->Fill(0.);
      break;
    case 5:
      orig = "_From_b";
      ((TH1F*)fOutputStudy->FindObject(Form("histOrig_D0_Bin%d",ptbin)))->Fill(1.);
      break;
    default:
      return;
  }

  Double_t highPt = 0; Double_t lead[3] = {0,0,0};  //infos for leading particle (pt,deltaphi)

  //loop over the tracks in the pool 
  Bool_t execPoolTr = fCorrelatorTr->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
  Bool_t execPoolKc = fCorrelatorKc->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
  Bool_t execPoolK0 = fCorrelatorK0->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
		
  Int_t NofEventsinPool = 1;
  if(fMixing) {
    NofEventsinPool = fCorrelatorTr->GetNofEventsInPool(); 
    if(!execPoolTr) {
      AliInfo("Mixed event analysis: track pool is not ready");
      NofEventsinPool = 0;
    }
  }

  //Charged tracks
  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++) {// loop on events in the pool; if it is SE analysis, stops at one (index not needed there)

    Bool_t analyzetracksTr = fCorrelatorTr->ProcessAssociatedTracks(jMix);// process all the tracks in the aodEvent, by applying the selection cuts
    if(!analyzetracksTr) {
      AliInfo("AliHFCorrelator::Cannot process the track array");
      continue;
    }
	
    for(Int_t iTrack = 0; iTrack<fCorrelatorTr->GetNofTracks(); iTrack++){ // looping on track candidates

      Bool_t runcorrelation = fCorrelatorTr->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* track = fCorrelatorTr->GetAssociatedParticle();
      if(track->GetLabel()<0) continue;
      if(track->Pt() < fPtThreshLow.at(ptbin) || track->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K
      if(track->Pt() < 0.3 || TMath::Abs(track->Eta())>0.8) continue; //discard tracks outside barrel (since it's kinematic MC and produces tracks all over rapidity region
      if(!fMixing) N_Charg++;

      AliAODMCParticle *trkMC = (AliAODMCParticle*)mcArray->At(track->GetLabel());
      if(!trkMC) continue;

      if (!trkMC->IsPhysicalPrimary()) {  //reject material budget, or other fake tracks
 	((TH1F*)fOutputStudy->FindObject(Form("hPhysPrim_Bin%d",ptbin)))->Fill(1.);  
  	continue;
      } else ((TH1F*)fOutputStudy->FindObject(Form("hPhysPrim_Bin%d",ptbin)))->Fill(0.);

      if (IsDDaughter(d,trkMC,mcArray)) continue;
      if (fSoftPiCut && IsSoftPion_MCKine(d,trkMC,mcArray)) continue; //remove soft pions (if requestes, e.g. for templates)

      FillSparsePlots(mcArray,mInv,origD0,PDGD0,track,ptbin,kTrack,0); //fills for charged tracks

      //retrieving leading info...
      if(track->Pt() > highPt) {
        lead[0] = fCorrelatorTr->GetDeltaPhi();
        lead[1] = fCorrelatorTr->GetDeltaEta();
        lead[2] = fReadMC ? CheckTrackOrigin(mcArray,trkMC) : 0;
        highPt = track->Pt();
      }

    } // end of tracks loop
  } //end of event loop for fCorrelatorTr

 if(fKaonCorr) { //loops for Kcharg and K0

  if(fMixing) {
    NofEventsinPool = fCorrelatorKc->GetNofEventsInPool(); 
    if(!execPoolKc) {
      AliInfo("Mixed event analysis: K+/- pool is not ready");
      NofEventsinPool = 0;
    }
  }

  //Charged Kaons loop
  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++) {// loop on events in the pool; if it is SE analysis, stops at one (index not needed there)
    Bool_t analyzetracksKc = fCorrelatorKc->ProcessAssociatedTracks(jMix);
    if(!analyzetracksKc) {
      AliInfo("AliHFCorrelator::Cannot process the K+/- array");
      continue;
    }  

    for(Int_t iTrack = 0; iTrack<fCorrelatorKc->GetNofTracks(); iTrack++){ // looping on charged kaons candidates

      Bool_t runcorrelation = fCorrelatorKc->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* kCharg = fCorrelatorKc->GetAssociatedParticle();
      if(kCharg->GetLabel()<1) continue;
      if(kCharg->Pt() < fPtThreshLow.at(ptbin) || kCharg->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K
      if(TMath::Abs(kCharg->Eta())>0.8) continue; //discard tracks outside barrel (since it's kinematic MC and produces tracks all over rapidity region
      if(!fMixing) N_KCharg++;

      AliAODMCParticle *kChargMC = (AliAODMCParticle*)mcArray->At(kCharg->GetLabel());
      if(!kChargMC) continue;

      if (IsDDaughter(d,kChargMC,mcArray)) continue;
      FillSparsePlots(mcArray,mInv,origD0,PDGD0,kCharg,ptbin,kKCharg,0); //fills for charged tracks

    } // end of charged kaons loop
  } //end of event loop for fCorrelatorKc

  if(fMixing) {
    NofEventsinPool = fCorrelatorK0->GetNofEventsInPool(); 
    if(!execPoolK0) {
      AliInfo("Mixed event analysis: K0 pool is not ready");
      NofEventsinPool = 0;
    }
  }

  //K0 loop
  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++) {// loop on events in the pool; if it is SE analysis, stops at one (index not needed there)
    Bool_t analyzetracksK0 = fCorrelatorK0->ProcessAssociatedTracks(jMix);
    if(!analyzetracksK0) {
      AliInfo("AliHFCorrelator::Cannot process the K0 array");
      continue;
    }  

    for(Int_t iTrack = 0; iTrack<fCorrelatorK0->GetNofTracks(); iTrack++){ // looping on k0 candidates

      Bool_t runcorrelation = fCorrelatorK0->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* k0 = fCorrelatorK0->GetAssociatedParticle();
      if(k0->GetLabel()<1) continue;
      if(k0->Pt() < fPtThreshLow.at(ptbin) || k0->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K
      if(TMath::Abs(k0->Eta())>0.8) continue; //discard tracks outside barrel (since it's kinematic MC and produces tracks all over rapidity region
  
      AliAODMCParticle *k0MC = (AliAODMCParticle*)mcArray->At(k0->GetLabel());
      if(!k0MC) continue;

      if (IsDDaughter(d,k0MC,mcArray)) continue;
      FillSparsePlots(mcArray,mInv,origD0,PDGD0,k0,ptbin,kK0,0); //fills for charged tracks

      if(!fMixing) N_Kaons++;

    } // end of charged kaons loop
  } //end of event loop for fCorrelatorK0

 } //end of 'if(fKaonCorr)'

  Double_t fillSpLeadMC[4] = {lead[0],mD0,lead[1],0.4}; //mD0 = mD0bar = 1.864

  //leading track correlations fill
  if(!fMixing && !fSpeed) {
    if(d->GetPdgCode()==421 && (fIsSelectedCandidate==1||fIsSelectedCandidate==3)) { //D0
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadMC); //c and b D0
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC); //c or b D0
      if(origD0==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC);  
      if(origD0==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC);  
      if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadMC);  //non HF
    }
    if(d->GetPdgCode()==-421 && fIsSelectedCandidate>1) { //D0bar
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadMC);
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC); //c or b D0
      if(origD0==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC);  
      if(origD0==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC); 
      if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadMC);  //non HF
    }
  }
    //Fill of count histograms
  if (!fAlreadyFilled && !fMixing) { 
    ((TH1F*)fOutputStudy->FindObject("hist_Count_Charg"))->Fill(N_Charg);
    ((TH1F*)fOutputStudy->FindObject("hist_Count_Kcharg"))->Fill(N_KCharg);
    ((TH1F*)fOutputStudy->FindObject("hist_Count_K0"))->Fill(N_Kaons);
  }

  fAlreadyFilled=kTRUE; //at least a D0 analyzed in the event; distribution plots already filled

}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::FillSparsePlots(TClonesArray* mcArray, Double_t mInv[], Int_t origD0, Int_t PdgD0, AliReducedParticle* track, Int_t ptbin, Int_t type, Int_t softpiME, Double_t wg) {
  //
  //fills the THnSparse for correlations, calculating the variables
  //

  //Initialization of variables
  Double_t mD0, mD0bar, deltaphi = 0., deltaeta = 0.;
  mD0 = mInv[0];
  mD0bar = mInv[1];

  if (fReadMC && track->GetLabel()<1) return;
  if (fReadMC && !(AliAODMCParticle*)mcArray->At(track->GetLabel())) return;
  Double_t ptTrack = track->Pt();
  Double_t d0Track = type!=kK0 ? track->GetImpPar() : 0.;
  Double_t phiTr = track->Phi();
  Double_t etaTr = track->Eta();
  Double_t origTr = fReadMC ? CheckTrackOrigin(mcArray,(AliAODMCParticle*)mcArray->At(track->GetLabel())) : 0;

  TString part = "", orig = "";

  switch (type) {
    case(kTrack): {
      part = "Charg";
      deltaphi = fCorrelatorTr->GetDeltaPhi();
      deltaeta = fCorrelatorTr->GetDeltaEta();
      break;
    }
    case(kKCharg): {
      part = "Kcharg";
      deltaphi = fCorrelatorKc->GetDeltaPhi();
      deltaeta = fCorrelatorKc->GetDeltaEta();
      break;
    }
    case(kK0): {
      part = "K0";
      deltaphi = fCorrelatorK0->GetDeltaPhi();
      deltaeta = fCorrelatorK0->GetDeltaEta();
      break;
    }
  }
  
  if(fMixing == kSE) {

    //Fixes limits; needed to include overflow into THnSparse projections!
    Double_t pTorig = track->Pt();
    Double_t d0orig = track->GetImpPar();
    Double_t ptLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d",ptbin,fPoolNum)))->GetAxis(2)->GetXmax(); //all plots have same axes...
    Double_t displLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d",ptbin,fPoolNum)))->GetAxis(3)->GetXmax();
    Double_t EtaLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d",ptbin,fPoolNum)))->GetAxis(4)->GetXmax();
    if(ptTrack > ptLim_Sparse) ptTrack = ptLim_Sparse-0.01;
    if(d0Track > displLim_Sparse) d0Track = (displLim_Sparse-0.001);
    if(deltaeta > EtaLim_Sparse) deltaeta = EtaLim_Sparse-0.01;
    if(deltaeta < -EtaLim_Sparse) deltaeta = -EtaLim_Sparse+0.01;
  
    //variables for filling histos
    Double_t fillSpPhiD0[5] = {deltaphi,mD0,ptTrack,d0Track,deltaeta};
    Double_t fillSpPhiD0bar[5] = {deltaphi,mD0bar,ptTrack,d0Track,deltaeta};
    Double_t fillSpWeigD0[5] = {deltaphi,mD0,deltaeta,ptTrack};
    Double_t fillSpWeigD0bar[5] = {deltaphi,mD0bar,deltaeta,ptTrack};

    Bool_t allowD0 = 0;
    Bool_t allowD0bar = 0;
    if(fSpeed) { //filling of sidebands in speed mode: 1 bin for LSB, 1 for RSB, no filling outside signal region and SB
      if(ptbin<9) {	    
        if(mD0 > fSignLeft_LowPt && mD0 < fSignRight_LowPt) allowD0 = 1;
        if(mD0bar > fSignLeft_LowPt && mD0bar < fSignRight_LowPt) allowD0bar = 1;
      } else {
        if(mD0 > fSignLeft_HighPt && mD0 < fSignRight_HighPt) allowD0 = 1;
        if(mD0bar > fSignLeft_HighPt && mD0bar < fSignRight_HighPt) allowD0bar = 1;
      }
      if(mD0 > fLSBLowLim.at(ptbin) && mD0 < fLSBUppLim.at(ptbin)) {allowD0 = 1; fillSpPhiD0[1] = 1.60; fillSpWeigD0[1] = 1.60;} //in LSB bin!
      if(mD0bar > fLSBLowLim.at(ptbin) && mD0bar < fLSBUppLim.at(ptbin)) {allowD0bar = 1; fillSpPhiD0bar[1] = 1.60; fillSpWeigD0bar[1] = 1.60;} //in LSB bin!
      if(mD0 > fRSBLowLim.at(ptbin) && mD0 < fRSBUppLim.at(ptbin)) {allowD0 = 1; fillSpPhiD0[1] = 2.18; fillSpWeigD0[1] = 2.18;} //in RSB bin!
      if(mD0bar > fRSBLowLim.at(ptbin) && mD0bar < fRSBUppLim.at(ptbin)) {allowD0bar = 1; fillSpPhiD0bar[1] = 2.18; fillSpWeigD0bar[1] = 2.18;} //in RSB bin!
    } //in this way if sidebands overlap with signal range in Mass axis, those overlapping bins will be void. But this creates no problems...
    else if(!fSpeed) { // Full Minv range in THnSparse!
      if((fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3)) allowD0 = 1;   
      if((fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3)) allowD0bar = 1;
    }
    
    if(fReadMC == 0) {
      //sparse fill for data (tracks, K+-, K0) + weighted
      if((fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) && allowD0) { //D0
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0,wg);
        if(!fSpeed) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigD0,pTorig*wg);
      }
      if((fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) && allowD0bar) { //D0bar
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0bar,wg);
        if(!fSpeed) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigD0bar,pTorig*wg);
      }
      if(!fAlreadyFilled) {
 	   ((TH1F*)fOutputStudy->FindObject(Form("hist_Pt_%s_Bin%d",part.Data(),ptbin)))->Fill(pTorig);
 	   ((TH1F*)fOutputStudy->FindObject(Form("hist_PhiDistr_%s",part.Data())))->Fill(phiTr);
 	   ((TH1F*)fOutputStudy->FindObject(Form("hist_EtaDistr_%s",part.Data())))->Fill(etaTr);
 	   ((TH2F*)fOutputStudy->FindObject(Form("hist_PhiVsEtaDistr_%s",part.Data())))->Fill(phiTr,etaTr);
      }
    }

    if(fReadMC) {

      if(origD0==4) {orig = "_From_c";} else {orig = "_From_b";}

      //sparse fill for data (tracks, K+-, K0) + weighted
      if(PdgD0==421 && (fIsSelectedCandidate==1||fIsSelectedCandidate==3)) { //D0 (from MCTruth)
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0,wg);
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0,wg);
         if(origD0==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0,wg);
         if(origD0==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0,wg);
         if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_NonHF_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0,wg);
         if(!fSpeed) {
           ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigD0,pTorig*wg);
           ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigD0,pTorig*wg);
           if(origD0==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigD0,pTorig*wg);
           if(origD0==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigD0,pTorig*wg);
           if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_NonHF_Bin%d",ptbin)))->Fill(fillSpWeigD0,pTorig*wg);
         }  
      }
      if(PdgD0==-421 && fIsSelectedCandidate>1) { //D0bar
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0bar,wg);
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0bar,wg);
         if(origD0==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0bar,wg);
         if(origD0==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0bar,wg); 
         if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_NonHF_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0bar,wg);
         if(!fSpeed) {
           ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigD0bar,pTorig*wg);
           ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigD0bar,pTorig*wg);
           if(origD0==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigD0bar,pTorig*wg);
           if(origD0==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigD0bar,pTorig*wg);
           if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_NonHF_Bin%d",ptbin)))->Fill(fillSpWeigD0bar,pTorig*wg);
         }
      } 
      if(!fAlreadyFilled) {
	((TH1F*)fOutputStudy->FindObject(Form("histDispl_%s_Bin%d",part.Data(),ptbin)))->Fill(d0orig); //Fills displacement histos
        if (origTr>=1&&origTr<=8) ((TH1F*)fOutputStudy->FindObject(Form("histDispl_%s_HF_Bin%d",part.Data(),ptbin)))->Fill(d0orig);
        if (origTr>=1&&origTr<=8) ((TH1F*)fOutputStudy->FindObject(Form("histDispl_%s_HF%s_Bin%d",part.Data(),orig.Data(),ptbin)))->Fill(d0orig);
        ((TH1F*)fOutputStudy->FindObject(Form("histDispl_%s%s_Bin%d",part.Data(),orig.Data(),ptbin)))->Fill(d0orig); //Fills displacement histos
        ((TH1F*)fOutputStudy->FindObject(Form("hist_Pt_%s_Bin%d",part.Data(),ptbin)))->Fill(pTorig);
        ((TH1F*)fOutputStudy->FindObject(Form("histOrig_%s_Bin%d",part.Data(),ptbin)))->Fill(origTr);
 	((TH1F*)fOutputStudy->FindObject(Form("hist_PhiDistr_%s",part.Data())))->Fill(phiTr);
 	((TH1F*)fOutputStudy->FindObject(Form("hist_PhiDistr_%s",part.Data())))->Fill(phiTr);
 	((TH1F*)fOutputStudy->FindObject(Form("hist_EtaDistr_%s",part.Data())))->Fill(etaTr);
 	((TH2F*)fOutputStudy->FindObject(Form("hist_PhiVsEtaDistr_%s",part.Data())))->Fill(phiTr,etaTr);
      }
    }//end MC case

  } //end of SE fill

  if(fMixing == kME) {

    //Fixes limits; needed to include overflow into THnSparse projections!
    Double_t EtaLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d_EvMix",ptbin,fPoolNum)))->GetAxis(2)->GetXmax();
    if(deltaeta > EtaLim_Sparse) deltaeta = EtaLim_Sparse-0.01;
    if(deltaeta < -EtaLim_Sparse) deltaeta = -EtaLim_Sparse+0.01;
    Double_t ptLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d_EvMix",ptbin,fPoolNum)))->GetAxis(3)->GetXmax(); //all plots have same axes...
    if(ptTrack > ptLim_Sparse) ptTrack = ptLim_Sparse-0.01;

    //variables for filling histos
    Double_t fillSpPhiD0[5] = {deltaphi,mD0,deltaeta,0.4,0}; //dummy for ME threshold! unless explicitly set by flag...
    Double_t fillSpPhiD0bar[5] = {deltaphi,mD0bar,deltaeta,0.4,0};
    if(fMEAxisThresh) {
      fillSpPhiD0[3] = ptTrack;
      fillSpPhiD0bar[3] = ptTrack;
    }
    if(softpiME==1) { //it's a softPi in the ME analysis! Fill it in the dedicated slice of ME THnSparse
      fillSpPhiD0[4] = 1;
      fillSpPhiD0bar[4] = 1;
    }

    Bool_t allowD0 = 0;
    Bool_t allowD0bar = 0;
    if(fSpeed) { //filling of sidebands in speed mode: 1 bin for LSB, 1 for RSB, no filling outside signal region and SB
      if(ptbin<9) {	    
        if(mD0 > fSignLeft_LowPt && mD0 < fSignRight_LowPt) allowD0 = 1;
        if(mD0bar > fSignLeft_LowPt && mD0bar < fSignRight_LowPt) allowD0bar = 1;
      } else {
        if(mD0 > fSignLeft_HighPt && mD0 < fSignRight_HighPt) allowD0 = 1;
        if(mD0bar > fSignLeft_HighPt && mD0bar < fSignRight_HighPt) allowD0bar = 1;
      }
      if(mD0 > fLSBLowLim.at(ptbin) && mD0 < fLSBUppLim.at(ptbin)) {allowD0 = 1; fillSpPhiD0[1] = 1.60;} //in LSB bin!
      if(mD0bar > fLSBLowLim.at(ptbin) && mD0bar < fLSBUppLim.at(ptbin)) {allowD0bar = 1; fillSpPhiD0bar[1] = 1.60;} //in LSB bin!
      if(mD0 > fRSBLowLim.at(ptbin) && mD0 < fRSBUppLim.at(ptbin)) {allowD0 = 1; fillSpPhiD0[1] = 2.18;} //in RSB bin!
      if(mD0bar > fRSBLowLim.at(ptbin) && mD0bar < fRSBUppLim.at(ptbin)) {allowD0bar = 1; fillSpPhiD0bar[1] = 2.18;} //in RSB bin!
    } //in this way if sidebands overlap with signal range in Mass axis, those overlapping bins will be void. But this creates no problems...
    else if(!fSpeed) { // Full Minv range in THnSparse!
      if((fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3)) allowD0 = 1;   
      if((fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3)) allowD0bar = 1;
    }    
    
    if(fReadMC == 0) {
      //sparse fill for data (tracks, K+-, K0)
      if((fIsSelectedCandidate == 1||fIsSelectedCandidate == 3) && allowD0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d_EvMix",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0,wg);
      if((fIsSelectedCandidate == 2||fIsSelectedCandidate == 3) && allowD0bar) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d_EvMix",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0bar,wg);
    }
    if(fReadMC == 1) {
      //sparse fill for data (tracks, K+-, K0)
      if(PdgD0==421 && (fIsSelectedCandidate==1||fIsSelectedCandidate==3))  ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d_EvMix",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0,wg);
      if(PdgD0==-421 && fIsSelectedCandidate>1) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d_EvMix",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiD0bar,wg);
    }//end MC case

  } //end of ME fill
  
  return;
}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSED0Correlations::CheckTrackOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  // checks on particle (#) origin:
  // 0) Not HF
  // 1) D->#
  // 2) D->X->#
  // 3) c hadronization
  // 4) B->#
  // 5) B->X-># (X!=D)
  // 6) B->D->#
  // 7) B->D->X->#
  // 8) b hadronization
  //
  if(fDebug>2) printf("AliAnalysisTaskSED0Correlations::CheckTrkOrigin() \n");
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isDdaugh=kFALSE;
  Bool_t isDchaindaugh=kFALSE;
  Bool_t isBdaugh=kFALSE;
  Bool_t isBchaindaugh=kFALSE;
  Bool_t isQuarkFound=kFALSE;

  if (mother<0) return -1;
  while (mother >= 0){
    istep++;
    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcMoth){
      pdgGranma = mcMoth->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isBchaindaugh=kTRUE;
        if(istep==1) isBdaugh=kTRUE;
      }
      if ((abspdgGranma > 400 && abspdgGranma < 500) || (abspdgGranma > 4000 && abspdgGranma < 5000)){
	isDchaindaugh=kTRUE;
        if(istep==1) isDdaugh=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) {isQuarkFound=kTRUE; if(abspdgGranma==5) isFromB = kTRUE;}
      mother = mcMoth->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      return -1;
    }
  }

  //decides what to return based on the flag status
  if(isQuarkFound) {
    if(!isFromB) {  //charm
      if(isDdaugh) return 1; //charm immediate
      else if(isDchaindaugh) return 2; //charm chain
      else return 3; //charm hadronization
    }
    else { //beauty
      if(isBdaugh) return 4; //b immediate
      else if(isBchaindaugh) { //b chain
        if(isDchaindaugh) {
          if(isDdaugh) return 6; //d immediate
          return 7; //d chain
          }
        else return 5; //b, not d
      }
      else return 8; //b hadronization
    }
  }
  else if(!isDdaugh && !isDchaindaugh && !isBdaugh && !isBchaindaugh) return 0; //no HF decay 
     //in this case, it's !isQuarkFound, but not in 100% cases it's a non HF particle!
     //rarely you can find a D/B meson which comes from a -1! It isn't a Non-HF, in that case! And I'll return -1...

  return -1; //some problem spotted
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSED0Correlations::IsDDaughter(AliAODMCParticle* d, AliAODMCParticle* track, TClonesArray* mcArray) const {

  //Daughter removal in MCKine case
  Bool_t isDaughter = kFALSE;
  Int_t labelD0 = d->GetLabel();

  Int_t mother = track->GetMother();

  //Loop on the mothers to find the D0 label (it must be the trigger D0, not a generic D0!)
  while (mother > 0){
    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(mcArray->At(mother)); //it's the mother of the track!
    if (mcMoth){
      if (mcMoth->GetLabel() == labelD0) isDaughter = kTRUE;
      mother = mcMoth->GetMother(); //goes back by one
    } else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }

  return isDaughter;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSED0Correlations::PtBinCorr(Double_t pt) const {
  //
  //give the pt bin where the pt lies.
  //
  Int_t ptbin=-1;
  if(pt<fBinLimsCorr.at(0)) return ptbin; //out of bounds
  
  Int_t i = 0;
  while(pt>fBinLimsCorr.at(i)) {ptbin=i; i++;}
  
  return ptbin;
}

//---------------------------------------------------------------------------
Bool_t AliAnalysisTaskSED0Correlations::SelectV0(AliAODv0* v0, AliAODVertex *vtx, Int_t opt, Int_t idArrayV0[][2]) const
{
  //
  // Selection for K0 hypotheses
  // options: 1 = selects mass invariant about 3 sigma inside the peak + threshold of 0.3 GeV
  // 	      2 = no previous selections

  if(!fCutsTracks->IsKZeroSelected(v0,vtx)) return kFALSE;

  AliAODTrack *v0Daug1 = (AliAODTrack*)v0->GetDaughter(0);
  if(!v0Daug1) {
    AliWarning("Error in casting to AOD track. Not a standard AOD?");
    return kFALSE;
  }
  AliAODTrack *v0Daug2 = (AliAODTrack*)v0->GetDaughter(1);
  if(!v0Daug2) {
    AliWarning("Error in casting to AOD track. Not a standard AOD?");
    return kFALSE;
  }

  if(opt==1) { //additional cuts for correlations (V0 has to be closer than 3 sigma from K0 mass)
    if(TMath::Abs(v0->MassK0Short()-0.4976) > 3*0.004) return kFALSE;
  }

  //This part removes double counting for swapped tracks!
  Int_t i = 0;  //while loop (until the last-written entry pair of ID!
  while(idArrayV0[i][0]!=-2 && idArrayV0[i][1]!=-2) {
    if((v0Daug1->GetID()==idArrayV0[i][0] && v0Daug2->GetID()==idArrayV0[i][1])||
       (v0Daug1->GetID()==idArrayV0[i][1] && v0Daug2->GetID()==idArrayV0[i][0])) return kFALSE;
    i++;
  }
  idArrayV0[i][0]=v0Daug1->GetID();
  idArrayV0[i][1]=v0Daug2->GetID();

  return kTRUE;
}

//---------------------------------------------------------------------------
Bool_t AliAnalysisTaskSED0Correlations::IsSoftPion_MCKine(AliAODMCParticle* d, AliAODMCParticle* track, TClonesArray* arrayMC) const
{
  //
  // Removes soft pions in Kine

  //Daughter removal in MCKine case
  Bool_t isSoftPi = kFALSE;
  Int_t labelD0 = d->GetLabel();

  Int_t mother = track->GetMother();
  if(mother<0) return isSoftPi; //safety check

  AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother)); //it's the mother of the track!
  if(!mcMoth){
    return isSoftPi;
  }
  if(TMath::Abs(mcMoth->GetPdgCode())==413 && mcMoth->GetNDaughters()==2) { //mother is D* with 2 daughs
    Int_t labdau1 = mcMoth->GetDaughter(0);
    Int_t labdau2 = mcMoth->GetDaughter(1);
    AliAODMCParticle* dau1 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(labdau1));
    AliAODMCParticle* dau2 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(labdau2));
    if(!dau1 || !dau2) return isSoftPi; //safety check
    if(dau1->GetLabel()==labelD0 || dau2->GetLabel()==labelD0) { //one of the daughs is the D0 trigger
      if((TMath::Abs(dau1->GetPdgCode())==421 && TMath::Abs(dau2->GetPdgCode())==211)||(TMath::Abs(dau1->GetPdgCode())==211 && TMath::Abs(dau2->GetPdgCode())==421)) {
	isSoftPi = kTRUE; //ok, soft pion was found
	return isSoftPi;
      }
    }
  } 

  return isSoftPi;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::FillTreeD0(AliAODRecoDecayHF2Prong* d, AliAODEvent* aod) {

  Int_t ptbin = PtBinCorr(d->Pt());
  if(ptbin < 0) return;

  Bool_t allowD0 = 0;
  Bool_t allowD0bar = 0;
  Double_t mD0, mD0bar;
  d->InvMassD0(mD0,mD0bar);
  fBranchD->invMass_D = 0;

  Float_t centEv = -9;
  if(fCutsD0->GetUseCentrality()) centEv = fCutsD0->GetCentrality(aod); //get event centrality with current estimator

  if(fSpeed) { //filling of sidebands in speed mode: 1 bin for LSB, 1 for RSB, no filling outside signal region and SB
    if(ptbin<9) {	    
      if(mD0 > fSignLeft_LowPt && mD0 < fSignRight_LowPt) allowD0 = 1;
      if(mD0bar > fSignLeft_LowPt && mD0bar < fSignRight_LowPt) allowD0bar = 1;
    } else {
      if(mD0 > fSignLeft_HighPt && mD0 < fSignRight_HighPt) allowD0 = 1;
      if(mD0bar > fSignLeft_HighPt && mD0bar < fSignRight_HighPt) allowD0bar = 1;
    }
    if(mD0 > fLSBLowLim.at(ptbin) && mD0 < fLSBUppLim.at(ptbin)) allowD0 = 1; //in LSB bin!
    if(mD0bar > fLSBLowLim.at(ptbin) && mD0bar < fLSBUppLim.at(ptbin)) allowD0bar = 1; //in LSB bin!
    if(mD0 > fRSBLowLim.at(ptbin) && mD0 < fRSBUppLim.at(ptbin)) allowD0 = 1; //in RSB bin!
    if(mD0bar > fRSBLowLim.at(ptbin) && mD0bar < fRSBUppLim.at(ptbin)) allowD0bar = 1; //in RSB bin!
  } //in this way if sidebands overlap with signal range in Mass axis, those overlapping bins will be void. But this creates no problems...
  else if(!fSpeed) { // Full Minv range in THnSparse!
    if((fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3)) allowD0 = 1;   
    if((fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3)) allowD0bar = 1;
  }   

  //Fill TTree for accepted candidates
  //**** NOTE: each candidate has a IDtrig, regardless of how many hypotheses are allowed (if both, the TTree is filled 2 times, but both times with the same IDtrig) ****//
  if(fIsSelectedCandidate>0 && (allowD0 || allowD0bar)) {
    ResetBranchD();
    fBranchD->phi_D = (Float_t)d->Phi();
    fBranchD->eta_D = (Float_t)d->Eta();
    fBranchD->pT_D = (Float_t)d->Pt();
    fBranchD->mult_D = (Float_t)fMultEv;
    fBranchD->zVtx_D = (Float_t)fzVtx;
    fBranchD->cent_D = (Float_t)centEv;
    fBranchD->period_D = (UInt_t)aod->GetPeriodNumber();
    fBranchD->orbit_D = (UInt_t)aod->GetOrbitNumber();
    fBranchD->BC_D = (UShort_t)aod->GetBunchCrossNumber();
    fBranchD->IDtrig_D = (Short_t)fNtrigD;
    fBranchD->sel_D = (Short_t)1; //******DUMMY FOR THE MOMENT******* - To be used for multiple selection fills (2^n = 1 if selction n is ok)
    fBranchD->pXdaug1_D = (Float_t)((AliVTrack*)d->GetDaughter(0))->Px();
    fBranchD->pXdaug2_D = (Float_t)((AliVTrack*)d->GetDaughter(1))->Px();
    fBranchD->pYdaug1_D = (Float_t)((AliVTrack*)d->GetDaughter(0))->Py();
    fBranchD->pYdaug2_D = (Float_t)((AliVTrack*)d->GetDaughter(1))->Py();
    fBranchD->pZdaug1_D = (Float_t)((AliVTrack*)d->GetDaughter(0))->Pz();
    fBranchD->pZdaug2_D = (Float_t)((AliVTrack*)d->GetDaughter(1))->Pz();
    fBranchD->hyp_D = (UShort_t)fIsSelectedCandidate;
    if((fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) && allowD0) {
      fBranchD->invMass_D = (Float_t)mD0;
      fBranchD->hyp_D = (UShort_t)1;
      fTreeD->Fill();
    }
    if((fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) && allowD0bar) {
      fBranchD->invMass_D = (Float_t)mD0bar;
      fBranchD->hyp_D = (UShort_t)2;
      fTreeD->Fill();
    }

    //store daughter ID, so that are not included in TTree of tracks, and tags the corresponding number of the D-meson trigger
    fDaughTrackID.push_back(((AliVTrack*)d->GetDaughter(0))->GetID());
    fDaughTrigNum.push_back(fNtrigD);
    fDaughTrackID.push_back(((AliVTrack*)d->GetDaughter(1))->GetID());
    fDaughTrigNum.push_back(fNtrigD);

    if(!fTrackArrayFilled) {
      fTrackArray = fCorrelatorTr->AcceptAndReduceTracks(aod); //track selection, needed for soft pion rejection (for the first trigger only)
      fTrackArrayFilled = kTRUE;
    }
 
    //soft pion rejection (for SE, in ME it's done in the AliHFOfflineCorrelator class
    for(Int_t iTrack = 0; iTrack<fTrackArray->GetEntriesFast(); iTrack++) { // looping on track candidates
      AliReducedParticle* track = (AliReducedParticle*)fTrackArray->At(iTrack);
      if(fSoftPiCut && !track->CheckSoftPi()) {  //identifies soft pions for the trigger under analysis and associate it to the track
        fSoftPiTrackID.push_back(track->GetID()); //tags tha track id as a soft pion
        fSoftPiTrigNum.push_back(fNtrigD); //identifies whose trigger the track is a softpion
        fNsoftPi++;
      }
    }

    fNtrigD++; //increase by 1 the index of D0 triggers in the event
  } //end of if for tree filling

  fAlreadyFilled=kTRUE; //A D0 was selected in the event (even if can be also outside mass range)! Enables saving the tracks of the event in the ME offline approach

  return;

}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::FillTreeTracks(AliAODEvent* aod) {

  if(!fTrackArrayFilled) {
    fTrackArray = fCorrelatorTr->AcceptAndReduceTracks(aod); //track selection, if not already done in FillTreeD0
    fTrackArrayFilled = kTRUE;
  }

  Float_t centEv = -9;
  if(fCutsD0->GetUseCentrality()) centEv = fCutsD0->GetCentrality(aod); //get event centrality with current estimator
  
  for(Int_t iTrack = 0; iTrack<fTrackArray->GetEntriesFast(); iTrack++){ // looping on track candidates
      
    AliReducedParticle* track = (AliReducedParticle*)fTrackArray->At(iTrack);
/*
      if(fReadMC) {
        AliAODMCParticle* trkKine = (AliAODMCParticle*)mcArray->At(track->GetLabel());
        if (!trkKine) continue;
        if (!trkKine->IsPhysicalPrimary()) {
 	  ((TH1F*)fOutputStudy->FindObject(Form("hPhysPrim_Bin%d",ptbin)))->Fill(1.);  
  	  continue; //reject the Reco track if correspondent Kine track is not primary
        } else ((TH1F*)fOutputStudy->FindObject(Form("hPhysPrim_Bin%d",ptbin)))->Fill(0.);
      }
*/

    //skip D-meson trigger daughters
    Short_t trigID = -1, trigID2 = -1, trigID3 = -1, trigID4 = -1;
    Int_t FoundTrig=0;
    Bool_t trackIsTrig=kFALSE;
    for(Int_t iID=0; iID<(int)fDaughTrackID.size(); iID++) {
      trackIsTrig=kFALSE; //reset flag to signal that the track is a trigger
      if(FoundTrig==0 && track->GetID() == fDaughTrackID.at(iID)) {
	trigID = fDaughTrigNum.at(iID); //associates corresponding trigID to daughters
	trackIsTrig=kTRUE;
      }
      if(FoundTrig==1 && track->GetID() == fDaughTrackID.at(iID)) {
	trigID2 = fDaughTrigNum.at(iID); //associates corresponding trigID to daughters
	trackIsTrig=kTRUE;
      }
      if(FoundTrig==2 && track->GetID() == fDaughTrackID.at(iID)) {
	trigID3 = fDaughTrigNum.at(iID); //associates corresponding trigID to daughters
	trackIsTrig=kTRUE;
      }	  
      if(FoundTrig==3 && track->GetID() == fDaughTrackID.at(iID)) {
	trigID4 = fDaughTrigNum.at(iID); //associates corresponding trigID to daughters
	trackIsTrig=kTRUE;
      }	
      if(trackIsTrig==kTRUE) FoundTrig++; //if track is a trigger, next track has to be stored in another position, for IDTrig!
    }
    
    //tags soft pions in the same way as daughter tracks (for the corresponding trigger)
    Bool_t trackIssoftPi=kFALSE;
    for(Int_t iID=0; iID<(int)fSoftPiTrackID.size(); iID++) {
      trackIssoftPi=kFALSE; //reset flag to signal that the track is a soft pion
      if(FoundTrig==0 && track->GetID() == fSoftPiTrackID.at(iID)) {
	trigID = fSoftPiTrigNum.at(iID); //associates corresponding trigID to daughters
	trackIssoftPi=kTRUE;
      }
      if(FoundTrig==1 && track->GetID() == fSoftPiTrackID.at(iID)) {
	trigID2 = fSoftPiTrigNum.at(iID); //associates corresponding trigID to daughters
	trackIssoftPi=kTRUE;
      }
      if(FoundTrig==2 && track->GetID() == fSoftPiTrackID.at(iID)) {
	trigID3 = fSoftPiTrigNum.at(iID); //associates corresponding trigID to daughters
	trackIssoftPi=kTRUE;
      } 
      if(FoundTrig==3 && track->GetID() == fSoftPiTrackID.at(iID)) {
	trigID4 = fSoftPiTrigNum.at(iID); //associates corresponding trigID to daughters
	trackIssoftPi=kTRUE;
      }	 
      if(trackIssoftPi==kTRUE) FoundTrig++;	//if track is a soft pion, next track has to be stored in another position, for IDTrig!
    }
    
    if(!AcceptTrackForMEOffline(track->Pt())) continue;

    //Fill TTree for accepted candidates
    ResetBranchTracks();
    fBranchTr->phi_Tr = (Float_t)track->Phi();
    fBranchTr->eta_Tr = (Float_t)track->Eta();
    fBranchTr->pT_Tr = (Float_t)track->Pt();
    fBranchTr->mult_Tr = (Float_t)fMultEv;
    fBranchTr->zVtx_Tr = (Float_t)fzVtx;
    fBranchTr->cent_Tr = (Float_t)centEv;
    fBranchTr->period_Tr = (UInt_t)aod->GetPeriodNumber();
    fBranchTr->orbit_Tr = (UInt_t)aod->GetOrbitNumber();
    fBranchTr->BC_Tr = (UShort_t)aod->GetBunchCrossNumber();
    fBranchTr->IDtrig_Tr = (Short_t)trigID;
    fBranchTr->IDtrig2_Tr = (Short_t)trigID2;
    fBranchTr->IDtrig3_Tr = (Short_t)trigID3;
    fBranchTr->IDtrig4_Tr = (Short_t)trigID4;
    fBranchTr->sel_Tr = (Short_t)1; //******DUMMY FOR THE MOMENT******* - To be used for multiple selection fills (2^n = 1 if selction n is ok)
    fTreeTr->Fill();

  } //end of track loop

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::FillTreeD0ForCutOptim(AliAODRecoDecayHF2Prong* d, AliAODEvent* aod) {

  Int_t ptbin = PtBinCorr(d->Pt());
  if(ptbin < 0) return;

  Float_t centEv = -9;
  if(fCutsD0->GetUseCentrality()) centEv = fCutsD0->GetCentrality(aod); //get event centrality with current estimator

  //recalculate vertex w/o daughters
  AliAODVertex *origownvtx=0x0;
  if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
  if(!fCutsD0->RecalcOwnPrimaryVtx(d,aod)) {
     fCutsD0->CleanOwnPrimaryVtx(d,aod,origownvtx);
     return;
  }

  //Preliminary vars
  Double_t mD0, mD0bar;
  d->InvMassD0(mD0,mD0bar);
  fBranchDCutVars->invMass = 0;
  Double_t ctsD0, ctsD0bar;
  d->CosThetaStarD0(ctsD0,ctsD0bar);

  //Topomatic
  Double_t dd0max=0.;
  for(Int_t ipr=0; ipr<2; ipr++) {
    Double_t diffIP, errdiffIP;
    d->Getd0MeasMinusExpProng(ipr,aod->GetMagneticField(),diffIP,errdiffIP);
    Double_t normdd0=0.;
    if(errdiffIP>0.) normdd0=diffIP/errdiffIP;
    if(ipr==0) dd0max=normdd0;
    else if(TMath::Abs(normdd0)>TMath::Abs(dd0max)) dd0max=normdd0;
  }

// printf("Centralità = %f, %f (ZNA), %f (V0M)\n",fCutsD0->GetCentrality(aod),fCutsD0->GetCentrality(aod,AliRDHFCuts::kCentZNA),fCutsD0->GetCentrality(aod,AliRDHFCuts::kCentV0M)); getchar();

  //Fill TTree for accepted candidates
  //if both hypotheses are ok, the TTree is filled 2 times, with the different cut values
  if(fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) {
    ResetBranchDForCutOptim();
    fBranchDCutVars->invMass = (Float_t)mD0;
    fBranchDCutVars->cent = (Float_t)centEv;
    fBranchDCutVars->pt = (Float_t)d->Pt();
    fBranchDCutVars->dca = (Float_t)d->GetDCA();
    fBranchDCutVars->cosThSt = (Float_t)TMath::Abs(ctsD0);
    fBranchDCutVars->pTk = (Float_t)d->Pt2Prong(1);
    fBranchDCutVars->pTpi = (Float_t)d->Pt2Prong(0);
    fBranchDCutVars->d0k = (Float_t)TMath::Abs(d->Getd0Prong(1)); 
    fBranchDCutVars->d0pi = (Float_t)TMath::Abs(d->Getd0Prong(0));
    fBranchDCutVars->d0xd0 = (Float_t)d->Prodd0d0();
    fBranchDCutVars->cosThPt = (Float_t)d->CosPointingAngle();
    fBranchDCutVars->normLxy = (Float_t)d->NormalizedDecayLengthXY();
    fBranchDCutVars->topom = (Float_t)TMath::Abs(dd0max);
    fTreeD->Fill();
  } //end of if for tree filling

  if(fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) {
    ResetBranchDForCutOptim();
    fBranchDCutVars->invMass = (Float_t)mD0bar;
    fBranchDCutVars->cent = (Float_t)centEv;
    fBranchDCutVars->pt = (Float_t)d->Pt();
    fBranchDCutVars->dca = (Float_t)d->GetDCA();
    fBranchDCutVars->cosThSt = (Float_t)TMath::Abs(ctsD0bar);
    fBranchDCutVars->pTk = (Float_t)d->Pt2Prong(0);
    fBranchDCutVars->pTpi = (Float_t)d->Pt2Prong(1);
    fBranchDCutVars->d0k = (Float_t)TMath::Abs(d->Getd0Prong(0)); 
    fBranchDCutVars->d0pi = (Float_t)TMath::Abs(d->Getd0Prong(1));
    fBranchDCutVars->d0xd0 = (Float_t)d->Prodd0d0();
    fBranchDCutVars->cosThPt = (Float_t)d->CosPointingAngle();
    fBranchDCutVars->normLxy = (Float_t)d->NormalizedDecayLengthXY();
    fBranchDCutVars->topom = (Float_t)TMath::Abs(dd0max);
    fTreeD->Fill();
  } //end of if for tree filling

  // unset recalculated primary vertex when not needed any more
  fCutsD0->CleanOwnPrimaryVtx(d,aod,origownvtx);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::ResetBranchD() {

  fBranchD->phi_D = 0.;
  fBranchD->eta_D = 0.;
  fBranchD->pT_D = 0.;
  fBranchD->invMass_D = 0.;
  fBranchD->mult_D = 0.;
  fBranchD->zVtx_D = 0.;
  fBranchD->cent_D = 0.;
  fBranchD->period_D = 0;
  fBranchD->orbit_D = 0;
  fBranchD->BC_D = 0;
  fBranchD->IDtrig_D = 0;
  fBranchD->sel_D = 0;
  fBranchD->pXdaug1_D = 0.;
  fBranchD->pXdaug2_D = 0.;
  fBranchD->pYdaug1_D = 0.;
  fBranchD->pYdaug2_D = 0.;
  fBranchD->pZdaug1_D = 0.;
  fBranchD->pZdaug2_D = 0.;
  fBranchD->hyp_D = 4;  //it can go from 0 to 2, so if you have a 4 in the TTree something's wrong in the filling...

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::ResetBranchTracks() {

  fBranchTr->phi_Tr = 0.;
  fBranchTr->eta_Tr = 0.;
  fBranchTr->pT_Tr = 0.;
  fBranchTr->mult_Tr = 0.;
  fBranchTr->zVtx_Tr = 0.;
  fBranchTr->cent_Tr = 0.;
  fBranchTr->period_Tr = 0;
  fBranchTr->orbit_Tr = 0;
  fBranchTr->BC_Tr = 0;
  fBranchTr->IDtrig_Tr = 0;
  fBranchTr->IDtrig2_Tr = 0;
  fBranchTr->IDtrig3_Tr = 0;
  fBranchTr->IDtrig4_Tr = 0; 
  fBranchTr->sel_Tr = 0;

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::ResetBranchDForCutOptim() {

  fBranchDCutVars->invMass = 0.;
  fBranchDCutVars->cent = 0.;
  fBranchDCutVars->pt = 0.;
  fBranchDCutVars->dca = 0.;
  fBranchDCutVars->cosThSt = 0.;
  fBranchDCutVars->pTk = 0.;
  fBranchDCutVars->pTpi = 0.;
  fBranchDCutVars->d0k = 0.;
  fBranchDCutVars->d0pi = 0.;
  fBranchDCutVars->d0xd0 = 0.;
  fBranchDCutVars->cosThPt = 0.;
  fBranchDCutVars->normLxy = 0;
  fBranchDCutVars->topom = 0.;
  
  return;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSED0Correlations::AcceptTrackForMEOffline(Double_t pt) {

  pt = pt*1000-(long)(pt*1000); //to extract decimals from 4° onwards - I.e. a.bcdefghi --> 0.efghi
  pt = (long)(pt*100); //to take only 4° and 5° decimal digits - I.e. 0.efghi --> ef

  if(pt > fFractAccME) return kFALSE; //reject track
  else return kTRUE; //accept track for offline ME

}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::PrintBinsAndLimits() {

  cout << "--------------------------\n";
  cout << "PtBins = " << fNPtBinsCorr << "\n";
  cout << "PtBin limits--------------\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << "Bin "<<i+1<<" = "<<fBinLimsCorr.at(i)<<" to "<<fBinLimsCorr.at(i+1)<<"\n";
  }
  cout << "\n--------------------------\n";
  cout << "PtBin tresh. tracks low---\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fPtThreshLow.at(i) << ", ";
  }
  cout << "\nPtBin tresh. tracks up----\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fPtThreshUp.at(i) << ", ";
  }
  cout << "\nSB limits (LSBLow)----\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fLSBLowLim.at(i) << ", ";
  }
  cout << "\nSB limits (LSBUpp)----\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fLSBUppLim.at(i) << ", ";
  }
  cout << "\nSB limits (RSBLow)----\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fRSBLowLim.at(i) << ", ";
  }
  cout << "\nSB limits (RSBUpp)----\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fRSBUppLim.at(i) << ", ";
  }

  cout << "\n--------------------------\n";
  cout << "D0 Eta cut for Correl = "<<fEtaForCorrel<<"\n";
  cout << "--------------------------\n";
  cout << "MC Truth = "<<fReadMC<<" - MC Reco: Trk = "<<fRecoTr<<", D0 = "<<fRecoD0<<"\n";
  cout << "--------------------------\n";
  cout << "Sel of Event tpye (PP,GS,FE,...)= "<<fSelEvType<<"\n";
  cout << "--------------------------\n";
  cout << "Ev Mixing = "<<fMixing<<"\n";
  cout << "--------------------------\n";
  cout << "ME thresh axis = "<<fMEAxisThresh<<"\n";
  cout << "--------------------------\n";
  cout << "Soft Pi Cut = "<<fSoftPiCut<<"\n";
  cout << "--------------------------\n";
  cout << "Speed (1 SBL/SBR bin) = "<<fSpeed<<"\n";
  cout << "--------------------------\n";
  cout << "All entries in Pool0 = "<<fMergePools<<"\n";
  cout << "--------------------------\n";  
  cout << "PtBin associated maximum edge = "<<fPtAssocLimit<<"\n";
  cout << "--------------------------\n";    
  cout << "Minimum D-meson pT = "<<fMinDPt<<"\n";
  cout << "--------------------------\n";  
  cout << "TTree filling = "<<fFillTrees<<"\n";
  cout << "--------------------------\n";  
}

