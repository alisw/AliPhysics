/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: ALICE Offline.                                                 *
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

//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    // 
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//        Copied from NetParticle Classes
//        Origin: Authors: Jochen Thaeder <jochen@thaeder.de>
//                         Michael Weber <m.weber@cern.ch>
//=========================================================================//

#include "TMath.h"
#include "TAxis.h"
#include "TSystem.h" 
#include "TFile.h" 
#include "TPRegexp.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliTracker.h"

#include "AliEbyEPidRatioHelper.h"

using namespace std;


ClassImp(AliEbyEPidRatioHelper)
//________________________________________________________________________
AliEbyEPidRatioHelper::AliEbyEPidRatioHelper() :
  fModeDistCreation(0),

  fInputEventHandler(NULL),
  fPIDResponse(NULL),
  fESD(NULL),
  fESDTrackCuts(NULL),
  fAOD(NULL),
  fAODtrackCutBit(1024),
  fIsMC(kFALSE),
  fMCEvent(NULL),
  fStack(NULL),

  fCentralityBin(-1),
  fCentralityPercentile(-1.),

  fCentralityBinMax(11),
  fVertexZMax(10.),
  fRapidityMax(0.9),
  fPhiMin(0.),
  fPhiMax(TMath::TwoPi()),
  fMinTrackLengthMC(70.),
  fNSigmaMaxCdd(3.),
  fNSigmaMaxCzz(3.),

  fPIDStrategy(0),
  fNSigmaMaxITS(4.),
  fNSigmaMaxTPC(4.),
  fNSigmaMaxTPClow(3.),
  fNSigmaMaxTOF(4.),
  fMinPtForTOFRequired(0.69),
  fMaxPtForTPClow(0.69),

  fHEventStat0(NULL),
  fHEventStat1(NULL),
  fHEventStatMax(6),

  fHTriggerStat(NULL),
  fNTriggers(5),

  fHCentralityStat(NULL),
  fHCentralityPer(NULL),
  fHCentralityPerAll(NULL),
  fNCentralityBins(11),
  
  fSubSamples(25),
  fRandom(NULL),
  fSubSampleIdx(1), 

  fIsRatio(kFALSE), 
  fIsPtBin(kFALSE), fIsDetectorWise(kFALSE) {
  // Constructor   
  
  AliLog::SetClassDebugLevel("AliEbyEPidRatioHelper",10);
}

const Float_t AliEbyEPidRatioHelper::fgkfHistBinWitdthRap = 0.1;
const Float_t AliEbyEPidRatioHelper::fgkfHistBinWitdthPt  = 0.3; // 0.08 // 300 MeV  // was 80 MeV

const Float_t AliEbyEPidRatioHelper::fgkfHistRangeCent[]  = {-0.5, 10.5};
const Int_t   AliEbyEPidRatioHelper::fgkfHistNBinsCent    = 11 ;

const Float_t AliEbyEPidRatioHelper::fgkfHistRangeEta[]   = {-0.9, 0.9};
const Int_t   AliEbyEPidRatioHelper::fgkfHistNBinsEta     = Int_t((AliEbyEPidRatioHelper::fgkfHistRangeEta[1] -
								   AliEbyEPidRatioHelper::fgkfHistRangeEta[0]) / 
								  AliEbyEPidRatioHelper::fgkfHistBinWitdthRap) +1;

const Float_t AliEbyEPidRatioHelper::fgkfHistRangeRap[]   = {-0.8, 0.8};
const Int_t   AliEbyEPidRatioHelper::fgkfHistNBinsRap     = Int_t((AliEbyEPidRatioHelper::fgkfHistRangeRap[1] - AliEbyEPidRatioHelper::fgkfHistRangeRap[0]) / AliEbyEPidRatioHelper::fgkfHistBinWitdthRap) +1;

const Float_t AliEbyEPidRatioHelper::fgkfHistRangePhi[]   = {0.0, static_cast<Float_t>(TMath::TwoPi())};
const Int_t   AliEbyEPidRatioHelper::fgkfHistNBinsPhi     = 42;

const Float_t AliEbyEPidRatioHelper::fgkfHistRangePt[]    = {0.2, 2.9}; // {0.2, 5.}; // was {0.3, 2.22}
const Int_t   AliEbyEPidRatioHelper::fgkfHistNBinsPt      = Int_t((AliEbyEPidRatioHelper::fgkfHistRangePt[1] - AliEbyEPidRatioHelper::fgkfHistRangePt[0]) / AliEbyEPidRatioHelper::fgkfHistBinWitdthPt); 

const Float_t AliEbyEPidRatioHelper::fgkfHistRangeSign[]  = {-1.5, 1.5};
const Int_t   AliEbyEPidRatioHelper::fgkfHistNBinsSign    =  3;

const Char_t* AliEbyEPidRatioHelper::fgkEventNames[]         = {"All", "IsTriggered", "HasVertex", "Vz<Vz_{Max}", "Centrality [0,100]%"};
const Char_t* AliEbyEPidRatioHelper::fgkCentralityMaxNames[] = {"5", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"};
const Char_t* AliEbyEPidRatioHelper::fgkTriggerNames[]       = {"kMB", "kCentral", "kSemiCentral", "kEMCEJE", "kEMCEGA" }; 
const Char_t* AliEbyEPidRatioHelper::fgkCentralityNames[]    = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%","50-60%", "60-70%", "70-80%", "80-90%", "90-100%"};

const Char_t* AliEbyEPidRatioHelper::fgkPidName[4]      = {"Nch","Npi","Nka","Npr"};
const Char_t* AliEbyEPidRatioHelper::fgkPidShLatex[4]      = {"N","#pi","K","p"};
const Char_t* AliEbyEPidRatioHelper::fgkPidLatex[4][2]  = {{"N_{-}","N_{+}"}, {"N_{#pi^{-}}","N_{#pi^{+}}"},{"N_{K^{-}}","N_{K^{+}}"}, {"N_{#bar{p}}","N_{p}"}};
const Char_t* AliEbyEPidRatioHelper::fgkPidTitles[4][2] = {{"Negative","Positive"},{"Anti-Pions","Pions"},{"Anti-Kaons","Kaons"}, {"Anti-Protons","Protons"}};


const Char_t* AliEbyEPidRatioHelper::fgkNetHistName[4]      = {"","Plus","Minus","Net"};
const Char_t* AliEbyEPidRatioHelper::fgkNetHistLatex[4]      = {"+ + +","+","-","+ - -"};
const Int_t AliEbyEPidRatioHelper::fgkfNetHistBin[4][4]  = {{3000,2400,1600,1200}, {1500,1200,800,600},{1500,1200,800,600},{600,600,600,600}};


//________________________________________________________________________
AliEbyEPidRatioHelper::~AliEbyEPidRatioHelper() {
  // Destructor

  if (fRandom)
    delete fRandom;

  return;
}

//________________________________________________________________________
void AliEbyEPidRatioHelper::SetPhiRange(Float_t f1, Float_t f2) {
  // -- Set phi range and adopt to phi-histogram
  
  fPhiMin = f1; 
  fPhiMax = (f1 < f2) ? f2 : f2+TMath::TwoPi();

  Float_t phiMin = fPhiMin;
  Float_t phiMax = fPhiMax;
  
  // -- Update Ranges
  Float_t binWidth = (AliEbyEPidRatioHelper::fgkfHistRangePhi[1] - AliEbyEPidRatioHelper::fgkfHistRangePhi[0]) / 
    Float_t(AliEbyEPidRatioHelper::fgkfHistNBinsPhi);

  Float_t lowEdge  = AliEbyEPidRatioHelper::fgkfHistRangePhi[0] - binWidth;
  Float_t highEdge = AliEbyEPidRatioHelper::fgkfHistRangePhi[0];

  for (Int_t ii = 1; ii <= AliEbyEPidRatioHelper::fgkfHistNBinsPhi; ++ii) {
    lowEdge += binWidth;
    highEdge += binWidth;

    if (phiMin >= lowEdge && phiMin < highEdge ) 
      phiMin = lowEdge;
    if (phiMax > lowEdge && phiMax <= highEdge ) 
      phiMax = highEdge;
  }
  
  printf(">>>> Update Phi Range : [%f,%f] -> [%f,%f]\n", fPhiMin, fPhiMax, phiMin, phiMax);
  fPhiMin = phiMin;
  fPhiMax = phiMax;
}



//________________________________________________________________________
Int_t AliEbyEPidRatioHelper::Initialize(AliESDtrackCuts *cuts, Bool_t isMC, Bool_t isRatio, Bool_t isPtBin, Bool_t isDetWise, Int_t trackCutBit, Int_t modeDistCreation) {
  // -- Initialize helper

  Int_t iResult = 0;

  // -- ESD track cuts
  fESDTrackCuts     = cuts;

  
  // -- Is MC
  fIsMC             = isMC;
  fIsRatio          = isRatio;
  fIsPtBin          = isPtBin;
  fIsDetectorWise   = isDetWise;

  

  // -- AOD track filter bit
  fAODtrackCutBit   = trackCutBit;
    
  // -- mode Distribution creation
  fModeDistCreation = modeDistCreation;

  // -- Setup event cut statistics 
  InitializeEventStats();

  // -- Setup trigger statistics 
  InitializeTriggerStats();

  // -- Setup centrality statistics 
  InitializeCentralityStats();

  // -- PRINT PID Strategy
  //    0 :   TPC(TPClow+TPCHigh)
  //    1 :   ITS
  //    2 :   TOF
  //    3 :   ITS+TPC(TPClow+TPCHigh)
  //    4 :   TPC(TPClow+TPCHigh)+TOF
  //    5 :   TPC(TPClow+TPCHigh)+TOF for pT >= fMinPtForTOFRequired TOF is required, below, only used if there
  //    6 :   TPC(TPClow+TPCHigh)+ITS+TOF with TOF only for those tracks which have TOF information
  //    7 :   TPC(TPClow+TPCHigh)+ITS+TOF for pT >= fMinPtForTOFRequired TOF is required, below, only used if there
  //    8 :   TPC(TPClow+TPCHigh)+ITS+TOF 
  printf(">>>>  PID STRATEGY: %d || sigmaMax: ITS %.2f TPC %.2f TOF %.2f \n", fPIDStrategy, fNSigmaMaxITS, fNSigmaMaxTPC, fNSigmaMaxTOF);            
                      
  // -- Initialize random number generator
  fRandom = new TRandom3();
  fRandom->SetSeed();
                      
  return iResult;
}

//________________________________________________________________________
Int_t AliEbyEPidRatioHelper::SetupEvent(AliESDInputHandler *esdHandler, AliAODInputHandler *aodHandler, AliMCEvent *mcEvent) {
  // -- Setup Event
  
  // -- Get ESD objects
  if(esdHandler){
    fInputEventHandler = static_cast<AliInputEventHandler*>(esdHandler);
    fESD               = dynamic_cast<AliESDEvent*>(fInputEventHandler->GetEvent());
    if (!fESD) {
      AliError("ESD event handler not available");
      return -1;
    }
  }

  // -- Get AOD objects
  else if(aodHandler){
    fInputEventHandler = static_cast<AliInputEventHandler*>(aodHandler);
    fAOD               = dynamic_cast<AliAODEvent*>(fInputEventHandler->GetEvent());
    if (!fAOD) {
      AliError("AOD event handler not available");
      return -1;
    }
  }

  // -- Get Common objects
  fPIDResponse = fInputEventHandler->GetPIDResponse();

  // -- Get MC objects
  fMCEvent     = mcEvent;
  if (fMCEvent)
    fStack     = fMCEvent->Stack();

  // -- Get event centrality
  // >  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90|90-100 --> 11 bins
  // >   0   1     2     3     4     5     6     7     8     9     10

  AliCentrality *centrality = NULL;

  if(esdHandler)
    centrality = fESD->GetCentrality();
  else if(aodHandler)
    centrality = ((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP();

  if (!centrality) {
    AliError("Centrality not available");
    return -1;
  }

  
  // Int_t a = centrality->GetCentralityClass5("V0M");
  // if (a < 0 || a >= 20 ) fCentralityBin = -2;
  // else if (a <= 1) fCentralityBin = a;
  // else fCentralityBin = 1 + centrality->GetCentralityClass10("V0M");
  
  
  
  Int_t centBin = centrality->GetCentralityClass10("V0M");
  if (centBin == 0) { fCentralityBin = centrality->GetCentralityClass5("V0M"); }
  else if (centBin == 11 || centBin == -1.)           { fCentralityBin = -1; }
  else if (centBin > 0 && centBin < fNCentralityBins) { fCentralityBin = centBin + 1; }
  else {  fCentralityBin = -2; }
  
  
  if (fCentralityBin >= fCentralityBinMax)
    fCentralityBin = -2;
    
  

  fCentralityPercentile = centrality->GetCentralityPercentile("V0M");
  
  fSubSampleIdx = fRandom->Integer(fSubSamples);

  return 0;
}
//________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::IsEventTriggered() {
  // -- Check if Event is triggered and fill Trigger Histogram
  
  Bool_t *aTriggerFired = new Bool_t[fNTriggers];
  for (Int_t ii = 0; ii < fNTriggers; ++ii)
    aTriggerFired[ii] = kFALSE;

  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kMB))          aTriggerFired[0] = kTRUE;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kCentral))     aTriggerFired[1] = kTRUE;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kSemiCentral)) aTriggerFired[2] = kTRUE;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kEMCEJE))      aTriggerFired[3] = kTRUE;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kEMCEGA))      aTriggerFired[4] = kTRUE;

  Bool_t isTriggered = kFALSE;

  for (Int_t ii=0; ii<fNTriggers; ++ii) {
    if(aTriggerFired[ii]) {
      isTriggered = kTRUE;
      fHTriggerStat->Fill(ii);
    }
  }
  
  delete[] aTriggerFired;

  return isTriggered;
}

//________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::IsEventRejected() {
  // -- Evaluate event statistics histograms
    
  Int_t *aEventCuts = new Int_t[fHEventStatMax];
  // set aEventCuts[ii] to 1 in case of reject
  
  for (Int_t ii=0;ii<fHEventStatMax; ++ii)
    aEventCuts[ii] = 0;

  Int_t iCut = 0;

  // -- 0 - Before Physics Selection   
  aEventCuts[iCut] = 0;

  // -- 1 - No Trigger fired
  ++iCut;
  if (!IsEventTriggered())
    aEventCuts[iCut] = 1;

  // -- 2 - No Vertex 
  ++iCut;
  const AliESDVertex* vtxESD = NULL;
  const AliAODVertex* vtxAOD = NULL;
  if (fESD){
    vtxESD = fESD->GetPrimaryVertexTracks();
    if (!vtxESD)
      aEventCuts[iCut] = 1;
  }
  else if (fAOD){
    vtxAOD = fAOD->GetPrimaryVertex();
    if (!vtxAOD)
      aEventCuts[iCut] = 1;
  }

  // -- 3 - Vertex z outside cut window
  ++iCut;
  if (vtxESD){
    if(TMath::Abs(vtxESD->GetZ()) > fVertexZMax) 
      aEventCuts[iCut] = 1;
  }
  else if(vtxAOD){
    if(TMath::Abs(vtxAOD->GetZ()) > fVertexZMax) 
      aEventCuts[iCut] = 1;
  }
  else
    aEventCuts[iCut] = 1;

  // -- 4 - Centrality = -1  (no centrality or not hadronic)
  ++iCut;
  if(fCentralityBin == -1.) 
    aEventCuts[iCut] = 1;

  // -- 5 - Centrality < fCentralityMax
  ++iCut;
  if(fCentralityBin == -2.) 
    aEventCuts[iCut] = 1;

  // -- Fill statistics / reject event
  Bool_t isRejected = FillEventStats(aEventCuts);

  // -- Cleanup 
  delete[] aEventCuts;

  //cout << isRejected << endl;

  return isRejected;
}

//________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::IsParticleAcceptedBasicCharged(AliVParticle *particle, Int_t idxMC) {
  // -- Check if MC particle is accepted for basic parameters
  
  if (!particle) 
    return kFALSE;

  // -- check if charged
  if (particle->Charge() == 0.0) 
    return kFALSE;
  
  // -- check if physical primary - ESD
  if (fESD) {
    if(!fStack->IsPhysicalPrimary(idxMC)) 
      return kFALSE;
  }
  // -- check if physical primary - AOD
  else {
    if(!(static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary()) 
      return kFALSE;
  }
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::IsParticleAcceptedBasicNeutral(AliVParticle *particle, Int_t idxMC) {
  // -- Check if MC particle is accepted for basic parameters
  
  if (!particle) 
    return kFALSE;

  // -- check if charged
  if (particle->Charge() != 0.0) 
    return kFALSE;
  
  // -- check if physical primary - ESD
  if (fESD) {
    if(!fStack->IsPhysicalPrimary(idxMC)) 
      return kFALSE;
  }
  // -- check if physical primary - AOD
  else {
    if(!(static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary()) 
      return kFALSE;
  }
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::IsParticleAcceptedRapidity(AliVParticle *particle, Double_t &yP, Int_t gCurPid) {
  
 if (gCurPid == 0) {
   yP = particle->Eta();
   return kTRUE;
 }
  
 Double_t mP = AliPID::ParticleMass(AliPID::kPion);
 if(gCurPid == 1) mP = AliPID::ParticleMass(AliPID::kPion);
 else if(gCurPid == 2) mP = AliPID::ParticleMass(AliPID::kKaon);
 else if(gCurPid == 3) mP = AliPID::ParticleMass(AliPID::kProton);
 
 // -- Calculate rapidities and kinematics
  Double_t p  = particle->P();
  Double_t pz = particle->Pz();

  Double_t eP = TMath::Sqrt(p*p + mP*mP);
  yP          = 0.5 * TMath::Log((eP + pz) / (eP - pz));  

  // -- Check Rapidity window
  if (TMath::Abs(yP) > fRapidityMax)
    return kFALSE;
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::IsParticleAcceptedPhi(AliVParticle *particle) {
  // -- Check if particle is accepted
  // > in phi
  // > return 0 if not accepted
  
  if (particle->Phi() > fPhiMin && particle->Phi() <= fPhiMax)
    return kTRUE;
  else if (particle->Phi() < fPhiMin && (particle->Phi() + TMath::TwoPi()) <= fPhiMax)
    return kTRUE;
  else
    return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::IsParticleFindable(Int_t label) {
  // -- Check if MC particle is findable tracks

  AliMCParticle *mcParticle = static_cast<AliMCParticle*>(fMCEvent->GetTrack(label));
  if(!mcParticle) 
    return kFALSE;
  
  Int_t counter; 
  Float_t tpcTrackLength = mcParticle->GetTPCTrackLength(AliTracker::GetBz(), 0.05, counter, 3.0); 

  return (tpcTrackLength > fMinTrackLengthMC);    
}
//________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::IsTrackAcceptedBasicCharged(AliVTrack* track) {
  // -- Check if track is accepted 
  // > for basic parameters

  if (!track)
    return kFALSE;
  
  if (track->Charge() == 0) 
    return kFALSE;
  
  return kTRUE;
} 
 
//________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::IsTrackAcceptedRapidity(AliVTrack *track, Double_t &yP, Int_t gCurPid) {
   if (gCurPid == 0) {
    yP = track->Eta();
    return kTRUE;
  }
  
  Double_t mP = AliPID::ParticleMass(AliPID::kPion);
  if(gCurPid == 1) mP = AliPID::ParticleMass(AliPID::kPion);
  else if(gCurPid == 2) mP = AliPID::ParticleMass(AliPID::kKaon);
  else if(gCurPid == 3) mP = AliPID::ParticleMass(AliPID::kProton);

  // -- Calculate rapidities and kinematics
  Double_t pvec[3];
  track->GetPxPyPz(pvec);

  Double_t p  = track->P();
  Double_t eP = TMath::Sqrt(p*p + mP*mP);
           yP = 0.5 * TMath::Log((eP + pvec[2]) / (eP - pvec[2]));
  
  // -- Check Rapidity window
  if (TMath::Abs(yP) > fRapidityMax)
    return kFALSE;
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::IsTrackAcceptedDCA(AliVTrack *vTrack) {
   Bool_t isAccepted = kTRUE;

  if (!fESD)
    return isAccepted;

  AliESDtrack* track = dynamic_cast<AliESDtrack*>(vTrack);
  if (!track)
    return kFALSE;
  
  // -- Get nHits SPD
  if (track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1)) {

    // -- Get DCA nSigmas
    Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
    track->GetImpactParameters(dca,cov);

    Float_t nSigmaCdd = (cov[0] != 0.) ? dca[0]/TMath::Sqrt(cov[0]) : -9.99; 
    Float_t nSigmaCzz = (cov[2] != 0.) ? dca[1]/TMath::Sqrt(cov[2]) : -9.99; 
    
    if (fNSigmaMaxCdd != 0.) {
      if (TMath::Abs(nSigmaCdd) > fNSigmaMaxCdd)
	isAccepted = kFALSE;
    }

    if (fNSigmaMaxCzz != 0.) {
      if (TMath::Abs(nSigmaCzz) > fNSigmaMaxCzz)
	isAccepted = kFALSE;
    }
  }

  return isAccepted;
}

//________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::IsTrackAcceptedPID(AliVTrack *track, Double_t* pid, AliPID::EParticleType gCurPid) {

  Bool_t isAcceptedITS    = kFALSE;
  Bool_t isAcceptedTPC    = kFALSE;
  Bool_t isAcceptedTPClow = kFALSE;
  Bool_t isAcceptedTOF    = kFALSE;
  Bool_t isAccepted       = kFALSE;

  // -- In case not PID is used
  if (gCurPid == AliPID::kUnknown) {
    pid[0] = 10.;
    pid[1] = 10.;
    pid[2] = 10.;
    return kTRUE;
  }
  
  if (fPIDResponse->NumberOfSigmas(AliPIDResponse::kITS, track, gCurPid, pid[0]) == AliPIDResponse::kDetPidOk) {
    if (TMath::Abs(pid[0]) < fNSigmaMaxITS) 
      isAcceptedITS = kTRUE;
  }

  if (fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, gCurPid, pid[1]) == AliPIDResponse::kDetPidOk) {
    if (TMath::Abs(pid[1]) < fNSigmaMaxTPC) 
      isAcceptedTPC = kTRUE;
    if (TMath::Abs(pid[1]) < fNSigmaMaxTPClow) 
      isAcceptedTPClow = kTRUE;
    if (track->Pt() < fMaxPtForTPClow)
      isAcceptedTPC = isAcceptedTPClow;
  }

  Bool_t hasPIDTOF = kFALSE;
  if (fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, gCurPid, pid[2]) == AliPIDResponse::kDetPidOk) {
    hasPIDTOF = kTRUE;
    if (TMath::Abs(pid[2]) < fNSigmaMaxTOF) 
      isAcceptedTOF = kTRUE;
  }
  // -- Check TOF missmatch for MC
  
  //if (ESD)
  if (fIsMC && isAcceptedTOF) {
    Int_t tofLabel[3];                                                                                                                                        
    //    AliESDtrack* track = dynamic_cast<AliESDtrack*>(vTrack);
    // TODO add code for AOD 

    (dynamic_cast<AliESDtrack*>(track))->GetTOFLabel(tofLabel);
   
    Bool_t hasMatchTOF = kTRUE;
    if (TMath::Abs(track->GetLabel()) != TMath::Abs(tofLabel[0]) || tofLabel[1] > 0) 
      hasMatchTOF = kFALSE;

    TParticle *matchedTrack = fStack->Particle(TMath::Abs(tofLabel[0]));
    if (TMath::Abs(matchedTrack->GetFirstMother()) == TMath::Abs(track->GetLabel())) 
      hasMatchTOF = kTRUE;

    isAcceptedTOF = hasMatchTOF;
  }

  //    0 :   TPC(TPClow+TPCHigh)
  //    1 :   ITS
  //    2 :   TOF
  //    3 :   ITS+TPC(TPClow+TPCHigh)
  //    4 :   TPC(TPClow+TPCHigh)+TOF
  //    5 :   TPC(TPClow+TPCHigh)+TOF for pT >= fMinPtForTOFRequired TOF is required, below, only used if there
  //    6 :   TPC(TPClow+TPCHigh)+ITS+TOF with TOF only for those tracks which have TOF information
  //    7 :   TPC(TPClow+TPCHigh)+ITS+TOF for pT >= fMinPtForTOFRequired TOF is required, below, only used if there
  //    8 :   TPC(TPClow+TPCHigh)+ITS+TOF 
  if (fPIDStrategy == 0) {             //  TPC PID
    isAccepted = isAcceptedTPC;
  }
  else if (fPIDStrategy == 1) {        //  ITS PID
    isAccepted = isAcceptedITS;
  }
  else if (fPIDStrategy == 2) {        //  TOF PID
    isAccepted = isAcceptedTOF;
  }
  else if (fPIDStrategy == 3) {        //  TPC+ITS PID
    isAccepted = isAcceptedTPC && isAcceptedITS;
  }
  else if (fPIDStrategy == 4) {        //  TPC+TOF PID
    isAccepted = isAcceptedTPC && isAcceptedTOF;
  }
  else if (fPIDStrategy == 5) {        //  TPC+TOF PID -- for pT >= fMinPtForTOFRequired TOF is required
    if (!hasPIDTOF && track->Pt() < fMinPtForTOFRequired) 
      isAcceptedTOF = kTRUE;
    isAccepted = isAcceptedTPC && isAcceptedTOF;
  }
  else if (fPIDStrategy == 6) {        //  ITS+TPC+TOF PID -- TOF only for those tracks which have TOF information
    isAccepted = isAcceptedTPC && isAcceptedITS;
    if (hasPIDTOF)
      isAccepted = isAccepted && isAcceptedTOF;
  }
  else if (fPIDStrategy == 7) {        //  ITS+TPC+TOF PID -- for pT >= fMinPtForTOFRequired TOF is required
    if (!hasPIDTOF && track->Pt() < fMinPtForTOFRequired)
      isAcceptedTOF = kTRUE;
    isAccepted = isAcceptedITS && isAcceptedTPC && isAcceptedTOF;
  }
  else if (fPIDStrategy == 8) {        //  ITS+TPC+TOF PID
    isAccepted = isAcceptedITS && isAcceptedTPC && isAcceptedTOF;
  }

  return isAccepted;
}

//________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::IsTrackAcceptedPhi(AliVTrack *track) {
  // -- Check if track is accepted
  // > in phi
  // > return 0 if not accepted
  
  if (track->Phi() > fPhiMin && track->Phi() <= fPhiMax)
    return kTRUE;
  else if (track->Phi() < fPhiMin && (track->Phi() + TMath::TwoPi()) <= fPhiMax)
    return kTRUE;
  else
    return kFALSE;
}

//________________________________________________________________________
void AliEbyEPidRatioHelper::BinLogAxis(const THnBase *hn, Int_t axisNumber, AliESDtrackCuts* cuts) {
    AliESDtrackCuts* esdTrackCuts = (cuts) ? cuts : fESDTrackCuts;

  // -- Make logarithmic binning 
  TAxis *axis = hn->GetAxis(axisNumber);
  Int_t  nBins = axis->GetNbins();

  Double_t from  = axis->GetXmin();
  Double_t to    = axis->GetXmax();
  Double_t *newBins = new Double_t[nBins + 1];
   
  newBins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./nBins);
  
  for (int ii = 1; ii <= nBins; ii++)
   newBins[ii] = factor * newBins[ii-1];
  
  axis->Set(nBins, newBins);

  delete [] newBins;

  // -- Update Ranges
  // ------------------
  Float_t ptRange[2];
  Float_t oldPtRange[2];
  esdTrackCuts->GetPtRange(ptRange[0],ptRange[1]);
  esdTrackCuts->GetPtRange(oldPtRange[0],oldPtRange[1]);

  Float_t minPtForTOFRequired = fMinPtForTOFRequired;
  Float_t maxPtForTPClow      = fMaxPtForTPClow;

  // -- Update minPt Cut
  Int_t bin = axis->FindBin(ptRange[0]+10e-7);
  ptRange[0] = axis->GetBinLowEdge(bin); 

  // -- Update maxPt Cut
  bin = axis->FindBin(ptRange[1]-10e-7);
  ptRange[1] = axis->GetBinUpEdge(bin); 

  if (ptRange[0] != oldPtRange[0] || ptRange[1] != oldPtRange[1]) {
    printf(">>>> Update Pt Range : [%f,%f] -> [%f,%f]\n", oldPtRange[0], oldPtRange[1], ptRange[0], ptRange[1]);
    esdTrackCuts->SetPtRange(ptRange[0],ptRange[1]);
  }

  // -- Update MinPtForTOFRequired
  bin = axis->FindBin(minPtForTOFRequired-10e-7);
  minPtForTOFRequired = axis->GetBinUpEdge(bin); 

  if (minPtForTOFRequired != fMinPtForTOFRequired) {
    printf(">>>> Update Min Pt for TOF : %f -> %f\n", fMinPtForTOFRequired, minPtForTOFRequired);
    fMinPtForTOFRequired = minPtForTOFRequired;
  }

  // -- Update MaxPtForTPClow
  bin = axis->FindBin(maxPtForTPClow-10e-7);
  maxPtForTPClow = axis->GetBinUpEdge(bin); 

  if (maxPtForTPClow != fMaxPtForTPClow) {
    printf(">>>> Update Max Pt for TPC Low : %f -> %f\n", fMaxPtForTPClow, maxPtForTPClow);
    fMaxPtForTPClow = maxPtForTPClow;
  }
}

//________________________________________________________________________
void AliEbyEPidRatioHelper::InitializeEventStats() {
  // -- Initialize event statistics histograms

  fHEventStat0 = new TH1F("hEventStat0","Event cut statistics 0;Event Cuts;Events", fHEventStatMax,-0.5,fHEventStatMax-0.5);
  fHEventStat1 = new TH1F("hEventStat1","Event cut statistics 1;Event Cuts;Events", fHEventStatMax,-0.5,fHEventStatMax-0.5);

  for ( Int_t ii=0; ii < fHEventStatMax-1; ii++ ) {
    fHEventStat0->GetXaxis()->SetBinLabel(ii+1, AliEbyEPidRatioHelper::fgkEventNames[ii]);
    fHEventStat1->GetXaxis()->SetBinLabel(ii+1, AliEbyEPidRatioHelper::fgkEventNames[ii]);
  }

  fHEventStat0->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", AliEbyEPidRatioHelper::fgkCentralityMaxNames[fCentralityBinMax-1]));
  fHEventStat1->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", AliEbyEPidRatioHelper::fgkCentralityMaxNames[fCentralityBinMax-1]));
}

//________________________________________________________________________
void AliEbyEPidRatioHelper::InitializeTriggerStats() {
  // -- Initialize trigger statistics histograms

  fHTriggerStat = new TH1F("hTriggerStat","Trigger statistics;Trigger;Events", fNTriggers,-0.5,fNTriggers-0.5);

  for ( Int_t ii=0; ii < fNTriggers; ii++ )
    fHTriggerStat->GetXaxis()->SetBinLabel(ii+1, AliEbyEPidRatioHelper::fgkTriggerNames[ii]);
}

//________________________________________________________________________
void AliEbyEPidRatioHelper::InitializeCentralityStats() {
  // -- Initialize trigger statistics histograms

  fHCentralityStat = new TH1F("hCentralityStat","Centrality statistics;Centrality Bins;Events", 
			      fNCentralityBins,-0.5,fNCentralityBins-0.5);
  
  for ( Int_t ii=0; ii < fNCentralityBins; ii++ )
    fHCentralityStat->GetXaxis()->SetBinLabel(ii+1, AliEbyEPidRatioHelper::fgkCentralityNames[ii]);
  
  fHCentralityPer = new TH1F("hCentralityPercentileAccepted","Centrality Percentile statistics;Centrality Bins;Events", 
			     100,-0.5,99.5);

  fHCentralityPerAll = new TH1F("hCentralityPercentileAll","Centrality Percentile statistics;Centrality Bins;Events", 
			     100,-0.5,99.5);


}
//________________________________________________________________________
Bool_t AliEbyEPidRatioHelper::FillEventStats(Int_t *aEventCuts) {
  // -- Fill event / centrality statistics 

  Bool_t isRejected = kFALSE;


  // -- Fill event statistics
  for (Int_t idx = 0; idx < fHEventStatMax ; ++idx) {

    if (aEventCuts[idx])
      isRejected = kTRUE;
    else
      fHEventStat0->Fill(idx);
  }
  
  for (Int_t idx = 0; idx < fHEventStatMax; ++idx) {
    if (aEventCuts[idx])
      break;
    fHEventStat1->Fill(idx);
  }

  // -- Fill centrality statistics of accepted events
  if (!isRejected) {
    fHCentralityStat->Fill(fCentralityBin);
    fHCentralityPer->Fill(fCentralityPercentile);
  }

  fHCentralityPerAll->Fill(fCentralityPercentile);



  return isRejected;
}
