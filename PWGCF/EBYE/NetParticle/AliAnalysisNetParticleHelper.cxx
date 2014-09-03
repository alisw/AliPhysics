//-*- Mode: C++ -*-

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

#include "AliAnalysisNetParticleHelper.h"

using namespace std;

/**
 * Class for NetParticle Distributions
 * -- Helper class for net particle istributions
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

ClassImp(AliAnalysisNetParticleHelper)

/*
 * ---------------------------------------------------------------------------------
 *                            particle names 
 * ---------------------------------------------------------------------------------
 */

  // MW make fgk ... static const
  const Char_t* aPartNames[AliPID::kSPECIES][2]             = {
    {"ele",     "posi"},
    {"mubar",   "mu"},
    {"pibar",   "pi"},
    {"kbar",    "k"},
    {"pbar",    "p"}
  };

  const Char_t* aPartTitles[AliPID::kSPECIES][2]            = {
    {"Electron",    "Positron"},
    {"Anti-Muon",   "Muon"},
    {"Anti-Pion",   "Proton"},
    {"Anti-Kaon",   "Kaon"},
    {"Anti-Proton", "Proton"}
  };

  const Char_t* aPartTitlesLatex[AliPID::kSPECIES][2]        = {
    {"e^{-}",   "e^{+}" },
    {"#mu^{-}", "#mu^{+}" },
    {"#pi^{-}", "#pi^{+}" },
    {"K^{-}",   "K^{+}" },
    {"#bar{p}", "p"}
  };

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
AliAnalysisNetParticleHelper::AliAnalysisNetParticleHelper() :
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

  fCentralityBinMax(7),
  fVertexZMax(10.),
  fRapidityMax(0.5),
  fPhiMin(0.),
  fPhiMax(TMath::TwoPi()),
  fMinTrackLengthMC(70.),
  fNSigmaMaxCdd(3.),
  fNSigmaMaxCzz(3.),

  fParticleSpecies(AliPID::kProton),

  fUsePID(kTRUE),
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
  fNCentralityBins(10),

  fNSubSamples(20),
  fSubSampleIdx(0),
  fRandom(NULL) {
  // Constructor   
  
  AliLog::SetClassDebugLevel("AliAnalysisNetParticleHelper",10);
}

const Float_t AliAnalysisNetParticleHelper::fgkfHistBinWitdthRap = 0.075;
const Float_t AliAnalysisNetParticleHelper::fgkfHistBinWitdthPt  = 0.3; // 0.08 // 300 MeV  // was 80 MeV

const Float_t AliAnalysisNetParticleHelper::fgkfHistRangeCent[]  = {-0.5, 8.5};
const Int_t   AliAnalysisNetParticleHelper::fgkfHistNBinsCent    = 9 ;

const Float_t AliAnalysisNetParticleHelper::fgkfHistRangeEta[]   = {-0.9, 0.9};
const Int_t   AliAnalysisNetParticleHelper::fgkfHistNBinsEta     = Int_t((AliAnalysisNetParticleHelper::fgkfHistRangeEta[1] -
									  AliAnalysisNetParticleHelper::fgkfHistRangeEta[0]) / 
									 AliAnalysisNetParticleHelper::fgkfHistBinWitdthRap) +1;

const Float_t AliAnalysisNetParticleHelper::fgkfHistRangeRap[]   = {-0.5, 0.5};
const Int_t   AliAnalysisNetParticleHelper::fgkfHistNBinsRap     = Int_t((AliAnalysisNetParticleHelper::fgkfHistRangeRap[1] -
									  AliAnalysisNetParticleHelper::fgkfHistRangeRap[0]) / 
									 AliAnalysisNetParticleHelper::fgkfHistBinWitdthRap) +1;

const Float_t AliAnalysisNetParticleHelper::fgkfHistRangePhi[]   = {0.0, static_cast<Float_t>(TMath::TwoPi())};
const Int_t   AliAnalysisNetParticleHelper::fgkfHistNBinsPhi     = 42 ;

const Float_t AliAnalysisNetParticleHelper::fgkfHistRangePt[]    = {0.2, 2.9}; // {0.2, 5.}; // was {0.3, 2.22}
const Int_t   AliAnalysisNetParticleHelper::fgkfHistNBinsPt      = Int_t((AliAnalysisNetParticleHelper::fgkfHistRangePt[1] -
									  AliAnalysisNetParticleHelper::fgkfHistRangePt[0]) / 
									 AliAnalysisNetParticleHelper::fgkfHistBinWitdthPt); 

const Float_t AliAnalysisNetParticleHelper::fgkfHistRangeSign[]  = {-1.5, 1.5};
const Int_t   AliAnalysisNetParticleHelper::fgkfHistNBinsSign    =  3;

const Char_t* AliAnalysisNetParticleHelper::fgkEventNames[]         = {"All", "IsTriggered", "HasVertex", "Vz<Vz_{Max}", "Centrality [0,90]%"};
const Char_t* AliAnalysisNetParticleHelper::fgkCentralityMaxNames[] = {"5", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"};
const Char_t* AliAnalysisNetParticleHelper::fgkTriggerNames[]       = {"kMB", "kCentral", "kSemiCentral", "kEMCEJE", "kEMCEGA" }; 
const Char_t* AliAnalysisNetParticleHelper::fgkCentralityNames[]    = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", 
								       "50-60%", "60-70%", "70-80%", "80-90%", "90-100%"};

//________________________________________________________________________
AliAnalysisNetParticleHelper::~AliAnalysisNetParticleHelper() {
  // Destructor

  if (fRandom)
    delete fRandom;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                    Setter
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleHelper::SetPhiRange(Float_t f1, Float_t f2) {
  // -- Set phi range and adopt to phi-histogram
  
  fPhiMin = f1; 
  fPhiMax = (f1 < f2) ? f2 : f2+TMath::TwoPi();

  Float_t phiMin = fPhiMin;
  Float_t phiMax = fPhiMax;
  
  // -- Update Ranges
  Float_t binWidth = (AliAnalysisNetParticleHelper::fgkfHistRangePhi[1] - AliAnalysisNetParticleHelper::fgkfHistRangePhi[0]) / 
    Float_t(AliAnalysisNetParticleHelper::fgkfHistNBinsPhi);

  Float_t lowEdge  = AliAnalysisNetParticleHelper::fgkfHistRangePhi[0] - binWidth;
  Float_t highEdge = AliAnalysisNetParticleHelper::fgkfHistRangePhi[0];

  for (Int_t ii = 1; ii <= AliAnalysisNetParticleHelper::fgkfHistNBinsPhi; ++ii) {
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
void AliAnalysisNetParticleHelper::SetParticleSpecies(AliPID::EParticleType pid) {
  // -- Set particle species (ID, Name, Title, Title LATEX)

  if ( Int_t(pid) < 0 || Int_t(pid) >= AliPID::kSPECIES) {  
    AliWarning("Particle ID not in AliPID::kSPECIES --> Set to protons");
    pid = AliPID::kProton;
  }  
  
  fParticleSpecies = pid;
  
  for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
    fPartName[idxPart]       = aPartNames[fParticleSpecies][idxPart];
    fPartTitle[idxPart]      = aPartTitles[fParticleSpecies][idxPart];
    fPartTitleLatex[idxPart] = aPartTitlesLatex[fParticleSpecies][idxPart];
  }
}

//________________________________________________________________________
void AliAnalysisNetParticleHelper::SetUsePID(Bool_t usePID) {
  // -- Set usage of PID
  //    > if turn off, set charge types (ID, Name, Title, Title LATEX)
  
  fUsePID = usePID;
  
  if (!usePID) {
    fParticleSpecies   = AliPID::kUnknown;

    fPartName[0]       = "neg";
    fPartName[1]       = "pos";
    fPartTitle[0]      = "Negative";
    fPartTitle[1]      = "Positive";
    fPartTitleLatex[0] = "Negative";
    fPartTitleLatex[1] = "Positive";
  }
}

/*
 * ---------------------------------------------------------------------------------
 *                                    Getter
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
TString AliAnalysisNetParticleHelper::GetParticleName(Int_t idxPart) {
  // -- Get particle Name

  if( idxPart != 0 && idxPart != 1){
    AliWarning("Particle type not known --> Set to antiparticles");
    idxPart = 0;
  }

  return fPartName[idxPart];
}

//________________________________________________________________________
TString AliAnalysisNetParticleHelper::GetParticleTitle(Int_t idxPart) {
  // -- Get particle Title

  if( idxPart != 0 && idxPart != 1){
    AliWarning("Particle type not known --> Set to antiparticles");
    idxPart = 0;
  }

  return fPartTitle[idxPart];
}

//________________________________________________________________________
TString AliAnalysisNetParticleHelper::GetParticleTitleLatex(Int_t idxPart) {
  // -- Get particle Title LATEX

  if( idxPart != 0 && idxPart != 1){
    AliWarning("Particle type not known --> Set to antiparticles");
    idxPart = 0;
  }

  return fPartTitleLatex[idxPart];
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Int_t AliAnalysisNetParticleHelper::Initialize(AliESDtrackCuts *cuts, Bool_t isMC, Int_t trackCutBit, Int_t modeDistCreation) {
  // -- Initialize helper

  Int_t iResult = 0;

  // -- ESD track cuts
  fESDTrackCuts     = cuts;

  // -- Is MC
  fIsMC             = isMC;

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
  printf(">>>> USE PID %d || PID STRATEGY: %d || sigmaMax: ITS %.2f TPC %.2f TOF %.2f \n", fUsePID, fPIDStrategy, fNSigmaMaxITS, fNSigmaMaxTPC, fNSigmaMaxTOF);            
                      
  // -- Initialize random number generator
  fRandom = new TRandom3();
  fRandom->SetSeed();
                      
  return iResult;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleHelper::SetupEvent(AliESDInputHandler *esdHandler, AliAODInputHandler *aodHandler, AliMCEvent *mcEvent) {
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
  // >  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  // >   0   1     2     3     4     5     6     7     8     9

  AliCentrality *centrality = NULL;

  if(esdHandler)
    centrality = fESD->GetCentrality();
  else if(aodHandler)
    centrality = ((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP();

  if (!centrality) {
    AliError("Centrality not available");
    return -1;
  }

  Int_t centBin = centrality->GetCentralityClass10("V0M");
  if (centBin == 0)
    fCentralityBin = centrality->GetCentralityClass5("V0M");
  else if (centBin == 10 || centBin == -1.)
    fCentralityBin = -1;
  else if (centBin > 0 && centBin < fNCentralityBins)
    fCentralityBin = centBin + 1;
  else
    fCentralityBin = -2;

  // -- Stay within the max centrality bin
  if (fCentralityBin >= fCentralityBinMax)
    fCentralityBin = -2;

  fCentralityPercentile = centrality->GetCentralityPercentile("V0M");

  // -- Get current subsample idx
  fSubSampleIdx = fRandom->Integer(fNSubSamples);

  return 0;
}

/*
 * ---------------------------------------------------------------------------------
 *                         Event / Trigger Statistics
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Bool_t AliAnalysisNetParticleHelper::IsEventTriggered() {
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
Bool_t AliAnalysisNetParticleHelper::IsEventRejected() {
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

  return isRejected;
}

/*
 * ---------------------------------------------------------------------------------
 *                         Accept Particle Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Bool_t AliAnalysisNetParticleHelper::IsParticleAcceptedBasicCharged(AliVParticle *particle, Int_t idxMC) {
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
Bool_t AliAnalysisNetParticleHelper::IsParticleAcceptedBasicNeutral(AliVParticle *particle, Int_t idxMC) {
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
Bool_t AliAnalysisNetParticleHelper::IsParticleAcceptedRapidity(AliVParticle *particle, Double_t &yP) {
  // -- Check if particle is accepted
  // > in rapidity
  // > if no pid : return kTRUE, yP = eta
  // > return 0 if not accepted

  if (!fUsePID) {
    yP = particle->Eta();
    return kTRUE;
  }

  Double_t mP = AliPID::ParticleMass(fParticleSpecies);

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
Bool_t AliAnalysisNetParticleHelper::IsParticleAcceptedPhi(AliVParticle *particle) {
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
Bool_t AliAnalysisNetParticleHelper::IsParticleFindable(Int_t label) {
  // -- Check if MC particle is findable tracks

  AliMCParticle *mcParticle = static_cast<AliMCParticle*>(fMCEvent->GetTrack(label));
  if(!mcParticle) 
    return kFALSE;
  
  Int_t counter; 
  Float_t tpcTrackLength = mcParticle->GetTPCTrackLength(AliTracker::GetBz(), 0.05, counter, 3.0); 

  return (tpcTrackLength > fMinTrackLengthMC);    
}

/*
 * ---------------------------------------------------------------------------------
 *                            Accept Track Methods - public
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Bool_t AliAnalysisNetParticleHelper::IsTrackAcceptedBasicCharged(AliVTrack* track) {
  // -- Check if track is accepted 
  // > for basic parameters

  if (!track)
    return kFALSE;
  
  if (track->Charge() == 0) 
    return kFALSE;
  
  return kTRUE;
} 
 
//________________________________________________________________________
Bool_t AliAnalysisNetParticleHelper::IsTrackAcceptedRapidity(AliVTrack *track, Double_t &yP) {
  // -- Check if track is accepted
  // > in rapidity
  // > if no pid : return kTRUE
  // > return 0 if not accepted

  if (!fUsePID) {
    yP = track->Eta();
    return kTRUE;
  }
  
  Double_t mP = AliPID::ParticleMass(fParticleSpecies);

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
Bool_t AliAnalysisNetParticleHelper::IsTrackAcceptedDCA(AliVTrack *vTrack) {
  // -- Check if track is accepted - ONLY FOR ESDs so far 
  // > for DCA, if both SPD layers have hits
  // > For now only Implemented for ESDs

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
Bool_t AliAnalysisNetParticleHelper::IsTrackAcceptedPID(AliVTrack *track, Double_t* pid) {
  // -- Check if track is accepted 
  // > provides ITS, TPC and TOF nSigmas to the argument

  Bool_t isAcceptedITS    = kFALSE;
  Bool_t isAcceptedTPC    = kFALSE;
  Bool_t isAcceptedTPClow = kFALSE;
  Bool_t isAcceptedTOF    = kFALSE;
  Bool_t isAccepted       = kFALSE;

  // -- In case not PID is used
  if (!fUsePID) {
    pid[0] = 10.;
    pid[1] = 10.;
    pid[2] = 10.;
    return kTRUE;
  }
  
  // -- Get PID with ITS and check
  if (fPIDResponse->NumberOfSigmas(AliPIDResponse::kITS, track, fParticleSpecies, pid[0]) == AliPIDResponse::kDetPidOk) {
    if (TMath::Abs(pid[0]) < fNSigmaMaxITS) 
      isAcceptedITS = kTRUE;
  }

  // -- Get PID with TPC and check
  if (fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, fParticleSpecies, pid[1]) == AliPIDResponse::kDetPidOk) {
    if (TMath::Abs(pid[1]) < fNSigmaMaxTPC) 
      isAcceptedTPC = kTRUE;
    if (TMath::Abs(pid[1]) < fNSigmaMaxTPClow) 
      isAcceptedTPClow = kTRUE;
    if (track->Pt() < fMaxPtForTPClow)
      isAcceptedTPC = isAcceptedTPClow;
  }

  // -- Get PID with TOF and check
  Bool_t hasPIDTOF = kFALSE;
  if (fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, fParticleSpecies, pid[2]) == AliPIDResponse::kDetPidOk) {
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
Bool_t AliAnalysisNetParticleHelper::IsTrackAcceptedPhi(AliVTrack *track) {
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

/*
 * ---------------------------------------------------------------------------------
 *                         Helper Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleHelper::BinLogAxis(const THnBase *hn, Int_t axisNumber, AliESDtrackCuts* cuts) {
  // -- Method for the correct logarithmic binning of histograms
  // -- and update fMinPtForTOFRequired using the logarithmic scale
  
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
  Int_t bin = axis->FindBin(ptRange[0]+10e-6);
  ptRange[0] = axis->GetBinLowEdge(bin); 

  // -- Update maxPt Cut
  bin = axis->FindBin(ptRange[1]-10e-6);
  ptRange[1] = axis->GetBinUpEdge(bin); 

  if (ptRange[0] != oldPtRange[0] || ptRange[1] != oldPtRange[1]) {
    printf(">>>> Update Pt Range : [%f,%f] -> [%f,%f]\n", oldPtRange[0], oldPtRange[1], ptRange[0], ptRange[1]);
    esdTrackCuts->SetPtRange(ptRange[0],ptRange[1]);
  }

  // -- Update MinPtForTOFRequired
  bin = axis->FindBin(minPtForTOFRequired-10e-6);
  minPtForTOFRequired = axis->GetBinUpEdge(bin); 

  if (minPtForTOFRequired != fMinPtForTOFRequired) {
    printf(">>>> Update Min Pt for TOF : %f -> %f\n", fMinPtForTOFRequired, minPtForTOFRequired);
    fMinPtForTOFRequired = minPtForTOFRequired;
  }

  // -- Update MaxPtForTPClow
  bin = axis->FindBin(maxPtForTPClow-10e-6);
  maxPtForTPClow = axis->GetBinUpEdge(bin); 

  if (maxPtForTPClow != fMaxPtForTPClow) {
    printf(">>>> Update Max Pt for TPC Low : %f -> %f\n", fMaxPtForTPClow, maxPtForTPClow);
    fMaxPtForTPClow = maxPtForTPClow;
  }
}

///////////////////////////////////////////////////////////////////////////////////

/*
 * ---------------------------------------------------------------------------------
 *                           Initialize - Private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleHelper::InitializeEventStats() {
  // -- Initialize event statistics histograms

  fHEventStat0 = new TH1F("hEventStat0","Event cut statistics 0;Event Cuts;Events", fHEventStatMax,-0.5,fHEventStatMax-0.5);
  fHEventStat1 = new TH1F("hEventStat1","Event cut statistics 1;Event Cuts;Events", fHEventStatMax,-0.5,fHEventStatMax-0.5);

  for ( Int_t ii=0; ii < fHEventStatMax-1; ii++ ) {
    fHEventStat0->GetXaxis()->SetBinLabel(ii+1, AliAnalysisNetParticleHelper::fgkEventNames[ii]);
    fHEventStat1->GetXaxis()->SetBinLabel(ii+1, AliAnalysisNetParticleHelper::fgkEventNames[ii]);
  }

  fHEventStat0->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", AliAnalysisNetParticleHelper::fgkCentralityMaxNames[fCentralityBinMax-1]));
  fHEventStat1->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", AliAnalysisNetParticleHelper::fgkCentralityMaxNames[fCentralityBinMax-1]));
}

//________________________________________________________________________
void AliAnalysisNetParticleHelper::InitializeTriggerStats() {
  // -- Initialize trigger statistics histograms

  fHTriggerStat = new TH1F("hTriggerStat","Trigger statistics;Trigger;Events", fNTriggers,-0.5,fNTriggers-0.5);

  for ( Int_t ii=0; ii < fNTriggers; ii++ )
    fHTriggerStat->GetXaxis()->SetBinLabel(ii+1, AliAnalysisNetParticleHelper::fgkTriggerNames[ii]);
}

//________________________________________________________________________
void AliAnalysisNetParticleHelper::InitializeCentralityStats() {
  // -- Initialize trigger statistics histograms

  fHCentralityStat = new TH1F("hCentralityStat","Centrality statistics;Centrality Bins;Events", 
			      fNCentralityBins,-0.5,fNCentralityBins-0.5);

  for ( Int_t ii=0; ii < fNCentralityBins; ii++ )
    fHCentralityStat->GetXaxis()->SetBinLabel(ii+1, AliAnalysisNetParticleHelper::fgkCentralityNames[ii]);
}

/*
 * ---------------------------------------------------------------------------------
 *                         Event / Trigger Statistics - private
 * ---------------------------------------------------------------------------------
 */
  
//________________________________________________________________________
Bool_t AliAnalysisNetParticleHelper::FillEventStats(Int_t *aEventCuts) {
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
  if (!isRejected)
    fHCentralityStat->Fill(fCentralityBin);

  return isRejected;
}
