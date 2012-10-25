//-*- Mode: C++ -*-

#include "TMath.h"
#include "TAxis.h"
#include "TSystem.h" 
#include "TFile.h" 
#include "TPRegexp.h"
#include "TDatabasePDG.h"

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

// Task for NetParticle checks
// Author: Jochen Thaeder <jochen@thaeder.de>

ClassImp(AliAnalysisNetParticleHelper)

/*
 * ---------------------------------------------------------------------------------
 *                            particle names 
 * ---------------------------------------------------------------------------------
 */

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

  const Char_t* aControlPartNames[2]      = {"lambdabar",     "lambda"};
  const Char_t* aControlPartTitles[2]     = {"Anti-Lambda",   "Lambda"};

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
AliAnalysisNetParticleHelper::AliAnalysisNetParticleHelper() :
  fESDHandler(NULL),
  fPIDResponse(NULL),
  fESD(NULL),
  fAODHandler(NULL),
  fAOD(NULL),
  fMCEvent(NULL),
  fStack(NULL),

  fCentralityBin(-1),
  fCentralityPercentile(-1.),

  fCentralityBinMax(7),
  fVertexZMax(10.),
  fRapidityMax(0.5),
  fMinTrackLengthMC(70.),
  fNSigmaMaxCdd(3.),
  fNSigmaMaxCzz(3.),

  fParticleSpecies(AliPID::kProton),
  fControlParticleCode(3122),
  fControlParticleIsNeutral(kTRUE),
  fControlParticleName("Lambda"),

  fNSigmaMaxTPC(2.5),
  fNSigmaMaxTOF(2.5),
  fMinPtForTOFRequired(0.8),

  fHEventStat0(NULL),
  fHEventStat1(NULL),
  fHEventStatMax(6),

  fHTriggerStat(NULL),
  fNTriggers(5),

  fHCentralityStat(NULL),
  fNCentralityBins(10),

  fEtaCorrFunc(NULL),
  fCorr0(NULL),
  fCorr1(NULL) {
  // Constructor   
  
  AliLog::SetClassDebugLevel("AliAnalysisNetParticleHelper",10);
}

//________________________________________________________________________
AliAnalysisNetParticleHelper::~AliAnalysisNetParticleHelper() {
  // Destructor

  for (Int_t jj = 0; jj < 2; ++jj) {
    if (fCorr0[jj]) delete[] fCorr0[jj];
    if (fCorr1[jj]) delete[] fCorr1[jj];
  }
  if (fCorr0) delete[] fCorr0;
  if (fCorr1) delete[] fCorr1;
  
  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                    Setter
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleHelper::SetParticleSpecies(AliPID::EParticleType pid) {
  // -- Set particle species (ID, Name, Title, Title LATEX)

  if( (Int_t)pid < 0 || (Int_t)pid > AliPID::kSPECIES){
    AliWarning("Particle ID not in AliPID::kSPECIES --> Set to protons");
    pid = AliPID::kProton;
  }  

  fParticleSpecies     = pid;

  for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
    fPartName[idxPart]       = aPartNames[fParticleSpecies][idxPart];
    fPartTitle[idxPart]      = aPartTitles[fParticleSpecies][idxPart];
    fPartTitleLatex[idxPart] = aPartTitlesLatex[fParticleSpecies][idxPart];
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

//________________________________________________________________________
TString AliAnalysisNetParticleHelper::GetControlParticleName(Int_t idxPart) {
  // -- Get particle Title LATEX

  if( idxPart != 0 && idxPart != 1){
    AliWarning("Particle type not known --> Set to antiparticles");
    idxPart = 0;
  }

  return aControlPartNames[idxPart];
}

//________________________________________________________________________
TString AliAnalysisNetParticleHelper::GetControlParticleTitle(Int_t idxPart) {
  // -- Get particle Title LATEX

  if( idxPart != 0 && idxPart != 1){
    AliWarning("Particle type not known --> Set to antiparticles");
    idxPart = 0;
  }

  return aControlPartTitles[idxPart];
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Int_t AliAnalysisNetParticleHelper::Initialize(Bool_t isMC) {
  // -- Initialize helper

  Int_t iResult = 0;

  // -- Setup event cut statistics 
  InitializeEventStats();

  // -- Setup trigger statistics 
  InitializeTriggerStats();

  // -- Setup centrality statistics 
  InitializeCentralityStats();

  // -- Load Eta correction function 
  iResult = InitializeEtaCorrection(isMC);

  // -- Load track by track correction function 
  iResult = InitializeTrackbyTrackCorrection();

  return iResult;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleHelper::SetupEvent(AliESDInputHandler *esdHandler, AliAODInputHandler *aodHandler, AliMCEvent *mcEvent) {
  // -- Setup Event

  // -- Get ESD objects
  if(esdHandler){
    fESDHandler  = esdHandler;
    fPIDResponse = esdHandler->GetPIDResponse();
    fESD         = fESDHandler->GetEvent();
  }

  // -- Get AOD objects
  else if(aodHandler){
    fAODHandler  = aodHandler;
    fPIDResponse = aodHandler->GetPIDResponse();
    fAOD         = fAODHandler->GetEvent();
  }

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
    centrality = fAOD->GetHeader()->GetCentralityP();

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

  if(esdHandler){
    // -- Update TPC pid with eta correction (only for ESDs?) XXX
    UpdateEtaCorrectedTPCPid();
  }

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

  if(fESDHandler){
    if ((fESDHandler->IsEventSelected() & AliVEvent::kMB))          aTriggerFired[0] = kTRUE;
    if ((fESDHandler->IsEventSelected() & AliVEvent::kCentral))     aTriggerFired[1] = kTRUE;
    if ((fESDHandler->IsEventSelected() & AliVEvent::kSemiCentral)) aTriggerFired[2] = kTRUE;
    if ((fESDHandler->IsEventSelected() & AliVEvent::kEMCEJE))      aTriggerFired[3] = kTRUE;
    if ((fESDHandler->IsEventSelected() & AliVEvent::kEMCEGA))      aTriggerFired[4] = kTRUE;
  }
  else if(fAODHandler){
    if ((fAODHandler->IsEventSelected() & AliVEvent::kMB))          aTriggerFired[0] = kTRUE;
    if ((fAODHandler->IsEventSelected() & AliVEvent::kCentral))     aTriggerFired[1] = kTRUE;
    if ((fAODHandler->IsEventSelected() & AliVEvent::kSemiCentral)) aTriggerFired[2] = kTRUE;
    if ((fAODHandler->IsEventSelected() & AliVEvent::kEMCEJE))      aTriggerFired[3] = kTRUE;
    if ((fAODHandler->IsEventSelected() & AliVEvent::kEMCEGA))      aTriggerFired[4] = kTRUE;
  }

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
    if(TMath::Abs(vtxESD->GetZv()) > fVertexZMax) 
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
Bool_t AliAnalysisNetParticleHelper::IsParticleAcceptedBasicCharged(TParticle *particle, Int_t idxMC) {
  // -- Check if MC particle is accepted for basic parameters
  
  if (!particle) 
    return kFALSE;

  // -- check if PDF code exists
  if (!particle->GetPDG()) 
    return kFALSE;
    
  // -- check if charged
  if (particle->GetPDG()->Charge() == 0.0) 
    return kFALSE;
      
  // -- check if physical primary
  if(!fStack->IsPhysicalPrimary(idxMC)) 
    return kFALSE;
        
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisNetParticleHelper::IsParticleAcceptedBasicCharged(AliAODMCParticle *particle) {
  // -- Check if MC particle is accepted for basic parameters
  
  if (!particle) 
    return kFALSE;

  // -- check if PDF code exists
  if (!(TDatabasePDG::Instance()->GetParticle(particle->PdgCode()))) 
    return kFALSE;
    
  // -- check if charged
  if ((TDatabasePDG::Instance()->GetParticle(particle->PdgCode()))->Charge() == 0.0) 
    return kFALSE;
      
  // -- check if physical primary
  if(!particle->IsPhysicalPrimary()) 
    return kFALSE;
        
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisNetParticleHelper::IsParticleAcceptedBasicNeutral(TParticle *particle, Int_t idxMC) {
  // -- Check if MC particle is accepted for basic parameters
  
  if (!particle) 
    return kFALSE;

  // -- check if PDF code exists
  if (!particle->GetPDG()) 
    return kFALSE;
    
  // -- check if neutral
  if (particle->GetPDG()->Charge() != 0.0) 
    return kFALSE;
      
  // -- check if physical primary
  if(!fStack->IsPhysicalPrimary(idxMC)) 
    return kFALSE;
        
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisNetParticleHelper::IsParticleAcceptedBasicNeutral(AliAODMCParticle *particle) {
  // -- Check if MC particle is accepted for basic parameters
  
  if (!particle) 
    return kFALSE;

  // -- check if PDF code exists
  if (!(TDatabasePDG::Instance()->GetParticle(particle->PdgCode()))) 
    return kFALSE;
    
  // -- check if charged
  if ((TDatabasePDG::Instance()->GetParticle(particle->PdgCode()))->Charge() != 0.0) 
    return kFALSE;
      
  // -- check if physical primary
  if(!particle->IsPhysicalPrimary()) 
    return kFALSE;
        
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisNetParticleHelper::IsParticleAcceptedRapidity(TParticle *particle, Double_t &yP) {
  // -- Check if particle is accepted
  // > in rapidity
  // > return 0 if not accepted

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
Bool_t AliAnalysisNetParticleHelper::IsTrackAcceptedBasicCharged(AliESDtrack *track) {
  // -- Check if track is accepted 
  // > for basic parameters

  if (!track)
    return kFALSE;
  
  if (track->Charge() == 0) 
    return kFALSE;
  
  // -- Get momentum for dEdx
  if (!track->GetInnerParam()) 
    return kFALSE;
  
  return kTRUE;
} 
 
//________________________________________________________________________
Bool_t AliAnalysisNetParticleHelper::IsTrackAcceptedBasicCharged(AliAODTrack *track) {
  // -- Check if track is accepted 
  // > for basic parameters

  if (!track)
    return kFALSE;
  
  if (track->Charge() == 0) 
    return kFALSE;
  
  //// -- Get momentum for dEdx --> returns always ZERO XXX
  //if (!track->GetInnerParam()) 
  //  return kFALSE;
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisNetParticleHelper::IsTrackAcceptedRapidity(AliVTrack *track, Double_t &yP) {
  // -- Check if track is accepted
  // > in rapidity
  // > return 0 if not accepted

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
Bool_t AliAnalysisNetParticleHelper::IsTrackAcceptedDCA(AliESDtrack *track) {  //ONLY FOR ESDs so far XXX
  // -- Check if track is accepted 
  // > for DCA, if both SPD layers have hits

  Bool_t isAccepted = kTRUE;

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
  // > provides TPC and TOF nSigmas to the argument

  Bool_t isAcceptedTPC = kFALSE;
  Bool_t isAcceptedTOF = kTRUE;
  
  // -- Get PID
  pid[0] = fPIDResponse->NumberOfSigmasTPC(track, fParticleSpecies);
  pid[1] = fPIDResponse->NumberOfSigmasTOF(track, fParticleSpecies);

  // -- Check PID with TPC
  if (TMath::Abs(pid[0]) < fNSigmaMaxTPC) 
    isAcceptedTPC = kTRUE;
  
  // -- Check PID with TOF
  if (isAcceptedTPC) {

    // Has TOF PID availible
    if (track->GetStatus() & AliVTrack::kTOFpid) {
      if (TMath::Abs(pid[1]) < fNSigmaMaxTOF) 
	isAcceptedTOF = kTRUE;
      else 
	isAcceptedTOF = kFALSE;
    }
    
    // Has no TOF PID availible
    else { 
      if (track->Pt() > fMinPtForTOFRequired)
	isAcceptedTOF = kFALSE;
      else
	isAcceptedTOF = kTRUE;
    }
  } // if (isAcceptedTPC) {

  return (isAcceptedTPC && isAcceptedTOF);
}

/*
 * ---------------------------------------------------------------------------------
 *                         Helper Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleHelper::UpdateEtaCorrectedTPCPid() {
  // -- Update eta corrected TPC pid 

  if (!fEtaCorrFunc)
    return;
  
  for (Int_t idxTrack = 0; idxTrack < fESD->GetNumberOfTracks(); ++idxTrack) {
    AliESDtrack *track = fESD->GetTrack(idxTrack); 

    // -- Check if charged track is accepted for basic parameters
    if (!IsTrackAcceptedBasicCharged(track))
      continue;
    
    Double_t newTPCSignal = track->GetTPCsignal() / fEtaCorrFunc->Eval(track->Eta());
    track->SetTPCsignal(newTPCSignal, track->GetTPCsignalSigma(), track->GetTPCsignalN());
  }
}

//________________________________________________________________________
Double_t AliAnalysisNetParticleHelper::GetTrackbyTrackCorrectionFactor(Double_t *aTrack,  Int_t flag) {
  // -- Get efficiency correctionf of particle dependent on (eta, phi, pt, centrality)

  Int_t idxPart = (aTrack[2] < 0) ? 0 : 1;
  THnSparseF* corrMatrix = (flag == 0) ? fCorr0[idxPart][fCentralityBin] : fCorr1[idxPart][fCentralityBin];
  
  Double_t dimBin[3] = {aTrack[3], aTrack[4], aTrack[1]}; // eta, phi, pt    

  Double_t corr = corrMatrix->GetBinContent(corrMatrix->GetBin(dimBin));
  if (corr == 0.) {
    AliError(Form("Should not happen : bin content = 0. >> eta: %.2f | phi : %.2f | pt : %.2f | cent %d", 
		  aTrack[3], aTrack[4], aTrack[1], fCentralityBin));
    corr = 1.;
  }
  
  return corr;
}

//________________________________________________________________________
void AliAnalysisNetParticleHelper::BinLogAxis(const THnSparseF *h, Int_t axisNumber) {
  // -- Method for the correct logarithmic binning of histograms
  
  TAxis *axis = h->GetAxis(axisNumber);
  Int_t  bins = axis->GetNbins();

  Double_t from  = axis->GetXmin();
  Double_t to    = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
   
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
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

  const Char_t* aEventNames[]      = {"All", "IsTriggered", "HasVertex", "Vz<Vz_{Max}", "Centrality [0,90]%"};
  const Char_t* aCentralityMaxNames[] = {"5", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"};

  fHEventStat0 = new TH1F("fHEventStat0","Event cut statistics 0;Event Cuts;Events", fHEventStatMax,-0.5,fHEventStatMax-0.5);
  fHEventStat1 = new TH1F("fHEventStat1","Event cut statistics 1;Event Cuts;Events", fHEventStatMax,-0.5,fHEventStatMax-0.5);

  for ( Int_t ii=0; ii < fHEventStatMax-1; ii++ ) {
    fHEventStat0->GetXaxis()->SetBinLabel(ii+1, aEventNames[ii]);
    fHEventStat1->GetXaxis()->SetBinLabel(ii+1, aEventNames[ii]);
  }
  fHEventStat0->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", aCentralityMaxNames[fCentralityBinMax-1]));
  fHEventStat1->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", aCentralityMaxNames[fCentralityBinMax-1]));
}

//________________________________________________________________________
void AliAnalysisNetParticleHelper::InitializeTriggerStats() {
  // -- Initialize trigger statistics histograms

  const Char_t* aTriggerNames[] = { "kMB", "kCentral", "kSemiCentral", "kEMCEJE", "kEMCEGA" };

  fHTriggerStat = new TH1F("fHTriggerStat","Trigger statistics;Trigger;Events", fNTriggers,-0.5,fNTriggers-0.5);

  for ( Int_t ii=0; ii < fNTriggers; ii++ )
    fHTriggerStat->GetXaxis()->SetBinLabel(ii+1, aTriggerNames[ii]);
}

//________________________________________________________________________
void AliAnalysisNetParticleHelper::InitializeCentralityStats() {
  // -- Initialize trigger statistics histograms

  const Char_t* aCentralityNames[] = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", 
				      "50-60%", "60-70%", "70-80%", "80-90%", "90-100%"};
 
  fHCentralityStat = new TH1F("fHCentralityStat","Centrality statistics;Centrality Bins;Events", 
			      fNCentralityBins,-0.5,fNCentralityBins-0.5);

  for ( Int_t ii=0; ii < fNCentralityBins; ii++ )
    fHCentralityStat->GetXaxis()->SetBinLabel(ii+1, aCentralityNames[ii]);
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleHelper::InitializeEtaCorrection(Bool_t isMC) {
  // -- Initialize eta correction maps for TPC pid
  
  TFile fileEtaCorrMaps("${ALICE_ROOT}/PWGDQ/dielectron/files/EtaCorrMaps.root");
  if (!fileEtaCorrMaps.IsOpen()) 
    return -1;

  TList *keyList = fileEtaCorrMaps.GetListOfKeys();
  
  TString sList("");
  sList = (isMC) ? "LHC11a10" :  "LHC10h.pass2";
  
  for (Int_t idx = 0; idx < keyList->GetEntries(); ++idx) {
    TString keyName = keyList->At(idx)->GetName();
    TPRegexp reg(keyName);
    if (reg.MatchB(sList)){
      AliInfo(Form("Using Eta Correction Function: %s",keyName.Data()));
      fEtaCorrFunc = static_cast<TF1*>(fileEtaCorrMaps.Get(keyName.Data()));
      return 0;
    }
  }

  return -2;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleHelper::InitializeTrackbyTrackCorrection() {
  // -- Initialize track by track correction matrices

  if( fParticleSpecies != AliPID::kProton){
    AliWarning("Only proton corrections (preliminary available) --> Will correct with these");
  }
  
  // In future --> !!!!!!!!!!!!!!!!!!!!!!!
  //TFile* corrFile  = TFile::Open("${ALICE_ROOT}/PWGCF/EBYE/NetParticle/eff/effectiveCorrection.root");
  TFile* corrFile  = TFile::Open("${ALICE_ROOT}/PWGCF/EBYE/NetParticle/eff/effectiveCorrectionProton.root");
  if (!corrFile) {
    AliError("Track-by-track correction file can not be opened!");
    return -1;
  }

  // -- correction - not cross section corrected
  fCorr0    = new THnSparseF**[2];
  fCorr0[0] = new THnSparseF*[fCentralityBinMax];
  fCorr0[1] = new THnSparseF*[fCentralityBinMax];

  for (Int_t idxCent = 0; idxCent < fCentralityBinMax; ++idxCent) {

    // In future --> !!!!!!!!!!!!!!!!!!!!!!!
    // THnSparseF *sp0 = static_cast<THnSparseF*>(corrFile->Get(Form("%s_Corr0_Cent_%d", fPartName[0].Data(), idxCent)));
    // THnSparseF *sp1 = static_cast<THnSparseF*>(corrFile->Get(Form("%s_Corr0_Cent_%d", fPartName[1].Data(), idxCent)));
    THnSparseF *sp0 = static_cast<THnSparseF*>(corrFile->Get(Form("pbar_Corr0_Cent_%d", idxCent)));
    THnSparseF *sp1 = static_cast<THnSparseF*>(corrFile->Get(Form("p_Corr0_Cent_%d", idxCent)));

    if (!sp0 || !sp1) {
      AliError(Form("Effective correction objects 0 - idxCent %d can not be retrieved!",idxCent));
      return -1;
    }
    
    fCorr0[0][idxCent] = static_cast<THnSparseF*>(sp0->Clone());
    fCorr0[1][idxCent] = static_cast<THnSparseF*>(sp1->Clone());
  }

  // -- correction - ross section corrected
  fCorr1    = new THnSparseF**[2];
  fCorr1[0] = new THnSparseF*[fCentralityBinMax];
  fCorr1[1] = new THnSparseF*[fCentralityBinMax];

  for (Int_t idxCent = 0; idxCent < fCentralityBinMax; ++idxCent) {

    // In future --> !!!!!!!!!!!!!!!!!!!!!!!
    //THnSparseF *sp0 = static_cast<THnSparseF*>(corrFile->Get(Form("%s_Corr1_Cent_%d", fPartName[0].Data(), idxCent)));
    //THnSparseF *sp1 = static_cast<THnSparseF*>(corrFile->Get(Form("%s_Corr1_Cent_%d", fPartName[1].Data(), idxCent)));
    THnSparseF *sp0 = static_cast<THnSparseF*>(corrFile->Get(Form("pbar_Corr1_Cent_%d", idxCent)));
    THnSparseF *sp1 = static_cast<THnSparseF*>(corrFile->Get(Form("p_Corr1_Cent_%d", idxCent)));

    if (!sp0 || !sp1) {
      AliError(Form("Effective correction objects 1 - idxCent %d can not be retrieved!",idxCent));
      return -1;
    }

    fCorr1[0][idxCent] = static_cast<THnSparseF*>(sp0->Clone());
    fCorr1[1][idxCent] = static_cast<THnSparseF*>(sp1->Clone());
  }

  return 0;
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
