#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TH1D.h"
#include "TFile.h"
#include "TParticle.h"
#include "TParticlePDG.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"

#include "AliAnalysisTaskFilterHe3.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"

class AliAnalysisTaskFilterHe3;

using namespace std;

ClassImp(AliAnalysisTaskFilterHe3)

//____________________________________________________________________________________//
AliAnalysisTaskFilterHe3::AliAnalysisTaskFilterHe3() :AliAnalysisTaskSE(), fESD(0), fOutputList(0),
  fListOfFiles(0), fESDtrackCuts(0), fESDtrackCutsPrimary(0), fPIDResponse(0), fMultSel(0), fUseMultTaskCentrality(0), ParticleType(AliPID::kHe3),
  fEventIdFile(0),
  fFileName(0),
  fHistZv(0),
  fHistdEdxData(0),
  fHistTof(0),
  fHistCent(0)
    
{
  //
  // default contstructor: do nothing
  //
}

//____________________________________________________________________________________//
AliAnalysisTaskFilterHe3::AliAnalysisTaskFilterHe3(const char *name) : AliAnalysisTaskSE(name), fESD(0), fOutputList(0),
										 fListOfFiles(0), fESDtrackCuts(0), fESDtrackCutsPrimary(0), fPIDResponse(0), fMultSel(0), fUseMultTaskCentrality(0), ParticleType(AliPID::kHe3),
										 fEventIdFile(0),
										 fFileName(0),
										 fHistZv(0),
										 fHistdEdxData(0),
										 fHistTof(0),
										 fHistCent(0)
  										 
{
  //
  // main constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  //
  Double_t paramDandTdata[5] = {6.70549, 6.11866, 8.86205e-15, 2.34059, 1.07029};
  Double_t paramHe3data[5] = {1.48718, 27.4992, 4.00313e-15, 2.48485, 8.31768}; // Z=2 needs slight re-tuning
  //
  Double_t paramDandTmc[5] = {20.1533, 2.58127, 0.00114169, 2.0373, 0.502123};
  Double_t paramHe3mc[5] = {20.1533, 2.58127, 0.00114169, 2.0373, 0.502123}; // Z=2 needs slight re-tuning
  //
  for (Int_t iP = 0; iP < 5; iP++)
  {
    fParamDeuteron[iP] = paramDandTdata[iP];
    fParamTriton[iP] = paramDandTdata[iP];
    fParamHe3[iP] = paramHe3data[iP];
    //
    fParamDeuteronMC[iP] = paramDandTmc[iP];
    fParamTritonMC[iP] = paramDandTmc[iP];
    fParamHe3MC[iP] = paramHe3mc[iP];
  }
}

//____________________________________________________________________________________//
AliAnalysisTaskFilterHe3::~AliAnalysisTaskFilterHe3()
{
  //
  // destructor
  //
  if (fOutputList)
  {
    delete fOutputList;
  }
}

//____________________________________________________________________________________//
void AliAnalysisTaskFilterHe3::UserCreateOutputObjects()
{
  //
  // create the output objects
  //
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  //
  fListOfFiles = new TTree("fListOfFiles", "NucleusCandidates");
  fListOfFiles->Branch("fEventIdFile", &fEventIdFile, "fEventIdFile/I");
  fListOfFiles->Branch("fFileName", &fFileName, 16000, 0);

  //
  // create QA histograms
  //
  fHistZv = new TH1F("fHistZv", "fHistZv; vertex z (cm)", 200, -40, 40);
  fHistCent = new TH1F("fHistCent", "fHistCent", 240, -10.0, 110.0);
  fHistdEdxData = new TH2F("fHistdEdxData", "fHistdEdxData; p/z (GeV/#it{c}); d#it{E}/d#it{x} in TPC (arb. units)", 1000, -5.0, 5.0, 1200, 0., 1200.);
  fHistdEdxDeuteronParam[0] = new TH3F("fHistdEdxDeuteronParam", "fHistdEdxDeuteronParam; P(Gev/c); dE/dx pull; TOF mass^{2}",
                                    100, -5, +5, 200, -10., 10., 500, -5.0, +5.0);
  fHistdEdxHe3Param[0] = new TH3F("fHistdEdxHe3Param", "fHistdEdxHe3Param; P(Gev/c); dE/dx pull; TOF mass^{2}",
                               100, -5, +5, 200, -10., 10., 500, -5.0, +5.0);
  fHistdEdxTritonParam[0] = new TH3F("fHistdEdxTritonParam", "fHistdEdxTritonParam; P(Gev/c); dE/dx pull; TOF mass^{2}",
                                  100, -5, +5, 200, -10., 10., 500, -5.0, +5.0);
  fHistTof = new TH2F("fHistTof", "all paritcles TOF; P(Gev/c); mass", 1000, -10, 10, 1000, 0, 5.0); // all particles TOF quality assurance

  //Once the trigger condition (on single track) is satisfied:
  //Control containers for a first quick check (to be optimized) 
  fHistdEdxDeuteronParam[1] = new TH3F("fHistdEdxDeuteronParam_triggerON", "fHistdEdxDeuteronParam_triggerON; P(Gev/c); dE/dx pull; TOF mass^{2}",
                                    100, -5, +5, 200, -10., 10., 500, -5.0, +5.0);
  fHistdEdxHe3Param[1] = new TH3F("fHistdEdxHe3Param_triggerON", "fHistdEdxHe3Param_triggerON; P(Gev/c); dE/dx pull; TOF mass^{2}",
                               100, -5, +5, 200, -10., 10., 500, -5.0, +5.0);
  fHistdEdxTritonParam[1] = new TH3F("fHistdEdxTritonParam_triggerON", "fHistdEdxTritonParam_triggerON; P(Gev/c); dE/dx pull; TOF mass^{2}",
                                  100, -5, +5, 200, -10., 10., 500, -5.0, +5.0);
  
  //
  //
  //add histograms to output
  //
  //
  fOutputList->Add(fHistZv);
  fOutputList->Add(fHistCent);
  fOutputList->Add(fHistdEdxData);
  fOutputList->Add(fHistdEdxDeuteronParam[0]);
  fOutputList->Add(fHistdEdxHe3Param[0]);
  fOutputList->Add(fHistdEdxTritonParam[0]);
  fOutputList->Add(fHistTof);
  fOutputList->Add(fHistdEdxDeuteronParam[1]);
  fOutputList->Add(fHistdEdxHe3Param[1]);
  fOutputList->Add(fHistdEdxTritonParam[1]);
  
  //
  // introduce track cuts
  //
  fESDtrackCuts = new AliESDtrackCuts;
  fESDtrackCuts->SetMinNClustersTPC(50);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetEtaRange(-0.8, 0.8);
  //
  fESDtrackCutsPrimary = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
  fESDtrackCutsPrimary->SetMaxDCAToVertexXY(0.5);
  fESDtrackCutsPrimary->SetMaxDCAToVertexZ(2.0);
  fESDtrackCutsPrimary->SetEtaRange(-0.8, 0.8);
  //
  //
  //
  PostData(1, fOutputList);
  PostData(2, fListOfFiles);
}

//____________________________________________________________________________________//
void AliAnalysisTaskFilterHe3::UserExec(Option_t *)
{
  //
  // loop over events
  //

  Bool_t isTriggered = kFALSE;
  //
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = 0x0;
  if (man)
  {
    inputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
    if (inputHandler)
    {
      fPIDResponse = inputHandler->GetPIDResponse();
    }
  }
  fESD = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESD)
    return;
  //
  // check if MC is available and if yes, initialise stack and McEvent
  //
  Bool_t mcTrue = kFALSE;
  AliMCEventHandler *eventHandlerMC = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (eventHandlerMC)
    mcTrue = kTRUE;
  AliMCEvent *mcEvent = 0x0;
  AliStack *stack = 0x0;
  if (mcTrue)
  {
    mcEvent = eventHandlerMC->MCEvent();
    stack = mcEvent->Stack();
  }
  //
  // define some useful constants
  //
  const Double_t massD = AliPID::ParticleMass(AliPID::kDeuteron);
  const Double_t massT = AliPID::ParticleMass(AliPID::kTriton);
  const Double_t massHe3 = AliPID::ParticleMass(AliPID::kHe3);
  //
  // primary vertex
  //
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1)
  {
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if (vertex->GetNContributors() < 1)
      vertex = 0x0;
  }
  if (!vertex)
    return;
  fHistZv->Fill(vertex->GetZ());
  if (TMath::Abs(vertex->GetZ()) > 10.0)
    return; // remove events with a vertex which is more than 10cm away
  //
  // RECONSTRUCTED PARTICLES
  //
  Int_t jTracks = fESD->GetNumberOfTracks();
  for (Int_t j = 0; j < jTracks; j++)
  {

    // start loop over tracks
    //
    AliESDtrack *track = static_cast<AliESDtrack *>(fESD->GetTrack(j));
    if (!track)
      continue;
    if (!fESDtrackCuts->AcceptTrack(track))
      continue; // check if track passes the cuts
    if (!track->GetInnerParam())
      continue; // check if track is a proper TPC track
    //
    // get MC information if it is available
    //
    Long64_t label = -9999;
    Long64_t pdgCode = 0;
    if (mcTrue)
    {
      label = track->GetLabel();
      TParticle *particleGen = ((AliMCParticle *)mcEvent->GetTrack(TMath::Abs(label)))->Particle();
      if (particleGen)
        pdgCode = particleGen->GetPdgCode();
    }
    //
    // basic track quantities
    //
    Double_t ptot = track->GetInnerParam()->GetP(); // momentum for dEdx determination
    Double_t sign = track->GetSign();               // charge
    Double_t tpcSignal = track->GetTPCsignal();
    Double_t pT = track->Pt();
    // DCA_XY
    //
    Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
    track->GetImpactParameters(dca, cov);
    //
    // fill the first QA dE/Dx histograms (Figure 1 left in paper draft)
    //
    if (track->GetTPCsignalN() > 50 && fESDtrackCutsPrimary->AcceptTrack(track))
      fHistdEdxData->Fill(ptot * sign, tpcSignal);
    //
    // get TOF signal if it is available
    //
    Bool_t hasTOF = kFALSE;
    UInt_t status = track->GetStatus();
    Bool_t hasTOFout = status & AliESDtrack::kTOFout;
    Bool_t hasTOFtime = status & AliESDtrack::kTIME;
    Double_t length = track->GetIntegratedLength();
    if (length > 350. && hasTOFout && hasTOFtime)
      hasTOF = kTRUE;
    //
    Double_t time0 = fPIDResponse->GetTOFResponse().GetStartTime(track->P());
    Double_t mass = -999.;
    Double_t time = -1;
    Double_t beta = 0;
    //
    if (hasTOF)
    {
      time = track->GetTOFsignal() - time0;
      if (time > 0)
      {
        beta = length / (TMath::C() * 1e-10 * time);
        Double_t gamma = 1 / TMath::Sqrt(1 - beta * beta);
        mass = ptot / TMath::Sqrt(gamma * gamma - 1); // using inner TPC mom. as approx.
      }
    }
    //
    // calculate TPC pid pulls and QA them including TOF
    //
    Float_t deutExp = -999;
    Float_t tritExp = -999;
    Float_t hel3Exp = -999;
    //
    if (ptot > 0.3)
    { // protection against numerical instabilities
      if (!mcTrue)
      {
        deutExp = AliExternalTrackParam::BetheBlochAleph(ptot / massD, fParamDeuteron[0], fParamDeuteron[1], fParamDeuteron[2], fParamDeuteron[3], fParamDeuteron[4]);
        tritExp = AliExternalTrackParam::BetheBlochAleph(ptot / massT, fParamTriton[0], fParamTriton[1], fParamTriton[2], fParamTriton[3], fParamTriton[4]);
        hel3Exp = 4.0 * AliExternalTrackParam::BetheBlochAleph(2.0 * ptot / massHe3, fParamHe3[0], fParamHe3[1], fParamHe3[2], fParamHe3[3], fParamHe3[4]);
      }
      else
      {
        deutExp = AliExternalTrackParam::BetheBlochAleph(ptot / massD, fParamDeuteronMC[0], fParamDeuteronMC[1], fParamDeuteronMC[2], fParamDeuteronMC[3], fParamDeuteronMC[4]);
        tritExp = AliExternalTrackParam::BetheBlochAleph(ptot / massT, fParamTritonMC[0], fParamTritonMC[1], fParamTritonMC[2], fParamTritonMC[3], fParamTritonMC[4]);
        hel3Exp = 4.0 * AliExternalTrackParam::BetheBlochAleph(2.0 * ptot / massHe3, fParamHe3MC[0], fParamHe3MC[1], fParamHe3MC[2], fParamHe3MC[3], fParamHe3MC[4]);
      }
    }
    //
    // fill QA and raw yield histograms
    //
    const Float_t avDeDxRes = 0.07; // approx dEdx resolution of 7%
    Double_t nSigmaDeut = (tpcSignal - deutExp) / (avDeDxRes * deutExp);
    Double_t nSigmaHe3 = (tpcSignal - hel3Exp) / (avDeDxRes * hel3Exp);
    Double_t nSigmaTrit = (tpcSignal - tritExp) / (avDeDxRes * tritExp);

    Double_t nSigmaA;
    switch((Int_t) ParticleType) {
    case (Int_t)AliPID::kHe3 : nSigmaA=nSigmaHe3; break;//printf("He3\n"); 
    case (Int_t)AliPID::kTriton : nSigmaA=nSigmaTrit; break;//printf("Triton\n"); 
    }

    //print cuts (as check)
    /*printf("Nsigma min = %f  max = %f\n", fMinNSigma, fMaxNSigma);
    printf("Pos Ptot min = %f max = %f\n", fMinPtotPos, fMaxPtotPos);
    printf("Neg Ptot min = %f  max = %f\n", fMinPtotNeg, fMaxPtotNeg);
    printf("Mass min = %f max = %f\n", fMinMass,fMaxMass);*/
    
    //
    // deuterons etc (positive and negative)
    //
    if (!mcTrue || pdgCode == 1000010020)
      fHistdEdxDeuteronParam[0]->Fill(ptot * sign, nSigmaDeut, mass * mass - massD * massD); // QA histogram
    if (!mcTrue || TMath::Abs(pdgCode) == 1000020030)
      fHistdEdxHe3Param[0]->Fill(ptot * sign, nSigmaHe3, mass * mass - massHe3 * massHe3); // QA histogram
    if (!mcTrue || pdgCode == 1000010030)
      fHistdEdxTritonParam[0]->Fill(ptot * sign, nSigmaTrit, mass * mass - massT * massT); // QA histogram
                                                                                 //
    // TRIGGER CONDITION
    //

    Bool_t isTriggeredCloneOnSingleTrack = kFALSE;
    
    Bool_t PosAndMomInRng = kFALSE;
    Bool_t NegAndMomInRng = kFALSE;
    if (sign < 0)
    {
      if (ptot > fMinPtotNeg && ptot < fMaxPtotNeg)
        NegAndMomInRng = kTRUE;
    }
    else if (ptot > fMinPtotPos && ptot < fMaxPtotPos)
      PosAndMomInRng = kTRUE;
    
    if (hasTOF && nSigmaA < fMaxNSigma && nSigmaA > fMinNSigma)
    {
      fHistTof->Fill(ptot * sign, mass);
      //
      if (fillSecifTOF && fMinMass < mass && mass < fMaxMass && track->GetTPCsignalN() > fMinNclsTPC)
	{
	  if (PosAndMomInRng || NegAndMomInRng) {
	    isTriggered = kTRUE;
	    isTriggeredCloneOnSingleTrack = kTRUE;
	  }
	}
    }
    
    if (PosAndMomInRng || NegAndMomInRng)
    {
      if (track->GetTPCsignalN() > fMinNclsTPC && fESDtrackCutsPrimary->AcceptTrack(track) && nSigmaA < fMaxNSigma && nSigmaA > fMinNSigma) {
        isTriggered = kTRUE;
	isTriggeredCloneOnSingleTrack = kTRUE;
      }
      if (fillSecifTOF == kFALSE && nSigmaA < fMaxNSigma && nSigmaA > fMinNSigma && fMinMass < mass && mass < fMaxMass && track->GetTPCsignalN() > fMinNclsTPC)
        isTriggered = kTRUE;
    }
    
    if(isTriggeredCloneOnSingleTrack) {//Check trigger
      fHistdEdxDeuteronParam[1]->Fill(ptot * sign, nSigmaDeut, mass * mass - massD * massD);
      fHistdEdxHe3Param[1]->Fill(ptot * sign, nSigmaHe3, mass * mass - massHe3 * massHe3);
      fHistdEdxTritonParam[1]->Fill(ptot * sign, nSigmaTrit, mass * mass - massT * massT);
    }
    
  } // end track loop
  //
  // get the file Name
  //
  if (isTriggered)
  {
    fFileName = inputHandler->GetTree()->GetCurrentFile()->GetName();
    fEventIdFile = fESD->GetHeader()->GetEventNumberInFile();
    fListOfFiles->Fill();
  }
  //
  // post the data and end the event loop
  //
  PostData(1, fOutputList);
  PostData(2, fListOfFiles);
}

//____________________________________________________________________________________//

void AliAnalysisTaskFilterHe3::Terminate(Option_t *)
{
}

//____________________________________________________________________________________//
