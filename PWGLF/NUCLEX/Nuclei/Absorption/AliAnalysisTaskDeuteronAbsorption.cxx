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

#include "AliAnalysisTaskDeuteronAbsorption.h"

#include <iostream>

#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>

#include <AliAnalysisManager.h>
#include <AliESDEvent.h>
#include <AliESDInputHandler.h>
#include <AliESDtrack.h>
#include <AliESDtrackCuts.h>
#include <AliMCEvent.h>
#include <AliMCEventHandler.h>
#include <AliPIDResponse.h>
#include <AliMultSelection.h>

class AliAnalysisTaskDeuteronAbsorption;

const AliPID::EParticleType AliAnalysisTaskDeuteronAbsorption::fgkSpecies[kNabsSpecies] = {AliPID::kKaon, AliPID::kProton, AliPID::kDeuteron, AliPID::kTriton, AliPID::kHe3};
const std::string AliAnalysisTaskDeuteronAbsorption::fgkParticleNames[kNabsSpecies] = {"Kaon", "Proton", "Deuteron", "Triton","3He"};
const double AliAnalysisTaskDeuteronAbsorption::fgkPhiParamPos[4][4] = {
      {1.38984e+00, -2.10187e+01, 5.81724e-02, 1.91938e+01},
      {2.02372e+00, -2.44456e+00, 8.99000e-01, 9.22399e-01},
      {4.21954e+00, -2.56555e+01, 4.17557e-02, 2.40301e+01},
      {5.17499e+00, -2.69241e+00, 6.97167e-01, 1.25974e+00}};
const double AliAnalysisTaskDeuteronAbsorption::fgkPhiParamNeg[4][4] = {
      {2.81984e+00, -1.81497e-01, -2.03494e+00, 2.64148e-01},
      {5.79322e+00, -5.44966e-02, -1.10803e+00, 1.29737e+00},
      {5.60000e+00, -2.06000e-01, -1.97130e+00, 2.67181e-01},
      {9.72180e+00, -4.35801e-02, -1.14550e+00, 1.49160e+00}};

ClassImp(AliAnalysisTaskDeuteronAbsorption); // classimp: necessary for root

AliAnalysisTaskDeuteronAbsorption::AliAnalysisTaskDeuteronAbsorption(const char *name) : AliAnalysisTaskSE(name),
                                                                                         fUseTRDboundariesCut{true},
                                                                                         fNtpcSigmas{5.},
                                                                                         fEventCuts{},
                                                                                         fUseTrackCuts{false},
                                                                                         fMindEdx{100.},
                                                                                         fMinTPCsignalN{50},
											 fUseTOFallClustersInfo{false},
											 ParticleType(AliPID::kHe3),
                                                                                         fPIDResponse{nullptr},
                                                                                         fESDtrackCuts{*AliESDtrackCuts::GetStandardITSTPCTrackCuts2011()},
                                                                                         fOutputList{nullptr},
                                                                                         fTreeTrack{nullptr},
											 tCentrality{-1},
											 tNglobalTracks{-1},	   
											 tPt{-999.},
                                                                                         tEta{-999.},
                                                                                         tPhi{-999.},
											 tSign{0.},
											 tdEdx{-999.},
											 tdEdxExp{-999.},
											 tdEdxExpSigma{-999.},
											 tnsigTPC{-999.},
                                                                                         tnsigTOF{-999.},
                                                                                         tmass2{-999.},
											 tITSchi2{-999.},   
                                                                                         tTOFsigDx{-999.},
                                                                                         tTOFsigDz{-999.},
                                                                                         tTOFchi2{-999.},
                                                                                         tTPCchi2{-999.},
                                                                                         tTPCxRows{-999.},
                                                                                         tTPCxRowsOverFindable{-999.},
                                                                                         tDCAxy{-999.},
                                                                                         tDCAz{-999.},
                                                                                         tT0res{-1.},
                                                                                         tT0mask{-1},
                                                                                         tTRDclsN{0},
                                                                                         tTRDntracklets{0},
                                                                                         tTRDNchamberdEdx{0},
                                                                                         tID{0},
                                                                                         tPdgCodeMc{0},
                                                                                         tTOFclsN{0},
											 tTOFchannel{-1},
                                                                                         tnPIDclsTPC{0},
                                                                                         tITSclsMap{0u},
                                                                                         tMCpt{0.f},
                                                                                         tMCabsMom{-1.},
                                                                                         tMCabsRadius{-1.},
											 tMCabsEta{-999.},
											 tMCabsPhi{-999.},
                                                                                         tMCtofMismatch{false},
											 tNdaughters{-99},
											 tIsReconstructed{false},
											 thasTOF{false},
											 tRunNumber{0},
											 tPIDforTracking{99},
                                                                                         fHistZv{nullptr},
                                                                                         fHist3TPCpid{nullptr},
                                                                                         fHist3TPCpidAll{nullptr},
                                                                                         fHist3TOFpid{nullptr},
                                                                                         fHist3TOFpidAll{nullptr},
                                                                                         fHist3TOFmass{nullptr},
                                                                                         fHist3TOFnsigma{nullptr},
                                                                                         fHist3TOFmassAll{nullptr}                                                                                        
{
  fESDtrackCuts.SetEtaRange(-0.8, 0.8);

  // constructor
  DefineInput(0, TChain::Class()); // define the input of the analysis: in this case we take a 'chain' of events
                                   // this chain is created by the analysis manager, so no need to worry about it,
                                   // it does its work automatically
  DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this case it's a list of histograms
                                   // you can add more output objects by calling DefineOutput(2, classname::Class())
                                   // if you add more output objects, make sure to call PostData for all of them, and to
//  if (fTreemode)
  DefineOutput(2, TTree::Class());
}

AliAnalysisTaskDeuteronAbsorption::~AliAnalysisTaskDeuteronAbsorption()
{
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  
  if (fOutputList)
    delete fOutputList; // at the end of your task, it is deleted from memory by calling this function

  if (fTreeTrack)
    delete fTreeTrack;
}

void AliAnalysisTaskDeuteronAbsorption::UserCreateOutputObjects()
{
  // create output objects
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man)
  {
    AliInputEventHandler *inputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
    if (inputHandler)
      fPIDResponse = inputHandler->GetPIDResponse();
  }

  // histograms used in the analysis
  // to an output file
  fOutputList = new TList();    // this is a list which will contain all of your histograms
  fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all objects it contains and will delete them
  //
  fHistZv = new TH1F("fHistZv", "fHistZv", 200, -40, 40); // histogram to monitor z-position of the primary vertex -- quality assurance
  fOutputList->Add(fHistZv);
  //

  std::string wTRD[2] = {"woTRD", "wTRD"};
  std::string wTOF[2] = {"woTOF", "wTOF"};
  std::string pos_neg[2] = {"neg", "pos"};

  for (Int_t iSpecies = 0; iSpecies < kNabsSpecies; ++iSpecies)
  {
    fHist3TPCpid[iSpecies] = new TH3F(Form("fHist3TPCpid%s", fgkParticleNames[iSpecies].data()), Form("%s; #it{p}/#it{z} (Gev/#it{c}); d#it{E}/d#it{x} (arb. units); #Phi (rad)", fgkParticleNames[iSpecies].data()), 400, -10, 10, 400, 0, 1000, 18, 0, TMath::TwoPi());
    fHist3TOFpid[iSpecies] = new TH3F(Form("fHist3TOFpid%s", fgkParticleNames[iSpecies].data()), Form("%s; #it{p}/#it{z} (Gev/#it{c}); #beta; #Phi (rad)", fgkParticleNames[iSpecies].data()), 400, -10, 10, 300, 0, 1.2, 18, 0, 2 * TMath::Pi());
    fHist3TOFmass[iSpecies] = new TH3F(Form("fHist3TOFmass%s", fgkParticleNames[iSpecies].data()), Form("%s; #it{p}/#it{z} (Gev/#it{c}); TOF m^{2} (GeV/#it{c}^{2})^{2}; #Phi (rad)", fgkParticleNames[iSpecies].data()), 400, -10, 10, 160, 0, 6.5, 18, 0, 2 * TMath::Pi());
    fHist3TOFnsigma[iSpecies] = new TH3F(Form("fHist3TOFnsigma%s", fgkParticleNames[iSpecies].data()), Form("%s; #it{p}/#it{z} (Gev/#it{c}); TOF n_{#sigma}; #Phi (rad)", fgkParticleNames[iSpecies].data()), 400, -10, 10, 300, -15, 15, 18, 0, 2 * TMath::Pi());
    fOutputList->Add(fHist3TPCpid[iSpecies]);
    fOutputList->Add(fHist3TOFpid[iSpecies]);
    fOutputList->Add(fHist3TOFmass[iSpecies]);
    fOutputList->Add(fHist3TOFnsigma[iSpecies]);
  }
  fHist3TPCpidAll = new TH3F("fHist3TPCpidAll", "; #it{p}/#it{z} (Gev/#it{c}); d#it{E}/d#it{x} (arb. units); #Phi (rad)", 400, -10, 10, 1000, 0, 1000, 18, 0, TMath::TwoPi());
  fHist3TOFpidAll = new TH3F("fHist3TOFpidAll", "; #it{p}/#it{z} (Gev/#it{c}); #beta; #Phi (rad)", 400, -10, 10, 300, 0, 1.2, 18, 0, 2 * TMath::Pi());
  fHist3TOFmassAll = new TH3F("fHist3TOFmassAll", "; #it{p}/#it{z} (Gev/#it{c}); TOF m^{2} (GeV/#it{c}^{2})^{2}; #Phi (rad)", 400, -10, 10, 160, 0, 6.5, 18, 0, 2 * TMath::Pi());
  fOutputList->Add(fHist3TPCpidAll);
  fOutputList->Add(fHist3TOFpidAll);
  fOutputList->Add(fHist3TOFmassAll);

  // Tree
  if (fTreemode)
  {
    OpenFile(2);
    fTreeTrack = new TTree("fTreeTrack", "Track Parameters");
    //fTreeTrack->Branch("tP", &tP, "tP/D");
    fTreeTrack->Branch("tCentrality", &tCentrality, "tCentrality/F");
    fTreeTrack->Branch("tNglobalTracks", &tNglobalTracks, "tNglobalTracks/I");
    fTreeTrack->Branch("tPt", &tPt, "tPt/F");
    fTreeTrack->Branch("tEta", &tEta, "tEta/F");
    fTreeTrack->Branch("tPhi", &tPhi, "tPhi/F");
    fTreeTrack->Branch("tSign", &tSign, "tSign/F");
    fTreeTrack->Branch("tdEdx", &tdEdx, "tdEdx/F");
    fTreeTrack->Branch("tdEdxExp", &tdEdxExp, "tdEdxExp/F");
    fTreeTrack->Branch("tdEdxExpSigma", &tdEdxExpSigma, "tdEdxExpSigma/F");
    fTreeTrack->Branch("tnsigTPC", &tnsigTPC, "tnsigTPC/F");
    fTreeTrack->Branch("tnsigTOF", &tnsigTOF, "tnsigTOF/F");
    fTreeTrack->Branch("tmass2", &tmass2, "tmass2/F");
    fTreeTrack->Branch("thasTOF", &thasTOF, "thasTOF/O");
    fTreeTrack->Branch("tnPIDclsTPC", &tnPIDclsTPC, "tnPIDclsTPC/b");
    fTreeTrack->Branch("tTOFsigDx", &tTOFsigDx, "tTOFsigDx/F");
    fTreeTrack->Branch("tTOFsigDz", &tTOFsigDz, "tTOFsigDz/F");
    fTreeTrack->Branch("tTOFchi2", &tTOFchi2, "tTOFchi2/F");
    fTreeTrack->Branch("tTOFclsN", &tTOFclsN, "tTOFclsN/b");
    fTreeTrack->Branch("tTOFchannel", &tTOFchannel, "tTOFchannel/I");
    fTreeTrack->Branch("tTRDclsN", &tTRDclsN, "tTRDclsN/I");
    fTreeTrack->Branch("tTRDntracklets", &tTRDntracklets, "tTRDntracklets/b");
    fTreeTrack->Branch("tTRDNchamberdEdx", &tTRDNchamberdEdx, "tTRDNchamberdEdx/b");
    fTreeTrack->Branch("tID", &tID, "tID/I");
    fTreeTrack->Branch("tPdgCodeMc", &tPdgCodeMc, "tPdgCodeMc/I");
    fTreeTrack->Branch("tTPCxRowsOverFindable", &tTPCxRowsOverFindable, "tTPCxRowsOverFindable/F");
    fTreeTrack->Branch("tITSchi2", &tITSchi2, "tITSchi2/F");         
    fTreeTrack->Branch("tTPCchi2", &tTPCchi2, "tTPCchi2/F");         
    fTreeTrack->Branch("tTPCxRows", &tTPCxRows, "tTPCxRows/F");        
    fTreeTrack->Branch("tDCAxy", &tDCAxy, "tDCAxy/F");           
    fTreeTrack->Branch("tDCAz", &tDCAz, "tDCAz/F");          
    fTreeTrack->Branch("tT0res", &tT0res, "tT0res/F");
    fTreeTrack->Branch("tT0mask", &tT0mask, "tT0mask/I");  
    fTreeTrack->Branch("tITSclsMap", &tITSclsMap, "tITSclsMap/b");
    fTreeTrack->Branch("tMCpt", &tMCpt, "tMCpt/F");
    fTreeTrack->Branch("tMCabsMom", &tMCabsMom, "tMCabsMom/F");
    fTreeTrack->Branch("tMCabsRadius", &tMCabsRadius, "tMCabsRadius/F");
    fTreeTrack->Branch("tMCabsEta", &tMCabsEta, "tMCabsEta/F");
    fTreeTrack->Branch("tMCabsPhi", &tMCabsPhi, "tMCabsPhi/F");
    fTreeTrack->Branch("tMCtofMismatch", &tMCtofMismatch, "tMCtofMismatch/O");
    fTreeTrack->Branch("tNdaughters", &tNdaughters, "tNdaughters/I");
    fTreeTrack->Branch("tIsReconstructed", &tIsReconstructed, "tIsReconstructed/O");
    fTreeTrack->Branch("tRunNumber", &tRunNumber, "tRunNumber/I");
    fTreeTrack->Branch("tPIDforTracking", &tPIDforTracking, "tPIDforTracking/b");

    fTreeTrack->Branch("tTOFclsN", &tTOFclsN, "tTOFclsN/b");
    if(fUseTOFallClustersInfo) {
      fTreeTrack->Branch("tNsigmaTOFarray", tNsigmaTOFarray, "[tTOFclsN]/F");
      fTreeTrack->Branch("tmass2array", tmass2array, "tmass2array[tTOFclsN]/F");
      fTreeTrack->Branch("tdxTOFarray", tdxTOFarray, "tdxTOFarray[tTOFclsN]/F");
      fTreeTrack->Branch("tdzTOFarray", tdzTOFarray, "tdzTOFarray[tTOFclsN]/F");
      fTreeTrack->Branch("tTOFarray", tTOFarray, "tTOFarray[tTOFclsN]/F");
      fTreeTrack->Branch("tLengtharray", tLengtharray, "tLengtharray[tTOFclsN]/F");
      fTreeTrack->Branch("texpTOFarray", texpTOFarray, "texpTOFarray[tTOFclsN]/F");
      fTreeTrack->Branch("tsigmaTOFarray", tsigmaTOFarray, "tsigmaTOFarray[tTOFclsN]/F");
      fTreeTrack->Branch("tTOFchannelarray", tTOFchannelarray, "tTOFchannelarray[tTOFclsN]/I");
      fTreeTrack->Branch("tMCtofMismatcharray", tMCtofMismatcharray, "tMCtofMismatcharray[tTOFclsN]/O");
      fTreeTrack->Branch("tNmatchableTracks", tNmatchableTracks, "tNmatchableTracks[tTOFclsN]/b");
    }
  }
  fEventCuts.AddQAplotsToList(fOutputList);

  PostData(1, fOutputList); // postdata will notify the analysis manager of changes / updates to the

  if (fTreemode)
    PostData(2, fTreeTrack);
}

void AliAnalysisTaskDeuteronAbsorption::UserExec(Option_t *)
{
  // main loop over events
  AliESDEvent* esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent)
    ::Fatal("AliAnalysisTaskDeuteronAbsorption::UserExec","No ESD event found.");  // if the pointer to the event is empty (getting it failed) skip this event
  Int_t nTracks = esdEvent->GetNumberOfTracks(); // see how many tracks there are in the event

  Bool_t isMC = false;
  AliMCEvent *mcEvent = nullptr;

  AliMCEventHandler * eventHandlerMC = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (eventHandlerMC)
  {
    mcEvent = eventHandlerMC->MCEvent();
    isMC = (mcEvent != nullptr);
  }

  if (!fEventCuts.AcceptEvent(esdEvent))
    return;

  // check for a proper primary vertex and monitor
  const AliESDVertex *vertex = esdEvent->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1)
  {
    // SPD vertex
    vertex = esdEvent->GetPrimaryVertexSPD();
    if (vertex->GetNContributors() < 1)
      vertex = nullptr;
  }
  if (vertex == nullptr) {
    PostData(1, fOutputList);
    return;
  }

  fHistZv->Fill(vertex->GetZ());
  if (TMath::Abs(vertex->GetZ()) > 10.0) {
    PostData(1, fOutputList);
    return; // remove events with a vertex which is more than 10cm away
  }

  tCentrality = ((AliMultSelection *) esdEvent->FindListObject("MultSelection"))->GetMultiplicityPercentile("V0M");
  // track loop
  std::vector<int> usedMC;

  tNglobalTracks = 0;
  for (Int_t i = 0; i < nTracks; i++)
  {
    AliESDtrack *track = static_cast<AliESDtrack *>(esdEvent->GetTrack(i)); 
    if (!track)
      continue;
    if (!(track->GetStatus() & AliVTrack::kTPCrefit))
      continue;
    if (!(track->GetStatus() & AliVTrack::kITSrefit))
      continue;
    if (!track->GetInnerParam())
      continue;                                   
    if (track->GetTPCsignalN() < fMinTPCsignalN)
      continue;
    if(TMath::Abs(track->GetInnerParam()->Eta()) > 0.8)
      continue;
    if(track->GetTPCsignalN() < 50)
      continue;
    if((track->GetTPCchi2() / track->GetTPCncls()) > 4)
      continue;

    UChar_t ITSclsMap = track->GetITSClusterMap();

    Int_t NclsITS=0;
    for(Int_t lay=0; lay < 6; lay++)
      if(TESTBIT(ITSclsMap, lay)) NclsITS++;

    Bool_t isSPD = TESTBIT(ITSclsMap, 0) || TESTBIT(ITSclsMap, 1);
    
    if((track->GetITSchi2()/Float_t(NclsITS)) > 36)
      continue;

    if(!isSPD)
      continue;

    tNglobalTracks++;
  }
  
  for (Int_t i = 0; i < nTracks; i++)
  {                                                                     
    AliESDtrack *track = static_cast<AliESDtrack *>(esdEvent->GetTrack(i)); // get a track (type AliESDDTrack) from the event
    if (!track)
      continue;
    if (!(track->GetStatus() & AliVTrack::kTPCrefit))
      continue;
    if (!(track->GetStatus() & AliVTrack::kITSrefit))
      continue;
    if (!fESDtrackCuts.AcceptTrack(track) && fUseTrackCuts)
      continue; // check if track passes the cuts
    if (!track->GetInnerParam())
      continue;                                       // check if track is a proper TPC track
    if (track->GetTPCsignalN() < fMinTPCsignalN)
      continue;

    // Process TOF information
    ULong_t status = (ULong_t)track->GetStatus();
    Bool_t hasTOFout  = status & AliVTrack::kTOFout;
    Bool_t hasTOFtime = status & AliVTrack::kTIME;
    const Float_t length = track->GetIntegratedLength();
    Bool_t hasTOF = hasTOFout && hasTOFtime;
    //
    Float_t ptot = track->GetTPCmomentum(); // momentum for dEdx determination
    Float_t tof = track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(track->P());
    //
    Float_t beta = -1.;
    Float_t mass2 = -1;
    //
    if (hasTOF) {
      beta = length / (TMath::C() * 1.e-10 * tof);
      if ((1 - beta * beta) > 0)
        mass2 = ptot * ptot * (1. / (beta * beta) - 1.);
    }

    Int_t pdgCodeTrackMc = 0;
    if (isMC)
    {
      AliVParticle *mcParticle = mcEvent->GetTrack(TMath::Abs(track->GetLabel()));
      pdgCodeTrackMc = mcParticle->PdgCode();
      tMCpt = mcParticle->Pt();
      usedMC.push_back(TMath::Abs(track->GetLabel()));
      Int_t tofL[3];
      track->GetTOFLabel(tofL);
      tMCtofMismatch = tofL[0] != TMath::Abs(track->GetLabel());
      if (std::abs(pdgCodeTrackMc) > 1000000)
      {

        Int_t counter = 0;
        Double_t totalMom[3]{0.};
        Double_t absVtx[3]{0., 0., 0.};
        Double_t absT{0.};
	tNdaughters = 0;
        for (Int_t c = mcParticle->GetDaughterFirst(); c <= mcParticle->GetDaughterLast(); c++)
        {
          AliVParticle *dPart = mcEvent->GetTrack(c);
	  if(!dPart) {
	    continue;
	  }
	  Double_t currentT{dPart->Tv()};
          if (counter == 0)
          {
            absT = currentT;
            dPart->XvYvZv(absVtx);
	    tMCabsEta = dPart->Eta();
	    tMCabsPhi = dPart->Phi();
	  }
          else if (std::abs(currentT - absT) > 1.e-10)
	    continue;
	  tNdaughters++;
          counter++;
          totalMom[0] += dPart->Px();
          totalMom[1] += dPart->Py();
          totalMom[2] += dPart->Pz();
	}
        tMCabsMom = std::sqrt(totalMom[0] * totalMom[0] + totalMom[1] * totalMom[1] + totalMom[2] * totalMom[2]);
        tMCabsRadius = std::hypot(absVtx[0], absVtx[1]);
      }
    }
    
    if (fTreemode && track->GetTPCsignal() > fMindEdx && fPIDResponse->NumberOfSigmasTPC(track, ParticleType) < fMaxNSigma && fPIDResponse->NumberOfSigmasTPC(track, ParticleType) > fMinNSigma)
    {

      //tP = track->GetInnerParam()->GetP();
      tPt = track->GetInnerParam()->GetSignedPt();
      tEta = track->GetInnerParam()->Eta();
      tPhi = track->GetInnerParam()->Phi();
      tSign = track->GetSign();
      tdEdx = track->GetTPCsignal();
      tdEdxExp = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, ParticleType);
      tdEdxExpSigma = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, ParticleType);
      tmass2 = mass2;
      tnPIDclsTPC = track->GetTPCsignalN();
      tTOFsigDx = track->GetTOFsignalDx();
      tTOFsigDz = track->GetTOFsignalDz();
      tTOFchi2 = track->GetTOFchi2();
      tTOFchannel = track->GetTOFCalChannel();
      tTRDclsN = track->GetTRDncls();
      tTRDntracklets = track->GetTRDntracklets();
      tTRDNchamberdEdx = track->GetTRDNchamberdEdx();
      tID = track->GetID();
      tnsigTPC = fPIDResponse->NumberOfSigmasTPC(track, ParticleType);
      tnsigTOF = fPIDResponse->NumberOfSigmasTOF(track, ParticleType);
      tPdgCodeMc = pdgCodeTrackMc;
      tITSchi2 = track->GetITSchi2();
      tTPCchi2 = track->GetTPCchi2() / track->GetTPCncls();
      tTPCxRows = track->GetTPCCrossedRows();
      tTPCxRowsOverFindable = tTPCxRows / track->GetTPCNclsF();
      track->GetImpactParameters(tDCAxy, tDCAz);
      tT0res = fPIDResponse->GetTOFResponse().GetStartTimeRes(track->GetTPCmomentum());
      tT0mask = fPIDResponse->GetTOFResponse().GetStartTimeMask(track->GetTPCmomentum());
      tITSclsMap = track->GetITSClusterMap();
      tIsReconstructed = true;
      thasTOF = hasTOF;
      tRunNumber = esdEvent->GetRunNumber();
      tPIDforTracking = track->GetPIDForTracking();

      tTOFclsN = track->GetTOFclusterN();
      if(fUseTOFallClustersInfo)
	FillTOFallMatchableHitsInfo(track);
      fTreeTrack->Fill();  
    }
    
    Float_t sign = track->GetSign();
    // fill QA histograms
    Float_t tpcNsigmas[kNabsSpecies];
    for (Int_t iSpecies = 0; iSpecies < kNabsSpecies; ++iSpecies) tpcNsigmas[iSpecies] = 999.;
    //
    fHist3TPCpidAll->Fill(ptot * sign, track->GetTPCsignal(), track->Phi());
    if (hasTOF && track->GetTPCsignal() > fMindEdx)
    {
      fHist3TOFpidAll->Fill(ptot * sign, beta, track->Phi());
      fHist3TOFmassAll->Fill(ptot * sign, mass2, track->Phi());
    }
    for (Int_t iSpecies = 0; iSpecies < kNabsSpecies; ++iSpecies)
    {
      tpcNsigmas[iSpecies] = fPIDResponse->NumberOfSigmasTPC(track, fgkSpecies[iSpecies]);
      if (std::abs(tpcNsigmas[iSpecies]) < fNtpcSigmas)
      {
        fHist3TPCpid[iSpecies]->Fill(ptot * sign, track->GetTPCsignal(), track->Phi());
        if (hasTOF && track->GetTPCsignal() > fMindEdx)
        {
          fHist3TOFpid[iSpecies]->Fill(ptot * sign, beta, track->Phi());
          fHist3TOFmass[iSpecies]->Fill(ptot * sign, mass2, track->Phi());
	  fHist3TOFnsigma[iSpecies]->Fill(ptot * sign, fPIDResponse->NumberOfSigmasTOF(track, fgkSpecies[iSpecies]), track->Phi());
        }
      }
    }
  } // end the track loop

  if(mcEvent && fTreemode) {
    for (Int_t iMC=0; iMC<fMCEvent->GetNumberOfTracks(); iMC++) {
      AliVParticle *mcParticle = mcEvent->GetTrack(iMC);
      tPdgCodeMc = mcParticle->PdgCode();
      if (std::abs(tPdgCodeMc) < 1000000000 || std::abs(tPdgCodeMc) > 1000020040) continue;

      if (std::find(usedMC.begin(), usedMC.end(), iMC) == usedMC.end()) {
        tMCpt = mcParticle->Pt();
        tIsReconstructed = false;
        fTreeTrack->Fill();
      } 
    }
  }
  
  // post the data
  PostData(1, fOutputList);
  PostData(2, fTreeTrack);
} // end the UserExec
//___________________________________________________________________________________________________________-
void AliAnalysisTaskDeuteronAbsorption::FillTOFallMatchableHitsInfo(AliESDtrack *track) {

  const Int_t fNtofClusters = track->GetTOFclusterN();
  Int_t *fTOFcluster = track->GetTOFclusterArray();// [fNtofClusters]
  
  Float_t NsigmaTOFarray[fNtofClusters], mass2array[fNtofClusters], dxTOFarray[fNtofClusters], dzTOFarray[fNtofClusters], TOFarray[fNtofClusters], Lengtharray[fNtofClusters], expTOFarray[fNtofClusters], sigmaTOFarray[fNtofClusters];
  Int_t NmatchableTracks[fNtofClusters], TOFchannelarray[fNtofClusters];
  Bool_t MCtofMismatcharray[fNtofClusters];

  for(Int_t i=0;i<fNtofClusters;i++) {
    NsigmaTOFarray[i] = -9e3;
    mass2array[i] = -9e3;
    dxTOFarray[i] = -9e3;
    dzTOFarray[i] = -9e3;
    TOFarray[i] = -9e3;
    Lengtharray[i] = -9e3;
    expTOFarray[i] = -9e12;
    sigmaTOFarray[i] = -9e1;
    TOFchannelarray[i] = -9;
    //MCtofMismatcharray[i] = true;
    NmatchableTracks[i] = -9;
  }  
  
  TClonesArray *tofclArray = track->GetESDEvent()->GetESDTOFClusters();
  
  for(Int_t i=0;i < fNtofClusters;i++) {//loop over tof clusters
    
    AliESDTOFCluster *tofcl = (AliESDTOFCluster *) tofclArray->At(fTOFcluster[i]);
    for(Int_t j=0;j < tofcl->GetNMatchableTracks();j++) {//loop over all matchable tracks
      
      if(tofcl->GetTrackIndex(j) == track->GetID()) {//I filter for the track I am interested in (i.e. that I started from)
	
	Float_t exptime = tofcl->GetIntegratedTime(ParticleType, j);

	Float_t time = tofcl->GetTime();
	Float_t startime = fPIDResponse->GetTOFResponse().GetStartTime(track->P());
	Float_t tof = time - startime;
	Float_t sigmaTOF = fPIDResponse->GetTOFResponse().GetExpectedSigma(track->P(), tof, ParticleType);
	
	NsigmaTOFarray[i] = (tof-exptime)/sigmaTOF;
	dxTOFarray[i] = tofcl->GetDx(j);
	dzTOFarray[i] = tofcl->GetDz(j);
	NmatchableTracks[i] = tofcl->GetNMatchableTracks();
	  
	Float_t beta = -1.;
	Float_t mass2 = -1;
	Float_t length = tofcl->GetLength(j);
	Float_t ptot = track->GetTPCmomentum();
	
	beta = length / (TMath::C() * 1.e-10 * tof);
	if ((1 - beta * beta) > 0)
	  mass2 = ptot * ptot * (1. / (beta * beta) - 1.);
	
	mass2array[i] = mass2;

	TOFchannelarray[i] = tofcl->GetTOFchannel();
	TOFarray[i] = tof;
	Lengtharray[i] = length;
	expTOFarray[i] = exptime;
	sigmaTOFarray[i] = sigmaTOF;

	MCtofMismatcharray[i] = ( TMath::Abs(track->GetLabel()) != tofcl->GetLabel() );

      }
    }
  }

  //To fill the tree global variables:
  for(Int_t i=0;i<fNtofClusters;i++) {
    tNsigmaTOFarray[i] = NsigmaTOFarray[i];
    tmass2array[i] = mass2array[i];
    tdxTOFarray[i] = dxTOFarray[i];
    tdzTOFarray[i] = dzTOFarray[i];
    tTOFarray[i] = TOFarray[i];
    tLengtharray[i] = Lengtharray[i];
    texpTOFarray[i] = expTOFarray[i];
    tsigmaTOFarray[i] = sigmaTOFarray[i];
    tTOFchannelarray[i] = TOFchannelarray[i];
    tMCtofMismatcharray[i] = MCtofMismatcharray[i];
    tNmatchableTracks[i] = NmatchableTracks[i];
  }  
  
  return;
}
