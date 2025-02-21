/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                       *
 * Author: Baldo Sahlmueller, Friederike Bock                     *
 * Version 1.0                                 *
 *                                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims    *
 * about the suitability of this software for any purpose. It is      *
 * provided "as is" without express or implied warranty.               *
 **************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis on conversion photons + calo photons
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPDGCode.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAnalysisTaskGammaJetsPureMC.h"
#include "AliVParticle.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include <algorithm>
#include <array>
#include <vector>
#include <map>

#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>

ClassImp(AliAnalysisTaskGammaJetsPureMC)

//________________________________________________________________________
AliAnalysisTaskGammaJetsPureMC::AliAnalysisTaskGammaJetsPureMC(): AliAnalysisTaskSE(),
  fOutputContainer(nullptr),
  fHistNEvents(nullptr),
  fHistXSection(nullptr),
  fHistPtHard(nullptr),
  fHistJetPtY_Std(nullptr),
  fHistJetEta_Std(nullptr),
  fHistJetPhi_Std(nullptr),
  fHistJetPtY_StdNN(nullptr),
  fHistJetEta_StdNN(nullptr),
  fHistJetPhi_StdNN(nullptr),
  fHistJetPtY_Det(nullptr),
  fHistJetEta_Det(nullptr),
  fHistJetPhi_Det(nullptr),
  fHistJetPtY_DetNN(nullptr),
  fHistJetEta_DetNN(nullptr),
  fHistJetPhi_DetNN(nullptr),
  fHistResponse_Std_Det(nullptr),
  fHistUnMatched_Std_Det(nullptr),
  fHistRecUnMatched_Std_Det(nullptr),
  fHistMultiMatched_Std_Det(nullptr),
  fHistResponse_Std_DetNN(nullptr),
  fHistUnMatched_Std_DetNN(nullptr),
  fHistRecUnMatched_Std_DetNN(nullptr),
  fHistMultiMatched_Std_DetNN(nullptr),
  fHistJetPtPartPtVsPart(nullptr),
  fHistJetPtPartFragVsPart(nullptr),
  fHistJetPtPartPtVsPartLead(nullptr),
  fHistJetPtPartFragVsPartLead(nullptr),
  fHistPrimaryParticles(nullptr),
  fHistEnergyFracParticleCat(nullptr),
  fJetRadius(0.4),
  fJetMinE(5.),
  fJetAccEta(0.8),
  fJetParticleAcc(0.4),
  fJetParticleAccFF(1.2),
  fJetAlgorithm(fastjet::antikt_algorithm),
  fJetStrategy(fastjet::Best),
  fJetAreaType(fastjet::active_area),
  fJetRecombScheme(fastjet::BIpt_scheme),
  fJetGhostArea(0.01),
  fGhostEtaMax(1.5),
  fActiveAreaRepeats(1),
  fAreaType(fastjet::active_area),
  fEffiNeutral(nullptr),
  fEffiCharged(nullptr),
  fVecNonMeasureable({-1}),
  fDoJetEnergyShift(false),
  fJetEnergyShift(1.),
  fRand(0)
{

}

//________________________________________________________________________
AliAnalysisTaskGammaJetsPureMC::AliAnalysisTaskGammaJetsPureMC(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(nullptr),
  fHistNEvents(nullptr),
  fHistXSection(nullptr),
  fHistPtHard(nullptr),
  fHistJetPtY_Std(nullptr),
  fHistJetEta_Std(nullptr),
  fHistJetPhi_Std(nullptr),
  fHistJetPtY_StdNN(nullptr),
  fHistJetEta_StdNN(nullptr),
  fHistJetPhi_StdNN(nullptr),
  fHistJetPtY_Det(nullptr),
  fHistJetEta_Det(nullptr),
  fHistJetPhi_Det(nullptr),
  fHistJetPtY_DetNN(nullptr),
  fHistJetEta_DetNN(nullptr),
  fHistJetPhi_DetNN(nullptr),
  fHistResponse_Std_Det(nullptr),
  fHistUnMatched_Std_Det(nullptr),
  fHistRecUnMatched_Std_Det(nullptr),
  fHistMultiMatched_Std_Det(nullptr),
  fHistResponse_Std_DetNN(nullptr),
  fHistUnMatched_Std_DetNN(nullptr),
  fHistRecUnMatched_Std_DetNN(nullptr),
  fHistMultiMatched_Std_DetNN(nullptr),
  fHistJetPtPartPtVsPart(nullptr),
  fHistJetPtPartFragVsPart(nullptr),
  fHistJetPtPartPtVsPartLead(nullptr),
  fHistJetPtPartFragVsPartLead(nullptr),
  fHistPrimaryParticles(nullptr),
  fHistEnergyFracParticleCat(nullptr),
  fJetRadius(0.4),
  fJetMinE(5.),
  fJetAccEta(0.8),
  fJetParticleAcc(0.4),
  fJetParticleAccFF(1.2),
  fJetAlgorithm(fastjet::antikt_algorithm),
  fJetStrategy(fastjet::Best),
  fJetAreaType(fastjet::active_area),
  fJetRecombScheme(fastjet::BIpt_scheme),
  fJetGhostArea(0.01),
  fGhostEtaMax(1.5),
  fActiveAreaRepeats(1),
  fAreaType(fastjet::active_area),
  fEffiNeutral(nullptr),
  fEffiCharged(nullptr),
  fVecNonMeasureable({-1}),
  fDoJetEnergyShift(false),
  fJetEnergyShift(1.),
  fRand(0)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaJetsPureMC::~AliAnalysisTaskGammaJetsPureMC()
{

}

//________________________________________________________________________
void AliAnalysisTaskGammaJetsPureMC::UserCreateOutputObjects(){

  // Create histograms
  if(fOutputContainer != nullptr){
    delete fOutputContainer;
    fOutputContainer          = nullptr;
  }
  if(fOutputContainer == nullptr){
    fOutputContainer          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  fHistNEvents                		= new TH1F("NEvents", "NEvents", 3, -0.5, 2.5);
  fHistNEvents->Sumw2();
  fOutputContainer->Add(fHistNEvents);

  fHistXSection               		= new TH1D("XSection", "XSection", 1000000, 0, 1e4);

  //   SetLogBinningXTH1(fHistXSection);
  fHistXSection->Sumw2();
  fOutputContainer->Add(fHistXSection);

  fHistPtHard                 		= new TH1F("PtHard", "PtHard", 400, 0, 200);
  fHistPtHard->Sumw2();
  fOutputContainer->Add(fHistPtHard);

  std::vector<double> vecJetPt;
  double jetPt = 0.;
  double maxJetPt = 500.;
  double epsilon = 1e-20;
  while(jetPt < maxJetPt){
    vecJetPt.push_back(jetPt);
    if(jetPt < 5. - epsilon) jetPt +=1.;
    else if(jetPt < 50. - epsilon) jetPt +=5.;
    else if(jetPt < 100. - epsilon) jetPt +=10.;
    else if(jetPt < 200. - epsilon) jetPt +=50.;
    else if(jetPt < 500. - epsilon) jetPt +=100.;
    else {
      vecJetPt.push_back(maxJetPt);
      break;
    }
  }

  std::vector<double> vecPartPt;
  double partPt = 0.;
  double maxPartPt = 500.;
  while(partPt < maxPartPt){
    vecPartPt.push_back(partPt);
    if(partPt < 10. - epsilon) partPt +=0.1;
    else if(partPt < 50. - epsilon) partPt +=1.;
    else if(partPt < 100. - epsilon) partPt +=5.;
    else if(partPt < 200. - epsilon) partPt +=10.;
    else {
      vecPartPt.push_back(maxPartPt);
      break;
    }
  }
  
  fHistJetPtY_Std = new TH2F("JetPtY_Std", "JetPtY_Std", vecJetPt.size()-1, vecJetPt.data(), 100, -2, 2);
  fHistJetPtY_Std->Sumw2();
  fOutputContainer->Add(fHistJetPtY_Std);
  fHistJetEta_Std = new TH1D("JetEta_Std", "JetEta_Std", 100, -1., 1.);
  fHistJetEta_Std->Sumw2();
  fOutputContainer->Add(fHistJetEta_Std);
  fHistJetPhi_Std = new TH1D("JetPhi_Std", "JetPhi_Std", 100, 0., 6.5);
  fHistJetPhi_Std->Sumw2();
  fOutputContainer->Add(fHistJetPhi_Std);


  fHistJetPtY_StdNN = new TH2F("JetPtY_StdNN", "JetPtY_StdNN", vecJetPt.size()-1, vecJetPt.data(), 100, -2, 2);
  fHistJetPtY_StdNN->Sumw2();
  fOutputContainer->Add(fHistJetPtY_StdNN);
  fHistJetEta_StdNN = new TH1D("JetEta_StdNN", "JetEta_StdNN", 100, -1., 1.);
  fHistJetEta_StdNN->Sumw2();
  fOutputContainer->Add(fHistJetEta_StdNN);
  fHistJetPhi_StdNN = new TH1D("JetPhi_StdNN", "JetPhi_StdNN", 100, 0., 6.5);
  fHistJetPhi_StdNN->Sumw2();
  fOutputContainer->Add(fHistJetPhi_StdNN);


  fHistJetPtY_Det = new TH2F("JetPtY_Det", "JetPtY_Det", vecJetPt.size()-1, vecJetPt.data(), 100, -2, 2);
  fHistJetPtY_Det->Sumw2();
  fOutputContainer->Add(fHistJetPtY_Det);
  fHistJetEta_Det = new TH1D("JetEta_Det", "JetEta_Det", 100, -1., 1.);
  fHistJetEta_Det->Sumw2();
  fOutputContainer->Add(fHistJetEta_Det);
  fHistJetPhi_Det = new TH1D("JetPhi_Det", "JetPhi_Det", 100, 0., 6.5);
  fHistJetPhi_Det->Sumw2();
  fOutputContainer->Add(fHistJetPhi_Det);


  fHistJetPtY_DetNN = new TH2F("JetPtY_DetNN", "JetPtY_DetNN", vecJetPt.size()-1, vecJetPt.data(), 100, -2, 2);
  fHistJetPtY_DetNN->Sumw2();
  fOutputContainer->Add(fHistJetPtY_DetNN);
  fHistJetEta_DetNN = new TH1D("JetEta_DetNN", "JetEta_DetNN", 100, -1., 1.);
  fHistJetEta_DetNN->Sumw2();
  fOutputContainer->Add(fHistJetEta_DetNN);
  fHistJetPhi_DetNN = new TH1D("JetPhi_DetNN", "JetPhi_DetNN", 100, 0., 6.5);
  fHistJetPhi_DetNN->Sumw2();
  fOutputContainer->Add(fHistJetPhi_DetNN);


  fHistResponse_Std_Det = new TH2F("Response_Std_Det", "Response_Std_Det", vecJetPt.size()-1, vecJetPt.data(), vecJetPt.size()-1, vecJetPt.data());
  fHistResponse_Std_Det->Sumw2();
  fOutputContainer->Add(fHistResponse_Std_Det);

  fHistUnMatched_Std_Det = new TH1D("JetPt_Unmatched_Std_Det", "JetPt_Unmatched_Std_Det", vecJetPt.size()-1, vecJetPt.data());
  fHistUnMatched_Std_Det->Sumw2();
  fOutputContainer->Add(fHistUnMatched_Std_Det);

  fHistRecUnMatched_Std_Det = new TH1D("JetPt_RecUnmatched_Std_Det", "JetPt_RecUnmatched_Std_Det", vecJetPt.size()-1, vecJetPt.data());
  fHistRecUnMatched_Std_Det->Sumw2();
  fOutputContainer->Add(fHistRecUnMatched_Std_Det);

  fHistMultiMatched_Std_Det = new TH1D("JetPt_Multimatched_Std_Det", "JetPt_Multimatched_Std_Det", vecJetPt.size()-1, vecJetPt.data());
  fHistMultiMatched_Std_Det->Sumw2();
  fOutputContainer->Add(fHistMultiMatched_Std_Det);

  fHistResponse_Std_DetNN = new TH2F("Response_Std_DetNN", "Response_Std_DetNN", vecJetPt.size()-1, vecJetPt.data(), vecJetPt.size()-1, vecJetPt.data());
  fHistResponse_Std_DetNN->Sumw2();
  fOutputContainer->Add(fHistResponse_Std_DetNN);

  fHistUnMatched_Std_DetNN = new TH1D("JetPt_Unmatched_Std_DetNN", "JetPt_Unmatched_Std_DetNN", vecJetPt.size()-1, vecJetPt.data());
  fHistUnMatched_Std_DetNN->Sumw2();
  fOutputContainer->Add(fHistUnMatched_Std_DetNN);

  fHistRecUnMatched_Std_DetNN = new TH1D("JetPt_RecUnmatched_Std_DetNN", "JetPt_RecUnmatched_Std_DetNN", vecJetPt.size()-1, vecJetPt.data());
  fHistRecUnMatched_Std_DetNN->Sumw2();
  fOutputContainer->Add(fHistRecUnMatched_Std_DetNN);

  fHistMultiMatched_Std_DetNN = new TH1D("JetPt_Multimatched_Std_DetNN", "JetPt_Multimatched_Std_DetNN", vecJetPt.size()-1, vecJetPt.data());
  fHistMultiMatched_Std_DetNN->Sumw2();
  fOutputContainer->Add(fHistMultiMatched_Std_DetNN);

  fHistPrimaryParticles = new TH1D("PrimPart_Pt", "PrimPart_Pt", 1000, 0, 100);
  fHistPrimaryParticles->Sumw2();
  fOutputContainer->Add(fHistPrimaryParticles);

  std::vector<double> vec001Binning;
  for(int i = 0; i <= 100; ++i){
    vec001Binning.push_back(i*0.01);
  }
  std::vector<double> vecEquidistant = {-0.5, 0.5, 1.5, 2.5};

  fHistEnergyFracParticleCat = new TH3F("HistEnergyFracParticleCat", "HistEnergyFracParticleCat", vec001Binning.size()-1, vec001Binning.data(), vecEquidistant.size()-1, vecEquidistant.data(), vecJetPt.size()-1, vecJetPt.data());
  fHistEnergyFracParticleCat->Sumw2();
  fOutputContainer->Add(fHistEnergyFracParticleCat);


  std::vector<double> vecParticleSpec;
  for(int i = 0; i < 18; ++i){
    vecParticleSpec.push_back(i-0.5);
  }
  fHistJetPtPartPtVsPart = new TH3F("JetPtPartPtVsPart", "JetPtPartPtVsPart", vecParticleSpec.size()-1, vecParticleSpec.data(), vecPartPt.size()-1, vecPartPt.data(), vecJetPt.size()-1, vecJetPt.data());
  fHistJetPtPartPtVsPart->Sumw2();
  fOutputContainer->Add(fHistJetPtPartPtVsPart);

  fHistJetPtPartPtVsPartLead = new TH3F("JetPtPartPtVsPartLeading", "JetPtPartPtVsPartLeading", vecParticleSpec.size()-1, vecParticleSpec.data(), vecPartPt.size()-1, vecPartPt.data(), vecJetPt.size()-1, vecJetPt.data());
  fHistJetPtPartPtVsPartLead->Sumw2();
  fOutputContainer->Add(fHistJetPtPartPtVsPartLead);

  //---------------------------
  // Fragmentation Binning
  //---------------------------
  std::vector<double> fVecBinsFragment;
  double valZ = 0;
  for (int i = 0; i < 1000; ++i) {
    fVecBinsFragment.push_back(valZ);
    if (valZ < 0.1 - epsilon)
      valZ += 0.01;
    else if (valZ < 0.2 - epsilon)
      valZ += 0.02;
    else if (valZ < 1 - epsilon)
      valZ += 0.05;
    else if (valZ < 1.2 - epsilon)
      valZ += 0.1;
    else
      break;
  }

  fHistJetPtPartFragVsPart = new TH3F("JetPtPartFragVsPart", "JetPtPartFragVsPart", vecParticleSpec.size()-1, vecParticleSpec.data(), fVecBinsFragment.size()-1, fVecBinsFragment.data(), vecJetPt.size()-1, vecJetPt.data());
  fHistJetPtPartFragVsPart->Sumw2();
  fOutputContainer->Add(fHistJetPtPartFragVsPart);

  fHistJetPtPartFragVsPartLead = new TH3F("JetPtPartFragVsPartLeading", "JetPtPartFragVsPartLeading", vecParticleSpec.size()-1, vecParticleSpec.data(), fVecBinsFragment.size()-1, fVecBinsFragment.data(), vecJetPt.size()-1, vecJetPt.data());
  fHistJetPtPartFragVsPartLead->Sumw2();
  fOutputContainer->Add(fHistJetPtPartFragVsPartLead);


  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaJetsPureMC::UserExec(Option_t *)
{

  fInputEvent = InputEvent();
    // std::cout << "I found an Event" << std::endl;
  fMCEvent = MCEvent();
  if(!fMCEvent){
    printf("fMCEvent is null\n");
    return;
  }
  //   cout << "I found an MC header" << endl;
  const AliVVertex* primVtxMC  = fMCEvent->GetPrimaryVertex();
  if(!primVtxMC){
    AliGenEventHeader* mcEH = fMCEvent->GenEventHeader();
    if(!mcEH){
      return;
    }
  }
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  if (TMath::Abs(mcProdVtxZ) < 10 ){
    fHistNEvents->Fill(0);
  } else {
    fHistNEvents->Fill(1);
  }
  

  AliGenEventHeader* mcEH = fMCEvent->GenEventHeader();
  AliGenPythiaEventHeader *pyH  = dynamic_cast<AliGenPythiaEventHeader*>(mcEH);
  AliGenHijingEventHeader *hiH  = 0;
  AliGenDPMjetEventHeader *dpmH = 0;

  // it can be only one save some casts
  // assuming PYTHIA and HIJING are the most likely ones...
  if(!pyH){
    hiH = dynamic_cast<AliGenHijingEventHeader*>(mcEH);
    if(!hiH){
      dpmH = dynamic_cast<AliGenDPMjetEventHeader*>(mcEH);
    }
  }

  // fetch the trials on a event by event basis, not from pyxsec.root otherwise
  // we will get a problem when running on proof since Notify may be called
  // more than once per file
  // consider storing this information in the AOD output via AliAODHandler
  Float_t ntrials = 0;
  if (!pyH || !hiH || dpmH) {
    AliGenCocktailEventHeader *ccEH = dynamic_cast<AliGenCocktailEventHeader *>(mcEH);
    if (ccEH) {
      TList *genHeaders = ccEH->GetHeaders();
      for (int imch=0; imch<genHeaders->GetEntries(); imch++) {
        if(!pyH)pyH = dynamic_cast<AliGenPythiaEventHeader*>(genHeaders->At(imch));
        if(!hiH)hiH = dynamic_cast<AliGenHijingEventHeader*>(genHeaders->At(imch));
        if(!dpmH)dpmH = dynamic_cast<AliGenDPMjetEventHeader*>(genHeaders->At(imch));
      }
    }
  }

  // take the trials from the p+p event
  if(hiH)ntrials = hiH->Trials();
  if(dpmH)ntrials = dpmH->Trials();
  if(pyH)ntrials = pyH->Trials();
  if(ntrials)fHistNEvents->Fill(2,ntrials);

  Double_t xSection = 0;
  Double_t ptHard = 0;
  if (pyH) xSection = pyH->GetXsection();
  if (pyH) ptHard = pyH->GetPtHard();
  if (xSection) fHistXSection->Fill(xSection);
  if (ptHard) fHistPtHard->Fill(ptHard);

  ProcessJets();

  PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskGammaJetsPureMC::Terminate(const Option_t *)
{

  //fOutputContainer->Print(); // Will crash on GRID
}


//_________________________________________________________________________________
void AliAnalysisTaskGammaJetsPureMC::SetLogBinningXTH1(TH1* histoRebin){
  TAxis *axisafter = histoRebin->GetXaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaJetsPureMC::SetLogBinningXTH2(TH2* histoRebin){
  TAxis *axisafter = histoRebin->GetXaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}

//_________________________________________________________________________________
bool AliAnalysisTaskGammaJetsPureMC::AcceptParticle(AliVParticle* particle){
  double pt = particle->Pt();
  if(particle->Charge() == 0){
    if(fRand.Rndm() > fEffiNeutral->Eval(pt)) return false;
  } else {
    if(fRand.Rndm() > fEffiCharged->Eval(pt)) return false;
  }
  return true;
}

//_________________________________________________________________________________
bool AliAnalysisTaskGammaJetsPureMC::IsNonMeasureable(int pdgCode, int charge) const {
  if(fVecNonMeasureable.size() > 0 && fVecNonMeasureable[0] == -1){ // If vector contains -1 in first element, set all neutral particles that are not photons to non measureable
    if(charge == 0 && pdgCode != 22) return true;
  }
  if(std::find(fVecNonMeasureable.begin(), fVecNonMeasureable.end(), pdgCode) != fVecNonMeasureable.end()) {
    return true;
  }
  return false;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaJetsPureMC::ProcessJets(){

  double JetR2 = fJetRadius*fJetRadius;
  std::vector<fastjet::PseudoJet> vecParticles;
  std::vector<fastjet::PseudoJet> vecParticlesNN;
  std::vector<fastjet::PseudoJet> vecParticlesDet;
  std::vector<fastjet::PseudoJet> vecParticlesNNDet;
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    // fill primary histograms
    AliVParticle* particle     = nullptr;
    particle                    = (AliVParticle *)fMCEvent->GetTrack(i);
    if (!particle) continue;

    // only select stable particles for jet finder
    if(!particle->IsPhysicalPrimary()) continue; // Only consider primary particles
    bool isAccepted = AcceptParticle(particle);
    int pdg = std::abs(particle->PdgCode());
    bool isNonMeas = IsNonMeasureable(pdg, particle->Charge());

    fastjet::PseudoJet jetPart(particle->Px(), particle->Py(), particle->Pz(), particle->P());
    jetPart.set_user_index(i);
    vecParticles.push_back(jetPart);
    if(!isNonMeas) vecParticlesNN.push_back(jetPart);
    if(isAccepted){
      vecParticlesDet.push_back(jetPart);
      if(!isNonMeas) vecParticlesNNDet.push_back(jetPart);
    }
  }

  std::vector<fastjet::PseudoJet> fVecJets_Std;
  std::vector<fastjet::PseudoJet> fVecJets_StdNN;
  std::vector<fastjet::PseudoJet> fVecJets_Det;
  std::vector<fastjet::PseudoJet> fVecJets_DetNN;

  fastjet::Selector sel_jets = fastjet::SelectorEMin(fJetMinE) * fastjet::SelectorEtaRange(-fJetAccEta, fJetAccEta);
  fastjet::JetDefinition jet_def(fJetAlgorithm, fJetRadius, fJetRecombScheme, fJetStrategy);
  fastjet::GhostedAreaSpec ghostSpec(fGhostEtaMax, fActiveAreaRepeats, fJetGhostArea, 0.1, 1e-100);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fAreaType,ghostSpec);

  double maxJetPt = 0.;
  // standard jets
  {
    fastjet::ClusterSequenceArea clust_seq_full(vecParticles, jet_def, area_def);
    fVecJets_Std = sel_jets(clust_seq_full.inclusive_jets());

    for(const auto & jet : fVecJets_Std){
      fHistJetPtY_Std->Fill(jet.pt(), jet.eta());
      fHistJetEta_Std->Fill(jet.eta());
      fHistJetPhi_Std->Fill(jet.phi());
      if(maxJetPt < jet.pt()){
        maxJetPt = jet.pt();
      }
    }
  }
  // standard jets NN
  {
    fastjet::ClusterSequenceArea clust_seq_full(vecParticlesNN, jet_def, area_def);
    fVecJets_StdNN = sel_jets(clust_seq_full.inclusive_jets());

    for(const auto & jet : fVecJets_StdNN){
      fHistJetPtY_StdNN->Fill(jet.pt(), jet.eta());
      fHistJetEta_StdNN->Fill(jet.eta());
      fHistJetPhi_StdNN->Fill(jet.phi());
    }
  }

  // standard jets Det
  {
    fastjet::ClusterSequenceArea clust_seq_full(vecParticlesDet, jet_def, area_def);
    fVecJets_Det = sel_jets(clust_seq_full.inclusive_jets());

    for(const auto & jet : fVecJets_Det){
      fHistJetPtY_Det->Fill(jet.pt(), jet.eta());
      fHistJetEta_Det->Fill(jet.eta());
      fHistJetPhi_Det->Fill(jet.phi());
    }
  }

  // standard jets Det NN
  {
    fastjet::ClusterSequenceArea clust_seq_full(vecParticlesNNDet, jet_def, area_def);
    fVecJets_DetNN = sel_jets(clust_seq_full.inclusive_jets());

    for(const auto & jet : fVecJets_DetNN){
      fHistJetPtY_DetNN->Fill(jet.pt()*fJetEnergyShift, jet.eta());
      fHistJetEta_DetNN->Fill(jet.eta());
      fHistJetPhi_DetNN->Fill(jet.phi());
    }
  }


  // Response Matrices
  FillResponseMatrixAndEffi(fVecJets_Std, fVecJets_Det, fHistResponse_Std_Det, fHistUnMatched_Std_Det, fHistMultiMatched_Std_Det, fHistRecUnMatched_Std_Det);
  FillResponseMatrixAndEffi(fVecJets_Std, fVecJets_DetNN, fHistResponse_Std_DetNN, fHistUnMatched_Std_DetNN, fHistMultiMatched_Std_DetNN, fHistRecUnMatched_Std_DetNN, fJetEnergyShift);



  // check for particles inside of jets
  std::vector<int> vecLeadingPDG(fVecJets_Std.size());
  for(auto & i : vecLeadingPDG) i = 0;
  std::vector<int> vecLeadingMotherPDG(fVecJets_Std.size());
  for(auto & i : vecLeadingMotherPDG) i = 0;
  std::vector<double> vecLeadingE(fVecJets_Std.size());
  for(auto & i : vecLeadingE) i = 0.;
  std::array<double, 4> arrEnergyFracPart = {0., 0., 0., 0.};
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    // fill primary histograms
    AliVParticle* particle     = nullptr;
    particle                    = (AliVParticle *)fMCEvent->GetTrack(i);
    if (!particle) continue;
    if(!particle->IsPhysicalPrimary()) continue; // Only consider primary particles
    if(std::abs(particle->Eta()) < 1.0){
      fHistPrimaryParticles->Fill(particle->Pt());

      arrEnergyFracPart[0]+=particle->Pt();
      int pdg = std::abs(particle->PdgCode());
      bool isNonMeas = IsNonMeasureable(pdg, particle->Charge());
      if(isNonMeas) arrEnergyFracPart[1] += particle->Pt();
      else if(particle->Charge() != 0) arrEnergyFracPart[2] += particle->Pt();
      else arrEnergyFracPart[3] += particle->Pt();

    }

    int iMother = particle->GetMother();
    auto particleMother                    = (AliVParticle *)fMCEvent->GetTrack(iMother);
    int pdgMother = -1;
    if(particleMother){
      pdgMother = particleMother->PdgCode();
    }

    int jetindex = -1;
    double minR = 10000;
    for(size_t j = 0; j < fVecJets_Std.size(); ++j){
      double dEta = std::abs(fVecJets_Std[j].eta() - particle->Eta());
      double dPhi = std::abs(fVecJets_Std[j].phi() - particle->Phi());
      if(dPhi > TMath::Pi()) dPhi-= 2*TMath::Pi();
      double R2 = dEta*dEta + dPhi*dPhi;
      if(R2 < minR && R2 < JetR2){
        minR = R2;
        jetindex = j;
      }
    }
    double jetPt = jetindex == -1 ? 0 : fVecJets_Std[jetindex].pt();
    int partIndex = GetParticleIndex(std::abs(particle->PdgCode()), std::abs(pdgMother));
    fHistJetPtPartPtVsPart->Fill(partIndex, particle->Pt(), jetPt);
    if(jetPt > 0)fHistJetPtPartFragVsPart->Fill(partIndex, particle->Pt()/jetPt, jetPt);

    if(jetindex >= 0){
      if(particle->P() > vecLeadingE[jetindex]){
        vecLeadingPDG[jetindex] = particle->PdgCode();
        vecLeadingMotherPDG[jetindex] = pdgMother;
        vecLeadingE[jetindex] = particle->Pt();
      }
    }
  }

  if(arrEnergyFracPart[0] > 0){
    fHistEnergyFracParticleCat->Fill(arrEnergyFracPart[1]/arrEnergyFracPart[0], 0., maxJetPt);
    fHistEnergyFracParticleCat->Fill(arrEnergyFracPart[2]/arrEnergyFracPart[0], 1., maxJetPt);
    fHistEnergyFracParticleCat->Fill(arrEnergyFracPart[3]/arrEnergyFracPart[0], 2., maxJetPt);
  }


  for(size_t j = 0; j < fVecJets_Std.size(); ++j){
    int partIndex = GetParticleIndex(std::abs(vecLeadingPDG[j]), std::abs(vecLeadingMotherPDG[j]));
    fHistJetPtPartPtVsPartLead->Fill(partIndex, vecLeadingE[j], fVecJets_Std[j].pt());
    if(fVecJets_Std[j].pt() > 0)fHistJetPtPartFragVsPartLead->Fill(partIndex, vecLeadingE[j]/fVecJets_Std[j].pt(), fVecJets_Std[j].pt());
  }

}


void AliAnalysisTaskGammaJetsPureMC::FillResponseMatrixAndEffi(std::vector<fastjet::PseudoJet> vecTrueJet, std::vector<fastjet::PseudoJet> vecRecJet, TH2F* hResp, TH1D* hUnMatched, TH1D* hMultiMatched, TH1D* hRecUnMatched, double energyShiftRec){

  
  std::vector<int> vecMatched;
  for(size_t irec = 0; irec < vecRecJet.size(); ++irec){
    double maxR = 10000;
    int index = -1;
    for(size_t itrue = 0; itrue < vecTrueJet.size(); ++itrue){
      double dEta = vecRecJet[irec].eta() - vecTrueJet[itrue].eta();
      double dPhi = vecRecJet[irec].phi() - vecTrueJet[itrue].phi();
      if(dPhi > TMath::Pi()) dPhi-= 2*TMath::Pi();
      double R2 = dEta*dEta + dPhi*dPhi;
      if(R2 < maxR && R2 < 0.5*fJetRadius*fJetRadius){
        // std::cout << "R2 " << R2 << "  maxR " << maxR << "  irec " << irec << "  itrue " << itrue << std::endl;
        maxR = R2;
        index = itrue;
      }
    }
    if(index != -1){
      // std::cout << "maxR " << maxR << "  vecRecJet[irec].pt() " << vecRecJet[irec].pt() << "  vecTrueJet[index].pt() " << vecTrueJet[index].pt() << "  irec " << irec << "  index " << index << std::endl;
      
      vecMatched.push_back(index);
      hResp->Fill(vecRecJet[irec].pt()*energyShiftRec, vecTrueJet[index].pt());
    } else {
      hRecUnMatched->Fill(vecRecJet[irec].pt()*energyShiftRec);
    }
  }
  for(size_t itrue = 0; itrue < vecTrueJet.size(); ++itrue){
    if(std::find(vecMatched.begin(), vecMatched.end(), itrue) == vecMatched.end()) {
      hUnMatched->Fill(vecTrueJet[itrue].pt());
    }
  }

  for(size_t i = 0; i < vecMatched.size(); ++i){
    for(size_t j = i+1; j < vecMatched.size(); ++j){
      if(vecMatched[i] == vecMatched[j]){
        hMultiMatched->Fill(vecTrueJet[vecMatched[i]].pt());
      }
    }
  }


}

//________________________________________________________________________
int AliAnalysisTaskGammaJetsPureMC::GetParticleIndex(const int pdgcode, const int motherpdg) const {
  if(pdgcode == 22){ // gamma
    if(motherpdg > 100 &&  motherpdg < 10000){ // decay photon
      return 0;
    } else { // direct photon
      return 1;
    }
  } else if(pdgcode == 211){ // pi+-
    return 2;
  } else if(pdgcode == 321){ // Kaon
    return 3;
  } else if(pdgcode == 130){ // K0s
    return 4;
  } else if(pdgcode == 310){ // K0l
    return 5;
  } else if(pdgcode == 2212){ // proton
    return 6;
  } else if(pdgcode == 2112){ // neutron
    return 7;
  } else if(pdgcode == 3122){  // Lambda
    return 8;
  } else if(pdgcode == 11){ // electron
    return 9;
  } else if(pdgcode == 13 ){ // Muons
    return 10;
  } else if(pdgcode == 3112 || pdgcode == 3222 ){ // charged Sigma baryons
    return 11;
  } else if(pdgcode == 3212 ){ // Sigma 0
    return 12;
  } else if(pdgcode == 3312){ // charged Cascade
    return 13;
  } else if(pdgcode == 3322){ // neutral Cascade
    return 14;
  } else if(pdgcode == 12 || pdgcode == 14 || pdgcode == 16 ){ // Neutrinos
    return 15;
  
  }
  return 16;
}


//________________________________________________________________________
void AliAnalysisTaskGammaJetsPureMC::SetEfficiency(TString strNeutral, TString strCharged){
  fEffiNeutral = new TF1("fEffiNeutral", strNeutral, 0, 10000);
  fEffiCharged = new TF1("fEffiCharged", strCharged, 0, 10000);
}

//________________________________________________________________________
void AliAnalysisTaskGammaJetsPureMC::SetParticlesNonMeas(TString strPart){
  TObjArray* tokens = strPart.Tokenize(",");
  if (!tokens) return; // Return if no tokens
  fVecNonMeasureable.clear();

  for (int i = 0; i < tokens->GetEntries(); ++i) {
      TObjString* objStr = dynamic_cast<TObjString*>(tokens->At(i));
      if (objStr) {
          TString token = objStr->GetString();
          token = token.Strip(TString::kBoth); // Trim whitespace
          fVecNonMeasureable.push_back(atoi(token.Data())); // Convert to int
      }
  }
  std::cout << "Printing vector of non measurable particles: " << std::endl;
  for(const auto & i : fVecNonMeasureable){
    std::cout << i << ", ";
  }
  std::cout << std::endl;

  delete tokens;
}

