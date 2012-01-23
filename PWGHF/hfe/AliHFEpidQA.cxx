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
//
// Class for PID QA
// Several studies done on clean samples of electrons, pions and kaons
// coming from V0 PID
// Compatible with both ESDs and AODs
//
// Autors:
//    Matus Kalisky <matus.kalisky@cern.ch>
//    Markus Heide <mheide@uni-muenster.de>
//    Markus Fasel <M.Fasel@gsi.de>
//


#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TString.h>
#include <TMultiLayerPerceptron.h>
#include <TFile.h>

#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVParticle.h"


#include "AliHFEcollection.h"
#include "AliHFEpidQA.h"
#include "AliHFEV0info.h"
#include "AliHFEV0pid.h"
#include "AliHFEV0pidMC.h"
#include "AliHFEV0cuts.h"
#include "AliHFEtrdPIDqa.h"


ClassImp(AliHFEpidQA)

  //__________________________________________
  AliHFEpidQA::AliHFEpidQA():
    fEvent(NULL),
    fMC(NULL), 
    fV0pid(NULL), 
    fV0pidMC(NULL), 
    fTRDpidQA(NULL), 
    fOutput(NULL), 
    fESDpid(NULL),
    fNNref(NULL),
    fTRDTotalChargeInSlice0(kFALSE)
{
  //
  // Default constructor
  //
  for(Int_t mom = 0; mom < 11; mom++){
    fNet[mom] = NULL;
  }
}

//__________________________________________
AliHFEpidQA::AliHFEpidQA(const AliHFEpidQA &ref):
  TObject(ref),
  fEvent(NULL), 
  fMC(NULL),
  fV0pid(NULL),
  fV0pidMC(NULL), 
  fTRDpidQA(NULL),
  fOutput(NULL), 
  fESDpid(NULL),
  fNNref(NULL),
  fTRDTotalChargeInSlice0(ref.fTRDTotalChargeInSlice0)
{
  //
  // Copy constructor
  //
  for(Int_t mom = 0; mom < 11; mom++){
    fNet[mom] = NULL;
  }
  ref.Copy(*this);
}

//__________________________________________
AliHFEpidQA &AliHFEpidQA::operator=(const AliHFEpidQA &ref){
  //
  // Assignment operator
  //
  if(this != &ref)
    ref.Copy(*this);
  return *this;
}

//__________________________________________
AliHFEpidQA::~AliHFEpidQA(){
  //
  // Destructor
  //

  // the pointers bellow are not dfeleted to prevend double-deleting of some of the content
  // these pointers are defined only once during the program call and should not cause a problem, 
  // but cleaner solution is necessary.
  //if(fV0pid) delete fV0pid;
  //if(fV0pidMC) delete fV0pidMC;
  //if(fOutput) delete fOutput;

  //  if(fTRDpidResponse) delete fTRDpidResponse; 
}

//__________________________________________
void AliHFEpidQA::Copy(TObject &o) const {
  //
  // Copy function
  //
  
  TObject::Copy(o);

  AliHFEpidQA &target = dynamic_cast<AliHFEpidQA &>(o);
  target.fMC = fMC;
  target.fTRDTotalChargeInSlice0 = fTRDTotalChargeInSlice0;

  if(target.fESDpid) delete target.fESDpid;
  target.fESDpid = fESDpid;
  if(target.fV0pid) delete target.fV0pid;
  if(fV0pid)
    target.fV0pid = dynamic_cast<AliHFEV0pid *>(fV0pid->Clone());
  else
    target.fV0pid = NULL;
  if(target.fV0pidMC) delete target.fV0pidMC;
  if(fV0pidMC) 
    target.fV0pidMC = dynamic_cast<AliHFEV0pidMC *>(fV0pidMC->Clone());
  else
    target.fV0pidMC = NULL;
  if(target.fTRDpidQA) delete target.fTRDpidQA;
  if(fTRDpidQA)
    target.fTRDpidQA = dynamic_cast<AliHFEtrdPIDqa *>(fTRDpidQA->Clone());
  else
    target.fTRDpidQA = NULL;
  if(target.fOutput) delete target.fOutput;
  if(fOutput)
    target.fOutput = dynamic_cast<AliHFEcollection *>(fOutput->Clone());
}

//__________________________________________
void AliHFEpidQA::Init(){
  //
  // Prepare task output
  //

  // load networks
  if(fNNref){
    for(Int_t mom = 0; mom < 11; mom++){                      
      fNet[mom] = (TMultiLayerPerceptron*) fNNref->Get(Form("NN_Mom%d", mom));
      if(!fNet[mom]){
	AliError(Form("No reference network for momentum bin %d!", mom));
      }
    }
  }
  
  fV0pid = new AliHFEV0pid("fV0pid");
  if(HasV0pidQA()) fV0pid->InitQA();
  fV0pidMC = new AliHFEV0pidMC();
  fV0pidMC->Init();

  fTRDpidQA = new AliHFEtrdPIDqa;
  if(fTRDTotalChargeInSlice0) fTRDpidQA->SetTotalChargeInSlice0();
  fTRDpidQA->Init();

  fOutput = new AliHFEcollection("pidQA", "PID QA output");

  const char *name[4] = {"Electron", "PionK0", "PionL", "Proton"};
  const char *title[4] = {"Electron", "K0 Pion", "Lambda Pion", "Proton"};
  const char *det[4] = {"ITS", "TPC", "TRD", "TOF"};
  const int effs[6] = {70, 75, 80, 85, 90, 95};
  for(Int_t i = 0; i < 4; i++){
    fOutput->CreateTH2F(Form("purity%s", name[i]), Form("%s Purity", title[i]), 2, -0.5, 1.5, 20, 0.1, 10, 1);

    for(Int_t idet = 0; idet < 4; idet++){
      // create all the histograms which all the detectors have in common
      if(idet != 2){ // no nSigma histogram for TRD
        fOutput->CreateTH2F(Form("h%s_nSigma_%s", det[idet], name[i]), Form("%s number of sigmas for %ss; p (GeV/c); number of sigmas", det[idet], title[i]), 20, 0.1, 10, 100, -7, 7, 0);
      }
      fOutput->CreateTH2F(Form("h%s_PID_p_%s", det[idet], name[i]), Form("%s PID for %ss; p (GeV/c); ITS PID", det[idet], title[i]), 100, 0.1, 10, 5, -0.5, 4.5, 0);
      fOutput->CreateTH2F(Form("h%s_El_like_%s", det[idet], name[i]), Form("%s Electron Likelihoods for %ss; p (GeV/c); likelihood", det[idet], title[i]), 25, 0.1, 10, 1000, 0., 1., 0);
      fOutput->CreateTH2F(Form("h%s_El_like_MC_%s", det[idet], name[i]), Form("%s Electron Likelihoods for MC %ss; p (GeV/c); likelihood", det[idet], title[i]), 25, 0.1, 10, 1000, 0., 1., 0);
    }

    //
    // ITS pid response
    //
    fOutput->CreateTH2F(Form("hITS_Signal_%s", name[i]), Form("ITS Signal org. for %ss", title[i]), 40, 0.1, 10, 400, 0, 1000, 0);
    fOutput->CreateTH2Fvector1(5, Form("hITS_dEdx_%s", name[i]), Form("ITS dEdx for %ss; p (GeV/c); dEdx (a.u.)", title[i]), 40, 0.1, 10, 400, 0, 1000, 0);
    
    //
    // TPC pid response
    //
    fOutput->CreateTH2F(Form("hTPC_dEdx_%s", name[i]), Form("TPC dEdx for %ss; p (GeV/c); dEdx (a.u.)", title[i]), 20, 0.1, 10, 200, 0, 200, 0);

    //
    // TRD pid response 
    //
    fOutput->CreateTH2F(Form("hTRD_trk_%s", name[i]), Form("%s PID tracklets; p (GeV/c); N TRD tracklets", title[i]), 100, 0.1, 10, 7, -0.5, 6.5, 0);
    // number of the non 0 slices per tracklet
    fOutput->CreateTH2F(Form("hTRD_Nslc_%s", name[i]), Form("%s PID slices > 0; p (GeV/c); N slc > 0", title[i]), 100, 0.1, 10, 9, -0.5, 8.5, 0);
    // location of the slices > 0 - where are the emtpy slices located ?
    fOutput->CreateTH2F(Form("hTRD_slc_%s", name[i]), Form("%s PID slices > 0 position; p (GeV/c); slice", title[i]), 100, 0.1, 10, 9, -0.5, 8.5, 0);
    fOutput->CreateTH2F(Form("hTRD_cls_%s", name[i]), Form("TRD clusters for %s candidates; p (GeV/c); N cls", title[i]), 25, 0.1, 10, 1000, 0, 1000);
    fOutput->CreateTH2F(Form("hTRD_dEdx_%s", name[i]), Form("TRD dEdx (trk) for %s candidates; p (GeV/c); tracklet dEdx (a.u.)", title[i]), 25, 0.1, 10, 1000, 0, 100000, 0);

    for(Int_t itl = 4; itl < 6; itl++){
      for(Int_t ieff = 0; ieff < 6; ieff++){
        fOutput->CreateTH2F(Form("hTRD_likeSel_%s_%dtls_%deff", name[i], itl, effs[ieff]), Form(" %s selected as electrons  with %d tracklets and %f electron efficiency", title[i], itl, static_cast<Double_t>(effs[ieff])/100.), 44, 0.1, 20., 100, 0., 1.);
        fOutput->CreateTH2F(Form("hTRD_likeRej_%s_%dtls_%deff", name[i], itl, effs[ieff]), Form(" %s rejected as electrons  with %d tracklets and %f electron efficiency", title[i], itl, static_cast<Double_t>(effs[ieff])/100.), 44, 0.1, 20., 100, 0., 1.);
      }
    }
    //
    // TOF pid response
    //
    fOutput->CreateTH2F(Form("hTOF_beta_%s", name[i]), Form("TOF beta for %s; p (GeV/c); #beta", title[i]),  50, 0.1, 5, 350, 0.4, 1.1, 0);
  }//.. loop over identified particle species

  // Global histograms
  fOutput->CreateTH1F("hITS_dEdx_nSamples", "ITS - number of non 0 dEdx samples", 4, 0.5, 4.5);
  fOutput->CreateTH2F("hTPC_PID", "TPC pid all tracks; tpc pid probability; species",100, 0, 1, 5, -0.5, 4.5 );
  fOutput->CreateTH2F("hTOF_PID", "TOF pid all tracks; tof pid probability; species",100, 0, 1,5,  -0.5, 4.5 );
  fOutput->CreateTH2F("hTOF_beta_all", "TOF beta for all nice single tracks; p (GeV/c); #beta", 100, 0.1, 10, 480, 0, 1.2, 0);
  fOutput->CreateTH2F("hTOF_qa_TmT0mT", "TOF (t - t0 - t[pion]) qa verus run",  10000, 114000, 124000, 200, -200, 200);
  fOutput->CreateTH2F("hTOF_qa_length", "TOF track length verus run",  10000, 114000, 112400, 200, 0, 1000);
  
  //
  // debug histograms
  //
  fOutput->CreateTH2F("hITS_kFlags", "ITS flags; flag; V0 candidates", 5, 0.5, 5.5, 5, -0.5, 4.5);  
  fOutput->CreateTH2F("hTPC_kFlags", "TPC flags; flag; V0 candidates", 5, 0.5, 5.5, 5, -0.5, 4.5);
  fOutput->CreateTH2F("hTRD_kFlags", "TRD flags; flag; V0 candidates", 5, 0.5, 5.5, 5, -0.5, 4.5);
  fOutput->CreateTH2F("hTOF_kFlags", "TOF flags; flag; V0 candidates", 5, 0.5, 5.5, 5, -0.5, 4.5);
  

  //
  // event based histograms
  //
  Int_t nT0[2] = {10000, 100};
  Double_t minT0[2] = {114500, -1000};
  Double_t maxT0[2] = {124500, 1000};
  fOutput->CreateTHnSparse("hEvent_T0", "T0 as a function of run number; run number; T0 (ns)", 2, nT0, minT0, maxT0);

  //
  // test the tender V0 supply
  //
  fOutput->CreateTH2F("h_tender_check_01", "tender -vs- HFEpidQA; HFEpidQA V0 candiadates; tender V0 candidates", 4, -1.5, 2.5, 4, -1.5, 2.5);
  fOutput->CreateTH2Fvector1(3, "h_tender_check_02", "tender -vs- HFEpidQA per type; AliHFEpidQA V0 ; tender V0", 4, -1.5, 2.5, 4, -1.5, 2.5);


  //
  // THnSpasre objects
  //

  // Create Illumination Plot
  // bins: species, pt, eta, phi, TPC status, TRD status
  {
    Int_t nbins[6] = {4, 40, 16, 72, 2, 2};
    Double_t min[6] = { 0, 0.1, -0.8, 0., 0., 0.};
    Double_t max[6] = { 4., 20., 0.8, 2*TMath::Pi(), 2., 2.};
    fOutput->CreateTHnSparse("hIllumination", "Illumination", 6, nbins, min, max);
    fOutput->BinLogAxis("hIllumination", 1);
  }

  // TPC THnSparse
  // bins: species, pt, n TPC clusters., TPC electron likelihood, TPC n sigmas, TPC signal
  {
    Int_t nbins[6] = { 5, 40, 20, 100, 100, 200};
    Double_t min[6] = { -0.5, 0.1, 0., 0., -5., 0.};
    Double_t max[6] = { 4.5, 20., 200, 1., 5., 120.};
    TString htitle = "TPC info sparse; VO identified species; p (GeV/c); n TPC clusters; TPC N sigma; TPC signal";
    fOutput->CreateTHnSparse("hTPCclusters",htitle,6, nbins, min, max);
    fOutput->BinLogAxis("hTPCclusters", 1);
  }
  // TRD THnSparse - entries per tracklet
  // species, p, tracklet position, number of PID tracklets, number of slices (non 0), number of clusters, electron likelihood, 

  {
    Int_t nbins[7] = {5, 20, 6, 7, 9, 45, 100};
    Double_t min[7] = {-0.5, 0.5, -0.5, -0.5, -0.5, -0.5, 0.};
    Double_t max[7] = {4.5, 10, 5.5, 6.5, 8.5, 179.5, 1.};
    TString htitle = "TRD tracklets; VO identified species; p (GeV/c); tracklet position; No. PID tacklets; No. slices; No. clusters; El. likelihood";
    fOutput->CreateTHnSparse("hTRDtracklets",htitle,7, nbins, min, max);
    fOutput->BinLogAxis("hTRDtracklets", 1);
  }
  

}
//__________________________________________
void AliHFEpidQA::Process(){
  //
  // Run PID QA
  //

  if(!fV0pid){
    AliError("V0pid not available! Forgotten to initialize?");
    return;
  }
  if(!fESDpid){
    AliError("fESDpid not initialized, I am leaving this!");
    return;
  }

  // to be udpated to AOD save mdoe
  if(!fEvent){
    AliError("AliVEvent not available, returning");
  }


  if(fMC) fV0pidMC->SetMCEvent(fMC);
  if(fMC) fV0pid->SetMCEvent(fMC);

  fV0pid->Process(fEvent);
  TObjArray *hfeelectrons = fV0pid->GetListOfElectrons();
  TObjArray *hfepionsK0 = fV0pid->GetListOfPionsK0();
  TObjArray *hfepionsL = fV0pid->GetListOfPionsL();
  TObjArray *hfeprotons = fV0pid->GetListOfProtons();

  // Get Track list for normal purpose
  TObjArray *electrons = MakeTrackList(hfeelectrons);
  TObjArray *pionsK0 = MakeTrackList(hfepionsK0);
  TObjArray *pionsL = MakeTrackList(hfepionsL);
  TObjArray *protons = MakeTrackList(hfeprotons);
  TObjArray *cleanElectrons = MakeCleanListElectrons(hfeelectrons);

  if(fMC){
    fV0pidMC->Process(electrons, AliHFEV0cuts::kRecoElectron);
    fV0pidMC->Process(pionsK0, AliHFEV0cuts::kRecoPionK0);
    fV0pidMC->Process(pionsL, AliHFEV0cuts::kRecoPionL);
    fV0pidMC->Process(protons, AliHFEV0cuts::kRecoProton);
  }

  AliDebug(2, Form("Number of Electrons      : %d", electrons->GetEntries()));
  AliDebug(2, Form("Number of K0 Pions       : %d", pionsK0->GetEntries()));
  AliDebug(2, Form("Number of Lambda Pions   : %d", pionsL->GetEntries()));
  AliDebug(2, Form("Number of Protons        : %d", protons->GetEntries()));
  if(fMC){
    AliDebug(2, "MC Information available. Doing Purity checks...");
    // Calculate the purity of the clean samples using MC 
    MakePurity(electrons, AliHFEV0cuts::kRecoElectron);
    MakePurity(pionsK0,  AliHFEV0cuts::kRecoPionK0);
    MakePurity(pionsL,  AliHFEV0cuts::kRecoPionL);
    MakePurity(protons,  AliHFEV0cuts::kRecoProton);
  }

  // some event wise checks
  CheckEvent();

  // Make Illumination Plot
  FillIllumination(electrons, AliHFEV0cuts::kRecoElectron);
  FillIllumination(pionsK0, AliHFEV0cuts::kRecoPionK0);
  FillIllumination(pionsL, AliHFEV0cuts::kRecoPionL);
  FillIllumination(protons, AliHFEV0cuts::kRecoProton);

  // Now we can do studies on the PID itself
  // For TRD use the TRD PID QA object
  fTRDpidQA->ProcessTracks(cleanElectrons, AliPID::kElectron);
  fTRDpidQA->ProcessTracks(pionsK0, AliPID::kPion);
  fTRDpidQA->ProcessTracks(pionsL, AliPID::kPion);
  fTRDpidQA->ProcessTracks(protons, AliPID::kProton);

  // Monitor TRD PID Response
  /*TestTRDResponse(cleanElectrons, AliHFEV0cuts::kRecoElectron);
  TestTRDResponse(pionsK0, AliHFEV0cuts::kRecoPionK0);
  TestTRDResponse(pionsL, AliHFEV0cuts::kRecoPionL);
  TestTRDResponse(protons, AliHFEV0cuts::kRecoProton);*/

  FillElectronLikelihoods(electrons,  AliHFEV0cuts::kRecoElectron); 
  FillElectronLikelihoods(pionsK0,  AliHFEV0cuts::kRecoPionK0); 
  FillElectronLikelihoods(pionsL,  AliHFEV0cuts::kRecoPionL); 
  FillElectronLikelihoods(protons,  AliHFEV0cuts::kRecoProton); 
  
  FillPIDresponse(electrons, AliHFEV0cuts::kRecoElectron);
  FillPIDresponse(pionsK0, AliHFEV0cuts::kRecoPionK0);
  FillPIDresponse(pionsL, AliHFEV0cuts::kRecoPionL);
  FillPIDresponse(protons, AliHFEV0cuts::kRecoProton);

  // check the tender V0s
  CheckTenderV0pid(electrons, AliHFEV0cuts::kRecoElectron);
  CheckTenderV0pid(pionsK0, AliHFEV0cuts::kRecoPionK0);
  CheckTenderV0pid(pionsL, AliHFEV0cuts::kRecoPionL);
  CheckTenderV0pid(protons, AliHFEV0cuts::kRecoProton);

  // Analysis done, flush the containers
  fV0pid->Flush();

  delete electrons;
  delete pionsL;
  delete pionsK0;
  delete protons;
  delete cleanElectrons;
}

//__________________________________________
void AliHFEpidQA::FillIllumination(const TObjArray * const tracks, Int_t species){
  //
  // Fill Illumination Plot
  //
  THnSparseF *hIllumination = dynamic_cast<THnSparseF *>(fOutput->Get("hIllumination"));
  if(!hIllumination) return;

  Double_t quantities[6]; memset(quantities, 0, sizeof(Double_t) *6);
  TIter trackIter(tracks);

  quantities[0] = species;
  TObject *o = NULL; AliESDtrack *esdtrack = NULL;
  while((o = trackIter())){
    if(!TString(o->IsA()->GetName()).CompareTo("AliESDtrack")){
      // work on local copy in order to not spoil others
      esdtrack = new AliESDtrack(*(static_cast<AliESDtrack *>(o)));  
      if(!esdtrack) continue;
    } else if(!TString(o->IsA()->GetName()).CompareTo("AliAODrack")){
      // Bad hack: Fill ESD track with AOD information
      esdtrack = new AliESDtrack(static_cast<AliAODTrack *>(o));
      if(!esdtrack) continue;
    } else {
      // Non usable
      continue;
    }

    // Propagate to the entrance of the TRD
    esdtrack->PropagateTo(300, fEvent->GetMagneticField());
    quantities[1] = esdtrack->Pt();
    quantities[2] = esdtrack->Eta();
    quantities[3] = esdtrack->Phi();
    quantities[4] = (esdtrack->GetStatus() & AliESDtrack::kTPCrefit) ? 1 : 0;
    quantities[5] = (esdtrack->GetStatus() & AliESDtrack::kTRDout) ? 1. : 0.;
    hIllumination->Fill(quantities);

    delete esdtrack;
  }
}
//__________________________________________
void AliHFEpidQA::FillTPCinfo(AliESDtrack *const esdtrack, Int_t species){
  //
  // Fill TPC Cluster Plots
  //
  THnSparseF *hTPCclusters = dynamic_cast<THnSparseF *>(fOutput->Get("hTPCclusters"));
  if(!hTPCclusters) return;

  Double_t quantities[6]; memset(quantities, 0, sizeof(Double_t) *6);
  
  Double_t pidProbs[5];
  const Int_t typePID[5] = {0, 2, 2, 3, 4};


  quantities[0] = species;

  
  esdtrack->GetTPCpid(pidProbs);   
    
  quantities[1] = esdtrack->P();
  quantities[2] = esdtrack->GetTPCNcls();
  quantities[3] = pidProbs[0];
  quantities[4] = fESDpid->NumberOfSigmasTPC(esdtrack,(AliPID::EParticleType)typePID[species]);
  quantities[5] = esdtrack->GetTPCsignal();
  hTPCclusters->Fill(quantities);

}

//__________________________________________
void AliHFEpidQA::TestTRDResponse(const TObjArray * const tracks, Int_t species){
  //
  // Test PID Response function of the TRD
  //
  Int_t effInt[6] = {70, 75, 80, 85, 90, 95};
  const char *sname[5] = {"Electron", "PionK0", "PionL", "Kaon", "Proton"};
  AliESDtrack *track = NULL;
  TIter trackIter(tracks);
  Float_t momenta[6], p;
  Int_t nmomenta;
  Double_t probs[5], likeEl;
  while((track = static_cast<AliESDtrack *>(trackIter()))){
    if(track->GetTRDntrackletsPID() < 4) continue;

    // calculate momentum at TRD level
    memset(momenta, 0, sizeof(Float_t) * 6); nmomenta = 0;
    for(Int_t ily = 0; ily < 6; ily++){
      if(track->GetTRDmomentum(ily) > 0.01) momenta[nmomenta++] = track->GetTRDmomentum(ily);
    }
    p = TMath::Mean(nmomenta, momenta);

    // Get Electron likelihood
    track->GetTRDpid(probs);
    likeEl = probs[0]/(probs[0] + probs[2]);

    for(Int_t ieff = 0; ieff < 6; ieff++){
      if(fESDpid->IdentifiedAsElectronTRD(track, static_cast<Double_t>(effInt[ieff])/100.)) fOutput->Fill(Form("hTRD_likeSel_%s_%dtls_%deff", sname[species], static_cast<Int_t>(track->GetTRDntrackletsPID()), effInt[ieff]), p, likeEl);
      else fOutput->Fill(Form("hTRD_likeRej_%s_%dtls_%deff", sname[species], static_cast<Int_t>(track->GetTRDntrackletsPID()), effInt[ieff]), p, likeEl);
    }
  }
}

//__________________________________________
void AliHFEpidQA::MakePurity(const TObjArray *tracks, Int_t species){
  //
  // Fill the QA histos for a given species
  //
  if(!fMC) return;
  AliDebug(3, Form("Doing Purity checks for species %d", species));
  Int_t pdg = GetPDG(species);
  TString hname;

  switch(species){
  case  AliHFEV0cuts::kRecoElectron:
    hname = "purityElectron";
    break;
  case  AliHFEV0cuts::kRecoPionK0:
    hname = "purityPionK0";
    break;
  case  AliHFEV0cuts::kRecoPionL:
    hname = "purityPionL";
    break;
  case  AliHFEV0cuts::kRecoProton:
    hname = "purityProton";
    break;
  default:  // non investigated species
    AliDebug(3, "Species not investigated");
    return;
  }  

  AliDebug(3, Form("Number of tracks: %d", tracks->GetEntries()));
  TIter trackIter(tracks);
  AliVParticle *recTrack = NULL, *mcTrack = NULL;
  while((recTrack = dynamic_cast<AliVParticle *>(trackIter()))){
    Int_t label = recTrack->GetLabel();
    AliDebug(4, Form("MC Label %d", label));
    mcTrack =fMC->GetTrack(TMath::Abs(label));
    if(!mcTrack){
      AliDebug(4, "MC track not available");
      continue; // we don't know
    }

    // Get the pdg code
    Int_t trackPdg = 0;
    if(!TString(mcTrack->IsA()->GetName()).CompareTo("AliMCParticle")){
      // case ESD
      AliMCParticle *mcp = dynamic_cast<AliMCParticle *>(mcTrack);
      if(!mcp) continue;
      trackPdg = TMath::Abs(mcp->Particle()->GetPdgCode());
    } else {
      // case AOD
      AliAODMCParticle *aodmcp = dynamic_cast<AliAODMCParticle *>(mcTrack);
      if(!aodmcp) continue;
      trackPdg = TMath::Abs(aodmcp->GetPdgCode());
    }
    if(trackPdg == pdg)    // Correct identification
      {
	fOutput->Fill(hname, 0., recTrack->Pt());
      }
    else  // Wrong identification
      fOutput->Fill(hname, 1., recTrack->Pt());
  }
}

//__________________________________________
void AliHFEpidQA::FillElectronLikelihoods(const TObjArray * const particles, Int_t species){
  //
  // Fill electron Likelihoods for the ITS, TPC and TOF
  // Required for the calculation of the electron efficiency, 
  // pion and proton efficiency and the thresholds
  //
  Long_t status = 0;
  const TString detname[4] = {"ITS", "TPC", "TRD", "TOF"};
  TString specname;

  switch(species){
  case  AliHFEV0cuts::kRecoElectron:
    specname = "Electron";
    break;
  case  AliHFEV0cuts::kRecoPionK0:
    specname = "PionK0";
    break;
  case  AliHFEV0cuts::kRecoPionL:
    specname = "PionL";
    break;
  case  AliHFEV0cuts::kRecoProton:
    specname = "Proton";
    break;
  default:
    AliDebug(2, Form("Species %d not investigated", species));
    return;
  };
  AliVParticle *recTrack = NULL;
  //  mcTrack =fMC->GetTrack(TMath::Abs(label));
  //   if(!mcTrack){
  //     AliDebug(4, "MC track not available");
  //     continue; // we don't know
  //   }

  TIter trackIter(particles);

  Double_t quantities[2];
  Double_t pidProbs[5];

  while((recTrack = dynamic_cast<AliVParticle *>(trackIter()))){
    if(!TString(recTrack->IsA()->GetName()).CompareTo("AliESDtrack")){
      // case ESD
      AliESDtrack *esdTrack = dynamic_cast<AliESDtrack *>(recTrack);
      if(!esdTrack) continue;
      status = esdTrack->GetStatus();

      //TPC momentum and likelihoods
      Double_t pTPC = 0.;
      pTPC = esdTrack->GetInnerParam() ? esdTrack->GetInnerParam()->P() : esdTrack->P();
      Bool_t mcFound = kFALSE;
      if(fMC){
	Int_t label = esdTrack->GetLabel();
	AliMCParticle *mcTrack = dynamic_cast<AliMCParticle *>(fMC->GetTrack(label));
	Int_t pdg = GetPDG(species);
	Int_t trackPdg = 0; 
	if(mcTrack){
	  trackPdg = TMath::Abs(mcTrack->Particle()->GetPdgCode());
	}
	if(pdg == trackPdg) mcFound = kTRUE;	
      }
      quantities[0] = pTPC;
      Bool_t detFlagSet = kFALSE;
      for(Int_t idet = 0; idet < 4; idet++){
        TString histname, histnameMC;
	histname = "h" + detname[idet] + "_El_like_" + specname;
	histnameMC = "h" + detname[idet] + "_El_like_MC_" + specname;

        switch(idet){
          case kITS:  esdTrack->GetITSpid(pidProbs);
                      detFlagSet = status & AliESDtrack::kITSpid;
                      break;
          case kTPC:  esdTrack->GetTPCpid(pidProbs);
                      detFlagSet = status & AliESDtrack::kTPCpid;
                      break;
          case kTRD:  esdTrack->GetTRDpid(pidProbs);
                      detFlagSet = status & AliESDtrack::kTRDpid;
                      break;
          case kTOF:  esdTrack->GetTOFpid(pidProbs);  
                      detFlagSet = status & AliESDtrack::kTOFpid;
                      break;
        };
        quantities[1] = pidProbs[AliPID::kElectron];
	// in case of TRD require 6 PID tracklets
	if(kTRD == idet && esdTrack->GetTRDntrackletsPID() != 6) continue;
        if(detFlagSet){
	  fOutput->Fill(histname, quantities[0], quantities[1]);
	  if(mcFound)
	    fOutput->Fill(histnameMC, quantities[0], quantities[1]);
        }
      }
    }//.. ESD
    else{
      //AOD
    }//.. aod 
  }//.. while tracks 
}
//__________________________________________
void AliHFEpidQA::FillPIDresponse(const TObjArray * const particles, Int_t species){
  //
  // Fill the PID response of different detectors to V0 daughter particles
  //
  TString hname, hname2, hname3;
   
  const TString typeName[5] = {"Electron", "PionK0", "PionL", "Kaon", "Proton"};
  const Int_t typePID[5] = {0, 2, 2, 3, 4};

  // PID THnSparse
  // axes:
  // 0) species,  1) momentum, 2) DCA xy, 3) DCA z 
  // 4) ITS signal 
  // 5) TPC Ncls 6) TPC signal 7) TPC nSigma, 
  // 8) TRD Ntrk, 9) TRD Ncls, 10) TRD dEdx, 

  //Double_t data[12];
 

  Int_t run = fEvent->GetRunNumber();
    
  AliVParticle *recTrack = NULL;
  TIter trackIter(particles); 
  while((recTrack = dynamic_cast<AliVParticle *>(trackIter()))){
    //for(Int_t i=0; i<12; ++i) data[i] = -99.;
    // ESD
    if(!TString(recTrack->IsA()->GetName()).CompareTo("AliESDtrack")){
      // case ESD
      AliESDtrack *esdTrack = dynamic_cast<AliESDtrack *>(recTrack);
      if(!esdTrack) continue;
      
      // for the PID THnSparse
      //data[0] = species;
      //data[1] = esdTrack->P();
      Float_t impactR = -1.;
      Float_t impactZ = -1.;
      esdTrack->GetImpactParameters(impactR, impactZ);
      //data[2] = impactR;
      //data[3] = impactZ;
      //data[11] = 0; // initialize the TOF pid cut on elecgrons to false
      // use ONLY tracks with PID flag TRUE
      ULong_t status = 0;
      status = esdTrack->GetStatus();

      //
      // DEBUG 
      //
      
      fOutput->Fill("hITS_kFlags", 5., species);
      if(status & AliESDtrack::kITSin)    fOutput->Fill("hITS_kFlags", 1., species);
      if(status & AliESDtrack::kITSout)   fOutput->Fill("hITS_kFlags", 2., species);
      if(status & AliESDtrack::kITSrefit) fOutput->Fill("hITS_kFlags", 3., species);
      if(status & AliESDtrack::kITSpid)   fOutput->Fill("hITS_kFlags", 4., species);
      
      fOutput->Fill("hTPC_kFlags", 5., species);
      if(status & AliESDtrack::kTPCin)    fOutput->Fill("hTPC_kFlags", 1., species);
      if(status & AliESDtrack::kTPCout)   fOutput->Fill("hTPC_kFlags", 2., species);
      if(status & AliESDtrack::kTPCrefit) fOutput->Fill("hTPC_kFlags", 3., species);
      if(status & AliESDtrack::kTPCpid)   fOutput->Fill("hTPC_kFlags", 4., species);

      fOutput->Fill("hTRD_kFlags", 5., species);
      if(status & AliESDtrack::kTRDin)    fOutput->Fill("hTRD_kFlags", 1., species);
      if(status & AliESDtrack::kTRDout)   fOutput->Fill("hTRD_kFlags", 2., species);
      if(status & AliESDtrack::kTRDrefit) fOutput->Fill("hTRD_kFlags", 3., species);
      if(status & AliESDtrack::kTRDpid)   fOutput->Fill("hTRD_kFlags", 4., species);

      fOutput->Fill("hTOF_kFlags", 5., species);
      if(status & AliESDtrack::kTOFin)    fOutput->Fill("hTOF_kFlags", 1., species);
      if(status & AliESDtrack::kTOFout)   fOutput->Fill("hTOF_kFlags", 2., species);
      if(status & AliESDtrack::kTOFrefit) fOutput->Fill("hTOF_kFlags", 3., species);
      if(status & AliESDtrack::kTOFpid)   fOutput->Fill("hTOF_kFlags", 4., species);

      
      //
      // ITS - 
      //
      if(status & AliESDtrack::kITSpid){
	Double_t p = esdTrack->P();

	// ITS signal
	//Double_t itsSignal = esdTrack->GetITSsignal();
	
	// ITS dEdx
	Double_t dEdxSamples[4];
	esdTrack->GetITSdEdxSamples(dEdxSamples);
	Int_t nSamples = 0;
	Double_t dEdxSum = 0.;
	hname = "hITS_dEdx_" + typeName[species];
	for(Int_t i=0; i<4; ++i){
	  if(dEdxSamples[i] > 0){
	    nSamples++;
	    fOutput->Fill(hname, i+1, p, dEdxSamples[i]);
	    dEdxSum += dEdxSamples[i];
	  }
	}
	if(4 == nSamples)fOutput->Fill(hname, 0, p, dEdxSum);
	fOutput->Fill("hITS_dEdx_nSamples", nSamples);

	Double_t signal = esdTrack->GetITSsignal();
	hname = "hITS_Signal_" + typeName[species];
	fOutput->Fill(hname, p, signal);
	//data[4] = signal;
	
	// ITS number of signas
	Double_t nsigma = fESDpid->NumberOfSigmasITS(esdTrack,(AliPID::EParticleType)typePID[species]);
	hname = "hITS_nSigma_" + typeName[species];
	fOutput->Fill(hname, p, nsigma);
	// ITS PID response
	Double_t itsPID[5] = {-1, -1, -1, -1, -1};
	esdTrack->GetITSpid(itsPID);
	Int_t ix = GetMaxPID(itsPID);
	hname = "hITS_PID_p_" + typeName[species];
	fOutput->Fill(hname, p, ix);
      }//.. kITSpid
      
      //
      // TPC
      //
      if(status & AliESDtrack::kTPCpid){
	// Make TPC clusters Plot
	//data[5] = esdTrack->GetTPCNcls();
	FillTPCinfo(esdTrack, species);

	Double_t p = esdTrack->GetInnerParam() ? esdTrack->GetInnerParam()->P() : esdTrack->P();
	// TPC dEdx
	Double_t dEdx = esdTrack->GetTPCsignal();
	hname = "hTPC_dEdx_" + typeName[species];
	fOutput->Fill(hname, p, dEdx);
	//data[6] = dEdx;
	
	//TPC number of sigmas
	Double_t nsigma = fESDpid->NumberOfSigmasTPC(esdTrack,(AliPID::EParticleType)typePID[species]);
	hname = "hTPC_nSigma_" + typeName[species];
	fOutput->Fill(hname, p, nsigma);
	//data[7] = nsigma;

	// TPC PID response
	hname = "hTPC_PID_p_" + typeName[species];
	Double_t tpcPID[5] = {-1, -1, -1, -1, -1};
	esdTrack->GetTPCpid(tpcPID);
	Int_t ix = GetMaxPID(tpcPID);
	fOutput->Fill(hname, p, ix);	
      }//.. kTPCpid
      
      //
      // TRD
      //
      
      if(status & AliESDtrack::kTRDpid){
	Double_t p = esdTrack->GetOuterParam() ? esdTrack->GetOuterParam()->P() : esdTrack->P();
	
	// TRD number of tracklets
	Int_t ntrk = esdTrack->GetTRDntrackletsPID();
	hname = "hTRD_trk_" + typeName[species];
	fOutput->Fill(hname, p, ntrk);
	//data[8] = ntrk;

	// TRD PID response
	hname = "hTRD_PID_p_" + typeName[species];
	Double_t trdPID[5] = {-1., -1., -1., -1., -1};
	esdTrack->GetTRDpid(trdPID);
	Int_t ix = GetMaxPID(trdPID);
	fOutput->Fill(hname, p, ix);
	// TRD n clusters
	Int_t ncls = esdTrack->GetTRDncls();
	hname = "hTRD_cls_" + typeName[species];
	fOutput->Fill(hname, p, ncls);
	//data[9] = ncls;

	// TRD - per tracklet - dEdx per, likelihood
	hname = "hTRD_Nslc_" + typeName[species];
	hname2 = "hTRD_slc_" + typeName[species];
	hname3 = "hTRD_dEdx_" + typeName[species];
	Int_t nSlices = esdTrack->GetNumberOfTRDslices();
	Double_t sumTot = 0.;
	Int_t not0Tot = 0;
	for(Int_t l=0; l< 6; ++l){
	  Double_t trkData[7] = {-1.,-1, -1, -1, -1, -1, -1};
	  trkData[0] = species;
	  trkData[1] = p;
	  trkData[2] = l;
	  trkData[3] = ntrk;
	  trkData[5] = ncls;	    
	  Double_t sum = 0.;
	  Int_t not0 = 0;
	  for(Int_t s=0; s<nSlices; ++s){
	    Double_t slice = esdTrack->GetTRDslice(l, s);
	    sum += slice;
	    if(slice > 0){
	      not0 += 1;
	      fOutput->Fill(hname2, p, s);
	    }
	  }//..slices

	  trkData[4] = not0;
	  fOutput->Fill(hname, p, not0);
	  fOutput->Fill(hname3, p, sum);
	  if(sum > 0){
	    sumTot += sum;
	    not0Tot += 1;
	  }
	  // lkelihoods per layer
	  if(not0Tot > 0 && fNNref){
	    Double_t likelihoods[5] = {-1., -1., -1., -1., -1};
	    TRDlikeTracklet(l, esdTrack, likelihoods);
	    trkData[6] = likelihoods[0];
	    //printf(" -D: species: %i, P; %f : %f, s: %f\n", species, p, likelihoods[0], s);
	  }
	  if(not0Tot) fOutput->Fill("hTRDtracklets", trkData);
	}//..layers
	// average dEx per number of tracklets
	//if(0 < not0Tot)
	//  data[10] = sumTot / not0Tot;
      }//.. kTRDpid
      
      
      //
      // TOF
      //
      if(status & AliESDtrack::kTOFpid){
	Double_t p = esdTrack->GetOuterParam() ? esdTrack->GetOuterParam()->P() : esdTrack->P();
	Double_t t0 = fESDpid->GetTOFResponse().GetStartTime(esdTrack->P());

	//TOF beta
	hname = "hTOF_beta_" + typeName[species];
	Float_t beta = TOFbeta(esdTrack);	
	fOutput->Fill(hname, p, beta);
	fOutput->Fill("hTOF_beta_all", p, beta);
	// TOF nSigma
	Double_t nsigma = fESDpid->NumberOfSigmasTOF(esdTrack,(AliPID::EParticleType)typePID[species]);
	hname = "hTOF_nSigma_" + typeName[species];
	fOutput->Fill(hname, p, nsigma);
	//if(beta > 0.97 && beta < 1.03){
	//  data[11] = 1;
	//}
	
	// TOF PID response
	hname = "hTOF_PID_p_" + typeName[species];
	Double_t tofPID[5] = {-1., -1., -1., -1., -1};
	esdTrack->GetTOFpid(tofPID);
	Int_t ix = GetMaxPID(tofPID);
	fOutput->Fill(hname, p, ix);
	
	// time of flight QA
	// - distribution of (time - t0 - pion_time)
	Double_t times[5];
	esdTrack->GetIntegratedTimes(times);
	Double_t tItrackL = esdTrack->GetIntegratedLength();
	Double_t tTOFsignal = esdTrack->GetTOFsignal();
	Double_t dT = tTOFsignal - t0 - times[2];
	fOutput->Fill("hTOF_qa_TmT0mT", run*1.0, dT);
	fOutput->Fill("hTOF_qa_length", run*1.0, tItrackL);

	
      }//.. kTOFpid
      // temporary - the PIDsparse needs rebuilding
      //fOutput->Fill("PIDsparse", data);
    }
    // AOD - comming soon
    else{
      continue;
    }
  }// .. tracks in TObjArray
  

}
//__________________________________________
void AliHFEpidQA:: CheckEvent(){
  //
  // check some event variables
  //

  // check the T0 as a function of run number (less than one bin per number
  Double_t t0 = fESDpid->GetTOFResponse().GetTimeZero();
  Int_t run = fEvent->GetRunNumber();
  Double_t data[2] = {run*1.0, t0*1000.};
  fOutput->Fill("hEvent_T0", data);
  

}
//__________________________________________
TList *AliHFEpidQA::GetOutput(){
  //
  // Getter for Output histograms
  //
  return fOutput->GetList();
}

//__________________________________________
TList *AliHFEpidQA::GetV0pidQA(){
  //
  // Getter for V0 PID QA histograms
  //
  return fV0pid->GetListOfQAhistograms();
}

//__________________________________________
TList *AliHFEpidQA::GetV0pidMC(){
  //
  // Getter for V0 PID QA histograms
  //
  if(fV0pidMC)
    return fV0pidMC->GetListOfQAhistograms();
  return NULL;
}

//__________________________________________
void AliHFEpidQA::RecalculateTRDpid(AliESDtrack * /*track*/, Double_t * /*pidProbs*/) const{
  //  fTRDpidResponse->MakePID(track);
  //  track->GetTRDpid(pidProbs);
}

//__________________________________________
void AliHFEpidQA::RecalculateTRDpid(AliAODTrack * /*track*/, Double_t * /*pidProbs*/) const{
  //  fTRDpidResponse->MakePID(track, pidProbs);
}
//___________________________________________________________________
Float_t AliHFEpidQA::TOFbeta(const AliESDtrack * const track) const {
  // computes the TOF beta
  Double_t l = track->GetIntegratedLength();  // cm
  Double_t t = track->GetTOFsignal();
  Double_t t0 = fESDpid->GetTOFResponse().GetTimeZero(); // ps

  //printf("-D: l: %f, t: %f, t0: %f\n", l, t, t0);

  if(l < 360. || l > 800.) return 0.;
  if(t <= 0.) return 0.;
  if(t0 >999990.0) return 0.;
  

  t -= t0; // subtract the T0

  l *= 0.01;  // cm ->m
  t *= 1e-12; //ps -> s

  
  Double_t v = l / t;
  Float_t beta = v / TMath::C();

  return beta;
}
//____________________________________________
Int_t AliHFEpidQA::GetMaxPID(const Double_t *pidProbs) const {
  //
  // return the index of maximal PID probability
  //
  Int_t ix = -1;
  Double_t tmp = 0.2;
  for(Int_t i=0; i<5; ++i){
    if(pidProbs[i] > tmp){
      ix = i;
      tmp = pidProbs[i];
    }
  }
  return ix;
}
//_____________________________________________
Int_t AliHFEpidQA::GetPDG(Int_t species){
  //
  // return the PDG particle code
  //

  Int_t pdg = 0;

  switch(species){
  case  AliHFEV0cuts::kRecoElectron:
    pdg = TMath::Abs(kElectron);
    break;
  case  AliHFEV0cuts::kRecoPionK0:
    pdg = TMath::Abs(kPiPlus);
    break;
  case  AliHFEV0cuts::kRecoPionL:
    pdg = TMath::Abs(kPiPlus);
    break;
  case  AliHFEV0cuts::kRecoProton:
    pdg = TMath::Abs(kProton);
    break;
  default:  // non investigated species
    AliDebug(3, "Species not recognised");
    return 0;
  }  

  return pdg;

}

//_____________________________________________
TObjArray * AliHFEpidQA::MakeTrackList(const TObjArray *tracks) const {
  //
  // convert list of AliHFEV0Info into a list of AliVParticle
  //
  TObjArray *output = new TObjArray;
  TIter trackInfos(tracks);
  AliHFEV0info *trackInfo = NULL;
  while((trackInfo = dynamic_cast<AliHFEV0info *>(trackInfos())))
    output->Add(trackInfo->GetTrack());

  return output;
}

//_____________________________________________
TObjArray * AliHFEpidQA::MakeCleanListElectrons(const TObjArray *electrons) const {
  //
  // Cleanup electron sample using TPC PID
  // PID requirement will allways be implemented to the pair
  // Strategy
  //
  TObjArray *tracks = new TObjArray;
  TIter candidates(electrons);
  AliESDEvent *esd; 
  //AliAODEvent *aod;
  AliHFEV0info *hfetrack;
  // const Double_t kSigmaTight = 1.;
  // const Double_t kSigmaLoose = 4.;
  const Double_t kSigmaTight = 2.;
  const Double_t kSigmaLoose = 2.;
  const Double_t shift = 0.57;
  if((esd = dynamic_cast<AliESDEvent *>(fEvent))){
    AliESDtrack *track = NULL, *partnerTrack = NULL;
    while((hfetrack = dynamic_cast<AliHFEV0info *>(candidates()))){
      track = dynamic_cast<AliESDtrack *>(hfetrack->GetTrack());
      if(!track) continue;
      partnerTrack = esd->GetTrack(hfetrack->GetPartnerID());
      Double_t nSigmaTrack = TMath::Abs(fESDpid->NumberOfSigmasTPC(track, AliPID::kElectron) - shift);
      Double_t nSigmaPartner = TMath::Abs(fESDpid->NumberOfSigmasTPC(partnerTrack, AliPID::kElectron) - shift);
      if((nSigmaTrack < kSigmaTight && nSigmaPartner < kSigmaLoose) || (nSigmaTrack < kSigmaLoose && nSigmaPartner < kSigmaTight))
        tracks->Add(track);
    }
  } 
  /*else {
    aod = dynamic_cast<AliAODEvent *>(fEvent);
    if(!aod) return NULL;
    //AliAODTrack *track = NULL, *partnerTrack = NULL;
    while((hfetrack = dynamic_cast<AliHFEV0info *>(candidates()))){
      if(!hfetrack) continue;
      //track = dynamic_cast<AliAODTrack *>(hfetrack->GetTrack());
      //partnerTrack = aod->GetTrack(hfetrack->GetPartnerID());
      // will be coming soon
    }
  }*/
  return tracks;
}
//___________________________________________________________
void AliHFEpidQA::CheckTenderV0pid(const TObjArray * const particles, Int_t species){

  //
  // retrieve the status bits from the TObject used to flag
  // the V0 daughter tracks and return the found PID
  // 0 - electron, 1 - pion, 2 - proton
  // 

  const Int_t id[5] = {0, 1, 1, -1, 2}; // convert AliHFEpid to simple 0, 1, 2

  AliVParticle *recTrack = NULL;
  TIter trackIter(particles); 
  while((recTrack = dynamic_cast<AliVParticle *>(trackIter()))){
    if(!TString(recTrack->IsA()->GetName()).CompareTo("AliESDtrack")){
      // case ESD
      AliESDtrack *track = dynamic_cast<AliESDtrack *>(recTrack);
      if(!track) continue;
      Int_t tPID = GetTenderV0pid(track);
      fOutput->Fill("h_tender_check_01", id[species]*1.0, tPID);
      fOutput->Fill("h_tender_check_02", id[species], id[species]*1.0, tPID);
      
    } //.. case ESD      
    
  }//.. iterate of tracks  
}
//___________________________________________________________
Int_t AliHFEpidQA::GetTenderV0pid(AliESDtrack * const track){
  //
  // retrieve the PID nformation stored in the status flags by the train tender
  //
  

  Int_t pid = -1;
  if(!track){
    return pid;
  }
  
  Int_t nTimes = 0;
  
  if(track->TestBit(2<<14)){
    pid = 0;
    nTimes++;
  }
  if(track->TestBit(2<<15)){
    pid = 1;
    nTimes++;
  }
  if(track->TestBit(2<<16)){
    pid = 2;
    nTimes++;
  }
  
  if(nTimes > 1){
    AliWarning("V0 track labeled multiple times by the V0 tender");
    pid = -1;
  }
    
  return pid;

}
//___________________________________________________________
Double_t AliHFEpidQA::TRDlikeTracklet(Int_t layer, AliESDtrack * const track, Double_t *likelihood){
  //
  // compute the TRD electron likelihoods for 1 tracklet
  // based on teh AliTRDpidRecalculator in train/until/tender
  // returns sum of the likelihoods (which should be 1)
  //

  const Double_t cScaleGain = 1./ 16000.;
  const Float_t pBins[11] ={0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0}; // momentum bins

  if(!track) return kFALSE;
  Float_t p = track->GetTRDmomentum(layer); // momentum for a tracklet in the ESDtrack
  if(p < 0) return kFALSE;

  Int_t mombin = TRDmomBin(p);  // momentum bin
  if(mombin < 0) return -1.;
  Float_t dEdxTRDsum = 0;              // dEdxTRDsum for checking if tracklet is available
  Float_t dEdxTRD[8];           // dEdx for a tracklet in the ESD slices
  Double_t ddEdxTRD[8];         // dEdx as Double_t for TMultiLayerPerceptron::Evaluate()
  
  Double_t prob[AliPID::kSPECIES];  // probabilities for all species in all layers

  for(Int_t is = 0; is < AliPID::kSPECIES; is++){
    likelihood[is] = 0.2;               // init probabilities
    prob[is] = 0.2; 
  }

  Double_t sum = 0.;

  for(Int_t islice = 0; islice<8; islice++){
    dEdxTRD[islice]=0.;                                    // init dE/dx
    ddEdxTRD[islice]=0.;                                   // init dE/dx
    dEdxTRD[islice]=track->GetTRDslice(layer,islice);       // get the dE/dx values
    dEdxTRDsum += dEdxTRD[islice];
    ddEdxTRD[islice]=(Double_t)dEdxTRD[islice]*cScaleGain;            // rescale dE/dx
    
  }
  for(Int_t is = 0; is < AliPID::kSPECIES; is++){
    Double_t probloc1, probloc2; 
    if(mombin == 0 && mombin < pBins[0]){          //  calculate PID for p > 0.6 GeV/c
      prob[is] = fNet[mombin]->Evaluate(is, ddEdxTRD);
    }
    else if(mombin == 10 && mombin >= pBins[10]){   //  calculate PID for p >= 10.0 GeV/c
      prob[is] = fNet[mombin]->Evaluate(is, ddEdxTRD);
    }
    else{                                                    //  standard calculation
      Int_t mombin1 = 0, mombin2 = 0;             // lower and upper momentum bin      
      if(p < pBins[mombin]) {mombin1 = mombin -1; mombin2 = mombin;}
      if(p >= pBins[mombin]) {mombin1 = mombin; mombin2 = mombin+1;}
      probloc1 = fNet[mombin1]->Evaluate(is, ddEdxTRD);
      probloc2 = fNet[mombin2]->Evaluate(is, ddEdxTRD);
      // weighting of the probability with regard to the track momentum
      prob[is] = probloc1 + (probloc2-probloc1)*(p-pBins[mombin1])/(pBins[mombin2]-pBins[mombin1]);
    }
    likelihood[is] = prob[is];
    sum += likelihood[is];
  }     
  
  return sum;
}
//__________________________________________________________________________
Int_t AliHFEpidQA::TRDmomBin(Double_t p) const {
  //
  // compute the momentum bin position
  // 

  const Float_t pBins[11] ={0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0}; // momentum bins
  
  Int_t pBin1 = -1;                                    // check bin1
  Int_t pBin2 = -1;                                    // check bin2

  if(p < 0) return -1;                                 // return -1 if momentum < 0
  if(p < pBins[0]) return 0;                                      // smallest momentum bin
  if(p >= pBins[10]) return 10; // largest momentum bin


  // calculate momentum bin for non extremal momenta
  for(Int_t iMomBin = 1; iMomBin < 11; iMomBin++){
    if(p < pBins[iMomBin]){
      pBin1 = iMomBin - 1;
      pBin2 = iMomBin;
    }
    else
      continue;

    if(p - pBins[pBin1] >= pBins[pBin2] - p){
      return pBin2;
    }
    else{
      return pBin1;
    }
  }

  return -1;


}
//__________________________________________________________________________


