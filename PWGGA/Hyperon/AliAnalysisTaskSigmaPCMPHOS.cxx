// @(#)root/base:$Id$
// Authors: Alexander Borissov, Sergei Solokhin    03/01/23

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

#include "AliRun.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSigmaPCMPHOS.h"
#include "AliTriggerAnalysis.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "AliAODMCParticle.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliPIDResponse.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliPID.h"
#include "AliAODInputHandler.h"
#include "TMath.h"
#include "TPDGCode.h"
#include "AliAODTrack.h"
#include "AliMCEventHandler.h"

// Analysis task to fill histograms with PHOS AOD clusters and cells
// Authors: Yuri Kharlov
// Date   : 28.05.2009

ClassImp(AliAnalysisTaskSigmaPCMPHOS)

    ////////////////////////////////////////////////////////////////////////////////
    //Default class constructor
    AliAnalysisTaskSigmaPCMPHOS::AliAnalysisTaskSigmaPCMPHOS() : AliAnalysisTaskSE(),
                                                                 fAODEvent(0), fGlobalEventID(0), fMCEvent(0x0),
                                                                 fPIDResponse(0),
                                                                 fOutputList(0),
                                                                 fStack(0x0),
                                                                 fPHOSGeo(0x0),
                                                                 fTriggerAnalysis(new AliTriggerAnalysis),
                                                                 fElectronMass(0),
                                                                 fGamma(0x0),
                                                                 fPCM(0x0),
                                                                 fPi0(0x0),
                                                                 fPi0Merged(0x0),
                                                                 fTracksPip(0x0),
                                                                 fTracksPim(0x0),
                                                                 fTracksPp(0x0),
                                                                 fTracksPm(0x0),
                                                                 fLambda(0x0),
                                                                 fMixGamma(0x0),
                                                                 fMixTracksPp(0x0),
                                                                 fMixTracksPm(0x0),
                                                                 fMixLambda(0x0),
                                                                 fAODMCTrackArray(0x0), fOnFlyVector(0x0), fFinderVector(0x0), fV0ParticleIDArray(0x0), fConvPhotonArray(0x0)
{
}

////////////////////////////////////////////////////////////////////////////////
// Class constructor for I/O operations
AliAnalysisTaskSigmaPCMPHOS::AliAnalysisTaskSigmaPCMPHOS(const char *name) : AliAnalysisTaskSE(name),
                                                                             fAODEvent(0), fGlobalEventID(0), fMCEvent(0x0),
                                                                             fPIDResponse(0),
                                                                             fOutputList(0),
                                                                             fStack(0x0),
                                                                             fPHOSGeo(0x0),
                                                                             fTriggerAnalysis(new AliTriggerAnalysis),
                                                                             fElectronMass(0),
                                                                             fGamma(0x0),
                                                                             fPCM(0x0),
                                                                             fPi0(0x0),
                                                                             fPi0Merged(0x0),
                                                                             fTracksPip(0x0),
                                                                             fTracksPim(0x0),
                                                                             fTracksPp(0x0),
                                                                             fTracksPm(0x0),
                                                                             fLambda(0x0),
                                                                             fMixGamma(0x0),
                                                                             fMixTracksPp(0x0),
                                                                             fMixTracksPm(0x0),
                                                                             fMixLambda(0x0),
                                                                             fAODMCTrackArray(0x0), fOnFlyVector(0x0), fFinderVector(0x0), fV0ParticleIDArray(0x0), fConvPhotonArray(0x0)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

////////////////////////////////////////////////////////////////////////////////
// Default class destructor
AliAnalysisTaskSigmaPCMPHOS::~AliAnalysisTaskSigmaPCMPHOS()
{
  if (fOutputList)
    delete fOutputList;
}

////////////////////////////////////////////////////////////////////////////////
// User-defined output objects, called once
void AliAnalysisTaskSigmaPCMPHOS::UserCreateOutputObjects()
{
  // Create AOD histograms

  //Save Particle Masses and other constants for later use
  fElectronMass = TDatabasePDG::Instance()->GetParticle(11)->Mass();

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fOutputList->Add(new TH1F("hClusterEnergy", "Cluster energy", 500, 0., 5.));
  fOutputList->Add(new TH2F("hClusterTOFvsE", "Cluster time vs energy", 100, 0., 10.e-7, 40, 0., 20.));

  fOutputList->Add(new TH2F("hLambdaMass", "gg mass", 50, 1.105, 1.125, 20, 0., 10.));
  fOutputList->Add(new TH2F("hLambdaBarMass", "gg mass", 50, 1.105, 1.125, 20, 0., 10.));

  //Gamma-gamma
  Int_t nMgg = 100;
  Double_t mggMax = 2.;
  Double_t ptggMax = 10.;
  Int_t nPtgg = 40;
  Double_t mHe3Max = 3.5;

  fOutputList->Add(new TH1F("hZvertex", "Z vertex", 200, -50., +50.));
  fOutputList->Add(new TH1F("hNvertexTracks", "N of primary tracks from the primary vertex", 150, 0., 150.));
  fOutputList->Add(new TH2F("hmHe3pGam", "hmHe3pGam", nMgg, 3.0, mHe3Max, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmHe3mGam", "hmHe3mGam", nMgg, 3.0, mHe3Max, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hGamGam", "hGamGam", nMgg, 0, 1., nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hGamGv0", "hGamGv0", nMgg, 0, 1., nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hGv0Gv0", "hGv0Gv0", nMgg, 0, 1., nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmGamGam", "hmGamGam", nMgg, 0, 1., nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hLamGam", "hLamGam", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hLamGamAlpha", "hLamGamAlpha", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmLamGam", "hmLamGam", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmixLamGam", "hmixLamGam", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmLamGamAlpha", "hmLamGamAlpha", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmixLamGamAlpha", "hmixLamGamAlpha", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmcLamGam", "hmcLamGam", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmcLamGamAlpha", "hmcLamGamAlpha", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmcLamGamAlpha1", "hmcLamGamAlpha1", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmcLamGamAlpha2", "hmcLamGamAlpha2", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmc4piLamGam", "hmc4piLamGam", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hLamGv0", "hLamGv0", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hLamGv0Alpha", "hLamGv0Alpha", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmLamGv0", "hmLamGv0", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));
  fOutputList->Add(new TH2F("hmLamGv0Alpha", "hmLamGv0Alpha", nMgg, 1.117, 1.220, nPtgg, 0., ptggMax));

  fOutputList->Add(new TH1F("hSelEvents", "Selected events", 12, 0.5, 12.5));
  fMixGamma = new TList();
  fMixGamma->SetOwner(kTRUE);
  fMixLambda = new TList();
  fMixLambda->SetOwner(kTRUE);
  PostData(1, fOutputList);
}

////////////////////////////////////////////////////////////////////////////////
// User-defined event analysis, called for each event
void AliAnalysisTaskSigmaPCMPHOS::UserExec(Option_t *option)
{
  // Main loop, called for each event  Analyze AOD

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr || !mgr->GetInputEventHandler())
    return;

  fAODEvent = dynamic_cast<AliAODEvent *>(InputEvent());
  if (!fAODEvent)
    return;

  AliAODHeader *aodHeader = dynamic_cast<AliAODHeader *>(fAODEvent->GetHeader());
  if (!aodHeader)
    return;

  AliVHeader *fAODEventHeader = fAODEvent->GetHeader(); //Get the Header from the event
  if (!fAODEventHeader)
    return;

  else
    fGlobalEventID = fAODEventHeader->GetEventIdAsLong(); //Get global ID of the event

  Int_t nTracks = fAODEvent->GetNumberOfTracks();
  fPIDResponse = mgr->GetInputEventHandler()->GetPIDResponse();

  // Event selection flags
  Bool_t fVerbose = kFALSE;
  Bool_t eventVtxExist = kFALSE;
  Bool_t eventVtxZ10cm = kFALSE;
  Bool_t eventPileup = kFALSE;
  Bool_t eventV0AND = kFALSE;
  Bool_t fAODEventVtxExist = kFALSE;

  FillHistogram("hSelEvents", 1);

  // Check the PID response
  if (!fPIDResponse)
    return;

  //Number of Tracks
  if (nTracks == 0)
    return; //No point in continuing if there are no tracks

  // Check primary vertex position
  Double_t primaryVtxPos[3] = {-999, -999, -999};
  const AliAODVertex *aodVtx = fAODEvent->GetPrimaryVertex();
  if (!aodVtx)
    return;

  aodVtx->GetXYZ(primaryVtxPos);

  
  //Fill V0 arrays
  fOnFlyVector.clear(); //clear the arrays
  fFinderVector.clear();
  fV0ParticleIDArray.clear();

  Int_t nV0 = fAODEvent->GetNumberOfV0s(); //Number of V0s in the event

  for (Int_t iV0 = 0; iV0 < nV0; iV0++)
  { //Loop over V0s in the event

    AliAODv0 *aodV0 = (AliAODv0 *)fAODEvent->GetV0(iV0); //Get V0 object
    if (!aodV0)
      continue;

    // Check basic V0 properties: 2 Daughters, opposite charge, total charge = 0
    if (aodV0->GetNDaughters() != 2)
      continue;
    if (aodV0->GetNProngs() != 2)
      continue;
    if (aodV0->GetCharge() != 0)
      continue;
    if (aodV0->ChargeProng(0) == aodV0->ChargeProng(1))
      continue;

    // Get daughter tracks
    AliAODTrack *trackN = dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(0));
    AliAODTrack *trackP = dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(1));
    if (!trackN || !trackP)
      continue;
    if (trackN->GetSign() == trackP->GetSign())
      continue;

    fFinderVector.push_back(trackN->GetID());
    fFinderVector.push_back(trackP->GetID());

    if (aodV0->GetOnFlyStatus())
    {
      fOnFlyVector.push_back(trackN->GetID());
      fOnFlyVector.push_back(trackP->GetID());
    }
  } //End of V0 Loop. Finished prearing the maps ==> To understand

  FillHistogram("hSelEvents", 3);

  /************************Start Event Processing**************************************/

  fMCEvent = MCEvent(); // Get MC event (called fMCEvent) from the input file

  if (fMCEvent)
  {
    fAODMCTrackArray = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fAODMCTrackArray)
      return;
    if (fMCEvent->Stack())
    {
      printf("-----------------------");
      fStack = static_cast<AliStack *>(fMCEvent->Stack());
      if (fStack)
        ProcessMC();
    }
  }

  FillV0PhotonArray();

  FillHistogram("hSelEvents", 4);

  if (!fPHOSGeo)
    fPHOSGeo = AliPHOSGeometry::GetInstance();

  // Checks if we have a primary vertex	// Get primary vertices form AOD
  if (fAODEvent->GetPrimaryVertexTracks()->GetNContributors() > 0)
    fAODEventVtxExist = kTRUE;
  else if (fAODEvent->GetPrimaryVertexSPD()->GetNContributors() > 0)
    fAODEventVtxExist = kTRUE;
  const AliAODVertex *esdVertex5 = fAODEvent->GetPrimaryVertex();

  Double_t vtx5[3];
  vtx5[0] = esdVertex5->GetX();
  vtx5[1] = esdVertex5->GetY();
  vtx5[2] = esdVertex5->GetZ();
  TVector3 vtx(vtx5);

  FillHistogram("hNvertexTracks", esdVertex5->GetNContributors());
  FillHistogram("hZvertex", esdVertex5->GetZ());
  if (TMath::Abs(esdVertex5->GetZ()) > 10.)
    return;

  FillHistogram("hSelEvents", 5);

  if (fAODEvent->IsPileupFromSPD())
    return;

  FillHistogram("hSelEvents", 6);

  // Fill event statistics for different selection criteria

  //Vtx class z-bin
  Int_t zvtx = (Int_t)((vtx5[2] + 10.) / 2.);
  if (zvtx < 0)
    zvtx = 0;
  if (zvtx > 9)
    zvtx = 9;
  FillHistogram("hSelEvents", 7);

  //        SelectHadrons() ;
  FillHistogram("hSelEvents", 8);
  SelectLambda();
  FillHistogram("hSelEvents", 9);
  SelectGamma();
  FillHistogram("hSelEvents", 10);

  //================REALs===========
  Int_t nGamma = fGamma->GetEntriesFast();
  Int_t nPCM = 1; // fPCM->GetEntriesFast() ;
  Int_t nLambda = fLambda->GetEntriesFast();
  Int_t nPp = 0; //fTracksPp->GetEntriesFast() ;
  Int_t nPm = 0; //fTracksPm->GetEntriesFast() ;
  FillHistogram("hSelEvents", 11);

  TLorentzVector pair;

  // fill gamma-gamma  ?repeat without tagged photons
  for (Int_t i = 0; i < nGamma; i++)
  {
    AliCaloPhoton *pv1 = (AliCaloPhoton *)fGamma->At(i);
    for (Int_t j = i + 1; j < nGamma; j++)
    {
      AliCaloPhoton *pv2 = (AliCaloPhoton *)fGamma->At(j);
      pair = *pv1 + *pv2;
      Double_t m = pair.M();
      Double_t pT = pair.Pt();
      FillHistogram("hGamGam", m, pT);
    }
  }

  FillHistogram("hSelEvents", 12);

  //Fill gamma-Lambda
  for (Int_t i = 0; i < nGamma; i++)
  {
    AliCaloPhoton *pv1 = (AliCaloPhoton *)fGamma->At(i);
    for (Int_t j = 0; j < nLambda; j++)
    {
      TLorentzVector *pv2 = (TLorentzVector *)fLambda->At(j);
      pair = *pv1 + *pv2;
      Double_t alpha = (pv1->E() - pv2->E()) / (pv1->E() + pv2->E());
      Double_t m = pair.M();
      Double_t pT = pair.Pt();
      FillHistogram("hLamGam", m, pT); //    -> Fill (m,pT);
      if (alpha < -0.7)
        FillHistogram("hLamGamAlpha", m, pT);
    }
  }

  // Fill gamma_v0-Lambda, idea to separate on Lambdap - particle and Lambdam - antiparticle
  TLorentzVector pvv1, pvv2;
  const Int_t nConvPhoton = fConvPhotonArray.size();
  for (Int_t i = 0; i < nConvPhoton - 1; i++)
  {
    AliAODv0 *v0_1 = (AliAODv0 *)fAODEvent->GetV0(fConvPhotonArray.at(i));
    if (!v0_1)
      continue;
    pvv1.SetXYZM(v0_1->Px(), v0_1->Py(), v0_1->Pz(), 0);
    for (Int_t j = 0; j < nLambda; j++)
    {
      TLorentzVector *pv2 = (TLorentzVector *)fLambda->At(j);
      pair = pvv1 + *pv2;
      Double_t alpha = (pvv1.E() - pv2->E()) / (pvv1.E() + pv2->E());
      Double_t pT = pair.Pt();
      Double_t m = pair.M();
      FillHistogram("hLamGv0", m, pT);
      if (alpha < -0.7)
        FillHistogram("hLamGv0Alpha", m, pT);
    }
    for (Int_t j = 0; j < nGamma; j++)
    {
      AliCaloPhoton *pv2 = (AliCaloPhoton *)fGamma->At(j);
      pair = pvv1 + *pv2;
      Double_t m = pair.M();
      Double_t pT = pair.Pt();
      FillHistogram("hGamGv0", m, pT); //   -> Fill (m,pT);
    }
    for (Int_t j = i + 1; j < nConvPhoton; j++)
    {

      AliAODv0 *v0_2 = (AliAODv0 *)fAODEvent->GetV0(fConvPhotonArray.at(j));
      if (!v0_2)
        continue;
      pvv2.SetXYZM(v0_2->Px(), v0_2->Py(), v0_2->Pz(), 0);
      pair = pvv1 + pvv2;
      Double_t m = pair.M();
      Double_t pT = pair.Pt();
      FillHistogram("hGv0Gv0", m, pT); //   -> Fill (m,pT);
    }
  } // end of converted photons

  //================MIXEDs - later, when REAL events will wi checked
  // * Fill gamma-gamma
  for (Int_t m = 0; m < fMixGamma->GetSize(); m++)
  {
    TClonesArray *tmp = (TClonesArray *)fMixGamma->At(m);
    for (Int_t j = 0; j < tmp->GetEntriesFast(); j++)
    {
      AliCaloPhoton *pv2 = (AliCaloPhoton *)tmp->At(j);
      for (Int_t i = 0; i < nGamma; i++)
      {
        AliCaloPhoton *pv1 = (AliCaloPhoton *)fGamma->At(i);
        pair = *pv1 + *pv2;
        Double_t m = pair.M();
        Double_t pT = pair.Pt();
        FillHistogram("hmGamGam", m, pT); //   -> Fill (m,pT);
      }
    }
  }

  for (Int_t m = 0; m < fMixGamma->GetSize(); m++)
  {
    TClonesArray *tmp = (TClonesArray *)fMixGamma->At(m);
    for (Int_t i = 0; i < tmp->GetEntriesFast(); i++)
    {
      AliCaloPhoton *pv1 = (AliCaloPhoton *)tmp->At(i);
      for (Int_t j = 0; j < nPp; j++)
      {
        AliCaloPhoton *pv2 = (AliCaloPhoton *)fTracksPp->At(j);
        pair = *pv1 + *pv2;
        Double_t m = pair.M();
        Double_t pT = pair.Pt();
        FillHistogram("hmHe3pGam", m, pT); //   -> Fill (m,pT);
      }
    }
  }

  for (Int_t m = 0; m < fMixGamma->GetSize(); m++)
  {
    TClonesArray *tmp = (TClonesArray *)fMixGamma->At(m);
    for (Int_t i = 0; i < tmp->GetEntriesFast(); i++)
    {
      AliCaloPhoton *pv1 = (AliCaloPhoton *)tmp->At(i);
      for (Int_t j = 0; j < nPm; j++)
      {
        AliCaloPhoton *pv2 = (AliCaloPhoton *)fTracksPm->At(j);
        pair = *pv1 + *pv2;
        Double_t m = pair.M();
        Double_t pT = pair.Pt();
        FillHistogram("hmHe3mGam", m, pT); //   -> Fill (m,pT);
      }
    }
  }

  // * Fill gamma-mix.Lambda
  for (Int_t m = 0; m < fMixLambda->GetSize(); m++)
  {
    TClonesArray *tmp = (TClonesArray *)fMixLambda->At(m);
    for (Int_t j = 0; j < tmp->GetEntriesFast(); j++)
    {
      TLorentzVector *pv2 = (TLorentzVector *)tmp->At(j);
      for (Int_t i = 0; i < nGamma; i++)
      {
        AliCaloPhoton *pv1 = (AliCaloPhoton *)fGamma->At(i);
        pair = *pv1 + *pv2;
        Double_t alpha = (pv1->E() - pv2->E()) / (pv1->E() + pv2->E());
        Double_t m = pair.M();
        Double_t pT = pair.Pt();
        FillHistogram("hmixLamGam", m, pT); //   -> Fill (m,pT);
        if (alpha < -0.7)
          FillHistogram("hmixLamGamAlpha", m, pT);
      }
    }
  }
  // * Fill mix.gamma-Lambda
  for (Int_t m = 0; m < fMixGamma->GetSize(); m++)
  {
    TClonesArray *tmp = (TClonesArray *)fMixGamma->At(m);
    for (Int_t i = 0; i < tmp->GetEntriesFast(); i++)
    {
      AliCaloPhoton *pv1 = (AliCaloPhoton *)tmp->At(i);
      for (Int_t j = 0; j < nLambda; j++)
      {
        TLorentzVector *pv2 = (TLorentzVector *)fLambda->At(j);
        pair = *pv1 + *pv2;
        Double_t alpha = (pv1->E() - pv2->E()) / (pv1->E() + pv2->E());
        Double_t m = pair.M();
        Double_t pT = pair.Pt();
        FillHistogram("hmLamGam", m, pT);
        if (alpha < -0.7)
          FillHistogram("hmLamGamAlpha", m, pT);
      }
    }
  }

  // * Fill gamma.v0-mix.Lambda
  for (Int_t m = 0; m < fMixLambda->GetSize(); m++)
  {
    TClonesArray *tmp = (TClonesArray *)fMixLambda->At(m);
    for (Int_t j = 0; j < tmp->GetEntriesFast(); j++)
    {
      TLorentzVector *pv2 = (TLorentzVector *)tmp->At(j);
      for (Int_t i = 0; i < nConvPhoton - 1; i++)
      {
        AliAODv0 *v0_1 = (AliAODv0 *)fAODEvent->GetV0(fConvPhotonArray.at(i));
        if (!v0_1)
          continue;
        pvv1.SetXYZM(v0_1->Px(), v0_1->Py(), v0_1->Pz(), 0);
        pair = pvv1 + *pv2;
        Double_t alpha = (pvv1.E() - pv2->E()) / (pvv1.E() + pv2->E());
        Double_t m = pair.M();
        Double_t pT = pair.Pt();
        FillHistogram("hmLamGv0", m, pT); //   -> Fill (m,pT);
        if (alpha < -0.7)
          FillHistogram("hmLamGv0Alpha", m, pT); //   -> Fill (m,pT);
      }
    }
  }

  //Now we either add current events to stack or remove ==> To check, skip for now
  //If no photons in current event - no need to add it to mixed
  const Int_t kMixEvents = 10;
  const Int_t kMixEventsHadr = 100;
  if (fGamma->GetEntriesFast() > 0)
  {
    fMixGamma->AddFirst(fGamma);
    fGamma = 0;
    if (fMixGamma->GetSize() > kMixEvents)
    { //Remove redundant events
      TClonesArray *tmp = static_cast<TClonesArray *>(fMixGamma->Last());
      fMixGamma->RemoveLast();
      delete tmp;
    }
  }

  // make mixed Lambda
  if (fLambda->GetEntriesFast() > 0)
  {
    fMixLambda->AddFirst(fLambda);
    fLambda = 0;
    if (fMixLambda->GetSize() > kMixEventsHadr)
    { //Remove redundant events
      TClonesArray *tmp = static_cast<TClonesArray *>(fMixLambda->Last());
      fMixLambda->RemoveLast();
      delete tmp;
    }
  }

  // Post output data.
  PostData(1, fOutputList);
}

////////////////////////////////////////////////////////////////////////////////
// Process simulated events if Monte-Carlo flag is true
void AliAnalysisTaskSigmaPCMPHOS::ProcessMC()
{

  //fill histograms for efficiensy etc. calculation

  const Double_t kRcut = 1.; //cut for primary particles
  Double_t vtx[3];
  vtx[0] = fAODEvent->GetPrimaryVertex()->GetX();
  vtx[1] = fAODEvent->GetPrimaryVertex()->GetY();
  vtx[2] = fAODEvent->GetPrimaryVertex()->GetZ();

  Int_t Daughter1 = 0;
  Int_t Daughter2 = 0;

  if (Int_t(fStack->GetNtrack()) < 2)
    return;

  Int_t nMCTracks = fMCEvent->GetNumberOfTracks();
  //loop over all MC tracks
  for (Int_t iMCtrack = 1; iMCtrack < nMCTracks; iMCtrack++)
  {

    AliAODMCParticle *mcPart = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(iMCtrack));
    if (!mcPart)
      continue;

    if (TMath::Abs(mcPart->Eta()) > 0.5)
      continue; //Acceptance Cut

    Int_t MCPartPDGCode = mcPart->PdgCode();

    if (!(MCPartPDGCode == 3122 || MCPartPDGCode == -3122))
      continue;

    Double_t pT = mcPart->Pt();
    Double_t m = 1.1926; // particle->M();

    FillHistogram("hmc4piLamGam", m, pT);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Fill array of photons' V0
void AliAnalysisTaskSigmaPCMPHOS::FillV0PhotonArray()
{

  // the same codes as for Sigma+ analysis
  Double_t primaryVtxPos[3] = {0, 0, 0};

  //Clear V0 Photon Array and reset counter
  fConvPhotonArray.clear();
  Int_t countPhotons = 0;

  Int_t nV0 = fAODEvent->GetNumberOfV0s(); //Number of V0s in the event
  if (nV0 == 0)
    return; //Return if there is no V0 to be processed

  Int_t non = 0;
  Int_t noff = 0;

  for (Int_t iV0 = 0; iV0 < nV0; iV0++)
  { //Loop over V0s in the event

    //Initialisation of the local Bools
    Bool_t isElectronTPC = kFALSE;
    Bool_t isPositronTPC = kFALSE;
    Bool_t isPhotonTPC = kFALSE;
    Bool_t isRealV0 = kFALSE;
    Bool_t isReallyPhoton = kFALSE;
    Bool_t isPhotonfromSigma = kFALSE;

    TVector3 vecN, vecP, vecM;                 //Momentum Vectors for V0 tracks
    TLorentzVector electron, positron, photon; //Lorentzvectors for invariant mass calculation

    // Daughter Track parameters for KF and ExtTrckPar initialization
    Double_t trackxyz[3];
    Double_t trackpxpypz[3];
    Double_t trackparams[6];
    Double_t covMatrix[21];

    AliAODv0 *aodV0 = (AliAODv0 *)fAODEvent->GetV0(iV0); //Get V0 object
    if (!aodV0)
      continue;

    // Check basic V0 properties: 2 Daughters, opposite charge, total charge = 0
    if (aodV0->GetNDaughters() != 2)
      continue;
    if (aodV0->GetNProngs() != 2)
      continue;
    if (aodV0->GetCharge() != 0)
      continue;
    if (aodV0->ChargeProng(0) == aodV0->ChargeProng(1))
      continue;

    // Get daughter tracks
    AliAODTrack *trackN = dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(0));
    AliAODTrack *trackP = dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(1));
    if (trackN->GetSign() == trackP->GetSign())
      continue;
    if (trackN->Charge() > 0)
    { //Check correct charge
      trackN = dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(1));
      trackP = dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(0));
    }

    if (!trackP || !trackN)
      continue;

    //If V0 is not On-fly, check the map if there is an equivalent On-fly V0
    if (!aodV0->GetOnFlyStatus())
    {
      noff++;

      //Check if the tracks are used by the on-the-fly finder
      Int_t nFound = fOnFlyVector.size();
      Bool_t isused = kFALSE;
      for (Int_t iID = 0; iID < nFound; iID++)
      {
        if (trackN->GetID() == fFinderVector[iID])
          isused = kTRUE;
        if (trackP->GetID() == fFinderVector[iID])
          isused = kTRUE;
      }
      if (isused)
        continue;
    }
    else
    {
      non++;
    }

    // Check track quality
    Int_t nTPCClustNeg = trackN->GetTPCNcls();
    Int_t nTPCClustPos = trackP->GetTPCNcls();

    // Daughter track PID using TPC
    Double_t nSigmaTPCelectron = fPIDResponse->NumberOfSigmasTPC(trackN, AliPID::kElectron);
    Double_t nSigmaTPCpositron = fPIDResponse->NumberOfSigmasTPC(trackP, AliPID::kElectron);
    if (TMath::Abs(nSigmaTPCelectron) < fMaxNsigDaughtTPC)
      isElectronTPC = kTRUE;
    if (TMath::Abs(nSigmaTPCpositron) < fMaxNsigDaughtTPC)
      isPositronTPC = kTRUE;
    if (isElectronTPC && isPositronTPC)
      isPhotonTPC = kTRUE;

    // Get topological values
    Double_t dcaV0Daughters = TMath::Abs(aodV0->DcaV0Daughters());
    Double_t dcaPosToPrimVtx = TMath::Abs(aodV0->DcaPosToPrimVertex());
    Double_t dcaNegToPrimVtx = TMath::Abs(aodV0->DcaNegToPrimVertex());
    Double_t cosPointAngle = aodV0->CosPointingAngle(primaryVtxPos);
    Double_t vtxPosV0[3];
    vtxPosV0[0] = aodV0->DecayVertexV0X();
    vtxPosV0[1] = aodV0->DecayVertexV0Y();
    vtxPosV0[2] = aodV0->DecayVertexV0Z();
    Double_t Vtxradius = TMath::Sqrt(vtxPosV0[0] * vtxPosV0[0] + vtxPosV0[1] * vtxPosV0[1]);

    //Calculating DCA of Photon to PV
    TVector3 CV(aodV0->DecayVertexV0X(), aodV0->DecayVertexV0Y(), aodV0->DecayVertexV0Z()); //Conv. Vertex
    TVector3 p(aodV0->Px(), aodV0->Py(), aodV0->Pz());                                      //Momentum vectors of the photons
    Double_t DCAPV = (p.Cross(CV)).Mag() / p.Mag();                                  //DCA to PV of Photons

    //Get reconstructed cartesian momentum
    vecN.SetXYZ(aodV0->MomNegX(), aodV0->MomNegY(), aodV0->MomNegZ()); //negative daughter
    vecP.SetXYZ(aodV0->MomPosX(), aodV0->MomPosY(), aodV0->MomPosZ()); //positive daughter
    vecM.SetXYZ(aodV0->MomV0X(), aodV0->MomV0Y(), aodV0->MomV0Z());    //mother

    //Custom Armenteros Podolanski calculation since V0 member functions are not reliable!
    Double_t pLNeg = vecN.Dot(vecM) / vecM.Mag(); //Momentum longitudinal
    Double_t pLPos = vecP.Dot(vecM) / vecM.Mag(); //to V0 momentum
    Double_t alpha = (pLPos - pLNeg) / (pLPos + pLNeg);
    Double_t qt = vecN.Perp(vecM);

    // Get kinematic values
    Double_t ptV0 = aodV0->Pt();
    Double_t pPos = trackP->P();
    Double_t pNeg = trackN->P();
    Double_t thetaPos = trackP->Theta();
    Double_t thetaNeg = trackN->Theta();
    Double_t totangle = aodV0->OpenAngleV0();

    //*******************AOD V0 MC treatment************************//
    // * AOD MC treatment
    if (fMCEvent)
    {
      AliAODMCParticle *V0Part = nullptr;
      AliAODMCParticle *NPart = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(TMath::Abs(trackN->GetLabel())));
      AliAODMCParticle *PPart = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(TMath::Abs(trackP->GetLabel())));
      AliAODMCParticle *V2Mother = nullptr;
      AliAODMCParticle *V1Mother = nullptr;
      AliAODMCParticle *V1Part = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(NPart->GetMother()));
      if (V1Part->GetMother() != -1)
      {
        V1Mother = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(TMath::Abs(V1Part->GetMother())));
        if (TMath::Abs(V1Mother->GetPdgCode()) == 3212)
        {
          Double_t m = V1Mother->M();
          Double_t pt = V1Mother->Pt();
          FillHistogram("hmcLamGamAlpha1", m, pt);
        }
      }

      AliAODMCParticle *V2Part = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(PPart->GetMother()));
      if (V2Part->GetMother() != -1)
      {
        V2Mother = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(TMath::Abs(V2Part->GetMother())));
        if (TMath::Abs(V2Mother->GetPdgCode()) == 3212)
        {
          Double_t m = V2Mother->M();
          Double_t pt = V2Mother->Pt();
          FillHistogram("hmcLamGamAlpha2", m, pt);
        }
      }

      if ((V1Part->GetMother() != -1) &&
          (V2Part->GetMother() != -1) &&
          (V1Mother->Pt() == V2Mother->Pt()))
      {
        if (TMath::Abs(V1Mother->GetPdgCode()) == 3212 &&
            TMath::Abs(V2Mother->GetPdgCode()) == 3212 &&
            V1Mother->Pt() == V2Mother->Pt())
        {
          Double_t m = V2Mother->M();
          Double_t pt = V2Mother->Pt();
          FillHistogram("hmcLamGamAlpha", m, pt);
        }
      }

      if (NPart && PPart)
      {

        if (NPart->GetMother() == PPart->GetMother() && NPart->GetMother() != -1)
        {
          V0Part = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(NPart->GetMother()));
          if (V0Part)
          {
            if (V0Part->GetPdgCode() == 22)
            {
              isReallyPhoton = kTRUE;

              AliAODMCParticle *V0Mother = nullptr;
              if (V0Part->GetMother() != -1)
                V0Mother = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(TMath::Abs(V0Part->GetMother())));
              if (V0Mother)
              {
                if (TMath::Abs(V0Mother->GetPdgCode()) == 3212)
                {
                  Double_t m = V0Mother->M();
                  Double_t pt = V0Mother->Pt();
                  FillHistogram("hmcLamGam", m, pt);
                }
              } //Mother of Photon exists and is a (anti)Sigma0
            }
          } //Mother exists and is a Photon
        }   //Both Tracks have a common Mother
      }     //Both Tracks have matched MC Particle
    }       //End of isMCEvent
    //  ***************End of AOD V0 MC treatment*********************** */

    //  Reconstruct photon with TLorentzVector
    electron.SetXYZM(vecN(0), vecN(1), vecN(2), fElectronMass);
    positron.SetXYZM(vecP(0), vecP(1), vecP(2), fElectronMass);
    photon = electron + positron;

    // Calculate photon invariant mass with TL
    Double_t photonmass = photon.M();

    // Angle calculation
    Double_t deltatheta = thetaPos - thetaNeg;

    // Acceptance Cut
    if (TMath::Abs(trackN->Eta()) > fMaxDaughtEta)
      continue;
    if (TMath::Abs(trackP->Eta()) > fMaxDaughtEta)
      continue;

    // Check Track quality and reject poor qualty tracks
    if (nTPCClustNeg < fMinTPCClustDaught)
      continue;
    if (nTPCClustPos < fMinTPCClustDaught)
      continue;

    // Armenteros-Podolanski Cuts
    if (TMath::Abs(alpha) > fMaxalpha)
      continue;

    if (TMath::Abs(qt) > fMaxqt)
      continue;

    // Angle Cut
    if (TMath::Abs(totangle) > fMaxopenangle)
      continue;

    if (TMath::Abs(deltatheta) > fMaxdeltatheta)
      continue;

    // CPA Cut
    if (cosPointAngle < fMinV0CPA)
      continue;

    // PID Cut
    if (!isPhotonTPC)
      continue;

    //Radius Cut
    if (Vtxradius < fMinV0Radius)
      continue;

    if (Vtxradius > fMaxV0Radius)
      continue;

    // Inv. Mass Cut
    if (photonmass > fMaxphotonmass)
      continue;

    // Store Photon candidates after selection
    fV0ParticleIDArray.push_back(trackN->GetID());
    fV0ParticleIDArray.push_back(trackP->GetID());

    fConvPhotonArray.push_back(iV0);
    countPhotons++;

  } //End of V0 Loop

  //return;
} //End of FillV0PhotonArray()

////////////////////////////////////////////////////////////////////////////////
// Template to fill 1D histogram
void AliAnalysisTaskSigmaPCMPHOS::FillHistogram(const char *key, Double_t x)
{
  //FillHistogram
  TH1 *hist = dynamic_cast<TH1 *>(fOutputList->FindObject(key));
  if (hist)
    hist->Fill(x);
}

////////////////////////////////////////////////////////////////////////////////
// Template to fill 2D histogram
void AliAnalysisTaskSigmaPCMPHOS::FillHistogram(const char *key, Double_t x, Double_t y)
{
  //FillHistogram
  TH1 *th1 = dynamic_cast<TH1 *>(fOutputList->FindObject(key));
  if (th1)
    th1->Fill(x, y);
}

////////////////////////////////////////////////////////////////////////////////
// Template to fill 3D histogram
void AliAnalysisTaskSigmaPCMPHOS::FillHistogram(const char *key, Double_t x, Double_t y, Double_t z)
{
  //Fills 1D histograms with key
  TObject *obj = fOutputList->FindObject(key);

  TH2 *th2 = dynamic_cast<TH2 *>(obj);
  if (th2)
  {
    th2->Fill(x, y, z);
    return;
  }

  TH3 *th3 = dynamic_cast<TH3 *>(obj);
  if (th3)
  {
    th3->Fill(x, y, z);
    return;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Select hadrons based on cuts specified in the header file
void AliAnalysisTaskSigmaPCMPHOS::SelectHadrons()
{
  //currently observe only He3 mucleus, have to fime 3_H^(a)Lambda mass peak, have to store also \pi+-
  // Multiplicity and momentum distribution of tracks
  const Double_t massPip = 0.13957;
  const Double_t massK = 0.493677;
  const Double_t massP = 0.938272;
  const Double_t massHe3 = 3.016029;

  Int_t nTracks = fAODEvent->GetNumberOfTracks();

  if (!fTracksPp)
    fTracksPp = new TClonesArray("AliCaloPhoton", fAODEvent->GetNumberOfTracks());
  else
    fTracksPp->Clear();
  if (!fTracksPm)
    fTracksPm = new TClonesArray("AliCaloPhoton", fAODEvent->GetNumberOfTracks());
  else
    fTracksPm->Clear();

  Int_t inPp = 0;
  Int_t inPm = 0;

  for (Int_t i = 0; i < nTracks; i++)
  {

    AliAODTrack *track = static_cast<AliAODTrack *>(fAODEvent->GetTrack(i)); // track (type AliAODTrack) from event
    if (!track->IsHybridGlobalConstrainedGlobal() || TMath::Abs(track->Eta()) > 0.8)
      continue;
    AliVParticle *inEvHMain = dynamic_cast<AliVParticle *>(track);
    if (!track || !track->TestFilterBit(1))
      continue; // if we failed, skip this track

    Bool_t pidProton = kFALSE;
    Bool_t pidHe3 = kFALSE;
    Double_t nsigmaProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);

    if (nsigmaProton < 10.67)
      pidProton = kTRUE;
    if (pidProton)
    {
      AliCaloPhoton *pr;
      if (track->Charge() > 0)
      {

        pr = new ((*fTracksPp)[inPp++]) AliCaloPhoton(track->Px(), track->Py(), track->Pz(), TMath::Sqrt(massP * massP + track->P() * track->P()));
      }
      else
      {

        pr = new ((*fTracksPm)[inPm++]) AliCaloPhoton(track->Px(), track->Py(), track->Pz(), TMath::Sqrt(massP * massP + track->P() * track->P()));
      }
      Double_t dca = track->DCA();
      pr->SetDispBit(dca > 0.8);
      pr->SetWeight(dca);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Select Lambda hyperons based on cuts specified in the header file
void AliAnalysisTaskSigmaPCMPHOS::SelectLambda()
{
  //Select Lmbdas from V0

  if (!fLambda)
    fLambda = new TClonesArray("TLorentzVector", 100);
  else
    fLambda->Clear();

  Int_t nv0 = fAODEvent->GetNumberOfV0s();
  Int_t inLambda = 0;
  const Double_t massLambda = 1.115683;

  while (nv0--)
  {
    AliAODv0 *v0 = fAODEvent->GetV0(nv0);
    if (!v0)
    {
      continue;
    }

    //Use onfly only
    if (!v0->GetOnFlyStatus())
      continue;

    const AliAODTrack *ntrack1 = (AliAODTrack *)v0->GetDaughter(1);
    if (!AcceptTrack(ntrack1))
      continue;

    const AliAODTrack *ptrack1 = (AliAODTrack *)v0->GetDaughter(0);
    if (!AcceptTrack(ptrack1))
      continue;

    // Remove like-sign
    if (ntrack1->Charge() == ptrack1->Charge())
    {
      continue;
    }

    if (v0->Pt() == 0)
    {
      continue;
    }

    if (ntrack1->GetKinkIndex(0) > 0 || ptrack1->GetKinkIndex(0) > 0)
      continue;

    Double_t lV0Position[3];
    v0->GetXYZ(lV0Position);

    Double_t lV0Radius = TMath::Sqrt(lV0Position[0] * lV0Position[0] + lV0Position[1] * lV0Position[1]);
    if (lV0Radius > 180.)
      continue;

    //DCA V0 daughters
    Double_t dca = v0->DcaV0Daughters();
    if (dca < 0.06)
      continue;

    Double_t cpa = v0->CosPointingAngle(fAODEvent->GetPrimaryVertex());
    if (cpa < 0.993)
      continue;

    Double_t nSigmaPosPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrack1, AliPID::kPion));
    Double_t nSigmaNegPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrack1, AliPID::kPion));
    Double_t nSigmaPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrack1, AliPID::kProton));
    Double_t nSigmaNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrack1, AliPID::kProton));

    Bool_t isLambda = 0;
    Bool_t isLambdaBar = 0;

    //DCA pi+- >0.02 cm
    //SCA p > 0.05 cm
    Float_t xyNeg = v0->DcaNegToPrimVertex();
    Float_t xyPos = v0->DcaPosToPrimVertex();
    if (ntrack1->Charge() > 0)
    { //Lambda and proton
      isLambda = (nSigmaNegProton < 3.7) && (nSigmaPosPion < 3.8) && (TMath::Abs(xyPos) > 0.02) && (TMath::Abs(xyNeg) > 0.05);
      isLambdaBar = (nSigmaPosProton < 3.9) && (nSigmaNegPion < 4.2) && (TMath::Abs(xyNeg) > 0.02) && (TMath::Abs(xyPos) > 0.05);
    }
    else
    {
      isLambda = (nSigmaPosProton < 3.7) && (nSigmaNegPion < 3.8) && (TMath::Abs(xyNeg) > 0.02) && (TMath::Abs(xyPos) > 0.05);
      isLambdaBar = (nSigmaNegProton < 3.9) && (nSigmaPosPion < 4.2) && (TMath::Abs(xyPos) > 0.02) && (TMath::Abs(xyNeg) > 0.05);
    }

    if (isLambda)
      FillHistogram("hLambdaMass", v0->MassLambda(), v0->Pt());
    if (isLambdaBar)
      FillHistogram("hLambdaBarMass", v0->MassAntiLambda(), v0->Pt());

    if (isLambda && TMath::Abs(v0->MassLambda() - 1.115) > 0.005)
      continue;
    if (isLambdaBar && TMath::Abs(v0->MassAntiLambda() - 1.115) > 0.005)
      continue;

    if (v0->Pt() < 0.5)
      continue;

    //So far combine Lambda and AntiLambda

    if (isLambda || isLambdaBar)
    {
      new ((*fLambda)[inLambda++]) TLorentzVector(v0->Px(), v0->Py(), v0->Pz(), TMath::Sqrt(massLambda * massLambda + v0->P() * v0->P()));
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Reject or accept track
Bool_t AliAnalysisTaskSigmaPCMPHOS::AcceptTrack(const AliAODTrack *t)
{
  if (!t->IsOn(AliAODTrack::kTPCrefit))
    return kFALSE;
  Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2, 1);
  if (nCrossedRowsTPC < 70)
    return kFALSE;
  Int_t findable = t->GetTPCNclsF();
  if (findable <= 0)
    return kFALSE;
  if (nCrossedRowsTPC / findable < 0.8)
    return kFALSE;

  return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
// Select photons based on cuts specified in the header file
void AliAnalysisTaskSigmaPCMPHOS::SelectGamma()
{

  //Select gamma in PHOS

  Int_t inPHOS = 0, iPi0Merged = 0;
  if (fGamma)
    fGamma->Clear();
  else
    fGamma = new TClonesArray("AliCaloPhoton", 100);

  const AliAODVertex *esdVertex5 = fAODEvent->GetPrimaryVertex();

  Double_t vtx5[3] = {esdVertex5->GetX(), esdVertex5->GetY(), esdVertex5->GetZ()};

  Int_t multClust = fAODEvent->GetNumberOfCaloClusters();

  for (Int_t i = 0; i < multClust; i++)
  {
    AliAODCaloCluster *clu = fAODEvent->GetCaloCluster(i);
    if (clu->GetType() != AliVCluster::kPHOSNeutral)
      continue; // always continue, why?
    if (clu->E() > 1.500)
      continue; // Ok 9sep22
    if (clu->GetM02() < 0.2)
      continue; // ok 9sep22, stong cut

    FillHistogram("hClusterEnergy", clu->E());
    FillHistogram("hClusterTOFvsE", clu->GetTOF(), clu->E());

    TLorentzVector pv1;
    clu->GetMomentum(pv1, vtx5);
    if (inPHOS >= fGamma->GetSize())
    {
      fGamma->Expand(inPHOS + 50);
    }
    AliCaloPhoton *p = new ((*fGamma)[inPHOS++]) AliCaloPhoton(pv1.X(), pv1.Py(), pv1.Z(), pv1.E());

    // what means Set*Bit after the cuts? - means that passed that particular cut?
    p->SetDispBit(clu->Chi2() < 2.5 * 2.5);
    p->SetTOFBit((clu->GetTOF() > -50.e-9) && (clu->GetTOF() < 50.e-9));
    p->SetCPVBit(clu->GetEmcCpvDistance() > 2.5);
    p->SetCPV2Bit(clu->GetEmcCpvDistance() > 1.);
    p->SetPrimary(0); //no matched partner yet
  }
}

////////////////////////////////////////////////////////////////////////////////
// Return Pi 0 mass
Double_t AliAnalysisTaskSigmaPCMPHOS::Pi0Mass(Double_t /*pt*/)
{
  return 0.137;
}

////////////////////////////////////////////////////////////////////////////////
// Return Pi 0 width
Double_t AliAnalysisTaskSigmaPCMPHOS::Pi0Width(Double_t /*pt*/)
{
  return 0.012; //2sigma
}

////////////////////////////////////////////////////////////////////////////////
// Return eta meson mass
Double_t AliAnalysisTaskSigmaPCMPHOS::EtaMass(Double_t /*pt*/)
{
  return 0.555;
}

////////////////////////////////////////////////////////////////////////////////
// Return eta meson width
Double_t AliAnalysisTaskSigmaPCMPHOS::EtaWidth(Double_t /*pt*/)
{
  return 0.030; //2sigma
}

////////////////////////////////////////////////////////////////////////////////
// Pion Dispersion cut
Double_t AliAnalysisTaskSigmaPCMPHOS::PionDispCut(Double_t m02, Double_t m20, Double_t E)
{
  //Returns ditance to pi0 peak center in sigmas
  //No Disp cut for soft energies
  if (E < 25.)
    return 999;

  //Parameterization using single pi0 simulation
  Double_t longMpi = 1.857398e+00 + 1.208331e+01 * TMath::Exp(-4.977723e-02 * E);
  Double_t longSpi = 3.820707e-01 + 1.000542e+00 * TMath::Exp(-3.877147e-02 * E);
  Double_t shortMpi = 1.152118e+00 - 4.076138e-01 * TMath::Exp(-2.372902e-02 * E);
  Double_t shortSpi = 1.517538e-01 + 9.382205e+00 * TMath::Exp(-1.563037e-01 * E);
  Double_t powerNpi = 2.055773e+00 + 9.616408e+03 * TMath::Exp(-2.664167e-01 * E);

  Double_t dx = (m02 - longMpi) / longSpi;
  Double_t dy = (m20 - shortMpi) / shortSpi;

  //we have non-gaussian power, so re-calculate in Gaussian sigmas
  return TMath::Sign(TMath::Sqrt(TMath::Power(TMath::Abs(dx), powerNpi) + TMath::Power(TMath::Abs(dy), powerNpi)), dx);
}

////////////////////////////////////////////////////////////////////////////////
// Class Termination
void AliAnalysisTaskSigmaPCMPHOS::Terminate(Option_t *option){};

////////////////////////////////////////////////////////////////////////////////
// Return Armenteros-Podolansky plot elements
void AliAnalysisTaskSigmaPCMPHOS::GetArPod(Double_t pos[3], Double_t neg[3], Double_t moth[3], Double_t arpod[2])
{

  //see header file for documentation

  TVector3 momentumVectorPositiveKF(pos[0], pos[1], pos[2]);
  TVector3 momentumVectorNegativeKF(neg[0], neg[1], neg[2]);
  TVector3 vecV0(moth[0], moth[1], moth[2]);

  Float_t thetaV0pos = TMath::ACos((momentumVectorPositiveKF * vecV0) / (momentumVectorPositiveKF.Mag() * vecV0.Mag()));
  Float_t thetaV0neg = TMath::ACos((momentumVectorNegativeKF * vecV0) / (momentumVectorNegativeKF.Mag() * vecV0.Mag()));

  Float_t alfa = ((momentumVectorPositiveKF.Mag()) * TMath::Cos(thetaV0pos) - (momentumVectorNegativeKF.Mag()) * TMath::Cos(thetaV0neg)) /
                 ((momentumVectorPositiveKF.Mag()) * TMath::Cos(thetaV0pos) + (momentumVectorNegativeKF.Mag()) * TMath::Cos(thetaV0neg));

  Float_t qt = momentumVectorPositiveKF.Mag() * TMath::Sin(thetaV0pos);

  arpod[0] = qt;
  arpod[1] = alfa;
}

////////////////////////////////////////////////////////////////////////////////
// Calculate rapidity
Double_t AliAnalysisTaskSigmaPCMPHOS::Rapidity(Double_t pt, Double_t pz, Double_t m)
{
  // calculates rapidity keeping the sign in case E == pz
  Double_t energy = TMath::Sqrt(pt * pt + pz * pz + m * m);
  if (energy != TMath::Abs(pz))
    return 0.5 * TMath::Log((energy + pz) / (energy - pz));
  return TMath::Sign(1.e30, pz);
}