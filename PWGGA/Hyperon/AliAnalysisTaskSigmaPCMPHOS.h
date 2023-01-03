// @(#)root/base:$Id$
// Authors: Alexander Borissov, Sergei Solokhin    03/01/23

#ifndef AliAnalysisTaskSigmaPCMPHOS_H
#define AliAnalysisTaskSigmaPCMPHOS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

class AliPHOSGeometry;
class AliTriggerAnalysis;
class AliPIDResponse;
class AliAODTrack;
class AliStack;
class AliMCEvent;

#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "AliStack.h"
#include "AliMCEvent.h"

class AliAnalysisTaskSigmaPCMPHOS : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskSigmaPCMPHOS();
  AliAnalysisTaskSigmaPCMPHOS(const char *name);
  virtual ~AliAnalysisTaskSigmaPCMPHOS();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void FillHistogram(const char *key, Double_t x);                         //Fill 1D histogram witn name key
  void FillHistogram(const char *key, Double_t x, Double_t y);             //Fill 2D histogram witn name key
  void FillHistogram(const char *key, Double_t x, Double_t y, Double_t z); //Fill 3D histogram witn name key

  //Functions for Event Processing
  void FillV0PhotonArray(); // Select V0 Photons and fill arrays (fConvPhotonArray)
  void SelectHadrons();
  void SelectGamma();
  void SelectLambda();
  void ProcessMC();
  void GetArPod(Double_t pos[3], Double_t neg[3], Double_t moth[3], Double_t arpod[2]);
  Double_t Rapidity(Double_t pt, Double_t pz, Double_t m);

  Bool_t AcceptTrack(const AliAODTrack *t);
  Double_t Pi0Mass(Double_t pt);
  Double_t Pi0Width(Double_t pt);
  Double_t EtaMass(Double_t pt);
  Double_t EtaWidth(Double_t pt);
  Double_t PionDispCut(Double_t m02, Double_t m20, Double_t E);

private:
  AliAnalysisTaskSigmaPCMPHOS(const AliAnalysisTaskSigmaPCMPHOS &);            // not implemented
  AliAnalysisTaskSigmaPCMPHOS &operator=(const AliAnalysisTaskSigmaPCMPHOS &); // not implemented

  ClassDef(AliAnalysisTaskSigmaPCMPHOS, 1); // PHOS analysis task

  TList *fOutputList;            //! Output list which contains all Histograms
  AliAODEvent *aodEvent;         //! AOD Event to be processed
  TClonesArray *AODMCTrackArray; //! TClonesArray containing the MC Particles
  AliPIDResponse *fPIDResponse;  //! PID response object
  AliMCEvent *mcEvent;           //! Monte-Carlo event flag
  AliStack *fStack;              //! Stack of simulated events

  std::vector<int> fOnFlyVector;       //! Save the track IDs used by the V0 Finders
  std::vector<int> fFinderVector;      //! to avoid double counting. Uses ESD IDs!
  std::vector<int> fV0ParticleIDArray; //! Save IDs of Particles found by any of the V0 Finders
  std::vector<int> fConvPhotonArray;   //! V0 Photon candidates

  Double_t cElectronMass;   //! Electron PDG mass
  Double_t cProtonMass;     //! Proton PDG mass
  Double_t cSigmaMass;      //! Sigma^0 PDG mass
  Double_t cPi0Mass;        //! Pi^0 PDG mass
  Double_t c;               //! speed of light
  Double_t primaryVtxPosX;  //! x position of the primary Vertex of the Event
  Double_t primaryVtxPosY;  //! y position of the primary Vertex of the Event
  Double_t primaryVtxPosZ;  //! z position of the primary Vertex of the Event
  UInt_t EventTriggers;     //! Triggers of the Event
  ULong64_t fGlobalEventID; //! Global ID of the Collision

  //FillProtonArray Cuts
  Double_t fMaxProtEta = 0.9;     //! maximum proton pseudorapidity cut
  Double_t fMinTPCClustProt = 60; //! minimum proton number of TPC clusters crossed cut
  Double_t fMaxNsigProtTPC = 3.5; //! maximum proton number of standard deviations in TPC cut
  Double_t fMaxNsigProtTOF = 5;   //! maximum proton number of standard deviations in TOF cut
  Double_t fMaxpOnlyTPCPID = 0.8;
  Double_t fMinProtpt = 0;  //! minimum proton transverse momentum cut
  Double_t fMaxProtpt = 15; //! maximum proton transverse momentum cut

  //ProcessMCParticles Cuts
  Double_t fMaxMCEta = 0.9; //! maximum pseudorapidity cut for simulated events

  //FillV0PhotonArray Cuts
  Double_t fMaxDaughtEta = 0.9;     //! maximum photon pseudorapidity cut
  Double_t fMinTPCClustDaught = 30; //! minimum photon number of TPC clusters crossed cut
  Double_t fMaxNsigDaughtTPC = 5;   //! maximum photon daughters number of standard deviations in TPC cut
  Double_t fMaxalpha = 1;           //! Armenteros-Podolansky plot assymmetry parameter cut
  Double_t fMaxqt = 0.05;           //! Armenteros-Podolansky plot momentum cut
  Double_t fMaxopenangle = 0.3;     //! maximum photon open angle cut
  Double_t fMaxdeltatheta = 0.1;    //! maximum photon delta change cut
  Double_t fMinV0CPA = 0.8;         //! minimum photon cosine of pointing angle cut
  Double_t fMinV0Radius = 3;        //! minimum photon V0 radius cut
  Double_t fMaxV0Radius = 220;      //! maximum photon V0 radius cut
  Double_t fMaxphotonmass = 0.1;    //! maximum electron-positron pair mass cut

  AliPHOSGeometry *fPHOSGeo;            //! PHOS geometry
  Int_t aodEventCounter;                //! number of analyzed events
  AliTriggerAnalysis *fTriggerAnalysis; //! Trigger Analysis for Normalisation

  TClonesArray *fGamma;     //! array of photons
  TClonesArray *fPCM;       //! array of photons from the photon conversion method
  TClonesArray *fPi0;       //! array of pi^0 mesons
  TClonesArray *fPi0Merged; //! array of merged pi^0 mesons
  TClonesArray *fTracksPim; //! array of negative pion tracks
  TClonesArray *fTracksPip; //! array of negative pion tracks
  TClonesArray *fTracksPp;  //! array of proton tracks
  TClonesArray *fTracksPm;  //! array of antiproton tracks
  TClonesArray *fLambda;    //! array of Lambda hyperons
  TList *fMixGamma;         //! list of mixed photon events
  TList *fMixTracksPp;      //! list of mixed proton tracks
  TList *fMixTracksPm;      //! list of mixed antiproton tracks
  TList *fMixLambda;        //! list of mixed Lambda hyperons
};
#endif