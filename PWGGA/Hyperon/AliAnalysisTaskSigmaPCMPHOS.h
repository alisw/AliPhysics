// @(#)root/base:$Id$
// Authors: Alexander Borissov, Sergei Solokhin, Nikita Gladin    03/01/23

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
class AliAODCaloCluster;
class AliAODv0;

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

  void FillHistogram(const char *key, Double_t x);                         //! Fill 1D histogram with name key
  void FillHistogram(const char *key, Double_t x, Double_t y);             //! Fill 2D histogram with name key
  void FillHistogram(const char *key, Double_t x, Double_t y, Double_t z); //! Fill 3D histogram with name key

  void FillCutsHistogram(const char *key, Double_t x);             //! Fill 1D parameter distribution histograms
  void FillCutsHistogram(const char *key, Double_t x, Double_t y); //! Fill 2D parameter distribution histograms

  void FillMCHistogram(const char *key, Double_t x);             //! Fill 1D MC-related histograms
  void FillMCHistogram(const char *key, Double_t x, Double_t y); //! Fill 2D MC-related histograms

  void FillPHOSCutHistograms(const AliAODCaloCluster *clu, const Int_t i); //! Fill PHOS parameter distribution histograms
  void FillPCMCutHistograms(const Int_t i, const Double_t params[]);       //! Fill PCM parameter distribution histograms
  void FillLambdaCutHistograms(const Int_t i);                             //! Fill Lambda parameter distribution histograms

  //Functions for Event Processing
  void FillV0PhotonArray();                 //! Select PCM photons and fill arrays
  void SelectGamma();                       //! Select PHOS Photons and fill arrays
  void SelectLambda();                      //! Select Lambda hyperons and fill arrays
  void ProcessMC();                         //! Process Monte-Carlo generated events
  Bool_t AcceptTrack(const AliAODTrack *t); //! decide whether or not to accept track

private:
  AliAnalysisTaskSigmaPCMPHOS(const AliAnalysisTaskSigmaPCMPHOS &);            //! not implemented
  AliAnalysisTaskSigmaPCMPHOS &operator=(const AliAnalysisTaskSigmaPCMPHOS &); //! not implemented

  ClassDef(AliAnalysisTaskSigmaPCMPHOS, 2); //! analysis task

  TList *fOutputList; //! Output list which contains histograms for analysis
  TList *fCutsList;   //! Output list which contains histograms for cuts tuning
  TList *fMCList;     //! Output list which contains histograms related to MC particles

  AliAODEvent *fAODEvent;         //! AOD Event to be processed
  AliAODv0 *fV0;                  //! photon conversion vertex
  TClonesArray *fAODMCTrackArray; //! TClonesArray containing the MC Particles
  AliPIDResponse *fPIDResponse;   //! PID response object
  AliMCEvent *fMCEvent;           //! Monte-Carlo event pointer
  AliAODv0 *fAODv0;               //! photon conversion vertex

  AliAODTrack *electronCandidate; //! electron candidate track
  AliAODTrack *positronCandidate; //! positron candidate track

  AliAODTrack *ptrack1; //! positive Lambda daughter candidate track
  AliAODTrack *ntrack1; //! negative Lambda daughter candidate track

  std::vector<int> fConvPhotonArray; //! V0 Photon candidates

  std::vector<TLorentzVector> fMCPhotonArray; //! Monte-Carlo generated PCM photons array
  std::vector<TLorentzVector> fMCLambdaArray; //! Monte-Carlo generated Lambda array
  std::vector<TLorentzVector> fMCPHOSArray;   //! Monte-Carlo generated PHOS photons array

  //Global cuts

  static constexpr Double_t fMaxZvertex = 10.0;          //! maximum Z vertex coordinate cut
  std::vector<Double_t> fEtaRangeCuts = {0.5, 0.8, 1.2}; //! pseudorapidity bounds used to fill histograms
  static constexpr Double_t fMaxDaughtEta = 1.2;         //! maximum daughters pseudorapidity cut

  //ProcessMCParticles Cuts
  static constexpr Double_t fMaxMCEta = 1.2; //! maximum pseudorapidity cut for simulated events

  //PCM Photon Cuts
  static constexpr Double_t fMinPhotonPt = 0.1; //! minimum electron-positron pair transverse momentum cut
  static constexpr Double_t fMaxPhotonPt = 1.5; //! maximum electron-positron pair transverse momentum cut

  static constexpr Double_t fMaxSigmaElectron = 4; //! maximum devation for electron signal
  static constexpr Double_t fMaxSigmaPositron = 4; //! maximum devation for positron signal

  static constexpr Double_t fMinTPCClustDaught = 30; //! minimum photon number of TPC clusters crossed cut
  static constexpr Double_t fMaxNsigDaughtTPC = 100; //! maximum photon daughters number of standard deviations in TPC cut
  static constexpr Double_t fMaxAlphaPhoton = 0.7;   //! Armenteros-Podolansky plot assymmetry parameter cut for photons
  static constexpr Double_t fMaxQtPhoton = 0.025;    //! Armenteros-Podolansky plot momentum cut for photons
  static constexpr Double_t fMaxopenangle = 10;      //! maximum photon open angle cut
  static constexpr Double_t fMaxdeltatheta = 10;     //! maximum photon delta change cut
  static constexpr Double_t fMinV0CPA = 0.8;         //! minimum photon cosine of pointing angle cut
  static constexpr Double_t fMinV0Radius = 3;        //! minimum photon V0 radius cut
  static constexpr Double_t fMaxV0Radius = 220;      //! maximum photon V0 radius cut
  static constexpr Double_t fMaxphotonMass = 0.1;    //! maximum electron-positron pair mass cut

  static constexpr Double_t fMaxPhotonDaughterTracksDCA = 1.5; //! maximum DCA between electron and positron tracks
  static constexpr Double_t fMinPhotonDaughterDCA = 0.0;       //! minimum DCA between photon daughter and decay vertex

  //PHOS cuts
  static constexpr Double_t fMaxClusterEnergy = 1.0;      //! maximum PHOS cluster energy
  static constexpr Double_t fMinClusterEnergy = 0.1;      //! minimum PHOS cluster energy
  static constexpr Double_t fMaxClusterTOF = 400.e-9;     //! maximum PHOS cluster time of flight
  static constexpr Double_t fMaxClusterAssymmetry = 10.0; //! maximum PHOS cluster assymmetry
  static constexpr Double_t fMinClusterAssymmetry = 0.1;  //! maximum PHOS cluster assymmetry

  //Select Lambda cuts

  static constexpr Double_t fMinLambdaPt = 0.5; //! minimum Lambda transverse momentum cut
  static constexpr Double_t fMaxLambdaPt = 10;  //! maximum Lambda transverse momentum cut

  static constexpr Double_t fMinCPA = 0.995;           //! minimum cosine of pointing angle for Lambda
  static constexpr Double_t fMaxlV0Radius = 180;       //! maximum Lambda decay radius
  static constexpr Double_t fMaxMassDeviation = 0.005; //! maximum deviation of Lambda candidate from PDG mass

  static constexpr Double_t fMaxSigmaNegProton = 4; //! maximum devation for antiproton signal
  static constexpr Double_t fMaxSigmaPosPion = 4;   //! maximum devation for positive pion signal
  static constexpr Double_t fMaxSigmaPosProton = 4; //! maximum devation for proton signal
  static constexpr Double_t fMaxSigmaNegPion = 4;   //! maximum devation for negative pion signal

  static constexpr Double_t fMaxLambdaDaughterTracksDCA = 1.5; //! maximum DCA between proton and pion tracks
  static constexpr Double_t fMinLambdaDaughterDCA = 0.02;      //! minimum DCA between Lambda daughter and decay vertex

  static constexpr Double_t fMaxAlphaLambda = 0.95; //! Armenteros-Podolansky plot assymmetry parameter cut for Lambdas
  static constexpr Double_t fMinAlphaLambda = 0.2;  //! Armenteros-Podolansky plot assymmetry parameter cut for Lambdas
  static constexpr Double_t fMaxQtLambda = 0.14;    //! Armenteros-Podolansky plot maximum momentum cut for Lambdas
  static constexpr Double_t fMinQtLambda = 0.02;    //! Armenteros-Podolansky plot minimum momentum cut for Lambdas

  //AcceptTrack cuts
  static constexpr Double_t fMinCrossedRowsTPC = 70;   //! minimum number of TPC clusters crossed by track
  static constexpr Double_t fMinTPCCLusterRatio = 0.8; //!  minimum ratio of number of TPC clusters crossed by track to number of findable tracks

  AliPHOSGeometry *fPHOSGeo;            //! PHOS geometry
  AliTriggerAnalysis *fTriggerAnalysis; //! Trigger Analysis for Normalisation

  TClonesArray *fGamma;  //! array of PHOS photons
  TClonesArray *fPCM;    //! array of photons from the photon conversion method
  TClonesArray *fLambda; //! array of Lambda hyperons
  TList *fMixGamma;      //! list of mixed photon events
  TList *fMixLambda;     //! list of mixed Lambda hyperons
};
#endif