/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskHypertritonKFTree_H
#define AliAnalysisTaskHypertritonKFTree_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

// includes added to play with KFParticle
#define HomogeneousField
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFPTrack.h"

class AliEventCuts;
class THnSparse;

class AliAnalysisTaskHypertritonKFTree : public AliAnalysisTaskSE  
{
public:
  AliAnalysisTaskHypertritonKFTree();
  AliAnalysisTaskHypertritonKFTree(const char *name);
  virtual ~AliAnalysisTaskHypertritonKFTree();
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  
  // Setter
  void SetRunAsData(Bool_t RunData) {kRunAsData=RunData;};
  void SetQA(Bool_t ActivateQA) {kDoQA=ActivateQA;};
  void SetIsMC(Bool_t UseMC) {kIsMC=UseMC;};
  void SetMCQA(Bool_t ActivateMCQA) {kDoMCQA=ActivateMCQA;};
  
private:
  void ProcessAOD();      // Track-by-track code without MC information
  void ProcessMC();       // Track-by-track code based on MC information
  Bool_t PassedBasicTrackQualityCuts (AliAODTrack* track); //Basic Track selection for the analysis
  Bool_t Helium3Selection (AliAODTrack* track); // Selection for 3He daughter
  Bool_t PionSelection (AliAODTrack* track); // Selection for π daughter
  Bool_t DaughterSelection (KFParticle kfpDaughter1, KFParticle kfpDaughter2);
  Bool_t HypertritonCandidateSelection (KFParticle kfpMother, KFParticle kfpDaughter1, KFParticle kfpDaughter2); // Selection for hypertriton candidate
  Double_t GetDCA (AliAODTrack *track , TString type);
  KFParticle CreateKFTrack(AliAODTrack *track, int pdgCode); // Converts AOD track to KFParticle
  KFVertex CreateKFVertex(const AliVVertex* vertex); // Converts AOD Vertex to KFVertex
  Double_t CalculatePointingAngle(KFParticle KFPart, KFVertex KFVtx); // calculates the pointing angle between the momentum of the particle and the vector connecting it to the given vertex
  
  void FillHe3Variables (AliAODTrack* track, KFParticle KFPart);
  void FillPionVariables (AliAODTrack* track, KFParticle KFPart);
  void FillDaughterVariables (KFParticle kfpDaughter1, KFParticle kfpDaughter2);
  
  AliEventCuts fEventCuts;        // Event cuts
  AliPIDResponse* fPIDResponse;   //! PID response
  KFVertex PrimVertex;
  Bool_t kRunAsData;              // True if analysis as data should be run
  Bool_t kDoQA;                   // True if QA should be stored (PID within ±5 sigma for pion and 3He)
  Bool_t kIsMC;                   // True if Monte Carlo simulation is analysed
  Bool_t kDoMCQA;                 // True if Monte Carlo QA should be stored (momentum correction, momentum ranges)
  
  AliAODEvent* fAODEvent;         //! input event
  AliMCEvent* fMCEvent;           //! corresponding mc event
  TList* fQAList;                 //! QA output list
  TList* fOutputList;             //! Output list
  
  TH1I* histoEventCentrality;     //! Count events within centrality percentile
  
  TTree* fCandidateTree;          //! Tree witth hypertriton candidates
  // Branches
  Float_t CentralityPercentile;                // Event centrality
  // Hypertrition candidate variables
  Float_t mass;                   // mass of the reconstructed candidate
  Float_t ErrorMass;              // mass of the reconstructed candidate
  Float_t p;                      // momentum of the reconstructed candidate
  Float_t pT;                     // transverse momentum of the reconstructed candidate
  Float_t Rapidity;               // Rapidity of the candidate
  Float_t massTopo;               // mass of the reconstructed candidate (with PV constraint)
  Float_t ErrorMassTopo;                // mass of the reconstructed candidate (with PV constraint)
  Float_t pTopo;                  // momentum of the reconstructed candidate (with PV constraint)
  Float_t pTTopo;                 // transverse momentum of the reconstructed candidate (with PV constraint)
  Float_t RapidityTopo;           // Rapidity of the candidate (with PV constraint)
  Int_t SignOfPair;               // Sign of the pair (LS or ULS) - for signal and background determination
  Float_t CosPointingAngle;       // Cosine of the pointing angle
  Float_t Chi2PerNDF;             // Chi2 per NDF for candidates
  Float_t Chi2PerNDFTopo;         // Chi2 per NDF for candidates (constrained to the primary vertex)
//  Float_t Chi2PerNDFMass;       // Chi2 per NDF for candidates (constrained to the mass 2.99131)
  Float_t DecayLength;            // Decay Length for candidates (constrained to the primary vertex)
  Float_t ErrorDecayLength;       // Error on the decay Length for candidates (constrained to the primary vertex)
  Float_t DecayLengthXY;          // Decay Length for candidates (constrained to the primary vertex) in xy
  Float_t ErrorDecayLengthXY;     // Error on the decay Length for candidates (constrained to the primary vertex) in xy
  Float_t DistanceToPV;           // Distance of closest approach to the primary vertex of the hypertriton candidate (KF information)
  Float_t DeviationFromPV;        // Deviation from the primary vertex of the hypertriton candidate (KF information)
  Float_t DistanceToPVXY;         // Distance of closest approach to the primary vertex of the hypertriton candidate (KF information) in xy
  Float_t DeviationFromPVXY;      // Deviation from the primary vertex of the hypertriton candidate (KF information) in xy
  
  // Daughter variables
  Float_t DistanceOfDaughters;    // Distance between the daughter tracks
  Float_t DeviationOfDaughters;   // chi2 deviation of the daughter tracks
  Float_t DistanceOfDaughtersXY;  // Distance between the daughter tracks in xy
  Float_t DeviationOfDaughtersXY; // chi2 deviation of the daughter tracks  in xy

//  pion
  Float_t pPion;                  // Momentum of the pion daughter
  Float_t pTPion;                 // Momentum of the pion daughter
  Float_t DCAPion;                // Distance of closest approach to the primary vertex of the pion daughter (AOD information)
  Float_t DistanceToPVPion;       // Distance of closest approach to the primary vertex of the pion daughter (KF information)
  Float_t DeviationFromPVPion;    // Deviation from the primary vertex of the pion daughter (KF information)
  Float_t DCAPionXY;              // Radial distance of closest approach to the primary vertex of the pion daughter (AOD information)
  Float_t DistanceToPVPionXY;     // Radial distance of closest approach to the primary vertex of the pion daughter (KF information)
  Float_t DeviationFromPVPionXY;  // Radial deviation from the primary vertex of the pion daughter (KF information)

  Int_t NClusterTPCPion;          // For Pion candidates
  Int_t NPIDClusterTPCPion;       // For Pion candidates
  
  Float_t TPCMomPion;             // For Pion candidates
  Float_t TPCnSigmaPion;          // For pion candidates

  Bool_t HasPointOnITSLayer0Pion;
  Bool_t HasPointOnITSLayer1Pion;
  Bool_t HasPointOnITSLayer2Pion;
  Bool_t HasPointOnITSLayer3Pion;
  Bool_t HasPointOnITSLayer4Pion;
  Bool_t HasPointOnITSLayer5Pion;
  
//  He3
  Float_t p3He;                   // Momentum of the 3He daughter
  Float_t pT3He;                  // Transverse momentum of the 3He daughter
  Float_t DCA3He;                 // Distance of closest approach to the primary vertex of the pion daughter (AOD information)
  Float_t DistanceToPV3He;        // Distance of closest approach to the primary vertex of the pion daughter (KF information)
  Float_t DeviationFromPV3He;     // Deviation from the primary vertex of the pion daughter (KF information)
  Float_t DCA3HeXY;               // Radial distance of closest approach to the primary vertex of the pion daughter (AOD information)
  Float_t DistanceToPV3HeXY;      // Radial distance of closest approach to the primary vertex of the pion daughter (KF information)
  Float_t DeviationFromPV3HeXY;   // Radial deviation from the primary vertex of the pion daughter (KF information)
  
  Int_t NClusterTPC3He;           // For 3He candidates
  Int_t NPIDClusterTPC3He;        // For 3He candidates
  
  Float_t TPCMom3He;              // For 3He candidates
  Float_t TPCnSigma3He;           // For 3He candidates
  Float_t TPCnSigma3H;            // For 3He candidates
  
  Bool_t HasPointOnITSLayer0He3;
  Bool_t HasPointOnITSLayer1He3;
  Bool_t HasPointOnITSLayer2He3;
  Bool_t HasPointOnITSLayer3He3;
  Bool_t HasPointOnITSLayer4He3;
  Bool_t HasPointOnITSLayer5He3;
    
  // MC output
  TTree* fGeneratedTreeMC;        //! Tree with generated hypertriton (denominator for efficiency)
  TTree* fCandidateTreeMC;        //! Tree with hypertriton candidates

  // Outputs for QA and corrections
  TH2F* fHistNsigmaTPCvsP3He;     //! Nsigma vs momentum (to check PID selection for 3He)
  TH2F* fHistNsigmaTPCvsPPion;    //! Nsigma vs momentum (to check PID selection for pions)

  TH2F* fHistPxTrueRecHe3;        //! MC p_x true vs reconstructed
  TH2F* fHistPyTrueRecHe3;        //! MC p_y true vs reconstructed
  TH2F* fHistPzTrueRecHe3;        //! MC p_z true vs reconstructed
  
  TH1F* fHistMomPion;             //! (MC) Momentum of pions from hypertriton
  TH1F* fHistMomHe3;              //! (MC) Momentum of He3 from hypertriton
  
  AliAnalysisTaskHypertritonKFTree(const AliAnalysisTaskHypertritonKFTree&); // not implemented
  AliAnalysisTaskHypertritonKFTree& operator=(const AliAnalysisTaskHypertritonKFTree&); // not implemented
  
  ClassDef(AliAnalysisTaskHypertritonKFTree, 1);
};

#endif
