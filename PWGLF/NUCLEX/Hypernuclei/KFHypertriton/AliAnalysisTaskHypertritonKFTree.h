/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskHypertritonKFTree_H
#define AliAnalysisTaskHypertritonKFTree_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

// includes added for KFParticle
#define HomogeneousField ///homogenous field in z direction
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFPTrack.h"

#include <vector>
using namespace std;            /// std namespace: so you can do things like 'cout'

class AliEventCuts;
class AliESDInputHandler;
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
  void SetQA(Bool_t ActivateQA) {kDoQA=ActivateQA;};
  void SetIsMC(Bool_t UseMC) {kIsMC=UseMC;};
  
  void SetRun2Body(Bool_t Activate2Body) {kRun2BodyDecay = Activate2Body;};
  void SetRun3Body(Bool_t Activate3Body) {kRun3BodyDecay = Activate3Body;};
  
private:
  void ProcessESD();      /// Track-by-track code without MC information
  void TwoBodyDecay(const vector<Int_t>& He3CandTrackId, const vector<Int_t>& PionCandTrackId);                         /// Implementation for 2-body decay
  void ThreeBodyDecay(const vector<Int_t>& DeuteronCandTrackId, const vector<Int_t>& ProtonCandTrackId, const vector<Int_t>& PionCandTrackId);  /// Implementation for 3-body decay
  
  void ProcessMC();       /// Track-by-track code based on MC information
  void TwoBodyDecayMC(const vector<Int_t>& ESDTrackId, Int_t Id3He, Int_t IdPion);                         /// MC implementation for 2-body decay
  void ThreeBodyDecayMC(const vector<Int_t>& ESDTrackId, Int_t IdDeuteron, Int_t IdProton, Int_t IdPion);  /// MC implementation for 3-body decay
  Bool_t PassedEventSelection();
  Bool_t PassedBasicTrackQualityCuts(AliESDtrack* track); //Basic Track selection for the analysis
  
  /// Two body decay daugthers
  Bool_t Helium3Selection(AliESDtrack* track);   /// Selection for 3He daughter
  Bool_t PionSelection(AliESDtrack* track);     /// Selection for π daughter (2- and 3-body)
  /// Three body decay daughters
  Bool_t DeuteronSelection(AliESDtrack* track); /// Selection for deuteron daughter
  Bool_t ProtonSelection(AliESDtrack* track);   /// Selection for proton daughter
  
  Bool_t DaughterSelection();
  Bool_t HypertritonCandidateSelection(KFParticle kfpMother, Bool_t TwoBody = true); // Selection for hypertriton candidate
  Double_t GetDCA(AliESDtrack *track , TString type);
  KFParticle CreateKFTrack(AliESDtrack *track, int pdgCode); /// Converts ESD track to KFParticle
  KFVertex CreateKFVertex(const AliVVertex* vertex);         /// Converts ESD Vertex to KFVertex
  Float_t CalculatePointingAngle(KFParticle KFPart, KFVertex KFVtx); // Calculates the pointing angle between the momentum of the particle and the vector connecting it to the given vertex
  
  void FillHe3Variables(AliESDtrack* track, KFParticle KFPart);
  void FillPionVariables(AliESDtrack* track, KFParticle KFPart);
  
  void FillDeuteronVariables(AliESDtrack* track, KFParticle KFPart);
  void FillProtonVariables(AliESDtrack* track, KFParticle KFPart);
  
  void FillDaughterVariables(KFParticle kfpDaughter1, KFParticle kfpDaughter2);                           /// 2-body decay version
  void FillDaughterVariables(KFParticle kfpDeuteron, KFParticle kfpProton, KFParticle kfpPion);  /// 3-body decay version
  
  void FillDistanceToSececondaryVertex(KFParticle kfpHelium, KFParticle kfpPion, KFParticle kfpMother);
  void FillDistanceToSececondaryVertex(KFParticle kfpDeuteron, KFParticle kfpProton, KFParticle kfpPion, KFParticle kfpMother);
  
  
  AliEventCuts fEventCuts;        /// Event cuts
  AliTimeRangeCut fTimeRangeCut;
  AliPIDResponse* fPIDResponse;   //! PID response
  KFVertex PrimVertex;
  
  /// use global sdt::vector container produced once
  /// container for found 3He candidates
  vector<Int_t> He3CandTrackIdPos;
  vector<Int_t> He3CandTrackIdNeg;
  /// container for found pion candidates
  vector<Int_t> PionCandTrackIdPos;
  vector<Int_t> PionCandTrackIdNeg;
  /// container for found deuteron candidates
  vector<Int_t> DeuteronCandTrackIdPos;
  vector<Int_t> DeuteronCandTrackIdNeg;
  /// container for found proton candidates
  vector<Int_t> ProtonCandTrackIdPos;
  vector<Int_t> ProtonCandTrackIdNeg;
  
  /// container to store in between
  vector<KFParticle> HyperTritonCandidates;
  vector<int> DeuteronDaughter;
  vector<int> ProtonDaughter;
  
  Bool_t kRun2BodyDecay;          /// Switch to turn  on/off 2-body decay part
  Bool_t kRun3BodyDecay;          /// Switch to turn  on/off 3-body decay part
  Bool_t kDoQA;                   /// True if QA should be stored
  Bool_t kIsMC;                   /// True if Monte Carlo simulation is analysed
  
  AliESDtrackCuts* fESDtrackCuts; //! track selection for ESD track candidates
  TList* fQAList;                 //! QA output list
  TList* fOutputList;             //! Output list
  AliMCEvent* fMCEvent;           //! corresponding mc event
  
  TH1F* hNumberOfEvents;          //! Count events at different track selection steps
  TH1I* histoEventCentrality;     //! Count events within centrality percentile
  
  TTree* CandidateTree;          //! Tree with hypertriton candidates from 2-body decay
  
  /// Branches
  Float_t CentralityPercentile;   /// Event centrality
  /// Hypertrition candidate variables
  Float_t mass;                   /// mass of the reconstructed candidate
  Float_t ErrorMass;              /// error on the mass of the reconstructed candidate
  Float_t px;                     /// momentum of the reconstructed candidate in x direction
  Float_t py;                     /// momentum of the reconstructed candidate in y direction
  Float_t pz;                     /// momentum of the reconstructed candidate in z direction
  Float_t Rapidity;               /// Rapidity of the candidate
  Int_t Charge;                   /// Charge of the candidate
  
  Float_t Chi2PerNDF;             /// Chi2 per NDF for candidates
  Float_t CosPointingAngle;       /// Cosine of the pointing angle
  Float_t DistanceToPV;           /// Distance of closest approach to the primary vertex of the hypertriton candidate (KF information)
  Float_t DeviationFromPV;        /// Deviation from the primary vertex of the hypertriton candidate (KF information)
  Float_t DistanceToPVXY;         /// Distance of closest approach to the primary vertex of the hypertriton candidate (KF information) in xy
  Float_t DeviationFromPVXY;      /// Deviation from the primary vertex of the hypertriton candidate (KF information) in xy
  
  /// properties after topological constraint to PV
  Float_t massTopo;               /// mass of the reconstructed candidate (with PV constraint)
  Float_t ErrorMassTopo;          /// error of the mass of the reconstructed candidate (with PV constraint)
  Float_t pxTopo;                 /// momentum of the reconstructed candidate in x direction (with PV constraint)
  Float_t pyTopo;                 /// momentum of the reconstructed candidate in y direction (with PV constraint)
  Float_t pzTopo;                 /// momentum of the reconstructed candidate in z direction (with PV constraint)
  Float_t RapidityTopo;           /// Rapidity of the candidate (with PV constraint)
  Float_t Chi2PerNDFTopo;         /// Chi2 per NDF for candidates (constrained to the primary vertex)
  Float_t CosPointingAngleTopo;   /// Cosine of the pointing angle  (with PV constraint)
  Float_t DistanceToPVTopo;       /// Distance of closest approach to the primary vertex of the hypertriton candidate (KF information)  (with PV constraint)
  Float_t DeviationFromPVTopo;    /// Deviation from the primary vertex of the hypertriton candidate (KF information)  (with PV constraint)
  Float_t DistanceToPVXYTopo;      /// Distance of closest approach to the primary vertex of the hypertriton candidate (KF information) in xy  (with PV constraint)
  Float_t DeviationFromPVXYTopo;  /// Deviation from the primary vertex of the hypertriton candidate (KF information) in xy  (with PV constraint)
  Float_t DecayLength;            /// Decay Length for candidates (constrained to the primary vertex)
  Float_t ErrorDecayLength;       /// Error on the decay Length for candidates (constrained to the primary vertex)
  Float_t DecayLengthXY;          /// Decay Length for candidates (constrained to the primary vertex) in xy
  Float_t ErrorDecayLengthXY;     /// Error on the decay Length for candidates (constrained to the primary vertex) in xy
  
  /// Daughter variables
  Float_t DistanceOfDaughters;    /// Distance between the daughter tracks
  Float_t DeviationOfDaughters;   /// chi2 deviation of the daughter tracks
  Float_t DistanceOfDaughtersXY;  /// Distance between the daughter tracks in xy
  Float_t DeviationOfDaughtersXY; /// chi2 deviation of the daughter tracks  in xy
  
  ///  pion
  Float_t pxPion;                 /// momentum of the pion dauther in x direction
  Float_t pyPion;                 /// momentum of the pion dauther in y direction
  Float_t pzPion;                 /// momentum of the pion dauther in z direction
  Int_t ChargePion;               /// Charge of pion
  Float_t DCAPion;                /// Distance of closest approach to the primary vertex of the pion daughter (ESD information)
  //  Float_t DistanceToPVPion;       /// Distance of closest approach to the primary vertex of the pion daughter (KF information)
  Float_t DeviationFromPVPion;    /// Deviation from the primary vertex of the pion daughter (KF information)
  Float_t DCAPionXY;              /// Radial distance of closest approach to the primary vertex of the pion daughter (ESD information)
  //  Float_t DistanceToPVPionXY;     /// Radial distance of closest approach to the primary vertex of the pion daughter (KF information)
  Float_t DeviationFromPVPionXY;  /// Radial deviation from the primary vertex of the pion daughter (KF information)
  Int_t NCrossedRowsTPCPion;          /// For Pion candidates
  Int_t NPIDClusterTPCPion;       /// For Pion candidates
  Float_t TPCMomPion;             /// For Pion candidates
  Float_t TPCnSigmaPion;          /// For pion candidates
  Bool_t HasPointOnITSLayer0Pion;
  Bool_t HasPointOnITSLayer1Pion;
  Bool_t HasPointOnITSLayer2Pion;
  Bool_t HasPointOnITSLayer3Pion;
  Bool_t HasPointOnITSLayer4Pion;
  Bool_t HasPointOnITSLayer5Pion;
  Int_t PIDForTrackingPion;       /// Pion would be 2 (see AliPID.h)
  Float_t DistanceToSecVertPion;  /// distance of the pion track to the secondary vertex (KF information)
  Float_t DeviationToSecVertPion; /// deviation of the pion track from the secondary vertex (KF information)
  
  ///  He3
  Float_t pxHe;                   /// momentum of the He dauther in x direction
  Float_t pyHe;                   /// momentum of the He dauther in y direction
  Float_t pzHe;                   /// momentum of the He dauther in z direction
  Int_t ChargeHe;                 /// Charge of helium-3
  Float_t DCA3He;                 /// Distance of closest approach to the primary vertex of the pion daughter (ESD information)
  //  Float_t DistanceToPV3He;        /// Distance of closest approach to the primary vertex of the pion daughter (KF information)
  Float_t DeviationFromPV3He;     /// Deviation from the primary vertex of the pion daughter (KF information)
  Float_t DCA3HeXY;               /// Radial distance of closest approach to the primary vertex of the pion daughter (ESD information)
  //  Float_t DistanceToPV3HeXY;      /// Radial distance of closest approach to the primary vertex of the pion daughter (KF information)
  Float_t DeviationFromPV3HeXY;   /// Radial deviation from the primary vertex of the pion daughter (KF information)
  Int_t NCrossedRowsTPC3He;           /// For 3He candidates
  Int_t NPIDClusterTPC3He;        /// For 3He candidates
  Float_t TPCMom3He;              /// For 3He candidates
  Float_t TPCnSigma3He;           /// For 3He candidates
  Float_t TPCnSigma3H;            /// For 3He candidates
  Bool_t HasPointOnITSLayer0He3;
  Bool_t HasPointOnITSLayer1He3;
  Bool_t HasPointOnITSLayer2He3;
  Bool_t HasPointOnITSLayer3He3;
  Bool_t HasPointOnITSLayer4He3;
  Bool_t HasPointOnITSLayer5He3;
  Int_t PIDForTrackingHe3;       /// Helium-3 would be 7 (see AliPID.h)
  Float_t DistanceToSecVertHe;   /// distance of the 3He track to the secondary vertex (KF information)
  Float_t DeviationToSecVertHe;  /// deviation of the 3He track from the secondary vertex (KF information)
  
  /// three body decay
  TTree* CandidateTree_3Body;          //! Tree with hypertriton candidates from 3-body decay
  
  /// additional branches for the 3-body reconstruction
  Float_t mass_Deuteron_Proton;    /// mass of the reconstructed candidate (only deuteron and proton daughter)
  Float_t mass_Proton_Pion;        /// mass of the reconstructed candidate build up from proton and pion daughter (rejection of Lambdas possible)
  Float_t CosPointingAngle_Deuteron_Proton;       /// Cosine of the pointing angle (only deuteron and proton daughter)
  
  /// Daughter variables
  Float_t DistanceOfDaughters_Deuteron_Proton;
  Float_t DeviationOfDaughters_Deuteron_Proton;
  Float_t DistanceOfDaughtersXY_Deuteron_Proton;
  Float_t DeviationOfDaughtersXY_Deuteron_Proton;
  
  Float_t DistanceOfDaughters_Deuteron_Pion;
  Float_t DeviationOfDaughters_Deuteron_Pion;
  Float_t DistanceOfDaughtersXY_Deuteron_Pion;
  Float_t DeviationOfDaughtersXY_Deuteron_Pion;
  
  Float_t DistanceOfDaughters_Proton_Pion;
  Float_t DeviationOfDaughters_Proton_Pion;
  Float_t DistanceOfDaughtersXY_Proton_Pion;
  Float_t DeviationOfDaughtersXY_Proton_Pion;
  
  Float_t OpeningAngle_Pion_Proton;
  Float_t OpeningAngle_Pion_Deuteron;
  Float_t OpeningAngle_Proton_Deuteron;
  
  ///  Deuteron
  Float_t pxDeuteron;                      /// momentum of the deuteron dauther in x direction
  Float_t pyDeuteron;                      /// momentum of the deuteron dauther in y direction
  Float_t pzDeuteron;                      /// momentum of the deuteron dauther in z direction
  Int_t ChargeDeuteron;                    /// Charge of deuteron
  Float_t DCADeuteron;                /// Distance of closest approach to the primary vertex of the Deuteron daughter (ESD information)
  //  Float_t DistanceToPVDeuteron;       /// Distance of closest approach to the primary vertex of the Deuteron daughter (KF information)
  Float_t DeviationFromPVDeuteron;    /// Deviation from the primary vertex of the Deuteron daughter (KF information)
  Float_t DCADeuteronXY;              /// Radial distance of closest approach to the primary vertex of the Deuteron daughter (ESD information)
  //  Float_t DistanceToPVDeuteronXY;     /// Radial distance of closest approach to the primary vertex of the Deuteron daughter (KF information)
  Float_t DeviationFromPVDeuteronXY;  /// Radial deviation from the primary vertex of the Deuteron daughter (KF information)
  Int_t NCrossedRowsTPCDeuteron;          /// For Deuteron candidates
  Int_t NPIDClusterTPCDeuteron;       /// For Deuteron candidates
  Float_t TPCMomDeuteron;             /// For Deuteron candidates
  Float_t TPCnSigmaDeuteron;          /// For Deuteron candidates
  Float_t TOFnSigmaDeuteron;          /// For Deuteron candidates
  Bool_t HasPointOnITSLayer0Deuteron;
  Bool_t HasPointOnITSLayer1Deuteron;
  Bool_t HasPointOnITSLayer2Deuteron;
  Bool_t HasPointOnITSLayer3Deuteron;
  Bool_t HasPointOnITSLayer4Deuteron;
  Bool_t HasPointOnITSLayer5Deuteron;
  Int_t PIDForTrackingDeuteron;       /// Deuteron would be 5 (see AliPID.h)
  Float_t DistanceToSecVertDeuteron;
  Float_t DeviationToSecVertDeuteron;
  
  /// Proton
  Float_t pxProton;                      /// momentum of the proton dauther in x direction
  Float_t pyProton;                      /// momentum of the  proton dauther in y direction
  Float_t pzProton;                      /// momentum of the proton dauther in z direction
  Int_t ChargeProton;                    /// Charge of proton
  Float_t DCAProton;                /// Distance of closest approach to the primary vertex of the Proton daughter (ESD information)
  //  Float_t DistanceToPVProton;       /// Distance of closest approach to the primary vertex of the Proton daughter (KF information)
  Float_t DeviationFromPVProton;    /// Deviation from the primary vertex of the Proton daughter (KF information)
  Float_t DCAProtonXY;              /// Radial distance of closest approach to the primary vertex of the Proton daughter (ESD information)
  //  Float_t DistanceToPVProtonXY;     /// Radial distance of closest approach to the primary vertex of the Proton daughter (KF information)
  Float_t DeviationFromPVProtonXY;  /// Radial deviation from the primary vertex of the Proton daughter (KF information)
  Int_t NCrossedRowsTPCProton;          /// For Proton candidates
  Int_t NPIDClusterTPCProton;       /// For Proton candidates
  Float_t TPCMomProton;             /// For Proton candidates
  Float_t TPCnSigmaProton;          /// For Proton candidates
  Float_t TOFnSigmaProton;          /// For Proton candidates
  Bool_t HasPointOnITSLayer0Proton;
  Bool_t HasPointOnITSLayer1Proton;
  Bool_t HasPointOnITSLayer2Proton;
  Bool_t HasPointOnITSLayer3Proton;
  Bool_t HasPointOnITSLayer4Proton;
  Bool_t HasPointOnITSLayer5Proton;
  Int_t PIDForTrackingProton;       /// Deuteron would be 4 (see AliPID.h)
  Float_t DistanceToSecVertProton;
  Float_t DeviationToSecVertProton;
  
  /// MC output
  Float_t pxMC;             /// momentum of the reconstructed candidate in x direction from MC
  Float_t pyMC;             /// momentum of the reconstructed candidate in y direction from MC
  Float_t pzMC;             /// momentum of the reconstructed candidate in z direction from MC
  
  Float_t pxHeMC;           /// momentum of the He dauther in x direction from MC
  Float_t pyHeMC;           /// momentum of the He dauther in y direction from MC
  Float_t pzHeMC;           /// momentum of the He dauther in z direction from MC
  
  Float_t pxPionMC;         /// momentum of the pion dauther in x direction from MC
  Float_t pyPionMC;         /// momentum of the pion dauther in y direction from MC
  Float_t pzPionMC;         /// momentum of the pion dauther in z direction from MC
  
  Float_t pxProtonMC;       /// momentum of the proton dauther in x direction from MC
  Float_t pyProtonMC;       /// momentum of the  proton dauther in y direction from MC
  Float_t pzProtonMC;       /// momentum of the proton dauther in z direction from MC
  
  Float_t pxDeuteronMC;     /// momentum of the deuteron dauther in x direction from MC
  Float_t pyDeuteronMC;     /// momentum of the deuteron dauther in y direction from MC
  Float_t pzDeuteronMC;     /// momentum of the deuteron dauther in z direction from MC
  
  /// Variance of the momenta (only for MC to calculate the pull)
  Float_t pxVariance;           /// momentum of the reconstructed candidate in x direction from MC
  Float_t pyVariance;           /// momentum of the reconstructed candidate in y direction from MC
  Float_t pzVariance;           /// momentum of the reconstructed candidate in z direction from MC
  
  Float_t pxHeVariance;         /// momentum of the He dauther in x direction from MC
  Float_t pyHeVariance;         /// momentum of the He dauther in y direction from MC
  Float_t pzHeVariance;         /// momentum of the He dauther in z direction from MC
  
  Float_t pxPionVariance;       /// momentum of the pion dauther in x direction from MC
  Float_t pyPionVariance;       /// momentum of the pion dauther in y direction from MC
  Float_t pzPionVariance;       /// momentum of the pion dauther in z direction from MC
  
  Float_t pxProtonVariance;     /// momentum of the proton dauther in x direction from MC
  Float_t pyProtonVariance;     /// momentum of the  proton dauther in y direction from MC
  Float_t pzProtonVariance;     /// momentum of the proton dauther in z direction from MC
  
  Float_t pxDeuteronVariance;   /// momentum of the deuteron dauther in x direction from MC
  Float_t pyDeuteronVariance;   /// momentum of the deuteron dauther in y direction from MC
  Float_t pzDeuteronVariance;   /// momentum of the deuteron dauther in z direction from MC
  
  Float_t pMC;                  /// momentum of the candidate  (generated tree to calculate efficiency)
  Float_t pTMC;                 /// transverse momentum of the  candidate (generated tree to calculate efficiency)
  Float_t RapidityMC;           /// Rapidity of the candidate
  Float_t DecayLengthMC;        /// Decay Length for candidates (constrained to the primary vertex)
  Float_t DecayLengthXYMC;      /// Decay Length for candidates in XY (constrained to the primary vertex)
  
  Float_t xSecVertexMC;         /// Position of the decay vertex of the hypertriton candidate (MC truth)
  Float_t ySecVertexMC;         /// Position of the decay vertex of the hypertriton candidate (MC truth)
  Float_t zSecVertexMC;         /// Position of the decay vertex of the hypertriton candidate (MC truth)
  
  Float_t xSecVertex;           /// Position of the decay vertex of the hypertriton candidate
  Float_t ySecVertex;           /// Position of the decay vertex of the hypertriton candidate
  Float_t zSecVertex;           /// Position of the decay vertex of the hypertriton candidate
  Float_t xSecVertexVariance;   /// Variance of the position of the decay vertex of the hypertriton candidate
  Float_t ySecVertexVariance;   /// Variance of the position of the decay vertex of the hypertriton candidate
  Float_t zSecVertexVariance;   /// Variance of the position of the decay vertex of the hypertriton candidate
  
  Int_t NumberOfDeuteronProtonCandidates; /// Number of valid deuteron proton candidates
  
  TTree* GeneratedTreeMC;           //! Tree with generated hypertriton (denominator for efficiency)
  
  TTree* GeneratedTreeMC_3Body;     //! Tree with hypertriton candidates from 3-body decay
  
  // Outputs for QA and corrections
  TH2F* fHistNsigmaTPCvsP3He;       //! Nsigma TPC vs momentum (to check PID selection for 3He)
  TH2F* fHistNsigmaTPCvsPPion;      //! Nsigma TPC vs momentum (to check PID selection for pions)
  
  TH2F* fHistNsigmaTPCvsPDeuteron;  //! Nsigma TPC vs momentum (to check PID selection for deuteron)
  TH2F* fHistNsigmaTPCvsPProton;    //! Nsigma TPC vs momentum (to check PID selection for Proton)
  
  TH2F* fHistNsigmaTOFvsPDeuteron;  //! Nsigma TOF vs momentum (to check PID selection for deuteron)
  TH2F* fHistNsigmaTOFvsPProton;    //! Nsigma TOF vs momentum (to check PID selection for Proton)
  
  TH2F* fHistpxTrueRecHe3;          //! MC p_x true vs reconstructed
  TH2F* fHistpyTrueRecHe3;          //! MC p_y true vs reconstructed
  TH2F* fHistpzTrueRecHe3;          //! MC p_z true vs reconstructed
  
  TH1F* fHistMomPion;               //! (MC) Momentum of pions from hypertriton
  TH1F* fHistMomHe3;                //! (MC) Momentum of He3 from hypertriton
  
  TH1F* fHistMomDeuteron;           //! (MC) Momentum of deuteron from hypertriton
  TH1F* fHistMomProton;             //! (MC) Momentum of proton from hypertriton
  TH1F* fHistMomPion_3Body;         //! (MC) Momentum of pions from hypertriton
  
  TH1I* fHistNDeuteronsPerEvent;    //! Number of deuterons per event to learn if we gain having a separate loop for pions
  TH1I* fHistNDeuteronProtonCandidatesPerEvent;    //! Number of deuteron+proton candidates per event to learn if we gain having a separate loop for pions
  
  AliAnalysisTaskHypertritonKFTree(const AliAnalysisTaskHypertritonKFTree&); // not implemented
  AliAnalysisTaskHypertritonKFTree& operator=(const AliAnalysisTaskHypertritonKFTree&); // not implemented
  
  ClassDef(AliAnalysisTaskHypertritonKFTree, 1);
};

#endif
