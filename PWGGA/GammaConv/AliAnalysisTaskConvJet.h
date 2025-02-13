#ifndef AliAnalysisTaskConvJet_H
#define AliAnalysisTaskConvJet_H
/**
 * \file AliAnalysisTaskConvJet.h
 *
 * \author Lizette Lamers <lizette.jacqueline.lamers@cern.ch>, Utrecht University
 * \author Mike Sas <mike.sas@cern.ch>, Utrecht University
 * \date Okt 4, 2018
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

/**
 * \class AliAnalysisTaskConvJet
 * \brief Implementation of a sample jet analysis task.
 *
 * This class in an implementation of a sample task for EMCal jet analysis.
 * It derives from AliAnalysisTaskEmcalJet.
 * It performs a simple analysis, producing jet spectra.
 */
class AliAnalysisTaskConvJet : public AliAnalysisTaskEmcalJet
{
 public:
  AliAnalysisTaskConvJet();
  AliAnalysisTaskConvJet(const char* name);
  virtual ~AliAnalysisTaskConvJet();

  void UserCreateOutputObjects();
  void Terminate(Option_t* option);

  static AliAnalysisTaskConvJet* AddTask_GammaConvJet(
    const char* ntracks = "usedefault",
    const char* nclusters = "usedefault",
    const char* ncells = "usedefault",
    const char* suffix = "",
    const double distToEMCBorder = 0.,
    const double distToSMEdges = 0.
    );

  Double_t GetNJets() { return fNJets; }
  std::vector<Double_t> GetVectorJetPt() { return fVectorJetPt; }
  std::vector<Double_t> GetVectorJetPx() { return fVectorJetPx; }
  std::vector<Double_t> GetVectorJetPy() { return fVectorJetPy; }
  std::vector<Double_t> GetVectorJetPz() { return fVectorJetPz; }
  std::vector<Double_t> GetVectorJetEta() { return fVectorJetEta; }
  std::vector<Double_t> GetVectorJetPhi() { return fVectorJetPhi; }
  std::vector<Double_t> GetVectorJetArea() { return fVectorJetR; }
  std::vector<Double_t> GetVectorJetNEF() { return fVectorJetNEF; }
  std::vector<Double_t> GetVectorJetNtracks() { return fVectorJetNCh; }
  std::vector<Double_t> GetVectorJetNclus() { return fVectorJetNClus; }

  // Jet constituents
  std::vector<AliVCluster*> GetJetClusters(int jet) { return fVecJetClusters[jet]; }
  std::vector<AliVParticle*> GetJetTracks(int jet) { return fVecJetTracks[jet]; }
  std::vector<AliVParticle*> GetTrueJetParticles(int jet) { return fVecTrueJetParticles[jet]; }
  double GetTrueJetLeadPartPt(int jet) const { return fVecTrueJetMaxPartPt[jet]; }
  int GetTrueJetLeadPartPDG(int jet) const { return fVecTrueJetMaxPartPDG[jet]; }

  Double_t GetTrueNJets() { return fTrueNJets; }
  std::vector<Double_t> GetTrueVectorJetPt() { return fTrueVectorJetPt; }
  std::vector<Double_t> GetTrueVectorJetPx() { return fTrueVectorJetPx; }
  std::vector<Double_t> GetTrueVectorJetPy() { return fTrueVectorJetPy; }
  std::vector<Double_t> GetTrueVectorJetPz() { return fTrueVectorJetPz; }
  std::vector<Double_t> GetTrueVectorJetEta() { return fTrueVectorJetEta; }
  std::vector<Double_t> GetTrueVectorJetPhi() { return fTrueVectorJetPhi; }
  std::vector<Double_t> GetTrueVectorJetArea() { return fTrueVectorJetR; }
  std::vector<Double_t> GetTrueVectorJetNPart() { return fTrueVectorJetNPart; }

  std::vector<int> GetTrueVectorJetParton() { return fTrueVectorJetParton; }
  std::vector<double> GetTrueVectorJetPartonPt() { return fTrueVectorJetPartonPt; }
  std::vector<double> GetTrueVectorJetPartonPx() { return fTrueVectorJetPartonPx; }
  std::vector<double> GetTrueVectorJetPartonPy() { return fTrueVectorJetPartonPy; }
  std::vector<double> GetTrueVectorJetPartonPz() { return fTrueVectorJetPartonPz; }

  // double GetTrueJetLeadingParticlePt(int index) {return fTrueVectorJetLeadingPartPt; }
  
  UInt_t GetAcceptanceType() { return fAccType; }
  UInt_t GetAcceptanceTypeMC() { return fAccTypeMC; }

  Double_t Get_Jet_Radius()
  {
    AliJetContainer* jetCont = 0;
    TIter next(&fJetCollArray);
    Double_t radius = -1;
    while ((jetCont = static_cast<AliJetContainer*>(next()))) {
      radius = jetCont->GetJetRadius();
    }
    return radius;
  }

  void FindPartonsJet(TClonesArray* arrMCPart);

  void SetDistToEMCBorder(const double tmp) {fDistToEMCBorder = tmp;}
  double GetDistToEMCBorder() const {return fDistToEMCBorder;}

  void SetDistToEMCSMEdge(const double tmp, const int mode) {fDistEMCSMEdge = tmp; fEMCSMEdgesMode = mode;}
  double GetDistToEMCSMEdge() const {return fDistEMCSMEdge;}

  void setWeightEnergyJets(const char * formula);

 protected:
  void ExecOnce();
  Bool_t FillHistograms();
  Bool_t Run();

  void DoJetLoop();
  std::tuple<double, int> GetLeadingPartPt(AliEmcalJet * jet, const bool isTrueJet = false);

  Double_t fNJets;                     // Number of reconstructed jets
  std::vector<Double_t> fVectorJetPt;  // Vector for the pt of the reconstructed jets
  std::vector<Double_t> fVectorJetPx;  // Vector for the px of the reconstructed jets
  std::vector<Double_t> fVectorJetPy;  // Vector for the py of the reconstructed jets
  std::vector<Double_t> fVectorJetPz;  // Vector for the pz of the reconstructed jets
  std::vector<Double_t> fVectorJetEta; // Vector for the eta of the reconstructed jets
  std::vector<Double_t> fVectorJetPhi; // Vector for the phi of the reconstructed jets
  std::vector<Double_t> fVectorJetR;   // Vector for the radius of the reconstructed jets
  std::vector<Double_t> fVectorJetNEF; // Vector for the neutral energy fraction of the reconstructed jets
  std::vector<Double_t> fVectorJetNClus;   // Vector for the number of calo clusters of the reconstructed jets
  std::vector<Double_t> fVectorJetNCh; // Vector for the number of tracks of the reconstructed jets

  Double_t fTrueNJets;                     // Number of true jets
  std::vector<Double_t> fTrueVectorJetPt;  // Vector for the pt of the true jets
  std::vector<Double_t> fTrueVectorJetPx;  // Vector for the px of the true jets
  std::vector<Double_t> fTrueVectorJetPy;  // Vector for the py of the true jets
  std::vector<Double_t> fTrueVectorJetPz;  // Vector for the pz of the true jets
  std::vector<Double_t> fTrueVectorJetEta; // Vector for the eta of the true jets
  std::vector<Double_t> fTrueVectorJetPhi; // Vector for the phi of the true jets
  std::vector<Double_t> fTrueVectorJetR;   // Vector for the radius of the true jets
  std::vector<Double_t> fTrueVectorJetNPart; // Vector for the number of particles contributing to this true jet

  std::vector<int> fTrueVectorJetParton;      // vector containing the mc stack id from the leading parton ("seed of the jet")
  std::vector<double> fTrueVectorJetPartonPt; // vector containing the pt of the leading parton ("seed of the jet")
  std::vector<double> fTrueVectorJetPartonPx; // vector containing the pt of the leading parton ("seed of the jet")
  std::vector<double> fTrueVectorJetPartonPy; // vector containing the pt of the leading parton ("seed of the jet")
  std::vector<double> fTrueVectorJetPartonPz; // vector containing the pt of the leading parton ("seed of the jet")

  // Jet constituents
  std::vector<std::vector<AliVCluster*>> fVecJetClusters;
  std::vector<std::vector<AliVParticle*>> fVecJetTracks;

  std::vector<std::vector<AliVParticle*>> fVecTrueJetParticles;

  std::vector<double> fVecTrueJetMaxPartPt;
  std::vector<double> fVecTrueJetMaxPartPDG;
  

  UInt_t fAccType;
  UInt_t fAccTypeMC;
  
  double fDistToEMCBorder;      // distance cut to the border of the EMCal, (if set to 0.4, jet with R=0.4 is fully on the EMCal)
  int fEMCSMEdgesMode;          // Mode 1 = cut away the edges, Mode 2 = keep the edges
  double fDistEMCSMEdge;        // distance cut to SM Edges. If jet axis is near that, cut it away (or the other way around)

  bool fApplyEnergyWeight;
  TF1* funcEnergyWeights;

 private:
  bool IsJetAccepted(const AliEmcalJet *jet);
  AliAnalysisTaskConvJet(const AliAnalysisTaskConvJet&);
  AliAnalysisTaskConvJet& operator=(const AliAnalysisTaskConvJet&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskConvJet, 20);
  /// \endcond
};
#endif
