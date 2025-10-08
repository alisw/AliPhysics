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
#include "TH3F.h"
#include "TLorentzVector.h"
#include "AliConvK0LambdaCuts.h"

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
    const double distToSMEdges = 0.,
    const bool addV0sToJet = false
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
  std::vector<TLorentzVector> GetJetClusters(int jet) { return fVecJetClusters[jet]; }
  std::vector<int> GetJetClustersLabel(int jet) { return fVecJetClusterLabel[jet]; }
  std::vector<float> GetJetClustersEFrac(int jet) { return fVecJetClusterEFrac[jet]; }
  std::vector<AliAODTrack*> GetJetTracks(int jet) { return fVecJetTracks[jet]; }
  std::vector<int> GetTrueJetParticles(int jet) { return fVecTrueJetParticles[jet]; }
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

  void setWeightEnergyJets(const char * formula, const int mode = 1);
  void SetMeasurablePart(TString str);

  void SetJetEnergyScaleAndResol(TString sJES, TString sJER, TString sJERCut){
    funcJES = new TF1("func_JES", sJES, 0, 1000);
    funcJER = new TF1("func_JER", sJER, 0, 1000);
    funcJERCut = new TF1("func_JERCut", sJERCut, 0, 1000);
  }
  double GetJES(double pt) const {
    return funcJES->Eval(pt);
  }
  double GetJER(double pt) const {
    return funcJER->Eval(pt);
  }
  double GetJERCut(double pt) const {
    if(!funcJERCut) return 1.e3;
    return funcJERCut->Eval(pt);
  }

  void SetMinRCutV0(double tmp) { fMinRV0 = tmp; }
  void SetMaxRCutV0(double tmp) { fMaxRV0 = tmp; }
  void SetMaxPtCutV0(double tmp) { fMaxPtCutV0 = tmp; }
  void SetMaxPtCutV0Leg(double tmp) { fMaxPtCutV0Leg = tmp; }
  void SetFracTPCClus(double tmp) { fMinFracTPCClusV0Leg = tmp; }
  void SetMaxAlphaV0(double tmp) { fMaxCutAlphaV0 = tmp; }
  void SetMaxQtV0(double tmp) { fMaxCutV0Qt = tmp; }
  void SetNSigmaCut(double tmp) { fnSigmaPID = tmp; }
  void SetCosPACut(double tmp) { fMinCosPA = tmp; }
  void SetMaxEtaCut(double tmp) { fMaxEtaCutVoLeg = tmp; }
  void SetSelectK0sMass(bool tmp) {fSelectK0sMass = tmp; }
  void SetSelectLambdaMass(bool tmp) {fSelectLambdaMass = tmp; }
  void SetSelectPhotonMass(bool tmp) {fSelectPhotonMass = tmp; }
  void SetV0FinderType(int tmp) {fV0FinderType = tmp; }
  void SetDoExtQA(bool tmp) { fDoFillExtendedHistos = tmp; }
  void SetK0LambdaCutsFromCutString(const TString &cutString, const bool SetDoLightOutput = false, const bool enableQAK0Lambda = false){
    if (fK0LambdaCuts) delete fK0LambdaCuts;
    fK0LambdaCuts = new AliConvK0LambdaCuts(cutString);
    if (SetDoLightOutput > 0)
      fK0LambdaCuts->SetLightOutput(SetDoLightOutput);
    fK0LambdaCuts->SetDoQA(enableQAK0Lambda); 
    fK0LambdaCuts->InitializeCutsFromCutString(cutString);
    fK0LambdaCuts->SetFillCutHistograms("K0LambdaCuts", true);
    fDoK0LambdaCutsAdvanced = true;
  }
  bool SetUseClusHadCorrEnergy(bool useHadCorr) {
    return fUseClusHadCorrEnergy = useHadCorr;
  }

  void SetAddV0sToJet(bool tmp) { fAddV0sToJet = tmp; }
  void AddV0sToJet(double weight = 1., const int isMC = 0);
  bool IsParticleInJet(const std::vector<double> &vecJetEta, const std::vector<double> &vecJetPhi, double partEta, double partPhi, int &matchedJet, double &RJetPi0Cand );
  TList *GetHistograms(){return fHistograms;}

 protected:
  void ExecOnce();
  Bool_t FillHistograms();
  Bool_t Run();

  void DoJetLoop();
  std::tuple<double, int> GetLeadingPartPt(AliEmcalJet * jet, const bool isTrueJet = false);
  bool IsNonMeasurable(const int pdg, const int charge);

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
  std::vector<std::vector<TLorentzVector>> fVecJetClusters;
  std::vector<std::vector<int>> fVecJetClusterLabel;
  std::vector<std::vector<float>> fVecJetClusterEFrac;
  std::vector<std::vector<AliAODTrack*>> fVecJetTracks;

  std::vector<std::vector<int>> fVecTrueJetParticles;

  std::vector<double> fVecTrueJetMaxPartPt;
  std::vector<double> fVecTrueJetMaxPartPDG;
  

  UInt_t fAccType;
  UInt_t fAccTypeMC;
  
  double fDistToEMCBorder;      // distance cut to the border of the EMCal, (if set to 0.4, jet with R=0.4 is fully on the EMCal)
  int fEMCSMEdgesMode;          // Mode 1 = cut away the edges, Mode 2 = keep the edges
  double fDistEMCSMEdge;        // distance cut to SM Edges. If jet axis is near that, cut it away (or the other way around)

  int fApplyEnergyWeight;       // 1: Apply constant energy shift for all jets, 2: Apply energy shift depending on the fraction of non-measurable particles
  TF1* funcEnergyWeights;
  std::vector<int> fVecMeasurable;

  TF1* funcJES;                 // function to parameterize jet energy scale
  TF1* funcJER;                 // function to parameterize jet energy resolution
  TF1* funcJERCut;              // function to parameterize jet energy resolution cut (used for jet matching for response)

  bool fAddV0sToJet;                     // flag if V0s should be added to the jet
  bool fV0sCurrEvtAdded;                 // flag to keep track if for the current event, the V0s were already added
  double fMinRV0;                        // minimum radius for V0
  double fMaxRV0;                        // maximum radius for V0
  double fMaxPtCutV0;                    // maximum Pt for V0
  double fMaxPtCutV0Leg;                 // maximum Pt for V0 leg
  double fMinFracTPCClusV0Leg;           // minimum fraction of TPC clusters for V0
  double fMaxCutAlphaV0;                 // maximum alpha for V0
  double fMaxCutV0Qt;                    // maximum Qt cut
  double fnSigmaPID;                     // nsigma cut on PID
  double fMinCosPA;                      // Cos PA cut
  double fMaxEtaCutVoLeg;                // Max Eta cut
  bool fSelectK0sMass;                   // Bool if V0s in K0s mass should be selected
  bool fSelectLambdaMass;                // Bool if V0s in Lambda mass should be selected
  bool fSelectPhotonMass;                // Bool if V0s in Photon mass should be selected
  int fV0FinderType;                     // Select V0finder type: 0-> both finders, 1->offlineFinder, 2, on-the-fly finder
  bool fDoFillExtendedHistos;            // Flag to fill extended QA histograms
  AliConvK0LambdaCuts* fK0LambdaCuts;    // Cuts for K0s and Lambda V0s
  bool fDoK0LambdaCutsAdvanced;          // Flag to use advanced K0s and Lambda cuts
  bool fUseClusHadCorrEnergy;            // Flag to use hadronic correction for cluster energy
  TList* fHistograms;                    //!  List containing all histograms
  TH2F* hJetEtaDiffWithV0;               //!  difference in eta between original and V0 included jets
  TH2F* hJetPhiDiffWithV0;               //!  difference in phi between original and V0 included jets
  TH2F* hJetPtDiffWithV0;                //!  difference in pt between original and V0 included jets
  TH1D* hV0Pt;                           //!  V0Pt that is added to the jet energy
  TH1D* hV0LegsPt;                       //!  Pt of V0 legs (each leg can only contribute once) that is added to the jet energy
  TH3F* fHistoV0TrueMotherLegPt;         //! Pt of V0 legs for different mother particles
  TH2F* fHistoV0TrueMotherLegPtTrue;     //! True Pt of V0 legs for different mother particles
  TH1F* fHistoV0LegHybrid;               //! Pt of V0 legs that are reconstructed as hybrid tracks
  TH2F* fHistoV0DaughterResolHybrid;     //! Resolution of V0 daughters that are reconstructed as hybrid tracks
  TH3F* fHistoV0TrueMotherLegPtHybrid;   //! Pt of V0 legs for different mother particles that are reconstructed as hybrid tracks
  TH2F* fHistoV0TrueMotherLegPtTrueHybrid; //! True Pt of V0 legs for different mother particles that are reconstructed as hybrid tracks
  TH2F* fHistoArmenterosV0;              //!  Armenteros-Podolanski plot for V0s added to jet
  TH3F* fHistoArmenterosV0Pt;            //!  Armenteros-Podolanski plot for V0s added to jet
  TH3F* fHistoArmenterosV0Source;        //!  Armenteros-Podolanski plot for V0s added to jet
  TH2F* fHistoMassPhoton;                //! Mass of photons with K0s expectation
  TH2F* fHistoV0DaughterResol;           //!  Resolution of V0 daughters
  TH2F* fHistoSourceV0Pt;                //! Reconstructed V0s for different sources (K0s, Lambda, photons, rest) vs pt
  TH2F* fHistoSourceV0DaughterPt;        //! Reconstructed V0s daughters for different sources (K0s, Lambda, photons, rest) vs pt
  TH2F* fHistoSourceV0DaughterTrackPt;   //! Reconstructed V0s daughter tracks that are reconstructed for different sources (K0s, Lambda, photons, rest) vs pt
  TH2F* fHistoV0Generated;               //! Generated V0 particles for different sources (K0s, Lambda, photons, rest) vs pt
  TH2F* fHistoV0GeneratedMeasDecay;      //! Generated V0 particles with corrected decay channel for different sources (K0s, Lambda, photons, rest) vs pt
  TH2F* fHistoV0DaughtersGenerated;      //! Generated V0 particle daughters with corrected decay channel for different sources (K0s, Lambda, photons, rest) vs pt
  TH3F* fHistoV0SourceVsPtVsJetPt;       //! Reconstructed V0s for different sources (K0s, Lambda, photons, rest) vs pt vs jet pt
  TH3F* fHistoV0SourceVsGenPtVsJetPt;    //! Reconstructed V0s for different sources (K0s, Lambda, photons, rest) vs gen pt vs jet pt
  TH3F* fHistoGenV0SourceVsPtVsJetPt;    //! Generated V0s for different sources (K0s, Lambda, photons, rest) vs pt vs jet pt
  TH3F* fHistoGenV0SourceMeasDecayVsPtVsJetPt; //! Generated V0s with corrected decay channel for different sources (K0s, Lambda, photons, rest) vs pt vs jet pt
  TH3F* fHistoGenV0SourceVsPtVsTrueJetPt;    //! Generated V0s for different sources (K0s, Lambda, photons, rest) vs pt vs true jet pt
  TH3F* fHistoGenV0SourceMeasDecayVsPtVsTrueJetPt; //! Generated V0s with corrected decay channel for different sources (K0s, Lambda, photons, rest) vs pt vs true jet pt
  TH3F* fHistoV0TrueMotherSourceVsPtVsJetPt;     //! Reconstructed V0s true mother for different sources (K0s, Lambda, photons, rest) vs pt vs jet pt
  TH3F* fHistoV0TrueMotherSourceVsPtVsJetPtHybrid; //! Reconstructed V0s true mother for different sources (K0s, Lambda, photons, rest) vs pt vs jet pt that are reconstructed as hybrid tracks
  TH2F* fHistoClusterFromPhotonInJetPt;     //! Pt of clusters that are identified as coming from photons in jet
  TH2F* fHistoClusterFromK0sInJetPt;          //! Pt of clusters that are identified as coming from K0s in jet
  TH2F* fHistoClusterFromK0sInJetPtK0sPt;    //! Pt of clusters that are identified as coming from K0s in jet vs K0s pt
  TH2F* fHistoClusterFromLambdaInJetPt;       //! Pt of clusters that are identified as coming from Lambda in jet
  TH2F* fHistoClusterFromLambdaInJetPtLambdaPt; //! Pt of clusters that are identified as coming from Lambda in jet vs Lambda pt
  TH2F* fHistoClusterFromPhotonInJetPtEFrac;   //! Energy fraction of clusters that are identified as coming from photons in jet
  TH2F* fHistoClusterFromK0sInJetPtEFrac;        //! Energy fraction of clusters that are identified as coming from K0s in jet
  TH2F* fHistoClusterFromLambdaInJetPtEFrac;  //! Energy fraction of clusters that are identified as coming from K0s in jet vs K0s pt


 private:
  bool IsJetAccepted(const AliEmcalJet *jet);
  AliAnalysisTaskConvJet(const AliAnalysisTaskConvJet&);
  AliAnalysisTaskConvJet& operator=(const AliAnalysisTaskConvJet&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskConvJet, 29);
  /// \endcond
};
#endif