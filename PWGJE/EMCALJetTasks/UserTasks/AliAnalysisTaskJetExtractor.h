#ifndef ALIANALYSISTASKJETEXTRACTOR_H
#define ALIANALYSISTASKJETEXTRACTOR_H

/* Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//###############################################################################################################################################3
class AliRDHFJetsCutsVertex;
class AliEmcalJetTree;
class AliHFJetsTaggingVertex;

/**
 * \class AliAnalysisTaskJetExtractor
 * \brief Analysis task that implements AliEmcalJetTree to extract jets to a tree
 *
 * This task extracts jets (including constituents) to a tree for further
 * analysis, e.g. to create a ML dataset.
 * Usage: Task will by default extract all jets with constituents.
 * To adjust extraction percentages, use task->GetJetTree().AddExtractionPercentage(Float_t minPt, Float_t maxPt, Float_t percentage);
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, Yale
 * \date Jul 27, 2018
 */
// 

class AliAnalysisTaskJetExtractor : public AliAnalysisTaskEmcalJet {
 public:
  // ## Define class to store and compare secondary vertices (needed for matching of jet splittings & sec. vertices)
  struct SimpleSecondaryVertex
  {
    Int_t           fIndex;         // Index of sec. vertices in the array that we save to the tree
    Double_t        fLxy;           // decay length
    AliVParticle*   fDaughter1;     // daughter particle
    AliVParticle*   fDaughter2;     // daughter particle
    AliVParticle*   fDaughter3;     // daughter particle
  };

  AliAnalysisTaskJetExtractor();
  AliAnalysisTaskJetExtractor(const char *name);
  virtual ~AliAnalysisTaskJetExtractor();
  static AliAnalysisTaskJetExtractor* AddTaskJetExtractor(TString trackArray, TString clusterArray, TString jetArray, TString rhoObject, Double_t jetRadius, AliRDHFJetsCutsVertex* vertexerCuts, const char* taskNameSuffix);
  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);
  AliEmcalJetTree*            GetJetTree() {return fJetTree;}

  void                        SetSaveConstituents(Bool_t val) {fSaveConstituents = val; fInitialized = kFALSE;}
  void                        SetSaveConstituentsIP(Bool_t val) {fSaveConstituentsIP = val; fInitialized = kFALSE;}
  void                        SetSaveConstituentPID(Bool_t val) {fSaveConstituentPID = val; fInitialized = kFALSE;}
  void                        SetSaveJetShapes(Bool_t val) {fSaveJetShapes = val; fInitialized = kFALSE;}
  void                        SetSaveJetSplittings(Bool_t val) {fSaveJetSplittings = val; fInitialized = kFALSE;}
  void                        SetSaveMCInformation(Bool_t val) {fSaveMCInformation = val; fInitialized = kFALSE;}
  void                        SetSaveSecondaryVertices(Bool_t val) {fSaveSecondaryVertices = val; fInitialized = kFALSE;}
  void                        SetSaveTriggerTracks(Bool_t val) {fSaveTriggerTracks = val; fInitialized = kFALSE;}
  void                        SetSaveCaloClusters(Bool_t val) {fSaveCaloClusters = val; fInitialized = kFALSE;}
  void                        SetMCParticleArrayName(const char* name)            { fMCParticleArrayName = name; }
  void                        SetHadronMatchingRadius(Double_t val)               { fHadronMatchingRadius  = val; }
  void                        SetJetMatchingRadius(Double_t val)                  { fJetMatchingRadius = val; }
  void                        SetJetMatchingSharedPtFraction(Double_t val)        { fJetMatchingSharedPtFraction = val; }
  void                        SetSecondaryVertexMaxChi2(Double_t val   )          { fSecondaryVertexMaxChi2 = val; }
  void                        SetSecondaryVertexMaxDispersion(Double_t val)       { fSecondaryVertexMaxDispersion = val; }
  void                        SetCustomStartupScript(const char* path)            { fCustomStartupScript = path; }
  void                        SetVertexerCuts(AliRDHFJetsCutsVertex* val)         { fVertexerCuts = val; }
  void                        SetSetEmcalJetFlavour(Bool_t val)                   { fSetEmcalJetFlavour = val; }
  void                        SetEventPercentage(Double_t val)                    { fEventPercentage  = val; }
  void                        SetTruthLabelRange(Int_t min, Int_t max)            { fTruthMinLabel = min; fTruthMaxLabel = max; }
  void                        SetSaveTrackPDGCode(Bool_t val)                     { fSaveTrackPDGCode = val; }
  void                        SetRandomSeed(ULong_t val)                          { fRandomSeed  = val; }
  void                        SetRandomSeedCones(ULong_t val)                     { fRandomSeedCones  = val; }
  void                        SetNeedEmbedClusterContainer(Bool_t val)            { fNeedEmbedClusterContainer = val;}
  
  void                        SetEventCutTriggerTrack(Double_t minPt, Double_t maxPt, Int_t minLabel=-9999999, Int_t maxLabel=+9999999)
                                { fEventCut_TriggerTrackMinPt = minPt; fEventCut_TriggerTrackMaxPt = maxPt; fEventCut_TriggerTrackMinLabel = minLabel;
                                  fEventCut_TriggerTrackMaxLabel = maxLabel;}

 protected:
  Bool_t                      CreateControlHistograms();
  void                        ExecOnce();
  Bool_t                      Run();
  void                        FillTrackControlHistograms(AliVTrack* track);
  void                        FillEventControlHistograms();
  void                        FillJetControlHistograms(AliEmcalJet* jet);
  void                        CalculateJetShapes(AliEmcalJet* jet, Double_t& leSub_noCorr, Double_t& angularity, Double_t& momentumDispersion, Double_t& trackPtMean, Double_t& trackPtMedian);
  void                        GetTrueJetPtFraction(AliEmcalJet* jet, Double_t& truePtFraction, Double_t& truePtFraction_mcparticles);
  bool                        PerformGeometricalJetMatching(AliJetContainer& contBase, AliJetContainer& contTag, double maxDist);
  void                        GetMatchedJetObservables(AliEmcalJet* jet,  Double_t& detJetPt, Double_t& partJetPt, Double_t& detJetDistance, Double_t& partJetDistance, Double_t& detJetMass, Double_t& partJetMass, Double_t& detJetAngularity, Double_t& partJetAngularity, Double_t& detJetpTD, Double_t& partJetpTD);
  void                        DoJetMatching(); 
  void                        GetJetType(AliEmcalJet* jet, Int_t& typeHM, Int_t& typePM, Int_t& typeIC);
  Bool_t                      IsTriggerTrackInEvent();
  Bool_t                      IsTrackInCone(const AliVParticle* track, Double_t eta, Double_t phi, Double_t radius);
  Bool_t                      IsClusterInCone( TLorentzVector clusterMomentum, Double_t eta, Double_t phi,Double_t radius);
  Bool_t                      IsStrangeJet(AliEmcalJet* jet);
  void                        PrintConfig();
  void                        AddPIDInformation(const AliVParticle* particle, Float_t& sigITS, Float_t& sigTPC, Float_t& sigTOF, Float_t& sigTRD, Short_t& recoPID, Int_t& truePID);
  void                        GetTrackImpactParameters(const AliVVertex* vtx, const AliAODTrack* track, Float_t& d0, Float_t& d0cov, Float_t& z0, Float_t& z0cov);
  void                        ReconstructSecondaryVertices(const AliVVertex* primVtx, const AliEmcalJet* jet, std::vector<Float_t>& secVtx_X, std::vector<Float_t>& secVtx_Y, std::vector<Float_t>& secVtx_Z, std::vector<Float_t>& secVtx_Mass, std::vector<Float_t>& secVtx_Lxy, std::vector<Float_t>& secVtx_SigmaLxy, std::vector<Float_t>& secVtx_Chi2, std::vector<Float_t>& secVtx_Dispersion);
  void                        GetJetSplittings(AliEmcalJet* jet, std::vector<Float_t>& splittings_radiatorE, std::vector<Float_t>& splittings_kT, std::vector<Float_t>& splittings_theta, std::vector<Int_t>& splittings_secVtx_rank, std::vector<Int_t>& splittings_secVtx_index);

  // ################## CUTS AND SETTINGS: Should be set during task initialization
  Bool_t                      fSaveConstituents;                        ///< save arrays of constituent basic properties
  Bool_t                      fSaveConstituentsIP;                      ///< save arrays of constituent impact parameters
  Bool_t                      fSaveConstituentPID;                      ///< save arrays of constituent PID parameters
  Bool_t                      fSaveJetShapes;                           ///< save jet shapes
  Bool_t                      fSaveJetSplittings;                       ///< save jet splittings from iterative CA reclustering
  Bool_t                      fSaveMCInformation;                       ///< save MC information
  Bool_t                      fSaveSecondaryVertices;                   ///< save reconstructed sec. vertex properties
  Bool_t                      fSaveTriggerTracks;                       ///< save event trigger track
  Bool_t                      fSaveCaloClusters;                        ///< save calorimeter clusters
  
  Double_t                    fEventPercentage;                         ///< percentage (0, 1] which will be extracted
  Double_t                    fEventCut_TriggerTrackMinPt;              ///< Event requirement, trigger track min pT
  Double_t                    fEventCut_TriggerTrackMaxPt;              ///< Event requirement, trigger track max pT
  Int_t                       fEventCut_TriggerTrackMinLabel;           ///< Event requirement, trigger track min label (can be used to selected toy particles)
  Int_t                       fEventCut_TriggerTrackMaxLabel;           ///< Event requirement, trigger track max label (can be used to selected toy particles)
  Int_t                       fEventCut_TriggerTrackOrigin;             ///< Event requirement, trigger track origin (0: no cut, 1: from embedding, 2: not from embedding)
  Int_t                       fTruthMinLabel;                           ///< min track label to consider it as true particle
  Int_t                       fTruthMaxLabel;                           ///< max track label to consider it as true particle
  Double_t                    fHadronMatchingRadius;                    ///< Matching radius to search for beauty/charm hadrons around jet
  Double_t                    fJetMatchingRadius;                       ///< Matching radius for geometrically matching jets
  Double_t                    fJetMatchingSharedPtFraction;             ///< Shared pT fraction required in matching

  TString                     fMCParticleArrayName;                     ///< Array name of MC particles in event (mcparticles)
  Bool_t                      fNeedEmbedClusterContainer;               ///< If we need to get embedded cluster container (true for hybrid event)

  ULong_t                     fRandomSeed;                              ///< random seed
  ULong_t                     fRandomSeedCones;                         ///< random seed
  AliRDHFJetsCutsVertex*      fVertexerCuts;                            ///< Cuts used for the vertexer (given in add task macro)

  Double_t                    fSecondaryVertexMaxChi2;                  ///< Max chi2 of secondary vertex (others will be discarded)
  Double_t                    fSecondaryVertexMaxDispersion;            ///< Max dispersion of secondary vertex (others will be discarded)
  TString                     fCustomStartupScript;                     ///< Path to custom shell script that will be executed
  Bool_t                      fSetEmcalJetFlavour;                      ///< if set, the flavour property of the AliEmcalJets will be set
  Bool_t                      fSaveTrackPDGCode;                        ///< save PDG code instead of code defined for AOD pid

  // ################## BUFFERS: Variables set during execution time
  AliEmcalJetTree*            fJetTree;
  Double_t                    fEventWeight;                             ///< event weight for each event (implemented for JEWEL)
  Double_t                    fImpactParameter;                         ///< IP for each event (implemented for JEWEL)
  Int_t                       fMultiplicity;                            ///< Multiplicity (number tracks, also for multiple containers)
  std::vector<Float_t>        fTriggerTracks_Pt;                        ///< found trigger track pT
  std::vector<Float_t>        fTriggerTracks_Eta;                       ///< found trigger track eta
  std::vector<Float_t>        fTriggerTracks_Phi;                       ///< found trigger track phi
  TClonesArray*               fMCParticleArray;                         //!<! Array of MC particles in event (usually mcparticles)
  TRandom3*                   fRandomGenerator;                         //!<! Random number generator, used for event + jet efficiency
  TRandom3*                   fRandomGeneratorCones;                    //!<! Random number generator, used for random cones
  AliHFJetsTaggingVertex*     fVtxTagger;                               //!<! class for sec. vertexing
  Bool_t                      fIsEmbeddedEvent;                         ///< Set to true if at least one embedding container is added to this task
  Bool_t                      fDoDetLevelMatching;                      ///< Whether or not we do det level matching                                 
  Bool_t                      fDoPartLevelMatching;                     ///< Whether or not we do particle level matching 

  std::vector<SimpleSecondaryVertex> fSimpleSecVertices;  ///< Vector of secondary vertices

  // ################## HELPER FUNCTIONS
  Double_t                    GetDistance(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2)
  {
    Double_t deltaPhi = TMath::Min(TMath::Abs(phi1-phi2),TMath::TwoPi() - TMath::Abs(phi1-phi2));
    return TMath::Sqrt((eta1-eta2)*(eta1-eta2) + deltaPhi*deltaPhi);
  }
  void                        FillHistogram(const char * key, Double_t x);
  void                        FillHistogram(const char * key, Double_t x, Double_t y);
  void                        FillHistogram(const char * key, Double_t x, Double_t y, Double_t add);
  void                        FillHistogram3D(const char * key, Double_t x, Double_t y, Double_t z, Double_t add = 0);
  template <class T> T*       AddHistogram1D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis");
  template <class T> T*       AddHistogram2D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0,  const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");
  template <class T> T*       AddHistogram3D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0, Int_t zBins = 100, Double_t zMin = 0.0, Double_t zMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");

 private:
  AliAnalysisTaskJetExtractor(const AliAnalysisTaskJetExtractor&);            // not implemented
  AliAnalysisTaskJetExtractor &operator=(const AliAnalysisTaskJetExtractor&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetExtractor, 9) // Jet extraction task
  /// \endcond
};


//###############################################################################################################################################3
const Int_t kMaxNumConstituents = 500;
/**
 * \class AliEmcalJetTree
 * \brief Class managing creation of a tree containing jets
 *
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, Yale
 * \date Jul 24, 2018
 */
// 
class AliEmcalJetTree : public TNamed
{
  public:
    AliEmcalJetTree();
    AliEmcalJetTree(const char* name);

    // ############ SETTERS
    void            AddExtractionPercentage(Float_t minPt, Float_t maxPt, Float_t percentage)
    {
      fExtractionPercentagePtBins.push_back(minPt);
      fExtractionPercentagePtBins.push_back(maxPt);
      fExtractionPercentages.push_back(percentage);
    }
    void            ResetExtractionPercentages() {fExtractionPercentages.clear(); fExtractionPercentagePtBins.clear();}
    void            AddExtractionJetTypeHM(Int_t type) {fExtractionJetTypes_HM.push_back(type);}
    void            AddExtractionJetTypePM(Int_t type) {fExtractionJetTypes_PM.push_back(type);}

    void            InitializeTree(Bool_t saveCaloClusters, Bool_t saveMCInformation, Bool_t saveMatchedJets_Det, Bool_t saveMatchedJets_Part, Bool_t saveConstituents, Bool_t saveConstituentsIP, Bool_t saveConstituentPID, Bool_t saveJetShapes, Bool_t saveSplittings, Bool_t saveSecondaryVertices, Bool_t saveTriggerTracks);

    // ######################################
    Bool_t          AddJetToTree(AliEmcalJet* jet, Bool_t saveConstituents, Bool_t saveConstituentsIP, Bool_t saveCaloClusters, Double_t* vertex, Float_t rho, Float_t rhoMass, Float_t centrality, Int_t multiplicity, Long64_t eventID, Float_t magField);
    void            FillBuffer_SecVertices(std::vector<Float_t>& secVtx_X, std::vector<Float_t>& secVtx_Y, std::vector<Float_t>& secVtx_Z, std::vector<Float_t>& secVtx_Mass, std::vector<Float_t>& secVtx_Lxy, std::vector<Float_t>& secVtx_SigmaLxy, std::vector<Float_t>& secVtx_Chi2, std::vector<Float_t>& secVtx_Dispersion);
    void            FillBuffer_JetShapes(AliEmcalJet* jet, Double_t leSub_noCorr, Double_t angularity, Double_t momentumDispersion, Double_t trackPtMean, Double_t trackPtMedian);
    void            FillBuffer_Splittings(std::vector<Float_t>& splittings_radiatorE, std::vector<Float_t>& splittings_kT, std::vector<Float_t>& splittings_theta, Bool_t saveSecondaryVertices, std::vector<Int_t>& splittings_secVtx_rank, std::vector<Int_t>& splittings_secVtx_index);
    void            FillBuffer_PID(std::vector<Float_t>& trackPID_ITS, std::vector<Float_t>& trackPID_TPC, std::vector<Float_t>& trackPID_TOF, std::vector<Float_t>& trackPID_TRD, std::vector<Short_t>& trackPID_Reco, std::vector<Int_t>& trackPID_Truth);
    void            FillBuffer_MonteCarlo(Int_t motherParton, Int_t motherHadron, Int_t partonInitialCollision,
                                    Float_t matchedJetDistance_Det, Float_t matchedJetPt_Det, Float_t matchedJetMass_Det, Float_t matchedJetAngularity_Det, Float_t matchedJetpTD_Det,
                                    Float_t matchedJetDistance_Part, Float_t matchedJetPt_Part, Float_t matchedJetMass_Part, Float_t matchedJetAngularity_Part, Float_t matchedJetpTD_Part,
                                    Float_t truePtFraction, Float_t truePtFraction_PartLevel, Float_t ptHard, Float_t eventWeight, Float_t impactParameter);
    void            FillBuffer_ImpactParameters(std::vector<Float_t>& trackIP_d0, std::vector<Float_t>& trackIP_z0, std::vector<Float_t>& trackIP_d0cov, std::vector<Float_t>& trackIP_z0cov);
    void            FillBuffer_TriggerTracks(std::vector<Float_t>& triggerTrackPt, std::vector<Float_t>& triggerTrackDeltaEta, std::vector<Float_t>& triggerTrackDeltaPhi);
    // ######################################

    void            SetRandomGenerator(TRandom3* gen) {fRandomGenerator = gen;}
    std::vector<Float_t> GetExtractionPercentagePtBins() {return fExtractionPercentagePtBins;}
    std::vector<Float_t> GetExtractionPercentages() {return fExtractionPercentages;}
    std::vector<Int_t> GetExtractionJetTypes_HM() {return fExtractionJetTypes_HM;}
    std::vector<Int_t> GetExtractionJetTypes_PM() {return fExtractionJetTypes_PM;}
    TTree*          GetTreePointer() {return fJetTree;}

  private:
    TTree*          fJetTree;                             //!<! tree structure
    Bool_t          fInitialized;                         ///< init state of tree
    TRandom3*       fRandomGenerator;                     //!<! random generator

    // Option flags
    std::vector<Float_t> fExtractionPercentages;          ///< Percentages which will be extracted for a given pT bin
    std::vector<Float_t> fExtractionPercentagePtBins;     ///< pT-bins associated with fExtractionPercentages
    std::vector<Int_t>   fExtractionJetTypes_HM;          ///< Jet-types that will be extracted with this tree (hadron matching)
    std::vector<Int_t>   fExtractionJetTypes_PM;          ///< Jet-types that will be extracted with this tree (parton matching)

    // Buffers that will be added to the tree
    Float_t         fBuffer_JetPt;                        //!<! array buffer
    Float_t         fBuffer_JetEta;                       //!<! array buffer
    Float_t         fBuffer_JetPhi;                       //!<! array buffer
    Float_t         fBuffer_JetArea;                      //!<! array buffer
    Int_t           fBuffer_NumTracks;                    //!<! array buffer
    Int_t           fBuffer_NumClusters;                  //!<! array buffer

    Float_t         fBuffer_Event_BackgroundDensityMass;  //!<! array buffer
    Float_t         fBuffer_Event_BackgroundDensity;      //!<! array buffer
    Float_t         fBuffer_Event_Vertex_X;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Y;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Z;               //!<! array buffer
    Float_t         fBuffer_Event_Centrality;             //!<! array buffer
    Int_t           fBuffer_Event_Multiplicity;           //!<! array buffer
    Long64_t        fBuffer_Event_ID;                     //!<! array buffer
    Float_t         fBuffer_Event_MagneticField;          //!<! array buffer
    Float_t         fBuffer_Event_PtHard;                 //!<! array buffer
    Float_t         fBuffer_Event_Weight;                 //!<! array buffer
    Float_t         fBuffer_Event_ImpactParameter;        //!<! array buffer

    Float_t*        fBuffer_Track_Pt;                     //!<! array buffer
    Float_t*        fBuffer_Track_Eta;                    //!<! array buffer
    Float_t*        fBuffer_Track_Phi;                    //!<! array buffer
    Float_t*        fBuffer_Track_Charge;                 //!<! array buffer
    Int_t*          fBuffer_Track_Label;                  //!<! array buffer
    Float_t*        fBuffer_Track_ProdVtx_X;              //!<! array buffer
    Float_t*        fBuffer_Track_ProdVtx_Y;              //!<! array buffer
    Float_t*        fBuffer_Track_ProdVtx_Z;              //!<! array buffer

    Float_t*        fBuffer_Cluster_Pt;                   //!<! array buffer
    Float_t*        fBuffer_Cluster_E;                    //!<! array buffer
    Float_t*        fBuffer_Cluster_Eta;                  //!<! array buffer
    Float_t*        fBuffer_Cluster_Phi;                  //!<! array buffer
    Float_t*        fBuffer_Cluster_M02;                  //!<! array buffer
    Float_t*        fBuffer_Cluster_Time;                 //!<! array buffer
    Int_t*          fBuffer_Cluster_Label;                //!<! array buffer

    Float_t         fBuffer_Shape_Mass_NoCorr;            //!<! array buffer
    Float_t         fBuffer_Shape_Mass_DerivCorr_1;       //!<! array buffer
    Float_t         fBuffer_Shape_Mass_DerivCorr_2;       //!<! array buffer
    Float_t         fBuffer_Shape_pTD_DerivCorr_1;        //!<! array buffer
    Float_t         fBuffer_Shape_pTD_DerivCorr_2;        //!<! array buffer
    Float_t         fBuffer_Shape_LeSub_NoCorr;           //!<! array buffer
    Float_t         fBuffer_Shape_LeSub_DerivCorr;        //!<! array buffer
    Float_t         fBuffer_Shape_Angularity_NoCorr;      //!<! array buffer
    Float_t         fBuffer_Shape_Angularity_DerivCorr_1; //!<! array buffer
    Float_t         fBuffer_Shape_Angularity_DerivCorr_2; //!<! array buffer
    Float_t         fBuffer_Shape_Circularity_DerivCorr_1;//!<! array buffer
    Float_t         fBuffer_Shape_Circularity_DerivCorr_2;//!<! array buffer
    Float_t         fBuffer_Shape_Sigma2_DerivCorr_1;     //!<! array buffer
    Float_t         fBuffer_Shape_Sigma2_DerivCorr_2;     //!<! array buffer
    Float_t         fBuffer_Shape_NumTracks_DerivCorr;    //!<! array buffer
    Float_t         fBuffer_Shape_MomentumDispersion;     //!<! array buffer
    Float_t         fBuffer_Shape_TrackPtMean;            //!<! array buffer
    Float_t         fBuffer_Shape_TrackPtMedian;          //!<! array buffer

    Int_t           fBuffer_Jet_MC_MotherParton;          //!<! array buffer
    Int_t           fBuffer_Jet_MC_MotherHadron;          //!<! array buffer
    Int_t           fBuffer_Jet_MC_MotherIC;              //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchedDetLevelJet_Distance;   //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchedDetLevelJet_Pt;         //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchedDetLevelJet_Mass;       //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchedDetLevelJet_Angularity; //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchedDetLevelJet_pTD;        //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchedPartLevelJet_Distance;   //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchedPartLevelJet_Pt;         //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchedPartLevelJet_Mass;       //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchedPartLevelJet_Angularity; //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchedPartLevelJet_pTD;        //!<! array buffer
    Float_t         fBuffer_Jet_MC_TruePtFraction;               //!<! array buffer
    Float_t         fBuffer_Jet_MC_TruePtFraction_PartLevel;     //!<! array buffer

  

    Int_t           fBuffer_NumTriggerTracks;
    Int_t           fBuffer_NumSecVertices;
    Int_t           fBuffer_NumSplittings;

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetTree, 12) // Jet tree class
    /// \endcond
};

#endif
