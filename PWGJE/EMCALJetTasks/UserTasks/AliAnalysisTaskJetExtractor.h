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

  AliAnalysisTaskJetExtractor();
  AliAnalysisTaskJetExtractor(const char *name);
  virtual ~AliAnalysisTaskJetExtractor();
  static AliAnalysisTaskJetExtractor* AddTaskJetExtractor(TString trackArray, TString jetArray, TString rhoObject, Double_t jetRadius, TString configFile, AliRDHFJetsCutsVertex* vertexerCuts, const char* taskNameSuffix);
  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);
  AliEmcalJetTree*            GetJetTree() {return fJetTree;}

  void                        ActivateTrueJetMatching(const char* arrayName, const char* rhoName) 
                                {fTruthJetsArrayName = arrayName; fTruthJetsRhoName = rhoName;}
  void                        SetTrueJetRhoMassName(const char* rhoName) 
                                {fTruthJetsRhoMassName = rhoName;}
  void                        SetTrueParticleArrayName(const char* partArrayName) 
                                {fTruthParticleArrayName = partArrayName;}
  void                        SetHadronMatchingRadius(Double_t val)               { fHadronMatchingRadius  = val; }
  void                        SetTrueJetMatchingRadius(Double_t val)              { fTrueJetMatchingRadius = val; }
  void                        SetSecondaryVertexMaxChi2(Double_t val   )          { fSecondaryVertexMaxChi2 = val; }
  void                        SetSecondaryVertexMaxDispersion(Double_t val)       { fSecondaryVertexMaxDispersion = val; }
  void                        SetCalculateSecondaryVertices(Bool_t val)           { fCalculateSecondaryVertices = val; }
  void                        SetVertexerCuts(AliRDHFJetsCutsVertex* val)         { fVertexerCuts = val; }
  void                        SetSetEmcalJetFlavour(Bool_t val)                   { fSetEmcalJetFlavour = val; }
  void                        SetEventPercentage(Double_t val)                    { fEventPercentage  = val; }
  void                        SetTruthLabelRange(Int_t min, Int_t max)            { fTruthMinLabel = min; fTruthMaxLabel = max; }
  void                        SetSaveTrackPDGCode(Bool_t val)                     { fSaveTrackPDGCode = val; }
  void                        SetRandomSeed(ULong_t val)                          { fRandomSeed  = val; }
  void                        SetRandomSeedCones(ULong_t val)                     { fRandomSeedCones  = val; }

  void                        SetEventCutTriggerTrack(Double_t minPt, Double_t maxPt, Int_t minLabel=-9999999, Int_t maxLabel=+9999999)
                                { fEventCut_TriggerTrackMinPt = minPt; fEventCut_TriggerTrackMaxPt = maxPt; fEventCut_TriggerTrackMinLabel = minLabel;
                                  fEventCut_TriggerTrackMaxLabel = maxLabel;}

 protected:
  Bool_t                      CreateControlHistograms();
  void                        ExecOnce();
  Bool_t                      Run();
  Bool_t                      IsTriggerTrackInEvent();
  void                        FillTrackControlHistograms(AliVTrack* track);
  void                        FillEventControlHistograms();
  void                        FillJetControlHistograms(AliEmcalJet* jet);
  void                        AddPIDInformation(AliVParticle* particle, Float_t& sigITS, Float_t& sigTPC, Float_t& sigTOF, Float_t& sigTRD, Short_t& recoPID, Int_t& truePID);
  Bool_t                      IsTrackInCone(AliVParticle* track, Double_t eta, Double_t phi, Double_t radius);
  Double_t                    GetTrueJetPtFraction(AliEmcalJet* jet);
  Double_t                    GetMatchedTrueJetObservables(AliEmcalJet* jet, Double_t& matchedJetPt, Double_t& matchedJetMass);
  void                        GetJetType(AliEmcalJet* jet, Int_t& typeHM, Int_t& typePM, Int_t& typeIC);
  Bool_t                      IsStrangeJet(AliEmcalJet* jet);
  void                        PrintConfig();

  void                        GetTrackImpactParameters(const AliVVertex* vtx, AliAODTrack* track, Float_t& d0, Float_t& d0cov, Float_t& z0, Float_t& z0cov);
  void                        ReconstructSecondaryVertices(const AliVVertex* primVtx, const AliEmcalJet* jet, std::vector<Float_t>& secVtx_X, std::vector<Float_t>& secVtx_Y, std::vector<Float_t>& secVtx_Z,
                                  std::vector<Float_t>& secVtx_Mass, std::vector<Float_t>& secVtx_Lxy, std::vector<Float_t>& secVtx_SigmaLxy, std::vector<Float_t>& secVtx_Chi2, std::vector<Float_t>& secVtx_Dispersion);

  AliEmcalJetTree*            fJetTree;
  AliHFJetsTaggingVertex*     fVtxTagger;                               //!<! class for sec. vertexing
  Double_t                    fHadronMatchingRadius;                    ///< Matching radius to search for beauty/charm hadrons around jet
  Double_t                    fTrueJetMatchingRadius;                   ///< Matching radius to true jet
  Double_t                    fSecondaryVertexMaxChi2;                  ///< Max chi2 of secondary vertex (others will be discarded)
  Double_t                    fSecondaryVertexMaxDispersion;            ///< Max dispersion of secondary vertex (others will be discarded)
  Bool_t                      fCalculateSecondaryVertices;              ///< Calculate the secondary vertices (instead of loading)
  AliRDHFJetsCutsVertex*      fVertexerCuts;                            ///< Cuts used for the vertexer (given in add task macro)
  Bool_t                      fSetEmcalJetFlavour;                      ///< if set, the flavour property of the AliEmcalJets will be set
  Double_t                    fEventPercentage;                         ///< percentage (0, 1] which will be extracted
  Int_t                       fTruthMinLabel;                           ///< min track label to consider it as true particle
  Int_t                       fTruthMaxLabel;                           ///< max track label to consider it as true particle
  Bool_t                      fSaveTrackPDGCode;                        ///< save PDG code instead of code defined for AOD pid
  ULong_t                     fRandomSeed;                              ///< random seed
  ULong_t                     fRandomSeedCones;                         ///< random seed

  // ################## EVENT CUTS
  Double_t                    fEventCut_TriggerTrackMinPt;              ///< Event requirement, trigger track min pT
  Double_t                    fEventCut_TriggerTrackMaxPt;              ///< Event requirement, trigger track max pT
  Int_t                       fEventCut_TriggerTrackMinLabel;           ///< Event requirement, trigger track min label (can be used to selected embedded particles)
  Int_t                       fEventCut_TriggerTrackMaxLabel;           ///< Event requirement, trigger track max label (can be used to selected embedded particles)

  // ################## BASIC EVENT VARIABLES
  AliJetContainer            *fJetsCont;                                //!<! Jets
  AliParticleContainer       *fTracksCont;                              //!<! Tracks
  TClonesArray*               fTruthParticleArray;                      //!<! Array of MC particles in event (mcparticles)
  TString                     fTruthJetsArrayName;                      ///< Array name for particle-level jets
  TString                     fTruthJetsRhoName;                        ///< Array name for particle-level rho
  TString                     fTruthJetsRhoMassName;                    ///< Array name for particle-level rho for mass
  TString                     fTruthParticleArrayName;                  ///< Array name of MC particles in event (mcparticles)
  TRandom3*                   fRandomGenerator;                         //!<! Random number generator, used for event + jet efficiency
  TRandom3*                   fRandomGeneratorCones;                    //!<! Random number generator, used for random cones
  Double_t                    fEventWeight;                             ///< event weight for each event (only works for JEWEL)
  Double_t                    fImpactParameter;                         ///< IP for each event (only works for JEWEL)


  std::vector<Float_t>        fTriggerTracks_Pt;                         ///< If trigger track requi. is used -> save pT
  std::vector<Float_t>        fTriggerTracks_Eta;                        ///< If trigger track requi. is used -> save eta
  std::vector<Float_t>        fTriggerTracks_Phi;                        ///< If trigger track requi. is used -> save phi


  // ################## HISTOGRAM HELPER FUNCTIONS
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
  ClassDef(AliAnalysisTaskJetExtractor, 1) // Jet extraction task
  /// \endcond
};


//###############################################################################################################################################3
static std::vector<Float_t> DEFAULT_VECTOR_FLOAT;
static std::vector<Short_t> DEFAULT_VECTOR_SHORT;
static std::vector<Int_t>   DEFAULT_VECTOR_INT;
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

    void            AddExtractionJetTypeHM(Int_t type) {fExtractionJetTypes_HM.push_back(type);}
    void            AddExtractionJetTypePM(Int_t type) {fExtractionJetTypes_PM.push_back(type);}

    void            ResetExtractionPercentages() {fExtractionPercentages.clear(); fExtractionPercentagePtBins.clear();}

    void            SetSaveEventProperties(Bool_t val) {fSaveEventProperties = val; fInitialized = kFALSE;}
    void            SetSaveConstituents(Bool_t val) {fSaveConstituents = val; fInitialized = kFALSE;}
    void            SetSaveConstituentsIP(Bool_t val) {fSaveConstituentsIP = val; fInitialized = kFALSE;}
    void            SetSaveConstituentPID(Bool_t val) {fSaveConstituentPID = val; fInitialized = kFALSE;}
    void            SetSaveJetShapes(Bool_t val) {fSaveJetShapes = val; fInitialized = kFALSE;}
    void            SetSaveMCInformation(Bool_t val) {fSaveMCInformation = val; fInitialized = kFALSE;}
    void            SetSaveSecondaryVertices(Bool_t val) {fSaveSecondaryVertices = val; fInitialized = kFALSE;}
    void            SetSaveTriggerTracks(Bool_t val) {fSaveTriggerTracks = val;}
    void            SetTrackContainer(AliParticleContainer* cont) {fTrackContainer = cont;}
    void            SetJetContainer(AliJetContainer* cont) {fJetContainer = cont;}
    void            SetRandomGenerator(TRandom3* gen) {fRandomGenerator = gen;}

    void            InitializeTree();

    // ############ GETTERS
    Bool_t          GetSaveEventProperties() {return fSaveEventProperties;}
    Bool_t          GetSaveConstituents() {return fSaveConstituents;}
    Bool_t          GetSaveConstituentsIP() {return fSaveConstituentsIP;}
    Bool_t          GetSaveConstituentPID() {return fSaveConstituentPID;}
    Bool_t          GetSaveJetShapes() {return fSaveJetShapes;}
    Bool_t          GetSaveMCInformation() {return fSaveMCInformation;}
    Bool_t          GetSaveSecondaryVertices() {return fSaveSecondaryVertices;}
    Bool_t          GetSaveTriggerTracks() {return fSaveTriggerTracks;}

    std::vector<Float_t> GetExtractionPercentagePtBins() {return fExtractionPercentagePtBins;}
    std::vector<Float_t> GetExtractionPercentages() {return fExtractionPercentages;}
    std::vector<Int_t> GetExtractionJetTypes_HM() {return fExtractionJetTypes_HM;}
    std::vector<Int_t> GetExtractionJetTypes_PM() {return fExtractionJetTypes_PM;}

    TTree*          GetTreePointer() {return fJetTree;}

    // ######################################
    Bool_t          AddJetToTree(AliEmcalJet* jet, Float_t vertexX = 0, Float_t vertexY = 0, Float_t vertexZ = 0, Float_t centrality = 0, Long64_t eventID = 0, Float_t magField = 0,
      Int_t motherParton = 0, Int_t motherHadron = 0, Int_t partonInitialCollision = 0, Float_t matchDistance = 0, Float_t matchedPt = 0, Float_t matchedMass = 0, Float_t truePtFraction = 0, Float_t ptHard = 0, Float_t eventWeight = 0, Float_t impactParameter = 0,
      std::vector<Float_t>& trackPID_ITS = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& trackPID_TPC = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& trackPID_TOF = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& trackPID_TRD = DEFAULT_VECTOR_FLOAT, std::vector<Short_t>& trackPID_Reco = DEFAULT_VECTOR_SHORT, std::vector<Int_t>& trackPID_Truth = DEFAULT_VECTOR_INT,
      std::vector<Float_t>& triggerTrackPt = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& triggerTrackDeltaEta = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& triggerTrackDeltaPhi = DEFAULT_VECTOR_FLOAT,
      std::vector<Float_t>& trackIP_d0 = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& trackIP_z0 = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& trackIP_d0cov = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& trackIP_z0cov = DEFAULT_VECTOR_FLOAT,
      std::vector<Float_t>& secVtx_X = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_Y = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_Z = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_Mass = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_Lxy = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_SigmaLxy = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_Chi2 = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_Dispersion = DEFAULT_VECTOR_FLOAT)
    {
      return AddJetToTree(jet, vertexX, vertexY, vertexZ, centrality, eventID, magField, motherParton, motherHadron, partonInitialCollision, matchDistance, matchedPt, matchedMass, truePtFraction, ptHard, eventWeight, impactParameter,
        trackPID_ITS.data(), trackPID_TPC.data(), trackPID_TOF.data(), trackPID_TRD.data(), trackPID_Reco.data(), trackPID_Truth.data(),
        triggerTrackPt.size(), triggerTrackPt.data(), triggerTrackDeltaEta.data(), triggerTrackDeltaPhi.data(),
        trackIP_d0.data(), trackIP_z0.data(), trackIP_d0cov.data(), trackIP_z0cov.data(),
        secVtx_X.size(), secVtx_X.data(), secVtx_Y.data(), secVtx_Z.data(), secVtx_Mass.data(), secVtx_Lxy.data(), secVtx_SigmaLxy.data(), secVtx_Chi2.data(), secVtx_Dispersion.data());
    }
    Bool_t          AddJetToTree(AliEmcalJet* jet, Float_t vertexX, Float_t vertexY, Float_t vertexZ, Float_t centrality, Long64_t eventID, Float_t magField,
      Int_t motherParton, Int_t motherHadron, Int_t partonInitialCollision, Float_t matchDistance, Float_t matchedPt, Float_t matchedMass, Float_t truePtFraction, Float_t ptHard, Float_t eventWeight,  Float_t impactParameter,
      Float_t* trackPID_ITS, Float_t* trackPID_TPC, Float_t* trackPID_TOF, Float_t* trackPID_TRD, Short_t* trackPID_Reco, Int_t* trackPID_Truth,
      Int_t numTriggerTracks, Float_t* triggerTrackPt, Float_t* triggerTrackDeltaEta, Float_t* triggerTrackDeltaPhi,
      Float_t* trackIP_d0, Float_t* trackIP_z0, Float_t* trackIP_d0cov, Float_t* trackIP_z0cov,
      Int_t numSecVertices, Float_t* secVtx_X, Float_t* secVtx_Y, Float_t* secVtx_Z, Float_t* secVtx_Mass, Float_t* secVtx_Lxy, Float_t* secVtx_SigmaLxy, Float_t* secVtx_Chi2, Float_t* secVtx_Dispersion);
    // ######################################

  private:
    TTree*          fJetTree;                             //!<! tree structure
    Bool_t          fInitialized;                         ///< init state of tree

    // Save flags
    Bool_t          fSaveEventProperties;                 ///< save general event properties (bgrd. etc.)
    Bool_t          fSaveConstituents;                    ///< save arrays of constituent basic properties
    Bool_t          fSaveConstituentsIP;                  ///< save arrays of constituent impact parameters
    Bool_t          fSaveConstituentPID;                  ///< save arrays of constituent PID parameters
    Bool_t          fSaveJetShapes;                       ///< save jet shapes
    Bool_t          fSaveMCInformation;                   ///< save MC information
    Bool_t          fSaveSecondaryVertices;               ///< save reconstructed sec. vertex properties
    Bool_t          fSaveTriggerTracks;                   ///< save event trigger track

    AliParticleContainer* fTrackContainer;                //!<! track container, used to work on tracks
    AliJetContainer*      fJetContainer;                  //!<! jet container, used to work on jets
    TRandom3*             fRandomGenerator;               //!<! random generator

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
    Int_t           fBuffer_NumConstituents;              //!<! array buffer

    Float_t         fBuffer_Event_BackgroundDensityMass;  //!<! array buffer
    Float_t         fBuffer_Event_BackgroundDensity;      //!<! array buffer
    Float_t         fBuffer_Event_Vertex_X;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Y;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Z;               //!<! array buffer
    Float_t         fBuffer_Event_Centrality;             //!<! array buffer
    Long64_t        fBuffer_Event_ID;                     //!<! array buffer
    Float_t         fBuffer_Event_MagneticField;          //!<! array buffer
    Float_t         fBuffer_Event_PtHard;                 //!<! array buffer
    Float_t         fBuffer_Event_Weight;                 //!<! array buffer
    Float_t         fBuffer_Event_ImpactParameter;        //!<! array buffer

    Float_t*        fBuffer_Const_Pt;                     //!<! array buffer
    Float_t*        fBuffer_Const_Eta;                    //!<! array buffer
    Float_t*        fBuffer_Const_Phi;                    //!<! array buffer
    Float_t*        fBuffer_Const_Charge;                 //!<! array buffer
    Int_t*          fBuffer_Const_Label;                  //!<! array buffer
    Float_t*        fBuffer_Const_ProdVtx_X;              //!<! array buffer
    Float_t*        fBuffer_Const_ProdVtx_Y;              //!<! array buffer
    Float_t*        fBuffer_Const_ProdVtx_Z;              //!<! array buffer

    Float_t         fBuffer_Shape_Mass;                   //!<! array buffer

    Int_t           fBuffer_Jet_MC_MotherParton;          //!<! array buffer
    Int_t           fBuffer_Jet_MC_MotherHadron;          //!<! array buffer
    Int_t           fBuffer_Jet_MC_MotherIC;              //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchDistance;         //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchedPt;             //!<! array buffer
    Float_t         fBuffer_Jet_MC_MatchedMass;           //!<! array buffer
    Float_t         fBuffer_Jet_MC_TruePtFraction;        //!<! array buffer

    Int_t           fBuffer_NumTriggerTracks;
    Int_t           fBuffer_NumSecVertices;


    /// \cond CLASSIMP
    ClassDef(AliEmcalJetTree, 4) // Jet tree class
    /// \endcond
};

#endif
