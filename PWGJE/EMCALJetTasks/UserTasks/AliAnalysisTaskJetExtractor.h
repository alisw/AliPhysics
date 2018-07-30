#ifndef ALIANALYSISTASKJETEXTRACTOR_H
#define ALIANALYSISTASKJETEXTRACTOR_H

/* Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//###############################################################################################################################################3
class AliEmcalJetTree;

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
  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  AliEmcalJetTree*            GetJetTree() {return fJetTree;}

 protected:
  Bool_t                      CreateControlHistograms();
  void                        ExecOnce();
  Bool_t                      Run();
  void                        FillTrackControlHistograms(AliVTrack* track);
  void                        FillEventControlHistograms();
  void                        FillJetControlHistograms(AliEmcalJet* jet);
  void                        AddPIDInformation(AliVParticle* particle, Float_t& sigITS, Float_t& sigTPC, Float_t& sigTOF, Float_t& sigTRD, Short_t& recoPID);
  Bool_t                      IsTrackInCone(AliVParticle* track, Double_t eta, Double_t phi, Double_t radius);

  AliEmcalJetTree*            fJetTree;

  // ################## BASIC EVENT VARIABLES
  AliJetContainer            *fJetsCont;                                //!<! Jets
  AliTrackContainer          *fTracksCont;                              //!<! Tracks


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

    void            SetSaveEventProperties(Bool_t val) {fSaveEventProperties = val; fInitialized = kFALSE;}
    void            SetSaveConstituents(Bool_t val) {fSaveConstituents = val; fInitialized = kFALSE;}
    void            SetSaveConstituentsIP(Bool_t val) {fSaveConstituentsIP = val; fInitialized = kFALSE;}
    void            SetSaveConstituentPID(Bool_t val) {fSaveConstituentPID = val; fInitialized = kFALSE;}
    void            SetSaveJetShapes(Bool_t val) {fSaveJetShapes = val; fInitialized = kFALSE;}
    void            SetSaveMCInformation(Bool_t val) {fSaveMCInformation = val; fInitialized = kFALSE;}
    void            SetSaveSecondaryVertices(Bool_t val) {fSaveSecondaryVertices = val; fInitialized = kFALSE;}

    void            InitializeTree();
    TTree*          GetTreePointer() {return fJetTree;}

    // ######################################
    Bool_t          AddJetToTree(AliEmcalJet* jet, Float_t bgrdDensity = 0, Float_t vertexX = 0, Float_t vertexY = 0, Float_t vertexZ = 0, Float_t centrality = 0, Long64_t eventID = 0, Float_t magField = 0,
      AliTrackContainer* trackCont = 0, Int_t motherParton = 0, Int_t motherHadron = 0, Float_t truePt = 0, Float_t ptHard = 0,
      std::vector<Float_t>& trackPID_ITS = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& trackPID_TPC = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& trackPID_TOF = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& trackPID_TRD = DEFAULT_VECTOR_FLOAT, std::vector<Short_t>& trackPID_Reco = DEFAULT_VECTOR_SHORT, std::vector<Short_t>& trackPID_Truth = DEFAULT_VECTOR_SHORT,
      std::vector<Float_t>& trackIP_d0 = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& trackIP_z0 = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& trackIP_d0cov = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& trackIP_z0cov = DEFAULT_VECTOR_FLOAT,
      std::vector<Float_t>& secVtx_X = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_Y = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_Z = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_Mass = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_Lxy = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_SigmaLxy = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_Chi2 = DEFAULT_VECTOR_FLOAT, std::vector<Float_t>& secVtx_Dispersion = DEFAULT_VECTOR_FLOAT)
    {
      return AddJetToTree(jet, bgrdDensity, vertexX, vertexY, vertexZ, centrality, eventID, magField, trackCont, motherParton, motherHadron, truePt, ptHard,
        trackPID_ITS.data(), trackPID_TPC.data(), trackPID_TOF.data(), trackPID_TRD.data(), trackPID_Reco.data(), trackPID_Truth.data(),
        trackIP_d0.data(), trackIP_z0.data(), trackIP_d0cov.data(), trackIP_z0cov.data(),
        secVtx_X.size(), secVtx_X.data(), secVtx_Y.data(), secVtx_Z.data(), secVtx_Mass.data(), secVtx_Lxy.data(), secVtx_SigmaLxy.data(), secVtx_Chi2.data(), secVtx_Dispersion.data());
    }
    Bool_t          AddJetToTree(AliEmcalJet* jet, Float_t bgrdDensity, Float_t vertexX, Float_t vertexY, Float_t vertexZ, Float_t centrality, Long64_t eventID, Float_t magField,
      AliTrackContainer* trackCont, Int_t motherParton, Int_t motherHadron, Float_t truePt, Float_t ptHard,
      Float_t* trackPID_ITS, Float_t* trackPID_TPC, Float_t* trackPID_TOF, Float_t* trackPID_TRD, Short_t* trackPID_Reco, Short_t* trackPID_Truth,
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

    // Option flags
    std::vector<Float_t> fExtractionPercentages;          ///< Percentages which will be extracted for a given pT bin
    std::vector<Float_t> fExtractionPercentagePtBins;     ///< pT-bins associated with fExtractionPercentages

    // Buffers that will be added to the tree
    Float_t         fBuffer_JetPt;                        //!<! array buffer
    Float_t         fBuffer_JetEta;                       //!<! array buffer
    Float_t         fBuffer_JetPhi;                       //!<! array buffer
    Float_t         fBuffer_JetArea;                      //!<! array buffer
    Int_t           fBuffer_NumConstituents;              //!<! array buffer

    Float_t         fBuffer_Event_BackgroundDensity;      //!<! array buffer
    Float_t         fBuffer_Event_Vertex_X;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Y;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Z;               //!<! array buffer
    Float_t         fBuffer_Event_Centrality;             //!<! array buffer
    Long64_t        fBuffer_Event_ID;                     //!<! array buffer
    Float_t         fBuffer_Event_MagneticField;          //!<! array buffer
    Float_t         fBuffer_Event_PtHard;                 //!<! array buffer

    Float_t*        fBuffer_Const_Pt;                     //!<! array buffer
    Float_t*        fBuffer_Const_Eta;                    //!<! array buffer
    Float_t*        fBuffer_Const_Phi;                    //!<! array buffer
    Float_t*        fBuffer_Const_Charge;                 //!<! array buffer
    Float_t*        fBuffer_Const_ProdVtx_X;              //!<! array buffer
    Float_t*        fBuffer_Const_ProdVtx_Y;              //!<! array buffer
    Float_t*        fBuffer_Const_ProdVtx_Z;              //!<! array buffer

    Float_t         fBuffer_Shape_Mass;                   //!<! array buffer

    Int_t           fBuffer_Jet_MC_MotherParton;          //!<! array buffer
    Int_t           fBuffer_Jet_MC_MotherHadron;          //!<! array buffer
    Float_t         fBuffer_Jet_MC_TruePt;                //!<! array buffer

    Int_t           fBuffer_NumSecVertices;

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetTree, 1) // Jet tree class
    /// \endcond
};

#endif
