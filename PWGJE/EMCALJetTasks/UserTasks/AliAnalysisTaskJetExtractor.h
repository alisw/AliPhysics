#ifndef ALIANALYSISTASKJETEXTRACTOR_H
#define ALIANALYSISTASKJETEXTRACTOR_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class THn;
class AliAODVertex;
class AliBasicJet;

//###############################################################################################################################################3

/**
 * \class AliAnalysisTaskJetExtractor
 * \brief Extract jets to tree
 *
 * This task extracts jets (including constituents) to a tree for further
 * analysis, e.g. to create a ML dataset
 * 
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date Dec 15, 2016
 */
// 
class AliAnalysisTaskJetExtractor : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskJetExtractor();
  AliAnalysisTaskJetExtractor(const char *name);
  virtual ~AliAnalysisTaskJetExtractor();
  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  // ################## SETTERS
  void                        SetExtractionCuts(Double_t percentage, Double_t minPt, Double_t maxPt, Int_t minCent = -1, Int_t maxCent = -1)
  {
    fExtractionPercentage = percentage;
    fExtractionCutMinPt = minPt;
    fExtractionCutMaxPt = maxPt;
    fExtractionCutMinCent = minCent;
    fExtractionCutMaxCent = maxCent;
  }

  void                        SetHadronMatchingRadius(Double_t val) { fHadronMatchingRadius = val; }
  void                        SetInitialCollisionMatchingRadius(Double_t val)     { fInitialCollisionMatchingRadius = val; }
  void                        SetTruthParticleArrayName(const char* val)     { fTruthParticleArrayName = val; }
  void                        SetSecondaryVertexUseThreeProng(Bool_t val)    { fSecondaryVertexUseThreeProng = val; }
  void                        SetSecondaryVertexMaxChi2(Int_t val)     { fSecondaryVertexMaxChi2 = val; }

  void                        SetExtractionCutListPIDHM(const char* val)
  { 
    fExtractionListPIDsHM.clear();
    TString tmpStr = val;
    TObjArray* tmpArr = tmpStr.Tokenize(",");
    for (Int_t i=0; i<tmpArr->GetEntries(); i++)
    {
      Int_t val = atoi((static_cast<TObjString*>(tmpArr->At(i)))->GetString().Data());
      fExtractionListPIDsHM.push_back(val);
    }
    fExtractionCutUseHM = kTRUE;
    fExtractionCutUseIC = kFALSE;
  }

  void                        SetExtractionCutListPIDIC(const char* val)
  { 
    fExtractionListPIDsIC.clear();
    TString tmpStr = val;
    TObjArray* tmpArr = tmpStr.Tokenize(",");
    for (Int_t i=0; i<tmpArr->GetEntries(); i++)
    {
      Int_t val = atoi((static_cast<TObjString*>(tmpArr->At(i)))->GetString().Data());
      fExtractionListPIDsIC.push_back(val);
    }
    fExtractionCutUseIC = kTRUE;
    fExtractionCutUseHM = kFALSE;
  }
  // ##################

 protected:
  Bool_t                      CreateControlHistograms();
  void                        ExecOnce();
  Bool_t                      Run();
  // ##################
  Bool_t                      IsJetSelected(AliEmcalJet* jet);
  void                        AddJetToTree(AliEmcalJet* jet);
  // ##################
  void                        FillEventControlHistograms();
  void                        FillJetControlHistograms(AliEmcalJet* jet);
  // ##################
  void                        CalculateEventProperties();
  void                        CalculateInitialCollisionJets();
  void                        GetLeadingJets(const char* opt, AliEmcalJet*& jetLeading, AliEmcalJet*& jetSubLeading);

  void                        CalculateJetProperties(AliEmcalJet* jet);
  void                        CalculateJetType(AliEmcalJet* jet, Int_t& typeIC, Int_t& typeHM);
  Int_t                       GetTrackPID(AliVParticle* part);
  Double_t                    GetTrackImpactParameter(const AliVVertex* vtx, AliAODTrack* track);
  void                        AddSecondaryVertices(const AliVVertex* vtx, AliEmcalJet* jet, AliBasicJet& basicJet);
  AliAODVertex*               GetSecondaryVertex(AliESDVertex* esdVtx, TObjArray* trkArray);


  // ################## BASIC EVENT VARIABLES
  AliJetContainer            *fJetsCont;                                //!<! Jets
  AliTrackContainer          *fTracksCont;                              //!<! Tracks
  TTree*                      fJetsTree;                                //!<! Jets that will be saved to a tree (optionally)
  void*                       fJetsTreeBuffer;                          //!<!  buffer for one jet (that will be saved to the tree)
  TRandom3*                   fRandom;                                  ///< random number generator

  // ################## CURRENT PROPERTIES
  Int_t                       fCurrentNJetsInEvents;                    ///< number of jets in event
  Int_t                       fCurrentJetTypeHM;                        ///< current jet type from hadron matching
  Int_t                       fCurrentJetTypeIC;                        ///< current jet type from initial collision (PYTHIA)
  AliEmcalJet*                fCurrentLeadingJet;                       //!<! leading jet (calculated event-by-event)
  AliEmcalJet*                fCurrentSubleadingJet;                    //!<! subleading jet (calculated event-by-event)
  AliEmcalJet*                fCurrentInitialParton1;                   //!<! jet that matches the initial parton 1 (PYTHIA)
  AliEmcalJet*                fCurrentInitialParton2;                   //!<! jet that matches the initial parton 2 (PYTHIA)
  Int_t                       fCurrentInitialParton1Type;               ///< type of initial parton 1
  Int_t                       fCurrentInitialParton2Type;               ///< type of initial parton 2
  Bool_t                      fFoundIC;                                 ///< status var showing that IC has been found

  // ################## CUTS
  Int_t                       fExtractionCutMinCent;                    ///< Extraction cut: minimum centrality
  Int_t                       fExtractionCutMaxCent;                    ///< Extraction cut: maximum centrality
  Bool_t                      fExtractionCutUseIC;                      ///< Extraction cut: Use IC PID
  Bool_t                      fExtractionCutUseHM;                      ///< Extraction cut: Use hadron matching PID
  Double_t                    fExtractionCutMinPt;                      ///< Extraction cut: minimum pT
  Double_t                    fExtractionCutMaxPt;                      ///< Extraction cut: minimum pT
  Double_t                    fExtractionPercentage;                    ///< Percentage of extracted jets
  std::vector<Int_t>          fExtractionListPIDsHM;                    ///< list of PIDs (hadron matching) that will be accepted
  std::vector<Int_t>          fExtractionListPIDsIC;                    ///< list of PIDs (initial collision) that will be accepted


  Double_t                    fHadronMatchingRadius;                    ///< Matching radius to search for beauty/charm hadrons around jet
  Double_t                    fInitialCollisionMatchingRadius;          ///< Matching radius to find a jet of the IC
  TString                     fTruthParticleArrayName;                  ///< Array name of MC particles in event (mcparticles)
  Bool_t                      fSecondaryVertexUseThreeProng;            ///< Reconstruct sec. vtx with 3 tracks (default: 2)
  Double_t                    fSecondaryVertexMaxChi2;                  ///< Max chi2 of secondary vertex (others will be discarded)

  // ######### HISTOGRAM FUNCTIONS
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

/**
 * \class AliBasicJetConstituent
 * \brief Simple class containing basic information for a constituent
 *
 * This class is used to save jet constituents with less overhead
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date Apr 21, 2016
 */
// 
class AliBasicJetConstituent
{
  public:
    AliBasicJetConstituent() : fEta(0), fPhi(0), fpT(0), fCharge(0), fPID(0), fVx(0), fVy(0), fVz(0), fImpactParameter(0) {}
    AliBasicJetConstituent(Float_t eta, Float_t phi, Float_t pt, Short_t charge, Short_t pid, Float_t vx, Float_t vy, Float_t vz, Float_t z)
    : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge), fPID(pid), fVx(vx), fVy(vy), fVz(vz), fImpactParameter(z)
    {
    }
    ~AliBasicJetConstituent();

    Float_t Pt()        { return fpT; }
    Float_t Phi()       { return fPhi; }
    Float_t Eta()       { return fEta; }
    Short_t Charge()    { return fCharge; }
    Short_t PID()       { return fPID; }

    Float_t Vx()        { return fVx; }
    Float_t Vy()        { return fVy; }
    Float_t Vz()        { return fVz; }

    Float_t ImpactParameter() { return fImpactParameter; }

  private:
    Float_t fEta;      ///< eta
    Float_t fPhi;      ///< phi
    Float_t fpT;       ///< pT
    Short_t fCharge;   ///< charge
    Short_t fPID;      ///< most probably PID code/PDG code

    Float_t fVx;       ///< production vertex X
    Float_t fVy;       ///< production vertex Y
    Float_t fVz;       ///< production vertex Z

    Float_t fImpactParameter; ///< impact parameter z
};

//###############################################################################################################################################3

/**
 * \class AliBasicJetSecondaryVertex
 * \brief Simple class containing basic information for a secondary
 * vertex
 * This class is used to save secondary vertices of a jet with less overhead
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date Dec 17, 2016
 */
// 
class AliBasicJetSecondaryVertex
{
  public:
    AliBasicJetSecondaryVertex() : fVx(0), fVy(0), fVz(0), fChi2(0) {}
    AliBasicJetSecondaryVertex(Float_t vx, Float_t vy, Float_t vz, Float_t chi2)
    : fVx(vx), fVy(vy), fVz(vz), fChi2(chi2)
    {
    }
    ~AliBasicJetSecondaryVertex();

    Float_t Vx()        { return fVx; }
    Float_t Vy()        { return fVy; }
    Float_t Vz()        { return fVz; }

    Float_t Chi2()      { return fChi2; }
  private:

    Float_t fVx;         ///< vertex X
    Float_t fVy;         ///< vertex Y
    Float_t fVz;         ///< vertex Z
    Float_t fChi2;       ///< Chi2/ndf for vertex
};

//###############################################################################################################################################3
/**
 * \class AliBasicJet
 * \brief Simple class containing basic information for a jet
 *
 * This class is used to save jets including constituents with minimal memory consumption.
 * Saved information focus on correlation analyses
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date Apr 21, 2016
 */
// 
class AliBasicJet
{
  public:
    AliBasicJet() : fEta(0), fPhi(0), fpT(0), fTruepT(0), fCharge(0), fRadius(0), fArea(0), fMotherInitialCollision(0), fMotherHadronMatching(0), fBackgroundDensity(0), fMagneticField(0), fVertexX(0), fVertexY(0), fVertexZ(0), fEventID(0), fCentrality(0), fConstituents() {}
    AliBasicJet(Float_t eta, Float_t phi, Float_t pt, Short_t charge, Float_t radius, Float_t area, Int_t partidIC, Int_t partidHM, Float_t bgrd, Float_t magfield, Float_t vtxX, Float_t vtxY, Float_t vtxZ, Long64_t id, Short_t cent)
    : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge), fRadius(radius), fArea(area), fMotherInitialCollision(partidIC), fMotherHadronMatching(partidHM), fBackgroundDensity(bgrd), fMagneticField(magfield), fVertexX(vtxX), fVertexY(vtxY), fVertexZ(vtxZ), fEventID(id), fCentrality(cent), fConstituents()
    {}
    ~AliBasicJet();

    // Basic jet properties
    Double_t                  Pt()       { return fpT; }
    Double_t                  TruePt()  { return fTruepT; }
    Double_t                  Phi()      { return fPhi; }
    Double_t                  Eta()      { return fEta; }
    Short_t                   Charge()   { return fCharge; }
    Double_t                  Radius() { return fRadius; }
    Double_t                  Area() { return fArea; }
    Int_t                     MotherInitialCollision() { return fMotherInitialCollision; }
    Int_t                     MotherHadronMatching() { return fMotherHadronMatching; }
    Double_t                  BackgroundDensity() { return fBackgroundDensity; }
    Double_t                  MagneticField() { return fMagneticField; }
    Double_t                  VertexX() { return fVertexX; }
    Double_t                  VertexY() { return fVertexY; }
    Double_t                  VertexZ() { return fVertexZ; }
    Long64_t                  EventID() { return fEventID; }
    Short_t                   Centrality() { return fCentrality; }
    Int_t                     GetNumbersOfConstituents() { return fConstituents.size(); }
    void                      SetTruePt(Double_t val)  {fTruepT = val;}

    // Basic constituent functions
    AliBasicJetConstituent*   GetJetConstituent(Int_t index) { return &fConstituents[index]; }
    void                      AddJetConstituent(Float_t eta, Float_t phi, Float_t pt, Short_t charge, Short_t pid=0, Float_t vx=0, Float_t vy=0, Float_t vz=0, Float_t z=0)
    {
      AliBasicJetConstituent c (eta, phi, pt, charge, pid, vx, vy, vz, z);
      AddJetConstituent(&c);
    }
    void                      AddJetConstituent(AliBasicJetConstituent* constituent) {fConstituents.push_back(*constituent); }

    // Basic secondary vertex functions
    AliBasicJetSecondaryVertex* GetSecondaryVertex(Int_t index) { return &fSecondaryVertices[index]; }
    void                      AddSecondaryVertex(Float_t vx, Float_t vy, Float_t vz, Float_t chi2=0)
    {
      AliBasicJetSecondaryVertex vtx (vx, vy, vz, chi2);
      AddSecondaryVertex(&vtx);
    }
    void                      AddSecondaryVertex(AliBasicJetSecondaryVertex* vtx) {fSecondaryVertices.push_back(*vtx);}

  private:
    Float_t   fEta;      ///< eta
    Float_t   fPhi;      ///< phi
    Float_t   fpT;       ///< pT
    Float_t   fTruepT;   ///< true pT (optional, e.g. from matching)
    Short_t   fCharge;   ///< charge
    Float_t   fRadius;   ///< jet radius
    Float_t   fArea;     ///< jet area
    Int_t     fMotherInitialCollision;  ///< PDG code of source particle (= initial collision quark/hadron)
    Int_t     fMotherHadronMatching;    ///< PDG code of source particle (according to matched hadrons around the jet)
    Float_t   fBackgroundDensity; ///< background
    Float_t   fMagneticField; ///< event magnetic field
    Float_t   fVertexX; ///< event vertex X
    Float_t   fVertexY; ///< event vertex Y
    Float_t   fVertexZ; ///< event vertex Z
    Long64_t  fEventID;  ///< Unique event id
    Short_t   fCentrality; ///< centrality

    std::vector<AliBasicJetSecondaryVertex> fSecondaryVertices; ///< vector of sec. vertices
    std::vector<AliBasicJetConstituent> fConstituents; ///< vector of constituents

};


#endif
