#ifndef ALIANALYSISTASKJETEXTRACTORHF_H
#define ALIANALYSISTASKJETEXTRACTORHF_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class THn;
class AliAODVertex;
class AliBasicJet;
class AliBasicJetConstituent;
class AliRDHFJetsCutsVertex;

//###############################################################################################################################################3

/**
 * \class AliAnalysisTaskJetExtractorHF
 * \brief Extract jets to tree
 *
 * This task extracts jets (including constituents) to a tree for further
 * analysis, e.g. to create a ML dataset. This task is written for HF jets
 * 
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date Jan 16, 2017
 */
// 
class AliAnalysisTaskJetExtractorHF : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskJetExtractorHF();
  AliAnalysisTaskJetExtractorHF(const char *name);
  virtual ~AliAnalysisTaskJetExtractorHF();
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

  void                        SetHadronMatchingRadius(Double_t val)               { fHadronMatchingRadius = val; }
  void                        SetTrueJetArrayName(const char* val)                { fTruthJetsArrayName = val; }
  void                        SetTrueRhoName(const char* val)                     { fTruthJetsRhoName = val; }
  void                        SetTruthParticleArrayName(const char* val)          { fTruthParticleArrayName = val; }
  void                        SetSecondaryVertexMaxChi2(Double_t val   )          { fSecondaryVertexMaxChi2 = val; }
  void                        SetSecondaryVertexMaxDispersion(Double_t val)       { fSecondaryVertexMaxDispersion = val; }
  void                        SetAddPIDSignal(Bool_t val)                         { fAddPIDSignal = val; }
  void                        SetCalculateSecondaryVertices(Bool_t val)           { fCalculateSecondaryVertices = val; }
  void                        SetVertexerCuts(AliRDHFJetsCutsVertex* val)         { fVertexerCuts = val; }
  void                        SetSetEmcalJetFlavour(Bool_t val)                   { fSetEmcalJetFlavour = val; }

  // This setter configures the extraction of jets below the extraction min pt cut
  void                        SetUnderflowBins(Int_t nBins, Double_t minPt, Double_t percentage)
  {
    fUnderflowNumBins     =  nBins;
    fUnderflowCutOff      =  minPt;
    fUnderflowPercentage  =  percentage;
  }

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
    fExtractionCutUsePM = kFALSE;
  }

  void                        SetExtractionCutListPIDPM(const char* val)
  { 
    fExtractionListPIDsPM.clear();
    TString tmpStr = val;
    TObjArray* tmpArr = tmpStr.Tokenize(",");
    for (Int_t i=0; i<tmpArr->GetEntries(); i++)
    {
      Int_t val = atoi((static_cast<TObjString*>(tmpArr->At(i)))->GetString().Data());
      fExtractionListPIDsPM.push_back(val);
    }
    fExtractionCutUsePM = kTRUE;
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
  void                        FillTrackControlHistograms(AliVTrack* track);
  void                        FillEventControlHistograms();
  void                        FillJetControlHistograms(AliEmcalJet* jet);
  // ##################
  void                        CalculateEventProperties();
  void                        GetLeadingJets(const char* opt, AliEmcalJet*& jetLeading, AliEmcalJet*& jetSubLeading);

  void                        FindMatchingJet(AliEmcalJet* jet);
  void                        CalculateJetProperties(AliEmcalJet* jet);
  void                        CalculateJetType(AliEmcalJet* jet, Int_t& typeHM, Int_t& typePM);
  Bool_t                      IsStrangeJet(AliEmcalJet* jet);
  void                        GetTrackImpactParameters(const AliVVertex* vtx, AliAODTrack* track, Double_t& d0, Double_t& d0cov, Double_t& z0, Double_t& z0cov);
  void                        AddSecondaryVertices(const AliVVertex* primVtx, const AliEmcalJet* jet, AliBasicJet& basicJet);
  void                        AddPIDInformation(AliVParticle* particle, AliBasicJetConstituent& constituent);
  Bool_t                      IsTrackInCone(AliVParticle* track, Double_t eta, Double_t phi, Double_t radius);


  // ################## BASIC EVENT VARIABLES
  TClonesArray*               fTruthParticleArray;                      //!<! Array of MC particles in event (mcparticles)
  AliJetContainer            *fJetsCont;                                //!<! Jets
  AliTrackContainer          *fTracksCont;                              //!<! Tracks
  TTree*                      fJetsTree;                                //!<! Jets that will be saved to a tree (optionally)
  void*                       fJetsTreeBuffer;                          //!<!  buffer for one jet (that will be saved to the tree)
  TRandom3*                   fRandom;                                  ///< random number generator
  AliHFJetsTaggingVertex*     fVtxTagger;                               //!<! class for sec. vertexing

  // ################## CURRENT PROPERTIES
  Int_t                       fCurrentNJetsInEvents;                    ///< number of tagged jets in event
  Int_t                       fCurrentNJetsInEventsTotal;               ///< number of inclusive jets in event
  Int_t                       fCurrentJetTypeHM;                        ///< current jet type from hadron matching
  Int_t                       fCurrentJetTypePM;                        ///< current jet type from parton matching
  AliEmcalJet*                fCurrentLeadingJet;                       //!<! leading jet (calculated event-by-event)
  AliEmcalJet*                fCurrentSubleadingJet;                    //!<! subleading jet (calculated event-by-event)
  Int_t                       fCurrentJetHFHadron;                      ///< PDG code for HF hadron in hadron-matched HF-jets
  Double_t                    fCurrentTrueJetPt;                        ///< truth jet pt buffer
  Int_t                       fUnderflowBinContents[100];               //!<! contents in "underflow" bins

  // ################## CUTS
  Int_t                       fExtractionCutMinCent;                    ///< Extraction cut: minimum centrality
  Int_t                       fExtractionCutMaxCent;                    ///< Extraction cut: maximum centrality
  Bool_t                      fExtractionCutUsePM;                      ///< Extraction cut: Use parton matching PID
  Bool_t                      fExtractionCutUseHM;                      ///< Extraction cut: Use hadron matching PID
  Double_t                    fExtractionCutMinPt;                      ///< Extraction cut: minimum pT
  Double_t                    fExtractionCutMaxPt;                      ///< Extraction cut: minimum pT
  Double_t                    fExtractionPercentage;                    ///< Percentage of extracted jets
  std::vector<Int_t>          fExtractionListPIDsHM;                    ///< list of PIDs (hadron matching) that will be accepted
  std::vector<Int_t>          fExtractionListPIDsPM;                    ///< list of PIDs (parton matching) that will be accepted
  Bool_t                      fSetEmcalJetFlavour;                      ///< if set, the flavour property of the AliEmcalJets will be set
  Int_t                       fUnderflowNumBins;                        ///< number of "underflow" low-pT bins
  Double_t                    fUnderflowCutOff;                         ///< "underflow" low-pT cut-off
  Double_t                    fUnderflowPercentage;                     ///< Percentage accepted in "underflow" low-pT bins

  Double_t                    fHadronMatchingRadius;                    ///< Matching radius to search for beauty/charm hadrons around jet
  TString                     fTruthJetsArrayName;                      ///< Array name for particle-level jets
  TString                     fTruthJetsRhoName;                        ///< Array name for particle-level rho
  TString                     fTruthParticleArrayName;                  ///< Array name of MC particles in event (mcparticles)
  Double_t                    fSecondaryVertexMaxChi2;                  ///< Max chi2 of secondary vertex (others will be discarded)
  Double_t                    fSecondaryVertexMaxDispersion;            ///< Max dispersion of secondary vertex (others will be discarded)
  Bool_t                      fAddPIDSignal;                            ///< Add pid signal to each jet constituent
  Bool_t                      fCalculateSecondaryVertices;              ///< Calculate the secondary vertices (instead of loading)
  AliRDHFJetsCutsVertex*      fVertexerCuts;                            ///< Cuts used for the vertexer (given in add task macro)


  // ######### HISTOGRAM FUNCTIONS
  void                        FillHistogram(const char * key, Double_t x);
  void                        FillHistogram(const char * key, Double_t x, Double_t y);
  void                        FillHistogram(const char * key, Double_t x, Double_t y, Double_t add);
  void                        FillHistogram3D(const char * key, Double_t x, Double_t y, Double_t z, Double_t add = 0);
  template <class T> T*       AddHistogram1D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis");
  template <class T> T*       AddHistogram2D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0,  const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");
  template <class T> T*       AddHistogram3D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0, Int_t zBins = 100, Double_t zMin = 0.0, Double_t zMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");


 private:
  AliAnalysisTaskJetExtractorHF(const AliAnalysisTaskJetExtractorHF&);            // not implemented
  AliAnalysisTaskJetExtractorHF &operator=(const AliAnalysisTaskJetExtractorHF&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetExtractorHF, 8) // Jet extraction task
  /// \endcond
};


//###############################################################################################################################################3

/**
 * \class AliBasicPID
 * \brief Simple class containing available (raw) information from PID detectors
 *
 * This class is used to save PID information for each constituents in a tree
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date Dec 20, 2016
 */
// 
class AliBasicPID
{
  public:
    AliBasicPID() : fSignalITS(0), fSignalTPC(0), fSignalTOF(0), fSignalTRD(0), fTruthPID(9), fRecoPID(9)  {}
    AliBasicPID(Float_t its, Float_t tpc, Float_t tof, Float_t trd, Short_t truthPID, Short_t recoPID)
    : fSignalITS(its), fSignalTPC(tpc), fSignalTOF(tof), fSignalTRD(trd), fTruthPID(truthPID), fRecoPID(recoPID)
    {
    }
    ~AliBasicPID();

    Float_t SignalITS()        { return fSignalITS; }
    Float_t SignalTPC()        { return fSignalTPC; }
    Float_t SignalTOF()        { return fSignalTOF; }
    Float_t SignalTRD()        { return fSignalTRD; }
    Short_t TruthPID()         { return fTruthPID; }
    Short_t RecoPID()          { return fRecoPID; }

  private:
    Float_t fSignalITS;      ///< ITS PID signal
    Float_t fSignalTPC;      ///< TPC PID signal
    Float_t fSignalTOF;      ///< TOF PID signal
    Float_t fSignalTRD;      ///< TRD PID signal
    Short_t fTruthPID;       ///< truth PID
    Short_t fRecoPID;        ///< reco PID
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
    AliBasicJetConstituent() : fEta(0), fPhi(0), fpT(0), fCharge(0), fVx(0), fVy(0), fVz(0), fImpactParameterD(0), fImpactParameterZ(0), fPID(0) {}
    AliBasicJetConstituent(Float_t eta, Float_t phi, Float_t pt, Short_t charge, Float_t vx, Float_t vy, Float_t vz, Float_t d0, Float_t d0cov, Float_t z0, Float_t z0cov)
    : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge), fVx(vx), fVy(vy), fVz(vz), fImpactParameterD(d0), fImpactParameterZ(z0), fCovImpactParameterD(d0cov), fCovImpactParameterZ(z0cov), fPID(0)
    {
    }
    ~AliBasicJetConstituent();

    Float_t Pt()        { return fpT; }
    Float_t Phi()       { return fPhi; }
    Float_t Eta()       { return fEta; }
    Short_t Charge()    { return fCharge; }

    Float_t Vx()        { return fVx; }
    Float_t Vy()        { return fVy; }
    Float_t Vz()        { return fVz; }

    Float_t SigImpactParameterD() { return (fCovImpactParameterD > 0) ? (fImpactParameterD/TMath::Sqrt(fCovImpactParameterD)) : 0; }
    Float_t SigImpactParameterZ() { return (fCovImpactParameterZ > 0) ? (fImpactParameterZ/TMath::Sqrt(fCovImpactParameterZ)) : 0; }

    Float_t ImpactParameterD() { return fImpactParameterD; }
    Float_t ImpactParameterZ() { return fImpactParameterZ; }
    Float_t CovImpactParameterD() { return fCovImpactParameterD; }
    Float_t CovImpactParameterZ() { return fCovImpactParameterZ; }

    AliBasicPID* PID()  { return &fPID.at(0); }

    void SetPIDSignal(Float_t its, Float_t tpc, Float_t tof, Float_t trd, Short_t truthPID, Short_t recoPID)
    { 
      if(fPID.size()) fPID.clear();
      AliBasicPID* newPID = new AliBasicPID(its,tpc,tof,trd,truthPID,recoPID);
      fPID.push_back(*newPID);
    }

  private:
    Float_t fEta;      ///< eta
    Float_t fPhi;      ///< phi
    Float_t fpT;       ///< pT
    Short_t fCharge;   ///< charge

    Float_t fVx;       ///< production vertex X
    Float_t fVy;       ///< production vertex Y
    Float_t fVz;       ///< production vertex Z

    Float_t fImpactParameterD; ///< impact parameter d (transversal IP)
    Float_t fImpactParameterZ; ///< impact parameter z (longitudional IP)
    Float_t fCovImpactParameterD; ///< cov D
    Float_t fCovImpactParameterZ; ///< cov Z

    std::vector<AliBasicPID> fPID; ///< PID 
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
    AliBasicJetSecondaryVertex() : fVx(0), fVy(0), fVz(0), fMass(0), fLxy(0), fSigmaLxy(0), fChi2(0), fDispersion(0) {}
    AliBasicJetSecondaryVertex(Float_t vx, Float_t vy, Float_t vz, Float_t chi2, Float_t dispersion, Float_t mass, Float_t lxy, Float_t sigmalxy)
    : fVx(vx), fVy(vy), fVz(vz), fMass(mass), fLxy(lxy), fSigmaLxy(sigmalxy), fChi2(chi2), fDispersion(dispersion)
    {
    }
    ~AliBasicJetSecondaryVertex();

    Float_t Vx()        { return fVx; }
    Float_t Vy()        { return fVy; }
    Float_t Vz()        { return fVz; }

    Float_t Mass()      { return fMass; }
    Float_t Lxy()       { return fLxy; }
    Float_t SigmaLxy()  { return fSigmaLxy; }

    Float_t Chi2()      { return fChi2; }
    Float_t Dispersion(){ return fDispersion; }
  private:

    Float_t fVx;         ///< vertex X
    Float_t fVy;         ///< vertex Y
    Float_t fVz;         ///< vertex Z
    Float_t fMass;       ///< Invariant mass of vertex daughters
    Float_t fLxy;        ///< Signed decay length in XY
    Float_t fSigmaLxy;   ///< Sigma of decay length in XY
    Float_t fChi2;       ///< Chi2/ndf for vertex
    Float_t fDispersion; ///< Dispersion of vertex
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
    AliBasicJet() : fEta(0), fPhi(0), fpT(0), fTruepT(0), fCharge(0), fMass(0), fRadius(0), fArea(0), fMotherPartonMatching(0), fMotherHadronMatching(0), fBackgroundDensity(0), fMagneticField(0), fVertexX(0), fVertexY(0), fVertexZ(0), fEventID(0), fCentrality(0), fEventPtHard(0), fConstituents() {}
    AliBasicJet(Float_t eta, Float_t phi, Float_t pt, Short_t charge, Float_t radius, Float_t area, Int_t partidPM, Int_t partidHM, Float_t bgrd, Float_t magfield, Float_t vtxX, Float_t vtxY, Float_t vtxZ, Long64_t id, Short_t cent, Float_t mass, Float_t ptHard)
    : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge), fMass(mass), fRadius(radius), fArea(area), fMotherPartonMatching(partidPM), fMotherHadronMatching(partidHM), fBackgroundDensity(bgrd), fMagneticField(magfield), fVertexX(vtxX), fVertexY(vtxY), fVertexZ(vtxZ), fEventID(id), fCentrality(cent), fEventPtHard(ptHard), fConstituents()
    {}
    ~AliBasicJet();

    // Basic jet properties
    Double_t                  Pt()       { return fpT; }
    Double_t                  TruePt()  { return fTruepT; }
    Double_t                  Phi()      { return fPhi; }
    Double_t                  Eta()      { return fEta; }
    Short_t                   Charge()   { return fCharge; }
    Double_t                  Mass()     { return fMass; }
    Double_t                  Radius() { return fRadius; }
    Double_t                  Area() { return fArea; }
    Int_t                     MotherPartonMatching() { return fMotherPartonMatching; }
    Int_t                     MotherHadronMatching() { return fMotherHadronMatching; }
    Double_t                  BackgroundDensity() { return fBackgroundDensity; }
    Double_t                  MagneticField() { return fMagneticField; }
    Double_t                  VertexX() { return fVertexX; }
    Double_t                  VertexY() { return fVertexY; }
    Double_t                  VertexZ() { return fVertexZ; }
    Long64_t                  EventID() { return fEventID; }
    Short_t                   Centrality() { return fCentrality; }
    Double_t                  EventPtHard() { return fEventPtHard; }
    Int_t                     GetNumbersOfConstituents() { return (Int_t)fConstituents.size(); }
    Int_t                     GetNumbersOfSecVertices() { return (Int_t)fSecondaryVertices.size(); }
    void                      SetTruePt(Double_t val)  {fTruepT = val;}

    // Basic constituent functions
    AliBasicJetConstituent*   GetJetConstituent(Int_t index) { return &fConstituents[index]; }
    void                      AddJetConstituent(Float_t eta, Float_t phi, Float_t pt, Short_t charge, Float_t vx=0, Float_t vy=0, Float_t vz=0, Float_t d0=0, Float_t d0cov=0, Float_t z0=0, Float_t z0cov=0)
    {
      AliBasicJetConstituent c (eta, phi, pt, charge, vx, vy, vz, d0, d0cov, z0, z0cov);
      AddJetConstituent(&c);
    }
    void                      AddJetConstituent(AliBasicJetConstituent* constituent) {fConstituents.push_back(*constituent); }

    // Basic secondary vertex functions
    AliBasicJetSecondaryVertex* GetSecondaryVertex(Int_t index) { return &fSecondaryVertices[index]; }
    void                      AddSecondaryVertex(Float_t vx, Float_t vy, Float_t vz, Float_t chi2=0, Float_t dispersion=0, Float_t mass=0, Float_t lxy=0, Float_t sigmalxy=0)
    {
      AliBasicJetSecondaryVertex vtx (vx, vy, vz, chi2, dispersion, mass, lxy, sigmalxy);
      AddSecondaryVertex(&vtx);
    }
    void                      AddSecondaryVertex(AliBasicJetSecondaryVertex* vtx) {fSecondaryVertices.push_back(*vtx);}

  private:
    Float_t   fEta;      ///< eta
    Float_t   fPhi;      ///< phi
    Float_t   fpT;       ///< pT
    Float_t   fTruepT;   ///< true pT (optional, e.g. from matching)
    Short_t   fCharge;   ///< charge
    Float_t   fMass;     ///< jet mass
    Float_t   fRadius;   ///< jet radius
    Float_t   fArea;     ///< jet area
    Int_t     fMotherPartonMatching; ///< PDG code of source particle (according to matched partons around the jet)
    Int_t     fMotherHadronMatching; ///< PDG code of source particle (according to matched hadrons around the jet)
    Float_t   fBackgroundDensity; ///< background
    Float_t   fMagneticField; ///< event magnetic field
    Float_t   fVertexX; ///< event vertex X
    Float_t   fVertexY; ///< event vertex Y
    Float_t   fVertexZ; ///< event vertex Z
    Long64_t  fEventID;  ///< Unique event id
    Short_t   fCentrality; ///< centrality
    Float_t   fEventPtHard; ///< pt hard of the corresponding event

    std::vector<AliBasicJetSecondaryVertex> fSecondaryVertices; ///< vector of sec. vertices
    std::vector<AliBasicJetConstituent> fConstituents; ///< vector of constituents

};


#endif
