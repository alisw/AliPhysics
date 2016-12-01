#ifndef ALIANALYSISTASKCHARGEDJETSHADRONCF_H
#define ALIANALYSISTASKCHARGEDJETSHADRONCF_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class THn;

//###############################################################################################################################################3
/**
 * \class AliChargedJetsHadronCFCuts
 * \brief Container class of cuts for AliAnalysisTaskChargedJetsHadronCF
 *
 *
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date Nov 20, 2016
 */
// 
class AliChargedJetsHadronCFCuts
{
 public:
  AliChargedJetsHadronCFCuts() {;}
  AliChargedJetsHadronCFCuts(const char* cutName, const char* outputName, Double_t minPt, Double_t maxPt, Double_t minCent, Double_t maxCent, Double_t lowerOffset, Double_t offset2D, Double_t slope2D, Double_t jetVetoPt)
  {
    fCutName = cutName;
    fOutputName = outputName;
    fPtMin = minPt;
    fPtMax = maxPt;
    fCentMin = minCent;
    fCentMax = maxCent;

    fCutLowerPercentageOffset = lowerOffset;
    fCut2DOffset = offset2D;
    fCut2DSlope = slope2D;

    fJetVetoPt = jetVetoPt;

    fAcceptedJets = 0;
    fArrayIndex = -1;
  }
  Bool_t IsCutFulfilled(Double_t pt, Double_t mcPt, Double_t cent, Double_t ptRatio, Double_t vetoPt)
  {
    // Simple kinematic cuts
    if(pt < fPtMin || pt >= fPtMax)
      return kFALSE;
    if(cent < fCentMin || cent >= fCentMax)
      return kFALSE;

    // jet veto
    if(vetoPt >= fJetVetoPt)
      return kFALSE;

    // Lower MC percentage cut
    if(ptRatio < fCutLowerPercentageOffset)
      return kFALSE;

    // Value of the 2D cut line. We cut everything above that
    Double_t cutLineValue = fCut2DSlope*mcPt + fCut2DOffset;
    if(ptRatio > cutLineValue)
      return kFALSE;

    return kTRUE;
  }

  ~AliChargedJetsHadronCFCuts();
  TString     fCutName;                              ///< name of this cut
  TString     fOutputName;                           ///< array name that is used to output objects associated with this cuts
  Double_t    fPtMin;                                ///< valid for jets above this pT
  Double_t    fPtMax;                                ///< valid for jets below this pT
  Double_t    fCentMin;                              ///< valid for centralities above
  Double_t    fCentMax;                              ///< valid for centralities below
  Double_t    fJetVetoPt;                            ///< jet veto pt

  Double_t    fCutLowerPercentageOffset;             ///< cut value
  Double_t    fCut2DOffset;                          ///< cut value
  Double_t    fCut2DSlope;                           ///< cut value


  Int_t       fAcceptedJets;                         ///< temporary var that holds how many jets passed
  Int_t       fArrayIndex;                           ///< array index that holds the output array index
};

//###############################################################################################################################################3

/**
 * \class AliAnalysisTaskChargedJetsHadronCF
 * \brief Support task for (charged) jet-hadron correlations
 *
 * This task selects the jets and optionally the particles which are saved as a TClonesArray to the event
 * and which are read by a correlation task. The tasks does a special event selection (e.g. dijet events)
 * and jet selection (e.g. fakefactor rejection). In addition, it fills auxiliary histograms
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date Apr 21, 2016
 */
// 
class AliAnalysisTaskChargedJetsHadronCF : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskChargedJetsHadronCF();
  AliAnalysisTaskChargedJetsHadronCF(const char *name);
  virtual ~AliAnalysisTaskChargedJetsHadronCF();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  // ######### SETTERS/GETTERS
  void                        ActivateJetEmbedding(const char* embArray, const char* embTrackArray, Double_t maxDistance, Int_t numMatchedJets)
  {
    fJetEmbeddingMaxDistance      = maxDistance;
    fJetEmbeddingArrayName        = embArray;
    fJetEmbeddingTrackArrayName   = embTrackArray;
    fJetEmbeddingNumMatchedJets   = numMatchedJets;
  }
  void                        AddJetEmbeddingCut(const char* cutName, const char* outputName, Double_t minPt, Double_t maxPt, Double_t minCent, Double_t maxCent, Double_t lowerOffset, Double_t offset2D, Double_t slope2D, Double_t jetVetoPt)
  {
    AliChargedJetsHadronCFCuts tmpCut(cutName, outputName, minPt, maxPt, minCent, maxCent, lowerOffset, offset2D, slope2D, jetVetoPt);
    fJetEmbeddingCuts.push_back(tmpCut);
  }

  void                        SetJetVetoArrayName(const char* name)   { fJetVetoArrayName = name; }
  void                        SetJetVetoJetByJet(Bool_t val)          { fJetVetoJetByJet = val; }
  void                        SetUsePerTrackMCPercentage(Bool_t val)   { fJetEmbeddingUsePerTrackMCPercentage = val; }
  void                        SetUseBgrdForMCPercentage(Bool_t val)   { fJetEmbeddingUseBgrdForMCPercentage = val; }
  void                        SetCreateEmbeddingPtPlotPerCut(Bool_t val)   { fJetEmbeddingCreatePtPlotPerCut = val; }
  void                        SetConstPtFilterBit(Int_t val)   { fConstPtFilterBit = val; }
  void                        SetNumberOfCentralityBins(Int_t val)   { fNumberOfCentralityBins = val; }
  void                        SetJetParticleArrayName(const char* name)   { fJetParticleArrayName = name; }
  void                        SetTrackParticleArrayName(const char* name) { fTrackParticleArrayName = name; }

  void                        SetJetOutputMode(Int_t mode) {fJetOutputMode = mode;}
  void                        SetPythiaExtractionMode(Int_t mode) {fPythiaExtractionMode = mode;}

  void                        ActivateJetExtraction(Double_t percentage, Double_t minPt, Double_t maxPt) {fExtractionPercentage = percentage; fExtractionMinPt = minPt; fExtractionMaxPt = maxPt;}
  void                        ActivateEventExtraction(Double_t percentage, Double_t minJetPt, Double_t maxJetPt) {fEventExtractionPercentage = percentage; fEventExtractionMinJetPt = minJetPt; fEventExtractionMaxJetPt = maxJetPt;}

 protected:
  void                        ExecOnce();
  Bool_t                      Run();

  // ######### META FUNCTIONS
  void                        BinLogAxis(const THn *h, Int_t axisNumber);
  void                        AddJetToTree(AliEmcalJet* jet);
  void                        AddEventToTree();
  void                        AddJetToOutputArray(AliEmcalJet* jet, Int_t arrayIndex, Int_t& jetsAlreadyInArray);
  void                        AddTrackToOutputArray(AliVTrack* track);
  void                        FillHistogramsTracks(AliVTrack* track);
  void                        FillHistogramsJets(AliEmcalJet* jet, const char* cutName);
  Bool_t                      IsJetSelected(AliEmcalJet* jet);

  AliJetContainer            *fJetsCont;                                //!<! Jets
  AliTrackContainer          *fTracksCont;                              //!<! Tracks
  TTree*                      fJetsTree;                                //!<! Jets that will be saved to a tree (optionally)
  void*                       fJetsTreeBuffer;                          //!<!  buffer for one jet (that will be saved to the tree)
  Double_t                    fExtractionPercentage;                    ///< percentage that is recorded
  Double_t                    fExtractionMinPt;                         ///< minimum pt of recorded jets
  Double_t                    fExtractionMaxPt;                         ///< maximum pt of recorded jets
  Double_t                    fEventExtractionPercentage;               ///< percentage of events that is recorded
  Double_t                    fEventExtractionMinJetPt;                 ///< minimum jet pt of recorded events
  Double_t                    fEventExtractionMaxJetPt;                 ///< maximum jet pt of recorded events
  
  Int_t                       fConstPtFilterBit;                        ///< For const pt plot, filter bit
  Int_t                       fNumberOfCentralityBins;                  ///< Number of centrality bins
  std::vector<TClonesArray*>  fJetsOutput;                              //!<! vector of arrays of basic correlation particles attached to the event (jets)
  TClonesArray               *fTracksOutput;                            //!<! Array of basic correlation particles attached to the event (tracks)
  TString                     fJetParticleArrayName;                    ///< Name of fJetsOutput array (if one uses only one)
  TString                     fTrackParticleArrayName;                  ///< Name of fTracksOutput array
  TClonesArray               *fJetEmbeddingArray;                       //!<! Array of generated jets imported into task (for embedding)
  TString                     fJetEmbeddingArrayName;                   ///< Name of array used to match jets
  TString                     fJetEmbeddingTrackArrayName;              ///< Name of array used to match tracks of jets
  Double_t                    fJetEmbeddingMaxDistance;                 ///< Max distance allowed to accept an embedded jet
  Int_t                       fJetEmbeddingNumMatchedJets;              ///< Number of matched leading jets that will be used
  Bool_t                      fJetEmbeddingUsePerTrackMCPercentage;     ///< When cutting on MC percentage, calculate it per track and not for all MC tracks
  Bool_t                      fJetEmbeddingUseBgrdForMCPercentage;      ///< When cutting on MC percentage, use bgrd. corr to calculate MC percentage
  Bool_t                      fJetEmbeddingCreatePtPlotPerCut;          ///< create TH3 per cut or only once
  std::vector<AliChargedJetsHadronCFCuts> fJetEmbeddingCuts;            ///< Cuts used in jet embedding

  TClonesArray               *fJetVetoArray;                            //!<! Array of jets imported into task used for veto a matching/embedding
  TString                     fJetVetoArrayName;                        ///< Name of array used for veto jets
  Bool_t                      fJetVetoJetByJet;                         ///< If true, the jet veto will be applied on a jet-by-jet basis
  std::vector<AliEmcalJet*>   fMatchedJets;                             ///< Jets matched in an event (embedded)
  std::vector<AliEmcalJet*>   fMatchedJetsReference;                    ///< Jets matched in an event (reference)
  TRandom3*                   fRandom;                                  ///< random number generator


  // Criteria for the selection of jets that are passed to the correlation task
  Int_t                       fJetOutputMode;                           ///< mode which jets are written to array (0: all accepted, 1: leading,  2: subleading, 3: leading+subleading)
  Int_t                       fPythiaExtractionMode;                    ///< Mode which PYTHIA-jets to extract for fJetOutputMode==6: 0: all, 1: quark-jets, 2: gluon jets

  // Event properties
  AliEmcalJet*                fLeadingJet;                              //!<!  leading jet (calculated event-by-event)
  AliEmcalJet*                fSubleadingJet;                           //!<!  subleading jet (calculated event-by-event)
  AliEmcalJet*                fInitialPartonMatchedJet1;                //!<!  On PYTHIA data and fJetOutputMode=6, this holds the PDG code of the initial collisions that was matched to this jet
  AliEmcalJet*                fInitialPartonMatchedJet2;                //!<!  On PYTHIA data and fJetOutputMode=6, this holds the PDG code of the initial collisions that was matched to this jet
  Int_t                       fAcceptedJets;                            //!<!  number accepted jets (calculated event-by-event)
  Int_t                       fAcceptedTracks;                          //!<!  number accepted tracks (calculated event-by-event)

  // ######### HISTOGRAM FUNCTIONS
  void                        FillHistogram(const char * key, Double_t x);
  void                        FillHistogram(const char * key, Double_t x, Double_t y);
  void                        FillHistogram(const char * key, Double_t x, Double_t y, Double_t add);
  void                        FillHistogram3D(const char * key, Double_t x, Double_t y, Double_t z, Double_t add = 0);

  template <class T> T*       AddHistogram1D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis");
  template <class T> T*       AddHistogram2D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0,  const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");
  template <class T> T*       AddHistogram3D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0, Int_t zBins = 100, Double_t zMin = 0.0, Double_t zMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");

  // ######### HELPER FUNCTIONS
  void                        GetTrackMCRatios(AliEmcalJet* jet, AliEmcalJet* mcJet, Double_t& trackRatio, Double_t& ptRatio);
  AliEmcalJet*                GetVetoJet(AliEmcalJet* jet);
  AliEmcalJet*                GetLeadingVetoJet();

  AliEmcalJet*                GetReferenceJet(AliEmcalJet* jet);
  Bool_t                      IsTrackInCone(AliVParticle* track, Double_t eta, Double_t phi, Double_t radius);
  void                        GetInitialCollisionJets();
  void                        GetMatchingJets();

  void                        GetLeadingJets(const char* opt, AliEmcalJet*& jetLeading, AliEmcalJet*& jetSubLeading);
  void                        GetLeadingJetsInArray(TClonesArray* arr, const char* opt, AliEmcalJet*& jetLeading, AliEmcalJet*& jetSubLeading);
  void                        CalculateEventProperties();


 private:
  AliAnalysisTaskChargedJetsHadronCF(const AliAnalysisTaskChargedJetsHadronCF&);            // not implemented
  AliAnalysisTaskChargedJetsHadronCF &operator=(const AliAnalysisTaskChargedJetsHadronCF&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskChargedJetsHadronCF, 9) // Charged jet+h analysis support task
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
    AliBasicJetConstituent() : fEta(0), fPhi(0), fpT(0), fCharge(0), fPID(0) {}
    AliBasicJetConstituent(Float_t eta, Float_t phi, Float_t pt, Short_t charge, Short_t pid)
    : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge), fPID(pid)
    {
    }
    ~AliBasicJetConstituent();

    Double_t Pt()       { return fpT; }
    Double_t Phi()      { return fPhi; }
    Double_t Eta()      { return fEta; }
    Short_t  Charge()   { return fCharge; }
    Short_t  PID()      { return fPID; }

  private:
    Float_t fEta;      ///< eta
    Float_t fPhi;      ///< phi
    Float_t fpT;       ///< pT
    Short_t fCharge;   ///< charge
    Short_t fPID;      ///< most probably PID code/PDG code
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
    AliBasicJet() : fEta(0), fPhi(0), fpT(0), fTruepT(0), fCharge(0), fRadius(0), fArea(0), fPDGCode(0), fBackgroundDensity(0), fMagneticField(0), fVertexX(0), fVertexY(0), fVertexZ(0), fEventID(0), fCentrality(0), fConstituents() {}
    AliBasicJet(Float_t eta, Float_t phi, Float_t pt, Short_t charge, Float_t radius, Float_t area, Float_t partid, Float_t bgrd, Float_t magfield, Float_t vtxX, Float_t vtxY, Float_t vtxZ, Long64_t id, Short_t cent)
    : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge), fRadius(radius), fArea(area), fPDGCode(partid), fBackgroundDensity(bgrd), fMagneticField(magfield), fVertexX(vtxX), fVertexY(vtxY), fVertexZ(vtxZ), fEventID(id), fCentrality(cent), fConstituents()
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
    Int_t                     PDGCode() { return fPDGCode; }
    Double_t                  BackgroundDensity() { return fBackgroundDensity; }
    Double_t                  MagneticField() { return fMagneticField; }
    Double_t                  VertexX() { return fVertexX; }
    Double_t                  VertexY() { return fVertexY; }
    Double_t                  VertexZ() { return fVertexZ; }
    Long64_t                  EventID() { return fEventID; }
    Short_t                   Centrality() { return fCentrality; }
    Int_t                     GetNumbersOfConstituents() { return fConstituents.size(); }

    // Basic constituent functions
    AliBasicJetConstituent*   GetJetConstituent(Int_t index) { return &fConstituents[index]; }
    void                      AddJetConstituent(Float_t eta, Float_t phi, Float_t pt, Short_t charge, Short_t pid=0)
    {
      AliBasicJetConstituent c (eta, phi, pt, charge, pid);
      AddJetConstituent(&c);
    }
    void                      AddJetConstituent(AliBasicJetConstituent* constituent) {fConstituents.push_back(*constituent); }
    void                      SetTruePt(Double_t val)  {fTruepT = val;}


  private:
    Float_t   fEta;      ///< eta
    Float_t   fPhi;      ///< phi
    Float_t   fpT;       ///< pT
    Float_t   fTruepT;   ///< true pT (optional, e.g. from matching)
    Short_t   fCharge;   ///< charge
    Float_t   fRadius;   ///< jet radius
    Float_t   fArea;     ///< jet area
    Int_t     fPDGCode;  ///< PDG code of source particle
    Float_t   fBackgroundDensity; ///< background
    Float_t   fMagneticField; ///< event magnetic field
    Float_t   fVertexX; ///< event vertex X
    Float_t   fVertexY; ///< event vertex Y
    Float_t   fVertexZ; ///< event vertex Z
    Long64_t  fEventID;  ///< Unique event id
    Short_t   fCentrality; ///< centrality

    std::vector<AliBasicJetConstituent> fConstituents; ///< vector of constituents

};


#endif
