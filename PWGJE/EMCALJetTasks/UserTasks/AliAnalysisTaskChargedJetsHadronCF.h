#ifndef ALIANALYSISTASKCHARGEDJETSHADRONCF_H
#define ALIANALYSISTASKCHARGEDJETSHADRONCF_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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

class THn;

class AliAnalysisTaskChargedJetsHadronCF : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskChargedJetsHadronCF();
  AliAnalysisTaskChargedJetsHadronCF(const char *name);
  virtual ~AliAnalysisTaskChargedJetsHadronCF();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  // ######### SETTERS/GETTERS
  void                        SetJetParticleArrayName(const char* name)   { fJetParticleArrayName = name; }
  void                        SetTrackParticleArrayName(const char* name) { fTrackParticleArrayName = name; }
  void                        SetJetMatchingArrayName(const char* name)   { fJetMatchingArrayName = name; }

  void                        ActivateRejectionJetConstituents(TF1* rejectionFunction)
                                { fRejectionFunction = rejectionFunction; }
  void                        SetJetOutputMode(Int_t mode) {fJetOutputMode = mode;}

  void                        SetEventCriteriumBackground(Double_t minValue, Double_t maxValue)   {fEventCriteriumMinBackground = minValue; fEventCriteriumMaxBackground = maxValue;}
  void                        SetEventCriteriumLeadingJets(Double_t leading, Double_t subleading) {fEventCriteriumMinLeadingJetPt = leading; fEventCriteriumMinSubleadingJetPt = subleading;}
  void                        SetEventCriteriumSelection(Int_t type);

  void                        ActivateJetExtraction(Double_t percentage, Double_t minPt) {fExtractionPercentage = percentage; fExtractionMinPt = minPt;}

 protected:
  void                        ExecOnce();
  Bool_t                      Run();

  // ######### META FUNCTIONS
  void                        BinLogAxis(const THn *h, Int_t axisNumber);
  void                        AddJetToTree(AliEmcalJet* jet);
  void                        AddJetToOutputArray(AliEmcalJet* jet);
  void                        AddTrackToOutputArray(AliVTrack* track);
  void                        FillHistogramsTracks(AliVTrack* track);
  void                        FillHistogramsJets(AliEmcalJet* jet);
  void                        FillHistogramsJetConstituents(AliEmcalJet* jet);
  Bool_t                      IsJetSelected(AliEmcalJet* jet);
  Bool_t                      IsEventSelected();

  AliJetContainer            *fJetsCont;                                //!<! Jets
  AliTrackContainer          *fTracksCont;                              //!<! Tracks
  TTree*                      fJetsTree;                                //!<! Jets that will be saved to a tree (optionally)
  void*                       fJetsTreeBuffer;                          //!<!  buffer for one jet (that will be saved to the tree)
  Double_t                    fExtractionPercentage;                    /// percentage that is recorded
  Double_t                    fExtractionMinPt;                         /// minimum pt of recorded jets
  Int_t                       fNumberOfCentralityBins;                  /// Number of centrality bins
  TClonesArray               *fJetsOutput;                              //!<! Array of basic correlation particles attached to the event (jets)
  TClonesArray               *fTracksOutput;                            //!<! Array of basic correlation particles attached to the event (tracks)
  TClonesArray               *fJetsInput;                               //!<! Array of generated jets imported into task (toy model)
  TString                     fJetParticleArrayName;                    /// Name of fJetsOutput array
  TString                     fTrackParticleArrayName;                  /// Name of fTracksOutput array
  TString                     fJetMatchingArrayName;                    /// Name of array used to match jets
  TRandom3*                   fRandom;                                  /// random number generator

  // Criteria for the selection of jets that are passed to the correlation task
  TF1*                        fRejectionFunction;                       /// Function describing the cut applied to jet const.
  Int_t                       fJetOutputMode;                           /// mode which jets are written to array (0: all accepted, 1: leading,  2: subleading, 3: leading+subleading)
  Double_t                    fMinFakeFactorPercentage;                 /// min fake factor (percentage relative to cut profile)
  Double_t                    fMaxFakeFactorPercentage;                 /// max fake factor (percentage relative to cut profile)
  Int_t                       fEventCriteriumMode;                      /// Mode of event selection
  Double_t                    fEventCriteriumMinBackground;             /// Minimum background
  Double_t                    fEventCriteriumMaxBackground;             /// Maximum background
  Double_t                    fEventCriteriumMinLeadingJetPt;           /// Min leading jet
  Double_t                    fEventCriteriumMinSubleadingJetPt;        /// Min subleading jet

  // Event properties
  AliEmcalJet*                fLeadingJet;                              //!<!  leading jet (calculated event-by-event)
  AliEmcalJet*                fSubleadingJet;                           //!<!  subleading jet (calculated event-by-event)
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
  Double_t                    CalculateFakeFactor(AliEmcalJet* jet);
  AliEmcalJet*                GetSubleadingJet(const char* opt);
  void                        CalculateEventProperties();


 private:
  AliAnalysisTaskChargedJetsHadronCF(const AliAnalysisTaskChargedJetsHadronCF&);            // not implemented
  AliAnalysisTaskChargedJetsHadronCF &operator=(const AliAnalysisTaskChargedJetsHadronCF&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskChargedJetsHadronCF, 6) // Charged jet+h analysis support task
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
class AliBasicJetConstituent : public TObject
{
  public:
    AliBasicJetConstituent() : fEta(0), fPhi(0), fpT(0), fCharge(0) {}
    AliBasicJetConstituent(Float_t eta, Float_t phi, Float_t pt, Short_t charge)
    : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge)
    {
    }
    ~AliBasicJetConstituent();

    Bool_t   IsEqual(const TObject* obj) { return (obj->GetUniqueID() == GetUniqueID()); }
    Double_t Pt()       { return fpT; }
    Double_t Phi()      { return fPhi; }
    Double_t Eta()      { return fEta; }
    Short_t  Charge()   { return fCharge; }

  private:
    Float_t fEta;      /// eta
    Float_t fPhi;      /// phi
    Float_t fpT;       /// pT
    Short_t fCharge;   /// charge

  /// \cond CLASSIMP
  ClassDef( AliBasicJetConstituent, 1); // very basic jet constituent object
  /// \endcond
};

//###############################################################################################################################################3
/**
 * \class AliBasicJet
 * \brief Simple class containing basic information for a jet
 *
 * This class is used to save jets including constituents with minimize memory consumption.
 * Saved information focus on correlation analyses
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date Apr 21, 2016
 */
// 
class AliBasicJet : public TObject
{
  public:
    AliBasicJet() : fEta(0), fPhi(0), fpT(0), fCharge(0), fRadius(0), fArea(0), fBackgroundDensity(0), fEventID(0), fCentrality(0), fConstituents() {}
    AliBasicJet(Float_t eta, Float_t phi, Float_t pt, Short_t charge, Float_t radius, Float_t area, Float_t bgrd, Long64_t id, Short_t cent)
    : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge), fRadius(radius), fArea(area), fBackgroundDensity(bgrd), fEventID(id), fCentrality(cent), fConstituents()
    {
    }
    ~AliBasicJet();

    // Basic jet properties
    Bool_t                    IsEqual(const TObject* obj) { return (obj->GetUniqueID() == GetUniqueID()); }
    Double_t                  Pt()       { return fpT; }
    Double_t                  Phi()      { return fPhi; }
    Double_t                  Eta()      { return fEta; }
    Short_t                   Charge()   { return fCharge; }
    Double_t                  Radius() { return fRadius; }
    Double_t                  Area() { return fArea; }
    Double_t                  BackgroundDensity() { return fBackgroundDensity; }
    Long64_t                  EventID() { return fEventID; }
    Short_t                   Centrality() { return fCentrality; }
    Int_t                     GetNumbersOfConstituents() { return fConstituents.size(); }

    // Basic constituent functions
    AliBasicJetConstituent*   GetJetConstituent(Int_t index) { return &fConstituents[index]; }
    void                      AddJetConstituent(Float_t eta, Float_t phi, Float_t pt, Short_t charge)
    {
      AliBasicJetConstituent c (eta, phi, pt, charge);
      AddJetConstituent(&c);
    }
    void                      AddJetConstituent(AliBasicJetConstituent* constituent) {fConstituents.push_back(*constituent); }


  private:
    Float_t   fEta;      /// eta
    Float_t   fPhi;      /// phi
    Float_t   fpT;       /// pT
    Short_t   fCharge;   /// charge
    Float_t   fRadius;   /// jet radius
    Float_t   fArea;     /// jet area
    Float_t   fBackgroundDensity; /// background
    Long64_t  fEventID;  /// Unique event id
    Short_t   fCentrality; /// centrality

    std::vector<AliBasicJetConstituent> fConstituents; /// vector of constituents

  /// \cond CLASSIMP
  ClassDef( AliBasicJet, 2); // very basic jet object
  /// \endcond
};
#endif
