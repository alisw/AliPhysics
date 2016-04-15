#ifndef AliAnalysisTaskEmcalJetExtractor_H
#define AliAnalysisTaskEmcalJetExtractor_H

// $Id$

class TH1;
class TH2;
class TH3;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliRhoParameter;

#include <vector>
#include "AliAnalysisTaskEmcalJet.h"


//###############################################################################################################################################3


class AliAnalysisTaskEmcalJetExtractor : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalJetExtractor();
  AliAnalysisTaskEmcalJetExtractor(const char *name);
  virtual ~AliAnalysisTaskEmcalJetExtractor();

  void                        UserCreateOutputObjects();
  void                        Initialize(Int_t modus, const char* treeName);
  void                        DefineExtraction(Int_t type, Int_t criterium, Double_t minPt, Double_t maxPt, Double_t percentage);

 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;

  AliJetContainer            *fJetsCont;                                //! Jets
  AliTrackContainer          *fTracksCont;                              //! Tracks
  void*                       fJetBuffer;                               //! buffer for one jet (that will be saved to the tree)
  TTree*                      fJetsOutput;                              //! Jets that will be saved to a tree
  Int_t                       fCounter;                                 // Event counter
  TRandom3*                   fRandom;                                  // Randomizer object
  Bool_t                      fIsExtractionDefined;                     // trigger if extraction is defined
  Int_t                       fExtractionType;                          // specifies how the jets are saved. 0-AliEmcalJet,1-AliBasicJet,2-AliBasicJet w/constituents
  Int_t                       fExtractionCriterium;                     // criterium for extraction: 0-MinBias
  Double_t                    fExtractionMinPt;                         // Jet min pt for extraction
  Double_t                    fExtractionMaxPt;                         // Jet max pt for extraction
  Double_t                    fExtractionPercentage;                    // extraction percentage, rest is thrown away

 private:
  AliAnalysisTaskEmcalJetExtractor(const AliAnalysisTaskEmcalJetExtractor&);            // not implemented
  AliAnalysisTaskEmcalJetExtractor &operator=(const AliAnalysisTaskEmcalJetExtractor&); // not implemented

  void                        PrintSettings();                          // Show all relevant settings (to be executed once after task is started)

  // ######### HISTOGRAM HELPER FUNCTIONS
  void                        FillHistogram(const char * key, Double_t x);
  void                        FillHistogram(const char * key, Double_t x, Double_t y);
  void                        FillHistogram(const char * key, Double_t x, Double_t y, Double_t add);
  template <class T> T*       AddHistogram1D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis");
  template <class T> T*       AddHistogram2D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0,  const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");


  ClassDef(AliAnalysisTaskEmcalJetExtractor, 1) // Jet extractor task
};

//###############################################################################################################################################3

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
    Float_t fEta;      // eta
    Float_t fPhi;      // phi
    Float_t fpT;       // pT
    Short_t fCharge;   // charge

    ClassDef( AliBasicJetConstituent, 1); // very basic jet constituent object
};

//###############################################################################################################################################3

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
    Float_t   fEta;      // eta
    Float_t   fPhi;      // phi
    Float_t   fpT;       // pT
    Short_t   fCharge;   // charge
    Float_t   fRadius;   // jet radius
    Float_t   fArea;     // jet area
    Float_t   fBackgroundDensity; // background
    Long64_t  fEventID;  // Unique event id
    Short_t   fCentrality; // centrality

    std::vector<AliBasicJetConstituent> fConstituents; // vector of constituents

    ClassDef( AliBasicJet, 2); // very basic jet object
};

#endif
