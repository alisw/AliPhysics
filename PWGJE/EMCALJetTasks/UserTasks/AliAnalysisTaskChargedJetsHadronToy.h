#ifndef ALIANALYSISTASKCHARGEDJETSHADRONTOY_H
#define ALIANALYSISTASKCHARGEDJETSHADRONTOY_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TClonesArray;
class TString;
class TRandom3;

// Toy model to create an event containing tracks

class AliAnalysisTaskChargedJetsHadronToy : public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskChargedJetsHadronToy();
  AliAnalysisTaskChargedJetsHadronToy(const char* name);
  virtual ~AliAnalysisTaskChargedJetsHadronToy();

  void  UserCreateOutputObjects();
  void  Terminate(Option_t *);

  // ### SETTERS/GETTERS
  void                        SetDistributionMultiplicity(TH1* val)           {fDistributionMultiplicity = val;}
  void                        SetDistributionPt(TH1* val)                     {fDistributionPt = val;}
  void                        SetDistributionEtaPhi(TH1* val)                 {fDistributionEtaPhi = val;}
  void                        SetDistributionV2(TH2* val)                     {fDistributionV2 = val;}
  void                        SetDistributionV3(TH2* val)                     {fDistributionV3 = val;}
  void                        SetDistributionV4(TH2* val)                     {fDistributionV4 = val;}
  void                        SetDistributionV5(TH2* val)                     {fDistributionV5 = val;}
  void                        SetRangeCentrality(Int_t min, Int_t max)        {fMinCentrality = min; fMaxCentrality = max;}

  void                        SetInputArrayName(const char* val)              {fInputArrayName = val;}
  void                        SetOutputArrayName(const char* val)             {fOutputArrayName = val;}

protected:
  // ### Settings
  TH1*                        fDistributionMultiplicity;          // histogram for multiplicity distribution
  TH1*                        fDistributionPt;                    // histogram for Pt distribution
  TH1*                        fDistributionEtaPhi;                // histogram for eta/phi distribution
  Int_t                       fMinCentrality;                     // minimum centrality
  Int_t                       fMaxCentrality;                     // maximum centrality

  TH2*                        fDistributionV2;                    /// Distribution for v2 in bins of pt and centrality
  TH2*                        fDistributionV3;                    /// Distribution for v3 in bins of pt and centrality
  TH2*                        fDistributionV4;                    /// Distribution for v4 in bins of pt and centrality
  TH2*                        fDistributionV5;                    /// Distribution for v5 in bins of pt and centrality

  // ### Input/output settings+arrays
  TString                     fInputArrayName;                    // Name of the TClonesArray that will be loaded
  TString                     fOutputArrayName;                   // Name of the TClonesArray that will be loaded
  TClonesArray*               fInputArray;                        //! input array containing tracks from events
  TClonesArray*               fOutputArray;                       //! input array containing tracks from events

  // ### Misc
  TRandom3*                   fRandom;                            //! random number generator
  Double_t                    fToyCent;                           /// eventwise centrality value
  Double_t                    fRandomPsi3;                        /// eventwise calculated psi 3
  Double_t                    fRandomPsi4;                        /// eventwise calculated psi 4
  Double_t                    fRandomPsi5;                        /// eventwise calculated psi 5

  Bool_t                      Run();
  void                        AssembleEvent();
  void                        CreateQAPlots();
  Double_t                    AddFlow(Double_t phi, Double_t pt);
  void                        ExecOnce();

  void                        FillHistogram(const char * key, Double_t x);
  void                        FillHistogram(const char * key, Double_t x, Double_t y);
  void                        FillHistogram(const char * key, Double_t x, Double_t y, Double_t add);
  template <class T> T*       AddHistogram1D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis");
  template <class T> T*       AddHistogram2D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0,  const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");

private:
  AliAnalysisTaskChargedJetsHadronToy(const AliAnalysisTaskChargedJetsHadronToy&);            // not implemented
  AliAnalysisTaskChargedJetsHadronToy &operator=(const AliAnalysisTaskChargedJetsHadronToy&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskChargedJetsHadronToy, 2); // Toy model
  /// \endcond
};

#endif
