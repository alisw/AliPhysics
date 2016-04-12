#ifndef ALIANALYSISTASKCHARGEDJETSHADRONCF_H
#define ALIANALYSISTASKCHARGEDJETSHADRONCF_H

// $Id$

class TH1;
class TH2;
class TH3;
class AliJetContainer;

#include "AliAnalysisTaskEmcalJet.h"

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


  void                        ActivateFakejetRejection(Double_t minFakeFactorPercentage, Double_t maxFakeFactorPercentage, TProfile* cutProfile)
                                { fMinFakeFactorPercentage = minFakeFactorPercentage; fMaxFakeFactorPercentage = maxFakeFactorPercentage; hFakeFactorCutProfile = cutProfile; fUseFakejetRejection = kTRUE; }
  void                        SetMinFakeFactorPercentage(Double_t value)  { fMinFakeFactorPercentage = value; }
  void                        SetMaxFakeFactorPercentage(Double_t value)  { fMaxFakeFactorPercentage = value; }
  void                        SetUseFakejetRejection(Bool_t value)  { fUseFakejetRejection = value; }
  void                        SetJetOutputMode(Int_t mode) {fJetOutputMode = mode;}

  void                        SetEventCriteriumBackground(Double_t minValue, Double_t maxValue)   {fEventCriteriumMinBackground = minValue; fEventCriteriumMaxBackground = maxValue;}
  void                        SetEventCriteriumLeadingJets(Double_t leading, Double_t subleading) {fEventCriteriumMinLeadingJetPt = leading; fEventCriteriumMinSubleadingJetPt = subleading;}
  void                        SetEventCriteriumSelection(Int_t type);

 protected:
  void                        ExecOnce();
  Bool_t                      Run();

  // ######### META FUNCTIONS
  void                        AddJetToOutputArray(AliEmcalJet* jet);
  void                        AddTrackToOutputArray(AliVTrack* track);
  void                        FillHistogramsTracks(AliVTrack* track);
  void                        FillHistogramsJets(AliEmcalJet* jet);
  void                        FillHistogramsJetConstituents(AliEmcalJet* jet);
  Bool_t                      IsJetSelected(AliEmcalJet* jet);
  Bool_t                      IsEventSelected();


  AliJetContainer            *fJetsCont;                                //!Jets
  AliTrackContainer          *fTracksCont;                              //!Tracks
  Int_t                       fNumberOfCentralityBins;                  // Number of centrality bins
  TClonesArray               *fJetsOutput;                              //!Array of basic correlation particles attached to the event (jets)
  TClonesArray               *fTracksOutput;                            //!Array of basic correlation particles attached to the event (tracks)
  TClonesArray               *fJetsInput;                               //!Array of generated jets imported into task (toy model)
  TString                     fJetParticleArrayName;                    // Name of fJetsOutput array
  TString                     fTrackParticleArrayName;                  // Name of fTracksOutput array
  TString                     fJetMatchingArrayName;                    // Name of array used to match jets
  TRandom3*                   fRandom;                                  // random number generator

  // Criteria for the selection of jets that are passed to the correlation task
  TProfile*                   hFakeFactorCutProfile;                    // profile of the fake factor cut distribution
  Int_t                       fJetOutputMode;                           // mode which jets are written to array (0: all accepted, 1: leading,  2: subleading, 3: leading+subleading)
  Bool_t                      fUseFakejetRejection;                     // Use fakejet rejection (a la ATLAS)
  Double_t                    fMinFakeFactorPercentage;                 // min fake factor (percentage relative to cut profile)
  Double_t                    fMaxFakeFactorPercentage;                 // max fake factor (percentage relative to cut profile)
  Int_t                       fEventCriteriumMode;                      // Mode of event selection
  Double_t                    fEventCriteriumMinBackground;             // Minimum background
  Double_t                    fEventCriteriumMaxBackground;             // Maximum background
  Double_t                    fEventCriteriumMinLeadingJetPt;           // Min leading jet
  Double_t                    fEventCriteriumMinSubleadingJetPt;        // Min subleading jet

  // Event properties
  AliEmcalJet*                fLeadingJet;                              //! leading jet (calculated event-by-event)
  AliEmcalJet*                fSubleadingJet;                           //! subleading jet (calculated event-by-event)
  Int_t                       fAcceptedJets;                            //! number accepted jets (calculated event-by-event)
  Int_t                       fAcceptedTracks;                          //! number accepted tracks (calculated event-by-event)

  // ######### HISTOGRAM FUNCTIONS
  void                        FillHistogram(const char * key, Double_t x);
  void                        FillHistogram(const char * key, Double_t x, Double_t y);
  void                        FillHistogram(const char * key, Double_t x, Double_t y, Double_t add);
  template <class T> T*       AddHistogram1D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis");
  template <class T> T*       AddHistogram2D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0,  const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");

  // ######### HELPER FUNCTIONS
  Double_t                    CalculateFakeFactor(AliEmcalJet* jet);
  AliEmcalJet*                GetSubleadingJet(const char* opt);
  void                        CalculateEventProperties();


 private:
  AliAnalysisTaskChargedJetsHadronCF(const AliAnalysisTaskChargedJetsHadronCF&);            // not implemented
  AliAnalysisTaskChargedJetsHadronCF &operator=(const AliAnalysisTaskChargedJetsHadronCF&); // not implemented

  ClassDef(AliAnalysisTaskChargedJetsHadronCF, 4) // Charged jet+h analysis support task
};
#endif
