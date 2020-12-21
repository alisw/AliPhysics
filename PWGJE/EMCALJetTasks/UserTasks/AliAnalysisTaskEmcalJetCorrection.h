#ifndef ALIANALYSISTASKEMCALJETCORRECTION_H
#define ALIANALYSISTASKEMCALJETCORRECTION_H

/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class AliAnalysisTaskEmcalJet;

//###############################################################################################################################################3
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  class TPython;
#endif

/**
 * \class AliAnalysisTaskEmcalJetCorrection
 * \brief Analysis task that corrects jet pT with scikit-learn models
 *
 * This task changes the pT of the input jet collection to a corrected value
 * which is calculated with a sklearn estimator
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, Yale
 * \date Feb 07, 2019
 */
// 
class AliAnalysisTaskEmcalJetCorrection : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalJetCorrection();
  AliAnalysisTaskEmcalJetCorrection(const char *name);
  virtual ~AliAnalysisTaskEmcalJetCorrection();
  static AliAnalysisTaskEmcalJetCorrection* AddTaskEmcalJetCorrection(TString modelName, TString trackArray, TString jetArray, TString rhoObject, Double_t jetRadius, const char* taskNameSuffix);
  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);
  void                        SetModelName(const char* val)                               {fModelName = val;}
  void                        SetBackgroundModelFileName(const char* val)                 {fBackgroundModelFileName = val;}
  void                        SetBackgroundModelInputParameters(const char* val)          {fBackgroundModelInputParameters = val;}
  void                        SetCustomPackages(const char* val)                          {fCustomPackages = val;}
  void                        SetPythonModulePath(const char* val)                        {fPythonModulePath = val;}

protected:
  void                        ExecOnce();
  Bool_t                      Run();
  void                        GetPtFromModel(AliEmcalJet* jet, Float_t& pt_ML);
  TString                     GetBackgroundModelArrayString(AliEmcalJet* jet);
  void                        CalculateJetShapes(AliEmcalJet* jet, Double_t& leSub_noCorr, Double_t& angularity, Double_t& momentumDispersion, Double_t& trackPtMean, Double_t& trackPtMedian);

  #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  TPython*                    fPythonCLI;
  #endif
  TString                     fCustomPackages;                          ///< Custom python modules to be installed locally
  TString                     fPythonModulePath;                        ///< The path of custom python modules (depends on local python version)
  AliJetContainer            *fJetsCont;                                //!<! Jets
  AliParticleContainer       *fTracksCont;                              //!<! Tracks
  TClonesArray*               fJetOutputArray;                          //!<! Array of corr. jets, attached to event
  TString                     fBackgroundModelFileName;                 ///< MVA model file name
  TString                     fBackgroundModelInputParameters;          ///< MVA model input parameters (comma-separated)
  TString                     fModelName;                               ///< Name of model (used for jet collection suffix)

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
  AliAnalysisTaskEmcalJetCorrection(const AliAnalysisTaskEmcalJetCorrection&);            // not implemented
  AliAnalysisTaskEmcalJetCorrection &operator=(const AliAnalysisTaskEmcalJetCorrection&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetCorrection, 2) // Jet correction task
  /// \endcond
};

#endif
