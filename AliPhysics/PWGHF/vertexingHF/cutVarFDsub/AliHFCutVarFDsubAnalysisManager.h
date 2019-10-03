#ifndef ALIHFCUTVARFDSUBANALYSISMANAGER_H
#define ALIHFCUTVARFDSUBANALYSISMANAGER_H
/// \class AliHFCutVarFDsubAnalysisManager
/// \brief Analysis manager for the cut variation feed down method analysis
///
///
///
///
/// \author Felix Reidt <felix.reidt@cern.ch>, CERN
/// \author Fabrizio Grosa <grosa@to.infn.it>, INFN Torino
/// \date Aug 17, 2015

#include "TObject.h"
#include "THnSparse.h"
#include "TString.h"

class TList;
class TH1F;

class AliHFCutVarFDsubAnalysisManager : public TObject {
protected:
  THnSparseF* fMCgenLevel[2]; //!<!
  THnSparseF* fMCafterCuts[2]; //!<!
  THnSparseF* fData; //!<!
  TList* fAxes; //!<!
  TList* fCuts; //!<!
  TList* fEffListVsCutSets[2]; //!<!
  TList* fEffListVsBins[2]; //!<!
  TList* fRawYields; //!<!
  TH1F* fCorrYieldPrompt; //!<!
  TH1F* fCorrYieldFD; //!<!
  TList* fResiduals; //!<!
  TList* fPulls; //!<!
  TList* fFprompt; //!<!
  TList* fFpromptRaw; //!<!
  TList* fIncDistError; //!<!

  TString fxAxisTitle; // title of the x-axis
  Double_t* fBinsX;    //!<! array containing the x-axis bins

  AliHFCutVarFDsubAnalysisManager(const AliHFCutVarFDsubAnalysisManager& analysisManager); /// Copy constructor
  AliHFCutVarFDsubAnalysisManager operator=(const AliHFCutVarFDsubAnalysisManager& analysisManager); // Assignment operator

public:
  enum { kPrompt=0, kFD };

  AliHFCutVarFDsubAnalysisManager(); /// Default constructor
  ~AliHFCutVarFDsubAnalysisManager(); /// Destructor


  void DrawDistributions(TString strOutputFolder="."); ///Draw the distributions of the cut variables
  void GetEfficiencies(TString strOutputFolder=".",TString strOutputFile="Efficiency.root",
                       Bool_t ptWeight=kFALSE, TF1* funcWeightsD=0x0, TF1* funcWeightsB=0x0); /// Obtain the efficiences from the MC THnSparses
  void DrawEfficiencies(TString strOutputFolder, TString prefix="eff", UInt_t xAxis=0);
  void GetRawYields(Bool_t drawFit=kFALSE, TString strOutputFolder=".", TString strOutputFile="RawYields.root"); /// Obtain the invariant mass distributions
  void GetXaxisInformation(); /// Obtain the x-axis information
  Bool_t Minimise(UInt_t method=0, UInt_t nIterations=10, Bool_t useWeights=kTRUE,
                  Double_t relSysteEffErr=0., Int_t nSim=1000); /// Obtain the corrected yields
  void DrawLines(TString strOutputFolder, Bool_t IsIncentre=kFALSE); ///Draw the lines correspondent to the cut sets

  /// \cond CLASSDEF
  ClassDef(AliHFCutVarFDsubAnalysisManager, 1);
  /// \endcond
};
#endif // ALIHFCUTVARFDSUBANALYSISMANAGER_H
