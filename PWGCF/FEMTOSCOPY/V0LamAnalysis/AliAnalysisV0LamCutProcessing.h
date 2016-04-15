#ifndef AliAnalysisV0LamCutProcessing_cxx
#define AliAnalysisV0LamCutProcessing_cxx

#include "AliAnalysisV0LamEventCollection.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"

class AliAnalysisV0LamCut
{
 public:
  AliAnalysisV0LamCut();
  AliAnalysisV0LamCut(std::vector<Double_t> variableCutValues, bool isUpperBound);
  ~AliAnalysisV0LamCut();
  int fNumberOfCutValues;
  bool fIsAnUpperLimit; //if this is false, the cut is a lower limit
  std::vector<Double_t> fCutValues;
};

class AliAnalysisV0LamCutProcessing
{
  //created and used in AliAnalysisV0Lam::init
public:
  enum CutType_t {kNoCut = 0,
		  kDCAProton = 1,
		  kDCAPion = 2,
		  kDCADaughter = 3,
		  kProperDecayLength = 4,
		  kEta = 5,
		  kCosPointing = 6,
		  kDCAV0 = 7,
		  kPt = 8,
		  kMassLam = 9,
		  kMassALam = 10};
  
  AliAnalysisV0LamCutProcessing(TList *const outputList, Int_t varCutIndex);
  ~AliAnalysisV0LamCutProcessing();
  void SetCentralityBin(int centBin) {fCurrentCentralityBin = centBin;}
  void CheckIfV0PassesCuts(AliReconstructedV0 *v0);
  void DoV0Histogramming(AliReconstructedV0 *v0);
  int GetNumberOfVariableCutValues() const {return fNumberOfVariableCutValues;}
  
 private:
  AliAnalysisV0LamCutProcessing(const AliAnalysisV0LamCutProcessing &processor); // Not implemented
  AliAnalysisV0LamCutProcessing operator=(const AliAnalysisV0LamCutProcessing &processor); // Not implemented
  void InitHistograms();
  void SortAndFillCutHistograms(AliReconstructedV0 *v0, bool isLambda);
  void ProcessCut(AliReconstructedV0 *v0, int index, bool isLambdaCandidate);
  void DetermineIfTrueV0(AliReconstructedV0 *v0, bool isLambda);
  void FillHist(AliReconstructedV0 *v0, int cutTypeIndex, int variableCutIndex, bool isLambda);

  std::vector<AliAnalysisV0LamCut> fCuts;
  int fNumberOfCutTypes;
  int fNumberOfVariableCutValues;
  int fCurrentCentralityBin;
  int fVariableCutType;
  TList *fOutputList; //! Compact output list where histograms are written
  //These hists are filled during the V0 reconstruction process
  TH2F *fHistDaughterPosDcaToPrimLam;
  TH2F *fHistDaughterPosDcaToPrimALam;
  TH2F *fHistDaughterNegDcaToPrimLam;
  TH2F *fHistDaughterNegDcaToPrimALam;
  TH2F *fHistDaughtersDcaLam;
  TH2F *fHistDaughtersDcaALam;
  TH2F *fHistDecayLengthLam;
  TH2F *fHistDecayLengthALam;
  TH2F *fHistProperDecayLengthLam;
  TH2F *fHistProperDecayLengthALam;
  TH2F *fHistEtaLam;
  TH2F *fHistEtaALam;
  TH2F *fHistCosPointingLam;
  TH2F *fHistCosPointingALam;
  TH2F *fHistDcaLam;
  TH2F *fHistDcaALam;
  TH2F *fHistPtLam;
  TH2F *fHistPtALam;
  TH2F *fHistMassLam; 
  TH2F *fHistMassALam;
  TH3F *fHistMassCentralityLam;
  TH3F *fHistMassCentralityALam;
};

#endif
