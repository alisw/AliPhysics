#include "TObject.h"

class AliCorrelationAnalysis : public TObject
{
public:
  AliCorrelationAnalysis();
  ~AliCorrelationAnalysis();

  void FillParentTHnSparse(TString fileName, Bool_t reduce = kFALSE, const TString &tag = "");
  void PlotDeltaPhiEtaGap(const TString &fileNamePbPb, TString fileNamePbPbMix = "",
			  const TString &fileNamepp = "", const TString &fileNamepp2 = "",
			  const TString &outputFile = "dphi_corr.root");
  void MergeDPhiFiles(const TString &fileName, const TString &fileName2, const TString &target);
  void RemoveWing(const TString &fileName, const TString &outputFile);

protected:
  void* GetUEHistogram(const TString &fileName, TList** listRef = 0, Bool_t mixed = kFALSE, const char* tag = "");
  void SetupRanges(void* obj);
  Double_t GetEtaCut(TTree *analysisSettings);
  void GetSumOfRatios(void* hVoid, void* hMixedVoid, TH1** hist, AliUEHist::CFStep step,
		      Int_t centralityBegin, Int_t centralityEnd, Float_t ptBegin, Float_t ptEnd,
		      Bool_t normalizePerTrigger = kTRUE, Bool_t useCentralityBinsDirectly = kFALSE);

  TString fCurrentFileName; //! current file name for caching

  void *fCacheSameEvent;    //! cache for same event
  void *fCacheMixedEvent;   //! cache for mixed event

  Float_t fPtMin;
  Float_t fPtMax;
  Float_t fZVtxRange;

  ClassDef(AliCorrelationAnalysis, 1);
};
