#ifndef AliPerfAnalyzeInvPt_h
#define AliPerfAnalyzeInvPt_h


class TNamed;
class TString;
class TList;
class TH1F;
class TH2F;
class TGraphErrors;

class AliPerfAnalyzeInvPt : public TNamed {

 public:
  AliPerfAnalyzeInvPt();
  AliPerfAnalyzeInvPt(Char_t* name, Char_t* title);
   virtual ~AliPerfAnalyzeInvPt(){;}
  void InitHistos(Double_t *binsXTheta,Double_t *fitParamTheta,Double_t *errFitParamTheta,Double_t *binsXPhi,Double_t *fitParamPhi,Double_t *errFitParamPhi);
   void InitFitFcn();
  void StartAnalysis(TH2F *histThetaInvPt, TH2F *histPhiInvPt, TObjArray* aFolderObj);
  void StartAnalysis(TObjArray *aFolderObj);
  void SetProjBinsPhi(const Double_t *pBins,Int_t sizep);
  void SetProjBinsTheta(const Double_t *tBins, Int_t sizet);
  void SetMakeFitOption(const Bool_t setGausFit, const Double_t exclusionR,const Double_t fitR );
 
  
   protected:
  Double_t fThetaBins[100];
  Double_t fPhiBins[100];
  
  Int_t fNThetaBins;
  Int_t fNPhiBins ;
  Double_t fRange;
  Double_t fExclRange ;
  Bool_t fFitGaus ;

 // projection histos
  TH1F *fHistFitTheta[100];
  TH1F *fHistFitPhi[100];

  TH2F *fHistH2InvPtTheta;
  TH2F *fHistH2InvPtPhi; 
  TGraphErrors *fGrMinPosTheta;
  TGraphErrors *fGrMinPosPhi;

 

 private:
  TF1 *fFitMinPos;
  TF1 *fFitMinPosRejP;
  TF1 *fFitInvGauss;
  TF1 *fFitInvGaussRejP;

 static Double_t Polynomial(Double_t *x, Double_t *par);
 static Double_t PolynomialRejP(Double_t *x, Double_t *par);
 static Double_t InvGauss(Double_t *x, Double_t *par);
 static Double_t InvGaussRejP(Double_t *x, Double_t *par);
 
  void MakeFit(TH1 *dmproy, TF1 * fitpb, Double_t &mean, Double_t &ErrMean,  Double_t &excl,Double_t &range);
  void MakeFitBetter(TH1 *dmproy, TF1 * fitpb2, Double_t &mean, Double_t &ErrMean, Double_t &f, Double_t &excl, Double_t &range);
  void MakeFitInvGauss(TH1 *dmproy, TF1 * fitpb2, Double_t &mean, Double_t &ErrMean,Double_t &excl , Double_t &range);
  void MakeFitInvGaussBetter(TH1 *dmproy, TF1 * fitpb2, Double_t &mean, Double_t &ErrMean, Double_t &f,Double_t &excl,  Double_t &range);
 

  

  AliPerfAnalyzeInvPt(const AliPerfAnalyzeInvPt&);            // not implemented 
  AliPerfAnalyzeInvPt& operator=(const AliPerfAnalyzeInvPt&); // not implemented 

  ClassDef(AliPerfAnalyzeInvPt, 1); 
};

#endif
