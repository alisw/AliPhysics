#ifndef ALIANALYSISTASKSEPPbCORRELATIONSJETV2KINE_H
#define ALIANALYSISTASKSEPPbCORRELATIONSJETV2KINE_H

#include "AliAnalysisTaskSE.h"

class TList;
class AliAnalysisTaskSEpPbCorrelationsJetV2Kine : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskSEpPbCorrelationsJetV2Kine();
  AliAnalysisTaskSEpPbCorrelationsJetV2Kine(const char *name);
  virtual ~AliAnalysisTaskSEpPbCorrelationsJetV2Kine();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *opt);
  
//virtual void Init();
  virtual void Terminate(Option_t *opt);
//virtual void NotifyRun();

 void SetTrigPtBinning(Int_t nBins, Double_t *limits);
 void SetAssocPtBinning(Int_t nBins, Double_t *limits);
 void SetCentBinning(Int_t nBins, Double_t *limits);
 void FillCorrelationTracksTPCTPC(TObjArray *triggerArray, TObjArray *associateArray, AliTHn *triggerHist, AliTHn *correlationHist, TObjArray *selected_TPC_Pairs, Double_t dcentrality);
 void FillCorrelationTracksTPCTPCMixed(Double_t centrality, TObjArray *triggerArray, TObjArray *associateArray, AliTHn *MixedCorrelation);
 void FillCorrelationTracksTPCTPCFMDA(TObjArray *triggerArray, TObjArray *associateArray, AliTHn *triggerHist, AliTHn *correlationHist);
 void FillCorrelationTracksTPCTPCFMDA_Mixed(Double_t centrality, TObjArray *triggerArray, TObjArray *associateArray, AliTHn *correlationHist);


 
 TObjArray* CloneTrack(TObjArray*selectedTrackArray);
 void SetAssoCut(Float_t mode) {fAsscoptCut=mode;}
 void SetCen1(Double_t d) {fCen1=d;}
 void SetCen2(Double_t d) {fCen2=d;}
 void SetPtMax(Double_t d) {fPtMax = d;}
 void SetPtMin(Double_t d) {fPtMin = d;}
 void SetAnaMode(TString s) {fMode = s;}

 
 Double_t RangePhi(Double_t DPhi);

 private:
 Double_t fCen1;
 Double_t fCen2;
 Double_t fPtMin;
 Double_t fPtMax;
 Double_t fVtxZ;
 Double_t fCentrality;
 Double_t fAsscoptCut;
 TString fMode;

 Int_t fNbinsPtTrig;
 Int_t fNbinsAssocPt;
 //Int_t fNbinsZvtx;

 AliTHn     *fHistTPCTPC_SS; //!
 AliTHn     *fHistTPCTPC_Mixed_SS; //!
 AliTHn     *fHistTPCTrig; //!
 AliTHn     *fHistTPCTPCFMDA; //!
 AliTHn     *fHistTPCTPCFMDA_Mixed; //!
 AliTHn     *fHistTPCTPCFMDATrig; //!

 AliTHn     *fHistTPCFMDA; //!
 AliTHn     *fHistTPCFMDC; //!
 AliTHn     *fHistFMDAFMDC; //!
 AliTHn     *fHistTPCFMDA_Mixed; //!
 AliTHn     *fHistTPCFMDC_Mixed; //!
 AliTHn     *fHistFMDAFMDC_Mixed; //!
 AliTHn     *fHistFMDTrig; //!


 TAxis *fPtTrigAxis;
 TAxis *fPtAssocAxis;
 //TAxis *fZvtxAxis; 

 AliEventPoolManager *fPoolMgr1; //! event pool manager for TPC-TPC
 AliEventPoolManager *fPoolMgr2; //! event pool manager for TPC-TPC-FMD

 
  AliAnalysisTaskSEpPbCorrelationsJetV2Kine(const AliAnalysisTaskSEpPbCorrelationsJetV2Kine&);
  AliAnalysisTaskSEpPbCorrelationsJetV2Kine& operator=(const AliAnalysisTaskSEpPbCorrelationsJetV2Kine&);

  TList *fListOutput; //!
  Int_t fT;
  ClassDef(AliAnalysisTaskSEpPbCorrelationsJetV2Kine, 1);
};

#endif
