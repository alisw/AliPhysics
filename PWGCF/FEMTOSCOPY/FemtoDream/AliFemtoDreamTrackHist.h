/*
 * AliFemtoDreamTrackHist.h
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMTRACKHIST_H_
#define ALIFEMTODREAMTRACKHIST_H_
#include "AliPIDResponse.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TProfile.h"

class AliFemtoDreamTrackHist {
 public:
  AliFemtoDreamTrackHist();
  AliFemtoDreamTrackHist(bool DCADist,bool CombSig);
  virtual ~AliFemtoDreamTrackHist();
  void FillConfig(int iBin, double val){fConfig->Fill(iBin, val);};
  void FillTrackCounter(int iBin){fCutCounter->Fill(iBin);};
  void FillpTCut(int i, double pT){fpTDist[i]->Fill(pT);};
  void FillpTPCCut(int i, double pTPC){fpTPCDist[i]->Fill(pTPC);};
  void FilletaCut(int i, double eta){fetaDist[i]->Fill(eta);};
  void FillphiCut(int i, double phi){fphiDist[i]->Fill(phi);};
  void FillTPCclsCut(int i, double nCls){fTPCCls[i]->Fill(nCls);};
  void FillDCAxyCut(int i, double pT, double dcaxy){
    fDCAxy[i]->Fill(pT, dcaxy);};
  void FillDCAzCut(int i, double pT, double dcaz){fDCAz[i]->Fill(pT, dcaz);};
  void FillTPCCrossedRowCut(int i, float Crossed){
    fTPCCrossedRows[i]->Fill(Crossed);
  };
  void FillTPCRatioCut(int i, float ratio){fTPCRatio[i]->Fill(ratio);};
  void FillTPCClsS(int i, double TPCClsS){fTPCClsS[i]->Fill(TPCClsS);};
  void FillHasSharedClsITS(int i,int layer,int yesno){
    fShrdClsITS[i]->Fill(layer,yesno);};
  void FillTPCdedx(int i, double mom, double dedx){
    fTPCdedx[i]->Fill(mom,dedx);
  };
  void FillTOFbeta(int i, double mom, double beta){
    fTOFbeta[i]->Fill(mom,beta);
  };
  void FillNSigTPC(int i, double mom, double nSigTPC){
    fNSigTPC[i]->Fill(mom, nSigTPC);
  };
  void FillNSigTOF(int i, double mom, double nSigTOF){
    fNSigTOF[i]->Fill(mom, nSigTOF);
  };
  void FillTPCStatus(int i, AliPIDResponse::EDetPidStatus statusTPC){
    fTPCStatus[i]->Fill(statusTPC);
  };
  void FillTOFStatus(int i, AliPIDResponse::EDetPidStatus statusTOF){
    fTOFStatus[i]->Fill(statusTOF);
  };
  void FillNSigComb(double pT, double nSigTPC, double nSigTOF);
  void FillDCAXYPtBins(double pT, double dcaxy);
  void FillTPCClsCPileUp(int i,int iCrit,double TPCClsC){
    fTPCClsCPiluUp[i]->Fill(iCrit,TPCClsC);
  }
  void FillITSSharedPileUp(int i,int iCrit,int yesno){
    fITShrdClsPileUp[i]->Fill(iCrit,yesno);
  }
  void SetName(TString name){fHistList->SetName(name.Data());};
  TList *GetHistList(){return fHistList;};
 private:
  TList *fHistList;         //!
  TList *fTrackCutQA[2];    //!
  TProfile *fConfig;        //!
  TH1F *fCutCounter;        //!
  TH1F *fpTDist[2];         //!
  TH1F *fpTPCDist[2];       //!
  TH1F *fetaDist[2];        //!
  TH1F *fphiDist[2];        //!
  TH1F *fTPCCls[2];         //!
  TH2F *fShrdClsITS[2];     //!
  TH2F *fDCAxy[2];          //!
  TH2F *fDCAz[2];           //!
  TH2F *fDCAXYPtBins;       //!
  TH1F *fTPCCrossedRows[2]; //!
  TH1F *fTPCRatio[2];       //!
  TH1F *fTPCClsS[2];        //!
  TH2F *fTPCdedx[2];        //!
  TH2F *fTOFbeta[2];        //!
  TH2F *fNSigTPC[2];        //!
  TH2F *fNSigTOF[2];        //!
  TH1F *fTPCStatus[2];      //!
  TH1F *fTOFStatus[2];      //!
  TH3F *fNSigCom;           //!
  TH2F *fTPCClsCPiluUp[2];  //!
  TH2F *fITShrdClsPileUp[2];//!
  ClassDef(AliFemtoDreamTrackHist,1);
};

#endif /* ALIFEMTODREAMTRACKHIST_H_ */
