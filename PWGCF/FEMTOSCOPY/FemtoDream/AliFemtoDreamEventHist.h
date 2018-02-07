/*
 * AliFemtoDreamEventHist.h
 *
 *  Created on: Nov 22, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMEVENTHIST_H_
#define ALIFEMTODREAMEVENTHIST_H_
#include "Rtypes.h"
#include "TList.h"
#include "TH1F.h"
#include "TProfile.h"
class AliFemtoDreamEventHist {
 public:
  AliFemtoDreamEventHist();
  virtual ~AliFemtoDreamEventHist();
  void FillEvtCounter(int iBin){fEvtCounter->Fill(iBin);};
  void FillCuts(int iBin,double val){fCutConfig->Fill(iBin,val);};
  void FillEvtNCont(int i, double val){fEvtNCont[i]->Fill(val);};
  void FillEvtVtxX(int i, double val){fEvtVtxX[i]->Fill(val);};
  void FillEvtVtxY(int i, double val){fEvtVtxY[i]->Fill(val);};
  void FillEvtVtxZ(int i, double val){fEvtVtxZ[i]->Fill(val);};
  void FillMultSPD(int i, double val){fMultDistSPD[i]->Fill(val);};
  void FillMultV0A(int i, double val){fMultDistV0A[i]->Fill(val);};
  void FillMultV0C(int i, double val){fMultDistV0C[i]->Fill(val);};
  void FillMultRef08(int i, double val){fMultDistRef08[i]->Fill(val);};
  TList *GetHistList() {return fEventCutList;};
 private:
  TList *fEventCutList;     //!
  TList *fEvtCutQA[2];      //!
  TH1F *fEvtCounter;        //!
  TProfile *fCutConfig;     //!
  TH1F *fEvtNCont[2];       //!
  TH1F *fEvtVtxX[2];        //!
  TH1F *fEvtVtxY[2];        //!
  TH1F *fEvtVtxZ[2];        //!
  TH1F *fMultDistSPD[2];    //!
  TH1F *fMultDistV0A[2];    //!
  TH1F *fMultDistV0C[2];    //!
  TH1F *fMultDistRef08[2];  //!
  ClassDef(AliFemtoDreamEventHist,1)
};

#endif /* ALIFEMTODREAMEVENTHIST_H_ */
