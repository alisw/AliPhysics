#ifndef ALILIGHTNEVENTHIST_H
#define ALILIGHTNEVENTHIST_H

/*
 * AliLightNEventHist.h
 *
 *  Created on: Nov 22, 2017
 *      Author: gu74req
 */

#include "Rtypes.h"
#include "TList.h"
#include "TH1F.h"
#include "TProfile.h"
class AliLightNEventHist {
public:
    AliLightNEventHist();
    virtual ~AliLightNEventHist();
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
    void FillV0Mpercentile(double val){fV0Mpercentile->Fill(val);};
    void FillV0MpercentileHM(double val){fV0MpercentileHM->Fill(val);};
    TList *GetHistList() {return fEventCutList;};
    void SetName(TString name){fEventCutList->SetName(name.Data());};
private:
    TList *fEventCutList;     //!
    TList *fEvtCutQA[2];      //!
    TH1F *fEvtCounter;        //!
    TProfile *fCutConfig;     //!
    TH1F *fV0Mpercentile;      //!
    TH1F *fV0MpercentileHM;      //!
    TH1F *fEvtNCont[2];       //!
    TH1F *fEvtVtxX[2];        //!
    TH1F *fEvtVtxY[2];        //!
    TH1F *fEvtVtxZ[2];        //!
    TH1F *fMultDistSPD[2];    //!
    TH1F *fMultDistV0A[2];    //!
    TH1F *fMultDistV0C[2];    //!
    TH1F *fMultDistRef08[2];  //!
    ClassDef(AliLightNEventHist,1)
};

#endif /* ALILIGHTNEVENTHIST_H */
