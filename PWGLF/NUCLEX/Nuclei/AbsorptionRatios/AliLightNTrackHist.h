#ifndef ALILIGHTNTRACKHIST_H
#define ALILIGHTNTRACKHIST_H

/*
 * AliLightNTrackHist.h
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#include "AliPIDResponse.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TProfile.h"
#include "TMath.h"
class AliLightNTrackHist {
public:
    AliLightNTrackHist();
    AliLightNTrackHist(bool DCADist,bool CombSig,bool PlotPmass2dca3D);
    virtual ~AliLightNTrackHist();
    void FillConfig(int iBin, double val){fConfig->Fill(iBin, val);};
    void FillTrackCounter(int iBin){fCutCounter->Fill(iBin);};
    void FillpCut(int i, double p){fpDist[i]->Fill(p);};
    void FillpTPCCut(int i, double pTPC){fpTPCDist[i]->Fill(pTPC);};
    void FillpDiff_p_pTPC(int i, double p, double difference){fDiff_p_pTPC[i]->Fill(p,difference);};
    void FilletaCut(int i, double eta){fetaDist[i]->Fill(eta);};
    void FillrapidityCut(int i, double rapidity){fRapidityDist[i]->Fill(rapidity);};
    void FillphiCut(int i, double phi){fphiDist[i]->Fill(phi);};
    void FillTPCclsCut(int i, double nCls){fTPCCls[i]->Fill(nCls);};
    void FillTPCclsHighPur(int i, double nCls){fTPCClsHighPur[i]->Fill(nCls);};
    void FillDCAxyCut(int i, double p, double dcaxy){
        fDCAxy[i]->Fill(p, dcaxy);};
    void FillDCAzCut(int i, double p, double dcaz){fDCAz[i]->Fill(p, dcaz);};
    void FillMass2sq(int i, double p, double mass2sq){fMass2sqHist[i]->Fill(p,mass2sq);};
    void FillEtaPhiTPConlyPID(int i, double eta, double phi){fEtaPhiTPConlyPIDHist[i]->Fill(eta,phi);};
    void FillEtaPhiTPCTOFPID(int i, double eta, double phi){fEtaPhiTPCTOFPIDHist[i]->Fill(eta,phi);};
    void FillTPCCrossedRowCut(int i, float Crossed){fTPCCrossedRows[i]->Fill(Crossed);};
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
    void FillITSdedx(int i, double mom, double dedx){
        fITSdedx[i]->Fill(mom,dedx);
    };
    void FillNSigTPC(int i, double mom, double nSigTPC){
        fNSigTPC[i]->Fill(mom, nSigTPC);
    };
    void FillNSigTOF(int i, double mom, double nSigTOF){
        fNSigTOF[i]->Fill(mom, nSigTOF);
    };
    void FillNSigITS(int i, double mom, double nSigITS){
        fNSigITS[i]->Fill(mom, nSigITS);
    };
    void FillTPCStatus(int i, AliPIDResponse::EDetPidStatus statusTPC){
        fTPCStatus[i]->Fill(statusTPC);
    };
    void FillTOFStatus(int i, AliPIDResponse::EDetPidStatus statusTOF){
        fTOFStatus[i]->Fill(statusTOF);
    };
    void FillITSStatus(int i, AliPIDResponse::EDetPidStatus statusITS){
        fITSStatus[i]->Fill(statusITS);
    };
    void FillNSigComb(double p, double nSigTPC, double nSigTOF);
    void FillDCAXYPBins(double p, double dcaxy);
    void FillP_mass2_DCAxy(double p, double mass2,double dcaxyCut);
    void FillTPCClsCPileUp(int i,int iCrit,double TPCClsC){
        fTPCClsCPiluUp[i]->Fill(iCrit,TPCClsC);
    }
    void FillITSSharedPileUp(int i,int iCrit,int yesno){
        fITShrdClsPileUp[i]->Fill(iCrit,yesno);
    }
    void FillTrackChi2perNDF(int i ,double Chi2perNDF){
        fTrackChi2perNDF[i]->Fill(Chi2perNDF);
    }
    void SetName(TString name){fHistList->SetName(name.Data());};
    TList *GetHistList(){return fHistList;};
private:
    TList *fHistList;         //!
    TList *fTrackCutQA[2];    //!
    TProfile *fConfig;        //!
    TH1F *fCutCounter;        //!
    TH1F *fpDist[2];         //!
    TH1F *fpTPCDist[2];       //!
    TH2F *fDiff_p_pTPC[2];        //!
    TH1F *fetaDist[2];        //!
    TH1F *fphiDist[2];        //!
    TH1F *fTPCCls[2];         //!
    TH1F *fTPCClsHighPur[2];  //!
    TH1F *fRapidityDist[2];   //!
    TH2F *fShrdClsITS[2];     //!
    TH2F *fDCAxy[2];          //!
    TH2F *fDCAz[2];           //!
    TH2F *fMass2sqHist[2];      //!
    TH2F *fEtaPhiTPConlyPIDHist[2];    //!
    TH2F *fEtaPhiTPCTOFPIDHist[2];    //!
    TH2F *fDCAXYPBins;       //!
    TH1F *fTPCCrossedRows[2]; //!
    TH1F *fTPCRatio[2];       //!
    TH1F *fTPCClsS[2];        //!
    TH2F *fTPCdedx[2];        //!
    TH2F *fITSdedx[2];        //!
    TH2F *fTOFbeta[2];        //!
    TH2F *fNSigTPC[2];        //!
    TH2F *fNSigTOF[2];        //!
    TH2F *fNSigITS[2];        //!
    TH1F *fTPCStatus[2];      //!
    TH1F *fTOFStatus[2];      //!
    TH1F *fITSStatus[2];      //!
    TH3F *fNSigCom;           //!
    TH3F *fP_mass2_DCAxyHist;    //!
    TH2F *fTPCClsCPiluUp[2];  //!
    TH2F *fITShrdClsPileUp[2];//!
    TH1F *fTrackChi2perNDF[2]; //!
    ClassDef(AliLightNTrackHist,1);
};

#endif /* ALILIGHTNTRACKHIST_H */
