#ifndef ALILIGHTNTRACKMCHIST_H
#define ALILIGHTNTRACKMCHIST_H

/*
 * AliLightNTrackMCHist.h
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#include "AliLightNBasePart.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TString.h"
#include "AliLog.h"
#include <iostream>

class AliLightNTrackMCHist {
public:
    AliLightNTrackMCHist();
    AliLightNTrackMCHist(bool contribSplitting,bool DCADist);
    virtual ~AliLightNTrackMCHist();
    void FillMCDCAXYPtBins(AliLightNBasePart::PartOrigin org,
                           int PDGCodeMoth,double pT,double dcaxy);
    void FillMCCorr(double pT){fMCCorrPt->Fill(pT);};
    void FillMCIdent(double pT){fMCIdentPt->Fill(pT);};
    void FillMCGen(double pT){fMCGenPt->Fill(pT);};
    void FillMCGenPrim(double pT){fMCGenPrimPt->Fill(pT);};
    void FillMCCorrPrim(double pT){fMCCorrPrimPt->Fill(pT);};
    void FillMCCorrPrimNoPID(double pT){fMCCorrPrimPnoPID->Fill(pT);};
    void FillMCCont(double pT){fMCContPt->Fill(pT);};
    void FillMCUnkn(double pT){fMCUnknownPt->Fill(pT);};
    void FillMCPrimary(double pT){fMCPrimaryPt->Fill(pT);};
    void FillMCMaterial(double pT){fMCMaterialPt->Fill(pT);};
    void FillMCFeeddown(double pT, double pdg){fMCFeeddownWeakPt->Fill(pT, pdg);};
    void FillStackGen(double p){fMCStackGen->Fill(p);};
    void FillStackGenPrimary(double p){fMCStackGenPrimary->Fill(p);};
    void FillMCpTPCCut(int i, double pTPC){fMCpTPCDist[i]->Fill(pTPC);};
    void FillMCetaCut(int i, double eta){fMCetaDist[i]->Fill(eta);};
    void FillMCphiCut(int i, double phi){fMCphiDist[i]->Fill(phi);};
    void FillMCTPCclsCut(int i, double pT, double nCls){
        fMCTPCCls[i]->Fill(pT, nCls);
    };
    void FillMCDCAxyCut(int i, double pT, double dcaxy){
        fMCDCAxy[i]->Fill(pT, dcaxy);
    };
    void FillMCDCAzCut(int i, double pT, double dcaz){
        fMCDCAz[i]->Fill(pT, dcaz);
    };
    void FillMCTPCCrossedRowCut(int i, double pT, float Crossed){
        fMCTPCCrossedRows[i]->Fill(pT,Crossed);};
    void FillMCTPCRatioCut(int i, double pT, float ratio){
        fMCTPCRatio[i]->Fill(pT,ratio);};
    void FillMCTPCdedx(int i, double mom, double dedx){
        fMCTPCdedx[i]->Fill(mom,dedx);
    };
    void FillMCTOFbeta(int i, double mom, double beta){
        fMCTOFbeta[i]->Fill(mom,beta);
    };
    void FillMCNSigTPC(int i, double mom, double nSigTPC){
        fMCNSigTPC[i]->Fill(mom, nSigTPC);
    };
    void FillMCNSigTOF(int i, double mom, double nSigTOF){
        fMCNSigTOF[i]->Fill(mom, nSigTOF);
    };
    TList *GetHistList() const {return fMCList;};
    TString ClassName() {return "AliLightNTrackMCHist";};
private:
    double fpTmin;                  //!
    double fpTmax;                  //!
    double fpTbins;                 //!
    bool fDoSplitting;              //!
    bool fDoDCAPlots;               //!
    
    TList *fMCList;                 //!
    TList *fMCStackList;            //!
    TList *fMCQAPlots[4];           //!
    TList *fDCAPlots;               //!
    TH1F *fMCCorrPt;                //!
    TH1F *fMCIdentPt;               //!
    TH1F *fMCGenPt;                 //!
    TH1F *fMCGenPrimPt;             //!
    TH1F *fMCCorrPrimPt;            //!
    TH1F *fMCCorrPrimPnoPID;        //!
    TH1F *fMCContPt;                //!
    TH1F *fMCUnknownPt;             //!
    
    TH1F *fMCPrimaryPt;             //!
    TH1F *fMCMaterialPt;            //!
    TH2F *fMCFeeddownWeakPt;        //!
    
    TH1F *fMCStackGen;              //!
    TH1F *fMCStackGenPrimary;       //!
    
    TH2F *fMCPrimDCAXYPtBins;       //!
    TH2F *fMCMaterialDCAXYPtBins;   //!
    TH2F *fMCSecondaryDCAXYPtBins;  //!
    TH2F *fMCSecLambdaDCAXYPtBins;  //!
    TH2F *fMCSecSigmaDCAXYPtBins;   //!
    
    TH1F *fMCpTPCDist[4];           //!
    TH1F *fMCetaDist[4];            //!
    TH1F *fMCphiDist[4];            //!
    TH2F *fMCTPCCls[4];             //!
    TH2F *fMCDCAxy[4];              //!
    TH2F *fMCDCAz[4];               //!
    
    TH2F *fMCTPCCrossedRows[4];     //!
    TH2F *fMCTPCRatio[4];           //!
    
    TH2F *fMCTPCdedx[4];            //!
    TH2F *fMCTOFbeta[4];            //!
    TH2F *fMCNSigTPC[4];            //!
    TH2F *fMCNSigTOF[4];            //!
    ClassDef(AliLightNTrackMCHist,1);
};

#endif /* ALILIGHTNTRACKMCHIST_H */
