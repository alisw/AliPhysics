/*
 * AliLightNTrackMCHist.cxx
 *
 *  Created on: Nov 14,2017
 *      Author: gu74req
 */

#include "AliLightNTrackMCHist.h"
AliLightNTrackMCHist::AliLightNTrackMCHist()
:fpTmin(0.0)
,fpTmax(4.5)
,fpTbins(45)
,fDoSplitting(false)
,fDoDCAPlots(false)
,fMCList(0)
,fMCStackList(0)
,fDCAPlots(0)
,fMCCorrPt(0)
,fMCIdentPt(0)
,fMCGenPt(0)
,fMCGenPrimPt(0)
,fMCCorrPrimPt(0)
,fMCContPt(0)
,fMCUnknownPt(0)
,fMCPrimaryPt(0)
,fMCMaterialPt(0)
,fMCFeeddownWeakPt(0)
,fMCStackGen(0)
,fMCStackGenPrimary(0)
,fMCPrimDCAXYPtBins(0)
,fMCMaterialDCAXYPtBins(0)
,fMCSecondaryDCAXYPtBins(0)
,fMCSecLambdaDCAXYPtBins(0)
,fMCSecSigmaDCAXYPtBins(0)
{
    for (int i=0;i<4;++i) {
        fMCQAPlots[i]=0;
        fMCpTPCDist[i]=0;
        fMCetaDist[i]=0;
        fMCphiDist[i]=0;
        fMCTPCCls[i]=0;
        fMCDCAxy[i]=0;
        fMCDCAz[i]=0;
        fMCTPCCrossedRows[i]=0;
        fMCTPCRatio[i]=0;
        fMCTPCdedx[i]=0;
        fMCTOFbeta[i]=0;
        fMCNSigTPC[i]=0;
        fMCNSigTOF[i]=0;
    }
}

AliLightNTrackMCHist::~AliLightNTrackMCHist()
{
    if (fMCList) {
        delete fMCList;
    }
}

AliLightNTrackMCHist::AliLightNTrackMCHist(bool contribSplitting,bool DCADist)
:fpTmin(0.0)
,fpTmax(4.5)
,fpTbins(45)
,fDoSplitting(contribSplitting)
,fDoDCAPlots(DCADist)
{
    double ptmin=0;
    double ptmax=4.5;
    double ptBins=45;
    double twoDBins=400;
    fMCList=new TList();
    fMCList->SetName("MonteCarlo");
    fMCList->SetOwner();
    
    fMCCorrPt=new TH1F("CorrParPt","Correct Particles Pt",ptBins,ptmin,ptmax);
    fMCCorrPt->Sumw2();
    fMCCorrPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCCorrPt);
    
    fMCIdentPt=new TH1F("IdentPartPt","Ident Particles Pt",ptBins,ptmin,ptmax);
    fMCIdentPt->Sumw2();
    fMCIdentPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCIdentPt);
    
    fMCGenPt=new TH1F("GenPartPt","Gen Particles Pt",ptBins,ptmin,ptmax);
    fMCGenPt->Sumw2();
    fMCGenPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCGenPt);
    
    fMCGenPrimPt=new TH1F("GenPrimPartP","Gen Primary Particles P",ptBins,ptmin,ptmax);
    fMCGenPrimPt->Sumw2();
    fMCGenPrimPt->GetXaxis()->SetTitle("#it{p}");
    fMCList->Add(fMCGenPrimPt);
    
    fMCCorrPrimPt=new TH1F("CorrPrimPartP","Corr Primary Particles P",ptBins,ptmin,ptmax);
    fMCCorrPrimPt->Sumw2();
    fMCCorrPrimPt->GetXaxis()->SetTitle("#it{p}");
    fMCList->Add(fMCCorrPrimPt);
    
    fMCCorrPrimPnoPID=new TH1F("CorrPrimNoPID","Corr Primary Particles NoPID",ptBins,ptmin,ptmax);
    fMCCorrPrimPnoPID->Sumw2();
    fMCCorrPrimPnoPID->GetXaxis()->SetTitle("#it{p}");
    fMCList->Add(fMCCorrPrimPnoPID);
    
    fMCContPt=new TH1F("ContPt","ContPt",ptBins,ptmin,ptmax);
    fMCContPt->Sumw2();
    fMCContPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCContPt);
    
    fMCUnknownPt=new TH1F("UnknPt","UnknPt",ptBins,ptmin,ptmax);
    fMCUnknownPt->Sumw2();
    fMCUnknownPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCUnknownPt);
    
    fMCPrimaryPt=new TH1F("PrimaryPt","PrimaryPt",ptBins,ptmin,ptmax);
    fMCPrimaryPt->Sumw2();
    fMCPrimaryPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCPrimaryPt);
    
    fMCMaterialPt=new TH1F("MatPt","MatPT",ptBins,ptmin,ptmax);
    fMCMaterialPt->Sumw2();
    fMCMaterialPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCMaterialPt);
    
    fMCFeeddownWeakPt=new TH2F("FeeddownPt","Feeddown Pt",ptBins,ptmin,ptmax,
                               213,3121,3334);
    fMCFeeddownWeakPt->Sumw2();
    fMCFeeddownWeakPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCFeeddownWeakPt);
    
    
    //MC stack for generated level
    fMCStackList=new TList();
    fMCStackList->SetName("MCStack_GenLevel");
    fMCStackList->SetOwner();
    fMCList->Add(fMCStackList);

    fMCStackGen= new TH1F("Stack_Generated","Stack_Generated",ptBins,ptmin,ptmax);
    fMCStackGen->Sumw2();
    fMCStackGen->GetXaxis()->SetTitle("p");
    fMCStackList->Add(fMCStackGen);
    
    fMCStackGenPrimary= new TH1F("Stack_Generated_primary","Stack_Generated_primary",ptBins,ptmin,ptmax);
    fMCStackGenPrimary->Sumw2();
    fMCStackGenPrimary->GetXaxis()->SetTitle("p");
    fMCStackList->Add(fMCStackGenPrimary);
    
    if (contribSplitting) {
        TString MCModes[4]={"Primary","Secondary","Material","Contamination"};
        for (int i=0;i<4;++i) {
            fMCQAPlots[i]=new TList();
            fMCQAPlots[i]->SetName(MCModes[i].Data());
            fMCQAPlots[i]->SetOwner();
            fMCList->Add(fMCQAPlots[i]);
            
            TString MCpTPCName=Form("MCpTPCDist%s",MCModes[i].Data());
            TString MCetaName=Form("MCEtaDist%s",MCModes[i].Data());
            TString MCphiName=Form("MCphiDist%s",MCModes[i].Data());
            TString MCTPCName=Form("MCTPCCls%s",MCModes[i].Data());
            TString MCDCAXYName=Form("MCMCDCAXY%s",MCModes[i].Data());
            TString MCDCAZName=Form("MCDCAZ%s",MCModes[i].Data());
            TString MCTPCCRName=Form("MCCrossedRows%s",MCModes[i].Data());
            TString MCTPCratioName=Form("MCTPCRatio%s",MCModes[i].Data());
            TString MCTPCdedxName=Form("MCTPCdedx%s",MCModes[i].Data());
            TString MCTOFbetaName=Form("MCTOFbeta%s",MCModes[i].Data());
            TString MCNSigTPCName=Form("MCNSigTPC%s",MCModes[i].Data());
            TString MCNSigTOFName=Form("MCNSigTOF%s",MCModes[i].Data());
            
            fMCpTPCDist[i]=new TH1F(MCpTPCName.Data(),MCpTPCName.Data(),ptBins,
                                    ptmin,ptmax);
            fMCpTPCDist[i]->Sumw2();
            fMCQAPlots[i]->Add(fMCpTPCDist[i]);
            
            fMCetaDist[i]=new TH1F(MCetaName.Data(),MCetaName.Data(),200,-10.,10.);
            fMCetaDist[i]->Sumw2();
            fMCQAPlots[i]->Add(fMCetaDist[i]);
            
            fMCphiDist[i]=new TH1F(MCphiName.Data(),MCphiName.Data(),200,0.,
                                   2*TMath::Pi());
            fMCphiDist[i]->Sumw2();
            fMCQAPlots[i]->Add(fMCphiDist[i]);
            
            fMCTPCCls[i]=new TH2F(MCTPCName.Data(),MCTPCName.Data(),ptBins,ptmin,
                                  ptmax,100,0,200.);
            fMCTPCCls[i]->Sumw2();
            fMCQAPlots[i]->Add(fMCTPCCls[i]);
            
            fMCDCAxy[i]=new TH2F(MCDCAXYName.Data(),MCDCAXYName.Data(),ptBins,ptmin,
                                 ptmax,2.5*twoDBins,-5.,5.);
            fMCDCAxy[i]->Sumw2();
            fMCQAPlots[i]->Add(fMCDCAxy[i]);
            
            fMCDCAz[i]=new TH2F(MCDCAZName.Data(),MCDCAZName.Data(),ptBins,ptmin,
                                ptmax,2.5*twoDBins,-5.,5.);
            fMCDCAz[i]->Sumw2();
            fMCQAPlots[i]->Add(fMCDCAz[i]);
            
            fMCTPCCrossedRows[i]=new TH2F(MCTPCCRName.Data(),MCTPCCRName.Data(),
                                          ptBins,ptmin,ptmax,100,0,200.);
            fMCTPCCrossedRows[i]->Sumw2();
            fMCQAPlots[i]->Add(fMCTPCCrossedRows[i]);
            
            fMCTPCRatio[i]=new TH2F(MCTPCratioName.Data(),MCTPCratioName.Data(),
                                    ptBins,ptmin,ptmax,100,0.,2.);
            fMCTPCRatio[i]->Sumw2();
            fMCQAPlots[i]->Add(fMCTPCRatio[i]);
            fMCTPCdedx[i]=new TH2F(MCTPCdedxName.Data(),MCTPCdedxName.Data(),
                                   ptBins,ptmin,ptmax,2*twoDBins,0.,400);
            fMCTPCdedx[i]->Sumw2();
            fMCQAPlots[i]->Add(fMCTPCdedx[i]);
            fMCTOFbeta[i]=new TH2F(MCTOFbetaName.Data(),MCTOFbetaName.Data(),
                                   ptBins,ptmin,ptmax,3.5*twoDBins,0.4,1.1);
            fMCTOFbeta[i]->Sumw2();
            fMCQAPlots[i]->Add(fMCTOFbeta[i]);
            fMCNSigTPC[i]=new TH2F(MCNSigTPCName.Data(),MCNSigTPCName.Data(),
                                   ptBins,ptmin,ptmax,twoDBins,-20.,20.);
            fMCNSigTPC[i]->Sumw2();
            fMCQAPlots[i]->Add(fMCNSigTPC[i]);
            fMCNSigTOF[i]=new TH2F(MCNSigTOFName.Data(),MCNSigTOFName.Data(),
                                   ptBins,ptmin,ptmax,twoDBins,-20.,20.);
            fMCNSigTOF[i]->Sumw2();
            fMCQAPlots[i]->Add(fMCNSigTOF[i]);
        }
    } else {
        for (int i=0;i<4;++i) {
            fMCQAPlots[i]=0;
            fMCpTPCDist[i]=0;
            fMCetaDist[i]=0;
            fMCphiDist[i]=0;
            fMCTPCCls[i]=0;
            fMCDCAxy[i]=0;
            fMCDCAz[i]=0;
            fMCTPCCrossedRows[i]=0;
            fMCTPCRatio[i]=0;
            fMCTPCdedx[i]=0;
            fMCTOFbeta[i]=0;
            fMCNSigTPC[i]=0;
            fMCNSigTOF[i]=0;
        }
    }
    
    if (DCADist) {
        fDCAPlots=new TList();
        fDCAPlots->SetName("DCAPtBinning");
        fDCAPlots->SetOwner();
        fMCList->Add(fDCAPlots);
        TString MCPridcaPtBinName=Form("DCAPtBinningPri");
        TString MCMatdcaPtBinName=Form("DCAPtBinningMat");
        TString MCSecdcaPtBinName=Form("DCAPtBinningSec");
        TString MCSecLamdcaPtBinName=Form("DCAPtBinningSecLam");
        TString MCSecSigdcaPtBinName=Form("DCAPtBinningSecSig");
        
        fMCPrimDCAXYPtBins=new TH2F(MCPridcaPtBinName.Data(),
                                    MCPridcaPtBinName.Data(),fpTbins,fpTmin,fpTmax
                                    ,1000,-5,5);
        fMCPrimDCAXYPtBins->Sumw2();
        fMCPrimDCAXYPtBins->GetXaxis()->SetTitle("p_{T}");
        fMCPrimDCAXYPtBins->GetYaxis()->SetTitle("dca_{XY}");
        fDCAPlots->Add(fMCPrimDCAXYPtBins);
        
        fMCMaterialDCAXYPtBins=new TH2F(MCMatdcaPtBinName.Data(),
                                        MCMatdcaPtBinName.Data(),fpTbins,fpTmin,
                                        fpTmax,1000,-5,5);
        fMCMaterialDCAXYPtBins->Sumw2();
        fMCMaterialDCAXYPtBins->GetXaxis()->SetTitle("p_{T}");
        fMCMaterialDCAXYPtBins->GetYaxis()->SetTitle("dca_{XY}");
        fDCAPlots->Add(fMCMaterialDCAXYPtBins);
        
        fMCSecondaryDCAXYPtBins=new TH2F(MCSecdcaPtBinName.Data(),
                                         MCSecdcaPtBinName.Data(),fpTbins,fpTmin,
                                         fpTmax,1000,-5,5);
        fMCSecondaryDCAXYPtBins->Sumw2();
        fMCSecondaryDCAXYPtBins->GetXaxis()->SetTitle("p_{T}");
        fMCSecondaryDCAXYPtBins->GetYaxis()->SetTitle("dca_{XY}");
        fDCAPlots->Add(fMCSecondaryDCAXYPtBins);
        
        fMCSecLambdaDCAXYPtBins=new TH2F(MCSecLamdcaPtBinName.Data(),
                                         MCSecLamdcaPtBinName.Data(),fpTbins,
                                         fpTmin,fpTmax,1000,-5,5);
        fMCSecLambdaDCAXYPtBins->Sumw2();
        fMCSecLambdaDCAXYPtBins->GetXaxis()->SetTitle("p_{T}");
        fMCSecLambdaDCAXYPtBins->GetYaxis()->SetTitle("dca_{XY}");
        fDCAPlots->Add(fMCSecLambdaDCAXYPtBins);
        
        fMCSecSigmaDCAXYPtBins=new TH2F(MCSecSigdcaPtBinName.Data(),
                                        MCSecSigdcaPtBinName.Data(),
                                        fpTbins,fpTmin,fpTmax,1000,-5,5);
        fMCSecSigmaDCAXYPtBins->Sumw2();
        fMCSecSigmaDCAXYPtBins->GetXaxis()->SetTitle("p_{T}");
        fMCSecSigmaDCAXYPtBins->GetYaxis()->SetTitle("dca_{XY}");
        fDCAPlots->Add(fMCSecSigmaDCAXYPtBins);
    } else {
        fDCAPlots=0;
        fMCPrimDCAXYPtBins=0;
        fMCMaterialDCAXYPtBins=0;
        fMCSecondaryDCAXYPtBins=0;
        fMCSecLambdaDCAXYPtBins=0;
        fMCSecSigmaDCAXYPtBins=0;
    }
}
void AliLightNTrackMCHist::FillMCDCAXYPtBins(
                                             AliLightNBasePart::PartOrigin org,int PDGCodeMoth,
                                             double pT,double dcaxy)
{
    if (!fDoDCAPlots) {
        AliFatal("FullBooking not set for SPCutHistograms! Cannot use this method");
    }
    if (org==AliLightNBasePart::kPhysPrimary) {
        fMCPrimDCAXYPtBins->Fill(pT,dcaxy);
    } else if (org==AliLightNBasePart::kWeak) {
        fMCSecondaryDCAXYPtBins->Fill(pT,dcaxy);
        if (TMath::Abs(PDGCodeMoth)==3222) {
            fMCSecSigmaDCAXYPtBins->Fill(pT,dcaxy);
        } else if (TMath::Abs(PDGCodeMoth)==3122) {
            fMCSecLambdaDCAXYPtBins->Fill(pT,dcaxy);
        } else {
            TString ErrHistSP=Form("Feeddown for %d not implemented",PDGCodeMoth);
            AliWarning(ErrHistSP.Data());
        }
    } else if (org==AliLightNBasePart::kMaterial) {
        fMCMaterialDCAXYPtBins->Fill(pT,dcaxy);
    } else {
        //For all particles with other origin i.e. missidentified by the analysis
        return;
    }
    return;
}
