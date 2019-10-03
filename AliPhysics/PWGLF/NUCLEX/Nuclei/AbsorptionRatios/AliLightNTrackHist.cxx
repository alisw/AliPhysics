/*
 * AliLightNTrackHist.cxx
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#include "AliLightNTrackHist.h"

ClassImp(AliLightNTrackHist)
AliLightNTrackHist::AliLightNTrackHist()
:fHistList(0)
,fConfig(0)
,fCutCounter(0)
,fDCAXYPBins(0)
,fNSigCom(0)
,fP_mass2_DCAxyHist(0)
{
    for (int i=0;i<2;++i) {
        fTrackCutQA[i]=0;
        fpDist[i]=0;
        fpTPCDist[i]=0;
        fDiff_p_pTPC[i]=0;
        fetaDist[i]=0;
        fphiDist[i]=0;
        fTPCCls[i]=0;
        fTPCClsS[i]=0;
        fShrdClsITS[i]=0;
        fDCAxy[i]=0;
        fDCAz[i]=0;
        fMass2sqHist[i]=0;
        fTPCCrossedRows[i]=0;
        fTPCRatio[i]=0;
        fTPCdedx[i]=0;
        fTOFbeta[i]=0;
        fNSigTPC[i]=0;
        fNSigTOF[i]=0;
        fTPCStatus[i]=0;
        fTOFStatus[i]=0;
        fTPCClsCPiluUp[i]=0;
        fITShrdClsPileUp[i]=0;
        fTrackChi2perNDF[i]=0;
    }
}
AliLightNTrackHist::AliLightNTrackHist(bool DCADist,bool CombSig,bool PlotPmass2dca3D) {
    TString sName[2] = {"before", "after"};
    double pmin = 0;
    double pmax = 4.5;
    int pBins = 90;
    int twoDBins = 400;
    fHistList = new TList();
    fHistList->SetName("TrackCuts");
    fHistList->SetOwner();
    
    fConfig = new TProfile("TrackCutConfig", "Track Cut Config", 22, 0, 22);
    fConfig->SetStats(0);
    fConfig->GetXaxis()->SetBinLabel(1, "p_{min}");
    fConfig->GetXaxis()->SetBinLabel(2, "p_{max}");
    fConfig->GetXaxis()->SetBinLabel(3, "#eta_{min}");
    fConfig->GetXaxis()->SetBinLabel(4, "#eta_{max}");
    fConfig->GetXaxis()->SetBinLabel(5, "Charge");
    fConfig->GetXaxis()->SetBinLabel(6, "ClsTPC_{min}");
    fConfig->GetXaxis()->SetBinLabel(7, "Filter Bit");
    fConfig->GetXaxis()->SetBinLabel(8, "Is MC");
    fConfig->GetXaxis()->SetBinLabel(9, "DCA_{xy}");
    fConfig->GetXaxis()->SetBinLabel(10, "DCA_{z}");
    fConfig->GetXaxis()->SetBinLabel(11, "Shared Cls");
    fConfig->GetXaxis()->SetBinLabel(12, "TPC Crossed Rows");
    fConfig->GetXaxis()->SetBinLabel(13, "PID Mom Thresh");
    fConfig->GetXaxis()->SetBinLabel(14, "PID nsig");
    fConfig->GetXaxis()->SetBinLabel(15, "Reject Pions");
    fConfig->GetXaxis()->SetBinLabel(16, "Smallest Sig");
    fConfig->GetXaxis()->SetBinLabel(17, "ITS Hit");
    fConfig->GetXaxis()->SetBinLabel(18, "TOF Timing");
    fConfig->GetXaxis()->SetBinLabel(19, "Pile Up Rej");
    fConfig->GetXaxis()->SetBinLabel(20, "TPC Refit");
    fConfig->GetXaxis()->SetBinLabel(21, "RapiMin");
    fConfig->GetXaxis()->SetBinLabel(22, "RapiMax");
    
    fHistList->Add(fConfig);
    
    fCutCounter = new TH1F("CutCounter", "Cut Counter", 23, 0, 23);
    fCutCounter->GetXaxis()->SetBinLabel(1, "Input");
    fCutCounter->GetXaxis()->SetBinLabel(2, "Filter Bit");
    fCutCounter->GetXaxis()->SetBinLabel(3, "#it{p} Cut");
    fCutCounter->GetXaxis()->SetBinLabel(4, "#eta Cut");
    fCutCounter->GetXaxis()->SetBinLabel(5, "Charge");
    fCutCounter->GetXaxis()->SetBinLabel(6, "PileUpITS");
    fCutCounter->GetXaxis()->SetBinLabel(7, "PileUpTOF");
    fCutCounter->GetXaxis()->SetBinLabel(8, "PileUp");
    fCutCounter->GetXaxis()->SetBinLabel(9, "nClsTPC");
    fCutCounter->GetXaxis()->SetBinLabel(10, "Shared Cls");
    fCutCounter->GetXaxis()->SetBinLabel(11, "TPC Refit");
    fCutCounter->GetXaxis()->SetBinLabel(12, "TPC Crossed Rows");
    fCutCounter->GetXaxis()->SetBinLabel(13, "TPC Row Ratio");
    fCutCounter->GetXaxis()->SetBinLabel(14, "Rapidity Cut");
    fCutCounter->GetXaxis()->SetBinLabel(15, "TPC OK");
    fCutCounter->GetXaxis()->SetBinLabel(16, "TPC PID");
    fCutCounter->GetXaxis()->SetBinLabel(17, "Reject Pions");
    fCutCounter->GetXaxis()->SetBinLabel(18, "TPC TOF OK");
    fCutCounter->GetXaxis()->SetBinLabel(19, "TPC TOF PID");
    fCutCounter->GetXaxis()->SetBinLabel(20, "Smallest Sig");
    fCutCounter->GetXaxis()->SetBinLabel(21, "Passes PID");
    fCutCounter->GetXaxis()->SetBinLabel(22, "DCA_{Z}");
    fCutCounter->GetXaxis()->SetBinLabel(23, "DCA_{XY}");
    
    fHistList->Add(fCutCounter);
    
    for (int i = 0;i<2;++i) {
        fTrackCutQA[i] = new TList();
        fTrackCutQA[i]->SetName(sName[i].Data());
        fTrackCutQA[i]->SetOwner();
        fHistList->Add(fTrackCutQA[i]);
        
        TString pName=Form("pDist_%s", sName[i].Data());
        fpDist[i] = new TH1F(pName.Data(), pName.Data(), 45, pmin, pmax);
        fpDist[i]->Sumw2();
        fpDist[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fTrackCutQA[i]->Add(fpDist[i]);
        
        TString pTPCName=Form("pTPCDist_%s", sName[i].Data());
        fpTPCDist[i] = new TH1F(pTPCName.Data(), pTPCName.Data(), 45, pmin, pmax);
        fpTPCDist[i]->Sumw2();
        fpTPCDist[i]->GetXaxis()->SetTitle("p_{TPC}");
        fTrackCutQA[i]->Add(fpTPCDist[i]);
        
        TString Diff_p_pTPCName = Form("Diff_p_pTPC_%s", sName[i].Data());
        fDiff_p_pTPC[i] = new TH2F(Diff_p_pTPCName.Data(), Diff_p_pTPCName.Data(), 45, pmin, pmax, 100, -0.3, 0.3);
        fDiff_p_pTPC[i]->Sumw2();
        fDiff_p_pTPC[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fDiff_p_pTPC[i]->GetYaxis()->SetTitle("#it{p} - #it{p}_{TPC} (GeV/#it{c})");
        fTrackCutQA[i]->Add(fDiff_p_pTPC[i]);
        
        TString etaName = Form("EtaDist_%s", sName[i].Data());
        fetaDist[i] = new TH1F(etaName.Data(), etaName.Data(), 200, -10., 10.);
        fetaDist[i]->Sumw2();
        fetaDist[i]->GetXaxis()->SetTitle("#eta");
        fTrackCutQA[i]->Add(fetaDist[i]);
        
        TString phiName = Form("phiDist_%s", sName[i].Data());
        fphiDist[i] = new TH1F(phiName.Data(), phiName.Data(), 200, 0., 2*TMath::Pi());
        fphiDist[i]->Sumw2();
        fphiDist[i]->GetXaxis()->SetTitle("#it{p}hi");
        fTrackCutQA[i]->Add(fphiDist[i]);
        
        TString RapidityName = Form("RapidityDist_%s", sName[i].Data());
        fRapidityDist[i] = new TH1F(RapidityName.Data(), RapidityName.Data(), 100, -5., 5.);
        fRapidityDist[i]->Sumw2();
        fRapidityDist[i]->GetXaxis()->SetTitle("#it{y}_{CMS}");
        fTrackCutQA[i]->Add(fRapidityDist[i]);
        
        TString TPCNClsName = Form("TPCCls_%s", sName[i].Data());
        fTPCCls[i] = new TH1F(TPCNClsName.Data(), TPCNClsName.Data(), 100, 0, 200.);
        fTPCCls[i]->Sumw2();
        fTPCCls[i]->GetXaxis()->SetTitle("# cls TPC");
        fTrackCutQA[i]->Add(fTPCCls[i]);
        
        TString TPCNClsHighPurName = Form("TPCCls_with_high_purity_%s", sName[i].Data());
        fTPCClsHighPur[i] = new TH1F(TPCNClsHighPurName.Data(), TPCNClsHighPurName.Data(), 100, 0, 200.);
        fTPCClsHighPur[i]->Sumw2();
        fTPCClsHighPur[i]->GetXaxis()->SetTitle("# cls TPC");
        fTrackCutQA[i]->Add(fTPCClsHighPur[i]);
        
        
        TString DCAXYName = Form("DCAXY_%s", sName[i].Data());
        fDCAxy[i] = new TH2F(DCAXYName.Data(), DCAXYName.Data(), pBins, pmin, pmax, 2.5*twoDBins, -5., 5.);
        fDCAxy[i]->Sumw2();
        fDCAxy[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fDCAxy[i]->GetYaxis()->SetTitle("DCA_{xy}");
        
        fTrackCutQA[i]->Add(fDCAxy[i]);
        
        TString DCAZName = Form("DCAZ_%s", sName[i].Data());
        fDCAz[i] = new TH2F(DCAZName.Data(), DCAZName.Data(), pBins, pmin, pmax, 2.5*twoDBins, -5., 5.);
        fDCAz[i]->Sumw2();
        fDCAz[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fDCAz[i]->GetYaxis()->SetTitle("DCA_{z}");
        fTrackCutQA[i]->Add(fDCAz[i]);
        
        TString Mass2sqName = Form("Mass2sq_%s", sName[i].Data());
        fMass2sqHist[i] = new TH2F(Mass2sqName.Data(), Mass2sqName.Data(), 45, pmin, pmax, 5*twoDBins, -10., 10.);
        fMass2sqHist[i]->Sumw2();
        fMass2sqHist[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fMass2sqHist[i]->GetYaxis()->SetTitle("m^{2}");
        fTrackCutQA[i]->Add(fMass2sqHist[i]);
        
        TString EtaPhiTPConlyPIDName = Form("EtaPhiTPConlyPID_%s", sName[i].Data());
        fEtaPhiTPConlyPIDHist[i] = new TH2F(EtaPhiTPConlyPIDName.Data(), EtaPhiTPConlyPIDName.Data(), 100, -2, 2, 200, 0., 2*TMath::Pi());
        fEtaPhiTPConlyPIDHist[i]->Sumw2();
        fEtaPhiTPConlyPIDHist[i]->GetXaxis()->SetTitle("#eta");
        fEtaPhiTPConlyPIDHist[i]->GetYaxis()->SetTitle("#phi");
        fTrackCutQA[i]->Add(fEtaPhiTPConlyPIDHist[i]);
        
        TString EtaPhiTPCTOFPIDName = Form("EtaPhiTPCTOFPID_%s", sName[i].Data());
        fEtaPhiTPCTOFPIDHist[i] = new TH2F(EtaPhiTPCTOFPIDName.Data(), EtaPhiTPCTOFPIDName.Data(), 100, -2, 2, 200, 0., 2*TMath::Pi());
        fEtaPhiTPCTOFPIDHist[i]->Sumw2();
        fEtaPhiTPCTOFPIDHist[i]->GetXaxis()->SetTitle("#eta");
        fEtaPhiTPCTOFPIDHist[i]->GetYaxis()->SetTitle("#phi");
        fTrackCutQA[i]->Add(fEtaPhiTPCTOFPIDHist[i]);
        
        TString TPCCRName = Form("CrossedRows_%s", sName[i].Data());
        fTPCCrossedRows[i] = new TH1F(TPCCRName.Data(), TPCCRName.Data(), 100, 0, 200.);
        fTPCCrossedRows[i]->Sumw2();
        fTPCCrossedRows[i]->GetXaxis()->SetTitle("# of crossed Rows");
        fTrackCutQA[i]->Add(fTPCCrossedRows[i]);
        
        TString TPCratioName = Form("TPCRatio_%s", sName[i].Data());
        fTPCRatio[i] = new TH1F(TPCratioName.Data(), TPCratioName.Data(), 100, 0., 2.);
        fTPCRatio[i]->Sumw2();
        fTPCRatio[i]->GetXaxis()->SetTitle("Ratio");
        fTrackCutQA[i]->Add(fTPCRatio[i]);
        
        TString TPCClsSName=Form("TPCClsS_%s",sName[i].Data());
        fTPCClsS[i]=new TH1F(TPCClsSName.Data(),TPCClsSName.Data(),200,0,200);
        fTPCClsS[i]->Sumw2();
        fTPCClsS[i]->GetXaxis()->SetTitle("TPC Cls S");
        fTrackCutQA[i]->Add(fTPCClsS[i]);
        
        TString ShrdClsITSName=Form("SharedClsITS_%s",sName[i].Data());
        fShrdClsITS[i]=new TH2F(ShrdClsITSName.Data(),ShrdClsITSName.Data(),
                                6,1,7,2,0,2);
        fShrdClsITS[i]->Sumw2();
        fShrdClsITS[i]->GetXaxis()->SetTitle("ITS Layer");
        fShrdClsITS[i]->GetYaxis()->SetBinLabel(1,"Has shared Cls");
        fShrdClsITS[i]->GetYaxis()->SetBinLabel(2,"Has no Shard Cls");
        fTrackCutQA[i]->Add(fShrdClsITS[i]);
        
        TString TPCdedxName = Form("TPCdedx_%s", sName[i].Data());
        fTPCdedx[i] = new TH2F(TPCdedxName.Data(), TPCdedxName.Data(), pBins, pmin, pmax, 4*twoDBins, 0., 800);
        fTPCdedx[i]->Sumw2();
        fTPCdedx[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fTPCdedx[i]->GetYaxis()->SetTitle("dEdx (arb. Units)");
        
        fTrackCutQA[i]->Add(fTPCdedx[i]);
        
        TString TOFbetaName = Form("TOFbeta_%s", sName[i].Data());
        fTOFbeta[i] = new TH2F(TOFbetaName.Data(), TOFbetaName.Data(), pBins, pmin, pmax, 3.5*twoDBins, 0.3, 1.1);
        fTOFbeta[i]->Sumw2();
        fTOFbeta[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fTOFbeta[i]->GetYaxis()->SetTitle("#beta_{TOF}");
        
        fTrackCutQA[i]->Add(fTOFbeta[i]);
        
        TString ITSdedxName = Form("ITSdedx_%s", sName[i].Data());
        fITSdedx[i] = new TH2F(ITSdedxName.Data(), ITSdedxName.Data(), pBins, pmin, pmax, 4*twoDBins, 0., 800);
        fITSdedx[i]->Sumw2();
        fITSdedx[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fITSdedx[i]->GetYaxis()->SetTitle("dEdx (arb. Units)");
        
        fTrackCutQA[i]->Add(fITSdedx[i]);
        
        TString NSigTPCName = Form("NSigTPC_%s", sName[i].Data());
        fNSigTPC[i] = new TH2F(NSigTPCName.Data(), NSigTPCName.Data(), pBins, pmin, pmax, 3.*twoDBins, -60., 60.);
        fNSigTPC[i]->Sumw2();
        fNSigTPC[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fNSigTPC[i]->GetYaxis()->SetTitle("n#sigma_{TPC}");
        
        fTrackCutQA[i]->Add(fNSigTPC[i]);
        
        TString NSigTOFName = Form("NSigTOF_%s", sName[i].Data());
        fNSigTOF[i] = new TH2F(NSigTOFName.Data(), NSigTOFName.Data(), pBins, pmin, pmax, 3.*twoDBins, -60., 60.);
        fNSigTOF[i]->Sumw2();
        fNSigTOF[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fNSigTOF[i]->GetYaxis()->SetTitle("n#sigma_{TOF}");
        
        fTrackCutQA[i]->Add(fNSigTOF[i]);
        
        TString NSigITSName = Form("NSigITS_%s", sName[i].Data());
        fNSigITS[i] = new TH2F(NSigITSName.Data(), NSigITSName.Data(), 45, pmin, pmax, 3.*twoDBins, -60., 60.);
        fNSigITS[i]->Sumw2();
        fNSigITS[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fNSigITS[i]->GetYaxis()->SetTitle("n#sigma_{ITS}");
        
        fTrackCutQA[i]->Add(fNSigITS[i]);
        
        TString TPCstatusname = Form("TPCStatus_%s",sName[i].Data());
        fTPCStatus[i] = new TH1F(TPCstatusname.Data(),TPCstatusname.Data(),4,0,4);
        fTPCStatus[i]->GetXaxis()->SetBinLabel(1, "kDetNoSignal");
        fTPCStatus[i]->GetXaxis()->SetBinLabel(2, "kDetPidOk");
        fTPCStatus[i]->GetXaxis()->SetBinLabel(3, "kDetMismatch");
        fTPCStatus[i]->GetXaxis()->SetBinLabel(4, "kDetNoParams");
        
        fTrackCutQA[i]->Add(fTPCStatus[i]);
        
        TString TOFstatusname = Form("TOFStatus_%s",sName[i].Data());
        fTOFStatus[i] = new TH1F(TOFstatusname.Data(),TOFstatusname.Data(),4,0,4);
        fTOFStatus[i]->GetXaxis()->SetBinLabel(1, "kDetNoSignal");
        fTOFStatus[i]->GetXaxis()->SetBinLabel(2, "kDetPidOk");
        fTOFStatus[i]->GetXaxis()->SetBinLabel(3, "kDetMismatch");
        fTOFStatus[i]->GetXaxis()->SetBinLabel(4, "kDetNoParams");
        
        fTrackCutQA[i]->Add(fTOFStatus[i]);
        
        TString ITSstatusname = Form("ITSStatus_%s",sName[i].Data());
        fITSStatus[i] = new TH1F(ITSstatusname.Data(),ITSstatusname.Data(),4,0,4);
        fITSStatus[i]->GetXaxis()->SetBinLabel(1, "kDetNoSignal");
        fITSStatus[i]->GetXaxis()->SetBinLabel(2, "kDetPidOk");
        fITSStatus[i]->GetXaxis()->SetBinLabel(3, "kDetMismatch");
        fITSStatus[i]->GetXaxis()->SetBinLabel(4, "kDetNoParams");
        
        fTrackCutQA[i]->Add(fITSStatus[i]);
        
        TString TrackChi2perNDFname = Form("TrackChi2perNDF_%s",sName[i].Data());
        fTrackChi2perNDF[i] = new TH1F(TrackChi2perNDFname.Data(),TrackChi2perNDFname.Data(),200,-1,10);
        fTrackChi2perNDF[i]->Sumw2();
        fTrackChi2perNDF[i]->GetXaxis()->SetTitle("#chi^{2}/NDF");
        fTrackCutQA[i]->Add(fTrackChi2perNDF[i]);
        
        TString TPCClsCPileUpName=Form("TPCClsCPileUp_%s",sName[i].Data());
        fTPCClsCPiluUp[i]=new TH2F(TPCClsCPileUpName.Data(),
                                   TPCClsCPileUpName.Data(),
                                   15,0,15,200,0,200);
        fTPCClsCPiluUp[i]->Sumw2();
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(1,"Hit ITS Layer 1");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(2,"Hit ITS Layer 2");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(3,"Hit ITS Layer 3");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(4,"Hit ITS Layer 4");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(5,"Hit ITS Layer 5");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(6,"Hit ITS Layer 6");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(7,"TOF Bunch.Xing");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(8,"no Hit ITS Layer 1");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(9,"no Hit ITS Layer 2");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(10,"no Hit ITS Layer 3");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(11,"no Hit ITS Layer 4");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(12,"no Hit ITS Layer 5");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(13,"no Hit ITS Layer 6");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(14,"no TOF Bunch.Xing");
        fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(15,"all ITS-TOF Req false");
        fTrackCutQA[i]->Add(fTPCClsCPiluUp[i]);
        
        TString ITSShrdClsPileUpName=Form("ITSShrdClsPileUp_%s",sName[i].Data());
        fITShrdClsPileUp[i]=new TH2F(ITSShrdClsPileUpName.Data(),
                                     ITSShrdClsPileUpName.Data(),
                                     15,0,15,12,0,12);
        fITShrdClsPileUp[i]->Sumw2();
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(1,"Hit ITS Layer 1");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(2,"Hit ITS Layer 2");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(3,"Hit ITS Layer 3");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(4,"Hit ITS Layer 4");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(5,"Hit ITS Layer 5");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(6,"Hit ITS Layer 6");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(7,"TOF Bunch.Xing");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(8,"no Hit ITS Layer 1");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(9,"no Hit ITS Layer 2");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(10,"no Hit ITS Layer 3");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(11,"no Hit ITS Layer 4");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(12,"no Hit ITS Layer 5");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(13,"no Hit ITS Layer 6");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(14,"no Hit TOF Bunch.Xing");
        fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(15,"all Hit ITS-TOF Req false");
        
        fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(1,"Has shared Cls Layer 1");
        fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(2,"Has shared Cls Layer 2");
        fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(3,"Has shared Cls Layer 3");
        fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(4,"Has shared Cls Layer 4");
        fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(5,"Has shared Cls Layer 5");
        fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(6,"Has shared Cls Layer 6");
        fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(7,"Has no Shard Cls 1");
        fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(8,"Has no Shard Cls 2");
        fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(9,"Has no Shard Cls 3");
        fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(10,"Has no Shard Cls 4");
        fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(11,"Has no Shard Cls 5");
        fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(12,"Has no Shard Cls 6");
        fTrackCutQA[i]->Add(fITShrdClsPileUp[i]);
        
    }
    if (CombSig) {
        fNSigCom = new TH3F("NSigComb","NSigComb",14,0.5,4.0,300,-30,30,300,-30,30);
        fHistList->Add(fNSigCom);
        fNSigCom->Sumw2();
        
    } else {
        fNSigCom=0;
    }
    
    if (DCADist) {
        TString dcaPBinName = Form("DCAXYPBinningTot");
        fDCAXYPBins = new TH2F(dcaPBinName.Data(), dcaPBinName.Data(),45,0.0,4.5,1000,-5,5);
        fDCAXYPBins->Sumw2();
        fDCAXYPBins->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fDCAXYPBins->GetYaxis()->SetTitle("dca_{XY}");
        fHistList->Add(fDCAXYPBins);
    } else {
        fDCAXYPBins=0;
    }
    
    if(PlotPmass2dca3D){
        fP_mass2_DCAxyHist = new TH3F("P_mass2_DCAxy", "P_mass2_DCAxy", 45, pmin, pmax, 5*twoDBins, -10., 10., 2.5*twoDBins, -5., 5.);
        fHistList->Add(fP_mass2_DCAxyHist);
        fP_mass2_DCAxyHist->Sumw2();
        fP_mass2_DCAxyHist->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
        fP_mass2_DCAxyHist->GetYaxis()->SetTitle("m^{2}");
        fP_mass2_DCAxyHist->GetZaxis()->SetTitle("DCAxy (cm)");
    }else{
        fP_mass2_DCAxyHist=0;
    }
}
AliLightNTrackHist::~AliLightNTrackHist() {
    for (int i=0;i<2;++i) {
        delete fTrackCutQA[i];
        delete fpDist[i];
        delete fpTPCDist[i];
        delete fDiff_p_pTPC[i];
        delete fetaDist[i];
        delete fphiDist[i];
        delete fTPCCls[i];
        delete fTPCClsS[i];
        delete fShrdClsITS[i];
        delete fDCAxy[i];
        delete fDCAz[i];
        delete fMass2sqHist[i];
        delete fTPCCrossedRows[i];
        delete fTPCRatio[i];
        delete fTPCdedx[i];
        delete fTOFbeta[i];
        delete fNSigTPC[i];
        delete fNSigTOF[i];
        delete fTPCStatus[i];
        delete fTOFStatus[i];
        delete fTPCClsCPiluUp[i];
        delete fITShrdClsPileUp[i];
        delete fTrackChi2perNDF[i];
    }
    delete fHistList;
    delete fP_mass2_DCAxyHist;
}

void AliLightNTrackHist::FillNSigComb(double p,double nSigTPC,double nSigTOF){
    fNSigCom->Fill(p,nSigTPC,nSigTOF);
}

void AliLightNTrackHist::FillDCAXYPBins(double p,double dcaxy){
    fDCAXYPBins->Fill(p,dcaxy);
}
void AliLightNTrackHist::FillP_mass2_DCAxy(double p, double mass2,double dcaxyCut){
    fP_mass2_DCAxyHist->Fill(p,mass2,dcaxyCut);
}
