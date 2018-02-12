/*
 * AliFemtoDreamCascadeHist.cxx
 *
 *  Created on: Jan 11, 2018
 *      Author: gu74req
 */

#include "AliFemtoDreamCascadeHist.h"

ClassImp(AliFemtoDreamCascadeHist)

AliFemtoDreamCascadeHist::AliFemtoDreamCascadeHist()
:fHistList()
,fCutCounter(0)
,fConfig(0)
,fInvMassPtXi(0)
,fInvMassPtv0(0)
{
  for (int i=0;i<2;++i) {
    fCascadeQA[i]=0;
    fInvMass[i]=0;
    fInvMassv0[i]=0;
    fXiPt[i]=0;;
    fP_Y_Xi[i]=0;;
    fDCAXiDaug[i]=0;;
    fMinDistVtxBach[i]=0;;
    fCPAXi[i]=0;;
    fTransRadiusXi[i]=0;;
    fv0MaxDCADaug[i]=0;;
    fCPAv0[i]=0;;
    fTransRadiusv0[i]=0;;
    fMinDistVtxv0[i]=0;;
    fMinDistVtxv0DaugPos[i]=0;;
    fMinDistVtxv0DaugNeg[i]=0;;
  }
}
AliFemtoDreamCascadeHist::AliFemtoDreamCascadeHist(double mass) {
  fHistList=new TList();
  fHistList->SetName("Cascade");
  fHistList->SetOwner();

  fCutCounter = new TH1F("CutCounter", "Cut Counter", 40, 0, 40);
  fCutCounter->GetXaxis()->SetBinLabel(1, "Input");
  fCutCounter->GetXaxis()->SetBinLabel(2, "Charge");
  fCutCounter->GetXaxis()->SetBinLabel(3, "Id Neg Pos");
  fCutCounter->GetXaxis()->SetBinLabel(4, "Id Neg Bach");
  fCutCounter->GetXaxis()->SetBinLabel(5, "Id Pos Bach");
  fCutCounter->GetXaxis()->SetBinLabel(6, "ReqITSTOFHit Pos");
  fCutCounter->GetXaxis()->SetBinLabel(7, "ReqITSTOFHit Neg");
  fCutCounter->GetXaxis()->SetBinLabel(8, "ReqITSTOFHit Bach");
  fCutCounter->GetXaxis()->SetBinLabel(9, "ReqTPCRefit Pos");
  fCutCounter->GetXaxis()->SetBinLabel(10, "ReqTPCRefit Neg");
  fCutCounter->GetXaxis()->SetBinLabel(11, "ReqTPCRefit Bach");
  fCutCounter->GetXaxis()->SetBinLabel(12, "Crossed Rows Pos");
  fCutCounter->GetXaxis()->SetBinLabel(13, "Crossed Rows Neg");
  fCutCounter->GetXaxis()->SetBinLabel(14, "Crossed Rows Bach");
  fCutCounter->GetXaxis()->SetBinLabel(15, "NClst TPC Pos");
  fCutCounter->GetXaxis()->SetBinLabel(16, "NClst TPC Neg");
  fCutCounter->GetXaxis()->SetBinLabel(17, "NClst TPC Bach");
  fCutCounter->GetXaxis()->SetBinLabel(18, "nSigmaTPC Pos");
  fCutCounter->GetXaxis()->SetBinLabel(19, "nSigmaTPC Neg");
  fCutCounter->GetXaxis()->SetBinLabel(20, "nSigmaTPC Bach");
  fCutCounter->GetXaxis()->SetBinLabel(21, "P_{T,min} Pos");
  fCutCounter->GetXaxis()->SetBinLabel(22, "P_{T,min} Neg");
  fCutCounter->GetXaxis()->SetBinLabel(23, "P_{T,min} Bach");
  fCutCounter->GetXaxis()->SetBinLabel(24, "DCAXiDaug");
  fCutCounter->GetXaxis()->SetBinLabel(25, "Min Dist Prim Vtx Bach");
  fCutCounter->GetXaxis()->SetBinLabel(26, "CPA Xi");
  fCutCounter->GetXaxis()->SetBinLabel(27, "Xi TransRadius");
  fCutCounter->GetXaxis()->SetBinLabel(28, "DCA v0 Vtx Daug");
  fCutCounter->GetXaxis()->SetBinLabel(29, "CPA v0");
  fCutCounter->GetXaxis()->SetBinLabel(30, "v0 TransRadius");
  fCutCounter->GetXaxis()->SetBinLabel(31, "Min Dist Prim Vtx v0");
  fCutCounter->GetXaxis()->SetBinLabel(32, "Min Dist Prim Vtx v0 Daughter");
  fCutCounter->GetXaxis()->SetBinLabel(33, "Inv Mass Lambda");
  fCutCounter->GetXaxis()->SetBinLabel(34, "Inv Mass Xi");

  fHistList->Add(fCutCounter);

  fConfig=new TProfile("Config","Config",22,0,22);
  fConfig->SetStats(0);
  fConfig->GetXaxis()->SetBinLabel(1,"Require ITSTOF Hit");
  fConfig->GetXaxis()->SetBinLabel(2,"Require TOF Refit Daug");
  fConfig->GetXaxis()->SetBinLabel(3,"Ratio CR");
  fConfig->GetXaxis()->SetBinLabel(4,"NCls TPC");
  fConfig->GetXaxis()->SetBinLabel(5,"N#Sigma_{TPC}");
  fConfig->GetXaxis()->SetBinLabel(6,"P_{T,Daug-min}");
  fConfig->GetXaxis()->SetBinLabel(7,"#Xi_{Mass}");
  fConfig->GetXaxis()->SetBinLabel(8,"#Xi_{Width}");
  fConfig->GetXaxis()->SetBinLabel(9,"#Xi Charge");
  fConfig->GetXaxis()->SetBinLabel(10,"DCA #Xi Daug");
  fConfig->GetXaxis()->SetBinLabel(11,"Min Dist. Bach To Vtx");
  fConfig->GetXaxis()->SetBinLabel(12,"CPA #Xi");
  fConfig->GetXaxis()->SetBinLabel(13,"Min Trans. Radius #Xi");
  fConfig->GetXaxis()->SetBinLabel(14,"Max Trans. Radius #Xi");
  fConfig->GetXaxis()->SetBinLabel(15,"v0_{Mass}");
  fConfig->GetXaxis()->SetBinLabel(16,"v0_{Width}");
  fConfig->GetXaxis()->SetBinLabel(17,"DCA v0 Daug");
  fConfig->GetXaxis()->SetBinLabel(18,"CPAv0");
  fConfig->GetXaxis()->SetBinLabel(19,"Min Trans. Radius v0");
  fConfig->GetXaxis()->SetBinLabel(20,"Max Trans. Radius v0");
  fConfig->GetXaxis()->SetBinLabel(21,"Min Dist. v0 To Vtx");
  fConfig->GetXaxis()->SetBinLabel(22,"Min Dist. v0 Daug to Vtx");
  fHistList->Add(fConfig);

  fInvMassPtXi=new TH2F("InvMassXiPt","InvMassXiPt",19,0.6,6.5,500,mass*0.9,mass*1.3);
  fInvMassPtXi->Sumw2();
  fInvMassPtXi->GetXaxis()->SetTitle("P_{T}");
  fInvMassPtXi->GetYaxis()->SetTitle("Inv Mass");
  fHistList->Add(fInvMassPtXi);

  fInvMassPtv0=new TH2F("InvMassv0Pt","InvMassv0Pt",19,0.6,6.5,400,0.9,1.2);
  fInvMassPtv0->Sumw2();
  fInvMassPtv0->GetXaxis()->SetTitle("P_{T}");
  fInvMassPtv0->GetYaxis()->SetTitle("Inv Mass");
  fHistList->Add(fInvMassPtv0);

  TString sName[2]={"before","after"};
  for (int i=0;i<2;++i) {
    fCascadeQA[i]=new TList();
    fCascadeQA[i]->SetOwner();
    fCascadeQA[i]->SetName(sName[i].Data());
    fHistList->Add(fCascadeQA[i]);

    TString InvMassXiName=Form("InvariantMassXi_%s",sName[i].Data());
    fInvMass[i]=new TH1F(InvMassXiName.Data(),InvMassXiName.Data(),500,mass*0.9,mass*1.3);
    fInvMass[i]->Sumw2();
    fCascadeQA[i]->Add(fInvMass[i]);

    TString InvMassv0Name=Form("InvariantMassv0_%s",sName[i].Data());
    fInvMassv0[i]=new TH1F(InvMassv0Name.Data(),InvMassv0Name.Data()
                           ,500,1.116*0.9,1.116*1.3);
    fInvMassv0[i]->Sumw2();
    fCascadeQA[i]->Add(fInvMassv0[i]);

    TString XiPtName=Form("XiPt_%s",sName[i].Data());
    fXiPt[i]=new TH1F(XiPtName.Data(),XiPtName.Data(),19,0.6,6.5);
    fXiPt[i]->Sumw2();
    fCascadeQA[i]->Add(fXiPt[i]);

    TString P_Y_XiName=Form("Xi_P_Y_%s",sName[i].Data());
    fP_Y_Xi[i]=new TH2F(P_Y_XiName.Data(),P_Y_XiName.Data(),50,-3.,3.,50,0.,5.);
    fP_Y_Xi[i]->GetXaxis()->SetName("Rapidity Y");
    fP_Y_Xi[i]->GetYaxis()->SetName("P");
    fP_Y_Xi[i]->Sumw2();
    fCascadeQA[i]->Add(fP_Y_Xi[i]);

    TString DCAXiDaugName=Form("DCAXiDaug_%s",sName[i].Data());
    fDCAXiDaug[i]=new TH1F(DCAXiDaugName.Data(),DCAXiDaugName.Data(),50,0,10);
    fDCAXiDaug[i]->Sumw2();
    fCascadeQA[i]->Add(fDCAXiDaug[i]);

    TString MinDistVtxBachName=Form("MinDistVtxBach_%s",sName[i].Data());
    fMinDistVtxBach[i]=new TH1F(MinDistVtxBachName.Data(),
                                MinDistVtxBachName.Data(),50,0,10);
    fMinDistVtxBach[i]->Sumw2();
    fCascadeQA[i]->Add(fMinDistVtxBach[i]);

    TString CPAXiName=Form("CPAXi_%s",sName[i].Data());
    fCPAXi[i]=new TH1F(CPAXiName.Data(),CPAXiName.Data(),100,0.97,1);
    fCPAXi[i]->Sumw2();
    fCascadeQA[i]->Add(fCPAXi[i]);

    TString TransRadiusXiName=Form("TransRadiusXi_%s",sName[i].Data());
    fTransRadiusXi[i]=new TH1F(TransRadiusXiName.Data(),
                               TransRadiusXiName.Data(),200,0,200);
    fTransRadiusXi[i]->Sumw2();
    fCascadeQA[i]->Add(fTransRadiusXi[i]);

    TString v0MaxDCADaugName=Form("v0MaxDCADaug_%s",sName[i].Data());
    fv0MaxDCADaug[i]=new TH1F(v0MaxDCADaugName.Data(),v0MaxDCADaugName.Data(),
                              50,0,10);
    fv0MaxDCADaug[i]->Sumw2();
    fCascadeQA[i]->Add(fv0MaxDCADaug[i]);

    TString CPAv0Name=Form("CPAv0_%s",sName[i].Data());
    fCPAv0[i]=new TH1F(CPAv0Name.Data(),CPAv0Name.Data(),100,0.97,1);
    fCPAv0[i]->Sumw2();
    fCascadeQA[i]->Add(fCPAv0[i]);

    TString TransRadiusv0Name=Form("TransRadiusv0_%s",sName[i].Data());
    fTransRadiusv0[i]=new TH1F(TransRadiusv0Name.Data(),
                               TransRadiusv0Name.Data(),200,0,200);
    fTransRadiusv0[i]->Sumw2();
    fCascadeQA[i]->Add(fTransRadiusv0[i]);

    TString MinDistVtxv0Name=Form("MinDistVtxv0_%s",sName[i].Data());
    fMinDistVtxv0[i]=new TH1F(MinDistVtxv0Name.Data(),MinDistVtxv0Name.Data(),
                              50,0,10);
    fMinDistVtxv0[i]->Sumw2();
    fCascadeQA[i]->Add(fMinDistVtxv0[i]);

    TString MinDistVtxv0DaugPosName=Form("MinDistVtxv0DaugPos_%s",sName[i].Data());
    fMinDistVtxv0DaugPos[i]=new TH1F(MinDistVtxv0DaugPosName.Data(),
                                  MinDistVtxv0DaugPosName.Data(),50,0,10);
    fMinDistVtxv0DaugPos[i]->Sumw2();
    fCascadeQA[i]->Add(fMinDistVtxv0DaugPos[i]);

    TString MinDistVtxv0DaugNameNeg=Form("MinDistVtxv0DaugNeg_%s",sName[i].Data());
    fMinDistVtxv0DaugNeg[i]=new TH1F(MinDistVtxv0DaugNameNeg.Data(),
                                  MinDistVtxv0DaugNameNeg.Data(),50,0,10);
    fMinDistVtxv0DaugNeg[i]->Sumw2();
    fCascadeQA[i]->Add(fMinDistVtxv0DaugNeg[i]);
  }
}

AliFemtoDreamCascadeHist::~AliFemtoDreamCascadeHist() {

}

