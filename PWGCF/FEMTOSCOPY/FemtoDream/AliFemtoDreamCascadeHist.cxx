/*
 * AliFemtoDreamCascadeHist.cxx
 *
 *  Created on: Jan 11, 2018
 *      Author: gu74req
 */

#include "AliFemtoDreamCascadeHist.h"

ClassImp(AliFemtoDreamCascadeHist)

AliFemtoDreamCascadeHist::AliFemtoDreamCascadeHist()
:fMinimalBooking(false)
,fHistList()
,fCutCounter(0)
,fConfig(0)
,fInvMassPtXi(0)
,fInvMassPtv0(0)
,fInvMassPerRunNumber()
{
  for (int i=0;i<2;++i) {
    fCascadeQA[i]=nullptr;
    fInvMass[i]=nullptr;
    fInvMassv0[i]=nullptr;
    fXiPt[i]=nullptr;
    fP_Y_Xi[i]=nullptr;
    fDCAXiDaug[i]=nullptr;
    fMinDistVtxBach[i]=nullptr;
    fCPAXi[i]=nullptr;
    fTransRadiusXi[i]=nullptr;
    fv0MaxDCADaug[i]=nullptr;
    fCPAv0[i]=nullptr;
    fv0Pt[i]=nullptr;
    fTransRadiusv0[i]=nullptr;
    fMinDistVtxv0[i]=nullptr;
    fMinDistVtxv0DaugPos[i]=nullptr;
    fMinDistVtxv0DaugNeg[i]=nullptr;
    fPodolandski[i]=nullptr;
  }
}

AliFemtoDreamCascadeHist::AliFemtoDreamCascadeHist(float mass,bool perRunnumber, int iRunMin=0, int iRunMax=0)
:fMinimalBooking(false)
{
  fHistList=new TList();
  fHistList->SetName("Cascade");
  fHistList->SetOwner();

  fCutCounter = new TH1F("CutCounter", "Cut Counter", 25, 0, 25);
  fCutCounter->GetXaxis()->SetBinLabel(1, "Input");
  fCutCounter->GetXaxis()->SetBinLabel(2, "Charge");
  fCutCounter->GetXaxis()->SetBinLabel(3, "Id Neg Pos");
  fCutCounter->GetXaxis()->SetBinLabel(4, "Id Neg Bach");
  fCutCounter->GetXaxis()->SetBinLabel(5, "Id Pos Bach");
  fCutCounter->GetXaxis()->SetBinLabel(6, "Neg Cuts Daug");
  fCutCounter->GetXaxis()->SetBinLabel(7, "Pos Cuts Daug");
  fCutCounter->GetXaxis()->SetBinLabel(8, "DCA v0 Vtx Daug");
  fCutCounter->GetXaxis()->SetBinLabel(9, "CPA v0");
  fCutCounter->GetXaxis()->SetBinLabel(10, "Pt v0");
  fCutCounter->GetXaxis()->SetBinLabel(11, "v0 TransRadius");
  fCutCounter->GetXaxis()->SetBinLabel(12, "Min Dist Prim Vtx v0");
  fCutCounter->GetXaxis()->SetBinLabel(13, "Min Dist Prim Vtx v0 Daughter");
  fCutCounter->GetXaxis()->SetBinLabel(14, "Inv Mass Lambda");
  fCutCounter->GetXaxis()->SetBinLabel(15, "Bach Cuts");
  fCutCounter->GetXaxis()->SetBinLabel(16, "DCAXiDaug");
  fCutCounter->GetXaxis()->SetBinLabel(17, "Min Dist Prim Vtx Bach");
  fCutCounter->GetXaxis()->SetBinLabel(18, "CPA Xi");
  fCutCounter->GetXaxis()->SetBinLabel(19, "Xi TransRadius");
  fCutCounter->GetXaxis()->SetBinLabel(20, "Rej Omega");
  fCutCounter->GetXaxis()->SetBinLabel(21, "PtCut");
  fCutCounter->GetXaxis()->SetBinLabel(22, "Inv Mass Xi");

  fHistList->Add(fCutCounter);

  fConfig=new TProfile("Config","Config",25,0,25);
  fConfig->SetStats(0);
  fConfig->GetXaxis()->SetBinLabel(1,"Xi Pt Min");
  fConfig->GetXaxis()->SetBinLabel(2,"Xi Pt Max");
  fConfig->GetXaxis()->SetBinLabel(3,"v0 Pt Min");
  fConfig->GetXaxis()->SetBinLabel(4,"v0 Pt Max");
  fConfig->GetXaxis()->SetBinLabel(5,"Xi Mass");
  fConfig->GetXaxis()->SetBinLabel(6,"Xi Width");
  fConfig->GetXaxis()->SetBinLabel(7,"Xi Charge");
  fConfig->GetXaxis()->SetBinLabel(8,"DCA Xi Daug Decay Vtx");
  fConfig->GetXaxis()->SetBinLabel(9,"Min DCA Bach");
  fConfig->GetXaxis()->SetBinLabel(10,"CPA Xi");
  fConfig->GetXaxis()->SetBinLabel(11,"Min Trans. Radius Xi");
  fConfig->GetXaxis()->SetBinLabel(12,"Max Trans. Radius Xi");
  fConfig->GetXaxis()->SetBinLabel(13,"v0 Mass");
  fConfig->GetXaxis()->SetBinLabel(14,"v0 Width");
  fConfig->GetXaxis()->SetBinLabel(15,"DCA V0 Daug Vtx");
  fConfig->GetXaxis()->SetBinLabel(16,"CPA v0");
  fConfig->GetXaxis()->SetBinLabel(17,"Min Trans. Radius v0");
  fConfig->GetXaxis()->SetBinLabel(18,"Max Trans. Radius v0");
  fConfig->GetXaxis()->SetBinLabel(19,"Min DCA v0");
  fConfig->GetXaxis()->SetBinLabel(20,"Min DCA v0 Daug");
  fConfig->GetXaxis()->SetBinLabel(21,"Rej. Omega Mass");
  fConfig->GetXaxis()->SetBinLabel(22,"Rej. Omega Width");
  fHistList->Add(fConfig);

  fInvMassPtXi=new TH2F("InvMassXiPt","InvMassXiPt",13,-0.2,6.3,500,mass*0.9,mass*1.3);
  fInvMassPtXi->Sumw2();
  fInvMassPtXi->GetXaxis()->SetTitle("P_{T}");
  fInvMassPtXi->GetYaxis()->SetTitle("Inv Mass");
  fHistList->Add(fInvMassPtXi);

  fInvMassPtv0=new TH2F("InvMassv0Pt","InvMassv0Pt",13,-0.2,6.3,400,0.9,1.2);
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
    fXiPt[i]=new TH1F(XiPtName.Data(),XiPtName.Data(),13,-0.2,6.3);
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

    TString v0PtName=Form("v0Pt_%s",sName[i].Data());
    fv0Pt[i]=new TH1F(v0PtName.Data(),v0PtName.Data(),100,0,10);
    fv0Pt[i]->Sumw2();
    fCascadeQA[i]->Add(fv0Pt[i]);

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

    TString PodoName=Form("Hodorlanski_%s",sName[i].Data());
    fPodolandski[i]=new TH2F(PodoName.Data(),PodoName.Data(),50,-1,1,50,0,1);
    fPodolandski[i]->Sumw2();
    fCascadeQA[i]->Add(fPodolandski[i]);
  }
  if (perRunnumber) {
    int nBins=iRunMax-iRunMin;
    TString InvMassRunNumbName="InvMassPerRunnumber";
    fInvMassPerRunNumber= new TH2F(InvMassRunNumbName.Data(),InvMassRunNumbName.Data(),
                                   nBins,iRunMin,iRunMax,500,mass*0.9,mass*1.3);
    fHistList->Add(fInvMassPerRunNumber);
  } else {
    fInvMassPerRunNumber=nullptr;
  }
}

AliFemtoDreamCascadeHist::AliFemtoDreamCascadeHist(TString minimalBooking,float mass)
:fMinimalBooking(true)
,fCutCounter(0)
,fConfig(0)
,fInvMassPtv0(0)
,fInvMassPerRunNumber()
{
  for (int i=0;i<2;++i) {
    fCascadeQA[i]=nullptr;
    fInvMass[i]=nullptr;
    fInvMassv0[i]=nullptr;
    fXiPt[i]=nullptr;
    fP_Y_Xi[i]=nullptr;
    fDCAXiDaug[i]=nullptr;
    fMinDistVtxBach[i]=nullptr;
    fCPAXi[i]=nullptr;
    fTransRadiusXi[i]=nullptr;
    fv0MaxDCADaug[i]=nullptr;
    fCPAv0[i]=nullptr;
    fv0Pt[i]=nullptr;
    fTransRadiusv0[i]=nullptr;
    fMinDistVtxv0[i]=nullptr;
    fMinDistVtxv0DaugPos[i]=nullptr;
    fMinDistVtxv0DaugNeg[i]=nullptr;
    fPodolandski[i]=nullptr;
  }
  fHistList=new TList();
  fHistList->SetOwner();
  fHistList->SetName(minimalBooking.Data());

  fInvMassPtXi=new TH2F("InvMassXiPt","InvMassXiPt",13,-0.2,6.3,500,mass*0.9,mass*1.3);
  fInvMassPtXi->Sumw2();
  fInvMassPtXi->GetXaxis()->SetTitle("P_{T}");
  fInvMassPtXi->GetYaxis()->SetTitle("Inv Mass");
  fHistList->Add(fInvMassPtXi);
}

AliFemtoDreamCascadeHist::~AliFemtoDreamCascadeHist() {

}

