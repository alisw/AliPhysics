/*
 * AliFemtoDreamCascadeCuts.cxx
 *
 *  Created on: Jan 11, 2018
 *      Author: gu74req
 */

#include "vector"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliLog.h"
ClassImp(AliFemtoDreamCascadeCuts)
AliFemtoDreamCascadeCuts::AliFemtoDreamCascadeCuts()
:fHist(0)
,fNegCuts(0)
,fPosCuts(0)
,fBachCuts(0)
,fHistList(0)
,fcutXiMass(false)
,fXiMass(0)
,fXiMassWidth(0)
,fcutXiCharge(false)
,fXiCharge(0)
,fcutDCAXiDaug(false)
,fMaxDCAXiDaug(0)
,fcutMinDistVtxBach(false)
,fMinDistVtxBach(0)
,fcutCPAXi(false)
,fCPAXi(0)
,fcutXiTransRadius(false)
,fMinXiTransRadius(0)
,fMaxXiTransRadius(0)
,fcutv0Mass(false)
,fv0Mass(0)
,fv0Width(0)
,fcutv0MaxDCADaug(false)
,fv0MaxDCADaug(0)
,fcutCPAv0(false)
,fCPAv0(0)
,fcutv0TransRadius(false)
,fMinv0TransRadius(0)
,fMaxv0TransRadius(0)
,fcutv0MinDistVtx(false)
,fv0MinDistVtx(0)
,fcutv0DaugMinDistVtx(false)
,fv0DaugMinDistVtx(0)
{
}

AliFemtoDreamCascadeCuts::~AliFemtoDreamCascadeCuts() {
  if (fNegCuts) {
    delete fNegCuts;
  }
  if (fPosCuts) {
    delete fPosCuts;
  }
  if (fBachCuts) {
    delete fBachCuts;
  }
}

AliFemtoDreamCascadeCuts* AliFemtoDreamCascadeCuts::XiCuts(bool isMC) {
  AliFemtoDreamCascadeCuts *XiCuts=new AliFemtoDreamCascadeCuts();
  XiCuts->SetXiMassRange(1.31486,0.05);
  XiCuts->SetCutXiDaughterDCA(1.6);
  XiCuts->SetCutXiMinDistBachToPrimVtx(0.05);
  XiCuts->SetCutXiCPA(0.97);
  XiCuts->SetCutXiTransverseRadius(0.8,200);
  XiCuts->Setv0MassRange(1.116,0.006);
  XiCuts->SetCutv0MaxDaughterDCA(1.6);
  XiCuts->SetCutv0CPA(0.97);
  XiCuts->SetCutv0TransverseRadius(1.4,200);
  XiCuts->SetCutv0MinDistToPrimVtx(0.07);
  XiCuts->SetCutv0MinDaugDistToPrimVtx(0.04);
  return XiCuts;
}

bool AliFemtoDreamCascadeCuts::isSelected(AliFemtoDreamCascade *casc) {
  bool pass=true;
  fHist->FillCutCounter(0);
  if (fcutXiCharge) {
    if (!(casc->GetCharge().at(0)==fXiCharge)) {
      pass=false;
    } else {
      fHist->FillCutCounter(1);
    }
  }
  if (pass) {
    std::vector<int> IDTracks=casc->GetIDTracks();
    if (IDTracks[0]==IDTracks[1]) {
      pass=false;
    } else {
      fHist->FillCutCounter(2);
    }
    if (IDTracks[0]==IDTracks[2]) {
      pass=false;
    } else {
      fHist->FillCutCounter(3);
    }
    if (IDTracks[1]==IDTracks[2]) {
      pass=false;
    } else {
      fHist->FillCutCounter(4);
    }
  }
  if (pass) {
    if (!fNegCuts->isSelected(casc->GetNegDaug())) {
      pass=false;
    } else {
      fHist->FillCutCounter(5);
    }
  }
  if (pass) {
    if (!fPosCuts->isSelected(casc->GetPosDaug())) {
      pass=false;
    } else {
      fHist->FillCutCounter(6);
    }
  }
  if (pass) {
    if (!fBachCuts->isSelected(casc->GetBach())) {
      pass=false;
    } else {
      fHist->FillCutCounter(7);
    }
  }
  if (pass&&fcutDCAXiDaug) {
    if (casc->GetXiDCADaug()>fMaxDCAXiDaug) {
      pass=false;
    } else {
      fHist->FillCutCounter(23);
    }
  }
  if (pass&&fcutMinDistVtxBach) {
    if (casc->BachDCAPrimVtx()<fMinDistVtxBach) {
      pass=false;
    } else {
      fHist->FillCutCounter(24);
    }
  }
  if (pass&&fcutCPAXi) {
    if (casc->GetCPA()<fCPAXi) {
      pass=false;
    } else {
      fHist->FillCutCounter(25);
    }
  }
  if (pass&&fcutXiTransRadius) {
    if ((casc->GetXiTransverseRadius()<fMinXiTransRadius)||
        (casc->GetXiTransverseRadius()>fMaxXiTransRadius)) {
      pass=false;
    } else {
      fHist->FillCutCounter(26);
    }
  }
  if (pass&&fcutv0MaxDCADaug) {
    if (casc->Getv0DCADaug()>fv0MaxDCADaug) {
      pass=false;
    } else {
      fHist->FillCutCounter(27);
    }
  }
  if (pass&&fcutCPAv0) {
    if (casc->Getv0CPA()<fCPAv0) {
      pass=false;
    } else {
      fHist->FillCutCounter(28);
    }
  }
  if (pass&&fcutv0TransRadius) {
    if ((casc->Getv0TransverseRadius()<fMinv0TransRadius)||
        (casc->Getv0TransverseRadius()>fMaxv0TransRadius)) {
      pass=false;
    } else {
      fHist->FillCutCounter(29);
    }
  }
  if (pass&&fcutv0MinDistVtx) {
    if (casc->Getv0DCAPrimVtx()<fv0MinDistVtx) {
      pass=false;
    } else {
      fHist->FillCutCounter(30);
    }
  }
  if (pass&&fcutv0DaugMinDistVtx) {
    if ((casc->Getv0PosToPrimVtx()<fv0DaugMinDistVtx)||
        (casc->Getv0NegToPrimVtx()<fv0DaugMinDistVtx)) {
      pass=false;
    } else {
      fHist->FillCutCounter(31);
    }
  }
  if (pass) {
    fHist->FillInvMassPtv0(casc->GetPt(),casc->Getv0Mass());
  }
  if (pass&&fcutv0Mass) {
    if ((casc->Getv0Mass()<(fv0Mass-fv0Width))||
        (casc->Getv0Mass()>(fv0Mass+fv0Width))) {
      pass=false;
    } else {
      fHist->FillCutCounter(32);
    }
  }
  if (pass) {
    fHist->FillInvMassPtXi(casc->GetPt(),casc->GetXiMass());
  }
  if (pass&&fcutXiMass) {
    if ((casc->GetXiMass()<(fXiMass-fXiMassWidth))||
        (casc->GetXiMass()>(fXiMass+fXiMassWidth))) {
      pass=false;
    } else {
      fHist->FillCutCounter(33);
    }
  }
  casc->SetUse(pass);
  casc->GetNegDaug()->SetUse(pass);
  casc->GetPosDaug()->SetUse(pass);
  casc->GetBach()->SetUse(pass);
  BookQA(casc);
  return pass;
}

void AliFemtoDreamCascadeCuts::Init() {
  fHist=new AliFemtoDreamCascadeHist(fXiMass);
  if (!(fNegCuts||fPosCuts||fBachCuts)) {
      AliFatal("Track Cuts Object Missing");
  }
  fNegCuts->Init();
  fPosCuts->Init();
  fBachCuts->Init();
  fHistList=new TList();
  fHistList->SetOwner();
  fHistList->SetName("CascadeCuts");
  fHistList->Add(fHist->GetHistList());
  fNegCuts->SetName("NegCuts");
  fHistList->Add(fNegCuts->GetQAHists());
  fPosCuts->SetName("PosCuts");
  fHistList->Add(fPosCuts->GetQAHists());
  fBachCuts->SetName("BachelorCuts");
  fHistList->Add(fBachCuts->GetQAHists());
}

void AliFemtoDreamCascadeCuts::BookQA(AliFemtoDreamCascade *casc) {
  for (int i=0;i<2;++i) {
    if (i==0||(i==1&&casc->UseParticle())) {
      fHist->FillInvMassXi(i,casc->GetXiMass());
      fHist->FillInvMassLambda(i,casc->Getv0Mass());
      fHist->FillXiPt(i,casc->GetMomentum().Pt());
      fHist->FillMomRapXi(i,casc->GetXiRapidity(),casc->GetMomentum().Mag());
      fHist->FillDCAXiDaug(i,casc->GetXiDCADaug());
      fHist->FillMinDistPrimVtxBach(i,casc->BachDCAPrimVtx());
      fHist->FillCPAXi(i,casc->GetCPA());
      fHist->FillTransverseRadiusXi(i,casc->GetXiTransverseRadius());
      fHist->FillMaxDCAv0Daug(i,casc->Getv0DCADaug());
      fHist->FillCPAv0(i,casc->Getv0CPA());
      fHist->FillTransverseRadiusv0(i,casc->Getv0TransverseRadius());
      fHist->FillMinDistPrimVtxv0(i,casc->Getv0DCAPrimVtx());
      fHist->FillMinDistPrimVtxv0DaugPos(i,casc->Getv0PosToPrimVtx());
      fHist->FillMinDistPrimVtxv0DaugNeg(i,casc->Getv0NegToPrimVtx());
    }
  }
  fNegCuts->BookQA(casc->GetNegDaug());
  fPosCuts->BookQA(casc->GetPosDaug());
  fBachCuts->BookQA(casc->GetBach());
  return;
}
