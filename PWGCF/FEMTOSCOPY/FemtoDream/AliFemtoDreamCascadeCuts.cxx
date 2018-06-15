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
,fMCHist(0)
,fNegCuts(0)
,fPosCuts(0)
,fBachCuts(0)
,fHistList(0)
,fMCHistList(0)
,fMinimalBooking(false)
,fMCData(false)
,fContribSplitting(false)
,fRunNumberQA(false)
,fMinRunNumber(0)
,fMaxRunNumber(0)
,fCutPt(false)
,fPtMin(0)
,fPtMax(0)
,fCutPtv0(false)
,fPtMinv0(0)
,fPtMaxv0(0)
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
,fRejOmega(false)
,fRejOmegaMass(0)
,fRejOmegaWidth(0)
,fPDGCasc(0)
,fPDGv0(0)
,fPDGPosDaug(0)
,fPDGNegDaug(0)
,fPDGBachDaug(0)
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

AliFemtoDreamCascadeCuts* AliFemtoDreamCascadeCuts::XiCuts(
    bool isMC,bool contribSplitting) {
  AliFemtoDreamCascadeCuts *XiCuts=new AliFemtoDreamCascadeCuts();
  XiCuts->SetIsMonteCarlo(isMC);
  XiCuts->SetContributionSplitting(contribSplitting);
  XiCuts->SetXiMassRange(1.322,0.005);
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
  XiCuts->SetRejectOmegas(1672,0.005);
  XiCuts->SetPtRangeXi(0.3,999.9);
  return XiCuts;
}

bool AliFemtoDreamCascadeCuts::isSelected(AliFemtoDreamCascade *casc) {
  bool pass=true;
  if (!fMinimalBooking)fHist->FillCutCounter(0);
  if (fcutXiCharge) {
    if (!(casc->GetCharge().at(0)==fXiCharge)) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(1);
    }
  }
  if (pass) {
    std::vector<int> IDTracks=casc->GetIDTracks();
    if (IDTracks[0]==IDTracks[1]) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(2);
    }
    if (IDTracks[0]==IDTracks[2]) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(3);
    }
    if (IDTracks[1]==IDTracks[2]) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(4);
    }
  }
  if (pass) {
    if (!fNegCuts->isSelected(casc->GetNegDaug())) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(5);
    }
  }
  if (pass) {
    if (!fPosCuts->isSelected(casc->GetPosDaug())) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(6);
    }
  }
  if (pass&&fcutv0MaxDCADaug) {
    if (casc->Getv0DCADaug()>fv0MaxDCADaug) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(7);
    }
  }
  if (pass&&fcutCPAv0) {
    if (casc->Getv0CPA()<fCPAv0) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(8);
    }
  }
  if (pass&&fCutPtv0) {
    if ((casc->Getv0Pt()<fPtMinv0) || (fPtMaxv0<casc->Getv0Pt())) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(9);
    }
  }
  if (pass&&fcutv0TransRadius) {
    if ((casc->Getv0TransverseRadius()<fMinv0TransRadius)||
        (casc->Getv0TransverseRadius()>fMaxv0TransRadius)) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(10);
    }
  }
  if (pass&&fcutv0MinDistVtx) {
    if (casc->Getv0DCAPrimVtx()<fv0MinDistVtx) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(11);
    }
  }
  if (pass&&fcutv0DaugMinDistVtx) {
    if ((casc->Getv0PosToPrimVtx()<fv0DaugMinDistVtx)||
        (casc->Getv0NegToPrimVtx()<fv0DaugMinDistVtx)) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(12);
    }
  }
  if (pass) {
    if (!fMinimalBooking)fHist->FillInvMassPtv0(casc->GetPt(),casc->Getv0Mass());
  }
  if (pass&&fcutv0Mass) {
    if ((casc->Getv0Mass()<(fv0Mass-fv0Width))||
        (casc->Getv0Mass()>(fv0Mass+fv0Width))) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(13);
    }
  }
  if (pass) {
    if (!fBachCuts->isSelected(casc->GetBach())) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(14);
    }
  }
  if (pass&&fcutDCAXiDaug) {
    if (casc->GetXiDCADaug()>fMaxDCAXiDaug) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(15);
    }
  }
  if (pass&&fcutMinDistVtxBach) {
    if (casc->BachDCAPrimVtx()<fMinDistVtxBach) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(16);
    }
  }
  if (pass&&fcutCPAXi) {
    if (casc->GetCPA()<fCPAXi) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(17);
    }
  }
  if (pass&&fcutXiTransRadius) {
    if ((casc->GetXiTransverseRadius()<fMinXiTransRadius)||
        (casc->GetXiTransverseRadius()>fMaxXiTransRadius)) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(18);
    }
  }
  if (pass&&fRejOmega) {
    if ((casc->GetOmegaMass()>(fRejOmegaMass-fRejOmegaWidth))&&
        (casc->GetOmegaMass()<(fRejOmegaMass+fRejOmegaWidth))) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(19);
    }
  }
  if (pass&&fCutPt) {
    if ((casc->GetPt()<fPtMin)||(fPtMax<casc->GetPt())) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(20);
    }
  }
  if (pass) {
    fHist->FillInvMassPt(casc->GetPt(),casc->GetMass());
    if (fRunNumberQA) {
      fHist->FillInvMassPerRunNumber(casc->GetEvtNumber(),casc->GetXiMass());
    }
  }
  if (pass&&fcutXiMass) {
    if ((casc->GetXiMass()<(fXiMass-fXiMassWidth))||
        (casc->GetXiMass()>(fXiMass+fXiMassWidth))) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(21);
    }
  }
  casc->SetUse(pass);
  casc->GetNegDaug()->SetUse(pass);
  casc->GetPosDaug()->SetUse(pass);
  casc->GetBach()->SetUse(pass);
  BookQA(casc);
  if (!fMinimalBooking) {
    if (fMCData) {
      BookMCQA(casc);
    }
  }
  return pass;
}

void AliFemtoDreamCascadeCuts::Init() {
  fNegCuts->SetMinimalBooking(fMinimalBooking);
  fNegCuts->Init();
  fPosCuts->SetMinimalBooking(fMinimalBooking);
  fPosCuts->Init();
  fBachCuts->SetMinimalBooking(fMinimalBooking);
  fBachCuts->Init();
  if (!fMinimalBooking) {
    fHist=new AliFemtoDreamCascadeHist(fXiMass,fRunNumberQA,fMinRunNumber,fMaxRunNumber);
    BookCuts();
    if (!(fNegCuts||fPosCuts||fBachCuts)) {
      AliFatal("Track Cuts Object Missing");
    }
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

    if (fMCData) {
      fMCHist=new AliFemtoDreamv0MCHist(
          400,1.0, 1.2,fContribSplitting,false);
      fMCHistList=new TList();
      fMCHistList->SetOwner();
      fMCHistList->SetName("CascadeMC");
      fMCHistList->Add(fMCHist->GetHistList());

      fNegCuts->SetMCName("NegCuts");
      fMCHistList->Add(fNegCuts->GetMCQAHists());
      fPosCuts->SetMCName("PosCuts");
      fMCHistList->Add(fPosCuts->GetMCQAHists());
      fBachCuts->SetMCName("BachCuts");
      fMCHistList->Add(fBachCuts->GetMCQAHists());
    }
  } else {
    fHist=new AliFemtoDreamCascadeHist("MinimalBooking",fXiMass);
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
}

void AliFemtoDreamCascadeCuts::BookQA(AliFemtoDreamCascade *casc) {
  if (!fMinimalBooking) {
    for (int i=0;i<2;++i) {
      if (i==0||(i==1&&casc->UseParticle())) {
        fHist->FillInvMass(i,casc->GetMass());
        fHist->FillInvMassLambda(i,casc->Getv0Mass());
        fHist->FillXiPt(i,casc->GetMomentum().Pt());
        fHist->FillMomRapXi(i,casc->GetXiRapidity(),casc->GetMomentum().Mag());
        fHist->FillMomRapOmega(i,casc->GetOmegaRapidity(),casc->GetMomentum().Mag());
        fHist->FillDCAXiPrimVtx(i,casc->GetDCAXiPrimVtx());
        fHist->FillDCAXiDaug(i,casc->GetXiDCADaug());
        fHist->FillMinDistPrimVtxBach(i,casc->BachDCAPrimVtx());
        fHist->FillCPAXi(i,casc->GetCPA());
        fHist->FillDecayLength(i,casc->GetDecayLength());
        fHist->Fillv0DecayLength(i,casc->Getv0DecayLength());
        fHist->FillTransverseRadiusXi(i,casc->GetXiTransverseRadius());
        fHist->FillMaxDCAv0Daug(i,casc->Getv0DCADaug());
        fHist->FillCPAv0(i,casc->Getv0CPA());
        fHist->FillCPAv0Xi(i,casc->Getv0XiPointingAngle());
        fHist->Fillv0Pt(i,casc->Getv0Pt());
        fHist->FillTransverseRadiusv0(i,casc->Getv0TransverseRadius());
        fHist->FillMinDistPrimVtxv0(i,casc->Getv0DCAPrimVtx());
        fHist->FillMinDistPrimVtxv0DaugPos(i,casc->Getv0PosToPrimVtx());
        fHist->FillMinDistPrimVtxv0DaugNeg(i,casc->Getv0NegToPrimVtx());
        fHist->FillPodolandski(i,casc->GetXiAlpha(),casc->GetPtArmXi());
      }
    }
  }
  if (casc->GetNegDaug()->IsSet()) fNegCuts->BookQA(casc->GetNegDaug());
  if (casc->GetPosDaug()->IsSet()) fPosCuts->BookQA(casc->GetPosDaug());
  if (casc->GetBach()->IsSet()) fBachCuts->BookQA(casc->GetBach());
  return;
}

void AliFemtoDreamCascadeCuts::BookMCQA(AliFemtoDreamCascade *casc) {
  if (!fMinimalBooking) {
    float pT=casc->GetPt();
    //  if (casc->GetHasDaughters()) {
    //    float etaNegDaug=casc->GetEta().at(1);
    //    float etaPosDaug=casc->GetEta().at(2);
    //    if (casc->GetMCPDGCode()==fPDGv0) {
    //      if (fpTmin<pT&&pT<fpTmax) {
    //        if (fPosCuts->GetEtaMin()<etaPosDaug&&etaPosDaug<fPosCuts->GetEtaMax()) {
    //          if (fNegCuts->GetEtaMin()<etaNegDaug&&etaNegDaug<fNegCuts->GetEtaMax()) {
    //            fMCHist->FillMCGen(pT);
    //          }
    //        }
    //      }
    //    }
    //  }
    if (casc->UseParticle()) {
      fMCHist->FillMCIdent(pT);
      AliFemtoDreamBasePart::PartOrigin tmpOrg=casc->GetParticleOrigin();
      if (casc->GetParticleOrigin()!=AliFemtoDreamBasePart::kFake) {
        if (casc->GetMCPDGCode()==fPDGCasc) {
          fMCHist->FillMCCorr(pT);
        } else {
          casc->SetParticleOrigin(AliFemtoDreamBasePart::kContamination);
        }
      }
      if (fContribSplitting) {
        FillMCContributions(casc);
      }
      casc->GetPosDaug()->SetParticleOrigin(casc->GetParticleOrigin());
      casc->GetNegDaug()->SetParticleOrigin(casc->GetParticleOrigin());
      casc->GetBach()->SetParticleOrigin(casc->GetParticleOrigin());
      fPosCuts->BookMC(casc->GetPosDaug());
      fNegCuts->BookMC(casc->GetNegDaug());
      fBachCuts->BookMC(casc->GetBach());
      casc->SetParticleOrigin(tmpOrg);
    }
  }
  return;
}

void AliFemtoDreamCascadeCuts::FillMCContributions(AliFemtoDreamCascade *casc) {
  if (!fMinimalBooking) {
    Double_t pT = casc->GetPt();
    Int_t iFill = -1;
    switch(casc->GetParticleOrigin()){
      case AliFemtoDreamBasePart::kPhysPrimary:
        fMCHist->FillMCPrimary(pT);
        iFill = 0;
        break;
      case AliFemtoDreamBasePart::kWeak:
        fMCHist->FillMCFeeddown(pT, TMath::Abs(casc->GetMotherWeak()));
        iFill = 1;
        break;
      case AliFemtoDreamBasePart::kMaterial:
        fMCHist->FillMCMaterial(pT);
        iFill = 2;
        break;
      case AliFemtoDreamBasePart::kContamination:
        fMCHist->FillMCCont(pT);
        iFill = 3;
        break;
      case AliFemtoDreamBasePart::kFake:
        fMCHist->FillMCCont(pT);
        iFill = 4;
        break;
      default:
        AliFatal("Type Not implemented");
        break;
    }
    if (iFill!=-1) {
      fMCHist->FillMCpT(iFill,pT);
      fMCHist->FillMCEta(iFill,casc->GetEta().at(0));
      fMCHist->FillMCPhi(iFill,casc->GetPhi().at(0));
      fMCHist->FillMCPodolanski(iFill,casc->GetPtArmXi(),casc->GetXiAlpha());
      fMCHist->FillMCCosPoint(iFill,pT,casc->GetCPA());
      fMCHist->FillMCBachDCAToPV(iFill,pT,casc->BachDCAPrimVtx());
      fMCHist->FillMCv0DecayLength(iFill,pT,casc->Getv0DecayLength());
      fMCHist->FillMCv0CPA(iFill,pT,casc->Getv0CPA());
      fMCHist->FillMCDecayLength(iFill,pT,casc->GetDecayLength());
      fMCHist->FillMCXiRapidity(iFill,casc->GetMomentum().Mag(),casc->GetXiRapidity());
      fMCHist->FillMCOmegaRapidity(iFill,casc->GetMomentum().Mag(),casc->GetOmegaRapidity());
      fMCHist->FillMCTransverseRadius(iFill,pT,casc->GetXiTransverseRadius());
      fMCHist->FillMCDCAPosDaugPrimVtx(iFill,pT,casc->Getv0PosToPrimVtx());
      fMCHist->FillMCDCANegDaugPrimVtx(iFill,pT,casc->Getv0NegToPrimVtx());
      fMCHist->FillMCDCADaugVtx(iFill,pT,casc->GetXiDCADaug());
      fMCHist->FillMCInvMass(iFill,casc->Getv0Mass());
      fMCHist->FillMCXiInvMass(iFill,pT,casc->GetXiMass());
      fMCHist->FillMCOmegaInvMass(iFill,pT,casc->GetOmegaMass());

    } else {
      std::cout << "this should not happen \n";
    }
  }
}

void AliFemtoDreamCascadeCuts::BookCuts() {
  if (fCutPt) {
    fHist->FillConfig(0,fPtMin);
    fHist->FillConfig(1,fPtMax);
  }
  if (fCutPtv0) {
    fHist->FillConfig(2,fPtMinv0);
    fHist->FillConfig(3,fPtMaxv0);
  }
  if (fcutXiMass) {
    fHist->FillConfig(4,fXiMass);
    fHist->FillConfig(5,fXiMassWidth);
  }
  if (fcutXiCharge) {
    fHist->FillConfig(6,fXiCharge);
  }
  if (fcutDCAXiDaug) {
    fHist->FillConfig(7,fMaxDCAXiDaug);
  }
  if (fcutMinDistVtxBach) {
    fHist->FillConfig(8,fMinDistVtxBach);
  }
  if (fcutCPAXi) {
    fHist->FillConfig(9,fCPAXi);
  }
  if (fcutXiTransRadius) {
    fHist->FillConfig(10,fMinXiTransRadius);
    fHist->FillConfig(11,fMaxXiTransRadius);
  }
  if (fcutv0Mass) {
    fHist->FillConfig(12,fv0Mass);
    fHist->FillConfig(13,fv0Width);
  }
  if (fcutv0MaxDCADaug) {
    fHist->FillConfig(14,fv0MaxDCADaug);
  }
  if (fcutCPAv0) {
    fHist->FillConfig(15,fCPAv0);
  }
  if (fcutv0TransRadius) {
    fHist->FillConfig(16,fMinv0TransRadius);
    fHist->FillConfig(17,fMaxv0TransRadius);
  }
  if (fcutv0MinDistVtx) {
    fHist->FillConfig(18,fv0MinDistVtx);
  }
  if (fcutv0DaugMinDistVtx) {
    fHist->FillConfig(19,fv0DaugMinDistVtx);
  }
  if (fRejOmega) {
    fHist->FillConfig(20,fRejOmegaMass);
    fHist->FillConfig(21,fRejOmegaWidth);
  }
}

