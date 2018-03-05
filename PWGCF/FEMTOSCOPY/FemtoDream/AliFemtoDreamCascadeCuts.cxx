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
  if (pass) {
    if (!fBachCuts->isSelected(casc->GetBach())) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(7);
    }
  }
  if (pass&&fcutDCAXiDaug) {
    if (casc->GetXiDCADaug()>fMaxDCAXiDaug) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(8);
    }
  }
  if (pass&&fcutMinDistVtxBach) {
    if (casc->BachDCAPrimVtx()<fMinDistVtxBach) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(9);
    }
  }
  if (pass&&fcutCPAXi) {
    if (casc->GetCPA()<fCPAXi) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(10);
    }
  }
  if (pass&&fcutXiTransRadius) {
    if ((casc->GetXiTransverseRadius()<fMinXiTransRadius)||
        (casc->GetXiTransverseRadius()>fMaxXiTransRadius)) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(11);
    }
  }
  if (pass&&fcutv0MaxDCADaug) {
    if (casc->Getv0DCADaug()>fv0MaxDCADaug) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(12);
    }
  }
  if (pass&&fcutCPAv0) {
    if (casc->Getv0CPA()<fCPAv0) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(13);
    }
  }
  if (pass&&fcutv0TransRadius) {
    if ((casc->Getv0TransverseRadius()<fMinv0TransRadius)||
        (casc->Getv0TransverseRadius()>fMaxv0TransRadius)) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(14);
    }
  }
  if (pass&&fcutv0MinDistVtx) {
    if (casc->Getv0DCAPrimVtx()<fv0MinDistVtx) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(15);
    }
  }
  if (pass&&fcutv0DaugMinDistVtx) {
    if ((casc->Getv0PosToPrimVtx()<fv0DaugMinDistVtx)||
        (casc->Getv0NegToPrimVtx()<fv0DaugMinDistVtx)) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(16);
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
      if (!fMinimalBooking)fHist->FillCutCounter(17);
    }
  }
  if (pass&&fRejOmega) {
    if ((casc->GetOmegaMass()>(fRejOmegaMass-fRejOmegaWidth))&&
        (casc->GetOmegaMass()<(fRejOmegaMass+fRejOmegaWidth))) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(18);
    }
  }
  if (pass) {
    if (!fMinimalBooking)fHist->FillInvMassPtXi(casc->GetPt(),casc->GetXiMass());
  }
  if (pass&&fcutXiMass) {
    if ((casc->GetXiMass()<(fXiMass-fXiMassWidth))||
        (casc->GetXiMass()>(fXiMass+fXiMassWidth))) {
      pass=false;
    } else {
      if (!fMinimalBooking)fHist->FillCutCounter(19);
    }
  }
  casc->SetUse(pass);
  casc->GetNegDaug()->SetUse(pass);
  casc->GetPosDaug()->SetUse(pass);
  casc->GetBach()->SetUse(pass);
  if (!fMinimalBooking) {
    BookQA(casc);
    if (fMCData) {
      BookMCQA(casc);
    }
  }
  return pass;
}

void AliFemtoDreamCascadeCuts::Init(bool MinimalBooking) {
  fMinimalBooking=MinimalBooking;
  fNegCuts->Init(fMinimalBooking);
  fPosCuts->Init(fMinimalBooking);
  fBachCuts->Init(fMinimalBooking);
  if (!fMinimalBooking) {
    fHist=new AliFemtoDreamCascadeHist(fXiMass);
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
  }
}

void AliFemtoDreamCascadeCuts::BookQA(AliFemtoDreamCascade *casc) {
  if (!fMinimalBooking) {
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
        fHist->FillPodolandski(i,casc->GetXiAlpha(),casc->GetPtArmXi());
      }
    }
    fNegCuts->BookQA(casc->GetNegDaug());
    fPosCuts->BookQA(casc->GetPosDaug());
    fBachCuts->BookQA(casc->GetBach());
  }
  return;
}

void AliFemtoDreamCascadeCuts::BookMCQA(AliFemtoDreamCascade *casc) {
  if (!fMinimalBooking) {
    double pT=casc->GetPt();
    //  if (casc->GetHasDaughters()) {
    //    double etaNegDaug=casc->GetEta().at(1);
    //    double etaPosDaug=casc->GetEta().at(2);
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
      fMCHist->FillMCXiDecayLength(iFill,pT,casc->GetXiDecayLength());
      fMCHist->FillMCOmegaDecayLength(iFill,pT,casc->GetOmegaDecayLength());
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

