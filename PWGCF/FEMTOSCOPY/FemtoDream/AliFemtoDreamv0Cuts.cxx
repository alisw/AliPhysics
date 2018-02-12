/*
 * AliFemtoDreamv0Cuts.cxx
 *
 *  Created on: Dec 12, 2017
 *      Author: gu74req
 */

#include "AliFemtoDreamv0Cuts.h"
#include "TDatabasePDG.h"
ClassImp(AliFemtoDreamv0Cuts)
AliFemtoDreamv0Cuts::AliFemtoDreamv0Cuts()
:fHistList()
,fMCHist(0)
,fHist(0)
,fPosCuts(0)
,fNegCuts(0)
,fMCData(false)
,fCPAPlots(false)
,fContribSplitting(false)
,fCutOnFlyStatus(false)
,fOnFlyStatus(false)
,fCutCharge(false)
,fCharge(0)
,fCutPt(false)
,fpTmin(0)
,fpTmax(0)
,fKaonRejection(false)
,fKaonRejLow(0)
,fKaonRejUp(0)
,fCutDecayVtxXYZ(false)
,fMaxDecayVtxXYZ(0)
,fCutTransRadius(false)
,fMinTransRadius(0)
,fMaxTransRadius(0)
,fCutMinDCADaugPrimVtx(false)
,fMinDCADaugToPrimVtx(0)
,fCutMaxDCADaugToDecayVtx(false)
,fMaxDCADaugToDecayVtx(0)
,fCutCPA(false)
,fMinCPA(0)
,fCutInvMass(false)
,fInvMassCutWidth(0)
,fAxisMinMass(0)
,fAxisMaxMass(1)
,fNumberXBins(1)
,fPDGv0(0)
,fPDGDaugP(0)
,fPDGDaugN(0)
{
}

AliFemtoDreamv0Cuts::~AliFemtoDreamv0Cuts() {
  if (fHist) {
    delete fHist;
  }
  if (fMCHist) {
    delete fMCHist;
  }
}
AliFemtoDreamv0Cuts* AliFemtoDreamv0Cuts::LambdaCuts(
    bool isMC,bool CPAPlots,bool SplitContrib) {
  AliFemtoDreamv0Cuts *LambdaCuts = new AliFemtoDreamv0Cuts();
  LambdaCuts->SetIsMonteCarlo(isMC);
  LambdaCuts->SetPlotCPADist(CPAPlots);
  LambdaCuts->SetPlotContrib(SplitContrib);

  LambdaCuts->SetCheckOnFlyStatus(false); //online = kTRUE, offline = kFALSE
  LambdaCuts->SetCutCharge(0);
  LambdaCuts->SetPtRange(0.3,999.);
  LambdaCuts->SetKaonRejection(0.48,0.515);
  LambdaCuts->SetCutMaxDecayVtx(100);
  LambdaCuts->SetCutTransverseRadius(0.2, 100);
  LambdaCuts->SetCutDCADaugToPrimVtx(0.05);
  LambdaCuts->SetCutDCADaugTov0Vtx(1.5);
  LambdaCuts->SetCutCPA(0.99);
  LambdaCuts->SetCutInvMass(0.004);
  LambdaCuts->SetAxisInvMassPlots(400,1.0, 1.2);

  return LambdaCuts;
}

bool AliFemtoDreamv0Cuts::isSelected(AliFemtoDreamv0 *v0) {
  bool pass=true;

  if (!v0->IsSet()) {
    pass=false;
  } else {
    fHist->FillTrackCounter(0);
  }
  if (pass) {
    if (!DaughtersPassCuts(v0)) {
      pass=false;
    }
  }
  if (pass) {
    if (!MotherPassCuts(v0)) {
      pass=false;
    }
  }
  if (pass) {
    if (fKaonRejection&&!RejectAsKaon(v0)) {
      pass=false;
    }
  }
  if (pass) {
    if (!CPAandMassCuts(v0)) {
      pass=false;
    }
  }
  v0->SetUse(pass);
  BookQA(v0);
  if (fMCData) {
    BookMC(v0);
  }
  return pass;
}

bool AliFemtoDreamv0Cuts::DaughtersPassCuts(AliFemtoDreamv0 *v0) {
  bool pass=true;
  bool passD1=true;
  bool passD2=true;
  if (!v0->GetHasDaughters()) {
    pass=false;
  } else {
    fHist->FillTrackCounter(1);
  }
  if (pass) {
    if (v0->GetCharge().at(1)<0&&v0->GetCharge().at(2)>0) {
      //at 1: Negative daughter, at 2: Positive Daughter should be the way it
      //was set, but sometimes it it stored in the wrong way
      fHist->FillTrackCounter(13);
      v0->Setv0Mass(CalculateInvMass(v0,fPDGDaugP,fPDGDaugN));
      passD1=fNegCuts->isSelected(v0->GetNegDaughter());
      passD2=fPosCuts->isSelected(v0->GetPosDaughter());
      if (passD1&&passD2) {
        fHist->FillTrackCounter(14);
      }
    } else {
      fHist->FillTrackCounter(15);
      v0->Setv0Mass(CalculateInvMass(v0,fPDGDaugN,fPDGDaugP));
      passD1=fPosCuts->isSelected(v0->GetNegDaughter());
      passD2=fNegCuts->isSelected(v0->GetPosDaughter());
      if (passD1&&passD2) {
        fHist->FillTrackCounter(16);
      }
    }
    if (!(passD1&&passD2)) {
      pass=false;
    } else {
      fHist->FillTrackCounter(2);
    }
  }
  pass=passD1&&passD2;
  return pass;
}
bool AliFemtoDreamv0Cuts::MotherPassCuts(AliFemtoDreamv0 *v0) {
  //all topological and kinematic cuts on the mother except for the CPA, this
  //will be filled later, in case also the CPA distributions are required.
  //Special CPA checks impelemented in the RejKaons later, to still ensure
  //proper Invariant Mass Plots.
  bool pass=true;
  if (fCutOnFlyStatus) {
    if (!(v0->GetOnlinev0()==fOnFlyStatus)) {
      pass=false;
    } else {
      fHist->FillTrackCounter(3);
    }
  }
  if (pass&&fCutCharge) {
    if (v0->GetCharge().at(0)!=fCharge) {
      pass=false;
    } else {
      fHist->FillTrackCounter(4);
    }
  }
  if (pass&&fCutPt) {
    if ((v0->GetPt()<fpTmin)||(v0->GetPt()>fpTmax)) {
      pass=false;
    } else {
      fHist->FillTrackCounter(5);
    }
  }
  if (pass&&fCutDecayVtxXYZ) {
    if ((v0->GetDCAv0Vtx(0)>fMaxDecayVtxXYZ)||
        (v0->GetDCAv0Vtx(1)>fMaxDecayVtxXYZ)||
        (v0->GetDCAv0Vtx(2)>fMaxDecayVtxXYZ)) {
      pass=false;
    } else {
      fHist->FillTrackCounter(6);
    }
  }
  if (pass&&fCutTransRadius) {
    if ((v0->GetTransverseRadius()<fMinTransRadius)||
        (v0->GetTransverseRadius()>fMaxTransRadius)) {
      pass=false;
    } else {
      fHist->FillTrackCounter(7);
    }
  }
  if (pass&&fCutMinDCADaugPrimVtx) {
    if ((v0->GetDCADaugPosVtx()<fMinDCADaugToPrimVtx)||
        (v0->GetDCADaugNegVtx()<fMinDCADaugToPrimVtx)) {
      pass=false;
    } else {
      fHist->FillTrackCounter(8);
    }
  }
  if (pass&&fCutMaxDCADaugToDecayVtx) {
    if (v0->GetDaugDCA()>fMaxDCADaugToDecayVtx) {
      pass=false;
    } else {
      fHist->FillTrackCounter(9);
    }
  }
  return pass;
}

bool AliFemtoDreamv0Cuts::RejectAsKaon(AliFemtoDreamv0 *v0) {
  bool pass=true;
  bool cpaPass=true;
  if (v0->GetCPA()<fMinCPA) {
    cpaPass=false;
  }
  double massKaon=CalculateInvMass(v0,211,211);
  if (cpaPass) {
    fHist->FillInvMassBefKaonRej(v0->Getv0Mass());
    fHist->FillInvMassKaon(massKaon);
  }
  if (fKaonRejLow<massKaon&&massKaon<fKaonRejUp) {
    pass=false;
  } else {
    fHist->FillTrackCounter(10);
  }
  return pass;
}

bool AliFemtoDreamv0Cuts::CPAandMassCuts(AliFemtoDreamv0 *v0) {
  //here we cut mass and cpa to fill the cpa distribution properly
  bool cpaPass=true;
  bool massPass=true;
  if (fCutCPA&&(v0->GetCPA()<fMinCPA)) {
    cpaPass=false;
  }
  if (fCutInvMass) {
    double massv0=TDatabasePDG::Instance()->GetParticle(fPDGv0)->Mass();
    if ((v0->Getv0Mass()<massv0-fInvMassCutWidth)||
        (massv0+fInvMassCutWidth<v0->Getv0Mass())) {
      massPass=false;
    }
  }
  //now with this information fill the histograms
  if (cpaPass) {
    fHist->FillInvMassPtBins(v0->GetPt(),v0->Getv0Mass());
    fHist->Fillv0MassDist(v0->Getv0Mass());
  }
  if (massPass&&fCPAPlots) {
    fHist->FillCPAPtBins(v0->GetPt(),v0->GetCPA());
    if (fMCData) {
      fMCHist->FillMCCPAPtBins(
          v0->GetParticleOrigin(),v0->GetPt(),v0->GetCPA());
    }
  }
  if (massPass) {
    fHist->FillTrackCounter(11);
  }
  if (massPass && cpaPass) {
    fHist->FillTrackCounter(12);
  }
  bool pass=massPass&&cpaPass;
  return pass;
}
void AliFemtoDreamv0Cuts::Init() {
  fHist=new AliFemtoDreamv0Hist(fNumberXBins,fAxisMinMass,fAxisMaxMass,
                                fCPAPlots);
  BookTrackCuts();
  if (fMCData) {
    fMCHist=new AliFemtoDreamv0MCHist(fNumberXBins,fAxisMinMass,fAxisMaxMass,
                                      fContribSplitting,fCPAPlots);
  }
  fPosCuts->Init();
  fNegCuts->Init();
  fHistList = new TList();
  fHistList->SetOwner();
  fHistList->SetName("v0Cuts");
  fHistList->Add(fHist->GetHistList());
  fPosCuts->SetName("PosCuts");
  fHistList->Add(fPosCuts->GetQAHists());
  fNegCuts->SetName("NegCuts");
  fHistList->Add(fNegCuts->GetQAHists());
}


void AliFemtoDreamv0Cuts::BookQA(AliFemtoDreamv0 *v0) {
  for (int i=0;i<2;++i) {
    if (i==0||(i==1&&v0->UseParticle())) {
      if (!v0->GetOnlinev0()) {
        fHist->FillOnFlyStatus(i,1);
      } else if (v0->GetOnlinev0()) {
        fHist->FillOnFlyStatus(i,0);
      }
      fHist->FillpTCut(i,v0->GetPt());
      fHist->FillEtaCut(i,v0->GetEta().at(0));
      fHist->Fillv0DecayVtxXCut(i,v0->GetDCAv0Vtx(0));
      fHist->Fillv0DecayVtxYCut(i,v0->GetDCAv0Vtx(1));
      fHist->Fillv0DecayVtxZCut(i,v0->GetDCAv0Vtx(2));
      fHist->FillTransverRadiusCut(i,v0->GetTransverseRadius());
      fHist->FillDCAPosDaugToPrimVtxCut(i,v0->GetDCADaugPosVtx());
      fHist->FillDCANegDaugToPrimVtxCut(i,v0->GetDCADaugNegVtx());
      fHist->FillDCADaugTov0VtxCut(i,v0->GetDaugDCA());
      fHist->FillCPACut(i,v0->GetCPA());
      fHist->FillInvMass(i,v0->Getv0Mass());
    }
  }
  v0->GetPosDaughter()->SetUse(v0->UseParticle());
  v0->GetNegDaughter()->SetUse(v0->UseParticle());
  if (fMCData) {
    v0->GetPosDaughter()->SetParticleOrigin(v0->GetParticleOrigin());
    v0->GetNegDaughter()->SetParticleOrigin(v0->GetParticleOrigin());
  }
  fPosCuts->BookQA(v0->GetPosDaughter());
  fNegCuts->BookQA(v0->GetNegDaughter());
}

void AliFemtoDreamv0Cuts::BookMC(AliFemtoDreamv0 *v0) {
  double pT=v0->GetPt();
  double etaNegDaug=v0->GetEta().at(1);
  double etaPosDaug=v0->GetEta().at(2);
  if (v0->GetMCPDGCode()==fPDGv0) {
    if (fpTmin<pT&&pT<fpTmax) {
      if (fPosCuts->GetEtaMin()<etaPosDaug&&etaPosDaug<fPosCuts->GetEtaMax()) {
        if (fNegCuts->GetEtaMin()<etaNegDaug&&etaNegDaug<fNegCuts->GetEtaMax()) {
          fMCHist->FillMCGen(pT);
        }
      }
    }
  }
  if (v0->UseParticle()) {
    fMCHist->FillMCIdent(pT);
    if (TMath::Abs(v0->GetMCPDGCode())==TMath::Abs(fPDGv0)) {
      fMCHist->FillMCCorr(pT);
    } else {
      v0->SetParticleOrigin(AliFemtoDreamBasePart::kContamination);
    }
    if (fContribSplitting) {
      FillMCContributions(v0);
    }
  }
}

void AliFemtoDreamv0Cuts::FillMCContributions(AliFemtoDreamv0 *v0) {
  Double_t pT = v0->GetPt();
  Int_t iFill = -1;
  switch(v0->GetParticleOrigin()){
    case AliFemtoDreamBasePart::kPhysPrimary:
      fMCHist->FillMCPrimary(pT);
      iFill = 0;
      break;
    case AliFemtoDreamBasePart::kWeak:
      fMCHist->FillMCFeeddown(pT, TMath::Abs(v0->GetMotherWeak()));
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
    default:
      AliFatal("Type Not implemented");
      break;
  }
  fMCHist->FillMCpT(iFill,pT);
  fMCHist->FillMCEta(iFill,v0->GetEta().at(0));
  fMCHist->FillMCPhi(iFill,v0->GetPhi().at(0));
  fMCHist->FillMCDCAVtxX(iFill,pT,v0->GetDCAv0Vtx(0));
  fMCHist->FillMCDCAVtxY(iFill,pT,v0->GetDCAv0Vtx(1));
  fMCHist->FillMCDCAVtxZ(iFill,pT,v0->GetDCAv0Vtx(2));
  fMCHist->FillMCTransverseRadius(iFill,pT,v0->GetTransverseRadius());
  fMCHist->FillMCDCAPosDaugPrimVtx(iFill,pT,v0->GetDCADaugPosVtx());
  fMCHist->FillMCDCANegDaugPrimVtx(iFill,pT,v0->GetDCADaugNegVtx());
  fMCHist->FillMCDCADaugVtx(iFill,pT,v0->GetDaugDCA());
  fMCHist->FillMCCosPoint(iFill,pT,v0->GetCPA());
  fMCHist->FillMCInvMass(iFill,v0->Getv0Mass());
}
void AliFemtoDreamv0Cuts::BookTrackCuts() {
  if (!fHist) {
    AliFatal("No Histograms available");
  }
  if (fCutOnFlyStatus) {
    fHist->FillConfig(0,1);
  }
  if (fCutCharge) {
    fHist->FillConfig(1,fCharge);
  }
  if (fCutPt) {
    fHist->FillConfig(2,fpTmin);
    fHist->FillConfig(3,fpTmax);
  }
  if (fKaonRejection) {
    fHist->FillConfig(4,1);
  }
  if (fCutDecayVtxXYZ) {
    fHist->FillConfig(5,fMaxDecayVtxXYZ);
  }
  if (fCutTransRadius) {
    fHist->FillConfig(6,fMinTransRadius);
    fHist->FillConfig(7,fMaxTransRadius);
  }
  if (fCutMinDCADaugPrimVtx) {
    fHist->FillConfig(8,fMinDCADaugToPrimVtx);
  }
  if (fCutMaxDCADaugToDecayVtx) {
    fHist->FillConfig(9,fMaxDCADaugToDecayVtx);
  }
  if (fCutInvMass) {
    fHist->FillConfig(10,fInvMassCutWidth);
  }
  if (fCutCPA) {
    fHist->FillConfig(11,fMinCPA);
  }
}

double AliFemtoDreamv0Cuts::CalculateInvMass(AliFemtoDreamv0 *v0,
                                             int PDGPosDaug,int PDGNegDaug) {
  Double_t invMass = 0;
  double massDP=TDatabasePDG::Instance()->GetParticle(PDGPosDaug)->Mass();
  double massDN=TDatabasePDG::Instance()->GetParticle(PDGNegDaug)->Mass();

  double EDaugP=TMath::Sqrt(
      massDP*massDP +
      v0->GetPosDaughter()->GetMomentum().X()*v0->GetPosDaughter()->GetMomentum().X()+
      v0->GetPosDaughter()->GetMomentum().Y()*v0->GetPosDaughter()->GetMomentum().Y()+
      v0->GetPosDaughter()->GetMomentum().Z()*v0->GetPosDaughter()->GetMomentum().Z());
  double EDaugN=TMath::Sqrt(
      massDN*massDN +
      v0->GetNegDaughter()->GetMomentum().X()*v0->GetNegDaughter()->GetMomentum().X()+
      v0->GetNegDaughter()->GetMomentum().Y()*v0->GetNegDaughter()->GetMomentum().Y()+
      v0->GetNegDaughter()->GetMomentum().Z()*v0->GetNegDaughter()->GetMomentum().Z());

  double energysum=EDaugP+EDaugN;
  double pSum2=
      (v0->GetNegDaughter()->GetMomentum().X()+v0->GetPosDaughter()->GetMomentum().X())*
      (v0->GetNegDaughter()->GetMomentum().X()+v0->GetPosDaughter()->GetMomentum().X())+

      (v0->GetNegDaughter()->GetMomentum().Y()+v0->GetPosDaughter()->GetMomentum().Y())*
      (v0->GetNegDaughter()->GetMomentum().Y()+v0->GetPosDaughter()->GetMomentum().Y())+

      (v0->GetNegDaughter()->GetMomentum().Z()+v0->GetPosDaughter()->GetMomentum().Z())*
      (v0->GetNegDaughter()->GetMomentum().Z()+v0->GetPosDaughter()->GetMomentum().Z());
  invMass = TMath::Sqrt(energysum*energysum - pSum2);
  return invMass;
}
