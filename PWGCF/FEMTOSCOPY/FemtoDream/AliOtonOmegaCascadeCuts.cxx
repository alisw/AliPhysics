/*
 * AliOtonOmegaCascadeCuts.cxx
 *
 *  Created on: Jan 11, 2018
 *      Author: gu74req
 */

#include "vector"
#include "AliOtonOmegaCascadeCuts.h"
#include "AliLog.h"
ClassImp(AliOtonOmegaCascadeCuts)
AliOtonOmegaCascadeCuts::AliOtonOmegaCascadeCuts()
    : fHist(nullptr),
      fMCHist(nullptr),
      fNegCuts(0),
      fPosCuts(0),
      fBachCuts(0),
      fHistList(0),
      fMCHistList(0),
      fMinimalBooking(false),
      fMCData(false),
      fContribSplitting(false),
      fCutFlag(0),
      fRunNumberQA(false),
      fMinRunNumber(0),
      fMaxRunNumber(0),
      fCutPt(false),
      fPtMin(0),
      fPtMax(0),
      fCutPtv0(false),
      fPtMinv0(0),
      fPtMaxv0(0),
      fcutXiMass(false),
      fXiMass(0),
      fXiMassWidth(0),
      fXiMassExcl(0),
      fXiMassWidthExcl(0),
      fcutXiCharge(false),
      fXiCharge(0),
      fcutDCAXiDaug(false),
      fMaxDCAXiDaug(0),
      fcutMinDistVtxBach(false),
      fMinDistVtxBach(0),
      fcutCPAXi(false),
      fCPAXi(0),
      fcutXiTransRadius(false),
      fMinXiTransRadius(0),
      fMaxXiTransRadius(0),
      fcutv0Mass(false),
      fv0Mass(0),
      fv0Width(0),
      fcutv0MaxDCADaug(false),
      fv0MaxDCADaug(0),
      fcutCPAv0(false),
      fCPAv0(0),
      fcutv0TransRadius(false),
      fMinv0TransRadius(0),
      fMaxv0TransRadius(0),
      fcutv0MinDistVtx(false),
      fv0MinDistVtx(0),
      fcutv0DaugMinDistVtx(false),
      fv0DaugMinDistVtx(0),
      fRejMassCut(false),
      fRejPDG(0),
      fRejMass(0),
      fRejWidth(0),
      fPDGCasc(0),
      fPDGv0(0),
      fPDGPosDaug(0),
      fPDGNegDaug(0),
      fPDGBachDaug(0) {
}

AliOtonOmegaCascadeCuts::AliOtonOmegaCascadeCuts(
    const AliOtonOmegaCascadeCuts& cuts)
    : fHist(cuts.fHist),
      fMCHist(cuts.fMCHist),
      fNegCuts(cuts.fNegCuts),
      fPosCuts(cuts.fPosCuts),
      fBachCuts(cuts.fBachCuts),
      fHistList(cuts.fHistList),
      fMCHistList(cuts.fMCHistList),
      fMinimalBooking(cuts.fMinimalBooking),
      fMCData(cuts.fMCData),
      fContribSplitting(cuts.fContribSplitting),
      fCutFlag(cuts.fCutFlag),
      fRunNumberQA(cuts.fRunNumberQA),
      fMinRunNumber(cuts.fMinRunNumber),
      fMaxRunNumber(cuts.fMaxRunNumber),
      fCutPt(cuts.fCutPt),
      fPtMin(cuts.fPtMin),
      fPtMax(cuts.fPtMax),
      fCutPtv0(cuts.fCutPtv0),
      fPtMinv0(cuts.fPtMinv0),
      fPtMaxv0(cuts.fPtMaxv0),
      fcutXiMass(cuts.fcutXiMass),
      fXiMass(cuts.fXiMass),
      fXiMassWidth(cuts.fXiMassWidth),
      fXiMassExcl(cuts.fXiMassExcl),
      fXiMassWidthExcl(cuts.fXiMassWidthExcl),
      fcutXiCharge(cuts.fcutXiCharge),
      fXiCharge(cuts.fXiCharge),
      fcutDCAXiDaug(cuts.fcutDCAXiDaug),
      fMaxDCAXiDaug(cuts.fMaxDCAXiDaug),
      fcutMinDistVtxBach(cuts.fcutMinDistVtxBach),
      fMinDistVtxBach(cuts.fMinDistVtxBach),
      fcutCPAXi(cuts.fcutCPAXi),
      fCPAXi(cuts.fCPAXi),
      fcutXiTransRadius(cuts.fcutXiTransRadius),
      fMinXiTransRadius(cuts.fMinXiTransRadius),
      fMaxXiTransRadius(cuts.fMaxXiTransRadius),
      fcutv0Mass(cuts.fcutv0Mass),
      fv0Mass(cuts.fv0Mass),
      fv0Width(cuts.fv0Width),
      fcutv0MaxDCADaug(cuts.fcutv0MaxDCADaug),
      fv0MaxDCADaug(cuts.fv0MaxDCADaug),
      fcutCPAv0(cuts.fcutCPAv0),
      fCPAv0(cuts.fCPAv0),
      fcutv0TransRadius(cuts.fcutv0TransRadius),
      fMinv0TransRadius(cuts.fMinv0TransRadius),
      fMaxv0TransRadius(cuts.fMaxv0TransRadius),
      fcutv0MinDistVtx(cuts.fcutv0MinDistVtx),
      fv0MinDistVtx(cuts.fv0MinDistVtx),
      fcutv0DaugMinDistVtx(cuts.fcutv0DaugMinDistVtx),
      fv0DaugMinDistVtx(cuts.fv0DaugMinDistVtx),
      fRejMassCut(cuts.fRejMassCut),
      fRejPDG(cuts.fRejPDG),
      fRejMass(cuts.fRejMass),
      fRejWidth(cuts.fRejWidth),
      fPDGCasc(cuts.fPDGCasc),
      fPDGv0(cuts.fPDGv0),
      fPDGPosDaug(cuts.fPDGPosDaug),
      fPDGNegDaug(cuts.fPDGNegDaug),
      fPDGBachDaug(cuts.fPDGBachDaug) {
}

AliOtonOmegaCascadeCuts& AliOtonOmegaCascadeCuts::operator=(
    const AliOtonOmegaCascadeCuts& cuts) {
  if (this == &cuts) {
    return *this;
  }

  this->fHist = cuts.fHist;
  this->fMCHist = cuts.fMCHist;
  this->fNegCuts = cuts.fNegCuts;
  this->fPosCuts = cuts.fPosCuts;
  this->fBachCuts = cuts.fBachCuts;
  this->fHistList = cuts.fHistList;
  this->fMCHistList = cuts.fMCHistList;
  this->fMinimalBooking = cuts.fMinimalBooking;
  this->fMCData = cuts.fMCData;
  this->fContribSplitting = cuts.fContribSplitting;
  this->fCutFlag = cuts.fCutFlag;
  this->fRunNumberQA = cuts.fRunNumberQA;
  this->fMinRunNumber = cuts.fMinRunNumber;
  this->fMaxRunNumber = cuts.fMaxRunNumber;
  this->fCutPt = cuts.fCutPt;
  this->fPtMin = cuts.fPtMin;
  this->fPtMax = cuts.fPtMax;
  this->fCutPtv0 = cuts.fCutPtv0;
  this->fPtMinv0 = cuts.fPtMinv0;
  this->fPtMaxv0 = cuts.fPtMaxv0;
  this->fcutXiMass = cuts.fcutXiMass;
  this->fXiMass = cuts.fXiMass;
  this->fXiMassWidth = cuts.fXiMassWidth;
  this->fXiMassExcl = cuts.fXiMassExcl;
  this->fXiMassWidthExcl = cuts.fXiMassWidthExcl;
  this->fcutXiCharge = cuts.fcutXiCharge;
  this->fXiCharge = cuts.fXiCharge;
  this->fcutDCAXiDaug = cuts.fcutDCAXiDaug;
  this->fMaxDCAXiDaug = cuts.fMaxDCAXiDaug;
  this->fcutMinDistVtxBach = cuts.fcutMinDistVtxBach;
  this->fMinDistVtxBach = cuts.fMinDistVtxBach;
  this->fcutCPAXi = cuts.fcutCPAXi;
  this->fCPAXi = cuts.fCPAXi;
  this->fcutXiTransRadius = cuts.fcutXiTransRadius;
  this->fMinXiTransRadius = cuts.fMinXiTransRadius;
  this->fMaxXiTransRadius = cuts.fMaxXiTransRadius;
  this->fcutv0Mass = cuts.fcutv0Mass;
  this->fv0Mass = cuts.fv0Mass;
  this->fv0Width = cuts.fv0Width;
  this->fcutv0MaxDCADaug = cuts.fcutv0MaxDCADaug;
  this->fv0MaxDCADaug = cuts.fv0MaxDCADaug;
  this->fcutCPAv0 = cuts.fcutCPAv0;
  this->fCPAv0 = cuts.fCPAv0;
  this->fcutv0TransRadius = cuts.fcutv0TransRadius;
  this->fMinv0TransRadius = cuts.fMinv0TransRadius;
  this->fMaxv0TransRadius = cuts.fMaxv0TransRadius;
  this->fcutv0MinDistVtx = cuts.fcutv0MinDistVtx;
  this->fv0MinDistVtx = cuts.fv0MinDistVtx;
  this->fcutv0DaugMinDistVtx = cuts.fcutv0DaugMinDistVtx;
  this->fv0DaugMinDistVtx = cuts.fv0DaugMinDistVtx;
  this->fRejMassCut = cuts.fRejMassCut;
  this->fRejPDG = cuts.fRejPDG;
  this->fRejMass = cuts.fRejMass;
  this->fRejWidth = cuts.fRejWidth;
  this->fPDGCasc = cuts.fPDGCasc;
  this->fPDGv0 = cuts.fPDGv0;
  this->fPDGPosDaug = cuts.fPDGPosDaug;
  this->fPDGNegDaug = cuts.fPDGNegDaug;
  this->fPDGBachDaug = cuts.fPDGBachDaug;

  return *this;
}

AliOtonOmegaCascadeCuts::~AliOtonOmegaCascadeCuts() {
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

AliOtonOmegaCascadeCuts* AliOtonOmegaCascadeCuts::XiCuts(
    bool isMC, bool contribSplitting, Int_t CutFlag) {
  AliOtonOmegaCascadeCuts *XiCuts = new AliOtonOmegaCascadeCuts();
  XiCuts->SetIsMonteCarlo(isMC);
  XiCuts->SetContributionSplitting(contribSplitting);
  
  if(CutFlag == 1){ // Cuts with TIME info for 2 tracks: proton and bachelor (its or tof time)
                    // sigmas TPC 5 for proton, pion, and bachelor
   XiCuts->SetXiMassRange(1.322, 0.005, 0., 0.);	//YES
   XiCuts->SetCutXiDaughterDCA(.9);		//YES
   XiCuts->SetCutXiMinDistBachToPrimVtx(0.03);	//as ESD new vertexer cut
   XiCuts->SetCutXiCPA(0.991);			//YES
   XiCuts->SetCutXiTransverseRadius(0.4, 200);	//as ESD new vertexer cut
   XiCuts->Setv0MassRange(1.116, 0.006);  	//YES
   XiCuts->SetCutv0MaxDaughterDCA(2.);   	//ESD new vertexer cut
   XiCuts->SetCutv0CPA(0.99);			//YES
   XiCuts->SetCutv0TransverseRadius(1., 200);	//as ESD new vertexer cut
   XiCuts->SetCutv0MinDistToPrimVtx(0.05);   	//as ESD new vertexer cut
   XiCuts->SetCutv0MinDaugDistToPrimVtx(0.02);	//as ESD new vertexer cut
   XiCuts->SetRejectMass(1.672, 0.005, 3334);	//YES
   XiCuts->SetPtRangeXi(0., 999.9);   		//as ESD new vertexer cut
  }else if(CutFlag == 2){ // Cuts with TIME info for 3 tracks: proton and bachelor (its or tof time)
                    // sigmas TPC 5 for proton, pion, and bachelor
   XiCuts->SetXiMassRange(1.322, 0.005, 0., 0.);        //YES
   XiCuts->SetCutXiDaughterDCA(.9);             //YES
   XiCuts->SetCutXiMinDistBachToPrimVtx(0.03);  //as ESD new vertexer cut
   XiCuts->SetCutXiCPA(0.991);                  //YES
   XiCuts->SetCutXiTransverseRadius(0.4, 200);  //as ESD new vertexer cut
   XiCuts->Setv0MassRange(1.116, 0.006);        //YES   
   XiCuts->SetCutv0MaxDaughterDCA(2.);         	//as ESD new vertexer cut
   XiCuts->SetCutv0CPA(0.99);                   //YES
   XiCuts->SetCutv0TransverseRadius(1., 200);  	//as ESD new vertexer cut
   XiCuts->SetCutv0MinDistToPrimVtx(0.05);      //as ESD new vertexer cut
   XiCuts->SetCutv0MinDaugDistToPrimVtx(0.02);  //as ESD new vertexer cut
   XiCuts->SetRejectMass(1.672, 0.005, 3334);   //YES
   XiCuts->SetPtRangeXi(0., 999.9);            	//as ESD new vertexer cut
  }else{ // CutFlag == 0) // Bernie&Andi cuts
   XiCuts->SetXiMassRange(1.322, 0.005, 0., 0.);
   XiCuts->SetCutXiDaughterDCA(1.6);
   XiCuts->SetCutXiMinDistBachToPrimVtx(0.05);
   XiCuts->SetCutXiCPA(0.98);
   XiCuts->SetCutXiTransverseRadius(0.8, 200);
   XiCuts->Setv0MassRange(1.116, 0.006);
   XiCuts->SetCutv0MaxDaughterDCA(1.5);
   XiCuts->SetCutv0CPA(0.97);
   XiCuts->SetCutv0TransverseRadius(1.4, 200);
   XiCuts->SetCutv0MinDistToPrimVtx(0.07);
   XiCuts->SetCutv0MinDaugDistToPrimVtx(0.05);
   XiCuts->SetRejectMass(1.672, 0.005, 3334);
   XiCuts->SetPtRangeXi(0.3, 999.9);
  }
  return XiCuts;
}

AliOtonOmegaCascadeCuts* AliOtonOmegaCascadeCuts::OmegaCuts(
    bool isMC, bool contribSplitting, Int_t CutFlag) {
  AliOtonOmegaCascadeCuts *OmegaCuts = new AliOtonOmegaCascadeCuts();
  OmegaCuts->SetIsMonteCarlo(isMC);
  OmegaCuts->SetContributionSplitting(contribSplitting);
  if(CutFlag == 1){ // Cuts with TIME info for 2 tracks: proton and bachelor (its or tof time)
                    // sigmas TPC 5 for proton, pion, and bachelor
   OmegaCuts->SetXiMassRange(1.672, 0.005, 0., 0.);		//YES
   OmegaCuts->SetCutXiDaughterDCA(.8);			//YES
   OmegaCuts->SetCutXiMinDistBachToPrimVtx(0.03);	//as ESD new vertexer cut
   OmegaCuts->SetCutXiCPA(0.995);			//YES
   OmegaCuts->SetCutXiTransverseRadius(0.4, 200);       //as ESD new vertexer cut
   OmegaCuts->Setv0MassRange(1.116, 0.006);		//YES
   OmegaCuts->SetCutv0MaxDaughterDCA(1.2);		//YES
   OmegaCuts->SetCutv0CPA(0.99);			//YES
   OmegaCuts->SetCutv0TransverseRadius(1., 200);       	//as ESD new vertexer cut
   OmegaCuts->SetCutv0MinDistToPrimVtx(0.06);		//YES
   OmegaCuts->SetCutv0MinDaugDistToPrimVtx(0.02);       //as ESD new vertexer cut
   OmegaCuts->SetRejectMass(1.322, 0.013, 3312);	//YES...but should be changed to asymmetric from 1.322+0.013
   OmegaCuts->SetPtRangeXi(0., 999.9);       		//as ESD new vertexer cut
  }else if(CutFlag == 2){ // Cuts with TIME info for 3 tracks: proton and bachelor (its or tof time)
                          // sigmas TPC 6 for proton, pion, and bachelor || SIX SIGMAS HERE!!!!
   OmegaCuts->SetXiMassRange(1.672, 0.005, 0., 0.);             //YES
   OmegaCuts->SetCutXiDaughterDCA(.8);                  //YES
   OmegaCuts->SetCutXiMinDistBachToPrimVtx(0.03);       //as ESD new vertexer cut
   OmegaCuts->SetCutXiCPA(0.995);                       //YES
   OmegaCuts->SetCutXiTransverseRadius(0.4, 200);       //as ESD new vertexer cut
   OmegaCuts->Setv0MassRange(1.116, 0.006);             //YES
   OmegaCuts->SetCutv0MaxDaughterDCA(1.2);              //YES
   OmegaCuts->SetCutv0CPA(0.99);                        //YES
   OmegaCuts->SetCutv0TransverseRadius(1., 200);       	//as ESD new vertexer cut
   OmegaCuts->SetCutv0MinDistToPrimVtx(0.05);           //YES...but I think this ESD new vertexer cut
   OmegaCuts->SetCutv0MinDaugDistToPrimVtx(0.02);       //as ESD new vertexer cut
   OmegaCuts->SetRejectMass(1.322, 0.013, 3312);        //YES...but should be changed to asymmetric from 1.322+0.013
   OmegaCuts->SetPtRangeXi(0., 999.9);                 	//as ESD new vertexer cut
  }else{ // CutFlag == 0) // Bernie&Andi cuts
   OmegaCuts->SetXiMassRange(1.672, 0.005, 0., 0.);
   OmegaCuts->SetCutXiDaughterDCA(2.);
   OmegaCuts->SetCutXiMinDistBachToPrimVtx(0.04);
   OmegaCuts->SetCutXiCPA(0.98);
   OmegaCuts->SetCutXiTransverseRadius(0.5, 200);
   OmegaCuts->Setv0MassRange(1.116, 0.006);
   OmegaCuts->SetCutv0MaxDaughterDCA(1.5);
   OmegaCuts->SetCutv0CPA(0.97);
   OmegaCuts->SetCutv0TransverseRadius(1.1, 200);
   OmegaCuts->SetCutv0MinDistToPrimVtx(0.06);
   OmegaCuts->SetCutv0MinDaugDistToPrimVtx(0.04);
   OmegaCuts->SetRejectMass(1.322, 0.005, 3312);
   OmegaCuts->SetPtRangeXi(0.3, 999.9);
  }
  return OmegaCuts;
}

bool AliOtonOmegaCascadeCuts::isSelected(AliOtonOmegaCascade *casc) {
  bool pass = true;
  if (casc->IsSet()) {
    if (!fMinimalBooking)
      fHist->FillCutCounter(0);
    if (fcutXiCharge) {
      if (!(casc->GetCharge().at(0) == fXiCharge)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(1);
      }
    }
    if (pass) {
      std::vector<int> IDTracks = casc->GetIDTracks();
      if (IDTracks[0] == IDTracks[1]) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(2);
      }
      if (IDTracks[0] == IDTracks[2]) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(3);
      }
      if (IDTracks[1] == IDTracks[2]) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(4);
      }
    }
    if (pass) {
      if (!fNegCuts->isSelected(casc->GetNegDaug())) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(5);
      }
    }
    if (pass) {
      if (!fPosCuts->isSelected(casc->GetPosDaug())) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(6);
      }
    }
    if (pass && fcutv0MaxDCADaug) {
      if (casc->Getv0DCADaug() > fv0MaxDCADaug) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(7);
      }
    }
    if (pass && fcutCPAv0) {
      if (casc->Getv0CPA() < fCPAv0) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(8);
      }
    }
    if (pass && fCutPtv0) {
      if ((casc->Getv0Pt() < fPtMinv0) || (fPtMaxv0 < casc->Getv0Pt())) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(9);
      }
    }
    if (pass && fcutv0TransRadius) {
      if ((casc->Getv0TransverseRadius() < fMinv0TransRadius)
          || (casc->Getv0TransverseRadius() > fMaxv0TransRadius)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(10);
      }
    }
    if (pass && fcutv0MinDistVtx) {
      if (casc->Getv0DCAPrimVtx() < fv0MinDistVtx) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(11);
      }
    }
    if (pass && fcutv0DaugMinDistVtx) {
      if ((casc->Getv0PosToPrimVtx() < fv0DaugMinDistVtx)
          || (casc->Getv0NegToPrimVtx() < fv0DaugMinDistVtx)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(12);
      }
    }
    if (pass) {
      if (!fMinimalBooking)
        fHist->FillInvMassPtv0(casc->GetPt(), casc->Getv0Mass());
    }
    if (pass && fcutv0Mass) {
      if ((casc->Getv0Mass() < (fv0Mass - fv0Width))
          || (casc->Getv0Mass() > (fv0Mass + fv0Width))) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(13);
      }
    }
    if (pass) {
      if (!fBachCuts->isSelected(casc->GetBach())) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(14);
      }
    }
    if (pass && fcutDCAXiDaug) {
      if (casc->GetXiDCADaug() > fMaxDCAXiDaug) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(15);
      }
    }
    if (pass && fcutMinDistVtxBach) {
      if (casc->BachDCAPrimVtx() < fMinDistVtxBach) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(16);
      }
    }
    if (pass && fcutCPAXi) {
      if (casc->GetCPA() < fCPAXi) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(17);
      }
    }
    if (pass && fcutXiTransRadius) {
      if ((casc->GetXiTransverseRadius() < fMinXiTransRadius)
          || (casc->GetXiTransverseRadius() > fMaxXiTransRadius)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(18);
      }
    }
    if (pass && fRejMassCut) {
      if ((casc->GetMass(fRejPDG) > (fRejMass - fRejWidth))
          && (casc->GetMass(fRejPDG) < (fRejMass + fRejWidth))) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(19);
      }
    }
    if (pass && fCutPt) {
      if ((casc->GetPt() < fPtMin) || (fPtMax < casc->GetPt())) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(20);
      }
    }
    if (pass) {
      fHist->FillInvMassPt(casc->GetPt(), casc->GetMass(fPDGCasc));
      if (fRunNumberQA) {
        fHist->FillInvMassPerRunNumber(casc->GetEvtNumber(), casc->GetMass(fPDGCasc));
      }
    }
//Kaon missidentification background removal cut, associated to the standard selection with XiCPA<0.995 cut
    if (pass && fCPAXi == 0.995  ) {
     //Fill some variables:
     Float_t pipx,pipy,pipz;
     Float_t prpx,prpy,prpz;
     Int_t piNcl;
     if(casc->GetBach()->GetCharge().at(0)==-1) {
      prpx = casc->GetPosDaug()->GetMomentum().Px();
      prpy = casc->GetPosDaug()->GetMomentum().Py();
      prpz = casc->GetPosDaug()->GetMomentum().Pz();
      pipx = casc->GetNegDaug()->GetMomentum().Px();
      pipy = casc->GetNegDaug()->GetMomentum().Py();
      pipz = casc->GetNegDaug()->GetMomentum().Pz();
      piNcl = casc->GetNegDaug()->GetNClsTPC();
     }else{
      pipx = casc->GetPosDaug()->GetMomentum().Px();
      pipy = casc->GetPosDaug()->GetMomentum().Py();
      pipz = casc->GetPosDaug()->GetMomentum().Pz();
      piNcl = casc->GetPosDaug()->GetNClsTPC();
      prpx = casc->GetNegDaug()->GetMomentum().Px();
      prpy = casc->GetNegDaug()->GetMomentum().Py();
      prpz = casc->GetNegDaug()->GetMomentum().Pz();
     }
     Float_t op = casc->GetMomentum().Mag();//omega momentum
     Float_t lpx = pipx+prpx;
     Float_t lpy = pipy+prpy;
     Float_t lpz = pipz+prpz;
     Float_t lp = sqrt(pow(lpx,2)+pow(lpy,2)+pow(lpz,2));
     Float_t pip = sqrt(pow(pipx,2)+pow(pipy,2)+pow(pipz,2));
     Float_t prp = sqrt(pow(prpx,2)+pow(prpy,2)+pow(prpz,2));
     Float_t kpx = casc->GetBach()->GetMomentum().Px();
     Float_t kpy = casc->GetBach()->GetMomentum().Py();
     Float_t kpz = casc->GetBach()->GetMomentum().Pz();
     Float_t kp = sqrt(pow(kpx,2)+pow(kpy,2)+pow(kpz,2));
     Float_t kPIDk = casc->GetBach()->GetnSigmaTPC(3);
     Float_t lr = casc->Getv0TransverseRadius();
     Float_t pi1E = sqrt(pow(0.13957018,2)+pow(pip,2));
     Float_t pi2E = sqrt(pow(0.13957018,2)+pow(prp,2));
     Float_t k0mass = sqrt(pow(pi1E+pi2E,2)-pow(lp,2));
     //do the actual cut inspired by ntuple cut:
     //momNEW = "op>2.&&op<4.&&lp>1.4&&lp<3.&&pip<.24&&kp>.56&&kp<1.24&&k0>.47&&kPIDk<2&&piNcl<150";
     //!momNEW&&"lr<100&&kPIDk>-2.5&&piNcl>60
     Bool_t isBkg = kFALSE;
     if(op>2.&&op<4.&&lp>1.4&&lp<3.&&pip<.24&&kp>.56&&kp<1.24&&k0mass>.47&&kPIDk<2&&piNcl<150) isBkg = kTRUE;
     Bool_t isSuperGoodOmega = kFALSE;
     if(lr<100&&kPIDk>-2.5&&piNcl>60) isSuperGoodOmega = kTRUE;
     if(isBkg) pass = false;
     if(!isSuperGoodOmega) pass = false;
    }
//
    if (pass && fcutXiMass) {
      if(
           (casc->GetMass(fPDGCasc) < (fXiMass - fXiMassWidth))||(casc->GetMass(fPDGCasc) > (fXiMass + fXiMassWidth))
           ||(fXiMassExcl>0.&&fXiMassWidthExcl>0.&&(  (casc->GetMass(fPDGCasc) > (fXiMassExcl - fXiMassWidthExcl))&&(casc->GetMass(fPDGCasc) < (fXiMassExcl + fXiMassWidthExcl))  ))    
      ){
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillCutCounter(21);
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
  } else {
    pass = false;
    if (!fMinimalBooking)
      fHist->FillCutCounter(22);
  }
  return pass;
}

void AliOtonOmegaCascadeCuts::Init() {
  fNegCuts->SetMinimalBooking(fMinimalBooking);
  fNegCuts->Init();
  fPosCuts->SetMinimalBooking(fMinimalBooking);
  fPosCuts->Init();
  fBachCuts->SetMinimalBooking(fMinimalBooking);
  fBachCuts->Init();
  if (!fMinimalBooking) {
    fHist = new AliFemtoDreamCascadeHist(fXiMass, fRunNumberQA, fMinRunNumber,
                                         fMaxRunNumber);
    BookCuts();
    if (!(fNegCuts || fPosCuts || fBachCuts)) {
      AliFatal("Track Cuts Object Missing");
    }
    fHistList = new TList();
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
      fMCHist = new AliFemtoDreamv0MCHist(400, 1.0, 1.2, fContribSplitting,
                                          false);
      fMCHistList = new TList();
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
    fHist = new AliFemtoDreamCascadeHist("MinimalBooking", fXiMass);
    fHistList = new TList();
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

void AliOtonOmegaCascadeCuts::BookQA(AliOtonOmegaCascade *casc) {
  if (!fMinimalBooking) {
    for (int i = 0; i < 2; ++i) {
      if (i == 0 || (i == 1 && casc->UseParticle())) {
        fHist->FillInvMass(i, casc->GetMass(fPDGCasc));
        fHist->FillInvMassLambda(i, casc->Getv0Mass());
        fHist->FillXiPt(i, casc->GetMomentum().Pt());
        fHist->FillMomRapXi(i, casc->GetXiRapidity(),
                            casc->GetMomentum().Mag());
        fHist->FillMomRapOmega(i, casc->GetOmegaRapidity(),
                               casc->GetMomentum().Mag());
        fHist->FillDCAXiPrimVtx(i, casc->GetDCAXiPrimVtx());
        fHist->FillDCAXiDaug(i, casc->GetXiDCADaug());
        fHist->FillMinDistPrimVtxBach(i, casc->BachDCAPrimVtx());
        fHist->FillCPAXi(i, casc->GetCPA());
        fHist->FillDecayLength(i, casc->GetDecayLength());
        fHist->Fillv0DecayLength(i, casc->Getv0DecayLength());
        fHist->FillTransverseRadiusXi(i, casc->GetXiTransverseRadius());
        fHist->FillMaxDCAv0Daug(i, casc->Getv0DCADaug());
        fHist->FillCPAv0(i, casc->Getv0CPA());
        fHist->FillCPAv0Xi(i, casc->Getv0XiPointingAngle());
        fHist->Fillv0Pt(i, casc->Getv0Pt());
        fHist->FillTransverseRadiusv0(i, casc->Getv0TransverseRadius());
        fHist->FillMinDistPrimVtxv0(i, casc->Getv0DCAPrimVtx());
        fHist->FillMinDistPrimVtxv0DaugPos(i, casc->Getv0PosToPrimVtx());
        fHist->FillMinDistPrimVtxv0DaugNeg(i, casc->Getv0NegToPrimVtx());
        fHist->FillPodolandski(i, casc->GetXiAlpha(), casc->GetPtArmXi());
      }
    }
  }
  if (casc->GetNegDaug()->IsSet())
    fNegCuts->BookQA(casc->GetNegDaug());
  if (casc->GetPosDaug()->IsSet())
    fPosCuts->BookQA(casc->GetPosDaug());
  if (casc->GetBach()->IsSet())
    fBachCuts->BookQA(casc->GetBach());
  return;
}

void AliOtonOmegaCascadeCuts::BookMCQA(AliOtonOmegaCascade *casc) {
  if (!fMinimalBooking) {
    float pT = casc->GetPt();
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
      fMCHist->FillMCIdent(pT);
      AliFemtoDreamBasePart::PartOrigin tmpOrg = casc->GetParticleOrigin();
      if (casc->GetParticleOrigin() != AliFemtoDreamBasePart::kFake) {
        if (casc->GetMCPDGCode() == fPDGCasc) {
          if (tmpOrg == AliFemtoDreamBasePart::kPhysPrimary) {
            fMCHist->FillMCCorr(pT);
            fMCHist->FillMCPtResolution(casc->GetMCPt(), casc->GetPt());
            fMCHist->FillMCThetaResolution(casc->GetMCTheta().at(0),
                                           casc->GetTheta().at(0),
                                           casc->GetMCPt());
            fMCHist->FillMCPhiResolution(casc->GetMCPhi().at(0),
                                         casc->GetPhi().at(0), casc->GetMCPt());
          }
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

void AliOtonOmegaCascadeCuts::FillMCContributions(AliOtonOmegaCascade *casc) {
  if (!fMinimalBooking) {
    Double_t pT = casc->GetPt();
    Int_t iFill = -1;
    switch (casc->GetParticleOrigin()) {
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
    if (iFill != -1) {
      fMCHist->FillMCpT(iFill, pT);
      fMCHist->FillMCEta(iFill, casc->GetEta().at(0));
      fMCHist->FillMCPhi(iFill, casc->GetPhi().at(0));
      fMCHist->FillMCPodolanski(iFill, casc->GetPtArmXi(), casc->GetXiAlpha());
      fMCHist->FillMCCosPoint(iFill, pT, casc->GetCPA());
      fMCHist->FillMCBachDCAToPV(iFill, pT, casc->BachDCAPrimVtx());
      fMCHist->FillMCv0DecayLength(iFill, pT, casc->Getv0DecayLength());
      fMCHist->FillMCv0CPA(iFill, pT, casc->Getv0CPA());
      fMCHist->FillMCDecayLength(iFill, pT, casc->GetDecayLength());
      fMCHist->FillMCXiRapidity(iFill, casc->GetMomentum().Mag(),
                                casc->GetXiRapidity());
      fMCHist->FillMCOmegaRapidity(iFill, casc->GetMomentum().Mag(),
                                   casc->GetOmegaRapidity());
      fMCHist->FillMCTransverseRadius(iFill, pT, casc->GetXiTransverseRadius());
      fMCHist->FillMCDCAPosDaugPrimVtx(iFill, pT, casc->Getv0PosToPrimVtx());
      fMCHist->FillMCDCANegDaugPrimVtx(iFill, pT, casc->Getv0NegToPrimVtx());
      fMCHist->FillMCDCADaugVtx(iFill, pT, casc->GetXiDCADaug());
      fMCHist->FillMCInvMass(iFill, casc->Getv0Mass());
      fMCHist->FillMCXiInvMass(iFill, pT, casc->GetXiMass());
      fMCHist->FillMCOmegaInvMass(iFill, pT, casc->GetOmegaMass());

    } else {
      std::cout << "this should not happen \n";
    }
  }
}

void AliOtonOmegaCascadeCuts::BookCuts() {
  if (fCutPt) {
    fHist->FillConfig(0, fPtMin);
    fHist->FillConfig(1, fPtMax);
  }
  if (fCutPtv0) {
    fHist->FillConfig(2, fPtMinv0);
    fHist->FillConfig(3, fPtMaxv0);
  }
  if (fcutXiMass) {
    fHist->FillConfig(4, fXiMass);
    fHist->FillConfig(5, fXiMassWidth);
  }
  if (fcutXiCharge) {
    fHist->FillConfig(6, fXiCharge);
  }
  if (fcutDCAXiDaug) {
    fHist->FillConfig(7, fMaxDCAXiDaug);
  }
  if (fcutMinDistVtxBach) {
    fHist->FillConfig(8, fMinDistVtxBach);
  }
  if (fcutCPAXi) {
    fHist->FillConfig(9, fCPAXi);
  }
  if (fcutXiTransRadius) {
    fHist->FillConfig(10, fMinXiTransRadius);
    fHist->FillConfig(11, fMaxXiTransRadius);
  }
  if (fcutv0Mass) {
    fHist->FillConfig(12, fv0Mass);
    fHist->FillConfig(13, fv0Width);
  }
  if (fcutv0MaxDCADaug) {
    fHist->FillConfig(14, fv0MaxDCADaug);
  }
  if (fcutCPAv0) {
    fHist->FillConfig(15, fCPAv0);
  }
  if (fcutv0TransRadius) {
    fHist->FillConfig(16, fMinv0TransRadius);
    fHist->FillConfig(17, fMaxv0TransRadius);
  }
  if (fcutv0MinDistVtx) {
    fHist->FillConfig(18, fv0MinDistVtx);
  }
  if (fcutv0DaugMinDistVtx) {
    fHist->FillConfig(19, fv0DaugMinDistVtx);
  }
  if (fRejMassCut) {
    fHist->FillConfig(20, fRejMass);
    fHist->FillConfig(21, fRejWidth);
  }
}

