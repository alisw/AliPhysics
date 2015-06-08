/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Class for storing and handling D0 meson candidates properties         //
//  for estimating the feed-down fraction using several sets of cuts      //
//     Andrea Rossi <andrea.rossi@cern.ch>                                //
//     Felix Reidt  <felix.reidt@cern.ch>                                 //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
#include <TNamed.h>
#include <THnSparse.h>
#include <TH1F.h>
#include <TClonesArray.h>
#include "AliAODRecoDecayHF2Prong.h"
#include "AliHFsubtractBFDcuts.h"

AliHFsubtractBFDcuts::AliHFsubtractBFDcuts() :
  TNamed(),
  fIsMC(kTRUE),
  fPtMCGenStep(0x0),
  fCutsData(0x0),
  fCutsMC(0x0)
{
  // default constructor
}

AliHFsubtractBFDcuts::AliHFsubtractBFDcuts(const char* name, const char* title) :
  TNamed(name,title),
  fIsMC(kTRUE),
  fPtMCGenStep(0x0),
  fCutsData(0x0),
  fCutsMC(0x0)
{
  // default constructor with name
}

void AliHFsubtractBFDcuts::InitHistos(){
  Int_t dimAxes[4]={500  ,24 ,30 ,400   };// mass, pt, normLXY, cosPointXY
  Double_t  min[4]={  1.7, 0., 0.,  0.96};
  Double_t  max[4]={  2.2,24.,30.,  1.  };
  fCutsData=new THnSparseF("fCutsDataFD","fCutsDataFD",4,dimAxes,min,max);
  fCutsData->GetAxis(0)->SetName("mass");
  fCutsData->GetAxis(0)->SetTitle("Mass (K,#pi) (GeV/#it{c^{2}})");
  fCutsData->GetAxis(1)->SetName("pt");
  fCutsData->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fCutsData->GetAxis(2)->SetName("NormDecLengthXY");
  fCutsData->GetAxis(2)->SetTitle("Normalized XY decay length");
  fCutsData->GetAxis(3)->SetName("CosPointXY");
  fCutsData->GetAxis(3)->SetTitle("Cos#theta_{point}^{XY}");

  Int_t dimAxesMC[5]={24 , 30,400   ,20 ,24 };// pt, normLXY, cosPointXY, #prongs, mother pt
  Double_t  minMC[5]={ 0., 0.,  0.96, 0., 0.};
  Double_t  maxMC[5]={24.,30.,  1.  ,20.,24.};
  fCutsMC=new THnSparseF("fCutsMCFD","fCutsMCFD",5,dimAxesMC,minMC,maxMC);
  fCutsMC->GetAxis(0)->SetName("pt");
  fCutsMC->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fCutsMC->GetAxis(1)->SetName("NormDecLengthXY");
  fCutsMC->GetAxis(1)->SetTitle("Normalized XY decay length");
  fCutsMC->GetAxis(2)->SetName("CosPointXY");
  fCutsMC->GetAxis(2)->SetTitle("Cos#theta_{point}^{XY}");
  fCutsMC->GetAxis(3)->SetName("nProngs");
  fCutsMC->GetAxis(3)->SetTitle("Number of Decay Prongs");
  fCutsMC->GetAxis(4)->SetName("Mother_pt");
  fCutsMC->GetAxis(4)->SetTitle("Mother #it{p}_{T} (GeV/#it{c})");

  fPtMCGenStep=new TH1F("fPtMCGenStep","fPtMCGenStep",20,0.,20);
  return;
}

void AliHFsubtractBFDcuts::FillGenStep(AliAODMCParticle *dzeropart,Double_t pt/*=-1.*/,Double_t weight/*=1.*/){
  if(pt<0){
    pt=dzeropart->Pt();
  }
  else  fPtMCGenStep->Fill(pt,weight);
  return;
}

void AliHFsubtractBFDcuts::FillSparses(AliAODRecoDecayHF2Prong *dzerocand,Int_t isSelected,Double_t pt,Double_t massD0,Double_t massD0bar,Double_t weight,TClonesArray* mcArray){

  if(isSelected<=0||isSelected>3){
    Printf("isSelected = %d", isSelected);
    return;
  }
  if(massD0<0){
    dzerocand->InvMassD0(massD0,massD0bar);
    Printf("mass D0 =%f, massD0bar = %f ", massD0, massD0bar);
  }
  if(pt<0){
    pt=dzerocand->Pt();
  }
  Double_t normalDecayLengXY=dzerocand->NormalizedDecayLengthXY();
  Double_t cptangXY=dzerocand->CosPointingAngleXY();
  Double_t pointData[4]={massD0,pt,normalDecayLengXY,cptangXY};
  Double_t pointMC[5]={pt,normalDecayLengXY,cptangXY,0.,pt};

  if(isSelected==1||isSelected==3){
    fCutsData->Fill(pointData,weight);
  }
  if(isSelected==2||isSelected==3){
    pointData[0]=massD0bar;
    fCutsData->Fill(pointData,weight);
  }

  if(fIsMC){
    Int_t labCand=GetCandidateLabel(dzerocand, mcArray);
    UInt_t nProngs=0;
    Bool_t decayChain=kFALSE;
    if (!DetermineDecayType(labCand, mcArray, nProngs, decayChain, pointMC[4])) {
      AliDebug(3, "Error during the decay type determination!");
    }
    else {
      pointMC[3]=(Double_t)nProngs;
    }
    fCutsMC->Fill(pointMC,weight);
    //Int_t* a =0x0;
    //*a=1;
  }

  return;
}

Int_t AliHFsubtractBFDcuts::GetCandidateLabel(AliAODRecoDecayHF2Prong *dzerocand,TClonesArray* mcArray) const {
  if (mcArray) {
    Int_t labDau0=((AliAODTrack*)dzerocand->GetDaughter(0))->GetLabel();
    AliAODMCParticle* firstDau=(AliAODMCParticle*)mcArray->UncheckedAt(TMath::Abs(labDau0));
    Int_t labCand = firstDau->GetMother();
    return labCand;
  }
  else return -1;
}


Bool_t AliHFsubtractBFDcuts::DetermineDecayType(Int_t labCand, TClonesArray* mcArray, UInt_t& nProngs, Bool_t& decayChain, Double_t& motherPt) const {
  if (mcArray) {
    AliAODMCParticle* cand=(AliAODMCParticle*)mcArray->UncheckedAt(labCand);
    Int_t labMother = cand->GetMother();
    AliAODMCParticle* mother = (AliAODMCParticle*)mcArray->UncheckedAt(labMother);
    Int_t pdgMother = mother->GetPdgCode();
    if (pdgMother<0) pdgMother*=-1; // treat particles and anti-particles the same way
    if (pdgMother==4 || pdgMother==2212 || pdgMother==2112) { // prompt production
      nProngs=1;
      return kTRUE;
    }
    if ((pdgMother%1000)/100==4) { // chained decay of charmed hadrons
      return DetermineDecayType(labMother, mcArray, nProngs, decayChain, motherPt);
    }
    if ((pdgMother%1000)/100!=5) {
      AliDebug(3, "Found strange decay!");
      nProngs=0;
      return kFALSE;
    }
    nProngs = 0;
    for (Int_t i=0; i<mcArray->GetEntriesFast(); ++i) {
      if (((AliAODMCParticle*)mcArray->UncheckedAt(i))->GetMother()==labMother) ++nProngs;
    }
    motherPt=mother->Pt();
    return kTRUE;
  }
  else return kFALSE;
}
