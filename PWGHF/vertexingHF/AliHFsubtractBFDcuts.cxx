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
#include "AliAODRecoDecayHF2Prong.h"
#include "AliHFsubtractBFDcuts.h"
AliHFsubtractBFDcuts::AliHFsubtractBFDcuts() :
  TNamed(),
  fisMC(kTRUE),
  fPtMCGenStep(0x0),
  fcutsData(0x0),
  fcutsMC(0x0)
{
  // default constructor
}

AliHFsubtractBFDcuts::AliHFsubtractBFDcuts(const char* name, const char* title) :
  TNamed(name,title),
  fisMC(kTRUE),
  fPtMCGenStep(0x0),
  fcutsData(0x0),
  fcutsMC(0x0)
{
  // default constructor with name
}

void AliHFsubtractBFDcuts::InitHistos(){
  Int_t dimAxes[4]={ 500, 24, 30, 400};// mass, pt, normLXY, cosPointXY
  Double_t min[4]={1.700, 0., 0.,0.96};
  Double_t max[4]={2.200,24.,30.,  1.};
  fcutsData=new THnSparseF("fcutsDataFD","fcutsDataFD",4,dimAxes,min,max);
  fcutsData->GetAxis(0)->SetName("mass");
  fcutsData->GetAxis(0)->SetTitle("Mass (K,#pi) (GeV/#it{c^{2}})");
  fcutsData->GetAxis(1)->SetName("pt");
  fcutsData->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fcutsData->GetAxis(2)->SetName("NormDecLengthXY");
  fcutsData->GetAxis(2)->SetTitle("Normalized XY decay length");
  fcutsData->GetAxis(3)->SetName("CosPointXY");
  fcutsData->GetAxis(3)->SetTitle("Cos#theta_{point}^{XY}");

  Int_t dimAxesMC[3]={24, 30, 400};// pt, normLXY, cosPointXY
  Double_t minMC[3]={ 0., 0.,0.96};
  Double_t maxMC[3]={24.,30.,  1.};
  fcutsMC=new THnSparseF("fcutsMCFD","fcutsMCFD",3,dimAxesMC,minMC,maxMC);
  fcutsMC->GetAxis(0)->SetName("pt");
  fcutsMC->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fcutsMC->GetAxis(1)->SetName("NormDecLengthXY");
  fcutsMC->GetAxis(1)->SetTitle("Normalized XY decay length");
  fcutsMC->GetAxis(2)->SetName("CosPointXY");
  fcutsMC->GetAxis(2)->SetTitle("Cos#theta_{point}^{XY}");

  fPtMCGenStep=new TH1F("fPtMCGenStep","fPtMCGenStep",20,0.,20);
  return;
}

void AliHFsubtractBFDcuts::FillGenStep(AliAODMCParticle *dzeropart,Double_t pt){
  if(pt<0){
    pt=dzeropart->Pt();
  }
  else  fPtMCGenStep->Fill(pt);
  return;
}

void AliHFsubtractBFDcuts::FillSparses(AliAODRecoDecayHF2Prong *dzerocand,Int_t isSelected,Double_t pt,Double_t massD0,Double_t massD0bar){

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
  Double_t pointdata[4]={massD0,pt,normalDecayLengXY,cptangXY};
  Double_t pointMC[3]={pt,normalDecayLengXY,cptangXY};

  if(isSelected==1||isSelected==3){
    fcutsData->Fill(pointdata);
  }
  if(isSelected==2||isSelected==3){
    pointdata[0]=massD0bar;
    fcutsData->Fill(pointdata);
  }

  Printf("isMC = %d ", fisMC);
  if(fisMC){
    Printf("filling!!!!!!!!!!!!!!");
    fcutsMC->Fill(pointMC);
  }

  return;
}
