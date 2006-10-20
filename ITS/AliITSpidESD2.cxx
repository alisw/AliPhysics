/**************************************************************************
 * Copyright(c) 2005-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------//
// ITS PID class --- method # 2                                          //
//                                                                       //
//                                                                       //
//The PID is based on the likelihood of all the four ITS' layers,        //
//without using the truncated mean for the dE/dx. The response           //
//functions for each layer are convoluted Landau-Gaussian functions.     // 
// Origin: Elena Bruna bruna@to.infn.it, Massimo Masera masera@to.infn.it//
//-----------------------------------------------------------------------//

#include "AliITSpidESD2.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliITStrackV2.h"
#include "AliITSclusterV2.h"
#include "AliITStrackerMI.h"
#include "AliITSLoader.h"
#include "AliITSPident.h"
#include "AliITSSteerPid.h"
#include "AliLog.h"

ClassImp(AliITSpidESD2)
//_________________________________________________________________________
  AliITSpidESD2::AliITSpidESD2():AliITSpidESD(),
fTracker(0),
fLoader(0),
fSp(0)
{ //
  //  The main constructor
}
//_________________________________________________________________________
AliITSpidESD2::AliITSpidESD2(AliITStrackerMI* tracker,AliITSLoader* loader):AliITSpidESD(),
fTracker(0),
fLoader(0),
fSp(0)
{ //
  //  The main constructor
  fTracker=tracker;
  fLoader=loader;
  fSp=new AliITSSteerPid();
  fSp->InitLayer();
}
//_________________________________________________________________________
AliITSpidESD2::~AliITSpidESD2(){
  //destructor

  if(fSp)delete fSp;

}
//______________________________________________________________________
AliITSpidESD2::AliITSpidESD2(const AliITSpidESD2 &ob) :AliITSpidESD(ob),
fTracker(ob.fTracker),
fLoader(ob.fLoader),
fSp(ob.fSp) 
{
  // Copy constructor
 
}

//______________________________________________________________________
AliITSpidESD2& AliITSpidESD2::operator=(const AliITSpidESD2& ob ){
  // Assignment operator
  this->~AliITSpidESD2();
  new(this) AliITSpidESD2(ob);
  return *this;
}

//_________________________________________________________________________
Int_t AliITSpidESD2::MakePID(AliESD *event)
{

  //
  //  This function calculates the "detector response" PID probabilities 
  //
  Double_t xr,par[5];
  AliITStrackV2* track=0;
  fLoader->LoadRecPoints();
  TTree *cTree=fLoader->TreeR();
  fTracker->LoadClusters(cTree);

  Int_t ntrk=event->GetNumberOfTracks();
  Double_t momits;
  // for (Int_t i=0; i<ntrk; i++) {
  for (Int_t i=0; i<ntrk; i++) {
    AliESDtrack *esdtr=event->GetTrack(i);
    if ((esdtr->GetStatus()&AliESDtrack::kITSin )==0)
      if ((esdtr->GetStatus()&AliESDtrack::kITSout)==0) continue;

    track = new AliITStrackV2(*esdtr);
    Double_t dEdxsignal=track->GetdEdx();
    track->GetExternalParameters(xr,par);
    if (par[4]!=0) {
      Float_t lamb=TMath::ATan(par[3]);
      momits=1/(TMath::Abs(par[4])*TMath::Cos(lamb));
    }
    else {
      AliWarning("Null particle momentum in ITS");
      momits = 0.;
    } 
    Double_t snp=track->GetSnp();
    Double_t tgl=track->GetTgl();
    const Int_t kns=AliPID::kSPECIES;
    Double_t condprobfun[kns];
    for(Int_t ii=0;ii<kns;ii++)condprobfun[ii]=0;
    Int_t cluindsdd1 = track->GetClusterIndex(3);
    Int_t cluindsdd2 = track->GetClusterIndex(2);
    Int_t cluindssd1 = track->GetClusterIndex(1);
    Int_t cluindssd2 = track->GetClusterIndex(0);
    Float_t q1,q1corr,q2,q2corr,q3,q3corr,q4,q4corr;
    AliITSclusterV2* clu1=(AliITSclusterV2*)fTracker->GetCluster(cluindsdd1);
    if(clu1!=0){
      q1=clu1->GetQ(); 
      q1corr=q1*TMath::Sqrt((1-snp*snp)/(1+tgl*tgl));
    }
    else{
      q1=-99;
      q1corr=-99;
    }
	
    AliITSclusterV2* clu2=(AliITSclusterV2*)fTracker->GetCluster(cluindsdd2);
    if(clu2!=0){
      q2=clu2->GetQ();
      q2corr=q2*TMath::Sqrt((1-snp*snp)/(1+tgl*tgl));
    }
    else{
      q2=-99;
      q2corr=-99;
    }
    
    AliITSclusterV2* clu3=(AliITSclusterV2*)fTracker->GetCluster(cluindssd1);
    if(clu3!=0){
      q3=clu3->GetQ();
      q3corr=q3*TMath::Sqrt((1-snp*snp)/(1+tgl*tgl));
    }
    else{
      q3=-99;
      q3corr=-99;
    }
    
    AliITSclusterV2* clu4=(AliITSclusterV2*)fTracker->GetCluster(cluindssd2);
    if(clu4!=0){
      q4=clu4->GetQ();
      q4corr=q4*TMath::Sqrt((1-snp*snp)/(1+tgl*tgl));
    }
    else{
      q4=-99;
      q4corr=-99;
    }
    Float_t qlay[4]={q1corr,q2corr,q3corr,q4corr};
    Float_t prip=0.33;
    Float_t prik=0.33;
    Float_t pripi=0.33;
    Float_t prie=0.;
    Double_t invPt=track->Get1Pt();
    AliITSPident mypid(momits,invPt,dEdxsignal,fSp,qlay,prip,prik,pripi,prie); 
    condprobfun[0]=0.;//el
    condprobfun[1]=0.;//mu
    condprobfun[2]=mypid.GetProdCondFunPi();//pi
    condprobfun[3]=mypid.GetProdCondFunK();//kaon
    condprobfun[4]=mypid.GetProdCondFunPro();//pro

    esdtr->SetITSpid(condprobfun);

    delete track;
   }
  fTracker->UnloadClusters();
  fLoader->UnloadRecPoints();
  return 0;
}
