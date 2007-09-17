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

/* $Id$ */

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
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliITStrackV2.h"
#include "AliITSRecPoint.h"
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
Int_t AliITSpidESD2::MakePID(AliESDEvent *event)
{

  //
  //  This function calculates the "detector response" PID probabilities 
  //
  Double_t xr,par[5];
  AliITStrackV2* track=0;
  fLoader->LoadRecPoints();
  TTree *cTree=fLoader->TreeR();
  fTracker->LoadClusters(cTree);
  printf("==== Landau Fit PID ITS ====== \n");
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
    Int_t cluind[12];
    for(Int_t ii=0;ii<12;ii++){
      cluind[ii]=track->GetClusterIndex(ii);
    }
    AliITSRecPoint* cluarr[12];
    Float_t qclu[8],qclucorr[8],nlay[8];
    for(Int_t ii=0;ii<8;ii++){
      qclu[ii]=0;
      qclucorr[ii]=0;
      nlay[ii]=0;
    }
    Int_t jj=0;
    cout<<"track = "<<i<<endl;
    for(Int_t ij=0;ij<12;ij++){
      cluind[ij]=track->GetClusterIndex(ij);
      cout<<cluind[ij]<<endl;
      if(cluind[ij]>0){
	cluarr[ij]=(AliITSRecPoint*)fTracker->GetCluster(cluind[ij]);
	Int_t lay=cluarr[ij]->GetLayer();
	cout<<"lay = "<<lay<<endl;
	if(lay>1){//sdd+ssd only
	  qclu[jj]=cluarr[ij]->GetQ(); 
	  qclucorr[jj]=qclu[jj]*TMath::Sqrt((1-snp*snp)/(1+tgl*tgl));
	  nlay[jj]=lay;
	  jj++;
	}
	else qclucorr[jj]=-1;
      }
    }

    Float_t prip=0.33;
    Float_t prik=0.33;
    Float_t pripi=0.33;
    Float_t prie=0.;
    AliITSPident mypid(momits,dEdxsignal,fSp,qclucorr,nlay,prip,prik,pripi,prie); 
    condprobfun[0]=mypid.GetProdCondFunPi();//el --PID in the ITS does not distinguish among Pi,el,mu
    condprobfun[1]=mypid.GetProdCondFunPi();//mu
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
