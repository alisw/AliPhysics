/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
/* $Id: AliTOFT0maker.cxx,v 1.8 2010/01/19 16:32:20 noferini Exp $ */

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  This class contains the basic functions for the time zero              //
//  evaluation with TOF detector informations.                             //
// Use case in an analysis task:                                           //
//                                                                         //
// Create the object in the task constructor (fTOFmaker is a private var)  //
// AliESDpid *extPID=new AliESDpid();                                      //
// fTOFmaker = new AliTOFT0maker(extPID);                                  //
// fTOFmaker->SetTimeResolution(100.0); // if you want set the TOF res     //
// 115 ps is the TOF default resolution value                              //
//                                                                         //
// Use the RemakePID method in the task::Exec                              //
// Double_t* calcolot0;                                                    //
// calcolot0=fTOFmaker->RemakePID(fESD);                                   //
// //calcolot0[0] = calculated event time                                  // 
// //calcolot0[1] = event time time resolution                             //
// //calcolot0[2] = average event time for the current fill                //
// //calcolot0[3] = tracks at TOF                                          // 
// //calcolot0[4] = calculated event time (only TOF)                       //
// //calcolot0[5] = event time time resolution (only TOF)                  //
// //calcolot0[6] = sigma t0 fill                                          //
// //calcolot0[7] = tracks at TOF really used in tht algorithm             // 
//                                                                         //
// Let consider that:                                                      //
// - the PIF is automatically recalculated with the event time subtrction  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "AliTOFT0v1.h"
#include "AliTOFT0maker.h"
#include "AliPID.h"
#include "AliESDpid.h"
#include "AliESDEvent.h"
#include "TFile.h"
#include "TH1F.h"
#include "AliTOFcalib.h"
#include "AliTOFRunParams.h"
#include "TRandom.h"

ClassImp(AliTOFT0maker)
           
//____________________________________________________________________________ 
AliTOFT0maker::AliTOFT0maker():
  TObject(),
  fT0TOF(NULL),
  fPIDesd(NULL),
  fExternalPIDFlag(kFALSE),
  fTOFcalib(NULL),
  fNoTOFT0(0),
  fTimeResolution(100),
  fT0sigma(1000),
  fHmapChannel(0),
  fKmask(0),
  fT0width(150.)
{
  // ctr
  fCalculated[0] = 0;
  fCalculated[1] = 0;
  fCalculated[2] = 0;
  fCalculated[3] = 0;

  fT0TOF = new AliTOFT0v1();
  if(AliPID::ParticleMass(0) == 0) new AliPID();

  fPIDesd = new AliESDpid();

  fPtCutMin[0] = 0.3;
  fPtCutMin[1] = 0.5;
  fPtCutMin[2] = 0.6;
  fPtCutMin[3] = 0.7;
  fPtCutMin[4] = 0.8;
  fPtCutMin[5] = 0.9;
  fPtCutMin[6] = 1;
  fPtCutMin[7] = 1.2;
  fPtCutMin[8] = 1.5;
  fPtCutMin[9] = 2;

  fPtCutMax[0] = 0.5;
  fPtCutMax[1] = 0.6;
  fPtCutMax[2] = 0.7;
  fPtCutMax[3] = 0.8;
  fPtCutMax[4] = 0.9;
  fPtCutMax[5] = 1;
  fPtCutMax[6] = 1.2;
  fPtCutMax[7] = 1.5;
  fPtCutMax[8] = 2;
  fPtCutMax[9] = 3;

  /* init arrays */
  for (Int_t i = 0; i < 10; i++) {
    fT0pt[i] = 0.;
    fT0ptSigma[i] = 0.;
  }
}
//____________________________________________________________________________ 
AliTOFT0maker::AliTOFT0maker(AliESDpid *externalPID, AliTOFcalib *tofCalib):
    TObject(),
    fT0TOF(NULL),
    fPIDesd(externalPID),
    fExternalPIDFlag(kTRUE),
    fTOFcalib(tofCalib),
    fNoTOFT0(0),
    fTimeResolution(100),
    fT0sigma(1000),
    fHmapChannel(0),
    fKmask(0),
    fT0width(150.)
{
  // ctr
  fCalculated[0] = 0;
  fCalculated[1] = 0;
  fCalculated[2] = 0;
  fCalculated[3] = 0;

  fT0TOF = new AliTOFT0v1();
  if(AliPID::ParticleMass(0) == 0) new AliPID();

  if(!fPIDesd){
    fPIDesd = new AliESDpid();
    fExternalPIDFlag = kFALSE;
  }

  fPtCutMin[0] = 0.3;
  fPtCutMin[1] = 0.5;
  fPtCutMin[2] = 0.6;
  fPtCutMin[3] = 0.7;
  fPtCutMin[4] = 0.8;
  fPtCutMin[5] = 0.9;
  fPtCutMin[6] = 1;
  fPtCutMin[7] = 1.2;
  fPtCutMin[8] = 1.5;
  fPtCutMin[9] = 2;

  fPtCutMax[0] = 0.5;
  fPtCutMax[1] = 0.6;
  fPtCutMax[2] = 0.7;
  fPtCutMax[3] = 0.8;
  fPtCutMax[4] = 0.9;
  fPtCutMax[5] = 1;
  fPtCutMax[6] = 1.2;
  fPtCutMax[7] = 1.5;
  fPtCutMax[8] = 2;
  fPtCutMax[9] = 3;

  /* init arrays */
  for (Int_t i = 0; i < 10; i++) {
    fT0pt[i] = 0.;
    fT0ptSigma[i] = 0.;
  }
}

/* copy-constructor and operator= suppressed


//____________________________________________________________________________ 
AliTOFT0maker::AliTOFT0maker(const AliTOFT0maker & t) :
  TObject(t),
  fT0TOF(t.fT0TOF),
  fPIDESD(t.fPIDESD),
  fNoTOFT0(t.fNoTOFT0),
  fTimeResolution(t.fTimeResolution),
  fT0sigma(t.fT0sigma),
  fHmapChannel(t.fHmapChannel),
  fKmask(t.fKmask)
  {
  // copy ctr
}

//____________________________________________________________________________ 
AliTOFT0maker& AliTOFT0maker::operator=(const AliTOFT0maker &t)
{
  //
  // assign. operator
  //

  if (this == &t)
    return *this;

  TObject::operator=(t);
  fTimeResolution = t.fTimeResolution;
  fT0sigma = t.fT0sigma;

  return *this;
}

*/

//____________________________________________________________________________ 
AliTOFT0maker::~AliTOFT0maker()
{
  // dtor

  delete fT0TOF;
  if (!fExternalPIDFlag) delete fPIDesd;
}
//____________________________________________________________________________ 
Double_t* AliTOFT0maker::ComputeT0TOF(AliESDEvent *esd,Double_t t0time,Double_t t0sigma){
  //
  // Remake TOF PID probabilities
  //

  Double_t t0tof[4];

  if(fKmask) ApplyMask(esd);

  /* get T0 spread from TOFcalib if available otherwise use default value */
  if (fTOFcalib && esd) {
    AliTOFRunParams *runParams = fTOFcalib->GetRunParams();
    if (runParams && runParams->GetTimestamp(0) != 0) {
      Float_t t0spread = runParams->EvalT0Spread(esd->GetTimeStamp());
      SetT0FillWidth(t0spread);
    }
  }

  fT0TOF->Init(esd);
  AliTOFT0v1* t0maker= fT0TOF;
  t0maker->SetTimeResolution(fTimeResolution*1e-12);

  t0maker->DefineT0("all",1.5,3.0);
  t0tof[0] = t0maker->GetResult(0);
  t0tof[1] = t0maker->GetResult(1);
  t0tof[2] = t0maker->GetResult(2);
  t0tof[3] = t0maker->GetResult(3);

  Float_t lT0Current=0.;
  fT0sigma=1000;

//   Int_t nrun = esd->GetRunNumber();
  Double_t t0fill = 0.;

  t0time += t0fill;

  Float_t sigmaFill = fT0width;

  if(sigmaFill < 20) sigmaFill = 140;

  fCalculated[0]=-1000*t0tof[0]; // best t0
  fCalculated[1]=1000*t0tof[1]; // sigma best t0
  fCalculated[2] = t0fill;    //t0 fill
  fCalculated[3] = t0tof[2];  // n TOF tracks
  fCalculated[4]=-1000*t0tof[0]; // TOF t0
  fCalculated[5]=1000*t0tof[1]; // TOF t0 sigma
  fCalculated[6]=sigmaFill; // sigma t0 fill
  fCalculated[7] = t0tof[3];  // n TOF tracks used for T0

  if(fCalculated[1] < sigmaFill && TMath::Abs(fCalculated[0] - t0fill) < 500 && fCalculated[1] < fTimeResolution*1.2){
    fT0sigma=fCalculated[1];
    lT0Current=fCalculated[0];
  }
  else{
    fCalculated[4] = t0fill;
    fCalculated[5] = sigmaFill;
  }

  if(fCalculated[1] < 1 || fT0sigma > sigmaFill || fCalculated[1] > fTimeResolution* 1.2){
    fT0sigma =1000;
    fCalculated[4] = t0fill;
    fCalculated[5] = sigmaFill;
  }

  if(t0sigma < 1000){
    if(fT0sigma < 1000){
      Double_t w1 = 1./t0sigma/t0sigma;
      Double_t w2 = 1./fCalculated[1]/fCalculated[1];

      Double_t wtot = w1+w2;

      lT0Current = (w1*t0time + w2*fCalculated[0]) / wtot;
      fT0sigma = TMath::Sqrt(1./wtot);
    }
    else{
      lT0Current=t0time;
      fT0sigma=t0sigma;
    }
  }

  if(fT0sigma < sigmaFill && TMath::Abs(lT0Current - t0fill) < 500){
    fCalculated[1]=fT0sigma;
    fCalculated[0]=lT0Current;
  }

  if(fT0sigma >= 1000 || fNoTOFT0){
    lT0Current = t0fill;
    fT0sigma = sigmaFill;

    fCalculated[0] = t0fill;
    fCalculated[1] = sigmaFill;
  }

  // T0 pt bin
  for(Int_t i=0;i<10;i++){
   t0maker->DefineT0("all",fPtCutMin[i],fPtCutMax[i]);
    t0tof[0] = t0maker->GetResult(0);
    t0tof[1] = t0maker->GetResult(1);
    t0tof[2] = t0maker->GetResult(2);
    t0tof[3] = t0maker->GetResult(3);
    fT0pt[i] =-1000*t0tof[0]; // best t0
    fT0ptSigma[i] =1000*t0tof[1]; // sigma best t0

    if(fT0ptSigma[i] < sigmaFill  && fT0ptSigma[i] < fTimeResolution * 1.2 && TMath::Abs(fT0pt[i] - t0fill) < 500){
      // Ok T0
    }
    else{
      fT0pt[i] = t0fill;
      fT0ptSigma[i] = sigmaFill;
    }
  }
  //----
  SetTOFResponse();

  return fCalculated;
}
//____________________________________________________________________________ 
Double_t  *AliTOFT0maker::GetT0p(Float_t p){// [0]=to -- [1] = sigma T0
  Int_t i=0;
  while(p > fPtCutMin[i] && i < 10) i++;
  if(i > 0) i--;
  
  fT0cur[0] = fT0pt[i];
  fT0cur[1] = fT0ptSigma[i];
  return fT0cur;
}
//____________________________________________________________________________ 
void AliTOFT0maker::SetTOFResponse(){
    fPIDesd->GetTOFResponse().SetTimeResolution(TMath::Sqrt(fT0sigma*fT0sigma + fTimeResolution*fTimeResolution));
}
//____________________________________________________________________________ 
Float_t AliTOFT0maker::GetExpectedSigma(Float_t mom, Float_t tof, Float_t mass){
  Double_t *sigmaT0 = GetT0p(mom);
  fPIDesd->GetTOFResponse().SetTimeResolution(TMath::Sqrt(sigmaT0[1]*sigmaT0[1] + fTimeResolution*fTimeResolution));
  Float_t sigma = fPIDesd->GetTOFResponse().GetExpectedSigma(mom,tof,mass);
  fPIDesd->GetTOFResponse().SetTimeResolution(TMath::Sqrt(fT0sigma*fT0sigma + fTimeResolution*fTimeResolution));

  return sigma;
}
//____________________________________________________________________________ 
void AliTOFT0maker::ApplyT0TOF(AliESDEvent *esd){
  //
  // Recalculate TOF PID probabilities
  //

  // subtruct t0 for each track
  Int_t ntracks = esd->GetNumberOfTracks();
  
  while (ntracks--) {
    AliESDtrack *t=esd->GetTrack(ntracks);
    
    if ((t->GetStatus()&AliESDtrack::kTOFout)==0) continue;
    
    Double_t time=t->GetTOFsignal();
    Float_t p = t->GetP();

    Double_t *t0=GetT0p(p);
    time -= t0[0];
    t->SetTOFsignal(time);
  }
  //
}
//____________________________________________________________________________ 
void  AliTOFT0maker::LoadChannelMap(char *filename){
  // Load the histo with the channel off map
  TFile *f= new TFile(filename);
  if(!f){
    printf("Cannot open the channel map file (%s)\n",filename);
    return;
  }
  
  fHmapChannel = (TH1F *) f->Get("hChEnabled");
  
  if(!fHmapChannel){
    printf("Cannot laod the channel map histo (from %s)\n",filename);
    return;
  }
    
}
//____________________________________________________________________________ 
void AliTOFT0maker::ApplyMask(AliESDEvent * const esd){
  // Switch off the disable channel
  if(!fHmapChannel){
    printf("Channel Map is not available\n");
    return;
  }
  
  Int_t ntracks = esd->GetNumberOfTracks();
  
  while (ntracks--) {
    AliESDtrack *t=esd->GetTrack(ntracks);    

    if ((t->GetStatus()&AliESDtrack::kTOFout)==0) continue;

    Int_t chan = t->GetTOFCalChannel();
 
    if(fHmapChannel->GetBinContent(chan) < 0.01){
      t->ResetStatus(AliESDtrack::kTOFout);
    }
  }
}

Float_t  
AliTOFT0maker::TuneForMC(AliESDEvent *esd){ // return true T0 event
  //
  // tune for MC data
  //

  Float_t TOFtimeResolutionDefault=80;

  Float_t t0 = gRandom->Gaus(0.,fT0width); 

  Float_t extraSmearing = 0;

  if(fTimeResolution > TOFtimeResolutionDefault){
    extraSmearing = TMath::Sqrt(fTimeResolution*fTimeResolution - TOFtimeResolutionDefault*TOFtimeResolutionDefault);
  }

  // subtruct t0 for each track
  Int_t ntracks = esd->GetNumberOfTracks();
  
  while (ntracks--) {
    AliESDtrack *t=esd->GetTrack(ntracks);
    
    if ((t->GetStatus()&AliESDtrack::kTOFout)==0) continue;
    
    /* check if channel is enabled */
    if (fTOFcalib && !fTOFcalib->IsChannelEnabled(t->GetTOFCalChannel())) {
      /* reset TOF status */
      t->ResetStatus(AliESDtrack::kTOFin);
      t->ResetStatus(AliESDtrack::kTOFout);
      t->ResetStatus(AliESDtrack::kTOFrefit);
      t->ResetStatus(AliESDtrack::kTOFpid);
    }

    Double_t time=t->GetTOFsignal();

    time += t0;

    if(extraSmearing>0){
      Float_t smearing = gRandom->Gaus(0.,extraSmearing);
      time += smearing;
    }

    t->SetTOFsignal(time);
  }
  //
  return t0;
}
