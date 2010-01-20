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
// fTOFmaker = new AliTOFT0maker();                                        //
// fTOFmaker->SetTimeResolution(115.0e-12); // if you want set the TOF res //
// 115 ps is the TOF default resolution value                              //
//                                                                         //
// Use the RemakePID method in the task::Exec                              //
// Double_t* calcolot0;                                                    //
// calcolot0=fTOFmaker->RemakePID(fESD);                                   //
// //calcolot0[0] = calculated event time                                  // 
// //calcolot0[1] = event time time resolution                             //
// //calcolot0[2] = average event time for the current fill                //
//                                                                         //
// Let consider that:                                                      //
// - the PIF is automatically recalculated with the event time subtrction  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <stdlib.h>

#include "AliTOFT0v1.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalibHisto.h"
#include "AliPID.h"
#include "AliESDpid.h"

ClassImp(AliTOFT0maker)
           
//____________________________________________________________________________ 
AliTOFT0maker::AliTOFT0maker() :
TObject(),
  fCalib(new AliTOFcalibHisto()),
  fESDswitch(0),
  fTimeResolution(115),
  fT0sigma(1000)
{
  //
  // ctr
  //
  
  fCalculated[0] = 0;
  fCalculated[1] = 0;
  fCalculated[2] = 0;

  fCalib->LoadCalibPar();

  if(AliPID::ParticleMass(0) == 0) new AliPID();
}
//____________________________________________________________________________ 
AliTOFT0maker::AliTOFT0maker(const AliTOFT0maker & t) :
TObject(),
  fCalib(t.fCalib),
  fESDswitch(t.fESDswitch),
  fTimeResolution(t.fTimeResolution),
  fT0sigma(t.fT0sigma)
{
}

//____________________________________________________________________________ 
AliTOFT0maker& AliTOFT0maker::operator=(const AliTOFT0maker &t)
{
  //
  // assign. operator
  //

  if (this == &t)
    return *this;
  fCalib = t.fCalib;
  fESDswitch = t.fESDswitch;
  fTimeResolution = t.fTimeResolution;
  fT0sigma = t.fT0sigma;

  return *this;
}
//____________________________________________________________________________ 
AliTOFT0maker::~AliTOFT0maker()
{
  // dtor
  if(fCalib) delete fCalib;
}
//____________________________________________________________________________ 
Double_t* AliTOFT0maker::RemakePID(AliESDEvent *esd,Double_t t0time,Double_t t0sigma){
  //
  // Remake TOF PID probabilities
  //

  Double_t *t0tof;

  AliTOFT0v1* t0maker=new AliTOFT0v1(esd);
  t0maker->SetCalib(fCalib);
  t0maker->SetTimeResolution(fTimeResolution*1e-12);

  if(! fESDswitch){
    t0tof=t0maker->DefineT0RawCorrection("all");
    TakeTimeRawCorrection(esd);
  }
  else t0tof=t0maker->DefineT0("all");

  Float_t lT0Current=0.;
  fT0sigma=1000;

  Int_t nrun = esd->GetRunNumber();
  Double_t t0fill = GetT0Fill(nrun);

  fCalculated[0]=-1000*t0tof[0];
  fCalculated[1]=1000*t0tof[1];
  fCalculated[2] = t0fill;

  if(fCalculated[1] < 150 && TMath::Abs(fCalculated[0] - t0fill) < 500){
    fT0sigma=fCalculated[1];
    lT0Current=fCalculated[0];
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

  if(fT0sigma >= 1000){
    lT0Current = t0fill;
    fT0sigma = 135;

    fCalculated[0] = t0fill;
    fCalculated[1] = 150;
  }

  RemakeTOFpid(esd,lT0Current);

  return fCalculated;
}
//____________________________________________________________________________ 
void AliTOFT0maker::TakeTimeRawCorrection(AliESDEvent * const esd){
  //
  // Take raw corrections for time measurements
  //

  Int_t ntracks = esd->GetNumberOfTracks();

  while (ntracks--) {
    AliESDtrack *t=esd->GetTrack(ntracks);
    
    if ((t->GetStatus()&AliESDtrack::kTOFout)==0) continue;
    
    Double_t time=t->GetTOFsignalRaw();
    Double_t tot = t->GetTOFsignalToT();
    Int_t chan = t->GetTOFCalChannel();
    Double_t corr = fCalib->GetFullCorrection(chan,tot) - fCalib->GetCorrection(AliTOFcalibHisto::kTimeSlewingCorr,chan,0);
    time -= corr*1000;

    Int_t crate = Int_t(fCalib->GetCalibMap(AliTOFcalibHisto::kDDL,chan));

    if(crate == 63 || crate == 62){
      time += 9200;
    }

    t->SetTOFsignal(time);
  }
}
//____________________________________________________________________________ 
void AliTOFT0maker::RemakeTOFpid(AliESDEvent *esd,Float_t timezero){
  //
  // Recalculate TOF PID probabilities
  //

  AliESDpid pidESD;
  pidESD.GetTOFResponse().SetTimeResolution(TMath::Sqrt(fT0sigma*fT0sigma + fTimeResolution*fTimeResolution));
  pidESD.MakePID(esd,kFALSE,timezero);
  
}
//____________________________________________________________________________ 
Double_t AliTOFT0maker::GetT0Fill(Int_t nrun) const {
  //
  // Return T0 of filling
  //

  Double_t t0;
  if(nrun==104065) t0= 1771614;
  else if(nrun==104068) t0= 1771603;
  else if(nrun==104070) t0= 1771594;
  else if(nrun==104073) t0= 1771610;
  else if(nrun==104080) t0= 1771305;
  else if(nrun==104083) t0= 1771613;
  else if(nrun==104157) t0= 1771665;
  else if(nrun==104159) t0= 1771679;
  else if(nrun==104160) t0= 1771633;
  else if(nrun==104316) t0= 1764344;
  else if(nrun==104320) t0= 1764342;
  else if(nrun==104321) t0= 1764371;
  else if(nrun==104439) t0= 1771750;
  else if(nrun==104792) t0= 1771755;
  else if(nrun==104793) t0= 1771762;
  else if(nrun==104799) t0= 1771828;
  else if(nrun==104800) t0= 1771788;
  else if(nrun==104801) t0= 1771796;
  else if(nrun==104802) t0= 1771775;
  else if(nrun==104803) t0= 1771795;
  else if(nrun==104824) t0= 1771751;
  else if(nrun==104825) t0= 1771763;
  else if(nrun==104845) t0= 1771792;
  else if(nrun==104852) t0= 1771817;
  else if(nrun==104864) t0= 1771825;
  else if(nrun==104865) t0= 1771827;
  else if(nrun==104867) t0= 1771841;
  else if(nrun==104876) t0= 1771856;
  else if(nrun==104878) t0= 1771847;
  else if(nrun==104879) t0= 1771830;
  else if(nrun==104892) t0= 1771837;
  else t0= 487;

  if(fESDswitch) t0 -= 487;
  
  return t0;
}
