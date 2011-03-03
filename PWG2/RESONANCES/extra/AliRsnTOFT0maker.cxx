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
/* $Id: AliRsnTOFT0maker.cxx,v 1.8 2010/01/19 16:32:20 noferini Exp $ */

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  This class contains the basic functions for the time zero              //
//  evaluation with TOF detector informations.                             //
// Use case in an analysis task:                                           //
//                                                                         //
// Create the object in the task constructor (fTOFmakerANA is a private var)  //
// fTOFmakerANA = new AliRsnTOFT0maker();                                        //
// fTOFmakerANA->SetTimeResolution(115.0e-12); // if you want set the TOF res //
// 115 ps is the TOF default resolution value                              //
//                                                                         //
// Use the RemakePID method in the task::Exec                              //
// Double_t* calcolot0;                                                    //
// calcolot0=fTOFmakerANA->RemakePID(fESD);                                   //
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
#include "AliTOFcalibHisto.h"
#include "AliPID.h"
#include "AliESDpid.h"
#include "TFile.h"
#include "AliRsnTOFT0maker.h"

ClassImp(AliRsnTOFT0maker)

//____________________________________________________________________________
AliRsnTOFT0maker::AliRsnTOFT0maker():
   fSettings(kNone),
   fCalib(new AliTOFcalibHisto()),
   fESDswitch(0),
   fTimeResolution(115),
   fT0sigma(1000),
   fHmapChannel(0),
   fKmask(0),
   fNoTOFT0(0)
{
   fCalculated[0] = 0;
   fCalculated[1] = 0;
   fCalculated[2] = 0;
   fCalculated[3] = 0;

   fCalib->LoadCalibPar();

   if (AliPID::ParticleMass(0) == 0) new AliPID();
   
   Int_t i;
   for (i = 0; i < 4; i++) fCalculated[i] = 0.0;
}
//____________________________________________________________________________
AliRsnTOFT0maker::AliRsnTOFT0maker(const AliRsnTOFT0maker & t) :
   TObject(),
   fSettings(t.fSettings),
   fCalib(t.fCalib),
   fESDswitch(t.fESDswitch),
   fTimeResolution(t.fTimeResolution),
   fT0sigma(t.fT0sigma),
   fHmapChannel(t.fHmapChannel),
   fKmask(t.fKmask),
   fNoTOFT0(t.fNoTOFT0)
{
}

//____________________________________________________________________________
AliRsnTOFT0maker& AliRsnTOFT0maker::operator=(const AliRsnTOFT0maker &t)
{
//
   // assign. operator
   //

   fSettings = t.fSettings;
   if (this == &t)
      return *this;
   fCalib = t.fCalib;
   fESDswitch = t.fESDswitch;
   fTimeResolution = t.fTimeResolution;
   fT0sigma = t.fT0sigma;

   return *this;
}
//____________________________________________________________________________
AliRsnTOFT0maker::~AliRsnTOFT0maker()
{
   // dtor
   if (fCalib) delete fCalib;
}
//____________________________________________________________________________
Double_t* AliRsnTOFT0maker::RemakePID(AliESDEvent *esd, Double_t t0time, Double_t t0sigma)
{
   //
   // Remake TOF PID probabilities
   //

   Double_t *t0tof;

   if (fKmask) ApplyMask(esd);

   AliTOFT0v1* t0makerANA = new AliTOFT0v1(esd);
//   t0makerANA->SetCalib(fCalib);
   t0makerANA->SetTimeResolution(fTimeResolution * 1e-12);

   if (! fESDswitch) {
      TakeTimeRawCorrection(esd);
   }

   t0tof = t0makerANA->DefineT0("all");

   Float_t lT0Current = 0.;
   fT0sigma = 1000;

   Int_t nrun = 0;//esd->GetRunNumber();
   Double_t t0fill = 175;//GetT0Fill(nrun);
   if (fSettings == kPass2 || fSettings == kPass4) {    // cambiato!!!!!!!!
      nrun = esd->GetRunNumber();
      t0fill = GetT0Fill(nrun);
   } else if (fSettings == kLHC09d10) { // e' il MC del pass2?
      t0fill = GetT0Fill(nrun);
   }

   fCalculated[0] = -1000 * t0tof[0];
   fCalculated[1] = 1000 * t0tof[1];
   fCalculated[2] = t0fill;
   fCalculated[3] = t0tof[2];

   if (fCalculated[1] < 150 && TMath::Abs(fCalculated[0] - t0fill) < 500) {
      fT0sigma = fCalculated[1];
      lT0Current = fCalculated[0];
   }

   if (t0sigma < 1000) {
      if (fT0sigma < 1000) {
         Double_t w1 = 1. / t0sigma / t0sigma;
         Double_t w2 = 1. / fCalculated[1] / fCalculated[1];

         Double_t wtot = w1 + w2;

         lT0Current = (w1 * t0time + w2 * fCalculated[0]) / wtot;
         fT0sigma = TMath::Sqrt(1. / wtot);
      } else {
         lT0Current = t0time;
         fT0sigma = t0sigma;
      }
   }

   if (fT0sigma >= 1000 || fNoTOFT0) {
      lT0Current = t0fill;
      fT0sigma = 135;

      fCalculated[0] = t0fill;
      fCalculated[1] = 150;
   }

   RemakeTOFpid(esd, lT0Current);

   return fCalculated;
}
//____________________________________________________________________________
void AliRsnTOFT0maker::TakeTimeRawCorrection(AliESDEvent * const esd)
{
   //
   // Take raw corrections for time measurements
   //

   Int_t ntracks = esd->GetNumberOfTracks();

   while (ntracks--) {
      AliESDtrack *t = esd->GetTrack(ntracks);

      if ((t->GetStatus()&AliESDtrack::kTOFout) == 0) continue;

      Double_t time = t->GetTOFsignalRaw();
      Double_t tot = t->GetTOFsignalToT();
      Int_t chan = t->GetTOFCalChannel();
      Double_t corr = fCalib->GetFullCorrection(chan, tot) - fCalib->GetCorrection(AliTOFcalibHisto::kTimeSlewingCorr, chan, 0);
      time -= corr * 1000;

      //Int_t crate = Int_t(fCalib->GetCalibMap(AliTOFcalibHisto::kDDL,chan));

//     if(crate == 63 || crate == 62){
//       time += 9200;

//    }

//     if(crate == 63 || crate == 62|| crate == 61){
//  printf("%i) %f\n",crate,time);
//     getchar();
//   }
      t->SetTOFsignal(time);
   }
}
//____________________________________________________________________________
void AliRsnTOFT0maker::RemakeTOFpid(AliESDEvent *esd, Float_t timezero)
{
   //
   // Recalculate TOF PID probabilities
   //
   AliESDpid pidESD;
   pidESD.GetTOFResponse().SetTimeResolution(TMath::Sqrt(fT0sigma * fT0sigma + fTimeResolution * fTimeResolution));
   pidESD.MakePID(esd, kFALSE, timezero);

}
//____________________________________________________________________________
Double_t AliRsnTOFT0maker::GetT0Fill(Int_t nrun) const
{
   //
   // Return T0 of filling
   //

   Double_t t0;
   if (nrun == 104065) t0 = 1771614;
   else if (nrun == 104068) t0 = 1771603;
   else if (nrun == 104070) t0 = 1771594;
   else if (nrun == 104073) t0 = 1771610;
   else if (nrun == 104080) t0 = 1771305;
   else if (nrun == 104083) t0 = 1771613;
   else if (nrun == 104157) t0 = 1771665;
   else if (nrun == 104159) t0 = 1771679;
   else if (nrun == 104160) t0 = 1771633;
   else if (nrun == 104316) t0 = 1764344;
   else if (nrun == 104320) t0 = 1764342;
   else if (nrun == 104321) t0 = 1764371;
   else if (nrun == 104439) t0 = 1771750;
   else if (nrun == 104792) t0 = 1771755;
   else if (nrun == 104793) t0 = 1771762;
   else if (nrun == 104799) t0 = 1771828;
   else if (nrun == 104800) t0 = 1771788;
   else if (nrun == 104801) t0 = 1771796;
   else if (nrun == 104802) t0 = 1771775;
   else if (nrun == 104803) t0 = 1771795;
   else if (nrun == 104824) t0 = 1771751;
   else if (nrun == 104825) t0 = 1771763;
   else if (nrun == 104845) t0 = 1771792;
   else if (nrun == 104852) t0 = 1771817;
   else if (nrun == 104864) t0 = 1771825;
   else if (nrun == 104865) t0 = 1771827;
   else if (nrun == 104867) t0 = 1771841;
   else if (nrun == 104876) t0 = 1771856;
   else if (nrun == 104878) t0 = 1771847;
   else if (nrun == 104879) t0 = 1771830;
   else if (nrun == 104892) t0 = 1771837;
   else t0 = 487;

   if (fESDswitch) t0 -= 487;
   else { if (fSettings == kPass4) t0 -= 37 * 1024 * 24.4; }

   return t0;
}

//____________________________________________________________________________
void  AliRsnTOFT0maker::LoadChannelMap(char *filename)
{
   // Load the histo with the channel off map
   TFile *f = new TFile(filename);
   if (!f) {
      printf("Cannot open the channel map file (%s)\n", filename);
      return;
   }

   fHmapChannel = (TH1F *) f->Get("hChEnabled");

   if (!fHmapChannel) {
      printf("Cannot laod the channel map histo (from %s)\n", filename);
      return;
   }

}
//____________________________________________________________________________
void AliRsnTOFT0maker::ApplyMask(AliESDEvent *esd)
{
   // Switch off the disable channel
   if (!fHmapChannel) {
      printf("Channel Map is not available\n");
      return;
   }

   Int_t ntracks = esd->GetNumberOfTracks();

   while (ntracks--) {
      AliESDtrack *t = esd->GetTrack(ntracks);

      if ((t->GetStatus()&AliESDtrack::kTOFout) == 0) continue;

      Int_t chan = t->GetTOFCalChannel();

      if (fHmapChannel->GetBinContent(chan) < 0.01) {
         t->ResetStatus(AliESDtrack::kTOFout);
      }
   }
}
