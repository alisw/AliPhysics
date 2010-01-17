#include <Riostream.h>
#include <stdlib.h>

#include "AliTOFT0v1.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalibHisto.h"
#include "AliPID.h"
#include "AliESDpid.h"

ClassImp(AliTOFT0maker)
           
//____________________________________________________________________________ 
  AliTOFT0maker::AliTOFT0maker(): fESDswitch(0), fTimeResolution(115),  fT0sigma(1000)
{
  fCalib = new AliTOFcalibHisto();
  fCalib->LoadCalibPar();

  if(AliPID::ParticleMass(0) == 0) new AliPID();
}
//____________________________________________________________________________ 
AliTOFT0maker::~AliTOFT0maker()
{
  // dtor
  if(fCalib) delete fCalib;
}
//____________________________________________________________________________ 
Double_t* AliTOFT0maker::RemakePID(AliESDEvent *esd,Double_t t0time,Double_t t0sigma){
  Double_t* calcolot0;

  AliTOFT0v1* t0maker=new AliTOFT0v1(esd);
  t0maker->SetCalib(fCalib);
  t0maker->SetTimeResolution(fTimeResolution*1e-12);

  if(! fESDswitch){
    calcolot0=t0maker->DefineT0RawCorrection("all");
    TakeTimeRawCorrection(esd);
  }
  else calcolot0=t0maker->DefineT0("all");

  calcolot0[0]*=-1000;
  calcolot0[1]*=1000;

  Float_t T0Current;
  fT0sigma=1000;

  if(calcolot0[1] < 300){
    fT0sigma=calcolot0[1];
    T0Current=calcolot0[0];
  }

  if(t0sigma < 1000){
    if(fT0sigma < 1000){
      Double_t w1 = 1./t0sigma/t0sigma;
      Double_t w2 = 1./calcolot0[1]/calcolot0[1];

      Double_t wtot = w1+w2;

      T0Current = (w1*t0time + w2*calcolot0[0]) / wtot;
      fT0sigma = TMath::Sqrt(1./wtot);
    }
    else{
      T0Current=t0time;
      fT0sigma=t0sigma;
    }
  }

  Int_t nrun = esd->GetRunNumber();
  Double_t t0fill = GetT0Fill(nrun);

  if(fT0sigma >= 1000){
    T0Current = t0fill;
    fT0sigma = 135;
  }

  RemakeTOFpid(esd,T0Current);

  return calcolot0;
}
//____________________________________________________________________________ 
void AliTOFT0maker::TakeTimeRawCorrection(AliESDEvent *esd){
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
  AliESDpid pidESD;
  pidESD.GetTOFResponse().SetTimeResolution(TMath::Sqrt(fT0sigma*fT0sigma + fTimeResolution*fTimeResolution));
  pidESD.MakePID(esd,kFALSE,timezero);
  
}
//____________________________________________________________________________ 
Double_t AliTOFT0maker::GetT0Fill(Int_t nrun){
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
  else t0= 1771837;

  if(fESDswitch) t0 -= 487;
  
  return t0;
}
