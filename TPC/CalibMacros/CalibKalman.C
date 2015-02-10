/// \file CalibKalman.C
///
/// ~~~
/// .x ~/rootlogon.C
/// gSystem->Load("libANALYSIS");
/// gSystem->Load("libTPCcalib");
/// 
/// .L $ALICE_ROOT/TPC/CalibMacros/CalibKalman.C+
/// ~~~



#include "TTreeStream.h"
#include "TFile.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "AliExternalTrackParam.h"
#include "AliTPCcalibTime.h"
#include "AliTPCkalmanTime.h"


void FitKalman(Int_t maxCount=50){
  TFile fcalib("CalibObjects.root");
  TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
  AliTPCcalibTime * calibTime = ( AliTPCcalibTime *)array->FindObject("calibTime");
  AliTPCkalmanTime kalman;
  //
  //

  TTreeSRedirector *pcstream = new TTreeSRedirector("debugFit.root");
  TH2D * hdz = calibTime->GetHistVdrift()->Projection(1,0);
  TH2D * hpress = calibTime->GetHistVdrift()->Projection(4,0);
  TH2D * htempA = calibTime->GetHistVdrift()->Projection(5,0);
  TH2D * htempC = calibTime->GetHistVdrift()->Projection(6,0);
  
  TH1 * hentries =  calibTime->GetHistVdrift()->Projection(0);
  Int_t nentries =hentries->GetNbinsX();
  TH1 *htemp=0;
  Int_t counter=0;
  for (Int_t ibin=0; ibin<nentries;ibin++){
    Double_t icon  = hentries->GetBinContent(ibin);
    if (icon<5) continue;
    if (counter>maxCount) break;
    Double_t time = hentries->GetBinCenter(ibin);
    //
    htemp   = hdz->ProjectionY("aaa",ibin,ibin+1);
    Double_t dzmean   = htemp->GetMean();
    Double_t dzrms    =  htemp->GetRMS();
    delete htemp;
    //
    htemp   = hpress->ProjectionY("aaa",ibin,ibin+1);
    Double_t dpress    = htemp->GetMean();
    Double_t rpress    =  htemp->GetRMS();
    delete htemp;
    htemp   = htempA->ProjectionY("aaa",ibin,ibin+1);
    Double_t mtemp    =  htemp->GetMean();
    Double_t rtemp    =  htemp->GetRMS();
    delete htemp;
    
    printf("ibin=%d\ttime=%f\tentries=%f\t%f\t%f\t%f\t%f\t%f\t%f\n",ibin,time,icon,dzmean,dzrms,dpress,rpress, mtemp,rtemp);
    if (counter==0) {
      kalman.Init(time,0,1,0.05,0.02);
    }
    Double_t sigma0t = 0.005/(24*60*60.);
    kalman.Propagate(time,sigma0t,pcstream);
    Float_t dvdriftrel = dzmean/500.;
    Float_t dverror    = (dzrms+0.1)/500.;
    Double_t  temperature = mtemp+273.15;
    Double_t  povertMeas = dpress/temperature;
    if (mtemp<15) continue;
    
    Float_t  nominalTemp = 19.03+273.15;
    Float_t  nominalPress= 973;
    Double_t povertNom = nominalPress/nominalTemp  ;
    Float_t dptratio   = (povertMeas-povertNom)/povertNom;
    kalman.Update(dvdriftrel,dverror,dptratio,pcstream);
    kalman.fState->Print();
    kalman.fCovariance->Print();
    (*pcstream)<<"fit"<<
      "time="<<time<<
      "dvdrift="<<dvdriftrel<<
      "dverror="<<dverror<<
      "dzrms="<<dzrms<<
      "dpt="<<dptratio<<
      "k.="<<&kalman<<
      "\n";
    counter++;
  }
  delete pcstream;
}

