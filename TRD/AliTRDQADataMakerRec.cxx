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

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Produces the data needed to calculate the quality assurance.          //
//  All data must be mergeable objects.                                   //
//                                                                        //
//  Author:                                                               //
//    Sylwester Radomski (radomski@physi.uni-heidelberg.de)               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1D.h> 
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TF1.h>
#include <TCanvas.h>

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliTRDcluster.h"
#include "AliTRDQADataMakerRec.h"
#include "AliTRDgeometry.h"
//#include "AliTRDdataArrayI.h"
#include "AliTRDrawStream.h"

#include "AliQAChecker.h"

ClassImp(AliTRDQADataMakerRec)

//____________________________________________________________________________ 
  AliTRDQADataMakerRec::AliTRDQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kTRD), "TRD Quality Assurance Data Maker")
{
  //
  // Default constructor
}

//____________________________________________________________________________ 
AliTRDQADataMakerRec::AliTRDQADataMakerRec(const AliTRDQADataMakerRec& qadm) :
  AliQADataMakerRec()
{
  //
  // Copy constructor 
  //

  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 

}

//__________________________________________________________________
AliTRDQADataMakerRec& AliTRDQADataMakerRec::operator=(const AliTRDQADataMakerRec& qadm)
{
  //
  // Equal operator.
  //

  this->~AliTRDQADataMakerRec();
  new(this) AliTRDQADataMakerRec(qadm);
  return *this;

}

//____________________________________________________________________________ 
void AliTRDQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //
  // Detector specific actions at end of cycle
  //
  //TStopwatch watch;
  //watch.Start();
  /**/
  AliDebug(AliQAv1::GetQADebugLevel(), "End of TRD cycle");
  
  if (task == AliQAv1::kRECPOINTS) {
    
    TH1D * hist = new TH1D("fitHist", "", 200, -0.5, 199.5);

    // loop over event types
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {

      //list[specie]->Print();
      
      // fill detector map;
      for(Int_t i = 0 ; i < 540 ; i++) {
        Double_t v = ((TH1D*)list[specie]->At(0))->GetBinContent(i+1);
        Int_t sm = i/30;
        Int_t det = i%30;
        TH2D *detMap = (TH2D*)list[specie]->At(87);
        Int_t bin = detMap->FindBin(sm, det);
        detMap->SetBinContent(bin, v);
      }
      
      // Rec points full chambers
      for (Int_t i = 0 ; i < 540 ; i++) {
        //AliDebug(AliQAv1::GetQADebugLevel(), Form("I = %d", i));
        //TH1D *h = ((TH2D*)list[specie]->At(1))->ProjectionY(Form("qaTRD_recPoints_amp_%d",i), i+1, i+1);
        hist->Reset();
	
	// project TH2D into TH1D 
        for(Int_t b = 1 ; b < hist->GetXaxis()->GetNbins()-1 ; b++) {
          Double_t xvalue = hist->GetBinCenter(b);
          Int_t bin = ((TH2D*)list[specie]->At(1))->FindBin(i,xvalue);
          Double_t value =  ((TH2D*)list[specie]->At(1))->GetBinContent(bin);
          hist->SetBinContent(b, value);
        }
        //AliDebug(AliQAv1::GetQADebugLevel(), Form("Sum = %d %f\n", i, hist->GetSum()));
        if (hist->GetSum() < 100) 
	  continue; // not enougth data in a chamber
	
        hist->Fit("landau", "q0", "goff", 10, 180);
        TF1 *fit = hist->GetFunction("landau");
        ((TH1D*)list[specie]->At(12))->Fill(fit->GetParameter(1));
        ((TH1D*)list[specie]->At(13))->Fill(fit->GetParameter(2));
      }


      // time-bin by time-bin sm by sm
      for(Int_t i=0; i<18; i++) { // loop over super-modules      
	for(Int_t j=0; j<kTimeBin; j++) { // loop over time bins
	  
          hist->Reset();
          for(Int_t b = 1 ; b < hist->GetXaxis()->GetNbins()-1 ; b++) {
            Double_t xvalue = hist->GetBinCenter(b);
            Double_t svalue = 0.0;
            for(Int_t det = i*30 ; det < (i+1)*30 ; det++) { // loop over detectors
              Int_t bin = ((TH3D*)list[specie]->At(10))->FindBin(det,j,xvalue);
              Double_t value =  ((TH3D*)list[specie]->At(10))->GetBinContent(bin);
              svalue += value;
            }
            //AliDebug(AliQAv1::GetQADebugLevel(), Form("v = %f\n", value));
            hist->SetBinContent(b, svalue);
          }
	  
          if (hist->GetSum() < 100) 
            continue;
	  
	  hist->Fit("landau", "q0", "goff", 10, 180);
	  TF1 *fit = hist->GetFunction("landau");
	  
	  TH1D *h1 = (TH1D*)list[specie]->At(14+18+i);
	  h1->SetMarkerStyle(20);
	  Int_t bin = h1->FindBin(j);
	  // printf("%d %d %d\n", det, j, bin);
	  
	  Double_t value = TMath::Abs(fit->GetParameter(1));
	  Double_t error = TMath::Abs(fit->GetParError(1));
	 
	  if (value/error < 3) continue; // insuficient statistics
	  
	  h1->SetBinContent(bin, value);
	  h1->SetBinError(bin, error);	
	}
      }
	
      // for numerical convergence
      TF1 *form = new TF1("formLandau", "landau", 0, 200);
	
      // time-bin by time-bin chamber by chamber
      for (Int_t i=0; i<540; i++) {
	for(Int_t j=0; j<kTimeBin; j++) {
	    
	  hist->Reset();
	  for(Int_t b = 1 ; b < hist->GetXaxis()->GetNbins()-1 ; b++) {
	    Double_t xvalue = hist->GetBinCenter(b);
	    Int_t bin = ((TH3D*)list[specie]->At(10))->FindBin(i,j,xvalue);
	    Double_t value =  ((TH3D*)list[specie]->At(10))->GetBinContent(bin);
	    //AliDebug(AliQAv1::GetQADebugLevel(), Form("v = %f\n", value));
	    hist->SetBinContent(b, value);
	  }
	  
	  if (hist->GetSum() < 100) 
	    continue;
	  
	  form->SetParameters(1000, 60, 20);
	  hist->Fit(form, "q0", "goff", 20, 180);
	  
	  Int_t sm = i/30;
	  Int_t det = i%30;
	  TH2D *h2 = (TH2D*)list[specie]->At(14+sm);
	  Int_t bin = h2->FindBin(det,j);
	  // printf("%d %d %d\n", det, j, bin);
	  
	  Double_t value = TMath::Abs(form->GetParameter(1));
	  Double_t error = TMath::Abs(form->GetParError(1));
	  
	  if (value/error < 3) continue;
	  
	  h2->SetBinContent(bin, value);
	  h2->SetBinError(bin, error);
	}
      }
      if (hist) 
	delete hist;
    }
  }
  //////////////////////////
  // const Int_t knbits = 6;
  // const char *suf[knbits] = {"TPCi", "TPCo", "TPCz", "TRDo", "TRDr", "TRDz"};
  //const char *sufRatio[4] = {"TRDrTRDo", "TRDoTPCo", "TRDrTPCo", "TRDzTPCo"};

  if (task == AliQAv1::kESDS) {
    
    const Int_t knRatio = 4;
    const Int_t kN[knRatio] = {4,3,4,5};
    const Int_t kD[knRatio] = {3,1,1,3}; 
    
    // create ratios
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      for(Int_t type = 0 ; type < 2 ; type++) {
        for(Int_t i = 0 ; i < knRatio ; i++) {
          TH1D *ratio = (TH1D*)list[specie]->At(19 + 2*i + type);
          TH1D *histN = (TH1D*)list[specie]->At(3 + 2*kN[i] + type);
          TH1D *histD = (TH1D*)list[specie]->At(3 + 2*kD[i] + type);
          BuildRatio(ratio, histN, histD);
          //ratio->Reset();
          //ratio->Add(histN);
          //ratio->Divide(histD);
        }
      }
      // ratio for the fraction of electrons per stack
      TH1D *histN = (TH1D*)list[specie]->At(33);
      TH1D *histD = (TH1D*)list[specie]->At(32);
      TH1D *ratio = (TH1D*)list[specie]->At(34);
      BuildRatio(ratio, histN, histD);
    }
  }
  // call the checker
  AliQAChecker::Instance()->Run(AliQAv1::kTRD, task, list) ;    

}

//____________________________________________________________________________ 
void AliTRDQADataMakerRec::InitESDs()
{
  //
  // Create ESDs histograms in ESDs subdir
  //
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  const Int_t kNhist = 36+5+4;

  TH1 *hist[kNhist];
  Int_t histoCounter = -1 ;

  hist[++histoCounter] = new TH1D("qaTRD_esd_ntracks", ";Number of tracks", 300, -0.5, 299.5);
  hist[++histoCounter] = new TH1D("qaTRD_esd_sector", ";Sector", 18, -0.5, 17.7);
  hist[++histoCounter] = new TH1D("qaTRD_esd_bits", ";Bits", 64, -0.5, 63.5);

  const Int_t knbits = 6;
  const char *suf[knbits] = {"TPCi", "TPCo", "TPCz", "TRDo", "TRDr", "TRDz"};

  // histo = 3
  for(Int_t i=0; i<knbits; i++) { 
    hist[++histoCounter] = new TH1D(Form("qaTRD_esd_pt%s",suf[i]), ";p_{T} (GeV/c);", 100, 0, 10);
    hist[++histoCounter] = new TH1D(Form("qaTRD_esd_trdz%s", suf[i]), ";z (cm)", 200, -400, 400); 
  }

  hist[++histoCounter] = new TH1D("qaTRD_esd_clsTRDo", "TRDo;number of clusters", 180, -0.5, 179.5);;
  hist[++histoCounter] = new TH1D("qaTRD_esd_clsTRDr", "TRDr;number of clusters", 180, -0.5, 179.5);;
  hist[++histoCounter] = new TH1D("qaTRD_esd_clsTRDz", "TRDz;number of clusters", 180, -0.5, 179.5);;
  //hist[++histoCounter] = new TH1D("qaTRD_esd_clsRatio", ";cluster ratio", 100, 0., 1.3);;

  hist[++histoCounter] = new TH2D("qaTRD_esd_sigMom", ";momentum (GeV/c);signal", 100, 0, 5, 200, 0, 1e3);

  // end of cycle plots (hist 19)
  const char *sufRatio[4] = {"TRDrTRDo", "TRDoTPCo", "TRDrTPCo", "TRDzTPCo"};

  for(int i=0; i<4; i++) {
    hist[++histoCounter] = new TH1D(Form("qaTRD_esd_pt%s",sufRatio[i]), 
				    Form("Efficiency in Pt %s;p_{T};eff", sufRatio[i]),
				    100, 0, 10);

    hist[++histoCounter] = new TH1D(Form("qaTRD_esd_trdz%s",sufRatio[i]), 
				    Form("Efficiency in Z %s;z (cm);eff", sufRatio[i]),
				    200, -400, 400);
  }

  // 27 - 31
  hist[27] = new TH1D("qaTRD_esd_quality", ";quality", 120, 0, 12);
  hist[28] = new TH1D("qaTRD_esd_budget", ";NN", 110, -1000, 100);
  hist[29] = new TH1D("qaTRD_esd_chi2", ";chi2", 200, 0, 100);
  hist[30] = new TH1D("qaTRD_esd_timeBin", ";time bin", 7, -0.5, 6.5);
  hist[31] = new TH1D("qaTRD_esd_pidQuality", "pid Quality", 7, -0.5, 6.5);

  // stack by stack electron identyfication
  hist[32] = new TH1D("qaTRD_esd_tracksStack", "number of all tracks;stack", 90, -0.5, 89.5);
  hist[33] = new TH1D("qaTRD_esd_electronStack", "number of electron tracks;stack", 90, -0.5, 89.5);
  hist[34] = new TH1D("qaTRD_esd_elRatioStack", "fraction of electron tracks;stack", 90, -0.5, 89.5);
  hist[35] = new TH1D("qaTRD_esd_thetaOut", ";tan(theta);", 100, -1, 1);
  
  const char *partType[5] = {"Electron", "Muon", "Pion", "Kaon", "Proton"}; 

  for(Int_t i=0; i<AliPID::kSPECIES; i++)
    hist[36+i] = new TH1D(Form("qaTRD_esd_pid%d",i),
			  Form("%s;probability",partType[i]), 100, 0, 1);
 
  // dE/dX vs momentum in three regions
  const char *zoneName[4] = {"total charge", "ampilification range", "plateau", "TR range"};
 
  // prepare the scale from 0.1 to 10 GeV
  const Int_t nscalex= 50;
  Double_t scalex[nscalex+1];
  Double_t dd = (TMath::Log(10) - TMath::Log(0.5)) / nscalex;
  for(Int_t ix=0; ix<nscalex+1; ix++) {
    scalex[ix] = 0.5 * TMath::Exp(dd * ix);
  }

  const Int_t nscaley = 50;
  Double_t scaley[nscaley+1];
  for(Int_t iy=0; iy<nscaley+1; iy++) {
    scaley[iy] = iy * (3e3/nscaley);
  }
    

  for(Int_t i=0; i<4; i++) {
    hist[41+i] = new TH2D(Form("qaTRD_esd_signalPzone_%d",i), 
			  Form("%s;momentum (GeV/c);singal (a.u.)", zoneName[i]),
			  nscalex, scalex, nscaley, scaley);
  }

  for(Int_t i=0; i<kNhist; i++) {
    //hist[i]->Sumw2();
    Add2ESDsList(hist[i], i, !expert, image);
  }

}

//____________________________________________________________________________ 
void AliTRDQADataMakerRec::InitRecPoints()
{
  //
  // Create Reconstructed Points histograms in RecPoints subdir
  //
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  //printf("Helo from Init rec points\n");

  const Int_t kNhist = 14 + 4 * 18 + 2 + 9;// + 540;
  TH1 *hist[kNhist];

  hist[0] = new TH1D("qaTRD_recPoints_det", ";Detector ID of the cluster", 540, -0.5, 539.5);
  hist[1] = new TH2D("qaTRD_recPoints_amp", ";Amplitude", 540, -0.5, 539, 200, -0.5, 199.5);
  hist[2] = new TH1D("qaTRD_recPoints_npad", ";Number of Pads", 12, -0.5, 11.5);

  hist[3] = new TH1D("qaTRD_recPoints_dist2", ";residuals [2pad]", 100, -1, 1);
  hist[4] = new TH1D("qaTRD_recPoints_dist3", ";residuals [3pad]", 100, -1, 1);
  hist[5] = new TH1D("qaTRD_recPoints_dist4", ";residuals [4pad]", 100, -1, 1);
  hist[6] = new TH1D("qaTRD_recPoints_dist5", ";residuals [5pad]", 100, -1, 1);

  hist[7] = new TH2D("qaTRD_recPoints_rowCol", ";row;col", 16, -0.5, 15.5, 145, -0.5, 144.5);
  hist[8] = new TH1D("qaTRD_recPoints_time", ";time bin", kTimeBin, -0.5, kTimeBin-0.5);
  hist[9] = new TH1D("qaTRD_recPoints_nCls", ";number of clusters", 500, -0.5, 499.5);

  hist[10] = new TH3D("qaTRD_recPoints_sigTime", ";chamber;time bin;signal", 
		      540, -0.5, 539.5, kTimeBin, -0.5, kTimeBin-0.5, 200, -0.5, 199.5);
  hist[11] = new TProfile("qaTRD_recPoints_prf", ";distance;center of gravity"
                         , 120, -0.6, 0.6, -1.2, 1.2, "");

  hist[12] = new TH1D("qaTRD_recPoints_ampMPV", ";amplitude MPV", 150, 0, 150);
  hist[13] = new TH1D("qaTRD_recPoints_ampSigma", ";amplitude Sigma", 200, 0, 200); 
  
  // chamber by chamber
  for(Int_t i=0; i<18; i++) {
    hist[14+i] = new TH2D(Form("qaTRD_recPoints_sigTime_sm%d",i), Form("sm%d;det;time bin"), 
			30, -0.5, 29.5, kTimeBin, -0.5, kTimeBin-0.5);
    hist[14+i]->SetMinimum(0);
    hist[14+i]->SetMaximum(150);
  }
 
  // time bin by time bin sm-by-sm
  for(Int_t i=0; i<18; i++) {
    hist[14+18+i] = new TH1D(Form("qaTRD_recPoints_sigTimeShape_sm%d", i), 
			     Form("sm%d;time bin;signal"),
			     kTimeBin, -0.5, kTimeBin-0.5);

    hist[14+18+i]->SetMaximum(150);    
  }

  // str = 50
  for(Int_t i=0; i<18; i++) {
    hist[50+i] = new TH1D(Form("qaTRD_recPoints_nCls_sm%d",i),
			  Form("sm%d;time bin;number of clusters",i),
			  kTimeBin, -0.5, kTimeBin-0.5);
  }

  // str = 68
  for(Int_t i=0; i<18; i++) {
    hist[68+i] = new TH1D(Form("qaTRD_recPoints_totalCharge_sm%d", i),
			  Form("sm%d;time bin;total charge", i),
			  kTimeBin, -0.5, kTimeBin-0.5);
  }

  hist[86] = new TH1D("qaTRD_recPoints_signal", ";amplitude", 400, -0.5, 399.5);
  hist[87] = new TH2D("qaTRD_recPoints_detMap", ";sm;chamber", 18, -0.5, 17.5, 30, -0.5, 29.5);

  
  // amplitude as a function of the pad size
  for(Int_t i=0; i<9; i++) {
    hist[88+i] = new TH1D(Form("qaTRD_recPoints_signalNpad_%d", i+2), ";amplitude, ADC", 400, -0.5, 399.5); 
  }
  
  // one 2D histogram per chamber
  //  for(Int_t i=0; i<540; i++) {
  //  hist[88+i] = new TH2D(Form("qaTRD_recPoints_map%d", i), ";col;row", 16, -0.5, 15.5, 144, -0.5, 143.5);
  //}


  for(Int_t i=0; i<kNhist; i++) {
    //hist[i]->Sumw2();
    Add2RecPointsList(hist[i], i, !expert, image);
  }
}

//____________________________________________________________________________ 
void AliTRDQADataMakerRec::InitRaws()
{
  //
  // create Raws histograms in Raws subdir
  //
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  AliInfo("Initialization of QA for Raw Data");

  //const Int_t kSM = 18;
  //const Int_t kNCh = 540;
  const Int_t kNhist = 6;
  TH1 *hist[kNhist];

  // link monitor
  const char *linkName[3] = {"smLink", "smBeaf", "smData"};
  for(Int_t i=0; i<3; i++) {
    hist[i] = new TH2D(Form("qaTRD_raws_%s", linkName[i]), 
		       ";super module;link", 
		       18, -0.5, 17.5, 60, -0.5, 59.5);
  }

  hist[3] = new TH1D("qaTRD_raws_errorHC", ";error ID;", 18, -3.5, 14.5);
  hist[4] = new TH1D("qaTRD_raws_errorMCM", ";error ID;", 18, -3.5, 14.5);
  hist[5] = new TH1D("qaTRD_raws_errorADC", ";error ID;", 18, -3.5, 14.5);

  /*
  // four histograms to be published
  hist[0] = new TH1D("qaTRD_raws_det", ";detector", 540, -0.5, 539.5);
  hist[1] = new TH1D("qaTRD_raws_sig", ";signal", 100, -0.5, 99.5);
  hist[2] = new TH1D("qaTRD_raws_timeBin", ";time bin", 40, -0.5, 39.5); 
  hist[3] = new TH1D("qaTRD_raws_smId", ";supermodule", 18, -0.5, 17.5);
  //
  
  // one double per MCM (not published)
  const Int_t kNMCM = 30 * 8 * 16;
  for(Int_t i=0; i<kSM; i++)
    hist[4+i] = new TH1D(Form("qaTRD_raws_sm%d",i),"",kNMCM, -0.5, kNMCM-0.5); 
  */

  // register
  for(Int_t i=0; i<kNhist; i++) {
    //hist[i]->Sumw2();
    Add2RawsList(hist[i], i, !expert, image, !saveCorr);
  }

}

//____________________________________________________________________________
void AliTRDQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //
  // Make QA data from ESDs
  //

  Int_t nTracks = esd->GetNumberOfTracks();
  GetESDsData(0)->Fill(nTracks);

  // track loop
  for (Int_t iTrack = 0; iTrack<nTracks; iTrack++) {

    AliESDtrack *track = esd->GetTrack(iTrack);
    const AliExternalTrackParam *paramOut = track->GetOuterParam();
    const AliExternalTrackParam *paramIn = track->GetInnerParam();

    // long track ..
    if (!paramIn) continue;
    if (!paramOut) continue;

    // not a kink
    if (track->GetKinkIndex(0) > 0) continue; 

    Double_t extZ = GetExtZ(paramIn);
    if (TMath::Abs(extZ) > 320) continue; // acceptance cut

    // .. in the acceptance
    Int_t sector = GetSector(paramOut->GetAlpha());
    Int_t stack = GetStack(paramOut);

    UInt_t u = 1;
    UInt_t status = track->GetStatus();
    for(Int_t bit=0; bit<32; bit++) 
      if (u<<bit & status) GetESDsData(2)->Fill(bit);

    const Int_t knbits = 6; 
    Int_t bit[6] = {0,0,0,0,0,0};    
    bit[0] = status & AliESDtrack::kTPCin;
    bit[1] = status & AliESDtrack::kTPCout;
    bit[2] = (status & AliESDtrack::kTPCout) && !(status & AliESDtrack::kTRDout);
    bit[3] = status & AliESDtrack::kTRDout;
    bit[4] = status & AliESDtrack::kTRDrefit;
    bit[5] = (status & AliESDtrack::kTRDout) && !(status & AliESDtrack::kTRDrefit);

    // transverse momentum
    //const Double_t *val = paramOut->GetParameter(); // parameters at the Outer plane
    Double_t pt = paramOut->Pt(); //1./TMath::Abs(val[4]);

    for(Int_t b=0; b<knbits; b++) {
      if (bit[b]) {
	GetESDsData(2*b+3)->Fill(pt); 
	GetESDsData(2*b+4)->Fill(extZ);
      }
    }

    // clusters
    for(Int_t b=0; b<3; b++) 
      if (bit[3+b]) GetESDsData(b+15)->Fill(track->GetTRDncls0());

    // refitted only
    if (!bit[4]) continue;

    //fQuality->Fill(track->GetTRDQuality());
    //fBudget->Fill(track->GetTRDBudget());
    //fSignal->Fill(track->GetTRDsignal());

    GetESDsData(1)->Fill(sector);
    GetESDsData(18)->Fill(track->GetP(), track->GetTRDsignal());

    GetESDsData(27)->Fill(track->GetTRDQuality());
    GetESDsData(28)->Fill(track->GetTRDBudget());
    GetESDsData(29)->Fill(track->GetTRDchi2());
    GetESDsData(30)->Fill(track->GetTRDTimBin(0));
    GetESDsData(31)->Fill(track->GetTRDntrackletsPID());
    
    
    // dedx
    for(Int_t k=0; k<4; ++k) {
      Double_t dedx = 0;
      for(Int_t j=0; j<6; j++) {
	dedx += track->GetTRDslice(j, k-1);
      }
      GetESDsData(41+k)->Fill(paramOut->GetP(), dedx/6.);
    }

    // probabilities
    if (status & AliESDtrack::kTRDpid) {
      for(Int_t k=0; k<AliPID::kSPECIES; ++k) 
	GetESDsData(36+k)->Fill(track->GetTRDpid(k));
    }

    // probabilities uniformity
    if (track->GetTRDntrackletsPID() < 6) continue;
    GetESDsData(35)->Fill(paramOut->GetZ()/paramOut->GetX());
    
    Int_t idx = 5 * sector + stack;
    GetESDsData(32)->Fill(idx); // all tracks
    if (track->GetTRDpid(AliPID::kElectron) > 0.9) 
      GetESDsData(33)->Fill(idx); // electrons only

    

    /*
    hist[27] = new TH1D("qaTRD_esd_quality", ";quality", 120, 0, 12);
    hist[28] = new TH1D("qaTRD_esd_budget", ";NN", 110, -1000, 100);
    hist[29] = new TH1D("qaTRD_esd_chi2", ";chi2", 300, 0, 100);
    hist[30] = new TH1D("qaTRD_esd_timeBin", 7, -0.5, 6.5);
    hist[31] = new TH1D("qaTRD_esd_pidQuality", 7, -0.5, 6.5);
    */

    /*
    // PID only
    if (status & AliESDtrack::kTRDpid) {

      for(Int_t l=0; l<6; l++) fTime->Fill(track->GetTRDTimBin(l));

      // fill pid histograms
      Double_t trdr0 = 0; //, tpcr0 = 0;
      Int_t trdBestPid = 5; //, tpcBestPid = 5;  // charged
      const Double_t kminPidValue = 0.9;

      //Double_t pp[5];
      //track->GetTPCpid(pp); // ESD inconsequence

      for(Int_t pid=0; pid<5; pid++) {
	
	trdr0 += track->GetTRDpid(pid);
	//tpcr0 += pp[pid];
	
	fTrdPID[pid]->Fill(track->GetTRDpid(pid));
	//fTpcPID[pid]->Fill(pp[pid]);
	
	if (track->GetTRDpid(pid) > kminPidValue) trdBestPid = pid;
	//if (pp[pid] > kminPidValue) tpcBestPid = pid;
      }

      fTrdPID[5]->Fill(trdr0); // check unitarity
      fTrdSigMomPID[trdBestPid]->Fill(track->GetP(), track->GetTRDsignal());

      //fTpcPID[5]->Fill(tpcr0); // check unitarity
      //fTpcSigMomPID[tpcBestPid]->Fill(track->GetP(), track->GetTPCsignal());
    }
    */

  }

}

//______________________________________________________________________________
Int_t AliTRDQADataMakerRec::GetSector(Double_t alpha) const 
{
  //
  // Gets the sector number 
  //

  Double_t size = TMath::DegToRad() * 20.; // shall use TRDgeo
  if (alpha < 0) alpha += 2*TMath::Pi();
  Int_t sector = (Int_t)(alpha/size);
  return sector;

}
//______________________________________________________________________________

Int_t AliTRDQADataMakerRec::GetStack(const AliExternalTrackParam *paramOut) const
{
  //
  // calculates the stack the track is in
  //
  
  const Double_t L = -0.9;
  const Double_t W = (2*L)/5;

  Double_t tan = paramOut->GetZ() / paramOut->GetX();
  Double_t pos = (tan - L) / W;
  return (Int_t) pos;
}

//______________________________________________________________________________
Double_t AliTRDQADataMakerRec::GetExtZ(const AliExternalTrackParam *in) const 
{
  //
  // Returns the Z position at the entry to TRD
  // using parameters from the TPC in
  //

  const Double_t kX0 = 300;

  Double_t x = in->GetX();
  const Double_t *par = in->GetParameter();
  Double_t theta = par[3];
  Double_t z = in->GetZ();

  Double_t zz = z + (kX0-x) * TMath::Tan(theta);
  return zz;

}

//____________________________________________________________________________
void AliTRDQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //
  // Makes QA data from raw data
  //
  
  AliInfo("Execution of QA for Raw data");
  //printf("Execution of QA for Raw data\n");

  // create raw reader TB
  rawReader->Reset();
  rawReader->SelectEquipment(0, 1024, 1041);
  rawReader->Select("TRD");
  
  // AliTRDrawStreamTB::AllowCorruptedData();
  // AliTRDrawStreamTB::DisableStackNumberChecker();
  // AliTRDrawStreamTB::DisableStackLinkNumberChecker();
  // AliTRDrawStreamTB::DisableSkipData();

  //AliTRDrawStreamTB *data = (AliTRDrawStreamTB*)AliTRDrawStreamBase::GetRawStream(rawReader);
  AliTRDrawStream *data = (AliTRDrawStream *)AliTRDrawStreamBase::GetRawStream(rawReader);  

  //if (raw->IsA()->GetName())


  // import the histograms
  TH2D *fSMLink[3];
  for(Int_t i=0; i<3; i++) {
    fSMLink[i] = (TH2D*)GetRawsData(i);
    //printf("address = %d\n", fSMLink[i]);
  }
    
  TH1D *fErrorHC =  (TH1D*)GetRawsData(3);
  TH1D *fErrorMCM = (TH1D*)GetRawsData(4);
  TH1D *fErrorADC = (TH1D*)GetRawsData(5);


  // loop over super-modules
  while (data->NextBuffer()) {
    
    //printf("processing next buffer\n");
    // check sm
    Int_t sm = rawReader->GetEquipmentId() - 1024;
    if (sm < 0 || sm > 18) return;
    
    // loop over links
    for (Int_t istack = 0; istack < 5; istack++) {	
      for (Int_t ilink = 0; ilink < 12; ilink++) {
	
	//Int_t det = sm * 30 + istack * 6 + ilink/2;	
      
	// check if data delivered
	if (!(data->IsLinkActiveInStack(istack, ilink))) continue;
	fSMLink[0]->Fill(sm, istack * 12 + ilink);
	
	// check if beaf-beaf
	if (data->GetLinkMonitorError(istack, ilink)) {
	  fSMLink[1]->Fill(sm, istack * 12 + ilink);
	  continue;
	}
	
	// fill histogram with HC header errors
	Int_t nErrHc = 0;
	
	nErrHc = FillBits(fErrorHC, data->GetH0ErrorCode(istack, ilink), 0);
	if (!nErrHc) fErrorHC->Fill(-3);
	      
	nErrHc = FillBits(fErrorHC, data->GetH1ErrorCode(istack, ilink), 2);
	if (!nErrHc) fErrorHC->Fill(-2);
	      
	nErrHc = FillBits(fErrorHC, data->GetHCErrorCode(istack, ilink), 4);
	if (!nErrHc) fErrorHC->Fill(-1);
	
	// data integrity protection 
	if (data->GetH0ErrorCode(istack, ilink) > 0) continue;
	if (data->GetH1ErrorCode(istack, ilink) > 0) continue;
	
	fSMLink[2]->Fill(sm, istack * 12 + ilink);	

	// loop over MCMs
	for (Int_t imcm = 0; imcm < data->GetHCMCMmax(istack, ilink); imcm++ ){
	  
	  Int_t nErrMcm = 0;
	
	  nErrMcm = FillBits(fErrorMCM, data->GetMCMhdErrorCode(istack, ilink, imcm), 0);
	  if (!nErrMcm) fErrorMCM->Fill(-3);
	
	  nErrMcm = FillBits(fErrorMCM, data->GetMCMADCMaskErrorCode(istack, ilink, imcm), 5);
	  if (!nErrMcm) fErrorMCM->Fill(-2);
	
	  nErrMcm = FillBits(fErrorMCM, data->GetMCMErrorCode(istack, ilink, imcm), 10);
	  if (!nErrMcm) fErrorMCM->Fill(-1);
	
	  // MCM protection
	  if ( (data->GetMCMhdErrorCode(istack,ilink,imcm)) & 2 ) continue;

	  // loop over ADC chanels 	
	  for (Int_t iadc=0; iadc < data->GetADCcount(istack, ilink, imcm); iadc++) {
	  
	    // fill ADC error bits
	    Int_t nErrAdc = FillBits(fErrorADC, data->GetADCErrorCode(), 0);
	    if (!nErrAdc) fErrorADC->Fill(-1);
	  }
	} // mcm 
      }
    } // link

    //printf("buffer analyzed\n");
  }
  
  // clean up

  AliInfo("sucessfull execution of QA for TRD raw data");
  //printf("sucessfull execution of QA for TRD raw data\n");
  

  /*
  // 157
  // T9 -- T10

  //const Int_t kSM = 18;
  //const Int_t kROC = 30;
  const Int_t kROB = 8;
  //const Int_t kLayer = 6;
  //const Int_t kStack = 5;
  const Int_t kMCM = 16;
  //  const Int_t kADC = 22;

  rawReader->Reset() ; 
  //AliTRDrawStreamBase::SetRawStreamVersion("TB");
  AliTRDrawStreamBase *raw = AliTRDrawStreamBase::GetRawStream(rawReader);
  AliDebug(2,Form("Stream version: %s", raw->IsA()->GetName()));

  while (raw->Next()) {

    GetRawsData(0)->Fill(raw->GetDet());

    // possibly needs changes with the new reader !!
    Int_t *sig = raw->GetSignals();
    for(Int_t i=0; i<3; i++) GetRawsData(1)->Fill(sig[i]);
    // ---

    GetRawsData(2)->Fill(raw->GetTimeBin());

    // calculate the index;
    Int_t sm = raw->GetSM();
    Int_t roc = raw->GetROC();
    Int_t rob = raw->GetROB();
    Int_t mcm = raw->GetMCM();
    //Int_t adc = raw->GetADC();

    //Int_t index = roc * (kROB*kMCM*kADC) + rob * (kMCM*kADC) + mcm * kADC + adc;
    Int_t  index = roc * (kROB*kMCM) + rob * kMCM + mcm;
    GetRawsData(3)->Fill(sm);
    GetRawsData(4+sm)->Fill(index);
  }

  delete raw;
  */
}

//____________________________________________________________________________
void AliTRDQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  //  
  // Makes data from RecPoints
  // 
  
  //  Info("MakeRecPoints", "making");

  Int_t nsize = Int_t(clustersTree->GetTotBytes() / (sizeof(AliTRDcluster))); 
  TObjArray *clusterArray = new TObjArray(nsize+1000); 

  TBranch *branch = clustersTree->GetBranch("TRDcluster");
  if (!branch) {
    AliError("Can't get the branch !");
    return;
  }
  branch->SetAddress(&clusterArray); 

  // Loop through all entries in the tree
  Int_t nEntries   = (Int_t)TMath::Ceil( clustersTree->GetEntries() );
  Int_t nbytes     = 0;
  AliTRDcluster *c = 0;
  Int_t nDet[540];
  for (Int_t i=0; i<540; i++) nDet[i] = 0;

  Int_t nCls = 0;
  
  //printf("nEntries = %d\n", nEntries);
  //nEntries++;
  
  /*
  // select the event 
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {    

    // Import the tree
    nbytes += clustersTree->GetEvent(iEntry);  
    Int_t nCluster = clusterArray->GetEntries();  
    
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) { 
      c = (AliTRDcluster *) clusterArray->At(iCluster);
      nCls++;
    }
  }

  if (nCls < 100) {
    delete clusterArray;
    return;
  }
  */

  /////

  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {    

    // Import the tree
    nbytes += clustersTree->GetEvent(iEntry);  

    // Get the number of points in the detector
    Int_t nCluster = clusterArray->GetEntries();  

  
    // Loop through all TRD digits
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) { 

      nCls++;
      c = (AliTRDcluster *) clusterArray->At(iCluster);

      Int_t iDet = c->GetDetector();
      Int_t nPads = c->GetNPads();

      nDet[iDet]++;
      GetRecPointsData(0)->Fill(iDet);
      GetRecPointsData(86)->Fill(c->GetQ());
      GetRecPointsData(1)->Fill(iDet, c->GetQ());
      GetRecPointsData(2)->Fill(nPads);
      if (nPads < 6)
	GetRecPointsData(1+c->GetNPads())->Fill(c->GetCenter());
      
      if (nPads < 10)
	GetRecPointsData(88+nPads-2)->Fill(c->GetQ());
      else GetRecPointsData(96)->Fill(c->GetQ());

      //if (c->GetPadTime() < 5)
      ((TH2D*)GetRecPointsData(7))->Fill(c->GetPadRow(), c->GetPadCol());
      GetRecPointsData(8)->Fill(c->GetPadTime());

      ((TH3D*)GetRecPointsData(10))->Fill(iDet, c->GetPadTime(), c->GetQ());
      
      Int_t iSM = iDet / 30;
      GetRecPointsData(50+iSM)->Fill(c->GetPadTime());
      GetRecPointsData(68+iSM)->Fill(c->GetPadTime(), c->GetQ());
      
      // total charge sm / det / timeBin
      //((TH2D*)GetRecPointsData(14+iSM))->Fill(iDet-iSM*30, c->GetPadTime(), c->GetQ());


      // PRF for 2pad
      //if (c->GetNPads() == 2) {
      Short_t *sig = c->GetSignals();
      Double_t frac = -10;

      if (sig[0] == 0 && sig[1] == 0 && sig[2] == 0 && sig[5] == 0 && sig[6] == 0) 
	frac = 1. * sig[4] / (sig[3] + sig[4]);

      if (sig[0] == 0 && sig[1] == 0 && sig[4] == 0 && sig[5] == 0 && sig[6] == 0)
	frac = -1. * sig[2] / (sig[2] + sig[3]);

      if (frac > -10)  ((TProfile*)GetRecPointsData(11))->Fill(c->GetCenter(), frac);
	
      //}
    }
  }

  /*
  for(Int_t i=0; i<540; i++) 
    if (nDet[i] > 0) GetRecPointsData(9)->Fill(nDet[i]);
  */
  GetRecPointsData(9)->Fill(nCls);
  
  delete clusterArray;

}

//____________________________________________________________________________ 
void AliTRDQADataMakerRec::StartOfDetectorCycle()
{
  //
  // Detector specific actions at start of cycle
  //

}

//__________________________________________________________________________
Int_t AliTRDQADataMakerRec::CheckPointer(TObject *obj, const char *name) 
{
  //
  // Checks initialization of pointers
  //

  if (!obj) AliWarning(Form("null pointer: %s", name));
  return !!obj;
}
//__________________________________________________________________________
void AliTRDQADataMakerRec::BuildRatio(TH1D *ratio, TH1D *histN, TH1D*histD) {
  //
  // Calculate the ratio of two histograms 
  // error are calculated assuming the histos have the same counts
  //

  // calclate

  Int_t nbins = histN->GetXaxis()->GetNbins();
  for(Int_t i=1; i<nbins+2; i++) {
    
    Double_t valueN = histN->GetBinContent(i);
    Double_t valueD = histD->GetBinContent(i);
    
    if (valueD < 1) {
      ratio->SetBinContent(i, 0);
      ratio->SetBinError(i, 0);
      continue;
    }

    Double_t eps = (valueN < valueD-valueN)? valueN : valueD-valueN;
    
    ratio->SetBinContent(i, valueN/valueD);
    ratio->SetBinError(i, TMath::Sqrt(eps)/valueD);
  }

  // style
  ratio->SetMinimum(-0.1);
  ratio->SetMaximum(1.1);
  ratio->SetMarkerStyle(20);
}
//__________________________________________________________________________

Int_t AliTRDQADataMakerRec::FillBits(TH1D *hist, Int_t code, Int_t offset) {

  Int_t nb = 0;
  UInt_t test = 1;
  for(Int_t i=0; i<8; i++) {
    if (code & test) {
      hist->Fill(i+offset);
      nb++;
    }
    test *= 2;       
  }
  
  return nb;
}

//__________________________________________________________________________
