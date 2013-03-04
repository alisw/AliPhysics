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


//Root includes
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TAxis.h>
//AliRoot includes
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCROC.h"
#include "AliMathBase.h"
#include "TTreeStream.h"

//date
#include "event.h"

//header file
#include "AliTPCCalibKr.h"

//----------------------------------------------------------------------------
// The AliTPCCalibKr class description (TPC Kr calibration).
//
//
// The AliTPCCalibKr keeps the array of TH3F histograms (TPC_max_padraw,TPC_max_pad,TPC_ADC_cluster),
// its data memebers and is filled by AliTPCCalibKrTask under conditions (Accept()).   
// 
// The ouput TH3F histograms are later used to determine the calibration parameters of TPC chambers.
// These calculations are done by using AliTPCCalibKr::Analyse() function. The ouput calibration 
// parameters (details in AliTPCCalibKr::Analyse()) are stored in the calibKr.root file for each TPC pad.
// In addition the debugCalibKr.root file with debug information is created.
//

/*
 
// Usage example:
//

// 1. Analyse output histograms:
TFile f("outHistFile.root");
AliTPCCalibKr *obj = (AliTPCCalibKr*) cOutput.FindObject("AliTPCCalibKr");
obj->SetRadius(0,0);
obj->SetStep(1,1);
obj->Analyse();

// 2. See calibration parameters e.g.:
TFile f("calibKr.root");
spectrMean->GetCalROC(70)->GetValue(40,40);
fitMean->GetCalROC(70)->GetValue(40,40);

// 3. See debug information e.g.:
TFile f("debugCalibKr.root");
.ls;

// -- Print calibKr TTree content 
calibKr->Print();

// -- Draw calibKr TTree variables
calibKr.Draw("fitMean");

*/


//
// Author: Jacek Otwinowski (J.Otwinowski@gsi.de) and Stafan Geartner (S.Gaertner@gsi.de)
//-----------------------------------------------------------------------------

ClassImp(AliTPCCalibKr)

AliTPCCalibKr::AliTPCCalibKr() : 
  TObject(),
  fASide(kTRUE),
  fCSide(kTRUE),
  fHistoKrArray(72),
  fADCOverClustSizeMin(0.0),
  fADCOverClustSizeMax(1.0e9),
  fMaxADCOverClustADCMin(0.0),
  fMaxADCOverClustADCMax(1.0e9),
  fTimeMin(0.0),
  fTimeMax(1.0e9),
  fClustSizeMin(0.0),
  fClustSizeMax(1.0e9),
  fTimebinRmsIrocMin(0.0),
  fPadRmsIrocMin(0.0),
  fRowRmsIrocMin(0.0),
  fClusterPadSize1DIrocMax(200),
  fCurveCoefficientIroc(1.0e9),
  fTimebinRmsOrocMin(0.0),
  fPadRmsOrocMin(0.0),
  fRowRmsOrocMin(0.0),
  fClusterPadSize1DOrocMax(200),
  fCurveCoefficientOroc(1.0e9),
  fIrocHistogramMin(100.),
  fIrocHistogramMax(6000.),
  fIrocHistogramNbins(200),
  fOrocHistogramMin(100.),
  fOrocHistogramMax(5500.),
  fOrocHistogramNbins(200),
  fRowRadius(0),
  fPadRadius(0),
  fRowStep(1),
  fPadStep(1)

{
  //
  // default constructor
  //

  // TObjArray with histograms
  fHistoKrArray.SetOwner(kTRUE); // is owner of histograms
  fHistoKrArray.Clear();
 
  // init cuts (by Stefan)
//  SetADCOverClustSizeRange(7,200);
//  SetMaxADCOverClustADCRange(0.01,0.4);
//  SetTimeRange(200,600);
//  SetClustSizeRange(6,200);

  //init  cuts (by Adam)
  SetTimebinRmsMin(1.6,0.8);
  SetPadRmsMin(0.825,0.55);
  SetRowRmsMin(0.1,0.1);
  SetClusterPadSize1DMax(15,11);
  SetCurveCoefficient(1.,2.);
 
  //set histograms settings
  SetIrocHistogram(200,100,6000);
  SetOrocHistogram(200,100,5500);

  // init histograms
  //Init();
}

//_____________________________________________________________________
AliTPCCalibKr::AliTPCCalibKr(const AliTPCCalibKr& pad) : 
  TObject(pad),
  
  fASide(pad.fASide),
  fCSide(pad.fCSide),
  fHistoKrArray(72),
  fADCOverClustSizeMin(pad.fADCOverClustSizeMin),
  fADCOverClustSizeMax(pad.fADCOverClustSizeMax),
  fMaxADCOverClustADCMin(pad.fMaxADCOverClustADCMin),
  fMaxADCOverClustADCMax(pad.fMaxADCOverClustADCMax),
  fTimeMin(pad.fTimeMin),
  fTimeMax(pad.fTimeMax),
  fClustSizeMin(pad.fClustSizeMin),
  fClustSizeMax(pad.fClustSizeMax),
  fTimebinRmsIrocMin(pad.fTimebinRmsIrocMin),
  fPadRmsIrocMin(pad.fPadRmsIrocMin),
  fRowRmsIrocMin(pad.fRowRmsIrocMin),
  fClusterPadSize1DIrocMax(pad.fClusterPadSize1DIrocMax),
  fCurveCoefficientIroc(pad.fCurveCoefficientIroc),
  fTimebinRmsOrocMin(pad.fTimebinRmsOrocMin),
  fPadRmsOrocMin(pad.fPadRmsOrocMin),
  fRowRmsOrocMin(pad.fRowRmsOrocMin),
  fClusterPadSize1DOrocMax(pad.fClusterPadSize1DOrocMax),
  fCurveCoefficientOroc(pad.fCurveCoefficientOroc),
  fIrocHistogramMin(pad.fIrocHistogramMin),
  fIrocHistogramMax(pad.fIrocHistogramMax),
  fIrocHistogramNbins(pad.fIrocHistogramNbins),
  fOrocHistogramMin(pad.fOrocHistogramMin),
  fOrocHistogramMax(pad.fOrocHistogramMax),
  fOrocHistogramNbins(pad.fOrocHistogramNbins),
  fRowRadius(pad.fRowRadius),
  fPadRadius(pad.fPadRadius),
  fRowStep(pad.fRowStep),
  fPadStep(pad.fPadStep)

{
  // copy constructor
 
  // TObjArray with histograms
  fHistoKrArray.SetOwner(kTRUE); // is owner of histograms
  fHistoKrArray.Clear();

  for (Int_t iSec = 0; iSec < 72; ++iSec) 
  {
    TH3F *hOld = pad.GetHistoKr(iSec);
	if(hOld) {
      TH3F *hNew = new TH3F( *pad.GetHistoKr(iSec) ); 
      fHistoKrArray.AddAt(hNew,iSec);
	}
  }
}

//_____________________________________________________________________
AliTPCCalibKr::~AliTPCCalibKr() 
{
  //
  // destructor
  //

  // for (Int_t iSec = 0; iSec < 72; ++iSec) {
  //     if (fHistoKrArray.At(iSec)) delete fHistoKrArray.RemoveAt(iSec);
  // }
  fHistoKrArray.Delete();
}

//_____________________________________________________________________
AliTPCCalibKr& AliTPCCalibKr::operator = (const  AliTPCCalibKr &source)
{
  // assignment operator

  if (&source == this) return *this;
  new (this) AliTPCCalibKr(source);

  return *this;
}

//_____________________________________________________________________
void AliTPCCalibKr::Init()
{
  // 
  // init output histograms 
  //

  // add histograms to the TObjArray
  for(Int_t i=0; i<72; ++i) {
    
    // C - side
    if(IsCSide(i) == kTRUE && fCSide == kTRUE) {
      TH3F *hist = CreateHisto(i);
      if(hist) fHistoKrArray.AddAt(hist,i);
    }
    
    // A - side
    if(IsCSide(i) == kFALSE && fASide == kTRUE) {
      TH3F *hist = CreateHisto(i);
      if(hist) fHistoKrArray.AddAt(hist,i);
    }
  }
}
 
//_____________________________________________________________________
Bool_t AliTPCCalibKr::Process(AliTPCclusterKr *cluster)
{
  //
  // process events 
  // call event by event
  //

  if(cluster) Update(cluster);
  else return kFALSE;
 
 return kTRUE;
}

//_____________________________________________________________________
TH3F* AliTPCCalibKr::CreateHisto(Int_t chamber)
{
    //
    // create new histogram
	//
    char name[256];
	TH3F *h;

	snprintf(name,256,"ADCcluster_ch%d",chamber);

    if( IsIROC(chamber) == kTRUE ) 
	{
	  h = new TH3F(name,name,63,0,63,108,0,108,fIrocHistogramNbins,fIrocHistogramMin,fIrocHistogramMax);
	} else {
	  h = new TH3F(name,name,96,0,96,140,0,140,fOrocHistogramNbins,fOrocHistogramMin,fOrocHistogramMax);
	}
    h->SetXTitle("padrow");
    h->SetYTitle("pad");
    h->SetZTitle("fADC");

return h;
}

//_____________________________________________________________________
Bool_t AliTPCCalibKr::IsIROC(Int_t chamber)
{
// check if IROCs
// returns kTRUE if IROCs and kFALSE if OROCs 

  if(chamber>=0 && chamber<36) return kTRUE;

return kFALSE;
}

//_____________________________________________________________________
Bool_t AliTPCCalibKr::IsCSide(Int_t chamber)
{
// check if C side
// returns kTRUE if C side and kFALSE if A side

  if((chamber>=18 && chamber<36) || (chamber>=54 && chamber<72)) return kTRUE;

return kFALSE;
}
 
//_____________________________________________________________________
Bool_t AliTPCCalibKr::Update(AliTPCclusterKr  *cl)
{
  //
  // fill existing histograms
  //

  if (!Accept(cl)) return kFALSE;
  TH3F *h = (TH3F*)fHistoKrArray.At(cl->GetSec());
  if(!h) return kFALSE;
  
  h->Fill(cl->GetMax().GetRow(),cl->GetMax().GetPad(),cl->GetADCcluster());
  
  return kTRUE;
}

//_____________________________________________________________________
Bool_t AliTPCCalibKr::Accept(AliTPCclusterKr  *cl){
  //
  // cuts
  //
  /*
    TCut cutR0("cutR0","fADCcluster/fSize<200");        // adjust it according v seetings - 
    TCut cutR1("cutR1","fADCcluster/fSize>7");          // cosmic tracks and noise removal
    TCut cutR2("cutR2","fMax.fAdc/fADCcluster<0.4");    // digital noise removal
    TCut cutR3("cutR3","fMax.fAdc/fADCcluster>0.01");   // noise removal
    TCut cutR4("cutR4","fMax.fTime>200");   // noise removal
    TCut cutR5("cutR5","fMax.fTime<600");   // noise removal
    TCut cutS1("cutS1","fSize<200");    // adjust it according v seetings - cosmic tracks
    TCut cutAll = cutR0+cutR1+cutR2+cutR3+cutR4+cutR5+cutS1;
  */
  /*
  //R0
  if ((float)cl->GetADCcluster()/ cl->GetSize() >200)        return kFALSE;
  // R1
  if ((float)cl->GetADCcluster()/ cl->GetSize() <7)          return kFALSE;
  //R2
  if ((float)cl->GetMax().GetAdc()/ cl->GetADCcluster() >0.4)  return kFALSE;
  //R3
  if ((float)cl->GetMax().GetAdc()/ cl->GetADCcluster() <0.01) return kFALSE;
  //R4
  if (cl->GetMax().GetTime() < 200) return kFALSE;
  //R5
  if (cl->GetMax().GetTime() > 600) return kFALSE;
  //S1
  if (cl->GetSize()>200) return kFALSE;
  if (cl->GetSize()<6)  return kFALSE;

  SetADCOverClustSizeRange(7,200);
  SetMaxADCOverClustADCRange(0.01,0.4);
  SetTimeRange(200,600);
  SetClustSizeRange(6,200);
  */

  //R0
  if ((float)cl->GetADCcluster()/ cl->GetSize() >fADCOverClustSizeMax)        return kFALSE;
  // R1
  if ((float)cl->GetADCcluster()/ cl->GetSize() <fADCOverClustSizeMin)          return kFALSE;
  //R2
  if ((float)cl->GetMax().GetAdc()/ cl->GetADCcluster() >fMaxADCOverClustADCMax)  return kFALSE;
  //R3
  if ((float)cl->GetMax().GetAdc()/ cl->GetADCcluster() <fMaxADCOverClustADCMin) return kFALSE;
  //R4
  if (cl->GetMax().GetTime() > fTimeMax) return kFALSE;
  //R5
  if (cl->GetMax().GetTime() < fTimeMin) return kFALSE;
  //S1
  if (cl->GetSize()>fClustSizeMax) return kFALSE;
  if (cl->GetSize()<fClustSizeMin)  return kFALSE;

  //
  // cuts by Adam
  //
  /*
    TCut cutAI0("cutAI0","fTimebinRMS>1.6");
    TCut cutAI1("cutAI1","fPadRMS>0.825");  
    TCut cutAI2("cutAI2","fRowRMS>0.1");    
    TCut cutAI3("cutAI3","fPads1D<15");  
    TCut cutAI4("cutAI4","fTimebinRMS+11.9-2.15*TMath::Log(1.*fADCcluster)<0"); 

    TCut cutAIAll = cutAI0+cutAI1+cutAI2+cutAI3+cutAI4;

    TCut cutAO0("cutAO0","fTimebinRMS>0.8");
    TCut cutAO1("cutAO1","fPadRMS>0.55");  
    TCut cutAO2("cutAO2","fRowRMS>0.1");    
    TCut cutAO3("cutAO3","fPads1D<11");  
    TCut cutAO4("cutAO4","fTimebinRMS+11.9-2.15*TMath::Log(2.*fADCcluster)<0"); 

    TCut cutAOAll = cutAO0+cutAO1+cutAO2+cutAO3+cutAO4;
    TCut cutAll("cutAll","(fSec<36&&cutAIAll)||(fSec>=36&&cutAOAll)");

  */
  if(cl->GetSec()<36){  //IROCs
    //AI0
    if((float)cl->GetTimebinRMS() <= fTimebinRmsIrocMin) return kFALSE;
    //AI1
    if((float)cl->GetPadRMS() <= fPadRmsIrocMin) return kFALSE;
    //AI2
    if((float)cl->GetRowRMS() <= fRowRmsIrocMin) return kFALSE;
    //AI3
    if(cl->GetPads1D() >= fClusterPadSize1DIrocMax) return kFALSE;
    //AI4
    if((float)cl->GetTimebinRMS()+11.9-2.15*TMath::Log(fCurveCoefficientIroc*(float)cl->GetADCcluster()) >= 0) return kFALSE;
  }else{  //OROCs
    //AO0
    if((float)cl->GetTimebinRMS() <= fTimebinRmsOrocMin) return kFALSE;
    //AO1
    if((float)cl->GetPadRMS() <= fPadRmsOrocMin) return kFALSE;
    //AO2
    if((float)cl->GetRowRMS() <= fRowRmsOrocMin) return kFALSE;
    //AO3
    if(cl->GetPads1D() >= fClusterPadSize1DOrocMax) return kFALSE;
    //AO4
    if((float)cl->GetTimebinRMS()+11.9-2.15*TMath::Log(fCurveCoefficientOroc*(float)cl->GetADCcluster()) >= 0) return kFALSE;
  }

  return kTRUE;

}

//_____________________________________________________________________
TH3F* AliTPCCalibKr::GetHistoKr(Int_t chamber) const
{
  // get histograms from fHistoKrArray
  return (TH3F*) fHistoKrArray.At(chamber);
}
 
//_____________________________________________________________________
void AliTPCCalibKr::Analyse() 
{
  //
  // analyse the histograms and extract krypton calibration parameters
  //

  // AliTPCCalPads that will contain the calibration parameters
  AliTPCCalPad* spectrMeanCalPad = new AliTPCCalPad("spectrMean", "spectrMean");
  AliTPCCalPad* spectrRMSCalPad = new AliTPCCalPad("spectrRMS", "spectrRMS");
  AliTPCCalPad* fitMeanCalPad = new AliTPCCalPad("fitMean", "fitMean");
  AliTPCCalPad* fitRMSCalPad = new AliTPCCalPad("fitRMS", "fitRMS");
  AliTPCCalPad* fitNormChi2CalPad = new AliTPCCalPad("fitNormChi2", "fitNormChi2");
  AliTPCCalPad* entriesCalPad = new AliTPCCalPad("entries", "entries");

  // file stream for debugging purposes
  TTreeSRedirector* debugStream = new TTreeSRedirector("debugCalibKr.root");

  // if entries in spectrum less than minEntries, then the fit won't be performed
  Int_t minEntries = 1; //300;

  Double_t windowFrac = 0.12;
  // the 3d histogram will be projected on the pads given by the following window size
  // set the numbers to 0 if you want to do a pad-by-pad calibration
  UInt_t rowRadius = fRowRadius;//4;
  UInt_t padRadius = fPadRadius;//4;
  
  // the step size by which pad and row are incremented is given by the following numbers
  // set them to 1 if you want the finest granularity
  UInt_t rowStep = fRowStep;//1;//2    // formerly: 2*rowRadius
  UInt_t padStep = fPadStep;//1;//4    // formerly: 2*padRadius

  for (Int_t chamber = 0; chamber < 72; chamber++) {
    //if (chamber != 71) continue;
    AliTPCCalROC roc(chamber);  // I need this only for GetNrows() and GetNPads()
    
    // Usually I would traverse each pad, take the spectrum for its neighbourhood and
    // obtain the calibration parameters. This takes very long, so instead I assign the same
    // calibration values to the whole neighbourhood and then go on to the next neighbourhood.
    UInt_t nRows = roc.GetNrows();
    for (UInt_t iRow = 0; iRow < nRows; iRow += rowStep) {
      UInt_t nPads = roc.GetNPads(iRow);
      //if (iRow >= 10) break;
      for (UInt_t iPad = 0; iPad < nPads; iPad += padStep) {
        //if (iPad >= 20) break;
        TH3F* h = GetHistoKr(chamber);
        if (!h) continue;
        
        // the 3d histogram will be projected on the pads given by the following bounds
        // for rows and pads
        Int_t rowLow = iRow - rowRadius;
        UInt_t rowUp = iRow + rowRadius + rowStep-1;
        Int_t padLow = iPad - padRadius;
        UInt_t padUp = iPad + padRadius + padStep-1;
        // if window goes out of chamber
        if (rowLow < 0) rowLow = 0;
        if (rowUp >= nRows) rowUp = nRows - 1;
        if (padLow < 0) padLow = 0;
        if (padUp >= nPads) padUp = nPads - 1;

        // project the histogram
        //TH1D* projH = h->ProjectionZ("projH", rowLow+1, rowUp+1, padLow+1, padUp+1); // SLOW
        TH1D* projH = ProjectHisto(h, "projH", rowLow+1, rowUp+1, padLow+1, padUp+1);
    
        // get the number of entries in the spectrum
        Double_t entries = projH->GetEntries();
        if (entries < minEntries) { delete projH; continue; }
        
        // get the two calibration parameters mean of spectrum and RMS of spectrum
        Double_t histMean = projH->GetMean();
        Double_t histRMS = (histMean != 0) ? projH->GetRMS() / histMean : 0.;
    
        // find maximum in spectrum to define a range (given by windowFrac) for which a Gauss is fitted
        Double_t maxEntries = projH->GetBinCenter(projH->GetMaximumBin());
		Int_t minBin = projH->FindBin((1.-windowFrac) * maxEntries);
		Int_t maxBin = projH->FindBin((1.+windowFrac) * maxEntries);
        Double_t integCharge =  projH->Integral(minBin,maxBin); 

        Int_t fitResult = projH->Fit("gaus", "Q0", "", (1.-windowFrac) * maxEntries, (1.+windowFrac) * maxEntries);

        if (fitResult != 0) {
          Error("Analyse", "Error while fitting spectrum for chamber %i, rows %i - %i, pads %i - %i, integrated charge %f.", chamber, rowLow, rowUp, padLow, padUp, integCharge);
          //printf("[ch%i r%i, p%i] entries = %f, maxEntries = %f,  %f\n", chamber, iRow, iPad, entries, maxEntries);
    
          delete projH;
          continue;
        }
    
        // get the two calibration parameters mean of gauss fit and sigma of gauss fit
        TF1* gausFit = projH->GetFunction("gaus");
        Double_t fitMean = gausFit->GetParameter(1);
        Double_t fitRMS = gausFit->GetParameter(2);
        Int_t numberFitPoints = gausFit->GetNumberFitPoints();
        if (numberFitPoints == 0) continue;
        Double_t fitNormChi2 = gausFit->GetChisquare() / numberFitPoints;
        delete gausFit;
        delete projH;
        if (fitMean <= 0) continue;
        //printf("[ch%i r%i, p%i] entries = %f, maxEntries = %f, fitMean = %f, fitRMS = %f\n", chamber, iRow, iPad, entries, maxEntries, fitMean, fitRMS);
    
        // write the calibration parameters for each pad that the 3d histogram was projected onto
        // (with considering the step size) to the CalPads
        // rowStep (padStep) odd: round down s/2 and fill this # of rows (pads) in both directions
        // rowStep (padStep) even: fill s/2 rows (pads) in ascending direction, s/2-1 in descending direction
        for (Int_t r = iRow - (rowStep/2 - (rowStep+1)%2); r <= (Int_t)(iRow + rowStep/2); r++) {
          if (r < 0 || r >= (Int_t)nRows) continue;
          UInt_t nPadsR = roc.GetNPads(r);
          for (Int_t p = iPad - (padStep/2 - (padStep+1)%2); p <= (Int_t)(iPad + padStep/2); p++) {
            if (p < 0 || p >= (Int_t)nPadsR) continue;
            spectrMeanCalPad->GetCalROC(chamber)->SetValue(r, p, histMean);
            spectrRMSCalPad->GetCalROC(chamber)->SetValue(r, p, histRMS);
            fitMeanCalPad->GetCalROC(chamber)->SetValue(r, p, fitMean);
            fitRMSCalPad->GetCalROC(chamber)->SetValue(r, p, fitRMS);
            fitNormChi2CalPad->GetCalROC(chamber)->SetValue(r, p, fitNormChi2);
            entriesCalPad->GetCalROC(chamber)->SetValue(r, p, entries);

            (*debugStream) << "calibKr" <<
              "sector=" << chamber <<          // chamber number
              "row=" << r <<                   // row number
              "pad=" << p <<                   // pad number
              "histMean=" << histMean <<       // mean of the spectrum
              "histRMS=" << histRMS <<         // RMS of the spectrum divided by the mean
              "fitMean=" << fitMean <<         // Gauss fitted mean of the 41.6 keV Kr peak
              "fitRMS=" << fitRMS <<           // Gauss fitted sigma of the 41.6 keV Kr peak
              "fitNormChi2" << fitNormChi2 <<  // normalized chi square of the Gauss fit
              "entries=" << entries <<         // number of entries for the spectrum
              "\n";
          }
        }
      }
    }
  }

  TFile f("calibKr.root", "recreate");
  spectrMeanCalPad->Write();
  spectrRMSCalPad->Write();
  fitMeanCalPad->Write();
  fitRMSCalPad->Write();
  fitNormChi2CalPad->Write();
  entriesCalPad->Write();
  f.Close();
  delete spectrMeanCalPad;
  delete spectrRMSCalPad;
  delete fitMeanCalPad;
  delete fitRMSCalPad;
  delete fitNormChi2CalPad;
  delete entriesCalPad;
  delete debugStream;
}

//_____________________________________________________________________
TH1D* AliTPCCalibKr::ProjectHisto(TH3F* histo3D, const char* name, Int_t xMin, Int_t xMax, Int_t yMin, Int_t yMax)
{
  // project the z-axis of a 3d histo to a specific range of the x- and y-axes,
  // replaces TH3F::ProjectZ() to gain more speed

  TAxis* xAxis = histo3D->GetXaxis();
  TAxis* yAxis = histo3D->GetYaxis();
  TAxis* zAxis = histo3D->GetZaxis();
  Double_t zMinVal = zAxis->GetXmin();
  Double_t zMaxVal = zAxis->GetXmax();
  
  Int_t nBinsZ = zAxis->GetNbins();
  TH1D* projH = new TH1D(name, name, nBinsZ, zMinVal, zMaxVal);

  Int_t nx = xAxis->GetNbins()+2;
  Int_t ny = yAxis->GetNbins()+2;
  Int_t bin = 0;
  Double_t entries = 0.;
  for (Int_t x = xMin; x <= xMax; x++) {
    for (Int_t y = yMin; y <= yMax; y++) {
      for (Int_t z = 0; z <= nBinsZ+1; z++) {
        bin = x + nx * (y + ny * z);
        Double_t val = histo3D->GetBinContent(bin);
        projH->Fill(zAxis->GetBinCenter(z), val);
        entries += val;
      }
    }
  }
  projH->SetEntries((Long64_t)entries);
  return projH;
}

//_____________________________________________________________________
Long64_t AliTPCCalibKr::Merge(TCollection* list) {
// merge function 
//

if (!list)
return 0;

if (list->IsEmpty())
return 1;

TIterator* iter = list->MakeIterator();
TObject* obj = 0;

    iter->Reset();
    Int_t count=0;
	while((obj = iter->Next()) != 0)
	{
	  AliTPCCalibKr* entry = dynamic_cast<AliTPCCalibKr*>(obj);
	  if (entry == 0) continue; 

		for(int i=0; i<72; ++i) { 
		  if(IsCSide(i) == kTRUE && fCSide == kTRUE) { 
		  ((TH3F*)fHistoKrArray.At(i))->Add( ((TH3F*)entry->fHistoKrArray.At(i)) );  
		  } 

		  if(IsCSide(i) == kFALSE && fASide == kTRUE) { 
		    ((TH3F*)fHistoKrArray.At(i))->Add( ((TH3F*)entry->fHistoKrArray.At(i)) );  
		  } 
		} 

	  count++;
	}

return count;
}
