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

/////////////////////////////////////////////////////////////////////////////////
//                                                                             
// AliTRDCalibraFit                                                               
//                                                                             
// This class is for the TRD calibration of the relative gain factor, the drift velocity,
// the time 0 and the pad response function. It fits the histos.       
// The 2D histograms or vectors (first converted in 1D histos) will be fitted  
// if they have enough entries, otherwise the (default) value of the choosen database 
// will be put. For the relative gain calibration the resulted factors will be globally 
// normalized to the gain factors of the choosen database. It unables to precise     
// previous calibration procedure.
// The function SetDebug enables the user to see:                                     
// _fDebug = 0: nothing, only the values are written in the tree if wanted
// _fDebug = 1: only the fit of the choosen calibration group fFitVoir (SetFitVoir)
// _fDebug = 2: a comparaison of the coefficients found and the default values 
//              in the choosen database.
//              fCoef , histogram of the coefs as function of the calibration group number
//              fDelta , histogram of the relative difference of the coef with the default
//                        value in the database as function of the calibration group number
//              fError , dirstribution of this relative difference
// _fDebug = 3: The coefficients in the choosen detector fDet (SetDet) as function of the
//              pad row and col number
// _fDebug = 4; The coeffcicients in the choosen detector fDet (SetDet) like in the 3 but with
//              also the comparaison histograms of the 1 for this detector
//
//                            
// Author:
//   R. Bailhache (R.Bailhache@gsi.de)
//                            
//////////////////////////////////////////////////////////////////////////////////////

#include <TLine.h>
#include <TH1I.h>
#include <TStyle.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TAxis.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TTreeStream.h>
#include <TLinearFitter.h>
#include <TVectorD.h>
#include <TROOT.h>
#include <TString.h>

#include "AliLog.h"
#include "AliMathBase.h"

#include "AliTRDCalibraFit.h"
#include "AliTRDCalibraMode.h"
#include "AliTRDCalibraVector.h"
#include "AliTRDCalibraVdriftLinearFit.h"
#include "AliTRDcalibDB.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDgeometry.h"
#include "AliTRDCommonParam.h"
#include "./Cal/AliTRDCalROC.h"
#include "./Cal/AliTRDCalPad.h"
#include "./Cal/AliTRDCalDet.h"


ClassImp(AliTRDCalibraFit)

AliTRDCalibraFit* AliTRDCalibraFit::fgInstance = 0;
Bool_t AliTRDCalibraFit::fgTerminated = kFALSE;

//_____________singleton implementation_________________________________________________
AliTRDCalibraFit *AliTRDCalibraFit::Instance()
{
  //
  // Singleton implementation
  //

  if (fgTerminated != kFALSE) {
    return 0;
  }

  if (fgInstance == 0) {
    fgInstance = new AliTRDCalibraFit();
  }

  return fgInstance;

}
//______________________________________________________________________________________
void AliTRDCalibraFit::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class
  //

  fgTerminated = kTRUE;

  if (fgInstance != 0) {
    delete fgInstance;
    fgInstance = 0;
  }

}
//______________________________________________________________________________________
AliTRDCalibraFit::AliTRDCalibraFit()
  :TObject()
  ,fGeo(0)
  ,fNumberOfBinsExpected(0)
  ,fMethod(0)
  ,fBeginFitCharge(3.5)
  ,fFitPHPeriode(1)
  ,fTakeTheMaxPH(kTRUE)
  ,fT0Shift0(0.124797)
  ,fT0Shift1(0.267451)
  ,fRangeFitPRF(1.0)
  ,fAccCDB(kFALSE)
  ,fMinEntries(800)
  ,fRebin(1)
  ,fScaleGain(0.02431)
  ,fNumberFit(0)
  ,fNumberFitSuccess(0)
  ,fNumberEnt(0)
  ,fStatisticMean(0.0)
  ,fDebugStreamer(0x0)
  ,fDebugLevel(0)
  ,fFitVoir(0)
  ,fMagneticField(0.5)
  ,fCalibraMode(new AliTRDCalibraMode())
  ,fCurrentCoefE(0.0)
  ,fCurrentCoefE2(0.0)
  ,fDect1(0)
  ,fDect2(0)
  ,fScaleFitFactor(0.0)
  ,fEntriesCurrent(0)
  ,fCountDet(0)
  ,fCount(0)
  ,fNbDet(0)
  ,fCalDet(0x0)
  ,fCalROC(0x0)
  ,fCalDet2(0x0)
  ,fCalROC2(0x0)
  ,fCalDetVdriftUsed(0x0)
  ,fCalDetExBUsed(0x0)
  ,fCurrentCoefDetector(0x0)
  ,fCurrentCoefDetector2(0x0)
  ,fVectorFit(0)
  ,fVectorFit2(0)
{
  //
  // Default constructor
  //

  fGeo         = new AliTRDgeometry();  
 
  // Current variables initialised
  for (Int_t k = 0; k < 2; k++) {
    fCurrentCoef[k]  = 0.0;
    fCurrentCoef2[k] = 0.0;
  }
  for (Int_t i = 0; i < 3; i++) {
    fPhd[i]          = 0.0;
    fDet[i]          = 0;
  }
 
}
//______________________________________________________________________________________
AliTRDCalibraFit::AliTRDCalibraFit(const AliTRDCalibraFit &c)
:TObject(c)
,fGeo(0)
,fNumberOfBinsExpected(c.fNumberOfBinsExpected)
,fMethod(c.fMethod)
,fBeginFitCharge(c.fBeginFitCharge)
,fFitPHPeriode(c.fFitPHPeriode)
,fTakeTheMaxPH(c.fTakeTheMaxPH)
,fT0Shift0(c.fT0Shift0)
,fT0Shift1(c.fT0Shift1)
,fRangeFitPRF(c.fRangeFitPRF)
,fAccCDB(c.fAccCDB)
,fMinEntries(c.fMinEntries)
,fRebin(c.fRebin)
,fScaleGain(c.fScaleGain)
,fNumberFit(c.fNumberFit)
,fNumberFitSuccess(c.fNumberFitSuccess)
,fNumberEnt(c.fNumberEnt)
,fStatisticMean(c.fStatisticMean)
,fDebugStreamer(0x0)
,fDebugLevel(c.fDebugLevel)
,fFitVoir(c.fFitVoir)
,fMagneticField(c.fMagneticField)
,fCalibraMode(0x0)
,fCurrentCoefE(c.fCurrentCoefE)
,fCurrentCoefE2(c.fCurrentCoefE2)
,fDect1(c.fDect1)
,fDect2(c.fDect2)
,fScaleFitFactor(c.fScaleFitFactor)
,fEntriesCurrent(c.fEntriesCurrent)
,fCountDet(c.fCountDet)
,fCount(c.fCount)
,fNbDet(c.fNbDet)
,fCalDet(0x0)
,fCalROC(0x0)
,fCalDet2(0x0)
,fCalROC2(0x0)
,fCalDetVdriftUsed(0x0)
,fCalDetExBUsed(0x0)
,fCurrentCoefDetector(0x0)
,fCurrentCoefDetector2(0x0)
,fVectorFit(0)
,fVectorFit2(0)
{
  //
  // Copy constructor
  //

  if(c.fCalibraMode)   fCalibraMode = new AliTRDCalibraMode(*c.fCalibraMode);
 
  //Current variables initialised
  for (Int_t k = 0; k < 2; k++) {
    fCurrentCoef[k]  = 0.0;
    fCurrentCoef2[k] = 0.0;
  }
  for (Int_t i = 0; i < 3; i++) {
    fPhd[i]          = 0.0;
    fDet[i]          = 0;
  }
  if(c.fCalDet) fCalDet   = new AliTRDCalDet(*c.fCalDet);
  if(c.fCalDet2) fCalDet2 = new AliTRDCalDet(*c.fCalDet2);
  
  if(c.fCalROC) fCalROC   = new AliTRDCalROC(*c.fCalROC);
  if(c.fCalROC2) fCalROC  = new AliTRDCalROC(*c.fCalROC2);

  if(c.fCalDetVdriftUsed) fCalDetVdriftUsed = new AliTRDCalDet(*c.fCalDetVdriftUsed);
  if(c.fCalDetExBUsed)    fCalDetExBUsed = new AliTRDCalDet(*c.fCalDetExBUsed);

  fVectorFit.SetName(c.fVectorFit.GetName());
  for(Int_t k = 0; k < c.fVectorFit.GetEntriesFast(); k++){
    AliTRDFitInfo *fitInfo = new AliTRDFitInfo();
    Int_t detector         = ((AliTRDFitInfo *)c.fVectorFit.UncheckedAt(k))->GetDetector();
    Int_t ntotal = 1;
    if (GetStack(detector) == 2) {
      ntotal = 1728;
    }
    else {
      ntotal = 2304;
    }
    Float_t *coef = new Float_t[ntotal];
    for (Int_t i = 0; i < ntotal; i++) {
      coef[i] = ((AliTRDFitInfo *)c.fVectorFit.UncheckedAt(k))->GetCoef()[i];
    }
    fitInfo->SetCoef(coef);
    fitInfo->SetDetector(detector);
    fVectorFit.Add((TObject *) fitInfo);
  }
  fVectorFit.SetName(c.fVectorFit.GetName());
  for(Int_t k = 0; k < c.fVectorFit2.GetEntriesFast(); k++){
    AliTRDFitInfo *fitInfo = new AliTRDFitInfo();
    Int_t detector         = ((AliTRDFitInfo *)c.fVectorFit2.UncheckedAt(k))->GetDetector();
    Int_t ntotal = 1;
    if (GetStack(detector) == 2) {
      ntotal = 1728;
    }
    else {
      ntotal = 2304;
    }
    Float_t *coef = new Float_t[ntotal];
    for (Int_t i = 0; i < ntotal; i++) {
      coef[i] = ((AliTRDFitInfo *)c.fVectorFit2.UncheckedAt(k))->GetCoef()[i];
    }
    fitInfo->SetCoef(coef);
    fitInfo->SetDetector(detector);
    fVectorFit2.Add((TObject *) fitInfo);
  }
  if (fGeo) {
    delete fGeo;
  }
  fGeo = new AliTRDgeometry();

}
//____________________________________________________________________________________
AliTRDCalibraFit::~AliTRDCalibraFit()
{
  //
  // AliTRDCalibraFit destructor
  //
  if ( fDebugStreamer ) delete fDebugStreamer;
  if ( fCalDet )  delete fCalDet;
  if ( fCalDet2 ) delete fCalDet2;
  if ( fCalROC )  delete fCalROC;
  if ( fCalROC2 ) delete fCalROC2;
  if ( fCalDetVdriftUsed)  delete fCalDetVdriftUsed;
  if ( fCalDetExBUsed) delete fCalDetExBUsed;
  if( fCurrentCoefDetector ) delete [] fCurrentCoefDetector;
  if( fCurrentCoefDetector2 ) delete [] fCurrentCoefDetector2; 
  fVectorFit.Delete();
  fVectorFit2.Delete();
  if (fGeo) {
    delete fGeo;
  }

}
//_____________________________________________________________________________
void AliTRDCalibraFit::Destroy() 
{
  //
  // Delete instance 
  //

  if (fgInstance) {
    delete fgInstance;
    fgInstance = 0x0;
  }

}
//_____________________________________________________________________________
void AliTRDCalibraFit::DestroyDebugStreamer() 
{
  //
  // Delete DebugStreamer
  //

  if ( fDebugStreamer ) delete fDebugStreamer;
  fDebugStreamer = 0x0;
 
}
//__________________________________________________________________________________
void AliTRDCalibraFit::RangeChargeIntegration(Float_t vdrift, Float_t t0, Int_t &begin, Int_t &peak, Int_t &end) const
{
  //
  // From the drift velocity and t0
  // return the position of the peak and maximum negative slope
  //
  
  const Float_t kDrWidth = AliTRDgeometry::DrThick();    // drift region
  Double_t widbins = 0.1;                                // 0.1 mus

  //peak and maxnegslope in mus
  Double_t begind = t0*widbins + fT0Shift0;
  Double_t peakd  = t0*widbins + fT0Shift1;
  Double_t maxnegslope = (kDrWidth + vdrift*peakd)/vdrift; 

  // peak and maxnegslope in timebin
  begin = TMath::Nint(begind*widbins);
  peak  = TMath::Nint(peakd*widbins);
  end   = TMath::Nint(maxnegslope*widbins); 

}
//____________Functions fit Online CH2d________________________________________
Bool_t AliTRDCalibraFit::AnalyseCH(const TH2I *ch)
{
  //
  // Fit the 1D histos, projections of the 2D ch on the Xaxis, for each
  // calibration group normalized the resulted coefficients (to 1 normally)
  //

  // Set the calibration mode
  //const char *name = ch->GetTitle();
  TString name = ch->GetTitle();
  if(!SetModeCalibration(name,0)) return kFALSE;

  // Number of Ybins (detectors or groups of pads)
  Int_t    nbins   = ch->GetNbinsX();// charge
  Int_t    nybins  = ch->GetNbinsY();// groups number
  if (!InitFit(nybins,0)) {
    return kFALSE;
  }
  if (!InitFitCH()) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;
  // Init fCountDet and fCount
  InitfCountDetAndfCount(0);
  // Beginning of the loop betwwen dect1 and dect2
  for (Int_t idect = fDect1; idect < fDect2; idect++) {
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi...
    UpdatefCountDetAndfCount(idect,0);
    ReconstructFitRowMinRowMax(idect, 0);
    // Take the histo
    TH1I *projch = (TH1I *) ch->ProjectionX("projch",idect+1,idect+1,(Option_t *)"e");
    projch->SetDirectory(0);
    // Number of entries for this calibration group
    Double_t nentries = 0.0;
    Double_t mean = 0.0;
    for (Int_t k = 0; k < nbins; k++) {
      Int_t binnb = (nbins+2)*(idect+1)+(k+1);
      nentries += ch->GetBinContent(binnb);
      mean += projch->GetBinCenter(k+1)*projch->GetBinContent(k+1);
      projch->SetBinError(k+1,TMath::Sqrt(projch->GetBinContent(k+1)));
    }
    projch->SetEntries(nentries);
    //printf("The number of entries for the group %d is %f\n",idect,nentries);
    if (nentries > 0) {
      fNumberEnt++;
      mean /= nentries;
    }
    // Rebin and statistic stuff
    if (fRebin > 1) {
      projch = ReBin((TH1I *) projch);
    }
    // This detector has not enough statistics or was off
    if (nentries <= fMinEntries) {
      NotEnoughStatisticCH(idect);
      if (fDebugLevel != 1) {
	delete projch;
      }
      continue;
    }
    // Statistics of the group fitted
    fStatisticMean += nentries;
    fNumberFit++;
    //Method choosen
    switch(fMethod)
      {
      case 0: FitMeanW((TH1 *) projch, nentries); break;
      case 1: FitMean((TH1 *) projch, nentries, mean); break;
      case 2: FitCH((TH1 *) projch, mean); break;
      case 3: FitBisCH((TH1 *) projch, mean); break;
      default: return kFALSE;
      }
    // Fill Infos Fit
    FillInfosFitCH(idect);
    // Memory!!!
    if (fDebugLevel != 1) {
      delete projch;
    }
  } // Boucle object
  // Normierungcharge
  if (fDebugLevel != 1) {
    NormierungCharge();
  }
  // Mean Statistic
  if (fNumberFit > 0) {
    AliInfo(Form("There are %d with at least one entries. %d fits have been proceeded (sucessfully or not...). There is a mean statistic of: %d over these fitted histograms and %d successfulled fits",fNumberEnt, fNumberFit, (Int_t) fStatisticMean/fNumberFit, fNumberFitSuccess));
    fStatisticMean = fStatisticMean / fNumberFit;
  }
  else {
    AliInfo(Form("There are %d with at least one entries. There is no fit!",fNumberEnt));
  }
  delete fDebugStreamer;
  fDebugStreamer = 0x0;

  return kTRUE;
}
//____________Functions fit Online CH2d________________________________________
Bool_t AliTRDCalibraFit::AnalyseCH(AliTRDCalibraVector *calvect)
{
  //
  // Reconstruct a 1D histo from the vectorCH for each calibration group,
  // fit the histo, normalized the resulted coefficients (to 1 normally)
  //

  // Set the calibraMode
  //const char *name = calvect->GetNameCH();
  TString name = calvect->GetNameCH();
  if(!SetModeCalibration(name,0)) return kFALSE;  

  // Number of Xbins (detectors or groups of pads)
  if (!InitFit((432*calvect->GetDetCha0(0)+108*calvect->GetDetCha2(0)),0)) {
    return kFALSE;
  }
  if (!InitFitCH()) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;
  // Init fCountDet and fCount
  InitfCountDetAndfCount(0);
  // Beginning of the loop between dect1 and dect2
  for (Int_t idect = fDect1; idect < fDect2; idect++) {
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi...........
    UpdatefCountDetAndfCount(idect,0);
    ReconstructFitRowMinRowMax(idect,0);
    // Take the histo
    Double_t nentries = 0.0;
    Double_t mean = 0.0;
    if(!calvect->GetCHEntries(fCountDet)) {
      NotEnoughStatisticCH(idect);
      continue;
    }
    
    TString tname("CH");
    tname += idect;
    TH1F *projch  = calvect->ConvertVectorCHHisto(fCountDet,(idect-(fCount-(fCalibraMode->GetNfragZ(0)*fCalibraMode->GetNfragRphi(0)))),(const char *) tname);
    projch->SetDirectory(0);
    for (Int_t k = 0; k < calvect->GetNumberBinCharge(); k++) {
      nentries += projch->GetBinContent(k+1);
      mean += projch->GetBinCenter(k+1)*projch->GetBinContent(k+1);
    }
    if (nentries > 0) {
      fNumberEnt++;
      mean /= nentries;
    }
    //printf("The number of entries for the group %d is %f\n",idect,nentries);
    // Rebin
    if (fRebin >  1) {
      projch = ReBin((TH1F *) projch);
    }
    // This detector has not enough statistics or was not found in VectorCH
    if (nentries <= fMinEntries) {
      NotEnoughStatisticCH(idect);
      continue;
    }
    // Statistic of the histos fitted
    fStatisticMean += nentries;
    fNumberFit++;
    //Method choosen
    switch(fMethod)
      {
      case 0: FitMeanW((TH1 *) projch, nentries); break;
      case 1: FitMean((TH1 *) projch, nentries, mean); break;
      case 2: FitCH((TH1 *) projch, mean); break;
      case 3: FitBisCH((TH1 *) projch, mean); break;
      default: return kFALSE;
      }
    // Fill Infos Fit
    FillInfosFitCH(idect); 
  } // Boucle object
  // Normierungcharge
  if (fDebugLevel != 1) {
    NormierungCharge();
  }
  // Mean Statistics
  if (fNumberFit > 0) {
    AliInfo(Form("There are %d with at least one entries. %d fits have been proceeded (sucessfully or not...). There is a mean statistic of: %d over these fitted histograms and %d successfulled fits",fNumberEnt, fNumberFit, (Int_t) fStatisticMean/fNumberFit, fNumberFitSuccess));
    fStatisticMean = fStatisticMean / fNumberFit;
  }
  else {
    AliInfo(Form("There are %d with at least one entries. There is no fit!",fNumberEnt));
  }
  delete fDebugStreamer;
  fDebugStreamer = 0x0;
  return kTRUE;
}
//____________Functions fit Online CH2d________________________________________
Double_t AliTRDCalibraFit::AnalyseCHAllTogether(const TH2I *ch)
{
  //
  // Fit the 1D histos, projections of the 2D ch on the Xaxis, for each
  // calibration group normalized the resulted coefficients (to 1 normally)
  //
  Int_t    nbins   = ch->GetNbinsX();// charge
  Int_t    nybins  = ch->GetNbinsY();// groups number
  // Take the histo
  TH1I *projch = (TH1I *) ch->ProjectionX("projch",1,nybins+1,(Option_t *)"e");
  projch->SetDirectory(0);
  // Number of entries for this calibration group
  Double_t nentries = 0.0;
  Double_t mean = 0.0;
  for (Int_t k = 0; k < nbins; k++) {
    nentries += projch->GetBinContent(k+1);
    mean += projch->GetBinCenter(k+1)*projch->GetBinContent(k+1);
    projch->SetBinError(k+1,TMath::Sqrt(projch->GetBinContent(k+1)));
  }
  projch->SetEntries(nentries);
  //printf("The number of entries for the group %d is %f\n",idect,nentries);
  if (nentries > 0) {
    mean /= nentries;
  }
  // This detector has not enough statistics or was off
  if (nentries <= fMinEntries) {
    delete projch;
    AliInfo(Form("There are %d entries for all together, it is not enough (%d)",(Int_t) nentries, fMinEntries));
    return -100.0;
  }
  //Method choosen
  switch(fMethod)
    {
    case 0: FitMeanW((TH1 *) projch, nentries); break;
    case 1: FitMean((TH1 *) projch, nentries, mean); break;
    case 2: FitCH((TH1 *) projch, mean); break;
    case 3: FitBisCH((TH1 *) projch, mean); break;
    default: return -100.0;
    }
  delete fDebugStreamer;
  fDebugStreamer = 0x0;
  
  if(fCurrentCoef[0] > 0.0) return fCurrentCoef[0];
  else return -100.0;
  
}
//________________functions fit Online PH2d____________________________________
Double_t AliTRDCalibraFit::AnalysePHAllTogether(const TProfile2D *ph)
{
  //
  // Take the 1D profiles (average pulse height), projections of the 2D PH
  // on the Xaxis, for each calibration group
  // Reconstruct a drift velocity
  // A first calibration of T0 is also made  using the same method
  //

  // Number of Xbins (detectors or groups of pads)
  Int_t    nbins   = ph->GetNbinsX();// time
  Int_t    nybins  = ph->GetNbinsY();// calibration group
 
  // Take the histo
  TH1D *projph = (TH1D *) ph->ProjectionX("projph",1,nybins+1,(Option_t *) "e");
  projph->SetDirectory(0); 
  // Number of entries for this calibration group
  Double_t nentries = 0;
  for(Int_t idect = 0; idect < nybins; idect++){
    for (Int_t k = 0; k < nbins; k++) {
      Int_t binnb = (nbins+2)*(idect+1)+(k+1);
      nentries += ph->GetBinEntries(binnb);
    }
  }
  //printf("AnalysePHAllTogether:: the number of entries is %f\n",nentries);
  // This detector has not enough statistics or was off
  if (nentries  <= fMinEntries) {
    AliInfo(Form("There are %d entries for all together, it is not enough (%d)",(Int_t) nentries, fMinEntries));
    if (fDebugLevel != 1) {
      delete projph;
    }
    return -100.0;
  }
  //Method choosen
  //printf("Method\n");
  switch(fMethod)
    {
    case 0: FitLagrangePoly((TH1 *) projph); break;
    case 1: FitPente((TH1 *) projph); break;
    case 2: FitPH((TH1 *) projph,0); break;
    default: return -100.0;
    }
  // Memory!!!
  if (fDebugLevel != 1) {
    delete projph;
  }
  delete fDebugStreamer;
  fDebugStreamer = 0x0;
  
  if(fCurrentCoef[0] > 0.0) return fCurrentCoef[0];
  else return -100.0;
  
}
//____________Functions fit Online PH2d________________________________________
Bool_t AliTRDCalibraFit::AnalysePH(AliTRDCalibraVector *calvect)
{
  //
  // Reconstruct the average pulse height from the vectorPH for each
  // calibration group
  // Reconstruct a drift velocity
  // A first calibration of T0 is also made  using the same method (slope method)
  //

  // Set the calibration mode
  //const char *name = calvect->GetNamePH();
  TString name = calvect->GetNamePH();
  if(!SetModeCalibration(name,1)) return kFALSE;

  // Number of Xbins (detectors or groups of pads)
  if (!InitFit((432*calvect->GetDetCha0(1)+108*calvect->GetDetCha2(1)),1)) {
    return kFALSE;
  }
  if (!InitFitPH()) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;
  // Init fCountDet and fCount
  InitfCountDetAndfCount(1);
  // Beginning of the loop
  for (Int_t idect = fDect1; idect < fDect2; idect++) {
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi...........
    UpdatefCountDetAndfCount(idect,1);
    ReconstructFitRowMinRowMax(idect,1);
    // Take the histo
    fEntriesCurrent = 0;
    if(!calvect->GetPHEntries(fCountDet)) {
      NotEnoughStatisticPH(idect,fEntriesCurrent);
      continue;
    }
    TString tname("PH");
    tname += idect;
    TH1F *projph  = calvect->CorrectTheError((TGraphErrors *) (calvect->ConvertVectorPHTGraphErrors(fCountDet,(idect-(fCount-(fCalibraMode->GetNfragZ(1)*fCalibraMode->GetNfragRphi(1)))),(const char *) tname)),fEntriesCurrent);
    projph->SetDirectory(0);
    if(fEntriesCurrent > 0) fNumberEnt++;
    //printf("The number of entries for the group %d is %d\n",idect,fEntriesCurrent);
    // This detector has not enough statistics or was off
    if (fEntriesCurrent <=  fMinEntries) {
      //printf("Not enough stat!\n");
      NotEnoughStatisticPH(idect,fEntriesCurrent);
      continue;
    }
    // Statistic of the histos fitted
    fNumberFit++;
    fStatisticMean += fEntriesCurrent;
    // Calcul of "real" coef
    CalculVdriftCoefMean();
    CalculT0CoefMean();
    //Method choosen
    switch(fMethod)
      {
      case 0: FitLagrangePoly((TH1 *) projph); break;
      case 1: FitPente((TH1 *) projph); break;
      case 2: FitPH((TH1 *) projph,(Int_t) (idect - fDect1)); break;
      default: return kFALSE;
      }
    // Fill the tree if end of a detector or only the pointer to the branch!!!
    FillInfosFitPH(idect,fEntriesCurrent);
  } // Boucle object
  
  // Mean Statistic
  if (fNumberFit > 0) {
    AliInfo(Form("There are %d with at least one entries. %d fits have been proceeded (sucessfully or not...). There is a mean statistic of: %d over these fitted histograms and %d successfulled fits",fNumberEnt, fNumberFit, (Int_t) fStatisticMean/fNumberFit,fNumberFitSuccess));
    fStatisticMean = fStatisticMean / fNumberFit;
  }
  else {
    AliInfo(Form("There are %d with at least one entries. There is no fit!",fNumberEnt));
  }
  delete fDebugStreamer;
  fDebugStreamer = 0x0;
  return kTRUE;
}
//________________functions fit Online PH2d____________________________________
Bool_t AliTRDCalibraFit::AnalysePH(const TProfile2D *ph)
{
  //
  // Take the 1D profiles (average pulse height), projections of the 2D PH
  // on the Xaxis, for each calibration group
  // Reconstruct a drift velocity
  // A first calibration of T0 is also made  using the same method
  //

  // Set the calibration mode
  //const char *name = ph->GetTitle();
  TString name = ph->GetTitle();
  if(!SetModeCalibration(name,1)) return kFALSE;
  
  //printf("Mode calibration set\n");

  // Number of Xbins (detectors or groups of pads)
  Int_t    nbins   = ph->GetNbinsX();// time
  Int_t    nybins  = ph->GetNbinsY();// calibration group
  if (!InitFit(nybins,1)) {
    return kFALSE;
  }

  //printf("Init fit\n");

  if (!InitFitPH()) {
    return kFALSE;
  }

  //printf("Init fit PH\n");

  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;
  // Init fCountDet and fCount
  InitfCountDetAndfCount(1);
  //printf("Init Count Det and fCount %d, %d\n",fDect1,fDect2);

  // Beginning of the loop
  for (Int_t idect = fDect1; idect < fDect2; idect++) {
    //printf("idect = %d\n",idect);
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi.......
    UpdatefCountDetAndfCount(idect,1);
    ReconstructFitRowMinRowMax(idect,1);
    // Take the histo
    TH1D *projph = (TH1D *) ph->ProjectionX("projph",idect+1,idect+1,(Option_t *) "e");
    projph->SetDirectory(0); 
    // Number of entries for this calibration group
    Double_t nentries = 0;
    for (Int_t k = 0; k < nbins; k++) {
      Int_t binnb = (nbins+2)*(idect+1)+(k+1);
      nentries += ph->GetBinEntries(binnb);
    }
    if (nentries > 0) {
      fNumberEnt++;
    }  
    //printf("The number of entries for the group %d is %f\n",idect,nentries);
    // This detector has not enough statistics or was off
    if (nentries  <= fMinEntries) {
      //printf("Not enough statistic!\n");
      NotEnoughStatisticPH(idect,nentries);     
      if (fDebugLevel != 1) {
	delete projph;
      }
      continue;
    }
    // Statistics of the histos fitted
    fNumberFit++;
    fStatisticMean += nentries;
    // Calcul of "real" coef
    CalculVdriftCoefMean();
    CalculT0CoefMean();
    //Method choosen
    //printf("Method\n");
    switch(fMethod)
      {
      case 0: FitLagrangePoly((TH1 *) projph); break;
      case 1: FitPente((TH1 *) projph); break;
      case 2: FitPH((TH1 *) projph,(Int_t) (idect - fDect1)); break;
      default: return kFALSE;
      }
    // Fill the tree if end of a detector or only the pointer to the branch!!!
    FillInfosFitPH(idect,nentries);
    // Memory!!!
    if (fDebugLevel != 1) {
      delete projph;
    }
  } // Boucle object
  // Mean Statistic
  if (fNumberFit > 0) {
    AliInfo(Form("There are %d with at least one entries. %d fits have been proceeded (sucessfully or not...). There is a mean statistic of: %d over these fitted histograms and %d successfulled fits",fNumberEnt, fNumberFit, (Int_t) fStatisticMean/fNumberFit,fNumberFitSuccess));
    fStatisticMean = fStatisticMean / fNumberFit;
  }
  else {
    AliInfo(Form("There are %d with at least one entries. There is no fit!",fNumberEnt));
  }
  delete fDebugStreamer;
  fDebugStreamer = 0x0;
  return kTRUE;
}
//____________Functions fit Online PRF2d_______________________________________
Bool_t AliTRDCalibraFit::AnalysePRF(const TProfile2D *prf)
{
  //
  // Take the 1D profiles (pad response function), projections of the 2D PRF
  // on the Xaxis, for each calibration group
  // Fit with a gaussian to reconstruct the sigma of the pad response function
  //

  // Set the calibration mode
  //const char *name = prf->GetTitle();
  TString name = prf->GetTitle();
  if(!SetModeCalibration(name,2)) return kFALSE;

  // Number of Ybins (detectors or groups of pads)
  Int_t    nybins  = prf->GetNbinsY();// calibration groups
  Int_t    nbins   = prf->GetNbinsX();// bins
  Int_t    nbg     = GetNumberOfGroupsPRF((const char *)prf->GetTitle());
  if((nbg > 0) || (nbg == -1)) return kFALSE;
  if (!InitFit(nybins,2)) {
    return kFALSE;
  }
  if (!InitFitPRF()) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;
  // Init fCountDet and fCount
  InitfCountDetAndfCount(2);
  // Beginning of the loop
  for (Int_t idect = fDect1; idect < fDect2; idect++) {
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi......
    UpdatefCountDetAndfCount(idect,2);
    ReconstructFitRowMinRowMax(idect,2);
    // Take the histo
    TH1D *projprf = (TH1D *) prf->ProjectionX("projprf",idect+1,idect+1,(Option_t *) "e");
    projprf->SetDirectory(0);
    // Number of entries for this calibration group
    Double_t nentries = 0;
    for (Int_t k = 0; k < nbins; k++) {
      Int_t binnb = (nbins+2)*(idect+1)+(k+1);
      nentries += prf->GetBinEntries(binnb);
    }
    if(nentries > 0) fNumberEnt++;
    // This detector has not enough statistics or was off
    if (nentries <= fMinEntries) {
      NotEnoughStatisticPRF(idect);
      if (fDebugLevel != 1) {
	delete projprf;
      }
      continue;
    }
    // Statistics of the histos fitted
    fNumberFit++;
    fStatisticMean += nentries;
    // Calcul of "real" coef
    CalculPRFCoefMean();
    //Method choosen
    switch(fMethod)
      {
      case 0: FitPRF((TH1 *) projprf); break;
      case 1: RmsPRF((TH1 *) projprf); break;
      default: return kFALSE;
      }
    // Fill the tree if end of a detector or only the pointer to the branch!!!
    FillInfosFitPRF(idect);
    // Memory!!!
    if (fDebugLevel != 1) {
      delete projprf;
    }
  } // Boucle object
  // Mean Statistic
  if (fNumberFit > 0) {
    AliInfo(Form("There are %d with at least one entries.",fNumberEnt));
    AliInfo(Form("%d fits have been proceeded (sucessfully or not...).",fNumberFit));
    AliInfo(Form("There is a mean statistic of: %d over these fitted histograms and %d successfulled fits"
                ,(Int_t) fStatisticMean/fNumberFit,fNumberFitSuccess));
    fStatisticMean = fStatisticMean / fNumberFit;
  }
  else {
    AliInfo(Form("There are %d with at least one entries. There is no fit!",fNumberEnt));
  }
  delete fDebugStreamer;
  fDebugStreamer = 0x0;
  return kTRUE;
}
//____________Functions fit Online PRF2d_______________________________________
Bool_t AliTRDCalibraFit::AnalysePRFMarianFit(const TProfile2D *prf)
{
  //
  // Take the 1D profiles (pad response function), projections of the 2D PRF
  // on the Xaxis, for each calibration group
  // Fit with a gaussian to reconstruct the sigma of the pad response function
  //

  // Set the calibration mode
  //const char *name = prf->GetTitle();
  TString name = prf->GetTitle();
  if(!SetModeCalibration(name,2)) return kFALSE;

  // Number of Ybins (detectors or groups of pads)
  TAxis   *xprf    = prf->GetXaxis();
  TAxis   *yprf    = prf->GetYaxis();
  Int_t    nybins  = yprf->GetNbins();// calibration groups
  Int_t    nbins   = xprf->GetNbins();// bins
  Float_t  lowedge = (Float_t) xprf->GetBinLowEdge(1);//lowedge in bins
  Float_t  upedge  = (Float_t) xprf->GetBinUpEdge(nbins);//upedge in bins
  Int_t    nbg     = GetNumberOfGroupsPRF((const char *)name);
  if(nbg == -1) return kFALSE;
  if(nbg > 0) fMethod = 1;
  else fMethod = 0;
  if (!InitFit(nybins,2)) {
    return kFALSE;
  }
  if (!InitFitPRF()) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;
  // Init fCountDet and fCount
  InitfCountDetAndfCount(2);
  // Beginning of the loop
  for (Int_t idect = fDect1; idect < fDect2; idect++) {
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi.......
    UpdatefCountDetAndfCount(idect,2);
    ReconstructFitRowMinRowMax(idect,2);
    // Build the array of entries and sum
    TArrayD arraye  = TArrayD(nbins);
    TArrayD arraym  = TArrayD(nbins);
    TArrayD arrayme = TArrayD(nbins);
    Double_t nentries = 0;
    //printf("nbins %d\n",nbins);
    for (Int_t k = 0; k < nbins; k++) {
      Int_t binnb = (nbins+2)*(idect+1)+(k+1);
      Double_t entries = (Double_t)prf->GetBinEntries(binnb);
      Double_t mean    = (Double_t)prf->GetBinContent(binnb);
      Double_t error   = (Double_t)prf->GetBinError(binnb); 
      //printf("for %d we have %f\n",k,entries);
      nentries += entries;
      arraye.AddAt(entries,k);
      arraym.AddAt(mean,k);
      arrayme.AddAt(error,k);
    }
    if(nentries > 0) fNumberEnt++;
    //printf("The number of entries for the group %d is %f\n",idect,nentries);
    // This detector has not enough statistics or was off
    if (nentries <= fMinEntries) {
      NotEnoughStatisticPRF(idect);
      continue;
    }
    // Statistics of the histos fitted
    fNumberFit++;
    fStatisticMean += nentries;
    // Calcul of "real" coef
    CalculPRFCoefMean();
    //Method choosen
    switch(fMethod)
      {
      case 0: FitPRFGausMI( arraye.GetArray(), arraym.GetArray(), arrayme.GetArray(), nbins, lowedge, upedge); break;
      case 1: FitTnpRange( arraye.GetArray(), arraym.GetArray(), arrayme.GetArray(), nbg, nbins); break;
      default: return kFALSE;
      }
    // Fill the tree if end of a detector or only the pointer to the branch!!!
    FillInfosFitPRF(idect);
  } // Boucle object
  // Mean Statistic
  if (fNumberFit > 0) {
    AliInfo(Form("There are %d with at least one entries.",fNumberEnt));
    AliInfo(Form("%d fits have been proceeded (sucessfully or not...).",fNumberFit));
    AliInfo(Form("There is a mean statistic of: %d over these fitted histograms and %d successfulled fits"
                ,(Int_t) fStatisticMean/fNumberFit,fNumberFitSuccess));
    fStatisticMean = fStatisticMean / fNumberFit;
  }
  else {
    AliInfo(Form("There are %d with at least one entries. There is no fit!",fNumberEnt));
  }
  delete fDebugStreamer;
  fDebugStreamer = 0x0;
  return kTRUE;
}
//____________Functions fit Online PRF2d_______________________________________
Bool_t AliTRDCalibraFit::AnalysePRF(AliTRDCalibraVector *calvect)
{
  //
  // Reconstruct the 1D histo (pad response function) from the vectorPRD for
  // each calibration group
  // Fit with a gaussian to reconstruct the sigma of the pad response function
  //

  // Set the calibra mode
  //const char *name = calvect->GetNamePRF();
  TString name = calvect->GetNamePRF();
  if(!SetModeCalibration(name,2)) return kFALSE;
  //printf("test0 %s\n",name);

  // Number of Xbins (detectors or groups of pads)
  if (!InitFit((432*calvect->GetDetCha0(2)+108*calvect->GetDetCha2(2)),2)) {
    //printf("test1\n");
    return kFALSE;
  }
  if (!InitFitPRF()) {
    ///printf("test2\n");
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;
  // Init fCountDet and fCount
  InitfCountDetAndfCount(2);
  // Beginning of the loop
  for (Int_t idect = fDect1; idect < fDect2; idect++) {
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi........
    UpdatefCountDetAndfCount(idect,2);
    ReconstructFitRowMinRowMax(idect,2);
    // Take the histo
    fEntriesCurrent = 0;
    if(!calvect->GetPRFEntries(fCountDet)) {
      NotEnoughStatisticPRF(idect);
      continue;
    }
    TString tname("PRF");
    tname += idect;
    TH1F *projprf  = calvect->CorrectTheError((TGraphErrors *) (calvect->ConvertVectorPRFTGraphErrors(fCountDet,(idect-(fCount-(fCalibraMode->GetNfragZ(1)*fCalibraMode->GetNfragRphi(1)))),(const char *) tname)),fEntriesCurrent);
    projprf->SetDirectory(0);
    if(fEntriesCurrent > 0) fNumberEnt++;
    // This detector has not enough statistics or was off
    if (fEntriesCurrent <= fMinEntries) {
      NotEnoughStatisticPRF(idect);
      continue;
    }
    // Statistic of the histos fitted
    fNumberFit++;
    fStatisticMean += fEntriesCurrent;
    // Calcul of "real" coef
    CalculPRFCoefMean();
    //Method choosen
    switch(fMethod)
      {
      case 1: FitPRF((TH1 *) projprf); break;
      case 2: RmsPRF((TH1 *) projprf); break;
      default: return kFALSE;
      }    
    // Fill the tree if end of a detector or only the pointer to the branch!!!
    FillInfosFitPRF(idect);
  } // Boucle object
  // Mean Statistics
  if (fNumberFit > 0) {
    AliInfo(Form("There are %d with at least one entries. %d fits have been proceeded (sucessfully or not...). There is a mean statistic of: %d over these fitted histograms and %d successfulled fits",fNumberEnt, fNumberFit, (Int_t) fStatisticMean/fNumberFit,fNumberFitSuccess));
  }
  else {
    AliInfo(Form("There are %d with at least one entries. There is no fit!",fNumberEnt));
  }
  delete fDebugStreamer;
  fDebugStreamer = 0x0;
  return kTRUE;
}
//____________Functions fit Online PRF2d_______________________________________
Bool_t AliTRDCalibraFit::AnalysePRFMarianFit(AliTRDCalibraVector *calvect)
{
  //
  // Reconstruct the 1D histo (pad response function) from the vectorPRD for
  // each calibration group
  // Fit with a gaussian to reconstruct the sigma of the pad response function
  //

  // Set the calibra mode
  //const char *name = calvect->GetNamePRF();
  TString name = calvect->GetNamePRF();
  if(!SetModeCalibration(name,2)) return kFALSE;
  //printf("test0 %s\n",name);
  Int_t    nbg     = GetNumberOfGroupsPRF((const char *)name);
  //printf("test1 %d\n",nbg);
  if(nbg == -1) return kFALSE;
  if(nbg > 0) fMethod = 1;
  else fMethod = 0;
  // Number of Xbins (detectors or groups of pads)
  if (!InitFit((432*calvect->GetDetCha0(2)+108*calvect->GetDetCha2(2)),2)) {
    //printf("test2\n");
    return kFALSE;
  }
  if (!InitFitPRF()) {
    //printf("test3\n");
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;
  // Variables
  Int_t nbins           = 0;
  Double_t *arrayx       = 0;
  Double_t *arraye       = 0;
  Double_t *arraym       = 0;
  Double_t *arrayme      = 0;
  Float_t lowedge       = 0.0;
  Float_t upedge        = 0.0;
  // Init fCountDet and fCount
  InitfCountDetAndfCount(2);
  // Beginning of the loop
  for (Int_t idect = fDect1; idect < fDect2; idect++) {
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi......
    UpdatefCountDetAndfCount(idect,2);
    ReconstructFitRowMinRowMax(idect,2);
    // Take the histo
    fEntriesCurrent  = 0;
    if(!calvect->GetPRFEntries(fCountDet)) {
      NotEnoughStatisticPRF(idect);
      continue;
    }
    TString tname("PRF");
    tname += idect;
    TGraphErrors *projprftree  = calvect->ConvertVectorPRFTGraphErrors(fCountDet,(idect-(fCount-(fCalibraMode->GetNfragZ(1)*fCalibraMode->GetNfragRphi(1)))),(const char *) tname);
    nbins   = projprftree->GetN();
    arrayx  = (Double_t *)projprftree->GetX();
    arraye  = (Double_t *)projprftree->GetEX();
    arraym  = (Double_t *)projprftree->GetY();
    arrayme = (Double_t *)projprftree->GetEY();
    Float_t step = arrayx[1]-arrayx[0];
    lowedge = arrayx[0] - step/2.0;
    upedge  = arrayx[(nbins-1)] + step/2.0;
    //printf("nbins est %d\n",nbins);
    for(Int_t k = 0; k < nbins; k++){
      fEntriesCurrent += (Int_t)arraye[k];
      //printf("for %d we have %f, %f\n",k,arraye[k],((projprftree->GetEX())[k]));
      if(arraye[k]>0.0) arrayme[k] = TMath::Sqrt(TMath::Abs(arrayme[k]-arraym[k]*arraym[k])/arraye[k]);
    }
    if(fEntriesCurrent > 0) fNumberEnt++;
    //printf("The number of entries for the group %d is %d\n",idect,fEntriesCurrent);
    // This detector has not enough statistics or was off
    if (fEntriesCurrent <= fMinEntries) {
      NotEnoughStatisticPRF(idect);
      continue;
    }
    // Statistic of the histos fitted
    fNumberFit++;
    fStatisticMean += fEntriesCurrent;
    // Calcul of "real" coef
    CalculPRFCoefMean();
    //Method choosen
    switch(fMethod)
      {
      case 0: FitPRFGausMI(arraye,arraym,arrayme,nbins,lowedge,upedge); break;
      case 1: FitTnpRange(arraye,arraym,arrayme,nbg,nbins); break;
      default: return kFALSE;
      }    
    // Fill the tree if end of a detector or only the pointer to the branch!!!
    FillInfosFitPRF(idect);
  } // Boucle object
  // Mean Statistics
  if (fNumberFit > 0) {
    AliInfo(Form("There are %d with at least one entries. %d fits have been proceeded (sucessfully or not...). There is a mean statistic of: %d over these fitted histograms and %d successfulled fits",fNumberEnt, fNumberFit, (Int_t) fStatisticMean/fNumberFit,fNumberFitSuccess));
  }
  else {
    AliInfo(Form("There are %d with at least one entries. There is no fit!",fNumberEnt));
  }
  delete fDebugStreamer;
  fDebugStreamer = 0x0;
  return kTRUE;
}
//____________Functions fit Online CH2d________________________________________
Bool_t AliTRDCalibraFit::AnalyseLinearFitters(AliTRDCalibraVdriftLinearFit *calivdli)
{
  //
  // The linear method
  //

  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;
  if(!InitFitLinearFitter()) return kFALSE;

  
  for(Int_t idet = 0; idet < 540; idet++){


    //printf("detector number %d\n",idet);

    // Take the result
    TVectorD param(2);
    TVectorD error(3);
    fEntriesCurrent = 0;
    fCountDet       = idet;
    Bool_t here     = calivdli->GetParam(idet,&param);
    Bool_t heree    = calivdli->GetError(idet,&error);
    //printf("here %d and heree %d\n",here, heree);
    if(heree) {
      fEntriesCurrent = (Int_t) error[2];
      fNumberEnt++;
    }
    //printf("Number of entries %d\n",fEntriesCurrent);
    // Nothing found or not enough statistic
    if((!heree) || (!here) || (fEntriesCurrent <= fMinEntries)) {
      NotEnoughStatisticLinearFitter();
      continue;
    }
    //param.Print();
    //error.Print();
    //Statistics
    fNumberFit++;
    fStatisticMean += fEntriesCurrent;     

    // Check the fit
    if((-(param[1])) <= 0.0) {
      NotEnoughStatisticLinearFitter();
      continue;
    }

    // CalculDatabaseVdriftandTan
    CalculVdriftLorentzCoef();
    //printf("AliTRDCalibraFit::AnalyzeVdriftLinearFit detector %d, vdrift %f and %f and exB %f and %f\n",idet,fCalDetVdriftUsed->GetValue(idet),fCurrentCoef[1],fCalDetExBUsed->GetValue(idet),fCurrentCoef2[1]);

    // Statistics   
    fNumberFitSuccess ++;

    // Put the fCurrentCoef
    fCurrentCoef[0]  = -param[1];
    // here the database must be the one of the reconstruction for the lorentz angle....
    fCurrentCoef2[0] = (param[0]+fCurrentCoef[1]*fCurrentCoef2[1])/fCurrentCoef[0];
    fCurrentCoefE    = error[1];
    fCurrentCoefE2   = error[0];
    if((TMath::Abs(fCurrentCoef2[0]) > 0.0000001) && (TMath::Abs(param[0]) > 0.0000001)){
      fCurrentCoefE2 = (fCurrentCoefE2/param[0]+fCurrentCoefE/fCurrentCoef[0])*fCurrentCoef2[0];
    }    

    // Fill
    FillInfosFitLinearFitter();

    
  }
  // Mean Statistics
  if (fNumberFit > 0) {
    AliInfo(Form("There are %d with at least one entries. %d fits have been proceeded (sucessfully or not...). There is a mean statistic of: %d over these fitted histograms and %d successfulled fits",fNumberEnt, fNumberFit, (Int_t) fStatisticMean/fNumberFit,fNumberFitSuccess));
  }
  else {
    AliInfo(Form("There are %d with at least one entries. There is no fit!",fNumberEnt));
  }
  delete fDebugStreamer;
  fDebugStreamer = 0x0;
  return kTRUE;
  
}
//____________Functions fit Online CH2d________________________________________
void AliTRDCalibraFit::AnalyseLinearFittersAllTogether(AliTRDCalibraVdriftLinearFit *calivdli, Double_t &vdriftoverall, Double_t &exboverall)
{
  //
  // The linear method
  //

  // Get the mean vdrift and exb used
  Double_t meanvdriftused = 0.0;
  Double_t meanexbused = 0.0;
  Double_t counterdet = 0.0;
  if((!fCalDetVdriftUsed) || (!fCalDetExBUsed)) {
    vdriftoverall = -100.0;
    exboverall = 100.0;
    return;
  }  

  // Add histos

  TH2S *linearfitterhisto = 0x0;
  
  for(Int_t idet = 0; idet < 540; idet++){
    
    TH2S * u = calivdli->GetLinearFitterHistoForce(idet);
    Double_t detectorentries = u->Integral();
    meanvdriftused += fCalDetVdriftUsed->GetValue(idet)*detectorentries;
    meanexbused += fCalDetExBUsed->GetValue(idet)*detectorentries;
    counterdet += detectorentries;

    //printf("detectorentries %f\n",detectorentries);
    
    //printf("AliTRDCalibraFit::AnalyzeVdriftLinearFitsAllTogether detector %d, vdrift %f and exB %f\n",idet,fCalDetVdriftUsed->GetValue(idet),fCalDetExBUsed->GetValue(idet));

    if(idet == 0) linearfitterhisto = u;
    else linearfitterhisto->Add(u);

  }
  if(counterdet > 0.0){
    meanvdriftused = meanvdriftused/counterdet;
    meanexbused = meanexbused/counterdet;    
  }
  else {
    vdriftoverall = -100.0;
    exboverall = 100.0;
    return;
  }
  
  
  //printf("AliTRDCalibraFit::AnalyzeVdriftLinearFitsAllTogether MEAN vdrift %f and exB %f\n",meanvdriftused,meanexbused);

  // Fit

  Int_t entries = 0;
  TAxis *xaxis = linearfitterhisto->GetXaxis();
  TAxis *yaxis = linearfitterhisto->GetYaxis();
  TLinearFitter linearfitter = TLinearFitter(2,"pol1");
  //printf("test\n");
  Double_t integral = linearfitterhisto->Integral();
  //printf("Integral is %f\n",integral);
  Bool_t securitybreaking = kFALSE;
  if(TMath::Abs(integral-1199) < 0.00001) securitybreaking = kTRUE;
  for(Int_t ibinx = 0; ibinx < linearfitterhisto->GetNbinsX(); ibinx++){
    for(Int_t ibiny = 0; ibiny < linearfitterhisto->GetNbinsY(); ibiny++){
      if(linearfitterhisto->GetBinContent(ibinx+1,ibiny+1)>0){
	Double_t x = xaxis->GetBinCenter(ibinx+1);
	Double_t y = yaxis->GetBinCenter(ibiny+1);
	
	for(Int_t k = 0; k < (Int_t)linearfitterhisto->GetBinContent(ibinx+1,ibiny+1); k++){
	  if(!securitybreaking){
	    linearfitter.AddPoint(&x,y);
	    entries++;
	  }
	  else {
	    if(entries< 1198){
	      linearfitter.AddPoint(&x,y);
	      entries++; 
	    }
	  }
	}
	
      }
    }
  }
      
  //printf("AnalyseLinearFittersAllTogether::Find %d entries\n",entries);
  //printf("Minstats %d\n",fMinEntries);

  

  // Eval the linear fitter
  if(entries > fMinEntries){
    TVectorD  par  = TVectorD(2);
    //printf("Fit\n");
    if((linearfitter.EvalRobust(0.8)==0)) {
      //printf("Take the param\n");
      linearfitter.GetParameters(par);
      //printf("Done\n");
      //par.Print();
      //printf("Finish\n");
      // Put the fCurrentCoef
      fCurrentCoef[0]  = -par[1];
      // here the database must be the one of the reconstruction for the lorentz angle....
      if(fCurrentCoef[0] > 0.0) fCurrentCoef2[0] = (par[0]+meanvdriftused*meanexbused)/fCurrentCoef[0];
      else fCurrentCoef2[0] = 100.0;      

    }
    else {
      
      fCurrentCoef[0] = -100.0;
      fCurrentCoef2[0] = 100.0;
      
    }
    
    
  }
  else {

    fCurrentCoef[0] = -100.0;
    fCurrentCoef2[0] = 100.0;
    
  }
  
  vdriftoverall = fCurrentCoef[0];
  exboverall = fCurrentCoef2[0];
  

  delete linearfitterhisto;
  delete fDebugStreamer;
  fDebugStreamer = 0x0;
  
}
//____________Functions for seeing if the pad is really okey___________________
//_____________________________________________________________________________
Int_t AliTRDCalibraFit::GetNumberOfGroupsPRF(TString nametitle)
{
  //
  // Get numberofgroupsprf
  //
  
  // Some patterns
  const Char_t *pattern0 = "Ngp0";
  const Char_t *pattern1 = "Ngp1";
  const Char_t *pattern2 = "Ngp2";
  const Char_t *pattern3 = "Ngp3";
  const Char_t *pattern4 = "Ngp4";
  const Char_t *pattern5 = "Ngp5";
  const Char_t *pattern6 = "Ngp6";

  // Nrphi mode
  if (strstr(nametitle.Data(),pattern0)) {
    return 0;
  }
  if (strstr(nametitle.Data(),pattern1)) {
    return 1;
  }
  if (strstr(nametitle.Data(),pattern2)) {
    return 2;
  }
  if (strstr(nametitle.Data(),pattern3)) {
    return 3;
  }
  if (strstr(nametitle.Data(),pattern4)) {
    return 4;
  }
  if (strstr(nametitle.Data(),pattern5)) {
    return 5;
  }
  if (strstr(nametitle.Data(),pattern6)){
    return 6;
  }
  else return -1;
 

}
//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::SetModeCalibration(TString name, Int_t i)
{
  //
  // Set fNz[i] and fNrphi[i] of the AliTRDCalibraFit::Instance()
  // corresponding to the given name
  //

  if(!SetNzFromTObject(name,i)) return kFALSE;
  if(!SetNrphiFromTObject(name,i)) return kFALSE;
  
  return kTRUE; 

}
//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::SetNrphiFromTObject(TString name, Int_t i)
{
  //
  // Set fNrphi[i] of the AliTRDCalibraFit::Instance()
  // corresponding to the given TObject
  //
  
  // Some patterns
  const Char_t *patternrphi0 = "Nrphi0";
  const Char_t *patternrphi1 = "Nrphi1";
  const Char_t *patternrphi2 = "Nrphi2";
  const Char_t *patternrphi3 = "Nrphi3";
  const Char_t *patternrphi4 = "Nrphi4";
  const Char_t *patternrphi5 = "Nrphi5";
  const Char_t *patternrphi6 = "Nrphi6";

  
  const Char_t *patternrphi10 = "Nrphi10";
  const Char_t *patternrphi100 = "Nrphi100";
  const Char_t *patternz10 = "Nz10";
  const Char_t *patternz100 = "Nz100";

  // Nrphi mode
  if ((strstr(name.Data(),patternrphi100)) && (strstr(name.Data(),patternz100))) {
    fCalibraMode->SetAllTogether(i);
    fNbDet = 540;
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 100",fNbDet));
    }
    return kTRUE;
  }
  if ((strstr(name.Data(),patternrphi10)) && (strstr(name.Data(),patternz10))) {
    fCalibraMode->SetPerSuperModule(i);
    fNbDet = 30;
    if (fDebugLevel > 1) {
      AliInfo(Form("fNDet %d and 100",fNbDet));
    }
    return kTRUE;
  }
  
  if (strstr(name.Data(),patternrphi0)) {
    fCalibraMode->SetNrphi(i ,0);
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 0",fNbDet));
    }
    return kTRUE;
  }
  if (strstr(name.Data(),patternrphi1)) {
    fCalibraMode->SetNrphi(i, 1);
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 1",fNbDet));
    }
    return kTRUE;
  }
  if (strstr(name.Data(),patternrphi2)) {
    fCalibraMode->SetNrphi(i, 2);
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 2",fNbDet));
    }    
    return kTRUE;
  }
  if (strstr(name.Data(),patternrphi3)) {
    fCalibraMode->SetNrphi(i, 3);
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 3",fNbDet));
    }   
    return kTRUE;
  }
  if (strstr(name.Data(),patternrphi4)) {
    fCalibraMode->SetNrphi(i, 4);
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 4",fNbDet));
    }   
    return kTRUE;
  }
  if (strstr(name.Data(),patternrphi5)) {
    fCalibraMode->SetNrphi(i, 5);
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 5",fNbDet));
    }
    return kTRUE;
  }
  if (strstr(name.Data(),patternrphi6)) {
    fCalibraMode->SetNrphi(i, 6);
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 6",fNbDet));
    }
    return kTRUE;
  }
  
  if (fDebugLevel > 1) {
    AliInfo(Form("fNbDet %d and rest",fNbDet));
  }
  fCalibraMode->SetNrphi(i ,0);
  return kFALSE;
  
}
//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::SetNzFromTObject(TString name, Int_t i)
{
  //
  // Set fNz[i] of the AliTRDCalibraFit::Instance()
  // corresponding to the given TObject
  //

  // Some patterns
  const Char_t *patternz0    = "Nz0";
  const Char_t *patternz1    = "Nz1";
  const Char_t *patternz2    = "Nz2";
  const Char_t *patternz3    = "Nz3";
  const Char_t *patternz4    = "Nz4";

  const Char_t *patternrphi10 = "Nrphi10";
  const Char_t *patternrphi100 = "Nrphi100";
  const Char_t *patternz10 = "Nz10";
  const Char_t *patternz100 = "Nz100";

  if ((strstr(name.Data(),patternrphi100)) && (strstr(name.Data(),patternz100))) {
    fCalibraMode->SetAllTogether(i);
    fNbDet = 540;
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 100",fNbDet));
    }
    return kTRUE;
  }
  if ((strstr(name.Data(),patternrphi10)) && (strstr(name.Data(),patternz10))) {
    fCalibraMode->SetPerSuperModule(i);
    fNbDet = 30;
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 10",fNbDet));
    }
    return kTRUE;
  }
  if (strstr(name.Data(),patternz0)) {
    fCalibraMode->SetNz(i, 0);
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 0",fNbDet));
    }
    return kTRUE;
  }
  if (strstr(name.Data(),patternz1)) {
    fCalibraMode->SetNz(i ,1);
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 1",fNbDet));
    }
    return kTRUE;
  }
  if (strstr(name.Data(),patternz2)) {
    fCalibraMode->SetNz(i ,2);
    if (fDebugLevel > 1) {    
      AliInfo(Form("fNbDet %d and 2",fNbDet));
    }
    return kTRUE;
  }
  if (strstr(name.Data(),patternz3)) {
    fCalibraMode->SetNz(i ,3);
    if (fDebugLevel > 1) {
      AliInfo(Form("fNbDet %d and 3",fNbDet));
    }
    return kTRUE;  
  }
  if (strstr(name.Data(),patternz4)) {
    fCalibraMode->SetNz(i ,4);
    if (fDebugLevel > 1) {    
      AliInfo(Form("fNbDet %d and 4",fNbDet));
    }
    return kTRUE;
  }
 
  if (fDebugLevel > 1) {
    AliInfo(Form("fNbDet %d and rest",fNbDet));
  }
  fCalibraMode->SetNz(i ,0);
  return kFALSE;
}
//______________________________________________________________________
void AliTRDCalibraFit::RemoveOutliers(Int_t type, Bool_t perdetector){
  //
  // Remove the results too far from the mean value and rms
  // type: 0 gain, 1 vdrift
  // perdetector
  //

  Int_t loop = (Int_t) fVectorFit.GetEntriesFast();
  if(loop != 540) {
    AliInfo("The Vector Fit is not complete!");
    return;
  }
  Int_t detector = -1;
  Int_t sector = -1;
  Float_t value  = 0.0;

  /////////////////////////////////
  // Calculate the mean values
  ////////////////////////////////
  // Initialisation
  ////////////////////////
  Double_t meanAll = 0.0;
   Double_t rmsAll = 0.0;
   Int_t countAll = 0;
   ////////////
  // compute
  ////////////
  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) fVectorFit.At(k))->GetDetector();
    sector = GetSector(detector);
    if(perdetector){
      value = ((AliTRDFitInfo *) fVectorFit.At(k))->GetCoef()[0];
      if(value > 0.0) {
	rmsAll += value*value;
	meanAll += value;
	countAll++;
      }
    }
    else {
      Int_t rowMax    = fGeo->GetRowMax(GetLayer(detector),GetStack(detector),GetSector(detector));
      Int_t colMax    = fGeo->GetColMax(GetLayer(detector));
      for (Int_t row = 0; row < rowMax; row++) {
	for (Int_t col = 0; col < colMax; col++) {
	  value = ((AliTRDFitInfo *) fVectorFit.At(k))->GetCoef()[(Int_t)(col*rowMax+row)];
	  if(value > 0.0) {
	    rmsAll += value*value;
	    meanAll += value;
	    countAll++;
	  }
	  
	} // Col
      } // Row
    }
  }  
  if(countAll > 0) {
    meanAll = meanAll/countAll;
    rmsAll = TMath::Sqrt(TMath::Abs(rmsAll/countAll - (meanAll*meanAll)));
  }
  //printf("RemoveOutliers: meanAll %f and rmsAll %f\n",meanAll,rmsAll);
  /////////////////////////////////////////////////
  // Remove outliers
  ////////////////////////////////////////////////
  Double_t defaultvalue = -1.0;
  if(type==1) defaultvalue = -1.5;
  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) fVectorFit.At(k))->GetDetector();
    sector = GetSector(detector);
    Int_t rowMax    = fGeo->GetRowMax(GetLayer(detector),GetStack(detector),GetSector(detector));
    Int_t colMax    = fGeo->GetColMax(GetLayer(detector));
    Float_t *coef = ((AliTRDFitInfo *) fVectorFit.At(k))->GetCoef();
  
    // remove the results too far away  
    for (Int_t row = 0; row < rowMax; row++) {
      for (Int_t col = 0; col < colMax; col++) {
	value = coef[(Int_t)(col*rowMax+row)];
	if((value > 0.0) && (rmsAll > 0.0) && (TMath::Abs(value-meanAll) > (2*rmsAll))) {
	  coef[(Int_t)(col*rowMax+row)] = defaultvalue;
	}
      } // Col
    } // Row
  }
}
//______________________________________________________________________
void AliTRDCalibraFit::RemoveOutliers2(Bool_t perdetector){
  //
  // Remove the results too far from the mean and rms
  // perdetector
  //

  Int_t loop = (Int_t) fVectorFit2.GetEntriesFast();
  if(loop != 540) {
    AliInfo("The Vector Fit is not complete!");
    return;
  }
  Int_t detector = -1;
  Int_t sector = -1;
  Float_t value  = 0.0;

  /////////////////////////////////
  // Calculate the mean values
  ////////////////////////////////
  // Initialisation
  ////////////////////////
  Double_t meanAll = 0.0;
  Double_t rmsAll = 0.0;
  Int_t countAll = 0;
  /////////////
  // compute
  ////////////
  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) fVectorFit2.At(k))->GetDetector();
    sector = GetSector(detector);
    if(perdetector){
      value = ((AliTRDFitInfo *) fVectorFit2.At(k))->GetCoef()[0];
      if(value < 70.0) {
	meanAll += value;
	rmsAll += value*value;
	countAll++;
      }
    }
    else {
      Int_t rowMax    = fGeo->GetRowMax(GetLayer(detector),GetStack(detector),GetSector(detector));
      Int_t colMax    = fGeo->GetColMax(GetLayer(detector));
      for (Int_t row = 0; row < rowMax; row++) {
	for (Int_t col = 0; col < colMax; col++) {
	  value = ((AliTRDFitInfo *) fVectorFit2.At(k))->GetCoef()[(Int_t)(col*rowMax+row)];
	  if(value < 70.0) {
	    rmsAll += value*value;
	    meanAll += value;
	    countAll++;
	  }	  
	} // Col
      } // Row
    }
  }  
  if(countAll > 0) {
    meanAll = meanAll/countAll;
    rmsAll = TMath::Sqrt(TMath::Abs(rmsAll/countAll - (meanAll*meanAll)));
  }
  //printf("Remove outliers 2: meanAll %f, rmsAll %f\n",meanAll,rmsAll);
  /////////////////////////////////////////////////
  // Remove outliers
  ////////////////////////////////////////////////
  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) fVectorFit2.At(k))->GetDetector();
    sector = GetSector(detector);
    Int_t rowMax    = fGeo->GetRowMax(GetLayer(detector),GetStack(detector),GetSector(detector));
    Int_t colMax    = fGeo->GetColMax(GetLayer(detector));
    Float_t *coef = ((AliTRDFitInfo *) fVectorFit2.At(k))->GetCoef();
  
    // remove the results too far away  
    for (Int_t row = 0; row < rowMax; row++) {
      for (Int_t col = 0; col < colMax; col++) {
	value = coef[(Int_t)(col*rowMax+row)];
	if((value < 70.0) && (rmsAll > 0.0) && (TMath::Abs(value-meanAll) > (2.5*rmsAll))) {
	  //printf("value outlier %f\n",value);
	  coef[(Int_t)(col*rowMax+row)] = 100.0;
	}
      } // Col
    } // Row
  }
}
//______________________________________________________________________
void AliTRDCalibraFit::PutMeanValueOtherVectorFit(Int_t ofwhat, Bool_t perdetector){
  //
  // ofwhat is equaled to 0: mean value of all passing detectors
  // ofwhat is equaled to 1: mean value of the detector, otherwise supermodule, otherwise all
  //

  Int_t loop = (Int_t) fVectorFit.GetEntriesFast();
  if(loop != 540) {
    AliInfo("The Vector Fit is not complete!");
    return;
  }
  Int_t detector = -1;
  Int_t sector = -1;
  Float_t value  = 0.0;

  /////////////////////////////////
  // Calculate the mean values
  ////////////////////////////////
  // Initialisation
  ////////////////////////
  Double_t meanAll = 0.0;
  Double_t meanSupermodule[18];
  Double_t meanDetector[540];
  Double_t rmsAll = 0.0;
  Double_t rmsSupermodule[18];
  Double_t rmsDetector[540];
  Int_t countAll = 0;
  Int_t countSupermodule[18];
  Int_t countDetector[540];
  for(Int_t sm = 0; sm < 18; sm++){
    rmsSupermodule[sm] = 0.0;
    meanSupermodule[sm] = 0.0;
    countSupermodule[sm] = 0;
  }
  for(Int_t det = 0; det < 540; det++){
    rmsDetector[det] = 0.0;
    meanDetector[det] = 0.0;
    countDetector[det] = 0;
  }
  ////////////
  // compute
  ////////////
  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) fVectorFit.At(k))->GetDetector();
    sector = GetSector(detector);
    if(perdetector){
      value = ((AliTRDFitInfo *) fVectorFit.At(k))->GetCoef()[0];
      if(value > 0.0) {
	rmsDetector[detector] += value*value;
	meanDetector[detector] += value;
	countDetector[detector]++;
	rmsSupermodule[sector] += value*value;
	meanSupermodule[sector] += value;
	countSupermodule[sector]++;
	rmsAll += value*value;
	meanAll += value;
	countAll++;
      }
    }
    else {
      Int_t rowMax    = fGeo->GetRowMax(GetLayer(detector),GetStack(detector),GetSector(detector));
      Int_t colMax    = fGeo->GetColMax(GetLayer(detector));
      for (Int_t row = 0; row < rowMax; row++) {
	for (Int_t col = 0; col < colMax; col++) {
	  value = ((AliTRDFitInfo *) fVectorFit.At(k))->GetCoef()[(Int_t)(col*rowMax+row)];
	  if(value > 0.0) {
	    rmsDetector[detector] += value*value;
	    meanDetector[detector] += value;
	    countDetector[detector]++;
	    rmsSupermodule[sector] += value*value;
	    meanSupermodule[sector] += value;
	    countSupermodule[sector]++;
	    rmsAll += value*value;
	    meanAll += value;
	    countAll++;
	  }
	  
	} // Col
      } // Row
    }
  }  
  if(countAll > 0) {
    meanAll = meanAll/countAll;
    rmsAll = TMath::Abs(rmsAll/countAll - (meanAll*meanAll));
  }
  for(Int_t sm = 0; sm < 18; sm++){
    if(countSupermodule[sm] > 0) {
      meanSupermodule[sm] = meanSupermodule[sm]/countSupermodule[sm];
      rmsSupermodule[sm] = TMath::Abs(rmsSupermodule[sm]/countSupermodule[sm] - (meanSupermodule[sm]*meanSupermodule[sm]));
    }
  }
  for(Int_t det = 0; det < 540; det++){
    if(countDetector[det] > 0) {
      meanDetector[det] = meanDetector[det]/countDetector[det];
      rmsDetector[det] = TMath::Abs(rmsDetector[det]/countDetector[det] - (meanDetector[det]*meanDetector[det]));
    }
  }
  //printf("Put mean value, meanAll %f, rmsAll %f\n",meanAll,rmsAll);
  ///////////////////////////////////////////////
  // Put the mean value for the no-fitted
  /////////////////////////////////////////////  
  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) fVectorFit.At(k))->GetDetector();
    sector = GetSector(detector);
    Int_t rowMax    = fGeo->GetRowMax(GetLayer(detector),GetStack(detector),GetSector(detector));
    Int_t colMax    = fGeo->GetColMax(GetLayer(detector));
    Float_t *coef = ((AliTRDFitInfo *) fVectorFit.At(k))->GetCoef();

    for (Int_t row = 0; row < rowMax; row++) {
      for (Int_t col = 0; col < colMax; col++) {
	value = coef[(Int_t)(col*rowMax+row)];
	if(value < 0.0) {
	  if((ofwhat == 0) && (meanAll > 0.0) && (countAll > 15)) coef[(Int_t)(col*rowMax+row)] = -TMath::Abs(meanAll);
	  if(ofwhat == 1){
	    if((meanDetector[detector] > 0.0) && (countDetector[detector] > 20)) coef[(Int_t)(col*rowMax+row)] = -TMath::Abs(meanDetector[detector]);
	    else if((meanSupermodule[sector] > 0.0) && (countSupermodule[sector] > 15)) coef[(Int_t)(col*rowMax+row)] = -TMath::Abs(meanSupermodule[sector]);
	    else if((meanAll > 0.0) && (countAll > 15)) coef[(Int_t)(col*rowMax+row)] = -TMath::Abs(meanAll);
	  }  
	}
	// Debug
	if(fDebugLevel > 1){
	  
	  if ( !fDebugStreamer ) {
	    //debug stream
	    TDirectory *backup = gDirectory;
	    fDebugStreamer = new TTreeSRedirector("TRDDebugFit.root");
	    if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
	  } 
	  
	  Float_t coefnow      = coef[(Int_t)(col*rowMax+row)]; 
	  
	  (* fDebugStreamer) << "PutMeanValueOtherVectorFit"<<
	    "detector="<<detector<<
	    "sector="<<sector<<
	    "row="<<row<<
	    "col="<<col<<
	    "before="<<value<<
	    "after="<<coefnow<<
	    "\n";  
	}
      } // Col
    } // Row
  }
}
//______________________________________________________________________
void AliTRDCalibraFit::PutMeanValueOtherVectorFit2(Int_t ofwhat, Bool_t perdetector){
  //
  // ofwhat is equaled to 0: mean value of all passing detectors
  // ofwhat is equaled to 1: mean value of the detector, otherwise supermodule, otherwise all
  //

  Int_t loop = (Int_t) fVectorFit2.GetEntriesFast();
  if(loop != 540) {
    AliInfo("The Vector Fit is not complete!");
    return;
  }
  Int_t detector = -1;
  Int_t sector = -1;
  Float_t value  = 0.0;

  /////////////////////////////////
  // Calculate the mean values
  ////////////////////////////////
  // Initialisation
  ////////////////////////
  Double_t meanAll = 0.0;
  Double_t rmsAll = 0.0;
  Double_t meanSupermodule[18];
  Double_t rmsSupermodule[18];
  Double_t meanDetector[540];
  Double_t rmsDetector[540];
  Int_t countAll = 0;
  Int_t countSupermodule[18];
  Int_t countDetector[540];
  for(Int_t sm = 0; sm < 18; sm++){
    rmsSupermodule[sm] = 0.0;
    meanSupermodule[sm] = 0.0;
    countSupermodule[sm] = 0;
  }
  for(Int_t det = 0; det < 540; det++){
    rmsDetector[det] = 0.0;
    meanDetector[det] = 0.0;
    countDetector[det] = 0;
  }
  // compute
  ////////////
  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) fVectorFit2.At(k))->GetDetector();
    sector = GetSector(detector);
    if(perdetector){
      value = ((AliTRDFitInfo *) fVectorFit2.At(k))->GetCoef()[0];
      if(value < 70.0) {
	rmsDetector[detector] += value*value;
	meanDetector[detector] += value;
	countDetector[detector]++;
	rmsSupermodule[sector] += value*value;
	meanSupermodule[sector] += value;
	countSupermodule[sector]++;
	meanAll += value;
	rmsAll += value*value;
	countAll++;
      }
    }
    else {
      Int_t rowMax    = fGeo->GetRowMax(GetLayer(detector),GetStack(detector),GetSector(detector));
      Int_t colMax    = fGeo->GetColMax(GetLayer(detector));
      for (Int_t row = 0; row < rowMax; row++) {
	for (Int_t col = 0; col < colMax; col++) {
	  value = ((AliTRDFitInfo *) fVectorFit2.At(k))->GetCoef()[(Int_t)(col*rowMax+row)];
	  if(value < 70.0) {
	    rmsDetector[detector] += value*value;
	    meanDetector[detector] += value;
	    countDetector[detector]++;
	    rmsSupermodule[sector] += value*value;
	    meanSupermodule[sector] += value;
	    countSupermodule[sector]++;
	    rmsAll += value*value;
	    meanAll += value;
	    countAll++;
	  }
	  
	} // Col
      } // Row
    }
  }  
  if(countAll > 0) {
    meanAll = meanAll/countAll;
    rmsAll = TMath::Abs(rmsAll/countAll - (meanAll*meanAll));
  }
  for(Int_t sm = 0; sm < 18; sm++){
    if(countSupermodule[sm] > 0) {
      meanSupermodule[sm] = meanSupermodule[sm]/countSupermodule[sm];
      rmsSupermodule[sm] = TMath::Abs(rmsSupermodule[sm]/countSupermodule[sm] - (meanSupermodule[sm]*meanSupermodule[sm]));
    }
  }
  for(Int_t det = 0; det < 540; det++){
    if(countDetector[det] > 0) {
      meanDetector[det] = meanDetector[det]/countDetector[det];
      rmsDetector[det] = TMath::Abs(rmsDetector[det]/countDetector[det] - (meanDetector[det]*meanDetector[det]));
    }
  }
  //printf("Put mean value 2: meanAll %f, rmsAll %f\n",meanAll,rmsAll);
  ////////////////////////////////////////////
  // Put the mean value for the no-fitted
  /////////////////////////////////////////////  
  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) fVectorFit2.At(k))->GetDetector();
    sector = GetSector(detector);
    Int_t rowMax    = fGeo->GetRowMax(GetLayer(detector),GetStack(detector),GetSector(detector));
    Int_t colMax    = fGeo->GetColMax(GetLayer(detector));
    Float_t *coef = ((AliTRDFitInfo *) fVectorFit2.At(k))->GetCoef();

    for (Int_t row = 0; row < rowMax; row++) {
      for (Int_t col = 0; col < colMax; col++) {
	value = coef[(Int_t)(col*rowMax+row)];
	if(value > 70.0) {
	  if((ofwhat == 0) && (meanAll > -3.0) && (countAll > 15)) coef[(Int_t)(col*rowMax+row)] = meanAll+100.0;
	  if(ofwhat == 1){
	    if((meanDetector[detector] > -3.0) && (countDetector[detector] > 20)) coef[(Int_t)(col*rowMax+row)] = meanDetector[detector]+100.0;
	    else if((meanSupermodule[sector] > -3.0) && (countSupermodule[sector] > 15)) coef[(Int_t)(col*rowMax+row)] = meanSupermodule[sector]+100.0;
	    else if((meanAll > -3.0) && (countAll > 15)) coef[(Int_t)(col*rowMax+row)] = meanAll+100.0;
	  }  
	}
	// Debug
	if(fDebugLevel > 1){
	  
	  if ( !fDebugStreamer ) {
	    //debug stream
	    TDirectory *backup = gDirectory;
	    fDebugStreamer = new TTreeSRedirector("TRDDebugFit.root");
	    if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
	  } 
	  
	  Float_t coefnow      = coef[(Int_t)(col*rowMax+row)]; 
	  
	  (* fDebugStreamer) << "PutMeanValueOtherVectorFit2"<<
	    "detector="<<detector<<
	    "sector="<<sector<<
	    "row="<<row<<
	    "col="<<col<<
	    "before="<<value<<
	    "after="<<coefnow<<
	    "\n";  
	}
      } // Col
    } // Row
  }
  
}
//_____________________________________________________________________________
AliTRDCalDet *AliTRDCalibraFit::CreateDetObjectVdrift(const TObjArray *vectorFit, Bool_t perdetector)
{
  //
  // It creates the AliTRDCalDet object from the AliTRDFitInfo
  // It takes the mean value of the coefficients per detector 
  // This object has to be written in the database
  //
  
  // Create the DetObject
  AliTRDCalDet *object = new AliTRDCalDet("ChamberVdrift","TRD drift velocities (detector value)");

  Int_t loop = (Int_t) vectorFit->GetEntriesFast();
  if(loop != 540) AliInfo("The Vector Fit is not complete!");
  Int_t detector = -1;
  Float_t value  = 0.0;
  
  //
  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) vectorFit->At(k))->GetDetector();
    Float_t mean  = 0.0;
    if(perdetector){
      mean = TMath::Abs(((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[0]);
    }
    else {
      Int_t   count = 0;
      Int_t rowMax    = fGeo->GetRowMax(GetLayer(detector),GetStack(detector),GetSector(detector));
      Int_t colMax    = fGeo->GetColMax(GetLayer(detector));
      for (Int_t row = 0; row < rowMax; row++) {
	for (Int_t col = 0; col < colMax; col++) {
	  value = ((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[(Int_t)(col*rowMax+row)];
	  mean += TMath::Abs(value);
	  count++;       
	} // Col
      } // Row
      if(count > 0) mean = mean/count;
    }
    object->SetValue(detector,mean);
  }
  
  return object;
}
//_____________________________________________________________________________
AliTRDCalDet *AliTRDCalibraFit::CreateDetObjectGain(const TObjArray *vectorFit, Bool_t meanOtherBefore, Double_t scaleFitFactor, Bool_t perdetector)
{
  //
  // It creates the AliTRDCalDet object from the AliTRDFitInfo
  // It takes the mean value of the coefficients per detector 
  // This object has to be written in the database
  //
  
  // Create the DetObject
  AliTRDCalDet *object = new AliTRDCalDet("ChamberGainFactor","GainFactor (detector value)");
  
  fScaleGain = scaleFitFactor;
 
  Int_t loop = (Int_t) vectorFit->GetEntriesFast();
  if(loop != 540) AliInfo("The Vector Fit is not complete!");
  Int_t detector = -1;
  Float_t value  = 0.0;

  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) vectorFit->At(k))->GetDetector();  
    Float_t mean  = 0.0;
    if(perdetector){
      value = ((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[0];
      if(!meanOtherBefore){
	if(value > 0) value = value*scaleFitFactor;
      }
      else value = value*scaleFitFactor;
      mean = TMath::Abs(value);
    }
    else{
      Int_t   count = 0;
      Int_t rowMax    = fGeo->GetRowMax(GetLayer(detector),GetStack(detector),GetSector(detector));
      Int_t colMax    = fGeo->GetColMax(GetLayer(detector));
      for (Int_t row = 0; row < rowMax; row++) {
	for (Int_t col = 0; col < colMax; col++) {
	  value = ((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[(Int_t)(col*rowMax+row)];
	  if(!meanOtherBefore) {
	    if(value > 0) value = value*scaleFitFactor;
	  }
	  else value = value*scaleFitFactor;
	  mean += TMath::Abs(value);
	  count++;       
	} // Col
      } // Row
      if(count > 0) mean = mean/count;
    }
    if(mean < 0.1) mean = 0.1;
    object->SetValue(detector,mean);
  }
 
  return object;
}
//_____________________________________________________________________________
AliTRDCalDet *AliTRDCalibraFit::CreateDetObjectT0(const TObjArray *vectorFit, Bool_t perdetector)
{
  //
  // It creates the AliTRDCalDet object from the AliTRDFitInfo2
  // It takes the min value of the coefficients per detector 
  // This object has to be written in the database
  //
  
  // Create the DetObject
  AliTRDCalDet *object = new AliTRDCalDet("ChamberT0","T0 (detector value)");
  
  Int_t loop = (Int_t) vectorFit->GetEntriesFast();
  if(loop != 540) AliInfo("The Vector Fit is not complete!");
  Int_t detector = -1;
  Float_t value  = 0.0;

  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) vectorFit->At(k))->GetDetector();   
    Float_t min  = 100.0;
    if(perdetector){
      value = ((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[0];
      //printf("Create det object %f for %d\n",value,k);
      // check successful
      if(value > 70.0) value = value-100.0;
      //
      min = value;
    }
    else{
      Int_t rowMax    = fGeo->GetRowMax(GetLayer(detector),GetStack(detector),GetSector(detector));
      Int_t colMax    = fGeo->GetColMax(GetLayer(detector));
      for (Int_t row = 0; row < rowMax; row++) {
	for (Int_t col = 0; col < colMax; col++) {
	  value = ((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[(Int_t)(col*rowMax+row)];
	  // check successful
	  if(value > 70.0) value = value-100.0;
	  //
	  if(min > value) min = value;
	} // Col
      } // Row
    }
    object->SetValue(detector,min);
  }

  return object;

}
//_____________________________________________________________________________
AliTRDCalDet *AliTRDCalibraFit::CreateDetObjectLorentzAngle(const TObjArray *vectorFit)
{
  //
  // It creates the AliTRDCalDet object from the AliTRDFitInfo2
  // It takes the min value of the coefficients per detector 
  // This object has to be written in the database
  //
  
  // Create the DetObject
  AliTRDCalDet *object = new AliTRDCalDet("tan(lorentzangle)","tan(lorentzangle) (detector value)");
  
  
  Int_t loop = (Int_t) vectorFit->GetEntriesFast();
  if(loop != 540) AliInfo("The Vector Fit is not complete!");
  Int_t detector = -1;
  Float_t value  = 0.0;

  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) vectorFit->At(k))->GetDetector();
    /*
      Int_t rowMax    = fGeo->GetRowMax(GetLayer(detector),GetStack(detector),GetSector(detector));
      Int_t colMax    = fGeo->GetColMax(GetLayer(detector));
      Float_t min  = 100.0;
      for (Int_t row = 0; row < rowMax; row++) {
      for (Int_t col = 0; col < colMax; col++) {
      value = ((AliTRDFitInfo *) fVectorFit2.At(k))->GetCoef()[(Int_t)(col*rowMax+row)];
      mean += -TMath::Abs(value);
      count++;       
      } // Col
      } // Row
      if(count > 0) mean = mean/count;
    */
    value = ((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[0];
    if(value > 70.0) value = value-100.0;
    object->SetValue(detector,value);
  }

  return object;
  
}
//_____________________________________________________________________________
TObject *AliTRDCalibraFit::CreatePadObjectGain(const TObjArray *vectorFit, Double_t scaleFitFactor, const AliTRDCalDet *detobject)
{
  //
  // It Creates the AliTRDCalPad object from AliTRDFitInfo
  // You need first to create the object for the detectors,
  // where the mean value is put.
  // This object has to be written in the database
  //
  
  // Create the DetObject
  AliTRDCalPad *object = new AliTRDCalPad("GainFactor","GainFactor (local variations)");
  
  if(!vectorFit){
    for(Int_t k = 0; k < 540; k++){
      AliTRDCalROC *calROC = object->GetCalROC(k);
      Int_t nchannels = calROC->GetNchannels();
      for(Int_t ch = 0; ch < nchannels; ch++){
	calROC->SetValue(ch,1.0);
      }
    }
  }
  else{

    Int_t loop = (Int_t) vectorFit->GetEntriesFast();
    if(loop != 540) AliInfo("The Vector Fit is not complete!");
    Int_t detector = -1;
    Float_t value  = 0.0;
    
    for (Int_t k = 0; k < loop; k++) {
      detector  = ((AliTRDFitInfo *) vectorFit->At(k))->GetDetector();
      AliTRDCalROC *calROC = object->GetCalROC(detector);
      Float_t mean         = detobject->GetValue(detector);
      if(TMath::Abs(mean) <= 0.0000000001) continue;
      Int_t rowMax    = calROC->GetNrows();
      Int_t colMax    = calROC->GetNcols();
      for (Int_t row = 0; row < rowMax; row++) {
	for (Int_t col = 0; col < colMax; col++) {
	  value = ((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[(Int_t)(col*rowMax+row)];
	  if(value > 0) value = value*scaleFitFactor;
	  calROC->SetValue(col,row,TMath::Abs(value)/mean);
	} // Col
      } // Row
    } 
  }

  return object;  
}
//_____________________________________________________________________________
TObject *AliTRDCalibraFit::CreatePadObjectVdrift(const TObjArray *vectorFit, const AliTRDCalDet *detobject)
{
  //
  // It Creates the AliTRDCalPad object from AliTRDFitInfo
  // You need first to create the object for the detectors,
  // where the mean value is put.
  // This object has to be written in the database
  //

  // Create the DetObject
  AliTRDCalPad *object = new AliTRDCalPad("LocalVdrift","TRD drift velocities (local variations)");

  if(!vectorFit){
    for(Int_t k = 0; k < 540; k++){
      AliTRDCalROC *calROC = object->GetCalROC(k);
      Int_t nchannels = calROC->GetNchannels();
      for(Int_t ch = 0; ch < nchannels; ch++){
	calROC->SetValue(ch,1.0);
      }
    }
  }
  else {
    
    Int_t loop = (Int_t) vectorFit->GetEntriesFast();
    if(loop != 540) AliInfo("The Vector Fit is not complete!");
    Int_t detector = -1;
    Float_t value  = 0.0;
    
    for (Int_t k = 0; k < loop; k++) {
      detector  = ((AliTRDFitInfo *) vectorFit->At(k))->GetDetector();
      AliTRDCalROC *calROC = object->GetCalROC(detector);
      Float_t mean         = detobject->GetValue(detector);
      if(mean == 0) continue;
      Int_t rowMax    = calROC->GetNrows();
      Int_t colMax    = calROC->GetNcols();
      for (Int_t row = 0; row < rowMax; row++) {
	for (Int_t col = 0; col < colMax; col++) {
	  value = ((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[(Int_t)(col*rowMax+row)];
	  calROC->SetValue(col,row,TMath::Abs(value)/mean);
	} // Col
      } // Row
    } 
  }
  return object;    

}
//_____________________________________________________________________________
TObject *AliTRDCalibraFit::CreatePadObjectT0(const TObjArray *vectorFit, const AliTRDCalDet *detobject)
{
  //
  // It Creates the AliTRDCalPad object from AliTRDFitInfo2
  // You need first to create the object for the detectors,
  // where the mean value is put.
  // This object has to be written in the database
  //
  
  // Create the DetObject
  AliTRDCalPad *object = new AliTRDCalPad("LocalT0","T0 (local variations)");

  if(!vectorFit){
    for(Int_t k = 0; k < 540; k++){
      AliTRDCalROC *calROC = object->GetCalROC(k);
      Int_t nchannels = calROC->GetNchannels();
      for(Int_t ch = 0; ch < nchannels; ch++){
	calROC->SetValue(ch,0.0);
      }
    }
  }
  else {
    
    Int_t loop = (Int_t) vectorFit->GetEntriesFast();
    if(loop != 540) AliInfo("The Vector Fit is not complete!");
    Int_t detector = -1;
    Float_t value  = 0.0;
    
    for (Int_t k = 0; k < loop; k++) {
      detector  = ((AliTRDFitInfo *) vectorFit->At(k))->GetDetector();
      AliTRDCalROC *calROC = object->GetCalROC(detector);
      Float_t min          = detobject->GetValue(detector);
      Int_t rowMax    = calROC->GetNrows();
      Int_t colMax    = calROC->GetNcols();
      for (Int_t row = 0; row < rowMax; row++) {
	for (Int_t col = 0; col < colMax; col++) {
	  value = ((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[(Int_t)(col*rowMax+row)];
	  // check successful
	  if(value > 70.0) value = value - 100.0;
	  //
	  calROC->SetValue(col,row,value-min);
	} // Col
      } // Row
    } 
  }
  return object;    

}
//_____________________________________________________________________________
TObject *AliTRDCalibraFit::CreatePadObjectPRF(const TObjArray *vectorFit)
{
  //
  // It Creates the AliTRDCalPad object from AliTRDFitInfo
  // This object has to be written in the database
  //
  
  // Create the DetObject
  AliTRDCalPad *object = new AliTRDCalPad("PRFWidth","PRFWidth");

  Int_t loop = (Int_t) vectorFit->GetEntriesFast();
  if(loop != 540) AliInfo("The Vector Fit is not complete!");
  Int_t detector = -1;
  Float_t value  = 0.0;

  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) vectorFit->At(k))->GetDetector();
    AliTRDCalROC *calROC = object->GetCalROC(detector);
    Int_t rowMax    = calROC->GetNrows();
    Int_t colMax    = calROC->GetNcols();
    for (Int_t row = 0; row < rowMax; row++) {
      for (Int_t col = 0; col < colMax; col++) {
	value = ((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[(Int_t)(col*rowMax+row)];
       	calROC->SetValue(col,row,TMath::Abs(value));
      } // Col
    } // Row
  } 

  return object;  

}
//_____________________________________________________________________________
AliTRDCalDet *AliTRDCalibraFit::MakeOutliersStatDet(const TObjArray *vectorFit, const char *name, Double_t &mean)
{
  //
  // It Creates the AliTRDCalDet object from AliTRDFitInfo
  // 0 successful fit 1 not successful fit
  // mean is the mean value over the successful fit
  // do not use it for t0: no meaning
  //
  
  // Create the CalObject
  AliTRDCalDet *object = new AliTRDCalDet(name,name);
  mean = 0.0;
  Int_t count = 0;
  
  Int_t loop = (Int_t) vectorFit->GetEntriesFast();
  if(loop != 540) {
    AliInfo("The Vector Fit is not complete! We initialise all outliers");
    for(Int_t k = 0; k < 540; k++){
      object->SetValue(k,1.0);
    }
  }
  Int_t detector = -1;
  Float_t value  = 0.0;
  
  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) vectorFit->At(k))->GetDetector();
    value = ((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[0];
    if(value <= 0) object->SetValue(detector,1.0);
    else {
      object->SetValue(detector,0.0);
      mean += value;
      count++;
    }
  }
  if(count > 0) mean /= count;
  return object;  
}
//_____________________________________________________________________________
TObject *AliTRDCalibraFit::MakeOutliersStatPad(const TObjArray *vectorFit, const char *name, Double_t &mean)
{
  //
  // It Creates the AliTRDCalPad object from AliTRDFitInfo
  // 0 not successful fit 1 successful fit
  // mean mean value over the successful fit
  //
  
  // Create the CalObject
  AliTRDCalPad *object = new AliTRDCalPad(name,name);
  mean = 0.0;
  Int_t count = 0;
  
  Int_t loop = (Int_t) vectorFit->GetEntriesFast();
  if(loop != 540) {
    AliInfo("The Vector Fit is not complete! We initialise all outliers");
    for(Int_t k = 0; k < 540; k++){
      AliTRDCalROC *calROC = object->GetCalROC(k);
      Int_t nchannels = calROC->GetNchannels();
      for(Int_t ch = 0; ch < nchannels; ch++){
	calROC->SetValue(ch,1.0);
      }
    }
  }
  Int_t detector = -1;
  Float_t value  = 0.0;
  
  for (Int_t k = 0; k < loop; k++) {
    detector  = ((AliTRDFitInfo *) vectorFit->At(k))->GetDetector();
    AliTRDCalROC *calROC = object->GetCalROC(detector);
    Int_t nchannels    = calROC->GetNchannels();
    for (Int_t ch = 0; ch < nchannels; ch++) {
      value = ((AliTRDFitInfo *) vectorFit->At(k))->GetCoef()[ch];
      if(value <= 0) calROC->SetValue(ch,1.0);
      else {
	calROC->SetValue(ch,0.0);
	mean += value;
	count++;
      }
    } // channels
  }
  if(count > 0) mean /= count;
  return object;  
}
//_____________________________________________________________________________
void AliTRDCalibraFit::SetPeriodeFitPH(Int_t periodeFitPH)
{ 
  //
  // Set FitPH if 1 then each detector will be fitted
  //

  if (periodeFitPH > 0) {
    fFitPHPeriode   = periodeFitPH; 
  }
  else {
    AliInfo("periodeFitPH must be higher than 0!");
  }

}
//_____________________________________________________________________________
void AliTRDCalibraFit::SetBeginFitCharge(Float_t beginFitCharge)
{ 
  //
  // The fit of the deposited charge distribution begins at
  // histo->Mean()/beginFitCharge
  // You can here set beginFitCharge
  //

  if (beginFitCharge > 0) {
    fBeginFitCharge = beginFitCharge; 
  }
  else {
    AliInfo("beginFitCharge must be strict positif!");
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::SetT0Shift0(Float_t t0Shift) 
{ 
  //
  // The t0 calculated with the maximum positif slope is shift from t0Shift0
  // You can here set t0Shift0
  //

  if (t0Shift > 0) {
    fT0Shift0 = t0Shift; 
  } 
  else {
    AliInfo("t0Shift0 must be strict positif!");
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::SetT0Shift1(Float_t t0Shift) 
{ 
  //
  // The t0 calculated with the maximum of the amplification region is shift from t0Shift1
  // You can here set t0Shift1
  //

  if (t0Shift > 0) {
    fT0Shift1 = t0Shift; 
  } 
  else {
    AliInfo("t0Shift must be strict positif!");
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::SetRangeFitPRF(Float_t rangeFitPRF)
{ 
  //
  // The fit of the PRF is from -rangeFitPRF to rangeFitPRF
  // You can here set rangeFitPRF
  //

  if ((rangeFitPRF >    0) && 
      (rangeFitPRF <= 1.5)) {
    fRangeFitPRF = rangeFitPRF;
  } 
  else {
    AliInfo("rangeFitPRF must be between 0 and 1.0");
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::SetMinEntries(Int_t minEntries)
{ 
  //
  // Minimum entries for fitting
  //

  if (minEntries >    0) {
    fMinEntries = minEntries;
  } 
  else {
    AliInfo("fMinEntries must be >= 0.");
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::SetRebin(Short_t rebin)
{ 
  //
  // Rebin with rebin time less bins the Ch histo
  // You can set here rebin that should divide the number of bins of CH histo
  //

  if (rebin > 0) {
    fRebin = rebin; 
    AliInfo("You have to be sure that fRebin divides fNumberBinCharge used!");
  } 
  else {
    AliInfo("You have to choose a positiv value!");
  }

}
//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::FillVectorFit()
{
  //
  // For the Fit functions fill the vector Fit
  //

  AliTRDFitInfo *fitInfo = new AliTRDFitInfo();

  Int_t ntotal = 1;
  if (GetStack(fCountDet) == 2) {
    ntotal = 1728;
  }
  else {
    ntotal = 2304;
  }

  //printf("For the detector %d , ntotal %d and fCoefCH[0] %f\n",countdet,ntotal,fCoefCH[0]);
  Float_t *coef = new Float_t[ntotal];
  for (Int_t i = 0; i < ntotal; i++) {
    coef[i] = fCurrentCoefDetector[i];
  }
  
  Int_t detector = fCountDet;
  // Set
  fitInfo->SetCoef(coef);
  fitInfo->SetDetector(detector);
  fVectorFit.Add((TObject *) fitInfo);

  return kTRUE;

}
//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::FillVectorFit2()
{
  //
  // For the Fit functions fill the vector Fit
  //

  AliTRDFitInfo *fitInfo = new AliTRDFitInfo();

  Int_t ntotal = 1;
  if (GetStack(fCountDet) == 2) {
    ntotal = 1728;
  }
  else {
    ntotal = 2304;
  }

  //printf("For the detector %d , ntotal %d and fCoefCH[0] %f\n",countdet,ntotal,fCoefCH[0]);
  Float_t *coef = new Float_t[ntotal];
  for (Int_t i = 0; i < ntotal; i++) {
    coef[i] = fCurrentCoefDetector2[i];
  }
  
  Int_t detector = fCountDet;
  // Set
  fitInfo->SetCoef(coef);
  fitInfo->SetDetector(detector);
  fVectorFit2.Add((TObject *) fitInfo);

  return kTRUE;

}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::InitFit(Int_t nbins, Int_t i)
{
  //
  // Init the number of expected bins and fDect1[i] fDect2[i] 
  //

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.01);
  
  // Mode groups of pads: the total number of bins!
  CalculNumberOfBinsExpected(i);
  
  // Quick verification that we have the good pad calibration mode!
  if (fNumberOfBinsExpected != nbins) {
    AliInfo(Form("It doesn't correspond to the mode of pad group calibration: expected %d and seen %d!",fNumberOfBinsExpected,nbins));
    return kFALSE;
  }
  
  // Security for fDebug 3 and 4
  if ((fDebugLevel >= 3) && 
      ((fDet[0] >  5) || 
       (fDet[1] >  4) || 
       (fDet[2] > 17))) {
    AliInfo("This detector doesn't exit!");
    return kFALSE;
  }

  // Determine fDet1 and fDet2 and set the fNfragZ and fNfragRphi for debug 3 and 4
  CalculDect1Dect2(i);

 
  return kTRUE;
}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::InitFitCH()
{
  //
  // Init the fVectorFitCH for normalisation
  // Init the histo for debugging 
  //

  gDirectory = gROOT;
 
  fScaleFitFactor = 0.0;
  fCurrentCoefDetector   = new Float_t[2304];
  for (Int_t k = 0; k < 2304; k++) {
    fCurrentCoefDetector[k] = 0.0;    
  }
  fVectorFit.SetName("gainfactorscoefficients");

  // fDebug == 0 nothing
  // fDebug == 1 and fFitVoir no histo
  if (fDebugLevel == 1) {
    if(!CheckFitVoir()) return kFALSE;
  }
  //Get the CalDet object
  if(fAccCDB){
    AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
    if (!cal) {
      AliInfo("Could not get calibDB");
      return kFALSE;
    }
    if(fCalDet) delete fCalDet;
    fCalDet = new AliTRDCalDet(*(cal->GetGainFactorDet()));
  }
  else{
    Float_t devalue = 1.0;
    if(fCalDet) delete fCalDet;
    fCalDet = new AliTRDCalDet("ChamberGainFactor","GainFactor (detector value)");
    for(Int_t k = 0; k < 540; k++){
      fCalDet->SetValue(k,devalue);
    }
  }
  return kTRUE;
  
}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::InitFitPH()
{
  //
  // Init the arrays of results 
  // Init the histos for debugging 
  //

  gDirectory = gROOT;
  fVectorFit.SetName("driftvelocitycoefficients");
  fVectorFit2.SetName("t0coefficients");

  fCurrentCoefDetector   = new Float_t[2304];
  for (Int_t k = 0; k < 2304; k++) {
    fCurrentCoefDetector[k] = 0.0;    
  }

  fCurrentCoefDetector2   = new Float_t[2304];
  for (Int_t k = 0; k < 2304; k++) {
    fCurrentCoefDetector2[k] = 0.0;    
  }
 
  //fDebug == 0 nothing
  // fDebug == 1 and fFitVoir no histo
  if (fDebugLevel == 1) {
    if(!CheckFitVoir()) return kFALSE;
  }
  //Get the CalDet object
  if(fAccCDB){
    AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
    if (!cal) {
      AliInfo("Could not get calibDB");
      return kFALSE;
    }
    if(fCalDet) delete fCalDet;
    if(fCalDet2) delete fCalDet2;
    fCalDet  = new AliTRDCalDet(*(cal->GetVdriftDet()));
    fCalDet2 = new AliTRDCalDet(*(cal->GetT0Det())); 
  }
  else{
    Float_t devalue  = 1.5;
    Float_t devalue2 = 0.0; 
    if(fCalDet) delete fCalDet;
    if(fCalDet2) delete fCalDet2;
    fCalDet  = new AliTRDCalDet("ChamberVdrift","TRD drift velocities (detector value)");
    fCalDet2 = new AliTRDCalDet("ChamberT0","T0 (detector value)");
    for(Int_t k = 0; k < 540; k++){
      fCalDet->SetValue(k,devalue);
      fCalDet2->SetValue(k,devalue2);
    }
  }
  return kTRUE;
}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::InitFitPRF()
{
  //
  // Init the calibration mode (Nz, Nrphi), the histograms for
  // debugging the fit methods if fDebug > 0, 
  //
  
  gDirectory = gROOT;
  fVectorFit.SetName("prfwidthcoefficients");
 
  fCurrentCoefDetector   = new Float_t[2304];
  for (Int_t k = 0; k < 2304; k++) {
    fCurrentCoefDetector[k] = 0.0;    
  }
  
  // fDebug == 0 nothing
  // fDebug == 1 and fFitVoir no histo
  if (fDebugLevel == 1) {
    if(!CheckFitVoir()) return kFALSE;
  }
  return kTRUE;
}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::InitFitLinearFitter()
{
  //
  // Init the fCalDet, fVectorFit fCurrentCoefDetector 
  //
  
  gDirectory = gROOT;
 
  fCurrentCoefDetector   = new Float_t[2304];
  fCurrentCoefDetector2  = new Float_t[2304];
  for (Int_t k = 0; k < 2304; k++) {
    fCurrentCoefDetector[k]  = 0.0;
    fCurrentCoefDetector2[k] = 0.0;    
  }

  if((!fCalDetVdriftUsed) || (!fCalDetExBUsed)) return kFALSE; 

  return kTRUE;
}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
void AliTRDCalibraFit::InitfCountDetAndfCount(Int_t i)
{
  //
  // Init the current detector where we are fCountDet and the
  // next fCount for the functions Fit... 
  //

  // Loop on the Xbins of ch!!
  fCountDet = -1; // Current detector
  fCount    =  0; // To find the next detector
  
  // If fDebug >= 3
  if (fDebugLevel >= 3) {
    // Set countdet to the detector
    fCountDet = AliTRDgeometry::GetDetector(fDet[0],fDet[1],fDet[2]);
    // Set counter to write at the end of the detector
    fCount = fDect2;
    // Get the right calib objects
    SetCalROC(i);
  }
  if(fDebugLevel == 1) {
    fCountDet = 0;
    fCalibraMode->CalculXBins(fCountDet,i);
    if((fCalibraMode->GetNz(i)!=100) && (fCalibraMode->GetNrphi(i)!=100)){
      while(fCalibraMode->GetXbins(i) <=fFitVoir){
	fCountDet++;
	fCalibraMode->CalculXBins(fCountDet,i);
	//printf("GetXBins %d\n",fCalibraMode->GetXbins(i));
      }      
    }
    else {
      fCountDet++;
    }
    fCount    = fCalibraMode->GetXbins(i);
    fCountDet--;
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    fCalibraMode->ModePadCalibration((Int_t) GetStack(fCountDet),i);
    fCalibraMode->ModePadFragmentation((Int_t) GetLayer(fCountDet)
				      ,(Int_t) GetStack(fCountDet)
				      ,(Int_t) GetSector(fCountDet),i);
  }
}
//_______________________________________________________________________________
void AliTRDCalibraFit::CalculNumberOfBinsExpected(Int_t i)
{
  //
  // Calculate the number of bins expected (calibration groups)
  //
  
  fNumberOfBinsExpected = 0;
  // All
  if((fCalibraMode->GetNz(i) == 100) && (fCalibraMode->GetNrphi(i) == 100)){
    fNumberOfBinsExpected = 1;
    return;
  }
  // Per supermodule
  if((fCalibraMode->GetNz(i) == 10) && (fCalibraMode->GetNrphi(i) == 10)){
    fNumberOfBinsExpected = 18;
    return;
  }
  // More
  fCalibraMode->ModePadCalibration(2,i);
  fCalibraMode->ModePadFragmentation(0,2,0,i);
  fCalibraMode->SetDetChamb2(i);
  if (fDebugLevel > 1) {
    AliInfo(Form("For the chamber 2: %d",fCalibraMode->GetDetChamb2(i)));
  }
  fNumberOfBinsExpected += 6 * 18 * fCalibraMode->GetDetChamb2(i);
  fCalibraMode->ModePadCalibration(0,i);
  fCalibraMode->ModePadFragmentation(0,0,0,i);
  fCalibraMode->SetDetChamb0(i);
  if (fDebugLevel > 1) {
    AliInfo(Form("For the other chamber 0: %d",fCalibraMode->GetDetChamb0(i)));
  }
  fNumberOfBinsExpected += 6 * 4 * 18 * fCalibraMode->GetDetChamb0(i);
 
}
//_______________________________________________________________________________
void AliTRDCalibraFit::CalculDect1Dect2(Int_t i)
{
  //
  // Calculate the range of fits
  //
  
  fDect1 = -1;
  fDect2 = -1;
  if (fDebugLevel == 1) {
    fDect1 = fFitVoir;
    fDect2 = fDect1 +1;
  }
  if ((fDebugLevel == 2) || (fDebugLevel == 0)) {
    fDect1 = 0;
    fDect2 = fNumberOfBinsExpected;
  }
  if (fDebugLevel >= 3) {
    fCountDet = AliTRDgeometry::GetDetector(fDet[0],fDet[1],fDet[2]);
    fCalibraMode->CalculXBins(fCountDet,i);
    fDect1 = fCalibraMode->GetXbins(i);
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    fCalibraMode->ModePadCalibration((Int_t) GetStack(fCountDet),i);
    fCalibraMode->ModePadFragmentation((Int_t) GetLayer(fCountDet)
				      ,(Int_t) GetStack(fCountDet)
				      ,(Int_t) GetSector(fCountDet),i);
    // Set for the next detector
    fDect2 = fDect1 + fCalibraMode->GetNfragZ(i)*fCalibraMode->GetNfragRphi(i);
  }
}
//_______________________________________________________________________________
Bool_t AliTRDCalibraFit::CheckFitVoir()
{
  //
  // Check if fFitVoir is in the range
  //
  
  if (fFitVoir < fNumberOfBinsExpected) {
    AliInfo(Form("We will see the fit of the object %d",fFitVoir));
  }
  else {
    AliInfo("fFitVoir is out of range of the histo!");
    return kFALSE;
  }
  return kTRUE;
}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
void AliTRDCalibraFit::UpdatefCountDetAndfCount(Int_t idect, Int_t i)
{
  //
  // See if we are in a new detector and update the
  // variables fNfragZ and fNfragRphi if yes 
  // Will never happen for only one detector (3 and 4)
  // Doesn't matter for 2
  //
  if (fCount == idect) {
    // On en est au detector (or first detector in the group)
    fCountDet += 1;
    AliDebug(2,Form("We are at the detector %d\n",fCountDet));
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    fCalibraMode->ModePadCalibration((Int_t) GetStack(fCountDet),i);
    fCalibraMode->ModePadFragmentation((Int_t) GetLayer(fCountDet)
			     	       ,(Int_t) GetStack(fCountDet)
                                       ,(Int_t) GetSector(fCountDet),i);
    // Set for the next detector
    fCount += fCalibraMode->GetNfragZ(i)*fCalibraMode->GetNfragRphi(i);
    // calib objects
    SetCalROC(i);
  }
}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
void AliTRDCalibraFit::ReconstructFitRowMinRowMax(Int_t idect, Int_t i)
{
  //
  // Reconstruct the min pad row, max pad row, min pad col and
  // max pad col of the calibration group for the Fit functions
  // idect is the calibration group inside the detector
  //
  if (fDebugLevel !=  1) {
    fCalibraMode->ReconstructionRowPadGroup((Int_t) (idect-(fCount-(fCalibraMode->GetNfragZ(i)*fCalibraMode->GetNfragRphi(i)))),i);
  }
  AliDebug(2,Form("AliTRDCalibraFit::ReconstructFitRowMinRowMax: the local calibration group is %d",idect-(fCount-(fCalibraMode->GetNfragZ(i)*fCalibraMode->GetNfragRphi(i)))));
  AliDebug(2,Form("AliTRDCalibraFit::ReconstructFitRowMinRowMax: the number of group per detector is %d",fCalibraMode->GetNfragZ(i)*fCalibraMode->GetNfragRphi(i)));
}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::NotEnoughStatisticCH(Int_t idect)
{
  //
  // For the case where there are not enough entries in the histograms
  // of the calibration group, the value present in the choosen database
  // will be put. A negativ sign enables to know that a fit was not possible.
  //
  
  if (fDebugLevel == 1) {
    AliInfo("The element has not enough statistic to be fitted");
  }
  else if (fNbDet > 0){
    Int_t firstdetector = fCountDet;
    Int_t lastdetector  = fCountDet+fNbDet;
    //AliInfo(Form("The element %d containing the detectors %d to %d has not enough statistic to be fitted",idect,firstdetector,lastdetector));
    // loop over detectors
    for(Int_t det = firstdetector; det < lastdetector; det++){

      //Set the calibration object again
      fCountDet = det;
      SetCalROC(0);   

      // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
      // Put them at 1
      fCalibraMode->ModePadCalibration((Int_t) GetStack(fCountDet),0);
      fCalibraMode->ModePadFragmentation((Int_t) GetLayer(fCountDet)
					 ,(Int_t) GetStack(fCountDet)
					 ,(Int_t) GetSector(fCountDet),0);
      // Reconstruct row min row max
      ReconstructFitRowMinRowMax(idect,0);      

      // Calcul the coef from the database choosen for the detector
      CalculChargeCoefMean(kFALSE);
      
      //stack 2, not stack 2
      Int_t factor = 0;
      if(GetStack(fCountDet) == 2) factor = 12;
      else factor = 16;
      
      // Fill the fCurrentCoefDetector with negative value to say: not fitted
      for (Int_t k = fCalibraMode->GetRowMin(0); k < fCalibraMode->GetRowMax(0); k++) {
	for (Int_t j = fCalibraMode->GetColMin(0); j < fCalibraMode->GetColMax(0); j++) {
	  fCurrentCoefDetector[(Int_t)(j*factor+k)] = -TMath::Abs(fCurrentCoef[1]);
	}
      }
      
      //Put default value negative
      fCurrentCoef[0] = -TMath::Abs(fCurrentCoef[1]);
      fCurrentCoefE   = 0.0;
      
      // Fill the stuff
      FillVectorFit();
      // Debug
      if(fDebugLevel > 1){ 
	
	if ( !fDebugStreamer ) {
	  //debug stream
	  TDirectory *backup = gDirectory;
	  fDebugStreamer = new TTreeSRedirector("TRDDebugFitCH.root");
	  if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
	} 
	
	Int_t   detector   = fCountDet;
	Int_t   caligroup  = idect;
	Short_t rowmin     = fCalibraMode->GetRowMin(0);
	Short_t rowmax     = fCalibraMode->GetRowMax(0);
	Short_t colmin     = fCalibraMode->GetColMin(0);
	Short_t colmax     = fCalibraMode->GetColMax(0);
	Float_t gf         = fCurrentCoef[0]; 
	Float_t gfs        = fCurrentCoef[1]; 
	Float_t gfE        = fCurrentCoefE;
	
	(*fDebugStreamer) << "FillFillCH" <<
	  "detector=" << detector <<
	  "caligroup=" << caligroup <<
	  "rowmin=" << rowmin <<
	  "rowmax=" << rowmax <<
	  "colmin=" << colmin <<
	  "colmax=" << colmax <<
	  "gf=" << gf <<
	  "gfs=" << gfs <<
	  "gfE=" << gfE <<
	  "\n"; 
	
      }
      // Reset
      for (Int_t k = 0; k < 2304; k++) {
	fCurrentCoefDetector[k] = 0.0;
      }
      
    }// loop detector
    AliDebug(2,Form("Check the count now: fCountDet %d",fCountDet));
  }
  else {

//AliInfo(Form("The element %d in this detector %d has not enough statistic to be fitted",idect-(fCount-(fCalibraMode->GetNfragZ(0)*fCalibraMode->GetNfragRphi(0))),fCountDet));
    
    // Calcul the coef from the database choosen
    CalculChargeCoefMean(kFALSE);

    //stack 2, not stack 2
    Int_t factor = 0;
    if(GetStack(fCountDet) == 2) factor = 12;
    else factor = 16;
    
    // Fill the fCurrentCoefDetector with negative value to say: not fitted
    for (Int_t k = fCalibraMode->GetRowMin(0); k < fCalibraMode->GetRowMax(0); k++) {
      for (Int_t j = fCalibraMode->GetColMin(0); j < fCalibraMode->GetColMax(0); j++) {
	fCurrentCoefDetector[(Int_t)(j*factor+k)] = -TMath::Abs(fCurrentCoef[1]);
      }
    }
    
    //Put default value negative
    fCurrentCoef[0] = -TMath::Abs(fCurrentCoef[1]);
    fCurrentCoefE   = 0.0;
   
    FillFillCH(idect);
  }
  
  return kTRUE;
}


//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::NotEnoughStatisticPH(Int_t idect,Double_t nentries)
{
  //
  // For the case where there are not enough entries in the histograms
  // of the calibration group, the value present in the choosen database
  // will be put. A negativ sign enables to know that a fit was not possible.
  //
  if (fDebugLevel == 1) {
    AliInfo("The element has not enough statistic to be fitted");
  }
  else if (fNbDet > 0) {

    Int_t firstdetector = fCountDet;
    Int_t lastdetector  = fCountDet+fNbDet;
//AliInfo(Form("The element %d containing the detectors %d to %d has not enough statistic to be fitted",idect,firstdetector,lastdetector));
    // loop over detectors
    for(Int_t det = firstdetector; det < lastdetector; det++){

      //Set the calibration object again
      fCountDet = det;
      SetCalROC(1);   

      // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
      // Put them at 1
      fCalibraMode->ModePadCalibration((Int_t) GetStack(fCountDet),1);
      fCalibraMode->ModePadFragmentation((Int_t) GetLayer(fCountDet)
					 ,(Int_t) GetStack(fCountDet)
					 ,(Int_t) GetSector(fCountDet),1);
      // Reconstruct row min row max
      ReconstructFitRowMinRowMax(idect,1);      

      // Calcul the coef from the database choosen for the detector
      CalculVdriftCoefMean();
      CalculT0CoefMean();
      
      //stack 2, not stack 2
      Int_t factor = 0;
      if(GetStack(fCountDet) == 2) factor = 12;
      else factor = 16;
      
      // Fill the fCurrentCoefDetector with negative value to say: not fitted
      for (Int_t k = fCalibraMode->GetRowMin(1); k < fCalibraMode->GetRowMax(1); k++) {
	for (Int_t j = fCalibraMode->GetColMin(1); j < fCalibraMode->GetColMax(1); j++) {
	  fCurrentCoefDetector[(Int_t)(j*factor+k)] = -TMath::Abs(fCurrentCoef[1]);
	  fCurrentCoefDetector2[(Int_t)(j*factor+k)] = fCurrentCoef2[1] + 100.0;
      	}
      }
      
      //Put default value negative
      fCurrentCoef[0]  = -TMath::Abs(fCurrentCoef[1]);
      fCurrentCoefE    = 0.0;
      fCurrentCoef2[0] = fCurrentCoef2[1] + 100.0;
      fCurrentCoefE2   = 0.0;
            
      // Fill the stuff
      FillVectorFit();
      FillVectorFit2();
      // Debug
      if(fDebugLevel > 1){ 

	if ( !fDebugStreamer ) {
	  //debug stream
	  TDirectory *backup = gDirectory;
	  fDebugStreamer = new TTreeSRedirector("TRDDebugFitPH.root");
	  if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
	} 
	
	
	Int_t   detector     = fCountDet;
	Int_t   caligroup    = idect;
	Short_t rowmin       = fCalibraMode->GetRowMin(1);
	Short_t rowmax       = fCalibraMode->GetRowMax(1);
	Short_t colmin       = fCalibraMode->GetColMin(1);
	Short_t colmax       = fCalibraMode->GetColMax(1);
	Float_t vf           = fCurrentCoef[0]; 
	Float_t vs           = fCurrentCoef[1]; 
	Float_t vfE          = fCurrentCoefE;
	Float_t t0f          = fCurrentCoef2[0]; 
	Float_t t0s          = fCurrentCoef2[1]; 
	Float_t t0E          = fCurrentCoefE2;
	
	
	
	(* fDebugStreamer) << "FillFillPH"<<
	"detector="<<detector<<
	  "nentries="<<nentries<<
	  "caligroup="<<caligroup<<
	  "rowmin="<<rowmin<<
	  "rowmax="<<rowmax<<
	  "colmin="<<colmin<<
	  "colmax="<<colmax<<
	  "vf="<<vf<<
	  "vs="<<vs<<
	  "vfE="<<vfE<<
	  "t0f="<<t0f<<
	  "t0s="<<t0s<<
	  "t0E="<<t0E<<
	  "\n";  
      }
      // Reset
      for (Int_t k = 0; k < 2304; k++) {
	fCurrentCoefDetector[k] = 0.0;
	fCurrentCoefDetector2[k] = 0.0;
      }
      
    }// loop detector
    AliDebug(2,Form("Check the count now: fCountDet %d",fCountDet));
  }    
  else {

//AliInfo(Form("The element %d in this detector %d has not enough statistic to be fitted",idect-(fCount-(fCalibraMode->GetNfragZ(1)*fCalibraMode->GetNfragRphi(1))),fCountDet));

    CalculVdriftCoefMean();
    CalculT0CoefMean();
  
    //stack 2 and not stack 2
    Int_t factor = 0;
    if(GetStack(fCountDet) == 2) factor = 12;
    else factor = 16;


    // Fill the fCurrentCoefDetector 2
    for (Int_t k = fCalibraMode->GetRowMin(1); k < fCalibraMode->GetRowMax(1); k++) {
      for (Int_t j = fCalibraMode->GetColMin(1); j < fCalibraMode->GetColMax(1); j++) {
	fCurrentCoefDetector[(Int_t)(j*factor+k)] = -TMath::Abs(fCurrentCoef[1]);
	fCurrentCoefDetector2[(Int_t)(j*factor+k)] = fCurrentCoef2[1] + 100.0;
      }
    }

    // Put the default value
    fCurrentCoef[0]  = -TMath::Abs(fCurrentCoef[1]);
    fCurrentCoefE    = 0.0;
    fCurrentCoef2[0] = fCurrentCoef2[1] + 100.0;
    fCurrentCoefE2   = 0.0;
     
    FillFillPH(idect,nentries);
    
  }
  
  return kTRUE;
  
}


//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::NotEnoughStatisticPRF(Int_t idect)
{
  //
  // For the case where there are not enough entries in the histograms
  // of the calibration group, the value present in the choosen database
  // will be put. A negativ sign enables to know that a fit was not possible.
  //
  
  if (fDebugLevel == 1) {
    AliInfo("The element has not enough statistic to be fitted");
  }
  else if (fNbDet > 0){
  
    Int_t firstdetector = fCountDet;
    Int_t lastdetector  = fCountDet+fNbDet;
//  AliInfo(Form("The element %d containing the detectors %d to %d has not enough statistic to be fitted",idect,firstdetector,lastdetector));
    
    // loop over detectors
    for(Int_t det = firstdetector; det < lastdetector; det++){

      //Set the calibration object again
      fCountDet = det;
      SetCalROC(2);   

      // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
      // Put them at 1
      fCalibraMode->ModePadCalibration((Int_t) GetStack(fCountDet),2);
      fCalibraMode->ModePadFragmentation((Int_t) GetLayer(fCountDet)
					 ,(Int_t) GetStack(fCountDet)
					 ,(Int_t) GetSector(fCountDet),2);
      // Reconstruct row min row max
      ReconstructFitRowMinRowMax(idect,2);      

      // Calcul the coef from the database choosen for the detector
      CalculPRFCoefMean();
      
      //stack 2, not stack 2
      Int_t factor = 0;
      if(GetStack(fCountDet) == 2) factor = 12;
      else factor = 16;
      
      // Fill the fCurrentCoefDetector with negative value to say: not fitted
      for (Int_t k = fCalibraMode->GetRowMin(2); k < fCalibraMode->GetRowMax(2); k++) {
	for (Int_t j = fCalibraMode->GetColMin(2); j < fCalibraMode->GetColMax(2); j++) {
	  fCurrentCoefDetector[(Int_t)(j*factor+k)] = -TMath::Abs(fCurrentCoef[1]);
	}
      }
      
      //Put default value negative
      fCurrentCoef[0] = -TMath::Abs(fCurrentCoef[1]);
      fCurrentCoefE   = 0.0;
      
      // Fill the stuff
      FillVectorFit();
      // Debug
      if(fDebugLevel > 1){
	
	if ( !fDebugStreamer ) {
	  //debug stream
	  TDirectory *backup = gDirectory;
	  fDebugStreamer = new TTreeSRedirector("TRDDebugFitPRF.root");
	  if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
	} 
	
	Int_t   detector     = fCountDet;
	Int_t   layer        = GetLayer(fCountDet);
	Int_t   caligroup    = idect;
	Short_t rowmin       = fCalibraMode->GetRowMin(2);
	Short_t rowmax       = fCalibraMode->GetRowMax(2);
	Short_t colmin       = fCalibraMode->GetColMin(2);
	Short_t colmax       = fCalibraMode->GetColMax(2);
	Float_t widf         = fCurrentCoef[0]; 
	Float_t wids         = fCurrentCoef[1]; 
	Float_t widfE        = fCurrentCoefE;
	
	(* fDebugStreamer) << "FillFillPRF"<<
	  "detector="<<detector<<
	  "layer="<<layer<<
	  "caligroup="<<caligroup<<
	  "rowmin="<<rowmin<<
	  "rowmax="<<rowmax<<
	  "colmin="<<colmin<<
	  "colmax="<<colmax<<
	  "widf="<<widf<<
	  "wids="<<wids<<
	  "widfE="<<widfE<<
	  "\n";  
      }
      // Reset
      for (Int_t k = 0; k < 2304; k++) {
	fCurrentCoefDetector[k] = 0.0;
      }
      
    }// loop detector
    AliDebug(2,Form("Check the count now: fCountDet %d",fCountDet));
  }
  else {
    
//  AliInfo(Form("The element %d in this detector %d has not enough statistic to be fitted",idect-(fCount-(fCalibraMode->GetNfragZ(2)*fCalibraMode->GetNfragRphi(2))),fCountDet));
    
    CalculPRFCoefMean();
    
    // stack 2 and not stack 2
    Int_t factor = 0;
    if(GetStack(fCountDet) == 2) factor = 12;
    else factor = 16;

    
    // Fill the fCurrentCoefDetector
    for (Int_t k = fCalibraMode->GetRowMin(2); k < fCalibraMode->GetRowMax(2); k++) {
      for (Int_t j = fCalibraMode->GetColMin(2); j < fCalibraMode->GetColMax(2); j++) {
	fCurrentCoefDetector[(Int_t)(j*factor+k)] = -TMath::Abs(fCurrentCoef[1]);
      }
    }

    // Put the default value
    fCurrentCoef[0] = -TMath::Abs(fCurrentCoef[1]);
    fCurrentCoefE   = 0.0;
    
    FillFillPRF(idect);
  }
  
  return kTRUE;

}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::NotEnoughStatisticLinearFitter()
{
  //
  // For the case where there are not enough entries in the histograms
  // of the calibration group, the value present in the choosen database
  // will be put. A negativ sign enables to know that a fit was not possible.
  //
  
  // Calcul the coef from the database choosen
  CalculVdriftLorentzCoef();

  Int_t factor = 0;
  if(GetStack(fCountDet) == 2) factor = 1728;
  else factor = 2304;
    
    
  // Fill the fCurrentCoefDetector
  for (Int_t k = 0; k < factor; k++) {
    fCurrentCoefDetector[k] = -TMath::Abs(fCurrentCoef[1]);
    // should be negative
    fCurrentCoefDetector2[k] = fCurrentCoef2[1]+100.0;
  }
   
  
  //Put default opposite sign only for vdrift
  fCurrentCoef[0]  = -TMath::Abs(fCurrentCoef[1]);
  fCurrentCoefE    = 0.0;
  fCurrentCoef2[0] = fCurrentCoef2[1]+100.0;
  fCurrentCoefE2 = 0.0; 
  
  FillFillLinearFitter();
    
  return kTRUE;
}

//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::FillInfosFitCH(Int_t idect)
{
  //
  // Fill the coefficients found with the fits or other
  // methods from the Fit functions
  //

  if (fDebugLevel != 1) {
    if (fNbDet > 0){
      Int_t firstdetector = fCountDet;
      Int_t lastdetector  = fCountDet+fNbDet;
      //    AliInfo(Form("The element %d containing the detectors %d to %d has been fitted",idect,firstdetector,lastdetector));
      // loop over detectors
      for(Int_t det = firstdetector; det < lastdetector; det++){
	
	//Set the calibration object again
	fCountDet = det;
	SetCalROC(0);   
	
	// Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
	// Put them at 1
	fCalibraMode->ModePadCalibration((Int_t) GetStack(fCountDet),0);
	fCalibraMode->ModePadFragmentation((Int_t) GetLayer(fCountDet)
					   ,(Int_t) GetStack(fCountDet)
					   ,(Int_t) GetSector(fCountDet),0);
	// Reconstruct row min row max
	ReconstructFitRowMinRowMax(idect,0);      
	
	// Calcul the coef from the database choosen for the detector
	if(fCurrentCoef[0] < 0.0) CalculChargeCoefMean(kFALSE);
	else CalculChargeCoefMean(kTRUE);
	
	//stack 2, not stack 2
	Int_t factor = 0;
	if(GetStack(fCountDet) == 2) factor = 12;
	else factor = 16;
	
	// Fill the fCurrentCoefDetector with negative value to say: not fitted
	Double_t coeftoput = 1.0;
	if(fCurrentCoef[0] < 0.0) coeftoput = - TMath::Abs(fCurrentCoef[1]);
	else coeftoput = fCurrentCoef[0];
	for (Int_t k = fCalibraMode->GetRowMin(0); k < fCalibraMode->GetRowMax(0); k++) {
	  for (Int_t j = fCalibraMode->GetColMin(0); j < fCalibraMode->GetColMax(0); j++) {
	    fCurrentCoefDetector[(Int_t)(j*factor+k)] = coeftoput;
	  }
	}
	
	// Fill the stuff
	FillVectorFit();
	// Debug
	if(fDebugLevel > 1){ 
	  
	  if ( !fDebugStreamer ) {
	    //debug stream
	    TDirectory *backup = gDirectory;
	    fDebugStreamer = new TTreeSRedirector("TRDDebugFitCH.root");
	    if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
	  } 
	  
	  Int_t   detector   = fCountDet;
	  Int_t   caligroup  = idect;
	  Short_t rowmin     = fCalibraMode->GetRowMin(0);
	  Short_t rowmax     = fCalibraMode->GetRowMax(0);
	  Short_t colmin     = fCalibraMode->GetColMin(0);
	  Short_t colmax     = fCalibraMode->GetColMax(0);
	  Float_t gf         = fCurrentCoef[0]; 
	  Float_t gfs        = fCurrentCoef[1]; 
	  Float_t gfE        = fCurrentCoefE;
	  
	  (*fDebugStreamer) << "FillFillCH" <<
	    "detector=" << detector <<
	    "caligroup=" << caligroup <<
	    "rowmin=" << rowmin <<
	    "rowmax=" << rowmax <<
	    "colmin=" << colmin <<
	    "colmax=" << colmax <<
	    "gf=" << gf <<
	    "gfs=" << gfs <<
	    "gfE=" << gfE <<
	    "\n"; 
	  
	}
	// Reset
	for (Int_t k = 0; k < 2304; k++) {
	  fCurrentCoefDetector[k] = 0.0;
	}
	
      }// loop detector
      //printf("Check the count now: fCountDet %d\n",fCountDet);
    }
    else{
      
      Int_t factor = 0;
      if(GetStack(fCountDet) == 2) factor = 12;
      else factor = 16; 
      
      for (Int_t k = fCalibraMode->GetRowMin(0); k < fCalibraMode->GetRowMax(0); k++) {
	for (Int_t j = fCalibraMode->GetColMin(0); j < fCalibraMode->GetColMax(0); j++) {
	  fCurrentCoefDetector[(Int_t)(j*factor+k)] = fCurrentCoef[0];
	}
      }
      
      FillFillCH(idect);
    }
  }

  return kTRUE;

}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::FillInfosFitPH(Int_t idect,Double_t nentries)
{
  //
  // Fill the coefficients found with the fits or other
  // methods from the Fit functions
  //

  if (fDebugLevel != 1) {
    if (fNbDet > 0){
      
      Int_t firstdetector = fCountDet;
      Int_t lastdetector  = fCountDet+fNbDet;
// AliInfo(Form("The element %d containing the detectors %d to %d has been fitted",idect,firstdetector,lastdetector));
      
      // loop over detectors
      for(Int_t det = firstdetector; det < lastdetector; det++){
	
	//Set the calibration object again
	fCountDet = det;
	SetCalROC(1);   
	
	// Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
	// Put them at 1
	fCalibraMode->ModePadCalibration((Int_t) GetStack(fCountDet),1);
	fCalibraMode->ModePadFragmentation((Int_t) GetLayer(fCountDet)
					   ,(Int_t) GetStack(fCountDet)
					   ,(Int_t) GetSector(fCountDet),1);
	// Reconstruct row min row max
	ReconstructFitRowMinRowMax(idect,1);      
	
	// Calcul the coef from the database choosen for the detector
	CalculVdriftCoefMean();
	CalculT0CoefMean();
	 	
	//stack 2, not stack 2
	Int_t factor = 0;
	if(GetStack(fCountDet) == 2) factor = 12;
	else factor = 16;
	
	// Fill the fCurrentCoefDetector with negative value to say: not fitted
	Double_t coeftoput  = 1.5;
	Double_t coeftoput2 = 0.0; 

	if(fCurrentCoef[0] < 0.0) coeftoput = - TMath::Abs(fCurrentCoef[1]);
	else coeftoput = fCurrentCoef[0];

	if(fCurrentCoef2[0] > 70.0) coeftoput2 = fCurrentCoef2[1] + 100.0;
	else coeftoput2 = fCurrentCoef2[0];

	for (Int_t k = fCalibraMode->GetRowMin(1); k < fCalibraMode->GetRowMax(1); k++) {
	  for (Int_t j = fCalibraMode->GetColMin(1); j < fCalibraMode->GetColMax(1); j++) {
	    fCurrentCoefDetector[(Int_t)(j*factor+k)]  = coeftoput;
	    fCurrentCoefDetector2[(Int_t)(j*factor+k)] = coeftoput2;
	  }
	}
	
	// Fill the stuff
	FillVectorFit();
	FillVectorFit2();
	// Debug
	if(fDebugLevel > 1){ 
	  
	  if ( !fDebugStreamer ) {
	    //debug stream
	    TDirectory *backup = gDirectory;
	    fDebugStreamer = new TTreeSRedirector("TRDDebugFitPH.root");
	    if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
	  } 
	  
	  
	  Int_t   detector     = fCountDet;
	  Int_t   caligroup    = idect;
	  Short_t rowmin       = fCalibraMode->GetRowMin(1);
	  Short_t rowmax       = fCalibraMode->GetRowMax(1);
	  Short_t colmin       = fCalibraMode->GetColMin(1);
	  Short_t colmax       = fCalibraMode->GetColMax(1);
	  Float_t vf           = fCurrentCoef[0]; 
	  Float_t vs           = fCurrentCoef[1]; 
	  Float_t vfE          = fCurrentCoefE;
	  Float_t t0f          = fCurrentCoef2[0]; 
	  Float_t t0s          = fCurrentCoef2[1]; 
	  Float_t t0E          = fCurrentCoefE2;
	  
	  
	  
	  (* fDebugStreamer) << "FillFillPH"<<
	    "detector="<<detector<<
	    "nentries="<<nentries<<
	    "caligroup="<<caligroup<<
	    "rowmin="<<rowmin<<
	    "rowmax="<<rowmax<<
	    "colmin="<<colmin<<
	    "colmax="<<colmax<<
	    "vf="<<vf<<
	    "vs="<<vs<<
	    "vfE="<<vfE<<
	    "t0f="<<t0f<<
	    "t0s="<<t0s<<
	    "t0E="<<t0E<<
	    "\n";  
	}
	// Reset
	for (Int_t k = 0; k < 2304; k++) {
	  fCurrentCoefDetector[k] = 0.0;
	  fCurrentCoefDetector2[k] = 0.0;
	}
	
      }// loop detector
      //printf("Check the count now: fCountDet %d\n",fCountDet);
    }
    else {
      
      Int_t factor = 0;
      if(GetStack(fCountDet) == 2) factor = 12;
      else factor = 16; 
      
      for (Int_t k = fCalibraMode->GetRowMin(1); k < fCalibraMode->GetRowMax(1); k++) {
	for (Int_t j = fCalibraMode->GetColMin(1); j < fCalibraMode->GetColMax(1); j++) {
	  fCurrentCoefDetector[(Int_t)(j*factor+k)]  = fCurrentCoef[0];
	  fCurrentCoefDetector2[(Int_t)(j*factor+k)] = fCurrentCoef2[0];
	}
      }  
      
      FillFillPH(idect,nentries);
    }
  }
  return kTRUE;
}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::FillInfosFitPRF(Int_t idect)
{
  //
  // Fill the coefficients found with the fits or other
  // methods from the Fit functions
  //
  
  if (fDebugLevel != 1) {
    if (fNbDet > 0){
    
      Int_t firstdetector = fCountDet;
      Int_t lastdetector  = fCountDet+fNbDet;
//    AliInfo(Form("The element %d containing the detectors %d to %d has been fitted",idect,firstdetector,lastdetector));
      
      // loop over detectors
      for(Int_t det = firstdetector; det < lastdetector; det++){
	
	//Set the calibration object again
	fCountDet = det;
	SetCalROC(2);   
	
	// Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
	// Put them at 1
	fCalibraMode->ModePadCalibration((Int_t) GetStack(fCountDet),2);
	fCalibraMode->ModePadFragmentation((Int_t) GetLayer(fCountDet)
					   ,(Int_t) GetStack(fCountDet)
					   ,(Int_t) GetSector(fCountDet),2);
	// Reconstruct row min row max
	ReconstructFitRowMinRowMax(idect,2);      
	
	// Calcul the coef from the database choosen for the detector
	CalculPRFCoefMean();
	
	//stack 2, not stack 2
	Int_t factor = 0;
	if(GetStack(fCountDet) == 2) factor = 12;
	else factor = 16;
	
	// Fill the fCurrentCoefDetector with negative value to say: not fitted
	Double_t coeftoput = 1.0;
	if(fCurrentCoef[0] < 0.0) coeftoput = - TMath::Abs(fCurrentCoef[1]);
	else coeftoput = fCurrentCoef[0];
	for (Int_t k = fCalibraMode->GetRowMin(2); k < fCalibraMode->GetRowMax(2); k++) {
	  for (Int_t j = fCalibraMode->GetColMin(2); j < fCalibraMode->GetColMax(2); j++) {
	    fCurrentCoefDetector[(Int_t)(j*factor+k)] = coeftoput;
	  }
	}
	
	// Fill the stuff
	FillVectorFit();
	// Debug
	if(fDebugLevel > 1){
	  
	  if ( !fDebugStreamer ) {
	    //debug stream
	    TDirectory *backup = gDirectory;
	    fDebugStreamer = new TTreeSRedirector("TRDDebugFitPRF.root");
	    if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
	  } 
	  
	  Int_t   detector     = fCountDet;
	  Int_t   layer        = GetLayer(fCountDet);
	  Int_t   caligroup    = idect;
	  Short_t rowmin       = fCalibraMode->GetRowMin(2);
	  Short_t rowmax       = fCalibraMode->GetRowMax(2);
	  Short_t colmin       = fCalibraMode->GetColMin(2);
	  Short_t colmax       = fCalibraMode->GetColMax(2);
	  Float_t widf         = fCurrentCoef[0]; 
	  Float_t wids         = fCurrentCoef[1]; 
	  Float_t widfE        = fCurrentCoefE;
	  
	  (* fDebugStreamer) << "FillFillPRF"<<
	    "detector="<<detector<<
	    "layer="<<layer<<
	    "caligroup="<<caligroup<<
	    "rowmin="<<rowmin<<
	    "rowmax="<<rowmax<<
	    "colmin="<<colmin<<
	    "colmax="<<colmax<<
	    "widf="<<widf<<
	    "wids="<<wids<<
	    "widfE="<<widfE<<
	    "\n";  
	}
	// Reset
	for (Int_t k = 0; k < 2304; k++) {
	  fCurrentCoefDetector[k] = 0.0;
	}
	
      }// loop detector
      //printf("Check the count now: fCountDet %d\n",fCountDet);
    }
    else {
      
      Int_t factor = 0;
      if(GetStack(fCountDet) == 2) factor = 12;
      else factor = 16; 
      
      // Pointer to the branch
      for (Int_t k = fCalibraMode->GetRowMin(2); k < fCalibraMode->GetRowMax(2); k++) {
	for (Int_t j = fCalibraMode->GetColMin(2); j < fCalibraMode->GetColMax(2); j++) {
	  fCurrentCoefDetector[(Int_t)(j*factor+k)] = fCurrentCoef[0];
	}
      }
      FillFillPRF(idect);   
    }
  }
  
  return kTRUE;

}
//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::FillInfosFitLinearFitter()
{
  //
  // Fill the coefficients found with the fits or other
  // methods from the Fit functions
  //
  
  Int_t factor = 0;
  if(GetStack(fCountDet) == 2) factor = 1728;
  else factor = 2304; 
  
  // Pointer to the branch
  for (Int_t k = 0; k < factor; k++) {
    fCurrentCoefDetector[k]  = fCurrentCoef[0];
    fCurrentCoefDetector2[k] = fCurrentCoef2[0];
  }
  
  FillFillLinearFitter();
  
  return kTRUE;

}
//________________________________________________________________________________
void AliTRDCalibraFit::FillFillCH(Int_t idect)
{
  //
  // DebugStream and fVectorFit
  //

  // End of one detector
  if ((idect == (fCount-1))) {
    FillVectorFit();
    // Reset
    for (Int_t k = 0; k < 2304; k++) {
      fCurrentCoefDetector[k] = 0.0;
    }
  }

  if(fDebugLevel > 1){ 

    if ( !fDebugStreamer ) {
      //debug stream
      TDirectory *backup = gDirectory;
      fDebugStreamer = new TTreeSRedirector("TRDDebugFitCH.root");
      if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
    } 
    
    Int_t   detector   = fCountDet;
    Int_t   caligroup  = idect;
    Short_t rowmin     = fCalibraMode->GetRowMin(0);
    Short_t rowmax     = fCalibraMode->GetRowMax(0);
    Short_t colmin     = fCalibraMode->GetColMin(0);
    Short_t colmax     = fCalibraMode->GetColMax(0);
    Float_t gf         = fCurrentCoef[0]; 
    Float_t gfs        = fCurrentCoef[1]; 
    Float_t gfE        = fCurrentCoefE;
    
    (*fDebugStreamer) << "FillFillCH" <<
      "detector=" << detector <<
      "caligroup=" << caligroup <<
      "rowmin=" << rowmin <<
      "rowmax=" << rowmax <<
      "colmin=" << colmin <<
      "colmax=" << colmax <<
      "gf=" << gf <<
      "gfs=" << gfs <<
      "gfE=" << gfE <<
      "\n"; 
    
  }
}
//________________________________________________________________________________
void AliTRDCalibraFit::FillFillPH(Int_t idect,Double_t nentries)
{
  //
  // DebugStream and fVectorFit and fVectorFit2
  //
  
  // End of one detector
    if ((idect == (fCount-1))) {
      FillVectorFit();
      FillVectorFit2();
      // Reset
      for (Int_t k = 0; k < 2304; k++) {
	fCurrentCoefDetector[k] = 0.0;
	fCurrentCoefDetector2[k] = 0.0;
      }
    }

    if(fDebugLevel > 1){ 

      if ( !fDebugStreamer ) {
	//debug stream
	TDirectory *backup = gDirectory;
	fDebugStreamer = new TTreeSRedirector("TRDDebugFitPH.root");
	if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
      } 
      
       
      Int_t   detector     = fCountDet;
      Int_t   caligroup    = idect;
      Short_t rowmin       = fCalibraMode->GetRowMin(1);
      Short_t rowmax       = fCalibraMode->GetRowMax(1);
      Short_t colmin       = fCalibraMode->GetColMin(1);
      Short_t colmax       = fCalibraMode->GetColMax(1);
      Float_t vf           = fCurrentCoef[0]; 
      Float_t vs           = fCurrentCoef[1]; 
      Float_t vfE          = fCurrentCoefE;
      Float_t t0f          = fCurrentCoef2[0]; 
      Float_t t0s          = fCurrentCoef2[1]; 
      Float_t t0E          = fCurrentCoefE2;
   


      (* fDebugStreamer) << "FillFillPH"<<
	"detector="<<detector<<
	"nentries="<<nentries<<
	"caligroup="<<caligroup<<
	"rowmin="<<rowmin<<
	"rowmax="<<rowmax<<
	"colmin="<<colmin<<
	"colmax="<<colmax<<
	"vf="<<vf<<
	"vs="<<vs<<
	"vfE="<<vfE<<
	"t0f="<<t0f<<
	"t0s="<<t0s<<
	"t0E="<<t0E<<
	"\n";  
    }

}
//________________________________________________________________________________
void AliTRDCalibraFit::FillFillPRF(Int_t idect)
{
  //
  // DebugStream and fVectorFit
  //

    // End of one detector
    if ((idect == (fCount-1))) {
      FillVectorFit();
      // Reset
      for (Int_t k = 0; k < 2304; k++) {
	fCurrentCoefDetector[k] = 0.0;
      }
    }

    
    if(fDebugLevel > 1){

      if ( !fDebugStreamer ) {
	//debug stream
	TDirectory *backup = gDirectory;
	fDebugStreamer = new TTreeSRedirector("TRDDebugFitPRF.root");
	if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
      } 
      
      Int_t   detector     = fCountDet;
      Int_t   layer        = GetLayer(fCountDet);
      Int_t   caligroup    = idect;
      Short_t rowmin       = fCalibraMode->GetRowMin(2);
      Short_t rowmax       = fCalibraMode->GetRowMax(2);
      Short_t colmin       = fCalibraMode->GetColMin(2);
      Short_t colmax       = fCalibraMode->GetColMax(2);
      Float_t widf         = fCurrentCoef[0]; 
      Float_t wids         = fCurrentCoef[1]; 
      Float_t widfE        = fCurrentCoefE;

      (* fDebugStreamer) << "FillFillPRF"<<
	"detector="<<detector<<
	"layer="<<layer<<
	"caligroup="<<caligroup<<
	"rowmin="<<rowmin<<
	"rowmax="<<rowmax<<
	"colmin="<<colmin<<
	"colmax="<<colmax<<
	"widf="<<widf<<
	"wids="<<wids<<
	"widfE="<<widfE<<
	"\n";  
    }

}
//________________________________________________________________________________
void AliTRDCalibraFit::FillFillLinearFitter()
{
  //
  // DebugStream and fVectorFit
  //

  // End of one detector
  FillVectorFit();
  FillVectorFit2();
  
  
  // Reset
  for (Int_t k = 0; k < 2304; k++) {
  fCurrentCoefDetector[k]  = 0.0;
  fCurrentCoefDetector2[k] = 0.0;
  }
  

  if(fDebugLevel > 1){

    if ( !fDebugStreamer ) {
      //debug stream
      TDirectory *backup = gDirectory;
      fDebugStreamer = new TTreeSRedirector("TRDDebugFitLinearFitter.root");
      if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
    } 
    
    //Debug: comparaison of the different methods (okey for first time but not for iterative procedure)
    AliTRDpadPlane *padplane = fGeo->GetPadPlane(GetLayer(fCountDet),GetStack(fCountDet));
    Float_t rowmd            = (padplane->GetRow0()+padplane->GetRowEnd())/2.;
    Float_t r                = AliTRDgeometry::GetTime0(GetLayer(fCountDet)); 
    Float_t tiltangle        = padplane->GetTiltingAngle();
    Int_t   detector         = fCountDet;
    Int_t   stack            = GetStack(fCountDet);
    Int_t   layer            = GetLayer(fCountDet);
    Float_t vf               = fCurrentCoef[0]; 
    Float_t vs               = fCurrentCoef[1]; 
    Float_t vfE              = fCurrentCoefE;
    Float_t lorentzangler    = fCurrentCoef2[0];
    Float_t elorentzangler   = fCurrentCoefE2;
    Float_t lorentzangles    = fCurrentCoef2[1];
   
    (* fDebugStreamer) << "FillFillLinearFitter"<<
      "detector="<<detector<<
      "stack="<<stack<<
      "layer="<<layer<<
      "rowmd="<<rowmd<<
      "r="<<r<<
      "tiltangle="<<tiltangle<<
      "vf="<<vf<<
      "vs="<<vs<<
      "vfE="<<vfE<<
      "lorentzangler="<<lorentzangler<<
      "Elorentzangler="<<elorentzangler<<
      "lorentzangles="<<lorentzangles<<
      "\n";  
  }
  
}
//
//____________Calcul Coef Mean_________________________________________________
//
//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::CalculT0CoefMean()
{
  //
  // For the detector Dect calcul the mean time 0
  // for the calibration group idect from the choosen database
  //

  fCurrentCoef2[1] = 0.0;
  if(fDebugLevel != 1){
    if(((fCalibraMode->GetNz(1) > 0) ||
       (fCalibraMode->GetNrphi(1) > 0)) && ((fCalibraMode->GetNz(1)    != 10) && (fCalibraMode->GetNz(1)    != 100))) {

      for (Int_t row = fCalibraMode->GetRowMin(1); row < fCalibraMode->GetRowMax(1); row++) {
	for (Int_t col = fCalibraMode->GetColMin(1); col < fCalibraMode->GetColMax(1); col++) {
	  fCurrentCoef2[1] += (Float_t) (fCalROC2->GetValue(col,row)+fCalDet2->GetValue(fCountDet));
	}
      }
      
      fCurrentCoef2[1] = fCurrentCoef2[1] / ((fCalibraMode->GetColMax(1)-fCalibraMode->GetColMin(1))*(fCalibraMode->GetRowMax(1)-fCalibraMode->GetRowMin(1)));
    
    }
    else {
     
      if(!fAccCDB){
	fCurrentCoef2[1] = fCalDet2->GetValue(fCountDet);
      }
      else{
	
	for(Int_t row = 0; row < fGeo->GetRowMax(GetLayer(fCountDet),GetStack(fCountDet),GetSector(fCountDet)); row++){
	  for(Int_t col = 0; col < fGeo->GetColMax(GetLayer(fCountDet)); col++){
	    fCurrentCoef2[1] += (Float_t) (fCalROC2->GetValue(col,row)+fCalDet2->GetValue(fCountDet));
	  }
	}
	fCurrentCoef2[1] = fCurrentCoef2[1] / ((fGeo->GetRowMax(GetLayer(fCountDet),GetStack(fCountDet),GetSector(fCountDet)))*(fGeo->GetColMax(GetLayer(fCountDet))));
      
      }
    }
  }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::CalculChargeCoefMean(Bool_t vrai)
{
  //
  // For the detector Dect calcul the mean gain factor
  // for the calibration group idect from the choosen database
  //

  fCurrentCoef[1] = 0.0;
  if(fDebugLevel != 1){
    if (((fCalibraMode->GetNz(0)    > 0) || 
	(fCalibraMode->GetNrphi(0) > 0)) && ((fCalibraMode->GetNz(0)    != 10) && (fCalibraMode->GetNz(0)    != 100))) {
      for (Int_t row = fCalibraMode->GetRowMin(0); row < fCalibraMode->GetRowMax(0); row++) {
	for (Int_t col = fCalibraMode->GetColMin(0); col < fCalibraMode->GetColMax(0); col++) {
	  fCurrentCoef[1] += (Float_t) (fCalROC->GetValue(col,row)*fCalDet->GetValue(fCountDet));
	  if (vrai) fScaleFitFactor += (Float_t) (fCalROC->GetValue(col,row)*fCalDet->GetValue(fCountDet));
	}
      }
      fCurrentCoef[1] = fCurrentCoef[1] / ((fCalibraMode->GetColMax(0)-fCalibraMode->GetColMin(0))*(fCalibraMode->GetRowMax(0)-fCalibraMode->GetRowMin(0)));
    }
    else {
      //Per detectors
      fCurrentCoef[1] = (Float_t) fCalDet->GetValue(fCountDet);
      if (vrai) fScaleFitFactor += ((Float_t) fCalDet->GetValue(fCountDet))*(fCalibraMode->GetColMax(0)-fCalibraMode->GetColMin(0))*(fCalibraMode->GetRowMax(0)-fCalibraMode->GetRowMin(0));
    }    
  }
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::CalculPRFCoefMean()
{
  //
  // For the detector Dect calcul the mean sigma of pad response
  // function for the calibration group idect from the choosen database
  //
  
  fCurrentCoef[1] = 0.0;
  if(fDebugLevel != 1){
    for (Int_t row = fCalibraMode->GetRowMin(2); row < fCalibraMode->GetRowMax(2); row++) {
      for (Int_t col = fCalibraMode->GetColMin(2); col < fCalibraMode->GetColMax(2); col++) {
	fCurrentCoef[1] += (Float_t) fCalROC->GetValue(col,row);
      }
    }
    fCurrentCoef[1] = fCurrentCoef[1] / ((fCalibraMode->GetColMax(2)-fCalibraMode->GetColMin(2))*(fCalibraMode->GetRowMax(2)-fCalibraMode->GetRowMin(2)));
  }
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::CalculVdriftCoefMean()
{
  //
  // For the detector dect calcul the mean drift velocity for the
  // calibration group idect from the choosen database
  //

  fCurrentCoef[1] = 0.0;
  if(fDebugLevel != 1){
    if (((fCalibraMode->GetNz(1)    > 0) || 
	(fCalibraMode->GetNrphi(1) > 0)) && ((fCalibraMode->GetNz(1)    != 10) && (fCalibraMode->GetNz(1)    != 100))) {
      
      for (Int_t row = fCalibraMode->GetRowMin(1); row < fCalibraMode->GetRowMax(1); row++) {
	for (Int_t col = fCalibraMode->GetColMin(1); col < fCalibraMode->GetColMax(1); col++) {
	  fCurrentCoef[1] += (Float_t) (fCalROC->GetValue(col,row)*fCalDet->GetValue(fCountDet));
	}
      }
      
      fCurrentCoef[1] = fCurrentCoef[1] / ((fCalibraMode->GetColMax(1)-fCalibraMode->GetColMin(1))*(fCalibraMode->GetRowMax(1)-fCalibraMode->GetRowMin(1)));
    
    }
    else {
      //per detectors
      fCurrentCoef[1] = (Float_t) fCalDet->GetValue(fCountDet);
    }  
  }
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::CalculVdriftLorentzCoef()
{
  //
  // For the detector fCountDet, mean drift velocity and tan lorentzangle
  //

  fCurrentCoef[1]  = fCalDetVdriftUsed->GetValue(fCountDet);
  fCurrentCoef2[1] = fCalDetExBUsed->GetValue(fCountDet); 

  return kTRUE;
}
//_____________________________________________________________________________
Float_t AliTRDCalibraFit::GetPRFDefault(Int_t layer) const
{
  //
  // Default width of the PRF if there is no database as reference
  //
  switch(layer)
    {
      // default database
      //case 0:  return 0.515;
      //case 1:  return 0.502;
      //case 2:  return 0.491;
      //case 3:  return 0.481;
      //case 4:  return 0.471;
      //case 5:  return 0.463;
      //default: return 0.0;

      // fit database
    case 0:  return 0.538429;
    case 1:  return 0.524302;
    case 2:  return 0.511591;
    case 3:  return 0.500140;
    case 4:  return 0.489821;
    case 5:  return 0.480524;
    default: return 0.0;
  }
}
//________________________________________________________________________________
void AliTRDCalibraFit::SetCalROC(Int_t i)
{
  //
  // Set the calib object for fCountDet
  //

  Float_t value = 0.0;
  
  //Get the CalDet object
  if(fAccCDB){
    AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
    if (!cal) {
      AliInfo("Could not get calibDB");
      return;
    }
    switch (i)
      {
      case 0: 
	if( fCalROC ){ 
	  fCalROC->~AliTRDCalROC();
	  new(fCalROC) AliTRDCalROC(*(cal->GetGainFactorROC(fCountDet)));
	}else fCalROC = new AliTRDCalROC(*(cal->GetGainFactorROC(fCountDet)));
	break;
      case 1:
	if( fCalROC ){ 
	  fCalROC->~AliTRDCalROC();
	  new(fCalROC) AliTRDCalROC(*(cal->GetVdriftROC(fCountDet)));
	}else fCalROC = new AliTRDCalROC(*(cal->GetVdriftROC(fCountDet)));
	if( fCalROC2 ){ 
	  fCalROC2->~AliTRDCalROC();
	  new(fCalROC2) AliTRDCalROC(*(cal->GetT0ROC(fCountDet)));
	}else fCalROC2 = new AliTRDCalROC(*(cal->GetT0ROC(fCountDet)));
	break;
      case 2:
	if( fCalROC ){ 
	  fCalROC->~AliTRDCalROC();
	  new(fCalROC) AliTRDCalROC(*(cal->GetPRFROC(fCountDet)));
	}else fCalROC = new AliTRDCalROC(*(cal->GetPRFROC(fCountDet)));
	break; 
      default: return;
      }
  }
  else{
    switch (i)
      {
      case 0:
	if(fCalROC) delete fCalROC;
	fCalROC = new AliTRDCalROC(GetLayer(fCountDet),GetStack(fCountDet)); 
	for(Int_t k = 0; k < fCalROC->GetNchannels(); k++){
	  fCalROC->SetValue(k,1.0);
	}
	break;
      case 1:
	if(fCalROC)  delete fCalROC;
	if(fCalROC2) delete fCalROC2;
	fCalROC  = new AliTRDCalROC(GetLayer(fCountDet),GetStack(fCountDet));
	fCalROC2 = new AliTRDCalROC(GetLayer(fCountDet),GetStack(fCountDet));
	for(Int_t k = 0; k < fCalROC->GetNchannels(); k++){
	  fCalROC->SetValue(k,1.0);
	  fCalROC2->SetValue(k,0.0);
	}
	break;
      case 2:
	if(fCalROC) delete fCalROC;
	value = GetPRFDefault(GetLayer(fCountDet));
	fCalROC = new AliTRDCalROC(GetLayer(fCountDet),GetStack(fCountDet)); 
	for(Int_t k = 0; k < fCalROC->GetNchannels(); k++){
	  fCalROC->SetValue(k,value);
	}
	break;
      default: return; 
      }
  }
  
}
//____________Fit Methods______________________________________________________

//_____________________________________________________________________________
void AliTRDCalibraFit::FitPente(TH1* projPH)
{
  //
  // Slope methode for the drift velocity
  //
  
  // Constants
  const Float_t kDrWidth = AliTRDgeometry::DrThick();
  Int_t binmax           = 0;
  Int_t binmin           = 0;
  fPhd[0]                = 0.0;
  fPhd[1]                = 0.0;
  fPhd[2]                = 0.0;
  Int_t ju               = 0;
  fCurrentCoefE          = 0.0;
  fCurrentCoefE2         = 0.0;
  fCurrentCoef[0]        = 0.0;
  fCurrentCoef2[0]       = 0.0;
  TLine *line            = new TLine();

  // Some variables
  TAxis   *xpph    = projPH->GetXaxis();
  Int_t    nbins   = xpph->GetNbins();
  Double_t lowedge = xpph->GetBinLowEdge(1);
  Double_t upedge  = xpph->GetBinUpEdge(xpph->GetNbins());
  Double_t widbins = (upedge - lowedge) / nbins;
  Double_t limit   = upedge + 0.5 * widbins; 
  Bool_t put = kTRUE;

  // Beginning of the signal
  TH1D *pentea = new TH1D("pentea","pentea",projPH->GetNbinsX(),0,(Float_t) limit);
  for (Int_t k = 1; k <  projPH->GetNbinsX(); k++) {
    pentea->SetBinContent(k,(Double_t) (projPH->GetBinContent(k+1) - projPH->GetBinContent(k)));
  }
  binmax = (Int_t) pentea->GetMaximumBin();
  if (binmax <= 1) {
    binmax = 2;
    AliInfo("Put the binmax from 1 to 2 to enable the fit");
  }
  if (binmax >= nbins) {
    binmax = nbins-1;
    put = kFALSE;
    AliInfo("Put the binmax from nbins-1 to nbins-2 to enable the fit");
  }
  pentea->Fit("pol2","0MR","",TMath::Max(pentea->GetBinCenter(binmax-1),0.0),pentea->GetBinCenter(binmax+1));
  Float_t l3P1am = pentea->GetFunction("pol2")->GetParameter(1);
  Float_t l3P2am = pentea->GetFunction("pol2")->GetParameter(2);
  Float_t l3P1amE = pentea->GetFunction("pol2")->GetParError(1);
  Float_t l3P2amE = pentea->GetFunction("pol2")->GetParError(2);
  if (TMath::Abs(l3P2am) > 0.00000001) {
    fPhd[0] = -(l3P1am / (2 * l3P2am));
  }
  if(!fTakeTheMaxPH){
    if((TMath::Abs(l3P1am) > 0.0000000001) && (TMath::Abs(l3P2am) > 0.00000000001)){
      fCurrentCoefE2 = (l3P1amE/l3P1am + l3P2amE/l3P2am)*fPhd[0];
    }
  }
  // Amplification region
  binmax = 0;
  ju     = 0;
  for (Int_t kbin = 1; kbin < projPH->GetNbinsX(); kbin ++) {
    if (((projPH->GetBinContent(kbin+1) - projPH->GetBinContent(kbin)) <= 0.0) && (ju == 0) && (kbin > (fPhd[0]/widbins))) {
      binmax = kbin;
      ju     = 1;
    }
  }
  if (binmax <= 1) {
    binmax = 2;
    AliInfo("Put the binmax from 1 to 2 to enable the fit");
  }
  if (binmax >= nbins) {
    binmax = nbins-1;
    put = kFALSE;
    AliInfo("Put the binmax from nbins-1 to nbins-2 to enable the fit");
  }
  projPH->Fit("pol2","0MR","",TMath::Max(projPH->GetBinCenter(binmax-1),0.0),projPH->GetBinCenter(binmax+1));
  Float_t l3P1amf = projPH->GetFunction("pol2")->GetParameter(1);
  Float_t l3P2amf = projPH->GetFunction("pol2")->GetParameter(2);
  Float_t l3P1amfE = projPH->GetFunction("pol2")->GetParError(1);
  Float_t l3P2amfE = projPH->GetFunction("pol2")->GetParError(2);
  if (TMath::Abs(l3P2amf) > 0.00000000001) {
    fPhd[1] = -(l3P1amf / (2 * l3P2amf));
  }
  if((TMath::Abs(l3P1amf) > 0.0000000001) && (TMath::Abs(l3P2amf) > 0.000000001)){
    fCurrentCoefE = (l3P1amfE/l3P1amf + l3P2amfE/l3P2amf)*fPhd[1];
  }
  if(fTakeTheMaxPH){
    fCurrentCoefE2 = fCurrentCoefE;
  }
  // Drift region
  TH1D *pente = new TH1D("pente","pente",projPH->GetNbinsX(),0,(Float_t) limit);
  for (Int_t k = TMath::Min(binmax+4,projPH->GetNbinsX()); k <  projPH->GetNbinsX(); k++) {
    pente->SetBinContent(k,(Double_t) (projPH->GetBinContent(k+1) - projPH->GetBinContent(k)));
  }
  binmin = 0;
  if(pente->GetEntries() > 0) binmin = (Int_t) pente->GetMinimumBin();
  if (binmin <= 1) {
    binmin = 2;
    AliInfo("Put the binmax from 1 to 2 to enable the fit");
  }
  if (binmin >= nbins) {
    binmin = nbins-1;
    put = kFALSE;
    AliInfo("Put the binmax from nbins-1 to nbins-2 to enable the fit");
  }
  pente->Fit("pol2"
            ,"0MR"
            ,""
            ,TMath::Max(pente->GetBinCenter(binmin-1),             0.0)
            ,TMath::Min(pente->GetBinCenter(binmin+1),(Double_t) limit));
  Float_t l3P1dr = pente->GetFunction("pol2")->GetParameter(1);
  Float_t l3P2dr = pente->GetFunction("pol2")->GetParameter(2);
  Float_t l3P1drE = pente->GetFunction("pol2")->GetParError(1);
  Float_t l3P2drE = pente->GetFunction("pol2")->GetParError(2);
  if (TMath::Abs(l3P2dr) > 0.00000001) {
    fPhd[2] = -(l3P1dr / (2 * l3P2dr));
  }
  if((TMath::Abs(l3P1dr) > 0.0000000001) && (TMath::Abs(l3P2dr) > 0.00000000001)){
    fCurrentCoefE += (l3P1drE/l3P1dr + l3P2drE/l3P2dr)*fPhd[2]; 
  }
  Float_t fPhdt0  = 0.0;
  Float_t t0Shift = 0.0;
  if(fTakeTheMaxPH) {
    fPhdt0 = fPhd[1];
    t0Shift = fT0Shift1;
  }
  else {
    fPhdt0 = fPhd[0];
    t0Shift = fT0Shift0;
  }

  if ((fPhd[2] > fPhd[0]) && 
      (fPhd[2] > fPhd[1]) && 
      (fPhd[1] > fPhd[0]) &&
      (put)) {
    fCurrentCoef[0] = (kDrWidth) / (fPhd[2]-fPhd[1]);
    fNumberFitSuccess++;

    if (fPhdt0 >= 0.0) {
      fCurrentCoef2[0] = (fPhdt0 - t0Shift) / widbins;
      if (fCurrentCoef2[0] < -3.0) {
        fCurrentCoef2[0] = fCurrentCoef2[1] + 100.0;
      }
    }
    else {
      fCurrentCoef2[0] = fCurrentCoef2[1] + 100.0;
    }

  }
  else {
    fCurrentCoef[0]  = -TMath::Abs(fCurrentCoef[1]);
    fCurrentCoef2[0] = fCurrentCoef2[1] + 100.0;
  }

  if (fDebugLevel == 1) {
    TCanvas *cpentei = new TCanvas("cpentei","cpentei",50,50,600,800);
    cpentei->cd();
    projPH->Draw();
    line->SetLineColor(2);
    line->DrawLine(fPhd[0],0,fPhd[0],projPH->GetMaximum());
    line->DrawLine(fPhd[1],0,fPhd[1],projPH->GetMaximum());
    line->DrawLine(fPhd[2],0,fPhd[2],projPH->GetMaximum());
    AliInfo(Form("fPhd[0] (beginning of the signal): %f"                  ,(Float_t) fPhd[0]));
    AliInfo(Form("fPhd[1] (end of the amplification region): %f"          ,(Float_t) fPhd[1]));
    AliInfo(Form("fPhd[2] (end of the drift region): %f"                  ,(Float_t) fPhd[2]));
    AliInfo(Form("fVriftCoef[1] (with only the drift region(default)): %f",(Float_t) fCurrentCoef[0]));
    TCanvas *cpentei2 = new TCanvas("cpentei2","cpentei2",50,50,600,800);
    cpentei2->cd();
    pentea->Draw();
    TCanvas *cpentei3 = new TCanvas("cpentei3","cpentei3",50,50,600,800);
    cpentei3->cd();
    pente->Draw();
  }
  else {
    delete pentea;
    delete pente;
  }
}
//_____________________________________________________________________________
void AliTRDCalibraFit::FitLagrangePoly(TH1* projPH)
{
  //
  // Slope methode but with polynomes de Lagrange
  //

  // Constants
  const Float_t kDrWidth = AliTRDgeometry::DrThick();
  Int_t binmax      = 0;
  Int_t binmin      = 0;
  //Double_t    *x    = new Double_t[5];
  //Double_t    *y    = new Double_t[5];
  Double_t x[5];
  Double_t y[5];
  x[0]              = 0.0;
  x[1]              = 0.0;
  x[2]              = 0.0;
  x[3]              = 0.0;
  x[4]              = 0.0;
  y[0]              = 0.0;
  y[1]              = 0.0;
  y[2]              = 0.0;
  y[3]              = 0.0;
  y[4]              = 0.0;
  fPhd[0]           = 0.0;
  fPhd[1]           = 0.0;
  fPhd[2]           = 0.0;
  Int_t ju          = 0;
  fCurrentCoefE     = 0.0;
  fCurrentCoefE2    = 1.0;
  fCurrentCoef[0]   = 0.0;
  fCurrentCoef2[0]  = 0.0;
  TLine *line = new TLine();
  TF1 * polynome = 0x0;
  TF1 * polynomea = 0x0;
  TF1 * polynomeb = 0x0;
  Double_t c0 = 0.0;
  Double_t c1 = 0.0;
  Double_t c2 = 0.0;
  Double_t c3 = 0.0;
  Double_t c4 = 0.0;
  
  // Some variables
  TAxis   *xpph    = projPH->GetXaxis();
  Int_t    nbins   = xpph->GetNbins();
  Double_t lowedge = xpph->GetBinLowEdge(1);
  Double_t upedge  = xpph->GetBinUpEdge(xpph->GetNbins());
  Double_t widbins = (upedge - lowedge) / nbins;
  Double_t limit   = upedge + 0.5 * widbins;

  
  Bool_t put = kTRUE;

  // Beginning of the signal
  TH1D *pentea = new TH1D("pentea","pentea",projPH->GetNbinsX(),0,(Float_t) limit);
  for (Int_t k = 1; k <  projPH->GetNbinsX(); k++) {
    pentea->SetBinContent(k,(Double_t) (projPH->GetBinContent(k+1) - projPH->GetBinContent(k)));
  }

  binmax = (Int_t) pentea->GetMaximumBin();

  Double_t minnn = 0.0;
  Double_t maxxx = 0.0;

  Int_t kase = nbins-binmax;
  
  switch(kase)
    {
    case 0:
      put = kFALSE;
      break;
    case 1:
      minnn = pentea->GetBinCenter(binmax-2);
      maxxx = pentea->GetBinCenter(binmax);
      x[0] = pentea->GetBinCenter(binmax-2);
      x[1] = pentea->GetBinCenter(binmax-1);
      x[2] = pentea->GetBinCenter(binmax);
      y[0] = pentea->GetBinContent(binmax-2);
      y[1] = pentea->GetBinContent(binmax-1);
      y[2] = pentea->GetBinContent(binmax);
      CalculPolynomeLagrange2(x,y,c0,c1,c2,c3,c4);
      //AliInfo("At the limit for beginning!");
      break;  
    case 2:
      minnn = pentea->GetBinCenter(binmax-2);
      maxxx = pentea->GetBinCenter(binmax+1);
      x[0] = pentea->GetBinCenter(binmax-2);
      x[1] = pentea->GetBinCenter(binmax-1);
      x[2] = pentea->GetBinCenter(binmax);
      x[3] = pentea->GetBinCenter(binmax+1);
      y[0] = pentea->GetBinContent(binmax-2);
      y[1] = pentea->GetBinContent(binmax-1);
      y[2] = pentea->GetBinContent(binmax);
      y[3] = pentea->GetBinContent(binmax+1);
      CalculPolynomeLagrange3(x,y,c0,c1,c2,c3,c4);
      break;
    default:
      switch(binmax){
      case 0:
	put = kFALSE;
	break;
      case 1:
	minnn = pentea->GetBinCenter(binmax);
	maxxx = pentea->GetBinCenter(binmax+2);
	x[0] = pentea->GetBinCenter(binmax);
	x[1] = pentea->GetBinCenter(binmax+1);
	x[2] = pentea->GetBinCenter(binmax+2);
	y[0] = pentea->GetBinContent(binmax);
	y[1] = pentea->GetBinContent(binmax+1);
	y[2] = pentea->GetBinContent(binmax+2);
	CalculPolynomeLagrange2(x,y,c0,c1,c2,c3,c4);
	break;
      case 2:
	minnn = pentea->GetBinCenter(binmax-1);
	maxxx = pentea->GetBinCenter(binmax+2);
	x[0] = pentea->GetBinCenter(binmax-1);
	x[1] = pentea->GetBinCenter(binmax);
	x[2] = pentea->GetBinCenter(binmax+1);
	x[3] = pentea->GetBinCenter(binmax+2);
	y[0] = pentea->GetBinContent(binmax-1);
	y[1] = pentea->GetBinContent(binmax);
	y[2] = pentea->GetBinContent(binmax+1);
	y[3] = pentea->GetBinContent(binmax+2);
       	CalculPolynomeLagrange3(x,y,c0,c1,c2,c3,c4);
	break;
      default:
     	minnn = pentea->GetBinCenter(binmax-2);
	maxxx = pentea->GetBinCenter(binmax+2);
	x[0] = pentea->GetBinCenter(binmax-2);
	x[1] = pentea->GetBinCenter(binmax-1);
	x[2] = pentea->GetBinCenter(binmax);
	x[3] = pentea->GetBinCenter(binmax+1);
	x[4] = pentea->GetBinCenter(binmax+2);
	y[0] = pentea->GetBinContent(binmax-2);
	y[1] = pentea->GetBinContent(binmax-1);
	y[2] = pentea->GetBinContent(binmax);
	y[3] = pentea->GetBinContent(binmax+1);
	y[4] = pentea->GetBinContent(binmax+2);
	CalculPolynomeLagrange4(x,y,c0,c1,c2,c3,c4);
	break;
      }
      break;
    }
  
  
  if(put) {
    polynomeb = new TF1("polb","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",minnn,maxxx);
    polynomeb->SetParameters(c0,c1,c2,c3,c4);
      
    Double_t step = (maxxx-minnn)/10000;
    Double_t l = minnn;
    Double_t maxvalue = 0.0;
    Double_t placemaximum = minnn;
    for(Int_t o = 0; o < 10000; o++){
      if(o == 0) maxvalue = polynomeb->Eval(l);
      if(maxvalue < (polynomeb->Eval(l))){
	maxvalue = polynomeb->Eval(l);
	placemaximum = l;
      }
      l += step;
    }
    fPhd[0] = placemaximum;
  }
  
  // Amplification region
  binmax = 0;
  ju     = 0;
  for (Int_t kbin = 1; kbin < projPH->GetNbinsX(); kbin ++) {
    if (((projPH->GetBinContent(kbin+1) - projPH->GetBinContent(kbin)) <= 0.0) && (ju == 0) && (kbin > (fPhd[0]/widbins))) {
      binmax = kbin;
      ju     = 1;
    }
  }
   
  Double_t minn = 0.0;
  Double_t maxx = 0.0;
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = 0.0;
  x[4] = 0.0;
  y[0] = 0.0;
  y[1] = 0.0;
  y[2] = 0.0;
  y[3] = 0.0;
  y[4] = 0.0;

  Int_t    kase1 = nbins - binmax;

  //Determination of minn and maxx
  //case binmax = nbins
  //pol2
  switch(kase1)
    {
    case 0:
      minn = projPH->GetBinCenter(binmax-2);
      maxx = projPH->GetBinCenter(binmax);
      x[0] = projPH->GetBinCenter(binmax-2);
      x[1] = projPH->GetBinCenter(binmax-1);
      x[2] = projPH->GetBinCenter(binmax);
      y[0] = projPH->GetBinContent(binmax-2);
      y[1] = projPH->GetBinContent(binmax-1);
      y[2] = projPH->GetBinContent(binmax);
      CalculPolynomeLagrange2(x,y,c0,c1,c2,c3,c4);
      //AliInfo("At the limit for the drift!");
      break;
    case 1:
      minn = projPH->GetBinCenter(binmax-2);
      maxx = projPH->GetBinCenter(binmax+1);
      x[0] = projPH->GetBinCenter(binmax-2);
      x[1] = projPH->GetBinCenter(binmax-1);
      x[2] = projPH->GetBinCenter(binmax);
      x[3] = projPH->GetBinCenter(binmax+1);
      y[0] = projPH->GetBinContent(binmax-2);
      y[1] = projPH->GetBinContent(binmax-1);
      y[2] = projPH->GetBinContent(binmax);
      y[3] = projPH->GetBinContent(binmax+1);
      CalculPolynomeLagrange3(x,y,c0,c1,c2,c3,c4);
      break;
    default:
      switch(binmax)
	{
	case 0:
	  put = kFALSE;
	  break;
	case 1:
	  minn = projPH->GetBinCenter(binmax);
	  maxx = projPH->GetBinCenter(binmax+2);
	  x[0] = projPH->GetBinCenter(binmax);
	  x[1] = projPH->GetBinCenter(binmax+1);
	  x[2] = projPH->GetBinCenter(binmax+2);
	  y[0] = projPH->GetBinContent(binmax);
	  y[1] = projPH->GetBinContent(binmax+1);
	  y[2] = projPH->GetBinContent(binmax+2);
	  CalculPolynomeLagrange2(x,y,c0,c1,c2,c3,c4);
	  break;
	case 2:
	  minn = projPH->GetBinCenter(binmax-1);
	  maxx = projPH->GetBinCenter(binmax+2);
	  x[0] = projPH->GetBinCenter(binmax-1);
	  x[1] = projPH->GetBinCenter(binmax);
	  x[2] = projPH->GetBinCenter(binmax+1);
	  x[3] = projPH->GetBinCenter(binmax+2);
	  y[0] = projPH->GetBinContent(binmax-1);
	  y[1] = projPH->GetBinContent(binmax);
	  y[2] = projPH->GetBinContent(binmax+1);
	  y[3] = projPH->GetBinContent(binmax+2);
	  CalculPolynomeLagrange3(x,y,c0,c1,c2,c3,c4);
	  break;
	default:
	  minn = projPH->GetBinCenter(binmax-2);
	  maxx = projPH->GetBinCenter(binmax+2);
	  x[0] = projPH->GetBinCenter(binmax-2);
	  x[1] = projPH->GetBinCenter(binmax-1);
	  x[2] = projPH->GetBinCenter(binmax);
	  x[3] = projPH->GetBinCenter(binmax+1);
	  x[4] = projPH->GetBinCenter(binmax+2);
	  y[0] = projPH->GetBinContent(binmax-2);
	  y[1] = projPH->GetBinContent(binmax-1);
	  y[2] = projPH->GetBinContent(binmax);
	  y[3] = projPH->GetBinContent(binmax+1);
	  y[4] = projPH->GetBinContent(binmax+2);
	  CalculPolynomeLagrange4(x,y,c0,c1,c2,c3,c4);
	  break;
	}
      break;
    }
  
  if(put) {
    polynomea = new TF1("pola","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",minn,maxx);
    polynomea->SetParameters(c0,c1,c2,c3,c4);
       
    Double_t step = (maxx-minn)/1000;
    Double_t l = minn;
    Double_t maxvalue = 0.0;
    Double_t placemaximum = minn;
    for(Int_t o = 0; o < 1000; o++){
      if(o == 0) maxvalue = polynomea->Eval(l);
      if(maxvalue < (polynomea->Eval(l))){
	maxvalue = polynomea->Eval(l);
	placemaximum = l;
      }
      l += step;
    }
    fPhd[1] = placemaximum;
  }
  
  // Drift region
  TH1D *pente = new TH1D("pente","pente", projPH->GetNbinsX(),0,(Float_t) limit);
  for (Int_t k = TMath::Min(binmax+4, projPH->GetNbinsX()); k <  projPH->GetNbinsX(); k++) {
    pente->SetBinContent(k,(Double_t) (projPH->GetBinContent(k+1) - projPH->GetBinContent(k)));
  }
  binmin = 0;
  if(pente->GetEntries() > 0) binmin = (Int_t) pente->GetMinimumBin();

  //should not happen
  if (binmin <= 1) {
    binmin = 2;
    put = 1;
    //AliInfo("Put the binmax from 1 to 2 to enable the fit");
  }
  
  //check
  if((projPH->GetBinContent(binmin)-projPH->GetBinError(binmin)) < (projPH->GetBinContent(binmin+1))) {
    //AliInfo("Too many fluctuations at the end!");
    put = kFALSE;
  }
  if((projPH->GetBinContent(binmin)+projPH->GetBinError(binmin)) > (projPH->GetBinContent(binmin-1))) {
    //AliInfo("Too many fluctuations at the end!");
    put = kFALSE;
  }
  if(TMath::Abs(pente->GetBinContent(binmin+1)) <= 0.0000000000001){
    //AliInfo("No entries for the next bin!");
    pente->SetBinContent(binmin,0);
    if(pente->GetEntries() > 0) binmin = (Int_t) pente->GetMinimumBin();
  }

  
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = 0.0;
  x[4] = 0.0;
  y[0] = 0.0;
  y[1] = 0.0;
  y[2] = 0.0;
  y[3] = 0.0;
  y[4] = 0.0;
  Double_t min = 0.0;
  Double_t max = 0.0;
  Bool_t case1 = kFALSE;
  Bool_t case2 = kFALSE;
  Bool_t case4 = kFALSE;

  //Determination of min and max
  //case binmin <= nbins-3
  //pol4 case 3
  if((binmin <= (nbins-3)) && ((binmin-2) >= TMath::Min(binmax+4, projPH->GetNbinsX()))){
    min = pente->GetBinCenter(binmin-2);
    max = pente->GetBinCenter(binmin+2);
    x[0] = pente->GetBinCenter(binmin-2);
    x[1] = pente->GetBinCenter(binmin-1);
    x[2] = pente->GetBinCenter(binmin);
    x[3] = pente->GetBinCenter(binmin+1);
    x[4] = pente->GetBinCenter(binmin+2);
    y[0] = pente->GetBinContent(binmin-2);
    y[1] = pente->GetBinContent(binmin-1);
    y[2] = pente->GetBinContent(binmin);
    y[3] = pente->GetBinContent(binmin+1);
    y[4] = pente->GetBinContent(binmin+2);
    //Calcul the polynome de Lagrange
    CalculPolynomeLagrange4(x,y,c0,c1,c2,c3,c4);
    //richtung +/-
    if((pente->GetBinContent(binmin+2) <= pente->GetBinContent(binmin+1)) &&
       (pente->GetBinContent(binmin-2) <= pente->GetBinContent(binmin-1))) {
      //AliInfo("polynome 4 false 1");
      put = kFALSE;
    }
    if(((binmin+3) <= (nbins-1)) &&
       (pente->GetBinContent(binmin+3) <= pente->GetBinContent(binmin+2)) &&
       ((binmin-3) >= TMath::Min(binmax+4, projPH->GetNbinsX())) &&
       (pente->GetBinContent(binmin-3) <= pente->GetBinContent(binmin-2))) {
      //AliInfo("polynome 4 false 2");
      put = kFALSE;
    }
    // poly 3
    if((pente->GetBinContent(binmin+2) <= pente->GetBinContent(binmin+1)) &&
       (pente->GetBinContent(binmin-2) > pente->GetBinContent(binmin-1))) {
      //AliInfo("polynome 4 case 1");
      case1 = kTRUE;
    }
    if((pente->GetBinContent(binmin+2) > pente->GetBinContent(binmin+1)) &&
       (pente->GetBinContent(binmin-2) <= pente->GetBinContent(binmin-1))) {
      //AliInfo("polynome 4 case 4");
      case4 = kTRUE;
    }
    
  }
  //case binmin = nbins-2
  //pol3 case 1
  if(((binmin == (nbins-2)) && ((binmin-2) >= TMath::Min(binmax+4, projPH->GetNbinsX()))) ||
     (case1)){
    min = pente->GetBinCenter(binmin-2);
    max = pente->GetBinCenter(binmin+1);
    x[0] = pente->GetBinCenter(binmin-2);
    x[1] = pente->GetBinCenter(binmin-1);
    x[2] = pente->GetBinCenter(binmin);
    x[3] = pente->GetBinCenter(binmin+1);
    y[0] = pente->GetBinContent(binmin-2);
    y[1] = pente->GetBinContent(binmin-1);
    y[2] = pente->GetBinContent(binmin);
    y[3] = pente->GetBinContent(binmin+1);
    //Calcul the polynome de Lagrange
    CalculPolynomeLagrange3(x,y,c0,c1,c2,c3,c4);
    //richtung +: nothing
    //richtung -
    if((pente->GetBinContent(binmin-2) <= pente->GetBinContent(binmin-1))) {
      //AliInfo("polynome 3- case 2");
      case2 = kTRUE;
    }
  }
  //pol3 case 4
  if(((binmin <= (nbins-3)) && ((binmin-1) == TMath::Min(binmax+4, projPH->GetNbinsX()))) ||
     (case4)){
    min = pente->GetBinCenter(binmin-1);
    max = pente->GetBinCenter(binmin+2);
    x[0] = pente->GetBinCenter(binmin-1);
    x[1] = pente->GetBinCenter(binmin);
    x[2] = pente->GetBinCenter(binmin+1);
    x[3] = pente->GetBinCenter(binmin+2);
    y[0] = pente->GetBinContent(binmin-1);
    y[1] = pente->GetBinContent(binmin);
    y[2] = pente->GetBinContent(binmin+1);
    y[3] = pente->GetBinContent(binmin+2);
    //Calcul the polynome de Lagrange
    CalculPolynomeLagrange3(x,y,c0,c1,c2,c3,c4);
    //richtung +
    if((pente->GetBinContent(binmin+2) <= pente->GetBinContent(binmin+1))) {
      //AliInfo("polynome 3+ case 2");      
      case2 = kTRUE;
    }
  }
  //pol2 case 5
  if((binmin <= (nbins-3)) && (binmin == TMath::Min(binmax+4, projPH->GetNbinsX()))){
    min = pente->GetBinCenter(binmin);
    max = pente->GetBinCenter(binmin+2);
    x[0] = pente->GetBinCenter(binmin);
    x[1] = pente->GetBinCenter(binmin+1);
    x[2] = pente->GetBinCenter(binmin+2);
    y[0] = pente->GetBinContent(binmin);
    y[1] = pente->GetBinContent(binmin+1);
    y[2] = pente->GetBinContent(binmin+2);
    //Calcul the polynome de Lagrange
    CalculPolynomeLagrange2(x,y,c0,c1,c2,c3,c4);
    //richtung +
    if((pente->GetBinContent(binmin+2) <= pente->GetBinContent(binmin+1))) {
      //AliInfo("polynome 2+ false");
      put = kFALSE;
    }
  }
  //pol2 case 2
  if(((binmin == (nbins-2)) && ((binmin-1) == TMath::Min(binmax+4, projPH->GetNbinsX()))) ||
     (case2)){
    min = pente->GetBinCenter(binmin-1);
    max = pente->GetBinCenter(binmin+1);
    x[0] = pente->GetBinCenter(binmin-1);
    x[1] = pente->GetBinCenter(binmin);
    x[2] = pente->GetBinCenter(binmin+1);
    y[0] = pente->GetBinContent(binmin-1);
    y[1] = pente->GetBinContent(binmin);
    y[2] = pente->GetBinContent(binmin+1);
    //Calcul the polynome de Lagrange
    CalculPolynomeLagrange2(x,y,c0,c1,c2,c3,c4);
    //richtung +: nothing
    //richtung -: nothing
  }
  //case binmin = nbins-1
  //pol2 case 0
  if((binmin == (nbins-1)) && ((binmin-2) >= TMath::Min(binmax+4, projPH->GetNbinsX()))){
    min = pente->GetBinCenter(binmin-2);
    max = pente->GetBinCenter(binmin);
    x[0] = pente->GetBinCenter(binmin-2);
    x[1] = pente->GetBinCenter(binmin-1);
    x[2] = pente->GetBinCenter(binmin);
    y[0] = pente->GetBinContent(binmin-2);
    y[1] = pente->GetBinContent(binmin-1);
    y[2] = pente->GetBinContent(binmin);
    //Calcul the polynome de Lagrange
    CalculPolynomeLagrange2(x,y,c0,c1,c2,c3,c4);
    //AliInfo("At the limit for the drift!");
    //fluctuation too big!
    //richtung +: nothing
    //richtung -
    if((pente->GetBinContent(binmin-2) <= pente->GetBinContent(binmin-1))) {
      //AliInfo("polynome 2- false ");
      put = kFALSE;
    }
  }
  if((binmin == (nbins-1)) && ((binmin-2) < TMath::Min(binmax+4, projPH->GetNbinsX()))) {
    put = kFALSE;
    //AliInfo("At the limit for the drift and not usable!");
  }

  //pass
  if((binmin == (nbins-2)) && ((binmin-1) < TMath::Min(binmax+4, projPH->GetNbinsX()))){
    put = kFALSE;
    //AliInfo("For the drift...problem!");
  }
  //pass but should not happen
  if((binmin <= (nbins-3)) && (binmin < TMath::Min(binmax+6, projPH->GetNbinsX()))){
    put = kFALSE;
    //AliInfo("For the drift...problem!");
  }
  
  if(put) {
    polynome = new TF1("pol","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",min,max);
    polynome->SetParameters(c0,c1,c2,c3,c4);
    //AliInfo(Form("GetMinimum of the function %f",polynome->GetMinimumX()));
    Double_t step = (max-min)/1000;
    Double_t l = min;
    Double_t minvalue = 0.0;
    Double_t placeminimum = min;
    for(Int_t o = 0; o < 1000; o++){
      if(o == 0) minvalue = polynome->Eval(l);
      if(minvalue > (polynome->Eval(l))){
	minvalue = polynome->Eval(l);
	placeminimum = l;
      }
      l += step;
    }
    fPhd[2] = placeminimum;
  }
  //printf("La fin %d\n",((Int_t)(fPhd[2]*10.0))+2);
  if((((Int_t)(fPhd[2]*10.0))+2) >= projPH->GetNbinsX()) fPhd[2] = 0.0;
  if(((((Int_t)(fPhd[2]*10.0))+2) < projPH->GetNbinsX()) && (projPH->GetBinContent(((Int_t)(fPhd[2]*10.0))+2)==0)) fPhd[2] = 0.0;
  
  Float_t fPhdt0  = 0.0;
  Float_t t0Shift = 0.0;
  if(fTakeTheMaxPH) {
    fPhdt0 = fPhd[1];
    t0Shift = fT0Shift1;
  }
  else {
    fPhdt0 = fPhd[0];
    t0Shift = fT0Shift0;
  }

  if ((fPhd[2] > fPhd[0]) && 
      (fPhd[2] > fPhd[1]) && 
      (fPhd[1] > fPhd[0]) &&
      (put)) {
    fCurrentCoef[0] = (kDrWidth) / (fPhd[2]-fPhd[1]);
    if(fCurrentCoef[0] > 2.5) fCurrentCoef[0] =  -TMath::Abs(fCurrentCoef[1]);
    else fNumberFitSuccess++;
    if (fPhdt0 >= 0.0) {
      fCurrentCoef2[0] = (fPhdt0 - t0Shift) / widbins;
      //printf("Value of timeoffset %f\n",fCurrentCoef2[0]);
      if (fCurrentCoef2[0] < -3.0) {
        fCurrentCoef2[0] = fCurrentCoef2[1] + 100.0;
      }
    }
    else {
      fCurrentCoef2[0] = fCurrentCoef2[1] + 100.0;
    }
  }
  else {
    ////printf("Put default %f\n",-TMath::Abs(fCurrentCoef[1]));
    fCurrentCoef[0]      = -TMath::Abs(fCurrentCoef[1]);
    
    if((fPhd[1] > fPhd[0]) &&
       (put)) {
      if (fPhdt0 >= 0.0) {
	fCurrentCoef2[0] = (fPhdt0 - t0Shift) / widbins;
	if (fCurrentCoef2[0] < -3.0) {
	  fCurrentCoef2[0] = fCurrentCoef2[1] + 100.0;
	}
      }
      else {
	fCurrentCoef2[0] = fCurrentCoef2[1] + 100.0;
      }
    }
    else{
      fCurrentCoef2[0]     = fCurrentCoef2[1] + 100.0;
      //printf("Fit failed!\n");
    }
  }
  
  if (fDebugLevel == 1) {
    TCanvas *cpentei = new TCanvas("cpentei","cpentei",50,50,600,800);
    cpentei->cd();
    projPH->Draw();
    line->SetLineColor(2);
    line->DrawLine(fPhd[0],0,fPhd[0],projPH->GetMaximum());
    line->DrawLine(fPhd[1],0,fPhd[1],projPH->GetMaximum());
    line->DrawLine(fPhd[2],0,fPhd[2],projPH->GetMaximum());
    AliInfo(Form("fPhd[0] (beginning of the signal): %f"                  ,(Float_t) fPhd[0]));
    AliInfo(Form("fPhd[1] (end of the amplification region): %f"          ,(Float_t) fPhd[1]));
    AliInfo(Form("fPhd[2] (end of the drift region): %f"                  ,(Float_t) fPhd[2]));
    AliInfo(Form("fVriftCoef[3] (with only the drift region(default)): %f",(Float_t) fCurrentCoef[0]));
    TCanvas *cpentei2 = new TCanvas("cpentei2","cpentei2",50,50,600,800);
    cpentei2->cd();
    pentea->Draw();
    TCanvas *cpentei3 = new TCanvas("cpentei3","cpentei3",50,50,600,800);
    cpentei3->cd();
    pente->Draw();
  }
  else {
    delete pentea;
    delete pente;
    if(polynome) delete polynome;
    if(polynomea) delete polynomea;
    if(polynomeb) delete polynomeb;
    //if(x) delete [] x;
    //if(y) delete [] y;
    if(line) delete line;

  }

  //Provisoire
  //if(fCurrentCoef[0] > 1.7) fCurrentCoef[0] = -TMath::Abs(fCurrentCoef[1]);
  //if((fCurrentCoef2[0] > 2.6) || (fCurrentCoef2[0] < 2.1)) fCurrentCoef2[0] = fCurrentCoef2[1] + 100.0;
  //printf("Value of timeoffset final %f\n",fCurrentCoef2[0]);
  projPH->SetDirectory(0);

}

//_____________________________________________________________________________
void AliTRDCalibraFit::FitPH(TH1* projPH, Int_t idect)
{
  //
  // Fit methode for the drift velocity
  //
  
  // Constants
  const Float_t kDrWidth = AliTRDgeometry::DrThick();  

  // Some variables
  TAxis   *xpph   = projPH->GetXaxis();
  Double_t upedge = xpph->GetBinUpEdge(xpph->GetNbins());

  TF1 *fPH = new TF1("fPH",AliTRDCalibraFit::PH,-0.05,3.2,6);
  fPH->SetParameter(0,0.469);     // Scaling
  fPH->SetParameter(1,0.18);      // Start 
  fPH->SetParameter(2,0.0857325); // AR
  fPH->SetParameter(3,1.89);      // DR
  fPH->SetParameter(4,0.08);      // QA/QD
  fPH->SetParameter(5,0.0);       // Baseline

  TLine *line = new TLine();

  fCurrentCoef[0]     = 0.0;
  fCurrentCoef2[0]    = 0.0;
  fCurrentCoefE       = 0.0;
  fCurrentCoefE2      = 0.0;
 
  if (idect%fFitPHPeriode == 0) {

    AliInfo(Form("The detector %d will be fitted",idect));
    fPH->SetParameter(0,(projPH->Integral()-(projPH->GetBinContent(1)*projPH->GetNbinsX())) * 0.00028); // Scaling
    fPH->SetParameter(1,fPhd[0] - 0.1);                                                                 // Start 
    fPH->SetParameter(2,fPhd[1] - fPhd[0]);                                                             // AR
    fPH->SetParameter(3,fPhd[2] - fPhd[1]);                                                             // DR
    fPH->SetParameter(4,0.225);                                                                         // QA/QD
    fPH->SetParameter(5,(Float_t) projPH->GetBinContent(1));
    
    if (fDebugLevel != 1) {
      projPH->Fit(fPH,"0M","",0.0,upedge);
    }
    else {
      TCanvas *cpente = new TCanvas("cpente","cpente",50,50,600,800);
      cpente->cd();
      projPH->Fit(fPH,"M+","",0.0,upedge);
      projPH->Draw("E0");
      line->SetLineColor(4);
      line->DrawLine(fPH->GetParameter(1)
                    ,0
                    ,fPH->GetParameter(1)
                    ,projPH->GetMaximum());
      line->DrawLine(fPH->GetParameter(1)+fPH->GetParameter(2)
                    ,0
                    ,fPH->GetParameter(1)+fPH->GetParameter(2)
                    ,projPH->GetMaximum());
      line->DrawLine(fPH->GetParameter(1)+fPH->GetParameter(2)+fPH->GetParameter(3)
                    ,0
                    ,fPH->GetParameter(1)+fPH->GetParameter(2)+fPH->GetParameter(3)
                    ,projPH->GetMaximum());
    }

    if (fPH->GetParameter(3) != 0) {
      fNumberFitSuccess++;
      fCurrentCoef[0]    = kDrWidth / (fPH->GetParameter(3));
      fCurrentCoefE      = (fPH->GetParError(3)/fPH->GetParameter(3))*fCurrentCoef[0];
      fCurrentCoef2[0]   = fPH->GetParameter(1);
      fCurrentCoefE2     = fPH->GetParError(1);
    } 
    else {
      fCurrentCoef[0]     = -TMath::Abs(fCurrentCoef[1]);
      fCurrentCoef2[0]    = fCurrentCoef2[1] + 100.0;
    }
 
  }
  else {

    // Put the default value
    fCurrentCoef[0]  = -TMath::Abs(fCurrentCoef[1]);
    fCurrentCoef2[0] = fCurrentCoef2[1] + 100.0;
  }

  if (fDebugLevel != 1) {
    delete fPH;
  }
  
}
//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::FitPRFGausMI(Double_t *arraye, Double_t *arraym, Double_t *arrayme, Int_t nBins, Float_t xMin, Float_t xMax)
{
  //
  // Fit methode for the sigma of the pad response function
  //

  TVectorD param(3);
  
  fCurrentCoef[0]  = 0.0;
  fCurrentCoefE = 0.0;

  Double_t ret = FitGausMI(arraye, arraym, arrayme, nBins, xMin, xMax,&param); 

  if(TMath::Abs(ret+4) <= 0.000000001){
    fCurrentCoef[0] = -fCurrentCoef[1];
    return kFALSE;
  }
  else {
    fNumberFitSuccess++;
    fCurrentCoef[0] = param[2];
    fCurrentCoefE   = ret;
    return kTRUE;
  }
}
//_____________________________________________________________________________
Double_t AliTRDCalibraFit::FitGausMI(Double_t *arraye, Double_t *arraym, Double_t *arrayme, Int_t nBins, Float_t xMin, Float_t xMax, TVectorD *param, Bool_t bError)
{
  //
  // Fit methode for the sigma of the pad response function
  //

  //We should have at least 3 points
  if(nBins <=3) return -4.0;

  TLinearFitter fitter(3,"pol2");
  fitter.StoreData(kFALSE);
  fitter.ClearPoints();
  TVectorD  par(3);
  Float_t binWidth = (xMax-xMin)/(Float_t)nBins;
  Float_t entries = 0;
  Int_t   nbbinwithentries = 0;
  for (Int_t i=0; i<nBins; i++){
    entries+=arraye[i];
    if(arraye[i] > 15) nbbinwithentries++;
    //printf("entries for i %d: %f\n",i,arraye[i]);
  }
  if ((entries<700) || (nbbinwithentries < ((Int_t)(nBins/2)))) return -4;
  //printf("entries %f\n",entries);
  //printf("nbbinwithentries %d\n",nbbinwithentries);  

  Int_t npoints=0;
  Float_t errorm = 0.0;
  Float_t errorn = 0.0;
  Float_t error  = 0.0;
  
  //
  for (Int_t ibin=0;ibin<nBins; ibin++){
      Float_t entriesI = arraye[ibin];
      Float_t valueI   = arraym[ibin];
      Double_t xcenter = 0.0;
      Float_t  val     = 0.0;
      if ((entriesI>15) && (valueI>0.0)){
      	xcenter = xMin+(ibin+0.5)*binWidth;
	errorm   = 0.0;
	errorn   = 0.0;
	error    = 0.0;
	if(!bError){
	  if((valueI + 0.01) > 0.0) errorm = TMath::Log((valueI + 0.01)/valueI);
	  if((valueI - 0.01) > 0.0) errorn = TMath::Log((valueI - 0.01)/valueI);
	  error = TMath::Max(TMath::Abs(errorm),TMath::Abs(errorn));
	}
	else{
	  if((valueI + arrayme[ibin]) > 0.0) errorm = TMath::Log((valueI + arrayme[ibin])/valueI);
	  if((valueI - arrayme[ibin]) > 0.0) errorn = TMath::Log((valueI - arrayme[ibin])/valueI);
	  error = TMath::Max(TMath::Abs(errorm),TMath::Abs(errorn));
	}
	if(TMath::Abs(error) < 0.000000001) continue;
	val      = TMath::Log(Float_t(valueI));
	fitter.AddPoint(&xcenter,val,error);
	npoints++;
      }

      if(fDebugLevel > 1){

      if ( !fDebugStreamer ) {
	//debug stream
	TDirectory *backup = gDirectory;
	fDebugStreamer = new TTreeSRedirector("TRDDebugFitPRF.root");
	if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
      } 
      
      Int_t    detector     = fCountDet;
      Int_t    layer        = GetLayer(fCountDet);
      Int_t    group        = ibin;    
     
      (* fDebugStreamer) << "FitGausMIFill"<<
	"detector="<<detector<<
	"layer="<<layer<<
	"nbins="<<nBins<<
	"group="<<group<<
	"entriesI="<<entriesI<<
	"valueI="<<valueI<<
	"val="<<val<<
	"xcenter="<<xcenter<<
	"errorm="<<errorm<<
	"errorn="<<errorn<<
	"error="<<error<<
	"bError="<<bError<<
	"\n";  
    }

  }

  if(npoints <=3) return -4.0;  

  Double_t chi2 = 0;
  if (npoints>3){
    fitter.Eval();
    fitter.GetParameters(par);
    chi2 = fitter.GetChisquare()/Float_t(npoints);
    
        
    if (!param)  param  = new TVectorD(3);
    if(TMath::Abs(par[2]) <= 0.000000001) return -4.0;
    Double_t  x      = TMath::Sqrt(TMath::Abs(-2*par[2])); 
    Double_t deltax = (fitter.GetParError(2))/x;
    Double_t errorparam2 = TMath::Abs(deltax)/(x*x);
    chi2 = errorparam2;
    
    (*param)[1] = par[1]/(-2.*par[2]);
    (*param)[2] = 1./TMath::Sqrt(TMath::Abs(-2.*par[2]));
    Double_t lnparam0 = par[0]+ par[1]* (*param)[1] +  par[2]*(*param)[1]*(*param)[1];
    if ( lnparam0>307 ) return -4;
    (*param)[0] = TMath::Exp(lnparam0);

    if(fDebugLevel > 1){

      if ( !fDebugStreamer ) {
	//debug stream
	TDirectory *backup = gDirectory;
	fDebugStreamer = new TTreeSRedirector("TRDDebugFitPRF.root");
	if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
      } 
      
      Int_t    detector     = fCountDet;
      Int_t    layer        = GetLayer(fCountDet);
           
     
      (* fDebugStreamer) << "FitGausMIFit"<<
	"detector="<<detector<<
	"layer="<<layer<<
	"nbins="<<nBins<<
	"errorsigma="<<chi2<<
	"mean="<<(*param)[1]<<
	"sigma="<<(*param)[2]<<
	"constant="<<(*param)[0]<<
	"\n";  
    }
  }

  if((chi2/(*param)[2]) > 0.1){
    if(bError){
      chi2 = FitGausMI(arraye,arraym,arrayme,nBins,xMin,xMax,param,kFALSE);
    }
    else return -4.0;
  }

  if(fDebugLevel == 1){
    TString name("PRF");
    name += (Int_t)xMin;
    name += (Int_t)xMax;  
    TCanvas *c1 = new TCanvas((const char *)name,(const char *)name,50,50,600,800);  
    c1->cd();
    name += "histo";
    TH1F *histo = new TH1F((const char *)name,(const char *)name,nBins,xMin,xMax);
    for(Int_t k = 0; k < nBins; k++){
      histo->SetBinContent(k+1,arraym[k]);
      histo->SetBinError(k+1,arrayme[k]);
    }
    histo->Draw();
    name += "functionf";
    TF1 *f1= new TF1((const char*)name,"[0]*exp(-(x-[1])^2/(2*[2]*[2]))",xMin,xMax);
    f1->SetParameter(0, (*param)[0]);
    f1->SetParameter(1, (*param)[1]);
    f1->SetParameter(2, (*param)[2]);
    f1->Draw("same");
  }

  
  return chi2;
 
}
//_____________________________________________________________________________
void AliTRDCalibraFit::FitPRF(TH1 *projPRF)
{
  //
  // Fit methode for the sigma of the pad response function
  //
  
  fCurrentCoef[0]  = 0.0;
  fCurrentCoefE = 0.0;

  if (fDebugLevel != 1) {
    projPRF->Fit("gaus","0M","",-fRangeFitPRF,fRangeFitPRF);
  }
  else {
    TCanvas *cfit = new TCanvas("cfit","cfit",50,50,600,800);
    cfit->cd();
    projPRF->Fit("gaus","M+","",-fRangeFitPRF,fRangeFitPRF);
    projPRF->Draw();
  }
  fCurrentCoef[0]  = projPRF->GetFunction("gaus")->GetParameter(2);
  fCurrentCoefE = projPRF->GetFunction("gaus")->GetParError(2);
  if(fCurrentCoef[0] <= 0.0) fCurrentCoef[0] = -fCurrentCoef[1];
  else {
    fNumberFitSuccess++;
  }
}
//_____________________________________________________________________________
void AliTRDCalibraFit::RmsPRF(TH1 *projPRF)
{
  //
  // Fit methode for the sigma of the pad response function
  //
  fCurrentCoef[0]   = 0.0;
  fCurrentCoefE  = 0.0;
  if (fDebugLevel == 1) {
    TCanvas *cfit = new TCanvas("cfit","cfit",50,50,600,800);
    cfit->cd();
    projPRF->Draw();
  }
  fCurrentCoef[0] = projPRF->GetRMS();
  if(fCurrentCoef[0] <= 0.0) fCurrentCoef[0] = -fCurrentCoef[1];
  else {
    fNumberFitSuccess++;
  }
}
//_____________________________________________________________________________
void AliTRDCalibraFit::FitTnpRange(Double_t *arraye, Double_t *arraym, Double_t *arrayme, Int_t nbg, Int_t nybins)
{
  //
  // Fit methode for the sigma of the pad response function with 2*nbg tan bins
  //
  
  TLinearFitter linearfitter = TLinearFitter(3,"pol2");
 

  Int_t   nbins    = (Int_t)(nybins/(2*nbg));
  Float_t lowedge  = -3.0*nbg;
  Float_t upedge   = lowedge + 3.0; 
  Int_t   offset   = 0;
  Int_t   npoints  = 0;
  Double_t xvalues = -0.2*nbg+0.1;
  Double_t y       = 0.0;
  Int_t   total    = 2*nbg;

  
  for(Int_t k = 0; k < total; k++){
    if(FitPRFGausMI(arraye+offset, arraym+offset, arrayme+offset, nbins, lowedge, upedge)){
      npoints++;
      y = fCurrentCoef[0]*fCurrentCoef[0];
      linearfitter.AddPoint(&xvalues,y,2*fCurrentCoefE*fCurrentCoef[0]);
    }
    
    if(fDebugLevel > 1){

      if ( !fDebugStreamer ) {
	//debug stream
	TDirectory *backup = gDirectory;
	fDebugStreamer = new TTreeSRedirector("TRDDebugFitPRF.root");
	if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
      } 
      
      Int_t    detector     = fCountDet;
      Int_t    layer        = GetLayer(fCountDet);
      Int_t    nbtotal      = total;  
      Int_t    group        = k;    
      Float_t  low          = lowedge;
      Float_t  up           = upedge;
      Float_t  tnp          = xvalues;
      Float_t  wid          = fCurrentCoef[0];
      Float_t  widfE        = fCurrentCoefE;

      (* fDebugStreamer) << "FitTnpRange0"<<
	"detector="<<detector<<
	"layer="<<layer<<
	"nbtotal="<<nbtotal<<
	"group="<<group<<
	"low="<<low<<
	"up="<<up<<
	"offset="<<offset<<
	"tnp="<<tnp<<
	"wid="<<wid<<
	"widfE="<<widfE<<
	"\n";  
    }
    
    offset  += nbins;
    lowedge += 3.0;
    upedge  += 3.0;
    xvalues += 0.2;

  }

  fCurrentCoefE = 0.0;
  fCurrentCoef[0] = 0.0;

  //printf("npoints\n",npoints);

  if(npoints < 3){
    fCurrentCoef[0] = -TMath::Abs(fCurrentCoef[1]);
  }
  else{
  
    TVectorD pars0;
    linearfitter.Eval();
    linearfitter.GetParameters(pars0);
    Double_t pointError0  =  TMath::Sqrt(linearfitter.GetChisquare()/npoints);
    Double_t errorsx0     =  linearfitter.GetParError(2)*pointError0;
    Double_t min0         = 0.0;
    Double_t ermin0       = 0.0;
    //Double_t prfe0      = 0.0;
    Double_t prf0         = 0.0;
    if((pars0[2] > 0.000000000001) && (TMath::Abs(pars0[1]) >= 0.000000000001)) {
      min0 = -pars0[1]/(2*pars0[2]);
      ermin0 = TMath::Abs(min0*(errorsx0/pars0[2]+linearfitter.GetParError(1)*pointError0/pars0[1]));
      prf0 = pars0[0]+pars0[1]*min0+pars0[2]*min0*min0;
      if(prf0 > 0.0) {
	/*
	  prfe0 = linearfitter->GetParError(0)*pointError0
	  +(linearfitter->GetParError(1)*pointError0/pars0[1]+ermin0/min0)*pars0[1]*min0
	  +(linearfitter->GetParError(2)*pointError0/pars0[2]+2*ermin0/min0)*pars0[2]*min0*min0;
	  prfe0 = prfe0/(2*TMath::Sqrt(prf0));
	  fCurrentCoefE   = (Float_t) prfe0;
	*/
	fCurrentCoef[0] = (Float_t) TMath::Sqrt(TMath::Abs(prf0));
      }
      else{
	fCurrentCoef[0] = -TMath::Abs(fCurrentCoef[1]);
      }
    }
    else {
      fCurrentCoef[0] = -TMath::Abs(fCurrentCoef[1]);
    }

    if(fDebugLevel > 1){

      if ( !fDebugStreamer ) {
	//debug stream
	TDirectory *backup = gDirectory;
	fDebugStreamer = new TTreeSRedirector("TRDDebugFitPRF.root");
	if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
      } 
      
      Int_t    detector     = fCountDet;
      Int_t    layer        = GetLayer(fCountDet);
      Int_t    nbtotal      = total;
      Double_t colsize[6]   = {0.635,0.665,0.695,0.725,0.755,0.785};  
      Double_t sigmax       = TMath::Sqrt(TMath::Abs(pars0[2]))*10000*colsize[layer];      

      (* fDebugStreamer) << "FitTnpRange1"<<
	"detector="<<detector<<
	"layer="<<layer<<
	"nbtotal="<<nbtotal<<
	"par0="<<pars0[0]<<
	"par1="<<pars0[1]<<
	"par2="<<pars0[2]<<
	"npoints="<<npoints<<
	"sigmax="<<sigmax<<
	"tan="<<min0<<
	"sigmaprf="<<fCurrentCoef[0]<<
	"sigprf="<<fCurrentCoef[1]<<
	"\n";  
    }
    
  }
  
}
//_____________________________________________________________________________
void AliTRDCalibraFit::FitMean(TH1 *projch, Double_t nentries, Double_t mean)
{
  //
  // Only mean methode for the gain factor
  //
 
  fCurrentCoef[0] = mean;
  fCurrentCoefE   = 0.0;
  if(nentries > 0) fCurrentCoefE = projch->GetRMS()/TMath::Sqrt(nentries);
  if (fDebugLevel == 1) {
    TCanvas *cpmean = new TCanvas("cpmean","cpmean",50,50,600,800);
    cpmean->cd();
    projch->Draw();
  }
  CalculChargeCoefMean(kTRUE);
  fNumberFitSuccess++;
}
//_____________________________________________________________________________
void AliTRDCalibraFit::FitMeanW(TH1 *projch, Double_t nentries)
{
  //
  // mean w methode for the gain factor
  //

  //Number of bins
  Int_t nybins = projch->GetNbinsX();
 
  //The weight function
  Double_t a = 0.00228515;
  Double_t b = -0.00231487;
  Double_t c = 0.00044298;
  Double_t d = -0.00379239;
  Double_t e = 0.00338349;

//   0 |0.00228515
//    1 |-0.00231487
//    2 |0.00044298
//    3 |-0.00379239
//    4 |0.00338349



  //A arbitrary error for the moment
  fCurrentCoefE = 0.0;
  fCurrentCoef[0] = 0.0;
  
  //Calcul 
  Double_t sumw = 0.0;
  Double_t sum = 0.0; 
  Float_t sumAll   = (Float_t) nentries;
  Int_t sumCurrent = 0;
  for(Int_t k = 0; k <nybins; k++){
    Double_t fraction = Float_t(sumCurrent)/Float_t(sumAll);
    if (fraction>0.95) break;
    Double_t weight = a + b*fraction + c*fraction*fraction + d *fraction*fraction*fraction+
      e*fraction*fraction*fraction*fraction;
    sumw += weight*projch->GetBinContent(k+1)*projch->GetBinCenter(k+1);
    sum  += weight*projch->GetBinContent(k+1); 
    sumCurrent += (Int_t) projch->GetBinContent(k+1);
    //printf("fraction %f, weight %f, bincontent %f\n",fraction,weight,projch->GetBinContent(k+1));   
  }
  if(sum > 0.0) fCurrentCoef[0] = (sumw/sum);

  if (fDebugLevel == 1) {
    TCanvas *cpmeanw = new TCanvas("cpmeanw","cpmeanw",50,50,600,800);
    cpmeanw->cd();
    projch->Draw();
  }
  fNumberFitSuccess++;
  CalculChargeCoefMean(kTRUE);
}
//_____________________________________________________________________________
void AliTRDCalibraFit::FitMeanWSm(TH1 *projch, Float_t sumAll)
{
  //
  // mean w methode for the gain factor
  //

  //Number of bins
  Int_t nybins = projch->GetNbinsX();
 
  //The weight function
  Double_t a = 0.00228515;
  Double_t b = -0.00231487;
  Double_t c = 0.00044298;
  Double_t d = -0.00379239;
  Double_t e = 0.00338349;

//   0 |0.00228515
//    1 |-0.00231487
//    2 |0.00044298
//    3 |-0.00379239
//    4 |0.00338349



  //A arbitrary error for the moment
  fCurrentCoefE = 0.0;
  fCurrentCoef[0] = 0.0;
  
  //Calcul 
  Double_t sumw = 0.0;
  Double_t sum = 0.0; 
  Int_t sumCurrent = 0;
  for(Int_t k = 0; k <nybins; k++){
    Double_t fraction = Float_t(sumCurrent)/Float_t(sumAll);
    if (fraction>0.95) break;
    Double_t weight = a + b*fraction + c*fraction*fraction + d *fraction*fraction*fraction+
      e*fraction*fraction*fraction*fraction;
    sumw += weight*projch->GetBinContent(k+1)*projch->GetBinCenter(k+1);
    sum  += weight*projch->GetBinContent(k+1); 
    sumCurrent += (Int_t) projch->GetBinContent(k+1);
    //printf("fraction %f, weight %f, bincontent %f\n",fraction,weight,projch->GetBinContent(k+1));   
  }
  if(sum > 0.0) fCurrentCoef[0] = (sumw/sum);

  if (fDebugLevel == 1) {
    TCanvas *cpmeanw = new TCanvas("cpmeanw","cpmeanw",50,50,600,800);
    cpmeanw->cd();
    projch->Draw();
  }
  fNumberFitSuccess++;
}
//_____________________________________________________________________________
void AliTRDCalibraFit::FitCH(TH1 *projch, Double_t mean)
{
  //
  // Fit methode for the gain factor
  //
 
  fCurrentCoef[0]  = 0.0;
  fCurrentCoefE    = 0.0;
  Double_t chisqrl = 0.0;
  Double_t chisqrg = 0.0;
  Double_t chisqr  = 0.0;
  TF1 *fLandauGaus = new TF1("fLandauGaus",FuncLandauGaus,0,300,5);

  projch->Fit("landau","0",""
             ,(Double_t) mean/fBeginFitCharge
             ,projch->GetBinCenter(projch->GetNbinsX()));
  Double_t l3P0         = projch->GetFunction("landau")->GetParameter(0);
  Double_t l3P1         = projch->GetFunction("landau")->GetParameter(1);
  Double_t l3P2         = projch->GetFunction("landau")->GetParameter(2);
  chisqrl = projch->GetFunction("landau")->GetChisquare();
    
  projch->Fit("gaus","0",""
	      ,(Double_t) mean/fBeginFitCharge
	      ,projch->GetBinCenter(projch->GetNbinsX()));
  Double_t g3P0         = projch->GetFunction("gaus")->GetParameter(0);
  Double_t g3P2         = projch->GetFunction("gaus")->GetParameter(2);
  chisqrg = projch->GetFunction("gaus")->GetChisquare();
        
  fLandauGaus->SetParameters(l3P0,l3P1,l3P2,g3P0,g3P2);
  if (fDebugLevel != 1) {
    projch->Fit("fLandauGaus","0",""
		,(Double_t) mean/fBeginFitCharge
		,projch->GetBinCenter(projch->GetNbinsX()));
    chisqr = projch->GetFunction("fLandauGaus")->GetChisquare();
  } 
  else  {
    TCanvas *cp = new TCanvas("cp","cp",50,50,600,800);
    cp->cd();
    projch->Fit("fLandauGaus","+",""
		,(Double_t) mean/fBeginFitCharge
		,projch->GetBinCenter(projch->GetNbinsX()));
    chisqr = projch->GetFunction("fLandauGaus")->GetChisquare();
    projch->Draw();
    fLandauGaus->Draw("same");
  }
  
  if ((projch->GetFunction("fLandauGaus")->GetParameter(1) > 0) && (projch->GetFunction("fLandauGaus")->GetParError(1) < (0.05*projch->GetFunction("fLandauGaus")->GetParameter(1))) && (chisqr < chisqrl) && (chisqr < chisqrg)) {
    //if ((projch->GetFunction("fLandauGaus")->GetParameter(1) > 0) && (chisqr < chisqrl) && (chisqr < chisqrg)) {
    fNumberFitSuccess++;
    CalculChargeCoefMean(kTRUE);
    fCurrentCoef[0]  = projch->GetFunction("fLandauGaus")->GetParameter(1);
    fCurrentCoefE    = projch->GetFunction("fLandauGaus")->GetParError(1);
  }
  else {
    CalculChargeCoefMean(kFALSE);
    fCurrentCoef[0] = -TMath::Abs(fCurrentCoef[1]);
  }
   
  if (fDebugLevel != 1) {
    delete fLandauGaus;
  }

}
//_____________________________________________________________________________
void AliTRDCalibraFit::FitBisCH(TH1* projch, Double_t mean)
{
  //
  // Fit methode for the gain factor more time consuming
  //


  //Some parameters to initialise
  Double_t widthLandau, widthGaus, mPV, integral;
  Double_t chisquarel = 0.0;
  Double_t chisquareg = 0.0;
  projch->Fit("landau","0M+",""
	      ,(Double_t) mean/6
	      ,projch->GetBinCenter(projch->GetNbinsX()));
  widthLandau  = projch->GetFunction("landau")->GetParameter(2);
  chisquarel = projch->GetFunction("landau")->GetChisquare();
  projch->Fit("gaus","0M+",""
	      ,(Double_t) mean/6
	      ,projch->GetBinCenter(projch->GetNbinsX()));
  widthGaus    = projch->GetFunction("gaus")->GetParameter(2);
  chisquareg = projch->GetFunction("gaus")->GetChisquare();
    
  mPV = (projch->GetFunction("landau")->GetParameter(1))/2;
  integral = (projch->GetFunction("gaus")->Integral(0.3*mean,3*mean)+projch->GetFunction("landau")->Integral(0.3*mean,3*mean))/2;
  
  // Setting fit range and start values
  Double_t fr[2];
  //Double_t sv[4] = { l3P2, fChargeCoef[1], projch->Integral("width"), fG3P2 };
  //Double_t sv[4]   = { fL3P2, fChargeCoef[1], fL3P0, fG3P2 };
  Double_t sv[4]   = { widthLandau, mPV, integral, widthGaus};
  Double_t pllo[4] = { 0.001, 0.001, projch->Integral()/3, 0.001};
  Double_t plhi[4] = { 300.0, 300.0, 30*projch->Integral(), 300.0};
  Double_t fp[4]   = { 1.0, 1.0, 1.0, 1.0 };
  Double_t fpe[4]  = { 1.0, 1.0, 1.0, 1.0 };
  fr[0]            = 0.3 * mean;
  fr[1]            = 3.0 * mean;
  fCurrentCoef[0]  = 0.0;
  fCurrentCoefE    = 0.0;

  Double_t chisqr;
  Int_t    ndf;
  TF1 *fitsnr = LanGauFit(projch,&fr[0],&sv[0]
                                ,&pllo[0],&plhi[0]
                                ,&fp[0],&fpe[0]
                                ,&chisqr,&ndf);
    
  Double_t projchPeak;
  Double_t projchFWHM;
  LanGauPro(fp,projchPeak,projchFWHM);

  if ((fp[1] > 0) && ((fpe[1] < (0.05*fp[1])) && (chisqr < chisquarel) && (chisqr < chisquareg))) {
    //if ((fp[1] > 0) && ((chisqr < chisquarel) && (chisqr < chisquareg))) {
    fNumberFitSuccess++;
    CalculChargeCoefMean(kTRUE);
    fCurrentCoef[0]  = fp[1];
    fCurrentCoefE = fpe[1];
    //chargeCoefE2 = chisqr;
  } 
  else {
    CalculChargeCoefMean(kFALSE);
    fCurrentCoef[0] = -TMath::Abs(fCurrentCoef[1]);
  }
  if (fDebugLevel == 1) {
    AliInfo(Form("fChargeCoef[0]: %f",(Float_t) fCurrentCoef[0]));
    TCanvas *cpy = new TCanvas("cpy","cpy",50,50,600,800);
    cpy->cd();
    projch->Draw();
    fitsnr->Draw("same");
  }
  else {
    delete fitsnr;
  }
} 
//_____________________________________________________________________________
void AliTRDCalibraFit::CalculPolynomeLagrange2(const Double_t *x, const Double_t *y, Double_t &c0, Double_t &c1, Double_t &c2, Double_t &c3, Double_t &c4) const
{
  //
  // Calcul the coefficients of the polynome passant par ces trois points de degre 2
  //
  Double_t x0 = y[0]/((x[0]-x[1])*(x[0]-x[2]));
  Double_t x1 = y[1]/((x[1]-x[0])*(x[1]-x[2]));
  Double_t x2 = y[2]/((x[2]-x[0])*(x[2]-x[1]));

  c4 = 0.0;
  c3 = 0.0;
  c2 = x0+x1+x2;
  c1 = -(x0*(x[1]+x[2])+x1*(x[0]+x[2])+x2*(x[0]+x[1]));
  c0 = x0*x[1]*x[2]+x1*x[0]*x[2]+x2*x[0]*x[1];

}

//_____________________________________________________________________________
void AliTRDCalibraFit::CalculPolynomeLagrange3(const Double_t *x, const Double_t *y, Double_t &c0, Double_t &c1, Double_t &c2, Double_t &c3, Double_t &c4) const
{
  //
  // Calcul the coefficients of the polynome passant par ces quatre points de degre 3
  //
  Double_t x0 = y[0]/((x[0]-x[1])*(x[0]-x[2])*(x[0]-x[3]));
  Double_t x1 = y[1]/((x[1]-x[0])*(x[1]-x[2])*(x[1]-x[3]));
  Double_t x2 = y[2]/((x[2]-x[0])*(x[2]-x[1])*(x[2]-x[3]));
  Double_t x3 = y[3]/((x[3]-x[0])*(x[3]-x[1])*(x[3]-x[2]));

  c4 = 0.0;
  c3 = x0+x1+x2+x3;
  c2 = -(x0*(x[1]+x[2]+x[3])
	   +x1*(x[0]+x[2]+x[3])
	   +x2*(x[0]+x[1]+x[3])
	   +x3*(x[0]+x[1]+x[2]));
  c1 = (x0*(x[1]*x[2]+x[1]*x[3]+x[2]*x[3])
	  +x1*(x[0]*x[2]+x[0]*x[3]+x[2]*x[3])
	  +x2*(x[0]*x[1]+x[0]*x[3]+x[1]*x[3])
	  +x3*(x[0]*x[1]+x[0]*x[2]+x[1]*x[2]));
  
  c0 = -(x0*x[1]*x[2]*x[3]
	  +x1*x[0]*x[2]*x[3]
	  +x2*x[0]*x[1]*x[3]
	  +x3*x[0]*x[1]*x[2]);  


}

//_____________________________________________________________________________
void AliTRDCalibraFit::CalculPolynomeLagrange4(const Double_t *x, const Double_t *y, Double_t &c0, Double_t &c1, Double_t &c2, Double_t &c3, Double_t &c4) const
{
  //
  // Calcul the coefficients of the polynome passant par ces cinqs points de degre 4
  //
  Double_t x0 = y[0]/((x[0]-x[1])*(x[0]-x[2])*(x[0]-x[3])*(x[0]-x[4]));
  Double_t x1 = y[1]/((x[1]-x[0])*(x[1]-x[2])*(x[1]-x[3])*(x[1]-x[4]));
  Double_t x2 = y[2]/((x[2]-x[0])*(x[2]-x[1])*(x[2]-x[3])*(x[2]-x[4]));
  Double_t x3 = y[3]/((x[3]-x[0])*(x[3]-x[1])*(x[3]-x[2])*(x[3]-x[4]));
  Double_t x4 = y[4]/((x[4]-x[0])*(x[4]-x[1])*(x[4]-x[2])*(x[4]-x[3]));
 

  c4 = x0+x1+x2+x3+x4;
  c3 = -(x0*(x[1]+x[2]+x[3]+x[4])
	   +x1*(x[0]+x[2]+x[3]+x[4])
	   +x2*(x[0]+x[1]+x[3]+x[4])
	   +x3*(x[0]+x[1]+x[2]+x[4])
	   +x4*(x[0]+x[1]+x[2]+x[3]));
  c2 = (x0*(x[1]*x[2]+x[1]*x[3]+x[1]*x[4]+x[2]*x[3]+x[2]*x[4]+x[3]*x[4])
	  +x1*(x[0]*x[2]+x[0]*x[3]+x[0]*x[4]+x[2]*x[3]+x[2]*x[4]+x[3]*x[4])
	  +x2*(x[0]*x[1]+x[0]*x[3]+x[0]*x[4]+x[1]*x[3]+x[1]*x[4]+x[3]*x[4])
	  +x3*(x[0]*x[1]+x[0]*x[2]+x[0]*x[4]+x[1]*x[2]+x[1]*x[4]+x[2]*x[4])
	  +x4*(x[0]*x[1]+x[0]*x[2]+x[0]*x[3]+x[1]*x[2]+x[1]*x[3]+x[2]*x[3]));

  c1 = -(x0*(x[1]*x[2]*x[3]+x[1]*x[2]*x[4]+x[1]*x[3]*x[4]+x[2]*x[3]*x[4])
	  +x1*(x[0]*x[2]*x[3]+x[0]*x[2]*x[4]+x[0]*x[3]*x[4]+x[2]*x[3]*x[4])
	  +x2*(x[0]*x[1]*x[3]+x[0]*x[1]*x[4]+x[0]*x[3]*x[4]+x[1]*x[3]*x[4])
	  +x3*(x[0]*x[1]*x[2]+x[0]*x[1]*x[4]+x[0]*x[2]*x[4]+x[1]*x[2]*x[4])
	  +x4*(x[0]*x[1]*x[2]+x[0]*x[1]*x[3]+x[0]*x[2]*x[3]+x[1]*x[2]*x[3]));

  c0 = (x0*x[1]*x[2]*x[3]*x[4]
	  +x1*x[0]*x[2]*x[3]*x[4]
	  +x2*x[0]*x[1]*x[3]*x[4]
	  +x3*x[0]*x[1]*x[2]*x[4]
	  +x4*x[0]*x[1]*x[2]*x[3]);

}
//_____________________________________________________________________________
void AliTRDCalibraFit::NormierungCharge()
{
  //
  // Normalisation of the gain factor resulting for the fits
  //
  
  // Calcul of the mean of choosen method by fFitChargeNDB
  Double_t sum         = 0.0;
  //printf("total number of entries %d\n",fVectorFitCH->GetEntriesFast());
  for (Int_t k = 0; k < (Int_t) fVectorFit.GetEntriesFast(); k++) {
    Int_t    total    = 0;
    Int_t    detector = ((AliTRDFitInfo *) fVectorFit.At(k))->GetDetector();
    Float_t *coef     = ((AliTRDFitInfo *) fVectorFit.At(k))->GetCoef();
    //printf("detector %d coef[0] %f\n",detector,coef[0]);
    if (GetStack(detector) == 2) {
      total = 1728;
    }
    if (GetStack(detector) != 2) {
      total = 2304;
    }
    for (Int_t j = 0; j < total; j++) {
      if (coef[j] >= 0) {
	sum += coef[j];
      }
    }
  }

  if (sum > 0) {
    fScaleFitFactor = fScaleFitFactor / sum;
  }
  else {
    fScaleFitFactor = 1.0;
  }  

  //methode de boeuf mais bon...
  Double_t scalefactor = fScaleFitFactor;
  
  if(fDebugLevel > 1){
    
    if ( !fDebugStreamer ) {
      //debug stream
      TDirectory *backup = gDirectory;
      fDebugStreamer = new TTreeSRedirector("TRDDebugFitCH.root");
      if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
    } 
    (* fDebugStreamer) << "NormierungCharge"<<
      "scalefactor="<<scalefactor<<
      "\n";  
    }
}
//_____________________________________________________________________________
TH1I *AliTRDCalibraFit::ReBin(const TH1I *hist) const
{
  //
  // Rebin of the 1D histo for the gain calibration if needed.
  // you have to choose fRebin, divider of fNumberBinCharge
  //

 TAxis *xhist  = hist->GetXaxis();
 TH1I  *rehist = new TH1I("projrebin","",(Int_t) xhist->GetNbins()/fRebin
                                        ,xhist->GetBinLowEdge(1)
                                        ,xhist->GetBinUpEdge(xhist->GetNbins()));

 AliInfo(Form("fRebin: %d",fRebin));
 Int_t i = 1;
 for (Int_t k = 1; k <= (Int_t) xhist->GetNbins()/fRebin; k++) {
   Double_t sum = 0.0;
   for (Int_t ji = i; ji < i+fRebin; ji++) {
     sum += hist->GetBinContent(ji);
   }
   sum = sum / fRebin;
   rehist->SetBinContent(k,sum);
   i += fRebin;
 }

 return rehist;

}

//_____________________________________________________________________________
TH1F *AliTRDCalibraFit::ReBin(const TH1F *hist) const
{
  //
  // Rebin of the 1D histo for the gain calibration if needed
  // you have to choose fRebin divider of fNumberBinCharge
  //

  TAxis *xhist  = hist->GetXaxis();
  TH1F  *rehist = new TH1F("projrebin","",(Int_t) xhist->GetNbins()/fRebin
                                         ,xhist->GetBinLowEdge(1)
                                         ,xhist->GetBinUpEdge(xhist->GetNbins()));

  AliInfo(Form("fRebin: %d",fRebin));
  Int_t i = 1;
  for (Int_t k = 1; k <= (Int_t) xhist->GetNbins()/fRebin; k++) {
    Double_t sum = 0.0;
    for (Int_t ji = i; ji < i+fRebin; ji++) {
      sum += hist->GetBinContent(ji);
    }
    sum = sum/fRebin;
    rehist->SetBinContent(k,sum);
    i += fRebin;
  }

  return rehist;
  
}
//
//____________Some basic geometry function_____________________________________
//

//_____________________________________________________________________________
Int_t AliTRDCalibraFit::GetLayer(Int_t d) const
{
  //
  // Reconstruct the plane number from the detector number
  //

  return ((Int_t) (d % 6));

}

//_____________________________________________________________________________
Int_t AliTRDCalibraFit::GetStack(Int_t d) const
{
  //
  // Reconstruct the stack number from the detector number
  //
  const Int_t kNlayer = 6;

  return ((Int_t) (d % 30) / kNlayer);

}

//_____________________________________________________________________________
Int_t AliTRDCalibraFit::GetSector(Int_t d) const
{
  //
  // Reconstruct the sector number from the detector number
  //
  Int_t fg = 30;

  return ((Int_t) (d / fg));

}

//
//____________Fill and Init tree Gain, PRF, Vdrift and T0______________________
//
//_______________________________________________________________________________
void AliTRDCalibraFit::ResetVectorFit()
{
  //
  // Reset the VectorFits
  //

  fVectorFit.SetOwner();
  fVectorFit.Clear();
  fVectorFit2.SetOwner();
  fVectorFit2.Clear();
  
}
//
//____________Private Functions________________________________________________
//

//_____________________________________________________________________________
Double_t AliTRDCalibraFit::PH(const Double_t *x, const Double_t *par) 
{
  //
  // Function for the fit
  //

  //TF1 *fAsymmGauss = new TF1("fAsymmGauss",AsymmGauss,0,4,6);

  //PARAMETERS FOR FIT PH
  // PASAv.4
  //fAsymmGauss->SetParameter(0,0.113755);
  //fAsymmGauss->SetParameter(1,0.350706);
  //fAsymmGauss->SetParameter(2,0.0604244);
  //fAsymmGauss->SetParameter(3,7.65596);
  //fAsymmGauss->SetParameter(4,1.00124);
  //fAsymmGauss->SetParameter(5,0.870597);  // No tail cancelation

  Double_t xx = x[0];
  
  if (xx < par[1]) {
    return par[5];
  }

  Double_t dx       = 0.005;
  Double_t xs       = par[1];
  Double_t ss       = 0.0;
  Double_t paras[2] = { 0.0, 0.0 };

  while (xs < xx) {
    if ((xs >= par[1]) &&
        (xs < (par[1]+par[2]))) {
      //fAsymmGauss->SetParameter(0,par[0]);
      //fAsymmGauss->SetParameter(1,xs);
      //ss += fAsymmGauss->Eval(xx);
      paras[0] = par[0];
      paras[1] = xs;
      ss += AsymmGauss(&xx,paras);
    }
    if ((xs >= (par[1]+par[2])) && 
        (xs <  (par[1]+par[2]+par[3]))) {
      //fAsymmGauss->SetParameter(0,par[0]*par[4]);
      //fAsymmGauss->SetParameter(1,xs);
      //ss += fAsymmGauss->Eval(xx);
      paras[0] = par[0]*par[4];
      paras[1] = xs;
      ss += AsymmGauss(&xx,paras);
    }
    xs += dx;
  }
  
  return ss + par[5];

}

//_____________________________________________________________________________
Double_t AliTRDCalibraFit::AsymmGauss(const Double_t *x, const Double_t *par)
{
  //
  // Function for the fit
  //

  //par[0] = normalization
  //par[1] = mean
  //par[2] = sigma
  //norm0  = 1
  //par[3] = lambda0
  //par[4] = norm1
  //par[5] = lambda1
  
  Double_t par1save = par[1];    
  //Double_t par2save = par[2];
  Double_t par2save = 0.0604244;
  //Double_t par3save = par[3];
  Double_t par3save = 7.65596;
  //Double_t par5save = par[5];
  Double_t par5save = 0.870597;
  Double_t dx       = x[0] - par1save;

  Double_t  sigma2  = par2save*par2save;
  Double_t  sqrt2   = TMath::Sqrt(2.0);
  Double_t  exp1    = par3save * TMath::Exp(-par3save * (dx - 0.5 * par3save * sigma2))
                               * (1.0 - AliMathBase::ErfFast((par3save * sigma2 - dx) / (sqrt2 * par2save)));
  Double_t  exp2    = par5save * TMath::Exp(-par5save * (dx - 0.5 * par5save * sigma2))
                               * (1.0 - AliMathBase::ErfFast((par5save * sigma2 - dx) / (sqrt2 * par2save)));

  //return par[0]*(exp1+par[4]*exp2);
  return par[0] * (exp1 + 1.00124 * exp2);

}

//_____________________________________________________________________________
Double_t AliTRDCalibraFit::FuncLandauGaus(const Double_t *x, const Double_t *par)
{
  //
  // Sum Landau + Gaus with identical mean
  //

  Double_t valLandau = par[0] * TMath::Landau(x[0],par[1],par[2]);
  //Double_t valGaus   = par[3] * TMath::Gaus(x[0],par[4],par[5]);
  Double_t valGaus   = par[3] * TMath::Gaus(x[0],par[1],par[4]);
  Double_t val       = valLandau + valGaus;

  return val;

}

//_____________________________________________________________________________
Double_t AliTRDCalibraFit::LanGauFun(const Double_t *x, const Double_t *par) 
{
  //
  // Function for the fit
  //
  // Fit parameters:
  // par[0]=Width (scale) parameter of Landau density
  // par[1]=Most Probable (MP, location) parameter of Landau density
  // par[2]=Total area (integral -inf to inf, normalization constant)
  // par[3]=Width (sigma) of convoluted Gaussian function
  //
  // In the Landau distribution (represented by the CERNLIB approximation), 
  // the maximum is located at x=-0.22278298 with the location parameter=0.
  // This shift is corrected within this function, so that the actual
  // maximum is identical to the MP parameter.
  //  

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np       = 100.0;             // Number of convolution steps
  Double_t sc       =   5.0;             // Convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow;
  Double_t xupp;
  Double_t step;
  Double_t i;
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp - xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for (i = 1.0; i <= np/2; i++) {

    xx    = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum  += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx    = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum  += fland * TMath::Gaus(x[0],xx,par[3]);

  }

  return (par[2] * step * sum * invsq2pi / par[3]);

}
//_____________________________________________________________________________
TF1 *AliTRDCalibraFit::LanGauFit(TH1 *his, const Double_t *fitrange, const Double_t *startvalues
                                      , const Double_t *parlimitslo, const Double_t *parlimitshi
                                      , Double_t *fitparams, Double_t *fiterrors
                                      , Double_t *chiSqr, Int_t *ndf) const
{
  //
  // Function for the fit
  //
  
  Int_t i;
  Char_t funname[100];
  
  TF1 *ffitold = (TF1 *) gROOT->GetListOfFunctions()->FindObject(funname);
  if (ffitold) {
    delete ffitold;
  }  

  TF1 *ffit    = new TF1(funname,LanGauFun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MP","Area","GSigma");
  
  for (i = 0; i < 4; i++) {
    ffit->SetParLimits(i,parlimitslo[i],parlimitshi[i]);
  }
  
  his->Fit(funname,"RB0");                   // Fit within specified range, use ParLimits, do not plot
  
  ffit->GetParameters(fitparams);            // Obtain fit parameters
  for (i = 0; i < 4; i++) {
    fiterrors[i] = ffit->GetParError(i);     // Obtain fit parameter errors
  }
  chiSqr[0] = ffit->GetChisquare();          // Obtain chi^2
  ndf[0]    = ffit->GetNDF();                // Obtain ndf

  return (ffit);                             // Return fit function
   
}

//_____________________________________________________________________________
Int_t AliTRDCalibraFit::LanGauPro(const Double_t *params, Double_t &maxx, Double_t &fwhm) 
{
  //
  // Function for the fit
  //

  Double_t p;
  Double_t x;
  Double_t fy;
  Double_t fxr;
  Double_t fxl;
  Double_t step;
  Double_t l;
  Double_t lold;

  Int_t    i        = 0;
  Int_t    maxcalls = 10000;
  
  // Search for maximum
  p    = params[1] - 0.1 * params[0];
  step = 0.05 * params[0];
  lold = -2.0;
  l    = -1.0;
  
  while ((l != lold) && (i < maxcalls)) {
    i++;
    lold = l;
    x    = p + step;
    l    = LanGauFun(&x,params);
    if (l < lold) {
      step = -step / 10.0;
    }
    p += step;
  }
  
  if (i == maxcalls) {
    return (-1);
  }
  maxx = x;
  fy = l / 2.0;

  // Search for right x location of fy  
  p    = maxx + params[0];
  step = params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;
  
  while ( (l != lold) && (i < maxcalls) ) {
    i++;
    
    lold = l;
    x = p + step;
    l = TMath::Abs(LanGauFun(&x,params) - fy);
    
    if (l > lold)
      step = -step/10;
 
    p += step;
  }
  
  if (i == maxcalls)
    return (-2);
  
  fxr = x;
  
  
  // Search for left x location of fy
  
  p = maxx - 0.5 * params[0];
  step = -params[0];
  lold = -2.0;
  l    = -1.0e300;
  i    = 0;
  
  while ((l != lold) && (i < maxcalls)) {
    i++;
    lold = l;
    x    = p + step;
    l    = TMath::Abs(LanGauFun(&x,params) - fy);
    if (l > lold) {
      step = -step / 10.0;
    }
    p += step;
  }
  
  if (i == maxcalls) {
    return (-3);
  }

  fxl  = x;
  fwhm = fxr - fxl;

  return (0);
}
//_____________________________________________________________________________
Double_t AliTRDCalibraFit::GausConstant(const Double_t *x, const Double_t *par)
{
  //
  // Gaus with identical mean
  //

  Double_t gauss   = par[0] * TMath::Gaus(x[0],0.0,par[1])+par[2];
 
  return gauss;

}



