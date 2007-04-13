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
// _fDebug = 1: a comparaison of the coefficients found and the default values 
//              in the choosen database.
//              fCoef , histogram of the coefs as function of the calibration group number
//              fDelta , histogram of the relative difference of the coef with the default
//                        value in the database as function of the calibration group number
//              fError , dirstribution of this relative difference
// _fDebug = 2: only the fit of the choosen calibration group fFitVoir (SetFitVoir)
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

#include <TTree.h>
#include <TLine.h>
#include <TH1I.h>
#include <TStyle.h>
#include <TProfile2D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TAxis.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TROOT.h>

#include "AliLog.h"
#include "AliCDBManager.h"

#include "AliTRDCalibraFit.h"
#include "AliTRDCalibraMode.h"
#include "AliTRDCalibraVector.h"
#include "AliTRDcalibDB.h"
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
  ,fWriteNameCoef(0)
  ,fFitPHOn(kFALSE)
  ,fFitPol2On(kFALSE)
  ,fFitLagrPolOn(kFALSE)
  ,fTakeTheMaxPH(kFALSE)
  ,fFitPHPeriode(0)
  ,fFitPHNDB(1)
  ,fBeginFitCharge(0.0)
  ,fT0Shift(0.0)
  ,fRangeFitPRF(0.0)
  ,fFitPRFOn(kFALSE)
  ,fRMSPRFOn(kFALSE)
  ,fFitPRFNDB(0)
  ,fMeanChargeOn(kFALSE)
  ,fFitChargeBisOn(kFALSE)
  ,fFitChargeOn(kFALSE)
  ,fFitMeanWOn(kFALSE)
  ,fFitChargeNDB(0)
  ,fAccCDB(kFALSE)
  ,fMinEntries(0)
  ,fRebin(0)
  ,fNumberFit(0)
  ,fNumberFitSuccess(0)
  ,fNumberEnt(0)
  ,fStatisticMean(0.0)
  ,fDebug(0)
  ,fFitVoir(0)
  ,fCalibraMode(0)
  ,fPRF(0)
  ,fGain(0)
  ,fT0(0)
  ,fVdrift(0)
  ,fVdriftDetector(0)
  ,fVdriftPad(0x0)
  ,fT0Detector(0)
  ,fT0Pad(0x0)
  ,fPRFDetector(0)
  ,fPRFPad(0x0)
  ,fCoefCH(0x0)
  ,fScaleFitFactor(0.0)
  ,fEntriesCurrent(0)
  ,fCalibraVector(0)
  ,fVectorFitCH(0) 
{
  //
  // Default constructor
  //

  fCalibraMode = new AliTRDCalibraMode();

  // Write
  for (Int_t i = 0; i < 3; i++) {
    fWriteCoef[i] = kFALSE;
  }

  // Debug Mode
  for (Int_t k = 0; k < 3; k++) {
    fDet[k] = 0;
  }

  for (Int_t i = 0; i < 3; i++) {
    fPhd[i] = 0.0;
  }

  // Init
  Init();
  
}

//______________________________________________________________________________________
AliTRDCalibraFit::AliTRDCalibraFit(const AliTRDCalibraFit &c)
  :TObject(c)
  ,fWriteNameCoef(0)
  ,fFitPHOn(kFALSE)
  ,fFitPol2On(kFALSE)
  ,fFitLagrPolOn(kFALSE)
  ,fTakeTheMaxPH(kFALSE)
  ,fFitPHPeriode(0)
  ,fFitPHNDB(1)
  ,fBeginFitCharge(0.0)
  ,fT0Shift(0.0)
  ,fRangeFitPRF(0.0)
  ,fFitPRFOn(kFALSE)
  ,fRMSPRFOn(kFALSE)
  ,fFitPRFNDB(0)
  ,fMeanChargeOn(kFALSE)
  ,fFitChargeBisOn(kFALSE)
  ,fFitChargeOn(kFALSE)
  ,fFitMeanWOn(kFALSE)
  ,fFitChargeNDB(0)
  ,fAccCDB(kFALSE)
  ,fMinEntries(0)
  ,fRebin(0) 
  ,fNumberFit(0)
  ,fNumberFitSuccess(0)
  ,fNumberEnt(0)
  ,fStatisticMean(0.0)
  ,fDebug(0)
  ,fFitVoir(0)
  ,fCalibraMode(0)
  ,fPRF(0)
  ,fGain(0)
  ,fT0(0)
  ,fVdrift(0)
  ,fVdriftDetector(0)
  ,fVdriftPad(0x0)
  ,fT0Detector(0)
  ,fT0Pad(0x0)
  ,fPRFDetector(0)
  ,fPRFPad(0x0)
  ,fCoefCH(0x0)
  ,fScaleFitFactor(0.0)
  ,fEntriesCurrent(0)
  ,fCalibraVector(0)
  ,fVectorFitCH(0)
{
  //
  // Copy constructor
  //

}

//____________________________________________________________________________________
AliTRDCalibraFit::~AliTRDCalibraFit()
{
  //
  // AliTRDCalibraFit destructor
  //

  ClearTree();

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
void AliTRDCalibraFit::ClearTree() 
{
  //
  // Delete the trees
  //
  
  if (fPRF) {
    delete fPRF;
    fPRF    = 0x0;
  }
  if (fGain) {
    delete fGain;
    fGain   = 0x0;
  }
  if (fT0) {
    delete fT0;
    fT0     = 0x0;
  }
  if (fVdrift) {
    delete fVdrift;
    fVdrift = 0x0;
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::Init() 
{
  //
  // Init some default values
  //

  // Write
  fWriteNameCoef        = "TRD.coefficient.root";
  
  // Fit
  fFitPHPeriode         = 1;
  fBeginFitCharge       = 3.5;
  fRangeFitPRF          = 1.0;
  fMinEntries           = 800;
  fT0Shift              = 0.126256;
  
  // Internal variables
  
  // Variables in the loop
  for (Int_t k = 0; k < 4; k++) {
    fChargeCoef[k] = 1.0;
    fVdriftCoef[k] = 1.5;
    fT0Coef[k]     = -1.0;
  }
  fChargeCoef[4] = 1.0;
  for (Int_t i = 0; i < 3; i++) {
    fPRFCoef[i]    = -1.0;
  }
    
  // Local database to be changed
  fRebin = 1;
  
}

//____________Functions fit Online CH2d________________________________________
Bool_t AliTRDCalibraFit::FitCHOnline(TH2I *ch)
{
  //
  // Fit the 1D histos, projections of the 2D ch on the Xaxis, for each
  // calibration group normalized the resulted coefficients (to 1 normally)
  // and write the results in a tree
  //

  //A small check
  if((fFitChargeNDB == 0) && (!fFitChargeOn)){
    AliInfo("You have choosen to write the default fit method but it is not on!");
    return kFALSE;
  }
  if((fFitChargeNDB == 1) && (!fMeanChargeOn)){
    AliInfo("You have choosen to write the mean method but it is not on!");
    return kFALSE;
  }
  if((fFitChargeNDB == 2) && (!fFitChargeBisOn)){
    AliInfo("You have choosen to write the second fit method but it is not on!");
    return kFALSE;
  }
  if((fFitChargeNDB == 4) && (!fFitMeanWOn)){
    AliInfo("You have choosen to write the mean w method but it is not on!");
    return kFALSE;
  }
 
   
  // Number of Xbins (detectors or groups of pads)
  TAxis   *xch     = ch->GetXaxis();
  Int_t    nbins   = xch->GetNbins();
  TAxis   *yph     = ch->GetYaxis();
  Int_t    nybins  = yph->GetNbins();
  if (!InitFit(nbins,0)) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;

  // Init fCountDet and fCount
  InitfCountDetAndfCount(0);
  
  // Beginning of the loop betwwen dect1 and dect2
  for (Int_t idect = fDect1[0]; idect < fDect2[0]; idect++) {

    TH1I *projch = (TH1I *) ch->ProjectionY("projch",idect+1,idect+1,(Option_t *)"e");
    projch->SetDirectory(0);

    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    UpdatefCountDetAndfCount(idect,0);
    
    // Reconstruction of the row and pad group: rowmin, row max ...
    ReconstructFitRowMinRowMax(idect, 0);
    
    // Number of entries for this calibration group
    Double_t nentries = 0.0;
    Double_t mean = 0.0;
    for (Int_t k = 0; k < nybins; k++) {
      nentries += ch->GetBinContent(ch->GetBin(idect+1,k+1));
      mean += projch->GetBinCenter(k+1)*projch->GetBinContent(k+1);
    }
    if (nentries > 0) {
      fNumberEnt++;
      mean /= nentries;
    }

   
    // Rebin and statistic stuff
    // Rebin
    if (fRebin > 1) {
      projch = ReBin((TH1I *) projch);
    }
    // This detector has not enough statistics or was off
    if (nentries < fMinEntries) {
      // Fill with the default infos
      NotEnoughStatistic(idect,0);
      // Memory!!!
      if (fDebug != 2) {
	delete projch;
      }
      continue;
    }
    
    // Statistics of the group fitted
    AliInfo(Form("For the group number %d there are %f stats",idect,nentries));
    fStatisticMean += nentries;
    fNumberFit++;
  
    
    // Method Mean and fit
    // idect is egal for fDebug = 0 and 2, only to fill the hist 
    fChargeCoef[1] = mean;
    if(fMeanChargeOn){
      FitMean((TH1 *) projch,(Int_t) (idect-fDect1[0]),nentries);
    }
    if(fFitChargeOn){
      FitCH((TH1 *) projch,(Int_t) (idect-fDect1[0]));
    }
    if(fFitChargeBisOn) { 
      FitBisCH((TH1 *) projch,(Int_t) (idect-fDect1[0]));
    }
    if(fFitMeanWOn){
      FitMeanW((TH1 *) projch,(Int_t) (idect-fDect1[0]));
    }
    
    // Visualise the detector for fDebug 3 or 4
    // Here is the reconstruction of the pad and row group is used!
    if (fDebug >= 3) {
      FillCoefChargeDB();
    }
    // Fill Infos Fit
    FillInfosFit(idect,0);
    
    // Memory!!!
    if (fDebug != 2) {
      delete projch;
    }
    
        
  } // Boucle object

  
  // Normierungcharge
  if (fDebug != 2) {
    NormierungCharge();
  }
  
  // Plot
  // 0 no plot, 1 and 4 error plot, 3 and 4 DB plot
  if ((fDebug == 1) || 
      (fDebug == 4)) {
    PlotWriteCH();
  }
  if ((fDebug == 4) || 
      (fDebug == 3)) {
    PlotCHDB();
  }

  // Mean Statistic
  if (fNumberFit > 0) {
    AliInfo(Form("There are %d with at least one entries.",fNumberEnt));
    AliInfo(Form("%d fits have been proceeded (sucessfully or not...).",fNumberFit));
    AliInfo(Form("There is a mean statistic of: %d over these fitted histograms and %d successfulled fits"
                , (Int_t) fStatisticMean/fNumberFit, fNumberFitSuccess));
    fStatisticMean = fStatisticMean / fNumberFit;
  }
  else {
    AliInfo(Form("There are %d with at least one entries. There is no fit!",fNumberEnt));
  }
  
  // Write the things!
  ConvertVectorFitCHTree();
  if (fWriteCoef[0]) {
    WriteFitInfos(0);       
  }
  
  return kTRUE;
  
}

//____________Functions fit Online CH2d________________________________________
Bool_t AliTRDCalibraFit::FitCHOnline()
{
  //
  // Reconstruct a 1D histo from the vectorCH for each calibration group,
  // fit the histo, normalized the resulted coefficients (to 1 normally)
  // and write the results in a tree
  //

  //A small check
  if((fFitChargeNDB == 0) && (!fFitChargeOn)){
    AliInfo("You have choosen to write the default fit method but it is not on!");
    return kFALSE;
  }
  if((fFitChargeNDB == 1) && (!fMeanChargeOn)){
    AliInfo("You have choosen to write the mean method but it is not on!");
    return kFALSE;
  }
  if((fFitChargeNDB == 2) && (!fFitChargeBisOn)){
    AliInfo("You have choosen to write the second fit method but it is not on!");
    return kFALSE;
  }
  if((fFitChargeNDB == 4) && (!fFitMeanWOn)){
    AliInfo("You have choosen to write the mean w method but it is not on!");
    return kFALSE;
  }
 

  //Warning
  if (!fCalibraVector) {
    AliError("You have first to set the calibravector before using this function!");
    return kFALSE;
  }

  // Number of Xbins (detectors or groups of pads)
  if (!InitFit(0,0)) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;
 
  // Init fCountDet and fCount
  InitfCountDetAndfCount(0);

  // Beginning of the loop between dect1 and dect2
  for (Int_t idect = fDect1[0]; idect < fDect2[0]; idect++) {
    
    // Search if the group is in the VectorCH
    Int_t place = fCalibraVector->SearchInVector(idect,0);
       
    // Is in
    TH1F *projch = 0x0;
    TString name("CH");
    name += idect;
    if (place != -1) {
      projch  = fCalibraVector->ConvertVectorCTHisto(place,(const char *) name);
      projch->SetDirectory(0);
    }

    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    UpdatefCountDetAndfCount(idect,0);
    
    // Reconstruction of the row and pad group: rowmin, row max ...
    ReconstructFitRowMinRowMax(idect,0);
    
    // Number of entries and mean
    Double_t nentries = 0.0;
    Double_t mean = 0.0;
    if (projch) {
      for (Int_t k = 0; k < fCalibraVector->GetNumberBinCharge(); k++) {
        nentries += projch->GetBinContent(k+1);
	mean += projch->GetBinCenter(k+1)*projch->GetBinContent(k+1);
      }
    }
    if (nentries > 0) {
      fNumberEnt++;
      mean /= nentries;
    }

    // Rebin and statistic stuff
    // Rebin
    if ((fRebin >  1) && 
        (place != -1)) {
      projch = ReBin((TH1F *) projch);
    }

    // This detector has not enough statistics or was not found in VectorCH
    if ((place == -1) || 
        ((place != -1) && 
         (nentries < fMinEntries))) {
      
      // Fill with the default infos
      NotEnoughStatistic(idect,0);

      // Memory!!!
      if (fDebug != 2) {
	delete projch;
      }     
      
      continue;

    }

    // Statistic of the histos fitted
    AliInfo(Form("For the group number %d there are %f stats",idect,nentries));
    fStatisticMean += nentries;
    fNumberFit++;

    
    // Method Mean and fit
    // idect is egal for fDebug = 0 and 2, only to fill the hist 
    fChargeCoef[1] = mean;
    if(fMeanChargeOn){
      FitMean((TH1 *) projch,(Int_t) (idect-fDect1[0]),nentries);
    }
    if(fFitChargeOn){
      FitCH((TH1 *) projch,(Int_t) (idect-fDect1[0]));
    }
    if(fFitChargeBisOn) { 
      FitBisCH((TH1 *) projch,(Int_t) (idect-fDect1[0]));
    }
    if(fFitMeanWOn){
      FitMeanW((TH1 *) projch,(Int_t) (idect-fDect1[0]));
    }
    
    // Visualise the detector for fDebug 3 or 4
    // Here is the reconstruction of the pad and row group is used!
    if (fDebug >= 3) {
      FillCoefChargeDB();
    }
    
    // Fill Infos Fit
    FillInfosFit(idect,0); 
    
    // Memory!!!
    if (fDebug != 2) {
      delete projch;
    }
    
  } // Boucle object

  
  // Normierungcharge
  if (fDebug != 2) {
    NormierungCharge();
  }

  
  // Plot
  // 0 no plot, 1 and 4 error plot, 3 and 4 DB plot
  if ((fDebug == 1) || 
      (fDebug == 4)){
    PlotWriteCH();
  }
  if((fDebug == 4) || 
     (fDebug == 3)){
    PlotCHDB();
  }
  
  // Mean Statistics
  if (fNumberFit > 0) {
    AliInfo(Form("There are %d with at least one entries.",fNumberEnt));
    AliInfo(Form("%d fits have been proceeded (sucessfully or not...).",fNumberFit));
    AliInfo(Form("There is a mean statistic of: %d over these fitted histograms and %d successfulled fits"
                ,(Int_t) fStatisticMean/fNumberFit, fNumberFitSuccess));
    fStatisticMean = fStatisticMean / fNumberFit;
  }
  else {
    AliInfo(Form("There are %d with at least one entries. There is no fit!",fNumberEnt));
  }
  
  // Write the things!
  ConvertVectorFitCHTree();
  if (fWriteCoef[0]) {
    WriteFitInfos(0);      
  }
  
  return kTRUE;
  
}

//____________Functions fit Online CH2d________________________________________
Bool_t AliTRDCalibraFit::FitCHOnline(TTree *tree)
{
  //
  // Look if the calibration group can be found in the tree, if yes take the
  // histo, fit it, normalized the resulted coefficients (to 1 normally) and
  // write the results in a tree
  //

  //A small check
  if((fFitChargeNDB == 0) && (!fFitChargeOn)){
    AliInfo("You have choosen to write the default fit method but it is not on!");
    return kFALSE;
  }
  if((fFitChargeNDB == 1) && (!fMeanChargeOn)){
    AliInfo("You have choosen to write the mean method but it is not on!");
    return kFALSE;
  }
  if((fFitChargeNDB == 2) && (!fFitChargeBisOn)){
    AliInfo("You have choosen to write the second fit method but it is not on!");
    return kFALSE;
  }
  if((fFitChargeNDB == 4) && (!fFitMeanWOn)){
    AliInfo("You have choosen to write the mean w method but it is not on!");
    return kFALSE;
  }
 
   
  // Number of Xbins (detectors or groups of pads)
  if (!InitFit(0,0)) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;
  
  // Initialise
  fCalibraVector = new AliTRDCalibraVector();
 
   
  // Init fCountDet and fCount
  InitfCountDetAndfCount(0);
  TH1F      *projch      = 0x0;
  tree->SetBranchAddress("histo",&projch);
  TObjArray *vectorplace = fCalibraVector->ConvertTreeVector(tree);

  // Beginning of the loop between dect1 and dect2
  for (Int_t idect = fDect1[0]; idect < fDect2[0]; idect++) {
    
    //Search if the group is in the VectorCH
    Int_t place = fCalibraVector->SearchInTreeVector(vectorplace,idect);
    
    // Is in
    if (place != -1) {
      // Variable
      tree->GetEntry(place);
    }
    
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    UpdatefCountDetAndfCount(idect,0);
    
    // Reconstruction of the row and pad group: rowmin, row max ...
    ReconstructFitRowMinRowMax(idect,0);
    
    // Number of entries and mean
    Double_t nentries = 0.0;
    Double_t mean = 0.0;
    if (projch) {
      for (Int_t k = 0; k < projch->GetXaxis()->GetNbins(); k++) {
        nentries += projch->GetBinContent(k+1);
	mean += projch->GetBinCenter(k+1)*projch->GetBinContent(k+1);
      }
    } 
    if (nentries > 0) {
      fNumberEnt++;   
      mean /= nentries;
    }

        
    // Rebin and statistic stuff
    // Rebin
    if ((fRebin  >  1) && 
        (place  != -1)) {
      projch = ReBin((TH1F *) projch);
    }

    // This detector has not enough statistics or was not found in VectorCH
    if((place == -1) || 
       ((place    !=          -1) && 
        (nentries  < fMinEntries))) {
      
      // Fill with the default infos
      NotEnoughStatistic(idect,0);
      
      continue;

    }
    
    // Statistics of the group fitted
    AliInfo(Form("For the group number %d there are %f stats",idect,nentries));
    fNumberFit++;
    fStatisticMean += nentries;
   
    // Method Mean and fit
    // idect is egal for fDebug = 0 and 2, only to fill the hist 
    fChargeCoef[1] = mean;
    if(fMeanChargeOn){
      FitMean((TH1 *) projch,(Int_t) (idect-fDect1[0]),nentries);
    }
    if(fFitChargeOn){
      FitCH((TH1 *) projch,(Int_t) (idect-fDect1[0]));
    }
    if(fFitChargeBisOn) { 
      FitBisCH((TH1 *) projch,(Int_t) (idect-fDect1[0]));
    }
    if(fFitMeanWOn){
      FitMeanW((TH1 *) projch,(Int_t) (idect-fDect1[0]));
    }

    // Visualise the detector for fDebug 3 or 4
    // Here is the reconstruction of the pad and row group is used!
    if (fDebug >= 3) {
      FillCoefChargeDB();
    }    

    // Fill Infos Fit
    FillInfosFit(idect,0); 
        
  } // Boucle object


  // Normierungcharge
  if (fDebug != 2) {
    NormierungCharge();
  }
  
 
  // Plot
  // 0 no plot, 1 and 4 error plot, 3 and 4 DB plot
  if ((fDebug == 1) || 
      (fDebug == 4)){
    PlotWriteCH();
  }
  if ((fDebug == 4) || 
      (fDebug == 3)){
    PlotCHDB();
  }

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

  // Write the things!
  ConvertVectorFitCHTree();
  if (fWriteCoef[0]) {
    WriteFitInfos(0);      
  }
  
  
  return kTRUE;
  
}

//________________functions fit Online PH2d____________________________________
Bool_t AliTRDCalibraFit::FitPHOnline(TProfile2D *ph)
{
  //
  // Take the 1D profiles (average pulse height), projections of the 2D PH
  // on the Xaxis, for each calibration group
  // Fit or use the slope of the average pulse height to reconstruct the
  // drift velocity write the results in a tree
  // A first calibration of T0 is also made  using the same method (slope method)
  //

  //A small check
  if((fFitPHNDB == 0) && (!fFitPHOn)){
    AliInfo("You have choosen to write the fit method but it is not on!");
    return kFALSE;
  }
  if((fFitPHNDB == 1) && (!fFitPol2On)){
    AliInfo("You have choosen to write the Pol2 method but it is not on!");
    return kFALSE;
  }
  if((fFitPHNDB == 3) && (!fFitLagrPolOn)){
    AliInfo("You have choosen to write the LagrPol2 method but it is not on!");
    return kFALSE;
  }
  
  // Number of Xbins (detectors or groups of pads)
  TAxis   *xph     = ph->GetXaxis();
  TAxis   *yph     = ph->GetYaxis();
  Int_t    nbins   = xph->GetNbins();
  Int_t    nybins  = yph->GetNbins();
  if (!InitFit(nbins,1)) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;

  // Init fCountDet and fCount
  InitfCountDetAndfCount(1);
  
  // Beginning of the loop
  for (Int_t idect = fDect1[1]; idect < fDect2[1]; idect++) {
      
    TH1D *projph = (TH1D *) ph->ProjectionY("projph",idect+1,idect+1,(Option_t *) "e");
    projph->SetDirectory(0); 

    // Number of entries for this calibration group
    Double_t nentries = 0;
    for (Int_t k = 0; k < nybins; k++) {
      nentries += ph->GetBinEntries(ph->GetBin(idect+1,k+1));
    }
    if (nentries > 0) {
      fNumberEnt++;
    }  

    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    UpdatefCountDetAndfCount(idect,1);

    // Reconstruction of the row and pad group: rowmin, row max ...
    ReconstructFitRowMinRowMax(idect,1);
    
    // Rebin and statistic stuff
    // This detector has not enough statistics or was off
    if (nentries  < fMinEntries) {
      
      // Fill with the default values
      NotEnoughStatistic(idect,1);     
      
      // Memory!!!
      if (fDebug != 2) {
	delete projph;
      }
           
      continue;

    }
    
    // Statistics of the histos fitted
    AliInfo(Form("For the group number %d there are %f stats",idect,nentries));
    fNumberFit++;
    fStatisticMean += nentries;
    
    // Calcul of "real" coef
    CalculVdriftCoefMean(fCountDet[1],(Int_t) (idect - fDect1[1]));
    CalculT0CoefMean(fCountDet[1],(Int_t) (idect - fDect1[1]));
    
    // Method Mean and fit
    // idect is egal for fDebug = 0 and 2, only to fill the hist 
    if(fFitPol2On){
      FitPente((TH1 *) projph,(Int_t) (idect - fDect1[1]));
    }
    if(fFitLagrPolOn){
      FitLagrangePoly((TH1 *) projph,(Int_t) (idect - fDect1[1]));
    }
    if (fFitPHOn) { 
      FitPH((TH1 *) projph,(Int_t) (idect - fDect1[1]));
    }
     
    // Visualise the detector for fDebug 3 or 4
    // Here is the reconstruction of the pad and row group is used!
    if (fDebug >= 3) {
      FillCoefVdriftDB();
      FillCoefT0DB();
    }
    // Fill the tree if end of a detector or only the pointer to the branch!!!
    FillInfosFit(idect,1);
    // Memory!!!
    if (fDebug != 2) {
      delete projph;
    }
    
  } // Boucle object

  
  // Plot
  // 0 no plot, 1 and 4 error plot, 3 and 4 DB plot
  if ((fDebug == 1) || 
      (fDebug == 4)) {
    PlotWritePH();
    PlotWriteT0();
  }
  if ((fDebug == 4) || 
      (fDebug == 3)) {
    PlotPHDB();
    PlotT0DB();
  }

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

  // Write the things!
  if(fWriteCoef[1]) {
    WriteFitInfos(1);
  }

  return kTRUE;
  
}

//____________Functions fit Online PH2d________________________________________
Bool_t AliTRDCalibraFit::FitPHOnline()
{
  //
  // Reconstruct the average pulse height from the vectorPH for each
  // calibration group
  // Fit or use the slope of the average pulse height to reconstruct the
  // drift velocity write the results in a tree
  // A first calibration of T0 is also made  using the same method (slope method)
  //

  //A small check
  if((fFitPHNDB == 0) && (!fFitPHOn)){
    AliInfo("You have choosen to write the fit method but it is not on!");
    return kFALSE;
  }
  if((fFitPHNDB == 1) && (!fFitPol2On)){
    AliInfo("You have choosen to write the Pol2 method but it is not on!");
    return kFALSE;
  }
  if((fFitPHNDB == 3) && (!fFitLagrPolOn)){
    AliInfo("You have choosen to write the LagrPol2 method but it is not on!");
    return kFALSE;
  }

  //Warning
  if (!fCalibraVector) {
    AliError("You have first to set the calibravector before using this function!");
    return kFALSE;
  }

   
  // Number of Xbins (detectors or groups of pads)
  if (!InitFit(0,1)) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;

  // Init fCountDet and fCount
  InitfCountDetAndfCount(1);

  // Beginning of the loop
  for (Int_t idect = fDect1[1]; idect < fDect2[1]; idect++) {

    // Search if the group is in the VectorCH
    Int_t place = fCalibraVector->SearchInVector(idect,1);
    
    // Is in
    TH1F    *projph = 0x0;
    TString name("PH");
    name += idect;
    if (place != -1) {
      //Entries
      fNumberEnt++;
      projph  = CorrectTheError((TGraphErrors *) (fCalibraVector->ConvertVectorPHisto(place,(const char *) name)));
      projph->SetDirectory(0);
    }

    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    UpdatefCountDetAndfCount(idect,1);
    
    // Reconstruction of the row and pad group: rowmin, row max ...
    ReconstructFitRowMinRowMax(idect,1);

    // Rebin and statistic stuff
    // This detector has not enough statistics or was off
    if ((place == -1) || 
        ((place           !=          -1) && 
         (fEntriesCurrent <  fMinEntries))) {
      
      // Fill with the default values
      NotEnoughStatistic(idect,1);

      // Memory!!!
      if (fDebug != 2) {
	delete projph;
      }

      continue;

    }

    // Statistic of the histos fitted
    AliInfo(Form("For the group number %d there are %d stats",idect,fEntriesCurrent));
    fNumberFit++;
    fStatisticMean += fEntriesCurrent;

    // Calcul of "real" coef
    CalculVdriftCoefMean(fCountDet[1],(Int_t) (idect - fDect1[1]));
    CalculT0CoefMean(fCountDet[1],(Int_t) (idect - fDect1[1]));

    // Method Mean and fit
    // idect is egal for fDebug = 0 and 2, only to fill the hist 
    if(fFitPol2On){
      FitPente((TH1 *) projph,(Int_t) (idect - fDect1[1]));
    }
    if(fFitLagrPolOn){
      FitLagrangePoly((TH1 *) projph,(Int_t) (idect - fDect1[1]));
    }
    if (fFitPHOn) { 
      FitPH((TH1 *) projph,(Int_t) (idect - fDect1[1]));
    }
  
    // Visualise the detector for fDebug 3 or 4
    // Here is the reconstruction of the pad and row group is used!
    if (fDebug >= 3) {
      FillCoefVdriftDB();
      FillCoefT0DB();
    }
    
    // Fill the tree if end of a detector or only the pointer to the branch!!!
    FillInfosFit(idect,1);
    
    // Memory!!!
    if (fDebug != 2) {
      delete projph;
    }
    
  } // Boucle object
  
  
  // Plot
  // 0 no plot, 1 and 4 error plot, 3 and 4 DB plot
  if ((fDebug == 1) || 
      (fDebug == 4)) {
    PlotWritePH();
    PlotWriteT0();
  }
  if ((fDebug == 4) || 
      (fDebug == 3)) {
    PlotPHDB();
    PlotT0DB();
  }
  
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
  
  // Write the things!
  if (fWriteCoef[1]) {
    WriteFitInfos(1);
  }

  return kTRUE;
  
}

//____________Functions fit Online PH2d________________________________________
Bool_t AliTRDCalibraFit::FitPHOnline(TTree *tree)
{
  //
  // Look if the calibration group can be found in the tree, if yes take the
  // histo, fit it, and write the results in a tree
  // A first calibration of T0 is also made  using the same method (slope method)
  //
   
  //A small check
  if ((fFitPHNDB == 0) && (!fFitPHOn)){
    AliInfo("You have choosen to write the fit method but it is not on!");
    return kFALSE;
  }
  if ((fFitPHNDB == 1) && (!fFitPol2On)){
    AliInfo("You have choosen to write the Pol2 method but it is not on!");
    return kFALSE;
  }
  if ((fFitPHNDB == 3) && (!fFitLagrPolOn)){
    AliInfo("You have choosen to write the LagrPol2 method but it is not on!");
    return kFALSE;
  }

  // Number of Xbins (detectors or groups of pads)
  if (!InitFit(0,1)) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;   
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;

  // Initialise
  fCalibraVector = new AliTRDCalibraVector();

  // Init fCountDet and fCount
  InitfCountDetAndfCount(1);
  TGraphErrors *projphtree  = 0x0;
  tree->SetBranchAddress("histo",&projphtree);
  TObjArray    *vectorplace = fCalibraVector->ConvertTreeVector(tree);
  
  // Beginning of the loop
  for (Int_t idect = fDect1[1]; idect < fDect2[1]; idect++) {

    // Search if the group is in the VectorCH
    Int_t place = fCalibraVector->SearchInTreeVector(vectorplace,idect);
    
    TH1F *projph = 0x0;
    // Is in
    if (place != -1) {
      //Entries
      fNumberEnt++;
      // Variable
      tree->GetEntry(place);
      projph = CorrectTheError(projphtree);
    }
    
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    UpdatefCountDetAndfCount(idect,1);

    // Reconstruction of the row and pad group: rowmin, row max ...
    ReconstructFitRowMinRowMax(idect,1);

    // Rebin and statistic stuff
    // This detector has not enough statistics or was off
    if((place == -1) || 
       ((place           !=         -1) && 
        (fEntriesCurrent < fMinEntries))) {
      
      // Fill with the default values
      NotEnoughStatistic(idect,1);
      
      // Memory!!!
      if (fDebug != 2) {
	delete projph;
      }
      
      continue;

    }
    
    // Statistics of the histos fitted
    AliInfo(Form("For the group number %d there are %d stats",idect,fEntriesCurrent));
    fNumberFit++;
    fStatisticMean += fEntriesCurrent;

    // Calcul of "real" coef
    CalculVdriftCoefMean(fCountDet[1],(Int_t) (idect - fDect1[1]));
    CalculT0CoefMean(fCountDet[1],(Int_t) (idect - fDect1[1]));
    
    // Method Mean and fit
    // idect is egal for fDebug = 0 and 2, only to fill the hist 
    if(fFitPol2On){
      FitPente((TH1 *) projph,(Int_t) (idect - fDect1[1]));
    }
    if(fFitLagrPolOn){
      FitLagrangePoly((TH1 *) projph,(Int_t) (idect - fDect1[1]));
    }
    if (fFitPHOn) { 
      FitPH((TH1 *) projph,(Int_t) (idect - fDect1[1]));
    }
  
    // Visualise the detector for fDebug 3 or 4
    // Here is the reconstruction of the pad and row group is used!
    if (fDebug >= 3) {
      FillCoefVdriftDB();
      FillCoefT0DB();
    }

    // Fill the tree if end of a detector or only the pointer to the branch!!!
    FillInfosFit(idect,1);

    // Memory!!!
    if (fDebug != 2) {
      delete projph;
    }
       
  } // Boucle object
 
  // Plot
  // 0 no plot, 1 and 4 error plot, 3 and 4 DB plot
  if ((fDebug == 1) || 
      (fDebug == 4)){
    PlotWritePH();
    PlotWriteT0();
  }
  if ((fDebug == 4) || 
      (fDebug == 3)){
    PlotPHDB();
    PlotT0DB();
  }

  // Mean Statistics
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
  
  // Write the things!
  if (fWriteCoef[1]) {
    WriteFitInfos(1);
  }

  return kTRUE;
  
}

//____________Functions fit Online PRF2d_______________________________________
Bool_t AliTRDCalibraFit::FitPRFOnline(TProfile2D *prf)
{
  //
  // Take the 1D profiles (pad response function), projections of the 2D PRF
  // on the Xaxis, for each calibration group
  // Fit with a gaussian to reconstruct the sigma of the pad response function
  // write the results in a tree
  //

  // A small check
  if ((fFitPRFNDB == 2) && (!fRMSPRFOn)){
    AliInfo("You have choosen to write the RMS method but it is not on!");
    return kFALSE;
  }
  if ((fFitPRFNDB == 0) && (!fFitPRFOn)){
    AliInfo("You have choosen to write the fit method but it is not on!");
    return kFALSE;
  }

  // Number of Xbins (detectors or groups of pads)
  TAxis   *xprf    = prf->GetXaxis();
  TAxis   *yprf    = prf->GetYaxis();
  Int_t    nybins  = yprf->GetNbins();
  Int_t    nbins   = xprf->GetNbins();
  if (!InitFit(nbins,2)) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;

  // Init fCountDet and fCount
  InitfCountDetAndfCount(2);
  
  // Beginning of the loop
  for (Int_t idect = fDect1[2]; idect < fDect2[2]; idect++) {

    TH1D *projprf = (TH1D *) prf->ProjectionY("projprf",idect+1,idect+1,(Option_t *) "e");
    projprf->SetDirectory(0);
    
    // Number of entries for this calibration group
    Double_t nentries = 0;
    for (Int_t k = 0; k < nybins; k++) {
      nentries += prf->GetBinEntries(prf->GetBin(idect+1,k+1));
    }
    if(nentries > 0) fNumberEnt++;
   
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    UpdatefCountDetAndfCount(idect,2);
    
    // Reconstruction of the row and pad group: rowmin, row max ...
    ReconstructFitRowMinRowMax(idect,2);
    
    // Rebin and statistic stuff
    // This detector has not enough statistics or was off
    if (nentries < fMinEntries) {
      
      // Fill with the default values
      NotEnoughStatistic(idect,2);
      
      // Memory!
      if (fDebug != 2) {
	delete projprf;
      }
      
      continue;

    }
    
    // Statistics of the histos fitted
    AliInfo(Form("For the group number %d there are %f stats",idect,nentries));
    fNumberFit++;
    fStatisticMean += nentries;
    
    // Calcul of "real" coef
    if ((fDebug == 1) || 
        (fDebug == 4)) {
      CalculPRFCoefMean(fCountDet[2],(Int_t) (idect - fDect1[2]));
    }
    
    // Method Mean and fit
    // idect is egal for fDebug = 0 and 2, only to fill the hist 
    if(fFitPRFOn){
      FitPRF((TH1 *) projprf,(Int_t) (idect - fDect1[2]));
    }
    if(fRMSPRFOn){
      RmsPRF((TH1 *) projprf,(Int_t) (idect - fDect1[2]));
    }
   
    
    // Visualise the detector for fDebug 3 or 4
    // Here is the reconstruction of the pad and row group is used!
    if (fDebug >= 3) {
      FillCoefPRFDB();
    }

    // Fill the tree if end of a detector or only the pointer to the branch!!!
    FillInfosFit(idect,2);
    
    // Memory!!!
    if (fDebug != 2) {
      delete projprf;
    }
    
  } // Boucle object
 
  // Plot
  // No plot, 1 and 4 error plot, 3 and 4 DB plot
  if ((fDebug == 1) || 
      (fDebug == 4)) {
    PlotWritePRF();
  }
  if ((fDebug == 4) || 
      (fDebug == 3)){
    PlotPRFDB();
  }

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

  // Write the things!
  if (fWriteCoef[2]) {
    WriteFitInfos(2);
  }

  return kTRUE;
  
}

//____________Functions fit Online PRF2d_______________________________________
Bool_t AliTRDCalibraFit::FitPRFOnline()
{
  //
  // Reconstruct the 1D histo (pad response function) from the vectorPRD for
  // each calibration group
  // Fit with a gaussian to reconstruct the sigma of the pad response function
  // write the results in a tree
  //

  // A small check
  if ((fFitPRFNDB == 2) && (!fRMSPRFOn)){
    AliInfo("You have choosen to write the RMS method but it is not on!");
    return kFALSE;
  }
  if ((fFitPRFNDB == 0) && (!fFitPRFOn)){
    AliInfo("You have choosen to write the fit method but it is not on!");
    return kFALSE;
  }

  // Warning
  if (!fCalibraVector) {
    AliError("You have first to set the calibravector before using this function!");
    return kFALSE;
  }

  // Number of Xbins (detectors or groups of pads)
  if (!InitFit(0,2)) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;

  // Init fCountDet and fCount
  InitfCountDetAndfCount(2);

  // Beginning of the loop
  for (Int_t idect = fDect1[2]; idect < fDect2[2]; idect++) {

    // Search if the group is in the VectorCH
    Int_t place = fCalibraVector->SearchInVector(idect,2);
    
    // Is in
    TH1F   *projprf = 0x0;
    TString name("PRF");
    name += idect;
    if (place != -1) {
      //Entries
      fNumberEnt++;
      projprf = CorrectTheError((TGraphErrors *) (fCalibraVector->ConvertVectorPHisto(place,(const char *)name)));
      projprf->SetDirectory(0);
    }

    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    UpdatefCountDetAndfCount(idect,2);

    // Reconstruction of the row and pad group: rowmin, row max ...
    ReconstructFitRowMinRowMax(idect,2);

    // Rebin and statistic stuff
    // This detector has not enough statistics or was off
    if ((place == -1) ||
        ((place           != -1) &&
         (fEntriesCurrent < fMinEntries))) {

      // Fill with the default values
      NotEnoughStatistic(idect,2);

      // Memory
      if (fDebug != 2) {
	delete projprf;
      }

      continue;

    }

    // Statistic of the histos fitted
    AliInfo(Form("For the group number %d there are %d stats", idect,fEntriesCurrent));
    fNumberFit++;
    fStatisticMean += fEntriesCurrent;

    // Calcul of "real" coef
    if ((fDebug == 1) || 
        (fDebug == 4)) {
      CalculPRFCoefMean(fCountDet[2],(Int_t) (idect-fDect1[2]));
    }

    // Method Mean and fit
    // idect is egal for fDebug = 0 and 2, only to fill the hist 
    if(fFitPRFOn){
      FitPRF((TH1 *) projprf,(Int_t) (idect-fDect1[2]));
    }
    if(fRMSPRFOn){
      RmsPRF((TH1 *) projprf,(Int_t) (idect-fDect1[2]));
    }
    
    // Visualise the detector for fDebug 3 or 4
    // Here is the reconstruction of the pad and row group is used!
    if (fDebug >= 3) {
      FillCoefPRFDB();
    }
    // Fill the tree if end of a detector or only the pointer to the branch!!!
    FillInfosFit(idect,2);
                    
    // Memory!!!
    if (fDebug != 2) {
      delete projprf;
    }
    
  } // Boucle object

  // Plot
  // No plot, 1 and 4 error plot, 3 and 4 DB plot
  if ((fDebug == 1) || 
      (fDebug == 4)) {
    PlotWritePRF();
  }
  if ((fDebug == 4) || 
      (fDebug == 3)) {
    PlotPRFDB();
  }

  // Mean Statistics
  if (fNumberFit > 0) {
    AliInfo(Form("There are %d with at least one entries.",fNumberEnt));
    AliInfo(Form("%d fits have been proceeded (sucessfully or not...).",fNumberFit));
    AliInfo(Form("There is a mean statistic of: %d over these fitted histograms and %d successfulled fits"
                ,(Int_t) fStatisticMean/fNumberFit,fNumberFitSuccess));
  }
  else {
    AliInfo(Form("There are %d with at least one entries. There is no fit!",fNumberEnt));
  }

  // Write the things!
  if (fWriteCoef[2]) {
    WriteFitInfos(2);
  }

  return kTRUE;
  
}

//____________Functions fit Online PRF2d_______________________________________
Bool_t AliTRDCalibraFit::FitPRFOnline(TTree *tree)
{
  //
  // Look if the calibration group can be found in the tree, if yes take
  // the histo, fit it, and write the results in a tree
  //

  // A small check
  if ((fFitPRFNDB == 2) && (!fRMSPRFOn)){
    AliInfo("You have choosen to write the RMS method but it is not on!");
    return kFALSE;
  }
  if ((fFitPRFNDB == 0) && (!fFitPRFOn)){
    AliInfo("You have choosen to write the fit method but it is not on!");
    return kFALSE;
  }

  // Number of Xbins (detectors or groups of pads)
  if (!InitFit(0,2)) {
    return kFALSE;
  }
  fStatisticMean        = 0.0;
  fNumberFit            = 0;
  fNumberFitSuccess     = 0;
  fNumberEnt            = 0;

  // Initialise
  fCalibraVector = new AliTRDCalibraVector();

  // Init fCountDet and fCount
  InitfCountDetAndfCount(2);
  TGraphErrors *projprftree = 0x0;
  tree->SetBranchAddress("histo",&projprftree);
  TObjArray    *vectorplace = fCalibraVector->ConvertTreeVector(tree);

  // Beginning of the loop
  for (Int_t idect = fDect1[2]; idect < fDect2[2]; idect++) {

    // Search if the group is in the VectorCH
    Int_t place = fCalibraVector->SearchInTreeVector(vectorplace,idect);
    
    // Is in   
    TH1F *projprf = 0x0;
    if (place != -1) {
      //Entries
      fNumberEnt++;
      // Variable
      tree->GetEntry(place);
      projprf = CorrectTheError(projprftree);
    }

    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    UpdatefCountDetAndfCount(idect,2);

    // Reconstruction of the row and pad group: rowmin, row max ...
    ReconstructFitRowMinRowMax(idect,2);

    // Rebin and statistic stuff
    // This detector has not enough statistics or was off
    if ((place == -1) ||
        ((place           !=          -1) &&
         (fEntriesCurrent  < fMinEntries))) {
      
      // Fill with the default values
      NotEnoughStatistic(idect,2);
      
      // Memory!!!
      if (fDebug != 2) {
	delete projprf;
      }
      
      continue;

    }

    // Statistics of the histos fitted
    AliInfo(Form("For the group number %d there are %d stats",idect,fEntriesCurrent));
    fNumberFit++;
    fStatisticMean += fEntriesCurrent;
        
    // Calcul of "real" coef
    if ((fDebug == 1) || 
        (fDebug == 4)){
      CalculPRFCoefMean(fCountDet[2],(Int_t) (idect - fDect1[2]));
    }
    
    // Method Mean and fit
    // idect is egal for fDebug = 0 and 2, only to fill the hist 
    if(fFitPRFOn){
      FitPRF((TH1 *) projprf,(Int_t) (idect - fDect1[2]));
    }
    if(fRMSPRFOn){
      RmsPRF((TH1 *) projprf,(Int_t) (idect - fDect1[2]));
    }
   
    // Visualise the detector for fDebug 3 or 4
    // Here is the reconstruction of the pad and row group is used!
    if (fDebug >= 3) {
      FillCoefPRFDB();
    }
    // Fill the tree if end of a detector or only the pointer to the branch!!!
    FillInfosFit(idect,2);

    // Memory!!!
    if (fDebug != 2) {
      delete projprf;
    }                
    
  } // Boucle object
  
  // Plot
  // No plot, 1 and 4 error plot, 3 and 4 DB plot
  if ((fDebug == 1) || 
      (fDebug == 4)){
    PlotWritePRF();
  }
  if ((fDebug == 4) || 
      (fDebug == 3)){
    PlotPRFDB();
  }
  
  // Mean Statistics
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

  // Write the things!
  if (fWriteCoef[2]) {
    WriteFitInfos(2);
  }

  return kTRUE;
  
}

//____________Functions for seeing if the pad is really okey___________________

//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::SetModeCalibrationFromTObject(TObject *object, Int_t i)
{
  //
  // Set fNz[i] and fNrphi[i] of the AliTRDCalibraFit::Instance()
  // corresponding to the given TObject
  //

  const char *nametitle = object->GetTitle();

  // Some patterns
  const Char_t *patternz0    = "Nz0";
  const Char_t *patternz1    = "Nz1";
  const Char_t *patternz2    = "Nz2";
  const Char_t *patternz3    = "Nz3";
  const Char_t *patternz4    = "Nz4";
  const Char_t *patternrphi0 = "Nrphi0";
  const Char_t *patternrphi1 = "Nrphi1";
  const Char_t *patternrphi2 = "Nrphi2";
  const Char_t *patternrphi3 = "Nrphi3";
  const Char_t *patternrphi4 = "Nrphi4";
  const Char_t *patternrphi5 = "Nrphi5";
  const Char_t *patternrphi6 = "Nrphi6";

  UShort_t testz    = 0;
  UShort_t testrphi = 0;

  // Nz mode
  if (strstr(nametitle,patternz0)) {
    testz++;
    fCalibraMode->SetNz(i, 0);
  }
  if (strstr(nametitle,patternz1)) {
    testz++;
    fCalibraMode->SetNz(i ,1);
  }
  if (strstr(nametitle,patternz2)) {
    testz++;
    fCalibraMode->SetNz(i ,2);
  }
  if (strstr(nametitle,patternz3)) {
    testz++;
    fCalibraMode->SetNz(i ,3);
  }
  if (strstr(nametitle,patternz4)) {
    testz++;
    fCalibraMode->SetNz(i ,4);
  }

  // Nrphi mode
  if (strstr(nametitle,patternrphi0)) {
    testrphi++;
    fCalibraMode->SetNrphi(i ,0);
  }
  if (strstr(nametitle,patternrphi1)) {
    testrphi++;
    fCalibraMode->SetNrphi(i, 1);
  }
  if (strstr(nametitle,patternrphi2)) {
    testrphi++;
    fCalibraMode->SetNrphi(i, 2);
  }
  if (strstr(nametitle,patternrphi3)) {
    testrphi++;
    fCalibraMode->SetNrphi(i, 3);
  }
  if (strstr(nametitle,patternrphi4)) {
    testrphi++;
    fCalibraMode->SetNrphi(i, 4);
  }
  if (strstr(nametitle,patternrphi5)) {
    testrphi++;
    fCalibraMode->SetNrphi(i, 5);
  }
  if (strstr(nametitle,patternrphi6)) {
    testrphi++;
    fCalibraMode->SetNrphi(i, 6);
  }
 
  // Look if all is okey
  if ((testz    == 1) && 
      (testrphi == 1)) {
    return kTRUE;
  }
  else {
    fCalibraMode->SetNrphi(i ,0);
    fCalibraMode->SetNz(i ,0);
    return kFALSE;
  }
  
}

//_____________________________________________________________________________
AliTRDCalDet *AliTRDCalibraFit::CreateDetObjectTree(TTree *tree, Int_t i)
{
  //
  // It creates the AliTRDCalDet object from the tree of the coefficient
  // for the calibration i (i != 2)
  // It takes the mean value of the coefficients per detector 
  // This object has to be written in the database
  //
  
  // Create the DetObject
  AliTRDCalDet *object = 0x0;
  if (i == 0) {
    object = new AliTRDCalDet("ChamberGainFactor","GainFactor (detector value)");
  }
  if (i == 1) {
    object = new AliTRDCalDet("ChamberVdrift","TRD drift velocities (detector value)");
  }
  else {
    object = new AliTRDCalDet("ChamberT0","T0 (detector value)");
  }
  
  // Read the Tree
  Int_t   detector = -1;
  Float_t values[2304];
  tree->SetBranchAddress("detector",&detector);
  if (i == 0) {
    tree->SetBranchAddress("gainPad",values);
  }
  if (i == 1) {
    tree->SetBranchAddress("vdrift" ,values);
  }
  if (i == 3) {
    tree->SetBranchAddress("t0"     ,values);
  }
  
  // For calculating the mean
  Float_t mean            = 0.0;
  Int_t   nto             = 0;
  Int_t   numberofentries = tree->GetEntries();
  
  if (numberofentries != 540) {
    AliInfo("The tree is not complete");
  }
  
  for (Int_t det = 0; det < numberofentries; ++det) {
    tree->GetEntry(det);
    if (GetChamber(detector) == 2) {
      nto = 1728;
    }
    else {
      nto = 2304;
    }
    mean = 0.0;
    if(i != 3){
      for (Int_t k = 0; k < nto; k++) {
	mean += TMath::Abs(values[k]) / nto;  
      }
    }
    else {
      for (Int_t k = 0; k < nto; k++) {
	if(k == 0) mean = values[k];
	if(mean > values[k]) mean = values[k];
      }
    }
    object->SetValue(detector,mean);
  }

  return object;

}

//_____________________________________________________________________________
TObject *AliTRDCalibraFit::CreatePadObjectTree(TTree *tree, Int_t i
                                          , AliTRDCalDet *detobject)
{
  //
  // It Creates the AliTRDCalPad object from the tree of the
  // coefficient for the calibration i (i != 2)
  // You need first to create the object for the detectors,
  // where the mean value is put.
  // This object has to be written in the database
  //
  
  // Create the DetObject
  AliTRDCalPad *object = 0x0;
  if (i == 0) {
    object = new AliTRDCalPad("GainFactor","GainFactor (local variations)");
  }
  if (i == 1) {
    object = new AliTRDCalPad("LocalVdrift","TRD drift velocities (local variations)");
  }
  else {
    object = new AliTRDCalPad("LocalT0","T0 (local variations)");
  }
  
  // Read the Tree
  Int_t   detector = -1;
  Float_t values[2304];
  tree->SetBranchAddress("detector",&detector);
  if (i == 0) {
    tree->SetBranchAddress("gainPad",values);
  }
  if (i == 1) {
    tree->SetBranchAddress("vdrift" ,values);
  }
  if (i == 3) {
    tree->SetBranchAddress("t0"     ,values);
  }
  
  // Variables
  Float_t mean            = 0.0;
  Int_t   numberofentries = tree->GetEntries();
  
  if (numberofentries != 540) {
    AliInfo("The tree is not complete");
  }
  
  for (Int_t det = 0; det < numberofentries; ++det) {
    tree->GetEntry(det);
    AliTRDCalROC *calROC = object->GetCalROC(detector);
    mean = detobject->GetValue(detector);
    if ((mean == 0) && (i != 3)) {
      continue;
    }
    Int_t rowMax = calROC->GetNrows();
    Int_t colMax = calROC->GetNcols();
    for (Int_t row = 0; row < rowMax; ++row) {
      for (Int_t col = 0; col < colMax; ++col) {
	if(i != 3) calROC->SetValue(col,row,TMath::Abs(values[(Int_t) (col*rowMax+row)])/mean);
	else calROC->SetValue(col,row,values[(Int_t) (col*rowMax+row)]-mean);
	
      } // Col
    } // Row
  }

  return object;

}

//_____________________________________________________________________________
TObject *AliTRDCalibraFit::CreatePadObjectTree(TTree *tree)
{
  //
  // It Creates the AliTRDCalPad object from the tree of the
  // coefficient for the calibration PRF (i = 2)
  // This object has to be written in the database
  //
  
  // Create the DetObject
  AliTRDCalPad *object = new AliTRDCalPad("PRFWidth","PRFWidth");

  // Read the Tree
  Int_t   detector = -1;
  Float_t values[2304];
  tree->SetBranchAddress("detector",&detector);
  tree->SetBranchAddress("width"   ,values);
   
  // Variables
  Int_t numberofentries = tree->GetEntries();

  if (numberofentries != 540) {
    AliInfo("The tree is not complete");
  }

  for (Int_t det = 0; det < numberofentries; ++det) {
    tree->GetEntry(det);
    AliTRDCalROC *calROC = object->GetCalROC(detector);
    Int_t rowMax = calROC->GetNrows();
    Int_t colMax = calROC->GetNcols();
    for (Int_t row = 0; row < rowMax; ++row) {
      for (Int_t col = 0; col < colMax; ++col) {
	calROC->SetValue(col,row,TMath::Abs(values[(Int_t) (col*rowMax+row)]));
      } // Col
    } // Row
  }

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
void AliTRDCalibraFit::SetFitPRFNDB(Int_t fitPRFNDB)
{ 
  //
  // TO choose the method that you write into the database
  //

  if ((fitPRFNDB >= 3) || (fitPRFNDB == 1)) {
    AliInfo("fitPRFNDB is not a correct number!"); 
  }
  else {
    fFitPRFNDB = fitPRFNDB;
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::SetFitChargeNDB(Int_t fitChargeNDB)
{ 
  //
  // To choose the method that you write into the database
  //
  if ((fitChargeNDB >= 5) || (fitChargeNDB == 3)) {
    AliInfo("fitChargeNDB is not a correct number!"); 
  }
  else {
    fFitChargeNDB = fitChargeNDB;
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::SetFitPHNDB(Int_t fitPHNDB)
{ 
  //
  // To choose the method that you write into the database
  //

  if ((fitPHNDB >= 4) || (fitPHNDB == 2)) {
    AliInfo("fitPHNDB is not a correct number!"); 
  }
  else {
    fFitPHNDB = fitPHNDB;
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
void AliTRDCalibraFit::SetT0Shift(Float_t t0Shift) 
{ 
  //
  // The t0 calculated with the maximum positif slope is shift from t0Shift
  // You can here set t0Shift
  //

  if (t0Shift > 0) {
    fT0Shift = t0Shift; 
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

//____________Pad Calibration Public___________________________________________

//____________Protected Functions______________________________________________
//____________Create the 2D histo to be filled online__________________________
//
//____________Fit______________________________________________________________
//____________Create histos if fDebug == 1 or fDebug >= 3______________________

//_____________________________________________________________________________
void AliTRDCalibraFit::InitArrayFitPH()
{
  //
  // Initialise fCoefVdrift[3] and fCoefVdriftE[2] to the right dimension
  //
  
  Int_t nbins = fDect2[1]-fDect1[1];

  fCoefVdrift[2] = new Double_t[nbins];

  // Init the pointer to nbins
  if (fFitPHOn) {
    fCoefVdrift[0] = new Double_t[nbins];
    fCoefVdriftE[0] = new Double_t[nbins];
    for(Int_t k = 0; k < nbins; k++){
      fCoefVdriftE[0][k] = 0.0;
    }
  }


  if (fFitPol2On){
    fCoefVdrift[1] = new Double_t[nbins];
    fCoefVdriftE[1] = new Double_t[nbins];
    for(Int_t k = 0; k < nbins; k++){
      fCoefVdriftE[1][k] = 0.0;
    }
  }
  if (fFitLagrPolOn){
    fCoefVdrift[3] = new Double_t[nbins];
    fCoefVdriftE[2] = new Double_t[nbins];
    for(Int_t k = 0; k < nbins; k++){
      fCoefVdriftE[2][k] = 0.0;
    }
  }
 
}

//_____________________________________________________________________________
void AliTRDCalibraFit::InitArrayFitT0()
{
  //
  // Initialise fCoefT0[3] and fCoefT0E[2] to the right dimension
  //
  
  Int_t nbins = fDect2[1]-fDect1[1];

  fCoefT0[2] = new Double_t[nbins];

  // Init the pointer to nbins
  if(fFitPHOn){
    fCoefT0[0] = new Double_t[nbins];
    fCoefT0E[0] = new Double_t[nbins];
    for(Int_t k = 0; k < nbins; k++){
      fCoefT0E[0][k] = 0.0;
    }
  }
  if(fFitPol2On){
    fCoefT0[1] = new Double_t[nbins];
    fCoefT0E[1] = new Double_t[nbins];
    for(Int_t k = 0; k < nbins; k++){
      fCoefT0E[1][k] = 0.0;
    }
  }
  if(fFitLagrPolOn){
    fCoefT0[3] = new Double_t[nbins];
    fCoefT0E[2] = new Double_t[nbins];
    for(Int_t k = 0; k < nbins; k++){
      fCoefT0E[2][k] = 0.0;
    }
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::InitArrayFitCH()
{
  //
  // Initialise fCoefCharge[4] and fCoefChargeE[3] to the right dimension
  //

  Int_t nbins = fDect2[0]-fDect1[0];

  //Init the pointer to nbins
  if(fMeanChargeOn){
    fCoefCharge[1] = new Double_t[nbins];
    fCoefChargeE[1] = new Double_t[nbins];
    for(Int_t k = 0; k < nbins; k++){
      fCoefChargeE[1][k] = 0.0;
    }
  }
  if(fFitMeanWOn){
    fCoefCharge[4] = new Double_t[nbins];
    fCoefChargeE[3] = new Double_t[nbins];
    for(Int_t k = 0; k < nbins; k++){
      fCoefChargeE[3][k] = 0.0;
    }
  }
  if(fFitChargeOn){
    fCoefCharge[0] = new Double_t[nbins];
    fCoefChargeE[0] = new Double_t[nbins];
    for(Int_t k = 0; k < nbins; k++){
      fCoefChargeE[0][k] = 0.0;
    }
  }

  if(fFitChargeBisOn){
    fCoefCharge[2] = new Double_t[nbins];
    fCoefChargeE[2] = new Double_t[nbins];
    for(Int_t k = 0; k < nbins; k++){
      fCoefChargeE[2][k] = 0.0;
    }
  }

  fCoefCharge[3] = new Double_t[nbins];
  
}

//_____________________________________________________________________________
void AliTRDCalibraFit::InitArrayFitPRF()
{
  //
  // Initialise fCoefPRF[2] and fCoefPRFE to the right dimension
  //

  Int_t nbins = fDect2[2]-fDect1[2];
  fCoefPRF[1] = new Double_t[nbins];

  //Init the pointer to nbins
  if(fFitPRFOn){
    fCoefPRF[0] = new Double_t[nbins];
    fCoefPRFE[0] = new Double_t[nbins];
    for(Int_t k = 0; k < nbins; k++){
      fCoefPRFE[0][k] = 0.0;
    }
  }
  if(fRMSPRFOn){
    fCoefPRF[2] = new Double_t[nbins];
    fCoefPRFE[1] = new Double_t[nbins];
    for(Int_t k = 0; k < nbins; k++){
      fCoefPRFE[1][k] = 0.0;
    }
  }
}

//_____________________________________________________________________________
void AliTRDCalibraFit::CreateFitHistoPRFDB(Int_t rowMax, Int_t colMax)
{
  //
  // Create the histos for fDebug = 3 and fDebug = 4 (Fit functions)
  //
  if(fFitPRFOn){
    fCoefPRFDB[0] = new TH2F("coefPRF0","",rowMax,0,rowMax,colMax,0,colMax);
    fCoefPRFDB[0]->SetStats(0);
    fCoefPRFDB[0]->SetXTitle("row Number");
    fCoefPRFDB[0]->SetYTitle("col Number");
    fCoefPRFDB[0]->SetZTitle("PRF width [pad width units]");
    fCoefPRFDB[0]->SetFillColor(6);
    fCoefPRFDB[0]->SetLineColor(6);
  }
  if(fRMSPRFOn){
    fCoefPRFDB[1] = new TH2F("coefPRF1","",rowMax,0,rowMax,colMax,0,colMax);
    fCoefPRFDB[1]->SetStats(0);
    fCoefPRFDB[1]->SetXTitle("row Number");
    fCoefPRFDB[1]->SetYTitle("col Number");
    fCoefPRFDB[1]->SetZTitle("PRF width [pad width units]");
    fCoefPRFDB[1]->SetFillColor(1);
    fCoefPRFDB[1]->SetLineColor(1);
  }
}

//_____________________________________________________________________________
void AliTRDCalibraFit::CreateFitHistoCHDB(Int_t rowMax, Int_t colMax)
{
  //
  // Create the histos for fDebug = 3 and fDebug = 4 (Fit functions)
  //

  if(fFitChargeOn){
  fCoefChargeDB[0] = new TH2F("coefchargedb0","",rowMax,0,rowMax,colMax,0,colMax);
  fCoefChargeDB[0]->SetStats(0);
  fCoefChargeDB[0]->SetXTitle("row Number");
  fCoefChargeDB[0]->SetYTitle("col Number");
  fCoefChargeDB[0]->SetZTitle("f_{g} Fit method");
  fCoefChargeDB[0]->SetFillColor(6);
  fCoefChargeDB[0]->SetLineColor(6);
  }
  if(fFitChargeBisOn){
    fCoefChargeDB[2] = new TH2F("coefchargedb2","",rowMax,0,rowMax,colMax,0,colMax);
    fCoefChargeDB[2]->SetStats(0);
    fCoefChargeDB[2]->SetXTitle("row Number");
    fCoefChargeDB[2]->SetYTitle("col Number");
    fCoefChargeDB[2]->SetZTitle("f_{g} Fitbis method");
    fCoefChargeDB[2]->SetFillColor(8);
    fCoefChargeDB[2]->SetLineColor(8);
  }
  if(fFitMeanWOn){
    fCoefChargeDB[3] = new TH2F("coefchargedb3","",rowMax,0,rowMax,colMax,0,colMax);
    fCoefChargeDB[3]->SetStats(0);
    fCoefChargeDB[3]->SetXTitle("row Number");
    fCoefChargeDB[3]->SetYTitle("col Number");
    fCoefChargeDB[3]->SetFillColor(1);
    fCoefChargeDB[3]->SetLineColor(1);
 
  }
  if(fMeanChargeOn){
  fCoefChargeDB[1] = new TH2F("coefchargedb1","",rowMax,0,rowMax,colMax,0,colMax);
  fCoefChargeDB[1]->SetStats(0);
  fCoefChargeDB[1]->SetXTitle("row Number");
  fCoefChargeDB[1]->SetYTitle("col Number");
  fCoefChargeDB[1]->SetZTitle("f_{g} Mean method");
  fCoefChargeDB[1]->SetFillColor(2);
  fCoefChargeDB[1]->SetLineColor(2);
  } 

}

//_____________________________________________________________________________
void AliTRDCalibraFit::CreateFitHistoPHDB(Int_t rowMax, Int_t colMax)
{
  //
  // Create the histos for fDebug = 3 and fDebug = 4 (Fit functions)
  //

  if(fFitPHOn){
  fCoefVdriftDB[0] = new TH2F("coefvdriftdb0","",rowMax,0,rowMax,colMax,0,colMax);
  fCoefVdriftDB[0]->SetStats(0);
  fCoefVdriftDB[0]->SetXTitle("row Number");
  fCoefVdriftDB[0]->SetYTitle("col Number");
  fCoefVdriftDB[0]->SetZTitle("v_{drift} Fit method");
  fCoefVdriftDB[0]->SetFillColor(6);
  fCoefVdriftDB[0]->SetLineColor(6);
  }

  if(fFitPol2On){
    fCoefVdriftDB[1] = new TH2F("coefvdriftdb1","",rowMax,0,rowMax,colMax,0,colMax);
    fCoefVdriftDB[1]->SetStats(0);
    fCoefVdriftDB[1]->SetXTitle("row Number");
    fCoefVdriftDB[1]->SetYTitle("col Number");
    fCoefVdriftDB[1]->SetZTitle("v_{drift} slope method");
    fCoefVdriftDB[1]->SetFillColor(2);
    fCoefVdriftDB[1]->SetLineColor(2);
  }
  if(fFitLagrPolOn){
    fCoefVdriftDB[2] = new TH2F("coefvdriftdb1","",rowMax,0,rowMax,colMax,0,colMax);
    fCoefVdriftDB[2]->SetStats(0);
    fCoefVdriftDB[2]->SetXTitle("row Number");
    fCoefVdriftDB[2]->SetYTitle("col Number");
    fCoefVdriftDB[2]->SetZTitle("v_{drift} slope method");
    fCoefVdriftDB[2]->SetFillColor(1);
    fCoefVdriftDB[2]->SetLineColor(1);
  }
 
}

//_____________________________________________________________________________
void AliTRDCalibraFit::CreateFitHistoT0DB(Int_t rowMax, Int_t colMax)
{
  //
  // Create the histos for fDebug = 3 and fDebug = 4 (Fit functions)
  //

  if(fFitPHOn){
    fCoefT0DB[0] = new TH2F("coefT0db0","",rowMax,0,rowMax,colMax,0,colMax);
    fCoefT0DB[0]->SetStats(0);
    fCoefT0DB[0]->SetXTitle("row Number");
    fCoefT0DB[0]->SetYTitle("col Number");
    fCoefT0DB[0]->SetZTitle("t0 Fit method");
    fCoefT0DB[0]->SetFillColor(6);
    fCoefT0DB[0]->SetLineColor(6);
  }
  if(fFitPol2On){
    fCoefT0DB[1] = new TH2F("coefT0db1","",rowMax,0,rowMax,colMax,0,colMax);
    fCoefT0DB[1]->SetStats(0);
    fCoefT0DB[1]->SetXTitle("row Number");
    fCoefT0DB[1]->SetYTitle("col Number");
    fCoefT0DB[1]->SetZTitle("t0 slope method");
    fCoefT0DB[1]->SetFillColor(2);
    fCoefT0DB[1]->SetLineColor(2);
  }
  if(fFitLagrPolOn){
    fCoefT0DB[2] = new TH2F("coefT0db1","",rowMax,0,rowMax,colMax,0,colMax);
    fCoefT0DB[2]->SetStats(0);
    fCoefT0DB[2]->SetXTitle("row Number");
    fCoefT0DB[2]->SetYTitle("col Number");
    fCoefT0DB[2]->SetZTitle("t0 slope method");
    fCoefT0DB[2]->SetFillColor(1);
    fCoefT0DB[2]->SetLineColor(1);
  }
  
}

//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::FillVectorFitCH(Int_t countdet)
{
  //
  // For the Fit functions fill the vector FitCH special for the gain calibration
  //

  AliTRDFitCHInfo *fitCHInfo = new AliTRDFitCHInfo();

  Int_t ntotal = 1;
  if (GetChamber(countdet) == 2) {
    ntotal = 1728;
  }
  else {
    ntotal = 2304;
  }

  //printf("For the detector %d , ntotal %d and fCoefCH[0] %f\n",countdet,ntotal,fCoefCH[0]);
  Float_t *coef = new Float_t[ntotal];
  for (Int_t i = 0; i < ntotal; i++) {
    coef[i] = fCoefCH[i];
  }
  
  Int_t detector = countdet;
  // Set
  fitCHInfo->SetCoef(coef);
  fitCHInfo->SetDetector(detector);
  fVectorFitCH->Add((TObject *) fitCHInfo);

  return kTRUE;

}

//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::InitFit(Int_t nbins, Int_t i)
{
  //
  // Init the calibration mode (Nz, Nrphi), the histograms for
  // debugging the fit methods if fDebug > 0, 
  //

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.01);

  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam");
    return kFALSE;
  }

  // Mode groups of pads: the total number of bins!
  Int_t numberofbinsexpected = 0;
  fCalibraMode->ModePadCalibration(2,i);
  fCalibraMode->ModePadFragmentation(0,2,0,i);
  fCalibraMode->SetDetChamb2(i);
  if (fDebug == 1) {
    AliInfo(Form("For the chamber 2: %d",fCalibraMode->GetDetChamb2(i)));
  }
  numberofbinsexpected += 6 * 18 * fCalibraMode->GetDetChamb2(i);
  fCalibraMode->ModePadCalibration(0,i);
  fCalibraMode->ModePadFragmentation(0,0,0,i);
  fCalibraMode->SetDetChamb0(i);
  if (fDebug == 1) {
    AliInfo(Form("For the other chamber 0: %d",fCalibraMode->GetDetChamb0(i)));
  }
  numberofbinsexpected += 6 * 4 * 18 * fCalibraMode->GetDetChamb0(i);
  
  // Quick verification that we have the good pad calibration mode if 2D histos!
  if (nbins != 0) {
    if (numberofbinsexpected != nbins) {
      AliInfo("It doesn't correspond to the mode of pad group calibration!");
      return kFALSE;
    }
  }

  // Security for fDebug 3 and 4
  if ((fDebug >= 3) && 
      ((fDet[0] >  5) || 
       (fDet[1] >  4) || 
       (fDet[2] > 17))) {
    AliInfo("This detector doesn't exit!");
    return kFALSE;
  }

  // Determine fDet1 and fDet2
  fDect1[i] = -1;
  fDect2[i] = -1;
  if (fDebug == 2) {
    fDect1[i] = fFitVoir;
    fDect2[i] = fDect1[i] +1;
  }
  if (fDebug <= 1) {
    fDect1[i] = 0;
    fDect2[i] = numberofbinsexpected;
  }
  if (fDebug >= 3) {
    fCalibraMode->CalculXBins(AliTRDgeometry::GetDetector(fDet[0],fDet[1],fDet[2]),i);
    fDect1[i] = fCalibraMode->GetXbins(i);
    fCalibraMode->CalculXBins((AliTRDgeometry::GetDetector(fDet[0],fDet[1],fDet[2])+1),i);
    fDect2[i] = fCalibraMode->GetXbins(i);
  }

  // Create the histos for debugging
  // CH
  if (i == 0) {
    
    gDirectory = gROOT;
    // Init the VectorFitCH
    fVectorFitCH = new TObjArray();
    fCoefCH      = new Float_t[2304];
    for (Int_t k = 0; k < 2304; k++) {
      fCoefCH[k] = 0.0;    
    }
    fScaleFitFactor = 0.0;

    // Number of Xbins(detectors or groups of pads) if Vector2d
    // Quick verification that we are not out of range!
    if (fCalibraVector) {
      if ((nbins                            == 0) && 
          (fCalibraVector->GetVectorCH()->GetEntriesFast()      >  0) && 
          ((Int_t) fCalibraVector->GetPlaCH()->GetEntriesFast() >  0)) {
	if ((Int_t) fCalibraVector->GetVectorCH()->GetEntriesFast() > numberofbinsexpected) {
	  AliInfo("ch doesn't correspond to the mode of pad group calibration!");
	  return kFALSE;
	}
	if ((Int_t) fCalibraVector->GetVectorCH()->GetEntriesFast() != 
            (Int_t) fCalibraVector->GetPlaCH()->GetEntriesFast()) {
	  AliInfo("VectorCH doesn't correspond to PlaCH!");
	  return kFALSE;
	}
      }
    }

    //
    // Debugging: Create the histos
    //

    // fDebug == 0 nothing
    
    // fDebug == 1 
    if (fDebug == 1) {
      InitArrayFitCH();
    }

    // fDebug == 2 and fFitVoir no histo
    if (fDebug == 2) {
      if (fFitVoir < numberofbinsexpected) {
	AliInfo(Form("We will see the fit of the object %d",fFitVoir));
      }
      else {
	AliInfo("fFitVoir is out of range of the histo!");
	return kFALSE;
      }
    }

    // fDebug == 3  or 4 and fDet
    if (fDebug >= 3) {
      if ((fCalibraMode->GetNz(0) == 0) && (fCalibraMode->GetNrphi(0) == 0)) {
	AliInfo("Do you really want to see one detector without pad groups?");
	return kFALSE;
      }
      else {
	AliInfo(Form("You will see the detector: iPlane %d, iChamb %d, iSect %d"
                    ,fDet[0],fDet[1],fDet[2]));
	// A little geometry:
	Int_t rowMax = parCom->GetRowMax(fDet[0],fDet[1],fDet[2]);
	Int_t colMax = parCom->GetColMax(fDet[0]);
	// Create the histos to visualise
	CreateFitHistoCHDB(rowMax,colMax);
	if (fDebug == 4) {
          InitArrayFitCH();
	}
      }
    }

  }
    
  // PH and T0
  if (i == 1) {
    
    // Number of Xbins (detectors or groups of pads) if vector2d
    // Quick verification that we are not out of range!
    if (fCalibraVector) {
      if ((nbins == 0) && 
          (fCalibraVector->GetVectorPH()->GetEntriesFast()      > 0) && 
          ((Int_t) fCalibraVector->GetPlaPH()->GetEntriesFast() > 0)) {
	if ((Int_t) fCalibraVector->GetVectorPH()->GetEntriesFast() > numberofbinsexpected) {
	  AliInfo("ph doesn't correspond to the mode of pad group calibration!");
	  return kFALSE;
	}
	if ((Int_t) fCalibraVector->GetVectorPH()->GetEntriesFast() != 
            (Int_t) fCalibraVector->GetPlaPH()->GetEntriesFast()) {
	  AliInfo("VectorPH doesn't correspond to PlaPH!");
	  return kFALSE;
	}
      }
    }
        
    // Init tree
    InitTreePH();
    InitTreeT0();    

    //
    // Debugging: Create the histos
    //

    // fDebug == 0 nothing
    
    // fDebug == 1 
    if (fDebug == 1) {
      // Create the histos replique de ph
      InitArrayFitPH();
      InitArrayFitT0();
    }

    // fDebug == 2 and fFitVoir no histo
    if (fDebug == 2) {
      if (fFitVoir < numberofbinsexpected) {
	AliInfo(Form("We will see the fit of the object %d",fFitVoir));
      }
      else {
	AliInfo("fFitVoir is out of range of the histo!");
	return kFALSE;
      }
    }

    // fDebug == 3  or 4 and fDet
    if (fDebug >= 3) {
      if ((fCalibraMode->GetNz(1)    == 0) && 
          (fCalibraMode->GetNrphi(1) == 0)) {
	AliInfo("Do you really want to see one detector without pad groups?");
	return kFALSE;
      }
      else  {
	AliInfo(Form("You will see the detector: iPlane %d, iChamb %d, iSect %d"
                    ,fDet[0],fDet[1],fDet[2]));
	// A little geometry:
	Int_t rowMax = parCom->GetRowMax(fDet[0],fDet[1],fDet[2]);
	Int_t colMax = parCom->GetColMax(fDet[0]);
	// Create the histos to visualise
	CreateFitHistoPHDB(rowMax,colMax);
	CreateFitHistoT0DB(rowMax,colMax);
	if (fDebug == 4) {
	  InitArrayFitPH();
	  InitArrayFitT0();
	}
      }
    }

  }

  // PRF
  if (i == 2) {
    
    // Number of Xbins(detectors or groups of pads) if vector2d
    if (fCalibraVector){
      if ((nbins == 0) && 
          (fCalibraVector->GetVectorPRF()->GetEntriesFast() > 0) && 
          (fCalibraVector->GetPlaPRF()->GetEntriesFast()    > 0)) {
	// Quick verification that we are not out of range!
	if ((Int_t) fCalibraVector->GetVectorPRF()->GetEntriesFast() > numberofbinsexpected) {
	  AliInfo("ch doesn't correspond to the mode of pad group calibration!");
	  return kFALSE;
	}
	if ((Int_t) fCalibraVector->GetVectorPRF()->GetEntriesFast() != 
            (Int_t) fCalibraVector->GetPlaPRF()->GetEntriesFast()) {
	  AliInfo("VectorPRF doesn't correspond to PlaCH!");
	  return kFALSE;
	}
      }
    }
    
    // Init tree
    InitTreePRF();

    //
    // Debugging: Create the histos
    //

    // fDebug == 0 nothing

    // fDebug == 1 
    if (fDebug == 1) {
      // Create the histos replique de ch
      InitArrayFitPRF();
    }
    
    // fDebug == 2 and fFitVoir no histo
    if (fDebug == 2) {
      if (fFitVoir < numberofbinsexpected) {
	AliInfo(Form("We will see the fit of the object %d",fFitVoir));
      }
      else {
	AliInfo("fFitVoir is out of range of the histo!");
	return kFALSE;
      }
    }

    // fDebug == 3  or 4 and fDet
    if (fDebug >= 3) {
      if ((fCalibraMode->GetNz(2)    == 0) && 
          (fCalibraMode->GetNrphi(2) == 0)) {
	AliInfo("Do you really want to see one detector without pad groups?");
	return kFALSE;
      }
      else {
	AliInfo(Form("You will see the detector: iPlane %d, iChamb %d, iSect %d"
                    ,fDet[0],fDet[1],fDet[2]));
	// A little geometry:
	Int_t rowMax = parCom->GetRowMax(fDet[0],fDet[1],fDet[2]);
	Int_t colMax = parCom->GetColMax(fDet[0]);
	// Create the histos to visualise
	CreateFitHistoPRFDB(rowMax,colMax);
	if (fDebug == 4) {
          InitArrayFitPRF();
	}
      }
    }

  }

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
  fCountDet[i] = -1; // Current detector
  fCount[i]    =  0; // To find the next detector
  
  // If fDebug >= 3
  if (fDebug >= 3) {

    // Set countdet to the detector
    fCountDet[i] = AliTRDgeometry::GetDetector(fDet[0],fDet[1],fDet[2]);
        
    // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
    fCalibraMode->ModePadCalibration(fDet[1],i);
    fCalibraMode->ModePadFragmentation(fDet[0],fDet[1],fDet[2],i);
       
    // Set counter to write at the end of the detector
    fCount[i] = fDect1[i] + fCalibraMode->GetNfragZ(i)*fCalibraMode->GetNfragRphi(i);

  }

}

//____________Functions for initialising the AliTRDCalibraFit in the code_________
void AliTRDCalibraFit::UpdatefCountDetAndfCount(Int_t idect, Int_t i)
{
  //
  // See if we are in a new detector and update the
  // variables fNfragZ and fNfragRphi if yes 
  //

  // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
  // If fDebug == 1 or 0
  if ((fDebug == 0) || 
      (fDebug == 1)) {

    if (fCount[i] == idect) {
      
      // On en est au detector
      fCountDet[i] += 1;
      
      // Determination of fNnZ, fNnRphi, fNfragZ and fNfragRphi
     fCalibraMode->ModePadCalibration((Int_t) GetChamber(fCountDet[i]),i);
     fCalibraMode->ModePadFragmentation((Int_t) GetPlane(fCountDet[i])
                          ,(Int_t) GetChamber(fCountDet[i])
                          ,(Int_t) GetSector(fCountDet[i]),i);

      // Set for the next detector
      fCount[i] += fCalibraMode->GetNfragZ(i)*fCalibraMode->GetNfragRphi(i);

    }

  }

}

//____________Functions for initialising the AliTRDCalibraFit in the code_________
void AliTRDCalibraFit::ReconstructFitRowMinRowMax(Int_t idect, Int_t i)
{
  //
  // Reconstruct the min pad row, max pad row, min pad col and
  // max pad col of the calibration group for the Fit functions
  //

  if (fDebug <  2) {
    fCalibraMode->ReconstructionRowPadGroup((Int_t) (idect-(fCount[i]-(fCalibraMode->GetNfragZ(i)
                                                                      *fCalibraMode->GetNfragRphi(i)))),i);
  }
  if (fDebug >= 3) {
    fCalibraMode->ReconstructionRowPadGroup((Int_t) (idect-fDect1[i]),i);
  }

}

//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::NotEnoughStatistic(Int_t idect, Int_t i)
{
  //
  // For the case where there are not enough entries in the histograms
  // of the calibration group, the value present in the choosen database
  // will be put. A negativ sign enables to know that a fit was not possible.
  //
  
  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }
  
  // Get cal
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

  if (fDebug != 2) {
    AliInfo(Form("The element %d in this detector %d has not enough statistic to be fitted"
                ,idect-(fCount[i]-(fCalibraMode->GetNfragZ(i)*fCalibraMode->GetNfragRphi(i))),fCountDet[i]));
  }
  if (fDebug == 2) {
    AliInfo("The element has not enough statistic to be fitted");
  }

  if ((i == 0) && (fDebug != 2)) {
    
    // Calcul the coef from the database choosen
    CalculChargeCoefMean(fCountDet[0],(Int_t) (idect-fDect1[0]),kFALSE);
    
    // Fill the coefCH[2304] with negative value to say: not fitted
    AliInfo(Form("The row min %d, the row max %d, the colmin %d and the col max %d"
                ,fCalibraMode->GetRowMin(0)
                ,fCalibraMode->GetRowMax(0)
                ,fCalibraMode->GetColMin(0)
                ,fCalibraMode->GetColMax(0)));
    for (Int_t k = fCalibraMode->GetRowMin(0); k < fCalibraMode->GetRowMax(0); k++) {
      for (Int_t j = fCalibraMode->GetColMin(0); j < fCalibraMode->GetColMax(0); j++) {
	if (GetChamber(fCountDet[0]) == 2) {
          fCoefCH[(Int_t)(j*12+k)] = -TMath::Abs(fChargeCoef[3]);
	}
	if (GetChamber(fCountDet[0]) != 2) {
          fCoefCH[(Int_t)(j*16+k)] = -TMath::Abs(fChargeCoef[3]);
	}
      }
    }

    // Put the default value negative 
    if ((fDebug == 1) || 
        (fDebug == 4)) {

      if (fFitChargeBisOn) {
	fCoefCharge[2][idect-fDect1[0]]=-TMath::Abs(fChargeCoef[3]);
      }
      if (fMeanChargeOn) {
	fCoefCharge[1][idect-fDect1[0]]=-TMath::Abs(fChargeCoef[3]);
      }
      if(fFitChargeOn){
	fCoefCharge[0][idect-fDect1[0]]=-TMath::Abs(fChargeCoef[3]);
      }

    }
      
    // End of one detector
    if ((idect == (fCount[0]-1))) {
      FillVectorFitCH((Int_t) fCountDet[0]);
      // Reset
      for (Int_t k = 0; k < 2304; k++) {
	fCoefCH[k] = 0.0;
      }
    }

  }
  
  if ((i == 1) && (fDebug != 2)) {

    CalculVdriftCoefMean(fCountDet[1],(Int_t) (idect-fDect1[1]));
    CalculT0CoefMean(fCountDet[1],(Int_t) (idect-fDect1[1]));

    // Put the default value (time0 can be negativ, so we stay with + )
    if ((fDebug == 1) || 
        (fDebug == 4)) {

      if (fFitPHOn) {
	fCoefVdrift[0][(idect-fDect1[1])] = -fVdriftCoef[2];
	fCoefT0[0][(idect-fDect1[1])] = fT0Coef[2];
      }

      if(fFitPol2On) {
	fCoefVdrift[1][(idect-fDect1[1])] = -fVdriftCoef[2];
	fCoefT0[1][(idect-fDect1[1])] = fT0Coef[2];
      }
      if(fFitLagrPolOn) {
	fCoefVdrift[3][(idect-fDect1[1])] = -fVdriftCoef[2];
	fCoefT0[3][(idect-fDect1[1])] = fT0Coef[2];
      }

    }
    
    // Put the default value
    if (fDebug >= 3) {
      if(fFitPHOn){
	fVdriftCoef[0] = fVdriftCoef[2];
	fT0Coef[0]     = fT0Coef[2];
      }
      if(fFitPol2On){
	fVdriftCoef[1] = fVdriftCoef[2];
	fT0Coef[1]     = fT0Coef[2];
      }
      if(fFitLagrPolOn){
	fVdriftCoef[3] = fVdriftCoef[2];
	fT0Coef[3]     = fT0Coef[2];
      }
      FillCoefVdriftDB();
      FillCoefT0DB();
    }

    // Fill the tree if end of a detector.
    // The pointer to the branch stays with the default value negative!!!
    // PH
    // Pointer to the branch
    for (Int_t k = fCalibraMode->GetRowMin(1); k < fCalibraMode->GetRowMax(1); k++) {
      for (Int_t j = fCalibraMode->GetColMin(1); j < fCalibraMode->GetColMax(1); j++) {
	if (GetChamber(fCountDet[1]) == 2) {
          fVdriftPad[(Int_t)(j*12+k)] = -TMath::Abs(fVdriftCoef[2]);
	}
	if (GetChamber(fCountDet[1]) != 2) {
          fVdriftPad[(Int_t)(j*16+k)] = -TMath::Abs(fVdriftCoef[2]);
	}
      }
    }

    // End of one detector
    if ((idect == (fCount[1]-1)) && (fDebug != 2)) {
      FillTreeVdrift((Int_t) fCountDet[1]);
    }

    // T0
    // Fill the tree if end of a detector.
    // The pointer to the branch stays with the default value positive!!!
    // Pointer to the branch
    for (Int_t k = fCalibraMode->GetRowMin(1); k < fCalibraMode->GetRowMax(1); k++) {
      for (Int_t j = fCalibraMode->GetColMin(1); j < fCalibraMode->GetColMax(1); j++) {
	if (GetChamber(fCountDet[1]) == 2) {
          fT0Pad[(Int_t)(j*12+k)] = fT0Coef[2];
	}
	if (GetChamber(fCountDet[1]) != 2) {
          fT0Pad[(Int_t)(j*16+k)] = fT0Coef[2];
	}
      }
    }

    // End of one detector
    if ((idect == (fCount[1]-1)) && (fDebug != 2)) {
      FillTreeT0((Int_t) fCountDet[1]);
    }

  }

  if ((i == 2) && (fDebug != 2)) {

    CalculPRFCoefMean(fCountDet[2],(Int_t) (idect-fDect1[2]));
      
    if ((fDebug == 1) || 
        (fDebug == 4)) {
      if(fFitPRFOn){
	fCoefPRF[0][(idect-fDect1[2])] = -fPRFCoef[1];
      }
      if(fRMSPRFOn){
	fCoefPRF[2][(idect-fDect1[2])] = -fPRFCoef[1];
      }
    }

    if (fDebug >= 3){
      if(fFitPRFOn){
	fPRFCoef[0] = fPRFCoef[1];
      }
      if(fRMSPRFOn){
	fPRFCoef[2] = fPRFCoef[1];
      }
      FillCoefPRFDB();
    }

    // Fill the tree if end of a detector.
    // The pointer to the branch stays with the default value 1.5!!!
    // Pointer to the branch
    for (Int_t k = fCalibraMode->GetRowMin(2); k < fCalibraMode->GetRowMax(2); k++) {
      for (Int_t j = fCalibraMode->GetColMin(2); j < fCalibraMode->GetColMax(2); j++) {
	if((parCom->GetColMax(GetPlane(fCountDet[2])) != (j+1)) && (j != 0)){
	  if (GetChamber(fCountDet[2]) == 2) {
            fPRFPad[(Int_t)(j*12+k)] = -fPRFCoef[1];
	  }
	  if (GetChamber(fCountDet[2]) != 2) {
            fPRFPad[(Int_t)(j*16+k)] = -fPRFCoef[1];
	  }
	}
	else {
	  if (fAccCDB) {
	    if (GetChamber(fCountDet[2]) == 2) {
              fPRFPad[(Int_t)(j*12+k)] = -((Float_t) cal->GetPRFWidth(fCountDet[2],j,k));
	    }
	    if (GetChamber(fCountDet[2]) != 2) {
              fPRFPad[(Int_t)(j*16+k)] = -((Float_t) cal->GetPRFWidth(fCountDet[2],j,k));
	    }
	  }
	  if (!fAccCDB) {
	    if (GetChamber(fCountDet[2]) == 2) {
              fPRFPad[(Int_t)(j*12+k)] = -((Float_t) GetPRFDefault(GetPlane(fCountDet[2])));
	     }
	    if (GetChamber(fCountDet[2]) != 2) {
              fPRFPad[(Int_t)(j*16+k)] = -((Float_t) GetPRFDefault(GetPlane(fCountDet[2])));
	    }
	  }
	}
      }
    }

    // End of one detector
    if ((idect == (fCount[2]-1)) && (fDebug != 2)) {
      FillTreePRF((Int_t) fCountDet[2]);
    }

  }
  
  return kTRUE;

}

//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::FillInfosFit(Int_t idect, Int_t i)
{
  //
  // Fill the coefficients found with the fits or other
  // methods from the Fit functions
  //

  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }

  // Get cal
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

  if ((i == 0) && (fDebug != 2)) {
    // Fill the coefCH[2304] with fChargeCoef[0]
    // that would be negativ only if the fit failed totally
    //printf("for fCountDet %d we have %f\n",fCountDet[0],fChargeCoef[fFitChargeNDB]);
    //printf("RowMin %d RowMax %d ColMin %d ColMax %d\n",fCalibraMode->GetRowMin(0),fCalibraMode->GetRowMax(0),fCalibraMode->GetColMin(0),fCalibraMode->GetColMax(0));
    for (Int_t k = fCalibraMode->GetRowMin(0); k < fCalibraMode->GetRowMax(0); k++) {
      for (Int_t j = fCalibraMode->GetColMin(0); j < fCalibraMode->GetColMax(0); j++) {
	if (GetChamber(fCountDet[0]) == 2) {
          fCoefCH[(Int_t)(j*12+k)] = fChargeCoef[fFitChargeNDB];
	}
	if (GetChamber(fCountDet[0]) != 2) {
          fCoefCH[(Int_t)(j*16+k)] = fChargeCoef[fFitChargeNDB];
	}
      }
    }                
    // End of one detector
    if ((idect == (fCount[0]-1))) {
      FillVectorFitCH((Int_t) fCountDet[0]);
      // Reset
      for (Int_t k = 0; k < 2304; k++) {
	fCoefCH[k] = 0.0;
      }
    }
  }

  if ((i == 1) && (fDebug != 2)) {

    // PH
    // Pointer to the branch: fVdriftCoef[1] will ne negativ only if the fit failed totally 
    for (Int_t k = fCalibraMode->GetRowMin(1); k < fCalibraMode->GetRowMax(1); k++) {
      for (Int_t j = fCalibraMode->GetColMin(1); j < fCalibraMode->GetColMax(1); j++) {
	if (GetChamber(fCountDet[1]) == 2) {
          fVdriftPad[(Int_t)(j*12+k)]=fVdriftCoef[fFitPHNDB];
	}
	if (GetChamber(fCountDet[1]) != 2) {
          fVdriftPad[(Int_t)(j*16+k)]=fVdriftCoef[fFitPHNDB];
	}
      }
    }                
    // End of one detector
    if ((idect == (fCount[1]-1)) && (fDebug != 2)) {
      FillTreeVdrift((Int_t) fCountDet[1]);
    }

    // T0
    // Pointer to the branch: fT0Coef[1] will ne negativ only if the fit failed totally 
    for (Int_t k = fCalibraMode->GetRowMin(1); k < fCalibraMode->GetRowMax(1); k++) {
      for (Int_t j = fCalibraMode->GetColMin(1); j < fCalibraMode->GetColMax(1); j++) {
	if (GetChamber(fCountDet[1]) == 2) {
          fT0Pad[(Int_t)(j*12+k)]=fT0Coef[fFitPHNDB];
	}
	if (GetChamber(fCountDet[1]) != 2) {
          fT0Pad[(Int_t)(j*16+k)]=fT0Coef[fFitPHNDB];
	}
      }
    }                
    // End of one detector
    if ((idect == (fCount[1]-1)) && (fDebug != 2)) {
      FillTreeT0((Int_t) fCountDet[1]);
    }

  }

  if ((i == 2) && (fDebug != 2)) {
    // Pointer to the branch
    for (Int_t k = fCalibraMode->GetRowMin(2); k < fCalibraMode->GetRowMax(2); k++) {
      for (Int_t j = fCalibraMode->GetColMin(2); j < fCalibraMode->GetColMax(2); j++) {
	if ((parCom->GetColMax(GetPlane(fCountDet[2])) != (j+1)) && (j != 0)) {
	  if (GetChamber(fCountDet[2]) == 2) {
            fPRFPad[(Int_t)(j*12+k)] = fPRFCoef[fFitPRFNDB];
	  }
	  if (GetChamber(fCountDet[2]) != 2) {
            fPRFPad[(Int_t)(j*16+k)] = fPRFCoef[fFitPRFNDB];
	  }
	}
	else {
	  if (fAccCDB) {
	    if (GetChamber(fCountDet[2]) == 2) {
              fPRFPad[(Int_t)(j*12+k)] = (Float_t) cal->GetPRFWidth(fCountDet[2],j,k);
	    }
	    if (GetChamber(fCountDet[2]) != 2) {
              fPRFPad[(Int_t)(j*16+k)] = (Float_t) cal->GetPRFWidth(fCountDet[2],j,k);
	    }
	  }
	  if (!fAccCDB) {
	    if (GetChamber(fCountDet[2]) == 2) {
              fPRFPad[(Int_t)(j*12+k)] = (Float_t) GetPRFDefault(GetPlane(fCountDet[2]));
	    }
	    if (GetChamber(fCountDet[2]) != 2) {
              fPRFPad[(Int_t)(j*16+k)] = (Float_t) GetPRFDefault(GetPlane(fCountDet[2])); 
	    }
	  }
	}
      }
    }
    // End of one detector
    if ((idect == (fCount[2]-1)) && (fDebug != 2)) {
      FillTreePRF((Int_t) fCountDet[2]);
    }
  }

  return kTRUE;

}

//____________Functions for initialising the AliTRDCalibraFit in the code_________
Bool_t AliTRDCalibraFit::WriteFitInfos(Int_t i)
{
  //
  // In the case the user wants to write a file with a tree of the found
  // coefficients for the calibration before putting them in the database
  //

  TFile *fout = TFile::Open(fWriteNameCoef,"UPDATE");
  // Check if the file could be opened
  if (!fout || !fout->IsOpen()) {
    AliInfo("No File found!");
    return kFALSE;
  }

  if ((i == 0) && (fDebug != 2)) {
    // The DB stuff
    if ((fDebug == 4) || 
        (fDebug == 3)) {
      WriteCHDB(fout);
    }
    // The tree
    fout->WriteTObject(fGain,fGain->GetName(),(Option_t *) "writedelete");
  }

  if ((i == 1) && (fDebug != 2)) {
    // The DB stuff
    if ((fDebug == 4) || 
        (fDebug == 3)) {
      WritePHDB(fout);
    }
    // The tree
    fout->WriteTObject(fVdrift,fVdrift->GetName(),(Option_t *) "writedelete");
    // The DB stuff
    if ((fDebug == 4) || 
        (fDebug == 3)) {
      WriteT0DB(fout);
    }
    // The tree
    fout->WriteTObject(fT0,fT0->GetName(),(Option_t *) "writedelete");
  }

  if ((i == 2) && (fDebug != 2)) {
    // The DB stuff
    if ((fDebug == 4) || 
        (fDebug == 3)) {
      WritePRFDB(fout);
    }
    // The tree
    fout->WriteTObject(fPRF,fPRF->GetName(),(Option_t *) "writedelete");
  }

  fout->Close();

  return kTRUE;

}

//
//____________Fill Coef DB in case of visualisation of one detector____________
//

//_____________________________________________________________________________
void AliTRDCalibraFit::FillCoefVdriftDB()
{
  //
  // Fill the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //
 
  for (Int_t row = fCalibraMode->GetRowMin(1); row < fCalibraMode->GetRowMax(1); row++) {
    for (Int_t col = fCalibraMode->GetColMin(1); col < fCalibraMode->GetColMax(1); col++) {
      if(fFitPol2On){
	fCoefVdriftDB[1]->SetBinContent(row+1,col+1,TMath::Abs(fVdriftCoef[1]));
      }
      if (fFitPHOn ) {
        fCoefVdriftDB[0]->SetBinContent(row+1,col+1,TMath::Abs(fVdriftCoef[0]));
      }
      if (fFitLagrPolOn ) {
        fCoefVdriftDB[2]->SetBinContent(row+1,col+1,TMath::Abs(fVdriftCoef[3]));
      }
    }
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::FillCoefT0DB()
{
  //
  // Fill the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //
 
  for (Int_t row = fCalibraMode->GetRowMin(1); row < fCalibraMode->GetRowMax(1); row++) {
    for (Int_t col = fCalibraMode->GetColMin(1); col < fCalibraMode->GetColMax(1); col++) {
      if(fFitPol2On){
	fCoefT0DB[1]->SetBinContent(row+1,col+1,TMath::Abs(fT0Coef[1]));
      }
      if (fFitPHOn) {
        fCoefT0DB[0]->SetBinContent(row+1,col+1,TMath::Abs(fT0Coef[0]));
      }
      if (fFitLagrPolOn) {
	fCoefT0DB[2]->SetBinContent(row+1,col+1,TMath::Abs(fT0Coef[3]));
      }
    }
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::FillCoefChargeDB()
{
  //
  // Fill the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //

  for (Int_t row = fCalibraMode->GetRowMin(0); row < fCalibraMode->GetRowMax(0); row++) {
    for (Int_t col = fCalibraMode->GetColMin(0); col < fCalibraMode->GetColMax(0); col++) {
      if (fMeanChargeOn) {
        fCoefChargeDB[1]->SetBinContent(row+1,col+1,TMath::Abs(fChargeCoef[1]));
      }
      if (fFitChargeBisOn) {
        fCoefChargeDB[2]->SetBinContent(row+1,col+1,TMath::Abs(fChargeCoef[2]));
      }
      if(fFitChargeOn){
	fCoefChargeDB[0]->SetBinContent(row+1,col+1,TMath::Abs(fChargeCoef[0]));
      }
      if(fFitMeanWOn){
	fCoefChargeDB[3]->SetBinContent(row+1,col+1,TMath::Abs(fChargeCoef[4]));
      }
    }
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::FillCoefPRFDB()
{
  //
  // Fill the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //
  if(fFitPRFOn){
    for (Int_t row = fCalibraMode->GetRowMin(2); row < fCalibraMode->GetRowMax(2); row++) {
      for (Int_t col = fCalibraMode->GetColMin(2); col < fCalibraMode->GetColMax(2); col++) {
	fCoefPRFDB[0]->SetBinContent(row+1,col+1,fPRFCoef[0]);
      }
    }
  }
  if(fRMSPRFOn){
    for (Int_t row = fCalibraMode->GetRowMin(2); row < fCalibraMode->GetRowMax(2); row++) {
      for (Int_t col = fCalibraMode->GetColMin(2); col < fCalibraMode->GetColMax(2); col++) {
	fCoefPRFDB[1]->SetBinContent(row+1,col+1,fPRFCoef[2]);
      }
    }
  }

}

//
//____________Plot histos CoefPRF....__________________________________________
//

//_____________________________________________________________________________
void AliTRDCalibraFit::PlotWriteCH()
{
  //
  // Scale the coefficients to one, create the graph errors and write them if wanted
  //

  //TObjArray of the grapherrors and so on
  TObjArray *listofgraphs = new TObjArray(); 

  Int_t nbins = fDect2[0]-fDect1[0];

  // Scale the coefs
  // We will check fScaleFitFactor for the fFitChargeNDB, otherwise we calculate and normalise to 1
  // It can be that fScaleFitFactor is different from scale if we have taken a no default database as reference
  //

  //counter
  Int_t counter[4];
  counter[0] = 0; //how many groups are fitted for 0
  counter[1] = 0; //how many groups are with mean for 1
  counter[2] = 0; //how many groups are fitted for 2
  counter[3] = 0; //how many groups are fitted for 4
  Double_t sum = 0.0;
  Double_t scale = 1.0; 

  // Scale the histo
  // Is -1 if no fit or mean, is 1 if fit or mean
  Double_t *xValuesFitted = new Double_t[nbins]; 
  Double_t *xValuesFittedMean = new Double_t[nbins];
  Double_t *xValuesFittedBis =  new Double_t[nbins];
  Double_t *xValuesFittedMeanW =  new Double_t[nbins];
  for(Int_t k = 0; k < nbins; k ++){
    xValuesFitted[k] = -1;
    xValuesFittedMean[k] = -1;
    xValuesFittedMeanW[k] = -1;
    xValuesFittedBis[k] = -1;
  }

  if(fFitChargeOn){
    sum = 0.0;
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefCharge[0][l] > 0){
	sum += fCoefCharge[0][l];
	xValuesFitted[counter[0]]= l;
	counter[0]++;
      }
    }
    scale = 1.0;
    if(sum > 0.0) scale = counter[0]/sum;
    if(fFitChargeNDB == 0){
      if(scale != fScaleFitFactor){
	AliInfo(Form("The normalisation is different from a nomalisation to one."));
	AliInfo(Form("For one we have %f and here %f",scale,fScaleFitFactor));
	if(!fAccCDB) {
	  AliInfo(Form("It is not normal because we didn't choose a reference database!"));
	}
      }
      scale = fScaleFitFactor;
    }
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefCharge[0][l] > 0){
	fCoefCharge[0][l]=fCoefCharge[0][l]*scale;
	fCoefChargeE[0][l]=fCoefChargeE[0][l]*scale;
      }
    }
  }

  if(fFitMeanWOn){
    sum = 0.0;
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefCharge[4][l] > 0){
	sum += fCoefCharge[4][l];
	xValuesFittedMeanW[counter[3]]= l;
	counter[3]++;
      }
    }
    scale = 1.0;
    if(sum > 0.0) scale = counter[3]/sum;
    if(fFitChargeNDB == 4){
      if(scale != fScaleFitFactor){
	AliInfo(Form("The normalisation is different from a nomalisation to one."));
	AliInfo(Form("For one we have %f and here %f",scale,fScaleFitFactor));
	if(!fAccCDB) {
	  AliInfo(Form("It is not normal because we didn't choose a reference database!"));
	}
      }
      scale = fScaleFitFactor;
    }
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefCharge[4][l] > 0){
	fCoefCharge[4][l]=fCoefCharge[4][l]*scale;
	fCoefChargeE[3][l]=fCoefChargeE[3][l]*scale;
      }
    }
  }

  if(fMeanChargeOn){
    sum = 0.0;
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefCharge[1][l] > 0){
	sum += fCoefCharge[1][l];
	xValuesFittedMean[counter[1]]= l;
	counter[1]++;
      }
    }
    scale = 1.0;
    if(sum > 0.0) scale = counter[1]/sum;
    if(fFitChargeNDB == 1){
      if(scale != fScaleFitFactor){
	AliInfo(Form("The normalisation is different from a nomalisation to one."));
	AliInfo(Form("For one we have %f and here %f",scale,fScaleFitFactor));
	if(!fAccCDB) {
	  AliInfo(Form("It is not normal because we didn't choose a reference database!"));
	}
      }
      scale = fScaleFitFactor;
    }
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefCharge[1][l] > 0){
	fCoefCharge[1][l]=fCoefCharge[1][l]*scale;
	fCoefChargeE[1][l]=fCoefChargeE[1][l]*scale;
      }
    }
  }

  if(fFitChargeBisOn){
    sum = 0.0;
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefCharge[2][l] > 0){
	sum += fCoefCharge[2][l];
	xValuesFittedBis[counter[2]]= l;
	counter[2]++;
      }
    }
    scale = 1.0;
    if(sum > 0.0) scale = counter[2]/sum;
    if(fFitChargeNDB == 0){
      if(scale != fScaleFitFactor){
	AliInfo(Form("The normalisation is different from a nomalisation to one."));
	AliInfo(Form("For one we have %f and here %f",scale,fScaleFitFactor));
	if(!fAccCDB) {
	  AliInfo(Form("It is not normal because we didn't choose a reference database!"));
	}
      }
      scale = fScaleFitFactor;
    }
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefCharge[2][l] > 0){
	fCoefCharge[2][l]=fCoefCharge[2][l]*scale;
	fCoefChargeE[2][l]=fCoefChargeE[2][l]*scale;
      }
    }
  }
  
  // Create the X and Xerror
  Double_t *xValues = new Double_t[nbins];
  Double_t *xValuesE = new Double_t[nbins];
  for(Int_t k = 0; k < nbins; k ++){
    xValues[k] = k;
    xValuesE[k] = 0.0;
  }

  // Create the graph erros and plot them
  TCanvas *cch1 = new TCanvas("cch1","",50,50,600,800);
  cch1->cd();
  TLegend *legch1 = new TLegend(0.4,0.6,0.89,0.89);

  TGraph *graphCharge3  = new TGraph(nbins,xValues,fCoefCharge[3]);
  graphCharge3->SetName("coefcharge3");
  graphCharge3->SetTitle("");
  graphCharge3->GetXaxis()->SetTitle("Det/Pad groups");
  graphCharge3->GetYaxis()->SetTitle("gain factor");
  graphCharge3->SetLineColor(4);
  graphCharge3->SetMarkerStyle(25);
  graphCharge3->SetMarkerColor(4);
  listofgraphs->Add((TObject *)graphCharge3);
  legch1->AddEntry(graphCharge3,"f_{g} simulated","p");
  graphCharge3->Draw("AP");

  if (fFitChargeOn) {
    TGraphErrors *graphCharge0  = new TGraphErrors(nbins,xValues,fCoefCharge[0],xValuesE,fCoefChargeE[0]);
    graphCharge0->SetName("coefcharge0");
    graphCharge0->SetTitle("");
    graphCharge0->GetXaxis()->SetTitle("Det/Pad groups");
    graphCharge0->GetYaxis()->SetTitle("gain factor");
    graphCharge0->SetMarkerColor(6);
    graphCharge0->SetLineColor(6);
    graphCharge0->SetMarkerStyle(26);
    listofgraphs->Add((TObject *)graphCharge0); 
    legch1->AddEntry(graphCharge0,"f_{g} fit","p"); 
    graphCharge0->Draw("P");
  }
  if (fFitMeanWOn) {
    TGraphErrors *graphCharge4  = new TGraphErrors(nbins,xValues,fCoefCharge[4],xValuesE,fCoefChargeE[3]);
    graphCharge4->SetName("coefcharge4");
    graphCharge4->SetTitle("");
    graphCharge4->GetXaxis()->SetTitle("Det/Pad groups");
    graphCharge4->GetYaxis()->SetTitle("gain factor");
    graphCharge4->SetMarkerColor(1);
    graphCharge4->SetLineColor(1);
    graphCharge4->SetMarkerStyle(30);
    listofgraphs->Add((TObject *)graphCharge4); 
    legch1->AddEntry(graphCharge4,"f_{g} Mean W","p"); 
    graphCharge4->Draw("P");
  }
  if (fMeanChargeOn) {
    TGraphErrors *graphCharge1  = new TGraphErrors(nbins,xValues,fCoefCharge[1],xValuesE,fCoefChargeE[1]);
    graphCharge1->SetName("coefcharge1");
    graphCharge1->GetXaxis()->SetTitle("Det/Pad groups");
    graphCharge1->GetYaxis()->SetTitle("gain factor");
    graphCharge1->SetTitle("");
    graphCharge1->SetMarkerColor(2);
    graphCharge1->SetLineColor(2);
    graphCharge1->SetMarkerStyle(24);
    legch1->AddEntry(graphCharge1,"f_{g} mean","p");
    graphCharge1->Draw("P");
    listofgraphs->Add((TObject *)graphCharge1);
  } 
  if (fFitChargeBisOn ) {
    TGraphErrors *graphCharge2  = new TGraphErrors(nbins,xValues,fCoefCharge[2],xValuesE,fCoefChargeE[2]);
    graphCharge2->SetName("coefcharge2");
    graphCharge2->SetTitle("");
    graphCharge2->GetXaxis()->SetTitle("Det/Pad groups");
    graphCharge2->GetYaxis()->SetTitle("gain factor");
    graphCharge2->SetMarkerColor(8);
    graphCharge2->SetLineColor(8);
    graphCharge2->SetMarkerStyle(25);
    legch1->AddEntry(graphCharge2,"f_{g} fitbis","p");
    graphCharge2->Draw("P");
    listofgraphs->Add((TObject *)graphCharge2);
  } 
  legch1->Draw("same");
  
  //Create the arrays and the graphs for the delta
  Int_t thefirst = 0;
  TCanvas *cch2 = new TCanvas("cch2","",50,50,600,800);
  cch2->Divide(2,1);
  TLegend *legch3 = new TLegend(0.4,0.6,0.89,0.89);
  TLegend *legch2 = new TLegend(0.4,0.6,0.89,0.89);

  if(fFitChargeOn){
    cch2->cd(2);
    Double_t *yValuesDelta = new Double_t[counter[0]];
    for(Int_t k = 0; k < counter[0]; k++){
      if (fCoefCharge[3][(Int_t)(xValuesFitted[k])] > 0.0) {
	yValuesDelta[k] = (fCoefCharge[0][(Int_t)xValuesFitted[k]]-fCoefCharge[3][(Int_t)xValuesFitted[k]])
                        / fCoefCharge[3][(Int_t)xValuesFitted[k]];
      }
      else {
        yValuesDelta[k] = 0.0;
      }
    }
    TGraph *graphDeltaCharge0 = new TGraph(counter[0],&xValuesFitted[0],yValuesDelta);
    graphDeltaCharge0->SetName("deltacharge0");
    graphDeltaCharge0->GetXaxis()->SetTitle("Det/Pad groups");
    graphDeltaCharge0->GetYaxis()->SetTitle("#Deltag/g_{sim}");
    graphDeltaCharge0->SetMarkerColor(6);
    graphDeltaCharge0->SetTitle("");
    graphDeltaCharge0->SetLineColor(6);
    graphDeltaCharge0->SetMarkerStyle(26);
    listofgraphs->Add((TObject *)graphDeltaCharge0);   
    legch3->AddEntry(graphDeltaCharge0,"fit","p");
    graphDeltaCharge0->Draw("AP");
   
    cch2->cd(1); 
    TH1I *histoErrorCharge0 = new TH1I("errorcharge0","",100  ,-0.10,0.10);
    histoErrorCharge0->SetXTitle("#Deltag/g_{sim}");
    histoErrorCharge0->SetYTitle("counts"); 
    histoErrorCharge0->SetLineColor(6);
    histoErrorCharge0->SetLineStyle(1);
    histoErrorCharge0->SetStats(0);
    Double_t maxvalue = 0.0;
    for(Int_t k = 0; k < counter[0]; k++){
      histoErrorCharge0->Fill(yValuesDelta[k]);
      if(k == 0) maxvalue = TMath::Abs(yValuesDelta[k]);
      if(maxvalue < (TMath::Abs(yValuesDelta[k]))) maxvalue = TMath::Abs(yValuesDelta[k]);
    } 
    AliInfo(Form("The maximum deviation found dor the fit method is %f",maxvalue));
    legch2->AddEntry(histoErrorCharge0,"f_{g} fit","l");
    histoErrorCharge0->Draw();
    listofgraphs->Add((TObject *)histoErrorCharge0); 
    thefirst =1;
  }

  if(fFitMeanWOn){
    cch2->cd(2);
    Double_t *yValuesDelta = new Double_t[counter[3]];
    for(Int_t k = 0; k < counter[3]; k++){
      if (fCoefCharge[3][(Int_t)(xValuesFittedMeanW[k])] > 0.0) {
	yValuesDelta[k] = (fCoefCharge[4][(Int_t)xValuesFittedMeanW[k]]-fCoefCharge[3][(Int_t)xValuesFittedMeanW[k]])
                       / fCoefCharge[3][(Int_t)xValuesFittedMeanW[k]];
      }
      else {
        yValuesDelta[k] = 0.0;
      }
    }
    TGraph *graphDeltaCharge4 = new TGraph(counter[3],&xValuesFittedMeanW[0],yValuesDelta);
    graphDeltaCharge4->SetName("deltacharge4");
    graphDeltaCharge4->GetXaxis()->SetTitle("Det/Pad groups");
    graphDeltaCharge4->GetYaxis()->SetTitle("#Deltag/g_{sim}");
    graphDeltaCharge4->SetMarkerColor(1);
    graphDeltaCharge4->SetTitle("");
    graphDeltaCharge4->SetLineColor(1);
    graphDeltaCharge4->SetMarkerStyle(30);
    listofgraphs->Add((TObject *)graphDeltaCharge4);   
    legch3->AddEntry(graphDeltaCharge4,"Mean W","p");
    if(thefirst == 0){
      graphDeltaCharge4->Draw("AP");
    }
    else {
      graphDeltaCharge4->Draw("P");
    }
   
    cch2->cd(1); 
    TH1I *histoErrorCharge4 = new TH1I("errorcharge4","",100  ,-0.10,0.10);
    histoErrorCharge4->SetXTitle("#Deltag/g_{sim}");
    histoErrorCharge4->SetYTitle("counts"); 
    histoErrorCharge4->SetLineColor(1);
    histoErrorCharge4->SetLineStyle(1);
    histoErrorCharge4->SetStats(0);
    Double_t maxvalue = 0.0;
    for(Int_t k = 0; k < counter[3]; k++){
      histoErrorCharge4->Fill(yValuesDelta[k]);
      if(k == 0) maxvalue = yValuesDelta[k];
      if(maxvalue < (TMath::Abs(yValuesDelta[k]))) maxvalue = TMath::Abs(yValuesDelta[k]);
    } 
    AliInfo(Form("The maximum deviation found for the meanW method is %f",maxvalue));
    legch2->AddEntry(histoErrorCharge4,"f_{g} Mean W","l");
    if(thefirst == 0){
      histoErrorCharge4->Draw();
    }
    else {
      histoErrorCharge4->Draw("same");
    }
    listofgraphs->Add((TObject *)histoErrorCharge4); 
    thefirst =1;
  }

  if (fMeanChargeOn) {
    cch2->cd(2);
    Double_t *yValuesDeltaMean = new Double_t[counter[1]];
    for (Int_t k = 0; k < counter[1]; k++){
      if (fCoefCharge[3][(Int_t)xValuesFittedMean[k]] > 0.0) {
	yValuesDeltaMean[k] = (fCoefCharge[1][(Int_t)xValuesFittedMean[k]]-fCoefCharge[3][(Int_t)xValuesFittedMean[k]])
                            / fCoefCharge[3][(Int_t)xValuesFittedMean[k]];
      }
      else {
        yValuesDeltaMean[k] = 0.0;
      }
    }
    TGraph *graphDeltaCharge1 = new TGraph(counter[1],&xValuesFittedMean[0],yValuesDeltaMean);
    graphDeltaCharge1->SetName("deltacharge1");
    graphDeltaCharge1->GetXaxis()->SetTitle("Det/Pad groups");
    graphDeltaCharge1->GetYaxis()->SetTitle("#Deltag/g_{sim}");
    graphDeltaCharge1->SetMarkerColor(2);
    graphDeltaCharge1->SetMarkerStyle(24);
    graphDeltaCharge1->SetLineColor(2);
    graphDeltaCharge1->SetTitle("");
    legch3->AddEntry(graphDeltaCharge1,"mean","p");
    if(thefirst == 0){
      graphDeltaCharge1->Draw("AP");
    }
    else {
      graphDeltaCharge1->Draw("P");
    }
    listofgraphs->Add((TObject *)graphDeltaCharge1);
    
    cch2->cd(1);
    TH1I *histoErrorCharge1 = new TH1I("errorcharge1","",100  ,-0.10,0.10);
    histoErrorCharge1->SetXTitle("#Deltag/g_{sim}");
    histoErrorCharge1->SetYTitle("counts");
    histoErrorCharge1->SetLineColor(2);
    histoErrorCharge1->SetLineStyle(2);
    histoErrorCharge1->SetStats(0); 
    Double_t maxvalue = 0.0;
    for(Int_t k = 0; k < counter[1]; k++){
      histoErrorCharge1->Fill(yValuesDeltaMean[k]);
      if(k == 0) maxvalue = TMath::Abs(yValuesDeltaMean[k]);
      if(maxvalue < (TMath::Abs(yValuesDeltaMean[k]))) maxvalue = TMath::Abs(yValuesDeltaMean[k]);
    }
    AliInfo(Form("The maximum deviation found for the mean method is %f",maxvalue));
    legch2->AddEntry(histoErrorCharge1,"f_{g} mean","l");
    if(thefirst == 0){
      histoErrorCharge1->Draw();
    }
    else {
      histoErrorCharge1->Draw("same");
    }
    listofgraphs->Add((TObject *)histoErrorCharge1);
    thefirst = 1;
  }
  
  if (fFitChargeBisOn) {
    cch2->cd(2);
    Double_t *yValuesDeltaBis = new Double_t[counter[2]];
    for(Int_t k = 0; k < counter[2]; k++){
      if (fCoefCharge[3][(Int_t)xValuesFittedBis[k]] > 0.0) {
	yValuesDeltaBis[k] = (fCoefCharge[2][(Int_t)xValuesFittedBis[k]]-fCoefCharge[3][(Int_t)xValuesFittedBis[k]])
                           / fCoefCharge[3][(Int_t)xValuesFittedBis[k]];
      }
      else {
        yValuesDeltaBis[k] = 0.0;
      }
    }
    TGraph *graphDeltaCharge2 = new TGraph(counter[2],&xValuesFittedBis[0],yValuesDeltaBis);
    graphDeltaCharge2->SetName("deltacharge2");
    graphDeltaCharge2->GetXaxis()->SetTitle("Det/Pad groups");
    graphDeltaCharge2->GetYaxis()->SetTitle("#Deltag/g_{sim}");
    graphDeltaCharge2->SetMarkerColor(8);
    graphDeltaCharge2->SetLineColor(8);
    graphDeltaCharge2->SetMarkerStyle(25);
    legch3->AddEntry(graphDeltaCharge2,"fit","p");
    graphDeltaCharge2->SetTitle("");
    if(thefirst == 0){
      graphDeltaCharge2->Draw("AP");
    }
    else {
      graphDeltaCharge2->Draw("P");
    }
    listofgraphs->Add((TObject *)graphDeltaCharge2);

    cch2->cd(1);
    TH1I *histoErrorCharge2 = new TH1I("errorcharge2","",100  ,-0.10, 0.10);
    histoErrorCharge2->SetXTitle("#Deltag/g_{sim}");
    histoErrorCharge2->SetYTitle("counts"); 
    histoErrorCharge2->SetLineColor(8);
    histoErrorCharge2->SetLineStyle(5);
    histoErrorCharge2->SetLineWidth(3);
    histoErrorCharge2->SetStats(0);
    Double_t maxvalue = 0.0;
    for(Int_t k = 0; k < counter[2]; k++){
      histoErrorCharge2->Fill(yValuesDeltaBis[k]);
      if(k == 0) maxvalue = TMath::Abs(yValuesDeltaBis[k]);
      if(maxvalue < (TMath::Abs(yValuesDeltaBis[k]))) maxvalue = TMath::Abs(yValuesDeltaBis[k]);
    }
    AliInfo(Form("The maximum deviation found for the fit bis method is %f",maxvalue));
    legch2->AddEntry(histoErrorCharge2,"f_{g} fitbis","l");
    if(thefirst == 0){
      histoErrorCharge2->Draw();
    }
    else {
      histoErrorCharge2->Draw("same");
    }
    listofgraphs->Add((TObject *)histoErrorCharge2);
    //it doesn't matter anymore but...
    thefirst = 1;
  }

  cch2->cd(2);
  legch3->Draw("same");
  cch2->cd(1);
  legch2->Draw("same");

  //Write if wanted
  if (fWriteCoef[0]){
    TFile *fout = TFile::Open(fWriteNameCoef,"UPDATE");
    // Check if the file could be opened
    if (!fout || !fout->IsOpen()) {
      AliInfo("No File found!");
    }
    
    else{
      for(Int_t k = 0; k < listofgraphs->GetEntriesFast(); k++){
	fout->WriteTObject((TObject *) listofgraphs->At(k),((TObject *)listofgraphs->At(k))->GetName()
                          ,(Option_t *) "OverWrite");
      }
    }
    fout->Close();
  }
    
}

//_____________________________________________________________________________
void AliTRDCalibraFit::PlotWritePH()
{
  //
  // create the graph errors and write them if wanted
  //

  //TObjArray of the grapherrors and so on
  TObjArray *listofgraphs = new TObjArray(); 

  Int_t nbins = fDect2[1]-fDect1[1];

  //See the number of fitted for delta
  
  //counter
  Int_t counter[3];
  counter[0] = 0;
  counter[1] = 0;
  counter[2] = 0;

  Double_t *xValuesFitted = new Double_t[nbins];
  Double_t *xValuesFittedPH = new Double_t[nbins];
  Double_t *xValuesFittedLP = new Double_t[nbins];
  for(Int_t k = 0; k < nbins; k ++){
    xValuesFitted[k] = -1;
    xValuesFittedPH[k] = -1;
    xValuesFittedLP[k] = -1;
  }

  if(fFitPol2On){
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefVdrift[1][l] > 0){
	xValuesFitted[counter[1]]=l;
	counter[1]++;
      }
    }
  }
  if(fFitLagrPolOn){
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefVdrift[3][l] > 0){
	xValuesFittedLP[counter[2]]=l;
	counter[2]++;
      }
    }
  }
  if(fFitPHOn){
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefVdrift[0][l] > 0){
	xValuesFittedPH[counter[0]]= l;
	counter[0]++;
      }
    }
  }

  //Create the X and Xerror
  Double_t *xValues = new Double_t[nbins];
  Double_t *xValuesE = new Double_t[nbins];
  for(Int_t k = 0; k < nbins; k ++){
    xValues[k] = k;
    xValuesE[k] = 0.0;
  }

  //Create the graph erros and plot them
  TCanvas *cph1 = new TCanvas("cph1","",50,50,600,800);
  cph1->cd();

  TGraph *graphVdrift2  = new TGraph(nbins,xValues,fCoefVdrift[2]);
  graphVdrift2->SetName("coefvdrift2");
  graphVdrift2->SetTitle("");
  graphVdrift2->GetXaxis()->SetTitle("Det/Pad groups");
  graphVdrift2->GetYaxis()->SetTitle("Vdrift [cm/#mus]");
  graphVdrift2->SetLineColor(4);
  listofgraphs->Add((TObject *)graphVdrift2);
  TLegend *legph1 = new TLegend(0.4,0.6,0.89,0.89);
  legph1->AddEntry(graphVdrift2,"Vdrift simulated","l");
  graphVdrift2->Draw("AL");

  if(fFitPol2On){
    TGraphErrors *graphVdrift1  = new TGraphErrors(nbins,xValues,fCoefVdrift[1],xValuesE,fCoefVdriftE[1]);
    graphVdrift1->SetName("coefvdrift1");
    graphVdrift1->SetTitle("");
    graphVdrift1->GetXaxis()->SetTitle("Det/Pad groups");
    graphVdrift1->GetYaxis()->SetTitle("Vdrift [cm/#mus]");
    graphVdrift1->SetMarkerColor(6);
    graphVdrift1->SetLineColor(6);
    graphVdrift1->SetMarkerStyle(26);
    listofgraphs->Add((TObject *)graphVdrift1);
    legph1->AddEntry(graphVdrift1,"Vdrift fit","p");
    graphVdrift1->Draw("P");
  }
  if (fFitPHOn) {
    TGraphErrors *graphVdrift0  = new TGraphErrors(nbins,xValues,fCoefVdrift[0],xValuesE,fCoefVdriftE[0]);
    graphVdrift0->SetName("coefVdrift0");
    graphVdrift0->GetXaxis()->SetTitle("Det/Pad groups");
    graphVdrift0->GetYaxis()->SetTitle("Vdrift [cm/#mus]");
    graphVdrift0->SetTitle("");
    graphVdrift0->SetMarkerColor(2);
    graphVdrift0->SetLineColor(2);
    graphVdrift0->SetMarkerStyle(24);
    legph1->AddEntry(graphVdrift0,"v_{fit PH}","p");
    graphVdrift0->Draw("P");
    listofgraphs->Add((TObject *)graphVdrift0);
  } 
  if (fFitLagrPolOn) {
    TGraphErrors *graphVdrift3  = new TGraphErrors(nbins,xValues,fCoefVdrift[3],xValuesE,fCoefVdriftE[2]);
    graphVdrift3->SetName("coefVdrift3");
    graphVdrift3->GetXaxis()->SetTitle("Det/Pad groups");
    graphVdrift3->GetYaxis()->SetTitle("Vdrift [cm/#mus]");
    graphVdrift3->SetTitle("");
    graphVdrift3->SetMarkerColor(1);
    graphVdrift3->SetLineColor(1);
    graphVdrift3->SetMarkerStyle(28);
    legph1->AddEntry(graphVdrift3,"v_{LagrPol}","p");
    graphVdrift3->Draw("P");
    listofgraphs->Add((TObject *)graphVdrift3);
  } 
  legph1->Draw("same");

  //Create the arrays and the graphs for the delta
  TCanvas *cph2 = new TCanvas("cph2","",50,50,600,800);
  cph2->Divide(2,1);
  TLegend *legph3 = new TLegend(0.4,0.6,0.89,0.89);
  TLegend *legph2 = new TLegend(0.4,0.6,0.89,0.89);
  Int_t first = 0;
 
  if(fFitPol2On){
    cph2->cd(2);
    Double_t *yValuesDelta = new Double_t[counter[1]];
    for (Int_t k = 0; k < counter[1]; k++){
      if (fCoefVdrift[2][(Int_t)(xValuesFitted[k])] > 0.0) {
	yValuesDelta[k] = (fCoefVdrift[1][(Int_t)xValuesFitted[k]]-fCoefVdrift[2][(Int_t)xValuesFitted[k]])
                        / fCoefVdrift[2][(Int_t)xValuesFitted[k]];
      }
      else {
        yValuesDelta[k] = 0.0;
      }
    }
    TGraph *graphDeltaVdrift1 = new TGraph(counter[1],&xValuesFitted[0],yValuesDelta);
    graphDeltaVdrift1->SetName("deltavdrift1");
    graphDeltaVdrift1->GetXaxis()->SetTitle("Det/Pad groups");
    graphDeltaVdrift1->GetYaxis()->SetTitle("#Deltav/v_{sim}");
    graphDeltaVdrift1->SetMarkerColor(6);
    graphDeltaVdrift1->SetTitle("");
    graphDeltaVdrift1->SetLineColor(6);
    graphDeltaVdrift1->SetMarkerStyle(26);
    listofgraphs->Add((TObject *)graphDeltaVdrift1);
    legph3->AddEntry(graphDeltaVdrift1,"v_{slope method}","p");
    graphDeltaVdrift1->Draw("AP");
  
    cph2->cd(1); 
    TH1I *histoErrorVdrift1 = new TH1I("errorvdrift1","",100  ,-0.2,0.2);
    histoErrorVdrift1->SetXTitle("#Deltav/v_{sim}");
    histoErrorVdrift1->SetYTitle("counts"); 
    histoErrorVdrift1->SetLineColor(6);
    histoErrorVdrift1->SetLineStyle(1);
    histoErrorVdrift1->SetStats(0);
    Double_t maxvalue = 0.0;
    for(Int_t k = 0; k < counter[1]; k++){
      histoErrorVdrift1->Fill(yValuesDelta[k]);
      if(k == 0) maxvalue = yValuesDelta[k];
      if(maxvalue < (TMath::Abs(yValuesDelta[k]))) maxvalue = TMath::Abs(yValuesDelta[k]);
    }
    AliInfo(Form("The maximum deviation found for the Pol2 method is %f",maxvalue));
    legph2->AddEntry(histoErrorVdrift1,"v_{slope method}","l");
    histoErrorVdrift1->Draw();
    listofgraphs->Add((TObject *)histoErrorVdrift1);
    first = 1;
  }

  if (fFitPHOn) {
    cph2->cd(2);
    Double_t *yValuesDeltaPH = new Double_t[counter[0]];
    for(Int_t k = 0; k < counter[0]; k++){
      if(fCoefVdrift[2][(Int_t)xValuesFittedPH[k]] > 0.0) {
	yValuesDeltaPH[k] = (fCoefVdrift[0][(Int_t)xValuesFittedPH[k]]-fCoefVdrift[2][(Int_t)xValuesFittedPH[k]])/fCoefVdrift[2][(Int_t)xValuesFittedPH[k]];
      }
      else yValuesDeltaPH[k] = 0.0;
    }
    TGraph *graphDeltaVdrift0 = new TGraph(counter[0],&xValuesFittedPH[0],yValuesDeltaPH);
    graphDeltaVdrift0->SetName("deltavdrift0");
    graphDeltaVdrift0->GetXaxis()->SetTitle("Det/Pad groups");
    graphDeltaVdrift0->GetYaxis()->SetTitle("#Deltav/v_{sim}");
    graphDeltaVdrift0->SetMarkerColor(2);
    graphDeltaVdrift0->SetMarkerStyle(24);
    graphDeltaVdrift0->SetLineColor(2);
    graphDeltaVdrift0->SetTitle("");
    legph3->AddEntry(graphDeltaVdrift0,"v_{fit PH}","p");
    if(first){
      graphDeltaVdrift0->Draw("P");
    }
    else {
      graphDeltaVdrift0->Draw("AP");
    }
    listofgraphs->Add((TObject *)graphDeltaVdrift0);
    cph2->cd(1);
    TH1I *histoErrorVdrift0 = new TH1I("errorvdrift0","",100  ,-0.2,0.2);
    histoErrorVdrift0->SetXTitle("#Deltav/v_{sim}");
    histoErrorVdrift0->SetYTitle("counts");
    histoErrorVdrift0->SetLineColor(2);
    histoErrorVdrift0->SetLineStyle(2);
    histoErrorVdrift0->SetStats(0); 
    Double_t maxvalue = 0.0;
    for(Int_t k = 0; k < counter[0]; k++){
      histoErrorVdrift0->Fill(yValuesDeltaPH[k]);
      if(k == 0) maxvalue = yValuesDeltaPH[k];
      if(maxvalue < (TMath::Abs(yValuesDeltaPH[k]))) maxvalue = TMath::Abs(yValuesDeltaPH[k]);
    }
    AliInfo(Form("The maximum deviation found for the fit method is %f",maxvalue));
    legph2->AddEntry(histoErrorVdrift0,"v_{fit PH}","l");
    if(first){
      histoErrorVdrift0->Draw("same");
    }
    else {
      histoErrorVdrift0->Draw();
    }
    listofgraphs->Add((TObject *)histoErrorVdrift0);
    first = 1;
  }

  if (fFitLagrPolOn) {
    cph2->cd(2);
    Double_t *yValuesDeltaPH = new Double_t[counter[2]];
    for (Int_t k = 0; k < counter[2]; k++){
      if (fCoefVdrift[2][(Int_t)xValuesFittedLP[k]] > 0.0) {
	yValuesDeltaPH[k] = (fCoefVdrift[3][(Int_t)xValuesFittedLP[k]]-fCoefVdrift[2][(Int_t)xValuesFittedLP[k]])
                          / fCoefVdrift[2][(Int_t)xValuesFittedLP[k]];
      }
      else {
        yValuesDeltaPH[k] = 0.0;
      }
    }
    TGraph *graphDeltaVdrift3 = new TGraph(counter[2],&xValuesFittedLP[0],yValuesDeltaPH);
    graphDeltaVdrift3->SetName("deltavdrift3");
    graphDeltaVdrift3->GetXaxis()->SetTitle("Det/Pad groups");
    graphDeltaVdrift3->GetYaxis()->SetTitle("#Deltav/v_{sim}");
    graphDeltaVdrift3->SetMarkerColor(1);
    graphDeltaVdrift3->SetMarkerStyle(28);
    graphDeltaVdrift3->SetLineColor(1);
    graphDeltaVdrift3->SetTitle("");
    legph3->AddEntry(graphDeltaVdrift3,"v_{LagrPol}","p");
    if(first){
      graphDeltaVdrift3->Draw("P");
    }
    else {
      graphDeltaVdrift3->Draw("AP");
    }
    listofgraphs->Add((TObject *)graphDeltaVdrift3);
    cph2->cd(1);
    TH1I *histoErrorVdrift3 = new TH1I("errorvdrift3","",100  ,-0.2,0.2);
    histoErrorVdrift3->SetXTitle("#Deltav/v_{sim}");
    histoErrorVdrift3->SetYTitle("counts");
    histoErrorVdrift3->SetLineColor(1);
    histoErrorVdrift3->SetLineStyle(1);
    histoErrorVdrift3->SetStats(0); 
    Double_t maxvalue = 0.0;
    for(Int_t k = 0; k < counter[2]; k++){
      histoErrorVdrift3->Fill(yValuesDeltaPH[k]);
      if(k == 0) maxvalue = yValuesDeltaPH[k];
      if(maxvalue < (TMath::Abs(yValuesDeltaPH[k]))) maxvalue = TMath::Abs(yValuesDeltaPH[k]);
    }
    AliInfo(Form("The maximum deviation found for the LagrPol method is %f",maxvalue));
    legph2->AddEntry(histoErrorVdrift3,"v_{LagrPol}","l");
    if(first){
      histoErrorVdrift3->Draw("same");
    }
    else {
      histoErrorVdrift3->Draw();
    }
    listofgraphs->Add((TObject *)histoErrorVdrift3);
    first = 1;
  }
  cph2->cd(2);
  legph3->Draw("same");
  cph2->cd(1);
  legph2->Draw("same");

  //Write if wanted
  if (fWriteCoef[1]){
    TFile *fout = TFile::Open(fWriteNameCoef,"UPDATE");
    // Check if the file could be opened
    if (!fout || !fout->IsOpen()) {
      AliInfo("No File found!");
    }
    
    else{
      for(Int_t k = 0; k < listofgraphs->GetEntriesFast(); k++){
	fout->WriteTObject((TObject *) listofgraphs->At(k),((TObject *)listofgraphs->At(k))->GetName()
                          ,(Option_t *) "OverWrite");
      }
    }
    fout->Close();
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::PlotWriteT0()
{
  //
  // create the graph errors and write them if wanted
  //

  //TObjArray of the grapherrors and so on
  TObjArray *listofgraphs = new TObjArray(); 

  Int_t nbins = fDect2[1]-fDect1[1];

  //See the number of fitted for delta: here T0 can be negative, we don't use the sign but the error
  //and the grapherrors of the coefficients contained the no fitted with error 0.0
  
  //counter
  Int_t counter[3];
  counter[0] = 0;
  counter[1] = 0;
  counter[2] = 0;

  Double_t *xValuesFitted = new Double_t[nbins];
  Double_t *xValuesFittedPH = new Double_t[nbins];
  Double_t *xValuesFittedLP = new Double_t[nbins];
  for(Int_t k = 0; k < nbins; k ++){
    xValuesFitted[k] = -1;
    xValuesFittedPH[k] = -1;
    xValuesFittedLP[k] = -1;
  }

  if(fFitPol2On){
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefT0E[1][l] != 0.0){
	xValuesFitted[counter[1]]=l;
	counter[1]++;
      }
    }
  }

  if(fFitPHOn){
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefT0E[0][l] != 0.0){
	xValuesFittedPH[counter[0]]= l;
	counter[0]++;
      }
    }
  }

  if(fFitLagrPolOn){
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefT0E[2][l] == 1.0){
	xValuesFittedLP[counter[2]]= l;
	counter[2]++;
      }
    }
  }

  //Create the X and Xerror
  Double_t *xValues = new Double_t[nbins];
  Double_t *xValuesE = new Double_t[nbins];
  for(Int_t k = 0; k < nbins; k ++){
    xValues[k] = k;
    xValuesE[k] = 0.0;
  }

  //Create the graph erros and plot them
  TCanvas *ct01 = new TCanvas("ct01","",50,50,600,800);
  ct01->cd();
  TLegend *legt01 = new TLegend(0.4,0.6,0.89,0.89);
 
  TGraph *graphT02  = new TGraph(nbins,xValues,fCoefT0[2]);
  graphT02->SetName("coeft02");
  graphT02->SetTitle("");
  graphT02->GetXaxis()->SetTitle("Det/Pad groups");
  graphT02->GetYaxis()->SetTitle("T0 [time bins]");
  graphT02->SetLineColor(4);
  listofgraphs->Add((TObject *)graphT02);
  legt01->AddEntry(graphT02,"T0 simulated","l");
  graphT02->Draw("AL");

  if(fFitPol2On){
    TGraphErrors *graphT01  = new TGraphErrors(nbins,xValues,fCoefT0[1],xValuesE,fCoefT0E[1]);
    graphT01->SetName("coeft01");
    graphT01->SetTitle("");
    graphT01->GetXaxis()->SetTitle("Det/Pad groups");
    graphT01->GetYaxis()->SetTitle("T0 [time bins]");
    graphT01->SetMarkerColor(6);
    graphT01->SetLineColor(6);
    graphT01->SetMarkerStyle(26);
    listofgraphs->Add((TObject *)graphT01);
    legt01->AddEntry(graphT01,"T0 slope method","p");
    graphT01->Draw("P");
  }
  if (fFitPHOn) {
    TGraphErrors *graphT00  = new TGraphErrors(nbins,xValues,fCoefT0[0],xValuesE,fCoefT0E[0]);
    graphT00->SetName("coeft00");
    graphT00->GetXaxis()->SetTitle("Det/Pad groups");
    graphT00->GetYaxis()->SetTitle("T0 [time bins]");
    graphT00->SetTitle("");
    graphT00->SetMarkerColor(2);
    graphT00->SetLineColor(2);
    graphT00->SetMarkerStyle(24);
    legt01->AddEntry(graphT00,"T0 fit","p");
    graphT00->Draw("P");
    listofgraphs->Add((TObject *)graphT00);
  } 
  if (fFitLagrPolOn) {
    TGraphErrors *graphT03  = new TGraphErrors(nbins,xValues,fCoefT0[3],xValuesE,xValuesE);
    graphT03->SetName("coeft03");
    graphT03->GetXaxis()->SetTitle("Det/Pad groups");
    graphT03->GetYaxis()->SetTitle("T0 [time bins]");
    graphT03->SetTitle("");
    graphT03->SetMarkerColor(1);
    graphT03->SetLineColor(1);
    graphT03->SetMarkerStyle(28);
    legt01->AddEntry(graphT03,"T0 LagrPol","p");
    graphT03->Draw("P");
    listofgraphs->Add((TObject *)graphT03);
  } 
  legt01->Draw("same");
  
  //Create the arrays and the graphs for the delta
  TCanvas *ct02 = new TCanvas("ct02","",50,50,600,800);
  ct02->Divide(2,1);
  TLegend *legt03 = new TLegend(0.4,0.6,0.89,0.89);
  TLegend *legt02 = new TLegend(0.4,0.6,0.89,0.89);
  Int_t first = 0;

  if(fFitPol2On){
    ct02->cd(2);
    Double_t *yValuesDelta = new Double_t[counter[1]];
    for(Int_t k = 0; k < counter[1]; k++){
      yValuesDelta[k] = (fCoefT0[1][(Int_t)xValuesFitted[k]]-fCoefT0[2][(Int_t)xValuesFitted[k]]);
    }
    TGraph *graphDeltaT01 = new TGraph(counter[1],&xValuesFitted[0],yValuesDelta);
    graphDeltaT01->SetName("deltat01");
    graphDeltaT01->GetXaxis()->SetTitle("Det/Pad groups");
    graphDeltaT01->GetYaxis()->SetTitle("#Deltat0 [time bins]");
    graphDeltaT01->SetMarkerColor(6);
    graphDeltaT01->SetTitle("");
    graphDeltaT01->SetLineColor(6);
    graphDeltaT01->SetMarkerStyle(26);
    listofgraphs->Add((TObject *)graphDeltaT01);
    legt03->AddEntry(graphDeltaT01,"T0_{slope method}","p");
    graphDeltaT01->Draw("AP");
 
    ct02->cd(1); 
    TH1I *histoErrorT01 = new TH1I("errort01","",100  ,-0.2,0.2);
    histoErrorT01->SetXTitle("#Deltat0 [time bins]");
    histoErrorT01->SetYTitle("counts"); 
    histoErrorT01->SetLineColor(6);
    histoErrorT01->SetLineStyle(1);
    histoErrorT01->SetStats(0);
    Double_t maxvalue = 0.0;
    for(Int_t k = 0; k < counter[1]; k++){
      histoErrorT01->Fill(yValuesDelta[k]);
      if(k == 0) maxvalue = yValuesDelta[k];
      if(maxvalue < (TMath::Abs(yValuesDelta[k]))) maxvalue = (TMath::Abs(yValuesDelta[k]));
    }
    AliInfo(Form("The maximum deviation found for the Pol2 method is %f",maxvalue));
    legt02->AddEntry(histoErrorT01,"T0_{slope method}","l");
    histoErrorT01->Draw();
    listofgraphs->Add((TObject *)histoErrorT01); 
    first = 1;
  }
  if (fFitPHOn) {
    ct02->cd(2);
    Double_t *yValuesDeltaPH = new Double_t[counter[0]];
    for(Int_t k = 0; k < counter[0]; k++){
      yValuesDeltaPH[k] = (fCoefT0[0][(Int_t)xValuesFittedPH[k]]-fCoefT0[2][(Int_t)xValuesFittedPH[k]]);
    }
    TGraph *graphDeltaT00 = new TGraph(counter[0],&xValuesFittedPH[0],yValuesDeltaPH);
    graphDeltaT00->SetName("deltat00");
    graphDeltaT00->GetXaxis()->SetTitle("Det/Pad groups");
    graphDeltaT00->GetYaxis()->SetTitle("#Deltat0 [time bins]");
    graphDeltaT00->SetMarkerColor(2);
    graphDeltaT00->SetMarkerStyle(24);
    graphDeltaT00->SetLineColor(2);
    graphDeltaT00->SetTitle("");
    legt03->AddEntry(graphDeltaT00,"T0_{fit PH}","p");
    if(first) {
      graphDeltaT00->Draw("P");
    }
    else{
      graphDeltaT00->Draw("AP");
    }
    listofgraphs->Add((TObject *)graphDeltaT00);
    ct02->cd(1);
    TH1I *histoErrorT00 = new TH1I("errort00","",100  ,-0.2,0.2);
    histoErrorT00->SetXTitle("#Deltat0 [time bins]");
    histoErrorT00->SetYTitle("counts");
    histoErrorT00->SetLineColor(2);
    histoErrorT00->SetLineStyle(2);
    histoErrorT00->SetStats(0); 
    Double_t maxvalue = 0.0;
    for(Int_t k = 0; k < counter[0]; k++){
      histoErrorT00->Fill(yValuesDeltaPH[k]);
      if(k == 0) maxvalue = yValuesDeltaPH[k];
      if(maxvalue < (TMath::Abs(yValuesDeltaPH[k]))) maxvalue = (TMath::Abs(yValuesDeltaPH[k]));
    }
    AliInfo(Form("The maximum deviation found for the fit method is %f",maxvalue));
    legt02->AddEntry(histoErrorT00,"T0_{fit PH}","l");
    if(first){
      histoErrorT00->Draw("same");
    }
    else{
      histoErrorT00->Draw();
    }
    listofgraphs->Add((TObject *)histoErrorT00);
    first = 1;
  }

  if (fFitLagrPolOn) {
    ct02->cd(2);
    Double_t *yValuesDeltaPH = new Double_t[counter[2]];
    for(Int_t k = 0; k < counter[2]; k++){
      yValuesDeltaPH[k] = (fCoefT0[3][(Int_t)xValuesFittedLP[k]]-fCoefT0[2][(Int_t)xValuesFittedLP[k]]);
    }
    TGraph *graphDeltaT03 = new TGraph(counter[2],&xValuesFittedLP[0],yValuesDeltaPH);
    graphDeltaT03->SetName("deltat03");
    graphDeltaT03->GetXaxis()->SetTitle("Det/Pad groups");
    graphDeltaT03->GetYaxis()->SetTitle("#Deltat0 [time bins]");
    graphDeltaT03->SetMarkerColor(1);
    graphDeltaT03->SetMarkerStyle(28);
    graphDeltaT03->SetLineColor(1);
    graphDeltaT03->SetTitle("");
    legt03->AddEntry(graphDeltaT03,"T0_{LagrPol}","p");
    if(first) {
      graphDeltaT03->Draw("P");
    }
    else{
      graphDeltaT03->Draw("AP");
    }
    listofgraphs->Add((TObject *)graphDeltaT03);
    ct02->cd(1);
    TH1I *histoErrorT03 = new TH1I("errort03","",100  ,-0.2,0.2);
    histoErrorT03->SetXTitle("#Deltat0 [time bins]");
    histoErrorT03->SetYTitle("counts");
    histoErrorT03->SetLineColor(1);
    histoErrorT03->SetLineStyle(1);
    histoErrorT03->SetStats(0); 
    Double_t maxvalue = 0.0;
    for(Int_t k = 0; k < counter[2]; k++){
      histoErrorT03->Fill(yValuesDeltaPH[k]);
      if(k == 0) maxvalue = yValuesDeltaPH[k];
      if(maxvalue < (TMath::Abs(yValuesDeltaPH[k]))) maxvalue = (TMath::Abs(yValuesDeltaPH[k]));
    }
    AliInfo(Form("The maximum deviation found for the LagrPol method is %f",maxvalue));
    legt02->AddEntry(histoErrorT03,"T0_{LagrPol}","l");
    if(first){
      histoErrorT03->Draw("same");
    }
    else{
      histoErrorT03->Draw();
    }
    listofgraphs->Add((TObject *)histoErrorT03);
    first = 1;
  }

  ct02->cd(2);
  legt03->Draw("same");
  ct02->cd(1);
  legt02->Draw("same");
  
  //Write if wanted
  if (fWriteCoef[1]){
    TFile *fout = TFile::Open(fWriteNameCoef,"UPDATE");
    // Check if the file could be opened
    if (!fout || !fout->IsOpen()) {
      AliInfo("No File found!");
    }
    
    else{
      for(Int_t k = 0; k < listofgraphs->GetEntriesFast(); k++){
	fout->WriteTObject((TObject *) listofgraphs->At(k),((TObject *)listofgraphs->At(k))->GetName()
                          ,(Option_t *) "OverWrite");
      }
    }
    fout->Close();
  }  

}

//_____________________________________________________________________________
void AliTRDCalibraFit::PlotWritePRF()
{
  //
  // create the graph errors and write them if wanted
  //

  //TObjArray of the grapherrors and so on
  TObjArray *listofgraphs = new TObjArray(); 

  Int_t nbins = fDect2[2]-fDect1[2];

  //See the number of fitted for delta
  
  //counter
  Int_t counter[2];
  counter[0] = 0;
  counter[1] = 0;

  Double_t *xValuesFitted = new Double_t[nbins];
  for(Int_t k = 0; k < nbins; k ++){
    xValuesFitted[k] = -1;
  }
  Double_t *xValuesRMS = new Double_t[nbins];
  for(Int_t k = 0; k < nbins; k ++){
    xValuesRMS[k] = -1;
  }

  if(fFitPRFOn){
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefPRF[0][l] > 0){
	xValuesFitted[counter[0]]=l;
	counter[0]++;
      }
    }
  }
  if(fRMSPRFOn){
    for(Int_t l = 0; l < nbins; l++){
      if(fCoefPRF[2][l] > 0){
	xValuesRMS[counter[1]]=l;
	counter[1]++;
      }
    }
  }
 

  //Create the X and Xerror
  Double_t *xValues = new Double_t[nbins];
  Double_t *xValuesE = new Double_t[nbins];
  for(Int_t k = 0; k < nbins; k ++){
    xValues[k] = k;
    xValuesE[k] = 0.0;
  }

  //Create the graph erros and plot them
  TCanvas *cprf1 = new TCanvas("cprf1","",50,50,600,800);
  cprf1->cd();
  TLegend *legprf1 = new TLegend(0.4,0.6,0.89,0.89);

  TGraph *graphPRF1  = new TGraph(nbins,xValues,fCoefPRF[1]);
  graphPRF1->SetName("coefprf1");
  graphPRF1->SetTitle("");
  graphPRF1->GetXaxis()->SetTitle("Det/Pad groups");
  graphPRF1->GetYaxis()->SetTitle("PRF width [p.u]");
  graphPRF1->SetLineColor(4);
  graphPRF1->SetMarkerColor(4);
  graphPRF1->SetMarkerStyle(25);
  graphPRF1->SetMarkerSize(0.7);
  listofgraphs->Add((TObject *)graphPRF1);
  legprf1->AddEntry(graphPRF1,"PRF width simulated","p");
  graphPRF1->Draw("AP");
  
  if(fFitPRFOn){
    TGraphErrors *graphPRF0  = new TGraphErrors(nbins,xValues,fCoefPRF[0],xValuesE,fCoefPRFE[0]);
    graphPRF0->SetName("coefprf0");
    graphPRF0->SetTitle("");
    graphPRF0->GetXaxis()->SetTitle("Det/Pad groups");
    graphPRF0->GetYaxis()->SetTitle("PRF Width [p.u]");
    graphPRF0->SetMarkerColor(6);
    graphPRF0->SetLineColor(6);
    graphPRF0->SetMarkerStyle(26);
    listofgraphs->Add((TObject *)graphPRF0);
    legprf1->AddEntry(graphPRF0,"PRF fit","p");
    graphPRF0->Draw("P");
  }
  if(fRMSPRFOn){
    TGraphErrors *graphPRF2  = new TGraphErrors(nbins,xValues,fCoefPRF[2],xValuesE,fCoefPRFE[1]);
    graphPRF2->SetName("coefprf2");
    graphPRF2->SetTitle("");
    graphPRF2->GetXaxis()->SetTitle("Det/Pad groups");
    graphPRF2->GetYaxis()->SetTitle("PRF Width [p.u]");
    graphPRF2->SetMarkerColor(1);
    graphPRF2->SetLineColor(1);
    graphPRF2->SetMarkerStyle(28);
    listofgraphs->Add((TObject *)graphPRF2);
    legprf1->AddEntry(graphPRF2,"PRF Rms","p");
    graphPRF2->Draw("P");
  }
  legprf1->Draw("same");

  
  //Create the arrays and the graphs for the delta
  TCanvas *cprf2 = new TCanvas("cprf2","",50,50,600,800);
  cprf2->Divide(2,1);
  Int_t first = 0;
  TLegend *legprf3 = new TLegend(0.4,0.6,0.89,0.89);
  TLegend *legprf2 = new TLegend(0.4,0.6,0.89,0.89);

  if(fFitPRFOn){
    cprf2->cd(2);
    Double_t *yValuesDelta = new Double_t[counter[0]];
    for(Int_t k = 0; k < counter[0]; k++){
      if(fCoefPRF[1][(Int_t)xValuesFitted[k]] > 0.0){
	yValuesDelta[k] = (fCoefPRF[0][(Int_t)xValuesFitted[k]]-fCoefPRF[1][(Int_t)xValuesFitted[k]])
                        / (fCoefPRF[1][(Int_t)xValuesFitted[k]]);
      }
    }
    TGraph *graphDeltaPRF0 = new TGraph(counter[0],&xValuesFitted[0],yValuesDelta);
    graphDeltaPRF0->SetName("deltaprf0");
    graphDeltaPRF0->GetXaxis()->SetTitle("Det/Pad groups");
    graphDeltaPRF0->GetYaxis()->SetTitle("#Delta#sigma/#sigma_{sim}");
    graphDeltaPRF0->SetMarkerColor(6);
    graphDeltaPRF0->SetTitle("");
    graphDeltaPRF0->SetLineColor(6);
    graphDeltaPRF0->SetMarkerStyle(26);
    listofgraphs->Add((TObject *)graphDeltaPRF0); 
    legprf3->AddEntry(graphDeltaPRF0,"#sigma_{fit}","p");
    graphDeltaPRF0->Draw("AP");
  
    cprf2->cd(1); 
    TH1I *histoErrorPRF0 = new TH1I("errorprf10","",100  ,-0.1,0.2);
    histoErrorPRF0->SetXTitle("#Delta#sigma/#sigma_{sim}");
    histoErrorPRF0->SetYTitle("counts"); 
    histoErrorPRF0->SetLineColor(6);
    histoErrorPRF0->SetLineStyle(1);
    histoErrorPRF0->SetStats(0);
    Double_t maxvalue = 0.0;
    for(Int_t k = 0; k < counter[0]; k++){
      histoErrorPRF0->Fill(yValuesDelta[k]);
      if(k == 0) maxvalue = yValuesDelta[k];
      if(maxvalue < (TMath::Abs(yValuesDelta[k]))) maxvalue = (TMath::Abs(yValuesDelta[k]));
    }  
    AliInfo(Form("The maximum deviation for the fit method is %f",maxvalue));
    legprf2->AddEntry(histoErrorPRF0,"#sigma_{fit}","l");
    histoErrorPRF0->Draw();
    listofgraphs->Add((TObject *)histoErrorPRF0); 
    first = 1;
  }

  if(fRMSPRFOn){
    cprf2->cd(2);
    Double_t *yValuesDelta = new Double_t[counter[1]];
    for(Int_t k = 0; k < counter[1]; k++){
      if(fCoefPRF[1][(Int_t)xValuesRMS[k]] > 0.0){
	yValuesDelta[k] = (fCoefPRF[2][(Int_t)xValuesRMS[k]]-fCoefPRF[1][(Int_t)xValuesRMS[k]])
                        / (fCoefPRF[1][(Int_t)xValuesRMS[k]]);
      }
    }
    TGraph *graphDeltaPRF2 = new TGraph(counter[1],&xValuesRMS[0],yValuesDelta);
    graphDeltaPRF2->SetName("deltaprf2");
    graphDeltaPRF2->GetXaxis()->SetTitle("Det/Pad groups");
    graphDeltaPRF2->GetYaxis()->SetTitle("#Delta#sigma/#sigma_{sim}");
    graphDeltaPRF2->SetMarkerColor(1);
    graphDeltaPRF2->SetTitle("");
    graphDeltaPRF2->SetLineColor(1);
    graphDeltaPRF2->SetMarkerStyle(28);
    listofgraphs->Add((TObject *)graphDeltaPRF2); 
    legprf3->AddEntry(graphDeltaPRF2,"#sigma_{rms}","p");
    if(first){
      graphDeltaPRF2->Draw("P");
    }
    else {
      graphDeltaPRF2->Draw("AP");
    }
  
    cprf2->cd(1); 
    TH1I *histoErrorPRF2 = new TH1I("errorprf12","",100  ,-0.1,0.2);
    histoErrorPRF2->SetXTitle("#Delta#sigma/#sigma_{sim}");
    histoErrorPRF2->SetYTitle("counts"); 
    histoErrorPRF2->SetLineColor(1);
    histoErrorPRF2->SetLineStyle(1);
    histoErrorPRF2->SetStats(0);
    Double_t maxvalue = 0.0;
    for(Int_t k = 0; k < counter[1]; k++){
      histoErrorPRF2->Fill(yValuesDelta[k]);
      if(k == 0) maxvalue = yValuesDelta[k];
      if(maxvalue < TMath::Abs(yValuesDelta[k])) maxvalue = TMath::Abs(yValuesDelta[k]);
    }
    AliInfo(Form("The maximum deviation for the rms is %f",maxvalue));  
    legprf2->AddEntry(histoErrorPRF2,"#sigma_{rms}","l");
    if(first){
      histoErrorPRF2->Draw("same");
    }
    else {
      histoErrorPRF2->Draw();
    }
    listofgraphs->Add((TObject *)histoErrorPRF2); 
    first = 1;
  }

  cprf2->cd(2);
  legprf3->Draw("same");
  cprf2->cd(1);
  legprf2->Draw("same");
  
  
  //Write if wanted
  if (fWriteCoef[2]){
    TFile *fout = TFile::Open(fWriteNameCoef,"UPDATE");
    // Check if the file could be opened
    if (!fout || !fout->IsOpen()) {
      AliInfo("No File found!");
    }
    
    else{
      for(Int_t k = 0; k < listofgraphs->GetEntriesFast(); k++){
	fout->WriteTObject((TObject *) listofgraphs->At(k),((TObject *)listofgraphs->At(k))->GetName()
                          ,(Option_t *) "OverWrite");
      }
    }
    fout->Close();
  }  
  
}

//
//____________Plot histos DB___________________________________________________
//

//_____________________________________________________________________________
void AliTRDCalibraFit::PlotCHDB()
{
  //
  // Plot the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //

  TCanvas *cchdb = new TCanvas("cchdb","",50,50,600,800);
  Int_t nb = 0;
  if(fFitChargeOn) nb++;
  if(fFitChargeBisOn) nb++;
  if(fMeanChargeOn) nb++;
  if(fFitMeanWOn) nb++;
  if(nb > 0){
    cchdb->Divide(nb,1);
    nb = 0;
    if(fMeanChargeOn){
      cchdb->cd(nb);
      fCoefChargeDB[1]->Draw("LEGO");
      nb++;
    }
    if(fFitChargeOn){
      cchdb->cd(nb);
      fCoefChargeDB[0]->Draw("LEGO");
      nb++;
    }
    if(fFitMeanWOn){
      cchdb->cd(nb);
      fCoefChargeDB[3]->Draw("LEGO");
      nb++;
    }
    if(fFitChargeBisOn){
      cchdb->cd(nb);
      fCoefChargeDB[2]->Draw("LEGO");
      //it doesn't matter anymore but....
      nb++;
    }
  }
}

//_____________________________________________________________________________
void AliTRDCalibraFit::PlotPHDB()
{
  //
  // Plot the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //

  TCanvas *cphdb = new TCanvas("cphdb","",50,50,600,800);
  Int_t nb = 0;
  if(fFitPol2On) nb++;
  if(fFitPHOn) nb++;
  if(fFitLagrPolOn) nb++;
  if(nb > 0){
    cphdb->Divide(nb,1);
    nb = 0;
    if(fFitPHOn){
      cphdb->cd(nb);
      fCoefVdriftDB[0]->Draw("LEGO");
      nb++;
    }
    if(fFitPol2On){
      cphdb->cd(nb);
      fCoefVdriftDB[1]->Draw("LEGO");
      nb++;
    }
    if(fFitLagrPolOn){
      cphdb->cd(nb);
      fCoefVdriftDB[2]->Draw("LEGO");
      nb++;
    }
  }
}

//_____________________________________________________________________________
void AliTRDCalibraFit::PlotT0DB()
{
  //
  // Plot the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //
  TCanvas *ct0db = new TCanvas("ct0db","",50,50,600,800);
  Int_t nb = 0;
  if(fFitPol2On) nb++;
  if(fFitPHOn) nb++;
  if(fFitLagrPolOn) nb++;
  if(nb > 0){
    ct0db->Divide(nb,1);
    nb = 0;
    if(fFitPHOn){
      ct0db->cd(nb);
      fCoefT0DB[0]->Draw("LEGO");
      nb++;
    }
    if(fFitPol2On){
      ct0db->cd(nb);
      fCoefT0DB[1]->Draw("LEGO");
      nb++;
    }
    if(fFitLagrPolOn){
      ct0db->cd(nb);
      fCoefT0DB[2]->Draw("LEGO");
      nb++;
    }
  }
}

//_____________________________________________________________________________
void AliTRDCalibraFit::PlotPRFDB()
{
  //
  // Plot the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //

  TCanvas *cprfdb = new TCanvas("cprfdb","",50,50,600,800);
  Int_t nb = 0;
  if(fFitPRFOn) nb++;
  if(fRMSPRFOn) nb++;
  if(nb > 0){
    cprfdb->Divide(nb,1);
    nb = 0;
    if(fFitPRFOn){
      cprfdb->cd(nb);
      fCoefPRFDB[0]->Draw("LEGO");
      nb++;
    }
    if(fRMSPRFOn){
      cprfdb->cd(nb);
      fCoefPRFDB[1]->Draw("LEGO");
      nb++;
    }
  }
}

//
//____________Write DB Histos__________________________________________________
//

//_____________________________________________________________________________
void AliTRDCalibraFit::WriteCHDB(TFile *fout)
{
  //
  // If wanted, write the debug histos for fDebug = 3 and fDebug = 4
  //
  if(fFitChargeOn){
    fout->WriteTObject(fCoefChargeDB[0],fCoefChargeDB[0]->GetName(),(Option_t *) "OverWrite");
  }
  if(fFitMeanWOn){
    fout->WriteTObject(fCoefChargeDB[3],fCoefChargeDB[0]->GetName(),(Option_t *) "OverWrite");
  }
  if (fMeanChargeOn) {
    fout->WriteTObject(fCoefChargeDB[1],fCoefChargeDB[1]->GetName(),(Option_t *) "OverWrite");
  }
  if (fFitChargeBisOn ) {
    fout->WriteTObject(fCoefChargeDB[2],fCoefChargeDB[2]->GetName(),(Option_t *) "OverWrite");
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::WritePHDB(TFile *fout)
{
  //
  // If wanted, write the debug histos for fDebug = 3 and fDebug = 4
  //

  if (fFitPHOn) {
    fout->WriteTObject(fCoefVdriftDB[0],fCoefVdriftDB[0]->GetName(),(Option_t *) "OverWrite");
  }
  if(fFitPol2On){
    fout->WriteTObject(fCoefVdriftDB[1],fCoefVdriftDB[1]->GetName(),(Option_t *) "OverWrite");
  }
  if(fFitLagrPolOn){
    fout->WriteTObject(fCoefVdriftDB[2],fCoefVdriftDB[2]->GetName(),(Option_t *) "OverWrite");
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::WriteT0DB(TFile *fout)
{
  //
  // If wanted, write the debug histos for fDebug = 3 and fDebug = 4
  //

  if (fFitPHOn) {
    fout->WriteTObject(fCoefT0DB[0],fCoefT0DB[0]->GetName(),(Option_t *) "OverWrite");
  }
  if(fFitPol2On){
    fout->WriteTObject(fCoefT0DB[1],fCoefT0DB[1]->GetName(),(Option_t *) "OverWrite");
  }
  if(fFitLagrPolOn){
    fout->WriteTObject(fCoefT0DB[2],fCoefT0DB[2]->GetName(),(Option_t *) "OverWrite");
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::WritePRFDB(TFile *fout)
{
  //
  // If wanted, write the debug histos for fDebug = 3 and fDebug = 4
  //
  if(fFitPRFOn){
    fout->WriteTObject(fCoefPRFDB[0],fCoefPRFDB[0]->GetName(),(Option_t *) "OverWrite");
  }
  if(fRMSPRFOn){
    fout->WriteTObject(fCoefPRFDB[1],fCoefPRFDB[1]->GetName(),(Option_t *) "OverWrite");
  }

}

//
//____________Calcul Coef Mean_________________________________________________
//

//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::CalculT0CoefMean(Int_t dect, Int_t idect)
{
  //
  // For the detector Dect calcul the mean time 0
  // for the calibration group idect from the choosen database
  //

  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB  Manager");
    return kFALSE;
  }

  fT0Coef[2] = 0.0;

  if ((fDebug != 2) && fAccCDB) {

    for (Int_t row = fCalibraMode->GetRowMin(1); row < fCalibraMode->GetRowMax(1); row++) {
      for (Int_t col = fCalibraMode->GetColMin(1); col < fCalibraMode->GetColMax(1); col++) {
	// Groups of pads
	if ((fCalibraMode->GetNz(1)    > 0)  || 
            (fCalibraMode->GetNrphi(1) > 0)) {
	  fT0Coef[2] += (Float_t) cal->GetT0(dect,col,row);
	}
	// Per detectors
	else {
          fT0Coef[2] += (Float_t) cal->GetT0Average(dect);
	}
      }
    }

    fT0Coef[2] = fT0Coef[2] / ((fCalibraMode->GetColMax(1)-fCalibraMode->GetColMin(1))
                             * (fCalibraMode->GetRowMax(1)-fCalibraMode->GetRowMin(1)));
    if ((fDebug == 1) || 
        (fDebug == 4)) {
      fCoefT0[2][idect] = fT0Coef[2];
    }

  }

  return kTRUE;
  
}

//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::CalculChargeCoefMean(Int_t dect, Int_t idect, Bool_t vrai)
{
  //
  // For the detector Dect calcul the mean gain factor
  // for the calibration group idect from the choosen database
  //

  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB  Manager");
    return kFALSE;
  }
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam  Manager");
    return kFALSE;
  }

  fChargeCoef[3] = 0.0;

  if (fDebug != 2) {

    for (Int_t row = fCalibraMode->GetRowMin(0); row < fCalibraMode->GetRowMax(0); row++) {
      for (Int_t col = fCalibraMode->GetColMin(0); col < fCalibraMode->GetColMax(0); col++) {
	// Groups of pads
	if ((fCalibraMode->GetNz(0)    > 0) || 
            (fCalibraMode->GetNrphi(0) > 0)) {
	  if (fAccCDB) {
            fChargeCoef[3] += (Float_t) cal->GetGainFactor(dect,col,row);
	  }
	  if (vrai && fAccCDB) {
            fScaleFitFactor += (Float_t) cal->GetGainFactor(dect,col,row);
	  }
	  if (!fAccCDB) {
            fChargeCoef[3] += 1.0;
	  }
	  if (vrai && (!fAccCDB)) {
            fScaleFitFactor += 1.0;
	  }
	}
	// Per detectors
	else {
	  if (fAccCDB) {
            fChargeCoef[3] += (Float_t) cal->GetGainFactorAverage(dect);
	  }
	  if (vrai && fAccCDB) {
            fScaleFitFactor += ((Float_t) cal->GetGainFactorAverage(dect));
	  }
	  if (!fAccCDB) {
            fChargeCoef[3] += 1.0;
	  }
	  if (vrai && (!fAccCDB)) {
            fScaleFitFactor += 1.0;
	  }
	}
      }
    }

    fChargeCoef[3] = fChargeCoef[3] / ((fCalibraMode->GetColMax(0)-fCalibraMode->GetColMin(0))
                                     * (fCalibraMode->GetRowMax(0)-fCalibraMode->GetRowMin(0)));
    if ((fDebug == 1) || 
        (fDebug == 4)) {
      fCoefCharge[3][idect]=fChargeCoef[3];
    }

  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::CalculPRFCoefMean(Int_t dect, Int_t idect)
{
  //
  // For the detector Dect calcul the mean sigma of pad response
  // function for the calibration group idect from the choosen database
  //

  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB  Manager");
    return kFALSE;
  }

  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam  Manager");
    return kFALSE;
  }

  fPRFCoef[1] = 0.0;
  Int_t cot = 0;

  if (fDebug != 2) {
    
    for (Int_t row = fCalibraMode->GetRowMin(2); row < fCalibraMode->GetRowMax(2); row++) {
      for (Int_t col = fCalibraMode->GetColMin(2); col < fCalibraMode->GetColMax(2); col++) {
	if ((parCom->GetColMax(GetPlane(dect)) != (col+1)) && (col != 0)) {
	  cot++;
	  if (fAccCDB) {
            fPRFCoef[1] += (Float_t) cal->GetPRFWidth(dect,col,row);
	  }
	  if (!fAccCDB) {
            fPRFCoef[1] += GetPRFDefault(GetPlane(dect));
	  }
	}
      }
    }

    if (cot > 0) {
      fPRFCoef[1] = fPRFCoef[1]/cot;
      if ((fDebug == 1) ||
          (fDebug == 4)) {
        fCoefPRF[1][idect] = fPRFCoef[1];
      }
    }
    if (cot <= 0) {
      if ((fDebug == 1) ||
          (fDebug == 4)) {
	if (fAccCDB) {
          fCoefPRF[1][idect] = cal->GetPRFWidth(dect,fCalibraMode->GetColMin(2),fCalibraMode->GetRowMin(2));
	}
	if (!fAccCDB) {
          fCoefPRF[1][idect] = GetPRFDefault(GetPlane(dect));
	}
      }
    }

  }

  return kTRUE;
  
}

//_____________________________________________________________________________
Bool_t AliTRDCalibraFit::CalculVdriftCoefMean(Int_t dect, Int_t idect)
{
  //
  // For the detector dect calcul the mean drift velocity for the
  // calibration group idect from the choosen database
  //

  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB  Manager");
    return kFALSE;
  }

  fVdriftCoef[2] = 0.0;

  if (fDebug != 2) {
    for (Int_t row = fCalibraMode->GetRowMin(1); row < fCalibraMode->GetRowMax(1); row++) {
      for (Int_t col = fCalibraMode->GetColMin(1); col < fCalibraMode->GetColMax(1); col++) {
	// Groups of pads
	if ((fCalibraMode->GetNz(1)    > 0) || 
            (fCalibraMode->GetNrphi(1) > 0)) {
	  if (fAccCDB) {
            fVdriftCoef[2] += (Float_t) cal->GetVdrift(dect,col,row);
	  }
	  if (!fAccCDB) {
            fVdriftCoef[2] += 1.5;
	  }
	}
	// Per detectors
	else {
	  if (fAccCDB) {
            fVdriftCoef[2] += (Float_t) cal->GetVdriftAverage(dect);
	  }
	  if (!fAccCDB) {
            fVdriftCoef[2] += 1.5;
	  }
	}
      }
    }
    fVdriftCoef[2] = fVdriftCoef[2] / ((fCalibraMode->GetColMax(1)-fCalibraMode->GetColMin(1))
                                     * (fCalibraMode->GetRowMax(1)-fCalibraMode->GetRowMin(1)));
    if ((fDebug == 1) || 
        (fDebug == 4)) {
      fCoefVdrift[2][idect] = fVdriftCoef[2];
    }
  }

  return kTRUE;
  
}

//_____________________________________________________________________________
Float_t AliTRDCalibraFit::GetPRFDefault(Int_t plane) const
{
  //
  // Default width of the PRF if there is no database as reference
  //

  if (plane == 0) {
    return 0.515;
  }
  if (plane == 1) {
    return 0.502;
  }
  if (plane == 2) {
    return 0.491;
  }
  if (plane == 3) {
    return 0.481;
  }
  if (plane == 4) {
    return 0.471;
  }
  if (plane == 5) {
    return 0.463;
  }
  else {
    return 0.0;
  }
  
}

//____________Fit Methods______________________________________________________

//_____________________________________________________________________________
void AliTRDCalibraFit::FitPente(TH1* projPH, Int_t idect)
{
  //
  // Slope methode for the drift velocity
  //
  
  // Constants
  const Float_t kDrWidth = AliTRDgeometry::DrThick();
  Int_t binmax   = 0;
  Int_t binmin   = 0;
  fPhd[0]        = 0.0;
  fPhd[1]        = 0.0;
  fPhd[2]        = 0.0;
  Int_t ju       = 0;
  Double_t vdriftCoefE = 0.0;
  Double_t t0CoefE = 0.0;
  fVdriftCoef[1] = 0.0;
  fT0Coef[1]     = 0.0;
  TLine *line = new TLine();

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
  if(fDebug == 2) AliInfo(Form("maximum positive bin for the positive slope %d",binmax));
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
  if (l3P2am != 0) {
    fPhd[0] = -(l3P1am / (2 * l3P2am));
  }
  if(!fTakeTheMaxPH){
    if((l3P1am != 0.0) && (l3P2am != 0.0)){
      t0CoefE = (l3P1amE/l3P1am + l3P2amE/l3P2am)*fPhd[0];
    }
  }
  if(fDebug == 2) AliInfo(Form("maximum extrapolated positive bin for the positive slope %f",fPhd[0]));
  
  // Amplification region
  binmax = 0;
  ju     = 0;
  for (Int_t kbin = 1; kbin < projPH->GetNbinsX(); kbin ++) {
    if (((projPH->GetBinContent(kbin+1) - projPH->GetBinContent(kbin)) <= 0.0) && 
         (ju == 0) && 
         (kbin > (fPhd[0]/widbins))) {
      binmax = kbin;
      ju     = 1;
    }
  }
  if(fDebug == 2) AliInfo(Form("For the amplification region binmax %d",binmax));
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
  
  if (l3P2amf != 0) {
    fPhd[1] = -(l3P1amf / (2 * l3P2amf));
  }
  if((l3P1amf != 0.0) && (l3P2amf != 0.0)){
    vdriftCoefE = (l3P1amfE/l3P1amf + l3P2amfE/l3P2amf)*fPhd[1];
  }
  if(fTakeTheMaxPH){
    t0CoefE = vdriftCoefE;
  }
  if(fDebug == 2) AliInfo(Form("For the amplification region extrapolated binmax %f",fPhd[1]));

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
  if(fDebug == 2) AliInfo(Form("For the drift region binmin %d",binmin));
  pente->Fit("pol2"
            ,"0MR"
            ,""
            ,TMath::Max(pente->GetBinCenter(binmin-1),             0.0)
            ,TMath::Min(pente->GetBinCenter(binmin+1),(Double_t) limit));
  Float_t l3P1dr = pente->GetFunction("pol2")->GetParameter(1);
  Float_t l3P2dr = pente->GetFunction("pol2")->GetParameter(2);
  Float_t l3P1drE = pente->GetFunction("pol2")->GetParError(1);
  Float_t l3P2drE = pente->GetFunction("pol2")->GetParError(2);
  if (l3P2dr != 0) {
    fPhd[2] = -(l3P1dr / (2 * l3P2dr));
  }
  if((l3P1dr != 0.0) && (l3P2dr != 0.0)){
    vdriftCoefE += (l3P1drE/l3P1dr + l3P2drE/l3P2dr)*fPhd[2]; 
  }
  if(fDebug == 2) AliInfo(Form("For the drift region extrapolated binmax %f",fPhd[2]));

  Float_t fPhdt0 = 0.0;
  if(fTakeTheMaxPH) fPhdt0 = fPhd[1];
  else fPhdt0 = fPhd[0];

  if ((fPhd[2] > fPhd[0]) && 
      (fPhd[2] > fPhd[1]) && 
      (fPhd[1] > fPhd[0]) &&
      (put)) {
    fVdriftCoef[1] = (kDrWidth) / (fPhd[2]-fPhd[1]);
    if(fFitPHNDB == 1) fNumberFitSuccess++;
    if (fPhdt0 >= 0.0) {
      fT0Coef[1] = (fPhdt0 - fT0Shift) / widbins;
      if (fT0Coef[1] < -1.0) {
        fT0Coef[1] = fT0Coef[2];
      }
    }
    else {
      fT0Coef[1] = fT0Coef[2];
    }
  }
  else {
    fVdriftCoef[1] = -TMath::Abs(fVdriftCoef[2]);
    fT0Coef[1]     = fT0Coef[2];
  }

  if ((fDebug == 1) || 
      (fDebug == 4)) {
    fCoefVdrift[1][idect] = fVdriftCoef[1];
    fCoefVdriftE[1] [idect] = vdriftCoefE;
    fCoefT0[1][idect] = fT0Coef[1];
    fCoefT0E[1][idect] = t0CoefE;
  }
  
  if (fDebug == 2) {
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
    AliInfo(Form("fVriftCoef[1] (with only the drift region(default)): %f",(Float_t) fVdriftCoef[1]));
    TCanvas *cpentei2 = new TCanvas("cpentei2","cpentei2",50,50,600,800);
    cpentei2->cd();
    pentea->Draw();
    TCanvas *cpentei3 = new TCanvas("cpentei3","cpentei3",50,50,600,800);
    cpentei3->cd();
    pente->Draw();
  }

  if (fDebug != 2) {
    delete pentea;
  }
  if (fDebug != 2) {
    delete pente;
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::FitLagrangePoly(TH1* projPH, Int_t idect)
{
  //
  // Slope methode but with polynomes de Lagrange
  //
  
  // Constants
  const Float_t kDrWidth = AliTRDgeometry::DrThick();
  Int_t binmax   = 0;
  Int_t binmin   = 0;
  Double_t    *x = new Double_t[5];
  Double_t    *y = new Double_t[5];
  x[0]           = 0.0;
  x[1]           = 0.0;
  x[2]           = 0.0;
  x[3]           = 0.0;
  x[4]           = 0.0;
  y[0]           = 0.0;
  y[1]           = 0.0;
  y[2]           = 0.0;
  y[3]           = 0.0;
  y[4]           = 0.0;
  fPhd[0]        = 0.0;
  fPhd[1]        = 0.0;
  fPhd[2]        = 0.0;
  Int_t ju       = 0;
  Double_t vdriftCoefE = 0.0;
  Double_t t0CoefE = 1.0;
  fVdriftCoef[3] = 0.0;
  fT0Coef[3]     = 0.0;
  TLine *line = new TLine();
  TF1 * polynome = 0x0;
  TF1 * polynomea = 0x0;
  TF1 * polynomeb = 0x0;
  Double_t *c = 0x0;
  
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
  if(fDebug == 2) AliInfo(Form("maximum positive bin for the positive slope %d",binmax));


  Double_t minnn = 0.0;
  Double_t maxxx = 0.0;

  //Determination of minnn and maxxx
  //case binmax = nbins -1
  //pol2
  if(binmax == (nbins-1)){
    minnn = pentea->GetBinCenter(binmax-2);
    maxxx = pentea->GetBinCenter(binmax);
    x[0] = pentea->GetBinCenter(binmax-2);
    x[1] = pentea->GetBinCenter(binmax-1);
    x[2] = pentea->GetBinCenter(binmax);
    y[0] = pentea->GetBinContent(binmax-2);
    y[1] = pentea->GetBinContent(binmax-1);
    y[2] = pentea->GetBinContent(binmax);
    //Calcul the polynome de Lagrange
    c = CalculPolynomeLagrange2(x,y);
    AliInfo("At the limit for beginning!");
  }
  //case binmax = nbins-2
  //pol3
  if(binmax == (nbins-2)){
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
    //Calcul the polynome de Lagrange
    c = CalculPolynomeLagrange3(x,y);
  }
  //case binmax <= nbins-3
  //pol4
  if(binmax <= (nbins-3)){
    if((binmax-2) >= 1){
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
    //Calcul the polynome de Lagrange
    c = CalculPolynomeLagrange4(x,y);
    }
    //pol3
    if((binmax-1) == 1){
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
    //Calcul the polynome de Lagrange
    c = CalculPolynomeLagrange3(x,y);
    }
    //pol2
    if(binmax == 1){
    minnn = pentea->GetBinCenter(binmax);
    maxxx = pentea->GetBinCenter(binmax+2);
    x[0] = pentea->GetBinCenter(binmax);
    x[1] = pentea->GetBinCenter(binmax+1);
    x[2] = pentea->GetBinCenter(binmax+2);
    y[0] = pentea->GetBinContent(binmax);
    y[1] = pentea->GetBinContent(binmax+1);
    y[2] = pentea->GetBinContent(binmax+2);
    //Calcul the polynome de Lagrange
    c = CalculPolynomeLagrange2(x,y);
    }
  }
  //pass but should not happen
  if((binmax <= (nbins-3)) && (binmax < 1)){
    put = kFALSE;
  }
     
  if(fDebug == 2) AliInfo(Form("For the beginning region binmax %d",binmax));

  if(put) {
    polynomeb = new TF1("polb","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",minnn,maxxx);
    polynomeb->SetParameters(c[0],c[1],c[2],c[3],c[4]);
    if(fDebug == 2) {
      AliInfo(Form("for the beginning: c[0] %f, c[1] %f, c[2] %f, c[3] %f, c[4] %f",c[0],c[1],c[2],c[3],c[4]));
    }
    
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

  if(fDebug == 2) AliInfo(Form("maximum extrapolated positive bin for the positive slope %f",fPhd[0]));
  
  // Amplification region
  binmax = 0;
  ju     = 0;
  for (Int_t kbin = 1; kbin < projPH->GetNbinsX(); kbin ++) {
    if (((projPH->GetBinContent(kbin+1) - projPH->GetBinContent(kbin)) <= 0.0) && 
         (ju == 0) && 
         (kbin > (fPhd[0]/widbins))) {
      binmax = kbin;
      ju     = 1;
    }
  }
  if(fDebug == 2) AliInfo(Form("For the amplification region binmax %d",binmax));
  
  Double_t minn = 0.0;
  Double_t maxx = 0.0;

  //Determination of minn and maxx
  //case binmax = nbins
  //pol2
  if(binmax == nbins){
    minn = projPH->GetBinCenter(binmax-2);
    maxx = projPH->GetBinCenter(binmax);
    x[0] = projPH->GetBinCenter(binmax-2);
    x[1] = projPH->GetBinCenter(binmax-1);
    x[2] = projPH->GetBinCenter(binmax);
    y[0] = projPH->GetBinContent(binmax-2);
    y[1] = projPH->GetBinContent(binmax-1);
    y[2] = projPH->GetBinContent(binmax);
    //Calcul the polynome de Lagrange
    c = CalculPolynomeLagrange2(x,y);
    AliInfo("At the limit for the drift!");
  }
  //case binmax = nbins-1
  //pol3
  if(binmax == (nbins-1)){
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
    //Calcul the polynome de Lagrange
    c = CalculPolynomeLagrange3(x,y);
  }
  //case binmax <= nbins-2
  //pol4
  if(binmax <= (nbins-2)){
    if((binmax-2) >= 1){
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
    //Calcul the polynome de Lagrange
    c = CalculPolynomeLagrange4(x,y);
    }
    //pol3
    if((binmax-1) == 1){
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
    //Calcul the polynome de Lagrange
    c = CalculPolynomeLagrange3(x,y);
    }
    //pol2
    if(binmax == 1){
    minn = projPH->GetBinCenter(binmax);
    maxx = projPH->GetBinCenter(binmax+2);
    x[0] = projPH->GetBinCenter(binmax);
    x[1] = projPH->GetBinCenter(binmax+1);
    x[2] = projPH->GetBinCenter(binmax+2);
    y[0] = projPH->GetBinContent(binmax);
    y[1] = projPH->GetBinContent(binmax+1);
    y[2] = projPH->GetBinContent(binmax+2);
    //Calcul the polynome de Lagrange
    c = CalculPolynomeLagrange2(x,y);
    }
  }
  //pass but should not happen
  if((binmax <= (nbins-2)) && (binmax < 1)){
    put = kFALSE;
  }
     
  if(fDebug == 2) AliInfo(Form("For the amplification region binmax %d",binmax));

  if(put) {
    polynomea = new TF1("pola","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",minn,maxx);
    polynomea->SetParameters(c[0],c[1],c[2],c[3],c[4]);
    if(fDebug == 2) {
      AliInfo(Form("for the amplification: c[0] %f, c[1] %f, c[2] %f, c[3] %f, c[4] %f",c[0],c[1],c[2],c[3],c[4]));
    }
    
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
  
  if(fDebug == 2) AliInfo(Form("For the amplification region extrapolated binmax %f",fPhd[1]));

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
    AliInfo("Put the binmax from 1 to 2 to enable the fit");
  }
  
  //check
  if((projPH->GetBinContent(binmin)-projPH->GetBinError(binmin)) < (projPH->GetBinContent(binmin+1))) put = kFALSE;
  if((projPH->GetBinContent(binmin)+projPH->GetBinError(binmin)) > (projPH->GetBinContent(binmin-1))) put = kFALSE;

  if(fDebug == 2) {
    AliInfo(Form("binmin %d BinContent %f BinError %f",binmin
                                                      ,projPH->GetBinContent(binmin)
                                                      ,projPH->GetBinError(binmin)));
    AliInfo(Form("binmin-1 %d BinContent %f BinError %f",binmin-1
                                                        ,projPH->GetBinContent(binmin-1)
                                                        ,projPH->GetBinError(binmin-1)));
    AliInfo(Form("binmin+1 %d BinContent %f BinError %f",binmin+1
                                                        ,projPH->GetBinContent(binmin+1)
                                                        ,projPH->GetBinError(binmin+1)));
  }
   
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
    c = CalculPolynomeLagrange4(x,y);
    //richtung +/-
    if((pente->GetBinContent(binmin+2) <= pente->GetBinContent(binmin+1)) &&
       (pente->GetBinContent(binmin-2) <= pente->GetBinContent(binmin-1))) put = kFALSE;
    if(((binmin+3) <= (nbins-1)) &&
       (pente->GetBinContent(binmin+3) <= pente->GetBinContent(binmin+2)) &&
       ((binmin-3) >= TMath::Min(binmax+4, projPH->GetNbinsX())) &&
       (pente->GetBinContent(binmin-3) <= pente->GetBinContent(binmin-2))) put = kFALSE;
    if((pente->GetBinContent(binmin+2) <= pente->GetBinContent(binmin+1)) &&
       (pente->GetBinContent(binmin-2) > pente->GetBinContent(binmin-1))) case1 = kTRUE;
    if((pente->GetBinContent(binmin+2) > pente->GetBinContent(binmin+1)) &&
       (pente->GetBinContent(binmin-2) <= pente->GetBinContent(binmin-1))) case4 = kTRUE;
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
    c = CalculPolynomeLagrange3(x,y);
    //richtung +: nothing
    //richtung -
    if((pente->GetBinContent(binmin-2) <= pente->GetBinContent(binmin-1))) case2 = kTRUE;
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
    c = CalculPolynomeLagrange3(x,y);
    //richtung +
    if((pente->GetBinContent(binmin+2) <= pente->GetBinContent(binmin+1))) case2 = kTRUE;
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
    c = CalculPolynomeLagrange2(x,y);
    //richtung +
    if((pente->GetBinContent(binmin+2) <= pente->GetBinContent(binmin+1))) put = kFALSE;
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
    c = CalculPolynomeLagrange2(x,y);
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
    c = CalculPolynomeLagrange2(x,y);
    AliInfo("At the limit for the drift!");
    //fluctuation too big!
    //richtung +: nothing
    //richtung -
    if((pente->GetBinContent(binmin-2) <= pente->GetBinContent(binmin-1))) put = kFALSE;
  }
  if((binmin == (nbins-1)) && ((binmin-2) < TMath::Min(binmax+4, projPH->GetNbinsX()))) {
    put = kFALSE;
    AliInfo("At the limit for the drift and not usable!");
  }

  //pass
  if((binmin == (nbins-2)) && ((binmin-1) < TMath::Min(binmax+4, projPH->GetNbinsX()))){
    put = kFALSE;
    AliInfo("For the drift...problem!");
  }
 
  //pass but should not happen
  if((binmin <= (nbins-3)) && (binmin < TMath::Min(binmax+4, projPH->GetNbinsX()))){
    put = kFALSE;
    AliInfo("For the drift...problem!");
  }

  if(fDebug == 2) AliInfo(Form("For the drift region binmax %d",binmin));

  if(put) {
    polynome = new TF1("pol","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",min,max);
    polynome->SetParameters(c[0],c[1],c[2],c[3],c[4]);
    if(fDebug == 2) {
      AliInfo(Form("c[0] %f, c[1] %f, c[2] %f, c[3] %f, c[4] %f",c[0],c[1],c[2],c[3],c[4]));
    }
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

  if(fDebug == 2) AliInfo(Form("For the drift region extrapolated binmax %f",fPhd[2]));

  Float_t fPhdt0 = 0.0;
  if(fTakeTheMaxPH) fPhdt0 = fPhd[1];
  else fPhdt0 = fPhd[0];

  if ((fPhd[2] > fPhd[0]) && 
      (fPhd[2] > fPhd[1]) && 
      (fPhd[1] > fPhd[0]) &&
      (put)) {
    fVdriftCoef[3] = (kDrWidth) / (fPhd[2]-fPhd[1]);
    if(fFitPHNDB == 3) fNumberFitSuccess++;
    if (fPhdt0 >= 0.0) {
      fT0Coef[3] = (fPhdt0 - fT0Shift) / widbins;
      if (fT0Coef[3] < -1.0) {
        fT0Coef[3] = fT0Coef[2];
      }
    }
    else {
      fT0Coef[3] = fT0Coef[2];
    }
  }
  else {
    fVdriftCoef[3] = -TMath::Abs(fVdriftCoef[2]);
    fT0Coef[3]     = fT0Coef[2];
  }

  if ((fDebug == 1) || 
      (fDebug == 4)) {
    fCoefVdrift[3][idect] = fVdriftCoef[3];
    fCoefVdriftE[2] [idect] = vdriftCoefE;
    fCoefT0[3][idect] = fT0Coef[3];
    fCoefT0E[2][idect] = t0CoefE;
  }
  
  if (fDebug == 2) {
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
    AliInfo(Form("fVriftCoef[3] (with only the drift region(default)): %f",(Float_t) fVdriftCoef[3]));
    TCanvas *cpentei2 = new TCanvas("cpentei2","cpentei2",50,50,600,800);
    cpentei2->cd();
    pentea->Draw();
    TCanvas *cpentei3 = new TCanvas("cpentei3","cpentei3",50,50,600,800);
    cpentei3->cd();
    pente->Draw();
  }

  if (fDebug != 2) {
    delete pentea;
    delete pente;
    delete polynome;
    delete polynomea;
    delete polynomeb;
  }

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

  fVdriftCoef[0] = 0.0;
  fT0Coef[0]     = 0.0;
  Double_t vdriftCoefE = 0.0;
  Double_t t0CoefE = 0.0;
 
  if (idect%fFitPHPeriode == 0) {

    AliInfo(Form("<AliTRDCalibraFit::FitPH> The detector %d will be fitted",idect));
    fPH->SetParameter(0,(projPH->Integral()-(projPH->GetBinContent(1)*projPH->GetNbinsX())) * 0.00028); // Scaling
    fPH->SetParameter(1,fPhd[0] - 0.1);                                                                 // Start 
    fPH->SetParameter(2,fPhd[1] - fPhd[0]);                                                             // AR
    fPH->SetParameter(3,fPhd[2] - fPhd[1]);                                                             // DR
    fPH->SetParameter(4,0.225);                                                                         // QA/QD
    fPH->SetParameter(5,(Float_t) projPH->GetBinContent(1));
    
    if (fDebug != 2) {
      projPH->Fit(fPH,"0M","",0.0,upedge);
    }

    if (fDebug == 2) {
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
      if(fFitPHNDB == 0) fNumberFitSuccess++;
      fVdriftCoef[0] = kDrWidth / (fPH->GetParameter(3));
      vdriftCoefE = (fPH->GetParError(3)/fPH->GetParameter(3))*fVdriftCoef[0];
      fT0Coef[0]     = fPH->GetParameter(1);
      t0CoefE = fPH->GetParError(1);
    } 
    else {
      fVdriftCoef[0] = -TMath::Abs(fVdriftCoef[2]);
      fT0Coef[0]     = fT0Coef[2];
    }

    if ((fDebug == 1) || 
        (fDebug == 4)) {
      fCoefVdrift[0][idect] = fVdriftCoef[0];
      fCoefVdriftE[0][idect] = vdriftCoefE;
      fCoefT0[0][idect] = fT0Coef[0];
      fCoefT0E[0][idect] = t0CoefE;
    }
    if (fDebug == 2) {
      AliInfo(Form("fVdriftCoef[0]: %f",(Float_t) fVdriftCoef[0]));
    }

  }

  else {

    // Put the default value 
    if ((fDebug <= 1) || 
        (fDebug == 4)) {
      fCoefVdrift[0][idect] = -TMath::Abs(fVdriftCoef[2]);
      fCoefT0[0][idect] = -TMath::Abs(fT0Coef[2]);
    }

  }

  if (fDebug != 2) {
    delete fPH;
  }
  
}

//_____________________________________________________________________________
void AliTRDCalibraFit::FitPRF(TH1 *projPRF, Int_t idect)
{
  //
  // Fit methode for the sigma of the pad response function
  //
  
  fPRFCoef[0] = 0.0;
  Double_t prfCoefE = 0.0;

  if (fDebug != 2) {
    
    projPRF->Fit("gaus","0M","",-fRangeFitPRF,fRangeFitPRF);
    
  }
  
  if (fDebug == 2) {
    TCanvas *cfit = new TCanvas("cfit","cfit",50,50,600,800);
    cfit->cd();
    projPRF->Fit("gaus","M+","",-fRangeFitPRF,fRangeFitPRF);
    projPRF->Draw();

  }

  fPRFCoef[0] = projPRF->GetFunction("gaus")->GetParameter(2);
  prfCoefE = projPRF->GetFunction("gaus")->GetParError(2);
  
  if(fPRFCoef[0] <= 0.0) fPRFCoef[0] = -fPRFCoef[1];
  else {
    if(fFitPRFNDB == 0) fNumberFitSuccess++;
  }

  if ((fDebug == 1) || 
      (fDebug == 4)) {
    fCoefPRF[0][idect] = fPRFCoef[0];
    fCoefPRFE[0][idect] = prfCoefE;
  }
  if (fDebug == 2) {
    AliInfo(Form("fPRFCoef[0]: %f",(Float_t) fPRFCoef[0]));
  }
  
}

//_____________________________________________________________________________
void AliTRDCalibraFit::RmsPRF(TH1 *projPRF, Int_t idect)
{
  //
  // Fit methode for the sigma of the pad response function
  //
  
  fPRFCoef[2] = 0.0;
  Double_t prfCoefE = 0.0;
  
 
  if (fDebug == 2) {
    TCanvas *cfit = new TCanvas("cfit","cfit",50,50,600,800);
    cfit->cd();
    projPRF->Draw();

  }

  fPRFCoef[2] = projPRF->GetRMS();
   
  if(fPRFCoef[2] <= 0.0) fPRFCoef[2] = -fPRFCoef[1];
  else {
    if(fFitPRFNDB == 2) fNumberFitSuccess++;
  }

  if ((fDebug == 1) || 
      (fDebug == 4)) {
    fCoefPRF[2][idect] = fPRFCoef[2];
    fCoefPRFE[1][idect] = prfCoefE;
  }
  if (fDebug == 2) {
    AliInfo(Form("fPRFCoef[2]: %f",(Float_t) fPRFCoef[2]));
  }
  
}

//_____________________________________________________________________________
void AliTRDCalibraFit::FitMean(TH1 *projch, Int_t idect, Double_t nentries)
{
  //
  // Only mean methode for the gain factor
  //
 
  Double_t chargeCoefE1 = 0.0;
  if(nentries > 0) chargeCoefE1 = projch->GetRMS()/TMath::Sqrt(nentries);
 
  if (fDebug == 2) {
    TCanvas *cpmean = new TCanvas("cpmean","cpmean",50,50,600,800);
    cpmean->cd();
    projch->Draw();
  }
  
  if(fFitChargeNDB == 1){
    CalculChargeCoefMean(fCountDet[0],(Int_t) idect,kTRUE);
    fNumberFitSuccess++;
  }
  if ((fDebug == 1) || 
      (fDebug == 4)) {
    fCoefCharge[1][idect]= fChargeCoef[1];
    fCoefChargeE[1][idect]= chargeCoefE1;
  }
}

//_____________________________________________________________________________
void AliTRDCalibraFit::FitMeanW(TH1 *projch, Int_t idect)
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

  //    0 |0.00228515
  //    1 |-0.00231487
  //    2 |0.00044298
  //    3 |-0.00379239
  //    4 |0.00338349

  //A arbitrary error for the moment
  Double_t chargeCoefE4 = 0.0;
  fChargeCoef[4] = 0.0;
  
  //Calcul 
  Double_t sumw = 0.0;
  Double_t sum = 0.0; 
  Int_t sumAll     = (Int_t) projch->GetEntries();
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
  if(sum > 0.0) fChargeCoef[4] = (sumw/sum);

  if (fDebug == 2) {
    AliInfo(Form("fChargeCoef[4] is %f for the dect %d",fChargeCoef[4],idect));
    TCanvas *cpmeanw = new TCanvas("cpmeanw","cpmeanw",50,50,600,800);
    cpmeanw->cd();
    projch->Draw();
  }
  
  if(fFitChargeNDB == 4){
    fNumberFitSuccess++;
    CalculChargeCoefMean(fCountDet[0],(Int_t) idect,kTRUE);
  }
  if ((fDebug == 1) || 
      (fDebug == 4)) {
    fCoefCharge[4][idect]= fChargeCoef[4];
    fCoefChargeE[3][idect]= chargeCoefE4;
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::FitCH(TH1 *projch, Int_t idect)
{
  //
  // Fit methode for the gain factor
  //
 
  fChargeCoef[0] = 0.0;
  Double_t chargeCoefE0 = 0.0;
  Double_t chisqrl = 0.0;
  Double_t chisqrg = 0.0;
  Double_t chisqr = 0.0;
  TF1 *fLandauGaus = new TF1("fLandauGaus",FuncLandauGaus,0,300,5);

  projch->Fit("landau","0",""
             ,(Float_t) fChargeCoef[1]/fBeginFitCharge
             ,projch->GetBinCenter(projch->GetNbinsX()));
  Double_t l3P0         = projch->GetFunction("landau")->GetParameter(0);
  Double_t l3P1         = projch->GetFunction("landau")->GetParameter(1);
  Double_t l3P2         = projch->GetFunction("landau")->GetParameter(2);
  chisqrl = projch->GetFunction("landau")->GetChisquare();
    
  projch->Fit("gaus","0",""
	      ,(Float_t) fChargeCoef[1]/fBeginFitCharge
	      ,projch->GetBinCenter(projch->GetNbinsX()));
  Double_t g3P0         = projch->GetFunction("gaus")->GetParameter(0);
  Double_t g3P2         = projch->GetFunction("gaus")->GetParameter(2);
  chisqrg = projch->GetFunction("gaus")->GetChisquare();
        
  fLandauGaus->SetParameters(l3P0,l3P1,l3P2,g3P0,g3P2);
  if ((fDebug <= 1) || 
      (fDebug >= 3)) {
    projch->Fit("fLandauGaus","0",""
		,(Float_t) fChargeCoef[1]/fBeginFitCharge
		,projch->GetBinCenter(projch->GetNbinsX()));
    chisqr = projch->GetFunction("fLandauGaus")->GetChisquare();
  }
  
  if (fDebug == 2) {
    TCanvas *cp = new TCanvas("cp","cp",50,50,600,800);
    cp->cd();
    projch->Fit("fLandauGaus","+",""
		,(Float_t) fChargeCoef[1]/fBeginFitCharge
		,projch->GetBinCenter(projch->GetNbinsX()));
    chisqr = projch->GetFunction("fLandauGaus")->GetChisquare();
    projch->Draw();
    fLandauGaus->Draw("same");
  }
  
  if ((projch->GetFunction("fLandauGaus")->GetParameter(1) > 0) && 
      (projch->GetFunction("fLandauGaus")->GetParError(1) < 
        (0.05*projch->GetFunction("fLandauGaus")->GetParameter(1))) && 
      (chisqr < chisqrl) && 
      (chisqr < chisqrg)) {
    // Calcul of "real" coef
    if(fFitChargeNDB == 0){
      fNumberFitSuccess++;
      CalculChargeCoefMean(fCountDet[0],(Int_t) idect,kTRUE);
    }
    else {
      CalculChargeCoefMean(fCountDet[0],(Int_t) idect,kFALSE);
    }
    fChargeCoef[0] = projch->GetFunction("fLandauGaus")->GetParameter(1);
    chargeCoefE0 = projch->GetFunction("fLandauGaus")->GetParError(1);
  }
  else {
    // Calcul of "real" coef
    CalculChargeCoefMean(fCountDet[0],(Int_t) idect,kFALSE);
    fChargeCoef[0] = -TMath::Abs(fChargeCoef[3]);
  }
  
  if (fDebug == 2) {
    AliInfo(Form("fChargeCoef[0]: %f",(Float_t) fChargeCoef[0]));
    AliInfo(Form("fChargeCoef[1]: %f",(Float_t) fChargeCoef[1]));
  }
  
  if ((fDebug == 1) || 
      (fDebug == 4)) {
    if (fChargeCoef[0] > 0.0) {
      fCoefCharge[0][idect]= fChargeCoef[0];
      fCoefChargeE[0][idect]= chargeCoefE0;
    }
  }
   
  if (fDebug != 2) {
    delete fLandauGaus;
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::FitBisCH(TH1* projch, Int_t idect)
{
  //
  // Fit methode for the gain factor more time consuming
  //

  //Some parameters to initialise
  Double_t widthLandau, widthGaus, MPV, Integral;
  Double_t chisquarel = 0.0;
  Double_t chisquareg = 0.0;
 
  projch->Fit("landau","0M+",""
	      ,(Float_t) fChargeCoef[1]/6
	      ,projch->GetBinCenter(projch->GetNbinsX()));
  widthLandau  = projch->GetFunction("landau")->GetParameter(2);
  chisquarel = projch->GetFunction("landau")->GetChisquare();
  
  projch->Fit("gaus","0M+",""
	      ,(Float_t) fChargeCoef[1]/6
	      ,projch->GetBinCenter(projch->GetNbinsX()));
  widthGaus    = projch->GetFunction("gaus")->GetParameter(2);
  chisquareg = projch->GetFunction("gaus")->GetChisquare();

  MPV      = (projch->GetFunction("landau")->GetParameter(1))/2;
  Integral = (projch->GetFunction("gaus")->Integral(0.3*fChargeCoef[1],3*fChargeCoef[1])
            + projch->GetFunction("landau")->Integral(0.3*fChargeCoef[1],3*fChargeCoef[1]))/2;

  // Setting fit range and start values
  Double_t fr[2];
  //Double_t sv[4] = { l3P2, fChargeCoef[1], projch->Integral("width"), fG3P2 };
  //Double_t sv[4]   = { fL3P2, fChargeCoef[1], fL3P0, fG3P2 };
  Double_t sv[4]   = { widthLandau, MPV, Integral, widthGaus};
  Double_t pllo[4] = { 0.001, 0.001, projch->Integral()/3, 0.001};
  Double_t plhi[4] = { 300.0, 300.0, 30*projch->Integral(), 300.0};
  Double_t fp[4]   = { 1.0, 1.0, 1.0, 1.0 };
  Double_t fpe[4]  = { 1.0, 1.0, 1.0, 1.0 };
  fr[0]            = 0.3 * fChargeCoef[1];
  fr[1]            = 3.0 * fChargeCoef[1];
  fChargeCoef[2]   = 0.0;
  Double_t chargeCoefE2 = 0.0;

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
    if(fFitChargeNDB == 2){
      fNumberFitSuccess++;
      CalculChargeCoefMean(fCountDet[0],(Int_t) idect,kTRUE);
    }
    else {
      CalculChargeCoefMean(fCountDet[0],(Int_t) idect,kFALSE);
    }
    fChargeCoef[2] = fp[1];
    chargeCoefE2 = fpe[1];
    //chargeCoefE2 = chisqr;
  } 
  else {
    CalculChargeCoefMean(fCountDet[0],(Int_t) idect,kFALSE);
    fChargeCoef[2] = -TMath::Abs(fChargeCoef[3]);
  }
  
  if (fDebug == 2) {
    AliInfo(Form("fChargeCoef[2]: %f",(Float_t) fChargeCoef[2]));
    TCanvas *cpy = new TCanvas("cpy","cpy",50,50,600,800);
    cpy->cd();
    projch->Draw();
    fitsnr->Draw("same");
  }

  if ((fDebug == 1) || 
      (fDebug == 4)) {
    if (fChargeCoef[2] > 0.0) {
      fCoefCharge[2][idect]= fChargeCoef[2];
      fCoefChargeE[2][idect]= chargeCoefE2;
    }
  }

  if (fDebug != 2) {
    delete fitsnr;
  }

}

//_____________________________________________________________________________
Double_t *AliTRDCalibraFit::CalculPolynomeLagrange2(Double_t *x, Double_t *y)
{
  //
  // Calcul the coefficients of the polynome passant par ces trois points de degre 2
  //
 
 Double_t *c = new Double_t[5];
  Double_t x0 = y[0]/((x[0]-x[1])*(x[0]-x[2]));
  Double_t x1 = y[1]/((x[1]-x[0])*(x[1]-x[2]));
  Double_t x2 = y[2]/((x[2]-x[0])*(x[2]-x[1]));

  c[4] = 0.0;
  c[3] = 0.0;
  c[2] = x0+x1+x2;
  c[1] = -(x0*(x[1]+x[2])+x1*(x[0]+x[2])+x2*(x[0]+x[1]));
  c[0] = x0*x[1]*x[2]+x1*x[0]*x[2]+x2*x[0]*x[1];

  return c;
  
}

//_____________________________________________________________________________
Double_t *AliTRDCalibraFit::CalculPolynomeLagrange3(Double_t *x, Double_t *y)
{
  //
  // Calcul the coefficients of the polynome passant par ces quatre points de degre 3
  //

  Double_t *c = new Double_t[5];
  Double_t x0 = y[0]/((x[0]-x[1])*(x[0]-x[2])*(x[0]-x[3]));
  Double_t x1 = y[1]/((x[1]-x[0])*(x[1]-x[2])*(x[1]-x[3]));
  Double_t x2 = y[2]/((x[2]-x[0])*(x[2]-x[1])*(x[2]-x[3]));
  Double_t x3 = y[3]/((x[3]-x[0])*(x[3]-x[1])*(x[3]-x[2]));

  c[4] = 0.0;
  c[3] = x0+x1+x2+x3;
  c[2] = -(x0*(x[1]+x[2]+x[3])
	   +x1*(x[0]+x[2]+x[3])
	   +x2*(x[0]+x[1]+x[3])
	   +x3*(x[0]+x[1]+x[2]));
  c[1] = (x0*(x[1]*x[2]+x[1]*x[3]+x[2]*x[3])
	  +x1*(x[0]*x[2]+x[0]*x[3]+x[2]*x[3])
	  +x2*(x[0]*x[1]+x[0]*x[3]+x[1]*x[3])
	  +x3*(x[0]*x[1]+x[0]*x[2]+x[1]*x[2]));
  
  c[0] = -(x0*x[1]*x[2]*x[3]
	  +x1*x[0]*x[2]*x[3]
	  +x2*x[0]*x[1]*x[3]
	  +x3*x[0]*x[1]*x[2]);  

  return c;

}

//_____________________________________________________________________________
Double_t *AliTRDCalibraFit::CalculPolynomeLagrange4(Double_t *x, Double_t *y)
{
  //
  // Calcul the coefficients of the polynome passant par ces cinqs points de degre 4
  //

  Double_t *c = new Double_t[5];
  Double_t x0 = y[0]/((x[0]-x[1])*(x[0]-x[2])*(x[0]-x[3])*(x[0]-x[4]));
  Double_t x1 = y[1]/((x[1]-x[0])*(x[1]-x[2])*(x[1]-x[3])*(x[1]-x[4]));
  Double_t x2 = y[2]/((x[2]-x[0])*(x[2]-x[1])*(x[2]-x[3])*(x[2]-x[4]));
  Double_t x3 = y[3]/((x[3]-x[0])*(x[3]-x[1])*(x[3]-x[2])*(x[3]-x[4]));
  Double_t x4 = y[4]/((x[4]-x[0])*(x[4]-x[1])*(x[4]-x[2])*(x[4]-x[3]));

  c[4] = x0+x1+x2+x3+x4;
  c[3] = -(x0*(x[1]+x[2]+x[3]+x[4])
	   +x1*(x[0]+x[2]+x[3]+x[4])
	   +x2*(x[0]+x[1]+x[3]+x[4])
	   +x3*(x[0]+x[1]+x[2]+x[4])
	   +x4*(x[0]+x[1]+x[2]+x[3]));
  c[2] = (x0*(x[1]*x[2]+x[1]*x[3]+x[1]*x[4]+x[2]*x[3]+x[2]*x[4]+x[3]*x[4])
	  +x1*(x[0]*x[2]+x[0]*x[3]+x[0]*x[4]+x[2]*x[3]+x[2]*x[4]+x[3]*x[4])
	  +x2*(x[0]*x[1]+x[0]*x[3]+x[0]*x[4]+x[1]*x[3]+x[1]*x[4]+x[3]*x[4])
	  +x3*(x[0]*x[1]+x[0]*x[2]+x[0]*x[4]+x[1]*x[2]+x[1]*x[4]+x[2]*x[4])
	  +x4*(x[0]*x[1]+x[0]*x[2]+x[0]*x[3]+x[1]*x[2]+x[1]*x[3]+x[2]*x[3]));

  c[1] = -(x0*(x[1]*x[2]*x[3]+x[1]*x[2]*x[4]+x[1]*x[3]*x[4]+x[2]*x[3]*x[4])
	  +x1*(x[0]*x[2]*x[3]+x[0]*x[2]*x[4]+x[0]*x[3]*x[4]+x[2]*x[3]*x[4])
	  +x2*(x[0]*x[1]*x[3]+x[0]*x[1]*x[4]+x[0]*x[3]*x[4]+x[1]*x[3]*x[4])
	  +x3*(x[0]*x[1]*x[2]+x[0]*x[1]*x[4]+x[0]*x[2]*x[4]+x[1]*x[2]*x[4])
	  +x4*(x[0]*x[1]*x[2]+x[0]*x[1]*x[3]+x[0]*x[2]*x[3]+x[1]*x[2]*x[3]));

  c[0] = (x0*x[1]*x[2]*x[3]*x[4]
	  +x1*x[0]*x[2]*x[3]*x[4]
	  +x2*x[0]*x[1]*x[3]*x[4]
	  +x3*x[0]*x[1]*x[2]*x[4]
	  +x4*x[0]*x[1]*x[2]*x[3]);

  return c;

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
  for (Int_t k = 0; k < (Int_t) fVectorFitCH->GetEntriesFast(); k++) {
    Int_t    total    = 0;
    Int_t    detector = ((AliTRDFitCHInfo *) fVectorFitCH->At(k))->GetDetector();
    Float_t *coef     = ((AliTRDFitCHInfo *) fVectorFitCH->At(k))->GetCoef();
    //printf("detector %d coef[0] %f\n",detector,coef[0]);
    if (GetChamber(detector) == 2) {
      total = 1728;
    }
    if (GetChamber(detector) != 2) {
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

  if ((fDebug == 3) || 
      (fDebug == 4)) {
    if ((fFitChargeOn) &&
	(fCoefChargeDB[0]->GetEntries()      > 0.0) && 
        (fCoefChargeDB[0]->GetSumOfWeights() > 0.0)) {
      fCoefChargeDB[0]->Scale(fCoefChargeDB[0]->GetEntries() / fCoefChargeDB[0]->GetSumOfWeights());
    }
    if ((fMeanChargeOn)  && 
        (fCoefChargeDB[1]->GetEntries()      > 0.0) && 
        (fCoefChargeDB[1]->GetSumOfWeights() > 0.0)) {
      fCoefChargeDB[1]->Scale(fCoefChargeDB[1]->GetEntries() / fCoefChargeDB[1]->GetSumOfWeights());
    }
    if ((fFitChargeBisOn) && 
        (fCoefChargeDB[2]->GetEntries()      > 0.0) && 
        (fCoefChargeDB[2]->GetSumOfWeights() > 0.0)) {
      fCoefChargeDB[2]->Scale(fCoefChargeDB[2]->GetEntries() / fCoefChargeDB[2]->GetSumOfWeights());
    }
    if ((fFitMeanWOn) && 
        (fCoefChargeDB[3]->GetEntries()      > 0.0) && 
        (fCoefChargeDB[3]->GetSumOfWeights() > 0.0)) {
      fCoefChargeDB[3]->Scale(fCoefChargeDB[3]->GetEntries() / fCoefChargeDB[3]->GetSumOfWeights());
    }
  }
  
}

//_____________________________________________________________________________
TH1I *AliTRDCalibraFit::ReBin(TH1I *hist) const
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

  if (fDebug == 2) {
    TCanvas *crebin = new TCanvas("crebin","",50,50,600,800);
    crebin->cd();
    rehist->Draw();
  }

  return rehist;

}

//_____________________________________________________________________________
TH1F *AliTRDCalibraFit::ReBin(TH1F *hist) const
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

  if (fDebug == 2) {
    TCanvas *crebin = new TCanvas("crebin","",50,50,600,800);
    crebin->cd();
    rehist->Draw();
  }
  
  return rehist;
  
}

//_____________________________________________________________________________
TH1F *AliTRDCalibraFit::CorrectTheError(TGraphErrors *hist)
{
  //
  // In the case of the vectors method the trees contains TGraphErrors for PH and PRF
  // to be able to add them after
  // We convert it to a TH1F to be able to applied the same fit function method
  // After having called this function you can not add the statistics anymore
  //

  TH1F *rehist = 0x0;

  Int_t     nbins   = hist->GetN();
  Double_t *x       = hist->GetX();
  Double_t *entries = hist->GetEX();
  Double_t *mean    = hist->GetY();
  Double_t *square  = hist->GetEY();
  fEntriesCurrent   = 0;

  if (nbins < 2) {
    return rehist; 
  }

  Double_t step     = x[1] - x[0]; 
  Double_t minvalue = x[0] - step/2;
  Double_t maxvalue = x[(nbins-1)] + step/2;

  rehist = new TH1F("projcorrecterror","",nbins,minvalue,maxvalue);

  for (Int_t k = 0; k < nbins; k++) {
    rehist->SetBinContent(k+1,mean[k]);
    if (entries[k] > 0.0) {
      fEntriesCurrent += (Int_t) entries[k];
      Double_t d = TMath::Abs(square[k] - (mean[k]*mean[k]));
      rehist->SetBinError(k+1,TMath::Sqrt(d/entries[k]));
    }
    else {
      rehist->SetBinError(k+1,0.0);
    }
  }

  return rehist;
 
}

//
//____________Some basic geometry function_____________________________________
//

//_____________________________________________________________________________
Int_t AliTRDCalibraFit::GetPlane(Int_t d) const
{
  //
  // Reconstruct the plane number from the detector number
  //

  return ((Int_t) (d % 6));

}

//_____________________________________________________________________________
Int_t AliTRDCalibraFit::GetChamber(Int_t d) const
{
  //
  // Reconstruct the chamber number from the detector number
  //
  Int_t fgkNplan = 6;

  return ((Int_t) (d % 30) / fgkNplan);

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

//_____________________________________________________________________________
void AliTRDCalibraFit::InitTreePRF()
{
  //
  // Init the tree where the coefficients from the fit methods can be stored
  //

  gDirectory = gROOT; 
  fPRFPad    = new Float_t[2304];
  fPRF       = new TTree("PRF","PRF");
  fPRF->Branch("detector",&fPRFDetector,"detector/I");
  fPRF->Branch("width"   ,fPRFPad      ,"width[2304]/F");

  // Set to default value for the plane 0 supposed to be the first one
  for (Int_t k = 0; k < 2304; k++) {
    fPRFPad[k] = 0.515;
  }
  fPRFDetector = -1;

}

//_____________________________________________________________________________
void AliTRDCalibraFit::FillTreePRF(Int_t countdet)
{
  //
  // Fill the tree with the sigma of the pad response function for the detector countdet
  //
  
  Int_t numberofgroup = 0; 
  fPRFDetector = countdet;
  fPRF->Fill();

  if (GetChamber((Int_t)(countdet+1)) == 2) {
    numberofgroup = 1728;
  }
  else {
    numberofgroup = 2304;
  }

  // Reset to default value for the next
  for (Int_t k = 0; k < numberofgroup; k++) {
    if (GetPlane((Int_t) (countdet+1)) == 0) {
      fPRFPad[k] = 0.515;
    }
    if (GetPlane((Int_t) (countdet+1)) == 1) {
      fPRFPad[k] = 0.502;
    }
    if (GetPlane((Int_t) (countdet+1)) == 2) {
      fPRFPad[k] = 0.491;
    }
    if (GetPlane((Int_t) (countdet+1)) == 3) {
      fPRFPad[k] = 0.481;
    }
    if (GetPlane((Int_t) (countdet+1)) == 4) {
      fPRFPad[k] = 0.471;
    }
    if (GetPlane((Int_t) (countdet+1)) == 5) {
      fPRFPad[k] = 0.463;
    }
  }

  fPRFDetector = -1;

}

//_____________________________________________________________________________
void AliTRDCalibraFit::ConvertVectorFitCHTree()
{
  //
  // Convert the vector stuff to a tree of 1D histos if the user
  // want to write it after the fill functions
  //

  Int_t detector      = -1;
  Int_t numberofgroup =  1;
  Float_t gainPad[2304];

  fGain = new TTree("Gain","Gain");
  fGain->Branch("detector",&detector,"detector/I");
  fGain->Branch("gainPad" ,gainPad  ,"gainPad[2304]/F");

  Int_t loop = (Int_t) fVectorFitCH->GetEntriesFast();
  for (Int_t k = 0; k < loop; k++) {
    detector = ((AliTRDFitCHInfo *) fVectorFitCH->At(k))->GetDetector();
    if (GetChamber((Int_t) ((AliTRDFitCHInfo *) fVectorFitCH->At(k))->GetDetector()) == 2) {
      numberofgroup = 1728;
    }
    else {
      numberofgroup = 2304;
    }
    for (Int_t i = 0; i < numberofgroup; i++) {
      if (((AliTRDFitCHInfo *) fVectorFitCH->At(k))->GetCoef()[i] >= 0) {
        gainPad[i] = ((AliTRDFitCHInfo *) fVectorFitCH->At(k))->GetCoef()[i] * fScaleFitFactor;
      }
      else {
	gainPad[i] = (Float_t) ((AliTRDFitCHInfo *) fVectorFitCH->At(k))->GetCoef()[i];
      }
    }
    //Finish the vector
    if(numberofgroup < 2304){
      for(Int_t i = numberofgroup; i < 2304; i++){
	gainPad[i] = -100.0;
      }
    }
    fGain->Fill();
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFit::FillTreeVdrift(Int_t countdet)
{
  //
  // Fill the tree with the drift velocities for the detector countdet
  //

  Int_t numberofgroup = 0;
  fVdriftDetector = countdet;
 
  fVdrift->Fill();
  if (GetChamber((Int_t)(countdet+1)) == 2) {
    numberofgroup = 1728;
  }
  else {
    numberofgroup = 2304;
  }
  // Reset to default value the gain coef
  for (Int_t k = 0; k < numberofgroup; k++) {
    fVdriftPad[k] = -1.5;
  }
  fVdriftDetector = -1;

}

//_____________________________________________________________________________
void AliTRDCalibraFit::InitTreePH()
{
  //
  // Init the tree where the coefficients from the fit methods can be stored
  //
  
  gDirectory = gROOT;
  fVdriftPad = new Float_t[2304];
  fVdrift    = new TTree("Vdrift","Vdrift");
  fVdrift->Branch("detector",&fVdriftDetector,"detector/I");
  fVdrift->Branch("vdrift"  ,fVdriftPad      ,"vdrift[2304]/F");
  // Set to default value for the plane 0 supposed to be the first one
  for (Int_t k = 0; k < 2304; k++) {
    fVdriftPad[k] = -1.5;
  }
  fVdriftDetector = -1;

}

//_____________________________________________________________________________
void AliTRDCalibraFit::FillTreeT0(Int_t countdet)
{
  //
  // Fill the tree with the t0 value for the detector countdet
  //

  Int_t numberofgroup = 0;

  fT0Detector = countdet;
 
  fT0->Fill();
  if (GetChamber((Int_t) (countdet+1)) == 2) {
    numberofgroup = 1728;
  }
  else {
    numberofgroup = 2304;
  }
  // Reset to default value 
  for (Int_t k = 0; k < numberofgroup; k++) {
    fT0Pad[k] = 0.0;
  }
  fT0Detector = -1;

}

//_____________________________________________________________________________
void AliTRDCalibraFit::InitTreeT0()
{
  //
  // Init the tree where the coefficients from the fit methods can be stored
  //
  
  gDirectory = gROOT;
  fT0Pad = new Float_t[2304];
  fT0 = new TTree("T0","T0");
  fT0->Branch("detector",&fT0Detector,"detector/I");
  fT0->Branch("t0",fT0Pad,"t0[2304]/F");
  //Set to default value for the plane 0 supposed to be the first one
  for(Int_t k = 0; k < 2304; k++){
    fT0Pad[k] = 0.0;
  }
  fT0Detector = -1;

}

//
//____________Private Functions________________________________________________
//

//_____________________________________________________________________________
Double_t AliTRDCalibraFit::PH(Double_t *x, Double_t *par) 
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
Double_t AliTRDCalibraFit::AsymmGauss(Double_t *x, Double_t *par) 
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
                               * (1.0 - TMath::Erf((par3save * sigma2 - dx) / (sqrt2 * par2save)));
  Double_t  exp2    = par5save * TMath::Exp(-par5save * (dx - 0.5 * par5save * sigma2))
                               * (1.0 - TMath::Erf((par5save * sigma2 - dx) / (sqrt2 * par2save)));

  //return par[0]*(exp1+par[4]*exp2);
  return par[0] * (exp1 + 1.00124 * exp2);

}

//_____________________________________________________________________________
Double_t AliTRDCalibraFit::FuncLandauGaus(Double_t *x, Double_t *par)
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
Double_t AliTRDCalibraFit::LanGauFun(Double_t *x, Double_t *par) 
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
TF1 *AliTRDCalibraFit::LanGauFit(TH1 *his, Double_t *fitrange, Double_t *startvalues
                                      , Double_t *parlimitslo, Double_t *parlimitshi
                                      , Double_t *fitparams, Double_t *fiterrors
                                      , Double_t *chiSqr, Int_t *ndf)
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
Int_t AliTRDCalibraFit::LanGauPro(Double_t *params, Double_t &maxx, Double_t &fwhm) 
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
Double_t AliTRDCalibraFit::GausConstant(Double_t *x, Double_t *par)
{
  //
  // Gaus with identical mean
  //

  Double_t gauss = par[0] * TMath::Gaus(x[0],0.0,par[1])+par[2];
 
  return gauss;

}
