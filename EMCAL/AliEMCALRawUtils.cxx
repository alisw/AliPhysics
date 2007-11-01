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
/* History of cvs commits:
 *
 * $Log$
 * Revision 1.5  2007/11/01 01:20:33  mvl
 * Further improvement of peak finding; more robust fit
 *
 * Revision 1.4  2007/10/31 17:15:24  mvl
 * Fixed bug in raw data unpacking; Added pedestal to signal fit; Added logic to deal with high/low gain
 *
 * Revision 1.3  2007/09/27 08:36:46  mvl
 * More robust setting of fit range in FitRawSignal (P. Hristov)
 *
 * Revision 1.2  2007/09/03 20:55:35  jklay
 * EMCAL e-by-e reconstruction methods from Cvetan
 *
 * Revision 1.1  2007/03/17 19:56:38  mvl
 * Moved signal shape routines from AliEMCAL to separate class AliEMCALRawUtils to streamline raw data reconstruction code.
 * */

//*-- Author: Marco van Leeuwen (LBL)
#include "AliEMCALRawUtils.h"

#include "TF1.h"
#include "TGraph.h"
#include "TSystem.h"

#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliCaloAltroMapping.h"
#include "AliAltroBuffer.h"
#include "AliRawReader.h"
#include "AliCaloRawStream.h"
#include "AliDAQ.h"

#include "AliEMCALLoader.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALDigit.h"


ClassImp(AliEMCALRawUtils)

// Signal shape parameters
Int_t    AliEMCALRawUtils::fgOrder         = 2 ;      // Order of gamma function 
Double_t AliEMCALRawUtils::fgTimeBinWidth  = 100E-9 ; // each sample is 100 ns
Double_t AliEMCALRawUtils::fgTau         = 235E-9 ;   // 235 ns (from CERN testbeam; not very accurate)
Double_t AliEMCALRawUtils::fgTimeTrigger = 1.5E-6 ;   // 15 time bins ~ 1.5 musec

// some digitization constants
Int_t    AliEMCALRawUtils::fgThreshold = 1;
Int_t    AliEMCALRawUtils::fgDDLPerSuperModule = 2;  // 2 ddls per SuperModule

AliEMCALRawUtils::AliEMCALRawUtils(): fHighLowGainFactor(0.) {
  fHighLowGainFactor = 16. ;          // adjusted for a low gain range of 82 GeV (10 bits) 
}
//____________________________________________________________________________
AliEMCALRawUtils::~AliEMCALRawUtils() {
}
//____________________________________________________________________________
void AliEMCALRawUtils::Digits2Raw()
{
  // convert digits of the current event to raw data
  
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliEMCALLoader *loader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));

  // get the digits
  loader->LoadDigits("EMCAL");
  loader->GetEvent();
  TClonesArray* digits = loader->Digits() ;
  
  if (!digits) {
    Warning("Digits2Raw", "no digits found !");
    return;
  }
    
  // get the geometry
  AliEMCALGeometry* geom = AliEMCALGeometry::GetInstance();
  if (!geom) {
    AliError(Form("No geometry found !"));
    return;
  }
  
  static const Int_t nDDL = 12*2; // 12 SM hardcoded for now. Buffers allocated dynamically, when needed, so just need an upper limit here
  AliAltroBuffer* buffers[nDDL];
  for (Int_t i=0; i < nDDL; i++)
    buffers[i] = 0;

  Int_t adcValuesLow[fgkTimeBins];
  Int_t adcValuesHigh[fgkTimeBins];

  //Load Mapping RCU files once
  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/EMCAL/mapping/RCU";
  TString path0 = path+"0.data";//This file will change in future
  TString path1 = path+"1.data";//This file will change in future
  AliAltroMapping * mapping[2] ; // For the moment only 2
  mapping[0] = new AliCaloAltroMapping(path0.Data());
  mapping[1] = new AliCaloAltroMapping(path1.Data());

  // loop over digits (assume ordered digits)
  for (Int_t iDigit = 0; iDigit < digits->GetEntries(); iDigit++) {
    AliEMCALDigit* digit = dynamic_cast<AliEMCALDigit *>(digits->At(iDigit)) ;
    if (digit->GetAmp() < fgThreshold) 
      continue;

    //get cell indices
    Int_t nSM = 0;
    Int_t nIphi = 0;
    Int_t nIeta = 0;
    Int_t iphi = 0;
    Int_t ieta = 0;
    Int_t nModule = 0;
    geom->GetCellIndex(digit->GetId(), nSM, nModule, nIphi, nIeta);
    geom->GetCellPhiEtaIndexInSModule(nSM, nModule, nIphi, nIeta,iphi, ieta) ;
    
    //Check which is the RCU of the cell.
    Int_t iRCU = -111;
    //RCU0
    if (0<=iphi&&iphi<8) iRCU=0; // first cable row
    else if (8<=iphi&&iphi<16 && 0<=ieta&&ieta<24) iRCU=0; // first half; 
    //second cable row
    //RCU1
    else if(8<=iphi&&iphi<16 && 24<=ieta&&ieta<48) iRCU=1; // second half; 
    //second cable row
    else if(16<=iphi&&iphi<24) iRCU=1; // third cable row
    
    //Which DDL?
    Int_t iDDL = fgDDLPerSuperModule* nSM + iRCU;
    if (iDDL >= nDDL)
      Fatal("Digits2Raw()","Non-existent DDL board number: %d", iDDL);

    if (buffers[iDDL] == 0) {      
      // open new file and write dummy header
      TString fileName = AliDAQ::DdlFileName("EMCAL",iDDL);
      buffers[iDDL] = new AliAltroBuffer(fileName.Data(),mapping[iRCU]);
      buffers[iDDL]->WriteDataHeader(kTRUE, kFALSE);  //Dummy;
    }
    
    // out of time range signal (?)
    if (digit->GetTimeR() > GetRawFormatTimeMax() ) {
      AliInfo("Signal is out of time range.\n");
      buffers[iDDL]->FillBuffer((Int_t)digit->GetAmp());
      buffers[iDDL]->FillBuffer(GetRawFormatTimeBins() );  // time bin
      buffers[iDDL]->FillBuffer(3);          // bunch length      
      buffers[iDDL]->WriteTrailer(3, ieta, iphi, nSM);  // trailer
      // calculate the time response function
    } else {
      Bool_t lowgain = RawSampledResponse(digit->GetTimeR(), digit->GetAmp(), adcValuesHigh, adcValuesLow) ; 
      if (lowgain) 
	buffers[iDDL]->WriteChannel(ieta, iphi, 0, GetRawFormatTimeBins(), adcValuesLow, fgThreshold);
      else 
	buffers[iDDL]->WriteChannel(ieta,iphi, 1, GetRawFormatTimeBins(), adcValuesHigh, fgThreshold);
    }
  }
  
  // write headers and close files
  for (Int_t i=0; i < nDDL; i++) {
    if (buffers[i]) {
      buffers[i]->Flush();
      buffers[i]->WriteDataHeader(kFALSE, kFALSE);
      delete buffers[i];
    }
  }
  mapping[0]->Delete();
  mapping[1]->Delete();
  loader->UnloadDigits();
}

//____________________________________________________________________________
void AliEMCALRawUtils::Raw2Digits(AliRawReader* reader,TClonesArray *digitsArr)
{
  // convert raw data of the current event to digits
  AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance();
  if (!geom) {
    AliError(Form("No geometry found !"));
    return;
  }

  digitsArr->Clear(); 

  if (!digitsArr) {
    Error("Raw2Digits", "no digits found !");
    return;
  }
  if (!reader) {
    Error("Raw2Digits", "no raw reader found !");
    return;
  }

  // Use AliAltroRawStream to read the ALTRO format.  No need to
  // reinvent the wheel :-) 
  AliCaloRawStream in(reader,"EMCAL");
  // Select EMCAL DDL's;
  reader->Select("EMCAL");
  //in.SetOldRCUFormat(kTRUE); // Needed for testbeam data
  
  cout << "Stream set up" << endl;

  // reading is from previously existing AliEMCALGetter.cxx
  // ReadRaw method
  TF1 * signalF = new TF1("signal", RawResponseFunction, 0, GetRawFormatTimeMax(), 4);
  
  Int_t id =  -1;
  Float_t time = 0. ; 
  Float_t amp = 0. ; 

  TGraph * gSig = new TGraph(GetRawFormatTimeBins()) ; 

  Int_t readOk = 1;
  Int_t lowGain = 0;

  while (readOk && in.GetModule() < 0) 
    readOk = in.Next();  // Go to first digit

  while (readOk) { 
    id =  geom->GetAbsCellIdFromCellIndexes(in.GetModule(), in.GetRow(), in.GetColumn()) ;
    lowGain = in.IsLowGain();
    Int_t maxTime = in.GetTime();  // timebins come in reverse order
    if (maxTime < 0 || maxTime >= GetRawFormatTimeBins()) {
      AliWarning(Form("Invalid time bin %d",maxTime));
      maxTime = GetRawFormatTimeBins();
    }
    gSig->Set(maxTime+1);
    // There is some kind of zero-suppression in the raw data, 
    // so set up the TGraph in advance
    for (Int_t i=0; i < maxTime; i++) {
      gSig->SetPoint(i, i * GetRawFormatTimeBinWidth(), 0);
    }

    Int_t iTime = 0;
    do {
      if (in.GetTime() >= gSig->GetN()) {
	  AliWarning("Too many time bins");
	  gSig->Set(in.GetTime());
      }
      gSig->SetPoint(in.GetTime(), 
		   in.GetTime() * GetRawFormatTimeBinWidth(), 
		   in.GetSignal()) ;
      if (in.GetTime() > maxTime)
        maxTime = in.GetTime();
      iTime++;
    } while ((readOk = in.Next()) && !in.IsNewHWAddress());
    signalF->SetRange(0,(Float_t)maxTime*GetRawFormatTimeBinWidth());

    FitRaw(gSig, signalF, amp, time) ; 
    
    if (amp > 0) {
      AliDebug(2,Form("id %d lowGain %d amp %g", id, lowGain, amp));
      AddDigit(digitsArr, id, lowGain, (Int_t)amp, time);
    }
	
    // Reset graph
    for (Int_t index = 0; index < gSig->GetN(); index++) {
      gSig->SetPoint(index, index * GetRawFormatTimeBinWidth(), 0) ;  
    } 
  }; // EMCAL entries loop
  
  delete signalF ; 
  delete gSig;
  
  return ; 
}

//____________________________________________________________________________ 
void AliEMCALRawUtils::AddDigit(TClonesArray *digitsArr, Int_t id, Int_t lowGain, Int_t amp, Float_t time) {
  //
  // Add a new digit. 
  // This routine checks whether a digit exists already for this tower 
  // and then decides whether to use the high or low gain info
  //
  // Called by Raw2Digits
  
  AliEMCALDigit *digit = 0, *tmpdigit = 0;
  
  TIter nextdigit(digitsArr);
  while (digit == 0 && (tmpdigit = (AliEMCALDigit*) nextdigit())) {
    if (tmpdigit->GetId() == id)
      digit = tmpdigit;
  }

  if (!digit) { // no digit existed for this tower; create one
    if (lowGain) 
      amp = Int_t(fHighLowGainFactor * amp); 
    Int_t idigit = digitsArr->GetEntries();
    new((*digitsArr)[idigit]) AliEMCALDigit( -1, -1, id, amp, time, idigit) ;	
  }
  else { // a digit already exists, check range 
         // (use high gain if signal < 800, otherwise low gain)
    if (lowGain) { // new digit is low gain
      if (digit->GetAmp() > 800) {  // use if stored digit is out of range
	digit->SetAmp(Int_t(fHighLowGainFactor * amp));
	digit->SetTime(time);
      }
    }
    else if (amp < 800) { // new digit is high gain; use if not out of range
      digit->SetAmp(amp);
      digit->SetTime(time);
    }
  }
}

//____________________________________________________________________________ 
void AliEMCALRawUtils::FitRaw(TGraph * gSig, TF1* signalF, Float_t & amp, Float_t & time)
{
  // Fits the raw signal time distribution; from AliEMCALGetter 

  const Int_t kNoiseThreshold = 5;
  const Int_t kNPedSamples = 10;
  amp = time = 0. ; 
  Double_t ped = 0;
  Int_t nPed = 0;

  for (Int_t index = 0; index < kNPedSamples; index++) {
    Double_t ttime, signal;
    gSig->GetPoint(index, ttime, signal) ; 
    if (signal > 0) {
      ped += signal;
      nPed++;
    }
  }

  if (nPed > 0)
    ped /= nPed;
  else {
    AliWarning("Could determine pedestal");	  
    ped = 10; // put some small value as first guess
  }

  Int_t max_found = 0;
  Int_t i_max = 0;
  Float_t max = -1;
  Float_t tmax = 0;
  Float_t max_fit = gSig->GetN()*GetRawFormatTimeBinWidth();
  Float_t min_after_sig = 9999;
  Int_t imin_after_sig = gSig->GetN();
  Float_t tmin_after_sig = gSig->GetN()*GetRawFormatTimeBinWidth();
  Int_t n_ped_after_sig = 0;

  for (Int_t i=kNPedSamples; i < gSig->GetN(); i++) {
    Double_t ttime, signal;
    gSig->GetPoint(i, ttime, signal) ; 
    if (!max_found && signal > max) {
      i_max = i;
      tmax = ttime;
      max = signal;
    }
    else if ( max > ped + kNoiseThreshold ) {
      max_found = 1;
      min_after_sig = signal;
      imin_after_sig = i;
      tmin_after_sig = ttime;
    }
    if (max_found) {
      if ( signal < min_after_sig) {
        min_after_sig = signal;
	imin_after_sig = i;
        tmin_after_sig = ttime;
      }
      if (i > tmin_after_sig + 5) {  // Two close peaks; end fit at minimum
        max_fit = tmin_after_sig;
        break;
      }
      if ( signal < ped + kNoiseThreshold)
        n_ped_after_sig++;
      if (n_ped_after_sig >= 5) {  // include 5 pedestal bins after peak
        max_fit = ttime;
        break;
      }
    }
  }

  if ( max - ped > kNoiseThreshold ) { // else its noise 
    AliDebug(2,Form("Fitting max %d ped %d", max, ped));
    signalF->SetParameter(0, ped) ; 
    signalF->SetParameter(1, max - ped) ; 
    signalF->SetParameter(2, tmax) ; 
    signalF->SetParLimits(2, 0, max_fit) ; 
    gSig->Fit(signalF, "QRON", "", 0., max_fit); //, "QRON") ; 
    amp = signalF->GetParameter(1); 
    time = signalF->GetParameter(2) - fgTimeTrigger;
  }
  return;
}
//__________________________________________________________________
Double_t AliEMCALRawUtils::RawResponseFunction(Double_t *x, Double_t *par)
{
  // Shape of the electronics raw reponse:
  // It is a semi-gaussian, 2nd order Gamma function of the general form
  //
  // t' = (t - t0 + tau) / tau
  // F = A * t**N * exp( N * ( 1 - t) )   for t >= 0
  // F = 0                                for t < 0 
  //
  // parameters:
  // ped: par[0]
  // A:   par[1]   // Amplitude = peak value
  // t0:  par[2]
  // tau: fgTau
  // N:   fgOrder
  //
  Double_t signal ;
  Double_t xx = ( x[0] - par[2] + fgTau ) / fgTau ; 

  if (xx <= 0) 
    signal = par[0] ;  
  else {  
    signal = par[0] + par[1] * TMath::Power(xx , fgOrder) * TMath::Exp(fgOrder * (1 - xx )) ; 

  }
  return signal ;  
}

//__________________________________________________________________
Bool_t AliEMCALRawUtils::RawSampledResponse(
const Double_t dtime, const Double_t damp, Int_t * adcH, Int_t * adcL) const 
{
  // for a start time dtime and an amplitude damp given by digit, 
  // calculates the raw sampled response AliEMCAL::RawResponseFunction

  const Int_t kRawSignalOverflow = 0x3FF ; 
  const Int_t pedVal = 32;
  Bool_t lowGain = kFALSE ; 

  TF1 signalF("signal", RawResponseFunction, 0, GetRawFormatTimeMax(), 4);
  signalF.SetParameter(0, pedVal) ; 
  signalF.SetParameter(1, damp) ; 
  signalF.SetParameter(2, dtime + fgTimeTrigger) ; 

  for (Int_t iTime = 0; iTime < GetRawFormatTimeBins(); iTime++) {
    Double_t time = iTime * GetRawFormatTimeBinWidth() ;
    Double_t signal = signalF.Eval(time) ;     
    adcH[iTime] =  static_cast<Int_t>(signal + 0.5) ;
    if ( adcH[iTime] > kRawSignalOverflow ){  // larger than 10 bits 
      adcH[iTime] = kRawSignalOverflow ;
      lowGain = kTRUE ; 
    }

    signal /= fHighLowGainFactor;

    adcL[iTime] =  static_cast<Int_t>(signal + 0.5) ;
    if ( adcL[iTime] > kRawSignalOverflow)  // larger than 10 bits 
      adcL[iTime] = kRawSignalOverflow ;
  }
  return lowGain ; 
}
