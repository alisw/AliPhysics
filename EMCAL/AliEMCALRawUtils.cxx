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
 * Revision 1.10  2007/12/06 13:58:11  hristov
 * Additional pritection. Do not delete the mapping, it is owned by another class
 *
 * Revision 1.9  2007/12/06 02:19:51  jklay
 * incorporated fitting procedure from testbeam analysis into AliRoot
 *
 * Revision 1.8  2007/12/05 02:30:51  jklay
 * modification to read Altro mappings into AliEMCALRecParam and pass to AliEMCALRawUtils from AliEMCALReconstructor; add option to AliEMCALRawUtils to set old RCU format (for testbeam) or not
 *
 * Revision 1.7  2007/11/14 15:51:46  gustavo
 * Take out few unnecessary prints
 *
 * Revision 1.6  2007/11/01 01:23:51  mvl
 * Removed call to SetOldRCUFormat, which is only needed for testbeam data
 *
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
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliCaloAltroMapping.h"
#include "AliAltroBuffer.h"
#include "AliRawReader.h"
#include "AliCaloRawStream.h"
#include "AliDAQ.h"
  
#include "AliEMCALRecParam.h"
#include "AliEMCALLoader.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALDigit.h"
#include "AliEMCAL.h"
  
ClassImp(AliEMCALRawUtils)
  
// Signal shape parameters
 Double_t AliEMCALRawUtils::fgTimeBinWidth  = 100E-9 ; // each sample is 100 ns
Double_t AliEMCALRawUtils::fgTimeTrigger = 1.5E-6 ;   // 15 time bins ~ 1.5 musec

// some digitization constants
Int_t    AliEMCALRawUtils::fgThreshold = 1;
Int_t    AliEMCALRawUtils::fgDDLPerSuperModule = 2;  // 2 ddls per SuperModule
Int_t    AliEMCALRawUtils::fgPedestalValue = 32;      // pedestal value for digits2raw
Double_t AliEMCALRawUtils::fgFEENoise = 3.;            // 3 ADC channels of noise (sampled)

AliEMCALRawUtils::AliEMCALRawUtils()
  : fHighLowGainFactor(0.), fOrder(0), fTau(0.), fNoiseThreshold(0),
    fNPedSamples(0), fGeom(0), fOption("")
{

  //These are default parameters.  
  //Can be re-set from without with setter functions
  fHighLowGainFactor = 16. ;          // adjusted for a low gain range of 82 GeV (10 bits) 
  fOrder = 2;                         // order of gamma fn
  fTau = 2.35;                        // in units of timebin, from CERN 2007 testbeam
  fNoiseThreshold = 3;
  fNPedSamples = 5;

  //Get Mapping RCU files from the AliEMCALRecParam                                 
  const TObjArray* maps = AliEMCALRecParam::GetMappings();
  if(!maps) AliFatal("Cannot retrieve ALTRO mappings!!");

  for(Int_t i = 0; i < 4; i++) {
    fMapping[i] = (AliAltroMapping*)maps->At(i);
  }

  //To make sure we match with the geometry in a simulation file,
  //let's try to get it first.  If not, take the default geometry
  AliRunLoader *rl = AliRunLoader::Instance();
  if(!rl) AliError("Cannot find RunLoader!");
  if (rl->GetAliRun() && rl->GetAliRun()->GetDetector("EMCAL")) {
    fGeom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
  } else {
    AliInfo(Form("Using default geometry in raw reco"));
    fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
  }

  if(!fGeom) AliFatal(Form("Could not get geometry!"));

}

//____________________________________________________________________________
AliEMCALRawUtils::AliEMCALRawUtils(AliEMCALGeometry *pGeometry)
  : fHighLowGainFactor(0.), fOrder(0), fTau(0.), fNoiseThreshold(0),
    fNPedSamples(0), fGeom(pGeometry), fOption("")
{
  //
  // Initialize with the given geometry - constructor required by HLT
  // HLT does not use/support AliRunLoader(s) instances
  // This is a minimum intervention solution
  // Comment by MPloskon@lbl.gov
  //

  //These are default parameters. 
  //Can be re-set from without with setter functions 
  fHighLowGainFactor = 16. ;          // adjusted for a low gain range of 82 GeV (10 bits)
  fOrder = 2;                         // order of gamma fn
  fTau = 2.35;                        // in units of timebin, from CERN 2007 testbeam
  fNoiseThreshold = 3;
  fNPedSamples = 5;

  //Get Mapping RCU files from the AliEMCALRecParam
  const TObjArray* maps = AliEMCALRecParam::GetMappings();
  if(!maps) AliFatal("Cannot retrieve ALTRO mappings!!");

  for(Int_t i = 0; i < 4; i++) {
    fMapping[i] = (AliAltroMapping*)maps->At(i);
  }

  if(!fGeom) AliFatal(Form("Could not get geometry!"));

}

//____________________________________________________________________________
AliEMCALRawUtils::AliEMCALRawUtils(const AliEMCALRawUtils& rawU)
  : TObject(),
    fHighLowGainFactor(rawU.fHighLowGainFactor), 
    fOrder(rawU.fOrder),
    fTau(rawU.fTau),
    fNoiseThreshold(rawU.fNoiseThreshold),
    fNPedSamples(rawU.fNPedSamples),
    fGeom(rawU.fGeom), 
    fOption(rawU.fOption)
{
  //copy ctor
  fMapping[0] = rawU.fMapping[0];
  fMapping[1] = rawU.fMapping[1];
  fMapping[2] = rawU.fMapping[2];
  fMapping[3] = rawU.fMapping[3];
}

//____________________________________________________________________________
AliEMCALRawUtils& AliEMCALRawUtils::operator =(const AliEMCALRawUtils &rawU)
{
  //assignment operator

  if(this != &rawU) {
    fHighLowGainFactor = rawU.fHighLowGainFactor;
    fOrder = rawU.fOrder;
    fTau = rawU.fTau;
    fNoiseThreshold = rawU.fNoiseThreshold;
    fNPedSamples = rawU.fNPedSamples;
    fGeom = rawU.fGeom;
    fOption = rawU.fOption;
    fMapping[0] = rawU.fMapping[0];
    fMapping[1] = rawU.fMapping[1];
    fMapping[2] = rawU.fMapping[2];
    fMapping[3] = rawU.fMapping[3];
  }

  return *this;

}

//____________________________________________________________________________
AliEMCALRawUtils::~AliEMCALRawUtils() {

}

//____________________________________________________________________________
void AliEMCALRawUtils::Digits2Raw()
{
  // convert digits of the current event to raw data
  
  AliRunLoader *rl = AliRunLoader::Instance();
  AliEMCALLoader *loader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));

  // get the digits
  loader->LoadDigits("EMCAL");
  loader->GetEvent();
  TClonesArray* digits = loader->Digits() ;
  
  if (!digits) {
    Warning("Digits2Raw", "no digits found !");
    return;
  }

  static const Int_t nDDL = 12*2; // 12 SM hardcoded for now. Buffers allocated dynamically, when needed, so just need an upper limit here
  AliAltroBuffer* buffers[nDDL];
  for (Int_t i=0; i < nDDL; i++)
    buffers[i] = 0;

  Int_t adcValuesLow[fgkTimeBins];
  Int_t adcValuesHigh[fgkTimeBins];

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
    fGeom->GetCellIndex(digit->GetId(), nSM, nModule, nIphi, nIeta);
    fGeom->GetCellPhiEtaIndexInSModule(nSM, nModule, nIphi, nIeta,iphi, ieta) ;
    
    //Check which is the RCU, 0 or 1, of the cell.
    Int_t iRCU = -111;
    //RCU0
    if (0<=iphi&&iphi<8) iRCU=0; // first cable row
    else if (8<=iphi&&iphi<16 && 0<=ieta&&ieta<24) iRCU=0; // first half; 
    //second cable row
    //RCU1
    else if(8<=iphi&&iphi<16 && 24<=ieta&&ieta<48) iRCU=1; // second half; 
    //second cable row
    else if(16<=iphi&&iphi<24) iRCU=1; // third cable row

    if (nSM%2==1) iRCU = 1 - iRCU; // swap for odd=C side, to allow us to cable both sides the same

    if (iRCU<0) 
      Fatal("Digits2Raw()","Non-existent RCU number: %d", iRCU);
    
    //Which DDL?
    Int_t iDDL = fgDDLPerSuperModule* nSM + iRCU;
    if (iDDL >= nDDL)
      Fatal("Digits2Raw()","Non-existent DDL board number: %d", iDDL);

    if (buffers[iDDL] == 0) {      
      // open new file and write dummy header
      TString fileName = AliDAQ::DdlFileName("EMCAL",iDDL);
      //Select mapping file RCU0A, RCU0C, RCU1A, RCU1C
      Int_t iRCUside=iRCU+(nSM%2)*2;
      //iRCU=0 and even (0) SM -> RCU0A.data   0
      //iRCU=1 and even (0) SM -> RCU1A.data   1
      //iRCU=0 and odd  (1) SM -> RCU0C.data   2
      //iRCU=1 and odd  (1) SM -> RCU1C.data   3
      //cout<<" nSM "<<nSM<<"; iRCU "<<iRCU<<"; iRCUside "<<iRCUside<<endl;
      buffers[iDDL] = new AliAltroBuffer(fileName.Data(),fMapping[iRCUside]);
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

  loader->UnloadDigits();
}

//____________________________________________________________________________
void AliEMCALRawUtils::Raw2Digits(AliRawReader* reader,TClonesArray *digitsArr)
{
  // convert raw data of the current event to digits                                                                                     

  digitsArr->Clear(); 

  if (!digitsArr) {
    Error("Raw2Digits", "no digits found !");
    return;
  }
  if (!reader) {
    Error("Raw2Digits", "no raw reader found !");
    return;
  }

  AliCaloRawStream in(reader,"EMCAL",fMapping);
  // Select EMCAL DDL's;
  reader->Select("EMCAL");

  //Updated fitting routine from 2007 beam test takes into account
  //possibility of two peaks in data and selects first one for fitting
  //Also sets some of the starting parameters based on the shape of the
  //given raw signal being fit

  TF1 * signalF = new TF1("signal", RawResponseFunction, 0, GetRawFormatTimeBins(), 5);
  signalF->SetParameters(10.,0.,fTau,fOrder,5.); //set all defaults once, just to be safe
  signalF->SetParNames("amp","t0","tau","N","ped");
  signalF->SetParameter(2,fTau); // tau in units of time bin
  signalF->SetParLimits(2,2,-1);
  signalF->SetParameter(3,fOrder); // order
  signalF->SetParLimits(3,2,-1);
  
  Int_t id =  -1;
  Float_t time = 0. ; 
  Float_t amp = 0. ; 

  //Graph to hold data we will fit (should be converted to an array
  //later to speed up processing
  TGraph * gSig = new TGraph(GetRawFormatTimeBins()); 

  Int_t readOk = 1;
  Int_t lowGain = 0;

  while (readOk && in.GetModule() < 0) 
    readOk = in.Next();  // Go to first digit

  Int_t col = 0;
  Int_t row = 0;

  while (readOk) { 

    id =  fGeom->GetAbsCellIdFromCellIndexes(in.GetModule(), in.GetRow(), in.GetColumn()) ;
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
      gSig->SetPoint(i, i , 0);
    }

    Int_t iTime = 0;
    do {
      if (in.GetTime() >= gSig->GetN()) {
	  AliWarning("Too many time bins");
	  gSig->Set(in.GetTime());
      }
      col = in.GetColumn();
      row = in.GetRow();
      
      gSig->SetPoint(in.GetTime(), in.GetTime(), in.GetSignal()) ;

      if (in.GetTime() > maxTime)
        maxTime = in.GetTime();
      iTime++;
    } while ((readOk = in.Next()) && !in.IsNewHWAddress());

    FitRaw(gSig, signalF, amp, time) ; 
    
    if (amp > 0 && amp < 2000) {  //check both high and low end of
				   //result, 2000 is somewhat arbitrary
      AliDebug(2,Form("id %d lowGain %d amp %g", id, lowGain, amp));

      AddDigit(digitsArr, id, lowGain, (Int_t)amp, time);
    }
	
    // Reset graph
    for (Int_t index = 0; index < gSig->GetN(); index++) {
      gSig->SetPoint(index, index, 0) ;  
    } 
    // Reset starting parameters for fit function
    signalF->SetParameters(10.,0.,fTau,fOrder,5.); //reset all defaults just to be safe

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
         // (use high gain if signal < cut value, otherwise low gain)
    if (lowGain) { // new digit is low gain
      if (digit->GetAmp() > fgkOverflowCut) {  // use if stored digit is out of range
	digit->SetAmp(Int_t(fHighLowGainFactor * amp));
	digit->SetTime(time);
      }
    }
    else if (amp < fgkOverflowCut) { // new digit is high gain; use if not out of range
      digit->SetAmp(amp);
      digit->SetTime(time);
    }
  }
}

//____________________________________________________________________________ 
void AliEMCALRawUtils::FitRaw(TGraph * gSig, TF1* signalF, Float_t & amp, Float_t & time)
{
  // Fits the raw signal time distribution; from AliEMCALGetter 

  amp = time = 0. ; 
  Double_t ped = 0;
  Int_t nPed = 0;

  for (Int_t index = 0; index < fNPedSamples; index++) {
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
    AliWarning("Could not determine pedestal");	  
    ped = 10; // put some small value as first guess
  }

  Int_t max_found = 0;
  Int_t i_max = 0;
  Float_t max = -1;
  Float_t max_fit = gSig->GetN();
  Float_t min_after_sig = 9999;
  Int_t tmin_after_sig = gSig->GetN();
  Int_t n_ped_after_sig = 0;
  Int_t plateau_width = 0;
  Int_t plateau_start = 9999;
  Float_t Cut = 0.3;

  for (Int_t i=fNPedSamples; i < gSig->GetN(); i++) {
    Double_t ttime, signal;
    gSig->GetPoint(i, ttime, signal) ; 
    if (!max_found && signal > max) {
      i_max = i;
      max = signal;
    }
    else if ( max > ped + fNoiseThreshold ) {
      max_found = 1;
      min_after_sig = signal;
      tmin_after_sig = i;
    }
    if (max_found) {
      if ( signal < min_after_sig) {
        min_after_sig = signal;
	tmin_after_sig = i;
      }
      if (i > tmin_after_sig + 5) {  // Two close peaks; end fit at minimum
        max_fit = tmin_after_sig;
        break;
      }
      if ( signal < Cut*max){   //stop fit at 30% amplitude(avoid the pulse shape falling edge)
        max_fit = i;
        break;
      }
      if ( signal < ped + fNoiseThreshold)
        n_ped_after_sig++;
      if (n_ped_after_sig >= 5) {  // include 5 pedestal bins after peak
        max_fit = i;
        break;
      }
    }
    //Add check on plateau
    if (signal >= fgkRawSignalOverflow - fNoiseThreshold) {
      if(plateau_width == 0) plateau_start = i;
      plateau_width++;
    }
  }

  if(plateau_width > 0) {
    for(int j = 0; j < plateau_width; j++) {
      //Note, have to remove the same point N times because after each
      //remove, the positions of all subsequent points have shifted down
      gSig->RemovePoint(plateau_start);
    }
  }

  if ( max - ped > fNoiseThreshold ) { // else its noise 
    AliDebug(2,Form("Fitting max %d ped %d", max, ped));
    signalF->SetRange(0,max_fit);

    if(max-ped > 50) 
      signalF->SetParLimits(2,1,3);

    signalF->SetParameter(4, ped) ; 
    signalF->SetParameter(1, i_max);
    signalF->SetParameter(0, max);
    
    gSig->Fit(signalF, "QROW"); // Note option 'W': equal errors on all points
    amp = signalF->GetParameter(0); 
    time = signalF->GetParameter(1)*GetRawFormatTimeBinWidth() - fgTimeTrigger;
  }
  return;
}
//__________________________________________________________________
Double_t AliEMCALRawUtils::RawResponseFunction(Double_t *x, Double_t *par)
{
  // Matches version used in 2007 beam test
  //
  // Shape of the electronics raw reponse:
  // It is a semi-gaussian, 2nd order Gamma function of the general form
  //
  // t' = (t - t0 + tau) / tau
  // F = A * t**N * exp( N * ( 1 - t) )   for t >= 0
  // F = 0                                for t < 0 
  //
  // parameters:
  // A:   par[0]   // Amplitude = peak value
  // t0:  par[1]
  // tau: par[2]
  // N:   par[3]
  // ped: par[4]
  //
  Double_t signal ;
  Double_t tau =par[2];
  Double_t N =par[3];
  Double_t ped = par[4];
  Double_t xx = ( x[0] - par[1] + tau ) / tau ;

  if (xx <= 0) 
    signal = ped ;  
  else {  
    signal = ped + par[0] * TMath::Power(xx , N) * TMath::Exp(N * (1 - xx )) ; 
  }
  return signal ;  
}

//__________________________________________________________________
Bool_t AliEMCALRawUtils::RawSampledResponse(
const Double_t dtime, const Double_t damp, Int_t * adcH, Int_t * adcL) const 
{
  // for a start time dtime and an amplitude damp given by digit, 
  // calculates the raw sampled response AliEMCAL::RawResponseFunction

  Bool_t lowGain = kFALSE ; 

  // A:   par[0]   // Amplitude = peak value
  // t0:  par[1]                            
  // tau: par[2]                            
  // N:   par[3]                            
  // ped: par[4]

  TF1 signalF("signal", RawResponseFunction, 0, GetRawFormatTimeBins(), 5);
  signalF.SetParameter(0, damp) ; 
  signalF.SetParameter(1, (dtime + fgTimeTrigger)/fgTimeBinWidth) ; 
  signalF.SetParameter(2, fTau) ; 
  signalF.SetParameter(3, fOrder);
  signalF.SetParameter(4, fgPedestalValue);

  for (Int_t iTime = 0; iTime < GetRawFormatTimeBins(); iTime++) {
    Double_t signal = signalF.Eval(iTime) ;     

    //According to Terry Awes, 13-Apr-2008
    //add gaussian noise in quadrature to each sample
    //Double_t noise = gRandom->Gaus(0.,fgFEENoise);
    //signal = sqrt(signal*signal + noise*noise);

    adcH[iTime] =  static_cast<Int_t>(signal + 0.5) ;
    if ( adcH[iTime] > fgkRawSignalOverflow ){  // larger than 10 bits 
      adcH[iTime] = fgkRawSignalOverflow ;
      lowGain = kTRUE ; 
    }

    signal /= fHighLowGainFactor;

    adcL[iTime] =  static_cast<Int_t>(signal + 0.5) ;
    if ( adcL[iTime] > fgkRawSignalOverflow)  // larger than 10 bits 
      adcL[iTime] = fgkRawSignalOverflow ;
  }
  return lowGain ; 
}
