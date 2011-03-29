// -*- mode: c++ -*-
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


//_________________________________________________________________________
//  Utility Class for handling Raw data
//  Does all transitions from Digits to Raw and vice versa, 
//  for simu and reconstruction
//
//  Note: the current version is still simplified. Only 
//    one raw signal per digit is generated; either high-gain or low-gain
//    Need to add concurrent high and low-gain info in the future
//    No pedestal is added to the raw signal.
//*-- Author: Marco van Leeuwen (LBL)


#include "AliEMCALRawUtils.h"
#include "TF1.h"
#include "TGraph.h"
#include <TRandom.h>
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliAltroBuffer.h"
#include "AliRawReader.h"
#include "AliCaloRawStreamV3.h"
#include "AliDAQ.h"
#include "AliEMCALRecParam.h"
#include "AliEMCALLoader.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRawDigit.h"
#include "AliEMCAL.h"
#include "AliCaloCalibPedestal.h"  
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliEMCALTriggerRawDigitMaker.h"
#include "AliEMCALTriggerSTURawStream.h"
#include "AliEMCALTriggerData.h"
#include "AliCaloConstants.h"
#include "AliCaloRawAnalyzer.h"
#include "AliCaloRawAnalyzerFactory.h"

using namespace CALO;
using namespace EMCAL;

Double_t AliEMCALRawUtils::fgTimeTrigger  = 600E-9 ;   // the time of the trigger as approximately seen in the data
Int_t    AliEMCALRawUtils::fgThreshold         = 1;
Int_t    AliEMCALRawUtils::fgPedestalValue     = 0;  // pedestal value for digits2raw, default generate ZS data
Double_t AliEMCALRawUtils::fgFEENoise          = 3.; // 3 ADC channels of noise (sampled)

ClassImp(AliEMCALRawUtils)


AliEMCALRawUtils::AliEMCALRawUtils( Algo::fitAlgorithm fitAlgo) : fNoiseThreshold(3),
								  fNPedSamples(4), 
								  fGeom(0), 
								  fOption(""),
								  fRemoveBadChannels(kFALSE),
								  fFittingAlgorithm(0),  
								  fTimeMin(-1.),
								  fTimeMax(1.),
								  fUseFALTRO(kTRUE),
								  fRawAnalyzer(0),
								  fTriggerRawDigitMaker(0x0)
{
  SetFittingAlgorithm(fitAlgo);
  // SetFittingAlgorithm(  Algo::kLMSOffline);

  //Get Mapping RCU files from the AliEMCALRecParam                                 
  const TObjArray* maps = AliEMCALRecParam::GetMappings();
  if(!maps) AliFatal("Cannot retrieve ALTRO mappings!!");

  for(Int_t i = 0; i < 4; i++) {
    fMapping[i] = (AliAltroMapping*)maps->At(i);
  }

  //To make sure we match with the geometry in a simulation file,
  //let's try to get it first.  If not, take the default geometry
  AliRunLoader *rl = AliRunLoader::Instance();
  if (rl && rl->GetAliRun()) {
    AliEMCAL * emcal = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"));
    if(emcal)fGeom = emcal->GetGeometry();
    else {
      AliDebug(1, Form("Using default geometry in raw reco"));
      fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
    }

  } else {
    AliDebug(1, Form("Using default geometry in raw reco"));
    fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
  }

  if(!fGeom) AliFatal(Form("Could not get geometry!"));
	
  fTriggerRawDigitMaker = new AliEMCALTriggerRawDigitMaker();

}


//____________________________________________________________________________
AliEMCALRawUtils::AliEMCALRawUtils(AliEMCALGeometry *pGeometry, Algo::fitAlgorithm fitAlgo) : //fHighLowGainFactor(16.), 
											      //fOrder(2), 
  //  fTau(2.35), 
  fNoiseThreshold(3),
  fNPedSamples(4), 
  fGeom(pGeometry), 
  fOption(""),
  fRemoveBadChannels(kFALSE),fFittingAlgorithm(0),
  fTimeMin(-1.),fTimeMax(1.),
  fUseFALTRO(kTRUE),fRawAnalyzer(0),
  fTriggerRawDigitMaker(0x0)
{
 

 // Initialize with the given geometry - constructor required by HLT
  // HLT does not use/support AliRunLoader(s) instances
  // This is a minimum intervention solution
  // Comment by MPloskon@lbl.gov
  SetFittingAlgorithm(fitAlgo);
  // SetFittingAlgorithm(  Algo::kLMSOffline);
  //Get Mapping RCU files from the AliEMCALRecParam
  const TObjArray* maps = AliEMCALRecParam::GetMappings();
  if(!maps) AliFatal("Cannot retrieve ALTRO mappings!!");
  
  for(Int_t i = 0; i < 4; i++) 
    {
      fMapping[i] = (AliAltroMapping*)maps->At(i);
    }

  if(!fGeom) AliFatal(Form("Could not get geometry!"));
  fTriggerRawDigitMaker = new AliEMCALTriggerRawDigitMaker();	
}


//____________________________________________________________________________
AliEMCALRawUtils::~AliEMCALRawUtils() 
{
  //dtor
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
  
  TArrayI adcValuesLow( TIMEBINS );
  TArrayI adcValuesHigh( TIMEBINS );
  
  // loop over digits (assume ordered digits)
  for (Int_t iDigit = 0; iDigit < digits->GetEntries(); iDigit++) 
    {
      AliEMCALDigit* digit = dynamic_cast<AliEMCALDigit *>(digits->At(iDigit)) ;
      if(!digit)
	{
	  AliFatal("NULL Digit");
	}
      else
	{
	  if (digit->GetAmplitude() < fgThreshold) 
	    {
	      continue;
	    }
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
      Int_t iDDL = NRCUSPERMODULE*nSM + iRCU;
      if (iDDL < 0 || iDDL >= nDDL){
        Fatal("Digits2Raw()","Non-existent DDL board number: %d", iDDL);
      }
      else{
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
        if (digit->GetTimeR() >  TIMEBINMAX  ) {
          AliInfo("Signal is out of time range.\n");
          buffers[iDDL]->FillBuffer((Int_t)digit->GetAmplitude());
          buffers[iDDL]->FillBuffer( TIMEBINS );  // time bin
          buffers[iDDL]->FillBuffer(3);          // bunch length      
          buffers[iDDL]->WriteTrailer(3, ieta, iphi, nSM);  // trailer
          // calculate the time response function
        } else {
          Bool_t lowgain = RawSampledResponse(digit->GetTimeR(), digit->GetAmplitude(), adcValuesHigh.GetArray(), adcValuesLow.GetArray()) ; 
          if (lowgain) 
            buffers[iDDL]->WriteChannel(ieta, iphi, 0, TIMEBINS, adcValuesLow.GetArray(), fgThreshold);
          else 
            buffers[iDDL]->WriteChannel(ieta,iphi, 1, TIMEBINS, adcValuesHigh.GetArray(), fgThreshold);
        }
      }// iDDL under the limits
    }//digit exists
  }//Digit loop
  
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
void AliEMCALRawUtils::AddDigit(TClonesArray *digitsArr, Int_t id, Int_t lowGain, Float_t amp, Float_t time, Float_t chi2, Int_t ndf) 
{
  AliEMCALDigit *digit = 0, *tmpdigit = 0;
  TIter nextdigit(digitsArr);
 
  while (digit == 0 && (tmpdigit = (AliEMCALDigit*) nextdigit())) 
    {
      if (tmpdigit->GetId() == id) digit = tmpdigit;
    }
  
  if (!digit) { // no digit existed for this tower; create one
    Int_t type = AliEMCALDigit::kHG; // use enum in AliEMCALDigit
    if (lowGain) 
      { 
	amp *= HGLGFACTOR;
	type = AliEMCALDigit::kLGnoHG;
      } 
    
    Int_t idigit = digitsArr->GetEntries();
    new((*digitsArr)[idigit]) AliEMCALDigit( -1, -1, id, amp, time, type, idigit, chi2, ndf); 
    AliDebug(2,Form("Add digit Id %d for the first time, type %d", id, type));
  }//digit added first time
  else 
    { // a digit already exists, check range 
		// (use high gain if signal < cut value, otherwise low gain)
      if (lowGain) 
	{ // new digit is low gain
	  if (digit->GetAmplitude() >  OVERFLOWCUT ) 
	    {  // use if previously stored (HG) digit is out of range
	      digit->SetAmplitude( HGLGFACTOR * amp);
	      digit->SetTime(time);
	      digit->SetType(AliEMCALDigit::kLG);
	      AliDebug(2,Form("Add LG digit ID %d for the second time, type %d", digit->GetId(), digit->GetType()));
	    }
	}//new low gain digit
      else { // new digit is high gain 
	
	if (amp <  OVERFLOWCUT  ) 
	  { // new digit is high gain; use if not out of range
	    digit->SetAmplitude(amp);
	    digit->SetTime(time);
	    digit->SetType(AliEMCALDigit::kHG);
	    AliDebug(2,Form("Add HG digit ID %d for the second time, type %d", digit->GetId(), digit->GetType()));
	  }
	else 
	  { // HG out of range, just change flag value to show that HG did exist
	    digit->SetType(AliEMCALDigit::kLG);
	    AliDebug(2,Form("Change LG digit to HG, ID %d, type %d", digit->GetId(), digit->GetType()));
	  }
      }//new high gain digit
    }//digit existed replace it
}


//____________________________________________________________________________
void AliEMCALRawUtils::Raw2Digits(AliRawReader* reader,TClonesArray *digitsArr, const AliCaloCalibPedestal* pedbadmap, TClonesArray *digitsTRG, AliEMCALTriggerData* trgData)
{

  if(digitsArr) digitsArr->Clear("C"); 
  if (!digitsArr) { Error("Raw2Digits", "no digits found !");return;}
  if (!reader) {Error("Raw2Digits", "no raw reader found !");return;}
  AliEMCALTriggerSTURawStream inSTU(reader);
  AliCaloRawStreamV3 in(reader,"EMCAL",fMapping);	
  reader->Select("EMCAL",0,43); // 43 = AliEMCALGeoParams::fgkLastAltroDDL
  fTriggerRawDigitMaker->Reset();	
  fTriggerRawDigitMaker->SetIO(reader, in, inSTU, digitsTRG, trgData);
  fRawAnalyzer->SetIsZeroSuppressed(true); // TMP - should use stream->IsZeroSuppressed(), or altro cfg registers later
    
  Int_t lowGain  = 0;
  Int_t caloFlag = 0; // low, high gain, or TRU, or LED ref.
  
  while (in.NextDDL()) 
    {
      //  fprintf(fp," TP1\n");
      while (in.NextChannel()) 
	{
	  //	  fprintf(fp," TP2\n");
    	  caloFlag = in.GetCaloFlag();
	  if (caloFlag > 2) continue; // Work with ALTRO and FALTRO 
    	  if(caloFlag < 2 && fRemoveBadChannels && pedbadmap->IsBadChannel(in.GetModule(),in.GetColumn(),in.GetRow()))
	    {
	      continue;
	    }  
      	  vector<AliCaloBunchInfo> bunchlist; 
	  while (in.NextBunch()) 
	    {
	      bunchlist.push_back( AliCaloBunchInfo(in.GetStartTimeBin(), in.GetBunchLength(), in.GetSignals() ) );
	    } 
	  if (bunchlist.size() == 0) continue;
      	  if ( caloFlag < 2 )
	    { // ALTRO
	      Int_t id = fGeom->GetAbsCellIdFromCellIndexes(in.GetModule(), in.GetRow(), in.GetColumn()) ;
	      lowGain  = in.IsLowGain();
	      fRawAnalyzer->SetL1Phase( in.GetL1Phase() );
	      AliCaloFitResults res =  fRawAnalyzer->Evaluate( bunchlist, in.GetAltroCFG1(), in.GetAltroCFG2());  
	      if(res.GetAmp() >= fNoiseThreshold )
		{
		  AddDigit(digitsArr, id, lowGain, res.GetAmp(),  res.GetTime(), res.GetChi2(),  res.GetNdf() ); 
		}
	    }//ALTRO
	  else if(fUseFALTRO)
	    {// Fake ALTRO
	      fTriggerRawDigitMaker->Add( bunchlist );
	    }//Fake ALTRO
	} // end while over channel   
    } //end while over DDL's, of input stream 
  fTriggerRawDigitMaker->PostProcess();	
  TrimDigits(digitsArr);
}


void AliEMCALRawUtils::TrimDigits(TClonesArray *digitsArr) 
{
  AliEMCALDigit *digit = 0;
  Int_t n = 0;
  Int_t nDigits = digitsArr->GetEntriesFast();
  TIter nextdigit(digitsArr);
  while ((digit = (AliEMCALDigit*) nextdigit())) {
    if (digit->GetType() == AliEMCALDigit::kLGnoHG) {
      AliDebug(1,Form("Remove digit with id %d, LGnoHG",digit->GetId()));
      digitsArr->Remove(digit);
    }
    else if(fTimeMin > digit->GetTime() || fTimeMax < digit->GetTime()) {
      digitsArr->Remove(digit);
      AliDebug(1,Form("Remove digit with id %d, Bad Time %e",digit->GetId(), digit->GetTime()));
    }
    else if (0 > digit->GetChi2()) {
      digitsArr->Remove(digit);
      AliDebug(1,Form("Remove digit with id %d, Bad Chi2 %e",digit->GetId(), digit->GetChi2()));
    }
    else {
      digit->SetIndexInList(n);	
      n++;
    }    
  }//while
  
  digitsArr->Compress();
  AliDebug(1,Form("N Digits before trimming : %d; after array compression %d",nDigits,digitsArr->GetEntriesFast()));
}


Double_t 
AliEMCALRawUtils::RawResponseFunction(Double_t *x, Double_t *par)
{
  Double_t signal = 0.;
  Double_t tau    = par[2];
  Double_t n      = par[3];
  Double_t ped    = par[4];
  Double_t xx     = ( x[0] - par[1] + tau ) / tau ;
  if (xx <= 0) 
    signal = ped ;  
  else 
    {  
      signal = ped + par[0] * TMath::Power(xx , n) * TMath::Exp(n * (1 - xx )) ; 
    }
  return signal ;  
}


Bool_t AliEMCALRawUtils::RawSampledResponse(const Double_t dtime, const Double_t damp, Int_t * adcH, 
					    Int_t * adcL, const Int_t keyErr) const 
{
  Bool_t lowGain = kFALSE ; 
  TF1 signalF("signal", RawResponseFunction, 0, TIMEBINS, 5);
  signalF.SetParameter(0, damp) ; 
  signalF.SetParameter(1, (dtime + fgTimeTrigger)/ TIMEBINWITH) ; 
  signalF.SetParameter(2, TAU) ; 
  signalF.SetParameter(3, ORDER);
  signalF.SetParameter(4, fgPedestalValue);
	
  Double_t signal=0.0, noise=0.0;
  for (Int_t iTime = 0; iTime <  TIMEBINS; iTime++) {
    signal = signalF.Eval(iTime) ;  
    
    if(keyErr>0) {
      noise = gRandom->Gaus(0.,fgFEENoise);
      signal += noise; 
    }
	  
    adcH[iTime] =  static_cast<Int_t>(signal + 0.5) ;
    if ( adcH[iTime] > MAXBINVALUE ){  // larger than 10 bits 
      adcH[iTime] = MAXBINVALUE ;
      lowGain = kTRUE ; 
    }
    signal /= HGLGFACTOR;
    adcL[iTime] =  static_cast<Int_t>(signal + 0.5) ;
    if ( adcL[iTime] > MAXBINVALUE )  // larger than 10 bits 
      adcL[iTime] = MAXBINVALUE ;
  }
  return lowGain ; 
}


//__________________________________________________________________
Double_t AliEMCALRawUtils::RawResponseFunctionLog(Double_t *x, Double_t *par)
{
  Double_t signal = 0. ;
  Double_t tau    = par[2];
  Double_t n      = par[3];
  Double_t xx     = ( x[0] - par[1] + tau ) / tau ;

  if (xx < 0)
    { 
      signal = par[0] - n*TMath::Log(TMath::Abs(xx)) + n * (1 - xx ) ;  
    }
  else 
    {  
      signal = par[0] + n*TMath::Log(xx) + n * (1 - xx ) ; 
    }
  return signal ;  
}



//__________________________________________________________________
void AliEMCALRawUtils::SetFittingAlgorithm(Int_t fitAlgo)              
{
  // fRawAnalyzer = AliCaloRawAnalyzerFactory::CreateAnalyzer(  Algo::kStandard );
  fRawAnalyzer = AliCaloRawAnalyzerFactory::CreateAnalyzer( fitAlgo );
  
  //fRawAnalyzer = AliCaloRawAnalyzerFactory::CreateAnalyzer( kStandard );

  fRawAnalyzer->SetNsampleCut(5); // requirement for fits to be done, for the new methods
  fRawAnalyzer->SetOverflowCut ( OVERFLOWCUT );
  fRawAnalyzer->SetAmpCut(fNoiseThreshold);
  fRawAnalyzer->SetFitArrayCut(fNoiseThreshold);
  //  return;
}



