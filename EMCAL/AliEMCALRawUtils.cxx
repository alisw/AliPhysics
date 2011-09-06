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
//*-- Major refactoring by Per Thomas Hille

#include "AliEMCALRawUtils.h"
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
#include "AliEMCALRawResponse.h"

using namespace CALO;
using namespace EMCAL;

ClassImp(AliEMCALRawUtils)


AliEMCALRawUtils::AliEMCALRawUtils( Algo::fitAlgorithm fitAlgo) : fNoiseThreshold(3),
								  fNPedSamples(4), 
								  fGeom(0), 
								  fOption(""),
								  fRemoveBadChannels(kFALSE),
								  fFittingAlgorithm(fitAlgo),  
								  fTimeMin(-1.),
								  fTimeMax(1.),
								  fUseFALTRO(kTRUE),
								  fRawAnalyzer(0),
								  fTriggerRawDigitMaker(0x0)
{
  SetFittingAlgorithm(fitAlgo);
  const TObjArray* maps = AliEMCALRecParam::GetMappings();
  if(!maps) AliFatal("Cannot retrieve ALTRO mappings!!");
  for(Int_t i = 0; i < 4; i++) 
    {
      fMapping[i] = (AliAltroMapping*)maps->At(i);
    }
  
  AliRunLoader *rl = AliRunLoader::Instance();
  if (rl && rl->GetAliRun()) 
    {
    AliEMCAL * emcal = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"));
    if(emcal)
      {
	fGeom = emcal->GetGeometry();
      }
    else 
      {
	AliDebug(1, Form("Using default geometry in raw reco"));
	fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
      }
    } 
  else 
    {
      AliDebug(1, Form("Using default geometry in raw reco"));
      fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
    }
  
  if(!fGeom) AliFatal(Form("Could not get geometry!"));
  fTriggerRawDigitMaker = new AliEMCALTriggerRawDigitMaker();
}


AliEMCALRawUtils::~AliEMCALRawUtils() 
{
  //dtor
  delete fRawAnalyzer;
  delete fTriggerRawDigitMaker;
}


void AliEMCALRawUtils::Digits2Raw()
{
  // convert digits of the current event to raw data
  AliRunLoader *rl = AliRunLoader::Instance();
  AliEMCALLoader *loader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));
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
	  if (digit->GetAmplitude() <  AliEMCALRawResponse::GetRawFormatThreshold() ) 
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
	  if (0<=iphi&&iphi<8) iRCU=0; // first cable row
	  else if (8<=iphi&&iphi<16 && 0<=ieta&&ieta<24) iRCU=0; // first half; 
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
	    if (buffers[iDDL] == 0) 
	      {      
		// open new file and write dummy header
		TString fileName = AliDAQ::DdlFileName("EMCAL",iDDL);
		//Select mapping file RCU0A, RCU0C, RCU1A, RCU1C
          Int_t iRCUside=iRCU+(nSM%2)*2;
          //iRCU=0 and even (0) SM -> RCU0A.data   0
          //iRCU=1 and even (0) SM -> RCU1A.data   1
          //iRCU=0 and odd  (1) SM -> RCU0C.data   2
          //iRCU=1 and odd  (1) SM -> RCU1C.data   3
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
          Bool_t lowgain = AliEMCALRawResponse::RawSampledResponse(digit->GetTimeR(), digit->GetAmplitude(), 
								   adcValuesHigh.GetArray(), adcValuesLow.GetArray()) ; 
	  
	  if (lowgain) 
            buffers[iDDL]->WriteChannel(ieta, iphi, 0, TIMEBINS, adcValuesLow.GetArray(),  AliEMCALRawResponse::GetRawFormatThreshold()  );
          else 
            buffers[iDDL]->WriteChannel(ieta,iphi, 1, TIMEBINS, adcValuesHigh.GetArray(),  AliEMCALRawResponse::GetRawFormatThreshold()  );
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



void AliEMCALRawUtils::AddDigit(TClonesArray *digitsArr, Int_t id, Int_t lowGain, Float_t amp, Float_t time, Float_t chi2, Int_t ndf) 
{
  // comment
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


void AliEMCALRawUtils::Raw2Digits(AliRawReader* reader,TClonesArray *digitsArr, const AliCaloCalibPedestal* pedbadmap, TClonesArray *digitsTRG, AliEMCALTriggerData* trgData)
{
  //conversion of raw data to digits
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
  
  Float_t bcTimePhaseCorr = 0; // for BC-based L1 phase correction
  Int_t bcMod4 = (reader->GetBCID() % 4); // LHC uses 40 MHz, EMCal uses 10 MHz clock
  if (bcMod4==0 || bcMod4==1) { 
    bcTimePhaseCorr = -1e-7; // subtract 100 ns for certain BC values
  } 

  while (in.NextDDL()) 
    {
      while (in.NextChannel()) 
	{
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
		  AddDigit(digitsArr, id, lowGain, res.GetAmp(),  res.GetTime()+bcTimePhaseCorr, res.GetChi2(),  res.GetNdf() ); 
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


void AliEMCALRawUtils::SetFittingAlgorithm(Int_t fitAlgo)              
{
  delete fRawAnalyzer; // delete doesn't do anything if the pointer is 0x0
  fRawAnalyzer = AliCaloRawAnalyzerFactory::CreateAnalyzer( fitAlgo );
  fRawAnalyzer->SetNsampleCut(5); // requirement for fits to be done, for the new methods
  fRawAnalyzer->SetOverflowCut ( OVERFLOWCUT );
  fRawAnalyzer->SetAmpCut(fNoiseThreshold);
  fRawAnalyzer->SetFitArrayCut(fNoiseThreshold);
}



