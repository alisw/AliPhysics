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
//This class produces PHOS trigger digits of one event.
//Authors: Henrik Qvigstad, Boris Polishchuk.

#include "AliPHOSTriggerRawDigiProducer.h"
#include "AliPHOSTriggerRawReader.h"
#include "AliPHOSTRURawReader.h"
#include "AliPHOSTriggerSTURawStream.h"
#include "AliPHOSTriggerParameters.h"
#include "AliPHOSTriggerRawDigit.h"
#include "AliPHOSGeometry.h"
#include "AliRawReader.h"
#include "AliCaloRawStreamV3.h"
#include "AliDAQ.h"

#include "TH1I.h"
#include "TH2I.h"

#include <iostream>
using namespace std;

ClassImp(AliPHOSTriggerRawDigiProducer)

AliPHOSTriggerRawDigiProducer::AliPHOSTriggerRawDigiProducer()
  :fModules(kNMods, false),
   fSaturationThreshold(950),
   fParameters(0),
   fRawReader(0),
   fRawStream(0),
   fTriggerReader(0)
{}

AliPHOSTriggerRawDigiProducer::AliPHOSTriggerRawDigiProducer(AliRawReader *rawReader)
  :fModules(kNMods, false),
   fSaturationThreshold(950),
   fParameters(0),
   fRawReader(rawReader),
   fRawStream(0),
   fTriggerReader(new AliPHOSTriggerRawReader)
{
  SetAnalyseModule(1);
  SetAnalyseModule(2);
  SetAnalyseModule(3);
  SetAnalyseModule(4);

  fRawStream = new AliCaloRawStreamV3(rawReader,"PHOS");
  // Select only data in ALTRO format and skip STU, the last PHOS DDL
  rawReader->Select("PHOS",0,AliDAQ::NumberOfDdls("PHOS")-2);
}

AliPHOSTriggerRawDigiProducer::~AliPHOSTriggerRawDigiProducer()
{
  delete fRawStream;
  delete fTriggerReader;
}

void AliPHOSTriggerRawDigiProducer::ProcessEvent(TClonesArray* tdigits)
{
  ProcessL0(tdigits);
  ProcessL1(tdigits);
}

void AliPHOSTriggerRawDigiProducer::ProcessL1(TClonesArray* tdigits)
{
  AliPHOSTriggerSTURawStream inPHOSSTU(fRawReader);
  Int_t iDigit = tdigits->GetEntries();
  
  fRawReader->Reset();
  fRawReader->Select("PHOS", 20,20);

  if(inPHOSSTU.ReadPayLoad()){
    
    for(Int_t iG = 0; iG<3; iG++) {// loop over trigger thresholds: high(0), medium(1), low(2).
      for(Int_t ith=0; ith<inPHOSSTU.GetNL1GammaPatch(iG); ith++){

  	Int_t itru, ieta, iphi;
  	inPHOSSTU.GetL1GammaPatch(ith,iG,itru,ieta,iphi);
	
  	Int_t x=-1, y=-1;
  	GetGammaPatchXY(itru,ieta,iphi,x,y);
	
  	Int_t module,xloc,zloc;
  	GetL1GammaPatchModuleXZ(itru,x,y,module,xloc,zloc);

	new((*tdigits)[iDigit]) AliPHOSTriggerRawDigit(module,xloc,zloc,iG,-1);
	iDigit++;
      }
    }
  }//if(inPHOSSTU.ReadPayLoad())
  
}

void AliPHOSTriggerRawDigiProducer::GetGammaPatchXY(Int_t itru, Int_t ieta, Int_t iphi, Int_t& x, Int_t& y)
{
  x =   ieta/*xpos in TRU*/ + (int)(itru%2) * 14/*xoff in Det*/ ;
  y =   iphi/*ypos in TRU*/ + (int)(itru/2) * 8 /*yoff in Det*/ ;
}

void AliPHOSTriggerRawDigiProducer::GetL1GammaPatchModuleXZ(Int_t itru, Int_t xglob, Int_t yglob, Int_t& module, Int_t& x, Int_t& z)
{
  //convert L1 gamma patch (xglob,yglob) in Detector Global system to local module (x,z).
  //Module numeration follows the "offline" agreement.
  
  // 0 <= xglob <= 27, 0 <= yglob <= 111.
  
  if(0<=itru && itru<4) module=1;
  if(4<=itru && itru<12 ) module=2;
  if(12<=itru && itru<20 ) module=3;
  if(20<=itru && itru<28 ) module=4;
  
  z = 56-2*(xglob+1);

  Int_t offset;
  Int_t mod = module; // online numeration
  
  if(mod==1) offset= -16;
  if(mod==2) offset = 16;
  if(mod==3) offset = 16+32;
  if(mod==4) offset = 16+32+32;
  
  x = (yglob - offset)*2; 
}

void AliPHOSTriggerRawDigiProducer::ProcessL0(TClonesArray* tdigits)
{
  
  fTriggerReader->Reset();

  tdigits->Clear();
  Int_t iDigit=0 ;

  while (fRawStream->NextDDL()) {
    // Skip STU DDL
    if (fRawStream->GetDDLNumber() == fgkSTUDDL) continue; 
    while (fRawStream->NextChannel()) {
      if (fRawStream->IsTRUData()) {
	fTriggerReader->ReadFromStream(fRawStream);
      }// IsTRUData
    }// NextChannel
  }//NextDDL
	
  // Loop over modules
  for(unsigned int mod = 0; mod < fModules.size(); ++mod) {
    if( fModules[mod] ) {
      
      // Loop over 4x4 cells
      for(int TRURow = 0; TRURow < kNTRURows; ++TRURow) {
	for(int branch = 0; branch < kNBranches; ++branch) {
	  
	  AliPHOSTRURawReader* truReader = fTriggerReader->GetTRU(mod, TRURow, branch);
	  if( truReader->IsActive() ) {
	    
	    for(int xIdx = 0; xIdx < kN4x4XPrTRURow; ++xIdx) {
	      for(int zIdx = 0; zIdx < kN4x4ZPrBranch; ++zIdx) {
	      
		// Determin if Trigger is flagged for any timeBin
		bool triggered = false;

		for(int timeBin = 0; timeBin < kNTRUTimeBins; ++timeBin){
		  if(truReader->IsActive(timeBin)) {
		    if( fTriggerReader->GetTRU(mod, TRURow, branch)->GetTriggerFlag(xIdx, zIdx, timeBin) ){
		      triggered = true;
		    } // end "if TriggerBit"
		  }
		}// end TimeBin loop
		
		if( triggered ){
		  // Get peak values
		  const int TSmax = Get4x4Max(fTriggerReader, fParameters, mod, TRURow, branch, xIdx, zIdx);
		  new((*tdigits)[iDigit]) AliPHOSTriggerRawDigit(mod,xIdx,zIdx,TRURow,branch,TSmax);
		  iDigit++;
		}// end  "if triggered"
	      
	      } // end zIdx loop
	    } // end xIdx loop
	  } // truReader->IsActive
	} // end branch loop
      } // end tru loop
    } // end "if module"
  } // end mod loop
  
}

int AliPHOSTriggerRawDigiProducer::Get2x2Max(AliPHOSTriggerRawReader* reader, AliPHOSTriggerParameters* params, int mod, int xIdx, int zIdx)
{
  int max = 0;
  for(int timeBin = 0; timeBin < kNTRUTimeBins; ++timeBin) {
    const int signal = Get2x2Signal(reader, params, mod, xIdx, zIdx, timeBin);
    if( max < signal ){
      max = signal;
    }
  }
  return max;
}


int AliPHOSTriggerRawDigiProducer::Get2x2Signal(AliPHOSTriggerRawReader* reader, AliPHOSTriggerParameters* parameters, int mod, int xIdx, int zIdx, int timeBin)
{
  const int TRURow = xIdx / kN2x2XPrTRURow;
  const int branch = zIdx / kN2x2ZPrBranch;
  const int TRUX = xIdx % kN2x2XPrTRURow; // 2x2 coordinates
  const int TRUZ = zIdx % kN2x2ZPrBranch; // 2x2 coordinates

  if( reader->GetTRU(mod, TRURow, branch)->IsActive() ){
    const int signal = reader->GetTRU(mod, TRURow, branch)->GetTriggerSignal( TRUX, TRUZ, timeBin);
    if( parameters )
      return signal - parameters->GetTRUPedestal(mod, TRURow, branch, TRUX, TRUZ);
    else
      return signal - AliPHOSTRURawReader::GetDefaultSignalValue();
  }
  else
    return 0;
}

int AliPHOSTriggerRawDigiProducer::Get4x4Max(AliPHOSTriggerRawReader* reader, AliPHOSTriggerParameters* params, int mod, int TRURow, int branch, int xIdx, int zIdx)
{
  int max = 0;
  for(int timeBin = 0; timeBin < kNTRUTimeBins; ++timeBin) {
    const int signal = Get4x4Signal(reader, params, mod, TRURow, branch, xIdx, zIdx, timeBin);
    if( max < signal ){
      max = signal;
    }
  }
  return max;
}


int AliPHOSTriggerRawDigiProducer::Get4x4Signal(AliPHOSTriggerRawReader* reader, AliPHOSTriggerParameters* params, int mod, int TRURow, int branch, int xIdx, int zIdx, int timeBin)
{
  const int modX = xIdx + TRURow * kN2x2XPrTRURow;
  const int modZ = zIdx + branch * kN2x2ZPrBranch;

  const int signal
    = Get2x2Signal(reader, params, mod, modX  , modZ  , timeBin)
      + Get2x2Signal(reader, params, mod, modX+1, modZ  , timeBin)
      + Get2x2Signal(reader, params, mod, modX  , modZ+1, timeBin)
      + Get2x2Signal(reader, params, mod, modX+1, modZ+1, timeBin);
  return signal;
}

bool AliPHOSTriggerRawDigiProducer::Is2x2Active(AliPHOSTriggerRawReader* reader, int mod, int xIdx, int zIdx)
{
  const int TRURow = xIdx / kN2x2XPrTRURow;
  const int branch = zIdx / kN2x2ZPrBranch;

  return reader->GetTRU(mod, TRURow, branch)->IsActive();
}

bool AliPHOSTriggerRawDigiProducer::Is2x2Active(AliPHOSTriggerRawReader* reader, int mod, int xIdx, int zIdx, int timeBin)
{
  const int TRURow = xIdx / kN2x2XPrTRURow;
  const int branch = zIdx / kN2x2ZPrBranch;

  return reader->GetTRU(mod, TRURow, branch)->IsActive(timeBin);
}


