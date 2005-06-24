/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
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

#include <TError.h>

#include "AliMUONTriggerCircuit.h"
#include "AliMUONTriggerDecision.h"
#include "AliMUONTriggerLut.h"
#include "AliMUON.h"
#include "AliMUONDigit.h"
#include "AliMUONConstants.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalTrigger.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliRawReader.h" // for raw data
#include "AliLog.h"


//----------------------------------------------------------------------
ClassImp(AliMUONTriggerDecision)

//----------------------------------------------------------------------
AliMUONTriggerDecision::AliMUONTriggerDecision(AliLoader* loader, Int_t iprint, AliMUONData* data)
  : TObject()
{
// Constructor 
  fDebug = iprint;            // print option
// iprint = 0 : don't print anything
// iprint = 1 : print Global Trigger Output
// iprint = 2 : print Local and Global Trigger Outputs
// iprint = 3 : iprint = 2 + detailed info on X strips
// iprint = 4 : iprint = 2 + detailed info on Y strip
// iprint = 5 : iprint = 2 + detailed info on X and Y strips
// Note : with iprint>2, the strips detailed info is given for all circuits

// Global Trigger information
  Int_t i;
  Int_t icirc;
  Int_t istrip;

  for (i=0; i<3; i++) {   // [0] : Low pt, [1] : High pt, [2] : All pt 
    fGlobalSinglePlus[i]=0;     // tot num of single plus 
    fGlobalSingleMinus[i]=0;    // tot num of single minus
    fGlobalSingleUndef[i]=0;    // tot num of single undefined
    fGlobalPairUnlike[i]=0;     // tot num of unlike-sign pairs
    fGlobalPairLike[i]=0;       // tot num of like-sign pairs
  }
  // Local Trigger information
  for (icirc=0; icirc<234; icirc++){
    fTrigger[icirc]=0;                   // trigger or not
    fStripX11[icirc]=0;                   // X strip in MC11 which triggers 
    fDev[icirc]=0;                        // deviation which triggers 
    fStripY11[icirc]=0;                   // Y strip in MC11 which triggers 
    for (i=0; i<2; i++) {           // pt information via LuT
      fLutLpt[icirc][i]=fLutHpt[icirc][i]=fLutApt[icirc][i]=0;    
    }
  }
  // bit pattern
  for (icirc=0; icirc<234; icirc++) {
    for (istrip=0; istrip<16; istrip++) {
      fXbit11[icirc][istrip]=fXbit12[icirc][istrip]=0;
      fYbit11[icirc][istrip]=fYbit12[icirc][istrip]=0;
      fYbit21[icirc][istrip]=fYbit22[icirc][istrip]=0;
      fYbit21U[icirc][istrip]=fYbit22U[icirc][istrip]=0;
      fYbit21D[icirc][istrip]=fYbit22D[icirc][istrip]=0;
    }
    for (istrip=0; istrip<32; istrip++) {
      fXbit21[icirc][istrip]=fXbit22[icirc][istrip]=0;
    }
  }

  fTriggerCircuit = new TObjArray(AliMUONConstants::NTriggerCircuit());

  // initialize loader's
  fLoader = loader;

  // initialize container
  if (data == 0){
    AliError("No MUONdata for trigger");
  }else{
    fMUONData = data;
  }

  // getting MUON
  fMUON = (AliMUON*) gAlice->GetDetector("MUON");

  // setting circuit
  for (icirc = 0; icirc < AliMUONConstants::NTriggerCircuit(); icirc++) {
    AliMUONTriggerCircuit* pCir = 0;
    pCir = &(fMUON->TriggerCircuit(icirc));
    fTriggerCircuit->AddAt(pCir, icirc);
  }

  // setting digits
  fDigits = new TObjArray(AliMUONConstants::NCh()); //NTriggerCh
  for (Int_t i=0; i<AliMUONConstants::NCh() ;i++) 
    fDigits->AddAt(new TClonesArray("AliMUONDigit",10000),i);
  fDigitIndices = new TArrayI[AliMUONConstants::NCh()];
}

//----------------------------------------------------------------------
AliMUONTriggerDecision::AliMUONTriggerDecision()
  : TObject(),
    fLoader(0),
    fTriggerCircuit(0),
    fMUONData(0),
    fMUON(0)
{
// Default constructor
  fDigitIndices = NULL;
}

//----------------------------------------------------------------------
AliMUONTriggerDecision::AliMUONTriggerDecision(const AliMUONTriggerDecision& rhs)
  : TObject(rhs) 
{
// Protected copy constructor

  AliFatal("Not implemented.");
}

//----------------------------------------------------------------------
void AliMUONTriggerDecision::ClearDigits()
{
  for ( int i=0;i<AliMUONConstants::NCh();i++ )
  {
    if ((*fDigits)[i]) ((TClonesArray*)fDigits->At(i))->Clear();
    fDigitIndices[i].Set(0);
  };
}

//----------------------------------------------------------------------
TClonesArray* AliMUONTriggerDecision::Digits(Int_t DetectionPlane)
{
  //Getting List of Digits
  if (fDigits)
    return ((TClonesArray*) fDigits->At(DetectionPlane));
  else
    return NULL;
}

//_____________________________________________________________________________
void AliMUONTriggerDecision::AddDigit(
		Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits,
		Int_t digitindex
	)
{
  //
  // Add a MUON digit to the list of Digits of the detection plane id
  // Also adds the digit index to the corresponding fDigitIndices arrays.
  //
  TClonesArray &ldigits = *Digits(id); 
  new(ldigits[ldigits.GetEntriesFast()]) AliMUONDigit(tracks,charges,digits);

  TArrayI& indices = fDigitIndices[id];
  indices.Set(indices.GetSize() + 1);
  indices[indices.GetSize() - 1] = digitindex;
}

//----------------------------------------------------------------------
AliMUONTriggerDecision::~AliMUONTriggerDecision()
{
// Destructor
  if (fTriggerCircuit){
    fTriggerCircuit->Clear();// Sets pointers to 0 since it is not the owner
    delete fTriggerCircuit;
  } 
//   if (fMUONData)
//     delete fMUONData;

  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
  }

  if (fDigitIndices)
    delete [] fDigitIndices;
}

//----------------------------------------------------------------------
AliMUONTriggerDecision& 
AliMUONTriggerDecision::operator=(const AliMUONTriggerDecision& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}    
          

//----------------------------------------------------------------------
void AliMUONTriggerDecision::Trigger(){
// main method of the class which calls the overall Trigger procedure

  ResetBit();
  SetBit();
  SetBitUpDownY();

  Int_t coinc44=0, resetMid=0; // initialize coincidence

  AliMUONTriggerCircuit* triggerCircuit;

  for (Int_t icirc=0; icirc<234; icirc++) {  // loop on circuits
    triggerCircuit = (AliMUONTriggerCircuit*)fTriggerCircuit->At(icirc);  	  
    //    Int_t idCircuit=triggerCircuit->GetIdCircuit(); 
    
    Int_t minDevStrip[5], minDev[5], coordY[5];
    for (Int_t i=0; i<5; i++) {
      minDevStrip[i]=minDev[i]=coordY[i]=0;
    }
    Int_t x2m=triggerCircuit->GetX2m();
    Int_t x2ud=triggerCircuit->GetX2ud();
    Int_t orMud[2]={0,0};
    triggerCircuit->GetOrMud(orMud);
        
// call triggerX
    TrigX(fXbit11[icirc],fXbit12[icirc],fXbit21[icirc],fXbit22[icirc], 
  	  coinc44, minDevStrip, minDev);
// call triggerY
    TrigY(fYbit11[icirc],fYbit12[icirc],fYbit21[icirc],fYbit22[icirc],
  	  fYbit21U[icirc],fYbit21D[icirc],fYbit22U[icirc],fYbit22D[icirc],
	  x2m,x2ud,orMud,resetMid,coinc44,coordY);
// call LocalTrigger     
    Int_t iTrigger=0;
    LocalTrigger(icirc, minDevStrip, minDev, coordY, iTrigger);

    if (iTrigger==1&&fDebug>1) { 
      PrintBitPatXInput(icirc);
      PrintBitPatYInput(icirc);
      PrintLocalOutput(minDevStrip, minDev, coordY);
    }      
  }  //  end loop on circuits

// call Global Trigger
  GlobalTrigger();
}

//----------------------------------------------------------------------
void AliMUONTriggerDecision::ResetBit(){
// reset bit pattern, global and local trigger output tables to 0
  
    Int_t i;
    Int_t icirc;
    Int_t istrip;

  for (icirc=0; icirc<234; icirc++) {
    for (istrip=0; istrip<16; istrip++) {
      fXbit11[icirc][istrip]=fXbit12[icirc][istrip]=0;
      fYbit11[icirc][istrip]=fYbit12[icirc][istrip]=0;
      fYbit21[icirc][istrip]=fYbit22[icirc][istrip]=0;
      fYbit21U[icirc][istrip]=fYbit22U[icirc][istrip]=0;
      fYbit21D[icirc][istrip]=fYbit22D[icirc][istrip]=0;
    }
    for (istrip=0; istrip<32; istrip++) {
      fXbit21[icirc][istrip]=fXbit22[icirc][istrip]=0;
    }
  }
  for (i=0; i<3; i++) { 
    fGlobalSinglePlus[i]=0;
    fGlobalSingleMinus[i]=0;
    fGlobalSingleUndef[i]=0;
    fGlobalPairLike[i]=0;
    fGlobalPairLike[i]=0;
  }
  for (icirc=0; icirc<234; icirc++){
    fTrigger[icirc]=0;
    fStripX11[icirc]=0;
    fDev[icirc]=0;                      
    fStripY11[icirc]=0;                 
    for (i=0; i<2; i++) {         
      fLutLpt[icirc][i]=fLutHpt[icirc][i]=fLutApt[icirc][i]=0;    
    }
  }
}

//----------------------------------------------------------------------
void AliMUONTriggerDecision::SetBit(){
// 1) loop over chambers and cathodes
// 2) load digits 
// 3) remove soft background
// 4) set the bit patterns

  Int_t cathode;
  AliMUONTriggerCircuit* triggerCircuit;

  for (Int_t chamber = 11; chamber < 15; chamber++){

      TClonesArray *muonDigits = Digits(chamber-1);
      Int_t ndigits = muonDigits->GetEntriesFast();
      AliDebug(3,Form("Found %d digits in %p %d", ndigits, (void*)muonDigits,chamber-1));

      AliMUONDigit  *mdig;
      
      for (Int_t digit = 0; digit < ndigits; digit++) {
	mdig    = (AliMUONDigit*)muonDigits->UncheckedAt(digit);
// get the center of the pad Id 
  	Int_t ix=mdig->PadX();
  	Int_t iy=mdig->PadY();
	cathode = mdig->Cathode() + 1;
	AliDebug(3,Form("cathode %d ix %d iy %d ",cathode,ix,iy));

// get the sum of the coded charge 
// see coding convention in AliMUONChamberTrigger::DisIntegration 	
	Int_t sumCharge=0;
	for (Int_t icharge=0; icharge<10; icharge++) {
	  sumCharge=sumCharge+mdig->TrackCharge(icharge);
	}
// apply condition on soft background	
	Int_t testCharge=sumCharge-(Int_t(sumCharge/10))*10;	
	if(sumCharge<=10||testCharge>0) {	  
// code pad
	  Int_t code=TMath::Abs(ix)*100+iy;
	  if (ix<0) { code=-code; }
	  
	  Int_t icirc;
	  Int_t istrip;
	  Int_t nStrip;

          // If I want to fetch the digits as in MUONCheck.C then I need to
          // know the correct digit index. These were stored in fDigitIndices
          // by the digitizer so we just need to fetch the correct value.
          Int_t digitindex = fDigitIndices[chamber-1][digit];

	  if (cathode==1) {
	    switch (chamber)
	      {
	      case 11:
		for (icirc=0; icirc<234; icirc++) {		  
		  triggerCircuit = (AliMUONTriggerCircuit*) fTriggerCircuit->At(icirc);  
		  for (istrip=0; istrip<16; istrip++) {
		    if (triggerCircuit->GetXcode(0,istrip)==code)
                    {
		      fXbit11[icirc][istrip]=1;
                      DigitFiredCircuit(icirc, cathode-1, chamber-1, digitindex);
                    };
		  }
		}
		break;
	      case 12:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = (AliMUONTriggerCircuit*)fTriggerCircuit->At(icirc);
		  for (istrip=0; istrip<16; istrip++) {
		    if (triggerCircuit->GetXcode(1,istrip)==code) 
                    {
		      fXbit12[icirc][istrip]=1;
                      DigitFiredCircuit(icirc, cathode-1, chamber-1, digitindex);
                    };
		  }
		}
		break;
	      case 13:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = (AliMUONTriggerCircuit*)fTriggerCircuit->At(icirc); 
		  for (istrip=0; istrip<32; istrip++) {
		    if (triggerCircuit->GetXcode(2,istrip)==code) 
                    {
		      fXbit21[icirc][istrip]=1;
                      DigitFiredCircuit(icirc, cathode-1, chamber-1, digitindex);
                    };
		  }
		}
		break;
	      case 14:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = (AliMUONTriggerCircuit*)fTriggerCircuit->At(icirc);
		  for (istrip=0; istrip<32; istrip++) {
		    if (triggerCircuit->GetXcode(3,istrip)==code) 
                    {
		      fXbit22[icirc][istrip]=1;		    
                      DigitFiredCircuit(icirc, cathode-1, chamber-1, digitindex);
                    };
		  }
		}		
		break;
	      }
	    
	  } else {                // Y plane 
	    switch (chamber)
	      {
	      case 11:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = (AliMUONTriggerCircuit*)fTriggerCircuit->At(icirc);
		  nStrip=triggerCircuit->GetNstripY();
		  for (istrip=0; istrip<nStrip; istrip++) {
		    if (triggerCircuit->GetYcode(0,istrip)==code) 
                    {
		      fYbit11[icirc][istrip]=1;
                      DigitFiredCircuit(icirc, cathode-1, chamber-1, digitindex);
                    };
		  }
		}
		break;
	      case 12:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = (AliMUONTriggerCircuit*)fTriggerCircuit->At(icirc);
		  nStrip=triggerCircuit->GetNstripY(); 
		  for (istrip=0; istrip<nStrip; istrip++) {
		    if (triggerCircuit->GetYcode(1,istrip)==code) 
                    {
		      fYbit12[icirc][istrip]=1;
                      DigitFiredCircuit(icirc, cathode-1, chamber-1, digitindex);
                    };
		  }
		}
		break;
	      case 13:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = (AliMUONTriggerCircuit*)fTriggerCircuit->At(icirc);
		  nStrip=triggerCircuit->GetNstripY();    
		  for (istrip=0; istrip<nStrip; istrip++) {
		    if (triggerCircuit->GetYcode(2,istrip)==code) 
                    {
		      fYbit21[icirc][istrip]=1;
                      DigitFiredCircuit(icirc, cathode-1, chamber-1, digitindex);
                    };
		  }
		}
		break;
	      case 14:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = (AliMUONTriggerCircuit*)fTriggerCircuit->At(icirc);
		  nStrip=triggerCircuit->GetNstripY();    
		  for (istrip=0; istrip<nStrip; istrip++) {
		    if (triggerCircuit->GetYcode(3,istrip)==code) 
                    {
		      fYbit22[icirc][istrip]=1;		      		 
                      DigitFiredCircuit(icirc, cathode-1, chamber-1, digitindex);
                    };
		  }
		}		
		break;
	      }
	  } // if cathode
	}  // remove soft background
      }   // end loop on digit
      fMUONData->ResetDigits();
//  }    // end loop on cathode
  }     // end loop on chamber
}  

//----------------------------------------------------------------------
void AliMUONTriggerDecision::SetBitUpDownY(){
// Set Y bit for up and down parts of circuits
  Int_t idModule, nStripX, nStripY, iPosCircuit;

  
  for (Int_t icirc=0; icirc<234; icirc++) {

    AliMUONTriggerCircuit* circuit;   // current circuit
    AliMUONTriggerCircuit* circuitD;  // circuit Down
    AliMUONTriggerCircuit* circuitU;  // circuit Up

    circuit = (AliMUONTriggerCircuit*)fTriggerCircuit->At(icirc);  
    idModule=circuit->GetIdModule();      // corresponding module Id.
    nStripX=circuit->GetNstripX();        // number of X strips
    nStripY=circuit->GetNstripY();        // number of Y strips
    iPosCircuit=circuit->GetPosCircuit(); // position of circuit in module

// fill lower part
    if (iPosCircuit==1) {               // need to scan lower module       
      if(idModule<91&&TMath::Abs(idModule)!=41&&idModule>-91) { 
	Int_t icircD=circuit->GetICircuitD();
	circuitD = (AliMUONTriggerCircuit*)fTriggerCircuit->At(icircD);  
	Int_t nStripD=circuitD->GetNstripY();
		
	if (TMath::Abs(idModule)==42) { // shift of +8 bits
	  for (Int_t istrip=0; istrip<nStripD; istrip++) {
	    fYbit21D[icirc][istrip+8]=fYbit21[icircD][istrip];
	    fYbit22D[icirc][istrip+8]=fYbit22[icircD][istrip];
	  }	  
	} else if (TMath::Abs(idModule)==52) { // shift of -8 bits
	  for (Int_t istrip=0; istrip<nStripD; istrip++) {
	    fYbit21D[icirc][istrip]=fYbit21[icircD][istrip+8];
	    fYbit22D[icirc][istrip]=fYbit22[icircD][istrip+8];
	  }
	} else {
	  for (Int_t istrip=0; istrip<nStripD; istrip++) {
	    fYbit21D[icirc][istrip]=fYbit21[icircD][istrip];
	    fYbit22D[icirc][istrip]=fYbit22[icircD][istrip];
	  }
	}
      }      
    } else {                         // lower strips within same module       
      for (Int_t istrip=0; istrip<nStripY; istrip++) { 
	fYbit21D[icirc][istrip]=fYbit21[icirc][istrip];
	fYbit22D[icirc][istrip]=fYbit22[icirc][istrip];
      }
    }    
    
// fill upper part
    if ((iPosCircuit==1&&nStripX==16)||(iPosCircuit==2&&nStripX==32)|| 
	(iPosCircuit==3&&nStripX==48)||(iPosCircuit==4&&nStripX==64)) {   
      if ((idModule>17||idModule<-17)&&TMath::Abs(idModule)!=61) {  
	Int_t icircU=circuit->GetICircuitU();
	circuitU = (AliMUONTriggerCircuit*)fTriggerCircuit->At(icircU);  
	Int_t nStripU=circuitU->GetNstripY();		
	
	if (TMath::Abs(idModule)==62) { // shift of +8 bits
	  for (Int_t istrip=0; istrip<nStripU; istrip++) {
	    fYbit21U[icirc][istrip+8]=fYbit21[icircU][istrip];
	    fYbit22U[icirc][istrip+8]=fYbit22[icircU][istrip];
	  }	  
	} else if (TMath::Abs(idModule)==52) { // shift of -8 bits
	  for (Int_t istrip=0; istrip<nStripU; istrip++) {
	    fYbit21U[icirc][istrip]=fYbit21[icircU][istrip+8];
	    fYbit22U[icirc][istrip]=fYbit22[icircU][istrip+8];
	  }
	} else {
	  for (Int_t istrip=0; istrip<nStripU; istrip++) {
	    fYbit21U[icirc][istrip]=fYbit21[icircU][istrip];
	    fYbit22U[icirc][istrip]=fYbit22[icircU][istrip];
	  }
	}
      }      
    } else {                       // upper strips within same module       
      for (Int_t istrip=0; istrip<nStripY; istrip++) { 
	fYbit21U[icirc][istrip]=fYbit21[icirc][istrip];
	fYbit22U[icirc][istrip]=fYbit22[icirc][istrip];
      }
    } 
  } // loop on circuit
}

//----------------------------------------------------------------------
// x part of trigger Algo
//----------------------------------------------------------------------
//----------------------------------------------------------------------
void AliMUONTriggerDecision::TrigX(Int_t ch1q[16], Int_t ch2q[16], 
				   Int_t ch3q[32], Int_t ch4q[32], 
				   Int_t coinc44, Int_t minDevStrip[5], 
				   Int_t minDev[5]){
// note : coinc44 = flag 0 or 1 (0 coincidence -> 3/4, 1 coincidence -> 4/4)
//---------------------------------------------------------
// step # 1 : declustering, reduction DS, calculate sgle & dble
//---------------------------------------------------------
  Int_t ch1e[19], ch2e[20], ch3e[35], ch4e[36]; 
  Int_t sgleHit1[31], sgleHit2[63];
  Int_t dbleHit1[31], dbleHit2[63];

  Int_t i;
  Int_t j;
  Int_t istrip;

  for (i=0; i<31; i++) {
    sgleHit1[i]=0;
    dbleHit1[i]=0;
  }
  for (i=0; i<63; i++) {
    sgleHit2[i]=0;
    dbleHit2[i]=0;
  }

//--- inititialize che using chq 
  for (i=0; i<19; i++) {
    if (i<1||i>16)  ch1e[i]=0; 
    else            ch1e[i]=ch1q[i-1]; 
  }
  for (i=0; i<20; i++) {
    if (i<2||i>17) ch2e[i]=0; 
    else           ch2e[i]=ch2q[i-2]; 
  }
  for (i=0; i<35; i++) {
    if (i<1||i>32) ch3e[i]=0; 
    else           ch3e[i]=ch3q[i-1];
  }
  for (i=0; i<36; i++) {
    if (i<2||i>33) ch4e[i]=0; 
    else           ch4e[i]=ch4q[i-2];
  }


//--- calculate dble & sgle first station
  for (i=0; i<=15; i++) {                   
    sgleHit1[2*i] = (!ch1e[i+1]|(ch1e[i]^ch1e[i+2])) & 
      (!ch2e[i+2] | (ch2e[i+1]^ch2e[i+3]));

    dbleHit1[2*i] = ch1e[i+1]&!(ch1e[i+2]^ch1e[i]) & 
      (ch2e[i+2] | (!ch2e[i]&ch2e[i+1]) | (ch2e[i+3]&!ch2e[i+4]));
  }

  for (i=0; i<=14; i++) {               
    sgleHit1[2*i+1] = (!ch1e[i+1]|!ch1e[i+2]|(ch1e[i]^ch1e[i+3])) & 
      (!ch2e[i+2] | !ch2e[i+3] | (ch2e[i+1]^ch2e[i+4]));
    dbleHit1[2*i+1] = ch1e[i+1]&ch1e[i+2]&!(ch1e[i]^ch1e[i+3]) & 
      (ch2e[i+2]&(!ch2e[i+1]|!ch2e[i]) | 
	      ch2e[i+3]&(ch2e[i+2]|!ch2e[i+4]|!ch2e[i+5]));
  }

//--- calculate dble & sgle second station
  for (i=0; i<=31; i++) {               
    sgleHit2[2*i] = (!ch3e[i+1]|(ch3e[i]^ch3e[i+2])) & 
      (!ch4e[i+2] | (ch4e[i+1]^ch4e[i+3]));
    dbleHit2[2*i] = ch3e[i+1]&!(ch3e[i+2]^ch3e[i]) & 
      (ch4e[i+2] | (!ch4e[i]&ch4e[i+1]) | (ch4e[i+3]&!ch4e[i+4]));
  }
  
  for (i=0; i<=30; i++) {               
    sgleHit2[2*i+1] = (!ch3e[i+1]|!ch3e[i+2]|(ch3e[i]^ch3e[i+3])) & 
      (!ch4e[i+2] | !ch4e[i+3] | (ch4e[i+1]^ch4e[i+4]));
    dbleHit2[2*i+1] = ch3e[i+1]&ch3e[i+2]&!(ch3e[i]^ch3e[i+3]) & 
      (ch4e[i+2]&(!ch4e[i+1]|!ch4e[i]) | 
       ch4e[i+3]&(ch4e[i+2]|!ch4e[i+4]|!ch4e[i+5]));
  }

//--- 
  if(fDebug==3||fDebug==5) {
    printf("===============================================================\n");
    printf(" X plane after sgle and dble \n");
    printf("                       0987654321098765432109876543210");
    printf("\n SGLE1                 ");
    for (istrip=30; istrip>=0; istrip--) printf("%i",(!sgleHit1[istrip]));
    printf("\n DBLE1                 ");
    for (istrip=30; istrip>=0; istrip--) printf("%i",dbleHit1[istrip]);
    printf("\n SGLE2 ");
    for (istrip=62; istrip>=0; istrip--) printf("%i",(!sgleHit2[istrip]));
    printf("\n DBLE2 ");
    for (istrip=62; istrip>=0; istrip--) printf("%i",dbleHit2[istrip]);
    printf("\n       210987654321098765432109876543210987654321098765432109876543210\n");
  }
  
//---------------------------------------------------------
// step # 2 : coincidence 3/4
//---------------------------------------------------------
  Int_t rearImage[31][31];
  for (i=0; i<31; i++) {
    for (j=0; j<31; j++) {
      rearImage[i][j]=0;
    }
  }

 Int_t notOr1=!dbleHit1[30] & !dbleHit1[29] & !dbleHit1[28] & !dbleHit1[27] & 
 !dbleHit1[26] & !dbleHit1[25] & !dbleHit1[24] & !dbleHit1[23] &
 !dbleHit1[22] & !dbleHit1[21] & !dbleHit1[20] & !dbleHit1[19] & 
 !dbleHit1[18] & !dbleHit1[17] & !dbleHit1[16] & !dbleHit1[15] & 
 !dbleHit1[14] & !dbleHit1[13] & !dbleHit1[12] & !dbleHit1[11] & 
 !dbleHit1[10] & !dbleHit1[9]  & !dbleHit1[8]  & !dbleHit1[7]  & 
 !dbleHit1[6]  & !dbleHit1[5]  & !dbleHit1[4]  & !dbleHit1[3]  & 
 !dbleHit1[2]  & !dbleHit1[1]  & !dbleHit1[0]  & !coinc44;

 Int_t notOr2= !dbleHit2[62] & !dbleHit2[61] & !dbleHit2[60] & !dbleHit2[59] & 
 !dbleHit2[58] & !dbleHit2[57] & !dbleHit2[56] & !dbleHit2[55] & 
 !dbleHit2[54] & !dbleHit2[53] & !dbleHit2[52] & !dbleHit2[51] & 
 !dbleHit2[50] & !dbleHit2[49] & !dbleHit2[48] & !dbleHit2[47] & 
 !dbleHit2[46] & !dbleHit2[45] & !dbleHit2[44] & !dbleHit2[43] & 
 !dbleHit2[42] & !dbleHit2[41] & !dbleHit2[40] & !dbleHit2[39] & 
 !dbleHit2[38] & !dbleHit2[37] & !dbleHit2[36] & !dbleHit2[35] & 
 !dbleHit2[34] & !dbleHit2[33] & !dbleHit2[32] & !dbleHit2[31] &
 !dbleHit2[30] & !dbleHit2[29] & !dbleHit2[28] & !dbleHit2[27] & 
 !dbleHit2[26] & !dbleHit2[25] & !dbleHit2[24] & !dbleHit2[23] & 
 !dbleHit2[22] & !dbleHit2[21] & !dbleHit2[20] & !dbleHit2[19] & 
 !dbleHit2[18] & !dbleHit2[17] & !dbleHit2[16] & !dbleHit2[15] & 
 !dbleHit2[14] & !dbleHit2[13] & !dbleHit2[12] & !dbleHit2[11] & 
 !dbleHit2[10] & !dbleHit2[9]  & !dbleHit2[8]  & !dbleHit2[7]  & 
 !dbleHit2[6]  & !dbleHit2[5]  & !dbleHit2[4]  & !dbleHit2[3]  & 
 !dbleHit2[2]  & !dbleHit2[1]  & !dbleHit2[0]  & !coinc44;	

// DS reduction
 for (i=0; i<31; i++) {
   sgleHit1[i] = !sgleHit1[i]&notOr1;
 }
 for (i=0; i<63; i++) {
   sgleHit2[i] = !sgleHit2[i]&notOr2;
 }

// extract rearImage
 for (i=0; i<31; i++){
   Int_t tmpSgleHit2[31];
   Int_t tmpDbleHit2[31];
   for (j=0; j<31; j++){
     tmpSgleHit2[j] = sgleHit2[i+j+1];
     tmpDbleHit2[j] = dbleHit2[i+j+1];
   }

   for (Int_t k=0; k<31; k++) {
     rearImage[i][k]=(sgleHit1[i]&tmpDbleHit2[k])|
       (dbleHit1[i]&(tmpSgleHit2[k]|tmpDbleHit2[k]));
   }
 }

  //-----------
 if(fDebug==3||fDebug==5) {
   printf("===============================================================\n");
   for (i=30; i>=0; i--) {
   printf("%i \t",i);
   for (istrip=31; istrip>=0; istrip--) printf("%i",rearImage[i][istrip]);
   printf("\n");   
   }
 }

//---------------------------------------------------------
// step # 3 : calculate deviation
//--------------------------------------------------------- 
 Int_t dev[31][6];
 for (i=0; i<31; i++) {
   for (j=0; j<6; j++) {
     dev[i][j]=0;
   }
 }

 for (i=0; i<31; i++){
   Int_t leftDev[5], rightDev[5]; 
   Int_t orL1, andL1, andL2, orR1, orR2, andR1, andR2, andR3;

// calculate Left deviation
 orL1=rearImage[i][16]|rearImage[i][18]|rearImage[i][20]|rearImage[i][22];
 andL1=!rearImage[i][17]&!rearImage[i][19]&!rearImage[i][21] & !orL1; 
 andL2=!rearImage[i][23]&!rearImage[i][24]&!rearImage[i][25]&!rearImage[i][26];
 
 leftDev[0] = (rearImage[i][16]|!rearImage[i][17]) & 
 (rearImage[i][16]|rearImage[i][18]|!rearImage[i][19]&
 (rearImage[i][20]|!rearImage[i][21])) &
 (orL1|!rearImage[i][23]&(rearImage[i][24]|!rearImage[i][25])) & 
 (orL1|rearImage[i][24]|rearImage[i][26]|!rearImage[i][27]&
 (rearImage[i][28]|!rearImage[i][29]));
				
 leftDev[1] = !rearImage[i][16] & 
 !(!rearImage[i][17]&!rearImage[i][18]&!rearImage[i][21]&!rearImage[i][22] & 
 (!rearImage[i][25]&!rearImage[i][26]&(rearImage[i][27]|rearImage[i][28]))) &
 (rearImage[i][17]|rearImage[i][18] | !rearImage[i][19]&!rearImage[i][20]) &
 (rearImage[i][17]|rearImage[i][18]|rearImage[i][21]|rearImage[i][22] | 
 !rearImage[i][23]&!rearImage[i][24]);
				
 leftDev[2] = (!rearImage[i][16]&!rearImage[i][17]&!rearImage[i][18]) & 
 (rearImage[i][19]|rearImage[i][20]|rearImage[i][21]|rearImage[i][22] | andL2);
		
 leftDev[3] = andL1;
		
 leftDev[4] = 
 !rearImage[i][27]&!rearImage[i][28]&!rearImage[i][29]&!rearImage[i][30] & 
 andL1 & andL2;

 // calculate Right deviation
 orR1=rearImage[i][8]|rearImage[i][10]|rearImage[i][12]|rearImage[i][14];
 orR2=rearImage[i][8]|rearImage[i][9]|rearImage[i][10]|rearImage[i][11];
 andR1=!rearImage[i][12]&!rearImage[i][13]&!rearImage[i][14]&!rearImage[i][15];
 andR2=
 !rearImage[i][8]&!rearImage[i][9]&!rearImage[i][10]&!rearImage[i][11] & andR1;
 andR3=!rearImage[i][4]&!rearImage[i][5]&!rearImage[i][6]&!rearImage[i][7]; 
		
 rightDev[0] = !rearImage[i][15]&(rearImage[i][14]|!rearImage[i][13]) & 
 ((rearImage[i][12]|rearImage[i][14]|!rearImage[i][11]&
 (rearImage[i][10]|!rearImage[i][9])) &
 ((orR1|!rearImage[i][7]&(rearImage[i][6]|!rearImage[i][5])) & 
 (orR1|rearImage[i][4]|rearImage[i][6]|!rearImage[i][3]&(rearImage[i][2]|
 !rearImage[i][1]))));
				
 rightDev[1] = !rearImage[i][15]&!rearImage[i][14] & 
 !(!rearImage[i][4]&!rearImage[i][5]&!rearImage[i][8]&!rearImage[i][9] &
 (!rearImage[i][12]&!rearImage[i][13]&(rearImage[i][2]|rearImage[i][3]))) &
 (rearImage[i][12]|rearImage[i][13] | !rearImage[i][10]&!rearImage[i][11]) & 
 (rearImage[i][8]|rearImage[i][9]|rearImage[i][12]|rearImage[i][13] | 
 !rearImage[i][6]&!rearImage[i][7]);
		
 rightDev[2] = andR1 & (orR2 | andR3); 
 rightDev[3] = andR2;		
 rightDev[4] = 
 !rearImage[i][0]&!rearImage[i][1]&!rearImage[i][2]&!rearImage[i][3] & 
 andR2 & andR3 ;

 // compare Left & Right deviations
 Int_t tmpLeftDev=0, tmpRightDev=0;
 for (j=0; j<5; j++){
   tmpLeftDev  = tmpLeftDev + Int_t(leftDev[j]*TMath::Power(2,j)); 
   tmpRightDev = tmpRightDev + Int_t(rightDev[j]*TMath::Power(2,j)); 
 }

 // assign mimimum deviation do dev[][]
 if (tmpLeftDev < tmpRightDev ){
   for (j=0; j<5; j++){ dev[i][j]=leftDev[j];}
   dev[i][5]=1;
 } else {
   for (j=0; j<5; j++){ dev[i][j]=rightDev[j];}
   dev[i][5]=0;
 }
  }
  
//---
 if(fDebug==3||fDebug==5) {
   printf("===============================================================\n");
   for (i=30; i>=0; i--) {
     printf("%i \t",i);
     for (istrip=5; istrip>=0; istrip--) printf("%i",dev[i][istrip]);
     printf(" \n");
   }
 }

//---------------------------------------------------------
// step # 4 : sort deviation
//--------------------------------------------------------- 
 Int_t bga1[16], bga2[8], bga3[4], bga4[2], bga5;
 Int_t tmpbga1[16][6], tmpbga2[8][6], tmpbga3[4][6], tmpbga4[2][6], tmpbga5[6];
 Int_t tmpMax[6]={1,1,1,1,1,0};

  for (i=0; i<15; i++) {
    Sort2x5(dev[2*i],dev[2*i+1],tmpbga1[i],bga1[i]);
  }  
    Sort2x5(dev[30],tmpMax,tmpbga1[15],bga1[15]);

//--    
  if(fDebug==3||fDebug==5) {
    printf("===============================================================\n");
    printf(" sorting : 1st level \n");
    for (i=15; i>=0; i--) {
      printf("\t %i \t",bga1[i]); 	
      for (j=5; j>=0; j--) printf("%i",tmpbga1[i][j]); 
     printf(" \n");
    }
  }

  for (i=0; i<8; i++) {  
    Sort2x5(tmpbga1[2*i],tmpbga1[2*i+1],tmpbga2[i],bga2[i]);
  }

//--    
  if(fDebug==3||fDebug==5) {
    printf("===============================================================\n");
    printf(" sorting : 2nd level \n");
    for (i=7; i>=0; i--) {
      printf("\t %i \t",bga2[i]); 	
      for (j=5; j>=0; j--) printf("%i",tmpbga1[i][j]); 	
      printf(" \n");
    }
  }
  
  for (i=0; i<4; i++) {  
    Sort2x5(tmpbga2[2*i],tmpbga2[2*i+1],tmpbga3[i],bga3[i]);
  }

//--    
  if(fDebug==3||fDebug==5) {
    printf("===============================================================\n");
    printf(" sorting : 3rd level \n");
    for (i=3; i>=0; i--) {
      printf("\t %i \t",bga3[i]); 	
      for (j=5; j>=0; j--) printf("%i",tmpbga3[i][j]); 
      printf(" \n");
    }
  }

  for (i=0; i<2; i++) {  
    Sort2x5(tmpbga3[2*i],tmpbga3[2*i+1],tmpbga4[i],bga4[i]);
  }

//--    
  if(fDebug==3||fDebug==5) {
    printf("===============================================================\n");
    printf(" sorting : 4th level \n");
    for (i=1; i>=0; i--) {
      printf("\t %i \t",bga4[i]); 	
      for (j=5; j>=0; j--) printf("%i",tmpbga4[i][j]);
      printf(" \n");
    }
  }
  
    Sort2x5(tmpbga4[0],tmpbga4[1],tmpbga5,bga5);

 // coding from 6 to 5 bits 
    minDev[4] = tmpbga5[5] | tmpbga5[4];
    for (i=0; i<4; i++) { 
      minDev[i]=tmpbga5[i] & !tmpbga5[4];
    }

 // find address of strip with minimum deviation 
    minDevStrip[4]=bga5;
    if (bga5<=1) minDevStrip[3]=bga4[bga5];

    Int_t tmpAd=minDevStrip[3]+minDevStrip[4]*2;
    if (tmpAd<=3) minDevStrip[2]=bga3[tmpAd];

    tmpAd=minDevStrip[2]+minDevStrip[3]*2+minDevStrip[4]*4;
    if (tmpAd<=7) minDevStrip[1]=bga2[tmpAd];

    tmpAd=minDevStrip[1]+minDevStrip[2]*2+minDevStrip[3]*4+minDevStrip[4]*8;
    if (tmpAd<=15) minDevStrip[0]=bga1[tmpAd];

    if(fDebug==3||fDebug==5) {
    printf("===============================================================\n");
    printf("minDevStrip = ");
    for  (i=4; i>=0; i--) printf("%i",minDevStrip[i]);
    printf(" minDev = ");
    for  (i=4; i>=0; i--) printf("%i",minDev[i]); 
    printf(" \n");
    printf("===============================================================\n");
  }

}

//---------------------------------------------
void AliMUONTriggerDecision::Sort2x5(Int_t dev1[6], Int_t dev2[6],
				     Int_t minDev[6], Int_t &dev1GTdev2){ 
// returns minimun between dev1 and dev2
 Int_t tmpDev1=0, tmpDev2=0;
 for (Int_t j=0; j<5; j++){
   tmpDev1 = tmpDev1 + Int_t(dev1[j]*TMath::Power(2,j)); 
   tmpDev2 = tmpDev2 + Int_t(dev2[j]*TMath::Power(2,j)); 
 }
 if (tmpDev1 <= tmpDev2 ){
   for (Int_t j=0; j<=5; j++) { minDev[j]=dev1[j];}
   dev1GTdev2=0;
 } else {
   for (Int_t j=0; j<=5; j++) { minDev[j]=dev2[j];}
   dev1GTdev2=1;   
 }
}

//----------------------------------------------------------------------
// y part of trigger Algo 
//----------------------------------------------------------------------
//----------------------------------------------------------------------
void AliMUONTriggerDecision::TrigY(Int_t y1[16], Int_t y2[16], 
				   Int_t y3[16], Int_t y4[16],
				   Int_t y3u[16], Int_t y3d[16], 
				   Int_t y4u[16], Int_t y4d[16],
				   Int_t x2m, Int_t x2ud, Int_t orMud[2], 
				   Int_t resetMid, Int_t coinc44, 
				   Int_t coordY[5]){
// note : resMid = 1 -> cancel 
//---------------------------------------------------------
// step # 1 : prehandling Y
//--------------------------------------------------------- 
    Int_t i;
    Int_t istrip;

  for (i=0; i<16; i++){
    y3[i]=y3[i]&!resetMid;
    y4[i]=y4[i]&!resetMid;
  }

  Int_t ch1[16], ch2[16], ch3[16], ch4[16];

  Int_t tmpy3to16[16], tmpy4to16[16];
  Int_t tmpy3uto16[16], tmpy3dto16[16], tmpy4uto16[16], tmpy4dto16[16];
  for (i=0; i<8; i++){
    ch1[2*i]   = y1[i]&x2m | y1[2*i]&!x2m;		
    ch1[2*i+1] = y1[i]&x2m | y1[2*i+1]&!x2m;

    ch2[2*i]   = y2[i]&x2m | y2[2*i]&!x2m;		
    ch2[2*i+1] = y2[i]&x2m | y2[2*i+1]&!x2m;

    tmpy3to16[2*i]   = y3[i]&x2m | y3[2*i]&!x2m;		
    tmpy3to16[2*i+1] = y3[i]&x2m | y3[2*i+1]&!x2m;

    tmpy4to16[2*i]   = y4[i]&x2m | y4[2*i]&!x2m;
    tmpy4to16[2*i+1] = y4[i]&x2m | y4[2*i+1]&!x2m;

    tmpy3uto16[2*i]   = y3u[i]&x2ud | y3u[2*i]&!x2ud; 
    tmpy3uto16[2*i+1] = y3u[i]&x2ud | y3u[2*i+1]&!x2ud;

    tmpy4uto16[2*i]   = y4u[i]&x2ud | y4u[2*i]&!x2ud; 
    tmpy4uto16[2*i+1] = y4u[i]&x2ud | y4u[2*i+1]&!x2ud;

    tmpy3dto16[2*i]   = y3d[i]&x2ud | y3d[2*i]&!x2ud; 
    tmpy3dto16[2*i+1] = y3d[i]&x2ud | y3d[2*i+1]&!x2ud;
    
    tmpy4dto16[2*i]   = y4d[i]&x2ud | y4d[2*i]&!x2ud; 
    tmpy4dto16[2*i+1] = y4d[i]&x2ud | y4d[2*i+1]&!x2ud;
  }
  
  if (orMud[0]==0&&orMud[1]==0){
    for (i=0; i<16; i++){
      ch3[i] = tmpy3to16[i];
      ch4[i] = tmpy4to16[i];
    }
  }
  if (orMud[0]==0&&orMud[1]==1){
      for (i=0; i<16; i++){
	ch3[i] = tmpy3uto16[i]|tmpy3to16[i];
	ch4[i] = tmpy4uto16[i]|tmpy4to16[i];
      }
  }
  if (orMud[0]==1&&orMud[1]==0){
      for (i=0; i<16; i++){
	ch3[i] = tmpy3dto16[i]|tmpy3to16[i];
	ch4[i] = tmpy4dto16[i]|tmpy4to16[i];
      }
  }
  if (orMud[0]==1&&orMud[1]==1){
      for (i=0; i<16; i++){
	ch3[i] = tmpy3dto16[i]|tmpy3to16[i]|tmpy3uto16[i];
	ch4[i] = tmpy4dto16[i]|tmpy4to16[i]|tmpy4uto16[i];
      }
  }

// debug
  if(fDebug==4||fDebug==5) {
    printf("===============================================================\n");  
    printf(" Y plane after PreHandling x2m x2ud orMud %i %i %i %i \n",
	   x2m,x2ud,orMud[0],orMud[1]);
    printf("                            ");
    for (istrip=15; istrip>=0; istrip--) {
      if (istrip>9)  printf("%i",istrip-10*Int_t(istrip/10));
      if (istrip<10) printf("%i",istrip);
    }  
    printf("\n YMC11                      ");
    for (istrip=15; istrip>=0; istrip--) printf("%i",ch1[istrip]); 
    printf("\n YMC12                      ");
    for (istrip=15; istrip>=0; istrip--) printf("%i",ch2[istrip]); 
    printf("\n YMC21                      ");
    for (istrip=15; istrip>=0; istrip--) printf("%i",ch3[istrip]); 
    printf("\n YMC22                      ");
    for (istrip=15; istrip>=0; istrip--) printf("%i",ch4[istrip]); 
    printf(" \n"); 
  }
//debug
  
//---------------------------------------------------------
// step # 2 : calculate sgle and dble, apply DS reduction
//--------------------------------------------------------- 
  Int_t sgle1[16], dble1[16];
  Int_t sgle2[16], dble2[16];

  // Calculate simple and double hits
  for (i=0; i<16; i++) {
    dble1[i] = ch1[i] & ch2[i];
    dble2[i] = ch3[i] & ch4[i];
    
    sgle1[i] = (ch1[i]|ch2[i]);
    sgle2[i] = (ch3[i]|ch4[i]);
  }

  //debug
  if(fDebug==4||fDebug==5) {
    printf("===============================================================\n");
    printf(" Y plane after sgle dble \n"); 
    printf("                            ");
    for (istrip=15; istrip>=0; istrip--) {
      if (istrip>9)  printf("%i",istrip-10*Int_t(istrip/10));
      if (istrip<10) printf("%i",istrip);
    }  
    printf("\n SGLE1                      ");
    for (istrip=15; istrip>=0; istrip--) printf("%i",sgle1[istrip]); 
    printf("\n DBLE1                      ");
    for (istrip=15; istrip>=0; istrip--) printf("%i",dble1[istrip]); 
    printf("\n SGLE2                      ");
    for (istrip=15; istrip>=0; istrip--) printf("%i",sgle2[istrip]); 
    printf("\n DBLE2                      ");
    for (istrip=15; istrip>=0; istrip--) printf("%i",dble2[istrip]); 
    printf(" \n"); 
  }
  //debug

  // DS Reduction 
  Int_t notOr1, notOr2;

  notOr1=!dble1[15] & !dble1[14] & !dble1[13] & !dble1[12] & 
	 !dble1[11] & !dble1[10] & !dble1[9]  & !dble1[8]  & 
	 !dble1[7]  & !dble1[6]  & !dble1[5]  & !dble1[4]  & 
	 !dble1[3]  & !dble1[2]  & !dble1[1]  & !dble1[0];

  notOr2=!dble2[15] & !dble2[14] & !dble2[13] & !dble2[12] & 
	 !dble2[11] & !dble2[10] & !dble2[9]  & !dble2[8]  & 
	 !dble2[7]  & !dble2[6]  & !dble2[5]  & !dble2[4]  & 
	 !dble2[3]  & !dble2[2]  & !dble2[1]  & !dble2[0];

  for (i=0; i<16; i++) {
    sgle1[i] = sgle1[i] & notOr1 & !coinc44;
    sgle2[i] = sgle2[i] & notOr2 & !coinc44;
  }

//---------------------------------------------------------
// step # 3 : 3/4 coincidence 
//--------------------------------------------------------- 
  Int_t frontImage[16];

  for (i=1; i<15; i++) {
  frontImage[i] = (dble1[i] | sgle1[i]) & 
    (dble2[i+1] | dble2[i] | dble2[i-1]) |
     dble1[i] & (sgle2[i+1] | sgle2[i] | sgle2[i-1]);
  }
  frontImage[0] = (dble1[0] | sgle1[0]) & 
    (dble2[1] | dble2[0]) | dble1[0] & (sgle2[1] | sgle2[0]);

  frontImage[15] = (dble1[15] | sgle1[15]) & 
    (dble2[15] | dble2[14]) | dble1[15] & (sgle2[15] | sgle2[14]);


//debug
  if(fDebug==4||fDebug==5) {
    printf("===============================================================\n");
    printf(" Y plane frontImage\n");
    printf("                            ");
  for (istrip=15; istrip>=0; istrip--) {
    if (istrip>9)  printf("%i",istrip-10*Int_t(istrip/10));
    if (istrip<10) printf("%i",istrip);
  }
  printf("\n                            ");
  for (istrip=15; istrip>=0; istrip--) printf("%i",frontImage[istrip]); 
  printf("\n");
  }
//debug

//---------------------------------------------------------
// step # 4 : Y position 
//--------------------------------------------------------- 
  Int_t or1, or2, and1, and2, and3;

 or1  = frontImage[7]|frontImage[5]|frontImage[3]|frontImage[1];
 or2  = frontImage[7]|frontImage[6]|frontImage[5]|frontImage[4];
 and1 = !frontImage[3]&!frontImage[2]&!frontImage[1]&!frontImage[0];
 and2 = !frontImage[7]&!frontImage[6]&!frontImage[5]&!frontImage[4] & and1;
 and3 = !frontImage[11]&!frontImage[10]&!frontImage[9]&!frontImage[8]; 
 
 coordY[0] = !frontImage[0]&(frontImage[1]|!frontImage[2]) & 
(frontImage[3]|frontImage[1]|!frontImage[4]&(frontImage[5]|!frontImage[6])) &
(or1|!frontImage[8]&(frontImage[9]|!frontImage[10])) & 
(or1|frontImage[11]|frontImage[9]|!frontImage[12]&(frontImage[13]|!frontImage[14]));
 
 coordY[1] = !frontImage[0]&!frontImage[1] & 
!(!frontImage[11]&!frontImage[10]&!frontImage[7]&!frontImage[6] & 
  !frontImage[3]&!frontImage[2]&(frontImage[13]|frontImage[12])) &
  (frontImage[3]|frontImage[2] | !frontImage[5]&!frontImage[4]) & 
  (frontImage[7]|frontImage[6]|frontImage[3]|frontImage[2] | 
!frontImage[9]&!frontImage[8]);
		
 coordY[2] = and1 & (or2 | and3);
		
 coordY[3] = and2;
		
 coordY[4] = !frontImage[15]&!frontImage[14]&!frontImage[13]&!frontImage[12] &
 and2 & and3 ;

}
//----------------------------------------------------------------------
// end of trigger Algo
//----------------------------------------------------------------------

//----------------------------------------------------------------------
void AliMUONTriggerDecision::LocalTrigger(Int_t icirc, 
					  Int_t minDevStrip[5], 
					  Int_t minDev[5], Int_t coordY[5], 
					  Int_t &iTrigger){
// returns local trigger answer for circuit icirc
  Int_t i;

  AliMUONTriggerCircuit* triggerCircuit;
  triggerCircuit = (AliMUONTriggerCircuit*) fTriggerCircuit->At(icirc);  	  
  Int_t idCircuit=triggerCircuit->GetIdCircuit();
  
  Int_t signDev=minDev[4];   
  Int_t deviation=0;
  for (i=0; i<4; i++) {          // extract deviation
    deviation = deviation+Int_t(minDev[i]*TMath::Power(2,i));   
  }
  
  Int_t istripX1Circ=0;
  for (i=0; i<5; i++) {          // extract X1 strip fired 
    istripX1Circ = istripX1Circ+Int_t(minDevStrip[i]*TMath::Power(2,i));   
  }
  
  Int_t iStripY=0;
  for (i=0; i<4; i++) {          // extract Y strip fired 
    iStripY = iStripY+Int_t(coordY[i]*TMath::Power(2,i));   
  }

// trigger or not 
  if (signDev==1&&deviation==0) {      // something in X ?
    iTrigger=0;    
  } else {
    if (coordY[4]==1&&iStripY==15) {   // something in Y ?
      iTrigger=0;
    } else {
      iTrigger=1;
    }
  }
  
  if (iTrigger==1) { 
// fill fTrigger fStripX11 fStripY11 
    fTrigger[icirc] = 1;
    fStripX11[icirc] = istripX1Circ;
    fStripY11[icirc] = iStripY;
    
// calculate deviation in [0+30]
    Int_t sign=0;
    if (signDev==0&&deviation!=0) sign=-1;
    if (signDev==0&&deviation==0) sign=0;
    if (signDev==1)               sign=1;    
    fDev[icirc] = sign * deviation + 15; // fill fDev 

// get Lut output for circuit/istripX/idev/istripY
    AliMUONTriggerLut* lut = new AliMUONTriggerLut;    
    //    lut->StartEvent();
    lut->GetLutOutput(icirc,fStripX11[icirc],fDev[icirc],fStripY11[icirc],
		      fLutLpt[icirc],fLutHpt[icirc],fLutApt[icirc]);
    //    lut->FinishEvent();
    delete lut;
    
    if (fDebug>1) {
      Float_t pt= // get ptCal corresponding to istripX1Circ/idev/iStripY
      triggerCircuit->PtCal(fStripX11[icirc],fDev[icirc],fStripY11[icirc]);
      printf("-------------------------------------------\n");
      printf(" Local Trigger info for circuit Id %i (number %i ) \n",
	     idCircuit,icirc);
      printf(" istripX1 signDev deviation istripY = %i %i %i %i \n", 
	     istripX1Circ,signDev,deviation,iStripY);      
      printf(" pt = %f  (GeV/c) \n",pt);
      printf("-------------------------------------------\n");
      printf(" Local Trigger Lut Output = Lpt : ");
      for (i=1; i>=0; i--) printf("%i",fLutLpt[icirc][i]);
      printf(" Hpt : ");
      for (i=1; i>=0; i--) printf("%i",fLutHpt[icirc][i]);
      printf(" Apt : ");
      for (i=1; i>=0; i--) printf("%i",fLutApt[icirc][i]);
      printf("\n");
      printf("-------------------------------------------\n");
    } // fDebug > 1    
  }  // local trigger = 1
}

//----------------------------------------------------------------------
void AliMUONTriggerDecision::GlobalTrigger(){
// loop on Lut[icirc] and give Global Trigger output
    Int_t i;

  for (Int_t icirc=0; icirc<234; icirc++){
    if (fLutLpt[icirc][0]==1&&fLutLpt[icirc][1]==1) 
      fGlobalSingleUndef[0] = fGlobalSingleUndef[0] + 1;
    if (fLutHpt[icirc][0]==1&&fLutHpt[icirc][1]==1) 
      fGlobalSingleUndef[1] = fGlobalSingleUndef[1] + 1;
    if (fLutApt[icirc][0]==1&&fLutApt[icirc][1]==1) 
      fGlobalSingleUndef[2] = fGlobalSingleUndef[2] + 1;
    
    if (fLutLpt[icirc][0]==0&&fLutLpt[icirc][1]==1) 
      fGlobalSinglePlus[0] = fGlobalSinglePlus[0] + 1;
    if (fLutHpt[icirc][0]==0&&fLutHpt[icirc][1]==1) 
      fGlobalSinglePlus[1] = fGlobalSinglePlus[1] + 1;
    if (fLutApt[icirc][0]==0&&fLutApt[icirc][1]==1) 
      fGlobalSinglePlus[2] = fGlobalSinglePlus[2] + 1;

    if (fLutLpt[icirc][0]==1&&fLutLpt[icirc][1]==0) 
      fGlobalSingleMinus[0] = fGlobalSingleMinus[0] + 1;
    if (fLutHpt[icirc][0]==1&&fLutHpt[icirc][1]==0) 
      fGlobalSingleMinus[1] = fGlobalSingleMinus[1] + 1;
    if (fLutApt[icirc][0]==1&&fLutApt[icirc][1]==0) 
      fGlobalSingleMinus[2] = fGlobalSingleMinus[2] + 1;
  }

  // like sign low, high and all pt
  for (i=0; i<3; i++) {
    fGlobalPairLike[i]=fGlobalSingleMinus[i]*(fGlobalSingleMinus[i]-1)/2 + 
      fGlobalSinglePlus[i]*(fGlobalSinglePlus[i]-1)/2 + 
      fGlobalSingleUndef[i]*(fGlobalSingleUndef[i]-1)/2 + 
      fGlobalSingleUndef[i]*fGlobalSinglePlus[i] + 
      fGlobalSingleUndef[i]*fGlobalSingleMinus[i];
  }

  // unlike sign low, high and all pt
  for (i=0; i<3; i++) {
    fGlobalPairUnlike[i]=fGlobalSingleMinus[i]*fGlobalSinglePlus[i] +
      fGlobalSingleUndef[i]*(fGlobalSingleUndef[i]-1)/2 + 
      fGlobalSingleUndef[i]*fGlobalSinglePlus[i] + 
      fGlobalSingleUndef[i]*fGlobalSingleMinus[i]; 
  }
  
  if (fDebug>=1) {
    printf("===================================================\n");
    printf(" Global Trigger output       Low pt  High pt   All\n");
    printf(" number of Single Plus      :\t");
    for (i=0; i<3; i++) printf("%i\t",fGlobalSinglePlus[i]);
    printf("\n");
    printf(" number of Single Minus     :\t");
    for (i=0; i<3; i++) printf("%i\t",fGlobalSingleMinus[i]);
    printf("\n");
    printf(" number of Single Undefined :\t"); 
    for (i=0; i<3; i++) printf("%i\t",fGlobalSingleUndef[i]);
    printf("\n");
    printf(" number of UnlikeSign pair  :\t"); 
    for (i=0; i<3; i++) printf("%i\t",fGlobalPairUnlike[i]);
    printf("\n");
    printf(" number of LikeSign pair    :\t");  
    for (i=0; i<3; i++) printf("%i\t",fGlobalPairLike[i]);
    printf("\n");
    printf("===================================================\n");
    printf("\n");
  }
}

//----------------------------------------------------------------------
void AliMUONTriggerDecision::PrintBitPatXInput(Int_t icirc){
// print bit pattern for X strips

    Int_t istrip;

  printf("-------- TRIGGER INPUT ---------\n");
  printf("===============================================================\n");
  printf("                            5432109876543210");
  printf("\n XMC11                      ");
  for (istrip=15; istrip>=0; istrip--) printf("%i",fXbit11[icirc][istrip]); 
  printf("\n XMC12                      ");
  for (istrip=15; istrip>=0; istrip--) printf("%i",fXbit12[icirc][istrip]);
  printf("\n XMC21              ");
  for (istrip=31; istrip>=0; istrip--) printf("%i",fXbit21[icirc][istrip]); 
  printf("\n XMC22              ");
  for (istrip=31; istrip>=0; istrip--) printf("%i",fXbit22[icirc][istrip]); 
  printf("\n                    ");
  printf("10987654321098765432109876543210\n");
}

//----------------------------------------------------------------------
void AliMUONTriggerDecision::PrintBitPatYInput(Int_t icirc){
// print bit pattern for Y strips

    Int_t istrip;

  AliMUONTriggerCircuit* triggerCircuit;
  triggerCircuit = (AliMUONTriggerCircuit*) fTriggerCircuit->At(icirc);  	  
  Int_t idCircuit=triggerCircuit->GetIdCircuit();
  Int_t nStrip=triggerCircuit->GetNstripY();

  printf("---------------------------------------------------------------\n");
  printf("                            ");
  for (istrip=nStrip-1; istrip>=0; istrip--) {
    if (istrip>9)  printf("%i",istrip-10*Int_t(istrip/10));
    if (istrip<10) printf("%i",istrip);
  }
  printf("\n YMC11                      ");
  for (istrip=nStrip-1; istrip>=0; istrip--) 
    printf("%i",fYbit11[icirc][istrip]); 
  printf("\n YMC12                      ");
  for (istrip=nStrip-1; istrip>=0; istrip--)
    printf("%i",fYbit12[icirc][istrip]); 
  printf("\n YMC21                      ");
  for (istrip=nStrip-1; istrip>=0; istrip--)
    printf("%i",fYbit21[icirc][istrip]); 
  printf("\n YMC22                      ");
  for (istrip=nStrip-1; istrip>=0; istrip--)
    printf("%i",fYbit22[icirc][istrip]); 
  printf("\n");
// tmp
  printf("---------------------------------------------------------------");
  printf("\n upper part of circuit %i",idCircuit);
  printf("\n UMC21                      ");
  for (istrip=15; istrip>=0; istrip--) printf("%i",fYbit21U[icirc][istrip]); 
  printf("\n UMC22                      ");
  for (istrip=15; istrip>=0; istrip--) printf("%i", fYbit22U[icirc][istrip]); 

  printf("\n lower part of circuit %i",idCircuit);
  printf("\n LMC21                      ");
  for (istrip=15; istrip>=0; istrip--) printf("%i",fYbit21D[icirc][istrip]);
  printf("\n LMC22                      ");
  for (istrip=15; istrip>=0; istrip--) printf("%i",fYbit22D[icirc][istrip]); 
  printf("\n");
  printf("===============================================================\n");
}

//----------------------------------------------------------------------
void AliMUONTriggerDecision::PrintLocalOutput(Int_t minDevStrip[5], 
					      Int_t minDev[5], 
					      Int_t coordY[5]){
// print Local trigger output before the LuT step

    Int_t i;

  printf("===============================================================\n");
  printf("-------- TRIGGER OUTPUT --------\n");
  printf("minDevStrip = ");
  for  (i=4; i>=0; i--) printf("%i",minDevStrip[i]);
  printf(" minDev = ");
  for  (i=4; i>=0; i--) printf("%i",minDev[i]);
  printf(" coordY = ");
  for  (i=4; i>=0; i--) printf("%i",coordY[i]); 
  printf(" \n");
}

//----------------------------------------------------------------------
//--- methods which return member data related info
//----------------------------------------------------------------------
Int_t AliMUONTriggerDecision::GetITrigger(Int_t icirc) const{
// returns Local Trigger Status
  return fTrigger[icirc];
}
//----------------------------------------------------------------------
Int_t AliMUONTriggerDecision::GetStripX11(Int_t icirc) const{
// returns fStripX11
  return fStripX11[icirc];
}
//----------------------------------------------------------------------
Int_t AliMUONTriggerDecision::GetDev(Int_t icirc) const{
// returns idev
  return fDev[icirc];
}
//----------------------------------------------------------------------
Int_t AliMUONTriggerDecision::GetStripY11(Int_t icirc) const{
// returns fStripY11;
   return fStripY11[icirc];
}
//----------------------------------------------------------------------
void AliMUONTriggerDecision::GetLutOutput(Int_t icirc, Int_t lpt[2], 
					  Int_t hpt[2], Int_t apt[2]) const {
// returns Look up Table output
  for (Int_t i=0; i<2; i++) {
    lpt[i]=fLutLpt[icirc][i];
    hpt[i]=fLutHpt[icirc][i];
    apt[i]=fLutApt[icirc][i];
  }
}
//----------------------------------------------------------------------
void AliMUONTriggerDecision::GetGlobalTrigger(Int_t singlePlus[3], 
					      Int_t singleMinus[3], 
					      Int_t singleUndef[3],
					      Int_t pairUnlike[3], 
					      Int_t pairLike[3]) const {
// returns Global Trigger information (0,1,2 : Lpt,Hpt,Apt)
// should not be used anymore.
  for (Int_t i=0; i<3; i++) { 
    singlePlus[i]  = fGlobalSinglePlus[i];
    singleMinus[i] = fGlobalSingleMinus[i];
    singleUndef[i] = fGlobalSingleUndef[i];
    pairUnlike[i]  = fGlobalPairUnlike[i];
    pairLike[i]    = fGlobalPairLike[i];    
  }
}

//_______________________________________________________________________
void AliMUONTriggerDecision::Digits2Trigger(){
// call the Trigger Algorithm and fill TreeD

  ClearDigitNumbers();

  fMUONData->ResetTrigger();
  Trigger();   
  AliMUONGlobalTrigger* pGloTrig = new AliMUONGlobalTrigger(fGlobalSinglePlus, fGlobalSingleMinus,
						       fGlobalSingleUndef, fGlobalPairUnlike, 
						       fGlobalPairLike);  
  // add a local trigger in the list 
  fMUONData->AddGlobalTrigger(*pGloTrig);
  
  for (Int_t icirc=0; icirc<AliMUONConstants::NTriggerCircuit(); icirc++) { 
    if(GetITrigger(icirc)==1) {
      Int_t localtr[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};      
      Int_t loLpt[2]={0,0}; Int_t loHpt[2]={0,0}; Int_t loApt[2]={0,0};
      GetLutOutput(icirc, loLpt, loHpt, loApt);
      localtr[0] = icirc;
      localtr[1] = GetStripX11(icirc);
      localtr[2] = GetDev(icirc);
      localtr[3] = GetStripY11(icirc);
      for (Int_t i = 0; i < 2; i++) {    // convert the Lut output in 1 digit 
	localtr[4] += Int_t(loLpt[i]*TMath::Power(2,i));
	localtr[5] += Int_t(loHpt[i]*TMath::Power(2,i));
	localtr[6] += Int_t(loApt[i]*TMath::Power(2,i));
      }

      for (Int_t i = 0; i < 16; i++) {    // convert X/Y bit in bit pattern
	localtr[7]  |= (fXbit11[icirc][i] << i);
	localtr[8]  |= (fXbit12[icirc][i] << i);

	// 8 first and last elts correspond to neighbouring cards
	localtr[9]  |= (fXbit21[icirc][i+8] << i);
	localtr[10] |= (fXbit22[icirc][i+8] << i);

	localtr[11] |= (fYbit11[icirc][i] << i);
	localtr[12] |= (fYbit12[icirc][i] << i);
	localtr[13] |= (fYbit21[icirc][i] << i);
	localtr[14] |= (fYbit22[icirc][i] << i);
      }

      AliMUONLocalTrigger* pLocTrig = new AliMUONLocalTrigger(localtr, fDigitNumbers[icirc]);
      fMUONData->AddLocalTrigger(*pLocTrig);  // add a local trigger in the list
    }
  }
}

//_______________________________________________________________________
void AliMUONTriggerDecision::ClearDigitNumbers()
{
// Clears the fDigitNumbers arrays so that they are all empty.

	for (Int_t i = 0; i < AliMUONConstants::NTriggerCircuit(); i++)
		fDigitNumbers[i].Set(0);
}

//_______________________________________________________________________
void AliMUONTriggerDecision::DigitFiredCircuit(
		Int_t circuit, Int_t cathode,
		Int_t chamber, Int_t digit
	)
{
// Registers that the specified digit fired the specified circuit.
// This digit gets added to an array which will be copied to
// AliMUONLocalTrigger when such an object is created for each circuit.

	Int_t digitnumber = AliMUONLocalTrigger::EncodeDigitNumber(chamber, cathode, digit);
	Int_t last = fDigitNumbers[circuit].GetSize();
	fDigitNumbers[circuit].Set(last + 1);
	fDigitNumbers[circuit][last] = digitnumber;
}

