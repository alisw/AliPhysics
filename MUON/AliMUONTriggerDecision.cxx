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
/*
$Log$
Revision 1.3  2000/06/25 17:02:19  pcrochet
scope problem on HP, i declared once, pow replaced by TMath::Power (PH)

Revision 1.2  2000/06/15 07:58:49  morsch
Code from MUON-dev joined

Revision 1.1.2.8  2000/06/14 14:54:34  morsch
Complete redesign, make use of TriggerCircuit and TriggerLut (PC)

Revision 1.1.2.5  2000/04/26 19:59:57  morsch
Constructor added.

Revision 1.1.2.4  2000/04/26 12:31:30  morsch
Modifications by P. Crochet:
- adapted to the new Trigger chamber geometry
- condition on soft background added
- contructor added in AliMUONTriggerDecision.h
- single-undefined taken into account in the output of GlobalTrigger()
- some bugs fixed

Revision 1.1.2.3  2000/03/21 09:29:58  morsch
Put back comments

Revision 1.1.2.2  2000/03/21 09:24:34  morsch
Author and responsible for the code: Philippe Crochet
*/


#include "AliMUONTriggerDecision.h"
#include "AliMUONTriggerLut.h"
#include "AliMUONHitMapA1.h"
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONPoints.h"
#include "AliSegmentation.h"
#include "AliMUONResponse.h"
#include "AliMUONChamber.h"
#include "AliMUONDigit.h"


#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TPad.h>
#include <TGraph.h> 
#include <TPostScript.h> 
#include <TMinuit.h> 
#include <iostream.h> 

//----------------------------------------------------------------------
ClassImp(AliMUONTriggerDecision)

//----------------------------------------------------------------------
AliMUONTriggerDecision::AliMUONTriggerDecision(Int_t iprint)
{
// Constructor 
  fiDebug = iprint;            // print option
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
    fiTrigger[icirc]=0;                   // trigger or not
    fStripX11[icirc]=0;                   // X strip in MC11 which triggers 
    fdev[icirc]=0;                        // deviation which triggers 
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
}

//----------------------------------------------------------------------
AliMUONTriggerDecision::~AliMUONTriggerDecision()
{
// Destructor
}

//----------------------------------------------------------------------
void AliMUONTriggerDecision::Trigger(){
// main method of the class which calls the overall Trigger procedure
//  cout << " In AliMUONTriggerDecision::Trigger " << "\n";

  ResetBit();
  SetBit();
  SetBitUpDownY();

  Int_t coinc44=0, resetMid=0; // initialize coincidence

  AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");  
  AliMUONTriggerCircuit* triggerCircuit;

  for (Int_t icirc=0; icirc<234; icirc++) {  // loop on circuits
    triggerCircuit = &(pMUON->TriggerCircuit(icirc));  	  
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

    if (iTrigger==1&&fiDebug>1) { 
      PrintBitPatXInput(icirc);
      PrintBitPatYInput(icirc);
      PrintLocalOutput(minDevStrip, minDev, coordY);
    }      
  }  //  end loop on circuits

// call Global Trigger
  GlobalTrigger();
  //  cout << " Leaving AliMUONTriggerDecision::Trigger " << "\n";
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
    fiTrigger[icirc]=0;
    fStripX11[icirc]=0;
    fdev[icirc]=0;                      
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

  AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");  
  AliMUONTriggerCircuit* triggerCircuit;

  for (Int_t chamber=11; chamber<15; chamber++){
    for (Int_t cathode=1; cathode<3; cathode++){
      
      AliMUONChamber*   iChamber;
      AliSegmentation*  segmentation;
      
      TClonesArray *muonDigits  = pMUON->DigitsAddress(chamber-1);
      if (muonDigits == 0) return;
      
      gAlice->ResetDigits();
      
      Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
      gAlice->TreeD()->GetEvent(nent-2+cathode-1);
      Int_t ndigits = muonDigits->GetEntriesFast();
      if (ndigits == 0) return;
      
      iChamber = &(pMUON->Chamber(chamber-1));
      segmentation=iChamber->SegmentationModel(cathode);
      AliMUONDigit  *mdig;
      
      for (Int_t digit=0; digit<ndigits; digit++) {
	mdig    = (AliMUONDigit*)muonDigits->UncheckedAt(digit);
// get the center of the pad Id 
  	Int_t ix=mdig->fPadX;
  	Int_t iy=mdig->fPadY;
// get the sum of the coded charge 
// see coding convention in AliMUONChamberTrigger::DisIntegration 	
	Int_t sumCharge=0;
	for (Int_t icharge=0; icharge<10; icharge++) {
	  sumCharge=sumCharge+mdig->fTcharges[icharge];
	}
// apply condition on soft background	
	Int_t testCharge=sumCharge-(Int_t(sumCharge/10))*10;
	testCharge=sumCharge-testCharge*10;
	if(sumCharge<=10||testCharge>0) {	  
// code pad
	  Int_t code=TMath::Abs(ix)*100+iy;
	  if (ix<0) { code=-code; }
	  
	  Int_t icirc;
	  Int_t istrip;
	  Int_t nStrip;

	  if (cathode==1) {
	    switch (chamber)
	      {
	      case 11:
		for (icirc=0; icirc<234; icirc++) {		  
		  triggerCircuit = &(pMUON->TriggerCircuit(icirc));  	  
		  for (istrip=0; istrip<16; istrip++) {
		    if (triggerCircuit->GetXcode(0,istrip)==code) 
		      fXbit11[icirc][istrip]=1;
		  }
		}
		break;
	      case 12:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = &(pMUON->TriggerCircuit(icirc));  	  
		  for (istrip=0; istrip<16; istrip++) {
		    if (triggerCircuit->GetXcode(1,istrip)==code) 
		      fXbit12[icirc][istrip]=1;
		  }
		}
		break;
	      case 13:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = &(pMUON->TriggerCircuit(icirc));  	  
		  for (istrip=0; istrip<32; istrip++) {
		    if (triggerCircuit->GetXcode(2,istrip)==code) 
		      fXbit21[icirc][istrip]=1;
		  }
		}
		break;
	      case 14:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = &(pMUON->TriggerCircuit(icirc));  	  
		  for (istrip=0; istrip<32; istrip++) {
		    if (triggerCircuit->GetXcode(3,istrip)==code) 
		      fXbit22[icirc][istrip]=1;		    
		  }
		}		
		break;
	      }
	    
	  } else {                // Y plane 
	    switch (chamber)
	      {
	      case 11:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = &(pMUON->TriggerCircuit(icirc));  	  
		  nStrip=triggerCircuit->GetNstripY();
		  for (istrip=0; istrip<nStrip; istrip++) {
		    if (triggerCircuit->GetYcode(0,istrip)==code) 
		      fYbit11[icirc][istrip]=1;
		  }
		}
		break;
	      case 12:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = &(pMUON->TriggerCircuit(icirc));  	  
		  nStrip=triggerCircuit->GetNstripY(); 
		  for (istrip=0; istrip<nStrip; istrip++) {
		    if (triggerCircuit->GetYcode(1,istrip)==code) 
		      fYbit12[icirc][istrip]=1;
		  }
		}
		break;
	      case 13:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = &(pMUON->TriggerCircuit(icirc));  	  
		  nStrip=triggerCircuit->GetNstripY();    
		  for (istrip=0; istrip<nStrip; istrip++) {
		    if (triggerCircuit->GetYcode(2,istrip)==code) 
		      fYbit21[icirc][istrip]=1;
		  }
		}
		break;
	      case 14:
		for (icirc=0; icirc<234; icirc++) {
		  triggerCircuit = &(pMUON->TriggerCircuit(icirc));  	  
		  nStrip=triggerCircuit->GetNstripY();    
		  for (istrip=0; istrip<nStrip; istrip++) {
		    if (triggerCircuit->GetYcode(3,istrip)==code) 
		      fYbit22[icirc][istrip]=1;		      		 
		  }
		}		
		break;
	      }
	  } // if cathode
	}  // remove soft background
      }   // end loop on digit
    }    // end loop on cathode
  }     // end loop on chamber
}  

//----------------------------------------------------------------------
void AliMUONTriggerDecision::SetBitUpDownY(){
// Set Y bit for up and down parts of circuits
  Int_t idModule, nStripX, nStripY, iPosCircuit;

  AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");
  
  for (Int_t icirc=0; icirc<234; icirc++) {

    AliMUONTriggerCircuit* circuit;   // current circuit
    AliMUONTriggerCircuit* circuitD;  // circuit Down
    AliMUONTriggerCircuit* circuitU;  // circuit Up

    circuit = &(pMUON->TriggerCircuit(icirc));  
    idModule=circuit->GetIdModule();      // corresponding module Id.
    nStripX=circuit->GetNstripX();        // number of X strips
    nStripY=circuit->GetNstripY();        // number of Y strips
    iPosCircuit=circuit->GetPosCircuit(); // position of circuit in module

// fill lower part
    if (iPosCircuit==1) {               // need to scan lower module       
      if(idModule<91&&TMath::Abs(idModule)!=41&&idModule>-91) { 
	Int_t icircD=circuit->GetICircuitD();
	circuitD = &(pMUON->TriggerCircuit(icircD));  
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
	circuitU = &(pMUON->TriggerCircuit(icircU));  
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
  if(fiDebug==3||fiDebug==5) {
    cout << "===============================================================" << "\n";
    cout << " X plane after sgle and dble " << " \n";
    cout << "                       0987654321098765432109876543210";
    cout << "\n SGLE1                 ";
    for (istrip=30; istrip>=0; istrip--) { cout << (!sgleHit1[istrip]); }
    cout << "\n DBLE1                 ";
    for (istrip=30; istrip>=0; istrip--) { cout << dbleHit1[istrip]; }
    cout << "\n SGLE2 ";
    for (istrip=62; istrip>=0; istrip--) { cout << (!sgleHit2[istrip]); }
    cout << "\n DBLE2 ";
    for (istrip=62; istrip>=0; istrip--) { cout << dbleHit2[istrip]; }
    cout << "\n       210987654321098765432109876543210987654321098765432109876543210" << "\n";
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
 if(fiDebug==3||fiDebug==5) {
    cout << "===============================================================" << "\n";
   for (i=30; i>=0; i--) {
     cout << i << "\t ";
     for (istrip=31; istrip>=0; istrip--) {
       cout << rearImage[i][istrip];
     }
     cout << " " << "\n";
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
  if(fiDebug==3||fiDebug==5) {
    cout << "===============================================================" << "\n";
    for (i=30; i>=0; i--) {
      cout << i << "\t ";
      for (istrip=5; istrip>=0; istrip--) { cout << dev[i][istrip]; }
      cout << " " << "\n";
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
  if(fiDebug==3||fiDebug==5) {
    cout << "===============================================================" << "\n";
    cout << " sorting : 1st level " << "\n";
    for (i=15; i>=0; i--) {
      cout << i << "\t " << bga1[i] << "\t"; 	
      for (j=5; j>=0; j--) {
	cout << tmpbga1[i][j]; 
      }
      cout << " " << "\n";
    }
  }

  for (i=0; i<8; i++) {  
    Sort2x5(tmpbga1[2*i],tmpbga1[2*i+1],tmpbga2[i],bga2[i]);
  }

//--    
  if(fiDebug==3||fiDebug==5) {
    cout << "===============================================================" << "\n";
    cout << " sorting : 2nd level " << "\n";
    for (i=7; i>=0; i--) {
      cout << i << "\t " << bga2[i] << "\t"; 	
      for (j=5; j>=0; j--) {
	cout << tmpbga2[i][j]; 
      }
      cout << " " << "\n";
    }
  }
    
  for (i=0; i<4; i++) {  
    Sort2x5(tmpbga2[2*i],tmpbga2[2*i+1],tmpbga3[i],bga3[i]);
  }

//--    
  if(fiDebug==3||fiDebug==5) {
    cout << "===============================================================" << "\n";
    cout << " sorting : 3rd level " << "\n";
    for (i=3; i>=0; i--) {
      cout << i << "\t " << bga3[i] << "\t"; 	
      for (j=5; j>=0; j--) {
	cout << tmpbga3[i][j]; 
      }
      cout << " " << "\n";
    }
  }

  for (i=0; i<2; i++) {  
    Sort2x5(tmpbga3[2*i],tmpbga3[2*i+1],tmpbga4[i],bga4[i]);
  }

//--    
  if(fiDebug==3||fiDebug==5) {
    cout << "===============================================================" << "\n";
    cout << " sorting : 4th level " << "\n";
    for (i=1; i>=0; i--) {
      cout << i << "\t " << bga4[i] << "\t"; 	
      for (j=5; j>=0; j--) {
	cout << tmpbga4[i][j]; 
      }
      cout << " " << "\n";
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

  if(fiDebug==3||fiDebug==5) {
    cout << "===============================================================" << "\n";
    cout << "minDevStrip = ";
    for  (i=4; i>=0; i--) {cout << minDevStrip[i];}
    cout << " minDev = ";
    for  (i=4; i>=0; i--) {cout << minDev[i];} 
    cout << " " << "\n";
    cout << "===============================================================" << "\n";
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
  if(fiDebug==4||fiDebug==5) {
    cout << "===============================================================" << "\n";  
    cout << " Y plane after PreHandling x2m x2ud orMud " 
	 << x2m << " , " << x2ud << " , " << orMud[0] << orMud[1] << "\n"; 
    cout << "                            ";
    for (istrip=15; istrip>=0; istrip--) {
      if (istrip>9) cout << istrip-10*Int_t(istrip/10);
      if (istrip<10) cout << istrip;
    }  
    cout << "\n YMC11                      ";
    for (istrip=15; istrip>=0; istrip--) {
      cout << ch1[istrip]; 
    }
    cout << "\n YMC12                      ";
    for (istrip=15; istrip>=0; istrip--) {
      cout << ch2[istrip]; 
    }
    cout << "\n YMC21                      ";
    for (istrip=15; istrip>=0; istrip--) {
      cout << ch3[istrip]; 
    }
    cout << "\n YMC22                      ";
    for (istrip=15; istrip>=0; istrip--) {
      cout << ch4[istrip]; 
    }
    cout << " \n"; 
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
  if(fiDebug==4||fiDebug==5) {
    cout << "===============================================================" << "\n";
    cout << " Y plane after sgle dble " << "\n"; 
    cout << "                            ";
    for (istrip=15; istrip>=0; istrip--) {
      if (istrip>9) { cout << istrip-10*Int_t(istrip/10);}
      if (istrip<10) { cout << istrip;}
    }  
    cout << "\n SGLE1                      ";
    for (istrip=15; istrip>=0; istrip--) {
      cout << sgle1[istrip]; 
    }
    cout << "\n DBLE1                      ";
    for (istrip=15; istrip>=0; istrip--) {
      cout << dble1[istrip]; 
    }
    cout << "\n SGLE2                      ";
    for (istrip=15; istrip>=0; istrip--) {
      cout << sgle2[istrip]; 
    }
    cout << "\n DBLE2                      ";
    for (istrip=15; istrip>=0; istrip--) {
      cout << dble2[istrip]; 
    }
    cout << " \n"; 
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
  if(fiDebug==4||fiDebug==5) {
    cout << "===============================================================" << "\n";
    cout << " Y plane frontImage\n";
    cout << "                            ";
  for (istrip=15; istrip>=0; istrip--) {
    if (istrip>9) cout << istrip-10*Int_t(istrip/10);
    if (istrip<10) cout << istrip;
  }
    cout << "\n                            ";
    for (istrip=15; istrip>=0; istrip--) {
      cout << frontImage[istrip]; 
    }
    cout << "\n";
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

  AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");  
  AliMUONTriggerCircuit* triggerCircuit;
  triggerCircuit = &(pMUON->TriggerCircuit(icirc));  	  
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
// fill fiTrigger fStripX11 fStripY11 
    fiTrigger[icirc] = 1;
    fStripX11[icirc] = istripX1Circ;
    fStripY11[icirc] = iStripY;
    
// calculate deviation in [0+30]
    Int_t sign=0;
    if (signDev==0&&deviation!=0) sign=-1;
    if (signDev==0&&deviation==0) sign=0;
    if (signDev==1)               sign=1;    
    fdev[icirc] = sign * deviation + 15; // fill fdev 

// get Lut output for circuit/istripX/idev/istripY
    AliMUONTriggerLut* lut = new AliMUONTriggerLut;    
    //    lut->StartEvent();
    lut->GetLutOutput(icirc,fStripX11[icirc],fdev[icirc],fStripY11[icirc],
		      fLutLpt[icirc],fLutHpt[icirc],fLutApt[icirc]);
    //    lut->FinishEvent();
    delete lut;
    
    if (fiDebug>1) {
      Float_t pt= // get ptCal corresponding to istripX1Circ/idev/iStripY
      triggerCircuit->PtCal(fStripX11[icirc],fdev[icirc],fStripY11[icirc]);
      cout << "-------------------------------------------" << "\n";
      cout << " Local Trigger info for circuit Id " << idCircuit 
	   << " (number " << icirc << ")" << "\n";
      cout << " istripX1 signDev deviation istripY = " 
	   << istripX1Circ << " , " << signDev 
	   << " , " << deviation << " , " << iStripY << "\n";      
      cout << " pt = " << pt << " (GeV/c) " << "\n";
      cout << "-------------------------------------------" << "\n";
      cout << " Local Trigger Lut Output = Lpt : " ;
      for (i=1; i>=0; i--) { cout << fLutLpt[icirc][i] ; }
      cout << " Hpt : ";
      for (i=1; i>=0; i--) { cout << fLutHpt[icirc][i] ; }
      cout << " Apt : ";
      for (i=1; i>=0; i--) { cout << fLutApt[icirc][i] ; }	  
      cout << "\n";
      cout << "-------------------------------------------" << "\n";
    } // fiDebug > 1    
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
  
  if (fiDebug>=1) {
    cout << "\n";
    cout << "===================================================" << "\n";
    cout << " Global Trigger output       " << "Low pt  High pt   All"  << "\n";
    cout << " number of Single Plus      :\t";
    for (i=0; i<3; i++) { cout << fGlobalSinglePlus[i] <<"\t";}
    cout << "\n";
    cout << " number of Single Minus     :\t";
    for (i=0; i<3; i++) { cout << fGlobalSingleMinus[i] <<"\t";}
    cout << "\n";
    cout << " number of Single Undefined :\t"; 
    for (i=0; i<3; i++) { cout << fGlobalSingleUndef[i] <<"\t";}
    cout << "\n";
    cout << " number of UnlikeSign pair  :\t"; 
    for (i=0; i<3; i++) { cout << fGlobalPairUnlike[i] <<"\t";}
    cout << "\n";
    cout << " number of LikeSign pair    :\t";  
    for (i=0; i<3; i++) { cout << fGlobalPairLike[i] <<"\t";}
    cout << "\n";
    cout << "===================================================" << "\n";
  }
}

//----------------------------------------------------------------------
void AliMUONTriggerDecision::PrintBitPatXInput(Int_t icirc){
// print bit pattern for X strips

    Int_t istrip;

  cout << "-------- TRIGGER INPUT ---------" << "\n";
  cout << "===============================================================" << "\n";
  cout << "                            5432109876543210";
  cout << "\n XMC11                      ";
  for (istrip=15; istrip>=0; istrip--) {
    cout << fXbit11[icirc][istrip]; 
  }
  cout << "\n XMC12                      ";
  for (istrip=15; istrip>=0; istrip--) {
    cout << fXbit12[icirc][istrip]; 
  }
  cout << "\n XMC21              ";
  for (istrip=31; istrip>=0; istrip--) {
    cout << fXbit21[icirc][istrip]; 
  }
  cout << "\n XMC22              ";
  for (istrip=31; istrip>=0; istrip--) {
    cout << fXbit22[icirc][istrip]; 
  }
  cout << "\n                    ";
  cout << "10987654321098765432109876543210" << "\n";
}

//----------------------------------------------------------------------
void AliMUONTriggerDecision::PrintBitPatYInput(Int_t icirc){
// print bit pattern for Y strips

    Int_t istrip;

  AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");  
  AliMUONTriggerCircuit* triggerCircuit;
  triggerCircuit = &(pMUON->TriggerCircuit(icirc));  	  
  Int_t idCircuit=triggerCircuit->GetIdCircuit();
  Int_t nStrip=triggerCircuit->GetNstripY();

  cout << "---------------------------------------------------------------" << "\n";
  cout << "                            ";
  for (istrip=nStrip-1; istrip>=0; istrip--) {
    if (istrip>9) { cout << istrip-10*Int_t(istrip/10);}
    if (istrip<10) { cout << istrip;}
  }
  cout << "\n YMC11                      ";
  for (istrip=nStrip-1; istrip>=0; istrip--) {
    cout << fYbit11[icirc][istrip]; 
  }
  cout << "\n YMC12                      ";
  for (istrip=nStrip-1; istrip>=0; istrip--) {
    cout << fYbit12[icirc][istrip]; 
  }
  cout << "\n YMC21                      ";
  for (istrip=nStrip-1; istrip>=0; istrip--) {
    cout << fYbit21[icirc][istrip]; 
  }
  cout << "\n YMC22                      ";
  for (istrip=nStrip-1; istrip>=0; istrip--) {
    cout << fYbit22[icirc][istrip]; 
  }
  cout << "\n";
// tmp
  cout << "---------------------------------------------------------------";
  cout << "\n upper part of circuit " << idCircuit ;
  cout << "\n UMC21                      ";
  for (istrip=15; istrip>=0; istrip--) {
    cout << fYbit21U[icirc][istrip]; 
  }
  cout << "\n UMC22                      ";
  for (istrip=15; istrip>=0; istrip--) {
    cout << fYbit22U[icirc][istrip]; 
  }

  cout << "\n lower part of circuit " << idCircuit ;
  cout << "\n LMC21                      ";
  for (istrip=15; istrip>=0; istrip--) {
    cout << fYbit21D[icirc][istrip]; 
  }
  cout << "\n LMC22                      ";
  for (istrip=15; istrip>=0; istrip--) {
    cout << fYbit22D[icirc][istrip]; 
  }
  cout << "\n";
  cout << "===============================================================" << "\n";
}
//----------------------------------------------------------------------
void AliMUONTriggerDecision::PrintLocalOutput(Int_t minDevStrip[5], 
					      Int_t minDev[5], 
					      Int_t coordY[5]){
// print Local trigger output before the LuT step

    Int_t i;

  cout << "===============================================================" << "\n";
  cout << "-------- TRIGGER OUTPUT --------" << "\n";
  cout << "minDevStrip = ";
  for  (i=4; i>=0; i--) {cout << minDevStrip[i];}
  cout << " minDev = ";
  for  (i=4; i>=0; i--) {cout << minDev[i];} 
  cout << " coordY = ";
  for  (i=4; i>=0; i--) {cout << coordY[i];} 
  cout << " " << "\n";  
}

//----------------------------------------------------------------------
//--- methods which return member data related info
//----------------------------------------------------------------------
Int_t AliMUONTriggerDecision::GetITrigger(Int_t icirc){
// returns Local Trigger Status
  return fiTrigger[icirc];
}
//----------------------------------------------------------------------
Int_t AliMUONTriggerDecision::GetStripX11(Int_t icirc){
// returns fStripX11
  return fStripX11[icirc];
}
//----------------------------------------------------------------------
Int_t AliMUONTriggerDecision::GetDev(Int_t icirc){
// returns idev
  return fdev[icirc];
}
//----------------------------------------------------------------------
Int_t AliMUONTriggerDecision::GetStripY11(Int_t icirc){
// returns fStripY11;
   return fStripY11[icirc];
}
//----------------------------------------------------------------------
void AliMUONTriggerDecision::GetLutOutput(Int_t icirc, Int_t lpt[2], 
					  Int_t hpt[2], Int_t apt[2]){
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
					      Int_t pairLike[3]){
// returns Global Trigger information (0,1,2 : Lpt,Hpt,Apt)
  for (Int_t i=0; i<3; i++) { 
    singlePlus[i]  = fGlobalSinglePlus[i];
    singleMinus[i] = fGlobalSingleMinus[i];
    singleUndef[i] = fGlobalSingleUndef[i];
    pairUnlike[i]  = fGlobalPairUnlike[i];
    pairLike[i]    = fGlobalPairLike[i];    
  }
}
//----------------------------------------------------------------------
//--- end of methods which return member data related info
//----------------------------------------------------------------------
//----------------------------------------------------------------------
/*
void AliMUONTriggerDecision::AddLocalTrigger(const AliMUONLocalTrigger c){
// Add a Local Trigger copy to the list
  AliMUON *MUON=(AliMUON*)gAlice->GetModule("MUON");
  MUON->AddLocalTrigger(c); 
  fNLocalTriggers++;
}
*/
