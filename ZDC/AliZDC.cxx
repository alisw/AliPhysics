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

/*
$Log$
Revision 1.12  2000/12/01 08:19:01  coppedis
Adding a message error if ZDC is constructed without DIPO

Revision 1.11  2000/11/30 17:21:03  coppedis
Introduce hit array fStHits reset only at the end of the event (for digitization)

Revision 1.10  2000/11/22 11:32:58  coppedis
Major code revision

Revision 1.9  2000/10/02 21:28:20  fca
Removal of useless dependecies via forward declarations

Revision 1.8  2000/07/10 13:58:01  fca
New version of ZDC from E.Scomparin & C.Oppedisano

Revision 1.7  2000/01/19 17:17:40  fca

Revision 1.6  1999/09/29 09:24:35  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Zero Degree Calorimeter                                                  //
//  This class contains the basic functions for the ZDCs                     //
//  Functions specific to one particular geometry are      	             //
//  contained in the derived classes                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

// --- ROOT system
#include <TBRIK.h>
#include <TNode.h>
#include "TGeometry.h"
#include "TTree.h"

// --- AliRoot header files
#include "AliZDC.h"
#include "AliZDCHit.h"
#include "AliRun.h"
#include "AliDetector.h"
#include "AliCallf77.h"
#include "AliConst.h"
#include "AliMC.h"

 
ClassImp(AliZDC)
 
//_____________________________________________________________________________
AliZDC::AliZDC()
{
  //
  // Default constructor for the Zero Degree Calorimeter base class
  //
  
  fIshunt = 1;

  fNhits = 0;
  fNStHits = 0;
  
  fStHits = new TClonesArray("AliZDCHit",1000);
  
  fNPrimaryHits = 0;
}
 
//_____________________________________________________________________________
AliZDC::AliZDC(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Standard constructor for the Zero Degree Calorimeter base class
  //

  // Check that DIPO is there (otherwise tracking is wrong!!!)
  
  AliModule* DIPO=gAlice->GetModule("DIPO");
  if(!DIPO) {
    Error("Constructor","ZDC needs DIPO!!!\n");
    exit(1);
  } 

  //
  // Allocate the array of hits

  fHits   = new TClonesArray("AliZDCHit",1000);
  gAlice->AddHitList(fHits);
  
  fStHits = new TClonesArray("AliZDCHit",1000);

  fNStHits = 0;

  fNPrimaryHits = 0;
  
  fIshunt =  1;
  
  fDimZN[0] = 3.52;
  fDimZN[1] = 3.52;
  fDimZN[2] = 50.;
  fDimZP[0] = 11.2;
  fDimZP[1] = 6.;
  fDimZP[2] = 75.;
  fPosZN[0] = 0.;
  fPosZN[1] = 0.;
  fPosZN[2] = 11650.;
  fPosZP[0] = -23.;
  fPosZP[1] = 0.;
  fPosZP[2] = 11600.;
  fFibZN[0] = 0.;
  fFibZN[1] = 0.01825;
  fFibZN[2] = 50.;
  fFibZP[0] = 0.;
  fFibZP[1] = 0.0275;
  fFibZP[2] = 75.;
  fGrvZN[0] = 0.03;
  fGrvZN[1] = 0.03;
  fGrvZN[2] = 50.;
  fGrvZP[0] = 0.04;
  fGrvZP[1] = 0.04;
  fGrvZP[2] = 75.;
  fDivZN[0] = 11;
  fDivZN[1] = 11;
  fDivZN[2] = 0;
  fDivZP[0] = 7;
  fDivZP[1] = 15;
  fDivZP[2] = 0;
  fTowZN[0] = 2;
  fTowZN[1] = 2;
  fTowZP[0] = 4;
  fTowZP[1] = 1;
  
  // EM Calorimeter
  fDimZEMPb  = 0.15*(TMath::Sqrt(2.));
  fDimZEMAir = 0.001;
  fFibRadZEM = 0.0315;
  fDivZEM[0] = 92;
  fDivZEM[1] = 0;
  fDivZEM[2] = 20;
  fDimZEM[0] = 2*fDivZEM[2]*(fDimZEMPb+fDimZEMAir+fFibRadZEM*(TMath::Sqrt(2.)));
  fDimZEM[1] = 3.5;
  fDimZEM[2] = 3.5;
  fDimZEM[3] = 45.;
  fDimZEM[4] = 0.;
  fDimZEM[5] = 0.;
  fFibZEM[0] = 0.;
  fFibZEM[1] = 0.0275;
  fFibZEM[2] = fDimZEM[2]/TMath::Sin(fDimZEM[3]*kDegrad)-fFibRadZEM;
  fPosZEM[0] = 0.;
  fPosZEM[1] = 5.8;
  fPosZEM[2] = 11600.;

}
//____________________________________________________________________________ 
AliZDC::~AliZDC()
{
  //
  // ZDC destructor
  //

  fIshunt   = 0;
//  delete fHits;
//  if(fStHits){
//    fStHits->Delete();
//    delete fStHits;
//    fNStHits = 0;
//  }
//  delete fDigits;
}
//_____________________________________________________________________________
void AliZDC::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // 		Add a ZDC hit to the hit list.
  // -> We make use of 2 array of hits:
  // [1]  fHits (the usual one) that contains hits for each PRIMARY
  // [2]  fStHits that contains hits for each EVENT and is used to
  //	  obtain digits at the end of each event
  //
  
  static Float_t primKinEn, xImpact, yImpact, sFlag;

  TClonesArray &lsthits = *fStHits;


  AliZDCHit *newquad, *curevquad, *curprimquad;
  newquad = new AliZDCHit(fIshunt, track, vol, hits);

  TClonesArray &lhits = *fHits;
   
  Int_t i,j,kStHit = 1;
  for(i=0; i<fNStHits; i++){
    // If the hits are equal (same track, same volume), sum them.
     curevquad = (AliZDCHit*) lsthits[i];
     kStHit = 1;
     if(*curevquad == *newquad){
        *curevquad = *curevquad+*newquad;
        kStHit = 0;
     } 
  }

  for(j=0; j<fNhits; j++){
    // If the hits are equal (same track, same volume), sum them.
     curprimquad = (AliZDCHit*) lhits[j];
     if(*curprimquad == *newquad){
        *curprimquad = *curprimquad+*newquad;
	delete newquad;
	return;
     } 
  }
  
  if(fNhits==0){
      // First hit -> setting flag for primary or secondary particle
      Int_t primary = gAlice->GetPrimary(track);     
      if(track != primary){
        newquad->fSFlag = 1;  // Hit created by secondary particle entering the ZDC
      }
      else if(track == primary){
        newquad->fSFlag = 0;  // Hit created by PRIMARY particle entering the ZDC
      }  
      fNPrimaryHits = fNPrimaryHits + 1;
      sFlag = newquad->fSFlag;
      primKinEn = newquad->fPrimKinEn;
      xImpact = newquad->fXImpact;
      yImpact = newquad->fYImpact;
   }
   else{       
      newquad->fPrimKinEn = primKinEn;
      newquad->fXImpact = xImpact;
      newquad->fYImpact = yImpact;
      newquad->fSFlag = sFlag;
   }

    //Otherwise create a new hit
    new(lhits[fNhits]) AliZDCHit(newquad);
    fNhits++;
    
    if(kStHit){
      new(lsthits[fNStHits]) AliZDCHit(newquad);
      fNStHits++;
    }
 
//    printf("\n  Primary Hits --------------------------------------------------------\n");
//    fHits->Print("");
//    printf("\n  Event Hits --------------------------------------------------------\n");
//    fStHits->Print("");

    delete newquad;
  }
  
//_____________________________________________________________________________
void AliZDC::ResetDigits()
{
  //
  // Reset number of digits and the digits array
  //
    
    AliDetector::ResetDigits();
//    fNStHits = 0;
//    if(fStHits) fStHits->Clear();
}

//_____________________________________________________________________________
void AliZDC::BuildGeometry()
{
  //
  // Build the ROOT TNode geometry for event display 
  // in the Zero Degree Calorimeter
  // This routine is dummy for the moment
  //

  TNode *node, *top;
  TBRIK *brik;
  const int kColorZDC  = kRed;
  
  //
  top=gAlice->GetGeometry()->GetNode("alice");
  
  // ZDC
    brik = new TBRIK("S_ZDC","ZDC box","void",300,300,5);
    top->cd();
    node = new TNode("ZDC","ZDC","S_ZDC",0,0,600,"");
    node->SetLineColor(kColorZDC);
    fNodes->Add(node);
}

//_____________________________________________________________________________
Int_t AliZDC::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Distance from the mouse to the Zero Degree Calorimeter
  // Dummy routine
  //
  return 9999;
}
 
//_____________________________________________________________________________
void AliZDC::StepManager()
{
  //
  // Routine called at every step in the Zero Degree Calorimeter
  //
}
