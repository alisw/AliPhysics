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
Revision 1.8  2000/07/10 13:58:01  fca
New version of ZDC from E.Scomparin & C.Oppedisano

Revision 1.7  2000/01/19 17:17:40  fca

Revision 1.6  1999/09/29 09:24:35  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Zero Degree Calorimeter                                                  //
//  This class contains the basic functions for the Time Of Flight           //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TBRIK.h>
#include <TNode.h>
#include "TGeometry.h"

#include "AliZDC.h"
#include "AliRun.h"
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
  fHits = 0;
}
 
//_____________________________________________________________________________
AliZDC::AliZDC(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Standard constructor for the Zero Degree Calorimeter base class
  //

  //
  // Allocate the array of hits
  fHits   = new TClonesArray("AliZDChit",  405);
  gAlice->AddHitList(fHits);
  
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
}
 
//_____________________________________________________________________________
void AliZDC::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a ZDC hit
  //
  static Float_t primKinEn, xImpact, yImpact, sFlag;

  TClonesArray &lhits = *fHits;

  AliZDChit *newquad, *curquad;
  newquad = new AliZDChit(fIshunt, track, vol, hits);
  Int_t i;
  for(i=0; i<fNhits; i++){
    // If the hits are equal (same track, same volume), sum them.
     curquad=(AliZDChit*) lhits[i];
     if(*curquad==*newquad){
        *curquad = *curquad+*newquad;
        delete newquad;
//        fHits->Print("");
	return;
      }
   }
   
   //Otherwise create a new hit.
   if(fNhits==0){
      // First hit -> setting flag for primary or secondary particle
      Int_t primary = gAlice->GetPrimary(track);     
      if(track != primary){
        newquad->fSFlag = 1;  // Hit created by secondary particle entering the ZDC
      }
      else if(track == primary){
        newquad->fSFlag = 0;  // Hit created by PRIMARY particle entering the ZDC
      }  
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
   new(lhits[fNhits++]) AliZDChit(newquad);
    if(fNhits==1) {
//      Int_t Curtrack = gAlice->CurrentTrack();
//      Int_t Prim = gAlice->GetPrimary(Curtrack);
//      printf ("		Primary track: %d, Current track: %d \n", 
//              Prim, Curtrack);
//      fHits->Print("");
    }
    delete newquad;
//    fHits->Print("");
  
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

ClassImp(AliZDChit)
  
//_____________________________________________________________________________
AliZDChit::AliZDChit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
  //
  // Add a ZDC hit
  //
  Int_t i;
  for(i=0; i<2; i++) fVolume[i] = vol[i];
  fX = hits[0];
  fY = hits[1];
  fZ = hits[2];
  fPrimKinEn = hits[3];
  fXImpact = hits[4];
  fYImpact = hits[5];
  fSFlag = hits[6];
  fLightPMQ = hits[7];
  fLightPMC = hits[8];
  fEnergy = hits[9]; 
}
