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
Revision 1.3  2000/06/28 15:16:35  morsch
(1) Client code adapted to new method signatures in AliMUONSegmentation (see comments there)
to allow development of slat-muon chamber simulation and reconstruction code in the MUON
framework. The changes should have no side effects (mostly dummy arguments).
(2) Hit disintegration uses 3-dim hit coordinates to allow simulation
of chambers with overlapping modules (MakePadHits, Disintegration).

Revision 1.2  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

Revision 1.1.2.3  2000/06/09 21:27:35  morsch
Most coding rule violations corrected.

Revision 1.1.2.2  2000/04/26 12:28:25  morsch
- flag pad hits with condition on ToF (CP)
- Tof included in the method DisIntegration (CP)

Revision 1.1.2.1  2000/02/17 14:30:54  morsch
Draft version

*/

#include "AliMUONChamberTrigger.h"
#include "AliMUONSegmentationTrigger.h"
#include "AliMUONResponseTrigger.h"
#include <TObjArray.h>
#include <TMath.h>
#include <iostream.h>

ClassImp(AliMUONChamberTrigger)

//-------------------------------------------
AliMUONChamberTrigger::AliMUONChamberTrigger() : AliMUONChamber()
{
// Default Constructor
}

//-------------------------------------------
void AliMUONChamberTrigger::DisIntegration(Float_t eloss, Float_t tof, 
					   Float_t xhit, Float_t yhit, Float_t zhit, 
					   Int_t& nnew,
					   Float_t newclust[6][500]) 
{
//    
//  Generates pad hits (simulated cluster) 
//  using the segmentation and the response model

  Int_t twentyNano;
  if (tof<75*TMath::Power(10,-9)) {
    twentyNano=1;
  } else {
    twentyNano=100;
  }

  //  cout << " time = " << tof << " , " << twentyNano << "\n";

  Float_t qp;
  nnew=0;
  for (Int_t i=1; i<=fnsec; i++) {
    AliMUONSegmentation * segmentation=
      (AliMUONSegmentation*) (*fSegmentation)[i-1];
    
// Find the module & strip Id. which has fired
    Int_t ix,iy;
    
    segmentation->GetPadIxy(xhit,yhit,0,ix,iy);
    segmentation->SetPad(ix,iy);

// treatment of GEANT hits w/o corresponding strip (due to the fact that
// the 2 geometries are computed in a very slightly different way) 
    if (ix==0&&iy==0) {
      cout << " AliMUONChamberTrigger hit w/o strip " << xhit << " , " << yhit << "\n";
    } else {          
      // --- store signal information for this strip
      newclust[0][nnew]=1.;                       // total charge
      newclust[1][nnew]=ix;                       // ix-position of pad
      newclust[2][nnew]=iy;                       // iy-position of pad
      newclust[3][nnew]=twentyNano;               // time of flight
      newclust[4][nnew]=segmentation->ISector();  // sector id
      newclust[5][nnew]=(Float_t) i;              // counter
      nnew++;
      // set hits
      segmentation->SetHit(xhit,yhit,zhit);
      // get the list of nearest neighbours
      Int_t nList, xList[2], yList[2];
      segmentation->Neighbours(ix,iy,&nList,xList,yList);
      
      for (Int_t j=0; j<nList; j++){
	
	// neighbour real coordinates (just for checks here)
	Float_t x,y,z;
	segmentation->GetPadCxy(xList[j],yList[j],x,y,z);
	// set pad (fx fy & fix fiy are the current pad coord. & Id.)
	segmentation->SetPad(xList[j],yList[j]);
	// get the chamber (i.e. current strip) response
	qp=fResponse->IntXY(segmentation);	  
	
	if (qp > 0.5) {
	  // --- store signal information for neighbours 
	  newclust[0][nnew]=qp;                       // total charge
	  newclust[1][nnew]=segmentation->Ix();       // ix-position of pad
	  newclust[2][nnew]=segmentation->Iy();       // iy-position of pad
	  newclust[3][nnew]=twentyNano;               // time of flight
	  newclust[4][nnew]=segmentation->ISector();  // sector id
	  newclust[5][nnew]=(Float_t) i;              // counter
	  nnew++;
	} // qp > 0.5  
      } // loop on neighbour
    } // endif hit w/o strip
  } // loop over planes
}








