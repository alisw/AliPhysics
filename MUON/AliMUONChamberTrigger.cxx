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
Revision 1.7  2000/10/02 16:58:29  egangler
Cleaning of the code :
-> coding conventions
-> void Streamers
-> some useless includes removed or replaced by "class" statement

Revision 1.6  2000/07/03 11:54:57  morsch
AliMUONSegmentation and AliMUONHitMap have been replaced by AliSegmentation and AliHitMap in STEER
The methods GetPadIxy and GetPadXxy of AliMUONSegmentation have changed name to GetPadI and GetPadC.

Revision 1.5  2000/06/29 12:34:09  morsch
AliMUONSegmentation class has been made independent of AliMUONChamber. This makes
it usable with any other geometry class. The link to the object to which it belongs is
established via an index. This assumes that there exists a global geometry manager
from which the pointer to the parent object can be obtained (in our case gAlice).

Revision 1.4  2000/06/29 06:52:02  pcrochet
pow changed to TMath::Power

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
#include "AliMUONResponseTriggerV1.h"
#include <TObjArray.h>
#include <TMath.h>
#include <iostream.h>

ClassImp(AliMUONChamberTrigger)

//-------------------------------------------

AliMUONChamberTrigger::AliMUONChamberTrigger()
{
// Default constructor
}


AliMUONChamberTrigger::AliMUONChamberTrigger(Int_t id) : AliMUONChamber(id)
{
// Constructor using chamber id
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
    AliSegmentation * segmentation=
      (AliSegmentation*) (*fSegmentation)[i-1];
    
// Find the module & strip Id. which has fired
    Int_t ix,iy;
    
    segmentation->GetPadI(xhit,yhit,0,ix,iy);
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

// cluster-size if AliMUONResponseTriggerV1, nothing if AliMUONResponseTrigger
      if (((AliMUONResponseTrigger*) fResponse)->SetGenerCluster()) {
  
	// set hits
	segmentation->SetHit(xhit,yhit,zhit);
	// get the list of nearest neighbours
	Int_t nList, xList[10], yList[10];
	segmentation->Neighbours(ix,iy,&nList,xList,yList);
	
	qp = 0;
	for (Int_t j=0; j<nList; j++){       // loop over neighbours	  
	  if (xList[j]!=0) {                 // existing neighbour	    
	    if (j==0||j==5||qp!=0) {         // built-up cluster-size
	      
	      // neighbour real coordinates (just for checks here)
	      Float_t x,y,z;
	      segmentation->GetPadC(xList[j],yList[j],x,y,z);
	      // set pad (fx fy & fix fiy are the current pad coord. & Id.)
	      segmentation->SetPad(xList[j],yList[j]);	  
	      // get the chamber (i.e. current strip) response
	      qp=fResponse->IntXY(segmentation);	  
	      
	      if (qp > 0.5) {		
		// --- store signal information for neighbours 
		newclust[0][nnew]=qp;                      // total charge
		newclust[1][nnew]=segmentation->Ix();      // ix-pos. of pad
		newclust[2][nnew]=segmentation->Iy();      // iy-pos. of pad
		newclust[3][nnew]=twentyNano;              // time of flight
		newclust[4][nnew]=segmentation->ISector(); // sector id
		newclust[5][nnew]=(Float_t) i;             // counter
		nnew++;
	      } // qp > 0.5 
	    } // built-up cluster-size
	  } // existing neighbour
	} // loop over neighbours
      } // endif hit w/o strip
    } // loop over planes
  } // if AliMUONResponseTriggerV1
}










