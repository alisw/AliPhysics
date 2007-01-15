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

#include <TMath.h>

#include "AliMUONChamberTrigger.h"
#include "AliMUONResponseTrigger.h"
#include "AliMUONHit.h"
#include "AliMUON.h"
#include "AliMUONSegmentation.h"
#include "AliMUONGeometrySegmentation.h"
#include "AliMUONGeometryTransformer.h"
#include "AliLog.h"

///
/// \class AliMUONChamberTrigger
///
/// Implementation of AliMUONChamber for the trigger
///
/// This class is to be deprecated.
///

/// \cond CLASSIMP
ClassImp(AliMUONChamberTrigger)
/// \endcond

//-------------------------------------------

AliMUONChamberTrigger::AliMUONChamberTrigger()
  : AliMUONChamber(),
    fkGeomTransformer(0)
{
/// Default constructor
}

//-------------------------------------------

AliMUONChamberTrigger:: ~AliMUONChamberTrigger()
{
/// Destructor
}

//-------------------------------------------

AliMUONChamberTrigger::AliMUONChamberTrigger(Int_t id,
                              const AliMUONGeometryTransformer* kGeometryTransformer) 
  : AliMUONChamber(id),
    fkGeomTransformer(kGeometryTransformer)
{
/// Constructor using chamber id
}

//-------------------------------------------
void AliMUONChamberTrigger::DisIntegration(AliMUONHit* hit,
					   Int_t& nnew,
					   Float_t newclust[6][500]) 
{
///  Generates pad hits (simulated cluster) 
///  using the segmentation and the response model


  Float_t   tof = hit->Age();
  Float_t  xhit = hit->X();
  Float_t  yhit = hit->Y();
  Float_t  zhit = hit->Z();
  Int_t      id = hit->DetElemId();

  Int_t twentyNano;
  if (tof<75*TMath::Power(10,-9)) {
    twentyNano=1;
  } else {
    twentyNano=100;
  }

  Float_t qp;
  nnew=0;
  for (Int_t i = 1; i <= 2; i++) {

    AliMUONGeometrySegmentation* segmentation=
      fMUON->GetSegmentation()->GetModuleSegmentationByDEId(id, i-1); 

    
// Find the module & strip Id. which has fired
    Int_t ix(-1);
    Int_t iy(-1);
    segmentation->GetPadI(id,xhit,yhit,0,ix,iy);
// treatment of GEANT hits w/o corresponding strip (due to the fact that
// geometry & segmentation are computed in a very slightly different way) 
    if ( ix<0 || iy<0 ) 
    {
      Float_t lx,ly,lz;
      fkGeomTransformer->Global2Local(id,xhit,yhit,0,lx,ly,lz);
      AliWarning(Form("AliMUONChamberTrigger hit w/o strip %i-%d %e %e "
                      "local %e %e %e ix,iy=%d,%d\n",id,i-1,xhit,yhit,lx,ly,lz,ix,iy));
    } else 
    {          
      segmentation->SetPad(id,ix,iy);
      if (xhit<0) ix = -ix;
      //    printf(" fId id fnsec xhit yhit zhit ix iy %i %i %i %f %f %f %i %i \n",fId,i,id,xhit,yhit,zhit,ix,iy);
      //     if (ix < 0 || ix > 10000) return;
      //     if (iy < 0 || iy > 10000) return;
      
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
	    segmentation->SetHit(id,xhit,yhit,zhit);
	    // get the list of nearest neighbours
	    Int_t nList, xList[10], yList[10];
	    segmentation->Neighbours(id,ix,iy,&nList,xList,yList);
	    
	    qp = 0;
	    for (Int_t j=0; j<nList; j++){       // loop over neighbours	  
		if (xList[j]!=0) {                 // existing neighbour	    
		    if (j==0||j==5||qp!=0) {         // built-up cluster-size
			
			// neighbour real coordinates (just for checks here)
			Float_t x,y,z;
			segmentation->GetPadC(id,xList[j],yList[j],x,y,z);
			// set pad (fx fy & fix fiy are the current pad coord. & Id.)
			segmentation->SetPad(id,xList[j],yList[j]);	  
			// get the chamber (i.e. current strip) response
			qp=fResponse->IntXY(id,segmentation);	  
			
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






