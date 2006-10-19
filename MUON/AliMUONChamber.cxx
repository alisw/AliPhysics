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

// -----------------------
// Class AliMUONChamber
// -----------------------
// MUON tracking chamber class
// now only providing DisIntegration function

// --- ROOT includes ---
#include <TRandom.h>
#include <TMath.h>
#include "AliRun.h"


// --- MUON includes ---
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONSegmentation.h"
#include "AliMUONHit.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONChamber)	
/// \endcond

//_______________________________________________________
AliMUONChamber::AliMUONChamber()
  : TObject(), 
    fId(0),
    fCurrentCorrel(1), // to avoid mistakes if ChargeCorrelInit is not called
    fResponse(0),
    fMUON(0)
{
/// Default constructor

  AliDebug(1, Form("default (empty) ctor this = %p", this));
}

//_______________________________________________________
AliMUONChamber::AliMUONChamber(Int_t id) 
  : TObject(), 
    fId(id),
    fCurrentCorrel(1), // to avoid mistakes if ChargeCorrelInit is not called
    fResponse(0),
    fMUON(0)
{
/// Standard constructor

    // muon
    fMUON = (AliMUON*)gAlice->GetModule("MUON");
    if (!fMUON) {
      AliFatal("MUON detector not defined.");
      return;
    }  

  AliDebug(1, Form("ctor this = %p", this) ); 
}

//_______________________________________________________
AliMUONChamber::~AliMUONChamber() 
{
/// Destructor

  AliDebug(1, Form("dtor this = %p", this));
}

//_____________________________________________________
void AliMUONChamber::ChargeCorrelationInit() 
{
/// Initialisation of charge correlation for current hit
/// the value is stored, and then used by Disintegration

  // exponential is here to avoid eventual problems in 0 
  // factor 2 because chargecorrel is q1/q2 and not q1/qtrue
    fCurrentCorrel = TMath::Exp(gRandom->Gaus(0,fResponse->ChargeCorrel()/2));
}

//_______________________________________________________
void AliMUONChamber::DisIntegration(AliMUONHit *hit, 
				    Int_t& nnew,Float_t newclust[6][500]) 
{
///  Generates pad hits (simulated cluster) 
///  using the segmentation and the response model 

  Float_t dx, dy;

  Float_t  xhit = hit->X();
  Float_t  yhit = hit->Y();
  Float_t  zhit = hit->Z();
  Int_t      id = hit->DetElemId();
  Float_t eloss = hit->Eloss();

  //
  // Width of the integration area
  //
  dx=fResponse->SigmaIntegration()*fResponse->ChargeSpreadX();
  dy=fResponse->SigmaIntegration()*fResponse->ChargeSpreadY();
  //
  // Get pulse height from energy loss
  Float_t qtot = fResponse->IntPH(eloss);
  //
  // Loop Over Pads
    
  Float_t qp; 
  nnew=0;
    
  // Cathode plane loop
  for (Int_t i = 1; i <= 2; i++) {
    Float_t qcath = qtot * (i==1? fCurrentCorrel : 1/fCurrentCorrel);
    
    AliMUONGeometrySegmentation* segmentation=
      fMUON->GetSegmentation()->GetModuleSegmentationByDEId(id, i-1); 

    for (segmentation->FirstPad(id, xhit, yhit, zhit, dx, dy); 
	 segmentation->MorePads(id); 
	 segmentation->NextPad(id)) 
      {
	qp=fResponse->IntXY(id, segmentation);
	qp=TMath::Abs(qp);
	//
	//
	if (qp > 1.e-4) 
	  {
	    if (nnew >= 500) // Perform a bounds check on nnew since it is assumed
	      // newclust only contains 500 elements.
	      {
		AliError("Limit of 500 pad responses reached.");
		return;
	      };
	    //
	    // --- store signal information
	    newclust[0][nnew]=qcath;                     // total charge
	    newclust[1][nnew]=segmentation->Ix();       // ix-position of pad
	    newclust[2][nnew]=segmentation->Iy();       // iy-position of pad
	    newclust[3][nnew]=qp * qcath;                // charge on pad
	    newclust[4][nnew]=segmentation->ISector();  // sector id
	    newclust[5][nnew]=(Float_t) i;              // counter
	    nnew++;
		
	  }
      } // Pad loop
  } // Cathode plane loop
}
