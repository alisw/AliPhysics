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
Revision 1.6  2000/10/09 14:01:12  morsch
Double inclusion of AliResponse removed.

Revision 1.5  2000/07/03 11:54:57  morsch
AliMUONSegmentation and AliMUONHitMap have been replaced by AliSegmentation and AliHitMap in STEER
The methods GetPadIxy and GetPadXxy of AliMUONSegmentation have changed name to GetPadI and GetPadC.

Revision 1.4  2000/06/29 12:34:09  morsch
AliMUONSegmentation class has been made independent of AliMUONChamber. This makes
it usable with any other geometry class. The link to the object to which it belongs is
established via an index. This assumes that there exists a global geometry manager
from which the pointer to the parent object can be obtained (in our case gAlice).

Revision 1.3  2000/06/28 15:16:35  morsch
(1) Client code adapted to new method signatures in AliMUONSegmentation (see comments there)
to allow development of slat-muon chamber simulation and reconstruction code in the MUON
framework. The changes should have no side effects (mostly dummy arguments).
(2) Hit disintegration uses 3-dim hit coordinates to allow simulation
of chambers with overlapping modules (MakePadHits, Disintegration).

Revision 1.2  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

Revision 1.1.2.5  2000/06/09 21:27:01  morsch
Most coding rule violations corrected.

Revision 1.1.2.4  2000/05/05 11:34:12  morsch
Log inside comments.

Revision 1.1.2.3  2000/05/05 10:09:52  morsch
Log messages included
*/

// --- MUON includes ---
#include "AliMUONChamber.h"

// --- ROOT includes ---

#include "TRandom.h"
#include "TMath.h"

ClassImp(AliMUONChamber)	

    AliMUONChamber::AliMUONChamber()
{
// Default constructor
    fSegmentation = new TObjArray(2);
    (*fSegmentation)[0] = 0;
    (*fSegmentation)[1] = 0;    
    fResponse=0;
    fnsec=1;
    fReconstruction=0;
    fId=0;
    // to avoid mistakes if ChargeCorrelInit is not called
    fChargeCorrel = 0;
    fCurrentCorrel =1;
}

    AliMUONChamber::AliMUONChamber(Int_t id) 
{
// Construtor with chamber id 
    fSegmentation = new TObjArray(2);
    (*fSegmentation)[0] = 0;
    (*fSegmentation)[1] = 0;    
    fResponse=0;
    fnsec=1;
    fReconstruction=0;
    fId=id;
    // to avoid mistakes if ChargeCorrelInit is not called
    fChargeCorrel = 0;
    fCurrentCorrel =1;
}

AliMUONChamber::~AliMUONChamber() 
{
// Destructor
    if (fSegmentation) delete fSegmentation;
}

AliMUONChamber::AliMUONChamber(const AliMUONChamber& rChamber)
{
// Dummy copy constructor
    ;
}


void AliMUONChamber::Init()
{
// Initalisation ..
//
// ... for chamber segmentation
    if ((*fSegmentation)[0]) 
    ((AliSegmentation *) (*fSegmentation)[0])->Init(fId);

    if (fnsec==2) {
	if ((*fSegmentation)[1])
	((AliSegmentation *) (*fSegmentation)[1])->Init(fId);
    }
}

Int_t   AliMUONChamber::SigGenCond(Float_t x, Float_t y, Float_t z)
{
// Ask segmentation if signal should be generated 
    if (fnsec==1) {
	return ((AliSegmentation*) (*fSegmentation)[0])
	    ->SigGenCond(x, y, z) ;
    } else {
	return (((AliSegmentation*) (*fSegmentation)[0])
		->SigGenCond(x, y, z)) ||
	    (((AliSegmentation*) (*fSegmentation)[1])
	     ->SigGenCond(x, y, z)) ;
    }
}


void    AliMUONChamber::SigGenInit(Float_t x, Float_t y, Float_t z)
{
//
// Initialisation of segmentation for hit
//  
    if (fnsec==1) {
	((AliSegmentation*) (*fSegmentation)[0])->SigGenInit(x, y, z) ;
    } else {
	((AliSegmentation*) (*fSegmentation)[0])->SigGenInit(x, y, z) ;
	((AliSegmentation*) (*fSegmentation)[1])->SigGenInit(x, y, z) ;
    }
}

void AliMUONChamber::ChargeCorrelationInit() {
// Initialisation of charge correlation for current hit
// the value is stored, and then used by Disintegration
if (fnsec==1) 
    fCurrentCorrel =1;
else 
    // exponential is here to avoid eventual problems in 0 
    fCurrentCorrel = TMath::Exp(gRandom->Gaus(0,fChargeCorrel/2));
}

void AliMUONChamber::DisIntegration(Float_t eloss, Float_t tof, 
				    Float_t xhit, Float_t yhit, Float_t zhit,
				    Int_t& nnew,Float_t newclust[6][500]) 
{
//    
//  Generates pad hits (simulated cluster) 
//  using the segmentation and the response model 
    Float_t dx, dy;
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
    
    Float_t qcheck=0, qp;
    nnew=0;
    
    // Cathode plane loop
    for (Int_t i=1; i<=fnsec; i++) {
	qcheck=0;
	Float_t qcath = qtot * (i==1? fCurrentCorrel : 1/fCurrentCorrel);
	AliSegmentation * segmentation=
	    (AliSegmentation *) (*fSegmentation)[i-1];
	for (segmentation->FirstPad(xhit, yhit, zhit, dx, dy); 
	     segmentation->MorePads(); 
	     segmentation->NextPad()) 
	{
	    qp=fResponse->IntXY(segmentation);
	    qp=TMath::Abs(qp);
//
//
	    if (qp > 1.e-4) {
		qcheck+=qp*qcath;
	    //
	    // --- store signal information
		newclust[0][nnew]=qcath;                     // total charge
		newclust[1][nnew]=segmentation->Ix();       // ix-position of pad
		newclust[2][nnew]=segmentation->Iy();       // iy-position of pad
		newclust[3][nnew]=qp * qcath;                // charge on pad
		newclust[4][nnew]=segmentation->ISector();  // sector id
		newclust[5][nnew]=(Float_t) i;              // counter
		nnew++;
//		if (i==2) printf("\n i, nnew, q %d %d %f", i, nnew, qp*qcath);
		
	    }
	} // Pad loop
	
    } // Cathode plane loop
}



void AliMUONChamber::InitGeo(Float_t zpos)
{
//    sensitive gas gap
      fdGas= 0.5;
//    3% radiation length of aluminum (X0=8.9 cm)      
      fdAlu= 3.0/100*8.9;
}


AliMUONChamber & AliMUONChamber::operator =(const AliMUONChamber& rhs)
{
// Dummy assignment operator
    return *this;
}
