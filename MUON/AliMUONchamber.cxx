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
*/

#include "AliMUONchamber.h"
#include "TMath.h"
ClassImp(AliMUONchamber)	
    AliMUONchamber::AliMUONchamber() 
{
    fSegmentation = new TObjArray(2);
    fResponse=0;
    fnsec=1;
}

void AliMUONchamber::Init()
{
    
    ((AliMUONsegmentation *) (*fSegmentation)[0])->Init(this);
    if (fnsec==2) {
	((AliMUONsegmentation *) (*fSegmentation)[1])->Init(this);
    }
    
}

void AliMUONchamber::DisIntegration(Float_t eloss, Float_t xhit, Float_t yhit,
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
    for (Int_t i=1; i<=fnsec; i++) {
	qcheck=0;
	AliMUONsegmentation * segmentation=(AliMUONsegmentation *) (*fSegmentation)[i-1];
	for (segmentation->FirstPad(xhit, yhit, dx, dy); 
	     segmentation->MorePads(); 
	     segmentation->NextPad()) 
	{
	    qp=fResponse->IntXY(segmentation);
	    qp=TMath::Abs(qp);
	    
//
//
	    if (qp > 1.e-4) {
		qcheck+=qp;
	    //
	    // --- store signal information
		newclust[0][nnew]=qtot;
		newclust[1][nnew]=segmentation->Ix();
		newclust[2][nnew]=segmentation->Iy();
		newclust[3][nnew]=qp * qtot;
		newclust[4][nnew]=segmentation->ISector();
		newclust[5][nnew]=(Float_t) i;
//		printf("\n pad hit %d %d %f %f \n",nnew,i,newclust[1][nnew],newclust[2][nnew]);
		nnew++;

		
	    }
	} // Pad loop
//	printf("\n check sum is %f %f %f %f %d \n",qcheck,qtot,xhit,yhit,fGid);
    } // Cathode plane loop
}



    void AliMUONchamber::InitGeo(Float_t)
{
//    sensitive gas gap
      fdGas= 0.5;
//    3% radiation length of aluminum (X0=8.9 cm)      
      fdAlu= 3.0/100*8.9;
}
