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

// --- ROOT includes ---
#include <TRandom.h>
#include <TMath.h>
#include "AliRun.h"


// --- MUON includes ---
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONHit.h"
#include "AliLog.h"

ClassImp(AliMUONChamber)	

AliMUONChamber::AliMUONChamber()
  : TObject(), 
    fId(0),
    fdGas(0.),
    fdAlu(0.),
    fZ(0.),
    fnsec(1),
    frMin(0.),
    frMax(0.),
    fCurrentCorrel(1), // to avoid mistakes if ChargeCorrelInit is not called
    fSegmentation2(0),
    fResponse(0),
    fGeometry(0),
    fMUON(0)
{
// Default constructor
}

//_______________________________________________________
AliMUONChamber::AliMUONChamber(Int_t id) 
  : TObject(), 
    fId(id),
    fdGas(0.),
    fdAlu(0.),
    fZ(0.),
    fnsec(1),
    frMin(0.),
    frMax(0.),
    fCurrentCorrel(1), // to avoid mistakes if ChargeCorrelInit is not called
    fSegmentation2(0),
    fResponse(0),
    fGeometry(0),
    fMUON(0)
{

    // muon
    fMUON = (AliMUON*)gAlice->GetModule("MUON");
    if (!fMUON) {
      AliFatal("MUON detector not defined.");
      return;
    }  

    // new segmentation
    fSegmentation2 = new TObjArray(2);
    fSegmentation2->AddAt(0,0);
    fSegmentation2->AddAt(0,1);
 
    fGeometry = new AliMUONGeometryModule(fId);

}

//_______________________________________________________
AliMUONChamber::AliMUONChamber(const AliMUONChamber& rChamber)
  : TObject(rChamber)
{
  // Protected copy constructor

  AliFatal("Not implemented.");
  // Dummy copy constructor
}

//_______________________________________________________
AliMUONChamber::~AliMUONChamber() 
{
  // Destructor

  if (fSegmentation2) {
    fSegmentation2->Delete();
    delete fSegmentation2;
  }
  
  delete fGeometry;
}

//_______________________________________________________
AliMUONChamber & AliMUONChamber::operator =(const AliMUONChamber& rhs)
{
  // Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}

//_______________________________________________________
Bool_t  AliMUONChamber::IsSensId(Int_t volId) const 
{
  // Returns true if the volume specified by volId is in the list
  // of sesitive volumes for this chamber

  return fGeometry->IsSensitiveVolume(volId);
}  


//_____________________________________________________
void AliMUONChamber::ChargeCorrelationInit() {
  // Initialisation of charge correlation for current hit
  // the value is stored, and then used by Disintegration
  if (fnsec==1) 
    fCurrentCorrel =1;
  else 
    // exponential is here to avoid eventual problems in 0 
    // factor 2 because chargecorrel is q1/q2 and not q1/qtrue
    fCurrentCorrel = TMath::Exp(gRandom->Gaus(0,fResponse->ChargeCorrel()/2));
}

//_______________________________________________________
void AliMUONChamber::InitGeo(Float_t /*zpos*/)
{
  //    sensitive gas gap
  fdGas= 0.5;
  //    3% radiation length of aluminum (X0=8.9 cm)      
  fdAlu= 3.0/100*8.9;
}
//_______________________________________________________
//
//                  NEW SEGMENTATION
//_______________________________________________________
void AliMUONChamber::Init(Int_t flag)
{
  // Initalisation ..
  //
  // ... for chamber segmentation

  if (!flag)    AliFatal("wrong segmentation type.");


  if (fSegmentation2->At(0)) 
    ((AliMUONGeometrySegmentation*) fSegmentation2->At(0))->Init(fId);

  if (fnsec==2) {
    if (fSegmentation2->At(1))
      ((AliMUONGeometrySegmentation*) fSegmentation2->At(1))->Init(fId);
  }
}

//_________________________________________________________________
Int_t   AliMUONChamber::SigGenCond(AliMUONHit *hit)
{
  // Ask segmentation if signal should be generated 

  Float_t x = hit->X();
  Float_t y = hit->Y();
  Float_t z = hit->Z();
  Int_t  id = hit->DetElemId();
  
  if (fnsec==1) {
    return  ((AliMUONGeometrySegmentation*)fSegmentation2->At(0))->SigGenCond(id, x, y, z);
  } else {
    return (((AliMUONGeometrySegmentation*) fSegmentation2->At(0))
	    ->SigGenCond(id, x, y, z)) ||
      (((AliMUONGeometrySegmentation*) fSegmentation2->At(1))
       ->SigGenCond(id, x, y, z)) ;
  }
}

//_________________________________________________________________
void    AliMUONChamber::SigGenInit(AliMUONHit *hit)
{
  //
  // Initialisation of segmentation for hit
  //  
  Float_t x = hit->X();
  Float_t y = hit->Y();
  Float_t z = hit->Z();
  Int_t  id = hit->DetElemId();

  if (fnsec==1) {
    ((AliMUONGeometrySegmentation*) fSegmentation2->At(0))->SigGenInit(id, x, y, z) ;
  } else {
    ((AliMUONGeometrySegmentation*) fSegmentation2->At(0))->SigGenInit(id, x, y, z) ;
    ((AliMUONGeometrySegmentation*) fSegmentation2->At(1))->SigGenInit(id, x, y, z) ;
  }
}

//_______________________________________________________
void AliMUONChamber::DisIntegration(AliMUONHit *hit, 
				    Int_t& nnew,Float_t newclust[6][500]) 
{
  //    
  //  Generates pad hits (simulated cluster) 
  //  using the segmentation and the response model 
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
  for (Int_t i=1; i<=fnsec; i++) {
    Float_t qcath = qtot * (i==1? fCurrentCorrel : 1/fCurrentCorrel);

    AliMUONGeometrySegmentation* segmentation=
      (AliMUONGeometrySegmentation*) fSegmentation2->At(i-1);

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
