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

/////////////////////////////////////////////////////////////////////////
//                                                                     //
//  Tracking sector                                                    //
//                                                                     //
/////////////////////////////////////////////////////////////////////////                       
                                
#include <TObject.h>

#include "AliRun.h"

#include "AliTRD.h" 
#include "AliTRDgeometry.h" 
#include "AliTRDcluster.h" 
#include "AliTRDtimeBin.h" 
#include "AliTRDtrackingSector.h" 
#include "AliTRDparameter.h"

ClassImp(AliTRDtrackingSector) 

//_______________________________________________________
AliTRDtrackingSector::~AliTRDtrackingSector()
{
  //
  // Destructor
  //

  delete fTimeBin;

}

//_______________________________________________________
AliTRDtimeBin &AliTRDtrackingSector::operator[](Int_t i)
{
  //
  // Index operator 
  //

  return *(fTimeBin+i);

}

//_______________________________________________________
void AliTRDtrackingSector::SetUp()
{ 
  //
  // Initialization
  //

  AliTRD *trd = (AliTRD*) gAlice->GetDetector("TRD");
  fGeom = trd->GetGeometry();

  fTimeBinSize = fGeom->GetTimeBinSize();

  fN = AliTRDgeometry::Nplan() * (Int_t(AliTRDgeometry::DrThick()
                                       /fTimeBinSize) + 1);

  fTimeBin = new AliTRDtimeBin[fN]; 

  if (!fPar) {
    fPar = new AliTRDparameter("TRDparameter","Standard TRD parameter");
  }

  fTimeBinSize = fPar->GetTimeBinSize();

}

//______________________________________________________
Double_t AliTRDtrackingSector::GetX(Int_t tb) const
{
  //
  // Get global x coordinate
  //

  if( (tb<0) || (tb>fN-1)) {
    fprintf(stderr,"AliTRDtrackingSector::GetX: TimeBin index is out of range !\n");
    return -99999.;
  }
  else { 
    
    Int_t tbPerPlane = fN/AliTRDgeometry::Nplan();
    Int_t localTb = tbPerPlane - tb%tbPerPlane - 1;

    Int_t plane = tb/tbPerPlane;
    Float_t t0 = fPar->GetTime0(plane);
    Double_t x = t0 - (localTb + 0.5) * fTimeBinSize;

    return x;

  }

}

//______________________________________________________
Int_t AliTRDtrackingSector::GetTimeBinNumber(Double_t x) const
{
  //
  // Returns the time bin number
  //

  Float_t rOut = fPar->GetTime0(AliTRDgeometry::Nplan()-1); 
  Float_t rIn  = fPar->GetTime0(0) - AliTRDgeometry::DrThick();


  if(x >= rOut) return fN-1;
  if(x <= rIn)  return 0;

  Int_t plane;
  for (plane = AliTRDgeometry::Nplan()-1; plane >= 0; plane--) {
    if(x > (fPar->GetTime0(plane) - AliTRDgeometry::DrThick())) break;
  }  
 
  Int_t tbPerPlane = fN/AliTRDgeometry::Nplan();
  Int_t localTb = Int_t((fPar->GetTime0(plane)-x)/fTimeBinSize);

  if((localTb < 0) || (localTb >= tbPerPlane)) {
    printf("AliTRDtrackingSector::GetTimeBinNumber: \n");
    printf("local time bin %d is out of bounds [0..%d]: x = %f \n",
	   localTb,tbPerPlane-1,x);
    return -1;
  }
      
  Int_t timeBin = (plane + 1) * tbPerPlane - 1 - localTb;

  return timeBin;
}

//______________________________________________________
Int_t AliTRDtrackingSector::GetTimeBin(Int_t det, Int_t localTb) const 
{
  //
  // Time bin
  //

  Int_t plane = fGeom->GetPlane(det);

  Int_t tbPerPlane = fN/AliTRDgeometry::Nplan();

  Int_t timeBin = (plane + 1) * tbPerPlane - 1 - localTb;

  return timeBin;

}


//______________________________________________________

Bool_t AliTRDtrackingSector::TECframe(Int_t tb, Double_t y, Double_t z) const
{
  //  
  // Returns <true> if point defined by <x(tb),y,z> is within 
  // the TEC G10 frame, otherwise returns <false>  
  //  

  if((tb > (fN-1)) || (tb < 0)) return kFALSE; 

  Int_t tbPerPlane = fN/AliTRDgeometry::Nplan();
  Int_t plane = tb/tbPerPlane;
  
  Double_t x = GetX(tb);
  y = TMath::Abs(y);

  if((y > fGeom->GetChamberWidth(plane)/2.) &&
     (y < x*TMath::Tan(0.5*AliTRDgeometry::GetAlpha()))) return kTRUE; 

  Double_t zmin, zmax;
  Float_t  fRowPadSize, fRow0;
  Int_t    nPadRows;

  for(Int_t iCha = 1; iCha < AliTRDgeometry::Ncham(); iCha++) {

    fRow0 = fPar->GetRow0(plane,iCha-1,0);
    fRowPadSize = fPar->GetRowPadSize(plane,iCha-1,0);
    nPadRows = fPar->GetRowMax(plane,iCha-1,0);
    zmin = fRow0 - fRowPadSize/2 + fRowPadSize * nPadRows;

    fRow0 = fPar->GetRow0(plane,iCha,0);
    fRowPadSize = fPar->GetRowPadSize(plane,iCha,0);
    zmax = fRow0 - fRowPadSize/2;

    if((z > zmin) && (z < zmax)) return kTRUE;     
  }

  return kFALSE;

}
