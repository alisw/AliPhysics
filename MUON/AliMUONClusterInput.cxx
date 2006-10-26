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

// ----------------------------
// Class AliMUONClusterInput
// ----------------------------
// Global data service for hit reconstruction
// Author: to be added

#include "AliMUONClusterInput.h"

#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometrySegmentation.h"
#include "AliMUONSegFactory.h"
#include "AliMUONSegmentation.h"
#include "AliMUONConstants.h"
#include "AliMUONMathieson.h"
#include "AliMUONRawCluster.h"
#include "AliMUONDigit.h"

#include "AliLog.h"

#include <TClonesArray.h>
#include <TMinuit.h>
#include <TGeoManager.h>

/// \cond CLASSIMP
ClassImp(AliMUONClusterInput)
/// \endcond

AliMUONClusterInput*        AliMUONClusterInput::fgClusterInput = 0; 
TMinuit*                    AliMUONClusterInput::fgMinuit = 0; 
AliMUONMathieson*           AliMUONClusterInput::fgMathieson = 0; 
AliMUONGeometryTransformer* AliMUONClusterInput::fgTransformer = 0; 
AliMUONSegmentation*        AliMUONClusterInput::fgSegmentation = 0; 

//______________________________________________________________________________
AliMUONClusterInput::AliMUONClusterInput()
  : TObject(),
    fNseg(0),
    fChamber(0),
    fCluster(0),
    fZ(0.),
    fChargeCorrel(1.),
    fDetElemId(0)
  
{
/// Default constructor

  fDigits[0]=0;
  fDigits[1]=0;
  fSegmentation2[0]=0;
  fSegmentation2[1]=0;
}

//______________________________________________________________________________
AliMUONClusterInput* AliMUONClusterInput::Instance()
{
/// return pointer to the singleton instance
    if (fgClusterInput == 0) {
	fgClusterInput = new AliMUONClusterInput();
	fgMinuit = new TMinuit(8);
	
	// Create segmentation with activated Root geometry  
	if ( ! gGeoManager ) {
	  AliFatalClass("Geometry not loaded.");
	  return fgClusterInput;
	}  
        fgTransformer = new AliMUONGeometryTransformer(true);
        fgTransformer->ReadGeometryData("volpath.dat", gGeoManager);
        AliMUONSegFactory factory(fgTransformer);
	fgSegmentation = factory.CreateSegmentation(); 
    }
    
    return fgClusterInput;
}

//______________________________________________________________________________
AliMUONClusterInput::~AliMUONClusterInput()
{
/// Destructor
    delete fgMinuit;
    delete fgMathieson;
    delete fgTransformer;
    delete fgSegmentation;
    fgMinuit = 0;
    fgMathieson = 0;
}

//______________________________________________________________________________
void AliMUONClusterInput::SetDigits(Int_t chamber, Int_t idDE, TClonesArray* dig1, TClonesArray* dig2)
{
  /// Set pointer to digits with corresponding segmentations and responses (two cathode planes)
    fChamber = chamber;
    fDetElemId = idDE;
    fDigits[0]  = dig1;
    fDigits[1]  = dig2; 
    fNDigits[0] = dig1->GetEntriesFast();
    fNDigits[1] = dig2->GetEntriesFast();
    
    delete fgMathieson;
    fgMathieson = new AliMUONMathieson();

    fSegmentation2[0]= fgSegmentation->GetModuleSegmentationByDEId(fDetElemId, 0);
    fSegmentation2[1]= fgSegmentation->GetModuleSegmentationByDEId(fDetElemId, 1);

    fNseg = 2;
    if (chamber < AliMUONConstants::NTrackingCh()) {
      if (chamber > 1 ) {
	fgMathieson->SetPitch(AliMUONConstants::Pitch());
	fgMathieson->SetSqrtKx3AndDeriveKx2Kx4(AliMUONConstants::SqrtKx3());
	fgMathieson->SetSqrtKy3AndDeriveKy2Ky4(AliMUONConstants::SqrtKy3());
	fChargeCorrel = AliMUONConstants::ChargeCorrel();
      } else {
	fgMathieson->SetPitch(AliMUONConstants::PitchSt1());
	fgMathieson->SetSqrtKx3AndDeriveKx2Kx4(AliMUONConstants::SqrtKx3St1());
	fgMathieson->SetSqrtKy3AndDeriveKy2Ky4(AliMUONConstants::SqrtKy3St1());
	fChargeCorrel = AliMUONConstants::ChargeCorrelSt1();
      }
    }
}

//______________________________________________________________________________
void AliMUONClusterInput::SetDigits(Int_t chamber, Int_t idDE, TClonesArray* dig)
{
/// Set pointer to digits with corresponding segmentations and responses (one cathode plane)

    fChamber = chamber;
    fDetElemId = idDE;
    fDigits[0] = dig;

    fSegmentation2[0]= fgSegmentation->GetModuleSegmentationByDEId(fDetElemId, 0);
    fNseg=1;
}

//______________________________________________________________________________
void  AliMUONClusterInput::SetCluster(AliMUONRawCluster* cluster)
{
/// Set the current cluster
  //PH printf("\n %p \n", cluster);
  fCluster=cluster;
  Float_t qtot;
  Int_t   i, cath, ix, iy;
  AliMUONDigit* digit;
  fNmul[0]=cluster->GetMultiplicity(0);
  fNmul[1]=cluster->GetMultiplicity(1);
  //PH printf("\n %p %p ", fDigits[0], fDigits[1]);
  
  for (cath=0; cath<2; cath++) {
    qtot=0;
    for (i=0; i<fNmul[cath]; i++) {
      // pointer to digit
      digit =(AliMUONDigit*)
		(fDigits[cath]->UncheckedAt(cluster->GetIndex(i,cath)));
	    // pad coordinates
	    ix = digit->PadX();
	    iy = digit->PadY();
	    // pad charge
	    fCharge[i][cath] = digit->Signal();
	    // pad centre coordinates
//	    fSegmentation[cath]->GetPadCxy(ix, iy, x, y);
            // globals kUsed in fitting functions
	    fix[i][cath]=ix;
	    fiy[i][cath]=iy;
	    // total charge per cluster
	    qtot+=fCharge[i][cath];
	    // Current z
	    Float_t xc, yc;
	    fSegmentation2[cath]->GetPadC(fDetElemId,ix,iy,xc,yc,fZ);
	} // loop over cluster digits
	fQtot[cath]=qtot;
	fChargeTot[cath]=Int_t(qtot);  
    }  // loop over cathodes
}

//______________________________________________________________________________
Float_t AliMUONClusterInput::DiscrChargeS1(Int_t i,Double_t *par) 
{
/// Compute the charge on first cathod only.
return DiscrChargeCombiS1(i,par,0);
}

//______________________________________________________________________________
Float_t AliMUONClusterInput::DiscrChargeCombiS1(Int_t i,Double_t *par, Int_t cath) 
{
/// \todo add comment
/// - par[0]    x-position of cluster
/// - param par[1]    y-position of cluster

   Float_t q1;
   fSegmentation2[cath]-> SetPad(fDetElemId, fix[i][cath], fiy[i][cath]);
   //  First Cluster
   fSegmentation2[cath]-> SetHit(fDetElemId, par[0],par[1],fZ);
   q1 = fgMathieson->IntXY(fDetElemId, fSegmentation2[cath]);
       
   Float_t value = fQtot[cath]*q1;
   return value;
}


//______________________________________________________________________________
Float_t AliMUONClusterInput::DiscrChargeS2(Int_t i,Double_t *par) 
{
/// \todo add comment
/// - par[0]    x-position of first  cluster
/// - par[1]    y-position of first  cluster
/// - par[2]    x-position of second cluster
/// - par[3]    y-position of second cluster
/// - par[4]    charge fraction of first  cluster
/// - 1-par[4]  charge fraction of second cluster

  Float_t q1, q2;
  
  fSegmentation2[0]->SetPad(fDetElemId, fix[i][0], fiy[i][0]);
  //  First Cluster
  fSegmentation2[0]->SetHit(fDetElemId, par[0],par[1],fZ);
  q1 = fgMathieson->IntXY(fDetElemId, fSegmentation2[0]);

  //  Second Cluster
  fSegmentation2[0]->SetHit(fDetElemId,par[2],par[3],fZ);
  q2 = fgMathieson->IntXY(fDetElemId, fSegmentation2[0]);
  
  Float_t value = fQtot[0]*(par[4]*q1+(1.-par[4])*q2);
  return value;
}

//______________________________________________________________________________
Float_t AliMUONClusterInput::DiscrChargeCombiS2(Int_t i,Double_t *par, Int_t cath) 
{
/// \todo add comment
/// - par[0]    x-position of first  cluster
/// - par[1]    y-position of first  cluster
/// - par[2]    x-position of second cluster
/// - par[3]    y-position of second cluster
/// - par[4]    charge fraction of first  cluster - first cathode
/// - 1-par[4]  charge fraction of second cluster 
/// - par[5]    charge fraction of first  cluster - second cathode

  Float_t q1, q2;

  fSegmentation2[cath]->SetPad(fDetElemId,fix[i][cath], fiy[i][cath]);
  //  First Cluster
  fSegmentation2[cath]->SetHit(fDetElemId,par[0],par[1],fZ);
  q1 = fgMathieson->IntXY(fDetElemId, fSegmentation2[cath]);
  
  //  Second Cluster
  fSegmentation2[cath]->SetHit(fDetElemId,par[2],par[3],fZ);
  q2 = fgMathieson->IntXY(fDetElemId, fSegmentation2[cath]);
  
  Float_t value;
  if (cath==0) {
    value = fQtot[0]*(par[4]*q1+(1.-par[4])*q2);
  } else {
    value = fQtot[1]*(par[5]*q1+(1.-par[5])*q2);
  }
  return value;
}
