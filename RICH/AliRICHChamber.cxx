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
  Revision 1.9  2001/02/13 20:15:34  jbarbosa
  Removed fNsec (number of cathodes - obsolete) related loops and calls.

  Revision 1.8  2000/12/18 17:45:43  jbarbosa
  Cleaned up PadHits object.

  Revision 1.7  2000/10/03 21:44:09  morsch
  Use AliSegmentation and AliHit abstract base classes.

  Revision 1.6  2000/10/02 15:44:37  jbarbosa
  Fixed forward declarations.

  Revision 1.5  2000/07/13 16:19:45  fca
  Mainly coding conventions + some small bug fixes

  Revision 1.4  2000/06/30 16:48:58  dibari
  New function GenerateTresholds() for pedestal simulation.

  Revision 1.3  2000/06/12 15:17:58  jbarbosa
  Cleaned up version.

  Revision 1.2  2000/05/18 13:45:57  jbarbosa
  Fixed feedback photon origin coordinates

  Revision 1.1  2000/04/19 12:57:20  morsch
  Newly structured and updated version (JB, AM)

*/


#include "AliRICHChamber.h"

#include <TLorentzVector.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TObjArray.h>
#include <TRotMatrix.h>
#include <AliRICHTresholdMap.h>
#include <AliSegmentation.h>
#include <AliRICHSegmentationV0.h>
#include <AliRICHGeometry.h>
#include <AliRICHResponse.h>

ClassImp(AliRICHChamber)	
    
AliRICHChamber::AliRICHChamber() 
{

//
// Chamber object constructor

    fSegmentation = 0;
    fResponse = 0;
    fGeometry = 0;
    fTresh = 0;
    frMin = 0.1;
    frMax = 140;
    for(Int_t i=0; i<50; ++i) fIndexMap[i] = 0;
}

AliRICHChamber::AliRICHChamber(const AliRICHChamber& Chamber)
{
// Copy Constructor
}


AliRICHResponse* AliRICHChamber::GetResponseModel()
{
//  
//  Get reference to response model
    return fResponse;
}

void   AliRICHChamber::ResponseModel(AliRICHResponse* thisResponse)
{
// Configure response model
    fResponse=thisResponse;
}

void AliRICHChamber::Init(Int_t id)
{
// Initialise chambers
    fSegmentation->Init(id);
}

void AliRICHChamber::LocaltoGlobal(Float_t pos[3],Float_t Globalpos[3])
{

// Local coordinates to global coordinates transformation

    Double_t *fMatrix;
    fMatrix =  fChamberMatrix->GetMatrix();
    Globalpos[0]=pos[0]*fMatrix[0]+pos[1]*fMatrix[3]+pos[2]*fMatrix[6];
    Globalpos[1]=pos[0]*fMatrix[1]+pos[1]*fMatrix[4]+pos[2]*fMatrix[7];
    Globalpos[2]=pos[0]*fMatrix[2]+pos[1]*fMatrix[5]+pos[2]*fMatrix[8];
    Globalpos[0]+=fChamberTrans[0];
    Globalpos[1]+=fChamberTrans[1];
    Globalpos[2]+=fChamberTrans[2];
}

void AliRICHChamber::GlobaltoLocal(Float_t pos[3],Float_t Localpos[3])
{

// Global coordinates to local coordinates transformation

    Double_t *fMatrixOrig;
    TMatrix fMatrixCopy(3,3);
    fMatrixOrig = fChamberMatrix->GetMatrix();
    for(Int_t i=0;i<3;i++)
      {
	for(Int_t j=0;j<3;j++)
	  fMatrixCopy(j,i)=fMatrixOrig[j+3*i];
      }
    fMatrixCopy.Invert();
    //Int_t elements=fMatrixCopy.GetNoElements();
    //printf("Elements:%d\n",elements);
    //fMatrixOrig= (Double_t*) fMatrixCopy;
    Localpos[0] = pos[0] - fChamberTrans[0];
    Localpos[1] = pos[1] - fChamberTrans[1];
    Localpos[2] = pos[2] - fChamberTrans[2];
    //printf("r1:%f, r2:%f, r3:%f\n",Localpos[0],Localpos[1],Localpos[2]);
    //printf("t1:%f t2:%f t3:%f\n",fChamberTrans[0],fChamberTrans[1],fChamberTrans[2]);
    Localpos[0]=Localpos[0]*fMatrixCopy(0,0)+Localpos[1]*fMatrixCopy(0,1)+Localpos[2]*fMatrixCopy(0,2);
    Localpos[1]=Localpos[0]*fMatrixCopy(1,0)+Localpos[1]*fMatrixCopy(1,1)+Localpos[2]*fMatrixCopy(1,2);
    Localpos[2]=Localpos[0]*fMatrixCopy(2,0)+Localpos[1]*fMatrixCopy(2,1)+Localpos[2]*fMatrixCopy(2,2);
    //Localpos[0]-=fChamberTrans[0];
    //Localpos[1]-=fChamberTrans[1];
    //Localpos[2]-=fChamberTrans[2];
} 


void AliRICHChamber::DisIntegration(Float_t eloss, Float_t xhit, Float_t yhit,
				    Int_t& nnew,Float_t newclust[5][500],ResponseType res) 
{
//    
//  Generates pad hits (simulated cluster) 
//  using the segmentation and the response model
    
    Float_t dx, dy;
    Float_t local[3];
    //Float_t source[3];
    Float_t global[3];
    //
    // Width of the integration area
    //
    dx=(fResponse->SigmaIntegration())*(fResponse->ChargeSpreadX());
    dy=(fResponse->SigmaIntegration())*(fResponse->ChargeSpreadY());
    //
    // Get pulse height from energy loss and generate feedback photons
    Float_t qtot=0;

    local[0]=xhit;
    // z-position of the wires relative to the RICH mother volume 
    // (2 mmm before CsI) old value: 6.076
    local[1]=1.276 + fGeometry->GetGapThickness()/2  - .2;
    //printf("AliRICHChamber feedback origin:%f",local[1]);
    local[2]=yhit;

    LocaltoGlobal(local,global);

    Int_t nFp=0;
    

    // To calculate wire sag, the origin of y-position must be the middle of the photcathode
    AliRICHSegmentationV0* segmentation = (AliRICHSegmentationV0*) GetSegmentationModel();
    Float_t newy;
    if (yhit>0)
      newy = yhit - segmentation->GetPadPlaneLength()/2;
    else
      newy = yhit + segmentation->GetPadPlaneLength()/2;
    
    if (res==kMip) {
	qtot = fResponse->IntPH(eloss, newy);
	nFp  = fResponse->FeedBackPhotons(global,qtot);
    } else if (res==kCerenkov) {
	qtot = fResponse->IntPH(newy);
	nFp  = fResponse->FeedBackPhotons(global,qtot);
    }

    //printf("Feedbacks:%d\n",nFp);
    
    //
    // Loop Over Pads
    
    Float_t qcheck=0, qp=0;
    
    nnew=0;
    for (fSegmentation->FirstPad(xhit, yhit, 0, dx, dy); 
	 fSegmentation->MorePads(); 
	 fSegmentation->NextPad()) 
      {
	qp= fResponse->IntXY(fSegmentation);
	qp= TMath::Abs(qp);
	
	//printf("Qp:%f Qtot %f\n",qp,qtot);
	
	if (qp > 1.e-4) {
	  qcheck+=qp;
	  //
	  // --- store signal information
	  newclust[0][nnew]=qp*qtot;
	  newclust[1][nnew]=fSegmentation->Ix();
	  newclust[2][nnew]=fSegmentation->Iy();
	  newclust[3][nnew]=fSegmentation->ISector();
	  nnew++;	
	  //printf("Newcluster:%d\n",i);
	}
      } // Pad loop
    //if (fSegmentation->ISector()==2)
      //printf("Nnew:%d\n\n\n\n",nnew);
}


AliRICHChamber& AliRICHChamber::operator=(const AliRICHChamber& rhs)
{
// Assignment operator
    return *this;
    
}


void AliRICHChamber::GenerateTresholds()
{

// Generates random treshold charges for all pads 

  //printf("Pads : %dx%d\n",fSegmentation->Npx(),fSegmentation->Npy());

  Int_t nx = fSegmentation->Npx();
  Int_t ny = fSegmentation->Npy();

  //Int_t size=nx*ny;

  //printf("Size:%d\n",size);

  fTresh = new AliRICHTresholdMap(fSegmentation);

  //printf("Generating tresholds...\n");

  for(Int_t i=-nx/2;i<nx/2;i++)
    {
      for(Int_t j=-ny/2;j<ny/2;j++)
	{
	  Int_t pedestal = (Int_t)(gRandom->Gaus(50, 10));
	  //Int_t pedestal =0;
	  fTresh->SetHit(i,j,pedestal);
	  //printf("Pad %d %d has pedestal %d.\n",i,j,pedestal);
	}
    }
      
}
