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
  Revision 1.1  2000/04/19 12:57:20  morsch
  Newly structured and updated version (JB, AM)

*/


#include "AliRICHChamber.h"
#include "AliRun.h"
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TRandom.h>

ClassImp(AliRICHChamber)	
    
AliRICHChamber::AliRICHChamber() 
{
    fSegmentation = 0;
    fResponse= 0;
    fGeometry= 0;
    frMin=0.1;
    frMax=140;
    fnsec=1;
}

//  
//  Get reference to response model
AliRICHResponse* AliRICHChamber::GetResponseModel()
{
    return fResponse;
}

// Configure response model
void   AliRICHChamber::ResponseModel(AliRICHResponse* thisResponse)
{
    fResponse=thisResponse;
}

void AliRICHChamber::Init()
{
    fSegmentation->Init(this);
}

void AliRICHChamber::LocaltoGlobal(Float_t pos[3],Float_t Globalpos[3])
{

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
				    Int_t& nnew,Float_t newclust[6][500],Response_t res) 
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

    Int_t Nfp=0;
    
    if (res==mip) {
	qtot = fResponse->IntPH(eloss);
	Nfp  = fResponse->FeedBackPhotons(global,qtot);
    } else if (res==cerenkov) {
	qtot = fResponse->IntPH();
	Nfp  = fResponse->FeedBackPhotons(global,qtot);
    }

    //printf("Feedbacks:%d\n",Nfp);
    
    //
    // Loop Over Pads
    
    Float_t qcheck=0, qp=0;
    
    nnew=0;
    for (Int_t i=1; i<=fnsec; i++) {
	qcheck=0;
	for (fSegmentation->FirstPad(xhit, yhit, dx, dy); 
	     fSegmentation->MorePads(); 
	     fSegmentation->NextPad()) 
	{
	    qp= fResponse->IntXY(fSegmentation);
	    qp= TMath::Abs(qp);

	    //printf("Qp:%f\n",qp);

	    if (qp > 1.e-4) {
		qcheck+=qp;
		//
		// --- store signal information
		newclust[0][nnew]=qtot;
		newclust[1][nnew]=fSegmentation->Ix();
		newclust[2][nnew]=fSegmentation->Iy();
		newclust[3][nnew]=qp * qtot;
		newclust[4][nnew]=fSegmentation->ISector();
		newclust[5][nnew]=(Float_t) i;
		nnew++;	
		//printf("Newcluster:%d\n",i);
	    }
	} // Pad loop
    } // Cathode plane loop
    //if (fSegmentation->ISector()==2)
      //printf("Nnew:%d\n\n\n\n",nnew);
}




