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

#include "AliRICHChamber.h"
#include "AliRICHConst.h" //for kR2d
#include "AliRICHParam.h"
#include <TRandom.h>
#include <TRotMatrix.h>
#include "AliRICHTresholdMap.h"
#include "AliSegmentation.h"
#include "AliRICHSegmentationV0.h"
#include "AliRICHGeometry.h"
#include "AliRICHResponse.h"

ClassImp(AliRICHChamber)	
//______________________________________________________________________________    
AliRICHChamber::AliRICHChamber() 
{//default ctor
  fpParam=0;    
  fpRotMatrix=0;
  
    fSegmentation = 0;
    fResponse = 0;
    fGeometry = 0;
    fTresh = 0;
    for(Int_t i=0; i<50; ++i) fIndexMap[i] = 0;
}
//______________________________________________________________________________
AliRICHChamber::AliRICHChamber(Int_t iModuleN,AliRICHParam *pParam)
{//named ctor. Defines all geometry parameters for the given module.
 //   6 7 |----> X
 // 3 4 5 | 
 // 1 2   V Z  
  SetCenter(0,pParam->Offset()-pParam->GapThickness()/2,0);//put to up position   
  switch(iModuleN){
    case 1:
      RotateX(-pParam->AngleYZ());   
      RotateZ(-pParam->AngleXY());      
      fName="RICHc1";fTitle="RICH chamber 1";
      break;      
    case 2:
      RotateZ(-pParam->AngleXY());      
      fName="RICHc2";fTitle="RICH chamber 2";
      break;      
    case 3:
      RotateX(-pParam->AngleYZ());
      fName="RICHc3";fTitle="RICH chamber 3";
      break;      
    case 4:          
      fName="RICHc4";fTitle="RICH chamber 4";  //no turns
      break;      
    case 5:
      RotateX(pParam->AngleYZ());
      fName="RICHc5";fTitle="RICH chamber 5";
      break;      
    case 6:
      RotateZ(pParam->AngleXY());      
      fName="RICHc6";fTitle="RICH chamber 6";
      break;      
    case 7:
      RotateX(pParam->AngleYZ());            
      RotateZ(pParam->AngleXY());      
      fName="RICHc7";fTitle="RICH chamber 7";
      break;      
    default:
      Fatal("named ctor","Wrong chamber number %i, check muster class ctor",iModuleN);
  }//switch(iModuleN)
  RotateZ(pParam->AngleRot());//apply common rotation  
  fpRotMatrix=new TRotMatrix("rot"+fName,"rot"+fName, Rot().ThetaX()*kR2d, Rot().PhiX()*kR2d,
                                                      Rot().ThetaY()*kR2d, Rot().PhiY()*kR2d,
                                                      Rot().ThetaZ()*kR2d, Rot().PhiZ()*kR2d);
  fpParam=pParam;
  fX=fCenterV3.X();fY=fCenterV3.Y();fZ=fCenterV3.Z();
}
//______________________________________________________________________________

void AliRICHChamber::LocaltoGlobal(Float_t local[3],Float_t global[3])
{//Local coordinates to global coordinates transformation

    Double_t *pMatrix;
    pMatrix =  fpRotMatrix->GetMatrix();
    global[0]=local[0]*pMatrix[0]+local[1]*pMatrix[3]+local[2]*pMatrix[6];
    global[1]=local[0]*pMatrix[1]+local[1]*pMatrix[4]+local[2]*pMatrix[7];
    global[2]=local[0]*pMatrix[2]+local[1]*pMatrix[5]+local[2]*pMatrix[8];
    global[0]+=fX;
    global[1]+=fY;
    global[2]+=fZ;
}

void AliRICHChamber::GlobaltoLocal(Float_t global[3],Float_t local[3])
{// Global coordinates to local coordinates transformation
    TMatrix matrixCopy(3,3);
    Double_t *pMatrixOrig = fpRotMatrix->GetMatrix();
    for(Int_t i=0;i<3;i++)
      {
	for(Int_t j=0;j<3;j++)
	  matrixCopy(j,i)=pMatrixOrig[j+3*i];
      }
    matrixCopy.Invert();
    local[0] = global[0] - fX;
    local[1] = global[1] - fY;
    local[2] = global[2] - fZ;
    local[0]=local[0]*matrixCopy(0,0)+local[1]*matrixCopy(0,1)+local[2]*matrixCopy(0,2);
    local[1]=local[0]*matrixCopy(1,0)+local[1]*matrixCopy(1,1)+local[2]*matrixCopy(1,2);
    local[2]=local[0]*matrixCopy(2,0)+local[1]*matrixCopy(2,1)+local[2]*matrixCopy(2,2);
} 

void AliRICHChamber::DisIntegration(Float_t eloss, Float_t xhit, Float_t yhit,
				    Int_t& iNpads,Float_t cluster[5][500],ResponseType res) 
{//Generates pad hits (simulated cluster) using the segmentation and the response model

  Float_t local[3],global[3];
// Width of the integration area
  Float_t dx=(fResponse->SigmaIntegration())*(fResponse->ChargeSpreadX());
  Float_t dy=(fResponse->SigmaIntegration())*(fResponse->ChargeSpreadY());
// Get pulse height from energy loss and generate feedback photons
  Float_t qtot=0;
  local[0]=xhit;
//z-position of the wires relative to the RICH mother volume (2 mm before CsI) old value: 6.076 ???????
  local[1]=1.276 + fGeometry->GetGapThickness()/2  - .2;
  local[2]=yhit;

  LocaltoGlobal(local,global);



//To calculate wire sag, the origin of y-position must be the middle of the photcathode
  AliRICHSegmentationV0* segmentation = (AliRICHSegmentationV0*) GetSegmentationModel();
  Float_t newy;
  if (yhit>0)
    newy = yhit - segmentation->GetPadPlaneLength()/2;
  else
    newy = yhit + segmentation->GetPadPlaneLength()/2;

  if(res==kMip){
    qtot = fResponse->IntPH(eloss, newy);
    fResponse->FeedBackPhotons(global,qtot);
  }else if(res==kPhoton){
    qtot = fResponse->IntPH(newy);
    fResponse->FeedBackPhotons(global,qtot);
  }

    // Loop Over Pads

  Float_t qcheck=0, qp=0;

  iNpads=0;
  for(fSegmentation->FirstPad(xhit, yhit, 0, dx, dy);
      fSegmentation->MorePads();
      fSegmentation->NextPad()) {
    qp= fResponse->IntXY(fSegmentation);
    qp= TMath::Abs(qp);
    if(qp >1.e-4){
      qcheck+=qp;
      cluster[0][iNpads]=qp*qtot;// --- store signal information
      cluster[1][iNpads]=fSegmentation->Ix();
      cluster[2][iNpads]=fSegmentation->Iy();
      cluster[3][iNpads]=fSegmentation->ISector();
      iNpads++;
    }
  }//pad loop
}//void AliRICHChamber::DisIntegration(...
//__________________________________________________________________________________________________
void AliRICHChamber::GenerateTresholds()
{//Generates random treshold charges for all pads
  Int_t nx = fSegmentation->Npx();
  Int_t ny = fSegmentation->Npy();

  fTresh = new AliRICHTresholdMap(fSegmentation);
  for(Int_t i=-nx/2;i<nx/2;i++){
    for(Int_t j=-ny/2;j<ny/2;j++){
      Int_t pedestal = (Int_t)(gRandom->Gaus(50, 10));
      fTresh->SetHit(i,j,pedestal);
    }
  }      
}//void AliRICHChamber::GenerateTresholds()
//__________________________________________________________________________________________________
void AliRICHChamber::Print(Option_t *) const
{
  printf("%s r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f  %6.2f,%6.2f %6.2f,%6.2f %6.2f,%6.2f\n",fName.Data(),
                     Rho(), ThetaD(),PhiD(),   X(),    Y(),    Z(),
                     ThetaXd(),PhiXd(),ThetaYd(),PhiYd(),ThetaZd(),PhiZd());
}//void AliRICHChamber::Print(Option_t *option)const
//__________________________________________________________________________________________________
void AliRICHChamber::SetChamberTransform(Float_t x,Float_t y,Float_t z,TRotMatrix *pRotMatrix) 
{
  SetCenter(x,y,z);
  fpRotMatrix=pRotMatrix;    
}
