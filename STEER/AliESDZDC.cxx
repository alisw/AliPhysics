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


//-------------------------------------------------------------------------
//                      Implementation of   Class AliESDZDC
//   This is a class that summarizes the ZDC data
//   for the ESD   
//   Origin: Christian Klein-Boesing, CERN, Christian.Klein-Boesing@cern.ch 
//-------------------------------------------------------------------------



#include "AliESDZDC.h"

ClassImp(AliESDZDC)

//______________________________________________________________________________
AliESDZDC::AliESDZDC() :
  TObject(),
  fZDCN1Energy(0),
  fZDCP1Energy(0),
  fZDCN2Energy(0),
  fZDCP2Energy(0),
  fZDCEMEnergy(0),
  fZDCEMEnergy1(0),
  fZDCParticipants(0),
  fZDCPartSideA(0),
  fZDCPartSideC(0),
  fImpactParameter(0),
  fImpactParamSideA(0),
  fImpactParamSideC(0),
  fESDQuality(0)
{
  for(Int_t i=0; i<5; i++){
    fZN1TowerEnergy[i] = fZN2TowerEnergy[i] = 0.;
    fZP1TowerEnergy[i] = fZP2TowerEnergy[i] = 0.;
    fZN1TowerEnergyLR[i] = fZN2TowerEnergyLR[i] = 0.;
    fZP1TowerEnergyLR[i] = fZP2TowerEnergyLR[i] = 0.;
  }
  for(Int_t i=0; i<2; i++){
    fZNACentrCoord[i] = fZNCCentrCoord[i] = 0.;
  }
  for(Int_t i=0; i<32; i++) fVMEScaler[i]=0;
}

//______________________________________________________________________________
AliESDZDC::AliESDZDC(const AliESDZDC& zdc) :
  TObject(zdc),
  fZDCN1Energy(zdc.fZDCN1Energy),
  fZDCP1Energy(zdc.fZDCP1Energy),
  fZDCN2Energy(zdc.fZDCN2Energy),
  fZDCP2Energy(zdc.fZDCP2Energy),
  fZDCEMEnergy(zdc.fZDCEMEnergy),
  fZDCEMEnergy1(zdc.fZDCEMEnergy1),
  fZDCParticipants(zdc.fZDCParticipants),
  fZDCPartSideA(zdc.fZDCPartSideA),
  fZDCPartSideC(zdc.fZDCPartSideC),
  fImpactParameter(zdc.fImpactParameter),
  fImpactParamSideA(zdc.fImpactParamSideA),
  fImpactParamSideC(zdc.fImpactParamSideC),
  fESDQuality(zdc.fESDQuality)
{
  // copy constructor
  for(Int_t i=0; i<5; i++){
     fZN1TowerEnergy[i] = zdc.fZN1TowerEnergy[i];
     fZN2TowerEnergy[i] = zdc.fZN2TowerEnergy[i];
     fZP1TowerEnergy[i] = zdc.fZP1TowerEnergy[i];
     fZP2TowerEnergy[i] = zdc.fZP2TowerEnergy[i];
     fZN1TowerEnergyLR[i] = zdc.fZN1TowerEnergyLR[i];
     fZN2TowerEnergyLR[i] = zdc.fZN2TowerEnergyLR[i];
     fZP1TowerEnergyLR[i] = zdc.fZP1TowerEnergyLR[i];
     fZP2TowerEnergyLR[i] = zdc.fZP2TowerEnergyLR[i];
  }
  for(Int_t i=0; i<2; i++){
    fZNACentrCoord[i] = zdc.fZNACentrCoord[i];
    fZNCCentrCoord[i] = zdc.fZNCCentrCoord[i];
  }
  for(Int_t i=0; i<32; i++) fVMEScaler[i] = zdc.fVMEScaler[i];
}

//______________________________________________________________________________
AliESDZDC& AliESDZDC::operator=(const AliESDZDC&zdc)
{
  // assigment operator
  if(this!=&zdc) {
    TObject::operator=(zdc);
    fZDCN1Energy = zdc.fZDCN1Energy;
    fZDCP1Energy = zdc.fZDCP1Energy;
    fZDCN2Energy = zdc.fZDCN2Energy;
    fZDCP2Energy = zdc.fZDCP2Energy;
    fZDCEMEnergy = zdc.fZDCEMEnergy;
    fZDCEMEnergy1 = zdc.fZDCEMEnergy1;
    for(Int_t i=0; i<5; i++){
       fZN1TowerEnergy[i] = zdc.fZN1TowerEnergy[i];
       fZN2TowerEnergy[i] = zdc.fZN2TowerEnergy[i];
       fZP1TowerEnergy[i] = zdc.fZP1TowerEnergy[i];
       fZP2TowerEnergy[i] = zdc.fZP2TowerEnergy[i];
       fZN1TowerEnergyLR[i] = zdc.fZN1TowerEnergyLR[i];
       fZN2TowerEnergyLR[i] = zdc.fZN2TowerEnergyLR[i];
       fZP1TowerEnergyLR[i] = zdc.fZP1TowerEnergyLR[i];
       fZP2TowerEnergyLR[i] = zdc.fZP2TowerEnergyLR[i];
    }
    //
    fZDCParticipants = zdc.fZDCParticipants;
    fZDCPartSideA = zdc.fZDCPartSideA;
    fZDCPartSideC = zdc.fZDCPartSideC;
    fImpactParameter = zdc.fImpactParameter;
    fImpactParamSideA = zdc.fImpactParamSideA;
    fImpactParamSideC = zdc.fImpactParamSideC;
    //
    for(Int_t i=0; i<2; i++){
         fZNACentrCoord[i] = zdc.fZNACentrCoord[i];
         fZNCCentrCoord[i] = zdc.fZNCCentrCoord[i];
    }
    //
    fESDQuality = zdc.fESDQuality;
    for(Int_t i=0; i<32; i++) fVMEScaler[i] = zdc.fVMEScaler[i];
  } 
  return *this;
}

//______________________________________________________________________________
void AliESDZDC::Copy(TObject &obj) const {
  
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if(this==&obj)return;
  AliESDZDC *robj = dynamic_cast<AliESDZDC*>(&obj);
  if(!robj)return; // not an AliESDZDC
  *robj = *this;

}


//______________________________________________________________________________
void AliESDZDC::Reset()
{
  // reset all data members
  fZDCN1Energy=0;
  fZDCP1Energy=0;
  fZDCN2Energy=0;
  fZDCP2Energy=0;
  fZDCEMEnergy=0;
  fZDCEMEnergy1=0;
  for(Int_t i=0; i<5; i++){
    fZN1TowerEnergy[i] = fZN2TowerEnergy[i] = 0.;
    fZP1TowerEnergy[i] = fZP2TowerEnergy[i] = 0.;
    fZN1TowerEnergyLR[i] = fZN2TowerEnergyLR[i] = 0.;
    fZP1TowerEnergyLR[i] = fZP2TowerEnergyLR[i] = 0.;
  }
  fZDCParticipants=0;  
  fZDCPartSideA=0;  
  fZDCPartSideC=0;  
  fImpactParameter=0;
  fImpactParamSideA=0;
  fImpactParamSideC=0;
  for(Int_t i=0; i<2; i++){
       fZNACentrCoord[i] = fZNCCentrCoord[i] = 0.;
  }
  fESDQuality=0;
  for(Int_t i=0; i<32; i++) fVMEScaler[i] = 0;
}

//______________________________________________________________________________
void AliESDZDC::Print(const Option_t *) const
{
  //  Print ESD for the ZDC
  printf("\n \t E_{ZNC} = %f TeV, E_{ZPC} = %f TeV, E_{ZNA} = %f TeV, E_{ZPA} = %f TeV,"
  " E_{ZEM} = %f GeV, Npart = %d, b = %1.2f fm\n",
  fZDCN1Energy/1000.,fZDCP1Energy/1000.,fZDCN2Energy/1000.,fZDCP2Energy/1000.,
  fZDCEMEnergy+fZDCEMEnergy1, fZDCParticipants,fImpactParameter);
  //
  printf(" ### fVMEScaler: \n");
  for(Int_t i=0; i<32; i++) printf("\t datum %d: %d \n",i,fVMEScaler[i]);
  printf("\n");
}

//______________________________________________________________________________
const Double_t * AliESDZDC::GetZNCCentroidInPbPb(Float_t beamEne) 
{
  // Provide coordinates of centroid over ZN (side C) front face
  const Float_t x[4] = {-1.75, 1.75, -1.75, 1.75};
  const Float_t y[4] = {-1.75, -1.75, 1.75, 1.75};
  Float_t numX=0., numY=0., den=0.;
  Float_t c, w; 
  const Float_t alpha=0.395;
  //
  for(Int_t i=0; i<4; i++)
    if(fZN1TowerEnergy[i+1]>0.) {
      w = TMath::Power(fZN1TowerEnergy[i+1], alpha);
      numX += x[i]*w;
      numY += y[i]*w;
      den += w;
    }
  if(den!=0){
    Float_t nSpecn = fZDCN1Energy/beamEne;
    c = 1.89358-0.71262/(nSpecn+0.71789);
    fZNCCentrCoord[0] = c*numX/den;
    fZNCCentrCoord[1] = c*numY/den;
  } else {
    fZNCCentrCoord[0] = fZNCCentrCoord[1] = 0;
  }
  return fZNCCentrCoord;
}

//______________________________________________________________________________
const Double_t * AliESDZDC::GetZNACentroidInPbPb(Float_t beamEne) 
{
  // Provide coordinates of centroid over ZN (side A) front face
  const Float_t x[4] = {-1.75, 1.75, -1.75, 1.75};
  const Float_t y[4] = {-1.75, -1.75, 1.75, 1.75};
  Float_t numX=0., numY=0., den=0.;
  Float_t c, w;
  const Float_t alpha=0.395;

  for(Int_t i=0; i<4; i++)
    if(fZN2TowerEnergy[i+1]>0.) {
      w = TMath::Power(fZN2TowerEnergy[i+1], alpha);
      numX += x[i]*w;
      numY += y[i]*w;
      den += w;
    }
  //
  if(den!=0){
    Float_t nSpecn = fZDCN2Energy/beamEne;
    c = 1.89358-0.71262/(nSpecn+0.71789);
    fZNACentrCoord[0] = c*numX/den;
    fZNACentrCoord[1] = c*numY/den;
  } else {
    fZNACentrCoord[0] = fZNACentrCoord[1] = 0;
  }
  return fZNACentrCoord;
}

//______________________________________________________________________________
const Double_t * AliESDZDC::GetZNCCentroidInpp() 
{
  // Provide coordinates of centroid over ZN (side C) front face
  const Float_t x[4] = {-1.75, 1.75, -1.75, 1.75};
  const Float_t y[4] = {-1.75, -1.75, 1.75, 1.75};
  Float_t numX=0., numY=0., den=0.;
  const Float_t alpha=0.5;
  Float_t w;
  //
  for(Int_t i=0; i<4; i++)
    if(fZN1TowerEnergy[i+1]>0.) {
      w = TMath::Power(fZN1TowerEnergy[i+1], alpha);
      numX += x[i]*w;
      numY += y[i]*w;
      den += w;
    }

  if(den!=0){
    fZNCCentrCoord[0] = numX/den;
    fZNCCentrCoord[1] = numY/den;
  } else {
    fZNCCentrCoord[0] = fZNCCentrCoord[1] = 0;
  }
  return fZNCCentrCoord;
}

//______________________________________________________________________________
const Double_t * AliESDZDC::GetZNACentroidInpp() 
{
  // Provide coordinates of centroid over ZN (side A) front face
  const Float_t x[4] = {-1.75, 1.75, -1.75, 1.75};
  const Float_t y[4] = {-1.75, -1.75, 1.75, 1.75};
  Float_t numX=0., numY=0., den=0.;
  const Float_t alpha=0.395;
  Float_t w;

  for(Int_t i=0; i<4; i++)
    if(fZN2TowerEnergy[i+1]>0.) {
      w = TMath::Power(fZN2TowerEnergy[i+1], alpha);
      numX += x[i]*w;
      numY += y[i]*w;
      den += w;
    }
  //
  if(den!=0){
    fZNACentrCoord[0] = numX/den;
    fZNACentrCoord[1] = numY/den;
  } else {
    fZNACentrCoord[0] = fZNACentrCoord[1] = 0;
  }
    
  return fZNACentrCoord;
}
