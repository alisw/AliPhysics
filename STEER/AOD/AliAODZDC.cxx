/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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
//     Class for AOD ZDC data
//     Author: Chiara Oppedisano
//     Chiara.Oppedisano@cern.ch March 2011
//-------------------------------------------------------------------------

#include <TMath.h>
#include "AliAODZDC.h"

ClassImp(AliAODZDC)

AliAODZDC::AliAODZDC() :
  AliVZDC(),
  fZNCEnergy(-999.),
  fZNAEnergy(-999.),
  fZPCEnergy(-999.),
  fZPAEnergy(-999.),
  fZEM1Energy(0.),
  fZEM2Energy(0.),
  fZDCParticipants(0),
  fZDCPartSideA(0),
  fZDCPartSideC(0),
  fImpactParameter(0),
  fImpactParamSideA(0),
  fImpactParamSideC(0),
  fZDCTDCSum(0),	 
  fZDCTDCDifference(0),
  fZNCTDC(-999.),
  fZNATDC(-999.),
  fZPCTDC(-999.),
  fZPATDC(-999.),
  fIsZNAfired(kFALSE),
  fIsZNCfired(kFALSE),
  fIsZPAfired(kFALSE),
  fIsZPCfired(kFALSE)
{
// Default constructor
  for(Int_t i=0; i<5; i++){
    fZNCTowerEnergy[i] = fZNATowerEnergy[i] = 0.;
    fZPCTowerEnergy[i] = fZPATowerEnergy[i] = 0.;
    fZNCTowerEnergyLR[i] = fZNATowerEnergyLR[i] = 0.;
    fZPCTowerEnergyLR[i] = fZPATowerEnergyLR[i] = 0.;
  }
  for(Int_t i=0; i<4; i++){
     fZNCTDCm[i] =  fZNATDCm[i] =  fZPCTDCm[i] = fZPATDCm[i] = -999.;
  }
}

//__________________________________________________________________________
AliAODZDC::AliAODZDC(const AliAODZDC &zdcAOD) :
  AliVZDC(zdcAOD),
  fZNCEnergy(zdcAOD.fZNCEnergy),
  fZNAEnergy(zdcAOD.fZNAEnergy),
  fZPCEnergy(zdcAOD.fZPCEnergy),
  fZPAEnergy(zdcAOD.fZPAEnergy),
  fZEM1Energy(zdcAOD.fZEM1Energy),
  fZEM2Energy(zdcAOD.fZEM2Energy),
  fZDCParticipants(zdcAOD.fZDCParticipants),
  fZDCPartSideA(zdcAOD.fZDCPartSideA),
  fZDCPartSideC(zdcAOD.fZDCPartSideC),
  fImpactParameter(zdcAOD.fImpactParameter),
  fImpactParamSideA(zdcAOD.fImpactParamSideA),
  fImpactParamSideC(zdcAOD.fImpactParamSideC),
  fZDCTDCSum(zdcAOD.fZDCTDCSum),	 
  fZDCTDCDifference(zdcAOD.fZDCTDCDifference),
  fZNCTDC(zdcAOD.fZNCTDC),
  fZNATDC(zdcAOD.fZNATDC),
  fZPCTDC(zdcAOD.fZPCTDC),
  fZPATDC(zdcAOD.fZPATDC),
  fIsZNAfired(zdcAOD.fIsZNAfired),
  fIsZNCfired(zdcAOD.fIsZNCfired),
  fIsZPAfired(zdcAOD.fIsZPAfired),
  fIsZPCfired(zdcAOD.fIsZPCfired)

{
// Constructor
  for(Int_t i=0; i<5; i++){
    fZNCTowerEnergy[i] = zdcAOD.fZNCTowerEnergy[i];
    fZNATowerEnergy[i] = zdcAOD.fZNATowerEnergy[i];
    fZPCTowerEnergy[i] = zdcAOD.fZPCTowerEnergy[i];
    fZPATowerEnergy[i] = zdcAOD.fZPATowerEnergy[i];
    fZNCTowerEnergyLR[i] = zdcAOD.fZNCTowerEnergyLR[i];
    fZNATowerEnergyLR[i] = zdcAOD.fZNATowerEnergyLR[i];
    fZPCTowerEnergyLR[i] = zdcAOD.fZPCTowerEnergyLR[i];
    fZPATowerEnergyLR[i] = zdcAOD.fZPATowerEnergyLR[i];
  }
  for(Int_t i=0; i<4; i++){
     fZNCTDCm[i] =  zdcAOD.fZNCTDCm[i];
     fZNATDCm[i] =  zdcAOD.fZNATDCm[i];
     fZPCTDCm[i] =  zdcAOD.fZPCTDCm[i];
     fZPATDCm[i] =  zdcAOD.fZPATDCm[i];
  }
}

//__________________________________________________________________________
AliAODZDC& AliAODZDC::operator=(const AliAODZDC& zdcAOD)
{
  // Assignment operator
  //
  if(this!=&zdcAOD) {
    TObject::operator=(zdcAOD);
    fZNCEnergy  = zdcAOD.fZNCEnergy;
    fZNAEnergy  = zdcAOD.fZNAEnergy;
    fZPCEnergy  = zdcAOD.fZPCEnergy;
    fZPAEnergy  = zdcAOD.fZPAEnergy;
    fZEM1Energy = zdcAOD.fZEM1Energy;
    fZEM2Energy = zdcAOD.fZEM2Energy;
    for(Int_t i=0; i<5; i++){
       fZNCTowerEnergy[i] = zdcAOD.fZNCTowerEnergy[i];
       fZNATowerEnergy[i] = zdcAOD.fZNATowerEnergy[i];
       fZPCTowerEnergy[i] = zdcAOD.fZPCTowerEnergy[i];
       fZPATowerEnergy[i] = zdcAOD.fZPATowerEnergy[i];
       fZNCTowerEnergyLR[i] = zdcAOD.fZNCTowerEnergyLR[i];
       fZNATowerEnergyLR[i] = zdcAOD.fZNATowerEnergyLR[i];
       fZPCTowerEnergyLR[i] = zdcAOD.fZPCTowerEnergyLR[i];
       fZPATowerEnergyLR[i] = zdcAOD.fZPATowerEnergyLR[i];
    }
    //
    fZDCParticipants = zdcAOD.fZDCParticipants;
    fZDCPartSideA = zdcAOD.fZDCPartSideA;
    fZDCPartSideC = zdcAOD.fZDCPartSideC;
    fImpactParameter = zdcAOD.fImpactParameter;
    fImpactParamSideA = zdcAOD.fImpactParamSideA;
    fImpactParamSideC = zdcAOD.fImpactParamSideC;
    //
    fZDCTDCSum = zdcAOD.fZDCTDCSum;    
    fZDCTDCDifference = zdcAOD.fZDCTDCDifference;
    fZNCTDC = zdcAOD.fZNCTDC;
    fZNATDC = zdcAOD.fZNATDC;
    fZPCTDC = zdcAOD.fZPCTDC;
    fZPATDC = zdcAOD.fZPATDC;
    for(Int_t i=0; i<4; i++){
       fZNCTDCm[i] =  zdcAOD.fZNCTDCm[i];
       fZNATDCm[i] =  zdcAOD.fZNATDCm[i];
       fZPCTDCm[i] =  zdcAOD.fZPCTDCm[i];
       fZPATDCm[i] =  zdcAOD.fZPATDCm[i];
    }
    fIsZNAfired = zdcAOD.fIsZNAfired;
    fIsZNCfired = zdcAOD.fIsZNCfired;
    fIsZPAfired = zdcAOD.fIsZPAfired;
    fIsZPCfired = zdcAOD.fIsZPCfired;
  } 
  return *this;
}

//______________________________________________________________________________
void  AliAODZDC::SetZNCTowers(const Double_t value[5], const Double_t valueLR[5])
{
  // Sets ZNC towers
  for(Int_t i=0; i<5; i++){
    fZNCTowerEnergy[i] = value[i];
    fZNCTowerEnergyLR[i] = valueLR[i];
  }
}

//______________________________________________________________________________
void  AliAODZDC::SetZNATowers(const Double_t value[5], const Double_t valueLR[5])
{
  // Sets ZNA towers
  for(Int_t i=0; i<5; i++){
    fZNATowerEnergy[i] = value[i];
    fZNATowerEnergyLR[i] = valueLR[i];
  }
}

//______________________________________________________________________________
void  AliAODZDC::SetZPCTowers(const Double_t value[5], const Double_t valueLR[5])
{
  // Sets ZPC towers
  for(Int_t i=0; i<5; i++){
    fZPCTowerEnergy[i] = value[i];
    fZPCTowerEnergyLR[i] = valueLR[i];
  }
}

//______________________________________________________________________________
void  AliAODZDC::SetZPATowers(const Double_t value[5], const Double_t valueLR[5])
{
  // Sets ZPA towers
  for(Int_t i=0; i<5; i++){
    fZPATowerEnergy[i] = value[i];
    fZPATowerEnergyLR[i] = valueLR[i];
  }
}

//______________________________________________________________________________
Bool_t AliAODZDC::GetZNCentroidInPbPb(Float_t beamEne, Double_t centrZNC[2], Double_t centrZNA[2]) 
{
  // Provides coordinates of centroid over ZN (side C) front face in PbPb
   if(beamEne==0){
    printf(" ZDC centroid in PbPb can't be calculated with E_beam = 0 !!!\n");
    //for(Int_t jj=0; jj<2; jj++) fZNCCentrCoord[jj] = 999.;
    return kFALSE;
  }

  const Float_t x[4] = {-1.75, 1.75, -1.75, 1.75};
  const Float_t y[4] = {-1.75, -1.75, 1.75, 1.75};
  const Float_t alpha=0.395;
  Float_t numXZNC=0., numYZNC=0., denZNC=0., cZNC, wZNC; 
  Float_t numXZNA=0., numYZNA=0., denZNA=0., cZNA, wZNA; 
  Float_t zncEnergy=0., znaEnergy=0.;
  //
  for(Int_t i=0; i<5; i++){
    zncEnergy += fZNCTowerEnergy[i];
    znaEnergy += fZNATowerEnergy[i];
  }
  for(Int_t i=0; i<4; i++){
    if(fZNCTowerEnergy[i+1]>0.) {
      wZNC = TMath::Power(fZNCTowerEnergy[i+1], alpha);
      numXZNC += x[i]*wZNC;
      numYZNC += y[i]*wZNC;
      denZNC += wZNC;
    }
    if(fZNATowerEnergy[i+1]>0.) {
      wZNA = TMath::Power(fZNATowerEnergy[i+1], alpha);
      numXZNA += x[i]*wZNA;
      numYZNA += y[i]*wZNA;
      denZNA += wZNA;
    }
  }
  //
  if(denZNC!=0){
    Float_t nSpecnC = zncEnergy/beamEne;
    cZNC = 1.89358-0.71262/(nSpecnC+0.71789);
    centrZNC[0] = cZNC*numXZNC/denZNC;
    centrZNC[1] = cZNC*numYZNC/denZNC;
  } 
  else{
    centrZNC[0] = centrZNC[1] = 999.;
  }
  if(denZNA!=0){
    Float_t nSpecnA = znaEnergy/beamEne;
    cZNA = 1.89358-0.71262/(nSpecnA+0.71789);
    centrZNA[0] = cZNA*numXZNA/denZNA;
    centrZNA[1] = cZNA*numYZNA/denZNA;
  } 
  else{
    centrZNA[0] = centrZNA[1] = 999.;
  }
 
  
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAODZDC::GetZNCentroidInpp(Double_t centrZNC[2], Double_t centrZNA[2]) 
{
  // Provides coordinates of centroid over ZN (side C) front face in pp
  const Float_t x[4] = {-1.75, 1.75, -1.75, 1.75};
  const Float_t y[4] = {-1.75, -1.75, 1.75, 1.75};
  const Float_t alpha=0.5;
  Float_t numXZNC=0., numYZNC=0., denZNC=0., wZNC; 
  Float_t numXZNA=0., numYZNA=0., denZNA=0., wZNA; 
  //
  for(Int_t i=0; i<4; i++){
    if(fZNCTowerEnergy[i+1]>0.) {
      wZNC = TMath::Power(fZNCTowerEnergy[i+1], alpha);
      numXZNC += x[i]*wZNC;
      numYZNC += y[i]*wZNC;
      denZNC += wZNC;
    }
    if(fZNATowerEnergy[i+1]>0.) {
      wZNA = TMath::Power(fZNATowerEnergy[i+1], alpha);
      numXZNA += x[i]*wZNA;
      numYZNA += y[i]*wZNA;
      denZNA += wZNA;
    }
  }
  //
  if(denZNC!=0){
    centrZNC[0] = numXZNC/denZNC;
    centrZNC[1] = numYZNC/denZNC;
  } 
  else{
    centrZNC[0] = centrZNC[1] = 999.;
  }
  if(denZNA!=0){
    centrZNA[0] = numXZNA/denZNA;
    centrZNA[1] = numYZNA/denZNA;
  } 
  else{
    centrZNA[0] = centrZNA[1] = 999.;
  }
  
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAODZDC::GetTDCSum(Float_t sum[4]) 
{
  // Provides value(s) != -999 if both ZN are fired
  for(int i=0; i<4; i++) sum[i] = -999.;
  int ind=0;
  if(IsZNANDfired()){
    for(Int_t i=0;i<4;++i){
      if(fZNCTDCm[i]>-999){
        for(Int_t j=0;j<4;++j){
          if(fZNATDCm[j]>-999){
	    sum[ind] = fZNCTDCm[i]+fZNATDCm[j];
	    ind++;
	  } // if ZNA[j] is hit
        } // Loop over ZNA TDC hits
      } // if ZNC[i] is hit
    } // Loop over ZNC TDC hits
    return kTRUE;
  } // ZN AND
  else return kFALSE;
}

//______________________________________________________________________________
Bool_t AliAODZDC::GetTDCDiff(Float_t diff[4]) 
{
  // Provides value(s) != -999 if both ZN are fired
  for(int i=0; i<4; i++) diff[i] = -999.;
  int ind=0;
  if(IsZNANDfired()){
    for(Int_t i=0;i<4;++i){
      if(fZNCTDCm[i]>-999){
        for(Int_t j=0;j<4;++j){
          if(fZNATDCm[j]>-999){
	    diff[ind] = fZNCTDCm[i]-fZNATDCm[j];
	    ind++;
	  } // if ZNA[j] is hit
        } // Loop over ZNA TDC hits
      } // if ZNC[i] is hit
    } // Loop over ZNC TDC hits
    return kTRUE;
  } // ZN AND
  else return kFALSE;
}
