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

/* $Id: AliEMCALRecoUtils.cxx 33808 2009-07-15 09:48:08Z gconesab $ */

///////////////////////////////////////////////////////////////////////////////
//
// Class AliEMCALRecoUtils
// Some utilities to recalculate the cluster position or energy linearity
//
//
// Author:  Gustavo Conesa (LPSC- Grenoble) 
///////////////////////////////////////////////////////////////////////////////

// --- standard c ---

// standard C++ includes
//#include <Riostream.h>

// ROOT includes

// STEER includes
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeoUtils.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliLog.h"

ClassImp(AliEMCALRecoUtils)
  
//______________________________________________
AliEMCALRecoUtils::AliEMCALRecoUtils():
  fNonLinearityFunction (kPi0GammaGamma)
{
//
  // Constructor.
  // Initialize all constant values which have to be used
  // during Reco algorithm execution
  //
  
  for(Int_t i = 0; i < 15 ; i++) {fMisalTransShift[i] = 0.; fMisalRotShift[i] = 0.; }
  for(Int_t i = 0; i < 6  ; i++) fNonLinearityParams[i] = 0.; 
  //By default kPi0GammaGamma case
  fNonLinearityParams[0] = 0.1457/0.1349766/1.038;
  fNonLinearityParams[1] = -0.02024/0.1349766/1.038;
  fNonLinearityParams[2] = 1.046;
  
}

//______________________________________________________________________
AliEMCALRecoUtils::AliEMCALRecoUtils(const AliEMCALRecoUtils & reco) 
: TNamed(reco), fNonLinearityFunction(reco.fNonLinearityFunction)
{
  //Copy ctor
  
  for(Int_t i = 0; i < 15 ; i++) {fMisalRotShift[i] = reco.fMisalRotShift[i]; fMisalTransShift[i] = reco.fMisalTransShift[i]; } 
  for(Int_t i = 0; i < 6  ; i++) fNonLinearityParams[i] = reco.fNonLinearityParams[i]; 
}


//______________________________________________________________________
AliEMCALRecoUtils & AliEMCALRecoUtils::operator = (const AliEMCALRecoUtils & reco) 
{
  //Assignment operator
  
  if(this == &reco)return *this;
  ((TNamed *)this)->operator=(reco);

  fNonLinearityFunction = reco.fNonLinearityFunction;
  for(Int_t i = 0; i < 15 ; i++) {fMisalTransShift[i] = reco.fMisalTransShift[i]; fMisalRotShift[i] = reco.fMisalRotShift[i];}
  for(Int_t i = 0; i < 6  ; i++) fNonLinearityParams[i] = reco.fNonLinearityParams[i]; 
  
  return *this;
}


//__________________________________________________
Float_t AliEMCALRecoUtils::CorrectClusterEnergyLinearity(AliVCluster* cluster){
// Correct cluster energy from non linearity functions
  Float_t energy = cluster->E();
  
  switch (fNonLinearityFunction) {
      
    case kPi0MC:
      //Non-Linearity correction (from MC with function ([0]*exp(-[1]/E))+(([2]/([3]*2.*TMath::Pi())*exp(-(E-[4])^2/(2.*[3]^2)))))
      //Double_t par0 = 1.001;
      //Double_t par1 = -0.01264;
      //Double_t par2 = -0.03632;
      //Double_t par3 = 0.1798;
      //Double_t par4 = -0.522;
       energy /= (fNonLinearityParams[0]*exp(-fNonLinearityParams[1]/energy))+
                  ((fNonLinearityParams[2]/(fNonLinearityParams[3]*2.*TMath::Pi())*
                    exp(-(energy-fNonLinearityParams[4])*(energy-fNonLinearityParams[4])/(2.*fNonLinearityParams[3]*fNonLinearityParams[3]))));
      break;
      
    case kPi0GammaGamma:

      //Non-Linearity correction (from Olga Data with function p0+p1*exp(-p2*E))
      //Double_t par0 = 0.1457;
      //Double_t par1 = -0.02024;
      //Double_t par2 = 1.046;
      energy /= (fNonLinearityParams[0]+fNonLinearityParams[1]*exp(-fNonLinearityParams[2]*energy)); //Olga function
      break;
      
    case kPi0GammaConversion:
      
      //Non-Linearity correction (Nicolas from Dimitri Data with function C*[1-a*exp(-b*E)])
      //Double_t C = 0.139393/0.1349766;
      //Double_t a = 0.0566186;
      //Double_t b = 0.982133;
      energy /= fNonLinearityParams[0]*(1-fNonLinearityParams[1]*exp(-fNonLinearityParams[2]*energy));
      
      break;
      
    case kNoCorrection:
      AliDebug(2,"No correction on the energy\n");
      break;
      
  }
  
  return energy;

}

//__________________________________________________
void AliEMCALRecoUtils::GetMaxEnergyCell(AliEMCALGeoUtils *geom, AliVCaloCells* cells, AliVCluster* clu, Int_t & absId,  Int_t& iSupMod, Int_t& ieta, Int_t& iphi)
{
  //For a given CaloCluster gets the absId of the cell 
  //with maximum energy deposit.
  
  Double_t eMax        = -1.;
  Double_t eCell       = -1.;
  Int_t    cellAbsId   = -1 ;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
	
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) {
     cellAbsId = clu->GetCellAbsId(iDig);
     eCell     = cells->GetCellAmplitude(cellAbsId);
     if(eCell > eMax)  { 
      eMax  = eCell; 
      absId = cellAbsId;
      //printf("\t new max: cell %d, e %f, ecell %f\n",maxId, eMax,eCell);
    }
  }// cell loop
  
  //Get from the absid the supermodule, tower and eta/phi numbers
  geom->GetCellIndex(absId,iSupMod,iTower,iIphi,iIeta); 
  //Gives SuperModule and Tower numbers
  geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
                                         iIphi, iIeta,iphi,ieta);    
  
}

//__________________________________________________
void AliEMCALRecoUtils::RecalculateClusterPosition(AliEMCALGeoUtils *geom, AliVCaloCells* cells, AliVCluster* clu, Int_t iParticle)
{
  //For a given CaloCluster recalculates the position for a given set of misalignment shifts and puts it again in the CaloCluster.
  
  Double_t eMax       = -1.;
  Double_t eCell      = -1.;
  Int_t    cellAbsId  = -1;
	
  Int_t maxId   = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
	Int_t iSupMod = -1;
  Int_t iphi = -1, ieta =-1;
  
  Float_t clEnergy = clu->E(); //Energy already recalibrated previously.
  Float_t weight = 0., weightedCol = 0., weightedRow = 0., totalWeight=0.;
  Bool_t  areInSameSM = kTRUE; //exclude clusters with cells in different SMs for now
  Int_t   startingSM = -1;
  
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) {
    cellAbsId = clu->GetCellAbsId(iDig);
    geom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta); 
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);			
    if     (iDig==0)  startingSM = iSupMod;
    else if(iSupMod != startingSM) areInSameSM = kFALSE;
    
    eCell  = cells->GetCellAmplitude(cellAbsId);
    
    weight = TMath::Log(eCell/clEnergy) + 4;
    if(weight < 0) weight = 0;
    totalWeight += weight;
    weightedCol += ieta*weight;
    weightedRow += iphi*weight;
    
    //printf("Max cell? cell %d, amplitude org %f, fraction %f, recalibration %f, amplitude new %f \n",cellAbsId, cells->GetCellAmplitude(cellAbsId), fraction, recalFactor, eCell) ;
    
    if(eCell > eMax)  { 
      eMax  = eCell; 
      maxId = cellAbsId;
      //printf("\t new max: cell %d, e %f, ecell %f\n",maxId, eMax,eCell);
    }
  }// cell loop
  
  //Get from the absid the supermodule, tower and eta/phi numbers  
  geom->GetCellIndex(maxId,iSupMod,iTower,iIphi,iIeta); 
  //Gives SuperModule and Tower numbers
  geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
                                         iIphi, iIeta,iphi,ieta);  
  
  Float_t xyzNew[]={0.,0.,0.};
  if(areInSameSM == kTRUE) {
    //printf("In Same SM\n");
    weightedCol = weightedCol/totalWeight;
    weightedRow = weightedRow/totalWeight;
    geom->RecalculateTowerPosition(weightedRow, weightedCol, iSupMod, clEnergy, iParticle, fMisalTransShift, fMisalRotShift, xyzNew); 
  }
  else {
    //printf("In Different SM\n");
    geom->RecalculateTowerPosition(iphi,        ieta,        iSupMod, clEnergy, iParticle, fMisalTransShift, fMisalRotShift, xyzNew); 
  }
  clu->SetPosition(xyzNew);

  
  //printf("\t Max : cell %d, iSupMod %d, ieta %d, iphi %d \n",maxId,iSupMod, ieta,iphi);
  
}

//__________________________________________________
void AliEMCALRecoUtils::Print(const Option_t *) const 
{
  // Print Parameters
  
  printf("AliEMCALRecoUtils Settings: \n");
  printf("Misalignment shifts\n");
  for(Int_t i=0; i<5; i++) printf("\t sector %d, traslation (x,y,z)=(%f,%f,%f), rotation (x,y,z)=(%f,%f,%f)\n",i, 
                                  fMisalTransShift[i*3],fMisalTransShift[i*3+1],fMisalTransShift[i*3+2],
                                  fMisalRotShift[i*3],  fMisalRotShift[i*3+1],  fMisalRotShift[i*3+2]   );
  printf("Non linearity function %d, parameters:\n", fNonLinearityFunction);
  for(Int_t i=0; i<6; i++) printf("param[%d]=%f\n",i, fNonLinearityParams[i]);
    
}
