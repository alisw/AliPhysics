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

//_________________________________________________________________________
//  AliAODCaloCluster extension for EMCAL to recalculate cluster 
//  parameters in case of recalibration.
//  Copy-paste from methods in AliEMCALRecPoint.
//*--
//*-- Author: Dmitri Peressounko (RRC KI) for PHOS
//*-- Adapted for EMCAL: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include "TVector3.h"
#include "TMath.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h" 
#include "AliEMCALGeometry.h" 
#include "AliEMCALPID.h" 
#include "AliEMCALAodCluster.h" 
#include "AliEMCALCalibData.h"
#include "AliAODCaloCells.h"

ClassImp(AliEMCALAodCluster)

//____________________________________________________________________________
AliEMCALAodCluster::AliEMCALAodCluster() : 
  AliAODCaloCluster(),fRecalibrated(0)
{
  // ctor
}
//____________________________________________________________________________
AliEMCALAodCluster::AliEMCALAodCluster(const AliAODCaloCluster & clu) : 
  AliAODCaloCluster(clu),fRecalibrated(0)
{
  // cpy ctor
}

//____________________________________________________________________________
AliEMCALAodCluster::~AliEMCALAodCluster()
{
  // dtor
}
//____________________________________________________________________________
void AliEMCALAodCluster::Recalibrate(AliEMCALCalibData * calibData, AliAODCaloCells *emcCells, TString emcalGeoName){
  //If not done yet, apply recalibration coefficients to energies list
  //NOTE that after recalibration fCellsAmpFraction contains not FRACTION but FULL energy 
  
  if(fRecalibrated)
   return ;
  
  if(!calibData)
    return ;

  AliEMCALGeometry * emcalgeo =  AliEMCALGeometry::GetInstance(emcalGeoName) ;
  if(!emcalgeo)
    AliFatal("AliEMCALGeometry was not constructed\n") ;
	
  Double32_t * cellsAmpFraction = GetCellsAmplitudeFraction(); 
  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  	
  for(Int_t i=0; i < GetNCells(); i++){
    
    //Get from the absid the supermodule, tower and eta/phi numbers
    emcalgeo->GetCellIndex(GetCellAbsId(i),iSupMod,iTower,iIphi,iIeta); 
    //Gives SuperModule and Tower numbers
    emcalgeo->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
					  iIphi, iIeta,iphi,ieta);   	
    
    Double_t energy = emcCells->GetCellAmplitude(GetCellAbsId(i)) ;
    AliDebug(2,Form("Recalibrate: cell %f, calib %f, fraction %f\n",energy,calibData->GetADCchannel(iSupMod,ieta,iphi),cellsAmpFraction[i]));
	cellsAmpFraction[i]*=energy*calibData->GetADCchannel(iSupMod,ieta,iphi);
  }
	
  SetCellsAmplitudeFraction(cellsAmpFraction);
  fRecalibrated=kTRUE; 
}
//____________________________________________________________________________
void  AliEMCALAodCluster::EvalAll(Float_t logWeight, TString geoname){
    //If recalibrated - recalculate all cluster parameters
  if(!fRecalibrated)
    return ;
  //printf("EvalAll e org %f\n",E());
  EvalEnergy() ; //Energy should be evaluated first
  //printf("EvalAll e2 %f\n",E());
  EvalPositionAndShowerShape(logWeight, geoname) ;
  //printf("EvalAll e3 %f\n",E());
  EvalPID() ; //Should be evaluated after energy and shower shape recalculation
  //printf("EvalAll e4 %f\n",E());
}
//____________________________________________________________________________
void AliEMCALAodCluster::EvalEnergy(){
  //Evaluate energy
  if(!fRecalibrated) // no need to recalibrate
    return ;
    
  Float_t energy=0. ;
  for(Int_t iDigit=0; iDigit < GetNCells(); iDigit++) {
    energy+=GetCellAmplitudeFraction(iDigit) ;
  }
  //printf("EvalEnergy: e %f\n", energy);
  SetE(energy);
  
   
}
////____________________________________________________________________________
//void AliEMCALAodCluster::EnergyCorrection(AliEMCALPID * pid){
//  //apply nonlinearity correction same as in AliEMCALPID.
//  SetE(pid->GetCalibratedEnergy(E())) ;
//}

//____________________________________________________________________________
void AliEMCALAodCluster::EvalPID(){           
	
  //re-evaluate identification parameters
//  pid->CalculatePID(E(),GetDispersion(),GetEmcCpvDistance(),GetTOF(),GetPID()) ;  
//  pid->CalculatePID(E(),GetDispersion(),GetM20(),GetM02(),GetEmcCpvDistance(),GetTOF(),GetPID()) ;

  //With bayesian
  AliEMCALPID *pid = new AliEMCALPID(kFALSE);
  pid->SetLowFluxParam(); // Need to be fixed
  Float_t pidlist[AliPID::kSPECIESN+1];
	for(Int_t i = 0; i < AliPID::kSPECIESN+1; i++) pidlist[i] = pid->GetPIDFinal(i);	
  SetPIDFromESD(pidlist);
	
}

//____________________________________________________________________________
void AliEMCALAodCluster::EvalPositionAndShowerShape(Float_t logWeight, TString emcalGeoName)
{
  // Calculates new center of gravity in the local EMCAL-module coordinates 
  // and tranfers into global ALICE coordinates
  // Calculates Dispersion and main axis
  if(!fRecalibrated) // no need to recalibrate
    return ;
  
  Int_t nstat  = 0;
  Float_t wtot = 0. ;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Double_t etai = -1.;
  Double_t phii = -1.;
	
  Double_t clXYZ[3] ={0.,0.,0.}; 
  Double_t xyzi[3]  ={0.,0.,0.};
 
  Double_t d     = 0.;
  Double_t dxx   = 0.;
  Double_t dzz   = 0.;
  Double_t dxz   = 0.;  
  Double_t xmean = 0.;
  Double_t zmean = 0.;

  AliEMCALGeometry * emcalgeo =  AliEMCALGeometry::GetInstance(emcalGeoName) ;
  if(!emcalgeo)
    AliFatal("AliEMCALGeometry was not constructed\n") ;
	
  Double_t dist  = TmaxInCm(Double_t(GetCellAmplitudeFraction(0)),0);
  for(Int_t iDigit=0; iDigit < GetNCells(); iDigit++) {
	
	//Get from the absid the supermodule, tower and eta/phi numbers
	emcalgeo->GetCellIndex(GetCellAbsId(iDigit),iSupMod,iTower,iIphi,iIeta); 
	emcalgeo->RelPosCellInSModule(GetCellAbsId(iDigit), dist, xyzi[0], xyzi[1], xyzi[2]);
	emcalgeo->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    Double_t ei = GetCellAmplitudeFraction(iDigit) ;
    if (E() > 0 && ei > 0) {
		Float_t w = ei;  
		if(logWeight > 0) w = TMath::Max( 0., logWeight + TMath::Log(ei/E()) ) ;
		etai=(Double_t)ieta;
		phii=(Double_t)iphi;		
		if(w > 0.0) {
			wtot += w ;
			nstat++;		
			for(Int_t i = 0; i < 3; i++ ) clXYZ[i]    += (w*xyzi[i]);
	
			//Shower shape
			dxx  += w * etai * etai ;
			xmean+= w * etai ;
			dzz  += w * phii * phii ;
			zmean+= w * phii ; 
			dxz  += w * etai * phii ; 
		}
    }
    else
      AliError(Form("Wrong energy %f and/or amplitude %f\n", ei, E()));
  }

  //Normalize to the weight	
  if (wtot > 0) {
		for(Int_t i=0; i<3; i++ ) clXYZ[i] /= wtot;
		xmean /= wtot ;
		zmean /= wtot ;
  }
  else
    AliError(Form("Wrong weight %f\n", wtot));

  //Put cluster position in the global system
  TVector3 gpos ;
  emcalgeo->GetGlobal(clXYZ, gpos, iSupMod);

  SetPosition(0, gpos[0]) ;
  SetPosition(1, gpos[1]) ;  
  SetPosition(2, gpos[2]) ;

  //Calculate dispersion	
  for(Int_t iDigit=0; iDigit < GetNCells(); iDigit++) {
		//Get from the absid the supermodule, tower and eta/phi numbers
		emcalgeo->GetCellIndex(GetCellAbsId(iDigit),iSupMod,iTower,iIphi,iIeta); 
		emcalgeo->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
		
		Double_t ei=GetCellAmplitudeFraction(iDigit) ;
		if (E() > 0 && ei > 0) {
			Float_t w = ei;  
			if(logWeight > 0) w = TMath::Max( 0., logWeight + TMath::Log(ei/E()) ) ;
			etai=(Double_t)ieta;
			phii=(Double_t)iphi;		
			if(w > 0.0)  d +=  w*((etai-xmean)*(etai-xmean)+(phii-zmean)*(phii-zmean)); 
		}
		else
			AliError(Form("Wrong energy %f and/or amplitude %f\n", ei, E()));
  }
	
  //Normalize to the weigth and set shower shape parameters
  if (wtot > 0 && nstat > 1) {
    d /= wtot ;
    dxx /= wtot ;
    dzz /= wtot ;
    dxz /= wtot ;
    dxx -= xmean * xmean ;
    dzz -= zmean * zmean ;
    dxz -= xmean * zmean ;
    SetM02(0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )) ;
    SetM20(0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz ));
  }
  else{
    d=0. ;
    SetM20(0.) ;
    SetM02(0.) ;
  }	
	
  if (d>=0)
	  SetDispersion(TMath::Sqrt(d)) ;
  else    
	  SetDispersion(0) ;

}

//_____________________________________________________________________
Double_t AliEMCALAodCluster::TmaxInCm(const Double_t e , const Int_t key) const
{ 
	// e energy in GeV)
	// key  =  0(gamma, default)
	//     !=  0(electron)
	static Double_t ca = 4.82;  // shower max parameter - first guess; ca=TMath::Log(1000./8.07)
	static Double_t x0 = 1.23;  // radiation lenght (cm)
	static Double_t tmax = 0.;   // position of electromagnetic shower max in cm
	
	tmax = 0.0;
	if(e>0.1) {
		tmax = TMath::Log(e) + ca;
		if      (key==0) tmax += 0.5; 
		else             tmax -= 0.5;
		tmax *= x0; // convert to cm
	}
	return tmax;
}

