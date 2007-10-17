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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.2  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Base class for reading data in order to do prompt gamma correlations
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

//---- ANALYSIS system ----
#include "Riostream.h"
#include "AliLog.h"
#include "AliGammaReader.h"

ClassImp(AliGammaReader)


//____________________________________________________________________________
AliGammaReader::AliGammaReader() : 
  TObject(), fDataType(0),
  fCTSEtaCut(0.), fEMCALEtaCut(0.), fPHOSEtaCut(0.),
  fNeutralPtCut(0.), fChargedPtCut(0.),
  fEMCALIPDistance(0.),   fPHOSIPDistance(0.), 
  fEMCALMinAngle(0.),   fPHOSMinAngle(0.), 
  fEMCALPID(0),fPHOSPID(0),
  fEMCALPhotonWeight(0.), fEMCALPi0Weight(0.),  fEMCALElectronWeight(0.),  
  fEMCALChargeWeight(0.),fEMCALNeutralWeight(0.),
  fPHOSPhotonWeight(0.), fPHOSPi0Weight(0.),  fPHOSElectronWeight(0.), 
  fPHOSChargeWeight(0.) , fPHOSNeutralWeight(0.) ,
  fPHOSWeightFormula(0), fPHOSPhotonWeightFormula(0x0), fPHOSPi0WeightFormula(0x0) 
{
  //Ctor

  fPhiEMCALCut[0]=0.;
  fPhiEMCALCut[1]=0.;
  fPhiPHOSCut[0]=0.;
  fPhiPHOSCut[1]=0.;
 
  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliGammaReader::AliGammaReader(const AliGammaReader & g) :   
  TObject(g), fDataType(g.fDataType),
  fCTSEtaCut(g.fCTSEtaCut),  fEMCALEtaCut(g.fEMCALEtaCut),  fPHOSEtaCut(g.fPHOSEtaCut),
  fNeutralPtCut(g.fNeutralPtCut), fChargedPtCut(g.fChargedPtCut),
  fEMCALIPDistance(g.fEMCALIPDistance),  fPHOSIPDistance(g.fPHOSIPDistance),  
  fEMCALMinAngle(g.fEMCALMinAngle),  fPHOSMinAngle(g.fPHOSMinAngle), 
  fEMCALPID(g.fEMCALPID), 
  fPHOSPID(g.fPHOSPID),
  fEMCALPhotonWeight(g.fEMCALPhotonWeight), 
  fEMCALPi0Weight(g.fEMCALPi0Weight), 
  fEMCALElectronWeight(g.fEMCALElectronWeight), 
  fEMCALChargeWeight(g.fEMCALChargeWeight), 
  fEMCALNeutralWeight(g.fEMCALNeutralWeight), 
  fPHOSPhotonWeight(g.fPHOSPhotonWeight),
  fPHOSPi0Weight(g.fPHOSPi0Weight),
  fPHOSElectronWeight(g.fPHOSElectronWeight), 
  fPHOSChargeWeight(g.fPHOSChargeWeight),
  fPHOSNeutralWeight(g.fPHOSNeutralWeight),
  fPHOSWeightFormula(g.fPHOSWeightFormula), 
  fPHOSPhotonWeightFormula(g.fPHOSPhotonWeightFormula), 
  fPHOSPi0WeightFormula(g.fPHOSPi0WeightFormula) 
{
  // cpy ctor

  fPhiEMCALCut[0]=g.fPhiEMCALCut[0];
  fPhiEMCALCut[1]=g.fPhiEMCALCut[1];
  fPhiPHOSCut[0]=g.fPhiPHOSCut[0];
  fPhiPHOSCut[1]=g.fPhiPHOSCut[1];
}

//_________________________________________________________________________
AliGammaReader & AliGammaReader::operator = (const AliGammaReader & source)
{
  // assignment operator

  if(&source == this) return *this;

  fDataType = source.fDataType ;
  fCTSEtaCut = source.fCTSEtaCut;  
  fEMCALEtaCut = source.fEMCALEtaCut;  
  fPHOSEtaCut = source.fPHOSEtaCut;
  fNeutralPtCut = source.fNeutralPtCut;
  fChargedPtCut = source.fChargedPtCut; 

  fPhiEMCALCut[0]=source.fPhiEMCALCut[0];
  fPhiEMCALCut[1]=source.fPhiEMCALCut[1];
  fPhiPHOSCut[0]=source.fPhiPHOSCut[0];
  fPhiPHOSCut[1]=source.fPhiPHOSCut[1];

  fEMCALIPDistance = source.fEMCALIPDistance; 
  fPHOSIPDistance = source.fPHOSIPDistance; 
  fEMCALMinAngle = source.fEMCALMinAngle; 
  fPHOSMinAngle = source.fPHOSMinAngle; 

  fEMCALPID = source.fEMCALPID ;
  fPHOSPID = source.fPHOSPID ;

  fEMCALPhotonWeight = source. fEMCALPhotonWeight ;
  fEMCALPi0Weight = source.fEMCALPi0Weight ;
  fEMCALElectronWeight = source.fEMCALElectronWeight; 
  fEMCALChargeWeight = source.fEMCALChargeWeight;
  fEMCALNeutralWeight = source.fEMCALNeutralWeight;

  fPHOSPhotonWeight = source.fPHOSPhotonWeight ;
  fPHOSPi0Weight = source.fPHOSPi0Weight ;
  fPHOSElectronWeight = source.fPHOSElectronWeight; 
  fPHOSChargeWeight = source.fPHOSChargeWeight;
  fPHOSNeutralWeight = source.fPHOSNeutralWeight;

  fPHOSWeightFormula       = source.fPHOSWeightFormula; 
  fPHOSPhotonWeightFormula = source.fPHOSPhotonWeightFormula; 
  fPHOSPi0WeightFormula    = source.fPHOSPi0WeightFormula;


  return *this;

}

//_______________________________________________________________
void AliGammaReader::InitParameters()
{
 
  //Initialize the parameters of the analysis.
  fDataType = kData ;
  fCTSEtaCut         = 0.7 ;  
  fEMCALEtaCut         = 0.7 ;  
  fPHOSEtaCut         = 0.12 ;
  fPhiEMCALCut[0] = 80 *TMath::DegToRad();
  fPhiEMCALCut[1] = 190*TMath::DegToRad();
  fPhiPHOSCut[0] = 220. *TMath::DegToRad();
  fPhiPHOSCut[1] = 320.*TMath::DegToRad();
  fNeutralPtCut   = 0.5 ;
  fChargedPtCut   = 0.5 ;

  fEMCALMinAngle    =   2.5 * TMath::DegToRad() ;  
  fPHOSMinAngle    = 0.45 * TMath::DegToRad() ; //3.6 ;
  fEMCALIPDistance    = 450. ;//cm  
  fPHOSIPDistance    = 445. ;//cm 460 (EMCA) - 15 (CPV)

  //pid, only for ESD data
  fEMCALPID = kFALSE;
  fPHOSPID = kFALSE;

  fEMCALPhotonWeight = 0.8 ;
  fEMCALPi0Weight = 0.5 ;
  fEMCALElectronWeight = 0.8 ;
  fEMCALChargeWeight = 0.5 ;
  fEMCALNeutralWeight = 0.5 ;

  fPHOSPhotonWeight = 0.75 ;
  fPHOSPi0Weight = 0.8 ;
  fPHOSElectronWeight = 0.5 ;
  fPHOSChargeWeight = 0.5 ;
  fPHOSNeutralWeight = 0.5 ;

  //Formula to set the PID weight threshold for photon or pi0
  fPHOSWeightFormula = kTRUE;
  fPHOSPhotonWeightFormula = 
    new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
  fPHOSPi0WeightFormula = 
    new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");

}


//________________________________________________________________
void AliGammaReader::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;

  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  printf("Data type           : %d\n", fDataType) ;
  printf("CTS Eta cut           : %f\n", fCTSEtaCut) ;
  printf("EMCAL Eta cut           : %f\n", fEMCALEtaCut) ;
  printf("PHOS Eta cut           : %f\n", fPHOSEtaCut) ;
  printf("Phi EMCAL cut           : [%f, %f]\n", fPhiEMCALCut[0],fPhiEMCALCut[1]) ;
  printf("Phi PHOS cut           : [%f, %f]\n", fPhiPHOSCut[0],fPhiPHOSCut[1]) ;
  printf("pT neutral cut           : %f GeV/c\n", fNeutralPtCut) ;
  printf("pT charged cut           : %f GeV/c\n", fChargedPtCut) ;

  if(fDataType == kMC || fDataType == kMCData){
    printf("IP distance to PHOS         : %f\n", fPHOSIPDistance) ;
    printf("IP distance to EMCAL         : %f\n", fEMCALIPDistance) ;
    printf("Min gamma decay aperture angle in PHOS         : %f\n", fPHOSMinAngle) ;
    printf("Min gamma decay aperture angle in EMCAL         : %f\n", fEMCALMinAngle) ;
  }
  
  if(fDataType != kMC){
    printf("PHOS PID on?               =     %d\n",  fPHOSPID) ; 
    printf("EMCAL PID  on?         =     %d\n",  fEMCALPID) ;
    printf("PHOS PID weight , photon %f, pi0 %f, e %f, charge %f, neutral %f \n",  
	   fPHOSPhotonWeight,  fPHOSPi0Weight, 
	   fPHOSElectronWeight,  fPHOSChargeWeight,   fPHOSNeutralWeight) ; 
    printf("EMCAL PID weight, photon %f, pi0 %f, e %f, charge %f, neutral %f\n",   
	   fEMCALPhotonWeight,  fEMCALPi0Weight, 
	   fEMCALElectronWeight,  fEMCALChargeWeight,  fEMCALNeutralWeight) ; 
    
    printf("PHOS Parametrized weight on?               =     %d\n",  fPHOSWeightFormula) ; 
    if(fPHOSWeightFormula){
      printf(">>>>>>>>>>> Photon weight formula<<<<<<<<<<<<\n");
      fPHOSPhotonWeightFormula->Print();
      printf(">>>>>>>>>>> Pi0    weight formula<<<<<<<<<<<<\n");
      fPHOSPhotonWeightFormula->Print();
    }
  }

} 


