#ifndef ALIEMCALPIDUTILS_H
#define ALIEMCALPIDUTILS_H

/* $Id: AliEMCALPIDUtils.h 33808 2009-07-15 09:48:08Z gconesab $ */

///////////////////////////////////////////////////////////////////////////////
// Class AliEMCALPIDUtils
// Compute PID weights for all the clusters
///////////////////////////////////////////////////////////////////////////////

//Root includes
#include "TNamed.h"
class TArrayD ;

//AliRoot includes
#include "AliPID.h" 

class AliEMCALPIDUtils : public TNamed {

public:
  
  AliEMCALPIDUtils();
/*   AliEMCALPIDUtils(Bool_t reconstructor); */
  //virtual ~AliEMCALPIDUtils() { }
  
  void     ComputePID(Double_t energy, Double_t lambda0); // give the PID of a cluster

  void     InitParameters();
  void     SetLowFluxParam();
  void     SetHighFluxParam();

  TArrayD  DistLambda0(const Double_t energy, const Int_t nature) ; // compute lambda0 distributions
  
  Double_t DistEnergy(const Double_t energy, const Int_t nature) ;

  Double_t GetPID(Int_t idx) const {if (idx>=0&&idx<3) return fPID[idx]; else return 0.;}
  Double_t GetPIDFinal(Int_t idx) const {if (idx>=0&&idx<AliPID::kSPECIESN) return fPIDFinal[idx]; else return 0.;}
  Double_t GetPIDWeight(Int_t idx) const {if (idx>=0&&idx<3) return fPIDWeight[idx]; else return 0.;}
  
  void    SetPID(Double_t val, Int_t idx) {if (idx>=0&&idx<3) fPID[idx] = val;}
  void    SetPIDFinal(Double_t val, Int_t idx) {if (idx>=0&&idx<AliPID::kSPECIESN) fPIDFinal[idx] = val;}
  void    SetPIDWeight(Double_t val, Int_t idx) {if (idx>=0&&idx<3) fPIDWeight[idx] = val;}
  void    SetPrintInfo(Bool_t yesno) {fPrintInfo = yesno;}
	
   	
	
 private:
  
  Double_t Polynomial(const Double_t x, const Double_t *params) const ;
  Double_t Polynomialinv(const Double_t x, const Double_t *params) const ;
  Double_t PolynomialMixed1(const Double_t x, const Double_t *params) const ;
  Double_t PolynomialMixed2(const Double_t x, const Double_t *params) const ;
  Double_t Polynomial0(const Double_t *params) const ;
  Double_t PowerExp(const Double_t x, const Double_t *params) const ;
	
protected: 	
  Bool_t   fPrintInfo;          // flag to decide if details about PID must be printed
  
  Double_t fGamma[6][6];            // Parameter to Compute PID for photons
  Double_t fGamma1to10[6][6];       // Parameter to Compute PID not used
  Double_t fHadron[6][6];	        // Parameter to Compute PID for hadrons, 1 to 10 GeV
  Double_t fHadron1to10[6][6];	    // Parameter to Compute PID for hadrons, 1 to 10 GeV
  Double_t fPiZero[6][6];           // Parameter to Compute PID for pi0
  Double_t fHadronEnergyProb[6]; 	// Parameter to Compute PID for energy ponderation for hadrons  	 
  Double_t fPiZeroEnergyProb[6]; 	// Parameter to Compute PID for energy ponderation for Pi0  	 
  Double_t fGammaEnergyProb[6]; 	// Parameter to Compute PID for energy ponderation for gamma  	 
   
  Float_t fPID[3];
  
  Float_t  fPIDFinal[AliPID::kSPECIESN+1]; // final PID format
  Float_t  fPIDWeight[3];                  // order: gamma, pi0, hadrons,
  Double_t fProbGamma;	                  // probility to be a Gamma
  Double_t fProbPiZero;	                  // probility to be a PiO
  Double_t fProbHadron;	                  // probility to be a Hadron
  Double_t fWeightHadronEnergy;	          // Weight for a  a Hadron to have a given energy  (parametr from a flat distrib from 0 to 100)
  Double_t fWeightGammaEnergy;	          // Weight for a  Gamma to have a given energy  (for the moment =1.)
  Double_t fWeightPiZeroEnergy;	          // Weight for a Pi0 Hadron to have a given energy (for the moment =1.)
  
  ClassDef(AliEMCALPIDUtils, 2)

};

#endif // ALIEMCALPIDUTILS_H


