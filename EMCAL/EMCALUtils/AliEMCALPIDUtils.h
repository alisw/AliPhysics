#ifndef ALIEMCALPIDUTILS_H
#define ALIEMCALPIDUTILS_H

//_________________________________________________________________________
/// \class AliEMCALPIDUtils
/// \brief Compute cluster PID weights 
///
///   Compute PID weights for all the clusters that are in AliESDs.root file
///   the AliESDs.root have to be in the same directory as the class.
///
///   It is called during reconstruction by EMCALrec/AliEMCALPID. 
///   It can be executed at analysis level.
///
///   Implementation based on simulations before 2009, executed in reconstruction but
///   not really used after 2009. To be revisited.
///
///   Do:    
///   AliEMCALPIDUtils *pid = new AliEMCALPIDUtils();
///   pid->SetPrintInfo(kTRUE);
///   pid->SetHighFluxParam(); //   pid->SetLowFluxParam(); 
///   
///   then in cluster loop do
///   pid->ComputePID(energy, lambda0);
///  	  
///   Compute PID Weight for all clusters in AliESDs.root file
///   keep this function for the moment for a simple verification, could be removed
///
///   Get back the probabilities with 
///   pid->GetPIDFinal(pidFinal) 
///
///   where Double_t pidFinal[AliPID::kSPECIESCN] is the standard PID for :
///
///	  kElectron :  fPIDFinal[0]
///	  kMuon     :  fPIDFinal[1]
///	  kPion	    :  fPIDFinal[2]
///	  kKaon	    :  fPIDFinal[3]
///	  kProton   :  fPIDFinal[4]
///	  kPhoton   :  fPIDFinal[5]
///	  kPi0	    :  fPIDFinal[6]
///	  kNeutron  :  fPIDFinal[7]
///	  kKaon0    :  fPIDFinal[8]
///	  kEleCon   :  fPIDFinal[9]
///	  kUnknown  :  fPIDFinal[10]
///
/// \author Genole Bourdaud, SUBATECH
/// First implementation, 2007
/// \author Marie Germain, <Marie.Germain@subatech.in2p3.fr>, SUBATECH 
/// New parametrization for low and high flux environment, 07/2009
/// \author Gustavo Conesa, <Gustavo.Conesa.Balbastre@cern.ch>, LNF, 08/2009 
/// Divide class in AliEMCALPID and AliEMCALPIDUtils, PIDUtils belong to library EMCALUtils 
///

// Root includes
#include "TNamed.h"
class TArrayD ;

// AliRoot includes
#include "AliPID.h" 

class AliEMCALPIDUtils : public TNamed {

public:
  
  AliEMCALPIDUtils();
    
  void     ComputePID(Double_t energy, Double_t lambda0); 

  void     InitParameters();
  void     SetLowFluxParam();
  void     SetHighFluxParam();

  TArrayD  DistLambda0(Double_t energy, Int_t nature) ; 
  
  Double_t DistEnergy (Double_t energy, Int_t nature) ;

  Double_t GetPID      (Int_t idx) const { if (idx>=0&&idx<3) return fPID[idx];       else return 0. ; }
  Double_t GetPIDFinal (Int_t idx) const { if (idx>=0&&idx<AliPID::kSPECIESCN) return fPIDFinal[idx]; else return 0. ; }
  Double_t GetPIDWeight(Int_t idx) const { if (idx>=0&&idx<3) return fPIDWeight[idx]; else return 0. ; }
  
  void     SetPID      (Double_t val, Int_t idx) { if (idx>=0&&idx<3) fPID[idx]       = val ; }
  void     SetPIDFinal (Double_t val, Int_t idx) { if (idx>=0&&idx<AliPID::kSPECIESCN) fPIDFinal[idx] = val ; }
  void     SetPIDWeight(Double_t val, Int_t idx) { if (idx>=0&&idx<3) fPIDWeight[idx] = val ; }

  void     SetPrintInfo(Bool_t yesno)            { fPrintInfo = yesno ; }
	
 private:
  
  Double_t Polynomial      (Double_t x, const Double_t * params) const ;
  Double_t Polynomialinv   (Double_t x, const Double_t * params) const ;
  Double_t PolynomialMixed1(Double_t x, const Double_t * params) const ;
  Double_t PolynomialMixed2(Double_t x, const Double_t * params) const ;
  Double_t Polynomial0     (            const Double_t * params) const ;
  Double_t PowerExp        (Double_t x, const Double_t * params) const ;
	
protected: 	
  Bool_t   fPrintInfo;            ///< Flag to decide if details about PID must be printed.
  
  Double_t fGamma[6][6];          ///< Parameter to Compute PID for photons.
  Double_t fGamma1to10[6][6];     ///< Parameter to Compute PID not used.
  Double_t fHadron[6][6];         ///< Parameter to Compute PID for hadrons, 1 to 10 GeV.
  Double_t fHadron1to10[6][6];	  ///< Parameter to Compute PID for hadrons, 1 to 10 GeV.
  Double_t fPiZero[6][6];         ///< Parameter to Compute PID for pi0.
  
  Double_t fHadronEnergyProb[6];  ///< Parameter to Compute PID for energy ponderation for hadrons. 	 
  Double_t fPiZeroEnergyProb[6];  ///< Parameter to Compute PID for energy ponderation for Pi0.  	 
  Double_t fGammaEnergyProb[6];   ///< Parameter to Compute PID for energy ponderation for gamma. 	 
  
  Float_t  fPID[3];               ///< Order: gamma, pi0, hadrons. Final container?
  Float_t  fPIDFinal[AliPID::kSPECIESCN+1]; ///< Final PID format, all species.
  Float_t  fPIDWeight[3];         ///< Order: gamma, pi0, hadrons. Temporal container?
  
  Double_t fProbGamma;	          ///< Probability to be a Gamma.
  Double_t fProbPiZero;	          ///< Probability to be a PiO.
  Double_t fProbHadron;	          ///< Probability to be a Hadron.
  
  Double_t fWeightHadronEnergy;	  ///< Weight for a Hadron to have a given energy  (parametr from a flat distrib from 0 to 100).
  Double_t fWeightGammaEnergy;	  ///< Weight for a Gamma to have a given energy  (for the moment =1.).
  Double_t fWeightPiZeroEnergy;	  ///< Weight for a Pi0 Hadron to have a given energy (for the moment =1.).
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALPIDUtils, 2) ;
  /// \endcond

};

#endif // ALIEMCALPIDUTILS_H


