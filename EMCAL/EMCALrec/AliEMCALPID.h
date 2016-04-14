#ifndef ALIEMCALPID_H
#define ALIEMCALPID_H

//_________________________________________________________________________
/// \class AliEMCALPID
/// \brief Compute cluster PID weights 
///
///   Compute PID weights for all the clusters that are in AliESDs.root file
///   the AliESDs.root have to be in the same directory as the class.
///   Wrapper class to be used during reconstruction, the main algorithm sits in 
///   EMCALUtils/AliEMCALPIDUtils.
///
///   Implementation based on simulations before 2009, executed in reconstruction but
///   not really used after 2009. To be revisited.
///
///   Do:    
///   AliEMCALPID *pid = new AliEMCALPID(kFALSE); // this calls the constructor which avoids the call to recparam 
///   pid->SetReconstructor(kFALSE);
///   pid->SetPrintInfo(kTRUE);
///   pid->SetHighFluxParam(); //   pid->SetLowFluxParam(); 
///   
///   then in cluster loop do
///   pid->ComputePID(energy, lambda0);
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
class TArrayD ;

// AliRoot includes
class AliESDEvent ;
#include "AliEMCALPIDUtils.h" 

class AliEMCALPID : public AliEMCALPIDUtils {

public:
  
  AliEMCALPID();
  
  AliEMCALPID(Bool_t reconstructor);
  
  //virtual ~AliEMCALPID() { }
  
  void    RunPID(AliESDEvent *esd);
  
  void    InitParameters();
  
  void    SetReconstructor(Bool_t yesno) { fReconstructor = yesno ; }
	
 private:
  
  Bool_t   fReconstructor ;  ///< Fill AliESDCaloCluster object when called from AliEMCALReconstructor
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALPID, 5) ;
  /// \endcond

};

#endif // ALIEMCALPID_H


