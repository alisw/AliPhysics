#ifndef ALIEMCALPID_H
#define ALIEMCALPID_H

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
// Class AliEMCALPID
// Compute PID weights for all the clusters
///////////////////////////////////////////////////////////////////////////////

//Root includes
class TArrayD ;

//AliRoot includes
class AliESDEvent ;
#include "AliEMCALPIDUtils.h" 

class AliEMCALPID : public AliEMCALPIDUtils {

public:
  
  AliEMCALPID();
  AliEMCALPID(Bool_t reconstructor);
  //virtual ~AliEMCALPID() { }
  
  void    RunPID(AliESDEvent *esd);
  void    InitParameters();
  void    SetReconstructor(Bool_t yesno) {fReconstructor = yesno;}
	
 private:
  
  Bool_t   fReconstructor;                // Fill esdcalocluster when called from EMCALReconstructor
  
  ClassDef(AliEMCALPID, 5)

};


#endif // ALIEMCALPID_H


