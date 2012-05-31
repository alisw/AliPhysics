
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Oystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#ifndef ALIHLTCALORECOPARAMHANDLER_H
#define ALIHLTCALORECOPARAMHANDLER_H

#include "AliHLTCaloConstantsHandler.h"
#include "AliHLTLogging.h"
#include "AliCDBPath.h"

class AliDetectorRecoParam;
class AliHLTCaloRecoParamHandler: public AliHLTCaloConstantsHandler, public AliHLTLogging
  {
  public:
    /** Constructor, */
      AliHLTCaloRecoParamHandler(TString det); // See above
      /** Destructor */
      virtual ~AliHLTCaloRecoParamHandler(); // See above
      /** Get the parameters from CDB */
      virtual Int_t GetParametersFromCDB(); // See above
      /** Get the corrected energy, should really be implemented by the child classes */
      virtual Float_t GetCorrectedEnergy(Float_t e) { return e; }
      /** Get the log weight */
      Float_t GetLogWeight() { return fLogWeight; } 
      /** Get rec point threshold */
      Float_t GetRecPointThreshold() { return fRecPointThreshold; }
      /** Get rec point member threshold */
      Float_t GetRecPointMemberThreshold() { return fRecPointMemberThreshold; }

     protected:
      /** Fill the parameters */
      virtual void FillParameters() = 0; //COMMENT
      /** The log weight used in calculating the cluster position */
      Float_t fLogWeight; //COMMENT
      /** The threshold for adding a digit to a recpoint */
      Float_t fRecPointMemberThreshold; //COMMENT
      /** The threshold for starting a recpoint */
      Float_t fRecPointThreshold; //COMMENT
      /** A reco param object */
      AliDetectorRecoParam *fRecoParamPtr; 	//! transient
      /** CDB path to the reco param object */
      AliCDBPath fRecoParamPath;    //COMMENT

     private:
	
      /** Default constructor, inhibited */
      AliHLTCaloRecoParamHandler(); // See above
      
      /** Copy constructor, not implemented */
      AliHLTCaloRecoParamHandler (const AliHLTCaloRecoParamHandler &);  //COMMENT
    
      /** Assignment operator, not implemented */
      AliHLTCaloRecoParamHandler & operator = (const AliHLTCaloRecoParamHandler &);  //COMMENT

      ClassDef(AliHLTCaloRecoParamHandler, 0);
	
  };

#endif // ALIHLTCALORECOPARAMHANDLER_H
