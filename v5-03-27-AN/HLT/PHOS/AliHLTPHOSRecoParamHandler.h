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

#ifndef ALIHLTPHOSRECOPARAMHANDLER_H
#define ALIHLTPHOSRECOPARAMHANDLER_H

#include "AliHLTCaloRecoParamHandler.h"

class AliPHOSRecoParam;
class AliPHOSPIDv1;

class AliHLTPHOSRecoParamHandler : public AliHLTCaloRecoParamHandler
  {
     public:
     
	 /** Constructor */
	 AliHLTPHOSRecoParamHandler();
	
	 /** Destructor */
	 virtual ~AliHLTPHOSRecoParamHandler();
	 
	 /** Get the energy corrected for non-linear effects etc. */
	 virtual Float_t GetCorrectedEnergy(Float_t e);

     protected:
	
	/** See base class for documentation */
	virtual void FillParameters(); //COMMENT
	 
     private:
	
      /** Copy constructor, not implemented */
      AliHLTPHOSRecoParamHandler (const AliHLTPHOSRecoParamHandler &); //COMMENT
    
	/** Assignment operator, not implemented */
      AliHLTPHOSRecoParamHandler & operator = (const AliHLTPHOSRecoParamHandler &); //COMMENT
      
  };

#endif // ALIHLTPHOSRECOPARAMHANDLER_H
