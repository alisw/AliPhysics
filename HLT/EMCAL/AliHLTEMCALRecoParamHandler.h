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

#ifndef ALIHLTEMCALRECOPARAMHANDLER_H
#define ALIHLTEMCALRECOPARAMHANDLER_H

#include "AliHLTCaloRecoParamHandler.h"

class AliHLTEMCALRecoParamHandler : public AliHLTCaloRecoParamHandler
  {
     public:
     
	 /** Constructor */
	 AliHLTEMCALRecoParamHandler();
	
	 /** Destructor */
	 virtual ~AliHLTEMCALRecoParamHandler();
	 
	 /** See base class for documentation */
	 virtual Int_t GetParametersFromCDB();
	 
	 /** Get the energy corrected for non-linear effects etc. */
	 virtual Float_t GetCorrectedEnergy(Float_t e);

     protected:
	
	/** Fill the parameters */
	virtual void FillParameters();
	 
     private:
	
      /** Copy constructor, not implemented */
      AliHLTEMCALRecoParamHandler (const AliHLTEMCALRecoParamHandler &); //COMMENT
    
	/** Assignment operator, not implemented */
      AliHLTEMCALRecoParamHandler & operator = (const AliHLTEMCALRecoParamHandler &); //COMMENT

  };

#endif // ALIHLTPHOSRECOPARAMHANDLER_H
