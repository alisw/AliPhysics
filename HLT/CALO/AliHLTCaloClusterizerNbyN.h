
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


#ifndef ALIHLTCALOCLUSTERIZERNBYN_H
#define ALIHLTCALOCLUSTERIZERNBYN_H

#include "AliHLTCaloClusterizer.h"


/**
 * Class does clusterization in for Calorimeters on an event basis. It is intended 
 * for use in HLT, but can also be used offline. It uses the static N by N method
 *
 * @file   AliHLTCaloClusterizerNbyN.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Clusterizer for CALO HLT
 */

class AliHLTCaloClusterizerNbyN : public AliHLTCaloClusterizer
{

public:
   
   /** Constructor */ 
   AliHLTCaloClusterizerNbyN(TString det);
    
   /** Destrucotr */
    virtual ~AliHLTCaloClusterizerNbyN();
   
   /** Clusterize the event */
    virtual Int_t ClusterizeEvent(Int_t nDigits);

    /** set the dimension of the grid */
    void SetGridDimension(Int_t n) { fN = n; }
    
   private:
      
      /** The "N" in N by N clustering */
      Int_t fN; // See above
     /** Default constructor, prohibited */
  AliHLTCaloClusterizerNbyN();                          // COMMENT
  
  /** Copy constructor, prohibited */
  AliHLTCaloClusterizerNbyN (const AliHLTCaloClusterizerNbyN &); //COMMENT
  
  /** Assignment operator, prohibited */
  AliHLTCaloClusterizerNbyN & operator = (const AliHLTCaloClusterizerNbyN &); //COMMENT


};

#endif // ALIHLTCALOCLUSTERIZERNBYN_H
