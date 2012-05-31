//-*- Mode: C++ -*-
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Albin Gaignette                                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTCALOHISTOPRODUCER_H
#define ALIHLTCALOHISTOPRODUCER_H

/** 
 * @file   AliHLTCaloHistoProducer
 * @author Svein Lindal slindal@fys.uio.no
 * @date 
 * @brief  Base class for calo physics histogram producers
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "Rtypes.h"
#include <vector>
#include "TObject.h"
#include "AliHLTLogging.h"

class TRefArray;
class AliHLTCaloClusterReader;
class AliHLTCaloClusterDataStruct;
class TObjArray;

/** 
 * @class AliHLTCaloHistoProducer
 *
 * Base class for calo physics histogram producers
 *
 *
 * @ingroup alihlt_phos
 */

using std::vector;

class AliHLTCaloHistoProducer : public TObject, public AliHLTLogging {

public:
  
  /** Constructor */
  AliHLTCaloHistoProducer();

  /** Destructor */
  virtual ~AliHLTCaloHistoProducer();

  /** Get a pointer to the TObjArray of histograms */
  TObjArray *GetHistograms();

  //** Loops of the calo clusters and fills histos
  virtual Int_t FillHistograms(Int_t nc, TRefArray * clusterArray ) = 0;
  virtual Int_t FillHistograms(Int_t nc, vector<AliHLTCaloClusterDataStruct*> &cVec) = 0;

protected:
  /** Cluster reader class   */
  AliHLTCaloClusterReader * fClusterReader;
  
  /** Pointer to the array of histograms */
  TObjArray *fHistArray;                     //!transient
  

private:

  /** Copy constructor */
  AliHLTCaloHistoProducer(const AliHLTCaloHistoProducer &);
  
  /** Assignment operator */
  AliHLTCaloHistoProducer & operator = (const AliHLTCaloHistoProducer &);
  
  ClassDef(AliHLTCaloHistoProducer, 0);

};
#endif
