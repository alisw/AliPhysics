//-*- Mode: C++ -*-
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Svein Lindal                                          *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** 
 * @file   AliHLTCaloHistoProducer
 * @author Svein Lindal <slindal@fys.uio.no>
 * @date 
 * @brief  Base class for Calo Physics histogram producers
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCaloHistoProducer.h"
#include "TObjArray.h"
#include "AliHLTCaloClusterReader.h"

AliHLTCaloHistoProducer::AliHLTCaloHistoProducer() :
  TObject(), 
  AliHLTLogging(),
  fClusterReader(NULL),
  fHistArray(NULL)
{
  // See header file for documentation
  fHistArray = new TObjArray;
  fClusterReader = new AliHLTCaloClusterReader();
}

AliHLTCaloHistoProducer::~AliHLTCaloHistoProducer()
{
  //destructor
  if(fClusterReader)
    delete fClusterReader;
  fClusterReader = NULL;

  if(fHistArray)
    delete fHistArray;
  fHistArray = NULL;

}

TObjArray* AliHLTCaloHistoProducer::GetHistograms()
{  
  // See header file for documentation
  return fHistArray;
}
