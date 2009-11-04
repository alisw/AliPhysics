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

#ifndef ALIHLTPHOSINVMASSHISTOGRAMPRODUCER_H
#define ALIHLTPHOSINVMASSHISTOGRAMPRODUCER_H

/** 
 * @file   AliHLTPHOSInvMassHistogramProducer
 * @author Albin Gaignette and Svein Lindal slindal@fys.uio.no
 * @date 
 * @brief  Produces Invariant mass histograms of PHOS clusters
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

//#include "AliHLTPHOSBase.h"

#include "Rtypes.h"
// #include "TClonesArray.h"

#include "AliHLTPHOSConstants.h" 

using namespace  PhosHLTConst;

class TObjArray;
class TH1F;
//class TH2F;
class AliHLTCaloClusterReader;
struct AliHLTCaloClusterHeaderStruct;


/** 
 * @class AliHLTPHOSInvMassHistogramProducer
 *
 * Class produces physics histograms for PHOS. It takes a TClonesArray
 * of AliESDCalocluster as input and fills several histograms
 *
 * Histograms (1D):
 * - Total number of clusters per event
 * - Energy distribution of clusters
 * - Total energy in event
 * - Invariant mass of two clusters
 * - Number of cells in clusters
 * - Fraction of cells with energy deposit
 * 
 * Histograms (2D):
 * - Number of cells in cluster vs cluster energy
 * - Number of clusters vs total energy
 *
 * @ingroup alihlt_phos
 */



//class AliHLTPHOSInvMassHistogramProducer : public AliHLTPHOSBase
class AliHLTPHOSInvMassHistogramProducer 
{
 public:

  /** Constructor */
  AliHLTPHOSInvMassHistogramProducer();

  /** Destructor */
  virtual ~AliHLTPHOSInvMassHistogramProducer();

  /** Copy constructor */
  AliHLTPHOSInvMassHistogramProducer(const AliHLTPHOSInvMassHistogramProducer &) :
    fClusterReader(NULL),
    fHistTwoClusterInvMass(0),
    fHistArrayPtr(0)
  {
    // Copy constructor not implemented
  }

  /** Assignment operator */
  AliHLTPHOSInvMassHistogramProducer & operator= (const AliHLTPHOSInvMassHistogramProducer)
  {
    // assignment
    return *this;
  }

  /** Analyse the clusters in the event */
  int DoEvent(AliHLTCaloClusterHeaderStruct* cHeader);

  /** Get a pointer to the TObjArray of histograms */
  TObjArray *GetHistograms();

  
  
 private:

  /** Cluster reader class   */
  AliHLTCaloClusterReader * fClusterReader;

  /** Histogram of the 2 cluster invariant mass */
  TH1F *fHistTwoClusterInvMass;                 //!transient

  /** Pointer to the array of histograms */
  TObjArray *fHistArrayPtr;                     //!transient

  ClassDef(AliHLTPHOSInvMassHistogramProducer, 0);

};
 
#endif
