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

#ifndef ALIHLTCALOHISTOINVMASS_H
#define ALIHLTCALOHISTOINVMASS_H

/** 
 * @file   AliHLTCaloHistoInvMass
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
//#include "PHOS/AliHLTPHOSConstants.h" 

class TObjArray;
class TH1F;
class AliHLTCaloClusterReader;
class AliESDEvent;
class TRefArray;
struct AliHLTCaloClusterHeaderStruct;
class TString;
class TRefArray;

/** 
 * @class AliHLTCaloHistoInvMass
 *
 * Class produces physics histograms for PHOS. It takes a TClonesArray
 * of AliESDCalocluster as input and fills several histograms
 *
 * Histograms (1D):
 *  * - Invariant mass of two clusters
 * 
 * @ingroup alihlt_phos
 */

class AliHLTCaloHistoInvMass {

public:
  
  /** Constructor */
  AliHLTCaloHistoInvMass(TString det);

  /** Destructor */
  virtual ~AliHLTCaloHistoInvMass();

  /** Get a pointer to the TObjArray of histograms */
  TObjArray *GetHistograms();

  //** Loops of the calo clusters and fills histos
  Int_t FillHistograms(Int_t nc, TRefArray * clustersArray);
  
private:

  /** Default constructor prohibited */
  AliHLTCaloHistoInvMass();
  
  /** Copy constructor */
  AliHLTCaloHistoInvMass(const AliHLTCaloHistoInvMass &);
  
  /** Assignment operator */
  AliHLTCaloHistoInvMass & operator= (const AliHLTCaloHistoInvMass);
  
  /** Cluster reader class   */
  AliHLTCaloClusterReader * fClusterReader;
  
  /** Histogram of the 2 cluster invariant mass */
  TH1F *fHistTwoClusterInvMass;                 //!transient
  
  /** Pointer to the array of histograms */
  TObjArray *fHistArrayPtr;                     //!transient
  
  ClassDef(AliHLTCaloHistoInvMass, 0);

};
 
#endif
