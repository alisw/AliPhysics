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

#ifndef ALIHLTPHOSHISTOPRODCELLENERGY
#define ALIHLTPHOSHISTOPRODCELLENERGY

/** 
 * @file   AliHLTPHOSHistoProdCellEnergy
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
#include "AliHLTPHOSConstants.h" 

class TObjArray;
class TH1F;
class TH2F;
class AliHLTCaloClusterReader;
struct AliHLTCaloClusterHeaderStruct;


using namespace  PhosHLTConst;


/** 
 * @class AliHLTPHOSHistoProdCellEnergy
 *
 * Class produces physics histograms for PHOS. It takes a TClonesArray
 * of AliESDCalocluster as input and fills several histograms
 *
 * Histograms (1D):
 *  * - Invariant mass of two clusters
 * 
 * @ingroup alihlt_phos
 */

class AliHLTPHOSHistoProdCellEnergy 
{
 public:
  
  /** Constructor */
  AliHLTPHOSHistoProdCellEnergy();

  /** Destructor */
  virtual ~AliHLTPHOSHistoProdCellEnergy();

  /** Copy constructor */
  AliHLTPHOSHistoProdCellEnergy(const AliHLTPHOSHistoProdCellEnergy &) :
    fClusterReader(NULL),
    fHistCellEnergy(NULL),
    fHistCellEnergyVsNCells(NULL),
    fHistArrayPtr(NULL)
  {
    // Copy constructor not implemented
  }

  /** Assignment operator */
  AliHLTPHOSHistoProdCellEnergy & operator= (const AliHLTPHOSHistoProdCellEnergy)
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
  TH1F *fHistCellEnergy;                 //!transient

  /** 2D histogram of cluster energy vs the number of cells in the cluster */
  TH2F *fHistCellEnergyVsNCells;
  
  /** Pointer to the array of histograms */
  TObjArray *fHistArrayPtr;                     //!transient

  ClassDef(AliHLTPHOSHistoProdCellEnergy, 0);

};
 
#endif
