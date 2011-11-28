//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTCALOHISTOCELLENERGY
#define ALIHLTCALOHISTOCELLENERGY
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** 
 * @file   AliHLTCaloHistoCellEnergy
 * @author Svein Lindal <slindal@fys.uio.no>
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

class TObjArray;
class TH1F;
class TH2F;
class TRefArray;
class TString;


/** 
 * @class AliHLTCaloHistoCellEnergy
 *
 * Class produces physics histograms for PHOS. It takes a TClonesArray
 * of AliESDCalocluster as input and fills several histograms
 *
 * Histograms (1D):
 *  * - Invariant mass of two clusters
 * 
 * @ingroup alihlt_phos
 */





class AliHLTCaloHistoCellEnergy 
{
 public:
  
  /** Constructor */
  AliHLTCaloHistoCellEnergy(TString det);

  /** Destructor */
  virtual ~AliHLTCaloHistoCellEnergy();

  /** Analyse the clusters in the event */
  Int_t FillHistograms(Int_t nc, TRefArray * clustersArray);
  
  /** Get a pointer to the TObjArray of histograms */
  TObjArray *GetHistograms();

  
  
 private:
  
  /** Default constructor prohibited */
  AliHLTCaloHistoCellEnergy();

  /** Copy constructor prohibited*/
  AliHLTCaloHistoCellEnergy(const AliHLTCaloHistoCellEnergy &);

  /** Assignment operator prohibited*/
  AliHLTCaloHistoCellEnergy & operator= (const AliHLTCaloHistoCellEnergy &);

  /** Histogram of the 2 cluster invariant mass */
  TH1F *fHistCellEnergy;                 //!transient

  /** 2D histogram of cluster energy vs the number of cells in the cluster */
  TH2F *fHistCellEnergyVsNCells;
  
  /** Pointer to the array of histograms */
  TObjArray *fHistArrayPtr;                     //!transient

  ClassDef(AliHLTCaloHistoCellEnergy, 1);

};
 
#endif
