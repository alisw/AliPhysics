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

#ifndef ALIHLTPHOSPHYSICSHISTOGRAMPRODUCER_H
#define ALIHLTPHOSPHYSICSHISTOGRAMPRODUCER_H

/** 
 * @file   AliHLTPHOSPhysicsHistogramProducer
 * @author Albin Gaignette
 * @date 
 * @brief  Histogram producer for PHOS HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

//#include "AliHLTPHOSBase.h"

#include "Rtypes.h"
#include "TClonesArray.h"
#include "AliHLTPHOSConstants.h" 

using namespace  PhosHLTConst;

class TObjArray;
class TH1F;
class TH2F;


/** 
 * @class AliHLTPHOSPhysicsHistogramProducer
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



//class AliHLTPHOSPhysicsHistogramProducer : public AliHLTPHOSBase
class AliHLTPHOSPhysicsHistogramProducer 
{
 public:

  /** Constructor */
  AliHLTPHOSPhysicsHistogramProducer();

  /** Destructor */
  virtual ~AliHLTPHOSPhysicsHistogramProducer();

  /** Copy constructor */
  AliHLTPHOSPhysicsHistogramProducer(const AliHLTPHOSPhysicsHistogramProducer &) :
  //    AliHLTPHOSBase(),
    fHistNcls(0),
    fHistEnergy(0),
    fHistTotEnergy(0),
    fHistTwoClusterInvMass(0),
    fHistNcells(0),
    fHistNcellsPercentage(0),
    fHistCellsEnergy(0),
    fHistNclusterTotE(0),
    fHistNcellsSumCells(0),
    fHistArrayPtr(0)
  {
    // Copy constructor not implemented
  }

  /** Assignment operator */
  AliHLTPHOSPhysicsHistogramProducer & operator= (const AliHLTPHOSPhysicsHistogramProducer)
  {
    // assignment
    return *this;
  }

  /** Analyse the clusters in the event */
  Int_t AnalyseClusters(TClonesArray* clusters);

  /** Get a pointer to the TObjArray of histograms */
  TObjArray *GetHistograms();
  
 private:

  /** Histogram of number of clusters */
  TH1F *fHistNcls;                              //!transient

  /** Histogram of the cluster energies */
  TH1F *fHistEnergy;                            //!transient

  /** Histogram of the total energy in PHOS */
  TH1F *fHistTotEnergy;                         //!transient

  /** Histogram of the 2 cluster invariant mass */
  TH1F *fHistTwoClusterInvMass;                 //!transient

  /** Histogram of the number of cells with energy */
  TH1F *fHistNcells;                            //!transient

  /** Histogram of the fraction of cells with energy deposit */
  TH1F *fHistNcellsPercentage;                  //!transient
  
  /** Histogram of number of cells in cluster vs cluster energy */
  TH2F *fHistCellsEnergy;                       //!transient

  /** Histogram of number of clusters vs total energy */
  TH2F *fHistNclusterTotE;                      //!transient
  
  /** Histogram of number of cells in the cluster vs total number of cells in PHOS */
  TH2F *fHistNcellsSumCells;                    //!transient

  /** Pointer to the array of histograms */
  TObjArray *fHistArrayPtr;                     //!transient

  ClassDef(AliHLTPHOSPhysicsHistogramProducer, 0);

};
 
#endif
