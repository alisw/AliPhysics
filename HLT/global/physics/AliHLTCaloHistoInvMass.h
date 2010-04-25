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


#include "AliHLTCaloHistoProducer.h"
#include "Rtypes.h"

class TObjArray;
class AliHLTCaloClusterDataStruct;
class TRefArray;
class TString;
class TH1F;

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

using std::vector;

class AliHLTCaloHistoInvMass : public AliHLTCaloHistoProducer {

public:
  
  /** Constructor */
  AliHLTCaloHistoInvMass(TString det);

  /** Destructor */
  virtual ~AliHLTCaloHistoInvMass();

  //** Loops of the calo clusters and fills histos
  Int_t FillHistograms(Int_t nc, TRefArray * clusterArray);
  Int_t FillHistograms(Int_t nc, vector<AliHLTCaloClusterDataStruct*> &cVec);
  //  template <class T> Int_t FillHistograms(Int_t nc, vector<T*> clusterVec);


  //   template <class T>
  //   Int_t FillHistograms(Int_t nc, vector<T*> &cVec);
  
private:

  /** Default constructor prohibited */
  AliHLTCaloHistoInvMass();
  
  /** Copy constructor */
  AliHLTCaloHistoInvMass(const AliHLTCaloHistoInvMass &);
  
  /** Assignment operator */
  AliHLTCaloHistoInvMass & operator= (const AliHLTCaloHistoInvMass &);
  
  /** Calculate 2 cluster inv mass and fill histograms */
  Int_t FillInvariantMassHistograms(Int_t nc, Float_t cPos[][3], Float_t cEnergy[]);

  /** Histogram of the 2 cluster invariant mass */
  TH1F *fHistTwoClusterInvMass0;                 //!transient
  TH1F *fHistTwoClusterInvMass1;                 //!transient
  TH1F *fHistTwoClusterInvMass2;                 //!transient
  TH1F *fHistTwoClusterInvMass3;                 //!transient
  TH1F *fHistTwoClusterInvMass4;                 //!transient

  ClassDef(AliHLTCaloHistoInvMass, 0);

};

#endif
