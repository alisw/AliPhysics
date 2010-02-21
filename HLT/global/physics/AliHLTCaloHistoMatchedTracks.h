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

#ifndef ALIHLTCALOHISTOMATCHEDTRACKS_H
#define ALIHLTCALOHISTOMATCHEDTRACKS_H

/** 
 * @file   AliHLTCaloHistoMatchedTracks
 * @author Albin Gaignette and Svein Lindal slindal@fys.uio.no
 * @date 
 * @brief  Produces Invariant mass histograms of PHOS clusters
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "Rtypes.h"
#include "AliHLTCaloHistoProducer.h"

class TRefArray;
class TObjArray;
class TH1F;


/** 
 * @class AliHLTCaloHistoMatchedTracks
 *
 * Class produces physics histograms for PHOS. It takes a TClonesArray
 * of AliESDCalocluster as input and fills several histograms
 *
 * @ingroup alihlt_phos
 */



//class AliHLTCaloHistoMatchedTracks : public AliHLTPHOSBase
class AliHLTCaloHistoMatchedTracks : public AliHLTCaloHistoProducer {

 public:

  /** Constructor */
  AliHLTCaloHistoMatchedTracks(TString det);

  /** Destructor */
  virtual ~AliHLTCaloHistoMatchedTracks();

  /** Loop over cluster data and fill histograms */
  Int_t FillHistograms(Int_t nc, TRefArray * clusterArray);
  Int_t FillHistograms(Int_t nc, vector<AliHLTCaloClusterDataStruct*> &cVec);

  /** Do the actual histogram filling, regardless of clustertype */
  template <class T>
  Int_t FillMatchedTracks(T*);
  
 private:
  
  /** Default constructor prohibited */
  AliHLTCaloHistoMatchedTracks();

  /** Copy constructor prohibited*/
  AliHLTCaloHistoMatchedTracks(const AliHLTCaloHistoMatchedTracks &);

  /** Assignment operator */
  AliHLTCaloHistoMatchedTracks & operator= (const AliHLTCaloHistoMatchedTracks);

  /** Histogram of the 2 cluster invariant mass */
  TH1F *fHistMatchDistance;                  //!transient
  
  /** Histograms of the energy distribution of mached and unmatched clusters */
  TH1F *fHistMatchedEnergy;                 //!transient
  TH1F *fHistUnMatchedEnergy;               //!transient

  ClassDef(AliHLTCaloHistoMatchedTracks, 0);

};
 
#endif
