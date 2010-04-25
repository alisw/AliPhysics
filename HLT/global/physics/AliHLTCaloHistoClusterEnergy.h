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

#ifndef ALIHLTCALOHISTOCLUSTERENERGY
#define ALIHLTCALOHISTOCLUSTERENERGY

/** 
 * @file   AliHLTCaloHistoClusterEnergy
 * @author Svein Lindal slindal@fys.uio.no
 * @date 
 * @brief  Produces histograms of cluster energy distributions
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

//#include "AliHLTPHOSBase.h"

#include "AliHLTCaloHistoProducer.h"
#include "Rtypes.h"
//class TObjArray;
class TH1F;
class TH2F;
//class AliESDEvent;
//class TRefArray;
class TString;



/** 
 * @class AliHLTCaloHistoClusterEnergy
 *
 * Class produces physics histograms for PHOS. It takes a TClonesArray
 * of AliESDCalocluster as input and fills several histograms
 *
 * 
 * 
 * 
 * @ingroup alihlt_phos
 */

class AliHLTCaloHistoClusterEnergy : public AliHLTCaloHistoProducer {
 
public:
  
  /** Constructor */
  AliHLTCaloHistoClusterEnergy(TString det);
  
  /** Destructor */
  virtual ~AliHLTCaloHistoClusterEnergy();

  /** Analyse the clusters in the event and fill histograms */
  Int_t FillHistograms(Int_t nc, TRefArray * clusterArray );
  Int_t FillHistograms(Int_t nc, vector<AliHLTCaloClusterDataStruct*> &cVec);

  template <class T>
  Int_t FillClusterEnergyHistos(T*);

 private:

  /** Default constructor prohibited*/
  AliHLTCaloHistoClusterEnergy();

  /** Copy constructor prohibited */
  AliHLTCaloHistoClusterEnergy(const AliHLTCaloHistoClusterEnergy &);

  /** Assignment operator prohibited */
  AliHLTCaloHistoClusterEnergy & operator= (const AliHLTCaloHistoClusterEnergy);

  /** Histogram of the 2 cluster invariant mass */
  TH1F * fHistClusterEnergy;                 //!transient

  /** 2D histogram of cluster energy vs the number of cells in the cluster */
  TH2F * fHistClusterEnergyVsNCells;         //!transient

/** 2D histogram of cluster energy deposit in eta vs phi */
  TH2F * fHistClusterEnergyDepositEtaPhi;         //!transient
  
  ClassDef(AliHLTCaloHistoClusterEnergy, 0);

};
 
#endif
