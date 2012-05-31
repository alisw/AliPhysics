/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Francesco Blanco                                      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#ifndef ALIHLTEMCALCLUSTERMONITOR_H
#define ALIHLTEMCALCLUSTERMONITOR_H

/**
 * Class makes histograms from information from raw data
 *
 * @file   AliHLTEMCALClusterMonitor.h
 * @author Francesco Blanco
 * @date
 * @brief  Histo maker for EMCAL HLT
 */


#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"


#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TObjArray.h"
#include "TString.h"

class AliHLTEMCALClusterMonitor : public TObject
{

public:
  /** Constructor */
  AliHLTEMCALClusterMonitor();

  /** Destructor */
  virtual ~AliHLTEMCALClusterMonitor();

  Int_t MakeHisto(AliHLTCaloClusterHeaderStruct *caloClusterHeaderPtr);

  TObjArray * GetHistograms();

private:

  AliHLTCaloClusterReader* fClusterReaderPtr; // !transient The reader

  TObjArray *hList;
	TH1F *hClusterEne;
	TH2F *hClusterEneVsTime;
	TH1I *hClusterCells;
	TH2F *hClusterEneVsCells;
	TH2F *hClusterEtaVsPhi;

  AliHLTEMCALClusterMonitor(const AliHLTEMCALClusterMonitor &);
  AliHLTEMCALClusterMonitor & operator = (const AliHLTEMCALClusterMonitor &);
  
  ClassDef(AliHLTEMCALClusterMonitor, 0); 

};


#endif
 
