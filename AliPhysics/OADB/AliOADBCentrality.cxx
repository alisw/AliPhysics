/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-------------------------------------------------------------------------
//     OADB class for run dependent centrality scaling and 
//     data for centrality determination
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include "AliOADBCentrality.h"
ClassImp(AliOADBCentrality);

//______________________________________________________________________________
AliOADBCentrality::AliOADBCentrality() : 
  TNamed(),
  fV0MScaleFactor(1.),
  fSPDScaleFactor(1.),
  fTPCScaleFactor(1.),
  fV0MScaleFactorMC(1.),
  fV0MSPDOutlierPar0(1.),
  fV0MSPDOutlierPar1(1.),
  fV0MTPCOutlierPar0(1.),
  fV0MTPCOutlierPar1(1.),
  fV0MSPDSigmaOutlierPar0(1.),
  fV0MSPDSigmaOutlierPar1(1.),
  fV0MSPDSigmaOutlierPar2(1.),
  fV0MTPCSigmaOutlierPar0(1.),
  fV0MTPCSigmaOutlierPar1(1.),
  fV0MTPCSigmaOutlierPar2(1.),
  fV0MZDCOutlierPar0(1.),
  fV0MZDCOutlierPar1(1.),
  fV0MZDCEcalOutlierPar0(1.),
  fV0MZDCEcalOutlierPar1(1.),
  fZVCut(10.),
  fOutliersCut(6.),
  fUseScaling(kFALSE),
  fUseCleaning(kTRUE),
  f1DHistos(),
  f2DHistos()
{
  // Default constructor
}
//______________________________________________________________________________
AliOADBCentrality::AliOADBCentrality(char* name) : 
  TNamed(name, "Centrality Scaling"),
  fV0MScaleFactor(1.),
  fSPDScaleFactor(1.),
  fTPCScaleFactor(1.),
  fV0MScaleFactorMC(1.),
  fV0MSPDOutlierPar0(1.),
  fV0MSPDOutlierPar1(1.),
  fV0MTPCOutlierPar0(1.),
  fV0MTPCOutlierPar1(1.),
  fV0MSPDSigmaOutlierPar0(1.),
  fV0MSPDSigmaOutlierPar1(1.),
  fV0MSPDSigmaOutlierPar2(1.),
  fV0MTPCSigmaOutlierPar0(1.),
  fV0MTPCSigmaOutlierPar1(1.),
  fV0MTPCSigmaOutlierPar2(1.),
  fV0MZDCOutlierPar0(1.),
  fV0MZDCOutlierPar1(1.),
  fV0MZDCEcalOutlierPar0(1.),
  fV0MZDCEcalOutlierPar1(1.),
  fZVCut(10.),
  fOutliersCut(6.),
  fUseScaling(kFALSE),
  fUseCleaning(kTRUE),
  f1DHistos(),
  f2DHistos()
{
  // Constructor
}
//______________________________________________________________________________
AliOADBCentrality::~AliOADBCentrality() 
{
  // destructor
}
