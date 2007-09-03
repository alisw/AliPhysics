 /**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
// --- AliRoot header files ---
#include "AliEMCALRecParam.h"
#include "AliLog.h"

ClassImp(AliEMCALRecParam)

//-----------------------------------------------------------------------------
// Container of EMCAL reconstruction parameters
// The purpose of this object is to store it to OCDB
// and retrieve it in AliEMCALClusterizerv1
// Author: Yuri Kharlov
//-----------------------------------------------------------------------------

AliEMCALRecParam::AliEMCALRecParam():
  fClusteringThreshold(0.5),fW0(4.5),fMinECut(0.45)
{
  // default reco values
}

//-----------------------------------------------------------------------------
void AliEMCALRecParam::Print(Option_t *) const
{
  printf("AliEMCALRecParam::Print()\n");
  // Print reconstruction parameters to stdout
  AliInfo(Form("Reconstruction parameters:\n fClusteringThreshold=%.3f,\n fW0=%.3f,\n fMinECut=%.3f",
	       fClusteringThreshold,fW0,fMinECut));

}
