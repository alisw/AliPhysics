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


//---------------------------------------------------------------------
// Service class for jet production data 
// Author: Andreas Morsch
// andreas.morsch@cern.ch
//---------------------------------------------------------------------


#include "AliJetProductionDataPDC2004.h"
ClassImp(AliJetProductionDataPDC2004)
 
 
////////////////////////////////////////////////////////////////////////

AliJetProductionDataPDC2004::AliJetProductionDataPDC2004() 
{
  // Constructor
    const Float_t kPthard[13] = 
	{20., 24., 29., 35., 42., 50., 60., 72., 86., 104., 125., 150., 180.};

    const TString kRunTitle[12] = 
	{"20-24GeV",  "24-29GeV",   "29-35GeV",   "35-42GeV", "42-50GeV", "50-60GeV", "60-72GeV", "72-86GeV",
	 "86-104GeV", "104-125GeV", "125-150GeV", "150-180GeV"};
    
    const Float_t kXsection[12]  =
	{1.10e-1, 6.00e-2, 2.00e-2, 1.16e-2, 5.33e-3, 2.68e-3, 1.22e-3, 5.26e-4, 2.53e-4, 1.00e-4, 4.44e-5, 1.84e-5};   

    const Float_t kGeneratedOverSelected[12] = 
	{2.67,   2.68,    2.47,    2.5,    2.3,    2.2,    2.1,    2.0,    1.9,    1.90,    1.7,    1.61};
    
	    
    fNbins     = 12;   
    fPtHardLimits    = new Float_t[fNbins + 1];
    fRunTitles       = new TString[fNbins];
    fWeights         = new Float_t[fNbins];
    Int_t i;
    
    for (i = 0; i < fNbins +1; i++) {
	fPtHardLimits[i] = kPthard[i];
    }

    for (i = 0; i < fNbins; i++) {
	fRunTitles[i] = kRunTitle[i];
	fWeights[i]   = kXsection[i] /  kGeneratedOverSelected[i];
    }
} 

////////////////////////////////////////////////////////////////////////

AliJetProductionDataPDC2004::~AliJetProductionDataPDC2004()
{
  // Destructor
}
