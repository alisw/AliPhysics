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


#include "AliJetProductionData.h"
#include "AliLog.h"

ClassImp(AliJetProductionData)
 
 
////////////////////////////////////////////////////////////////////////

AliJetProductionData::AliJetProductionData() 
{
  // Constructor
    fNbins           = 0;   
    fPtHardLimits    = 0x0;
    fRunTitles       = 0x0;
} 

////////////////////////////////////////////////////////////////////////

AliJetProductionData::~AliJetProductionData()
{
  // Destructor
    delete fPtHardLimits;
    delete fRunTitles;
    
}

void AliJetProductionData::GetPtHardLimits(Int_t bin, Float_t& ptmin, Float_t& ptmax)
{
// Get pt_hard limits for given bin
    if (bin >= 0 && bin < fNbins) {
	ptmin = fPtHardLimits[bin];
	ptmax = fPtHardLimits[bin + 1];	
    } else {
	AliFatal("Bin out of range !");
    }
}

TString AliJetProductionData::GetRunTitle(Int_t bin)
{
    // Get run title for given bin
    
    if (bin < 0 || bin >= fNbins) 
	AliFatal("Bin out of range !");
    
    return fRunTitles[bin];
}

Float_t  AliJetProductionData::GetWeight(Int_t bin)
{
    // Get weight for given bin
      if (bin < 0 || bin >= fNbins) 
	  AliFatal("Bin out of range !");
      return fWeights[bin];
}
