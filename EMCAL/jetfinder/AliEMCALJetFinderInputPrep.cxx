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

//_________________________________________________________________________
//  Base Class for JetFinder Input Preparation
// --
//*-- Author: Mark Horner (LBL/UCT)
// --



#include <stdio.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TTree.h>
#include "Riostream.h"
//.....................
#include "AliEMCALJet.h"
#include "AliEMCALParton.h"
#include "AliEMCALJetFinderInputPrep.h"
#include "AliEMCALJetFinderInput.h"

ClassImp(AliEMCALJetFinderInputPrep)

//________________________________________________________________________
AliEMCALJetFinderInputPrep::AliEMCALJetFinderInputPrep() 
  : fDebug(0),fInputObject(),fPythiaComparison(0)
{
	// Default constructor
if (fDebug > 0) Info("AliEMCALJetFinderInputPrep","Beginning Constructor");	
  fDebug = 0;
  fPythiaComparison = 0; // This requires lots of checks 
  fInputObject.SetDebug(0);
}
AliEMCALJetFinderInputPrep::~AliEMCALJetFinderInputPrep()
{

}

Int_t AliEMCALJetFinderInputPrep::FillFromFile(TString * /*filename*/,
					       AliEMCALJetFinderFileType_t /*fileType*/,
					       Int_t /*EventNumber*/)
{
  return 0;
}
void AliEMCALJetFinderInputPrep::Reset(AliEMCALJetFinderResetType_t resettype)
{ 
	//  Reset data
if (fDebug > 1) Info("Reset","Beginning Reset");
	
	switch (resettype){
		case kResetData:
			fInputObject.Reset(resettype);
			break;
		case kResetTracks:
			fInputObject.Reset(kResetTracks);
			break;
		case kResetDigits:
			fInputObject.Reset(kResetDigits);
			break;
		case kResetParameters:
			break;
		case kResetAll:
			fInputObject.Reset(kResetAll);
			break;
	        case kResetPartons:
		  Warning("FillFromFile", "kResetPartons not implemented") ; 
		  break;
	       case kResetParticles:
		 Warning("FillFromFile", "kResetParticles not implemented") ; 
		 break;
	       case kResetJets:
		 Warning("FillFromFile", "kResetJets not implemented") ; 
		 break;
	}// end switch

}



