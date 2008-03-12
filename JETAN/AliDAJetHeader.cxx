//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************

//----------------------------------------------------------------------------
// Deterministic Annealing Jet header class
// Stores parameters of DA jet algorithm
// Author: Davide Perrino (davide.perrino@ba.infn.it, davide.perrino@cern.ch)
//----------------------------------------------------------------------------

#include "AliDAJetHeader.h"

ClassImp(AliDAJetHeader)


//---------------------------------------------------------------------
AliDAJetHeader::AliDAJetHeader():
	fDirectory("/home/perrino/events"),
	fFileOut("jets.root"),
	fPytOnly(kTRUE),
	fPtCut(0.),
	fEtaCut(.9),
	fChgOnly(kTRUE),
	fSelectJets(kTRUE),
	fNclustMax(16),
	fFixedCl(kFALSE),
	fEtMin(10.)
{
    // Constructor
}
