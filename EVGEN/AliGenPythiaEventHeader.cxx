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

/*
$Log$
Revision 1.1  2001/07/13 09:34:53  morsch
Event header class for Pythia added.

*/

#include "AliGenPythiaEventHeader.h"
ClassImp(AliGenPythiaEventHeader)


void AliGenPythiaEventHeader::AddJet(Float_t px, Float_t py, Float_t pz, Float_t e)
{
//
//  Add a jet 
//
    if (fNJets < 10) {
	fJets[0][fNJets] = px;
	fJets[1][fNJets] = py;
	fJets[2][fNJets] = pz;
	fJets[3][fNJets] = e;
	fNJets++;
    } else {
	printf("\nWarning: More than 10 jets triggered !!\n");
    }
}

void AliGenPythiaEventHeader::TriggerJet(Int_t i, Float_t p[4])
{
//
// Give back jet #i
//
    if (i >= fNJets) {
	printf("\nWarning: TriggerJet, index out of Range!!\n");
    } else {
	p[0] = fJets[0][i];
	p[1] = fJets[1][i];
	p[2] = fJets[2][i];
	p[3] = fJets[3][i];
    }
}
