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

/* $Id$*/

#include "TFlukaScoringOption.h"
ClassImp(TFlukaScoringOption);


TFlukaScoringOption::TFlukaScoringOption()
{
    // Default constructor
}


TFlukaScoringOption::TFlukaScoringOption(const char* name, const char* sdum, Int_t npar,  Float_t what[12])
    : TNamed(name, sdum)
{
    // Constructor
    fNpar = npar;
    for (Int_t i = 0; i < 12; i++) fWhat[i] = what[i];
}
