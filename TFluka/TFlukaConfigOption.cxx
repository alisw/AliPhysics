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

#include "TFlukaConfigOption.h"
ClassImp(TFlukaConfigOption)


TFlukaConfigOption::TFlukaConfigOption()
{
    // Default constructor
}


TFlukaConfigOption::TFlukaConfigOption(const char* cutName, Double_t cut)
    : TNamed(cutName, "Cut")
{
    // Constructor
    fCutValue        = cut;
    fMedium          = -1;
}

TFlukaConfigOption::TFlukaConfigOption(const char* cutName, Double_t cut, Int_t imed)
    : TNamed(cutName, "Cut")
{
    // Constructor
    fCutValue       = cut;
    fMedium         = imed;
}


TFlukaConfigOption::TFlukaConfigOption(const char* procName, Int_t flag)
    : TNamed(procName, "Process")
{
    // Constructor
    fProcessFlag     = flag;
    fMedium          = -1;
}

TFlukaConfigOption::TFlukaConfigOption(const char* procName, Int_t flag, Int_t imed)
    : TNamed(procName, "Process")
{
    // Constructor
    fProcessFlag    = flag;
    fMedium         = imed;
}
