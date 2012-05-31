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

/* $Id: AliAnalysisTaskBaseLine.cxx 46301 2011-01-06 14:25:27Z agheata $ */

//
//
// This empty task is used for the analysis train to estimate the memory and CPU consumption without any user code
//
//

#include "AliAnalysisTaskBaseLine.h"

ClassImp(AliAnalysisTaskBaseLine)

AliAnalysisTaskBaseLine::AliAnalysisTaskBaseLine() 
   :AliAnalysisTaskSE()
{
}

AliAnalysisTaskBaseLine::AliAnalysisTaskBaseLine(const char *name)
   :AliAnalysisTaskSE(name)
{
}

AliAnalysisTaskBaseLine::~AliAnalysisTaskBaseLine()
{
}

void AliAnalysisTaskBaseLine::UserExec(Option_t *) 
{
}
