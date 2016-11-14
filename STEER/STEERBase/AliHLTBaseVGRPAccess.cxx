// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: David Rohr <drohr@cern.ch> for The ALICE HLT Project.   *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************

#include "AliHLTBaseVGRPAccess.h"
#include <stddef.h>
#include <stdio.h>

AliHLTBaseVGRPAccess* AliHLTBaseVGRPAccess::fgPtr; //Auto initialized to NULL, will be overriden by STEER:AliHLTBaseGRPAccess

int AliHLTBaseVGRPAccess::GetStartTime() {
    return(fgPtr ? fgPtr->GetStartTimeInternal() : 0);
}
int AliHLTBaseVGRPAccess::GetEndTime() {
    return(fgPtr ? fgPtr->GetEndTimeInternal() : 0);
}
