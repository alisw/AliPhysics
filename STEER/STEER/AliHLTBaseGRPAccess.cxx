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

#include "AliHLTBaseGRPAccess.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"

bool AliHLTBaseGRPAccess::fgGrpInitialized = false;
const AliGRPObject* AliHLTBaseGRPAccess::fgGrpObj = NULL;

void AliHLTBaseGRPAccess::initGRP()
{
    if (fgGrpInitialized == false)
    {
        static AliGRPManager grp;
        grp.ReadGRPEntry();
        const AliGRPObject *grpObjTmp = grp.GetGRPData();
        fgGrpObj = grpObjTmp;
        fgGrpInitialized = true;
    }
}

int AliHLTBaseGRPAccess::GetStartTimeInternal()
{
    initGRP();
    int retVal = fgGrpObj->GetTimeStart();
    return(retVal);
}

int AliHLTBaseGRPAccess::GetEndTimeInternal()
{
    initGRP();
    int retVal = fgGrpObj->GetTimeEnd();
    return(retVal);
}

static AliHLTBaseVGRPAccess* gAliHLTBaseGRPAccessGlobal = AliHLTBaseGRPAccess::Instance();
