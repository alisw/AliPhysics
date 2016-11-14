// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIHLTBASEGRPACCESS_H
#define ALIHLTBASEGRPACCESS_H

#include "AliHLTBaseVGRPAccess.h"
#include <stdio.h>

class AliGRPObject;

class AliHLTBaseGRPAccess : public AliHLTBaseVGRPAccess
{
private:
    virtual int GetStartTimeInternal();
    virtual int GetEndTimeInternal();

    AliHLTBaseGRPAccess()
    {
        AliHLTBaseVGRPAccess::SetGRPAccess(this);
    }

    static bool fgGrpInitialized;
    static const AliGRPObject *fgGrpObj;
    static void initGRP();

protected:
    virtual ~AliHLTBaseGRPAccess() {};

public:
    static AliHLTBaseGRPAccess* Instance()
    {
        static AliHLTBaseGRPAccess instance;
        return(&instance);
    }
};

#endif
