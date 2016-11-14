// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIHLTBASEVGRPACCESS_H
#define ALIHLTBASEVGRPACCESS_H

class AliHLTBaseVGRPAccess
{
protected:
    static void SetGRPAccess(AliHLTBaseVGRPAccess* ptr) {fgPtr = ptr;}
    virtual int GetStartTimeInternal() = 0;
    virtual int GetEndTimeInternal() = 0;
    virtual ~AliHLTBaseVGRPAccess() {};

public:
    static int GetStartTime();
    static int GetEndTime();

private:
    static AliHLTBaseVGRPAccess* fgPtr;
};

#endif
