#ifndef ALI_TRD_PREPROCESSOR_H
#define ALI_TRD_PREPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TRD preprocessor for the database SHUTTLE                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliPreprocessor.h"

class AliTRDPreprocessor : public AliPreprocessor
{

  public:

    AliTRDPreprocessor(AliShuttleInterface *shuttle);
    virtual ~AliTRDPreprocessor();

  protected:

    virtual void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* /*dcsAliasMap*/);

  private:
    
    ClassDef(AliTRDPreprocessor,0);

};

#endif
