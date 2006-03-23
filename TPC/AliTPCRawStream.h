#ifndef ALITPCRAWSTREAM_H
#define ALITPCRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to TPC digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliAltroRawStream.h"

class AliRawReader;
class AliAltroMapping;

class AliTPCRawStream: public AliAltroRawStream {
  public :
    AliTPCRawStream(AliRawReader* rawReader);
    virtual ~AliTPCRawStream();

    virtual void             Reset();

  protected :
    AliTPCRawStream(const AliTPCRawStream& stream);
    AliTPCRawStream& operator = (const AliTPCRawStream& stream);

    virtual void ApplyAltroMapping();

    AliAltroMapping *fMapping[6];   // Pointers to ALTRO mapping

    ClassDef(AliTPCRawStream, 0)    // base class for reading TPC raw digits
};

#endif
