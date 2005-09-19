/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $Id$

/// \ingroup motif
/// \class AliMpMotifReader
/// \brief Class that takes care of reading the motifs data.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_READER_H
#define ALI_MP_MOTIF_READER_H

#include <fstream>

#include <TObject.h>
#include <TString.h>
#include <TVector2.h>
#include <Riostream.h>

#include "AliMpStationType.h"
#include "AliMpPlaneType.h"
#include "AliMpIntPair.h"
#include "AliMpContainers.h"

class AliMpMotifMap;
class AliMpVMotif;
class AliMpMotifSpecial;
class AliMpMotifType;

class AliMpMotifReader : public TObject
{
  public:
    AliMpMotifReader(AliMpStationType station, AliMpPlaneType plane);
    AliMpMotifReader();
    virtual ~AliMpMotifReader();
  
    // methods   
    AliMpMotifType*     BuildMotifType(const TString& motifTypeId);
    AliMpMotifSpecial*  BuildMotifSpecial(const TString& motifID,
                                          AliMpMotifType* motifType);

    // set methods
    void SetVerboseLevel(Int_t verboseLevel); 
    
  protected:
    AliMpMotifReader(const AliMpMotifReader& right);
    AliMpMotifReader&  operator = (const AliMpMotifReader& right);

  private:
#ifdef WITH_ROOT
    static const Int_t   fgkSeparator;  // the separator used for conversion
                                        // of string to Int_t
    
    // methods
    Int_t   GetIndex(const string& s) const;
    Int_t   GetIndex(const AliMpIntPair& pair) const;
    string  GetString(Int_t index) const;
    AliMpIntPair  GetPair(Int_t index) const;
#endif

    // data members  
    AliMpStationType  fStationType; // station type 
    AliMpPlaneType    fPlaneType;   // plane type 
    Int_t             fVerboseLevel;// verbose level

  ClassDef(AliMpMotifReader,1)  // Data reader
};

#endif //ALI_MP_MOTIF_READER_H
