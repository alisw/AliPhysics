/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifReader.h,v 1.5 2005/09/26 16:10:46 ivana Exp $

/// \ingroup motif
/// \class AliMpMotifReader
/// \brief Class that takes care of reading the motifs data.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_READER_H
#define ALI_MP_MOTIF_READER_H

#include "AliMpContainers.h"

#ifdef WITH_STL
#include <map>
#endif    

#ifdef WITH_ROOT
#include <TExMap.h>
#endif    

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
#ifdef WITH_STL
    typedef std::map<std::string, std::pair<Int_t,Int_t> > PadMapType;
    typedef PadMapType::iterator PadMapTypeIterator;
#endif    
#ifdef WITH_ROOT
    typedef TExMap PadMapType;
#endif    

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
    // data members  
    AliMpStationType  fStationType; // station type 
    AliMpPlaneType    fPlaneType;   // plane type 
    Int_t             fVerboseLevel;// verbose level

  ClassDef(AliMpMotifReader,1)  // Data reader
};

#endif //ALI_MP_MOTIF_READER_H
