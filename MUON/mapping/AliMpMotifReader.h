/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifReader.h,v 1.9 2006/05/24 13:58:18 ivana Exp $

/// \ingroup motif
/// \class AliMpMotifReader
/// \brief Class that takes care of reading the motifs data.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_READER_H
#define ALI_MP_MOTIF_READER_H

#include <TObject.h>

#include "AliMpContainers.h"

#include "AliMpStationType.h"
#include "AliMpPlaneType.h"
#include "AliMpIntPair.h"
#include "AliMpContainers.h"

#ifdef WITH_ROOT
#include <TExMap.h>
#endif    
#include <TString.h>
#include <TVector2.h>
#include <Riostream.h>

#include <fstream>
#ifdef WITH_STL
#include <map>
#endif    

class AliMpMotifMap;
class AliMpVMotif;
class AliMpMotifSpecial;
class AliMpMotifType;

class AliMpMotifReader : public TObject
{
  public:
#ifdef WITH_STL
    /// Map of int pair to string
    typedef std::map<std::string, std::pair<Int_t,Int_t> > PadMapType;
    /// Map of int pair to string iterator
    typedef PadMapType::iterator PadMapTypeIterator;
#endif    
#ifdef WITH_ROOT
    /// Map of int pair to string
    typedef TExMap PadMapType;
#endif    

  public:
    AliMpMotifReader(AliMp::StationType station, AliMp::PlaneType plane);
    AliMpMotifReader();
    virtual ~AliMpMotifReader();
  
    // methods   
    AliMpMotifType*     BuildMotifType(const TString& motifTypeId);
    AliMpMotifSpecial*  BuildMotifSpecial(const TString& motifID,
                                          AliMpMotifType* motifType,
                                          Double_t scale=1.0);
    TString MotifSpecialName(const TString& motifID, Double_t scale);
    
  private:
    /// Not implemented
    AliMpMotifReader(const AliMpMotifReader& right);
    /// Not implemented
    AliMpMotifReader&  operator = (const AliMpMotifReader& right);

    // data members  
    AliMp::StationType  fStationType; ///< station type 
    AliMp::PlaneType    fPlaneType;   ///< plane type 

  ClassDef(AliMpMotifReader,0)  // Data reader
};

#endif //ALI_MP_MOTIF_READER_H
