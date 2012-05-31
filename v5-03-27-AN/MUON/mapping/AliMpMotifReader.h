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

#include "AliMpStationType.h"
#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"

#include <TExMap.h>
#include <TString.h>
#include <Riostream.h>

#include <fstream>

class AliMpMotifMap;
class AliMpVMotif;
class AliMpMotifSpecial;
class AliMpMotifType;
class AliMpDataStreams;

class AliMpMotifReader : public TObject
{
  public:
    AliMpMotifReader(const AliMpDataStreams& dataStreams,
                     AliMp::StationType station, 
                     AliMq::Station12Type station12,
                     AliMp::PlaneType plane);
    virtual ~AliMpMotifReader();
  
    // methods   
    AliMpMotifType*     BuildMotifType(const TString& motifTypeId);
    AliMpMotifSpecial*  BuildMotifSpecial(const TString& motifID,
                                          AliMpMotifType* motifType,
                                          Double_t scale=1.0);
    TString MotifSpecialName(const TString& motifID, Double_t scale);
    
  private:
    /// Not implemented
    AliMpMotifReader();
    /// Not implemented
    AliMpMotifReader(const AliMpMotifReader& right);
    /// Not implemented
    AliMpMotifReader&  operator = (const AliMpMotifReader& right);

    // data members  
    const AliMpDataStreams& fkDataStreams;///< data streams
    AliMp::StationType    fStationType;   ///< station type 
    AliMq::Station12Type  fStation12Type; ///< station12 type 
    AliMp::PlaneType      fPlaneType;     ///< plane type 

  ClassDef(AliMpMotifReader,0)  // Data reader
};

#endif //ALI_MP_MOTIF_READER_H
