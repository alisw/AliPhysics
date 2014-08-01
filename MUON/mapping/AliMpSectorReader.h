/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSectorReader.h,v 1.7 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpSectorReader
/// \brief Class that takes care of reading the sector data.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SECTOR_READER_H
#define ALI_MP_SECTOR_READER_H

#include <TObject.h>

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpXDirection.h"

#include <TString.h>

#include <fstream>

class AliMpSector;
class AliMpZone;
class AliMpSubZone;
class AliMpRow;
class AliMpVRowSegmentSpecial;
class AliMpMotifReader;
class AliMpVMotif;
class AliMpMotifSpecial;
class AliMpMotifType;
class AliMpDataStreams;

using std::istream;

class AliMpSectorReader : public TObject
{
  public:
    AliMpSectorReader(AliMq::Station12Type station, AliMp::PlaneType plane);
    virtual ~AliMpSectorReader();
  
    // methods   
    AliMpSector*  BuildSector(const AliMpDataStreams& dataStreams);
    
  private:  
    /// Not implemented
    AliMpSectorReader();
    /// Not implemented
    AliMpSectorReader(const AliMpSectorReader& right);
    /// Not implemented
    AliMpSectorReader& operator = (const AliMpSectorReader& right);

     // static methods
    static const TString&  GetSectorKeyword();       
    static const TString&  GetZoneKeyword();         
    static const TString&  GetSubZoneKeyword();      
    static const TString&  GetRowKeyword();          
    static const TString&  GetSectorSpecialKeyword();
    static const TString&  GetMotifKeyword();        
    static const TString&  GetRowSpecialKeyword();   
    static const TString&  GetPadRowsKeyword();      
    static const TString&  GetPadRowSegmentKeyword();
  
   // methods
    void  ReadSectorData(const AliMpDataStreams& dataStreams,
                         istream& in);
    void  ReadZoneData(const AliMpDataStreams& dataStreams,
                         istream& in);
    void  ReadSubZoneData(const AliMpDataStreams& dataStreams,
                         istream& in, AliMpZone* zone);
    void  ReadRowSegmentsData(const AliMpDataStreams& dataStreams,
                         istream& in,
                         AliMpZone* zone, AliMpSubZone* subZone);

    AliMpVMotif*  ReadMotifData(const AliMpDataStreams& dataStreams,
                         istream& in, AliMpZone* zone);
    void  ReadSectorSpecialData(const AliMpDataStreams& dataStreams,
                         istream& in, AliMp::XDirection direction);
    void  ReadMotifsSpecialData(const AliMpDataStreams& dataStreams,
                         istream& in);
    void  ReadRowSpecialData(const AliMpDataStreams& dataStreams,
                         istream& in, AliMp::XDirection direction);
    void  ReadRowSegmentSpecialData(const AliMpDataStreams& dataStreams,
                         istream& in,
                         AliMpVRowSegmentSpecial* segment,
                         AliMp::XDirection direction);

    // data members  
    AliMq::Station12Type  fStationType; ///< station type 
    AliMp::PlaneType      fPlaneType;   ///< plane type 
    AliMpSector*          fSector;      ///< sector
    AliMpMotifReader*     fMotifReader; ///< motif reader

  ClassDef(AliMpSectorReader,0)  // Data reader
};

#endif //ALI_MP_READER_H
