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

#include "AliMpStationType.h"
#include "AliMpPlaneType.h"
#include "AliMpXDirection.h"
#include "AliMpIntPair.h"

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

class AliMpSectorReader : public TObject
{
  public:
    AliMpSectorReader(AliMp::StationType station, AliMp::PlaneType plane);
    AliMpSectorReader();
    virtual ~AliMpSectorReader();
  
    // methods   
    AliMpSector*  BuildSector();
    
  private:  
    AliMpSectorReader(const AliMpSectorReader& right);
    AliMpSectorReader&  operator = (const AliMpSectorReader& right);

    // methods
    void  ReadSectorData(ifstream& in);
    void  ReadZoneData(ifstream& in);
    void  ReadSubZoneData(ifstream& in, AliMpZone* zone);
    void  ReadRowSegmentsData(ifstream& in,
                          AliMpZone* zone, AliMpSubZone* subZone);
    AliMpVMotif*  ReadMotifData(ifstream& in, AliMpZone* zone);
    void  ReadSectorSpecialData(ifstream& in, AliMp::XDirection direction);
    void  ReadMotifsSpecialData(ifstream& in);
    void  ReadRowSpecialData(ifstream& in, AliMp::XDirection direction);
    void  ReadRowSegmentSpecialData(ifstream& in,
                          AliMpVRowSegmentSpecial* segment,
			  AliMp::XDirection direction);
    
    // static data members
    static const TString  fgkSectorKeyword;        ///< sector keyword
    static const TString  fgkZoneKeyword;          ///< zone keyword
    static const TString  fgkSubZoneKeyword;       ///< subzone keyword
    static const TString  fgkRowKeyword;           ///< row keyword
    static const TString  fgkEofKeyword;           ///< eof keyword
    static const TString  fgkSectorSpecialKeyword; ///< sector special keyword
    static const TString  fgkMotifKeyword;         ///< motif keyword
    static const TString  fgkRowSpecialKeyword;    ///< row special keyword
    static const TString  fgkPadRowsKeyword;       ///< pad rows keyword
    static const TString  fgkPadRowSegmentKeyword; ///< pad row segment keyword
  
    // data members  
    AliMp::StationType  fStationType; ///< station type 
    AliMp::PlaneType    fPlaneType;   ///< plane type 
    AliMpSector*        fSector;      ///< sector
    AliMpMotifReader*   fMotifReader; ///< motif reader

  ClassDef(AliMpSectorReader,1)  // Data reader
};

#endif //ALI_MP_READER_H
