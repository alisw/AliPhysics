// $Id$
// Category: sector
//
// Class AliMpReader
// -------------------
// Class that takes care of reading the sector data.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_READER_H
#define ALI_MP_READER_H

#include <fstream>

#include <TObject.h>
#include <TString.h>

#include "AliMpSectorTypes.h"
#include "AliMpStationType.h"
#include "AliMpPlaneType.h"
#include "AliMpXDirection.h"
#include "AliMpIntPair.h"

class AliMpSector;
class AliMpZone;
class AliMpSubZone;
class AliMpRow;
class AliMpVRowSegmentSpecial;
class AliMpVMotif;
class AliMpMotifSpecial;
class AliMpMotifType;

class AliMpReader : public TObject
{
  public:
    AliMpReader(AliMpStationType station, AliMpPlaneType plane);
    AliMpReader();
    virtual ~AliMpReader();
  
    // methods   
    AliMpSector*        BuildSector();
    AliMpMotifType*     BuildMotifType(const TString& motifTypeId);
    AliMpMotifSpecial*  BuildMotifSpecial(const TString& motifID, 
                                          AliMpMotifType* motifType);

    // set methods
    void SetVerboseLevel(Int_t verboseLevel); 
    
  private:
#ifdef WITH_ROOT
    static const Int_t   fgkSeparator;  // the separator used for conversion
                                        // of string to Int_t
    
    // methods
    Int_t  GetIndex(const std::string& s) const;
    Int_t  GetIndex(const AliMpIntPair& pair) const;
    std::string  GetString(Int_t index) const;
    AliMpIntPair  GetPair(Int_t index) const;
#endif
  
    // methods
    void  ReadSectorData(ifstream& in);
    void  ReadZoneData(ifstream& in);
    void  ReadSubZoneData(ifstream& in, AliMpZone* zone);
    void  ReadRowSegmentsData(ifstream& in,
                          AliMpZone* zone, AliMpSubZone* subZone);
    AliMpVMotif*  ReadMotifData(ifstream& in, AliMpZone* zone);
    void  ReadSectorSpecialData(ifstream& in, AliMpXDirection direction);
    void  ReadMotifsSpecialData(ifstream& in);
    void  ReadRowSpecialData(ifstream& in, AliMpXDirection direction);
    void  ReadRowSegmentSpecialData(ifstream& in,
                          AliMpVRowSegmentSpecial* segment,
			  AliMpXDirection direction);
    
    // static data members
    static const TString  fgkSectorKeyword;        // sector keyword
    static const TString  fgkZoneKeyword;          // zone keyword
    static const TString  fgkSubZoneKeyword;       // subzone keyword
    static const TString  fgkRowKeyword;           // row keyword
    static const TString  fgkEofKeyword;           // eof keyword
    static const TString  fgkSectorSpecialKeyword; // sector special keyword
    static const TString  fgkMotifKeyword;         // motif keyword
    static const TString  fgkRowSpecialKeyword;    // row special keyword
    static const TString  fgkPadRowsKeyword;       // pad rows keyword
    static const TString  fgkPadRowSegmentKeyword; // pad row segment keyword
  
    // data members  
    AliMpStationType  fStationType; // station type 
    AliMpPlaneType    fPlaneType;   // plane type 
    AliMpSector*      fSector;      // sector
    Int_t             fVerboseLevel;// verbose level

  ClassDef(AliMpReader,1)  // Data reader
};

#endif //ALI_MP_READER_H
