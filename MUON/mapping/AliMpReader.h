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
#include <TVector2.h>

#include "AliMpPlaneType.h"

class AliMpSector;
class AliMpZone;
class AliMpSubZone;
class AliMpRow;
class AliMpRowSegmentSpecial;
class AliMpVMotif;
class AliMpMotifSpecial;
class AliMpMotifType;

class AliMpReader : public TObject
{
  public:
    AliMpReader(AliMpPlaneType plane);
    AliMpReader();
    virtual ~AliMpReader();
  
    // methods   
    AliMpSector*        BuildSector();
    AliMpMotifType*     BuildMotifType(TString motifTypeId);
    AliMpMotifSpecial*  BuildMotifSpecial(TString motifID, 
                                          AliMpMotifType* motifType);

    // set methods
    void SetVerboseLevel(Int_t verboseLevel); 
    
  private:
    // methods
    void  ReadSectorData(ifstream& in);
    void  ReadZoneData(ifstream& in);
    void  ReadSubZoneData(ifstream& in, AliMpZone* zone);
    void  ReadRowSegmentsData(ifstream& in,
                          AliMpZone* zone, AliMpSubZone* subZone);
    AliMpVMotif*  ReadMotifData(ifstream& in, AliMpZone* zone);
    void  ReadSectorSpecialData(ifstream& in);
    void  ReadMotifsSpecialData(ifstream& in);
    void  ReadRowSpecialData(ifstream& in);
    void  ReadRowSegmentSpecialData(ifstream& in,
                          AliMpRowSegmentSpecial* segment);
    
    // static data members
    static const TString  fgkSectorKeyword;
    static const TString  fgkZoneKeyword;
    static const TString  fgkSubZoneKeyword;
    static const TString  fgkRowKeyword;
    static const TString  fgkEofKeyword;
    static const TString  fgkSectorSpecialKeyword;
    static const TString  fgkMotifKeyword;
    static const TString  fgkRowSpecialKeyword;
    static const TString  fgkPadRowsKeyword;
    static const TString  fgkPadRowSegmentKeyword;
  
    // data members    
    AliMpPlaneType  fPlaneType;   // plane type 
    AliMpSector*    fSector;      // sector
    Int_t           fVerboseLevel;// verbose level

  ClassDef(AliMpReader,1)  // Data reader
};

#endif //ALI_MP_READER_H

