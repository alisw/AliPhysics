// $Id$
// Category: sector
//
// Class AliMpReader
// -------------------
// Class that takes care of reading the sector data.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#if !defined(__HP_aCC) && !defined(__alpha)
  #include <sstream>
#endif

#include <Riostream.h>
#include <Rstrstream.h>
#include <TSystem.h>
#include <TError.h>
#include <TMath.h>

#include "AliMpReader.h"
#include "AliMpSector.h"
#include "AliMpFiles.h"
#include "AliMpZone.h"
#include "AliMpSubZone.h"
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpRowSegment.h"
#include "AliMpRowSegmentLSpecial.h"
#include "AliMpRowSegmentRSpecial.h"
#include "AliMpPadRow.h"
#include "AliMpMotifMap.h"
#include "AliMpMotif.h"
#include "AliMpMotifSpecial.h"
#include "AliMpMotifType.h"
#include "AliMpConnection.h"
#include "AliMpIntPair.h"
#include "AliMpDirection.h"

ClassImp(AliMpReader)

#ifdef WITH_ROOT
const Int_t    AliMpReader::fgkSeparator = 100;
#endif

const TString  AliMpReader::fgkSectorKeyword  = "SECTOR_DATA";
const TString  AliMpReader::fgkZoneKeyword    = "ZONE";
const TString  AliMpReader::fgkSubZoneKeyword = "SUBZONE";
const TString  AliMpReader::fgkRowKeyword     = "ROW_SEGMENT";
const TString  AliMpReader::fgkEofKeyword     = "EOF";
const TString  AliMpReader::fgkSectorSpecialKeyword  = "SECTOR_SPECIAL_DATA";
const TString  AliMpReader::fgkMotifKeyword          = "MOTIF";
const TString  AliMpReader::fgkRowSpecialKeyword     = "ROW";
const TString  AliMpReader::fgkPadRowsKeyword        = "PAD_ROWS";
const TString  AliMpReader::fgkPadRowSegmentKeyword  = "PAD_ROW_SEGMENT";

//_____________________________________________________________________________
AliMpReader::AliMpReader(AliMpStationType station, AliMpPlaneType plane) 
  : TObject(),
    fStationType(station),
    fPlaneType(plane),
    fSector(0),
    fVerboseLevel(0)
{
//
}

//_____________________________________________________________________________
AliMpReader::AliMpReader() 
  : TObject(),
    fStationType(kStation1),
    fPlaneType(kBendingPlane),
    fSector(0),
    fVerboseLevel(0)
{
//
}

//_____________________________________________________________________________
AliMpReader::AliMpReader(const AliMpReader& right) 
  : TObject(right) {
// 
  Fatal("AliMpReader", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpReader::~AliMpReader() {
//  
}

//
// operators
//

//_____________________________________________________________________________
AliMpReader& AliMpReader::operator=(const AliMpReader& right)
{
  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//
// private methods
//

#ifdef WITH_ROOT
//_____________________________________________________________________________
Int_t  AliMpReader::GetIndex(const string& s) const 
{
// Converts the TString to integer.
// ---

  if (s.length() > 5) {
    Fatal("GetIndex", "String too long.");
    return 0;
  }  

  Int_t index = 0;
  for (Int_t i=s.length(); i>=0; --i)  index = index*100 + int(s[i]);
  
  return index;
}

//______________________________________________________________________________
Int_t  AliMpReader::GetIndex(const AliMpIntPair& pair) const
{
// Converts the pair of integers to integer.
// ---

  if (pair.GetFirst() >= fgkSeparator || pair.GetSecond() >= fgkSeparator)
    Fatal("GetIndex", "Index out of limit.");
      
  return pair.GetFirst()*fgkSeparator + pair.GetSecond() + 1;
}  

//_____________________________________________________________________________
string  AliMpReader::GetString(Int_t index) const
{
// Converts the integer index to the string.
// ---

  string s;
  while (index >0) {
    Char_t c = index%100;
    s += c;
    index = index/100;
  }
  
  return s;
  
}

//______________________________________________________________________________
AliMpIntPair  AliMpReader::GetPair(Int_t index) const
{
// Converts the integer index to the pair of integers.
// ---

  return AliMpIntPair((index-1)/fgkSeparator, (index-1)%fgkSeparator);
}  
#endif

//_____________________________________________________________________________
void  AliMpReader::ReadSectorData(ifstream& in)
{
// Reads sector input data;
// prepares zones and rows vectors to be filled in.
// ---

  TString keyword;
  in >> keyword;
  
  if (fVerboseLevel>0) 
    cout << keyword << endl;

  if (keyword != fgkSectorKeyword) {
     Fatal("ReadSectorData", "Wrong file format.");
     return;
  }   
    
  Int_t nofZones, nofRows;
  TString directionStr;
  in >> nofZones;
  in >> nofRows;
  in >> directionStr;
  
  AliMpDirection direction;
  direction = (directionStr == "Y") ? kY  :  kX;
  if (fVerboseLevel>0) 
     cout << nofZones << " " <<  nofRows << endl;

  fSector = new AliMpSector("Not defined", nofZones, nofRows,direction);
  
  TString nextKeyword;
  in >> nextKeyword;
    
  if (nextKeyword != fgkZoneKeyword) {
     Fatal("ReadSectorData", "Wrong file format.");
     return;
  }      
  
  ReadZoneData(in);
}  

//_____________________________________________________________________________
void AliMpReader::ReadZoneData(ifstream& in)
{
// Reads zone input data;
// creates zone and adds it to zones vector.
// ---

  Int_t zoneID;
  Double_t  sizex, sizey; 
  in >> zoneID;    
  in >> sizex;
  in >> sizey;
  if (fVerboseLevel>0) 
     cout << fgkZoneKeyword << " " <<  zoneID << "  "        
          << sizex << " " << sizey << endl;
  
  AliMpZone* zone =  fSector->GetZone(zoneID);
  zone->SetPadDimensions(TVector2(sizex/2.,sizey/2.)); 
      
  TString nextKeyword;
  in >> nextKeyword;
    
  if (nextKeyword != fgkSubZoneKeyword) {
     Fatal("ReadZoneData", "Wrong file format.");
     return;
  }  
    
  ReadSubZoneData(in, zone);
}

//_____________________________________________________________________________
void AliMpReader::ReadSubZoneData(ifstream& in, AliMpZone* zone)
{
// Reads subzone input data;
// creates subzone and its to the specified zone.
// ---

  if (fVerboseLevel>0) 
    cout << fgkSubZoneKeyword << " ";

  AliMpVMotif* motif = ReadMotifData(in, zone);
  AliMpSubZone* subZone = new AliMpSubZone(motif); 
  zone->AddSubZone(subZone); 

  TString nextKeyword;
  in >> nextKeyword;
    
  if (nextKeyword != fgkRowKeyword) {
     Fatal("ReadSubZoneData", "Wrong file format.");
     return;
  }  
    
  ReadRowSegmentsData(in, zone, subZone);
}   

//_____________________________________________________________________________
AliMpVMotif*  AliMpReader::ReadMotifData(ifstream& in, AliMpZone* zone)
{
// Reads the motif input data.
// ---

  TString  motifID;
  TString  motifTypeID;
  in >> motifID;
  in >> motifTypeID;
  if (fVerboseLevel>0) {
    cout << motifID << " " 
         << motifTypeID << endl;
  }	 
  
  AliMpMotifMap* motifMap = fSector->GetMotifMap();

  AliMpMotifType* motifType = 0;
  AliMpVMotif* motif 
    = motifMap->FindMotif(motifID, motifTypeID, zone->GetPadDimensions());
  if (!motif) {    
    motifType = motifMap->FindMotifType(motifTypeID);
    if (!motifType) {
      motifType = BuildMotifType(motifTypeID);     
      motifMap->AddMotifType(motifType);
    }
    
    if (zone->GetPadDimensions().X() != 0. && zone->GetPadDimensions().Y() != 0.) 
      motif = new AliMpMotif(motifID, motifType, zone->GetPadDimensions());
    else 
      motif = BuildMotifSpecial(motifID, motifType);
      
    if (motif) 
      motifMap->AddMotif(motif);

  }  

  return motif;   
}  

//_____________________________________________________________________________
void AliMpReader::ReadRowSegmentsData(ifstream& in, 
                                      AliMpZone* zone, AliMpSubZone* subZone)
{
// Reads row segments input data of a specified zone and subzone;
// creates row segment and adds it to the specified subzone
// and a corresponding row in the rows vector.
// ---

  TString nextKeyword;
  do {
    // 
    // Read data from file
    //
    Int_t offX, offY, inRow, nofMotifs, firstMotifPositionId, firstMotifPositionDId;
    in >> offX;
    in >> offY;
    in >> inRow;
    in >> nofMotifs;
    in >> firstMotifPositionId;
    in >> firstMotifPositionDId;
    if (fVerboseLevel>0) 
       cout << fgkRowKeyword << " " 
            << offX << " " << offY << " " << inRow << " " << nofMotifs << " " 
	    << firstMotifPositionId << " " << firstMotifPositionDId
	    << endl;

    in >> nextKeyword;

    //
    // Process data
    //
    AliMpRow* row = fSector->GetRow(inRow);
    AliMpVMotif* motif = subZone->GetMotif();
    
    // Create row segment and add it to its zone, row   
    AliMpVRowSegment* rowSegment
      = new AliMpRowSegment(row, motif, AliMpIntPair(offX, offY), nofMotifs, 
                        firstMotifPositionId, firstMotifPositionDId);
			
    subZone->AddRowSegment(rowSegment);
    row->AddRowSegment(rowSegment);
  }
  while (!in.eof() && (nextKeyword == fgkRowKeyword));

  if (in.eof()) return;

  if (nextKeyword == fgkZoneKeyword) {
    ReadZoneData(in);
  }
  else if (nextKeyword == fgkSubZoneKeyword) {
    ReadSubZoneData(in, zone);
  }   
  else {
    Fatal("ReadRowSegmentsData", "Wrong file format.");
  } 
}   

//_____________________________________________________________________________
void AliMpReader::ReadSectorSpecialData(ifstream& in, AliMpXDirection direction)
{
// Reads sector input data
// with a special (irregular) motifs.
// ---

  TString keyword;
  in >> keyword;
  if (fVerboseLevel>0) 
    cout << keyword << endl;

  if (keyword != fgkSectorSpecialKeyword) {
     Fatal("ReadSectorSpecialData", "Wrong file format.");
     return;
  }   

  TString nextKeyword;
  in >> nextKeyword;
  if (fVerboseLevel>0) 
    cout << keyword << endl;
    
  if (nextKeyword != fgkMotifKeyword) {
    Fatal("ReadSectorSpecialData", "Wrong file format.");
    return;
  }  

  ReadMotifsSpecialData(in);     
  ReadRowSpecialData(in, direction); 
}  

//_____________________________________________________________________________
void AliMpReader::ReadMotifsSpecialData(ifstream& in)
{
// Reads the special (irregular) motifs input data.
// ---

  if (fVerboseLevel>0) 
    cout << fgkMotifKeyword << " ";

  TString nextKeyword;
  do {
    Int_t zone;
    in >> zone;
    AliMpVMotif* motif =  ReadMotifData(in, fSector->GetZone(zone));
    AliMpSubZone* subZone = new AliMpSubZone(motif); 
    fSector->GetZone(zone)->AddSubZone(subZone); 
  
    in >> nextKeyword;
    if (fVerboseLevel>0) 
      cout << nextKeyword << " ";      
  }
  while (nextKeyword == fgkMotifKeyword);
    
  if (nextKeyword != fgkRowSpecialKeyword) {
     Fatal("ReadMotifSpecialData", "Wrong file format.");
     return;
  }      
}  

//_____________________________________________________________________________
void AliMpReader::ReadRowSpecialData(ifstream& in, AliMpXDirection direction)
{
// Reads row input data
// with a special (irregular) motifs.
// ---

  Int_t id;
  in >> id;
  if (fVerboseLevel>0) 
      cout << id << endl;      
  
  // Get the row and its border
  AliMpRow* row = fSector->GetRow(id);

  AliMpVRowSegmentSpecial* segment = 0;
  if (direction == kLeft) {
    AliMpVRowSegment* firstNormalSeg = row->GetRowSegment(0);
    Double_t offsetX = firstNormalSeg->LeftBorderX();
  
    // Create a special row segment
    segment = new AliMpRowSegmentLSpecial(row, offsetX);
    row->AddRowSegmentInFront(segment);
  }
  else { 
    AliMpVRowSegment* precedentNormalSeg 
      = row->GetRowSegment(row->GetNofRowSegments()-1);
    Double_t offsetX = precedentNormalSeg->RightBorderX();
  
    // Create a special row segment
    segment = new AliMpRowSegmentRSpecial(row, offsetX);
    row->AddRowSegment(segment);
  }  
      
  TString nextKeyword;
  in >> nextKeyword;
  if (fVerboseLevel>0) 
    cout << nextKeyword << " ";
    
  if (nextKeyword != fgkPadRowsKeyword) {
     Fatal("ReadRowSpecialData", "Wrong file format.");
     return;
  }  
    
  ReadRowSegmentSpecialData(in, segment, direction); 
  
  // Update row segment and set it to all subzones associated with
  // contained motifs
  
  segment->UpdateMotifVector();
  segment->UpdatePadsOffset();
  
  for (Int_t i=0; i<segment->GetNofMotifs(); i++) {
    AliMpSubZone* subZone = 0;
    Int_t j = 0;
    while (!subZone && j<fSector->GetNofZones())
      subZone = fSector->GetZone(++j)->FindSubZone(segment->GetMotif(i));
    
    subZone->AddRowSegment(segment);
  }  
}  

//_____________________________________________________________________________
void AliMpReader::ReadRowSegmentSpecialData(ifstream& in, 
                                            AliMpVRowSegmentSpecial* segment,
					    AliMpXDirection direction)
{
// Reads row segment input data
// with a special (irregular) motifs.
// ---

  Int_t nofPadRows;
  in >> nofPadRows;
  if (fVerboseLevel>0) 
    cout << nofPadRows << endl;
  
  TString keyword;
  in >> keyword;
  if (fVerboseLevel>0) 
    cout << keyword << " ";
    
  if (keyword != fgkPadRowSegmentKeyword) {
     Fatal("ReadRowSegmentSpecialData", "Wrong file format.");
     return;
  }  
  
  //
  // Process data
  //
    
  PadRowVector  newPadRows;
  for (Int_t i=0; i<nofPadRows; i++) {
    
     // Create pad row
     AliMpPadRow* padRow = new AliMpPadRow(direction);
     segment->AddPadRow(padRow);
     
     // Keep the new rows in a temporary vector
#ifdef WITH_STL
     newPadRows.push_back(padRow);
#endif
#ifdef WITH_ROOT
     newPadRows.Add(padRow);
#endif
  }   
      
  TString nextKeyword;
  do {
    // 
    // Read data from file
    //
    Int_t    nofPadsInRow, motifPositionId;
    TString   motifId, motifTypeId;
    in >> nofPadsInRow;
    in >> motifId;
    in >> motifPositionId; 
  
    if (fVerboseLevel>0) 
       cout << nofPadsInRow << " " << motifId << " " << motifPositionId << endl;

    in >> nextKeyword;
    if (fVerboseLevel>0) 
      cout << nextKeyword << " ";

    //
    // Process data
    //
    
    for (Int_t i=0; i<nofPadRows; i++) {
    
      // Get pad row from the temporary vector
#ifdef WITH_STL
      AliMpPadRow* padRow = newPadRows[i];
#endif
#ifdef WITH_ROOT
      AliMpPadRow* padRow = (AliMpPadRow*)newPadRows[i];
#endif
      
      // Find motif
      AliMpVMotif* motif = fSector->GetMotifMap()->FindMotif(motifId);
      
      if (!motif) {
        Fatal("ReadRowSegmentSpecialData", "Unknown motif.");
	return;
      }

      // Create pad row segment
      padRow->AddPadRowSegment(dynamic_cast<AliMpMotif *>(motif), 
                               motifPositionId, nofPadsInRow);
    }  
  }
  while (!in.eof() && (nextKeyword == fgkPadRowSegmentKeyword));
  
  if (in.eof()) return;

  if (nextKeyword == fgkPadRowsKeyword) {
    ReadRowSegmentSpecialData(in, segment, direction);
  }
  else if (nextKeyword == fgkRowSpecialKeyword) {
    ReadRowSpecialData(in, direction);
  }   
  else {
    Fatal("ReadRowSegmentSpecialData", "Wrong file format.");
  } 
}  

//
// public methods
//

//_____________________________________________________________________________
AliMpSector* AliMpReader::BuildSector()
{
// Reads the mapping data from ascii file
// $MINSTALL/data/fileName and creates the basic objects:
// zones, subzones, rows, row segments, motifs.
// ---

  // Open input file
  ifstream in(AliMpFiles::Instance()
              ->SectorFilePath(fStationType, fPlaneType).Data(), ios::in);
  if (!in) {
     cerr << AliMpFiles::Instance()
              ->SectorFilePath(fStationType, fPlaneType) << endl;	
     Error("BuildSector", "File not found.");
     return 0;
  }
  
  ReadSectorData(in);
  fSector->SetRowSegmentOffsets();

  // Open input file for special inner zone
  TString sectorSpecialFileName 
    = AliMpFiles::Instance()->SectorSpecialFilePath(fStationType, fPlaneType);
  if (!gSystem->AccessPathName(sectorSpecialFileName.Data())) {
    ifstream in2(sectorSpecialFileName.Data(), ios::in);
    if (!in2) {	
       cerr << AliMpFiles::Instance()
              ->SectorSpecialFilePath(fStationType, fPlaneType) << endl;	
       Error("BuildSector", "File not found.");
       return 0;
    }
    
    ReadSectorSpecialData(in2, kLeft);
  }   

  // Open input file for special outer zone
  TString sectorSpecialFileName2 
    = AliMpFiles::Instance()->SectorSpecialFilePath2(fStationType, fPlaneType);
  if (!gSystem->AccessPathName(sectorSpecialFileName2.Data())) {
    ifstream in3(sectorSpecialFileName2.Data(), ios::in);
    if (!in3) {	
       cerr << AliMpFiles::Instance()
              ->SectorSpecialFilePath2(fStationType, fPlaneType) << endl;	
       Error("Build", "File not found.");
       return 0;
    }
    
    ReadSectorSpecialData(in3, kRight);
  }   

  fSector->Initialize();
  
  return fSector;
}  

//_____________________________________________________________________________
AliMpMotifType* AliMpReader::BuildMotifType(const TString& motifTypeId)
{

  // Read the files describing a motif in the "$MINSTALL/data" directory
  // and fill the AliMpMotifType structure with.
  // The files mentioned are are named padPos<maskName>.dat
  // and connect<maskName>.dat

  AliMpMotifType*  motifType = new AliMpMotifType(motifTypeId);	

  TString strPadPos 
    = AliMpFiles::Instance()
      ->PadPosFilePath(fStationType, fPlaneType, motifTypeId);
  ifstream padPos(strPadPos.Data());
  if (fVerboseLevel>0) cout<<"Opening file "<<strPadPos<<endl;

  PadMapType positions;

  char line[256];
  do {
    padPos.getline(line,255);
    if (!padPos) break;

#if defined (__HP_aCC) || (__alpha)
    strstream strline;
    strline << line;
#else
    istringstream strline(line);
#endif    
    string key;

    strline>>key;
    if ((key=="#") || (key=="") ) continue;

    int i,j;
    strline>>i>>j;
#ifdef WITH_STL
    positions[key].first=i;
    positions[key].second=j;
#endif
#ifdef WITH_ROOT
    positions.Add(GetIndex(key), GetIndex(AliMpIntPair(i,j))); 
#endif
  } while (!padPos.eof());

  padPos.close();

  if (fVerboseLevel>0) 
    cout << "Opening file "
         << AliMpFiles::Instance()->BergToGCFilePath(fStationType)
         << endl;

  ifstream bergToGCFile(AliMpFiles::Instance()->BergToGCFilePath(fStationType));
  Int_t gassiChannel[80];
  while(1) {
    Int_t bergNum;
    TString gcStr;
    bergToGCFile>>bergNum>>gcStr;
    if (!bergToGCFile.good()) break;
    if (gcStr=="GND") continue;
    if (bergNum>80) {
        Fatal("BuildMotifType","Berg number > 80 ...");
        continue;
    }
    gassiChannel[bergNum-1]= atoi(gcStr);
  }
  bergToGCFile.close();
  
  TString strMotif 
    = AliMpFiles::Instance()
      ->MotifFilePath(fStationType, fPlaneType, motifTypeId);
  ifstream motif(strMotif);
  if (fVerboseLevel>0) cout<<"Opening file "<<strMotif<<endl;


  Int_t nofPadsX=0;
  Int_t nofPadsY=0;

  do {
  
    Int_t ix,iy,numBerg,numKapton,padNum,gassiNum;

    TString lineStr,token;
    lineStr.ReadLine(motif);
    if (!motif.good()) break;
#if defined (__HP_aCC) || (__alpha)
    strstream tokenList;
    tokenList << lineStr.Data();
#else
    istringstream tokenList(lineStr.Data());
#endif    
    
    token.ReadToken(tokenList);
    if (!tokenList.good()) continue; // column is missing...
    if ( (token.Length()>0) && (token[0]=='#') ) continue; // this is a comment line
    
    numBerg = atoi(token.Data());
    if (numBerg==0) {
      Warning("BuildMotifType","Berg number invalid");
      continue;
    }
    
    token.ReadToken(tokenList);
    if (!tokenList.good()) continue; // column is missing...
    numKapton = atoi(token.Data());
    if (numKapton==0) continue;

    
    token.ReadToken(tokenList);
    if (!tokenList.good()) continue; // column is missing...
    if (token=="GND") continue;
    string padName = token.Data();
    padNum = motifType->PadNum(token);
    
     token.ReadToken(tokenList);
     if (token.IsNull() ) continue; // column is missing...
//     if (token[0]!='E') {
//       cerr<<"Problem : gassinumber isn't begining with E:"<<token<<endl;
//       continue;
//     }  else {
//        gassiNum = atoi(token.Data() +1 )-1;
//     }
    if ( (numBerg<1) || (numBerg>80) ) {
        Warning("BuildMotifType","Berg number outside range");
        continue;
    }
    
    gassiNum  = gassiChannel[numBerg-1];

#ifdef WITH_STL
    PadMapTypeIterator iter = positions.find(padName);
    if (iter==positions.end()) {
      cerr<<"Problem: Pad number "<<padNum<<" found in the file "<<strMotif
	  <<" but not in the file"<<strPadPos<<endl;
      continue;
    }

    ix= iter->second.first;
    iy= iter->second.second;
#endif

#ifdef WITH_ROOT
    Long_t value = positions.GetValue(GetIndex(padName));
    if (!value) {
      cerr<<"Problem: Pad number "<<padNum<<" found in the file "<<strMotif
	  <<" but not in the file"<<strPadPos<<endl;
      continue;
    }

    ix = GetPair(value).GetFirst();
    iy = GetPair(value).GetSecond();
#endif

    motifType->AddConnection(AliMpIntPair(ix,iy),
                  new AliMpConnection(padNum,numBerg,numKapton,gassiNum));

    if (ix>=nofPadsX) nofPadsX=ix+1;
    if (iy>=nofPadsY) nofPadsY=iy+1;

  } while (!motif.eof());    


  motifType->SetNofPads(nofPadsX, nofPadsY);

  motif.close();

  return motifType;
}


//_____________________________________________________________________________
AliMpMotifSpecial*  
AliMpReader::BuildMotifSpecial(const TString& motifID,
                               AliMpMotifType* motifType)
{
// Build a special motif by reading the file motifSpecial<motifId>.dat
// in the data directory
// ---

  // Open the input file
  ifstream in(AliMpFiles::Instance()
              ->MotifSpecialFilePath(fStationType, fPlaneType, motifID).Data(), 
	      ios::in);
  if (!in) {	
     Error("BuildMotifSpecial", "File not found.");
     return 0;
  }

  AliMpMotifSpecial* res = new AliMpMotifSpecial(motifID,motifType);
  Int_t i,j;
  Double_t x,y;
  in >> i;
  while (!in.eof()){
    in >>j >>x >> y;
    res->SetPadDimensions(AliMpIntPair(i,j),TVector2(x/2.,y/2.));
    in >> i;
  }
  
  in.close();
  return res;
}


//_____________________________________________________________________________
void AliMpReader::SetVerboseLevel(Int_t verboseLevel)
{
// Sets verbose level.
// ---

  fVerboseLevel = verboseLevel;
}

