// $Id$
// Category: sector
//
// Class AliMpReader
// -------------------
// Class that takes care of reading the sector data.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <string>
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
#include "AliMpRowSegmentSpecial.h"
#include "AliMpPadRow.h"
#include "AliMpPadRowSegment.h"
#include "AliMpMotifMap.h"
#include "AliMpMotif.h"
#include "AliMpMotifSpecial.h"
#include "AliMpMotifType.h"
#include "AliMpConnection.h"
#include "AliMpIntPair.h"
#include "AliMpDirection.h"

ClassImp(AliMpReader)

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
AliMpReader::AliMpReader(AliMpPlaneType plane) 
  : TObject(),
    fPlaneType(plane),
    fSector(0),
    fVerboseLevel(0)
{
//
}

//_____________________________________________________________________________
AliMpReader::AliMpReader() 
  : TObject(),
    fPlaneType(kBendingPlane),
    fSector(0),
    fVerboseLevel(0)
{
//
}

//_____________________________________________________________________________
AliMpReader::~AliMpReader() {
//  
}

//
// private methods
//

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
void AliMpReader::ReadSectorSpecialData(ifstream& in)
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
  ReadRowSpecialData(in); 
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

    AliMpVMotif* motif =  ReadMotifData(in, fSector->GetZone(1));
    AliMpSubZone* subZone = new AliMpSubZone(motif); 
    fSector->GetZone(1)->AddSubZone(subZone); 
             // special row segments are always in zone 1
  
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
void AliMpReader::ReadRowSpecialData(ifstream& in)
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
  AliMpVRowSegment* firstNormalSeg = row->GetRowSegment(0);
  Double_t offsetX = firstNormalSeg->LeftBorderX();
  
  // Create a special row segment
  AliMpRowSegmentSpecial* segment = new AliMpRowSegmentSpecial(row, offsetX);
  row->AddRowSegmentInFront(segment);
  
  TString nextKeyword;
  in >> nextKeyword;
  if (fVerboseLevel>0) 
    cout << nextKeyword << " ";
    
  if (nextKeyword != fgkPadRowsKeyword) {
     Fatal("ReadRowSpecialData", "Wrong file format.");
     return;
  }  
    
  ReadRowSegmentSpecialData(in, segment); 
  
  // Update row segment and set it to all subzones associated with
  // contained motifs
  
  segment->UpdateMotifVector();
  segment->UpdatePadsOffset();
  
  for (Int_t i=0; i<segment->GetNofMotifs(); i++)
    fSector->GetZone(1)->FindSubZone(segment->GetMotif(i))->AddRowSegment(segment);
}  

//_____________________________________________________________________________
void AliMpReader::ReadRowSegmentSpecialData(ifstream& in, 
                                            AliMpRowSegmentSpecial* segment)
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
     AliMpPadRow* padRow = new AliMpPadRow();
     segment->AddPadRow(padRow);
     
     // Keep the new rows in a temporary vector
     newPadRows.push_back(padRow);
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
      AliMpPadRow* padRow = newPadRows[i];
      
      // Find motif
      AliMpVMotif* motif = fSector->GetMotifMap()->FindMotif(motifId);
      
      if (!motif) {
        Fatal("ReadRowSegmentSpecialData", "Unknown motif.");
	return;
      }

      // Create pad row segment
      AliMpPadRowSegment* padRowSegment 
        = new AliMpPadRowSegment(padRow, dynamic_cast<AliMpMotif *>(motif), 
                             motifPositionId, nofPadsInRow);
      padRow->AddPadRowSegment(padRowSegment);	 
    }  
  }
  while (!in.eof() && (nextKeyword == fgkPadRowSegmentKeyword));
  
  if (in.eof()) return;

  if (nextKeyword == fgkPadRowsKeyword) {
    ReadRowSegmentSpecialData(in, segment);
  }
  else if (nextKeyword == fgkRowSpecialKeyword) {
    ReadRowSpecialData(in);
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
  ifstream in(AliMpFiles::Instance()->SectorFilePath(fPlaneType).Data(), ios::in);
  if (!in) {	
     Error("Read", "File not found.");
     return 0;
  }
  
  ReadSectorData(in);
  fSector->SetRowSegmentOffsets();

  // Open second input file
  TString sectorSpecialFileName 
    = AliMpFiles::Instance()->SectorSpecialFilePath(fPlaneType);
  if (!gSystem->AccessPathName(sectorSpecialFileName.Data())) {
    ifstream in2(sectorSpecialFileName.Data(), ios::in);
    if (!in2) {	
       Error("Read", "File not found.");
       return 0;
    }
    
    ReadSectorSpecialData(in2);
  }   

  fSector->Initialize();
  
  return fSector;
}  

//_____________________________________________________________________________
AliMpMotifType* AliMpReader::BuildMotifType(TString motifTypeId)
{

  // Read the files describing a motif in the "$MINSTALL/data" directory
  // and fill the AliMpMotifType structure with.
  // The files mentioned are are named padPos<maskName>.dat
  // and connect<maskName>.dat

  AliMpMotifType*  motifType = new AliMpMotifType(motifTypeId);	

  TString strPadPos 
    = AliMpFiles::Instance()->PadPosFilePath(fPlaneType,motifTypeId);
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
    positions[key].first=i;
    positions[key].second=j;
  } while (!padPos.eof());

  padPos.close();

  if (fVerboseLevel>0) cout <<"Opening file "
                            <<AliMpFiles::Instance()->BergToGCFilePath()
                            <<endl;

  ifstream bergToGCFile(AliMpFiles::Instance()->BergToGCFilePath());
  Int_t GassiChannel[80];
  while(1) {
    Int_t bergNum;
    TString GCStr;
    bergToGCFile>>bergNum>>GCStr;
    if (!bergToGCFile.good()) break;
    if (GCStr=="GND") continue;
    if (bergNum>80) {
        Fatal("BuildMotifType","Berg number > 80 ...");
        continue;
    }
    GassiChannel[bergNum-1]= atoi(GCStr);
  }
  bergToGCFile.close();
  
  TString strMotif 
    = AliMpFiles::Instance()->MotifFilePath(fPlaneType,motifTypeId);
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
    
    gassiNum  = GassiChannel[numBerg-1];

    PadMapTypeIterator iter = positions.find(padName);
    if (iter==positions.end()) {
      cerr<<"Problem: Pad number "<<padNum<<" found in the file "<<strMotif
	  <<" but not in the file"<<strPadPos<<endl;
      continue;
    }

    ix= iter->second.first;
    iy= iter->second.second;

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
AliMpReader::BuildMotifSpecial(TString motifID,AliMpMotifType* motifType)
{
// Build a special motif by reading the file motifSpecial<motifId>.dat
// in the data directory
// ---

  // Open the input file
  ifstream in(AliMpFiles::Instance()
                ->MotifSpecialFilePath(fPlaneType, motifID).Data(), ios::in);
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

