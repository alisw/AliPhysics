/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
// $MpId: AliMpSectorReader.cxx,v 1.9 2006/05/24 13:58:46 ivana Exp $
// Category: sector
//
// Class AliMpSectorReader
// -----------------------
// Class that takes care of reading the sector data.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpSectorReader.h"
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
#include "AliMpMotifReader.h"
#include "AliMpMotifMap.h"
#include "AliMpMotif.h"
#include "AliMpMotifSpecial.h"
#include "AliMpMotifType.h"
#include "AliMpConnection.h"
#include "AliMpIntPair.h"
#include "AliMpDirection.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include <Riostream.h>
#include <Rstrstream.h>
#include <TSystem.h>
#include <TMath.h>

#if !defined(__HP_aCC) && !defined(__alpha)
  #include <sstream>
#endif

/// \cond CLASSIMP
ClassImp(AliMpSectorReader)
/// \endcond

const TString  AliMpSectorReader::fgkSectorKeyword  = "SECTOR_DATA";
const TString  AliMpSectorReader::fgkZoneKeyword    = "ZONE";
const TString  AliMpSectorReader::fgkSubZoneKeyword = "SUBZONE";
const TString  AliMpSectorReader::fgkRowKeyword     = "ROW_SEGMENT";
const TString  AliMpSectorReader::fgkEofKeyword     = "EOF";
const TString  AliMpSectorReader::fgkSectorSpecialKeyword  = "SECTOR_SPECIAL_DATA";
const TString  AliMpSectorReader::fgkMotifKeyword          = "MOTIF";
const TString  AliMpSectorReader::fgkRowSpecialKeyword     = "ROW";
const TString  AliMpSectorReader::fgkPadRowsKeyword        = "PAD_ROWS";
const TString  AliMpSectorReader::fgkPadRowSegmentKeyword  = "PAD_ROW_SEGMENT";

//_____________________________________________________________________________
AliMpSectorReader::AliMpSectorReader(AliMpStationType station, 
                                     AliMpPlaneType plane) 
  : TObject(),
    fStationType(station),
    fPlaneType(plane),
    fSector(0),
    fMotifReader(new AliMpMotifReader(station, plane))
{
// Standard constructor
}

//_____________________________________________________________________________
AliMpSectorReader::AliMpSectorReader() 
  : TObject(),
    fStationType(kStation1),
    fPlaneType(kBendingPlane),
    fSector(0),
    fMotifReader(0)
{
// Default constructor
}

//_____________________________________________________________________________
AliMpSectorReader::~AliMpSectorReader() 
{
/// Destructor  

  delete fMotifReader;
}

//
// private methods
//

//_____________________________________________________________________________
void  AliMpSectorReader::ReadSectorData(ifstream& in)
{
/// Read sector input data;
/// prepare zones and rows vectors to be filled in.

  TString keyword;
  in >> keyword;
  
  AliDebugStream(2) << keyword << endl;

  if (keyword != fgkSectorKeyword) {
     Fatal("ReadSectorData", "Wrong file format.");
     return;
  }   
    
  Int_t nofZones, nofRows;
  TString directionStr;
  Double_t offsetX, offsetY;
  in >> nofZones;
  in >> nofRows;
  in >> directionStr;
  in >> offsetX;
  in >> offsetY;
  
  AliMpDirection direction;
  direction = (directionStr == "Y") ? kY  :  kX;

  AliDebugStream(2) << nofZones << " " <<  nofRows << endl;

  fSector = new AliMpSector("Not defined", nofZones, nofRows,direction,
                            TVector2(offsetX, offsetY));
  
  TString nextKeyword;
  in >> nextKeyword;
    
  if (nextKeyword != fgkZoneKeyword) {
     Fatal("ReadSectorData", "Wrong file format.");
     return;
  }      
  
  ReadZoneData(in);
}  

//_____________________________________________________________________________
void AliMpSectorReader::ReadZoneData(ifstream& in)
{
/// Read zone input data;
/// create zone and adds it to zones vector.

  Int_t zoneID;
  Double_t  sizex, sizey; 
  in >> zoneID;    
  in >> sizex;
  in >> sizey;
  AliDebugStream(2)
     << fgkZoneKeyword << " " <<  zoneID << "  " 
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
void AliMpSectorReader::ReadSubZoneData(ifstream& in, AliMpZone* zone)
{
/// Read subzone input data;
/// create subzone and its to the specified zone.

  AliDebugStream(2) << fgkSubZoneKeyword << endl;

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
AliMpVMotif*  AliMpSectorReader::ReadMotifData(ifstream& in, AliMpZone* zone)
{
/// Read the motif input data.

  TString  motifID;
  TString  motifTypeID;
  in >> motifID;
  in >> motifTypeID;

  AliDebugStream(2) << motifID << " " << motifTypeID << endl;
  
  AliMpMotifMap* motifMap = fSector->GetMotifMap();

  AliMpMotifType* motifType = 0;
  AliMpVMotif* motif 
    = motifMap->FindMotif(motifID, motifTypeID, zone->GetPadDimensions());
  if (!motif) {    
    motifType = motifMap->FindMotifType(motifTypeID);
    if (!motifType) {
      motifType = fMotifReader->BuildMotifType(motifTypeID);     
      motifMap->AddMotifType(motifType);
    }
    
    if (zone->GetPadDimensions().X() != 0. && zone->GetPadDimensions().Y() != 0.) 
      motif = new AliMpMotif(motifID, motifType, zone->GetPadDimensions());
    else 
      motif = fMotifReader->BuildMotifSpecial(motifID, motifType);
      
    if (motif) 
      motifMap->AddMotif(motif);

  }  

  return motif;   
}  

//_____________________________________________________________________________
void AliMpSectorReader::ReadRowSegmentsData(ifstream& in, 
                                      AliMpZone* zone, AliMpSubZone* subZone)
{
/// Read row segments input data of a specified zone and subzone;
/// creates row segment and add it to the specified subzone
/// and a corresponding row in the rows vector.

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
    
    firstMotifPositionId |= AliMpConstants::ManuMask(fPlaneType);
    
    AliDebugStream(2)
      << fgkRowKeyword << " " 
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
void AliMpSectorReader::ReadSectorSpecialData(ifstream& in, AliMpXDirection direction)
{
/// Read sector input data
/// with a special (irregular) motifs.

  TString keyword;
  in >> keyword;

  AliDebugStream(2) << keyword << endl;

  if (keyword != fgkSectorSpecialKeyword) {
     Fatal("ReadSectorSpecialData", "Wrong file format.");
     return;
  }   

  TString nextKeyword;
  in >> nextKeyword;

  AliDebugStream(2) << keyword << endl;
    
  if (nextKeyword != fgkMotifKeyword) {
    Fatal("ReadSectorSpecialData", "Wrong file format.");
    return;
  }  

  ReadMotifsSpecialData(in);     
  ReadRowSpecialData(in, direction); 
}  

//_____________________________________________________________________________
void AliMpSectorReader::ReadMotifsSpecialData(ifstream& in)
{
/// Read the special (irregular) motifs input data.

  AliDebugStream(2) << fgkMotifKeyword << endl;

  TString nextKeyword;
  do {
    Int_t zone;
    in >> zone;
    AliMpVMotif* motif =  ReadMotifData(in, fSector->GetZone(zone));
    AliMpSubZone* subZone = new AliMpSubZone(motif); 
    fSector->GetZone(zone)->AddSubZone(subZone); 
  
    in >> nextKeyword;

    AliDebugStream(2) << nextKeyword << endl;      
  }
  while (nextKeyword == fgkMotifKeyword);
    
  if (nextKeyword != fgkRowSpecialKeyword) {
     Fatal("ReadMotifSpecialData", "Wrong file format.");
     return;
  }      
}  

//_____________________________________________________________________________
void AliMpSectorReader::ReadRowSpecialData(ifstream& in, AliMpXDirection direction)
{
/// Read row input data
/// with a special (irregular) motifs.

  Int_t id;
  in >> id;

  AliDebugStream(2) << id << endl;      
  
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
  
  AliDebugStream(2) << nextKeyword << endl;
    
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
void AliMpSectorReader::ReadRowSegmentSpecialData(ifstream& in, 
                                            AliMpVRowSegmentSpecial* segment,
					    AliMpXDirection direction)
{
/// Read row segment input data
/// with a special (irregular) motifs.

  Int_t nofPadRows;
  in >> nofPadRows;
  
  AliDebugStream(2) << nofPadRows << endl;
  
  TString keyword;
  in >> keyword;

  AliDebugStream(2) << keyword << endl;
    
  if (keyword != fgkPadRowSegmentKeyword) {
     Fatal("ReadRowSegmentSpecialData", "Wrong file format.");
     return;
  }  
  
  //
  // Process data
  //
    
  AliMpVRowSegmentSpecial::PadRowVector  newPadRows;
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
  
    motifPositionId |= AliMpConstants::ManuMask(fPlaneType);

    AliDebugStream(2)
      << nofPadsInRow << " " << motifId << " " << motifPositionId << endl;

    in >> nextKeyword;

    AliDebugStream(2) << nextKeyword << endl;

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
AliMpSector* AliMpSectorReader::BuildSector()
{
/// Read the mapping data from ascii data file
/// and create the basic objects:                                            \n
/// zones, subzones, rows, row segments, motifs.

  // Open input file
  ifstream in(AliMpFiles::SectorFilePath(fStationType, fPlaneType).Data(), ios::in);
  if (!in) {
     AliErrorStream()
       << "File " << AliMpFiles::SectorFilePath(fStationType, fPlaneType) 
       << " not found." << endl;
     return 0;
  }
  
  ReadSectorData(in);
  fSector->SetRowSegmentOffsets();

  // Open input file for special inner zone
  TString sectorSpecialFileName 
    = AliMpFiles::SectorSpecialFilePath(fStationType, fPlaneType);
  if (!gSystem->AccessPathName(sectorSpecialFileName.Data())) {
    ifstream in2(sectorSpecialFileName.Data(), ios::in);
    if (!in2) {	
       AliErrorStream()
         << "File " << AliMpFiles::SectorSpecialFilePath(fStationType, fPlaneType) 
	 << " not found." << endl;
       return 0;
    }
    
    ReadSectorSpecialData(in2, kLeft);
  }   

  // Open input file for special outer zone
  TString sectorSpecialFileName2 
    = AliMpFiles::SectorSpecialFilePath2(fStationType, fPlaneType);
  if (!gSystem->AccessPathName(sectorSpecialFileName2.Data())) {
    ifstream in3(sectorSpecialFileName2.Data(), ios::in);
    if (!in3) {	
       AliErrorStream()
         << "File " << AliMpFiles::SectorSpecialFilePath2(fStationType, fPlaneType) 
	 << " not found."<< endl;	
       return 0;
    }
    
    ReadSectorSpecialData(in3, kRight);
  }   

  fSector->Initialize();
  
  return fSector;
}  

