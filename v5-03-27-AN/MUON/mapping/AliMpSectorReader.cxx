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

//-----------------------------------------------------------------------------
// Class AliMpSectorReader
// -----------------------
// Class that takes care of reading the sector data.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpSectorReader.h"
#include "AliMpSector.h"
#include "AliMpFiles.h"
#include "AliMpDataStreams.h"
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
#include "AliMpDirection.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include <Riostream.h>
#include <Rstrstream.h>
#include <TSystem.h>
#include <TMath.h>

#include <limits>

#if !defined(__HP_aCC) && !defined(__alpha)
  #include <sstream>
#endif

/// \cond CLASSIMP
ClassImp(AliMpSectorReader)
/// \endcond

//
// static private methods
//

//_____________________________________________________________________________
const TString& AliMpSectorReader::GetSectorKeyword()
{
  /// sector keyword
  static const TString kSectorKeyword  = "SECTOR_DATA";
  return kSectorKeyword;
}  

//_____________________________________________________________________________
const TString& AliMpSectorReader::GetZoneKeyword()
{
  /// zone keyword
  static const TString kZoneKeyword = "ZONE";
  return kZoneKeyword;
}  

//_____________________________________________________________________________
const TString& AliMpSectorReader::GetSubZoneKeyword() 
{
  /// subzone keyword
  static const TString kSubZoneKeyword = "SUBZONE";
  return kSubZoneKeyword;
}  

//_____________________________________________________________________________
const TString& AliMpSectorReader::GetRowKeyword()
{
  /// row keyword
  static const TString kRowKeyword = "ROW_SEGMENT";
  return kRowKeyword;
}  

//_____________________________________________________________________________
const TString& AliMpSectorReader::GetSectorSpecialKeyword()
{
  /// sector special keyword
  static const TString kSectorSpecialKeyword = "SECTOR_SPECIAL_DATA";
  return kSectorSpecialKeyword;
}  

//_____________________________________________________________________________
const TString& AliMpSectorReader::GetMotifKeyword()
{
  /// motif keyword
  static const TString kMotifKeyword = "MOTIF";
  return kMotifKeyword;
}  

//_____________________________________________________________________________
const TString& AliMpSectorReader::GetRowSpecialKeyword()
{
  /// row special keyword
  static const TString kRowSpecialKeyword = "ROW";
  return kRowSpecialKeyword;
}  

//_____________________________________________________________________________
const TString& AliMpSectorReader::GetPadRowsKeyword()
{
  /// pad rows keyword
  static const TString kPadRowsKeyword = "PAD_ROWS";
  return kPadRowsKeyword;
}  

//_____________________________________________________________________________
const TString& AliMpSectorReader::GetPadRowSegmentKeyword()
{
  /// pad row segment keyword
  static const TString kPadRowSegmentKeyword = "PAD_ROW_SEGMENT";
  return kPadRowSegmentKeyword;
}  

//
// ctors, dtor
//
  
//_____________________________________________________________________________
AliMpSectorReader::AliMpSectorReader(const AliMpDataStreams& dataStreams,
                                     AliMq::Station12Type station, 
                                     AliMp::PlaneType plane) 
  : TObject(),
    fkDataStreams(dataStreams),
    fStationType(station),
    fPlaneType(plane),
    fSector(0),
    fMotifReader(new AliMpMotifReader(dataStreams, AliMp::kStation12, station, plane))
 
{
/// Standard constructor
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
void  AliMpSectorReader::ReadSectorData(istream& in)
{
/// Read sector input data;
/// prepare zones and rows vectors to be filled in.

  TString keyword;
  in >> keyword;
  
  AliDebugStream(2) << keyword << endl;

  if (keyword != GetSectorKeyword()) {
     AliErrorStream() << "Wrong file format." << endl;
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
  
  AliMp::Direction direction;
  direction = (directionStr == "Y") ? AliMp::kY  :  AliMp::kX;

  AliDebugStream(2) << nofZones << " " <<  nofRows << endl;
  
  if ( nofZones < 0 || nofZones >= std::numeric_limits<Int_t>::max() ||
       nofRows < 0  || nofRows >= std::numeric_limits<Int_t>::max() ) {
    AliErrorStream() << "Wrong nofZones/nofRows value." << endl;
    return;
  }         

  fSector = new AliMpSector("Not defined", nofZones, nofRows,direction,
                            offsetX, offsetY);
  
  TString nextKeyword;
  in >> nextKeyword;
    
  if (nextKeyword != GetZoneKeyword()) {
     AliErrorStream() << "Wrong file format." << endl;
     return;
  }      
  
  ReadZoneData(in);
}  

//_____________________________________________________________________________
void AliMpSectorReader::ReadZoneData(istream& in)
{
/// Read zone input data;
/// create zone and adds it to zones vector.

  Int_t zoneID;
  Double_t  sizex, sizey; 
  in >> zoneID;    
  in >> sizex;
  in >> sizey;
  AliDebugStream(2)
     << GetZoneKeyword() << " " <<  zoneID << "  " 
     << sizex << " " << sizey << endl;
  
  AliMpZone* zone =  fSector->GetZone(zoneID);
  zone->SetPadDimensions(sizex/2.,sizey/2.); 
      
  TString nextKeyword;
  in >> nextKeyword;
    
  if (nextKeyword != GetSubZoneKeyword()) {
     AliErrorStream() << "Wrong file format." << endl;
     return;
  }  
    
  ReadSubZoneData(in, zone);
}

//_____________________________________________________________________________
void AliMpSectorReader::ReadSubZoneData(istream& in, AliMpZone* zone)
{
/// Read subzone input data;
/// create subzone and its to the specified zone.

  AliDebugStream(2) << GetSubZoneKeyword() << endl;

  AliMpVMotif* motif = ReadMotifData(in, zone);
  AliMpSubZone* subZone = new AliMpSubZone(motif); 
  zone->AddSubZone(subZone); 

  TString nextKeyword;
  in >> nextKeyword;
    
  if (nextKeyword != GetRowKeyword()) {
     AliErrorStream() << "Wrong file format." << endl;
     return;
  }  
    
  ReadRowSegmentsData(in, zone, subZone);
}   

//_____________________________________________________________________________
AliMpVMotif*  AliMpSectorReader::ReadMotifData(istream& in, AliMpZone* zone)
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
    = motifMap->FindMotif(motifID, motifTypeID, 
                          zone->GetPadDimensionX(), zone->GetPadDimensionY());
  if (!motif) {    
    motifType = motifMap->FindMotifType(motifTypeID);
    if (!motifType) {
      motifType = fMotifReader->BuildMotifType(motifTypeID);     
      motifMap->AddMotifType(motifType);
    }
    
    if (zone->GetPadDimensionX() != 0. && zone->GetPadDimensionY() != 0.) 
      motif = new AliMpMotif(motifID, motifType, 
                     zone->GetPadDimensionX(), zone->GetPadDimensionY());
    else 
      motif = fMotifReader->BuildMotifSpecial(motifID, motifType);
      
    if (motif) 
      motifMap->AddMotif(motif);

  }  

  return motif;   
}  

//_____________________________________________________________________________
void AliMpSectorReader::ReadRowSegmentsData(istream& in, 
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
      << GetRowKeyword() << " " 
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
      = new AliMpRowSegment(row, motif, offX, offY, nofMotifs, 
                        firstMotifPositionId, firstMotifPositionDId);
			
    subZone->AddRowSegment(rowSegment);
    row->AddRowSegment(rowSegment);
  }
  while (!in.eof() && (nextKeyword == GetRowKeyword()));

  if (in.eof()) return;

  if (nextKeyword == GetZoneKeyword()) {
    ReadZoneData(in);
  }
  else if (nextKeyword == GetSubZoneKeyword()) {
    ReadSubZoneData(in, zone);
  }   
  else {
    AliErrorStream() << "Wrong file format." << endl;
  } 
}   

//_____________________________________________________________________________
void AliMpSectorReader::ReadSectorSpecialData(istream& in, AliMp::XDirection direction)
{
/// Read sector input data
/// with a special (irregular) motifs.

  TString keyword;
  in >> keyword;

  AliDebugStream(2) << keyword << endl;

  if (keyword != GetSectorSpecialKeyword()) {
     AliErrorStream() << "Wrong file format." << endl;
     return;
  }   

  TString nextKeyword;
  in >> nextKeyword;

  AliDebugStream(2) << keyword << endl;
    
  if (nextKeyword != GetMotifKeyword()) {
    AliErrorStream() << "Wrong file format." << endl;
    return;
  }  

  ReadMotifsSpecialData(in);     
  ReadRowSpecialData(in, direction); 
}  

//_____________________________________________________________________________
void AliMpSectorReader::ReadMotifsSpecialData(istream& in)
{
/// Read the special (irregular) motifs input data.

  AliDebugStream(2) << GetMotifKeyword() << endl;

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
  while (nextKeyword == GetMotifKeyword());
    
  if (nextKeyword != GetRowSpecialKeyword()) {
     AliErrorStream() << "Wrong file format." << endl;
     return;
  }      
}  

//_____________________________________________________________________________
void AliMpSectorReader::ReadRowSpecialData(istream& in, AliMp::XDirection direction)
{
/// Read row input data
/// with a special (irregular) motifs.

  Int_t id;
  in >> id;

  AliDebugStream(2) << id << endl;      
  
  // Get the row and its border
  AliMpRow* row = fSector->GetRow(id);

  AliMpVRowSegmentSpecial* segment = 0;
  if (direction == AliMp::kLeft) {
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
    
  if (nextKeyword != GetPadRowsKeyword()) {
     AliErrorStream() << "Wrong file format." << endl;
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
    
    if (subZone) subZone->AddRowSegment(segment);
  }  
}  

//_____________________________________________________________________________
void AliMpSectorReader::ReadRowSegmentSpecialData(istream& in, 
                                            AliMpVRowSegmentSpecial* segment,
					    AliMp::XDirection direction)
{
/// Read row segment input data
/// with a special (irregular) motifs.

  Int_t nofPadRows;
  in >> nofPadRows;
  
  AliDebugStream(2) << nofPadRows << endl;
  
  if ( nofPadRows < 0 || nofPadRows >= std::numeric_limits<Int_t>::max()) {
    AliErrorStream() << "Wrong nofPadRows value." << endl;
    return;
  }         

  TString keyword;
  in >> keyword;

  AliDebugStream(2) << keyword << endl;
    
  if (keyword != GetPadRowSegmentKeyword()) {
     AliErrorStream() << "Wrong file format." << endl;
     return;
  }  
  
  //
  // Process data
  //
    
  TObjArray newPadRows;
  for (Int_t i=0; i<nofPadRows; i++) {
    
     // Create pad row
     AliMpPadRow* padRow = new AliMpPadRow(direction);
     segment->AddPadRow(padRow);
     
     // Keep the new rows in a temporary vector
     newPadRows.Add(padRow);
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
      AliMpPadRow* padRow = (AliMpPadRow*)newPadRows[i];
      
      // Find motif
      AliMpVMotif* motif = fSector->GetMotifMap()->FindMotif(motifId);
      
      if (!motif) {
        AliErrorStream() << "Unknown motif" << endl;
	return;
      }

      // Create pad row segment
      padRow->AddPadRowSegment(dynamic_cast<AliMpMotif *>(motif), 
                               motifPositionId, nofPadsInRow);
    }  
  }
  while (!in.eof() && (nextKeyword == GetPadRowSegmentKeyword()));
  
  if (in.eof()) return;

  if (nextKeyword == GetPadRowsKeyword()) {
    ReadRowSegmentSpecialData(in, segment, direction);
  }
  else if (nextKeyword == GetRowSpecialKeyword()) {
    ReadRowSpecialData(in, direction);
  }   
  else {
    AliErrorStream() << "Wrong file format." << endl;
  } 
}  

//
// public methods
//

//_____________________________________________________________________________
AliMpSector* AliMpSectorReader::BuildSector()
{
/// Read the mapping data from stream and create the basic objects:         \n
/// zones, subzones, rows, row segments, motifs.

  // Open input stream
  //
  istream& in 
    = fkDataStreams.
        CreateDataStream(AliMpFiles::SectorFilePath(fStationType,fPlaneType));

  ReadSectorData(in);
  delete &in;
  
  fSector->SetRowSegmentOffsets();

  // Open input stream for special inner zone
  
  // add is data function
  
  TString sectorSpecialFileName 
    = AliMpFiles::SectorSpecialFilePath(fStationType, fPlaneType);
  if ( fkDataStreams.IsDataStream(sectorSpecialFileName) ) {
    istream& in2 
      = fkDataStreams.
          CreateDataStream(sectorSpecialFileName);
  
    ReadSectorSpecialData(in2, AliMp::kLeft);
    
    delete &in2;
  }  

  // Open input file for special outer zone
  TString sectorSpecialFileName2 
    = AliMpFiles::SectorSpecialFilePath2(fStationType, fPlaneType);
  if ( fkDataStreams.IsDataStream(sectorSpecialFileName2) ) {
    istream& in3
      = fkDataStreams.
          CreateDataStream(sectorSpecialFileName2);
    
    ReadSectorSpecialData(in3, AliMp::kRight);
    
    delete &in3;
  }   

  fSector->Initialize();
  
  return fSector;
}  

