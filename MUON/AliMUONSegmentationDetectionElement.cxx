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

/* $Id$ */
//===========================================================
//  Segmentation classes for MUON Detection Elements      
//        Gines MARTINEZ, SUBATECH July 04                
//  This class interfaces with the mapping and segmentation
//  files MUON.
//  This files are placed by default in 
//  $ALICE_ROOT/MUON/mapping/data/Stationxxx/yyy_plane/
//  There are in tracking 23 types of detection elements
//  8 SectorSt1, 8 SectorSt2, 2 122000SR1, 2 122000NR1, 4 112200SR2, 4 112200NR2 
//  4 122200S, 4 122200N, 8 222000N,8 220000N,  8 330000N, 4 122300N, 8 112230NR3 
//  8 112230N, 8 222330N, 8 223300N, 16 333000N, 4  122330N, 8 112233NR3, 8 112233N 
//  8 222333N, 8 223330N, 8 333300N 
//  Detailed information in Alice Technical Note xxxxxxxx (2004)
//===========================================================

#include <Riostream.h>

#include "TClonesArray.h"
#include "TMap.h"

#include "AliMUONSegmentationDetectionElement.h"
#include "AliMUONSegmentManuIndex.h"
#include "AliMUONSegmentPosition.h"
#include "AliMUONSegmentIndex.h"


//___________________________________________
ClassImp(AliMUONSegmentationDetectionElement)

//static data member
const TString AliMUONSegmentationDetectionElement::fgkDefaultTop =getenv("ALICE_ROOT") ;
const TString AliMUONSegmentationDetectionElement::fgkStationDir = "/MUON/mapping/data/station";
const TString AliMUONSegmentationDetectionElement::fgkBendingDir= "/bending_plane";  
const TString AliMUONSegmentationDetectionElement::fgkNonBendingDir= "/non-bending_plane";
const TString AliMUONSegmentationDetectionElement::fgkFileExt = ".map";
const TString AliMUONSegmentationDetectionElement::fgkBendingExt = ".Bending";
const TString AliMUONSegmentationDetectionElement::fgkNonBendingExt = ".NonBending";

//___________________________________________

AliMUONSegmentationDetectionElement::AliMUONSegmentationDetectionElement() : TObject()
{
  fMapManuIndexIndex= 0x0;
  fMapIndexManuIndex= 0x0;
  fMapIndexPosition= 0x0;
  fXlocalSegmentPositions= 0x0;
  fYlocalSegmentPositions= 0x0;
}


// AliMUONSegmentationDetectionElement::AliMUONSegmentationDetectionElement(const char* ElementType)
// {
  
//   fMapManuIndexIndex= 0x0;
//   fMapIndexManuIndex= 0x0;
//   fMapIndexPosition= 0x0;
//   fXlocalSegmentPositions= 0x0;
//   fYlocalSegmentPositions= 0x0;
// }

//________________________________________________
AliMUONSegmentationDetectionElement::AliMUONSegmentationDetectionElement(const AliMUONSegmentationDetectionElement& rhs): TObject(rhs) 
{
// Protected copy constructor

  Fatal("AliMUONSegmentationDetectionElementModule", "Not implemented.");
}
//_________________________________________________
AliMUONSegmentationDetectionElement::~AliMUONSegmentationDetectionElement(){
  fListOfIndexes->Delete();
  fListOfManuIndexes->Delete();
  fListOfPositions->Delete();
  fMapManuIndexIndex->Clear();
  fMapIndexManuIndex->Clear();
  fMapIndexPosition->Clear();
}
//_________________________________________________
AliMUONSegmentIndex * AliMUONSegmentationDetectionElement::GetIndex( const char * SegmentManuIndexName)
{
  if (fMapManuIndexIndex) return  (AliMUONSegmentIndex*)  fMapManuIndexIndex->GetValue(SegmentManuIndexName);
  else {
    Warning("GetIndex","SegmentManuIndex %s out of DetectionElement Mapping %s",
	    SegmentManuIndexName,fDetectionElementType.Data());
    return 0x0;
  }
}

//_________________________________________________
AliMUONSegmentManuIndex * AliMUONSegmentationDetectionElement::GetManuIndex( const char * SegmentIndexName)
{
  if (fMapIndexManuIndex) return (AliMUONSegmentManuIndex*) fMapIndexManuIndex->GetValue(SegmentIndexName);
  else {
    Warning("GetManuIndex","SegmentIndex %s out of Detection Element mapping %s",
	    SegmentIndexName,fDetectionElementType.Data());
    return 0x0;
  }
}

//_________________________________________________
AliMUONSegmentPosition  * AliMUONSegmentationDetectionElement::GetPosition( const char * SegmentIndexName)
{
  if (fMapIndexPosition) return (AliMUONSegmentPosition*) fMapIndexPosition->GetValue(SegmentIndexName);
  else {
    Warning("GetPosition","SegmentIndex %s out of DetectionElement mapping %s",
	    SegmentIndexName, fDetectionElementType.Data());
    return 0x0;
  }
}

//_________________________________________________
AliMUONSegmentManuIndex * AliMUONSegmentationDetectionElement::FindManuIndex(const char* ManuIndexName)
{
  if (fMapManuIndexIndex) return (AliMUONSegmentManuIndex*) fMapManuIndexIndex->FindObject(ManuIndexName);
  else  {
    Warning("FindManuIndex","SegmentManuIndex %s out of DetectionElement mapping %s",
	    ManuIndexName,fDetectionElementType.Data());
    return 0x0;
  }
}

//_________________________________________________
AliMUONSegmentManuIndex * AliMUONSegmentationDetectionElement::FindIndex(const char* IndexName)
{
  if (fMapIndexPosition) return (AliMUONSegmentManuIndex *) fMapIndexPosition->FindObject(IndexName);
  else {
    Warning("FindIndex","SegmentIndex %s out of DetectionElement mapping %s",
	    IndexName,fDetectionElementType.Data());
    return 0x0;
  }
}

//_________________________________________________
void    AliMUONSegmentationDetectionElement::Init(const char * DetectionElementType)
{
  TString ElementType(DetectionElementType);
  fSegmentationMappingFile_Bending = fgkDefaultTop+fgkStationDir+"345"
    +fgkBendingDir+"/"+ElementType+fgkBendingExt+fgkFileExt;
  printf("file is %s\n", fSegmentationMappingFile_Bending.Data());
  fSegmentationMappingFile_NonBending = fgkDefaultTop+fgkStationDir+"345"
    +fgkNonBendingDir+"/"+ElementType+fgkNonBendingExt+fgkFileExt;

  if (fMapManuIndexIndex==0x0) { 
    fListOfIndexes = new TObjArray(10000);
    fListOfManuIndexes =new TObjArray(10000);
    fListOfPositions =new TObjArray(10000);
    fMapManuIndexIndex= new TMap();
    fMapIndexManuIndex = new TMap();
    fMapIndexPosition = new TMap();
  }
  else {
    fListOfIndexes->Delete();
    fListOfManuIndexes->Delete();
    fListOfPositions->Delete();
    fMapManuIndexIndex->Clear();
    fMapIndexManuIndex->Clear();
    fMapIndexPosition->Clear();
  }



  Int_t icathode;
  //Bendingplane
  icathode=0;
  ReadingSegmentationMappingFile(fSegmentationMappingFile_Bending ,icathode);
  //NonBendingplane
  icathode=1;
  ReadingSegmentationMappingFile(fSegmentationMappingFile_NonBending,icathode);
  
}
//_______________________________________________________________
void AliMUONSegmentationDetectionElement::ReadingSegmentationMappingFile(TString infile, Int_t cathode)
{ 
  ifstream in( infile,  ios::in);
  if (!in) {
    //    Error("ReadingSegmentationMappingFile", "File not found.");
  }
  else {
    Int_t id, ix, iy, idmanu, idchannel;
    Float_t x, y;
    while ( !in.eof()) {
      in >> id >> ix >> iy >> x >> y >> idmanu >> idchannel;
      char name[10];
      sprintf(name,"%d%d",ix,iy);
      printf("%s id=%d ix=%d iy=%d x=%f y=%f idmanu=%d and idchannel=%d\n",name,id, ix, iy,  x, y,idmanu, idchannel);
      
      fListOfIndexes->AddAt(     new AliMUONSegmentIndex(id,ix,iy,cathode), id);
      fListOfManuIndexes->AddAt( new AliMUONSegmentManuIndex(id,idmanu,0,idchannel), id);;
      fListOfPositions->AddAt(   new AliMUONSegmentPosition(id, x, y,cathode), id);;

      ( (AliMUONSegmentIndex* )fListOfIndexes->At(id))->Print();
      ( (AliMUONSegmentManuIndex*)fListOfManuIndexes->At(id))->Print();
      ( (AliMUONSegmentPosition*)fListOfPositions->At(id))->Print();
     
      fMapManuIndexIndex->Add( fListOfManuIndexes->At(id), fListOfIndexes->At(id));
      fMapIndexManuIndex->Add( fListOfIndexes->At(id),     fListOfManuIndexes->At(id));
      fMapIndexPosition->Add(  fListOfIndexes->At(id),     fListOfPositions->At(id));
    } 
  }
  in.close();
}
