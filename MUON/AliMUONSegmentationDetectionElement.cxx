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
  //Default constructor
  fMapManuIndexIndex= 0x0;
  fMapIndexManuIndex= 0x0;
  fMapIndexPosition= 0x0;
  fMapPositionIndex=0X0;
}


// AliMUONSegmentationDetectionElement::AliMUONSegmentationDetectionElement(const char* ElementType)
// {
  
//   fMapManuIndexIndex= 0x0;
//   fMapIndexManuIndex= 0x0;
//   fMapIndexPosition= 0x0;
// }

//________________________________________________
AliMUONSegmentationDetectionElement::AliMUONSegmentationDetectionElement(const AliMUONSegmentationDetectionElement& rhs): TObject(rhs) 
{
// Protected copy constructor

  Fatal("AliMUONSegmentationDetectionElementModule", "Not implemented.");
}
//_________________________________________________
AliMUONSegmentationDetectionElement::~AliMUONSegmentationDetectionElement(){
  //Class destructor
  fListOfIndexes->Delete();
  fListOfManuIndexes->Delete();
  fListOfPositions->Delete();
  fMapManuIndexIndex->Clear();
  fMapIndexManuIndex->Clear();
  fMapIndexPosition->Clear();
  fMapPositionIndex->Clear();
}

//_________________________________________________
AliMUONSegmentIndex * AliMUONSegmentationDetectionElement::GetIndex(Int_t manu, Int_t channel) const
{
  // Getting AliMUONSegmentIndex from ManuIndex
  return GetIndex( AliMUONSegmentManuIndex::Name(manu, channel).Data() ) ;
}
//_________________________________________________
AliMUONSegmentIndex * AliMUONSegmentationDetectionElement::GetIndex( const char * SegmentManuIndexName) const
{
  // Getting AliMUONSegmentIndex from name of AliMUONSegmentManuIndex
  if (fMapManuIndexIndex) return  (AliMUONSegmentIndex*)  fMapManuIndexIndex->GetValue(SegmentManuIndexName);
  else {
    Warning("GetIndex","SegmentManuIndex %s out of DetectionElement Mapping %s",
	    SegmentManuIndexName,fDetectionElementType.Data());
    return 0x0;
  }
}
//_________________________________________________
AliMUONSegmentManuIndex * AliMUONSegmentationDetectionElement::GetManuIndex(Int_t padx, Int_t pady, Int_t cathode ) const
{
  // Getting ManuIndex from Index
  return GetManuIndex( AliMUONSegmentIndex::Name(padx, pady, cathode).Data() ); 
}
//_________________________________________________
AliMUONSegmentManuIndex * AliMUONSegmentationDetectionElement::GetManuIndex( const char * SegmentIndexName) const
{
  // Getting ManuIndex from manuname
  if (fMapIndexManuIndex) return (AliMUONSegmentManuIndex*) fMapIndexManuIndex->GetValue(SegmentIndexName);
  else {
    Warning("GetManuIndex","SegmentIndex %s out of Detection Element mapping %s",
	    SegmentIndexName,fDetectionElementType.Data());
    return 0x0;
  }
}
//_________________________________________________
AliMUONSegmentPosition  * AliMUONSegmentationDetectionElement::GetPosition(Int_t padx, Int_t pady, Int_t cathode ) const
{
  //Getting position from index
  return GetPosition(  AliMUONSegmentIndex::Name(padx, pady, cathode).Data() ); 
}
//_________________________________________________
AliMUONSegmentPosition  * AliMUONSegmentationDetectionElement::GetPosition( const char * SegmentIndexName) const
{
  // Getting position from indexname
  if (fMapIndexPosition) return (AliMUONSegmentPosition*) fMapIndexPosition->GetValue(SegmentIndexName);
  else {
    Warning("GetPosition","SegmentIndex %s out of DetectionElement mapping %s",
	    SegmentIndexName, fDetectionElementType.Data());
    return 0x0;
  }
}
//_________________________________________________
AliMUONSegmentIndex * AliMUONSegmentationDetectionElement::GetIndexFromPosition(Float_t x, Float_t y, Int_t cathode) const
{
  // Getting Index from position if position is a center pad position
  return GetIndexFromPosition( AliMUONSegmentPosition::Name(x,y, cathode).Data() );
}
//_________________________________________________
AliMUONSegmentIndex * AliMUONSegmentationDetectionElement::GetIndexFromPosition( const char * PositionName) const
{
  // Getting index form positionname
  if (fMapPositionIndex) return  (AliMUONSegmentIndex*)  fMapPositionIndex->GetValue(PositionName);
  else {
    Warning("GetIndexFromPosition","SegmentPosition %s out of DetectionElement Mapping %s",
	    PositionName,fDetectionElementType.Data());
    return 0x0;
  }
}

//_________________________________________________
AliMUONSegmentManuIndex * AliMUONSegmentationDetectionElement::FindManuIndex( const char* ManuIndexName)const
{
  // Getting AliMUONSegmentManuIndex objecto from manu index
  if (fMapManuIndexIndex) return (AliMUONSegmentManuIndex*) fMapManuIndexIndex->FindObject(ManuIndexName);
  else  {
    Warning("FindManuIndex","SegmentManuIndex %s out of DetectionElement mapping %s",
	    ManuIndexName,fDetectionElementType.Data());
    return 0x0;
  }
}

//_________________________________________________
AliMUONSegmentIndex * AliMUONSegmentationDetectionElement::FindIndex(const char* IndexName) const
{
  // Getting 
  if (fMapIndexPosition) return (AliMUONSegmentIndex *) fMapIndexPosition->FindObject(IndexName);
  else {
    Warning("FindIndex","SegmentIndex %s out of DetectionElement mapping %s",
	    IndexName,fDetectionElementType.Data());
    return 0x0;
  }
}
//_________________________________________________
AliMUONSegmentIndex * AliMUONSegmentationDetectionElement::FindIndexFromPosition(Float_t /*x*/, Float_t /*y*/, Int_t /*cathode*/ ) const
{

 //  char * name = AliMUONSegmentPosition::Name(x,y);

//   if (GetIndexFromPosition( AliMUONSegmentPosition::Name(x,y)) ) 
//     return GetIndexFromPosition( AliMUONSegmentPosition::Name(x,y) );
  
//   Float_t xl= ((Int_t) x*10 )/10.;
//   Float_t yl= ((Int_t) y*10 )/10.;
//   Int_t ix,iy, ixp;
  
//   for(ix=1; ix<4; ix++) {
//     xl = ((Int_t) 0.5+10.*(x +  ((Float_t) ix )*0.1))/10.; 
//     for(iy=-ix; iy<ix+1; iy++) {
//       printf("A %d and %d and %f and %f \n",ix,iy,xl,yl);
//       yl = ((Int_t) 10.*(y +  ((Float_t) iy )*0.1))/10. ;
//       sprintf(name,"%5.1f-%5.1f",xl, yl);
//       if (GetIndexFromPosition(name)) break;
//     }
//     if (GetIndexFromPosition(name)) break;
    

//     for(ixp=ix-1; ixp>-ix-1; ixp--) {
//       xl = ((Int_t) 0.5+10.*(x +  ((Float_t) ixp )*0.1))/10. ;
//       printf("B %d and %d and %f and %f \n",ixp,ix, xl, yl);
//       sprintf(name,"%5.1f-%5.1f",xl, yl);
//       if (GetIndexFromPosition(name)) break;
//     }
//     if (GetIndexFromPosition(name)) break;
    
//     for(iy=ix-1; iy>-ix-1; iy--) {
//       yl = ((Int_t) 0.5+10.*(y +  ((Float_t) iy )*0.1))/10. ;
//       printf("C %d and %d and %f and %f \n",-ix,iy,xl,yl);
//       sprintf(name,"%5.1f-%5.1f",xl, yl);
//       if (GetIndexFromPosition(name)) break;
//     }
//     if (GetIndexFromPosition(name)) break;
    
//     for(ixp=-ix+1; ixp<ix+1; ixp++) {
//       xl = ((Int_t) 0.5+10.*(x +  ((Float_t) ixp )*0.1))/10. ;
//       printf("D %d and %d and %f and %f \n",ixp,-ix,xl,yl);
//       sprintf(name,"%5.1f-%5.1f",xl, yl);
//       if (GetIndexFromPosition(name)) break;
//     }
//     if (GetIndexFromPosition(name)) break;
//   }
//   return GetIndexFromPosition(name);
  return 0x0;
}
//_________________________________________________
void    AliMUONSegmentationDetectionElement::Init(const char * DetectionElementType)
{
  TString elementType(DetectionElementType);
  fSegmentationMappingFileBending = fgkDefaultTop+fgkStationDir+"345"
    +fgkBendingDir+"/"+elementType+fgkBendingExt+fgkFileExt;
  printf("file is %s\n", fSegmentationMappingFileBending.Data());
  fSegmentationMappingFileNonBending = fgkDefaultTop+fgkStationDir+"345"
    +fgkNonBendingDir+"/"+elementType+fgkNonBendingExt+fgkFileExt;

  if (fMapManuIndexIndex==0x0) { 
    fListOfIndexes = new TObjArray(15000);
    fListOfManuIndexes =new TObjArray(15000);
    fListOfPositions =new TObjArray(15000);
    fMapManuIndexIndex= new TMap();
    fMapIndexManuIndex = new TMap();
    fMapIndexPosition = new TMap();
    fMapPositionIndex = new TMap();
  }
  else {
    fListOfIndexes->Delete();
    fListOfManuIndexes->Delete();
    fListOfPositions->Delete();
    fMapManuIndexIndex->Clear();
    fMapIndexManuIndex->Clear();
    fMapIndexPosition->Clear();
    fMapPositionIndex->Clear();
  }
  Int_t icathode;
  //Bendingplane
  icathode=0;
  ReadingSegmentationMappingFile(fSegmentationMappingFileBending ,icathode);
  //NonBendingplane
  icathode=1;
  ReadingSegmentationMappingFile(fSegmentationMappingFileNonBending,icathode);
  
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
      printf("id=%d ix=%d iy=%d x=%f y=%f idmanu=%d and idchannel=%d\n",id, ix, iy,  x, y,idmanu, idchannel);
      
      fListOfIndexes->AddAt(     new AliMUONSegmentIndex(id,ix,iy,cathode), id);
      fListOfManuIndexes->AddAt( new AliMUONSegmentManuIndex(id,idmanu,0,idchannel), id);;
      fListOfPositions->AddAt(   new AliMUONSegmentPosition(id, x, y,cathode), id);;

      ( (AliMUONSegmentIndex* )fListOfIndexes->At(id))->Print();
      ( (AliMUONSegmentManuIndex*)fListOfManuIndexes->At(id))->Print();
      ( (AliMUONSegmentPosition*)fListOfPositions->At(id))->Print();
     
      fMapManuIndexIndex->Add( fListOfManuIndexes->At(id), fListOfIndexes->At(id));
      fMapIndexManuIndex->Add( fListOfIndexes->At(id),     fListOfManuIndexes->At(id));
      fMapIndexPosition->Add(  fListOfIndexes->At(id),     fListOfPositions->At(id));
      fMapPositionIndex->Add(  fListOfPositions->At(id), fListOfIndexes->At(id) );
    } 
  }
  in.close();
}
