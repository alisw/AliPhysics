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
#include "TMath.h"

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
  fWireD  = 0.25; // wire pitch in cm
  fWireX0 = 0.;  // X0 position of the 1st wire
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

//__________________________________________________
AliMUONSegmentIndex *  AliMUONSegmentationDetectionElement::FindIndexFromPosition(Float_t x, Float_t y, Int_t cathode)
{
  // Finding x_wire corresponding to x
  Float_t x_wire = GetAnod(x);

  //Finding pad corresponding to the position (x_wire, y) in a zone of size 3cm x 10cm or 10cm x 3cm depending on cathode plane
  // 
  Int_t ix_max = (cathode==0) ? TMath::Nint(5./AliMUONSegmentPosition::GetUnit())  :  TMath::Nint(1.5/AliMUONSegmentPosition::GetUnit());
  Int_t iy_max = (cathode==0) ? TMath::Nint(1.5/AliMUONSegmentPosition::GetUnit()) :  TMath::Nint(5./AliMUONSegmentPosition::GetUnit());

  AliMUONSegmentIndex * segmentindex =0x0;
  AliMUONSegmentIndex * foundsegmentindex =0x0;
  AliMUONSegmentPosition * segmentposition=0x0;
  Int_t   ix, iy;
  Float_t xt,yt;
  Float_t distance = 99999999.;
  //printf("%d %d \n",ix_max, iy_max);

  for(ix=-ix_max; ix<ix_max; ix++) {
    xt = x_wire + ((Float_t)ix)*AliMUONSegmentPosition::GetUnit();
    for(iy=-iy_max; iy<iy_max; iy++) {
      yt = y + ((Float_t)iy)*AliMUONSegmentPosition::GetUnit();
      segmentindex = GetIndexFromPosition( xt, yt, cathode);
      if ( segmentindex ) {
	// segmentindex->Print();
	segmentposition = GetPosition(segmentindex->GetName());
	if ( segmentposition->Distance(x_wire, y) < distance ) {     
	  //printf("adfafa xt %f yt %f distance %f \n", xt, yt, segmentposition->Distance(xt,yt) );	
	  distance = segmentposition->Distance(x_wire,y);
	  foundsegmentindex = segmentindex;
	}
      }
    }
  }
 if (!foundsegmentindex) {
   Warning("FindIndexFromPosition","Not found Index for position x=%5.2f y=%5.2f \n",x,y);
 }    
 return foundsegmentindex;
}

//____________________________________________________-
Float_t AliMUONSegmentationDetectionElement::GetAnod(Float_t xhit) const
{
  // Returns for a hit position xhit the position of the nearest anode wire    
  Float_t wire= ( (xhit- fWireX0)>0 ) ? 
    Int_t( (xhit- fWireX0)/fWireD )+0.5 : 
    Int_t( (xhit- fWireX0)/fWireD )-0.5;
  return fWireD*wire+fWireX0; 
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
//__________________________________________________
void  AliMUONSegmentationDetectionElement::GetPadC(Int_t ix, Int_t iy, Int_t cathode, Float_t &x, Float_t &y )
{
  x = GetPosition(ix,iy,cathode)->GetXlocal();
  y = GetPosition(ix,iy,cathode)->GetYlocal();
}
//__________________________________________________
void  AliMUONSegmentationDetectionElement::GetPadI(Float_t x, Float_t y, Int_t cathode, Int_t &padx, Int_t &pady)
{

  AliMUONSegmentIndex * segmentindex = FindIndexFromPosition(x,y,cathode);

  if (segmentindex) {
    padx = segmentindex->GetPadX();
    pady = segmentindex->GetPadX();
  }
  else {
    Warning("GetPadI","Not found Index for position x=%5.2f y=%5.2f \n",x,y);
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
  return GetIndexFromPosition( AliMUONSegmentPosition::Name(x, y, cathode).Data() );
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
void    AliMUONSegmentationDetectionElement::Init(const char * DetectionElementType)
{
  TString elementType(DetectionElementType);
  fSegmentationMappingFileBending = fgkDefaultTop+fgkStationDir+"345"
    +fgkBendingDir+"/"+elementType+fgkBendingExt+fgkFileExt;
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
  Info("ReadingSegmentationMappingFile","%s", fSegmentationMappingFileBending.Data());
  ReadingSegmentationMappingFile(fSegmentationMappingFileBending ,icathode);
  //NonBendingplane
  icathode=1;
  Info("Init","Reading mapping file is %s\n", fSegmentationMappingFileNonBending.Data());
  ReadingSegmentationMappingFile(fSegmentationMappingFileNonBending,icathode);
  
}
//_______________________________________________________________
void AliMUONSegmentationDetectionElement::ReadingSegmentationMappingFile(TString infile, Int_t cathode)
{ 
  ifstream in( infile,  ios::in);
  if (!in) {
    Error("ReadingSegmentationMappingFile", "File not found.");
  }
  else {
    Int_t id, ix, iy, idmanu, idchannel;
    Float_t x, y;
    do {
      in >> id >> ix >> iy >> x >> y >> idmanu >> idchannel;
      //     printf("id=%d ix=%d iy=%d x=%f y=%f idmanu=%d and idchannel=%d\n",id, ix, iy,  x, y,idmanu, idchannel);
      
      fListOfIndexes->AddAt(     new AliMUONSegmentIndex(id,ix,iy,cathode), id);
      fListOfManuIndexes->AddAt( new AliMUONSegmentManuIndex(id,idmanu,0,idchannel), id);;
      fListOfPositions->AddAt(   new AliMUONSegmentPosition(id, x, y,cathode), id);;

      //( (AliMUONSegmentIndex* )fListOfIndexes->At(id))->Print();
      //( (AliMUONSegmentManuIndex*)fListOfManuIndexes->At(id))->Print();
      //( (AliMUONSegmentPosition*)fListOfPositions->At(id))->Print();
     
      fMapManuIndexIndex->Add( fListOfManuIndexes->At(id), fListOfIndexes->At(id));
      fMapIndexManuIndex->Add( fListOfIndexes->At(id),     fListOfManuIndexes->At(id));
      fMapIndexPosition->Add(  fListOfIndexes->At(id),     fListOfPositions->At(id));
      fMapPositionIndex->Add(  fListOfPositions->At(id),   fListOfIndexes->At(id) );
    }  
    while ( !in.eof()); 
  }
  in.close();
}
