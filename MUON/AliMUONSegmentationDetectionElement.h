#ifndef ALIMUONSEGMENTATIONDETECTIONELEMENT_H
#define ALIMUONSEGMENTATIONDETECTIONELEMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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


#include <TObject.h>
#include <TString.h>
#include <TArrayF.h>

class TClonesArray;
class TMap;

class AliMUONSegmentManuIndex;
class AliMUONSegmentPosition;
class AliMUONSegmentIndex;

class AliMUONSegmentationDetectionElement : public TObject {
 public:
  AliMUONSegmentationDetectionElement();
  //AliMUONSegmentationDetectionElement(const char* ElementType="");
  virtual ~AliMUONSegmentationDetectionElement();

  // User functions
  AliMUONSegmentIndex     * GetIndex(Int_t manu, Int_t channel) const;
  AliMUONSegmentIndex     * GetIndexFromPosition(Float_t x, Float_t y, Int_t icathode) const;
  AliMUONSegmentManuIndex * GetManuIndex( Int_t padx, Int_t pady, Int_t cathode) const ;
  AliMUONSegmentPosition  * GetPosition(Int_t padx, Int_t pady, Int_t cathode) const ;

  
  AliMUONSegmentIndex     * GetIndex( const char * SegmentManuIndexName)const;
  AliMUONSegmentIndex     * GetIndexFromPosition( const char * PositionName)const;
  AliMUONSegmentManuIndex * GetManuIndex( const char * SegmentIndexName) const;
  AliMUONSegmentPosition  * GetPosition( const char * SegmentIndexName) const;
  TObjArray *            ListOfIndexes() {return fListOfIndexes;}
  TObjArray *            ListOfManuIndexes() {return fListOfIndexes;}
  TObjArray *            ListOfPositions() {return fListOfIndexes;}
 
  AliMUONSegmentManuIndex * FindManuIndex(const char* ManuIndexName="") const;
  AliMUONSegmentIndex * FindIndex(const char* IndexName="") const;

  AliMUONSegmentIndex * FindIndexFromPosition(Float_t x, Float_t y, Int_t cathode) const;
  
  void     Init(const char * DetectionElementType="slat220000N");

  void     ReadingSegmentationMappingFile(TString infile, Int_t cathode);
  
 protected:
  AliMUONSegmentationDetectionElement(const AliMUONSegmentationDetectionElement& rhs);
  
 private:
  // static data members  
  static const TString fgkDefaultTop;  // Top directory of $Alice_ROOT/MUON/mapping
  static const TString fgkStationDir;  // Directory for station station1, station2, station345
  static const TString fgkBendingDir;    //bending plane directory
  static const TString fgkNonBendingDir; //non-bending plane directory
  static const TString fgkFileExt;  // File extention
  static const TString fgkBendingExt;  // bending file extention
  static const TString fgkNonBendingExt;  // bending file extention

  // data members
  TString   fDetectionElementType;               //  Type of detection element St1Sector, slat220000N, etc ....
  TString   fSegmentationMappingFileBending;    //  Segmentation & mapping file for bending plane
  TString   fSegmentationMappingFileNonBending; //  Segmentation & mapping file for non bending plane
  TObjArray * fListOfIndexes;        // TObject Array fo AliMUONSegmentIndex
  TObjArray * fListOfManuIndexes;   // TObject Array fo AliMUONSegmentManuIndex
  TObjArray * fListOfPositions;  // TObject Array fo AliMUONSegmentPositions
  TMap *    fMapManuIndexIndex;  // Map with key ManuIndex and value = Index
  TMap *    fMapIndexManuIndex;// Map with key ManuIndexIndex and value = ManuIndex
  TMap *    fMapIndexPosition;// Map with key Index and value = Position
  TMap *    fMapPositionIndex;// Map with key Index and value = Position
  

  ClassDef(AliMUONSegmentationDetectionElement,1) // Segmentation for MUON detection elements
    
    };
#endif






