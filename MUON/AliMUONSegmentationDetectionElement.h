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

#include <TString.h>

#include "AliSegmentation.h"

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

  // User common functions

  AliMUONSegmentIndex     * FindIndexFromPosition(Float_t x, Float_t y, Int_t cathode);
  AliMUONSegmentIndex     * GetIndex(Int_t manu, Int_t channel) const;
  AliMUONSegmentIndex     * GetIndexFromPosition(Float_t x, Float_t y, Int_t cathode) const;
  AliMUONSegmentManuIndex * GetManuIndex( Int_t padx, Int_t pady, Int_t cathode) const ;
  AliMUONSegmentPosition  * GetPosition(Int_t padx, Int_t pady, Int_t cathode) const ;
  TObjArray *               ListOfIndexes()     {return fListOfIndexes;}
  TObjArray *               ListOfManuIndexes() {return fListOfManuIndexes;}
  TObjArray *               ListOfPositions()   {return fListOfPositions;}
  
  void                      Init(const char * DetectionElementType="slat220000N");
  void                      ReadingSegmentationMappingFile(TString infile, Int_t cathode);

  // virtual functions from AliSegmentation. In future this class should derive from AliSegmentation
  Float_t      GetAnod(Float_t xhit) const; // Anod wire coordinate closest to xhit
  void         SetDAnod(Float_t D) {fWireD = D;};  // Wire Pitch
  void         GetPadI(Float_t x, Float_t y , Int_t cathode, Int_t &padx, Int_t &pady); // Transform from Position to closest Index 
  void         GetPadC(Int_t ix, Int_t iy, Int_t cathode, Float_t &x, Float_t &y );  // Transform from Index to Position 
/*   void      FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy);// Initialiser  */
/*   void      NextPad();  // Stepper  */
/*   Int_t     MorePads(); // Condition  */
/*   void      Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]); // Get next neighbours  */ 

 protected:  
  AliMUONSegmentIndex     * GetIndex( const char * SegmentManuIndexName)const;
  AliMUONSegmentIndex     * GetIndexFromPosition( const char * PositionName)const;
  AliMUONSegmentManuIndex * GetManuIndex( const char * SegmentIndexName) const;
  AliMUONSegmentPosition  * GetPosition( const char * SegmentIndexName) const;
  AliMUONSegmentManuIndex * FindManuIndex(const char* ManuIndexName="") const;
  AliMUONSegmentIndex     * FindIndex(const char* IndexName="") const;

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
  Float_t     fWireD;         // Wire pitch in cm
  Float_t     fWireX0;        // Initial X0 position in local coordinates of the first wire
  Int_t       fCurrentSegment;// Index of the current segment
  TString     fDetectionElementType;               //  Type of detection element St1Sector, slat220000N, etc ....
  TString     fSegmentationMappingFileBending;    //  Segmentation & mapping file for bending plane
  TString     fSegmentationMappingFileNonBending; //  Segmentation & mapping file for non bending plane
  TObjArray * fListOfIndexes;        // TObject Array fo AliMUONSegmentIndex
  TObjArray * fListOfManuIndexes;   // TObject Array fo AliMUONSegmentManuIndex
  TObjArray * fListOfPositions;  // TObject Array fo AliMUONSegmentPositions
  TMap *      fMapManuIndexIndex;  // Map with key ManuIndex and value = Index
  TMap *      fMapIndexManuIndex;// Map with key ManuIndexIndex and value = ManuIndex
  TMap *      fMapIndexPosition;// Map with key Index and value = Position
  TMap *      fMapPositionIndex;// Map with key Index and value = Position

  ClassDef(AliMUONSegmentationDetectionElement,1) // Segmentation for MUON detection elements
    
    };
#endif






