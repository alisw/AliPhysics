// -*- mode: c++ -*-
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ALIFMDRING_H
#define ALIFMDRING_H
#ifndef ALIFMDPOLYGON_H
# include <AliFMDPolygon.h>
#endif
#ifndef ROOT_TArrayI
# include <TArrayI.h>
#endif

class TBrowser;
class TNode;
class TObjArray;
class TShape;
class TList;


//__________________________________________________________________
struct AliFMDRing : public TObject
{
  Char_t   fId;			 // ID
  Bool_t   fDetailed;
  Int_t    fActiveId;		 // Active volume 
  Int_t    fPrintboardBottomId;  // Print board bottom volume
  Int_t    fPrintboardTopId;     // Print board top volume
  Int_t    fRingId;		 // Ring volume
  Int_t    fSectionId;		 // Section volumes 
  Int_t    fStripId;		 // Strip volumes 
  Int_t    fVirtualBackId;	 // Virtual Back volume
  Int_t    fVirtualFrontId;	 // Virtual Front volume

  Double_t fBondingWidth;	 // With of bonding pad on sensor
  Double_t fWaferRadius;	 // Size of wafer the sensor was made from 
  Double_t fSiThickness;	 // Thickness of sensor
  Double_t fLowR;		 // Lower radius of ring
  Double_t fHighR;		 // Upper radius of ring
  Double_t fTheta;		 // Opening angle of the silicon wafers
  Int_t    fNStrips;		 // Number of strips 
  Double_t fRingDepth;           // How far the ring extends beyond
				 // the z value given. 
  Double_t fLegRadius;		 // Radius of support legs 
  Double_t fLegLength;		 // Radius of support legs 
  Double_t fLegOffset;		 // Radius of support legs 

  Double_t fModuleSpacing;	 // Staggering offset 
  Double_t fPrintboardThickness; // Thickness of print board

  TArrayI    fRotations;	 // Array of rotations
  TShape*    fShape;             // Shape used for event display
  TObjArray* fRotMatricies;      // Matricies used for event display

  AliFMDPolygon  fPolygon;		 // Polygon shape 
public:
  //----------------------------------------------------------------
  AliFMDRing(Char_t id='\0', Bool_t detailed=kTRUE);
  virtual ~AliFMDRing();
  void   Init();
  bool   IsWithin(size_t moduleNo, double x, double y) const;
  void   SetupCoordinates();  
  void   SetupGeometry(Int_t vacuumId, Int_t siId, Int_t pcbId, 
		       Int_t pbRotId, Int_t idRotId);
  void   Geometry(const char* mother, Int_t baseId, Double_t z, Int_t pbRotId,
		  Int_t idRotId);
  void   SimpleGeometry(TList* nodes, 
			TNode* mother, 
			Int_t colour, 
			Double_t z, 
			Int_t n);
  void   Gsatt();
  void   Draw(Option_t* opt="HOL") const; //*MENU*
  void   Browse(TBrowser* b);
  Bool_t IsFolder() const { return kTRUE; }
  
  Char_t   GetId()                  const { return fId; }
  Int_t    GetActiveId()	    const { return fActiveId; }
  Int_t    GetPrintboardBottomId()  const { return fPrintboardBottomId; }
  Int_t    GetPrintboardTopId()     const { return fPrintboardTopId; }
  Int_t    GetRingId()		    const { return fRingId; }
  Int_t    GetSectionId()           const { return fSectionId; }
  Int_t    GetStripId()		    const { return fStripId; }
  Int_t    GetVirtualBackId()	    const { return fVirtualBackId; }
  Int_t    GetVirtualFrontId()	    const { return fVirtualFrontId; }
  Double_t GetBondingWidth()	    const { return fBondingWidth; }
  Double_t GetWaferRadius()	    const { return fWaferRadius; }
  Double_t GetSiThickness()	    const { return fSiThickness; }
  Double_t GetLowR()		    const { return fLowR; }
  Double_t GetHighR()		    const { return fHighR; }
  Double_t GetTheta()		    const { return fTheta; }
  Int_t    GetNStrips()		    const { return fNStrips; }
  Int_t    GetNSectors()	    const { return Int_t(360 / fTheta); }
  Double_t GetLegRadius()	    const { return fLegRadius; }
  Double_t GetLegLength()	    const { return fLegLength; }
  Double_t GetModuleSpacing()	    const { return fModuleSpacing; }
  Double_t GetPrintboardThickness() const { return fPrintboardThickness; }
  Double_t GetRingDepth()           const { return fRingDepth; }

  void SetBondingWidth(Double_t width)	          { fBondingWidth = width; }
  void SetWaferRadius(Double_t radius)	          { fWaferRadius = radius; } 
  void SetSiThickness(Double_t thickness)	  { fSiThickness = thickness; }
  void SetLowR(Double_t lowR)		          { fLowR = lowR; }
  void SetHighR(Double_t highR)		          { fHighR = highR; }
  void SetTheta(Double_t theta)		          { fTheta = theta; }
  void SetNStrips(Int_t nStrips)		  { fNStrips = nStrips; }
  void SetLegRadius(Double_t radius)	          { fLegRadius = radius; }
  void SetLegLength(Double_t length)	          { fLegLength = length; }
  void SetLegOffset(Double_t offset)	          { fLegOffset = offset; }

  void SetModuleSpacing(Double_t       spacing)	  { fModuleSpacing = spacing; }
  void SetPrintboardThickness(Double_t thickness) { fPrintboardThickness = thickness; }

  ClassDef(AliFMDRing, 1) // Ring volume parameters 
};
#endif 
//
// EOF
//
