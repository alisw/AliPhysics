#ifndef ALIFMDSUBDETECTOR_H
#define ALIFMDSUBDETECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ROOT_TObject
# include <TObject.h>
#endif

class TNode;
class TList;
class AliFMDRing;


//__________________________________________________________________
struct AliFMDSubDetector : public TObject
{
protected:
  Int_t       fId;                 // Detector number 
  Double_t    fInnerZ;             // Position of outer ring along z
  Double_t    fOuterZ;             // Position of outer ring along z 
  Double_t    fHoneycombThickness; // Thickness of honeycomb plate 
  Double_t    fKaptionThickness;   // Thickness of kaption of honeycomb
  Double_t    fInnerHoneyLowR;     // Inner radius of inner honeycomb
  Double_t    fInnerHoneyHighR;    // Outer radius of inner honeycomb
  Int_t       fInnerHoneyTopId;    // Volume ID of top of inner honeycomb
  Int_t       fInnerHoneyBottomId; // Volume ID of bottom of inner honeycomb
  Double_t    fOuterHoneyLowR;     // Inner radius of outer honeycomb
  Double_t    fOuterHoneyHighR;    // Outer radius of outer honeycomb
  Int_t       fOuterHoneyTopId;    // Volume ID of top of outer honeycomb   
  Int_t       fOuterHoneyBottomId; // Volume ID of bottom of outer honeycomb

  Int_t       fRotationId;         // The ID of the sub-detector rotation
  
  AliFMDRing* fInner;              // Reference to inner ring description
  AliFMDRing* fOuter;              // Reference to outer ring description
public:  
  AliFMDSubDetector(Int_t n);
  virtual ~AliFMDSubDetector() {}
  virtual void   SetupGeometry(Int_t airId, Int_t kaptionId);  
  virtual void   Geometry(const char* mother, Int_t pbRotId, 
			  Int_t idRotId, Double_t z=0);
  virtual void   SimpleGeometry(TList* nodes, 
				TNode* mother, 
				Int_t colour, 
				Double_t zMother);
  
  virtual void   Gsatt();
  virtual void   Draw(Option_t* option="BIOL0") const; //*MENU*
  virtual Bool_t CheckHit(Char_t ring, Int_t module, Double_t x, Double_t y);

  void   SetInner(AliFMDRing* r)             { fInner = r; }
  void   SetOuter(AliFMDRing* r)             { fOuter = r; }
  void   SetInnerZ(Double_t z)               { fInnerZ = z; }
  void   SetOuterZ(Double_t z)               { fOuterZ = z; }
  void   SetHoneycombThickness(Double_t t=1) { fHoneycombThickness = t; }
  void   SetInnerHoneyLowR(Double_t r)       { fInnerHoneyLowR = r; }
  void   SetInnerHoneyHighR(Double_t r)      { fInnerHoneyHighR = r; }
  void   SetOuterHoneyLowR(Double_t r)       { fOuterHoneyLowR = r; }
  void   SetOuterHoneyHighR(Double_t r)      { fOuterHoneyHighR = r; }
  void   SetKaptionThickness(Double_t t=.1)  { fKaptionThickness = t; }
     
  Double_t    GetInnerZ()             const { return fInnerZ; }
  Double_t    GetOuterZ()             const { return fOuterZ; }
  AliFMDRing* GetInner()              const { return fInner; }
  AliFMDRing* GetOuter()              const { return fOuter; }
  Double_t    GetHoneycombThickness() const { return fHoneycombThickness; }
  Double_t    GetInnerHoneyLowR()     const { return fInnerHoneyLowR; }
  Double_t    GetInnerHoneyHighR()    const { return fInnerHoneyHighR; }
  Double_t    GetOuterHoneyLowR()     const { return fOuterHoneyLowR; }
  Double_t    GetOuterHoneyHighR()    const { return fOuterHoneyHighR; }
  Double_t    GetKaptionThickness()   const { return fKaptionThickness; }
  Int_t       GetId()                 const { return fId; }     
  Bool_t      IsFolder()              const { return kTRUE; }

  ClassDef(AliFMDSubDetector, 1) // FMD Sub detector base class
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//
