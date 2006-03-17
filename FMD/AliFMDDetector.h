#ifndef ALIFMDDETECTOR_H
#define ALIFMDDETECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ROOT_TNamed
# include <TNamed.h>
#endif
class AliFMDRing;
class TGeoMatrix;


//__________________________________________________________________
/** Base class for the geometry description and parameters of the FMD
    sub detectors FMD1, FMD2, and FMD3.

    This class hold common parameters of the specific FMD detectors.
*/
class AliFMDDetector : public TNamed 
{
public:
  AliFMDDetector(Int_t id, AliFMDRing* inner, AliFMDRing* outer);
  AliFMDDetector(const AliFMDDetector& other);
  AliFMDDetector& operator=(const AliFMDDetector& other);
  virtual ~AliFMDDetector() {}
  /** Initialize the geometry */
  virtual void Init();
  /** Find the transformations that correspond to modules of this
      detector, and store them in the arrays. */
  virtual void InitTransformations();
  
  /** @param x Detector number */
  void SetId(Int_t x) { fId = x; }
  /** @param x Position of outer ring along z */
  void SetInnerZ(Double_t x) { fInnerZ = x; }
  /** @param x Position of outer ring along z */
  void SetOuterZ(Double_t x) { fOuterZ = x; }
  /** @param x Thickness of honeycomb plate */
  void SetHoneycombThickness(Double_t x=1) { fHoneycombThickness = x; }
  /** @param x Thickness of aluminium of honeycomb */
  void SetAlThickness(Double_t x=.1) { fAlThickness = x; }
  /** @param x Inner radius of inner honeycomb */
  void SetInnerHoneyLowR(Double_t x) { fInnerHoneyLowR = x; }
  /** @param x Outer radius of inner honeycomb */
  void SetInnerHoneyHighR(Double_t x) { fInnerHoneyHighR = x; }
  /** @param x Inner radius of outer honeycomb */
  void SetOuterHoneyLowR(Double_t x) { fOuterHoneyLowR = x; }
  /** @param x Outer radius of outer honeycomb */
  void SetOuterHoneyHighR(Double_t x) { fOuterHoneyHighR = x; }
    
  /** @return Detector number */
  Int_t GetId() const { return fId; }
  /** @return Position of outer ring along z */
  Double_t GetInnerZ() const { return fInnerZ; }
  /** @return Position of outer ring along z */
  Double_t GetOuterZ() const { return fOuterZ; }
  /** @return Thickness of honeycomb plate */
  Double_t GetHoneycombThickness() const { return fHoneycombThickness; }
  /** @return Thickness of aluminium of honeycomb */
  Double_t GetAlThickness() const { return fAlThickness; }
  /** @return Inner radius of inner honeycomb */
  Double_t GetInnerHoneyLowR() const { return fInnerHoneyLowR; }
  /** @return Outer radius of inner honeycomb */
  Double_t GetInnerHoneyHighR() const { return fInnerHoneyHighR; }
  /** @return Inner radius of outer honeycomb */
  Double_t GetOuterHoneyLowR() const { return fOuterHoneyLowR; }
  /** @return Outer radius of outer honeycomb */
  Double_t GetOuterHoneyHighR() const { return fOuterHoneyHighR; }
    
  /** @return Inner ring information */
  AliFMDRing* GetInner() const { return fInner; }
  /** @return Outer ring information */
  AliFMDRing* GetOuter() const { return fOuter; }
  /** @param id Id of ring to get 
      @return Pointer to ring, 0 on failure */
  AliFMDRing* GetRing(Char_t id) const;
  /** @param id Id of ring to get 
      @return Z position of ring or 0 on failure */
  Double_t GetRingZ(Char_t id) const;
  
  void Detector2XYZ(Char_t ring, UShort_t sector, UShort_t strip, 
		    Double_t& x, Double_t& y, Double_t& z) const;
  Bool_t XYZ2Detector(Double_t x, Double_t y, Double_t z, 
		      Char_t& ring, UShort_t& sector, UShort_t& strip) const;
protected:
  Bool_t        HasAllTransforms(Char_t ring) const;
  TGeoMatrix*   FindTransform(Char_t ring, UShort_t sector) const;
  Int_t		fId;			// Detector number
  Double_t	fInnerZ;		// Position of outer ring along z
  Double_t	fOuterZ;		// Position of outer ring along z
  Double_t	fHoneycombThickness;	// Thickness of honeycomb plate
  Double_t	fAlThickness;		// Thickness of aluminium of honeycomb
  Double_t	fInnerHoneyLowR;	// Inner radius of inner honeycomb
  Double_t	fInnerHoneyHighR;	// Outer radius of inner honeycomb
  Double_t	fOuterHoneyLowR;	// Inner radius of outer honeycomb
  Double_t	fOuterHoneyHighR;	// Outer radius of outer honeycomb
  AliFMDRing*	fInner;			// Pointer to inner ring information
  AliFMDRing*	fOuter;			// Pointer to outer ring information
  TObjArray*    fInnerTransforms;       // List of inner module <-> global
  TObjArray*    fOuterTransforms;       // List of outer module <-> global

  ClassDef(AliFMDDetector, 1); // 
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
