#ifndef ALIFMDDETECTOR_H
#define ALIFMDDETECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
//__________________________________________________________________
//
// Utility class to help implement the FMD geometry.  This provides
// the interface for the concrete geometry implementations of the FMD
// sub-detectors. 
/** @file    AliFMDDetector.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:36:27 2006
    @brief   Sub-detector base class declaration
    @ingroup FMD_base
*/
#ifndef ROOT_TNamed
# include <TNamed.h>
#endif
class AliFMDRing;
class TGeoMatrix;

/** @defgroup FMD_base Basic classes */
//__________________________________________________________________
/** @brief Base class for the geometry description and parameters of
    the FMD sub detectors FMD1, FMD2, and FMD3.

    This class hold common parameters of the specific FMD detectors.
    @ingroup FMD_base
*/
class AliFMDDetector : public TNamed 
{
public:
  /** Constructor
      @param id    Detector number
      @param inner Pointer to inner ring geometry
      @param outer Pointer to inner outer geometry
      @return  */
  AliFMDDetector(Int_t id, AliFMDRing* inner, AliFMDRing* outer);
  /** Copy CTOR 
      @param other Object to copy from. */
  AliFMDDetector(const AliFMDDetector& other);
  /** Assignment operator 
      @param other Object to assign from
      @return reference to this object  */
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
  
  /** Translate detector coordinates (detector, ring, sector, strip)
      to spatial coordinates (x, y, z) in the master reference frame
      of ALICE.  The member function uses the transformations
      previously obtained from the TGeoManager.
      @param ring   Ring id
      @param sector Sector number
      @param strip  Strip number
      @param x      On return, X coordinate 
      @param y      On return, Y coordinate 
      @param z      On return, Z coordinate  */
  void Detector2XYZ(Char_t ring, UShort_t sector, UShort_t strip, 
		    Double_t& x, Double_t& y, Double_t& z) const;
  /** Translate spatial coordinates (x,y,z) in the master reference
      frame of ALICE to the detector coordinates (detector, ring,
      sector, strip).  Note, that if this method is to be used in
      reconstruction or the like, then the input z-coordinate should
      be corrected for the events interactions points z-coordinate,
      like  
      @code 
      geom->XYZ2Detector(x,y,z-ipz,d,r,s,t);
      @endcode
      @param x       X coordinate
      @param y 	     Y coordinate
      @param z 	     Z coordinate
      @param ring    On return, Ring id		   
      @param sector  On return, Sector number	   
      @param strip   On return, Strip number	   
      @return @c  false of (@a x, @a y, @a z) is not within this
      detector.  */
  Bool_t XYZ2Detector(Double_t x, Double_t y, Double_t z, 
		      Char_t& ring, UShort_t& sector, UShort_t& strip) const;

  /** Declare alignable volumes */
  virtual void SetAlignableVolumes() const;
protected:
  /** Check if we have all transformations for a ring 
      @param ring Ring to check for 
      @return @c true if we got all transforms  */
  Bool_t        HasAllTransforms(Char_t ring) const;
  /** Get transformation matrix for a sector in a ring 
      @param ring   Ring id
      @param sector Sector numberr 
      @return Matrix on success, 0 otherwise */
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
  TObjArray*    fInnerTransforms;       // List of inner module global
  TObjArray*    fOuterTransforms;       // List of outer module global

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
