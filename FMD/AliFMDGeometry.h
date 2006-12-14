#ifndef ALIFMDGEOMETRY_H
#define ALIFMDGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDGeometry.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:40:37 2006
    @brief   Geometry mananger for the FMD
*/
//____________________________________________________________________
//                                                                          
// Forward Multiplicity Detector based on Silicon wafers. 
//
// This class is a singleton that handles the geometry parameters of
// the FMD detectors.  
// The actual code is done by various separate classes.
//                                                       
#ifndef ALIGEOMETRY_H
# include <AliGeometry.h>
#endif
#ifndef ROOT_TArrayI
# include <TArrayI.h>
#endif
#ifndef ROOT_TMatrixFfwd
# include <TMatrixFfwd.h>
#endif
class TVector3;
class TParticle;
class AliRecPoint;
class AliFMDRing;
class AliFMDDetector;
class AliFMD1;
class AliFMD2;
class AliFMD3;
class AliFMDGeometryBuilder;


//__________________________________________________________________
/** @brief Singleton object of FMD geometry descriptions and parameters.
    This class is a singleton that handles the geometry parameters of
    the FMD detectors.  
                                                          
    The actual code is done by various separate classes.   Below is
    diagram showing the relationship between the various FMD classes
    that handles the geometry 
    @verbatim
                                  +------------+ 
                               +- | AliFMDRing |
    			   2   |  +------------+
         +----------------+<>--+        |				
         | AliFMDGeometry |             ^                       	
         +----------------+<>--+        V 1..2                     	
              		   3   | +----------------+ 		
               		       +-| AliFMDDetector | 		
                		 +----------------+		
                                        ^
                                        |
                          +-------------+-------------+
                          |             |             |	      
                     +---------+   +---------+   +---------+
                     | AliFMD1 |   | AliFMD2 |   | AliFMD3 |
                     +---------+   +---------+   +---------+
         
    @endverbatim
    -  AliFMDRing 
       This class contains all stuff needed to do with a ring.  It's
       used by the AliFMDDetector objects to instantise inner and
       outer rings.  The AliFMDRing objects are shared by the
       AliFMDDetector objects, and owned by the AliFMDv1 object. 
    -  AliFMD1, AliFMD2, and AliFMD3 
       These are specialisation of AliFMDDetector, that contains the
       particularities of each of the sub-detector system.  It is
       envisioned that the classes should also define the support
       volumes and material for each of the detectors.

    @ingroup FMD_base
 */
class AliFMDGeometry : public AliGeometry
{
public:
  /** @return Singleton */
  static AliFMDGeometry* Instance();
  /** Initialize */
  virtual void Init();
  /** Initialize transforms */
  virtual void InitTransformations();
  /** @return Get inner description */
  AliFMDRing*     GetInner() const { return fInner; }
  /** @return Get outer description */
  AliFMDRing*     GetOuter() const { return fOuter; }
  /** @return Get FMD1 description */
  AliFMD1*        GetFMD1()  const { return (fUseFMD1 ? fFMD1 : 0); }
  /** @return Get FMD2 description */
  AliFMD2*        GetFMD2()  const { return (fUseFMD2 ? fFMD2 : 0); }
  /** @return Get FMD3 description */
  AliFMD3*        GetFMD3()  const { return (fUseFMD3 ? fFMD3 : 0); }
  /** Get description of a sub-detector
      @param i Sub-detector #
      @return Description of sub-detector, or 0 */
  AliFMDDetector* GetDetector(Int_t i) const;
  /** Get description of a ring
      @param i Ring id
      @return Description of ring, or 0 */
  AliFMDRing*     GetRing(Char_t i) const;
  /** @param i IF true, disable sub-detector @a i */
  void            Disable(Int_t i);
  /** @param i IF true, enable sub-detector @a i */
  void            Enable(Int_t i);
  /** @return Density @f$ \rho@f$ of silicon */
  Double_t        GetSiDensity() const { return 2.33; }
  /** Translate detector coordinates (detector, ring, sector, strip)
      to spatial coordinates (x, y, z) in the master reference frame
      of ALICE.  The member function uses the transformations
      previously obtained from the TGeoManager.
      @param detector Detector number
      @param ring     Ring id
      @param sector   Sector number
      @param strip    Strip number
      @param x        On return, X coordinate 
      @param y        On return, Y coordinate 
      @param z        On return, Z coordinate  */
  void            Detector2XYZ(UShort_t detector, Char_t ring, 
			       UShort_t sector, UShort_t strip, 
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
      @param x        X coordinate
      @param y 	      Y coordinate
      @param z 	      Z coordinate
      @param detector On return, Detector number
      @param ring     On return, Ring id		   
      @param sector   On return, Sector number	   
      @param strip    On return, Strip number	   
      @return @c  false of (@a x, @a y, @a z) is not within this
      detector.  */
  Bool_t          XYZ2Detector(Double_t x, Double_t y, Double_t z, 
			       UShort_t& detector, Char_t& ring, 
			       UShort_t& sector, UShort_t& strip) const;
  /** Make the geometry.  This delegates to AliFMDGeometryBuilder */
  void   Build();
  /** @return Get detector offset in paths */
  Int_t  GetDetectorOff() const    { return fDetectorOff; }
  /** @return Get sensor offset in paths */
  Int_t  GetModuleOff() const      { return fModuleOff;   }
  /** @return Get ring offset in paths */
  Int_t  GetRingOff() const        { return fRingOff;     }
  /** @return Get ring sector in paths */
  Int_t  GetSectorOff() const      { return fSectorOff;   }
  /** @param off Detector off-set set in geometry path */
  void   SetDetectorOff(Int_t off) { fDetectorOff = off; }
  /** @param off Module off-set set in geometry path */
  void   SetModuleOff(Int_t off)   { fModuleOff   = off; }
  /** @param off Ring off-set set in geometry path */
  void   SetRingOff(Int_t off)     { fRingOff     = off; }
  /** @param off Sectord off-set set in geometry path */
  void   SetSectorOff(Int_t off)   { fSectorOff   = off; }
  /** Check if volume @a vol is marked as active 
      @param vol Volume ID
      @return  @c true if @a vol is declared active */
  Bool_t IsActive(Int_t vol) const;
  /** Set active volumes 
      @param active Active volume id array 
      @param n elements of @a active */
  void   SetActive(Int_t* active, Int_t n);
  /** @param id Register volume @a id to be active */
  void   AddActive(Int_t id);
  /** Set an external geometry builder
      @param b Geometry builder */
  void   SetBuilder(AliFMDGeometryBuilder* b) { fBuilder = b; }
  /** Extract informaton from TGeoManager */
  void   ExtractGeomInfo();
  /** Whether we are to use a detailed geometry or not
      @param det if @c true, make a detailed geometry. */
  void   SetDetailed(Bool_t det) { fDetailed = det; }
  /** @return @c true if geometry is detailed */
  Bool_t IsDetailed() const { return fDetailed; }
  /** @param ass Whether to use assemblies or not */
  void   UseAssembly(Bool_t ass)  { fUseAssembly = ass; }

  // AliGeometry member functions 
  /** Get global coordinates cooresponding to a rec point. 
      @param p   Reconstructed point.
      @param pos On return, the position
      @param mat On return, the material at @a post */
  virtual void GetGlobal(const AliRecPoint* p, TVector3& pos, 
			 TMatrixF& mat) const;
  /** Get global coordinates cooresponding to a rec point. 
      @param p   Reconstructed point.
      @param pos On return, the position */
  virtual void GetGlobal(const AliRecPoint* p, TVector3& pos) const;
  /** Check if particle will hit an active detector element.  Note
      done yet. 
      @param particle Track 
      @return @c true if @a particle will hit this detector */
  virtual Bool_t Impact(const TParticle* particle) const;
  /** Declare alignable volumes */
  virtual void SetAlignableVolumes() const;
protected:
  Bool_t        fIsInitialized; // Whether singleton is initalized
  AliFMDRing*	fInner;		// Inner ring geometry information
  AliFMDRing*	fOuter;		// Outer ring geometry information
  AliFMD1*	fFMD1;		// FMD1 geometry information
  AliFMD2*	fFMD2;		// FMD2 geometry information
  AliFMD3*	fFMD3;		// FMD3 geometry information
  Bool_t	fUseFMD1;	// Wheter to Use FMD1 or not
  Bool_t	fUseFMD2;	// Wheter to Use FMD2 or not
  Bool_t	fUseFMD3;	// Wheter to Use FMD3 or not
  static AliFMDGeometry* fgInstance; // Singleton instance 
  /** CTOR */
  AliFMDGeometry();
  /** Copy CTOR
      @param other To copy from  */
  AliFMDGeometry(const AliFMDGeometry& other);
  /** Assignment operator 
      @param other To assig from
      @return reference to this.  */
  AliFMDGeometry& operator=(const AliFMDGeometry& other);
  virtual ~AliFMDGeometry() {}
  
  AliFMDGeometryBuilder* fBuilder; // Geometry builder 
  Int_t fDetectorOff;              // Detector off-set 
  Int_t fModuleOff;                // Module off-set 
  Int_t fRingOff;                  // ring offset
  Int_t fSectorOff;                // Sector offset    
  TArrayI fActive;                 // Active volumes
  Bool_t fDetailed;                // Whether to make detailed geom
  Bool_t fUseAssembly;             // Whther to use assemblies 

  ClassDef(AliFMDGeometry,1); // Geometry parameters and manager 
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
