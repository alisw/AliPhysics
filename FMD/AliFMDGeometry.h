#ifndef ALIFMDGEOMETRY_H
#define ALIFMDGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
//____________________________________________________________________ 
//  
// Forward Multiplicity Detector based on Silicon wafers. 
//
// This class is a singleton that handles the geometry parameters of
// the FMD detectors.  
//                                                       
#ifndef ALIGEOMETRY_H
# include <AliGeometry.h>
#endif
class TVector3;
class TMatrix;
class TParticle;
class AliRecPoint;
class AliFMDRing;
class AliFMDDetector;
class AliFMD1;
class AliFMD2;
class AliFMD3;
#ifndef USE_PRE_MOVE
#ifndef ROOT_TArrayI
# include <TArrayI.h>
#endif
class AliFMDGeometryBuilder;
class TArrayI;
#endif


//__________________________________________________________________
/** Singleton object of FMD geometry descriptions and parameters.
 */
class AliFMDGeometry : public AliGeometry
{
public:
  static AliFMDGeometry* Instance();
  virtual void Init();
  AliFMDRing*     GetInner() const { return fInner; }
  AliFMDRing*     GetOuter() const { return fOuter; }
  AliFMD1*        GetFMD1()  const { return (fUseFMD1 ? fFMD1 : 0); }
  AliFMD2*        GetFMD2()  const { return (fUseFMD2 ? fFMD2 : 0); }
  AliFMD3*        GetFMD3()  const { return (fUseFMD3 ? fFMD3 : 0); }
  AliFMDDetector* GetDetector(Int_t i) const;
  AliFMDRing*     GetRing(Char_t i) const;
  void            Disable(Int_t i);
  void            Enable(Int_t i);
  Double_t        GetSiDensity() const { return 2.33; }
  void            Detector2XYZ(UShort_t detector, Char_t ring, 
			       UShort_t sector, UShort_t strip, 
			       Double_t& x, Double_t& y, Double_t& z) const;
  Bool_t          XYZ2Detector(Double_t x, Double_t y, Double_t z, 
			       UShort_t& detector, Char_t& ring, 
			       UShort_t& sector, UShort_t& strip) const;
#ifndef USE_PRE_MOVE
  void   Build();
  Int_t  GetDetectorOff() const    { return fDetectorOff; }
  Int_t  GetModuleOff() const      { return fModuleOff;   }
  Int_t  GetRingOff() const        { return fRingOff;     }
  Int_t  GetSectorOff() const      { return fSectorOff;   }
  void   SetDetectorOff(Int_t off) { fDetectorOff = off; }
  void   SetModuleOff(Int_t off)   { fModuleOff   = off; }
  void   SetRingOff(Int_t off)     { fRingOff     = off; }
  void   SetSectorOff(Int_t off)   { fSectorOff   = off; }
  Bool_t IsActive(Int_t vol) const;
  void   SetActive(Int_t* active, Int_t n);
  void   AddActive(Int_t id);
  void   SetBuilder(AliFMDGeometryBuilder* b) { fBuilder = b; }
  void   ExtractGeomInfo();
  void   SetDetailed(Bool_t det) { fDetailed = det; }
  Bool_t IsDetailed() const { return fDetailed; }
  void   UseAssembly(Bool_t ass)  { fUseAssembly = ass; }
#endif  

  // AliGeometry member functions 
  virtual void GetGlobal(const AliRecPoint* p, TVector3& pos, 
			 TMatrix& mat) const;
  virtual void GetGlobal(const AliRecPoint* p, TVector3& pos) const;
  virtual Bool_t Impact(const TParticle* particle) const;
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
  AliFMDGeometry();
  AliFMDGeometry(const AliFMDGeometry& other);
  AliFMDGeometry& operator=(const AliFMDGeometry& other);
  virtual ~AliFMDGeometry() {}
  
#ifndef USE_PRE_MOVE
  AliFMDGeometryBuilder* fBuilder;
  Int_t fDetectorOff;
  Int_t fModuleOff;  
  Int_t fRingOff;    
  Int_t fSectorOff;  
  TArrayI fActive;
  Bool_t fDetailed;
  Bool_t fUseAssembly;
#endif

  ClassDef(AliFMDGeometry,1); //
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
