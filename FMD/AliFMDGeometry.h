#ifndef ALIFMDGEOMETRY_H
#define ALIFMDGEOMETRY_H
//____________________________________________________________________ 
//  
// $Id$ 
//
#ifndef ALIGEOMETRY_H
# include <AliGeometry.h>
#endif
#ifndef ROOT_TObjArray
# include <TObjArray.h>
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
  void Detector2XYZ(UShort_t detector, 
		    Char_t ring, UShort_t sector, UShort_t strip, 
		    Double_t& x, Double_t& y, Double_t& z) const;

  // AliGeometry member functions 
  virtual void GetGlobal(const AliRecPoint* p, TVector3& pos, 
			 TMatrix& mat) const;
  virtual void GetGlobal(const AliRecPoint* p, TVector3& pos) const;
  virtual Bool_t Impact(const TParticle* particle) const;
protected:
  Bool_t        fIsInitialized;
  AliFMDRing*	fInner;		// Inner ring geometry information
  AliFMDRing*	fOuter;		// Outer ring geometry information
  AliFMD1*	fFMD1;		// FMD1 geometry information
  AliFMD2*	fFMD2;		// FMD2 geometry information
  AliFMD3*	fFMD3;		// FMD3 geometry information
  Bool_t	fUseFMD1;	// Wheter to Use FMD1 or not
  Bool_t	fUseFMD2;	// Wheter to Use FMD2 or not
  Bool_t	fUseFMD3;	// Wheter to Use FMD3 or not
  static AliFMDGeometry* fgInstance;
  AliFMDGeometry();
  virtual ~AliFMDGeometry() {}

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
