//
// $Id$
//
#ifndef ALIFMD3_H
#define ALIFMD3_H
/** @file    AliFMD3.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:23:48 2006
    @brief   Geometry parameters of the FMD3 detector. 
*/
#ifndef ALIFMDDETECTOR_H
# include "AliFMDDetector.h"
#endif

/** @class AliFMD3 AliFMD3.h <FMD/AliFMD3.h> 
    @brief Geometry parameters of the FMD3 detector. 
    FMD3 has a fairly complicated support structure.  The cone also
    supports the beam-pipe.
    @image html FMD3.png 
    @ingroup FMD_base
*/
class AliFMD3 : public AliFMDDetector 
{
public: 
  /** Constructor 
      @param inner Pointer to inner ring description 
      @param outer Pointer to outer ring description */
  AliFMD3(AliFMDRing* inner, AliFMDRing* outer);
  /** Destructor */
  virtual ~AliFMD3(){}

  /** Initialize the geometry */
  virtual void Init();

  /** @param z Z position of front of nose */
  void SetNoseZ(Double_t z=-46) { fNoseZ = z; }
  /** @param r Nose inner radius */
  void SetNoseLowR(Double_t r=5.5) { fNoseLowR = r; }
  /** @param r Nose outer radius */
  void SetNoseHighR(Double_t r=6.7) { fNoseHighR = r; }
  /** @param l Length of nose in Z */
  void SetNoseLength(Double_t l=2.8) { fNoseLength = l; }
  /** @param r Inner radius of base of cone */
  void SetBackLowR(Double_t r=61./2) { fBackLowR = r; }
  /** @param r Outer radius of base of cone */
  void SetBackHighR(Double_t r=66.8/2) { fBackHighR = r; }
  /** @param l Length of base of cone in Z */
  void SetBackLength(Double_t l=1.4) { fBackLength = l; }
  /** @param t Thickness of support beams */
  void SetBeamThickness(Double_t t=.5) { fBeamThickness = t; }
  /** @param w Width of support beams */
  void SetBeamWidth(Double_t w=6) { fBeamWidth = w; }
  /** @param l Length of the cone in Z */
  void SetConeLength(Double_t l=30.9) { fConeLength = l; }
  /** @param r Outer radius of flanges */
  void SetFlangeR(Double_t r=49.25) { fFlangeR = r; }
  /** @param n Number of support beams */
  void SetNBeam(Int_t n=8) { fNBeam = n; }
  /** @param n Number of support flanges */
  void SetNFlange(Int_t n=4) { fNFlange = n; }

  /** @return Z position of front of nose */
  Double_t GetNoseZ() const { return fNoseZ; }
  /** @return Nose inner radius */
  Double_t GetNoseLowR() const { return fNoseLowR; }
  /** @return Nose outer radius */
  Double_t GetNoseHighR() const { return fNoseHighR; }
  /** @return Length of nose in Z */
  Double_t GetNoseLength() const { return fNoseLength; }
  /** @return Inner radius of base of cone */
  Double_t GetBackLowR() const { return fBackLowR; }
  /** @return Outer radius of base of cone */
  Double_t GetBackHighR() const { return fBackHighR; }
  /** @return Length of base of cone in Z */
  Double_t GetBackLength() const { return fBackLength; }
  /** @return Thickness of support beams */
  Double_t GetBeamThickness() const { return fBeamThickness; }
  /** @return Width of support beams */
  Double_t GetBeamWidth() const { return fBeamWidth; }
  /** @return Length of the cone in Z */
  Double_t GetConeLength() const { return fConeLength; }
  /** @return Outer radius of flanges */
  Double_t GetFlangeR() const { return fFlangeR; }
  /** @return Midpoint of mother volume */
  Double_t GetZ() const { return fZ; }
  /** @return Slope of cone */
  Double_t GetAlpha() const { return fAlpha; }
  /** @return Number of support beams */
  Int_t GetNBeam() const { return fNBeam; }
  /** @return Number of support flanges */
  Int_t GetNFlange() const { return fNFlange; }

  /** Get the cone radii at @a z. 
      @param z Point to evaulate at 
      @param opt If @c "O" get the outer radii, if @c "I" get the
      inner radii. 
      @return the radius of the cone */
  Double_t ConeR(Double_t z, Option_t* opt="O") const;

protected: 
  Double_t	fNoseZ;			// Z position of front of nose
  Double_t	fNoseLowR;		// Nose inner radius
  Double_t	fNoseHighR;		// Nose outer radius
  Double_t	fNoseLength;		// Length of nose in Z
  Double_t	fBackLowR;		// Inner radius of base of cone
  Double_t	fBackHighR;		// Outer radius of base of cone
  Double_t	fBackLength;		// Length of base of cone in Z
  Double_t	fBeamThickness;		// Thickness of support beams
  Double_t	fBeamWidth;		// Width of support beams
  Double_t	fConeLength;		// Length of the cone in Z
  Double_t	fFlangeR;		// Outer radius of flanges
  Double_t	fZ;			// Midpoint of mother volume
  Double_t	fAlpha;			// Slope of cone
  Int_t		fNBeam;			// Number of support beams
  Int_t		fNFlange;		// Number of support flangesy
  ClassDef(AliFMD3, 1);
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
//
// EOF
//
