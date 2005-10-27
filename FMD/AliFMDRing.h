#ifndef ALIFMDRING_H
#define ALIFMDRING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
//__________________________________________________________________
//
// Parameters of the FMD rings. 
// 
// This class is responsible to make the (common) rings of the three
// sub-detectors. 
//
#ifndef ROOT_TNamed
# include <TNamed.h>
#endif
#ifndef ROOT_TObjArray
# include <TObjArray.h>
#endif

class TBrowser;
class TVector2;

/** Geometry description and parameters of a ring in the FMD
    detector. 
    
    As there are only 2 kinds of rings @e Inner (@c 'I') and @e
    Outer (@c 'O') the two objects of this class is owned by the
    Geometry::FMD singleton object.  The 3 Geometry::FMDDetector
    objects shares these two instances as needed. 
*/
class AliFMDRing : public TNamed
{
public:
  AliFMDRing(Char_t fId);
  virtual ~AliFMDRing() {}
  /** Initialize the ring geometry */
  virtual void Init();
  
  /** @param x Value of The Id of this ring type */
  void SetId(Char_t x) { fId = x; }
  /** @param x Value of With of bonding pad on sensor */
  void SetBondingWidth(Double_t x=.5) { fBondingWidth = x; }
  /** @param x Value of Size of wafer the sensor was made from */
  void SetWaferRadius(Double_t x=13.4/2) { fWaferRadius = x; }
  /** @param x Value of Thickness of sensor */
  void SetSiThickness(Double_t x=.03) { fSiThickness = x; }
  /** @param x Value of Lower radius of ring */
  void SetLowR(Double_t x) { fLowR = x; }
  /** @param x Value of Upper radius of ring */
  void SetHighR(Double_t x) { fHighR = x; }
  /** @param x Value of Opening angle of the silicon wafers */
  void SetTheta(Double_t x) { fTheta = x; }
  /** @param x Value of Number of strips */
  void SetNStrips(Int_t x) { fNStrips = x; }
  /** @param x Value of How far the ring extends beyond the z value given. */
  void SetRingDepth(Double_t x) { fRingDepth = x; }
  /** @param x Value of Radius of support legs */
  void SetLegRadius(Double_t x=.5) { fLegRadius = x; }
  /** @param x Value of Radius of support legs */
  void SetLegLength(Double_t x=1) { fLegLength = x; }
  /** @param x Value of Radius of support legs */
  void SetLegOffset(Double_t x=2) { fLegOffset = x; }
  /** @param x Value of Staggering offset */
  void SetModuleSpacing(Double_t x=1) { fModuleSpacing = x; }
  /** @param x Value of Thickness of print board */
  void SetPrintboardThickness(Double_t x=.1) { fPrintboardThickness = x; }
  /** @param x Value of Thickness of copper on print board */
  void SetCopperThickness(Double_t x=.01) { fCopperThickness = x; }
  /** @param x Value of Thickness of chip on print board */
  void SetChipThickness(Double_t x=.01) { fChipThickness = x; }
  /** @param x Value of spacing between si and PCB */
  void SetSpacing(Double_t x=.05) { fSpacing = x; }

  /** @return The Id of this ring type */
  Char_t GetId() const { return fId; }
  /** @return With of bonding pad on sensor */
  Double_t GetBondingWidth() const { return fBondingWidth; }
  /** @return Size of wafer the sensor was made from */
  Double_t GetWaferRadius() const { return fWaferRadius; }
  /** @return Thickness of sensor */
  Double_t GetSiThickness() const { return fSiThickness; }
  /** @return Lower radius of ring */
  Double_t GetLowR() const { return fLowR; }
  /** @return Upper radius of ring */
  Double_t GetHighR() const { return fHighR; }
  /** @return Opening angle of the sector (half that of silicon wafers) */
  Double_t GetTheta() const { return fTheta; }
  /** @return Number of strips */
  Int_t GetNStrips() const { return fNStrips; }
  /** @return Number of sectors */
  Int_t GetNSectors() const { return Int_t(360. / fTheta); }
  /** @return Number of modules (2 sectors per module) */
  Int_t GetNModules() const { return GetNSectors() / 2; }
  /** @return How far the ring extends beyond the z value given. */
  Double_t GetRingDepth() const { return fRingDepth; }
  /** @return Radius of support legs */
  Double_t GetLegRadius() const { return fLegRadius; }
  /** @return Radius of support legs */
  Double_t GetLegLength() const { return fLegLength; }
  /** @return Radius of support legs */
  Double_t GetLegOffset() const { return fLegOffset; }
  /** @return Staggering offset */
  Double_t GetModuleSpacing() const { return fModuleSpacing; }
  /** @return Thickness of print board */
  Double_t GetPrintboardThickness() const { return fPrintboardThickness; }
  /** @return Thickness copper of print board */
  Double_t GetCopperThickness() const { return fCopperThickness; }
  /** @return Thickness chip of print board */
  Double_t GetChipThickness() const { return fChipThickness; }
  /** @return Value of spacing between si and PCB */
  Double_t GetSpacing() const { return fSpacing; }

  /** @return List of verticies */
  const TObjArray& GetVerticies() const { return fVerticies; }
  /** @return Number of verticies */
  Int_t GetNVerticies() const { return fVerticies.GetEntries(); }
  /** @param i Vertex number 
      @return the ith vertex */
  TVector2* GetVertex(Int_t i) const;
     
  void Detector2XYZ(UShort_t sector, UShort_t strip, 
		    Double_t& x, Double_t& y, Double_t& z) const;
  
private: 
  Char_t	fId;			// The Id of this ring type
  Double_t	fBondingWidth;		// With of bonding pad on sensor
  Double_t	fWaferRadius;		// Size of wafer sensor was made from
  Double_t	fSiThickness;		// Thickness of sensor
  Double_t	fLowR;			// Lower radius of ring
  Double_t	fHighR;			// Upper radius of ring
  Double_t	fTheta;			// Opening angle of the silicon wafers
  Int_t		fNStrips;		// Number of strips
  Double_t	fRingDepth;		// How far the ring extends beyond z
  Double_t	fLegRadius;		// Radius of support legs
  Double_t 	fLegLength;		// Radius of support legs
  Double_t	fLegOffset;		// Radius of support legs
  Double_t	fModuleSpacing;		// Staggering offset
  Double_t	fPrintboardThickness;	// Thickness of print board
  Double_t	fCopperThickness;	// Thickness of Cu on print board
  Double_t	fChipThickness;		// Thickness of chip on print board
  Double_t      fSpacing;               // Spacing between si and PCB
  
  TObjArray	fVerticies;		// List of verticies

  ClassDef(AliFMDRing, 0);
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
