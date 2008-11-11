#ifndef ALIFMDRING_H
#define ALIFMDRING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDRing.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:47:43 2006
    @brief   FMD ring geometry parameters 
*/
//__________________________________________________________________
//
// Parameters of the FMD rings. 
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

/** 
 * @brief Geometry description and parameters of a ring in the FMD 
 * detector. 
 *
 * As there are only 2 kinds of rings @e Inner (@c 'I') and @e Outer
 * (@c 'O') the two objects of this class is owned by the
 * AliFMDGeometry singleton object.  The 3 AliFMDDetector objects
 * shares these two instances as needed.
 * @ingroup FMD_base
*/
class AliFMDRing : public TNamed
{
public:
  /** 
   * CTOR
   *
   * @param fId Ring ID  
   */
  AliFMDRing(Char_t fId);
  /** 
   * DTOR  
   */
  virtual ~AliFMDRing() {}
  /** 
   * Initialize the ring geometry 
   */
  virtual void Init();
  
  /** 
   *
   * @param x Value of The Id of this ring type
   */
  void SetId(Char_t x) { fId = x; }
  /** 
   *
   * @param x Value of With of bonding pad on sensor
   */
  void SetBondingWidth(Double_t x=.5) { fBondingWidth = x; }
  /** 
   *
   * @param x Value of Size of wafer the sensor was made from
   */
  void SetWaferRadius(Double_t x=13.4/2) { fWaferRadius = x; }
  /** 
   *
   * @param x Value of Thickness of sensor
   */
  void SetSiThickness(Double_t x=.032) { fSiThickness = x; }
  /** 
   *
   * @param x Value of Lower radius of ring
   */
  void SetLowR(Double_t x) { fLowR = x; }
  /** 
   *
   * @param x Value of Upper radius of ring
   */
  void SetHighR(Double_t x) { fHighR = x; }
  /** 
   *
   * @param x Value of Opening angle of the silicon wafers
   */
  void SetTheta(Double_t x) { fTheta = x; }
  /** 
   *
   * @param x Value of Number of strips
   */
  void SetNStrips(Int_t x) { fNStrips = x; }
  /** 
   *
   * @param x Value of How far the ring extends beyond the z value given.
   */
  void SetRingDepth(Double_t x) { fRingDepth = x; }
  /** 
   *
   * @param x Value of Radius of support legs
   */
  void SetLegRadius(Double_t x=.25) { fLegRadius = x; }
  /** 
   *
   * @param x Value of Radius of support legs
   */
  void SetLegLength(Double_t x=.9) { fLegLength = x; }
  /** 
   *
   * @param x Value of Radius of support legs
   */
  void SetLegOffset(Double_t x=2) { fLegOffset = x; }
  /** 
   *
   * @param x Value of Staggering offset
   */
  void SetModuleSpacing(Double_t x=.5) { fModuleSpacing = x; }
  /** 
   *
   * @param x Value of Thickness of print board
   */
  void SetPrintboardThickness(Double_t x=.08) { fPrintboardThickness = x; }
  /** 
   *
   * @param x Value of Thickness of copper on print board
   */
  void SetCopperThickness(Double_t x=.01) { fCopperThickness = x; }
  /** 
   *
   * @param x Value of Thickness of chip on print board
   */
  void SetChipThickness(Double_t x=.01) { fChipThickness = x; }
  /** 
   *
   * @param x Value of spacing between si and PCB
   */
  void SetSpacing(Double_t x=.05) { fSpacing = x; }
  /** 
   *
   * @param x Thickness of honeycomb plate
   */
  void SetHoneycombThickness(Double_t x=0.65) { fHoneycombThickness = x; }
  /** 
   *
   * @param x Thickness of aluminium of honeycomb
   */
  void SetAlThickness(Double_t x=.1) { fAlThickness = x; }

  /** 
   *
   * @return The Id of this ring type
   */
  Char_t GetId() const { return fId; }
  /** 
   *
   * @return With of bonding pad on sensor
   */
  Double_t GetBondingWidth() const { return fBondingWidth; }
  /** 
   *
   * @return Size of wafer the sensor was made from
   */
  Double_t GetWaferRadius() const { return fWaferRadius; }
  /** 
   *
   * @return Thickness of sensor
   */
  Double_t GetSiThickness() const { return fSiThickness; }
  /** 
   *
   * @return Minimum r for an active strip
   */
  Double_t GetMinR() const { return fMinR; }
  /** 
   *
   * @return Maximum r for an active strip
   */
  Double_t GetMaxR() const { return fMaxR; }
  /** 
   *
   * @return Lower radius of ring
   */
  Double_t GetLowR() const { return fLowR; }
  /** 
   *
   * @return Upper radius of ring
   */
  Double_t GetHighR() const { return fHighR; }
  /** 
   *
   * @return Opening angle of the sector (half that of silicon wafers)
   */
  Double_t GetTheta() const { return fTheta; }
  /** 
   *
   * @return Number of strips
   */
  Int_t GetNStrips() const { return fNStrips; }
  /** 
   *
   * @return Number of sectors
   */
  Int_t GetNSectors() const { return Int_t(360. / fTheta); }
  /** 
   *
   * @return Number of modules (2 sectors per module)
   */
  Int_t GetNModules() const { return GetNSectors() / 2; }
  /** 
   *
   * @return How far the ring extends beyond the z value given.
   */
  Double_t GetRingDepth() const { return fRingDepth; }
  /** 
   *
   * @return Radius of support legs
   */
  Double_t GetLegRadius() const { return fLegRadius; }
  /** 
   *
   * @return Radius of support legs
   */
  Double_t GetLegLength() const { return fLegLength; }
  /** 
   *
   * @return Radius of support legs
   */
  Double_t GetLegOffset() const { return fLegOffset; }
  /** 
   *
   * @return Staggering offset
   */
  Double_t GetModuleSpacing() const { return fModuleSpacing; }
  /** 
   *
   * @return Thickness of print board
   */
  Double_t GetPrintboardThickness() const { return fPrintboardThickness; }
  /** 
   *
   * @return Thickness copper of print board
   */
  Double_t GetCopperThickness() const { return fCopperThickness; }
  /** 
   *
   * @return Thickness chip of print board
   */
  Double_t GetChipThickness() const { return fChipThickness; }
  /** 
   *
   * @return Value of spacing between si and PCB
   */
  Double_t GetSpacing() const { return fSpacing; }
  /** 
   *
   * @return Thickness of honeycomb plate
   */
  Double_t GetHoneycombThickness() const { return fHoneycombThickness; }
  /** 
   *
   * @return Thickness of aluminium of honeycomb
   */
  Double_t GetAlThickness() const { return fAlThickness; }
  /** 
   *
   * @return The strip pitch
   */ 
  Double_t GetPitch() const { return (fMaxR - fMinR) / fNStrips; }
  /** 
   *
   * @return Radius (in cm) correspondig to strip @a strip
   */
  Double_t GetStripRadius(UShort_t strip) const;
  /** 
   *
   * @return Full depth of (low) modules in this (half) ring
   */
  Double_t GetModuleDepth() const;
  /** 
   *
   * @return Full depth of this (half) ring
   */
  Double_t GetFullDepth() const;
  /** 
   * Get the inner radius of the digitizer cards 
   *
   * @return The inner radius of the digitizer cards 
   */
  Double_t GetFMDDLowR() const { return 1.2*GetLowR(); }
  /** 
   * Get the outer radius of the digitizer cards 
   *
   * @return The outer radius of the digitizer cards 
   */
  Double_t GetFMDDHighR() const { return .95*GetHighR(); }
  /** 
   *
   * @return Thickness of print board
   */
  Double_t GetFMDDPrintboardThickness() const { return 2*fPrintboardThickness; }
  /** 
   *
   * @return Thickness copper of print board
   */
  Double_t GetFMDDCopperThickness() const { return 2*fCopperThickness; }
  /** 
   *
   * @return Thickness chip of print board
   */
  Double_t GetFMDDChipThickness() const { return 2*fChipThickness; }

  /** 
   *
   * @return List of verticies
   */
  const TObjArray& GetVerticies() const { return fVerticies; }
  /** 
   *
   * @return Number of verticies
   */
  Int_t GetNVerticies() const { return fVerticies.GetEntries(); }
  /** 
   * @param i Vertex number 
   *
   * @return the ith vertex 
   */
  TVector2* GetVertex(Int_t i) const;
  /** 
   *
   * @return List of verticies
   */
  const TObjArray& GetSensorVerticies() const { return fSensorVerticies; }
  /** 
   * @param i Vertex number 
   *
   * @return the ith vertex 
   */
  TVector2* GetSensorVertex(Int_t i) const;
  /** 
   *
   * @return List of verticies
   */
  const TObjArray& GetHybridVerticies() const { return fHybridVerticies; }
  /** 
   * @param i Vertex number 
   *
   * @return the ith vertex 
   */
  TVector2* GetHybridVertex(Int_t i) const;

  /** 
   * Get a list of feet positions 
   * 
   * 
   * @return List TVector2 of feet positions on hybrid card
   */  
  const TObjArray& GetFeetPositions() const { return fFeetPositions; }
  /** 
   * Get the number of feet positions 
   * 
   * 
   * @return Number of feet positions 
   */
  Int_t GetNFeetPositions() const { return fFeetPositions.GetEntries(); }
  /** 
   * Get the @a i feet position
   * 
   * @param i Index
   * 
   * @return The foot position of stand-off @a i 
   */  
  TVector2* GetFootPosition(Int_t i) const;
    
     
  /** Not used */
  void Detector2XYZ(UShort_t sector, UShort_t strip, 
		    Double_t& x, Double_t& y, Double_t& z) const;
  /** Not used */
  Bool_t XYZ2Detector(Double_t x, Double_t y, Double_t z, 
		      UShort_t& sector, UShort_t& strip) const;
private: 
  Char_t	fId;			// The Id of this ring type
  Double_t	fBondingWidth;		// With of bonding pad on sensor
  Double_t	fWaferRadius;		// Size of wafer sensor was made from
  Double_t	fSiThickness;		// Thickness of sensor
  Double_t	fLowR;			// Lower radius of ring
  Double_t	fHighR;			// Upper radius of ring
  Double_t	fMinR;			// Lower radius of active strips
  Double_t	fMaxR;			// Upper radius of active strips
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
  Double_t	fHoneycombThickness;	// Thickness of honeycomb plate
  Double_t	fAlThickness;		// Thickness of aluminium of honeycomb
  
  TObjArray	fVerticies;		// List of active sensor verticies
  TObjArray     fSensorVerticies;       // List of physical sensor verticies
  TObjArray     fHybridVerticies;       // List of hybrid card verticies
  TObjArray     fFeetPositions;         // List of feet positions
  
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
