//
// $Id$
//
#ifndef ALIFMD3SUPPORT_H
#define ALIFMD3SUPPORT_H

#ifndef ROOT_TObject
# include <TObject.h>
#endif
#ifndef ROOT_TArrayI
# include <TArrayI.h>
#endif

class AliFMD3Support : public TObject
{
public:
  AliFMD3Support();
  virtual ~AliFMD3Support();
  virtual void   SetupGeometry(Int_t airId, Int_t cId,
			       Double_t innerZl, Double_t innerZh, 
			       Double_t innerRl, Double_t outerZl, 
			       Double_t outerZh, Double_t outerRl);
  virtual void   Geometry(const char* mother, Int_t idRotId, Double_t z=0);
  virtual void   Gsatt() const;
  
  void SetNoseZ(Double_t         x=-46)		{ fNoseZ = x; }
  void SetNoseLowR(Double_t      x=5.5)		{ fNoseLowR = x; }
  void SetNoseHighR(Double_t     x=6.7)		{ fNoseHighR = x; }
  void SetNoseLength(Double_t    x=2.8)		{ fNoseLength = x; }
  void SetBackLowR(Double_t      x=61./2)	{ fBackLowR = x; }
  void SetBackHighR(Double_t     x=66.8/2)	{ fBackHighR = x; }
  void SetBackLength(Double_t    x=1.4)		{ fBackLength = x; }
  void SetBeamThickness(Double_t x=.5)		{ fBeamThickness = x; }
  void SetBeamWidth(Double_t     x=6)		{ fBeamWidth = x; }
  void SetConeLength(Double_t    x=30.9)	{ fConeLength = x; }
  void SetFlangeR(Double_t       x=49.25)	{ fFlangeR = x; }
  void SetNBeam(Int_t            n=8)           { fNBeam = n; }
  void SetNFlange(Int_t          n=4)           { fNFlange = n; }
  
  Double_t GetNoseZ() const		{ return fNoseZ; }
  Double_t GetNoseLowR() const		{ return fNoseLowR; }
  Double_t GetNoseHighR() const		{ return fNoseHighR; }
  Double_t GetNoseLength() const	{ return fNoseLength; }
  Double_t GetBackLowR() const		{ return fBackLowR; }
  Double_t GetBackHighR() const		{ return fBackHighR; }
  Double_t GetBackLength() const	{ return fBackLength; }
  Double_t GetBeamThickness() const	{ return fBeamThickness; }
  Double_t GetBeamWidth() const		{ return fBeamWidth; }
  Double_t GetConeLength() const	{ return fConeLength; }
  Double_t GetFlangeR() const	        { return fFlangeR; }
  Int_t    GetNBeam() const             { return fNBeam; }
  Int_t    GetNFlange() const           { return fNFlange; }
  Double_t GetZ() const                 { return fZ; }
  Double_t ConeR(Double_t z, Option_t* opt="O") const;
protected:
  Double_t fNoseZ;		// Z position of front of nose
  Double_t fNoseLowR;		// Nose inner radius
  Double_t fNoseHighR;		// Nose outer radius
  Double_t fNoseLength;		// Length of nose in Z
  Double_t fBackLowR;		// Inner radius of base of cone
  Double_t fBackHighR;		// Outer radius of base of cone
  Double_t fBackLength;		// Length of base of cone in Z
  Double_t fBeamThickness;	// Thickness of support beams
  Double_t fBeamWidth;		// Width of support beams
  Double_t fConeLength;		// Length of the cone in Z
  Double_t fFlangeR;            // Outer radius of flanges
  Double_t fZ;                  // Midpoint of mother volume 
  Double_t fAlpha;              // Slope of cone 
  
  Int_t    fNBeam;              // Number of support beams 
  Int_t    fNFlange;            // Number of support flanges
  Int_t    fNoseId;             // Id of nose volume
  Int_t    fBeamId;             // Id of beam volumes 
  Int_t    fBackId;             // Id of base volume
  Int_t    fFlangeId;           // Id of flange volume
  TArrayI  fRotations;
  
  static const Char_t* fgkNoseName; // Name of nose volume 
  static const Char_t* fgkBeamName; // Name of beam volumes
  static const Char_t* fgkBackName; // Name of base volume 
  static const Char_t* fgkFlangeName; // Name of flange volume 
  
  ClassDef(AliFMD3Support,1); // Geometry of Support for the FMD3 
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
