#ifndef ALITRDTRACKONLINE_H
#define ALITRDTRACKONLINE_H

#include "TObject.h"
#include "TList.h"

#include "Math/IFunction.h"
#include "Math/Minimizer.h"

#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"

class AliVTrdTracklet;
class AliTRDtrackPosition;
class AliTRDtrackOnline;
class AliTRDgeometry;

class AliTRDtrackParametrization : public TNamed
{
 public:
  AliTRDtrackParametrization(const char* name = "", const char* title = "");
  ~AliTRDtrackParametrization() {}

  virtual void SetParams(ROOT::Math::Minimizer * minim) = 0;
  virtual void GetParams(ROOT::Math::Minimizer * minim) = 0;
  virtual void SetValues(const Double_t *par) = 0;
  virtual Int_t GetNDim() const = 0;

  virtual void UpdateTitle() {}

  virtual AliTRDtrackPosition ExtrapolateToLayer(Int_t layer) = 0;
  virtual AliTRDtrackPosition ExtrapolateToX(Float_t x) = 0;

  Bool_t IsFitGood() const { return fFitGood; }

 protected:
  Bool_t fFitGood;

  ClassDef(AliTRDtrackParametrization, 1);
};


class AliTRDtrackParametrizationStraightLine : public AliTRDtrackParametrization
{
 public:
  AliTRDtrackParametrizationStraightLine();
  AliTRDtrackParametrizationStraightLine(Double_t offsetY, Double_t slopeY,
					 Double_t offsetZ, Double_t slopeZ);

  virtual void SetParams(ROOT::Math::Minimizer * minim);
  virtual void GetParams(ROOT::Math::Minimizer * minim);
  virtual void SetValues(const Double_t *par);
  virtual Int_t GetNDim() const { return 4; }

  AliTRDtrackPosition ExtrapolateToLayer(Int_t layer);
  AliTRDtrackPosition ExtrapolateToX(Float_t x);

  Double_t GetOffsetY() const { return fOffsetY; }
  Double_t GetOffsetZ() const { return fOffsetZ; }
  Double_t GetSlopeY()  const { return fSlopeY; }
  Double_t GetSlopeZ()  const { return fSlopeZ; }

  void Print(Option_t *option = "") const;

 protected:
  Double_t fOffsetY;
  Double_t fSlopeY;
  Double_t fOffsetZ;
  Double_t fSlopeZ;

  ClassDef(AliTRDtrackParametrizationStraightLine, 1);
};


class AliTRDtrackParametrizationCurved : public AliTRDtrackParametrization
{
 public:
  AliTRDtrackParametrizationCurved();

  virtual void SetParams(ROOT::Math::Minimizer * minim);
  virtual void GetParams(ROOT::Math::Minimizer * minim);
  virtual void SetValues(const Double_t *par);
  virtual Int_t GetNDim() const { return 4; }

  AliTRDtrackPosition ExtrapolateToLayer(Int_t layer);
  AliTRDtrackPosition ExtrapolateToX(Float_t x);

  Float_t GetY(Float_t x);

  void Print(Option_t *option = "") const;

 protected:
  // parameters
  Double_t fRadiusInv;
  Double_t fOffsetY;
  Double_t fOffsetZ;
  Double_t fSlopeZ;

  // fixed values
  Double_t fOffsetX;

  ClassDef(AliTRDtrackParametrizationCurved, 1);
};


class AliTRDtrackOnline : public TObject
{
 public:
  AliTRDtrackOnline();
  ~AliTRDtrackOnline();

  void AddTracklet(AliVTrdTracklet *trkl);

  Bool_t Fit(ROOT::Math::Minimizer *minim);

  Int_t GetNTracklets() const { return fNTracklets; }
  AliVTrdTracklet* GetTracklet(Int_t i) const { return i < fNTracklets ? (AliVTrdTracklet*) fTracklets[i] : 0x0; }

  AliTRDtrackPosition ExtrapolateToLayer(Int_t layer);

  void AddParametrization(AliTRDtrackParametrization *param) { fTrackParametrizations.Add(param); }
  const TList& GetParametrizations() const { return fTrackParametrizations; }

  void Print(Option_t *option = "") const;

  static Float_t GetX(AliVTrdTracklet *trkl) { return fgGeometry->GetTime0(trkl->GetDetector() % 6); }
  static Float_t GetZ(AliVTrdTracklet *trkl) { return fgGeometry->GetPadPlane((trkl->GetDetector() % 6), (trkl->GetDetector()/6) % 5)->GetRowPos(trkl->GetBinZ()) -
      fgGeometry->GetPadPlane((trkl->GetDetector() % 6), (trkl->GetDetector()/6) % 5)->GetRowSize(trkl->GetBinZ()); }
  static AliTRDgeometry *fgGeometry;

 protected:
  static const Int_t fgkMaxTracklets = 10;

  Int_t fNTracklets;
  TObjArray fTracklets;

  TList fTrackParametrizations;

  ClassDef(AliTRDtrackOnline, 1);
};


class AliTRDtrackPosition : public TObject
{
 public:
  AliTRDtrackPosition(Float_t y, Float_t z, Float_t dy = 0.);
  ~AliTRDtrackPosition();

  Float_t GetY()  const { return fY; }
  Float_t GetZ()  const { return fZ; }
  Float_t GetdY() const { return fDy; }

  Float_t Distance(AliVTrdTracklet *trkl) const;

 protected:
  Float_t fY;
  Float_t fZ;
  Float_t fDy;

  ClassDef(AliTRDtrackPosition, 1);
};


class AliTRDtrackResiduals : public ROOT::Math::IBaseFunctionMultiDim
{
public:
  AliTRDtrackResiduals(const AliTRDtrackOnline *track, AliTRDtrackParametrization *param);
  AliTRDtrackResiduals(const AliTRDtrackResiduals &rhs);
  AliTRDtrackResiduals& operator=(const AliTRDtrackResiduals &rhs);
  ~AliTRDtrackResiduals() {}

  AliTRDtrackResiduals* Clone() const;
  UInt_t NDim() const { return fParam->GetNDim(); }
  Double_t DoEval(const Double_t *par) const;

protected:
  const AliTRDtrackOnline *fTrack; // reference to track being fitted
  AliTRDtrackParametrization *fParam; // reference to the used parametrization

  static AliTRDgeometry *fgGeometry;
};

#endif
