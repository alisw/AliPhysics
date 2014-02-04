#include <limits>

#include "TObject.h"
#include "TList.h"
#include "TMath.h"
#include "Math/Minimizer.h"

#include "AliLog.h"
#include "AliVTrdTracklet.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"

#include "AliTRDtrackOnline.h"

AliTRDgeometry *AliTRDtrackOnline::fgGeometry = new AliTRDgeometry();

AliTRDtrackOnline::AliTRDtrackOnline() :
  TObject(),
  fNTracklets(0),
  fTracklets(),
  fTrackParametrizations()
{

}


AliTRDtrackOnline::~AliTRDtrackOnline()
{

}


void AliTRDtrackOnline::AddTracklet(AliVTrdTracklet *trkl)
{
  if (fNTracklets == fgkMaxTracklets)
    return;
  else
    fTracklets[fNTracklets++] = trkl;
}


Bool_t AliTRDtrackOnline::Fit(ROOT::Math::Minimizer *minim)
{
  // fit all attached parametrizations

  Bool_t minSuccess = kFALSE;

  if (minim) {
    TIter param(&fTrackParametrizations);

    while (AliTRDtrackParametrization *par = (AliTRDtrackParametrization*) param()) {

      AliTRDtrackResiduals res(this, par);
      minim->Clear();
      minim->SetFunction(res);
      par->SetParams(minim);
      minSuccess = minim->Minimize();
      par->GetParams(minim);
    }
  }

  return minSuccess;
}


AliTRDtrackPosition AliTRDtrackOnline::ExtrapolateToLayer(Int_t /* layer */)
{
  Int_t maxLayer = -1;
  AliVTrdTracklet *trklBest = 0x0;
  for (Int_t iTracklet = fNTracklets-1; iTracklet > -1; iTracklet--) {
    AliVTrdTracklet *trkl = (AliVTrdTracklet*) fTracklets[iTracklet];
    if (trkl->GetDetector() % 6 >= maxLayer) {
      maxLayer = trkl->GetDetector() % 6;
      trklBest = trkl;
    }
  }
  if (trklBest)
    return AliTRDtrackPosition(trklBest->GetLocalY(), GetZ(trklBest));
  else {
    AliFatal("No tracklet in this track");
    return AliTRDtrackPosition(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
  }
}


void AliTRDtrackOnline::Print(Option_t* /* option */) const
{
  printf("track with %i tracklets:\n", GetNTracklets());
  for (Int_t iTracklet = 0; iTracklet < fNTracklets; iTracklet++) {
    printf("  0x%08x %i %4.1f %4.1f\n",
	   ((AliVTrdTracklet*) fTracklets[iTracklet])->GetTrackletWord(),
	   ((AliVTrdTracklet*) fTracklets[iTracklet])->GetDetector() % 6,
	   ((AliVTrdTracklet*) fTracklets[iTracklet])->GetLocalY(),
	   GetZ((AliVTrdTracklet*) fTracklets[iTracklet]));
  }
  TIter next(&fTrackParametrizations);
  while (AliTRDtrackParametrization *param = (AliTRDtrackParametrization*) next()) {
    param->Print();
  }
}

AliTRDtrackPosition::AliTRDtrackPosition(Float_t y, Float_t z, Float_t dy) :
  TObject(),
  fY(y),
  fZ(z),
  fDy(dy)
{

}

AliTRDtrackPosition::~AliTRDtrackPosition()
{

}

Float_t AliTRDtrackPosition::Distance(AliVTrdTracklet *trkl) const
{
  return TMath::Hypot(trkl->GetLocalY() - fY, AliTRDtrackOnline::GetZ(trkl) - fZ);
}


AliTRDtrackParametrization::AliTRDtrackParametrization(const char* name, const char* title) :
  TNamed(name, title),
  fFitGood(kFALSE)
{

}

AliTRDtrackParametrizationStraightLine::AliTRDtrackParametrizationStraightLine() :
  AliTRDtrackParametrization("straight line", "straight line"),
  fOffsetY(0),
  fSlopeY(0),
  fOffsetZ(0),
  fSlopeZ(0)
{

}

AliTRDtrackParametrizationStraightLine::AliTRDtrackParametrizationStraightLine(Double_t offsetY, Double_t slopeY,
									       Double_t offsetZ, Double_t slopeZ) :
  AliTRDtrackParametrization("straight line", Form("straight line: y = %4.2f + %4.2f * x, z = %4.2f + %4.2f *x",
						   offsetY, slopeY, offsetZ, slopeZ)),
  fOffsetY(offsetY),
  fSlopeY(slopeY),
  fOffsetZ(offsetZ),
  fSlopeZ(slopeZ)
{

}

void AliTRDtrackParametrizationStraightLine::SetParams(ROOT::Math::Minimizer * minim)
{
  minim->SetVariable(0, "offsety", 0., 0.1);
  minim->SetVariable(1, "slopey", 0., 0.1);
  // minim->SetVariable(2, "offsetz", 0., 0.1);
  minim->SetFixedVariable(2, "offsetz", 0.);
  minim->SetVariable(3, "slopez", 0., 0.1);
}

void AliTRDtrackParametrizationStraightLine::GetParams(ROOT::Math::Minimizer * minim)
{
  fOffsetY = minim->X()[0];
  fSlopeY  = minim->X()[1];
  fOffsetZ = minim->X()[2];
  fSlopeZ  = minim->X()[3];
}

void AliTRDtrackParametrizationStraightLine::SetValues(const Double_t *par)
{
  fOffsetY = par[0];
  fSlopeY  = par[1];
  fOffsetZ = par[2];
  fSlopeZ  = par[3];
}

AliTRDtrackPosition AliTRDtrackParametrizationStraightLine::ExtrapolateToLayer(Int_t layer)
{
  Float_t y = fOffsetY + fSlopeY * AliTRDtrackOnline::fgGeometry->GetTime0(layer);
  Float_t z = fOffsetZ + fSlopeZ * AliTRDtrackOnline::fgGeometry->GetTime0(layer);
  return AliTRDtrackPosition(y, z, fSlopeY*3.);
}

AliTRDtrackPosition AliTRDtrackParametrizationStraightLine::ExtrapolateToX(Float_t x)
{
  Float_t y = fOffsetY + fSlopeY * x;
  Float_t z = fOffsetZ + fSlopeZ * x;
  return AliTRDtrackPosition(y, z, fSlopeY*3.);
}

void AliTRDtrackParametrizationStraightLine::Print(Option_t * /* option */) const
{
  printf("straight line: offsetY = %4.1f, slopeY = %4.1f; offsetZ = %4.1f, slopeZ = %4.1f\n",
	 fOffsetY, fSlopeY, fOffsetZ, fSlopeZ);
}

AliTRDtrackParametrizationCurved::AliTRDtrackParametrizationCurved() :
  AliTRDtrackParametrization("helix", "helix"),
  fRadiusInv(0.),
  fOffsetY(0.),
  fOffsetZ(0.),
  fSlopeZ(0.),
  fOffsetX(300.)
{

}


void AliTRDtrackParametrizationCurved::SetParams(ROOT::Math::Minimizer * minim)
{
  minim->SetVariable(0, "offsety", 0., 0.1);
  minim->SetVariable(1, "invradius", 0., 0.1);
  // minim->SetVariable(2, "offsetz", 1., 0.1);
  minim->SetFixedVariable(2, "offsetz", 0.);
  minim->SetVariable(3, "slopez", 0., 0.1);
}


void AliTRDtrackParametrizationCurved::GetParams(ROOT::Math::Minimizer * minim)
{
  this->SetValues(minim->X());
}


void AliTRDtrackParametrizationCurved::SetValues(const Double_t *par)
{
  fOffsetY    = par[0];
  fRadiusInv  = par[1];
  fOffsetZ    = par[2];
  fSlopeZ     = par[3];
}

AliTRDtrackPosition AliTRDtrackParametrizationCurved::ExtrapolateToLayer(Int_t layer)
{
  return ExtrapolateToX(AliTRDtrackOnline::fgGeometry->GetTime0(layer));
}

AliTRDtrackPosition AliTRDtrackParametrizationCurved::ExtrapolateToX(Float_t x)
{
  Double_t yext1 = GetY(x);
  Double_t yext2 = GetY(x + 3.);

  Double_t zext = fOffsetZ + fSlopeZ * x;

  return AliTRDtrackPosition(yext1, zext, yext2-yext1);
}

Float_t AliTRDtrackParametrizationCurved::GetY(Float_t x)
{
 Double_t yext = 0.;
  // use Taylor expansion for small 1/R
  if (TMath::Abs(fRadiusInv) < 1.) {
    // offset
    yext  = fOffsetY * x/fOffsetX;
    // linear term
    yext += - (fOffsetX - x) * x * fRadiusInv /
      (2 * (fOffsetX*1./TMath::Sqrt(fOffsetX*fOffsetX + fOffsetY*fOffsetY)) *
       (fOffsetX*1./TMath::Sqrt(fOffsetX*fOffsetX + fOffsetY*fOffsetY)) *
       (fOffsetX*1./TMath::Sqrt(fOffsetX*fOffsetX + fOffsetY*fOffsetY)));
  }
  else {
    Double_t disc = 1./(fOffsetX*fOffsetX + fOffsetY*fOffsetY) - fRadiusInv*fRadiusInv/4.;
    if (disc < 0) {
      AliError("Discriminant < 0");
      return 1000.;
    }
    yext = TMath::Sqrt(disc) -
      TMath::Sqrt((fRadiusInv*fOffsetY/2. + fOffsetX * TMath::Sqrt(disc)) *
		  (fRadiusInv*fOffsetY/2. + fOffsetX * TMath::Sqrt(disc)) /
		  (fOffsetX*fOffsetX) -
		  fRadiusInv*fRadiusInv/(fOffsetX*fOffsetX)* x*x +
		  fRadiusInv*fRadiusInv/fOffsetX * x +
		  2 * fRadiusInv * fOffsetY * (x - fOffsetX)/(fOffsetX*fOffsetX) * TMath::Sqrt(disc));
    yext = fOffsetY/2. - fOffsetX * yext / fRadiusInv;
  }

  return yext;
}

void AliTRDtrackParametrizationCurved::Print(Option_t * /* option */) const
{
  printf("helix curve: 1/R = %f, y = %4.1f\n", fRadiusInv, fOffsetY);
}


AliTRDtrackResiduals::AliTRDtrackResiduals(const AliTRDtrackOnline *track, AliTRDtrackParametrization *param) :
  ROOT::Math::IBaseFunctionMultiDim(),
  fTrack(track),
  fParam(param)
{

}

AliTRDtrackResiduals::AliTRDtrackResiduals(const AliTRDtrackResiduals &rhs) :
  ROOT::Math::IBaseFunctionMultiDim(rhs),
  fTrack(rhs.fTrack),
  fParam(rhs.fParam)
{

}

AliTRDtrackResiduals& AliTRDtrackResiduals::operator=(const AliTRDtrackResiduals &rhs)
{
  if (&rhs != this) {
    ROOT::Math::IBaseFunctionMultiDim::operator=(rhs);
    fTrack = rhs.fTrack;
    fParam = rhs.fParam;
  }

  return *this;
}

AliTRDtrackResiduals* AliTRDtrackResiduals::Clone() const
{
  return new AliTRDtrackResiduals(*this);
}

Double_t AliTRDtrackResiduals::DoEval(const Double_t *par) const
{
  // calculate chi2 for the given values for the parametrization

  // initialisation
  Float_t deltaY = 0.;
  Float_t deltaZ = 0.;
  Float_t chi2 = 0.;

  // actually set the values for the parametrization
  fParam->SetValues(par);

  // loop over all contributing tracklets
  for (Int_t iTracklet = 0; iTracklet < fTrack->GetNTracklets(); iTracklet++) {
    AliVTrdTracklet *trkl = fTrack->GetTracklet(iTracklet);

    // Int_t layer = trkl->GetDetector() % 6;

    AliTRDtrackPosition pos = fParam->ExtrapolateToX(AliTRDtrackOnline::GetX(trkl));
    Float_t yext = pos.GetY();
    Float_t zext = pos.GetZ();

    AliTRDpadPlane *pp = fgGeometry->GetPadPlane(trkl->GetDetector());
    Float_t zlen = 0.5 * pp->GetRowSize(trkl->GetBinZ());
    Float_t zpad = pp->GetRowPos(trkl->GetBinZ()) - zlen;
    zpad = AliTRDtrackOnline::GetZ(trkl);
    Float_t zrel = zext - zpad;
    if (zrel > zlen)
      zrel = zlen;
    else if (zrel < -zlen)
      zrel = -zlen;

    Float_t ycorr = trkl->GetLocalY() + TMath::Tan(TMath::Pi()/180.*pp->GetTiltingAngle()) * zrel;

    deltaY = ycorr        - yext;
    deltaZ = AliTRDtrackOnline::GetZ(trkl) - zext;
    deltaY /= 0.3;
    deltaZ /= 3.;
//     printf("in layer %i: deltaY = %f, deltaZ = %f\n", layer, deltaY, deltaZ);

    chi2 += deltaY*deltaY + deltaZ*deltaZ;
  }

//   printf("chi2 = %f\n", chi2);
  return chi2;
}
