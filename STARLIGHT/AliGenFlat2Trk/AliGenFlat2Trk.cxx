// -*- C++ -*-
// $Id$

#include <Riostream.h>
#include <TFile.h>
#include <TMath.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>

#include <TF1.h>

#include <TVector2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

#include "AliGenEventHeader.h"
#include "AliRun.h"
#include "AliLog.h"

#include "AliGenFlat2Trk.h"

ClassImp(AliGenFlat2Trk);

//------------------------------------------------------------

AliGenFlat2Trk::AliGenFlat2Trk()
  : fMass(0)
  , fMinvMin(0)
  , fMinvMax(0)
  , fPtPairMin(0)
  , fPtPairMax(0)
  , fYPairMin(0)
  , fYPairMax(0)
  , fTries(0)
  , fRunAsMCGenerator(kTRUE)
  , fDebug(0)
  , fEvent(0)
  , fNpartProd(0)
  , fTheta(NULL)
  , fHeader(NULL) {
  fPid[0] = fPid[1] = 0;
}

AliGenFlat2Trk::AliGenFlat2Trk(Int_t    pid,
			       Double_t mInvMax,
			       Double_t ptPairMin,
			       Double_t ptPairMax,
			       Double_t yPairMin,
			       Double_t yPairMax,
			       Int_t    tries,
			       Bool_t   runAsMCGenerator,
			       TString  thetaFunction)
  : fMass(0)
  , fMinvMin(0) // default lower bound is at threshold
  , fMinvMax(mInvMax)
  , fPtPairMin(ptPairMin)
  , fPtPairMax(ptPairMax)
  , fYPairMin(yPairMin)
  , fYPairMax(yPairMax)
  , fTries(tries)
  , fRunAsMCGenerator(runAsMCGenerator)
  , fDebug(0)
  , fEvent(0)
  , fNpartProd(0)
  , fTheta(NULL)
  , fHeader(NULL)
{
  SetPid(pid);

  // Avoid zero pt
  if (fPtMin == 0.0) fPtMin = 1.E-04;

  // decay angle distribution (default: spin1 -> pi+pi-)
  fTheta = new TF1("fTheta",
		   (thetaFunction == ""
		    ? "sin(x)*(1-cos(x)*cos(x))"
		    : thetaFunction),
		   0.0, TMath::Pi());
}

AliGenFlat2Trk::~AliGenFlat2Trk()
{
  SafeDelete(fTheta);
  SafeDelete(fHeader);
}

Double_t AliGenFlat2Trk::SetMinvMin(Double_t m) {
  fMinvMin = TMath::Max(2*fMass, m);
  return fMinvMin;
}

void AliGenFlat2Trk::SetPid(Int_t pid) {
  TParticlePDG *pPlus  = TDatabasePDG::Instance()->GetParticle(+pid);
  TParticlePDG *pMinus = TDatabasePDG::Instance()->GetParticle(-pid);
  fPid[0] = pid;
  fPid[1] = (pMinus ? -pid : pid);
  if (!pPlus)
    AliFatalF("particle with pdg code=%d not found", pid);

  fMass    = TDatabasePDG::Instance()->GetParticle(fPid[0])->Mass();
  fMinvMin = TMath::Max(2*fMass, fMinvMin);
}

void AliGenFlat2Trk::Init()
{
  // print configuration
  TParticlePDG *pPlus  = TDatabasePDG::Instance()->GetParticle(fPid[0]);
  TParticlePDG *pMinus = TDatabasePDG::Instance()->GetParticle(fPid[1]);
  AliInfoF("produced particle species: %d,%d (%s,%s)",
	   fPid[0], fPid[1], pPlus->GetName(), pMinus->GetName());

  TString form = (fTheta ? fTheta->GetExpFormula() : "theta function is not set");
  form.ReplaceAll("x", "theta");
  AliInfoF("decay angle function: '%s'", form.Data());

  AliInfoF("Minv     = [%6.3f,%6.3f]", fMinvMin,   fMinvMax);
  AliInfoF("pair-pT  = [%6.3f,%6.3f]", fPtPairMin, fPtPairMax);
  AliInfoF("pair-Y   = [%6.3f,%6.3f]", fYPairMin,  fYPairMax);

  if (!fTheta)
    AliFatal("NULL == fTheta");
}

AliGenFlat2Trk::PPair
AliGenFlat2Trk::GeneratePair(Double_t y,         // rapidity
			     Double_t pt,        // p_T
			     Double_t m) const { // minv
  const TVector2 vpt(TVector2(pt, 0.0).Rotate(gRandom->Uniform(0.0, TMath::TwoPi())));
  const Double_t mt(TMath::Sqrt(m*m + pt*pt));
  const TLorentzVector v(vpt.X(),
			 vpt.Y(),
			 mt*TMath::SinH(y),
			 mt*TMath::CosH(y));

  TVector3 vp3;
  vp3.SetMagThetaPhi(TMath::Sqrt(0.25*m*m - fMass*fMass),
		     fTheta->GetRandom(),
		     TMath::TwoPi() * gRandom->Rndm());

  TLorentzVector vp4_A( vp3, 0.5*m);
  TLorentzVector vp4_B(-vp3, 0.5*m);

  return PPair(vp4_A.Transform(TLorentzRotation(v.BoostVector())),
	       vp4_B.Transform(TLorentzRotation(v.BoostVector())),
	       1.);
}

AliGenFlat2Trk::PPair
AliGenFlat2Trk::GeneratePair() const {
  const Double_t y  = gRandom->Uniform(fYPairMin,  fYPairMax);
  const Double_t pt = gRandom->Uniform(fPtPairMin, fPtPairMax);
  const Double_t m  = gRandom->Uniform(fMinvMin,   fMinvMax);

  PPair ppResult;
  Int_t counter=0;
  for (Int_t i=0; i<fTries; ++i) {
    const PPair pp(GeneratePair(y, pt, m));
    if (   pp.fV0.Rapidity() > fYMin && pp.fV0.Rapidity() < fYMax
	&& pp.fV1.Rapidity() > fYMin && pp.fV1.Rapidity() < fYMax
	&& pp.fV0.Perp() > fPtMin && pp.fV0.Perp() < fPtMax
	&& pp.fV1.Perp() > fPtMin && pp.fV1.Perp() < fPtMax) {
      ppResult=pp;
      ++counter;
    }
  }
  ppResult.fWeight = Double_t(counter)/Double_t(fTries);

  return ppResult;
}

void AliGenFlat2Trk::Generate()
{
  Float_t random[6] = { 0.,0.,0., 0.,0.,0. };
  Float_t origin[3] = { 0.,0.,0. };

  for (Int_t j=0; j<3; ++j)
    origin[j] = fOrigin[j];

  if (fVertexSmear == kPerEvent) {
    Rndm(random, 6);
    for (Int_t j=0; j<3; ++j) {
      origin[j] +=
	fOsigma[j]
	*TMath::Cos(2*random[2*j]*TMath::Pi())
	*TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
    }
  }

  PPair pp = GeneratePair();
  while (pp.fWeight == 0.0)
    pp = GeneratePair();

  AliInfoF("Generate (M,pT,Y)=(%6.3f,%6.3f,%6.3f) 1st,2nd: Y=(%5.2f,%5.2f), Pt=(%6.3f,%6.3f)",
	   pp.M(), pp.Perp(), pp.Y(),
	   pp.fV0.Rapidity(), pp.fV1.Rapidity(),
	   pp.fV0.Perp(),     pp.fV1.Perp());
  pp.fV0.Print();
  pp.fV1.Print();

  if (fRunAsMCGenerator) {
    Int_t ntr = 0;
    PushTrack(fTrackIt, -1, fPid[0],
	      pp.fV0.Px(), pp.fV0.Py(), pp.fV0.Pz(), pp.fV0.E(), // 4-momentun
	      origin[0], origin[1], origin[2], 0.0,      // vertex+tof
	      0.0, 0.0, 0.0,                             // polarization
	      kPPrimary, ntr,
	      pp.fWeight);
    ++fNpartProd;

    PushTrack(fTrackIt, -1, fPid[1],
	      pp.fV1.Px(), pp.fV1.Py(), pp.fV1.Pz(), pp.fV1.E(), // 4-momentun
	      origin[0], origin[1], origin[2], 0.0,      // vertex+tof
	      0.0, 0.0, 0.0,                             // polarization
	      kPPrimary, ntr,
	      pp.fWeight);
    ++fNpartProd;

    ++fEvent;

    // vertex
    const TArrayF eventVertex(3, origin);

    // Header
    SafeDelete(fHeader);
    fHeader = new AliGenEventHeader("AliGenFlat2Trk");

    // Event Vertex
    fHeader->SetPrimaryVertex(eventVertex);
    fHeader->SetNProduced(fNpartProd);
    AddHeader(fHeader);
  }
}
