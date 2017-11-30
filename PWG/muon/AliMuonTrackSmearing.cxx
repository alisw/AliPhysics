/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliMuonTrackSmearing.h"

#include <math.h>

#include "TF1.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TDatabasePDG.h"

#include "AliLog.h"
//#include "AliVParticle.h"
//#include "AliMCEvent.h"
//#include "AliAODTrack.h"
//#include "AliPID.h"

/// \cond CLASSIMP
ClassImp(AliMuonTrackSmearing) // Class implementation in ROOT context
/// \endcond

//________________________________________________________________________
AliMuonTrackSmearing::AliMuonTrackSmearing ( TRootIOCtor* /*ioCtor*/ ) :
TObject(),
fChosenFunc(0),
fSigmaTrk(0.),
fSigmaTrkCut(0.),
fSigmaxChSt1(0.),
fSigmayChSt1(0.),
fSigmayCh(0.),
fNSigmaShift(0.),
fZB02(0.),
fZB23(0.),
fZB310(0.),
fMuMass(0.),
fTuneKalman(0),
fCrystalBallTails(),
fCrystalBall(NULL),
fRecoCharge(0.),
fRAbs(-1.),
fRecoTrack()
//,
//fRecoTrackList()
{
  /// Ctr.
}


//________________________________________________________________________
AliMuonTrackSmearing::AliMuonTrackSmearing ( Int_t chosenFunc ) :
TObject(),
fChosenFunc(chosenFunc),
fSigmaTrk(0.002),
fSigmaTrkCut(4.),
fSigmaxChSt1(0.),
fSigmayChSt1(0.),
fSigmayCh(0.),
fNSigmaShift(0.),
fZB02(3.82),
fZB23(4.72),
fZB310(4.46),
fMuMass(TDatabasePDG::Instance()->GetParticle("mu-")->Mass()),
fTuneKalman(kFALSE),
fCrystalBallTails(),
fCrystalBall(NULL),
fRecoCharge(0.),
fRAbs(-1.),
fRecoTrack()
//,
//fRecoTrackList()
{
  /// Ctr.
  SetupDefaultValues();
}

//________________________________________________________________________
AliMuonTrackSmearing::~AliMuonTrackSmearing()
{
  /// Dtr
  delete fCrystalBall;
//  ClearRecoTrackList();
}

//________________________________________________________________________
void AliMuonTrackSmearing::Print ( Option_t* /*option*/ ) const
{
  /// Print information

  TString funcNames[3] = {"CrystalBall","BreitWigner","Gaus"};
  printf("\n*** Muon track smearing parameters ***\n");
  printf("Chosen function: %s\n",funcNames[fChosenFunc].Data());
  printf("Sigma x cluster station 1: %g\n",fSigmaxChSt1);
  if ( fChosenFunc == kCrystalBall ) printf("  tail param: %g  %g\n",fCrystalBallTails[0],fCrystalBallTails[1]);
  printf("Sigma y cluster station 1: %g\n",fSigmayChSt1);
  if ( fChosenFunc == kCrystalBall ) printf("  tail param: %g  %g\n",fCrystalBallTails[2],fCrystalBallTails[3]);
  printf("Sigma y slope: %g\n",fSigmayCh);
  if ( fChosenFunc == kCrystalBall ) printf("  tail param: %g  %g\n",fCrystalBallTails[4],fCrystalBallTails[5]);
  printf("Sigma for tracker: %g\n",fSigmaTrk);
  printf("Sigma cut for tracker: %g\n",fSigmaTrkCut);
  Double_t sigmaThetaRes = TMath::Sqrt(SigmaThetaDevFromRes2());
  printf("Sigma systematic shift %g * %g = %g\n",fNSigmaShift,sigmaThetaRes,fNSigmaShift*sigmaThetaRes);
  printf("*****************************\n");
}



//________________________________________________________________________
void AliMuonTrackSmearing::SetupDefaultValues ( )
{
  /// Setup default values for parameters
  switch ( fChosenFunc ) {
    case kBreitWigner:
      fSigmaxChSt1 = 0.000550;
      fSigmayChSt1 = 0.000190;
      fSigmayCh = 0.000280;
      break;
    case kCrystalBall:
      fSigmaxChSt1 = 0.000401;
      fSigmayChSt1 = 0.000110;
      fSigmayCh = 0.000153;
      SetCrystalBallParams ( 2.017293, 1.890778, 1.514588, 1.938707, 1.280116, 2.239019 );
      break;
    default:
      fSigmaxChSt1 = 0.000419;
      fSigmayChSt1 = 0.000130;
      fSigmayCh = 0.000207;
  }
}

//________________________________________________________________________
void AliMuonTrackSmearing::SetCrystalBallParams ( Double_t xChSt1Par1, Double_t xChSt1Par2,
                                                   Double_t yChSt1Par1, Double_t yChSt1Par2,
                                                   Double_t yChPar1, Double_t yChPar2 )
{
  /// Setup tail parameters of CrystalBall
  fCrystalBallTails.clear();
  fCrystalBallTails.push_back(xChSt1Par1);
  fCrystalBallTails.push_back(xChSt1Par2);
  fCrystalBallTails.push_back(yChSt1Par1);
  fCrystalBallTails.push_back(yChSt1Par2);
  fCrystalBallTails.push_back(yChPar1);
  fCrystalBallTails.push_back(yChPar2);
}

//________________________________________________________________________
//void AliMuonTrackSmearing::ClearRecoTrackList ()
//{
//  /// Clear list of reconstructed tracks
//  /// This needs to be called at the end of each event to avoid memory leaks
//  /// FIXME: in principle one could simply use unique_ptr, defined in C++11
//  /// but the interpreter of CINT has some problems with it
//  for ( AliVParticle* obj : fRecoTrackList ) delete obj;
//  fRecoTrackList.clear();
//}

//________________________________________________________________________
Double_t AliMuonTrackSmearing::CrystalBallSymmetric ( Double_t *xx,Double_t *par )
{
  /// Crystal Ball definition

  ///par[0] = Normalization
  ///par[1] = mean
  ///par[2] = sigma
  ///par[3] = alpha = alpha'
  ///par[4] = n = n'

  Double_t tp = fabs((xx[0]-par[1])/par[2]);

  Double_t absAlpha = fabs(par[3]);
  Double_t ap = pow(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
  Double_t bp = par[4]/absAlpha - absAlpha;

  if (tp < absAlpha) return par[0]*(exp(-0.5*tp*tp)); // gaussian core
  else return par[0]*(ap/pow(bp + tp, par[4])); //left and right tails
}

//________________________________________________________________________
void AliMuonTrackSmearing::ComputeRecoTrack ( Double_t pGen, Double_t etaGen,
                                                Double_t phiGen, Double_t chargeGen )
{
  /// Given the generated MC track parameters,
  /// compute the track with smeared parameters according to resolution

  fRecoCharge = chargeGen;

  Double_t sign = ( chargeGen > 0. ) ? 1. : -1.;

  Double_t eLoss02 = ELoss(pGen,1.5); // = 110cm eLoss_Common + 305 cm eLoss_tungsten
  Double_t eLoss23 = ELoss(pGen,2.5); // = 378cm eLoss_Common + 37 cm eLoss_tungsten
  Double_t eLoss310 = ELoss(pGen,6.); // = 378cm eLoss_Common + 37 cm eLoss_steel
  Double_t eLossC = (eLoss02 - 305./37.*eLoss23) / (110 - 378.*305./37);
  Double_t eLossW = (eLoss02 - 110.*eLossC) / 305.;
  Double_t eLossS = (eLoss310 - 378.*eLossC) / 37.;
  Double_t fwhmELoss02 = TMath::Sqrt(FWHMELoss2(pGen, 1.5));
  Double_t fwhmELoss23 = TMath::Sqrt(FWHMELoss2(pGen, 2.5));
  Double_t fwhmELoss310 = TMath::Sqrt(FWHMELoss2(pGen, 6.));
  Double_t fwhmELossC = (fwhmELoss02 - 305./37.*fwhmELoss23) / (110 - 378.*305./37);
  Double_t fwhmELossW = (fwhmELoss02 - 110.*fwhmELossC) / 305.;
  Double_t fwhmELossS = (fwhmELoss310 - 378.*fwhmELossC) / 37.;

  // get mean energy loss and Branson plane for this eta value
  Double_t theta = 2.*TMath::ATan(TMath::Exp(etaGen))*TMath::RadToDeg();
  Double_t eLoss, fwhmELoss, zB;
  if (theta < 2.) {
    eLoss = eLoss02;
    fwhmELoss = fwhmELoss02;
    zB = fZB02;
  } else if (theta < 3.) {
    eLoss = eLoss23;
    fwhmELoss = fwhmELoss23;
    zB = fZB23;
  } else {
    eLoss = eLoss310;
    fwhmELoss = fwhmELoss310;
    zB = fZB310;
  }

  // MCS in the absorber
  Double_t pHalfLoss = pGen-0.5*TMath::Min(pGen-1.e-6,eLoss);
  Double_t slopeX = TMath::Cos(phiGen)/TMath::SinH(etaGen);
  Double_t slopeY = TMath::Sin(phiGen)/TMath::SinH(etaGen);
  Double_t sigmaMCSAbs = TMath::Sqrt(SigmaSlopeFromMCSInAbs2((theta < 2.) ? pGen : pHalfLoss, theta));
  Double_t slopeXMCS = GenRndGaus(slopeX,sigmaMCSAbs);
  Double_t slopeYMCS = GenRndGaus(slopeY,sigmaMCSAbs);

  // track position at the end of the absorber
  Double_t slopeXAbs, slopeYAbs;
  if (theta < 2) { // 305 cm of tungsten from 200 cm to 505 cm
    slopeXAbs = slopeXMCS + 2./TMath::Sqrt(3.) * (slopeXMCS - slopeX) * fZB23 / (fZB23 - 2.) * 3.05 / 5.05;
    slopeYAbs = slopeYMCS + 2./TMath::Sqrt(3.) * (slopeYMCS - slopeY) * fZB23 / (fZB23 - 2.) * 3.05 / 5.05;
  } else {
    slopeXAbs = slopeXMCS + 2./TMath::Sqrt(3.) * (slopeXMCS - slopeX) * zB / (zB - 0.9) * 4.15 / 5.05;
    slopeYAbs = slopeYMCS + 2./TMath::Sqrt(3.) * (slopeYMCS - slopeY) * zB / (zB - 0.9) * 4.15 / 5.05;
  }
  Double_t slopeRAbs = TMath::Sqrt(slopeXAbs*slopeXAbs + slopeYAbs*slopeYAbs);
  Double_t thetaAbs = TMath::ATan(slopeRAbs) * TMath::RadToDeg();
  fRAbs = slopeRAbs * 505.;

  // energy loss according to the "real" path in the absorber (should be done in 3D...)
  if ((theta > 2. && thetaAbs < 2.) || (theta > 3. && thetaAbs < 3.)) {

    Double_t lAt2out = (thetaAbs < 2.) ? 415.*TMath::Power((2.-theta)/(thetaAbs-theta),2./3.) : 415.;
    Double_t lAt2 = lAt2out*(505.-415.)/(505.-lAt2out);
    Double_t lAt3out = (theta > 3.) ? 415.*TMath::Power((3.-theta)/(thetaAbs-theta),2./3.) : 0.;
    Double_t lAt3 = lAt3out*(505.-415.)/(505.-lAt3out);

    Double_t lInC = TMath::Max(TMath::Min(lAt2,378.),110.);
    Double_t lInS = TMath::Max(lAt3-378.,0.);

    eLoss = (415.-lInC-lInS)*eLossW + lInS*eLossS + lInC*eLossC;
    fwhmELoss = (415.-lInC-lInS)*fwhmELossW + lInS*fwhmELossS + lInC*fwhmELossC;

  } else if ((theta < 2. && thetaAbs > 2.) || (theta < 3. && thetaAbs > 3.)) {

    Double_t lInMaterial = (theta < 2.) ? 305. : 415.;
    Double_t lAt2out = (theta < 2.) ? lInMaterial*TMath::Power((2.-theta)/(thetaAbs-theta),2./3.) : 0.;
    Double_t lAt2 = (415.-lInMaterial) + lAt2out*(505.-lInMaterial)/(505.-lAt2out);
    Double_t lAt3out = (thetaAbs > 3.) ? lInMaterial*TMath::Power((3.-theta)/(thetaAbs-theta),2./3.) : lInMaterial;
    Double_t lAt3 = (415.-lInMaterial) + lAt3out*(505.-lInMaterial)/(505.-lAt3out);

    Double_t lInC = 110. + TMath::Max(378.-TMath::Max(lAt2,110.),0.);
    Double_t lInS = 415.-TMath::Max(lAt3,378.);

    eLoss = (415.-lInC-lInS)*eLossW + lInS*eLossS + lInC*eLossC;
    fwhmELoss = (415.-lInC-lInS)*fwhmELossW + lInS*fwhmELossS + lInC*fwhmELossC;

  }

  // energy loss in the absorber
  Double_t dp = -1.;
  while (dp < 0.) dp = gRandom->Landau(eLoss+0.22278298*0.25*fwhmELoss,0.25*fwhmELoss);
  Double_t pAbsEnd = pGen - dp;
  if (pAbsEnd < 0.5) pAbsEnd = 0.5; // below ~0.5 GeV/c the track should not even exist

  // Branson plane according to thetaAbs
  Double_t zBRec;
  if (thetaAbs < 2.) zBRec = fZB02;
  else if (thetaAbs < 3.) zBRec = fZB23;
  else zBRec = fZB310;

  // compute reconstructed slopes
  Double_t sigmaMCSCh2 = SigmaSlopeFromMCSInCh2(pAbsEnd, kFALSE, zBRec);
  Double_t sigmaResX2 = SigmaSlopeFromRes2(kFALSE, kFALSE, zBRec);
  Double_t sigmaResY2 = SigmaSlopeFromRes2(kTRUE, kFALSE, zBRec);
  Double_t slopeXRec = 0., slopeYRec = 0.;
  if ( fChosenFunc == kGaus ) {
    slopeXRec = GenRndGaus(slopeXMCS,TMath::Sqrt(sigmaMCSCh2+sigmaResX2));
    slopeYRec = GenRndGaus(slopeYMCS,TMath::Sqrt(sigmaMCSCh2+sigmaResY2));
  } else {
    Double_t slopeXMCSCh = GenRndGaus(slopeXMCS,TMath::Sqrt(sigmaMCSCh2));
    Double_t slopeYMCSCh = GenRndGaus(slopeYMCS,TMath::Sqrt(sigmaMCSCh2));
    Double_t sigmaResX = TMath::Sqrt(sigmaResX2);
    Double_t sigmaResY = TMath::Sqrt(sigmaResY2);
    if ( fChosenFunc == kBreitWigner ) {
      slopeXRec = GenRndBreitWigner(slopeXMCSCh, sigmaResX, fSigmaTrkCut * sigmaResX / fSigmaxChSt1 * fSigmaTrk);
      slopeYRec = GenRndBreitWigner(slopeYMCSCh, sigmaResY, fSigmaTrkCut * sigmaResY / fSigmayChSt1 * fSigmaTrk);
    } else if ( fChosenFunc == kCrystalBall ) {
      slopeXRec = GenRndCrystalBall(slopeXMCSCh, sigmaResX, fCrystalBallTails[0], fCrystalBallTails[1], fSigmaTrkCut * sigmaResX / fSigmaxChSt1 * fSigmaTrk);
      slopeYRec = GenRndCrystalBall(slopeYMCSCh, sigmaResY, fCrystalBallTails[2], fCrystalBallTails[3], fSigmaTrkCut * sigmaResY / fSigmayChSt1 * fSigmaTrk);
    }
  }

  // compute reconstructed momentum at first cluster
  Double_t sigmaThetaRes2 = SigmaThetaDevFromRes2();
  Double_t sigmaThetaRes = TMath::Sqrt(sigmaThetaRes2);
  Double_t mean = PToThetaDev(pAbsEnd) + sign * fNSigmaShift * sigmaThetaRes;
  Double_t pAbsEndRec = 0.;
  if ( fChosenFunc == kGaus ) {
    Double_t sigmaThetaDev = TMath::Sqrt(SigmaThetaDevFromMCS2(pAbsEnd) + sigmaThetaRes2);
    pAbsEndRec = ThetaDevToP(GenRndGaus(mean,sigmaThetaDev));
  } else {
    Double_t thetaDevMCS = GenRndGaus(mean,TMath::Sqrt(SigmaThetaDevFromMCS2(pAbsEnd)));
    if ( fChosenFunc == kBreitWigner ) {
      pAbsEndRec = ThetaDevToP(GenRndBreitWigner(thetaDevMCS, sigmaThetaRes, fSigmaTrkCut * sigmaThetaRes / fSigmayCh * fSigmaTrk));
    } else if ( fChosenFunc == kCrystalBall ) {
      pAbsEndRec = ThetaDevToP(GenRndCrystalBall(thetaDevMCS, sigmaThetaRes, fCrystalBallTails[4], fCrystalBallTails[5], fSigmaTrkCut * sigmaThetaRes / fSigmayCh * fSigmaTrk));
    }
  }
  if (pAbsEndRec < 0.) {
    fRecoCharge *= -1.;
    pAbsEndRec = TMath::Abs(pAbsEndRec);
  }

  // compute reconstructed kinematics
  Double_t phiRec = TMath::ATan2(slopeYRec,slopeXRec) + TMath::Pi();
  Double_t etaRec = -TMath::ASinH(1./TMath::Sqrt(slopeXRec*slopeXRec + slopeYRec*slopeYRec));
  Double_t pRec = pAbsEndRec + ELoss(pAbsEndRec, thetaAbs);
  Double_t pTRec = pRec/TMath::CosH(etaRec);
//  Double_t pxRec = pTRec*TMath::Cos(phiRec);
//  Double_t pyRec = pTRec*TMath::Sin(phiRec);
//  Double_t pZRec = pTRec*TMath::SinH(etaRec);

  fRecoTrack.SetPtEtaPhiM(pTRec,etaRec,phiRec,fMuMass);

//  AliDebug(1,Form("Track smearing: p: %g -> %g  eta %g -> %g  phi %g -> %g  charge %g -> %g",pGen,pRec,etaGen,etaRec,phiGen,phiRec,chargeGen,fRecoCharge));
  AliDebug(1,Form("Track smearing: p: %g -> %g  eta %g -> %g  phi %g -> %g  charge %g -> %g",pGen,fRecoTrack.P(),etaGen,fRecoTrack.Eta(),phiGen,fRecoTrack.Phi(),chargeGen,fRecoCharge));

//
//  Double_t mTRec = TMath::Sqrt(fMuMass*fMuMass + pTRec*pTRec);
//  Double_t yRec = TMath::ASinH(pZRec/mTRec);
}


//________________________________________________________________________
Double_t AliMuonTrackSmearing::ELoss ( Double_t momentum, Double_t theta ) const
{
  // Returns the total momentum energy loss of muon in the absorber
  if (theta < 2.) return 9.3 + 0.05*momentum;
  else if (theta < 3.) return 2.9 + 0.0134*momentum;
  return 2.4 + 0.007*momentum;
}

//________________________________________________________________________
Double_t AliMuonTrackSmearing::ELossFluctuation2 ( Double_t momentum, Double_t rhoZoverA ) const
{
  /// Returns the total momentum energy loss fluctuation of muon in the absorber
  Double_t k = 0.307075e-3; // GeV.g^-1.cm^2
  Double_t pathLength = 415.;
  Double_t p2=momentum*momentum;
  Double_t beta2=p2/(p2 + fMuMass*fMuMass);
  Double_t fwhm = 2. * k * rhoZoverA * pathLength / beta2; // FWHM of the energy loss Landau distribution
  return fwhm*fwhm;
}

//________________________________________________________________________
Double_t AliMuonTrackSmearing::FWHMELoss2 ( Double_t momentum, Double_t theta ) const
{
  /// energy loss dispersion of muon in the absorber
  Double_t rhoZoverA; // effective rho * Z / A of the absorber (Carbone = 2.21*0.5; Be = 0.44*1.85)
  if (fTuneKalman) {
    Double_t corr = 2.; // FWHM are summed quadratically instead of linearly in AliMUONTrackExtrap --> bug to be corrected
    if (theta < 2.) rhoZoverA = 3.3*corr;
    else if (theta < 3.) rhoZoverA = 0.81*corr;
    else rhoZoverA = 0.62*corr;
  } else {
    if (theta < 2.) rhoZoverA = 11.;
    else if (theta < 3.) rhoZoverA = 2.7;
    else rhoZoverA = 1.7;
  }
  return ELossFluctuation2(momentum, rhoZoverA);
}


//________________________________________________________________________
Double_t AliMuonTrackSmearing::GenRndBreitWigner ( Double_t mean, Double_t sigma, Double_t max ) const
{
  /// generate a random number following a Breit-Wigner distribution in the range mean +- max
  Double_t val = 0.;
  do val = gRandom->BreitWigner(mean,sigma);
  while ( TMath::Abs(val-mean) > max );
  return val;
}


//________________________________________________________________________
Double_t AliMuonTrackSmearing::GenRndCrystalBall ( Double_t mean, Double_t sigma, Double_t tail1, Double_t tail2, Double_t max )
{
  /// generate a random number following a Crystal Ball distribution in the range mean +- max
  if ( ! fCrystalBall ) {
    fCrystalBall = new TF1("CrystalBall2", this, &AliMuonTrackSmearing::CrystalBallSymmetric, -2., 2., 5, "AliMuonTrackSmearing", "CrystalBallSymmetric");
    fCrystalBall->SetNpx(1000);
  }

  fCrystalBall->SetParameters(1.,mean,sigma,tail1,tail2);
  fCrystalBall->SetRange(mean-max,mean+max);
  return fCrystalBall->GetRandom();
}


//________________________________________________________________________
Double_t AliMuonTrackSmearing::GenRndGaus ( Double_t mean, Double_t sigma ) const
{
  /// generate a random number following a gaussian distribution
  return gRandom->Gaus(mean,sigma);
}

//________________________________________________________________________
//AliVParticle* AliMuonTrackSmearing::GetRecoTrack ( const AliVParticle* recoTrack, const AliMCEvent* mcEvent )
//{
//  /// Given the reconstructed MC track, recovers the generated particle
//  /// and returns the reconstructed track according to resolution
//
//  if ( recoTrack->GetLabel() < 0 ) return 0x0;
//  AliVParticle* genParticle = mcEvent->GetTrack(recoTrack->GetLabel());
//  return GetRecoTrack(genParticle,recoTrack);
//}

//________________________________________________________________________
TLorentzVector AliMuonTrackSmearing::GetRecoTrack ( Double_t pGen, Double_t etaGen, Double_t phiGen, Double_t chargeGen, Double_t &recoCharge, Double_t &rAbs )
{
  /// Given the generated MC track parameters,
  /// returns the track with smeared parameters according to resolution

  ComputeRecoTrack(pGen,etaGen,phiGen,chargeGen);
  recoCharge = fRecoCharge;
  rAbs = fRAbs;
  return fRecoTrack;
}

////________________________________________________________________________
//AliVParticle* AliMuonTrackSmearing::GetRecoTrack ( const AliVParticle* genParticle )
//{
//  /// Given the generated MC Track,
//  /// returns the track with smeared parameters according to resolution
//  return GetRecoTrack(genParticle,(AliVParticle*)0x0);
//}
//
////________________________________________________________________________
//AliVParticle* AliMuonTrackSmearing::GetRecoTrack ( const AliVParticle* genParticle, const AliVParticle* recoTrack )
//{
//  /// Given the generated MC track,
//  /// returns the track with smeared parameters according to resolution
//  /// If the reconstructed track from the standard MC is provided as well
//  /// the returned track is a clone of it with modified momentum parameters
//  /// otherwise create a new AliAODTrack
//
//  ComputeRecoTrack(genParticle->P(),genParticle->Eta(),genParticle->Phi(),genParticle->Charge());
//  AliAODTrack* aodTrack = 0x0;
//  if ( recoTrack && ( recoTrack->IsA() == AliAODTrack::Class() ) ) {
//    aodTrack = static_cast<AliAODTrack*>(recoTrack->Clone());
//    aodTrack->SetPt(fRecoTrack.Pt());
//    aodTrack->SetPhi(fRecoTrack.Phi());
//    aodTrack->SetTheta(fRecoTrack.Theta());
//    aodTrack->SetCharge(fRecoCharge);
//  }
//  else {
//    Double_t pVec[3] = {fRecoTrack.Px(),fRecoTrack.Py(),fRecoTrack.Pz()};
//    Double_t pos[3] = {genParticle->Xv(),genParticle->Yv(),genParticle->Zv()};
//    UChar_t itsClusMap(0);
//    UInt_t selectInfo(0);
//    aodTrack = new AliAODTrack (genParticle->GetUniqueID(), // ID
//                                0,                          // label
//                                pVec,                       // momentum
//                                kTRUE,                      // cartesian coordinate system
//                                pos,                        // position
//                                kFALSE,                     // isDCA
//                                0x0,                        // covariance matrix
//                                fRecoCharge,                // charge
//                                itsClusMap,                 // ITSClusterMap
//                                0x0,                        // origin vertex
//                                kFALSE,                     // used for vertex fit?
//                                kFALSE,                     // used for primary vertex fit?
//                                AliAODTrack::kPrimary,      // track type
//                                selectInfo);
//    aodTrack->SetPIDForTracking(AliPID::kMuon);
//  }
//  AliDebug(1,Form("px %g -> %g  py %g -> %g  pz %g -> %g",genParticle->Px(),aodTrack->Px(),genParticle->Py(),aodTrack->Py(),genParticle->Pz(),aodTrack->Pz()));
//  fRecoTrackList.push_back(aodTrack);
//  return aodTrack;
//}

//________________________________________________________________________
Double_t AliMuonTrackSmearing::MCS2 ( Double_t momentum, Double_t dZ, Double_t x0) const
{
  /// Return the angular dispersion square due to multiple Coulomb scattering
  /// through a material of thickness "dZ" and of radiation length "x0"
  Double_t p2=momentum*momentum;
  Double_t betap=p2/TMath::Sqrt(p2 + fMuMass*fMuMass);
  Double_t theta02 = 0.0136 / betap * (1 + 0.038 * TMath::Log(dZ/x0));
  return theta02 * theta02 * dZ / x0;
}


//________________________________________________________________________
Double_t AliMuonTrackSmearing::PToThetaDev ( Double_t momentum ) const
{
  /// deviation angle for a track of a given momentum in the spectrometer
  return 3./momentum*0.3; // BL = 3Tm; GeV/c = 10^9/3.10^8
}


//________________________________________________________________________
Double_t AliMuonTrackSmearing::SigmaSlopeFromMCSInAbs2 ( Double_t momentum, Double_t theta ) const
{
  //// slope dispersion due to the MCS in the front absorber
  Double_t x0;
  //  if (theta < 2.) x0 = (tuneKalman) ? 1. : 0.5;
  //  else x0 = (tuneKalman) ? 35.4 : 19.3; // X0 Carbone = 42.7/2.21; Be = 65.2/1.85
  if (theta < 2.) x0 = 1.;
  else x0 = 35.4; // X0 Carbone = 42.7/2.21; Be = 65.2/1.85
  return MCS2(momentum, 415., x0) / 4.; // sigma reduced by ~2 by the Branson correction
}


//________________________________________________________________________
Double_t AliMuonTrackSmearing::SigmaSlopeFromMCSInCh2 ( Double_t momentum, Bool_t at1stCl, Double_t zB ) const
{
  /// slope dispersion due to the MCS in chambers
  Double_t sigmaMCS2 = (fTuneKalman) ? 2.*MCS2(momentum, 0.065, 1.) : 2.*MCS2(momentum, 0.065, 1.) + 2.*MCS2(momentum, 0.075, 1.);
  if (at1stCl) return sigmaMCS2;
  return sigmaMCS2*(5.36-zB)*(5.36-zB)/zB/zB; // propagate to the Branson plane then to the vertex
}


//________________________________________________________________________
Double_t AliMuonTrackSmearing::SigmaSlopeFromRes2 ( Bool_t bendingDir, Bool_t at1stCl, Double_t zB ) const
{
  /// slope dispersion due to chamber resolution
  static Double_t magicFactorX = 1./23.; // to get angular resolution^2 from station 1 x-resolution^2
  static Double_t magicFactorY = 1./3.; // to get angular resolution^2 from station 1 y-resolution^2
  Double_t sigmaSt12, sigmaSlope2;
  if (bendingDir) {
    sigmaSt12 = fSigmayChSt1*fSigmayChSt1/2.;
    sigmaSlope2 = magicFactorY*sigmaSt12;
  } else {
    sigmaSt12 = fSigmaxChSt1*fSigmaxChSt1/2.;
    sigmaSlope2 = magicFactorX*sigmaSt12;
  }
  if (at1stCl) return sigmaSlope2;
  else return (sigmaSlope2*(5.36-zB)*(5.36-zB) + sigmaSt12)/zB/zB; // propagate to the Branson plane then to the vertex
}


//________________________________________________________________________
Double_t AliMuonTrackSmearing::SigmaThetaDevFromMCS2 ( Double_t momentum ) const
{
  /// resolution of the deviation angle for a track of a given momentum in the spectrometer
  return fTuneKalman ? 4.*MCS2(momentum, 0.035, 1.) : 2*MCS2(momentum, 0.075, 1.) + 4.*MCS2(momentum, 0.035, 1.);
}

//________________________________________________________________________
Double_t AliMuonTrackSmearing::SigmaThetaDevFromRes2() const
{
  /// resolution of the deviation angle for a track of a given momentum in the spectrometer
  static Double_t magicFactor = 1./2.2; // to get thetaDev resolution^2 from averaged chamber resolution^2
  return magicFactor*fSigmayCh*fSigmayCh;
}

//________________________________________________________________________
Double_t AliMuonTrackSmearing::ThetaDevToP ( Double_t thetaDev ) const
{
  /// deviation angle for a track of a given momentum in the spectrometer
  return 3./thetaDev*0.3; // BL = 3Tm; GeV/c = 10^9/3.10^8
}

//________________________________________________________________________
Double_t AliMuonTrackSmearing::PResVsP( const Double_t *x, const Double_t *par )
{
  /// expected momentum resolution versus p
  // par[0] = angular position in the absorber
  // par[1] < 0. = resolution at vertex; par[1] > 0. = resolution at first cluster
  Double_t p = *x;
  Double_t theta = par[0];
  Bool_t atFirstCluster = (par[1] > 0.);
  Double_t dp = ELoss(p, theta);
  Double_t pCorr = TMath::Max(atFirstCluster ? p : p - dp, 0.5);
  Double_t thetaDev = PToThetaDev(pCorr);
  Double_t sigmaThetaDev2 = SigmaThetaDevFromRes2();
  if (fChosenFunc == kBreitWigner) sigmaThetaDev2 /= 2.*log(2.); // FWHM = 2 * gamma = 2 * srqt(2*ln(2)) * sigma
  sigmaThetaDev2 += SigmaThetaDevFromMCS2(pCorr);
  if (atFirstCluster) return 100.*TMath::Sqrt(sigmaThetaDev2)/thetaDev;
  else {
    Double_t sigmaELoss2 = FWHMELoss2(p, theta) / (8.*log(2.)); // gaussian: fwmh = 2 * srqt(2*ln(2)) * sigma
    return 100.*TMath::Sqrt(sigmaELoss2 + sigmaThetaDev2/thetaDev/thetaDev*pCorr*pCorr)/p;
  }
}

//________________________________________________________________________
Double_t AliMuonTrackSmearing::SlopeResVsP( const Double_t *x, const Double_t *par )
{
  /// expected slope resolution versus p
  // par[0] = angular position in the absorber
  // par[1] < 0. = resolution along X; par[1] > 0. = resolution along Y
  // par[2] < 0. = resolution at vertex; par[1] > 0. = resolution at first cluster
  Double_t p = *x;
  Double_t theta = par[0];
  Bool_t atFirstCluster = (par[2] > 0.);
  Double_t dp = ELoss(p, theta);
  Double_t zB;
  if (theta < 2.) zB = fZB02;
  else if (theta < 3.) zB = fZB23;
  else zB = fZB310;
  Double_t pCorr = TMath::Max(atFirstCluster ? p : p - dp, 0.5);
  Double_t sigmaMCSCh2 = SigmaSlopeFromMCSInCh2(pCorr, atFirstCluster, zB);
  Double_t sigmaRes2 = SigmaSlopeFromRes2((par[1] > 0.), atFirstCluster, zB);
  if (fChosenFunc == kBreitWigner) sigmaRes2 /= 2.*log(2.); // FWHM = 2 * gamma = 2 * srqt(2*ln(2)) * sigma
  if (atFirstCluster) return TMath::Sqrt(sigmaMCSCh2+sigmaRes2);
  else {
    Double_t pHalfCorr = p-0.5*TMath::Min(p-1.e-6,dp);
    Double_t sigmaMCSAbs2 = SigmaSlopeFromMCSInAbs2(pHalfCorr,theta);
    return TMath::Sqrt(sigmaMCSAbs2+sigmaMCSCh2+sigmaRes2);
  }
}

