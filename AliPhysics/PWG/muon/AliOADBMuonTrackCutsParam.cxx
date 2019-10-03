/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliOADBMuonTrackCutsParam.h"

#include "TVector3.h"

#include "AliLog.h"

using namespace std;

ClassImp(AliOADBMuonTrackCutsParam)


//________________________________________________________________________
AliOADBMuonTrackCutsParam::AliOADBMuonTrackCutsParam () :
TNamed("AliOADBMuonTrackCutsParam", "OADB object for Muon track cuts"),
fMeanDcaX(0.),
fMeanDcaY(0.),
fMeanDcaZ(0.),
fMeanPCorr23(0.),
fMeanPCorr310(0.),
fSigmaPdca23(0.),
fSigmaPdca310(0.),
fNSigmaPdcaCut(0.),
fChi2NormCut(0.),
fRelPResolution(0.),
fSlopeResolution(0.),
fSharpPtApt(0.),
fSharpPtLpt(0.),
fSharpPtHpt(0.)
{
  // default ctor
}


//________________________________________________________________________
AliOADBMuonTrackCutsParam::AliOADBMuonTrackCutsParam ( const char* name ) :
TNamed(name, "OADB object for Muon track cuts"),
fMeanDcaX(0.),
fMeanDcaY(0.),
fMeanDcaZ(0.),
fMeanPCorr23(0.),
fMeanPCorr310(0.),
fSigmaPdca23(0.),
fSigmaPdca310(0.),
fNSigmaPdcaCut(0.),
fChi2NormCut(0.),
fRelPResolution(0.),
fSlopeResolution(0.),
fSharpPtApt(0.),
fSharpPtLpt(0.),
fSharpPtHpt(0.)
{
  // ctor, better use this one
}


//________________________________________________________________________
AliOADBMuonTrackCutsParam::~AliOADBMuonTrackCutsParam()
{
  // dtor
}

//________________________________________________________________________
AliOADBMuonTrackCutsParam::AliOADBMuonTrackCutsParam ( const AliOADBMuonTrackCutsParam& other ) :
TNamed ( other ),
fMeanDcaX ( other.fMeanDcaX ),
fMeanDcaY ( other.fMeanDcaY ),
fMeanDcaZ ( other.fMeanDcaZ ),
fMeanPCorr23 ( other.fMeanPCorr23 ),
fMeanPCorr310 ( other.fMeanPCorr310 ),
fSigmaPdca23 ( other.fSigmaPdca23 ),
fSigmaPdca310 ( other.fSigmaPdca310 ),
fNSigmaPdcaCut ( other.fNSigmaPdcaCut ),
fChi2NormCut ( other.fChi2NormCut ),
fRelPResolution ( other.fRelPResolution ),
fSlopeResolution ( other.fSlopeResolution ),
fSharpPtApt ( other.fSharpPtApt ),
fSharpPtLpt ( other.fSharpPtLpt ),
fSharpPtHpt ( other.fSharpPtHpt )
{
// Copy ctor  
}


//________________________________________________________________________
AliOADBMuonTrackCutsParam& AliOADBMuonTrackCutsParam::operator=(const AliOADBMuonTrackCutsParam& other)
{
  //Assignment operator
  if ( &other == this ) return *this;
  TNamed::operator=(other);

  fMeanDcaX = other.fMeanDcaX;
  fMeanDcaY = other.fMeanDcaY;
  fMeanDcaZ = other.fMeanDcaZ;
  fMeanPCorr23 = other.fMeanPCorr23;
  fMeanPCorr310 = other.fMeanPCorr310;
  fSigmaPdca23 = other.fSigmaPdca23;
  fSigmaPdca310 = other.fSigmaPdca310;
  fNSigmaPdcaCut = other.fNSigmaPdcaCut;
  fChi2NormCut = other.fChi2NormCut;
  fRelPResolution = other.fRelPResolution;
  fSlopeResolution = other.fSlopeResolution;
  fSharpPtApt = other.fSharpPtApt;
  fSharpPtLpt = other.fSharpPtLpt;
  fSharpPtHpt = other.fSharpPtHpt;

  return *this;
}

//________________________________________________________________________
void AliOADBMuonTrackCutsParam::SetMeanDCA ( Double_t xAtDca, Double_t yAtDca, Double_t zAtDca )
{
  /// Set mean DCA from track
  fMeanDcaX = xAtDca;
  fMeanDcaY = yAtDca;
  fMeanDcaZ = zAtDca;
}

//________________________________________________________________________
TVector3 AliOADBMuonTrackCutsParam::GetMeanDCA () const
{ 
  /// Get mean DCA from track
  return TVector3 ( fMeanDcaX, fMeanDcaY, fMeanDcaZ );
}


//________________________________________________________________________
void AliOADBMuonTrackCutsParam::SetMeanPCorr ( Double_t pCorrThetaAbs23, Double_t pCorrThetaAbs310 )
{
  /// Set mean p correction
  fMeanPCorr23 = pCorrThetaAbs23;
  fMeanPCorr310 = pCorrThetaAbs310;
}

//________________________________________________________________________
Double_t AliOADBMuonTrackCutsParam::GetMeanPCorr23 ( ) const
{
  /// Get mean p correction in 2<theta_abs<3 deg
  return fMeanPCorr23;
}

//________________________________________________________________________
Double_t AliOADBMuonTrackCutsParam::GetMeanPCorr310 ( ) const
{
  /// Get mean p correction in 3<theta_abs<10 deg
  return fMeanPCorr310;
}

//________________________________________________________________________
void AliOADBMuonTrackCutsParam::SetSigmaPdca ( Double_t sigmaThetaAbs23, Double_t sigmaThetaAbs310 )
{ 
  /// Set sigma pdca
  fSigmaPdca23 = sigmaThetaAbs23;
  fSigmaPdca310 = sigmaThetaAbs310;
}

//________________________________________________________________________
Double_t AliOADBMuonTrackCutsParam::GetSigmaPdca23 ( ) const
{ 
  /// Get mean pdca in 2<theta_abs<3 deg
  return fSigmaPdca23;
}

//________________________________________________________________________
Double_t AliOADBMuonTrackCutsParam::GetSigmaPdca310 ( ) const
{ 
  /// Get mean pdca in 3<theta_abs<10 deg
  return fSigmaPdca310;
}

//________________________________________________________________________
void AliOADBMuonTrackCutsParam::SetNSigmaPdca ( Double_t nSigmas )
{ 
  /// Set N sigma pdca cut
  fNSigmaPdcaCut = nSigmas;
}

//________________________________________________________________________
Double_t AliOADBMuonTrackCutsParam::GetNSigmaPdca () const
{
  /// Get N sigma pdca cut
  return fNSigmaPdcaCut;
}

//________________________________________________________________________
void AliOADBMuonTrackCutsParam::SetChi2NormCut ( Double_t chi2normCut )
{
  /// Set cut on normalized chi2 of tracks
  fChi2NormCut = chi2normCut;
}

//________________________________________________________________________
Double_t AliOADBMuonTrackCutsParam::GetChi2NormCut () const
{
  /// Get cut on normalized chi2 of tracks
  return fChi2NormCut;
}

//________________________________________________________________________
void AliOADBMuonTrackCutsParam::SetRelPResolution ( Double_t relPResolution )
{
  /// Set relative momentum resolution
  fRelPResolution = relPResolution;
}

//________________________________________________________________________
Double_t AliOADBMuonTrackCutsParam::GetRelPResolution () const
{
  /// Get relative momentum resolution
  return fRelPResolution;
}


//________________________________________________________________________
void AliOADBMuonTrackCutsParam::SetSlopeResolution ( Double_t slopeResolution )
{
  /// Set slope resolution
  fSlopeResolution = slopeResolution;
}

//________________________________________________________________________
Double_t AliOADBMuonTrackCutsParam::GetSlopeResolution () const
{
  /// Get slope resolution
  return fSlopeResolution;
}

//________________________________________________________________________
void AliOADBMuonTrackCutsParam::SetSharpPtCut ( Double_t valueApt, Double_t valueLpt, Double_t valueHpt  )
{
  /// Set sharp tracker cut matching the trigger level
  
  fSharpPtApt = valueApt;
  fSharpPtLpt = valueLpt;
  fSharpPtHpt = valueHpt;
}

//________________________________________________________________________
Double_t AliOADBMuonTrackCutsParam::GetSharpPtCut ( Int_t trigPtCut, Bool_t warn ) const
{
  /// Get sharp tracker cut matching the trigger level
  /// trigPtCut can be 0 (Apt), 1 (Lpt) or 2 (Hpt)
  switch ( trigPtCut ) {
  case 0:
    return fSharpPtApt;
  case 1:
    return fSharpPtLpt;
  case 2:
    return fSharpPtHpt;
  }
  
  if ( warn ) AliError("Allowed values for trigPtCut are 0 (Apt), 1 (Lpt), 2 (Hpt)");
  return 0.;
}

//________________________________________________________________________
void AliOADBMuonTrackCutsParam::Print ( Option_t* /*option*/ ) const
{
  /// Print info
  printf(" *** Muon track parameter summary: ***\n");
  printf("  Mean vertex DCA: (%g, %g, %g)\n", fMeanDcaX, fMeanDcaY, fMeanDcaZ);
  printf("  Mean p correction (GeV/c): theta2-3 = %g  theta3-10 = %g\n", fMeanPCorr23, fMeanPCorr310);
  printf("  Sigma p x DCA (cm x GeV/c): theta2-3 = %g  theta3-10 = %g\n", fSigmaPdca23, fSigmaPdca310);
  printf("  Cut p x DCA in units of sigma: %g\n", fNSigmaPdcaCut);
  printf("  Cut on track chi2/NDF: %g\n", fChi2NormCut);
  printf("  Momentum resolution: %g\n", fRelPResolution);
  printf("  Slope resolution: %g\n", fSlopeResolution);
  printf("  Sharp pt cut: %g (Apt)  %g (Lpt)  %g (Hpt)\n", fSharpPtApt, fSharpPtLpt, fSharpPtHpt);
  printf(" ********************************\n");
}
