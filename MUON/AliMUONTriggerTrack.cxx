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

// $Id$

//-----------------------------------------------------------------------------
// Class AliMUONTriggerTrack
//---------------------------
// Reconstructed Trigger track in ALICE dimuon spectrometer
// Note: equivalent to AliMUONTriggerTrack for tracking,
// No need for a AliMUONTriggerTrackParam
// Author: Philippe Crochet
//-----------------------------------------------------------------------------

#include "AliMUONTriggerTrack.h"
#include "AliMUONTrackReconstructor.h" 
#include "TString.h"
#include <Riostream.h>
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerTrack)
/// \endcond

//__________________________________________________________________________
AliMUONTriggerTrack::AliMUONTriggerTrack()
  : TObject(),
    fx11(0),
    fy11(0),
    fz11(0.),
    fz21(0.),
    fSlopeX(0),
    fSlopeY(0),
    floTrgNum(0),
    fGTPattern(0),
    fHitsPatternInTrigCh(0),
    fCovariances(0x0)
{
  /// default ctr
      AliDebug(5,Form("this=%p",this));
}
//__________________________________________________________________________
AliMUONTriggerTrack::AliMUONTriggerTrack(Float_t x11, Float_t y11, Float_t z11, Float_t z21, Float_t slopeX, Float_t slopeY, Int_t loTrgNum, Long_t theGTPattern, UShort_t hitsPatternInTrigCh)
    : TObject(),
      fx11(x11),
      fy11(y11),
      fz11(z11),
      fz21(z21),
      fSlopeX(slopeX),
      fSlopeY(slopeY),
      floTrgNum(loTrgNum),
      fGTPattern(theGTPattern),
      fHitsPatternInTrigCh(hitsPatternInTrigCh),
      fCovariances(0x0)
{
/// ctor from local trigger output
        AliDebug(5,Form("this=%p x11=%f y11=%f z11=%f z21=%f slopeX=%f slopeY=%f loTrgNum=%d GTPattern=%ld HitsPatternInTrigCh %i",
                        this,x11,y11,z11,z21,slopeX,slopeY,loTrgNum,theGTPattern,fHitsPatternInTrigCh));

}

//__________________________________________________________________________
AliMUONTriggerTrack::~AliMUONTriggerTrack()
{
  /// Destructor
  AliDebug(5,Form("this=%p",this));
  if (fCovariances) {
    delete fCovariances;
    fCovariances = 0x0;
  }
}

//__________________________________________________________________________
AliMUONTriggerTrack::AliMUONTriggerTrack (const AliMUONTriggerTrack& theMUONTriggerTrack)
    : TObject(theMUONTriggerTrack),
      fx11(theMUONTriggerTrack.fx11),
      fy11(theMUONTriggerTrack.fy11),
      fz11(theMUONTriggerTrack.fz11),
      fz21(theMUONTriggerTrack.fz21),
      fSlopeX(theMUONTriggerTrack.fSlopeX),
      fSlopeY(theMUONTriggerTrack.fSlopeY),
      floTrgNum(theMUONTriggerTrack.floTrgNum),
      fGTPattern(theMUONTriggerTrack.fGTPattern),
      fHitsPatternInTrigCh(theMUONTriggerTrack.fHitsPatternInTrigCh),
      fCovariances(0x0)
{
///
/// copy ctor
///
  if (theMUONTriggerTrack.fCovariances) fCovariances = new TMatrixD(*(theMUONTriggerTrack.fCovariances));
  AliDebug(5,Form("this=%p copy ctor",this));

}
      
//__________________________________________________________________________
AliMUONTriggerTrack & AliMUONTriggerTrack::operator=(const AliMUONTriggerTrack&
theMUONTriggerTrack)
{
/// Assignment operator

    // check assignement to self
    if (this == &theMUONTriggerTrack)
	return *this;
    
    /// base class assignement
    TObject::operator=(theMUONTriggerTrack);

    fx11 = theMUONTriggerTrack.fx11;
    fy11 = theMUONTriggerTrack.fy11;
    fz11 = theMUONTriggerTrack.fz11;
    fz21 = theMUONTriggerTrack.fz21;
    fSlopeX = theMUONTriggerTrack.fSlopeX;
    fSlopeY = theMUONTriggerTrack.fSlopeY;
    floTrgNum = theMUONTriggerTrack.floTrgNum;
    fGTPattern = theMUONTriggerTrack.fGTPattern;
    fHitsPatternInTrigCh = theMUONTriggerTrack.fHitsPatternInTrigCh;

    if (theMUONTriggerTrack.fCovariances) {
      if (fCovariances) *fCovariances = *(theMUONTriggerTrack.fCovariances);
      else fCovariances = new TMatrixD(*(theMUONTriggerTrack.fCovariances));
    } else {
      delete fCovariances;
      fCovariances = 0x0;
    }

    return *this;
}

//__________________________________________________________________________
void
AliMUONTriggerTrack::Print(Option_t* opt) const
{
/// Printing
  TString optString(opt);
  optString.ToUpper();
  if ( optString.Contains("FULL") ) optString = "PARAM COV";

  if ( optString.Contains("PARAM"))
    cout << Form("(X,Y,Z)11=(%7.2f,%7.2f,%7.2f) Z21=%7.2f Slope(X,Y)=(%7.2f,%7.2f) LocalBoard #%3d GlobalTriggerPattern %x HitsPatternInTrigCh %x",
		 fx11,fy11,fz11,fz21,fSlopeX,fSlopeY,floTrgNum,fGTPattern,fHitsPatternInTrigCh) << endl;

  if ( optString.Contains("COV") ){
    if ( ! fCovariances ) cout << "Covariances not initialized " << endl;
    else fCovariances->Print();
  }
}

//__________________________________________________________________________
void AliMUONTriggerTrack::SetCovariances(const TMatrixD& covariances)
{
  /// Set the covariance matrix
  if (fCovariances) *fCovariances = covariances;
  else fCovariances = new TMatrixD(covariances);
}

//__________________________________________________________________________
void AliMUONTriggerTrack::SetCovariances(const Double_t matrix[3][3])
{
  /// Set the covariance matrix
  if (fCovariances) fCovariances->SetMatrixArray(&(matrix[0][0]));
  else fCovariances = new TMatrixD(3,3,&(matrix[0][0]));
}

//__________________________________________________________________________
const TMatrixD& AliMUONTriggerTrack::GetCovariances() const
{
  /// Return the covariance matrix (create it before if needed)
  if (!fCovariances) {
    fCovariances = new TMatrixD(3,3);
    fCovariances->Zero();
  }
  return *fCovariances;
}

//__________________________________________________________________________
Bool_t AliMUONTriggerTrack::Match(AliMUONTriggerTrack &track,
				  Double_t sigmaCut) const
{
  /// Try to match this track with the given track. Matching conditions:
  /// - x, y position and y slope within sigmaCut
  
  // Find the track with the covariances correctly set
  // Extrapolate to the z of the other track
  Bool_t hasCov1 = ( GetCovariances().NonZeros() != 0 );
  Bool_t hasCov2 = ( track.GetCovariances().NonZeros() != 0 );

  const AliMUONTriggerTrack* trackToExtrap = ( hasCov2 ) ? &track : this;
  const AliMUONTriggerTrack* fixedTrack = ( hasCov2 ) ? this : &track;

  TMatrixD paramDiff(3,1);
  Double_t deltaZ = fixedTrack->GetZ11() - trackToExtrap->GetZ11();
  paramDiff(0,0) = fixedTrack->GetX11() - trackToExtrap->GetX11();
  paramDiff(1,0) = fixedTrack->GetY11() - ( trackToExtrap->GetY11() + trackToExtrap->GetSlopeY() * deltaZ );
  paramDiff(2,0) = fixedTrack->GetSlopeY() - trackToExtrap->GetSlopeY();
  Double_t chi2 = 0.;

  TMatrixD cov1(fixedTrack->GetCovariances());
  TMatrixD cov2(trackToExtrap->GetCovariances());

  // Extrapolate covariances to z
  if ( deltaZ != 0 ) {
    if ( hasCov1 || hasCov2 ){
      TMatrixD jacob(3,3);
      jacob.UnitMatrix();
      jacob(1,2) = deltaZ;
      TMatrixD tmp(trackToExtrap->GetCovariances(),TMatrixD::kMultTranspose,jacob);
      TMatrixD tmp2(jacob,TMatrixD::kMult,tmp);
      cov2 = tmp2;
    }
  }

  AliDebug(3, Form("track1 Y11 %f  track2 Y11: %f (Z11 %f)  -> %f (Z11 %f)", fixedTrack->GetY11(), trackToExtrap->GetY11(), trackToExtrap->GetZ11(), trackToExtrap->GetY11() + trackToExtrap->GetSlopeY() * deltaZ, fixedTrack->GetZ11()));

  TMatrixD sumCov(cov1,TMatrixD::kPlus,cov2);
  if (sumCov.Determinant() != 0) {
    sumCov.Invert();      
    TMatrixD tmp(sumCov,TMatrixD::kMult,paramDiff);
    TMatrixD chi2M(paramDiff,TMatrixD::kTransposeMult,tmp);
    chi2 = chi2M(0,0);
  } else {
    AliWarning(" Determinant = 0");
    Double_t sigma2 = 0.;
    for (Int_t iVar = 0; iVar < 3; iVar++) {
      sigma2 = cov1(iVar,iVar) + cov2(iVar,iVar);
      chi2 += paramDiff(iVar,0) * paramDiff(iVar,0) / sigma2;
    }
  }

  if ( chi2/3 > sigmaCut * sigmaCut )
    return kFALSE;
  
  return kTRUE;
}
