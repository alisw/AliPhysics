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

/* $Id$ */

 
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD cluster                                                              //
//                                                                           //
/////////////////////////////////////////////////////////////////////////////// 

#include "TMath.h"

#include "AliLog.h"
#include "AliTRDcluster.h"
#include "AliTRDgeometry.h"
#include "AliTRDCommonParam.h"

ClassImp(AliTRDcluster)

//___________________________________________________________________________
AliTRDcluster::AliTRDcluster() 
  :AliCluster() 
  ,fPadCol(0)
  ,fPadRow(0)
  ,fPadTime(0)
  ,fLocalTimeBin(0)
  ,fNPads(0)
  ,fClusterMasking(0)
  ,fDetector(0)
  ,fQ(0)
  ,fCenter(0)
{ 
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 7; i++) {
    fSignals[i] = 0;
  }

}

//___________________________________________________________________________
AliTRDcluster::AliTRDcluster(Int_t det, Float_t q
                           , Float_t *pos, Float_t *sig
                           , Int_t *tracks, Char_t npads, Short_t *signals
                           , UChar_t col, UChar_t row, UChar_t time
                           , Char_t timebin, Float_t center, UShort_t volid)
  :AliCluster(volid,pos[0],pos[1],pos[2],sig[0],sig[1],0.0,0x0) 
  ,fPadCol(col)
  ,fPadRow(row)
  ,fPadTime(time)
  ,fLocalTimeBin(timebin)
  ,fNPads(npads)
  ,fClusterMasking(0)
  ,fDetector(det)
  ,fQ(q)
  ,fCenter(center)
{ 
  //
  // Constructor
  //

  for (Int_t i = 0; i < 7; i++) {
    fSignals[i] = signals[i];
  }

  if (tracks) {
    AddTrackIndex(tracks);
  }

}

//_____________________________________________________________________________
AliTRDcluster::AliTRDcluster(const AliTRDcluster &c)
  :AliCluster(c)
  ,fPadCol(c.fPadCol)
  ,fPadRow(c.fPadRow)
  ,fPadTime(c.fPadTime)
  ,fLocalTimeBin(c.fLocalTimeBin)
  ,fNPads(c.fNPads)
  ,fClusterMasking(c.fClusterMasking)
  ,fDetector(c.fDetector)
  ,fQ(c.fQ)
  ,fCenter(c.fCenter)
{
  //
  // Copy constructor 
  //

  SetBit(kInChamber, c.IsInChamber());
  SetLabel(c.GetLabel(0),0);
  SetLabel(c.GetLabel(1),1);
  SetLabel(c.GetLabel(2),2);

  SetY(c.GetY());
  SetZ(c.GetZ());
  SetSigmaY2(c.GetSigmaY2());
  SetSigmaZ2(c.GetSigmaZ2());  

  for (Int_t i = 0; i < 7; i++) {
    fSignals[i] = c.fSignals[i];
  }

}

//_____________________________________________________________________________
void AliTRDcluster::AddTrackIndex(Int_t *track)
{
  //
  // Adds track index. Currently assumed that track is an array of
  // size 9, and up to 3 track indexes are stored in fTracks[3].
  // Indexes are sorted according to:
  //  1) index of max number of appearances is stored first
  //  2) if two or more indexes appear equal number of times, the lowest
  //     ones are stored first;
  //

  const Int_t kSize = 9;
  Int_t  entries[kSize][2];

  Int_t  i = 0;
  Int_t  j = 0;
  Int_t  k = 0;
  Int_t  index;
  Bool_t indexAdded;

  for (i = 0; i < kSize; i++) {
    entries[i][0] = -1;
    entries[i][1] =  0;
  }                                 

  for (k = 0; k < kSize; k++) {

    index      = track[k];
    indexAdded = kFALSE; 

    j = 0;
    if (index >= 0) {
      while ((!indexAdded) && (j < kSize)) {
        if ((entries[j][0] == index) || 
            (entries[j][1] ==     0)) {
          entries[j][0] = index;
          entries[j][1] = entries[j][1] + 1;
          indexAdded    = kTRUE;
        }
        j++;
      }
    }

  }

  // Sort by number of appearances and index value
  Int_t swap = 1;
  Int_t tmp0;
  Int_t tmp1;
  while (swap > 0) {
    swap = 0;
    for (i = 0; i < (kSize - 1); i++) {
      if ((entries[i][0]   >= 0) && 
          (entries[i+1][0] >= 0)) {
        if ((entries[i][1] < entries[i+1][1]) ||
            ((entries[i][1] == entries[i+1][1]) &&
             (entries[i][0] >  entries[i+1][0]))) {
          tmp0            = entries[i][0];
          tmp1            = entries[i][1];
          entries[i][0]   = entries[i+1][0];
          entries[i][1]   = entries[i+1][1];
          entries[i+1][0] = tmp0;
          entries[i+1][1] = tmp1;
          swap++;
        }
      }
    }
  }               

  // Set track indexes
  for (i = 0; i < 3; i++) {
    SetLabel(entries[i][0],i);
  }

  return;

}          

//_____________________________________________________________________________
void AliTRDcluster::Clear(Option_t *)
{
  //
  // Reset all member to the default value
  //
  fPadCol=0;
  fPadRow=0;
  fPadTime=0;
  fLocalTimeBin=0;
  fNPads=0;
  fClusterMasking=0;
  fDetector=0;
  for (Int_t i=0; i < 7; i++) fSignals[i]=0;
  fQ = 0;
  fCenter = 0;
  for (Int_t i = 0; i < 3; i++) SetLabel(0,i);
  SetX(0);
  SetY(0);
  SetZ(0);
  SetSigmaY2(0);
  SetSigmaZ2(0);
  SetVolumeId(0);
}

//_____________________________________________________________________________
Float_t AliTRDcluster::GetSumS() const
{
  //
  // Returns the total charge from a not unfolded cluster
  //

  Float_t sum = 0.0;
  for (Int_t i = 0; i < 7; i++) {
    sum += fSignals[i];
  }

  return sum;

}

//_____________________________________________________________________________
Float_t AliTRDcluster::GetXloc(Double_t t0, Double_t vd, Double_t *const q, Double_t *const xq, Double_t z)
{
//
// (Re)Calculate cluster position in the x direction in local chamber coordinates (with respect to the anode wire 
// position) using all available information from tracking.
// Input parameters:
//   t0 - calibration aware trigger delay [us]
//   vd - drift velocity in the region of the cluster [cm/us]
//   z  - distance to the anode wire [cm]. By default 0.2 !!
//   q & xq - array of charges and cluster positions from previous clusters in the tracklet [a.u.]
// Output values :
//   return x position of the cluster with respect to the 
//   anode wire using all tracking information
//
// The estimation of the radial position is based on calculating the drift time and the drift velocity at the point of 
// estimation. The drift time can be estimated according to the expression:
// BEGIN_LATEX
// t_{drift} = t_{bin} - t_{0} - t_{cause}(x) - t_{TC}(q_{i-1}, q_{i-2}, ...)
// END_LATEX
// where t_0 is the delay of the trigger signal. t_cause is the causality delay between ionisation electrons hitting 
// the anode and the registration of maximum signal by the electronics - it is due to the rising time of the TRF 
// convoluted with the diffusion width. t_TC is the residual charge from previous bins due to residual tails after tail 
// cancellation.
//
// The drift velocity is considered to vary linearly with the drift length (independent of the distance to the anode wire 
// in the z direction). Thus one can write the calculate iteratively the drift length from the expression:
// BEGIN_LATEX
// x = t_{drift}(x)*v_{drfit}(x)
// END_LATEX
//
// Authors
// Alex Bercuci <A.Bercuci@gsi.de>
//

  AliTRDCommonParam *cp = AliTRDCommonParam::Instance(); 
  Double_t fFreq = cp->GetSamplingFrequency();
  //drift time corresponding to the center of the time bin
  Double_t td = (fPadTime + .5)/fFreq; // [us] 
  if(td < t0+2.5) return 0.; // do not calculate radial posion of clusters in the amplification region

  // correction for t0
  td -= t0;

  Double_t x = vd*td, xold=0.;
  Float_t tc0   = 0.244, // TRF rising time 0.2us
          dtcdx = 0.009, // diffusion contribution to the rising time of the signal
          kTC   = 0.;    // tail cancellation residual
  while(TMath::Abs(x-xold)>1.e-3){ // convergence on 10um level 
    xold = x;
    Float_t tc  = tc0 - dtcdx*x; 
    Float_t tq  = 0.;
    if(q && xq){
      for(Int_t iq=0; iq<3; iq++) tq += q[iq]*TMath::Exp(-kTC*(x - xq[iq]));
    }
    Float_t vdcorr = x/cp->TimeStruct(vd, x+.5*AliTRDgeometry::CamHght(), z);
    x    = (td - tc - tq) * vdcorr;
  }
  return x;
}

//_____________________________________________________________________________
Float_t AliTRDcluster::GetYloc(Double_t s2, Double_t W, Double_t xd, Double_t wt, Double_t *const y1, Double_t *const y2)
{
//
// (Re)Calculate cluster position in the y direction in local chamber coordinates using all available information from tracking.
//
// Input parameters:
//   s2 - sigma of gaussian parameterization (see bellow for the exact parameterization)
//   W  - pad width
//   xd - drift length (with respect to the anode wire) [cm]
//   wt - omega*tau = tg(a_L)
// Output values :
//   y1 and y2 - partial positions based on 2 pads clusters
//   return y position of the cluster from all information
//
// Estimation of y coordinate is based on the gaussian approximation of the PRF. Thus one may
// calculate the y position knowing the signals q_i-1, q_i and q_i+1 in the 3 adiacent pads by:
// BEGIN_LATEX
// y = #frac{1}{w_{1}+w_{2}}#[]{w_{1}#(){y_{0}-#frac{W}{2}+#frac{s^{2}}{W}ln#frac{q_{i}}{q_{i-1}}}+w_{2}#(){y_{0}+ #frac{W}{2}+#frac{s^{2}}{W}ln#frac{q_{i+1}}{q_{i}}}}
// END_LATEX
// where W is the pad width, y_0 is the position of the center pad and s^2 is given by
// BEGIN_LATEX
// s^{2} = s^{2}_{0} + s^{2}_{diff} (x,B) + #frac{tg^{2}(#phi-#alpha_{L})*l^{2}}{12}
// END_LATEX
// with s_0 being the PRF for 0 drift and track incidence phi equal to the lorentz angle a_L and the diffusion term 
// being described by:
// BEGIN_LATEX
// s_{diff} (x,B) = #frac{D_{L}#sqrt{x}}{1+#(){#omega#tau}^{2}}
// END_LATEX
// with x being the drift length. The weights w_1 and w_2 are taken to be q_i-1^2 and q_i+1^2 respectively
// 
// Authors
// Alex Bercuci <A.Bercuci@gsi.de>
// Theodor Rascanu <trascanu@stud.uni-frankfurt.de>
//
  Float_t y0 = GetY()-W*fCenter;
  Double_t w1 = fSignals[2]*fSignals[2];
  Double_t w2 = fSignals[4]*fSignals[4];
  Float_t y1r = fSignals[2]>0 ? (y0 - .5*W + s2*TMath::Log(fSignals[3]/(Float_t)fSignals[2])/W) : 0.;
  Float_t y2r = fSignals[4]>0 ? (y0 + .5*W + s2*TMath::Log(fSignals[4]/(Float_t)fSignals[3])/W) : 0.;

  if(y1) (*y1) = y1r;
  if(y2) (*y2) = y2r;

  Double_t ld  = TMath::Max(xd - 0.*AliTRDgeometry::CamHght(), 0.);

  return (w1*y1r+w2*y2r)/(w1+w2) - ld*wt;
}

//_____________________________________________________________________________
Bool_t AliTRDcluster::IsEqual(const TObject *o) const
{
  //
  // Compare relevant information of this cluster with another one
  //
  
  const AliTRDcluster *inCluster = dynamic_cast<const AliTRDcluster*>(o);
  if (!o || !inCluster) return kFALSE;

  if ( AliCluster::GetX() != inCluster->GetX() ) return kFALSE;
  if ( AliCluster::GetY() != inCluster->GetY() ) return kFALSE;
  if ( AliCluster::GetZ() != inCluster->GetZ() ) return kFALSE;
  if ( fQ != inCluster->fQ ) return kFALSE;
  if ( fDetector != inCluster->fDetector ) return kFALSE;
  if ( fPadCol != inCluster->fPadCol ) return kFALSE;
  if ( fPadRow != inCluster->fPadRow ) return kFALSE;
  if ( fPadTime != inCluster->fPadTime ) return kFALSE;
  if ( fClusterMasking != inCluster->fClusterMasking ) return kFALSE;
  if ( IsInChamber() != inCluster->IsInChamber() ) return kFALSE;
  if ( IsShared() != inCluster->IsShared() ) return kFALSE;
  if ( IsUsed() != inCluster->IsUsed() ) return kFALSE;
  
  return kTRUE;
}

//_____________________________________________________________________________
void AliTRDcluster::Print(Option_t *o) const
{
  AliInfo(Form("Det[%3d] LTrC[%7.2f %7.2f %7.2f] Q[%f] Stat[in(%c) use(%c) sh(%c)]", 
    fDetector, GetX(), GetY(), GetZ(), fQ, 
    IsInChamber() ? 'y' : 'n', IsUsed() ? 'y' : 'n', IsShared() ? 'y' : 'n'));

  if(strcmp(o, "a")!=0) return;
  AliInfo(Form("LChC[c(%3d) r(%2d) t(%2d)] t-t0[%2d] Npad[%d] cen[%5.3f] mask[%d]", fPadCol, fPadRow, fPadTime, fLocalTimeBin, fNPads, fCenter, fClusterMasking)); 
  AliInfo(Form("Signals[%3d %3d %3d %3d %3d %3d %3d]", fSignals[0], fSignals[1], fSignals[2], fSignals[3], fSignals[4], fSignals[5], fSignals[6]));
}


//_____________________________________________________________________________
void AliTRDcluster::SetPadMaskedPosition(UChar_t position)
{
  //
  // store the pad corruption position code
  // 
  // Code: 1 = left cluster
  //       2 = middle cluster;
  //       4 = right cluster
  //
  for(Int_t ipos = 0; ipos < 3; ipos++)
    if(TESTBIT(position, ipos))
      SETBIT(fClusterMasking, ipos);
}

//_____________________________________________________________________________
void AliTRDcluster::SetPadMaskedStatus(UChar_t status)
{
  //
  // store the status of the corrupted pad
  //
  // Code: 2 = noisy
  //       4 = Bridged Left
  //       8 = Bridged Right
  //      32 = Not Connected
  for(Int_t ipos = 0; ipos < 5; ipos++)
    if(TESTBIT(status, ipos))
      SETBIT(fClusterMasking, ipos + 3);
}
