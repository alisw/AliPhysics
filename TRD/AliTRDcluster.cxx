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

const Int_t AliTRDcluster::fgkNlut = 128;
Double_t *AliTRDcluster::fgLUT = 0x0;

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
AliTRDcluster::AliTRDcluster(Int_t det, UChar_t col, UChar_t row, UChar_t time, const Short_t *sig, UShort_t vid) 
  :AliCluster() 
  ,fPadCol(col)
  ,fPadRow(row)
  ,fPadTime(time)
  ,fLocalTimeBin(0)
  ,fNPads(0)
  ,fClusterMasking(0)
  ,fDetector(det)
  ,fQ(0.)
  ,fCenter(0.)
{ 
  //
  // Constructor for self constructing cluster. In this approach the information is inserted gradualy into the 
  // cluster and all dependencies are (re)calculated inside the cluster itself.
  //
  // A.Bercuci <A.Bercuci@gsi.de>

  memcpy(&fSignals, sig, 7*sizeof(Short_t));
  fQ = fSignals[2]+fSignals[3]+fSignals[4];
  SetVolumeId(vid);
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

//___________________________________________________________________________
Double_t AliTRDcluster::GetSX(Int_t tb, Double_t z)
{
  if(tb<1 || tb>=24) return 10.; // return huge [10cm]
  const Double_t sx[24][10]={
    {0.000e+00, 9.352e-01, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 2.309e+00},
    {8.387e-02, 8.718e-02, 8.816e-02, 9.444e-02, 9.993e-02, 1.083e-01, 1.161e-01, 1.280e-01, 1.417e-01, 1.406e-01},
    {1.097e-01, 1.105e-01, 1.127e-01, 1.151e-01, 1.186e-01, 1.223e-01, 1.272e-01, 1.323e-01, 1.389e-01, 1.490e-01},
    {1.407e-01, 1.404e-01, 1.414e-01, 1.430e-01, 1.429e-01, 1.449e-01, 1.476e-01, 1.494e-01, 1.515e-01, 1.589e-01},
    {1.681e-01, 1.679e-01, 1.666e-01, 1.657e-01, 1.656e-01, 1.649e-01, 1.652e-01, 1.662e-01, 1.671e-01, 1.694e-01},
    {1.745e-01, 1.737e-01, 1.707e-01, 1.690e-01, 1.643e-01, 1.610e-01, 1.612e-01, 1.628e-01, 1.638e-01, 1.659e-01},
    {1.583e-01, 1.558e-01, 1.535e-01, 1.488e-01, 1.445e-01, 1.419e-01, 1.428e-01, 1.451e-01, 1.462e-01, 1.494e-01},
    {1.414e-01, 1.391e-01, 1.368e-01, 1.300e-01, 1.256e-01, 1.259e-01, 1.285e-01, 1.326e-01, 1.358e-01, 1.406e-01},
    {1.307e-01, 1.289e-01, 1.261e-01, 1.216e-01, 1.193e-01, 1.165e-01, 1.201e-01, 1.241e-01, 1.274e-01, 1.344e-01},
    {1.251e-01, 1.227e-01, 1.208e-01, 1.155e-01, 1.110e-01, 1.116e-01, 1.133e-01, 1.187e-01, 1.229e-01, 1.308e-01},
    {1.234e-01, 1.209e-01, 1.175e-01, 1.127e-01, 1.094e-01, 1.093e-01, 1.109e-01, 1.155e-01, 1.210e-01, 1.275e-01},
    {1.215e-01, 1.187e-01, 1.156e-01, 1.108e-01, 1.070e-01, 1.065e-01, 1.090e-01, 1.134e-01, 1.196e-01, 1.251e-01},
    {1.202e-01, 1.180e-01, 1.151e-01, 1.108e-01, 1.070e-01, 1.058e-01, 1.089e-01, 1.127e-01, 1.183e-01, 1.256e-01},
    {1.207e-01, 1.176e-01, 1.142e-01, 1.109e-01, 1.072e-01, 1.069e-01, 1.088e-01, 1.122e-01, 1.182e-01, 1.252e-01},
    {1.213e-01, 1.182e-01, 1.156e-01, 1.102e-01, 1.076e-01, 1.063e-01, 1.091e-01, 1.132e-01, 1.181e-01, 1.243e-01},
    {1.205e-01, 1.180e-01, 1.150e-01, 1.104e-01, 1.072e-01, 1.063e-01, 1.083e-01, 1.132e-01, 1.183e-01, 1.243e-01},
    {1.212e-01, 1.195e-01, 1.135e-01, 1.107e-01, 1.070e-01, 1.065e-01, 1.097e-01, 1.126e-01, 1.185e-01, 1.238e-01},
    {1.201e-01, 1.184e-01, 1.155e-01, 1.111e-01, 1.088e-01, 1.075e-01, 1.089e-01, 1.131e-01, 1.189e-01, 1.237e-01},
    {1.197e-01, 1.186e-01, 1.147e-01, 1.113e-01, 1.085e-01, 1.077e-01, 1.105e-01, 1.137e-01, 1.188e-01, 1.245e-01},
    {1.213e-01, 1.194e-01, 1.154e-01, 1.114e-01, 1.091e-01, 1.082e-01, 1.098e-01, 1.140e-01, 1.194e-01, 1.247e-01},
    {1.210e-01, 1.189e-01, 1.155e-01, 1.119e-01, 1.088e-01, 1.080e-01, 1.105e-01, 1.141e-01, 1.195e-01, 1.244e-01},
    {1.196e-01, 1.189e-01, 1.145e-01, 1.105e-01, 1.095e-01, 1.083e-01, 1.087e-01, 1.121e-01, 1.173e-01, 1.208e-01},
    {1.123e-01, 1.129e-01, 1.108e-01, 1.110e-01, 1.080e-01, 1.065e-01, 1.056e-01, 1.066e-01, 1.071e-01, 1.095e-01},
    {1.136e-01, 1.135e-01, 1.130e-01, 1.122e-01, 1.113e-01, 1.071e-01, 1.041e-01, 1.025e-01, 1.014e-01, 9.973e-02}
  };
  if(z>=0. && z<.25) return sx[tb][Int_t(z/.025)];
  
  Double_t m = 1.e-8; for(Int_t id=10; id--;) if(sx[tb][id]>m) m=sx[tb][id];
  return m;
}

//___________________________________________________________________________
Double_t AliTRDcluster::GetSY(Int_t tb, Double_t z)
{
  if(tb<1 || tb>=24) return 10.; // return huge [10cm]
  const Double_t sy[24][10]={
    {0.000e+00, 2.610e-01, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 4.680e-01},
    {3.019e-02, 3.036e-02, 3.131e-02, 3.203e-02, 3.294e-02, 3.407e-02, 3.555e-02, 3.682e-02, 3.766e-02, 3.824e-02},
    {1.773e-02, 1.778e-02, 1.772e-02, 1.790e-02, 1.807e-02, 1.833e-02, 1.873e-02, 1.905e-02, 1.958e-02, 2.029e-02},
    {1.774e-02, 1.772e-02, 1.746e-02, 1.738e-02, 1.756e-02, 1.756e-02, 1.739e-02, 1.720e-02, 1.743e-02, 1.769e-02},
    {2.064e-02, 2.078e-02, 2.069e-02, 2.060e-02, 2.033e-02, 2.024e-02, 2.022e-02, 1.961e-02, 1.922e-02, 1.901e-02},
    {2.382e-02, 2.379e-02, 2.371e-02, 2.333e-02, 2.318e-02, 2.285e-02, 2.255e-02, 2.244e-02, 2.174e-02, 2.132e-02},
    {2.615e-02, 2.589e-02, 2.539e-02, 2.493e-02, 2.420e-02, 2.396e-02, 2.362e-02, 2.342e-02, 2.321e-02, 2.330e-02},
    {2.640e-02, 2.638e-02, 2.577e-02, 2.548e-02, 2.477e-02, 2.436e-02, 2.416e-02, 2.401e-02, 2.399e-02, 2.402e-02},
    {2.647e-02, 2.632e-02, 2.587e-02, 2.546e-02, 2.465e-02, 2.447e-02, 2.429e-02, 2.415e-02, 2.429e-02, 2.475e-02},
    {2.657e-02, 2.637e-02, 2.580e-02, 2.525e-02, 2.492e-02, 2.441e-02, 2.446e-02, 2.441e-02, 2.478e-02, 2.491e-02},
    {2.640e-02, 2.608e-02, 2.583e-02, 2.539e-02, 2.478e-02, 2.440e-02, 2.456e-02, 2.464e-02, 2.486e-02, 2.533e-02},
    {2.636e-02, 2.630e-02, 2.584e-02, 2.542e-02, 2.483e-02, 2.451e-02, 2.449e-02, 2.467e-02, 2.496e-02, 2.554e-02},
    {2.634e-02, 2.629e-02, 2.583e-02, 2.526e-02, 2.480e-02, 2.460e-02, 2.458e-02, 2.472e-02, 2.518e-02, 2.549e-02},
    {2.629e-02, 2.621e-02, 2.581e-02, 2.527e-02, 2.480e-02, 2.458e-02, 2.451e-02, 2.485e-02, 2.516e-02, 2.547e-02},
    {2.629e-02, 2.607e-02, 2.573e-02, 2.543e-02, 2.485e-02, 2.464e-02, 2.452e-02, 2.476e-02, 2.505e-02, 2.550e-02},
    {2.635e-02, 2.613e-02, 2.578e-02, 2.523e-02, 2.491e-02, 2.465e-02, 2.470e-02, 2.467e-02, 2.515e-02, 2.564e-02},
    {2.613e-02, 2.602e-02, 2.587e-02, 2.526e-02, 2.507e-02, 2.482e-02, 2.456e-02, 2.486e-02, 2.509e-02, 2.572e-02},
    {2.620e-02, 2.599e-02, 2.563e-02, 2.528e-02, 2.484e-02, 2.462e-02, 2.464e-02, 2.476e-02, 2.513e-02, 2.571e-02},
    {2.634e-02, 2.596e-02, 2.565e-02, 2.519e-02, 2.497e-02, 2.457e-02, 2.450e-02, 2.481e-02, 2.511e-02, 2.540e-02},
    {2.593e-02, 2.589e-02, 2.563e-02, 2.511e-02, 2.472e-02, 2.453e-02, 2.452e-02, 2.474e-02, 2.501e-02, 2.543e-02},
    {2.576e-02, 2.582e-02, 2.526e-02, 2.505e-02, 2.462e-02, 2.446e-02, 2.445e-02, 2.466e-02, 2.486e-02, 2.510e-02},
    {2.571e-02, 2.549e-02, 2.533e-02, 2.501e-02, 2.453e-02, 2.443e-02, 2.445e-02, 2.450e-02, 2.448e-02, 2.469e-02},
    {2.812e-02, 2.786e-02, 2.776e-02, 2.723e-02, 2.695e-02, 2.650e-02, 2.642e-02, 2.617e-02, 2.612e-02, 2.610e-02},
    {3.251e-02, 3.267e-02, 3.223e-02, 3.183e-02, 3.125e-02, 3.106e-02, 3.067e-02, 3.010e-02, 2.936e-02, 2.927e-02}
  };
  if(z>=0. && z<.25) return sy[tb][Int_t(z/.025)];

  Double_t m = 1.e-8; for(Int_t id=10; id--;) if(sy[tb][id]>m) m=sy[tb][id];

  return m;
}

//___________________________________________________________________________
Double_t AliTRDcluster::GetXcorr(Int_t tb, Double_t z)
{
  // drift length correction [cm]
  // TODO to be parametrized in term of drift velocity
  // A.Bercuci (Mar 28 2009)

  if(tb<0 || tb>=24) return 0.;
  const Int_t nd = 5;
  const Double_t dx[24][nd]={
    {+1.747e-01,+3.195e-01,+1.641e-01,+1.607e-01,+6.002e-01},
    {+5.468e-02,+5.760e-02,+6.365e-02,+8.003e-02,+1.067e-01},
    {-6.327e-02,-6.339e-02,-6.423e-02,-6.900e-02,-7.949e-02},
    {-1.417e-01,-1.424e-01,-1.450e-01,-1.465e-01,-1.514e-01},
    {-1.637e-01,-1.619e-01,-1.622e-01,-1.613e-01,-1.648e-01},
    {-1.386e-01,-1.334e-01,-1.261e-01,-1.276e-01,-1.314e-01},
    {-8.799e-02,-8.299e-02,-7.861e-02,-8.038e-02,-8.436e-02},
    {-5.139e-02,-4.849e-02,-4.641e-02,-4.965e-02,-5.286e-02},
    {-2.927e-02,-2.773e-02,-2.807e-02,-3.021e-02,-3.378e-02},
    {-1.380e-02,-1.229e-02,-1.335e-02,-1.547e-02,-1.984e-02},
    {-4.168e-03,-4.601e-03,-5.462e-03,-8.164e-03,-1.035e-02},
    {+2.044e-03,+1.889e-03,+9.603e-04,-1.342e-03,-3.736e-03},
    {+3.568e-03,+3.581e-03,+2.391e-03,+2.942e-05,-1.585e-03},
    {+4.403e-03,+4.571e-03,+3.509e-03,+8.703e-04,-1.425e-03},
    {+4.941e-03,+4.808e-03,+3.284e-03,+1.105e-03,-1.208e-03},
    {+5.124e-03,+5.022e-03,+4.305e-03,+2.023e-03,-1.145e-03},
    {+4.882e-03,+4.008e-03,+3.408e-03,+7.886e-04,-1.356e-03},
    {+3.852e-03,+3.539e-03,+2.057e-03,+1.670e-04,-1.993e-03},
    {+2.154e-03,+2.111e-03,+5.723e-04,-1.254e-03,-3.256e-03},
    {+1.755e-03,+2.101e-03,+9.516e-04,-1.649e-03,-3.394e-03},
    {+1.617e-03,+1.662e-03,+4.169e-04,-9.843e-04,-4.309e-03},
    {-9.204e-03,-9.069e-03,-1.182e-02,-1.458e-02,-1.880e-02},
    {-6.727e-02,-6.820e-02,-6.804e-02,-7.134e-02,-7.615e-02},
    {-1.802e-01,-1.733e-01,-1.633e-01,-1.601e-01,-1.632e-01}
  };
//   const Double_t dx[24][nd]={
//     {+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00},
//     {-2.763e-04,-2.380e-04,-6.286e-04,-9.424e-04,+1.046e-03,+1.932e-03,+1.620e-03,+1.951e-03,-1.321e-03,-1.115e-03},
//     {-1.825e-03,-9.245e-04,-1.012e-03,-8.215e-04,+2.703e-05,+1.403e-03,+2.340e-03,+2.577e-03,+2.017e-03,+8.006e-04},
//     {-3.070e-03,-8.563e-04,-1.257e-03,+8.491e-05,+4.503e-04,-2.467e-05,-1.793e-04,+5.085e-04,+1.321e-03,+4.056e-04},
//     {-3.637e-03,-2.857e-03,-3.098e-03,-2.304e-03,-1.467e-03,-1.755e-03,+4.585e-04,+2.757e-03,+3.184e-03,+3.525e-03},
//     {-9.884e-03,-7.695e-03,-7.290e-03,-3.990e-03,-9.982e-04,+2.226e-03,+3.375e-03,+6.419e-03,+7.209e-03,+6.891e-03},
//     {-6.844e-03,-5.807e-03,-4.012e-03,-1.566e-03,+5.035e-04,+2.024e-03,+3.225e-03,+3.918e-03,+5.942e-03,+6.024e-03},
//     {-2.628e-03,-2.201e-03,-4.562e-04,+9.832e-04,+3.411e-03,+2.062e-03,+1.526e-03,+9.350e-04,+8.842e-04,+1.007e-03},
//     {+6.603e-04,+1.545e-03,+1.681e-03,+1.918e-03,+2.165e-03,+1.825e-03,+1.691e-03,-1.923e-04,+1.835e-04,-1.284e-03},
//     {+1.895e-03,+1.586e-03,+2.000e-03,+3.537e-03,+2.526e-03,+1.316e-03,+8.229e-04,-7.671e-05,-2.175e-03,-3.528e-03},
//     {+2.927e-03,+3.369e-03,+3.603e-03,+2.675e-03,+2.737e-03,+1.133e-03,+4.318e-04,-1.215e-03,-2.443e-03,-3.116e-03},
//     {+3.072e-03,+3.564e-03,+3.612e-03,+3.149e-03,+2.768e-03,+1.186e-03,+3.083e-04,-1.447e-03,-2.480e-03,-3.263e-03},
//     {+2.697e-03,+3.565e-03,+3.759e-03,+2.855e-03,+2.909e-03,+6.564e-04,-5.224e-04,-3.309e-04,-1.636e-03,-3.739e-03},
//     {+3.633e-03,+3.232e-03,+3.727e-03,+3.024e-03,+3.365e-03,+1.598e-03,-6.903e-04,-1.039e-03,-3.176e-03,-4.472e-03},
//     {+2.999e-03,+3.942e-03,+3.322e-03,+3.162e-03,+1.978e-03,+1.657e-03,-4.760e-04,-8.343e-04,-2.346e-03,-3.281e-03},
//     {+3.734e-03,+3.098e-03,+3.435e-03,+2.512e-03,+2.651e-03,+1.745e-03,+9.424e-04,-1.404e-03,-3.177e-03,-4.444e-03},
//     {+3.204e-03,+4.003e-03,+3.068e-03,+2.697e-03,+3.187e-03,+3.878e-04,-1.124e-04,-1.855e-03,-2.584e-03,-3.807e-03},
//     {+2.653e-03,+3.631e-03,+2.327e-03,+3.460e-03,+1.810e-03,+1.244e-03,-3.651e-04,-2.664e-04,-2.307e-03,-3.642e-03},
//     {+2.538e-03,+3.208e-03,+2.390e-03,+3.519e-03,+1.763e-03,+1.330e-04,+1.669e-04,-1.422e-03,-1.685e-03,-3.519e-03},
//     {+2.605e-03,+2.465e-03,+2.771e-03,+2.966e-03,+2.361e-03,+6.029e-04,-4.435e-04,-1.876e-03,-1.694e-03,-3.757e-03},
//     {+2.866e-03,+3.315e-03,+3.146e-03,+2.117e-03,+1.933e-03,+9.339e-04,+9.556e-04,-1.314e-03,-3.615e-03,-3.558e-03},
//     {+4.002e-03,+3.543e-03,+3.631e-03,+4.127e-03,+1.919e-03,-2.852e-04,-9.484e-04,-2.060e-03,-4.477e-03,-5.491e-03},
//     {+6.029e-03,+5.147e-03,+4.286e-03,+2.215e-03,+9.240e-04,-1.554e-03,-2.366e-03,-3.635e-03,-5.372e-03,-6.467e-03},
//     {+3.941e-03,+3.995e-03,+5.638e-04,-3.332e-04,-2.539e-03,-3.764e-03,-3.647e-03,-4.900e-03,-5.414e-03,-5.202e-03}
//   };
  if(z>=0. && z<.25) return dx[tb][Int_t(z/.025)];

  Double_t m = 0.; for(Int_t id=nd; id--;) m+=dx[tb][id];
  return m/nd;
}

//___________________________________________________________________________
Double_t AliTRDcluster::GetYcorr(Int_t ly, Float_t y)
{
// PRF correction TODO to be replaced by the gaussian 
// approximation with full error parametrization
  const Float_t cy[AliTRDgeometry::kNlayer][3] = {
    { 4.014e-04, 8.605e-03, -6.880e+00},
    {-3.061e-04, 9.663e-03, -6.789e+00},
    { 1.124e-03, 1.105e-02, -6.825e+00},
    {-1.527e-03, 1.231e-02, -6.777e+00},
    { 2.150e-03, 1.387e-02, -6.783e+00},
    {-1.296e-03, 1.486e-02, -6.825e+00}
  }; 

  return cy[ly][0] + cy[ly][1] * TMath::Sin(cy[ly][2] * y);
}

//_____________________________________________________________________________
Float_t AliTRDcluster::GetXloc(Double_t t0, Double_t vd, Double_t *const /*q*/, Double_t *const /*xq*/, Double_t /*z*/)
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

  // calculate t0 corrected time bin
  Double_t td = fPadTime - t0;
  fLocalTimeBin = TMath::Nint(td);
  //drift time corresponding to the center of the time bin
  td = (td + .5)/fFreq; // [us] 
  // correction for t0
  td -= t0;
  // calculate radial posion of clusters in the drift region
  if(td < .2) return 0.;
  // TRF rising time (fitted)
  // It should be absorbed by the t0. For the moment t0 is 0 for simulations.
  // A.Bercuci (Mar 26 2009)
  td -= 0.189;

  // apply fitted correction 
  Float_t x = td*vd + GetXcorr(fLocalTimeBin);
  if(x>0.&&x<.5*AliTRDgeometry::CamHght()+AliTRDgeometry::CdrHght()) SetInChamber();

  return x;

/*
  // invert drift time function
  Double_t xM= AliTRDgeometry::CamHght()+AliTRDgeometry::CdrHght(),
           x = vd*td + .5*AliTRDgeometry::CamHght(), 
           t = cp->TimeStruct(vd, x, z), dx1=0.,dx2;
  while(TMath::Abs(td-t)>1.e-4){ // convergence on 100ps
    dx2 = vd*(td-t);
    if(TMath::Abs(TMath::Abs(dx2)-TMath::Abs(dx1))<1.e-6){
      x+=.5*dx2;
      break;
    } else x+=dx2;

    if(x<0. || x>xM) return 0.;
    t = cp->TimeStruct(vd, x, z);
    dx1 = dx2;
  }

  return x-.5*AliTRDgeometry::CamHght();
*/
}

//_____________________________________________________________________________
Float_t AliTRDcluster::GetYloc(Double_t y0, Double_t s2, Double_t W, Double_t *const y1, Double_t *const y2)
{

  //printf("  s[%3d %3d %3d] w[%f %f] yr[%f %f]\n", fSignals[2], fSignals[3], fSignals[4], w1/(w1+w2), w2/(w1+w2), y1r*W, y2r*W);
  if(IsRPhiMethod(kCOG)) GetDYcog();
  else if(IsRPhiMethod(kLUT)) GetDYlut();
  if(IsRPhiMethod(kGAUS)) GetDYgauss(s2/W/W, y1, y2);

  if(y1) (*y1)*=W;
  if(y2) (*y2)*=W;

  return y0+fCenter*W+(IsRPhiMethod(kLUT)?GetYcorr(AliTRDgeometry::GetLayer(fDetector), fCenter):0.);
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
  AliInfo(Form("Det[%3d] LTrC[%+6.2f %+6.2f %+6.2f] Q[%5.1f] FLAG[in(%c) use(%c) sh(%c)] Y[%s]", 
    fDetector, GetX(), GetY(), GetZ(), fQ, 
    IsInChamber() ? 'y' : 'n', 
    IsUsed() ? 'y' : 'n', 
    IsShared() ? 'y' : 'n',
    IsRPhiMethod(kGAUS)?"GAUS":(IsRPhiMethod(kLUT)?"LUT":"COG")
  ));

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

//___________________________________________________________________________
Float_t AliTRDcluster::GetDYcog(Double_t *const, Double_t *const)
{
//
// Get COG position
// Used for clusters with more than 3 pads - where LUT not applicable
//
  Double_t sum = fSignals[1]
                +fSignals[2]
                +fSignals[3] 
                +fSignals[4]
                +fSignals[5];

  // ???????????? CBL
  // Go to 3 pad COG ????
  // ???????????? CBL
  fCenter = (0.0 * (-fSignals[1] + fSignals[5])
                      + (-fSignals[2] + fSignals[4])) / sum;

  return fCenter;
}

//___________________________________________________________________________
Float_t AliTRDcluster::GetDYlut(Double_t *const, Double_t *const)
{
  //
  // Calculates the cluster position using the lookup table.
  // Method provided by Bogdan Vulpescu.
  //

  if(!fgLUT) FillLUT();

  Double_t ampL = fSignals[2],
           ampC = fSignals[3],
           ampR = fSignals[4];
  Int_t ilayer = AliTRDgeometry::GetLayer(fDetector);

  Double_t x    = 0.0;
  Double_t xmin, xmax, xwid;

  Int_t    side = 0;
  Int_t    ix;

  Double_t xMin[AliTRDgeometry::kNlayer] = { 
    0.006492, 0.006377, 0.006258, 0.006144, 0.006030, 0.005980 
  };
  Double_t xMax[AliTRDgeometry::kNlayer] = { 
    0.960351, 0.965870, 0.970445, 0.974352, 0.977667, 0.996101 
  };

  if      (ampL > ampR) {
    x    = (ampL - ampR) / ampC;
    side = -1;
  } 
  else if (ampL < ampR) {
    x    = (ampR - ampL) / ampC;
    side = +1;
  }

  if (ampL != ampR) {

    xmin = xMin[ilayer] + 0.000005;
    xmax = xMax[ilayer] - 0.000005;
    xwid = (xmax - xmin) / 127.0;

    if      (x < xmin) fCenter = 0.0000;
    else if (x > xmax) fCenter = side * 0.5000;
    else {
      ix      = (Int_t) ((x - xmin) / xwid);
      fCenter = side * fgLUT[ilayer*fgkNlut+ix];
    }
  } else fCenter = 0.0;

  return fCenter;
}

//___________________________________________________________________________
Float_t AliTRDcluster::GetDYgauss(Double_t s2w, Double_t *const y1, Double_t *const y2)
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
  Double_t w1 = fSignals[2]*fSignals[2];
  Double_t w2 = fSignals[4]*fSignals[4];
  Double_t w = w1+w2;
  if(w<1.){
    AliError("Missing side signals for cluster.");
    Print("a");
    return 0.;
  }  

  //Double_t s2w = s2/W/W;
  Float_t y1r  = fSignals[2]>0 ? (-0.5 + s2w*TMath::Log(fSignals[3]/(Float_t)fSignals[2])) : 0.;
  Float_t y2r  = fSignals[4]>0 ? (0.5 + s2w*TMath::Log(fSignals[4]/(Float_t)fSignals[3])) : 0.;

  if(y1) (*y1) = y1r;
  if(y2) (*y2) = y2r;

  return fCenter      = (w1*y1r+w2*y2r)/w;
}



//_____________________________________________________________________________
void AliTRDcluster::FillLUT()
{
  //
  // Create the LUT
  //

  // The lookup table from Bogdan
  Float_t lut[AliTRDgeometry::kNlayer][fgkNlut] = {  
    {
      0.0070, 0.0150, 0.0224, 0.0298, 0.0374, 0.0454, 0.0533, 0.0611, 
      0.0684, 0.0755, 0.0827, 0.0900, 0.0975, 0.1049, 0.1120, 0.1187, 
      0.1253, 0.1318, 0.1385, 0.1453, 0.1519, 0.1584, 0.1646, 0.1704, 
      0.1762, 0.1821, 0.1879, 0.1938, 0.1996, 0.2053, 0.2108, 0.2160, 
      0.2210, 0.2260, 0.2310, 0.2361, 0.2411, 0.2461, 0.2509, 0.2557, 
      0.2602, 0.2646, 0.2689, 0.2732, 0.2774, 0.2816, 0.2859, 0.2901, 
      0.2942, 0.2983, 0.3022, 0.3061, 0.3099, 0.3136, 0.3172, 0.3207, 
      0.3242, 0.3278, 0.3312, 0.3347, 0.3382, 0.3416, 0.3450, 0.3483, 
      0.3515, 0.3547, 0.3579, 0.3609, 0.3639, 0.3669, 0.3698, 0.3727, 
      0.3756, 0.3785, 0.3813, 0.3842, 0.3870, 0.3898, 0.3926, 0.3952, 
      0.3979, 0.4005, 0.4032, 0.4057, 0.4082, 0.4108, 0.4132, 0.4157, 
      0.4181, 0.4205, 0.4228, 0.4252, 0.4275, 0.4299, 0.4322, 0.4345, 
      0.4367, 0.4390, 0.4412, 0.4434, 0.4456, 0.4478, 0.4499, 0.4520, 
      0.4541, 0.4562, 0.4583, 0.4603, 0.4623, 0.4643, 0.4663, 0.4683, 
      0.4702, 0.4722, 0.4741, 0.4758, 0.4774, 0.4790, 0.4805, 0.4824, 
      0.4844, 0.4863, 0.4883, 0.4902, 0.4921, 0.4940, 0.4959, 0.4978 
    },
    {
      0.0072, 0.0156, 0.0235, 0.0313, 0.0394, 0.0478, 0.0561, 0.0642, 
      0.0718, 0.0792, 0.0868, 0.0947, 0.1025, 0.1101, 0.1172, 0.1241, 
      0.1309, 0.1378, 0.1449, 0.1518, 0.1586, 0.1650, 0.1710, 0.1770, 
      0.1830, 0.1891, 0.1952, 0.2011, 0.2070, 0.2125, 0.2177, 0.2229, 
      0.2280, 0.2332, 0.2383, 0.2435, 0.2484, 0.2533, 0.2581, 0.2627, 
      0.2670, 0.2714, 0.2757, 0.2799, 0.2842, 0.2884, 0.2927, 0.2968, 
      0.3008, 0.3048, 0.3086, 0.3123, 0.3159, 0.3195, 0.3231, 0.3266, 
      0.3301, 0.3335, 0.3370, 0.3404, 0.3438, 0.3471, 0.3504, 0.3536, 
      0.3567, 0.3598, 0.3628, 0.3657, 0.3686, 0.3715, 0.3744, 0.3772, 
      0.3800, 0.3828, 0.3856, 0.3884, 0.3911, 0.3938, 0.3965, 0.3991, 
      0.4016, 0.4042, 0.4067, 0.4092, 0.4116, 0.4140, 0.4164, 0.4187, 
      0.4211, 0.4234, 0.4257, 0.4280, 0.4302, 0.4325, 0.4347, 0.4369, 
      0.4391, 0.4413, 0.4434, 0.4456, 0.4477, 0.4497, 0.4518, 0.4538, 
      0.4558, 0.4578, 0.4598, 0.4618, 0.4637, 0.4656, 0.4675, 0.4694, 
      0.4713, 0.4732, 0.4750, 0.4766, 0.4781, 0.4797, 0.4813, 0.4832, 
      0.4851, 0.4870, 0.4888, 0.4906, 0.4925, 0.4942, 0.4960, 0.4978
    },
    {
      0.0075, 0.0163, 0.0246, 0.0328, 0.0415, 0.0504, 0.0592, 0.0674, 
      0.0753, 0.0832, 0.0914, 0.0996, 0.1077, 0.1154, 0.1225, 0.1296, 
      0.1369, 0.1442, 0.1515, 0.1585, 0.1652, 0.1714, 0.1776, 0.1839, 
      0.1902, 0.1965, 0.2025, 0.2085, 0.2141, 0.2194, 0.2247, 0.2299, 
      0.2352, 0.2405, 0.2457, 0.2507, 0.2557, 0.2604, 0.2649, 0.2693, 
      0.2737, 0.2780, 0.2823, 0.2867, 0.2909, 0.2951, 0.2992, 0.3033, 
      0.3072, 0.3110, 0.3146, 0.3182, 0.3218, 0.3253, 0.3288, 0.3323, 
      0.3357, 0.3392, 0.3426, 0.3459, 0.3492, 0.3524, 0.3555, 0.3586, 
      0.3616, 0.3645, 0.3674, 0.3703, 0.3731, 0.3759, 0.3787, 0.3815, 
      0.3843, 0.3870, 0.3897, 0.3925, 0.3950, 0.3976, 0.4002, 0.4027, 
      0.4052, 0.4076, 0.4101, 0.4124, 0.4148, 0.4171, 0.4194, 0.4217, 
      0.4239, 0.4262, 0.4284, 0.4306, 0.4328, 0.4350, 0.4371, 0.4393, 
      0.4414, 0.4435, 0.4455, 0.4476, 0.4496, 0.4516, 0.4536, 0.4555, 
      0.4575, 0.4594, 0.4613, 0.4632, 0.4650, 0.4669, 0.4687, 0.4705, 
      0.4723, 0.4741, 0.4758, 0.4773, 0.4789, 0.4804, 0.4821, 0.4839, 
      0.4857, 0.4875, 0.4893, 0.4910, 0.4928, 0.4945, 0.4961, 0.4978
    },
    {
      0.0078, 0.0171, 0.0258, 0.0345, 0.0438, 0.0532, 0.0624, 0.0708, 
      0.0791, 0.0875, 0.0962, 0.1048, 0.1130, 0.1206, 0.1281, 0.1356, 
      0.1432, 0.1508, 0.1582, 0.1651, 0.1716, 0.1780, 0.1845, 0.1910, 
      0.1975, 0.2038, 0.2099, 0.2155, 0.2210, 0.2263, 0.2317, 0.2371, 
      0.2425, 0.2477, 0.2528, 0.2578, 0.2626, 0.2671, 0.2715, 0.2759, 
      0.2803, 0.2846, 0.2890, 0.2933, 0.2975, 0.3016, 0.3056, 0.3095, 
      0.3132, 0.3168, 0.3204, 0.3239, 0.3274, 0.3309, 0.3344, 0.3378, 
      0.3412, 0.3446, 0.3479, 0.3511, 0.3543, 0.3574, 0.3603, 0.3633, 
      0.3662, 0.3690, 0.3718, 0.3747, 0.3774, 0.3802, 0.3829, 0.3857, 
      0.3883, 0.3910, 0.3936, 0.3962, 0.3987, 0.4012, 0.4037, 0.4061, 
      0.4085, 0.4109, 0.4132, 0.4155, 0.4177, 0.4200, 0.4222, 0.4244, 
      0.4266, 0.4288, 0.4309, 0.4331, 0.4352, 0.4373, 0.4394, 0.4414, 
      0.4435, 0.4455, 0.4475, 0.4494, 0.4514, 0.4533, 0.4552, 0.4571, 
      0.4590, 0.4608, 0.4626, 0.4645, 0.4662, 0.4680, 0.4698, 0.4715, 
      0.4733, 0.4750, 0.4766, 0.4781, 0.4796, 0.4812, 0.4829, 0.4846, 
      0.4863, 0.4880, 0.4897, 0.4914, 0.4930, 0.4946, 0.4963, 0.4979
    },
    {
      0.0081, 0.0178, 0.0270, 0.0364, 0.0463, 0.0562, 0.0656, 0.0744, 
      0.0831, 0.0921, 0.1013, 0.1102, 0.1183, 0.1261, 0.1339, 0.1419, 
      0.1499, 0.1576, 0.1648, 0.1715, 0.1782, 0.1849, 0.1917, 0.1984, 
      0.2048, 0.2110, 0.2167, 0.2223, 0.2278, 0.2333, 0.2389, 0.2444, 
      0.2497, 0.2548, 0.2598, 0.2645, 0.2691, 0.2735, 0.2780, 0.2824, 
      0.2868, 0.2912, 0.2955, 0.2997, 0.3038, 0.3078, 0.3116, 0.3152, 
      0.3188, 0.3224, 0.3259, 0.3294, 0.3329, 0.3364, 0.3398, 0.3432, 
      0.3465, 0.3497, 0.3529, 0.3561, 0.3591, 0.3620, 0.3649, 0.3677, 
      0.3705, 0.3733, 0.3761, 0.3788, 0.3816, 0.3843, 0.3869, 0.3896, 
      0.3922, 0.3948, 0.3973, 0.3998, 0.4022, 0.4047, 0.4070, 0.4094, 
      0.4117, 0.4139, 0.4162, 0.4184, 0.4206, 0.4227, 0.4249, 0.4270, 
      0.4291, 0.4313, 0.4334, 0.4354, 0.4375, 0.4395, 0.4415, 0.4435, 
      0.4455, 0.4474, 0.4493, 0.4512, 0.4531, 0.4550, 0.4568, 0.4586, 
      0.4604, 0.4622, 0.4639, 0.4657, 0.4674, 0.4691, 0.4708, 0.4725, 
      0.4742, 0.4758, 0.4773, 0.4788, 0.4803, 0.4819, 0.4836, 0.4852, 
      0.4869, 0.4885, 0.4901, 0.4917, 0.4933, 0.4948, 0.4964, 0.4979
    },
    {
      0.0085, 0.0189, 0.0288, 0.0389, 0.0497, 0.0603, 0.0699, 0.0792, 
      0.0887, 0.0985, 0.1082, 0.1170, 0.1253, 0.1336, 0.1421, 0.1505, 
      0.1587, 0.1662, 0.1733, 0.1803, 0.1874, 0.1945, 0.2014, 0.2081, 
      0.2143, 0.2201, 0.2259, 0.2316, 0.2374, 0.2431, 0.2487, 0.2541, 
      0.2593, 0.2642, 0.2689, 0.2735, 0.2781, 0.2826, 0.2872, 0.2917, 
      0.2961, 0.3003, 0.3045, 0.3086, 0.3125, 0.3162, 0.3198, 0.3235, 
      0.3270, 0.3306, 0.3342, 0.3377, 0.3411, 0.3446, 0.3479, 0.3511, 
      0.3543, 0.3575, 0.3605, 0.3634, 0.3663, 0.3691, 0.3720, 0.3748, 
      0.3775, 0.3803, 0.3830, 0.3857, 0.3884, 0.3911, 0.3937, 0.3962, 
      0.3987, 0.4012, 0.4036, 0.4060, 0.4084, 0.4107, 0.4129, 0.4152, 
      0.4174, 0.4196, 0.4218, 0.4239, 0.4261, 0.4282, 0.4303, 0.4324, 
      0.4344, 0.4365, 0.4385, 0.4405, 0.4425, 0.4445, 0.4464, 0.4483, 
      0.4502, 0.4521, 0.4539, 0.4558, 0.4576, 0.4593, 0.4611, 0.4629, 
      0.4646, 0.4663, 0.4680, 0.4697, 0.4714, 0.4730, 0.4747, 0.4759, 
      0.4769, 0.4780, 0.4790, 0.4800, 0.4811, 0.4827, 0.4843, 0.4859, 
      0.4874, 0.4889, 0.4905, 0.4920, 0.4935, 0.4950, 0.4965, 0.4979
    }
  }; 

  if(!fgLUT) fgLUT = new Double_t[AliTRDgeometry::kNlayer*fgkNlut];

  for (Int_t ilayer = 0; ilayer < AliTRDgeometry::kNlayer; ilayer++) {
    for (Int_t ilut  = 0; ilut  < fgkNlut; ilut++  ) {
      fgLUT[ilayer*fgkNlut+ilut] = lut[ilayer][ilut];
    }
  }
}

