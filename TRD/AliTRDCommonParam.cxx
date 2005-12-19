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
// Class containing constant common parameters                           //
//                                                                           //
// Request an instance with AliTRDCommonParam::Instance()                 //
// Then request the needed values                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliRun.h"

#include "AliTRDCommonParam.h"
#include "AliTRDpadPlane.h"


ClassImp(AliTRDCommonParam)

AliTRDCommonParam* AliTRDCommonParam::fgInstance = 0;
Bool_t AliTRDCommonParam::fgTerminated = kFALSE;

//_ singleton implementation __________________________________________________
AliTRDCommonParam* AliTRDCommonParam::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  // 
  
  if (fgTerminated != kFALSE)
    return 0;

  if (fgInstance == 0)
    fgInstance = new AliTRDCommonParam();
  
  return fgInstance;
}

void AliTRDCommonParam::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class and sets the terminated flag, instances cannot be requested anymore
  // This function can be called several times.
  //
  
  fgTerminated = kTRUE;
  
  if (fgInstance != 0)
  {
    delete fgInstance;
    fgInstance = 0;
  }
}

//_____________________________________________________________________________
AliTRDCommonParam::AliTRDCommonParam()
{
  //
  // constructor
  //
  
  fField              = 0.0;
  fPRFbin             = 0;
  fPRFlo              = 0.0;
  fPRFhi              = 0.0;
  fPRFwid             = 0.0;
  fPRFpad             = 0;

  fPRFsmp             = 0;
  
  fExBOn              = kFALSE;
  fPRFOn              = kFALSE;
  
  fPadPlaneArray      = 0;
  
  Init();
}

//_____________________________________________________________________________
void AliTRDCommonParam::Init()
{
  //
  // constructor helper
  //
  
  // E x B effects
  fExBOn          = kTRUE;

  // The pad response function
  fPRFOn          = kTRUE;

  // The magnetic field strength in Tesla
  Double_t x[3] = { 0.0, 0.0, 0.0 };
  Double_t b[3];
  gAlice->Field(x,b);  // b[] is in kilo Gauss
  fField = b[2] * 0.1; // Tesla
  
  // Create the sampled PRF
  SamplePRF();
  
  // ----------------------------------------------------------------------------
  // The pad planes
  // ----------------------------------------------------------------------------
  
  fPadPlaneArray = new TObjArray(kNplan * kNcham);
  
  for (Int_t iplan = 0; iplan < kNplan; iplan++) {
    for (Int_t icham = 0; icham < kNcham; icham++) {
      Int_t ipp = iplan + icham * kNplan;
      fPadPlaneArray->AddAt(new AliTRDpadPlane(iplan,icham),ipp);
    }
  }
}

//_____________________________________________________________________________
AliTRDCommonParam::~AliTRDCommonParam() 
{
  //
  // destructor
  //
  
  if (fPRFsmp) {
    delete [] fPRFsmp;
    fPRFsmp = 0;
  }

  if (fPadPlaneArray) {
    fPadPlaneArray->Delete();
    delete fPadPlaneArray;
    fPadPlaneArray = 0;
  }
}

//_____________________________________________________________________________
AliTRDCommonParam::AliTRDCommonParam(const AliTRDCommonParam &p):TObject(p)
{
  //
  // copy constructor
  //

  ((AliTRDCommonParam &) p).Copy(*this);
}


//_____________________________________________________________________________
AliTRDCommonParam &AliTRDCommonParam::operator=(const AliTRDCommonParam &p)
{
  //
  // Assignment operator
  //

  if (this != &p) ((AliTRDCommonParam &) p).Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void AliTRDCommonParam::Copy(TObject &p) const
{
  //
  // Copy function
  //
  
  AliTRDCommonParam* target = dynamic_cast<AliTRDCommonParam*> (&p);
  if (!target)
    return;
  
  target->fExBOn              = fExBOn;
  target->fPRFOn              = fPRFOn;
  target->fPRFbin             = fPRFbin;
  target->fPRFlo              = fPRFlo;
  target->fPRFhi              = fPRFhi;
  target->fPRFwid             = fPRFwid;
  target->fPRFpad             = fPRFpad;
  if (target->fPRFsmp) delete [] target->fPRFsmp;
  target->fPRFsmp = new Float_t[fPRFbin];
  for (Int_t iBin = 0; iBin < fPRFbin; iBin++) {
    target->fPRFsmp[iBin] = fPRFsmp[iBin];
  }
}

//_____________________________________________________________________________
void AliTRDCommonParam::SamplePRF()
{
  //
  // Samples the pad response function
  //

  const Int_t kPRFbin = 61;

  Float_t prf[kNplan][kPRFbin] = { {2.9037e-02, 3.3608e-02, 3.9020e-02, 4.5292e-02,
                    5.2694e-02, 6.1362e-02, 7.1461e-02, 8.3362e-02,
                    9.7063e-02, 1.1307e-01, 1.3140e-01, 1.5235e-01,
                    1.7623e-01, 2.0290e-01, 2.3294e-01, 2.6586e-01,
                    3.0177e-01, 3.4028e-01, 3.8077e-01, 4.2267e-01,
                    4.6493e-01, 5.0657e-01, 5.4655e-01, 5.8397e-01,
                    6.1767e-01, 6.4744e-01, 6.7212e-01, 6.9188e-01,
                    7.0627e-01, 7.1499e-01, 7.1851e-01, 7.1499e-01,
                    7.0627e-01, 6.9188e-01, 6.7212e-01, 6.4744e-01,
                    6.1767e-01, 5.8397e-01, 5.4655e-01, 5.0657e-01,
                    4.6493e-01, 4.2267e-01, 3.8077e-01, 3.4028e-01,
                    3.0177e-01, 2.6586e-01, 2.3294e-01, 2.0290e-01,
                    1.7623e-01, 1.5235e-01, 1.3140e-01, 1.1307e-01,
                    9.7063e-02, 8.3362e-02, 7.1461e-02, 6.1362e-02,
                    5.2694e-02, 4.5292e-02, 3.9020e-02, 3.3608e-02,
                    2.9037e-02},
                   {2.5478e-02, 2.9695e-02, 3.4655e-02, 4.0454e-02,
                    4.7342e-02, 5.5487e-02, 6.5038e-02, 7.6378e-02,
                    8.9696e-02, 1.0516e-01, 1.2327e-01, 1.4415e-01,
                    1.6794e-01, 1.9516e-01, 2.2573e-01, 2.5959e-01,
                    2.9694e-01, 3.3719e-01, 3.7978e-01, 4.2407e-01,
                    4.6889e-01, 5.1322e-01, 5.5569e-01, 5.9535e-01,
                    6.3141e-01, 6.6259e-01, 6.8882e-01, 7.0983e-01,
                    7.2471e-01, 7.3398e-01, 7.3761e-01, 7.3398e-01,
                    7.2471e-01, 7.0983e-01, 6.8882e-01, 6.6259e-01,
                    6.3141e-01, 5.9535e-01, 5.5569e-01, 5.1322e-01,
                    4.6889e-01, 4.2407e-01, 3.7978e-01, 3.3719e-01,
                    2.9694e-01, 2.5959e-01, 2.2573e-01, 1.9516e-01,
                    1.6794e-01, 1.4415e-01, 1.2327e-01, 1.0516e-01,
                    8.9696e-02, 7.6378e-02, 6.5038e-02, 5.5487e-02,
                    4.7342e-02, 4.0454e-02, 3.4655e-02, 2.9695e-02,
                    2.5478e-02},
                   {2.2363e-02, 2.6233e-02, 3.0782e-02, 3.6140e-02,
                    4.2535e-02, 5.0157e-02, 5.9197e-02, 6.9900e-02,
                    8.2707e-02, 9.7811e-02, 1.1548e-01, 1.3601e-01,
                    1.5998e-01, 1.8739e-01, 2.1840e-01, 2.5318e-01,
                    2.9182e-01, 3.3373e-01, 3.7837e-01, 4.2498e-01,
                    4.7235e-01, 5.1918e-01, 5.6426e-01, 6.0621e-01,
                    6.4399e-01, 6.7700e-01, 7.0472e-01, 7.2637e-01,
                    7.4206e-01, 7.5179e-01, 7.5551e-01, 7.5179e-01,
                    7.4206e-01, 7.2637e-01, 7.0472e-01, 6.7700e-01,
                    6.4399e-01, 6.0621e-01, 5.6426e-01, 5.1918e-01,
                    4.7235e-01, 4.2498e-01, 3.7837e-01, 3.3373e-01,
                    2.9182e-01, 2.5318e-01, 2.1840e-01, 1.8739e-01,
                    1.5998e-01, 1.3601e-01, 1.1548e-01, 9.7811e-02,
                    8.2707e-02, 6.9900e-02, 5.9197e-02, 5.0157e-02,
                    4.2535e-02, 3.6140e-02, 3.0782e-02, 2.6233e-02,
                    2.2363e-02},
                   {1.9635e-02, 2.3167e-02, 2.7343e-02, 3.2293e-02,
                    3.8224e-02, 4.5335e-02, 5.3849e-02, 6.4039e-02,
                    7.6210e-02, 9.0739e-02, 1.0805e-01, 1.2841e-01,
                    1.5216e-01, 1.7960e-01, 2.1099e-01, 2.4671e-01,
                    2.8647e-01, 3.2996e-01, 3.7660e-01, 4.2547e-01,
                    4.7536e-01, 5.2473e-01, 5.7215e-01, 6.1632e-01,
                    6.5616e-01, 6.9075e-01, 7.1939e-01, 7.4199e-01,
                    7.5838e-01, 7.6848e-01, 7.7227e-01, 7.6848e-01,
                    7.5838e-01, 7.4199e-01, 7.1939e-01, 6.9075e-01,
                    6.5616e-01, 6.1632e-01, 5.7215e-01, 5.2473e-01,
                    4.7536e-01, 4.2547e-01, 3.7660e-01, 3.2996e-01,
                    2.8647e-01, 2.4671e-01, 2.1099e-01, 1.7960e-01,
                    1.5216e-01, 1.2841e-01, 1.0805e-01, 9.0739e-02,
                    7.6210e-02, 6.4039e-02, 5.3849e-02, 4.5335e-02,
                    3.8224e-02, 3.2293e-02, 2.7343e-02, 2.3167e-02,
                    1.9635e-02},
                   {1.7224e-02, 2.0450e-02, 2.4286e-02, 2.8860e-02,
                    3.4357e-02, 4.0979e-02, 4.8966e-02, 5.8612e-02,
                    7.0253e-02, 8.4257e-02, 1.0102e-01, 1.2094e-01,
                    1.4442e-01, 1.7196e-01, 2.0381e-01, 2.4013e-01,
                    2.8093e-01, 3.2594e-01, 3.7450e-01, 4.2563e-01,
                    4.7796e-01, 5.2991e-01, 5.7974e-01, 6.2599e-01,
                    6.6750e-01, 7.0344e-01, 7.3329e-01, 7.5676e-01,
                    7.7371e-01, 7.8410e-01, 7.8793e-01, 7.8410e-01,
                    7.7371e-01, 7.5676e-01, 7.3329e-01, 7.0344e-01,
                    6.6750e-01, 6.2599e-01, 5.7974e-01, 5.2991e-01,
                    4.7796e-01, 4.2563e-01, 3.7450e-01, 3.2594e-01,
                    2.8093e-01, 2.4013e-01, 2.0381e-01, 1.7196e-01,
                    1.4442e-01, 1.2094e-01, 1.0102e-01, 8.4257e-02,
                    7.0253e-02, 5.8612e-02, 4.8966e-02, 4.0979e-02,
                    3.4357e-02, 2.8860e-02, 2.4286e-02, 2.0450e-02,
                    1.7224e-02},
                   {1.5096e-02, 1.8041e-02, 2.1566e-02, 2.5793e-02,
                    3.0886e-02, 3.7044e-02, 4.4515e-02, 5.3604e-02,
                    6.4668e-02, 7.8109e-02, 9.4364e-02, 1.1389e-01,
                    1.3716e-01, 1.6461e-01, 1.9663e-01, 2.3350e-01,
                    2.7527e-01, 3.2170e-01, 3.7214e-01, 4.2549e-01,
                    4.8024e-01, 5.3460e-01, 5.8677e-01, 6.3512e-01,
                    6.7838e-01, 7.1569e-01, 7.4655e-01, 7.7071e-01,
                    7.8810e-01, 7.9871e-01, 8.0255e-01, 7.9871e-01,
                    7.8810e-01, 7.7071e-01, 7.4655e-01, 7.1569e-01,
                    6.7838e-01, 6.3512e-01, 5.8677e-01, 5.3460e-01,
                    4.8024e-01, 4.2549e-01, 3.7214e-01, 3.2170e-01,
                    2.7527e-01, 2.3350e-01, 1.9663e-01, 1.6461e-01,
                    1.3716e-01, 1.1389e-01, 9.4364e-02, 7.8109e-02,
                    6.4668e-02, 5.3604e-02, 4.4515e-02, 3.7044e-02,
                    3.0886e-02, 2.5793e-02, 2.1566e-02, 1.8041e-02,
                    1.5096e-02}};

  // More sampling precision with linear interpolation
  fPRFlo  = -1.5;
  fPRFhi  =  1.5;
  Float_t pad[kPRFbin];
  Int_t   sPRFbin = kPRFbin;  
  Float_t sPRFwid = (fPRFhi - fPRFlo) / ((Float_t) sPRFbin);
  for (Int_t iPad = 0; iPad < sPRFbin; iPad++) {
    pad[iPad] = ((Float_t) iPad + 0.5) * sPRFwid + fPRFlo;
  }
  fPRFbin = 500;  
  fPRFwid = (fPRFhi - fPRFlo) / ((Float_t) fPRFbin);
  fPRFpad = ((Int_t) (1.0 / fPRFwid));

  if (fPRFsmp) delete [] fPRFsmp;
  fPRFsmp = new Float_t[kNplan*fPRFbin];

  Int_t   ipos1;
  Int_t   ipos2;
  Float_t diff;

  for (Int_t iPla = 0; iPla < kNplan; iPla++) {

    for (Int_t iBin = 0; iBin < fPRFbin; iBin++) {

      Float_t bin = (((Float_t) iBin) + 0.5) * fPRFwid + fPRFlo;
      ipos1 = ipos2 = 0;
      diff  = 0;
      do {
        diff = bin - pad[ipos2++];
      } while ((diff > 0) && (ipos2 < kPRFbin));
      if      (ipos2 == kPRFbin) {
        fPRFsmp[iPla*fPRFbin+iBin] = prf[iPla][ipos2-1];
      }
      else if (ipos2 == 1) {
        fPRFsmp[iPla*fPRFbin+iBin] = prf[iPla][ipos2-1];
      }
      else {
        ipos2--;
        if (ipos2 >= kPRFbin) ipos2 = kPRFbin - 1;
        ipos1 = ipos2 - 1;
        fPRFsmp[iPla*fPRFbin+iBin] = prf[iPla][ipos2] 
                                   + diff * (prf[iPla][ipos2] - prf[iPla][ipos1]) 
                                          / sPRFwid;
      }

    }
  } 

}

//_____________________________________________________________________________
Int_t AliTRDCommonParam::PadResponse(Double_t signal, Double_t dist
    , Int_t plane, Double_t *pad) const
{
  //
  // Applies the pad response
  //

  Int_t iBin  = ((Int_t) (( - dist - fPRFlo) / fPRFwid));
  Int_t iOff  = plane * fPRFbin;

  Int_t iBin0 = iBin - fPRFpad + iOff;
  Int_t iBin1 = iBin           + iOff;
  Int_t iBin2 = iBin + fPRFpad + iOff;

  pad[0] = 0.0;
  pad[1] = 0.0;
  pad[2] = 0.0;
  if ((iBin1 >= 0) && (iBin1 < (fPRFbin*kNplan))) {

    if (iBin0 >= 0) {
      pad[0] = signal * fPRFsmp[iBin0];
    }
    pad[1] = signal * fPRFsmp[iBin1];
    if (iBin2 < (fPRFbin*kNplan)) {
      pad[2] = signal * fPRFsmp[iBin2];
    }

    return 1;

  }
  else {

    return 0;

  }

}

//_____________________________________________________________________________
AliTRDpadPlane *AliTRDCommonParam::GetPadPlane(Int_t p, Int_t c) const
{
  //
  // Returns the pad plane for a given plane <p> and chamber <c> number
  //

  Int_t ipp = p + c * kNplan;
  return ((AliTRDpadPlane *) fPadPlaneArray->At(ipp));

}

//_____________________________________________________________________________
Int_t AliTRDCommonParam::GetRowMax(Int_t p, Int_t c, Int_t /*s*/) const
{
  //
  // Returns the number of rows on the pad plane
  //

  return GetPadPlane(p,c)->GetNrows();

}

//_____________________________________________________________________________
Int_t AliTRDCommonParam::GetColMax(Int_t p) const
{
  //
  // Returns the number of rows on the pad plane
  //

  return GetPadPlane(p,0)->GetNcols();

}

//_____________________________________________________________________________
Double_t AliTRDCommonParam::GetRow0(Int_t p, Int_t c, Int_t /*s*/) const
{
  //
  // Returns the position of the border of the first pad in a row
  //

  return GetPadPlane(p,c)->GetRow0();

}

//_____________________________________________________________________________
Double_t AliTRDCommonParam::GetCol0(Int_t p) const
{
  //
  // Returns the position of the border of the first pad in a column
  //

  return GetPadPlane(p,0)->GetCol0();

}
