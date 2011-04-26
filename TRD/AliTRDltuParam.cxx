#include <stdio.h>

#include "TMath.h"

#include "AliTRDltuParam.h"

// definition of geometry constants
Float_t AliTRDltuParam::fgZrow[6][5] = {
  {301, 177, 53, -57, -181},
  {301, 177, 53, -57, -181},
  {315, 184, 53, -57, -188},
  {329, 191, 53, -57, -195},
  {343, 198, 53, -57, -202},
  {347, 200, 53, -57, -204}};
Float_t AliTRDltuParam::fgX[6] =
  {300.65, 313.25, 325.85, 338.45, 351.05, 363.65};
Float_t AliTRDltuParam::fgTiltingAngle[6] =
  {-2., 2., -2., 2., -2., 2.};
Int_t   AliTRDltuParam::fgDyMax =  63;
Int_t   AliTRDltuParam::fgDyMin = -64;
Float_t AliTRDltuParam::fgBinDy = 140e-6;
Float_t AliTRDltuParam::fgWidthPad[6] =
  {0.635, 0.665, 0.695, 0.725, 0.755, 0.785};
Float_t AliTRDltuParam::fgLengthInnerPadC1[6] =
  {7.5, 7.5, 8.0, 8.5, 9.0, 9.0};
Float_t AliTRDltuParam::fgLengthOuterPadC1[6] =
  {7.5, 7.5, 7.5, 7.5, 7.5, 8.5};
Float_t AliTRDltuParam::fgLengthInnerPadC0 = 9.0;
Float_t AliTRDltuParam::fgLengthOuterPadC0 = 8.0;
Float_t AliTRDltuParam::fgScalePad = 256. * 32.;
Float_t AliTRDltuParam::fgDriftLength = 3.e-2;

AliTRDltuParam::AliTRDltuParam() :
  TObject(),
  fMagField(0.),
  fOmegaTau(0.),
  fPtMin(0.1),
  fNtimebins(20 << 5),
  fScaleQ0(0),
  fScaleQ1(0),
  fPidTracklengthCorr(kFALSE),
  fTiltCorr(kFALSE)
{

}

AliTRDltuParam::~AliTRDltuParam()
{

}

Int_t AliTRDltuParam::GetDyCorrection(Int_t det, Int_t rob, Int_t mcm) const
{
  // calculate the correction of the deflection
  // i.e. Lorentz angle and tilt correction (if active)

  Int_t layer = det % 6;

  Float_t dyTilt = ( fgDriftLength * TMath::Tan(fgTiltingAngle[layer]) *
  		     GetLocalZ(det, rob, mcm) / fgX[layer] );

  // calculate Lorentz correction
  Float_t dyCorr = - fOmegaTau * fgDriftLength;

  if(fTiltCorr)
    dyCorr += dyTilt; // add tilt correction

  return (int) TMath::Nint(dyCorr * fgScalePad / fgWidthPad[layer]);
}

void AliTRDltuParam::GetDyRange(Int_t det, Int_t rob, Int_t mcm, Int_t ch,
				 Int_t &dyMinInt, Int_t &dyMaxInt) const
{
  // calculate the deflection range in which tracklets are accepted

  dyMinInt = fgDyMin;
  dyMaxInt = fgDyMax;

  if (TMath::Abs(fMagField) < 0.1)
    return;

  Float_t e = 0.30;

  Float_t maxDeflTemp = GetPerp(det, rob, mcm, ch)/2. *         // Sekante/2
    (e * TMath::Abs(fMagField) / fPtMin);   // 1/R

  Float_t maxDeflAngle = 0.;

  if (maxDeflTemp < 1.) {
    maxDeflAngle = TMath::ASin(maxDeflTemp);

    Float_t dyMin = ( fgDriftLength *
		      tan(GetPhi(det, rob, mcm, ch) - maxDeflAngle) );

    dyMinInt = Int_t(dyMin / fgBinDy);
    if (dyMinInt < fgDyMin)
      dyMinInt = fgDyMin;

    Float_t dyMax = ( fgDriftLength *
  		      TMath::Tan(GetPhi(det, rob, mcm, ch) + maxDeflAngle) );

    dyMaxInt = Int_t(dyMax / fgBinDy);
    if (dyMaxInt > fgDyMax)
      dyMaxInt = fgDyMax;
  }
  if ((dyMaxInt - dyMinInt) <= 0) {
    printf("strange dy range: [%i,%i]\n", dyMinInt, dyMaxInt);
  }
}

Float_t AliTRDltuParam::GetElongation(Int_t det, Int_t rob, Int_t mcm, Int_t ch) const
{
  Int_t layer = det % 6;

  Float_t elongation = TMath::Abs(GetDist(det, rob, mcm, ch) / fgX[layer]);

  // sanity check
  if(elongation<0.001) {
    elongation=1.;
  }
  return elongation;
}

void AliTRDltuParam::GetCorrectionFactors(Int_t det, Int_t rob, Int_t mcm, Int_t ch,
					  UInt_t &cor0, UInt_t &cor1) const
{
  if (fPidTracklengthCorr == kTRUE ) {
    cor0 = Int_t ((1.0*fScaleQ0* (1./GetElongation(det, rob, mcm, ch)) ));
    cor1 = Int_t ((1.0*fScaleQ1* (1./GetElongation(det, rob, mcm, ch)) ));
  }
  else {
    cor0 = fScaleQ0;
    cor1 = fScaleQ1;
  }
}

Int_t AliTRDltuParam::GetNtimebins() const
{
  return fNtimebins;
}

Float_t AliTRDltuParam::GetX(Int_t det, Int_t /* rob */, Int_t /* mcm */) const
{
  Int_t layer = det%6;
  return fgX[layer];
}

Float_t AliTRDltuParam::GetLocalY(Int_t det, Int_t rob, Int_t mcm, Int_t ch) const
{
  Int_t layer = det%6;
  // calculate the pad position as in the TRAP
  Float_t ypos = (-4 + 1 + (rob&0x1) * 4 + (mcm&0x3)) * 18 - ch - 0.5; // y position in bins of pad widths
  return ypos*fgWidthPad[layer];
}

Float_t AliTRDltuParam::GetLocalZ(Int_t det, Int_t rob, Int_t mcm) const
{
  Int_t stack = (det%30) / 6;
  Int_t layer = det % 6;
  Int_t row   = (rob/2) * 4 + mcm/4;

  if (stack == 2) {
    if (row == 0)
      return (fgZrow[layer][stack] - 0.5 * fgLengthOuterPadC0);
    else if (row == 11)
      return (fgZrow[layer][stack] - 1.5 * fgLengthOuterPadC0 - (row - 1) * fgLengthInnerPadC0);
    else
      return (fgZrow[layer][stack] - fgLengthOuterPadC0 - (row - 0.5) * fgLengthInnerPadC0);
  }
  else {
    if (row == 0)
      return (fgZrow[layer][stack] - 0.5 * fgLengthOuterPadC1[layer]);
    else if (row == 15)
      return (fgZrow[layer][stack] - 1.5 * fgLengthOuterPadC1[layer] - (row - 1) * fgLengthInnerPadC1[layer]);
    else
      return (fgZrow[layer][stack] - fgLengthOuterPadC1[layer] - (row - 0.5) * fgLengthInnerPadC1[layer]);
  }
}

Float_t AliTRDltuParam::GetPerp(Int_t det, Int_t rob, Int_t mcm, Int_t ch) const
{
  return TMath::Sqrt(GetLocalY(det, rob, mcm, ch)*GetLocalY(det, rob, mcm, ch) +
		     GetX(det, rob, mcm)*GetX(det, rob, mcm) );
}

Float_t AliTRDltuParam::GetPhi(Int_t det, Int_t rob, Int_t mcm, Int_t ch) const
{
  return TMath::ATan2(GetLocalY(det, rob, mcm, ch), GetX(det, rob, mcm));
}

Float_t AliTRDltuParam::GetDist(Int_t det, Int_t rob, Int_t mcm, Int_t ch) const
{
  return TMath::Sqrt(GetLocalY(det, rob, mcm, ch)*GetLocalY(det, rob, mcm, ch) +
		     GetX(det, rob, mcm)*GetX(det, rob, mcm) +
		     GetLocalZ(det, rob, mcm)*GetLocalZ(det, rob, mcm) );
}
