#include "AliDisplacedVertexSelectionAD.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDAD.h"
#include "AliESDADfriend.h"
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliADCalibData.h>
#include <TList.h>
#include <TH2.h>
#include <TH1.h>
#include <TSpline.h>
#include <TMath.h>

//--------------------------------------------------------------------
//
// See
//
//  https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AD_Offlline
// 

namespace {
  const Bool_t kCminusA = true;
}
//--------------------------------------------------------------------
void 
AliDisplacedVertexSelectionAD::SetupForData(TList*      l,
					    const char* name,
					    Bool_t      mc)
{
  // --- Get Calibrations from OCDB ----------------------------------
  AliCDBManager* cdb   = AliCDBManager::Instance();
  AliCDBEntry*   slew  = cdb->Get("AD/Calib/TimeSlewing");
  AliCDBEntry*   calib = cdb->Get("AD/Calib/Data");
  if (slew)  fSlewing = static_cast<TList*>(slew->GetObject());
  if (calib) fCalib   = static_cast<AliADCalibData*>(calib->GetObject());
  
  // --- Make diagnostics --------------------------------------------
  TList* out = new TList;
  out->SetName("displacedVertexAD");
  out->SetOwner();
  l->Add(out);

  Int_t    nT        = 1000;
  Double_t minT      = -100;
  Double_t maxT      = +100;
  Int_t    nZ        = Int_t(2*fMaxBunch*37.5)*2; // half-centimeter bins
  Double_t minZ      = -fMaxBunch*37.5;
  Double_t maxZ      = +fMaxBunch*37.5;

  fAC       = new TH2D("ac", "t_{A} vs t_{C}", nT, minT, maxT, nT, minT, maxT);
  fAC->SetDirectory(0);
  fAC->SetXTitle("t_{A}");
  fAC->SetYTitle("t_{C}");

  const char* diffText = (kCminusA ?
			  "#Deltat = t_{C}-t_{A}" :
			  "#Deltat = t_{A}-t_{C}");
  const char* sumText  = (kCminusA ?
			  "#sumt = t_{C}+t_{A}" :
			  "#sumt = t_{A}+t_{C}");
  
  fSumDelta =  new TH2D("sumDelta","#Deltat vs #sumt",
			nT, minT, maxT, 1.5*nT, minT, 2*maxT);
  fSumDelta->SetXTitle(diffText);
  fSumDelta->SetYTitle(sumText);
  fSumDelta->SetDirectory(0);

  fSumDeltaSatA =  static_cast<TH2*>(fSumDelta->Clone("sumDeltaSatA"));
  fSumDeltaSatA->SetTitle("#Deltat vs #sumt (A)");
  fSumDeltaSatA->SetDirectory(0);
  fSumDeltaSatA->SetMarkerColor(kRed-3);
  fSumDeltaSatA->SetMarkerStyle(20);
  fSumDeltaSatA->SetMarkerSize(.2);
  
  fSumDeltaSatC =  static_cast<TH2*>(fSumDelta->Clone("sumDeltaSatC"));
  fSumDeltaSatC->SetTitle("#Deltat vs #sumt (C)");
  fSumDeltaSatC->SetDirectory(0);
  fSumDeltaSatC->SetMarkerColor(kBlue-3);
  fSumDeltaSatC->SetMarkerStyle(21);
  fSumDeltaSatC->SetMarkerSize(.2);
  
  fIPzAll = new TH1D("ipZ", "AD IP_{#it{z}}",
		     nZ, minZ, maxZ);
  fIPzAll->SetDirectory(0);
  fIPzAll->SetXTitle("IP_{#it{z}} [cm]");
  fIPzAll->SetYTitle("Events");
  fIPzAll->SetFillStyle(3001);
  fIPzAll->SetFillColor(kGray+2);
  fIPzAll->SetLineColor(kGray+2);

  fIPzSatA = static_cast<TH1*>(fIPzAll->Clone("ipZsatA"));
  fIPzSatA->SetTitle("Main-Satellite IP_{#it{z}} (A)");
  fIPzSatA->SetFillColor(kRed-3);
  fIPzSatA->SetLineColor(kRed-3);
  fIPzSatA->SetDirectory(0);

  fIPzSatC = static_cast<TH1*>(fIPzAll->Clone("ipZsatC"));
  fIPzSatC->SetTitle("Satellite-Main IP_{#it{z}} (C)");
  fIPzSatC->SetFillColor(kBlue-3);
  fIPzSatC->SetLineColor(kBlue-3);
  fIPzSatC->SetDirectory(0);

  fIPzMain = static_cast<TH1*>(fIPzAll->Clone("ipZmain"));
  fIPzMain->SetTitle("Main-Main IP_{#it{z}}");
  fIPzMain->SetFillColor(kGreen-3);
  fIPzMain->SetLineColor(kGreen-3);
  fIPzMain->SetDirectory(0);
  
  fIPzDelta = new TH2D("deltaIPz", "#Deltat vs IP_{#it{z}}",
		       nT/5, minT/2.5, maxT/2.5,  nZ, minZ, maxZ);
  fIPzDelta->SetDirectory(0);
  fIPzDelta->SetXTitle(diffText);
  fIPzDelta->SetYTitle("IP_{#it{z}} [cm]");

  fIPzSum = new TH2D("sumIPz", "#sumt vs IP_{#it{z}}",
		     1.5*nT, minT, 2*maxT, nZ, minZ, maxZ);
  fIPzSum->SetDirectory(0);
  fIPzSum->SetXTitle(sumText);
  fIPzSum->SetYTitle("IP_{#it{z}} [cm]");

  fIPzBunch = new TH2D("bunchIpz", "De-bunch # vs IP_{#it{z}}",
		       2*fMaxBunch+1, -fMaxBunch-.5, +fMaxBunch+.5,
		       nZ, minZ, maxZ);
  fIPzBunch->SetDirectory(0);
  fIPzBunch->SetXTitle("De-bunch #");
  fIPzBunch->SetYTitle("IP_{#it{z}} [cm]");
  
  out->Add(fAC);
  out->Add(fSumDelta);
  out->Add(fSumDeltaSatA);
  out->Add(fSumDeltaSatC);
  out->Add(fIPzAll);
  out->Add(fIPzSatA);
  out->Add(fIPzSatC);
  out->Add(fIPzMain);
  out->Add(fIPzDelta);
  out->Add(fIPzSum);
  out->Add(fIPzBunch);
}

//--------------------------------------------------------------------
void 
AliDisplacedVertexSelectionAD::Print(Option_t* option) const
{
  Printf("Displaced vertex from AD signal");
}

//--------------------------------------------------------------------
Bool_t
AliDisplacedVertexSelectionAD::Process(const AliESDEvent* esd)
{
  fIPz                      = kInvalidVtxZ;
  fEventType                = kUnknown;
  AliESDAD*       adESD     = esd->GetADData();
  AliESDfriend*   esdFriend = esd->FindFriend();
  AliESDADfriend* adFriend  = esdFriend->GetADfriend();
  
  Double_t tA   = MeanTime(adESD, adFriend, true);
  Double_t tC   = MeanTime(adESD, adFriend, false);
  if (tA  == -kInvalidTime || tC == -kInvalidTime) return true;

  fAC->Fill(tA, tC);
  
  Double_t sumT   = tA + tC;
  Double_t deltaT = tC - tA;

  fSumDelta->Fill(deltaT, sumT);

  // C [m/s] * deltaT [ns] / 2 * 10^-9 [s/ns] * 10^2 [cm/m]
  fIPz = TMath::C() * deltaT / 2 * 1e-9 * 1e2;
  fIPzAll->Fill(fIPz);

  fIPzDelta->Fill(deltaT, fIPz);
  fIPzSum->Fill(sumT, fIPz);
  
  Double_t r2 = 2*fSpacing/2*fSpacing/2;
  // Double_t r  = TMath::Sqrt(r2);
  Int_t    b  = -100;
  for (Int_t i = -fMaxBunch; i <= fMaxBunch; i++) {
    Double_t ir = i*fSpacing;
    Double_t x2 = TMath::Power(deltaT-ir,2);
    Double_t y2 = TMath::Power(sumT  -ir,2);
    Double_t u2 = TMath::Power(-sumT -ir,2);

    if ((x2+y2) < r2 || (x2+u2) < r2) {
      // Got a collision - which kind?
      b = i;
      if (i == 0)
	// Main-main 
	fEventType = kMain;      
      // else if (TMath::SignBit(deltaT) == TMath::SignBit(sumT))
      else if (std::signbit(deltaT) == std::signbit(sumT))
	// If difference and sum have same signs
	fEventType = (kCminusA ? kSatelliteA : kSatelliteC);
      else
	// If difference and sum have opposite signs
	fEventType = (kCminusA ? kSatelliteC : kSatelliteA);
      break;
    }
  }
  switch (fEventType) {
  case kMain:       fIPzMain->Fill(fIPz); break;
  case kSatelliteA: fIPzSatA->Fill(fIPz);fSumDeltaSatA->Fill(deltaT,sumT);break;
  case kSatelliteC: fIPzSatC->Fill(fIPz);fSumDeltaSatC->Fill(deltaT,sumT);break;
  default:                                        break;
  }

  if (b > -100) fIPzBunch->Fill(b, fIPz);
  
  return true;
}
//--------------------------------------------------------------------
Double_t
AliDisplacedVertexSelectionAD::MeanTime(AliESDAD*       adESD,
					AliESDADfriend* adFriend,
					Bool_t          aNotC) const
{
  if (!fCalib || !fSlewing) return -kInvalidTime;
  if (!adESD || !adFriend)  return -kInvalidTime;

  const Int_t chMin = (aNotC ?  8 : 0);
  const Int_t chMax = (aNotC ? 16 : 8);

  const Double_t tRes = fCalib->GetTimeResolution(aNotC);

  Int_t    cnt   = 0;
  Double_t sumW  = 0;
  Double_t sumWT = 0;
  for (Int_t ch = chMin; ch < chMax; ++ch) {
    const Double_t charge = adESD->GetAdc(ch);
    if (adFriend->GetWidth(ch) <= 1e-6 || charge <= 10.)  continue;

    TSpline3* sts = static_cast<TSpline3*>(fSlewing->At(ch));
    Double_t  q2  = TMath::Power(charge, 2);
    cnt++;
    sumW        += q2;
    sumWT       += q2 * (adFriend->GetTime(ch) -
			 tRes * sts->Eval(TMath::Log10(1./charge)));
  }

  return (cnt > 0 ? sumWT/sumW : -kInvalidTime);    
}
//--------------------------------------------------------------------
//
// EOF
// 
