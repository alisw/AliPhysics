#include "AliFMDESDFixer.h"
#include "AliESDFMD.h"
#include "AliFMDStripIndex.h"
#include "AliForwardUtil.h"
#include "AliForwardCorrectionManager.h"
#include "AliFMDCorrNoiseGain.h"
#include "AliLog.h"
#include <TH1.h>
#include <TList.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TMath.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TVector3.h>
#include <iostream>
#include <iomanip>


//____________________________________________________________________
AliFMDESDFixer::AliFMDESDFixer() 
  : TObject(), 
    fRecoFactor(1),
    fMaxNoiseCorr(0.05),
    fRecalculateEta(true),
    fXtraDead(0),
    fHasXtraDead(false),
    fInvalidIsEmpty(false),
    fNoiseChange(0),
    fEtaChange(0),
    fDeadChange(0)
{}

//____________________________________________________________________
AliFMDESDFixer::AliFMDESDFixer(const char*) 
  : TObject(), 
    fRecoFactor(1),
    fMaxNoiseCorr(0.05),
    fRecalculateEta(false),
    fXtraDead(AliFMDStripIndex::Pack(3,'O',19,511)+1),
    fHasXtraDead(false),
    fInvalidIsEmpty(false),
    fNoiseChange(0),
    fEtaChange(0),
    fDeadChange(0)
{}


//____________________________________________________________________
AliFMDESDFixer::AliFMDESDFixer(const AliFMDESDFixer& o) 
  : TObject(), 
    fRecoFactor(o.fRecoFactor),
    fMaxNoiseCorr(o.fMaxNoiseCorr),
    fRecalculateEta(o.fRecalculateEta),
    fXtraDead(o.fXtraDead),
    fHasXtraDead(o.fHasXtraDead),
    fInvalidIsEmpty(o.fInvalidIsEmpty),
    fNoiseChange(0),
    fEtaChange(0),
    fDeadChange(0)
{}

//____________________________________________________________________
AliFMDESDFixer&
AliFMDESDFixer::operator=(const AliFMDESDFixer& o) 
{ 
  if (&o == this) return *this;

  fRecoFactor = o.fRecoFactor;
  fMaxNoiseCorr = o.fMaxNoiseCorr;
  fRecalculateEta = o.fRecalculateEta;
  fXtraDead = o.fXtraDead;
  fHasXtraDead = o.fHasXtraDead;
  fInvalidIsEmpty = o.fInvalidIsEmpty;
  fNoiseChange = 0;
  fEtaChange  = 0;
  fDeadChange = 0;

  return *this; 
};

//____________________________________________________________________
void
AliFMDESDFixer::CreateOutputObjects(TList* l)
{
  TList* d = new TList;
  d->SetOwner();
  d->SetName(GetName());
  l->Add(d);

  d->Add(AliForwardUtil::MakeParameter("recoFactor", fRecoFactor));
  d->Add(AliForwardUtil::MakeParameter("recalcEta",  fRecalculateEta));
  d->Add(AliForwardUtil::MakeParameter("invalidIsEmpty",  fInvalidIsEmpty));
  
  fNoiseChange = new TH1D("noiseChange", "#delta#Delta due to noise",30,0,.3);
  fNoiseChange->SetDirectory(0);
  fNoiseChange->SetFillColor(kYellow+1);
  fNoiseChange->SetFillStyle(3001);
  d->Add(fNoiseChange);

  fEtaChange = new TH1D("etaChange", "#delta#eta",40,-1,1);
  fEtaChange->SetDirectory(0);
  fEtaChange->SetFillColor(kCyan+1);
  fEtaChange->SetFillStyle(3001);
  d->Add(fEtaChange);

  fDeadChange = new TH1D("deadChange", "#deltaN_{dead}",100,0,51200);
  fDeadChange->SetDirectory(0);
  fDeadChange->SetFillColor(kMagenta+1);
  fDeadChange->SetFillStyle(3001);
  d->Add(fDeadChange);

  TObjArray* extraDead = new TObjArray;
  extraDead->SetOwner();
  extraDead->SetName("extraDead");
  fXtraDead.Compact();
  UInt_t firstBit = fXtraDead.FirstSetBit();
  UInt_t nBits    = fXtraDead.GetNbits();
  for (UInt_t i = firstBit; i < nBits; i++) {
    if (!fXtraDead.TestBitNumber(i)) continue;
    UShort_t dd, s, t;
    Char_t   r;
    AliFMDStripIndex::Unpack(i, dd, r, s, t);
    extraDead->Add(AliForwardUtil::MakeParameter(Form("FMD%d%c[%02d,%03d]",
						      dd, r, s, t), 
						 UShort_t(i)));
  }
  d->Add(extraDead);
  fHasXtraDead = nBits > 0;

}

//____________________________________________________________________
void
AliFMDESDFixer::AddDead(UShort_t d, Char_t r, UShort_t s, UShort_t t)
{
  if (d < 1 || d > 3) {
    Warning("AddDead", "Invalid detector FMD%d", d);
    return;
  }
  Bool_t inner = (r == 'I' || r == 'i');
  if (d == 1 && !inner) { 
    Warning("AddDead", "Invalid ring FMD%d%c", d, r);
    return;
  }
  if ((inner && s >= 20) || (!inner && s >= 40)) { 
    Warning("AddDead", "Invalid sector FMD%d%c[%02d]", d, r, s);
    return;
  }
  if ((inner && t >= 512) || (!inner && t >= 256)) { 
    Warning("AddDead", "Invalid strip FMD%d%c[%02d,%03d]", d, r, s, t);
    return;
  }
    
  Int_t id = AliFMDStripIndex::Pack(d, r, s, t);
  // Int_t i  = 0;
  fXtraDead.SetBitNumber(id, true);
}
//____________________________________________________________________
void
AliFMDESDFixer::AddDeadRegion(UShort_t d,  Char_t r, 
				   UShort_t s1, UShort_t s2, 
				   UShort_t t1, UShort_t t2)
{
  // Add a dead region spanning from FMD<d><r>[<s1>,<t1>] to 
  // FMD<d><r>[<s2>,<t2>] (both inclusive)
  for (Int_t s = s1; s <= s2; s++) 
    for (Int_t t = t1; t <= t2; t++) 
      AddDead(d, r, s, t);
}
//____________________________________________________________________
void
AliFMDESDFixer::AddDead(const Char_t* script)
{
  if (!script || script[0] == '\0') return;
  
  const char* scr = gSystem->Which(gROOT->GetMacroPath(), script);
  if (!scr) {
    AliWarningF("%s not found in %s", script, gROOT->GetMacroPath());
    return;
  }
  AliInfoF("Reading additional dead strips from %s", scr);
  
  gROOT->Macro(Form("%s((AliFMDESDFixer*)%p);", scr, this));

  gInterpreter->UnloadFile(scr);
  delete scr;
}

//____________________________________________________________________
Bool_t
AliFMDESDFixer::IsDead(UShort_t d, Char_t r, UShort_t s, UShort_t t) const
{
  Int_t id = AliFMDStripIndex::Pack(d, r, s, t);
  return fXtraDead.TestBitNumber(id); 
}

//____________________________________________________________________
Int_t
AliFMDESDFixer::FindTargetNoiseFactor(const AliESDFMD& esd, Bool_t check) const
{
  if (!IsUseNoiseCorrection()) 
    // If the reconstruction factor was high (4 or more), do nothing 
    return 0;

  Int_t target = 0;
  if (AliESDFMD::Class_Version() < 4) {
    // IF we running with older STEER - we fix it here
    target = 4;					
  } else {
#if 1
    if (!esd.TestBit(1 << 14)) { 
      // If the bit isn't set, do nothing
      return 0;
    }
#else 
    // Uncommented until Peter commits patch to STEER/ESD
    if (!esd.NeedNoiseFix()) { 
      // If the bit isn't set, do nothing
      return 0;
    }
#endif
    target = Int_t(esd.GetNoiseFactor());
  }
  // Get the target factor - even thought the method below returns a
  // floating point value, we know that the noise factor is always
  // integer, so we coerce it to be the same here. 
  target -= fRecoFactor;

  // If the target factor is the same or smaller than the assumed
  // factor, we have nothing to do here, and we return immediately
  if (target <= 0) return 0;

  // Get the scaled noise from the correction mananger 
  if (check && !AliForwardCorrectionManager::Instance().GetNoiseGain()) 
    return 0;

  return target;
}

//____________________________________________________________________
#define ETA2COS(ETA)						\
  TMath::Cos(2*TMath::ATan(TMath::Exp(-TMath::Abs(ETA))))

//____________________________________________________________________
void
AliFMDESDFixer::Fix(AliESDFMD& esd, const TVector3& ip)
{

  const AliFMDCorrNoiseGain* ng  = 0;
  Int_t tgtFactor = FindTargetNoiseFactor(esd, false);
  if (tgtFactor > 0) 
    ng  = AliForwardCorrectionManager::Instance().GetNoiseGain();

  if (!ng && !fHasXtraDead && !fRecalculateEta && !fInvalidIsEmpty) 
    // We have nothing to do!
    return;

  UShort_t nDead = 0;
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nQ = d == 1 ? 1 : 2;

    for (UShort_t q = 0; q < nQ; q++) { 
      Char_t   r  = (q == 0 ? 'I' : 'O');
      UShort_t nS = (q == 0 ?  20 :  40);
      UShort_t nT = (q == 0 ? 512 : 256);
      
      for (UShort_t s = 0; s < nS; s++) { 
	for (UShort_t t = 0; t < nT; t++) { 
	  Double_t mult     = esd.Multiplicity(d,r,s,t);
	  Double_t eta      = esd.Eta(d,r,s,t);
	  Double_t cosTheta = 0;

	  if (CheckDead(d,r,s,t,mult)) nDead++;
	  
	  // Possibly re-calculate eta 
	  if (fRecalculateEta) RecalculateEta(d,r,s,t,ip,eta,mult,cosTheta);

	  // Possibly correct for poor treatment of ZS in reconstruction. 
	  if (ng && mult != AliESDFMD::kInvalidMult) {
	    if (cosTheta <= 0) cosTheta = ETA2COS(eta);
	    if (!NoiseCorrect(tgtFactor,ng->Get(d,r,s,t), cosTheta, mult))
	      nDead++;
	  } 

	  // Write out final values to object 
	  if (mult >= AliESDFMD::kInvalidMult) mult = AliESDFMD::kInvalidMult;
	  esd.SetMultiplicity(d,r,s,t,mult);
	  esd.SetEta(d,r,s,t,eta);
	} // for t
      } // for s
    } // for q
  } // for d
  fDeadChange->Fill(nDead);
}

//____________________________________________________________________
Bool_t
AliFMDESDFixer::CheckDead(UShort_t d, Char_t r, UShort_t s, UShort_t t,
			  Double_t& mult)
{
  // Correct for zero's being flagged as invalid 
  if (mult == AliESDFMD::kInvalidMult && fInvalidIsEmpty) mult = 0;

  // Take into account what we're defined as dead 
  if (IsDead(d,r,s,t)) {
    mult = AliESDFMD::kInvalidMult;
    return true;
  }
  return false;
}

//____________________________________________________________________
void
AliFMDESDFixer::RecalculateEta(UShort_t d, Char_t r, UShort_t s, UShort_t t,
			       const TVector3& ip, Double_t& eta,
			       Double_t& mult, 
			       Double_t& cosTheta)
{
  Double_t oldEta = eta, newEta = eta, newPhi=0;
  // Double_t newEta = AliForwardUtil::GetEtaFromStrip(d,r,s,t, zvtx);
  // eta             = newEta;
  cosTheta = ETA2COS(eta);
  if (!AliForwardUtil::GetEtaPhi(d, r, s, t, ip, newEta, newPhi)) return;
  if (TMath::Abs(eta) < 1) {
    ::Warning("RecalculateEta",
	      "FMD%d%c[%2d,%3d] (%f,%f,%f) eta=%f phi=%f (was %f)",
	      d, r, s, t, ip.X(), ip.Y(), ip.Z(), newEta, newPhi, oldEta);
    return;
  }
  
  eta = newEta;
  fEtaChange->Fill(newEta-oldEta);

  if (mult == AliESDFMD::kInvalidMult) return;

  Double_t newCos = ETA2COS(newEta);
  Double_t oldCos = ETA2COS(oldEta);
  Double_t corr   = newCos / oldCos;
  cosTheta        = newCos;
  mult            *= corr;
}
//____________________________________________________________________
Bool_t
AliFMDESDFixer::NoiseCorrect(Int_t target, Double_t corr, Double_t cosTheta, 
			     Double_t& mult)
{
  if (corr > fMaxNoiseCorr || corr <= 0) { 
    mult = AliESDFMD::kInvalidMult;
    return false;
  }
  Double_t add = corr * target * cosTheta;
  fNoiseChange->Fill(add);

  mult += add;
  return true;
}

#define PF(N,V,...)					\
  AliForwardUtil::PrintField(N,V, ## __VA_ARGS__)
#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)
#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)

//____________________________________________________________________
void
AliFMDESDFixer::Print(Option_t*) const
{
  AliForwardUtil::PrintTask(*this);
  gROOT->IncreaseDirLevel();
  PFB("Consider invalid null", fInvalidIsEmpty);
  PFB("Has extra dead", fHasXtraDead || fXtraDead.GetNbits() > 0);
  PFV("Reco noise factor", fRecoFactor);
  PFV("Max noise corr", fMaxNoiseCorr);
  PFB("Recalc. eta", fRecalculateEta);
  gROOT->DecreaseDirLevel();
}

//____________________________________________________________________
// 
// EOF
// 
