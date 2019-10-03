// 
// Utilities used in the forward multiplcity analysis 
// 
//
#include "AliForwardUtil.h"
//#include <ARVersion.h>
#include <AliAnalysisManager.h>
#include "AliAODForwardMult.h"
#include <AliLog.h>
#include <AliInputEventHandler.h>
#include <AliAODInputHandler.h>
#include <AliAODHandler.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliAnalysisTaskSE.h>
#include <AliPhysicsSelection.h>
#include <AliTriggerAnalysis.h>
#include <AliMultiplicity.h>
#include "AliMultEstimator.h" // <-- Sigh, needed by AliMultSelection
#include "AliMultVariable.h"  // <-- Sigh, needed by AliMultSelection
#include "AliMultInput.h"     // <-- Sigh, needed by AliMultSelection
#include "AliMultSelection.h"
#include "AliMultSelectionCuts.h"
#include <TParameter.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TMath.h>
#include <TError.h>
#include <TROOT.h>
#include <TVector3.h>
#define FIT_OPTIONS "RNS"

//====================================================================
ULong_t AliForwardUtil::AliROOTRevision()
{
#ifdef ALIROOT_SVN_REVISION
  return ALIROOT_SVN_REVISION;
#elif defined(ALIROOT_REVISION)
  static ULong_t ret = 0;
  if (ret != 0) return ret;

  // Select first 32bits of the 40byte long check-sum
  TString rev(ALIROOT_REVISION, 8);
  for (ULong_t i = 0; i < 8; i++) {
    ULong_t p = 0;
    switch (rev[i]) { 
    case '0': p = 0; break;
    case '1': p = 1; break;
    case '2': p = 2; break;
    case '3': p = 3; break;
    case '4': p = 4; break;
    case '5': p = 5; break;
    case '6': p = 6; break;
    case '7': p = 7; break;
    case '8': p = 8; break;
    case '9': p = 9; break;
    case 'a': case 'A': p = 10; break;
    case 'b': case 'B': p = 11; break;
    case 'c': case 'C': p = 12; break;
    case 'd': case 'D': p = 13; break;
    case 'e': case 'E': p = 14; break;
    case 'f': case 'F': p = 15; break;
    }
    ret |= (p << (32-4*(i+1)));
  }
  return ret;
#else
  return 0;
#endif
}
//____________________________________________________________________
ULong_t AliForwardUtil::AliROOTBranch()
{
  // Do something here when we switch to git - sigh!
#if !defined(ALIROOT_SVN_BRANCH) && !defined(ALIROOT_BRANCH) 
  return 0;
#endif
  static ULong_t ret = 0;
  if (ret != 0) return ret;
  TString str;
  TString top;
#ifdef ALIROOT_SVN_BRANCH
  str = ALIROOT_SVN_BRANCH;
  top = "trunk";
#elif defined(ALIROOT_BRANCH)
  str = ALIROOT_BRANCH;
  top = "master";
#endif
  if (str.IsNull()) return 0xFFFFFFFF;
  if (str[0] == 'v') str.Remove(0,1);
  if (str.EqualTo(top)) return ret = 0xFFFFFFFF;

  TObjArray*   tokens = str.Tokenize("-");
  TObjString*  pMajor = tokens->GetEntries()>0 ? 
    (static_cast<TObjString*>(tokens->At(0))) : 0;
  TObjString*  pMinor = tokens->GetEntries()>1 ? 
    (static_cast<TObjString*>(tokens->At(1))) : 0;
  TObjString*  pRelea = tokens->GetEntries() > 2 ? 
    static_cast<TObjString*>(tokens->At(2)) : 0;
  TObjString* pAn     = tokens->GetEntries() > 3 ? 
    static_cast<TObjString*>(tokens->At(3)) : 0;
  TString sMajor,sMinor,sRelea;
  if (pMajor) sMajor = pMajor->String().Strip(TString::kLeading, '0'); 
  if (pMinor) sMinor = pMinor->String().Strip(TString::kLeading, '0');
  if (pRelea) sRelea = pRelea->String().Strip(TString::kLeading, '0');
  //
  ret = (((sMajor.Atoi() & 0xFF) << 12) |
    ((sMinor.Atoi() & 0xFF) <<  8) |
    ((sRelea.Atoi() & 0xFF) <<  4) |
    (pAn ? 0xAA : 0));
  
  return ret;
}

//====================================================================
UShort_t
AliForwardUtil::ParseCollisionSystem(const char* sys)
{
  // 
  // Parse a collision system spec given in a string.   Known values are 
  // 
  //  - "ppb", "p-pb", "pa", "p-a"  which returns kPPb
  //  - "pp", "p-p"                 which returns kPP 
  //  - "PbPb", "Pb-Pb", "A-A",     which returns kPbPb 
  //  - Everything else gives kUnknown 
  // 
  // Parameters:
  //    sys Collision system spec 
  // 
  // Return:
  //    Collision system id 
  //
  TString s(sys);
  s.ToLower();
  // we do pA first to avoid pp catch on ppb string (AH)
  if (s.Contains("p-pb")  || s.Contains("ppb"))   return AliForwardUtil::kPPb;
  if (s.Contains("p-a")   || s.Contains("pa"))    return AliForwardUtil::kPPb;
  if (s.Contains("pb-p")  || s.Contains("pbp"))   return AliForwardUtil::kPbp;
  if (s.Contains("a-p")   || s.Contains("ap"))    return AliForwardUtil::kPbp;
  if (s.Contains("p-p")   || s.Contains("pp"))    return AliForwardUtil::kPP; 
  if (s.Contains("pb-pb") || s.Contains("pbpb"))  return AliForwardUtil::kPbPb;
  if (s.Contains("xe-xe") || s.Contains("xexe"))  return AliForwardUtil::kXeXe;
  if (s.Contains("xe54-xe54"))                    return AliForwardUtil::kXeXe;
  if (s.Contains("a-a")   || s.Contains("aa"))    return AliForwardUtil::kPbPb;
  return AliForwardUtil::kUnknown;
}
//____________________________________________________________________
UShort_t
AliForwardUtil::ParseCollisionSystem(Int_t b1a, Int_t b1z,
				     Int_t b2a, Int_t b2z,
				     const char* sys)
{
  if (b1a == 129 && b1z == 54 && b2a == 129 && b2z == 54)
    return AliForwardUtil::kXeXe;
  if (b1a == 208 && b1z == 82 && b2a == 208 && b2z == 82)
    return AliForwardUtil::kPbPb;
  if (b1a == 208 && b1z == 82 && b2a ==   1 && b2z ==  1)
    return AliForwardUtil::kPbp;
  if (b1a ==   1 && b1z ==  1 && b2a == 208 && b2z == 82)
    return AliForwardUtil::kPPb;
  if (b1a ==   1 && b1z ==  1 && b2a ==   1 && b2z ==  1)
    return AliForwardUtil::kPP;
  return ParseCollisionSystem(sys);
}

//____________________________________________________________________
const char*
AliForwardUtil::CollisionSystemString(UShort_t sys)
{
  // 
  // Get a string representation of the collision system 
  // 
  // Parameters:
  //    sys  Collision system 
  // - kPP   -> "pp"
  // - kPbPb -> "PbPb"
  // - kPPb  -> "pPb"
  // - kPbp  -> "Pbp"
  // - anything else gives "unknown"
  // 
  // Return:
  //    String representation of the collision system 
  //
  switch (sys) { 
  case AliForwardUtil::kPP:   return "pp";
  case AliForwardUtil::kPbPb: return "PbPb";
  case AliForwardUtil::kPPb:  return "pPb";
  case AliForwardUtil::kPbp:  return "Pbp";
  case AliForwardUtil::kXeXe: return "XeXe";
  }
  return "unknown";
}
//____________________________________________________________________
Float_t
AliForwardUtil::BeamRapidity(Float_t beam, UShort_t z, UShort_t a)
{
  const Double_t pMass = 9.38271999999999995e-01;
  const Double_t nMass = 9.39564999999999984e-01;
  Double_t       beamE = z * beam / 2;
  Double_t       beamM = z * pMass + (a - z) * nMass;
  Double_t       beamP = TMath::Sqrt(beamE * beamE - beamM * beamM);
  Double_t       beamY = .5* TMath::Log((beamE+beamP) / (beamE-beamP));
  return beamY;
}
//____________________________________________________________________
Float_t
AliForwardUtil::CenterOfMassEnergy(Float_t beam, 
				   UShort_t z1, 
				   UShort_t a1, 
				   Short_t z2, 
				   Short_t a2) 
{
  // Calculate the center of mass energy given target/projectile 
  // mass and charge numbers
  if (z2 < 0) z2 = z1;
  if (a2 < 0) a2 = a1;
  return TMath::Sqrt(Float_t(z1*z2)/a1/a2) * beam;
}
//____________________________________________________________________
Float_t
AliForwardUtil::CenterOfMassRapidity(UShort_t z1, 
				     UShort_t a1, 
				     Short_t z2, 
				     Short_t a2) 
{
  // Calculate the center of mass rapidity (shift) given target/projectile 
  // mass and charge numbers
  if (z2 < 0) z2 = z1;
  if (a2 < 0) a2 = a1;
  if (z2 == z1 && a2 == a1) return 0;
  return .5 * TMath::Log(Float_t(z1*a2)/z2/a1);
}

namespace {
  UShort_t CheckSNN(Float_t energy)
  {
    if (TMath::Abs(energy - 900.)   < 10)  return 900;
    if (TMath::Abs(energy - 2400.)  < 10)  return 2400;
    if (TMath::Abs(energy - 2760.)  < 20)  return 2760;
    if (TMath::Abs(energy - 4400.)  < 10)  return 4400;
    if (TMath::Abs(energy - 5000.)  < 10)  return 5000;
    if (TMath::Abs(energy - 5022.)  < 10)  return 5023;
    if (TMath::Abs(energy - 5125.)  < 30)  return 5100;
    if (TMath::Abs(energy - 5200.)  < 50)  return 5200;
    if (TMath::Abs(energy - 5440.)  < 20)  return 5440;
    if (TMath::Abs(energy - 5500.)  < 40)  return 5500;
    if (TMath::Abs(energy - 7000.)  < 10)  return 7000;
    if (TMath::Abs(energy - 8000.)  < 10)  return 8000;
    if (TMath::Abs(energy - 10000.) < 10)  return 10000;
    if (TMath::Abs(energy - 13000.) < 10)  return 13000;
    if (TMath::Abs(energy - 14000.) < 10)  return 14000;
    return 0;
  }
}
//____________________________________________________________________
UShort_t
AliForwardUtil::ParseCenterOfMassEnergy(UShort_t sys, Float_t beam)
{
  // 
  // Parse the center of mass energy given as a float and return known 
  // values as a unsigned integer
  // 
  // Parameters:
  //    sys   Collision system (needed for AA)
  //    beam  Center of mass energy * total charge 
  // 
  // Return:
  //    Center of mass energy per nucleon
  //
  Float_t energy = beam; 
  // Below no longer needed apparently
  // if (sys == AliForwardUtil::kPbPb) energy = energy / 208 * 82;
  if (sys == AliForwardUtil::kPPb) 
    energy = CenterOfMassEnergy(beam, 82, 208, 1, 1);
  if (sys == AliForwardUtil::kPbp) 
    energy = CenterOfMassEnergy(beam, 1, 1, 82, 208);
  else if (sys == AliForwardUtil::kPbPb) 
    energy = CenterOfMassEnergy(beam, 82, 208, 82, 208);
  else if (sys == AliForwardUtil::kXeXe) 
    energy = CenterOfMassEnergy(beam, 54, 129, 54, 129);
  UShort_t ret = CheckSNN(energy);
  if (ret > 1) return ret;
  if (sys == AliForwardUtil::kPbPb ||
      sys == AliForwardUtil::kPPb  ||
      sys == AliForwardUtil::kPbp) {
    ret = CheckSNN(beam);
  }
  return ret;
}
//____________________________________________________________________
const char* 
AliForwardUtil::CenterOfMassEnergyString(UShort_t cms)
{
  // 
  // Get a string representation of the center of mass energy per nuclean
  // 
  // Parameters:
  //    cms  Center of mass energy per nucleon
  // 
  // Return:
  //    String representation of the center of mass energy per nuclean
  //
  return Form("%04dGeV", cms);
}
//____________________________________________________________________
Short_t
AliForwardUtil::ParseMagneticField(Float_t v)
{
  // 
  // Parse the magnetic field (in kG) as given by a floating point number
  // 
  // Parameters:
  //    field  Magnetic field in kG 
  // 
  // Return:
  //    Short integer value of magnetic field in kG 
  //
  if (TMath::Abs(v - 5.)  < 1 ) return +5;
  if (TMath::Abs(v - 2.5) < 1 ) return +2;
  if (TMath::Abs(v - 2.5) < 1 ) return -2;
  if (TMath::Abs(v + 5.)  < 1 ) return -5;
  if (TMath::Abs(v) < 1)        return 0;
  ::Warning("ParseMagneticfield", "Magnetic field value: %f not known", v);
  return 999;
}
//____________________________________________________________________
const char* 
AliForwardUtil::MagneticFieldString(Short_t f)
{
  // 
  // Get a string representation of the magnetic field
  // 
  // Parameters:
  //    field Magnetic field in kG
  // 
  // Return:
  //    String representation of the magnetic field
  //
  return Form("%01dkG", f);
}
//_____________________________________________________________________
AliAODEvent* AliForwardUtil::GetAODEvent(AliAnalysisTaskSE* task)
{
  // Check if AOD is the output event
  if (!task) ::Fatal("GetAODEvent", "Null task given, cannot do that");

  AliAODEvent* ret = task->AODEvent();
  if (ret) return ret; 
  
  // Check if AOD is the input event 
  ret = dynamic_cast<AliAODEvent*>(task->InputEvent());
  if (!ret) ::Warning("GetAODEvent", "No AOD event found");
  
  return ret; 
}
//_____________________________________________________________________
UShort_t AliForwardUtil::CheckForAOD()
{
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  if (dynamic_cast<AliAODInputHandler*>(am->GetInputEventHandler())) {
    // ::Info("CheckForAOD", "Found AOD Input handler");
    return 1;
  }
  if (dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler())) {
    // ::Info("CheckForAOD", "Found AOD Output handler");
    return 2;
  }

  ::Warning("CheckForAOD", 
	    "Neither and input nor output AOD handler is specified");
  return 0;
}
//_____________________________________________________________________
Bool_t AliForwardUtil::CheckForTask(const char* clsOrName, Bool_t cls)
{
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  if (!cls) { 
    AliAnalysisTask* t = am->GetTask(clsOrName);
    if (!t) { 
      ::Warning("CheckForTask", "Task %s not found in manager", clsOrName);
      return false;
    }
    ::Info("CheckForTask", "Found task %s", clsOrName);
    return true;
  }
  TClass* dep = gROOT->GetClass(clsOrName);
  if (!dep) { 
    ::Warning("CheckForTask", "Unknown class %s for needed task", clsOrName);
    return false;
  }
  TIter next(am->GetTasks());
  TObject* o = 0;
  while ((o = next())) { 
    if (o->IsA()->InheritsFrom(dep)) {
      ::Info("CheckForTask", "Found task of class %s: %s", 
	     clsOrName, o->GetName());
      return true;
    }
  }
  ::Warning("CheckForTask", "No task of class %s was found", clsOrName);
  return false;
}

//_____________________________________________________________________
TObject* AliForwardUtil::MakeParameter(const char* name, UShort_t value)
{
  TParameter<int>* ret = new TParameter<int>(name, value);
  ret->SetMergeMode('f');
  ret->SetUniqueID(value);
  return ret;
}
//_____________________________________________________________________
TObject* AliForwardUtil::MakeParameter(const char* name, Int_t value)
{
  TParameter<int>* ret = new TParameter<int>(name, value);
  ret->SetMergeMode('f');
  ret->SetUniqueID(value);
  return ret;
}
//_____________________________________________________________________
TObject* AliForwardUtil::MakeParameter(const char* name, ULong_t value)
{
  TParameter<Long_t>* ret = new TParameter<Long_t>(name, value);
  ret->SetMergeMode('f');
  ret->SetUniqueID(value);
  return ret;
}
//_____________________________________________________________________
TObject* AliForwardUtil::MakeParameter(const char* name, Double_t value)
{
  TParameter<double>* ret = new TParameter<double>(name, value);
  // Float_t v = value;
  // UInt_t* tmp = reinterpret_cast<UInt_t*>(&v);
  ret->SetMergeMode('f');
  // ret->SetUniqueID(*tmp);
  return ret;
}
//_____________________________________________________________________
TObject* AliForwardUtil::MakeParameter(const char* name, Bool_t value)
{
  TParameter<bool>* ret = new TParameter<bool>(name, value);
  ret->SetMergeMode('f');
  ret->SetUniqueID(value);
  return ret;
}

//_____________________________________________________________________
void AliForwardUtil::GetParameter(TObject* o, UShort_t& value)
{
  if (!o) return;
  TParameter<int>* p = static_cast<TParameter<int>*>(o);
  if (p->TestBit(BIT(19)))
    value = p->GetVal(); 
  else
    value = o->GetUniqueID();
}
//_____________________________________________________________________
void AliForwardUtil::GetParameter(TObject* o, Int_t& value)
{
  if (!o) return;
  TParameter<int>* p = static_cast<TParameter<int>*>(o);
  if (p->TestBit(BIT(19)))
    value = p->GetVal(); 
  else
    value = o->GetUniqueID();
}
//_____________________________________________________________________
void AliForwardUtil::GetParameter(TObject* o, ULong_t& value)
{
  if (!o) return;
  TParameter<Long_t>* p = static_cast<TParameter<Long_t>*>(o);
  if (p->TestBit(BIT(19)))
    value = p->GetVal(); 
  else
    value = o->GetUniqueID();
}
//_____________________________________________________________________
void AliForwardUtil::GetParameter(TObject* o, Double_t& value)
{
  if (!o) return;
  TParameter<double>* p = static_cast<TParameter<double>*>(o);
  if (p->TestBit(BIT(19)))
    value = p->GetVal(); // o->GetUniqueID();
  else {
    UInt_t  i = o->GetUniqueID();
    Float_t v = *reinterpret_cast<Float_t*>(&i);
    value = v;
  }
}
//_____________________________________________________________________
void AliForwardUtil::GetParameter(TObject* o, Bool_t& value)
{
  if (!o) return;
  TParameter<bool>* p = static_cast<TParameter<bool>*>(o);
  if (p->TestBit(BIT(19)))
    value = p->GetVal(); // o->GetUniqueID();
  else
    value = o->GetUniqueID();
}
  
//_____________________________________________________________________
Double_t AliForwardUtil::GetStripR(Char_t ring, UShort_t strip)
{
  // Get max R of ring
  // 
  // New implementation has only one branch
  const Double_t minR[] = {  4.5213, 15.4 };
  const Double_t maxR[] = { 17.2,    28.0 };
  const Int_t    nStr[] = { 512,     256  };

  Int_t      q       = (ring == 'I' || ring == 'i') ? 0 : 1;  
  Double_t   rad     = maxR[q] - minR[q];
  Double_t   segment = rad / nStr[q];
  Double_t   r       = minR[q] + segment*strip;

  return r;
}

//_____________________________________________________________________
Double_t AliForwardUtil::GetSectorZ(UShort_t det, Char_t ring, UShort_t sec)
{
  Int_t          hybrid  = sec / 2;
  Int_t          q       = (ring == 'I' || ring == 'i') ? 0 : 1;
  Int_t          r       = q == 0 ? 1 : 0;
  const Double_t zs[][2] = { { 320.266+1.5, kInvalidValue }, 
			     {  83.666,     74.966+.5 },
			     { -63.066+.5, -74.966 } };
  if (det > 3 || zs[det-1][q] == kInvalidValue) {
    ::Warning("GetSectorZ", "Unknown sub-detector FMD%d%c", det, ring);
    return kInvalidValue;
  }

  Double_t z = zs[det-1][q];
  switch (det) {
  case 1: if ((hybrid % 2) == 1) z -= .5; break;
  case 2: if ((hybrid % 2) == r) z -= .5; break;
  case 3: if ((hybrid % 2) == q) z -= .5; break;
  }

  return z;
}

//_____________________________________________________________________
Double_t AliForwardUtil::GetSectorPhi(UShort_t d, Char_t ring, UShort_t sec)
{
  UShort_t nSec = (ring == 'I' || ring == 'i') ? 20 : 40;
  Double_t base = float(sec+.5) / nSec * TMath::TwoPi();
  switch (d) {
  case 1:  base += TMath::Pi()/2;  break;
  case 2:  break; 
  case 3:  base = TMath::Pi() - base; break;
  default:
    ::Warning("GetSectorPhi", "Unknown detector %d", d);
    return kInvalidValue;
  }
  if (base < 0)              base += TMath::TwoPi();
  if (base > TMath::TwoPi()) base -= TMath::TwoPi();
  return base;
}
  
//_____________________________________________________________________
Double_t AliForwardUtil::GetEtaFromStrip(UShort_t det, Char_t ring, 
					 UShort_t sec, UShort_t strip, 
					 Double_t zvtx)
{
  // Calculate eta from strip with vertex (redundant with
  // AliESDFMD::Eta but support displaced vertices)
  //
  // Slightly more optimized version that uses less branching 
  
  // Get R of the strip
  Double_t   r     = GetStripR(ring, strip);
  Double_t   z     = GetSectorZ(det, ring, sec);   
  Double_t   theta = TMath::ATan2(r,z-zvtx);
  Double_t   eta   = -1*TMath::Log(TMath::Tan(0.5*theta));
  
  return eta;
}

//_____________________________________________________________________
Bool_t AliForwardUtil::GetXYZ(UShort_t det, Char_t ring,
			      UShort_t sec, UShort_t strip,
			      const TVector3& ip,
			      TVector3& pos)
{
  Double_t   rD      = GetStripR(ring, strip);
  Double_t   phiD    = GetSectorPhi(det,ring, sec);
  Double_t   zD      = GetSectorZ(det, ring, sec);
  if (phiD == kInvalidValue || zD == kInvalidValue) {
    pos.SetXYZ(kInvalidValue,kInvalidValue,kInvalidValue);
    return false;
  }
  Double_t   xD      = rD*TMath::Cos(phiD);
  Double_t   yD      = rD*TMath::Sin(phiD);
  Double_t   iX      = ip.X(); if (iX > 100) iX = 0; // No X
  Double_t   iY      = ip.Y(); if (iY > 100) iY = 0; // No Y
  Double_t   dX      = xD-iX;
  Double_t   dY      = yD-iY;
  Double_t   dZ      = zD-ip.Z();
  pos.SetXYZ(dX, dY, dZ);
  return true;
}
//_____________________________________________________________________
Bool_t AliForwardUtil::GetEtaPhi(UShort_t det, Char_t ring,
				 UShort_t sec, UShort_t strip,
				 const TVector3& ip,
				 Double_t& eta, Double_t& phi)
{
  TVector3 pos;
  if (!GetXYZ(det, ring, sec, strip, ip, pos)) {
    ::Warning("GetEtaPhi", "Invalid position for FMD%d%c[%2d,%3d]=(%f,%f,%f)",
	      det, ring, sec, strip, pos.X(), pos.Y(), pos.Z());    
    eta = kInvalidValue;
    phi = kInvalidValue;
    return false;
  }
  Double_t   r       = TMath::Sqrt(TMath::Power(pos.X(),2)+
				   TMath::Power(pos.Y(),2));
  Double_t   theta   = TMath::ATan2(r, pos.Z());
  Double_t   tant    = TMath::Tan(theta/2);
  if (TMath::Abs(theta) < 1e-9) {
    ::Warning("GetEtaPhi","tan(theta/2)=%f very small", tant);
    eta = kInvalidValue;
    phi = kInvalidValue;
    return false;
  }
  phi = TMath::ATan2(pos.Y(), pos.X());
  eta = -TMath::Log(tant);
  if (phi < 0)              phi += TMath::TwoPi();
  if (phi > TMath::TwoPi()) phi -= TMath::TwoPi();

  return true;
}

//_____________________________________________________________________
Bool_t AliForwardUtil::GetEtaPhiFromStrip(Char_t    r,
					  UShort_t  strip,
					  Double_t& eta, Double_t& phi , 
					  Double_t  ipX, Double_t  ipY)
{
  Double_t rs  = GetStripR(r, strip);
  Double_t sx  = rs*TMath::Cos(phi);
  Double_t sy  = rs*TMath::Sin(phi);
  Double_t dx  = sx-ipX;
  Double_t dy  = sy-ipY;
  Double_t rv  = TMath::Sqrt(TMath::Power(dx,2) + TMath::Power(dy,2));
  Double_t the = 2*TMath::ATan(TMath::Exp(-eta));
  Double_t tth = TMath::Tan(the);
  if (TMath::Abs(tth) < 1e-9) {
    ::Warning("GetEtaPhiFromStrip",
	      "eta=%f -> theta=%f tan(theta)=%f invalid (no change)",
	      eta, the, tth);
    return false;
  }
  Double_t z   = rs / tth;
  // Printf("IP(x,y)=%f,%f S(x,y)=%f,%f D(x,y)=%f,%f R=%f theta=%f "
  //        "tan(theta)=%f z=%f",
  //        ipX, ipY, sx, sy, dx, dy, rv, the, TMath::Tan(the), z);
  eta          = -TMath::Log(TMath::Tan(TMath::ATan2(rv,z)/2));
  phi          = TMath::ATan2(dy,dx);
  if (phi < 0) phi += TMath::TwoPi();

  return true;
}


//_____________________________________________________________________
Double_t AliForwardUtil::GetPhiFromStrip(Char_t ring, UShort_t strip, 
					 Double_t phi,
					 Double_t xvtx, Double_t yvtx)
{
  // Calculate eta from strip with vertex (redundant with
  // AliESDFMD::Eta but support displaced vertices)

  // Unknown x,y -> no change
  if (yvtx > 999 || xvtx > 999) return phi;

  //Get max R of ring
  Double_t r   = GetStripR(ring, strip);
  Double_t amp = TMath::Sqrt(xvtx*xvtx+yvtx*yvtx) / r;
  Double_t pha = (TMath::Abs(yvtx) < 1e-12  ? 0 : TMath::ATan2(xvtx, yvtx));
  Double_t cha = amp * TMath::Cos(phi+pha);
  phi += cha;
  if (phi < 0)              phi += TMath::TwoPi();
  if (phi > TMath::TwoPi()) phi -= TMath::TwoPi();
  return phi;
}
//====================================================================
TAxis*
AliForwardUtil::MakeFullIpZAxis(Int_t nCenter)
{
  TArrayD bins;
  MakeFullIpZAxis(nCenter, bins);
  TAxis* a = new TAxis(bins.GetSize()-1,bins.GetArray());
  return a;
}
//____________________________________________________________________
void
AliForwardUtil::MakeFullIpZAxis(Int_t nCenter, TArrayD& bins)
{
  // Custom vertex axis that will include satellite vertices 
  // Satellite vertices are at k*37.5 where k=-10,-9,...,9,10 
  // Nominal vertices are usually in -10 to 10 and we should have 
  // 10 bins in that range.  That gives us a total of 
  //
  //   10+10+10=30 bins 
  // 
  // or 31 bin boundaries 
  if (nCenter % 2 == 1) 
    // Number of central bins is odd - make it even
    nCenter--;
  const Double_t mCenter = 20;
  const Int_t    nSat    = 10;
  const Int_t    nBins   = 2*nSat + nCenter;
  const Int_t    mBin    = nBins / 2;
  Double_t       dCenter = 2*mCenter / nCenter;
  bins.Set(nBins+1);
  bins[mBin] = 0;
  for (Int_t i = 1; i <= nCenter/2; i++) { 
    // Assign from the middle out 
    Double_t  v  = i * dCenter;
    // Printf("Assigning +/-%7.2f to %3d/%3d", v,mBin-i,mBin+i);
    bins[mBin-i] = -v;
    bins[mBin+i] = +v;
  }
  for (Int_t i = 1; i <= nSat; i++) { 
    Double_t v = (i+.5) * 37.5;
    Int_t    o = nCenter/2+i;
    // Printf("Assigning +/-%7.2f to %3d/%3d", v,mBin-o,mBin+o);
    bins[mBin-o] = -v;
    bins[mBin+o] = +v;
  }
}
//____________________________________________________________________
void 
AliForwardUtil::MakeLogScale(Int_t    nBins, 
			     Int_t    minOrder, 
			     Int_t    maxOrder, 
			     TArrayD& bins)
{
  Double_t dO = Double_t(maxOrder-minOrder) / nBins; 
  bins.Set(nBins+1);
  for (Int_t i = 0; i <= nBins; i++) bins[i] = TMath::Power(10, i * dO+minOrder);
}

//____________________________________________________________________
void 
AliForwardUtil::PrintTask(const TObject& o)
{
  Int_t ind = gROOT->GetDirLevel();
  if (ind > 0) 
    // Print indention 
    std::cout << std::setfill(' ') << std::setw(ind) << " " << std::flush;

  TString t = TString::Format("%s %s", o.GetName(), o.ClassName());
  const Int_t maxN = 75;
  std::cout << "--- " << t << " " << std::setfill('-') 
	    << std::setw(maxN-ind-5-t.Length()) << "-" << std::endl;
}
//____________________________________________________________________
void
AliForwardUtil::PrintName(const char* name)
{
  Int_t ind = gROOT->GetDirLevel();
  if (ind > 0) 
    // Print indention 
    std::cout << std::setfill(' ') << std::setw(ind) << " " << std::flush;
    
  // Now print field name 
  const Int_t maxN  = 29;
  Int_t       width = maxN - ind;
  TString     n(name);
  if (n.Length() > width-1) {
    // Truncate the string, and put in "..."
    n.Remove(width-4);
    n.Append("...");
  }
  n.Append(":");
  std::cout << std::setfill(' ') << std::left << std::setw(width) 
	    << n << std::right << std::flush;
}
//____________________________________________________________________
void
AliForwardUtil::PrintField(const char* name, const char* value, ...)
{
  PrintName(name);

  // Now format the field value 
  va_list ap;
  va_start(ap, value);
  static char buf[512];
  vsnprintf(buf, 511, value, ap);
  buf[511] = '\0';
  va_end(ap);

  std::cout << buf << std::endl;
}


//====================================================================
Float_t AliForwardUtil::GetCentrality(const AliVEvent& event, 
				      const TString&   method, 
				      Int_t&           qual, 
				      Bool_t           verbose)
{
  qual = 0xFFFF;
  Float_t cent = GetCentralityMult(event,method,qual,verbose);
  if (qual < 0xFFF) return cent;
  
  // No selection object found, try compat
  return GetCentralityCompat(event,method,qual,verbose);  
}
//____________________________________________________________________
Float_t AliForwardUtil::GetCentralityMult(const AliVEvent& event, 
					  const TString&   method, 
					  Int_t&           qual, 
					  Bool_t           verbose)
{
  TObject* o = event.FindListObject("MultSelection");
  if (!o) {
    if (verbose) 
      ::Warning("AliForwardUtil::GetCentralityMult",
		"No MultSelection object found in event");
    return -1;
  }
  AliMultSelection* sel = static_cast<AliMultSelection*>(o);  
  if (!sel->GetEstimatorList() ||
      sel->GetEstimatorList()->GetEntries() <= 0){
    if (verbose) {
      ::Warning("AliForwardUtil::GetCentralityMult",
		"No list of estimators, falling back to compat");
      sel->PrintInfo();
    }
    return -1;
  }
  AliMultEstimator* est = sel->GetEstimator(method);
  // if (verbose) sel->GetEstimatorList()->ls();
  if (!est) {
    if (verbose) {
      ::Warning("AliForwardUtil::GetCentralityMult",
		"Unknown estimator: %s", method.Data());
      sel->GetEstimatorList()->Print();
    }
    return -1;
  }
  // 198 -> 1: beyond anchor 
  // 199 -> 2: no calib 
  // 200 -> 3: not desired trigger
  // 201 -> 4: not INEL>0 with tracklets
  // 202 -> 5: vertex Z not within 10cm
  // 203 -> 6: tagged as pileup (SPD)
  // 204 -> 7: inconsistent SPD/tracking vertex
  // 205 -> 8: rejected by tracklets-vs-clusters
  Float_t cent = est->GetPercentile();
  qual         = sel->GetEvSelCode();
  if (qual == AliMultSelectionCuts::kNoCalib) cent = -1;
  if (verbose)
    ::Info("AliForwardUtil::GetCentralityMult",
	   "Got centrality %5.1f%% (%d)", /*" - old %5.1f%% (%d)",*/
	   cent, qual/*, old, oldQual*/);
  return cent;

}

//____________________________________________________________________
Float_t AliForwardUtil::GetCentralityCompat(const AliVEvent& event,
					    const TString&   method,
					    Int_t&           qual,
					    Bool_t           verbose)
{
  if (event.IsA()->InheritsFrom(AliESDEvent::Class()))
    return GetCentralityCompat(static_cast<const AliESDEvent&>(event),
			       method,qual,verbose);
  if (event.IsA()->InheritsFrom(AliAODEvent::Class()))
    return GetCentralityCompat(static_cast<const AliAODEvent&>(event),
			       method,qual,verbose);
  return -1;
}
//____________________________________________________________________
Float_t AliForwardUtil::GetCentralityCompat(const AliESDEvent& event,
					    const TString&     method,
					    Int_t&             qual,
					    Bool_t             verbose)
{
  AliCentrality* centObj = const_cast<AliESDEvent&>(event).GetCentrality();
  if (!centObj) { 
    if (verbose) 
      ::Warning("AliForwardUtil::GetCentralityCompat",
		"No centrality object found in ESD");
    return -1;
  }
  Float_t cent = centObj->GetCentralityPercentileUnchecked(method);  
  qual         = centObj->GetQuality();
  if (verbose)
    ::Info("AliForwardUtil::GetCentralityCompat<ESD>",
	   "Got centrality %5.1f%% (%d)", cent, qual);  
  if (qual & 0x1)
    qual = AliMultSelectionCuts::kRejVzCut;
  if (qual & 0x2)
    qual = AliMultSelectionCuts::kRejTrackletsVsClusters;
  if (qual & 0x4)
    qual = AliMultSelectionCuts::kRejConsistencySPDandTrackVertices;
  if (qual & 0x8)
    qual = AliMultSelectionCuts::kRejTrackletsVsClusters;
  return cent;
}
//____________________________________________________________________
Float_t AliForwardUtil::GetCentralityCompat(const AliAODEvent& event,
					    const TString&     method,
					    Int_t&             qual,
					    Bool_t             verbose)
{
  AliAODHeader* hdr = dynamic_cast<AliAODHeader*>(event.GetHeader());
  if (!hdr) {
    if (verbose)
      ::Warning("AliForwardUtil::GetCentralityCompat","Not a standard AOD");
    return -1;
  }
  AliCentrality* cP = hdr->GetCentralityP();
  if (!cP) {
    if (verbose)
      ::Warning("AliForwardUtil::GetCentralityCompat",
		"No centrality found in AOD");
    return -1;
  }
  Float_t cent = cP->GetCentralityPercentile(method);
  qual         = cP->GetQuality();
  if (qual & 0x1) qual = AliMultSelectionCuts::kRejVzCut;
  if (qual & 0x2) qual = AliMultSelectionCuts::kRejTrackletsVsClusters;
  if (qual & 0x4) qual = AliMultSelectionCuts::kRejConsistencySPDandTrackVertices;
  if (qual & 0x8) qual = AliMultSelectionCuts::kRejTrackletsVsClusters;
		    
  if (verbose)
    ::Info("AliForwardUtil::GetCentralityCompat<AOD>",
	   "Got centrality %5.1f%% (%d)", cent, qual);
  return cent;
}
//====================================================================
AliForwardUtil::Histos::~Histos()
{
  // 
  // Destructor
  //
}

//____________________________________________________________________
void
AliForwardUtil::Histos::Delete(Option_t* opt)
{
  if (fFMD1i) delete fFMD1i;
  if (fFMD2i) delete fFMD2i;
  if (fFMD2o) delete fFMD2o;
  if (fFMD3i) delete fFMD3i;
  if (fFMD3o) delete fFMD3o;
  fFMD1i = 0;
  fFMD2i = 0;
  fFMD2o = 0;
  fFMD3i = 0;
  fFMD3o = 0;
  TObject::Delete(opt);
}

//____________________________________________________________________
TH2D*
AliForwardUtil::Histos::Make(UShort_t d, Char_t r, const TAxis& etaAxis)
{
  // 
  // Make a histogram 
  // 
  // Parameters:
  //    d        Detector
  //    r        Ring 
  //    etaAxis  Eta axis to use
  // 
  // Return:
  //    Newly allocated histogram 
  //
  Int_t ns = (r == 'I' || r == 'i') ? 20 : 40;
  TH2D* hist = 0;
  if (etaAxis.GetXbins() && etaAxis.GetXbins()->GetArray())
    hist = new TH2D(Form("FMD%d%c_cache", d, r), 
		    Form("FMD%d%c cache", d, r),
		    etaAxis.GetNbins(), etaAxis.GetXbins()->GetArray(), 
		    ns, 0, TMath::TwoPi());
  else
    hist = new TH2D(Form("FMD%d%c_cache", d, r), 
		    Form("FMD%d%c cache", d, r),
		    etaAxis.GetNbins(), etaAxis.GetXmin(), 
		    etaAxis.GetXmax(), ns, 0, TMath::TwoPi());
  hist->SetXTitle("#eta");
  hist->SetYTitle("#phi [radians]");
  hist->SetZTitle("d^{2}N_{ch}/d#etad#phi");
  hist->Sumw2();
  hist->SetDirectory(0);

  return hist;
}
//____________________________________________________________________
void
AliForwardUtil::Histos::RebinEta(TH2D* hist, const TAxis& etaAxis)
{
  TAxis* xAxis = hist->GetXaxis();
  if (etaAxis.GetXbins() && etaAxis.GetXbins()->GetArray())
    xAxis->Set(etaAxis.GetNbins(), etaAxis.GetXbins()->GetArray());
  else
    xAxis->Set(etaAxis.GetNbins(), etaAxis.GetXmin(), etaAxis.GetXmax());
  hist->Rebuild();
}
  

//____________________________________________________________________
void
AliForwardUtil::Histos::Init(const TAxis& etaAxis)
{
  // 
  // Initialize the object 
  // 
  // Parameters:
  //    etaAxis Eta axis to use 
  //
  fFMD1i = Make(1, 'I', etaAxis);
  fFMD2i = Make(2, 'I', etaAxis);
  fFMD2o = Make(2, 'O', etaAxis);
  fFMD3i = Make(3, 'I', etaAxis);
  fFMD3o = Make(3, 'O', etaAxis);
}
//____________________________________________________________________
void
AliForwardUtil::Histos::ReInit(const TAxis& etaAxis)
{
  // 
  // Initialize the object 
  // 
  // Parameters:
  //    etaAxis Eta axis to use 
  //
  if (!fFMD1i) fFMD1i = Make(1, 'i', etaAxis); else RebinEta(fFMD1i, etaAxis);
  if (!fFMD2i) fFMD2i = Make(2, 'i', etaAxis); else RebinEta(fFMD2i, etaAxis);
  if (!fFMD2o) fFMD2o = Make(2, 'o', etaAxis); else RebinEta(fFMD2o, etaAxis);
  if (!fFMD3i) fFMD3i = Make(3, 'i', etaAxis); else RebinEta(fFMD3i, etaAxis);
  if (!fFMD3o) fFMD3o = Make(3, 'o', etaAxis); else RebinEta(fFMD3o, etaAxis);
}

//____________________________________________________________________
void
AliForwardUtil::Histos::Clear(Option_t* option)
{
  // 
  // Clear data 
  // 
  // Parameters:
  //    option Not used 
  //
  if (fFMD1i) { fFMD1i->Reset(option); fFMD1i->ResetBit(kSkipRing); }
  if (fFMD2i) { fFMD2i->Reset(option); fFMD2i->ResetBit(kSkipRing); }
  if (fFMD2o) { fFMD2o->Reset(option); fFMD2o->ResetBit(kSkipRing); }
  if (fFMD3i) { fFMD3i->Reset(option); fFMD3i->ResetBit(kSkipRing); }
  if (fFMD3o) { fFMD3o->Reset(option); fFMD3o->ResetBit(kSkipRing); }
}

//____________________________________________________________________
TH2D*
AliForwardUtil::Histos::Get(UShort_t d, Char_t r) const
{
  // 
  // Get the histogram for a particular detector,ring
  // 
  // Parameters:
  //    d Detector 
  //    r Ring 
  // 
  // Return:
  //    Histogram for detector,ring or nul 
  //
  switch (d) { 
  case 1: return fFMD1i;
  case 2: return (r == 'I' || r == 'i' ? fFMD2i : fFMD2o);
  case 3: return (r == 'I' || r == 'i' ? fFMD3i : fFMD3o);
  }
  return 0;
}
//====================================================================
TList*
AliForwardUtil::RingHistos::DefineOutputList(TList* d) const
{
  // 
  // Define the outout list in @a d
  // 
  // Parameters:
  //    d Where to put the output list
  // 
  // Return:
  //    Newly allocated TList object or null
  //
  if (!d) return 0;
  TList* list = new TList;
  list->SetOwner();
  list->SetName(fName.Data());
  d->Add(list);
  return list;
}
//____________________________________________________________________
TList*
AliForwardUtil::RingHistos::GetOutputList(const TList* d) const
{
  // 
  // Get our output list from the container @a d
  // 
  // Parameters:
  //    d where to get the output list from 
  // 
  // Return:
  //    The found TList or null
  //
  if (!d) return 0;
  TList* list = static_cast<TList*>(d->FindObject(fName.Data()));
  return list;
}

//____________________________________________________________________
TH1*
AliForwardUtil::RingHistos::GetOutputHist(const TList* d, const char* name) const
{
  // 
  // Find a specific histogram in the source list @a d
  // 
  // Parameters:
  //    d     (top)-container 
  //    name  Name of histogram
  // 
  // Return:
  //    Found histogram or null
  //
  return static_cast<TH1*>(d->FindObject(name));
}


//====================================================================
AliForwardUtil::DebugGuard::DebugGuard(Int_t lvl, Int_t msgLvl, 
				       const char* format, ...)
  : fMsg("")
{
  if (lvl < msgLvl) return; 
  va_list ap;
  va_start(ap, format);
  Format(fMsg, format, ap);
  va_end(ap);
  Output(+1, fMsg);
}
//____________________________________________________________________
AliForwardUtil::DebugGuard::~DebugGuard()
{
  if (fMsg.IsNull()) return;
  Output(-1, fMsg);
}
//____________________________________________________________________
void
AliForwardUtil::DebugGuard::Message(Int_t lvl, Int_t msgLvl, 
				    const char* format, ...)
{
  if (lvl < msgLvl) return; 
  TString msg;
  va_list ap;
  va_start(ap, format);
  Format(msg, format, ap);
  va_end(ap);
  Output(0, msg);
}

//____________________________________________________________________
void
AliForwardUtil::DebugGuard::Format(TString& out, const char* format, va_list ap)
{
  static char buf[512];
  Int_t n = gROOT->GetDirLevel() + 2;
  for (Int_t i = 0; i < n; i++) buf[i] = ' ';
  vsnprintf(&(buf[n]), 511-n, format, ap);
  buf[511] = '\0';
  out = buf;  
}
//____________________________________________________________________
void
AliForwardUtil::DebugGuard::Output(int in, TString& msg)
{
  msg[0] = (in > 0 ? '>' :  in < 0 ? '<' : '=');
  AliLog::Message(AliLog::kInfo, msg, 0, 0, "PWGLF/forward", 0, 0);
  if      (in > 0) gROOT->IncreaseDirLevel();
  else if (in < 0) gROOT->DecreaseDirLevel();
}
//====================================================================
AliForwardUtil::SuppressGuard::SuppressGuard(Int_t lvl)
  : save(gErrorIgnoreLevel)
{
  gErrorIgnoreLevel = lvl;
  AliLog::SetModuleDebugLevel("ROOT", lvl);
  AliLog::SetGlobalDebugLevel(lvl);
}
//____________________________________________________________________
AliForwardUtil::SuppressGuard::~SuppressGuard()
{
  gErrorIgnoreLevel = save;
  AliLog::SetModuleDebugLevel("ROOT", save);
  AliLog::SetGlobalDebugLevel(save);
}

//
// EOF
//
