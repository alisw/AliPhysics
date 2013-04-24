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
#include <TParameter.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TMath.h>
#include <TError.h>
#include <TROOT.h>

//====================================================================
ULong_t AliForwardUtil::AliROOTRevision()
{
#ifdef ALIROOT_SVN_REVISION
  return ALIROOT_SVN_REVISION;
#else 
  return 0;
#endif
}
//____________________________________________________________________
ULong_t AliForwardUtil::AliROOTBranch()
{
#ifdef ALIROOT_SVN_BRANCH
  static ULong_t ret = 0;
  if (ret != 0) return ret;
  
  TString str(ALIROOT_SVN_BRANCH);
  if (str[0] == 'v') str.Remove(0,1);
  if (str.EqualTo("trunk")) return ret = 0xFFFFFFFF;

  TObjArray*   tokens = str.Tokenize("-");
  TObjString*  pMajor = static_cast<TObjString*>(tokens->At(0));
  TObjString*  pMinor = static_cast<TObjString*>(tokens->At(1));
  TObjString*  pRelea = static_cast<TObjString*>(tokens->At(2));
  TObjString* pAn     = (tokens->GetEntries() > 3 ? 
    static_cast<TObjString*>(tokens->At(3)) : 0);
  TString sMajor = pMajor->String().Strip(TString::kLeading, '0');
  TString sMinor = pMinor->String().Strip(TString::kLeading, '0');
  TString sRelea = pRelea->String().Strip(TString::kLeading, '0');
  
  ret = (((sMajor.Atoi() & 0xFF) << 12) |
    ((sMinor.Atoi() & 0xFF) <<  8) |
    ((sRelea.Atoi() & 0xFF) <<  4) |
    (pAn ? 0xAA : 0));
  
  return ret;
#else 
  return 0;
#endif
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
  if (s.Contains("a-p")   || s.Contains("ap"))    return AliForwardUtil::kPPb;
  if (s.Contains("p-p")   || s.Contains("pp"))    return AliForwardUtil::kPP; 
  if (s.Contains("pb-pb") || s.Contains("pbpb"))  return AliForwardUtil::kPbPb;
  if (s.Contains("a-a")   || s.Contains("aa"))    return AliForwardUtil::kPbPb;
  return AliForwardUtil::kUnknown;
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
  // - kPP -> "pp"
  // - kPbPb -> "PbPb" 
  // - anything else gives "unknown"
  // 
  // Return:
  //    String representation of the collision system 
  //
  switch (sys) { 
  case AliForwardUtil::kPP:   return "pp";
  case AliForwardUtil::kPbPb: return "PbPb";
  case AliForwardUtil::kPPb:  return "pPb";
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
  if (TMath::Abs(energy - 900.)   < 10)  return 900;
  if (TMath::Abs(energy - 2400.)  < 10)  return 2400;
  if (TMath::Abs(energy - 2750.)  < 20)  return 2750;
  if (TMath::Abs(energy - 4400.)  < 10)  return 4400;
  if (TMath::Abs(energy - 5022.)  < 10)  return 5023;
  if (TMath::Abs(energy - 5500.)  < 40)  return 5500;
  if (TMath::Abs(energy - 7000.)  < 10)  return 7000;
  if (TMath::Abs(energy - 8000.)  < 10)  return 8000;
  if (TMath::Abs(energy - 10000.) < 10)  return 10000;
  if (TMath::Abs(energy - 14000.) < 10)  return 14000;
  return 0;
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
  if (TMath::Abs(v - 5.) < 1 ) return +5;
  if (TMath::Abs(v + 5.) < 1 ) return -5;
  if (TMath::Abs(v) < 1)       return 0;
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
    ::Info("CheckForAOD", "Found AOD Input handler");
    return 1;
  }
  if (dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler())) {
    ::Info("CheckForAOD", "Found AOD Output handler");
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
TObject* AliForwardUtil::MakeParameter(const Char_t* name, UShort_t value)
{
  TParameter<int>* ret = new TParameter<int>(name, value);
  ret->SetUniqueID(value);
  return ret;
}
//_____________________________________________________________________
TObject* AliForwardUtil::MakeParameter(const Char_t* name, Int_t value)
{
  TParameter<int>* ret = new TParameter<int>(name, value);
  ret->SetUniqueID(value);
  return ret;
}
//_____________________________________________________________________
TObject* AliForwardUtil::MakeParameter(const Char_t* name, ULong_t value)
{
  TParameter<Long_t>* ret = new TParameter<Long_t>(name, value);
  ret->SetUniqueID(value);
  return ret;
}
//_____________________________________________________________________
TObject* AliForwardUtil::MakeParameter(const Char_t* name, Double_t value)
{
  TParameter<double>* ret = new TParameter<double>(name, value);
  Float_t v = value;
  UInt_t* tmp = reinterpret_cast<UInt_t*>(&v);
  ret->SetUniqueID(*tmp);
  return ret;
}
//_____________________________________________________________________
TObject* AliForwardUtil::MakeParameter(const Char_t* name, Bool_t value)
{
  TParameter<bool>* ret = new TParameter<bool>(name, value);
  ret->SetUniqueID(value);
  return ret;
}

//_____________________________________________________________________
void AliForwardUtil::GetParameter(TObject* o, UShort_t& value)
{
  if (!o) return;
  value = o->GetUniqueID();
}
//_____________________________________________________________________
void AliForwardUtil::GetParameter(TObject* o, Int_t& value)
{
  if (!o) return;
  value = o->GetUniqueID();
}
//_____________________________________________________________________
void AliForwardUtil::GetParameter(TObject* o, ULong_t& value)
{
  if (!o) return;
  value = o->GetUniqueID();
}
//_____________________________________________________________________
void AliForwardUtil::GetParameter(TObject* o, Double_t& value)
{
  if (!o) return;
  UInt_t  i = o->GetUniqueID();
  Float_t v = *reinterpret_cast<Float_t*>(&i);
  value = v;
}
//_____________________________________________________________________
void AliForwardUtil::GetParameter(TObject* o, Bool_t& value)
{
  if (!o) return;
  value = o->GetUniqueID();
}
  
//_____________________________________________________________________
Double_t AliForwardUtil::GetStripR(Char_t ring, UShort_t strip)
{
  //Get max R of ring
  const Double_t iMinR = 4.5213;
  const Double_t iMaxR = 17.2;
  const Double_t oMinR = 15.4;
  const Double_t oMaxR = 28.0;
  
  Double_t   minR    = (ring == 'I' || ring == 'i') ? iMinR : oMinR;
  Double_t   maxR    = (ring == 'I' || ring == 'i') ? iMaxR : oMaxR;
  Double_t   nStrips = (ring == 'I' || ring == 'i') ? 512   : 256;
  Double_t   rad     =  maxR - minR;
  Double_t   segment = rad / nStrips;
  Double_t   r       =  minR + segment*strip;

  return r;
}

//_____________________________________________________________________
Double_t AliForwardUtil::GetEtaFromStrip(UShort_t det, Char_t ring, 
					 UShort_t sec, UShort_t strip, 
					 Double_t zvtx)
{
  // Calculate eta from strip with vertex (redundant with
  // AliESDFMD::Eta but support displaced vertices)
  
  //Get max R of ring
  Double_t   r         = GetStripR(ring, strip);
  Int_t      hybrid    = sec / 2;
  Bool_t     inner     = (ring == 'I' || ring == 'i');
  Double_t   z         = 0;
  switch (det) { 
  case 1: z = 320.266;                     break;
  case 2: z = (inner ?  83.666 :  74.966); break;
  case 3: z = (inner ? -63.066 : -74.966); break; 
  default: return -999999;
  }
  if ((hybrid % 2) == 0) z -= .5;
  
  Double_t   theta = TMath::ATan2(r,z-zvtx);
  Double_t   eta   = -1*TMath::Log(TMath::Tan(0.5*theta));
  
  return eta;
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
Int_t    AliForwardUtil::fgConvolutionSteps  = 100;
Double_t AliForwardUtil::fgConvolutionNSigma = 5;
namespace {
  // 
  // The shift of the most probable value for the ROOT function TMath::Landau 
  //
  const Double_t  mpshift  = -0.22278298;
  // 
  // Integration normalisation 
  //
  const Double_t  invSq2pi = 1. / TMath::Sqrt(2*TMath::Pi());

  // 
  // Utility function to use in TF1 defintition 
  //
  Double_t landauGaus1(Double_t* xp, Double_t* pp) 
  {
    Double_t x        = xp[0];
    Double_t constant = pp[AliForwardUtil::ELossFitter::kC];
    Double_t delta    = pp[AliForwardUtil::ELossFitter::kDelta];
    Double_t xi       = pp[AliForwardUtil::ELossFitter::kXi];
    Double_t sigma    = pp[AliForwardUtil::ELossFitter::kSigma];
    Double_t sigmaN   = pp[AliForwardUtil::ELossFitter::kSigmaN];

    return constant * AliForwardUtil::LandauGaus(x, delta, xi, sigma, sigmaN);
  }

  Double_t landauGausComposite(Double_t* xp, Double_t* pp)
  {
    Double_t x           = xp[0];
    Double_t cP          = pp[AliForwardUtil::ELossFitter::kC];
    Double_t deltaP      = pp[AliForwardUtil::ELossFitter::kDelta];
    Double_t xiP         = pp[AliForwardUtil::ELossFitter::kXi];
    Double_t sigmaP      = pp[AliForwardUtil::ELossFitter::kSigma];
    Double_t cS          = pp[AliForwardUtil::ELossFitter::kSigma+1];
    Double_t deltaS      = pp[AliForwardUtil::ELossFitter::kSigma+2];
    Double_t xiS         = pp[AliForwardUtil::ELossFitter::kSigma+3];
    Double_t sigmaS      = pp[AliForwardUtil::ELossFitter::kSigma+4];

    return (cP * AliForwardUtil::LandauGaus(x,deltaP,xiP,sigmaP,0) + 
	    cS * AliForwardUtil::LandauGaus(x,deltaS,xiS,sigmaS,0));
  }
    
  // 
  // Utility function to use in TF1 defintition 
  //
  Double_t landauGausN(Double_t* xp, Double_t* pp) 
  {
    Double_t  x        = xp[0];
    Double_t constant  = pp[AliForwardUtil::ELossFitter::kC];
    Double_t delta     = pp[AliForwardUtil::ELossFitter::kDelta];
    Double_t xi        = pp[AliForwardUtil::ELossFitter::kXi];
    Double_t sigma     = pp[AliForwardUtil::ELossFitter::kSigma];
    Double_t sigmaN    = pp[AliForwardUtil::ELossFitter::kSigmaN];
    Int_t     n        = Int_t(pp[AliForwardUtil::ELossFitter::kN]);
    Double_t* a        = &(pp[AliForwardUtil::ELossFitter::kA]);

    return constant * AliForwardUtil::NLandauGaus(x, delta, xi, sigma, sigmaN,
						  n, a);
  }
  // 
  // Utility function to use in TF1 defintition 
  //
  Double_t landauGausI(Double_t* xp, Double_t* pp) 
  {
    Double_t x         = xp[0];
    Double_t constant  = pp[AliForwardUtil::ELossFitter::kC];
    Double_t delta     = pp[AliForwardUtil::ELossFitter::kDelta];
    Double_t xi        = pp[AliForwardUtil::ELossFitter::kXi];
    Double_t sigma     = pp[AliForwardUtil::ELossFitter::kSigma];
    Double_t sigmaN    = pp[AliForwardUtil::ELossFitter::kSigmaN];
    Int_t    i         = Int_t(pp[AliForwardUtil::ELossFitter::kN]);

    return constant * AliForwardUtil::ILandauGaus(x,delta,xi,sigma,sigmaN,i);
  }


}
//____________________________________________________________________
Double_t 
AliForwardUtil::Landau(Double_t x, Double_t delta, Double_t xi)
{
  // 
  // Calculate the shifted Landau
  // @f[
  //    f'_{L}(x;\Delta,\xi) = f_L(x;\Delta+0.22278298\xi)
  // @f]
  //
  // where @f$ f_{L}@f$ is the ROOT implementation of the Landau
  // distribution (known to have @f$ \Delta_{p}=-0.22278298@f$ for
  // @f$\Delta=0,\xi=1@f$. 
  // 
  // Parameters:
  //    x      Where to evaluate @f$ f'_{L}@f$ 
  //    delta  Most probable value 
  //    xi     The 'width' of the distribution 
  //
  // Return:
  //    @f$ f'_{L}(x;\Delta,\xi) @f$
  //
  return TMath::Landau(x, delta - xi * mpshift, xi);
}
//____________________________________________________________________
Double_t 
AliForwardUtil::LandauGaus(Double_t x, Double_t delta, Double_t xi,
			   Double_t sigma, Double_t sigmaN)
{
  // 
  // Calculate the value of a Landau convolved with a Gaussian 
  // 
  // @f[ 
  // f(x;\Delta,\xi,\sigma') = \frac{1}{\sigma' \sqrt{2 \pi}}
  //    \int_{-\infty}^{+\infty} d\Delta' f'_{L}(x;\Delta',\xi)
  //    \exp{-\frac{(\Delta-\Delta')^2}{2\sigma'^2}}
  // @f]
  // 
  // where @f$ f'_{L}@f$ is the Landau distribution, @f$ \Delta@f$ the
  // energy loss, @f$ \xi@f$ the width of the Landau, and 
  // @f$ \sigma'^2=\sigma^2-\sigma_n^2 @f$.  Here, @f$\sigma@f$ is the
  // variance of the Gaussian, and @f$\sigma_n@f$ is a parameter modelling 
  // noise in the detector.  
  //
  // Note that this function uses the constants fgConvolutionSteps and
  // fgConvolutionNSigma
  // 
  // References: 
  //  - <a href="http://dx.doi.org/10.1016/0168-583X(84)90472-5">Nucl.Instrum.Meth.B1:16</a>
  //  - <a href="http://dx.doi.org/10.1103/PhysRevA.28.615">Phys.Rev.A28:615</a>
  //  - <a href="http://root.cern.ch/root/htmldoc/tutorials/fit/langaus.C.html">ROOT implementation</a>
  // 
  // Parameters:
  //    x         where to evaluate @f$ f@f$
  //    delta     @f$ \Delta@f$ of @f$ f(x;\Delta,\xi,\sigma')@f$
  //    xi        @f$ \xi@f$ of @f$ f(x;\Delta,\xi,\sigma')@f$
  //    sigma     @f$ \sigma@f$ of @f$\sigma'^2=\sigma^2-\sigma_n^2 @f$
  //    sigma_n   @f$ \sigma_n@f$ of @f$\sigma'^2=\sigma^2-\sigma_n^2 @f$
  // 
  // Return:
  //    @f$ f@f$ evaluated at @f$ x@f$.  
  //
  Double_t deltap = delta - xi * mpshift;
  Double_t sigma2 = sigmaN*sigmaN + sigma*sigma;
  Double_t sigma1 = sigmaN == 0 ? sigma : TMath::Sqrt(sigma2);
  Double_t xlow   = x - fgConvolutionNSigma * sigma1;
  Double_t xhigh  = x + fgConvolutionNSigma * sigma1;
  Double_t step   = (xhigh - xlow) / fgConvolutionSteps;
  Double_t sum    = 0;
  
  for (Int_t i = 0; i <= fgConvolutionSteps/2; i++) { 
    Double_t x1 = xlow  + (i - .5) * step;
    Double_t x2 = xhigh - (i - .5) * step;
    
    sum += TMath::Landau(x1, deltap, xi, kTRUE) * TMath::Gaus(x, x1, sigma1);
    sum += TMath::Landau(x2, deltap, xi, kTRUE) * TMath::Gaus(x, x2, sigma1);
  }
  return step * sum * invSq2pi / sigma1;
}

//____________________________________________________________________
Double_t 
AliForwardUtil::ILandauGaus(Double_t x, Double_t delta, Double_t xi, 
			    Double_t sigma, Double_t sigmaN, Int_t i)
{
  // 
  // Evaluate 
  // @f[ 
  //    f_i(x;\Delta,\xi,\sigma') = f(x;\Delta_i,\xi_i,\sigma_i')
  // @f] 
  // corresponding to @f$ i@f$ particles i.e., with the substitutions 
  // @f{eqnarray*}{ 
  //    \Delta    \rightarrow \Delta_i    &=& i(\Delta + \xi\log(i))
  //    \xi       \rightarrow \xi_i       &=& i \xi
  //    \sigma    \rightarrow \sigma_i    &=& \sqrt{i}\sigma
  //    \sigma'^2 \rightarrow \sigma_i'^2 &=& \sigma_n^2 + \sigma_i^2
  // @f} 
  // 
  // Parameters:
  //    x        Where to evaluate 
  //    delta    @f$ \Delta@f$ 
  //    xi       @f$ \xi@f$ 
  //    sigma    @f$ \sigma@f$ 
  //    sigma_n  @f$ \sigma_n@f$
  //    i        @f$ i @f$
  // 
  // Return:
  //    @f$ f_i @f$ evaluated
  //  
  Double_t deltaI =  (i == 1 ? delta : i * (delta + xi * TMath::Log(i)));
  Double_t xiI    =  i * xi;
  Double_t sigmaI =  (i == 1 ? sigma : TMath::Sqrt(Double_t(i))*sigma);
  if (sigmaI < 1e-10) { 
    // Fall back to landau 
    return AliForwardUtil::Landau(x, deltaI, xiI);
  }
  return AliForwardUtil::LandauGaus(x, deltaI, xiI, sigmaI, sigmaN);
}

//____________________________________________________________________
Double_t 
AliForwardUtil::IdLandauGausdPar(Double_t x, 
				 UShort_t par,   Double_t dPar, 
				 Double_t delta, Double_t xi, 
				 Double_t sigma, Double_t sigmaN, 
				 Int_t    i)
{
  // 
  // Numerically evaluate 
  // @f[ 
  //    \left.\frac{\partial f_i}{\partial p_i}\right|_{x}
  // @f] 
  // where @f$ p_i@f$ is the @f$ i^{\mbox{th}}@f$ parameter.  The mapping 
  // of the parameters is given by 
  //
  // - 0: @f$\Delta@f$ 
  // - 1: @f$\xi@f$ 
  // - 2: @f$\sigma@f$ 
  // - 3: @f$\sigma_n@f$ 
  //
  // This is the partial derivative with respect to the parameter of
  // the response function corresponding to @f$ i@f$ particles i.e.,
  // with the substitutions
  // @f[ 
  //    \Delta    \rightarrow \Delta_i    = i(\Delta + \xi\log(i))
  //    \xi       \rightarrow \xi_i       = i \xi
  //    \sigma    \rightarrow \sigma_i    = \sqrt{i}\sigma
  //    \sigma'^2 \rightarrow \sigma_i'^2 = \sigma_n^2 + \sigma_i^2
  // @f] 
  // 
  // Parameters:
  //    x        Where to evaluate 
  //    ipar     Parameter number 
  //    dp       @f$ \epsilon\delta p_i@f$ for some value of @f$\epsilon@f$
  //    delta    @f$ \Delta@f$ 
  //    xi       @f$ \xi@f$ 
  //    sigma    @f$ \sigma@f$ 
  //    sigma_n  @f$ \sigma_n@f$
  //    i        @f$ i@f$
  // 
  // Return:
  //    @f$ f_i@f$ evaluated
  //  
  if (dPar == 0) return 0;
  Double_t dp      = dPar;
  Double_t d2      = dPar / 2;
  Double_t deltaI  =  i * (delta + xi * TMath::Log(i));
  Double_t xiI     =  i * xi;
  Double_t si      =  TMath::Sqrt(Double_t(i));
  Double_t sigmaI  =  si*sigma;
  Double_t y1      = 0;
  Double_t y2      = 0;
  Double_t y3      = 0;
  Double_t y4      = 0;
  switch (par) {
  case 0: 
    y1 = ILandauGaus(x, deltaI+i*dp, xiI, sigmaI, sigmaN, i);
    y2 = ILandauGaus(x, deltaI+i*d2, xiI, sigmaI, sigmaN, i);
    y3 = ILandauGaus(x, deltaI-i*d2, xiI, sigmaI, sigmaN, i);
    y4 = ILandauGaus(x, deltaI-i*dp, xiI, sigmaI, sigmaN, i);
    break;
  case 1: 
    y1 = ILandauGaus(x, deltaI, xiI+i*dp, sigmaI, sigmaN, i);
    y2 = ILandauGaus(x, deltaI, xiI+i*d2, sigmaI, sigmaN, i);
    y3 = ILandauGaus(x, deltaI, xiI-i*d2, sigmaI, sigmaN, i);
    y4 = ILandauGaus(x, deltaI, xiI-i*dp, sigmaI, sigmaN, i);
    break;
  case 2: 
    y1 = ILandauGaus(x, deltaI, xiI, sigmaI+si*dp, sigmaN, i);
    y2 = ILandauGaus(x, deltaI, xiI, sigmaI+si*d2, sigmaN, i);
    y3 = ILandauGaus(x, deltaI, xiI, sigmaI-si*d2, sigmaN, i);
    y4 = ILandauGaus(x, deltaI, xiI, sigmaI-si*dp, sigmaN, i);
    break;
  case 3: 
    y1 = ILandauGaus(x, deltaI, xiI, sigmaI, sigmaN+dp, i);
    y2 = ILandauGaus(x, deltaI, xiI, sigmaI, sigmaN+d2, i);
    y3 = ILandauGaus(x, deltaI, xiI, sigmaI, sigmaN-d2, i);
    y4 = ILandauGaus(x, deltaI, xiI, sigmaI, sigmaN-dp, i);
    break;
  default:
    return 0;
  } 
  
  Double_t d0  = y1 - y4;
  Double_t d1  = 2 * (y2 - y3);
  
  Double_t g   = 1/(2*dp) * (4*d1 - d0) / 3;
   
  return g;
}

//____________________________________________________________________
Double_t 
AliForwardUtil::NLandauGaus(Double_t x, Double_t delta, Double_t xi, 
			    Double_t sigma, Double_t sigmaN, Int_t n, 
			    const Double_t* a)
{
  // 
  // Evaluate 
  // @f[ 
  //   f_N(x;\Delta,\xi,\sigma') = \sum_{i=1}^N a_i f_i(x;\Delta,\xi,\sigma'a)
  // @f] 
  // 
  // where @f$ f(x;\Delta,\xi,\sigma')@f$ is the convolution of a
  // Landau with a Gaussian (see LandauGaus).  Note that 
  // @f$ a_1 = 1@f$, @f$\Delta_i = i(\Delta_1 + \xi\log(i))@f$, 
  // @f$\xi_i=i\xi_1@f$, and @f$\sigma_i'^2 = \sigma_n^2 + i\sigma_1^2@f$. 
  //  
  // References: 
  //  - <a href="http://dx.doi.org/10.1016/0168-583X(84)90472-5">Nucl.Instrum.Meth.B1:16</a>
  //  - <a href="http://dx.doi.org/10.1103/PhysRevA.28.615">Phys.Rev.A28:615</a>
  //  - <a href="http://root.cern.ch/root/htmldoc/tutorials/fit/langaus.C.html">ROOT implementation</a>
  // 
  // Parameters:
  //    x        Where to evaluate @f$ f_N@f$
  //    delta    @f$ \Delta_1@f$ 
  //    xi       @f$ \xi_1@f$
  //    sigma    @f$ \sigma_1@f$ 
  //    sigma_n  @f$ \sigma_n@f$ 
  //    n        @f$ N@f$ in the sum above.
  //    a        Array of size @f$ N-1@f$ of the weights @f$ a_i@f$ for 
  //                 @f$ i > 1@f$ 
  // 
  // Return:
  //    @f$ f_N(x;\Delta,\xi,\sigma')@f$ 
  //
  Double_t result = ILandauGaus(x, delta, xi, sigma, sigmaN, 1);
  for (Int_t i = 2; i <= n; i++) 
    result += a[i-2] * AliForwardUtil::ILandauGaus(x,delta,xi,sigma,sigmaN,i);
  return result;
}
namespace { 
  const Int_t kColors[] = { kRed+1, 
			    kPink+3, 
			    kMagenta+2, 
			    kViolet+2, 
			    kBlue+1, 
			    kAzure+3, 
			    kCyan+1, 
			    kTeal+2, 
			    kGreen+2, 
			    kSpring+3, 
			    kYellow+2, 
			    kOrange+2 };
}

//____________________________________________________________________
TF1*
AliForwardUtil::MakeNLandauGaus(Double_t  c, 
				Double_t  delta, Double_t xi, 
				Double_t  sigma, Double_t sigmaN, Int_t n, 
				const Double_t* a, 
				Double_t  xmin, Double_t xmax)
{
  // 
  // Generate a TF1 object of @f$ f_N@f$ 
  // 
  // Parameters:
  //    c         Constant			       
  //    delta     @f$ \Delta@f$ 		       
  //    xi 	      @f$ \xi_1@f$	       	       
  //    sigma     @f$ \sigma_1@f$ 	       	       
  //    sigma_n   @f$ \sigma_n@f$ 	       	       
  //    n 	      @f$ N@f$ - how many particles to sum to
  //    a         Array of size @f$ N-1@f$ of the weights @f$ a_i@f$ for 
  //                  @f$ i > 1@f$ 
  //    xmin      Least value of range  
  //    xmax      Largest value of range
  // 
  // Return:
  //    Newly allocated TF1 object
  //
  Int_t npar       = AliForwardUtil::ELossFitter::kN+n;
  TF1* landaun     = new TF1(Form("nlandau%d", n), &landauGausN,xmin,xmax,npar);
  // landaun->SetLineStyle(((n-2) % 10)+2); // start at dashed
  landaun->SetLineColor(kColors[((n-1) % 12)]); // start at red
  landaun->SetLineWidth(2);
  landaun->SetNpx(500);
  landaun->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}", "N");

  // Set the initial parameters from the seed fit 
  landaun->SetParameter(AliForwardUtil::ELossFitter::kC,      c);       
  landaun->SetParameter(AliForwardUtil::ELossFitter::kDelta,  delta);   
  landaun->SetParameter(AliForwardUtil::ELossFitter::kXi,     xi);      
  landaun->SetParameter(AliForwardUtil::ELossFitter::kSigma,  sigma);   
  landaun->SetParameter(AliForwardUtil::ELossFitter::kSigmaN, sigmaN); 
  landaun->FixParameter(AliForwardUtil::ELossFitter::kN,      n);       

  // Set the range and name of the scale parameters 
  for (UShort_t i = 2; i <= n; i++) {// Take parameters from last fit 
    landaun->SetParameter(AliForwardUtil::ELossFitter::kA+i-2, a[i-2]);
    landaun->SetParName(AliForwardUtil::ELossFitter::kA+i-2, Form("a_{%d}", i));
  }
  return landaun;
}
//____________________________________________________________________
TF1*
AliForwardUtil::MakeILandauGaus(Double_t  c, 
				Double_t  delta, Double_t xi, 
				Double_t  sigma, Double_t sigmaN, Int_t i, 
				Double_t  xmin, Double_t xmax)
{
  // 
  // Generate a TF1 object of @f$ f_I@f$ 
  // 
  // Parameters:
  //    c        Constant
  //    delta    @f$ \Delta@f$ 
  //    xi       @f$ \xi_1@f$	       
  //    sigma    @f$ \sigma_1@f$ 	       
  //    sigma_n  @f$ \sigma_n@f$ 	       
  //    i 	     @f$ i@f$ - the number of particles
  //    xmin     Least value of range
  //    xmax     Largest value of range
  // 
  // Return:
  //    Newly allocated TF1 object
  //
  Int_t npar       = AliForwardUtil::ELossFitter::kN+1;
  TF1* landaui     = new TF1(Form("ilandau%d", i), &landauGausI,xmin,xmax,npar);
  // landaui->SetLineStyle(((i-2) % 10)+2); // start at dashed
  landaui->SetLineColor(kColors[((i-1) % 12)]); // start at red
  landaui->SetLineWidth(1);
  landaui->SetNpx(500);
  landaui->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}", "i");

  // Set the initial parameters from the seed fit 
  landaui->SetParameter(AliForwardUtil::ELossFitter::kC,      c);       
  landaui->SetParameter(AliForwardUtil::ELossFitter::kDelta,  delta);   
  landaui->SetParameter(AliForwardUtil::ELossFitter::kXi,     xi);      
  landaui->SetParameter(AliForwardUtil::ELossFitter::kSigma,  sigma);   
  landaui->SetParameter(AliForwardUtil::ELossFitter::kSigmaN, sigmaN); 
  landaui->FixParameter(AliForwardUtil::ELossFitter::kN,      i);       

  return landaui;
}

//====================================================================
AliForwardUtil::ELossFitter::ELossFitter(Double_t lowCut, 
					 Double_t maxRange, 
					 UShort_t minusBins) 
  : fLowCut(lowCut), fMaxRange(maxRange), fMinusBins(minusBins), 
    fFitResults(0), fFunctions(0)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    lowCut     Lower cut of spectrum - data below this cuts is ignored
  //    maxRange   Maximum range to fit to 
  //    minusBins  The number of bins below maximum to use 
  //
  fFitResults.SetOwner();
  fFunctions.SetOwner();
}
//____________________________________________________________________
AliForwardUtil::ELossFitter::~ELossFitter()
{
  // 
  // Destructor
  // 
  //
  fFitResults.Delete();
  fFunctions.Delete();
}
//____________________________________________________________________
void
AliForwardUtil::ELossFitter::Clear()
{
  // 
  // Clear internal arrays 
  // 
  //
  fFitResults.Clear();
  fFunctions.Clear();
}
//____________________________________________________________________
TF1*
AliForwardUtil::ELossFitter::Fit1Particle(TH1* dist, Double_t sigman)
{
  // 
  // Fit a 1-particle signal to the passed energy loss distribution 
  // 
  // Note that this function clears the internal arrays first 
  // 
  // Parameters:
  //    dist    Data to fit the function to 
  //    sigman If larger than zero, the initial guess of the
  //               detector induced noise. If zero or less, then this 
  //               parameter is ignored in the fit (fixed at 0)
  // 
  // Return:
  //    The function fitted to the data 
  //

  // Clear the cache 
  Clear();
  
  // Find the fit range 
  dist->GetXaxis()->SetRangeUser(fLowCut, fMaxRange);
  
  // Get the bin with maximum 
  Int_t    peakBin = dist->GetMaximumBin();
  Double_t peakE   = dist->GetBinLowEdge(peakBin);
  
  // Get the low edge 
  dist->GetXaxis()->SetRangeUser(fLowCut, peakE);
  Int_t    minBin = peakBin - fMinusBins; // dist->GetMinimumBin();
  Double_t minE   = TMath::Max(dist->GetBinCenter(minBin),fLowCut);
  Double_t maxE   = dist->GetBinCenter(peakBin+2*fMinusBins);

  Int_t    minEb = dist->GetXaxis()->FindBin(minE);
  Int_t    maxEb = dist->GetXaxis()->FindBin(maxE);
  Double_t intg  = dist->Integral(minEb, maxEb);
  if (intg <= 0) {
    ::Warning("Fit1Particle", 
	      "Integral of %s between [%f,%f] [%03d,%03d] = %f < 0", 
	      dist->GetName(), minE, maxE, minEb, maxEb, intg);
    return 0;
  }
    
  // Restore the range 
  dist->GetXaxis()->SetRangeUser(0, fMaxRange);
  
  // Define the function to fit 
  TF1* landau1 = new TF1("landau1", landauGaus1, minE,maxE,kSigmaN+1);

  // Set initial guesses, parameter names, and limits  
  landau1->SetParameters(1,peakE,peakE/10,peakE/5,sigman);
  landau1->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}");
  landau1->SetNpx(500);
  landau1->SetParLimits(kDelta, minE, fMaxRange);
  landau1->SetParLimits(kXi,    0.00, fMaxRange);
  landau1->SetParLimits(kSigma, 1e-5, fMaxRange);
  if (sigman <= 0)  landau1->FixParameter(kSigmaN, 0);
  else              landau1->SetParLimits(kSigmaN, 0, fMaxRange);

  // Do the fit, getting the result object 
  ::Info("Fit1Particle", "Fitting in the range %f,%f", minE, maxE);
  TFitResultPtr r = dist->Fit(landau1, "RNQS", "", minE, maxE);
  // landau1->SetRange(minE, fMaxRange);
  fFitResults.AddAtAndExpand(new TFitResult(*r), 0);
  fFunctions.AddAtAndExpand(landau1, 0);

  return landau1;
}
//____________________________________________________________________
TF1*
AliForwardUtil::ELossFitter::FitNParticle(TH1* dist, UShort_t n, 
					  Double_t sigman)
{
  // 
  // Fit a N-particle signal to the passed energy loss distribution 
  //
  // If there's no 1-particle fit present, it does that first 
  // 
  // Parameters:
  //    dist   Data to fit the function to 
  //    n      Number of particle signals to fit 
  //    sigman If larger than zero, the initial guess of the
  //               detector induced noise. If zero or less, then this 
  //               parameter is ignored in the fit (fixed at 0)
  // 
  // Return:
  //    The function fitted to the data 
  //

  // Get the seed fit result 
  TFitResult* r = static_cast<TFitResult*>(fFitResults.At(0));
  TF1*        f = static_cast<TF1*>(fFunctions.At(0));
  if (!r || !f) { 
    f = Fit1Particle(dist, sigman);
    r = static_cast<TFitResult*>(fFitResults.At(0));
    if (!r || !f) { 
      ::Warning("FitNLandau", "No first shot at landau fit");
      return 0;
    }
  }

  // Get some parameters from seed fit 
  Double_t delta1  = r->Parameter(kDelta);
  Double_t xi1     = r->Parameter(kXi);
  Double_t maxEi   = n * (delta1 + xi1 * TMath::Log(n)) + 2 * n * xi1;
  Double_t minE    = f->GetXmin();

  Int_t    minEb = dist->GetXaxis()->FindBin(minE);
  Int_t    maxEb = dist->GetXaxis()->FindBin(maxEi);
  Double_t intg  = dist->Integral(minEb, maxEb);
  if (intg <= 0) {
    ::Warning("FitNParticle",
	      "Integral of %s between [%f,%f] [%03d,%03d] = %f < 0", 
	      dist->GetName(), minE, maxEi, minEb, maxEb, intg);
    return 0;
  }

  // Array of weights 
  TArrayD a(n-1);
  for (UShort_t i = 2; i <= n; i++) 
    a.fArray[i-2] = (n == 2 ? 0.05 : 0.000001);
  // Make the fit function 
  TF1* landaun = MakeNLandauGaus(r->Parameter(kC),
				 r->Parameter(kDelta),
				 r->Parameter(kXi),
				 r->Parameter(kSigma),
				 r->Parameter(kSigmaN),
				 n, a.fArray, minE, maxEi);
  landaun->SetParLimits(kDelta,  minE, fMaxRange);       // Delta
  landaun->SetParLimits(kXi,     0.00, fMaxRange);       // xi
  landaun->SetParLimits(kSigma,  1e-5, fMaxRange);       // sigma
  // Check if we're using the noise sigma 
  if (sigman <= 0)  landaun->FixParameter(kSigmaN, 0);
  else              landaun->SetParLimits(kSigmaN, 0, fMaxRange);

  // Set the range and name of the scale parameters 
  for (UShort_t i = 2; i <= n; i++) {// Take parameters from last fit 
    landaun->SetParLimits(kA+i-2, 0,1);
  }

  // Do the fit 
  ::Info("Fit1Particle", "Fitting in the range %f,%f", minE, maxEi);
  TFitResultPtr tr = dist->Fit(landaun, "RSQN", "", minE, maxEi);
  
  // landaun->SetRange(minE, fMaxRange);
  fFitResults.AddAtAndExpand(new TFitResult(*tr), n-1);
  fFunctions.AddAtAndExpand(landaun, n-1);
  
  return landaun;
}  
//____________________________________________________________________
TF1*
AliForwardUtil::ELossFitter::FitComposite(TH1* dist, Double_t sigman)
{
  // 
  // Fit a composite particle signal to the passed energy loss
  // distribution
  // 
  // Parameters:
  //    dist    Data to fit the function to 
  //    sigman If larger than zero, the initial guess of the
  //               detector induced noise. If zero or less, then this 
  //               parameter is ignored in the fit (fixed at 0)
  // 
  // Return:
  //    The function fitted to the data 
  //

  // Find the fit range 
  dist->GetXaxis()->SetRangeUser(fLowCut, fMaxRange);
  
  // Get the bin with maximum 
  Int_t    peakBin = dist->GetMaximumBin();
  Double_t peakE   = dist->GetBinLowEdge(peakBin);
  
  // Get the low edge 
  dist->GetXaxis()->SetRangeUser(fLowCut, peakE);
  Int_t    minBin = peakBin - fMinusBins; // dist->GetMinimumBin();
  Double_t minE   = TMath::Max(dist->GetBinCenter(minBin),fLowCut);
  Double_t maxE   = dist->GetBinCenter(peakBin+2*fMinusBins);

  // Get the range in bins and the integral of that range 
  Int_t    minEb = dist->GetXaxis()->FindBin(minE);
  Int_t    maxEb = dist->GetXaxis()->FindBin(maxE);
  Double_t intg  = dist->Integral(minEb, maxEb);
  if (intg <= 0) {
    ::Warning("Fit1Particle", 
	      "Integral of %s between [%f,%f] [%03d,%03d] = %f < 0", 
	      dist->GetName(), minE, maxE, minEb, maxEb, intg);
    return 0;
  }
    
  // Restore the range 
  dist->GetXaxis()->SetRangeUser(0, fMaxRange);
  
  // Define the function to fit 
  TF1* seed = new TF1("landauSeed", landauGaus1, minE,maxE,kSigmaN+1);

  // Set initial guesses, parameter names, and limits  
  seed->SetParameters(1,peakE,peakE/10,peakE/5,sigman);
  seed->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}");
  seed->SetNpx(500);
  seed->SetParLimits(kDelta, minE, fMaxRange);
  seed->SetParLimits(kXi,    0.00, fMaxRange);
  seed->SetParLimits(kSigma, 1e-5, fMaxRange);
  if (sigman <= 0)  seed->FixParameter(kSigmaN, 0);
  else              seed->SetParLimits(kSigmaN, 0, fMaxRange);

  // Do the fit, getting the result object 
  ::Info("FitComposite", "Fitting seed in the range %f,%f", minE, maxE);
  /* TFitResultPtr r = */ dist->Fit(seed, "RNQS", "", minE, maxE);

  maxE = dist->GetXaxis()->GetXmax();
  TF1* comp = new TF1("composite", landauGausComposite, 
		      minE, maxE, kSigma+1+4);
  comp->SetParNames("C",       "#Delta_{p}",       "#xi",       "#sigma",
		    "C#prime", "#Delta_{p}#prime", "#xi#prime", "#sigma#prim");
  comp->SetParameters(0.8 * seed->GetParameter(kC),  // 0 Primary weight 
		      seed->GetParameter(kDelta),    // 1 Primary Delta
		      seed->GetParameter(kDelta)/10, // 2 primary Xi
		      seed->GetParameter(kDelta)/5,  // 3 primary sigma
		      1.20 * seed->GetParameter(kC), // 5 Secondary weight
		      seed->GetParameter(kDelta),    // 6 secondary Delta
		      seed->GetParameter(kXi),       // 7 secondary Xi
		      seed->GetParameter(kSigma));   // 8 secondary sigma
		      
  // comp->SetParLimits(kC,       minE, fMaxRange); // C
  comp->SetParLimits(kDelta,      minE, fMaxRange); // Delta
  comp->SetParLimits(kXi,         0.00, fMaxRange); // Xi 
  comp->SetParLimits(kSigma,      1e-5, fMaxRange); // Sigma
  // comp->SetParLimits(kSigma+1, minE, fMaxRange); // C
  comp->SetParLimits(kSigma+2,    minE/10, fMaxRange); // Delta
  comp->SetParLimits(kSigma+3,    0.00,    fMaxRange); // Xi 
  comp->SetParLimits(kSigma+4,    1e-6,    fMaxRange); // Sigma
  comp->SetLineColor(kRed+1);
  comp->SetLineWidth(3);
  
  // Do the fit, getting the result object 
  ::Info("FitComposite", "Fitting composite in the range %f,%f", minE, maxE);
  /* TFitResultPtr r = */ dist->Fit(comp, "RNQS", "", minE, maxE);

#if 0
  TF1* part1 = static_cast<TF1*>(seed->Clone("part1"));
  part1->SetLineColor(kGreen+1);
  part1->SetLineWidth(4);
  part1->SetRange(minE, maxE);
  part1->SetParameters(comp->GetParameter(0), // C 
		       comp->GetParameter(1), // Delta
		       comp->GetParameter(2), // Xi
		       comp->GetParameter(3), // sigma
		       0);
  part1->Save(minE,maxE,0,0,0,0);
  dist->GetListOfFunctions()->Add(part1);

  TF1* part2 = static_cast<TF1*>(seed->Clone("part2"));
  part2->SetLineColor(kBlue+1);
  part2->SetLineWidth(4);
  part2->SetRange(minE, maxE);
  part2->SetParameters(comp->GetParameter(4), // C 
		       comp->GetParameter(5), // Delta
		       comp->GetParameter(6), // Xi
		       comp->GetParameter(7), // sigma
		       0);
  part2->Save(minE,maxE,0,0,0,0);
  dist->GetListOfFunctions()->Add(part2);
#endif
  return comp;
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
AliForwardUtil::Histos::Make(UShort_t d, Char_t r, 
			     const TAxis& etaAxis) const
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
  TH2D* hist = new TH2D(Form("FMD%d%c_cache", d, r), 
			Form("FMD%d%c cache", d, r),
			etaAxis.GetNbins(), etaAxis.GetXmin(), 
			etaAxis.GetXmax(), ns, 0, 2*TMath::Pi());
  hist->SetXTitle("#eta");
  hist->SetYTitle("#phi [radians]");
  hist->SetZTitle("d^{2}N_{ch}/d#etad#phi");
  hist->Sumw2();
  hist->SetDirectory(0);

  return hist;
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
AliForwardUtil::Histos::Clear(Option_t* option)
{
  // 
  // Clear data 
  // 
  // Parameters:
  //    option Not used 
  //
  if (fFMD1i) fFMD1i->Reset(option);
  if (fFMD2i) fFMD2i->Reset(option);
  if (fFMD2o) fFMD2o->Reset(option);
  if (fFMD3i) fFMD3i->Reset(option);
  if (fFMD3o) fFMD3o->Reset(option);
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



//
// EOF
//
