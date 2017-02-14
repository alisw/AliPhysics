#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TRandom.h"
#include "TSystem.h"
#include "AliGenPythia.h"
#include "AliLog.h"
#endif

AliGenPythia* GenPythia6_Powheg(){

  AliGenPythia *gener = new AliGenPythia(1);
  gener->SetProcess(VAR_PYTHIA_POWHEG_PROCESS);
  gener->SetStrucFunc(VAR_LHAPDF_STRUCFUNC_SET);
  gener->SetReadLHEF("pwgevents.lhe");
  gener->UseNewMultipleInteractionsScenario(); // pt ordering is better when coupling with POWHEG
  gener->SetProjectile(VAR_PROJECTILE_NAME,VAR_PROJECTILE_A,VAR_PROJECTILE_Z);
  gener->SetTarget(VAR_TARGET_NAME,VAR_TARGET_A,VAR_TARGET_Z);
  if ( VAR_PROJECTILE_A > 1 || VAR_TARGET_A > 1 ) {
    gener->SetNuclearPDF(VAR_NPDF_SET);
    gener->SetUseNuclearPDF(kTRUE);
  }
  if ( VAR_PROJECTILE_A != VAR_TARGET_A ) {
    gener->SetUseLorentzBoost(kTRUE);
  }
  if ( VAR_PYTHIA_POWHEG_PROCESS == kPyBeautyPWHG || VAR_PYTHIA_POWHEG_PROCESS == kPyCharmPWHG ) {
    gener->SetForceDecay(kSemiMuonic);
    gener->SwitchHFOff();
  }
  gener->SetPhiRange(0., 360.);
  gener->SetCutOnChild(1);
  gener->SetChildThetaRange(168.0,178.5);
  gener->SetChildPtRange(VAR_CHILD_PT_MIN, 1.e10);
  gener->SetNumberOfAcceptedParticles(1);
  gener->SetPdgCodeParticleforAcceptanceCut(13);
  gener->SetTrackingFlag(1);

  return gener;
}
