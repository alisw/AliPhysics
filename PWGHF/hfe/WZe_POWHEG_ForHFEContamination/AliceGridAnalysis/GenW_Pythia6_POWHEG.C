#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TRandom.h"
#include "TSystem.h"
#include "AliGenPythia.h"
#include "AliLog.h"
#endif

AliGenPythia* GenW_Pythia6_POWHEG(){

  gSystem->Setenv("LHAPATH",gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets")); // Needed to run lhapdf-5.9.1 on grid
  AliLog::SetClassDebugLevel("AliGenPythia",1);

  AliGenPythia *gener = new AliGenPythia(1);
  gener->SetProcess(kPyWPWHG);
  gener->SetStrucFunc(kCT10nlo);
  gener->SetReadLHEF("pwgevents.lhe");
  gener->UseNewMultipleInteractionsScenario(); // pt ordering is better when coupling with POWHEG
  gener->SetProjectile("p",1,1); // pp
  gener->SetTarget("p",1,1);
  gener->SetNuclearPDF(19);    // 0: ESK08, 8: EPS08, 9: EPS09lo, 19: EPS09nlo
  gener->SetUseNuclearPDF(kTRUE);
  gener->SetUseLorentzBoost(kTRUE);
  gener->SetPhiRange(0., 360.);
  gener->SetCutOnChild(1);
  //gener->SetChildThetaRange(168.0,178.5);
  gener->SetChildYRange(-1,1);
  gener->SetNumberOfAcceptedParticles(1);
  gener->SetPdgCodeParticleforAcceptanceCut(11);
  gener->SetTrackingFlag(1);

  return gener;
}
