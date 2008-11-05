 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>

#include "AliFMDAnalysisTaskSharing.h"
#include "AliAnalysisManager.h"
#include "AliESDFMD.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliESDVertex.h"
#include "AliFMDAnaParameters.h"

ClassImp(AliFMDAnalysisTaskSharing)

//_____________________________________________________________________
AliFMDAnalysisTaskSharing::AliFMDAnalysisTaskSharing()
: fDebug(0),
  fESD(0x0),
  fOutputESD(),
  foutputESDFMD(),
  fSharedThis(kFALSE),
  fSharedPrev(kFALSE)
{
  // Default constructor
  DefineInput (0, AliESDEvent::Class());
  DefineOutput(0, AliESDEvent::Class());
}
//_____________________________________________________________________
AliFMDAnalysisTaskSharing::AliFMDAnalysisTaskSharing(const char* name):
    AliAnalysisTask(name, "AnalysisTaskFMD"),
    fDebug(0),
    fESD(0x0),
    fOutputESD(),
    foutputESDFMD(),
    fSharedThis(kFALSE),
    fSharedPrev(kFALSE)

{
  DefineInput (0, AliESDEvent::Class());
  DefineOutput(0, AliESDEvent::Class());
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::CreateOutputObjects()
{
  fOutputESD.CreateStdContent();
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::ConnectInputData(Option_t */*option*/)
{
  fESD = (AliESDEvent*)GetInputData(0);
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::Exec(Option_t */*option*/)
{
  AliESD* old = fESD->GetAliESDOld();
  if (old) {
    fESD->CopyFromOldESD();
  }
  
  foutputESDFMD.Clear();
  
  fOutputESD.SetPrimaryVertexSPD(fESD->GetPrimaryVertexSPD());
  
  AliESDFMD* fmd = fESD->GetFMDData();
  
  if (!fmd) return;
  
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      for(UShort_t sec =0; sec < nsec;  sec++) {
	fSharedThis      = kFALSE;
	fSharedPrev      = kFALSE;
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  foutputESDFMD.SetMultiplicity(det,ring,sec,strip,0.);
	  Float_t mult = fmd->Multiplicity(det,ring,sec,strip);
	  if(mult == AliESDFMD::kInvalidMult || mult == 0) continue;
	  	  
	  Float_t Eprev = 0;
	  Float_t Enext = 0;
	  if(strip != 0)
	    if(fmd->Multiplicity(det,ring,sec,strip-1) != AliESDFMD::kInvalidMult)
	      Eprev = fmd->Multiplicity(det,ring,sec,strip-1);
	  if(strip != nstr - 1)
	    if(fmd->Multiplicity(det,ring,sec,strip+1) != AliESDFMD::kInvalidMult)
	    Enext = fmd->Multiplicity(det,ring,sec,strip+1);
	  
	  Float_t nParticles = GetMultiplicityOfStrip(mult,Eprev,Enext,det,ring);
	  foutputESDFMD.SetMultiplicity(det,ring,sec,strip,nParticles);
	  foutputESDFMD.SetEta(det,ring,sec,strip,fmd->Eta(det,ring,sec,strip));
	  Float_t eta = fmd->Eta(det,ring,sec,strip);
	  
	}
      }
    }
  }
  fOutputESD.SetFMDData(&foutputESDFMD);
    
  PostData(0, &fOutputESD); 
  
}
//_____________________________________________________________________
Float_t AliFMDAnalysisTaskSharing::GetMultiplicityOfStrip(Float_t mult,
							  Float_t Eprev,
							  Float_t Enext,
							  Int_t   det,
							  Char_t  ring) {
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  Float_t nParticles = 0;
  Float_t cutLow  = 0.2;
  Float_t cutHigh = pars->GetMPV(det,ring) - 2*pars->GetSigma(det,ring);
  Float_t Etotal = mult;
  /* 
  if(mult > 3*pars->GetMPV(det,ring) && 
     (Enext > 3*pars->GetMPV(det,ring) || (Enext > 3*pars->GetMPV(det,ring))))
    return 0;
  
  if(mult > 5*pars->GetMPV(det,ring))
    return 0;
  */
  if(fSharedThis) {
    fSharedThis      = kFALSE;
    fSharedPrev      = kTRUE;
    return 0.;
  }
  
  if(Etotal < 0.33*pars->GetMPV(det,ring)) {
    fSharedThis      = kFALSE;
    fSharedPrev      = kFALSE;
    return 0.; 
  }
  
  if(Eprev > cutLow && Eprev < cutHigh && !fSharedPrev ) {
    Etotal += Eprev;
  }
  
  if(Enext > cutLow && Enext < cutHigh ) {
    Etotal += Enext;
    fSharedThis      = kTRUE;
  }
  
  if(Etotal > cutHigh ) {
    nParticles = 1;
    fSharedPrev      = kTRUE;
  }
  else {
    fSharedThis      = kFALSE;
    fSharedPrev      = kFALSE;
  }
  
  return nParticles;
}

//_____________________________________________________________________
//
// EOF
//
