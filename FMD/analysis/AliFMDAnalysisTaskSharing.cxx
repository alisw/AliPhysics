 
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
  foutputESDFMD(),
  fSharedThis(kFALSE),
  fSharedPrev(kFALSE),
  fDiagList(),
  fStandalone(kTRUE),
  fEsdVertex(0)
{
  // Default constructor
  DefineInput (0, AliESDEvent::Class());
  DefineOutput(0, AliESDFMD::Class());
  DefineOutput(1, AliESDVertex::Class());
  DefineOutput(2, AliESDEvent::Class());
  DefineOutput(3, TList::Class());
}
//_____________________________________________________________________
AliFMDAnalysisTaskSharing::AliFMDAnalysisTaskSharing(const char* name, Bool_t SE):
    AliAnalysisTask(name, "AnalysisTaskFMD"),
    fDebug(0),
    fESD(0x0),
    foutputESDFMD(),
    fSharedThis(kFALSE),
    fSharedPrev(kFALSE),
    fDiagList(),
    fStandalone(kTRUE),
    fEsdVertex(0)
{
  fStandalone = SE;
  if(fStandalone) {
    DefineInput (0, AliESDEvent::Class());
    DefineOutput(0, AliESDFMD::Class());
    DefineOutput(1, AliESDVertex::Class());
    DefineOutput(2, AliESDEvent::Class());
    DefineOutput(3, TList::Class());
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::CreateOutputObjects()
{
  if(!foutputESDFMD)
    foutputESDFMD = new AliESDFMD();
  
  if(!fEsdVertex)
    fEsdVertex    = new AliESDVertex();
  //Diagnostics
  fDiagList.SetName("Sharing diagnostics");
  for(Int_t det = 1; det<=3; det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    
    for(Int_t iring = 0;iring<nRings; iring++) {
      Char_t ringChar = (iring == 0 ? 'I' : 'O');
      TH1F* hEdist        = new TH1F(Form("Edist_before_sharing_FMD%d%c", det, ringChar),
				     Form("Edist_before_sharing_FMD%d%c", det, ringChar),
				     200,0,5);
      TH1F* hEdist_after  = new TH1F(Form("Edist_after_sharing_FMD%d%c", det, ringChar),
				     Form("Edist_after_sharing_FMD%d%c", det, ringChar),
				     200,0,5);
      fDiagList.Add(hEdist);
      fDiagList.Add(hEdist_after);

    }
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::ConnectInputData(Option_t */*option*/)
{
  if(fStandalone)
    fESD = (AliESDEvent*)GetInputData(0);
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::Exec(Option_t */*option*/)
{
  
  AliESD* old = fESD->GetAliESDOld();
  if (old) {
    fESD->CopyFromOldESD();
  }
  
  foutputESDFMD->Clear();
  
  AliESDFMD* fmd = fESD->GetFMDData();
  
  if (!fmd) return;
  
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      
      TH1F* hEdist = (TH1F*)fDiagList.FindObject(Form("Edist_before_sharing_FMD%d%c",det,ring));
      
      for(UShort_t sec =0; sec < nsec;  sec++) {
	fSharedThis      = kFALSE;
	fSharedPrev      = kFALSE;
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  foutputESDFMD->SetMultiplicity(det,ring,sec,strip,0.);
	  Float_t mult = fmd->Multiplicity(det,ring,sec,strip);

	  if(mult == AliESDFMD::kInvalidMult || mult == 0) continue;
	  
	  hEdist->Fill(mult);
	  Float_t Eprev = 0;
	  Float_t Enext = 0;
	  if(strip != 0)
	    if(fmd->Multiplicity(det,ring,sec,strip-1) != AliESDFMD::kInvalidMult)
	      Eprev = fmd->Multiplicity(det,ring,sec,strip-1);
	  if(strip != nstr - 1)
	    if(fmd->Multiplicity(det,ring,sec,strip+1) != AliESDFMD::kInvalidMult)
	    Enext = fmd->Multiplicity(det,ring,sec,strip+1);
	  
	  Float_t nParticles = GetMultiplicityOfStrip(mult,Eprev,Enext,det,ring);
	  foutputESDFMD->SetMultiplicity(det,ring,sec,strip,nParticles);
	  foutputESDFMD->SetEta(det,ring,sec,strip,fmd->Eta(det,ring,sec,strip));
	  
	}
      }
    }
  }
  
  Double_t vertex[3];
  GetVertex(vertex);
  fEsdVertex->SetXYZ(vertex);
  if(fStandalone) {
    PostData(0, foutputESDFMD); 
    PostData(1, fEsdVertex); 
    PostData(2, fESD); 
    PostData(3, &fDiagList); 
  }
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
  TH1F* hEdist = (TH1F*)fDiagList.FindObject(Form("Edist_after_sharing_FMD%d%c",det,ring));
  hEdist->Fill(Etotal);
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
void AliFMDAnalysisTaskSharing::GetVertex(Double_t* vertexXYZ) 
{
  const AliESDVertex* vertex = 0;
  vertex = fESD->GetPrimaryVertex();
  if (!vertex)        vertex = fESD->GetPrimaryVertexSPD();
  if (!vertex)        vertex = fESD->GetPrimaryVertexTPC();
  if (!vertex)        vertex = fESD->GetVertex();
  
  if (vertex) {
    vertex->GetXYZ(vertexXYZ);
    return;
  }
  else if (fESD->GetESDTZERO()) { 
    vertexXYZ[0] = 0;
    vertexXYZ[1] = 0;
    vertexXYZ[2] = fESD->GetT0zVertex();

    return;
  }
  
  return;
}
//_____________________________________________________________________
//
// EOF
//
