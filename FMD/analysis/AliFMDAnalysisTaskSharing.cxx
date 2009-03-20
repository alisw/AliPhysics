 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>
#include <TMath.h>
#include "AliFMDDebug.h"
#include "AliFMDAnalysisTaskSharing.h"
#include "AliAnalysisManager.h"
#include "AliESDFMD.h"
#include "AliFMDGeometry.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliFMDAnaParameters.h"
#include "AliFMDParameters.h"

ClassImp(AliFMDAnalysisTaskSharing)

//_____________________________________________________________________
AliFMDAnalysisTaskSharing::AliFMDAnalysisTaskSharing()
: fDebug(0),
  fESD(0x0),
  foutputESDFMD(),
  fEnergy(0),
  fNstrips(0),
  fSharedThis(kFALSE),
  fSharedPrev(kFALSE),
  fDiagList(),
  fStandalone(kTRUE),
  fEsdVertex(0),
  fStatus(kTRUE)
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
    fEnergy(0),
    fNstrips(0),
    fSharedThis(kFALSE),
    fSharedPrev(kFALSE),
    fDiagList(),
    fStandalone(kTRUE),
    fEsdVertex(0),
    fStatus(kTRUE)
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
				     1000,0,25);
      TH1F* hEdist_after  = new TH1F(Form("Edist_after_sharing_FMD%d%c", det, ringChar),
				     Form("Edist_after_sharing_FMD%d%c", det, ringChar),
				     1000,0,25);
      
      
      TH1F* hNstripsHit    = new TH1F(Form("N_strips_hit_FMD%d%c",det,ringChar),
				     Form("N_strips_hit_FMD%d%c",det,ringChar),
				     25,0,25);
      fDiagList.Add(hEdist);
      fDiagList.Add(hEdist_after);
      fDiagList.Add(hNstripsHit);

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
  
  Double_t vertex[3];
  GetVertex(vertex);
  fEsdVertex->SetXYZ(vertex);
  
  if(vertex[0] == 0 && vertex[1] == 0 && vertex[2] == 0) {
  
    fStatus = kFALSE;
    return;
  }
  else
    fStatus = kTRUE;
  const AliMultiplicity* testmult = fESD->GetMultiplicity();
  
  Int_t nTrackLets = testmult->GetNumberOfTracklets();
  
  if(nTrackLets < 1000) foutputESDFMD->SetUniqueID(kTRUE);
  else foutputESDFMD->SetUniqueID(kFALSE);
  
  AliESDFMD* fmd = fESD->GetFMDData();
  
  if (!fmd) return;
  Int_t nHits = 0;
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
	fEnergy = 0;
	fNstrips = 0;
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  foutputESDFMD->SetMultiplicity(det,ring,sec,strip,0.);
	  Float_t mult = fmd->Multiplicity(det,ring,sec,strip);
	  
	  if(mult == AliESDFMD::kInvalidMult || mult == 0) continue;
	  
	  Double_t eta  = fmd->Eta(det,ring,sec,strip);//EtaFromStrip(det,ring,sec,strip,vertex[2]);
	  //std::cout<<EtaFromStrip(det,ring,sec,strip,vertex[2]) <<"    "<<fmd->Eta(det,ring,sec,strip)<<std::endl;
	  
	  hEdist->Fill(mult);
	  if(fmd->IsAngleCorrected())
	    mult = mult/TMath::Cos(Eta2Theta(fmd->Eta(det,ring,sec,strip)));
	  Float_t Eprev = 0;
	  Float_t Enext = 0;
	  if(strip != 0)
	    if(fmd->Multiplicity(det,ring,sec,strip-1) != AliESDFMD::kInvalidMult) {
	      Eprev = fmd->Multiplicity(det,ring,sec,strip-1);
	      if(fmd->IsAngleCorrected())
		Eprev = Eprev/TMath::Cos(Eta2Theta(fmd->Eta(det,ring,sec,strip-1)));
	    }
	  if(strip != nstr - 1)
	    if(fmd->Multiplicity(det,ring,sec,strip+1) != AliESDFMD::kInvalidMult) {
	      Enext = fmd->Multiplicity(det,ring,sec,strip+1);
	      if(fmd->IsAngleCorrected())
		Enext = Enext/TMath::Cos(Eta2Theta(fmd->Eta(det,ring,sec,strip+1)));
	    }
	  
	  Float_t merged_energy = GetMultiplicityOfStrip(mult,eta,Eprev,Enext,det,ring,sec,strip);

	  if(merged_energy > 0 )
	    nHits++;
	  foutputESDFMD->SetMultiplicity(det,ring,sec,strip,merged_energy);
	  foutputESDFMD->SetEta(det,ring,sec,strip,eta);
	  
	}
      }
    }
  }
  
  if(fStandalone) {
    PostData(0, foutputESDFMD); 
    PostData(1, fEsdVertex); 
    PostData(2, fESD); 
    PostData(3, &fDiagList); 
  }
}
//_____________________________________________________________________
Float_t AliFMDAnalysisTaskSharing::GetMultiplicityOfStrip(Float_t mult,
							  Float_t eta,
							  Float_t Eprev,
							  Float_t Enext,
							  UShort_t   det,
							  Char_t  ring,
							  UShort_t /*sec*/,
							  UShort_t /*strip*/) {
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
 
  Float_t merged_energy = 0;
  Float_t nParticles = 0;
  Float_t cutLow  = 0.2;
  Float_t cutHigh = pars->GetMPV(det,ring,eta) -1*pars->GetSigma(det,ring,eta);
  // Float_t cutPart = pars->GetMPV(det,ring,eta) - 5*pars->GetSigma(det,ring,eta);
  Float_t Etotal  = mult;
  
  //if(mult > 5)
  //  std::cout<<mult<<"    "<<det<<"    "<<ring<<"   "<<sec<<"    "<<strip<<std::endl;
  
  if(foutputESDFMD->GetUniqueID() == kTRUE) {
    
    if(mult > cutLow ) {
      fEnergy = fEnergy + mult;
      fNstrips++;
    }
  if((Enext <0.01 && fEnergy >0) || fNstrips >2 ) {
    
          
    //if((fEnergy*TMath::Cos(Eta2Theta(eta))) > cutPart || fNstrips > 1) {
      nParticles = 1;
      merged_energy = fEnergy*TMath::Cos(Eta2Theta(eta));
      TH1F* hEdist = (TH1F*)fDiagList.FindObject(Form("Edist_after_sharing_FMD%d%c",det,ring));
      hEdist->Fill(fEnergy);
      TH1F* hNstrips = (TH1F*)fDiagList.FindObject(Form("N_strips_hit_FMD%d%c",det,ring));
      hNstrips->Fill(fNstrips);
      //  std::cout<<Form("Merged signals %f %f %f into %f , %f in strip %d, sec %d, ring %c, det %d",Eprev, mult, Enext, fEnergy/TMath::Cos(Eta2Theta(eta)),fEnergy,strip,sec,ring,det )<<std::endl;
      
      // }
    // else
    //std::cout<<Form("NO HIT  for  %f %f %f into %f , %f in strip %d, sec %d, ring %c, det %d, cuts %f , %f",Eprev, mult, Enext, fEnergy/TMath::Cos(Eta2Theta(eta)),fEnergy,strip,sec,ring,det,cutPart,cutHigh )<<std::endl;
    
    fEnergy  = 0;
    fNstrips = 0;
    return merged_energy;
  }
  
  return 0;
  
  }
  else {
     
  if(fSharedThis) {
    fSharedThis      = kFALSE;
    fSharedPrev      = kTRUE;
    return 0.;
  }
  
  if(mult < cutLow) {
    fSharedThis      = kFALSE;
    fSharedPrev      = kFALSE;
    return 0;
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
  
  Etotal = Etotal*TMath::Cos(Eta2Theta(eta));
  if(Etotal > 0) {
    merged_energy = Etotal;
    fSharedPrev      = kTRUE;
    //if(det == 3 && ring =='I')
    //  std::cout<<Form("Merged signals %f %f %f into %f , %f in strip %d, sec %d, ring %c, det %d",Eprev, mult, Enext, Etotal/TMath::Cos(Eta2Theta(eta)),Etotal,strip,sec,ring,det )<<std::endl;
  }
    else{// if(Etotal > 0) {
      //if(det == 3 && ring =='I')
      //	std::cout<<Form("NO HIT  for  %f %f %f into %f , %f in strip %d, sec %d, ring %c, det %d, cuts %f , %f",Eprev, mult, Enext, Etotal/TMath::Cos(Eta2Theta(eta)),Etotal,strip,sec,ring,det,cutPart,cutHigh )<<std::endl;
    fSharedThis      = kFALSE;
    fSharedPrev      = kFALSE;
  }
  // merged_energy = mult;
  
  return merged_energy; 
  }  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::GetVertex(Double_t* vertexXYZ) 
{
  const AliESDVertex* vertex = 0;
  vertex = fESD->GetPrimaryVertex();
  if(!vertex || (vertexXYZ[0] == 0 && vertexXYZ[1] == 0 && vertexXYZ[2] == 0))        
    vertex = fESD->GetPrimaryVertexSPD();
  if(!vertex || (vertexXYZ[0] == 0 && vertexXYZ[1] == 0 && vertexXYZ[2] == 0))        
    vertex = fESD->GetPrimaryVertexTPC();
  if(!vertex || (vertexXYZ[0] == 0 && vertexXYZ[1] == 0 && vertexXYZ[2] == 0))    
    vertex = fESD->GetVertex();
  if (vertex && (vertexXYZ[0] != 0 || vertexXYZ[1] != 0 || vertexXYZ[2] != 0)) {
    vertex->GetXYZ(vertexXYZ);
    //std::cout<<vertex->GetName()<<"   "<< vertex->GetTitle() <<"   "<< vertex->GetZv()<<std::endl;
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
Float_t AliFMDAnalysisTaskSharing::Eta2Theta(Float_t eta) {

  Float_t theta = 2*TMath::ATan(TMath::Exp(-1*eta));
  
  if(eta < 0)
    theta = theta-TMath::Pi();
  
  std::cout<<"From eta2Theta: "<<theta<<"   "<<eta<<std::endl;
  return theta;
  


}
//_____________________________________________________________________
Double_t AliFMDAnalysisTaskSharing::EtaFromStrip(UShort_t det, 
						Char_t ring, 
						UShort_t sector, 
						UShort_t strip, 
						Double_t zvtx)
{
   
  AliFMDGeometry* geo = AliFMDGeometry::Instance();
  
  Double_t x,y,z;
  geo->Detector2XYZ(det,ring,sector,strip,x,y,z);
  
  Double_t r = TMath::Sqrt(x*x+y*y);
  
  Double_t z_real      = z-zvtx;
  Double_t theta       = TMath::ATan2(r,z_real);
  // std::cout<<"From EtaFromStrip "<<theta<<std::endl;
  Double_t eta         =  -1*TMath::Log(TMath::Tan(0.5*theta));
  
  // std::cout<<det<<"   "<<ring<<"   "<<sector<<"   "<<strip<<"   "<<r<<"    "<<z_real<<"   "<<theta<<"    "<<eta<<std::endl;
  
  return eta;
}
//_____________________________________________________________________
//
// EOF
//
