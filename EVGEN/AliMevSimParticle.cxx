
#include "AliMevSimParticle.h"


ClassImp(AliMevSimParticle)

   
///////////////////////////////////////////////////////////////////////////////////////

AliMevSimParticle::AliMevSimParticle()
  : TMevSimPartTypeParams() {

}

///////////////////////////////////////////////////////////////////////////////////////

AliMevSimParticle::AliMevSimParticle(PDG_t pdg, Int_t multmean, Int_t multvc, 
		  Float_t tempmean, Float_t tempstdev, Float_t sigmamean,
		  Float_t sigmastdev, Float_t expvelmean, Float_t expvelstdev)

  : TMevSimPartTypeParams(0, multmean, multvc, tempmean, tempstdev, 
			  sigmamean, sigmastdev, expvelmean, expvelstdev)  {


  // Calculate geant ID from pdg
  fConv = new TMevSimConverter();
  fPdg = pdg;
  if (fConv) fGPid = fConv->IdFromPDG(pdg);  

}

///////////////////////////////////////////////////////////////////////////////////////

AliMevSimParticle::~AliMevSimParticle() {
}

///////////////////////////////////////////////////////////////////////////////////////

void  AliMevSimParticle::SetPDG(PDG_t pdg) {

  fPdg = pdg;
  fGPid = fConv->IdFromPDG(pdg);
}

///////////////////////////////////////////////////////////////////////////////////////

PDG_t AliMevSimParticle::GetPDG() {
  
  fPdg = (PDG_t)fConv->PDGFromId(fGPid);
  return fPdg;
}

///////////////////////////////////////////////////////////////////////////////////////


  

