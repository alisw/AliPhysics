/*************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Class for electrons from beauty study
// Counting electrons from beauty
// by DCA cuts, background subtraction 
//
// Authors:
//   Hongyan Yang <hongyan@physi.uni-heidelberg.de>
//

#include "TMath.h"
#include "TList.h"
#include "AliLog.h"

#include <TParticle.h>
#include <TDatabasePDG.h>

#include "THnSparse.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliMCParticle.h"
#include "AliESDtrack.h"

#include "AliStack.h" 

#include "AliHFEdisplacedElectrons.h"

ClassImp(AliHFEdisplacedElectrons)

//__________________________________________________________
const Float_t AliHFEdisplacedElectrons::fgkDcaMinPtIntv[13] = {
  // 
  // define DCA low limits for single electrons cut in each Pt bin 
  // preliminary numbers 
  // these numbers should be replaced by the best numbers determined by 
  // cut efficiency study: signal/background  (not used right now)
  //
  450, // 0.0 - 0.5
  450, // 0.5 - 1.0
  450, // 1.0 - 1.5
  500, // 1.5 - 2.0
  400, // 2.0 - 2.5
  300, // 2.5 - 3.0
  300, // 3.0 - 4.0
  300, // 4.0 - 5.0
  200, // 5.0 - 7.0
  150, // 7.0 - 9.0
  150, // 9.0 - 12.0
  100, //12.0 - 16.0
  50}; //16.0 - 20.0
//__________________________________________________________
const Float_t AliHFEdisplacedElectrons::fgkPtIntv[14] = {
  //
  // define pT bins for spectra of single electrons 
  //
  0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.0,9.0,12.0,16.0,20.0};

//__________________________________________________________
const Char_t *AliHFEdisplacedElectrons::fgkKineVar[3] = {
  "y", "phi", "pt"};

//__________________________________________________________
const Char_t *AliHFEdisplacedElectrons::fgkKineVarTitle[3] ={
  "rapidity;y;dN/dy",   "azimuthal;#phi;dN/d#phi",   "transverse momentum;p_{T};dN/dp_{T}", 
};

//__________________________________________________________
AliHFEdisplacedElectrons::AliHFEdisplacedElectrons():
  fDebugLevel(0)
  , fTHnSparseDcaMcPionInfo(NULL)
  , fTHnSparseDcaMcEleInfo(NULL)
  , fTHnSparseDcaDataEleInfo(NULL)
  , fOutputList(0x0)
{
  //
  // default constructor
  //

}

//__________________________________________________________
AliHFEdisplacedElectrons::AliHFEdisplacedElectrons(const AliHFEdisplacedElectrons &ref):
  TObject(ref)
  , fDebugLevel(ref.fDebugLevel)
  , fTHnSparseDcaMcPionInfo(ref.fTHnSparseDcaMcPionInfo)
  , fTHnSparseDcaMcEleInfo(ref.fTHnSparseDcaMcEleInfo)
  , fTHnSparseDcaDataEleInfo(ref.fTHnSparseDcaDataEleInfo)
  , fOutputList(ref.fOutputList)
  
{
  //
  // copy constructor
  //
}


//__________________________________________________________
AliHFEdisplacedElectrons&AliHFEdisplacedElectrons::operator=(const AliHFEdisplacedElectrons &ref)
{
  //
  // Assignment operator
  //

  if(this == &ref) return *this;
  AliHFEdisplacedElectrons::operator=(ref);

  fDebugLevel = ref.fDebugLevel;

  fTHnSparseDcaMcPionInfo = ref.fTHnSparseDcaMcPionInfo;
  fTHnSparseDcaMcEleInfo = ref.fTHnSparseDcaMcEleInfo;
  fTHnSparseDcaDataEleInfo = ref.fTHnSparseDcaDataEleInfo;
  fOutputList = ref.fOutputList;

  return *this;
}

//__________________________________________________________
AliHFEdisplacedElectrons::~AliHFEdisplacedElectrons()
{
  //
  // default constructor
  //

  if(fTHnSparseDcaMcPionInfo) 
    delete fTHnSparseDcaMcPionInfo;
  if(fTHnSparseDcaMcEleInfo) 
    delete fTHnSparseDcaMcEleInfo;
  if(fTHnSparseDcaDataEleInfo) 
    delete fTHnSparseDcaDataEleInfo;
    
  if(fOutputList){
    fOutputList->Clear();
    delete fOutputList;
    
  }

  Printf("analysis done\n");
}

//__________________________________________________________
void AliHFEdisplacedElectrons::InitAnalysis(){
  //
  // init analysis (no intialization yet)
  //
  

  Printf("initialize analysis\n");

}


//__________________________________________________________
void AliHFEdisplacedElectrons::CreateOutputs(TList* const displacedList){

  //
  //  create output fOutputList
  //

  // THnSparseF
  // 8 interested electron sources: others, photon conv, direct photon,  pi0, eta, b, b->c, c
  // 21 possible DCA cuts XY
  // 21 possible DCA cuts Z
  // 13 pT bins
  // 10 rapidity bins
  // 10 azimuthal angle phi bins


  if(!displacedList) return;

  fOutputList = displacedList;
  fOutputList -> SetName("information");

 
  // electron source 
  Int_t nBinsEleSource = 8;
  Double_t minEleSource = -1.5; 
  Double_t maxEleSource = 6.5;  
  Double_t *binLimEleSource = new Double_t[nBinsEleSource+1];
  for(Int_t i=0; i<=nBinsEleSource; i++)
    binLimEleSource[i] = minEleSource + i*(maxEleSource-minEleSource)/nBinsEleSource;

  // dca bins: the same as XY and Z
  // these should be variable bins as well
  Int_t nBinsDca = kNDcaMin-1;  // 12 bins
  Double_t minDca = -0.5; 
  Double_t maxDca = 20.5; 
  Double_t dcaBinWidth = (maxDca-minDca)/nBinsDca;
  Double_t *binLimDca = new Double_t[nBinsDca+1];
  for(Int_t i=0; i<=nBinsDca; i++)
    binLimDca[i] = minDca + i*dcaBinWidth;
  

 // pt bins
  Int_t nBinsPt = kNPtIntv-1;
  Double_t *binLimPt = new Double_t[nBinsPt+1];
  for(Int_t i=0; i<=nBinsPt; i++)
    binLimPt[i] = fgkPtIntv[i];  // variable bins

  // rapidity bins
  Int_t nBinsRap = 10;
  Double_t minRap = -1.0; 
  Double_t maxRap = 1.0;
  Double_t *binLimRap = new Double_t[nBinsRap+1];
  for(Int_t i=0; i<=nBinsRap; i++)
    binLimRap[i] = minRap + i*(maxRap-minRap)/nBinsRap;

  // azumuthal phi angle
  Int_t nBinsPhi = 10;
  Double_t minPhi = 0; 
  Double_t maxPhi = 2*TMath::Pi();
  Double_t *binLimPhi = new Double_t[nBinsPhi+1];
  for(Int_t i=0; i<=nBinsPhi; i++)
    binLimPhi[i] = minPhi + i*(maxPhi-minPhi)/nBinsPhi;



  //for MC pionss
  const Int_t nVarPion = 5;
  Int_t iBinPion[nVarPion] = {nBinsDca, nBinsDca, nBinsPt, nBinsRap, nBinsPhi};
   
  //  THnSparseF *fTHnSparseDcaMcPionInfo = NULL; // empty for the moment 
  if(HasMCData()){
    fTHnSparseDcaMcPionInfo = new THnSparseF("dcaMcPionInfo", 
					     "MC info:;dcaXY [50 #mum];dcaZ [50 #mum];pT [GeV/c];y [rapidity];#phi [rad];", 	
					     nVarPion, iBinPion);
    fTHnSparseDcaMcPionInfo->SetBinEdges(0, binLimDca);  // dca xy cut
    fTHnSparseDcaMcPionInfo->SetBinEdges(1, binLimDca);  // dca z cut
    fTHnSparseDcaMcPionInfo->SetBinEdges(2, binLimPt);   // pt
    fTHnSparseDcaMcPionInfo->SetBinEdges(3, binLimRap);  // rapidity
    fTHnSparseDcaMcPionInfo->SetBinEdges(4, binLimPhi);  // phi
    fTHnSparseDcaMcPionInfo->Sumw2();

    fOutputList -> AddAt(fTHnSparseDcaMcPionInfo, kMCpion);
  }

  // for MC electrons
  const Int_t nVar = 6;
  Int_t iBin[nVar] = {nBinsEleSource,nBinsDca, nBinsDca, nBinsPt, nBinsRap, nBinsPhi};
  
  //  THnSparseF *fTHnSparseDcaMcEleInfo = NULL; // empty for the moment 
  if(HasMCData()){
    fTHnSparseDcaMcEleInfo = new THnSparseF("dcaMcElectronInfo", 
					 "MC info:;ID [electron source id];dcaXY [50 #mum];dcaZ [50 #mum];pT [GeV/c];y [rapidity];#phi [rad];",		
					 nVar, iBin);
    fTHnSparseDcaMcEleInfo->SetBinEdges(0, binLimEleSource); // electron source
    fTHnSparseDcaMcEleInfo->SetBinEdges(1, binLimDca);  // dca xy cut
    fTHnSparseDcaMcEleInfo->SetBinEdges(2, binLimDca);  // dca z cut
    fTHnSparseDcaMcEleInfo->SetBinEdges(3, binLimPt);   // pt
    fTHnSparseDcaMcEleInfo->SetBinEdges(4, binLimRap);  // rapidity
    fTHnSparseDcaMcEleInfo->SetBinEdges(5, binLimPhi);  // phi
    fTHnSparseDcaMcEleInfo->Sumw2();
    fOutputList -> AddAt(fTHnSparseDcaMcEleInfo, kMCelectron);
  }

  
  // for ESD: HFE pid

  //  THnSparseF *fTHnSparseDcaDataEleInfo = NULL;  // empty for the moment
  const Int_t nVarData = 5;
  Int_t iBinData[nVarData] = {nBinsDca, nBinsDca, nBinsPt, nBinsRap, nBinsPhi};  
  
  fTHnSparseDcaDataEleInfo = new THnSparseF("dcaDataElectronInfo", 
					 "Data info:;dcaXY [50 #mum];dcaZ [50 #mum];pT [GeV/c];y [rapidity];#phi [rad];",		
					 nVarData, iBinData);    
  fTHnSparseDcaDataEleInfo->SetBinEdges(0, binLimDca);  // dca xy cut
  fTHnSparseDcaDataEleInfo->SetBinEdges(1, binLimDca);  // dca z cut
  fTHnSparseDcaDataEleInfo->SetBinEdges(2, binLimPt);   // pt
  fTHnSparseDcaDataEleInfo->SetBinEdges(3, binLimRap);  // rapidity
  fTHnSparseDcaDataEleInfo->SetBinEdges(4, binLimPhi);  // phi
  fTHnSparseDcaDataEleInfo->Sumw2();
  
  fOutputList -> AddAt(fTHnSparseDcaDataEleInfo, kData);
  
  AliInfo("THnSparse histograms are created\n");
  fOutputList->Print();

}


//__________________________________________________________
void AliHFEdisplacedElectrons::FillMCOutput(AliESDEvent * const fESDEvent, AliESDtrack* const esdTrack, AliStack * const stack)
{

  // fill output
  if(!esdTrack) return;
  AliESDtrack *track = esdTrack;
  Double_t pt   = track->Pt();
  Double_t eta = track->Eta();
  Double_t phi = track->Phi();
  
  Int_t label = track->GetLabel();
  if(label<0 || label>stack->GetNtrack()) return;  

  TParticle *particle = stack->Particle(label);
  if(!particle) return;
 

  // obtain impact parameters in xy and y
  const AliESDVertex *primVtx = fESDEvent->GetPrimaryVertex();      

  Float_t magneticField = 5;  // initialized as 5kG
  magneticField = fESDEvent->GetMagneticField();  // in kG 
  

  Double_t dz[2];   // error of dca in cm
  Double_t covardz[3];
  track->PropagateToDCA(primVtx,magneticField, 1000., dz, covardz); 


  Double_t dcaXY = TMath::Abs(dz[0])*1.0e4;  // conv dca in cm to dca in micron 
  Double_t dcaZ = TMath::Abs(dz[1])*1.0e4;

  // do PID with MC

  if(HasMCData() && TMath::Abs(GetMCpid(stack, label))==kPDGelectron){
  
    Int_t eleLabel = label; 
    
    const Int_t nvarMC=6;
    Double_t var[nvarMC];
    var[0] = -1;
    Int_t sourcePdg = ElectronFromSource(stack, eleLabel);     

    if(sourcePdg==kPDGgamma){ // check direct photon or not  // fixme      
      if(ElePhotonDirect(stack, eleLabel)!=-1)	
	var[0] = kEleDirectPhotonConv;    
      else 	
	var[0] = kElePhotonConv;     
      if(fDebugLevel>=10) 
	printf("photonconv=====this electron %d is from %d, source id=%.1f\n", eleLabel, sourcePdg, var[0]);   
    }   // photon or direct photon -> e
    
    if(sourcePdg==kPDGpi0){      
      var[0] = kElePi0;   
      if(fDebugLevel>=10) printf("pi0======this electron %d is from %d, source id=%.1f\n", eleLabel, sourcePdg, var[0]);    
    } 

    if(sourcePdg==kPDGeta){ 
      var[0] = kEleEta;     
      if(fDebugLevel>=10) 
	printf("eta=====this electron %d is from %d, source id=%.1f\n", eleLabel, sourcePdg, var[0]);    
    }   // from eta -> e

    if(TMath::Abs(sourcePdg%10000)/100==kPDGbeauty ||    // for special intermediate meson states: like 100553       
       TMath::Abs(sourcePdg)/1000==kPDGbeauty ||      
       TMath::Abs(sourcePdg)/100==kPDGbeauty ||     
       TMath::Abs(sourcePdg)==kPDGbeauty){     
      var[0]=kEleB;       
      if(fDebugLevel>=10) 
	printf("beauty======electron %d is from %d with id %.1f\n", eleLabel, ElectronFromSource(stack, eleLabel), var[0]);   
    } // direct beauty  -> e      
    if(TMath::Abs(sourcePdg%10000)/100==kPDGcharm ||    // for special intermediate meson states: like 1004**      
       TMath::Abs(sourcePdg)/1000==kPDGcharm ||       
       TMath::Abs(sourcePdg)/100==kPDGcharm ||        
       TMath::Abs(sourcePdg)==kPDGcharm){           
      // two cases: 
      //            electron from b->c->e     
      //            electron from c->e     
      if(ElectronFromCharm(stack, eleLabel)!=-1){
	var[0] = ElectronFromCharm(stack, eleLabel);
	if(fDebugLevel>=10)
	  printf("charm----->electron %d has mother %d is from %.1f\n", eleLabel, ElectronFromSource(stack, eleLabel), var[0]);
      } 
    }  // charm electrons: b->c->e or c->e

    if(fDebugLevel>=10) printf("others----->electron %d has mother %d is from %.1f\n", eleLabel, ElectronFromSource(stack, eleLabel), var[0]);
    
    if(dcaXY<1000) var[1] = Int_t(dcaXY)/50;   // larger than 1mm should go to the last bin
    else
      var[1] = 20;
    if(dcaZ<1000) var[2] = Int_t(dcaZ)/50;
    else
      var[2] = 20;

    var[3] = pt; // pt 
    var[4] = eta; // eta
    var[5] = phi; // phi

    (dynamic_cast<THnSparseF *>(fOutputList->At(kMCelectron)))->Fill(var);
  }

  if(HasMCData() && TMath::Abs(GetMCpid(stack, label))==kPDGpion){
      
    const Int_t nvarPionMC=5;
    Double_t varPion[nvarPionMC];
    if(dcaXY<1000)
      varPion[0] = Int_t(dcaXY)/50; // dca xy 
    else
      varPion[0] = 20;
    if(dcaZ<1000)
      varPion[1] = Int_t(dcaZ)/50; // dca Z 
    else
      varPion[1] = 20;
    
//     varPion[0] = TMath::Abs(dcaXY); // dca xy 
//     varPion[1] = TMath::Abs(dcaZ); // dca z    
    varPion[2] = pt; // pt 
    varPion[3] = eta; //eta
    varPion[4] = phi; // phi
    
    (dynamic_cast<THnSparseF *>(fOutputList->At(kMCpion)))->Fill(varPion);
  }
  
}



//__________________________________________________________
void AliHFEdisplacedElectrons::FillESDOutput(AliESDEvent * const fESDEvent, AliESDtrack* const esdTrack)
{

  // fill output
  if(!esdTrack) return;
  AliESDtrack *track = esdTrack;

  Double_t pt   = track->Pt();
  Double_t eta = track->Eta();
  Double_t phi = track->Phi();
  
  // obtain impact parameters in xy and y
  const AliESDVertex *primVtx = fESDEvent->GetPrimaryVertex();    
  
  Float_t magneticField = 5;  // initialized as 5kG
  magneticField = fESDEvent->GetMagneticField();  // in kG 
  Double_t dz[2];   // error of dca in cm
  Double_t covardz[3];
  track->PropagateToDCA(primVtx,magneticField, 1000., dz, covardz); 
  Double_t dcaXY = TMath::Abs(dz[0]*1.0e4);  // conv dca in cm to dca in micron 
  Double_t dcaZ = TMath::Abs(dz[1]*1.0e4);
 
  
  const Int_t nvarData = 5;
  Double_t varData[nvarData]; 
  // fixme
  if(dcaXY<1000)
    varData[0] = Int_t(dcaXY)/50; // dca xy 
  else 
    varData[0] = 20;
  if(dcaZ<1000)
    varData[1] = Int_t(dcaZ)/50; // dca Z 
  else
    varData[1] = 20;

  varData[2] = pt; // pt 
  varData[3] = eta; //eta
  varData[4] = phi; // phi
  
  (dynamic_cast<THnSparseF *>(fOutputList->At(kData)))->Fill(varData);

}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::ElectronFromSource(AliStack * const stack, Int_t label) const
{

  //
  //  return electron source label via electron label
  //

  // kEleSource is supposed to be either photon conv, direct photon conv, pi0, eta, beauty --> 0, 1, 2, 3, 4

  Int_t eleLabel = label;

  if(eleLabel<0 || eleLabel>stack->GetNtrack()) return -1;

  TParticle *particle = stack->Particle(eleLabel);
  if(!particle) return -1;

  Int_t motherLabel = particle->GetFirstMother();
  if(motherLabel<0 || motherLabel>stack->GetNtrack()) return -1;
  TParticle *motherPart = stack->Particle(motherLabel);
  if(!motherPart) return -1;

  Int_t pdgCode = TMath::Abs(motherPart->GetPdgCode());

  if(pdgCode==kPDGelectron) {
    if(fDebugLevel>=10) printf("particle label: %d...(motherLabel=%d : motherPdg=%d)  grandmother's pdg code was returned...%d \n",
			       label, motherLabel, pdgCode, ElectronFromSource(stack, motherLabel));
    return ElectronFromSource(stack, motherLabel);
  }

  return pdgCode;    

}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::ElectronFromCharm(AliStack * const stack, Int_t label) const
{
  //
  //  separate electron: kEleC from c->eX, kEleBC from b->c->eX
  //

  Int_t motherLabel = label;
  TParticle *motherParticle = stack->Particle(motherLabel);  // mother part
  if(!motherParticle) return -1;
  Int_t gMotherLabel = motherParticle->GetFirstMother();  // grand mother
  if(gMotherLabel<0 || gMotherLabel>stack->GetNtrack()) return -1;
  
  TParticle *gMotherPart = stack->Particle(gMotherLabel);
  if(!gMotherPart) return -1;
  
  Int_t pdgCode = gMotherPart->GetPdgCode();
  if(TMath::Abs(pdgCode%10000)/100==kPDGbeauty ||    // for special intermediate meson states: like 100553
     TMath::Abs(pdgCode)/1000==kPDGbeauty || TMath::Abs(pdgCode)/100==kPDGbeauty || TMath::Abs(pdgCode)==kPDGbeauty) {
    if(fDebugLevel>=10)  printf("this electron label %d is from mother %d, and finally from %d\n", label, ElectronFromSource(stack, label), pdgCode);
    return kEleBC;
  }  // for sure it is from BC
  
  else 
    if(TMath::Abs(pdgCode%10000)/100==kPDGcharm ||  // for special intermediate meson states: like 100443, 10443, 204*3 (*=1, 2, 3, 4)
       TMath::Abs(pdgCode)/1000==kPDGcharm || 
       TMath::Abs(pdgCode)/100==kPDGcharm || 
       TMath::Abs(pdgCode)==kPDGcharm){
      
      if(CharmFromBeauty(stack, gMotherLabel)!=-1)
	return kEleBC;
      else
	return kEleC;
      
    }
  
    else
      return -1;
     
}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::CharmFromBeauty(AliStack * const stack, Int_t hfLabel) const
{

  //
  //  check if charm meson/hadron is from beauty decay
  //  -1, not from beauty  
  //  return the label of the mother of this charm hadron
  //

  Int_t charmLabel = hfLabel;
  
  //  AliStack *stack = fMC->Stack();
  if(charmLabel<0 || charmLabel>stack->GetNtrack()) return -1;
  TParticle *particle = stack->Particle(charmLabel);
  if(!particle) return -1;

  Int_t motherCharmLabel = particle->GetFirstMother();
  if(motherCharmLabel<0 || motherCharmLabel>stack->GetNtrack()) return -1;

  TParticle *motherCharmPart = stack->Particle(motherCharmLabel);
  if(!motherCharmPart) return -1;

  Int_t pdgCode = motherCharmPart->GetPdgCode();
  
  if(TMath::Abs(pdgCode%10000)/100==kPDGbeauty ||   // for special intermediate meson states: like 100553
     TMath::Abs(pdgCode)/1000==kPDGbeauty || 
     TMath::Abs(pdgCode)/100==kPDGbeauty || 
     TMath::Abs(pdgCode)==kPDGbeauty) 
    return motherCharmLabel;  
  else 
    if(TMath::Abs(pdgCode%10000)/100==kPDGcharm ||    // for special intermediate meson states: like 100443, 10443, 204*3 (*=1, 2, 3, 4)
       TMath::Abs(pdgCode)/1000==kPDGcharm || 
       TMath::Abs(pdgCode)/100==kPDGcharm) 
      return CharmFromBeauty(stack, motherCharmLabel);   // loop over to see if charm is from beauty -- if yes, return the label of mother of the charm
  
    else
      return -1;
}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::ElePhotonDirect(AliStack * const stack, Int_t label) const 
{
  //
  // electron is from photon, and check if this photon is direct
  //

  Int_t eleLabel = label;
  
  TParticle *particle = stack->Particle(eleLabel);
  if(particle->GetFirstMother()<0 || particle->GetFirstMother()>stack->GetNtrack())
    return -1;
  
  TParticle *photonPart = stack->Particle(particle->GetFirstMother());
  if(!photonPart) 
    return -1;
  Int_t motherPhotonLabel = photonPart->GetFirstMother();
  if(motherPhotonLabel<0 || motherPhotonLabel>stack->GetNtrack())
    return -1;

  TParticle *motherPhotonPart = stack->Particle(motherPhotonLabel);
  if(!motherPhotonPart) 
    return -1;
  
  Int_t pdgMotherPhoton = motherPhotonPart->GetPdgCode();
  if(TMath::Abs(pdgMotherPhoton)<=10 || TMath::Abs(pdgMotherPhoton)==21)
    return 1;

  else 
    return -1;


}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::GetMCpid(AliStack* const stack, Int_t partLabel) const
{

  //
  // Simply pdg code
  //
  
  Int_t label = partLabel;
//   AliStack* stack = fMC->Stack();
  if((label < 0) || (label >= stack->GetNtrack())) return -1;  

  // MC Information
  TParticle * particle = stack->Particle(label);
  if(!particle) return -1;
  Int_t pdg = particle->GetPdgCode();

  return pdg;
 
}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::GetMotherLabel(AliStack *const stack, Int_t eleLabel) const 
{
  //
  // Simply label of mother
  //


  Int_t label = eleLabel;
  //  AliStack* stack = fMC->Stack();
  if((label < 0) || (label >= stack->GetNtrack())) return -1;  

  // MC Information
  TParticle * particle = stack->Particle(label);
  if(!particle) return -1;

  return particle->GetFirstMother();
 
}


//__________________________________________________________
Float_t AliHFEdisplacedElectrons::GetRapidity(TParticle *part) const
{
  // return rapidity
  
  Float_t rapidity = -999;        
  if((part->Energy() - part->Pz())*(part->Energy() + part->Pz())>0)
    rapidity = 0.5*(TMath::Log((part->Energy()+part->Pz()) / (part->Energy()-part->Pz()))); 
  
  return rapidity;
}


//__________________________________________________________
Float_t AliHFEdisplacedElectrons::GetTrackRapidity(AliESDtrack * const track) const
{
  // return rapidity of electron
  
  Float_t px = track->Px();
  Float_t py = track->Py();
  Float_t pz = track->Pz();

  // electron mass 0.00051099906 GeV/c2
  TParticlePDG* electron = TDatabasePDG::Instance()->GetParticle(kPDGelectron);
  Double_t mass = electron->Mass();
  
  Float_t en = TMath::Sqrt(px*px + py*py + pz*pz + mass*mass);

  Float_t rapidity = -999;        
  if((en - pz)*(en + pz)>0)
    rapidity = 0.5*(TMath::Log((en - pz)*(en + pz))); 
  
  return rapidity;
}
