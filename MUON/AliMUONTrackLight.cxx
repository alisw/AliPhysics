/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *      SigmaEffect_thetadegrees                                          *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------------------------------
// Compact information for the muon generated tracks in the MUON arm 
// useful at the last stage of the analysis chain
// provides a link between the reconstructed track and the generated particle 
// stores kinematical information at gen. and rec. level and 
// the decay history of the muon, allowing the identification of the 
// mother process 
// 
// To be used together with AliMUONPairLight
//
// This class was prepared by INFN Cagliari, July 2006
// (authors: H.Woehri, A.de Falco)
//-----------------------------------------------------------------------------

// 13 Nov 2007:
// Added a temporary fix to FindRefTrack to be able to handle reconstructed tracks
// generated from ESD muon track information. The problem is that the ESD data at
// the moment only contains the first hit on chamber 1. Hopefully in the near future
// this will be fixed and all hit information will be available.
//  - Artur Szostak <artursz@iafrica.com>

#include "AliMUONTrackLight.h"
#include "AliMUONTrack.h"
#include "AliMUONConstants.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTrackParam.h"

#include "AliESDMuonTrack.h"
#include "AliStack.h"
#include "AliLog.h"

#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TString.h"

#include <cstdio>

ClassImp(AliMUONTrackLight) 

//===================================================================

AliMUONTrackLight::AliMUONTrackLight() 
  : TObject(), 
    fPrec(), 
    fIsTriggered(kFALSE),
    fCharge(-999), 
    fChi2(-1), 
    fCentr(-1),
    fPgen(), 
    fTrackPythiaLine(-999),
    fTrackPDGCode(-999),
    fOscillation(kFALSE), 
    fNParents(0),
    fWeight(1)    
{
  /// default constructor
  fPgen.SetPxPyPzE(0.,0.,0.,0.); 
  fPrec.SetPxPyPzE(0.,0.,0.,0.); 
  for (Int_t i=0; i<3; i++) fXYZ[i]=-999; 
  for (Int_t npar = 0; npar < fgkNParentsMax; npar++){
    fParentPDGCode[npar] = -1; 
    fParentPythiaLine[npar] = -1;
  }
  for (Int_t i = 0; i < 4; i++){
    fQuarkPDGCode[i] = -1; 
    fQuarkPythiaLine[i] = -1; 
  }
}

//============================================
AliMUONTrackLight::AliMUONTrackLight(const AliMUONTrackLight &muonCopy) 
  : TObject(muonCopy), 
    fPrec(muonCopy.fPrec), 
    fIsTriggered(muonCopy.fIsTriggered),
    fCharge(muonCopy.fCharge), 
    fChi2(muonCopy.fChi2), 
    fCentr(muonCopy.fCentr),
    fPgen(muonCopy.fPgen), 
    fTrackPythiaLine(muonCopy.fTrackPythiaLine),
    fTrackPDGCode(muonCopy.fTrackPDGCode),
    fOscillation(muonCopy.fOscillation), 
    fNParents(muonCopy.fNParents),
    fWeight(muonCopy.fWeight)
{
  /// copy constructor
  for (Int_t i=0; i<3; i++) fXYZ[i]=muonCopy.fXYZ[i]; 
  for (Int_t npar = 0; npar < fgkNParentsMax; npar++){
    fParentPDGCode[npar] = muonCopy.fParentPDGCode[npar]; 
    fParentPythiaLine[npar] = muonCopy.fParentPythiaLine[npar];
  }
  for (Int_t i = 0; i < 4; i++){
    fQuarkPDGCode[i] = muonCopy.fQuarkPDGCode[i]; 
    fQuarkPythiaLine[i] = muonCopy.fQuarkPythiaLine[i]; 
  }
}

//============================================
AliMUONTrackLight::AliMUONTrackLight(AliESDMuonTrack* muonTrack)
  : TObject(), 
    fPrec(), 
    fIsTriggered(kFALSE),
    fCharge(-999), 
    fChi2(-1),
    fCentr(-1),
    fPgen(), 
    fTrackPythiaLine(-999),
    fTrackPDGCode(-999),
    fOscillation(kFALSE), 
    fNParents(0),
    fWeight(1)
{ 
  /// constructor
  fPgen.SetPxPyPzE(0.,0.,0.,0.); 
  for (Int_t npar = 0; npar < fgkNParentsMax; npar++){
    fParentPDGCode[npar] = -1; 
    fParentPythiaLine[npar] = -1;
  }
  for (Int_t i = 0; i < 4; i++){
    fQuarkPDGCode[i] = -1; 
    fQuarkPythiaLine[i] = -1; 
  }
  FillFromESD(muonTrack);
}

//============================================
AliMUONTrackLight::~AliMUONTrackLight()
{
/// Destructor
} 

//============================================
AliMUONTrackLight& AliMUONTrackLight::operator=(const AliMUONTrackLight& muonCopy)
{
  // check assignment to self
  if (this == &muonCopy) return *this;

  // base class assignment
  TObject::operator=(muonCopy);

  // assignment operator
  fPrec = muonCopy.fPrec; 
  fIsTriggered = muonCopy.fIsTriggered;
  fCharge = muonCopy.fCharge; 
  fChi2 = muonCopy.fChi2; 
  fCentr = muonCopy.fCentr;
  fPgen = muonCopy.fPgen; 
  fTrackPythiaLine = muonCopy.fTrackPythiaLine;
  fTrackPDGCode = muonCopy.fTrackPDGCode;
  fOscillation = muonCopy.fOscillation; 
  fNParents = muonCopy.fNParents;
  fWeight = muonCopy.fWeight;
  
  for (Int_t i=0; i<3; i++) fXYZ[i]=muonCopy.fXYZ[i]; 
  for (Int_t npar = 0; npar < fgkNParentsMax; npar++){
    fParentPDGCode[npar] = muonCopy.fParentPDGCode[npar]; 
    fParentPythiaLine[npar] = muonCopy.fParentPythiaLine[npar];
  }
  for (Int_t i = 0; i < 4; i++){
    fQuarkPDGCode[i] = muonCopy.fQuarkPDGCode[i]; 
    fQuarkPythiaLine[i] = muonCopy.fQuarkPythiaLine[i]; 
  }
}    

//============================================

void AliMUONTrackLight::FillFromAliMUONTrack(AliMUONTrack *trackReco,Double_t zvert){
  /// this method sets the muon reconstructed momentum according to the value given by AliMUONTrack
  AliMUONTrackParam* trPar = trackReco->GetTrackParamAtVertex();
  if (!trPar) {
    AliError("The track must contain the parameters at vertex");
    return;
  }
  this->SetCharge(Int_t(TMath::Sign(1.,trPar->GetInverseBendingMomentum())));
  this->SetPxPyPz(trPar->Px(),trPar->Py(), trPar->Pz()); 
  this->SetTriggered(trackReco->GetMatchTrigger()); 
  
  Double_t xyz[3] = { trPar->GetNonBendingCoor(), 
		      trPar->GetBendingCoor(),
		      trPar->GetZ()};
  if (zvert!=-9999) xyz[2] = zvert;
  this->SetVertex(xyz); 
}

//============================================
void AliMUONTrackLight::FillFromESD(AliESDMuonTrack* muonTrack,Double_t zvert){
  /// computes prec and charge from ESD track
  Double_t mumass = TDatabasePDG::Instance()->GetParticle(13)->Mass(); 
  Double_t thetaX = muonTrack->GetThetaX();
  Double_t thetaY = muonTrack->GetThetaY();
  Double_t tanthx = TMath::Tan(thetaX);
  Double_t tanthy = TMath::Tan(thetaY);
  Double_t pYZ    =  1./TMath::Abs(muonTrack->GetInverseBendingMomentum());
  Double_t pz     = - pYZ / TMath::Sqrt(1.0 + tanthy * tanthy);
  Double_t px     = pz * tanthx;
  Double_t py     = pz * tanthy;
  fCharge   = Int_t(TMath::Sign(1.,muonTrack->GetInverseBendingMomentum()));
  Double_t energy = TMath::Sqrt(mumass * mumass + px*px + py*py + pz*pz);
  fPrec.SetPxPyPzE(px,py,pz,energy);
  // get the position
  fXYZ[0] = muonTrack->GetNonBendingCoor();
  fXYZ[1] = muonTrack->GetBendingCoor();
  if (zvert==-9999) fXYZ[2] = muonTrack->GetZ();
  else fXYZ[2] = zvert;
  // get the chi2 per d.o.f.
  fChi2 = muonTrack->GetChi2()/ (2.0 * muonTrack->GetNHit() - 5);
  fIsTriggered = muonTrack->GetMatchTrigger();
}

//============================================
void AliMUONTrackLight::SetPxPyPz(Double_t px, Double_t py, Double_t pz){ 
  /// set the reconstructed 4-momentum, assuming the particle is a muon
  Double_t mumass = TDatabasePDG::Instance()->GetParticle(13)->Mass(); 
  Double_t energy = TMath::Sqrt(mumass * mumass + px*px + py*py + pz*pz);
  fPrec.SetPxPyPzE(px,py,pz,energy);
}

//============================================
void AliMUONTrackLight::FillMuonHistory(AliStack *stack, TParticle *part){
  /// scans the muon history to determine parents pdg code and pythia line
  Int_t countP = -1;
  Int_t parents[10], parLine[10];
  Int_t lineM = part->GetFirstMother();//line in the Pythia output of the particle's mother

  TParticle *mother;
  Int_t status=-1, pdg=-1;
  while(lineM >= 0){
    
    mother = stack->Particle(lineM); //direct mother of rec. track
    pdg = mother->GetPdgCode();//store PDG code of first mother
    // break if a string, gluon, quark or diquark is found 
    if(pdg == 92 || pdg == 21 || TMath::Abs(pdg) < 10 || IsDiquark(pdg)) break;
    parents[++countP] = pdg;
    parLine[countP] = lineM;
    status = mother->GetStatusCode();//get its status code to check if oscillation occured
    if(IsB0(parents[countP]) && status == 12) this->SetOscillation(kTRUE);
    lineM = mother->GetFirstMother();
  }
  //store all the fragmented parents in an array:
  for(int i = 0; i <= countP; i++){
    this->SetParentPDGCode(i,parents[countP-i]);
    this->SetParentPythiaLine(i,parLine[countP-i]);
  }
  fNParents = countP+1;
  countP = -1;

  //and store the lines of the string and further quarks in another array:
  while(lineM >= 0){
    mother = stack->Particle(lineM);
    pdg = mother->GetPdgCode();
    //now, get information before the fragmentation
    this->SetQuarkPythiaLine(++countP, lineM);//store the line of the string in index 0
    this->SetQuarkPDGCode(countP, pdg);//store the pdg of the quarks in index 1,2
    lineM = mother->GetFirstMother();
  }
  
  //check if in case of HF production, the string points to the correct end
  //and correct it in case of need:
  countP = 1;
  for(int par = 0; par < 4; par++){
    if(TMath::Abs(this->GetQuarkPDGCode(par)) < 6){
      countP = par; //get the quark just before hadronisation
      break;
    }
  }
  if(this->GetQuarkPythiaLine(countP) > -1 && (this->GetParentFlavour(0)==4 || this->GetParentFlavour(0)==5)){
    if(this->GetParentFlavour(0) != TMath::Abs(this->GetQuarkPDGCode(countP))){

      AliWarning(Form("quark flavour of parent and that of quark do not correspond: %d %d --> correcting\n",
          this->GetParentFlavour(0), TMath::Abs(this->GetQuarkPDGCode(countP)))
        );
      
      pdg = this->GetQuarkPDGCode(countP);
      Int_t line = this->GetQuarkPythiaLine(countP);
      this->ResetQuarkInfo();
      while(TMath::Abs(pdg) != this->GetParentFlavour(0)){//pdg of q,g in Pythia listing following the wrong string end
                                                        //must coincide with the flavour of the last fragmented mother

	pdg = stack->Particle(++line)->GetPdgCode();    
      }
      //now, we have the correct string end and the correct line
      //continue to fill again all parents and corresponding lines
      while(line >= 0){
	mother = stack->Particle(line);//get again the mother
	pdg = mother->GetPdgCode();
	this->SetQuarkPythiaLine(countP, line);
	this->SetQuarkPDGCode(countP++, pdg);
	line = mother->GetFirstMother();
      }
      this->PrintInfo("h");
    }//mismatch
  }
}

//====================================
void AliMUONTrackLight::ResetQuarkInfo(){
  /// resets parton information
  for(int pos = 1; pos < 4; pos++){//[0] is the string
    this->SetQuarkPDGCode(pos,-1);
    this->SetQuarkPythiaLine(pos,-1);
  }
}

//====================================
Bool_t AliMUONTrackLight::IsB0(Int_t intTest) const {
  /// checks if the particle is a B0 
  Int_t bMes0[2] = {511,531};//flavour code of B0d and B0s
  Bool_t answer = kFALSE;
  for(int i = 0; i < 2; i++){
    if(TMath::Abs(intTest) == bMes0[i]){
      answer = kTRUE;
      break;
    }
  }
  return answer;
}
//====================================
Bool_t AliMUONTrackLight::IsMotherAResonance(Int_t index) const {
  /// checks if mother is a resonance
  Int_t intTest = GetParentPDGCode(index); 
  // the resonance pdg code is built this way
  // x00ffn where x=0,1,.. (1S,2S... states), f=quark flavour 
  Int_t id=intTest%100000; 
  return (!((id-id%10)%110));
}
//====================================
Int_t AliMUONTrackLight::GetParentFlavour(Int_t idParent) const {
  /// returns the flavour of parent idParent (idParent=0 is the oldest 
  /// hadronized parent)
  Int_t pdg = GetParentPDGCode(idParent); 
  Int_t quark = TMath::Abs(pdg/100);
  if(quark > 9) quark = quark/10;
  return quark;
}

//====================================
void AliMUONTrackLight::PrintInfo(const Option_t* opt){
  /// prints information about the track: 
  /// - "H" muon's decay history
  /// - "K" muon kinematics
  /// - "A" all variables
  TString options(opt);
  options.ToUpper();

  if(options.Contains("H") || options.Contains("A")){ //muon decay history
    char *name= new char[100];
    TString pdg = "", line = "";
    for(int i = 3; i >= 0; i--){
      if(this->GetQuarkPythiaLine(i)>= 0){
	snprintf(name, 100, "%4d --> ", this->GetQuarkPythiaLine(i));
	line += name;
	snprintf(name, 100, "%4d --> ", this->GetQuarkPDGCode(i));
	pdg += name;
      }
    }
    for(int i = 0; i < fNParents; i++){ 
      if(this->GetParentPythiaLine(i)>= 0){
	snprintf(name, 100, "%7d --> ", this->GetParentPythiaLine(i));
	line += name;
	snprintf(name, 100, "%7d --> ", this->GetParentPDGCode(i));
	pdg += name;
      }
    }
    snprintf(name, 100, "%4d", this->GetTrackPythiaLine()); line += name;
    snprintf(name, 100, "%4d", this->GetTrackPDGCode()); pdg += name;

    printf("\nmuon's decay history:\n");
    printf(" PDG: %s\n", pdg.Data());
    printf("line: %s\n", line.Data());
  }
  if(options.Contains("K") || options.Contains("A")){ //muon kinematic

    Int_t charge = this->GetCharge();
    Double_t *vtx = this->GetVertex();
    TLorentzVector momRec = this->GetPRec();
    TLorentzVector momGen = this->GetPGen();
    printf("the track's charge is %d\n", charge);
    printf("Primary vertex: Vx = %1.3f, Vy = %1.3f, Vz = %1.3f\n", vtx[0], vtx[1], vtx[2]);
    printf("Generated:     Px = %1.3f, Py = %1.3f, Pz = %1.3f\n", momGen.Px(), momGen.Py(), momGen.Pz());
    printf("Reconstructed: Px = %1.3f, Py = %1.3f, Pz = %1.3f\n", momRec.Px(), momRec.Py(), momRec.Pz());
    printf("Rec. variables: pT %1.3f, pseudo-rapidity %1.3f, theta %1.3f (%1.3f degree), phi %1.3f (%1.3f degree)\n", 
	   momRec.Pt(), momRec.Eta(), momRec.Theta(), 180./TMath::Pi() * momRec.Theta(), 
	   momRec.Phi(), 180./TMath::Pi() * momRec.Phi());
  }
}
//====================================
Bool_t AliMUONTrackLight::IsParentPionOrKaon(Int_t idparent){
  /// checks if a muon comes from a pion or kaon or a particle that decays into one of these two
  Int_t pdg = this->GetParentPDGCode(idparent); 
  if (TMath::Abs(pdg)==211 || //pi+
      TMath::Abs(pdg)==321 || //K+
      TMath::Abs(pdg)==213 || //rho+
      TMath::Abs(pdg)==311 || //K0
      TMath::Abs(pdg)==313 || //K*0
      TMath::Abs(pdg)==323    //K*+
      ) { 
    return kTRUE;
  }
  else return kFALSE;
}
//====================================
Bool_t AliMUONTrackLight::IsDiquark(Int_t pdg) const{
  /// check if the provided pdg code corresponds to a diquark 
  pdg = TMath::Abs(pdg);
  if((pdg > 1000) && (pdg%100 < 10)) return kTRUE;
  else return kFALSE;
}
