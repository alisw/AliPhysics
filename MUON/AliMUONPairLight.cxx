/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *      SigmaEffect_thetadegrees                                                                  *
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

//===================================================================
//This class was prepared by INFN Cagliari, July 2006
//(authors: H.Woehri, A.de Falco)
// 
// Compact information for the generated muon pairs in the MUON arm 
// useful at the last stage of the analysis chain
// Pairs are built with two AliMUONTrackLight objects 
// Using the class AliMUONTrackLight this class combines the decay
// information ("history") of the reconstructed tracks and fills
// a series of flags for the formed reconstructed dimuon:
// fIsCorrelated, fCreationProcess, fIsFeedDown, ...
// for information about the dimuon, use PrintInfo with the appropriate
// printflag
// To be used together with AliMUONTrackLight
//===================================================================


//MUON classes
#include "AliMUONPairLight.h"
//Root classes
#include "TString.h"

ClassImp(AliMUONPairLight) 

//====================================
AliMUONPairLight::AliMUONPairLight() : 
  TObject(), 
  fMu0(),
  fMu1(), 
  fCreationProcess(-999),
  fIsCorrelated(kFALSE), 
  fCauseOfCorrelation (-1),
  fIsFeedDown(kFALSE)
{
  /// default constructor
  ; 
}

//====================================

AliMUONPairLight::AliMUONPairLight(AliMUONPairLight &dimuCopy) 
  : TObject(dimuCopy),
    fMu0(dimuCopy.fMu0),
    fMu1(dimuCopy.fMu1), 
    fCreationProcess(dimuCopy.fCreationProcess),
    fIsCorrelated(dimuCopy.fIsCorrelated), 
    fCauseOfCorrelation (dimuCopy.fCauseOfCorrelation),
    fIsFeedDown(dimuCopy.fIsFeedDown)
{ 
/// copy constructor
///   fMu0 = AliMUONTrackLight(dimuCopy.fMu0); 
///   fMu1 = AliMUONTrackLight(dimuCopy.fMu1); 
///   fIsCorrelated = dimuCopy.fIsCorrelated;
///   fCauseOfCorrelation = dimuCopy.fCauseOfCorrelation;
///   fCreationProcess = dimuCopy.fCreationProcess;
///   fIsFeedDown = dimuCopy.fIsFeedDown;
  ;
}

//====================================

AliMUONPairLight::~AliMUONPairLight(){
  /// destructor
}

//====================================

Bool_t AliMUONPairLight::IsAResonance(){
  /// checks if muon pair comes from a resonance decay  
  if (!fIsCorrelated) return kFALSE;   //if muons not correlated, cannot be a resonance
  //if muons are correlated, check if the PDG of the
  //common mother is a resonance
  Int_t pdg = GetCauseOfCorrelation();
  if(pdg < 10) return kFALSE;
  Int_t id=pdg%100000;
  if(((id-id%10)%110)) return kFALSE;
  else return kTRUE;
  printf("<AliMUONPairLight::IsAResonance> arriving after this piece of code\n");

}

//====================================

AliMUONTrackLight* AliMUONPairLight::GetMuon(Int_t index)  { 
  /// return muon 0 or 1
   if (index==0) return &fMu0;
   else if (index==1) return &fMu1; 
   else{ printf ("Index can be either 0 or 1\n"); return 0;}
   //   else return &fMu1; 
}

//====================================

Int_t AliMUONPairLight::GetMuonMotherPDG(Int_t imuon, Int_t mother) { 
  /// return muon mother pdg code
  if (imuon==0) return fMu0.GetParentPDGCode(mother); 
  else if (imuon==1) return fMu1.GetParentPDGCode(mother); 
  else { printf ("Index must be only 0 or 1\n"); return -999; } 
}

//====================================
void AliMUONPairLight::SetProcess(){
  /// finds the process related to the muon pair (open charm/beauty, resonance, 
  /// uncorrelated...) 
  AliMUONTrackLight *mu1 = &fMu0;
  AliMUONTrackLight *mu2 = &fMu1;

  //1.) check if one of the tracks is not a muon
  if(IsOneTrackNotAMuon()) { 
    this->SetCorrelated(kFALSE);
    return;
  }

  // check if the two muons are correlated
  // first check if they come from the same hadron (resonance or beauty/charm meson)
  Int_t npar1 = mu1->GetNParents(); 
  Int_t npar2 = mu2->GetNParents(); 
  //  for (Int_t imoth1 = 0; imoth1<npar1; imoth1++) { 
  for (Int_t imoth1 = npar1-1; imoth1>=0; imoth1--) { 
    Int_t lineMo1 = mu1->GetParentPythiaLine(imoth1);
    //    for (Int_t imoth2 = 0; imoth2<npar2; imoth2++) { 
    for (Int_t imoth2 = npar2-1; imoth2>=0; imoth2--) { 
      Int_t lineMo2 = mu2->GetParentPythiaLine(imoth2);
      if(lineMo1 == lineMo2) { 
	this->SetCorrelated(kTRUE); 
	this->SetCauseOfCorrelation(mu1->GetParentPDGCode(imoth1));
	if(!IsAResonance()) fCreationProcess = 3; 
	else fCreationProcess = -1;
	return;   // <<<<---------------- RETURN? 
      }
    }
  }  
  
  // if Open Beauty/Charm we can have 3 creation processes 
  // (pair creation [0], gluon splitting [1] or flavour excitation [2])
  //
  // 1.) gluon splitting: gluon (stored with index 2, id=21) must be the same 
  if (mu1->GetQuarkPythiaLine(2) == mu2->GetQuarkPythiaLine(2) && mu1->GetQuarkPDGCode(2) == 21) { 
    this->SetCorrelated(kTRUE); 
    if(GetCauseOfCorrelation() == -1){
      this->SetCauseOfCorrelation(mu1->GetQuarkPDGCode(2));
    }
    fCreationProcess = 1; 
    return;
  }
  
  // 2.) pair creation: if pythia line 6 of one track *and* pythia line 7 of second track
  // are filled with a Q and Qbar
  Int_t line1 = mu1->GetQuarkPythiaLine(2); //[2] ... very first quark
  Int_t line2 = mu2->GetQuarkPythiaLine(2); 
  Int_t flavour1 = TMath::Abs(mu1->GetQuarkPDGCode(2)); //[2] ... very first quark
  Int_t flavour2 = TMath::Abs(mu2->GetQuarkPDGCode(2)); 
  if ((line1 == 6 || line1 == 7) && (line2 == 6 || line2 == 7)) { 
    if((flavour1 == 4 || flavour1 == 5) && (flavour2 == 4 || flavour2 == 5)){
      this->SetCorrelated(kTRUE);
      fCreationProcess = 0; 
      return;
    }
  }

  // 3.)flavour excitation:
  if((((line1 == 6 || line1 == 7) && (flavour1 == 4 || flavour1 == 5)) && //first track has Q in line 6 or 7
      ((line2 == 2 || line2 == 3) && (flavour2 == 21 || flavour2 < 10))) || //second track has a g/q in line 2 or 3
      (((line2 == 6 || line2 == 7) && (flavour2 == 4 || flavour2 == 5)) && //or the same,
      ((line1 == 2 || line1 == 3) && (flavour1 == 21 || flavour1 < 10)))){ // swapping the track's indices

//     printf("candidate for flavour excitation\n");

    //Hermine: is it needed to check the lines 4-7???

    //now, we have a candidate for flavour excitation
    //we must also verify if in the Pythia listing
    //the "incoming" lines 4 and 5 and "outgoing" lines 6 and 7
    //are filled with g Q each: e.g. 4(g),5(Q),5(g),6(Q)
//     Int_t pdg4 = TMath::Abs(stack()->Particle(4)->GetPdgCode());
//     Int_t pdg5 = TMath::Abs(stack()->Particle(5)->GetPdgCode());
//     Int_t pdg6 = TMath::Abs(stack()->Particle(6)->GetPdgCode());
//     Int_t pdg7 = TMath::Abs(stack()->Particle(7)->GetPdgCode());
//     if((pdg4 == 21 && pdg5 < 10 || pdg5 == 21 && pdg4 < 10 ) &&
//        (pdg6 == 21 && pdg7 < 10 || pdg7 == 21 && pdg6 < 10 )){
//       this->PrintInfo("H");
      this->SetCorrelated(kTRUE);
      fCreationProcess = 2;
      return;
//     }
  }
}

//====================================
void AliMUONPairLight::SetMuons(AliMUONTrackLight mu0, AliMUONTrackLight mu1){
  /// set the two muons 
  fMu0 = mu0; 
  fMu1 = mu1; 
  this->SetProcess();
} 

//====================================
void AliMUONPairLight::PrintInfo(Option_t* opt){
  /// print information about muon pairs
  /// Options: 
  /// - "H" single muons' decay histories
  /// - "K" dimuon kinematics
  /// - "F" dimuon flags
  /// - "A" all variables
  TString options(opt);
  options.ToUpper();

  if(options.Contains("H") || options.Contains("A")){//muon decay histories

    AliMUONTrackLight *mu1 = &fMu0;
    AliMUONTrackLight *mu2 = &fMu1;

    printf("========= History =======================\n");
    printf("first muon");
    mu1->PrintInfo("H");
    printf("second muon");
    mu2->PrintInfo("H");
    printf("=========================================\n");
  }
  if(options.Contains("F") || options.Contains("A")){//flags
    printf("the flags set for this muon pair are:\n");
    printf("=====================================\n");
    if(this->IsOneTrackNotAMuon()) printf("(*) one rec. track is not a muon\n");
    fIsCorrelated ? printf("(*) it is a correlated pair\n") : printf("(*) it is not a correlated pair\n");
    if(IsOpenCharm()) printf("(*) correlated open charm: ");
    if(IsOpenBeauty()) printf("(*) correlated open beauty: ");
    if(IsOpenCharm() || IsOpenBeauty()){
      switch(fCreationProcess){
      case 0:
	printf("pair creation");
	break;
      case 1:
	printf("gluon splitting");
	break;
      case 2:
	printf("flavour excitation");
	break;
      case 3:
	printf("both muons come from same fragmented mother");
	break;
      }
      if(this->GetMuon(0)->GetOscillation() || this->GetMuon(1)->GetOscillation()) 
	printf("... where oscillation occured\n");
      else{
	if(IsOpenBeauty())
	  printf(" (no oscillation)\n");
	else
	  printf("\n");
      }
    }
    IsAResonance() ? printf("(*) it is a resonance: %d\n", this->GetMuonMotherPDG(fIsFeedDown)) : printf("(*) it is not a resonance\n");
    fIsFeedDown ? printf("(*) mother has feed-down: %d --> %d\n", this->GetMuonMotherPDG(1), this->GetMuonMotherPDG(0)) : printf("(*) no feed-down\n");
    printf("=====================================\n");
  }
  if(options.Contains("K") || options.Contains("A")){//dimuon kinematics
    Double_t *vtx = this->GetMuon(0)->GetVertex();
    TLorentzVector momRec = this->GetPRec();
    TLorentzVector momGen = this->GetPGen();
    printf("the dimuon charge is %d\n", this->GetCharge());
    printf("primary Vertex: Vx = %1.3f, Vy = %1.3f, Vz = %1.3f\n", vtx[0], vtx[1], vtx[2]);
    printf("Generated:     Px = %1.3f, Py = %1.3f, Pz = %1.3f\n", momGen.Px(), momGen.Py(), momGen.Pz());
    printf("Reconstructed: Px = %1.3f, Py = %1.3f, Pz = %1.3f\n", momRec.Px(), momRec.Py(), momRec.Pz());
    //rapidity, pT, angles, ...
    printf("Rec. variables: pT %1.3f, pseudo-rapidity %1.3f, openingAngle %1.3f (%1.3f degree), theta %1.3f (%1.3f degree), phi %1.3f (%1.3f degree)\n", 
	   momRec.Pt(), momRec.Eta(), 
	   TMath::Pi()/180.*this->GetOpeningAngle(), this->GetOpeningAngle(), 
	   momRec.Theta(), 180./TMath::Pi() * momRec.Theta(), 
	   momRec.Phi(), 180./TMath::Pi() * momRec.Phi());
  }
}

Double_t AliMUONPairLight::GetOpeningAngle() { 
  // opening angle between the two muons in the lab frame (in degrees)
  TLorentzVector pRecMu0 =  fMu0.GetPRec();
  TLorentzVector pRecMu1 =  fMu1.GetPRec();
  TVector3 pRecMu03 = pRecMu0.Vect();
  TVector3 pRecMu13 = pRecMu1.Vect();
  Double_t scalar = pRecMu03.Dot(pRecMu13);
  Double_t modMu0 = pRecMu03.Mag();
  Double_t modMu1 = pRecMu13.Mag();
  Double_t theta = (TMath::ACos(scalar/(modMu0*modMu1)))*(180./TMath::Pi());
  return theta; 
}
