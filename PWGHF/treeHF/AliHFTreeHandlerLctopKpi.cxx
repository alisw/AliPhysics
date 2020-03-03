/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerLctopKpi
// \brief helper class to handle a tree for Lc->pKpi cut optimisation and MVA analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
// F. Prino, prino@to.infn.it
// L. Vermunt, luuk.vermunt@cern.ch
// L. van Doremalen, lennart.van.doremalen@cern.ch
// J. Norman, jaime.norman@cern.ch
// G. Luparello, grazia.luparello@cern.ch
/////////////////////////////////////////////////////////////

#include <TString.h>
#include "AliHFTreeHandlerLctopKpi.h"
#include "AliAODRecoDecayHF3Prong.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerLctopKpi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerLctopKpi::AliHFTreeHandlerLctopKpi():
  AliHFTreeHandler(),
  fSigmaVertex(-9999.),
  fDist12toPrim(-9999.),
  fDist23toPrim(-9999.),
  fNormd0MeasMinusExp(-9999.),
  fSumImpParProngs(-9999.),
  fResonantDecayType(0),
  fResonantDecayTypeMC(0)
{
  //
  // Default constructor
  //

  fNProngs=3; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fImpParProng[iProng] = -9999.;
    fDCAProng[iProng] = -9999.;
  }
}

//________________________________________________________________
AliHFTreeHandlerLctopKpi::AliHFTreeHandlerLctopKpi(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fSigmaVertex(-9999.),
  fDist12toPrim(-9999.),
  fDist23toPrim(-9999.),
  fNormd0MeasMinusExp(-9999.),
  fSumImpParProngs(-9999.),
  fResonantDecayType(0),
  fResonantDecayTypeMC(0)
{
  //
  // Standard constructor
  //

  fNProngs=3; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fImpParProng[iProng] = -9999.;
    fDCAProng[iProng] = -9999.;
  }
}

//________________________________________________________________
AliHFTreeHandlerLctopKpi::~AliHFTreeHandlerLctopKpi()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerLctopKpi::BuildTree(TString name, TString title) 
{
  fIsMCGenTree=false;

  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  //set common variables
  AddCommonDmesonVarBranches();

  //set Lc variables
  fTreeVar->Branch("sig_vert",&fSigmaVertex);
  fTreeVar->Branch("dist_12",&fDist12toPrim);
  fTreeVar->Branch("dist_23",&fDist23toPrim);
  fTreeVar->Branch("max_norm_d0d0exp",&fNormd0MeasMinusExp);
  fTreeVar->Branch("sum_d0d0_prongs",&fSumImpParProngs);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
    fTreeVar->Branch(Form("dca_prong%d",iProng),&fDCAProng[iProng]);
  }
  fTreeVar->Branch("resonant_decay_mc",&fResonantDecayType);

  //set single-track variables
  AddSingleTrackBranches();
  if (fFillJets) AddJetBranches();

  //set PID variables
  if(fPidOpt!=kNoPID) AddPidBranches(true,true,true,true,true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerLctopKpi::SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo, AliPIDResponse *pidrespo)
{
  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandType&kSignal)) return true;
  }
  fRunNumber=runnumber;
  fEvID=eventID;
  fEvIDExt=eventID_Ext;
  fEvIDLong=eventID_Long;
  fPtGen=ptgen;
  
  //topological variables
  //common
  fPt=cand->Pt();
  fY=cand->Y(4122);
  fEta=cand->Eta();
  fPhi=cand->Phi();
  fDecayLength=cand->DecayLength();
  fDecayLengthXY=cand->DecayLengthXY();
  fNormDecayLengthXY=cand->NormalizedDecayLengthXY();
  fCosP=cand->CosPointingAngle();
  fCosPXY=cand->CosPointingAngleXY();
  fImpParXY=cand->ImpParXY();
  fDCA=cand->GetDCA();
  fNormd0MeasMinusExp=ComputeMaxd0MeasMinusExp(cand,bfield);
  fSumImpParProngs=cand->Getd0Prong(0)*cand->Getd0Prong(0)+cand->Getd0Prong(1)*cand->Getd0Prong(1)+cand->Getd0Prong(2)*cand->Getd0Prong(2);
  
  //Lc -> pKpi variables
  if(masshypo==1){ //pKpi
    fInvMass=((AliAODRecoDecayHF3Prong*)cand)->InvMassLcpKpi();
  }
  else if(masshypo==2){ //piKp
    fInvMass=((AliAODRecoDecayHF3Prong*)cand)->InvMassLcpiKp();
  }
  else return false;
  fSigmaVertex=((AliAODRecoDecayHF3Prong*)cand)->GetSigmaVert();
  fDist12toPrim=((AliAODRecoDecayHF3Prong*)cand)->GetDist12toPrim();
  fDist23toPrim=((AliAODRecoDecayHF3Prong*)cand)->GetDist23toPrim();
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    fImpParProng[iProng]=cand->Getd0Prong(iProng);
    fDCAProng[iProng]=cand->GetDCA(iProng);
  }
    
  //single track variables
  AliAODTrack* prongtracks[3];
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) prongtracks[iProng] = (AliAODTrack*)cand->GetDaughter(iProng);
  bool setsingletrack = SetSingleTrackVars(prongtracks);  
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt==kNoPID) return true;

  bool setpid = SetPidVars(prongtracks,pidrespo,true,true,true,true,true);
  if(!setpid) return false;

  return true;
}

void AliHFTreeHandlerLctopKpi::AddBranchResonantDecay(TTree *t) {
  t->Branch("resonant_decay_mc",&fResonantDecayTypeMC);
}


int AliHFTreeHandlerLctopKpi::GetLcResonantDecay(TClonesArray *arrMC, AliAODMCParticle *mcPart) 
{
  int numberOfLambdac=0;
  if(TMath::Abs(mcPart->GetPdgCode())!=4122) return -1;
  Int_t nDaugh = (Int_t)mcPart->GetNDaughters();
  if(nDaugh<2) return -1;
  if(nDaugh>3) return -1;
  AliAODMCParticle* pdaugh1 = (AliAODMCParticle*)arrMC->At(mcPart->GetDaughterLabel(0));
  if(!pdaugh1) {return -1;}
  Int_t number1 = TMath::Abs(pdaugh1->GetPdgCode());
  AliAODMCParticle* pdaugh2 = (AliAODMCParticle*)arrMC->At(mcPart->GetDaughterLabel(1));
  if(!pdaugh2) {return -1;}
  Int_t number2 = TMath::Abs(pdaugh2->GetPdgCode());
  if(nDaugh==3){
    Int_t thirdDaugh=mcPart->GetDaughterLabel(1)-1;
    AliAODMCParticle* pdaugh3 = (AliAODMCParticle*)arrMC->At(thirdDaugh);
    if(!pdaugh3) return -1;
    Int_t number3 = TMath::Abs(pdaugh3->GetPdgCode());
    if((number1==321 && number2==211 && number3==2212) || (number1==211 && number2==321 && number3==2212) || (number1==211 && number2==2212 && number3==321) || (number1==321 && number2==2212 && number3==211) || (number1==2212 && number2==321 && number3==211) || (number1==2212 && number2==211 && number3==321)) numberOfLambdac = kNonResonant;
  }
  if(nDaugh==2){
    Int_t nfiglieK=0;
    if((number1==2212 && number2==313)){
      nfiglieK=pdaugh2->GetNDaughters();
      if(nfiglieK!=2) return -1;
      AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)arrMC->At(pdaugh2->GetDaughterLabel(0));
      AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)arrMC->At(pdaugh2->GetDaughterLabel(1));
      if(!pdaughK1) return -1;
      if(!pdaughK2) return -1;
      if((TMath::Abs(pdaughK1->GetPdgCode())==211 && TMath::Abs(pdaughK2->GetPdgCode())==321) || (TMath::Abs(pdaughK1->GetPdgCode())==321 && TMath::Abs(pdaughK2->GetPdgCode())==211)) numberOfLambdac = kKstar;
    }
    if((number1==313 && number2==2212)){
      nfiglieK=pdaugh1->GetNDaughters();
      if(nfiglieK!=2) return -1;
      AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)arrMC->At(pdaugh1->GetDaughterLabel(0));
      AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)arrMC->At(pdaugh1->GetDaughterLabel(1));
      if(!pdaughK1) return -1;
      if(!pdaughK2) return -1;
      if((TMath::Abs(pdaughK1->GetPdgCode())==211 && TMath::Abs(pdaughK2->GetPdgCode())==321) || (TMath::Abs(pdaughK1->GetPdgCode())==321 && TMath::Abs(pdaughK2->GetPdgCode())==211)) numberOfLambdac = kKstar;
    }
    Int_t nfiglieDelta=0;
    if(number1==321 && number2==2224){
      nfiglieDelta=pdaugh2->GetNDaughters();
      if(nfiglieDelta!=2) return -1;
      AliAODMCParticle *pdaughD1=(AliAODMCParticle*)arrMC->At(pdaugh2->GetDaughterLabel(0));
      AliAODMCParticle *pdaughD2=(AliAODMCParticle*)arrMC->At(pdaugh2->GetDaughterLabel(1));
      if(!pdaughD1) return -1;
      if(!pdaughD2) return -1;
      if((TMath::Abs(pdaughD1->GetPdgCode())==211 && TMath::Abs(pdaughD2->GetPdgCode())==2212) || (TMath::Abs(pdaughD1->GetPdgCode())==2212 && TMath::Abs(pdaughD2->GetPdgCode())==211)) numberOfLambdac = kDelta;
    }
    if(number1==2224 && number2==321){
      nfiglieDelta=pdaugh1->GetNDaughters();
      if(nfiglieDelta!=2) return -1;
      AliAODMCParticle* pdaughD1 = (AliAODMCParticle*)arrMC->At(pdaugh1->GetDaughterLabel(0));
      AliAODMCParticle* pdaughD2 = (AliAODMCParticle*)arrMC->At(pdaugh1->GetDaughterLabel(1));
      if(!pdaughD1) return -1;
      if(!pdaughD2) return -1;
      if((TMath::Abs(pdaughD1->GetPdgCode())==211 && TMath::Abs(pdaughD2->GetPdgCode())==2212) || (TMath::Abs(pdaughD1->GetPdgCode())==2212 && TMath::Abs(pdaughD2->GetPdgCode())==211)) numberOfLambdac = kDelta;
    }

    Int_t nfiglieLa=0;
    if(number1==3124 && number2==211){
      nfiglieLa=pdaugh1->GetNDaughters();
      if(nfiglieLa!=2) return -1;
      AliAODMCParticle *pdaughL1=(AliAODMCParticle*)arrMC->At(pdaugh1->GetDaughterLabel(0));
      AliAODMCParticle *pdaughL2=(AliAODMCParticle*)arrMC->At(pdaugh1->GetDaughterLabel(1));
      if(!pdaughL1) return -1;
      if(!pdaughL2) return -1;
      if((TMath::Abs(pdaughL1->GetPdgCode())==321 && TMath::Abs(pdaughL2->GetPdgCode())==2212) || (TMath::Abs(pdaughL1->GetPdgCode())==2212 && TMath::Abs(pdaughL2->GetPdgCode())==321)) numberOfLambdac = kL1520;
    }
    if(number1==211 && number2==3124){
      nfiglieLa=pdaugh2->GetNDaughters();
      if(nfiglieLa!=2) return -1;
      AliAODMCParticle *pdaughL1=(AliAODMCParticle*)arrMC->At(pdaugh2->GetDaughterLabel(0));
      AliAODMCParticle *pdaughL2=(AliAODMCParticle*)arrMC->At(pdaugh2->GetDaughterLabel(1));
      if(!pdaughL1) return -1;
      if(!pdaughL2) return -1;
      if((TMath::Abs(pdaughL1->GetPdgCode())==321 && TMath::Abs(pdaughL2->GetPdgCode())==2212) || (TMath::Abs(pdaughL1->GetPdgCode())==2212 && TMath::Abs(pdaughL2->GetPdgCode())==321)) numberOfLambdac = kL1520;

    }
  }
  return numberOfLambdac;
}
