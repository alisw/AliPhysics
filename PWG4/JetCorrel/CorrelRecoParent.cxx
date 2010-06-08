/**************************************************************************
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
/* $Id: $ */

//______________________________________________________________________________
// Container class for reconstructed parents. Reconstruction method uses
// AliKFParticle or TLorentzVector as indicated by kUseAliKF
//-- Author: Paul Constantin

#include "CorrelRecoParent.h"

using namespace std;

CorrelRecoParent_t::CorrelRecoParent_t() :
  CorrelParticle_t(), fAssym(-999.), fOpenAng(-999.), fjcESD(NULL){
  // default constructor
}

CorrelRecoParent_t& CorrelRecoParent_t::operator=(const CorrelRecoParent_t& rhs){
  fPt      = rhs.Pt()*rhs.Q();
  fPhi     = rhs.Phi();
  fEta     = rhs.Eta();
  fMass    = rhs.M();
  fID      = rhs.ID();
  fAssym   = rhs.Assym();
  fOpenAng = rhs.OpenAng();
  fjcESD     = rhs.Evt();
  return *this;
}

CorrelRecoParent_t* CorrelRecoParent_t::Copy(){
  // creates and returns deep object copy
  CorrelRecoParent_t *copy = new CorrelRecoParent_t;
  copy->fPt      = this->Pt()*this->Q();
  copy->fPhi     = this->Phi();
  copy->fEta     = this->Eta();
  copy->fMass    = this->M();
  copy->fID      = this->ID();
  copy->fAssym   = this->Assym();
  copy->fOpenAng = this->OpenAng();
  copy->fjcESD     = this->Evt();
  return copy;
}

Bool_t CorrelRecoParent_t::Reconstruct(CorrelParticle_t *p1, CorrelParticle_t *p2, Bool_t kUseAliKF){
  // main method for parent reconstruction

  if(p1->ID()==t_photon && p2->ID()==t_photon) fID = t_diphoton;
  else if(p1->ID()==t_electron && p2->ID()==t_electron) fID = t_dielectron;
  else if(p1->ID()==t_jet && p2->ID()==t_jet) fID = t_dijet;
  else fID = t_dihadron;
  
  if(fID==t_dielectron && kUseAliKF){
    // code for parent reconstruction based on AliKFParticle:
    if(!fjcESD)
      {std::cerr<<"CorrelRecoParent_t::Reconstruct - undefined event"<<std::endl; exit(-1);}
    CorrelKFTrack_t* e1 = dynamic_cast<CorrelKFTrack_t*>(p1);
    CorrelKFTrack_t* e2 = dynamic_cast<CorrelKFTrack_t*>(p2);
    if(!e1 || !e2)
      {std::cerr<<"CorrelRecoParent_t::Reconstruct - failed particle casting"<<std::endl; exit(-1);}

    AliKFVertex primVtx(*(fjcESD->GetPrimaryVertex()));
    if(primVtx.GetNContributors()<=0) return kFALSE;
    AliKFParticle::SetField(fjcESD->GetMagneticField());
    AliKFParticle* pKF1 = new AliKFParticle; AliKFParticle* pKF2 = new AliKFParticle;
    if(e1->Param()==NULL || e1->Covar()==NULL || e2->Param()==NULL || e2->Covar()==NULL)
      {std::cerr<<"CorrelRecoParent_t::Reconstruct - null parameter pointer"<<std::endl; exit(-1);}      
    pKF1->Create(e1->Param(), e1->Covar(), e1->Q(), 11*e1->Q());  // AliKFParticle uses Pythia PIDs
    pKF2->Create(e2->Param(), e2->Covar(), e2->Q(), 11*e2->Q());  // so electron=11 and positron=-11

    AliKFParticle DiElectron(*pKF1,*pKF2);
    DiElectron.SetMassConstraint(91,3); // the dielectron mass cut (1st=mass,2nd=sigma)
    primVtx += DiElectron;
    DiElectron.SetProductionVertex(primVtx);
    Double_t chi2=100000.;
    if(DiElectron.GetNDF()!=0) chi2=DiElectron.GetChi2()/DiElectron.GetNDF();
    if(chi2>3.5) return kFALSE; // the dielectron vertex cut

    Float_t px = DiElectron.GetPx();
    Float_t py = DiElectron.GetPy();
    Float_t pz = DiElectron.GetPz();
    Float_t ener = DiElectron.GetE();
    Float_t theta = TMath::ATan(TMath::Sqrt(px*px+py*py)/pz);
    fPt = TMath::Sqrt(px*px+py*py)*Float_t(p1->Q()*p2->Q());
    fPhi = TMath::ATan(py/px);
    fEta = -TMath::Log(TMath::Tan(theta/2.));
    fMass = ener - TMath::Sqrt(px*px+py*py+pz*pz);
    fAssym = TMath::Abs((e1->P()-e2->P())/(e1->P()+e2->P()));
    fOpenAng = -999.; // compute opening angle

    delete pKF1; delete pKF2;
    
  } else {

    // code for parent reconstruction based on TLorentzVector:
    TLorentzVector part1, part2, fPartSum;
    part1.SetPxPyPzE(p1->Px(), p1->Py(), p1->Pz(), p1->E());
    part2.SetPxPyPzE(p2->Px(), p2->Py(), p2->Pz(), p2->E());

    fPartSum.SetPxPyPzE(0.,0.,0.,0.);
    fPartSum = part1 + part2;
    Float_t pT = fPartSum.Pt();
    if(TMath::Abs(p1->Q())>0 && TMath::Abs(p2->Q())>0)
      fPt = pT*Float_t(p1->Q()*p2->Q()); // fPt stores charge as its sign
    else
      fPt = pT;
    fPhi = fPartSum.Phi();
    fEta = fPartSum.Eta();
    fMass = fPartSum.M();
    fAssym = TMath::Abs((p1->P()-p2->P())/(p1->P()+p2->P()));
    fOpenAng = part1.Angle(part2.Vect());
    if(NotInMass(fID,fMass)) return kFALSE;
  }
  
  return kTRUE;
}

Bool_t CorrelRecoParent_t::NotInMass(cPartType_t pID, Float_t mass){
  // THE MASS RANGES SHOULD PROBABLY BE MOMENTUM AND CENTRALITY DEPENDENT!!!
  if(pID!=t_dielectron && pID!=t_diphoton) return kFALSE;

  const Float_t kZ0MassMin = 82;    // 91-3*3
  const Float_t kZ0MassMax = 100;   // 91+3*3
  const Float_t kPi0MassMin = 0.11; // 0.14-3*0.1
  const Float_t kPi0MassMax = 0.17; // 0.14+3*0.1
  if(pID==t_dielectron && mass>kZ0MassMin && mass<kZ0MassMax) return kFALSE;
  if(pID==t_diphoton && mass>kPi0MassMin && mass<kPi0MassMax) return kFALSE;

  return kTRUE;
}

void CorrelRecoParent_t::Show() const {
  // printout method
  std::cout<<" RecoParent pT="<<Pt()<<" phi="<<Phi()<<" eta="<<Eta()<<" m="<<M()<<" id="<<ID()
	   <<" assym="<<Assym()<<" open="<<OpenAng()<<std::endl;
}
