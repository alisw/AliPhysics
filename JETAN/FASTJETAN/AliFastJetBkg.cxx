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

/* $Id$ */

//---------------------------------------------------------------------
// Class to calculate the background per unit area
// manages the search for jets
// Authors: Elena Bruna elena.bruna@yale.edu
//          Sevil Salur ssalur@lbl.gov
//
// 2011 :
// renamed from AliJetBkg to AliFastJetBkg as this class uses only FASTJET based algos.
//---------------------------------------------------------------------

#include <Riostream.h> 
#include <TClonesArray.h>
#include <TF1.h>
#include <TString.h>

#include "AliJetHeader.h"
#include "AliFastJetHeaderV1.h"
#include "AliAODJet.h"
#include "AliFastJetInput.h"
#include "AliFastJetBkg.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"

#include<vector> 

using namespace std;

ClassImp(AliFastJetBkg)

////////////////////////////////////////////////////////////////////////

AliFastJetBkg::AliFastJetBkg():
  TObject(),
  fHeader(0),
  fInputFJ(0)
{
  // Default constructor
}

//______________________________________________________________________
AliFastJetBkg::AliFastJetBkg(const AliFastJetBkg& input):
  TObject(input),
  fHeader(input.fHeader),
  fInputFJ(input.fInputFJ)
{
  // copy constructor
}

//______________________________________________________________________
AliFastJetBkg& AliFastJetBkg::operator=(const AliFastJetBkg& source)
{
  // Assignment operator. 
  if(this!=&source){ 
   TObject::operator=(source);
   fHeader = source.fHeader;
   fInputFJ = source.fInputFJ;
  }
  
  return *this;

}

//___________________________________________________________________
void AliFastJetBkg::BkgFastJetb(Double_t& rho,Double_t& sigma, 
				Double_t& meanarea)
{
  // Bkg estimation
   
  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader; 
  Int_t debug  = header->GetDebug();     // debug option
  if(debug>0) cout<<"===============  AliFastJetBkg::BkgFastJetb()  =========== "<<endl;
  vector<fastjet::PseudoJet> inputParticles=fInputFJ->GetInputParticles();
  
  double rParamBkg = header->GetRparamBkg(); //Radius for background calculation

  Double_t medianb,sigmab,meanareab;
  CalcRhob(medianb,sigmab,meanareab,inputParticles,rParamBkg,"All");
  rho=medianb;
  sigma=sigmab;
  meanarea=meanareab;
 
}

//_________________________________________________________________
void AliFastJetBkg::BkgFastJetWoHardest(Double_t& rho,Double_t& sigma, 
					Double_t& meanarea)
{

  // Bkg estimation without hardest jet

  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader; 
  Int_t debug  = header->GetDebug();     // debug option
  if(debug) cout<<"===============  AliFastJetBkg::BkgWoHardest()  =========== "<<endl;
  vector<fastjet::PseudoJet> inputParticles=fInputFJ->GetInputParticles();
  
  double rParamBkg = header->GetRparamBkg(); //Radius for background calculation  
  Double_t medianb,sigmab,meanareab;
  CalcRhoWoHardest(medianb,sigmab,meanareab,inputParticles,rParamBkg,"All");
  rho=medianb;
  sigma=sigmab;
  meanarea=meanareab;

}

//____________________________________________________________________
void AliFastJetBkg::CalcRhob(Double_t& median,Double_t& 
			     sigma,Double_t& 
			     meanarea,vector<fastjet::PseudoJet> inputParticles,Double_t 
			     rParamBkg,TString method)
{
  // calculate rho using the fastjet method

  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Int_t debug  = header->GetDebug();     // debug option

  fastjet::Strategy strategy = header->GetStrategy();
  fastjet::RecombinationScheme recombScheme = header->GetRecombScheme();
  fastjet::JetAlgorithm algorithm = header->GetBGAlgorithm();
  fastjet::JetDefinition jetDef(algorithm, rParamBkg, recombScheme, strategy);

  // create an object that specifies how we to define the area
  fastjet::AreaDefinition areaDef;
  double ghostEtamax = header->GetGhostEtaMax(); 
  double ghostArea   = header->GetGhostArea(); 
  int    activeAreaRepeats = header->GetActiveAreaRepeats(); 

  // now create the object that holds info about ghosts

  fastjet::GhostedAreaSpec ghost_spec(ghostEtamax, activeAreaRepeats, ghostArea);
  // and from that get an area definition
  fastjet::AreaType areaType = header->GetAreaType();
  areaDef = fastjet::AreaDefinition(areaType,ghost_spec);
  
  //fastjet::ClusterSequenceArea clust_seq(inputParticles, jetDef);
  fastjet::ClusterSequenceArea clust_seq(inputParticles, jetDef,areaDef);
  TString comment = "Running FastJet algorithm for BKG calculation with the following setup. ";
  comment+= "Jet definition: ";
  comment+= TString(jetDef.description());
  // comment+= ". Area definition: ";
  // comment+= TString(areaDef.description());
  comment+= ". Strategy adopted by FastJet: ";
  comment+= TString(clust_seq.strategy_string());
  comment+= Form("Method: %s",method.Data());
  header->SetComment(comment);
  if(debug>0){
    cout << "--------------------------------------------------------" << endl;
    cout << comment << endl;
    cout << "--------------------------------------------------------" << endl;
  }

  vector<fastjet::PseudoJet> inclusiveJets = clust_seq.inclusive_jets(0.);
  vector<fastjet::PseudoJet> jets = sorted_by_pt(inclusiveJets); 

  double phiMin = 0, phiMax = 0, rapMin = 0, rapMax = 0;

  phiMin = 0;
  phiMax = 2*TMath::Pi();
 
  rapMax = ghostEtamax - rParamBkg;
  rapMin = - ghostEtamax + rParamBkg;

  fastjet::RangeDefinition range(rapMin, rapMax, phiMin, phiMax);

  double medianb, sigmab, meanareab;
  clust_seq.get_median_rho_and_sigma(inclusiveJets, range, false, medianb, sigmab, meanareab, false);
  median=medianb;
  sigma=sigmab;
  meanarea=meanareab; 
 
}

//____________________________________________________________________
void AliFastJetBkg::CalcRhoWoHardest(Double_t& median,Double_t& 
				     sigma,Double_t& meanarea,vector<fastjet::PseudoJet> inputParticles,Double_t 
				     rParamBkg,TString method)
{
  // calculate rho (without the hardest jet) using the fastjet method

  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Int_t debug  = header->GetDebug();     // debug option

  fastjet::Strategy strategy = header->GetStrategy();
  fastjet::RecombinationScheme recombScheme = header->GetRecombScheme();
  fastjet::JetAlgorithm algorithm = header->GetBGAlgorithm();
  fastjet::JetDefinition jetDef(algorithm, rParamBkg, recombScheme, strategy);

  // create an object that specifies how we to define the area
  fastjet::AreaDefinition areaDef;
  double ghostEtamax = header->GetGhostEtaMax(); 
  double ghostArea   = header->GetGhostArea(); 
  int    activeAreaRepeats = header->GetActiveAreaRepeats(); 

  // now create the object that holds info about ghosts

  fastjet::GhostedAreaSpec ghost_spec(ghostEtamax, activeAreaRepeats, ghostArea);
  // and from that get an area definition
  fastjet::AreaType areaType = header->GetAreaType();
  areaDef = fastjet::AreaDefinition(areaType,ghost_spec);
  //fastjet::ClusterSequenceArea clust_seq(inputParticles, jetDef);
  fastjet::ClusterSequenceArea clust_seq(inputParticles,jetDef,areaDef);
  TString comment = "Running FastJet algorithm for BKG calculation with the following setup. ";
  comment+= "Jet definition: ";
  comment+= TString(jetDef.description());
  // comment+= ". Area definition: ";
  // comment+= TString(areaDef.description());
  comment+= ". Strategy adopted by FastJet: ";
  comment+= TString(clust_seq.strategy_string());
  comment+= Form("Method: %s",method.Data());
  header->SetComment(comment);
if(debug>0){
    cout << "--------------------------------------------------------" << endl;
    cout << comment << endl;
    cout << "--------------------------------------------------------" << endl;
  }

  vector<fastjet::PseudoJet> inclusiveJets = clust_seq.inclusive_jets(0.);
  vector<fastjet::PseudoJet> jets = sorted_by_pt(inclusiveJets); 
  vector<fastjet::PseudoJet> jets2=sorted_by_pt(inclusiveJets);
  if(jets2.size()>=2) jets2.erase(jets2.begin(),jets2.begin()+1);
    
  double phiMin = 0, phiMax = 0, rapMin = 0, rapMax = 0;

  phiMin = 0;
  phiMax = 2*TMath::Pi();
 
  rapMax = ghostEtamax - rParamBkg;
  rapMin = - ghostEtamax + rParamBkg;

  fastjet::RangeDefinition range(rapMin, rapMax, phiMin, phiMax);

  double medianb, sigmab, meanareab;
  clust_seq.get_median_rho_and_sigma(jets2, range, false, medianb, sigmab, 
				     meanareab, false);
  median=medianb;
  sigma=sigmab;
  meanarea=meanareab; 

}

//___________________________________________________________________
Float_t AliFastJetBkg::BkgFastJet()
{
  // Return background  
 
  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Int_t debug  = header->GetDebug();     // debug option

  if(debug>0) cout<<"===============  AliFastJetBkg::BkgFastJet()  =========== "<<endl;
  vector<fastjet::PseudoJet> inputParticles=fInputFJ->GetInputParticles();
  
  if(debug>0) cout<<"printing inputParticles for BKG "<<inputParticles.size()<<endl;
  
  double rParamBkg = header->GetRparamBkg(); //Radius for background calculation
  Double_t rho=CalcRho(inputParticles,rParamBkg,"All");
  if(debug) cout<<"-------- rho (from all part)="<<rho<<endl; 
  return rho;
 
}

//___________________________________________________________________
Float_t AliFastJetBkg::BkgChargedFastJet()
{
  // Background for charged jets
  
  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Int_t debug  = header->GetDebug();     // debug option

  if(debug>0) cout<<"===============  AliFastJetBkg::BkgChargedFastJet()  =========== "<<endl;

  vector<fastjet::PseudoJet> inputParticlesCharged=fInputFJ->GetInputParticlesCh();
  
  if(debug>0) cout<<"printing CHARGED inputParticles for BKG "<<inputParticlesCharged.size()<<endl;

  double rParam = header->GetRparam();

  Double_t rho=CalcRho(inputParticlesCharged,rParam,"Charg");

  if(debug>0) cout<<"-------- rho (from CHARGED part)="<<rho<<endl; 
  return rho;

}

//___________________________________________________________________
Float_t AliFastJetBkg::BkgStat()
{
  // background subtraction using statistical method
 
  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Int_t debug  = header->GetDebug();     // debug option

  if(debug>0) cout<<"==============AliFastJetBkg::BkgStat()============="<<endl;
  //TO BE IMPLEMENTED 
  Int_t nTracks= 0;
  TF1 fun("fun",BkgFunction,0,800,1);
  Double_t enTot=fun.Eval(nTracks);
  Double_t accEMCal=2*0.7*110./180*TMath::Pi();//2.68 area of EMCal
  return enTot/accEMCal;

}

//___________________________________________________________________
Float_t AliFastJetBkg::BkgFastJetCone(TClonesArray* fAODJets)
{
  // Cone background subtraction method applied on the fastjet: REmove the particles of the
  // two largest jets with the given R from the estimation of new rho. 

  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Int_t debug  = header->GetDebug();     // debug option

  if(debug>0) cout<<"==============AliFastJetBkg::SubtractFastJetBackgCone()============="<<endl;

  Float_t rc= header->GetRparam();

  //Hard wired Calorimeter area (get it later from the AliJetReaderHeader.h)
  Double_t accEMCal=2*0.7*110./180*TMath::Pi();//2.68 area of EMCal

  Int_t nJ=fAODJets->GetEntries(); //this must be the # of jets... 
  if(debug>0) cout<<"nJets:  "<<nJ<<endl;
  
  // Information extracted from fInputParticle
  // load input vectors and calculate total energy in array
  Float_t pt,eta,phi;
  Float_t jeteta = 0,jetphi = 0,jeteta1 = 0, jetphi1 = 0;
  Float_t rhoback=0.0;

  Float_t ptallback=0.0; //particles without the jet
  Float_t restarea=accEMCal; //initial area set 
  Bool_t acc=0;
  Bool_t acc1=0;
  Float_t rCone=0.4;
  
  if(nJ==1) { 
    AliAODJet *jettmp = (AliAODJet*)(fAODJets->At(0));
    jeteta=jettmp->Eta();
    jetphi=jettmp->Phi();
    acc=EmcalAcceptance(jeteta,jetphi,rCone);
    if(acc==1)restarea= accEMCal-TMath::Pi()*rc*rc;
    if(debug) cout<<" acc  "<<acc<<endl;
  }
  
  if(nJ>=2) { 
    AliAODJet *jettmp = (AliAODJet*)(fAODJets->At(0));
    AliAODJet *jettmp1 = (AliAODJet*)(fAODJets->At(1));
    jeteta=jettmp->Eta();
    jetphi=jettmp->Phi();
    jeteta1=jettmp1->Eta();
    jetphi1=jettmp1->Phi(); 
    acc=EmcalAcceptance(jeteta,jetphi,rCone);
    acc1=EmcalAcceptance(jeteta1,jetphi1,rCone);
    if(acc1==1 && acc==1)restarea= accEMCal-2*TMath::Pi()*rc*rc;
    if(acc1==1 && acc==0)restarea= accEMCal-TMath::Pi()*rc*rc;
    if(acc1==0 && acc==1)restarea= accEMCal-TMath::Pi()*rc*rc;

    if(debug) cout<<" acc1="<<acc<<"  acc2="<<acc1<<"  restarea="<<restarea<<endl;

  }
  
  // cout<<" nIn = "<<nIn<<endl;
  Float_t sumpt=0;
  vector<fastjet::PseudoJet> inputParticles=fInputFJ->GetInputParticles();
  for(UInt_t i=0; i<inputParticles.size(); i++)
    { // Loop over input list of particles
      pt    = inputParticles[i].perp();
      eta   = inputParticles[i].eta();
      phi   = inputParticles[i].phi();

      // To be updated
      //cout<<"test emcal acceptance for particles "<<EmcalAcceptance(eta,phi,0.)<<endl;
	
      Float_t deta=0.0, dphi=0.0, dr=100.0;
      Float_t deta1=0.0, dphi1=0.0, dr1=100.0;

      //cout<<i<<"  pt="<<pt<<"  eta="<<eta<<"  phi="<<phi<<endl;
      if(phi>1.396 && phi<3.316 && eta>-0.7 && eta<0.7){
	sumpt+=pt;
	//if(i<30)cout<<i<<" pt = "<<pt<<endl;

	if(nJ==1 && acc==1) { 
	  deta = eta - jeteta;
	  dphi = phi - jetphi;
	  if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	  if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	  dr = TMath::Sqrt(deta * deta + dphi * dphi);
	  if(dr<=rc)sumpt-=pt;
	}
	
	if(nJ>=2) { 
	  if(acc==1){
	    deta = eta - jeteta;
	    dphi = phi - jetphi;
	    if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	    if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	    dr = TMath::Sqrt(deta * deta + dphi * dphi);
	    if(dr<=rc)sumpt-=pt;
	  }
	  if(acc1==1){
	    deta1 = eta - jeteta1;
	    dphi1 = phi - jetphi1;
	    if (dphi1 < -TMath::Pi()) dphi1= -dphi1 - 2.0 * TMath::Pi();
	    if (dphi1 > TMath::Pi()) dphi1 = 2.0 * TMath::Pi() - dphi1;
	    dr1 = TMath::Sqrt(deta1 * deta1 + dphi1 * dphi1);
	    if(dr1<=rc)sumpt-=pt;
	  }
	}

	if(dr >= rc && dr1 >=rc) { 
	  // particles outside both cones
	  if(debug>1) cout<<" out of the cone  "<<dr<<"    "<<deta<<"  deltaeta  "<<dphi<<"  dphi "<<i<<"  particle  "<<endl;
	  if(debug>1) cout<<" out of the cone  "<<dr1<<"    "<<deta1<<"  deltaeta1  "<<dphi1<<"  dphi1 "<<i<<"  particle  "<<endl;
	  ptallback+=pt;
	}
      }
    } // End loop on input list of particles 
  
  if(debug>0) cout<<"total area left "<<restarea<<endl;
  if(debug>0) cout<<"sumpt="<<sumpt<<endl;
  // if(acc==1 || acc1==1) rhoback= ptallback/restarea;
  // else rhoback=ptallback;

  rhoback= ptallback/restarea;
  if(debug)cout<<"rhoback    "<<rhoback<<"     "<<nJ<<"   "<<endl;

  return rhoback;
   
}

//___________________________________________________________________
Double_t AliFastJetBkg::CalcRho(vector<fastjet::PseudoJet> inputParticles,Double_t rParamBkg,TString method)
{
  // calculate rho using the fastjet method

  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Int_t debug  = header->GetDebug();     // debug option

  fastjet::Strategy strategy = header->GetStrategy();
  fastjet::RecombinationScheme recombScheme = header->GetRecombScheme();
  fastjet::JetAlgorithm algorithm = header->GetBGAlgorithm();
  fastjet::JetDefinition jetDef(algorithm, rParamBkg, recombScheme, strategy);

  // create an object that specifies how we to define the area
  fastjet::AreaDefinition areaDef;
  double ghostEtamax       = header->GetGhostEtaMax(); 
  double ghostArea         = header->GetGhostArea(); 
  int    activeAreaRepeats = header->GetActiveAreaRepeats(); 
  
  // now create the object that holds info about ghosts

  if (method.Contains("Charg"))ghostEtamax=0.9;

  fastjet::GhostedAreaSpec ghost_spec(ghostEtamax, activeAreaRepeats, ghostArea);
  // and from that get an area definition
  fastjet::AreaType areaType = header->GetAreaType();
  areaDef = fastjet::AreaDefinition(areaType,ghost_spec);
  if(debug>0) cout<<"rParamBkg="<<rParamBkg<<"  ghostEtamax="<<ghostEtamax<<"  ghostArea="<<ghostArea<<" areadef="<<TString(areaDef.description())<< endl;
  //fastjet::ClusterSequenceArea clust_seq(inputParticles, jetDef);
  fastjet::ClusterSequenceArea clust_seq(inputParticles, jetDef,areaDef);
  TString comment = "Running FastJet algorithm for BKG calculation with the following setup. ";
  comment+= "Jet definition: ";
  comment+= TString(jetDef.description());
  // comment+= ". Area definition: ";
  // comment+= TString(areaDef.description());
  comment+= ". Strategy adopted by FastJet: ";
  comment+= TString(clust_seq.strategy_string());
  header->SetComment(comment);
  if(debug>0){
    cout << "--------------------------------------------------------" << endl;
    cout << comment << endl;
    cout << "--------------------------------------------------------" << endl;
  }

  double ptmin = header->GetPtMin(); 
  vector<fastjet::PseudoJet> inclusiveJets = clust_seq.inclusive_jets(ptmin);
  vector<fastjet::PseudoJet> jets = sorted_by_pt(inclusiveJets); 

  if (debug>0) {
    cout<<"# of BKG jets = "<<jets.size()<<endl;
    for (size_t j = 0; j < jets.size(); j++) { // loop for jets   
      printf("BKG Jet found %5d %9.5f %8.5f %10.3f %4.4f \n",(Int_t)j,jets[j].rap(),jets[j].phi(),jets[j].perp(),clust_seq.area(jets[j]));
    }
  }
  
  double phiMin = 0, phiMax = 0, rapMin = 0, rapMax = 0;

  if (method.Contains("All")){
    phiMin = 80.*TMath::Pi()/180+rParamBkg;
    phiMax = 190.*TMath::Pi()/180-rParamBkg;
  }
  if (method.Contains("Charg")){
    phiMin = 0;
    phiMax = 2*TMath::Pi();
  }
  rapMax = ghostEtamax - rParamBkg;
  rapMin = - ghostEtamax + rParamBkg;

  fastjet::RangeDefinition range(rapMin, rapMax, phiMin, phiMax);

  Double_t rho=clust_seq.median_pt_per_unit_area(range);
  // double median, sigma, meanArea;
  // clust_seq.get_median_rho_and_sigma(inclusiveJets, range, false, median, sigma, meanArea, true);
  // fastjet::ActiveAreaSpec area_spec(ghostEtamax,activeAreaRepeats,ghostArea);

  // fastjet::ClusterSequenceActiveArea clust_seq_bkg(inputParticles, jetDef,area_spec);

  if(debug>0) cout<<"bkg in R="<<rParamBkg<<"  : "<<rho<<" range: Rap="<<rapMin<<","<<rapMax<<" --  phi="<<phiMin<<","<<phiMax<<endl;

  return rho;

}

//___________________________________________________________________
Double_t  AliFastJetBkg::BkgFunction(Double_t */*x*/,Double_t */*par*/)
{
  // to be implemented--- (pT + Energy in EMCal Acceptance vs Multiplicity)
  return 1;

}

//___________________________________________________________________
Bool_t AliFastJetBkg::EmcalAcceptance(Float_t eta, Float_t phi, Float_t radius) const
{
  // Apply emcal acceptance cuts
  // To be updated

  Float_t meanPhi=190./180.*TMath::Pi()-110./180.*TMath::Pi()/2;
  Float_t deltaphi=110./180.*TMath::Pi();
  Float_t phicut=deltaphi/2.-radius;
  Float_t etacut=0.7-radius;
  //cout<<"  eta    "<<eta<<"  phi    "<<phi<<endl;
  //cout<<etacut<<"    "<<phicut<<"    "<<meanPhi<<"    "<<endl;
  if(TMath::Abs(eta)<etacut && TMath::Abs(phi-meanPhi)<phicut) return 1;
  else return 0; 

}
