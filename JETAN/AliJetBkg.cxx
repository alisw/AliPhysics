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

#include <Riostream.h> 
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TArrayF.h>
#include <TClonesArray.h>
#include <TF1.h>
#include <TString.h>

#include "AliJetHeader.h"
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliFastJetFinder.h"
#include "AliFastJetHeaderV1.h"
#include "AliJetReaderHeader.h"
#include "AliJetReader.h"
#include "AliJetUnitArray.h"
#include "AliFastJetInput.h"
#include "AliESDEvent.h"


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
// get info on how fastjet was configured
#include "fastjet/config.h"

#ifdef ENABLE_PLUGIN_SISCONE
#include "fastjet/SISConePlugin.hh"
#endif

#include<sstream>  // needed for internal io
#include<vector> 
#include <cmath> 

using namespace std;
#include "AliJetBkg.h"

ClassImp(AliJetBkg)

////////////////////////////////////////////////////////////////////////

AliJetBkg::AliJetBkg():
  TObject(),
  fReader(0),
  fHeader(0),
  fInputFJ(0)
{
  // Default constructor

}
//______________________________________________________________________
AliJetBkg::AliJetBkg(const AliJetBkg& input):
  TObject(input),
    fReader(input.fReader),
    fHeader(input.fHeader),
    fInputFJ(input.fInputFJ)
{
  // copy constructor

}
//______________________________________________________________________
AliJetBkg& AliJetBkg::operator=(const AliJetBkg& source){
    // Assignment operator. 
    this->~AliJetBkg();
    new(this) AliJetBkg(source);
    return *this;

}


//___________________________________________________________________
  void AliJetBkg::BkgFastJetb(Double_t& rho,Double_t& sigma, 
Double_t& meanarea){

    AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader; 
    Bool_t debug  = header->GetDebug();     // debug option 
    if(debug)cout<<"===============  AliJetBkg::BkgFastJetb()  =========== "<<endl;
  vector<fastjet::PseudoJet> inputParticles=fInputFJ->GetInputParticles();
  
  //cout<<"printing inputParticles for BKG "<<inputParticles.size()<<endl;
  

   

  double rParamBkg = header->GetRparamBkg(); //Radius for background calculation

  Double_t medianb,sigmab,meanareab;
  CalcRhob(medianb,sigmab,meanareab,inputParticles,rParamBkg,"All");
  rho=medianb;
  sigma=sigmab;
  meanarea=meanareab;
 
}
//_________________________________________________________________

  void AliJetBkg::BkgFastJetWoHardest(Double_t& rho,Double_t& sigma, 
Double_t& meanarea){
    AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader; 
    Bool_t debug  = header->GetDebug();     // debug option 
    if(debug)cout<<"===============  AliJetBkg::BkgWoHardest()  =========== "<<endl;
  vector<fastjet::PseudoJet> inputParticles=fInputFJ->GetInputParticles();
  
  //cout<<"printing inputParticles for BKG "<<inputParticles.size()<<endl;   
  double rParamBkg = header->GetRparamBkg(); //Radius for background calculation  
  Double_t medianb,sigmab,meanareab;
  CalcRhoWoHardest(medianb,sigmab,meanareab,inputParticles,rParamBkg,"All");
  rho=medianb;
  sigma=sigmab;
  meanarea=meanareab;

  }

//____________________________________________________________________

 void AliJetBkg::CalcRhob(Double_t& median,Double_t& 
sigma,Double_t& 
meanarea,vector<fastjet::PseudoJet> inputParticles,Double_t 
rParamBkg,TString method){
  //calculate rho using the fastjet method

  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Bool_t debug  = header->GetDebug();     // debug option

  fastjet::Strategy strategy = header->GetStrategy();
  fastjet::RecombinationScheme recombScheme = header->GetRecombScheme();
  fastjet::JetAlgorithm algorithm = header->GetAlgorithm(); 
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
  //  comment+= ". Area definition: ";
  //comment+= TString(areaDef.description());
  comment+= ". Strategy adopted by FastJet: ";
  comment+= TString(clust_seq.strategy_string());
  header->SetComment(comment);
  if(debug){
    cout << "--------------------------------------------------------" << endl;
    cout << comment << endl;
    cout << "--------------------------------------------------------" << endl;
  }

  vector<fastjet::PseudoJet> inclusiveJets = clust_seq.inclusive_jets(0.);
  vector<fastjet::PseudoJet> jets = sorted_by_pt(inclusiveJets); 

 
  //cout<<"# of BKG jets = "<<jets.size()<<endl;
      //  for (size_t j = 0; j < jets.size(); j++) { // loop for jets   
      //printf("BKG Jet found %5d %9.5f %8.5f %10.3f %4.4f 
      //\n",(Int_t)j,jets[j].rap(),jets[j].phi(),jets[j].perp(),
      //clust_seq.area(jets[j]));
      // }
  
  
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


 void AliJetBkg::CalcRhoWoHardest(Double_t& median,Double_t& 
sigma,Double_t& 
meanarea,vector<fastjet::PseudoJet> inputParticles,Double_t 
rParamBkg,TString method){
  //calculate rho using the fastjet method

  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Bool_t debug  = header->GetDebug();     // debug option

  fastjet::Strategy strategy = header->GetStrategy();
  fastjet::RecombinationScheme recombScheme = header->GetRecombScheme();
  fastjet::JetAlgorithm algorithm = header->GetAlgorithm(); 
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
  //cout<<"rParamBkg="<<rParamBkg<<"  ghostEtamax="<<ghostEtamax<<"  ghostArea="<<ghostArea<<" areadef="<<TString(areaDef.description())<<endl;
  //fastjet::ClusterSequenceArea clust_seq(inputParticles, jetDef);
  fastjet::ClusterSequenceArea clust_seq(inputParticles,jetDef,areaDef);
  TString comment = "Running FastJet algorithm for BKG calculation with the following setup. ";
  comment+= "Jet definition: ";
  comment+= TString(jetDef.description());
  //  comment+= ". Area definition: ";
  //comment+= TString(areaDef.description());
  comment+= ". Strategy adopted by FastJet: ";
  comment+= TString(clust_seq.strategy_string());
  header->SetComment(comment);
  if(debug){
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
Float_t AliJetBkg::BkgFastJet(){
  
  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Bool_t debug  = header->GetDebug();     // debug option

  if(debug)cout<<"===============  AliJetBkg::BkgFastJet()  =========== "<<endl;
  vector<fastjet::PseudoJet> inputParticles=fInputFJ->GetInputParticles();
  
  if(debug)cout<<"printing inputParticles for BKG "<<inputParticles.size()<<endl;
  
  for(UInt_t i=0;i<inputParticles.size();i++){
    //  cout<<"                "<<inputParticles[i].px()<<" "<<inputParticles[i].py()<<" "<<inputParticles[i].pz()<<endl;
    
  }
   
  double rParamBkg = header->GetRparamBkg(); //Radius for background calculation
  Double_t rho=CalcRho(inputParticles,rParamBkg,"All");
  if(debug)cout<<"-------- rho (from all part)="<<rho<<endl; 
  return rho;
 
}
//___________________________________________________________________
Float_t  AliJetBkg::BkgChargedFastJet(){

  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Bool_t debug  = header->GetDebug();     // debug option

  if(debug)cout<<"===============  AliJetBkg::BkgChargedFastJet()  =========== "<<endl;

  vector<fastjet::PseudoJet> inputParticlesCharged=fInputFJ->GetInputParticlesCh();
  
  if(debug)cout<<"printing CHARGED inputParticles for BKG "<<inputParticlesCharged.size()<<endl;

  for(UInt_t i=0;i<inputParticlesCharged.size();i++){
    // cout<<"                "<<inputParticlesCharged[i].px()<<" "<<inputParticlesCharged[i].py()<<" "<<inputParticlesCharged[i].pz()<<endl;

  } 

  double rParam = header->GetRparam();

  Double_t rho=CalcRho(inputParticlesCharged,rParam,"Charg");

  if(debug)cout<<"-------- rho (from CHARGED part)="<<rho<<endl; 
  return rho;
}



Float_t AliJetBkg::BkgStat()
{
  //background subtraction using statistical method
 
  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Bool_t debug  = header->GetDebug();     // debug option

  if(debug)cout<<"==============AliJetBkg::BkgStat()============="<<endl;
  //TO BE IMPLEMENTED 
  // Int_t nTracks= fReader->GetESD()->GetNumberOfTracks();
  Int_t nTracks= 0;
  TF1 fun("fun",BkgFunction,0,800,1);
  Double_t enTot=fun.Eval(nTracks);
  Double_t accEMCal=2*0.7*110./180*TMath::Pi();//2.68 area of EMCal
  return enTot/accEMCal;

}

////////////////////////////////////////////////////////////////////////
Float_t AliJetBkg::BkgFastJetCone(TClonesArray* fAODJets)
{

  // Cone background subtraction method applied on the fastjet: REmove the particles of the
  // two largest jets with the given R from the estimation of new rho. 

  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Bool_t debug  = header->GetDebug();     // debug option

  if(debug)cout<<"==============AliJetBkg::SubtractFastJetBackgCone()============="<<endl;

  Float_t rc= header->GetRparam();

  //Hard wired Calorimeter area (get it later from the AliJetReaderHeader.h)
  Double_t accEMCal=2*0.7*110./180*TMath::Pi();//2.68 area of EMCal

  Int_t nJ=fAODJets->GetEntries(); //this must be the # of jets... 
  if(debug)cout<<"nJets:  "<<nJ<<endl;
  

  //begin unit array 
  TClonesArray* fUnit = fReader->GetUnitArray();
  if(fUnit == 0) { cout << "Could not get the momentum array" << endl; return -99; }

  Int_t nIn = fUnit->GetEntries();
  if(nIn == 0) { cout << "entries = 0 ; Event empty !!!" << endl ; return -99; }
  
  // Information extracted from fUnitArray
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
    AliAODJet *jettmp = dynamic_cast<AliAODJet*>(fAODJets->At(0));
    jeteta=jettmp->Eta();
    jetphi=jettmp->Phi();
    acc=EmcalAcceptance(jeteta,jetphi,rCone);
    if(acc==1)restarea= accEMCal-TMath::Pi()*rc*rc;
    if(debug)cout<<" acc  "<<acc<<endl;
  }

  
  if(nJ>=2) { 
    AliAODJet *jettmp = dynamic_cast<AliAODJet*>(fAODJets->At(0));
    AliAODJet *jettmp1 = dynamic_cast<AliAODJet*>(fAODJets->At(1));
    jeteta=jettmp->Eta();
    jetphi=jettmp->Phi();
    jeteta1=jettmp1->Eta();
    jetphi1=jettmp1->Phi(); 
    acc=EmcalAcceptance(jeteta,jetphi,rCone);
    acc1=EmcalAcceptance(jeteta1,jetphi1,rCone);
    if(acc1==1 && acc==1)restarea= accEMCal-2*TMath::Pi()*rc*rc;
    if(acc1==1 && acc==0)restarea= accEMCal-TMath::Pi()*rc*rc;
    if(acc1==0 && acc==1)restarea= accEMCal-TMath::Pi()*rc*rc;

    if(debug)cout<<" acc1="<<acc<<"  acc2="<<acc1<<"  restarea="<<restarea<<endl;

  }
  
  // cout<<" nIn = "<<nIn<<endl;
  Float_t sumpt=0;
  for(Int_t i=0; i<nIn; i++) 
    { //Unit Array Loop
      AliJetUnitArray *uArray = (AliJetUnitArray*)fUnit->At(i);
      if(uArray->GetUnitEnergy()>0.){
	
	pt    = uArray->GetUnitEnergy();
       	eta   = uArray->GetUnitEta();
	phi   = uArray->GetUnitPhi();

	//cout<<"test emcal acceptance for particles "<<EmcalAcceptance(eta,phi,0.)<<endl;
	
	Float_t deta=0.0, dphi=0.0, dr=100.0;
	Float_t deta1=0.0, dphi1=0.0, dr1=100.0;
	
	//cout<<i<<"  pt="<<pt<<"  eta="<<eta<<"  phi="<<phi<<endl;
	if(phi>1.396 && phi<3.316 && eta>-0.7 && eta<0.7){
	  sumpt+=pt;
	  //if(i<30)cout<<i<<" pt unit = "<<pt<<endl;

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
	  
	    //cout<<" out of the cone  "<<dr<<"    "<<deta<<"  deltaeta  "<<dphi<<"  dphi "<<i<<"  particle  "<<endl;
	    //cout<<" out of the cone  "<<dr1<<"    "<<deta1<<"  deltaeta1  "<<dphi1<<"  dphi1 "<<i<<"  particle  "<<endl;
	    ptallback+=pt;
	  }
	}
	//cout<<"  ipart  "<<ipart<<"  rhointegral  "<<rhoback <<endl;
      }		
    } // End loop on UnitArray 
  
  if(debug)cout<<"total area left "<<restarea<<endl;
  if(debug)cout<<"sumpt="<<sumpt<<endl;
  // if(acc==1 || acc1==1) rhoback= ptallback/restarea;
  //else rhoback=ptallback;

  rhoback= ptallback/restarea;
  if(debug)cout<<"rhoback    "<<rhoback<<"     "<<nJ<<"   "<<endl;

  return rhoback;
   
}


Double_t AliJetBkg::CalcRho(vector<fastjet::PseudoJet> inputParticles,Double_t rParamBkg,TString method){
  //calculate rho using the fastjet method

  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Bool_t debug  = header->GetDebug();     // debug option

  fastjet::Strategy strategy = header->GetStrategy();
  fastjet::RecombinationScheme recombScheme = header->GetRecombScheme();
  fastjet::JetAlgorithm algorithm = header->GetAlgorithm(); 
  fastjet::JetDefinition jetDef(algorithm, rParamBkg, recombScheme, strategy);

  // create an object that specifies how we to define the area
  fastjet::AreaDefinition areaDef;
  double ghostEtamax = header->GetGhostEtaMax(); 
  double ghostArea   = header->GetGhostArea(); 
  int    activeAreaRepeats = header->GetActiveAreaRepeats(); 
  
  // now create the object that holds info about ghosts

  if (method.Contains("Charg"))ghostEtamax=0.9;

  fastjet::GhostedAreaSpec ghost_spec(ghostEtamax, activeAreaRepeats, ghostArea);
  // and from that get an area definition
  fastjet::AreaType areaType = header->GetAreaType();
  areaDef = fastjet::AreaDefinition(areaType,ghost_spec);
  if(debug)cout<<"rParamBkg="<<rParamBkg<<"  ghostEtamax="<<ghostEtamax<<"  ghostArea="<<ghostArea<<" areadef="<<TString(areaDef.description())<<endl;
  //fastjet::ClusterSequenceArea clust_seq(inputParticles, jetDef);
  fastjet::ClusterSequenceArea clust_seq(inputParticles, jetDef,areaDef);
  TString comment = "Running FastJet algorithm for BKG calculation with the following setup. ";
  comment+= "Jet definition: ";
  comment+= TString(jetDef.description());
  //  comment+= ". Area definition: ";
  //comment+= TString(areaDef.description());
  comment+= ". Strategy adopted by FastJet: ";
  comment+= TString(clust_seq.strategy_string());
  header->SetComment(comment);
  if(debug){
    cout << "--------------------------------------------------------" << endl;
    cout << comment << endl;
    cout << "--------------------------------------------------------" << endl;
  }

  double ptmin = header->GetPtMin(); 
  vector<fastjet::PseudoJet> inclusiveJets = clust_seq.inclusive_jets(ptmin);
  vector<fastjet::PseudoJet> jets = sorted_by_pt(inclusiveJets); 

  if (debug) {
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
  //    double median, sigma, meanArea;
  //clust_seq.get_median_rho_and_sigma(inclusiveJets, range, false, median, sigma, meanArea, true);
  //fastjet::ActiveAreaSpec area_spec(ghostEtamax,activeAreaRepeats,ghostArea);

  // fastjet::ClusterSequenceActiveArea clust_seq_bkg(inputParticles, jetDef,area_spec);


  if(debug)cout<<"bkg in R="<<rParamBkg<<"  : "<<rho<<" range: Rap="<<rapMin<<","<<rapMax<<" --  phi="<<phiMin<<","<<phiMax<<endl;
  return rho;
}

Float_t  AliJetBkg::EtaToTheta(Float_t arg)
{
  //  return (180./TMath::Pi())*2.*atan(exp(-arg));
  return 2.*atan(exp(-arg));
}


Double_t  AliJetBkg::BkgFunction(Double_t */*x*/,Double_t */*par*/){
  //to be implemented--- (pT + Energy in EMCal Acceptance vs Multiplicity)
  return 1;
}


Bool_t AliJetBkg::EmcalAcceptance(const Float_t eta, const Float_t phi, const Float_t radius){
 
  Float_t meanPhi=190./180.*TMath::Pi()-110./180.*TMath::Pi()/2;
  Float_t deltaphi=110./180.*TMath::Pi();
  Float_t phicut=deltaphi/2.-radius;
  Float_t etacut=0.7-radius;
  //cout<<"  eta    "<<eta<<"  phi    "<<phi<<endl;
  //cout<<etacut<<"    "<<phicut<<"    "<<meanPhi<<"    "<<endl;
  if(TMath::Abs(eta)<etacut && TMath::Abs(phi-meanPhi)<phicut) return 1;
  else return 0; 
}
