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

/* History of cvs commits:
 *
 * $Log$
 *
 */

//_________________________________________________________________________
// Class for the analysis of gamma correlations (gamma-jet, 
// gamma-hadron and isolation cut.
// This class contains 3 main methods: one to fill lists of particles (ESDs) comming 
//  from the CTS (ITS+TPC) and the calorimeters;  the other one tags a candidate 
//  cluster as isolated;  the last one search in the 
//  corresponing calorimeter for the highest energy cluster, identified it as 
//  prompt photon;
//
//  Class created from old AliPHOSGammaJet 
//  (see AliRoot versions previous Release 4-09)
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include <TFile.h>
#include <TParticle.h>
#include <TH2.h>

#include "AliAnaGammaDirect.h" 
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "Riostream.h"
#include "AliLog.h"

ClassImp(AliAnaGammaDirect)

//____________________________________________________________________________
  AliAnaGammaDirect::AliAnaGammaDirect(const char *name) : 
    AliAnalysisTask(name,""), fChain(0), fESD(0),
    fOutputContainer(new TObjArray(100)), 
    fPrintInfo(0), fMinGammaPt(0.),
    fCalorimeter(""), fEMCALPID(0),fPHOSPID(0),
    fConeSize(0.),fPtThreshold(0.),fPtSumThreshold(0), 
    fNCones(0),fNPtThres(0),fMakeICMethod(0)
{
  //Ctor        
  TList * list = gDirectory->GetListOfKeys() ; 
  TIter next(list) ; 
  TH2F * h = 0 ;
  Int_t index ; 
  for (index = 0 ; index < list->GetSize()-1 ; index++) { 
    //-1 to avoid GammaJet Task
    h = dynamic_cast<TH2F*>(gDirectory->Get(list->At(index)->GetName())) ; 
    fOutputContainer->Add(h) ; 
  }
  
  for(Int_t i = 0; i < 10 ; i++){
    fConeSizes[i]=0;
    fPtThresholds[i]=0;
  }
  
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
  
}


//____________________________________________________________________________
AliAnaGammaDirect::AliAnaGammaDirect(const AliAnaGammaDirect & g) : 
  AliAnalysisTask(g), fChain(g.fChain), fESD(g.fESD),
  fOutputContainer(g. fOutputContainer),  fPrintInfo(g.fPrintInfo),
  fMinGammaPt(g.fMinGammaPt), fCalorimeter(g.fCalorimeter),
  fEMCALPID(g.fEMCALPID),fPHOSPID(g.fPHOSPID),
  fConeSize(g.fConeSize),
  fPtThreshold(g.fPtThreshold),
  fPtSumThreshold(g.fPtSumThreshold), 
  fNCones(g.fNCones),fNPtThres(g.fNPtThres),
  fMakeICMethod(g.fMakeICMethod)
{
  // cpy ctor
  SetName (g.GetName()) ; 
  SetTitle(g.GetTitle()) ; 

  for(Int_t i = 0; i < 10 ; i++){
    fConeSizes[i]=  g.fConeSizes[i];
    fPtThresholds[i]=   g.fPtThresholds[i];
  }
}

//____________________________________________________________________________
AliAnaGammaDirect::~AliAnaGammaDirect() 
{
  // Remove all pointers
  fOutputContainer->Clear() ; 
  delete fOutputContainer ;
  delete fhNGamma    ;  
  delete fhPhiGamma  ; 
  delete fhEtaGamma   ;  
  delete fhPtCandidate ;
  delete [] fhPtThresIsolated ;
  delete [] fhPtSumIsolated ;

}

//____________________________________________________________________________
void AliAnaGammaDirect::CreateParticleList(TClonesArray * pl, 
					TClonesArray * plCTS, 
					TClonesArray * plEMCAL,  
					TClonesArray * plPHOS){
  
  //Create a list of particles from the ESD. These particles have been measured 
  //by the Central Tracking system (TPC+ITS), PHOS and EMCAL 
  
  Int_t index = pl->GetEntries() ; 
  Int_t npar  = 0 ;
  Float_t *pid = new Float_t[AliPID::kSPECIESN];  
  AliDebug(3,"Fill particle lists");
  if(fPrintInfo)
    AliInfo(Form("fCalorimeter %s",fCalorimeter.Data()));

  Double_t v[3] ; //vertex ;
  fESD->GetVertex()->GetXYZ(v) ; 

  //########### PHOS ##############
  if(fCalorimeter == "PHOS"){

    Int_t begphos = fESD->GetFirstPHOSCluster();  
    Int_t endphos = fESD->GetFirstPHOSCluster() + 
      fESD->GetNumberOfPHOSClusters() ;  
    Int_t indexNePHOS = plPHOS->GetEntries() ;
    AliDebug(3,Form("First PHOS particle %d, last particle %d", begphos,endphos));

    if(fCalorimeter == "PHOS"){
      for (npar =  begphos; npar <  endphos; npar++) {//////////////PHOS track loop
	AliESDCaloCluster * clus = fESD->GetCaloCluster(npar) ; // retrieve track from esd

	//Create a TParticle to fill the particle list
	
	Float_t en = clus->GetClusterEnergy() ;
	Float_t *p = new Float_t();
	clus->GetGlobalPosition(p) ;
	TVector3 pos(p[0],p[1],p[2]) ; 
	Double_t phi  = pos.Phi();
	Double_t theta= pos.Theta();
	Double_t px = en*TMath::Cos(phi)*TMath::Sin(theta);;
	Double_t py = en*TMath::Sin(phi)*TMath::Sin(theta);
	Double_t pz = en*TMath::Cos(theta);

	TParticle * particle = new TParticle() ;
	particle->SetMomentum(px,py,pz,en) ;
	AliDebug(4,Form("PHOS clusters: pt %f, phi %f, eta %f", particle->Pt(),particle->Phi(),particle->Eta()));
	
	//Select only photons
	
	pid=clus->GetPid();
	//cout<<"pid "<<pid[AliPID::kPhoton]<<endl ;
	if( !fPHOSPID)
	  new((*plPHOS)[indexNePHOS++])   TParticle(*particle) ;
	else if( pid[AliPID::kPhoton] > 0.75)
	  new((*plPHOS)[indexNePHOS++])   TParticle(*particle) ;
      }
    }
  }

  //########### CTS (TPC+ITS) #####################
  Int_t begtpc   = 0 ;  
  Int_t endtpc   = fESD->GetNumberOfTracks() ;
  Int_t indexCh  = plCTS->GetEntries() ;
  AliDebug(3,Form("First CTS particle %d, last particle %d", begtpc,endtpc));
  
  for (npar =  begtpc; npar <  endtpc; npar++) {////////////// track loop
    AliESDtrack * track = fESD->GetTrack(npar) ; // retrieve track from esd
    //We want tracks fitted in the detectors:
    ULong_t status=AliESDtrack::kTPCrefit;
    status|=AliESDtrack::kITSrefit;
   
    //We want tracks whose PID bit is set:
    //     ULong_t status =AliESDtrack::kITSpid;
    //     status|=AliESDtrack::kTPCpid;
  
    if ( (track->GetStatus() & status) == status) {//Check if the bits we want are set
      // Do something with the tracks which were successfully
      // re-fitted 
      Double_t en = 0; //track ->GetTPCsignal() ;
      Double_t mom[3];
      track->GetPxPyPz(mom) ;
      Double_t px = mom[0];
      Double_t py = mom[1];
      Double_t pz = mom[2]; //Check with TPC people if this is correct.
      Int_t pdg = 11; //Give any charged PDG code, in this case electron.
      //I just want to tag the particle as charged
       TParticle * particle = new TParticle(pdg, 1, -1, -1, -1, -1, 
						 px, py, pz, en, v[0], v[1], v[2], 0);
  
      //TParticle * particle = new TParticle() ;
      //particle->SetMomentum(px,py,pz,en) ;
 
      new((*plCTS)[indexCh++])       TParticle(*particle) ;    
      new((*pl)[index++])           TParticle(*particle) ;
    }
  }
  
  //################ EMCAL ##############
  
  Int_t begem = fESD->GetFirstEMCALCluster();  
  Int_t endem = fESD->GetFirstEMCALCluster() + 
    fESD->GetNumberOfEMCALClusters() ;  
  Int_t indexNe  = plEMCAL->GetEntries() ; 
  
  AliDebug(3,Form("First EMCAL particle %d, last particle %d",begem,endem));
    
    for (npar =  begem; npar <  endem; npar++) {//////////////EMCAL track loop
      AliESDCaloCluster * clus = fESD->GetCaloCluster(npar) ; // retrieve track from esd
      Int_t clustertype= clus->GetClusterType();
      if(clustertype == AliESDCaloCluster::kClusterv1){
	Float_t en = clus->GetClusterEnergy() ;
	Float_t *p = new Float_t();
	clus->GetGlobalPosition(p) ;
	TVector3 pos(p[0],p[1],p[2]) ;
	Double_t phi  = pos.Phi();
	Double_t theta= pos.Theta();
	Double_t px = en*TMath::Cos(phi)*TMath::Sin(theta);;
	Double_t py = en*TMath::Sin(phi)*TMath::Sin(theta);
	Double_t pz = en*TMath::Cos(theta);

	pid=clus->GetPid();
	if(fCalorimeter == "EMCAL")
	  {
	    TParticle * particle = new TParticle() ;
	    particle->SetMomentum(px,py,pz,en) ;
	    AliDebug(4,Form("EMCAL clusters: pt %f, phi %f, eta %f", particle->Pt(),particle->Phi(),particle->Eta()));
	    if(!fEMCALPID) //Only identified particles
	      new((*plEMCAL)[indexNe++])       TParticle(*particle) ;
	    else if(pid[AliPID::kPhoton] > 0.75)
	      new((*plEMCAL)[indexNe++])       TParticle(*particle) ;	    
	  }
	else
	  {
	    Int_t pdg = 0;
	    if(fEMCALPID) 
	      {
		if( pid[AliPID::kPhoton] > 0.75) //This has to be fixen.
		  pdg = 22;
		else if( pid[AliPID::kPi0] > 0.75)
		  pdg = 111;
	      }
	    else
	      pdg = 22; //No PID, assume all photons
	    
	    TParticle * particle = new TParticle(pdg, 1, -1, -1, -1, -1, 
						 px, py, pz, en, v[0], v[1], v[2], 0);
	    AliDebug(4,Form("EMCAL clusters: pt %f, phi %f, eta %f", particle->Pt(),particle->Phi(),particle->Eta()));
	    
	    new((*plEMCAL)[indexNe++])       TParticle(*particle) ; 
	    new((*pl)[index++])           TParticle(*particle) ;
	  }
      }
    }

    AliDebug(3,"Particle lists filled");
    
}



//____________________________________________________________________________
void AliAnaGammaDirect::Exec(Option_t *) 
{
  
  // Processing of one event
    
  //Get ESDs
  Long64_t entry = fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if (fPrintInfo) 
     AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 

  //CreateTLists with arrays of TParticles. Filled with particles only relevant for the analysis.

  TClonesArray * particleList = new TClonesArray("TParticle",1000); // All particles refitted in CTS and detected in EMCAL (jet)
  TClonesArray * plCTS         = new TClonesArray("TParticle",1000); // All particles refitted in Central Tracking System (ITS+TPC)
  TClonesArray * plNe         = new TClonesArray("TParticle",1000);   // All particles measured in Jet Calorimeter (EMCAL)
  TClonesArray * plPHOS     = new TClonesArray("TParticle",1000);  // All particles measured in PHOS as Gamma calorimeter
  TClonesArray * plEMCAL   = new TClonesArray("TParticle",1000);  // All particles measured in EMCAL as Gamma calorimeter

  TParticle *pGamma = new TParticle(); //It will contain the kinematics of the found prompt gamma
  TParticle *pLeading = new TParticle(); //It will contain the kinematics of the found leading particle
  
  //Fill lists with photons, neutral particles and charged particles
  //look for the highest energy photon in the event inside fCalorimeter

  AliDebug(2, "Fill particle lists, get prompt gamma");

  //Fill particle lists 
  CreateParticleList(particleList, plCTS,plEMCAL,plPHOS); 

  if(fCalorimeter == "PHOS")
    plNe = plPHOS;
  if(fCalorimeter == "EMCAL")
    plNe = plEMCAL;


  //_______________Analysis 1__________________________
  //Look for the prompt photon in the selected calorimeter
  if(fMakeICMethod < 3){
    
    Bool_t iIsInPHOSorEMCAL = kFALSE ; //To check if Gamma was in any calorimeter
    //Search highest energy prompt gamma in calorimeter
    GetPromptGamma(plNe,  plCTS, pGamma, iIsInPHOSorEMCAL) ; 
    
    AliDebug(1, Form("Is Gamma in %s? %d",fCalorimeter.Data(),iIsInPHOSorEMCAL));
    
    //If there is any photon candidate in fCalorimeter
    if(iIsInPHOSorEMCAL){
      if (fPrintInfo)
	AliInfo(Form("Prompt Gamma: pt %f, phi %f, eta %f", pGamma->Pt(),pGamma->Phi(),pGamma->Eta())) ;
      
    }//Gamma in Calo
    
  }//Analysis 1


  //_______________Analysis 2__________________________
  //Look for the prompt photon in the selected calorimeter
  //Isolation Cut Analysis for both methods and different pt cuts and cones
  if(fMakeICMethod == 3){
    
    for(Int_t ipr = 0; ipr < plNe->GetEntries() ; ipr ++ ){
      TParticle * pCandidate = dynamic_cast<TParticle *>(plNe->At(ipr)) ;
      
      if(pCandidate->Pt() > fMinGammaPt){
	
	Bool_t  icPtThres   = kFALSE;
	Bool_t  icPtSum     = kFALSE;
	
	Float_t ptC             = pCandidate->Pt() ;
	
	fhPtCandidate->Fill(ptC);
	
	for(Int_t icone = 0; icone<fNCones; icone++){
	  fConeSize  = fConeSizes[icone] ;
	  Float_t coneptsum = 0 ;
	  for(Int_t ipt = 0; ipt<fNPtThres;ipt++){ 
	    fPtThreshold   =  fPtThresholds[ipt] ;
	    MakeIsolationCut(plCTS,plNe, pCandidate, ipr, icPtThres, icPtSum,coneptsum);
	    AliDebug(4,Form("Candidate pt %f, pt in cone %f, Isolated? ICPt %d, ICSum %d",
			    pCandidate->Pt(), coneptsum, icPtThres, icPtSum));
	    
	    fhPtThresIsolated[icone][ipt]->Fill(ptC); 
	  }//pt thresh loop
	  fhPtSumIsolated[icone]->Fill(ptC,coneptsum) ;
	}//cone size loop
      }//min pt candidate
    }//candidate loop
  }//Analysis 2
  
  AliDebug(2, "End of analysis, delete pointers");
  
  particleList->Delete() ; 
  plCTS->Delete() ;
  plNe->Delete() ;
  plPHOS->Delete() ;
  pLeading->Delete();
  pGamma->Delete();

  delete plNe ;
  delete plCTS ;
  delete particleList ;
  //  delete pLeading;
  //  delete pGamma;

  PostData(0, fOutputContainer);
}    


//____________________________________________________________________________
void AliAnaGammaDirect::GetPromptGamma(TClonesArray * pl, TClonesArray * plCTS, TParticle *pGamma, Bool_t &Is) const 
{
  //Search for the prompt photon in Calorimeter with pt > fMinGammaPt

  Double_t pt = 0;
  Int_t index = -1; 
  for(Int_t ipr = 0;ipr < pl->GetEntries() ; ipr ++ ){
    TParticle * particle = dynamic_cast<TParticle *>(pl->At(ipr)) ;

    if((particle->Pt() > fMinGammaPt) && (particle->Pt() > pt)){
      index = ipr ;
      pt = particle->Pt();
      pGamma->SetMomentum(particle->Px(),particle->Py(),particle->Pz(),particle->Energy());
      AliDebug(4,Form("Cluster in calo: pt %f, phi %f, eta %f", pGamma->Pt(),pGamma->Phi(),pGamma->Eta())) ;
      Is  = kTRUE;
    }
  }

  //Do Isolation?
  if(fMakeICMethod && Is)
    {
      Float_t coneptsum = 0 ;
      Bool_t  icPtThres   = kFALSE;
      Bool_t  icPtSum     = kFALSE;
      MakeIsolationCut(plCTS,pl, pGamma, index, 
		       icPtThres, icPtSum,coneptsum);
      if(fMakeICMethod == 1) //Pt thres method
	Is = icPtThres ;
      if(fMakeICMethod == 2) //Pt cone sum method
	Is = icPtSum ;
    }
  
  if(Is){
    AliDebug(3,Form("Cluster with p_{T} larger than %f found in calorimeter ", fMinGammaPt)) ;
    AliDebug(3,Form("Gamma: pt %f, phi %f, eta %f", pGamma->Pt(),pGamma->Phi(),pGamma->Eta())) ;
    //Fill prompt gamma histograms
    fhNGamma  ->Fill(pGamma->Pt());
    fhPhiGamma->Fill( pGamma->Pt(),pGamma->Phi());
    fhEtaGamma->Fill(pGamma->Pt(),pGamma->Eta());
  }
  else
    AliDebug(1,Form("NO Cluster with pT larger than %f found in calorimeter ", fMinGammaPt)) ;
}

  //____________________________________________________________________________
void AliAnaGammaDirect::Init(const Option_t * )
{
  // Initialisation of branch container 
 
  AliDebug(2,Form("*** Initialization of %s", GetName())) ; 
  
  // Get input data
  fChain = dynamic_cast<TChain *>(GetInputData(0)) ;
  if (!fChain) {
    AliError(Form("Input 0 for %s not found\n", GetName()));
    return ;
  }
  
  if (!fESD) {
    // One should first check if the branch address was taken by some other task
    char ** address = (char **)GetBranchAddress(0, "ESD") ;
    if (address) 
      fESD = (AliESD *)(*address) ; 
    if (!fESD) 
      fChain->SetBranchAddress("ESD", &fESD) ;  
  }
  // The output objects will be written to 
  TDirectory * cdir = gDirectory ; 
  // Open a file for output #0
  char outputName[1024] ; 
  sprintf(outputName, "%s.root", GetName() ) ; 
  OpenFile(0, outputName , "RECREATE") ; 
  if (cdir) 
    cdir->cd() ; 

 //  //Initialize the parameters of the analysis.
  fCalorimeter="PHOS";
  fPrintInfo           = kTRUE;
  fMinGammaPt  = 10. ;

  //Fill particle lists when PID is ok
  fEMCALPID = kFALSE;
  fPHOSPID = kFALSE;

  fConeSize             = 0.2 ; 
  fPtThreshold         = 2.0; 
  fPtSumThreshold  = 1.; 

  fNCones           = 4 ; 
  fNPtThres         = 4 ; 
  fConeSizes[0] = 0.1; fConeSizes[0] = 0.2; fConeSizes[2] = 0.3; fConeSizes[3] = 0.4;
  fPtThresholds[0]=1.; fPtThresholds[0]=2.; fPtThresholds[0]=3.; fPtThresholds[0]=4.;

  fMakeICMethod = 1; // 0 don't isolate, 1 pt thresh method, 2 cone pt sum method, 3 make isolation study

  //Initialization of histograms 
  MakeHistos() ;
}

//__________________________________________________________________
void  AliAnaGammaDirect::MakeIsolationCut(TClonesArray * plCTS, 
					   TClonesArray * plNe, 
					   TParticle * pCandidate, 
					   Int_t index, 
					   Bool_t  &icmpt,  Bool_t  &icms, 
					   Float_t &coneptsum) const 
{  
  //Search in cone around a candidate particle if it is isolated 
  Float_t phiC  = pCandidate->Phi() ;
  Float_t etaC = pCandidate->Eta() ;
  Float_t pt     = -100. ;
  Float_t eta   = -100.  ;
  Float_t phi    = -100.  ;
  Float_t rad   = -100 ;
  Int_t    n        = 0 ;
  TParticle * particle  = new TParticle();

  coneptsum = 0; 
  icmpt = kFALSE;
  icms  = kFALSE;

  //Check charged particles in cone.
  for(Int_t ipr = 0;ipr < plCTS->GetEntries() ; ipr ++ ){
    particle = dynamic_cast<TParticle *>(plCTS->At(ipr)) ;
    pt    = particle->Pt();
    eta  = particle->Eta();
    phi  = particle->Phi() ;
    
    //Check if there is any particle inside cone with pt larger than  fPtThreshold
    rad = TMath::Sqrt(TMath::Power((eta-etaC),2) +
		      TMath::Power((phi-phiC),2));
    if(rad<fConeSize){
      AliDebug(3,Form("charged in cone pt %f, phi %f, eta %f, R %f ",pt,phi,eta,rad));
      coneptsum+=pt;
      if(pt > fPtThreshold ) n++;
    }
  }// charged particle loop
  
  //Check neutral particles in cone.
  for(Int_t ipr = 0;ipr < plNe->GetEntries() ; ipr ++ ){
    if(ipr != index){//Do not count the candidate
      particle = dynamic_cast<TParticle *>(plNe->At(ipr)) ;
      pt    = particle->Pt();
      eta  = particle->Eta();
      phi  = particle->Phi() ;
      
      //Check if there is any particle inside cone with pt larger than  fPtThreshold
      rad = TMath::Sqrt(TMath::Power((eta-etaC),2) +
			TMath::Power((phi-phiC),2));
      if(rad<fConeSize){
	AliDebug(3,Form("charged in cone pt %f, phi %f, eta %f, R %f ",pt,phi,eta,rad));
	coneptsum+=pt;
	if(pt > fPtThreshold ) n++;
      }
    }
  }// neutral particle loop
  
  if(n == 0) 
    icmpt =  kTRUE ;
  if(coneptsum < fPtSumThreshold)
    icms  =  kTRUE ;

}

//___________________________________________________________________
void AliAnaGammaDirect::MakeHistos()
{
  // Create histograms to be saved in output file and 
  // stores them in fOutputContainer
  
  fOutputContainer = new TObjArray(10000) ;

  //Histograms of highest gamma identified in Event
  fhNGamma  = new TH1F("NGamma","Number of #gamma over PHOS",240,0,120); 
  fhNGamma->SetYTitle("N");
  fhNGamma->SetXTitle("p_{T #gamma}(GeV/c)");
  fOutputContainer->Add(fhNGamma) ; 
  
  fhPhiGamma  = new TH2F
    ("PhiGamma","#phi_{#gamma}",200,0,120,200,0,7); 
  fhPhiGamma->SetYTitle("#phi");
  fhPhiGamma->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhPhiGamma) ; 
  
  fhEtaGamma  = new TH2F
    ("EtaGamma","#phi_{#gamma}",200,0,120,200,-0.8,0.8); 
  fhEtaGamma->SetYTitle("#eta");
  fhEtaGamma->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhEtaGamma) ;

  if( fMakeICMethod== 3 ){

    //Isolation cut histograms
    fhPtCandidate  = new TH1F
      ("PtCandidate","p_{T} of candidate particles for isolation",240,0,120); 
    fhPtCandidate->SetXTitle("p_{T} (GeV/c)");
    fOutputContainer->Add(fhPtCandidate) ;
    
    char name[128];
    char title[128];
    for(Int_t icone = 0; icone<fNCones; icone++){
      sprintf(name,"PtSumIsolated_Cone_%d",icone);
      sprintf(title,"Candidate cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
      fhPtSumIsolated[icone]  = new TH2F(name, title,240,0,120,120,0,10);
      fhPtSumIsolated[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
      fhPtSumIsolated[icone]->SetXTitle("p_{T} (GeV/c)");
      fOutputContainer->Add(fhPtSumIsolated[icone]) ; 
      
      for(Int_t ipt = 0; ipt<fNPtThres;ipt++){ 
	sprintf(name,"PtThresIsol_Cone_%d_Pt%d",icone,ipt);
	sprintf(title,"Isolated candidate p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	fhPtThresIsolated[icone][ipt]  = new TH1F(name, title,240,0,120);
	fhPtThresIsolated[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	fOutputContainer->Add(fhPtThresIsolated[icone][ipt]) ; 
      }//icone loop
  }//ipt loop
  }
}

void AliAnaGammaDirect::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;

  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  printf("Cone Size               =     %f\n", fConeSize) ; 
  printf("pT threshold           =     %f\n", fPtThreshold) ;
  printf("pT sum threshold   =     %f\n", fPtSumThreshold) ; 
  printf("Min Gamma pT      =     %f\n", fMinGammaPt) ; 
  printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ; 
} 

void AliAnaGammaDirect::Terminate(Option_t *)
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
    

}
