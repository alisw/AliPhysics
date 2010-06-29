/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes iGetEntriesFast(s hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id: AliAnaPi0EbE.cxx 28688 2008-09-11 15:04:07Z gconesab $ */

//_________________________________________________________________________
// Class for the analysis of high pT pi0 event by event
// Pi0 identified by one of the following:
//  -Invariant mass of 2 cluster in calorimeter
//  -Shower shape analysis in calorimeter
//  -Invariant mass of one cluster in calorimeter and one photon reconstructed in TPC (in near future)
//
// -- Author: Gustavo Conesa (LNF-INFN) &  Raphaelle Ichou (SUBATECH)
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TList.h>
#include <TClonesArray.h>
//#include <TObjString.h>
#include <TH2F.h>
//#include "Riostream.h"

// --- Analysis system --- 
#include "AliAnaPi0EbE.h" 
#include "AliCaloTrackReader.h"
#include "AliIsolationCut.h"
#include "AliNeutralMesonSelection.h"
#include "AliCaloPID.h"
#include "AliMCAnalysisUtils.h"
#include "AliStack.h"
#include "AliFiducialCut.h"
#include "TParticle.h"
#include "AliAODCaloCluster.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnaPi0EbE)
  
//____________________________________________________________________________
AliAnaPi0EbE::AliAnaPi0EbE() : 
AliAnaPartCorrBaseClass(),  fAnaType(kIMCalo),fCalorimeter(""),
fMinDist(0.),fMinDist2(0.),fMinDist3(0.),	
fInputAODGammaConv(0x0),fInputAODGammaConvName(""),
fhPtPi0(0),fhPhiPi0(0),fhEtaPi0(0), 
fhPtMCNoPi0(0),fhPhiMCNoPi0(0),fhEtaMCNoPi0(0), 
fhPtMCPi0(0),fhPhiMCPi0(0),fhEtaMCPi0(0)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();
  
}
/*
//____________________________________________________________________________
AliAnaPi0EbE::AliAnaPi0EbE(const AliAnaPi0EbE & p) : 
AliAnaPartCorrBaseClass(p),  fAnaType(p.fAnaType), fCalorimeter(p.fCalorimeter),
fMinDist(p.fMinDist),fMinDist2(p.fMinDist2), fMinDist3(p.fMinDist3),
fInputAODGammaConv(new TClonesArray(*p.fInputAODGammaConv)),
fInputAODGammaConvName(p.fInputAODGammaConvName),
fhPtPi0(p.fhPtPi0),fhPhiPi0(p.fhPhiPi0),fhEtaPi0(p.fhEtaPi0), 
fhPtMCNoPi0(p.fhPtMCNoPi0),fhPhiMCNoPi0(p.fhPhiMCNoPi0),fhEtaMCNoPi0(p.fhEtaMCNoPi0), 
fhPtMCPi0(p.fhPtMCPi0),fhPhiMCPi0(p.fhPhiMCPi0),fhEtaMCPi0(p.fhEtaMCPi0) 
{
  // cpy ctor
}

//_________________________________________________________________________
AliAnaPi0EbE & AliAnaPi0EbE::operator = (const AliAnaPi0EbE & p)
{
  // assignment operator
  
  if(&p == this) return *this;
  
  fAnaType     = p.fAnaType ;
  fCalorimeter = p.fCalorimeter ;
  fMinDist     = p.fMinDist;
  fMinDist2    = p.fMinDist2;
  fMinDist3    = p.fMinDist3;
  
  fInputAODGammaConv     = new TClonesArray(*p.fInputAODGammaConv) ;
  fInputAODGammaConvName = p.fInputAODGammaConvName ;
  
  fhPtPi0 = p.fhPtPi0; fhPhiPi0 = p.fhPhiPi0; fhEtaPi0 = p.fhEtaPi0;  
  fhPtMCNoPi0 = p.fhPtMCNoPi0; fhPhiMCNoPi0 = p.fhPhiMCNoPi0; fhEtaMCNoPi0 = p.fhEtaMCPi0;  
  fhPtMCPi0 = p.fhPtMCPi0; fhPhiMCPi0 = p.fhPhiMCPi0; fhEtaMCPi0 = p.fhEtaMCPi0;  
  
  return *this;
  
}
*/
//____________________________________________________________________________
AliAnaPi0EbE::~AliAnaPi0EbE() 
{
  //dtor
  if(fInputAODGammaConv){
    fInputAODGammaConv->Clear() ; 
    delete fInputAODGammaConv ;
  }
}

//________________________________________________________________________
TList *  AliAnaPi0EbE::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("Pi0EbEHistos") ; 
  
  Int_t nptbins  = GetHistoPtBins();
  Int_t nphibins = GetHistoPhiBins();
  Int_t netabins = GetHistoEtaBins();
  Float_t ptmax  = GetHistoPtMax();
  Float_t phimax = GetHistoPhiMax();
  Float_t etamax = GetHistoEtaMax();
  Float_t ptmin  = GetHistoPtMin();
  Float_t phimin = GetHistoPhiMin();
  Float_t etamin = GetHistoEtaMin();	
  
  fhPtPi0  = new TH1F("hPtPi0","Number of identified  #pi^{0} decay",nptbins,ptmin,ptmax); 
  fhPtPi0->SetYTitle("N");
  fhPtPi0->SetXTitle("p_{T #pi^{0}}(GeV/c)");
  outputContainer->Add(fhPtPi0) ; 
  
  fhPhiPi0  = new TH2F
    ("hPhiPi0","#phi_{#pi^{0}}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
  fhPhiPi0->SetYTitle("#phi");
  fhPhiPi0->SetXTitle("p_{T #pi^{0}} (GeV/c)");
  outputContainer->Add(fhPhiPi0) ; 
  
  fhEtaPi0  = new TH2F
    ("hEtaPi0","#phi_{#pi^{0}}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
  fhEtaPi0->SetYTitle("#eta");
  fhEtaPi0->SetXTitle("p_{T #pi^{0}} (GeV/c)");
  outputContainer->Add(fhEtaPi0) ;
  
  if(IsDataMC()) {
    if((GetReader()->GetDataType() == AliCaloTrackReader::kMC && fAnaType!=kSSCalo) || 
       GetReader()->GetDataType() != AliCaloTrackReader::kMC){
      
      fhPtMCPi0  = new TH1F("hPtMCPi0","Identified pi0 from pi0",nptbins,ptmin,ptmax); 
      fhPtMCPi0->SetYTitle("N");
      fhPtMCPi0->SetXTitle("p_{T #pi^{0}}(GeV/c)");
      outputContainer->Add(fhPtMCPi0) ; 
      
      fhPhiMCPi0  = new TH2F
	("hPhiMCPi0","Identified pi0 from pi0",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiMCPi0->SetYTitle("#phi");
      fhPhiMCPi0->SetXTitle("p_{T #pi^{0}} (GeV/c)");
      outputContainer->Add(fhPhiMCPi0) ; 
      
      fhEtaMCPi0  = new TH2F
	("hEtaMCPi0","Identified pi0 from pi0",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaMCPi0->SetYTitle("#eta");
      fhEtaMCPi0->SetXTitle("p_{T #pi^{0}} (GeV/c)");
      outputContainer->Add(fhEtaMCPi0) ;
      
      fhPtMCNoPi0  = new TH1F("hPtMCNoPi0","Identified pi0 not from pi0",nptbins,ptmin,ptmax); 
      fhPtMCNoPi0->SetYTitle("N");
      fhPtMCNoPi0->SetXTitle("p_{T #pi^{0}}(GeV/c)");
      outputContainer->Add(fhPtMCNoPi0) ; 
      
      fhPhiMCNoPi0  = new TH2F
	("hPhiMCNoPi0","Identified pi0 not from pi0",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiMCNoPi0->SetYTitle("#phi");
      fhPhiMCNoPi0->SetXTitle("p_{T #pi^{0}} (GeV/c)");
      outputContainer->Add(fhPhiMCNoPi0) ; 
      
      fhEtaMCNoPi0  = new TH2F
	("hEtaMCNoPi0","Identified pi0 not from pi0",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaMCNoPi0->SetYTitle("#eta");
      fhEtaMCNoPi0->SetXTitle("p_{T #pi^{0}} (GeV/c)");
      outputContainer->Add(fhEtaMCNoPi0) ;
      
    }
  }//Histos with MC
  
  
  //Keep neutral meson selection histograms if requiered
  //Setting done in AliNeutralMesonSelection
  
  if(fAnaType!=kSSCalo && GetNeutralMesonSelection()){
    
    TList * nmsHistos = GetNeutralMesonSelection()->GetCreateOutputObjects() ;
    if(GetNeutralMesonSelection()->AreNeutralMesonSelectionHistosKept())
      for(Int_t i = 0; i < nmsHistos->GetEntries(); i++) outputContainer->Add(nmsHistos->At(i)) ;
	delete nmsHistos;
	  
  }
  
  //Save parameters used for analysis
//  TString parList ; //this will be list of parameters used for this analysis.
//  char onePar[255] ;
//  
//  sprintf(onePar,"--- AliAnaPi0EbE ---\n") ;
//  parList+=onePar ;	
//  sprintf(onePar,"fAnaType=%d (Pi0 selection type) \n",fAnaType) ;
//  parList+=onePar ;
//  
//  if(fAnaType == kSSCalo){
//    sprintf(onePar,"Calorimeter: %s\n",fCalorimeter.Data()) ;
//    parList+=onePar ;
//    sprintf(onePar,"fMinDist =%2.2f (Minimal distance to bad channel to accept cluster) \n",fMinDist) ;
//    parList+=onePar ;
//    sprintf(onePar,"fMinDist2=%2.2f (Cuts on Minimal distance to study acceptance evaluation) \n",fMinDist2) ;
//    parList+=onePar ;
//    sprintf(onePar,"fMinDist3=%2.2f (One more cut on distance used for acceptance-efficiency study) \n",fMinDist3) ;
//    parList+=onePar ;
//  }
//  
//  //Get parameters set in base class.
//  parList += GetBaseParametersList() ;
//  
//  //Get parameters set in PID class.
//  if(fAnaType == kSSCalo) parList += GetCaloPID()->GetPIDParametersList() ;
//  
//  TObjString *oString= new TObjString(parList) ;
//  outputContainer->Add(oString);
  
  return outputContainer ;
  
}

//__________________________________________________________________
void  AliAnaPi0EbE::MakeAnalysisFillAOD() 
{
  //Do analysis and fill aods
  
  switch(fAnaType) 
    {
    case kIMCalo:
      MakeInvMassInCalorimeter();
      break;
      
    case kSSCalo:
      MakeShowerShapeIdentification();
      break;
      
    case kIMCaloTracks:
      MakeInvMassInCalorimeterAndCTS();
      break;
      
    }
}

//__________________________________________________________________
void  AliAnaPi0EbE::MakeInvMassInCalorimeter() 
{
  //Do analysis and fill aods
  //Search for the photon decay in calorimeters
  //Read photon list from AOD, produced in class AliAnaPhoton
  //Check if 2 photons have the mass of the pi0.
  
  TLorentzVector mom1;
  TLorentzVector mom2;
  TLorentzVector mom ;
  Int_t tag1 = 0;
  Int_t tag2 = 0;
  Int_t tag  = 0;
  
  if(!GetInputAODBranch()){
    printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - No input calo photons in AOD with name branch < %s >, STOP \n",GetInputAODName().Data());
    abort();
  }
  for(Int_t iphoton = 0; iphoton < GetInputAODBranch()->GetEntriesFast(); iphoton++){
    AliAODPWG4Particle * photon1 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(iphoton));
    mom1 = *(photon1->Momentum());
    
    
    for(Int_t jphoton = iphoton+1; jphoton < GetInputAODBranch()->GetEntriesFast()-1; jphoton++){
      
      AliAODPWG4Particle * photon2 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(jphoton));
      mom2 = *(photon2->Momentum());
	  Int_t input = -1;	//if -1 photons come from different files, not a pi0
	  if(photon1->GetInputFileIndex() == photon2->GetInputFileIndex()) input = photon1->GetInputFileIndex();
	  
      //Select good pair (good phi, pt cuts, aperture and invariant mass)
      if(GetNeutralMesonSelection()->SelectPair(mom1, mom2))
	{
	  if(GetDebug()>1) 
	    printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - Selected gamma pair: pt %f, phi %f, eta%f \n",(mom1+mom2).Pt(), (mom1+mom2).Phi()*180./3.1416, (mom1+mom2).Eta());
	  
	  //Play with the MC stack if available
	  if(IsDataMC()){
	    //Check origin of the candidates
		Int_t  label1 = photon1->GetLabel();
		Int_t  label2 = photon2->GetLabel();
	    tag1 = GetMCAnalysisUtils()->CheckOrigin(label1, GetReader(), photon1->GetInputFileIndex());
	    tag2 = GetMCAnalysisUtils()->CheckOrigin(label2, GetReader(), photon2->GetInputFileIndex());
	    
	    if(GetDebug() > 0) printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - Origin of: photon1 %d; photon2 %d \n",tag1, tag2);
	    if(GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCPi0Decay) && GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCPi0Decay)){
	      
	      //Check if pi0 mother is the same
		  if(GetReader()->ReadStack()){ 
				TParticle * mother1 = GetMCStack()->Particle(label1);//photon in kine tree
				label1 = mother1->GetFirstMother();
				//mother1 = GetMCStack()->Particle(label1);//pi0
				
				TParticle * mother2 = GetMCStack()->Particle(label2);//photon in kine tree
				label2 = mother2->GetFirstMother();
				//mother2 = GetMCStack()->Particle(label2);//pi0
		  }
		  else if(GetReader()->ReadAODMCParticles() && (input > -1)){
				AliAODMCParticle * mother1 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles(photon1->GetInputFileIndex()))->At(label1);//photon in kine tree
				label1 = mother1->GetMother();
				//mother1 = GetMCStack()->Particle(label1);//pi0
				AliAODMCParticle * mother2 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles(photon2->GetInputFileIndex()))->At(label2);//photon in kine tree
				label2 = mother2->GetMother();
				//mother2 = GetMCStack()->Particle(label2);//pi0
		  }
			
	      //printf("mother1 %d, mother2 %d\n",label1,label2);
	      if(label1 == label2)
			  GetMCAnalysisUtils()->SetTagBit(tag,AliMCAnalysisUtils::kMCPi0);
	    }
	  }//Work with stack also   
	  
	  //Create AOD for analysis
	  mom = mom1+mom2;
	  AliAODPWG4Particle pi0 = AliAODPWG4Particle(mom);
	  //pi0.SetLabel(calo->GetLabel(0));
	  pi0.SetPdg(AliCaloPID::kPi0);
	  pi0.SetDetector(photon1->GetDetector());
	  pi0.SetTag(tag);  
	  //Set the indeces of the original caloclusters  
	  pi0.SetCaloLabel(photon1->GetCaloLabel(0), photon2->GetCaloLabel(0));
	  pi0.SetInputFileIndex(input);
	  AddAODParticle(pi0);
	}//pi0
    }//2n photon loop
    
  }//1st photon loop
  
  if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - End fill AODs \n");  
  
}

//__________________________________________________________________
void  AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() 
{
  //Do analysis and fill aods
  //Search for the photon decay in calorimeters
  //Read photon list from AOD, produced in class AliAnaPhoton and AliGammaConversion
  //Check if 2 photons have the mass of the pi0.
  
  TLorentzVector mom1;
  TLorentzVector mom2;
  TLorentzVector mom ;
  Int_t tag1 = 0;
  Int_t tag2 = 0;
  Int_t tag  = 0;
  
  if(!GetInputAODBranch()){
    printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - No input calo photons in AOD branch with name < %s > , STOP\n",GetInputAODName().Data());
    abort();
  }
  for(Int_t iphoton = 0; iphoton < GetInputAODBranch()->GetEntriesFast(); iphoton++){
    AliAODPWG4Particle * photon1 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(iphoton));
    mom1 = *(photon1->Momentum());
    
    //Play with the MC stack if available
    fInputAODGammaConv = (TClonesArray *) GetReader()->GetOutputEvent()->FindListObject(fInputAODGammaConvName);
    if(!fInputAODGammaConv) {
      printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - No input gamma conversions in AOD branch with name < %s >, STOP\n",fInputAODGammaConvName.Data());
      abort();	
    }
    for(Int_t jphoton = iphoton+1; jphoton < fInputAODGammaConv->GetEntriesFast()-1; jphoton++){
      AliAODPWG4Particle * photon2 =  (AliAODPWG4Particle*) (fInputAODGammaConv->At(jphoton));
      mom2 = *(photon2->Momentum());
		
	  Int_t input = -1;	//if -1 photons come from different files, not a pi0
	  if(photon1->GetInputFileIndex() == photon2->GetInputFileIndex()) input = photon1->GetInputFileIndex();

      //Select good pair (good phi, pt cuts, aperture and invariant mass)
      if(GetNeutralMesonSelection()->SelectPair(mom1, mom2)){
	if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - Selected gamma pair: pt %f, phi %f, eta%f\n",(mom1+mom2).Pt(), (mom1+mom2).Phi()*180./3.1416, (mom1+mom2).Eta());
	
	if(IsDataMC()){
	  Int_t	label1 = photon1->GetLabel();
	  Int_t	label2 = photon2->GetLabel();
	  tag1 = GetMCAnalysisUtils()->CheckOrigin(label1, GetReader(), photon1->GetInputFileIndex());
	  tag2 = GetMCAnalysisUtils()->CheckOrigin(label2, GetReader(), photon2->GetInputFileIndex());
	  if(GetDebug() > 0) printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - Origin of: photon1 %d; photon2 %d \n",tag1, tag2);
	  if(GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCPi0Decay) && GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCPi0Decay)){
	    //Check if pi0 mother is the same
	  
		if(GetReader()->ReadStack()){ 
			TParticle * mother1 = GetMCStack()->Particle(label1);//photon in kine tree
			label1 = mother1->GetFirstMother();
			//mother1 = GetMCStack()->Particle(label1);//pi0
	    
			TParticle * mother2 = GetMCStack()->Particle(label2);//photon in kine tree
			label2 = mother2->GetFirstMother();
			//mother2 = GetMCStack()->Particle(label2);//pi0
	    }
		else if(GetReader()->ReadAODMCParticles() && (input > -1)){
			AliAODMCParticle * mother1 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles(photon1->GetInputFileIndex()))->At(label1);//photon in kine tree
			label1 = mother1->GetMother();
			//mother1 = GetMCStack()->Particle(label1);//pi0
			AliAODMCParticle * mother2 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles(photon2->GetInputFileIndex()))->At(label2);//photon in kine tree
			label2 = mother2->GetMother();
			//mother2 = GetMCStack()->Particle(label2);//pi0
		}
		  
		//printf("mother1 %d, mother2 %d\n",label1,label2);
	    if(label1 == label2)
	      GetMCAnalysisUtils()->SetTagBit(tag,AliMCAnalysisUtils::kMCPi0);
	  }
	}//Work with stack also   
	
	//Create AOD for analysis
	mom = mom1+mom2;
	AliAODPWG4Particle pi0 = AliAODPWG4Particle(mom);
	//pi0.SetLabel(calo->GetLabel(0));
	pi0.SetPdg(AliCaloPID::kPi0);
	pi0.SetDetector(photon1->GetDetector());
	pi0.SetTag(tag);
	//Set the indeces of the original tracks or caloclusters  
	pi0.SetCaloLabel(photon1->GetCaloLabel(0), -1);
	pi0.SetTrackLabel(photon2->GetTrackLabel(0), photon2->GetTrackLabel(1));
	pi0.SetInputFileIndex(input);
	AddAODParticle(pi0);
      }//pi0
    }//2n photon loop
    
  }//1st photon loop
  
  if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - End fill AODs \n");  
  
}


//__________________________________________________________________
void  AliAnaPi0EbE::MakeShowerShapeIdentification() 
{
  //Search for pi0 in fCalorimeter with shower shape analysis 
  
  TObjArray * pl = 0x0; 
  
  //Get vertex for photon momentum calculation
  Double_t vertex[]  = {0,0,0} ; //vertex 
  Double_t vertex2[] = {0,0,0} ; //vertex from second aod input
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) 
  {
	  GetReader()->GetVertex(vertex);
	  if(GetReader()->GetSecondInputAODTree()) GetReader()->GetSecondInputAODVertex(vertex2);
  }
	
  //Select the Calorimeter of the photon
  if(fCalorimeter == "PHOS")
    pl = GetAODPHOS();
  else if (fCalorimeter == "EMCAL")
    pl = GetAODEMCAL();
  //Fill AODCaloClusters and AODParticle with PHOS aods
  TLorentzVector mom ;
  for(Int_t icalo = 0; icalo < pl->GetEntriesFast(); icalo++){
    AliAODCaloCluster * calo = (AliAODCaloCluster*) (pl->At(icalo));	
    
    //Cluster selection, not charged, with pi0 id and in fiducial cut
	  
	//Input from second AOD?
	Int_t input = 0;
	if     (fCalorimeter == "EMCAL" && GetReader()->GetAODEMCALNormalInputEntries() <= icalo) input = 1 ;
	else if(fCalorimeter == "PHOS"  && GetReader()->GetAODPHOSNormalInputEntries()  <= icalo) input = 1;
	  
	//Get Momentum vector, 
	if     (input == 0) calo->GetMomentum(mom,vertex) ;//Assume that come from vertex in straight line
	else if(input == 1) calo->GetMomentum(mom,vertex2);//Assume that come from vertex in straight line  
	  
	//If too small or big pt, skip it
    if(mom.Pt() < GetMinPt() || mom.Pt() > GetMaxPt() ) continue ; 
    //Check acceptance selection
    if(IsFiducialCutOn()){
      Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
      if(! in ) continue ;
    }
    
    //Create AOD for analysis
    AliAODPWG4Particle aodpi0 = AliAODPWG4Particle(mom);
    aodpi0.SetLabel(calo->GetLabel(0));
    //Set the indeces of the original caloclusters  
    aodpi0.SetCaloLabel(calo->GetID(),-1);
    aodpi0.SetDetector(fCalorimeter);
    if(GetDebug() > 1) 
      printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - FillAOD: Min pt cut and fiducial cut passed: pt %3.2f, phi %2.2f, eta %1.2f\n",aodpi0.Pt(),aodpi0.Phi(),aodpi0.Eta());	
    
    //Check Distance to Bad channel, set bit.
    Double_t distBad=calo->GetDistToBadChannel() ; //Distance to bad channel
    if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
    if(distBad < fMinDist) //In bad channel (PHOS cristal size 2.2x2.2 cm)
      continue ;
    
    if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - FillAOD: Bad channel cut passed %4.2f\n",distBad);
    
    if(distBad > fMinDist3) aodpi0.SetDistToBad(2) ;
    else if(distBad > fMinDist2) aodpi0.SetDistToBad(1) ; 
    else aodpi0.SetDistToBad(0) ;
    
    //Check PID
    //PID selection or bit setting
    if(GetReader()->GetDataType() == AliCaloTrackReader::kMC){
      //Get most probable PID, check PID weights (in MC this option is mandatory)
      aodpi0.SetPdg(GetCaloPID()->GetPdg(fCalorimeter,calo->PID(),mom.E()));//PID with weights
      if(GetDebug() > 1) 
	printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - FillAOD: PDG of identified particle %d\n",aodpi0.GetPdg());
      //If primary is not pi0, skip it.
      if(aodpi0.GetPdg() != AliCaloPID::kPi0) continue ;
    }					
    else if(IsCaloPIDOn()){
      //Skip matched clusters with tracks
      if(calo->GetNTracksMatched() > 0) continue ;
      
      //Get most probable PID, 2 options check PID weights 
      //or redo PID, recommended option for EMCal.		
      if(!IsCaloPIDRecalculationOn())
	aodpi0.SetPdg(GetCaloPID()->GetPdg(fCalorimeter,calo->PID(),mom.E()));//PID with weights
      else
	aodpi0.SetPdg(GetCaloPID()->GetPdg(fCalorimeter,mom,calo));//PID recalculated
      
      if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - PDG of identified particle %d\n",aodpi0.GetPdg());
      
      //If cluster does not pass pid, not pi0, skip it.
      if(aodpi0.GetPdg() != AliCaloPID::kPi0) continue ;			
      
    }
    else{
      //Set PID bits for later selection (AliAnaPi0 for example)
      //GetPDG already called in SetPIDBits.
      GetCaloPID()->SetPIDBits(fCalorimeter,calo,&aodpi0);
      if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - PID Bits set \n");		
    }
    
    if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - Pi0 selection cuts passed: pT %3.2f, pdg %d\n",aodpi0.Pt(), aodpi0.GetPdg());
    
    //Play with the MC stack if available
    //Check origin of the candidates
    if(IsDataMC()){
      if((GetReader()->GetDataType() == AliCaloTrackReader::kMC && fAnaType!=kSSCalo) || 
		 GetReader()->GetDataType() != AliCaloTrackReader::kMC){
		  aodpi0.SetInputFileIndex(input);
		  Int_t tag	=0;
		  tag = GetMCAnalysisUtils()->CheckOrigin(calo->GetLabel(0),GetReader(), aodpi0.GetInputFileIndex());
		  //GetMCAnalysisUtils()->CheckMultipleOrigin(calo->GetLabels(),calo->GetNLabel(), GetReader(), aodpi0.GetInputFileIndex(), tag);
		  aodpi0.SetTag(tag);
		  if(GetDebug() > 0) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - Origin of candidate %d\n",aodpi0.GetTag());
      }
    }//Work with stack also   
    
    //Add AOD with pi0 object to aod branch
    AddAODParticle(aodpi0);
    
  }//loop
  
  if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - End fill AODs \n");  
  
}
//__________________________________________________________________
void  AliAnaPi0EbE::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms
  
  if(!GetOutputAODBranch()){
    printf("AliAnaPi0EbE::MakeAnalysisFillHistograms()  - No output pi0 in AOD branch with name < %s >,STOP \n",GetOutputAODName().Data());
    abort();
  }
  //Loop on stored AOD pi0
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaPi0EbE::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    
    AliAODPWG4Particle* pi0 =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = pi0->GetPdg();
	  
    if(IsCaloPIDOn() && pdg != AliCaloPID::kPi0) continue;              
    
    //Fill pi0 histograms 
    Float_t pt  = pi0->Pt();
    Float_t phi = pi0->Phi();
    Float_t eta = pi0->Eta();
    
    fhPtPi0  ->Fill(pt);
    fhPhiPi0 ->Fill(pt,phi);
    fhEtaPi0 ->Fill(pt,eta);
    
    if(IsDataMC()){
      if((GetReader()->GetDataType() == AliCaloTrackReader::kMC && fAnaType!=kSSCalo) || 
	 GetReader()->GetDataType() != AliCaloTrackReader::kMC){
	if(GetMCAnalysisUtils()->CheckTagBit(pi0->GetTag(), AliMCAnalysisUtils::kMCPi0)){
	  fhPtMCPi0  ->Fill(pt);
	  fhPhiMCPi0 ->Fill(pt,phi);
	  fhEtaMCPi0 ->Fill(pt,eta);
	}
	else{
	  fhPtMCNoPi0  ->Fill(pt);
	  fhPhiMCNoPi0 ->Fill(pt,phi);
	  fhEtaMCNoPi0 ->Fill(pt,eta);
	}
      }
    }//Histograms with MC
    
  }// aod loop
  
}


//____________________________________________________________________________
void AliAnaPi0EbE::Init()
{ 
  //Init
  //Do some checks
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD()){
    printf("AliAnaPi0EbE::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD()){
    printf("AliAnaPi0EbE::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
}

//____________________________________________________________________________
void AliAnaPi0EbE::InitParameters()
{
  //Initialize the parameters of the analysis.  
  AddToHistogramsName("AnaPi0EbE_");

  fInputAODGammaConvName = "gammaconv" ;
  fAnaType = kIMCalo ;
  fCalorimeter = "EMCAL" ;
  fMinDist  = 2.;
  fMinDist2 = 4.;
  fMinDist3 = 5.;
  
}

//__________________________________________________________________
void AliAnaPi0EbE::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print("");
  printf("Analysis Type = %d \n",  fAnaType) ;
  if(fAnaType == kSSCalo){     
    printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ;
    printf("Min Distance to Bad Channel   = %2.1f\n",fMinDist);
    printf("Min Distance to Bad Channel 2 = %2.1f\n",fMinDist2);
    printf("Min Distance to Bad Channel 3 = %2.1f\n",fMinDist3); 
  } 
  printf("    \n") ;
  
} 
