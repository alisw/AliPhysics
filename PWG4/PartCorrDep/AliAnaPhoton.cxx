/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes hereby granted      *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id: AliAnaPhoton.cxx 28688 2008-09-11 15:04:07Z gconesab $ */

//_________________________________________________________________________
//
// Class for the photon identification.
// Clusters from calorimeters are identified as photons
// and kept in the AOD. Few histograms produced.
//
// -- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TH2F.h>

// --- Analysis system --- 
#include "AliAnaPhoton.h" 
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"

ClassImp(AliAnaPhoton)
  
//____________________________________________________________________________
  AliAnaPhoton::AliAnaPhoton() : 
    AliAnaPartCorrBaseClass(), fCalorimeter(""), 
	fMinDist(0.),fMinDist2(0.),fMinDist3(0.),
	fhPtPhoton(0),fhPhiPhoton(0),fhEtaPhoton(0),
	//MC
    fhPtPrompt(0),fhPhiPrompt(0),fhEtaPrompt(0), 
    fhPtFragmentation(0),fhPhiFragmentation(0),fhEtaFragmentation(0), 
    fhPtPi0Decay(0),fhPhiPi0Decay(0),fhEtaPi0Decay(0), 
    fhPtOtherDecay(0),fhPhiOtherDecay(0),fhEtaOtherDecay(0), 
    fhPtConversion(0),fhPhiConversion(0),fhEtaConversion(0), 
    fhPtUnknown(0),fhPhiUnknown(0),fhEtaUnknown(0)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();

}

//____________________________________________________________________________
AliAnaPhoton::AliAnaPhoton(const AliAnaPhoton & g) : 
  AliAnaPartCorrBaseClass(g), fCalorimeter(g.fCalorimeter),
   fMinDist(g.fMinDist),fMinDist2(g.fMinDist2), fMinDist3(g.fMinDist3),
   fhPtPhoton(g.fhPtPhoton),fhPhiPhoton(g.fhPhiPhoton),fhEtaPhoton(g.fhEtaPhoton), 
  //MC
  fhPtPrompt(g.fhPtPrompt),fhPhiPrompt(g.fhPhiPrompt),fhEtaPrompt(g.fhEtaPrompt), 
  fhPtFragmentation(g.fhPtFragmentation),fhPhiFragmentation(g.fhPhiFragmentation),fhEtaFragmentation(g.fhEtaFragmentation), 
  fhPtPi0Decay(g.fhPtPi0Decay),fhPhiPi0Decay(g.fhPhiPi0Decay),fhEtaPi0Decay(g.fhEtaPi0Decay), 
  fhPtOtherDecay(g.fhPtOtherDecay),fhPhiOtherDecay(g.fhPhiOtherDecay),fhEtaOtherDecay(g.fhEtaOtherDecay), 
  fhPtConversion(g. fhPtConversion),fhPhiConversion(g.fhPhiConversion),fhEtaConversion(g.fhEtaConversion), 
  fhPtUnknown(g.fhPtUnknown),fhPhiUnknown(g.fhPhiUnknown),fhEtaUnknown(g.fhEtaUnknown) 
{
  // cpy ctor
  
}

//_________________________________________________________________________
AliAnaPhoton & AliAnaPhoton::operator = (const AliAnaPhoton & g)
{
  // assignment operator
  
  if(&g == this) return *this;

  fCalorimeter = g.fCalorimeter ;
  fMinDist  = g.fMinDist;
  fMinDist2 = g.fMinDist2;
  fMinDist3 = g.fMinDist3;
  
  fhPtPhoton = g.fhPtPhoton ; 
  fhPhiPhoton = g.fhPhiPhoton ;
  fhEtaPhoton = g.fhEtaPhoton ;
 
  fhPtPrompt = g.fhPtPrompt;
  fhPhiPrompt = g.fhPhiPrompt;
  fhEtaPrompt = g.fhEtaPrompt; 
  fhPtFragmentation = g.fhPtFragmentation;
  fhPhiFragmentation = g.fhPhiFragmentation;
  fhEtaFragmentation = g.fhEtaFragmentation; 
  fhPtPi0Decay = g.fhPtPi0Decay;
  fhPhiPi0Decay = g.fhPhiPi0Decay;
  fhEtaPi0Decay = g.fhEtaPi0Decay; 
  fhPtOtherDecay = g.fhPtOtherDecay;
  fhPhiOtherDecay = g.fhPhiOtherDecay;
  fhEtaOtherDecay = g.fhEtaOtherDecay; 
  fhPtConversion = g. fhPtConversion;
  fhPhiConversion = g.fhPhiConversion;
  fhEtaConversion = g.fhEtaConversion; 
  fhPtUnknown = g.fhPtUnknown;
  fhPhiUnknown = g.fhPhiUnknown;
  fhEtaUnknown = g.fhEtaUnknown; 

  return *this;
  
}

//____________________________________________________________________________
AliAnaPhoton::~AliAnaPhoton() 
{
  //dtor

}


//________________________________________________________________________
TList *  AliAnaPhoton::GetCreateOutputObjects()
{  
	// Create histograms to be saved in output file and 
	// store them in outputContainer
	TList * outputContainer = new TList() ; 
	outputContainer->SetName("PhotonHistos") ; 
	
	Int_t nptbins  = GetHistoNPtBins();
	Int_t nphibins = GetHistoNPhiBins();
	Int_t netabins = GetHistoNEtaBins();
	Float_t ptmax  = GetHistoPtMax();
	Float_t phimax = GetHistoPhiMax();
	Float_t etamax = GetHistoEtaMax();
	Float_t ptmin  = GetHistoPtMin();
	Float_t phimin = GetHistoPhiMin();
	Float_t etamin = GetHistoEtaMin();	

	//Histograms of highest Photon identified in Event
	fhPtPhoton  = new TH1F("hPtPhoton","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
	fhPtPhoton->SetYTitle("N");
	fhPtPhoton->SetXTitle("p_{T #gamma}(GeV/c)");
	outputContainer->Add(fhPtPhoton) ; 
	
	fhPhiPhoton  = new TH2F
    ("hPhiPhoton","#phi_{#gamma}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
	fhPhiPhoton->SetYTitle("#phi");
	fhPhiPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
	outputContainer->Add(fhPhiPhoton) ; 
	
	fhEtaPhoton  = new TH2F
    ("hEtaPhoton","#phi_{#gamma}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
	fhEtaPhoton->SetYTitle("#eta");
	fhEtaPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
	outputContainer->Add(fhEtaPhoton) ;
	
	if(IsDataMC()){
		
		fhPtPrompt  = new TH1F("hPtPrompt","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
		fhPtPrompt->SetYTitle("N");
		fhPtPrompt->SetXTitle("p_{T #gamma}(GeV/c)");
		outputContainer->Add(fhPtPrompt) ; 
		
		fhPhiPrompt  = new TH2F
		("hPhiPrompt","#phi_{#gamma}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
		fhPhiPrompt->SetYTitle("#phi");
		fhPhiPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
		outputContainer->Add(fhPhiPrompt) ; 
		
		fhEtaPrompt  = new TH2F
		("hEtaPrompt","#phi_{#gamma}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
		fhEtaPrompt->SetYTitle("#eta");
		fhEtaPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
		outputContainer->Add(fhEtaPrompt) ;
		
		fhPtFragmentation  = new TH1F("hPtFragmentation","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
		fhPtFragmentation->SetYTitle("N");
		fhPtFragmentation->SetXTitle("p_{T #gamma}(GeV/c)");
		outputContainer->Add(fhPtFragmentation) ; 
		
		fhPhiFragmentation  = new TH2F
		("hPhiFragmentation","#phi_{#gamma}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
		fhPhiFragmentation->SetYTitle("#phi");
		fhPhiFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
		outputContainer->Add(fhPhiFragmentation) ; 
		
		fhEtaFragmentation  = new TH2F
		("hEtaFragmentation","#phi_{#gamma}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
		fhEtaFragmentation->SetYTitle("#eta");
		fhEtaFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
		outputContainer->Add(fhEtaFragmentation) ;
		
		fhPtPi0Decay  = new TH1F("hPtPi0Decay","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
		fhPtPi0Decay->SetYTitle("N");
		fhPtPi0Decay->SetXTitle("p_{T #gamma}(GeV/c)");
		outputContainer->Add(fhPtPi0Decay) ; 
		
		fhPhiPi0Decay  = new TH2F
		("hPhiPi0Decay","#phi_{#gamma}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
		fhPhiPi0Decay->SetYTitle("#phi");
		fhPhiPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
		outputContainer->Add(fhPhiPi0Decay) ; 
		
		fhEtaPi0Decay  = new TH2F
		("hEtaPi0Decay","#phi_{#gamma}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
		fhEtaPi0Decay->SetYTitle("#eta");
		fhEtaPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
		outputContainer->Add(fhEtaPi0Decay) ;
		
		fhPtOtherDecay  = new TH1F("hPtOtherDecay","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
		fhPtOtherDecay->SetYTitle("N");
		fhPtOtherDecay->SetXTitle("p_{T #gamma}(GeV/c)");
		outputContainer->Add(fhPtOtherDecay) ; 
		
		fhPhiOtherDecay  = new TH2F
		("hPhiOtherDecay","#phi_{#gamma}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
		fhPhiOtherDecay->SetYTitle("#phi");
		fhPhiOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
		outputContainer->Add(fhPhiOtherDecay) ; 
		
		fhEtaOtherDecay  = new TH2F
		("hEtaOtherDecay","#phi_{#gamma}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
		fhEtaOtherDecay->SetYTitle("#eta");
		fhEtaOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
		outputContainer->Add(fhEtaOtherDecay) ;
		
		fhPtConversion  = new TH1F("hPtConversion","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
		fhPtConversion->SetYTitle("N");
		fhPtConversion->SetXTitle("p_{T #gamma}(GeV/c)");
		outputContainer->Add(fhPtConversion) ; 
		
		fhPhiConversion  = new TH2F
		("hPhiConversion","#phi_{#gamma}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
		fhPhiConversion->SetYTitle("#phi");
		fhPhiConversion->SetXTitle("p_{T #gamma} (GeV/c)");
		outputContainer->Add(fhPhiConversion) ; 
		
		fhEtaConversion  = new TH2F
		("hEtaConversion","#phi_{#gamma}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
		fhEtaConversion->SetYTitle("#eta");
		fhEtaConversion->SetXTitle("p_{T #gamma} (GeV/c)");
		outputContainer->Add(fhEtaConversion) ;
		
		fhPtUnknown  = new TH1F("hPtUnknown","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
		fhPtUnknown->SetYTitle("N");
		fhPtUnknown->SetXTitle("p_{T #gamma}(GeV/c)");
		outputContainer->Add(fhPtUnknown) ; 
		
		fhPhiUnknown  = new TH2F
		("hPhiUnknown","#phi_{#gamma}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
		fhPhiUnknown->SetYTitle("#phi");
		fhPhiUnknown->SetXTitle("p_{T #gamma} (GeV/c)");
		outputContainer->Add(fhPhiUnknown) ; 
		
		fhEtaUnknown  = new TH2F
		("hEtaUnknown","#phi_{#gamma}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
		fhEtaUnknown->SetYTitle("#eta");
		fhEtaUnknown->SetXTitle("p_{T #gamma} (GeV/c)");
		outputContainer->Add(fhEtaUnknown) ;
	}//Histos with MC
	
 	//Save parameters used for analysis
	TString parList ; //this will be list of parameters used for this analysis.
	char onePar[255] ;
	
	sprintf(onePar,"--- AliAnaPhoton ---\n") ;
	parList+=onePar ;	
	sprintf(onePar,"Calorimeter: %s\n",fCalorimeter.Data()) ;
	parList+=onePar ;
	sprintf(onePar,"fMinDist =%2.2f (Minimal distance to bad channel to accept cluster) \n",fMinDist) ;
	parList+=onePar ;
	sprintf(onePar,"fMinDist2=%2.2f (Cuts on Minimal distance to study acceptance evaluation) \n",fMinDist2) ;
	parList+=onePar ;
	sprintf(onePar,"fMinDist3=%2.2f (One more cut on distance used for acceptance-efficiency study) \n",fMinDist3) ;
	parList+=onePar ;
	
 	//Get parameters set in base class.
	parList += GetBaseParametersList() ;
	
	//Get parameters set in PID class.
	parList += GetCaloPID()->GetPIDParametersList() ;
	
 	//Get parameters set in FidutialCut class (not available yet)
	//parlist += GetCaloPID()->GetFidCutParametersList() 
	
 	TObjString *oString= new TObjString(parList) ;
	outputContainer->Add(oString);
	
	return outputContainer ;
	
}

//____________________________________________________________________________
void AliAnaPhoton::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  SetOutputAODClassName("AliAODPWG4Particle");
  SetOutputAODName("photons");
  fCalorimeter = "PHOS" ;
  fMinDist  = 2.;
  fMinDist2 = 4.;
  fMinDist3 = 5.;
  
}

//__________________________________________________________________
void  AliAnaPhoton::MakeAnalysisFillAOD() 
{
	//Do analysis and fill aods
	//Search for photons in fCalorimeter 
	
	TClonesArray * pl = new TClonesArray; 

	//Get vertex for photon momentum calculation
	Double_t vertex[]={0,0,0} ; //vertex ;
	if(!GetReader()->GetDataType()== AliCaloTrackReader::kMC) GetReader()->GetVertex(vertex);

	//Select the Calorimeter of the photon
	if(fCalorimeter == "PHOS")
		pl = GetAODPHOS();
	else if (fCalorimeter == "EMCAL")
		pl = GetAODEMCAL();
	//Fill AODCaloClusters and AODParticle with PHOS aods
	TLorentzVector mom ;
	for(Int_t icalo = 0; icalo < pl->GetEntriesFast(); icalo++){
		AliAODCaloCluster * calo =  (AliAODCaloCluster*) (pl->At(icalo));	

		//Cluster selection, not charged, with photon id and in fidutial cut
		//Get Momentum vector, 
		calo->GetMomentum(mom,vertex);//Assume that come from vertex in straight line
		//If too small or big pt, skip it
		if(mom.Pt() < GetMinPt() || mom.Pt() > GetMaxPt() ) continue ; 
		//Check acceptance selection
		if(IsFidutialCutOn()){
			Bool_t in = GetFidutialCut()->IsInFidutialCut(mom,fCalorimeter) ;
			if(! in ) continue ;
		}
		
		//Create AOD for analysis
		AliAODPWG4Particle aodph = AliAODPWG4Particle(mom);
		aodph.SetLabel(calo->GetLabel(0));
		//printf("Index %d, Id %d\n",icalo, calo->GetID());
		//Set the indeces of the original caloclusters  
		aodph.SetCaloLabel(calo->GetID(),-1);
		aodph.SetDetector(fCalorimeter);
		if(GetDebug() > 1) printf("AliAnaPhoton::FillAOD: Min pt cut and fidutial cut passed: pt %3.2f, phi %2.2f, eta %1.2f\n",aodph.Pt(),aodph.Phi(),aodph.Eta());	
			
		//Check Distance to Bad channel, set bit.
		Double_t distBad=calo->GetDistToBadChannel() ; //Distance to bad channel
		if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
		if(distBad < fMinDist) //In bad channel (PHOS cristal size 2.2x2.2 cm)
		continue ;
		
		if(GetDebug() > 1) printf("AliAnaPhoton::FillAOD: Bad channel cut passed %4.2f\n",distBad);
				
		if(distBad > fMinDist3) aodph.SetDistToBad(2) ;
		else if(distBad > fMinDist2) aodph.SetDistToBad(1) ; 
		else aodph.SetDistToBad(0) ;
		
		//Check PID
		//PID selection or bit setting
		if(GetReader()->GetDataType() == AliCaloTrackReader::kMC){
			//Get most probable PID, check PID weights (in MC this option is mandatory)
			aodph.SetPdg(GetCaloPID()->GetPdg(fCalorimeter,calo->PID(),mom.E()));//PID with weights
			if(GetDebug() > 1) printf("AliAnaPhoton::FillAOD: PDG of identified particle %d\n",aodph.GetPdg());
			//If primary is not photon, skip it.
			if(aodph.GetPdg() != AliCaloPID::kPhoton) continue ;
		}					
		else if(IsCaloPIDOn()){
			//Skip matched clusters with tracks
			if(calo->GetNTracksMatched() > 0) continue ;
		
			//Get most probable PID, 2 options check PID weights 
			//or redo PID, recommended option for EMCal.		
			if(!IsCaloPIDRecalculationOn())
				aodph.SetPdg(GetCaloPID()->GetPdg(fCalorimeter,calo->PID(),mom.E()));//PID with weights
			else
				aodph.SetPdg(GetCaloPID()->GetPdg(fCalorimeter,mom,calo));//PID recalculated
				
			if(GetDebug() > 1) printf("AliAnaPhoton::FillAOD: PDG of identified particle %d\n",aodph.GetPdg());
			
			//If cluster does not pass pid, not photon, skip it.
			if(aodph.GetPdg() != AliCaloPID::kPhoton) continue ;			
				
		}
		else{
			//Set PID bits for later selection (AliAnaPi0 for example)
			//GetPDG already called in SetPIDBits.
			GetCaloPID()->SetPIDBits(fCalorimeter,calo,&aodph);
			if(GetDebug() > 1) printf("AliAnaPhoton::FillAOD: PID Bits set \n");		
		}
		
		if(GetDebug() > 1) printf("AliAnaPhoton::FillAOD: Photon selection cuts passed: pT %3.2f, pdg %d\n",aodph.Pt(), aodph.GetPdg());
		
		//Play with the MC stack if available
		//Check origin of the candidates
		if(IsDataMC()){
			aodph.SetTag(GetCaloPID()->CheckOrigin(calo->GetLabel(0),GetMCStack()));
			if(GetDebug() > 0) printf("AliAnaPhoton::FillAOD: Origin of candidate %d\n",aodph.GetTag());
		}//Work with stack also   
		
		//Add AOD with photon object to aod branch
		AddAODParticle(aodph);
	
	}//loop

if(GetDebug() > 1) printf("AliAnaPhoton::FillAOD: End fill AODs \n");  

}

//__________________________________________________________________
void  AliAnaPhoton::MakeAnalysisFillHistograms() 
{
	//Do analysis and fill histograms

	//Get vertex for photon momentum calculation
	Double_t v[]={0,0,0} ; //vertex ;
	//if(!GetReader()->GetDataType()== AliCaloTrackReader::kMC) 
	GetReader()->GetVertex(v);

	//Loop on stored AOD photons
	Int_t naod = GetOutputAODBranch()->GetEntriesFast();
	if(GetDebug() > 0) printf("AliAnaPhoton::FillHistos: aod branch entries %d\n", naod);

	for(Int_t iaod = 0; iaod < naod ; iaod++){
		AliAODPWG4Particle* ph =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
		Int_t pdg = ph->GetPdg();

		if(GetDebug() > 3) printf("AliAnaPhoton::FillHistos: PDG %d, MC TAG %d, Calorimeter %s\n", ph->GetPdg(),ph->GetTag(), (ph->GetDetector()).Data()) ;
		
		//If PID used, fill histos with photons in Calorimeter fCalorimeter
		if(IsCaloPIDOn())  
			if(pdg != AliCaloPID::kPhoton) continue; 
		if(ph->GetDetector() != fCalorimeter) continue;
		
		if(GetDebug() > 1) printf("AliAnaPhoton::FillHistos: ID Photon: pt %f, phi %f, eta %f\n", ph->Pt(),ph->Phi(),ph->Eta()) ;
		
		//Fill photon histograms 
		Float_t ptcluster = ph->Pt();
		Float_t phicluster = ph->Phi();
		Float_t etacluster = ph->Eta();
		
		fhPtPhoton   ->Fill(ptcluster);
		fhPhiPhoton ->Fill(ptcluster,phicluster);
		fhEtaPhoton ->Fill(ptcluster,etacluster);
		
		if(IsDataMC()){
			Int_t tag =ph->GetTag();
			if(tag == AliCaloPID::kMCPrompt){
				fhPtPrompt   ->Fill(ptcluster);
				fhPhiPrompt ->Fill(ptcluster,phicluster);
				fhEtaPrompt ->Fill(ptcluster,etacluster);
			}
			else if(tag==AliCaloPID::kMCFragmentation)
			{
				fhPtFragmentation   ->Fill(ptcluster);
				fhPhiFragmentation ->Fill(ptcluster,phicluster);
				fhEtaFragmentation ->Fill(ptcluster,etacluster);
			}
			else if(tag==AliCaloPID::kMCPi0Decay)
			{
				fhPtPi0Decay   ->Fill(ptcluster);
				fhPhiPi0Decay ->Fill(ptcluster,phicluster);
				fhEtaPi0Decay ->Fill(ptcluster,etacluster);
			}
			else if(tag==AliCaloPID::kMCEtaDecay || tag==AliCaloPID::kMCOtherDecay)
			{
				fhPtOtherDecay   ->Fill(ptcluster);
				fhPhiOtherDecay ->Fill(ptcluster,phicluster);
				fhEtaOtherDecay ->Fill(ptcluster,etacluster);
			}
			else if(tag==AliCaloPID::kMCConversion)
			{
				fhPtConversion   ->Fill(ptcluster);
				fhPhiConversion ->Fill(ptcluster,phicluster);
				fhEtaConversion ->Fill(ptcluster,etacluster);
			}
			else{
				
				fhPtUnknown   ->Fill(ptcluster);
				fhPhiUnknown ->Fill(ptcluster,phicluster);
				fhEtaUnknown ->Fill(ptcluster,etacluster);
			}
		}//Histograms with MC
		
	}// aod loop

}


//__________________________________________________________________
void AliAnaPhoton::Print(const Option_t * opt) const
{
	//Print some relevant parameters set for the analysis
	
	if(! opt)
		return;

	printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
	AliAnaPartCorrBaseClass::Print(" ");
    printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ;
	printf("Min Distance to Bad Channel   = %2.1f\n",fMinDist);
	printf("Min Distance to Bad Channel 2 = %2.1f\n",fMinDist2);
	printf("Min Distance to Bad Channel 3 = %2.1f\n",fMinDist3);
	printf("    \n") ;
	
} 
