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
 * Revision 1.4  2007/10/29 13:48:42  gustavo
 * Corrected coding violations
 *
 * Revision 1.2  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Base class for gamma and correlation analysis
// It is called by the task class AliAnalysisGammaTask and it connects the input (ESD/AOD/MonteCarlo)
// got with AliGammaReader (produces TClonesArrays of TParticles), with the analysis classes 
// AliAnaGammaDirect, AliAnaGammaCorrelation ....
//
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---

#include <TParticle.h>
#include <TH2.h>

//---- AliRoot system ---- 
#include "AliAnaGamma.h" 
#include "AliGammaReader.h" 
#include "AliAnaGammaDirect.h" 
#include "AliAnaGammaCorrelation.h" 
#include "AliAnaGammaSelection.h" 
#include "AliNeutralMesonSelection.h"
#include "AliAODCaloCluster.h"
#include "AliAODEvent.h"
#include "Riostream.h"
#include "AliLog.h"

ClassImp(AliAnaGamma)


//____________________________________________________________________________
  AliAnaGamma::AliAnaGamma() : 
    TObject(),
    fOutputContainer(0x0), 
    fAnaType(0),  fCalorimeter(0), fData(0x0), fKine(0x0), 
    fReader(0x0), fGammaDirect(0x0), fGammaCorrelation(0x0), fGammaSelection(0x0),
    fNeutralMesonSelection(0x0), fAODclusters(0x0), fNAODclusters(0)
{
  //Default Ctor

  //Initialize parameters, pointers and histograms
  if(!fReader)
    fReader = new AliGammaReader();
  if(!fGammaDirect)
    fGammaDirect = new AliAnaGammaDirect();
  if(!fGammaCorrelation)
    fGammaCorrelation = new AliAnaGammaCorrelation();
  if(!fGammaSelection)
    fGammaSelection = new AliAnaGammaSelection();
  if(!fNeutralMesonSelection)
    fNeutralMesonSelection = new AliNeutralMesonSelection();

  fAODclusters = 0;

  InitParameters();
  
}

//____________________________________________________________________________
AliAnaGamma::AliAnaGamma(const AliAnaGamma & g) :   
  TObject(),
  fOutputContainer(g. fOutputContainer), 
  fAnaType(g.fAnaType),  fCalorimeter(g.fCalorimeter), 
  fData(g.fData), fKine(g.fKine),fReader(g.fReader),
  fGammaDirect(g.fGammaDirect), fGammaCorrelation(g.fGammaCorrelation),
  fGammaSelection(g.fGammaSelection),
  fNeutralMesonSelection(g.fNeutralMesonSelection),  
  fAODclusters(g. fAODclusters), fNAODclusters(g.fNAODclusters)
{
  // cpy ctor
  
}

//_________________________________________________________________________
AliAnaGamma & AliAnaGamma::operator = (const AliAnaGamma & source)
{
  // assignment operator

  if(this == &source)return *this;
  ((TObject *)this)->operator=(source);

  fOutputContainer = source.fOutputContainer ;
  fAnaType = source.fAnaType;
  fCalorimeter = source.fCalorimeter ;
  fData = source.fData ; 
  fKine = source.fKine ;
  fReader = source.fReader ;
  fGammaDirect = source.fGammaDirect ;
  fGammaCorrelation = source.fGammaCorrelation ;
  fGammaSelection = source.fGammaSelection ;
  fNeutralMesonSelection = source.fNeutralMesonSelection ;

  return *this;

}

//____________________________________________________________________________
AliAnaGamma::~AliAnaGamma() 
{
  // Remove all pointers.

  // Protection added in case of NULL pointers (MG)
  if (fOutputContainer) {
     fOutputContainer->Clear();
     delete fOutputContainer ;
  }   

  if (fData) delete fData ; 
  if (fKine) delete fKine ;
  if (fReader) delete fReader ;
  if (fGammaDirect) delete fGammaDirect ;
  if (fGammaCorrelation) delete fGammaCorrelation ;
  if (fGammaSelection) delete fGammaSelection ;
  if (fNeutralMesonSelection) delete fNeutralMesonSelection ;

}

//________________________________________________________________________
void AliAnaGamma::Init()
{  

  //Init container histograms and other common variables

  //Histograms container
  fOutputContainer = new TList ;
  
  //Fill container with appropriate histograms
  
  //Set selection  analysis histograms
  TList * selectcontainer =  fGammaSelection->GetCreateOutputObjects(); 
  for(Int_t i = 0; i < selectcontainer->GetEntries(); i++){
    Bool_t  add = kTRUE ;
    TString name = (selectcontainer->At(i))->GetName();   
    if(!fReader->IsEMCALOn() && name.Contains("EMCAL")) add = kFALSE;
    if(!fReader->IsPHOSOn() && name.Contains("PHOS"))   add = kFALSE;
    if(!fReader->IsCTSOn() &&  !fGammaSelection->FillCTS() && name.Contains("CTS"))   add = kFALSE;
    if(add) fOutputContainer->Add(selectcontainer->At(i)) ;
  }  //Set selection  analysis histograms
 
  
  //Set prompt photon analysis histograms
  TList * promptcontainer =  fGammaDirect->GetCreateOutputObjects(); 
  for(Int_t i = 0; i < promptcontainer->GetEntries(); i++)
    fOutputContainer->Add(promptcontainer->At(i)) ;
  
  //Check if selected options are correct or set them when necessary
  if(fReader->GetDataType() == AliGammaReader::kMCData){
    fGammaDirect->SetMC();//Only useful with AliGammaMCDataReader, by default kFALSE
    fGammaSelection->SetMC();//Only useful with AliGammaMCDataReader, by default kFALSE
  }
  
  if(fAnaType == kCorrelation){
    
    //Check if selected options are correct
    if (fGammaDirect->GetICMethod()==AliAnaGammaDirect::kSeveralIC)
      AliFatal("Correlation not allowed with multiple isolation cuts, kCorrelation and kSeveralIC do not go together");
    
    if(fGammaCorrelation->GetCorrelationType() == AliAnaGammaCorrelation::kParton && 
       fReader->GetDataType() != AliGammaReader::kMC)
      AliFatal("kParton must be analyzed with data kMC");
    
    //Set the parameters for the neutral pair selection depending on the analysis, 
    fNeutralMesonSelection->SetDeltaPhiCutRange(fGammaCorrelation->GetDeltaPhiMinCut(), 
						fGammaCorrelation->GetDeltaPhiMaxCut());  
    

    if(fGammaCorrelation->GetCorrelationType() == AliAnaGammaCorrelation::kHadron){
      fNeutralMesonSelection->SetPhiPtSelection(AliNeutralMesonSelection::kSelectPhiMinPt);
      fNeutralMesonSelection->SetMinPt(fGammaCorrelation->GetMinPtHadron());
      
    }
    
    if(fGammaCorrelation->GetCorrelationType() == AliAnaGammaCorrelation::kJetLeadCone){
      fNeutralMesonSelection->SetPhiPtSelection(AliNeutralMesonSelection::kSelectPhiPtRatio);
      fNeutralMesonSelection->SetRatioCutRange(fGammaCorrelation->GetRatioMinCut(), 
					       fGammaCorrelation->GetRatioMaxCut());
    }

    //Set the neutral mesosn selection histograms
    TList * neutralmesoncontainer =  fNeutralMesonSelection->GetCreateOutputObjects();
    if(fNeutralMesonSelection->AreNeutralMesonSelectionHistosKept()){
      for(Int_t i = 0; i < neutralmesoncontainer->GetEntries(); i++)
	fOutputContainer->Add(neutralmesoncontainer->At(i)) ;
    }
    
    //Set correlation histograms
    TList * correlationcontainer =  fGammaCorrelation->GetCreateOutputObjects();
    for(Int_t i = 0; i < correlationcontainer->GetEntries(); i++)
      fOutputContainer->Add(correlationcontainer->At(i)) ;
    fGammaCorrelation->SetOutputContainer(fOutputContainer);
    fGammaCorrelation->SetNeutralMesonSelection(fNeutralMesonSelection);
 
  }//kCorrelation  
  
}

//____________________________________________________________________________
void AliAnaGamma::InitParameters()
{

  //Init data members
 
  fAnaType = kPrompt;
  fCalorimeter = "EMCAL";

}

//__________________________________________________________________
void AliAnaGamma::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  printf("Analysis type           =     %d\n", fAnaType) ;
  printf("Calorimeter           =     %s\n", fCalorimeter.Data()) ;

  switch(fAnaType)
    {
    case kPrompt:
      {
	fGammaDirect->Print("");
      }// case kIsolationCut
      break;
   
    case kCorrelation:
      {
	fGammaCorrelation->Print("");
      }//  case kCorrelation
      break;
      
    }//switch
} 


//____________________________________________________________________________
Bool_t AliAnaGamma::ProcessEvent(Long64_t entry){

  AliDebug(1,Form("Entry %d",entry));
  //cout<<"Event >>>>>>>>>>>>> "<<entry<<endl;
  if(!fOutputContainer)
    AliFatal("Histograms not initialized");

  //CreateTLists with arrays of TParticles. Filled with particles only relevant for the analysis.

  TClonesArray * plCTS      = new TClonesArray("TParticle",1000); // All particles refitted in Central Tracking System (ITS+TPC)
  TClonesArray * plEMCAL    = new TClonesArray("TParticle",1000);   // All particles measured in Jet Calorimeter (EMCAL)
  TClonesArray * plPHOS     = new TClonesArray("TParticle",1000);  // All particles measured  Gamma calorimeter
  TClonesArray * plPrimCTS      = new TClonesArray("TParticle",1000); // primary tracks
  TClonesArray * plPrimEMCAL    = new TClonesArray("TParticle",1000);   // primary emcal clusters
  TClonesArray * plPrimPHOS     = new TClonesArray("TParticle",1000);  // primary phos clusters
  TClonesArray * plParton   = new TClonesArray("TParticle",1000);  // All partons
  //Fill lists with photons, neutral particles and charged particles
  //look for the highest energy photon in the event inside fCalorimeter
    
  //Fill particle lists 
  if(fReader->GetDataType() == AliGammaReader::kData){
    AliDebug(1,"Data analysis");
    fReader->CreateParticleList(fData, NULL,plCTS,plEMCAL,plPHOS,NULL,NULL,NULL); 
  }
  else if( fReader->GetDataType()== AliGammaReader::kMC){
    AliDebug(1,"Kinematics analysis");
    fReader->CreateParticleList(fKine, NULL,plCTS,plEMCAL,plPHOS,plParton,NULL,NULL); 
  }
  else if(fReader->GetDataType() == AliGammaReader::kMCData) {
   AliDebug(1,"Data + Kinematics analysis");
   fReader->CreateParticleList(fData, fKine,plCTS,plEMCAL,plPHOS,plPrimCTS,plPrimEMCAL,plPrimPHOS); 
  }
  else
    AliError("Option not implemented");

  //Fill AOD with calorimeter
  //Temporal solution, just for testing
  //FillAODs(plPHOS,plEMCAL);  

  //Select particles to do the final analysis.
  if(fReader->IsEMCALOn()) fGammaSelection->Selection("EMCAL",plEMCAL,plPrimEMCAL);
  if(fReader->IsPHOSOn())  fGammaSelection->Selection("PHOS",plPHOS,plPrimPHOS);
  if(fReader->IsCTSOn() &&  fGammaSelection->FillCTS()) fGammaSelection->Selection("CTS",plCTS,plPrimCTS);


  //Search highest energy prompt gamma in calorimeter
  if(fCalorimeter == "PHOS")
    MakeAnalysis(plPHOS, plEMCAL, plCTS, plParton, plPrimPHOS) ; 
  else if (fCalorimeter == "EMCAL")
    MakeAnalysis(plEMCAL, plPHOS, plCTS,plParton, plPrimEMCAL) ; 
  else
    AliFatal("Wrong calorimeter name");

  plCTS->Clear() ;
  plEMCAL->Clear() ;
  plPHOS->Clear() ;
  plParton->Clear() ;
  plPrimCTS->Clear() ;
  plPrimEMCAL->Clear() ;
  plPrimPHOS->Clear() ;

  delete plCTS ;
  delete plPHOS ;
  delete plEMCAL ;
  delete plParton ;
  delete plPrimCTS ;
  delete plPrimPHOS ;
  delete plPrimEMCAL ;

  return kTRUE;

}

//____________________________________________________________________________
void AliAnaGamma::MakeAnalysis(TClonesArray * plCalo, TClonesArray * plNe, TClonesArray * plCTS, TClonesArray * plParton, TClonesArray * plPrimCalo)  {
  
  TParticle * pGamma = new TParticle ;
  Bool_t isInCalo = kFALSE ;
  
  switch(fAnaType)
    {
      
      //Only Prompt photon analysis
    case kPrompt:
      {	  
	AliDebug(1,"kPrompt analysis");
	switch(fGammaDirect->GetICMethod())
	  {
	    
	  case AliAnaGammaDirect::kSeveralIC:
	    {
	      fGammaDirect->MakeSeveralICAnalysis(plCalo, plCTS);
	      AliDebug(1,"kSeveralIC analysis");
	    }
	    break;
	    
	  default :
	    {
	      fGammaDirect->GetPromptGamma(plCalo, plCTS,plPrimCalo, pGamma,isInCalo);
	      if(!isInCalo)
		AliDebug(1,"Prompt gamma not found");
	    }
	    break;
	  }//IC method
      }// case kPrompt:
      break;
      
      //Correlate prompt photon with something: parton, hadron, jet.
    case kCorrelation:
      {
	AliDebug(1,"kCorrelation analysis");
	//Find prompt photon	
	fGammaDirect->GetPromptGamma(plCalo, plCTS,plPrimCalo, pGamma,isInCalo);

	if(isInCalo){//If prompt photon found, do correlation
	  
	  switch(fGammaCorrelation->GetCorrelationType())
	    {
	    case AliAnaGammaCorrelation::kParton:
	      {
		AliDebug(1,"kParton correlation");
		fGammaCorrelation->MakeGammaCorrelation(pGamma, plParton, NULL);
	      }//  case kParton
	      break;
      
	    case AliAnaGammaCorrelation::kHadron:
	      {
		AliDebug(1,"kHadron correlation");
		fGammaCorrelation->MakeGammaCorrelation(pGamma, plCTS, plNe);
	      }//  case kHadron
	      break;

	    case AliAnaGammaCorrelation::kJetLeadCone:
	      {		
		AliDebug(1,"kJetLeadCone correlation");
		fGammaCorrelation->MakeGammaCorrelation(pGamma, plCTS, plNe);
	      }//  case kJetLeadCone
	      break;
	      
	    case AliAnaGammaCorrelation::kJetFinder:
	      {	
		AliDebug(1,"kJetFinder correlation");
		printf("Analysis not implemented \n");
	      }//  case kJetFinder
	      break;
	    }// switch correlation
	}// is in calo
	else  AliDebug(2,"Prompt gamma not found");
      }// case kCorrelation
      break;

    } //switch(fAnaType)

  delete pGamma ; 
  
}

//____________________________________________________
void AliAnaGamma::AddCluster(AliAODCaloCluster p)
{
// Add new jet to the list
  cout<<"AOD list pointer "<<fAODclusters<<" nAODclusters "<<fNAODclusters<<endl;
  cout<<"list entries "<<fAODclusters->GetEntries()<<endl;
  new ((*fAODclusters)[fNAODclusters++]) AliAODCaloCluster(p);
}

//___________________________________________________
void AliAnaGamma::ConnectAOD(AliAODEvent* aod)
{
// Connect to the AOD
  fAODclusters = aod->GetCaloClusters();
}

//____________________________________________________________________________
void AliAnaGamma::FillAODs(TClonesArray * plPHOS, TClonesArray * plEMCAL){
  //Fill AOD caloClusters
  //Temporal method, just for testing AOD creation
  
  //Fill AOD with PHOS clusters
  Int_t nphos =  plPHOS->GetEntries() ;
  cout<<"PHOS entries "<<nphos<<endl;
  for(Int_t ipr = 0;ipr < nphos ; ipr ++ ){
    TParticle * particle = dynamic_cast<TParticle *>(plPHOS->At(ipr)) ;
    Float_t pos []= {0,0,0};
    AliAODCaloCluster phos(ipr,0,0x0,particle->Energy(),pos,0x0,AliAODCluster::kPHOSNeutral,0);
    AddCluster(phos);
  }
  
  //Fill AOD with EMCAL clusters
  cout<<"EMCAL entries "<<plEMCAL->GetEntries()<<endl;
  for(Int_t ipr = 0;ipr < plEMCAL->GetEntries() ; ipr ++ ){
    TParticle * particle = dynamic_cast<TParticle *>(plEMCAL->At(ipr)) ;
    Float_t pos []= {0,0,0};
    AliAODCaloCluster emcal(ipr+nphos,0,0x0,particle->Energy(),pos,0x0,AliAODCluster::kEMCALClusterv1,0);
    AddCluster(emcal);
  }
}
