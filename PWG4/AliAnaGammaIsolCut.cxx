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
 * Revision 1.1  2007/01/23 17:17:29  schutz
 * New Gamma package
 *
 *
 */

//_________________________________________________________________________
// Class for the analysis of gamma correlations (gamma-jet, 
// gamma-hadron and isolation cut.
// This class makes isolation cut analysis for 2 IC methods 
//(cone pt sum and particle pt threshold), for different cone sizes 
//and pt thresholds
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

#include "AliAnaGammaIsolCut.h" 
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "Riostream.h"
#include "AliLog.h"

ClassImp(AliAnaGammaIsolCut)

//____________________________________________________________________________
  AliAnaGammaIsolCut::AliAnaGammaIsolCut(const char *name) : 
    AliAnaGammaDirect(name),
    fOutputContainer(new TObjArray(100)),  
    fNCones(0),fNPtThres(0)
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
AliAnaGammaIsolCut::AliAnaGammaIsolCut(const AliAnaGammaIsolCut & ic) : 
  AliAnaGammaDirect(ic),
  fOutputContainer(ic. fOutputContainer), 
  fNCones(ic.fNCones),fNPtThres(ic.fNPtThres)
{
  // cpy ctor
  SetName (ic.GetName()) ; 
  SetTitle(ic.GetTitle()) ; 

  for(Int_t i = 0; i < 10 ; i++){
    fConeSizes[i]=  ic.fConeSizes[i];
    fPtThresholds[i]=   ic.fPtThresholds[i];
  }
}

//____________________________________________________________________________
AliAnaGammaIsolCut::~AliAnaGammaIsolCut() 
{
  // Remove all pointers
  fOutputContainer->Clear() ; 
  delete fOutputContainer ;

  delete fhPtCandidate ;
  delete [] fhPtThresIsolated ;
  delete [] fhPtSumIsolated ;

}



//____________________________________________________________________________
void AliAnaGammaIsolCut::Exec(Option_t *) 
{
  
  // Processing of one event
    
  //Get ESDs
  Long64_t entry = GetChain()->GetReadEntry() ;
  
  if (!GetESD()) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if (GetPrintInfo()) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(GetChain()))->GetFile()->GetName(), entry)) ; 

  //CreateTLists with arrays of TParticles. Filled with particles only relevant for the analysis.

  TClonesArray * particleList = new TClonesArray("TParticle",1000); // All particles refitted in CTS and detected in EMCAL (jet)
  TClonesArray * plCTS         = new TClonesArray("TParticle",1000); // All particles refitted in Central Tracking System (ITS+TPC)
  TClonesArray * plNe         = new TClonesArray("TParticle",1000);   // All particles measured in Jet Calorimeter (EMCAL)
  TClonesArray * plPHOS     = new TClonesArray("TParticle",1000);  // All particles measured in PHOS as Gamma calorimeter
  TClonesArray * plEMCAL   = new TClonesArray("TParticle",1000);  // All particles measured in EMCAL as Gamma calorimeter
  
  
  AliDebug(2, "Fill particle lists");
  
  //Fill particle lists 
  CreateParticleList(particleList, plCTS,plEMCAL,plPHOS); 
  
  if(GetCalorimeter() == "PHOS")
    plNe = plPHOS;
  if(GetCalorimeter() == "EMCAL")
    plNe = plEMCAL;
  
  //Isolation Cut Analysis for both methods and different pt cuts and cones
  
  for(Int_t ipr = 0; ipr < plNe->GetEntries() ; ipr ++ ){
    TParticle * pCandidate = dynamic_cast<TParticle *>(plNe->At(ipr)) ;
    
    if(pCandidate->Pt() > GetMinGammaPt()){
      
      Bool_t  icPtThres   = kFALSE;
      Bool_t  icPtSum     = kFALSE;
      
      Float_t ptC             = pCandidate->Pt() ;
      
      fhPtCandidate->Fill(ptC);
      
      for(Int_t icone = 0; icone<fNCones; icone++){
	SetConeSize(fConeSizes[icone]) ;
	Float_t coneptsum = 0 ;
	for(Int_t ipt = 0; ipt<fNPtThres;ipt++){ 
	  SetPtThreshold(fPtThresholds[ipt]) ;
	  MakeIsolationCut(plCTS,plNe, pCandidate, ipr, icPtThres, icPtSum,coneptsum);
	  AliDebug(4,Form("Candidate pt %f, pt in cone %f, Isolated? ICPt %d, ICSum %d",
			  pCandidate->Pt(), coneptsum, icPtThres, icPtSum));
	  
	  fhPtThresIsolated[icone][ipt]->Fill(ptC); 
	}//pt thresh loop
	fhPtSumIsolated[icone]->Fill(ptC,coneptsum) ;
      }//cone size loop
    }//min pt candidate
  }//candidate loop
  
  AliDebug(2, "End of analysis, delete pointers");
  
  particleList->Delete() ; 
  plCTS->Delete() ;
  plNe->Delete() ;
  plPHOS->Delete() ;
  plEMCAL->Delete() ;

  delete plNe ;
  delete plCTS ;
  //delete plPHOS ;
  //delete plEMCAL ;
  delete particleList ;

  PostData(0, fOutputContainer);
}    

  //____________________________________________________________________________
void AliAnaGammaIsolCut::Init(const Option_t * )
{
  // Initialisation of branch container 
  AliAnaGammaDirect::Init();

  fNCones           = 4 ; 
  fNPtThres         = 4 ; 
  fConeSizes[0] = 0.1; fConeSizes[0] = 0.2; fConeSizes[2] = 0.3; fConeSizes[3] = 0.4;
  fPtThresholds[0]=1.; fPtThresholds[0]=2.; fPtThresholds[0]=3.; fPtThresholds[0]=4.;

  //Initialization of histograms 
  MakeHistos() ;
}

//___________________________________________________________________
void AliAnaGammaIsolCut::MakeHistos()
{
  // Create histograms to be saved in output file and 
  // stores them in fOutputContainer
  
  fOutputContainer = new TObjArray(10000) ; 

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


void AliAnaGammaIsolCut::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;

  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  printf("N Cone Sizes               =     %d\n", fNCones) ; 
  printf("N pT thresholds           =     %d\n", fNPtThres) ;
  printf("Cone Sizes                  =    \n") ;
  for(Int_t i = 0; i < fNCones; i++)
    printf("   %f;",  fConeSizes[i]) ;
  printf("    \n") ;
  for(Int_t i = 0; i < fNPtThres; i++)
    printf("   %f;",  fPtThresholds[i]) ;
  printf("    \n") ;
  
} 

void AliAnaGammaIsolCut::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
 
}
