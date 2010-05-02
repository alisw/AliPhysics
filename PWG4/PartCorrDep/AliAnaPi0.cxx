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

//_________________________________________________________________________
// Class to collect two-photon invariant mass distributions for
// extractin raw pi0 yield.
//
//-- Author: Dmitri Peressounko (RRC "KI") 
//-- Adapted to PartCorr frame by Lamia Benhabib (SUBATECH)
//-- and Gustavo Conesa (INFN-Frascati)
//_________________________________________________________________________


// --- ROOT system ---
#include "TH3.h"
#include "TH2D.h"
//#include "Riostream.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TROOT.h"
#include "TClonesArray.h"
//#include "TObjString.h"

//---- AliRoot system ----
#include "AliAnaPi0.h"
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"
#include "AliStack.h"
#include "AliFiducialCut.h"
#include "TParticle.h"
#include "AliAODCaloCluster.h"
#include "AliVEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliNeutralMesonSelection.h"

ClassImp(AliAnaPi0)

//________________________________________________________________________________________________________________________________________________  
AliAnaPi0::AliAnaPi0() : AliAnaPartCorrBaseClass(),
fNCentrBin(0),fNZvertBin(0),fNrpBin(0),
fNPID(0),fNmaxMixEv(0), fZvtxCut(0.),fCalorimeter(""),
fNModules(12), fUseAngleCut(kFALSE), fEventsList(0x0), //fhEtalon(0x0), 
fhReMod(0x0), fhRe1(0x0),fhMi1(0x0),fhRe2(0x0),fhMi2(0x0),fhRe3(0x0),fhMi3(0x0),fhEvents(0x0),
fhRealOpeningAngle(0x0),fhRealCosOpeningAngle(0x0),
fhPrimPt(0x0), fhPrimAccPt(0x0), fhPrimY(0x0), fhPrimAccY(0x0), fhPrimPhi(0x0), fhPrimAccPhi(0x0),
fhPrimOpeningAngle(0x0),fhPrimCosOpeningAngle(0x0)
{
//Default Ctor
 InitParameters();
 
}

//________________________________________________________________________________________________________________________________________________
AliAnaPi0::AliAnaPi0(const AliAnaPi0 & ex) : AliAnaPartCorrBaseClass(ex),  
fNCentrBin(ex.fNCentrBin),fNZvertBin(ex.fNZvertBin),fNrpBin(ex.fNrpBin),
fNPID(ex.fNPID),fNmaxMixEv(ex.fNmaxMixEv),fZvtxCut(ex.fZvtxCut), fCalorimeter(ex.fCalorimeter),
fNModules(ex.fNModules), fUseAngleCut(ex.fUseAngleCut), fEventsList(0x0), //fhEtalon(ex.fhEtalon), 
fhReMod(ex.fhReMod), fhRe1(ex.fhRe1),fhMi1(ex.fhMi1),fhRe2(ex.fhRe2),fhMi2(ex.fhMi2),
fhRe3(ex.fhRe3),fhMi3(ex.fhMi3),fhEvents(ex.fhEvents),
fhRealOpeningAngle(ex.fhRealOpeningAngle),fhRealCosOpeningAngle(ex.fhRealCosOpeningAngle),
fhPrimPt(ex.fhPrimPt), fhPrimAccPt(ex.fhPrimAccPt), fhPrimY(ex.fhPrimY), 
fhPrimAccY(ex.fhPrimAccY), fhPrimPhi(ex.fhPrimPhi), fhPrimAccPhi(ex.fhPrimAccPhi),
fhPrimOpeningAngle(ex.fhPrimOpeningAngle),fhPrimCosOpeningAngle(ex.fhPrimCosOpeningAngle)
{
  // cpy ctor
  //Do not need it
}
//
////________________________________________________________________________________________________________________________________________________
//AliAnaPi0 & AliAnaPi0::operator = (const AliAnaPi0 & ex)
//{
//  // assignment operator
//  
//  if(this == &ex)return *this;
//  ((AliAnaPartCorrBaseClass *)this)->operator=(ex);
//  
//  fNCentrBin = ex.fNCentrBin  ; fNZvertBin = ex.fNZvertBin  ; fNrpBin = ex.fNrpBin  ; 
//  fNPID = ex.fNPID  ; fNmaxMixEv = ex.fNmaxMixEv  ; fZvtxCut = ex.fZvtxCut  ;  fCalorimeter = ex.fCalorimeter  ;  
//  fNModules = ex.fNModules; fEventsList = ex.fEventsList  ;  //fhEtalon = ex.fhEtalon  ; 
//  fhRe1 = ex.fhRe1  ; fhMi1 = ex.fhMi1  ; fhRe2 = ex.fhRe2  ; fhMi2 = ex.fhMi2  ; fhReMod = ex.fhReMod; 
//  fhRe3 = ex.fhRe3  ; fhMi3 = ex.fhMi3  ; fhEvents = ex.fhEvents  ; fUseAngleCut = ex.fUseAngleCut;
//  fhPrimPt = ex.fhPrimPt  ;  fhPrimAccPt = ex.fhPrimAccPt  ;  fhPrimY = ex.fhPrimY  ;  
//  fhPrimAccY = ex.fhPrimAccY  ;  fhPrimPhi = ex.fhPrimPhi  ;  fhPrimAccPhi = ex.fhPrimAccPhi ;
//  fhRealOpeningAngle = ex.fhRealOpeningAngle; fhRealCosOpeningAngle = ex.fhRealCosOpeningAngle;
//  fhPrimOpeningAngle = ex.fhPrimOpeningAngle; fhPrimCosOpeningAngle = ex.fhPrimCosOpeningAngle;
//
//  return *this;
//  
//}

//________________________________________________________________________________________________________________________________________________
AliAnaPi0::~AliAnaPi0() {
  // Remove event containers
  if(fEventsList){
    for(Int_t ic=0; ic<fNCentrBin; ic++){
      for(Int_t iz=0; iz<fNZvertBin; iz++){
	for(Int_t irp=0; irp<fNrpBin; irp++){
	  fEventsList[ic*fNZvertBin*fNrpBin+iz*fNrpBin+irp]->Delete() ;
	  delete fEventsList[ic*fNZvertBin*fNrpBin+iz*fNrpBin+irp] ;
	}
      }
    }
    delete[] fEventsList; 
    fEventsList=0 ;
  }
	
}

//________________________________________________________________________________________________________________________________________________
void AliAnaPi0::InitParameters()
{
//Init parameters when first called the analysis
//Set default parameters
  SetInputAODName("PWG4Particle");
  
  AddToHistogramsName("AnaPi0_");
  fNModules = 12; // set maximum to maximum number of EMCAL modules
  fNCentrBin = 1;
  fNZvertBin = 1;
  fNrpBin    = 1;
  fNPID      = 9;
  fNmaxMixEv = 10;
  fZvtxCut   = 40;
  fCalorimeter  = "PHOS";
  fUseAngleCut = kFALSE;
	
}
//________________________________________________________________________________________________________________________________________________
//void AliAnaPi0::Init()
//{  
  //Init some data members needed in analysis
  
  //Histograms binning and range
 // if(!fhEtalon){                                                   //  p_T      alpha   d m_gg    
 //   fhEtalon = new TH3D("hEtalon","Histo with binning parameters",50,0.,25.,10,0.,1.,200,0.,1.) ; 
 //   fhEtalon->SetXTitle("P_{T} (GeV)") ;
 //   fhEtalon->SetYTitle("#alpha") ;
 //   fhEtalon->SetZTitle("m_{#gamma#gamma} (GeV)") ;
 // }
  
//}

//________________________________________________________________________________________________________________________________________________
TList * AliAnaPi0::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer
  
  //create event containers
  fEventsList = new TList*[fNCentrBin*fNZvertBin*fNrpBin] ;
	
  for(Int_t ic=0; ic<fNCentrBin; ic++){
    for(Int_t iz=0; iz<fNZvertBin; iz++){
      for(Int_t irp=0; irp<fNrpBin; irp++){
	fEventsList[ic*fNZvertBin*fNrpBin+iz*fNrpBin+irp] = new TList() ;
      }
    }
  }
  	
  TList * outputContainer = new TList() ; 
  outputContainer->SetName(GetName()); 
	
  fhReMod = new TH3D*[fNModules] ;
  fhRe1 = new TH3D*[fNCentrBin*fNPID] ;
  fhRe2 = new TH3D*[fNCentrBin*fNPID] ;
  fhRe3 = new TH3D*[fNCentrBin*fNPID] ;
  fhMi1 = new TH3D*[fNCentrBin*fNPID] ;
  fhMi2 = new TH3D*[fNCentrBin*fNPID] ;
  fhMi3 = new TH3D*[fNCentrBin*fNPID] ;
  
  char key[255] ;
  char title[255] ;
  Int_t nptbins   = GetHistoPtBins();
  Int_t nphibins  = GetHistoPhiBins();
  Int_t netabins  = GetHistoEtaBins();
  Float_t ptmax   = GetHistoPtMax();
  Float_t phimax  = GetHistoPhiMax();
  Float_t etamax  = GetHistoEtaMax();
  Float_t ptmin   = GetHistoPtMin();
  Float_t phimin  = GetHistoPhiMin();
  Float_t etamin  = GetHistoEtaMin();	
	
  Int_t nmassbins = GetHistoMassBins();
  Int_t nasymbins = GetHistoAsymmetryBins();
  Float_t massmax = GetHistoMassMax();
  Float_t asymmax = GetHistoAsymmetryMax();
  Float_t massmin = GetHistoMassMin();
  Float_t asymmin = GetHistoAsymmetryMin();
	
  for(Int_t ic=0; ic<fNCentrBin; ic++){
    for(Int_t ipid=0; ipid<fNPID; ipid++){
		
      //Distance to bad module 1
      sprintf(key,"hRe_cen%d_pid%d_dist1",ic,ipid) ;
      sprintf(title,"Real m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      
      //fhEtalon->Clone(key);
      //fhRe1[ic*fNPID+ipid]=(TH3D*)fhEtalon->Clone(key) ;
      //fhRe1[ic*fNPID+ipid]->SetName(key) ;
      //fhRe1[ic*fNPID+ipid]->SetTitle(title) ;
	  fhRe1[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhRe1[ic*fNPID+ipid]) ;
      
      sprintf(key,"hMi_cen%d_pid%d_dist1",ic,ipid) ;
      sprintf(title,"Mixed m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      //fhMi1[ic*fNPID+ipid]=(TH3D*)fhEtalon->Clone(key) ;
      //fhMi1[ic*fNPID+ipid]->SetName(key) ;
      //fhMi1[ic*fNPID+ipid]->SetTitle(title) ;
	  fhMi1[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhMi1[ic*fNPID+ipid]) ;
      
      //Distance to bad module 2
      sprintf(key,"hRe_cen%d_pid%d_dist2",ic,ipid) ;
      sprintf(title,"Real m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      //fhRe2[ic*fNPID+ipid]=(TH3D*)fhEtalon->Clone(key) ;
      //fhRe2[ic*fNPID+ipid]->SetName(key) ;
      //fhRe2[ic*fNPID+ipid]->SetTitle(title) ;
	  fhRe2[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhRe2[ic*fNPID+ipid]) ;
      
      sprintf(key,"hMi_cen%d_pid%d_dist2",ic,ipid) ;
      sprintf(title,"Mixed m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      //fhMi2[ic*fNPID+ipid]=(TH3D*)fhEtalon->Clone(key) ;
      //fhMi2[ic*fNPID+ipid]->SetName(key) ;
      //fhMi2[ic*fNPID+ipid]->SetTitle(title) ;
	  fhMi2[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhMi2[ic*fNPID+ipid]) ;
      
      //Distance to bad module 3
      sprintf(key,"hRe_cen%d_pid%d_dist3",ic,ipid) ;
      sprintf(title,"Real m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      //fhRe3[ic*fNPID+ipid]=(TH3D*)fhEtalon->Clone(key) ;
      //fhRe3[ic*fNPID+ipid]->SetName(key) ; 
      //fhRe3[ic*fNPID+ipid]->SetTitle(title) ;
	  fhRe3[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhRe3[ic*fNPID+ipid]) ;
      
      sprintf(key,"hMi_cen%d_pid%d_dist3",ic,ipid) ;
      sprintf(title,"Mixed m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      //fhMi3[ic*fNPID+ipid]=(TH3D*)fhEtalon->Clone(key) ;
      //fhMi3[ic*fNPID+ipid]->SetName(key) ;
      //fhMi3[ic*fNPID+ipid]->SetTitle(title) ;
	  fhMi3[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhMi3[ic*fNPID+ipid]) ;
    }
  }
  
  
  fhEvents=new TH3D("hEvents","Number of events",fNCentrBin,0.,1.*fNCentrBin,
		    fNZvertBin,0.,1.*fNZvertBin,fNrpBin,0.,1.*fNrpBin) ;
  outputContainer->Add(fhEvents) ;
  
	
  fhRealOpeningAngle  = new TH2D
  ("hRealOpeningAngle","Angle between all #gamma pair vs E_{#pi^{0}}",nptbins,ptmin,ptmax,200,0,0.5); 
  fhRealOpeningAngle->SetYTitle("#theta(rad)");
  fhRealOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
  outputContainer->Add(fhRealOpeningAngle) ;

  fhRealCosOpeningAngle  = new TH2D
  ("hRealCosOpeningAngle","Cosinus of angle between all #gamma pair vs E_{#pi^{0}}",nptbins,ptmin,ptmax,200,-1,1); 
  fhRealCosOpeningAngle->SetYTitle("cos (#theta) ");
  fhRealCosOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
  outputContainer->Add(fhRealCosOpeningAngle) ;
	
	
  //Histograms filled only if MC data is requested 	
  if(IsDataMC() || (GetReader()->GetDataType() == AliCaloTrackReader::kMC) ){
   // if(fhEtalon->GetXaxis()->GetXbins() && fhEtalon->GetXaxis()->GetXbins()->GetSize()){ //Variable bin size
   //   fhPrimPt = new TH1D("hPrimPt","Primary pi0 pt",fhEtalon->GetXaxis()->GetNbins(),fhEtalon->GetXaxis()->GetXbins()->GetArray()) ;
   //   fhPrimAccPt = new TH1D("hPrimAccPt","Primary pi0 pt with both photons in acceptance",fhEtalon->GetXaxis()->GetNbins(),
	//		     fhEtalon->GetXaxis()->GetXbins()->GetArray()) ;
   // }
   // else{
   //   fhPrimPt = new TH1D("hPrimPt","Primary pi0 pt",fhEtalon->GetXaxis()->GetNbins(),fhEtalon->GetXaxis()->GetXmin(),fhEtalon->GetXaxis()->GetXmax()) ;
   //   fhPrimAccPt = new TH1D("hPrimAccPt","Primary pi0 pt with both photons in acceptance",
	//		     fhEtalon->GetXaxis()->GetNbins(),fhEtalon->GetXaxis()->GetXmin(),fhEtalon->GetXaxis()->GetXmax()) ;
   // }

	fhPrimPt     = new TH1D("hPrimPt","Primary pi0 pt",nptbins,ptmin,ptmax) ;
    fhPrimAccPt  = new TH1D("hPrimAccPt","Primary pi0 pt with both photons in acceptance",nptbins,ptmin,ptmax) ;
    outputContainer->Add(fhPrimPt) ;
    outputContainer->Add(fhPrimAccPt) ;
    
    fhPrimY      = new TH1D("hPrimaryRapidity","Rapidity of primary pi0",netabins,etamin,etamax) ; 
    outputContainer->Add(fhPrimY) ;
    
    fhPrimAccY   = new TH1D("hPrimAccRapidity","Rapidity of primary pi0",netabins,etamin,etamax) ; 
    outputContainer->Add(fhPrimAccY) ;
    
	fhPrimPhi    = new TH1D("hPrimaryPhi","Azimithal of primary pi0",nphibins,phimin*TMath::RadToDeg(),phimax*TMath::RadToDeg()) ; 
    outputContainer->Add(fhPrimPhi) ;
    
    fhPrimAccPhi = new TH1D("hPrimAccPhi","Azimithal of primary pi0 with accepted daughters",nphibins,phimin*TMath::RadToDeg(),phimax*TMath::RadToDeg()) ; 
    outputContainer->Add(fhPrimAccPhi) ;
    
    
    fhPrimOpeningAngle  = new TH2D
      ("hPrimOpeningAngle","Angle between all primary #gamma pair vs E_{#pi^{0}}",nptbins,ptmin,ptmax,200,0,0.5); 
    fhPrimOpeningAngle->SetYTitle("#theta(rad)");
    fhPrimOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
    outputContainer->Add(fhPrimOpeningAngle) ;
    
    fhPrimCosOpeningAngle  = new TH2D
      ("hPrimCosOpeningAngle","Cosinus of angle between all primary #gamma pair vs E_{#pi^{0}}",nptbins,ptmin,ptmax,200,-1,1); 
    fhPrimCosOpeningAngle->SetYTitle("cos (#theta) ");
    fhPrimCosOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
    outputContainer->Add(fhPrimCosOpeningAngle) ;
    
  }
  
  for(Int_t imod=0; imod<fNModules; imod++){
    //Module dependent invariant mass
    sprintf(key,"hReMod_%d",imod) ;
    sprintf(title,"Real m_{#gamma#gamma} distr. for Module %d",imod) ;
    //fhEtalon->Clone(key);
    //fhReMod[imod]=(TH3D*)fhEtalon->Clone(key) ;
    //fhReMod[imod]->SetName(key) ;
    //fhReMod[imod]->SetTitle(title) ;
    fhReMod[imod]  = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
    outputContainer->Add(fhReMod[imod]) ;
  }
  
  
//  //Save parameters used for analysis
//  TString parList ; //this will be list of parameters used for this analysis.
//  char onePar[255] ;
//  sprintf(onePar,"--- AliAnaPi0 ---\n") ;
//  parList+=onePar ;	
//  sprintf(onePar,"Number of bins in Centrality:  %d \n",fNCentrBin) ;
//  parList+=onePar ;
//  sprintf(onePar,"Number of bins in Z vert. pos: %d \n",fNZvertBin) ;
//  parList+=onePar ;
//  sprintf(onePar,"Number of bins in Reac. Plain: %d \n",fNrpBin) ;
//  parList+=onePar ;
//  sprintf(onePar,"Depth of event buffer: %d \n",fNmaxMixEv) ;
//  parList+=onePar ;
//  sprintf(onePar,"Number of different PID used:  %d \n",fNPID) ;
//  parList+=onePar ;
//  sprintf(onePar,"Cuts: \n") ;
//  parList+=onePar ;
//  sprintf(onePar,"Z vertex position: -%f < z < %f \n",fZvtxCut,fZvtxCut) ;
//  parList+=onePar ;
//  sprintf(onePar,"Calorimeter: %s \n",fCalorimeter.Data()) ;
//  parList+=onePar ;
//  sprintf(onePar,"Number of modules: %d \n",fNModules) ;
//  parList+=onePar ;
//  
//  TObjString *oString= new TObjString(parList) ;
//  outputContainer->Add(oString);
  
  return outputContainer;
}

//_________________________________________________________________________________________________________________________________________________
void AliAnaPi0::Print(const Option_t * /*opt*/) const
{
  //Print some relevant parameters set for the analysis
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");

  printf("Number of bins in Centrality:  %d \n",fNCentrBin) ;
  printf("Number of bins in Z vert. pos: %d \n",fNZvertBin) ;
  printf("Number of bins in Reac. Plain: %d \n",fNrpBin) ;
  printf("Depth of event buffer: %d \n",fNmaxMixEv) ;
  printf("Number of different PID used:  %d \n",fNPID) ;
  printf("Cuts: \n") ;
  printf("Z vertex position: -%2.3f < z < %2.3f \n",fZvtxCut,fZvtxCut) ;
  printf("Number of modules:             %d \n",fNModules) ;
  printf("Select pairs with their angle: %d \n",fUseAngleCut) ;
  printf("------------------------------------------------------\n") ;
} 


//____________________________________________________________________________________________________________________________________________________
void AliAnaPi0::MakeAnalysisFillHistograms() 
{
  //Process one event and extract photons from AOD branch 
  // filled with AliAnaPhoton and fill histos with invariant mass
  
  //Apply some cuts on event: vertex position and centrality range  
  Int_t iRun=(GetReader()->GetInputEvent())->GetRunNumber() ;
  if(IsBadRun(iRun)) return ;	
  
  Double_t vert[]={0,0,0} ; //vertex ;
  GetReader()->GetVertex(vert);
  if(vert[2]<-fZvtxCut || vert[2]> fZvtxCut) return ; //Event can not be used (vertex, centrality,... cuts not fulfilled)
  
  //Get Centrality and calculate centrality bin
  //Does not exist in ESD yet???????
  Int_t curCentrBin=0 ;
  
  //Get Reaction Plain position and calculate RP bin
  //does not exist in ESD yet????
  Int_t curRPBin=0 ;	
  
  Int_t curZvertBin=(Int_t)(0.5*fNZvertBin*(vert[2]+fZvtxCut)/fZvtxCut) ;
  
  fhEvents->Fill(curCentrBin+0.5,curZvertBin+0.5,curRPBin+0.5) ;
  
  Int_t nPhot = GetInputAODBranch()->GetEntriesFast() ;
  if(GetDebug() > 1) printf("AliAnaPi0::MakeAnalysisFillHistograms() - Photon entries %d\n", nPhot);
  Int_t module1 = -1;
  Int_t module2 = -1;
  for(Int_t i1=0; i1<nPhot-1; i1++){
    AliAODPWG4Particle * p1 = (AliAODPWG4Particle*) (GetInputAODBranch()->At(i1)) ;
    TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());
    //Get Module number
    module1 = GetModuleNumber(p1);
    for(Int_t i2=i1+1; i2<nPhot; i2++){
      AliAODPWG4Particle * p2 = (AliAODPWG4Particle*) (GetInputAODBranch()->At(i2)) ;
      TLorentzVector photon2(p2->Px(),p2->Py(),p2->Pz(),p2->E());
      //Get module number
      module2 = GetModuleNumber(p2);
      Double_t m  = (photon1 + photon2).M() ;
      Double_t pt = (photon1 + photon2).Pt();
      Double_t a  = TMath::Abs(p1->E()-p2->E())/(p1->E()+p2->E()) ;
      if(GetDebug() > 2)
	printf("AliAnaPi0::MakeAnalysisFillHistograms() - Current Event: pT: photon1 %2.2f, photon2 %2.2f; Pair: pT %2.2f, mass %2.3f, a %f2.3\n",
	       p1->Pt(), p2->Pt(), pt,m,a);
				
      //Check if opening angle is too large or too small compared to what is expected	
      Double_t angle   = photon1.Angle(photon2.Vect());
      //if(fUseAngleCut && !GetNeutralMesonSelection()->IsAngleInWindow((photon1+photon2).E(),angle)) continue;
      //printf("angle %f\n",angle);
      if(fUseAngleCut && angle < 0.1) continue;
      fhRealOpeningAngle   ->Fill(pt,angle);
      fhRealCosOpeningAngle->Fill(pt,TMath::Cos(angle));
      
      //Fill module dependent histograms
      //if(module1==module2) printf("mod1 %d\n",module1);
      if(module1==module2 && module1 >=0 && module1<fNModules)
	fhReMod[module1]->Fill(pt,a,m) ;
      
      for(Int_t ipid=0; ipid<fNPID; ipid++)
	{
	  if((p1->IsPIDOK(ipid,AliCaloPID::kPhoton)) && (p2->IsPIDOK(ipid,AliCaloPID::kPhoton))){ 
	    fhRe1[curCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
	    if(p1->DistToBad()>0 && p2->DistToBad()>0){
	      fhRe2[curCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
	      if(p1->DistToBad()>1 && p2->DistToBad()>1){
		fhRe3[curCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
	      }
	    }
	  }
	} 
    }
  }
  
  //Fill mixed
  
  TList * evMixList=fEventsList[curCentrBin*fNZvertBin*fNrpBin+curZvertBin*fNrpBin+curRPBin] ;
  Int_t nMixed = evMixList->GetSize() ;
  for(Int_t ii=0; ii<nMixed; ii++){  
    TClonesArray* ev2= (TClonesArray*) (evMixList->At(ii));
    Int_t nPhot2=ev2->GetEntriesFast() ;
    Double_t m = -999;
    if(GetDebug() > 1) printf("AliAnaPi0::MakeAnalysisFillHistograms() - Mixed event %d photon entries %d\n", ii, nPhot);
    
    for(Int_t i1=0; i1<nPhot; i1++){
      AliAODPWG4Particle * p1 = (AliAODPWG4Particle*) (GetInputAODBranch()->At(i1)) ;
      TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());
      for(Int_t i2=0; i2<nPhot2; i2++){
	AliAODPWG4Particle * p2 = (AliAODPWG4Particle*) (ev2->At(i2)) ;
	
	TLorentzVector photon2(p2->Px(),p2->Py(),p2->Pz(),p2->E());
	m =           (photon1+photon2).M() ; 
	Double_t pt = (photon1 + photon2).Pt();
	Double_t a  = TMath::Abs(p1->E()-p2->E())/(p1->E()+p2->E()) ;
	
	//Check if opening angle is too large or too small compared to what is expected
	Double_t angle   = photon1.Angle(photon2.Vect());
	//if(fUseAngleCut && !GetNeutralMesonSelection()->IsAngleInWindow((photon1+photon2).E(),angle)) continue;
	if(fUseAngleCut && angle < 0.1) continue;  
	
	if(GetDebug() > 2)
	  printf("AliAnaPi0::MakeAnalysisFillHistograms() - Mixed Event: pT: photon1 %2.2f, photon2 %2.2f; Pair: pT %2.2f, mass %2.3f, a %f2.3\n",
		 p1->Pt(), p2->Pt(), pt,m,a);			
	for(Int_t ipid=0; ipid<fNPID; ipid++){ 
	  if((p1->IsPIDOK(ipid,AliCaloPID::kPhoton)) && (p2->IsPIDOK(ipid,AliCaloPID::kPhoton))){ 
	    fhMi1[curCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
	    if(p1->DistToBad()>0 && p2->DistToBad()>0){
	      fhMi2[curCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
	      if(p1->DistToBad()>1 && p2->DistToBad()>1){
		fhMi3[curCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
	      }
	      
	    }
	  }
	}
      }
    }
  }
  
  TClonesArray *currentEvent = new TClonesArray(*GetInputAODBranch());
  //Add current event to buffer and Remove redandant events 
  if(currentEvent->GetEntriesFast()>0){
    evMixList->AddFirst(currentEvent) ;
    currentEvent=0 ; //Now list of particles belongs to buffer and it will be deleted with buffer
    if(evMixList->GetSize()>=fNmaxMixEv)
      {
	TClonesArray * tmp = (TClonesArray*) (evMixList->Last()) ;
	evMixList->RemoveLast() ;
	delete tmp ;
      }
  } 
  else{ //empty event
    delete currentEvent ;
    currentEvent=0 ; 
  }
  
  //Acceptance
  if(IsDataMC() && GetReader()->ReadStack()){	
    AliStack * stack = GetMCStack();
    if(stack && (IsDataMC() || (GetReader()->GetDataType() == AliCaloTrackReader::kMC)) ){
      for(Int_t i=0 ; i<stack->GetNprimary(); i++){
	TParticle * prim = stack->Particle(i) ;
	if(prim->GetPdgCode() == 111){
	  Double_t pi0Pt = prim->Pt() ;
	  //printf("pi0, pt %2.2f\n",pi0Pt);
	  if(prim->Energy() == TMath::Abs(prim->Pz()))  continue ; //Protection against floating point exception	  
	  Double_t pi0Y  = 0.5*TMath::Log((prim->Energy()-prim->Pz())/(prim->Energy()+prim->Pz())) ;
	  Double_t phi   = TMath::RadToDeg()*prim->Phi() ;
	  if(TMath::Abs(pi0Y) < 0.5){
	    fhPrimPt->Fill(pi0Pt) ;
	  }
	  fhPrimY  ->Fill(pi0Y) ;
	  fhPrimPhi->Fill(phi) ;
	  
	  //Check if both photons hit Calorimeter
	  Int_t iphot1=prim->GetFirstDaughter() ;
	  Int_t iphot2=prim->GetLastDaughter() ;
	  if(iphot1>-1 && iphot1<stack->GetNtrack() && iphot2>-1 && iphot2<stack->GetNtrack()){
	    TParticle * phot1 = stack->Particle(iphot1) ;
	    TParticle * phot2 = stack->Particle(iphot2) ;
	    if(phot1 && phot2 && phot1->GetPdgCode()==22 && phot2->GetPdgCode()==22){
	      //printf("2 photons: photon 1: pt %2.2f, phi %3.2f, eta %1.2f; photon 2: pt %2.2f, phi %3.2f, eta %1.2f\n",
	      //	phot1->Pt(), phot1->Phi()*180./3.1415, phot1->Eta(), phot2->Pt(), phot2->Phi()*180./3.1415, phot2->Eta());
	      
	      TLorentzVector lv1, lv2;
	      phot1->Momentum(lv1);
	      phot2->Momentum(lv2);
	      
	      Bool_t inacceptance = kFALSE;
	      if(fCalorimeter == "PHOS"){
		if(GetPHOSGeometry() && GetCaloUtils()->IsPHOSGeoMatrixSet()){
		  Int_t mod ;
		  Double_t x,z ;
		  if(GetPHOSGeometry()->ImpactOnEmc(phot1,mod,z,x) && GetPHOSGeometry()->ImpactOnEmc(phot2,mod,z,x)) 
		    inacceptance = kTRUE;
		  if(GetDebug() > 2) printf("In %s Real acceptance? %d\n",fCalorimeter.Data(),inacceptance);
		}
		else{
		  
		  if(GetFiducialCut()->IsInFiducialCut(lv1,fCalorimeter) && GetFiducialCut()->IsInFiducialCut(lv2,fCalorimeter)) 
		    inacceptance = kTRUE ;
		  if(GetDebug() > 2) printf("In %s fiducial cut acceptance? %d\n",fCalorimeter.Data(),inacceptance);
		}
		
	      }	   
	      else if(fCalorimeter == "EMCAL" && GetCaloUtils()->IsEMCALGeoMatrixSet()){
		if(GetEMCALGeometry()){
		  if(GetEMCALGeometry()->Impact(phot1) && GetEMCALGeometry()->Impact(phot2)) 
		    inacceptance = kTRUE;
		  if(GetDebug() > 2) printf("In %s Real acceptance? %d\n",fCalorimeter.Data(),inacceptance);
		}
		else{
		  if(GetFiducialCut()->IsInFiducialCut(lv1,fCalorimeter) && GetFiducialCut()->IsInFiducialCut(lv2,fCalorimeter)) 
		    inacceptance = kTRUE ;
		  if(GetDebug() > 2) printf("In %s fiducial cut acceptance? %d\n",fCalorimeter.Data(),inacceptance);
		}
	      }	  
	      
	      if(inacceptance){
		
		fhPrimAccPt->Fill(pi0Pt) ;
		fhPrimAccPhi->Fill(phi) ;
		fhPrimAccY->Fill(pi0Y) ;
		Double_t angle  = lv1.Angle(lv2.Vect());
		fhPrimOpeningAngle   ->Fill(pi0Pt,angle);
		fhPrimCosOpeningAngle->Fill(pi0Pt,TMath::Cos(angle));
		
	      }//Accepted
	    }// 2 photons      
	  }//Check daughters exist
	}// Primary pi0
      }//loop on primaries	
    }//stack exists and data is MC
  }//read stack
  else if(GetReader()->ReadAODMCParticles()){
    if(GetDebug() >= 0)  printf("AliAnaPi0::MakeAnalysisFillHistograms() - Acceptance calculation with MCParticles not implemented yet\n");
  }	
  
}	

//________________________________________________________________________
void AliAnaPi0::ReadHistograms(TList* outputList)
{
  // Needed when Terminate is executed in distributed environment
  // Refill analysis histograms of this class with corresponding histograms in output list. 
  
  // Histograms of this analsys are kept in the same list as other analysis, recover the position of
  // the first one and then add the next.
  Int_t index = outputList->IndexOf(outputList->FindObject(GetAddedHistogramsStringToName()+"hRe_cen0_pid0_dist1"));
  
  if(!fhRe1) fhRe1 = new TH3D*[fNCentrBin*fNPID] ;
  if(!fhRe2) fhRe2 = new TH3D*[fNCentrBin*fNPID] ;
  if(!fhRe3) fhRe3 = new TH3D*[fNCentrBin*fNPID] ;
  if(!fhMi1) fhMi1 = new TH3D*[fNCentrBin*fNPID] ;
  if(!fhMi2) fhMi2 = new TH3D*[fNCentrBin*fNPID] ;
  if(!fhMi3) fhMi3 = new TH3D*[fNCentrBin*fNPID] ;	
  if(!fhReMod) fhReMod = new TH3D*[fNModules] ;	
  
  for(Int_t ic=0; ic<fNCentrBin; ic++){
    for(Int_t ipid=0; ipid<fNPID; ipid++){
      fhRe1[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhMi1[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhRe2[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhMi2[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhRe3[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhMi3[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
    }
  }
  
  fhEvents = (TH3D *) outputList->At(index++); 
  
  //Histograms filled only if MC data is requested 	
  if(IsDataMC() || (GetReader()->GetDataType() == AliCaloTrackReader::kMC) ){
    fhPrimPt     = (TH1D*)  outputList->At(index++);
    fhPrimAccPt  = (TH1D*)  outputList->At(index++);
    fhPrimY      = (TH1D*)  outputList->At(index++);
    fhPrimAccY   = (TH1D*)  outputList->At(index++);
    fhPrimPhi    = (TH1D*)  outputList->At(index++);
    fhPrimAccPhi = (TH1D*)  outputList->At(index++);
  }
  
  for(Int_t imod=0; imod < fNModules; imod++)
    fhReMod[imod] = (TH3D*) outputList->At(index++);
  
}


//____________________________________________________________________________________________________________________________________________________
void AliAnaPi0::Terminate(TList* outputList) 
{
  //Do some calculations and plots from the final histograms.
  
  printf(" *** %s Terminate:\n", GetName()) ; 
  
  //Recover histograms from output histograms list, needed for distributed analysis.    
  ReadHistograms(outputList);
  
  if (!fhRe1) {
    printf("AliAnaPi0::Terminate() - Error: Remote output histograms not imported in AliAnaPi0 object");
    return;
  }
  
  printf("AliAnaPi0::Terminate()         Mgg Real        : %5.3f , RMS : %5.3f \n", fhRe1[0]->GetMean(),   fhRe1[0]->GetRMS() ) ;
  
  char nameIM[128];
  sprintf(nameIM,"AliAnaPi0_%s_cPt",fCalorimeter.Data());
  TCanvas  * cIM = new TCanvas(nameIM, "", 400, 10, 600, 700) ;
  cIM->Divide(2, 2);
  
  cIM->cd(1) ; 
  //gPad->SetLogy();
  TH1D * hIMAllPt = (TH1D*) fhRe1[0]->ProjectionZ(Form("IMPtAll_%s",fCalorimeter.Data()));
  hIMAllPt->SetLineColor(2);
  hIMAllPt->SetTitle("No cut on  p_{T, #gamma#gamma} ");
  hIMAllPt->Draw();

  cIM->cd(2) ; 
  TH3F * hRe1Pt5 = (TH3F*)fhRe1[0]->Clone(Form("IMPt5_%s",fCalorimeter.Data()));
  hRe1Pt5->GetXaxis()->SetRangeUser(0,5);
  TH1D * hIMPt5 = (TH1D*) hRe1Pt5->Project3D(Form("IMPt5_%s_pz",fCalorimeter.Data()));
  hIMPt5->SetLineColor(2);  
  hIMPt5->SetTitle("0 < p_{T, #gamma#gamma} < 5 GeV/c");
  hIMPt5->Draw();
  
  cIM->cd(3) ; 
  TH3F * hRe1Pt10 =  (TH3F*)fhRe1[0]->Clone(Form("IMPt10_%s",fCalorimeter.Data()));
  hRe1Pt10->GetXaxis()->SetRangeUser(5,10);
  TH1D * hIMPt10 = (TH1D*) hRe1Pt10->Project3D(Form("IMPt10_%s_pz",fCalorimeter.Data()));
  hIMPt10->SetLineColor(2);  
  hIMPt10->SetTitle("5 < p_{T, #gamma#gamma} < 10 GeV/c");
  hIMPt10->Draw();
  
  cIM->cd(4) ; 
  TH3F * hRe1Pt20 =  (TH3F*)fhRe1[0]->Clone(Form("IMPt20_%s",fCalorimeter.Data()));
  hRe1Pt20->GetXaxis()->SetRangeUser(10,20);
  TH1D * hIMPt20 = (TH1D*) hRe1Pt20->Project3D(Form("IMPt20_%s_pz",fCalorimeter.Data()));
  hIMPt20->SetLineColor(2);  
  hIMPt20->SetTitle("10 < p_{T, #gamma#gamma} < 20 GeV/c");
  hIMPt20->Draw();
   
  char nameIMF[128];
  sprintf(nameIMF,"AliAnaPi0_%s_Mgg.eps",fCalorimeter.Data());
  cIM->Print(nameIMF);

  char namePt[128];
  sprintf(namePt,"AliAnaPi0_%s_cPt",fCalorimeter.Data());
  TCanvas  * cPt = new TCanvas(namePt, "", 400, 10, 600, 700) ;
  cPt->Divide(2, 2);

  cPt->cd(1) ; 
  //gPad->SetLogy();
  TH1D * hPt = (TH1D*) fhRe1[0]->Project3D("x");
  hPt->SetLineColor(2);
  hPt->SetTitle("No cut on  M_{#gamma#gamma} ");
  hPt->Draw();

  cPt->cd(2) ; 
  TH3F * hRe1IM1 = (TH3F*)fhRe1[0]->Clone(Form("Pt1_%s",fCalorimeter.Data()));
  hRe1IM1->GetZaxis()->SetRangeUser(0.05,0.21);
  TH1D * hPtIM1 = (TH1D*) hRe1IM1->Project3D("x");
  hPtIM1->SetLineColor(2);  
  hPtIM1->SetTitle("0.05 < M_{#gamma#gamma} < 0.21 GeV/c^{2}");
  hPtIM1->Draw();
  
  cPt->cd(3) ; 
  TH3F * hRe1IM2 = (TH3F*)fhRe1[0]->Clone(Form("Pt2_%s",fCalorimeter.Data()));
  hRe1IM2->GetZaxis()->SetRangeUser(0.09,0.17);
  TH1D * hPtIM2 = (TH1D*) hRe1IM2->Project3D("x");
  hPtIM2->SetLineColor(2);  
  hPtIM2->SetTitle("0.09 < M_{#gamma#gamma} < 0.17 GeV/c^{2}");
  hPtIM2->Draw();

  cPt->cd(4) ; 
  TH3F * hRe1IM3 = (TH3F*)fhRe1[0]->Clone(Form("Pt3_%s",fCalorimeter.Data()));
  hRe1IM3->GetZaxis()->SetRangeUser(0.11,0.15);
  TH1D * hPtIM3 = (TH1D*) hRe1IM1->Project3D("x");
  hPtIM3->SetLineColor(2);  
  hPtIM3->SetTitle("0.11 < M_{#gamma#gamma} < 0.15 GeV/c^{2}");
  hPtIM3->Draw();
   
  char namePtF[128];
  sprintf(namePtF,"AliAnaPi0_%s_Pt.eps",fCalorimeter.Data());
  cPt->Print(namePtF);

 
  char line[1024] ; 
  sprintf(line, ".!tar -zcf %s_%s.tar.gz *.eps", GetName(),fCalorimeter.Data()) ; 
  gROOT->ProcessLine(line);
  sprintf(line, ".!rm -fR AliAnaPi0_%s*.eps",fCalorimeter.Data()); 
  gROOT->ProcessLine(line);
 
  printf(" AliAnaPi0::Terminate() - !! All the eps files are in %s_%s.tar.gz !!!\n", GetName(), fCalorimeter.Data());

}





