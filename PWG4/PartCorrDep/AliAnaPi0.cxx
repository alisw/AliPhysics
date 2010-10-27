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
#include "TObjString.h"

//---- AliRoot system ----
#include "AliAnaPi0.h"
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"
#include "AliStack.h"
#include "AliFiducialCut.h"
#include "TParticle.h"
#include "AliVEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliNeutralMesonSelection.h"
#include "AliMixedEvent.h"


ClassImp(AliAnaPi0)

//________________________________________________________________________________________________________________________________________________  
AliAnaPi0::AliAnaPi0() : AliAnaPartCorrBaseClass(),
fDoOwnMix(kFALSE),fNCentrBin(0),//fNZvertBin(0),fNrpBin(0),
fNPID(0),fNmaxMixEv(0), fCalorimeter(""),
fNModules(12), fUseAngleCut(kFALSE), fEventsList(0x0), fMultiCutAna(kFALSE),
fNPtCuts(0),fPtCuts(0x0),fNAsymCuts(0),fAsymCuts(0x0),
fNCellNCuts(0),fCellNCuts(0x0),fNPIDBits(0),fPIDBits(0x0),fhReMod(0x0),fhReDiffMod(0x0),
fhRe1(0x0),      fhMi1(0x0),      fhRe2(0x0),      fhMi2(0x0),      fhRe3(0x0),      fhMi3(0x0),
fhReInvPt1(0x0), fhMiInvPt1(0x0), fhReInvPt2(0x0), fhMiInvPt2(0x0), fhReInvPt3(0x0), fhMiInvPt3(0x0),
fhRePtNCellAsymCuts(0x0), fhRePIDBits(0x0),
fhRePtMult(0x0), fhRePtMultAsy07(0x0)  ,   
fhEvents(0x0), fhRealOpeningAngle(0x0),fhRealCosOpeningAngle(0x0),
fhPrimPt(0x0), fhPrimAccPt(0x0), fhPrimY(0x0), fhPrimAccY(0x0), fhPrimPhi(0x0), fhPrimAccPhi(0x0),
fhPrimOpeningAngle(0x0),fhPrimCosOpeningAngle(0x0)
{
//Default Ctor
 InitParameters();
 
}

//________________________________________________________________________________________________________________________________________________
AliAnaPi0::~AliAnaPi0() {
  // Remove event containers
  
  if(fDoOwnMix && fEventsList){
    for(Int_t ic=0; ic<fNCentrBin; ic++){
      for(Int_t iz=0; iz<GetNZvertBin(); iz++){
        for(Int_t irp=0; irp<GetNRPBin(); irp++){
          fEventsList[ic*GetNZvertBin()*GetNRPBin()+iz*GetNRPBin()+irp]->Delete() ;
          delete fEventsList[ic*GetNZvertBin()*GetNRPBin()+iz*GetNRPBin()+irp] ;
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
//  fNZvertBin = 1;
//  fNrpBin    = 1;
  fNPID      = 9;
  fNmaxMixEv = 10;
 
  fCalorimeter  = "PHOS";
  fUseAngleCut = kFALSE;
  
  fMultiCutAna = kFALSE;
  
  fNPtCuts = 3;
  fPtCuts  = new Float_t[fNPtCuts];
  fPtCuts[0] = 0.; fPtCuts[1] = 0.2;   fPtCuts[2] = 0.3;

  fNAsymCuts = 3;
  fAsymCuts  = new Float_t[fNAsymCuts];
  fAsymCuts[0] = 0.7;  fAsymCuts[1] = 0.8;   fAsymCuts[2] = 1.;   
  
  fNCellNCuts = 3;
  fCellNCuts  = new Int_t[fNCellNCuts];
  fCellNCuts[0] = 1; fCellNCuts[1] = 2;   fCellNCuts[2] = 3;   
  
  fNPIDBits = 3;
  fPIDBits  = new Int_t[fNPIDBits];
  fPIDBits[0] = 2;   fPIDBits[1] = 4;   fPIDBits[2] = 6; // check dispersion, neutral, dispersion&&neutral
  
}


//________________________________________________________________________________________________________________________________________________
TObjString * AliAnaPi0::GetAnalysisCuts()
{  
	 //Save parameters used for analysis
	 TString parList ; //this will be list of parameters used for this analysis.
   const Int_t buffersize = 255;
	 char onePar[buffersize] ;
	 snprintf(onePar,buffersize,"--- AliAnaPi0 ---\n") ;
	 parList+=onePar ;	
	 snprintf(onePar,buffersize,"Number of bins in Centrality:  %d \n",fNCentrBin) ;
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"Number of bins in Z vert. pos: %d \n",GetNZvertBin()) ;
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"Number of bins in Reac. Plain: %d \n",GetNRPBin()) ;
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"Depth of event buffer: %d \n",fNmaxMixEv) ;
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"Number of different PID used:  %d \n",fNPID) ;
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"Cuts: \n") ;
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"Z vertex position: -%f < z < %f \n",GetZvertexCut(),GetZvertexCut()) ;
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"Calorimeter: %s \n",fCalorimeter.Data()) ;
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"Number of modules: %d \n",fNModules) ;
	 parList+=onePar ;
  if(fMultiCutAna){
    snprintf(onePar, buffersize," pT cuts: n = %d, pt > ",fNPtCuts) ;
    for(Int_t i = 0; i < fNPtCuts; i++) snprintf(onePar,buffersize,"%s %2.2f;",onePar,fPtCuts[i]);
    parList+=onePar ;
    
    snprintf(onePar,buffersize, " N cell in cluster cuts: n = %d, nCell > ",fNCellNCuts) ;
    for(Int_t i = 0; i < fNCellNCuts; i++) snprintf(onePar,buffersize,"%s %d;",onePar,fCellNCuts[i]);
    parList+=onePar ;
    
    snprintf(onePar,buffersize," Asymmetry cuts: n = %d, asymmetry < ",fNAsymCuts) ;
    for(Int_t i = 0; i < fNAsymCuts; i++) snprintf(onePar,buffersize,"%s %2.2f;",onePar,fAsymCuts[i]);
    parList+=onePar ;
    
    snprintf(onePar,buffersize," PID selection bits: n = %d, PID bit =\n",fNPIDBits) ;
    for(Int_t i = 0; i < fNPIDBits; i++) snprintf(onePar,buffersize,"%s %d;",onePar,fPIDBits[i]);
    parList+=onePar ;
  }
  
	 return new TObjString(parList) ;	
}

//________________________________________________________________________________________________________________________________________________
TList * AliAnaPi0::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer
  
  //create event containers
  fEventsList = new TList*[fNCentrBin*GetNZvertBin()*GetNRPBin()] ;
	
  for(Int_t ic=0; ic<fNCentrBin; ic++){
    for(Int_t iz=0; iz<GetNZvertBin(); iz++){
      for(Int_t irp=0; irp<GetNRPBin(); irp++){
        fEventsList[ic*GetNZvertBin()*GetNRPBin()+iz*GetNRPBin()+irp] = new TList() ;
      }
    }
  }
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName(GetName()); 
	
  fhReMod = new TH3D*[fNModules] ;
  fhReDiffMod = new TH3D*[fNModules+1] ;
  
  fhRe1 = new TH3D*[fNCentrBin*fNPID] ;
  fhRe2 = new TH3D*[fNCentrBin*fNPID] ;
  fhRe3 = new TH3D*[fNCentrBin*fNPID] ;
  fhMi1 = new TH3D*[fNCentrBin*fNPID] ;
  fhMi2 = new TH3D*[fNCentrBin*fNPID] ;
  fhMi3 = new TH3D*[fNCentrBin*fNPID] ;
    
  fhReInvPt1 = new TH3D*[fNCentrBin*fNPID] ;
  fhReInvPt2 = new TH3D*[fNCentrBin*fNPID] ;
  fhReInvPt3 = new TH3D*[fNCentrBin*fNPID] ;
  fhMiInvPt1 = new TH3D*[fNCentrBin*fNPID] ;
  fhMiInvPt2 = new TH3D*[fNCentrBin*fNPID] ;
  fhMiInvPt3 = new TH3D*[fNCentrBin*fNPID] ;
    
  const Int_t buffersize = 255;
  char key[buffersize] ;
  char title[buffersize] ;
  
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
      snprintf(key, buffersize,"hRe_cen%d_pid%d_dist1",ic,ipid) ;
      snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      fhRe1[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhRe1[ic*fNPID+ipid]) ;
            
      //Distance to bad module 2
      snprintf(key, buffersize,"hRe_cen%d_pid%d_dist2",ic,ipid) ;
      snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      fhRe2[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhRe2[ic*fNPID+ipid]) ;
      
      //Distance to bad module 3
      snprintf(key, buffersize,"hRe_cen%d_pid%d_dist3",ic,ipid) ;
      snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      fhRe3[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhRe3[ic*fNPID+ipid]) ;
      
      //Inverse pT 
      //Distance to bad module 1
      snprintf(key, buffersize,"hReInvPt_cen%d_pid%d_dist1",ic,ipid) ;
      snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      fhReInvPt1[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhReInvPt1[ic*fNPID+ipid]) ;
      
      //Distance to bad module 2
      snprintf(key, buffersize,"hReInvPt_cen%d_pid%d_dist2",ic,ipid) ;
      snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      fhReInvPt2[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhReInvPt2[ic*fNPID+ipid]) ;
      
      //Distance to bad module 3
      snprintf(key, buffersize,"hReInvPt_cen%d_pid%d_dist3",ic,ipid) ;
      snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      fhReInvPt3[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhReInvPt3[ic*fNPID+ipid]) ;
      
      if(fDoOwnMix){
        //Distance to bad module 1
        snprintf(key, buffersize,"hMi_cen%d_pid%d_dist1",ic,ipid) ;
        snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
        fhMi1[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
        outputContainer->Add(fhMi1[ic*fNPID+ipid]) ;
        
        //Distance to bad module 2
        snprintf(key, buffersize,"hMi_cen%d_pid%d_dist2",ic,ipid) ;
        snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
        fhMi2[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
        outputContainer->Add(fhMi2[ic*fNPID+ipid]) ;
        
        //Distance to bad module 3
        snprintf(key, buffersize,"hMi_cen%d_pid%d_dist3",ic,ipid) ;
        snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
        fhMi3[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
        outputContainer->Add(fhMi3[ic*fNPID+ipid]) ;
        
        //Inverse pT
        //Distance to bad module 1
        snprintf(key, buffersize,"hMiInvPt_cen%d_pid%d_dist1",ic,ipid) ;
        snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
        fhMiInvPt1[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
        outputContainer->Add(fhMiInvPt1[ic*fNPID+ipid]) ;
        
        //Distance to bad module 2
        snprintf(key, buffersize,"hMiInvPt_cen%d_pid%d_dist2",ic,ipid) ;
        snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
        fhMiInvPt2[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
        outputContainer->Add(fhMiInvPt2[ic*fNPID+ipid]) ;
        
        //Distance to bad module 3
        snprintf(key, buffersize,"hMiInvPt_cen%d_pid%d_dist3",ic,ipid) ;
        snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
        fhMiInvPt3[ic*fNPID+ipid] = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
        outputContainer->Add(fhMiInvPt3[ic*fNPID+ipid]) ;
        
        
      }
    }
  }
  
  if(fMultiCutAna){
    
    fhRePIDBits         = new TH2D*[fNPIDBits];
    for(Int_t ipid=0; ipid<fNPIDBits; ipid++){
      snprintf(key,   buffersize,"hRe_pidbit%d",ipid) ;
      snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for PIDBit=%d",fPIDBits[ipid]) ;
      fhRePIDBits[ipid] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
      outputContainer->Add(fhRePIDBits[ipid]) ;
    }// pid bit loop
    
    fhRePtNCellAsymCuts = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
    for(Int_t ipt=0; ipt<fNPtCuts; ipt++){
      for(Int_t icell=0; icell<fNCellNCuts; icell++){
        for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
          snprintf(key,   buffersize,"hRe_pt%d_cell%d_asym%d",ipt,icell,iasym) ;
          snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for pt >%2.2f, ncell>%d and asym >%2.2f ",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]) ;
          Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
          //printf("ipt %d, icell %d, iassym %d, index %d\n",ipt, icell, iasym, index);
          fhRePtNCellAsymCuts[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          outputContainer->Add(fhRePtNCellAsymCuts[index]) ;
        }
      }
    }
    
    fhRePtMult      = new TH3D("hRePtMult","(p_{T},C,M)_{#gamma#gamma}, 0<A<1.0",      nptbins,ptmin,ptmax,40,0.,160.,nmassbins,massmin,massmax);
    fhRePtMultAsy07 = new TH3D("hRePtMultAsy07","(p_{T},C,M)_{#gamma#gamma}, 0<A<0.7",nptbins,ptmin,ptmax,40,0.,160.,nmassbins,massmin,massmax);
    outputContainer->Add(fhRePtMult) ;
    outputContainer->Add(fhRePtMultAsy07) ;

  }// multi cuts analysis
  
  fhEvents=new TH3D("hEvents","Number of events",fNCentrBin,0.,1.*fNCentrBin,
                    GetNZvertBin(),0.,1.*GetNZvertBin(),GetNRPBin(),0.,1.*GetNRPBin()) ;
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
  if(IsDataMC()){

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
    ("hPrimOpeningAngle","Angle between all primary #gamma pair vs E_{#pi^{0}}",nptbins,ptmin,ptmax,100,0,0.5); 
    fhPrimOpeningAngle->SetYTitle("#theta(rad)");
    fhPrimOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
    outputContainer->Add(fhPrimOpeningAngle) ;
    
    fhPrimCosOpeningAngle  = new TH2D
    ("hPrimCosOpeningAngle","Cosinus of angle between all primary #gamma pair vs E_{#pi^{0}}",nptbins,ptmin,ptmax,100,-1,1); 
    fhPrimCosOpeningAngle->SetYTitle("cos (#theta) ");
    fhPrimCosOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
    outputContainer->Add(fhPrimCosOpeningAngle) ;
    
  }
  
  TString * pairname = new TString[fNModules];
  if(fCalorimeter=="EMCAL"){
    pairname[0]="A side (0-2)"; 
    pairname[1]="C side (1-3)";
    pairname[2]="Sector 0 (0-1)"; 
    pairname[3]="Sector 1 (2-3)";
    for(Int_t i = 4 ; i < fNModules ; i++) pairname[i]="";}
  if(fCalorimeter=="PHOS") {
    pairname[0]="(0-1)"; 
    pairname[1]="(0-2)";
    pairname[2]="(1-3)";
    for(Int_t i = 3 ; i < fNModules ; i++) pairname[i]="";}

  for(Int_t imod=0; imod<fNModules; imod++){
    //Module dependent invariant mass
    snprintf(key, buffersize,"hReMod_%d",imod) ;
    snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for Module %d",imod) ;
    fhReMod[imod]  = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
    outputContainer->Add(fhReMod[imod]) ;

    snprintf(key, buffersize,"hReDiffMod_%d",imod) ;
    snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for Different Modules: %s",(pairname[imod]).Data()) ;
    fhReDiffMod[imod]  = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
    outputContainer->Add(fhReDiffMod[imod]) ;
  }
  
  delete [] pairname;
  
  snprintf(key, buffersize,"hReDiffMod_%d",fNModules) ;
  snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for all Modules Combination") ;
  fhReDiffMod[fNModules]  = new TH3D(key,title,nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax,nmassbins,massmin,massmax) ;
  outputContainer->Add(fhReDiffMod[fNModules]) ;
  
  
//  for(Int_t i = 0; i < outputContainer->GetEntries() ; i++){
//  
//    printf("Histogram %d, name: %s\n ",i, outputContainer->At(i)->GetName());
//  
//  }
  
  return outputContainer;
}

//_________________________________________________________________________________________________________________________________________________
void AliAnaPi0::Print(const Option_t * /*opt*/) const
{
  //Print some relevant parameters set for the analysis
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");

  printf("Number of bins in Centrality:  %d \n",fNCentrBin) ;
  printf("Number of bins in Z vert. pos: %d \n",GetNZvertBin()) ;
  printf("Number of bins in Reac. Plain: %d \n",GetNRPBin()) ;
  printf("Depth of event buffer: %d \n",fNmaxMixEv) ;
  printf("Number of different PID used:  %d \n",fNPID) ;
  printf("Cuts: \n") ;
  printf("Z vertex position: -%2.3f < z < %2.3f \n",GetZvertexCut(),GetZvertexCut()) ;
  printf("Number of modules:             %d \n",fNModules) ;
  printf("Select pairs with their angle: %d \n",fUseAngleCut) ;
  if(fMultiCutAna){
    printf("pT cuts: n = %d, \n",fNPtCuts) ;
    printf("\tpT > ");
    for(Int_t i = 0; i < fNPtCuts; i++) printf("%2.2f ",fPtCuts[i]);
    printf("GeV/c\n");
    
    printf("N cell in cluster cuts: n = %d, \n",fNCellNCuts) ;
    printf("\tnCell > ");
    for(Int_t i = 0; i < fNCellNCuts; i++) printf("%d ",fCellNCuts[i]);
    printf("\n");

    printf("Asymmetry cuts: n = %d, \n",fNAsymCuts) ;
    printf("\tasymmetry < ");
    for(Int_t i = 0; i < fNAsymCuts; i++) printf("%2.2f ",fAsymCuts[i]);
    printf("\n");

    printf("PID selection bits: n = %d, \n",fNPIDBits) ;
    printf("\tPID bit = ");
    for(Int_t i = 0; i < fNPIDBits; i++) printf("%d ",fPIDBits[i]);
    printf("\n");
  
  }
  printf("------------------------------------------------------\n") ;
} 

//_____________________________________________________________
void AliAnaPi0::FillAcceptanceHistograms(){
  //Fill acceptance histograms if MC data is available
  
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
    if(GetDebug() >= 0)  printf("AliAnaPi0::FillAcceptanceHistograms() - Acceptance calculation with MCParticles not implemented yet\n");
  }	
}

//____________________________________________________________________________________________________________________________________________________
void AliAnaPi0::MakeAnalysisFillHistograms() 
{
  //Process one event and extract photons from AOD branch 
  // filled with AliAnaPhoton and fill histos with invariant mass
  
  //In case of MC data, fill acceptance histograms
  FillAcceptanceHistograms();
  
  //Apply some cuts on event: vertex position and centrality range  
  Int_t iRun=(GetReader()->GetInputEvent())->GetRunNumber() ;
  if(IsBadRun(iRun)) return ;	
  
  Int_t nPhot = GetInputAODBranch()->GetEntriesFast() ;
  if(GetDebug() > 1) 
    printf("AliAnaPi0::MakeAnalysisFillHistograms() - Photon entries %d\n", nPhot);
  if(nPhot < 2 )
    return ; 
  Int_t module1 = -1;
  Int_t module2 = -1;
  Double_t vert[] = {0.0, 0.0, 0.0} ; //vertex 
  Int_t evtIndex1 = 0 ; 
  Int_t currentEvtIndex = -1 ; 
  Int_t curCentrBin = 0 ; 
  Int_t curRPBin    = 0 ; 
  Int_t curZvertBin = 0 ;
  
  for(Int_t i1=0; i1<nPhot-1; i1++){
    AliAODPWG4Particle * p1 = (AliAODPWG4Particle*) (GetInputAODBranch()->At(i1)) ;
    // get the event index in the mixed buffer where the photon comes from 
    // in case of mixing with analysis frame, not own mixing
    evtIndex1 = GetEventIndex(p1, vert) ; 
    //printf("charge = %d\n", track->Charge());
    if ( evtIndex1 == -1 )
      return ; 
    if ( evtIndex1 == -2 )
      continue ; 
    if (evtIndex1 != currentEvtIndex) {
      //Get Reaction Plan position and calculate RP bin
      //does not exist in ESD yet????
      curCentrBin = 0 ; 
      curRPBin    = 0 ;
      curZvertBin = (Int_t)(0.5*GetNZvertBin()*(vert[2]+GetZvertexCut())/GetZvertexCut()) ;
      fhEvents->Fill(curCentrBin+0.5,curZvertBin+0.5,curRPBin+0.5) ;
      currentEvtIndex = evtIndex1 ; 
    }
    
    //printf("AliAnaPi0::MakeAnalysisFillHistograms(): Photon 1 Evt %d  Vertex : %f,%f,%f\n",evtIndex1, GetVertex(evtIndex1)[0] ,GetVertex(evtIndex1)[1],GetVertex(evtIndex1)[2]);

    TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());
    //Get Module number
    module1 = GetModuleNumber(p1);
    for(Int_t i2=i1+1; i2<nPhot; i2++){
      AliAODPWG4Particle * p2 = (AliAODPWG4Particle*) (GetInputAODBranch()->At(i2)) ;
      Int_t evtIndex2 = GetEventIndex(p2, vert) ; 
      if ( evtIndex2 == -1 )
        return ; 
      if ( evtIndex2 == -2 )
        continue ;    
      if (GetMixedEvent() && (evtIndex1 == evtIndex2))
        continue ;
      //printf("AliAnaPi0::MakeAnalysisFillHistograms(): Photon 2 Evt %d  Vertex : %f,%f,%f\n",evtIndex2, GetVertex(evtIndex2)[0] ,GetVertex(evtIndex2)[1],GetVertex(evtIndex2)[2]);
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
      if(fUseAngleCut && angle < 0.1) 
        continue;
      fhRealOpeningAngle   ->Fill(pt,angle);
      fhRealCosOpeningAngle->Fill(pt,TMath::Cos(angle));
      
      //Fill module dependent histograms
      //if(module1==module2) printf("mod1 %d\n",module1);
      if(module1==module2 && module1 >=0 && module1<fNModules)
        fhReMod[module1]->Fill(pt,a,m) ;
      else  
        fhReDiffMod[fNModules]->Fill(pt,a,m) ;
      
      if(fCalorimeter=="EMCAL"){
        if((module1==0 && module2==2) || (module1==2 && module2==0)) fhReDiffMod[0]->Fill(pt,a,m) ; 
        if((module1==1 && module2==3) || (module1==3 && module2==1)) fhReDiffMod[1]->Fill(pt,a,m) ; 
        if((module1==0 && module2==1) || (module1==1 && module2==0)) fhReDiffMod[2]->Fill(pt,a,m) ;
        if((module1==2 && module2==3) || (module1==3 && module2==2)) fhReDiffMod[3]->Fill(pt,a,m) ; 
      }
      else {
        if((module1==0 && module2==1) || (module1==1 && module2==0)) fhReDiffMod[0]->Fill(pt,a,m) ; 
        if((module1==0 && module2==2) || (module1==2 && module2==0)) fhReDiffMod[1]->Fill(pt,a,m) ; 
        if((module1==1 && module2==2) || (module1==2 && module2==1)) fhReDiffMod[2]->Fill(pt,a,m) ;
      }
      
      for(Int_t ipid=0; ipid<fNPID; ipid++){
        if((p1->IsPIDOK(ipid,AliCaloPID::kPhoton)) && (p2->IsPIDOK(ipid,AliCaloPID::kPhoton))){ 
          fhRe1     [curCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
          fhReInvPt1[curCentrBin*fNPID+ipid]->Fill(pt,a,m,1./pt) ;
          if(p1->DistToBad()>0 && p2->DistToBad()>0){
            fhRe2     [curCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
            fhReInvPt2[curCentrBin*fNPID+ipid]->Fill(pt,a,m,1./pt) ;
            if(p1->DistToBad()>1 && p2->DistToBad()>1){
              fhRe3     [curCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
              fhReInvPt3[curCentrBin*fNPID+ipid]->Fill(pt,a,m,1./pt) ;
            }// bad 3
          }// bad2
        }// bad 1
      }// pid loop
      
      //Multi cuts analysis 
      if(fMultiCutAna){
        //Histograms for different PID bits selection
        for(Int_t ipid=0; ipid<fNPIDBits; ipid++){
          
          if(p1->IsPIDOK(fPIDBits[ipid],AliCaloPID::kPhoton)    && 
             p2->IsPIDOK(fPIDBits[ipid],AliCaloPID::kPhoton))   fhRePIDBits[ipid]->Fill(pt,m) ;
          
          //printf("ipt %d, ipid%d, name %s\n",ipt, ipid, fhRePtPIDCuts[ipt*fNPIDBits+ipid]->GetName());
        } // pid bit cut loop
        
        //Several pt,ncell and asymetry cuts
        //Get the number of cells
        Int_t ncell1 = 0;
        Int_t ncell2 = 0;
        AliVEvent * event = GetReader()->GetInputEvent();
        if(event){
          for(Int_t iclus = 0; iclus < event->GetNumberOfCaloClusters(); iclus++){
            AliVCluster *cluster = event->GetCaloCluster(iclus);
            
            Bool_t is = kFALSE;
            if     (fCalorimeter == "EMCAL" && GetReader()->IsEMCALCluster(cluster)) is = kTRUE;
            else if(fCalorimeter == "PHOS"  && GetReader()->IsPHOSCluster (cluster)) is = kTRUE;
            
            if(is){
              if      (p1->GetCaloLabel(0) == cluster->GetID()) ncell1 = cluster->GetNCells();
              else if (p2->GetCaloLabel(0) == cluster->GetID()) ncell2 = cluster->GetNCells();
            } // PHOS or EMCAL cluster as requested in analysis
            
            if(ncell2 > 0 &&  ncell1 > 0) break; // No need to continue the iteration
            
          }
          //printf("e 1: %2.2f, e 2: %2.2f, ncells: n1 %d, n2 %d\n", p1->E(), p2->E(),ncell1,ncell2);
        }
        for(Int_t ipt=0; ipt<fNPtCuts; ipt++){          
          for(Int_t icell=0; icell<fNCellNCuts; icell++){
            for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
              
              if(p1->Pt() >   fPtCuts[ipt]      && p2->Pt() > fPtCuts[ipt]        && 
                 a        <   fAsymCuts[iasym]                                    && 
                 ncell1   >=  fCellNCuts[icell] && ncell2   >= fCellNCuts[icell]) fhRePtNCellAsymCuts[((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym]->Fill(pt,m) ;
              
              //printf("ipt %d, icell%d, iasym %d, name %s\n",ipt, icell, iasym,  fhRePtNCellAsymCuts[((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym]->GetName());
            }// pid bit cut loop
          }// icell loop
        }// pt cut loop
        
        fhRePtMult->Fill(pt,GetTrackMultiplicity(),m) ;
        if(a < 0.7)
          fhRePtMultAsy07->Fill(pt,GetTrackMultiplicity(),m) ;

      }// multiple cuts analysis
    }// second same event particle
  }// first cluster
  
  if(fDoOwnMix){
    //Fill mixed
    TList * evMixList=fEventsList[curCentrBin*GetNZvertBin()*GetNRPBin()+curZvertBin*GetNRPBin()+curRPBin] ;
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
              fhMi1     [curCentrBin*fNPID+ipid]->Fill(pt,   a,m) ;
              fhMiInvPt1[curCentrBin*fNPID+ipid]->Fill(1./pt,a,m) ;
              if(p1->DistToBad()>0 && p2->DistToBad()>0){
                fhMi2     [curCentrBin*fNPID+ipid]->Fill(pt,   a,m) ;
                fhMiInvPt2[curCentrBin*fNPID+ipid]->Fill(1./pt,a,m) ;
                if(p1->DistToBad()>1 && p2->DistToBad()>1){
                  fhMi3     [curCentrBin*fNPID+ipid]->Fill(pt,   a,m) ;
                  fhMiInvPt3[curCentrBin*fNPID+ipid]->Fill(1./pt,a,m) ;
                }
                
              }
            }
          }//loop for histograms
        }// second cluster loop
      }//first cluster loop
    }//loop on mixed events
    
    TClonesArray *currentEvent = new TClonesArray(*GetInputAODBranch());
    //Ad  d current event to buffer and Remove redundant events 
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
  }// DoOwnMix
  
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
  if(!fhReInvPt1) fhReInvPt1  = new TH3D*[fNCentrBin*fNPID] ;
  if(!fhReInvPt2) fhReInvPt2  = new TH3D*[fNCentrBin*fNPID] ;
  if(!fhReInvPt3) fhReInvPt3  = new TH3D*[fNCentrBin*fNPID] ;
  if(!fhMiInvPt1) fhMiInvPt1  = new TH3D*[fNCentrBin*fNPID] ;
  if(!fhMiInvPt2) fhMiInvPt2  = new TH3D*[fNCentrBin*fNPID] ;
  if(!fhMiInvPt3) fhMiInvPt3  = new TH3D*[fNCentrBin*fNPID] ;	
  if(!fhReMod)    fhReMod     = new TH3D*[fNModules] ;	
  if(!fhReDiffMod)fhReDiffMod = new TH3D*[fNModules+1] ;	
  

  for(Int_t ic=0; ic<fNCentrBin; ic++){
    for(Int_t ipid=0; ipid<fNPID; ipid++){
      fhRe1[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhRe2[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhRe3[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      
      fhReInvPt1[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhReInvPt2[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhReInvPt3[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      
      fhMi1[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhMi2[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhMi3[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      
      fhMiInvPt1[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhMiInvPt2[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);
      fhMiInvPt3[ic*fNPID+ipid] = (TH3D*) outputList->At(index++);      
    }
  }
  
  if(fMultiCutAna){
    
    if(!fhRePtNCellAsymCuts) fhRePtNCellAsymCuts = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
    if(!fhRePIDBits)         fhRePIDBits         = new TH2D*[fNPIDBits];

    for(Int_t ipid=0; ipid<fNPIDBits; ipid++){
      fhRePIDBits[ipid] = (TH2D*) outputList->At(index++);
    }// ipid loop
    
    for(Int_t ipt=0; ipt<fNPtCuts; ipt++){
      for(Int_t icell=0; icell<fNCellNCuts; icell++){
        for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
          fhRePtNCellAsymCuts[((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym] = (TH2D*) outputList->At(index++);
        }// iasym
      }// icell loop
    }// ipt loop
    fhRePtMult      = (TH3D*) outputList->At(index++);
    fhRePtMultAsy07 = (TH3D*) outputList->At(index++);
  }// multi cut analysis 
  
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
    
  const Int_t buffersize = 255;

  char nameIM[buffersize];
  snprintf(nameIM, buffersize,"AliAnaPi0_%s_cPt",fCalorimeter.Data());
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
   
  char nameIMF[buffersize];
  snprintf(nameIMF,buffersize,"AliAnaPi0_%s_Mgg.eps",fCalorimeter.Data());
  cIM->Print(nameIMF);

  char namePt[buffersize];
  snprintf(namePt,buffersize,"AliAnaPi0_%s_cPt",fCalorimeter.Data());
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
   
  char namePtF[buffersize];
  snprintf(namePtF,buffersize,"AliAnaPi0_%s_Pt.eps",fCalorimeter.Data());
  cPt->Print(namePtF);

  char line[buffersize] ; 
  snprintf(line,buffersize,".!tar -zcf %s_%s.tar.gz *.eps", GetName(),fCalorimeter.Data()) ; 
  gROOT->ProcessLine(line);
  snprintf(line, buffersize,".!rm -fR AliAnaPi0_%s*.eps",fCalorimeter.Data()); 
  gROOT->ProcessLine(line);
 
  printf(" AliAnaPi0::Terminate() - !! All the eps files are in %s_%s.tar.gz !!!\n", GetName(), fCalorimeter.Data());

}
  //____________________________________________________________________________________________________________________________________________________
Int_t AliAnaPi0::GetEventIndex(AliAODPWG4Particle * part, Double_t * vert)  
{
  // retieves the event index and checks the vertex
  //    in the mixed buffer returns -2 if vertex NOK
  //    for normal events   returns 0 if vertex OK and -1 if vertex NOK
  
  Int_t evtIndex = -1 ; 
  if(GetReader()->GetDataType()!=AliCaloTrackReader::kMC){
    
    if (GetMixedEvent()){
      
      evtIndex = GetMixedEvent()->EventIndexForCaloCluster(part->GetCaloLabel(0)) ;
      GetVertex(vert,evtIndex); 
      
      if(TMath::Abs(vert[2])> GetZvertexCut())
        evtIndex = -2 ; //Event can not be used (vertex, centrality,... cuts not fulfilled)
    } else {// Single event
      
      GetVertex(vert);
      
      if(TMath::Abs(vert[2])> GetZvertexCut())
        evtIndex = -1 ; //Event can not be used (vertex, centrality,... cuts not fulfilled)
      else 
        evtIndex = 0 ;
    }
  }//No MC reader
  else {
    evtIndex = 0;
    vert[0] = 0. ; 
    vert[1] = 0. ; 
    vert[2] = 0. ; 
  }
  
  return evtIndex ; 
}

