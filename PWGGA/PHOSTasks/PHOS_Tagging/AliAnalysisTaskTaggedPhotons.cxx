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

/* $Id: */

//_________________________________________________________________________
// Analysis for Tagged Photons
// Prepares all necessary histograms for calculation of 
// the yield of pi0 decay photon in calorimeter:
// Marks photons which makes pi0 with some other and
// fill invariant mass distributions for estimate background below pi0
// peak so that estimate proportion of fake pairs. 
// Fills as well controll MC histograms with clasification of the photon origin 
// and check of the proportion of truly tagged photons.
// 
//
//*-- Dmitry Blau, Dmitri Peresunko 
//////////////////////////////////////////////////////////////////////////////

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THashList.h>
#include <TFile.h>
#include <TROOT.h>

#include "AliAnalysisTaskTaggedPhotons.h" 
#include "AliAnalysisManager.h"
#include "AliAODEvent.h" 
#include "AliAODEvent.h" 
#include "AliVCluster.h" 
#include "AliAODPWG4Particle.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "TGeoManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliPHOSGeometry.h"
#include "AliTriggerAnalysis.h"
#include "AliEMCALGeometry.h"


ClassImp(AliAnalysisTaskTaggedPhotons)

//______________________________________________________________________________
AliAnalysisTaskTaggedPhotons::AliAnalysisTaskTaggedPhotons() : 
  AliAnalysisTaskSE(),
  fPHOSgeom(0x0),
  fOutputContainer(0x0),
  fStack(0x0),
  fTrackEvent(0x0),
  fPHOSEvent(0x0),
  fCurrentMixedList(0x0),
  fTriggerAnalysis(0x0),
  fZmax(0.),
  fZmin(0.),
  fPhimax(0.),
  fPhimin(0.),
  fCentrality(0.),
  fCentBin(0) 
{
  //Deafult constructor
  //no memory allocations
  for(Int_t i=0;i<10;i++)
    for(Int_t j=0;j<2;j++)
      fPHOSEvents[i][j]=0x0 ;    //Container for PHOS photons
  for(Int_t mod=0; mod<6; mod++)
    fPHOSBadMap[mod]=0x0 ;
}
//______________________________________________________________________________
AliAnalysisTaskTaggedPhotons::AliAnalysisTaskTaggedPhotons(const char *name) : 
  AliAnalysisTaskSE(name),
  fPHOSgeom(0x0),
  fOutputContainer(0x0),
  fStack(0x0),
  fTrackEvent(0x0),
  fPHOSEvent(0x0),
  fCurrentMixedList(0x0),
  fTriggerAnalysis(new AliTriggerAnalysis),
  fZmax(-60.),
  fZmin(60.),
  fPhimax(250.),
  fPhimin(320.),
  fCentrality(0.),
  fCentBin(0)
{
  // Constructor.

  // Output slots #0 write into a TH1 container
  DefineOutput(1,THashList::Class());
  // Set bad channel map (empty so far)
  for(Int_t i=0;i<1;i++)
    for(Int_t j=0;j<5;j++)
      fPHOSEvents[i][j]=0x0 ;    //Container for PHOS photons
  // Create AOD track cut

  for(Int_t mod=0; mod<6; mod++)
    fPHOSBadMap[mod]=0x0 ;
  
  
}

//____________________________________________________________________________
AliAnalysisTaskTaggedPhotons::AliAnalysisTaskTaggedPhotons(const AliAnalysisTaskTaggedPhotons& ap) :
  AliAnalysisTaskSE(ap.GetName()),  
  fPHOSgeom(0x0),
  fOutputContainer(0x0),
  fStack(0x0),
  fTrackEvent(0x0),
  fPHOSEvent(0x0),
  fCurrentMixedList(0x0),
  fTriggerAnalysis(new AliTriggerAnalysis),  
  fZmax(-60.),
  fZmin(60.),
  fPhimax(250.),
  fPhimin(320.),
  fCentrality(0.),
  fCentBin(0)
{
  // cpy ctor
  fZmax=ap.fZmax ;
  fZmin=ap.fZmin ;
  fPhimax=ap.fPhimax ;
  fPhimin=ap.fPhimin ;
  for(Int_t i=0;i<1;i++)
    for(Int_t j=0;j<5;j++)
      fPHOSEvents[i][j]=0x0 ;    //Container for PHOS photons
  for(Int_t mod=0; mod<6; mod++)
    fPHOSBadMap[mod]=0x0 ;

}

//_____________________________________________________________________________
AliAnalysisTaskTaggedPhotons& AliAnalysisTaskTaggedPhotons::operator = (const AliAnalysisTaskTaggedPhotons& ap)
{
// assignment operator

  this->~AliAnalysisTaskTaggedPhotons();
  new(this) AliAnalysisTaskTaggedPhotons(ap);
  return *this;
}

//______________________________________________________________________________
AliAnalysisTaskTaggedPhotons::~AliAnalysisTaskTaggedPhotons()
{
  // dtor
  if(fOutputContainer && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    fOutputContainer->Clear() ; 
    delete fOutputContainer ;
  }
  for(Int_t i=0;i<10;i++)
    for(Int_t j=0;j<2;j++)
      if(fPHOSEvents[i][j]){
        delete fPHOSEvents[i][j] ;
	fPHOSEvents[i][j]=0x0 ;   
      }
}
//________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::UserCreateOutputObjects()
{ 


  // Create the outputs containers
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);
  fOutputContainer->SetName(GetName()) ; 

  //QA histograms
  fOutputContainer->Add(new TH1F("hSelEvents","Event selection", 10,0.,10.)) ;
  
  //vertex distribution
  fOutputContainer->Add(new TH1F("hNvertexTracks","N of primary tracks from the primary vertex",150,0.,150.));
  fOutputContainer->Add(new TH1F("hZvertex","Z vertex",200,-50.,+50.));
  fOutputContainer->Add(new TH2F("hTrackMult","Charged track multiplicity",100,0.,100.,250,0.,500.));
  
  //centrality
  fOutputContainer->Add(new TH1F("hCentrality","Ccentrality",100,0.,100.));
 
  fOutputContainer->Add(new TH2F("hTOF","cluster TOF",200,0.,20.,300,-3.e-6,6.e-6));
  
  
  fOutputContainer->Add(new TH2F("hCluNXZM1","Clu (X,Z), M1"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM2","Clu (X,Z), M2"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM3","Clu (X,Z), M3"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM1","Clu E(X,Z), M1"  ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM2","Clu E(X,Z), M2"  ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM3","Clu E(X,Z), M3"  ,64,0.5,64.5, 56,0.5,56.5));

  fOutputContainer->Add(new TH2F("hCluArea2M1","Clu (X,Z), M1"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluArea2M2","Clu (X,Z), M1"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluArea2M3","Clu (X,Z), M1"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluArea3M1","Clu (X,Z), M1"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluArea3M2","Clu (X,Z), M1"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluArea3M3","Clu (X,Z), M1"   ,64,0.5,64.5, 56,0.5,56.5));
  
  char cPID[4][5] ;
  snprintf(cPID[0],5,"All") ;
  snprintf(cPID[1],5,"Disp");
  snprintf(cPID[2],5,"CPV") ;
  snprintf(cPID[3],5,"Both"); 
  
 
  const Int_t nPt=400 ;
  const Double_t ptMax=40. ;
  const Int_t nM=400 ;
  const Double_t mMax=1. ;

  const Int_t nCenBin=5;
  for(Int_t cen=0; cen<nCenBin; cen++){

  for(Int_t iPID=0; iPID<4; iPID++){
    
    //Inclusive spectra
    fOutputContainer->Add(new TH1F(Form("hPhot_%s_cent%d",cPID[iPID],cen),"Spectrum of all reconstructed particles",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hPhot_Dist2_%s_cent%d",cPID[iPID],cen),"Spectrum of all reconstructed particles",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hPhot_Dist3_%s_cent%d",cPID[iPID],cen),"Spectrum of all reconstructed particles",nPt,0.,ptMax)) ;
  
    fOutputContainer->Add(new TH1F(Form("hPhot_Area1_%s_cent%d",cPID[iPID],cen),"Spectrum of all reconstructed particles",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hPhot_Area2_%s_cent%d",cPID[iPID],cen),"Spectrum of all reconstructed particles",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hPhot_Area3_%s_cent%d",cPID[iPID],cen),"Spectrum of all reconstructed particles",nPt,0.,ptMax)) ;

    fOutputContainer->Add(new TH1F(Form("hPhot_nStrictTagged_Area1_%s_cent%d",cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;

    for(Int_t itag=0; itag<18; itag++){
      fOutputContainer->Add(new TH1F(Form("hPhot_nTagged%d_Area1_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_nTagged%d_Area2_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_nTagged%d_Area3_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;
    }    
    for(Int_t kind=1; kind<33; kind*=2){
      fOutputContainer->Add(new TH1F(Form("hPhot_Isolation%d_%s_cent%d",kind,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_Isolation%d_Area1_%s_cent%d",kind,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_nStrictTagged_Isolation%d_Area1_%s_cent%d",kind,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;    
      fOutputContainer->Add(new TH1F(Form("hPhot_nTagged_Isolation%d_Area1_%s_cent%d",kind,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;
    }
  }
    
  
  fOutputContainer->Add(new TH1F(Form("hTaggedMult_cent%d",cen),"Spectrum of multiply tagged photons",nPt,0.,ptMax)) ;

  //Invariant mass distributions for fake corrections
  
  for(Int_t iPID=0; iPID<4; iPID++){
    fOutputContainer->Add(new TH2F(Form("hInvM_Re_Emin1_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Re_Emin2_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Re_Emin3_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Re_Emin1_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Re_Emin2_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Re_Emin3_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;

    fOutputContainer->Add(new TH2F(Form("hInvM_Mi_Emin1_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Mi_Emin2_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Mi_Emin3_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Mi_Emin1_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Mi_Emin2_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Mi_Emin3_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
  }  

 
  fOutputContainer->Add(new TH2F(Form("QA_Cone1_Tracks_cent%d",cen),"Two-photon inv. mass vs first photon pt",50,0.,50.,200,0.,100.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Cone2_Tracks_cent%d",cen),"Two-photon inv. mass vs first photon pt",50,0.,50.,200,0.,100.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Cone3_Tracks_cent%d",cen),"Two-photon inv. mass vs first photon pt",50,0.,50.,200,0.,100.)) ;
  }//centrality
  
  //MC  
  char partName[4][10] ;
  snprintf(partName[0],10,"gamma") ;
  snprintf(partName[1],10,"pi0");
  snprintf(partName[2],10,"eta") ;
  snprintf(partName[3],10,"omega"); 
  
  if(AliAnalysisManager::GetAnalysisManager()){
    AliMCEventHandler* mctruth = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    if(mctruth){
      
      fOutputContainer->Add(new TH1F("hMCConversionRadius","Clusters without label",600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCRecPi0Vtx","Secondary pi0s",100,0.,10.,600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCRecEtaVtx","Secondary etas",100,0.,10.,600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCRecOmegaVtx","Secondary etas",100,0.,10.,600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCRecEtaprVtx","Secondary etas",100,0.,10.,600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCRecK0sVtx","Secondary K0s",100,0.,10.,600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCRecK0lVtx","Secondary K0l",100,0.,10.,600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCGammaPi0MisConvR","Converted photons",400,0.,40.,600,0.,600.)) ;
 
  for(Int_t cen=0; cen<nCenBin; cen++){
  for(Int_t ipart=0; ipart<4; ipart++){  
    fOutputContainer->Add(new TH2F(Form("hMC_%s_vertex",partName[ipart]),"vertex",nPt,0.,ptMax,150,0.,600.)) ;
    fOutputContainer->Add(new TH1F(Form("hMC_all_%s",partName[ipart]),"Spectum (full rapifity)",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMC_unitEta_%s",partName[ipart]),"Spectum, |y|<0.15",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMC_rap_%s",partName[ipart]),"Rapidity",100,-5.,5.)) ;
    fOutputContainer->Add(new TH1F(Form("hMC_phi_%s",partName[ipart]),"Azimuthal angle",100,0.,TMath::TwoPi())) ;
  }

  
  fOutputContainer->Add(new TH2F(Form("LabelsNPrim_cent%d",cen),"vertex",nPt,0.,ptMax,20,0.,20.)) ;
  fOutputContainer->Add(new TH1F(Form("LabelsGamma_cent%d",cen),"vertex",nPt,0.,ptMax)) ;
  fOutputContainer->Add(new TH2F(Form("LabelsGammaE_cent%d",cen),"vertex",nPt,0.,ptMax,100,0.,2.)) ;
   
  fOutputContainer->Add(new TH1F(Form("hMCRecNoLabel_cent%d",cen),"Clusters without label",nPt,0.,ptMax)) ;
  fOutputContainer->Add(new TH1F(Form("hMCConversionRadius_cent%d",cen),"Clusters without label",600,0.,600.)) ;
  fOutputContainer->Add(new TH2F(Form("hMCRecPi0Vtx_cent%d",cen),"Secondary pi0s",100,0.,10.,600,0.,600.)) ;
  fOutputContainer->Add(new TH2F(Form("hMCRecEtaVtx_cent%d",cen),"Secondary etas",100,0.,10.,600,0.,600.)) ;
  fOutputContainer->Add(new TH2F(Form("hMCRecOmegaVtx_cent%d",cen),"Secondary etas",100,0.,10.,600,0.,600.)) ;
  fOutputContainer->Add(new TH2F(Form("hMCRecEtaprVtx_cent%d",cen),"Secondary etas",100,0.,10.,600,0.,600.)) ;
  fOutputContainer->Add(new TH2F(Form("hMCRecK0sVtx_cent%d",cen),"Secondary K0s",100,0.,10.,600,0.,600.)) ;
  fOutputContainer->Add(new TH2F(Form("hMCRecK0lVtx_cent%d",cen),"Secondary K0l",100,0.,10.,600,0.,600.)) ;
  fOutputContainer->Add(new TH2F(Form("hMCGammaPi0MisConvR_cent%d",cen),"Converted photons",400,0.,40.,600,0.,600.)) ;
  
  
  fOutputContainer->Add(new TH2F(Form("hMCGammaPi0PrimMgg_cent%d",cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
  fOutputContainer->Add(new TH2F(Form("hMCGammaPi0RecMgg_cent%d",cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
 
  for(Int_t iPID=0; iPID<4; iPID++){    
    fOutputContainer->Add(new TH1F(Form("hMCRecNoPrim_%s_cent%d",cPID[iPID],cen),"Reconstructed particles withour primary",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecGamma_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: Gamma",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecPhotConv_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: e+-",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecPi0_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecEta_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: eta",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecOmega_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: Gamma",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecEtapr_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: eta prime",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecPbar_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: bar(p)",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecNbar_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: bar(n)",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecPipm_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: pipm",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecN_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: n",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecP_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: p",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecKpm_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: K+-",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecK0s_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: K0s",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecK0l_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: K0l",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecUnknownCh_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: K0l",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecUnknownNeu_%s_cent%d",cPID[iPID],cen),"Reconstructed particles with primary: K0l",nPt,0.,ptMax)) ;

    //Decay photons
    fOutputContainer->Add(new TH1F(Form("hMCRecGammaDir_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, no primary",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecGammaPi0_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecGammaEta_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from eta",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecGammaOmega_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from omega",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecGammaOther_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from other",nPt,0.,ptMax)) ;
    
    //Pi0 decay photons
    fOutputContainer->Add(new TH1F(Form("hMCRecGammaPi0Dalitz_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCRecGammaPi0NoStack_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi02Gamma_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0Rec_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0RecSoft_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0RecCells_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;

    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0Tagged_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0FakeTagged_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0TrueTagged_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0MultyTagged_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0MisGeo_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0MisGeoFA1_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0MisGeoFA2_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0MisGeoFA3_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0MisConv_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0MisEmin_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0MisOther_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0MisFakePrim_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0MisFoundPrim_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hMCGammaPi0MisNPhot_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;  

    //all clusters fake tagged
    fOutputContainer->Add(new TH1F(Form("hMCAllFakeTagged_%s_cent%d",cPID[iPID],cen),"Reconstructed gammas, from pi0",nPt,0.,ptMax)) ;

    fOutputContainer->Add(new TH2F(Form("hMC_InvM_Re_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hMC_InvM_Re_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
   
  }
  fOutputContainer->Add(new TH1F(Form("hMCGammaPi0MisPartner_cent%d",cen),"Spectrum of missed partners",nPt,0.,ptMax)) ;
  fOutputContainer->Add(new TH2F(Form("hMCGammaPi0MisPartnerEtaPhi_cent%d",cen),"Spectrum of missed partners",100,-0.2,0.2,100,4.5,5.6)) ;
    }
    }
  }//If MC handler exists...
  
  for(Int_t i=0;i<1;i++)
    for(Int_t j=0;j<5;j++)
      fPHOSEvents[i][j]=0x0 ;    //Container for PHOS photons
  

  PostData(1, fOutputContainer);


}

//______________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::UserExec(Option_t *) 
{
  //Select events
  //Select photons
  //Fill QA histograms
  //Fill Tagging histogsms
  
  // Event selection flags

  FillHistogram("hSelEvents",1) ;
    
  AliVEvent* event = (AliVEvent*)InputEvent();
  if(!event){
    AliDebug(1,"No event") ;
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents",2) ;
  
  //read geometry if not read yet
  if(fPHOSgeom==0)
    InitGeometry() ;
  
  //MC stack init
  fStack=0x0 ;
  if(AliAnalysisManager::GetAnalysisManager()){
    AliMCEventHandler* mctruth = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    if(mctruth)
      fStack = mctruth->MCEvent()->Stack();
  }
  
  
  // Checks if we have a primary vertex
  // Get primary vertices form AOD

/*  
  if(event->GetPrimaryVertexTracks()->GetNContributors()>0)
    eventVtxExist    = kTRUE;
  else 
    if(event->GetPrimaryVertexSPD()->GetNContributors()>0)
      eventVtxExist    = kTRUE;
*/

  Double_t vtx5[3];
  vtx5[0] = event->GetPrimaryVertex()->GetX();
  vtx5[1] = event->GetPrimaryVertex()->GetY();
  vtx5[2] = event->GetPrimaryVertex()->GetZ();

  FillHistogram("hNvertexTracks",event->GetPrimaryVertex()->GetNContributors());
  FillHistogram("hZvertex"      ,vtx5[2]);
  if (TMath::Abs(vtx5[2]) > 10. ){
    PostData(1, fOutputContainer);
    return ;
  }
    
  FillHistogram("hSelEvents",3) ;
  
  
//  if (event->IsPileupFromSPD()){
//    PostData(1, fOutputContainer);
//    return ;
//  }

  FillHistogram("hSelEvents",4) ;

  //centrality
  AliCentrality *centrality = event->GetCentrality();
  if( centrality )
    fCentrality=centrality->GetCentralityPercentile("V0M");
  else {
    AliError("Event has 0x0 centrality");
    fCentrality = -1.;
  }
  FillHistogram("hCentrality",fCentrality) ;

  if(fCentrality<0. || fCentrality>=100.){
    PostData(1, fOutputContainer);
    return ;
  }
  fCentBin = (Int_t)(fCentrality/20.) ; 

  FillHistogram("hSelEvents",5) ;
  
  
  //Vtx class z-bin
  Int_t zvtx = 0 ; 

  //Calculate charged multiplicity
  Int_t trackMult = 0;
  if(fTrackEvent)
    fTrackEvent->Clear() ;
  else
    fTrackEvent = new TClonesArray("AliAODPWG4Particle",event->GetNumberOfTracks()) ;
  for (Int_t i=0;i<event->GetNumberOfTracks();++i) {
    AliVParticle *track = event->GetTrack(i) ;
    if(TMath::Abs(track->Eta())< 0.9){
      if(trackMult>=fTrackEvent->GetSize())
	fTrackEvent->Expand(2*trackMult) ;
      new ((*fTrackEvent)[trackMult]) AliAODPWG4Particle(track->Px(),track->Py(),track->Pz(),track->P());
      trackMult++;
    }
  }
  FillHistogram("hTrackMult",fCentrality,trackMult+0.5) ;

  if(!fPHOSEvents[zvtx][fCentBin]) 
    fPHOSEvents[zvtx][fCentBin]=new TList() ;
  fCurrentMixedList = fPHOSEvents[zvtx][fCentBin] ;

   const Double_t rcut=1. ; //cut on vertex to consider particle as "primary" 
 
  //---------Select photons-------------------
  Int_t multClust = event->GetNumberOfCaloClusters();
  if(!fPHOSEvent)
    fPHOSEvent   = new TClonesArray("AliAODPWG4Particle",multClust);
  else
    fPHOSEvent->Clear() ;
  Int_t inList = 0; //counter of caloClusters

  for (Int_t i=0; i<multClust; i++) {
    AliVCluster * clu = event->GetCaloCluster(i);
    
    if(!clu->IsPHOS())
      continue ; 
    
    if(clu->E()<0.1) 
      continue;

    if(clu->GetNCells()<3)
      continue ;          
    
    if(clu->GetM02()<0.2) 
      continue ;          
    
    FillHistogram("hTOF",clu->E(),clu->GetTOF()) ;
    if(TMath::Abs(clu->GetTOF())>100.e-9) 
      continue ;          
    
    if(clu->GetDistanceToBadChannel()<2.5)
      continue ;

    Float_t pos[3] ;
    clu->GetPosition(pos) ;
    Int_t fidArea=GetFiducialArea(pos) ;
    if(fidArea==0) //Bad cell
      continue; 

    
    TVector3 global1(pos) ;
    Int_t relId[4] ;
    fPHOSgeom->GlobalPos2RelId(global1,relId) ;
    Int_t mod  = relId[0] ;
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;
    
    FillHistogram(Form("hCluNXZM%d",mod),cellX,cellZ,1.);
    FillHistogram(Form("hCluEXZM%d",mod),cellX,cellZ,clu->E());
    if(fidArea>1){
      FillHistogram(Form("hCluArea2M%d",mod),cellX,cellZ,1.);
      if(fidArea>2){
         FillHistogram(Form("hCluArea3M%d",mod),cellX,cellZ,1.);
      }
    }
    
    TLorentzVector momentum ;
    clu->GetMomentum(momentum, vtx5);
    AliAODPWG4Particle *p = new ((*fPHOSEvent)[inList]) AliAODPWG4Particle(momentum.Px(),momentum.Py(),momentum.Pz(),clu->E() );
    inList++;

    Int_t isolation = EvalIsolation(&momentum) ;
    p->SetBtag(isolation) ;
    
    p->SetDistToBad((Int_t)(1.+clu->GetDistanceToBadChannel()/2.2));
    
    p->SetTag(0); //Strict PID pi0 partner not found
    p->SetTagged(kFALSE);   //Reconstructed pairs found
    Bool_t sure = 0;
    p->SetFiducialArea(fidArea) ;

    if(fStack){    
     //This is primary entered PHOS
     Int_t primLabel=FindPrimary(clu,sure) ;
     //Look what particle left vertex
     if(primLabel>-1){
       TParticle * prim = fStack->Particle(primLabel) ;
       Int_t iparent=primLabel;
       TParticle * parent = prim;
       while((parent->R() > rcut) && (iparent>-1)){
         iparent=parent->GetFirstMother();
         parent=fStack->Particle(iparent);
      }
       p->SetCaloLabel(primLabel,iparent); //This and partner cluster
     }
     else{
       p->SetCaloLabel(-1,-1); //This and partner cluster
     }
       
    }
    else{  
     //This is primary at vertex(R<1cm)
     p->SetCaloLabel(-1,-1); //This and partner cluster
    }
    p->SetLabel(i); //Cluster index    
    //PID criteria
//    p->SetDispBit(clu->Chi2()) ;
    p->SetDispBit(TestLambda(clu->E(),clu->GetM20(),clu->GetM02())) ;
    p->SetTOFBit(TestTOF(clu->GetTOF(),clu->E())) ;
    p->SetChargedBit(clu->GetEmcCpvDistance()>10.) ;      
  }
  
  FillMCHistos() ;
  FillTaggingHistos() ;

  //Remove old events
  fCurrentMixedList->AddFirst(fPHOSEvent);
  fPHOSEvent=0x0 ;
  if(fCurrentMixedList->GetSize() > 10){
    TClonesArray *tmp = static_cast <TClonesArray*> (fCurrentMixedList->Last());
    fCurrentMixedList->Remove(tmp);
    delete tmp;
  }
  
  PostData(1, fOutputContainer);

}
//________________________________________________
void AliAnalysisTaskTaggedPhotons::FillMCHistos(){
   
  //MC info about this particle
  if(!fStack)
    return ;
  
  const Double_t rcut=1. ; //cut on vertex to consider particle as "primary" 
  AliAODEvent* event = (AliAODEvent*)InputEvent();
  
  //Fill spectra of primary particles
  char partName[10] ;
  for(Int_t i=0;i<fStack->GetNtrack();i++){
    TParticle* particle =  fStack->Particle(i);
    //Primary particle
    Double_t r=particle->R() ;

    switch(particle->GetPdgCode()){
      case 111:  snprintf(partName,10,"pi0") ; 
                 break ;
      case 221:  snprintf(partName,10,"eta") ;
                 break ;
      case 223:  snprintf(partName,10,"omega") ;
                 break ;
      case 22:   snprintf(partName,10,"gamma") ;
                 break ;
      default:   continue ;
    }
    
    Double_t pt = particle->Pt() ;
    //Distribution over vertex
    FillHistogram(Form("hMC_%s_vertex",partName),pt,r) ;

    if(r >rcut)
      continue ;
   
    //Total number of pi0 with creation radius <1 cm
    Double_t w = PrimaryParticleWeight(particle) ;  
    FillHistogram(Form("hMC_all_%s",partName),pt,w) ;
    if(TMath::Abs(particle->Y())<0.15){
      FillHistogram(Form("hMC_unitEta_%s",partName),pt,w) ;
    }

    FillHistogram(Form("hMC_rap_%s",partName),particle->Y(),w) ;
    
    Double_t phi=particle->Phi() ;
    while(phi<0.)phi+=TMath::TwoPi() ;
    while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
    FillHistogram(Form("hMC_phi_%s",partName),phi,w) ;
  }
 
  
  //Clussify reconstructed clusters
  //First - photons (from vertex) and contaminations
  //Second - photons from different sources
  //Third - photons from pi0s - missed for different reasons
  
  const Int_t n=fPHOSEvent->GetEntriesFast() ;
  for(Int_t i=0;i<n;i++){
    AliAODPWG4Particle *p = static_cast<AliAODPWG4Particle*>(fPHOSEvent->At(i));
    Int_t iprim=p->GetCaloLabel(0) ; //Particle entered PHOS
    if(iprim<0){ //No label!
      FillHistogram("hMCRecNoLabel",p->Pt());
      continue ;
    }     
    
    TParticle * prim = fStack->Particle(iprim) ;
      
    //Look what particle left vertex
    Int_t iparent=p->GetCaloLabel(1); //Particle left vertex
    TParticle * parent =fStack->Particle(iparent); 

    Int_t grandpaPDG=-1 ;
    if(parent){
      grandpaPDG=parent->GetPdgCode() ;
    }
    switch(grandpaPDG){
    case -1: //no primary
	  FillPIDHistograms("hMCRecNoPrim",p);
	  break ;	  
    case 22:  
	  FillPIDHistograms("hMCRecGamma",p);
	  break ;
    case  11:
    case -11: //electron/positron conversion
	  FillPIDHistograms("hMCRecPhotConv",p);  //Reconstructed with photon from conversion primary
	  FillHistogram("hMCConversionRadius",prim->R());
	  break ;
    case 111: //Pi0 decay           //Primary decay photon (as in MC)
 	  FillPIDHistograms("hMCRecPi0",p);
          //Strange, vertex of pi0?
 	  FillHistogram("hMCRecPi0Vtx",p->Pt(),parent->R());
	  break ;
    case 221: //eta decay
	  FillPIDHistograms("hMCRecEta",p);
          //Strange, vertex of eta?
 	  FillHistogram("hMCRecEtaVtx",p->Pt(),parent->R());
	  break ;  
    case 223: //omega meson decay
	  FillPIDHistograms("hMCRecOmega",p);
          //Strange, vertex of omega?
 	  FillHistogram("hMCRecOmegaVtx",p->Pt(),parent->R());
	  break ;
    case 331: //eta' decay
	  FillPIDHistograms("hMCRecEtapr",p);
          //Strange, vertex of eta'?
 	  FillHistogram("hMCRecEtaprVtx",p->Pt(),parent->R());
	  break ;	  
    case -2212:
	  FillPIDHistograms("hMCRecPbar",p);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
    case -2112: //antineutron & antiproton 
	  FillPIDHistograms("hMCRecNbar",p);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
    case  211:
    case -211:
	  FillPIDHistograms("hMCRecPipm",p);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
    case 2112:
	  FillPIDHistograms("hMCRecN",p);  //Reconstructed with photon from antibaryon annihilation
	  break ;	        
    case 2212:
	  FillPIDHistograms("hMCRecP",p);  //Reconstructed with photon from antibaryon annihilation
	  break ;	        
    case  321:
    case -321:
	  FillPIDHistograms("hMCRecKpm",p);  //Reconstructed with photon from conversion primary
	  break ;
    case 130:
	  FillHistogram("hMCRecK0sVtx",p->Pt(),prim->R());
	  FillPIDHistograms("hMCRecK0s",p);  //Reconstructed with photon from conversion primary
	  break ;
    case 310: 	  
	  FillHistogram("hMCRecK0lVtx",p->Pt(),prim->R());
	  FillPIDHistograms("hMCRecK0l",p);  //Reconstructed with photon from conversion primary
	  break ;
    default:  
          if(parent->GetPDG()->Charge()!=0)
            FillPIDHistograms("hMCRecUnknownCh",p);
	  else
            FillPIDHistograms("hMCRecUnknownNeu",p);
	  printf("Unknown PDG code: %d \n",grandpaPDG) ;
    }  
//Put here filling histograms vs PID...    
    
        
    //Now classify decay photons
    if(grandpaPDG==22){
      //Separate gammas from pi0s, eta, omega, eta', e+-, other, no primary
      Int_t igrandpa = -1 ;
      if(parent)
        igrandpa=parent->GetFirstMother();
      Int_t decayPDG=-1 ; 
      TParticle * decay=0x0 ;
      if(igrandpa>-1){
	decay=fStack->Particle(igrandpa) ;
	decayPDG=decay->GetPdgCode() ;
      }
      switch(decayPDG){
	case -1 :// Direct gamma?
            FillPIDHistograms("hMCRecGammaDir",p);
	    break ;
	case 111 :// pi0 decays
            FillPIDHistograms("hMCRecGammaPi0",p);
	    break ;
	case 221 :// eta decays
            FillPIDHistograms("hMCRecGammaEta",p);
	    break ;
	case 223 :// omega decays
            FillPIDHistograms("hMCRecGammaOmega",p);
	    break ;
	default:
            FillPIDHistograms("hMCRecGammaOther",p);
      }	    


      //Classify photons from pi0 decay
      if(decayPDG==111){
	//Dalitz decay
	if(decay->GetNDaughters()>2){
          FillPIDHistograms("hMCRecGammaPi0Dalitz",p);
	  continue ;
	}
	Int_t ipartner = decay->GetFirstDaughter(); 
	if(ipartner==iparent){ //look other
          ipartner = decay->GetLastDaughter(); 
	}
	//No partner in Stack
	if(ipartner==-1){
          FillPIDHistograms("hMCRecGammaPi0NoStack",p) ;
	}
        else{
	  //	  
          FillPIDHistograms("hMCGammaPi02Gamma",p) ;
	  TParticle * partner = fStack->Particle(ipartner);
	  //Check if partner is registered and made correct mass
	  Double_t mPrim=(parent->Energy()+partner->Energy())*(parent->Energy()+partner->Energy())-
	                 (parent->Px()+partner->Px())*(parent->Px()+partner->Px()) -
	                 (parent->Py()+partner->Py())*(parent->Py()+partner->Py()) -
	                 (parent->Pz()+partner->Pz())*(parent->Pz()+partner->Pz()) ;
          if(mPrim>0.)mPrim=TMath::Sqrt(mPrim) ;
	  FillHistogram(Form("hMCGammaPi0PrimMgg_cent%d",fCentBin),mPrim,p->Pt()) ;	  
	  
	  //find corresponding PWG4particle
	  Bool_t isPartnerRec=kFALSE ;
	  Bool_t isTagged=kFALSE ;
          for(Int_t j=0;j<n;j++){
	    if(j==i) 
	      continue ;
            AliAODPWG4Particle *p2 = static_cast<AliAODPWG4Particle*>(fPHOSEvent->At(j));
	    Double_t invMass = p->GetPairMass(p2);      
            if(p2->GetCaloLabel(1)==ipartner){ //partner
              isPartnerRec=kTRUE ;
  	      FillHistogram("hMCGammaPi0RecMgg",invMass,p->Pt()) ;
              FillPIDHistograms("hMCGammaPi0Rec",p) ;
	  
	      //estimate proportion passed PI0 cut
              if(IsInPi0Band(invMass,p->Pt(),2)){ //Check type!!!!!!!!!!!!!!
                FillPIDHistograms("hMCGammaPi0TrueTagged",p) ;		  
                if(isTagged)
                  FillPIDHistograms("hMCGammaPi0MultyTagged",p) ;		  
		isTagged=kTRUE ;
	      }
	    }
	    else{
	      if(IsInPi0Band(invMass,p->Pt(),2)){//Check type!!!!!!!!!!!!!!
                FillPIDHistograms("hMCGammaPi0FakeTagged",p) ;
                if(isTagged)
                  FillPIDHistograms("hMCGammaPi0MultyTagged",p) ;		  
		isTagged=kTRUE ;		
	      }
	    }
	  }
          if(isTagged)          
            FillPIDHistograms("hMCGammaPi0Tagged",p) ;
	  
	  //If no PWG4particle trace the reason
	  if(!isPartnerRec){
             //conversion
	      if(partner->GetPdgCode()==22){ 
		Bool_t isPartnerLost=kFALSE; //If partner is lost for some reason
		
	        //did not hit PHOS	
		Int_t modulenum ;
		Double_t ztmp=0.,xtmp=0. ;
		Bool_t impact=fPHOSgeom->ImpactOnEmc(partner,modulenum,ztmp,xtmp) ;
		//Check Bad Map
		Int_t fidArea=0 ;                  
		TVector3 partnerGlobaPos ;
		if(impact){
                  fPHOSgeom->Local2Global(modulenum,xtmp, ztmp,partnerGlobaPos) ;
                  Float_t pos[3] ;
		  partnerGlobaPos.GetXYZ(pos);
                  fidArea=GetFiducialArea(pos) ;
                }

		if(!impact || (fidArea==0) ){ //this photon cannot hit PHOS		  
		  FillPIDHistograms("hMCGammaPi0MisGeo",p);
		  Int_t iFidArea = p->GetFiducialArea(); 
		  if(iFidArea>0){
		    FillPIDHistograms("hMCGammaPi0MisGeoFA1",p) ;  //Spectrum of tagged with missed partner
		    if(iFidArea>1){
		      FillPIDHistograms("hMCGammaPi0MisGeoFA2",p) ;  //Spectrum of tagged with missed partner
		      if(iFidArea>2){
			FillPIDHistograms("hMCGammaPi0MisGeoFA3",p) ;  //Spectrum of tagged with missed partner
		      }
		    }
		  }
		  isPartnerLost=kTRUE;
		}

		//this photon is converted before it is registered
		if(!isPartnerLost){
		  if(partner->GetNDaughters()>0 && 
		     fStack->Particle(partner->GetFirstDaughter())->R()<450.){  
		    FillPIDHistograms("hMCGammaPi0MisConv",p);
		    FillHistogram("hMCGammaPi0MisConvR",p->Pt(),fStack->Particle(partner->GetFirstDaughter())->R()) ;  //Spectrum of tagged with missed partner
		    isPartnerLost=kTRUE;
		  }
		}
		
        	// too soft	
		if(!isPartnerLost && 
		   partner->Energy()<0.3){ //energy is not enough to be registered by PHOS
		  FillPIDHistograms("hMCGammaPi0MisEmin",p);  
		  isPartnerLost=kTRUE;
		}
		//Other reason
		if(!isPartnerLost){
		  FillHistogram("hMCGammaPi0MisPartner",partner->Pt());
		  FillHistogram("hMCGammaPi0MisPartnerEtaPhi",partner->Eta(),partner->Phi());	
		  //May be overlap?
		  Bool_t fakePrimary=kFALSE ;
                  Int_t multClust = event->GetNumberOfCaloClusters();
                  for (Int_t iclu=0; iclu<multClust; iclu++) {
                    AliAODCaloCluster * clu = event->GetCaloCluster(iclu);
                    Float_t pos[3] ;
                    clu->GetPosition(pos) ;
                    TVector3 global1(pos) ;
		    Double_t d=(partnerGlobaPos-global1).Mag();
		    if( d<4 ){//same cluster
                      //check is partner is in list?
                      Int_t nlab=clu->GetNLabels() ;
                      for(Int_t ilab=0;  ilab<nlab;  ilab++){
		        Int_t labelA = clu->GetLabelAt(ilab) ;
		        while(labelA>-1){
			  if(labelA==ipartner){
			    if(clu->E()>0.3)
         	              FillPIDHistograms("hMCGammaPi0RecSoft",p);
			    else
			      if(clu->GetNCells()<3)
         	                FillPIDHistograms("hMCGammaPi0RecCells",p);
			      else
       	                        FillPIDHistograms("hMCGammaPi0MisFoundPrim",p);
			    break ; 
			  }
                          labelA=fStack->Particle(labelA)->GetFirstMother() ;
		        }
		      }                      
                      fakePrimary=kTRUE ;
                      break ;
		    }
		  }//end of loop over clusters
                  //No cluster in this region
                  if(fakePrimary)
		    FillPIDHistograms("hMCGammaPi0MisFakePrim",p);
		  else
		    FillPIDHistograms("hMCGammaPi0MisOther",p);
		  
		}
	      }//Partner - photon
	      else{//partner not photon
		FillPIDHistograms("hMCGammaPi0MisNPhot",p);                
	      }
	    }
	  }
	}
      }
    } //PHOS clusters 
   
   
   //Fill histograms with all clusters fake tagging
   for(Int_t i=0;i<n-1;i++){
     AliAODPWG4Particle *p = static_cast<AliAODPWG4Particle*>(fPHOSEvent->At(i));

     Int_t ipartner=-999 ; // this will be partner if this cluster from pi0 decay
     //Look what particle left vertex
     Int_t iparent=p->GetCaloLabel(1); //Particle left vertex
     TParticle * parent =fStack->Particle(iparent); 
     if(parent){	
       Int_t grandpaPDG=parent->GetPdgCode() ;
       if(grandpaPDG==22){
         Int_t igrandpa =parent->GetFirstMother();
         if(igrandpa>-1){
           TParticle * decay=fStack->Particle(igrandpa) ;
	   Int_t decayPDG=decay->GetPdgCode() ;
           if(decayPDG==111){
	     //Dalitz decay
	     if(decay->GetNDaughters()==2){     
	       ipartner = decay->GetFirstDaughter(); 
	       if(ipartner==iparent){ //look other
                 ipartner = decay->GetLastDaughter(); 
	       }
	       if(ipartner==-1) //no partner
                 ipartner=-999 ;
   	     }
	   }
	 }
       }
     }
     for(Int_t j=i+1;j<n;j++){
       AliAODPWG4Particle *p2 = static_cast<AliAODPWG4Particle*>(fPHOSEvent->At(j));
       if(p2->GetCaloLabel(1)!=ipartner){ //partner
         Double_t invMass = p->GetPairMass(p2);      
	 if(IsInPi0Band(invMass,p->Pt(),2)){
             FillPIDHistograms("hMCAllFakeTagged",p) ;
	 }
	 if(IsInPi0Band(invMass,p2->Pt(),2)){
             FillPIDHistograms("hMCAllFakeTagged",p2) ;
         }
       }
     }
   } 
   
}
//________________________________________________
void AliAnalysisTaskTaggedPhotons::FillTaggingHistos(){
  //Fill all necessary histograms


 //First fill invariant mass distrubutions and mark tagged photons
  //Invariant Mass analysis
  const Int_t n=fPHOSEvent->GetEntriesFast() ;
  for(Int_t i=0;i<n-1;i++){
    AliAODPWG4Particle *p1 = static_cast<AliAODPWG4Particle*>(fPHOSEvent->At(i));
    for(Int_t j = i+1 ; j < n ; j++) {
      AliAODPWG4Particle * p2 = static_cast<AliAODPWG4Particle*>(fPHOSEvent->At(j));
      
      Double_t invMass = p1->GetPairMass(p2);   
      
      if(p2->E()>0.1){
        FillPIDHistograms("hInvM_Re_Emin1",p1,invMass) ;
        if(p2->E()>0.2){
          FillPIDHistograms("hInvM_Re_Emin2",p1,invMass) ;
          if(p2->E()>0.3){
            FillPIDHistograms("hInvM_Re_Emin3",p1,invMass) ;
	  }
	}
      }
	
      if(p1->E()>0.1){
        FillPIDHistograms("hInvM_Re_Emin1",p2,invMass) ;
        if(p1->E()>0.2){
          FillPIDHistograms("hInvM_Re_Emin2",p2,invMass) ;
          if(p1->E()>0.3){
            FillPIDHistograms("hInvM_Re_Emin3",p2,invMass) ;
	  }
	}
      }
      if(TestPID(3, p2)){
        if(p2->E()>0.1){
          FillPIDHistograms("hInvM_Re_Emin1_Strict",p1,invMass) ;
          if(p2->E()>0.2){
            FillPIDHistograms("hInvM_Re_Emin2_Strict",p1,invMass) ;
            if(p2->E()>0.3){
              FillPIDHistograms("hInvM_Re_Emin3_Strict",p1,invMass) ;
	    }
	  }
	}
      }
      if(TestPID(3, p1)){
        if(p1->E()>0.1){
          FillPIDHistograms("hInvM_Re_Emin1_Strict",p2,invMass) ;
          if(p1->E()>0.1){
            FillPIDHistograms("hInvM_Re_Emin2_Strict",p2,invMass) ;
            if(p1->E()>0.3){
              FillPIDHistograms("hInvM_Re_Emin3_Strict",p2,invMass) ;
	    }
	  }
	}
      }
      if(IsSamePi0(p1,p2)){
        FillPIDHistograms("hMC_InvM_Re",p1,invMass) ;
        FillPIDHistograms("hMC_InvM_Re",p2,invMass) ;
        if(TestPID(3, p2)){
          FillPIDHistograms("hMC_InvM_Re_Strict",p1,invMass) ;
        }
        if(TestPID(3, p1)){
          FillPIDHistograms("hMC_InvM_Re_Strict",p2,invMass) ;
        }
      }

      //Tagging: 1,2,3 sigma
      //Emin=100,200,300 Mev
      //IsInPi0Band(mass, Ptphoton, type Emin cut
      Int_t tag1=0 ;
      for(Int_t eminType=0; eminType<3; eminType++){
        if(p2->E()>0.1*(eminType+1)){
	  //Set corresponding bit
	  Double_t nsigma = IsInPi0Band(invMass,p1->Pt(),eminType) ; //in band with n sigmas
	  for(Int_t isigma=0; isigma<3; isigma++){
  	    if(nsigma<1+isigma){
   	      tag1|= (1 << (3*eminType+isigma)) ;
	      if(TestPID(3, p2))
   	        tag1|= (1 << (3*eminType+isigma+9)) ;
	    }
	  }
	}
      }
      p1->SetTag(tag1) ;
      Int_t tag2=0 ;
      for(Int_t eminType=0; eminType<3; eminType++){
        if(p1->E()>0.1*(eminType+1)){
	  //Set corresponding bit
	  Double_t nsigma = IsInPi0Band(invMass,p2->Pt(),eminType) ; //in band with n sigmas
	  for(Int_t isigma=0; isigma<3; isigma++){
  	    if(nsigma<1+isigma){
   	      tag2|= (1 << (3*eminType+isigma)) ;
	      if(TestPID(3, p2))
   	        tag2|= (1 << (3*eminType+isigma+9)) ;
	    }
	  }
	}
      }
      p2->SetTag(tag2) ;
      
      if(tag1 & (1<<7)){ //2 sigma, Emin=0.3: default tagging
        if(p1->IsTagged()){//Multiple tagging
          FillHistogram(Form("hTaggedMult_cent%d",fCentBin),p1->Pt());
        }  
        p1->SetTagged(kTRUE) ;
      }
      if(tag2 & (1<<7)){ //2 sigma, Emin=0.3: default tagging
        if(p2->IsTagged()){//Multiple tagging
          FillHistogram(Form("hTaggedMult_cent%d",fCentBin),p2->Pt());
        }  
        p2->SetTagged(kTRUE) ;
      }      
    }
  }
  
  
  
  //Single particle histgams
  for(Int_t i=0;i<n;i++){
    AliAODPWG4Particle *p = static_cast<AliAODPWG4Particle*>(fPHOSEvent->At(i));

    Int_t isolation = p->GetBtag();

    //Inclusive spectra
    FillPIDHistograms("hPhot",p) ;
      
    if(p->DistToBad()>1){
      FillPIDHistograms("hPhot_Dist2",p) ;
      if(p->DistToBad()>2){
        FillPIDHistograms("hPhot_Dist3",p) ;
      }
    }
      
      
    for(Int_t kind=1; kind<33; kind*=2){
      if((isolation&kind)){
        FillPIDHistograms(Form("hPhot_Isolation%d",kind),p) ;
      }
    }
      
    Int_t iFidArea = p->GetFiducialArea(); 
    if(iFidArea>0){
      FillPIDHistograms("hPhot_Area1",p) ;
      for(Int_t kind=1; kind<33; kind*=2){
        if((isolation&kind)){
          FillPIDHistograms(Form("hPhot_Isolation%d_Area1",kind),p) ;
	}
      }

      //Fill taggings with 
      //3 Emin cuts
      //Default Emin, 1,2,3 sigmas
      //strict and loose PID cut on partner
      Int_t tag=p->GetTag() ;
      for(Int_t ibit=0; ibit<18; ibit++){
        if((tag & (1<<ibit))==0){ 
          FillPIDHistograms(Form("hPhot_nTagged%d_Area1",ibit),p) ;
//        for(Int_t kind=1; kind<33; kind*=2){
//          if((isolation&kind)){
//            FillPIDHistograms(Form("hPhot_nStrictTagged_Isolation%d_Area1",kind),p) ;
//          }
//	  }
	}
      }

      if(iFidArea>1){
        FillPIDHistograms("hPhot_Area2",p) ;
        for(Int_t ibit=0; ibit<18; ibit++){
          if((tag & (1<<ibit))==0){ 
            FillPIDHistograms(Form("hPhot_nTagged%d_Area2",ibit),p) ;
	  }
        }
	if(iFidArea>2){
          FillPIDHistograms("hPhot_Area3",p) ;
          for(Int_t ibit=0; ibit<18; ibit++){
            if((tag & (1<<ibit))==0){ 
              FillPIDHistograms(Form("hPhot_nTagged%d_Area3",ibit),p) ;
	    }
	  }
	}
      }
    }
  } 
  
   //Fill Mixed InvMass distributions:
  TIter nextEv(fCurrentMixedList) ;
  for(Int_t i=0;i<n;i++){
    AliAODPWG4Particle *p1 = static_cast<AliAODPWG4Particle*>(fPHOSEvent->At(i));
    while(TClonesArray * event2 = static_cast<TClonesArray*>(nextEv())){
      Int_t nPhotons2 = event2->GetEntriesFast() ;
      for(Int_t j=0; j < nPhotons2 ; j++){
        AliAODPWG4Particle * p2 = static_cast<AliAODPWG4Particle*>(event2->At(j)) ;
        Double_t invMass = p1->GetPairMass(p2);

	if(p2->E()>0.1){
          FillPIDHistograms("hInvM_Mi_Emin1",p1,invMass) ;
	  if(p2->E()>0.2){
            FillPIDHistograms("hInvM_Mi_Emin2",p1,invMass) ;
	    if(p2->E()>0.3){
              FillPIDHistograms("hInvM_Mi_Emin3",p1,invMass) ;
	    }
	  }
	}
        if(TestPID(3, p2)){
  	  if(p2->E()>0.1){
            FillPIDHistograms("hInvM_Mi_Emin1_Strict",p1,invMass) ;
    	    if(p2->E()>0.2){
              FillPIDHistograms("hInvM_Mi_Emin2_Strict",p1,invMass) ;
  	      if(p2->E()>0.3){
                FillPIDHistograms("hInvM_Mi_Emin3_Strict",p1,invMass) ;
	      }
	    }
	  }
	}
	
	if(p1->E()>0.1){
          FillPIDHistograms("hInvM_Mi_Emin1",p2,invMass) ;
  	  if(p1->E()>0.2){
            FillPIDHistograms("hInvM_Mi_Emin2",p2,invMass) ;
	    if(p1->E()>0.3){
              FillPIDHistograms("hInvM_Mi_Emin3",p2,invMass) ;
	    }
	  }
	}
        if(TestPID(3, p1)){
	  if(p1->E()>0.1){
            FillPIDHistograms("hInvM_Mi_Emin1_Strict",p2,invMass) ;
	    if(p1->E()>0.2){
              FillPIDHistograms("hInvM_Mi_Emin2_Strict",p2,invMass) ;
	      if(p1->E()>0.3){
                FillPIDHistograms("hInvM_Mi_Emin3_Strict",p2,invMass) ;
	      }
	    }
	  }
	}
      }
    }
  } 
 
  
}

//______________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::Init()
{
  // Intialisation of parameters
}

//______________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  if (fDebug > 1) Printf("Terminate()");
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskTaggedPhotons::IsInPi0Band(Double_t m, Double_t pt, Int_t /*type*/)const
{
  //Parameterization of the pi0 peak region
  Double_t mpi0mean =  1.33259e-01 -  2.910e-03 * TMath::Exp(-pt/3.616) ;
//  Double_t mpi0mean = 0.136 ; //fPi0MeanP0 + fPi0MeanP1 * pt + fPi0MeanP2 * pt*pt + fPi0MeanP3 * pt*pt*pt;

    Double_t mpi0sigma=TMath::Sqrt(4.17927e-03*4.17927e-03+2.81581e-03*2.81581e-03/pt+3.59218e-04*3.59218e-04*pt*pt) ;
//  Double_t mpi0sigma= 0.005; //TMath::Sqrt(fPi0SigmaP0 * fPi0SigmaP0 / pt + fPi0SigmaP1 * fPi0SigmaP1 + fPi0SigmaP2 * fPi0SigmaP2 / pt / pt);
 
  return (m>mpi0mean-2*mpi0sigma && m<mpi0mean+2*mpi0sigma) ;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskTaggedPhotons::IsSamePi0(const AliAODPWG4Particle *p1, const AliAODPWG4Particle *p2)const{
  //Looks through parents and finds if there was commont pi0 among ancestors

  if(!fStack)
    return kFALSE ; //can not say anything

  Int_t prim1 = p1->GetCaloLabel(0);
  while(prim1!=-1){ 
    Int_t prim2 = p2->GetCaloLabel(0);
    while(prim2!=-1){ 
      if(prim1==prim2){
        if(fStack->Particle(prim1)->GetPdgCode()==111)
          return kTRUE ;
        else
          return kFALSE ;
      }
      prim2=fStack->Particle(prim2)->GetFirstMother() ;
    }
    prim1=fStack->Particle(prim1)->GetFirstMother() ;
  }
  return kFALSE ;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskTaggedPhotons::GetFiducialArea(const Float_t * position)const{
  //calculates in which kind of fiducial area photon hit

  TVector3 global1(position) ;
  Int_t relId[4] ;
  fPHOSgeom->GlobalPos2RelId(global1,relId) ;
  Int_t mod  = relId[0] ;
  Int_t cellX = relId[2];
  Int_t cellZ = relId[3] ;

  if(fPHOSBadMap[mod] && fPHOSBadMap[mod]->GetBinContent(cellX,cellZ)>0)
    return 0 ;

  //New we are in good channel, 
  //calculate distance to the closest group of bad channels
  const Int_t cut1=10;
  const Int_t cut2=15;
//For 3-mod configuration
//  if((mod==3 && cellX<=cut1) || (mod==1 && cellX>=65-cut1) || cellZ<=cut1 || cellZ>=57-cut1)
//For 1+3 configuraion
  if( cellX<=cut1 ||  cellX>=65-cut1 || cellZ<=cut1 || cellZ>=57-cut1)
    return 1;
//  //and from large dead area
//  if(fPHOSBadMap[mod]->Integral(cellX-cut1,cellX+cut1,cellZ-cut1,cellZ+cut1)>0.2*cut1*cut1)
//    return 1 ;
//Full configuration
//    if((mod==3 && cellX<=cut2) || (mod==1 && cellX>=65-cut2) || cellZ<=cut2 || cellZ>=57-cut2)
//1+3 configuration
  if( cellX<=cut2 || cellX>=65-cut2 || cellZ<=cut2 || cellZ>=57-cut2)
    return 2;
//  if(fPHOSBadMap[mod]->Integral(cellX-cut2,cellX+cut2,cellZ-cut2,cellZ+cut2)>0.8*cut2*cut2)
//    return 2 ;
  //Very good channel
  return 3 ;

}
//______________________________________________________________________________^M
void  AliAnalysisTaskTaggedPhotons::InitGeometry(){
  //Rotation matrixes are set with Tender
  
  if(fPHOSgeom) return ;
  else
    fPHOSgeom = AliPHOSGeometry::GetInstance() ;
  
}
//_____________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::FillHistogram(const char * key,Double_t x)const{
  //FillHistogram
  TH1I * tmpI = dynamic_cast<TH1I*>(fOutputContainer->FindObject(key)) ;
  if(tmpI){
    tmpI->Fill(x) ;
    return ;
  }
  TH1F * tmpF = dynamic_cast<TH1F*>(fOutputContainer->FindObject(key)) ;
  if(tmpF){
    tmpF->Fill(x) ;
    return ;
  }
  TH1D * tmpD = dynamic_cast<TH1D*>(fOutputContainer->FindObject(key)) ;
  if(tmpD){
    tmpD->Fill(x) ;
    return ;
  }
  AliInfo(Form("can not find histogram <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::FillHistogram(const char * key,Double_t x,Double_t y)const{
  //FillHistogram
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH1F")){
    ((TH1F*)tmp)->Fill(x,y) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y) ;
    return ;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s",key,tmp->IsA()->GetName())) ;
}

//_____________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
  //Fills 1D histograms with key
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y,z) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F")){
    ((TH3F*)tmp)->Fill(x,y,z) ;
    return ;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::FillPIDHistograms(const char * name, const AliAODPWG4Particle * p) const{

  FillHistogram(Form("%s_All_cent%d",name,fCentBin),p->Pt()) ;
  if(p->GetDispBit())
    FillHistogram(Form("%s_Disp_cent%d",name,fCentBin),p->Pt()) ;
  if(p->GetChargedBit())
    FillHistogram(Form("%s_CPV_cent%d",name,fCentBin),p->Pt()) ;
  if(p->GetDispBit() && p->GetChargedBit()) 
    FillHistogram(Form("%s_Both_cent%d",name,fCentBin),p->Pt()) ;
  
}
//_____________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::FillPIDHistograms(const char * name, const AliAODPWG4Particle * p,Double_t x) const{

  FillHistogram(Form("%s_All_cent%d",name,fCentBin),x,p->Pt()) ;
  if(p->GetDispBit())
    FillHistogram(Form("%s_Disp_cent%d",name,fCentBin),x,p->Pt()) ;
  if(p->GetChargedBit())
    FillHistogram(Form("%s_CPV_cent%d",name,fCentBin),x,p->Pt()) ;
  if(p->GetDispBit() && p->GetChargedBit()) 
    FillHistogram(Form("%s_Both_cent%d",name,fCentBin),x,p->Pt()) ;
  
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskTaggedPhotons::TestLambda(Double_t pt,Double_t l1,Double_t l2){
  
  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;
  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
              0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
              0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return (R2<2.5*2.5) ;
  
}
//_________________________________________________________________________________
Int_t AliAnalysisTaskTaggedPhotons::EvalIsolation(TLorentzVector * ph){

   // Check if this particle is isolated. 
   //We use several cone radii and epsilons.
   //As well we look at charged particles and EMCAL ones

   const Double_t coneR1=0.3 ;
   const Double_t coneR2=0.4 ;
   const Double_t coneR3=0.5 ;

   const Double_t epsilon1=0.1 ;
   const Double_t epsilon2=0.05 ;

   if(!ph) return 0 ;

   //Sum of energies in cones, tracks and clusters in EMCAL
   Double_t eCone1 = 0;
   Double_t eCone2 = 0;
   Double_t eCone3 = 0;

   Double_t  phiTrig = ph->Phi();
   Double_t  etaTrig = ph->Eta();

   Int_t n=fTrackEvent->GetEntriesFast() ;
   for(Int_t itr=0; itr<n; itr++){
     AliAODPWG4Particle * track = (AliAODPWG4Particle*)fTrackEvent->At(itr) ;
         
     Double_t deleta = etaTrig - track->Eta() ;
     Double_t delphi = phiTrig - track->Phi() ;      
     while(delphi<-TMath::Pi()) delphi+=TMath::TwoPi() ;
     while(delphi>TMath::Pi()) delphi-=TMath::TwoPi() ;
   
     Double_t dr    = TMath::Sqrt(deleta * deleta + delphi * delphi);

     if(dr<coneR3){
       eCone3+=track->Pt() ;
       if(dr<coneR2){
	 eCone2+=track->Pt() ;
         if(dr<coneR1){
	   eCone1+=track->Pt() ;
	 }
	}
      }	
    }   
         
    //Fill QA histgams
    Double_t ptTrig=ph->Pt() ;
    FillHistogram(Form("QA_Cone1_Tracks_cent%d",fCentBin),ptTrig,eCone1) ;
    FillHistogram(Form("QA_Cone2_Tracks_cent%d",fCentBin),ptTrig,eCone2) ;
    FillHistogram(Form("QA_Cone3_Tracks_cent%d",fCentBin),ptTrig,eCone3) ;
    
    
    //Fill Bits
    Int_t iCone1E1 = (epsilon1*ptTrig > eCone1) ;
    Int_t iCone2E1 = (epsilon1*ptTrig > eCone2) ;
    Int_t iCone3E1 = (epsilon1*ptTrig > eCone3) ;
    
    Int_t iCone1E2 = (epsilon2*ptTrig > eCone1) ;
    Int_t iCone2E2 = (epsilon2*ptTrig > eCone2) ;
    Int_t iCone3E2 = (epsilon2*ptTrig > eCone3) ;
    
    
    Int_t isolation=   iCone1E1+  2*iCone2E1   +4*iCone3E1+
                    8*iCone1E2 +16*iCone2E2 +32*iCone3E2 ;
    return isolation ;		    
}
//_________________________________________________________________________________
Bool_t AliAnalysisTaskTaggedPhotons::TestPID(Int_t iPID, AliAODPWG4Particle* part){
//   //Checks PID of particle
  
  if(!part) return kFALSE ;
  if(iPID==0) return kTRUE ;
  if(iPID==1) return part->GetDispBit() ;
  if(iPID==2) return part->GetChargedBit() ;
  if(iPID==3) return part->GetDispBit() && part->GetChargedBit() ;
  return kFALSE ;
    
}
//_________________________________________________________________________________
Double_t AliAnalysisTaskTaggedPhotons::PrimaryParticleWeight(TParticle * /*particle*/){
  return 1;
  
}
//___________________________________________________________________________
Int_t AliAnalysisTaskTaggedPhotons::FindPrimary(AliVCluster*clu,  Bool_t&sure){
  //Finds primary and estimates if it unique one?
  //First check can it be photon/electron
  const Double_t emFraction=0.9; //part of energy of cluster to be assigned to EM particle
  Int_t n=clu->GetNLabels() ;
  FillHistogram(Form("LabelsNPrim_cent%d",fCentBin),clu->E(),n) ;
  if(n==1){
    sure=kTRUE;
    return clu->GetLabelAt(0) ;
  }
  //Number of clusters with one or more photons
  Bool_t hasGamma=kFALSE ;
  Double_t eMax=0. ;
  for(Int_t i=0;  i<n;  i++){
    TParticle*  p=  fStack->Particle(clu->GetLabelAt(i)) ;
    Int_t pdg = p->GetPdgCode() ;
    if(pdg==22){
      hasGamma=kTRUE ;
      if(p->Energy()>eMax){
	eMax=p->Energy();
      }
    }
  }
  if(hasGamma){
    FillHistogram(Form("LabelsGamma_cent%d",fCentBin),clu->E()) ;
    FillHistogram(Form("LabelsGammaE_cent%d",fCentBin),clu->E(),eMax/clu->E()) ;
  }  
  
  for(Int_t i=0;  i<n;  i++){
    TParticle*  p=  fStack->Particle(clu->GetLabelAt(i)) ;
    Int_t pdg = p->GetPdgCode() ;
    if(pdg==22  ||  pdg==11 || pdg == -11){
      if(p->Energy()>emFraction*clu->E()){
	sure=kTRUE ;
	return clu->GetLabelAt(i);
      }
    }
  }

  Double_t*  Ekin=  new  Double_t[n] ;
  for(Int_t i=0;  i<n;  i++){
    TParticle*  p=  fStack->Particle(clu->GetLabelAt(i)) ;
    Ekin[i]=p->P() ;  // estimate of kinetic energy
    if(p->GetPdgCode()==-2212  ||  p->GetPdgCode()==-2112){
      Ekin[i]+=1.8  ;  //due to annihilation
    }
  }
  Int_t iMax=0;
  Double_t eSubMax=0. ;
  eMax=0.;
  for(Int_t i=0;  i<n;  i++){
    if(Ekin[i]>eMax){
      eSubMax=eMax;
      eMax=Ekin[i];
      iMax=i;
    }
  }
  if(eSubMax>0.8*eMax)//not obvious primary
    sure=kFALSE;
  else
    sure=kTRUE;
  delete[]  Ekin;
  return  clu->GetLabelAt(iMax) ;
}


