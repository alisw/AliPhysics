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
#include "AliCaloPhoton.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "TGeoManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliPHOSGeometry.h"
#include "AliTriggerAnalysis.h"
#include "AliPHOSTriggerUtils.h"
#include "AliEMCALGeometry.h"
#include "AliAnalysisUtils.h"
#include "AliOADBContainer.h"
#include "AliAODMCHeader.h"

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
  fUtils(0x0),
  fPHOSTrigUtils(0x0),
  fZmax(0.),
  fZmin(0.),
  fPhimax(0.),
  fPhimin(0.),
  fMinBCDistance(0.),
  fCentrality(0.),
  fCentBin(0), 
  fIsMB(kTRUE),
  fIsMC(0),
  fIsFastMC(0)
{
  //Deafult constructor
  //no memory allocations
  for(Int_t i=0;i<10;i++)
    for(Int_t j=0;j<2;j++)
      fPHOSEvents[i][j]=0x0 ;    //Container for PHOS photons
  for(Int_t i=0;i<6;i++)
    fPHOSBadMap[i]=0x0;
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
  fUtils(0x0),
  fPHOSTrigUtils(0x0),
  fZmax(-60.),
  fZmin(60.),
  fPhimax(250.),
  fPhimin(320.),
  fMinBCDistance(0.),
  fCentrality(0.),
  fCentBin(0),
  fIsMB(0),
  fIsMC(0),
  fIsFastMC(0)
{
  // Constructor.

  // Output slots #0 write into a TH1 container
  DefineOutput(1,THashList::Class());
  // Set bad channel map (empty so far)
  for(Int_t i=0;i<1;i++)
    for(Int_t j=0;j<5;j++)
      fPHOSEvents[i][j]=0x0 ;    //Container for PHOS photons  
  for(Int_t i=0;i<6;i++)
    fPHOSBadMap[i]=0x0;
  
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
  fUtils(0x0),
  fPHOSTrigUtils(0x0),
  fZmax(-60.),
  fZmin(60.),
  fPhimax(250.),
  fPhimin(320.),
  fMinBCDistance(0.),  
  fCentrality(0.),
  fCentBin(0),
  fIsMB(kTRUE),
  fIsMC(0),
  fIsFastMC(0)
{
  // cpy ctor
  fZmax=ap.fZmax ;
  fZmin=ap.fZmin ;
  fPhimax=ap.fPhimax ;
  fPhimin=ap.fPhimin ;
  for(Int_t i=0;i<1;i++)
    for(Int_t j=0;j<5;j++)
      fPHOSEvents[i][j]=0x0 ;    //Container for PHOS photons
  for(Int_t i=0;i<6;i++)
    fPHOSBadMap[i]=0x0;

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
  fOutputContainer->Add(new TH2F("hTrackEtaPhi","Charged track eta vs phi distribution",200,-2.,2.,200,0.,TMath::TwoPi()));
  fOutputContainer->Add(new TH2F("hTrackEtaPt","Charged track eta vs pt distribution",200,-2.,2.,200,0.,50.));
  
  //centrality
  fOutputContainer->Add(new TH1F("hCentrality","Ccentrality",100,0.,100.));
  fOutputContainer->Add(new TH2F("hPHOSCentrality","PHOS vs centrality",100,0.,100.,100,0.,100.)); 
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
  
  fOutputContainer->Add(new TH2F("hTofM1","TOF in mod1",200,-1.e-6,1.e-6,200,0.,20.)) ;
  fOutputContainer->Add(new TH2F("hTofM2","TOF in mod2",200,-1.e-6,1.e-6,200,0.,20.)) ;
  fOutputContainer->Add(new TH2F("hTofM3","TOF in mod3",200,-1.e-6,1.e-6,200,0.,20.)) ;

  fOutputContainer->Add(new TH3F("hDispE","Dispersion vs E",120,0.,60.,100,0.,10.,100,0.,10.)) ;
  
  char cPID[4][5] ;
  snprintf(cPID[0],5,"All") ;
  snprintf(cPID[1],5,"Disp");
  snprintf(cPID[2],5,"CPV") ;
  snprintf(cPID[3],5,"Both"); 
  
 
  const Int_t nPt=650 ;
  const Double_t ptMax=65. ;
  const Int_t nM=400 ;
  const Double_t mMax=1. ;

  //QA histograms
   for(Int_t iPID=0; iPID<4; iPID++){
    fOutputContainer->Add(new TH1F(Form("hPhotM1_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M1",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hPhotM2_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M2",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hPhotM3_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M3",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hPhotM1A2_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M1",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hPhotM2A2_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M2",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hPhotM3A2_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M3",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hPhotM1A3_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M1",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hPhotM2A3_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M2",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH1F(Form("hPhotM3A3_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M3",nPt,0.,ptMax)) ;
   }

   //Trigger QA
   fOutputContainer->Add(new TH2F("hInvM_Re_mod1","Two-photon inv. mass (both in m1)",nM,0.,mMax,nPt,0.,ptMax)) ;
   fOutputContainer->Add(new TH2F("hInvM_Re_mod3","Two-photon inv. mass (both in m3)",nM,0.,mMax,nPt,0.,ptMax)) ;
   fOutputContainer->Add(new TH2F("hInvM_Mi_mod1","Two-photon inv. mass (both in m1)",nM,0.,mMax,nPt,0.,ptMax)) ;
   fOutputContainer->Add(new TH2F("hInvM_Mi_mod3","Two-photon inv. mass (both in m3)",nM,0.,mMax,nPt,0.,ptMax)) ;

  
  
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

    for(Int_t itag=0; itag<18; itag++){
      fOutputContainer->Add(new TH1F(Form("hPhot_TaggedMult%d_Area1_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of multiply tagged photons",nPt,0.,ptMax)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_TaggedMult%d_Area2_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of multiply tagged photons",nPt,0.,ptMax)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_TaggedMult%d_Area3_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of multiply tagged photons",nPt,0.,ptMax)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_Tagged%d_Area1_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_Tagged%d_Area2_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_Tagged%d_Area3_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;

      fOutputContainer->Add(new TH1F(Form("hPhot_TrueTagged%d_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of all tagged particles",nPt,0.,ptMax)) ;      
    }    
    for(Int_t kind=1; kind<33; kind*=2){
      fOutputContainer->Add(new TH1F(Form("hPhot_Isolation%d_%s_cent%d",kind,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_Isolation%d_Area1_%s_cent%d",kind,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_nTagged_Isolation%d_Area1_%s_cent%d",kind,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;
    }
  }
  for(Int_t kind=1; kind<33; kind*=2){
    fOutputContainer->Add(new TH1F(Form("hPi_Isolation%d_cent%d",kind,cen),"Spectrum of all reconstructed particles, no PID",nPt,0.,ptMax)) ;
  }

  

  //Invariant mass distributions for fake corrections
  
  for(Int_t iPID=0; iPID<4; iPID++){


    fOutputContainer->Add(new TH2F(Form("hInvM_Re_Emin1_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Re_Emin2_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Re_Emin3_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Mi_Emin1_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Mi_Emin2_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hInvM_Mi_Emin3_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Re_Emin1_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Re_Emin2_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Re_Emin3_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Re_Emin1_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Re_Emin2_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Re_Emin3_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;

    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Mi_Emin1_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Mi_Emin2_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Mi_Emin3_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Mi_Emin1_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Mi_Emin2_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Mi_Emin3_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;

    
    
  }  

 
  fOutputContainer->Add(new TH2F(Form("QA_Cone1_Tracks_cent%d",cen),"Cone energy",50,0.,50.,200,0.,50.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Cone2_Tracks_cent%d",cen),"Cone energy",50,0.,50.,200,0.,50.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Cone3_Tracks_cent%d",cen),"Cone energy",50,0.,50.,200,0.,50.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_PCone1_Tracks_cent%d",cen),"Cone energy",50,0.,50.,200,0.,50.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_PCone2_Tracks_cent%d",cen),"Cone energy",50,0.,50.,200,0.,50.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_PCone3_Tracks_cent%d",cen),"Cone energy",50,0.,50.,200,0.,50.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Pi0Cone1_Tracks_cent%d",cen),"Cone energy",50,0.,50.,200,0.,50.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Pi0Cone2_Tracks_cent%d",cen),"Cone energy",50,0.,50.,200,0.,50.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Pi0Cone3_Tracks_cent%d",cen),"Cone energy",50,0.,50.,200,0.,50.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Pi0PCone1_Tracks_cent%d",cen),"Cone energy",50,0.,50.,200,0.,50.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Pi0PCone2_Tracks_cent%d",cen),"Cone energy",50,0.,50.,200,0.,50.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Pi0PCone3_Tracks_cent%d",cen),"Cone energy",50,0.,50.,200,0.,50.)) ;
  }//centrality
  
  //MC  
  char partName[15][10] ;
  snprintf(partName[0],10,"gamma") ;
  snprintf(partName[1],10,"pi0");
  snprintf(partName[2],10,"eta") ;
  snprintf(partName[3],10,"omega"); 
  snprintf(partName[4],10,"K0s"); 
  snprintf(partName[5],10,"Kpm"); 
  snprintf(partName[6],10,"pipm"); 
  snprintf(partName[7],10,"n"); 
  snprintf(partName[8],10,"nbar"); 
  snprintf(partName[9],10,"p"); 
  snprintf(partName[10],10,"pbar"); 
  
  
  if(fIsMC){
      
      fOutputContainer->Add(new TH1F("hMCConversionRadius","Clusters without label",600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCRecPi0Vtx","Secondary pi0s",100,0.,10.,600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCRecEtaVtx","Secondary etas",100,0.,10.,600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCRecOmegaVtx","Secondary etas",100,0.,10.,600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCRecEtaprVtx","Secondary etas",100,0.,10.,600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCRecK0sVtx","Secondary K0s",100,0.,10.,600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCRecK0lVtx","Secondary K0l",100,0.,10.,600,0.,600.)) ;
      fOutputContainer->Add(new TH2F("hMCGammaPi0MisConvR","Converted photons",400,0.,40.,600,0.,600.)) ;
      
      for(Int_t mod=1; mod<5; mod++){
        fOutputContainer->Add(new TH1F(Form("hMCMinBiasPhot%d",mod),"MinBias photons",500,0.,50.)) ;
        fOutputContainer->Add(new TH1F(Form("hMCMinBiasPhotMap%d",mod),"MinBias photons in trigger area",500,0.,50.)) ;
        fOutputContainer->Add(new TH1F(Form("hMCTrigPhot%d",mod),"Triggered photons",500,0.,50.)) ;
      }
 
  for(Int_t ipart=0; ipart<11; ipart++){  
    fOutputContainer->Add(new TH2F(Form("hMC%s_ptrap",partName[ipart]),"Spectrum of primary photons",100,0.,10.,200,-1.,1.)) ;
    fOutputContainer->Add(new TH2F(Form("hMC%s_ptphi",partName[ipart]),"Spectrum of primary photons",100,0.,10.,100,0.,TMath::TwoPi())) ;
    fOutputContainer->Add(new TH2F(Form("hMC_%s_vertex",partName[ipart]),"vertex",nPt,0.,ptMax,150,0.,600.)) ;
    for(Int_t cen=0; cen<nCenBin; cen++){
       fOutputContainer->Add(new TH1F(Form("hMC_all_%s_cent%d",partName[ipart],cen),"Spectum (full rapifity)",nPt,0.,ptMax)) ;
       fOutputContainer->Add(new TH1F(Form("hMC_unitEta_%s_cent%d",partName[ipart],cen),"Spectum, |y|<0.15",nPt,0.,ptMax)) ;
       fOutputContainer->Add(new TH1F(Form("hMC_rap_%s_cent%d",partName[ipart],cen),"Rapidity",100,-5.,5.)) ;
       fOutputContainer->Add(new TH1F(Form("hMC_phi_%s_cent%d",partName[ipart],cen),"Azimuthal angle",100,0.,TMath::TwoPi())) ;
    }
  }
  for(Int_t cen=0; cen<nCenBin; cen++){
       fOutputContainer->Add(new TH1F(Form("hMC_unitEta_GammaPi0_cent%d",cen),"Spectum, |y|<0.15",nPt,0.,ptMax)) ;
  }
  
  for(Int_t cen=0; cen<nCenBin; cen++){
  
    fOutputContainer->Add(new TH2F(Form("LabelsNPrim_cent%d",cen),"vertex",nPt,0.,ptMax,20,0.,20.)) ;
    fOutputContainer->Add(new TH1F(Form("LabelsGamma_cent%d",cen),"vertex",nPt,0.,ptMax)) ;
    fOutputContainer->Add(new TH2F(Form("LabelsGammaE_cent%d",cen),"vertex",nPt,0.,ptMax,100,0.,2.)) ;
    
    
       //Sort registered particles spectra according MC information
       for(Int_t iPID=0; iPID<4; iPID++){
         fOutputContainer->Add(new TH1F(Form("hMCRecPhoton_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPhoton_Area1_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPhoton_Area2_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPhoton_Area3_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecE_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPbar_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecNbar_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecP_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPipm_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecKpm_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecN_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecCharg_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecNeutral_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecK0s_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecNoPRim_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,0.,ptMax )) ;

	 //Decay photons	 
         fOutputContainer->Add(new TH1F(Form("hMCRecPhotPi0_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPhotEta_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPhotOmega_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPhotOther_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPhotNoPrim_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,0.,ptMax )) ;
       
   
    
         //MC tagging: reasons of partner loss etc.
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnStack_%s_cent%d",cPID[iPID],cen),"Decay photons with partner not in Stack", nPt,0.,ptMax )) ;
         for(Int_t iType=0; iType<9; iType++)
  	   fOutputContainer->Add(new TH1F(Form("hMCDecWithFoundPartnType%d_%s_cent%d",iType,cPID[iPID],cen),"Decay photon with found partner", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWithWrongMass_%s_cent%d",cPID[iPID],cen),"Decay photon with wrong mass", nPt,0.,ptMax )) ;       
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnAccept_%s_cent%d",cPID[iPID],cen),"Decay photon with parttner not in PHOS acc", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnAcceptFA1_%s_cent%d",cPID[iPID],cen),"Decay photons with partner missed due geometry Fid. area. 1", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnAcceptFA2_%s_cent%d",cPID[iPID],cen),"Decay photons with partner missed due geometry Fid. area. 2", nPt,0.,ptMax)) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnAcceptFA3_%s_cent%d",cPID[iPID],cen),"Decay photons with partner missed due geometry Fid. area. 3", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnConv_%s_cent%d",cPID[iPID],cen),"Decay photons with partner missed due to conversion", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnEmin_%s_cent%d",cPID[iPID],cen),"Decay photons with partner missed due to low energy", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnOther_%s_cent%d",cPID[iPID],cen),"Decay photons with partner missed due unknown reason", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnAll_%s_cent%d",cPID[iPID],cen),"Decay photons with partner missed due to any reason", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnNPhot_%s_cent%d",cPID[iPID],cen),"pi0 decay photon with non-photon partner", nPt,0.,ptMax )) ;

         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnCutEmin_%s_cent%d",cPID[iPID],cen),"Decay photons with rec. partner but failed Emin cut", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnCutNcell_%s_cent%d",cPID[iPID],cen),"Decay photons with rec. partner but failed Ncell cut", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnCutEcross_%s_cent%d",cPID[iPID],cen),"Decay photons with rec. partner but failed Ecross cut", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnCutM02_%s_cent%d",cPID[iPID],cen),"Decay photons with rec. partner but failed M02 cut", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnDefCuts_%s_cent%d",cPID[iPID],cen),"Decay photons with rec. partner but failed default cuts", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWRecPartn_%s_cent%d",cPID[iPID],cen),"Decay photons with rec partner", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWRecUniqPartn_%s_cent%d",cPID[iPID],cen),"Decay photons with rec partner", nPt,0.,ptMax )) ;

         fOutputContainer->Add(new TH1F(Form("hMCDecMerged_%s_cent%d",cPID[iPID],cen),"Decay photons with rec partner", nPt,0.,ptMax )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecUnfolded_%s_cent%d",cPID[iPID],cen),"Decay photons with rec partner", nPt,0.,ptMax )) ;
 
	 
	 fOutputContainer->Add(new TH2F(Form("hMC_InvM_Re_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
         fOutputContainer->Add(new TH2F(Form("hMC_InvM_Re_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,0.,ptMax)) ;
       }    
       fOutputContainer->Add(new TH2F(Form("hMCmass_cent%d",cen),"Mass with reconstructed decay partner",nM,0.,mMax,nPt,0.,ptMax )) ;			
    }    
  }

  //If we work with MC, need to set Sumw2 - we will use weights
  if(fIsMC){
      for(Int_t i=0; i<fOutputContainer->GetSize();i++){
        ((TH1*)fOutputContainer->At(i))->Sumw2() ; 
      }
  }
    
  
  
  for(Int_t i=0;i<10;i++)
    for(Int_t j=0;j<5;j++)
      fPHOSEvents[i][j]=0x0 ;    //Container for PHOS photons
  

  //Prepare PHOS trigger utils if necessary
  if(!fIsMB )
   fPHOSTrigUtils = new AliPHOSTriggerUtils("PHOSTrig") ; 
  
  
  PostData(1, fOutputContainer);


}

//______________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::UserExec(Option_t *) 
{
  //Select events
  //Select photons
  //Fill QA histograms
  //Fill Tagging histogsms
  
  const Double_t kEcrossCut=0.98 ;
  const Double_t kTOFMaxCut= 100.e-9 ;  
  const Double_t kTOFMinCut=-100.e-9 ;  
  
  // Event selection flags
  //  FillHistogram("hSelEvents",0) ;
    
  AliVEvent* event = (AliVEvent*)InputEvent();
  if(!event){
    AliDebug(1,"No event") ;
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents",1) ;

  //MC stack init
  fStack = (TClonesArray*)event->FindListObject(AliAODMCParticle::StdBranchName());
  
  //read geometry if not read yet
  if(fPHOSgeom==0){
    InitGeometry() ;
  }

  if(!fUtils) 
      fUtils = new AliAnalysisUtils();
  
  if((!fIsFastMC) && (!fIsMC)){

    Bool_t isMB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7)  ; 
    Bool_t isPHI7 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kPHI7);

    if((fIsMB && !isMB) || (!fIsMB && !isPHI7)){
      PostData(1, fOutputContainer);
      return;    
    }
  }
  
  //If we work with PHOS trigger, prepapre TriggerUtils
  if(!fIsMB)
   fPHOSTrigUtils->SetEvent(event) ;
  
  FillHistogram("hSelEvents",2) ;
  
  // Checks if we have a primary vertex
  // Get primary vertices form AOD

  Double_t vtx5[3];
  vtx5[0] = event->GetPrimaryVertex()->GetX();
  vtx5[1] = event->GetPrimaryVertex()->GetY();
  vtx5[2] = event->GetPrimaryVertex()->GetZ();

  FillHistogram("hNvertexTracks",event->GetPrimaryVertex()->GetNContributors());
  FillHistogram("hZvertex"      ,vtx5[2]);
  if(fIsFastMC){ //vertex from header
    AliAODMCHeader *cHeaderAOD = dynamic_cast<AliAODMCHeader*>(event->FindListObject(AliAODMCHeader::StdBranchName()));
    if(!cHeaderAOD){
      PostData(1, fOutputContainer);
      return ;      
    }
    cHeaderAOD->GetVertex(vtx5);
  }
  if (TMath::Abs(vtx5[2]) > 10. ){
    PostData(1, fOutputContainer);
    return ;
  }
  FillHistogram("hSelEvents",3) ;
  //Vtx class z-bin
  Int_t zvtx = TMath::Min(9,Int_t((vtx5[2]+10.)/2.)) ; 

  
  if(!fIsFastMC){
//    if (event->IsPileupFromSPD()){
//      PostData(1, fOutputContainer);
//      return ;
//    }
  
    if(!fUtils->IsVertexSelected2013pA(event)){
      PostData(1, fOutputContainer);
      return ;
    }
  
    FillHistogram("hSelEvents",4) ;
  
    if(fUtils->IsPileUpEvent(event)){
      PostData(1, fOutputContainer);
      return ;
    }
    FillHistogram("hSelEvents",5) ;
  }
  
  //centrality
  if(!fIsFastMC){  
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
  }
  else{
    fCentrality=1.;
  }
  fCentBin = (Int_t)(fCentrality/20.) ; 

  FillHistogram("hSelEvents",6) ;
  
  
  //Calculate charged multiplicity
  Int_t trackMult = 0;
  if(fTrackEvent)
    fTrackEvent->Clear() ;
  else
    fTrackEvent = new TClonesArray("AliCaloPhoton",event->GetNumberOfTracks()) ;

  for (Int_t i=0;i<event->GetNumberOfTracks();++i) {
    AliAODTrack *track = (AliAODTrack*)event->GetTrack(i) ;
    if(!track->IsHybridGlobalConstrainedGlobal())
      continue ;
    if(TMath::Abs(track->Eta())< 0.9){
      if(trackMult>=fTrackEvent->GetSize())
	fTrackEvent->Expand(2*trackMult) ;
      new ((*fTrackEvent)[trackMult]) AliCaloPhoton(track->Px(),track->Py(),track->Pz(),track->P());
      trackMult++;
     FillHistogram("hTrackEtaPhi",track->Eta(),track->Phi()) ;
     FillHistogram("hTrackEtaPt",track->Eta(),track->Pt()) ;
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
    fPHOSEvent   = new TClonesArray("AliCaloPhoton",multClust);
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
  
    if(clu->GetMCEnergyFraction()>kEcrossCut) //Ecross cut, should be filled with Tender
     continue ;    
    
    if(clu->GetDistanceToBadChannel()<fMinBCDistance)
      continue ;

    Float_t pos[3] ;
    clu->GetPosition(pos) ;
    Int_t fidArea=GetFiducialArea(pos) ;
//    if(fidArea==0) //Bad cell
//      continue; 
    
    TVector3 global1(pos) ;
    Int_t relId[4] ;
    fPHOSgeom->GlobalPos2RelId(global1,relId) ;
    Int_t mod  = relId[0] ;
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;

    FillHistogram("hTOF",clu->E(),clu->GetTOF()) ;
    FillHistogram(Form("hTofM%d",mod),clu->GetTOF(),clu->E()) ;
    if((!fIsMC) && (clu->GetTOF() < kTOFMinCut || clu->GetTOF() > kTOFMaxCut))
      continue ;          
    
   FillHistogram("hDispE", clu->E(),clu->GetM02(),clu->GetM20()) ;
    
    TLorentzVector momentum ;
    clu->GetMomentum(momentum, vtx5);
    
    AliCaloPhoton *p = new ((*fPHOSEvent)[inList]) AliCaloPhoton(momentum.Px(),momentum.Py(),momentum.Pz(),clu->E() );
//    momentum*= clu->GetCoreEnergy()/momentum.E() ;
//    AliCaloPhoton *p = new ((*fPHOSEvent)[inList]) AliCaloPhoton(momentum.Px(),momentum.Py(),momentum.Pz(),clu->GetCoreEnergy() );
    inList++;
    
    p->SetModule(mod) ;
    Int_t isolation = EvalIsolation(&momentum,kTRUE) ;
    p->SetIsolationTag(isolation) ;
    
    p->SetDistToBad((Int_t)(1.+clu->GetDistanceToBadChannel()/2.2));
    p->SetBC(i) ; //reference to CaloCluster
    p->SetTagInfo(0); //No pi0 partners found so far
    p->SetTagged(kFALSE);   //Reconstructed pairs found
    
    p->SetFiducialArea(fidArea) ;

    //Mark photons fired trigger
    if(!fIsMB){   
      if(fIsMC){
        p->SetTrig(fPHOSTrigUtils->IsFiredTriggerMC(clu)) ;    
      }
      else
        p->SetTrig(fPHOSTrigUtils->IsFiredTrigger(clu)) ;    
    }
    
    if(fIsMB || ((!fIsMB) && p->IsTrig()) ){
      FillHistogram(Form("hCluNXZM%d",mod),cellX,cellZ,1.);
      FillHistogram(Form("hCluEXZM%d",mod),cellX,cellZ,clu->E());
      if(fidArea>1){
        FillHistogram(Form("hCluArea2M%d",mod),cellX,cellZ,1.);
        if(fidArea>2){
          FillHistogram(Form("hCluArea3M%d",mod),cellX,cellZ,1.);
        }
      }
    }
    
    if(fIsMC){    
       //Check trigger efficiency
       FillHistogram(Form("hMCMinBiasPhot%d",mod),clu->E()) ;
       //CheckTrigBadMap
       if((!fIsMB) && fPHOSTrigUtils->TestBadMap(mod, cellX,cellZ))
         FillHistogram(Form("hMCMinBiasPhotMap%d",mod),clu->E()) ;
       
       if(p->IsTrig())
          FillHistogram(Form("hMCTrigPhot%d",mod),clu->E()) ;
       //Look for MC particle entered PHOS
       FillHistogram(Form("LabelsNPrim_cent%d",fCentBin),clu->E(),float(clu->GetNLabels())) ;
       Int_t primLabel=clu->GetLabelAt(0) ; //FindPrimary(clu,sure) ;
       //Look what particle left vertex
       if(primLabel>-1){
         AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(primLabel) ;
         Int_t iparent=primLabel;
         AliAODMCParticle * parent = prim;
         Double_t r2=prim->Xv()*prim->Xv()+prim->Yv()*prim->Yv() ;
         while((r2 > rcut*rcut) && (iparent>-1)){
           iparent=parent->GetMother();
           parent=(AliAODMCParticle*)fStack->At(iparent);
           r2=parent->Xv()*parent->Xv()+parent->Yv()*parent->Yv() ;
         }
         p->SetPrimary(primLabel) ;
         p->SetPrimaryAtVertex(iparent) ;
	 p->SetWeight(PrimaryParticleWeight(parent)) ;
       }
       else{
         p->SetPrimary(-1); //Primary index    
         p->SetPrimaryAtVertex(-1) ;
	 p->SetWeight(1.) ;
       }
    }
    else{  
      p->SetPrimary(-1); //Primary index    
      p->SetPrimaryAtVertex(-1) ;
      p->SetWeight(1.) ;
    }
    //PID criteria
//  Cut on Core Lambdas    
    p->SetDispBit(clu->Chi2()<2.5*2.5) ;
//  Cut on FullLamdas
//    p->SetDispBit(clu->GetDispersion()<2.5*2.5) ;
    p->SetTOFBit(TestTOF(clu->GetTOF(),clu->E())) ;
    p->SetCPVBit(clu->GetEmcCpvDistance()>2.5) ;   
    
    if(fIsMB || (!fIsMB && p->IsTrig())){ 
      FillHistogram(Form("hPhotM%d_All",mod),p->Pt()) ;
      if(fidArea>1){
       FillHistogram(Form("hPhotM%dA2_All",mod),p->Pt()) ;
        if(fidArea>2){
          FillHistogram(Form("hPhotM%dA3_All",mod),p->Pt()) ;
        }
      }
      if(p->IsDispOK()){
        FillHistogram(Form("hPhotM%d_Disp",mod),p->Pt()) ;
        if(fidArea>1){
         FillHistogram(Form("hPhotM%dA2_Disp",mod),p->Pt()) ;
          if(fidArea>2){
            FillHistogram(Form("hPhotM%dA3_Disp",mod),p->Pt()) ;
          }
        }
        if(p->IsCPVOK()){
          FillHistogram(Form("hPhotM%d_Both",mod),p->Pt()) ;
          if(fidArea>1){
            FillHistogram(Form("hPhotM%dA2_Both",mod),p->Pt()) ;
            if(fidArea>2){
              FillHistogram(Form("hPhotM%dA3_Both",mod),p->Pt()) ;
            }
          }
        }
      } 
      if(p->IsCPVOK()){
        FillHistogram(Form("hPhotM%d_CPV",mod),p->Pt()) ;
        if(fidArea>1){
         FillHistogram(Form("hPhotM%dA2_CPV",mod),p->Pt()) ;
          if(fidArea>2){
            FillHistogram(Form("hPhotM%dA3_CPV",mod),p->Pt()) ;
          }
        }
      }
    }
  }
  FillHistogram("hPHOSCentrality",fCentrality,inList+0.5) ;
  
  
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
  if(!fIsMC)
    return ;
  const Double_t rcut=1. ; //cut on vertex to consider particle as "primary" 
  const Double_t phiMin=260.*TMath::Pi()/180. ;
  const Double_t phiMax=320.*TMath::Pi()/180. ;

  AliVEvent* event = (AliVEvent*)InputEvent();
  
  Int_t nPrim = fStack->GetEntriesFast() ;
  //Fill Primary particl yields
  
  for(Int_t i=0;i<nPrim;i++){
    AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(i) ;
    Double_t r2=prim->Xv()*prim->Xv()+prim->Yv()*prim->Yv() ;
    if(r2>rcut*rcut)
      continue ;

    Int_t pdg=prim->GetPdgCode() ;    
    char partName[30] ;
    if(pdg == 111)
      snprintf(partName,30,"pi0") ;
    else
      if(pdg == 221)
        snprintf(partName,30,"eta") ;
      else
        if(pdg == 22)
           snprintf(partName,30,"gamma") ;
	else
          if(pdg == 310)
             snprintf(partName,30,"K0s") ;
	  else
            if(abs(pdg) == 321)
              snprintf(partName,30,"Kpm") ;
  	    else
              if(abs(pdg) == 211)
                snprintf(partName,30,"pipm") ;
	      else  
                if(abs(pdg) == 2212)
                  snprintf(partName,30,"p") ;
 	        else  
                  if(abs(pdg) ==-2212)
                    snprintf(partName,30,"pbar") ;
   	          else  
                    if(abs(pdg) == 2112)
                      snprintf(partName,30,"n") ;
 	            else  
                      if(abs(pdg) ==-2112)
                        snprintf(partName,30,"nbar") ;
 	              else  
                        continue ;      	      

    //Primary particle
    Double_t phi=prim->Phi() ;
    while(phi<0.)phi+=TMath::TwoPi() ;
    while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
    Double_t pt=prim->Pt() ;

    //Total number of pi0 with creation radius <1 cm
    Double_t w = PrimaryParticleWeight(prim) ;  
    FillHistogram(Form("hMC_all_%s_cent%d",partName,fCentBin),pt,w) ;
    if(TMath::Abs(prim->Y())<0.13){
      FillHistogram(Form("hMC_phi_%s_cent%d",partName,fCentBin),phi,w) ;
      if(phi>phiMin && phi<phiMax){
        FillHistogram(Form("hMC_unitEta_%s_cent%d",partName,fCentBin),pt,w) ;
	if(prim->GetMother()>-1){
          AliAODMCParticle * primGM = (AliAODMCParticle*)fStack->At(prim->GetMother()) ;
	  if(primGM->GetPdgCode()==111)
            FillHistogram(Form("hMC_unitEta_GammaPi0_cent%d",fCentBin),pt,w) ;
	}
      }
    }

    FillHistogram(Form("hMC_rap_%s_cent%d",partName,fCentBin),prim->Y(),w) ;
    //Some additional QA
    if(pdg == 111){
       FillHistogram("hMCpi0_ptrap",pt,prim->Y(),w) ;   
       FillHistogram("hMCpi0_ptphi",pt,phi,w) ;   
    }
    if(pdg == 22){
       FillHistogram("hMCgamma_ptrap",pt,prim->Y(),w) ;   
       FillHistogram("hMCgamma_ptphi",pt,phi,w) ;   
    }
    
  }
  
 
  
  //Clussify reconstructed clusters
  //First - photons (from vertex) and contaminations
  //Second - photons from different sources
  //Third - photons from pi0s - missed for different reasons
  
  const Int_t n=fPHOSEvent->GetEntriesFast() ;
  for(Int_t i=0;i<n;i++){
    AliCaloPhoton *p = static_cast<AliCaloPhoton*>(fPHOSEvent->At(i));
    Int_t label=p->GetPrimary() ;
    if(label<0){ //No label!
      FillHistogram("hMCRecNoLabel",p->Pt(),p->GetWeight());
      continue ;
    }     

    
    AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(p->GetPrimary()) ;
    //Look what particle left virtex
    Int_t iparent=p->GetPrimary();
    AliAODMCParticle * parent = prim;
    while(parent->Xv()*parent->Xv()+parent->Yv()*parent->Yv() > rcut*rcut){
	iparent=parent->GetMother();
	if(iparent<0)
	  break ;
	parent = (AliAODMCParticle*)fStack->At(iparent) ;	
      }
      Int_t parentPDG=parent->GetPdgCode() ;    
      Int_t iFidArea=p->GetFiducialArea(); 
      switch(parentPDG){
	case 22: //electron/positron conversion
	  FillPIDHistograms("hMCRecPhoton",p);  //Reconstructed with photon from conversion primary
          if(iFidArea>0){
	    FillPIDHistograms("hMCRecPhoton_Area1",p);
            if(iFidArea>1){
	      FillPIDHistograms("hMCRecPhoton_Area2",p);
              if(iFidArea>1){
	        FillPIDHistograms("hMCRecPhoton_Area3",p);
	      }
	    }
	  }
	  break ;
	case  11:
	case -11: //electron/positron conversion
	  FillPIDHistograms("hMCRecE",p);  //Reconstructed with photon from conversion primary
	  break ;
	case -2212:
	  FillPIDHistograms("hMCRecPbar",p);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
	case -2112: //antineutron & antiproton conversion
	  FillPIDHistograms("hMCRecNbar",p);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
	case  211:
	case -211:
	  FillPIDHistograms("hMCRecPipm",p);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
	case 2212:
	  FillPIDHistograms("hMCRecP",p);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
	case  321:
	case -321:
	  FillPIDHistograms("hMCRecKpm",p);  //Reconstructed with photon from conversion primary
	  break ;
	case 310:
	  FillPIDHistograms("hMCRecK0s",p);  //Reconstructed with photon from conversion primary
	  break ;
	case 2112: //antineutron & antiproton conversion
	  FillPIDHistograms("hMCRecN",p);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
	case -1: //direct photon or no primary
	  FillPIDHistograms("hMCRecNoPRim",p);
	  break ;	  
	default:  
	  if(parent->Charge()!=0)
	    FillPIDHistograms("hMCRecCharg",p);  //Reconstructed with photon from antibaryon annihilation
	  else 
	    FillPIDHistograms("hMCRecNeutral",p);  //Reconstructed with photon from antibaryon annihilation
      }  
    
    
      //Now classify decay photon
      if(parentPDG==22){
	Int_t iGrandParent=parent->GetMother();
	if(iGrandParent<0 || iGrandParent>=fStack->GetEntriesFast()){
	  FillPIDHistograms("hMCRecPhotNoPrim",p);
          continue ;	  
	}
	AliAODMCParticle * grandParent = (AliAODMCParticle*)fStack->At(iGrandParent) ;	
        Int_t grandParentPDG=grandParent->GetPdgCode() ;     
        switch(grandParentPDG){
	case 111: //pi0
	  FillPIDHistograms("hMCRecPhotPi0",p);
	  break ;  		
	case 221: //eta decay
	  FillPIDHistograms("hMCRecPhotEta",p);
	  break ;  
	case 223: //omega meson decay
	  FillPIDHistograms("hMCRecPhotOmega",p);
	  break ;
	default:
	  FillPIDHistograms("hMCRecPhotOther",p);
	}
	//--------consider pi0 decays--------------------
	if(grandParentPDG==111){
	  //First find which daughter is our cluster
          //iparent - index of curent photon	  
	  Int_t ipartner=grandParent->GetDaughter(0) ;
	  if(ipartner==iparent){//look for other
  	    if(grandParent->GetNDaughters()>1){
	      ipartner=grandParent->GetDaughter(1);  
	    }
	    else{
	      ipartner=-1 ;
	    }
	  }
 	  //There is no partner in stack
	  if(ipartner==-1){
            FillPIDHistograms("hMCDecWMisPartnStack",p) ;
	  }
          else{
	    AliAODMCParticle * partner = (AliAODMCParticle *)fStack->At(ipartner);
	    //Check if partner is registered and made correct mass
	    //If not - trace the reason
	    AliCaloPhoton *pp = 0x0 ;
	  
	    for(Int_t ii=0;ii<n;ii++){
	      if(i==ii) continue; 
	      AliCaloPhoton * tmp = static_cast<AliCaloPhoton*>(fPHOSEvent->At(ii));
	      Int_t ipartnPrim = tmp->GetPrimary() ;
	      while(ipartnPrim>-1){
                if(ipartnPrim==ipartner){
		  break ;
		}
	        ipartnPrim = ((AliAODMCParticle *)fStack->At(ipartnPrim))->GetMother();
	      }
	      if(ipartnPrim==ipartner){
	        pp=tmp ;
	        break ;
	      }
	    }
	    if(pp){
	      //Partner reconstructed, but did not pass cuts
                FillPIDHistograms("hMCDecWRecUniqPartn",p) ;	
	    }
 	    //Partner not found. Check if it is not dominant contributor?
	    if(!pp){
  	      for(Int_t ii=0;(ii<n) && (!pp);ii++){
	        if(i==ii) continue; 
	        AliCaloPhoton * tmp = static_cast<AliCaloPhoton*>(fPHOSEvent->At(ii));
	        Int_t iCaloCluster=tmp->GetBC();//index of AODCaloCluster
                AliVCluster * clu = event->GetCaloCluster(iCaloCluster);
		Int_t nCluPrimaries = clu->GetNLabels() ;
		for(Int_t iAODLabel=0; (iAODLabel<nCluPrimaries) && (!pp); iAODLabel++){
		  Int_t ipartnPrim = clu->GetLabelAt(iAODLabel) ;
	          while(ipartnPrim>-1){
                    if(ipartnPrim==ipartner){
		      break ;
		    }
	            ipartnPrim = ((AliAODMCParticle *)fStack->At(ipartnPrim))->GetMother();
	          }
	          if(ipartnPrim==ipartner){
	            pp=tmp ;
	            break ;
	          }
		}
	      }
	    }
	    //If partner still not found, check if it is in same cluster
	    if(!pp){
	      Int_t iCaloCluster=p->GetBC();//index of AODCaloCluster
              AliVCluster * clu = event->GetCaloCluster(iCaloCluster);
	      Int_t nCluPrimaries = clu->GetNLabels() ;
	      for(Int_t iAODLabel=0; iAODLabel<nCluPrimaries ; iAODLabel++){
		Int_t ipartnPrim = clu->GetLabelAt(iAODLabel) ;
	        while(ipartnPrim>-1){
                  if(ipartnPrim==ipartner){
		      break ;
		  }
	          ipartnPrim = ((AliAODMCParticle *)fStack->At(ipartnPrim))->GetMother();
	        }
	        if(ipartnPrim==ipartner){ //yes, this cluster contains both primary
                  if(clu->GetNExMax()<2){ //was not unfolded
	            FillPIDHistograms("hMCDecMerged",p) ;
		  }
		  else{
	            FillPIDHistograms("hMCDecUnfolded",p) ;
		  }
                }
	      }
	    }	    

	    if(pp){
	      //Partner reconstructed, but did not pass cuts
                FillPIDHistograms("hMCDecWRecPartn",p) ;	
    	        Double_t invMass=(*p+ *pp).M() ;
	        FillHistogram(Form("hMCmass_cent%d",fCentBin),invMass,p->Pt(),p->GetWeight()) ;
		Double_t nSigma=InPi0Band(invMass,p->Pt()) ;
		// analog to Tag
                for(Int_t eminType=0; eminType<3; eminType++){
                  if(pp->E()>0.1*(eminType+1)){
  	            for(Int_t isigma=0; isigma<3; isigma++){
  	              if(nSigma<1.+isigma){
			 Int_t iType=3*eminType+isigma ;
	                 FillPIDHistograms(Form("hMCDecWithFoundPartnType%d",iType),p) ;
		      }
		    }
		  }
	        }
	        if(nSigma>3.){
	          FillPIDHistograms("hMCDecWithWrongMass",p) ;
	        }
	    }
	    else{//Partner not reconstructed
	      if(partner->GetPdgCode()==22){
		Bool_t isPartnerLost=kFALSE; //If partner is lost for some reason
		
		//Check if partner miss PHOS
  	        Int_t modulenum ;
		Double_t ztmp=0.,xtmp=0. ;
		Double_t vtx[3]={partner->Xv(),partner->Yv(),partner->Zv()} ;
		Bool_t impact=fPHOSgeom->ImpactOnEmc(vtx,partner->Theta(),partner->Phi(),modulenum,ztmp,xtmp) ;
		
		if(impact){//still check bad map
                  Int_t relId[4] ;
                  fPHOSgeom->RelPosToRelId(modulenum,xtmp,ztmp,relId) ;
                  if ( !IsGoodChannel(modulenum,relId[2],relId[3]) ) {
                     impact=kFALSE ;		    
		  }  
		}
 
		if(!impact){ //this photon cannot hit PHOS		  
		  FillPIDHistograms("hMCDecWMisPartnAccept",p) ;  //Spectrum of tagged with missed partner
		  if(iFidArea>0){
		    FillPIDHistograms("hMCDecWMisPartnAcceptFA1",p) ;  //Spectrum of tagged with missed partner
		    if(iFidArea>1){
		      FillPIDHistograms("hMCDecWMisPartnAcceptFA2",p) ;  //Spectrum of tagged with missed partner
		      if(iFidArea>2){
			FillPIDHistograms("hMCDecWMisPartnAcceptFA3",p) ;  //Spectrum of tagged with missed partner
		      }
		    }
		  }
		  isPartnerLost=kTRUE;
		}
		
                if(!isPartnerLost){
		  //this photon is converted before it is registered
		  if(partner->GetNDaughters()>0){
		    AliAODMCParticle* tmpP=(AliAODMCParticle*)fStack->At(partner->GetDaughter(0));
		    if(tmpP->Xv()*tmpP->Xv()+tmpP->Yv()*tmpP->Yv()<450.*450.){  
		      FillPIDHistograms("hMCDecWMisPartnConv",p) ;  //Spectrum of tagged with missed partner
		      isPartnerLost=kTRUE;
		    }
		  }
 		}
		if(!isPartnerLost && 
		   partner->E()<0.3){ //energy is not enough to be registered by PHOS
		  FillPIDHistograms("hMCDecWMisPartnEmin",p) ;  //Spectrum of tagged with missed partner
		  isPartnerLost=kTRUE;
		}
		if(!isPartnerLost){ //Reason not found!!!!!                  		  
		  FillPIDHistograms("hMCDecWMisPartnOther",p);
                  Int_t multClust = event->GetNumberOfCaloClusters();
                  for (Int_t iclu=0; (iclu<multClust) && (!isPartnerLost); iclu++) {
                    AliVCluster * clu = event->GetCaloCluster(iclu);
                    if(!clu->IsPHOS())
                      continue ; 
		    if(clu->E()==p->E()) //same cluster as current
		      continue ;
	      	    Int_t nCluPrimaries = clu->GetNLabels() ;
		    for(Int_t iAODLabel=0; (iAODLabel<nCluPrimaries) && (!isPartnerLost); iAODLabel++){
		      Int_t ipartnPrim = clu->GetLabelAt(iAODLabel) ;
	              while(ipartnPrim>-1){
                        if(ipartnPrim==ipartner)
		          break ;
	                ipartnPrim = ((AliAODMCParticle *)fStack->At(ipartnPrim))->GetMother();
	              }
	              if(ipartnPrim==ipartner){
		        isPartnerLost=kTRUE;
	                break ;
	              }
		    }
		  }
		  if(isPartnerLost){//Did not pass default cuts                 		  
		    FillPIDHistograms("hMCDecWMisPartnDefCuts",p);
		  }		  
		}
		else{//Sum of all missed partners
		  FillPIDHistograms("hMCDecWMisPartnAll",p);
		}
	      }//Partner - photon
	      else{//partner not photon
		FillPIDHistograms("hMCDecWMisPartnNPhot",p);                
	      }
	      
	    }//Partner not reconstructed
	  }//Partner in stack
	}//photon from pi0 decay
      }//photon
    } //PHOS clusters    
   
}
    
//________________________________________________
void AliAnalysisTaskTaggedPhotons::FillTaggingHistos(){
  //Fill all necessary histograms


 //First fill invariant mass distrubutions and mark tagged photons
  //Invariant Mass analysis
  const Int_t n=fPHOSEvent->GetEntriesFast() ;
  for(Int_t i=0;i<n-1;i++){
    AliCaloPhoton *p1 = static_cast<AliCaloPhoton*>(fPHOSEvent->At(i));
    for(Int_t j = i+1 ; j < n ; j++) {
      AliCaloPhoton * p2 = static_cast<AliCaloPhoton*>(fPHOSEvent->At(j));
      
      //At least one photon should be trigger in PHOS triggered events
      if(!fIsMB && !p1->IsTrig() &&  !p2->IsTrig() ) 
        continue ;
      
      Double_t invMass = (*p1 + *p2).M();   

      Double_t nsigma1 = InPi0Band(invMass,p1->Pt()) ; //in band with n sigmas
      Double_t nsigma2 = InPi0Band(invMass,p2->Pt()) ; //in band with n sigmas

      if((p1->E()>0.1) && (p2->E()>0.1)){
        FillPIDHistograms("hInvM_Re_Emin1",p1,p2,invMass,kTRUE) ;
        if((p1->E()>0.2) && (p2->E()>0.2)){
          FillPIDHistograms("hInvM_Re_Emin2",p1,p2,invMass,kTRUE) ;
          if((p1->E()>0.3) && (p2->E()>0.3)){
            FillPIDHistograms("hInvM_Re_Emin3",p1,p2,invMass,kTRUE) ; 
	    
	    if(p1->Module()==1 && p2->Module()==1){
	       Double_t w=TMath::Sqrt(p1->GetWeight()*p2->GetWeight()) ;
               FillHistogram("hInvM_Re_mod1",invMass,(*p1 + *p2).Pt(),w) ; 
	    }
	    if(p1->Module()==3 && p2->Module()==3){
	       Double_t w=TMath::Sqrt(p1->GetWeight()*p2->GetWeight()) ;
               FillHistogram("hInvM_Re_mod3",invMass,(*p1 + *p2).Pt(),w) ; 
	    }
	    //Fill izolated pi0s
	    if(nsigma1<2 || nsigma2<2){ //2 sigma band
	      TLorentzVector pi0=*p1+*p2 ;
              Int_t isolation=EvalIsolation(&pi0,0) ;
 	      for(Int_t kind=1; kind<33; kind*=2){
                 if((isolation&kind)){
                   FillHistogram(Form("hPi_Isolation%d_cent%d",kind,fCentBin),pi0.Pt()) ;
	         }
              }
	    }
	    
	  }
	}
      }
      
      
      if(p2->E()>0.1){
        FillPIDHistograms("hSingleInvM_Re_Emin1",p1,invMass) ;
        if(p2->E()>0.2){
          FillPIDHistograms("hSingleInvM_Re_Emin2",p1,invMass) ;
          if(p2->E()>0.3){
            FillPIDHistograms("hSingleInvM_Re_Emin3",p1,invMass) ;
	  }
	}
      }
	
      if(p1->E()>0.1){
        FillPIDHistograms("hSingleInvM_Re_Emin1",p2,invMass) ;
        if(p1->E()>0.2){
          FillPIDHistograms("hSingleInvM_Re_Emin2",p2,invMass) ;
          if(p1->E()>0.3){
            FillPIDHistograms("hSingleInvM_Re_Emin3",p2,invMass) ;
	  }
	}
      }
      if(TestPID(3, p2)){
        if(p2->E()>0.1){
          FillPIDHistograms("hSingleInvM_Re_Emin1_Strict",p1,invMass) ;
          if(p2->E()>0.2){
            FillPIDHistograms("hSingleInvM_Re_Emin2_Strict",p1,invMass) ;
            if(p2->E()>0.3){
              FillPIDHistograms("hSingleInvM_Re_Emin3_Strict",p1,invMass) ;
	    }
	  }
	}
      }
      if(TestPID(3, p1)){
        if(p1->E()>0.1){
          FillPIDHistograms("hSingleInvM_Re_Emin1_Strict",p2,invMass) ;
          if(p1->E()>0.1){
            FillPIDHistograms("hSingleInvM_Re_Emin2_Strict",p2,invMass) ;
            if(p1->E()>0.3){
              FillPIDHistograms("hSingleInvM_Re_Emin3_Strict",p2,invMass) ;
	    }
	  }
	}
      }
      if(IsSameParent(p1,p2)==111){
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
      //InPi0Band(mass, Ptphoton, type Emin cut
      //Set Bits to 1 if tagged
      Int_t tag1=0 ;
      for(Int_t eminType=0; eminType<3; eminType++){
        if(p2->E()>0.1*(eminType+1)){
	  //Set corresponding bit
	  for(Int_t isigma=0; isigma<3; isigma++){
  	    if(nsigma1<1+isigma){
   	      tag1|= (1 << (3*eminType+isigma)) ;
	      if(p2->IsPIDOK(3))
   	        tag1|= (1 << (3*eminType+isigma+9)) ;	  
	    }
	  }
	}
      }
      Int_t oldTag1=p1->GetTagInfo() ;
      for(Int_t ibit=0; ibit<18; ibit++){
        if(((oldTag1 & (1<<ibit))!=0) && //Already tagged 
           ((tag1 & (1<<ibit))!=0)){//Multiple tagging
         Int_t iFidArea = p1->GetFiducialArea(); 
         if(iFidArea>0){
           FillPIDHistograms(Form("hPhot_TaggedMult%d_Area1",ibit),p1) ;
           if(iFidArea>1){
             FillPIDHistograms(Form("hPhot_TaggedMult%d_Area2",ibit),p1) ;
             if(iFidArea>2){
               FillPIDHistograms(Form("hPhot_TaggedMult%d_Area3",ibit),p1) ;
             }
	   }
	 }
	}
      }
      tag1=tag1|oldTag1 ;
      p1->SetTagInfo(tag1) ;
      Int_t tag2=0 ;
      for(Int_t eminType=0; eminType<3; eminType++){
        if(p1->E()>0.1*(eminType+1)){
	  //Set corresponding bit
	  for(Int_t isigma=0; isigma<3; isigma++){
  	    if(nsigma2<1+isigma){
   	      tag2|= (1 << (3*eminType+isigma)) ;
  	      if(p1->IsPIDOK(3))
   	        tag2|= (1 << (3*eminType+isigma+9)) ;	  
	    }
	  }
	}
      }
      Int_t oldTag2=p2->GetTagInfo() ;
      for(Int_t ibit=0; ibit<18; ibit++){
        if(((oldTag2 & (1<<ibit))!=0) && //Already tagged 
           ((tag2 & (1<<ibit))!=0)){//Multiple tagging
         Int_t iFidArea = p2->GetFiducialArea(); 
         if(iFidArea>0){
           FillPIDHistograms(Form("hPhot_TaggedMult%d_Area1",ibit),p2) ;
           if(iFidArea>1){
             FillPIDHistograms(Form("hPhot_TaggedMult%d_Area2",ibit),p2) ;
             if(iFidArea>2){
               FillPIDHistograms(Form("hPhot_TaggedMult%d_Area3",ibit),p2) ;
             }
	   }
	 }
	}
      }
      tag2=tag2|oldTag2 ;
      p2->SetTagInfo(tag2) ;
            
      Bool_t trueTag=(IsSameParent(p1,p2)==111); //true tag
      p1->SetUnfolded(p1->IsntUnfolded()|trueTag) ;
      p2->SetUnfolded(p2->IsntUnfolded()|trueTag) ;
      
    }
  }
  
  
  
  //Single particle histgams
  for(Int_t i=0;i<n;i++){
    AliCaloPhoton *p = static_cast<AliCaloPhoton*>(fPHOSEvent->At(i));
    
   //photon should be trigger in PHOS triggered events
   if(!fIsMB && !p->IsTrig() ) 
     continue ;

    Int_t isolation = p->GetIsolationTag();

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
      Int_t tag=p->GetTagInfo() ;
      for(Int_t ibit=0; ibit<18; ibit++){
        if((tag & (1<<ibit))!=0){ 
          FillPIDHistograms(Form("hPhot_Tagged%d_Area1",ibit),p) ;
	  if(p->IsntUnfolded()) //true tag
             FillPIDHistograms(Form("hPhot_TrueTagged%d",ibit),p) ;
	  
	}
      }

      if(iFidArea>1){
        FillPIDHistograms("hPhot_Area2",p) ;
        for(Int_t ibit=0; ibit<18; ibit++){
          if((tag & (1<<ibit))!=0){ 
            FillPIDHistograms(Form("hPhot_Tagged%d_Area2",ibit),p) ;
	  }
        }
	if(iFidArea>2){
          FillPIDHistograms("hPhot_Area3",p) ;
          for(Int_t ibit=0; ibit<18; ibit++){
            if((tag & (1<<ibit))!=0){ 
              FillPIDHistograms(Form("hPhot_Tagged%d_Area3",ibit),p) ;
	    }
	  }
	}
      }
    }
  } 
  
   //Fill Mixed InvMass distributions:
  TIter nextEv(fCurrentMixedList) ;
  for(Int_t i=0;i<n;i++){
    AliCaloPhoton *p1 = static_cast<AliCaloPhoton*>(fPHOSEvent->At(i));
    while(TClonesArray * event2 = static_cast<TClonesArray*>(nextEv())){
      Int_t nPhotons2 = event2->GetEntriesFast() ;
      for(Int_t j=0; j < nPhotons2 ; j++){
        AliCaloPhoton * p2 = static_cast<AliCaloPhoton*>(event2->At(j)) ;
        Double_t invMass = (*p1 + *p2).M();

        //At least one photon should be trigger in PHOS triggered events
        if(!fIsMB && !p1->IsTrig() &&  !p2->IsTrig() ) 
          continue ;	
	
        if((p1->E()>0.1) && (p2->E()>0.1)){
          FillPIDHistograms("hInvM_Mi_Emin1",p1,p2,invMass,kFALSE) ;
          if((p1->E())>0.2 && (p2->E()>0.2)){
            FillPIDHistograms("hInvM_Mi_Emin2",p1,p2,invMass,kFALSE) ;
            if((p1->E())>0.3 && (p2->E()>0.3)){
              FillPIDHistograms("hInvM_Mi_Emin3",p1,p2,invMass,kFALSE) ;

	      if(p1->Module()==1 && p2->Module()==1){
	        Double_t w=TMath::Sqrt(p1->GetWeight()*p2->GetWeight()) ;
                FillHistogram("hInvM_Mi_mod1",invMass,(*p1 + *p2).Pt(),w) ;
	      }
	      if(p1->Module()==3 && p2->Module()==3){
	        Double_t w=TMath::Sqrt(p1->GetWeight()*p2->GetWeight()) ;
                FillHistogram("hInvM_Mi_mod3",invMass,(*p1 + *p2).Pt(),w) ; 
	      }
	    
  	    }
 	  }
        }
	
	
	if(p2->E()>0.1){
          FillPIDHistograms("hSingleInvM_Mi_Emin1",p1,invMass) ;
	  if(p2->E()>0.2){
            FillPIDHistograms("hSingleInvM_Mi_Emin2",p1,invMass) ;
	    if(p2->E()>0.3){
              FillPIDHistograms("hSingleInvM_Mi_Emin3",p1,invMass) ;
	    }
	  }
	}
        if(TestPID(3, p2)){
  	  if(p2->E()>0.1){
            FillPIDHistograms("hSingleInvM_Mi_Emin1_Strict",p1,invMass) ;
    	    if(p2->E()>0.2){
              FillPIDHistograms("hSingleInvM_Mi_Emin2_Strict",p1,invMass) ;
  	      if(p2->E()>0.3){
                FillPIDHistograms("hSingleInvM_Mi_Emin3_Strict",p1,invMass) ;
	      }
	    }
	  }
	}
	
	if(p1->E()>0.1){
          FillPIDHistograms("hSingleInvM_Mi_Emin1",p2,invMass) ;
  	  if(p1->E()>0.2){
            FillPIDHistograms("hSingleInvM_Mi_Emin2",p2,invMass) ;
	    if(p1->E()>0.3){
              FillPIDHistograms("hSingleInvM_Mi_Emin3",p2,invMass) ;
	    }
	  }
	}
        if(TestPID(3, p1)){
	  if(p1->E()>0.1){
            FillPIDHistograms("hSingleInvM_Mi_Emin1_Strict",p2,invMass) ;
	    if(p1->E()>0.2){
              FillPIDHistograms("hSingleInvM_Mi_Emin2_Strict",p2,invMass) ;
	      if(p1->E()>0.3){
                FillPIDHistograms("hSingleInvM_Mi_Emin3_Strict",p2,invMass) ;
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
Double_t AliAnalysisTaskTaggedPhotons::InPi0Band(Double_t m, Double_t pt)const
{
  //Parameterization of the pi0 peak region
  //for LHC13bcdef
//  Double_t mpi0mean =  0.13447 - 1.41259e-03 * TMath::Exp(-pt/1.30044) ;  
  //Parameterization of data 30.08.2014
  Double_t mpi0mean =  0.135 ;  

  //Double_t mpi0sigma=TMath::Sqrt(5.22245e-03*5.22245e-03 +2.86851e-03*2.86851e-03/pt) + 9.09932e-05*pt ;
  //Parameterization of data 30.08.2014
  Double_t mpi0sigma=TMath::Sqrt(4.67491e-03*4.67491e-03 +3.42884e-03*3.42884e-03/pt) + 1.24324e-04*pt ;

  return TMath::Abs(m-mpi0mean)/mpi0sigma ;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskTaggedPhotons::IsSameParent(const AliCaloPhoton *p1, const AliCaloPhoton *p2)const{
  //Looks through parents and finds if there was commont pi0 among ancestors

  if(!fIsMC)
    return 0 ; //can not say anything

  Int_t prim1 = p1->GetPrimary();
  while(prim1!=-1){ 
    Int_t prim2 = p2->GetPrimary();
  
    while(prim2!=-1){       
      if(prim1==prim2){
	return ((AliAODMCParticle*)fStack->At(prim1))->GetPdgCode() ;
      }
      prim2=((AliAODMCParticle*)fStack->At(prim2))->GetMother() ;
    }
    prim1=((AliAODMCParticle*)fStack->At(prim1))->GetMother() ;
  }
  return 0 ;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskTaggedPhotons::GetFiducialArea(const Float_t * position)const{
  //calculates in which kind of fiducial area photon hit

  TVector3 global1(position) ;
  Int_t relId[4] ;
  fPHOSgeom->GlobalPos2RelId(global1,relId) ;
//  Int_t mod  = relId[0] ;
  Int_t cellX = relId[2];
  Int_t cellZ = relId[3] ;

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
//Full configuration
//    if((mod==3 && cellX<=cut2) || (mod==1 && cellX>=65-cut2) || cellZ<=cut2 || cellZ>=57-cut2)
//1+3 configuration
  if( cellX<=cut2 || cellX>=65-cut2 || cellZ<=cut2 || cellZ>=57-cut2)
    return 2;
  //Very good channel
  return 3 ;

}
//______________________________________________________________________________^M
void  AliAnalysisTaskTaggedPhotons::InitGeometry(){
  //Rotation matrixes are set with Tender
  
  if(fPHOSgeom) return ;
  
  
  fPHOSgeom = AliPHOSGeometry::GetInstance() ;
 
  if(!fPHOSgeom){ //Geometry not yet constructed with Tender
    fPHOSgeom = AliPHOSGeometry::GetInstance("IHEP","");

    AliOADBContainer geomContainer("phosGeo");
    geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
    TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(170000,"PHOSRotationMatrixes");
    for(Int_t mod=0; mod<5; mod++) {
      if(!matrixes->At(mod)) continue;
      fPHOSgeom->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;   
    }
  }
  
  //Read BadMap for MC simulations
  Int_t runNumber=196208 ; //LHC13bcdef
  AliOADBContainer badmapContainer(Form("phosBadMap"));
  badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root","phosBadMap");
  TObjArray *maps = (TObjArray*)badmapContainer.GetObject(runNumber,"phosBadMap");
  if(!maps){
      AliError("TaggedPhotons: Can not read Bad map\n") ;    
  }
  else{
    AliInfo(Form("TaggedPhotons: Setting PHOS bad map with name %s \n",maps->GetName())) ;
    for(Int_t mod=0; mod<5;mod++){
      if(fPHOSBadMap[mod]) 
        delete fPHOSBadMap[mod] ;
      TH2I * h = (TH2I*)maps->At(mod) ;      
      if(h)
        fPHOSBadMap[mod]=new TH2I(*h) ;
    }
  }        
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
void AliAnalysisTaskTaggedPhotons::FillPIDHistograms(const char * name, const AliCaloPhoton * p) const{

  FillHistogram(Form("%s_All_cent%d",name,fCentBin),p->Pt(),p->GetWeight()) ;
  if(p->IsDispOK())
    FillHistogram(Form("%s_Disp_cent%d",name,fCentBin),p->Pt(),p->GetWeight()) ;
  if(p->IsCPVOK())
    FillHistogram(Form("%s_CPV_cent%d",name,fCentBin),p->Pt(),p->GetWeight()) ;
  if(p->IsDispOK() && p->IsCPVOK()) 
    FillHistogram(Form("%s_Both_cent%d",name,fCentBin),p->Pt(),p->GetWeight()) ;
  
}
//_____________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::FillPIDHistograms(const char * name, const AliCaloPhoton * p,Double_t x) const{

  FillHistogram(Form("%s_All_cent%d",name,fCentBin),x,p->Pt(),p->GetWeight()) ;
  if(p->IsDispOK())
    FillHistogram(Form("%s_Disp_cent%d",name,fCentBin),x,p->Pt(),p->GetWeight()) ;
  if(p->IsCPVOK())
    FillHistogram(Form("%s_CPV_cent%d",name,fCentBin),x,p->Pt(),p->GetWeight()) ;
  if(p->IsDispOK() && p->IsCPVOK()) 
    FillHistogram(Form("%s_Both_cent%d",name,fCentBin),x,p->Pt(),p->GetWeight()) ;
  
}
//_____________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::FillPIDHistograms(const char * name, const AliCaloPhoton * p1,const AliCaloPhoton * p2,Double_t x, Bool_t /*isRe*/) const{

  Double_t ptPi = (*p1 + *p2).Pt() ;
  Double_t w=TMath::Sqrt(p1->GetWeight()*p2->GetWeight()) ;
//  if(isRe){
//    Int_t pdg=IsSameParent(p1,p2) ;
//    if(pdg!=0)
//      w=p1->GetWeight() ;
//  }
  FillHistogram(Form("%s_All_cent%d",name,fCentBin),x,ptPi,w) ;
  if(p1->IsDispOK() && p2->IsDispOK())
    FillHistogram(Form("%s_Disp_cent%d",name,fCentBin),x,ptPi,w) ;
  if(p1->IsCPVOK() && p2->IsCPVOK())
    FillHistogram(Form("%s_CPV_cent%d",name,fCentBin),x,ptPi,w) ;
  if(p1->IsDispOK() && p1->IsCPVOK() && p2->IsDispOK() && p2->IsCPVOK()) 
    FillHistogram(Form("%s_Both_cent%d",name,fCentBin),x,ptPi,w) ;
  
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
Int_t AliAnalysisTaskTaggedPhotons::EvalIsolation(TLorentzVector * ph, Bool_t isPhoton){

   // Check if this particle is isolated. 
   //We use several cone radii and epsilons.
   //As well we look at charged particles and EMCAL ones

   const Double_t coneR1=0.3 ;
   const Double_t coneR2=0.4 ;
   const Double_t coneR3=0.5 ;

   const Double_t cutEcone1=1.5+0.8 ;
   const Double_t cutEcone2=2.+0.8 ;
   const Double_t cutEcone3=2.5+0.8 ;

   if(!ph) return 0 ;

   //Sum of energies in cones, tracks and clusters in EMCAL
   Double_t eCone1 = 0;
   Double_t eCone2 = 0;
   Double_t eCone3 = 0;
   //Cross-check, energy in cone perpr to photon
   Double_t eP1Cone1 = 0;
   Double_t eP1Cone2 = 0;
   Double_t eP1Cone3 = 0;
   Double_t eP2Cone1 = 0;
   Double_t eP2Cone2 = 0;
   Double_t eP2Cone3 = 0;
   

   Double_t  phiTrig = ph->Phi();
   Double_t  etaTrig = ph->Eta();

   Int_t n=fTrackEvent->GetEntriesFast() ;
   for(Int_t itr=0; itr<n; itr++){
     AliCaloPhoton * track = (AliCaloPhoton*)fTrackEvent->At(itr) ;
         
     Double_t deleta = etaTrig - track->Eta() ;
     Double_t delphi = phiTrig - track->Phi() ;      
     while(delphi<-TMath::Pi()) delphi+=TMath::TwoPi() ;
     while(delphi>TMath::Pi()) delphi-=TMath::TwoPi() ;
   
     //Perp cones
     Double_t delphiP1 = phiTrig +TMath::PiOver2()- track->Phi() ;      
     while(delphiP1<-TMath::Pi()) delphiP1+=TMath::TwoPi() ;
     while(delphiP1>TMath::Pi()) delphiP1-=TMath::TwoPi() ;
     Double_t delphiP2 = phiTrig -TMath::PiOver2()- track->Phi() ;      
     while(delphiP2<-TMath::Pi()) delphiP2+=TMath::TwoPi() ;
     while(delphiP2>TMath::Pi()) delphiP2-=TMath::TwoPi() ;
     
     
     Double_t dr    = TMath::Sqrt(deleta * deleta + delphi * delphi);
     Double_t drP1    = TMath::Sqrt(deleta * deleta + delphiP1 * delphiP1);
     Double_t drP2    = TMath::Sqrt(deleta * deleta + delphiP2 * delphiP2);

     if(dr<coneR3){
       eCone3+=track->Pt() ;
       if(dr<coneR2){
	 eCone2+=track->Pt() ;
         if(dr<coneR1){
	   eCone1+=track->Pt() ;
	 }
	}
      }	
       
    if(drP1<coneR3){
       eP1Cone3+=track->Pt() ;
       if(drP1<coneR2){
	 eP1Cone2+=track->Pt() ;
         if(drP1<coneR1){
	   eP1Cone1+=track->Pt() ;
	 }
       }
    }	
       
    if(drP2<coneR3){
       eP2Cone3+=track->Pt() ;
       if(drP2<coneR2){
	 eP2Cone2+=track->Pt() ;
         if(drP2<coneR1){
	   eP2Cone1+=track->Pt() ;
	 }
	}
      }	
    }   
         
    //Fill QA histgams
    Double_t ptTrig=ph->Pt() ;
    if(isPhoton){
    FillHistogram(Form("QA_Cone1_Tracks_cent%d",fCentBin),ptTrig,eCone1) ;
    FillHistogram(Form("QA_Cone2_Tracks_cent%d",fCentBin),ptTrig,eCone2) ;
    FillHistogram(Form("QA_Cone3_Tracks_cent%d",fCentBin),ptTrig,eCone3) ;
    FillHistogram(Form("QA_PCone1_Tracks_cent%d",fCentBin),ptTrig,eP1Cone1) ;
    FillHistogram(Form("QA_PCone2_Tracks_cent%d",fCentBin),ptTrig,eP1Cone2) ;
    FillHistogram(Form("QA_PCone3_Tracks_cent%d",fCentBin),ptTrig,eP1Cone3) ;
    FillHistogram(Form("QA_PCone1_Tracks_cent%d",fCentBin),ptTrig,eP2Cone1) ;
    FillHistogram(Form("QA_PCone2_Tracks_cent%d",fCentBin),ptTrig,eP2Cone2) ;
    FillHistogram(Form("QA_PCone3_Tracks_cent%d",fCentBin),ptTrig,eP2Cone3) ;
    }
    else{
    FillHistogram(Form("QA_Pi0Cone1_Tracks_cent%d",fCentBin),ptTrig,eCone1) ;
    FillHistogram(Form("QA_Pi0Cone2_Tracks_cent%d",fCentBin),ptTrig,eCone2) ;
    FillHistogram(Form("QA_Pi0Cone3_Tracks_cent%d",fCentBin),ptTrig,eCone3) ;
    FillHistogram(Form("QA_Pi0PCone1_Tracks_cent%d",fCentBin),ptTrig,eP1Cone1) ;
    FillHistogram(Form("QA_Pi0PCone2_Tracks_cent%d",fCentBin),ptTrig,eP1Cone2) ;
    FillHistogram(Form("QA_Pi0PCone3_Tracks_cent%d",fCentBin),ptTrig,eP1Cone3) ;
    FillHistogram(Form("QA_Pi0PCone1_Tracks_cent%d",fCentBin),ptTrig,eP2Cone1) ;
    FillHistogram(Form("QA_Pi0PCone2_Tracks_cent%d",fCentBin),ptTrig,eP2Cone2) ;
    FillHistogram(Form("QA_Pi0PCone3_Tracks_cent%d",fCentBin),ptTrig,eP2Cone3) ;
      
    }
    
    //Fill Bits
    Int_t iCone1E1 = (cutEcone1 > eCone1) ;
    Int_t iCone2E1 = (cutEcone2 > eCone2) ;
    Int_t iCone3E1 = (cutEcone3 > eCone3) ;
    
    Int_t iCone1E2 = (1.5*cutEcone1 > eCone1) ;
    Int_t iCone2E2 = (1.5*cutEcone2 > eCone2) ;
    Int_t iCone3E2 = (1.5*cutEcone3 > eCone3) ;
    
    
    Int_t isolation=   iCone1E1+  2*iCone2E1   +4*iCone3E1+
                    8*iCone1E2 +16*iCone2E2 +32*iCone3E2 ;
    return isolation ;		    
}
//_________________________________________________________________________________
Bool_t AliAnalysisTaskTaggedPhotons::TestPID(Int_t iPID, AliCaloPhoton* part){
//   //Checks PID of particle
  
  if(!part) return kFALSE ;
  if(iPID==0) return kTRUE ;
  if(iPID==1) return part->IsDispOK() ;
  if(iPID==2) return part->IsCPVOK() ;
  if(iPID==3) return part->IsDispOK() && part->IsCPVOK() ;
  return kFALSE ;
    
}
//_________________________________________________________________________________
Double_t AliAnalysisTaskTaggedPhotons::PrimaryParticleWeight(AliAODMCParticle * particle){
  if(!fIsMC)  //This is real data
     return 1;
  //Classify parent at vertex
  //Introduce for eta and pi0 weights   
     
  Double_t r2=particle->Xv()*particle->Xv()+particle->Yv()*particle->Yv() ;
  Int_t mother = particle->GetMother() ;
  while(mother>-1){
    if(r2<1.)
      break ;
    particle = (AliAODMCParticle*) fStack->At(mother);
    mother = particle->GetMother() ;
    r2=particle->Xv()*particle->Xv()+particle->Yv()*particle->Yv() ;
  }
  
  //Particle within 1 cm from the virtex
  Int_t pdg = particle->GetPdgCode() ;
  Double_t x = particle->Pt() ;
  mother = particle->GetMother() ;
  while(TMath::Abs(pdg)<100){//gamma, electrons, muons 
    if(mother>-1){
      AliAODMCParticle * tmpP=(AliAODMCParticle*)fStack->At(mother) ;
      pdg=tmpP->GetPdgCode() ;
      x = tmpP->Pt() ;
      mother = tmpP->GetMother() ;
    }
    else{ //direct photon/electron....
      return 1.; 
    }
  } 
  if(pdg == 111){
  //Pi0
     if(x<1) return 1. ;
     else return fWeightParamPi0[0]*(TMath::Power(x,fWeightParamPi0[1])*
       (1.+fWeightParamPi0[2]*x+fWeightParamPi0[3]*x*x)/
       (1.+fWeightParamPi0[4]*x+fWeightParamPi0[5]*x*x) +fWeightParamPi0[6])  ;
  }
  if(pdg == 221){
  //Eta - same same shape, but yield 0.48 and Br(eta->2gamma)
     Double_t norm=0.48 * 0.3943;
     if(x<1) return norm ;
     else return norm*fWeightParamPi0[0]*(TMath::Power(x,fWeightParamPi0[1])*
       (1.+fWeightParamPi0[2]*x+fWeightParamPi0[3]*x*x)/
       (1.+fWeightParamPi0[4]*x+fWeightParamPi0[5]*x*x) +fWeightParamPi0[6]) ;
  }
  return 1. ;
}
//_________________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::SetPi0WeightParameters(TArrayD * ar){
 
  for(Int_t i=0; i<7; i++){ //Array range
    if(ar->GetSize()>i) fWeightParamPi0[i]=ar->At(i) ;
    else fWeightParamPi0[i]=0.;
  }
  //normalize to obtain smooth transition at 1 GeV
  Double_t x=1. ;
  fWeightParamPi0[0]=1./(TMath::Power(x,fWeightParamPi0[1])*
       (1.+fWeightParamPi0[2]*x+fWeightParamPi0[3]*x*x)/
       (1.+fWeightParamPi0[4]*x+fWeightParamPi0[5]*x*x) +fWeightParamPi0[6]) ;
  
  
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
    AliAODMCParticle*  p=  (AliAODMCParticle*)fStack->At(clu->GetLabelAt(i)) ;
    Int_t pdg = p->GetPdgCode() ;
    if(pdg==22){
      hasGamma=kTRUE ;
      if(p->E()>eMax){
	eMax=p->E();
      }
    }
  }
  if(hasGamma){
    FillHistogram(Form("LabelsGamma_cent%d",fCentBin),clu->E()) ;
    FillHistogram(Form("LabelsGammaE_cent%d",fCentBin),clu->E(),eMax/clu->E()) ;
  }  
  
  for(Int_t i=0;  i<n;  i++){
    AliAODMCParticle*  p= (AliAODMCParticle*) fStack->At(clu->GetLabelAt(i)) ;
    Int_t pdg = p->GetPdgCode() ;
    if(pdg==22  ||  pdg==11 || pdg == -11){
      if(p->E()>emFraction*clu->E()){
	sure=kTRUE ;
	return clu->GetLabelAt(i);
      }
    }
  }

  Double_t*  Ekin=  new  Double_t[n] ;
  for(Int_t i=0;  i<n;  i++){
    AliAODMCParticle*  p=(AliAODMCParticle*) fStack->At(clu->GetLabelAt(i)) ;
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
//___________________________________________________________________________
Bool_t AliAnalysisTaskTaggedPhotons::IsGoodChannel(Int_t mod, Int_t ix, Int_t iz)
{
  //Check if this channel belogs to the good ones
  //Used only for estimate of photon impact in MC simulations
  //For true clusters a bad map in Tender should be used
  
  if(mod>4 || mod<1){
    return kFALSE ;
  }
  if(!fPHOSBadMap[mod]){
     return kTRUE ;
  }
  if(fPHOSBadMap[mod]->GetBinContent(ix,iz)>0)
    return kFALSE ;
  else
    return kTRUE ;
}
