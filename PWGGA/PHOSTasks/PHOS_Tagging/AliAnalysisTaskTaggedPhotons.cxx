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
#include <TProfile.h>
#include <TH2.h>
#include <TH3.h>
#include <TMatrixD.h>
#include <TMatrixDEigen.h>
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
#include "AliMultSelection.h"
#include "AliMCAnalysisUtils.h"
#include "AliGenPythiaEventHeader.h"

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
  fMCGenerEventHeaderToAccept(""),
  fCentEstimator(1),
  fNCenBin(5),   
  fCentrality(0),
  fCentWeight(1.),
  fCentBin(0), 
  fRunNumber(0),
  fForseRun(kFALSE),
  fIsMB(kTRUE),
  fUseCaloFastTr(kFALSE),
  fIsMC(0),
  fIsFastMC(0),
  fIsMCWithPileup(0),
  fDefTof(true),
  fRP(0.),
  fJetPtHardFactor(2.5),
  fZmax(0.),
  fZmin(0.),
  fPhimax(0.),
  fPhimin(0.),
  fMinBCDistance(0.),
  fTimeCut(12.5e-9),
  fEminCut(0.1),
  fNonlinA(1.),
  fNonlinB(0.),
  fNonlinC(1.),
  fNonlinD(0.),
  fNonlinE(1.),
  fUnsmA(0.5),
  fUnsmB(1.8),
  fNPID(4),
  fMCType(kFullMC),
  fCutType(kDefCut),
  fPHOSTrigger(kPHOSAny),
  fTrackSelection(kLHC13x)
{
  //Deafult constructor
  //no memory allocations
  for(Int_t i=0;i<10;i++)
    for(Int_t j=0;j<5;j++)
      fPHOSEvents[i][j]=0x0 ;    //Container for PHOS photons
      
  for(Int_t i=0;i<6;i++)
    fPHOSBadMap[i]=0x0;
  //Centrality binning
  fCenBinEdges.Set(fNCenBin);
  for(int cen=1; cen<=fNCenBin; cen++)
    fCenBinEdges.AddAt(int(100.*cen/fNCenBin),cen-1) ; 

  fWeightParamPi0[0]=1.;  

  for(Int_t i=1; i<7; i++)
    fWeightParamPi0[i]=0.;  

  for(Int_t i=0; i<6; i++)
   fCentralityWeights[i]= 0x0; //!Weights to correct centrality non-flatness

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
  fMCGenerEventHeaderToAccept(""),
  fCentEstimator(1),
  fNCenBin(5),   
  fCentrality(0),
  fCentWeight(1.),
  fCentBin(0), 
  fRunNumber(0),
  fForseRun(kFALSE),
  fIsMB(kTRUE),
  fUseCaloFastTr(kFALSE),
  fIsMC(0),
  fIsFastMC(0),
  fIsMCWithPileup(0),
  fDefTof(true),
  fRP(0.),
  fJetPtHardFactor(2.5),
  fZmax(-60.),
  fZmin(60.),
  fPhimax(250.),
  fPhimin(320.),
  fMinBCDistance(0.),
  fTimeCut(12.5e-9),
  fEminCut(0.1),
  fNonlinA(1.),
  fNonlinB(0.),
  fNonlinC(1.),
  fNonlinD(0.),
  fNonlinE(1.),
  fUnsmA(0.5),
  fUnsmB(1.8),
  fNPID(4),
  fMCType(kFullMC),
  fCutType(kDefCut),
  fPHOSTrigger(kPHOSAny),
  fTrackSelection(kLHC13x)
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
  //Centrality binning
  fCenBinEdges.Set(fNCenBin);
  for(int cen=1; cen<=fNCenBin; cen++)
    fCenBinEdges.AddAt(int(100.*cen/fNCenBin),cen-1) ; 
  
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
  fMCGenerEventHeaderToAccept(""),
  fCentEstimator(1),
  fNCenBin(5),   
  fCentrality(0),
  fCentWeight(1.),
  fCentBin(0), 
  fRunNumber(0),
  fForseRun(kFALSE),
  fIsMB(kTRUE),
  fUseCaloFastTr(kFALSE),
  fIsMC(0),
  fIsFastMC(0),
  fIsMCWithPileup(0),
  fDefTof(true),
  fRP(0.),
  fJetPtHardFactor(2.5),
  fZmax(-60.),
  fZmin(60.),
  fPhimax(250.),
  fPhimin(320.),
  fMinBCDistance(0.),
  fTimeCut(12.5e-9),
  fEminCut(0.1),
  fNonlinA(1.),
  fNonlinB(0.),
  fNonlinC(1.),
  fNonlinD(0.),
  fNonlinE(1.),
  fUnsmA(0.5),
  fUnsmB(1.8),
  fNPID(4),
  fMCType(kFullMC),
  fCutType(kDefCut),
  fPHOSTrigger(kPHOSAny),
  fTrackSelection(kLHC13x)
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
  fCenBinEdges.Set(ap.fCenBinEdges.GetSize(),ap.fCenBinEdges.GetArray());

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
  fOutputContainer->Add(new TH1F("hSelEvents","Event selection", 15,0.,15.)) ;
  
  //vertex distribution
  fOutputContainer->Add(new TH1F("hNvertexTracks","N of primary tracks from the primary vertex",150,0.,150.));
  fOutputContainer->Add(new TH1F("hZvertex","Z vertex",200,-50.,+50.));
  fOutputContainer->Add(new TH2F("hTrackMult","Charged track multiplicity",100,0.,100.,250,0.,500.));
  fOutputContainer->Add(new TH2F("hTrackEtaPhi","Charged track eta vs phi distribution",200,-2.,2.,200,0.,TMath::TwoPi()));
  fOutputContainer->Add(new TH2F("hTrackEtaPt","Charged track eta vs pt distribution",200,-2.,2.,200,0.,50.));
  
  //centrality
  fOutputContainer->Add(new TH1F("hCentrality","Ccentrality",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentralityRaw","Centrality",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentralityV0A","Centrality V0A",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentralityV0C","Centrality V0C",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentralityV0M","Centrality V0M",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentralityCL1","Centrality CL1",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentralityZNA","Centrality ZNA",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentralityZNC","Centrality ZNC",105,0.,105.));
  fOutputContainer->Add(new TH2F("hCentralityCorr","Centrality",105,0.,105.,105,0.,105.));
 
  fOutputContainer->Add(new TH1F("hCentrality1V0A","Centrality V0A",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1V0C","Centrality V0C",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1V0M","Centrality V0M",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1CL1","Centrality CL1",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1ZNA","Centrality ZNA",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1ZNC","Centrality ZNC",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1OnlineV0M","Centrality OnlineV0M",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1OnlineV0A","Centrality OnlineV0A",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1OnlineV0C","Centrality OnlineV0C",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1ADM","Centrality ADM",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1ADA","Centrality ADA",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1ADC","Centrality ADC",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1SPDClusters","Centrality SPDClusters",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1SPDTracklets","Centrality SPDTracklets",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1RefMult05","Centrality RefMult05",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality1RefMult08","Centrality RefMult08",105,0.,105.));
  fOutputContainer->Add(new TH3F("hSpheriSpheroMult","SSM",22,-0.2,1.,22,-0.2,1.,100,0.,100.)) ;


  fOutputContainer->Add(new TH1F("hCentrality2V0A","Centrality V0A",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2V0C","Centrality V0C",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2V0M","Centrality V0M",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2CL1","Centrality CL1",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2ZNA","Centrality ZNA",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2ZNC","Centrality ZNC",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2OnlineV0M","Centrality OnlineV0M",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2OnlineV0A","Centrality OnlineV0A",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2OnlineV0C","Centrality OnlineV0C",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2ADM","Centrality ADM",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2ADA","Centrality ADA",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2ADC","Centrality ADC",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2SPDClusters","Centrality SPDClusters",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2SPDTracklets","Centrality SPDTracklets",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2RefMult05","Centrality RefMult05",105,0.,105.));
  fOutputContainer->Add(new TH1F("hCentrality2RefMult08","Centrality RefMult08",105,0.,105.));  

  fOutputContainer->Add(new TH2F("hPHOSCentrality","PHOS vs centrality",105,0.,105.,100,0.,100.)); 
  fOutputContainer->Add(new TH2F("hTOF","cluster TOF",200,0.,20.,300,-3.e-6,6.e-6));
  
  //Reaction Plane
  fOutputContainer->Add(new TH2F("hRPCentrality","Reaction plane",100,0.,100.,50,0.,TMath::Pi()));
  
  
  fOutputContainer->Add(new TH2F("hCluNXZM1","Clu (X,Z), M1"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM2","Clu (X,Z), M2"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM3","Clu (X,Z), M3"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM4","Clu (X,Z), M4"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM1","Clu E(X,Z), M1"  ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM2","Clu E(X,Z), M2"  ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM3","Clu E(X,Z), M3"  ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM4","Clu E(X,Z), M4"  ,64,0.5,64.5, 56,0.5,56.5));
  
  fOutputContainer->Add(new TH2F("hTofM1","TOF in mod1",200,-1.e-7,1.e-7,200,0.,20.)) ;
  fOutputContainer->Add(new TH2F("hTofM2","TOF in mod2",200,-1.e-7,1.e-7,200,0.,20.)) ;
  fOutputContainer->Add(new TH2F("hTofM3","TOF in mod3",200,-1.e-7,1.e-7,200,0.,20.)) ;
  fOutputContainer->Add(new TH2F("hTofM4","TOF in mod4",200,-1.e-7,1.e-7,200,0.,20.)) ;

  fOutputContainer->Add(new TH3F("hDispE","Dispersion vs E",130,0.,65.,100,0.,10.,100,0.,10.)) ;
  fOutputContainer->Add(new TH3F("hDispEneu","Dispersion vs E",130,0.,65.,100,0.,10.,100,0.,10.)) ;
  
  char cPID[8][15] ;
  snprintf(cPID[0],5,"All") ;
  snprintf(cPID[1],5,"Disp");
  snprintf(cPID[2],5,"CPV") ;
  snprintf(cPID[3],5,"Both"); 
  snprintf(cPID[4],6,"Disp2");
  snprintf(cPID[5],6,"Disp3");
  snprintf(cPID[6],6,"Both2"); 
  snprintf(cPID[7],6,"Both3"); 
  
  fNPID=4 ;  //Extend to 8 to look at PID cuts systematics
  
  Int_t nPt=73;
  Double_t ptBins[74]={0.,0.1,0.2,0.3,0.4, 0.5,0.6,0.7,0.8,0.9, 1.0,1.1,1.2,1.3,1.4, 1.5,1.6,1.7,1.8,1.9, 2.0,2.2,2.4,2.5,2.6, 
                       2.8,3.,3.2,3.4,3.6, 3.8,4.0,4.5,4.8,5., 5.5,5.6,6.0,6.4,6.5, 7.0,7.2,7.5,8.0,8.5, 9.0,9.5,10.,11.,12., 13.,14.,15.,16.,17., 18.,19.,20.,22.,24., 25.,26.,28.,30.,32., 35.,40.,45.,50.,55., 60.,65.,70.,80.};

                     
                     
  const Int_t nM=500 ;
  const Double_t mMax=1. ;

  //QA histograms
   for(Int_t iPID=0; iPID<fNPID; iPID++){
    fOutputContainer->Add(new TH1F(Form("hPhotM1_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M1",nPt,ptBins)) ;
    fOutputContainer->Add(new TH1F(Form("hPhotM2_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M2",nPt,ptBins)) ;
    fOutputContainer->Add(new TH1F(Form("hPhotM3_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M3",nPt,ptBins)) ;
    fOutputContainer->Add(new TH1F(Form("hPhotM4_%s",cPID[iPID]),"Spectrum of all reconstructed particles, M3",nPt,ptBins)) ;
   }

   //Module QA
   for(Int_t mod=1; mod<5; mod++){
     fhReMod[mod] = new TH2F(Form("hInvM_Re_mod%d",mod),"Two-photon inv. mass (both in m1)",nM,0.,mMax,nPt,ptBins) ;       //!
     fOutputContainer->Add(fhReMod[mod]) ;
     fhMiMod[mod] = new TH2F(Form("hInvM_Mi_mod%d",mod),"Two-photon inv. mass (both in m1)",nM,0.,mMax,nPt,ptBins);       //!
     fOutputContainer->Add(fhMiMod[mod]) ;
   }
  
  
  for(Int_t cen=0; cen<fNCenBin; cen++){

  for(Int_t iPID=0; iPID<fNPID; iPID++){
    
    //Inclusive spectra
    fOutputContainer->Add(new TH1F(Form("hPhot_%s_cent%d",cPID[iPID],cen),"Spectrum of all reconstructed particles",nPt,ptBins)) ;
    fOutputContainer->Add(new TH1F(Form("hPhotWcorr_%s_cent%d",cPID[iPID],cen),"Spectrum of all reconstructed particles",nPt,ptBins)) ;
    fOutputContainer->Add(new TH1F(Form("hPhot_Dist2_%s_cent%d",cPID[iPID],cen),"Spectrum of all reconstructed particles",nPt,ptBins)) ;
    fOutputContainer->Add(new TH1F(Form("hPhot_Dist3_%s_cent%d",cPID[iPID],cen),"Spectrum of all reconstructed particles",nPt,ptBins)) ;
  
    for(Int_t itag=0; itag<18; itag++){
      fOutputContainer->Add(new TH1F(Form("hPhot_TaggedMult%d_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of multiply tagged photons",nPt,ptBins)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_TaggedMult%d_Isolation2_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of multiply tagged isolated photons",nPt,ptBins)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_Tagged%d_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of tagged photons",nPt,ptBins)) ;
      fOutputContainer->Add(new TH1F(Form("hPhotWcorr_Tagged%d_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of tagged photons",nPt,ptBins)) ;      
      fOutputContainer->Add(new TH1F(Form("hPhot_Tagged%d_Isolation2_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of tagged and isolated photons",nPt,ptBins)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_TrueTagged%d_%s_cent%d",itag,cPID[iPID],cen),"Spectrum of all tagged particles",nPt,ptBins)) ;      
    }    
    for(Int_t kind=0; kind<20; kind++){
      fOutputContainer->Add(new TH1F(Form("hPhot_Isolation%d_%s_cent%d",kind,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,ptBins)) ;
      fOutputContainer->Add(new TH1F(Form("hPhot_Tagged_Isolation%d_%s_cent%d",kind,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,ptBins)) ;
    }
    for(Int_t kind=0; kind<20; kind++){
      fOutputContainer->Add(new TH1F(Form("hPhot_nTagged_Isolation%d_%s_cent%d",kind,cPID[iPID],cen),"Spectrum of all reconstructed particles, no PID",nPt,ptBins)) ;
    }
  }
  for(Int_t kind=0; kind<20; kind++){
    fhPiIsolation[kind][cen]=new TH1F(Form("hPi_Isolation%d_cent%d",kind,cen),"Spectrum of all reconstructed particles, no PID",nPt,ptBins) ;
    fOutputContainer->Add(fhPiIsolation[kind][cen]) ;
  }

  

  //Invariant mass distributions for fake corrections
  
  for(Int_t iPID=0; iPID<fNPID; iPID++){
    for(Int_t iEmin=0; iEmin<3; iEmin++){
       fhRe[iEmin][cen][iPID] = new TH2F(Form("hInvM_Re_Emin%d_%s_cent%d",iEmin+1,cPID[iPID],cen),
                                         "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ; 
       fOutputContainer->Add(fhRe[iEmin][cen][iPID] ) ;
       fhMi[iEmin][cen][iPID] = new TH2F(Form("hInvM_Mi_Emin%d_%s_cent%d",iEmin+1,cPID[iPID],cen),
                                         "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ;
       fOutputContainer->Add(fhMi[iEmin][cen][iPID]) ;

       fhReSingle[iEmin][cen][iPID]= new TH2F(Form("hSingleInvM_Re_Emin%d_%s_cent%d",iEmin+1,cPID[iPID],cen),
                                               "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ;
       fOutputContainer->Add(fhReSingle[iEmin][cen][iPID]) ;
       
       fhMiSingle[iEmin][cen][iPID]= new TH2F(Form("hSingleInvM_Mi_Emin%d_%s_cent%d",iEmin+1,cPID[iPID],cen),
                                               "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ;
       fOutputContainer->Add(fhMiSingle[iEmin][cen][iPID]) ;
       
       fhReSingleIso[iEmin][cen][iPID]= new TH2F(Form("hSingleInvM_Re_Emin%d_Iso_%s_cent%d",iEmin+1,cPID[iPID],cen),
                                               "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ;
       fOutputContainer->Add(fhReSingleIso[iEmin][cen][iPID]) ;

       fhReTruePi0[iEmin][cen][iPID] = new TH2F(Form("hInvM_ReTruePi0_Emin%d_%s_cent%d",iEmin+1,cPID[iPID],cen),
                                         "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ; 
       fOutputContainer->Add(fhReTruePi0[iEmin][cen][iPID] ) ;

       fhReTrueEta[iEmin][cen][iPID] = new TH2F(Form("hInvM_ReTrueEta_Emin%d_%s_cent%d",iEmin+1,cPID[iPID],cen),
                                         "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ; 
       fOutputContainer->Add(fhReTrueEta[iEmin][cen][iPID] ) ;
       
       fhMiSingleIso[iEmin][cen][iPID]= new TH2F(Form("hSingleInvM_Mi_Emin%d_Iso_%s_cent%d",iEmin+1,cPID[iPID],cen),
                                               "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ;
       fOutputContainer->Add(fhMiSingleIso[iEmin][cen][iPID]) ;
    } 
        
    
    
    
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Re_Emin1_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Re_Emin2_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Re_Emin3_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins)) ;

    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Mi_Emin1_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Mi_Emin2_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F(Form("hSingleInvM_Mi_Emin3_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins)) ;

  }  

 
  fOutputContainer->Add(new TH2F(Form("QA_Cone1_Tracks_cent%d",cen),"Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Cone2_Tracks_cent%d",cen),"Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Cone3_Tracks_cent%d",cen),"Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_TCone1_Tracks_cent%d",cen),"Tagged Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_TCone2_Tracks_cent%d",cen),"Tagged Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_TCone3_Tracks_cent%d",cen),"Tagged Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_PCone1_Tracks_cent%d",cen),"Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_PCone2_Tracks_cent%d",cen),"Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_PCone3_Tracks_cent%d",cen),"Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Pi0Cone1_Tracks_cent%d",cen),"Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Pi0Cone2_Tracks_cent%d",cen),"Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Pi0Cone3_Tracks_cent%d",cen),"Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Pi0PCone1_Tracks_cent%d",cen),"Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Pi0PCone2_Tracks_cent%d",cen),"Cone energy",nPt,ptBins,400,0.,10.)) ;
  fOutputContainer->Add(new TH2F(Form("QA_Pi0PCone3_Tracks_cent%d",cen),"Cone energy",nPt,ptBins,400,0.,10.)) ;
  
  fOutputContainer->Add(new TProfile(Form("QA_profCone1_Tracks_cent%d",cen),"Cone energy",nPt,ptBins)) ;
  fOutputContainer->Add(new TProfile(Form("QA_profCone2_Tracks_cent%d",cen),"Cone energy",nPt,ptBins)) ;
  fOutputContainer->Add(new TProfile(Form("QA_profCone3_Tracks_cent%d",cen),"Cone energy",nPt,ptBins)) ;
  fOutputContainer->Add(new TProfile(Form("QA_profTCone1_Tracks_cent%d",cen),"Cone energy",nPt,ptBins)) ;
  fOutputContainer->Add(new TProfile(Form("QA_profTCone2_Tracks_cent%d",cen),"Cone energy",nPt,ptBins)) ;
  fOutputContainer->Add(new TProfile(Form("QA_profTCone3_Tracks_cent%d",cen),"Cone energy",nPt,ptBins)) ;
  fOutputContainer->Add(new TProfile(Form("QA_profPCone1_Tracks_cent%d",cen),"Cone energy",nPt,ptBins)) ;
  fOutputContainer->Add(new TProfile(Form("QA_profPCone2_Tracks_cent%d",cen),"Cone energy",nPt,ptBins)) ;
  fOutputContainer->Add(new TProfile(Form("QA_profPCone3_Tracks_cent%d",cen),"Cone energy",nPt,ptBins)) ;
  
  fOutputContainer->Add(new TH3F(Form("QA_PCone12_Tracks_cent%d",cen),"Cone1 vs Cone2",65,0.,65.,100,0.,25.,100,0.,25.)) ;
  fOutputContainer->Add(new TH3F(Form("QA_PConeRP_Tracks_cent%d",cen),"Cone1 vs RP",65,0.,65.,100,0.,25.,6,0.,TMath::Pi())) ;
  
  }//centrality
 
  //pi0 partners energy (asimmetry)
  fhQAIsoEpartn = new TH2F("QA_Iso_Epartn","Partner energy distribuion",60,0.,60.,100,0.,1.) ;
  fOutputContainer->Add(fhQAIsoEpartn) ;
  fhQAAllEpartn = new TH2F("QA_All_Epartn","Partner energy distribuion",60,0.,60.,100,0.,1.) ;
  fOutputContainer->Add(fhQAAllEpartn) ;
  fhQAIsoEpartnBg = new TH2F("QA_Iso_Epartn_Bg","Partner energy distribuion, Bg region",60,0.,60.,100,0.,1.) ;
  fOutputContainer->Add(fhQAIsoEpartnBg) ;
  fhQAAllEpartnBg = new TH2F("QA_All_Epartn_Bg","Partner energy distribuion, Bg region",60,0.,60.,100,0.,1.) ;
  fOutputContainer->Add(fhQAAllEpartnBg) ;

  fhQAIsozpartn = new TH2F("QA_Iso_zpartn","Partner dZ distribuion",60,0.,60.,100,0.,100.) ;
  fOutputContainer->Add(fhQAIsozpartn) ;
  fhQAAllzpartn = new TH2F("QA_All_zpartn","Partner dZ distribuion",60,0.,60.,100,0.,100.) ;
  fOutputContainer->Add(fhQAAllzpartn) ;
  fhQAIsozpartnBg = new TH2F("QA_Iso_zpartn_Bg","Partner dZ distribuion, Bg region",60,0.,60.,100,0.,100.) ;
  fOutputContainer->Add(fhQAIsozpartnBg) ;
  fhQAAllzpartnBg = new TH2F("QA_All_zpartn_Bg","Partner dZ distribuion, Bg region",60,0.,60.,100,0.,100.) ;
  fOutputContainer->Add(fhQAAllzpartnBg) ;
 
  fhQAIsoxpartn = new TH2F("QA_Iso_xpartn","Partner dX distribuion",60,0.,60.,100,0.,100.) ; 
  fOutputContainer->Add(fhQAIsoxpartn) ;
  fhQAAllxpartn = new TH2F("QA_All_xpartn","Partner dX distribuion",60,0.,60.,100,0.,100.) ;
  fOutputContainer->Add(fhQAAllxpartn) ;
  fhQAIsoxpartnBg = new TH2F("QA_Iso_xpartn_Bg","Partner dX distribuion, Bg region",60,0.,60.,100,0.,100.) ;
  fOutputContainer->Add(fhQAIsoxpartnBg) ;
  fhQAAllxpartnBg = new TH2F("QA_All_xpartn_Bg","Partner dX distribuion, Bg region",60,0.,60.,100,0.,100.) ;
  fOutputContainer->Add(fhQAAllxpartnBg) ;
 
 
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
            
      for(Int_t mod=1; mod<5; mod++){
        fOutputContainer->Add(new TH1F(Form("hMCMinBiasPhot%d",mod),"MinBias photons",500,0.,50.)) ;
        fOutputContainer->Add(new TH1F(Form("hMCMinBiasPhotMap%d",mod),"MinBias photons in trigger area",500,0.,50.)) ;
        fOutputContainer->Add(new TH1F(Form("hMCTrigPhot%d",mod),"Triggered photons",500,0.,50.)) ;
      }
 
  for(Int_t ipart=0; ipart<11; ipart++){  
    fOutputContainer->Add(new TH2F(Form("hMC%s_ptrap",partName[ipart]),"Spectrum of primary photons",100,0.,10.,200,-1.,1.)) ;
    fOutputContainer->Add(new TH2F(Form("hMC%s_ptphi",partName[ipart]),"Spectrum of primary photons",100,0.,10.,100,0.,TMath::TwoPi())) ;
    fOutputContainer->Add(new TH2F(Form("hMC_%s_vertex",partName[ipart]),"vertex",nPt,ptBins,150,0.,600.)) ;
    for(Int_t cen=0; cen<fNCenBin; cen++){
       fOutputContainer->Add(new TH1F(Form("hMC_all_%s_cent%d",partName[ipart],cen),"Spectum (full rapifity)",nPt,ptBins)) ;
       fOutputContainer->Add(new TH1F(Form("hMC_unitEta_%s_cent%d",partName[ipart],cen),"Spectum, |y|<0.15",nPt,ptBins)) ;
       fOutputContainer->Add(new TH1F(Form("hMC_rap_%s_cent%d",partName[ipart],cen),"Rapidity",100,-1.,1.)) ;
       fOutputContainer->Add(new TH1F(Form("hMC_phi_%s_cent%d",partName[ipart],cen),"Azimuthal angle",100,0.,TMath::TwoPi())) ;
    }
  }
  for(Int_t cen=0; cen<fNCenBin; cen++){
    fOutputContainer->Add(new TH1F(Form("hMC_unitEta_GammaPi0_cent%d",cen),"Spectum, |y|<0.15",nPt,ptBins)) ;
  }
  for(Int_t iPID=0; iPID<fNPID; iPID++){
    for(Int_t iEmin=0; iEmin<3; iEmin++){
      fhReSingleTruePi0[iEmin][iPID] = new TH2F(Form("hSingleTruePi0InvM_Re_Emin%d_%s",iEmin+1,cPID[iPID]),
                                                       "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ;
      fOutputContainer->Add(fhReSingleTruePi0[iEmin][iPID]) ;
        
      fhReSingleTruePi0K0s[iEmin][iPID] = new TH2F(Form("hSingleTruePi0K0sInvM_Re_Emin%d_%s",
                                                               iEmin+1,cPID[iPID]),
                                            "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ;
      fOutputContainer->Add(fhReSingleTruePi0K0s[iEmin][iPID]) ;
      fhReSingleTruePi0Kpm[iEmin][iPID] = new TH2F(Form("hSingleTruePi0KpmInvM_Re_Emin%d_%s",iEmin+1,cPID[iPID]),
                                                       "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ;
      fOutputContainer->Add(fhReSingleTruePi0Kpm[iEmin][iPID]) ;
      fhReSingleTruePi0Nbar[iEmin][iPID] = new TH2F(Form("hSingleTruePi0InvM_Re_Emin%d_%s",iEmin+1,cPID[iPID]),
                                                       "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ;
      fOutputContainer->Add(fhReSingleTruePi0Nbar[iEmin][iPID]) ;
      
        
      fhReSingleTruePi0K0sC[iEmin][iPID] = new TH2F(Form("hSingleTruePi0K0sCorrInvM_Re_Emin%d_%s",
                                                               iEmin+1,cPID[iPID]),
                                            "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ;
      fOutputContainer->Add(fhReSingleTruePi0K0sC[iEmin][iPID]) ;
      fhReSingleTruePi0KpmC[iEmin][iPID] = new TH2F(Form("hSingleTruePi0KpmCorrInvM_Re_Emin%d_%s",iEmin+1,cPID[iPID]),
                                                       "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ;
      fOutputContainer->Add(fhReSingleTruePi0KpmC[iEmin][iPID]) ;
      fhReSingleTruePi0NbarC[iEmin][iPID] = new TH2F(Form("hSingleTruePi0NbarCorrInvM_Re_Emin%d_%s",iEmin+1,cPID[iPID]),
                                                       "Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins) ;
      fOutputContainer->Add(fhReSingleTruePi0NbarC[iEmin][iPID]) ;
      
    }
  }
  
  for(Int_t cen=0; cen<fNCenBin; cen++){
  
    fOutputContainer->Add(new TH2F(Form("LabelsNPrim_cent%d",cen),"vertex",nPt,ptBins,20,0.,20.)) ;
    fOutputContainer->Add(new TH1F(Form("LabelsGamma_cent%d",cen),"vertex",nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F(Form("LabelsGammaE_cent%d",cen),"vertex",nPt,ptBins,100,0.,2.)) ;
    
    
       //Sort registered particles spectra according MC information
       for(Int_t iPID=0; iPID<fNPID; iPID++){
         fOutputContainer->Add(new TH1F(Form("hMCRecPhoton_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPhotonOnly_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecE_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPbar_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecNbar_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecP_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPipm_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecKpm_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecN_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecCharg_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecNeutral_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecK0s_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecK0l_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecNoPRim_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. electrons", nPt,ptBins )) ;

	 //Decay photons	 
         fOutputContainer->Add(new TH1F(Form("hMCRecPhotPi0_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPhotEta_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPhotOmega_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPhotOther_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCRecPhotNoPrim_%s_cent%d",cPID[iPID],cen),"Spectrum of rec. photons", nPt,ptBins )) ;
       
   
    
         //MC tagging: reasons of partner loss etc.
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnStack_%s_cent%d",cPID[iPID],cen),"Decay photons with partner not in Stack", nPt,ptBins )) ;
         for(Int_t iType=0; iType<9; iType++)
  	   fOutputContainer->Add(new TH1F(Form("hMCDecWithFoundPartnType%d_%s_cent%d",iType,cPID[iPID],cen),"Decay photon with found partner", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWithWrongMass_%s_cent%d",cPID[iPID],cen),"Decay photon with wrong mass", nPt,ptBins )) ;       
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnAccept_%s_cent%d",cPID[iPID],cen),"Decay photon with parttner not in PHOS acc", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnConv_%s_cent%d",cPID[iPID],cen),"Decay photons with partner missed due to conversion", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnEmin_%s_cent%d",cPID[iPID],cen),"Decay photons with partner missed due to low energy", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnOther_%s_cent%d",cPID[iPID],cen),"Decay photons with partner missed due unknown reason", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnAll_%s_cent%d",cPID[iPID],cen),"Decay photons with partner missed due to any reason", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnNPhot_%s_cent%d",cPID[iPID],cen),"pi0 decay photon with non-photon partner", nPt,ptBins )) ;

         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnCutEmin_%s_cent%d",cPID[iPID],cen),"Decay photons with rec. partner but failed Emin cut", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnCutNcell_%s_cent%d",cPID[iPID],cen),"Decay photons with rec. partner but failed Ncell cut", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnCutEcross_%s_cent%d",cPID[iPID],cen),"Decay photons with rec. partner but failed Ecross cut", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnCutM02_%s_cent%d",cPID[iPID],cen),"Decay photons with rec. partner but failed M02 cut", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWMisPartnDefCuts_%s_cent%d",cPID[iPID],cen),"Decay photons with rec. partner but failed default cuts", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWRecPartn_%s_cent%d",cPID[iPID],cen),"Decay photons with rec partner", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecWRecUniqPartn_%s_cent%d",cPID[iPID],cen),"Decay photons with rec partner", nPt,ptBins )) ;

         fOutputContainer->Add(new TH1F(Form("hMCDecMerged_%s_cent%d",cPID[iPID],cen),"Decay photons with rec partner", nPt,ptBins )) ;
         fOutputContainer->Add(new TH1F(Form("hMCDecUnfolded_%s_cent%d",cPID[iPID],cen),"Decay photons with rec partner", nPt,ptBins )) ;
 
	 
	 fOutputContainer->Add(new TH2F(Form("hMC_InvM_Re_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins)) ;
         fOutputContainer->Add(new TH2F(Form("hMC_InvM_Re_Strict_%s_cent%d",cPID[iPID],cen),"Two-photon inv. mass vs first photon pt",nM,0.,mMax,nPt,ptBins)) ;
       }    
       fOutputContainer->Add(new TH2F(Form("hMCmass0_cent%d",cen),"Mass with reconstructed decay partner Emin=0.1",nM,0.,mMax,nPt,ptBins )) ;			
       fOutputContainer->Add(new TH2F(Form("hMCmass1_cent%d",cen),"Mass with reconstructed decay partner Emin=0.2",nM,0.,mMax,nPt,ptBins )) ;			
       fOutputContainer->Add(new TH2F(Form("hMCmass2_cent%d",cen),"Mass with reconstructed decay partner Emin=0.3",nM,0.,mMax,nPt,ptBins )) ;			
    }    
  }

  //If we work with MC or use trigger, need to set Sumw2 - we will use weights
//   if(fIsMC){
      for(Int_t i=0; i<fOutputContainer->GetSize();i++){
        ((TH1*)fOutputContainer->At(i))->Sumw2() ; 
      }
//   }
    
  
  
  for(Int_t i=0;i<10;i++)
    for(Int_t j=0;j<fNCenBin;j++)
      fPHOSEvents[i][j]=0x0 ;    //Container for PHOS photons
  

  //Prepare PHOS trigger utils if necessary
  if(!fIsMB ){
    fPHOSTrigUtils = new AliPHOSTriggerUtils("PHOSTrig") ; 
    if(fForseRun){
      printf("Forse run %d \n", fRunNumber) ;       
      fPHOSTrigUtils->ForseUsingRun(fRunNumber) ; 
    }
  } 
   
  if(fMCType==kJetJet){
    fMCAnalysisUtils = new AliMCAnalysisUtils() ; 
  }
  
  PostData(1, fOutputContainer);


}

//______________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::UserExec(Option_t *) 
{
  //Select events
  //Select photons
  //Fill QA histograms
  //Fill Tagging histogsms
        
  AliAODEvent* event = (AliAODEvent*)InputEvent();
  if(!event){
    AliDebug(1,"No event") ;
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents",0) ;
  if(!fForseRun)
    fRunNumber=event->GetRunNumber() ;

  //MC stack init
  fStack = (TClonesArray*)event->FindListObject(AliAODMCParticle::StdBranchName());
  
  //read geometry if not read yet
  if(fPHOSgeom==0){
    InitGeometry() ;
  }

  if(!fUtils) 
      fUtils = new AliAnalysisUtils();
  
  if((!fIsFastMC) && (!fIsMC)){
    TString trigClasses = ((AliAODHeader*)event->GetHeader())->GetFiredTriggerClasses();
    if(fUseCaloFastTr){
      if(fIsMB){  //Select only INT7 events
        if( !trigClasses.Contains("CINT7-B-NOPF-") ){
          PostData(1, fOutputContainer);
          return;    
        }
      }
      else{ //PHI7 
        if(!trigClasses.Contains("CPHI7-B-NOPF-CALOFAST")){ //only calofast to avoid double counting
          PostData(1, fOutputContainer);
          return;    
        }
      }  
    } else{
        
      Bool_t isMB = (fInputHandler->IsEventSelected() & AliVEvent::kINT7)  ; 
      Bool_t isPHI7 = (fInputHandler->IsEventSelected() & AliVEvent::kPHI7);
      
      switch(fPHOSTrigger){
      case kPHOSL0: isPHI7 = (trigClasses.Contains("CPHI7-B-NOPF-")) ;
                    break ;
      case kPHOSL1low : isPHI7 =( trigClasses.Contains("CPHI7PHL-B-NOPF-") &&
                                 !trigClasses.Contains("CPHI7-B-NOPF-") ) ;
                        break ;
      case kPHOSL1med: isPHI7 =( trigClasses.Contains("CPHI7PHM-B-NOPF-") &&
                                !trigClasses.Contains("CPHI7PHL-B-NOPF-")) ;
                       break ;
      case kPHOSL1high: isPHI7=( trigClasses.Contains("CPHI7PHH-B-NOPF-") &&
                                !trigClasses.Contains("CPHI7PHM-B-NOPF-")) ;
                        break ;
      case kPHOSAny: ; //do nothing
      default: ;
      }
        
      if((fIsMB && !isMB) || (!fIsMB && !isPHI7)){
        PostData(1, fOutputContainer);
        return;    
      }
    }
  }

  //If we work with PHOS trigger, prepapre TriggerUtils
  if(!fIsMB)
   fPHOSTrigUtils->SetEvent(event) ;
  
  FillHistogram("hSelEvents",1) ;
  
  // Checks if we have a primary vertex
  // Get primary vertices form AOD

  Double_t vtx5[3];
  vtx5[0] = event->GetPrimaryVertex()->GetX();
  vtx5[1] = event->GetPrimaryVertex()->GetY();
  vtx5[2] = event->GetPrimaryVertex()->GetZ();

  FillHistogram("hNvertexTracks",event->GetPrimaryVertex()->GetNContributors());
  if(fIsFastMC){ //vertex from header
    AliAODMCHeader *cHeaderAOD = dynamic_cast<AliAODMCHeader*>(event->FindListObject(AliAODMCHeader::StdBranchName()));
    if(!cHeaderAOD){
      PostData(1, fOutputContainer);
      return ;      
    }
    cHeaderAOD->GetVertex(vtx5);
  }
  FillHistogram("hZvertex"      ,vtx5[2]);
  if (TMath::Abs(vtx5[2]) > 10. ){
    PostData(1, fOutputContainer);
    return ;
  }
  FillHistogram("hSelEvents",2) ;
  //Vtx class z-bin
  Int_t zvtx = TMath::Min(9,Int_t((vtx5[2]+10.)/2.)) ; 
  

  Double_t mcVtxX=0.,mcVtxY=0.,mcVtxZ=0.;
  if(fStack){
    const AliVVertex* primVtxMC   = MCEvent()->GetPrimaryVertex();
    mcVtxX   = primVtxMC->GetX();
    mcVtxY   = primVtxMC->GetY();
    mcVtxZ   = primVtxMC->GetZ();
  }
  
  //Number of contributors
  Int_t nPrimContributors = event->GetPrimaryVertex()->GetNContributors() ;
  if (!fIsMC && (nPrimContributors < 1)){
    PostData(1, fOutputContainer);
    return;   
  }
  FillHistogram("hSelEvents",3) ;
  
  if(fMCType==kJetJet){
    //reject events with too hard jets compared to ptHardBin  
    if(!AcceptJJevent()){
      FillHistogram("hSelEvents",13) ;
      PostData(1, fOutputContainer);
      return;   
    }
  }
  
  

  if(!SelectCentrality(event)){
    PostData(1, fOutputContainer);
    return ;
  }
  FillHistogram("hCentrality",fCentrality,fCentWeight) ;
  FillHistogram("hCentralityRaw",fCentrality) ;
  FillHistogram("hSelEvents",6) ;
    
//   AliEventplane *eventPlane = event->GetEventplane();
//   if( ! eventPlane ) { //Event has no event plane
//     PostData(1, fOutputContainer);
//     return;
//   }
//   //V0A
//   const Int_t harmonics = 2; 
//   Double_t qx=0., qy=0.;  
//   //Whole V0
//   fRP = eventPlane->CalculateVZEROEventPlane(event,10, harmonics,qx,qy);
//   
//   FillHistogram("hRPCentrality",fCentrality,fRP,fCentWeight) ;
  FillHistogram("hSelEvents",7) ;
  
  
  //Calculate charged multiplicity
  Int_t trackMult = 0;
  if(fTrackEvent)
    fTrackEvent->Clear() ;
  else
    fTrackEvent = new TClonesArray("AliCaloPhoton",event->GetNumberOfTracks()) ;

  for (Int_t i=0;i<event->GetNumberOfTracks();++i) {
    AliAODTrack *track = (AliAODTrack*)event->GetTrack(i) ;
    if(fTrackSelection==kLHC13x){
      if(!track->IsOn(AliVTrack::kITSpureSA))
        continue ;
      if(TMath::Abs(track->Eta())> 0.8)
        continue ;
      if(track->Pt() < 0.15 || track->Pt()>10.) continue ;
      if(track->GetITSNcls()<4) continue ; 
      if(track->GetITSchi2()> 36.*track->GetITSNcls()) continue ; 
      float dr, dz;
      track->GetImpactParameters(dr, dz);
      // Check pointing to the primary vertex
      if(TMath::Abs(dr) > 3.2) continue;
      if(TMath::Abs(dz) > 2.4) continue;
    }
    if(fTrackSelection==kFAST || fTrackSelection ==kCENTwoSSD || fTrackSelection==kCENTwSSD){
      UInt_t bit = 1 << 5;
      if (!track->TestFilterBit(bit))
         continue ;
      if(TMath::Abs(track->Eta())> 0.8)
        continue ;
      if(track->Pt() < 0.15 || track->Pt()>10.) continue ;  
    }
    
    if(trackMult>=fTrackEvent->GetSize())
	fTrackEvent->Expand(2*trackMult) ;
    new ((*fTrackEvent)[trackMult]) AliCaloPhoton(track->Px(),track->Py(),track->Pz(),track->P());
      trackMult++;
    FillHistogram("hTrackEtaPhi",track->Eta(),track->Phi(),fCentWeight) ;
    FillHistogram("hTrackEtaPt",track->Eta(),track->Pt(),fCentWeight) ;
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
    
    if (clu->GetType() !=AliVCluster::kPHOSNeutral ) 
      continue ; 
    
    if(clu->E()<fEminCut) 
      continue;
    
    if(fCutType ==kDefCut){
      if(clu->GetNCells()<3)
        continue ;          
    
      if(clu->GetM02()<0.1) 
        continue ;          
        
    }
    if(fCutType ==kLowECut){
      if(clu->E()>2. && clu->GetNCells()<3)
        continue ;          
    
      if(clu->E()>2. && clu->GetM02()<0.1) 
        continue ;          
    }
    
//     if(clu->GetMCEnergyFraction()>kEcrossCut) //Ecross cut, should be filled with Tender
//      continue ;    
//     
//     if(clu->GetDistanceToBadChannel()<fMinBCDistance)
//       continue ;

    Float_t pos[3] ;
    clu->GetPosition(pos) ;
    Int_t fidArea=GetFiducialArea(pos) ;
    
    TVector3 global1(pos) ;
    Int_t relId[4] ;
    fPHOSgeom->GlobalPos2RelId(global1,relId) ;
    Int_t mod  = relId[0] ;
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;
    
    TVector3 local ;
    fPHOSgeom->Global2Local(local,global1,mod);
    
    FillHistogram("hTOF",clu->E(),clu->GetTOF()) ;
    FillHistogram(Form("hTofM%d",mod),clu->GetTOF(),clu->E()) ;
//     if((!fIsMC) && (clu->GetTOF() < kTOFMinCut || clu->GetTOF() > kTOFMaxCut))
//       continue ;          
    if((!fIsMC || (fIsMC && fIsMCWithPileup)) && (TMath::Abs(clu->GetTOF()) > fTimeCut))
      continue ;          
    
    
    TLorentzVector momentum ;
    clu->GetMomentum(momentum, vtx5);
    Double_t cluE = NonLinearity(clu->E()); 
    if(fIsMC){
       Int_t primLabel=clu->GetLabelAt(0) ; //FindPrimary(clu,sure) ;
       if(primLabel>-1){
         AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(primLabel) ;
         cluE = Unsmear(cluE,prim->E()) ;
       }
    }
    Double_t sc = cluE/clu->E(); 
    AliCaloPhoton *p = new ((*fPHOSEvent)[inList]) AliCaloPhoton(sc*momentum.Px(),sc*momentum.Py(),sc*momentum.Pz(),sc*clu->E() );
    inList++;
    
    p->SetModule(mod) ;
//Isolation should be done after tagging
//     Int_t isolation = EvalIsolation(&momentum,kTRUE) ;
//     p->SetIsolationTag(isolation) ;
    
    p->SetDistToBad((Int_t)(1.+clu->GetDistanceToBadChannel()/2.2));
    p->SetBC(i) ; //reference to CaloCluster
    p->SetTagInfo(0); //No pi0 partners found so far
    p->SetTagged(kFALSE);   //Reconstructed pairs found
    p->SetEMCx(local.X()) ;
    p->SetEMCz(local.Z()) ;
    p->SetCluster(clu); 
    
    p->SetFiducialArea(fidArea) ;

    //Mark photons fired trigger
    if(!fIsMB){   
      if(fIsMC){
        if(fPHOSTrigger==kPHOSAny || fPHOSTrigger==kPHOSL0)    
          p->SetTrig(fPHOSTrigUtils->IsFiredTriggerMC(clu)&1) ; 
        else
          p->SetTrig(fPHOSTrigUtils->IsFiredTriggerMC(clu)&(1<<(fPHOSTrigger-1))) ; 
      }
      else
        p->SetTrig(fPHOSTrigUtils->IsFiredTrigger(clu)) ;    
    }
    if(fIsMB || ((!fIsMB) && p->IsTrig()) ){
      FillHistogram(Form("hCluNXZM%d",mod),cellX,cellZ,1.);
      FillHistogram(Form("hCluEXZM%d",mod),cellX,cellZ,cluE);
    }
    
    if(fIsMC){    
       //Check trigger efficiency
       FillHistogram(Form("hMCMinBiasPhot%d",mod),cluE,fCentWeight) ;
       //CheckTrigBadMap
       if((!fIsMB) && fPHOSTrigUtils->TestBadMap(mod, cellX,cellZ))
         FillHistogram(Form("hMCMinBiasPhotMap%d",mod),cluE,fCentWeight) ;
       
       if(p->IsTrig())
          FillHistogram(Form("hMCTrigPhot%d",mod),cluE,fCentWeight) ;
       //Look for MC particle entered PHOS
       FillHistogram(Form("LabelsNPrim_cent%d",fCentBin),cluE,float(clu->GetNLabels())) ;
       Int_t primLabel=clu->GetLabelAt(0) ; //FindPrimary(clu,sure) ;
       //Look what particle left vertex
       if(primLabel>-1){
         AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(primLabel) ;
         Int_t iparent=primLabel;
         AliAODMCParticle * parent = prim;
         Double_t r2=(prim->Xv()-mcVtxX)*(prim->Xv()-mcVtxX)+
                (prim->Yv()-mcVtxY)*(prim->Yv()-mcVtxY)+
                (prim->Zv()-mcVtxZ)*(prim->Zv()-mcVtxZ) ;
         while((r2 > rcut*rcut) && (iparent>-1)){
           if(parent->GetMother()<0)
               break;
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
//    p->SetDispBit(clu->Chi2()<2.5*2.5) ;
//  Cut on FullLamdas
    p->SetDispBit(clu->GetDispersion()<2.5*2.5) ;
    p->SetNsigmaFullDisp(TMath::Sqrt(clu->GetDispersion())) ;
// printf("Disp=%e \n",p->GetNsigmaFullDisp()) ;    
    p->SetTOFBit(TestTOF(clu->GetTOF(),cluE)) ;
    p->SetCPVBit(clu->GetEmcCpvDistance()>2.5) ;  
    p->SetNsigmaCPV(clu->GetEmcCpvDistance()) ;

    FillHistogram("hDispE", cluE,clu->GetM02(),clu->GetM20()) ;
    if(p->IsCPVOK())
      FillHistogram("hDispEneu", cluE,clu->GetM02(),clu->GetM20()) ;

    
    if(fIsMB || (!fIsMB && p->IsTrig())){ 
      FillHistogram(Form("hPhotM%d_All",mod),p->Pt(),fCentWeight) ;
      if(p->IsDispOK()){
        FillHistogram(Form("hPhotM%d_Disp",mod),p->Pt(),fCentWeight) ;
        if(p->IsCPVOK()){
          FillHistogram(Form("hPhotM%d_Both",mod),p->Pt(),fCentWeight) ;
        }
      } 
      if(p->IsCPVOK()){
        FillHistogram(Form("hPhotM%d_CPV",mod),p->Pt(),fCentWeight) ;
      }
    }
  }
  FillHistogram("hPHOSCentrality",fCentrality,inList+0.5) ;
  
  
  FillMCHistos() ;
  FillTaggingHistos() ;

  //Remove old events
  fCurrentMixedList->AddFirst(fPHOSEvent);
  fPHOSEvent=0x0 ;
  if(fCurrentMixedList->GetSize() > 20){
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
  const AliVVertex* primVtxMC   = MCEvent()->GetPrimaryVertex();
  Double_t mcVtxX   = primVtxMC->GetX();
  Double_t mcVtxY   = primVtxMC->GetY();
  Double_t mcVtxZ   = primVtxMC->GetZ();
  
  Int_t nPrim = fStack->GetEntriesFast() ;
  //Fill Primary particl yields
  
  for(Int_t i=0;i<nPrim;i++){
    AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(i) ;
    Double_t r2=(prim->Xv()-mcVtxX)*(prim->Xv()-mcVtxX)+
                (prim->Yv()-mcVtxY)*(prim->Yv()-mcVtxY)+
                (prim->Zv()-mcVtxZ)*(prim->Zv()-mcVtxZ) ;
    if(r2>rcut*rcut){
      continue ;      
    }
    if(fIsMCWithPileup){ // remove pileup
      if(TMath::Abs(prim->T()) > fTimeCut)
        continue ;
    }

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
    Double_t w = fCentWeight*PrimaryParticleWeight(prim) ;
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
    
   //photon should be trigger in PHOS triggered events
   if(!fIsMB && !p->IsTrig() ) 
     continue ;
        
    Double_t w1TOF = 1.; 
    if(fIsMC){ //simulate TOF cut efficiency
      w1TOF=TOFCutEff(p->Pt()) ; 
    }
    
    Int_t label=p->GetPrimary() ;
    if(label<0){ //No label!
      FillHistogram("hMCRecNoLabel",p->Pt(),w1TOF*p->GetWeight());
      continue ;
    }     

    
    AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(label) ;
    //Look what particle left virtex
    Int_t iparent=label;
    AliAODMCParticle * parent = prim;
    while((parent->Xv()-mcVtxX)*(parent->Xv()-mcVtxX)+
          (parent->Yv()-mcVtxY)*(parent->Yv()-mcVtxY)+
          (parent->Zv()-mcVtxZ)*(parent->Zv()-mcVtxZ)> rcut*rcut){
	iparent=parent->GetMother();
	if(iparent<0)
	  break ;
        AliAODMCParticle* tmp = (AliAODMCParticle*)fStack->At(iparent) ;
        int pdgtmp= abs(tmp->GetPdgCode()) ;
        if(pdgtmp==3122 || pdgtmp==3112 || pdgtmp==3212 || pdgtmp==3222) //Lambda, Sigmas
          break ; //keep proton/pion as primary  
	parent = tmp ;	
    }
      Int_t parentPDG=parent->GetPdgCode() ;    
      switch(parentPDG){
	case 22: //electron/positron conversion
        case 111: //Bug in assigning label to cluster
        case 221: 
	  FillPIDHistogramsW("hMCRecPhoton",p,w1TOF);  //Reconstructed with photon from conversion primary
          {
          bool isOnly=kTRUE;
          //Check if other contributions from same primary?
          UInt_t iL=1;
          while(iL<p->GetCluster()->GetNLabels()){
            AliAODMCParticle * testPrim1=prim ;
            bool isSame=kFALSE;
            while(testPrim1){
              int lpr2= p->GetCluster()->GetLabelAt(iL) ;
              AliAODMCParticle * prim2 = nullptr;
              if(lpr2>-1)
                prim2 = (AliAODMCParticle*)fStack->At(lpr2) ;
              while(prim2){
                if(testPrim1==prim2){//same parent  
                   isSame=kTRUE;
                   break ;
                }
                lpr2 = prim2->GetMother();
                if(lpr2>-1)
                  prim2 = (AliAODMCParticle*)fStack->At(lpr2) ;
                else
                  prim2 = nullptr;  
                
              } 
              if(isSame)
                break ;
              Int_t lpr1 = testPrim1->GetMother();
                if(lpr1>-1)
                  testPrim1 = (AliAODMCParticle*)fStack->At(lpr1) ;
                else
                  testPrim1 = nullptr;  
            }
            if(isSame){
              if(testPrim1){ //same but not photon?
                int pdgPT=testPrim1->GetPdgCode() ;
//                 if(pdgPT!=22 && abs(pdgPT)!=11 && pdgPT!=111 && pdgPT!=221){
                if(pdgPT!=22 && abs(pdgPT)!=11 ){
                  isOnly=kFALSE;
                  break ;
                }
              }
              iL++ ;
            }
            else{
              isOnly=kFALSE;
              break ;
            }
          }
          if(isOnly){  
	    FillPIDHistogramsW("hMCRecPhotonOnly",p,w1TOF);  //single photon
          }
          }
	  break ;
	case  11:
	case -11: //electron/positron conversion
	  FillPIDHistogramsW("hMCRecE",p,w1TOF);  //Reconstructed with photon from conversion primary
	  break ;
	case -2212:
	  FillPIDHistogramsW("hMCRecPbar",p,w1TOF);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
	case -2112: //antineutron & antiproton conversion
	  FillPIDHistogramsW("hMCRecNbar",p,w1TOF);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
	case  211:
	case -211:
	  FillPIDHistogramsW("hMCRecPipm",p,w1TOF);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
	case 2212:
	  FillPIDHistogramsW("hMCRecP",p,w1TOF);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
	case  321:
	case -321:
	  FillPIDHistogramsW("hMCRecKpm",p,w1TOF);  //Reconstructed with photon from conversion primary
	  break ;
	case 310:
	  FillPIDHistogramsW("hMCRecK0s",p,w1TOF);  //Reconstructed with photon from conversion primary
	  break ;
	case 130:
	  FillPIDHistogramsW("hMCRecK0l",p,w1TOF);  //Reconstructed with photon from conversion primary
	  break ;
	case 2112: //antineutron & antiproton conversion
	  FillPIDHistogramsW("hMCRecN",p,w1TOF);  //Reconstructed with photon from antibaryon annihilation
	  break ;	  
	case -1: //direct photon or no primary
	  FillPIDHistogramsW("hMCRecNoPRim",p,w1TOF);
	  break ;	  
	default:  
	  if(parent->Charge()!=0)
	    FillPIDHistogramsW("hMCRecCharg",p,w1TOF);  //Reconstructed with photon from antibaryon annihilation
	  else 
	    FillPIDHistogramsW("hMCRecNeutral",p,w1TOF);  //Reconstructed with photon from antibaryon annihilation
      }  
    
    
      //Now classify decay photon
      if(parentPDG==22 || parentPDG==111 || parentPDG==221){
	Int_t iGrandParent=parent->GetMother();
	if(iGrandParent<0 || iGrandParent>=fStack->GetEntriesFast()){
	  FillPIDHistogramsW("hMCRecPhotNoPrim",p,w1TOF);
          continue ;	  
	}
	AliAODMCParticle * grandParent = (AliAODMCParticle*)fStack->At(iGrandParent) ;	
        Int_t grandParentPDG=grandParent->GetPdgCode() ;     
        switch(grandParentPDG){
	case 111: //pi0
	  FillPIDHistogramsW("hMCRecPhotPi0",p,w1TOF);
	  break ;  		
	case 221: //eta decay
	  FillPIDHistogramsW("hMCRecPhotEta",p,w1TOF);
	  break ;  
	case 223: //omega meson decay
	  FillPIDHistogramsW("hMCRecPhotOmega",p,w1TOF);
	  break ;
	default:
	  FillPIDHistogramsW("hMCRecPhotOther",p,w1TOF);
	}
	//--------consider pi0 decays--------------------
	if(grandParentPDG==111){
	  //First find which daughter is our cluster
          //iparent - index of curent photon	  
	  Int_t ipartner=grandParent->GetDaughterLabel(0) ;
	  if(ipartner==iparent){//look for other
  	    if(grandParent->GetNDaughters()>1){
	      ipartner=grandParent->GetDaughterLabel(1);  
	    }
	    else{
	      ipartner=-1 ;
	    }
	  }
 	  //There is no partner in stack
	  if(ipartner==-1){
            FillPIDHistogramsW("hMCDecWMisPartnStack",p,w1TOF) ;
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
                FillPIDHistogramsW("hMCDecWRecUniqPartn",p,w1TOF) ;	
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
	            FillPIDHistogramsW("hMCDecMerged",p,w1TOF) ;
		  }
		  else{
	            FillPIDHistogramsW("hMCDecUnfolded",p,w1TOF) ;
		  }
                }
	      }
	    }	    

	    if(pp){
	      //Partner reconstructed, but did not pass cuts
                FillPIDHistogramsW("hMCDecWRecPartn",p,w1TOF) ;	
    	        Double_t invMass=(*p+ *pp).M() ;
                for(Int_t eminType=0; eminType<3; eminType++){
                  if(pp->E()>fEminCut+0.1*eminType){
	             FillHistogram(Form("hMCmass%d_cent%d",eminType,fCentBin),invMass,p->Pt(),p->GetWeight()*w1TOF*TOFCutEff(pp->Pt())) ;
                  }
                }
		Double_t nSigma=InPi0Band(invMass,p->Pt()) ;
		// analog to Tag
                for(Int_t eminType=0; eminType<3; eminType++){
                  if(pp->E()>fEminCut+0.1*eminType){
  	            for(Int_t isigma=0; isigma<3; isigma++){
  	              if(nSigma<1.+isigma){
			 Int_t iType=3*eminType+isigma ;
	                 FillPIDHistogramsW(Form("hMCDecWithFoundPartnType%d",iType),p,w1TOF*TOFCutEff(pp->Pt())) ;
		      }
		    }
		  }
	        }
	        if(nSigma>3.){
	          FillPIDHistogramsW("hMCDecWithWrongMass",p,w1TOF) ;
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
		  FillPIDHistogramsW("hMCDecWMisPartnAccept",p,w1TOF) ;  //Spectrum of tagged with missed partner
		  isPartnerLost=kTRUE;
		}
		
                if(!isPartnerLost){
		  //this photon is converted before it is registered
		  if(partner->GetNDaughters()>0){
		    AliAODMCParticle* tmpP=(AliAODMCParticle*)fStack->At(partner->GetDaughterLabel(0));
		    if(tmpP->Xv()*tmpP->Xv()+tmpP->Yv()*tmpP->Yv()<450.*450.){  
		      FillPIDHistogramsW("hMCDecWMisPartnConv",p,w1TOF) ;  //Spectrum of tagged with missed partner
		      isPartnerLost=kTRUE;
		    }
		  }
 		}
		if(!isPartnerLost && 
		   partner->E()<0.3){ //energy is not enough to be registered by PHOS
		  FillPIDHistogramsW("hMCDecWMisPartnEmin",p,w1TOF) ;  //Spectrum of tagged with missed partner
		  isPartnerLost=kTRUE;
		}
		if(!isPartnerLost){ //Reason not found!!!!!                  		  
		  FillPIDHistogramsW("hMCDecWMisPartnOther",p,w1TOF);
                  Int_t multClust = event->GetNumberOfCaloClusters();
                  for (Int_t iclu=0; (iclu<multClust) && (!isPartnerLost); iclu++) {
                    AliVCluster * clu = event->GetCaloCluster(iclu);
                    if(!clu->IsPHOS())
                      continue ; 
		    if(NonLinearity(clu->E())==p->E()) //same cluster as current
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
		    FillPIDHistogramsW("hMCDecWMisPartnDefCuts",p,w1TOF);
		  }		  
		}
		else{//Sum of all missed partners
		  FillPIDHistogramsW("hMCDecWMisPartnAll",p,w1TOF);
		}
	      }//Partner - photon
	      else{//partner not photon
		FillPIDHistogramsW("hMCDecWMisPartnNPhot",p,w1TOF);                
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

  //Default isolation cut  
  const Int_t kDefISolation =2;  
  //Minimal opening angle cut  
  const Double_t minAngle=TMath::ATan(2.2*sqrt(2.)/460.) ;   

  const Int_t n=fPHOSEvent->GetEntriesFast() ;
  for(Int_t i=0;i<n;i++){
    AliCaloPhoton *p = static_cast<AliCaloPhoton*>(fPHOSEvent->At(i));
    
    //Isolation tags for all photons
    Int_t isolation = EvalIsolation(p,kTRUE) ;
    p->SetIsolationTag(isolation) ;
  }
     
  //First fill invariant mass distrubutions and mark tagged photons
  //Invariant Mass analysis
  for(Int_t i=0;i<n-1;i++){
    AliCaloPhoton *p1 = static_cast<AliCaloPhoton*>(fPHOSEvent->At(i));
    Double_t ptP1 = p1->Pt() ;
//     Double_t w1=fCentWeight*p1->GetWeight() ;
    Double_t w1TOF = 1.; 
    if(fIsMC){ //simulate TOF cut efficiency
      w1TOF=TOFCutEff(ptP1) ; 
    }
    for(Int_t j = i+1 ; j < n ; j++) {
      AliCaloPhoton * p2 = static_cast<AliCaloPhoton*>(fPHOSEvent->At(j));
      
      if(p1->Angle(p2->Vect())< minAngle)
          continue ;
          
      //At least one photon should be trigger in PHOS triggered events
      if(!fIsMB && !p1->IsTrig() &&  !p2->IsTrig() ) 
        continue ;
      
      Double_t invMass = (*p1 + *p2).M();   
      Double_t ptPi = (*p1 + *p2).Pt() ;
      Double_t ptP2 = p2->Pt() ;
//       Double_t w2=fCentWeight*p2->GetWeight() ;
      Double_t w2TOF=1.;
      Double_t w=TMath::Sqrt(p1->GetWeight()*p2->GetWeight()) ;
      Int_t grandparent,grandParentPDG;
      Double_t ptGrandParent=0;
      Int_t commonParent = IsSameParent(p1,p2,grandparent);
      if(fIsMC ){ //simulate TOF cut efficiency
        w2TOF=TOFCutEff(ptP2); 
        w*=w1TOF*w2TOF; 
        if(grandparent>-1){
          AliAODMCParticle * tmp = ((AliAODMCParticle*)fStack->At(grandparent)) ;
          grandParentPDG = tmp->GetPdgCode() ;
          ptGrandParent = tmp->Pt() ;
        }
      }

      Double_t nsigma1 = InPi0Band(invMass,p1->Pt()) ; //in band with n sigmas
      Double_t nsigma2 = InPi0Band(invMass,p2->Pt()) ; //in band with n sigmas

      if((p1->E()>fEminCut) && (p2->E()>fEminCut)){
        for(Int_t iPID=0; iPID<fNPID; iPID++){  
          if(TestPID(iPID, p1,p2)){
            fhRe[0][fCentBin][iPID]->Fill(invMass,ptPi,w) ;  
            if(commonParent==111){
              fhReTruePi0[0][fCentBin][iPID]->Fill(invMass,ptPi,w) ;  
              //grandparent
              //K0s,Kpm,Nbar
              if(p2->E()>fEminCut){
                fhReSingleTruePi0[0][iPID]->Fill(invMass,ptP1,w) ;  
                if(grandParentPDG==130){
                  fhReSingleTruePi0K0s[0][iPID]->Fill(invMass,ptP1,w) ; 
                  double wcorr = DPMJETKCorr(ptGrandParent) ;
                  fhReSingleTruePi0K0sC[0][iPID]->Fill(invMass,ptP1,w*wcorr) ; 
                }
                if(TMath::Abs(grandParentPDG)==321){
                  fhReSingleTruePi0Kpm[0][iPID]->Fill(invMass,ptP1,w) ;                      
                  double wcorr = DPMJETKCorr(ptGrandParent) ;
                  fhReSingleTruePi0KpmC[0][iPID]->Fill(invMass,ptP1,w*wcorr) ; 
                }
                if(grandParentPDG==-2212 || grandParentPDG == -2112){ //nbar, pbar
                  fhReSingleTruePi0Nbar[0][iPID]->Fill(invMass,ptP1,w) ;                      
                  double wcorr = DPMJETpCorr(ptGrandParent) ;
                  fhReSingleTruePi0NbarC[0][iPID]->Fill(invMass,ptP1,w*wcorr) ; 
                }
                if(p2->E()>fEminCut+0.1){
                  fhReSingleTruePi0[1][iPID]->Fill(invMass,ptP1,w) ;  
                  if(grandParentPDG==130){
                    fhReSingleTruePi0K0s[1][iPID]->Fill(invMass,ptP1,w) ;                      
                    double wcorr = DPMJETKCorr(ptGrandParent) ;
                    fhReSingleTruePi0K0sC[1][iPID]->Fill(invMass,ptP1,w*wcorr) ;                      
                  }
                  if(TMath::Abs(grandParentPDG)==321){
                    fhReSingleTruePi0Kpm[1][iPID]->Fill(invMass,ptP1,w) ;                      
                    double wcorr = DPMJETKCorr(ptGrandParent) ;
                    fhReSingleTruePi0KpmC[1][iPID]->Fill(invMass,ptP1,w*wcorr) ;                      
                  }
                  if(grandParentPDG==-2212 || grandParentPDG == -2112){ //nbar, pbar
                    fhReSingleTruePi0Nbar[1][iPID]->Fill(invMass,ptP1,w) ;                      
                    double wcorr = DPMJETpCorr(ptGrandParent) ;
                    fhReSingleTruePi0NbarC[1][iPID]->Fill(invMass,ptP1,w*wcorr) ;                      
                  }
                  if(p2->E()>fEminCut+0.2){
                    fhReSingleTruePi0[2][iPID]->Fill(invMass,ptP1,w) ;  
                    if(grandParentPDG==130){
                      fhReSingleTruePi0K0s[2][iPID]->Fill(invMass,ptP1,w) ;                      
                      double wcorr = DPMJETKCorr(ptGrandParent) ;
                      fhReSingleTruePi0K0sC[2][iPID]->Fill(invMass,ptP1,w*wcorr) ;                      
                    }
                    if(TMath::Abs(grandParentPDG)==321){
                      fhReSingleTruePi0Kpm[2][iPID]->Fill(invMass,ptP1,w) ;                      
                      double wcorr = DPMJETKCorr(ptGrandParent) ;
                      fhReSingleTruePi0KpmC[2][iPID]->Fill(invMass,ptP1,w*wcorr) ;                      
                    }
                    if(grandParentPDG==-2212 || grandParentPDG == -2112){ //nbar, pbar
                      fhReSingleTruePi0Nbar[2][iPID]->Fill(invMass,ptP1,w) ;                      
                      double wcorr = DPMJETpCorr(ptGrandParent) ;
                      fhReSingleTruePi0NbarC[2][iPID]->Fill(invMass,ptP1,w*wcorr) ;                      
                    }
                  }
                }
              }
              if(p1->E()>fEminCut){
                fhReSingleTruePi0[0][iPID]->Fill(invMass,ptP2,w) ;  
                if(grandParentPDG==130){
                  fhReSingleTruePi0K0s[0][iPID]->Fill(invMass,ptP2,w) ;                      
                  double wcorr = DPMJETKCorr(ptGrandParent) ;
                  fhReSingleTruePi0K0sC[0][iPID]->Fill(invMass,ptP2,w*wcorr) ;                      
                }
                if(TMath::Abs(grandParentPDG)==321){
                  fhReSingleTruePi0Kpm[0][iPID]->Fill(invMass,ptP2,w) ;                      
                  double wcorr = DPMJETKCorr(ptGrandParent) ;
                  fhReSingleTruePi0KpmC[0][iPID]->Fill(invMass,ptP2,w*wcorr) ;                      
                }
                if(grandParentPDG==-2212 || grandParentPDG == -2112){ //nbar, pbar
                  fhReSingleTruePi0Nbar[0][iPID]->Fill(invMass,ptP2,w) ;                      
                  double wcorr = DPMJETpCorr(ptGrandParent) ;
                  fhReSingleTruePi0NbarC[0][iPID]->Fill(invMass,ptP2,w*wcorr) ;                      
                }
                if(p1->E()>fEminCut+0.1){
                  fhReSingleTruePi0[1][iPID]->Fill(invMass,ptP2,w) ;  
                  if(grandParentPDG==130){
                    fhReSingleTruePi0K0s[1][iPID]->Fill(invMass,ptP2,w) ;                      
                    double wcorr = DPMJETKCorr(ptGrandParent) ;
                    fhReSingleTruePi0K0sC[1][iPID]->Fill(invMass,ptP2,w*wcorr) ;                      
                  }
                  if(TMath::Abs(grandParentPDG)==321){
                    fhReSingleTruePi0Kpm[1][iPID]->Fill(invMass,ptP2,w) ;                      
                    double wcorr = DPMJETKCorr(ptGrandParent) ;
                    fhReSingleTruePi0KpmC[1][iPID]->Fill(invMass,ptP2,w*wcorr) ;                      
                  }
                  if(grandParentPDG==-2212 || grandParentPDG == -2112){ //nbar, pbar
                    fhReSingleTruePi0Nbar[1][iPID]->Fill(invMass,ptP2,w) ;                      
                    double wcorr = DPMJETpCorr(ptGrandParent) ;
                    fhReSingleTruePi0NbarC[1][iPID]->Fill(invMass,ptP2,w*wcorr) ;                      
                  }
                  if(p1->E()>fEminCut+0.2){
                    fhReSingleTruePi0[2][iPID]->Fill(invMass,ptP2,w) ;  
                    if(grandParentPDG==130){
                      fhReSingleTruePi0K0s[2][iPID]->Fill(invMass,ptP2,w) ;                      
                      double wcorr = DPMJETKCorr(ptGrandParent) ;
                      fhReSingleTruePi0K0sC[2][iPID]->Fill(invMass,ptP2,w*wcorr) ;                      
                    }
                    if(TMath::Abs(grandParentPDG)==321){
                      fhReSingleTruePi0Kpm[2][iPID]->Fill(invMass,ptP2,w) ;                      
                      double wcorr = DPMJETKCorr(ptGrandParent) ;
                      fhReSingleTruePi0KpmC[2][iPID]->Fill(invMass,ptP2,w*wcorr) ;                      
                    }
                    if(grandParentPDG==-2212 || grandParentPDG == -2112){ //nbar, pbar
                      fhReSingleTruePi0Nbar[2][iPID]->Fill(invMass,ptP2,w) ;                      
                      double wcorr = DPMJETpCorr(ptGrandParent) ;
                      fhReSingleTruePi0NbarC[2][iPID]->Fill(invMass,ptP2,w*wcorr) ;                      
                    }
                  }
                }
              }              
            }
            if(commonParent==221)
              fhReTrueEta[0][fCentBin][iPID]->Fill(invMass,ptPi,w) ;  
            if((p1->E()>fEminCut+0.1) && (p2->E()>fEminCut+0.1)){
              fhRe[1][fCentBin][iPID]->Fill(invMass,ptPi,w) ;  
              if(commonParent==111)
                fhReTruePi0[1][fCentBin][iPID]->Fill(invMass,ptPi,w) ;  
              if(commonParent==221)
                fhReTrueEta[1][fCentBin][iPID]->Fill(invMass,ptPi,w) ;  
              if((p1->E()>fEminCut+0.2) && (p2->E()>fEminCut+0.2)){
                fhRe[2][fCentBin][iPID]->Fill(invMass,ptPi,w) ;  
                if(commonParent==111)
                  fhReTruePi0[2][fCentBin][iPID]->Fill(invMass,ptPi,w) ;  
                if(commonParent==221)
                  fhReTrueEta[2][fCentBin][iPID]->Fill(invMass,ptPi,w) ;  
              }
            }
          }
        }
        if((p1->E()>fEminCut+0.2) && (p2->E()>fEminCut+0.2)){
          if(p1->Module()==p2->Module()){
	       Double_t ww=fCentWeight*TMath::Sqrt(p1->GetWeight()*p2->GetWeight()) ;
               fhReMod[p1->Module()]->Fill(invMass,ptPi,ww) ; 
	  }
        }
      }
      
      if((p1->E()>fEminCut+0.2) && (p2->E()>fEminCut+0.2)){
        //Fill izolated pi0s
        if(nsigma1<2 || nsigma2<2){ //2 sigma band
          TLorentzVector pi0=*p1+*p2 ;
          Int_t isolation=EvalIsolation(&pi0,0) ;
          for(Int_t kind=0; kind<5; kind++){
             if(isolation&(1<<kind)){
               fhPiIsolation[kind][fCentBin]->Fill(ptPi,fCentWeight) ;  
             }
          }
        }
	    
	double alpha=TMath::Abs(p1->E()-p2->E())/(p1->E()+p2->E()) ;
        double dz = TMath::Abs(p1->EMCz()-p2->EMCz()) ;
        double dx = TMath::Abs(p1->EMCx()-p2->EMCx()) ;
	if(nsigma1<2 || nsigma2<2){ //2 sigma band
	  fhQAAllEpartn->Fill(p1->Pt(),alpha) ;
	  fhQAAllzpartn->Fill(p1->Pt(),dz) ;
	  fhQAAllxpartn->Fill(p1->Pt(),dx) ;
          if(p1->GetIsolationTag()&kDefISolation){//default isolation
	    fhQAIsoEpartn->Fill(p1->Pt(),alpha) ; 
	    fhQAIsozpartn->Fill(p1->Pt(),dz) ; 
	    fhQAIsoxpartn->Fill(p1->Pt(),dx) ; 
          }
	  fhQAAllEpartn->Fill(p2->Pt(),alpha) ;
	  fhQAAllzpartn->Fill(p2->Pt(),dz) ;
	  fhQAAllxpartn->Fill(p2->Pt(),dx) ;
          if(p2->GetIsolationTag()&kDefISolation){
	    fhQAIsoEpartn->Fill(p2->Pt(),alpha) ; //default isolation
	    fhQAIsozpartn->Fill(p2->Pt(),dz) ; 
	    fhQAIsoxpartn->Fill(p2->Pt(),dx) ; 
          }
        }
        else{
	  if(nsigma1<4 || nsigma2<4){ //2 sigma band
	    fhQAAllEpartnBg->Fill(p1->Pt(),alpha) ;
	    fhQAAllzpartnBg->Fill(p1->Pt(),dz) ;
	    fhQAAllxpartnBg->Fill(p1->Pt(),dx) ;
            if(p1->GetIsolationTag()&kDefISolation){
	      fhQAIsoEpartnBg->Fill(p1->Pt(),alpha) ;
	      fhQAIsozpartnBg->Fill(p1->Pt(),dz) ; 
	      fhQAIsoxpartnBg->Fill(p1->Pt(),dx) ; 
            }
	    fhQAAllEpartnBg->Fill(p2->Pt(),alpha) ;
	    fhQAAllzpartnBg->Fill(p2->Pt(),dz) ;
	    fhQAAllxpartnBg->Fill(p2->Pt(),dx) ;
            if(p2->GetIsolationTag()&kDefISolation){
	      fhQAIsoEpartnBg->Fill(p2->Pt(),alpha) ;
	      fhQAIsozpartnBg->Fill(p2->Pt(),dz) ; 
	      fhQAIsoxpartnBg->Fill(p2->Pt(),dx) ; 
            }
          }
        }
      }
      
      
      if(p2->E()>fEminCut){
        for(Int_t iPID=0; iPID<fNPID; iPID++){  
          if(TestPID(iPID,p1)){
            fhReSingle[0][fCentBin][iPID]->Fill(invMass,ptP1,w) ;  
            if(p1->GetIsolationTag()&kDefISolation){
              fhReSingleIso[0][fCentBin][iPID]->Fill(invMass,ptP1,w) ;  
            }
            if(p2->E()>fEminCut+0.1){
              fhReSingle[1][fCentBin][iPID]->Fill(invMass,ptP1,w) ;  
              if(p1->GetIsolationTag()&kDefISolation){
                fhReSingleIso[1][fCentBin][iPID]->Fill(invMass,ptP1,w) ;  
              }
              if(p2->E()>fEminCut+0.2){
                fhReSingle[2][fCentBin][iPID]->Fill(invMass,ptP1,w) ;  
                if(p1->GetIsolationTag()&kDefISolation){
                  fhReSingleIso[2][fCentBin][iPID]->Fill(invMass,ptP1,w) ;  
                }
              }
            }
          }
        }
      }

      if(p1->E()>fEminCut){
        for(Int_t iPID=0; iPID<fNPID; iPID++){  
          if(TestPID(iPID,p2)){
            fhReSingle[0][fCentBin][iPID]->Fill(invMass,ptP2,w) ;  
            if(p2->GetIsolationTag()&kDefISolation){
              fhReSingleIso[0][fCentBin][iPID]->Fill(invMass,ptP2,w) ;  
            }
            if(p1->E()>fEminCut+0.1){
              fhReSingle[1][fCentBin][iPID]->Fill(invMass,ptP2,w) ;  
              if(p2->GetIsolationTag()&kDefISolation){
                fhReSingleIso[1][fCentBin][iPID]->Fill(invMass,ptP2,w) ;  
              }
              if(p1->E()>fEminCut+0.2){
                fhReSingle[2][fCentBin][iPID]->Fill(invMass,ptP2,w) ;  
                if(p2->GetIsolationTag()&kDefISolation){
                  fhReSingleIso[2][fCentBin][iPID]->Fill(invMass,ptP2,w) ;  
                }
              }
            }
          }
        }
      }

      
      //Tagging: 1,2,3 sigma
      //Emin=100,200,300 Mev
      //InPi0Band(mass, Ptphoton, type Emin cut
      //Set Bits to 1 if tagged
      Int_t tag1=0 ;
      for(Int_t eminType=0; eminType<3; eminType++){
        if(p2->E()>fEminCut+0.1*eminType){
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
           //Calculate combined probability to pass
           Float_t wtofOld = p1->GetTagWeight(ibit) ;
//            wtofOld=1.-(1.-wtofOld)*(1.-w2TOF);
           wtofOld+=w2TOF;
           p1->SetTagWeight(ibit,wtofOld) ;  //Add weight of TOF cut of the partner
               
           FillPIDHistogramsW(Form("hPhot_TaggedMult%d",ibit),p1,w1TOF*w2TOF) ;
           if(p1->GetIsolationTag()&kDefISolation){               
             FillPIDHistogramsW(Form("hPhot_TaggedMult%d_Isolation2",ibit),p1,w1TOF*w2TOF) ;
           }
	}
	else{ //single tag so far, just add weight
          if((tag1 & (1<<ibit))!=0){
            p1->SetTagWeight(ibit,w2TOF) ;  //Add weight of TOF cut of the partner
          }
        }
      }
      tag1=tag1|oldTag1 ;
      p1->SetTagInfo(tag1) ;
      
      
      Int_t tag2=0 ;
      for(Int_t eminType=0; eminType<3; eminType++){
        if(p1->E()>fEminCut+0.1*eminType){
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
           Float_t wtofOld = p2->GetTagWeight(ibit) ;
//            wtofOld=1.-(1.-wtofOld)*(1.-w1TOF);
           wtofOld+=w1TOF;
           p2->SetTagWeight(ibit,wtofOld) ;  //Add weight of TOF cut of the partner
           FillPIDHistogramsW(Form("hPhot_TaggedMult%d",ibit),p2,w1TOF*w2TOF) ;
           if(p2->GetIsolationTag()&kDefISolation){               
             FillPIDHistogramsW(Form("hPhot_TaggedMult%d_Isolation2",ibit),p2,w1TOF*w2TOF) ;
           }
	}
	else{
          if((tag2 & (1<<ibit))!=0){
            p2->SetTagWeight(ibit,w1TOF) ;  
          }
        }
      }
      tag2=tag2|oldTag2 ;
      p2->SetTagInfo(tag2) ;
            
      Int_t gp;
      Bool_t trueTag=(IsSameParent(p1,p2,gp)==111); //true tag
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

    Double_t w1TOF = 1.; 
    if(fIsMC){ //simulate TOF cut efficiency
      w1TOF=TOFCutEff(p->Pt()) ; 
    }

    Int_t isolation=p->GetIsolationTag() ;
   
    //Inclusive spectra
    FillPIDHistogramsW("hPhot",p,w1TOF) ;
    if(fIsMC){
      double wcorr = 1;
      int l = p->GetPrimaryAtVertex() ;
      if(l>-1){
        AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(l) ;
        if(prim->GetPdgCode()==130 || prim->GetPdgCode()==310 || TMath::Abs(prim->GetPdgCode())==321){
          wcorr=DPMJETKCorr(prim->Pt()) ;
        }
        else{
          if(TMath::Abs(prim->GetPdgCode())==2212 || TMath::Abs(prim->GetPdgCode())==2112 || 
             TMath::Abs(prim->GetPdgCode())==3122){
               wcorr=DPMJETpCorr(prim->Pt()) ;
          }
        }
      }
      FillPIDHistogramsW("hPhotWcorr",p,w1TOF*wcorr) ;
    }    
      
    if(p->DistToBad()>1){
      FillPIDHistogramsW("hPhot_Dist2",p,w1TOF) ;
      if(p->DistToBad()>2){
        FillPIDHistogramsW("hPhot_Dist3",p,w1TOF) ;
      }
    }
      
    //Normal and random isolations      
    Int_t tag=p->GetTagInfo() ;
    for(Int_t kind=0; kind<20; kind++){
      if((isolation&(1<<kind))){
        FillPIDHistogramsW(Form("hPhot_Isolation%d",kind),p,w1TOF) ;
        if((tag & (1<<6))!=0){ //bit6: Emin=300 MeV+1sigma+all partners
          FillPIDHistogramsW(Form("hPhot_Tagged_Isolation%d",kind),p,p->GetTagWeight(6)*w1TOF) ;
	}
	else{
          FillPIDHistogramsW(Form("hPhot_nTagged_Isolation%d",kind),p,w1TOF) ;     
        }
      }
    }
    
   for(Int_t ibit=0; ibit<18; ibit++){
     if((tag & (1<<ibit))!=0){ 
       FillPIDHistogramsW(Form("hPhot_Tagged%d",ibit),p,p->GetTagWeight(ibit)*w1TOF) ;
       if(fIsMC){
         double wcorr = 1;
         int l = p->GetPrimaryAtVertex() ;
         if(l>-1){
           AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(l) ;
           if(prim->GetPdgCode()==130 || prim->GetPdgCode()==310 || TMath::Abs(prim->GetPdgCode())==321){
             wcorr=DPMJETKCorr(prim->Pt()) ;
           }
           else{
             if(TMath::Abs(prim->GetPdgCode())==2212 || TMath::Abs(prim->GetPdgCode())==2112 || 
                TMath::Abs(prim->GetPdgCode())==3122){
               wcorr=DPMJETpCorr(prim->Pt()) ;
             }
           }
         }
         FillPIDHistogramsW(Form("hPhotWcorr_Tagged%d",ibit),p,p->GetTagWeight(ibit)*w1TOF*wcorr) ;           
       }
       if(isolation&kDefISolation){
          FillPIDHistogramsW(Form("hPhot_Tagged%d_Isolation2",ibit),p,p->GetTagWeight(ibit)*w1TOF) ;           
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
        
       if(p1->Angle(p2->Vect())< minAngle)
          continue ;       
        
        Double_t invMass = (*p1 + *p2).M();
        Double_t ptPi = (*p1 + *p2).Pt() ;
        Double_t ptP1 = p1->Pt() ;
        Double_t ptP2 = p2->Pt() ;
        Double_t w=fCentWeight*TMath::Sqrt(p1->GetWeight()*p2->GetWeight()) ;
//         Double_t w1=fCentWeight*p1->GetWeight() ;
//         Double_t w2=fCentWeight*p2->GetWeight() ;
        Double_t w1TOF = 1.; 
        Double_t w2TOF = 1.; 
        if(fIsMC ){ //simulate TOF cut efficiency
          w1TOF=TOFCutEff(ptP1) ; 
          w2TOF=TOFCutEff(ptP2) ; 
          w*=w1TOF*w2TOF ;
        }
 
        //At least one photon should be trigger in PHOS triggered events
        if(!fIsMB && !p1->IsTrig() &&  !p2->IsTrig() ) 
          continue ;	
	
        if((p1->E()>fEminCut) && (p2->E()>fEminCut)){
          for(Int_t iPID=0; iPID<fNPID; iPID++){  
            if(TestPID(iPID,p1,p2)){
              fhMi[0][fCentBin][iPID]->Fill(invMass,ptPi,w) ;  
              if((p1->E()>fEminCut+0.1) && (p2->E()>fEminCut+0.1)){
                fhMi[1][fCentBin][iPID]->Fill(invMass,ptPi,w) ;  
                if((p1->E()>fEminCut+0.2) && (p2->E()>fEminCut+0.2)){
                  fhMi[2][fCentBin][iPID]->Fill(invMass,ptPi,w) ;  

                }
              }
            }
          }
          if((p1->E()>fEminCut+0.2) && (p2->E()>fEminCut+0.2)){
            if(p1->Module()==p2->Module()){
              fhMiMod[p1->Module()]->Fill(invMass,(*p1 + *p2).Pt(),w) ; 
            }
          }
        }
	
	if(p2->E()>fEminCut){
          for(Int_t iPID=0; iPID<fNPID; iPID++){  
            if(TestPID(iPID,p1)){
              fhMiSingle[0][fCentBin][iPID]->Fill(invMass,ptP1,w) ;  
              if(p1->GetIsolationTag()&kDefISolation){
                fhMiSingleIso[0][fCentBin][iPID]->Fill(invMass,ptP1,w) ;  
              }
	      if(p2->E()>fEminCut+0.1){
                fhMiSingle[1][fCentBin][iPID]->Fill(invMass,ptP1,w) ;  
                if(p1->GetIsolationTag()&kDefISolation){
                  fhMiSingleIso[1][fCentBin][iPID]->Fill(invMass,ptP1,w) ;  
                }
	        if(p2->E()>fEminCut+0.2){
                  fhMiSingle[2][fCentBin][iPID]->Fill(invMass,ptP1,w) ;  
                  if(p1->GetIsolationTag()&kDefISolation){
                    fhMiSingleIso[2][fCentBin][iPID]->Fill(invMass,ptP1,w) ;  
                  }
                }
              }
            }
          }
        }
	
	if(p1->E()>fEminCut){
          for(Int_t iPID=0; iPID<fNPID; iPID++){  
            if(TestPID(iPID,p2)){
              fhMiSingle[0][fCentBin][iPID]->Fill(invMass,ptP2,w) ;  
              if(p2->GetIsolationTag()&kDefISolation){
                fhMiSingleIso[0][fCentBin][iPID]->Fill(invMass,ptP2,w) ;  
              }
	      if(p1->E()>fEminCut+0.1){
                fhMiSingle[1][fCentBin][iPID]->Fill(invMass,ptP2,w) ;  
                if(p2->GetIsolationTag()&kDefISolation){
                  fhMiSingleIso[1][fCentBin][iPID]->Fill(invMass,ptP2,w) ;  
                }
	        if(p1->E()>fEminCut+0.2){
                  fhMiSingle[2][fCentBin][iPID]->Fill(invMass,ptP2,w) ;  
                  if(p2->GetIsolationTag()&kDefISolation){
                    fhMiSingleIso[2][fCentBin][iPID]->Fill(invMass,ptP2,w) ;  
                  }
                }
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
//  Double_t mpi0mean =  0.135 ;  

//   //Parameterization 13.10.2018
//   Double_t mpi0mean =1.34693e-01-3.68195e-04*TMath::TanH((pt-5.00834e+00)/1.81347e+00) ;  
  
  Double_t mpi0mean = 0; 
  Double_t mpi0sigma=1.;
  if(fRunNumber>=265015 && fRunNumber<=267166){ //LHC16qrst, tune 19.06.2020
     mpi0mean = -8.63864e-01+(5.32935e-02+1.18250e-01*pt+2.01664e-01*pt*pt+2.87553e-01*pt*pt*pt+pt*pt*pt*pt)/(5.47622e-02+1.11935e-01*pt+2.11192e-01*pt*pt+2.82832e-01*pt*pt*pt+pt*pt*pt*pt) ;
     mpi0sigma = 3.94524e-04/pt/pt-3.06055e-04/pt+8.07149e-03-2.51764e-03*sqrt(pt)+5.32026e-04*pt ;
  }
  else{
     //Parameterization 21.08.2018 with updated NonLin Run2TuneMC
     mpi0mean = 1.36269e-01-1.81643456e-05/((pt-4.81920e-01)*(pt-4.81920e-01)+3.662247e-02)-2.15520e-04*exp(-pt/1.72016e+00) ;  
     //Parameterization 13.10.2018 with updated NonLin Run2TuneMC
     mpi0sigma=TMath::Sqrt(2.59195e-05+1.101556186e-05/pt+2.e-8*pt*pt) ;
  }   
  
  //Double_t mpi0sigma=TMath::Sqrt(5.22245e-03*5.22245e-03 +2.86851e-03*2.86851e-03/pt) + 9.09932e-05*pt ;
  //Parameterization of data 30.08.2014
//   Double_t mpi0sigma=TMath::Sqrt(4.67491e-03*4.67491e-03 +3.42884e-03*3.42884e-03/pt) + 1.24324e-04*pt ;

//   //Parameterization 13.10.2018
//   Double_t mpi0sigma=TMath::Sqrt(3.79261e-03*3.79261e-03/pt+4.76506e-03*4.76506e-03+4.87152e-05*4.87152e-05*pt*pt*pt) ;
  
  return TMath::Abs(m-mpi0mean)/mpi0sigma ;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskTaggedPhotons::IsSameParent(const AliCaloPhoton *p1, const AliCaloPhoton *p2, Int_t &grandparent)const{
  //Looks through parents and finds if there was commont pi0 among ancestors

  if(!fIsMC){
    grandparent=-1;  
    return 0 ; //can not say anything
  }
  Int_t prim1 = p1->GetPrimary();
  while(prim1!=-1){ 
    Int_t prim2 = p2->GetPrimary();
  
    while(prim2!=-1){       
      if(prim1==prim2){
        AliAODMCParticle * tmp = ((AliAODMCParticle*)fStack->At(prim1)) ;
        grandparent = tmp->GetMother() ;
        return tmp->GetPdgCode() ;
      }
      prim2=((AliAODMCParticle*)fStack->At(prim2))->GetMother() ;
    }
    prim1=((AliAODMCParticle*)fStack->At(prim1))->GetMother() ;
  }
  grandparent=-1;  
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
    if(fRunNumber<209122){ //Run1  
       AliError("TaggedPhotons: Can not get geometry from TENDER, creating PHOS geometry for Run1\n") ;    
       fPHOSgeom = AliPHOSGeometry::GetInstance("IHEP","");
    }
    else{
       AliError("TaggedPhotons: Can not get geometry from TENDER, creating PHOS geometry for Run2\n") ;    
       fPHOSgeom = AliPHOSGeometry::GetInstance("Run2","");
    }
    AliOADBContainer geomContainer("phosGeo");
    geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
    TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
    for(Int_t mod=0; mod<5; mod++) {
      if(!matrixes->At(mod)) continue;
      fPHOSgeom->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;   
    }
  }
  
  //Read BadMap for MC simulations
  AliOADBContainer badmapContainer(Form("phosBadMap"));
  badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root","phosBadMap");
  TObjArray *maps = (TObjArray*)badmapContainer.GetObject(fRunNumber,"phosBadMap");
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
  if(tmp->IsA() == TClass::GetClass("TProfile")){
    ((TProfile*)tmp)->Fill(x,y) ;
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
void AliAnalysisTaskTaggedPhotons::FillPIDHistograms(const char * name,  AliCaloPhoton * p) const{

  double pt= p->Pt(); 
  double w=fCentWeight*p->GetWeight() ;
  FillHistogram(Form("%s_All_cent%d",name,fCentBin),pt,w) ;
  if(p->IsDispOK())
    FillHistogram(Form("%s_Disp_cent%d",name,fCentBin),pt,w) ;
  if(fNPID>4){
      if(p->GetNsigmaFullDisp()<3.){
        FillHistogram(Form("%s_Disp3_cent%d",name,fCentBin),pt,w) ;
        if(p->GetNsigmaFullDisp()<2.)
          FillHistogram(Form("%s_Disp2_cent%d",name,fCentBin),pt,w) ;
      }
  }
  if(p->IsCPVOK()){
    FillHistogram(Form("%s_CPV_cent%d",name,fCentBin),pt,w) ;
    if(p->IsDispOK()) 
      FillHistogram(Form("%s_Both_cent%d",name,fCentBin),pt,w) ;
    if(fNPID>4){
      if(p->GetNsigmaFullDisp()<3.){
        FillHistogram(Form("%s_Both3_cent%d",name,fCentBin),pt,w) ;
        if(p->GetNsigmaFullDisp()<2.)
          FillHistogram(Form("%s_Both2_cent%d",name,fCentBin),pt,w) ;
      }
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::FillPIDHistogramsW(const char * name,  AliCaloPhoton * p, Double_t w) const{

  double pt= p->Pt(); 
  double ww=w*fCentWeight*p->GetWeight() ;
  FillHistogram(Form("%s_All_cent%d",name,fCentBin),pt,ww) ;
  if(p->IsDispOK())
    FillHistogram(Form("%s_Disp_cent%d",name,fCentBin),pt,ww) ;
  if(fNPID>4){
      if(p->GetNsigmaFullDisp()<3.){
        FillHistogram(Form("%s_Disp3_cent%d",name,fCentBin),pt,ww) ;
        if(p->GetNsigmaFullDisp()<2.)
          FillHistogram(Form("%s_Disp2_cent%d",name,fCentBin),pt,ww) ;
      }
  }
  if(p->IsCPVOK()){
    FillHistogram(Form("%s_CPV_cent%d",name,fCentBin),pt,ww) ;
    if(p->IsDispOK()) 
      FillHistogram(Form("%s_Both_cent%d",name,fCentBin),pt,ww) ;
    if(fNPID>4){
      if(p->GetNsigmaFullDisp()<3.){
        FillHistogram(Form("%s_Both3_cent%d",name,fCentBin),pt,ww) ;
        if(p->GetNsigmaFullDisp()<2.)
          FillHistogram(Form("%s_Both2_cent%d",name,fCentBin),pt,ww) ;
      }
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::FillPIDHistograms(const char * name,  AliCaloPhoton * p,Double_t x) const{

  double pt= p->Pt(); 
  double w = fCentWeight*p->GetWeight() ;
  FillHistogram(Form("%s_All_cent%d",name,fCentBin),x,pt,w) ;
  if(p->IsDispOK())
    FillHistogram(Form("%s_Disp_cent%d",name,fCentBin),x,pt,w) ;
  if(fNPID>4){
      if(p->GetNsigmaFullDisp()<3.){
        FillHistogram(Form("%s_Disp3_cent%d",name,fCentBin),pt,w) ;
        if(p->GetNsigmaFullDisp()<2.)
          FillHistogram(Form("%s_Disp2_cent%d",name,fCentBin),pt,w) ;
      }
  }
  if(p->IsCPVOK()){
    FillHistogram(Form("%s_CPV_cent%d",name,fCentBin),x,pt,w) ;
    if(p->IsDispOK() ) 
      FillHistogram(Form("%s_Both_cent%d",name,fCentBin),x,pt,w) ;
    if(fNPID>4){
      if(p->GetNsigmaFullDisp()<3.){
        FillHistogram(Form("%s_Both3_cent%d",name,fCentBin),pt,w) ;
        if(p->GetNsigmaFullDisp()<2.)
          FillHistogram(Form("%s_Both2_cent%d",name,fCentBin),pt,w) ;
      }
    }  
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::FillPIDHistograms(const char * name,  AliCaloPhoton * p1, AliCaloPhoton * p2,Double_t x, Bool_t /*isRe*/) const{

  Double_t ptPi = (*p1 + *p2).Pt() ;
  Double_t w=fCentWeight*TMath::Sqrt(p1->GetWeight()*p2->GetWeight()) ;
//  if(isRe){
//    Int_t pdg=IsSameParent(p1,p2) ;
//    if(pdg!=0)
//      w=p1->GetWeight() ;
//  }
  FillHistogram(Form("%s_All_cent%d",name,fCentBin),x,ptPi,w) ;
  if(p1->IsDispOK() && p2->IsDispOK())
    FillHistogram(Form("%s_Disp_cent%d",name,fCentBin),x,ptPi,w) ;
  if(fNPID>4){
    if(p1->GetNsigmaFullDisp()<3. && p2->GetNsigmaFullDisp()<3.){
      FillHistogram(Form("%s_Disp3_cent%d",name,fCentBin),x,ptPi,w) ;
      if(p1->GetNsigmaFullDisp()<2. && p2->GetNsigmaFullDisp()<2.)
        FillHistogram(Form("%s_Disp2_cent%d",name,fCentBin),x,ptPi,w) ;
    }
  }
  if(p1->IsCPVOK() && p2->IsCPVOK()){
    FillHistogram(Form("%s_CPV_cent%d",name,fCentBin),x,ptPi,w) ;
    if(p1->IsDispOK()  && p2->IsDispOK() ) 
      FillHistogram(Form("%s_Both_cent%d",name,fCentBin),x,ptPi,w) ;
    if(fNPID>4){
      if(p1->GetNsigmaFullDisp()<3. && p2->GetNsigmaFullDisp()<3.){
        FillHistogram(Form("%s_Both3_cent%d",name,fCentBin),x,ptPi,w) ;
        if(p1->GetNsigmaFullDisp()<2. && p2->GetNsigmaFullDisp()<2.)
          FillHistogram(Form("%s_Both2_cent%d",name,fCentBin),x,ptPi,w) ;
      }
    }
  } 
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
   if(!ph) return 0 ;

   const Double_t coneR1=0.3 ;
   const Double_t coneR2=0.4 ;
   const Double_t coneR3=0.5 ;
   
   const Double_t eCut=2.; //GeV

   //Calculate UE cut + threshold
   Double_t a[5]={3.40313,2.14956,1.51481,8.94272e-01,3.96717e-01} ;

   //TODO!!! cut for pp
//   Double_t a[6]={0.49,0.49,0.49,0.49,0.49,0.49} ;  //pp

   const Double_t cutEcone1=eCut+coneR1*coneR1*TMath::Pi()*(a[fCentBin]+0.008*ph->Pt()) ;
   const Double_t cutEcone2=eCut+coneR2*coneR2*TMath::Pi()*(a[fCentBin]+0.008*ph->Pt()) ;
   const Double_t cutEcone3=eCut+coneR3*coneR3*TMath::Pi()*(a[fCentBin]+0.008*ph->Pt()) ;


   //Sum of energies in cones, tracks and clusters in PHOS
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
         
    //Orientation of the perp. cone v.r.t. Reaction Plane     
    Double_t dphiRP= phiTrig+TMath::PiOver2()-fRP ;     
    while(dphiRP<0) dphiRP+=TMath::Pi() ;
    while(dphiRP>TMath::Pi()) dphiRP-=TMath::Pi() ;
     
    
    //Fill QA histgams
    Double_t ptTrig=ph->Pt() ;
    if(isPhoton){
    FillHistogram(Form("QA_Cone1_Tracks_cent%d",fCentBin),ptTrig,eCone1) ;
    FillHistogram(Form("QA_Cone2_Tracks_cent%d",fCentBin),ptTrig,eCone2) ;
    FillHistogram(Form("QA_Cone3_Tracks_cent%d",fCentBin),ptTrig,eCone3) ;
    FillHistogram(Form("QA_profCone1_Tracks_cent%d",fCentBin),ptTrig,eCone1) ;
    FillHistogram(Form("QA_profCone2_Tracks_cent%d",fCentBin),ptTrig,eCone2) ;
    FillHistogram(Form("QA_profCone3_Tracks_cent%d",fCentBin),ptTrig,eCone3) ;
    AliCaloPhoton * photon = (AliCaloPhoton*)ph ;
    if((photon->GetTagInfo()& (1<<6))!=0){ //bit6: Emin=300 MeV+1sigma+all partners
      FillHistogram(Form("QA_TCone1_Tracks_cent%d",fCentBin),ptTrig,eCone1) ;
      FillHistogram(Form("QA_TCone2_Tracks_cent%d",fCentBin),ptTrig,eCone2) ;
      FillHistogram(Form("QA_TCone3_Tracks_cent%d",fCentBin),ptTrig,eCone3) ;      
      FillHistogram(Form("QA_profTCone1_Tracks_cent%d",fCentBin),ptTrig,eCone1) ;
      FillHistogram(Form("QA_profTCone2_Tracks_cent%d",fCentBin),ptTrig,eCone2) ;
      FillHistogram(Form("QA_profTCone3_Tracks_cent%d",fCentBin),ptTrig,eCone3) ;      
    }
    FillHistogram(Form("QA_PCone1_Tracks_cent%d",fCentBin),ptTrig,eP1Cone1) ;
    FillHistogram(Form("QA_PCone2_Tracks_cent%d",fCentBin),ptTrig,eP1Cone2) ;
    FillHistogram(Form("QA_PCone3_Tracks_cent%d",fCentBin),ptTrig,eP1Cone3) ;
    FillHistogram(Form("QA_PCone1_Tracks_cent%d",fCentBin),ptTrig,eP2Cone1) ;
    FillHistogram(Form("QA_PCone2_Tracks_cent%d",fCentBin),ptTrig,eP2Cone2) ;
    FillHistogram(Form("QA_PCone3_Tracks_cent%d",fCentBin),ptTrig,eP2Cone3) ;
    FillHistogram(Form("QA_profPCone1_Tracks_cent%d",fCentBin),ptTrig,eP1Cone1) ;
    FillHistogram(Form("QA_profPCone2_Tracks_cent%d",fCentBin),ptTrig,eP1Cone2) ;
    FillHistogram(Form("QA_profPCone3_Tracks_cent%d",fCentBin),ptTrig,eP1Cone3) ;
    FillHistogram(Form("QA_profPCone1_Tracks_cent%d",fCentBin),ptTrig,eP2Cone1) ;
    FillHistogram(Form("QA_profPCone2_Tracks_cent%d",fCentBin),ptTrig,eP2Cone2) ;
    FillHistogram(Form("QA_profPCone3_Tracks_cent%d",fCentBin),ptTrig,eP2Cone3) ;
    
    FillHistogram(Form("QA_PCone12_Tracks_cent%d",fCentBin),ptTrig,eP1Cone2,eP2Cone2) ;
    FillHistogram(Form("QA_PConeRP_Tracks_cent%d",fCentBin),ptTrig,eP1Cone2,dphiRP) ;
    FillHistogram(Form("QA_PConeRP_Tracks_cent%d",fCentBin),ptTrig,eP2Cone2,dphiRP) ;
    
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
    
    Int_t iCone1E2 = (1.+cutEcone1 > eCone1) ;
    Int_t iCone2E2 = (1.+cutEcone2 > eCone2) ;
    Int_t iCone3E2 = (1.+cutEcone3 > eCone3) ;
    
    //PerpCone1
    Int_t iP1Cone1 = (cutEcone1 > eP1Cone1) ;
    Int_t iP1Cone2 = (cutEcone2 > eP1Cone2) ;
    Int_t iP1Cone3 = (cutEcone3 > eP1Cone3) ;
    
    Int_t iP2Cone1 = (cutEcone1 > eP2Cone1) ;
    Int_t iP2Cone2 = (cutEcone2 > eP2Cone2) ;
    Int_t iP2Cone3 = (cutEcone3 > eP2Cone3) ;
    
    Int_t isolation=     iCone1E1 +   2*iCone2E1 +   4*iCone3E1
                    +  8*iCone1E2 +  16*iCone2E2 +  32*iCone3E2
                    + 64*iP1Cone1 + 128*iP1Cone2 + 256*iP1Cone3
                    +512*iP2Cone1 +1024*iP2Cone2 +2048*iP2Cone3;
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
Bool_t AliAnalysisTaskTaggedPhotons::TestPID(Int_t iPID, AliCaloPhoton* p1, AliCaloPhoton* p2){
//   //Checks PID of particle
  
  if(!p1 || !p2) return kFALSE ;
  if(iPID==0) return kTRUE ;
  if(iPID==1) return p1->IsDispOK() && p2->IsDispOK() ;
  if(iPID==2) return p1->IsCPVOK()  && p2->IsCPVOK() ;
  if(iPID==3) return p1->IsDispOK() && p1->IsCPVOK() && p2->IsDispOK() && p2->IsCPVOK() ;
  return kFALSE ;
    
}
//_________________________________________________________________________________
Double_t AliAnalysisTaskTaggedPhotons::PrimaryParticleWeight(AliAODMCParticle * particle){
  if(!fIsMC)  //This is real data
     return 1;
  //Classify parent at vertex
  //Introduce for eta and pi0 weights   
    
  if(fMCType==kFullMC || fMCType==kJetJet){
    return 1.; //For full MC scan
  }
  Int_t mother = particle->GetMother() ;
  Int_t pdg = particle->GetPdgCode() ;
  Int_t parentPDG=0; 
  if(fMCType==kSingleGamma){
    parentPDG=22; 
  }
  if(fMCType==kSinglePi0){
    parentPDG=111; 
  }
  if(fMCType==kSingleEta){
    parentPDG=221; 
  }
  while(mother>-1 && pdg!=parentPDG){
    particle = (AliAODMCParticle*) fStack->At(mother);
    mother = particle->GetMother() ;
    pdg = particle->GetPdgCode() ;
  }
  if(pdg!=parentPDG){
    return 0. ;
  }
  
  //Particle within 1 cm from the virtex
  Double_t x = particle->Pt() ;
  if(x<0.4) x=0.4 ;
  return fWeightParamPi0[0]*(TMath::Power(x,fWeightParamPi0[1])*
       (fWeightParamPi0[2]+x*fWeightParamPi0[3]+x*x)/
       (fWeightParamPi0[4]+x*fWeightParamPi0[5]+x*x) +fWeightParamPi0[6])  ;
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
    FillHistogram(Form("LabelsGamma_cent%d",fCentBin),NonLinearity(clu->E())) ;
    FillHistogram(Form("LabelsGammaE_cent%d",fCentBin),NonLinearity(clu->E()),eMax/NonLinearity(clu->E())) ;
  }  
  
  for(Int_t i=0;  i<n;  i++){
    AliAODMCParticle*  p= (AliAODMCParticle*) fStack->At(clu->GetLabelAt(i)) ;
    Int_t pdg = p->GetPdgCode() ;
    if(pdg==22  ||  pdg==11 || pdg == -11){
      if(p->E()>emFraction*NonLinearity(clu->E())){
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
//___________________________________________________________________________
Double_t AliAnalysisTaskTaggedPhotons::TrigCentralityWeight(Double_t x){
   Double_t w=1.; 
   if(fRunNumber >=195344 && fRunNumber <= 196311){ //LHC13bcde  
   switch(fCentEstimator){
    case 4 : w = 10877.2 - 567.661*x + 18.1458*x*x -0.31884*x*x*x + 0.00269787*x*x*x*x -8.65923e-06*x*x*x*x*x; //CL1   
             break;
    case 3 : w = 5593.67 -75.5495*x + 1.76812*x*x -0.0514218*x*x*x + 0.000585638*x*x*x*x-2.24744e-06*x*x*x*x*x; //ZNA       
             break;         
    case 2 : w = 8329.24 - 300.106*x + 9.27297*x*x - 0.210205*x*x*x +0.00275734*x*x*x*x-1.91452e-05*x*x*x*x*x+ 5.49793e-08*x*x*x*x*x*x ;//V0M 
             break;
    case 1 : 
    default: w= 2.896558e+04+x*(-2.470793e+03+x*(2.318660e+02+x*(-1.403869e+01+x*(5.216917e-01+x*(-1.218504e-02+    
                x*(1.792322e-04+x*(-1.608286e-06+x*(8.026836e-09-1.70525e-11*x)))))))) ;  //V0A
   }
   }
    if(fRunNumber >=196433 && fRunNumber <= 197388 ){ //LHC13f  
    switch(fCentEstimator){
    case 4 : w=  3.41752 - 0.152255*x+0.00432026*x*x-7.03916e-05*x*x*x +  5.55297e-07*x*x*x*x -1.65227e-09*x*x*x*x*x ; //CL1
             break;
    case 3 : w= 1.88202 - 0.0206302*x + 0.000366297*x*x -1.26288e-05*x*x*x + 1.49687e-07*x*x*x*x -5.6476e-10*x*x*x*x*x ;//ZNC
             break;         
    case 2 : w= 3.13185 - 0.12687*x + 0.00354685*x*x - 5.91787e-05*x*x*x + 4.79696e-07*x*x*x*x - 1.46955e-09*x*x*x*x*x; //V0M
             break;
    case 1 : 
    default:  w= 1.813239e+05+x*(-1.153753e+04+x*(7.810399e+02+x*(-3.481358e+01+x*(9.353658e-01+x*(-1.534967e-02 +     
                        x*(1.513630e-04+x*(-8.353821e-07 +x*(2.110175e-09 -9.88144e-13*x)))))))) ;  //V0C
    }
    }
    if(w<=0.)
       return 1.;
    else
       return 1./w ; 
}
//___________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::SetCentralityWeights(TString filename){
 
  TFile * fin = new TFile(filename.Data()) ;
  if(!fin->IsOpen()){
    return ;  
  }

  TH1F * tmp = (TH1F*)fin->Get("hCentrality2V0A") ;
  if(tmp){  
    gROOT->cd() ;  
    fCentralityWeights[0] = new TH1F(*tmp) ;  
    fCentralityWeights[0]->SetName("WeightV0A") ;
  }
  tmp = (TH1F*)fin->Get("hCentrality2V0C") ;
  if(tmp){  
    gROOT->cd() ;  
    fCentralityWeights[1] = new TH1F(*tmp) ;  
    fCentralityWeights[1]->SetName("WeightV0C") ;  
  }
  tmp = (TH1F*)fin->Get("hCentrality2V0M") ;
  if(tmp){  
    gROOT->cd() ;  
    fCentralityWeights[2] = new TH1F(*tmp) ;  
    fCentralityWeights[2]->SetName("WeightV0M") ;  
  }
  tmp = (TH1F*)fin->Get("hCentrality2ZNA") ;
  if(tmp){  
    gROOT->cd() ;  
    fCentralityWeights[3] = new TH1F(*tmp) ;  
    fCentralityWeights[3]->SetName("WeightZNA") ;  
  }
  tmp = (TH1F*)fin->Get("hCentrality2ZNC") ;
  if(tmp){  
    gROOT->cd() ;  
    fCentralityWeights[4] = new TH1F(*tmp) ;  
    fCentralityWeights[4]->SetName("WeightZNC") ;  
  }
  tmp = (TH1F*)fin->Get("hCentrality2CL1") ;
  if(tmp){  
    gROOT->cd() ;  
    fCentralityWeights[5] = new TH1F(*tmp) ;  
    fCentralityWeights[5]->SetName("WeightCL1") ;  
  }
  fin->Close() ;
  delete fin ;
  
}
//___________________________________________________________________________
Double_t AliAnalysisTaskTaggedPhotons::MBCentralityWeight(Double_t x){   //Correction for Pileup cut centrality bias

  Double_t w=1.;
  if(fRunNumber >=195344 && fRunNumber <= 196311){ //LHC13bcde  
    switch(fCentEstimator){
    case 4 : if(fCentralityWeights[5]) w=fCentralityWeights[5]->GetBinContent(fCentralityWeights[5]->FindBin(x)); //CL1
             break;
    case 3 : if(fCentralityWeights[3]) w=fCentralityWeights[3]->GetBinContent(fCentralityWeights[3]->FindBin(x)); //ZNA
             break;         
    case 2 : if(fCentralityWeights[2]) w=fCentralityWeights[2]->GetBinContent(fCentralityWeights[2]->FindBin(x)); //V0M
             break;
    case 1 : 
    default: if(fCentralityWeights[0]) w=fCentralityWeights[0]->GetBinContent(fCentralityWeights[0]->FindBin(x)); //V0A
    }
  }
  else if(fRunNumber >=196433 && fRunNumber <= 197388 ){ //LHC13f  
    switch(fCentEstimator){
    case 4 : if(fCentralityWeights[5]) w=fCentralityWeights[5]->GetBinContent(fCentralityWeights[5]->FindBin(x)); //CL1
             break;
    case 3 : if(fCentralityWeights[4]) w=fCentralityWeights[4]->GetBinContent(fCentralityWeights[4]->FindBin(x)); //ZNC
             break;         
    case 2 : if(fCentralityWeights[2]) w=fCentralityWeights[2]->GetBinContent(fCentralityWeights[2]->FindBin(x)); //V0M
             break;
    case 1 : 
    default: if(fCentralityWeights[1]) w=fCentralityWeights[1]->GetBinContent(fCentralityWeights[1]->FindBin(x)); //V0C
    }
  }
 
  if(w==0.){
    w=1.;    
  }
  return 1./w ;  
    
}
//___________________________________________________________________________
Bool_t AliAnalysisTaskTaggedPhotons::SelectCentrality(AliVEvent * event){
  //Select centrality either in p-Pb or in pp collisions
  //returns true if no errors occured, false in case of errors  
  
  if(fIsFastMC){ 
    fCentWeight=1.; 
    fCentrality=1.;
    return kTRUE ;
  }
        
    
  //In case of p-Pb data  
  if(fRunNumber >=195344 && fRunNumber <= 197388 ){ //LHC13bcdef: pPb 5.02
    
    //Fill Centrality before vertex/pileup cuts
    AliCentrality *centrality = event->GetCentrality();
    if(!centrality)
      return kFALSE ;  

    Double_t v0A =centrality->GetCentralityPercentile("V0A");
    FillHistogram("hCentrality1V0A",v0A) ;
      
    Double_t v0C =centrality->GetCentralityPercentile("V0C");
    FillHistogram("hCentrality1V0C",v0C) ;

    Double_t v0M =centrality->GetCentralityPercentile("V0M");
    FillHistogram("hCentrality1V0M",v0M) ;
    
    Double_t cl1 =centrality->GetCentralityPercentile("CL1");
    FillHistogram("hCentrality1CL1",cl1) ;
      
    Double_t zna =centrality->GetCentralityPercentile("ZNA");
    FillHistogram("hCentrality1ZNA",zna) ;

    Double_t znc =centrality->GetCentralityPercentile("ZNC");
    FillHistogram("hCentrality1ZNC",znc) ;
    
//       if(fRunNumber >=195344 && fRunNumber <= 197388){ //LHC13bcdef  
    if(!fUtils->IsVertexSelected2013pA(event)){
      return kFALSE ;
    }
//       }
    
    FillHistogram("hSelEvents",4) ;
  
    if(fUtils->IsPileUpEvent(event)){
          return kFALSE ;
    }
    FillHistogram("hSelEvents",5) ;
    
    
    //Centrality distribution after event selection
    v0A =centrality->GetCentralityPercentile("V0A");
    FillHistogram("hCentrality2V0A",v0A) ;
      
    v0C =centrality->GetCentralityPercentile("V0C");
    FillHistogram("hCentrality2V0C",v0C) ;

    v0M =centrality->GetCentralityPercentile("V0M");
    FillHistogram("hCentrality2V0M",v0M) ;
    
    cl1 =centrality->GetCentralityPercentile("CL1");
    FillHistogram("hCentrality2CL1",cl1) ;
      
    zna =centrality->GetCentralityPercentile("ZNA");
    FillHistogram("hCentrality2ZNA",zna) ;

    znc =centrality->GetCentralityPercentile("ZNC");
    FillHistogram("hCentrality2ZNC",znc) ;
    
    //centrality
    fCentrality=1.;
    if(fRunNumber >=195344 && fRunNumber <= 196311){ //LHC13bcde  
        switch(fCentEstimator){
            case 4 : fCentrality=centrality->GetCentralityPercentile("CL1");
                     break;
            case 3 : fCentrality=centrality->GetCentralityPercentile("ZNA");
                     break;         
            case 2 : fCentrality=centrality->GetCentralityPercentile("V0M");
                     break;
            case 1 : 
            default: fCentrality=centrality->GetCentralityPercentile("V0A");
        }
     }
     else if(fRunNumber >=196433 && fRunNumber <= 197388 ){ //LHC13f  
        switch(fCentEstimator){
            case 4 : fCentrality=centrality->GetCentralityPercentile("CL1");
                     break;
            case 3 : fCentrality=centrality->GetCentralityPercentile("ZNC");
                     break;         
            case 2 : fCentrality=centrality->GetCentralityPercentile("V0M");
                     break;
            case 1 : 
            default: fCentrality=centrality->GetCentralityPercentile("V0C");
        }
     }

     fCentBin=0;
     while(fCentBin<fNCenBin && fCentrality>fCenBinEdges.At(fCentBin))
        fCentBin++ ;
     if(fCentBin>=fNCenBin) fCentBin=fNCenBin-1; 
      
     if(fIsMB)
        fCentWeight=MBCentralityWeight(fCentrality); 
     else
        fCentWeight=TrigCentralityWeight(fCentrality); 
     return kTRUE ;
  } //end of p-Pb case
    
  
  
  //In case of p-Pb data Run2  
  if(fRunNumber >=265015 	 && fRunNumber <= 267166 ){ //LHC16qrst
    
    //Fill Centrality before vertex/pileup cuts
    AliMultSelection *multSelection = (AliMultSelection*) event -> FindListObject("MultSelection");
    if(!multSelection)
      return kFALSE ;  

    Double_t v0A =multSelection->GetMultiplicityPercentile("V0A");
    FillHistogram("hCentrality1V0A",v0A) ;
      
    Double_t v0C =multSelection->GetMultiplicityPercentile("V0C");
    FillHistogram("hCentrality1V0C",v0C) ;

    Double_t v0M =multSelection->GetMultiplicityPercentile("V0M");
    FillHistogram("hCentrality1V0M",v0M) ;
    
    Double_t cl1 =multSelection->GetMultiplicityPercentile("CL1");
    FillHistogram("hCentrality1CL1",cl1) ;
      
    Double_t zna =multSelection->GetMultiplicityPercentile("ZNA");
    FillHistogram("hCentrality1ZNA",zna) ;

    Double_t znc =multSelection->GetMultiplicityPercentile("ZNC");
    FillHistogram("hCentrality1ZNC",znc) ;
           
    //multSelection
    fCentrality=1.;
    switch(fCentEstimator){
      case 4 : fCentrality=cl1;
               break;
      case 3 : fCentrality=zna;
               break;         
      case 2 : fCentrality=v0M;
               break;
      case 1 : 
      default: fCentrality=v0A;
    }

     fCentBin=0;
     while(fCentBin<fNCenBin && fCentrality>fCenBinEdges.At(fCentBin))
        fCentBin++ ;
     if(fCentBin>=fNCenBin) fCentBin=fNCenBin-1; 
           
     if(fIsMB)
        fCentWeight=MBCentralityWeight(fCentrality); 
     else
        fCentWeight=TrigCentralityWeight(fCentrality); 
     return kTRUE ;
  } //end of p-Pb Run2 case
    
  
  
  //pp collisions
  //Fill Multiplicity before vertex/pileup cuts
  AliMultSelection *multSelection = (AliMultSelection*) event -> FindListObject("MultSelection");
    
  if (!multSelection ){
    return kFALSE ;  
  }
  else{
      

    Double_t v0M = multSelection->GetMultiplicityPercentile("V0M");
    FillHistogram("hCentrality1V0M",v0M) ;
    if(v0M==199. && !fIsMB){
      v0M=0.;  // No calibration
    }
    
    Double_t v0A =multSelection->GetMultiplicityPercentile("V0A");
    FillHistogram("hCentrality1V0A",v0A) ;
      
    Double_t v0C =multSelection->GetMultiplicityPercentile("V0C");
    FillHistogram("hCentrality1V0C",v0C) ;

    Double_t olv0M =multSelection->GetMultiplicityPercentile("OnlineV0M");
    FillHistogram("hCentrality1OnlineV0M",olv0M) ;
    
    Double_t olv0A =multSelection->GetMultiplicityPercentile("OnlineV0A");
    FillHistogram("hCentrality1OnlineV0A",olv0A) ;
    
    Double_t olv0C =multSelection->GetMultiplicityPercentile("OnlineV0C");
    FillHistogram("hCentrality1OnlineV0C",olv0C) ;
    
    Double_t adM = multSelection->GetMultiplicityPercentile("ADM");
    FillHistogram("hCentrality1ADM",adM) ;

    Double_t adA =multSelection->GetMultiplicityPercentile("ADA");
    FillHistogram("hCentrality1ADA",adA) ;
      
    Double_t adC =multSelection->GetMultiplicityPercentile("ADC");
    FillHistogram("hCentrality1ADC",adC) ;

    Double_t cl1 =multSelection->GetMultiplicityPercentile("SPDClusters");
    FillHistogram("hCentrality1SPDClusters",cl1) ;
      
    Double_t spd =multSelection->GetMultiplicityPercentile("SPDTracklets");
    FillHistogram("hCentrality1SPDTracklets",spd) ;

    Double_t refMu05 =multSelection->GetMultiplicityPercentile("RefMult05");
    FillHistogram("hCentrality1RefMult05",refMu05) ;
    
    Double_t refMu08 =multSelection->GetMultiplicityPercentile("RefMult08");
    FillHistogram("hCentrality1RefMult08",refMu08) ;
    
    Double_t sphericity = CalculateSphericity() ;
    
    Double_t spherocity = CalculateSpherocity() ;
    FillHistogram("hSpheriSpheroMult",sphericity,spherocity,v0M) ;
    
        
    FillHistogram("hSelEvents",4) ;
  
    if(fUtils->IsPileUpEvent(event)){
          return kFALSE ;
    }
    FillHistogram("hSelEvents",5) ;
    
    FillHistogram("hCentrality2V0M",v0M) ;
    FillHistogram("hCentrality2V0A",v0A) ;
    FillHistogram("hCentrality2V0C",v0C) ;
    FillHistogram("hCentrality2OnlineV0M",olv0M) ;
    FillHistogram("hCentrality2OnlineV0A",olv0A) ;
    FillHistogram("hCentrality2OnlineV0C",olv0C) ;
    FillHistogram("hCentrality2ADM",adM) ;
    FillHistogram("hCentrality2ADA",adA) ;
    FillHistogram("hCentrality2ADC",adC) ;
    FillHistogram("hCentrality2SPDClusters",cl1) ;
    FillHistogram("hCentrality2SPDTracklets",spd) ;
    FillHistogram("hCentrality2RefMult05",refMu05) ;
    FillHistogram("hCentrality2RefMult08",refMu08) ;
    
    
    //centrality
    fCentrality=1.;
    switch(fCentEstimator){
            case 15 : fCentrality=100*sphericity ; //
                     break;
            case 14 : fCentrality=100*spherocity ; //
                     break;
            case 13 : fCentrality=v0M;
                     break;
            case 12 : fCentrality=v0A;
                     break;         
            case 11 : fCentrality=v0C;
                     break;         
            case 10 : fCentrality=olv0M;
                     break;         
            case 9 : fCentrality=olv0A;
                     break;         
            case 8 : fCentrality=olv0C;
                     break;         
            case 7 : fCentrality=multSelection->GetMultiplicityPercentile("ADM");
                     break;         
            case 6 : fCentrality=multSelection->GetMultiplicityPercentile("ADA");
                     break;         
            case 5 : fCentrality=multSelection->GetMultiplicityPercentile("ADC");
                     break;         
            case 4 : fCentrality=multSelection->GetMultiplicityPercentile("SPDClusters");
                     break;         
            case 3 : fCentrality=multSelection->GetMultiplicityPercentile("SPDTracklets");
                     break;         
            case 2 : fCentrality=multSelection->GetMultiplicityPercentile("RefMult05");
                     break;         
            case 1 : 
            default: fCentrality=multSelection->GetMultiplicityPercentile("RefMult08");
                     break;         
      }
    }
    
     fCentBin=0;
     if(fCentrality<0){
       FillHistogram("hSelEvents",11) ;
       return kFALSE ;  
     }
     
     if(fCentrality>100){
       FillHistogram("hSelEvents",12) ;
       return kFALSE ;  
     }
     
     while(fCentBin<fNCenBin && fCentrality>fCenBinEdges.At(fCentBin))
        fCentBin++ ;
     if(fCentBin>=fNCenBin) fCentBin=fNCenBin-1; 
      
     if(fIsMB)
        fCentWeight=1. ;// CentralityWeightPP(fCentrality); 
     else
        fCentWeight=1. ; //TrigCentralityWeightPP(fCentrality); 
     return kTRUE ;

}
//___________________________________________________________________________
Double_t AliAnalysisTaskTaggedPhotons::TOFCutEff(Double_t x ){
  //Time cut efficiency in MC  
  //time cut efficiency in 15n+17pq, update 02.10.2019
  //12.5 ns cut
  if(TMath::Abs(fTimeCut-12.5e-9)<0.01*12.5e-9){  
     if(x<4.6){
       return TMath::Exp(3.64952e-03*(-5.80032e+01-1.53442e+02*x+1.30994e+02*x*x+-3.53094e+01*x*x*x+x*x*x*x)/
                        (-7.75638e-02+8.64761e-01*x+1.22320e+00*x*x-1.00177e+00*x*x*x+x*x*x*x)) ;
     }
     else{
       return 0.63922783*(1.-1.63273e-01*TMath::TanH((x-7.94528e+00)/1.28997e+00))*
                         (-4.39257e+00*x+2.25503e+00*x*x+x*x*x)/(2.37160e+00*x-6.93786e-01*x*x+x*x*x) ;
     }
  }
  if(TMath::Abs(fTimeCut-30.e-9)<0.01*30.e-9){ 
    if(fDefTof){  
      //Improved param for E<0.2 GeV (07.07.2020)  
      if(x>1.6)x=1.6;
      if(x<0.14)x=0.14;
      return TMath::Exp((6.24104e+05-6.40577e+06*x+1.25640e+07*x*x-9.59211e+06*x*x*x+2.53582e+06*x*x*x*x)/(1.-2.34722e+06*x+1.68667e+07*x*x)) ;
    }
    else{ //alternatice param        
      //another param 11.03.2021
      if(x>1.6)x=1.6;
      return 1./(1.+TMath::Exp((4.83230e+01-8.89758e+01*x+1.10897e+03*x*x-5.73755e+03*x*x*x-1.43777e+03*x*x*x*x)/
            (1.-1.23667e+02*x+1.07255e+03*x*x+5.87221e+02*x*x*x))) ;
    }
  }
  if(TMath::Abs(fTimeCut-100.e-9)<0.01*100.e-9){
    //parameterization 01.08.2020   
    if(x>1.1)x=1.1;
    if(x<0.11)x=0.11;
    return TMath::Exp((-1.15295e+05+2.26754e+05*x-1.26063e+05*x*x+x*x*x)/(1.-3.16443e+05*x+3.68044e+06*x*x+x*x*x)) ;
  }
  //no other parameterizations so far
  return 1.; 
}
//___________________________________________________________________________
double AliAnalysisTaskTaggedPhotons::CalculateSphericity(){
  TMatrixD Si(2,2);
  TMatrixD S(2,2);
  Double_t p(0);

  Int_t numberOfRecTracks = 0;
  Int_t nTr = fInputEvent->GetNumberOfTracks() ;
  for(Int_t iTracks = 0; iTracks<nTr; iTracks++){
      AliAODTrack* curTrack = (AliAODTrack*) fInputEvent->GetTrack(iTracks);
      if(fIsMB){
        if(curTrack->GetID()<0) continue; // Avoid double counting of tracks
        if(!curTrack->IsHybridGlobalConstrainedGlobal()) continue;
      }
      else{
        if(!curTrack->IsOn(AliVTrack::kITSpureSA)) continue ;
        if(!curTrack->TestFilterBit(AliAODTrack::kTrkITSsa)) continue;
        if(curTrack->GetITSNcls()<4) continue ;
      }
      if(TMath::Abs(curTrack->Eta())>0.8) continue;
      if(curTrack->Pt()<0.5) continue;
      numberOfRecTracks++;

      Si(0,0)=pow(curTrack->Px(),2);
      Si(0,1)=curTrack->Px()*curTrack->Py();
      Si(1,0)=curTrack->Py()*curTrack->Px();
      Si(1,1)=pow(curTrack->Py(),2);
      S += (1./curTrack->Pt())*Si;
      p += curTrack->Pt();
  }
  if(p==0 || numberOfRecTracks<3){
      return -0.1;
  }else{
      TMatrixDEigen eg((1./p)*S); 
      return (2*eg.GetEigenValues()(1,1))/(eg.GetEigenValues()(0,0)+eg.GetEigenValues()(1,1));
  }
}
//________________________________________________________________________
Double_t AliAnalysisTaskTaggedPhotons::CalculateSpherocity(){
  const int minMulti=3;
  const double minPt=0.15;
  Int_t nTr = fInputEvent->GetNumberOfTracks() ;
  if(nTr<minMulti) return -0.1;//if not enough tracks (1st check), return -0.1
  int trackMulti=0;
  
  Double_t sumpt=0;
  std::vector<Double_t> nx;
  std::vector<Double_t> ny;
  std::vector<Double_t> px;
  std::vector<Double_t> py;
  for(Int_t itrk=0;itrk<nTr;++itrk) {
    AliAODTrack *track=(AliAODTrack*)fInputEvent->GetTrack(itrk);
    if(!track) continue;

    if(fIsMB){
      if(track->GetID()<0) continue; // Avoid double counting of tracks
        if(!track->IsHybridGlobalConstrainedGlobal()) continue;
    }
    else{
      if(!track->IsOn(AliVTrack::kITSpureSA)) continue ;
      if(track->GetITSNcls()<4) continue ;
    }
    if(TMath::Abs(track->Eta())>0.8) continue;
//       if(!((AliAODTrack*)track)->TestFilterBit(fAODFilterBit)) continue; //Check filter bit
//       //Manually check for ITS & TPC refit
//       if((track->GetStatus()&AliVTrack::kITSrefit)!=AliVTrack::kITSrefit) continue;
//       if((track->GetStatus()&AliVTrack::kTPCrefit)!=AliVTrack::kTPCrefit) continue;
  
    Double_t pt = track->Pt();
    if(pt<minPt) continue;
    sumpt+=pt;
    Double_t phi = track->Phi();
    nx.push_back(TMath::Cos(phi));
    ny.push_back(TMath::Sin(phi));
    px.push_back(pt*nx[trackMulti]);
    py.push_back(pt*ny[trackMulti]);
    ++trackMulti;
  };
  if(trackMulti<minMulti) {
    return -0.1; //if not enought tracks
  };
  //Calculating spherocity now
  Double_t retval=2; //Return value should be the minimal value in range 0..1, so 2 is a good starting point.
  for(Int_t i=0; i<trackMulti; i++) {
    Double_t num=0;
    for(Int_t j=0;j<trackMulti;j++)
      num+=TMath::Abs(ny[i]*px[j] -nx[i]*py[j]);
    Double_t sFull = TMath::Power((num/sumpt),2);
    if(sFull < retval)
      retval = sFull;
  }
  if(retval>1) return -0.2; //If sph>1, something went wrong
  retval=retval*TMath::Pi()*TMath::Pi()/4; //normalization
  return retval;

}
Bool_t AliAnalysisTaskTaggedPhotons::AcceptJJevent(){
  //Accept JJ event only if Jet pt is smaller than fJetPtHardFactor * ptHardBin
  //to avoid fluctuations.
    
    // Pythia header
    TString pyGenName       = ""; 
    TString pyProcessName   = "";  
    Int_t   pyProcess       = 0;
    Int_t   pyFirstGenPart  = 0; 
    Int_t   pythiaVersion   = 0;
    
    // Init it first to 0 to tell the method to recover it.
    AliGenPythiaEventHeader * genPythiaEventHeader = 
    fMCAnalysisUtils->GetPythiaEventHeader(MCEvent(),fMCGenerEventHeaderToAccept,
                                           pyGenName,pyProcessName,pyProcess,pyFirstGenPart,pythiaVersion);    
    if ( genPythiaEventHeader )
  {      
    Int_t nTriggerJets =  genPythiaEventHeader->NTriggerJets();
    Float_t ptHard = genPythiaEventHeader->GetPtHard();
    TParticle * jet =  0;    
    Float_t tmpjet[]={0,0,0,0};
    for(Int_t ijet = 0; ijet< nTriggerJets; ijet++)
    {
      genPythiaEventHeader->TriggerJet(ijet, tmpjet);
      jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
      
      //Compare jet pT and pt Hard
      if(jet->Pt() > fJetPtHardFactor * ptHard)
      {
        return kFALSE;
      }
    } // jet loop
    
    if(jet) delete jet;
  } // pythia header
  
  return kTRUE ;
}
Double_t AliAnalysisTaskTaggedPhotons::Unsmear(Double_t e, Double_t eMC){

  // Too big difference btw measured and MC energy: wrong MC? do not correct e  
  if(TMath::Abs(e-eMC)>0.15*e) return e ;
    
  //Correct for NonLinearity
  Double_t corr =1.00356*(1.-3.17857e-01*exp(-eMC*5.51728))*(1-0.0973089*exp(-eMC*1.58757))*(1.+2.63610e-01*exp(-eMC*2.63638)) ;
  eMC/=corr;  
  Double_t scale = 1. - fUnsmA*TMath::Exp(-e*fUnsmB) ;
  return (eMC + (e-eMC)*scale) ; 
    
}
Double_t AliAnalysisTaskTaggedPhotons::DPMJETKCorr(Double_t x){
  //kaon yield correction for DPMJET vs pT
  return (1.23546+0.613098*exp(-x/4.81805))*(1.-0.796292*exp(-x/0.660332));
}
Double_t AliAnalysisTaskTaggedPhotons::DPMJETpCorr(Double_t x){
  //baryon yield correction for DPMJET vs pT
  return (1.18027+16.6995*exp(-x/2.42366))*(1.-1.00768*exp(-x/9.11308));
}


