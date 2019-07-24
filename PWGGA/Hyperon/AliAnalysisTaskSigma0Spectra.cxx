/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author:  Alexander Borissov (PNU)   Dmitri Peressounko (RRC KI)        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 * updated for Sigma^0 analysis Alexander Borissov, Jihye Song, Dec. 2014 *
 **************************************************************************/
////////////////////////////////////////////////
//---------------------------------------------
// Class used to do analysis on conversion + calorimeter pairs
// A lot of code cut-and-pasted from AliV0Reader and GammaConversion classes
//
//---------------------------------------------
////////////////////////////////////////////////

// root
#include "TChain.h"
#include "TH3.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "TRandom.h"


// analysis
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskSigma0Spectra.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrackCuts.h"
#include "AliESDCaloCluster.h"
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliESDv0.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliCaloParticle.h"

class Riostream;
class TFile;

ClassImp(AliAnalysisTaskSigma0Spectra)


//______________________________________________________________________________
AliAnalysisTaskSigma0Spectra::AliAnalysisTaskSigma0Spectra():
AliAnalysisTaskSigma0(),
fPHOSPi0Event(0x0),
fEMCALPi0Event(0x0),
fConvPi0Event(0x0)
{
    Int_t nBin=10 ;
    for(Int_t i=0;i<nBin;i++){
        fPHOSPi0Events[i]=0 ;
    }
    for(Int_t i=0;i<nBin;i++){
        fEMCALPi0Events[i]=0 ;
    }
    for(Int_t i=0;i<nBin;i++){
        fConvPi0Events[i]=0 ;
    }
}


//______________________________________________________________________________
AliAnalysisTaskSigma0Spectra::AliAnalysisTaskSigma0Spectra(const char* name):
AliAnalysisTaskSigma0(name),
fPHOSPi0Event(0x0),
fEMCALPi0Event(0x0),
fConvPi0Event(0x0)
{
    Int_t nBin=10 ;
    for(Int_t i=0;i<nBin;i++){
        fPHOSPi0Events[i]=new TList()  ;
    }
    
    for(Int_t i=0;i<nBin;i++){
        fEMCALPi0Events[i]=new TList()  ;
    }
    
    for(Int_t i=0;i<nBin;i++){
        fConvPi0Events[i]=new TList()  ;
    }
}
//_____________________________________________________
AliAnalysisTaskSigma0Spectra::~AliAnalysisTaskSigma0Spectra()
{
    // Remove all pointers
    if(fPHOSPi0Event){
        delete fPHOSPi0Event ;
        fPHOSPi0Event=0x0 ;
    }
    Int_t nBin=10 ;
    for(Int_t i=0;i<nBin;i++){
        if(fPHOSPi0Events[i]){
            delete fPHOSPi0Events[i] ;
            fPHOSPi0Events[i]=0x0 ;
        }
    }
    
    if(fEMCALPi0Event){
        delete fEMCALPi0Event ;
        fEMCALPi0Event=0x0 ;
    }
    for(Int_t i=0;i<nBin;i++){
        if(fEMCALPi0Events[i]){
            delete fEMCALPi0Events[i] ;
            fEMCALPi0Events[i]=0x0 ;
        }
    }
    
    if(fConvPi0Event){
        delete fConvPi0Event ;
        fConvPi0Event=0x0 ;
    }
    for(Int_t i=0;i<nBin;i++){
        if(fConvPi0Events[i]){
            delete fConvPi0Events[i] ;
            fConvPi0Events[i]=0x0 ;
        }
    }
    
}
//_____________________________________________________
void AliAnalysisTaskSigma0Spectra::Init()
{
    AliAnalysisTaskSigma0::Init() ;
}
//____________________________________________________________
void AliAnalysisTaskSigma0Spectra::UserCreateOutputObjects()
{
    AliAnalysisTaskSigma0::UserCreateOutputObjects();
    
    //Single spectra
    // Int_t nchxt = 500 ; unused variable
    Int_t npt=150 ;   //was 500 for pi0 spectra
    Double_t ptmin = 0. ;  Double_t ptmax=15. ;   // was 50.
    //  Int_t nchMgg = 400 ;     // Double_t mggmax=0.7;
    //    Double_t mSigmax=1.4;  Double_t mSigmin=1.1;
    // Double_t mLammax=1.20;  Double_t mLammin=1.05;
    //    Int_t nchMgg = 280;
    
    Double_t mSigmax=1.230;  Double_t mSigmin=1.130;                                                                                                              
    Double_t mLammax=1.122;  Double_t mLammin=1.110;                                                                                                               
    Int_t nchMgg = 100;                                                                                                                                       

    //    fOutputContainer->Add(new TH1F("hPhosGamPt","PHOS gamma  Pt", npt*2, ptmin,ptmax )) ;
    fOutputContainer->Add(new TH1F("hCaloGamPt","Conv gamma Pt",  npt*2, ptmin,ptmax )) ;
    fOutputContainer->Add(new TH2F("hV0LamMvsPt","V0 Lam ALL Mass vs pt",120, mLammin,mLammax, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hV0PLamMvsPt","V0 Lam Part Mass vs pt",120, mLammin,mLammax, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hV0ALamMvsPt","V0 Lam APart Mass vs pt",120, mLammin,mLammax, 100,0.,ptmax)) ;
    
    
    //    fOutputContainer->Add(new TH2F("hCCReMvsPt1","CC Mass vs ArPod1",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    // fOutputContainer->Add(new TH2F("hCCReMvsPtnArPod0","CC Mass vs ptnArPod0",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    // fOutputContainer->Add(new TH2F("hCCReMvsPtnArPod","CC Mass vs ptnArPod",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    fOutputContainer->Add(new TH2F("hCCReMvsPt","CC Mass vs pt",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    fOutputContainer->Add(new TH2F("hCCRenegMvsPt","CC Mass vs pt",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    fOutputContainer->Add(new TH2F("hCCReposMvsPt","CC Mass vs pt",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    
    
    // fOutputContainer->Add(new TH2F("hCCReMvsPt0","PP Mass vs pt0",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    //   fOutputContainer->Add(new TH2F("hCCRenegMvsPt0","PP Mass vs pt0",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    //  fOutputContainer->Add(new TH2F("hCCReposMvsPt0","PP Mass vs pt0",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    
    // fOutputContainer->Add(new TH2F("hCCMiMvsPt","PP Mass vs pt Mix",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    // fOutputContainer->Add(new TH2F("hCCMinegMvsPt","PP Mass vs pt Mix",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    // fOutputContainer->Add(new TH2F("hCCMiposMvsPt","PP Mass vs pt Mix",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    
    //    fOutputContainer->Add(new TH2F("hCCMiMvsPt0","PP Mass vs pt Mix0",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    //fOutputContainer->Add(new TH2F("hCCMinegMvsPt0","PP Mass vs pt Mix0",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
    // fOutputContainer->Add(new TH2F("hCCMiposMvsPt0","PP Mass vs pt Mix0",nchMgg, mSigmin,mSigmax ,npt,0.,ptmax));
        
    //  fOutputContainer->Add(new TH1F("hSig0PhosRPhiEta","Sig0 P Rapidity", npt, -3, 3 )) ;
    // fOutputContainer->Add(new TH1F("hSig0CaloRPhiEta","Sig0 C Rapidity",  npt, -3, 3 )) ;
    
    //  fOutputContainer->Add(new TH1F("hSig0PhosRapALL","Sig0 P Rapidity", npt, -3, 3 )) ;
    // fOutputContainer->Add(new TH1F("hSig0CaloRapALL","Sig0 C Rapidity",  npt, -3, 3 )) ;
    
    // fOutputContainer->Add(new TH2F("hPhosGamPhiEta","PHOS gamma  PhiEta", 32, -3.2, 3.2, 20, -2, 2  )) ;
    fOutputContainer->Add(new TH2F("hCaloGamPhiEta","Conv gamma PhiETA",  32, -3.2, 3.2, 20, -2, 2  )) ;
    
    //    fOutputContainer->Add(new TH2F("hSig0CaloRapvsPt","C Rap vs pt",100, -2., 2., npt,0.,ptmax));
    //   fOutputContainer->Add(new TH2F("hSig0PhosRapvsPt","P Rap vs pt",100, -2., 2., npt,0.,ptmax));
    
    // fOutputContainer->Add(new TH2F("ArPodSig0P","Armenteros-Podolanski Sig0P ;#alpha;p_{t}",100,-1.0,1.0, 80,0,0.8) );
    // fOutputContainer->Add(new TH2F("ArPodSig0C","Armenteros-Podolanski Sig0C ;#alpha;p_{t}",100,-1.0,1.0, 80,0,0.8) );
    //fOutputContainer->Add(new TH2F("ArPodSig0P2","Armenteros-Podolanski2 Sig0P ;#alpha;p_{t}",100,-1.0,1.0, 80,0,0.8) );
    // fOutputContainer->Add(new TH2F("ArPodSig0C2","Armenteros-Podolanski2 Sig0C ;#alpha;p_{t}",100,-1.0,1.0, 80,0,0.8) );
    
    
    char key[55] ;
    for(Int_t is=0;   is<=5 ; is++){        
        sprintf(key,"hDtPSig0MvsPtRec%d",is) ;
        fOutputContainer->Add(new TH2F(key,"pT reconstructed of PSigma0",nchMgg,mSigmin,mSigmax, npt,0.,ptmax)) ;
        sprintf(key,"hDtASig0MvsPtRec%d",is) ;
        fOutputContainer->Add(new TH2F(key,"pT reconstructed of ASigma0",nchMgg,mSigmin,mSigmax, npt,0.,ptmax)) ;
        
        sprintf(key,"hDtPSig0MvsPtMix%d",is) ;
        fOutputContainer->Add(new TH2F(key,"mix BG pT reconstructed of PSigma0",nchMgg,mSigmin,mSigmax, npt,0.,ptmax)) ;
        sprintf(key,"hDtASig0MvsPtMix%d",is) ;
        fOutputContainer->Add(new TH2F(key,"mix BG pT reconstructed of ASigma0",nchMgg,mSigmin,mSigmax, npt,0.,ptmax)) ;
        

        sprintf(key,"hDtPLamMvsPtRec%d",is) ;
        fOutputContainer->Add(new TH2F(key,"pT reconstructed of PLam",nchMgg,mLammin,mLammax, npt,0.,ptmax)) ;
        sprintf(key,"hDtALamMvsPtRec%d",is) ;
        fOutputContainer->Add(new TH2F(key,"pT reconstructed of ALam",nchMgg,mLammin,mLammax, npt,0.,ptmax)) ;        
    }

    for(Int_t is=0;   is<=5 ; is++){        
        sprintf(key,"hPHOS_Sig0MvsPtRec%d",is) ;
        fOutputContainer->Add(new TH2F(key,"pT reconstructed of Sigma0",nchMgg,mSigmin,mSigmax, npt,0.,ptmax)) ;
        sprintf(key,"hcPHOS_Sig0MvsPtRec%d",is) ;
        fOutputContainer->Add(new TH2F(key,"pTc reconstructed of Sigma0",nchMgg,mSigmin,mSigmax, npt,0.,ptmax)) ;

        sprintf(key,"hPHOS_PSig0MvsPtRec%d",is) ;
        fOutputContainer->Add(new TH2F(key,"pT reconstructed of PSigma0",nchMgg,mSigmin,mSigmax, npt,0.,ptmax)) ;
        sprintf(key,"hPHOS_ASig0MvsPtRec%d",is) ;
        fOutputContainer->Add(new TH2F(key,"pT reconstructed of ASigma0",nchMgg,mSigmin,mSigmax, npt,0.,ptmax)) ;
        
	/*	sprintf(key,"hPHOS_Sig0MvsPtMix%d",is) ;
        fOutputContainer->Add(new TH2F(key,"mix BG pT reconstructed of Sigma0",nchMgg,mSigmin,mSigmax, npt,0.,ptmax)) ;
	sprintf(key,"hPHOS_PSig0MvsPtMix%d",is) ;
        fOutputContainer->Add(new TH2F(key,"mix BG pT reconstructed of PSigma0",nchMgg,mSigmin,mSigmax, npt,0.,ptmax)) ;
        sprintf(key,"hPHOS_ASig0MvsPtMix%d",is) ;
        fOutputContainer->Add(new TH2F(key,"mix BG pT reconstructed of ASigma0",nchMgg,mSigmin,mSigmax, npt,0.,ptmax)) ;
	*/
  }
     
    PostData(1, fOutputContainer);
}


//_____________________________________________________
void AliAnalysisTaskSigma0Spectra::UserExec(Option_t *option)
{
    
    //   printf("Sigma0Spectra Execute analysis for current event \n");
    const Int_t nEventsToKeep=100 ;
    if(fPHOSPi0Event || fEMCALPi0Event  ){
        if(   fPHOSPi0Event->GetEntriesFast()>0
           || fEMCALPi0Event->GetEntriesFast()>0
           || fConvPi0Event->GetEntriesFast()>0   )   {
            //valirable with Vertex position should not change since last event yet...
            //If it is not defined so far, there should not be photons
            Int_t izvtx = (Int_t)((fZvtx+10.)/2.) ;
            if(izvtx<0)izvtx=0 ;
            if(izvtx>9)izvtx=9 ;
            
            TList * prevConv = fConvPi0Events[izvtx] ;

           prevConv->AddFirst(fConvPi0Event) ;
            fConvPi0Event=0;
            
            if(prevConv->GetSize()>nEventsToKeep){//Remove redundant events
                TClonesArray * tmp = static_cast<TClonesArray*>(prevConv->Last()) ;
                prevConv->RemoveLast() ;
                delete tmp ;
            }
        }
    }
    if(fPHOSPi0Event)
    fPHOSPi0Event->Clear() ;
    else
    fPHOSPi0Event = new TClonesArray("AliCaloParticle",10) ;
    //  FillHistogram("hQA_EventsTrig",5.5) ;
    
    if(fEMCALPi0Event)
    fEMCALPi0Event->Clear() ;
    else
    fEMCALPi0Event = new TClonesArray("AliCaloParticle",10) ;
    if(fConvPi0Event)
    fConvPi0Event->Clear() ;
    else
    fConvPi0Event = new TClonesArray("AliCaloParticle",10) ;
    //  FillHistogram("hQA_EventsTrig",6.5) ;
    
    AliAnalysisTaskSigma0::UserExec(option) ;
    
    FillRealMixed() ;
    PostData(1, fOutputContainer);
    
}
//_____________________________________________________
void AliAnalysisTaskSigma0Spectra::FillSpectrum(){ }

//_____________________________________________________
void AliAnalysisTaskSigma0Spectra::FillRealMixed(){
    
    Int_t izvtx = (Int_t)((fZvtx+10.)/2.) ;
    if(izvtx<0)izvtx=0 ;
    if(izvtx>9)izvtx=9 ;
    //    TList * prevPHOS = fPHOSEvents[izvtx] ; 
    TList * prevConv = fConvEvents[izvtx] ;

    //    Double_t fCentCuts[5] = {0, 20, 40, 60, 100 }; 
    
    //    if( fPHOSEvent->GetEntriesFast() )  printf("PHOS clusters = %d \n",fPHOSEvent->GetEntriesFast()) ;
    char key[55] ;
    // Fill Real distributions  PHOS + PHOS
    Int_t nConv=fConvEvent->GetEntriesFast() ;
    Int_t nPHOS=fPHOSEvent->GetEntriesFast() ; 
    //  Int_t nEMCAL=fEMCALEvent->GetEntriesFast() ;   // can use EMCAL later
    Int_t nLam=fGenpi0Event->GetEntriesFast() ;
    
    Double_t MassLam = 1.115683 ;
    //  Double_t Rmin[5]={ 1.4, 1.2, 1.0, 0.8, 0.6, 0.5 };     // Double_t dRphieta = 0; unused variable
 

  //========================================================================================================== 
  // PHOS
  for(Int_t iPHOS=0; iPHOS<nPHOS; iPHOS++){
    AliCaloParticle * ph1 = static_cast<AliCaloParticle*>(fPHOSEvent->At(iPHOS)) ;

    //    FillHistogram("P_gamPt", ph1->Pt() ) ;
    //    FillHistogram("P_gamPhiEta", ph1->Phi(), ph1->Eta() ) ;

    for(Int_t iConv2=0; iConv2<nLam; iConv2++){
      AliCaloParticle * ph2c = static_cast<AliCaloParticle*>(fGenpi0Event->At(iConv2)) ;     
      
      //      if (  ph2c->M() < 1 ) continue;   // consider only (bar)Lambda from conversions
      //      printf("M-PC  M %f PT %f  E %f nconv %d \n", ph2c->M(), ph2c->Pt(), ph2c->E(), nConv );
 
      Double_t ModLam =  ph2c->GetXt(); 
      //      printf("-----PHOS ModLam %f MassLam %f evt %d \n",ModLam,ph2c->M(),fESDEvent->GetEventNumberInFile());

      Double_t Mlamrec = ph2c->M();		
      Double_t corr_factor = MassLam / Mlamrec;		
      TLorentzVector ph2corr;
      ph2corr = corr_factor * *ph2c;

      //      dRphieta=sqrt( (ph1->Phi()-ph2c->Phi())*(ph1->Phi()-ph2c->Phi()) + (ph1->Eta()-ph2c->Eta())*(ph1->Eta()-ph2c->Eta()) );
      TLorentzVector pi0p=*ph1 + ph2corr ; 
      Double_t mphos=pi0p.M() ;
      Double_t ptphos=pi0p.Pt() ;
      Double_t pos[3]= { ph1->Px(),  ph1->Py(),  ph1->Pz() };
      Double_t neg[3]= { ph2corr.Px(),  ph2corr.Py(),  ph2corr.Pz() }; 
      Double_t moth[3]= { pi0p.Px(),  pi0p.Py(),  pi0p.Pz() };
      Double_t arpod[2]= {0,0};

      GetArPod( pos, neg, moth, arpod );
      if( ( abs(  arpod[1] ) > 0.6 && arpod[0]< 0.12 ) ) {
	Double_t rapphos = Rapidity(  ptphos,  pi0p.Pz(),  mphos );
     
	for (Int_t irap = 0; irap < 6; irap++) {
	  if( fabs(rapphos) < fEtaCuts[irap] ) {
	    sprintf(key,"hcPHOS_Sig0MvsPtRec%d",irap) ;
	    FillHistogram(key, mphos, ptphos );
	    
	    if( ModLam > 0 ) {
	      sprintf(key,"hPHOS_PSig0MvsPtRec%d",irap) ;
	      FillHistogram(key, mphos, ptphos );
	    }
	    else  if( ModLam < 0 ) {
	      sprintf(key,"hPHOS_ASig0MvsPtRec%d",irap) ;
	      FillHistogram(key, mphos, ptphos );
	    }
	  }
	}
      }

      TLorentzVector sig0prec=*ph1 + *ph2c ; 
      Double_t mrec=sig0prec.M() ;
      Double_t ptrec=sig0prec.Pt() ;
      Double_t posr[3]= { ph1->Px(),  ph1->Py(),  ph1->Pz() };
      Double_t negr[3]= { ph2c->Px(),  ph2c->Py(),  ph2c->Pz() }; 
      Double_t mothr[3]= { sig0prec.Px(),  sig0prec.Py(), sig0prec.Pz() };
      Double_t arpodr[2]= {0,0};
      GetArPod( posr, negr, mothr, arpodr );
      if( ( abs(  arpodr[1] ) > 0.6 && arpodr[0]< 0.12 ) ) {
	/* FillHistogram("PC_Re_mvsPt_0", mrec, ptrec ) ;
	if( ModLam < 0 ) 	FillHistogram("PC_Reneg_mvsPt_0",mrec,ptrec ) ;
	else if(   ModLam > 0 )	FillHistogram("PC_Repos_mvsPt_0", mrec, ptrec ) ;
	if(  m > 1.18 && m < 1.2 ) {
	  FillHistogram("ArPodSig0P2",  arpod[1], arpod[0] ) ; 
	  } */
	Double_t rap = Rapidity(  ptrec,  sig0prec.Pz(),  mrec );
	for (Int_t irap = 0; irap < 6; irap++) {
	  if( fabs(rap) < fEtaCuts[irap] ) {
	    sprintf(key,"hPHOS_Sig0MvsPtRec%d",irap) ;
	    FillHistogram(key, mrec, ptrec );	    
	  }
	}
      }    
    }   
  }        // end Lam+Phos  
  /* Lam+PHOS mixed
  for(Int_t icon2=0; icon2<nLam; icon2++){
    AliCaloParticle * ph1 = static_cast<AliCaloParticle*>(fGenpi0Event->At(icon2)) ;
    Double_t ModLam =  ph1->GetXt() ;
    //    printf("MixPC MOdLam %f \n", ModLam);

    Double_t Mlamrec = ph1->M();		
    Double_t corr_factor = MassLam / Mlamrec;		
    TLorentzVector ph1corr;
    ph1corr = corr_factor * *ph1;
    
    for(Int_t ev=0; ev<prevPHOS->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;    
      Int_t nMixPHOS = mixPHOS->GetEntriesFast();
 
      //      printf("--------------mixed PHOS prev i %d  nMIX %d \n", ev, nMixPHOS );
      for(Int_t iPHOS2 = 0; iPHOS2<mixPHOS->GetEntriesFast(); iPHOS2++){
	AliCaloParticle * ph2 = static_cast<AliCaloParticle*>(mixPHOS->At(iPHOS2)) ;
      
	TLorentzVector sig0mix= ph1corr + *ph2 ; 
	Double_t mmix=  sig0mix.M()  ;
	Double_t ptmix=sig0mix.Pt() ;       
	Double_t rapmix = Rapidity(  ptmix,  sig0mix.Pz(),  mmix );
	for (Int_t irap = 0; irap < 6; irap++) {
	  if( fabs(rapmix) < fEtaCuts[irap] ) {
	    sprintf(key,"hPHOS_Sig0MvsPtMix%d",irap) ;
	    FillHistogram(key, mmix, ptmix );
	
	    if( ModLam > 0 ) {
	      sprintf(key,"hPHOS_PSig0MvsPtMix%d",irap) ;
	      FillHistogram(key, mmix, ptmix );
	    }
	    else  if( ModLam < 0 ) {
	      sprintf(key,"hPHOS_ASig0MvsPtMix%d",irap) ;
	      FillHistogram(key, mmix, ptmix );
	    }
	  }
	}
      }
    }
    } */   // END of Lam+PHOS mixed data


        

    
       
    //==========================================================================================================
    // V0s of  Lambda  and  gammma-conv   CONVERSION starts
    if( nConv >= 1 ){
      //    printf("CONV nConv %d nLam %d \n", nConv, nLam );
        
      for(Int_t icon=0; icon<nConv; icon++){
	AliCaloParticle * ph1 = static_cast<AliCaloParticle*>(fConvEvent->At(icon)) ;
        
	FillHistogram("hCaloGamPt", ph1->Pt() ) ;
	FillHistogram("hCaloGamPhiEta", ph1->Phi(), ph1->Eta() ) ;
      }
    }
    
    if( nLam >= 1 ){
      for(Int_t icon=0; icon<nLam; icon++){
	AliCaloParticle * ph1 = static_cast<AliCaloParticle*>(fGenpi0Event->At(icon)) ;
            
	Double_t rap = Rapidity(  ph1->Pt(), ph1->Pz(), ph1->M() );
            
	FillHistogram("hV0LamMvsPt", ph1->M(),  ph1->Pt() ) ;
	Double_t ModLam =  ph1->GetXt();
	//      printf("V0 Lam MOdLam %f \n", ModLam);
	if( ModLam > 0 ) 	FillHistogram("hV0PLamMvsPt", ph1->M(),  ph1->Pt() ) ;
	else if(   ModLam < 0 )	FillHistogram("hV0ALamMvsPt", ph1->M(),  ph1->Pt() ) ;
            
            
	if( ph1->M() > MassLam - 0.5100 &&  ph1->M() < MassLam + 0.5100 ) { // see rapidity distributions from data
	  //	  char key[55] ;
	  for (Int_t irap = 0; irap < 6; irap++) {
	    if( fabs(rap) < fEtaCuts[irap] ) {
	      if( ModLam > 0 ) {
		sprintf(key,"hDtPLamMvsPtRec%d",irap) ;
		FillHistogram(key, ph1->M(), ph1->Pt() );
	      }
	      else  if( ModLam < 0 ) {
		sprintf(key,"hDtALamMvsPtRec%d",irap) ;
		FillHistogram(key, ph1->M(),  ph1->Pt() );
	      }
	    }
	  }
	}
      }
        
      //==========================================================================================================
      // Conversions: Sigma0 = Lambda + gammma-conv
      for(Int_t icon=0; icon<nConv; icon++){
	AliCaloParticle * ph1 = static_cast<AliCaloParticle*>(fConvEvent->At(icon)) ;
        
	for(Int_t icon2=0; icon2<nLam; icon2++){
	  AliCaloParticle * ph2 = static_cast<AliCaloParticle*>(fGenpi0Event->At(icon2)) ;
	  Double_t ModLam =  ph2->GetXt();
	  // printf("-----CONV ModLam %f MassLam %f evt %d tCentr \n",ModLam,ph2->M(),fESDEvent->GetEventNumberInFile(), fCentr );
	  
	  Double_t Mlamrec = ph2->M();
	  Double_t corr_factor = 1; // MassLam / Mlamrec;
	  TLorentzVector	ph2corr;
	  ph2corr = corr_factor * *ph2;
	  
	  TLorentzVector sig0=*ph1 + ph2corr ;
		
	  Double_t pt=sig0.Pt() ;
	  Double_t pos[3]= { ph1->Px(),  ph1->Py(),  ph1->Pz() };
	  Double_t neg[3]= { ph2corr.Px(),  ph2corr.Py(),  ph2corr.Pz() };
	  Double_t moth[3]= { sig0.Px(), sig0.Py(),  sig0.Pz() };
                
	  Double_t arpod[2]= {0,0};
	  GetArPod( pos, neg, moth, arpod );
	  if( (fabs(arpod[1]) > 0.6 && arpod[0]< 0.12) ) {
	    
	    Double_t m=  sig0.M();                                               
	    //                 if (  m > 1.25  || !( Mlamrec > 1.110 &&  Mlamrec < 1.120) ) continue;                
	    if (  m > 1.230  || !( Mlamrec > 1.110 &&  Mlamrec < 1.122) ) continue;  // sa in SelectLambda          
	    // Store Tree variables
	    fLambdaMod = ph2->GetXt();
	    fLambdaPx = ph2->Px();
	    fLambdaPy = ph2->Py();
	    fLambdaPz = ph2->Pz();
	    fLambdaMass = ph2 -> M();
	    fLambdaCosPointingAngle = ph2->GetV0CosPointingAngle();
	    fLambdaDCADaughters = ph2->GetV0DCADaughters();
	    fLambdaDCAtoPVNeg = ph2->GetV0IPNeg();
	    fLambdaDCAtoPVPos = ph2->GetV0IPPos();
	    fLambdaRadius = ph2->GetV0Radius();
	    fLambdaArmPt = ph2->GetArmPt();
	    fLambdaArmAlpha = ph2->GetArmAlpha();
	    fLambdaEta = ph2->Eta();

	    //	fLambdaChi2 = ph2->GetV0Chi2();
	    
	    //	fLambdaNumberofCrossedRows = ;
	    //	fLambdaRatioofCrossedRows = ;
	    //	fLambdaRunNumber = fESDEvent->GetRunNumber();
	    
	    /* ABB 8jan15 - fill Tree variables in the selection routine! 
            // have tune variables to V0Reader !!! abb-2dec15
	    fGammaPx = ph1->Px();
	    fGammaPy = ph1->Py();
	    fGammaPz = ph1->Pz();
	    fGammaMass = ph1 -> M();
	    fGammaCosPointingAngle = ph1->GetV0CosPointingAngle();
	    fGammaDCADaughters = ph1->GetV0DCADaughters();
	    fGammaDCAtoPVNeg = ph1->GetV0IPNeg();
	    fGammaDCAtoPVPos = ph1->GetV0IPPos();
	    fGammaRadius = ph1->GetV0Radius();
	    fGammaArmPt = ph1->GetArmPt();
	    fGammaArmAlpha = ph1->GetArmAlpha();
	    fGammaEta = ph1->Eta();
	    fGammaChi2 = ph1->GetV0Chi2();
	    fGammaZConv = ph1->GetZConv(); 
	    */

	    fSigmaMass = sig0.M();
	    fSigmaPt = sig0.Pt();
	    fSigmaArmPt = arpod[0];
	    fSigmaArmAlpha = arpod[1];
	    fCentr = fCentrality;             

	    //	    fTreeVariablePt = sqrt( fLambdaPx*fLambdaPx + fLambdaPy*fLambdaPy );

	    tTreeEvent->Fill();
	    fTreeV0->Fill();
	    //	    fTreeGammaV0->Fill();
	                    
	    printf("Sig0-CC  M %f PT %f  E %f nconv %d Centr %f \n", sig0.M(), sig0.Pt(), sig0.E(), nConv, fCentr );
                            
	    FillHistogram("hCCReMvsPtnArPod",m,pt ) ;
                
                    
	    FillHistogram("hCCReMvsPt",m,pt ) ;
	    if( ModLam < 0 ) 	FillHistogram("hCCRenegMvsPt",m,pt ) ;
	    else if(   ModLam > 0 )	FillHistogram("hCCReposMvsPt", m, pt ) ;
                 
                    
	    Double_t rap = Rapidity(  pt,  sig0.Pz(),  m );
	    char key[55] ;
	    for (Int_t irap = 0; irap < 6; irap++) {
	      if( fabs(rap) < fEtaCuts[irap] ) {
		if( ModLam > 0 ) {
		  sprintf(key,"hDtPSig0MvsPtRec%d",irap) ;
		  FillHistogram(key, m, pt );
		}
		else  if( ModLam < 0 ) {
		  sprintf(key,"hDtASig0MvsPtRec%d",irap) ;
		  FillHistogram(key, m, pt );
		}
	      }
	    }

	    /*     for (Int_t icent = 0; icent < 4; icent++) {
	      if(  fabs(fCentr) >= fCentCuts[icent] && fabs(fCentr) < fCentCuts[icent+1] ) {
		if( ModLam > 0 ) {
		  sprintf(key,"hCentDtPSig0MvsPtRec%d",icent) ;
		  FillHistogram(key, m, pt );
		}
		else  if( ModLam < 0 ) {
		  sprintf(key,"hCentDtASig0MvsPtRec%d",icent) ;
		  FillHistogram(key, m, pt );
		}
	      }
	      }                      
	    FillHistogram("hSig0CaloRapALL", rap  ) ;
	    FillHistogram("hSig0CaloRapvsPt", rap, pt  ) ;
	    
	    FillHistogram("hSig0CaloRPhiEta", rap  ) ;
	    //                    FillHistogram("rec_Sig0", pt  ) ;
	    //  if( fabs(rap) < 1.4 )  FillHistogram("rec_Sig0eta1", pt) ;
	    FillHistogram("ArPodSig0C",  arpod[1], arpod[0] ) ;
	    */	    

	    TLorentzVector sig0prec=*ph1 + *ph2 ;
	    Double_t mrec=sig0prec.M() ;
	    Double_t ptrec=sig0prec.Pt() ;
	    Double_t posr[3]= { ph1->Px(),  ph1->Py(),  ph1->Pz() };
	    Double_t negr[3]= { ph2->Px(),  ph2->Py(),  ph2->Pz() };
	    Double_t mothr[3]= { sig0prec.Px(),  sig0prec.Py(), sig0prec.Pz() };
	    Double_t arpodr[2]= {0,0};

	    GetArPod( posr, negr, mothr, arpodr );
	    if( !( fabs(  arpodr[1] ) > 0.6 && arpodr[0]< 0.12 ) ) continue ;	    
	    //	    Double_t rapmix = Rapidity(  ptrec,  sig0prec.Pz(),  mrec );
	    // char key[55] ;
	    //printf("MIX Sig0mix-CC  M %f PT %f E %f nLam %d Centr %f \n",sig0mix.M(),sig0mix.Pt(),sig0mix.E(),nLam,fCentr);


	    /*   FillHistogram("hCCReMvsPtnArPod0",m,pt ) ;
	    if( ( fabs(  arpodr[1] ) > 0.5 && arpodr[0]< 0.14 ) )  FillHistogram("hCCReMvsPt1", mrec, ptrec ) ;
	    if( ( fabs(  arpodr[1] ) > 0.6 && arpodr[0]< 0.12 ) ) {
	      FillHistogram("hCCReMvsPt0", mrec, ptrec ) ;
	      if( ModLam < 0 ) 	FillHistogram("hCCRenegMvsPt0",mrec,ptrec ) ;
	      else if(   ModLam > 0 )	FillHistogram("hCCReposMvsPt0", mrec, ptrec ) ;
	      
	      if(  m > 1.18 && m < 1.2 ) {
		FillHistogram("ArPodSig0C2",  arpod[1], arpod[0] ) ;
	      }
	      } */
	  }
	}
      }
       
      //Lam + mixed Conv Photon
      for(Int_t icon2=0; icon2<nLam; icon2++){
	AliCaloParticle * ph2 = static_cast<AliCaloParticle*>(fGenpi0Event->At(icon2)) ;
	Double_t ModLam =  ph2->GetXt() ;
        
	Double_t corr_factor = 1.0 ;
	TLorentzVector	ph2corr;
	ph2corr = corr_factor * *ph2;
        
	for(Int_t ev=0; ev< prevConv->GetSize();ev++){
	  TClonesArray * mixConv = static_cast<TClonesArray*>(prevConv->At(ev)) ;
	  
	  for(Int_t icon=0; icon<mixConv->GetEntriesFast() ; icon++){
	    
	    AliCaloParticle * ph1 = static_cast<AliCaloParticle*>(mixConv->At(icon)) ;
	    
	    TLorentzVector sig0mix=*ph1 + ph2corr ;
	    Double_t mmix= fabs( sig0mix.M() ) ;
	    Double_t ptmix=sig0mix.Pt() ;
	    // FillHistogram("hCCMiMvsPt",mmix,ptmix ) ;
                    
	    // if( ModLam < 0 )  FillHistogram("hCCMinegMvsPt",mmix,ptmix ) ;
	    // else if( ModLam > 0 )  FillHistogram("hCCMiposMvsPt",mmix,ptmix ) ;
	    
	    Double_t posr[3]= { ph1->Px(),  ph1->Py(),  ph1->Pz() };
	    Double_t negr[3]= { ph2corr.Px(),  ph2corr.Py(),  ph2corr.Pz() };
	    Double_t mothr[3]= { sig0mix.Px(),  sig0mix.Py(), sig0mix.Pz() };
	    Double_t arpodr[2]= {0,0};
	    
	    GetArPod( posr, negr, mothr, arpodr );
	    if( !( fabs(  arpodr[1] ) > 0.6 && arpodr[0]< 0.12 ) ) continue ;
	    
	    Double_t rapmix = Rapidity(  ptmix,  sig0mix.Pz(),  mmix );
	    char key[55] ;
	    //printf("MIX Sig0mix-CC  M %f PT %f E %f nLam %d Centr %f \n",sig0mix.M(),sig0mix.Pt(),sig0mix.E(),nLam,fCentr);
	    for (Int_t irap = 0; irap < 6; irap++) {
	      if( fabs(rapmix) < fEtaCuts[irap] ) {
		if( ModLam > 0 ) {
		  sprintf(key,"hDtPSig0MvsPtMix%d",irap) ;
		  FillHistogram(key, mmix, ptmix );
		}
		else  if( ModLam < 0 ) {
		  sprintf(key,"hDtASig0MvsPtMix%d",irap) ;
		  FillHistogram(key, mmix, ptmix );
		}
	      }
	    }
	    /*  for (Int_t icent = 0; icent < 4; icent++) {
	      if(  fabs(fCentr) >= fCentCuts[icent] && fabs(fCentr) < fCentCuts[icent+1] ) {
	      //	  printf(" cmin %f cmax %f \n", fCentCuts[icent], fCentCuts[icent+1]);
		
		if( ModLam > 0 ) {
		  sprintf(key,"hCentDtPSig0MvsPtMix%d",icent) ;
		  FillHistogram(key, mmix, ptmix );
		}
		else  if( ModLam < 0 ) {
		  sprintf(key,"hCentDtASig0MvsPtMix%d",icent) ;
		  FillHistogram(key, mmix, ptmix );
		}
	      }
	      } */                    
	  }
	}
      }    // END of Conv mixed data
    }     // END of if for Conv and Lam     
}      //END RealMixed

//_____________________________________________________
void AliAnalysisTaskSigma0Spectra::FillCorr200(TLorentzVector * trig, const Int_t itype  )
{  
    Double_t xesum = 0;
    Double_t yesum = 0;
    
    Int_t n=fTrackEvent->GetEntriesFast() ;
    Int_t itype2 = itype*10;
    
    char key0[155]; char key1[155]; char key2[155]; char key3[155]; char key4[155]; char key5[55]; 
  
    for(Int_t i=0; i<n;i++){
        
        Int_t hnum = -1;
        
        Double_t ptTrig = trig->Pt();
        printf("FillCorr2 itype2 %d ptTrig %f \n",itype2,ptTrig );    
        
        AliCaloParticle * partn = static_cast<AliCaloParticle*>(fTrackEvent->At(i)) ;
        if(trig==partn) //for the case of track trigger
        continue ;
        
        Double_t dphi=trig->Phi()-partn->Phi() ;
        
        if(dphi <= -TMath::PiOver2()) dphi+=TMath::TwoPi();
        if(dphi > 3*TMath::PiOver2()) dphi-=TMath::TwoPi();
        
        
        Double_t xe=-partn->Pt()*TMath::Cos(dphi)/trig->Pt() ;     
        Double_t ye=-partn->Pt()*TMath::Sin(dphi)/trig->Pt() ;  
        
        xesum += xe ;  // if( xe > 0 ) xesumpos += xe;  else if( xe < 0 ) xesumneg +=  fabs( xe };
        yesum += ye ;
        
        Double_t pttr = partn->Pt();
        
        Double_t cosdelphi = TMath::Cos(dphi) ;
        if( cosdelphi == 0 ) cosdelphi = 1.e-6;
        Double_t rat= -pttr/ptTrig * cosdelphi / fabs( cosdelphi )  ;
        
        if( (ptTrig > 1  && ptTrig <= 5  ) && (pttr > 1 && pttr <= 5  )  )  hnum  = itype2 ;
        if( (ptTrig > 5 && ptTrig < 10 ) && (pttr > 1  && pttr <= 5  )  )  hnum  = itype2 + 1 ;
        if( (ptTrig > 5 && ptTrig < 10 ) && (pttr > 5  && pttr <10 )  )  hnum  = itype2 + 2 ;
        if( (ptTrig > 10 && ptTrig < 20 ) && (pttr > 1  && pttr <= 5  )  )  hnum  = itype2 + 3 ;
        if( (ptTrig > 10 && ptTrig < 20 ) && (pttr > 5  && pttr <= 10 )  )  hnum  = itype2 + 4 ;
        if( (ptTrig > 10 && ptTrig < 20 ) && (pttr > 10 && pttr < 20 )  )  hnum  = itype2 + 5 ;
        
        sprintf(key0,"hdelphi%d",hnum); sprintf(key1,"hrat%d",hnum); 
        sprintf(key2,"hxe%d",hnum); sprintf(key4,"hye%d",hnum);   sprintf(key3,"hmgg%d",hnum);
        sprintf(key5,"Multiplicity%d",hnum);
        
        FillHistogram(key0,dphi);
        FillHistogram(key1,rat); FillHistogram(key2,xe); 
        FillHistogram(key3, trig->M() );   
        FillHistogram(key4,ye);
        FillHistogram(key5,fMultiplicity);
        
    }
}


//______________________________________________________________________________
void AliAnalysisTaskSigma0Spectra::GetArPod( Double_t pos[3], Double_t neg[3], Double_t moth[3],  Double_t arpod[2] ){
    
    //see header file for documentation
    
    TVector3 momentumVectorPositiveKF(pos[0],pos[1],pos[2]);
    TVector3 momentumVectorNegativeKF(neg[0],neg[1],neg[2]);
    TVector3 vecV0(moth[0],moth[1],moth[2]);
    
    if( momentumVectorPositiveKF.Mag()*vecV0.Mag() == 0. || momentumVectorNegativeKF.Mag() * vecV0.Mag() == 0 ) {
        arpod[0] = 0;
        arpod[1] = 0;
        return;
    }
    
    
    Float_t costhetaV0pos=( momentumVectorPositiveKF* vecV0)/(momentumVectorPositiveKF.Mag() * vecV0.Mag()) ;
    Float_t costhetaV0neg=( momentumVectorNegativeKF* vecV0)/(momentumVectorNegativeKF.Mag() * vecV0.Mag());
    
    if ( fabs( costhetaV0pos ) > 1 ||  fabs( costhetaV0neg ) > 1 ) {
        arpod[0] = 0;
        arpod[1] = 0;
        return;
    }
    
    Float_t thetaV0pos=TMath::ACos(( momentumVectorPositiveKF* vecV0)/(momentumVectorPositiveKF.Mag() * vecV0.Mag()));
    Float_t thetaV0neg=TMath::ACos(( momentumVectorNegativeKF* vecV0)/(momentumVectorNegativeKF.Mag() * vecV0.Mag()));
    
    Float_t alfa =((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)-(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg))/
    ((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;
    
    Float_t qt = momentumVectorPositiveKF.Mag()*TMath::Sin(thetaV0pos);
    
    arpod[0]=qt;
    arpod[1]=alfa;
    
}




