 /**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Dmitri Peressounko (RRC KI), Alexander Borissov (WSU)          *
 *                                                                        *  
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class used to do analysis on conversion + calorimeter pairs
// A lot of code cut-and-pasted from AliV0Reader and GammaConversion classes
// 
//---------------------------------------------
////////////////////////////////////////////////
// root
#include "TH3.h"
#include "TClonesArray.h"
#include "AliAnalysisTaskSigma0Corr.h"
//#include "AliStack.h"
//#include "AliLog.h"
//#include "AliPHOSGeoUtils.h" 
//#include "AliEMCALGeometry.h" 
#include "AliCaloParticle.h"


ClassImp(AliAnalysisTaskSigma0Corr)

//______________________________________________________________________________
AliAnalysisTaskSigma0Corr::AliAnalysisTaskSigma0Corr():
AliAnalysisTaskSigma0Spectra()
{
}
//______________________________________________________________________________
AliAnalysisTaskSigma0Corr::AliAnalysisTaskSigma0Corr(const char* name):
  AliAnalysisTaskSigma0Spectra(name)
{
}
//_____________________________________________________
AliAnalysisTaskSigma0Corr::~AliAnalysisTaskSigma0Corr() 
{
  // Remove all pointers
}
//_____________________________________________________
void AliAnalysisTaskSigma0Corr::Init()
{
  AliAnalysisTaskSigma0Spectra::Init() ;
}
//____________________________________________________________
void AliAnalysisTaskSigma0Corr::UserCreateOutputObjects()
{
    AliAnalysisTaskSigma0Spectra::UserCreateOutputObjects();

//  char key[55] ;
//  Int_t nTrig=20 ;
//  Double_t ptTrigMax=40. ;	      
//  Int_t nPartn=20 ;
//  Double_t ptPartnMax=10 ;
  // Int_t nPhi=50 ;
  // Double_t phiMin=-TMath::PiOver2() ;
  //Double_t phiMax=3*TMath::PiOver2() ;

  //    fOutputContainer->Add(new TH2F("hrat-hxe-0","ratio=pT-trig/pT vx xE", 10,-1.2,1.2, 10,-1.5,1.5  )) ;	      
  //  fOutputContainer->Add(new TH2F("hrat-hxe-1","ratio=pT-trig/pT vx xE", 40,-1.2,1.2, 30,-1.5,1.5  )) ;	      
  //  fOutputContainer->Add(new TH2F("hrat-hxe-2","ratio=pT-trig/pT vx xE", 40,-1.2,1.2, 30,-1.5,1.5  )) ;	      
 //  fOutputContainer->Add(new TH2F("hrat-hxe-3","ratio=pT-trig/pT vx xE", 40,-1.2,1.2, 30,-1.5,1.5  )) ;	      
  //  fOutputContainer->Add(new TH2F("hrat0-hxe-4","ratio=pT-trig/pT vx xE", 40,-1.2,1.2, 30,-1.5,1.5  )) ;	      

  /*   fOutputContainer->Add(new TH3F("PHOS_single_track_all","Correlation with PHOS cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ; 
  fOutputContainer->Add(new TH3F("PHOS_single_track_DefIsol","Correlation with PHOS cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("PHOS_single_track_Direct","Correlation with PHOS cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("PHOS_single_track_DirectIsolated","Correlation with PHOS cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ; */

  /*
  //  fOutputContainer->Add(new TH3F("EMCAL_single_track_all","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;

  fOutputContainer->Add(new TH3F("EMCAL_single_track_all","Correlation with PHOS cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_single_track_DefIsol","Correlation with PHOS cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_single_track_Direct","Correlation with PHOS cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_single_track_DirectIsolated","Correlation with PHOS cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;

  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_DefIsol","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_Direct","Correlation with  EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_DirectIsolated","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  */

  /* xE distributions
  Int_t nXe=100 ;

  fOutputContainer->Add(new TH3F("PHOS_single_track_all_Xe","Correlation with PHOS cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("PHOS_single_track_DefIsol_Xe","Correlation with PHOS cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("PHOS_single_track_Direct_Xe","Correlation with PHOS cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("PHOS_single_track_DirectIsolated_Xe","Correlation with PHOS cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;  */

  /*
  fOutputContainer->Add(new TH3F("EMCAL_single_track_all_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_single_track_DefIsol_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_single_track_Direct_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_single_track_DirectIsolated_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  */

  /*  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_DefIsol_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_Direct_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_DirectIsolated_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  */

  /*
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_isol","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_isol_peak","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_isol_side","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  // fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_incl","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_incl_peak","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_incl_side","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_isol","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_isol_peak","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_isol_side","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  // fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_incl","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_incl_peak","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_incl_side","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nPhi,phiMin,phiMax)) ;
  


  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_isol_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_isol_peak_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_isol_side_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  // fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_incl","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_incl_peak_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_all_incl_side_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_isol_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_isol_peak_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_isol_side_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  // fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_incl","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_incl_peak_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  fOutputContainer->Add(new TH3F("EMCAL_pi0_track_lead_incl_side_Xe","Correlation with EMCAL cluster",nTrig,0.,ptTrigMax,nPartn,0.,ptPartnMax,nXe,-1.5,1.5)) ;
  */
  
  
	      	     
  //  fOutputContainer->Add(new TH2F("hxyesum","xyEsum", 30,-1.5,1.5, 30,-1.5,1.5  )) ; // NO TH2F allowed - memory???
  //    fOutputContainer->Add(new TH2F("hxyesum","PP Mass vs pt", 5,0.,0.7, 5,0.,10 )) ;


  //  if( 1> 0 ) return;

    //  char key1[155] ;  //  char key2[155] ;  

//  for(Int_t m3=1; m3<= 7; m3++){
//   for(Int_t m0=1; m0<= 2; m0++){
//     //  for(Int_t m1=0; m1< 2; m1++){
 
//     for(Int_t m2=0; m2< 1; m2++){   
//       // Int_t  m = m0*10+m1*100+m2;
//       Int_t  m = m3*100 + m0*10+m2;

//       snprintf(key1,155,"hdelphi-%d",m) ;
//       fOutputContainer->Add(new TH1F(key1,"delphi",94,-3.14,6.28)) ;
//       snprintf(key1,155,"hrat-%d",m) ;
//       fOutputContainer->Add(new TH1F(key1,"ratio=pTa/pTtrig",150,-3., 3.)) ;
//       snprintf(key1,155,"hxe-%d",m) ;
//       fOutputContainer->Add(new TH1F(key1,"xE",150,-3.,3. )) ;
//       snprintf(key1,155,"hye-%d",m) ;
//       fOutputContainer->Add(new TH1F(key1,"yE",150,-3.,3. )) ;
//       snprintf(key1,155,"hmgg-%d",m) ;
//       fOutputContainer->Add(new TH1F(key1,"Mass of gamma+gamma",70, 0., 0.7 )) ;
 
//       snprintf(key1,155,"Multiplicity-%d",m) ;
//       fOutputContainer->Add(new TH1F(key1,"Multiplicity",71, -0.5, 70.5 )) ;
      
//       //    snprintf(key1,155,"hrat-hxe-%d",m) ;
//       //    fOutputContainer->Add(new TH2F(key1,"ratio=pT-trig/pT vx xE", 40,-1.2,1.2, 30,-1.5,1.5  )) ;
//     }
   
//     /*    Int_t msum =  m0*10+m1*100;
//     snprintf(key1,155,"hxesum-%d",msum) ; 
//     fOutputContainer->Add(new TH1F( key1,"xEsum",150,-3.,3. )) ;
//     snprintf(key1,155,"hyesum-%d",msum) ; 
//     fOutputContainer->Add(new TH1F( key1,"yEsum",150,-3.,3. )) ;       
//     */
//   }
//  }
  // }
  PostData(1, fOutputContainer); 
}

//_____________________________________________________
//_____________________________________________________
//_____________________________________________________
void AliAnalysisTaskSigma0Corr::UserExec(Option_t *option)
{

  // Execute analysis for current event
  

   AliAnalysisTaskSigma0Spectra::UserExec(option) ;

   //    if ( 1>0 ) return;



   if( 1>0 ){
     PostData(1, fOutputContainer);
     return;
   }

  //Fill dphi correlations
  //First - single PHOS cluster vs tracks
  /*  SKIP   if(fPHOSEvent->GetEntriesFast()>0 && fLeadingPHOS>=0){ //it exist 
    AliCaloParticle * trig = static_cast<AliCaloParticle*>(fPHOSEvent->At(fLeadingPHOS)) ;
    // printf("PHOS trig2 (%d)=%p \n",fLeadingPHOS,trig) ;
    FillCorrelation(trig,"PHOS_single_track_all") ; 

    Int_t isolation = trig->IsIsolated();
    Int_t kind = 2 ; //Cone=0.4, epsilon=0.1, tracks only

    //    printf("PP isolation2= %d, kind=%d, res=%d \n",isolation,kind,(isolation&kind)) ;

     if((isolation&kind !=0))
       FillCorrelation(trig,"PHOS_single_track_DefIsol") ;   
    if(!trig->IsTagged())
       FillCorrelation(trig,"PHOS_single_track_Direct") ; 
    if(!trig->IsTagged() && (isolation&kind) )
    FillCorrelation(trig,"PHOS_single_track_DirectIsolated") ; 
    } */


  // printf(".. .. .. .\n");



  //  if( fTrackEvent->GetEntriesFast()>0 && fLeadingTrack >= 0 ){
    //    AliCaloParticle * ltrk = static_cast<AliCaloParticle*>(fTrackEvent->At(fLeadingTrack));
    //    printf("Lead Track num %d Pt %f   \n",  fLeadingTrack, fELeadingTrack );  
  //  }


  //----------   pi^0 --------PHOS-------------------------------------
  if(fPHOSPi0Event->GetEntriesFast()>0 && fLeadingPi0PHOS>=0   &&  fELeadingPi0PHOS  > fELeadingTrack ){ 

    //    AliCaloParticle * trig = static_cast<AliCaloParticle*>(fPHOSEvent->At(fLeadingPHOS)) ;

    printf("Lead pi^0 PHOS %d Pt %f n-pi0 %d  \n",  fLeadingPi0PHOS, fELeadingPi0PHOS, fPHOSPi0Event->GetEntriesFast() ); 

    for(Int_t i=0; i< fPHOSPi0Event->GetEntriesFast(); i++){ 
      AliCaloParticle * pi0 = static_cast<AliCaloParticle*>(fPHOSPi0Event->At(i)) ;

      FillDeltaPhi(pi0, 1 ) ;
    }
  }


  //  if ( 1>0 ) return;

  /*
  //2nd - single EMCAL cluster vs tracks
  if(fEMCALEvent->GetEntriesFast()>0 && fLeadingEMCAL>=0){ //it exist 
    AliCaloParticle * trig = static_cast<AliCaloParticle*>(fEMCALEvent->At(fLeadingEMCAL)) ;

    //       FillCorrelation(trig,"EMCAL_single_track_all") ; 

    //    Int_t itype = 0 ;
    Int_t isolation = trig->IsIsolated();
    Int_t kind = 2048 ; //Cone=0.4, epsilon=0.1, tracks only

    // printf("Lead photon EMCAL num %d Pt %f \n",  fLeadingEMCAL, fELeadingEMCAL );

    //    printf("EE isolation4= %d, kind=%d, res=%d \n \n",isolation,kind,(isolation&kind)) ;

	//    kind = 1 ;

    //    printf("EE isolation1= %d, kind=%d, res=%d \n",isolation,kind,(isolation&kind)) ;
    if((isolation&kind !=0))
       FillCorrelation(trig,"EMCAL_single_track_DefIsol") ;   
    if(!trig->IsTagged())
       FillCorrelation(trig,"EMCAL_single_track_Direct") ; 
    if(!trig->IsTagged() && (isolation&kind) )
       FillCorrelation(trig,"EMCAL_single_track_DirectIsolated") ;  
  */

  
  //3d - pi0 from EMCAL cluster vs tracks -----------------------------------------
  // abb - jan12

  if( fEMCALPi0Event->GetEntriesFast() <=0 ) return;

  if(fEMCALPi0Event->GetEntriesFast()>0 && fLeadingPi0EMCAL>=0   &&  fELeadingPi0EMCAL  > fELeadingTrack ){

    //    printf("Lead pi^0 EMCAL %d Pt %f   \n",  fLeadingPi0EMCAL, fELeadingPi0EMCAL ); 
    for(Int_t i=0; i< fEMCALPi0Event->GetEntriesFast(); i++){ 
      AliCaloParticle * pi0 = static_cast<AliCaloParticle*>(fEMCALPi0Event->At(i)) ;

      //- all pi0 = incl + lead
      // FillCorrelation(pi0,"EMCAL_pi0_track_all") ;  
      FillDeltaPhi(pi0, 2 ) ; 
      // Int_t isolation = pi0->IsIsolated();
      // Int_t kind = 2048 ;
    }
  }

    //   printf("EE-CORR isolation1= %d, kind=%d, res=%d \n",isolation,kind,(isolation&kind)) ;
  
    /*
    //    if((isolation == 4095 ) )  {
    if((isolation&kind !=0))  {
      // FillCorrelation(pi0,"EMCAL_pi0_track_all_isol") ;  
      FillDeltaPhi(pi0, 2 ) ; 
      if( pi0->M() > 0.12 && pi0->M()<0.16 ) {
	// FillCorrelation(pi0,"EMCAL_pi0_track_all_isol_peak") ;  
	FillDeltaPhi(pi0, 3 ) ; 
      }
      else if(  (pi0->M() > 0.06 && pi0->M()<0.09) ||  (pi0->M() > 0.190 && pi0->M()<0.25) ) {
	// FillCorrelation(pi0,"EMCAL_pi0_track_all_isol_side") ;  
	FillDeltaPhi(pi0, 4 ) ; 
      }
    }
    */  

    /*  
    // if( 1>0  )  {             // inclusive - all data
    // FillCorrelation(pi0,"EMCAL_pi0_track_all_incl") ;  // FillDeltaPhi(pi0, 5 ) ; 
    if( pi0->M() > 0.12 && pi0->M()<0.16 ) {
      // FillCorrelation(pi0,"EMCAL_pi0_track_all_incl_peak") ;   
      FillDeltaPhi(pi0, 5 ) ; 
    }
    else if(  (pi0->M() > 0.06 && pi0->M()<0.09) ||  (pi0->M() > 0.190 && pi0->M()<0.25) ) {
      // FillCorrelation(pi0,"EMCAL_pi0_track_all_incl_side") ;   
      FillDeltaPhi(pi0, 6 ) ; 
    }

    if((isolation == 0))  {
      // FillCorrelation(pi0,"EMCAL_pi0_track_all_isol") ;  
      FillDeltaPhi(pi0, 7 ) ; 
      if( pi0->M() > 0.12 && pi0->M()<0.16 ) {
	// FillCorrelation(pi0,"EMCAL_pi0_track_all_isol_peak") ;  
	FillDeltaPhi(pi0, 8 ) ; 
      }
      else if(  (pi0->M() > 0.06 && pi0->M()<0.09) ||  (pi0->M() > 0.190 && pi0->M()<0.25) ) {
	// FillCorrelation(pi0,"EMCAL_pi0_track_all_isol_side") ;  
	FillDeltaPhi(pi0, 9 ) ; 
      }
    }
    */

    // only leading pi0
    /*    if( i == fLeadingPi0EMCAL ){
      //      FillCorrelation(pi0,"EMCAL_pi0_track_lead") ;   

      FillDeltaPhi(pi0, 11 ) ; 
      Int_t isolation = pi0->IsIsolated();

      //      if((isolation == 4095 ) )  {
      Int_t kind = 2048 ;
      if((isolation&kind !=0))  {
	// FillCorrelation(pi0,"EMCAL_pi0_track_lead_isol") ;   
	FillDeltaPhi(pi0, 12 ) ; 
	if( pi0->M() > 0.12 && pi0->M()<0.16 ) {
	  //  FillCorrelation(pi0,"EMCAL_pi0_track_lead_isol_peak") ;   
	  FillDeltaPhi(pi0, 13 ) ; 
	}
	else if(  (pi0->M() > 0.06 && pi0->M()<0.09) ||  (pi0->M() > 0.190 && pi0->M()<0.25) ) {
	  //  FillCorrelation(pi0,"EMCAL_pi0_track_lead_isol_side") ;  
	  FillDeltaPhi(pi0, 14 ) ; 
	}
      }
      */
      //  if( 1>0 )  {        // inclusive - all leading
      //	FillCorrelation(pi0,"EMCAL_pi0_track_lead_incl") ;  //	FillDeltaPhi(pi0, 15 ) ; 

      /*
      if( pi0->M() > 0.12 && pi0->M()<0.16 ) {
	//  FillCorrelation(pi0,"EMCAL_pi0_track_lead_incl_peak") ;  
	FillDeltaPhi(pi0, 15 ) ; 

	//   printf("Lead pi0 EMCAL num %d Pt %f \n", i, fELeadingPi0EMCAL );

      }
      else if(  (pi0->M() > 0.06 && pi0->M()<0.09) ||  (pi0->M() > 0.190 && pi0->M()<0.25) ) {
	//	FillCorrelation(pi0,"EMCAL_pi0_track_lead_incl_side") ;  
	FillDeltaPhi(pi0, 16 ) ; 
      }
      */

      /*
      if((isolation  == 0))  {
	// FillCorrelation(pi0,"EMCAL_pi0_track_lead_isol") ;   
	FillDeltaPhi(pi0, 17 ) ; 
	if( pi0->M() > 0.12 && pi0->M()<0.16 ) {
	  //  FillCorrelation(pi0,"EMCAL_pi0_track_lead_isol_peak") ;   
	  FillDeltaPhi(pi0, 18 ) ; 
	}
	else if(  (pi0->M() > 0.06 && pi0->M()<0.09) ||  (pi0->M() > 0.190 && pi0->M()<0.25) ) {
	  //  FillCorrelation(pi0,"EMCAL_pi0_track_lead_isol_side") ;  
	  FillDeltaPhi(pi0, 19 ) ; 
	}
	}      
    } 
      */  // end of leading pi0 
    //  }   // end of EMCAL pi0 events


    /* DIMAP - Dec11:
  if( fEMCALPi0Event->GetEntriesFast()>0 && fLeadingPi0EMCAL>=0){ //it exist 
    AliCaloParticle * trig = static_cast<AliCaloParticle*>(fEMCALPi0Event->At(fLeadingPi0EMCAL)) ;
    //    printf("EMCAL trig2 (%d)=%p \n",fLeadingPi0EMCAL,trig) ;
    
    FillCorrelation(trig,"EMCAL_pi0_track_all") ; 

    printf("Corr-EE-LEADING  E %f numb %d   \n", fELeadingPi0EMCAL, fLeadingPi0EMCAL  );  
    //    printf("Corr-EE-trig  m %f p %f    pt %f \n", trig->M(), trig->P(), trig->Pt() );
    //    TLorentzVector * pi = static_cast<TLorentzVector*>(fEMCALPi0Event->At(fLeadingPi0EMCAL)) ;
    //  printf("Corr-EMCAL-pi flead %d  M %f E  %f p %f  \n \n ",fLeadingPi0EMCAL, pi->M(), pi->E() , pi->P() ) ;

    Int_t itype = 0 ;
    FillDeltaPhi(pi0, itype) ; 

    Int_t isolation = trig->IsIsolated();
    Int_t kind = 3 ; 

    //    printf("EE isolation2= %d, kind=%d, res=%d mpi %f \n",isolation,kind,(isolation&kind), trig.M() ) ;
    //    printf("EE isolation-all = %d, kind=%d, res=%d \n",isolation,kind,(isolation&kind)) ;

    if((isolation&kind !=0))
       FillCorrelation(trig,"EMCAL_pi0_track_DefIsol") ;   
    if(!trig->IsTagged())
       FillCorrelation(trig,"EMCAL_pi0_track_Direct") ;     

    if(!trig->IsTagged() &&  (isolation&kind) ){
       FillCorrelation(trig,"EMCAL_pi0_track_DirectIsolated") ; 
       itype = 5 ;
       FillDeltaPhi(pi0, itype) ; 
    }

    kind = 2; //Cone=0.4, epsilon=0.1, tracks only

    if(!trig->IsTagged() && (isolation&kind) ){


      // FillCorrelation(trig,"EMCAL_pi0_track_DirectIsolated") ; 
       itype = 10 ;
       FillDeltaPhi(pi0, itype) ; 
    }
  }
    */

  //4th - pi0 from EMCAL cluster + Conv  vs tracks -----------------------------------------
  //  if( fConvPi0Event->GetEntriesFast()>0 && fLeadingPi0EMConv>=0){ //it exist 
  //  printf("Corr-EC-LEADING  E %f numb %d   \n", fELeadingPi0EMConv, fLeadingPi0EMConv  );  
  // }
    
//  ProcessMC();
  PostData(1, fOutputContainer);
}
//_____________________________________________________
void AliAnalysisTaskSigma0Corr::FillCorrelation(AliCaloParticle * trig, const char * kind)
{
  
   char kindXe[55] ;
   sprintf(kindXe,"%s_Xe",kind) ;
   Int_t n=fTrackEvent->GetEntriesFast() ;

   //   printf( " kind %s \n", kind);
   // if( kind == "EMCAL_single_track_DirectIsolated")    printf( "ISOLATED kind %s \n", kind);

   for(Int_t i=0; i<n;i++){
     AliCaloParticle * partn = static_cast<AliCaloParticle*>(fTrackEvent->At(i)) ;
     if(trig==partn) //for the case of track trigger
       continue ;
     Double_t dphi=trig->Phi() - partn->Phi() ;
     //     while(dphi<-TMath::PiOver2())dphi+=TMath::TwoPi() ;
     //  while(dphi>3*TMath::PiOver2())dphi-=TMath::TwoPi() ;

     if(dphi <= -TMath::PiOver2()) dphi+=TMath::TwoPi();
     if(dphi > 3*TMath::PiOver2()) dphi-=TMath::TwoPi();

     FillHistogram(kind,trig->Pt(),partn->Pt(),dphi) ;
   
     Double_t xe=-partn->Pt()*TMath::Cos(dphi)/trig->Pt() ;
     FillHistogram(kindXe,trig->Pt(),partn->Pt(),xe) ;
   }
}
//_____________________________________________________
void AliAnalysisTaskSigma0Corr::FillDeltaPhi(AliCaloParticle * trig, const Int_t itype  )
{
  
  Double_t xesum = 0; //  Double_t xesumpos = 0;   Double_t xesumneg = 0; 
  Double_t yesum = 0;

  Int_t n=fTrackEvent->GetEntriesFast() ;

  Int_t itype2 = itype*10;

  char key0[155]; char key1[155]; char key2[155]; char key3[155]; char key4[155]; char key5[55]; 
   //   printf( " itype %d \n", itype );   //   sprintf(key0,"hdelphi-",itype) ;   ///  sprintf(key1,"hrat-",itype) ;
   //   sprintf(key2,"hxe-",itype) ;   ///  sprintf(key3,"hrat-hxe-",itype) ;   //   sprintf(key4,"hmvspt-",itype) ;   
   //   printf("%s %s ",key0, key1 ) ;

  for(Int_t i=0; i<n;i++){

    Int_t hnum = -1;
    

    AliCaloParticle * partn = static_cast<AliCaloParticle*>(fTrackEvent->At(i)) ;
    if(trig==partn) //for the case of track trigger
      continue ;

    Double_t dphi=trig->Phi()-partn->Phi() ;
    //     while(dphi<-TMath::PiOver2())dphi+=TMath::TwoPi() ;     //  while(dphi>3*TMath::PiOver2())dphi-=TMath::TwoPi() ;

    if(dphi <= -TMath::PiOver2()) dphi+=TMath::TwoPi();

    if(dphi > 3*TMath::PiOver2()) dphi-=TMath::TwoPi();


    // if( dphi > 2/3*TMath::Pi() &&  dphi < 5/3*TMath::Pi() ) {

    Double_t xe=-partn->Pt()*TMath::Cos(dphi)/trig->Pt() ;     
    Double_t ye=-partn->Pt()*TMath::Sin(dphi)/trig->Pt() ;  

    xesum += xe ;  // if( xe > 0 ) xesumpos += xe;  else if( xe < 0 ) xesumneg +=  fabs( xe };
    yesum += ye ;

    Double_t ptTrig = trig->Pt();
    Double_t pttr = partn->Pt();

    Double_t cosdelphi = TMath::Cos(dphi) ;
    if( cosdelphi == 0 ) cosdelphi = 1.e-6;
    Double_t rat= -pttr/ptTrig * cosdelphi / fabs( cosdelphi )  ;

    if ( dphi > 0.66*TMath::Pi() &&  dphi < 1.33*TMath::Pi() && pttr > 0.3   ){
  
         if( (ptTrig > 6 && ptTrig <= 7  )   )  hnum  = itype2 ;
    else if( (ptTrig > 7 && ptTrig <= 8  )   )  hnum  = itype2 + 1 ;
    else if( (ptTrig > 8 && ptTrig <= 9  )   )  hnum  = itype2 + 2 ;
    else if( (ptTrig > 9 && ptTrig <= 10 )   )  hnum  = itype2 + 3 ;
    else if( (ptTrig > 10 && ptTrig <= 12)   )  hnum  = itype2 + 4 ;
    else if( (ptTrig > 12 && ptTrig <= 16 )  )  hnum  = itype2 + 5 ;
    else if( (ptTrig > 16 && ptTrig <= 30 )  )  hnum  = itype2 + 6 ;
    }

    //     printf("itype %d hnum %d trig %f ass %f \n", itype, hnum, ptTrig, pttr );
    //  key0 += hnum ;
    //     key1 += hnum ;
    //     key2 += hnum ;
    // key3 += hnum ;
  
    if( hnum == -1 ) return;

    sprintf(key0,"hdelphi-%d",hnum); sprintf(key1,"hrat-%d",hnum); 
    sprintf(key2,"hxe-%d",hnum); sprintf(key4,"hye-%d",hnum);   sprintf(key3,"hmgg-%d",hnum);
    sprintf(key5,"Multiplicity-%d",hnum);

    // printf("%s %s Multi %f ",key0, key1, fMultiplicity ) ;

    FillHistogram(key0,dphi); 
    FillHistogram(key1,rat); FillHistogram(key2,xe); 
    FillHistogram(key3, trig->M() );   
    FillHistogram(key4,ye);
    FillHistogram(key5,fMultiplicity);
  }
  /*  Int_t  hnum2 = itype2 ;
  sprintf(key0,"hxesum-%d",hnum2);
  FillHistogram( key0,xesum );
  sprintf(key0,"hyesum-%d",hnum2);
  FillHistogram( key0,yesum );
  */
  //  FillHistogram( "hxyesum",xesum,yesum );

 
}
