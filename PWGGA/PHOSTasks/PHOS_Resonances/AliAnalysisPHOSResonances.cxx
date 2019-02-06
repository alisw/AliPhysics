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
 
// Analysis task for resonanse measurement in PHOS
// Authors: Dmitri Peresunko

#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "THashList.h"
#include "TLorentzVector.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisPHOSResonances.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliTriggerAnalysis.h"
#include "AliAODInputHandler.h"
#include "AliPIDResponse.h"


ClassImp(AliAnalysisPHOSResonances) ;

//________________________________________________________________________
AliAnalysisPHOSResonances::AliAnalysisPHOSResonances(const char *name) 
: AliAnalysisTaskSE(name),
  fOutputContainer(0),
  //fTriggerAnalysis(new AliTriggerAnalysis),
  fPIDResponse(0),
  fGamma(0x0), 
  fPi0(0x0), 
  fPi0Merged(0x0), 
  fTracksElm(0x0), 
  fTracksElp(0x0), 
  fTracksPim(0x0), 
  fTracksPip(0x0), 
  fTracksKm(0x0), 
  fTracksKp(0x0), 
  fTracksPm(0x0), 
  fTracksPp(0x0), 
  fLambda(0x0),  
  fMixGamma(0x0),
  fMixPi0(0x0),
  fMixPi0Merged(0x0),
  fMixElm(0x0),
  fMixElp(0x0),
  fMixTracksPim(0x0),
  fMixTracksPip(0x0),
  fMixTracksKm(0x0),
  fMixTracksKp(0x0),
  fMixTracksPm(0x0),
  fMixTracksPp(0x0),
  fMixLambda(0x0),
  fhPr(0x0)
{
  // Constructor
  for(Int_t i=0; i<300; i++){
    fhHistos[i] = 0x0 ;
  }
	
  // Output slots #0 write into a TH1 container
  DefineOutput(1,THashList::Class());

}
AliAnalysisPHOSResonances:: AliAnalysisPHOSResonances(const AliAnalysisPHOSResonances& rh):
  AliAnalysisTaskSE(rh.GetName()),
  fOutputContainer(0x0),
  fPIDResponse(0x0),
  fEvent(0x0),
  fGamma(0x0), 
  fPi0(0x0), 
  fPi0Merged(0x0), 
  fTracksElm(0x0), 
  fTracksElp(0x0), 
  fTracksPim(0x0), 
  fTracksPip(0x0), 
  fTracksKm(0x0), 
  fTracksKp(0x0), 
  fTracksPm(0x0), 
  fTracksPp(0x0), 
  fLambda(0x0),  
  fMixGamma(0x0),
  fMixPi0(0x0),
  fMixPi0Merged(0x0),
  fMixElm(0x0),
  fMixElp(0x0),
  fMixTracksPim(0x0),
  fMixTracksPip(0x0),
  fMixTracksKm(0x0),
  fMixTracksKp(0x0),
  fMixTracksPm(0x0),
  fMixTracksPp(0x0),
  fMixLambda(0x0)
{
  if(fOutputContainer)
    delete fOutputContainer ;  
  fOutputContainer = new THashList() ; 
  fListOfChannels = rh.fListOfChannels ;
  fnPID = rh.fnPID ;
  for(Int_t i=0; i<300; i++){
     if(rh.fhHistos[i]){ 
         fhHistos[i] = new TH3F(*(rh.fhHistos[i])) ;
         fOutputContainer->Add(fhHistos[i]);
     }
     else{
         fhHistos[i] =0x0 ;
     }
  }
  if(rh.fhPr){
    fhPr = new TH2F(*(rh.fhPr)) ;
    fOutputContainer->Add(fhPr) ;
  }
  else{
      fhPr=0x0 ;
  }
      
}

//________________________________________________________________________
AliAnalysisPHOSResonances& AliAnalysisPHOSResonances::operator=(const AliAnalysisPHOSResonances& rh){
 // assignment operator

  this->~AliAnalysisPHOSResonances();
  new(this) AliAnalysisPHOSResonances(rh);
  return *this;
}
//________________________________________________________________________
AliAnalysisPHOSResonances::~AliAnalysisPHOSResonances(){
  //Destructor
  if(fOutputContainer){
     delete fOutputContainer;
     fOutputContainer=0x0 ;
  }
  if(fPIDResponse){
     delete fPIDResponse ;
     fPIDResponse=0x0 ;
  }
 
  if(fGamma){ delete fGamma; fGamma=0x0 ;}
  if(fPi0){delete fPi0; fPi0=0x0 ;} 
  if(fPi0Merged){delete fPi0Merged; fPi0Merged=0x0 ;}
  if(fTracksElm){delete fTracksElm; fTracksElm=0x0;} 
  if(fTracksElp){delete fTracksElp; fTracksElp=0x0;} 
  if(fTracksPim){delete fTracksPim; fTracksPim=0x0;} 
  if(fTracksPip){delete fTracksPip; fTracksPip=0x0;} 
  if(fTracksKm){delete fTracksKm; fTracksKm=0x0;} 
  if(fTracksKp){delete fTracksKp; fTracksKp=0x0;} 
  if(fTracksPm){delete fTracksPm; fTracksPm=0x0;} 
  if(fTracksPp){delete fTracksPp; fTracksPp=0x0;} 
  if(fLambda){delete fLambda; fLambda=0x0;} 
  
  if(fMixGamma){delete fMixGamma; fMixGamma=0x0;} 
  if(fMixPi0){delete fMixPi0; fMixPi0=0x0;} 
  if(fMixPi0Merged){delete fMixPi0Merged; fMixPi0Merged=0x0;} 
  if(fMixElm){delete fMixElm; fMixElm=0x0;} 
  if(fMixElp){delete fMixElp; fMixElp=0x0;} 
  if(fMixTracksPim){delete fMixTracksPim; fMixTracksPim=0x0;} 
  if(fMixTracksPip){delete fMixTracksPip; fMixTracksPip=0x0;} 
  if(fMixTracksKm){delete fMixTracksKm; fMixTracksKm=0x0;} 
  if(fMixTracksKp){delete fMixTracksKp; fMixTracksKp=0x0;} 
  if(fMixTracksPm){delete fMixTracksPm; fMixTracksPm=0x0;} 
  if(fMixTracksPp){delete fMixTracksPp; fMixTracksPp=0x0;} 
 
  if(fMixLambda){delete fMixLambda; fMixLambda=0x0;} 
  //No need to delete histograms in array fhHistos[]!
  //They are deleted as content of fOutputContainer
    
}
//________________________________________________________________________
void AliAnalysisPHOSResonances::UserCreateOutputObjects() {
	// Create histograms
	// Called once

	// AOD histograms
	if(fOutputContainer != NULL){
		delete fOutputContainer;
	}
	fOutputContainer = new THashList();
	fOutputContainer->SetOwner(kTRUE);

        fOutputContainer->Add(new TH1F("hSelEvents","Events selected",10,0.,10.)) ;
        
	fOutputContainer->Add(new TH1F("hCellMultEvent", "PHOS cell multiplicity per event",2000,0,2000));
	fOutputContainer->Add(new TH1F("hClusterMult", "CaloCluster multiplicity", 100,0,100));
	fOutputContainer->Add(new TH1F("hPHOSClusterMult","PHOS cluster multiplicity",100,0,100));
	fOutputContainer->Add(new TH1F("hCellEnergy", "Cell energy", 5000,0.,50.));
	fOutputContainer->Add(new TH1F("hClusterEnergy", "Cluster energy", 5000,0.,50.));
	fOutputContainer->Add(new TH2F("hClusterEvsN", "Cluster energy vs digit multiplicity",     5000,0.,50.,40,0.,40.));
	fOutputContainer->Add(new TH1F("hCellMultClu","Cell multiplicity per cluster",200,0,200));
	fOutputContainer->Add(new TH1F("hModule","Module events",5,0.,5.));
	fOutputContainer->Add(new TH2F("hClusterTOFvsE", "Cluster time vs energy", 500,-250.e-9,250.e-9,40,0.,40.));
        for(Int_t module=1; module<5; module++){
	  fOutputContainer->Add(new TH2F(Form("hModule%d",module), Form("Cluster occupancy in module %d",module), 64,0.,64,56,0.,56.));
        }

	fOutputContainer->Add(new TH1F("hZvertex","Z vertex",200,-50.,+50.));
	fOutputContainer->Add(new TH1F("hNvertexTracks",   "N of primary tracks from the primary vertex",150,0.,150.));
	fOutputContainer->Add(new TH1F("hT0TOF","T0 time (s)",2000,-2.0e-9,2.0e-9));
	
	fOutputContainer->Add(new TH1F("hTrackMult" ,"Charged track multiplicity",150,0.,150.));
	fOutputContainer->Add(new TH1F("hPionMult"  ,"#pi^{#pm} multiplicity"    ,150,0.,150.));
	fOutputContainer->Add(new TH1F("hKaonMult"  ,"K^{#pm} multiplicity"      ,150,0.,150.));
	fOutputContainer->Add(new TH1F("hProtonMult","p,#bar{p} multiplicity"    ,150,0.,150.));
	fOutputContainer->Add(new TH1F("hUndefMult" ,"Undefined multiplicity"    ,150,0.,150.));
	

 	fOutputContainer->Add(new TH2F("hLambdaMass" ,"#Lambda mass",100,0.5,1.5,10,0.,5.));
 	fOutputContainer->Add(new TH2F("hLambdaBarMass" ,"#bar{#Lambda} mass",100,0.5,1.5,10,0.,5.));

  	fOutputContainer->Add(new TH2F("hEp" ,"Ep",100,0.,2.,100,0.,50.));
        
        char cPID[14][15] ;
        fnPID=7; 
        snprintf(cPID[0],15,"All") ;
        snprintf(cPID[1],15,"TOF") ;
        snprintf(cPID[2],15,"Disp");
        snprintf(cPID[3],15,"CPV") ;
        snprintf(cPID[4],15,"Both"); 
        snprintf(cPID[5],15,"NoPi0"); 
        snprintf(cPID[6],15,"NoPi0Eta"); 
        
        
        
       
        //Gamma-gamma
        Int_t nMgg=400;
        Double_t mggMax=2.;
        Int_t nPtgg=80;
        Double_t ptggMax=80.;
        

        for(Int_t i=0; i<300;i++){
          fhHistos[i]=0x0;
        }

        for(Int_t iPID=0; iPID<fnPID; iPID++){
          if(fListOfChannels.At(gg)){  
            fhHistos[gg*fnPID*2+iPID] = new TH3F(Form("hRegg%s",cPID[iPID]) ,"Real (#gamma #gamma)", 1000,0.,4.,40,0.,20.,10,0.,1.);
	    fhHistos[gg*fnPID*2+iPID+fnPID] = new TH3F(Form("hMigg%s",cPID[iPID]) ,"Mixed (#gamma #gamma)",1000,0.,4.,40,0.,20.,10,0.,1.);
          }          
          
          if(fListOfChannels.At(gpp)){  
 	    fhHistos[gpp*fnPID*2+iPID] =new TH3F(Form("hRegpp%s",cPID[iPID]) ,"Real (#gamma p)",nMgg,1.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
	    fhHistos[gpp*fnPID*2+iPID+fnPID] = new TH3F(Form("hMigpp%s",cPID[iPID]) ,"Mixed (#gamma p)",nMgg,1.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
          }
          
          if(fListOfChannels.At(gpm)){  
	    fhHistos[gpm*fnPID*2+iPID] =new TH3F(Form("hRegpm%s",cPID[iPID]) ,"Real (#gamma #bar{p})",nMgg,1.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
	    fhHistos[gpm*fnPID*2+iPID+fnPID] =new TH3F(Form("hMigpm%s",cPID[iPID]) ,"Mixed (#gamma #bar{p})",nMgg,1.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
          }
          
          if(fListOfChannels.At(gpp)){ //DCA comes along gpp  
	    fhHistos[gppDCA*fnPID*2+iPID] =new TH3F(Form("hRegppDCA%s",cPID[iPID]) ,"Real (#gamma p)",nMgg,1.,mggMax,nPtgg,0.,ptggMax,100,0.,10.);
	    fhHistos[gppDCA*fnPID*2+iPID+fnPID] = new TH3F(Form("hMigppDCA%s",cPID[iPID]) ,"Mixed (#gamma p)",nMgg,1.,mggMax,nPtgg,0.,ptggMax,100,0.,10.);
          }
          
          if(fListOfChannels.At(gpm)){ //DCA comes along gpm 
	    fhHistos[gpmDCA*fnPID*2+iPID] =new TH3F(Form("hRegpmDCA%s",cPID[iPID]) ,"Real (#gamma #bar{p})",nMgg,1.,mggMax,nPtgg,0.,ptggMax,100,0.,10.);
	    fhHistos[gpmDCA*fnPID*2+iPID+fnPID] =new TH3F(Form("hMigpmDCA%s",cPID[iPID]) ,"Mixed (#gamma #bar{p})",nMgg,1.,mggMax,nPtgg,0.,ptggMax,100,0.,10.);
          }
          
          
          //Gamma-pi0
          if(fListOfChannels.At(gpi0)){  
	    fhHistos[gpi0*fnPID*2+iPID] =new TH3F(Form("hRegpi0%s",cPID[iPID]) ,"Real (#gamma #pi^{0})",nMgg,0.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
	    fhHistos[gpi0*fnPID*2+iPID+fnPID] =new TH3F(Form("hMigpi0%s",cPID[iPID]) ,"Mixed (#gamma #pi^{0})",nMgg,0.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
          }

          //Gamma-pi0
          if(fListOfChannels.At(gpi0Merg)){  
	    fhHistos[gpi0Merg*fnPID*2+iPID] =new TH3F(Form("hRegpi0Merg%s",cPID[iPID]) ,"Real (#gamma #pi^{0})",nMgg,0.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
	    fhHistos[gpi0Merg*fnPID*2+iPID+fnPID] =new TH3F(Form("hMigpi0Merg%s",cPID[iPID]) ,"Mixed (#gamma #pi^{0})",nMgg,0.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
          }


          //pi0-pi0
          if(fListOfChannels.At(pi0pi0)){  
	    fhHistos[pi0pi0*fnPID*2+iPID] =new TH3F(Form("hRepi0pi0%s",cPID[iPID]) ,"Real (#pi^{0} #pi^{0})",nMgg,0.,mggMax,nPtgg,0.,ptggMax,10,0.,1.);
	    fhHistos[pi0pi0*fnPID*2+iPID+fnPID] =new TH3F(Form("hMipi0pi0%s",cPID[iPID]) ,"Mixed (#pi^{0} #pi^{0})",nMgg,0.,mggMax,nPtgg,0.,ptggMax,10,0.,1.);
          }

          //pi0-pi0Merg
          if(fListOfChannels.At(pi0pi0M)){  
	    fhHistos[pi0pi0M*fnPID*2+iPID] =new TH3F(Form("hRepi0pi0Merg%s",cPID[iPID]) ,"Real (#pi^{0} #pi^{0})",nMgg,0.,mggMax,nPtgg,0.,ptggMax,10,0.,1.);
	    fhHistos[pi0pi0M*fnPID*2+iPID+fnPID] =new TH3F(Form("hMipi0pi0Merg%s",cPID[iPID]) ,"Mixed (#pi^{0} #pi^{0})",nMgg,0.,mggMax,nPtgg,0.,ptggMax,10,0.,1.);
          }
          
          //pi0-p
          if(fListOfChannels.At(pi0pp)){  
	    fhHistos[pi0pp*fnPID*2+iPID] =new TH3F(Form("hRepi0pp%s",cPID[iPID]) ,"Real (#pi^{0} p)",nMgg,1.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
	    fhHistos[pi0pp*fnPID*2+iPID+fnPID] =new TH3F(Form("hMipi0pp%s",cPID[iPID]) ,"Mixed (#pi^{0} p)",nMgg,1.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
          }
          
          if(fListOfChannels.At(pi0pm)){  
	    fhHistos[pi0pm*fnPID*2+iPID] =new TH3F(Form("hRepi0pm%s",cPID[iPID]) ,"Real (#pi^{0} #bar{p})",nMgg,1.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
	    fhHistos[pi0pm*fnPID*2+iPID+fnPID] =new TH3F(Form("hMipi0pm%s",cPID[iPID]) ,"Mixed (#pi^{0} #bar{p})",nMgg,1.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
          }
          
          
          //pi0-p DCA
          if(fListOfChannels.At(pi0pp)){  
	    fhHistos[pi0ppDCA*fnPID*2+iPID] =new TH3F(Form("hRepi0ppDCA%s",cPID[iPID]) ,"Real (#pi^{0} p)",nMgg,1.,mggMax,nPtgg,0.,ptggMax,100,0.,10.);
	    fhHistos[pi0ppDCA*fnPID*2+iPID+fnPID] =new TH3F(Form("hMipi0ppDCA%s",cPID[iPID]) ,"Mixed (#pi^{0} p)",nMgg,1.,mggMax,nPtgg,0.,ptggMax,100,0.,10.);
          }
          
          if(fListOfChannels.At(pi0pm)){  
	    fhHistos[pi0pmDCA*fnPID*2+iPID] =new TH3F(Form("hRepi0pmDCA%s",cPID[iPID]) ,"Real (#pi^{0} #bar{p})",nMgg,1.,mggMax,nPtgg,0.,ptggMax,100,0.,10.);
	    fhHistos[pi0pmDCA*fnPID*2+iPID+fnPID] =new TH3F(Form("hMipi0pmDCA%s",cPID[iPID]) ,"Mixed (#pi^{0} #bar{p})",nMgg,1.,mggMax,nPtgg,0.,ptggMax,100,0.,10.);
          }
          

          //lambda-gamma
          if(fListOfChannels.At(gLam)){  
	    fhHistos[gLam*fnPID*2+iPID] =new TH3F(Form("hRegLam%s",cPID[iPID]) ,"Real (#gamma #Lambda)",nMgg,1.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
	    fhHistos[gLam*fnPID*2+iPID+fnPID] =new TH3F(Form("hMigLam%s",cPID[iPID]) ,"Mixed (#gamma #Lambda)",nMgg,1.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
          }
        
          //lambda-pi0
          if(fListOfChannels.At(pi0Lam)){  
	    fhHistos[pi0Lam*fnPID*2+iPID] =new TH3F(Form("hRepi0Lam%s",cPID[iPID]) ,"Real (#pi^{0} #Lambda)",nMgg,1.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
	    fhHistos[pi0Lam*fnPID*2+iPID+fnPID] =new TH3F(Form("hMipi0Lam%s",cPID[iPID]) ,"Mixed (#pi^{0} #Lambda)",nMgg,1.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
          }
          
          //e+e-
          if(fListOfChannels.At(ee)){  
	    fhHistos[ee*fnPID*2+iPID] =new TH3F(Form("hReee%s",cPID[iPID]) ,"Real (e^{+}e^{-})",400,0.,6,40,0.,50.,10,-1.,1.);
	    fhHistos[ee*fnPID*2+iPID+fnPID] =new TH3F(Form("hMiee%s",cPID[iPID]) ,"Mixed (e^{+}e^{-})",400,0.,6,40,0.,50.,10,-1.,1.);
          }
          
          
          //gamma pi+ pi-
          if(fListOfChannels.At(gpippim)){  
	    fhHistos[gpippim*fnPID*2+iPID] =new TH3F(Form("hRegpippim%s",cPID[iPID]) ,"Real (#gamma#pi^{+}#pi^{-})",nMgg,0.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
	    fhHistos[gpippim*fnPID*2+iPID+fnPID] =new TH3F(Form("hMigpippim%s",cPID[iPID]) ,"Mixed (#gamma#pi^{+}#pi^{-})",nMgg,0.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
          }
          
          //pi0 pi+ pi-
          if(fListOfChannels.At(pi0pippim)){  
	    fhHistos[pi0pippim*fnPID*2+iPID] =new TH3F(Form("hRepi0pippim%s",cPID[iPID]) ,"Real (#pi^{0}#pi^{+}#pi^{-})",nMgg,0.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
	    fhHistos[pi0pippim*fnPID*2+iPID+fnPID] =new TH3F(Form("hMipi0pippim%s",cPID[iPID]) ,"Mixed (#pi^{0}#pi^{+}#pi^{-})",nMgg,0.,mggMax,nPtgg,0.,ptggMax,10,-1.,1.);
          }
         }
        
        
        fhPr = new TH2F("PrDCA","DCA",200,0.,10.,100,0.,10.);
        fOutputContainer->Add(fhPr) ;
        
                 
        for(Int_t i=0; i<300;i++){
          if( fhHistos[i])
          fOutputContainer->Add(fhHistos[i]);
        }


	fMixGamma = new TList() ;
	fMixGamma->SetOwner() ;
        fMixPi0= new TList() ;
	fMixPi0->SetOwner() ;
        fMixPi0Merged= new TList() ;
	fMixPi0Merged->SetOwner() ;

	fMixTracksPp = new TList() ;
	fMixTracksPp->SetOwner() ;
	fMixTracksPm = new TList() ;
	fMixTracksPm->SetOwner() ;
	fMixTracksPip = new TList() ;
	fMixTracksPip->SetOwner() ;
	fMixTracksPim = new TList() ;
	fMixTracksPim->SetOwner() ;
	fMixLambda = new TList() ;
	fMixLambda->SetOwner() ;

	PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisPHOSResonances::UserExec(Option_t *) {
	// Main loop, called for each event
	// Analyze AOD
        
	FillHistogram("hSelEvents",1) ;

	fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
	if (!fEvent) {
		Printf("ERROR: Could not retrieve event");
		return;
	}
	FillHistogram("hSelEvents",2) ;

	AliAODInputHandler *esdH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	if (esdH) 
		fPIDResponse = esdH->GetPIDResponse();
	else {
		Printf("ERROR: Could not get AODInputHandler");
		return;
	}
	FillHistogram("hSelEvents",3) ;

	// Checks if we have a primary vertex
	// Get primary vertices form AOD
/*
	if      (fEvent->GetPrimaryVertexTracks()->GetNContributors()>0)
		fEventVtxExist    = kTRUE;
	else if (fEvent->GetPrimaryVertexSPD()   ->GetNContributors()>0)
		fEventVtxExist    = kTRUE;
*/
	const AliAODVertex *esdVertex5 = fEvent->GetPrimaryVertex();

	Double_t vtx5[3];
	vtx5[0] = esdVertex5->GetX();
	vtx5[1] = esdVertex5->GetY();
	vtx5[2] = esdVertex5->GetZ();
        TVector3 vtx(vtx5) ;
	
	FillHistogram("hNvertexTracks",esdVertex5->GetNContributors());
	FillHistogram("hZvertex"      ,esdVertex5->GetZ());
	if (TMath::Abs(esdVertex5->GetZ()) > 10. )
          return ;

	FillHistogram("hSelEvents",4) ;

	if (fEvent->IsPileupFromSPD())
          return ;

	FillHistogram("hSelEvents",5) ;

	// Fill event statistics for different selection criteria

	//Vtx class z-bin
	Int_t zvtx = (Int_t)((vtx5[2]+10.)/2.) ;
	if(zvtx<0)zvtx=0 ;
	if(zvtx>9)zvtx=9 ;

        if(fPi0)fPi0->Clear() ;
        else fPi0 = new TClonesArray("AliCaloPhoton",100) ;
        if(fPi0Merged)fPi0Merged->Clear() ;
        else fPi0Merged = new TClonesArray("AliCaloPhoton",10) ;
        
        //
        SelectHadrons() ;
        if(fListOfChannels.At(gLam) || fListOfChannels.At(pi0Lam)){ 
          SelectLambda() ;
        }
        SelectGamma() ;
        if(fListOfChannels.At(ee)){  
          SelectElectrons() ;  
        }
        
        
        //================REALs===========
        const Double_t massPi0=0.1349770;
        
        //Fill gamma-gamma and select pi0s
        Int_t nGamma = 0;
        if(fGamma) nGamma = fGamma->GetEntriesFast() ;
        Int_t nPi0Merged = 0;
        if(fPi0Merged) nPi0Merged = fPi0Merged->GetEntriesFast() ;
        Int_t nPp = 0;
        if(fTracksPp) nPp = fTracksPp->GetEntriesFast() ;
        Int_t nPm = 0;
        if(fTracksPm) nPm = fTracksPm->GetEntriesFast() ;
        Int_t nPip = 0;
        if(fTracksPip)nPip = fTracksPip->GetEntriesFast() ;
        Int_t nPim = 0 ;
        if(fTracksPim) nPim = fTracksPim->GetEntriesFast() ;
        Int_t nLambda=0;
        if(fLambda) nLambda=fLambda->GetEntriesFast() ;
        Int_t nPi0=0;
        Int_t nEl=0;
        if(fTracksElm) nEl=fTracksElm->GetEntriesFast() ;
        Int_t nPo=0;
        if(fTracksElp) nPo=fTracksElp->GetEntriesFast() ;
        
        TLorentzVector pair;
        
        //PID cuts: Dispersion, CPV, time, combinations
        //Anti-CPV (R<1)+ time
        for(Int_t i=0; i<nGamma; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
          
          for(Int_t j=i+1; j<nGamma; j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)fGamma->At(j) ;
            pair=*pv1 + *pv2 ;
 	    Double_t alpha = TMath::Abs(pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
            Double_t m=pair.M();
            Double_t pT=pair.Pt() ;
            FillHistogram(0,gg,All,m,pT,alpha) ;
            if(!pv1->IsTOFOK()|| !pv2->IsTOFOK())continue ;
            
            FillHistogram(0,gg,TOF,m,pT,alpha) ;
            
            if(pv1->IsDispOK()&& pv2->IsDispOK()){
              FillHistogram(0,gg,Disp,m,pT,alpha) ;
              if(pv1->IsCPVOK()&& pv2->IsCPVOK()){
                FillHistogram(0,gg,Both,m,pT,alpha) ;
              }
            }
            if(pv1->IsCPVOK()&& pv2->IsCPVOK()){
              FillHistogram(0,gg,CPV,m,pT,alpha) ;
            }
            
            
            if(TMath::Abs(m-Pi0Mass(pT))<Pi0Width(pT)){ //pi0 candidate
//               AliCaloPhoton * pi0 = new((*fPi0)[nPi0++]) AliCaloPhoton(pair.Px(),pair.Py(),pair.Pz(),pair.E()) ;
              AliCaloPhoton * pi0 = new((*fPi0)[nPi0++]) AliCaloPhoton(pair.Px(),pair.Py(),pair.Pz(),
                                                                       TMath::Sqrt(pair.P()*pair.P()+massPi0*massPi0) ) ;
              //Remember photons
              pi0->SetPrimary(i) ;
              pi0->SetPrimaryAtVertex(j) ;
              pi0->SetDispBit(pv1->IsDispOK()&& pv2->IsDispOK()) ;
              pi0->SetCPVBit(pv1->IsCPVOK()&& pv2->IsCPVOK());
              pi0->SetWeight(alpha) ;
              //Mark photons
              pv1->SetPrimary(111) ;
              pv2->SetPrimary(111) ;
            }
            if(TMath::Abs(m-EtaMass(pT))<EtaWidth(pT)){ //Eta candidate
              if(pv1->GetPrimary()==0)   
                pv1->SetPrimary(221) ;
              if(pv2->GetPrimary()==0)   
                pv2->SetPrimary(221) ;
            }
          }
        }
        
        //repeat without tagged photons
        for(Int_t i=0; i<nGamma; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
          if(pv1->GetPrimary()==111 ||!pv1->IsTOFOK() || !pv1->IsDispOK() || !pv1->IsCPVOK())
             continue ; 
          for(Int_t j=i+1; j<nGamma; j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)fGamma->At(j) ;
            if(pv2->GetPrimary()==111 || !pv2->IsTOFOK() ||!pv2->IsDispOK() || !pv2->IsCPVOK())
              continue ; 
            pair=*pv1 + *pv2 ;
 	    Double_t alpha = TMath::Abs(pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
            Double_t m=pair.M();
            Double_t pT=pair.Pt() ;
            FillHistogram(0,gg,NoPi0,m,pT,alpha) ;
            if(pv1->GetPrimary()==0 && pv2->GetPrimary()==0)
              FillHistogram(0,gg,NoPi0Eta,m,pT,alpha) ;
          }
        }
        
        //Fill gamma-pi0 
        if(fListOfChannels.At(gpi0)){  
        for(Int_t i=0; i<nGamma; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
          if(!pv1->IsTOFOK())continue ;
          for(Int_t j=0; j<nPi0; j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)fPi0->At(j) ;
            if(i!=pv2->GetPrimary() && i!=pv2->GetPrimaryAtVertex()){ //photon not from pi0
              pair=*pv1 + *pv2 ;
 	      Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
              Double_t m=pair.M();
              Double_t pT=pair.Pt() ;
              FillHistogram(0,gpi0,All,m,pT,alpha) ;
              if(pv2->GetWeight()>0.4)continue ;
              FillHistogram(0,gpi0,TOF,m,pT,alpha) ;
              if(pv1->IsDispOK()&& pv2->IsDispOK()){
                FillHistogram(0,gpi0,Disp,m,pT,alpha) ;
                if(pv1->IsCPVOK()&& pv2->IsCPVOK()){
                   FillHistogram(0,gpi0,Both,m,pT,alpha) ;
                   if(pv1->GetPrimary()!=111){
                     FillHistogram(0,gpi0,NoPi0,m,pT,alpha) ;
                     if(pv1->GetPrimary()==0){
                       FillHistogram(0,gpi0,NoPi0Eta,m,pT,alpha) ;
                     }
                   }
                }
              }
              if(pv1->IsCPVOK()&& pv2->IsCPVOK()){
                FillHistogram(0,gpi0,CPV,m,pT,alpha) ;
              }
              
            }
          }
        }
        }

        //Fill gamma-mergedpi0 
        if(fListOfChannels.At(gpi0Merg)){  
        for(Int_t i=0; i<nGamma; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
          if(!pv1->IsTOFOK())continue ;
          for(Int_t j=0; j<nPi0Merged; j++){
              AliCaloPhoton * pv2 = (AliCaloPhoton*)fPi0Merged->At(j) ;
              pair=*pv1 + *pv2 ;
 	      Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
              Double_t m=pair.M();
              Double_t pT=pair.Pt() ;
              FillHistogram(0,gpi0Merg,All,m,pT,alpha) ;
              if(!pv2->IsTOFOK())continue ;
              FillHistogram(0,gpi0Merg,TOF,m,pT,alpha) ;
              if(pv1->IsDispOK()&& pv2->IsDispOK()){
                FillHistogram(0,gpi0Merg,Disp,m,pT,alpha) ;
                if(pv1->IsCPVOK()&& pv2->IsCPVOK()){
                   FillHistogram(0,gpi0Merg,Both,m,pT,alpha) ;
                   if(pv1->GetPrimary()!=111){
                     FillHistogram(0,gpi0Merg,NoPi0,m,pT,alpha) ;
                     if(pv1->GetPrimary()==0){
                       FillHistogram(0,gpi0Merg,NoPi0Eta,m,pT,alpha) ;
                     }
                   }
                }
              }
              if(pv1->IsCPVOK()&& pv2->IsCPVOK()){
                FillHistogram(0,gpi0Merg,CPV,m,pT,alpha) ;
              }
          }
        }
        }
        
        //Fill pi0-pi0 
        if(fListOfChannels.At(pi0pi0)){  
        for(Int_t i=0; i<nPi0-1; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fPi0->At(i) ;
          if(pv1->GetWeight()>0.6)continue ; //alpha cut
          for(Int_t j=i+1; j<nPi0; j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)fPi0->At(j) ;
            if(pv2->GetWeight()>0.6)continue ; //alpha cut
            pair=*pv1 + *pv2 ;
            Double_t alpha = TMath::Abs(pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
            Double_t m=pair.M();
            Double_t pT=pair.Pt() ;
            FillHistogram(0,pi0pi0,All,m,pT,alpha) ;
            if(pv1->IsDispOK()&& pv2->IsDispOK()){
              FillHistogram(0,pi0pi0,Disp,m,pT,alpha) ;
              if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
                FillHistogram(0,pi0pi0,Both,m,pT,alpha) ;
              }
            }
            if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
              FillHistogram(0,pi0pi0,CPV,m,pT,alpha) ;
            }
          }
        }
        }
            
        //Fill pi0M-pi0M 
        if(fListOfChannels.At(pi0pi0M)){  
        for(Int_t i=0; i<nPi0; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fPi0->At(i) ;
          for(Int_t j=0; j<nPi0Merged; j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)fPi0Merged->At(j) ;
            pair=*pv1 + *pv2 ;
            Double_t alpha = TMath::Abs(pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
            Double_t m=pair.M();
            Double_t pT=pair.Pt() ;
            FillHistogram(0,pi0pi0M,All,m,pT,alpha) ;
            if(pv1->IsDispOK()&&pv2->IsDispOK()){
              FillHistogram(0,pi0pi0M,Disp,m,pT,alpha) ;
              if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
                FillHistogram(0,pi0pi0M,Both,m,pT,alpha) ;
              }
            }
            if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
              FillHistogram(0,pi0pi0M,CPV,m,pT,alpha) ;
            }
          }
        }
        }
            
            
        
        //Fill pi0-p 
        if(fListOfChannels.At(pi0pp)){  
        for(Int_t i=0; i<nPi0; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fPi0->At(i) ;
          if(pv1->GetWeight()>0.4)continue ; //alpha cut
          for(Int_t j=0; j<nPp; j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)fTracksPp->At(j) ;
            pair=*pv1 + *pv2 ;
            Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
            Double_t m=pair.M();
            Double_t pT=pair.Pt() ;
            Double_t dca= pv2->GetWeight() ; 
            FillHistogram(0,pi0pp,All,m,pT,alpha) ;
            FillHistogram(0,pi0ppDCA,All,m,pT,dca) ;
            if(pv1->IsDispOK()){
              FillHistogram(0,pi0pp,Disp,m,pT,alpha) ;
              FillHistogram(0,pi0ppDCA,Disp,m,pT,dca) ;
              if(pv1->IsCPVOK()){
                FillHistogram(0,pi0pp,Both,m,pT,alpha) ;
                FillHistogram(0,pi0ppDCA,Both,m,pT,dca) ;
              }
            }
            if(pv1->IsCPVOK()){
              FillHistogram(0,pi0pp,CPV,m,pT,alpha) ;
              FillHistogram(0,pi0ppDCA,CPV,m,pT,dca) ;
            }
          }
        }
        }

        //Fill pi0-pbar 
        if(fListOfChannels.At(pi0pm)){  
        for(Int_t i=0; i<nPi0; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fPi0->At(i) ;
          if(pv1->GetWeight()>0.4)continue ; //alpha cut
          for(Int_t j=0; j<nPm; j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)fTracksPm->At(j) ;
            pair=*pv1 + *pv2 ;
            Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
            Double_t m=pair.M();
            Double_t pT=pair.Pt() ;
            Double_t dca= pv2->GetWeight() ; 
            FillHistogram(0,pi0pm,All,m,pT,alpha) ;
            FillHistogram(0,pi0pmDCA,All,m,pT,dca) ;
            if(pv1->IsDispOK()){
              FillHistogram(0,pi0pm,Disp,m,pT,alpha) ;
              FillHistogram(0,pi0pmDCA,Disp,m,pT,dca) ;
              if(pv1->IsCPVOK()){
                FillHistogram(0,pi0pm,Both,m,pT,alpha) ;
                FillHistogram(0,pi0pmDCA,Both,m,pT,dca) ;
              }
            }
            if(pv1->IsCPVOK()){
              FillHistogram(0,pi0pm,CPV,m,pT,alpha) ;
              FillHistogram(0,pi0pmDCA,CPV,m,pT,dca) ;
            }
          }
        }
        }
        
        
        //Fill gamma-p 
        if(fListOfChannels.At(gpp)){  
        for(Int_t i=0; i<nGamma; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
          for(Int_t j=0; j<nPp; j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)fTracksPp->At(j) ;
            pair=*pv1 + *pv2 ;
 	    Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
            Double_t m=pair.M();
            Double_t pT=pair.Pt() ;
            Double_t dca= pv2->GetWeight() ; 
            FillHistogram(0,gpp,All,m,pT,alpha) ;
            FillHistogram(0,gppDCA,All,m,pT,dca) ;
            if(pv1->IsDispOK()){
              FillHistogram(0,gpp,Disp,m,pT,alpha) ;
              FillHistogram(0,gppDCA,Disp,m,pT,dca) ;
              if(pv1->IsCPVOK()){
                FillHistogram(0,gpp,Both,m,pT,alpha) ;
                FillHistogram(0,gppDCA,Both,m,pT,dca) ;
                if(pv1->GetPrimary()!=111){
                  FillHistogram(0,gpp,NoPi0,m,pT,alpha) ;
                  FillHistogram(0,gppDCA,NoPi0,m,pT,dca) ;
                  if(pv1->GetPrimary()==0){
                    FillHistogram(0,gpp,NoPi0Eta,m,pT,alpha) ;
                    FillHistogram(0,gppDCA,NoPi0Eta,m,pT,dca) ;
                  }
                }
              }
            }
            if(pv1->IsCPVOK()){
              FillHistogram(0,gpp,CPV,m,pT,alpha) ;
              FillHistogram(0,gppDCA,CPV,m,pT,dca) ;
            }            
          }
        }
        }

        //Fill gamma-pbar 
        if(fListOfChannels.At(gpm)){  
        for(Int_t i=0; i<nGamma; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
          for(Int_t j=0; j<nPm; j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)fTracksPm->At(j) ;
            pair=*pv1 + *pv2 ;
 	    Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
            Double_t m=pair.M();
            Double_t pT=pair.Pt() ;
            Double_t dca= pv2->GetWeight() ; 
            FillHistogram(0,gpm,All,m,pT,alpha) ;
            FillHistogram(0,gpmDCA,All,m,pT,dca) ;
            if(pv1->IsDispOK()){
              FillHistogram(0,gpm,Disp,m,pT,alpha) ;
              FillHistogram(0,gpmDCA,Disp,m,pT,dca) ;
              if(pv1->IsCPVOK()){
                FillHistogram(0,gpm,Both,m,pT,alpha) ;
                FillHistogram(0,gpm,Both,m,pT,alpha) ;
                if(pv1->GetPrimary()!=111){
                  FillHistogram(0,gpm,NoPi0,m,pT,alpha) ;
                  FillHistogram(0,gpmDCA,NoPi0,m,pT,dca) ;
                  if(pv1->GetPrimary()==0){
                    FillHistogram(0,gpm,NoPi0Eta,m,pT,alpha) ;
                    FillHistogram(0,gpmDCA,NoPi0Eta,m,pT,dca) ;
                  }
                }
              }
            }
            if(pv1->IsCPVOK()){
              FillHistogram(0,gpm,CPV,m,pT,alpha) ;
              FillHistogram(0,gpmDCA,CPV,m,pT,dca) ;
            }            
          }
        }
        }
        
         
        //Fill gamma-Lambda 
        if(fListOfChannels.At(gLam)){  
        for(Int_t i=0; i<nGamma; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
          for(Int_t j=0; j<nLambda; j++){
            TLorentzVector * pv2 = (TLorentzVector*)fLambda->At(j) ;
            pair=*pv1 + *pv2 ;
 	    Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
            Double_t m=pair.M();
            Double_t pT=pair.Pt() ;
            FillHistogram(0,gLam,All,m,pT,alpha) ;
            if(!pv1->IsTOFOK())continue ;
            FillHistogram(0,gLam,TOF,m,pT,alpha) ;
            if(pv1->IsDispOK()){
              FillHistogram(0,gLam,Disp,m,pT,alpha) ;
              if(pv1->IsCPVOK()){
                FillHistogram(0,gLam,Both,m,pT,alpha) ;
                if(pv1->GetPrimary()!=111){
                  FillHistogram(0,gLam,NoPi0,m,pT,alpha) ;
                  if(pv1->GetPrimary()==0){
                    FillHistogram(0,gLam,NoPi0Eta,m,pT,alpha) ;
                  }
                }
              }
            }
            if(pv1->IsCPVOK()){
              FillHistogram(0,gLam,CPV,m,pT,alpha) ;
            }
          }
        }
        }
       
         //Fill pi0-Lambda 
        if(fListOfChannels.At(pi0Lam)){  
        for(Int_t i=0; i<nPi0; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fPi0->At(i) ;
          if(pv1->GetWeight()>0.4)continue ; //alpha cut
          for(Int_t j=0; j<nLambda; j++){
            TLorentzVector * pv2 = (TLorentzVector*)fLambda->At(j) ;
            pair=*pv1 + *pv2 ;
 	    Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
            Double_t m=pair.M();
            Double_t pT=pair.Pt() ;
            FillHistogram(0,pi0Lam,All,m,pT,alpha) ;
            if(pv1->IsDispOK()){
              FillHistogram(0,pi0Lam,Disp,m,pT,alpha) ;
              if(pv1->IsCPVOK()){
                FillHistogram(0,pi0Lam,Both,m,pT,alpha) ;
              }
            }
            if(pv1->IsCPVOK()){
              FillHistogram(0,pi0Lam,CPV,m,pT,alpha) ;
            }
          }
        }
        }
      
        //Fill ee 
        if(fListOfChannels.At(ee)){  
        for(Int_t i=0; i<nEl; i++){
          AliCaloPhoton * pv1 = (AliCaloPhoton*)fTracksElm->At(i) ;
          for(Int_t j=0; j<nPo; j++){
            TLorentzVector * pv2 = (TLorentzVector*)fTracksElp->At(j) ;
            pair=*pv1 + *pv2 ;
 	    Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
            Double_t m=pair.M();
            Double_t pT=pair.Pt() ;
            FillHistogram(0,ee,All,m,pT,alpha) ;
          }
        }
        }

        
        //Fill gamma pi+ pi- 
        if(fListOfChannels.At(gpippim)){  
          for(Int_t i=0; i<nGamma; i++){
            AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
            for(Int_t j=0; j<nPip; j++){
              TLorentzVector * pv2 = (TLorentzVector*)fTracksPip->At(j) ;
              for(Int_t k=0; k<nPim; k++){
                TLorentzVector * pv3 = (TLorentzVector*)fTracksPim->At(k) ;
                pair=*pv1 + *pv2 + *pv3 ;
 	        Double_t alpha2 = (pv1->E()-pv2->E()-pv3->E())/(pv1->E()+pv2->E()+pv3->E()) ;
 	        Double_t alpha = (pv2->E()-pv3->E())/(pv2->E()+pv3->E()) ;
                Double_t m=pair.M();
                Double_t pT=pair.Pt() ;
                FillHistogram(0,gpippim,All,m,pT,alpha) ;
                if(!pv1->IsTOFOK())continue ;
                FillHistogram(0,gpippim,TOF,m,pT,alpha) ;
                if(pv1->IsDispOK()){
                  FillHistogram(0,gpippim,Disp,m,pT,alpha) ;
                  if(pv1->IsCPVOK()){
                    FillHistogram(0,gpippim,Both,m,pT,alpha) ;
                    FillHistogram(0,gpippim,NoPi0Eta,m,pT,alpha2) ;
                    if(pv1->GetPrimary()!=111){
                      FillHistogram(0,gpippim,NoPi0,m,pT,alpha) ;
//                       if(pv1->GetPrimary()==0){
//                         FillHistogram(0,gpippim,NoPi0Eta,m,pT,alpha) ;
//                       }
                    }
                  }
                }
                if(pv1->IsCPVOK()){
                  FillHistogram(0,gpippim,CPV,m,pT,alpha) ;
                }
              }
            }
          }
        }
 
 
        //Fill pi0 pi+ pi- 
        if(fListOfChannels.At(pi0pippim)){  
          for(Int_t i=0; i<nGamma; i++){
            AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
            for(Int_t j=0; j<nPip; j++){
              TLorentzVector * pv2 = (TLorentzVector*)fTracksPip->At(j) ;
              for(Int_t k=0; k<nPim; k++){
                TLorentzVector * pv3 = (TLorentzVector*)fTracksPim->At(k) ;
                pair=*pv1 + *pv2 + *pv3 ;
 	        Double_t alpha2 = (pv1->E()-pv2->E()-pv3->E())/(pv1->E()+pv2->E()+pv3->E()) ;
 	        Double_t alpha = (pv2->E()-pv3->E())/(pv2->E()+pv3->E()) ;
                Double_t m=pair.M();
                Double_t pT=pair.Pt() ;
                FillHistogram(0,pi0pippim,All,m,pT,alpha) ;
                if(!pv1->IsTOFOK())continue ;
                FillHistogram(0,pi0pippim,TOF,m,pT,alpha) ;
                if(pv1->IsDispOK()){
                  FillHistogram(0,pi0pippim,Disp,m,pT,alpha) ;
                  if(pv1->IsCPVOK()){
                    FillHistogram(0,pi0pippim,Both,m,pT,alpha) ;
                    FillHistogram(0,pi0pippim,NoPi0Eta,m,pT,alpha2) ;
                    if(pv1->GetPrimary()!=111){
                      FillHistogram(0,pi0pippim,NoPi0,m,pT,alpha) ;
//                       if(pv1->GetPrimary()==0){
//                         FillHistogram(0,pi0pippim,NoPi0Eta,m,pT,alpha) ;
//                       }
                    }
                  }
                }
                if(pv1->IsCPVOK()){
                  FillHistogram(0,pi0pippim,CPV,m,pT,alpha) ;
                }
              }
            }
          }
        }
 
        //================MIXEDs==================================

       //Fill gamma-gamma 
       if(fListOfChannels.At(gg)){  
       for(Int_t m=0; m<fMixGamma->GetSize(); m++){
         TClonesArray * tmp = (TClonesArray*)fMixGamma->At(m);		
         for(Int_t j=0; j<tmp->GetEntriesFast(); j++){
           AliCaloPhoton * pv2 = (AliCaloPhoton*)tmp->At(j) ;
           for(Int_t i=0; i<nGamma; i++){
               AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
               pair=*pv1 + *pv2 ;
 	       Double_t alpha = TMath::Abs(pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
               Double_t m=pair.M();
               Double_t pT=pair.Pt() ;
               FillHistogram(1,gg,All,m,pT,alpha) ;
               if(!pv1->IsTOFOK()|| !pv2->IsTOFOK())continue ;
               FillHistogram(1,gg,TOF,m,pT,alpha) ;
               
               if(pv1->IsDispOK()){
                 FillHistogram(1,gg,Disp,m,pT,alpha) ;
                 if(pv1->IsCPVOK()){
                   FillHistogram(1,gg,Both,m,pT,alpha) ;
                   if(pv1->GetPrimary()!=111 && pv2->GetPrimary()!=111){
                     FillHistogram(1,gg,NoPi0,m,pT,alpha) ;
                     if(pv1->GetPrimary()==0 && pv2->GetPrimary()==0){
                       FillHistogram(1,gg,NoPi0Eta,m,pT,alpha) ;
                     }
                   }
                 }
               }
               if(pv1->IsCPVOK()){
                 FillHistogram(1,gg,CPV,m,pT,alpha) ;
               }
             }
           }
        }
        }
       
        //Fill gamma-pi0 
       if(fListOfChannels.At(gpi0)){  
        for(Int_t m=0; m<fMixPi0->GetSize(); m++){
           TClonesArray * tmp = (TClonesArray*)fMixPi0->At(m);		
           for(Int_t j=0; j<tmp->GetEntriesFast(); j++){
             AliCaloPhoton * pv2 = (AliCaloPhoton*)tmp->At(j) ;
             for(Int_t i=0; i<nGamma; i++){
               AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
               if(!pv1->IsTOFOK())continue ;
               pair=*pv1 + *pv2 ;
               Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
               Double_t m=pair.M();
               Double_t pT=pair.Pt() ;
               FillHistogram(1,gpi0,All,m,pT,alpha) ;
               if(pv2->GetWeight()>0.4)continue ;
               FillHistogram(1,gpi0,TOF,m,pT,alpha) ;
               if(pv1->IsDispOK()&&pv2->IsDispOK()){
                 FillHistogram(1,gpi0,Disp,m,pT,alpha) ;
                 if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
                   FillHistogram(1,gpi0,Both,m,pT,alpha) ;
                   if(pv1->GetPrimary()!=111){
                     FillHistogram(1,gpi0,NoPi0,m,pT,alpha) ;
                     if(pv1->GetPrimary()==0){
                       FillHistogram(1,gpi0,NoPi0Eta,m,pT,alpha) ;
                     }
                   }
                   
                 }
               }
               if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
                 FillHistogram(1,gpi0,CPV,m,pT,alpha) ;
               }
             }
          }
        }
        for(Int_t m=0; m<fMixGamma->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixGamma->At(m);		
          for(Int_t i=0; i<tmp->GetEntriesFast(); i++){
             AliCaloPhoton * pv1 = (AliCaloPhoton*)tmp->At(i) ;
             if(!pv1->IsTOFOK())continue ;
             for(Int_t j=0; j<nPi0; j++){
               AliCaloPhoton * pv2 = (AliCaloPhoton*)fPi0->At(j) ;
               pair=*pv1 + *pv2 ;
               Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
               Double_t m=pair.M();
               Double_t pT=pair.Pt() ;
               FillHistogram(1,gpi0,All,m,pT,alpha) ;
               if(pv2->GetWeight()>0.4)continue ; //alpha for pi0
               FillHistogram(1,gpi0,TOF,m,pT,alpha) ;
               if(pv1->IsDispOK()&&pv2->IsDispOK()){
                 FillHistogram(1,gpi0,Disp,m,pT,alpha) ;
                 if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
                   FillHistogram(1,gpi0,Both,m,pT,alpha) ;
                   if(pv1->GetPrimary()!=111){
                     FillHistogram(1,gpi0,NoPi0,m,pT,alpha) ;
                     if(pv1->GetPrimary()==0){
                       FillHistogram(1,gpi0,NoPi0Eta,m,pT,alpha) ;
                     }
                   }
                   
                 }
               }
               if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
                 FillHistogram(1,gpi0,CPV,m,pT,alpha) ;
               }
             }
          }
        }
        }
        
        //Fill gamma-pi0 
        if(fListOfChannels.At(gpi0Merg)){  
        for(Int_t m=0; m<fMixPi0Merged->GetSize(); m++){
           TClonesArray * tmp = (TClonesArray*)fMixPi0Merged->At(m);		
           for(Int_t j=0; j<tmp->GetEntriesFast(); j++){
             AliCaloPhoton * pv2 = (AliCaloPhoton*)tmp->At(j) ;
             for(Int_t i=0; i<nGamma; i++){
               AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
               if(!pv1->IsTOFOK()||!pv2->IsTOFOK())continue ;
               pair=*pv1 + *pv2 ;
               Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
               Double_t m=pair.M();
               Double_t pT=pair.Pt() ;
               FillHistogram(1,gpi0Merg,All,m,pT,alpha) ;
//                if(pv2->GetWeight()>0.4)continue ;
               FillHistogram(1,gpi0Merg,TOF,m,pT,alpha) ;
               if(pv1->IsDispOK()&&pv2->IsDispOK()){
                 FillHistogram(1,gpi0Merg,Disp,m,pT,alpha) ;
                 if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
                   FillHistogram(1,gpi0Merg,Both,m,pT,alpha) ;
                   if(pv1->GetPrimary()!=111){
                     FillHistogram(1,gpi0Merg,NoPi0,m,pT,alpha) ;
                     if(pv1->GetPrimary()==0){
                       FillHistogram(1,gpi0Merg,NoPi0Eta,m,pT,alpha) ;
                     }
                   }
                   
                 }
               }
               if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
                 FillHistogram(1,gpi0Merg,CPV,m,pT,alpha) ;
               }
             }
          }
        }
        for(Int_t m=0; m<fMixGamma->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixGamma->At(m);		
          for(Int_t i=0; i<tmp->GetEntriesFast(); i++){
             AliCaloPhoton * pv1 = (AliCaloPhoton*)tmp->At(i) ;
             if(!pv1->IsTOFOK())continue ;
             for(Int_t j=0; j<nPi0Merged; j++){
               AliCaloPhoton * pv2 = (AliCaloPhoton*)fPi0Merged->At(j) ;
               pair=*pv1 + *pv2 ;
               Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
               Double_t m=pair.M();
               Double_t pT=pair.Pt() ;
               FillHistogram(1,gpi0Merg,All,m,pT,alpha) ;
               if(!pv2->IsTOFOK())continue ;
               FillHistogram(1,gpi0Merg,TOF,m,pT,alpha) ;
               if(pv1->IsDispOK()&&pv2->IsDispOK()){
                 FillHistogram(1,gpi0Merg,Disp,m,pT,alpha) ;
                 if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
                   FillHistogram(1,gpi0Merg,Both,m,pT,alpha) ;
                   if(pv1->GetPrimary()!=111){
                     FillHistogram(1,gpi0Merg,NoPi0,m,pT,alpha) ;
                     if(pv1->GetPrimary()==0){
                       FillHistogram(1,gpi0Merg,NoPi0Eta,m,pT,alpha) ;
                     }
                   }
                   
                 }
               }
               if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
                 FillHistogram(1,gpi0Merg,CPV,m,pT,alpha) ;
               }
             }
          }
        }
        }
        
    
       //Fill pi0-pi0 
       if(fListOfChannels.At(pi0pi0)){  
       for(Int_t i=0; i<nPi0; i++){
         AliCaloPhoton * pv1 = (AliCaloPhoton*)fPi0->At(i) ;
         if(pv1->GetWeight()>0.6)continue ;
         for(Int_t m=0; m<fMixPi0->GetSize(); m++){
           TClonesArray * tmp = (TClonesArray*)fMixPi0->At(m);		
           for(Int_t j=0; j<tmp->GetEntriesFast(); j++){
             AliCaloPhoton * pv2 = (AliCaloPhoton*)tmp->At(j) ;
             if(pv2->GetWeight()>0.6)continue ;
             pair=*pv1 + *pv2 ;
             Double_t alpha = TMath::Abs((pv1->E()-pv2->E()))/(pv1->E()+pv2->E()) ;
             Double_t m=pair.M();
             Double_t pT=pair.Pt() ;
             FillHistogram(1,pi0pi0,All,m,pT,alpha) ;
             if(pv1->IsDispOK()&&pv2->IsDispOK()){
               FillHistogram(1,pi0pi0,Disp,m,pT,alpha) ;
               if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
                   FillHistogram(1,pi0pi0,Both,m,pT,alpha) ;
               }
             }
             if(pv1->IsCPVOK()&&pv2->IsCPVOK()){
                 FillHistogram(1,pi0pi0,CPV,m,pT,alpha) ;
             }
           }
         }
       }
       }
        
        //Fill pi0-p 
       if(fListOfChannels.At(pi0pp)){  
        for(Int_t m=0; m<fMixTracksPp->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixTracksPp->At(m);		
          for(Int_t j=0; j<tmp->GetEntriesFast(); j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)tmp->At(j) ;
            for(Int_t i=0; i<nPi0; i++){
              AliCaloPhoton * pv1 = (AliCaloPhoton*)fPi0->At(i) ;
              if(pv1->GetWeight()>0.4)continue ; //alpha cut
              pair=*pv1 + *pv2 ;
              Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
              Double_t m=pair.M();
              Double_t pT=pair.Pt() ;
              Double_t dca= pv2->GetWeight() ; 
              FillHistogram(1,pi0pp,All,m,pT,alpha) ;
              FillHistogram(1,pi0ppDCA,All,m,pT,dca) ;
              if(pv1->IsDispOK()){
                FillHistogram(1,pi0pp,Disp,m,pT,alpha) ;
                FillHistogram(1,pi0ppDCA,Disp,m,pT,dca) ;
                if(pv1->IsCPVOK()){
                  FillHistogram(1,pi0pp,Both,m,pT,alpha) ;
                  FillHistogram(1,pi0ppDCA,Both,m,pT,dca) ;
                }
              }
              if(pv1->IsCPVOK()){
                FillHistogram(1,pi0pp,CPV,m,pT,alpha) ;
                FillHistogram(1,pi0ppDCA,CPV,m,pT,dca) ;
              }
            }
          }
        }
        for(Int_t m=0; m<fMixPi0->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixPi0->At(m);		
          for(Int_t i=0; i<tmp->GetEntriesFast(); i++){
            AliCaloPhoton * pv1 = (AliCaloPhoton*)tmp->At(i) ;
            if(pv1->GetWeight()>0.4)continue ; //alpha cut
            for(Int_t j=0; j<nPp; j++){
              AliCaloPhoton * pv2 = (AliCaloPhoton*)fTracksPp->At(j) ;
              pair=*pv1 + *pv2 ;
              Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
              Double_t m=pair.M();
              Double_t pT=pair.Pt() ;
              Double_t dca= pv2->GetWeight() ; 
              FillHistogram(1,pi0pp,All,m,pT,alpha) ;
              FillHistogram(1,pi0ppDCA,All,m,pT,dca) ;
              if(pv1->IsDispOK()){
                FillHistogram(1,pi0pp,Disp,m,pT,alpha) ;
                FillHistogram(1,pi0ppDCA,Disp,m,pT,dca) ;
                if(pv1->IsCPVOK()){
                  FillHistogram(1,pi0pp,Both,m,pT,alpha) ;
                  FillHistogram(1,pi0ppDCA,Both,m,pT,dca) ;
                }
              }
              if(pv1->IsCPVOK()){
                FillHistogram(1,pi0pp,CPV,m,pT,alpha) ;
                FillHistogram(1,pi0ppDCA,CPV,m,pT,dca) ;
              }
            }
          }
        }
        }

        //Fill pi0-pbar 
       if(fListOfChannels.At(pi0pm)){  
        for(Int_t m=0; m<fMixTracksPm->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixTracksPm->At(m);		
          for(Int_t j=0; j<tmp->GetEntriesFast(); j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)tmp->At(j) ;
            for(Int_t i=0; i<nPi0; i++){
              AliCaloPhoton * pv1 = (AliCaloPhoton*)fPi0->At(i) ;
              if(pv1->GetWeight()>0.4)continue ; //alpha cut
              pair=*pv1 + *pv2 ;
              Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
              Double_t m=pair.M();
              Double_t pT=pair.Pt() ;
              Double_t dca= pv2->GetWeight() ; 
              FillHistogram(1,pi0pm,All,m,pT,alpha) ;
              FillHistogram(1,pi0pmDCA,All,m,pT,dca) ;
              if(pv1->IsDispOK()){
                FillHistogram(1,pi0pm,Disp,m,pT,alpha) ;
                FillHistogram(1,pi0pmDCA,Disp,m,pT,dca) ;
                if(pv1->IsCPVOK()){
                  FillHistogram(1,pi0pm,Both,m,pT,alpha) ;
                  FillHistogram(1,pi0pmDCA,Both,m,pT,dca) ;
                }
              }
              if(pv1->IsCPVOK()){
                FillHistogram(1,pi0pm,CPV,m,pT,alpha) ;
                FillHistogram(1,pi0pmDCA,CPV,m,pT,dca) ;
              }
            }
          }
        }
        for(Int_t m=0; m<fMixPi0->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixPi0->At(m);		
          for(Int_t i=0; i<tmp->GetEntriesFast(); i++){
            AliCaloPhoton * pv1 = (AliCaloPhoton*)tmp->At(i) ;
            if(pv1->GetWeight()>0.4)continue ; //alpha cut
            for(Int_t j=0; j<nPm; j++){
              AliCaloPhoton * pv2 = (AliCaloPhoton*)fTracksPm->At(j) ;
              pair=*pv1 + *pv2 ;
              Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
              Double_t m=pair.M();
              Double_t pT=pair.Pt() ;
              Double_t dca= pv2->GetWeight() ; 
              FillHistogram(1,pi0pm,All,m,pT,alpha) ;
              FillHistogram(1,pi0pmDCA,All,m,pT,dca) ;
              if(pv1->IsDispOK()){
                FillHistogram(1,pi0pm,Disp,m,pT,alpha) ;
                FillHistogram(1,pi0pmDCA,Disp,m,pT,dca) ;
                if(pv1->IsCPVOK()){
                  FillHistogram(1,pi0pm,Both,m,pT,alpha) ;
                  FillHistogram(1,pi0pmDCA,Both,m,pT,dca) ;
                }
              }
              if(pv1->IsCPVOK()){
                FillHistogram(1,pi0pm,CPV,m,pT,alpha) ;
                FillHistogram(1,pi0pmDCA,CPV,m,pT,dca) ;
              }
            }
          }
        }
        }
       
        //Fill gamma-p 
        if(fListOfChannels.At(gpp)){  
        for(Int_t m=0; m<fMixTracksPp->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixTracksPp->At(m);		
          for(Int_t j=0; j<tmp->GetEntriesFast(); j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)tmp->At(j) ;
            for(Int_t i=0; i<nGamma; i++){
              AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
              pair=*pv1 + *pv2 ;
 	      Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
              Double_t m=pair.M();
              Double_t pT=pair.Pt() ;
              Double_t dca= pv2->GetWeight() ; 
              FillHistogram(1,gpp,All,m,pT,alpha) ;
              FillHistogram(1,gppDCA,All,m,pT,dca) ;
              if(!pv1->IsTOFOK())continue ;
              FillHistogram(1,gpp,TOF,m,pT,alpha) ;
              FillHistogram(1,gppDCA,TOF,m,pT,dca) ;
              if(pv1->IsDispOK()){
                FillHistogram(1,gpp,Disp,m,pT,alpha) ;
                FillHistogram(1,gppDCA,Disp,m,pT,dca) ;
                if(pv1->IsCPVOK()){
                  FillHistogram(1,gpp,Both,m,pT,alpha) ;
                  FillHistogram(1,gppDCA,Both,m,pT,alpha) ;
                  if(pv1->GetPrimary()!=111){
                    FillHistogram(1,gpp,NoPi0,m,pT,alpha) ;
                    FillHistogram(1,gppDCA,NoPi0,m,pT,dca) ;
                    if(pv1->GetPrimary()==0){
                      FillHistogram(1,gpp,NoPi0Eta,m,pT,alpha) ;
                      FillHistogram(1,gppDCA,NoPi0Eta,m,pT,dca) ;
                    }
                  }
                }
              }
              if(pv1->IsCPVOK()){
                FillHistogram(1,gpp,CPV,m,pT,alpha) ;
                FillHistogram(1,gppDCA,CPV,m,pT,dca) ;
              }            
            }
          }
        }
        for(Int_t m=0; m<fMixGamma->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixGamma->At(m);		
          for(Int_t i=0; i<tmp->GetEntriesFast(); i++){
             AliCaloPhoton * pv1 = (AliCaloPhoton*)tmp->At(i) ;
             for(Int_t j=0; j<nPp; j++){
               AliCaloPhoton * pv2 = (AliCaloPhoton*)fTracksPp->At(j) ;
               pair=*pv1 + *pv2 ;
 	       Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
               Double_t m=pair.M();
               Double_t pT=pair.Pt() ;
               Double_t dca= pv2->GetWeight() ; 
               FillHistogram(1,gpp,All,m,pT,alpha) ;
               FillHistogram(1,gppDCA,All,m,pT,dca) ;
               if(!pv1->IsTOFOK())continue ;
               FillHistogram(1,gpp,TOF,m,pT,alpha) ;
               FillHistogram(1,gppDCA,TOF,m,pT,dca) ;
               if(pv1->IsDispOK()){
                 FillHistogram(1,gpp,Disp,m,pT,alpha) ;
                 FillHistogram(1,gppDCA,Disp,m,pT,dca) ;
                 if(pv1->IsCPVOK()){
                   FillHistogram(1,gpp,Both,m,pT,alpha) ;
                   FillHistogram(1,gppDCA,Both,m,pT,dca) ;
                   if(pv1->GetPrimary()!=111){
                     FillHistogram(1,gpp,NoPi0,m,pT,alpha) ;
                     FillHistogram(1,gppDCA,NoPi0,m,pT,dca) ;
                     if(pv1->GetPrimary()==0){
                       FillHistogram(1,gpp,NoPi0Eta,m,pT,alpha) ;
                       FillHistogram(1,gppDCA,NoPi0Eta,m,pT,dca) ;
                     }
                   }
                 }
               }
               if(pv1->IsCPVOK()){
                 FillHistogram(1,gpp,CPV,m,pT,alpha) ;
                 FillHistogram(1,gppDCA,CPV,m,pT,dca) ;
               }
             }
          }
        }
        }

        //Fill gamma-pbar 
        if(fListOfChannels.At(gpm)){  
        for(Int_t m=0; m<fMixTracksPm->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixTracksPm->At(m);		
          for(Int_t j=0; j<tmp->GetEntriesFast(); j++){
            AliCaloPhoton * pv2 = (AliCaloPhoton*)tmp->At(j) ;
            for(Int_t i=0; i<nGamma; i++){
              AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
              pair=*pv1 + *pv2 ;
 	      Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
              Double_t m=pair.M();
              Double_t pT=pair.Pt() ;
               Double_t dca= pv2->GetWeight() ; 
               FillHistogram(1,gpm,All,m,pT,alpha) ;
               FillHistogram(1,gpmDCA,All,m,pT,dca) ;
               if(!pv1->IsTOFOK())continue ;
               FillHistogram(1,gpm,TOF,m,pT,alpha) ;
               FillHistogram(1,gpmDCA,TOF,m,pT,dca) ;
               if(pv1->IsDispOK()){
                 FillHistogram(1,gpm,Disp,m,pT,alpha) ;
                 FillHistogram(1,gpmDCA,Disp,m,pT,dca) ;
                 if(pv1->IsCPVOK()){
                   FillHistogram(1,gpm,Both,m,pT,alpha) ;
                   FillHistogram(1,gpmDCA,Both,m,pT,dca) ;
                   if(pv1->GetPrimary()!=111){
                     FillHistogram(1,gpm,NoPi0,m,pT,alpha) ;
                     FillHistogram(1,gpmDCA,NoPi0,m,pT,dca) ;
                     if(pv1->GetPrimary()==0){
                       FillHistogram(1,gpm,NoPi0Eta,m,pT,alpha) ;
                       FillHistogram(1,gpmDCA,NoPi0Eta,m,pT,dca) ;
                     }
                   }
                 }
               }
               if(pv1->IsCPVOK()){
                 FillHistogram(1,gpm,CPV,m,pT,alpha) ;
                 FillHistogram(1,gpmDCA,CPV,m,pT,dca) ;
               }
            }
          }
        }
        for(Int_t m=0; m<fMixGamma->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixGamma->At(m);		
          for(Int_t i=0; i<tmp->GetEntriesFast(); i++){
             AliCaloPhoton * pv1 = (AliCaloPhoton*)tmp->At(i) ;
             for(Int_t j=0; j<nPm; j++){
               AliCaloPhoton * pv2 = (AliCaloPhoton*)fTracksPm->At(j) ;
               pair=*pv1 + *pv2 ;
               Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
               Double_t m=pair.M();
               Double_t pT=pair.Pt() ;
               Double_t dca= pv2->GetWeight() ; 
               FillHistogram(1,gpm,All,m,pT,alpha) ;
               FillHistogram(1,gpmDCA,All,m,pT,dca) ;
               if(!pv1->IsTOFOK())continue ;
               FillHistogram(1,gpm,TOF,m,pT,alpha) ;
               FillHistogram(1,gpmDCA,TOF,m,pT,dca) ;
               if(pv1->IsDispOK()){
                 FillHistogram(1,gpm,Disp,m,pT,alpha) ;
                 FillHistogram(1,gpmDCA,Disp,m,pT,dca) ;
                 if(pv1->IsCPVOK()){
                   FillHistogram(1,gpm,Both,m,pT,alpha) ;
                   FillHistogram(1,gpmDCA,Both,m,pT,dca) ;
                   if(pv1->GetPrimary()!=111){
                     FillHistogram(1,gpm,NoPi0,m,pT,alpha) ;
                     FillHistogram(1,gpmDCA,NoPi0,m,pT,dca) ;
                     if(pv1->GetPrimary()==0){
                       FillHistogram(1,gpm,NoPi0Eta,m,pT,alpha) ;
                       FillHistogram(1,gpmDCA,NoPi0Eta,m,pT,dca) ;
                     }
                   }
                 }
               }
               if(pv1->IsCPVOK()){
                 FillHistogram(1,gpm,CPV,m,pT,alpha) ;
                 FillHistogram(1,gpmDCA,CPV,m,pT,dca) ;
               }
            }
          }
        }
        }
         
        //Fill gamma-Lambda 
        if(fListOfChannels.At(gLam)){  
        for(Int_t m=0; m<fMixLambda->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixLambda->At(m);		
          for(Int_t j=0; j<tmp->GetEntriesFast(); j++){
            TLorentzVector * pv2 = (TLorentzVector*)tmp->At(j) ;
            for(Int_t i=0; i<nGamma; i++){
              AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
              pair=*pv1 + *pv2 ;
              Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
              Double_t m=pair.M();
              Double_t pT=pair.Pt() ;
              FillHistogram(1,gLam,All,m,pT,alpha) ;
              if(!pv1->IsTOFOK())continue ;
              FillHistogram(1,gLam,TOF,m,pT,alpha) ;
              if(pv1->IsDispOK()){
                 FillHistogram(1,gLam,Disp,m,pT,alpha) ;
                 if(pv1->IsCPVOK()){
                   FillHistogram(1,gLam,Both,m,pT,alpha) ;
                  if(pv1->GetPrimary()!=111){
                    FillHistogram(1,gLam,NoPi0,m,pT,alpha) ;
                    if(pv1->GetPrimary()==0){
                      FillHistogram(1,gLam,NoPi0Eta,m,pT,alpha) ;
                    }
                  }
                 }
               }
               if(pv1->IsCPVOK()){
                 FillHistogram(1,gLam,CPV,m,pT,alpha) ;
               }
            }
          }
        }
        for(Int_t m=0; m<fMixGamma->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixGamma->At(m);		
          for(Int_t i=0; i<tmp->GetEntriesFast(); i++){
             AliCaloPhoton * pv1 = (AliCaloPhoton*)tmp->At(i) ;
             for(Int_t j=0; j<nLambda; j++){
               TLorentzVector * pv2 = (TLorentzVector*)fLambda->At(j) ;
               pair=*pv1 + *pv2 ;
 	       Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
               Double_t m=pair.M();
               Double_t pT=pair.Pt() ;
               FillHistogram(1,gLam,All,m,pT,alpha) ;
               if(!pv1->IsTOFOK())continue ;
               FillHistogram(1,gLam,TOF,m,pT,alpha) ;
                if(pv1->IsDispOK()){
                 FillHistogram(1,gLam,Disp,m,pT,alpha) ;
                 if(pv1->IsCPVOK()){
                   FillHistogram(1,gLam,Both,m,pT,alpha) ;
                   if(pv1->GetPrimary()!=111){
                     FillHistogram(1,gLam,NoPi0,m,pT,alpha) ;
                     if(pv1->GetPrimary()==0){
                       FillHistogram(1,gLam,NoPi0Eta,m,pT,alpha) ;
                     }
                   }
                 }
               }
               if(pv1->IsCPVOK()){
                 FillHistogram(1,gLam,CPV,m,pT,alpha) ;
               }
             }
          }
        }
        }
       
        //Fill pi0-Lambda 
        if(fListOfChannels.At(pi0Lam)){  
        for(Int_t m=0; m<fMixLambda->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixLambda->At(m);		
          for(Int_t j=0; j<tmp->GetEntriesFast(); j++){
            TLorentzVector * pv2 = (TLorentzVector*)tmp->At(j) ;
            for(Int_t i=0; i<nPi0; i++){
              AliCaloPhoton * pv1 = (AliCaloPhoton*)fPi0->At(i) ;
              if(pv1->GetWeight()>0.4)continue ; //alpha cut
              pair=*pv1 + *pv2 ;
              Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
              Double_t m=pair.M();
              Double_t pT=pair.Pt() ;
              FillHistogram(1,pi0Lam,All,m,pT,alpha) ;
               if(pv1->IsDispOK()){
                 FillHistogram(1,pi0Lam,Disp,m,pT,alpha) ;
                 if(pv1->IsCPVOK()){
                   FillHistogram(1,pi0Lam,Both,m,pT,alpha) ;
                 }
               }
               if(pv1->IsCPVOK()){
                 FillHistogram(1,pi0Lam,CPV,m,pT,alpha) ;
               }
            }
          }
        }
        for(Int_t m=0; m<fMixPi0->GetSize(); m++){
          TClonesArray * tmp = (TClonesArray*)fMixPi0->At(m);		
          for(Int_t i=0; i<tmp->GetEntriesFast(); i++){
             AliCaloPhoton * pv1 = (AliCaloPhoton*)tmp->At(i) ;
             if(pv1->GetWeight()>0.4)continue ; //alpha cut
             for(Int_t j=0; j<nLambda; j++){
               TLorentzVector * pv2 = (TLorentzVector*)fLambda->At(j) ;
               pair=*pv1 + *pv2 ;
 	       Double_t alpha = (pv1->E()-pv2->E())/(pv1->E()+pv2->E()) ;
               Double_t m=pair.M();
               Double_t pT=pair.Pt() ;
               FillHistogram(1,pi0Lam,All,m,pT,alpha) ;
               if(pv1->IsDispOK()){
                 FillHistogram(1,pi0Lam,Disp,m,pT,alpha) ;
                 if(pv1->IsCPVOK()){
                   FillHistogram(1,pi0Lam,Both,m,pT,alpha) ;
                 }
               }
               if(pv1->IsCPVOK()){
                 FillHistogram(1,pi0Lam,CPV,m,pT,alpha) ;
               }
             }
          }
        }
        }
        
        
        if(fListOfChannels.At(gpippim)){  
          for(Int_t m=0; m<fMixTracksPip->GetSize(); m++){
            TClonesArray * tmp = (TClonesArray*)fMixTracksPip->At(m);	
            for(Int_t j=0; j<tmp->GetEntriesFast(); j++){
              AliCaloPhoton * pv2 = (AliCaloPhoton*)tmp->At(j) ;
              for(Int_t n=0; n<fMixTracksPim->GetSize(); n++){
                if(m==n) continue ;  
                TClonesArray * tmp2 = (TClonesArray*)fMixTracksPim->At(n);	
                for(Int_t k=0; k<tmp2->GetEntriesFast(); k++){
                  AliCaloPhoton * pv3 = (AliCaloPhoton*)tmp2->At(k) ;
                  for(Int_t i=0; i<nGamma; i++){
                    AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;
                    pair=*pv1 + *pv2 + *pv3 ;
 	            Double_t alpha = (pv2->E()-pv3->E())/(pv2->E()+pv3->E()) ;
 	            Double_t alpha2 = (pv1->E()-pv2->E()-pv3->E())/(pv1->E()+pv2->E()+pv3->E()) ;
                    Double_t m=pair.M();
                    Double_t pT=pair.Pt() ;
                    FillHistogram(1,gpippim,All,m,pT,alpha) ;
                    if(!pv1->IsTOFOK())continue ;
                    FillHistogram(1,gpippim,TOF,m,pT,alpha) ;
                    if(pv1->IsDispOK()){
                      FillHistogram(1,gpippim,Disp,m,pT,alpha) ;
                      if(pv1->IsCPVOK()){
                        FillHistogram(1,gpippim,Both,m,pT,alpha) ;
                        FillHistogram(1,gpippim,NoPi0Eta,m,pT,alpha2) ;
                        if(pv1->GetPrimary()!=111){
                          FillHistogram(1,gpippim,NoPi0,m,pT,alpha) ;
                        }
                      }
                    }
                    if(pv1->IsCPVOK()){
                      FillHistogram(1,gpippim,CPV,m,pT,alpha) ;
                    }
                  }
                }
              }            
            }
          }
        }
        
        if(fListOfChannels.At(pi0pippim)){  
          for(Int_t m=0; m<fMixTracksPip->GetSize(); m++){
            TClonesArray * tmp = (TClonesArray*)fMixTracksPip->At(m);	
            for(Int_t j=0; j<tmp->GetEntriesFast(); j++){
              AliCaloPhoton * pv2 = (AliCaloPhoton*)tmp->At(j) ;
              for(Int_t n=0; n<fMixTracksPim->GetSize(); n++){
                if(m==n) continue ;  
                TClonesArray * tmp2 = (TClonesArray*)fMixTracksPim->At(n);	
                for(Int_t k=0; k<tmp2->GetEntriesFast(); k++){
                  AliCaloPhoton * pv3 = (AliCaloPhoton*)tmp2->At(k) ;
                  for(Int_t i=0; i<nPi0; i++){
                    AliCaloPhoton * pv1 = (AliCaloPhoton*)fPi0->At(i) ;
                    pair=*pv1 + *pv2 + *pv3 ;
 	            Double_t alpha = (pv2->E()-pv3->E())/(pv2->E()+pv3->E()) ;
 	            Double_t alpha2 = (pv1->E()-pv2->E()-pv3->E())/(pv1->E()+pv2->E()+pv3->E()) ;
                    Double_t m=pair.M();
                    Double_t pT=pair.Pt() ;
                    FillHistogram(1,pi0pippim,All,m,pT,alpha) ;
                    if(!pv1->IsTOFOK())continue ;
                    FillHistogram(1,pi0pippim,TOF,m,pT,alpha) ;
                    if(pv1->IsDispOK()){
                      FillHistogram(1,pi0pippim,Disp,m,pT,alpha) ;
                      if(pv1->IsCPVOK()){
                        FillHistogram(1,pi0pippim,Both,m,pT,alpha) ;
                        FillHistogram(1,pi0pippim,NoPi0Eta,m,pT,alpha2) ;
                        if(pv1->GetPrimary()!=111){
                          FillHistogram(1,pi0pippim,NoPi0,m,pT,alpha) ;
                        }
                      }
                    }
                    if(pv1->IsCPVOK()){
                      FillHistogram(1,pi0pippim,CPV,m,pT,alpha) ;
                    }
                  }
                }
              }            
            }
          }
            
        }
        
    
	//=====================================================
	//Now we either add current events to stack or remove
	//If no photons in current event - no need to add it to mixed
	const Int_t kMixEvents=100;
	const Int_t kMixEventsHadr=4;
	if(fGamma && fGamma->GetEntriesFast()>0){
		fMixGamma->AddFirst(fGamma) ;
		fGamma=0;
		if(fMixGamma->GetSize()>kMixEvents){//Remove redundant events
			TClonesArray * tmp = static_cast<TClonesArray*>(fMixGamma->Last()) ;
			fMixGamma->RemoveLast() ;
			delete tmp ;
		}
	}
	if(fPi0 && fPi0->GetEntriesFast()>0){
		fMixPi0->AddFirst(fPi0) ;
		fPi0=0;
		if(fMixPi0->GetSize()>kMixEvents){//Remove redundant events
			TClonesArray * tmp = static_cast<TClonesArray*>(fMixPi0->Last()) ;
			fMixPi0->RemoveLast() ;
			delete tmp ;
		}
	}
	if(fPi0Merged && fPi0Merged->GetEntriesFast()>0){
		fMixPi0Merged->AddFirst(fPi0Merged) ;
		fPi0Merged=0;
		if(fMixPi0Merged->GetSize()>kMixEvents){//Remove redundant events
			TClonesArray * tmp = static_cast<TClonesArray*>(fMixPi0Merged->Last()) ;
			fMixPi0Merged->RemoveLast() ;
			delete tmp ;
		}
	}
	
	if(fTracksPp && fTracksPp->GetEntriesFast()>0){
		fMixTracksPp->AddFirst(fTracksPp) ;
			fTracksPp=0;
		if(fMixTracksPp->GetSize()>kMixEventsHadr){//Remove redundant events
			TClonesArray * tmp = static_cast<TClonesArray*>(fMixTracksPp->Last()) ;
			fMixTracksPp->RemoveLast() ;
			delete tmp;
		}
	}
	if(fTracksPm && fTracksPm->GetEntriesFast()>0){
		fMixTracksPm->AddFirst(fTracksPm) ;
			fTracksPm=0;
		if(fMixTracksPm->GetSize()>kMixEventsHadr){//Remove redundant events
			TClonesArray * tmp = static_cast<TClonesArray*>(fMixTracksPm->Last()) ;
			fMixTracksPm->RemoveLast() ;
			delete tmp;
		}
	}

	
	if(fTracksPip && fTracksPip->GetEntriesFast()>0){
		fMixTracksPip->AddFirst(fTracksPip) ;
			fTracksPip=0;
		if(fMixTracksPip->GetSize()>kMixEventsHadr){//Remove redundant events
			TClonesArray * tmp = static_cast<TClonesArray*>(fMixTracksPip->Last()) ;
			fMixTracksPip->RemoveLast() ;
			delete tmp;
		}
	}
	if(fTracksPim && fTracksPim->GetEntriesFast()>0){
		fMixTracksPim->AddFirst(fTracksPim) ;
			fTracksPim=0;
		if(fMixTracksPim->GetSize()>kMixEventsHadr){//Remove redundant events
			TClonesArray * tmp = static_cast<TClonesArray*>(fMixTracksPim->Last()) ;
			fMixTracksPim->RemoveLast() ;
			delete tmp;
		}
	}
	
	
	if(fLambda && fLambda->GetEntriesFast()>0){
		fMixLambda->AddFirst(fLambda) ;
		fLambda=0;
		if(fMixLambda->GetSize()>kMixEventsHadr){//Remove redundant events
			TClonesArray * tmp = static_cast<TClonesArray*>(fMixLambda->Last()) ;
			fMixLambda->RemoveLast() ;
			delete tmp ;
		}
	}


	// Post output data.
	PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisPHOSResonances::Terminate(Option_t *)
{
}
//_____________________________________________________________________________
void AliAnalysisPHOSResonances::FillHistogram(const char * key,Double_t x)const{
	//FillHistogram
	TH1 * hist = dynamic_cast<TH1*>(fOutputContainer->FindObject(key)) ;
	if(hist)
		hist->Fill(x) ;
	else
		AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisPHOSResonances::FillHistogram(const char * key,Double_t x,Double_t y)const{
	//FillHistogram
	TH1 * th1 = dynamic_cast<TH1*> (fOutputContainer->FindObject(key));
	if(th1)
		th1->Fill(x, y) ;
	else
		AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisPHOSResonances::FillHistogram(Int_t remi,Channels_t ch,PID_t pid,Double_t x, Double_t y,Double_t z) const{
    //Fill 2D histogram witn name key
    
    TH3F * tmp = fhHistos[ch*fnPID*2+pid+fnPID*remi] ;
    if(tmp)
      tmp->Fill(x,y,z) ;
    else
      printf("No histogram: Re=%d, channel=%d, PID=%d \n",remi,ch,pid);  
    
}

//_____________________________________________________________________________
void AliAnalysisPHOSResonances::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
	//Fills 1D histograms with key
	TObject * obj = fOutputContainer->FindObject(key);

	TH2 * th2 = dynamic_cast<TH2*> (obj);
	if(th2) {
		th2->Fill(x, y, z) ;
		return;
	}

	TH3 * th3 = dynamic_cast<TH3*> (obj);
	if(th3) {
		th3->Fill(x, y, z) ;
		return;
	}

	AliError(Form("can not find histogram (of instance TH2) <%s> ",key)) ;
}

//_____________________________________________________________________________
void AliAnalysisPHOSResonances::SelectHadrons() {

        
	// Multiplicity and momentum distribution of tracks
	const Double_t massPip = 0.13957;
// 	const Double_t massK = 0.493677;
	const Double_t massP = 0.938272;

	Int_t nTracks = fEvent->GetNumberOfTracks();
	if(!fTracksPip)
		fTracksPip = new TClonesArray("TLorentzVector",fEvent->GetNumberOfTracks()) ;
        else
            fTracksPip->Clear() ;
	if(!fTracksPim)
		fTracksPim = new TClonesArray("TLorentzVector",fEvent->GetNumberOfTracks()) ;
        else
            fTracksPim->Clear() ;
// 	if(!fTracksKp)
// 		fTracksKp = new TClonesArray("TLorentzVector",fEvent->GetNumberOfTracks()) ;
//         else
//             fTracksKp->Clear() ;
// 	if(!fTracksKm)
// 		fTracksKm = new TClonesArray("TLorentzVector",fEvent->GetNumberOfTracks()) ;
//         else
//             fTracksKm->Clear() ;
	if(!fTracksPp)
		fTracksPp = new TClonesArray("AliCaloPhoton",fEvent->GetNumberOfTracks()) ;
        else
            fTracksPp->Clear() ;
	if(!fTracksPm)
		fTracksPm = new TClonesArray("AliCaloPhoton",fEvent->GetNumberOfTracks()) ;
        else
            fTracksPm->Clear() ;
	
	Int_t inPip=0 ;
	Int_t inPim=0 ;
// 	Int_t inKp=0 ;
// 	Int_t inKm=0 ;
	Int_t inPp=0 ;
	Int_t inPm=0 ;

	for (Int_t i=0;i<nTracks;i++) {
            AliAODTrack *track = (AliAODTrack*)fEvent->GetTrack(i) ;
            if(!track->IsHybridGlobalConstrainedGlobal() || TMath::Abs(track->Eta())> 0.8) 
                continue ;
            AliVParticle *inEvHMain = dynamic_cast<AliVParticle*>(track);
            Bool_t pidPion=kFALSE, pidKaon=kFALSE, pidProton=kFALSE, pidUndef=kFALSE;
	    Double_t nsigmaProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton));
	    Double_t nsigmaKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon)); 
	    Double_t nsigmaPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion)); 
   	    // guess the particle based on the smaller nsigma
	    if((nsigmaKaon  <nsigmaPion) && (nsigmaKaon  <nsigmaProton) && (nsigmaKaon  <3)) pidKaon   = kTRUE;
	    if((nsigmaPion  <nsigmaKaon) && (nsigmaPion  <nsigmaProton) && (nsigmaPion  <3)) pidPion   = kTRUE;
	    if((nsigmaProton<nsigmaKaon) && (nsigmaProton<nsigmaPion  ) && (nsigmaProton<3)) pidProton = kTRUE;
	    if (!pidPion && !pidKaon && !pidProton) pidUndef = kTRUE;

   	    if (pidPion) {
              if(track->Charge()>0){
                new((*fTracksPip)[inPip++]) TLorentzVector(track->Px(),track->Py(),track->Pz(),TMath::Sqrt(massPip*massPip+track->P()*track->P())) ;
              }
              else{
                new((*fTracksPim)[inPim++]) TLorentzVector(track->Px(),track->Py(),track->Pz(),TMath::Sqrt(massPip*massPip+track->P()*track->P())) ;
              }
            }
//             if (pidKaon){
//               if(track->Charge()>0){
//                 new((*fTracksKp)[inKp++]) TLorentzVector(track->Px(),track->Py(), track->Pz(), TMath::Sqrt(massK*massK+track->P()*track->P())) ;
//               }
//               else{
//                 new((*fTracksKm)[inKm++]) TLorentzVector(track->Px(),track->Py(), track->Pz(), TMath::Sqrt(massK*massK+track->P()*track->P())) ;
//               }
// 
//             }
            if (pidProton) {
              AliCaloPhoton * pr;  
              if(track->Charge()>0){ 
                pr = new((*fTracksPp)[inPp++]) AliCaloPhoton(track->Px(),track->Py(), track->Pz(), TMath::Sqrt(massP*massP+track->P()*track->P())) ;
              }
              else{
                pr =new((*fTracksPm)[inPm++]) AliCaloPhoton(track->Px(),track->Py(), track->Pz(), TMath::Sqrt(massP*massP+track->P()*track->P())) ;
              }
              Double_t dca=track->DCA() ;
              pr->SetDispBit(dca>0.8) ;
              pr->SetWeight(dca) ;
              fhPr->Fill(dca,track->Pt()) ;
            }
	}
}


//_____________________________________________________________________________
void AliAnalysisPHOSResonances::SelectLambda() {
  //Select Lmbdas from V0
    
    if(!fLambda)
      fLambda = new TClonesArray("TLorentzVector",100) ;
    else
      fLambda->Clear() ;  
    
    Int_t nv0 = fEvent->GetNumberOfV0s();
    Int_t inLambda=0 ;
    const Double_t massLambda=1.115683;
    while (nv0--) {
        AliAODv0 *v0=fEvent->GetV0(nv0);
        if (!v0) {continue;} 
        
        //Use onfly only
        if(!v0->GetOnFlyStatus()) continue ;

        const AliAODTrack *ntrack1=(AliAODTrack *)v0->GetDaughter(1);
        if (!AcceptTrack(ntrack1)) continue;
        
        const AliAODTrack *ptrack1=(AliAODTrack *)v0->GetDaughter(0);
        if (!AcceptTrack(ptrack1)) continue;
        
	// Remove like-sign
	if (ntrack1->Charge() == ptrack1->Charge()){
	  continue;
	} 
            
	if (v0->Pt()==0) 	{continue;}
        
        if( ntrack1->GetKinkIndex(0) > 0 || ptrack1->GetKinkIndex(0) > 0) continue ;        
        
        Double_t  lV0Position[3];
        v0->GetXYZ(lV0Position);
        
        Double_t  lV0Radius      = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
        if(lV0Radius>180.)
           continue ; 
      
        
 /*      // DCA between daughter and Primary Vertex:
       Double_t lDcaPosToPrimVertex=TMath::Abs(pos->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField));
       Double_t lDcaNegToPrimVertex=TMath::Abs(pos->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField));
       if(lDcaPosToPrimVertex<0.06 || lDcaNegToPrimVertex<0.06)
         continue ; 
  */      
        //DCA V0 daughters
        Double_t dca=v0->DcaV0Daughters();
        if (dca<0.06) continue;
        

        Double_t cpa=v0->CosPointingAngle(fEvent->GetPrimaryVertex());
        if (cpa<0.993) continue;

//         Double_t xyz[3]; v0->GetSecondaryVtx(xyz);
//         Double_t r2=xyz[0]*xyz[0] + xyz[1]*xyz[1];
//         if (r2<0.9*0.9) return kFALSE;
//         if (r2>100*100) return kFALSE;

    
 	Double_t  nSigmaPosPion =	TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrack1, AliPID::kPion));
	Double_t  nSigmaNegPion =	TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrack1, AliPID::kPion));
	Double_t  nSigmaPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrack1, AliPID::kProton));
	Double_t  nSigmaNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrack1, AliPID::kProton));
        
        Bool_t isLambda=0;
        Bool_t isLambdaBar=0;
        
        //DCA pi+- >0.02 cm
        //SCA p > 0.05 cm
        Float_t xyNeg=v0->DcaNegToPrimVertex();
        Float_t xyPos=v0->DcaPosToPrimVertex();
        if(ntrack1->Charge()>0){ //Lambda and proton
          isLambda= (nSigmaNegProton<3.7) && (nSigmaPosPion<3.8) && (TMath::Abs(xyPos)>0.02) && (TMath::Abs(xyNeg)>0.05) ;
          isLambdaBar=(nSigmaPosProton<3.9) && (nSigmaNegPion<4.2) && (TMath::Abs(xyNeg)>0.02) && (TMath::Abs(xyPos)>0.05) ;
        }
        else{
          isLambda= (nSigmaPosProton<3.7) && (nSigmaNegPion<3.8) && (TMath::Abs(xyNeg)>0.02) && (TMath::Abs(xyPos)>0.05) ;
          isLambdaBar=(nSigmaNegProton<3.9) && (nSigmaPosPion<4.2) && (TMath::Abs(xyPos)>0.02) && (TMath::Abs(xyNeg)>0.05) ;
        }
 
        if(isLambda)
            FillHistogram("hLambdaMass",v0->MassLambda(),v0->Pt()) ;
        if(isLambdaBar)
            FillHistogram("hLambdaBarMass",v0->MassAntiLambda(),v0->Pt()) ;
        
        if(isLambda && TMath::Abs(v0->MassLambda()-1.115)>0.005)
            continue ;
        if(isLambdaBar && TMath::Abs(v0->MassAntiLambda()-1.115)>0.005)
            continue ;
        
        if(v0->Pt()<0.5)
            continue ;
        
        //So far combine Lambda and AntiLambda
	   
	if(isLambda || isLambdaBar){
          new((*fLambda)[inLambda++]) TLorentzVector(v0->Px(),v0->Py(),v0->Pz(),TMath::Sqrt(massLambda*massLambda+v0->P()*v0->P())) ;
            
        }
    }
}   
//_____________________________________________________________________________
Bool_t AliAnalysisPHOSResonances::AcceptTrack(const AliAODTrack *t) {
  if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  //if (t->GetKinkIndex(0)>0) return kFALSE;                                                                                                                

  Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < 70) return kFALSE;
  Int_t findable=t->GetTPCNclsF();
  if (findable <= 0) return kFALSE;
  if (nCrossedRowsTPC/findable < 0.8) return kFALSE;

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisPHOSResonances::SelectElectrons(){
 
    
  //Select electrons   
    
  Int_t inPHOSp=0;
  Int_t inPHOSm=0;
  if(fTracksElm)
      fTracksElm->Clear() ;
  else
      fTracksElm = new TClonesArray("TLorentzVector",10) ;
  if(fTracksElp)
      fTracksElp->Clear() ;
  else
      fTracksElp = new TClonesArray("TLorentzVector",10) ;
   
  Int_t multClust = fEvent->GetNumberOfCaloClusters();
 
  for (Int_t i=0; i<multClust; i++) {
    AliAODCaloCluster *clu = fEvent->GetCaloCluster(i);
    if (clu->GetType() !=AliVCluster::kPHOSNeutral ) continue;
    if (clu->E()<0.1) continue;
    if(clu->GetNCells()<3) continue ;
    if(clu->GetM02()<0.2)   continue ;    
//     if(clu->GetMCEnergyFraction()>0.98) //Ecross fCuts, should be filled with Tender
//      continue ;    
    if(clu->Chi2()>2.5*2.5)continue ;
    if(clu->GetEmcCpvDistance()>2.)continue ;
//     if((clu->GetTOF()>25.e-9) || (clu->GetTOF() <-25.e-9) )
//       continue ;
   
    Int_t ntrM = clu->GetNTracksMatched();
    AliAODTrack * track =0x0 ;
    if(ntrM>0)
      track = (AliAODTrack*)clu->GetTrackMatched(0);
    
    if(!track)
      continue ;  
    
    Double_t Ep=clu->E()/track->P();
    FillHistogram("hEp",Ep,clu->E()) ;
    if(Ep<0.9 || Ep>1.1)
      continue;  
           
    TLorentzVector pv1(Ep*track->Px(),Ep*track->Py(),Ep*track->Pz(),clu->E()) ; 
    Short_t charge =  track->Charge() ;
    if(charge>0){
      if(inPHOSp>=fTracksElp->GetSize()){
        fTracksElp->Expand(inPHOSp+5) ;
      }
      new((*fTracksElp)[inPHOSp++]) TLorentzVector(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    }
    else{
      if(inPHOSm>=fTracksElm->GetSize()){
        fTracksElm->Expand(inPHOSm+5) ;
      }
      new((*fTracksElm)[inPHOSm++]) TLorentzVector(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    }
  }
    
    
    
}
//_____________________________________________________________________________
void AliAnalysisPHOSResonances::SelectGamma(){
 
    
  //Select electrons   
    
  Int_t inPHOS=0, iPi0Merged=0;
  if(fGamma)
      fGamma->Clear() ;
  else
      fGamma = new TClonesArray("AliCaloPhoton",100) ;
   
  const AliAODVertex *esdVertex5 = fEvent->GetPrimaryVertex();

  Double_t vtx5[3] ={esdVertex5->GetX(),esdVertex5->GetY(),esdVertex5->GetZ()};

  Int_t multClust = fEvent->GetNumberOfCaloClusters();
 
  for (Int_t i=0; i<multClust; i++) {
    AliAODCaloCluster *clu = fEvent->GetCaloCluster(i);
    if (clu->GetType() !=AliVCluster::kPHOSNeutral ) continue;
    if (clu->E()<0.04) continue;
//     if(clu->GetNCells()<3) continue ;
    if(clu->GetM02()<0.2)   continue ;    
//     if(clu->GetMCEnergyFraction()>0.98) //Ecross fCuts, should be filled with Tender
//      continue ;    
//     if(clu->Chi2()>2.5*2.5)continue ;
//     if(clu->GetEmcCpvDistance()<2.5)continue ;
//     if((clu->GetTOF()>25.e-9) || (clu->GetTOF() <-25.e-9) )
//       continue ;
              
    FillHistogram("hClusterTOFvsE",clu->GetTOF(),clu->E()) ;
    
    TLorentzVector pv1 ;
    clu->GetMomentum(pv1 ,vtx5);
    if(inPHOS>=fGamma->GetSize()){
      fGamma->Expand(inPHOS+50) ;
    }
    AliCaloPhoton * p = new((*fGamma)[inPHOS++]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    p->SetDispBit(clu->Chi2()<2.5*2.5) ;
    p->SetTOFBit((clu->GetTOF()>-50.e-9) && (clu->GetTOF() <50.e-9) ) ;
    p->SetCPVBit(clu->GetEmcCpvDistance()>2.5) ;   
    p->SetCPV2Bit(clu->GetEmcCpvDistance()>1.) ;  
    p->SetPrimary(0) ; //no matched partner yet
    
    Double_t nSigmaPi = PionDispCut(clu->GetM02(), clu->GetM20(), clu->E()); 
    if(nSigmaPi>-1.5 && nSigmaPi<3){ //this is pi0
      AliCaloPhoton * pM = new((*fPi0Merged)[iPi0Merged++]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),TMath::Sqrt(TMath::Power(pv1.E(),2)+0.135*0.135)) ;
      pM->SetTOFBit((clu->GetTOF()>-50.e-9) && (clu->GetTOF() <50.e-9) ) ;
      pM->SetDispBit(TMath::Abs(nSigmaPi)<1.) ;
      pM->SetCPVBit(clu->GetEmcCpvDistance()>2.5) ;   
      pM->SetCPV2Bit(clu->GetEmcCpvDistance()>1.) ;  
      pM->SetPrimary(0) ; //no matched partner yet
      
    }
    
    
  }
    
    
    
}


//_____________________________________________________________________________
Double_t AliAnalysisPHOSResonances::Pi0Mass(Double_t /*pt*/){
  return 0.137 ;  
}
//_____________________________________________________________________________
Double_t AliAnalysisPHOSResonances::Pi0Width(Double_t /*pt*/){
  return 0.012;   //2sigma
}
//_____________________________________________________________________________
Double_t AliAnalysisPHOSResonances::EtaMass(Double_t /*pt*/){
  return 0.555 ;  
}
//_____________________________________________________________________________
Double_t AliAnalysisPHOSResonances::EtaWidth(Double_t /*pt*/){
  return 0.030;   //2sigma
}
//______________________________________________________
Double_t AliAnalysisPHOSResonances::PionDispCut(Double_t m02, Double_t m20, Double_t E){
   //Returns ditance to pi0 peak center in sigmas

   //No Disp cut for soft energies 
   if(E<25.)
     return 999;  
    
   //Parameterization using single pi0 simulation 
   Double_t longMpi = 1.857398e+00 + 1.208331e+01*TMath::Exp(-4.977723e-02*E) ; 
   Double_t longSpi = 3.820707e-01 + 1.000542e+00*TMath::Exp(-3.877147e-02*E) ; 
   Double_t shortMpi = 1.152118e+00 - 4.076138e-01*TMath::Exp(-2.372902e-02*E) ; 
   Double_t shortSpi = 1.517538e-01 + 9.382205e+00*TMath::Exp(-1.563037e-01*E) ; 
   Double_t powerNpi = 2.055773e+00 + 9.616408e+03*TMath::Exp(-2.664167e-01*E) ; 
   
   Double_t dx = (m02-longMpi)/longSpi ;
   Double_t dy = (m20-shortMpi)/shortSpi ;
   
   //we have non-gaussian power, so re-calculate in Gaussian sigmas 
   return TMath::Sign(TMath::Sqrt(TMath::Power(TMath::Abs(dx),powerNpi)+TMath::Power(TMath::Abs(dy),powerNpi)),dx) ;
    
    
}

