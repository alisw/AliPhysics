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
 
// Analysis task for gamma-hadron fluctuation study
// Authors: Dmitri Peresunko and Evgenia Nekgasova

#include "TH2.h"
#include "TH3.h"
#include "THashList.h"
#include "TGeoGlobalMagField.h"
#include "TList.h"
#include "TClonesArray.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisPHOSFluctuations.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliAODInputHandler.h"
#include "AliPIDResponse.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCParticle.h"
#include "AliPPVsMultUtils.h"
#include "AliMultSelection.h"

#include "AliCaloPhoton.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include "AliAODTrack.h"

#include <vector>

ClassImp(AliAnalysisPHOSFluctuations)

//________________________________________________________________________
AliAnalysisPHOSFluctuations::AliAnalysisPHOSFluctuations(const char *name) 
: AliAnalysisTaskSE(name),
  fPhPtMin(0.3), 
  fPhPtMax(1.),
  fPi0PtMin(0.3), 
  fPi0PtMax(1.),
  fChPtMin(0.3), 
  fChPtMax(1.),
  fChEtaCutMax(0.8), 
  fPhEtaCutMax(0.8), 
  fPi0EtaCutMax(0.8), 
  fChPhiMin(0.), 
  fChPhiMax(6.2831853 ),
  fPhPhiMin(0.), 
  fPhPhiMax(6.2831853 ),
  fPi0PhiMin(0.), 
  fPi0PhiMax(6.2831853 ),
  fRunNumber(0),                           
  fCentrality(0),                         
  fCentBin(0),
  fRunType(kpp),
  fPPUtils(nullptr),          
  fOutputContainer(nullptr),         
  fCurrentMixedList(nullptr),          
  fStack(nullptr),           
  fPIDResponse (nullptr),   
  fUtils(nullptr),           
  fArrGamma (nullptr),      
  fRecPipm(0),
  fRecPipmTrue(0),
  fhSelEvents(nullptr),     
  fhTrackPt(nullptr),      
  fhPionPt(nullptr),       
  fhKaonPt(nullptr),       
  fhProtonPt(nullptr),     
  fhUndefPt(nullptr),      
  fhPhotonPt(nullptr),     
  fhMCPhotonPt(nullptr),   
  fhMCPrimPi0N(nullptr),        
  fhMCPrimPi01N (nullptr),       
  fhMCPrimPi0NoresN(nullptr),   
  fhMCPrimPi0Nores1N(nullptr),  
  fhMCPrimGammaN(nullptr),      
  fhMCPrimGamma1N(nullptr),      
  fhMCPrimGammaPi0N(nullptr),   
  fhMCPrimGammaPi01N(nullptr),  
  fhMCPrimGammaPi0SingleN(nullptr),        
  fhMCPrimGammaPi0Single1N(nullptr),       
  fhMCPrimGammaAllSingleN(nullptr),        
  fhMCPrimGammaAllSingle1N(nullptr),       
  fhMCPrimGammaPi0SingleNoresN(nullptr),   
  fhMCPrimGammaPi0SingleNores1N(nullptr),  
  fhMCPrimPipmN(nullptr),                      
  fhMCPrimPipm1N (nullptr),                   
  fhMCPrimPipmNoresNa(nullptr),              
  fhMCPrimPipmNores1Na (nullptr),             
  fhMCPrimPipmNoresNb(nullptr),              
  fhMCPrimPipmNores1Nb (nullptr),             
  fhMCPrimPipmPi0(nullptr),                  
  fhMCPrimPipmPi0Nores(nullptr),             
  fhMCPrimPipmGamma(nullptr),                
  fhMCPrimPipmGammaPi0 (nullptr),             
  fhMCPrimPipmGammaPi0Single (nullptr),       
  fhMCPrimPipmGammaAllSingle (nullptr),       
  fhMCPrimPipmGammaPi0SingleNores(nullptr),  
  fhPipmN(nullptr),                           
  fhPipm1N (nullptr),                          
  fhPipmTrueN(nullptr),                       
  fhPipmTrue1N(nullptr),                      
  fEgammaEpi0(nullptr)                        
{
    
  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());

  for(Int_t i=0;i<20;i++)
    for(Int_t j=0;j<NCENT;j++)
      fPHOSEvents[i][j]=nullptr ;    //Container for PHOS photons
  for(Int_t i=0;i<NPID;i++){
    fRecPhot[i] =0;
    fRecPhotTrue[i] =0;
    fRecPhotTruePi0[i] =0;
    fRecPhotTruePi0Single[i] =0;
    fhGammaN[i]=nullptr ;                  
    fhGamma1N[i]=nullptr ;                 
    fhGammaTrueN[i]=nullptr ;              
    fhGammaTrue1N[i]=nullptr ;             
    fhGammaTruePi0N[i]=nullptr ;           
    fhGammaTruePi01N[i]=nullptr ;          
    fhGammaTruePi0SingleN[i]=nullptr ;     
    fhGammaTruePi0Single1N[i]=nullptr ;    

    fhGammaPipm[i]=nullptr ;               
    fhGammaPipmTrue[i]=nullptr ;           
    fhGammaPipmTruePi0[i]=nullptr ;        
    fhGammaPipmTruePi0Single[i]=nullptr ;  
    fhReal[i]=nullptr ;                     
    fhRealTrue[i]=nullptr ;                  
    fhRealCommon[i]=nullptr ;              
    fhMixed[i]=nullptr ;                     
  }
}
//________________________________________________________________________
void AliAnalysisPHOSFluctuations::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  // AOD histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  const int nCenBin=20;
  //General QA
  fhSelEvents = new TH1F("hSelEvents","Selected events",20,0.,20.) ;
  fOutputContainer->Add(fhSelEvents);
  fOutputContainer->Add(new TH1F("hZvertex","Z vertex",200,-50.,+50.));
  fOutputContainer->Add(new TH1F("hNvertexTracks","N of primary tracks from the primary vertex",150,0.,150.));
  fOutputContainer->Add(new TH1F("hClusterMult","hClusterMult",200,0.,200.));
  fOutputContainer->Add(new TH1F("hCentralityV0M","hCentralityV0M",100,0.,100.));
  fOutputContainer->Add(new TH2F("hTrackCentralityV0M","hTrackCentralityV0M",100,0.,100.,100,0.,2000.));


  fhTrackPt = new TH1F("hTrackPt","hTrackPt",100,0.,10.) ;
  fOutputContainer->Add(fhTrackPt);
  fhPionPt = new TH1F("hPionPt","hTrackPt",100,0.,10.) ;
  fOutputContainer->Add(fhPionPt);  

  //MC histos
  fhMCPrimPi0N = new TH1F("hMCPrimPi0N",Form("MC #pi^{0} <N> %5.3f<pT<%5.3f",fPi0PtMin,fPi0PtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPi0N);
  fhMCPrimPi01N = new TH1F("hMCPrimPi01N",Form("MC #pi^{0} <N(N-1)> %5.3f<pT<%5.3f",fPi0PtMin,fPi0PtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPi01N);
  fhMCPrimPi0NoresN = new TH1F("hMCPrimPi0NoresN",Form("MC #pi^{0} no res <N> %5.3f<pT<%5.3f",fPi0PtMin,fPi0PtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPi0NoresN);
  fhMCPrimPi0Nores1N = new TH1F("hMCPrimPi0Nores1N",Form("MC #pi^{0} no res <N(N-1)> %5.3f<pT<%5.3f",fPi0PtMin,fPi0PtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPi0Nores1N);

  fhMCPrimGammaN = new TH1F("hMCPrimGammaN",Form("MC #gamma <N> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimGammaN);
  fhMCPrimGamma1N = new TH1F("hMCPrimGamma1N",Form("MC #gamma <N(N-1)> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimGamma1N);
  fhMCPrimGammaPi0N = new TH1F("hMCPrimGammaPi0N",Form("MC #gamma #pi^{0} <N> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimGammaPi0N);
  fhMCPrimGammaPi01N = new TH1F("hMCPrimGammaPi01N",Form("MC #gamma #pi^{0} <N(N-1)> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimGammaPi01N);
  fhMCPrimGammaPi0SingleN = new TH1F("hMCPrimGammaPi0SingleN",Form("MC #gamma_{s} <N> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimGammaPi0SingleN);
  fhMCPrimGammaPi0Single1N = new TH1F("hMCPrimGammaPi0Single1N",Form("MC #gamma_{s} <N(N-1)> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimGammaPi0Single1N);
  fhMCPrimGammaAllSingleN = new TH1F("hMCPrimGammaAllSingleN",Form("MC #gamma_{s} <N> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimGammaAllSingleN);
  fhMCPrimGammaAllSingle1N = new TH1F("hMCPrimGammaAllSingle1N",Form("MC #gamma_{s} <N(N-1)> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimGammaAllSingle1N);

  fhMCPrimGammaPi0SingleNoresN = new TH1F("hMCPrimGammaPi0SingleNoresN",Form("MC #gamma_{s} <N> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimGammaPi0SingleNoresN);
  fhMCPrimGammaPi0SingleNores1N = new TH1F("hMCPrimGammaPi0SingleNores1N",Form("MC #gamma_{s} <N(N-1)> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimGammaPi0SingleNores1N);


  fhMCPrimPipmN = new TH1F("hMCPrimPipmN",Form("MC #pi^{#pm} <N> %5.3f<pT<%5.3f",fChPtMin,fChPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipmN);
  fhMCPrimPipm1N = new TH1F("hMCPrimPipm1N",Form("MC #pi^{0} <N(N-1)> %5.3f<pT<%5.3f",fChPtMin,fChPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipm1N);

  fhMCPrimPipmPi0 = new TH1F("hMCPrimPipmPi0",Form("MC <N_{#pi}N_{#pi}> %5.3f<pT<%5.3f",fChPtMin,fChPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipmPi0);
  fhMCPrimPipmGamma = new TH1F("hMCPrimPipmGamma",Form("MC <N_{#gamma}N_{#pi}> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipmGamma);
  fhMCPrimPipmGammaPi0 = new TH1F("hMCPrimPipmGammaPi0",Form("MC <N_{#gamma}N_{#pi}> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipmGammaPi0);
  fhMCPrimPipmGammaPi0Single = new TH1F("hMCPrimPipmGammaPi0Single",Form("MC <N_{#gamma}N_{#pi}> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipmGammaPi0Single);
  fhMCPrimPipmGammaAllSingle = new TH1F("hMCPrimPipmGammaAllSingle",Form("MC <N_{#gamma}N_{#pi}> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipmGammaAllSingle);
  fhMCPrimPipmGammaPi0SingleNores = new TH1F("hMCPrimPipmGammaPi0SingleNores",Form("MC <N_{#gamma}N_{#pi}> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipmGammaPi0SingleNores);


  fhMCPrimPipmNoresNa = new TH1F("hMCPrimPipmNoresNp",Form("MC #pi^{#pm} no res <N> %5.3f<pT<%5.3f",fChPtMin,fChPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipmNoresNa);
  fhMCPrimPipmNores1Na = new TH1F("hMCPrimPipmNores1Np",Form("MC #pi^{#pm} no res <N(N-1)> %5.3f<pT<%5.3f",fChPtMin,fChPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipmNores1Na);
  fhMCPrimPipmNoresNb = new TH1F("hMCPrimPipmNoresNg",Form("MC #pi^{#pm} no res <N> %5.3f<pT<%5.3f",fChPtMin,fChPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipmNoresNb);
  fhMCPrimPipmNores1Nb = new TH1F("hMCPrimPipmNores1Ng",Form("MC #pi^{#pm} no res <N(N-1)> %5.3f<pT<%5.3f",fChPtMin,fChPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipmNores1Nb);
  fhMCPrimPipmPi0Nores = new TH1F("hMCPrimPipmPi0Nores",Form("MC #pi^{#pm} no res <N_{#pi}N_{#pi}> %5.3f<pT<%5.3f",fChPtMin,fChPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhMCPrimPipmPi0Nores);

  //reco histos
  fhPipmN = new TH1F("hPipmN",Form("#pi^{#pm} <N> %5.3f<pT<%5.3f",fChPtMin,fChPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhPipmN);
  fhPipm1N = new TH1F("hPipm1N",Form("#pi^{#pm} <N(N-1)> %5.3f<pT<%5.3f",fChPtMin,fChPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhPipm1N);
  fhPipmTrueN = new TH1F("hPipmTrueN",Form("true #pi^{#pm} <N> %5.3f<pT<%5.3f",fChPtMin,fChPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhPipmTrueN);
  fhPipmTrue1N = new TH1F("hPipmTrue1N",Form("true #pi^{#pm} <N(N-1)> %5.3f<pT<%5.3f",fChPtMin,fChPtMax),nCenBin,0.,100.) ;
  fOutputContainer->Add(fhPipmTrue1N);

  for(int iPID=0; iPID<NPID; iPID++){
    fhGammaN[iPID] = new TH1F(Form("hGammaNpid%d",iPID),Form("#gamma <N> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
    fOutputContainer->Add(fhGammaN[iPID]);
    fhGamma1N[iPID] = new TH1F(Form("hGamma1Npid%d",iPID),Form("#gamma <N(N-1)> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
    fOutputContainer->Add(fhGamma1N[iPID]);
    fhGammaTrueN[iPID] = new TH1F(Form("hGammaTrueNpid%d",iPID),Form("#gamma true <N> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
    fOutputContainer->Add(fhGammaTrueN[iPID]);
    fhGammaTrue1N[iPID] = new TH1F(Form("hGammaTrue1Npid%d",iPID),Form("#gamma true <N(N-1)> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
    fOutputContainer->Add(fhGammaTrue1N[iPID]);

    fhGammaTruePi0N[iPID] = new TH1F(Form("hGammaTruePi0Npid%d",iPID),Form("#gamma true pi0 <N> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
    fOutputContainer->Add(fhGammaTruePi0N[iPID]);
    fhGammaTruePi01N[iPID] = new TH1F(Form("hGammaTruePi01Npid%d",iPID),Form("#gamma true pi0 <N(N-1)> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
    fOutputContainer->Add(fhGammaTruePi01N[iPID]);
    fhGammaTruePi0SingleN[iPID] = new TH1F(Form("hGammaTruePi0SingleNpid%d",iPID),Form("#gamma_{s} true pi0 <N(N-1)> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
    fOutputContainer->Add(fhGammaTruePi0SingleN[iPID]);
    fhGammaTruePi0Single1N[iPID] = new TH1F(Form("hGammaTruePi0Single1Npid%d",iPID),Form("#gamma_{s} true pi0 <N(N-1)> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
    fOutputContainer->Add(fhGammaTruePi0Single1N[iPID]);

    fhGammaPipm[iPID] = new TH1F(Form("hGammaPipmpid%d",iPID),Form("<N_{#gamma}N_{#pi}> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
    fOutputContainer->Add(fhGammaPipm[iPID]);
    fhGammaPipmTrue[iPID] = new TH1F(Form("hGammaPipmTruepid%d",iPID),Form("true <N_{#gamma}N_{#pi}> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
    fOutputContainer->Add(fhGammaPipmTrue[iPID]);
    fhGammaPipmTruePi0[iPID] = new TH1F(Form("hGammaPipmTruePi0pid%d",iPID),Form("true #pi^{0} <N_{#gamma}N_{#pi}> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
    fOutputContainer->Add(fhGammaPipmTruePi0[iPID]);
    fhGammaPipmTruePi0Single[iPID] = new TH1F(Form("hGammaPipmTruePi0Singlepid%d",iPID),Form("true single <N_{#gamma}N_{#pi}> %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),nCenBin,0.,100.) ;
    fOutputContainer->Add(fhGammaPipmTruePi0Single[iPID]);
  }

  //Inv mass
  for(int iPID=0; iPID<NPID; iPID++){
    fhReal[iPID] = new TH2F(Form("hRealpid%d",iPID),Form("real %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),500,0.,1.,nCenBin,0.,100.) ;
    fOutputContainer->Add(fhReal[iPID]);
    fhRealTrue[iPID] = new TH2F(Form("hRealTruepid%d",iPID),Form("true gamma real %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),500,0.,1.,nCenBin,0.,100.) ;
    fOutputContainer->Add(fhRealTrue[iPID]);
    fhRealCommon[iPID] = new TH2F(Form("hRealCommonpid%d",iPID),Form("common parent %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),500,0.,1.,nCenBin,0.,100.) ;
    fOutputContainer->Add(fhRealCommon[iPID]);
    fhMixed[iPID] = new TH2F(Form("hMixedpid%d",iPID),Form("mixed %5.3f<pT<%5.3f",fPhPtMin,fPhPtMax),500,0.,1.,nCenBin,0.,100.) ;
    fOutputContainer->Add(fhMixed[iPID]);
  }  
  fEgammaEpi0 = new TH2F("EgammavsEpi0","Pt gamma vs Pt pi0",100,0.,2.,100,0.,2.) ;
  fOutputContainer->Add(fEgammaEpi0) ;

  fPPUtils = new AliPPVsMultUtils();

  PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisPHOSFluctuations::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze AOD
  AliAODEvent* event = (AliAODEvent*)InputEvent();
  if(!event){
    AliDebug(1,"No event") ;
    PostData(1, fOutputContainer);
    return;
  }
  fhSelEvents->Fill(1) ;
  fRunNumber=event->GetRunNumber() ;

  if(!fUtils) 
      fUtils = new AliAnalysisUtils();
  
  Bool_t isMB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7)  ; 
           
  if(!isMB ){
    PostData(1, fOutputContainer);
    return ;
  }
    
  fhSelEvents->Fill(2) ;
  
  // Checks if we have a primary vertex
  // Get primary vertices form AOD

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
  fhSelEvents->Fill(3) ;

  if(!fUtils->IsVertexSelected2013pA(event)){
    PostData(1, fOutputContainer);
    return ;
  }
  
  fhSelEvents->Fill(4) ;
  
  if(fUtils->IsPileUpEvent(event)){
    PostData(1, fOutputContainer);
    return ;
  }
  fhSelEvents->Fill(5) ;

  AliMultSelection *multSelection = (AliMultSelection*) event -> FindListObject("MultSelection");

  fCentrality =multSelection->GetMultiplicityPercentile("V0M");
  FillHistogram("hCentralityV0M",fCentrality) ;

  // if ((fRunType==1)||(fRunType==2)){

  //   AliCentrality *centrality = event->GetCentrality();
  //   fCentralityCL1=centrality->GetCentralityPercentile("CL1");
  //   fCentralityV0A=centrality->GetCentralityPercentile("V0A");
  //   fCentralityZNA=centrality->GetCentralityPercentile("ZNA");

  //   FillHistogram("hCentralityV0A",fCentralityV0A) ;
  //   FillHistogram("hCentralityCL1",fCentralityCL1) ;
  //   FillHistogram("hCentralityZNA",fCentralityZNA) ;

  // }

  fhSelEvents->Fill(6) ;
  
  if(!fPIDResponse){
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
  }
    
  fStack = (TClonesArray*)event->FindListObject(AliAODMCParticle::StdBranchName());
  
  Int_t zvtx = TMath::Min(9,Int_t((vtx5[2]+10.)/2.)) ; 
  fCentBin=fCentrality/20.;
  if(fCentBin>=NCENT)fCentBin=NCENT-1;

  if(!fPHOSEvents[zvtx][fCentBin]) 
    fPHOSEvents[zvtx][fCentBin]=new TList() ;
  fCurrentMixedList = fPHOSEvents[zvtx][fCentBin] ;

  ProcessMC() ;

  SelectPhotons() ;
  
  SelectTracks() ;

  //Fill histograms
  fhPipmN->Fill(fCentrality,fRecPipm) ;
  fhPipm1N->Fill(fCentrality,fRecPipm*(fRecPipm-1)) ;
  fhPipmTrueN->Fill(fCentrality,fRecPipmTrue) ;
  fhPipmTrue1N->Fill(fCentrality,fRecPipmTrue*(fRecPipmTrue-1)) ;
  
  
  for(Int_t iPID=0; iPID<NPID; iPID++){
    fhGammaN[iPID]->Fill(fCentrality,fRecPhot[iPID]) ;
    fhGamma1N[iPID]->Fill(fCentrality,(fRecPhot[iPID]*(fRecPhot[iPID]-1)) ) ;
    fhGammaTrueN[iPID]->Fill(fCentrality,fRecPhotTrue[iPID]) ;
    fhGammaTrue1N[iPID]->Fill(fCentrality,(fRecPhotTrue[iPID]*(fRecPhotTrue[iPID]-1)) ) ;
    fhGammaTruePi0N[iPID]->Fill(fCentrality,fRecPhotTruePi0[iPID]) ;
    fhGammaTruePi01N[iPID]->Fill(fCentrality,(fRecPhotTruePi0[iPID]*(fRecPhotTruePi0[iPID]-1)) ) ;
    fhGammaTruePi0SingleN[iPID]->Fill(fCentrality,fRecPhotTruePi0Single[iPID]) ;
    fhGammaTruePi0Single1N[iPID]->Fill(fCentrality,(fRecPhotTruePi0Single[iPID]*(fRecPhotTruePi0Single[iPID]-1)) ) ;

    fhGammaPipm[iPID]->Fill(fCentrality,fRecPhot[iPID]*fRecPipm);
    fhGammaPipmTrue[iPID]->Fill(fCentrality,fRecPhotTrue[iPID]*fRecPipmTrue);
    fhGammaPipmTruePi0[iPID]->Fill(fCentrality,fRecPhotTruePi0[iPID]*fRecPipmTrue);
    fhGammaPipmTruePi0Single[iPID]->Fill(fCentrality,fRecPhotTruePi0Single[iPID]*fRecPipmTrue);
  }

  // Post output data.
  PostData(1, fOutputContainer);
}  
//________________________________________________________________________
void AliAnalysisPHOSFluctuations::SelectTracks() 
{  
  //Select charged pions
  //Calculate charged multiplicity
  fRecPipm=0 ;
  fRecPipmTrue=0 ;
  
  int nTracks= fInputEvent->GetNumberOfTracks() ;
  
  FillHistogram("hTrackCentralityV0M",fCentrality,nTracks) ;
  for (Int_t i=0;i<nTracks;++i) {
    AliAODTrack *track = (AliAODTrack*)fInputEvent->GetTrack(i) ;
    if(!track->IsHybridGlobalConstrainedGlobal())
      continue ;
    double pt = track->Pt();
    if(pt<fChPtMin || pt >fChPtMax)
      continue ;

    if(track->GetType()!= AliAODTrack::kPrimary ){
     continue;
    }
    double eta = TMath::Abs(track->Eta()) ;
    if(eta>fChEtaCutMax){
      continue ;  
    }
    double phi=track->Phi() ;
    while(phi<0.)phi+=TMath::TwoPi() ;
    while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
    if(phi<fChPhiMin || phi >fChPhiMax)
      continue ;
    fhTrackPt->Fill(pt);
    bool pidPion=kFALSE, pidKaon=kFALSE, pidProton=kFALSE, pidUndef=kFALSE;
    if(fPIDResponse) {
	    double nsigmaProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
	    double nsigmaKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)); 
	    double nsigmaPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)); 
        
    	// guess the particle based on the smaller nsigma
	    if((nsigmaKaon  <nsigmaPion) && (nsigmaKaon  <nsigmaProton) && (nsigmaKaon  <3)) pidKaon   = kTRUE;
	    if((nsigmaPion  <nsigmaKaon) && (nsigmaPion  <nsigmaProton) && (nsigmaPion  <3)) pidPion   = kTRUE;
	    if((nsigmaProton<nsigmaKaon) && (nsigmaProton<nsigmaPion  ) && (nsigmaProton<3)) pidProton = kTRUE;
	    if (!pidPion && !pidKaon && !pidProton) pidUndef = kTRUE;
    }

	  if (pidPion) {
      fhPionPt->Fill(pt);
      fRecPipm++ ;
      if(fStack){
        int iprim=track->GetLabel() ;
        if(iprim>-1){
          AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(iprim) ;
          if(TMath::Abs(prim->GetPdgCode())==211 ){
            fRecPipmTrue++ ;
          }
        }
      }
    }
	}
}
//________________________________________________________________________
void AliAnalysisPHOSFluctuations:: ProcessMC(){
    
  //Scan MC if available and construct dynamic variable from primary particles  
    
  if(!fStack)
    return ;
  
  int nMCPrimPi0 = 0 ;       //number of primary pi0 
  int nMCPrimPi0Nores = 0;   //number of primary pi0 
  int nMCPrimGamma = 0;      //number of primary gamma 
  int nMCPrimGammaPi0 = 0 ;  //number of primary gamma from pi0
  int nMCPrimGammaPi0Single = 0;      //number of primary gamma from pi0
  int nMCPrimGammaAllSingle = 0;      //number of primary gamma single from pi0 and other hadr
  int nMCPrimGammaPi0SingleNores = 0; //number of primary gamma from pi0
  int nMCPrimPipm = 0;         //number of primary pi+- 
  int nMCPrimPipmNores1 = 0 ;   //number of primary pi+- no common parent with pi0
  int nMCPrimPipmNores2 = 0 ;   //number of primary pi+- no common parent with gamma

  fChPrimaryList1.clear() ;
  fChPrimaryList2.clear() ;
  fPhPrimaryList.clear() ;
  fPi0PrimaryList.clear() ;
  fpi0list.clear() ;

  int nPrim = fStack->GetEntriesFast() ;

  for(Int_t i=0;i<nPrim;i++){
    AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(i) ;
    double r2=prim->Xv()*prim->Xv()+prim->Yv()*prim->Yv() ;
    if(r2>0.25)
      continue ;
    Int_t pdg=prim->GetPdgCode();

    if(pdg==111){ //pi0
      double pt = prim->Pt() ;
      if(pt<fPi0PtMin || pt>fPi0PtMax)
        continue ;
      if(TMath::Abs(prim->Eta())>fPi0EtaCutMax)
        continue; 
      double phi=prim->Phi() ;
      while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
      while(phi<0)phi+=TMath::TwoPi() ;
      if(phi<fPi0PhiMin || phi>fPi0PhiMax)
        continue ;

      nMCPrimPi0++ ;
      fPi0PrimaryList.push_back(i) ;
    }
    if(pdg==22){ //gamma
      double pt = prim->Pt() ;
      if(pt<fPhPtMin || pt>fPhPtMax)
        continue ;
      if(TMath::Abs(prim->Eta())>fPhEtaCutMax)
        continue; 
      double phi=prim->Phi() ;
      while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
      while(phi<0)phi+=TMath::TwoPi() ;
      if(phi<fPhPhiMin || phi>fPhPhiMax)
        continue ;

      nMCPrimGamma++ ;
      int parent = prim->GetMother() ;
      if(parent>-1){
        int pdgParent = ((AliAODMCParticle*)fStack->At(parent))->GetPdgCode() ;
        if(pdgParent==111){
          nMCPrimGammaPi0++;
          fEgammaEpi0->Fill(pt,((AliAODMCParticle*)fStack->At(parent))->Pt());

          bool found = false ;
          for(int iparpi0:fpi0list){
            if(iparpi0==parent){
              found=true;
              break;
            }
          }
          if(!found){
            nMCPrimGammaPi0Single++;
            nMCPrimGammaAllSingle++;
            fpi0list.push_back(parent) ;
            fPhPrimaryList.push_back(i) ;
          }
        }
        else{
          nMCPrimGammaAllSingle++;
        }
      }
      else{
        nMCPrimGammaAllSingle++;       
      }
    }
    if(pdg==211 || pdg ==-211){ //pipm
      double pt = prim->Pt() ;
      if(pt<fChPtMin || pt>fChPtMax)
        continue ;
      if(TMath::Abs(prim->Eta())>fChEtaCutMax)
        continue; 
      double phi=prim->Phi() ;
      while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
      while(phi<0)phi+=TMath::TwoPi() ;
      if(phi<fChPhiMin || phi>fChPhiMax)
        continue ;
      nMCPrimPipm++ ;
      fChPrimaryList1.push_back(i) ;
      fChPrimaryList2.push_back(i) ;
    }
  }

  //Remove resonance decays
  for(int i=0; i<nMCPrimPi0; i++){
    int prim1=fPi0PrimaryList[i];
    for(int j=0; j<nMCPrimPipm; j++){
      int prim2=fChPrimaryList1[j];
      int pdg = CommonAnsestor(prim1,prim2) ;
      if(pdg!=0 && pdg!=2212 && pdg!=2112){ //EPOS assignment
        // printf("Common pi0+pi = %d\n",pdg) ;
        fPi0PrimaryList[i]=-1;
        fChPrimaryList1[j]=-1;
      }
    }
  }
  for(int i=0; i<nMCPrimGammaPi0Single; i++){
    int prim1=fPhPrimaryList[i];
    for(int j=0; j<nMCPrimPipm; j++){
      int prim2=fChPrimaryList2[j];
      int pdg = CommonAnsestor(prim1,prim2) ;
      if(pdg!=0 && pdg!=2212 && pdg!=2112){
        // printf("Common gamma+pi = %d\n",pdg) ;
        fPhPrimaryList[i]=-1;
        fChPrimaryList2[j]=-1;
      }
    }
  }
  for(int i=0; i<nMCPrimPi0; i++){
    if(fPi0PrimaryList[i]>=0){
      nMCPrimPi0Nores++;
    }
  }
  for(int i=0; i<nMCPrimGammaPi0Single; i++){
    if(fPhPrimaryList[i]>=0){
      nMCPrimGammaPi0SingleNores++;
    }
  }
  for(int i=0; i<nMCPrimPipm; i++){
    if(fChPrimaryList1[i]>=0){
      nMCPrimPipmNores1++;
    }
    if(fChPrimaryList2[i]>=0){
      nMCPrimPipmNores2++;
    }
  }

  //Fill histograms
  fhMCPrimPi0N->Fill(fCentrality,nMCPrimPi0) ;
  fhMCPrimPi01N->Fill(fCentrality,nMCPrimPi0*(nMCPrimPi0-1)) ;
  fhMCPrimPi0NoresN->Fill(fCentrality,nMCPrimPi0Nores) ;
  fhMCPrimPi0Nores1N->Fill(fCentrality,nMCPrimPi0Nores*(nMCPrimPi0Nores-1)) ;

  fhMCPrimGammaN->Fill(fCentrality,nMCPrimGamma) ;
  fhMCPrimGamma1N->Fill(fCentrality,nMCPrimGamma*(nMCPrimGamma-1)) ;
  fhMCPrimGammaPi0N->Fill(fCentrality,nMCPrimGammaPi0) ;
  fhMCPrimGammaPi01N->Fill(fCentrality,nMCPrimGammaPi0*(nMCPrimGammaPi0-1)) ;
  fhMCPrimGammaPi0SingleN->Fill(fCentrality,nMCPrimGammaPi0Single) ;
  fhMCPrimGammaPi0Single1N->Fill(fCentrality,nMCPrimGammaPi0Single*(nMCPrimGammaPi0Single-1)) ;
  fhMCPrimGammaAllSingleN->Fill(fCentrality,nMCPrimGammaAllSingle) ;
  fhMCPrimGammaAllSingle1N->Fill(fCentrality,nMCPrimGammaAllSingle*(nMCPrimGammaAllSingle-1)) ;
  fhMCPrimGammaPi0SingleNoresN->Fill(fCentrality,nMCPrimGammaPi0SingleNores) ;
  fhMCPrimGammaPi0SingleNores1N->Fill(fCentrality,nMCPrimGammaPi0SingleNores*(nMCPrimGammaPi0SingleNores-1)) ;

  fhMCPrimPipmN->Fill(fCentrality,nMCPrimPipm) ;
  fhMCPrimPipm1N->Fill(fCentrality,nMCPrimPipm*(nMCPrimPipm-1)) ;
  fhMCPrimPipmPi0->Fill(fCentrality,nMCPrimPipm*nMCPrimPi0) ;
  fhMCPrimPipmNoresNa->Fill(fCentrality,nMCPrimPipmNores1) ;
  fhMCPrimPipmNores1Na->Fill(fCentrality,nMCPrimPipmNores1*(nMCPrimPipmNores1-1)) ;
  fhMCPrimPipmPi0Nores->Fill(fCentrality,nMCPrimPipmNores1*nMCPrimPi0Nores) ;

  fhMCPrimPipmGamma->Fill(fCentrality,nMCPrimPipm*nMCPrimGamma) ;
  fhMCPrimPipmGammaPi0->Fill(fCentrality,nMCPrimPipm*nMCPrimGammaPi0) ;
  fhMCPrimPipmGammaPi0Single->Fill(fCentrality,nMCPrimPipm*nMCPrimGammaPi0Single) ;
  fhMCPrimPipmGammaAllSingle->Fill(fCentrality,nMCPrimPipm*nMCPrimGammaAllSingle) ;
  fhMCPrimPipmNoresNb->Fill(fCentrality,nMCPrimPipmNores2) ;
  fhMCPrimPipmNores1Nb->Fill(fCentrality,nMCPrimPipmNores2*(nMCPrimPipmNores2-1)) ;
  fhMCPrimPipmGammaPi0SingleNores->Fill(fCentrality,nMCPrimPipmNores2*nMCPrimGammaPi0SingleNores) ;

}
//________________________________________________________________________
void AliAnalysisPHOSFluctuations::Terminate(Option_t *)
{
}

//_____________________________________________________________________________
void AliAnalysisPHOSFluctuations::FillHistogram(const char * key,Double_t x)const{
  //FillHistogram
  TH1 * hist = dynamic_cast<TH1*>(fOutputContainer->FindObject(key)) ;
  if(hist)
    hist->Fill(x) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisPHOSFluctuations::FillHistogram(const char * key,Double_t x,Double_t y)const{
  //FillHistogram
  TH1 * th1 = dynamic_cast<TH1*> (fOutputContainer->FindObject(key));
  if(th1)
    th1->Fill(x, y) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}

//_____________________________________________________________________________
void AliAnalysisPHOSFluctuations::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
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
  
  AliError(Form("can not find histogram (of instance TH3) <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisPHOSFluctuations::SelectPhotons(){

  if(!fArrGamma)
    fArrGamma = new TClonesArray("AliCaloPhoton",100);
  else
    fArrGamma->Clear();

  frecpi0.clear() ;
  const Double_t kTOFMaxCut= 30.e-9 ;  
  const Double_t kTOFMinCut=-30.e-9 ;  
  const Double_t rcut = 1. ; // cut to define primary photons

  Int_t nGamma=0;
  for(Int_t iPID=0; iPID<NPID; iPID++){  
    fRecPhot[iPID] =0 ;
    fRecPhotTrue[iPID] =0 ;
    fRecPhotTruePi0[iPID] =0 ;
    fRecPhotTruePi0Single[iPID] =0 ;
  }
    
  Int_t multClust = fInputEvent->GetNumberOfCaloClusters();
  FillHistogram("hClusterMult",multClust);

  Double_t vtx5[3];
  vtx5[0] = fInputEvent->GetPrimaryVertex()->GetX();
  vtx5[1] = fInputEvent->GetPrimaryVertex()->GetY();
  vtx5[2] = fInputEvent->GetPrimaryVertex()->GetZ();
  TVector3 vertex(vtx5);
  
  for (Int_t i=0; i<multClust; i++) {
    AliAODCaloCluster *clu = (AliAODCaloCluster*)fInputEvent->GetCaloCluster(i);
    if ( clu->GetType() != AliVCluster::kPHOSNeutral ) continue;    
     
//    if(clu->GetNCells()<3) continue ;
//    if(clu->GetM02()<0.2)continue ;
   if(clu->GetTOF() < kTOFMinCut || clu->GetTOF() > kTOFMaxCut)
     continue ;          
  
    TLorentzVector p;
    clu->GetMomentum(p,vtx5);
    double pt=p.Pt();
    if(pt<fPhPtMin || pt>fPhPtMax)
        continue ;

    bool dispBit = (clu->GetDispersion()<2.5*2.5) ;
    bool cpvBit = clu->GetEmcCpvDistance()>2.5 ;  
    
    bool isPhoton = kFALSE ;
    bool isFromPi0 = kFALSE ;
    bool isSingle = true ;
    if(fStack){ //Find primary
      int primLabel=clu->GetLabelAt(0) ; //FindPrimary(clu,sure) ;
      //Look what particle left vertex
      if(primLabel>-1){
        AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(primLabel) ;
        int iparent=primLabel;
        AliAODMCParticle * parent = prim;
        double r2=(prim->Xv()-vtx5[0])*(prim->Xv()-vtx5[0])+(prim->Yv()-vtx5[1])*(prim->Yv()-vtx5[1]) ;
        while((r2 > rcut*rcut) && (iparent>-1)){
          iparent=parent->GetMother();
          parent=(AliAODMCParticle*)fStack->At(iparent);
          r2=(parent->Xv()-vtx5[0])*(parent->Xv()-vtx5[0])+(parent->Yv()-vtx5[1])*(parent->Yv()-vtx5[1]) ;
        }
        if(parent->GetPdgCode()==111){
          isFromPi0=kTRUE;
        }
        if((r2 < rcut*rcut) && (parent->GetPdgCode()==22)){
          isPhoton=true ;
          int iMother = parent->GetMother() ;
          //find parent pi0
          while(iMother>-1){
            AliAODMCParticle * mthr = (AliAODMCParticle*)fStack->At(iMother) ;
            if(mthr->GetPdgCode()==111){
              isFromPi0=kTRUE;
              for(int ipi0: frecpi0){
                if(iMother==ipi0){
                  isSingle=false ;
                  break ;
                }
              }
              break ;
            }
            iMother = mthr->GetMother() ;
          }
          if(isFromPi0 && isSingle){
            frecpi0.push_back(iMother) ;
          }
        }
      }
    }
    
    for(int iPID=0; iPID<4; iPID++){
      if(iPID==1 && !dispBit) continue;
      if(iPID==2 && !cpvBit) continue;
      if(iPID==3 && (!dispBit || !cpvBit)) continue;

      fRecPhot[iPID]++ ;
      if(isPhoton){
        fRecPhotTrue[iPID]++ ;
        if(isFromPi0){
          fRecPhotTruePi0[iPID]++ ;
          if(isSingle){
            fRecPhotTruePi0Single[iPID]++ ;
          }
        }
      }
    }

    AliCaloPhoton* Gamma = new((*fArrGamma)[nGamma++]) AliCaloPhoton(p.Px(),p.Py(),p.Pz(),clu->E());
    Gamma->SetPrimary(clu->GetLabelAt(0)) ;
    Gamma->SetDispBit(dispBit) ;
    Gamma->SetCPVBit(cpvBit) ;
    if(isPhoton)
      Gamma->SetPhoton(1);
    else 
      Gamma->SetPhoton(0);
  }

  //Fill Real and mixed distributions
  for(Int_t iGamma=0; iGamma<nGamma-1; iGamma++){
    AliCaloPhoton * gamma1 = (AliCaloPhoton*)fArrGamma->At(iGamma);
    for(Int_t jGamma=iGamma+1; jGamma<nGamma; jGamma++){
      AliCaloPhoton * gamma2 = (AliCaloPhoton*)fArrGamma->At(jGamma); 
      Double_t m12=(*gamma1 + *gamma2).M();
      for(int iPID=0; iPID<4; iPID++){
        if(iPID==1 && !(gamma1->IsDispOK() && gamma2->IsDispOK())) continue;
        if(iPID==2 && !(gamma1->IsCPVOK() && gamma2->IsDispOK())) continue;
        if(iPID==3 && !(gamma1->IsDispOK() && gamma2->IsDispOK() &&  gamma1->IsCPVOK() && gamma2->IsDispOK())) continue;
        fhReal[iPID]->Fill(m12,fCentrality) ;

        if(CommonAnsestor(gamma1->GetPrimary(),gamma2->GetPrimary())!=0){
          fhRealCommon[iPID]->Fill(m12,fCentrality);
        }          
        if(gamma1->IsPhoton() && gamma2->IsPhoton() )
          fhRealTrue[iPID]->Fill(m12,fCentrality);
      }
    }
  }


  // Mixed inv masses
  for(Int_t iGamma=0; iGamma<nGamma; iGamma++){
    AliCaloPhoton * gamma1 = (AliCaloPhoton*)fArrGamma->At(iGamma);
    TIter nextEv(fCurrentMixedList) ;
    while(TClonesArray * event2 = static_cast<TClonesArray*>(nextEv())){
      Int_t nGamma2 = event2->GetEntriesFast() ;
      for(Int_t jGamma=0; jGamma < nGamma2 ; jGamma++){
        AliCaloPhoton * gamma2 = static_cast<AliCaloPhoton*>(event2->At(jGamma)) ;
        double m12=(*gamma1 + *gamma2).M();
        for(int iPID=0; iPID<4; iPID++){
          if(iPID==1 && !(gamma1->IsDispOK() && gamma2->IsDispOK())) continue;
          if(iPID==2 && !(gamma1->IsCPVOK() && gamma2->IsDispOK())) continue;
          if(iPID==3 && !(gamma1->IsDispOK() && gamma2->IsDispOK() &&  gamma1->IsCPVOK() && gamma2->IsDispOK())) continue;
          fhMixed[iPID]->Fill(m12,fCentrality) ;
        }
      }
    }
  }

  fCurrentMixedList->AddFirst(fArrGamma); 
  fArrGamma=0x0;
  if(fCurrentMixedList->GetSize() > 20){
    TClonesArray *tmp = static_cast <TClonesArray*> (fCurrentMixedList->Last());
    fCurrentMixedList->Remove(tmp);
    delete tmp;
  }
}
//_____________________________________________________________________________
int AliAnalysisPHOSFluctuations::CommonAnsestor(int prim1,int prim2){

  //Find common ansestor except quarks and gluons

  if(!fStack) return 0 ;
  if(prim1==-1 || prim2==-1) return 0 ;

  while(prim1>-1){
    Int_t pr2=prim2 ; 
    while(pr2>-1){
      if(prim1==pr2){ 
        if((abs(((AliAODMCParticle*)fStack->At(prim1))->GetPdgCode()) <= 8) || 
           ((AliAODMCParticle*)fStack->At(prim1))->GetPdgCode() == 21)
          return 0;
        else
          return ((AliAODMCParticle*)fStack->At(prim1))->GetPdgCode() ;
      }
      AliAODMCParticle * part2 = (AliAODMCParticle*)fStack->At(pr2) ;
      pr2=part2->GetMother();
    }
    AliAODMCParticle * part1 = (AliAODMCParticle*)fStack->At(prim1) ;
    prim1=part1->GetMother();
  }
  return 0 ;
}
