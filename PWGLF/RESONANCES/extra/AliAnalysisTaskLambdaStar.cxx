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
#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TH1F.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TExMap.h>
#include <TVector3.h>
#include <TAxis.h>
#include <TRandom.h>
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliAODEvent.h>
#include <AliAODVertex.h>
#include <AliAODv0.h>
#include <AliAODInputHandler.h>
#include "AliAnalysisTaskLambdaStar.h"
#include <AliCentrality.h>
#include <AliPID.h>
#include <AliPIDResponse.h>
#include <AliExternalTrackParam.h>
#include <AliAODTrack.h>
#include <AliVTrack.h>


ClassImp(AliAnalysisTaskLambdaStar)

#endif 

//________________________________________________________________________
AliAnalysisTaskLambdaStar::AliAnalysisTaskLambdaStar() 
  : AliAnalysisTaskSE(),
  fkAbsZvertexCut(10.0),
  fkCentCut(80.0),
  fkKaonMass(0.49366),
  fkProMass(0.9382720),
  fPIDResponse(0),
  fCirc(0), 
  fCentMin(0), 
  fCentMax(0), 
  fNSigma(0),
  fRejNSigma(0),
  fClusterTPC(0),
  fDCAxy(0),
  fFilterBit(0),
  fNMix(0),
  fPatch(0),
  fCentPerPatch(0),
  fResoBuffer(0),
  fAOD(0), 
  fPrimaryVtx(0), 
  fOutputList(0), 
  fOutputPrimaries(0),
  fGTI(0),
  fTrackBuffSize(18000),
  fHistGoodEvent(0),
  fHistEvent(0),
  fHistZVertexCent(0),
  fPriHistShare(0),
  fHistMassPtPKMin(0),   
  fHistMassPtPbarKPlus(0),               
  fHistMassPtPKMinMix(0),   
  fHistMassPtPbarKPlusMix(0),               
  fHistMassPtPKPlusLS(0),   
  fHistMassPtPbarKMinLS(0),
  fPriHistTPCnsigmakaon_nocut(0),
  fPriHistTPCnsigmaproton_nocut(0),
  fPriHistTOFnsigmakaon_nocut(0),
  fPriHistTOFnsigmaproton_nocut(0),
  fPriHistTPCTOFnsigmakaon_nocut(0),
  fPriHistTPCTOFnsigmaproton_nocut(0),
  fPriHistTPCnsigmakaon(0),
  fPriHistTPCnsigmaproton(0),
  fPriHistTOFnsigmakaon(0),
  fPriHistTOFnsigmaproton(0),
  fPriHistTPCTOFnsigmakaon(0),
  fPriHistTPCTOFnsigmaproton(0),
  fPriHistTPCsignal(0),
  fPriHistTPCsignalproton(0),
  fPriHistTPCsignalkaon(0),
  fPriHistDCAxyYPtPro(0),           
  fPriHistDCAxyYPtAPro(0),            
  fPriHistDCAxyYPtKPlus(0),           
  fPriHistDCAxyYPtKMinus(0)   
  
{
  // Dummy constructor
  fPrimaryVtxPosition[0]=0;
  fPrimaryVtxPosition[1]=0;
  fPrimaryVtxPosition[2]=0;
}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::AliAnalysisTaskLambdaStar(const char *name) 
  : AliAnalysisTaskSE(name),
    fkAbsZvertexCut(10.0),
    fkCentCut(80.0),
    fkKaonMass(0.493),
    fkProMass(0.9382720),
    fPIDResponse(0), 
    fCirc(0), 
    fCentMin(0), 
    fCentMax(0), 
    fNSigma(0),
    fRejNSigma(0),
    fNMix(0),
    fClusterTPC(0),
    fDCAxy(0),
    fFilterBit(0),
    fPatch(0),
    fCentPerPatch(0),
    fResoBuffer(0),
    fAOD(0), 
    fPrimaryVtx(0), 
    fOutputList(0), 
    fOutputPrimaries(0),
    fGTI(0),
    fTrackBuffSize(18000),
    fHistGoodEvent(0),
    fHistEvent(0),
    fHistZVertexCent(0),
    fPriHistShare(0),
    fHistMassPtPKMin(0),   
    fHistMassPtPbarKPlus(0),               
    fHistMassPtPKMinMix(0),   
    fHistMassPtPbarKPlusMix(0),               
    fHistMassPtPKPlusLS(0),   
    fHistMassPtPbarKMinLS(0),
    fPriHistTPCnsigmakaon_nocut(0),
    fPriHistTPCnsigmaproton_nocut(0),
    fPriHistTOFnsigmakaon_nocut(0),
    fPriHistTOFnsigmaproton_nocut(0),
    fPriHistTPCTOFnsigmakaon_nocut(0),
    fPriHistTPCTOFnsigmaproton_nocut(0),
    fPriHistTPCnsigmakaon(0),
    fPriHistTPCnsigmaproton(0),
    fPriHistTOFnsigmakaon(0),
    fPriHistTOFnsigmaproton(0),
    fPriHistTPCTOFnsigmakaon(0),
    fPriHistTPCTOFnsigmaproton(0),
    fPriHistTPCsignal(0),
    fPriHistTPCsignalproton(0),
    fPriHistTPCsignalkaon(0),
    fPriHistDCAxyYPtPro(0),           
    fPriHistDCAxyYPtAPro(0),            
    fPriHistDCAxyYPtKPlus(0),           
    fPriHistDCAxyYPtKMinus(0)   
    
{
  // Constructor
  fPrimaryVtxPosition[0]=0;
  fPrimaryVtxPosition[1]=0;
  fPrimaryVtxPosition[2]=0;
  
  // Define output slots only here
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  
}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::~AliAnalysisTaskLambdaStar() {
  // Destructor, go through the data member and delete them
  // fPIDResponse is just a pointer to the pid response task,
  // we don't create it so we don't delete it. It comes from 
  // the AliInputEventHandler
  
  if(fResoBuffer){
    delete fResoBuffer;
    fResoBuffer=0;
  }
    
  // The lists containing the histograms
  if (fOutputList){
    fOutputList->Delete();
    delete fOutputList;
    fOutputList=0;
  }
  if (fOutputPrimaries){
    fOutputPrimaries->Delete();
    delete fOutputPrimaries;
    fOutputPrimaries=0;
  }
    
  // Array, note the [] with the delete
  if (fGTI)
    delete[] fGTI;
  fGTI=0;
  
}
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::UserCreateOutputObjects()
{
  // Create histograms and other objects and variables
  // Called once
  
  // Get the PID response object
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(!man){AliError("Couldn't get the analysis manager!");}
  
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(!inputHandler){AliError("Couldn't get the input handler!");}
  fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse){AliError("Couldn't get the PID response task!");}
  
  // Create the buffer for event mixing
  // Standard values are
  //  fkZvertexBins(10),
  //  fkCentBins(10),
  //  fkMixBuff(5),
  //  fkPriTrackLim(500),
  //  fkV0Lim(50),
 
  fResoBuffer = new ResoBuffer(10,10,fNMix,1200,fkAbsZvertexCut,fkCentCut);

  // In AODs, TPC only tracks don't have the pid information stored.
  // Also, the TPC only tracks don't have any resolution in the DCAxy
  // to distinguish between primaries and secondaries so we need the
  // corresponding global track for every TPC only track. The way to do 
  // this is to just store the pointer to the global track for every id.
  fGTI = new AliAODTrack *[fTrackBuffSize]; // Array of pointers 

  // Create the output list
  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputPrimaries = new TList();
  fOutputPrimaries->SetOwner();
  
  // Invariant mass binning for lambdas
  const Int_t nMinvBins = 300;
  const Float_t minvLowEdge=1.4, minvHiEdge=1.7;
  
  // Control hist for event cuts
  fHistGoodEvent = new TH1F("h1GoodEvent","No of events passing the cuts.",10,-.5,9.5);
  fHistEvent = new TH1F("hEvent", "Accepted Event Vs. Centrality", 100,0.0,100);
  fHistZVertexCent = new TH2D("fHistZVertexCent"," Vz  Vs Centrality; V_{Z} {cm} ; Centrality {%}", 60, -14.5, 14.5,100,0,100);
  // 3d y pt mass
  fHistMassPtPKMin = new TH2D ("InvMassPtPKMin","M_{inv}pt Vs mass",nMinvBins,minvLowEdge,minvHiEdge,300,0.,30.);  
  fHistMassPtPbarKPlus = new TH2D ("InvMassPtPKPlus","M_{inv}pt Vs mass",nMinvBins,minvLowEdge,minvHiEdge ,300, 0., 30.0);
  fHistMassPtPKMinMix = new TH2D ("InvMassPtPKMinMix","M_{inv}pt Vs mass",nMinvBins,minvLowEdge,minvHiEdge ,300, 0., 30.);
  fHistMassPtPbarKPlusMix = new TH2D ("InvMassPtPKPlusMix","M_{inv}pt Vs mass",nMinvBins,minvLowEdge,minvHiEdge , 300, 0., 30.0);
  fHistMassPtPKPlusLS = new TH2D ("InvMassPtPKMinLS","M_{inv}pt Vs mass",nMinvBins,minvLowEdge,minvHiEdge ,300, 0.0, 30.);
  fHistMassPtPbarKMinLS = new TH2D ("InvMassPtPKPlusLS","M_{inv}pt Vs mass",nMinvBins,minvLowEdge,minvHiEdge ,300,0.0,30.0);
 
  fOutputList->Add(fHistGoodEvent);
  fOutputList->Add(fHistEvent);
  fOutputList->Add(fHistZVertexCent);
  fOutputList->Add(fHistMassPtPKMin);
  fOutputList->Add(fHistMassPtPbarKPlus);
  fOutputList->Add(fHistMassPtPKMinMix);
  fOutputList->Add(fHistMassPtPbarKPlusMix);
  fOutputList->Add(fHistMassPtPKPlusLS);
  fOutputList->Add(fHistMassPtPbarKMinLS);
  
  //nsigma plots before nsigma selection
  fPriHistTPCnsigmakaon_nocut = new TH2D("nsigmaTPCkaonVsPt_nocut", "nsigmaTPC for kaon before nsigma selection",300,0.0,10.0, 200,-10.0,10.0);
  fPriHistTPCnsigmaproton_nocut = new TH2D("nsigmaTPCprotonVsPt_nocut", "nsigmaTPC for proton before nsigma selection",300,0.0,10.0, 200,-10.0,10.0);
  fPriHistTOFnsigmakaon_nocut = new TH2D("nsigmaTOFkaonVsPt_nocut", "nsigmaTOF for kaon before nsigma selection",300,0.0,10.0, 200,-10.0,10.0);
  fPriHistTOFnsigmaproton_nocut = new TH2D("nsigmaTOFprotonVsPt_nocut", "nsigmaTOF for proton before nsigma selection",300,0.0,10.0, 200,-10.0,10.0);
  fPriHistTPCTOFnsigmakaon_nocut = new TH2D("nsigmaTPCTOFkaonVsPt_nocut", "nsigmaTPCTOF for kaon before nsigma selection",200,-10.0,10.0, 200,-10.0,10.0);
  fPriHistTPCTOFnsigmaproton_nocut = new TH2D("nsigmaTPCTOFprotonVsPt_nocut", "nsigmaTPCTOF for proton before nsigma selection",200,-10.0,10.0, 200,-10.0,10.0 );
  
  //nsigma plots
  fPriHistTPCnsigmakaon = new TH2D("nsigmaTPCkaonVsPt", "nsigmaTPC for kaon",500,0.0,10.0, 200,-10.0,10.0);
  fPriHistTPCnsigmaproton = new TH2D("nsigmaTPCprotonVsPt", "nsigmaTPC for proton",500,0.0,10.0, 200,-10.0,10.0);
  fPriHistTOFnsigmakaon = new TH2D("nsigmaTOFkaonVsPt", "nsigmaTOF for kaon ",500,0.0,10.0, 200,-10.0,10.0);
  fPriHistTOFnsigmaproton = new TH2D("nsigmaTOFprotonVsPt", "nsigmaTOF for proton ",100,0.0,10.0, 200,-10.0,10.0);
  fPriHistTPCTOFnsigmakaon = new TH2D("nsigmaTPCTOFkaonVsPt", "nsigmaTPCTOF for kaon",200,-10.0,10.0, 200,-10.0,10.0);
  fPriHistTPCTOFnsigmaproton = new TH2D("nsigmaTPCTOFprotonVsPt", "nsigmaTPCTOF for proton",200,-10.0,10.0, 200,-10.0,10.0);
  
  // Shared clusters
  fPriHistShare = new TH1F ("h1PriShare","Shared clusters, primaries;#shared clusters;counts", 160,0,160);
  // dEdx analysis
  fPriHistTPCsignal = new TH2F ("h2TPCsignal","TPC signal for all tracks ;p_{tot};dEdx",500,0.0,10.0,500,0.0,1000);
  fPriHistTPCsignalkaon = new TH2F ("h2TPCsignalkaon","TPC signal for kaons;p_{tot};dEdx",500,0.0,10.0,500,0.0,1000);
  fPriHistTPCsignalproton = new TH2F ("h2TPCsignalproton","TPC signal for protons;p_{tot};dEdx",500,0.0,10.0,500,0.0,1000);
    
    
  //  Common for all protons - DCA xy distribution to determine primaries, secondaries from weak decay and secondaries from material
  fPriHistDCAxyYPtPro = new TH3F ("h3DCAxyYPtPro","DCAxy vs (y,pt) protons",100,-3.,3.,30,-1.5,1.5,100,0.,30);
  fPriHistDCAxyYPtAPro = new TH3F ("h3DCAxyYPtAPro","DCAxy vs (y,pt) anti-protons",100,-3.,3.,30,-1.5,1.5,100,0.,30);
  fPriHistDCAxyYPtKPlus = new TH3F ("h3DCAxyYPtKPlus","DCAxy vs (y,pt) kplus",100,-3.,3.,30,-1.5,1.5,100,0.,30.);
  fPriHistDCAxyYPtKMinus = new TH3F ("h3DCAxyYPtKMinus","DCAxy vs (y,pt) kminus",100,-3.,3.,30,-1.5,1.5,100,0.,30.);  
  
  fOutputPrimaries->Add(fPriHistTPCnsigmakaon_nocut);
  fOutputPrimaries->Add(fPriHistTPCnsigmaproton_nocut);
  fOutputPrimaries->Add(fPriHistTOFnsigmakaon_nocut);
  fOutputPrimaries->Add(fPriHistTOFnsigmaproton_nocut);
  fOutputPrimaries->Add(fPriHistTPCTOFnsigmakaon_nocut);
  fOutputPrimaries->Add(fPriHistTPCTOFnsigmaproton_nocut);
  fOutputPrimaries->Add(fPriHistTPCnsigmakaon);
  fOutputPrimaries->Add(fPriHistTPCnsigmaproton);
  fOutputPrimaries->Add(fPriHistTOFnsigmakaon);
  fOutputPrimaries->Add(fPriHistTOFnsigmaproton);
  fOutputPrimaries->Add(fPriHistTPCTOFnsigmakaon);
  fOutputPrimaries->Add(fPriHistTPCTOFnsigmaproton);
  fOutputPrimaries->Add(fPriHistShare);
  fOutputPrimaries->Add(fPriHistTPCsignal);
  fOutputPrimaries->Add(fPriHistTPCsignalkaon);
  fOutputPrimaries->Add(fPriHistTPCsignalproton);
  fOutputPrimaries->Add(fPriHistDCAxyYPtPro);
  fOutputPrimaries->Add(fPriHistDCAxyYPtAPro);
  fOutputPrimaries->Add(fPriHistDCAxyYPtKPlus);
  fOutputPrimaries->Add(fPriHistDCAxyYPtKMinus);
  
  // Post the data
  PostData(1, fOutputList);
  PostData(2, fOutputPrimaries);
  
}

//________________________________________________________________________
void AliAnalysisTaskLambdaStar::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event and Fill a control histogram
  fHistGoodEvent->Fill(0.0);
  // Get the event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    return;
  }
  
  // Fill a control histogram
  fHistGoodEvent->Fill(1.0);  
  
  // Get the centrality selection
  AliCentrality *centrality=NULL;
  centrality = fAOD->GetCentrality();
  if (!centrality) {
    printf ("ERROR: couldn't get the AliCentrality\n");
    return;
  }
  
  // Fill a control histogram
  fHistGoodEvent->Fill(2.0);
  
  if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected())) return;
  fHistGoodEvent->Fill(3.0);
  
  // Analyze only events using multiplicity in V0 detector (standard)
  Float_t centralityPercentile = centrality->GetCentralityPercentileUnchecked("V0M");
  if(!centrality->GetCentralityPercentileUnchecked("V0M"))
    {return ;}
  fHistGoodEvent->Fill(4.0);
  
  if ( centralityPercentile <= fCentMin || centralityPercentile > fCentMax){return;}
  if ( centralityPercentile >= 5 && centralityPercentile <= 20){  
    ApplyCentralityPatchPbPb2011(centrality);  }
  
  // Fill a control histogram
  fHistGoodEvent->Fill(5.0);
  
  // Primary vertex, GetPrimaryVertex() returns the "best" reconstructed vertex
  fPrimaryVtx = fAOD->GetPrimaryVertex();
  if (!fPrimaryVtx){
    printf ("ERROR: no primary vertex\n");
    return;
  }
  // Fill a control histogram
  fHistGoodEvent->Fill(6.0);
  
  fPrimaryVtx->GetXYZ(fPrimaryVtxPosition);  
  if (TMath::Abs(fPrimaryVtxPosition[2]) > fkAbsZvertexCut)
    return;
  // Fill a control histogram
  fHistGoodEvent->Fill(7.0);
  //fill centrality histogram
  fHistEvent->Fill(centralityPercentile);
  
  //Fill  Z-vertex histo
  fHistZVertexCent->Fill(fPrimaryVtxPosition[2], centralityPercentile);
  
  // Multiplicity
  if (!(fAOD->GetNumberOfTracks())) {
    return;
  }
  
  // Fill a control histogram
  fHistGoodEvent->Fill(8.0);
  
  // Set up the event buffer to store this event
  fResoBuffer->ShiftAndAdd(fAOD);
  
  // Reset the reference array to the global tracks..
  ResetGlobalTrackReference();
  
  //  AliAODTrack *track=NULL;
   // cout<<"I am In .....the event loop."<<endl;
   // cout<<"rejected nsigma = "<<fRejNSigma<<endl;
  for (Int_t iTrack=0;iTrack<fAOD->GetNumberOfTracks();iTrack++){
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if (!track) continue;
    // Store the reference of the global tracks
    StoreGlobalTrackReference(track);  
  }
  
  for (Int_t iTrack=0;iTrack<fAOD->GetNumberOfTracks();iTrack++){
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if (!track) continue;
    if(!track->TestFilterBit(fFilterBit))  
      continue;
    if(!AcceptTrack(track))
      continue;
    // Reject tracks with shared clusters
    if(!GoodTPCFitMapSharedMap(track))
      continue;
    
    // Visualization of TPC dE/dx
    FillDedxHist(track);
    
    // Depending on momentum choose pid method
    if (track->P() < 0.5){
      ProcessTPC(track,fNSigma);
    }
      else if (track->P() < 0.7){
       ProcessHybridPro(track,fCirc,fNSigma);
       }
   else if (track->P() < 100.0) {
      ProcessHybrid(track,fCirc,fNSigma);
    }
    
  } // End of loop over primary tracks
  
  
  // Process real events
  ProcessReal( );
  ProcessLikeSignBkg(); 
  // Process mixed events
  ProcessMixed();
  
  // Post output data.
  PostData(1, fOutputList);
  PostData(2, fOutputPrimaries);

}


//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ProcessTPC(AliAODTrack* track, Double_t nsig){
  
  Float_t TMom = track->Pt();
  Double_t nsigmapion = 999, nsigmakaon=999,nsigmaproton=999,nsigmaelectron=999;
  
  nsigmaelectron = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron) ;
  nsigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon) ;
  nsigmapion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion) ;
  nsigmaproton =fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton) ;
  
  fPriHistTPCnsigmakaon_nocut->Fill(TMom,nsigmakaon);
  fPriHistTPCnsigmaproton_nocut->Fill(TMom,nsigmaproton);
  
//  if(nsigmaelectron  < 3.0  &&  nsigmapion   > 3.0  && nsigmakaon   > 3.0
//     && nsigmaproton > 3.0 ) return; // DANGEROUS LINE
    
//  if( ( nsigmakaon==nsigmapion ) && ( nsigmakaon==nsigmaproton )) return ; // USELESS LINE never happens
  
    // KAON IDENTIFICATION
//  if( ( nsigmakaon   <  nsigmapion ) && ( nsigmakaon < nsigmaproton ) && (nsigmakaon <= nsig) ){ // NOT OKAY, only (nsigmakaon <= nsig) or if you want to reject pions and protons do (nsigmakaon <= nsig) && (nsigmapion >= rejnsig) && (nsigmaproton >= rejnsig)
    
    
  if( (TMath::Abs(nsigmakaon) <= nsig) && (TMath::Abs(nsigmapion) >= fRejNSigma) && (TMath::Abs(nsigmaproton) >= fRejNSigma) ) {  // CLEAN KAON IDENTIFICATION
    
    fPriHistTPCsignalkaon->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
    fPriHistTPCnsigmakaon->Fill(TMom,nsigmakaon);
      
    if (track->Charge()>0){  
      // Cut .1 cm on DCAxy and fill a histogram
      if(goodDCAKaon(track)){
        // Add to the  event
	fResoBuffer->GetEvt(0)->AddKPlus(track);
      }
    }
    else{
      // Cut .1 cm on DCAxy and fill a histogram    
      if(goodDCAKaon(track)){
	// Add to the  event
	fResoBuffer->GetEvt(0)->AddKMin(track);
      }
    }
  }

 // if( ( nsigmaproton==nsigmapion ) && ( nsigmaproton==nsigmakaon )) return ; // USELESS LINE
  
    // PROTON IDENTIFICATION
//  if( ( nsigmaproton   < nsigmapion ) && ( nsigmaproton < nsigmakaon ) && (nsigmaproton   <= nsig )){ // NOT OK, see above
    
    if( (TMath::Abs(nsigmaproton) <= nsig ) && (TMath::Abs(nsigmakaon) >= fRejNSigma) && (TMath::Abs(nsigmapion) >= fRejNSigma) ){ // CLEAN PROTON IDENTIFICATION

        fPriHistTPCsignalproton->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
        fPriHistTPCnsigmaproton->Fill(TMom,nsigmaproton);

    if (track->Charge()>0){
      // Cut .1 cm on DCAxy and fill a histogram
      if(goodDCA(track)){
	// Add to the  event
	fResoBuffer->GetEvt(0)->AddPro(track);
      }
    }
    else{
      // Cut .1 cm on DCAxy and fill a histogram
      if(goodDCA(track)){
	// Add to the  event
	fResoBuffer->GetEvt(0)->AddAPro(track);
      }
    }
  }
  
} // End of void ProcessTPC


//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ProcessHybridPro(AliAODTrack *track, Bool_t circ, Double_t nsig){
  Float_t TMom = track->Pt();
  Double_t nsigmapion = 999, nsigmakaon=999,nsigmaproton=999,nsigmaelectron=999;
  
  nsigmaelectron = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron) ;
  nsigmapion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion) ;
  nsigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon) ;
  nsigmaproton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton) ;
  
  fPriHistTPCnsigmakaon_nocut->Fill(TMom,nsigmakaon);
  fPriHistTPCnsigmaproton_nocut->Fill(TMom,nsigmaproton);
  
  Double_t nsigmaprotonTOF=999.,nsigmakaonTOF=999.,nsigmapionTOF=999.;
  Double_t nsigmaTPCTOFkProton=999.,nsigmaTPCTOFkKaon=999.,nsigmaTPCTOFkPion=999.;
  Bool_t fHasTOFPID;
  
  //if( nsigmaelectron  < 3.0  &&  nsigmapion  > 3.0  && nsigmakaon  > 3.0  && nsigmaproton > 3.0 )  return; //   DANGEROUS LINE
  
  if(track->GetStatus() & AliVTrack::kTOFpid){ fHasTOFPID=kTRUE; }
  else{ fHasTOFPID=kFALSE; }
  
  
  if (fHasTOFPID){
    nsigmapionTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion) ;
    nsigmakaonTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon) ;
    nsigmaprotonTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton) ;
    
    fPriHistTOFnsigmaproton_nocut->Fill(TMom,nsigmaprotonTOF);
    fPriHistTPCTOFnsigmaproton_nocut->Fill(nsigmaprotonTOF,nsigmaproton);
    fPriHistTOFnsigmakaon_nocut->Fill(TMom,nsigmakaonTOF);
    fPriHistTPCTOFnsigmakaon_nocut->Fill(nsigmakaonTOF,nsigmakaon);
    
    
    if(circ){   
      Double_t d2Proton=nsigmaproton * nsigmaproton + nsigmaprotonTOF * nsigmaprotonTOF;
      Double_t d2Kaon=nsigmakaon * nsigmakaon + nsigmakaonTOF * nsigmakaonTOF;
      Double_t d2Pion=nsigmapion * nsigmapion + nsigmapionTOF * nsigmapionTOF;
      
      nsigmaTPCTOFkProton  =  TMath::Sqrt(d2Proton);
      nsigmaTPCTOFkKaon    =  TMath::Sqrt(d2Kaon);
      nsigmaTPCTOFkPion    =  TMath::Sqrt(d2Pion);
    }
    else{  
      nsigmaTPCTOFkProton  =  nsigmaprotonTOF;
      nsigmaTPCTOFkKaon    =  nsigmakaonTOF;
      nsigmaTPCTOFkPion    =  nsigmapionTOF;
    }
  }
  else{
    nsigmaTPCTOFkProton = nsigmaproton;
    nsigmaTPCTOFkKaon   = nsigmakaon;
    nsigmaTPCTOFkPion   = nsigmapion;
  }
  
  
  /*if(circ){ if( (nsigmaTPCTOFkKaon >= nsigmaTPCTOFkPion ) && ( nsigmaTPCTOFkKaon >= nsigmaTPCTOFkProton )) return ;  } // RANDOM Rejection not OkAY
  else{ 
    if( (nsigmakaon >= nsigmapion ) && ( nsigmakaon >= nsigmaproton )) return ; RANDOM Rejection not OkAY
    if(fHasTOFPID) {
      if( (nsigmakaonTOF >= nsigmapionTOF ) && ( nsigmakaonTOF >= nsigmaprotonTOF )) return ; }	
      }*/ //USELESS line
  
  
  if(!circ){
    if(fHasTOFPID){ //PURE Kaon selection by TPC and TOF
        if( ( TMath::Abs(nsigmaTPCTOFkKaon)   <= nsig ) && (TMath::Abs(nsigmaTPCTOFkPion) >= fRejNSigma) && (TMath::Abs(nsigmaTPCTOFkProton) >= fRejNSigma)  && ( TMath::Abs(nsigmakaon) <= 5.0) ) {
            
	fPriHistTPCnsigmakaon->Fill(TMom,nsigmakaon);
	fPriHistTOFnsigmakaon->Fill(TMom,nsigmakaonTOF);
	fPriHistTPCTOFnsigmakaon->Fill(nsigmakaonTOF,nsigmakaon);
    fPriHistTPCsignalkaon->Fill(track->GetTPCmomentum(),track->GetTPCsignal());

	// Distinguish between charges
	if (track->Charge() > 0) {
	  // Cut .1 cm on DCAxy and fill a histogram
	  if(goodDCAKaon(track)){
	    // Add to the  event
	    fResoBuffer->GetEvt(0)->AddKPlus(track);
	  }
    }
	else {
	  // Cut .1 cm on DCAxy and fill a histogram
	  if(goodDCAKaon(track)){
	    // add to the  event
	    fResoBuffer->GetEvt(0)->AddKMin(track);
	  }
	}
      }
    }
    else  { //PURE Kaon selection by TPC if TOF is not present
        if( (TMath::Abs(nsigmakaon) <= nsig) && (TMath::Abs(nsigmaproton) >= fRejNSigma) && (TMath::Abs(nsigmapion) >= fRejNSigma)){
            
            fPriHistTPCsignalkaon->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
            fPriHistTPCnsigmakaon->Fill(TMom,nsigmakaon);

	// Distinguish between charges
	if (track->Charge() > 0) {
	  // Cut .1 cm on DCAxy and fill a histogram
	  if(goodDCAKaon(track)){
	    // Add to the  event
	    fResoBuffer->GetEvt(0)->AddKPlus(track);
	  }
	}
	else {
	  // Cut .1 cm on DCAxy and fill a histogram
	  if(goodDCAKaon(track)){
	    // add to the  event
	    fResoBuffer->GetEvt(0)->AddKMin(track);
	  }
	}
      }
    }    
  } //end of  if(!circ)
  
  else { // circular pid PURE Kaon selection by TPC or TPC||TOF by circular nsigma cuts
      if( (TMath::Abs(nsigmaTPCTOFkKaon)  <= nsig ) && (TMath::Abs(nsigmaTPCTOFkProton) >= fRejNSigma) && (TMath::Abs(nsigmaTPCTOFkPion) >= fRejNSigma) ){
          
          fPriHistTPCsignalkaon->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
          fPriHistTPCnsigmakaon->Fill(TMom,nsigmakaon);
         fPriHistTOFnsigmakaon->Fill(TMom,nsigmakaonTOF);
         fPriHistTPCTOFnsigmakaon->Fill(nsigmakaonTOF,nsigmakaon);
      if (track->Charge() > 0) {
	// Cut .1 cm on DCAxy and fill a histogram
	if(goodDCAKaon(track)){
	  // Add to the  event
	  fResoBuffer->GetEvt(0)->AddKPlus(track);
	}
      }
      else {
	// Cut .1 cm on DCAxy and fill a histogram
	if(goodDCAKaon(track)){
	  // add to the  event
	  fResoBuffer->GetEvt(0)->AddKMin(track);
	}
      }
    }
  }
  
  //Proton selection
  // if( ( nsigmaproton == nsigmapion ) && ( nsigmaproton == nsigmakaon )) return ; // USELESS LINE
  // if( (nsigmaproton < nsigmapion) && ( nsigmaproton < nsigmakaon) && ( nsigmaproton   <= nsig) ){ //// USELESS LINE Random Rejection

    if(  ( TMath::Abs(nsigmaproton)   <= nsig) &&( TMath::Abs(nsigmapion) >= fRejNSigma) && ( TMath::Abs(nsigmakaon) >= fRejNSigma)  ){  // CLEAN Proton Selection by TPC
    fPriHistTPCnsigmaproton->Fill(TMom,nsigmaproton);
        fPriHistTPCsignalproton->Fill(track->GetTPCmomentum(),track->GetTPCsignal());

    // Distinguish between charges
    if (track->Charge() > 0) {
      // Cut .1 cm on DCAxy and fill a histogram
      if(goodDCA(track)){	
	// Add to the  event
	fResoBuffer->GetEvt(0)->AddPro(track);
      }
    }
    else {
      // Cut .1 cm on DCAxy and fill a histogram
      if(goodDCA(track)){
	// add to the  event
	fResoBuffer->GetEvt(0)->AddAPro(track);
      }
    }  
  }
  
} // End of ProcessHybridPro

//________________________________________________________________________

void AliAnalysisTaskLambdaStar::ProcessHybrid(AliAODTrack *track, Bool_t circ,Double_t nsig){
  Float_t TMom = track->Pt();
  Double_t nsigmapion = 999, nsigmakaon=999,nsigmaproton=999,nsigmaelectron=999;
  
  nsigmaelectron = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron) ;
  nsigmapion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion) ;
  nsigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon) ;
  nsigmaproton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton) ;
  
  fPriHistTPCnsigmakaon_nocut->Fill(TMom,nsigmakaon);
  fPriHistTPCnsigmaproton_nocut->Fill(TMom,nsigmaproton);
  
  Double_t nsigmaprotonTOF=999.,nsigmakaonTOF=999.,nsigmapionTOF=999.;
  Double_t nsigmaTPCTOFkProton=999.,nsigmaTPCTOFkKaon=999.,nsigmaTPCTOFkPion=999.;
  Bool_t fHasTOFPID;
  
  //if(nsigmaelectron  < 3.0 && nsigmapion   > 3.0  && nsigmakaon   > 3.0 && nsigmaproton > 3.0 ) return ; // DANGEROUS LINE
  
  if(track->GetStatus() & AliVTrack::kTOFpid){ fHasTOFPID=kTRUE; }
  else{ fHasTOFPID=kFALSE; }
   
  if (fHasTOFPID){  
    nsigmakaonTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon) ;
    nsigmaprotonTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton) ;
    nsigmapionTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion) ;
    
    fPriHistTOFnsigmaproton_nocut->Fill(TMom,nsigmaprotonTOF);
    fPriHistTPCTOFnsigmaproton_nocut->Fill(nsigmaprotonTOF,nsigmaproton);
    fPriHistTOFnsigmakaon_nocut->Fill(TMom,nsigmakaonTOF);
    fPriHistTPCTOFnsigmakaon_nocut->Fill(nsigmakaonTOF,nsigmakaon);
    
    if(circ){  
      Double_t d2Proton=nsigmaproton * nsigmaproton + nsigmaprotonTOF * nsigmaprotonTOF;
      Double_t d2Kaon=nsigmakaon * nsigmakaon + nsigmakaonTOF * nsigmakaonTOF;
      Double_t d2Pion=nsigmapion * nsigmapion + nsigmapionTOF * nsigmapionTOF;
      
      nsigmaTPCTOFkProton  =  TMath::Sqrt(d2Proton);
      nsigmaTPCTOFkKaon    =  TMath::Sqrt(d2Kaon);
      nsigmaTPCTOFkPion    =  TMath::Sqrt(d2Pion);
    }
    else{
      nsigmaTPCTOFkProton  =  nsigmaprotonTOF;
      nsigmaTPCTOFkKaon    =  nsigmakaonTOF;
      nsigmaTPCTOFkPion    =  nsigmapionTOF;
    }
  }
  else {
    nsigmaTPCTOFkProton = nsigmaproton;
    nsigmaTPCTOFkKaon   = nsigmakaon;
    nsigmaTPCTOFkPion   = nsigmapion;
  }
  

  
  if(!circ){
    if(fHasTOFPID){
      if( ( TMath::Abs(nsigmaTPCTOFkKaon)   <= nsig ) && (TMath::Abs(nsigmaTPCTOFkProton) >= fRejNSigma) && (TMath::Abs(nsigmaTPCTOFkPion) >= fRejNSigma)  && ( TMath::Abs(nsigmakaon) <= 5.0) ){
          
    fPriHistTPCsignalkaon->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
	fPriHistTPCnsigmakaon->Fill(TMom,nsigmakaon);
	fPriHistTOFnsigmakaon->Fill(TMom,nsigmakaonTOF);
	fPriHistTPCTOFnsigmakaon->Fill(nsigmakaonTOF,nsigmakaon);

	// Distinguish between charges
	if (track->Charge() > 0) {
	  // Cut .1 cm on DCAxy and fill a histogram
	  if(goodDCAKaon(track)){
	    // Add to the  event
	    fResoBuffer->GetEvt(0)->AddKPlus(track);
	  }
	}
	else {
	  // Cut .1 cm on DCAxy and fill a histogram
	  if(goodDCAKaon(track)){
	    // add to the  event
	    fResoBuffer->GetEvt(0)->AddKMin(track);
	  }
	}
      }
    }
    else if (track->P() < 0.8) {
        if( (TMath::Abs(nsigmakaon) <= nsig) && (TMath::Abs(nsigmaproton) >= fRejNSigma) && (TMath::Abs(nsigmapion) >= fRejNSigma)){
           
            fPriHistTPCsignalkaon->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
            fPriHistTPCnsigmakaon->Fill(TMom,nsigmakaon);

	// Distinguish between charges
	if (track->Charge() > 0) {
	  // Cut .1 cm on DCAxy and fill a histogram
	  if(goodDCAKaon(track)){
	    // Add to the  event
	    fResoBuffer->GetEvt(0)->AddKPlus(track);
	  }
    }
	else {
	  // Cut .1 cm on DCAxy and fill a histogram
	  if(goodDCAKaon(track)){
	    // add to the  event
	    fResoBuffer->GetEvt(0)->AddKMin(track);
	  }
	}
      }
    }
    
  } //end of if(!circ)
  
  else  { // for circular
      if(  (TMath::Abs(nsigmaTPCTOFkKaon)  <= nsig ) && (TMath::Abs(nsigmaTPCTOFkProton) >= fRejNSigma) && (TMath::Abs(nsigmaTPCTOFkPion) >= fRejNSigma) ){
          
      fPriHistTPCsignalkaon->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
      fPriHistTPCnsigmakaon->Fill(TMom,nsigmakaon);

      if (fHasTOFPID){
	fPriHistTOFnsigmakaon->Fill(TMom,nsigmakaonTOF);
	fPriHistTPCTOFnsigmakaon->Fill(nsigmakaonTOF,nsigmakaon);
      }
      if (track->Charge() > 0) {
	if(goodDCAKaon(track)){
	  fResoBuffer->GetEvt(0)->AddKPlus(track);
	}
      }
      else {
	// Cut .1 cm on DCAxy and fill a histogram
	if(goodDCAKaon(track)){
	  // add to the  event
	  fResoBuffer->GetEvt(0)->AddKMin(track);
	}
      }
    }
  }
  
    
  //  PROTON selection
  if(!circ){
    if(fHasTOFPID){
        if( ( TMath::Abs(nsigmaTPCTOFkProton) <= nsig ) && (TMath::Abs(nsigmaTPCTOFkPion) >= fRejNSigma )&& (TMath::Abs(nsigmaTPCTOFkKaon)  >= fRejNSigma )  && ( TMath::Abs(nsigmaproton) <= 5.0 ) ){ //Clean Proton selection by TPC and TOF
            
    fPriHistTPCsignalproton->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
	fPriHistTPCnsigmaproton->Fill(TMom,nsigmaproton);
	fPriHistTOFnsigmaproton->Fill(TMom,nsigmaprotonTOF);
	fPriHistTPCTOFnsigmaproton->Fill(nsigmaprotonTOF,nsigmaproton);
	if (track->Charge() > 0) {
	  if(goodDCA(track)){
	    fResoBuffer->GetEvt(0)->AddPro(track);
	  }
	}
	else {
	  if(goodDCA(track)){
	    fResoBuffer->GetEvt(0)->AddAPro(track);
	  }
	}
      }
    }
    else if(track->P() < 1.4){
      if( (TMath::Abs(nsigmaproton) <= nsig)  &&  (TMath::Abs(nsigmapion) >= fRejNSigma) && (TMath::Abs(nsigmakaon) >= fRejNSigma) ){
	
          fPriHistTPCsignalproton->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
          fPriHistTPCnsigmaproton->Fill(TMom,nsigmaproton);

	if (track->Charge() > 0) {
	  if(goodDCA(track)){
	    fResoBuffer->GetEvt(0)->AddPro(track);
	  }
	}
	else {	  
	  if(goodDCA(track)){
	    fResoBuffer->GetEvt(0)->AddAPro(track);
	  }
	}  
      }
    }
  }
  
  else
    {
      if( (TMath::Abs(nsigmaTPCTOFkProton)  <= nsig) && (TMath::Abs(nsigmaTPCTOFkPion) >= fRejNSigma) && (TMath::Abs(nsigmaTPCTOFkKaon)  >= fRejNSigma )){
          
    fPriHistTPCsignalproton->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
	fPriHistTPCnsigmaproton->Fill(TMom,nsigmaproton);

	if (fHasTOFPID){
	  fPriHistTOFnsigmaproton->Fill(TMom,nsigmaprotonTOF);
	  fPriHistTPCTOFnsigmaproton->Fill(nsigmaprotonTOF,nsigmaproton);
	}
	if (track->Charge() > 0) {
	  if(goodDCA(track)){
	    fResoBuffer->GetEvt(0)->AddPro(track);
	  }
	}
	else {
	  if(goodDCA(track)){
	    fResoBuffer->GetEvt(0)->AddAPro(track);	    
	  }
	}
      }
    }
  
  
} // End of ProcessHybrid


Double_t AliAnalysisTaskLambdaStar::ApplyCentralityPatchPbPb2011( AliCentrality *central){
  //This part rejects randomly events such that the centrality gets flat for LHC11h Pb-Pb data
  //for 0-5% and 10-20% centrality bin
  Double_t cent = (Float_t)(central->GetCentralityPercentile("V0M"));
  Double_t rnd_hc, testf, ff, N1, N2;
  
  if(fCentPerPatch==510){
    
    N1 = 1.9404e+06;
    N2 = 1.56435e+06;
    ff = 5.04167e+06 - 1.49885e+06*cent + 2.35998e+05*cent*cent -1.22873e+04*cent*cent*cent;
  }
  
  if(fCentPerPatch==1020){
    N1 = 1.56435e+06;
    N2 = 4.20e+05;
    ff = 1.68062e+08 - 5.19673e+07*cent + 6.4068e+06*cent*cent + 6.4068e+06*cent*cent*cent - 392687*cent*cent*cent*cent - 145.07*cent*cent*cent*cent*cent;
  }
  
  testf = ( N2 + (N1-ff) ) / N1;
  rnd_hc = gRandom->Rndm();
  
  
  if (rnd_hc < 0 || rnd_hc > 1 )
    {
      AliWarning("Wrong Random number generated");
      return -999.0;
    }
  
  if (rnd_hc < testf)
    return cent;
  else
    return -999.0;
}



//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ProcessReal() {
  // Process real events
  
  Int_t iPro,iKminus,iAPro,iKplus;
  
  // Proton K- loop
  Int_t nproton = fResoBuffer->GetEvt(0)->GetNPro();
  Int_t nkmin = fResoBuffer->GetEvt(0)->GetNKMin();
  
  for (iPro = 0; iPro < nproton; iPro++){
    // Skip if unUseIt() entry
    if (!fResoBuffer->GetEvt(0)->fProTracks[iPro].UseIt())
      continue;
    // Kminus loop
    for (iKminus=0;iKminus < nkmin;iKminus++){
      
      // Skip if unUseIt() entry
      if (!fResoBuffer->GetEvt(0)->fKMinTracks[iKminus].UseIt())
  	continue;
      
      
      Double_t  pairrap =  Rapidity(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(0)->fKMinTracks[iKminus]);
      Double_t  invmass = MInv(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(0)->fKMinTracks[iKminus]);
      Double_t  pairpt  = Pt(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(0)->fKMinTracks[iKminus]);
      //      Double_t  ctheta1  = Costheta(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(0)->fKMinTracks[iKminus]);
      // Double_t  ctheta2  = Costheta1(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(0)->fKMinTracks[iKminus]);
      //Double_t  openang =  OpeningAngle(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(0)->fKMinTracks[iKminus]);
      
      if(TMath::Abs(pairrap) > 0.5) continue;
      
      //if(openang < 0.4) continue;
      //  if(TMath::Abs(ctheta1) > 0.8 && TMath::Abs(ctheta2 > 0.8)) continue;	       
      //      cout<<"*****openang"<<openang<<"********"<<ctheta1<<"****"<<ctheta2<<endl;
      //cout<<"rap"<<pairrap<<"invmass"<<invmass<<"pairpt"<<pairpt<<endl;

      fHistMassPtPKMin->Fill(invmass,pairpt);
      
    }// Kaon loop
    
  }// Proton loop
  
  Int_t npbar   = fResoBuffer->GetEvt(0)->GetNAPro();
  Int_t nkplus  = fResoBuffer->GetEvt(0)->GetNKPlus();

  
  for (iAPro = 0; iAPro<npbar; iAPro++){
    
    // Skip if unUseIt() entry
    
    if (!fResoBuffer->GetEvt(0)->fAProTracks[iAPro].UseIt())  continue;
    
    // Kplus loop
    
    for (iKplus=0;iKplus< nkplus;iKplus++){
      
      
      if (!fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus].UseIt()) continue;
      
      Double_t  pairrap =  Rapidity(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus]);
      Double_t  invmass = MInv(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus]);
      Double_t  pairpt  = Pt(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus]);
      //   Double_t  openang  = OpeningAngle(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus]);
      //Double_t  ctheta1  = Costheta(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus]);
      //Double_t  ctheta2  = Costheta1(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus]);
      
      
      if(TMath::Abs(pairrap) > 0.5) continue;
      
      
      //    if(openang < 0.4) continue;
      //   if(TMath::Abs(ctheta1) > 0.8 && TMath::Abs(ctheta2 > 0.8)) continue;
      //	    cout<<"*****openang"<<openang<<"********"<<ctheta1<<"****"<<ctheta2<<endl;
      //cout<<" rap "<<pairrap<<"invmass     "<<invmass<<"    pairpt    "<<pairpt<<endl;
      
      fHistMassPtPbarKPlus->Fill(invmass,pairpt);
      // Fill the ThnSparse     
      
      
    }// Kplus loop
    
  }//A Proton loop
  
  
}

//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ProcessMixed() {
  // Process mixed events
  
  Int_t iPro, iKminus, iAPro, iKplus;
  
  // Int_t nmixed = fResoBuffer->GetMixBuffSize();
  
  //  cout<<"nmixed*******"<<nmixed<<endl;
  
  // Loop over the event buffer
  for (UChar_t iMix = 1;iMix<fResoBuffer->GetMixBuffSize();iMix++){
    
    Int_t nproton = fResoBuffer->GetEvt(0)->GetNPro();
    Int_t nkmin = fResoBuffer->GetEvt(iMix)->GetNKMin();
    
    for (iPro = 0; iPro < nproton; iPro++){
      
      // Skip if unUseIt() entry
      if (!fResoBuffer->GetEvt(0)->fProTracks[iPro].UseIt())
  	continue;
      
      // Proton loop
    for (iKminus=0;iKminus< nkmin;iKminus++){
	
      // Skip if unUseIt() entry
      if (!(fResoBuffer->GetEvt(iMix))->fKMinTracks[iKminus].UseIt())
	continue;
      
      
	Double_t  pairrap =  Rapidity(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(iMix)->fKMinTracks[iKminus]);
	Double_t  invmass = MInv(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(iMix)->fKMinTracks[iKminus]);
	Double_t  pairpt  = Pt(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(iMix)->fKMinTracks[iKminus]);
	//	Double_t  openang  = OpeningAngle(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(iMix)->fKMinTracks[iKminus]);
	
	//cout<<invmass<<"****"<<pairpt<<"**********mixed"<<endl;
	
	if(TMath::Abs(pairrap) > 0.5) continue;
	
	
	//if(openang < 0.4) continue;
	
	fHistMassPtPKMinMix->Fill(invmass,pairpt);
  	
    }// Kmin loop
    
    }// Proton loop
    
    
    Int_t npbar   = fResoBuffer->GetEvt(0)->GetNAPro();
    
    Int_t nkplus  = fResoBuffer->GetEvt(iMix)->GetNKPlus();
    
    //    for (iAPro = 0; iAPro < fResoBuffer->GetEvt(0)->GetNAPro(); iAPro++){
    
    for (iAPro = 0; iAPro < npbar; iAPro++){
      // Skip if unUseIt() entry
      if (!fResoBuffer->GetEvt(0)->fAProTracks[iAPro].UseIt())
	continue;
      // Kplus loop
      for (iKplus=0;iKplus< nkplus ;iKplus++){
	
	// Skip if unUseIt() entry
	if (!fResoBuffer->GetEvt(iMix)->fKPlusTracks[iKplus].UseIt())
	  continue;
	
	Double_t  pairrap =  Rapidity(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(iMix)->fKPlusTracks[iKplus]);
	Double_t  invmass = MInv(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(iMix)->fKPlusTracks[iKplus]);
	Double_t  pairpt  = Pt(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(iMix)->fKPlusTracks[iKplus]);
	
	// 	Double_t  openang  = OpeningAngle(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(iMix)->fKPlusTracks[iKplus]);
	
    
	if(TMath::Abs(pairrap) > 0.5) continue;
	    	    
	//if(openang < 0.4) continue;

       //cout<<invmass<<"****"<<pairpt<<"**********mixed"<<endl;
	fHistMassPtPbarKPlusMix->Fill(invmass,pairpt);
          	
    }// AProton loop
    }// kplus loop
    
  }// Event buffer loop
  
}// End of void ProcessMixed 
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ProcessLikeSignBkg() {
  
  Int_t iPro, iKminus, iAPro, iKplus ;
  
  // Proton K+ loop
  Int_t nproton = fResoBuffer->GetEvt(0)->GetNPro();
  Int_t nkplus = fResoBuffer->GetEvt(0)->GetNKPlus();
  
  
  for (iPro = 0; iPro < nproton; iPro++){
    // Skip if unUseIt() entry
    if (!fResoBuffer->GetEvt(0)->fProTracks[iPro].UseIt())
      continue;
    // Kplus loop
    for (iKplus=0;iKplus< nkplus;iKplus++){

      // Skip if unUseIt() entry
      if (!fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus].UseIt())
  	continue;

      Double_t  pairrap =  Rapidity(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus]);
      Double_t  invmass = MInv(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus]);
      Double_t  pairpt  = Pt(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus]);
      //  Double_t  openang  = OpeningAngle(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus]);
      //Double_t  ctheta1  = Costheta(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus]);
      //Double_t  ctheta2  = Costheta1(fResoBuffer->GetEvt(0)->fProTracks[iPro], fResoBuffer->GetEvt(0)->fKPlusTracks[iKplus]);

      if(TMath::Abs(pairrap) > 0.5) continue;

      // if(openang < 0.4) continue;
      //cout<<"rap"<<pairrap<<"invmass"<<invmass<<"pairpt"<<pairpt<<endl;

      fHistMassPtPKPlusLS->Fill(invmass,pairpt);
            
    }// Kaon loop

  }// Proton loop

  Int_t npbar = fResoBuffer->GetEvt(0)->GetNAPro();
  Int_t nkmin = fResoBuffer->GetEvt(0)->GetNKMin();

  for (iAPro = 0; iAPro < npbar; iAPro++){
    // Skip if unUseIt() entry
    if (!fResoBuffer->GetEvt(0)->fAProTracks[iAPro].UseIt())
      continue;
    // Kminus loop
    for (iKminus=0;iKminus<nkmin;iKminus++){

      // Skip if unUseIt() entry
      if (!fResoBuffer->GetEvt(0)->fKMinTracks[iKminus].UseIt())
  	continue;
      
       Double_t  pairrap =  Rapidity(fResoBuffer->GetEvt(0)->fAProTracks[iPro], fResoBuffer->GetEvt(0)->fKMinTracks[iKminus]);
       Double_t  invmass = MInv(fResoBuffer->GetEvt(0)->fAProTracks[iPro], fResoBuffer->GetEvt(0)->fKMinTracks[iKminus]);
       Double_t  pairpt  = Pt(fResoBuffer->GetEvt(0)->fAProTracks[iPro], fResoBuffer->GetEvt(0)->fKMinTracks[iKminus]);
       //  Double_t  openang  = OpeningAngle(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(0)->fKMinTracks[iKminus]);
       //Double_t  ctheta1  = Costheta(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(0)->fKMinTracks[iKminus]);
       //Double_t  ctheta2  = Costheta1(fResoBuffer->GetEvt(0)->fAProTracks[iAPro], fResoBuffer->GetEvt(0)->fKMinTracks[iKminus]);

       
       
       if(TMath::Abs(pairrap) > 0.5) continue;
       
       //       if(openang < 0.4) continue;
       //cout<<"rap"<<pairrap<<"invmass"<<invmass<<"pairpt"<<pairpt<<endl;

       fHistMassPtPbarKMinLS->Fill(invmass,pairpt);
       
       // Fill the ThnSparse     
      
                 
    }// Kaon loop

  }// Proton loop

   
    
    
} // 



Float_t AliAnalysisTaskLambdaStar::MInv(ResoBufferTrack track1, ResoBufferTrack track2){
  

  Float_t e1 = TMath::Sqrt(track1.fP[0]*track1.fP[0] + track1.fP[1]*track1.fP[1] + track1.fP[2]*track1.fP[2] + fkProMass*fkProMass);
  Float_t e2 = TMath::Sqrt(track2.fP[0]*track2.fP[0] + track2.fP[1]*track2.fP[1] + track2.fP[2]*track2.fP[2] + fkKaonMass*fkKaonMass);

  Double_t massinv =  TMath::Sqrt(((e1+e2) * (e1+e2)) - (((track1.fP[0]+track2.fP[0])*(track1.fP[0]+track2.fP[0])) + ((track1.fP[1]+track2.fP[1])*(track1.fP[1]+track2.fP[1])) + ((track1.fP[2]+track2.fP[2])*(track1.fP[2]+track2.fP[2]))));

 return massinv;
}

Float_t AliAnalysisTaskLambdaStar::Costheta(ResoBufferTrack track1, ResoBufferTrack track2){
  

  Double_t massvtx = 1.520 ;
  Double_t massp[2];
  massp[0] = fkProMass;
  massp[1] = fkKaonMass;

  Double_t pStar = TMath::Sqrt((massvtx*massvtx-(massp[0]*massp[0])-(massp[1]*massp[1]))*(massvtx*massvtx-(massp[0]*massp[0])-(massp[1]*massp[1]))- (4.* massp[0]*massp[0]*massp[1]*massp[1]))/(2.*massvtx);

  Double_t ener =  TMath::Sqrt((massvtx * massvtx) + (((track1.fP[0]+track2.fP[0])*(track1.fP[0]+track2.fP[0])) + ((track1.fP[1]+track2.fP[1])*(track1.fP[1]+track2.fP[1])) + ((track1.fP[2]+track2.fP[2])*(track1.fP[2]+track2.fP[2]))));


  Double_t momlam = TMath::Sqrt(((track1.fP[0]+track2.fP[0])*(track1.fP[0]+track2.fP[0])) + ((track1.fP[1]+track2.fP[1])*(track1.fP[1]+track2.fP[1])) + ((track1.fP[2]+track2.fP[2])*(track1.fP[2]+track2.fP[2])));

  //  Double_t e=E(pdgvtx);

  Double_t beta = momlam/ener;
  Double_t gamma = ener/massvtx;
  TVector3 mom(track1.fP[0],track1.fP[1],track1.fP[2]);
  TVector3 momTot((track1.fP[0] + track2.fP[0]), (track1.fP[1] + track2.fP[1]), (track1.fP[2] + track2.fP[2]));
  Double_t QlProng =  mom.Dot(momTot)/momTot.Mag();
  Double_t cts = (QlProng/gamma-beta*TMath::Sqrt(pStar*pStar+massp[0]*massp[0]))/pStar;

   return  TMath::Cos(cts);  

}

Float_t AliAnalysisTaskLambdaStar::Costheta1(ResoBufferTrack track1, ResoBufferTrack track2){
  

  Double_t massvtx = 1.520 ;
  Double_t massp[2];
  massp[0] = fkProMass;
  massp[1] = fkKaonMass;

  Double_t pStar = TMath::Sqrt((massvtx*massvtx-(massp[0]*massp[0])-(massp[1]*massp[1]))*(massvtx*massvtx-(massp[0]*massp[0])-(massp[1]*massp[1]))- (4.* massp[0]*massp[0]*massp[1]*massp[1]))/(2.*massvtx);

  Double_t ener =  TMath::Sqrt((massvtx * massvtx) + (((track1.fP[0]+track2.fP[0])*(track1.fP[0]+track2.fP[0])) + ((track1.fP[1]+track2.fP[1])*(track1.fP[1]+track2.fP[1])) + ((track1.fP[2]+track2.fP[2])*(track1.fP[2]+track2.fP[2]))));


  Double_t momlam = TMath::Sqrt(((track1.fP[0]+track2.fP[0])*(track1.fP[0]+track2.fP[0])) + ((track1.fP[1]+track2.fP[1])*(track1.fP[1]+track2.fP[1])) + ((track1.fP[2]+track2.fP[2])*(track1.fP[2]+track2.fP[2])));

  //  Double_t e=E(pdgvtx);

  Double_t beta = momlam/ener;
  Double_t gamma = ener/massvtx;
  TVector3 mom(track2.fP[0],track2.fP[1],track2.fP[2]);
  TVector3 momTot((track1.fP[0] + track2.fP[0]), (track1.fP[1] + track2.fP[1]), (track1.fP[2] + track2.fP[2]));
  Double_t QlProng =  mom.Dot(momTot)/momTot.Mag();
  Double_t cts = (QlProng/gamma-beta*TMath::Sqrt(pStar*pStar+massp[1]*massp[1]))/pStar;

   return  TMath::Cos(cts);  

}




//________________________________________________________________________
Float_t AliAnalysisTaskLambdaStar::Rapidity(ResoBufferTrack track1, ResoBufferTrack track2){

  Float_t e1 = TMath::Sqrt((track1.fP[0]*track1.fP[0]) + (track1.fP[1]*track1.fP[1]) + (track1.fP[2]*track1.fP[2]) + (fkProMass*fkProMass));
  Float_t e2 = TMath::Sqrt((track2.fP[0]*track2.fP[0]) + (track2.fP[1]*track2.fP[1]) + (track2.fP[2]*track2.fP[2]) + (fkKaonMass*fkKaonMass));
  Double_t e = e1 + e2;
  Double_t pz = track1.fP[2]+ track2.fP[2];

  if (e != TMath::Abs(pz)) { // energy was not equal to pz
    return 0.5*TMath::Log((e+pz)/(e-pz));
  } else { // energy was equal to pz
    return -999.;
  }

}


Double_t AliAnalysisTaskLambdaStar::OpeningAngle(ResoBufferTrack track1, ResoBufferTrack track2){


  Double_t e1e2 = (track1.fP[0]*track2.fP[0]) + (track1.fP[1]*track2.fP[1]) + (track1.fP[2]*track2.fP[2]);
  Double_t e2 = TMath::Sqrt((track2.fP[0]*track2.fP[0]) + (track2.fP[1]*track2.fP[1]) + (track2.fP[2]*track2.fP[2]));
  Double_t e1 = TMath::Sqrt((track1.fP[0]*track1.fP[0]) + (track1.fP[1]*track1.fP[1]) + (track1.fP[2]*track1.fP[2]));
  
  return TMath::Cos(e1e2/(e1*e2));
  //  return   57.296*(TMath::ACos(e1e2/e1*e2));

}



//________________________________________________________________________
Bool_t AliAnalysisTaskLambdaStar::goodDCA(AliAODTrack *track) {
  
  // Get the DCAxy and DCAz. There also exists a TPC only 
  // impact parameter, but this has not enough resolution 
  // to discriminate between primaries, secondaries and material
 
  Float_t xy=0.,rap=RapidityProton(track),pt=track->Pt();
  
      xy = DCAxy(track, fAOD);
    
  // Fill the DCAxy histograms
  if (track->Charge() > 0){
    fPriHistDCAxyYPtPro->Fill(xy,rap,pt);
  }
  else{
    fPriHistDCAxyYPtAPro->Fill(xy,rap,pt);
  }
  // Do a cut. 0.1 cm shows highest significance for primaries
  if (xy>fDCAxy)
    return kFALSE;
  return kTRUE;
}


Bool_t AliAnalysisTaskLambdaStar::goodDCAKaon(AliAODTrack *track) {
  
  // Get the DCAxy and DCAz. There also exists a TPC only 
  // impact parameter, but this has not enough resolution 
  // to discriminate between primaries, secondaries and material
 
  Float_t xy=0.,rap=RapidityKaon(track),pt=track->Pt();
  
	xy = DCAxy(track, fAOD);
      
   
  // Fill the DCAxy histograms
  if (track->Charge() > 0){
    fPriHistDCAxyYPtKPlus->Fill(xy,rap,pt);
  }
  else{
    fPriHistDCAxyYPtKMinus->Fill(xy,rap,pt);
  }
  // Do a cut. 0.1 cm shows highest significance for primaries
  if (xy>fDCAxy)
    return kFALSE;
  return kTRUE;
}
//_______________________________________________________________
Float_t AliAnalysisTaskLambdaStar::RapidityProton(AliAODTrack *track){
  // Can't find how to set the assumed mass for the AliAODTrack.
  // Same stuff as in AliAODTrack::Y() just with proton mass
  Double_t e = TMath::Sqrt(track->P()*track->P() + fkProMass*fkProMass);
  Double_t pz = track->Pz();
  if (e != TMath::Abs(pz)) { // energy was not equal to pz
    return 0.5*TMath::Log((e+pz)/(e-pz));
  } else { // energy was equal to pz
    return -999.;
  }
}


Float_t AliAnalysisTaskLambdaStar::RapidityKaon(AliAODTrack *track){
  // Can't find how to set the assumed mass for the AliAODTrack.
  // Same stuff as in AliAODTrack::Y() just with proton mass
  Double_t e = TMath::Sqrt(track->P()*track->P() + fkKaonMass*fkKaonMass);
  Double_t pz = track->Pz();
  if (e != TMath::Abs(pz)) { // energy was not equal to pz
    return 0.5*TMath::Log((e+pz)/(e-pz));
  } else { // energy was equal to pz
    return -999.;
  }
}


//________________________________________________________________________

//________________________________________________________________________
Float_t AliAnalysisTaskLambdaStar::DCAxy(const AliAODTrack *track, const AliVEvent *evt){
  // Note that AliAODTrack::PropagateToDCA() changes the track. 
  // Don't know whether this is what one wants?
  if(!track){
    printf("Pointer to track is zero!\n");
    return -9999.;
  }

  // Create an external parameter from the AODtrack
  AliExternalTrackParam etp; etp.CopyFromVTrack(track);
  // Propagation through the beam pipe would need a correction 
  // for material, I guess.
  if(etp.GetX()>3.) {
    printf("This method can be used only for propagation inside the beam pipe\n");
    printf("  id: %d, filtermap: %d\n",track->GetID(),track->GetFilterMap());
    return -9999.; 
  }
  // Do the propagation
  Double_t dca[2]={-9999.,-9999.},covar[3]={0.,0.,0.};
  if(!etp.PropagateToDCA(evt->GetPrimaryVertex(),evt->GetMagneticField(),10.,dca,covar)) return -9999.;
  // return the DCAxy
  return dca[0];
}
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::FillDedxHist(const AliVTrack *track){
  // This is for visualization. Fill the the dE/dx histograms
  // for all tracks, not only for those, where only the TPC
  // is used for PID. Thus avoiding the sharp cut off at a 
  // momentum of 0.75 GeV/c.
  
  // Positive tracks

    fPriHistTPCsignal->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
  
  

}
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::StoreGlobalTrackReference(AliAODTrack *track){
  
  // Check that the id is positive
  if(track->GetID()<0){
    //    printf("Warning: track has negative ID: %d\n",track->GetID());
    return;
  }

  // Check id is not too big for buffer
  if(track->GetID()>=fTrackBuffSize){
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n"
	   ,track->GetID(),fTrackBuffSize);
    return;
  }

  // Warn if we overwrite a track
  if(fGTI[track->GetID()]){
    // Seems like there are FilterMap 0 tracks
    // that have zero TPCNcls, don't store these!
    if( (!track->GetFilterMap()) &&
	(!track->GetTPCNcls())   )
      return;

    // Imagine the other way around,
    // the zero map zero clusters track
    // is stored and the good one wants 
    // to be added. We ommit the warning
    // and just overwrite the 'bad' track
    if( fGTI[track->GetID()]->GetFilterMap() ||
	fGTI[track->GetID()]->GetTPCNcls()   ){
      // If we come here, there's a problem
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
	     (fGTI[track->GetID()])->GetTPCNcls(),track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
	     (fGTI[track->GetID()])->GetFilterMap(),track->GetFilterMap());
    }
  } // Two tracks same id

 
  // Assign the pointer
  (fGTI[track->GetID()]) = track;
}
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ResetGlobalTrackReference(){
  // Sets all the pointers to zero. To be called at
  // the beginning or end of an event
  for(UShort_t i=0;i<fTrackBuffSize;i++){
    fGTI[i]=0;
  }
}
//________________________________________________________________________
Bool_t AliAnalysisTaskLambdaStar::AcceptTrack(const AliAODTrack *track){
  // Apply additional track cuts

  
  Float_t nCrossed = track->GetTPCClusterInfo(2, 1);
  if(nCrossed<fClusterTPC)
    return kFALSE;
  //if(!track->GetTPCNclsF())
  //  return kFALSE; // Note that the AliESDtrackCuts would here return kTRUE
  //if((nCrossed/track->GetTPCNclsF()) < .8)
  //  return kFALSE;

  if(TMath::Abs(track->Eta()) > 0.8) return kFALSE;
  if(TMath::Abs(track->Pt()) < 0.2) return kFALSE;

  return kTRUE;

}
//________________________________________________________________________

//________________________________________________________________________
Bool_t AliAnalysisTaskLambdaStar::GoodTPCFitMapSharedMap(const AliAODTrack *track){
  // Rejects tracks with shared clusters after filling a control histogram
  // This overload is used for primaries

  // Get the shared maps
  const TBits sharedMap = track->GetTPCSharedMap();
  // Fill a control histogram
  // Reject shared clusters
  if((sharedMap.CountBits()) >= 1){
      // Fill a control histogram
      fPriHistShare->Fill(sharedMap.CountBits());

    // Bad track, has too many shared clusters!
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________

//________________________________________________________________________
Float_t AliAnalysisTaskLambdaStar::Pt(ResoBufferTrack track1, ResoBufferTrack track2){
 
  return  TMath::Sqrt((track1.fP[0] + track2.fP[0])*(track1.fP[0] + track2.fP[0])  + (track1.fP[1] + track2.fP[1])*(track1.fP[1] + track2.fP[1]));

}


//________________________________________________________________________
AliAnalysisTaskLambdaStar& AliAnalysisTaskLambdaStar::operator=(const AliAnalysisTaskLambdaStar& atpl)
{
  if(this!=&atpl){
  // One operation with the atpl to get rid of the warning unused parameter
  fPrimaryVtxPosition[0]=atpl.fPrimaryVtxPosition[0];
  printf("Assignment operator not implemented\n");
  }
  return *this;
}
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}
//________________________________________________________________________
//
//
//     Classes in the class AliAnalysisTaskLambdaStar
//         ResoBuffer, ResoBufferEvent, ResoBufferV0 and ResoBufferTrack
//
//________________________________________________________________________
//
//                        ResoBufferTrack
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBufferTrack::ResoBufferTrack():
  fID(65535)

{
  // Standard constructor, initialize everything with values indicating 
  // a track that should not be used
  
  // No idea how to initialize the arrays nicely like the fID(65535)..
  for (UChar_t i=0;i<3;i++){
    fP[i]=-9999.;
  
    for (UChar_t j=0;j<9;j++){
      //      fXglobal[j][i]=-9999.;
      fXshifted[j][i]=-9999.;
    }
  }
}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBufferTrack::ResoBufferTrack(const AliAODTrack *track,const Float_t bfield,const Float_t priVtx[3]):
  fID(65535)  
{
 
  // Use the function to have the code in one place
  Set(track,bfield,priVtx);
}
//________________________________________________________________________

//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ResoBufferTrack::GetShiftedPositionAtShiftedRadii(const AliAODTrack *track, const Float_t bfield, const Float_t priVtx[3]){
  // Gets the global position of the track at nine different radii in the TPC
  // track is the track you want to propagate
  // bfield is the magnetic field of your event
  // globalPositionsAtRadii is the array of global positions in the radii and xyz
  
  // Initialize the array to something indicating there was no propagation

  for(Int_t i=0;i<9;i++){
    for(Int_t j=0;j<3;j++){
      fXshifted[i][j]=-9999.;
    }
  }

   // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp; etp.CopyFromVTrack(track);
  //  printf("\nAfter CopyFromVTrack\n");
  //  etp.Print();
 
  // The global position of the the track
  Double_t xyz[3]={-9999.,-9999.,-9999.};  

  // Counter for which radius we want
  Int_t iR=0; 
  // The radii at which we get the global positions
  // IROC (OROC) from 84.1 cm to 132.1 cm (134.6 cm to 246.6 cm)
  // Compare squared radii for faster code
  Float_t RSquaredWanted[9]={85.*85.,105.*105.,125.*125.,145.*145.,165.*165.,
			     185.*185.,205.*205.,225.*225.,245.*245.}; 
  // The shifted radius we are at, squared. Compare squared radii for faster code
  Float_t shiftedRadiusSquared=0;

  // Propagation is done in local x of the track
 
  for (Float_t x = 58.;x<247.;x+=1.){
  
  // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
    // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
    // the track is straight, i.e. has inifinite pt and doesn't get bent. 
    // If the track's momentum is smaller than infinite, it will develop a y-component,
    // which adds to the global radius

    // Stop if the propagation was not succesful. This can happen for low pt tracks
    // that don't reach outer radii
   
    if(!etp.PropagateTo(x,bfield))break;
    
    etp.GetXYZ(xyz); // GetXYZ returns global coordinates

    // Without shifting the primary vertex to (0.,0.,0.) the next line would just be
    // WRONG: globalRadiusSquared = xyz[0]*xyz[0]+xyz[1]*xyz[1];
    // but as we shift the primary vertex we want to compare positions at shifted radii.
    // I can't draw in ASCII but please take a piece of paper and just visualize it once.

    // Changing plus to minus on July10th2012

    shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
                         + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);

    // Roughly reached the radius we want
    if(shiftedRadiusSquared > RSquaredWanted[iR]){
      
      // Bigger loop has bad precision, we're nearly one centimeter too far, 
      // go back in small steps.
      while (shiftedRadiusSquared>RSquaredWanted[iR]){
	x-=.1;
	//	printf("propagating to x %5.2f\n",x);
	if(!etp.PropagateTo(x,bfield))break;
	etp.GetXYZ(xyz); // GetXYZ returns global coordinates
	// Added the shifting also here on July11th2012
	shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
	                     + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);
      }
      //      printf("At Radius:%05.2f (local x %5.2f). Setting position to x %4.1f y %4.1f z %4.1f\n",TMath::Sqrt(globalRadiusSquared),x,xyz[0],xyz[1],xyz[2]);
      fXshifted[iR][0]=xyz[0]-priVtx[0];
      fXshifted[iR][1]=xyz[1]-priVtx[1];
      fXshifted[iR][2]=xyz[2]-priVtx[2];
      // Indicate we want the next radius    
      iR+=1;
    }
    if(iR>=8){
      // TPC edge reached
      return;
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ResoBufferTrack::Set(const AliAODTrack *track,const Float_t bfield,const Double_t priVtx[3]){
  // Overloaded function
  Float_t priVtxPos[3]={(Float_t)priVtx[0],(Float_t)priVtx[1],(Float_t)priVtx[2]};
  Set(track,bfield,priVtxPos);
}


void AliAnalysisTaskLambdaStar::ResoBufferTrack::SetP(const AliAODTrack *track){
  fP[0] = track->Px();
  fP[1] = track->Py();
  fP[2] = track->Pz();

}
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ResoBufferTrack::Set(const AliAODTrack *track,const Float_t bfield,const Float_t priVtx[3]){
 
  // Set the ID, a good ID also indicates to use the track
  if(track->GetID() >=0){
    // global tracks, i.e. v0 daughters
    fID = track->GetID();
  }
  else {
    // e.g. tpc only tracks, i.e. primary protons
    fID = -track->GetID()-1;

  }
  // Set the momentum
  
   track->PxPyPz(fP);

  
  GetShiftedPositionAtShiftedRadii(track,bfield,priVtx);

}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBufferTrack::ResoBufferTrack(const ResoBufferTrack& fbt):
  fID(fbt.fID)
 {
  // Copy constructor

  for (UChar_t i=0;i<3;i++){
    fP[i]=fbt.fP[i];
    for (UChar_t j=0;j<9;j++){
      //      fXglobal[j][i]=fbt.fXglobal[j][i];
      fXshifted[j][i]=fbt.fXshifted[j][i];
    }
  }
 }
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBufferTrack& AliAnalysisTaskLambdaStar::ResoBufferTrack::operator=(const ResoBufferTrack& fbt){
  // Assignment operator, from wikipedia :)
  
  // Protect against self-assignment
  if(this != &fbt){
    fID = fbt.fID;
    for (UChar_t i=0;i<3;i++){
      fP[i]=fbt.fP[i];
      for (UChar_t j=0;j<9;j++){
	//	fXglobal[j][i]=fbt.fXglobal[j][i];
	fXshifted[j][i]=fbt.fXshifted[j][i];
      }
    }
  }
  // By convention, always return *this (Could it be the convention is called
  return *this;
}
//________________________________________________________________________






//                        ResoBufferEvent
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBufferEvent::ResoBufferEvent():
  fPriTrackLim(0)
  ,fProTracks(0),fAProTracks(0)
  ,fKPlusTracks(0),fKMinTracks(0)
  ,fNProTracks(0),fNAProTracks(0),fNKPlusTracks(0),fNKMinTracks(0)
  ,fBfield(-9999.)
{
  // Standard constructor, all pointer to zero
  fPriVtxPos[0]=-9999.;
  fPriVtxPos[1]=-9999.;
  fPriVtxPos[2]=-9999.;

  printf("This constructor has zero size in the arrays!\n");
}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBufferEvent::ResoBufferEvent(const UShort_t priTrackBuff,const Double_t bfield,const Double_t priVtxPos[3]):
  fPriTrackLim(priTrackBuff)
  ,fProTracks(new ResoBufferTrack[fPriTrackLim])
  ,fAProTracks(new ResoBufferTrack[fPriTrackLim])
  ,fKPlusTracks(new ResoBufferTrack[fPriTrackLim])
  ,fKMinTracks(new ResoBufferTrack[fPriTrackLim])
  ,fNProTracks(0),fNAProTracks(0),fNKPlusTracks(0),fNKMinTracks(0)
  ,fBfield(-bfield)

{
  // Constructor.
  fPriVtxPos[0] = priVtxPos[0]; // This is some old C++
  fPriVtxPos[1] = priVtxPos[1];
  fPriVtxPos[2] = priVtxPos[2];
}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBufferEvent::ResoBufferEvent(const UShort_t priTrackBuff):
  fPriTrackLim(priTrackBuff)
  ,fProTracks(new ResoBufferTrack[fPriTrackLim])
  ,fAProTracks(new ResoBufferTrack[fPriTrackLim])
  ,fKPlusTracks(new ResoBufferTrack[fPriTrackLim])
  ,fKMinTracks(new ResoBufferTrack[fPriTrackLim])
  ,fNProTracks(0),fNAProTracks(0),fNKPlusTracks(0),fNKMinTracks(0)
  ,fBfield(-9999.)
 
{  
  // Constructor. fBfield and fPriVtxPos not needed yet, can be set later.
  fPriVtxPos[0] = -9999.; // This is C++03
  fPriVtxPos[1] = -9999.;
  fPriVtxPos[2] = -9999.;

  //  printf("constructed eventwith NBgLam: %u\n",fNBgLamTracks);
}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBufferEvent::ResoBufferEvent(const ResoBufferEvent &fbe):
  fPriTrackLim(fbe.GetPriTrackLim())
  ,fProTracks(new ResoBufferTrack[fPriTrackLim])
  ,fAProTracks(new ResoBufferTrack[fPriTrackLim])
  ,fKPlusTracks(new ResoBufferTrack[fPriTrackLim])
  ,fKMinTracks(new ResoBufferTrack[fPriTrackLim])
  ,fNProTracks(fbe.GetNPro()),fNAProTracks(fbe.GetNAPro())
  ,fNKPlusTracks(fbe.GetNKPlus()),fNKMinTracks(fbe.GetNKMin())
  ,fBfield(fbe.GetBfield())
{
  // Copy constructor
  fbe.GetVtxPos(fPriVtxPos);
  // Avoid to much creation and deletion of objects
  UShort_t i;
  // Copy the primary tracks
  for (i=0;i<fPriTrackLim;i++){
    fProTracks[i]=fbe.fProTracks[i];
    fAProTracks[i]=fbe.fAProTracks[i];
    fKPlusTracks[i]=fbe.fKPlusTracks[i];
    fKMinTracks[i]=fbe.fKMinTracks[i];
  }

}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBufferEvent& AliAnalysisTaskLambdaStar::ResoBufferEvent::operator=(const ResoBufferEvent &fbe){
  // Assignment operator

  // Protect against self-assignment
  if(this!=&fbe){
    // Well, we use arrays of a constant size to avoid
    // excessive memory allocation and won't give this up.
    // So we'll only copy as much as fits on the left side
    // from the right side.
    // DON'T COPY THE ARRAY SIZES fV0Lim AND fPriTrackLim !!!
    if(fPriTrackLim < fbe.GetPriTrackLim() ){
      // AliWarning(Form("Trying to assign too big event (buffer %d/%d) to"
      // 		    " this (buffer %d/%d). Only partially copying.",
      // 		    fbe.GetPriTrackLim(),
      // 		    fPriTrackLim,));
      printf("Trying to assign too big event (buffer %d) to"
    		    " this (buffer %d. Only partially copying.\n",
	     fbe.GetPriTrackLim(),
	     fPriTrackLim);
    }
    // Always start with the easy stuff :)
    fbe.GetVtxPos(fPriVtxPos);
    fBfield = fbe.GetBfield();
    // Number of tracks is minimum of array size of 'this'
    // and the number of tracks from the right side
    fNProTracks = TMath::Min(fPriTrackLim,fbe.GetNPro());
    fNAProTracks = TMath::Min(fPriTrackLim,fbe.GetNAPro());

    fNKPlusTracks = TMath::Min(fPriTrackLim,fbe.GetNKPlus());
    fNKMinTracks = TMath::Min(fPriTrackLim,fbe.GetNKMin());
    
    // Avoid creation and deletion of 'i' for every loop
    UShort_t i;
    // Copy primary tracks. No need to set a 'bad track'
    // flag for the entries above GetNPro() (...) as
    // above everything is bad by definition.
    // Protons
    for (i=0;i<GetNPro();i++)
      fProTracks[i]=fbe.fProTracks[i];
    // Anti-protons
    for (i=0;i<GetNAPro();i++)
      fAProTracks[i]=fbe.fAProTracks[i];

    // Kplus
    for (i=0;i<GetNPro();i++){
      fKPlusTracks[i]=fbe.fKPlusTracks[i];
    }
    // Anti-lambdas
    for (i=0;i<GetNKMin();i++){
      fKMinTracks[i]=fbe.fKMinTracks[i];
    }
    
  }
  return *this;
}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBufferEvent::~ResoBufferEvent(){
  // Destructor

  // Delete the arrays of tracks,
  // note the [] with the delete
  if(fProTracks){
    delete[] fProTracks;
    fProTracks=0;
  }
  if(fAProTracks){
    delete[] fAProTracks;
    fAProTracks=0;
  }
  if(fKPlusTracks){
    delete[] fKPlusTracks;
    fKPlusTracks=0;
  }
  if(fKMinTracks){
    delete[] fKMinTracks;
    fKMinTracks=0;
  }

}

//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ResoBufferEvent::Reset(const Double_t bfield, const Double_t priVtxPos[3]){
 
 // Reset the old event, i.e., make clear 'here is no info'
  // by setting the 'number of stored ...' to zero
  fNProTracks=0;
  fNAProTracks=0;
  fNKPlusTracks=0;
  fNKMinTracks=0;
  
  // And set the new event properties 
  fBfield = bfield;
  fPriVtxPos[0]=priVtxPos[0];
  fPriVtxPos[1]=priVtxPos[1];
  fPriVtxPos[2]=priVtxPos[2];
}
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ResoBufferEvent::AddPro(const AliAODTrack *track){
  // Add a proton to this event

  // Check whether there is still space in the array
  if(fNProTracks > fPriTrackLim-1){
    // AliWarning(Form("Cannot add proton, array size (%d) too small"
    // 		    ,fPriTrackLim));
    printf("Cannot add proton, array size (%d) too small\n"
    		    ,fPriTrackLim);
    return;
  }

  fProTracks[fNProTracks].Set(track,fBfield,fPriVtxPos);
  //  fProTracks[fNProTracks].SetP(track);

  fNProTracks++;
  //  printf("Added proton %d/%d\n",fNProTracks,fPriTrackLim);

}  
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ResoBufferEvent::AddAPro(const AliAODTrack *track){
  // Add a anti-proton to this event

  
  if(fNAProTracks > fPriTrackLim-1){
  
    printf("Cannot add anti-proton, array size (%d) too small\n"
		    ,fPriTrackLim);
    return;
  }
  // Add the V0 at the end of the array

  fAProTracks[fNAProTracks].Set(track,fBfield,fPriVtxPos);

  //fAProTracks[fNAProTracks].SetP(track);

  fNAProTracks++;
}  
//________________________________________________________________________

void AliAnalysisTaskLambdaStar::ResoBufferEvent::AddKPlus(const AliAODTrack *track){
  // Add a kplus to this event

  // Check whether there is still space in the array
  if(fNKPlusTracks > fPriTrackLim-1){

    printf("Cannot add kplus, array size (%d) too small\n"
    		    ,fPriTrackLim);
    return;
  }

  fKPlusTracks[fNKPlusTracks].Set(track,fBfield,fPriVtxPos);
  //fKPlusTracks[fNKPlusTracks].SetP(track);
  fNKPlusTracks++;


}  
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ResoBufferEvent::AddKMin(const AliAODTrack *track){
  // Add a anti-proton to this event

  
  if(fNKMinTracks > fPriTrackLim-1){
  
    printf("Cannot add kminus, array size (%d) too small\n"
		    ,fPriTrackLim);
    return;
  }
  // Add the V0 at the end of the array
  fKMinTracks[fNKMinTracks].Set(track,fBfield,fPriVtxPos);
  //fKMinTracks[fNKMinTracks].SetP(track);
  fNKMinTracks++;
}  

//________________________________________________________________________


//________________________________________________________________________
//
//                        ResoBuffer
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBuffer::ResoBuffer() :
  fkZvertexBins(0),
  fkCentBins(0),
  fkMixBuffSize(0),
  fkPriTrackLim(0),
  fZvertexAxis(0),
  fCentAxis(0),
  fCurEvt(0),
  fEC(0)
{
  // Dummy constructor, create arrays with zero size
  // Note that some data member are constant, you
  // won't be able to create the ResoBuffer first with this
  // constructor and then set the appropiate size.
  
}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBuffer::ResoBuffer(const UChar_t ZvertexBins,const UChar_t CentBins,const UChar_t MixBuff,const UShort_t PriTrackLim, const Float_t AbsZvertexCut,const Float_t CentCut) :
  fkZvertexBins(ZvertexBins),
  fkCentBins(CentBins),
  fkMixBuffSize(MixBuff),
  fkPriTrackLim(PriTrackLim),
  fZvertexAxis(new TAxis(fkZvertexBins,-AbsZvertexCut,AbsZvertexCut)),
  fCentAxis(new TAxis (fkCentBins,0.0,CentCut)),
  fCurEvt(new ResoBufferEvent *[fkMixBuffSize]),
  fEC(new ResoBufferEvent ***[fkZvertexBins])
{
  // Constructor, creates at once all events with all tracks
  //  printf ("Creating with pritracklim %d and v0lim %d\n",fkPriTrackLim,fkV0Lim);

  // Create the array step by step
  // Bins in z of the primary vertex position. Do this as
  // the detector looks different from a different z coordinate
  for (UChar_t iZBin=0;iZBin<fkZvertexBins;iZBin++){
    fEC[iZBin] = new ResoBufferEvent **[fkCentBins];
    // Bins in centrality
    for (UChar_t iCentBin=0;iCentBin<fkCentBins;iCentBin++){
      fEC[iZBin][iCentBin] = new ResoBufferEvent *[fkMixBuffSize];
      // The number of events to keep for one mixing class
      for(UChar_t iMixBuff=0;iMixBuff<fkMixBuffSize;iMixBuff++){
	// Create an event to hold the info for mixing
	fEC[iZBin][iCentBin][iMixBuff] = new ResoBufferEvent(fkPriTrackLim);
      }
    }
  }
}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBuffer::ResoBuffer(const AliAnalysisTaskLambdaStar::ResoBuffer &fb) :
  fkZvertexBins(fb.fkZvertexBins),
  fkCentBins(fb.fkCentBins),
  fkMixBuffSize(fb.fkMixBuffSize),
  fkPriTrackLim(fb.fkPriTrackLim),
  fZvertexAxis(new TAxis(*(fb.fZvertexAxis))),
  fCentAxis(new TAxis (*(fb.fCentAxis))),
  fCurEvt(new ResoBufferEvent *[fkMixBuffSize]),
  fEC(new ResoBufferEvent ***[fkZvertexBins])
{
  // Copy constructor. Linux complains not having this and 
  // compiling this task with aliroot

  printf("ResoBuffer ctor not tested yet, be cautious\n");
  
  // Create the array step by step
  // Bins in z of the primary vertex position. Do this as
  // the detector looks different from a different z coordinate
  for (UChar_t iZBin=0;iZBin<fkZvertexBins;iZBin++){
    fEC[iZBin] = new ResoBufferEvent **[fkCentBins];
    // Bins in centrality
    for (UChar_t iCentBin=0;iCentBin<fkCentBins;iCentBin++){
      fEC[iZBin][iCentBin] = new ResoBufferEvent *[fkMixBuffSize];
      // The number of events to keep for one mixing class
      for(UChar_t iMixBuff=0;iMixBuff<fkMixBuffSize;iMixBuff++){
	// Create an event to hold the info for mixing
	fEC[iZBin][iCentBin][iMixBuff] = new ResoBufferEvent(*(fb.fEC[iZBin][iCentBin][iMixBuff]));
      }
    }
  }
}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBuffer& AliAnalysisTaskLambdaStar::ResoBuffer::operator=(const AliAnalysisTaskLambdaStar::ResoBuffer& fb){
  //Assignment operator
  if(this!=&fb){
    printf("ResoBuffer assignment operator not implemented\n");
  }
  return *this;
  
}
//________________________________________________________________________
AliAnalysisTaskLambdaStar::ResoBuffer::~ResoBuffer(){
  // Destructor
  // The axes to fin the correct bins
  if(fZvertexAxis){
    delete fZvertexAxis;
    fZvertexAxis=0;
  }
  if(fCentAxis){
    delete fCentAxis;
    fCentAxis=0;
  }
  // fCurEvt is an array of pointer
  if(fCurEvt){
    delete[] fCurEvt;
    fCurEvt=0;
  }
  // Delete all the events and the pointer to them
  for (UChar_t iZBin=0;iZBin<fkZvertexBins;iZBin++){
    for (UChar_t iCentBin=0;iCentBin<fkCentBins;iCentBin++){
      for(UChar_t iMixBuff=0;iMixBuff<fkMixBuffSize;iMixBuff++){
	if(fEC[iZBin][iCentBin][iMixBuff]){
	  delete fEC[iZBin][iCentBin][iMixBuff];
	  fEC[iZBin][iCentBin][iMixBuff]=0;
	}
      }
      if(fEC[iZBin][iCentBin]){
	delete fEC[iZBin][iCentBin];
	fEC[iZBin][iCentBin]=0;
      }
    }
    if(fEC[iZBin]){
      delete fEC[iZBin];
      fEC[iZBin]=0;
    }
  }
  delete[] fEC;
  fEC=0;
}
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ResoBuffer::ShiftAndAdd(AliAODEvent *evt){
  // Shift the events in the appropiate centrality / zvertex bin and set the 
  // current event pointer correctly
  Double_t priVtxPos[3];
  evt->GetPrimaryVertex()->GetXYZ(priVtxPos);
  //  printf("Mag field: %f\n",evt->GetMagneticField());
  ShiftAndAdd(evt->GetMagneticField(),
	      priVtxPos,
	      evt->GetCentrality()->GetCentralityPercentileUnchecked("V0M"));
}
//________________________________________________________________________
void AliAnalysisTaskLambdaStar::ResoBuffer::ShiftAndAdd(const Double_t bfield,const Double_t priVtxPos[3],const Float_t centrality){

  // Shift the events in the appropiate centrality / zvertex bin and set the 
  // current event pointer correctly

  // Find the correct centrality/zvertex bin 

  const UChar_t ZvertexBin = fZvertexAxis->FindFixBin(priVtxPos[2]) - 1; // -1 for array starting at 0

  const UChar_t CentBin = fCentAxis->FindFixBin(centrality) - 1;// -1 for array starting at 0

  // The new current event is the old last event
  fCurEvt[0] = fEC[ZvertexBin][CentBin][fkMixBuffSize-1];

  // Shift the pointer, starting from the back
  UChar_t iMix;
  for(iMix=fkMixBuffSize-1;iMix>0;iMix--){
    fEC[ZvertexBin][CentBin][iMix] = fEC[ZvertexBin][CentBin][iMix-1];
  }
  // And reset the zero'th one
  fEC[ZvertexBin][CentBin][0] = fCurEvt[0];
  fEC[ZvertexBin][CentBin][0]->Reset(bfield,priVtxPos);
  // Also set the pointer to the other events..
  for (iMix=1;iMix<fkMixBuffSize;iMix++){
    fCurEvt[iMix] = fEC[ZvertexBin][CentBin][iMix];
  }
}
