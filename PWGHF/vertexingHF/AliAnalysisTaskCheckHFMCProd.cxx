#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliGenerator.h"
#include "AliVertexingHFUtils.h"
#include "AliMultiplicity.h"
#include <TParticle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TChain.h>
#include "AliESDInputHandlerRP.h"
#include "AliAODInputHandler.h"
#include "AliAODTracklets.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskCheckHFMCProd.h"

/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */ 

//*************************************************************************
// Implementation of class AliAnalysisTaskCheckHFMCProd
// AliAnalysisTask to check MC production at ESD+Kine level
// 
//
// Authors: F. Prino, prino@to.infn.it
//          
//*************************************************************************

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCheckHFMCProd);
/// \endcond

//______________________________________________________________________________
AliAnalysisTaskCheckHFMCProd::AliAnalysisTaskCheckHFMCProd() : 
  AliAnalysisTaskSE("HFMCChecks"), 
  fOutput(0),
  fHistoNEvents(0),
  fHistoPhysPrim(0),
  fHistoPhysPrimPiKPi09(0),
  fHistoPhysPrimPiKPi09vsb(0),
  fHistoTracks(0),
  fHistoSelTracks(0),
  fHistoTracklets(0),
  fHistoTrackletsEta1(0),
  fHistoPtPhysPrim(0),
  fHistoEtaPhysPrim(0),
  fHistoSPD3DVtxX(0),
  fHistoSPD3DVtxY(0),
  fHistoSPD3DVtxZ(0),
  fHistoSPDZVtxX(0),
  fHistoSPDZVtxY(0),
  fHistoSPDZVtxZ(0),
  fHistoTRKVtxX(0),
  fHistoTRKVtxY(0),
  fHistoTRKVtxZ(0),
  fHistoNcharmed(0),
  fHistoNbVsNc(0),
  fHistOriginPrompt(0),
  fHistOriginFeeddown(0),
  fHistPtBDecLenBXYFeeddown(0),
  fHistMotherID(0),
  fHistDSpecies(0),
  fHistBSpecies(0),
  fHistLcDecayChan(0),
  fHistNcollHFtype(0),
  fHistNinjectedvsb(0),
  fHistEtaPhiPtGenEle(0),
  fHistEtaPhiPtGenPi(0),
  fHistEtaPhiPtGenK(0),
  fHistEtaPhiPtGenPro(0),
  fHistEtaPhiPtRecEle(0),
  fHistEtaPhiPtRecPi(0),
  fHistEtaPhiPtRecK(0),
  fHistEtaPhiPtRecPro(0),
  fHistPtRecVsPtGen(0),
  fHistPhiRecVsPhiGen(0),
  fHistEtaRecVsEtaGen(0),
  fHistPtRecGood(0),
  fHistPtRecFake(0),
  fSearchUpToQuark(kFALSE),
  fSystem(0),
  fESDtrackCuts(0x0),
  fReadMC(kTRUE),
  fPtMin(0.),
  fPtMax(40.),
  fNPtBins(40),
  fPtMinB(0.),
  fPtMaxB(40.),
  fNPtBinsB(40),
  fYMin(-2.),
  fYMax(2.),
  fNYBins(40),
  fEvent(nullptr)
{
  //
  for(Int_t i=0; i<5; i++){
    fHistBYPtAllDecay[i]=0x0;
    fHistYPtAllDecay[i]=0x0;
    fHistYPtPromptAllDecay[i]=0x0;
    fHistYPtFeeddownAllDecay[i]=0x0;
    fHistYPtPrompt[i]=0x0;
    fHistYPtFeeddown[i]=0x0;
    fHistPtDDecLenPrompt[i]=0x0;
    fHistPtDDecLenXYPrompt[i]=0x0;
    fHistPtDPtBDecLenFeeddown[i]=0x0;
    fHistPtDPtBDecLenXYFeeddown[i]=0x0;
  }
  for(Int_t i=0; i<2; i++){
    fHistYPtD0byDecChannel[i]=0x0;
    fHistYPtDplusbyDecChannel[i]=0x0;
    fHistYPtDsbyDecChannel[i]=0x0;
  }
  fHistYPtDplusbyDecChannel[2]=0x0;

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskCheckHFMCProd::~AliAnalysisTaskCheckHFMCProd(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  delete fESDtrackCuts;

}
   
//___________________________________________________________________________
void AliAnalysisTaskCheckHFMCProd::UserCreateOutputObjects() {
  /// create output histos

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistoNEvents = new TH1F("hNEvents", "Number of processed events",3,-0.5,2.5);
  fHistoNEvents->SetMinimum(0);
  fOutput->Add(fHistoNEvents);

  Double_t maxMult=100.;
  if(fSystem==1) maxMult=10000.;
  if(fSystem==2) maxMult=500.;
  fHistoPhysPrim = new TH1F("hPhysPrim","",100,-0.5,maxMult-0.5);
  fOutput->Add(fHistoPhysPrim);
  fHistoPhysPrimPiKPi09 = new TH1F("hPhysPrimPiKPi09","",100,-0.5,maxMult-0.5);
  fOutput->Add(fHistoPhysPrimPiKPi09);
  fHistoPhysPrimPiKPi09vsb = new TH2F("hPhysPrimPiKPi09vsb","",50,0.,20.,100,-0.5,maxMult-0.5);
  fOutput->Add(fHistoPhysPrimPiKPi09vsb);

  fHistoTracks = new TH1F("hTracks","",100,-0.5,maxMult*2-0.5);
  fOutput->Add(fHistoTracks);
  fHistoSelTracks = new TH1F("hSelTracks","",100,-0.5,maxMult-0.5);
  fOutput->Add(fHistoSelTracks);
  fHistoTracklets = new TH1F("hTracklets","",100,-0.5,maxMult-0.5);
  fOutput->Add(fHistoTracklets);
  fHistoTrackletsEta1 = new TH1F("hTrackletsEta1","",100,-0.5,maxMult-0.5);
  fOutput->Add(fHistoTrackletsEta1);
  fHistoPtPhysPrim = new TH1F("hPtPhysPrim","",100,0.,20.);
  fOutput->Add(fHistoPtPhysPrim);
  fHistoEtaPhysPrim = new TH1F("hEtaPhysPrim","",100,-10.,10.);
  fOutput->Add(fHistoEtaPhysPrim);

  fHistoSPD3DVtxX = new TH1F("hSPD3DvX","",100,-1.,1.);
  fOutput->Add(fHistoSPD3DVtxX);
  fHistoSPD3DVtxY = new TH1F("hSPD3DvY","",100,-1.,1.);
  fOutput->Add(fHistoSPD3DVtxY);
  fHistoSPD3DVtxZ = new TH1F("hSPD3DvZ","",100,-15.,15.);
  fOutput->Add(fHistoSPD3DVtxZ);

  fHistoSPDZVtxX = new TH1F("hSPDZvX","",100,-1.,1.);
  fOutput->Add(fHistoSPDZVtxX);
  fHistoSPDZVtxY = new TH1F("hSPDZvY","",100,-1.,1.);
  fOutput->Add(fHistoSPDZVtxY);
  fHistoSPDZVtxZ = new TH1F("hSPDZvZ","",100,-15.,15.);
  fOutput->Add(fHistoSPDZVtxZ);


  fHistoTRKVtxX = new TH1F("hTRKvX","",100,-1.,1.);
  fOutput->Add(fHistoTRKVtxX);
  fHistoTRKVtxY = new TH1F("hTRKvY","",100,-1.,1.);
  fOutput->Add(fHistoTRKVtxY);
  fHistoTRKVtxZ = new TH1F("hTRKvZ","",100,-15.,15.);
  fOutput->Add(fHistoTRKVtxZ);

  Int_t nBinscb=11;
  if(fSystem==1) nBinscb=200;
  if(fSystem==2) nBinscb=21;
  Double_t maxncn=nBinscb-0.5;
  fHistoNcharmed = new TH2F("hncharmed","",100,-0.5,maxMult-0.5,nBinscb,-0.5,maxncn);
  fOutput->Add(fHistoNcharmed);
  fHistoNbVsNc = new TH2F("hnbvsnc","",nBinscb,-0.5,maxncn,nBinscb,-0.5,maxncn);
  fOutput->Add(fHistoNbVsNc);

  fHistYPtPrompt[0] = new TH2F("hyptD0prompt","D0 - Prompt",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtPrompt[1] = new TH2F("hyptDplusprompt","Dplus - Prompt",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtPrompt[2] = new TH2F("hyptDstarprompt","Dstar - Prompt",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtPrompt[3] = new TH2F("hyptDsprompt","Ds - Prompt",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtPrompt[4] = new TH2F("hyptLcprompt","Lc - Prompt",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);

  fHistBYPtAllDecay[0] = new TH2F("hyptB0AllDecay","B0 - All",fNPtBinsB, fPtMinB, fPtMaxB, fNYBins, fYMin, fYMax);
  fHistBYPtAllDecay[1] = new TH2F("hyptBplusAllDecay","Bplus - All",fNPtBinsB, fPtMinB, fPtMaxB, fNYBins, fYMin, fYMax);
  fHistBYPtAllDecay[2] = new TH2F("hyptBstarAllDecay","Bstar - All",fNPtBinsB, fPtMinB, fPtMaxB, fNYBins, fYMin, fYMax);
  fHistBYPtAllDecay[3] = new TH2F("hyptBsAllDecay","Bs - All",fNPtBinsB, fPtMinB, fPtMaxB, fNYBins, fYMin, fYMax);
  fHistBYPtAllDecay[4] = new TH2F("hyptLbAllDecay","LB - All",fNPtBinsB, fPtMinB, fPtMaxB, fNYBins, fYMin, fYMax);

  fHistYPtAllDecay[0] = new TH2F("hyptD0AllDecay","D0 - All",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtAllDecay[1] = new TH2F("hyptDplusAllDecay","Dplus - All",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtAllDecay[2] = new TH2F("hyptDstarAllDecay","Dstar - All",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtAllDecay[3] = new TH2F("hyptDsAllDecay","Ds - All",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtAllDecay[4] = new TH2F("hyptLcAllDecay","Lc - All",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);

  fHistYPtPromptAllDecay[0] = new TH2F("hyptD0promptAllDecay","D0 - Prompt",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtPromptAllDecay[1] = new TH2F("hyptDpluspromptAllDecay","Dplus - Prompt",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtPromptAllDecay[2] = new TH2F("hyptDstarpromptAllDecay","Dstar - Prompt",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtPromptAllDecay[3] = new TH2F("hyptDspromptAllDecay","Ds - Prompt",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtPromptAllDecay[4] = new TH2F("hyptLcpromptAllDecay","Lc - Prompt",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);

  fHistYPtFeeddownAllDecay[0] = new TH2F("hyptD0feeddownAllDecay","D0 - FromB",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtFeeddownAllDecay[1] = new TH2F("hyptDplusfeeddownAllDecay","Dplus - FromB",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtFeeddownAllDecay[2] = new TH2F("hyptDstarfeeddownAllDecay","Dstar - FromB",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtFeeddownAllDecay[3] = new TH2F("hyptDsfeeddownAllDecay","Ds - FromB",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtFeeddownAllDecay[4] = new TH2F("hyptLcfeeddownAllDecay","Lc - FromB",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);


  fHistYPtFeeddown[0] = new TH2F("hyptD0feeddown","D0 - Feeddown",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtFeeddown[1] = new TH2F("hyptDplusfeeddown","Dplus - Feeddown",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtFeeddown[2] = new TH2F("hyptDstarfeedown","Dstar - Feeddown",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtFeeddown[3] = new TH2F("hyptDsfeedown","Ds - Feeddown",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtFeeddown[4] = new TH2F("hyptLcfeedown","Lc - Feeddown",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);

  TString pnams[5]={"D0","Dplus","Dstar","Ds","Lc"};
  for(Int_t j=0; j<5; j++){
    fHistPtDDecLenPrompt[j] = new TH2F(Form("h%sPtDDecLenPrompt",pnams[j].Data()),"; p_{T}(D) ; Dec Len (cm)",fNPtBins, fPtMin, fPtMax,500,0.,5.);
    fHistPtDDecLenXYPrompt[j] = new TH2F(Form("h%sPtDDecLenXYPrompt",pnams[j].Data()),"; p_{T}(D) ; Dec Len XY (cm)",fNPtBins, fPtMin, fPtMax,500,0.,5.);
    fHistPtDPtBDecLenFeeddown[j] = new TH3F(Form("h%sPtDPtBDecLenFeeddown",pnams[j].Data()),"; p_{T}(D) ; p_{T} (B) ; Dec Len (cm)",fNPtBins, fPtMin, fPtMax,100,0.,100.,500,0.,5.);
    fHistPtDPtBDecLenXYFeeddown[j] = new TH3F(Form("h%sPtDPtBDecLenXYFeeddown",pnams[j].Data()),"; p_{T}(D) ; p_{T} (B) ; Dec Len XY (cm)",fNPtBins, fPtMin, fPtMax,100,0.,100.,500,0.,5.);
  }
  
  for(Int_t ih=0; ih<5; ih++){
    fHistBYPtAllDecay[ih]->SetMinimum(0);
    fOutput->Add(fHistBYPtAllDecay[ih]);
    fHistYPtAllDecay[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtAllDecay[ih]);
    fHistYPtPromptAllDecay[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtPromptAllDecay[ih]);
    fHistYPtFeeddownAllDecay[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtFeeddownAllDecay[ih]);
    fHistYPtPrompt[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtPrompt[ih]);
    fHistYPtFeeddown[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtFeeddown[ih]);
    fOutput->Add(fHistPtDDecLenPrompt[ih]);
    fOutput->Add(fHistPtDDecLenXYPrompt[ih]);
    fOutput->Add(fHistPtDPtBDecLenFeeddown[ih]);
    fOutput->Add(fHistPtDPtBDecLenXYFeeddown[ih]);    
  }

  fHistYPtD0byDecChannel[0] = new TH2F("hyptD02","D0 - 2prong",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtD0byDecChannel[1] = new TH2F("hyptD04","D0 - 4prong",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtDplusbyDecChannel[0] = new TH2F("hyptDplusnonreson","Dplus - non reson",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtDplusbyDecChannel[1] = new TH2F("hyptDplusreson","Dplus - reson via K0*",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtDplusbyDecChannel[2] = new TH2F("hyptDplusKKpi","Dplus -> KKpi",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtDsbyDecChannel[0] = new TH2F("hyptDsphi","Ds - vis Phi",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);
  fHistYPtDsbyDecChannel[1] = new TH2F("hyptDsk0st","Ds - via k0*",fNPtBins, fPtMin, fPtMax, fNYBins, fYMin, fYMax);

  for(Int_t ih=0; ih<2; ih++){

    fHistYPtD0byDecChannel[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtD0byDecChannel[ih]);
    fHistYPtDplusbyDecChannel[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtDplusbyDecChannel[ih]);
    fHistYPtDsbyDecChannel[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtDsbyDecChannel[ih]);
  }
  fHistYPtDplusbyDecChannel[2]->SetMinimum(0);
  fOutput->Add(fHistYPtDplusbyDecChannel[2]);

  fHistOriginPrompt=new TH1F("hOriginPrompt","",100,0.,0.5);
  fHistOriginPrompt->SetMinimum(0);
  fOutput->Add(fHistOriginPrompt);
  fHistOriginFeeddown=new TH1F("hOriginFeeddown","",100,0.,0.5);
  fHistOriginFeeddown->SetMinimum(0);
  fOutput->Add(fHistOriginFeeddown);
  fHistPtBDecLenBXYFeeddown = new TH2F("hPtBDecLenBXYFeeddown","; p_{T} (B) ; Dec Len XY (cm)",100,0.,100.,500,0.,5.);
  fOutput->Add(fHistPtBDecLenBXYFeeddown);
  fHistMotherID=new TH1F("hMotherID","",1000,-1.5,998.5);
  fHistMotherID->SetMinimum(0);
  fOutput->Add(fHistMotherID);
  fHistDSpecies=new TH1F("hDSpecies","",10,-0.5,9.5);
  fHistDSpecies->GetXaxis()->SetBinLabel(1,"D0");
  fHistDSpecies->GetXaxis()->SetBinLabel(2,"D0bar");
  fHistDSpecies->GetXaxis()->SetBinLabel(3,"D+");
  fHistDSpecies->GetXaxis()->SetBinLabel(4,"D-");
  fHistDSpecies->GetXaxis()->SetBinLabel(5,"D*+");
  fHistDSpecies->GetXaxis()->SetBinLabel(6,"D*-");
  fHistDSpecies->GetXaxis()->SetBinLabel(7,"Ds+");
  fHistDSpecies->GetXaxis()->SetBinLabel(8,"Ds-");
  fHistDSpecies->GetXaxis()->SetBinLabel(9,"Lc+");
  fHistDSpecies->GetXaxis()->SetBinLabel(10,"Lc-");
  fHistDSpecies->SetMinimum(0);
  fOutput->Add(fHistDSpecies);
  fHistBSpecies=new TH1F("hBSpecies","",10,-0.5,9.5);
  fHistBSpecies->GetXaxis()->SetBinLabel(1,"B0");
  fHistBSpecies->GetXaxis()->SetBinLabel(2,"B0bar");
  fHistBSpecies->GetXaxis()->SetBinLabel(3,"B+");
  fHistBSpecies->GetXaxis()->SetBinLabel(4,"B-");
  fHistBSpecies->GetXaxis()->SetBinLabel(5,"B*+");
  fHistBSpecies->GetXaxis()->SetBinLabel(6,"B*-");
  fHistBSpecies->GetXaxis()->SetBinLabel(7,"Bs+");
  fHistBSpecies->GetXaxis()->SetBinLabel(8,"Bs-");
  fHistBSpecies->GetXaxis()->SetBinLabel(9,"Lb+");
  fHistBSpecies->GetXaxis()->SetBinLabel(10,"Lb-");
  fHistBSpecies->SetMinimum(0);
  fOutput->Add(fHistBSpecies);
  fHistLcDecayChan=new TH1F("hLcDecayChan","",9,-2.5,6.5);
  fHistLcDecayChan->GetXaxis()->SetBinLabel(1,"Violates p cons");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(2,"Other decay");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(3,"Error");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(4,"pK#pi non res");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(5,"pK#pi via K*0");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(6,"pK#pi via #Delta++");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(7,"pK#pi via #Lambda1520");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(8,"pK0s#rightarrowp#pi#pi");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(9,"#pi#Lambda#rightarrowp#pi#pi");
  fHistLcDecayChan->SetMinimum(0);
  fOutput->Add(fHistLcDecayChan);

  fHistNcollHFtype=new TH2F("hNcollHFtype","",5,-1.5,3.5,30,-0.5,29.5);
  fOutput->Add(fHistNcollHFtype);
  fHistNinjectedvsb=new TH2F("hNinjectedvsb","",50,0.,20.,140,-0.5,139.5);
  fOutput->Add(fHistNinjectedvsb);

  Double_t binseta[11]={-1.0,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.0};
  const Int_t nBinsPhi=40;
  Double_t binsphi[nBinsPhi+1];
  for(Int_t ib=0; ib<=nBinsPhi; ib++) binsphi[ib]=ib*TMath::Pi()/20.;
  const Int_t nBinsPt=24;  
  Double_t binspt[nBinsPt+1]={0.,0.10,0.15,0.2,0.25,
			      0.3,0.4,0.5,0.6,0.7,
			      0.8,0.9,1.,1.25,1.5,
			      1.75,2.,2.5,3.,4.,
			      5.,7.5,10.,15.,20.};

   fHistEtaPhiPtGenEle=new TH3F("hEtaPhiPtGenEle","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtGenEle);
  fHistEtaPhiPtGenPi=new TH3F("hEtaPhiPtGenPi","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtGenPi);
  fHistEtaPhiPtGenK=new TH3F("hEtaPhiPtGenK","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtGenK);
  fHistEtaPhiPtGenPro=new TH3F("hEtaPhiPtGenPro","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtGenPro);


  fHistEtaPhiPtRecEle=new TH3F("hEtaPhiPtRecEle","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtRecEle);
  fHistEtaPhiPtRecPi=new TH3F("hEtaPhiPtRecPi","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtRecPi);
  fHistEtaPhiPtRecK=new TH3F("hEtaPhiPtRecK","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtRecK);
  fHistEtaPhiPtRecPro=new TH3F("hEtaPhiPtRecPro","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtRecPro);

  fHistPtRecVsPtGen=new TH2F("hPtRecVsPtGen","  ; particle p_{T}  (GeV/c); track p_{T} (GeV/c)",100,0.,10.,100,0.,10.);
  fOutput->Add(fHistPtRecVsPtGen);
  fHistPhiRecVsPhiGen=new TH2F("hPhiRecVsPhiGen","  ; particle #varphi ; track #varphi",100,0.,2.*TMath::Pi(),100,0.,2.*TMath::Pi());
  fOutput->Add(fHistPhiRecVsPhiGen);
  fHistEtaRecVsEtaGen=new TH2F("hEtaRecVsEtaGen","  ; particle #eta ; track #eta",100,0.,2.*TMath::Pi(),100,0.,2.*TMath::Pi());
  fOutput->Add(fHistEtaRecVsEtaGen);

  fHistPtRecGood=new TH1F("hPtRecGood"," ; track p_{T} (GeV/c)",100,0.,10.);
  fOutput->Add(fHistPtRecGood);
  fHistPtRecFake=new TH1F("hPtRecFake"," ; track p_{T} (GeV/c)",100,0.,10.);
  fOutput->Add(fHistPtRecFake);

  PostData(1,fOutput);

}
//______________________________________________________________________________
void AliAnalysisTaskCheckHFMCProd::UserExec(Option_t *)
{
  //

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD = man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  Bool_t isAOD = man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  if(isESD)
    fEvent = dynamic_cast<AliESDEvent*>(InputEvent());
  else if(isAOD)
    fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  else
    fEvent = nullptr;

  if(!fEvent) {
    printf("AliAnalysisTaskCheckHFMCProd::Exec(): bad event\n");
    return;
  } 

  fHistoNEvents->Fill(0);

  if(!fESDtrackCuts){
    Int_t year=2011;
    if(fEvent->GetRunNumber()<=139517) year=2010;
    if(year==2010) fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE);
    else fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
    fESDtrackCuts->SetMaxDCAToVertexXY(2.4);
    fESDtrackCuts->SetMaxDCAToVertexZ(3.2);
    fESDtrackCuts->SetDCAToVertex2D(kTRUE);
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  }

  Int_t nTracks=fEvent->GetNumberOfTracks();
  fHistoTracks->Fill(nTracks);
  Int_t nSelTracks=0;
  AliVTrack* tr = nullptr;
  for(Int_t it=0; it<nTracks; it++){
    if(isESD)
      tr = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(it));
    else
      tr = dynamic_cast<AliAODTrack*>(fEvent->GetTrack(it));

    UInt_t status=tr->GetStatus();
    if(!(status&AliESDtrack::kITSrefit)) continue;
    if(!(status&AliESDtrack::kTPCin)) continue;
    nSelTracks++;
  }
  fHistoSelTracks->Fill(nSelTracks);

  Int_t nTracklets=0;
  Int_t nTracklets1=0;
  if(isESD)
  {
    const AliMultiplicity* mult=dynamic_cast<AliESDEvent*>(fEvent)->GetMultiplicity();
    nTracklets=mult->GetNumberOfTracklets();
    for(Int_t it=0; it<nTracklets; it++){
      Double_t eta=TMath::Abs(mult->GetEta(it));
      if(eta<1) nTracklets1++;
    }
  }
  else
  {
    AliAODTracklets* tracklets=dynamic_cast<AliAODEvent*>(fEvent)->GetTracklets();
    nTracklets=tracklets->GetNumberOfTracklets();
    nTracklets1=AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(dynamic_cast<AliAODEvent*>(fEvent), -1., 1.);
  }
  fHistoTracklets->Fill(nTracklets);
  fHistoTrackletsEta1->Fill(nTracklets1);
  
  AliVVertex *spdv = nullptr;
  if(isESD)
    spdv=(AliESDVertex*)(dynamic_cast<AliESDEvent*>(fEvent)->GetVertex());
  else
    spdv=dynamic_cast<AliAODEvent*>(fEvent)->GetPrimaryVertexSPD();
  if(spdv && spdv->IsFromVertexer3D()){
    fHistoSPD3DVtxX->Fill(spdv->GetX());
    fHistoSPD3DVtxY->Fill(spdv->GetY());
    fHistoSPD3DVtxZ->Fill(spdv->GetZ());
  }
  if(spdv && spdv->IsFromVertexerZ()){
    fHistoSPDZVtxX->Fill(spdv->GetX());
    fHistoSPDZVtxY->Fill(spdv->GetY());
    fHistoSPDZVtxZ->Fill(spdv->GetZ());
  }
  AliVVertex *trkv = nullptr;
  if(isESD)
    trkv = (AliESDVertex*)(dynamic_cast<AliESDEvent*>(fEvent)->GetPrimaryVertex());
  else
    trkv = dynamic_cast<AliAODEvent*>(fEvent)->GetPrimaryVertex();
  if(trkv && trkv->GetNContributors()>1){
    fHistoTRKVtxX->Fill(trkv->GetX());
    fHistoTRKVtxY->Fill(trkv->GetY());
    fHistoTRKVtxZ->Fill(trkv->GetZ());
  }

  AliMCEvent* mcEvent = nullptr;
  AliAODMCHeader* mcHeader = nullptr;
  TClonesArray* arrayMC = nullptr;
  if (fReadMC)
  {
    Double_t mcVtx[3] = {-999., -999., -999.};
    if(isESD)
    {
      AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler *>(man->GetMCtruthEventHandler());
      if (!eventHandler)
      {
        Printf("ERROR: Could not retrieve MC event handler");
        return;
      }
      mcEvent = eventHandler->MCEvent();
      if (!mcEvent)
      {
        Printf("ERROR: Could not retrieve MC event");
        return;
      }

      const AliVVertex *mcVert = mcEvent->GetPrimaryVertex();
      if (!mcVert)
      {
        Printf("ERROR: generated vertex not available");
        return;
      }

      mcVtx[0] = mcVert->GetX();
      mcVtx[1] = mcVert->GetY();
      mcVtx[2] = mcVert->GetZ();
    }
    else
    {
      arrayMC = (TClonesArray *)(dynamic_cast<AliAODEvent*>(fEvent)->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
      if (!arrayMC)
      {
        Printf("ERROR: MC particles branch not found");
        return;
      }

      // load MC header
      mcHeader = (AliAODMCHeader *)(dynamic_cast<AliAODEvent*>(fEvent)->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
      if (!mcHeader)
      {
        printf("ERROR: MC header branch not found");
        return;
      }
      mcHeader->GetVertex(mcVtx);
    }
  
    if (TMath::Abs(mcVtx[2]) > 10)
      return;

    TString genname = "";
    if(isESD)
      genname = mcEvent->GenEventHeader()->ClassName();

    Int_t nColl = -1;
    Double_t imppar = -999.;
    Int_t nInjected = 0;
    Int_t typeHF = -1;
    TList *lgen = nullptr;
    if ((isESD && genname.Contains("CocktailEventHeader")) || isAOD)
    {
      if(isESD)
      {
        AliGenCocktailEventHeader *cockhead = (AliGenCocktailEventHeader *)mcEvent->GenEventHeader();
        lgen = cockhead->GetHeaders();
      }
      else
      {
        lgen = mcHeader->GetCocktailHeaders();
      }

      for (Int_t ig = 0; ig < lgen->GetEntries(); ig++)
      {
        AliGenerator *gen = (AliGenerator *)lgen->At(ig);
        TString title = gen->GetName();
        if (title.Contains("bchadr"))
        {
          typeHF = 1;
          nInjected++;
        }
        else if (title.Contains("chadr"))
        {
          typeHF = 0;
          nInjected++;
        }
        else if (title.Contains("bele"))
        {
          typeHF = 3;
          nInjected++;
        }
        else if (title.Contains("cele"))
        {
          typeHF = 2;
          nInjected++;
        }
        else if (title.Contains("pythiaHF"))
        {
          nInjected++;
        }
        else if (title.Contains("hijing") || title.Contains("Hijing"))
        {
          AliGenHijingEventHeader *hijh = (AliGenHijingEventHeader *)lgen->At(ig);
          imppar = hijh->ImpactParameter();
        }
      }
      nColl = lgen->GetEntries();
      fHistNcollHFtype->Fill(typeHF, nColl);
      fHistNinjectedvsb->Fill(imppar, nInjected);
    }
    else if (genname.Contains("HijingEventHeader"))
    {
      AliGenHijingEventHeader *hijh = (AliGenHijingEventHeader *)mcEvent->GenEventHeader();
      imppar = hijh->ImpactParameter();
    }
    else
    {
      TString genTitle = mcEvent->GenEventHeader()->GetTitle();
      if (genTitle.Contains("bchadr"))
        typeHF = 1;
      else if (genTitle.Contains("chadr"))
        typeHF = 0;
      else if (genTitle.Contains("bele"))
        typeHF = 3;
      else if (genTitle.Contains("cele"))
        typeHF = 2;
      fHistNcollHFtype->Fill(typeHF, 1.);
    }

    Int_t nParticles = 0;
    if(isESD)
      nParticles = mcEvent->GetNumberOfTracks();
    else
      nParticles = arrayMC->GetEntriesFast();

    Double_t dNchdy = 0.;
    Int_t nb = 0, nc = 0;
    Int_t nCharmed = 0;
    Int_t nPhysPrim = 0;
    Int_t nPiKPeta09 = 0;
    for (Int_t i = 0; i < nParticles; i++)
    {
      AliVParticle *mcPart = nullptr;
      if(isESD)
      {
        mcPart = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(i));
        TParticle *part = (TParticle *)mcEvent->Particle(i);
        if (!part || !mcPart)
          continue;
      }
      else
      {
        mcPart = dynamic_cast<AliAODMCParticle *>(arrayMC->At(i));
      }

      Int_t pdg = mcPart->PdgCode();
      Int_t absPdg = TMath::Abs(pdg);
      if (absPdg == 4)
        nc++;
      if (absPdg == 5)
        nb++;
      
      Double_t pt = mcPart->Pt();
      Double_t eta = mcPart->Eta();
      Double_t phi = mcPart->Phi();
      Double_t energy = mcPart->E();
      Double_t pz = mcPart->Pz();

      if ((isESD && mcEvent->IsPhysicalPrimary(i)) || (isAOD && dynamic_cast<AliAODMCParticle *>(mcPart)->IsPhysicalPrimary()))
      {
        fHistoEtaPhysPrim->Fill(eta);
        if (absPdg == 11)
          fHistEtaPhiPtGenEle->Fill(eta, phi, pt);
        else if (absPdg == 211)
          fHistEtaPhiPtGenPi->Fill(eta, phi, pt);
        else if (absPdg == 321)
          fHistEtaPhiPtGenK->Fill(eta, phi, pt);
        else if (absPdg == 2212)
          fHistEtaPhiPtGenPro->Fill(eta, phi, pt);

        if (TMath::Abs(eta) < 0.5)
        {
          dNchdy += 0.6666; // 2/3 for the ratio charged/all
          nPhysPrim++;
        }
        if (TMath::Abs(eta) < 0.9)
        {
          fHistoPtPhysPrim->Fill(pt);
          if (absPdg == 211 || absPdg == 321 || absPdg == 2212)
          {
            nPiKPeta09++;
          }
        }
      }
      Float_t rapid = -999.;
      if (TMath::Abs(energy-TMath::Abs(pz))>0.001)
      {
        rapid = 0.5 * TMath::Log((energy + pz) / (energy - pz));
      }
      Int_t iPart = -1;
      Int_t iType = 0;
      Int_t iSpecies = -1;
      Int_t dummy[4];
      if (absPdg == 421)
      {
        iSpecies = 0;
        iType = isESD ? AliVertexingHFUtils::CheckD0Decay(mcEvent, i, dummy) : AliVertexingHFUtils::CheckD0Decay(arrayMC, dynamic_cast<AliAODMCParticle*>(mcPart), dummy);
        if (iType > 0)
          iPart = 0;
      }
      else if (absPdg == 411)
      {
        iSpecies = 1;
        iType = isESD ? AliVertexingHFUtils::CheckDplusDecay(mcEvent, i, dummy) : AliVertexingHFUtils::CheckDplusDecay(arrayMC, dynamic_cast<AliAODMCParticle*>(mcPart), dummy);
        if (iType < 0)
        {
          Int_t iTypeKKpi = isESD ? AliVertexingHFUtils::CheckDplusKKpiDecay(mcEvent, i, dummy) : AliVertexingHFUtils::CheckDplusKKpiDecay(arrayMC, dynamic_cast<AliAODMCParticle*>(mcPart), dummy);
          if (iTypeKKpi > 0)
            iType = 3;
        }
        if (iType > 0)
          iPart = 1;
      }
      else if (absPdg == 413)
      {
        iSpecies = 2;
        iType = isESD ? AliVertexingHFUtils::CheckDstarDecay(mcEvent, i, dummy) : AliVertexingHFUtils::CheckDstarDecay(arrayMC, dynamic_cast<AliAODMCParticle*>(mcPart), dummy);
        if (iType > 0)
          iPart = 2;
      }
      else if (absPdg == 431)
      {
        iSpecies = 3;
        iType = isESD ? AliVertexingHFUtils::CheckDsDecay(mcEvent, i, dummy) : AliVertexingHFUtils::CheckDsDecay(arrayMC, dynamic_cast<AliAODMCParticle*>(mcPart), dummy);
        if (iType == 1 || iType == 2)
          iPart = 3;
      }
      else if (absPdg == 4122)
      {
        iSpecies = 4;
        iType = isESD ? AliVertexingHFUtils::CheckLcpKpiDecay(mcEvent, i, dummy) : AliVertexingHFUtils::CheckLcpKpiDecay(arrayMC, dynamic_cast<AliAODMCParticle*>(mcPart), dummy);
        if (iType < 0)
        {
          Int_t iTypeV0 = isESD ? AliVertexingHFUtils::CheckLcV0bachelorDecay(mcEvent, i, dummy) : AliVertexingHFUtils::CheckLcV0bachelorDecay(arrayMC, dynamic_cast<AliAODMCParticle*>(mcPart), dummy);
          if (iTypeV0 == 1)
            iType = 5;
          if (iTypeV0 == 2)
            iType = 6;
        }
        fHistLcDecayChan->Fill(iType);
        if (iType >= 0)
          iPart = 4;
      }
      if (iSpecies >= 0)
        fHistYPtAllDecay[iSpecies]->Fill(pt, rapid);

      // check beauty mesons
      if (absPdg == 511)
        fHistBYPtAllDecay[0]->Fill(pt, rapid);
      else if (absPdg == 521)
        fHistBYPtAllDecay[1]->Fill(pt, rapid);
      else if (absPdg == 513)
        fHistBYPtAllDecay[2]->Fill(pt, rapid);
      else if (absPdg == 531)
        fHistBYPtAllDecay[3]->Fill(pt, rapid);
      else if (absPdg == 5122)
        fHistBYPtAllDecay[4]->Fill(pt, rapid);

      if (pdg == 511)
        fHistBSpecies->Fill(0);
      else if (pdg == -511)
        fHistBSpecies->Fill(1);
      else if (pdg == 521)
        fHistBSpecies->Fill(2);
      else if (pdg == -521)
        fHistBSpecies->Fill(3);
      else if (pdg == 513)
        fHistBSpecies->Fill(4);
      else if (pdg == -513)
        fHistBSpecies->Fill(5);
      else if (pdg == 531)
        fHistBSpecies->Fill(6);
      else if (pdg == -531)
        fHistBSpecies->Fill(7);
      else if (pdg == 5122)
        fHistBSpecies->Fill(8);
      else if (pdg == -5122)
        fHistBSpecies->Fill(9);

      if (iSpecies < 0)
        continue; // not a charm meson

      if (pdg == 421)
        fHistDSpecies->Fill(0);
      else if (pdg == -421)
        fHistDSpecies->Fill(1);
      else if (pdg == 411)
        fHistDSpecies->Fill(2);
      else if (pdg == -411)
        fHistDSpecies->Fill(3);
      else if (pdg == 413)
        fHistDSpecies->Fill(4);
      else if (pdg == -413)
        fHistDSpecies->Fill(5);
      else if (pdg == 431)
        fHistDSpecies->Fill(6);
      else if (pdg == -431)
        fHistDSpecies->Fill(7);
      else if (pdg == 4122)
        fHistDSpecies->Fill(8);
      else if (pdg == -4122)
        fHistDSpecies->Fill(9);

      Double_t distx = mcPart->Xv() - mcVtx[0];
      Double_t disty = mcPart->Yv() - mcVtx[1];
      Double_t distz = mcPart->Zv() - mcVtx[2];
      Double_t distToVert = TMath::Sqrt(distx * distx + disty * disty + distz * distz);
      Double_t distToVertXY = TMath::Sqrt(distx * distx + disty * disty);
      AliVParticle* mcDau0 = nullptr;
      Int_t iDau0=mcPart->GetDaughterFirst();
      if(iDau0>=0){
	if(isESD) mcDau0 = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(iDau0));
	else mcDau0 = dynamic_cast<AliAODMCParticle *>(arrayMC->At(iDau0));
      }
      Double_t declen3D=-9999.;
      Double_t declenXY=-9999.;
      if(mcDau0){
	Double_t decdistx= mcDau0->Xv() - mcVtx[0];
	Double_t decdisty  = mcDau0->Yv() - mcVtx[1];
	Double_t decdistz = mcDau0->Zv() - mcVtx[2];
	declen3D=TMath::Sqrt(decdistx * decdistx + decdisty * decdisty + decdistz * decdistz);
	declenXY=TMath::Sqrt(decdistx * decdistx + decdisty * decdisty);
      }
      fHistMotherID->Fill(mcPart->GetMother());
      Int_t iFromB = isESD ? AliVertexingHFUtils::CheckOrigin(mcEvent, dynamic_cast<AliMCParticle*>(mcPart), fSearchUpToQuark) : AliVertexingHFUtils::CheckOrigin(arrayMC, dynamic_cast<AliAODMCParticle*>(mcPart), fSearchUpToQuark);
      if (iFromB == 4)
      {
        fHistYPtPromptAllDecay[iSpecies]->Fill(pt, rapid);
        fHistOriginPrompt->Fill(distToVert);
	if(mcDau0){
	  fHistPtDDecLenPrompt[iSpecies]->Fill(pt,declen3D);
	  fHistPtDDecLenXYPrompt[iSpecies]->Fill(pt,declenXY);
	}
      }
      else if (iFromB == 5)
      {
	Double_t ptB= isESD ? AliVertexingHFUtils::GetBeautyMotherPt(mcEvent, dynamic_cast<AliMCParticle*>(mcPart)) : AliVertexingHFUtils::GetBeautyMotherPt(arrayMC, dynamic_cast<AliAODMCParticle*>(mcPart));
        fHistYPtFeeddownAllDecay[iSpecies]->Fill(pt, rapid);
        fHistOriginFeeddown->Fill(distToVert);
	fHistPtBDecLenBXYFeeddown->Fill(ptB,distToVertXY);
	if(mcDau0){
	  fHistPtDPtBDecLenFeeddown[iSpecies]->Fill(pt,ptB,declen3D);
	  fHistPtDPtBDecLenXYFeeddown[iSpecies]->Fill(pt,ptB,declenXY);
	}
      }

      if (iPart < 0)
        continue;
      if (iType < 0)
        continue;
      nCharmed++;
      if (iPart == 0 && iType > 0 && iType <= 2)
        fHistYPtD0byDecChannel[iType - 1]->Fill(pt, rapid);
      else if (iPart == 1 && iType > 0 && iType <= 3)
        fHistYPtDplusbyDecChannel[iType - 1]->Fill(pt, rapid);
      else if (iPart == 3 && iType > 0 && iType <= 2)
        fHistYPtDsbyDecChannel[iType - 1]->Fill(pt, rapid);

      if (iFromB == 4 && iPart >= 0 && iPart < 5)
        fHistYPtPrompt[iPart]->Fill(pt, rapid);
      else if (iFromB == 5 && iPart >= 0 && iPart < 5)
        fHistYPtFeeddown[iPart]->Fill(pt, rapid);
    }

    AliESDtrack* track = nullptr;
    for(Int_t i=0; i<nTracks; i++){
      if(isESD)
        track = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(i));
      else
      {
        AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(fEvent->GetTrack(i));
        // convert to ESD track here
        track = new AliESDtrack(aodTrack);
        // set the TPC cluster info
        track->SetTPCClusterMap(aodTrack->GetTPCClusterMap());
        track->SetTPCSharedMap(aodTrack->GetTPCSharedMap());
        track->SetTPCPointsF(aodTrack->GetTPCNclsF());
      }

      if(fESDtrackCuts->AcceptTrack(track)){
        if(track->GetLabel()>0) fHistPtRecGood->Fill(track->Pt());
        else fHistPtRecFake->Fill(track->Pt());
        Int_t label=TMath::Abs(track->GetLabel());

        AliVParticle *mcPart = nullptr;
        if(isESD)
          mcPart = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(label));
        else
          mcPart = dynamic_cast<AliAODMCParticle *>(arrayMC->At(label));

        if(!mcPart)
          continue;

        if((isESD && mcEvent->IsPhysicalPrimary(label)) || (isAOD && (dynamic_cast<AliAODMCParticle *>(mcPart)->IsPhysicalPrimary()))){
          Int_t absPdg=TMath::Abs(mcPart->PdgCode());
          Double_t pt = mcPart->Pt();
          Double_t eta = mcPart->Eta();
          Double_t phi = mcPart->Phi();
          if(absPdg==11) fHistEtaPhiPtRecEle->Fill(eta,phi,pt);
          else if(absPdg==211) fHistEtaPhiPtRecPi->Fill(eta,phi,pt);
          else if(absPdg==321) fHistEtaPhiPtRecK->Fill(eta,phi,pt);
          else if(absPdg==2212) fHistEtaPhiPtRecPro->Fill(eta,phi,pt);
          fHistPtRecVsPtGen->Fill(pt,track->Pt());
          fHistPhiRecVsPhiGen->Fill(phi,track->Phi());
          fHistEtaRecVsEtaGen->Fill(eta,track->Eta());
        }
      }
      if(isAOD)
      {
        delete track;
        track=nullptr;
      }
    }
    fHistoNcharmed->Fill(dNchdy,nCharmed);
    fHistoNbVsNc->Fill(nc,nb);
    fHistoPhysPrim->Fill(nPhysPrim);
    fHistoPhysPrimPiKPi09->Fill(nPiKPeta09);
    fHistoPhysPrimPiKPi09vsb->Fill(imppar,nPiKPeta09);
  }

  PostData(1,fOutput); 
}

//______________________________________________________________________________
void AliAnalysisTaskCheckHFMCProd::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  return;
}
