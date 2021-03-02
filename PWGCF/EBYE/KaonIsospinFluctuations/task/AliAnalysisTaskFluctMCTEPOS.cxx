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



#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TList.h"
#include "TPDGCode.h"
#include "THnSparse.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliVTrack.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliMultSelectionBase.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliAODPid.h"




#include "AliMCParticle.h"

#include "AliInputEventHandler.h"
#include "AliAODpidUtil.h"
#include "AliAnalysisTaskFluctMCTEPOS.h"

ClassImp(AliAnalysisTaskFluctMCTEPOS)

//________________________________________________________________________
AliAnalysisTaskFluctMCTEPOS::AliAnalysisTaskFluctMCTEPOS() // All data members should be initialised here
:AliAnalysisTaskSE(),fIsMonteCarlo(true), fCentrality(0),fMCImpactParameter(0),
fOutput(0),
fPIDResponse(0),
fESD(0), 
fPrimaryVtx(0), 
fTrackBuffSize(18000),
fHistGoodEvent(0),
fHistPt(0), 
fHistEta(0),
fHistNV0(0),
fHistBBK0Pos(0),
fHistBBK0Neg(0),
fHistBBPion(0),
fHistZVertexCent(0),
fHistCosPaMK0(0),	
fHistcTauMK0(0),		
fHistDcaMK0(0),	
fHistRapMK0(0),	
fHistArmPodK0(0),
fHistoCorrelation(0),
fKshortSparse(0)



{
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskFluctMCTEPOS::AliAnalysisTaskFluctMCTEPOS(const char *name) // All data members should be initialised here
  :AliAnalysisTaskSE(name), fIsMonteCarlo(true),  fCentrality(0),fMCImpactParameter(0),
fOutput(0),
fPIDResponse(0),
fESD(0), 
fPrimaryVtx(0), 
fTrackBuffSize(18000),
fHistGoodEvent(0),
fHistPt(0), 
fHistEta(0),
fHistNV0(0),
fHistBBK0Pos(0),
fHistBBK0Neg(0),
fHistBBPion(0),
fHistZVertexCent(0), 
fHistCosPaMK0(0),	
fHistcTauMK0(0),		
fHistDcaMK0(0),	
fHistRapMK0(0),	
fHistArmPodK0(0),
fHistoCorrelation(0),
fKshortSparse(0)
{

    DefineOutput(1, TList::Class());                  // for output list
}

//________________________________________________________________________
AliAnalysisTaskFluctMCTEPOS::~AliAnalysisTaskFluctMCTEPOS()
{

    if (fOutput) delete fOutput;    

}

//________________________________________________________________________
void AliAnalysisTaskFluctMCTEPOS::UserCreateOutputObjects()
{
   
  fOutput = new TList();
  fOutput->SetOwner();  
	

  fHistGoodEvent = new TH1F("hGoodEvent","No of events passing the cuts.",10,   .5,9.5);
    
      
    Int_t ptbins = 100;
    Float_t ptlow = 0.1, ptup = 40.0;
    fHistPt = new TH1F("fHistPt", "P_{T} distribution for reconstructed K0s", ptbins, ptlow, ptup);
    fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fHistPt->SetMarkerStyle(kFullCircle);
	
    Int_t etabins = 40;
    Float_t etalow = -1.0, etaup = 1.0;
    fHistEta = new TH1F("fHistEta","#eta distribution for reconstructed K0s",etabins, etalow, etaup);
    fHistEta->GetXaxis()->SetTitle("#eta");
    fHistEta->GetYaxis()->SetTitle("counts");
	
        
    Int_t K0Bins[4]   = {100, 300,  140, 24 };
    Double_t K0Min[4] = {-0.5,  0.,   0.414, -1.2};
    Double_t K0Max[4] = { 99.5 ,30.,0.582,1.2};
    
    const Char_t *K0Title[] = {"centrality","Pt","Mass","Eta"};
    
    fKshortSparse = new THnSparseD("fThnKshort", "", 4, K0Bins, K0Min, K0Max);
    for (Int_t iaxis = 0; iaxis < 4; iaxis++)
    fKshortSparse->GetAxis(iaxis)->SetTitle(K0Title[iaxis]);

    Int_t fgSparseDataBins[5]   = {2000,  2000,   2000,   2000, 100};
    Double_t fgSparseDataMin[5] = {-0.5, -0.5,   -0.5,   -0.5, -0.5};
    Double_t fgSparseDataMax[5] = {1999.5, 1999.5, 1999.5, 1999.5, 99.5};
	
    const Char_t *fgkSparseDataTitle[] = {"K_{s}^{0}", "K+", "K-", "kaons","centrality"};
    
	 
    fHistoCorrelation = new THnSparseD("fThnCorr", "", 5, fgSparseDataBins, fgSparseDataMin, fgSparseDataMax);
    for (Int_t iaxis = 0; iaxis < 5; iaxis++)
    fHistoCorrelation->GetAxis(iaxis)->SetTitle(fgkSparseDataTitle[iaxis]);
	
	

    //Fill histos here	
   
    fOutput->Add(fHistPt);
    //fOutput->Add(fHistEta);
    fOutput->Add(fHistoCorrelation);
    
    PostData(1, fOutput); 
}



//________________________________________________________________________

  static Bool_t AcceptMCTrack(AliMCParticle *mctrk){
  if(!mctrk) return kFALSE;
 
 
  if (mctrk->Eta() < -0.5 || mctrk->Eta() > 0.5) return kFALSE;
  if (mctrk->Pt() < 0.2 || mctrk->Pt() > 1.5 )    return kFALSE;
  
  return kTRUE;
  }  



//________________________________________________________________________
void AliAnalysisTaskFluctMCTEPOS::UserExec(Option_t *) 
{
 
  AliESDEvent *lESDevent = 0x0;
  AliMCEvent  *lMCevent  = 0x0;
  AliStack    *lMCstack  = 0x0;
    
  AliVEvent* event = InputEvent();

  // Appropriate for ESD analysis!
  lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!lESDevent) {
    AliWarning("ERROR: lESDevent not available \n");
    return;
  }                                                   
  
  lMCevent = dynamic_cast<AliMCEvent *>(MCEvent());
  if (!lMCevent) {
    Printf("ERROR: Could not retrieve MC event \n");

    return;
  }
    
  lMCstack = lMCevent->Stack();
  if (!lMCstack) {
    Printf("ERROR: Could not retrieve MC stack \n");

    return;
  }
  Bool_t isINT7selected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
  if(!isINT7selected)return ;

  const AliESDVertex *fPrimaryVtx = lESDevent->GetPrimaryVertex();
  if (!fPrimaryVtx){
    printf ("ERROR: no primary vertex\n");
    return;
  }

  Double_t zv=fPrimaryVtx->GetZ();
  //  cout<<zv<<"   = vertex z "<<endl;    
  if (TMath::Abs(zv) > 10.0) return;
   
  //===========================Centrality calculation =====================
  
  AliMultSelection *MultSelection = (AliMultSelection*)lESDevent->FindListObject("MultSelection");
  if(!MultSelection) return;
  double centPercentile = MultSelection->GetMultiplicityPercentile("V0M");
    
  if(centPercentile < 0 || centPercentile >= 80) return;
  //cout<<centPercentile<<"   = Centrality "<<endl;  

    Int_t nkaons = 0;
    Int_t nkshort = 0;
    Int_t nkaonsp = 0;
    Int_t nkaonsn = 0;
    
       // next are the typical event/track loops:
        for(Int_t i = 0; i < lMCevent->GetNumberOfTracks(); i++)
        {
	  //  cout << ", trackId = " << i << endl;
	  AliMCParticle* mctrk = dynamic_cast<AliMCParticle*>(lMCevent->GetTrack(i));
	  if(!mctrk)
	    continue;
	  //cout << "pt = " << mctrk->Pt() << endl;
	  
	  TParticle * part = lMCstack->Particle(i);
	  if( !part )
	    continue;
            

	  if( !lMCstack->IsPhysicalPrimary(i) )
	    continue;

      Int_t code=part->GetPdgCode();

      if (TMath::Abs(code) == 321){
	nkaons++; }
      if (code == 321){
	nkaonsp++;}
      if( code == -321){
	nkaonsn++;}
      if(code == 310){
      nkshort++;}
       }
    
   Double_t vsparse[5];

   
   vsparse[0]   = nkshort;
   vsparse[1]   = nkaonsp;
   vsparse[2]   = nkaonsn;
   vsparse[3]   = nkaons;
   vsparse[4]   = centPercentile;
   
   fHistoCorrelation->Fill(vsparse);
   
	PostData(1, fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskFluctMCTEPOS::Terminate(Option_t *) 
{
	// Draw result to screen, or perform fitting, normalizations
	// Called once at the end of the query
	
	fOutput = dynamic_cast<TList*> (GetOutputData(1));
	if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
	
	fHistPt = dynamic_cast<TH1F*> (fOutput->FindObject("fHistPt"));
	if (!fHistPt) { Printf("ERROR: could not retrieve fHistPt"); return;}
	
	
	
	TCanvas *c = new TCanvas("AliAnalysisTaskFluctMCTEPOS","P_{T} & #eta",10,10,1020,510);
	c->Divide(2,1);
	c->cd(1)->SetLogy();
	fHistPt->DrawCopy("E");

}




