/**************************************************************************
* Copyright(c) 1998-2008,ALICE Experiment at CERN,All rights reserved.    *
*                                                                         *
* Author: The ALICE Off-line Project.                                     *
* Contributors are mentioned in the code where appropriate.               *
*                                                                         *
* Permission to use,copy,modify and distribute this software and its      *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee,provided that the above copyright notice appears in all     *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

/////////////////////////////////////////////////////
// AliAnalysisTaskFlowCascade:
// Analysis task to select Xi and Omega candidates for flow analysis.
//
// Author: Zhong-Bao.Yin@cern.ch
//////////////////////////////////////////////////////

#include "TChain.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TVector3.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliVParticle.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliVVertex.h"
#include "AliESDVZERO.h"
#include "AliESDUtils.h"

#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliPIDResponse.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAODVZERO.h"
#include "AliAODcascade.h"

#include "TMath.h"
#include "TObjArray.h"
#include "AliFlowCandidateTrack.h"

#include "AliFlowTrackCuts.h"
#include "AliFlowEventCuts.h"
#include "AliFlowEvent.h"
#include "AliFlowCommonConstants.h"

#include "AliAnalysisTaskFlowCascade.h"

ClassImp(AliAnalysisTaskFlowCascade)

//_____________________________________________________________________________
  AliAnalysisTaskFlowCascade::AliAnalysisTaskFlowCascade() :
    AliAnalysisTaskSE(),
    //    fMinCent(0), fMaxCent(0),
    fSpecie(0),
    fMassBins(0),
    fMinMass(0.0),
    fMaxMass(0.0),
    fCutsEvent(NULL),
    fCutsRPTPC(NULL),
    fCutsRPVZE(NULL),
    fCutsPOI(NULL),
    fCutsDau(NULL),
    fPIDResponse(NULL),
    fFlowEventTPC(NULL),
    fFlowEventVZE(NULL),
    fCandidates(NULL),
    fQAList(NULL)
{
  //ctor                                                                       
  for (Int_t i=0; i!=8; ++i)
    fCascadeCuts[i] = 0;
  
}

//_____________________________________________________________________________
AliAnalysisTaskFlowCascade
::AliAnalysisTaskFlowCascade(const char *name,
			     AliFlowEventCuts *cutsEvent, 
			     AliFlowTrackCuts *cutsRPTPC,
			     AliFlowTrackCuts *cutsRPVZE,
			     /* AliESDtrackCuts */ AliFlowTrackCuts *cutsDau ) :
  AliAnalysisTaskSE(name),
  //fMinCent(minCent), fMaxCent(maxCent),
  fSpecie(0),
  fMassBins(0),
  fMinMass(0.0),
  fMaxMass(0.0),
  fCutsEvent(cutsEvent),
  fCutsRPTPC(cutsRPTPC),
  fCutsRPVZE(cutsRPVZE),
  fCutsPOI(NULL),
  fCutsDau(cutsDau),
  fPIDResponse(NULL),
  fFlowEventTPC(NULL),
  fFlowEventVZE(NULL),
  fCandidates(NULL),
  fQAList(NULL)
{
  //ctor                                                                       
  for (Int_t i=0; i!=8; ++i)
    fCascadeCuts[i] = 0;

  DefineInput( 0,TChain::Class());
  DefineOutput(1,AliFlowEventSimple::Class()); // TPC object
  DefineOutput(2,AliFlowEventSimple::Class()); // VZE object
  DefineOutput(3,TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskFlowCascade::~AliAnalysisTaskFlowCascade()
{
  if(fQAList) delete fQAList;
  if (fFlowEventTPC) delete fFlowEventTPC;
  if (fFlowEventVZE) delete fFlowEventVZE;
  if (fCandidates)   delete fCandidates;
  if (fCutsDau)      delete fCutsDau;
  if (fCutsPOI)      delete fCutsPOI;
  if (fCutsRPTPC)   delete fCutsRPTPC;
  if (fCutsRPVZE)   delete fCutsRPVZE;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowCascade::UserCreateOutputObjects()
{

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler
    = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  
  fQAList = new TList();
  fQAList->SetOwner();
  AddQAEvents();
  AddQACandidates();
  //  PostData(3,fQAList);

  AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
  cc->SetNbinsMult(1);
  cc->SetMultMin(0);
  cc->SetMultMax(1);

  cc->SetNbinsPt(20);
  cc->SetPtMin(0.0);
  cc->SetPtMax(10.0);

  cc->SetNbinsPhi(1);
  cc->SetPhiMin(0.0);
  cc->SetPhiMax(TMath::TwoPi());

  cc->SetNbinsEta(1);
  cc->SetEtaMin(-2.0);
  cc->SetEtaMax(+2.0);

  cc->SetNbinsQ(3);
  cc->SetQMin(0.0);
  cc->SetQMax(3.0);

  cc->SetNbinsMass(fMassBins);
  cc->SetMassMin(fMinMass);
  cc->SetMassMax(fMaxMass);

  fCutsPOI = new AliFlowTrackCuts("null_cuts");
  fCutsPOI->SetParamType(fCutsRPTPC->GetParamType());
  fCutsPOI->SetPtRange(+1,-1); // select nothing QUICK   
  fCutsPOI->SetEtaRange(+1,-1); // select nothing VZERO  

  fFlowEventTPC = new AliFlowEvent(3000);
  fFlowEventVZE = new AliFlowEvent(1000);
  fCandidates = new TObjArray(100);
  fCandidates->SetOwner();

  PostData(1,fFlowEventTPC);
  PostData(2,fFlowEventVZE);
  PostData(3,fQAList);
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowCascade::AddQAEvents()
{
  TList *tQAEvents = new TList();
  tQAEvents->SetName("Events");
  tQAEvents->SetOwner();
  TH1I* tEvent = new TH1I("Event","Number of Events",   3,0,3);
  tQAEvents->Add(tEvent);
  
  TH1D *tTPCRFP = new TH1D("RFPTPC",
			   "TPC Reference Flow Particles;multiplicity",
			   100, 0, 3000); 
  tQAEvents->Add(tTPCRFP);
  TH1D *tVZERFP = new TH1D("RFPVZE", 
			   "VZERO Reference Flow Particles;multiplicity",
			   100, 0, 30000); 
  tQAEvents->Add(tVZERFP);

  TProfile *tCuts = new TProfile("Cuts","Analysis Cuts",10,0,10);
  tCuts->Fill(0.5,fCascadeCuts[0],1); 
  tCuts->GetXaxis()->SetBinLabel(1,"dcaXiDau");
  tCuts->Fill(1.5,fCascadeCuts[1],1); 
  tCuts->GetXaxis()->SetBinLabel(2,"XiCPA");
  tCuts->Fill(2.5,fCascadeCuts[2],1); 
  tCuts->GetXaxis()->SetBinLabel(3,"dcaV0Vtx");
  tCuts->Fill(3.5,fCascadeCuts[3],1); 
  tCuts->GetXaxis()->SetBinLabel(4,"dcaBachVtx");
  tCuts->Fill(4.5,fCascadeCuts[4],1); 
  tCuts->GetXaxis()->SetBinLabel(5,"dcaV0Dau");
  tCuts->Fill(5.5,fCascadeCuts[5],1); 
  tCuts->GetXaxis()->SetBinLabel(6,"V0CPA");
  tCuts->Fill(6.5,fCascadeCuts[6],1); 
  tCuts->GetXaxis()->SetBinLabel(7,"dcaV0DauVtx");
  tCuts->Fill(7.5,fCascadeCuts[7],1); 
  tCuts->GetXaxis()->SetBinLabel(8,"V0Mass");
  tQAEvents->Add(tCuts);

  fQAList->Add(tQAEvents);
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowCascade::AddQACandidates()
{
  TList *tQACandidates;
  
  tQACandidates = new TList();
  tQACandidates->SetOwner();
  tQACandidates->SetName("Candidates");

  TH1F* tChi2Xi = new TH1F("Chi2Xi", 
		       "Cascade #chi^{2}; #chi^{2}; Number of Cascades", 
		       160, 0, 160);
  tQACandidates->Add(tChi2Xi);

  TH1F* tDCAXiDaughters 
    = new TH1F( "DcaXiDaughters",  
		"DCA between Xi Daughters; DCA (cm); Number of Cascades", 
		100, 0., 0.5);
  tQACandidates->Add(tDCAXiDaughters);

  TH1F * tDCABachToPrimVertex
    = new TH1F("DcaBachToPrimVertex", 
	       "DCA of Bach. to Prim. Vertex; DCA (cm);Number of Cascades", 
	       250, 0., 2.5);
  tQACandidates->Add(tDCABachToPrimVertex);
  
  TH1F * tXiCosOfPointingAngle
    = new TH1F("XiCosOfPointingAngle",
	       "Cos of Xi Pointing Angle; Cos (Xi Point.Angl);Number of Xis", 
	       200, 0.99, 1.0);
  tQACandidates->Add(tXiCosOfPointingAngle);

  TH1F * tXiRadius = new TH1F("XiRadius",  
			      "Casc. decay transv. radius; r (cm); Counts" , 
			      1050, 0., 105.0 );
  tQACandidates->Add(tXiRadius);
  
  TH1F *tMassLambda 
    = new TH1F("MassLambdaAsCascDghter",
	       "#Lambda assoc. to Casc. candidates; Eff. Mass (GeV/c^{2}); Counts", 
	       300,1.00,1.3);
  tQACandidates->Add(tMassLambda);

  TH1F *tV0Chi2 = new TH1F("V0Chi2Xi", 
		       "V0 #chi^{2}, in cascade; #chi^{2};Counts", 
		       160, 0, 40);
  tQACandidates->Add(tV0Chi2);
  
  TH1F * tV0CosOfPointingAngle 
    = new TH1F("V0CosOfPointingAngleXi", 
	       "Cos of V0 Pointing Angle, in cascade;Cos(V0 Point. Angl); Counts", 
	       200, 0.98, 1.0);
  tQACandidates->Add(tV0CosOfPointingAngle);

  TH1F *tV0Radius  = new TH1F("V0RadiusXi", 
			  "V0 decay radius, in cascade; radius (cm); Counts", 
			  1050, 0., 105.0);
  tQACandidates->Add(tV0Radius);
  
  TH1F * tDcaV0DaughtersXi 
    = new TH1F("DcaV0DaughtersXi", 
	       "DCA between V0 daughters, in cascade;DCA (cm);Number of V0s", 
	       120, 0., 0.6);
  tQACandidates->Add(tDcaV0DaughtersXi);

  TH1F * tDcaV0ToPrimVertex 
    = new TH1F("DcaV0ToPrimVertexXi", 
	       "DCA of V0 to Prim. Vertex, in cascade;DCA (cm);Number of Cascades", 200, 0., 1.);
  tQACandidates->Add(tDcaV0ToPrimVertex);

  TH1F * tDCAPosToPrimVertex =
    new TH1F("DcaPosToPrimVertexXi", 
	     "DCA of V0 pos daughter to Prim. Vertex;DCA (cm);Counts", 
	     300, 0, 3);
  tQACandidates->Add(tDCAPosToPrimVertex);

  TH1F * tDCANegToPrimVertex 
    =  new TH1F("DcaNegToPrimVertexXi", 
		"DCA of V0 neg daughter to Prim. Vertex;DCA (cm);Counts", 
		300, 0, 3);
  tQACandidates->Add(tDCANegToPrimVertex);

  TH1F *tV0toXiCosOfPointingAngle 
    = new TH1F("V0toXiCosOfPointingAngle", 
	       "Cos. of V0 Ptng Angl Xi vtx; Cos(V0 Point. Angl / Xi vtx); Counts", 
	       100, 0.99, 1.0);
  tQACandidates->Add(tV0toXiCosOfPointingAngle);

  TH2F *th2Armenteros 
    = new TH2F("Armenteros", 
	       "#alpha_{Arm}(casc. cand.) Vs Pt_{Arm}(casc. cand.); #alpha_{Arm} ; Pt_{Arm} (GeV/c)", 
	       140, -1.2, 1.2, 300, 0., 0.3);
  tQACandidates->Add(th2Armenteros);

  TH2F *th2TPCdEdxOfCascDghters
    = new TH2F( "TPCdEdxOfCascDghters",
		"TPC dE/dx of the cascade daughters; charge x || #vec{p}_{TPC inner wall}(Casc. daughter) || (GeV/c); TPC signal (ADC) ", 
		200, -10.0, 10.0, 450, 0., 900.);
  tQACandidates->Add(th2TPCdEdxOfCascDghters);

  TH2F *th2MassVsPtAll 
    = new TH2F("MassVsPtAll", 
	       "M_{candidates} vs Pt; Pt (GeV/c); M (GeV/c^{2})", 
	       100, 0., 10., fMassBins, fMinMass, fMaxMass);
  tQACandidates->Add(th2MassVsPtAll);

  TH1F *tDistToVtxZAfter 
    = new TH1F("DistToVtxZAfter", 
	       "Distance to vtx z after propagation to vtx; z [cm]", 
	       100, -5., 5.);
  tQACandidates->Add(tDistToVtxZAfter);

  TH1F * tDistToVtxXYAfter 
    = new TH1F("DistToVtxXYAfter", 
	       "Distance to vtx xy after propagation to vtx",
	       500, 0., 50.);
  tQACandidates->Add(tDistToVtxXYAfter);

  TH2F *th2DistToVtxZBeforeVsAfter 
    = new TH2F("DistToVtxZBeforeVsAfter", 
	       "Distance to vtx z before vs after propagation; Distance before [cm]; Distance after [cm]", 
	       500, -50., 50., 100, -5., 5.);
  tQACandidates->Add(th2DistToVtxZBeforeVsAfter);

  TH2F *th2DistToVtxXYBeforeVsAfter
    = new TH2F("DistToVtxXYBeforeVsAfter",
	       "Distance to vtx xy before vs after propagation; Distance before [cm]; Distance after [cm]", 
	       500, 0., 50, 500, 0., 50);
  tQACandidates->Add(th2DistToVtxXYBeforeVsAfter);

  TH2F * th2PxBeforeVsAfter 
    = new TH2F("PxBeforeVsAfter", 
	       "Px before vs after propagation; Px [GeV/c]; Px' [GeV/c]", 
	       200, -10., 10, 200, -10., 10.);
  tQACandidates->Add(th2PxBeforeVsAfter);

  TH2F * th2PyBeforeVsAfter
    = new TH2F("PyBeforeVsAfter",
               "Py before vs after propagation; Py [GeV/c]; Py' [GeV/c]",
               200, -10., 10, 200, -10., 10.);
  tQACandidates->Add(th2PyBeforeVsAfter);
  
  TH2F * th2PhiPosBeforeVsAfter
    = new TH2F("PhiPosBeforeVsAfter", 
	       "Phi for positively charged candidates before vs after propagation; #phi; #phi'", 
	       360, 0., 2.0*TMath::Pi(), 360, 0., 2.0*TMath::Pi());
  tQACandidates->Add(th2PhiPosBeforeVsAfter);

  TH2F *th2PhiNegBeforeVsAfter
    = new TH2F("PhiNegBeforeVsAfter",
               "Phi for negatively charged candidates before vs after propagation; #phi; #phi'",
               360, 0., 2.0*TMath::Pi(), 360, 0., 2.0*TMath::Pi());
  tQACandidates->Add(th2PhiNegBeforeVsAfter);

  fQAList->Add(tQACandidates);
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowCascade::NotifyRun()
{
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowCascade::UserExec(Option_t *)
{
  AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  Bool_t acceptEvent=kFALSE; 
  fCandidates->SetLast(-1);

  if(fESD) {
    // recorrecting VZERO (for pass 1 only)
    /*
    Float_t *vChCorr = new Float_t[64];
    Float_t dummy;
    AliESDUtils::GetCorrV0(fESD,dummy,NULL,vChCorr);
    AliESDVZERO *vzero = (AliESDVZERO*) fESD->GetVZEROData();
    vzero->SetMultiplicity( vChCorr );
    delete [] vChCorr;
    */
    //

    ((TH1I*)((TList*)fQAList->FindObject("Events"))->FindObject("Event"))->Fill(0);
    
    const AliVVertex *vtxGlb = fESD->GetPrimaryVertexTracks();
    const AliVVertex *vtxSPD = fESD->GetPrimaryVertexSPD();
    if(!vtxGlb || !vtxSPD) return;
    if( fCutsEvent->IsSelected(fESD, 0) && (TMath::Abs(vtxSPD->GetZ()-vtxGlb->GetZ()) <= 0.5) ) {
      acceptEvent = kTRUE;
    ((TH1I*)((TList*)fQAList->FindObject("Events"))->FindObject("Event"))->Fill(2);
      ReadFromESDv0(fESD);
    }
  } else if(fAOD) {
    const AliVVertex *vtxGlb = fAOD->GetPrimaryVertex();
    const AliVVertex *vtxSPD = fAOD->GetPrimaryVertexSPD();
    if(!vtxGlb || !vtxSPD) return;

    ((TH1I*)((TList*)fQAList->FindObject("Events"))->FindObject("Event"))->Fill(0);
        
    if(fCutsEvent->IsSelected(fAOD, 0) && (TMath::Abs(vtxSPD->GetZ()-vtxGlb->GetZ()) <= 0.5) ) {
      acceptEvent = kTRUE;
      ((TH1I*)((TList*)fQAList->FindObject("Events"))->FindObject("Event"))->Fill(2);
      ReadFromAODv0(fAOD);
    }

    
    /*

    AliAODHeader *aodHeader = fAOD->GetHeader();
    if(!aodHeader) return;
    AliCentrality *centrality = aodHeader->GetCentralityP();
    if(!centrality) return;
    Double_t cent = centrality->GetCentralityPercentile("V0M" );
    Double_t cent1 = centrality->GetCentralityPercentile("TRK" );
    if(TMath::Abs(cent-cent1) >= 5.) return;
    
    if(cent<fMinCent||cent>=fMaxCent) return; //centrality cut
    
    Double_t zvtx = fAOD->GetPrimaryVertex()->GetZ();
    if(TMath::Abs(zvtx)>10.) return; //vertex cut
    
    ((TH1I*)((TList*)fQAList->FindObject("Events"))->FindObject("Event"))->Fill(2);
    ReadFromAODv0(fAOD);
    */  
  }
  
  if(!acceptEvent) return;

  ((TH1D*)((TList*)fQAList->FindObject("Events"))->FindObject("RFPTPC"))
    ->Fill(fFlowEventTPC->GetNumberOfRPs() );
  Double_t mult=0;
  for(Int_t i=0; i != fFlowEventVZE->GetNumberOfRPs(); ++i) {
    AliFlowTrackSimple *pTrack = fFlowEventVZE->GetTrack(i);
    mult += pTrack->Weight();
  }
  ((TH1D*)((TList*)fQAList->FindObject("Events"))->FindObject("RFPVZE"))
    ->Fill( mult );

  //  if(fDebug) printf("TPCevent %d | VZEevent %d\n",
  //                  fFlowEventTPC->NumberOfTracks(),
  //                  fFlowEventVZE->NumberOfTracks() );
  AddCandidates();

  PostData(1,fFlowEventTPC);
  PostData(2,fFlowEventVZE);
  PostData(3,fQAList);

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowCascade::AddCandidates(){
  
  //  if(fDebug) printf("I received %d candidates\n",
  //		    fCandidates->GetEntriesFast());
  for(int iCand=0; iCand!=fCandidates->GetEntriesFast(); ++iCand ) {
    AliFlowCandidateTrack *cand 
      = dynamic_cast<AliFlowCandidateTrack*>(fCandidates->At(iCand));
    if(!cand) continue;
    // if(fDebug) 
    //  printf(" >Checking at candidate %d with %d daughters: mass %f\n",
    //	     iCand, cand->GetNDaughters(), cand->Mass());
    
    // untagging ===>                      
    for(int iDau=0; iDau != cand->GetNDaughters(); ++iDau) {
      // if(fDebug) 
      // printf(" >Daughter %d with fID %d", iDau, cand->GetIDDaughter(iDau));
      for(int iRPs=0; iRPs != fFlowEventTPC->NumberOfTracks(); ++iRPs ) {
        AliFlowTrack *iRP 
	  = dynamic_cast<AliFlowTrack*>(fFlowEventTPC->GetTrack( iRPs ));
        if (!iRP) continue;
        if( !iRP->InRPSelection() ) continue;
        if( cand->GetIDDaughter(iDau) == iRP->GetID() ) {
          //if(fDebug) printf(" was in RP set");
          iRP->SetForRPSelection(kFALSE);
          fFlowEventTPC->SetNumberOfRPs( fFlowEventTPC->GetNumberOfRPs() -1 );
        }
      }
      //if(fDebug) printf("\n");
    }
    // <=== untagging
    cand->SetForPOISelection(kTRUE);
    fFlowEventTPC->InsertTrack( ((AliFlowTrack*) cand) );
    fFlowEventVZE->InsertTrack( ((AliFlowTrack*) cand) );
  }

  //  if(fDebug) printf("TPCevent %d | VZEevent %d\n",
  //                fFlowEventTPC->NumberOfTracks(),
  //                fFlowEventVZE->NumberOfTracks() );

}

//______________________________________________________________________________
void AliAnalysisTaskFlowCascade::ReadFromESDv0(AliESDEvent *fESD)
{
  //AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts("null_cuts");
  //cutsPOI->SetParamType( fCutsRP->GetParamType() );
  //cutsPOI->SetParamType( AliFlowTrackCuts::kGlobal );
  // cutsPOI->SetPtRange(+1,-1); // select nothing
  //cutsPOI->SetEtaRange(+1,-1); // select nothing VZERO

  fCutsRPTPC->SetEvent(fESD,MCEvent());
  fCutsRPVZE->SetEvent(fESD,MCEvent());

  fCutsPOI->SetEvent(fESD,MCEvent());

  fFlowEventTPC->Fill(fCutsRPTPC,fCutsPOI);
  fFlowEventVZE->Fill(fCutsRPVZE,fCutsPOI);
  
  Double_t trkPrimaryVtxPos[3] = {-100., -100., -100.};
  Double_t bestPrimaryVtxPos[3] = {-100., -100., -100.};
  int nCascades=fESD->GetNumberOfCascades();
  
  const AliESDVertex *primaryTrackingESDVtx = fESD->GetPrimaryVertexTracks();
  primaryTrackingESDVtx->GetXYZ(trkPrimaryVtxPos);
  
  const AliESDVertex *primaryBestESDVtx = fESD->GetPrimaryVertex();
  primaryBestESDVtx->GetXYZ(bestPrimaryVtxPos);

  Double_t b = fESD->GetMagneticField();

  for(int i = 0; i != nCascades; ++i) {
    
    // Double_t trkPrimaryVtxRadius3D = -500.;
    // Double_t bestPrimaryVtxRadius3D = -500.;
    Double_t effMassXi = 0.;
    Double_t chi2Xi = -1.;
    Double_t dcaXiDaughters = -1.;
    Double_t XiCosOfPointingAngle = -1.;
    Double_t posXi[3] = {-1000., -1000., -1000.};
    Double_t XiRadius = -1000.;

    // Double_t innerWallMomCascDghters[3] = {-100., -100., -100.};
    //Double_t tpcSignalCascDghters[3] = {-100., -100., -100.};

    Double_t invMassLambdaAsCascDghter = 0.;
    Double_t V0Chi2Xi = -1.;
    Double_t dcaV0DaughtersXi = -1.;
    
    Double_t dcaBachToPrimaryVtxXi = -1.;
    Double_t dcaV0ToPrimaryVtxXi = -1.;
    Double_t dcaPosToPrimaryVtxXi = -1.;
    Double_t dcaNegToPrimaryVtxXi = -1.;
    Double_t V0CosOfPointingAngleXi = -1.;
    Double_t posV0Xi[3] = {-1000., -1000., -1000.};
    Double_t V0RadiusXi = -1000.;
    Double_t V0quality = 0.;

    Double_t invMassXiMinus = 0.;
    Double_t invMassXiPlus = 0.;
    Double_t invMassOmegaMinus = 0.;
    Double_t invMassOmegaPlus = 0.;

    /*
    Bool_t isPosInXiProton = kFALSE;
    Bool_t isPosInXiPion = kFALSE;
    Bool_t isPosInOmegaProton = kFALSE;
    Bool_t isPosInOmegaPion = kFALSE;
    
    Bool_t isNegInXiProton = kFALSE;
    Bool_t isNegInXiPion = kFALSE;
    Bool_t isNegInOmegaProton = kFALSE;
    Bool_t isNegInOmegaPion = kFALSE;

    Bool_t isBachelorKaon = kFALSE;
    Bool_t isBachelorPion = kFALSE;
    */

    Bool_t isBachelorKaonForTPC = kFALSE;
    Bool_t isBachelorPionForTPC = kFALSE;
    Bool_t isNegPionForTPC = kFALSE;
    Bool_t isPosPionForTPC = kFALSE;
    Bool_t isNegProtonForTPC = kFALSE;
    Bool_t isPosProtonForTPC = kFALSE;

    Double_t XiPx = 0., XiPy = 0., XiPz = 0.;
    Double_t XiPt = 0.;
    Double_t XiPtot = 0.;
    
    Double_t bachPx = 0., bachPy = 0., bachPz = 0.;
    Double_t bachPt = 0.;
    Double_t bachPtot = 0.;
    
    //Short_t chargeXi = -2;
    Double_t V0toXiCosOfPointingAngle = 0.;
    
    Double_t rapXi = -20.;
    Double_t rapOmega = -20.;
    Double_t phi = 6.3;
    Double_t alphaXi = -200.;
    Double_t ptArmXi = -200.;
    //    TLorentzVector lv1, lv2, lv3, lv12, lvXi;

    Double_t distToVtxZBefore = -999.;
    Double_t distToVtxZAfter = -999.;
    Double_t distToVtxXYBefore = -999.;
    Double_t distToVtxXYAfter = -999.;
    Double_t XiPAfter[3] = {-999., -999., -999.};
    Double_t phiAfter = -999.;
    
    AliESDcascade *xi = fESD->GetCascade(i);
    if(!xi) continue;
    
    if(xi->Charge()<0)
      xi->ChangeMassHypothesis(V0quality, 3312); // Xi- hypothesis
    else if(xi->Charge() > 0)
      xi->ChangeMassHypothesis(V0quality, -3312);
    else continue;

    effMassXi = xi->GetEffMassXi();
    chi2Xi = xi->GetChi2Xi();
    dcaXiDaughters = xi->GetDcaXiDaughters();
    XiCosOfPointingAngle 
      = xi->GetCascadeCosineOfPointingAngle(bestPrimaryVtxPos[0],
					    bestPrimaryVtxPos[1],
					    bestPrimaryVtxPos[2]);
    xi->GetXYZcascade(posXi[0], posXi[1], posXi[2]);
    XiRadius = TMath::Sqrt(posXi[0]*posXi[0]
			   +posXi[1]*posXi[1]
			   +posXi[2]*posXi[2]);
    
    UInt_t idxPosXi = (UInt_t)TMath::Abs(xi->GetPindex());
    UInt_t idxNegXi = (UInt_t)TMath::Abs(xi->GetNindex());
    UInt_t idxBach = (UInt_t)TMath::Abs(xi->GetBindex());

    if(idxBach == idxPosXi || idxBach == idxNegXi) continue;

    AliESDtrack *pTrkXi = fESD->GetTrack(idxPosXi);
    AliESDtrack *nTrkXi = fESD->GetTrack(idxNegXi);
    AliESDtrack *bTrkXi = fESD->GetTrack(idxBach);

    if( !pTrkXi || !nTrkXi || !bTrkXi ) continue;

    if( !fCutsDau->IsSelected(pTrkXi) 
	|| !fCutsDau->IsSelected(nTrkXi)
	|| !fCutsDau->IsSelected(bTrkXi) ) continue;

    const AliExternalTrackParam *pExtTrk = pTrkXi->GetInnerParam();
    const AliExternalTrackParam *nExtTrk = nTrkXi->GetInnerParam();
    const AliExternalTrackParam *bExtTrk = bTrkXi->GetInnerParam();
    
    if(pExtTrk && pTrkXi->IsOn(AliESDtrack::kTPCin)){
      ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("TPCdEdxOfCascDghters"))->Fill(pExtTrk->GetP()*pExtTrk->Charge(), pTrkXi->GetTPCsignal());
    }
    if(nExtTrk && nTrkXi->IsOn(AliESDtrack::kTPCin)){
      ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("TPCdEdxOfCascDghters"))->Fill(nExtTrk->GetP()*nExtTrk->Charge(), nTrkXi->GetTPCsignal());
    }
    if(bExtTrk && bTrkXi->IsOn(AliESDtrack::kTPCin)){
      ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("TPCdEdxOfCascDghters"))->Fill(bExtTrk->GetP()*bExtTrk->Charge(), bTrkXi->GetTPCsignal());
    }
    
    invMassLambdaAsCascDghter = xi->GetEffMass(); // from V0
    dcaV0DaughtersXi = xi->GetDcaV0Daughters();
    V0Chi2Xi = xi->GetChi2V0();
    V0CosOfPointingAngleXi 
      = xi->GetV0CosineOfPointingAngle(bestPrimaryVtxPos[0],
				       bestPrimaryVtxPos[1],
				       bestPrimaryVtxPos[2]);
    dcaV0ToPrimaryVtxXi = xi->GetD(bestPrimaryVtxPos[0], 
				   bestPrimaryVtxPos[1],
				   bestPrimaryVtxPos[2]);
    dcaBachToPrimaryVtxXi = TMath::Abs(bTrkXi->GetD(bestPrimaryVtxPos[0],
						    bestPrimaryVtxPos[1],
						    bestPrimaryVtxPos[2]));
    
    //V0
    xi->GetXYZ(posV0Xi[0], posV0Xi[1], posV0Xi[2]);
    V0RadiusXi = TMath::Sqrt(posV0Xi[0]*posV0Xi[0]
			     +posV0Xi[1]*posV0Xi[1]
			     +posV0Xi[2]*posV0Xi[2]);
    dcaPosToPrimaryVtxXi = TMath::Abs(pTrkXi->GetD(bestPrimaryVtxPos[0],
						   bestPrimaryVtxPos[1],
						   bestPrimaryVtxPos[2]));
    dcaNegToPrimaryVtxXi = TMath::Abs(nTrkXi->GetD(bestPrimaryVtxPos[0],
						   bestPrimaryVtxPos[1],
						   bestPrimaryVtxPos[2]));

    //apply cuts
    //if(XiRadius < 0.9 || XiRadius > 100.) continue;
    //if(dcaXiDaughters > 0.2) continue;
    //if(XiCosOfPointingAngle < 0.99) continue;
    //if(dcaV0ToPrimaryVtxXi < 0.03) continue;
    //if(dcaBachToPrimaryVtxXi < 0.01) continue;
    
    if(dcaXiDaughters > fCascadeCuts[0]) continue;
    if(XiCosOfPointingAngle < fCascadeCuts[1]) continue;
    if(dcaV0ToPrimaryVtxXi < fCascadeCuts[2]) continue;
    if(dcaBachToPrimaryVtxXi < fCascadeCuts[3]) continue;

    //V0 mass cut?
    //  if(TMath::Abs(invMassLambdaAsCascDghter-1.11568) > 0.01) continue;
    if(TMath::Abs(invMassLambdaAsCascDghter-1.11568) > fCascadeCuts[7]) 
      continue;
    
    //if(dcaV0DaughtersXi > 1.) continue;
    //if(V0CosOfPointingAngleXi > 0.9999) continue;
    //if(dcaPosToPrimaryVtxXi < 0.1) continue;
    //if(dcaNegToPrimaryVtxXi < 0.1) continue;

    if(dcaV0DaughtersXi > fCascadeCuts[4]) continue;
    if(V0CosOfPointingAngleXi > fCascadeCuts[5]) continue;
    if(dcaPosToPrimaryVtxXi < fCascadeCuts[6]) continue;
    if(dcaNegToPrimaryVtxXi < fCascadeCuts[6]) continue;

    //if(V0RadiusXi < 1.0 || V0RadiusXi > 100) continue;
    
    //other cuts?
    // change mass hypothesis to cover all the possibilities
    if(bTrkXi->Charge()<0){
      V0quality = 0.;
      xi->ChangeMassHypothesis(V0quality, 3312); //Xi- hyp.
      invMassXiMinus = xi->GetEffMassXi();

      V0quality = 0.;
      xi->ChangeMassHypothesis(V0quality, 3334); //Omega- hyp.
      invMassOmegaMinus = xi->GetEffMassXi();

      V0quality = 0.;
      xi->ChangeMassHypothesis(V0quality, 3312); //back to default hyp.
    }

    if(bTrkXi->Charge() > 0){
      V0quality = 0.;
      xi->ChangeMassHypothesis(V0quality, -3312); //anti-Xi- hyp.
      invMassXiPlus = xi->GetEffMassXi();

      V0quality = 0.;
      xi->ChangeMassHypothesis(V0quality, -3334); //anti-Omega- hyp.
      invMassOmegaPlus = xi->GetEffMassXi();

      V0quality = 0.;
      xi->ChangeMassHypothesis(V0quality, -3312); //back to default hyp.
    }

    //PID on the daughter tracks
    /*
    //A - Combined PID
    //Resonable guess the priors for the cascade track sample
    //(e, mu, pi, K, p)
    Double_t priorsGuessXi[5] = {0, 0, 2, 0, 1};
    Double_t priorsGuessOmega[5] = {0, 0, 1, 1, 1};

    //Combined bachelor-daughter PID
    AliPID pidXi;
    pidXi.SetPriors(priorsGuessXi);
    AliPID pidOmega;
    pidOmega.SetPriors(priorsGuessOmega);

    if(pTrkXi->IsOn(AliESDtrack::kESDpid)){// combined PID exists
      Double_t r[10] = {0.};
      pTrkXi->GetESDpid(r);
      pidXi.SetProbabilities(r);
      pidOmega.SetProbabilities(r);

      //Check if the V0 postive track is proton (case for Xi-)
      Double_t pProton = pidXi.GetProbability(AliPID::kProton);
      if(pProton > pidXi.GetProbability(AliPID::kElectron)
         && pProton > pidXi.GetProbability(AliPID::kMuon)
         && pProton > pidXi.GetProbability(AliPID::kPion)
         && pProton > pidXi.GetProbability(AliPID::kKaon))
        isPosInXiProton = kTRUE;
        
      //Check if the V0 postive track is a pi+ (case for Xi+)
      Double_t pPion = pidXi.GetProbability(AliPID::kPion);
      if(pPion > pidXi.GetProbability(AliPID::kElectron)
         && pPion > pidXi.GetProbability(AliPID::kMuon)
         && pPion > pidXi.GetProbability(AliPID::kKaon)
         && pPion > pidXi.GetProbability(AliPID::kProton))
        isPosInXiPion = kTRUE;
      // Check if the V0 positive track is a proton (case for Omega-)
      pProton = pidOmega.GetProbability(AliPID::kProton);
      if(pProton > pidOmega.GetProbability(AliPID::kElectron)
         && pProton > pidOmega.GetProbability(AliPID::kMuon)
         && pProton > pidOmega.GetProbability(AliPID::kPion)
         && pProton > pidOmega.GetProbability(AliPID::kKaon))
        isPosInOmegaProton = kTRUE;
    
      // Check if the V0 positive track is a pi+ (case for Omega+)
      pPion =  pidOmega.GetProbability(AliPID::kPion);
      if(pPion > pidOmega.GetProbability(AliPID::kElectron)
         && pPion > pidOmega.GetProbability(AliPID::kMuon)
         && pPion > pidOmega.GetProbability(AliPID::kKaon)
         && pPion > pidOmega.GetProbability(AliPID::kProton))
        isPosInOmegaPion = kTRUE;
    }

    //Combined V0-negative-daughter PID
    pidXi.SetPriors(priorsGuessXi);
    pidOmega.SetPriors(priorsGuessOmega);
    if(nTrkXi->IsOn(AliESDtrack::kESDpid)){
      Double_t r[10] = {0.};
      nTrkXi->GetESDpid(r);
      pidXi.SetProbabilities(r);
      pidOmega.SetProbabilities(r);
      
      // Check if the V0 negative track is a pi- (case for Xi-)
      Double_t pPion = pidXi.GetProbability(AliPID::kPion);
      if(pPion > pidXi.GetProbability(AliPID::kElectron)
         && pPion > pidXi.GetProbability(AliPID::kMuon)
         && pPion >  pidXi.GetProbability(AliPID::kKaon)
         && pPion > pidXi.GetProbability(AliPID::kProton))
        isNegInXiPion = kTRUE;

      // Check if the V0 negative track is an anti-proton (case for Xi+)
      Double_t pProton = pidXi.GetProbability(AliPID::kProton);
      if(pProton > pidXi.GetProbability(AliPID::kElectron)
         && pProton > pidXi.GetProbability(AliPID::kMuon)
         && pProton > pidXi.GetProbability(AliPID::kPion) 
         && pProton > pidXi.GetProbability(AliPID::kKaon))
        isNegInXiProton = kTRUE;
      
      // Check if the V0 negative track is a pi- (case for Omega-)
      pPion = pidOmega.GetProbability(AliPID::kPion);
      if(pPion > pidOmega.GetProbability(AliPID::kElectron)
         && pPion > pidOmega.GetProbability(AliPID::kMuon)
         && pPion > pidOmega.GetProbability(AliPID::kKaon)
         && pPion > pidOmega.GetProbability(AliPID::kProton))
        isNegInOmegaPion = kTRUE;
      
      // Check if the V0 negative track is an anti-proton (case for Omega+)   
      pProton =  pidOmega.GetProbability(AliPID::kProton);
      if(pProton > pidOmega.GetProbability(AliPID::kElectron)
         && pProton > pidOmega.GetProbability(AliPID::kMuon)
         && pProton > pidOmega.GetProbability(AliPID::kPion)
         && pProton > pidOmega.GetProbability(AliPID::kKaon))
        isNegInOmegaProton = kTRUE;
      
    }

    // Combined bachelor PID
    pidXi.SetPriors(priorsGuessXi);
    pidOmega.SetPriors(priorsGuessOmega);
    if(bTrkXi->IsOn(AliESDtrack::kESDpid)){//Combined PID exists
      Double_t r[10] = {0.};
      bTrkXi->GetESDpid(r);
      pidXi.SetProbabilities(r);
      pidOmega.SetProbabilities(r);
      
      //Check if the bachelor track is a pion
      Double_t pPion = pidXi.GetProbability(AliPID::kPion);
      if(pPion >  pidXi.GetProbability(AliPID::kElectron)
	 && pPion >  pidXi.GetProbability(AliPID::kMuon)
	 && pPion > pidXi.GetProbability(AliPID::kKaon)
	 && pPion > pidXi.GetProbability(AliPID::kProton))
	isBachelorPion = kTRUE;
      
      // Check if the bachelor track is a kaon
      Double_t pKaon = pidOmega.GetProbability(AliPID::kKaon);
      if(pKaon > pidOmega.GetProbability(AliPID::kElectron)
	 && pKaon > pidOmega.GetProbability(AliPID::kMuon)
	 && pKaon > pidOmega.GetProbability(AliPID::kPion)
	 && pKaon > pidOmega.GetProbability(AliPID::kProton))
	isBachelorKaon = kTRUE;
    }
    */


    //B - TPC PID: 3-sigma bands on Bethe-Bloch curve
    //Bachelor
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bTrkXi, AliPID::kKaon))<3.)
      isBachelorKaonForTPC = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bTrkXi, AliPID::kPion))<3.)
      isBachelorPionForTPC = kTRUE;

    //Negative V0 daughter
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrkXi, AliPID::kPion))<3.)
      isNegPionForTPC = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrkXi, AliPID::kProton))<3.)
      isNegProtonForTPC = kTRUE;
    
    //Positive V0 daughter
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrkXi, AliPID::kPion))<3.)
      isPosPionForTPC = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrkXi, AliPID::kProton))<3.)
      isPosProtonForTPC = kTRUE;
 
   
    //Extra QA information
    xi->GetPxPyPz(XiPx, XiPy, XiPz);
    XiPt = TMath::Sqrt(XiPx*XiPx + XiPy*XiPy);
    XiPtot= TMath::Sqrt(XiPx*XiPx + XiPy*XiPy + XiPz*XiPz);
    
    XiPAfter[0] = XiPx;
    XiPAfter[1] = XiPy;
    XiPAfter[2] = XiPz;

    xi->GetBPxPyPz(bachPx, bachPy, bachPz);
    bachPt = TMath::Sqrt(bachPx*bachPx + bachPy*bachPy);
    bachPtot = TMath::Sqrt(bachPx*bachPx + bachPy*bachPy + bachPz*bachPz);

    //chargeXi = xi->Charge();
    
    V0toXiCosOfPointingAngle 
      = xi->GetV0CosineOfPointingAngle(posXi[0], posXi[1], posXi[2]);
    rapXi = xi->RapXi();
    rapOmega = xi->RapOmega();
    phi = xi->Phi();
    alphaXi = xi->AlphaXi();
    ptArmXi = xi->PtArmXi();

    distToVtxZBefore = posXi[2]-bestPrimaryVtxPos[2];
    distToVtxXYBefore 
      = TMath::Sqrt((posXi[0] - bestPrimaryVtxPos[0])
		    *(posXi[0] - bestPrimaryVtxPos[0])
		    +(posXi[1] - bestPrimaryVtxPos[1])
                    *(posXi[1] - bestPrimaryVtxPos[1]));
    
    //propagation to the best primary vertex to determine the momentum
    Propagate(bestPrimaryVtxPos, posXi, XiPAfter, b, xi->Charge());
    distToVtxZAfter = posXi[2] - bestPrimaryVtxPos[2];
    distToVtxXYAfter = TMath::Sqrt((posXi[0] - bestPrimaryVtxPos[0])
				   *(posXi[0] - bestPrimaryVtxPos[0])
				   +(posXi[1] - bestPrimaryVtxPos[1])
				   *(posXi[1] - bestPrimaryVtxPos[1]));
    phiAfter = TMath::Pi() + TMath::ATan2(-XiPAfter[1],-XiPAfter[0]);
    
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DistToVtxZAfter"))->Fill(distToVtxZAfter);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DistToVtxXYAfter"))->Fill(distToVtxXYAfter);
    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DistToVtxZBeforeVsAfter"))->Fill(distToVtxZBefore, distToVtxZAfter);
    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DistToVtxXYBeforeVsAfter"))->Fill(distToVtxXYBefore, distToVtxXYAfter);
    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("PxBeforeVsAfter"))->Fill(XiPx, XiPAfter[0]);
    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("PyBeforeVsAfter"))->Fill(XiPy, XiPAfter[1]);
    if(xi->Charge()>0)
      ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("PhiPosBeforeVsAfter"))->Fill(phi, phiAfter);
    else if(xi->Charge()<0)
      ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("PhiNegBeforeVsAfter"))->Fill(phi, phiAfter);

    
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("Chi2Xi"))->Fill(chi2Xi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DcaXiDaughters"))->Fill(dcaXiDaughters);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DcaBachToPrimVertex"))->Fill(dcaBachToPrimaryVtxXi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("XiCosOfPointingAngle"))->Fill(XiCosOfPointingAngle);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("XiRadius"))->Fill(XiRadius);
    
    //V0
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("MassLambdaAsCascDghter"))->Fill(invMassLambdaAsCascDghter);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("V0Chi2Xi"))->Fill(V0Chi2Xi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DcaV0DaughtersXi"))->Fill(dcaV0DaughtersXi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("V0CosOfPointingAngleXi"))->Fill(V0CosOfPointingAngleXi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("V0RadiusXi"))->Fill(V0RadiusXi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DcaV0ToPrimVertexXi"))->Fill(dcaV0ToPrimaryVtxXi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DcaPosToPrimVertexXi"))->Fill(dcaPosToPrimaryVtxXi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DcaNegToPrimVertexXi"))->Fill(dcaNegToPrimaryVtxXi);
    
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("V0toXiCosOfPointingAngle"))->Fill(V0toXiCosOfPointingAngle);
    
    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("Armenteros"))->Fill(alphaXi, ptArmXi);
  
    //PID cuts with TPC cuts
    if(xi->Charge() < 0){
      if(isPosProtonForTPC
         && isNegPionForTPC){
	
	switch(fSpecie) {
	case 0:
	  if(isBachelorPionForTPC && TMath::Abs(rapXi) < 0.8){
	    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("MassVsPtAll"))->Fill(XiPt, invMassXiMinus);
	  
	  //candidate inserting
	    MakeTrack(invMassXiMinus, XiPt, /*xi->Phi()*/
		      phiAfter, xi->Eta(),  pTrkXi->GetID(),
		      nTrkXi->GetID(), bTrkXi->GetID());
	  }
	  break;
	  
	case 1:
	  if(isBachelorKaonForTPC && TMath::Abs(rapOmega) < 0.8){
	    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("MassVsPtAll"))->Fill(XiPt, invMassOmegaMinus);
	    MakeTrack(invMassOmegaMinus, XiPt, /*xi->Phi()*/
		      phiAfter, xi->Eta(),
		      pTrkXi->GetID(), nTrkXi->GetID(), bTrkXi->GetID());
	  }
	  break;
	}
      }
    }

    if(xi->Charge() > 0){
      if(isNegProtonForTPC
         && isPosPionForTPC){
	
	switch (fSpecie){
	case 0:
	  if(isBachelorPionForTPC && TMath::Abs(rapXi) < 0.8){
	    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("MassVsPtAll"))->Fill(XiPt, invMassXiPlus);
	  
	  //candidate inserting                                              
	    MakeTrack(invMassXiPlus, XiPt, /*xi->Phi()*/
		      phiAfter, xi->Eta(),
		      pTrkXi->GetID(), nTrkXi->GetID(), bTrkXi->GetID());
	  }
	  break;

	case 1:
	  if(isBachelorKaonForTPC && TMath::Abs(rapOmega) < 0.8){
	    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("MassVsPtAll"))->Fill(XiPt, invMassOmegaPlus);
	    MakeTrack(invMassOmegaPlus, XiPt, /*xi->Phi()*/
		      phiAfter, xi->Eta(),
		      pTrkXi->GetID(), nTrkXi->GetID(), bTrkXi->GetID());
	  }
	  break;
	}
      }
    }
  
  }

  return;
}

//______________________________________________________________________________
void AliAnalysisTaskFlowCascade::ReadFromAODv0(AliAODEvent *fAOD)
{
  
  fCutsRPTPC->SetEvent(fAOD, MCEvent());
  fCutsRPVZE->SetEvent(fAOD, MCEvent());
  fCutsPOI->SetEvent(fAOD, MCEvent());
  fFlowEventTPC->Fill(fCutsRPTPC, fCutsPOI);
  fFlowEventVZE->Fill(fCutsRPVZE, fCutsPOI);

  //  Double_t trkPrimaryVtxPos[3] = {-100., -100., -100.};
  Double_t bestPrimaryVtxPos[3] = {-100., -100., -100.};

  Double_t b = fAOD->GetMagneticField();

  int nCascades=fAOD->GetNumberOfCascades();
  const AliAODVertex *primaryBestAODVtx = fAOD->GetPrimaryVertex();
  primaryBestAODVtx->GetXYZ(bestPrimaryVtxPos);
  
  // calculation part dedicated to Xi vertices
  for(Int_t iXi = 0; iXi < nCascades; iXi++){
    Double_t effMassXi = 0.;
    Double_t chi2Xi = -1.;
    Double_t dcaXiDaughters = -1.;
    Double_t XiCosOfPointingAngle = -1.;
    Double_t posXi[3] = {-1000., -1000., -1000.};
    Double_t XiRadius = -1000.;
    
    Double_t invMassLambdaAsCascDghter = 0.;
    Double_t V0Chi2Xi = -1.;
    Double_t dcaV0DaughtersXi = -1.;
    
    Double_t dcaBachToPrimaryVtxXi = -1.;
    Double_t dcaV0ToPrimaryVtxXi = -1.;
    Double_t dcaPosToPrimaryVtxXi = -1.;
    Double_t dcaNegToPrimaryVtxXi = -1.;
    Double_t V0CosOfPointingAngleXi = -1.;
    Double_t posV0Xi[3] = {-1000., -1000., -1000.};
    Double_t V0RadiusXi = -1000.;
    //    Double_t V0quality = 0.;

    Double_t invMassXiMinus = 0.;
    Double_t invMassXiPlus = 0.;
    Double_t invMassOmegaMinus = 0.;
    Double_t invMassOmegaPlus = 0.;
    
    /*
    Bool_t isPosInXiProton = kFALSE;
    Bool_t isPosInXiPion = kFALSE;
    Bool_t isPosInOmegaProton = kFALSE;
    Bool_t isPosInOmegaPion = kFALSE;
    
    Bool_t isNegInXiProton = kFALSE;
    Bool_t isNegInXiPion = kFALSE;
    Bool_t isNegInOmegaProton = kFALSE;
    Bool_t isNegInOmegaPion = kFALSE;

    Bool_t isBachelorKaon = kFALSE;
    Bool_t isBachelorPion = kFALSE;
    */


    Bool_t isBachelorKaonForTPC = kFALSE;
    Bool_t isBachelorPionForTPC = kFALSE;
    Bool_t isNegPionForTPC = kFALSE;
    Bool_t isPosPionForTPC = kFALSE;
    Bool_t isNegProtonForTPC = kFALSE;
    Bool_t isPosProtonForTPC = kFALSE;
    
    Double_t XiPx = 0., XiPy = 0., XiPz = 0.;
    Double_t XiPt = 0.;
    Double_t XiPtot = 0.;
    
    Double_t bachPx = 0., bachPy = 0., bachPz = 0.;
    Double_t bachPt = 0.;
    Double_t bachPtot = 0.;
    
    //Short_t chargeXi = -2;
    Double_t V0toXiCosOfPointingAngle = 0.;
    
    Double_t rapXi = -20.;
    Double_t rapOmega = -20.;
    Double_t phi = 6.3;
    Double_t alphaXi = -200.;
    Double_t ptArmXi = -200.;

    Double_t distToVtxZBefore = -999.;
    Double_t distToVtxZAfter = -999.;
    Double_t distToVtxXYBefore = -999.;
    Double_t distToVtxXYAfter = -999.;
    Double_t XiPAfter[3] = {-999., -999., -999.};
    Double_t phiAfter = -999.;

    const AliAODcascade *xi = fAOD->GetCascade(iXi);
    if (!xi) continue;

    effMassXi = xi->MassXi(); //default working hypothesis: Xi- decay
    chi2Xi = xi->Chi2Xi();
    dcaXiDaughters = xi->DcaXiDaughters();
    XiCosOfPointingAngle = xi->CosPointingAngleXi(bestPrimaryVtxPos[0],
						  bestPrimaryVtxPos[1],
						  bestPrimaryVtxPos[2]);
    posXi[0] = xi->DecayVertexXiX();
    posXi[1] = xi->DecayVertexXiY();
    posXi[2] = xi->DecayVertexXiZ();
    XiRadius = TMath::Sqrt(posXi[0]*posXi[0]
                           +posXi[1]*posXi[1]
                           +posXi[2]*posXi[2]);

    AliAODTrack *pTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
    AliAODTrack *nTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
    AliAODTrack *bTrkXi 
      = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) );

    if(!pTrkXi || !nTrkXi || !bTrkXi) continue;

    UInt_t idxPosXi  = (UInt_t) TMath::Abs( pTrkXi->GetID() );
    UInt_t idxNegXi  = (UInt_t) TMath::Abs( nTrkXi->GetID() );
    UInt_t idxBach   = (UInt_t) TMath::Abs( bTrkXi->GetID() );

    if(idxBach == idxNegXi || idxBach == idxPosXi) continue;

    if( !fCutsDau->IsSelected(pTrkXi) 
        || !fCutsDau->IsSelected(nTrkXi)
        || !fCutsDau->IsSelected(bTrkXi) ) continue;

    
    if(pTrkXi->IsOn(AliESDtrack::kTPCin)){
      ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("TPCdEdxOfCascDghters"))->Fill(pTrkXi->P()*pTrkXi->Charge(), pTrkXi->GetTPCsignal());
    }
    if( nTrkXi->IsOn(AliESDtrack::kTPCin) ){
      ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("TPCdEdxOfCascDghters"))->Fill(nTrkXi->P()*nTrkXi->Charge(), nTrkXi->GetTPCsignal());
    }
    if(bTrkXi->IsOn(AliESDtrack::kTPCin)){
      ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("TPCdEdxOfCascDghters"))->Fill(bTrkXi->P()*bTrkXi->Charge(), bTrkXi->GetTPCsignal());
    }
    
    if(xi->ChargeXi() < 0)
      invMassLambdaAsCascDghter = xi->MassLambda();
    else
      invMassLambdaAsCascDghter = xi->MassAntiLambda();
    
    dcaV0DaughtersXi = xi->DcaV0Daughters();
    V0Chi2Xi = xi->Chi2V0();
    V0CosOfPointingAngleXi 
      = xi->CosPointingAngle(bestPrimaryVtxPos);
    dcaV0ToPrimaryVtxXi = xi->DcaV0ToPrimVertex();
    dcaBachToPrimaryVtxXi = xi->DcaBachToPrimVertex();
    
    //V0
    posV0Xi[0] = xi->DecayVertexV0X();
    posV0Xi[1] = xi->DecayVertexV0Y();
    posV0Xi[2] = xi->DecayVertexV0Z();
    V0RadiusXi = TMath::Sqrt(posV0Xi[0]*posV0Xi[0]
                             +posV0Xi[1]*posV0Xi[1]
                             +posV0Xi[2]*posV0Xi[2]);
    dcaPosToPrimaryVtxXi = xi->DcaPosToPrimVertex();
    dcaNegToPrimaryVtxXi = xi->DcaNegToPrimVertex();

    //apply cuts ?
    // if(XiRadius < 1. || XiRadius > 100.) continue;
    //if(dcaXiDaughters > 0.1) continue;
    //if(XiCosOfPointingAngle < 0.999) continue;
    //if(dcaV0ToPrimaryVtxXi < 0.05) continue;
    //if(dcaBachToPrimaryVtxXi < 0.03) continue;

    if(dcaXiDaughters > fCascadeCuts[0]) continue;
    if(XiCosOfPointingAngle < fCascadeCuts[1]) continue;
    if(dcaV0ToPrimaryVtxXi < fCascadeCuts[2]) continue;
    if(dcaBachToPrimaryVtxXi < fCascadeCuts[3]) continue;
    
    //V0 mass cut?
    //if(TMath::Abs(invMassLambdaAsCascDghter-1.11568) > 0.006) continue;
    if(TMath::Abs(invMassLambdaAsCascDghter-1.11568) > fCascadeCuts[7]) 
      continue;

    //if(dcaV0DaughtersXi > 1.) continue;
    //if(V0CosOfPointingAngleXi > 0.9999) continue;
    //if(dcaPosToPrimaryVtxXi < 0.1) continue;
    //if(dcaNegToPrimaryVtxXi < 0.1) continue;
    if(dcaV0DaughtersXi > fCascadeCuts[4]) continue;
    if(V0CosOfPointingAngleXi > fCascadeCuts[5]) continue;
    if(dcaPosToPrimaryVtxXi < fCascadeCuts[6]) continue;
    if(dcaNegToPrimaryVtxXi < fCascadeCuts[6]) continue;
    
    // if(V0RadiusXi < 1.0 || V0RadiusXi > 100) continue;
    
    //other cuts?


    //???
    if(xi->ChargeXi()<0){
      invMassXiMinus = xi->MassXi();
      invMassOmegaMinus = xi->MassOmega();
    }else{
      invMassXiPlus = xi->MassXi();
      invMassOmegaPlus = xi->MassOmega();
    }

    /*
    if(pTrkXi->GetMostProbablePID() == AliAODTrack::kProton) {
      isPosInXiProton = kTRUE;
      isPosInOmegaProton = kTRUE;
    }
    if(pTrkXi->GetMostProbablePID() == AliAODTrack::kPion){
      isPosInXiPion = kTRUE;
      isPosInOmegaPion = kTRUE;
    }
    
    if(nTrkXi->GetMostProbablePID() == AliAODTrack::kPion){
      isNegInXiPion = kTRUE;
      isNegInOmegaPion = kTRUE;
    }
    if(nTrkXi->GetMostProbablePID() == AliAODTrack::kProton){
      isNegInXiProton = kTRUE;
      isNegInOmegaProton = kTRUE;
    }

    if(bTrkXi->GetMostProbablePID() == AliAODTrack::kPion)
      isBachelorPion = kTRUE;
    if(bTrkXi->GetMostProbablePID() == AliAODTrack::kKaon)
      isBachelorKaon = kTRUE;
    */

    //PID with TPC only: 
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bTrkXi, AliPID::kKaon))<3.)
      isBachelorKaonForTPC = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bTrkXi, AliPID::kPion))<3.)
      isBachelorPionForTPC = kTRUE;

    //Negative V0 daughter
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrkXi, AliPID::kPion))<3.)
      isNegPionForTPC = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrkXi, AliPID::kProton))<3.)
      isNegProtonForTPC = kTRUE;
    
    //Positive V0 daughter
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrkXi, AliPID::kPion))<3.)
      isPosPionForTPC = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrkXi, AliPID::kProton))<3.)
      isPosProtonForTPC = kTRUE;

    //Extra QA information
    XiPx = xi->MomXiX();
    XiPy = xi->MomXiY();
    XiPz = xi->MomXiZ();
    XiPt = TMath::Sqrt(XiPx*XiPx + XiPy*XiPy);
    XiPtot= TMath::Sqrt(XiPx*XiPx + XiPy*XiPy + XiPz*XiPz);
    
    bachPx = xi->MomBachX();
    bachPy = xi->MomBachY();
    bachPz = xi->MomBachZ();
    
    bachPt = TMath::Sqrt(bachPx*bachPx + bachPy*bachPy);
    bachPtot = TMath::Sqrt(bachPx*bachPx + bachPy*bachPy + bachPz*bachPz);
  
    V0toXiCosOfPointingAngle = xi->CosPointingAngle( xi->GetDecayVertexXi() );
    
    rapXi = xi->RapXi();
    rapOmega = xi->RapOmega();
    phi = xi->Phi();
    alphaXi = xi->AlphaXi();
    ptArmXi = xi->PtArmXi();

    distToVtxZBefore = posXi[2]-bestPrimaryVtxPos[2];
    distToVtxXYBefore
      = TMath::Sqrt((posXi[0] - bestPrimaryVtxPos[0])
                    *(posXi[0] - bestPrimaryVtxPos[0])
                    +(posXi[1] - bestPrimaryVtxPos[1])
                    *(posXi[1] - bestPrimaryVtxPos[1]));


    XiPAfter[0] = XiPx;
    XiPAfter[1] = XiPy;
    XiPAfter[2] = XiPz;
    //propagation to the best primary vertex to determine the momentum
    Propagate(bestPrimaryVtxPos, posXi, XiPAfter, b, xi->ChargeXi());
    distToVtxZAfter = posXi[2] - bestPrimaryVtxPos[2];
    distToVtxXYAfter = TMath::Sqrt((posXi[0] - bestPrimaryVtxPos[0])
                                   *(posXi[0] - bestPrimaryVtxPos[0])
                                   +(posXi[1] - bestPrimaryVtxPos[1])
                                   *(posXi[1] - bestPrimaryVtxPos[1]));
    phiAfter = TMath::Pi() + TMath::ATan2(-XiPAfter[1],-XiPAfter[0]);

    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DistToVtxZAfter"))->Fill(distToVtxZAfter);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DistToVtxXYAfter"))->Fill(distToVtxXYAfter);
    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DistToVtxZBeforeVsAfter"))->Fill(distToVtxZBefore, distToVtxZAfter);
    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DistToVtxXYBeforeVsAfter"))->Fill(distToVtxXYBefore, distToVtxXYAfter);
    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("PxBeforeVsAfter"))->Fill(XiPx, XiPAfter[0]);
    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("PyBeforeVsAfter"))->Fill(XiPy, XiPAfter[1]);
    if(xi->ChargeXi()>0)
      ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("PhiPosBeforeVsAfter"))->Fill(phi, phiAfter);
    else if(xi->ChargeXi()<0)
      ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("PhiNegBeforeVsAfter"))->Fill(phi, phiAfter);
    
    //for default hypothesis
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("Chi2Xi"))->Fill(chi2Xi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DcaXiDaughters"))->Fill(dcaXiDaughters);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DcaBachToPrimVertex"))->Fill(dcaBachToPrimaryVtxXi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("XiCosOfPointingAngle"))->Fill(XiCosOfPointingAngle);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("XiRadius"))->Fill(XiRadius);
    
    //V0
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("MassLambdaAsCascDghter"))->Fill(invMassLambdaAsCascDghter);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("V0Chi2Xi"))->Fill(V0Chi2Xi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DcaV0DaughtersXi"))->Fill(dcaV0DaughtersXi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("V0CosOfPointingAngleXi"))->Fill(V0CosOfPointingAngleXi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("V0RadiusXi"))->Fill(V0RadiusXi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DcaV0ToPrimVertexXi"))->Fill(dcaV0ToPrimaryVtxXi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DcaPosToPrimVertexXi"))->Fill(dcaPosToPrimaryVtxXi);
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("DcaNegToPrimVertexXi"))->Fill(dcaNegToPrimaryVtxXi);
    
    ((TH1F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("V0toXiCosOfPointingAngle"))->Fill(V0toXiCosOfPointingAngle);
    
    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("Armenteros"))->Fill(alphaXi, ptArmXi);
  
    //with PID cuts
    if(xi->ChargeXi()<0){
      if(isPosProtonForTPC && isNegPionForTPC){
	switch (fSpecie){
	case 0:
	  if( isBachelorPionForTPC && TMath::Abs(rapXi) < 0.8){
	    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("MassVsPtAll"))->Fill(XiPt, invMassXiMinus);
	    MakeTrack(invMassXiMinus, XiPt, /*xi->Phi(),*/
		      phiAfter, xi->Eta(), 
		      pTrkXi->GetID(), nTrkXi->GetID(), bTrkXi->GetID());
	  }// endif
      
	  break;

	case 1:
	  if(isBachelorKaonForTPC && TMath::Abs(rapOmega) < 0.8
	     && (invMassXiMinus > 1.32486 || invMassXiMinus < 1.30486)){
	    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("MassVsPtAll"))->Fill(XiPt, invMassOmegaMinus);      
	    
	    MakeTrack(invMassOmegaMinus, XiPt, /* xi->Phi(),*/
		      phiAfter,   xi->Eta(),
		      pTrkXi->GetID(), nTrkXi->GetID(), bTrkXi->GetID());
	  }//endif
	  break;
	}
      }
    }//end if ChargeXi()<0
    
    if(xi->ChargeXi() > 0){
      if(isNegProtonForTPC && isPosPionForTPC){ 
	switch(fSpecie){
	case 0:
	  if (isBachelorPionForTPC  && TMath::Abs(rapXi) < 0.8){
	    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("MassVsPtAll"))->Fill(XiPt, invMassXiPlus);
	  
	    //candidate inserting                                              
	    MakeTrack(invMassXiPlus, XiPt, /* xi->Phi(),*/
		      phiAfter, xi->Eta(),
		      pTrkXi->GetID(), nTrkXi->GetID(), bTrkXi->GetID());
	  }//endif particle id
	  break;
	  
	case 1:
	  if(isBachelorKaonForTPC && TMath::Abs(rapOmega) < 0.8
	     && (invMassXiPlus > 1.32486 || invMassXiPlus < 1.30486)){
	    ((TH2F*)((TList*)fQAList->FindObject("Candidates"))->FindObject("MassVsPtAll"))->Fill(XiPt, invMassOmegaPlus);
	    MakeTrack(invMassOmegaPlus, XiPt, /* xi->Phi(),*/
		      phiAfter, xi->Eta(),
		      pTrkXi->GetID(), nTrkXi->GetID(), bTrkXi->GetID());
	  }//endif particle id
	}
      }
    }//endif ChargeXi()>0  
  }//for Xi candidate loop

  return;

}


void AliAnalysisTaskFlowCascade::MakeTrack( double mass, 
					    double pt, 
					    double phi, 
					    double eta, 
					    int iid, 
					    int jid,
					    int kid) {
  // create track for flow tasks        
  if(fCandidates->GetLast()+1>=fCandidates->GetSize()) {
    fCandidates->Expand( 2*fCandidates->GetSize() );
  }
  Bool_t overwrite = kTRUE;

  AliFlowCandidateTrack *sTrack 
    = (static_cast<AliFlowCandidateTrack*> (fCandidates->At( fCandidates->GetLast()+1 )));
  if( !sTrack ) { // creates new
    sTrack = new AliFlowCandidateTrack();
    overwrite = kFALSE;
  } else { // overwrites
    sTrack->ClearMe();
  }


  sTrack->SetMass(mass);
  sTrack->SetPt(pt);
  sTrack->SetPhi(phi);
  sTrack->SetEta(eta);
  sTrack->AddDaughter(iid);
  sTrack->AddDaughter(jid);
  sTrack->AddDaughter(kid);
  sTrack->SetForPOISelection(kTRUE);
  sTrack->SetForRPSelection(kFALSE);
 
  if(overwrite) {
    fCandidates->SetLast( fCandidates->GetLast()+1 );
  } else {
    fCandidates->AddLast(sTrack);
  }
  
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowCascade::Terminate(Option_t *)
{

}

void AliAnalysisTaskFlowCascade::Propagate(Double_t vv[3], 
					   Double_t x[3], 
					   Double_t p[3], 
					   Double_t bz, 
					   Short_t sign){
  //Propagation to the primary vertex to determine the px and py
  //x, p are the position and momentum as input and output
  //bz is the magnetic field along z direction
  //sign is the charge of particle for propagation

  Double_t pp = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  Double_t len = (vv[2]-x[2])*pp/p[2];
  Double_t a = -kB2C*bz*sign;  

  Double_t rho = a/pp;
  x[0] += p[0]*TMath::Sin(rho*len)/a - p[1]*(1-TMath::Cos(rho*len))/a;
  x[1] += p[1]*TMath::Sin(rho*len)/a + p[0]*(1-TMath::Cos(rho*len))/a;
  x[2] += p[2]*len/pp;

  Double_t p0=p[0];
  p[0] = p0  *TMath::Cos(rho*len) - p[1]*TMath::Sin(rho*len);
  p[1] = p[1]*TMath::Cos(rho*len) + p0  *TMath::Sin(rho*len);
}


//=====================================================================       
void AliAnalysisTaskFlowCascade::SetCommonConstants(Int_t massBins, 
						    Double_t minMass, 
						    Double_t maxMass)
{
  // setter for mass bins                          
                          
  fMassBins = massBins;
  fMinMass = minMass;
  fMaxMass = maxMass;
}


//====================================================================         
void AliAnalysisTaskFlowCascade::SetCuts2010(int set) {

  // fCascadeCuts[0]: DcaXiDaughter; fCascadeCuts[1]: XiCosOfPointingAngle
  // fCascadeCuts[2]: DcaV0ToPrimaryVtxXi; fCascadeCuts[3]: DcaBachToPrimaryVtxXi
  // fCascadeCuts[4]: DcaV0DaughersXi; fCascadeCuts[5]: V0CosOfPointingAngleXi
  // fCascadeCuts[6]: DcaV0DaughterToPrimaryVtxXi; fCascadeCuts[7]: V0MassWidth
  
  switch(set){
  
  case 0: //tighter
    fCascadeCuts[0] = 0.2; fCascadeCuts[1] = 0.999;
    fCascadeCuts[2] = 0.03; fCascadeCuts[3] = 0.05;
    fCascadeCuts[4] = .5; fCascadeCuts[5] = 0.9998;
    fCascadeCuts[6] = 0.1; fCascadeCuts[7] = 0.006;
    break;
    
  case 1: //middle
    fCascadeCuts[0] = 0.3; fCascadeCuts[1] = 0.99;
    fCascadeCuts[2] = 0.01; fCascadeCuts[3] = 0.03;
    fCascadeCuts[4] = .6; fCascadeCuts[5] = 0.9999;
    fCascadeCuts[6] = 0.1; fCascadeCuts[7] = 0.008;
    break;

  case 2: //looser
    fCascadeCuts[0] = 0.3; fCascadeCuts[1] = 0.99;
    fCascadeCuts[2] = 0.01; fCascadeCuts[3] = 0.03;
    fCascadeCuts[4] = 1.; fCascadeCuts[5] = 1.;
    fCascadeCuts[6] = 0.1; fCascadeCuts[7] = 0.01;
    break;
  }
  
}
