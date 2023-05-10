#include <Riostream.h>
#include "TTree.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDatabasePDG.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliVEventHandler.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"

#include "AliAnalysisTaskfn.h"
#include "AliPIDResponse.h"
#include "AliPhysicsSelection.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"
#include "AliPPVsMultUtils.h"
#include "AliAnalysisUtils.h"

//_____ Event pool includes
#include "AliEventPoolManager.h"
//#include "AliPool.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskfn)
ClassImp(AliCompactTrack)

AliAnalysisTaskfn::AliAnalysisTaskfn(): AliAnalysisTaskSE(), fOutput(0), fPoolMgr(0X0), fpionreduced(0), fPIDResponse(0), fVevent(0), lESDevent(0), fEventCuts(0), fESDtrackCuts(0), fHistZVertex(0),fHistCentralityEvtCount(0), fHisteventsummary(0),f1Unlike(0),f1Like(0),f1Mix(0)
{

}

AliAnalysisTaskfn::AliAnalysisTaskfn(const char *name): AliAnalysisTaskSE(name), fOutput(0), fPoolMgr(0X0), fpionreduced(0), fPIDResponse(0), fVevent(0), lESDevent(0), fEventCuts(0), fESDtrackCuts(0), fHistZVertex(0),fHistCentralityEvtCount(0), fHisteventsummary(0),f1Unlike(0),f1Like(0),f1Mix(0)
{
  DefineInput(0, TChain::Class()); 
  DefineOutput(1, TList::Class()); 
}


AliAnalysisTaskfn::~AliAnalysisTaskfn()
{
  //------------------------------------------------
  // DESTRUCTOR
  //------------------------------------------------
    
  if (fOutput){
    delete fOutput;
    fOutput = 0x0;
  }
  
  if(fVevent) {
    delete fVevent;
    fVevent=0x0;
  }

  if(lESDevent) {
    delete lESDevent;
    lESDevent=0x0;
  }

  if(fESDtrackCuts) {
    delete fESDtrackCuts;
    fESDtrackCuts=0x0;
  }


  if (fPoolMgr) { delete fPoolMgr; fPoolMgr = 0x0; }
  if (fPIDResponse) { delete fPIDResponse; fPIDResponse = 0x0; }
}

//________________________________________________________________________
void AliAnalysisTaskfn::UserCreateOutputObjects()
{
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  if(!fESDtrackCuts ){
    fESDtrackCuts = new AliESDtrackCuts();
    fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1); 
    fESDtrackCuts->SetPtRange(0.15,20.0);
    fESDtrackCuts->SetEtaRange(-0.8,0.8);
  }


  SetupForMixing();

  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!
  OpenFile(1);
  fHisteventsummary            = new TH1F("fHisteventsummary", "Event cut summary", 5,0.0,5.0);       
  fHistZVertex            = new TH1F("fHistZVertex", "Z vertex distribution", 100,-10,10);       
  fHistCentralityEvtCount = new TH1F("fHistCentralityEvtCount", "Centrality Event Count", 100,0,100);

  Int_t bins[2]={250, 200};
  Double_t xmin[2]={1.0, 0.0};
  Double_t xmax[2]={2.0, 20.0};
  f1Unlike = new THnSparseD("f1Unlike", "Unlike histogram", 2, bins, xmin, xmax);
  f1Like = new THnSparseD("f1Like", "Like histogram", 2, bins, xmin, xmax);
  f1Mix = new THnSparseD("f1Mix", "Mix histogram", 2, bins, xmin, xmax);


  f1Unlike->Sumw2();
  f1Like->Sumw2();
  f1Mix->Sumw2();


  fOutput->Add(fHisteventsummary);
  fOutput->Add(fHistZVertex);
  fOutput->Add(fHistCentralityEvtCount);
  fOutput->Add(f1Unlike);
  fOutput->Add(f1Like);
  fOutput->Add(f1Mix);
  PostData(1, fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskfn::UserExec(Option_t *)
{
     
  // Main loop
  // Called for each event
  //AliAODEvent *lAODevent = 0x0;
  
  lESDevent = dynamic_cast <AliESDEvent*> (InputEvent());

  if (!(lESDevent)) {
    AliWarning("ERROR: event not available \n");
    PostData(1, fOutput);
    return;
  }
  
  ////trigger/////////////
  /*
  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isSelected = 0;
  isSelected = (maskIsSelected & AliVEvent::kINT7);
  //isSelected = (maskIsSelected & (AliVEvent::kHighMultV0 | AliVEvent::kINT7));
  */
  
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t maskIsSelected = inputHandler->IsEventSelected();
  Bool_t isSelected = kFALSE;
  isSelected = maskIsSelected& AliVEvent::kINT7;;



  if ( ! isSelected)
    {
      
      PostData(1, fOutput);
      return;
    }

  fHisteventsummary->Fill(0.5);

  //AliVEvent *fVevent=0x0;
  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fVevent){
    printf("ERROR: fVevent not available\n");
    return;
  }


  //event selection  
  const AliVVertex *vertex = fVevent->GetPrimaryVertex();
  if (!GoodEvent(vertex)) return; 
  double zv=vertex->GetZ();
  //double yv=vertex->GetY();
  //double xv=vertex->GetX();
  fHistZVertex->Fill(zv);
  ////////////////////////////////


  //////centrality selection/////////
  Float_t lV0M;
  Int_t lEvSelCode = 300;
  AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
  if( !MultSelection)
    {
      AliWarning("AliMultSelection object not found!");
      
      PostData(1, fOutput);
      return;
    }
  else
    {
      lV0M = MultSelection->GetMultiplicityPercentile("V0M"); 
    }
    
  

  fHistCentralityEvtCount->Fill(lV0M);
  fHisteventsummary->Fill(1.5);
  
  std::vector<Alikks0Container> kks0Candidates;
  std::vector<AlipionContainer> pionCandidates;
  Alikks0Container kks0;
  AlipionContainer pion;
    
  TLorentzVector daughterk(0.0,0.0,0.0,0.0);
  TLorentzVector daughterks0(0.0,0.0,0.0,0.0);
  TLorentzVector motherkks0(0.0,0.0,0.0,0.0);
  TLorentzVector mother(0.0,0.0,0.0,0.0);
  TLorentzVector mothermix(0.0,0.0,0.0,0.0);
    


  Int_t ntracks = lESDevent->GetNumberOfTracks();
  Int_t nv0s = lESDevent->GetNumberOfV0s();
    


  TObjArray* piontracksAccepted= new TObjArray;
  piontracksAccepted -> SetOwner(kTRUE);

  Int_t nproton=0;
  Int_t nkks0pair=0;
  Int_t npion=0;
  int mult=0;
  Double_t tV0mom[3] = {0.0, 0.0, 0.0};
  Double_t trkPt=0.0;

  for(Int_t itr = 0; itr < ntracks; itr++)
    {
      AliVTrack   *track = (AliVTrack*)lESDevent->GetTrack(itr);
      if(!track)      continue;
      AliESDtrack *esdtrack  = dynamic_cast<AliESDtrack*>(track);
      if(!esdtrack)      continue;
      if(!fESDtrackCuts->AcceptTrack(esdtrack))  continue;
      mult=mult+1;

      if(IsPion(track))
	{
	  piontracksAccepted->Add(new AliCompactTrack(track->Px(),track->Py(),track->Pz(),esdtrack->Charge()));
	  pion.particle.SetXYZM(track->Px(), track->Py(), track->Pz(), 0.1396);
	  pion.charge=esdtrack->Charge();
	  pion.tracknumber=itr;
	  pionCandidates.push_back(pion);
	  npion=npion+1;

	  
	}
 

      if(IsKaon(track)){
	daughterk.SetXYZM(track->Px(), track->Py(), track->Pz(), 0.4937);

	for(Int_t itrv0 = 0; itrv0 < nv0s; itrv0++)
	  {
	    AliESDv0 *v0=lESDevent->GetV0(itrv0);
	    if(!v0) continue;
	    v0->GetPxPyPz(tV0mom[0],tV0mom[1],tV0mom[2] );
	    if(!CheckESDV0(v0, lESDevent))continue;
	    if(v0->GetEffMass()<0.48 || v0->GetEffMass()>0.51) continue;
            v0->ChangeMassHypothesis(310);
	    daughterks0.SetXYZM(tV0mom[0], tV0mom[1], tV0mom[2], 0.4976);
	    motherkks0=daughterk+daughterks0;
	    kks0.particle.SetXYZM(motherkks0.Px(),motherkks0.Py(),motherkks0.Pz(),motherkks0.M());
	    kks0.charge=esdtrack->Charge();
	    kks0.tracknumber=itr;
	    kks0Candidates.push_back(kks0);
	    nkks0pair=nkks0pair+1;

	  }
      }
    }

  if(pionCandidates.size()<1 || kks0Candidates.size()<1) 
    {
      PostData(1, fOutput);
      return;
    }

  fHisteventsummary->Fill(2.5);


  for(int j=0;j<kks0Candidates.size();j++)
    {

      kks0=kks0Candidates[j];
      if (kks0.particle.M() > 1.04) continue;
      /////////////make same event pair///////////////////////
      for(int i=0;i<pionCandidates.size();i++)
	{
	  pion=pionCandidates[i];
	  if(pion.tracknumber==kks0.tracknumber) continue;
	  mother.SetPxPyPzE(pion.particle.Px()+kks0.particle.Px(), pion.particle.Py()+kks0.particle.Py(), pion.particle.Pz()+kks0.particle.Pz(), pion.particle.E()+kks0.particle.E());
	  if(mother.M()>2.0) continue;	  
	  if(pion.charge*kks0.charge<0)
	    {
	      /////////fill unlike histo//////////
	      f1Unlike->Fill(mother.M(),mother.Pt());
	      
	    }
	  if(pion.charge*kks0.charge>0)
	    {
	      /////////fill like histo//////////
	      f1Like->Fill(mother.M(),mother.Pt());
	    }
	}

    }
 fHisteventsummary->Fill(3.5);

  ///////////////////////////////////////////////////////////////////
      
  //////////////mixed event///////////////////
  double pionmixenergy=0.0;
  AliEventPool* pool = fPoolMgr->GetEventPool(mult,zv); 
  if (!pool) AliInfo("No pool found for centrality ");
  else
    {
      if (pool->IsReady())
	{
	  //Int_t nMix = pool->GetCurrentNEvents();
	  Int_t nMix=10;
	  for (Int_t jMix=0; jMix<nMix; jMix++)
	    {
	      TObjArray* bgTracks = pool->GetRandomEvent();
	      Int_t numTracks = bgTracks->GetEntriesFast();

	      for(int ipion = 0; ipion < numTracks; ipion++)
		{
		  AliCompactTrack *piontrackmix = (AliCompactTrack*) bgTracks->At(ipion);
		  if(!piontrackmix)
		    {
		      AliFatal(Form("ERROR: Could not receive mix pool track %d\n",ipion));
		      continue;
		    }
		  AliESDtrack *esdtrackmix  = dynamic_cast<AliESDtrack*>(piontrackmix);
		  for(int j=0;j<kks0Candidates.size();j++)
		    {
		      kks0=kks0Candidates[j];
		      if (kks0.particle.M() > 1.04) continue;
		      if(piontrackmix->Charge()*kks0.charge>0)continue;
		      pionmixenergy=TMath::Sqrt(0.139*0.139+piontrackmix->Px()*piontrackmix->Px()+piontrackmix->Py()*piontrackmix->Py()+piontrackmix->Pz()*piontrackmix->Pz());
		      mothermix.SetPxPyPzE(piontrackmix->Px()+kks0.particle.Px(), piontrackmix->Py()+kks0.particle.Py(), piontrackmix->Pz()+kks0.particle.Pz(), pionmixenergy+kks0.particle.E());
		      if(mothermix.M()>2.0) continue;
		      //////////fill mix histo//////////
		      f1Mix->Fill(mothermix.M(),mothermix.Pt());
		    }
		}

	    }
	}
      pool->UpdatePool(piontracksAccepted);
    }


  /////////////////////////////////////


    
  PostData(1, fOutput);
}
//________________________________________________________________________
void AliAnalysisTaskfn::Terminate(Option_t *)
{
  
}

//________________________________________________________________________

Bool_t AliAnalysisTaskfn::GoodEvent(const AliVVertex *vertex) //all cuts taken from alice code github pp analysis from strangeness and correlation analysis
{

  Bool_t EventAccepted;
  EventAccepted = fEventCuts.AcceptEvent(fVevent);
  if (!EventAccepted)
    {
      PostData(1, fOutput);
      return kFALSE;
    }   


  if (!vertex)
    {
      PostData(1, fOutput);
      return kFALSE;
    }  

  if (vertex->GetNContributors() < 1)
    {
      PostData(1, fOutput);
      return kFALSE;
    }
  
  if (TMath::Abs(vertex->GetZ()) > 10.0)
    {
      PostData(1, fOutput);
      return kFALSE;
    }

  
  
  AliPPVsMultUtils *multUtils = new AliPPVsMultUtils();
  Bool_t spdpileup=multUtils->IsNotPileupSPDInMultBins(fVevent);
  Bool_t inel0=multUtils->IsINELgtZERO(fVevent);
  Bool_t inconstspdvertex=multUtils->HasNoInconsistentSPDandTrackVertices(fVevent);
  
  if(!spdpileup) 
    {
      PostData(1, fOutput);
      return kFALSE;
    }
  if(!inel0)
    {
      PostData(1, fOutput);
      return kFALSE;
    }
  if(!inconstspdvertex) 
    {
      PostData(1, fOutput);
      return kFALSE;
    }
  if( lESDevent->IsIncompleteDAQ() ) 
    {
      PostData(1, fOutput);
      return kFALSE;
    }

    AliAnalysisUtils *fUtils = new AliAnalysisUtils();
    if(fUtils->IsSPDClusterVsTrackletBG(fVevent))
    {
      PostData(1, fOutput);
      return kFALSE;
    }


  return kTRUE;

}
//-----------------------------------------------
Bool_t AliAnalysisTaskfn::IsPion(AliVTrack *vtrack)
{

  
  AliESDtrack *esdtrack  = dynamic_cast<AliESDtrack*>(vtrack);
  Double_t nsigmatpcpion=TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kPion));
  Double_t nsigmatofpion=TMath::Abs(fPIDResponse->NumberOfSigmasTOF(esdtrack, AliPID::kPion));
  Bool_t TOFHIT=kFALSE;
  
  
  if ((vtrack->GetStatus() & AliESDtrack::kTOFout)) TOFHIT=kTRUE;
  if ((vtrack->GetStatus() & AliESDtrack::kTIME  )) TOFHIT=kTRUE;

  if(!TOFHIT)
    {
      if(nsigmatpcpion>3.0)return kFALSE;
    }

  if(TOFHIT)
    {
      //if(l<350)return kFALSE;                                                                                                               
      if(nsigmatofpion>3.0)return kFALSE;
    }



  return kTRUE;
}

Bool_t AliAnalysisTaskfn::IsKaon(AliVTrack *vtrack)
{

  AliESDtrack *esdtrack  = dynamic_cast<AliESDtrack*>(vtrack);
  Double_t nsigmatpckaon=TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kKaon));
  Double_t nsigmatofkaon=TMath::Abs(fPIDResponse->NumberOfSigmasTOF(esdtrack, AliPID::kKaon));
  Bool_t TOFHIT=kFALSE;

  
  if ((vtrack->GetStatus() & AliESDtrack::kTOFout)) TOFHIT=kTRUE;
  if ((vtrack->GetStatus() & AliESDtrack::kTIME  )) TOFHIT=kTRUE;

  if(!TOFHIT)
    {
      if(nsigmatpckaon>3.0)return kFALSE;
    }

  if(TOFHIT)
    {
      //if(l<350)return kFALSE;                                                                                                               
      if(nsigmatofkaon>3.0)return kFALSE;
    }


    
  return kTRUE;
}



//_________________________________________________________________________________________________
Bool_t AliAnalysisTaskfn::CheckESDV0(AliESDv0 *v0, AliESDEvent *lESDEvent)
{


  if (v0->GetOnFlyStatus()) {
    return kFALSE; 
  }
  
 
  Double_t xPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetX();
  Double_t yPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetY();
  Double_t zPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetZ();


  // retrieve the V0 daughters                                                                                                               
  UInt_t lIdxPos      = (UInt_t) TMath::Abs(v0->GetPindex());
  UInt_t lIdxNeg      = (UInt_t) TMath::Abs(v0->GetNindex());
  AliESDtrack *pTrack = lESDEvent->GetTrack(lIdxPos);
  AliESDtrack *nTrack = lESDEvent->GetTrack(lIdxNeg);

 

  // filter like-sign V0
  if ( TMath::Abs(((pTrack->GetSign()) - (nTrack->GetSign())) ) < 0.1) {
    return kFALSE;
  }
  

 
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterK0s");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.15,30);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0); //                                                                                                 
  esdTrackCuts->SetMinNCrossedRowsTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetMinDCAToVertexXY(0.06);
  esdTrackCuts->SetMinDCAToVertexZ(0.06);

  if (!esdTrackCuts->IsSelected(pTrack)) {
    return kFALSE;
  }

  if (!esdTrackCuts->IsSelected(nTrack)) {
    return kFALSE;
  }
  


  // topological checks
  if (TMath::Abs(v0->GetDcaV0Daughters()) > 1.0) {
    return kFALSE;
  }

  if (TMath::Abs(v0->GetD(xPrimaryVertex, yPrimaryVertex, zPrimaryVertex)) > 0.3) {
    return kFALSE;
  }

 
  Bool_t fCheckOOBPileup=kTRUE;
  if (fCheckOOBPileup) {
    Double_t bfield = lESDEvent->GetMagneticField();
    if(!TrackPassesOOBPileupCut(pTrack, bfield) &&
       !TrackPassesOOBPileupCut(nTrack, bfield)) return kFALSE;
  }
   
  if ( (TMath::Abs(v0->GetV0CosineOfPointingAngle()) < 0.97) || (TMath::Abs(v0->GetV0CosineOfPointingAngle()) >= 1 ) ) {
    return kFALSE;
  }

  
  if (TMath::Abs(v0->Eta())> 0.8) {
    return kFALSE;
  }

  if (TMath::Abs(v0->RapK0Short())>0.5){
    return kFALSE;
  }  

  Double_t v0Position[3];
  v0->GetXYZ(v0Position[0],v0Position[1],v0Position[2]);
  Double_t radius = TMath::Sqrt(TMath::Power(v0Position[0],2) + TMath::Power(v0Position[1],2));
  if ( ( radius < 0.5 ) || ( radius > 100 ) ) {
    return kFALSE;
  }

  
  Double_t tV0mom[3];
  v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] );
  Double_t lV0TotalMomentum =  TMath::Sqrt(tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );
  Double_t fLength = TMath::Sqrt(TMath::Power(v0Position[0]- xPrimaryVertex,2) + TMath::Power(v0Position[1] - yPrimaryVertex,2)+ TMath::Power(v0Position[2]- zPrimaryVertex,2));
  if( TMath::Abs(0.497*fLength/lV0TotalMomentum) > 15)
    {
      AliDebugClass(2, "Failed Lifetime Cut on positive track V0");
      return kFALSE;
    }


  
  // Competing v0rejection for K0s vs Lambda
  Double_t altmass=0.0;
  Double_t fToleranceVeto=0.004;
  //Double_t fToleranceVeto=0.01;
  
  v0->ChangeMassHypothesis(kLambda0);
  if ((TMath::Abs(v0->GetEffMass() - 1.115683)) < fToleranceVeto) {
    AliDebugClass(2, "Failed competing V0 rejection check");
    return kFALSE;
  }

  v0->ChangeMassHypothesis(kK0Short);

  // check PID
  Double_t posnsTPC   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion));
  Double_t negnsTPC   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion));
  //if(((negnsTPC > 4) && (posnsTPC > 4))) {
  if(! ((negnsTPC <= 4) && (posnsTPC <= 4)) ) {  
  return kFALSE;
  }

  return kTRUE;
  

}

//_________________________________________________________________________________________________
void AliAnalysisTaskfn::SetupForMixing()
{

  ////////////////////////////
  // Set-up for Mixed Event //
  ////////////////////////////
  const Int_t trackDepth = 10000;
  const Int_t poolsize = 100;
  const Int_t nmix = 10;
  Double_t centralityBins[] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,500.0};
  Double_t vertexBins[] = {-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0};
  const Int_t nCentralityBins = 11;
  const Int_t nVertexBins = 10;
  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centralityBins, nVertexBins, vertexBins);
  fPoolMgr->SetTargetValues(trackDepth,0.1,nmix);
}
//___________________________________________________________


Bool_t AliAnalysisTaskfn::TrackPassesOOBPileupCut(AliESDtrack* t, Double_t b){
  if (!t) return true;
  if ((t->GetStatus() & AliESDtrack::kITSrefit) == AliESDtrack::kITSrefit) return true;
  if (t->GetTOFExpTDiff(b, true) + 2500 > 1e-6) return true;
  return false;
}


