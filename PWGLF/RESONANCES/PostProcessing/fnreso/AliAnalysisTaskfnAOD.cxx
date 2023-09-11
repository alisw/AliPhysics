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

#include "AliAnalysisTaskfnAOD.h"
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

ClassImp(AliAnalysisTaskfnAOD)
ClassImp(AliReducedTrack)

AliAnalysisTaskfnAOD::AliAnalysisTaskfnAOD(): 
AliAnalysisTaskSE(), 
fOutput(0), 
fPoolMgr(0X0), 
fPIDResponse(0), 
fVevent(0), 
lESDevent(0), 
lAODevent(0), 
fEventCuts(0), 
fESDtrackCuts(0), 
fHistVz(0),
fHistCentrality(0), 
fHisteventsummary(0),
f1Unlike(0),
f1Like(0),
f1Mix(0),
hist1(0),
hist2(0),
hist3(0),
hist4(0),
hist5(0),
hist6(0),
hist7(0),
fFilterBit(32),
kkshmasscut(1.04),
nsigtpcpion(2),
nsigtofpion(3),
nsigtpckaon(2),
nsigtofkaon(3),
dcaxypos(0.06),
dcaxyneg(0.06),
dcav0daugh(1.0),
dcav0pv(0.3),
cospa(0.97),
lowrad(0.5),
lifetime(15),
pidpion(4),
nCRcut(70.),
ratiocrfccut(0.8),
chi2globalcut(36.0),
chi2cut(36.0)  
{

}

AliAnalysisTaskfnAOD::AliAnalysisTaskfnAOD(const char *name): 
AliAnalysisTaskSE(name), 
fOutput(0), 
fPoolMgr(0X0), 
fPIDResponse(0),
fVevent(0), 
lESDevent(0), 
lAODevent(0), 
fEventCuts(0), 
fESDtrackCuts(0), 
fHistVz(0),
fHistCentrality(0), 
fHisteventsummary(0),
f1Unlike(0),
f1Like(0),
f1Mix(0),
hist1(0),
hist2(0),
hist3(0),
hist4(0),
hist5(0),
hist6(0),
hist7(0),
fFilterBit(32),
kkshmasscut(1.04),
nsigtpcpion(2),
nsigtofpion(3),
nsigtpckaon(2),
nsigtofkaon(3),
dcaxypos(0.06),
dcaxyneg(0.06),
dcav0daugh(1.0),
dcav0pv(0.3),
cospa(0.97),
lowrad(0.5),
lifetime(15),
pidpion(4),
nCRcut(70.),
ratiocrfccut(0.8),
chi2globalcut(36.0),
chi2cut(36.0)  
{
  DefineInput(0, TChain::Class()); 
  DefineOutput(1, TList::Class()); 
}


AliAnalysisTaskfnAOD::~AliAnalysisTaskfnAOD()
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

  if(lAODevent) {
    delete lAODevent;
    lAODevent=0x0;
  }

  if(fESDtrackCuts) {
    delete fESDtrackCuts;
    fESDtrackCuts=0x0;
  }


  if (fPoolMgr) { delete fPoolMgr; fPoolMgr = 0x0; }
  if (fPIDResponse) { delete fPIDResponse; fPIDResponse = 0x0; }
}

//________________________________________________________________________
void AliAnalysisTaskfnAOD::UserCreateOutputObjects()
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


  EventMixing();

  fOutput = new TList();
  fOutput->SetOwner(); 
  OpenFile(1);
  fHisteventsummary            = new TH1F("fHisteventsummary", "Event cut summary", 5,0.0,5.0);       
  fHistVz            = new TH1F("fHistZVertex", "Z vertex distribution", 100,-10,10);       
  fHistCentrality = new TH1F("fHistCentrality", "Centrality distribution", 100,0,100);

  Int_t bins[2]={250, 200};
  Double_t xmin[2]={1.0, 0.0};
  Double_t xmax[2]={2.0, 20.0};
  f1Unlike = new THnSparseD("f1Unlike", "Unlike histogram", 2, bins, xmin, xmax);
  f1Like = new THnSparseD("f1Like", "Like histogram", 2, bins, xmin, xmax);
  f1Mix = new THnSparseD("f1Mix", "Mix histogram", 2, bins, xmin, xmax);
  hist1=new TH1D("hist1", "hist1", 85, 0, 170);
  hist2=new TH1D("hist2", "hist2", 20, 0, 2);
  hist3=new TH1D("hist3", "hist3", 10, 0, 10);
  hist4=new TH1D("hist4", "hist4", 25, 0, 50);
  hist5=new TH1D("hist5", "hist5", 50, -5, 5);
  hist6=new TH1D("hist6", "hist6", 50, -5, 5);
  hist7=new TH1D("hist7", "hist7", 25, 0, 50);


  f1Unlike->Sumw2();
  f1Like->Sumw2();
  f1Mix->Sumw2();


  fOutput->Add(fHisteventsummary);
  fOutput->Add(fHistVz);
  fOutput->Add(fHistCentrality);
  fOutput->Add(f1Unlike);
  fOutput->Add(f1Like);
  fOutput->Add(f1Mix);
  fOutput->Add(hist1);
  fOutput->Add(hist2);
  fOutput->Add(hist3);
  fOutput->Add(hist4);
  fOutput->Add(hist5);
  fOutput->Add(hist6);
  fOutput->Add(hist7);

  PostData(1, fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskfnAOD::UserExec(Option_t *)
{
     
  // Main loop
  // Called for each event
  //AliAODEvent *lAODevent = 0x0;
  


  lAODevent = dynamic_cast <AliAODEvent*> (InputEvent());

  if (!(lAODevent)) {
    AliWarning("ERROR: AOD event not available \n");
    PostData(1, fOutput);
    return;
  }
  
  ////trigger/////////////
 
  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isSelected = 0;
  isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
 

  if ( ! isSelected)
    {
      
      PostData(1, fOutput);
      return;
    }

  fHisteventsummary->Fill(0.5);

  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fVevent){
    printf("ERROR: fVevent not available\n");
    return;
  }


  //event selection  
  const AliVVertex *vertex = fVevent->GetPrimaryVertex();
  if (!GoodEvent(vertex)) return; 
  double zv=vertex->GetZ();
  ////////////////////////////////


  //////centrality selection/////////
  Float_t lV0M;
  Int_t lEvSelCode = 300;
  AliMultSelection *MultSelection = (AliMultSelection*) lAODevent -> FindListObject("MultSelection");
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
    
  

  fHistVz->Fill(zv);
  fHistCentrality->Fill(lV0M);
  fHisteventsummary->Fill(1.5);
  


  std::vector<AlikkshPair> kks0Candidates;
  std::vector<Alipi> pionCandidates;
  AlikkshPair kks0;
  Alipi pion;
    
  TLorentzVector kaon(0.0,0.0,0.0,0.0);
  TLorentzVector kshort(0.0,0.0,0.0,0.0);
  TLorentzVector kaonkshort(0.0,0.0,0.0,0.0);
  TLorentzVector motherf1(0.0,0.0,0.0,0.0);
  TLorentzVector motherf1mix(0.0,0.0,0.0,0.0);
    


  Int_t ntracks = lAODevent->GetNumberOfTracks();
  Int_t nv0s = lAODevent->GetNumberOfV0s();
  Double_t pionmass=0.139;     
  Double_t kaonmass=0.4937;     


  TObjArray* selectedpiontracks= new TObjArray;
  selectedpiontracks -> SetOwner(kTRUE);


  Int_t nkks0pair=0;
  Int_t npion=0;
  int mult=0;
  Double_t trkPt=0.0, trkEta=0.0;
  Float_t nCrossedRowsTPC=0.0, ratiocrfindcls=0.0;
  Int_t findablecls=0;
  Double_t chi2perclsTPC=0.0, chi2perclsITS=0.0, chi2constglobal=0.0;
  Float_t dca[2]{0., 0.};



  for(Int_t itr = 0; itr < ntracks; itr++)
    {
      AliVTrack   *track = (AliVTrack*)lAODevent->GetTrack(itr);
      if(!track)      continue;
      //AliAODTrack* aodtrack = dynamic_cast <AliAODTrack*> (fVevent->GetTrack(itr));
      AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(track);
      if(!aodtrack)      continue;
      if(!aodtrack->TestFilterBit(fFilterBit))  continue;
      trkPt=aodtrack->Pt();
      trkEta=aodtrack->Eta();
      if (trkPt < 0.15) continue;      
      if (TMath::Abs(trkEta) > 0.8) continue; 
      mult=mult+1;

      //track cuts fill
      nCrossedRowsTPC = aodtrack->GetTPCClusterInfo(2,1);
      findablecls = aodtrack->GetTPCNclsF();
      ratiocrfindcls = nCrossedRowsTPC/findablecls;
      chi2perclsTPC = aodtrack->GetTPCchi2perCluster();
      AliAODVertex *vertex = aodtrack->GetProdVertex();
      if (vertex)
        {
          if (vertex->GetType() == AliAODVertex::kKink) continue;
	}

      if (!(aodtrack->GetStatus() & AliVTrack::kTPCrefit)) continue;
      if (!(aodtrack->GetStatus() & AliVTrack::kITSrefit)) continue;
      chi2constglobal = aodtrack->GetChi2TPCConstrainedVsGlobal();
      aodtrack->GetImpactParameters(dca[0], dca[1]);
      chi2perclsITS= aodtrack->GetITSchi2()/aodtrack->GetITSNcls();

      if (nCrossedRowsTPC > nCRcut  && ratiocrfindcls > ratiocrfccut && chi2constglobal < chi2globalcut && chi2perclsITS < chi2cut)
	{
	
      hist1->Fill(nCrossedRowsTPC);
      hist2->Fill(ratiocrfindcls);
      hist3->Fill(chi2perclsTPC);
      hist4->Fill(chi2constglobal);
      hist5->Fill(dca[0]);
      hist6->Fill(dca[1]);
      hist7->Fill(chi2perclsITS);
      /////////////////////////

      if(IsKaon(track)){
	kaon.SetXYZM(track->Px(), track->Py(), track->Pz(), kaonmass);

	for(Int_t itrv0 = 0; itrv0 < nv0s; itrv0++)
	  {
	    AliAODv0 *v0=lAODevent->GetV0(itrv0);
	    if(!v0) continue;
	   
	    if(!IsV0(v0, lAODevent)) continue;
	    kshort.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), 0.4976);
	    kaonkshort=kaon+kshort;
	    kks0.charge=aodtrack->Charge();
	    kks0.trkid=itr;
	    kks0.particle.SetXYZM(kaonkshort.Px(),kaonkshort.Py(),kaonkshort.Pz(),kaonkshort.M());
	    kks0Candidates.push_back(kks0);
	    nkks0pair=nkks0pair+1;
	    

	  }
      }


      if(IsPion(track))
	{
	  selectedpiontracks->Add(new AliReducedTrack(track->Px(),track->Py(),track->Pz(),aodtrack->Charge()));
	  pion.charge=aodtrack->Charge();
	  pion.trkid=itr;
	  pion.particle.SetXYZM(track->Px(), track->Py(), track->Pz(), pionmass);
	  pionCandidates.push_back(pion);
	  npion=npion+1;
	}
      
	}

    }
  
  Int_t piontracksize = pionCandidates.size();
  Int_t kkshorttracksize = kks0Candidates.size();


  if(piontracksize<1 || kkshorttracksize<1) 
    {
      PostData(1, fOutput);
      return;
    }

  fHisteventsummary->Fill(2.5);


  for (const auto& kks0 : kks0Candidates)
    {
      if (kks0.particle.M() > kkshmasscut)
	continue;
      // Creating same event pair
      for (const auto& pion : pionCandidates)
	{
	  if (pion.trkid == kks0.trkid)
            continue;

	  motherf1 = pion.particle + kks0.particle;
	  if (motherf1.M() > 2.0)
            continue;
	  if (TMath::Abs(motherf1.Rapidity())>=0.5)                                                                                         
	    continue;       

	  if (pion.charge * kks0.charge < 0)
	    {
	      // Fill unlike histo
	      f1Unlike->Fill(motherf1.M(), motherf1.Pt());
	      
	    }
	  else if (pion.charge * kks0.charge > 0)
	    {
	      // Fill like histo
	      f1Like->Fill(motherf1.M(), motherf1.Pt());
	      
	    }
	}
    }



 fHisteventsummary->Fill(3.5);

  ///////////////////////////////////////////////////////////////////
      
 // Mixed event
 double pionmixenergy = 0.0;
 AliEventPool* pool = fPoolMgr->GetEventPool(lV0M, zv);
 if (pool && pool->IsReady())
   {
     Int_t nMix = 5; // Set the number of mixed events

     for (Int_t jMix = 0; jMix < nMix; jMix++)
       {
	 TObjArray* bgTracks = pool->GetRandomEvent();
	 Int_t numTracks = bgTracks->GetEntriesFast();

	 for (const auto& bgTrack : *bgTracks)
	   {
	     AliReducedTrack* piontrackmix = dynamic_cast<AliReducedTrack*>(bgTrack);
	     if (!piontrackmix)
	       {
		 AliFatal(Form("ERROR: Could not receive mix pool track %d\n", bgTrack->GetUniqueID()));
		 continue;
	       }

	     AliAODTrack* aodtrackmix = dynamic_cast<AliAODTrack*>(piontrackmix);

	     for (const auto& kks0 : kks0Candidates)
	       {
		 if (kks0.particle.M() > kkshmasscut)
		   continue;

		 if (piontrackmix->Charge() * kks0.charge > 0)
		   continue;

		 pionmixenergy = TMath::Sqrt(pionmass * pionmass + piontrackmix->Px() * piontrackmix->Px() + piontrackmix->Py() * piontrackmix->Py() + piontrackmix->Pz() * piontrackmix->Pz());
		 motherf1mix.SetPxPyPzE(piontrackmix->Px()+kks0.particle.Px(), piontrackmix->Py()+kks0.particle.Py(), piontrackmix->Pz()+kks0.particle.Pz(), pionmixenergy+kks0.particle.E());
		 if (motherf1mix.M() > 2.0)
		   continue;
		 if (TMath::Abs(motherf1mix.Rapidity())>=0.5)                                         
		   continue;
		 // Fill mix histo
		 f1Mix->Fill(motherf1mix.M(), motherf1mix.Pt());
		 
	       }
	   }
       }
   }

 if (pool)
   pool->UpdatePool(selectedpiontracks);


  /////////////////////////////////////

  

    
  PostData(1, fOutput);
}
//________________________________________________________________________
void AliAnalysisTaskfnAOD::Terminate(Option_t *)
{
  
}

//________________________________________________________________________

Bool_t AliAnalysisTaskfnAOD::GoodEvent(const AliVVertex *vertex) //all cuts taken from alice code github pp analysis from strangeness and correlation analysis
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

  
  /*
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
  */
  if( lAODevent->IsIncompleteDAQ() ) 
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

    fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);


    return kTRUE;

}
//-----------------------------------------------
Bool_t AliAnalysisTaskfnAOD::IsPion(AliVTrack *vtrack)
{

  
  AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(vtrack);
  Double_t nsigmatpcpion=TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodtrack, AliPID::kPion));
  Double_t nsigmatofpion=TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodtrack, AliPID::kPion));
  Bool_t TOFHIT=kFALSE;
  
  TOFHIT = HasTOF(aodtrack);

  
  if(!TOFHIT)
    {
      if(nsigmatpcpion>nsigtpcpion)return kFALSE;
    }

  else
    {
      if(nsigmatofpion>nsigtofpion)return kFALSE;
    }

  return kTRUE;
}

Bool_t AliAnalysisTaskfnAOD::IsKaon(AliVTrack *vtrack)
{

  AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(vtrack);
  Double_t nsigmatpckaon=TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodtrack, AliPID::kKaon));
  Double_t nsigmatofkaon=TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodtrack, AliPID::kKaon));
  Bool_t TOFHIT=kFALSE;
  
  
  TOFHIT = HasTOF(aodtrack);

  if(!TOFHIT)
    {
      if(nsigmatpckaon>nsigtpckaon)return kFALSE;
    }
  else
    {
      if(nsigmatofkaon>nsigtofkaon)return kFALSE;
    }

    
  return kTRUE;
}

//--------------------------------------------


Bool_t AliAnalysisTaskfnAOD::HasTOF(AliAODTrack *track)
{
  bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  //const float len = track->GetIntegratedLength();
  //bool hasTOF = hasTOFout && hasTOFtime && (len > 350.);
  bool hasTOF = hasTOFout && hasTOFtime;
  return hasTOF;
}



//_________________________________________________________________________________________________
Bool_t AliAnalysisTaskfnAOD::IsV0(AliAODv0 *v0, AliAODEvent *lAODEvent)
{

  
  if (v0->GetOnFlyStatus()) {
    return kFALSE; 
  }
  
  Double_t xPrimaryVertex = lAODEvent->GetPrimaryVertex()->GetX();
  Double_t yPrimaryVertex = lAODEvent->GetPrimaryVertex()->GetY();
  Double_t zPrimaryVertex = lAODEvent->GetPrimaryVertex()->GetZ();
  Double_t primVtx[3] = {xPrimaryVertex, yPrimaryVertex, zPrimaryVertex};


  
  // retrieve the V0 daughters
  AliAODTrack *pTrack = (AliAODTrack *)(v0->GetSecondaryVtx()->GetDaughter(0));
  AliAODTrack *nTrack = (AliAODTrack *)(v0->GetSecondaryVtx()->GetDaughter(1));;



  // filter like-sign V0
  if ( TMath::Abs(((pTrack->Charge()) - (nTrack->Charge())) ) < 0.1) {
    return kFALSE;
  }  

  if (!AcceptAODtracks(pTrack, nTrack)) return kFALSE;


  // topological checks
  if (TMath::Abs(v0->DcaPosToPrimVertex()) < dcaxypos) {
    return kFALSE;
  }
  if (TMath::Abs(v0->DcaNegToPrimVertex()) < dcaxyneg) {
    return kFALSE;
  }
  if (TMath::Abs(v0->DcaV0Daughters()) > dcav0daugh) {
    return kFALSE;
  }
  if (TMath::Abs(v0->DcaV0ToPrimVertex()) > dcav0pv) {
    return kFALSE;
  }


  Bool_t fCheckOOBPileup=kTRUE;
  if (fCheckOOBPileup) {
    Double_t bfield = lAODEvent->GetMagneticField();
    if(!TrackPassesOOBPileupCut(pTrack, bfield) &&
       !TrackPassesOOBPileupCut(nTrack, bfield)) return kFALSE;
  }

   
  if ((TMath::Abs(v0->CosPointingAngle(primVtx)) < cospa) || (TMath::Abs(v0->CosPointingAngle(primVtx)) >= 1 ) ) {
    return kFALSE;
  }
  
  if (TMath::Abs(v0->Eta())> 0.8) {
    return kFALSE;
  }

  if (TMath::Abs(v0->RapK0Short()) > 0.8) {
    AliDebugClass(2, "Failed check on V0 rapidity");
    return kFALSE;
  }

 
  Double_t radius = v0->RadiusV0();
  if (( radius < lowrad ) || ( radius > 200 ) ) {
    return kFALSE;
  }

  Double_t lV0TotalMomentum  = TMath::Sqrt(v0->Ptot2V0());
  //Double_t fLength = v0->DecayLength(primVtx);
  Double_t fLength = TMath::Sqrt(TMath::Power(v0->DecayVertexV0X() - xPrimaryVertex,2) +
				 TMath::Power(v0->DecayVertexV0Y() - yPrimaryVertex,2) +
				 TMath::Power(v0->DecayVertexV0Z() - zPrimaryVertex,2));

  if( TMath::Abs(0.497*fLength/lV0TotalMomentum) > lifetime)
    {
      AliDebugClass(2, "Failed Lifetime Cut on positive track V0");
      return kFALSE;
    }


  Double_t altmass=0.0;
  Double_t fMass = 0.497614;
  Double_t fTolerance = 0.03;
  altmass = v0->MassK0Short();


  if ((TMath::Abs(altmass - fMass)) > fTolerance) {
    AliDebugClass(2,"V0 is not in the expected inv mass range  Mass");
    return kFALSE;
  }


  // Competing v0rejection for K0s vs Lambda
  Double_t fToleranceVeto=0.004;
  
  altmass = v0->MassLambda();
  if ((TMath::Abs(altmass - 1.115683)) < fToleranceVeto) {
    AliDebugClass(2, "Failed competing V0 rejection check");
    return kFALSE;
  }
  altmass = v0->MassAntiLambda();
  if ((TMath::Abs(altmass - 1.115683)) < fToleranceVeto) {
    AliDebugClass(2, "Failed competing anti-V0 rejection check");
    return kFALSE;
  }


  altmass = v0->MassK0Short();

  // check PID
  Double_t posnsTPC   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion));
  Double_t negnsTPC   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion));
  if(! ((negnsTPC <= pidpion) && (posnsTPC <= pidpion)) ) {  
    return kFALSE;
  }

  return kTRUE;
  
  

  
}

//_________________________________________________________________________________________________
void AliAnalysisTaskfnAOD::EventMixing()
{

  ////////////////////////////
  // Set-up for Mixed Event //
  ////////////////////////////
  const Int_t trackDepth = 10000;
  const Int_t poolsize = 100;
  const Int_t nmix = 5;
  Double_t centralityBins[] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0};
  Double_t vertexBins[] = {-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0};
  const Int_t nCentralityBins = 11;
  const Int_t nVertexBins = 10;
  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centralityBins, nVertexBins, vertexBins);
  fPoolMgr->SetTargetValues(trackDepth,0.1,nmix);
}
//___________________________________________________________


Bool_t AliAnalysisTaskfnAOD::TrackPassesOOBPileupCut(AliAODTrack* t, Double_t b){
  if (!t) return true;
  if ((t->GetStatus() & AliVTrack::kITSrefit) == AliVTrack::kITSrefit) return true;
  if (t->GetTOFExpTDiff(b, true) + 2500 > 1e-6) return true;
  return false;
}


//___________________________________________________________________

Bool_t AliAnalysisTaskfnAOD::AcceptAODtracks(AliAODTrack *pTrack, AliAODTrack *nTrack)
{

  Double_t peta = pTrack->Eta();
  Double_t neta = nTrack->Eta();
  Double_t ppT = pTrack->Pt();
  Double_t npT = nTrack->Pt();

  if (TMath::Abs(peta) > 0.8 || TMath::Abs(neta) > 0.8) 
    {
      return kFALSE;
    }
  
  if (ppT < 0.15 || npT < 0.15) 
    {
      return kFALSE;
    }
  
  if( pTrack->GetKinkIndex(0)>0 || nTrack->GetKinkIndex(0)>0 )
    {
      return kFALSE;
    }

  if( !(pTrack->GetStatus() & AliAODTrack::kTPCrefit)) return kFALSE;
  if( !(nTrack->GetStatus() & AliAODTrack::kTPCrefit)) return kFALSE;

  if (pTrack->GetTPCchi2perCluster()>100) return kFALSE;
  if (nTrack->GetTPCchi2perCluster()>100) return kFALSE;


  Double_t nCrossedRowsTPCpos = pTrack->GetTPCClusterInfo(2,1);
  Double_t nCrossedRowsTPCneg = nTrack->GetTPCClusterInfo(2,1);
  Double_t findablepos = pTrack->GetTPCNclsF();
  Double_t findableneg = nTrack->GetTPCNclsF();
  if (nCrossedRowsTPCpos < 70 || nCrossedRowsTPCneg < 70) return kFALSE;
  if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) return kFALSE;  
  if (nCrossedRowsTPCpos/findablepos < 0.8) return kFALSE;
  if (nCrossedRowsTPCneg/findableneg < 0.8) return kFALSE;

  return kTRUE;
}


