
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

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"

#include "Aliresonance.h"
#include "AliPIDResponse.h"
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

ClassImp(Aliresonance)
ClassImp(AliCompactTrack)

Aliresonance::Aliresonance(): AliAnalysisTaskSE(), fOutput(0), fPoolMgr(0X0), fPIDResponse(0), fVevent(0), lAODevent(0), fEventCuts(0), fHistZVertex(0),fHistCentralityEvtCount(0), fHisteventsummary(0),f1Unlike(0),f1Like(0),f1Mix(0)
{

}

Aliresonance::Aliresonance(const char *name): AliAnalysisTaskSE(name), fOutput(0), fPoolMgr(0X0), fPIDResponse(0), fVevent(0), lAODevent(0), fEventCuts(0), fHistZVertex(0),fHistCentralityEvtCount(0), fHisteventsummary(0),f1Unlike(0),f1Like(0),f1Mix(0)
{
  DefineInput(0, TChain::Class()); 
  DefineOutput(1, TList::Class()); 
}


Aliresonance::~Aliresonance()
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

  if(lAODevent) {
    delete lAODevent;
    lAODevent=0x0;
  }

  if (fPoolMgr) { 
    delete fPoolMgr; 
    fPoolMgr = 0x0; 
  }
  
  if (fPIDResponse) { 
    delete fPIDResponse; 
    fPIDResponse = 0x0; 
  }
}

//________________________________________________________________________
void Aliresonance::UserCreateOutputObjects()
{
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  MixingEvents();

  fOutput = new TList();
  fOutput->SetOwner(); 
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
void Aliresonance::UserExec(Option_t *)
{
     
  // Main loop
  // Called for each event
  //AliAODEvent *lAODevent = 0x0;
  
  lAODevent = dynamic_cast <AliAODEvent*> (InputEvent());

  if (!(lAODevent)) {
    AliWarning("ERROR: event not available \n");
    PostData(1, fOutput);
    return;
  }
  
  ////trigger selection/////////////
  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isSelected = 0;
  //isSelected = (maskIsSelected & AliVEvent::kINT7);
  isSelected = (maskIsSelected & (AliVEvent::kHighMultV0 | AliVEvent::kINT7));
  

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
    
  

  fHistZVertex->Fill(zv);
  fHistCentralityEvtCount->Fill(lV0M);

  fHisteventsummary->Fill(1.5);
  


  ////////////////////////////////////////////////////////////////


  TLorentzVector daughterkaon(0.0,0.0,0.0,0.0);
  TLorentzVector daughterksh(0.0,0.0,0.0,0.0);
  TLorentzVector motherkksh(0.0,0.0,0.0,0.0);
  TLorentzVector motherf1(0.0,0.0,0.0,0.0);
  TLorentzVector mothermixf1(0.0,0.0,0.0,0.0);



  std::vector<AlikkshContainer> kkshCandidates;
  std::vector<AlipiContainer> piCandidates;
  AlikkshContainer kksh;
  AlipiContainer pion;
    
  TObjArray* Acceptedpiontracks= new TObjArray;
  Acceptedpiontracks -> SetOwner(kTRUE);   


  Int_t ntracks = lAODevent->GetNumberOfTracks();
  Int_t nv0s = lAODevent->GetNumberOfV0s();
  Double_t pionmass=0.139; 
  Int_t nkksh=0;
  Int_t npion=0;
  int mult=0;
  Double_t tV0mom[3] = {0.0, 0.0, 0.0};
  Double_t trkPt=0.0, trkEta=0.0;


  for(Int_t i = 0; i < ntracks; i++)
    {
      AliVTrack   *track = (AliVTrack*)lAODevent->GetTrack(i);
      if(!track)      continue;
      AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(track);
      if(!aodtrack)      continue;
      if(!aodtrack->TestFilterBit(96))  continue;
      trkPt=aodtrack->Pt();
      trkEta=aodtrack->Eta();
      if (trkPt < 0.15) continue;      
      if (TMath::Abs(trkEta) > 0.8) continue;     
      mult=mult+1;

         
      if(IsKaon(track)){
	daughterkaon.SetXYZM(track->Px(), track->Py(), track->Pz(), 0.4937);

	for(Int_t iv0 = 0; iv0 < nv0s; iv0++)
	  {
	    AliAODv0 *v0=lAODevent->GetV0(iv0);
	    if(!v0) continue;
	    v0->GetPxPyPz(tV0mom);
	    if(!IsV0(v0, lAODevent)) continue;
	    daughterksh.SetXYZM(tV0mom[0], tV0mom[1], tV0mom[2], 0.4976);
	    motherkksh=daughterkaon+daughterksh;
	    kksh.particle.SetXYZM(motherkksh.Px(),motherkksh.Py(),motherkksh.Pz(),motherkksh.M());
	    kksh.charge=aodtrack->Charge();
	    kksh.trackid=i;
	    kkshCandidates.push_back(kksh);
	    nkksh=nkksh+1;
	  }
      }

      if(IsPion(track))
	{
	  Acceptedpiontracks->Add(new AliCompactTrack(track->Px(),track->Py(),track->Pz(),aodtrack->Charge()));
	  pion.particle.SetXYZM(track->Px(), track->Py(), track->Pz(), pionmass);
	  pion.charge=aodtrack->Charge();
	  pion.trackid=i;
	  piCandidates.push_back(pion);
	  npion=npion+1;
	}
 

    }

  if(piCandidates.size()<1 || kkshCandidates.size()<1) 
    {
      PostData(1, fOutput);
      return;
    }

  fHisteventsummary->Fill(2.5);


    /////////////unlike and like sign dist///////////////////////

  for(int j=0;j<kkshCandidates.size();j++)
    {

      kksh=kkshCandidates[j];
      if ((kksh.particle.M()) > 1.04) continue;  

      for(int i=0;i<piCandidates.size();i++)
	{
	  pion=piCandidates[i];
	  if(pion.trackid==kksh.trackid) continue; //to prevent same track counting
	  motherf1.SetPxPyPzE(pion.particle.Px()+kksh.particle.Px(), pion.particle.Py()+kksh.particle.Py(), pion.particle.Pz()+kksh.particle.Pz(), pion.particle.E()+kksh.particle.E());
	  
	  if(motherf1.M()>2.0) continue;
	 
	  if(pion.charge*kksh.charge<0)
	    {
	      
	      f1Unlike->Fill(motherf1.M(),motherf1.Pt());
	      
	    }
	  if(pion.charge*kksh.charge>0)
	    {
	      
	      f1Like->Fill(motherf1.M(),motherf1.Pt());
	    }
	}

    }
 fHisteventsummary->Fill(3.5);

  ///////////////////////////////////////////////////////////////////
      
  //////////////mixed event dist/////////////////// from github code
  double pionmixenergy=0.0;
  AliEventPool* pool = fPoolMgr->GetEventPool(mult,zv); 
  if (!pool) AliInfo("ERROR: Pool unavailable for given centrality ");
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

	      for(int ipi = 0; ipi < numTracks; ipi++)
		{
		  AliCompactTrack *piontrackmix = (AliCompactTrack*) bgTracks->At(ipi);
		  if(!piontrackmix)
		    {
		      AliFatal(Form("ERROR: No mix pool track %d\n",ipi));
		      continue;
		    }
		  
		  AliAODTrack *aodtrackmix  = dynamic_cast<AliAODTrack*>(piontrackmix);
		  
		  for(int j=0;j<kkshCandidates.size();j++)
		    {
		      kksh=kkshCandidates[j];
		      if ((kksh.particle.M()) > 1.04) continue;  
		      if(piontrackmix->Charge()*kksh.charge>0)continue;
		      
		      pionmixenergy=TMath::Sqrt(pionmass*pionmass+piontrackmix->Px()*piontrackmix->Px()+piontrackmix->Py()*piontrackmix->Py()+piontrackmix->Pz()*piontrackmix->Pz());
		      mothermixf1.SetPxPyPzE(piontrackmix->Px()+kksh.particle.Px(), piontrackmix->Py()+kksh.particle.Py(), piontrackmix->Pz()+kksh.particle.Pz(), pionmixenergy+kksh.particle.E());
		      
		      if(mothermixf1.M()>2.0) continue;
		   
		      f1Mix->Fill(mothermixf1.M(),mothermixf1.Pt());
		    }
		}

	    }
	}
      pool->UpdatePool(Acceptedpiontracks);
    }


  /////////////////////////////////////


  fHisteventsummary->Fill(4.5);    
  PostData(1, fOutput);
}
//________________________________________________________________________
void Aliresonance::Terminate(Option_t *)
{
  
}

//________________________________________________________________________

Bool_t Aliresonance::GoodEvent(const AliVVertex *vertex) //all cuts taken from alice code github pp analysis from strangeness and correlation analysis
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


  return kTRUE;

}
//-----------------------------------------------
Bool_t Aliresonance::IsPion(AliVTrack *vtrack)
{

  
  AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(vtrack);
  Double_t nsigmatpcpion=TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodtrack, AliPID::kPion));
  Double_t nsigmatofpion=TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodtrack, AliPID::kPion));
  
  Double_t trkPt=aodtrack->Pt();
  Double_t nsigmacircularcut=sqrt(nsigmatpcpion*nsigmatpcpion + nsigmatofpion*nsigmatofpion);
  Double_t tofsig=0.0;

    if (trkPt<0.5)
    {
      if(nsigmatpcpion>3.0) return kFALSE;      
    }
    else 
    {
      if(nsigmacircularcut>3.0) return kFALSE;
    }
  
  return kTRUE;
}

Bool_t Aliresonance::IsKaon(AliVTrack *vtrack)
{

  AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(vtrack);
  Double_t nsigmatpckaon=TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodtrack, AliPID::kKaon));
  Double_t nsigmatofkaon=TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodtrack, AliPID::kKaon));
  
  Double_t trkPt=aodtrack->Pt();
  Double_t nsigmacircularcut=sqrt(nsigmatpckaon*nsigmatpckaon + nsigmatofkaon*nsigmatofkaon);
  Double_t tofsig=0.0;
  
  if (trkPt<0.45)
    {
      if(nsigmatpckaon>3.0) return kFALSE;      
    }
  else
    {
      if(nsigmacircularcut>3.0) return kFALSE;
    }
    
  return kTRUE;
}



//_________________________________________________________________________________________________
Bool_t Aliresonance::IsV0(AliAODv0 *v0, AliAODEvent *lAODEvent)
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
  if ( TMath::Abs(((pTrack->GetSign()) - (nTrack->GetSign())) ) < 0.1) {
    return kFALSE;
  }
  

  // check track quality cuts
  if (!AcceptAODtracks(pTrack, nTrack)) return kFALSE;


  // topological checks
  if (TMath::Abs(v0->DcaPosToPrimVertex()) < 0.06) {
    return kFALSE;
  }
  if (TMath::Abs(v0->DcaNegToPrimVertex()) < 0.06) {
    return kFALSE;
  }
  if (TMath::Abs(v0->DcaV0Daughters()) > 1.0 ) {
    return kFALSE;
  }
  if (TMath::Abs(v0->DcaV0ToPrimVertex()) > 0.3) {
    return kFALSE;
  }
 
  Bool_t fCheckOOBPileup=kTRUE;
  if (fCheckOOBPileup) {
    Double_t bfield = lAODEvent->GetMagneticField();
    if(!TrackPassesOOBPileupCut(pTrack, bfield) &&
       !TrackPassesOOBPileupCut(nTrack, bfield)) return kFALSE;
  }
   
  if ((TMath::Abs(v0->CosPointingAngle(primVtx)) < 0.97) || (TMath::Abs(v0->CosPointingAngle(primVtx)) >= 1 ) ) {
    return kFALSE;
  }
  
  if (TMath::Abs(v0->RapK0Short()) > 0.5) {
    return kFALSE;
  }
  
  if (TMath::Abs(v0->Eta())> 0.8) {
    return kFALSE;
  }
 
  Double_t radius = v0->RadiusV0();
  if (( radius < 0.5 ) || ( radius > 200 ) ) {
    return kFALSE;
  }

  Double_t lV0TotalMomentum  = TMath::Sqrt(v0->Ptot2V0());
  Double_t fLength = v0->DecayLength(primVtx);
  if( TMath::Abs(0.497*fLength/lV0TotalMomentum) > 20)
    {
      AliDebugClass(2, "Failed Lifetime Cut on positive track V0");
      return kFALSE;
    }

  
  // Competing v0rejection for K0s vs Lambda
  Double_t altmass=0.0;
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
  if(! ((negnsTPC <= 4) && (posnsTPC <= 4)) ) {  
  return kFALSE;
  }


  if(v0->MassK0Short()<0.482 || v0->MassK0Short()>0.512)
   {
     return kFALSE;
   }

  return kTRUE;
  

}

//_________________________________________________________________________________________________
void Aliresonance::MixingEvents() //from git PWGCF correlations code
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



Bool_t Aliresonance::TrackPassesOOBPileupCut(AliAODTrack* t, Double_t b){
  if (!t) return true;
  if ((t->GetStatus() & AliVTrack::kITSrefit) == AliVTrack::kITSrefit) return true;
  if (t->GetTOFExpTDiff(b, true) + 2500 > 1e-6) return true;
  return false;
}


Bool_t Aliresonance::AcceptAODtracks(AliAODTrack *pTrack, AliAODTrack *nTrack)
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

