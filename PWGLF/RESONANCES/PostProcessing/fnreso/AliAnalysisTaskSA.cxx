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

#include "AliAnalysisTaskSA.h"
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

ClassImp(AliAnalysisTaskSA)
ClassImp(AliCompTrack)

AliAnalysisTaskSA::AliAnalysisTaskSA(): 
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
fHisteventmult(0),
kstarUnlike(0),
kstarLike(0),
//kstarposLike(0),
//kstarnegLike(0),
kstarMix(0),
//fHistpionpt(0),
//fHistkaonpt(0),
/*fHistnsigtpcpion(0),
fHistnsigtpckaon(0),
fHistnsigtofpion(0),
fHistnsigtofkaon(0),*/
frame(2),
fListTRKCorr(NULL),
fListNUACorr(NULL),
fListV0MCorr(NULL),
fHCorrectMCposChrg(NULL),
fHCorrectMCnegChrg(NULL),
fHCorrectNUAposChrg(NULL),  
fHCorrectNUAnegChrg(NULL), 
fHCorrectEVNTWGTChrg(NULL)
{

}

AliAnalysisTaskSA::AliAnalysisTaskSA(const char *name): 
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
fHisteventmult(0),
kstarUnlike(0),
kstarLike(0),
//kstarposLike(0),
//kstarnegLike(0),
kstarMix(0),
//fHistpionpt(0),
//fHistkaonpt(0),
/*fHistnsigtpcpion(0),
fHistnsigtpckaon(0),
fHistnsigtofpion(0),
fHistnsigtofkaon(0),*/
frame(2),
fListTRKCorr(NULL),
fListNUACorr(NULL),
fListV0MCorr(NULL),
fHCorrectMCposChrg(NULL),
fHCorrectMCnegChrg(NULL),
fHCorrectNUAposChrg(NULL),  
fHCorrectNUAnegChrg(NULL), 
fHCorrectEVNTWGTChrg(NULL)
{
  DefineInput(0, TChain::Class()); 
  DefineOutput(1, TList::Class()); 
}


AliAnalysisTaskSA::~AliAnalysisTaskSA()
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


  if(fListTRKCorr)  delete fListTRKCorr;
  if(fListNUACorr)  delete fListNUACorr;
  if(fListV0MCorr)  delete fListV0MCorr;



  if (fPoolMgr) { delete fPoolMgr; fPoolMgr = 0x0; }
  if (fPIDResponse) { delete fPIDResponse; fPIDResponse = 0x0; }
}

//________________________________________________________________________
void AliAnalysisTaskSA::UserCreateOutputObjects()
{
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  if(!fESDtrackCuts ){
    fESDtrackCuts = new AliESDtrackCuts();
    fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1); 
    fESDtrackCuts->SetPtRange(0.15,20.0);
    fESDtrackCuts->SetEtaRange(-0.8,0.8);
    //fESDtrackCuts->SetRapRange(-0.5,0.5);
  }


  EventMixing();

  fOutput = new TList();
  fOutput->SetOwner(); 
  OpenFile(1);
  fHisteventsummary            = new TH1F("fHisteventsummary", "Event cut summary", 5,0.0,5.0);       
  fHisteventmult            = new TH1F("fHisteventmult", "Event multiplicity", 10,0.0,10.0);       
  fHistVz            = new TH1F("fHistZVertex", "Z vertex distribution", 100,-10,10);       
  fHistCentrality = new TH1F("fHistCentrality", "Centrality distribution", 100,0,100);

  /*
  Int_t bins[4]={90, 18, 200, 20};
  Double_t xmin[4]={0.6, 0.0, 0.0, -1.0};
  Double_t xmax[4]={1.5, 90.0, 20.0, 1.0};
  */
  Int_t bins[4]={90, 100, 200, 20};
  Double_t xmin[4]={0.6, 0.0, 0.0, -1.0};
  Double_t xmax[4]={1.5, 100.0, 20.0, 1.0};
  kstarUnlike = new THnSparseD("kstarUnlike", "Unlike histogram", 4, bins, xmin, xmax);
  kstarLike = new THnSparseD("kstarLike", "Like histogram", 4, bins, xmin, xmax);
  //kstarposLike = new THnSparseD("kstarposLike", "Pos Like histogram", 3, bins, xmin, xmax);
  //kstarnegLike = new THnSparseD("kstarnegLike", "Neg Like histogram", 3, bins, xmin, xmax);
  kstarMix = new THnSparseD("kstarMix", "Mix histogram", 4, bins, xmin, xmax);

  //fHistpionpt = new TH1D("pionpt", "pionpt", 200, 0.0, 20.0);
  //fHistkaonpt = new TH1D("kaonpt", "kaonpt", 200, 0.0, 20.0);
  /*fHistnsigtpcpion = new TH1D ("nsigmatpcpion","nsigmatpcpion",100,-5.0,5.0);
  fHistnsigtofpion = new TH1D ("nsigmatofpion","nsigmatofpion",100,-5.0,5.0);
  fHistnsigtpckaon = new TH1D ("nsigmatpckaon","nsigmatpckaon",100,-5.0,5.0);
  fHistnsigtofkaon = new TH1D ("nsigmatofkaon","nsigmatofkaon",100,-5.0,5.0);
  */

  kstarUnlike->Sumw2();
  kstarLike->Sumw2();
  //kstarposLike->Sumw2();
  //kstarnegLike->Sumw2();
  kstarMix->Sumw2();
  //fHistpionpt->Sumw2();
  //fHistkaonpt->Sumw2();
  /*fHistnsigtpcpion->Sumw2();
  fHistnsigtpckaon->Sumw2();
  fHistnsigtofpion->Sumw2();
  fHistnsigtofkaon->Sumw2();
  */
  if(fListTRKCorr){
    std::cout<<"\n UserCreateOutputObject::Info() Tlist for MC tracking Efficiency Found.!!\n"<<std::endl;
  }
  else{
    std::cout<<"\n\n ******* WARNING No TList for Trk Efficiency Correction!!\n using TrkWgt = 1.0 \n "<<std::endl;
  }

  if(fListNUACorr){
    std::cout<<"\n UserCreateOutputObject::Info() Tlist for NUA Correction Found.!!\n"<<std::endl;
    //fListNUACorr->ls(); 
  }
  else{
    std::cout<<"\n\n ******* WARNING No TList NUA Correction!!\n using NUAWgt = 1.0 \n "<<std::endl;
  }

  if(fListV0MCorr){
    std::cout<<"\n UserCreateOutputObject::Info() Tlist for EVNTWGT Correction Found.!!\n"<<std::endl;
    //fListV0MCorr->ls(); 
  }
  else{
    std::cout<<"\n\n ******* WARNING No TList EVNTWGT Correction!!\n using EVNTWGTWgt = 1.0 \n "<<std::endl;
  }
 





  fOutput->Add(fHisteventsummary);
  fOutput->Add(fHisteventmult);
  fOutput->Add(fHistVz);
  fOutput->Add(fHistCentrality);
  fOutput->Add(kstarUnlike);
  fOutput->Add(kstarLike);
  //fOutput->Add(kstarposLike);
  //fOutput->Add(kstarnegLike);
  fOutput->Add(kstarMix);
  //fOutput->Add(fHistpionpt);
  //fOutput->Add(fHistkaonpt);
  /*fOutput->Add(fHistnsigtpcpion);
  fOutput->Add(fHistnsigtpckaon);
  fOutput->Add(fHistnsigtofpion);
  fOutput->Add(fHistnsigtofkaon);
  */
  PostData(1, fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskSA::UserExec(Option_t *)
{
     
  // Main loop
  
  
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
 
  //getting correction histograms for eff and NUA
  Int_t runNumber = lAODevent->GetRunNumber();
  if(fListTRKCorr) GetMCCorrectionHist(runNumber,lV0M);  //use centrality dependent MC efficiency (Temporary!!)

  

  if(fListNUACorr){
    GetNUACorrectionHist(runNumber,0);         //Charge
  }
  
  if(fListV0MCorr){
    GetV0MCorrectionHist(runNumber,0);         //Charge
  }


  
  std::vector<Alikaon> kaonCandidates;
  std::vector<Alipion> pionCandidates;
  Alikaon kaon;
  Alipion pion;
    
  TLorentzVector daughterkaon(0.0,0.0,0.0,0.0);
  TLorentzVector daughterpion(0.0,0.0,0.0,0.0);
  TLorentzVector motherkstar(0.0,0.0,0.0,0.0);
  TLorentzVector motherkstarmix(0.0,0.0,0.0,0.0);
    


  Int_t ntracks = lAODevent->GetNumberOfTracks();
  Int_t nv0s = lAODevent->GetNumberOfV0s();
  Double_t pionmass=0.139;     
  Double_t kaonmass=0.4937;     


  TObjArray* selectedpiontracks= new TObjArray;
  selectedpiontracks -> SetOwner(kTRUE);


  Int_t nkaon=0;
  Int_t npion=0;
  int mult=0;
  Double_t trackpt=0.0;
  Double_t trkPt=0.0, trkEta=0.0;


  for(Int_t itr = 0; itr < ntracks; itr++)
    {
      //AliVTrack   *track = (AliVTrack*)lESDevent->GetTrack(itr);
      AliVTrack   *track = (AliVTrack*)lAODevent->GetTrack(itr);
      if(!track)      continue;
      //AliESDtrack *esdtrack  = dynamic_cast<AliESDtrack*>(track);
      AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(track);
      if(!aodtrack)      continue;
      //if(!fESDtrackCuts->AcceptTrack(esdtrack))  continue;
      if(!aodtrack->TestFilterBit(32))  continue;
      trkPt=aodtrack->Pt();
      trkEta=aodtrack->Eta();
      if (trkPt < 0.15) continue;
      if (TMath::Abs(trkEta) > 0.8) continue;
      mult=mult+1;
      

      if(IsKaon(track))
	{
	  daughterkaon.SetXYZM(track->Px(), track->Py(), track->Pz(), kaonmass);
	  kaon.charge=aodtrack->Charge();
	  kaon.trkid=itr;
	  kaon.particle.SetXYZM(daughterkaon.Px(),daughterkaon.Py(),daughterkaon.Pz(),daughterkaon.M());
	  kaonCandidates.push_back(kaon);
	  nkaon=nkaon+1;
	  fHisteventmult->Fill(0.5);
	  //fHistkaonpt->Fill(trkPt);	  
	}
      


      if(IsPion(track))
	{

	  selectedpiontracks->Add(new AliCompTrack(track->Px(),track->Py(),track->Pz(),aodtrack->Charge()));
	  pion.charge=aodtrack->Charge();
	  pion.trkid=itr;
	  pion.particle.SetXYZM(track->Px(), track->Py(), track->Pz(), pionmass);
	  pionCandidates.push_back(pion);
	  npion=npion+1;
	  fHisteventmult->Fill(1.5);	  
	  //fHistpionpt->Fill(trkPt);

	}
 
    }

  Int_t piontracksize = pionCandidates.size();
  Int_t kaontracksize = kaonCandidates.size();


  if(piontracksize<1 || kaontracksize<1) 
    {
      PostData(1, fOutput);
      return;
    }

  fHisteventsummary->Fill(2.5);

  Double_t Costhetastar=0.0;

  for (const auto& kaon : kaonCandidates)
    {
      fHisteventmult->Fill(2.5);	  
      // Creating same event pair
      for (const auto& pion : pionCandidates)
	{
	  fHisteventmult->Fill(3.5);	  
	  if (pion.trkid == kaon.trkid)
            continue;
	  fHisteventmult->Fill(4.5);	  
	  motherkstar = pion.particle + kaon.particle;
	  if (frame==0)
	    Costhetastar = CosThetaStar(motherkstar, pion.particle, kaon.particle);
	  else if (frame==1)
	    Costhetastar = CosThetaStarHel(motherkstar, pion.particle, kaon.particle);

	  if (motherkstar.M() > 1.5)
            continue;

	  if (TMath::Abs(motherkstar.Rapidity())>=0.5)
	  continue;

	  fHisteventmult->Fill(5.5);	  

	  if (pion.charge * kaon.charge < 0)
	    {
	      fHisteventmult->Fill(6.5);	  
	      // Fill unlike histo
	      kstarUnlike->Fill(motherkstar.M(), lV0M, motherkstar.Pt(), Costhetastar);
	      //kstarUnlike->Fill(motherkstar.M(), lV0M, motherkstar.Pt());
	      
	    }
	  else if (pion.charge * kaon.charge > 0)
	    {
	      fHisteventmult->Fill(7.5);	  
	      // Fill like histo
	      kstarLike->Fill(motherkstar.M(), lV0M, motherkstar.Pt(), Costhetastar);
	      //kstarLike->Fill(motherkstar.M(), lV0M, motherkstar.Pt());
	    }
	  /*
	   if ((pion.charge > 0) && (kaon.charge > 0))
            {
              kstarposLike->Fill(motherkstar.M(), lV0M, motherkstar.Pt());
            }
          else if ((pion.charge < 0) && (kaon.charge < 0))
            { 
              kstarnegLike->Fill(motherkstar.M(), lV0M, motherkstar.Pt());
            }
	  */
	}
    }




 fHisteventsummary->Fill(3.5);

  ///////////////////////////////////////////////////////////////////
      
 Double_t Costhetastarmix=0.0;
 
 // Mixed event
 double pionmixenergy = 0.0;
 TLorentzVector pionmixVector(0.0,0.0,0.0,0.0);
 AliEventPool* pool = fPoolMgr->GetEventPool(lV0M, zv);
 if (pool && pool->IsReady())
   {
     //Int_t nMix = 10; // Set the number of mixed events
     Int_t nMix = 5; // Set the number of mixed events

     for (Int_t jMix = 0; jMix < nMix; jMix++)
       {
	 TObjArray* bgTracks = pool->GetRandomEvent();
	 Int_t numTracks = bgTracks->GetEntriesFast();

	 for (const auto& bgTrack : *bgTracks)
	   {
	     AliCompTrack* piontrackmix = dynamic_cast<AliCompTrack*>(bgTrack);
	     if (!piontrackmix)
	       {
		 AliFatal(Form("ERROR: Could not receive mix pool track %d\n", bgTrack->GetUniqueID()));
		 continue;
	       }

	     AliAODTrack* aodtrackmix = dynamic_cast<AliAODTrack*>(piontrackmix);

	     for (const auto& kaon : kaonCandidates)
	       {
	
		 if (piontrackmix->Charge() * kaon.charge > 0)
		   continue;

		 pionmixenergy = TMath::Sqrt(pionmass * pionmass + piontrackmix->Px() * piontrackmix->Px() + piontrackmix->Py() * piontrackmix->Py() + piontrackmix->Pz() * piontrackmix->Pz());		 
		 pionmixVector.SetPxPyPzE(piontrackmix->Px(), piontrackmix->Py(), piontrackmix->Pz(), pionmixenergy);
		 motherkstarmix.SetPxPyPzE(piontrackmix->Px()+kaon.particle.Px(), piontrackmix->Py()+kaon.particle.Py(), piontrackmix->Pz()+kaon.particle.Pz(), pionmixenergy+kaon.particle.E());
		 if (motherkstarmix.M() > 1.5)
		   continue;
		 if (TMath::Abs(motherkstarmix.Rapidity())>=0.5)
		 continue;

		 if (frame==0)
		   Costhetastarmix = CosThetaStar(motherkstarmix, pionmixVector, kaon.particle);
		 else if (frame==1)
		   Costhetastarmix = CosThetaStarHel(motherkstarmix, pionmixVector, kaon.particle);
		 // Fill mix histo
		 kstarMix->Fill(motherkstarmix.M(), lV0M, motherkstarmix.Pt(), Costhetastarmix);
		 //kstarMix->Fill(motherkstarmix.M(), lV0M, motherkstarmix.Pt());
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
void AliAnalysisTaskSA::Terminate(Option_t *)
{
  
}

//________________________________________________________________________

Bool_t AliAnalysisTaskSA::GoodEvent(const AliVVertex *vertex) 

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
  
  
  if( lAODevent->IsIncompleteDAQ() ) 
    {
      PostData(1, fOutput);
      return kFALSE;
    }

    AliAnalysisUtils *fUtils = new AliAnalysisUtils();
    /* if(fUtils->IsSPDClusterVsTrackletBG(fVevent))
    {
      PostData(1, fOutput);
      return kFALSE;
    }
    */

    fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE); 


  return kTRUE;

}
//-----------------------------------------------
Bool_t AliAnalysisTaskSA::IsPion(AliVTrack *vtrack)
{

  
  AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(vtrack);
  Double_t nsigmatpcpion=TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodtrack, AliPID::kPion));
  Double_t nsigmatofpion=TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodtrack, AliPID::kPion));
  Bool_t TOFHIT=kFALSE;

  TOFHIT = HasTOF(aodtrack);

  
  if(!TOFHIT)
    {
      if(nsigmatpcpion>2.0)return kFALSE;
    }

  else
    {
      if(nsigmatofpion>3.0)return kFALSE;
    }


  /*
  Double_t trkPt=aodtrack->Pt();
  Double_t nsigmacircularcut=TMath::Abs(sqrt(nsigmatpcpion*nsigmatpcpion + nsigmatofpion*nsigmatofpion));
  Double_t tofsig=0.0;

  if (trkPt<0.5)
    {
      if(nsigmatpcpion>3.0) return kFALSE;
    }
  else
    {
      //if(nsigmacircularcut>3.0) return kFALSE;
      if(nsigmatofpion>3.0) return kFALSE;                                                                                                  
    }
  */  
  return kTRUE;
}

Bool_t AliAnalysisTaskSA::IsKaon(AliVTrack *vtrack)
{

  AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(vtrack);
  Double_t nsigmatpckaon=TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodtrack, AliPID::kKaon));
  Double_t nsigmatofkaon=TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodtrack, AliPID::kKaon));
  Bool_t TOFHIT=kFALSE;

  TOFHIT = HasTOF(aodtrack);

  if(!TOFHIT)
    {
      if(nsigmatpckaon>2.0)return kFALSE;
    }
  else
    {
      if(nsigmatofkaon>3.0)return kFALSE;
    }

  /*
  Double_t trkPt=aodtrack->Pt();
  Double_t nsigmacircularcut=TMath::Abs(sqrt(nsigmatpckaon*nsigmatpckaon + nsigmatofkaon*nsigmatofkaon));
  Double_t tofsig=0.0;

  if (trkPt<0.45)
  {
    if(nsigmatpckaon>3.0) return kFALSE;
      }
  else
    {
      //if(nsigmacircularcut>3.0) return kFALSE;
      if(nsigmatofkaon>3.0) return kFALSE;                                                                                                  
    }
  */  
    
  return kTRUE;
}


Bool_t AliAnalysisTaskSA::HasTOF(AliAODTrack *track)
{
  bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len = track->GetIntegratedLength();
  bool hasTOF = hasTOFout && hasTOFtime && (len > 350.);
  return hasTOF;
}



//_________________________________________________________________________________________________
void AliAnalysisTaskSA::EventMixing()
{

  ////////////////////////////
  // Set-up for Mixed Event //
  ////////////////////////////
  const Int_t trackDepth = 1000000;
  const Int_t poolsize = 100;
  const Int_t nmix = 5;
  Double_t centralityBins[] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0};
  Double_t vertexBins[] = {-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0};
  const Int_t nCentralityBins = 10;
  const Int_t nVertexBins = 10;
  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centralityBins, nVertexBins, vertexBins);
  fPoolMgr->SetTargetValues(trackDepth,0.1,nmix);
}
//___________________________________________________________



Double_t AliAnalysisTaskSA::CosThetaStar(TLorentzVector mother, TLorentzVector daughter0, TLorentzVector daughter1)
{

  TVector3 BeamVect;
  BeamVect.SetXYZ(0,0,1);
  TVector3 momentumM(mother.Vect());
  TVector3 normal(mother.Y() / momentumM.Mag(), -mother.X() / momentumM.Mag(), 0.0);

  // Computes components of beta
  Double_t betaX = -mother.X() / mother.E();
  Double_t betaY = -mother.Y() / mother.E();
  Double_t betaZ = -mother.Z() / mother.E();

  // Computes Lorentz transformation of the momentum of the first daughter
  // into the rest frame of the mother and theta*
  daughter0.Boost(betaX, betaY, betaZ);
  TVector3 momentumD = daughter0.Vect();

  Double_t cosThetaStar = normal.Dot(momentumD) / momentumD.Mag();

  return cosThetaStar;

}




Double_t AliAnalysisTaskSA::CosThetaStarHel(TLorentzVector mother, TLorentzVector daughter0, TLorentzVector daughter1)

{

  // Computes components of beta
  Double_t betaX = -mother.X() / mother.E();
  Double_t betaY = -mother.Y() / mother.E();
  Double_t betaZ = -mother.Z() / mother.E();

  daughter0.Boost(betaX, betaY, betaZ);

  TVector3 zAxisHE = (mother.Vect()).Unit();
  Double_t thetaHE= zAxisHE.Dot((daughter0.Vect()).Unit());

  return thetaHE;

}


/*
Double_t AliAnalysisTaskSA::CosThetaStarEP(TLorentzVector mother, TLorentzVector daughter0, TLorentzVector daughter1)

{

  
  // Computes components of beta
  Double_t betaX = -mother.X() / mother.E();
  Double_t betaY = -mother.Y() / mother.E();
  Double_t betaZ = -mother.Z() / mother.E();

  daughter0.Boost(betaX, betaY, betaZ);

  TVector3 zAxisHE = (mother.Vect()).Unit();
  Double_t thetaHE= zAxisHE.Dot((daughter0.Vect()).Unit());

  return thetaHE;

}
*/





void AliAnalysisTaskSA::GetNUACorrectionHist(Int_t run, Int_t kParticleID)
{

  if(fListNUACorr){

    if(kParticleID==0){ //charge
      fHCorrectNUAposChrg = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Pos_Cent%d_Run%d",0,run)); 
      fHCorrectNUAnegChrg = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Neg_Cent%d_Run%d",0,run));
    }
    else{
      fHCorrectNUAposChrg=NULL;
      fHCorrectNUAnegChrg=NULL;
    }
    
  }//------> if list Exist

  // else {
  //   printf("\n ******** Warning: No NUA Correction File/List...!! \n Run= %d, Use NUA Wgt = 1.0 ********* \n",run);
  // }
}







void AliAnalysisTaskSA::GetV0MCorrectionHist(Int_t run, Int_t kParticleID)
{

  if(fListV0MCorr){
    //cout<<"*****************AT LAST I GOT INTO EVNTWGT CORRECTION FUNCTION******************************"<<endl;
    if(kParticleID==0){ //charge
      fHCorrectEVNTWGTChrg = (TH1F *) fListV0MCorr->FindObject(Form("hwgtCharge_Run%d",run)); 
    
    }
  }
  else{
    // cout<<"*****************AT LAST I GOT INTO EVNTWGT CORRECTION FUNCTION BUT LIST NOT PRESENT******************************"<<endl;
    fHCorrectEVNTWGTChrg=NULL;
  }
}
    



void AliAnalysisTaskSA::GetMCCorrectionHist(Int_t run,Float_t centr){

  if(fListTRKCorr) {
    //cout<<"\n =========> Info: Found TList with MC Tracking Corr Histograms <=========== "<<endl;
    /// Default: Centrality Independent MC efficiency:
    
    fHCorrectMCposChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgPos");
    fHCorrectMCnegChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgNeg");


    /// Centrality dependent MC efficiency: (Temporary)
    /*if(centr>5.0){
      fHCorrectMCposChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgPos");
      fHCorrectMCposPion =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyPionPos");
      fHCorrectMCposKaon =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyKaonPos");
      fHCorrectMCposProt =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyProtPos");

      fHCorrectMCnegChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgNeg");
      fHCorrectMCnegPion =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyPionNeg");
      fHCorrectMCnegKaon =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyKaonNeg");
      fHCorrectMCnegProt =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyProtNeg");
    }
    else{
      fHCorrectMCposChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgPosCent0");
      fHCorrectMCposPion =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyPionPosCent0");
      fHCorrectMCposKaon =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyKaonPosCent0");
      fHCorrectMCposProt =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyProtPosCent0");

      fHCorrectMCnegChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgNegCent0");
      fHCorrectMCnegPion =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyPionNegCent0");
      fHCorrectMCnegKaon =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyKaonNegCent0");
      fHCorrectMCnegProt =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyProtNegCent0");
    }
    */
    
    //for(int i=0;i<10;i++) {
    //fFB_Efficiency_Cent[i] = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_%d",i));
    //}
  }
  // else if(!fListTRKCorr){
  //   std::cout<<"\n\n !!!!**** Warning : No FB Efficiency Correction, run = "<<run<<"..!!!**** \n using MC TrkWgt = 1.0 \n"<<std::endl;
  // }
}










