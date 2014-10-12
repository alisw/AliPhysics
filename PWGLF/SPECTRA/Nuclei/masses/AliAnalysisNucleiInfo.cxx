#include "AliAnalysisNucleiInfo.h"

// ROOT includes
#include <TMath.h>
#include "TChain.h"

// AliRoot includes
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVEvent.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliCentrality.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TProfile.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisManager.h"
#include "TFile.h"

ClassImp(AliAnalysisNucleiInfo)

//_____________________________________________________________________________
AliAnalysisNucleiInfo::AliAnalysisNucleiInfo():
  AliAnalysisTaskSE(),
  FilterBit(16),                                
  EtaLimit(),                           
  DCAxyCut(1000.),                               
  DCAzCut(1000.),                                
  NsigmaTpcCut(2.0),
  iBconf(0),                           
  kTOF(0),
//iTriggerSel(-99),
  fAOD(NULL), 
  fESD(NULL),
  fEvent(NULL),
  fPIDResponse(NULL)
{
  EtaLimit[0]=-99.0;
  EtaLimit[1]=99.0;
  
  fList[0]=new TList();
  fList[0]->SetName("results");
  
  fList[1]=new TList();
  fList[1]->SetName("results2");
}
//______________________________________________________________________________
AliAnalysisNucleiInfo::AliAnalysisNucleiInfo(const char *name):
  AliAnalysisTaskSE(name),
  FilterBit(16),                                
  EtaLimit(),                           
  DCAxyCut(1000.),                               
  DCAzCut(1000.),                                
  NsigmaTpcCut(2.0),                           
  iBconf(0),                                  
  kTOF(0),
  //iTriggerSel(-99),
  fAOD(NULL), 
  fESD(NULL),
  fEvent(NULL),
  fPIDResponse(NULL)
{

  EtaLimit[0]=-99.0;
  EtaLimit[1]=99.0;

  fList[0]=new TList();
  DefineOutput(1, TList::Class());
  fList[0]->SetName("results");
  
  fList[1]=new TList();
  DefineOutput(2, TList::Class());
  fList[1]->SetName("results2");
}
//_____________________________________________________________________________
AliAnalysisNucleiInfo::~AliAnalysisNucleiInfo()
{
  if(fList[0]) delete fList[0];
  if(fList[1]) delete fList[1];
}
//______________________________________________________________________________
void AliAnalysisNucleiInfo::UserCreateOutputObjects()
{
  Char_t namePart[nPart][30];
  snprintf(namePart[0],30,"e");
  snprintf(namePart[1],30,"mu");
  snprintf(namePart[2],30,"pi");
  snprintf(namePart[3],30,"K");
  snprintf(namePart[4],30,"p");
  snprintf(namePart[5],30,"d");
  snprintf(namePart[6],30,"t");
  snprintf(namePart[7],30,"He3");
  snprintf(namePart[8],30,"He4");
  
  Char_t name[nSpec][30];
  snprintf(name[0],20,"e_plus");
  snprintf(name[1],20,"mu_plus");
  snprintf(name[2],20,"pi_plus");
  snprintf(name[3],20,"K_plus");
  snprintf(name[4],20,"p");
  snprintf(name[5],20,"d");
  snprintf(name[6],20,"t");
  snprintf(name[7],20,"He3");
  snprintf(name[8],20,"He4");
  snprintf(name[9],20,"e_minus");
  snprintf(name[10],20,"mu_plus");
  snprintf(name[11],20,"pi_plus");
  snprintf(name[12],20,"K_plus");
  snprintf(name[13],20,"p_bar");
  snprintf(name[14],20,"d_bar");
  snprintf(name[15],20,"t_bar");
  snprintf(name[16],20,"He3_bar");
  snprintf(name[17],20,"He4_bar");
  
  Int_t hbins[2];

  for(Int_t iB=0;iB<nBconf;iB++) {
    
    htemp[iB] = new TH1F("htemp","htemp (avoid the problem with the empty list...);B field",20,-10,10);

    //htriggerbits[iB] = new TH1I("htriggerbits","htriggerbits; bits",10,-5,5);
    htriggerbits[iB][0] = new TH1I("htriggerbits_0","trigger mask; bits",45,-5,40);
    htriggerbits[iB][1] = new TH1I("htriggerbits_1","trigger bits (exclusive); bits",45,-5,40);
    
    hZvertex[iB][0] = new TH1F("hZvertex_Selected","Vertex distribution of selected events;z vertex (cm)",240,-30,30);
    hZvertex[iB][1] = new TH1F("hZvertex_Analyzed","Vertex distribution of analyzed events;z vertex (cm)",240,-30,30);
    
    hEta[iB] = new TH1F("hEta_Analyzed","|#eta| distribution after the track cuts;#eta",200,-1.0,1.0);
    
    hPhi[iB] = new TH1F("hPhi_Analyzed","#phi distribution after the track cuts;#phi (rad.)",90,0,6.3);//Each TRD supermodule is divided for 5 (DeltaPhi(TRD)=0.35 theoretical)
    
    //hbins[0]=500; hbins[1]=2000;
    hbins[0]=500; hbins[1]=2000;//hbins[0]=100; hbins[1]=500
    fdEdxVSp[iB][0] = new TH2F("fdEdxVSp_pos","dE/dx vs p (positive charge); p_{TPC}/z (GeV/c); dE/dx_{TPC} (a.u.)",hbins[0],0,10,hbins[1],0,1000);
    fdEdxVSp[iB][1] = new TH2F("fdEdxVSp_neg","dE/dx vs p (negative charge); p_{TPC}/z (GeV/c); dE/dx_{TPC} (a.u.)",hbins[0],0,10,hbins[1],0,1000);

    Char_t name_hDeDxExp[nPart][200];
    Char_t title_hDeDxExp[nPart][200];
    for(Int_t i=0;i<nPart;i++) {
      snprintf(name_hDeDxExp[i],200,"hDeDxExp_%s",namePart[i]);
      snprintf(title_hDeDxExp[i],200,"Expected dE/dx of %s in the TPC;p_{TPC}/z (GeV/c);dE/dx_{TPC} (a.u.)",namePart[i]);
      hDeDxExp[iB][i] = new TProfile(name_hDeDxExp[i],title_hDeDxExp[i],200,0,10,0,1000,"");//,500,0,5,0,1000,""); toram
    }

    Char_t name_fNsigmaTpc[nSpec][200];
    Char_t title_fNsigmaTpc[nSpec][200];
    hbins[0]=200; hbins[1]=200;
    for(Int_t i=0;i<nSpec;i++) {
      snprintf(name_fNsigmaTpc[i],200,"NsigmaTpc_%s",name[i]);
      snprintf(title_fNsigmaTpc[i],200,"NsigmaTpc_%s;p_{TPC}/z (GeV/c);n_{#sigma_{TPC}}^{%s}",name[i],name[i]);
      fNsigmaTpc[iB][i] = new TH2F(name_fNsigmaTpc[i],title_fNsigmaTpc[i],hbins[0],0,10,hbins[1],-10,10);
    }
    
    hbins[0]=500; hbins[1]=520;//hbins[0]=200; hbins[1]=260;
    fBetaTofVSp[iB][0] = new TH2F("fBetaTofVSp_pos","#beta_{TOF} vs p/z (positive charge);p(GeV/c);#beta_{TOF}",hbins[0],0,10,hbins[1],0.4,1.05);
    fBetaTofVSp[iB][1] = new TH2F("fBetaTofVSp_neg","#beta_{TOF} vs p/z (negative charge);p(GeV/c);#beta_{TOF}",hbins[0],0,10,hbins[1],0.4,1.05);
    
    Char_t name_hBetaExp[nPart][200];
    Char_t title_hBetaExp[nPart][200];
    for(Int_t i=0;i<nPart;i++) {
      snprintf(name_hBetaExp[i],200,"hBetaTofVsP_Exp_%s",namePart[i]);
      snprintf(title_hBetaExp[i],200,"Expected #beta_{TOF} vs p/z of %s;p/z (GeV/c); #beta_{TOF}",namePart[i]);
      hBetaExp[iB][i] = new TProfile(name_hBetaExp[i],title_hBetaExp[i],200,0,10,0.4,1.05,"");
    }
    
    Char_t name_fNsigmaTof[nSpec][200];
    Char_t title_fNsigmaTof[nSpec][200];    
    hbins[0]=200; hbins[1]=200;
    for(Int_t i=0;i<nSpec;i++) {
      snprintf(name_fNsigmaTof[i],200,"NsigmaTof_%s",name[i]);
      snprintf(title_fNsigmaTof[i],200,"NsigmaTof_%s;p_{T}/z (GeV/c);n_{#sigma_{TOF}}^{%s}",name[i],name[i]);
      fNsigmaTof[iB][i] = new TH2F(name_fNsigmaTof[i],title_fNsigmaTof[i],hbins[0],0,10,hbins[1],-10,10);
    }

    Char_t name_fTofMinusExp[2][nSpec][200];
    Char_t title_fTofMinusExp[2][nSpec][200];    
    hbins[0]=200; hbins[1]=2000;
    
    for(Int_t i=0;i<nSpec;i++) {
      snprintf(name_fTofMinusExp[0][i],200,"TofMinusExp_%s",name[i]);
      snprintf(title_fTofMinusExp[0][i],200,"TofMinusExp_%s;p_{T}/z (GeV/c);tof-t_{exp}^{%s} (ps)",name[i],name[i]);
      
      snprintf(name_fTofMinusExp[1][i],200,"TofMinusExpWtpc_%s",name[i]);
      snprintf(title_fTofMinusExp[1][i],200,"TofMinusExpWtpc_%s (with a 2#sigma TPC cut);p_{T}/z (GeV/c);tof-t_{exp}^{%s} (ps)",name[i],name[i]);
    }
    for(Int_t it=0;it<2;it++) {
      for(Int_t i=0;i<nSpec;i++) {
	fTofMinusExp[iB][it][i] = new TH2F(name_fTofMinusExp[it][i],title_fTofMinusExp[it][i],hbins[0],0,10,hbins[1],-20000,20000);
      }
    }

    Char_t name_h2DCAap[18][200];
    Char_t title_h2DCAap[18][200];
    
    for(Int_t iS=0;iS<nSpec;iS++) {
      snprintf(name_h2DCAap[iS],200,"h2DCAap_%s",name[iS]);
      snprintf(title_h2DCAap[iS],200,"h2DCA_%s in for p_{T}/z<1.5GeV/c;DCA_{xy} (cm);DCA_{z} (cm)",name[iS]);
      if(iS==5 || iS==5+9) h2DCAap[iB][iS] = new TH2F(name_h2DCAap[iS],title_h2DCAap[iS],1750,-3.5,3.5,1750,-3.5,3.5);
      else h2DCAap[iB][iS] = new TH2F(name_h2DCAap[iS],title_h2DCAap[iS],1,-3.5,3.5,1,-3.5,3.5);
    }
    
    fList[iB]->Add(htemp[iB]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(htriggerbits[iB][i]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(hZvertex[iB][i]);
    fList[iB]->Add(hEta[iB]);
    fList[iB]->Add(hPhi[iB]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(fdEdxVSp[iB][i]);
    for(Int_t i=0;i<nPart;i++) {
      fList[iB]->Add(hDeDxExp[iB][i]);
    }    
    for(Int_t i=0;i<nSpec;i++) {
      fList[iB]->Add(fNsigmaTpc[iB][i]);
    }
    for(Int_t i=0;i<2;i++) fList[iB]->Add(fBetaTofVSp[iB][i]);
    for(Int_t i=0;i<nPart;i++) {
      fList[iB]->Add(hBetaExp[iB][i]);
    }
    for(Int_t i=0;i<nSpec;i++) {
      fList[iB]->Add(fNsigmaTof[iB][i]);
    }  
    for(Int_t i=0;i<nSpec;i++) {
      fList[iB]->Add(fTofMinusExp[iB][0][i]);
    } 
    for(Int_t i=0;i<nSpec;i++) {
      fList[iB]->Add(fTofMinusExp[iB][1][i]);
    } 
    for(Int_t iS=0;iS<nSpec;iS++) {
      if(iS==5 || iS==5+9) fList[iB]->Add(h2DCAap[iB][iS]);
    }
    
    // Post output data.
    PostData(1, fList[0]);
    PostData(2, fList[1]);
        
  }//end iB loop
}
//______________________________________________________________________________
void AliAnalysisNucleiInfo::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fAOD && !fESD){
    Printf("%s:%d AODEvent and ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    return;
  }
  
  if(fESD) fEvent = fESD;
  else fEvent = fAOD;

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse();
  
  //--------------------------Magnetic field polarity--------------------
  Double_t fBfield=fEvent->GetMagneticField();
  if(fBfield<0.0) iBconf=0;//B--
  else iBconf=1;//B++
  for(Int_t i=0;i<nBconf;i++) htemp[i]->Fill(fBfield);
    
  //-------------------------zVertex determination of event----------------
  Double_t zvtx = 9999.9;
  const AliVVertex* vtxEVENT = fEvent->GetPrimaryVertex();
  if(vtxEVENT->GetNContributors()>0) zvtx = vtxEVENT->GetZ();
  
  hZvertex[iBconf][0]->Fill(zvtx);
  
  //---------------------------EVENT CUTS-----------------------------
  if(TMath::Abs(zvtx) < 10.0){

    //TRIGGER SELECTION
    Int_t iTrigger=-2;

    if(inputHandler->IsEventSelected() & AliVEvent::kMB) iTrigger = 0;
    if(inputHandler->IsEventSelected() & AliVEvent::kCentral) iTrigger = 16;
    if(inputHandler->IsEventSelected() & AliVEvent::kSemiCentral) iTrigger = 17;
    //if((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()) & AliVEvent::kAny) iTrigger = 35;
    
    if(iTriggerSel!=-99) {//if a dedicated trigger is required
      if(iTrigger!=iTriggerSel) return;
    }
    
    for(Int_t i=0;i<32;i++) {
      Int_t bit=(1<<i);
      if(inputHandler->IsEventSelected() & bit) htriggerbits[iBconf][0]->Fill(i);
    }
    if(inputHandler->IsEventSelected() & AliVEvent::kAny) htriggerbits[iBconf][0]->Fill(35);
    if(inputHandler->IsEventSelected() & AliVEvent::kAnyINT) htriggerbits[iBconf][0]->Fill(36);
    
    htriggerbits[iBconf][1]->Fill(iTrigger);
    
    hZvertex[iBconf][1]->Fill(zvtx);
    
    Int_t nTracks = fEvent->GetNumberOfTracks();
    
    //----------------------loop on the TRACKS-----------------------------
    for(Int_t iT = 0; iT < nTracks; iT++) { 
      AliVTrack* track = (AliVTrack *) fEvent->GetTrack(iT);
      
      if (!track){
	continue;
      }
      
     //For the geometrical cuts
      Double_t eta = track->Eta();
      
      Bool_t trkFlag = 0;
      trkFlag = ((AliAODTrack *) track)->TestFilterBit(FilterBit);
      //TestFilterBit(16) -- Standard Cuts with very loose DCA: GetStandardITSTPCTrackCuts2011(kFALSE) && SetMaxDCAToVertexXY(2.4) && SetMaxDCAToVertexZ(3.2) && SetDCaToVertex2D(kTRUE)
      //TestFilterBit(32) (STARDARD) -- Standard Cuts with very tight DCA cut ( 7sigma^primaries: 7*(0.0015+0.0050/pt^1.1) ) : GetStandardITSTPCTrackCuts2011(). 
      
      //-------------------------------------start TRACK CUTS--------------------------------
      if ((track->Pt() < 0.2) || (eta<EtaLimit[0]) || (eta>EtaLimit[1]) || !trkFlag)
	continue;
        
       //Vertex determination
      Double_t b[2] = {-99., -99.};
      Double_t bCov[3] = {-99., -99., -99.};
      if (!track->PropagateToDCA(fEvent->GetPrimaryVertex(), fEvent->GetMagneticField(), 100., b, bCov))
	continue;
      
      Double_t DCAxy = b[0];
      Double_t DCAz = b[1];
      
      //Cut on the DCAxy
      Bool_t isDCAxyCut=kFALSE;
      if(TMath::Abs(DCAxy)<DCAxyCut) isDCAxyCut=kTRUE;
      
      //Cut on the DCAz
      Bool_t isDCAzCut=kFALSE;
      if(TMath::Abs(DCAz)<DCAzCut) isDCAzCut=kTRUE;

      if (!isDCAxyCut || !isDCAzCut)
	continue;

      //-------------------------------------end TRACK CUTS----------------------------------
     
      //-------------------------------------Track info--------------------------------------
      Double_t phi= track->Phi();
      Double_t charge = (Double_t)track->Charge();
      Double_t p = track->P();
      Double_t pt = track->Pt();
      Double_t dedx = track->GetTPCsignal();
      Double_t pTPC = track->GetTPCmomentum();
      Double_t tof  = track->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);
      Double_t beta = 0.0;
      //Double_t M2 = 999.9;
      //Double_t Z2 = 999.9;
      
      kTOF = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
      
      //-----------------------------TPC info------------------------------
      Double_t nsigmaTPC[nPart];
      Double_t expdedx[nPart];
      
      Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,mu,pi,K,p,d,t,He3,He4
      Int_t FlagPid = 0;
      
      for(Int_t iS=0;iS<9;iS++){
	nsigmaTPC[iS] = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) iS);
	//TPC identification:
	if(TMath::Abs(nsigmaTPC[iS])<NsigmaTpcCut) {
	  FlagPid += ((Int_t)TMath::Power(2,iS));
	}
      }
      
      hEta[iBconf]->Fill(eta);
      hPhi[iBconf]->Fill(phi);
      
      //More TPC info:
      for(Int_t iS=0;iS<9;iS++){
	expdedx[iS] = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, (AliPID::EParticleType) iS, AliTPCPIDResponse::kdEdxDefault, kTRUE);
	hDeDxExp[iBconf][iS]->Fill(pTPC,expdedx[iS]);
	nsigmaTPC[iS] = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) iS);
	//fNsigmaTpc[iBconf][iS]->Fill(pTPC,nsigmaTPC[iS]);
	if(charge>0) {//positive particle
	  fNsigmaTpc[iBconf][iS]->Fill(pTPC,nsigmaTPC[iS]);
	}
	else {//negative particle
	  fNsigmaTpc[iBconf][iS+nPart]->Fill(pTPC,nsigmaTPC[iS]);
	}
      }
          
      if(charge>0) fdEdxVSp[iBconf][0]->Fill(pTPC,dedx);
      else fdEdxVSp[iBconf][1]->Fill(pTPC,dedx);
      
      //ITS info
      for(Int_t iS=0;iS<9;iS++){
	if(FlagPid & stdFlagPid[iS]) {
	  if(pt<1.5) {
	    if(charge>0) {
	      h2DCAap[iBconf][iS]->Fill(DCAxy,DCAz);
	      h2DCAap[iBconf][iS]->Fill(-DCAxy,-DCAz);
	    }
	    else {
	      h2DCAap[iBconf][iS+nPart]->Fill(DCAxy,DCAz);
	      h2DCAap[iBconf][iS+nPart]->Fill(-DCAxy,-DCAz);
	    }
	  }
	}
      }
      
      //-----------------------------TOF info------------------------------
      
      Double_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.875612859,2.808921005,1.404195741,1.863689620};

      //----------------------------------------kTOF available-----------------------------
           
      if(kTOF) {
	
	Double_t exptimes[9];
	track->GetIntegratedTimes(exptimes);
	//Integrated times of the Nuclei:
	for(Int_t iN=5;iN<9;iN++) {
	  exptimes[iN] = exptimes[4]*exptimes[4]*(massOverZ[iN]*massOverZ[iN]/p/p+1)/(massOverZ[4]*massOverZ[4]/p/p+1);
	  exptimes[iN] = TMath::Sqrt(exptimes[iN]);
	}  
	
	beta=exptimes[0];
	beta=beta/tof;//beta = L/tof/c = t_e/tof
	
	Double_t nsigmaTOF[9];
	for(Int_t iS=0;iS<9;iS++){
	  nsigmaTOF[iS] = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType) iS);
	  if(charge>0) {
	    hBetaExp[iBconf][iS]->Fill(p,exptimes[0]/exptimes[iS]);
	    fNsigmaTof[iBconf][iS]->Fill(pt,nsigmaTOF[iS]);
	  }
	  else {
	    hBetaExp[iBconf][iS+nPart]->Fill(p,exptimes[0]/exptimes[iS]);
	    fNsigmaTof[iBconf][iS+nPart]->Fill(pt,nsigmaTOF[iS]);
	  }
	}
	
	if(charge>0) fBetaTofVSp[iBconf][0]->Fill(p,beta);
	else fBetaTofVSp[iBconf][1]->Fill(p,beta);
	
	for(Int_t iS=0;iS<9;iS++){
	  if(charge>0) {
	    fTofMinusExp[iBconf][0][iS]->Fill(pt,tof-exptimes[iS]);
	    if(FlagPid & stdFlagPid[iS]) fTofMinusExp[iBconf][1][iS]->Fill(pt,tof-exptimes[iS]);
	  }
	  else {
	    fTofMinusExp[iBconf][0][iS+nPart]->Fill(pt,tof-exptimes[iS]);
	    if(FlagPid & stdFlagPid[iS]) fTofMinusExp[iBconf][1][iS+nPart]->Fill(pt,tof-exptimes[iS]);
	  }
	}

      }//end kTOF available
    }//end track loop
  }//end loop on the events
}

//_____________________________________________________________________________
void AliAnalysisNucleiInfo::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
