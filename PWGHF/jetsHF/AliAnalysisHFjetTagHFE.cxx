// $Id$
//
// HF jet via HFE
// Author: S.Sakai

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliAODEvent.h"  // added
#include "AliAODHandler.h" // added 
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"

#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliAnalysisHFjetTagHFE.h"

ClassImp(AliAnalysisHFjetTagHFE)

//________________________________________________________________________
AliAnalysisHFjetTagHFE::AliAnalysisHFjetTagHFE() : 
  AliAnalysisTaskEmcalJet("AliAnalysisHFjetTagHFE", kTRUE),
  fVevent(0),
  //ftrack(0),
  fpidResponse(0),
  fHistTracksPt(0),
  fHistClustersPt(0),
  fHistLeadingJetPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsPtLeadHad(0),
  fHistJetsCorrPtArea(0),
  fHistPtDEtaDPhiTrackClus(0),
  fHistPtDEtaDPhiClusTrack(0),
  fHistClustDx(0),
  fHistClustDz(0),
  fHistTPCnSigma(0), //my
  fHistEop(0),
  fHistIncEle(0),
  fHistIncjet(0),
  fHistHFjet(0), 
  fInvmassULS(0), // my
  fInvmassLS(0), // my
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0)

{
  // Default constructor.

  fHistTracksPt       = new TH1*[fNcentBins];
  fHistClustersPt     = new TH1*[fNcentBins];
  fHistLeadingJetPt   = new TH1*[fNcentBins];
  fHistJetsPhiEta     = new TH2*[fNcentBins];
  fHistJetsPtArea     = new TH2*[fNcentBins];
  fHistJetsPtLeadHad  = new TH2*[fNcentBins];
  fHistJetsCorrPtArea = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistTracksPt[i] = 0;
    fHistClustersPt[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsPtLeadHad[i] = 0;
    fHistJetsCorrPtArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisHFjetTagHFE::AliAnalysisHFjetTagHFE(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistTracksPt(0),
  fHistClustersPt(0),
  fHistLeadingJetPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsPtLeadHad(0),
  fHistJetsCorrPtArea(0),
  fHistPtDEtaDPhiTrackClus(0),
  fHistPtDEtaDPhiClusTrack(0),
  fHistClustDx(0),
  fHistClustDz(0),
  fHistTPCnSigma(0),//my
  fHistEop(0),
  fHistIncEle(0),
  fHistIncjet(0),
  fHistHFjet(0),
  fInvmassULS(0),
  fInvmassLS(0),//my
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0)
{
  // Standard constructor.

  fHistTracksPt       = new TH1*[fNcentBins];
  fHistClustersPt     = new TH1*[fNcentBins];
  fHistLeadingJetPt   = new TH1*[fNcentBins];
  fHistJetsPhiEta     = new TH2*[fNcentBins];
  fHistJetsPtArea     = new TH2*[fNcentBins];
  fHistJetsPtLeadHad  = new TH2*[fNcentBins];
  fHistJetsCorrPtArea = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistTracksPt[i] = 0;
    fHistClustersPt[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsPtLeadHad[i] = 0;
    fHistJetsCorrPtArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisHFjetTagHFE::~AliAnalysisHFjetTagHFE()
{
  // Destructor.
  //delete ftrack;
  //delete fCaloClusters;
}

//________________________________________________________________________
void AliAnalysisHFjetTagHFE::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fJetsCont           = GetJetContainer(0);
  if(fJetsCont) { //get particles and clusters connected to jets
    fTracksCont       = fJetsCont->GetParticleContainer();
    fCaloClustersCont = fJetsCont->GetClusterContainer();
  } else {        //no jets, just analysis tracks and clusters
    fTracksCont       = GetParticleContainer(0);
    fCaloClustersCont = GetClusterContainer(0);
  }
  if(fTracksCont) fTracksCont->SetClassName("AliVTrack");
  if(fCaloClustersCont) fCaloClustersCont->SetClassName("AliVCluster");

  TString histname;

  for (Int_t i = 0; i < fNcentBins; i++) {
    if (fParticleCollArray.GetEntriesFast()>0) {
      histname = "fHistTracksPt_";
      histname += i;
      fHistTracksPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistTracksPt[i]->GetXaxis()->SetTitle("p_{T,track} (GeV/c)");
      fHistTracksPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistTracksPt[i]);
    }

    if (fClusterCollArray.GetEntriesFast()>0) {
      histname = "fHistClustersPt_";
      histname += i;
      fHistClustersPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistClustersPt[i]->GetXaxis()->SetTitle("p_{T,clus} (GeV/c)");
      fHistClustersPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistClustersPt[i]);
    }

    if (fJetCollArray.GetEntriesFast()>0) {
      histname = "fHistLeadingJetPt_";
      histname += i;
      fHistLeadingJetPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
      fHistLeadingJetPt[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistLeadingJetPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistLeadingJetPt[i]);
      
      histname = "fHistJetsPhiEta_";
      histname += i;
      fHistJetsPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 50, -1, 1, 101, 0, TMath::Pi()*2 + TMath::Pi()/200);
      fHistJetsPhiEta[i]->GetXaxis()->SetTitle("#eta");
      fHistJetsPhiEta[i]->GetYaxis()->SetTitle("#phi");
      fOutput->Add(fHistJetsPhiEta[i]);
      
      histname = "fHistJetsPtArea_";
      histname += i;
      fHistJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 30, 0, 3);
      fHistJetsPtArea[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistJetsPtArea[i]->GetYaxis()->SetTitle("area");
      fOutput->Add(fHistJetsPtArea[i]);

      histname = "fHistJetsPtLeadHad_";
      histname += i;
      fHistJetsPtLeadHad[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistJetsPtLeadHad[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistJetsPtLeadHad[i]->GetYaxis()->SetTitle("p_{T,lead} (GeV/c)");
      fHistJetsPtLeadHad[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetsPtLeadHad[i]);
    
      if (!(GetJetContainer()->GetRhoName().IsNull())) {
	histname = "fHistJetsCorrPtArea_";
	histname += i;
	fHistJetsCorrPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 30, 0, 3);
	fHistJetsCorrPtArea[i]->GetXaxis()->SetTitle("p_{T}^{corr} [GeV/c]");
	fHistJetsCorrPtArea[i]->GetYaxis()->SetTitle("area");
	fOutput->Add(fHistJetsCorrPtArea[i]);
      }
    }
  }

  histname = "fHistPtDEtaDPhiTrackClus";
  fHistPtDEtaDPhiTrackClus = new TH3F(histname.Data(),Form("%s;#it{p}_{T}^{track};#Delta#eta;#Delta#varphi",histname.Data()),100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiTrackClus);

  histname = "fHistPtDEtaDPhiClusTrack";
  fHistPtDEtaDPhiClusTrack = new TH3F(histname.Data(),Form("%s;#it{p}_{T}^{clus};#Delta#eta;#Delta#varphi",histname.Data()),100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiClusTrack);

  fHistClustDx = new TH1F("fHistClustDx","fHistClustDx;Dx",1000,0.,1.);
  fOutput->Add(fHistClustDx);

  fHistClustDz = new TH1F("fHistClustDz","fHistClustDz;Dz",1000,0.,1.);
  fOutput->Add(fHistClustDz);

  fHistTPCnSigma = new TH2F("fHistTPCnSigma","TPC nSigma;p_{T}(GeV/c);n#sigms",100,0.,20.,250,-5.,5.);
  fOutput->Add(fHistTPCnSigma);

  fHistEop = new TH2F("fHistEop","E/p;p_{T}(GeV/c);E/p",100,0.,20.,200,0.,4.);
  fOutput->Add(fHistEop);

  fHistIncEle = new TH1F("fHistIncEle","Inclusive electron;p_{T}",100,0.,20.);
  fOutput->Add(fHistIncEle);

  fHistIncjet = new TH2F("fHistIncjet","Inc jet;p_{T}",20,0,20,150,0.,150.);
  fOutput->Add(fHistIncjet);

  fHistHFjet = new TH2F("fHistHFjet","HF jet;p_{T}",20,0,20,150,0.,150.);
  fOutput->Add(fHistHFjet);

  fInvmassULS = new TH2F("fInvmassULS","ULS mass;p_{T};mass",20,0,20,150,0.,0.3);
  fOutput->Add(fInvmassULS);

  fInvmassLS = new TH2F("fInvmassLS","LS mass;p_{T};mass",20,0,20,150,0.,0.3);
  fOutput->Add(fInvmassLS);


  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisHFjetTagHFE::FillHistograms()
{
  // Fill histograms.

  if (fTracksCont) {
    AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle(0)); 
    while(track) {
      fHistTracksPt[fCentBin]->Fill(track->Pt()); 
      track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    }
  }

  if (fCaloClustersCont) {
    AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(0); 
    while(cluster) {
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);
      fHistClustersPt[fCentBin]->Fill(nPart.Pt());
      Double_t dx = cluster->GetTrackDx();
      Double_t dz = cluster->GetTrackDz();
      fHistClustDx->Fill(dx);
      fHistClustDz->Fill(dz);
      cluster = fCaloClustersCont->GetNextAcceptCluster();
    }
  }

  if (fJetsCont) {
    AliEmcalJet *jet = fJetsCont->GetNextAcceptJet(0); 
    while(jet) {

      fHistJetsPtArea[fCentBin]->Fill(jet->Pt(), jet->Area());
      fHistJetsPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi());

      Float_t ptLeading = fJetsCont->GetLeadingHadronPt(jet);
      fHistJetsPtLeadHad[fCentBin]->Fill(jet->Pt(), ptLeading);

      if (fHistJetsCorrPtArea[fCentBin]) {
	Float_t corrPt = jet->Pt() - fJetsCont->GetRhoVal() * jet->Area();
	fHistJetsCorrPtArea[fCentBin]->Fill(corrPt, jet->Area());
      }
      jet = fJetsCont->GetNextAcceptJet(); 
    }
    
    jet = fJetsCont->GetLeadingJet();
    if(jet) fHistLeadingJetPt[fCentBin]->Fill(jet->Pt());
  }

  CheckClusTrackMatching();

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisHFjetTagHFE::CheckClusTrackMatching()
{
  
  if(!fTracksCont || !fCaloClustersCont)
    return;

  Double_t deta = 999;
  Double_t dphi = 999;

  //Get closest cluster to track
  AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle(0)); 
  while(track) {
    //Get matched cluster
    Int_t emc1 = track->GetEMCALcluster();
    if(fCaloClustersCont && emc1>=0) {
      AliVCluster *clusMatch = fCaloClustersCont->GetCluster(emc1);
      if(clusMatch) {
	AliPicoTrack::GetEtaPhiDiff(track, clusMatch, dphi, deta);
	fHistPtDEtaDPhiTrackClus->Fill(track->Pt(),deta,dphi);
      }
    }
    track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
  }
  
  //Get closest track to cluster
  AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(0); 
  while(cluster) {
    TLorentzVector nPart;
    cluster->GetMomentum(nPart, fVertex);
    fHistClustersPt[fCentBin]->Fill(nPart.Pt());
    
    //Get matched track
    AliVTrack *mt = NULL;      
    AliAODCaloCluster *acl = dynamic_cast<AliAODCaloCluster*>(cluster);
    if(acl) {
      if(acl->GetNTracksMatched()>1)
	mt = static_cast<AliVTrack*>(acl->GetTrackMatched(0));
    }
    else {
      AliESDCaloCluster *ecl = dynamic_cast<AliESDCaloCluster*>(cluster);
      Int_t im = ecl->GetTrackMatchedIndex();
      if(fTracksCont && im>=0) {
	mt = static_cast<AliVTrack*>(fTracksCont->GetParticle(im));
      }
    }
    if(mt) {
      AliPicoTrack::GetEtaPhiDiff(mt, cluster, dphi, deta);
      fHistPtDEtaDPhiClusTrack->Fill(nPart.Pt(),deta,dphi);
      
      /* //debugging
	 if(mt->IsEMCAL()) {
	 Int_t emc1 = mt->GetEMCALcluster();
	 Printf("current id: %d  emc1: %d",fCaloClustersCont->GetCurrentID(),emc1);
	 AliVCluster *clm = fCaloClustersCont->GetCluster(emc1);
	 AliPicoTrack::GetEtaPhiDiff(mt, clm, dphi, deta);
	 Printf("deta: %f dphi: %f",deta,dphi);
	 }
      */
    }
    cluster = fCaloClustersCont->GetNextAcceptCluster();
  }
}

//________________________________________________________________________
void AliAnalysisHFjetTagHFE::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;

}

//________________________________________________________________________
Bool_t AliAnalysisHFjetTagHFE::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
  double Zvertex = pVtx->GetZ();  
  double Yvertex = pVtx->GetY();  
  double Xvertex = pVtx->GetX();  

  TClonesArray* ftrack = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
  //ftrack = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
  int ntracks = ftrack->GetEntries();

  TClonesArray* fCaloClusters = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters")); 

    ///////////////////
    //PID initialised//
    //////////////////
    //AliPIDResponse *fpidResponse = fInputHandler->GetPIDResponse();
    fpidResponse = fInputHandler->GetPIDResponse();

  //AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle(0)); 
  //while(track) {

    for (Int_t itrack = 0; itrack < ntracks; itrack++) {

        AliVParticle* ptrack = dynamic_cast<AliVTrack*>(ftrack->At(itrack));
        AliVTrack *track = dynamic_cast<AliVTrack*>(ptrack);
        AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(track);  // to apply cuts

        Bool_t isElectron = kFALSE;
        Bool_t fFlagNonHFE=kFALSE;
        //Bool_t fFlagPhotonicElec = kFALSE;
        Double_t epTarray[3];
        for(int i=0; i<3; i++)epTarray[i] = 0.0;

        // get track information
        Double_t pt = track->Pt(); 
        Double_t px = track->Px(); 
        Double_t py = track->Py(); 
        Double_t pz = track->Pz(); 
        Double_t eta = track->Eta(); 

        if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; // AOD track level
        if(pt<2.0)continue;
        if(fabs(eta>0.6))continue;
        if(track->GetTPCNcls() < 80) continue;
        if(atrack->GetITSNcls() < 3) continue;   // AOD track level
        if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) continue;

        // Get TPC nSigma
        Double_t dEdx =-999, fTPCnSigma=-999;
        dEdx = track->GetTPCsignal();
        fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);

        if(fTPCnSigma<-1 || fTPCnSigma>3)continue;
        fHistTPCnSigma->Fill(pt,fTPCnSigma);

        
        // Get E/p
        Int_t EMCalIndex = -1;
        EMCalIndex = track->GetEMCALcluster();
        if(EMCalIndex < 0) continue;
        
        AliVCluster *clustMatch=0x0;
        //if(!fUseTender) clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex);
        //clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex);
        clustMatch = dynamic_cast<AliVCluster*>(fCaloClusters->At(EMCalIndex));  // here
        
        if(clustMatch && clustMatch->IsEMCAL())
        {
            /////////////////////////////////////////////
            //Properties of tracks matched to the EMCAL//
            /////////////////////////////////////////////
            if(TMath::Abs(clustMatch->GetTrackDx())>0.05 || TMath::Abs(clustMatch->GetTrackDz())>0.05) continue;
            
            Double_t clustMatchE = clustMatch->E();
            
            //EMCAL EID info
            Double_t eop = -1.0;
            if(track->P()>0)eop = clustMatchE/track->P();
            cout << "eop = " << eop << endl;
            fHistEop->Fill(pt,eop);
            if(eop>0.85 && eop<1.3)isElectron = kTRUE;  
                 

            if(isElectron)
              {
               epTarray[0] = px;
               epTarray[1] = py;
               epTarray[2] = py;

               SelectPhotonicElectron(itrack, ftrack, track, fFlagNonHFE);

              }
          }

     if(isElectron)fHistIncEle->Fill(pt);

    // ++++++ find e in jet

    if (fJetsCont) 
       {
        AliEmcalJet *jet = fJetsCont->GetNextAcceptJet(0);  // full or charge ?
        while(jet) 
           {
            //jet->Pt() 
            //jet->Area()
            //jet->Eta() 
            //jet->Phi()
            //Float_t ptLeading = fJetsCont->GetLeadingHadronPt(jet);

            Bool_t iTagHFjet = tagHFjet( jet, epTarray, 0, pt);
	    Float_t corrPt = jet->Pt() - fJetsCont->GetRhoVal() * jet->Area();
            if(iTagHFjet)fHistIncjet->Fill(pt,corrPt);
            if(iTagHFjet && !fFlagNonHFE)fHistHFjet->Fill(pt,corrPt);
           }
         jet = fJetsCont->GetNextAcceptJet(); 
       }   
    
   }  // end of Track loop

  delete ftrack;
  delete fCaloClusters;

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.


}


Bool_t AliAnalysisHFjetTagHFE::tagHFjet(AliEmcalJet* jet, double *epT, int MCpid, double &maxpT_e)
{
  Bool_t HFjetTag = kFALSE;

  for (unsigned j = 0; j< jet->GetNumberOfTracks(); j++) 
      {
       AliVParticle *jetcont;
       jetcont = static_cast<AliVParticle*>(jet->TrackAt(j, fTracks));
       if(jetcont) continue;
       
       if(epT[0] == jetcont->Px() && epT[1] == jetcont->Py() && epT[2] == jetcont->Pz()) // electron in jet
         {
          HFjetTag = kTRUE;
          cout << "jet tag by HFE" << endl;
         }
     
      }

 return HFjetTag;
}


void AliAnalysisHFjetTagHFE::SelectPhotonicElectron(Int_t itrack, TClonesArray *ftrackCL, AliVTrack *track, Bool_t &fFlagPhotonicElec)
//void AliAnalysisHFjetTagHFE::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec)
{
    ///////////////////////////////////////////
    //////Non-HFE - Invariant mass method//////
    ///////////////////////////////////////////
    
    Bool_t flagPhotonicElec = kFALSE;
    
    Int_t ntracks = -999;
    //if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
    //ntracks = ftrack->GetEntries();
    ntracks = ftrackCL->GetEntries();
    
    for (Int_t jtrack = 0; jtrack < ntracks; jtrack++) {
        AliVParticle* VAssotrack = 0x0;
        //if(!fUseTender) VAssotrack  = fVevent->GetTrack(jtrack);
        //VAssotrack = dynamic_cast<AliVTrack*>(ftrack->At(jtrack)); //take tracks from Tender list
        VAssotrack = dynamic_cast<AliVTrack*>(ftrackCL->At(jtrack)); //take tracks from Tender list
        
        if (!VAssotrack) {
            printf("ERROR: Could not receive track %d\n", jtrack);
            continue;
        }

        AliVTrack *Assotrack = dynamic_cast<AliVTrack*>(VAssotrack);
        AliAODTrack *aAssotrack = dynamic_cast<AliAODTrack*>(VAssotrack);

        //------reject same track
        if(jtrack==itrack) continue;

        Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
        Double_t ptAsso=-999., nsigma=-999.0, mass=-999., width = -999;
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;

        nsigma = fpidResponse->NumberOfSigmasTPC(Assotrack, AliPID::kElectron);
        ptAsso = Assotrack->Pt();
        Int_t chargeAsso = Assotrack->Charge();
        Int_t charge = track->Charge();
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        //------track cuts applied
        if(!aAssotrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
        if(aAssotrack->GetTPCNcls() < 70) continue;
        if((!(aAssotrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aAssotrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        
        //-------loose cut on partner electron
        if(ptAsso <0.2) continue;
        if(aAssotrack->Eta()<-0.9 || aAssotrack->Eta()>0.9) continue;
        if(nsigma < -3 || nsigma > 3) continue;
        
        //-------define KFParticle to get mass
        AliKFParticle::SetField(fVevent->GetMagneticField());
        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*Assotrack, fPDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        //-------Get mass
        Int_t MassCorrect;
        MassCorrect = recg.GetMass(mass,width);

        if(fFlagLS)
            if(track->Pt()>1) fInvmassLS->Fill(track->Pt(),mass);
        if(fFlagULS)
            if(track->Pt()>1) fInvmassULS->Fill(track->Pt(),mass);
        
        if(mass<100 && fFlagULS && !flagPhotonicElec) flagPhotonicElec = kTRUE; //Tag Non-HFE (random mass cut, not optimised)
    }
    fFlagPhotonicElec = flagPhotonicElec;
}


//________________________________________________________________________
void AliAnalysisHFjetTagHFE::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
