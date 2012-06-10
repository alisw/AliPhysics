// $Id$
//
// Simulation EMCal task.
//
// Author: Saehanseul Oh

#include "AliAnalysisTaskSOH.h"

#include <TH1F.h>
#include <TH2F.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALTrack.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
#include "AliMCEvent.h"

ClassImp(AliAnalysisTaskSOH)

//________________________________________________________________________
AliAnalysisTaskSOH::AliAnalysisTaskSOH() :
  AliAnalysisTaskSE(), 
  fESD(0), 
  fMC(0), 
  fEsdTrackCuts(0x0),
  fHybridTrackCuts1(0x0),
  fHybridTrackCuts2(0x0),
  fTrackIndices(0x0),
  fOutputList(0x0),        
  fHEventStat(0),      
  fHTrkEffParGenPt(0), 
  fHTrkEffDetGenPt(0), 
  fHTrkEffDetRecPt(0), 
  fHEOverPVsPt(0x0),
  fHEMCalResponsePion(0x0), 
  fHEMCalResponseElec(0x0),
  fHEMCalResponseProton(0x0), 
  fHEMCalRecPhiEtaClus(0x0),
  fHEMCalRecPhiEtaTrk(0x0), 
  fHClsPhiEta(0x0),
  fHEMCalRecdPhidEta(0x0),
  fHEMCalRecdPhidEtaP(0x0),
  fHEMCalRecdPhidEtaM(0x0),
  fHEMCalRecdPhidEta_Truth(0x0),
  fHEMCalRecdPhidEtaP_Truth(0x0),
  fHEMCalRecdPhidEtaM_Truth(0x0),
  fHEMCalRecdPhidEtaposEta(0x0),
  fHEMCalRecdPhidEtanegEta(0x0)
{
  // Constructor

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}


//________________________________________________________________________
AliAnalysisTaskSOH::AliAnalysisTaskSOH(const char *name) :
  AliAnalysisTaskSE(name), 
  fESD(0), 
  fMC(0), 
  fEsdTrackCuts(0x0),
  fHybridTrackCuts1(0x0),
  fHybridTrackCuts2(0x0),
  fTrackIndices(0x0),
  fOutputList(0x0),        
  fHEventStat(0),      
  fHTrkEffParGenPt(0), 
  fHTrkEffDetGenPt(0), 
  fHTrkEffDetRecPt(0), 
  fHEOverPVsPt(0x0),
  fHEMCalResponsePion(0x0), 
  fHEMCalResponseElec(0x0),
  fHEMCalResponseProton(0x0), 
  fHEMCalRecPhiEtaClus(0x0),
  fHEMCalRecPhiEtaTrk(0x0), 
  fHClsPhiEta(0x0),
  fHEMCalRecdPhidEta(0x0),
  fHEMCalRecdPhidEtaP(0x0),
  fHEMCalRecdPhidEtaM(0x0),
  fHEMCalRecdPhidEta_Truth(0x0),
  fHEMCalRecdPhidEtaP_Truth(0x0),
  fHEMCalRecdPhidEtaM_Truth(0x0),
  fHEMCalRecdPhidEtaposEta(0x0),
  fHEMCalRecdPhidEtanegEta(0x0)
{
  // Constructor

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}


//________________________________________________________________________
AliAnalysisTaskSOH::~AliAnalysisTaskSOH()
{
  // Destructor.

  delete fEsdTrackCuts;
  delete fHybridTrackCuts1;
  delete fHybridTrackCuts2;
  delete fTrackIndices;
}


//________________________________________________________________________
void AliAnalysisTaskSOH::UserCreateOutputObjects()
{
  // Create histograms, called once.

  OpenFile(1);
  
  fOutputList = new TList();
  fOutputList->SetOwner(1);

  fHEventStat = new TH1F("fHEventStat","Event statistics for analysis",1,0,1);
  fHEventStat->GetXaxis()->SetBinLabel(1,"Event");
  fOutputList->Add(fHEventStat);

  fHTrkEffParGenPt = new TH1F("fHTrkEffParGenPt", "Particle level truth p_{T} distribution of generated primary charged pions;p_{T}^{gen} (GeV/c)",15,0.15,3.1);
  fOutputList->Add(fHTrkEffParGenPt);

  fHTrkEffDetGenPt = new TH1F("fHTrkEffDetGenPt", "Detector level truth p_{T} distribution of primary charged pions;p_{T}^{gen} (GeV/c)",15,0.15,3.1);
  fOutputList->Add(fHTrkEffDetGenPt);

  fHTrkEffDetRecPt = new TH1F("fHTrkEffDetRecPt", "Reconstructed track p_{T} distribution of primary charged pions;p_{T}^{rec} (GeV/c)",15,0.15,3.1);
  fOutputList->Add(fHTrkEffDetRecPt);

  fHEOverPVsPt = new TH2F("fHEOverPVsPt", "E/P vs track p_{T}; p_{T} (GeV/c); E/P", 20 , 0, 4, 40, 0, 3.2);
  fOutputList->Add(fHEOverPVsPt);

  fHEMCalResponsePion = new TH2F("fHEMCalResponsePion", "Pion E/P vs track p_{T}; p_{T} (GeV/c); E/P", 100 , 0, 4, 100, 0, 3.2);
  fOutputList->Add(fHEMCalResponsePion);

  fHEMCalResponseElec = new TH2F("fHEMCalResponseElec", "Electron E/P vs track p_{T};  p_{T} (GeV/c); E/P", 100 , 0, 4, 100, 0, 3.2);
  fOutputList->Add(fHEMCalResponseElec);

  fHEMCalResponseProton = new TH2F("fHEMCalResponseProton", "Proton E/P vs track p_{T};  p_{T} (GeV/c); E/P", 100 , 0, 4, 100, 0, 3.2);
  fOutputList->Add(fHEMCalResponseElec);

  fHEMCalRecPhiEtaClus = new TH2F("fHEMCalRecPhiEtaClus","EMCAL Cluster#phi-#eta; #eta; #phi",1000,-3,3,1000,-4,4);
  fOutputList->Add(fHEMCalRecPhiEtaClus);

  fHEMCalRecPhiEtaTrk = new TH2F("fHEMCalRecPhiEtaTrk","EMCAL Track #phi-#eta; #eta; #phi",1000,-3,3,1000,-4,4);
  fOutputList->Add(fHEMCalRecPhiEtaTrk);

  fHClsPhiEta = new TH2F("fHClsPhiEta", "cluster #phi-#eta; #eta; #phi",1000,-3,3,1000,-4,4);
  fOutputList->Add(fHClsPhiEta);

  fHEMCalRecdPhidEta = new TH2F("fHEMCalRecdPhidEta","EMCAL Cluster-Track #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEta);

  fHEMCalRecdPhidEtaP = new TH2F("fHEMCalRecdPhidEtaP","EMCAL Charge+ Cluster-Track #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEtaP);

  fHEMCalRecdPhidEtaM = new TH2F("fHEMCalRecdPhidEtaM","EMCAL Charge- Cluster-Track #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEtaM);

  fHEMCalRecdPhidEta_Truth = new TH2F("fHEMCalRecdPhidEta_Truth","EMCAL Cluster-Track(Truth matched) #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEta_Truth);

  fHEMCalRecdPhidEtaP_Truth = new TH2F("fHEMCalRecdPhidEtaP_Truth","EMCAL Charge+ Cluster-Track(Truth matched) #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEtaP_Truth);

  fHEMCalRecdPhidEtaM_Truth = new TH2F("fHEMCalRecdPhidEtaM_Truth","EMCAL Charge- Cluster-Track(Truth matched) #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEtaM_Truth);

  fHEMCalRecdPhidEtaposEta = new TH2F("fHEMCalRecdPhidEtaposEta","(+eta track) EMCAL Cluster-Track #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEtaposEta);

  fHEMCalRecdPhidEtanegEta = new TH2F("fHEMCalRecdPhidEtanegEta","(-eta track) EMCAL Cluster-Track #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEtanegEta);
  fTrackIndices = new TArrayI();

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskSOH::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  // Post output data.
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }

  fMC = MCEvent();
  if (!fMC) {
    printf("ERROR: fMC not available\n");
    return;
  }

  fHEventStat->Fill(0.5);

  if(fTrackIndices) 
    fTrackIndices->Reset();

  ProcessTrack();
  ProcessMc();
  
  PostData(1, fOutputList);
}   
   
//________________________________________________________________________
void  AliAnalysisTaskSOH::ProcessTrack()
{
  // Process track.

  fTrackIndices->Set(fESD->GetNumberOfTracks());
  AliDebug(3,Form("%s:%d Selecting tracks",(char*)__FILE__,__LINE__));

  Int_t isMth = 0;
  Int_t nTracks = 0;

  Float_t ClsPos[3] = {-999,-999,-999};
  Double_t emcTrkpos[3] = {-999,-999,-999};

  for(Int_t itr=0; itr<fESD->GetNumberOfCaloClusters(); itr++)
  {
    AliESDCaloCluster *cluster = fESD->GetCaloCluster(itr);
    if(!cluster) continue;
    cluster->GetPosition(ClsPos);
    TVector3 VClsPos(ClsPos[0], ClsPos[1], ClsPos[2]);
    fHClsPhiEta->Fill(VClsPos.Eta(), VClsPos.Phi());
  }

  for(Int_t itr=0; itr<fESD->GetNumberOfTracks(); itr++)
  {
    AliESDtrack *esdtrack = fESD->GetTrack(itr);
    if(!esdtrack)continue;
    AliESDtrack *newTrack = GetAcceptTrack(esdtrack);
    if(!newTrack) continue;
  
    Double_t clsE = -1;
    Int_t clsIndex = newTrack->GetEMCALcluster();
    if(newTrack->GetEMCALcluster()>-1)
    {
      AliESDCaloCluster *cluster = fESD->GetCaloCluster(clsIndex);
      if(IsGoodCluster(cluster))
      {
        isMth=1;
	      
        cluster->GetPosition(ClsPos);
        TVector3 VClsPos(ClsPos[0], ClsPos[1], ClsPos[2]);
	      
        AliEMCALTrack EMCTrk(*newTrack);
        if(!EMCTrk.PropagateToGlobal(ClsPos[0], ClsPos[1], ClsPos[2], 0.0, 0.0)) {
          continue;
        }
        EMCTrk.GetXYZ(emcTrkpos);
        TVector3 VemcTrkPos(emcTrkpos[0],emcTrkpos[1],emcTrkpos[2]);
	      
        fHEMCalRecPhiEtaClus->Fill(VClsPos.Eta(), VClsPos.Phi());
        fHEMCalRecPhiEtaTrk->Fill(VemcTrkPos.Eta(), VemcTrkPos.Phi());
	      
        Double_t dPhi = VClsPos.Phi() - VemcTrkPos.Phi();
        if (dPhi < -1*TMath::Pi()) dPhi += (2*TMath::Pi());
        else if (dPhi > TMath::Pi()) dPhi -= (2*TMath::Pi());
	      
        Double_t dEta = VClsPos.Eta() - VemcTrkPos.Eta();
	      
        fHEMCalRecdPhidEta->Fill(dEta, dPhi);	

        if((newTrack->GetLabel())>-1 && (newTrack->GetLabel()) < fMC->GetNumberOfTracks())
        {
          AliVParticle *vParticle = fMC->GetTrack(newTrack->GetLabel());
          if(IsGoodMcParticle(vParticle, newTrack->GetLabel()))
          {
            fHEMCalRecdPhidEta_Truth->Fill(dEta, dPhi);
            if(vParticle->Charge() > 0) fHEMCalRecdPhidEtaP_Truth->Fill(dEta, dPhi);
            if(vParticle->Charge() < 0) fHEMCalRecdPhidEtaM_Truth->Fill(dEta, dPhi);
          }
        }

        if(esdtrack->Charge() > 0) {fHEMCalRecdPhidEtaP->Fill(dEta, dPhi);}
        if(esdtrack->Charge() < 0) {fHEMCalRecdPhidEtaM->Fill(dEta, dPhi);}

        if(VemcTrkPos.Eta() > 0) fHEMCalRecdPhidEtaposEta->Fill(dEta, dPhi);
        if(VemcTrkPos.Eta() < 0) fHEMCalRecdPhidEtanegEta->Fill(dEta, dPhi);
	      
        clsE = cluster->E();
        if(newTrack->P()>0) fHEOverPVsPt->Fill(newTrack->Pt(),clsE/newTrack->P());
      }
	  
      Int_t ipart = newTrack->GetLabel();
      if(ipart>-1 && ipart<fMC->GetNumberOfTracks())
      {
        AliVParticle* vParticle = fMC->GetTrack(ipart);
        if(isMth && vParticle)
        {
          if(TMath::Abs(vParticle->PdgCode())==211)
          {
            fHEMCalResponsePion->Fill(newTrack->Pt(),clsE/newTrack->P());
          }
          if(TMath::Abs(vParticle->PdgCode())==11)
          {
            fHEMCalResponseElec->Fill(newTrack->Pt(),clsE/newTrack->P());
          }
          if(TMath::Abs(vParticle->PdgCode())==2212)
          {
            fHEMCalResponseProton->Fill(newTrack->Pt(),clsE/newTrack->P());
          }
        }
      }
    }
    fTrackIndices->AddAt(itr,nTracks);
    nTracks++;
  }
  fTrackIndices->Set(nTracks);
}

//________________________________________________________________________
void AliAnalysisTaskSOH::ProcessMc()
{
  // Process MC.

  for(Int_t ipart=0; ipart<fMC->GetNumberOfTracks(); ipart++)
  {
    AliVParticle* vParticle = fMC->GetTrack(ipart);
    if(!IsGoodMcParticle(vParticle, ipart)) continue;
    Int_t pdgCode = vParticle->PdgCode();
      
    //tracking effciency
    if(TMath::Abs(vParticle->Eta())<0.9 && TMath::Abs(pdgCode==211))
    {
      fHTrkEffParGenPt->Fill(vParticle->Pt());
      for(Int_t j=0; j<fTrackIndices->GetSize(); j++)
      {
        AliESDtrack *esdtrack = fESD->GetTrack(fTrackIndices->At(j));
        if(esdtrack && esdtrack->GetLabel()==ipart)
        {
          fHTrkEffDetGenPt->Fill(vParticle->Pt());
          fHTrkEffDetRecPt->Fill(esdtrack->Pt());
          break;
        }
      }
    }
  } 
}

//________________________________________________________________________
AliESDtrack *AliAnalysisTaskSOH::GetAcceptTrack(AliESDtrack *esdtrack)
{
  // Get accepted track.

  static AliESDtrack newTrack;
  if(fEsdTrackCuts->AcceptTrack(esdtrack))
  {
    esdtrack->Copy(newTrack);
    newTrack.SetTRDQuality(0);
  }
  else if(fHybridTrackCuts1->AcceptTrack(esdtrack))
  {
    if(esdtrack->GetConstrainedParam())
    {
      esdtrack->Copy(newTrack);
      const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
      newTrack.Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
      newTrack.SetTRDQuality(1);		
    }
    else 
      return 0x0;
  }
  else if(fHybridTrackCuts2->AcceptTrack(esdtrack))
  {
    if(esdtrack->GetConstrainedParam())
    {
      esdtrack->Copy(newTrack);
      const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
      newTrack.Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
      newTrack.SetTRDQuality(2);		
    }
    else 
      return 0x0;
  }
  else
  {
    return 0x0;
  }

  return &newTrack;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSOH::IsGoodMcParticle(AliVParticle* vParticle, Int_t ipart)
{
  // Return true if good MC particle.

  if(!vParticle) return kFALSE;
  if(!fMC->IsPhysicalPrimary(ipart)) return kFALSE;
  if (TMath::Abs(vParticle->Eta())>2) return kFALSE;
  if(vParticle->Pt()<0.15) return kFALSE;
  return kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskSOH::IsGoodCluster(AliESDCaloCluster *cluster)
{
  // Return true if good cluster.

  if(!cluster) return kFALSE;
  if (!cluster->IsEMCAL()) return kFALSE;
  return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskSOH::Terminate(Option_t *) 
{
  // Terminate analysis.
}
