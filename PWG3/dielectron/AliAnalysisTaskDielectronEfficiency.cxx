/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//#####################################################
//#                                                   #
//#  Simple efficiency study for dielectrons          #
//#  Author: Jens Wiechula Jens.Wiechula@cern.ch      #
//#                                                   #
//#####################################################

#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TROOT.h>
#include "TChain.h"
#include <TCanvas.h>
// #include "TTree.h"
#include <TH1.h>
#include <TH2F.h>
#include <THashList.h>

#include "AliAnalysisManager.h"

#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliKFParticle.h"
#include "AliESDtrackCuts.h"
#include "AliKineTrackCuts.h"
#include "AliLog.h"

#include "AliDielectronHistos.h"
#include "AliAnalysisTaskDielectronEfficiency.h"

ClassImp(AliAnalysisTaskDielectronEfficiency)


//=================================================================================
AliAnalysisTaskDielectronEfficiency::AliAnalysisTaskDielectronEfficiency() :
  AliAnalysisTask(),
  fInputEvent(0),
  fHist(0),
  fESDtrackCuts(0),
  fKineCutsLegs(0),
  fKineCutsMother(0),
  fIdMCMother(443),
  fIdMCDaughterP(-11),
  fIdMCDaughterN(11),
  fPDG(TDatabasePDG::Instance())
{
}

//=================================================================================
AliAnalysisTaskDielectronEfficiency::AliAnalysisTaskDielectronEfficiency(const char *name) :
  AliAnalysisTask(name,name),
  fInputEvent(0),
  fHist(0),
  fESDtrackCuts(new AliESDtrackCuts),
  fKineCutsLegs(new AliKineTrackCuts),
  fKineCutsMother(new AliKineTrackCuts),
  fIdMCMother(443),
  fIdMCDaughterP(-11),
  fIdMCDaughterN(11),
  fPDG(TDatabasePDG::Instance())
{
  //
  // named constructor. This is the one that should be used by the user, oterwise the
  // essential objects are not created!
  //
  DefineInput(0,TChain::Class());
  DefineOutput(0, THashList::Class());
}

//=================================================================================
AliAnalysisTaskDielectronEfficiency::~AliAnalysisTaskDielectronEfficiency()
{
  //
  // dtor
  //
  if (fESDtrackCuts)   delete fESDtrackCuts;
  if (fKineCutsLegs)   delete fKineCutsLegs;
  if (fKineCutsMother) delete fKineCutsMother;
}
//=================================================================================
void AliAnalysisTaskDielectronEfficiency::CreateOutputObjects() {
  //
  // Create histograms
  // Called once
  //

  //-------------------
  // MC truth produced
  fHist=new AliDielectronHistos;
  fHist->AddClass("MC;MCcut;DataSameMother;Event;DataCuts;DataTRDCuts");

  fHist->UserHistogram("MC",    "JpsiMCPt"  ,"MC Jpsi;Pt [GeV]"         ,100,0,10);
  fHist->UserHistogram("MC",    "mass"    ,"MC Jpsi; Inv.Mass [GeV]" ,100,0,4);
  fHist->UserHistogram("MC",    "e+Pt"    ,"MC e+ from JPsi;Pt [GeV]" ,100,0,10);
  fHist->UserHistogram("MC",    "e-Pt"    ,"MC e- from JPsi;Pt [GeV]" ,100,0,10);
  fHist->UserHistogram("MC", "dndyPt"      ,"MC Jpsi procution; Rapidity;Pt",100,-4,4,100,0,10);
  fHist->UserHistogram("MC","dndy"        ,"MC dNdy Jpsi;Rapidity;Entries/event"    ,100,-4,4);
  fHist->GetHistogram("MC","dndy")->Sumw2();

  //-------------------
  //MC truth after cuts
  fHist->UserHistogram("MCcut",    "JpsiMCPt"  ,"MC Jpsi;Pt [GeV]"         ,100,0,10);
  fHist->UserHistogram("MCcut",    "mass"    ,"MC Jpsi; Inv.Mass [GeV]" ,100,0,4);
  fHist->UserHistogram("MCcut",    "e+Pt"    ,"MC e+ from JPsi;Pt [GeV]" ,100,0,10);
  fHist->UserHistogram("MCcut",    "e-Pt"    ,"MC e- from JPsi;Pt [GeV]" ,100,0,10);
  fHist->UserHistogram("MCcut", "dndyPt"      ,"MC Jpsi procution; Rapidity;Pt",100,-4,4,100,0,10);
  fHist->UserHistogram("MCcut","dndy"        ,"MC dNdy Jpsi;Rapidity;Entries/event"    ,100,-4,4);
  fHist->GetHistogram("MCcut","dndy")->Sumw2();
  
  //-----------------
  //reconstructed data with cuts on MC truth
  fHist->UserHistogram("DataSameMother","JpsiMCPt","Rec Jpsi; MC Pt [GeV]", 100,0,10);
  fHist->UserHistogram("DataSameMother","dndyPtMC" ,"Rec Jpsi procution; Rapidity;Pt",100,-4,4,100,0,10);
  fHist->UserHistogram("DataSameMother", "e+Pt"   ,"Rec e+ from JPsi;MC Pt [GeV]" ,100,0,10);
  fHist->UserHistogram("DataSameMother", "e-Pt"   ,"Rec e- from JPsi;MC Pt [GeV]" ,100,0,10);
  fHist->UserHistogram("DataSameMother","dndy"    ,"Rec Jpsi;Rapidity;Entries/event"             ,100,-4,4);
  fHist->GetHistogram("DataSameMother","dndy")->Sumw2();
  
  fHist->UserHistogram("DataSameMother","mass"    ,"Rec Jpsi (KF); Inv.Mass [GeV]" ,100,0,4);
  fHist->UserHistogram("DataSameMother","JpsiPt"  ,"Rec Jpsi (KF); Pt [GeV]"       ,100,0,10);
  fHist->UserHistogram("DataSameMother","Chi2"    ,"Rec Jpsi (KF); #Chi^{2}"       ,100,0,50);
  fHist->UserHistogram("DataSameMother","dndyPt"   ,"Rec Jpsi procution (KF); Rapidity;Pt",100,-4,4,100,0,10);
  //------------------
  // reconstructed data after ESD track cuts and cuts on MC truth
  //------------------
  fHist->UserHistogram("DataCuts","JpsiMCPt","Rec Jpsi; MC Pt [GeV]", 100,0,10);
  fHist->UserHistogram("DataCuts","dndyPtMC" ,"Rec Jpsi procution; Rapidity;Pt",100,-4,4,100,0,10);
  fHist->UserHistogram("DataCuts", "e+Pt"   ,"Rec e+ from JPsi;MC Pt [GeV]" ,100,0,10);
  fHist->UserHistogram("DataCuts", "e-Pt"   ,"Rec e- from JPsi;MC Pt [GeV]" ,100,0,10);
  fHist->UserHistogram("DataCuts","dndy"    ,"Rec Jpsi;Rapidity;Entries/event"             ,100,-4,4);
  fHist->GetHistogram("DataCuts","dndy")->Sumw2();
  
  fHist->UserHistogram("DataCuts","mass"    ,"Rec Jpsi (KF); Inv.Mass [GeV]" ,100,0,4);
  fHist->UserHistogram("DataCuts","JpsiPt"  ,"Rec Jpsi (KF); Pt [GeV]"       ,100,0,10);
  fHist->UserHistogram("DataCuts","Chi2"    ,"Rec Jpsi (KF); #Chi^{2}"       ,100,0,50);
  fHist->UserHistogram("DataCuts","dndyPt"   ,"Rec Jpsi procution (KF); Rapidity;Pt",100,-4,4,100,0,10);
  //------------------
  // after ESD track cuts + TRD cuts + MC cuts
  //------------------
  fHist->UserHistogram("DataTRDCuts","JpsiMCPt","Rec Jpsi; MC Pt [GeV]", 100,0,10);
  fHist->UserHistogram("DataTRDCuts","dndyPtMC" ,"Rec Jpsi procution; Rapidity;Pt",100,-4,4,100,0,10);
  fHist->UserHistogram("DataTRDCuts", "e+Pt"   ,"Rec e+ from JPsi;MC Pt [GeV]" ,100,0,10);
  fHist->UserHistogram("DataTRDCuts", "e-Pt"   ,"Rec e- from JPsi;MC Pt [GeV]" ,100,0,10);
  fHist->UserHistogram("DataTRDCuts","dndy"    ,"Rec Jpsi;Rapidity;Entries/event"             ,100,-4,4);
  fHist->GetHistogram("DataTRDCuts","dndy")->Sumw2();
  
  fHist->UserHistogram("DataTRDCuts","mass"    ,"Rec Jpsi (KF); Inv.Mass [GeV]" ,100,0,4);
  fHist->UserHistogram("DataTRDCuts","JpsiPt"  ,"Rec Jpsi (KF); Pt [GeV]"       ,100,0,10);
  fHist->UserHistogram("DataTRDCuts","Chi2"    ,"Rec Jpsi (KF); #Chi^{2}"       ,100,0,50);
  fHist->UserHistogram("DataTRDCuts","dndyPt"   ,"Rec Jpsi procution (KF); Rapidity;Pt",100,-4,4,100,0,10);

  //-----------------
  //Event information
  //-----------------
  fHist->UserHistogram("Event","NEvents","Number of events",1,0,1);

}

// //____________________________________________________________
void AliAnalysisTaskDielectronEfficiency::ConnectInputData(Option_t *) {
  //
  // connect the input data
  //
  fInputEvent=0;
  TTree* tree=dynamic_cast<TTree*>(GetInputData(0));
  if (!tree) {
    printf("ERROR: Could not read chain from input slot 0\n");
  } else {
    AliInputEventHandler *eventH = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!eventH) {
      AliError("Could not get ESDInputHandler");
    } else {
      fInputEvent = eventH->GetEvent();
      AliInfo("*** CONNECTED NEW EVENT ****");
    }
  }
}

//=================================================================================
void AliAnalysisTaskDielectronEfficiency::Exec(Option_t *) {
  //
  // Main loop. Called for every event
  // Process the event in FillPlots and post the data afterwards
  //
  if (!fInputEvent) {
    Printf("ERROR: Could not get input event\n");
    return;
  }

  FillPlots(fInputEvent);
  
  PostData(0, const_cast<THashList*>(fHist->GetHistogramList()));
}


//====================================================================================
void AliAnalysisTaskDielectronEfficiency::FillPlots(AliVEvent *event)
{
  //
  // Fill histograms
  //
  AliESDEvent *esd=dynamic_cast<AliESDEvent*>(event);
  if (!esd) return;
  Int_t ntrack=esd->GetNumberOfTracks();

  // Fetch Stack 
  AliMCEventHandler *mcH = (AliMCEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
  if(!mcH) {
    AliError("No MC handler found\n")
    return;
  }
  
  AliMCEvent *mcev=mcH->MCEvent();
  if (!mcev){
    AliError("No MC event found\n")
    return;
  }
  AliStack*  pStack = mcev->Stack();
  
  if (!pStack) return;

  //fill event info
  fHist->Fill("Event","NEvents",0);

  //fill MC histograms
  FillMCInfo(pStack);

  //
  Float_t massMother=0;
  if (fIdMCMother>-1) fPDG->GetParticle(fIdMCMother)->Mass();
  
  TLorentzVector v;
  //loop over all tracks
  for (Int_t itrack=0; itrack<ntrack; ++itrack){
    //negative particles only in this loop
    AliESDtrack *trackN=esd->GetTrack(itrack);
    if (trackN->Charge()!=-1) continue;

    //MC truth
    Int_t labelN=TMath::Abs(trackN->GetLabel());
    if (labelN<0) continue;
    TParticle *pN=pStack->Particle(labelN);
    Int_t pdgN=pN->GetPdgCode();

    //MC mother
    Int_t idMotherN=pN->GetFirstMother();
    TParticle *motherN=0;
    Int_t pdgMotherN=-1;
    if (fIdMCMother>-1&&idMotherN>-1){
      motherN=pStack->Particle(idMotherN);
      pdgMotherN=motherN->GetPdgCode();
    }
    
    for (Int_t itrack2=0; itrack2<ntrack; ++itrack2){
      //positive particles only in this loop
      AliESDtrack *trackP=esd->GetTrack(itrack2);
      if (trackP->Charge()!=1) continue;
      
      //MC truth
      Int_t labelP=TMath::Abs(trackP->GetLabel());
      if (labelP<0) continue;
      TParticle *pP=pStack->Particle(labelP);
      //       Int_t pdgP=pP->GetPdgCode();

      //MC mother
      Int_t idMotherP=pP->GetFirstMother();
//       TParticle *motherP=0;
      //       Int_t pdgMotherP=0;
      if (idMotherP>-1){
//         motherP=pStack->Particle(idMotherP);
      //       pdgMotherP=motherP->GetPdgCode();
      }
      //===============
      //Fill histograms
      //===============
      Bool_t motherOK=kFALSE;
      if (fIdMCMother==-1) motherOK=kTRUE;
      else if (pdgMotherN==fIdMCMother) motherOK=kTRUE;

      if (pdgN==fIdMCDaughterN && motherOK){ //electron and mother is fIdMCMother
        AliKFParticle electron(*trackN,11);
        AliKFParticle positron(*trackP,-11);
        AliKFParticle jpsi(electron);
        jpsi+=positron;

        Bool_t sameMother=kFALSE;
        Bool_t motherCutOK=kFALSE;

        if (fIdMCMother==-1) {
          //accept as same mother if we don't requite
          sameMother=kTRUE;
          motherCutOK=kTRUE;
        }
        else
        {
          if (idMotherN==idMotherP) sameMother=kTRUE;
          if (fKineCutsMother->IsSelected(motherN)) motherCutOK=kTRUE;
        }

        if ( sameMother  &&  // same mother
             motherCutOK && //cuts mother MC truth
            fKineCutsLegs->IsSelected(pN) && fKineCutsLegs->IsSelected(pP)){ //cuts legs MC truth

          //MC only data
          if (motherN){
            fHist->Fill("DataSameMother","JpsiMCPt",motherN->Pt());
            fHist->Fill("DataSameMother","e+Pt",pP->Pt());
            fHist->Fill("DataSameMother","e-Pt",pN->Pt());
            v.SetPxPyPzE(motherN->Px(),motherN->Py(),motherN->Pz(),motherN->Energy());
            fHist->Fill("DataSameMother","dndy",v.Rapidity());
            fHist->Fill("DataSameMother","dndyPtMC",v.Rapidity(),motherN->Pt());
          }

          //reconstructed data
          v.SetPtEtaPhiM(jpsi.GetPt(),jpsi.GetEta(),jpsi.GetPhi(),massMother);
//           printf("Jpsi: %f,%f,%f,%f\n",jpsi.GetPt(),jpsi.GetEta(),jpsi.GetPhi(),massMother);
          fHist->Fill("DataSameMother","JpsiPt",jpsi.GetPt());
          fHist->Fill("DataSameMother","Chi2",jpsi.GetChi2()/jpsi.GetNDF());
          fHist->Fill("DataSameMother","mass",jpsi.GetMass());
          fHist->Fill("DataSameMother","dndyPt",v.Rapidity(),jpsi.GetPt());

          //histograms after ESD cuts
          if (fESDtrackCuts->IsSelected(trackN)&&fESDtrackCuts->IsSelected(trackP)){
            //MC only
            if (motherN) {
              fHist->Fill("DataCuts","JpsiMCPt",motherN->Pt());
              fHist->Fill("DataCuts","e+Pt",pP->Pt());
              fHist->Fill("DataCuts","e-Pt",pN->Pt());
              v.SetPxPyPzE(motherN->Px(),motherN->Py(),motherN->Pz(),motherN->Energy());
              fHist->Fill("DataCuts","dndy",v.Rapidity());
              fHist->Fill("DataCuts","dndyPtMC",motherN->Eta(),motherN->Pt());
            }

          //reconstructed data
            v.SetPtEtaPhiM(jpsi.GetPt(),jpsi.GetEta(),jpsi.GetPhi(),massMother);
            fHist->Fill("DataCuts","JpsiPt",jpsi.GetPt());
            fHist->Fill("DataCuts","Chi2",jpsi.GetChi2()/jpsi.GetNDF());
            fHist->Fill("DataCuts","mass",jpsi.GetMass());
            fHist->Fill("DataCuts","dndyPt",v.Rapidity(),jpsi.GetPt());

            //Additional TRD cuts
            if ( ((trackN->GetStatus()&AliESDtrack::kTRDrefit)!=0) && trackN->GetTRDntrackletsPID()>4 ){
              if ( ((trackP->GetStatus()&AliESDtrack::kTRDrefit)!=0) && trackP->GetTRDntrackletsPID()>4 ){
                if (motherN){
                  fHist->Fill("DataTRDCuts","JpsiMCPt",motherN->Pt());
                  fHist->Fill("DataTRDCuts","e+Pt",pP->Pt());
                  fHist->Fill("DataTRDCuts","e-Pt",pN->Pt());
                  v.SetPxPyPzE(motherN->Px(),motherN->Py(),motherN->Pz(),motherN->Energy());
                  fHist->Fill("DataTRDCuts","dndy",v.Rapidity());
                  fHist->Fill("DataTRDCuts","dndyPtMC",v.Rapidity(),motherN->Pt());
                }
                
                //reconstructed data
                v.SetPtEtaPhiM(jpsi.GetPt(),jpsi.GetEta(),jpsi.GetPhi(),massMother);
                fHist->Fill("DataTRDCuts","JpsiPt",jpsi.GetPt());
                fHist->Fill("DataTRDCuts","Chi2",jpsi.GetChi2()/jpsi.GetNDF());
                fHist->Fill("DataTRDCuts","mass",jpsi.GetMass());
                fHist->Fill("DataTRDCuts","dndyPt",v.Rapidity(),jpsi.GetPt());
              }
            }
          }
        }
      }
    }
  }
}

//____________________________________________________________
// Int_t AliAnalysisTaskDielectronEfficiency::Merge(TList *list)
// {
//   //
//   // Merge function
//   //
//   if (!list) return 0;
//   if (list->IsEmpty()) return 1;
// 
//   TIter next(list);
//   while ( (TObject *o=next()) ){
//     AliAnalysisTaskDielectronEfficiency *task=dynamic_cast<AliAnalysisTaskDielectronEfficiency*>o;
//     if (!o) continue;
//     fNev+=task->fNev;
//   }
// }


void AliAnalysisTaskDielectronEfficiency::FillMCInfo(AliStack * const pStack)
{
  //
  // fill pure MC histograms
  //
  
  TLorentzVector v;
  //Fill MC info
  for (Int_t ipart=0; ipart<pStack->GetNtrack(); ++ipart){
    TParticle *part=pStack->Particle(ipart);
//     printf("Particle %d\n",part->GetPdgCode());
    if (part->GetPdgCode()!=fIdMCMother || part->GetNDaughters()!=2) continue;
    TParticle *d1=pStack->Particle(part->GetFirstDaughter());
    TParticle *d2=pStack->Particle(part->GetLastDaughter());
    TParticle *dP=0;
    TParticle *dN=0;
    if (fPDG->GetParticle(d1->GetPdgCode())->Charge()>0){
      dP=d1;
      dN=d2;
    }else{
      dP=d2;
      dN=d1;
    }
    if ( dP->GetPdgCode()!=fIdMCDaughterP || dN->GetPdgCode()!=fIdMCDaughterN ) continue;
    v.SetPxPyPzE(part->Px(),part->Py(),part->Pz(),part->Energy());
    fHist->Fill("MC","JpsiMCPt",part->Pt());
    fHist->Fill("MC","dndy",v.Rapidity());
    fHist->Fill("MC","dndyPt",v.Rapidity(),part->Pt());
    fHist->Fill("MC","e-Pt",dN->Pt());
    fHist->Fill("MC","e+Pt",dP->Pt());
    //e+ e- inv mass
    TLorentzVector vE;
    vE.SetPxPyPzE(dN->Px(),dN->Py(),dN->Pz(),dN->Energy());
    TLorentzVector vP;
    vP.SetPxPyPzE(dP->Px(),dP->Py(),dP->Pz(),dP->Energy());
    fHist->Fill("MC","mass",(vE+vP).M());


    //cuts
    if (!fKineCutsMother->IsSelected(part)) continue;
    if (!fKineCutsLegs->IsSelected(d1) || !fKineCutsLegs->IsSelected(d2) ) continue;
      
    v.SetPxPyPzE(part->Px(),part->Py(),part->Pz(),part->Energy());
    fHist->Fill("MCcut","JpsiMCPt",part->Pt());
    fHist->Fill("MCcut","dndy",v.Rapidity());
    fHist->Fill("MCcut","dndyPt",v.Rapidity(),part->Pt());
    fHist->Fill("MCcut","e-Pt",dN->Pt());
    fHist->Fill("MCcut","e+Pt",dP->Pt());
  }
}
void AliAnalysisTaskDielectronEfficiency::SetupDefaultCuts(Int_t type)
{
  //
  // setup standard ESD track cuts
  //

  if (type==0){
    //ESD cuts
    fESDtrackCuts->SetMaxDCAToVertexZ(3.0);
    fESDtrackCuts->SetMaxDCAToVertexXY(3.0);
    fESDtrackCuts->SetRequireTPCRefit(kTRUE);
    fESDtrackCuts->SetRequireITSRefit(kTRUE);
    fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    
    fESDtrackCuts->SetMinNClustersTPC(50);
    fESDtrackCuts->SetMaxChi2PerClusterTPC(4);

    //MC cuts
    fKineCutsLegs->SetEtaRange(-0.9,0.9);
    fKineCutsMother->SetRapRange(-0.9,0.9);
  } else if (type==1) {
//     fESDtrackCuts->SetMaxCovDiagonalElements(2, 2, .5, .5, 2);
    fESDtrackCuts->SetMaxDCAToVertexZ(3.0);
    fESDtrackCuts->SetMaxDCAToVertexXY(3.0);
    fESDtrackCuts->SetRequireTPCRefit(kTRUE);
    fESDtrackCuts->SetRequireITSRefit(kTRUE);
    fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    
    fESDtrackCuts->SetMinNClustersTPC(50);
    fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
    
    //MC cuts
//     fKineCutsLegs->SetEtaRange(-0.9,0.9);
//     fKineCutsMother->SetRapRange(-0.9,0.9);
    
  }
}

//===================================================================================
void AliAnalysisTaskDielectronEfficiency::Terminate(Option_t *) {
  //
  // Called once at the end of the query
  //

  
  AliDielectronHistos *hist=new AliDielectronHistos;
  hist->SetHistogramList(*(THashList*)GetOutputData(0));
  
  if (hist->GetHistogram("Event","NEvents")){
    //get number of events
    Double_t nev=hist->GetHistogram("Event","NEvents")->GetBinContent(1);
    //
    //normalise dndy histograms
    //
    hist->GetHistogram("DataSameMother","dndy")->Scale(1./nev);
    hist->GetHistogram("MC","dndy")->Scale(1./nev);
    hist->GetHistogram("DataCuts","dndy")->Scale(1./nev);
    hist->GetHistogram("DataTRDCuts","dndy")->Scale(1./nev);

    //
    // create the efficiency histograms
    //
    // hEffTracking/2D only tracking effects, no esd cuts
    // hEffESDCuts  tracking plus ESD cuts
    // hEffTRDCuts  tracking plus ESD plus TRD cuts
    //
    TH1F *hEffTracking=(TH1F*)hist->GetHistogram("DataSameMother","JpsiMCPt")->Clone("Efficiency");
    hEffTracking->Divide(hist->GetHistogram("MCcut","JpsiMCPt"));
    hEffTracking->SetTitle("Efficiencies");
    hist->UserHistogram("DataSameMother",hEffTracking);
    
    TH1F *hEffESDCuts=(TH1F*)hist->GetHistogram("DataCuts","JpsiMCPt")->Clone("Efficiency");
    hEffESDCuts->Divide(hist->GetHistogram("MCcut","JpsiMCPt"));
    hEffTracking->SetTitle("Efficiencies");
    hist->UserHistogram("DataCuts",hEffESDCuts);
    
    TH1F *hEffTRDCuts=(TH1F*)hist->GetHistogram("DataTRDCuts","JpsiMCPt")->Clone("Efficiency");
    hEffTRDCuts->Divide(hist->GetHistogram("MCcut","JpsiMCPt"));
    hEffTRDCuts->SetTitle("Efficiencies");
    hist->UserHistogram("DataTRDCuts",hEffTRDCuts);
    
    hist->DrawSame("Efficiency");

    //2D efficiencies
    TH2F *hEffTracking2D=(TH2F*)hist->GetHistogram("DataSameMother","dndyPtMC")->Clone("2DEfficiency");
    hEffTracking2D->Divide(hist->GetHistogram("MCcut","dndyPt"));
    hEffTracking2D->SetTitle("2D Efficiency - tracking");
    hist->UserHistogram("DataSameMother",hEffTracking2D);

    TH2F *hEffESDCuts2D=(TH2F*)hist->GetHistogram("DataCuts","dndyPtMC")->Clone("2DEfficiency");
    hEffESDCuts2D->Divide(hist->GetHistogram("MCcut","dndyPt"));
    hEffESDCuts2D->SetTitle("2D Efficiency - quality cuts");
    hist->UserHistogram("DataCuts",hEffESDCuts2D);
    
    TH2F *hEffTRDCuts2D=(TH2F*)hist->GetHistogram("DataTRDCuts","dndyPtMC")->Clone("2DEfficiency");
    hEffTRDCuts2D->Divide(hist->GetHistogram("MCcut","dndyPt"));
    hEffTRDCuts2D->SetTitle("2D Efficiency - quality+TRD cuts");
    hist->UserHistogram("DataTRDCuts",hEffTRDCuts2D);
    
//
    // Draw all histograms of all histogram classes
    // Use the Draw functionality of AliDielectronHistos
    //
    hist->Draw();
    
    //
    // Draw all histograms with the same name of all classes into one canvas
    // Use the Draw functionality of AliDielectronHistos
    //
    hist->DrawSame("JpsiMCPt");
  }

  PostData(0, const_cast<THashList*>(hist->GetHistogramList()));
}

