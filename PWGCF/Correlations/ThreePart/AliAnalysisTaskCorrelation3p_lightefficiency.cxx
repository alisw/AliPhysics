/*************************************************************************
* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliAnalysisTaskCorrelation3p_lightefficiency.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliEventPoolManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliOADBContainer.h"
#include "AliCentrality.h"
#include "AliVVZERO.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TChain.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "THashList.h"
#include "TMath.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "TCanvas.h"
#include "TObjectTable.h"
ClassImp(AliAnalysisTaskCorrelation3p_lightefficiency)
//
//Task to create three particle correlations.
//Authors:
// Matthias Richter
// Paul Baetzing || pbatzing@cern.ch
//


AliAnalysisTaskCorrelation3p_lightefficiency::AliAnalysisTaskCorrelation3p_lightefficiency()
  : AliAnalysisTaskSE("AliAnalysisTaskCorrelation3p")
  , fOutput(NULL)
  , fTextBox(NULL)
  , fOption("")
  , fCentrality(NULL)
  , fVertexobj(NULL)
  , fVertex()
  , fCollisionType(AliAnalysisTaskCorrelation3p_lightefficiency::PbPb)
  , fisESD(kFALSE)
  , fisAOD(kFALSE)
  , fRemoveSignals(kFALSE)
  , fMcArray(NULL)
  , fMBinEdges(TArrayD())  
  , fZBinEdges(TArrayD())  
  , fMaxNEventMix(100)
  , fMinNofTracksMix(10)
  , fCentralityEstimator("V0M")
  , fCentralityPercentile(0)
  , fMultiplicity(0)
  , fMaxVz(10.0)
  , fMaxMult(0)
  , fMaxNumberOfTracksInPPConsidered(200)
  , fNTriggers(0)
  , fNAssociated(0)  
  , fAcceptancecut(0.8)
  , fMinPt(3.0)
  , fMaxPt(8.0)
  , fMinNClustersTPC(70)
{
  // default constructor
  // 
  Double_t Medges[6] = {0,20,40,60,80,90};
  TArrayD MBEdges(6, Medges);
  Double_t Zedges[6] = {-10,-5,-2.5,2.5,5,10};
  TArrayD ZBedges(6, Zedges);
  fMBinEdges = MBEdges;
  fZBinEdges = ZBedges;
  DefineSlots();
}

AliAnalysisTaskCorrelation3p_lightefficiency::AliAnalysisTaskCorrelation3p_lightefficiency(const char *name, const char* opt)
  : AliAnalysisTaskSE(name)
  , fOutput(NULL)
  , fTextBox(NULL)
  , fOption(opt)
  , fCentrality(NULL)
  , fVertexobj(NULL)
  , fVertex()
  , fCollisionType(AliAnalysisTaskCorrelation3p_lightefficiency::PbPb)
  , fisESD(kFALSE)
  , fisAOD(kFALSE)
  , fRemoveSignals(kFALSE)
  , fMcArray(NULL)
  , fMBinEdges(TArrayD())  
  , fZBinEdges(TArrayD())  
  , fMaxNEventMix(100)
  , fMinNofTracksMix(10)
  , fCentralityEstimator("V0M")
  , fCentralityPercentile(0)
  , fMultiplicity(0)
  , fMaxVz(10.0)
  , fMaxMult(0)
  , fMaxNumberOfTracksInPPConsidered(200)
  , fNTriggers(0)
  , fNAssociated(0)  
  , fAcceptancecut(0.8)
  , fMinPt(3.0)
  , fMaxPt(8.0)
  , fMinNClustersTPC(70)
  {
  // constructor with options
  //
  //
  Double_t Medges[6] = {0,20,40,60,80,90};
  TArrayD MBEdges(6, Medges);
  Double_t Zedges[6] = {-10,-7.5,-2.5,2.5,7.5,10};
  TArrayD ZBedges(6, Zedges);
  fMBinEdges = MBEdges;
  fZBinEdges = ZBedges;
  DefineSlots();
  }

int AliAnalysisTaskCorrelation3p_lightefficiency::DefineSlots()
{
  // define the data slots
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  return 0;
}

AliAnalysisTaskCorrelation3p_lightefficiency::~AliAnalysisTaskCorrelation3p_lightefficiency()
{
  // destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    fOutput = 0;
  }
}

void AliAnalysisTaskCorrelation3p_lightefficiency::UserCreateOutputObjects()
{
  // create result objects and add to output list
  TH1::SetDefaultSumw2(kTRUE);//want the collection of weights on all histograms.
  fOutput = new THashList;
  fOutput->SetOwner();
    if(fCollisionType==pp)for(int i=0;i<fMBinEdges.GetSize();i++){
      fMBinEdges.AddAt(fMaxNumberOfTracksInPPConsidered*fMBinEdges.At(i)/fMBinEdges.At(fMBinEdges.GetSize()-1),i);
    }
  
  InitializeQAhistograms();
  InitializeEffHistograms();
  // all tasks must post data once for all outputs
  PostData(1, fOutput);
}

void AliAnalysisTaskCorrelation3p_lightefficiency::UserExec(Option_t* /*option*/)
{
  FillHistogram("Check",0.0);
  // process the event
  TObject* pInput=InputEvent();
  if (!pInput) {AliError("failed to get input");return;}
  AliVEvent *pEvent = dynamic_cast<AliVEvent*>(pInput);
  if(!pEvent){AliError(Form("input of wrong class type %s, expecting AliVEvent", pInput->ClassName()));return;}

 //Get the Array.
  GetMCArray();

  //if it is not found, return without doing anything:
  if(!fMcArray) return;

  //Find out if it is AOD or ESD.  
  fisESD=pEvent->IsA()==AliESDEvent::Class();
  fisAOD=pEvent->IsA()==AliAODEvent::Class();
  //Get the runnumber and find which bin and fill value this corresponds to.
  GetCentralityAndVertex();
  if(!SelectEvent()) return;//events are rejected.

  if(fCollisionType==AliAnalysisTaskCorrelation3p_lightefficiency::pp) FillHistogram("centVsZVertex",fMultiplicity,fVertex[2]);//only fill with selected events.
  if(fCollisionType==AliAnalysisTaskCorrelation3p_lightefficiency::PbPb) FillHistogram("centVsZVertex",fCentralityPercentile,fVertex[2]);
  //To fill with tracks and pions:
  TObjArray allrelevantParticles;
  //Fill all the tracks
  FillHistogram("centVsNofTracks",fCentralityPercentile, GetTracks(&allrelevantParticles, pEvent));
  
  if(fNTriggers>=1)FillHistogram("NAssociatedETriggered",fNAssociated);
  //If VZERO data, fill the Multiplicity histograms.
  PostData(1, fOutput);  
}

void AliAnalysisTaskCorrelation3p_lightefficiency::FinishTaskOutput()
{
  // end of the processing
    TH1 * hist = dynamic_cast<TH1*>(fOutput->FindObject("trackCount")) ;
    if (hist) cout << "FinishTaskOutput: " << hist->GetEntries() << " events(s)" << endl;
}

void AliAnalysisTaskCorrelation3p_lightefficiency::Terminate(Option_t *)
{
  // last action on the client
//   gObjectTable->Print();

}

Int_t AliAnalysisTaskCorrelation3p_lightefficiency::GetTracks(TObjArray* allrelevantParticles, AliVEvent *pEvent)
{
  Int_t nofTracks = 0;
  nofTracks=pEvent->GetNumberOfTracks();
  FillHistogram("trackCount",nofTracks);
  for (int i=0; i<nofTracks; i++) {
    AliVParticle* t=pEvent->GetTrack(i);
    if (!t) continue;
    FillHistogram("trackUnselectedPt",t->Pt());
    FillHistogram("trackUnselectedPhi",t->Phi());
    FillHistogram("trackUnselectedTheta",t->Theta());
    if (!IsSelected(t)) continue;
    if(fRemoveSignals){

      if(IsAddedSignal(t)){ 	FillHistogram("ptrejected",t->Pt()); continue;}
      else{			FillHistogram("ptaccepted",t->Pt());}
    }

    allrelevantParticles->Add(t);
    if(fCollisionType==AliAnalysisTaskCorrelation3p_lightefficiency::pp){
      FillHistogram("hnTracksinBins",fMultiplicity,fVertex[2],t->Phi(),t->Eta(),t->Pt());
      if(fMcArray&&t->GetLabel()>=0){if(dynamic_cast<AliAODMCParticle*>(fMcArray->At(t->GetLabel()))->IsPhysicalPrimary())FillHistogram("hnTracksinBinsRecPP",fMultiplicity,fVertex[2],t->Phi(),t->Eta(),t->Pt());}
    }
    if(fCollisionType==AliAnalysisTaskCorrelation3p_lightefficiency::PbPb){
      FillHistogram("hnTracksinBins",fCentralityPercentile,fVertex[2],t->Phi(),t->Eta(),t->Pt());
      if(fMcArray&&t->GetLabel()>=0){if(dynamic_cast<AliAODMCParticle*>(fMcArray->At(t->GetLabel()))->IsPhysicalPrimary())FillHistogram("hnTracksinBinsRecPP",fCentralityPercentile,fVertex[2],t->Phi(),t->Eta(),t->Pt());}
    }
    
    FillHistogram("trackPt",t->Pt());
    FillHistogram("trackPhi",t->Phi());
    FillHistogram("trackTheta",t->Theta());
  }
  

  if(fMcArray){
    int nofMCParticles = fMcArray->GetEntriesFast();
    for (int i=0;i<nofMCParticles;i++){
      AliVParticle* t =  (AliVParticle *) fMcArray->At(i);
      if (!t) continue;
//       if(dynamic_cast<AliAODMCParticle*>(t)->GetGeneratorIndex()==1&&dynamic_cast<AliAODMCParticle*>(t)->GetMother()==-1)cout <<" In MCarray loop"<< i <<endl;
      GeneratorStat(t);
      if( t->Charge()!=0){//check if they are physical primary particles
	if (!IsSelected(t)) continue;
// 	if (IsAddedSignal(t)) continue;
	if(fCollisionType==AliAnalysisTaskCorrelation3p_lightefficiency::pp){FillHistogram("hnTracksinBinsMC",fMultiplicity,fVertex[2],t->Phi(),t->Eta(),t->Pt());}
	if(fCollisionType==AliAnalysisTaskCorrelation3p_lightefficiency::PbPb){FillHistogram("hnTracksinBinsMC",fCentralityPercentile,fVertex[2],t->Phi(),t->Eta(),t->Pt());}
	
      }
    }
  }

  return nofTracks;
}

Bool_t AliAnalysisTaskCorrelation3p_lightefficiency::IsAddedSignal(AliVParticle* p){
  if(fMcArray){
    AliAODMCParticle* p_mc = dynamic_cast<AliAODMCParticle*>(p);

    if(!p_mc){
//       if(p->GetLabel()<=-1) return kFALSE;//no mc particle corresponding to the reconstruced track.
      p_mc = dynamic_cast<AliAODMCParticle*>(fMcArray->At(abs(p->GetLabel())));
    }
//     p_mc->Print();
    if(p_mc->GetMother()==-1){//in the case where p_mc has no mother, check the origin
      if(p_mc->GetGeneratorIndex()==1) return kTRUE;
    }
    while(p_mc->GetMother()>-1){//if the particle has a mother, follow the mother to the origin, and check the generator
//       p_mc->Print();
//     if(p_mc->GetMother()<0){p_mc->Print(); cout << p_mc->GetMother()<<endl;}

     AliAODMCParticle* p_mc_mother = dynamic_cast<AliAODMCParticle*>(fMcArray->At(p_mc->GetMother()));
//      p_mc_mother->Print();
     if(p_mc_mother->GetMother()==-1){
//        p_mc_mother->Print();
//        cout << p_mc_mother->GetGeneratorIndex()<<endl;
//        if(p_mc_mother->GetGeneratorIndex()==0){
// 	 //is from hijing
// 	 cout << "Hijing:"<<endl;
// 	 p_mc_mother->Print();
// 	 return kFALSE;
//        }
       if(p_mc_mother->GetGeneratorIndex()!=0){
	 //is from added signal
// 	 cout << "Added signal:"<<endl;
// 	 p_mc_mother->Print();
	 return kTRUE;
       }
       else return kFALSE;//
     }
     else p_mc = p_mc_mother;
    }
    return kFALSE;
  } 
  return kTRUE;
}


Bool_t AliAnalysisTaskCorrelation3p_lightefficiency::IsSelected(AliVParticle* p)
{
  //Performs selection cuts for tracks and triggers
  if (p->IsA()==AliESDtrack::Class() && IsSelectedTrackESD(p)) return IsSelectedTrack(p);
  if (p->IsA()==AliAODTrack::Class() && IsSelectedTrackAOD(p)) return IsSelectedTrack(p);
  if (p->IsA()==AliAODMCParticle::Class() && dynamic_cast<AliAODMCParticle*>(p)->IsPhysicalPrimary()) return IsSelectedTrack(p);
  return kFALSE;
}

Bool_t AliAnalysisTaskCorrelation3p_lightefficiency::IsSelectedTrack(AliVParticle* p)
{
  if (p->Pt()<=fMinPt) return kFALSE;
  if (fMaxPt>fMinPt && p->Pt()>fMaxPt) return kFALSE;
  float etatrigger=p->Eta();
  if (etatrigger<=-fAcceptancecut || etatrigger>=fAcceptancecut) return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskCorrelation3p_lightefficiency::IsSelectedTrackAOD(AliVParticle* t)
{
  AliAODTrack *AODt = dynamic_cast<AliAODTrack*>(t);
  Bool_t isselected = kTRUE;
  Double_t DCAtang=-999.0;
  Double_t DCAlong=-999.0;
  GetDCA(DCAtang,DCAlong,AODt);
  if((AODt->HasPointOnITSLayer(1)||AODt->HasPointOnITSLayer(2)))   FillHistogram("TrackDCAandonITS",DCAtang,DCAlong,1);
  else FillHistogram("TrackDCAandonITS",DCAtang,DCAlong,0);
//   isselected = isselected&&(AODt->TestFilterBit(BIT(4)));  //filter bits: BIT(4) = standard cuts, loose DCA; BIT(5) standard cuts, tight DCA 
// //  if(isselected) cout << "FilterBitPassed"<<endl;
//   isselected = isselected&&(AODt->GetFilterMap()&AliVTrack::kITSrefit);
// //   if(isselected) cout << "ITSrefit passed"<<endl;
// //   isselected = isselected&&((AODt->GetFilterMap()&AliAODTrack::kTPCrefit)||fCollisionType==PbPb);//in the PbPb AODs it seems this is not set.
// //   if(isselected) cout << "TPCrefit passed"<<endl;
//   isselected = isselected&&(AODt->HasPointOnITSLayer(1)||AODt->HasPointOnITSLayer(2));//in first or second ITS layer
// //   if(isselected) cout << "ITS any passed"<<endl;
// //   isselected = isselected&&(AODt->HasPointOnITSLayer(2));//in second ITS layer
// //   if(isselected) cout << "ITS layer 2 passed"<<endl;
//   isselected = isselected&&(AODt->GetTPCNcls()>70);
// //   if(isselected) cout << "more than 70 TPC clusters"<<endl;
//   isselected = isselected&&(abs(DCAtang)<0.5);//cm. DCA less then 0.5 cm in transverse direction.
// //   if(isselected) cout << "DCA tang passed"<<endl;
//   isselected = isselected&&(abs(DCAlong)<3);//cm. DCA less then 3 cm in the longitudinal direction.
// //   if(isselected) cout << "DCA long passed"<<endl;
  //Hybrid tracks give flat distributions
  isselected = AODt->IsHybridGlobalConstrainedGlobal();
  if( (AODt->HasPointOnITSLayer(1)||AODt->HasPointOnITSLayer(2))&&isselected)   FillHistogram("TrackDCAandonITSselected",DCAtang,DCAlong,1);
  if(!(AODt->HasPointOnITSLayer(1)||AODt->HasPointOnITSLayer(2))&&isselected)   FillHistogram("TrackDCAandonITSselected",DCAtang,DCAlong,0);
  return isselected; 
}

Bool_t AliAnalysisTaskCorrelation3p_lightefficiency::IsSelectedTrackESD(AliVParticle* t)
{
  if(t)return kFALSE;//ESD is currently not supported
  else return kFALSE;
}

void AliAnalysisTaskCorrelation3p_lightefficiency::GeneratorStat(AliVParticle* p)
{
  //Code to fill the generator stats histogram for all particles that have no mother
  AliAODMCParticle* mcpart = dynamic_cast<AliAODMCParticle*>(p);
  if(mcpart){
    if(mcpart->GetMother()==-1){
      FillHistogramGenPar("genparticle",FillPDG(mcpart->GetPdgCode()),FillGen(mcpart->GetGeneratorIndex()),mcpart->Pt());
    }
  }
}

Int_t AliAnalysisTaskCorrelation3p_lightefficiency::FillGen(Int_t GetGeneratorIndex)
{
//   cout << GetGeneratorIndex <<endl;
  switch(GetGeneratorIndex){
    case 0 : return 1;//"Generator 1 (Hijing?)"
    case 1 : return 2;//"Generator 2 (Pythia/embedded?)"
  }
  return 3;
}

const char* AliAnalysisTaskCorrelation3p_lightefficiency::NameGen(Int_t BinGenIndex)
{
  switch(BinGenIndex){
    case 1 : return "Generator 1 (Hijing?)";
    case 2 : return "Generator 2 (Pythia/embedded?)";
    case 3 : return "other generators (GEANT and co?)";
  }
  return "";
}



Int_t AliAnalysisTaskCorrelation3p_lightefficiency::FillPDG(Int_t GetPDG)
{
  Int_t abspdg = TMath::Abs(GetPDG);
  if(abspdg<=9){//quarks
    if(abspdg<=3) return 1;//up, down or strange quark
    if(abspdg==4) return 2;//charm quark
    if(abspdg==5) return 3;//bottom quark
    if(abspdg==6) return 4;//top quark
    if(abspdg>6&&abspdg<9) return 99;//other particles
    if(abspdg==9) return 5;//gluon
  }
  if(abspdg>9&&abspdg<11) return 99;//other particles
  if(abspdg>=11&&abspdg<19){//leptons
    if(abspdg==11||abspdg==13||abspdg==15) return 6;//electrons, muons and tauons
    if(abspdg>15||abspdg==12||abspdg==14) return 99;//other particles
  }
  if(abspdg>=19&&abspdg<21) return 99;//other particles
  if(abspdg>=21&&abspdg<38){//Gauge particles
    if(abspdg == 21) return 5;//gluon
    if(abspdg == 22||abspdg==23||abspdg==25)return 7;//neutral gauge boson
    if(abspdg==24)return 8;//W boson
    if(abspdg>25) return 99;//any other boson
  }
  if(abspdg>=38&&abspdg<111) return 99; //any other particle
  if(abspdg>=111&&abspdg<411){//light mesons, I = 1/0, strange mesons
    if(abspdg==111) return 9;//pi0
    if(abspdg==211) return 10;//pi+-
    if(abspdg==113||abspdg==213) return 11;//rho mesons
    if(abspdg==115||abspdg==215||abspdg==117||abspdg==217||abspdg==119||abspdg==219)return 12;//a/b/rho3 mesons
    if(abspdg==221||abspdg==331)return 13;//eta meson
    if(abspdg==223||abspdg==333||abspdg==225||abspdg==335||abspdg==227||abspdg==337||abspdg==229)return 14;//eta meson
    if(abspdg==130||abspdg==310||abspdg==311||abspdg==313||abspdg==315||abspdg==317||abspdg==319)return 15;//neutral kaon
    if(abspdg==321||abspdg==323||abspdg==325||abspdg==327||abspdg==329)return 16;//charged kaon
    return 99;//if none other is chosen, return other particle
  }
 if(abspdg>=411&&abspdg<511){//charmed mesons
   if(abspdg==411||abspdg==413||abspdg==415||abspdg==431||abspdg==433||abspdg==435)return 18;//charged D meson
   if(abspdg==421||abspdg==423||abspdg==425)return 17;//neutral D meson
   if(abspdg==443)return 21;//J/Psi
   if(abspdg==441||abspdg==445)return 22;//other ccbar mesons.
   return 99;//if none other is chosen, return other particle
 }
 if(abspdg>=511&&abspdg<1114){//bottom mesons, other particles up to proton.
   if(abspdg==511||abspdg==513||abspdg==515||abspdg==531||abspdg==533||abspdg==535)return 19;//neutral B meson
   if(abspdg==521||abspdg==523||abspdg==525||abspdg==541||abspdg==543||abspdg==545)return 20;//charged B meson
   if(abspdg==553)return 23;//Ypsilon
   if(abspdg==551||abspdg==555||abspdg==557)return 22;//other bbbar mesons.
   return 99;//if none other is chosen, return other particle
 }
 if(abspdg>=1114&&abspdg<3122){//light baryons
   if(abspdg==1114||abspdg==2114||abspdg==2214||abspdg==2224)return 27;//Delta baryons
   if(abspdg==2212)return 25;//proton
   if(abspdg==2112)return 26;//neutron
   return 99;//if none other is chosen, return other particle
 }
 if(abspdg>=3122&&abspdg<3335){
   return 28;//strange baryon
 }
 if(abspdg>=3335&&abspdg<4112)return 99;
 if(abspdg>=4122&&abspdg<4444)return 29;//charmed baryon
 if(abspdg>=4444&&abspdg<5122)return 99;
 if(abspdg>=5122&&abspdg<5554)return 30;//bottom baryon
 
}

const char* AliAnalysisTaskCorrelation3p_lightefficiency::NamePDG(Int_t BinPDG)
{
  //switch the bin number to the name:
  switch(BinPDG){
    case 1 : return "u/d/s";//up , down or strange quark
    case 2 : return "c";//charm quark
    case 3 : return "b";//bottom quark
    case 4 : return "t";//top quark
    case 5 : return "g";//gluon
    case 6 : return "e/mu/tau";//electrons, muons and tauons
    case 7 : return "G0";//neutral gauge boson
    case 8 : return "W";//W boson
    case 9 : return "pi0";//pi0
    case 10 : return "pi+-";//pi+-
    case 11 : return "rho";//rho mesons
    case 12 : return "a/b/rho3";//b/a/rho3 mesons
    case 13 : return "eta";//eta meson
    case 14 : return "omega/phi/f";//omega, phi or f mesons
    case 15 : return "K0";//neutral kaon
    case 16 : return "K+-";//charged kaon
    case 17 : return "D0";//neutral D meson
    case 18 : return "D+-";//charged D meson
    case 19 : return "B0";//neutral B meson
    case 20 : return "B+-";//charged B meson
    case 21 : return "J/Psi";//J/Psi
    case 22 : return "other ccbar";//Other ccbar mesons
    case 23 : return "Ypsilon";//Ypsilon meson
    case 24 : return "other bbbar";//other bbbar mesons.
    case 25 : return "p";//proton
    case 26 : return "n";//neutron
    case 27 : return "Delta";//delta baryon
    case 28 : return "Strange baryons";//strange baryons
    case 29 : return "Charmed baryons";//charmed baryons
    case 30 : return "Bottom baryons";//bottom baryons
    case 99 : return "other particles";//other
  }
  return "";
}

void AliAnalysisTaskCorrelation3p_lightefficiency::GetMCArray()
{//Gets the MCarray if we are in AOD.
    fMcArray = 0;
    AliAODInputHandler* aodHandler=dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (aodHandler){
      AliAODEvent *aod=aodHandler->GetEvent();
      if (aod ) {
	fMcArray = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
	if (!fMcArray) AliError("Could not retrieve MC array!");
      }
      else AliError("Could not retrieve AOD event! MC is only supported in AOD.");
    }
}

void AliAnalysisTaskCorrelation3p_lightefficiency::GetDCA(Double_t& DCAtang, Double_t& DCAlong, AliAODTrack* AODt)
{
if(AODt->TestBit(AliAODTrack::kIsDCA)){
  DCAtang = AODt->DCA();
  DCAlong = AODt->ZAtDCA();
}
else{
  if(fVertex){
    Double_t fBzkg = dynamic_cast<AliAODEvent*>(InputEvent())->GetMagneticField();
    Double_t* dca = new Double_t[2];
    Double_t* dcacov = new Double_t[3];

    Double_t kVeryBigno = 1000000;
    if(AODt->PropagateToDCA(fVertexobj,fBzkg,kVeryBigno,dca,dcacov)){DCAtang=dca[0];DCAlong = dca[1];}
    else{DCAtang = -999;DCAlong=-999;}
    delete dca;
    delete dcacov;
    }
  }
} 


void AliAnalysisTaskCorrelation3p_lightefficiency::GetCentralityAndVertex()
{
  fVertexobj = InputEvent()->GetPrimaryVertex();
  if(fisAOD)
  {
    // Fill AliAODEvent interface specific information
    AliAODHeader *header = dynamic_cast<AliAODHeader*>(dynamic_cast<AliAODEvent*>(InputEvent())->GetHeader());
    fCentrality =  dynamic_cast<AliCentrality*>(header->GetCentralityP());
    fMultiplicity = header->GetRefMultiplicity();
  }  
  if(fisESD)
  {
    fCentrality = dynamic_cast<AliESDEvent*>(InputEvent())->GetCentrality();
    fMultiplicity = 0.0 ;//not implemented atm
  }
  if(fCentrality)fCentralityPercentile = fCentrality->GetCentralityPercentile(fCentralityEstimator);

  if(fCentrality)FillHistogram("centrality",fCentralityPercentile);
  FillHistogram("multiplicity",fMultiplicity);
  //Get the primary Vertex
  if( fVertexobj ) {
    fVertex[0] = fVertexobj->GetX();
    fVertex[1] = fVertexobj->GetY();
    fVertex[2] = fVertexobj->GetZ();}
  else return;
}
Bool_t AliAnalysisTaskCorrelation3p_lightefficiency::SelectEvent()
{//This function provides the event cuts for this class.
  if(fCollisionType==pp){//With pp, the following cuts are applied:
    FillHistogram("Eventbeforeselection",fVertex[2],fMultiplicity,0);
    if(!fVertexobj){AliError("Vertex object not found.");return kFALSE;}//Runs only after GetCentralityAndVertex().
    if(fVertexobj->GetNContributors()<1) return kFALSE; // no tracks go into reconstructed vertex
    if(abs(fVertex[2])>fMaxVz) return kFALSE;//Vertex is too far out
    if(fMultiplicity>fMaxNumberOfTracksInPPConsidered) return kFALSE;//Out of multiplicity bounds in pp, no histograms will be filled.
    if(InputEvent()->IsPileupFromSPD(3,0.8,3.,2.,5.))return kFALSE;  //reject for pileup.
    FillHistogram("Eventafterselection",fVertex[2],fMultiplicity,0);
    
  }
  if(fCollisionType==PbPb){
    FillHistogram("Eventbeforeselection",fVertex[2],fMultiplicity,fCentralityPercentile);
    if(!fVertexobj){AliError("Vertex object not found.");return kFALSE;}//Runs only after GetCentralityAndVertex().
    if(fVertexobj->GetNContributors()<1) return kFALSE; // no tracks go into reconstructed vertex
    if(abs(fVertex[2])>fMaxVz) return kFALSE;//Vertex is too far out
    if(!fCentrality){AliError("Centrality object not found.");return kFALSE;}//Centrality must be defined in the PbPb case.
    if(fCentrality->GetQuality()!=0)return kFALSE;//bad centrality.
    if(fCentralityPercentile<0) return kFALSE;//centrality is not defined
    if(fCentralityPercentile>fMaxMult) return kFALSE;//Out of centrality bounds in PbPb, will not fill any histogram.
    FillHistogram("Eventafterselection",fVertex[2],fMultiplicity,fCentralityPercentile);

  }
  return kTRUE;
}



void AliAnalysisTaskCorrelation3p_lightefficiency::InitializeQAhistograms()
{
  //Function that initializes the QA histograms 
  if (!fOutput) return;
  //QA histograms
  fOutput->Add(new TH1D("Check","Check",1,-0.5,0.5));

  fOutput->Add(new TH3D("TrackDCAandonITS","DCA tangential vs DCA longitudinal vs Is in the first two ITS layers",50,-2,2,50,-5,5,2,-0.5,1.5));
  fOutput->Add(new TH3D("TrackDCAandonITSselected","DCA tangential vs DCA longitudinal vs Is in the first two ITS layers for selected events",50,-2,2,50,-5,5,2,-0.5,1.5));
  fOutput->Add(new TH3D("Eventbeforeselection","Vertex vs Multiplicity vs Centrality before event selection.", 50,-15,15,50,0,4000,50,0,100));
  fOutput->Add(new TH3D("Eventafterselection","Vertex vs Multiplicity vs Centrality after event selection.", 50,-15,15,50,0,4000,50,0,100));
  fOutput->Add(new TH1D("trackCount", "trackCount", 1000,  0, 4000));
  fOutput->Add(new TH1D("trackUnselectedPt"   , "trackPt"   , 1000,  0, 20));
  fOutput->Add(new TH1D("trackPt"   , "trackPt"   , 100,  0, 20));
  fOutput->Add(new TH1D("trackUnselectedPhi"  , "trackPhi"  ,  180,  0., 2*TMath::Pi()));
  fOutput->Add(new TH1D("trackPhi"  , "trackPhi"  ,  180,  0., 2*TMath::Pi()));
  fOutput->Add(new TH1D("trackUnselectedTheta", "trackTheta",  180, 0, TMath::Pi()));
  fOutput->Add(new TH1D("trackTheta", "trackTheta",  180, 0.0, TMath::Pi()));
  fOutput->Add(new TH1D("Ntriggers","Number of triggers per event",50,-0.5,49.5));
  
  fOutput->Add(new TH1D("centrality", "Centrality",  100,  0, 100));
  fOutput->Add(new TH1D("multiplicity", "Multiplicity of tracklets",  100,  0, fMaxNumberOfTracksInPPConsidered));
  fOutput->Add(new TH2D("centVsNofTracks", "centVsNofTracks", 100, 0, 100, 100, 0, 2000));
  if(fCollisionType==PbPb)fOutput->Add(new TH2D("centVsZVertex", "centvszvertex", 100, 0, 100, 100, -10, 10));
  if(fCollisionType==pp)fOutput->Add(new TH2D("centVsZVertex", "centvszvertex", 100, 0, fMaxNumberOfTracksInPPConsidered, 100, -10, 10));
  if(fRemoveSignals){ 
    fOutput->Add(new TH1D("ptrejected","pt distribution of primary particles that are added signals.",100,  0, 20));
    fOutput->Add(new TH1D("ptaccepted","pt distribution of particles that are not added signals.",100,  0, 20));
    TH3D* genparticlehist = new TH3D("genparticle","Particle types vs generator type vs pT",99,0,1,3,0,1,100,  0, 20);
    for(int i = 1;i<=99;i++){
      TString label = NamePDG(i);
      genparticlehist->GetXaxis()->SetBinLabel(i,label);
    }
    for(int i = 1;i<=3;i++){
      TString label = NameGen(i);
      genparticlehist->GetYaxis()->SetBinLabel(i,label);
    }
    fOutput->Add(genparticlehist);
  }
}

void AliAnalysisTaskCorrelation3p_lightefficiency::InitializeEffHistograms()
{
  //Function that initializes the Efficiency histogram
  Int_t    nofMBins=2*(fMBinEdges.GetSize()-1);
  Double_t MBinmined=fMBinEdges.At(0);
  Double_t MBinmaxed=fMBinEdges.At(fMBinEdges.GetSize()-1);
  Double_t Mbins[nofMBins+1];
  for (int i = 0;i<=nofMBins;i++){
    if(i%2==0){
      Mbins[i] = fMBinEdges.At(i/2);
    }
    else{
      Mbins[i] = fMBinEdges.At((i-1)/2) + (fMBinEdges.At((i+1)/2) - fMBinEdges.At((i-1)/2))/2.0;
    }
  }
  Int_t    nofZBins=2*(fZBinEdges.GetSize()-1);
  Double_t ZBinmined=fZBinEdges.At(0);
  Double_t ZBinmaxed=fZBinEdges.At(fZBinEdges.GetSize()-1);
  Double_t Zbins[nofZBins+1];
  for (int i = 0;i<=nofZBins;i++){
    if(i%2==0){
      Zbins[i] = fZBinEdges.At(i/2);
    }
    else{
      Zbins[i] = fZBinEdges.At((i-1)/2) + (fZBinEdges.At((i+1)/2) - fZBinEdges.At((i-1)/2))/2.0;
    }
  }

  Int_t nphi = 1;//36;//72;
  Double_t phimin = 0.0;
  Double_t phimax = 2.0*TMath::Pi();
  Int_t nEta = 63;//126;
  Double_t EtaMin = -0.9;
  Double_t EtaMax =  0.9;
  Double_t pTmin = fMinPt;
  Double_t pTmax = fMaxPt;
  Int_t npT = (pTmax-pTmin)/0.1 - 0.5;//rounding up

  Int_t  bins[5]   = {nofMBins, nofZBins,nphi,nEta,npT};
  Double_t xmin[5] = {MBinmined,ZBinmined,phimin,EtaMin,pTmin};
  Double_t xmax[5] = {MBinmaxed,ZBinmaxed,phimax,EtaMax,pTmax};
  fOutput->Add(new THnF("hnTracksinBins","Tracks in different bins.",5,bins,xmin,xmax));
  dynamic_cast<THnF*>(fOutput->FindObject("hnTracksinBins"))->GetAxis(0)->Set(nofMBins,Mbins);
  dynamic_cast<THnF*>(fOutput->FindObject("hnTracksinBins"))->GetAxis(1)->Set(nofZBins,Zbins);
  fOutput->Add(new THnF("hnTracksinBinsRecPP","Tracks in different bins from reconstruction that originate from a PhysicalPrimary MC particle.",5,bins,xmin,xmax));
  dynamic_cast<THnF*>(fOutput->FindObject("hnTracksinBinsRecPP"))->GetAxis(0)->Set(nofMBins,Mbins);
  dynamic_cast<THnF*>(fOutput->FindObject("hnTracksinBinsRecPP"))->GetAxis(1)->Set(nofZBins,Zbins);
  fOutput->Add(new THnF("hnTracksinBinsMC","Tracks in different bins, MC truth for charged particles.",5,bins,xmin,xmax));
  dynamic_cast<THnF*>(fOutput->FindObject("hnTracksinBinsMC"))->GetAxis(0)->Set(nofMBins,Mbins);
  dynamic_cast<THnF*>(fOutput->FindObject("hnTracksinBinsMC"))->GetAxis(1)->Set(nofZBins,Zbins);
}



void AliAnalysisTaskCorrelation3p_lightefficiency::SetMixingScheme(Int_t MaxNEventMix, Int_t MinNofTracksMix, TArrayD MBinEdges, TArrayD ZBinEdges)
{
  fMaxNEventMix= MaxNEventMix;
  fMinNofTracksMix = MinNofTracksMix;
  cout << MBinEdges.GetSize()<<endl;
  for(int i=0; i<MBinEdges.GetSize()-1; ++i)
    if(MBinEdges.At(i) > MBinEdges.At(i+1)) AliFatal("edges are not sorted");
  for(int i=0; i<ZBinEdges.GetSize()-1; ++i)
    if(ZBinEdges.At(i) > ZBinEdges.At(i+1)) AliFatal("edges are not sorted");  
  fMBinEdges = MBinEdges;
  fZBinEdges = ZBinEdges;
  fMaxMult = fMBinEdges.At(fMBinEdges.GetSize()-1);
}


void AliAnalysisTaskCorrelation3p_lightefficiency::FillHistogram(const char* key, Double_t x)
{
  TH1 * hist = dynamic_cast<TH1*>(fOutput->FindObject(key)) ;
  if(hist)
    hist->Fill(x) ;
  else AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}

void AliAnalysisTaskCorrelation3p_lightefficiency::FillHistogram(const char* key, Double_t x, Double_t y)
{
  TH2 * hist = dynamic_cast<TH2*>(fOutput->FindObject(key)) ;
  if(hist)
    hist->Fill(x,y) ;
  else AliError(Form("can not find histogram (of instance TH2) <%s> ",key)) ;
}

void AliAnalysisTaskCorrelation3p_lightefficiency::FillHistogram(const char* key, Double_t x, Double_t y, Double_t z)
{
  TH3 * hist = dynamic_cast<TH3*>(fOutput->FindObject(key)) ;
  if(hist)
    hist->Fill(x,y,z) ;
  else AliError(Form("can not find histogram (of instance TH3) <%s> ",key)) ;
}

void AliAnalysisTaskCorrelation3p_lightefficiency::FillHistogram(const char* key, Double_t x, Double_t y, Double_t z,Double_t a, Double_t b)
{
  THnF * hist = dynamic_cast<THnF*>(fOutput->FindObject(key)) ;
  if(hist){
    Double_t s[5] = {x,y,z,a,b};
    hist->Fill(s) ;
  }
  else AliError(Form("can not find histogram (of instance THNf) <%s> ",key)) ;
}

void AliAnalysisTaskCorrelation3p_lightefficiency::FillHistogramGenPar(const char* key, Int_t x, Int_t y,Double_t z)
{
  TH3 * hist = dynamic_cast<TH3*>(fOutput->FindObject(key)) ;
  if(hist){
    hist->Fill(hist->GetXaxis()->GetBinCenter(x),hist->GetYaxis()->GetBinCenter(y),z) ;}
  else AliError(Form("can not find histogram (of instance TH3) <%s> ",key)) ;
}

