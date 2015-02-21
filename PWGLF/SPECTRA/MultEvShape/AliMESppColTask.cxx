#include "AliLog.h"

#include "TROOT.h"
#include "TString.h"
#include "TH1.h"
#include "THnSparse.h"
#include "TTreeStream.h"
#include "TGeoGlobalMagField.h"

#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "TParticlePDG.h"
#include "TParticle.h"

#include "AliMESppColTask.h"
#include "AliMESeventInfo.h"
#include "AliMEStrackInfo.h"
#include "AliMCEvent.h"

#include "AliEventPoolManager.h"


ClassImp(AliMESppColTask)
ClassImp(AliMESppColTask::AliMESppColTaskExchange)
ClassImp(AliMESppColTask::AliMESppColMixEvent)

AliMESppColTask::AliMESppColTaskExchange::AliMESppColTaskExchange()
  : TObject()
  ,fN(0)
{
  //
  // Constructor
  //
  for(Int_t i(0); i<NMAXMULT; i++){
    fDEta[i] = -999.;
    fDPhi[i] = -999.;
  }
}

//________________________________________________________________________
void AliMESppColTask::AliMESppColTaskExchange::Add(Float_t de, Float_t dp)
{
  if(fN>=NMAXMULT){
    AliWarning("Correlation bufer exhausted");
    return;
  }
  //AliInfo(Form("Add correl @ %d [%f %f]", n, de, dp));
  fDEta[fN] = de;
  fDPhi[fN] = dp;
  fN++;
}

//________________________________________________________________________
AliMESppColTask::AliMESppColMixEvent::AliMESppColMixEvent()
  : TNamed()
  ,fParent(NULL)
  ,fPoolManager(NULL)
  ,fPool(NULL)
//   ,fEvInfo(NULL)
  ,fAssociatedTracks(NULL)
//   ,fMixedTrack(NULL)
  ,fmixing(kFALSE)
  ,fNTracks(0)
  ,fPoolContent(0)  
  ,fMult(-1)
   ,fPhiMin(0)
  ,fPhiMax(0)
  ,fPtTrigg(0)
  ,fEtaTrigg(0)
  ,fPhiTrigg(0) 
  ,fidTrigg(0)
  ,fDeltaEta(0)
  ,fDeltaPhi(0)
{
  //
  // Constructor
  //
}

//________________________________________________________________________
AliMESppColTask::AliMESppColMixEvent::AliMESppColMixEvent(AliMESppColTask *p)
  : TNamed()
  ,fParent(p)
  ,fPoolManager(NULL)
  ,fPool(NULL)
//   ,fEvInfo(NULL)
  ,fAssociatedTracks(NULL)
//   ,fMixedTrack(NULL)
  ,fmixing(kFALSE)
  ,fNTracks(0)
  ,fPoolContent(0)  
  ,fMult(-1)
  ,fPhiMin(0)
  ,fPhiMax(0)
  ,fPtTrigg(0)
  ,fEtaTrigg(0)
  ,fPhiTrigg(0)  
  ,fidTrigg(0)
  ,fDeltaEta(0)
  ,fDeltaPhi(0)
{
  //
  // Constructor
  //
}

AliMESppColTask::AliMESppColMixEvent::~AliMESppColMixEvent()
{
  //
  // Destructor
  //
//   printf("fPoolManager[%x]\n", (void*)fPoolManager);
  if(fPoolManager) {delete fPoolManager; fPoolManager=0;}
//   printf("fPool[%x]\n", (void*)fPool);
  if(fPool) {delete fPool; fPool=0;}
//   if(fEvInfo) {delete fEvInfo; fEvInfo=0;}
//   printf("fAssociatedTracks[%x]\n", (void*)fAssociatedTracks);
  if(fAssociatedTracks) {
    fAssociatedTracks->Clear(); 
//     printf("fAssociatedTracks->Clear()\n");
    delete fAssociatedTracks; 
    fAssociatedTracks=0;
  }
//   printf("fMixedTrack[%x]\n", (void*)fMixedTrack);
// //   if(fMixedTrack) {delete fMixedTrack; fMixedTrack=0;}
//     if(fMixedTrack) {
//     fMixedTrack->Clear(); 
//     printf("fMixedTrack->Clear()\n");
//     delete fMixedTrack; 
//     fMixedTrack=0;
//   }
  if(fNTracks) fNTracks=0;
  
  if(fPhiMin) fPhiMin=0;
  if(fPhiMax) fPhiMax=0;
  
  if(fPtTrigg)  fPtTrigg=0;
  if(fPhiTrigg) fPhiTrigg=0;
  if(fEtaTrigg) fEtaTrigg=0;
  if(fidTrigg) fidTrigg=0;
  
  
  if(fDeltaPhi) fDeltaPhi=0;
  if(fDeltaEta) fDeltaEta=0;
  
}
//________________________________________________________________________
AliMESppColTask::AliMESppColTask()
  : AliMESbaseTask()
  , fCorrelator(0x0)
//   , fmixing(kFALSE)
{
  //
  // Constructor
  //
}

//________________________________________________________________________
AliMESppColTask::AliMESppColTask(const char *name)
  : AliMESbaseTask(name)
  , fCorrelator(0x0)
//   , fmixing(kTRUE)
{
  //
  // Constructor
  //
  DefineOutput(AliMESbaseTask::kQA,          TList::Class());
}

//________________________________________________________________________
AliMESppColTask::~AliMESppColTask()
{
  //deconstructor
  if(fCorrelator) {delete fCorrelator; fCorrelator = 0;}
}

//________________________________________________________________________
void AliMESppColTask::UserCreateOutputObjects()
{
  //define user data containers
  AliMESbaseTask::UserCreateOutputObjects();  
//   fmixing=kTRUE;
  fCorrelator = new AliMESppColMixEvent(this);
  fCorrelator->SetEventMixing(kTRUE);
  
  Double_t Pi = TMath::Pi();
  fCorrelator->SetDeltaPhiInterval( -0.5*Pi, 1.5*Pi);
  
  Bool_t pooldef = fCorrelator->DefineEventPool();
  if(!pooldef) AliInfo("Warning:: Event pool not defined properly");
  
  
}

//________________________________________________________________________
void AliMESppColTask::UserExec(Option_t *opt)
{
  // Run user analysis. The following objects are allocated after calling AliMESbaseTask::UserExec(opt)
  // fEvInfo  -  reconstructed event information (class AliMESeventInfo)
  // fTracks  -  reconstructed array of tracks (class TObjArray of AliMEStrackInfo)
  // fMCevInfo-  MC event information (class AliMESeventInfo)
  // fMCtracks-  MC array of tracks (class TObjArray of AliMEStrackInfo)
  AliMESbaseTask::UserExec(opt);
  Bool_t UseEM=kTRUE;
  fCorrelator->SetEventMixing(kTRUE);
  
  if(!(fEvInfo   = dynamic_cast<AliMESeventInfo*>(GetInputData(AliMESbaseTask::kEventInfo)))){ 
          AliInfo("Ev Info Not defined.");
          return;
    }
  if(!fEvInfo->HasTriggerMB() && !fEvInfo->HasTriggerHM()) return;
  if(!fEvInfo->HasVertex()) return;  
  if(TMath::Abs(fEvInfo->GetVertexZ())>10.) return; 
  if((fEvInfo->GetMultiplicity(AliMESeventInfo::kComb))<0) return;

  if(!(fTracks   = dynamic_cast<TObjArray*>(GetInputData(AliMESbaseTask::kTracks)))){
    AliError("REC track array missing. Processing skipped");
    return;
  }
  AliMEStrackInfo *t(NULL); 
  Double_t lead[6] = {0., 0., 0., 0., 0., 0.};
  Double_t lem[6] = {0., 0., 0., 0., 0., 0.};
  THnSparse *L(NULL), *Lem(NULL);
  L=(THnSparse*)fHistosQA->At(0);
  Lem=(THnSparse*)fHistosQA->At(1);
  
  Bool_t correlatorON = fCorrelator->InitializeEventPool();
  if(!correlatorON) {
    AliInfo("AliMESppColMixEvent didn't initialize the pool correctly or processed a bad event");
    return;
  }

  Double_t etaTrigg(-999.), phiTrigg(-999.); Int_t idTrigg(-999);   Double_t ptTrigg(0.1);
 
  for(Int_t iTrack = 0; iTrack<fTracks->GetEntries(); iTrack++){
    if(!(t = (AliMEStrackInfo*)fTracks ->At(iTrack))) continue;
    if((t->Pt())>ptTrigg){
      ptTrigg=t->Pt();
      etaTrigg=t->Eta();
      phiTrigg=t->Phi();
      idTrigg = iTrack;
    }
  }
  
  lead[0] = fTracks->GetEntries(); 
  lead[1] = ptTrigg; 
  lead[2] = etaTrigg;
  lead[3] = phiTrigg;
  
   for (Int_t iTracks = 0; iTracks < fTracks->GetEntries(); iTracks++) {
    if(!(t = (AliMEStrackInfo*)fTracks ->At(iTracks))) continue;
    Double_t dEta=-999., dPhi=-999.;
  //------celelalte trakuri inafara de leading  
    if(iTracks!=idTrigg){
      dEta = t->Eta()- etaTrigg; lead[4] = t->Eta()- etaTrigg;
      dPhi = t->Phi()- phiTrigg; lead[5] = t->Phi()- phiTrigg;     
      if(dPhi<-TMath::PiOver2()) { dPhi+=TMath::TwoPi(); lead[5]+=TMath::TwoPi();}
      else if(dPhi>TMath::Pi()+TMath::PiOver2()) { dPhi-=TMath::TwoPi(); lead[5]-=TMath::TwoPi();}
      
      if(L) L->Fill(lead);
      if(DebugLevel()>0){

      (*AliMESbaseTask::DebugStream()) << "leadingSE"
        <<"ptTrigg=" << ptTrigg
        <<"etaTrigg=" << etaTrigg
        <<"phiTrigg=" << phiTrigg
        <<"dEta=" << dEta
        <<"dPhi=" << dPhi        
        << "\n";
      }
    }
   }

//   phiTrigg = fCorrelator->SetCorrectPhiRange(phiTrigg);
//   fCorrelator->SetTriggerParticleProperties(ptTrigg, phiTrigg, etaTrigg, idTrigg);
  
//   cout<<"111111111111111111111111111111"<< endl;
  Bool_t execPool = fCorrelator->ProcessEventPool();
  
  if(UseEM && !execPool) AliWarning("Mixed event analysis: pool is not ready"); 
  else{
    Int_t NofEventsinPool = 1;

    if(UseEM) NofEventsinPool = fCorrelator->GetNofEventsInPool();
    
    for (Int_t jMix =0; jMix < NofEventsinPool; jMix++){// loop on events in the pool; if it is SE analysis, stops at one -> to be done
      Bool_t analyzetracks = fCorrelator->ProcessAssociatedTracks(jMix);
      if(!analyzetracks) {
        AliInfo("AliMESppColMixEvent::Cannot process the track array");
        continue;
      }
      Int_t NofTracks = fCorrelator->GetNofTracks();
      
      Double_t etaTriggPool(-999.), phiTriggPool(-999.); Int_t idTriggPool(-999);   Double_t ptTriggPool(0.1);
 
      for(Int_t iTrackPool = 0; iTrackPool<NofTracks; iTrackPool++){
        if(!(t = (AliMEStrackInfo*)fTracks ->At(iTrackPool))) continue;
        if((t->Pt())>ptTriggPool){
          ptTriggPool=t->Pt();
          etaTriggPool=t->Eta();
          phiTriggPool=t->Phi();
          idTriggPool = iTrackPool;
        }
      }
      lem[0] = fTracks->GetEntries(); 
      lem[1] = ptTriggPool; 
      lem[2] = etaTriggPool;
      lem[3] = phiTriggPool;
      
//       phiTriggPool = fCorrelator->SetCorrectPhiRange(phiTriggPool);
      fCorrelator->SetTriggerParticleProperties(ptTriggPool, phiTriggPool, etaTriggPool, idTriggPool);
        
        for(Int_t iTrack = 0; iTrack<NofTracks; iTrack++){//second loop on track candidates
        if(idTriggPool!=iTrack){
          Bool_t runcorrelation = fCorrelator->Correlate(iTrack);
          if(!runcorrelation) continue;
          Double_t DeltaPhi = fCorrelator->GetDeltaPhi(); lem[5] = fCorrelator->GetDeltaPhi();
          Double_t DeltaEta = fCorrelator->GetDeltaEta(); lem[4] = fCorrelator->GetDeltaEta();
//           cout<<"33333333333333333322222222222222222222"<<endl;
          if(Lem) Lem->Fill(lem);
          if(DebugLevel()>0){
//             cout<<"5555555555555555555555555555555"<<endl;
          (*AliMESbaseTask::DebugStream()) << "Correlate"
//             <<"ptTrigg=" << ptTrigg
//             <<"etaTrigg=" << etaTrigg
//             <<"phiTrigg=" << phiTrigg
            <<"ptTriggPool=" << ptTriggPool
            <<"etaTriggPool=" << etaTriggPool
            <<"phiTriggPool=" << phiTriggPool
            <<"DeltaEta=" << DeltaEta
            <<"DeltaPhi=" << DeltaPhi
            << "\n";
          }
        }
        }
      }
  }
  Bool_t updated = fCorrelator->PoolUpdate();
  if(!updated) AliInfo("Pool was not updated");
}
//________________________________________________________________________
Bool_t AliMESppColTask::PostProcess()
{
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliMESppColTask::AliMESppColMixEvent::DefineEventPool()
{
  
  Int_t MaxNofEvInfos = 200;//1000;
  Int_t MinNofTraks = 1000;//2000;
  Int_t nMultiplicityBins = 9;
  Double_t multiplicityBins[] = { 0, 6, 12, 19, 28, 39, 49, 59, 71, 82 };
  Int_t nZvtxBins = 7;
  Double_t vertexBins[] = {  -10, -5, -2.5,0, 2.5, 5, 10 };
  
  fPoolManager = new AliEventPoolManager(MaxNofEvInfos,MinNofTraks,nMultiplicityBins,(Double_t*)multiplicityBins,nZvtxBins,(Double_t*) vertexBins);
  
  if(!fPoolManager) return kFALSE;
  
  return kTRUE;

}

Bool_t AliMESppColTask::AliMESppColMixEvent::InitializeEventPool()
{
  
  Int_t mult = fParent->fEvInfo->GetMultiplicity(AliMESeventInfo::kComb); 
  Double_t zvertex = fParent->fEvInfo->GetVertexZ(); // zvertex
  fPool = fPoolManager->GetEventPool(mult, zvertex);
  fPool->SetTargetEvents(8); //set the minimum number of events to mix
  
  if(TMath::Abs(zvertex)>=10 || mult>82 || mult<0){
    AliInfo(Form("Event with zvertex = %5.2f cm and multiplicty = %d is out of pool bins, SKIPPING", zvertex, mult));
    return kFALSE;
  }
  
  fPool = fPoolManager->GetEventPool(mult, zvertex);
  
  if(!fPool){
    AliInfo(Form("No pool found for multiplicity = %d and zvertex = %f cm", mult, zvertex));
    return kFALSE;
  }
  
  fPool->PrintInfo();
  
  return kTRUE;
  
}

Bool_t AliMESppColTask::AliMESppColMixEvent::ProcessEventPool()
{
  
  if(!fmixing) return kFALSE;
  if(!fPool->IsReady()) return kFALSE;
//   AliInfo("Processing Event Pool");
  fPoolContent = fPool->GetCurrentNEvents();
  fPool->PrintInfo();
  return kTRUE;
  
}

Bool_t AliMESppColTask::AliMESppColMixEvent::Correlate(Int_t loopindex)
{
  if(loopindex >= fNTracks) return kFALSE;
  if(!fAssociatedTracks) return kFALSE;
  
  AliMEStrackInfo * fMixedTrack = (AliMEStrackInfo*)fAssociatedTracks->At(loopindex);
  fDeltaPhi = SetCorrectPhiRange(fPhiTrigg - fMixedTrack->Phi());        
  fDeltaEta = fEtaTrigg - fMixedTrack->Eta();
//   AliInfo("Corelare...");
  return kTRUE;

}

Bool_t AliMESppColTask::AliMESppColMixEvent::ProcessAssociatedTracks(Int_t EvLoopIndex)
{
  //associatedTracks should be deleted in the user task
  fAssociatedTracks = new TObjArray();
  
  if(!fmixing){ //Single Event analysis
    fAssociatedTracks = AcceptAndReduceTracks();
  }
  
  
  if(fmixing){ //Mixed Event analysis
    fAssociatedTracks = fPool->GetEvent(EvLoopIndex);
//     AliInfo("Getting event from pool...");
  }
  if(!fAssociatedTracks) return kFALSE;
  fNTracks = fAssociatedTracks->GetEntriesFast();
  
  return kTRUE;
  
}

Bool_t AliMESppColTask::AliMESppColMixEvent::PoolUpdate()
{
 
  if(!fmixing) return kFALSE;
  if(!fPool) return kFALSE;
  if(fmixing){ //the pool will be updated for Event Mixing procedure
//       AliInfo("Pool is updating...");
    TObjArray* objArr = NULL;
    objArr = (TObjArray*)AcceptAndReduceTracks();
    if(objArr->GetEntriesFast()>0) fPool->UpdatePool(objArr); //this ensures that the pool is updated only if there are entries in the array
  }
  
  return kTRUE;
  
}


Double_t AliMESppColTask::AliMESppColMixEvent::SetCorrectPhiRange(Double_t phi){
  
  Double_t pi = TMath::Pi();
  
  if(phi<-0.5*pi) phi+=2*pi;
  if(phi>pi+0.5*pi) phi-=2*pi;
  
  return phi;
  
}

TObjArray*  AliMESppColTask::AliMESppColMixEvent::AcceptAndReduceTracks()
{
  
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  AliMEStrackInfo *track(NULL);
   
  for (Int_t iTrack=0; iTrack<fParent->fTracks->GetEntries(); ++iTrack) {
      if(!(track = (AliMEStrackInfo*)fParent->fTracks->At(iTrack))) continue;
      if(TMath::Abs(track->Y()) > 1.) continue;
    
      
      tracksClone->Add(new AliMEStrackInfo(*track));
  }
  return tracksClone;
}

// TH2D* AliMESppColTask::AliMESppColMixEvent::NormToPeak(TH2D * inputHisto)
// {
// // function that normalizes the SE distribution to the bin centered at (delta phi, delta eta) = (0,0)
// Int_t* centralBins = FindCentralBin(inputHisto);
// inputHisto->Scale(1./inputHisto->GetBinContent(centralBins[0], centralBins[1]));
// 
// return inputHisto;
// 
// }

// TH2D* AliMESppColTask::AliMESppColMixEvent::ApplyEventMixingCorrection(TH2D * SEhisto, TH2D * MEhisto)
// //function that applis the event mixing correction
// {
//   TH2D * outputHisto;
//   NormToPeak(MEhisto);
//   
//   outputHisto = Divide2DHistos(SEhisto, MEhisto);
//   
//   return outputHisto;
// }
// 
// TH2D* AliMESppColTask::AliMESppColMixEvent::Divide2DHistos(TH2D * NumHisto, TH2D * DenomHisto)
// {
//   Int_t binsX = NumHisto->GetNbinsX();
//   Int_t binsY = NumHisto->GetNbinsY();
//   
//   Double_t XlowerBin = NumHisto->GetXaxis()->GetBinLowEdge(1);
//   Double_t XupperBin = NumHisto->GetXaxis()->GetBinLowEdge(binsX + 1);
//   
//   Double_t YlowerBin = NumHisto->GetYaxis()->GetBinLowEdge(1);
//   Double_t YupperBin = NumHisto->GetYaxis()->GetBinLowEdge(binsY + 1);
//   
//   if((binsX != DenomHisto->GetNbinsX()) || (binsY != DenomHisto->GetNbinsY())) printf(">>> Warning! Dividing two histos with different binning!");
//   
//   if((XlowerBin != DenomHisto->GetXaxis()->GetBinLowEdge(1)) || (XupperBin != DenomHisto->GetXaxis()->GetBinLowEdge(binsX + 1))) printf(">>> Warning! Dividing two histos with different X axis ranges!");
// 
//   if((YlowerBin != DenomHisto->GetYaxis()->GetBinLowEdge(1)) || (YupperBin != DenomHisto->GetYaxis()->GetBinLowEdge(binsY+1))) printf(">>> Warning! Dividing two histos with different Y axis ranges!");
// 
//   TString name = NumHisto->GetName();
//   name += "_ratio";
// // define the output histo
//   TH2D * outputHisto = new TH2D(name.Data(),name.Data(),binsX,XlowerBin,XupperBin,binsY,YlowerBin,YupperBin);
//   
//   Double_t ratio = 0;
//   Double_t ratioerr = 0;
//   Double_t numvalue = 0; Double_t denomvalue = 0;
//   Double_t numvalerr = 0; Double_t denomvalerr = 0;
//   
//   for (Int_t x =1; x<binsX+1; x++){ // loop on delta phi
//     for (Int_t y=1; y<binsY+1; y++){ // loop on delta eta
//       
//       numvalue = NumHisto->GetBinContent(x,y);
//       denomvalue = DenomHisto->GetBinContent(x,y);
//       numvalerr = NumHisto->GetBinError(x,y);
//       denomvalerr = DenomHisto->GetBinError(x,y);
//       
//       if(!denomvalue) {
//         printf("Error: Dividing by zero - cannot divide histos\n"); return NULL;
//       }
//       
//       ratio = numvalue/denomvalue;
//       ratioerr = TMath::Sqrt((numvalerr/denomvalue)*(numvalerr/denomvalue) + ratio*ratio * (denomvalerr/denomvalue)*(denomvalerr/denomvalue));
//       
//       outputHisto->SetBinContent(x,y,ratio);
//       outputHisto->SetBinError(x,y,ratioerr);
//     }
//   }
//   
//   return outputHisto;
//   
// }


//________________________________________________________
Bool_t AliMESppColTask::BuildQAHistos()
{
  fHistosQA = new TList(); fHistosQA->SetOwner(kTRUE);
  
  TString st1, st2;
  THnSparseI *L(NULL);
  if(!(L = (THnSparseI*)gROOT->FindObject("LeadingInfo"))){
    const Int_t ndim(6);
    const Char_t *cldTitle[ndim] = {"ntrk", "p_{t}", "#eta", "#phi", "#Delta#eta", "#Delta#phi"};
    const Int_t cldNbins[ndim]   = {  30, 100,     100,    180,     100,        120};
    const Double_t cldMin[ndim]  = {   0.5, 0.,  -1.0,     0.,    -1.0,     -TMath::Pi()/2.},
                   cldMax[ndim]  = {   30.5, 10.,     1., (2*TMath::Pi()), 1.,  TMath::Pi()*1.5};
    st1 = "Leading Info;";
    for(Int_t idim(0); idim<ndim; idim++){ st1 += cldTitle[idim]; st1+=";";}
    L = new THnSparseI("LeadingInfo", st1.Data(), ndim, cldNbins, cldMin, cldMax);
  } else L->Reset();
  fHistosQA->AddAt(L, 0);
  
  THnSparseI *Lem(NULL);
  if(!(Lem = (THnSparseI*)gROOT->FindObject("EventMixing"))){
    const Int_t ndim2(6);
    const Char_t *cldTitle2[ndim2] = {"ntrk", "p_{t}", "#eta", "#phi", "#Delta#eta", "#Delta#phi"};
    const Int_t cldNbins2[ndim2]   = {30, 100,     100,    180,     100,        120};
    const Double_t cldMin2[ndim2]  = {0.5, 0.,  -1.0,     0.,    -1.0,     -TMath::Pi()/2.},
                   cldMax2[ndim2]  = {30.5, 10.,     1., (2*TMath::Pi()), 1.,  TMath::Pi()*1.5};
    st2 = "Event Mixing Info;";
    for(Int_t idim(0); idim<ndim2; idim++){ st2 += cldTitle2[idim]; st2+=";";}
    Lem = new THnSparseI("EventMixing", st2.Data(), ndim2, cldNbins2, cldMin2, cldMax2);
  } else Lem->Reset();
  fHistosQA->AddAt(Lem, 1);
   
  return kTRUE;
}

