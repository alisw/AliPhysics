//
// This class is the task for connecting together 
// MC information and the RC information 
//
// The task is a wrapper over two components
// AliGenInfoMaker
// AliESDRecInfoMaker.h

// ROOT includes
#include <TChain.h>
#include <TMath.h>

// ALIROOT includes
#include <TTreeStream.h>
#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"

#include <AliESD.h>
#include "AliGenInfoTask.h"
#include "AliGenInfoMaker.h"
#include "AliHelix.h"

//
#include "AliMCInfo.h"
#include "AliComparisonObject.h"
#include "AliESDRecInfo.h"
#include "AliTPCParamSR.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TFile.h"
#include "AliTPCseed.h"

// STL includes
#include <iostream>

using namespace std;

ClassImp(AliGenInfoTask)

//________________________________________________________________________
AliGenInfoTask::AliGenInfoTask() : 
  AliAnalysisTask(), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fESDfriend(0),
  fCompList(0),         //array of comparison objects
  fGenTracksArray(0),  //clones array with filtered particles
  fGenKinkArray(0),    //clones array with filtered Kinks
  fGenV0Array(0),      //clones array with filtered V0s
  fRecTracksArray(0),  //clones array with filtered particles
  fDebugStreamer(0),
  fStreamLevel(0),
  fDebugLevel(0),
  fDebugOutputPath()
{
  //
  // Default constructor (should not be used)
  //
}

AliGenInfoTask::AliGenInfoTask(const AliGenInfoTask& /*info*/) : 
  AliAnalysisTask(), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fESDfriend(0),
  fCompList(0),
  fGenTracksArray(0),  //clones array with filtered particles
  fGenKinkArray(0),    //clones array with filtered Kinks
  fGenV0Array(0),      //clones array with filtered V0s
  fRecTracksArray(0),  //clones array with filtered particles
  //
  fDebugStreamer(0),
  fStreamLevel(0),
  fDebugLevel(),
  fDebugOutputPath()
{
  //
  // Default constructor 
  //
}



//________________________________________________________________________
AliGenInfoTask::AliGenInfoTask(const char *name) : 
  AliAnalysisTask(name, "AliGenInfoTask"), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fESDfriend(0),
  fCompList(0),
  fGenTracksArray(0),  //clones array with filtered particles
  fGenKinkArray(0),    //clones array with filtered Kinks
  fGenV0Array(0),      //clones array with filtered V0s
  fRecTracksArray(0),  //clones array with filtered particles
  fDebugStreamer(0),
  fStreamLevel(0),
  fDebugLevel(0),
  fDebugOutputPath()
{
  //
  // Normal constructor
  //
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList
  DefineOutput(0, TObjArray::Class());
  //
  //
  fCompList = new TObjArray;
}

AliGenInfoTask::~AliGenInfoTask(){
  //
  //
  //
  if (fDebugLevel>0)  printf("AliGenInfoTask::~AliGenInfoTask\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer=0;
  if(fCompList)   delete fCompList;  
  fCompList =0; 
}


//________________________________________________________________________
void AliGenInfoTask::ConnectInputData(Option_t *) 
{
  //
  // Connect the input data
  //
  if(fDebugLevel>3)
    cout << "AnalysisTaskTPCCluster::ConnectInputData()" << endl;

  TTree* tree=dynamic_cast<TTree*>(GetInputData(0));
  if (!tree) {
    //Printf("ERROR: Could not read chain from input slot 0");
  }
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!esdH) {
      //Printf("ERROR: Could not get ESDInputHandler");
    }
    else {
      fESD = esdH->GetEvent();
      fESDfriend=static_cast<AliESDfriend*>(fESD->FindListObject("AliESDfriend"));
      fESD->SetESDfriend(fESDfriend);
      //Printf("*** CONNECTED NEW EVENT ****");
    }  
  }
  AliMCEventHandler* mcinfo = (AliMCEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());  
  mcinfo->SetReadTR(kTRUE);
  
  fMCinfo = mcinfo->MCEvent();


}


//_____________________________________________________________________________
Bool_t AliGenInfoTask::AddComparisonObject(AliComparisonObject *pObj) 
{
  // add comparison object to the list
  if(pObj == 0) {
      Printf("ERROR: Could not add comparison object");
	  return kFALSE;
  }
  // add object to the list
  fCompList->AddLast(pObj);       
  return kTRUE;
}




//________________________________________________________________________
void AliGenInfoTask::CreateOutputObjects() 
{
  //
  // Connect the output objects
  //
  if(fDebugLevel>3)
    cout << "AnalysisTaskTPCCluster::CreateOutputObjects()" << endl;
  
}


//________________________________________________________________________
void AliGenInfoTask::Exec(Option_t *) {
  //
  // Execute analysis for current event 
  //

  if(fDebugLevel>3)
    cout << "AliGenInfoTask::Exec()" << endl;
    

  // If MC has been connected   
  if (fGenTracksArray) fGenTracksArray->Delete();
  if (fRecTracksArray) fRecTracksArray->Delete();

  if (!fMCinfo){
    cout << "Not MC info\n" << endl;
  }else{
    //mcinfo->Print();
    ProcessMCInfo();
    ProcessESDInfo();
    DumpInfo();
    ProcessComparison();
  }
  //
}      




//________________________________________________________________________
void AliGenInfoTask::Terminate(Option_t *) {
    //
    // Terminate loop
    //
  if(fDebugLevel>3)
    printf("AliGenInfoTask: Terminate() \n");  
  //
  if (fDebugLevel>0) printf("AliTPCcalibBase::Terminate\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer = 0;
  return;
}





TTreeSRedirector *AliGenInfoTask::GetDebugStreamer(){
  //
  // Get Debug streamer
  // In case debug streamer not yet initialized and StreamLevel>0 create new one
  //
  if (fStreamLevel==0) return 0;
  if (fDebugStreamer) return fDebugStreamer;
  TString dsName;
  dsName=GetName();
  dsName+="Debug.root";
  dsName.ReplaceAll(" ","");
  fDebugStreamer = new TTreeSRedirector(dsName.Data());
  return fDebugStreamer;
}



AliMCInfo*  AliGenInfoTask::GetTrack(Int_t index, Bool_t force){
  //
  // Get the MC info for given track
  //
  if (!fGenTracksArray) fGenTracksArray = new TClonesArray("AliMCInfo",1000);
  if (index>fGenTracksArray->GetEntriesFast()) fGenTracksArray->Expand(index*2+10);
  AliMCInfo * info = (AliMCInfo*)fGenTracksArray->At(index);
  if (!force) return info;
  if (!info){
    info = new ((*fGenTracksArray)[index]) AliMCInfo;
  }
  return info;
}

AliESDRecInfo*  AliGenInfoTask::GetRecTrack(Int_t index, Bool_t force){
  //
  // Get the MC info for given track
  //
  if (!fRecTracksArray) fRecTracksArray = new TClonesArray("AliESDRecInfo",1000);
  if (index>fRecTracksArray->GetEntriesFast()) fRecTracksArray->Expand(index*2+10);
  AliESDRecInfo * info = (AliESDRecInfo*)fRecTracksArray->At(index);
  if (!force) return info;
  if (!info){
    info = new ((*fRecTracksArray)[index]) AliESDRecInfo;
  }
  return info;
}




void  AliGenInfoTask::ProcessMCInfo(){
  //
  // Dump information from MC to the array
  //
  //
  TParticle * particle= new TParticle;
  TClonesArray * trefs = new TClonesArray("AliTrackReference");
  //
  //
  // Process tracks
  //
  Int_t npart = fMCinfo->GetNumberOfTracks();
  if (npart==0) return;
  Double_t vertex[4]={0,0,0,0};
  fMCinfo->GetParticleAndTR(0, particle, trefs);
  if (particle){
    vertex[0]=particle->Vx();
    vertex[1]=particle->Vy();
    vertex[2]=particle->Vz();
    vertex[3]=particle->R();
  }

  for (Int_t ipart=0;ipart<npart;ipart++){
    Int_t status = fMCinfo->GetParticleAndTR(ipart, particle, trefs);
    if (status<0) continue;
    if (!particle) continue;
    if (!trefs) continue;
    if (!AcceptParticle(particle)) continue;
    //if (trefs->GetEntries()<1) continue;
    AliMCInfo * mcinfo = GetTrack(ipart,kTRUE);
    mcinfo->Update(particle,trefs,vertex,ipart);
    //
    TTreeSRedirector *pcstream = GetDebugStreamer();
    if (pcstream){
      (*pcstream)<<"MC"<<
	"p.="<<particle<<
	"MC.="<<mcinfo<<
	"\n";
    }
  }
}

void AliGenInfoTask::ProcessESDInfo(){
  //
  //
  //
  if (!fESD) return;
  static AliTPCParamSR param;
  fESDfriend=static_cast<AliESDfriend*>(fESD->FindListObject("AliESDfriend"));
  if (!fESDfriend) {
    //Printf("ERROR: fESDfriend not available");
    return;
  }
  fESD->SetESDfriend(fESDfriend);
  //
  //
  Int_t ntracks = fESD->GetNumberOfTracks();
  for (Int_t itrack=0; itrack<ntracks; itrack++){
    AliESDtrack *track = fESD->GetTrack(itrack);
    Int_t label = TMath::Abs(track->GetLabel());
    AliMCInfo * mcinfo = GetTrack(label,kFALSE);
    if (!mcinfo) continue;
    AliESDRecInfo *recInfo= GetRecTrack(label,kTRUE);
    recInfo->AddESDtrack(track,mcinfo);
    recInfo->Update(mcinfo,&param,kTRUE);
  }
  //


} 


void AliGenInfoTask::ProcessComparison(){
  //
  //
  //
  static AliESDRecInfo dummy;
  Int_t npart = fMCinfo->GetNumberOfTracks();
  for (Int_t ipart=0;ipart<npart;ipart++){
    AliMCInfo * mcinfo = GetTrack(ipart,kFALSE);
    if (!mcinfo) continue;
    AliESDRecInfo *recInfo= GetRecTrack(ipart,kFALSE);
    if (!recInfo) recInfo=&dummy; 
    //
    for (Int_t icomp = 0; icomp<fCompList->GetEntries(); icomp++){
      AliComparisonObject *pObj= (AliComparisonObject *)fCompList->At(icomp);
      if (pObj){
	pObj->Exec(mcinfo,recInfo);
      }
    }    
  }
  PostData(0, fCompList);
}

void AliGenInfoTask::DumpInfo(){
  //
  //
  //
  TParticle * particle= new TParticle;
  TClonesArray * trefs = new TClonesArray("AliTrackReference");
  //

  static AliESDRecInfo dummy;
  Int_t npart = fMCinfo->GetNumberOfTracks();
  for (Int_t ipart=0;ipart<npart;ipart++){
    AliMCInfo * mcinfo = GetTrack(ipart,kFALSE);
    if (!mcinfo) continue;
    AliESDRecInfo *recInfo= GetRecTrack(ipart,kFALSE);
    if (!recInfo) recInfo=&dummy; 
    TTreeSRedirector *pcstream = GetDebugStreamer();
    
    fMCinfo->GetParticleAndTR(ipart, particle, trefs);
    Int_t counter=0;
    Float_t length=0;
    if (trefs!=0 && particle!=0){
      length = GetTPCTrackLength(*trefs, particle , fESD->GetMagneticField(), counter, 3.0);
    }

    if (pcstream){
      (*pcstream)<<"RC"<<	
	"MC.="<<mcinfo<<
	"RC.="<<recInfo<<
	"length="<<length<<
	"counter="<<counter<<
	//	"tr.="<<tpctrack<<
	"\n";
    }    
  }
}



Bool_t AliGenInfoTask::AcceptParticle(TParticle *part){
  //
  /*
    MC cuts
    TCut cutPt("p.Pt()>0.1");
    TCut cutZ("abs(p.Vz())<250");
    TCut cutR("abs(p.R())<250");
    //
  */
  //
  //
  if (part->Pt()<0.1) return kFALSE;
  if (TMath::Abs(part->Vz())>250) return kFALSE;
  if (part->R()>360)  return kFALSE;
  if (part->GetPDG()){
    if (part->GetPDG()->Charge()==0) return kFALSE;
  }
  return kTRUE;
}



void AliGenInfoTask::FinishTaskOutput()
{
  //
  // According description in AliAnalisysTask this method is call
  // on the slaves before sending data
  //
  Terminate("slave");
  RegisterDebugOutput(fDebugOutputPath.Data());

}



void AliGenInfoTask::RegisterDebugOutput(const char */*path*/){
  //
  // store  - copy debug output to the destination position
  // currently ONLY for local copy
  TString dsName;
  dsName=GetName();
  dsName+="Debug.root";
  dsName.ReplaceAll(" ","");
  TString dsName2=fDebugOutputPath.Data();
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+="/";;
  dsName2+=gSystem->HostName();
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+="/";
  dsName2+=gSystem->BaseName(gSystem->pwd());
  dsName2+="/";
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+=dsName;
  AliInfo(Form("copy %s\t%s\n",dsName.Data(),dsName2.Data()));
  printf("copy %s\t%s\n",dsName.Data(),dsName2.Data());
  TFile::Cp(dsName.Data(),dsName2.Data());
}


Float_t  AliGenInfoTask::GetTPCTrackLength(const TClonesArray& trackRefs, TParticle*part, Float_t bz,  Int_t &counter, Float_t deadWidth){
  //
  // return track length in geometrically active volume of TPC.
  // z nad rphi acceptance is included
  // doesn't take into account dead channel and ExB  
  // Intput:
  // trackRefs
  // bz - magnetic field
  // deadWidth - dead zone in r-phi
  // Additional output:
  // counter   - number of circles
  const Float_t kRMin = 90;
  const Float_t kRMax = 245;
  const Float_t kZMax = 250;
  const Float_t kMinPt= 0.01; 
  Float_t length =0;
  Int_t nrefs = trackRefs.GetEntriesFast();
  AliExternalTrackParam param;
  Double_t cv[21];
  for (Int_t i=0; i<21;i++) cv[i]=0;
  counter=0;
  //
  //
  AliTrackReference *ref0 = (AliTrackReference*)trackRefs.At(0);
  Float_t direction = 0;
  //
  for (Int_t iref=1; iref<nrefs;iref++){
    AliTrackReference *ref = (AliTrackReference*)trackRefs.At(iref);
    if (!ref) continue;
    if (!ref0 || ref0->DetectorId()!= AliTrackReference::kTPC){
      ref0 = ref;
      direction = ((ref0->X()*ref0->Px()+ref0->Y()*ref0->Py())>0)? 1.:-1.;
      continue;
    }
    Float_t newdirection = ((ref->X()*ref->Px()+ref->Y()*ref->Py())>0)? 1.:-1.;
    if (newdirection*direction<0) {
      counter++;  //circle counter 
      direction = newdirection;
      continue;
    }
    if (counter>0) continue;
    if (ref0->Pt()<kMinPt) break;
    Float_t radius0 = TMath::Max(TMath::Min(ref0->R(),kRMax),kRMin);;
    Float_t radius1 = TMath::Max(TMath::Min(ref->R(),kRMax),kRMin);
    Double_t xyz[3]={ref0->X(),ref0->Y(),ref0->Z()};
    Double_t pxyz[3]={ref0->Px(),ref0->Py(),ref0->Pz()};
    Double_t alpha;
       param.Set(xyz,pxyz,cv,TMath::Nint(part->GetPDG()->Charge()/3.));
    
    for (Float_t radius = radius0; radius<radius1; radius+=1){
      param.GetXYZAt(radius, bz, xyz);
      if (TMath::Abs(xyz[2])>kZMax) continue;
      Float_t gradius = TMath::Sqrt(xyz[1]*xyz[1]+xyz[0]*xyz[0]);
      if (gradius>kRMax) continue;
      alpha = TMath::ATan2(xyz[1],xyz[0]);
      if (alpha<0) alpha+=TMath::TwoPi();
      //
      Int_t sector  = Int_t(9*alpha/TMath::Pi());
      Float_t lalpha = alpha-((sector+0.5)*TMath::Pi()/9.);
      Float_t dedge  = (TMath::Tan(TMath::Pi()/18.)-TMath::Abs(TMath::Tan(lalpha)))*gradius; 
      if (dedge>deadWidth) length++;
    }
    if (ref->DetectorId()!= AliTrackReference::kTPC) break; 
    ref0 = ref;
  }
  return length;
}
