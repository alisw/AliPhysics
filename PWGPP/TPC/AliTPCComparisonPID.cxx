//
//
// 

//
// ROOT includes
#include <TChain.h>
#include <TMath.h>
#include <TVectorD.h>
#include <TSystem.h>
#include <TFile.h>
#include <TParticle.h>

// ALIROOT includes
#include <TTreeStream.h>
#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMathBase.h"

#include <AliESD.h>
#include "AliExternalTrackParam.h"
#include "AliTracker.h"
#include "AliTPCseed.h"
//
#include "AliTPCComparisonPID.h"
//
#include <THnSparse.h>

//

// STL includes
#include <iostream>

using namespace std;

ClassImp(AliTPCComparisonPID)

//________________________________________________________________________
AliTPCComparisonPID::AliTPCComparisonPID() : 
  AliAnalysisTask(), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fTPCsignal(0),
  fTPCsignalNorm(0),
  fDebugStreamer(0),
  fStreamLevel(0),
  fDebugLevel(0),
  fDebugOutputPath()
{
  //
  // Default constructor (should not be used)
  //
}

AliTPCComparisonPID::AliTPCComparisonPID(const AliTPCComparisonPID& info) : 
  AliAnalysisTask(info), 
  fMCinfo(info.fMCinfo),     //! MC event handler
  fESD(info.fESD),        //!
  fTPCsignal(0),
  fTPCsignalNorm(0),
  //
  fDebugStreamer(0),
  fStreamLevel(0),
  fDebugLevel(),
  fDebugOutputPath()
{
  //
  // Dummy Copy  constructor - no copy constructor for THnSparse 
  //
}



//________________________________________________________________________
AliTPCComparisonPID::AliTPCComparisonPID(const char *name) : 
  AliAnalysisTask(name, "AliTPCComparisonPID"), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fTPCsignal(0),
  fTPCsignalNorm(0),
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
  DefineOutput(0, AliTPCComparisonPID::Class());
  //
  //make histos
  Init(); 
}

void AliTPCComparisonPID::Init(){
  //
  // Init dEdx histogram
  // Dimensions
  // 0 - particle specie as defined in the AliPID - negatives+5 <0,9>
  // 1 - momenta - at the entrance of the TPC
  // 2 - tan lambda- fP[3]
  // 3 - betagamma
  // 4 - measurement - dEdx or dEdx/BB
  //
  Double_t xmin[5],  xmax[5];
  Int_t    nbins[5];
  // pid
  nbins[0]=10;
  xmin[0]=0; xmax[0]=10;
  // momenta
  nbins[1]=30;
  xmin[1]=0.1; xmax[1]=3;
  //P3
  nbins[2]=20;
  xmin[2]=-1.5; xmax[2]=1.5;
  //
  // log (betagamma)
  //
  nbins[3]=50;
  xmin[3]=0.1; xmax[3]=100;
  //
  // 
  nbins[4]=400;
  xmin[4]=20; xmax[4]=400;
  fTPCsignal = new THnSparseF("TPC signal","TPC signal",5,nbins,xmin,xmax);
  nbins[4]=100;
  xmin[4]=25; xmax[4]=75;
  fTPCsignal = new THnSparseF("TPC signal Norm","TPC signal Norm",5,nbins,xmin,xmax);
  //
}





AliTPCComparisonPID::~AliTPCComparisonPID(){
  //
  //
  //
  if (fDebugLevel>0)  printf("AliTPCComparisonPID::~AliTPCComparisonPID\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer=0;
  delete fTPCsignal;
  delete fTPCsignalNorm;
}


//________________________________________________________________________
void AliTPCComparisonPID::ConnectInputData(Option_t *) 
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
      //Printf("*** CONNECTED NEW EVENT ****");
    }  
  }
  AliMCEventHandler* mcinfo = (AliMCEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());  
  mcinfo->SetReadTR(kTRUE);
  
  fMCinfo = mcinfo->MCEvent();


}






//________________________________________________________________________
void AliTPCComparisonPID::CreateOutputObjects() 
{
  //
  // Connect the output objects
  //
  if(fDebugLevel>3)
    cout << "AnalysisTaskTPCCluster::CreateOutputObjects()" << endl;

}


//________________________________________________________________________
void AliTPCComparisonPID::Exec(Option_t *) {
  //
  // Execute analysis for current event 
  //

  if(fDebugLevel>3)
    cout << "AliTPCComparisonPID::Exec()" << endl;
    

  // If MC has been connected   

  if (!fMCinfo){
    cout << "Not MC info\n" << endl;
  }else{
    ProcessMCInfo();
    //mcinfo->Print();
    //DumpInfo();
  }
  //
  PostData(0, this);
}      




//________________________________________________________________________
void AliTPCComparisonPID::Terminate(Option_t *) {
    //
    // Terminate loop
    //
  if(fDebugLevel>3)
    printf("AliTPCComparisonPID: Terminate() \n");  
  //
  if (fDebugLevel>0) printf("AliMCtrackingTestTask::Terminate\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer = 0;
  return;
}



TTreeSRedirector *AliTPCComparisonPID::GetDebugStreamer(){
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




void  AliTPCComparisonPID::ProcessMCInfo(){
  //
  //
  //
  //
  Int_t npart   = fMCinfo->GetNumberOfTracks();
  Int_t ntracks = fESD->GetNumberOfTracks(); 
  if (npart<=0) return;
  if (ntracks<=0) return;
  //
  //
  TParticle * particle= new TParticle;
  TClonesArray * trefs = new TClonesArray("AliTrackReference");
  
  for (Int_t itrack=0;itrack<ntracks;itrack++){
    AliESDtrack *track = fESD->GetTrack(itrack);
    const AliExternalTrackParam *in=track->GetInnerParam();
    if (!in) continue;
    Int_t ipart = TMath::Abs(track->GetLabel());
    //
    Int_t status = fMCinfo->GetParticleAndTR(ipart, particle, trefs);
    if (status<0) continue;
    if (!particle) continue;
    if (!trefs) continue;
    //
    //
    Double_t mom = in->GetP();
    Double_t dedx=track->GetTPCsignal();
    Double_t mass = particle->GetMass();
    Double_t bg  =mom/mass;
    Double_t betheBloch = AliMathBase::BetheBlochAleph(bg);
    //
    // Fill histos
    //
    Double_t x[5];
    //PID
    Int_t pdg = particle->GetPdgCode();
    for (Int_t iType=0;iType<5;iType++) {
      if (AliPID::ParticleCode(iType)==TMath::Abs(pdg)){
	x[0]=iType;
	if (pdg<0) x[0]+=5;
      }
    }
    x[1]= mom;
    x[2]= track->GetTgl();
    x[3]= TMath::Log(bg);
    x[4]= dedx;
    fTPCsignal->Fill(x);
    x[4]=dedx/betheBloch;
    fTPCsignalNorm->Fill(x);    
  }
}




void AliTPCComparisonPID::FinishTaskOutput()
{
  //
  // According description in AliAnalisysTask this method is call
  // on the slaves before sending data
  //
  Terminate("slave");
  gSystem->Exec("pwd");
  RegisterDebugOutput();

}


void AliTPCComparisonPID::RegisterDebugOutput(){
  //
  //
  //
  //
  // store  - copy debug output to the destination position
  // currently ONLY for local copy
  TString dsName;
  dsName=GetName();
  dsName+="Debug.root";
  dsName.ReplaceAll(" ","");
  TString dsName2=fDebugOutputPath.Data();
  gSystem->MakeDirectory(dsName2.Data());
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
