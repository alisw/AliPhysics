//
// This class is the task to check the 
// Propagation and Update method used in the 
//               1. AliExternalTrackParam 
//               2. AliTracker
//
// Pure Monte-Carlo data used, not influence of detectors

// Input - TParticle + Array of track references - (Points alogn track trajectories)
// Output - Trees with track references - no histograms
//          MC tree -  test for material budget correction 
//                     see function ProcessRefTracker
//          MCupdate tree - test for correctness of propagation and update
//                     see function AliMCTrackingTestTask::FitTrackRefs
//
// Principle - Creates AliExternalTrackParam form 1 Track Refernece - 
//             Propagate it to other
// Magnetic field and the geometry has to be created before using it 

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

// ALIROOT includes
#include <TTreeStream.h>
#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"

#include <AliESD.h>
#include "AliMCTrackingTestTask.h"
#include "AliGenInfoMaker.h"
#include "AliHelix.h"

//
#include "AliMCInfo.h"
#include "AliComparisonObject.h"
#include "AliESDRecInfo.h"
#include "AliTPCParamSR.h"
#include "AliTracker.h"
#include "AliTPCseed.h"

// STL includes
#include <iostream>

using namespace std;

ClassImp(AliMCTrackingTestTask)

//________________________________________________________________________
AliMCTrackingTestTask::AliMCTrackingTestTask() : 
  AliAnalysisTask(), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fDebugStreamer(0),
  fStreamLevel(0),
  fDebugLevel(0),
  fDebugOutputPath()
{
  //
  // Default constructor (should not be used)
  //
}

AliMCTrackingTestTask::AliMCTrackingTestTask(const AliMCTrackingTestTask& /*info*/) : 
  AliAnalysisTask(), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
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
AliMCTrackingTestTask::AliMCTrackingTestTask(const char *name) : 
  AliAnalysisTask(name, "AliMCTrackingTestTask"), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
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
}

AliMCTrackingTestTask::~AliMCTrackingTestTask(){
  //
  //
  //
  if (fDebugLevel>0)  printf("AliMCTrackingTestTask::~AliMCTrackingTestTask\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer=0;
}


//________________________________________________________________________
void AliMCTrackingTestTask::ConnectInputData(Option_t *) 
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
void AliMCTrackingTestTask::CreateOutputObjects() 
{
  //
  // Connect the output objects
  //
  if(fDebugLevel>3)
    cout << "AnalysisTaskTPCCluster::CreateOutputObjects()" << endl;
  
}


//________________________________________________________________________
void AliMCTrackingTestTask::Exec(Option_t *) {
  //
  // Execute analysis for current event 
  //

  if(fDebugLevel>3)
    cout << "AliMCTrackingTestTask::Exec()" << endl;
    

  // If MC has been connected   

  if (!fMCinfo){
    cout << "Not MC info\n" << endl;
  }else{
    ProcessMCInfo();
    //mcinfo->Print();
    //DumpInfo();
  }
  //
}      




//________________________________________________________________________
void AliMCTrackingTestTask::Terminate(Option_t *) {
    //
    // Terminate loop
    //
  if(fDebugLevel>3)
    printf("AliMCTrackingTestTask: Terminate() \n");  
  //
  if (fDebugLevel>0) printf("AliMCtrackingTestTask::Terminate\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer = 0;
  return;
}



TTreeSRedirector *AliMCTrackingTestTask::GetDebugStreamer(){
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






AliExternalTrackParam * AliMCTrackingTestTask::MakeTrack(const AliTrackReference* ref, TParticle*part)
{
  //
  // Make track out of the track ref
  // part - TParticle used to determine chargr
  // the covariance matrix - equal 0 - starting from ideal MC position
  Double_t xyz[3]={ref->X(),ref->Y(),ref->Z()};
  Double_t pxyz[3]={ref->Px(),ref->Py(),ref->Pz()};
  Double_t cv[21];
  for (Int_t i=0; i<21;i++) cv[i]=0;
  if (!part->GetPDG()) return 0;
  AliExternalTrackParam * param = new AliExternalTrackParam(xyz,pxyz,cv,TMath::Nint(part->GetPDG()->Charge()/3.));
  return param;
}

Bool_t  AliMCTrackingTestTask::PropagateToPoint(AliExternalTrackParam *param, Double_t *xyz, Double_t mass, Float_t step){
  // 
  // Propagate track to point xyz using 
  // AliTracker::PropagateTo functionality
  //
  //  param - track parameters
  //  xyz   - position to propagate
  //  mass  - particle mass
  //  step  - step to be used
  Double_t radius=TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  AliTracker::PropagateTrackTo(param, radius+step, mass, step, kTRUE,0.99);
  AliTracker::PropagateTrackTo(param, radius+0.5, mass, step*0.1, kTRUE,0.99);
  Double_t sxyz[3]={0,0,0};
  AliESDVertex vertex(xyz,sxyz);
  Bool_t isOK = param->PropagateToDCA(&vertex,AliTracker::GetBz(),10);
  return isOK;
}


void  AliMCTrackingTestTask::ProcessMCInfo(){
  //
  //
  //
   //
  TParticle * particle= new TParticle;
  TClonesArray * trefs = new TClonesArray("AliTrackReference");
  const Double_t kPcut=0.1;
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
  //
  //

  for (Int_t ipart=0;ipart<npart;ipart++){
    Int_t status = fMCinfo->GetParticleAndTR(ipart, particle, trefs);
    if (status<0) continue;
    if (!particle) continue;
    if (!trefs) continue;
    Int_t nref = trefs->GetEntries();
    if (nref<5) continue;
    FitTrackRefs(particle,trefs);
    AliTrackReference * tpcIn=0;
    AliTrackReference * tpcOut=0;
    AliTrackReference * trdIn=0;
    AliTrackReference * trdOut=0;
    AliTrackReference * itsIn=0;
    AliTrackReference * itsOut=0;
    Double_t rmax=0;
    Double_t rmin=1000;
    for (Int_t iref=0;iref<nref;iref++){
      AliTrackReference * ref = (AliTrackReference*)trefs->At(iref);
      if (!ref) continue;
      
      Float_t dir = ref->X()*ref->Px()+ref->Y()*ref->Py();
      
      if (dir<0) break; // oposite direction - looping track - return back
      if (ref->P()<kPcut) continue;
      if (ref->R()<rmax) break;
      //if (ref->R()<rmin)  break; 
      //TPC
      if (ref->DetectorId()==AliTrackReference::kTPC){
	if (!tpcIn) {
	  tpcIn  = ref;
	}else{
	  if (ref->R()>tpcIn->R()) tpcOut = ref;
	}	
      }
      //ITS
      if (ref->DetectorId()==AliTrackReference::kITS){
	if (!itsIn) {
	  itsIn  = ref;
	}else{
	  if (ref->R()>itsIn->R()) itsOut = ref;
	}	
      }
      //TRD
      if (ref->DetectorId()==AliTrackReference::kTRD){
	if (!trdIn) {
	  trdIn  = ref;
	}else{
	  if (ref->R()>trdIn->R()) trdOut = ref;
	}	
      }      
      if (ref->R()<rmin) rmin=ref->R();
      if (ref->R()>rmax) rmax=ref->R();
    }
    if (tpcIn && tpcOut) {
      ProcessRefTracker(tpcIn,tpcOut,particle,1);
      ProcessRefTracker(tpcIn,tpcOut,particle,3);
    }
    if (itsIn && itsOut) ProcessRefTracker(itsIn,itsOut,particle,0);
    if (trdIn && trdOut) ProcessRefTracker(trdIn,trdOut,particle,2);
  }
}



void AliMCTrackingTestTask::ProcessRefTracker(AliTrackReference* refIn,  AliTrackReference* refOut, TParticle*part,Int_t type){
  //
  // Test propagation from In to out
  //
  AliExternalTrackParam *param = 0;
  AliExternalTrackParam *paramMC = 0;
  Double_t xyzIn[3]={refIn->X(),refIn->Y(), refIn->Z()};
  Double_t mass = part->GetMass();
  Double_t step=1;
  //
  param=MakeTrack(refOut,part);
  paramMC=MakeTrack(refOut,part);
  if (!param) return;
  if (type<3) PropagateToPoint(param,xyzIn, mass, step);
  if (type==3) {
    AliTPCseed seed;
    seed.Set(param->GetX(),param->GetAlpha(),param->GetParameter(),param->GetCovariance());
    Float_t alpha= TMath::ATan2(refIn->Y(),refIn->X());
    if(seed.Rotate(alpha-seed.GetAlpha())==kFALSE) return;
    seed.SetMass(mass);
    for (Float_t xlayer= seed.GetX(); xlayer>refIn->R(); xlayer-=step){
      seed.PropagateTo(xlayer);
    }
    seed.PropagateTo(refIn->R());
    param->Set(seed.GetX(),seed.GetAlpha(),seed.GetParameter(),seed.GetCovariance());
  }
  TTreeSRedirector *pcstream = GetDebugStreamer();
  TVectorD gpos(3);
  TVectorD gmom(3);
  param->GetXYZ(gpos.GetMatrixArray());
  param->GetPxPyPz(gmom.GetMatrixArray());
  if (pcstream){
    (*pcstream)<<"MC"<<
      "type="<<type<<
      "step="<<step<<
      "refIn.="<<refIn<<
      "refOut.="<<refOut<<
      "p.="<<part<<
      "par.="<<param<<   
      "parMC.="<<paramMC<<   
      "gpos.="<<&gpos<<
      "gmom.="<<&gmom<<
      "\n";
  }
}


void  AliMCTrackingTestTask::FitTrackRefs(TParticle * part, TClonesArray * trefs){
  //
  //
  //
  //
  const Int_t kMinRefs=6;
  Int_t nrefs = trefs->GetEntries();
  if (nrefs<kMinRefs) return; // we should have enough references
  Int_t iref0 =-1;
  Int_t iref1 =-1;
  
  for (Int_t iref=0; iref<nrefs; iref++){
    AliTrackReference * ref = (AliTrackReference*)trefs->At(iref);
    if (!ref) continue;    
    Float_t dir = ref->X()*ref->Px()+ref->Y()*ref->Py();
    if (dir<0) break;
    if (ref->DetectorId()!=AliTrackReference::kTPC) continue;
    if (iref0<0) iref0 = iref;
    iref1 = iref;    
  }
  if (iref1-iref0<kMinRefs) return;
  Double_t covar[15];
  for (Int_t icov=0; icov<15; icov++) covar[icov]=0;
  covar[0]=1; 
  covar[2]=1; 
  covar[5]=1;
  covar[9]=1;
  covar[14]=1;

  AliTrackReference * refIn = (AliTrackReference*)trefs->At(iref0);
  AliTrackReference * refOut = (AliTrackReference*)trefs->At(iref1);
  AliExternalTrackParam *paramPropagate= MakeTrack(refIn,part);
  AliExternalTrackParam *paramUpdate   = MakeTrack(refIn,part);
  paramUpdate->AddCovariance(covar);
  Double_t mass = part->GetMass();
  Double_t charge = part->GetPDG()->Charge()/3.;
/*
  Float_t alphaIn= TMath::ATan2(refIn->Y(),refIn->X());
  Float_t radiusIn= refIn->R();
  Float_t alphaOut= TMath::ATan2(refOut->Y(),refOut->X());
  Float_t radiusOut= refOut->R();
*/
  Bool_t isOKP=kTRUE;
  Bool_t isOKU=kTRUE;
  AliMagF * field = (AliMagF*) TGeoGlobalMagField::Instance()->GetField();
  for (Int_t iref = iref0; iref<=iref1; iref++){
    AliTrackReference * ref = (AliTrackReference*)trefs->At(iref);
    Float_t alphaC= TMath::ATan2(ref->Y(),ref->X());
    Double_t pos[3] = {ref->X(), ref->Y(), ref->Z()};
    Double_t mag[3];
    field->Field(pos,mag);
    isOKP&=paramPropagate->Rotate(alphaC);
    isOKU&=paramUpdate->Rotate(alphaC);
    for (Float_t xref= paramPropagate->GetX(); xref<ref->R(); xref++){
      isOKP&=paramPropagate->PropagateTo(xref, mag[2]);
      isOKU&=paramUpdate->PropagateTo(xref, mag[2]);
    }
    isOKP&=paramPropagate->PropagateTo(ref->R(), mag[2]);
    isOKU&=paramUpdate->PropagateTo(ref->R(), mag[2]);
    Double_t clpos[2] = {0, ref->Z()};
    Double_t clcov[3] = { 0.005,0,0.005};
    isOKU&= paramUpdate->Update(clpos, clcov);  
  }
  TTreeSRedirector *pcstream = GetDebugStreamer();
  if (pcstream){
    TVectorD gposU(3);
    TVectorD gmomU(3);
    TVectorD gposP(3);
    TVectorD gmomP(3);
    paramUpdate->GetXYZ(gposU.GetMatrixArray());
    paramUpdate->GetPxPyPz(gmomU.GetMatrixArray());
    paramPropagate->GetXYZ(gposP.GetMatrixArray());
    paramPropagate->GetPxPyPz(gmomP.GetMatrixArray());

     (*pcstream)<<"MCupdate"<<
       "isOKU="<<isOKU<<
       "isOKP="<<isOKP<<
       "m="<<mass<<
       "q="<<charge<<
       "part.="<<part<<
       "refIn.="<<refIn<<
       "refOut.="<<refOut<<
       "pP.="<<paramPropagate<<
       "pU.="<<paramUpdate<<
       "gposU.="<<&gposU<<
       "gmomU.="<<&gmomU<<
       "gposP.="<<&gposP<<
       "gmomP.="<<&gmomP<<
       "\n";
   }
}




void AliMCTrackingTestTask::FinishTaskOutput()
{
  //
  // According description in AliAnalisysTask this method is call
  // on the slaves before sending data
  //
  Terminate("slave");
  gSystem->Exec("pwd");
  RegisterDebugOutput();

}


void AliMCTrackingTestTask::RegisterDebugOutput(){
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


/*
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool; 
  TChain * chain = tool.MakeChain("mctracking.txt","MC",0,100);
  chain->Lookup();
  //
  //
  chain->SetAlias("pdg","(p.fPdgCode)");
  chain->SetAlias("dPRec","(refOut.P()-par.P())/refIn.P()");
  chain->SetAlias("dPMC","(refOut.P()-refIn.P())/refIn->P()");
  chain->SetAlias("dPtRec","(refOut.Pt()-par.Pt())/refIn.Pt()");
  chain->SetAlias("dPtMC","(refOut.Pt()-refIn.Pt())/refIn->Pt()");


  // ITS
  chain->Draw("-sqrt(dPRec):-sqrt(dPMC)","abs(pdg)!=11&&type==0&&p.Pt()<0.5&&abs(p.R())<1&&abs(refOut.fZ)<220");
  htemp->SetYTitle("#sqrt{#DeltaP_{rec}/P}");
  htemp->SetXTitle("#sqrt{#DeltaP_{mc}/P}");
  gPad->SaveAs("picLoss/dPcorr_ITS_step1.gif");
  gPad->SaveAs("picLoss/dPcorr_ITS_step1.eps");
  // TPC
  chain->Draw("-sqrt(dPRec):-sqrt(dPMC)","abs(pdg)!=11&&type==1&&p.Pt()<0.5&&abs(p.R())<1&&abs(refOut.fZ)<220");
  htemp->SetYTitle("#sqrt{#DeltaP_{rec}/P}");
  htemp->SetXTitle("#sqrt{#DeltaP_{mc}/P}");
  gPad->SaveAs("picLoss/dPcorr_TPC_step1.gif");
  gPad->SaveAs("picLoss/dPcorr_TPC_step1.eps");
  //
   // TPC
  chain->Draw("-sqrt(dPRec):-sqrt(dPMC)","abs(pdg)!=11&&type==3&&p.Pt()<0.5&&abs(p.R())<1&&abs(refOut.fZ)<220");
  htemp->SetYTitle("#sqrt{#DeltaP_{rec}/P}");
  htemp->SetXTitle("#sqrt{#DeltaP_{mc}/P}");
  gPad->SaveAs("picLoss/dPcorr_TPCseed_step1.gif");
  gPad->SaveAs("picLoss/dPcorr_TPCseed_step1.eps");


  // TRD
  chain->Draw("-sqrt(dPRec):-sqrt(dPMC)","abs(pdg)!=11&&type==2&&p.Pt()<0.5&&abs(p.R())<1&&abs(refOut.fZ)<220");
  htemp->SetYTitle("#sqrt{#DeltaP_{rec}/P}");
  htemp->SetXTitle("#sqrt{#DeltaP_{mc}/P}");
  gPad->SaveAs("picLoss/dPcorr_TRD_step1.gif");
  gPad->SaveAs("picLoss/dPcorr_TRD_step1.eps");

  //
  //
  //
  chain->Draw("(par.Pt()-refIn.Pt())/refIn.Pt()>>his(100,-0.02,0.02)","abs(pdg)!=11&&type==3&&p.Pt()<0.5&&abs(p.R())<1&&abs(refOut.fZ)<220");
  his->SetXTitle("(P_{trec}-P_{tmc})/P_{tmc}");
  gPad->SaveAs("picLoss/dPtcorr_TPCseed_step1_1D.eps");
  gPad->SaveAs("picLoss/dPtcorr_TPCseed_step1_1D.gif");

  chain->Draw("(par.P()-refIn.P())/refIn.P()>>his(100,-0.02,0.02)","abs(pdg)!=11&&type==3&&p.Pt()<0.5&&abs(p.R())<1&&abs(refOut.fZ)<220");
  his->SetXTitle("(P_{rec}-P_{mc})/P_{mc}");
  gPad->SaveAs("picLoss/dPcorr_TPCseed_step1_1D.eps");
  gPad->SaveAs("picLoss/dPcorr_TPCseed_step1_1D.gif");
*/
