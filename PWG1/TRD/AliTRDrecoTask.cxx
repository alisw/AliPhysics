///////////////////////////////////////////////////////
//
// Basic class for Performance/Calibration TRD tasks
// 
// It performs generic tasks like :
//   - data file manegment
//   - reference container management
//   - debug container management
//   - interaction with AliAnalysisManager
//   - Plot functor loop
//
// Author: Alexandru Bercuci <A.Bercuci@gsi.de>, 10/09/2008
//
/////////////////////////////////////////////////////////

#include "TClass.h"
#include "TMethod.h"
#include "TMethodCall.h"
#include "TMethodArg.h"
#include "TFile.h"
#include "TChain.h"
#include "TList.h"
#include "TMap.h"
#include "TH1.h"
#include "TF1.h"
#include "TObjArray.h"
#include "TDirectory.h"
#include "TTreeStream.h"
#include "TBox.h"
#include "TVectorT.h"

#include "AliLog.h"
#include "AliAnalysisTask.h"
#include "AliExternalTrackParam.h"

#include "info/AliTRDeventInfo.h"
#include "AliTRDrecoTask.h"
#include "AliTRDtrackV1.h"
#include "AliTRDpidUtil.h"

ClassImp(AliTRDrecoTask)

TList* AliTRDrecoTask::fgTrendPoint(NULL);
TTreeSRedirector* AliTRDrecoTask::fgDebugStream(NULL);
//_______________________________________________________
AliTRDrecoTask::AliTRDrecoTask()
  : AliAnalysisTaskSE()
  ,fNRefFigures(0)
  ,fDets(NULL)
  ,fContainer(NULL)
  ,fEvent(NULL)
  ,fTracks(NULL)
  ,fkTrack(NULL)
  ,fkMC(NULL)
  ,fkESD(NULL)
  ,fSpecies(-6)
  ,fPt(-1.)
  ,fPhi(0.)
  ,fEta(0.)
  ,fPlotFuncList(NULL)
  ,fRunTerminate(kFALSE)
{
// Default constructor
  snprintf(fNameId, 10, "no name");
}

//_______________________________________________________
AliTRDrecoTask::AliTRDrecoTask(const char *name, const char *title)
  : AliAnalysisTaskSE(name)
  ,fNRefFigures(0)
  ,fDets(NULL)
  ,fContainer(NULL)
  ,fEvent(NULL)
  ,fTracks(NULL)
  ,fkTrack(NULL)
  ,fkMC(NULL)
  ,fkESD(NULL)
  ,fSpecies(-6)
  ,fPt(-1.)
  ,fPhi(0.)
  ,fEta(0.)
  ,fPlotFuncList(NULL)
  ,fRunTerminate(kFALSE)
{
// Constructor for all derived performance tasks

  SetTitle(title);
  snprintf(fNameId, 10, "no name");
  DefineInput (1, TObjArray::Class()); // track list
  DefineInput (2, AliTRDeventInfo::Class()); // event info object
  DefineOutput(1, TObjArray::Class()); // histogram list
}

//_______________________________________________________
AliTRDrecoTask::~AliTRDrecoTask() 
{

  // Generic task destructor

  AliDebug(2, Form(" Ending task %s[%s]", GetName(), GetTitle()));
  if(fgDebugStream){ 
    delete fgDebugStream;
    fgDebugStream = NULL;
  }

  if(fPlotFuncList){
    fPlotFuncList->Delete();
    delete fPlotFuncList;
    fPlotFuncList = NULL;
  }
  
  if(fDets){
    if(fDets->IsOwner()) fDets->Delete();
    delete fDets;
    fDets = NULL;
  }

  if(fContainer && !(AliAnalysisManager::GetAnalysisManager() && AliAnalysisManager::GetAnalysisManager()->IsProofMode())){
    if(fContainer->IsOwner()) fContainer->Delete();
    delete fContainer;
    fContainer = NULL;
  }

  if(fgTrendPoint){
    TFile::Open("TRD.PerformanceTrend.root", "UPDATE");
    fgTrendPoint->Write();
    delete fgTrendPoint;
    fgTrendPoint=NULL;
    gFile->Close();
  }
}

//_______________________________________________________
Int_t AliTRDrecoTask::GetNRefFigures() const  
{ 
  if(!fNRefFigures) AliWarning("No reference plots available.");
  return fNRefFigures; 
} 

//_______________________________________________________
void AliTRDrecoTask::UserCreateOutputObjects()
{
  if(!HasFunctorList()) InitFunctorList();
  fContainer = Histos();
  PostData(1, fContainer);
}

//_______________________________________________________
void AliTRDrecoTask::UserExec(Option_t *)
{
// Loop over Plot functors published by particular tasks

  fTracks = dynamic_cast<TObjArray *>(GetInputData(1));
  fEvent  = dynamic_cast<AliTRDeventInfo *>(GetInputData(2));

  if(!fPlotFuncList){
    AliWarning("No functor list defined for the reference plots");
    return;
  }
  if(!fTracks) return;
  if(!fTracks->GetEntriesFast()) return;
  else AliDebug(2, Form("Tracks[%d] for %s", fTracks->GetEntriesFast(), GetName()));

  AliTRDtrackInfo *trackInfo = NULL;
  TIter plotIter(fPlotFuncList);
  TObjArrayIter trackIter(fTracks);
  while((trackInfo = dynamic_cast<AliTRDtrackInfo*>(trackIter()))){
    fkTrack = trackInfo->GetTrack();
    fkMC    = trackInfo->GetMCinfo();
    fkESD   = trackInfo->GetESDinfo();
    // cache properties of the track at TRD entrance
    // check input track status
    AliExternalTrackParam *tin(NULL);
    fPt=-1; fEta=0.; fPhi=0.; fSpecies=-6;
    if(!fkTrack || !(tin = fkTrack->GetTrackIn())) AliDebug(2, "Track did not entered TRD fiducial volume.");
    else {
      fPt   = tin->Pt();
      fEta  = tin->Eta();
      Double_t xyz[3];
      if(!tin->GetXYZ(xyz)) AliDebug(2, "Failed getting global track postion");
      else fPhi  = TMath::ATan2(xyz[1], xyz[0]);
      fSpecies= fkTrack->Charge()*(AliTRDpidUtil::Mass2Pid(fkTrack->GetMass())+1);
    }

    TMethodCall *plot = NULL;
    plotIter.Reset();
    while((plot=dynamic_cast<TMethodCall*>(plotIter()))){
      plot->Execute(this);
    }
  }
}

//_______________________________________________________
Bool_t AliTRDrecoTask::GetRefFigure(Int_t /*ifig*/)
{
  AliWarning("Retrieving reference figures not implemented.");
  return kFALSE;
}

//_______________________________________________________
Bool_t AliTRDrecoTask::PutTrendValue(const Char_t *name, Double_t val)
{
// Generic publisher for trend values

  if(!fgTrendPoint){
    fgTrendPoint = new TList();
    fgTrendPoint->SetOwner();
  }
  fgTrendPoint->AddLast(new TNamed(Form("%s_%s", GetName(), name), Form("%f", val)));
  return kTRUE;
}

//_______________________________________________________
void AliTRDrecoTask::InitFunctorList()
{
// Initialize list of functors

  TClass *c = this->IsA();
  if(fPlotFuncList) fPlotFuncList->Clear();

  TMethod *m = NULL;
  TIter methIter(c->GetListOfMethods());
  while((m=dynamic_cast<TMethod*>(methIter()))){
    TString name(m->GetName());
    if(!name.BeginsWith("Plot")) continue;
    if(!fPlotFuncList) fPlotFuncList = new TList();
    fPlotFuncList->AddLast(new TMethodCall(c, (const char*)name, ""));
  }
}

//_______________________________________________________
Bool_t AliTRDrecoTask::Load(const Char_t *file, const Char_t *dir)
{
// Generic container loader

  if(!TFile::Open(file)){
    AliWarning(Form("Couldn't open file %s.", file));
    return kFALSE;
  }
  if(!gFile->cd(dir)){
    AliWarning(Form("Couldn't cd to %s in %s.", dir, file));
    return kFALSE;
  }
  TObjArray *o = NULL;
  if(!(o = (TObjArray*)gDirectory->Get(GetName()))){
    AliWarning("Missing histogram container.");
    return kFALSE;
  }
  fContainer = (TObjArray*)o->Clone(GetName());
  gFile->Close();
  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDrecoTask::LoadDetectorMap(const Char_t *file, const Char_t *dir)
{
// Load detector map.

  if(!TFile::Open(file)){
    AliWarning(Form("Couldn't open file %s.", file));
    return kFALSE;
  }
  if(!gFile->cd(dir)){
    AliWarning(Form("Couldn't cd to %s in %s.", dir, file));
    return kFALSE;
  }
  TObjArray *info = NULL;
  if(!(info = (TObjArray*)gDirectory->Get("TRDinfoGen"))){
    AliWarning("Missing TRDinfoGen container.");
    return kFALSE;
  }
  TObjArray *dets = (TObjArray*)info->FindObject("Chambers");
  if(!dets){
    AliWarning("Missing detector map from TRDinfoGen results.");
    info->ls();
    return kFALSE;
  }
  fDets = (TObjArray*)dets->Clone("Chambers");
  gFile->Close();
  return kTRUE;
}


//________________________________________________________
Bool_t AliTRDrecoTask::Save(TObjArray * const results){
  //
  // Store the output graphs in a ROOT file
  // Input TObject array will not be written as Key to the file,
  // only content itself
  //

  TDirectory *cwd = gDirectory;
  if(!TFile::Open(Form("TRD.Result%s.root", GetName()), "RECREATE")) return kFALSE;

  TIterator *iter = results->MakeIterator();
  TObject *inObject = NULL, *outObject = NULL;
  while((inObject = iter->Next())){
    outObject = inObject->Clone();
    outObject->Write(NULL, TObject::kSingleKey);
  }
  delete iter;
  gFile->Close(); delete gFile;
  cwd->cd(); 
  return kTRUE;
}

//_______________________________________________________
Bool_t AliTRDrecoTask::PostProcess()
{
// To be implemented by particular tasks

  AliWarning("Post processing of reference histograms not implemented.");
  return kTRUE;
}

//_______________________________________________________
void AliTRDrecoTask::MakeDetectorPlot(Int_t ly)
{
// Draw chamber boundaries in eta/phi plots with misalignments
// based on info collected by AliTRDinfoGen

  if(!fDets){
    AliWarning("Detector map and status not available.");
    return;
  }

  TBox *gdet = new TBox();
  gdet->SetLineWidth(kBlack);gdet->SetFillColor(kBlack);
  Int_t style[] = {0, 3003};
  for(Int_t idet(0); idet<540; idet++){
    if(idet%6 != ly) continue;
    TVectorF *det((TVectorF*)fDets->At(idet));
    if(!det) continue;
    AliDebug(2, Form("det[%03d] 0[%+4.1f %+4.1f] 1[%+4.1f %+4.1f]", idet, (*det)[0], (*det)[1], (*det)[2], (*det)[3]));
    Int_t iopt = Int_t((*det)[4]);
    gdet->SetFillStyle(style[iopt]);
    gdet->DrawBox((*det)[0], (*det)[1], (*det)[2], (*det)[3]);
  }
}


//_______________________________________________________
void AliTRDrecoTask::MakeSummary()
{
// To be implemented by particular tasks
  AliWarning("Summary not available");
}

//_______________________________________________________
void AliTRDrecoTask::SetDebugLevel(Int_t level)
{
// Generic debug handler

  AliAnalysisTaskSE::SetDebugLevel(level);
  if(DebugLevel()>=1){
    AliInfo(Form("Debug Level for Task %s set to %d", GetName(), level));
    TDirectory *savedir = gDirectory;
    fgDebugStream = new TTreeSRedirector("TRD.DebugPerformance.root");
    savedir->cd();
  }
}

//____________________________________________________________________
void AliTRDrecoTask::Terminate(Option_t *)
{
  //
  // Terminate
  //

  if(fgDebugStream){ 
    delete fgDebugStream;
    fgDebugStream = NULL;
  }
  fContainer = dynamic_cast<TObjArray *>(GetOutputData(1));
  if(fContainer && fRunTerminate){
    PostProcess();
    MakeSummary();
  }
}

//________________________________________________________
void AliTRDrecoTask::Adjust(TF1 *f, TH1 * const h)
{
// Helper function to avoid duplication of code
// Make first guesses on the fit parameters

  // find the intial parameters of the fit !! (thanks George)
  Int_t nbinsy = Int_t(.5*h->GetNbinsX());
  Double_t sum = 0.;
  for(Int_t jbin=nbinsy-4; jbin<=nbinsy+4; jbin++) sum+=h->GetBinContent(jbin); sum/=9.;
  f->SetParLimits(0, 0., 3.*sum);
  f->SetParameter(0, .9*sum);

  f->SetParLimits(1, -.2, .2);
  f->SetParameter(1, -0.1);

  f->SetParLimits(2, 0., 4.e-1);
  f->SetParameter(2, 2.e-2);
  if(f->GetNpar() <= 4) return;

  f->SetParLimits(3, 0., sum);
  f->SetParameter(3, .1*sum);

  f->SetParLimits(4, -.3, .3);
  f->SetParameter(4, 0.);

  f->SetParLimits(5, 0., 1.e2);
  f->SetParameter(5, 2.e-1);
}
