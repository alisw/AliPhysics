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

#include "AliLog.h"
#include "AliAnalysisTask.h"

#include "AliTRDrecoTask.h"

ClassImp(AliTRDrecoTask)

TList* AliTRDrecoTask::fgTrendPoint(NULL);
TTreeSRedirector* AliTRDrecoTask::fgDebugStream(NULL);
//_______________________________________________________
AliTRDrecoTask::AliTRDrecoTask()
  : AliAnalysisTaskSE()
  ,fNRefFigures(0)
  ,fContainer(NULL)
  ,fTracks(NULL)
  ,fkTrack(NULL)
  ,fkMC(NULL)
  ,fkESD(NULL)
  ,fPlotFuncList(NULL)
{
// Default constructor  
}

//_______________________________________________________
AliTRDrecoTask::AliTRDrecoTask(const char *name, const char *title)
  : AliAnalysisTaskSE(name)
  ,fNRefFigures(0)
  ,fContainer(NULL)
  ,fTracks(NULL)
  ,fkTrack(NULL)
  ,fkMC(NULL)
  ,fkESD(NULL)
  ,fPlotFuncList(NULL)
{
// Constructor for all derived performance tasks

  SetTitle(title);
  DefineInput (1, TObjArray::Class()); // track list
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
  
  if(fContainer){
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
void AliTRDrecoTask::UserExec(Option_t *)
{
// Loop over Plot functors published by particular tasks

  fTracks = dynamic_cast<TObjArray *>(GetInputData(1));
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

    TMethodCall *plot = NULL;
    plotIter.Reset();
    while((plot=dynamic_cast<TMethodCall*>(plotIter()))){
      plot->Execute(this);
    }
  }
  PostData(1, fContainer);
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
  return kFALSE;
}

//_______________________________________________________
void AliTRDrecoTask::MakeSummary(){
// To be implemented by particular tasks
  AliWarning("Summary not available");
}

//_______________________________________________________
void AliTRDrecoTask::SetDebugLevel(Int_t level)
{
// Generic debug handler

  AliAnalysisTaskSE::SetDebugLevel(level);
  if(DebugLevel()>=1){
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
  if(HasPostProcess()) PostProcess();
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
