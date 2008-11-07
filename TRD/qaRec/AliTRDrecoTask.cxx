#include "TClass.h"
#include "TMethod.h"
#include "TMethodCall.h"
#include "TMethodArg.h"
#include "TFile.h"
#include "TList.h"
#include "TObjArray.h"
#include "TDirectory.h"
#include "TTreeStream.h"

#include "AliLog.h"
#include "AliAnalysisTask.h"

#include "AliTRDrecoTask.h"

ClassImp(AliTRDrecoTask)

//_______________________________________________________
AliTRDrecoTask::AliTRDrecoTask(const char *name, const char *title)
  : AliAnalysisTask(name, title)
  ,fNRefFigures(0)
  ,fDebugLevel(0)
  ,fPlotFuncList(0x0)
  ,fContainer(0x0)
  ,fTracks(0x0)
  ,fTrack(0x0)
  ,fMC(0x0)
  ,fESD(0x0)
  ,fDebugStream(0x0)
{
  DefineInput(0, TObjArray::Class());
  DefineOutput(0, TObjArray::Class());
}

//_______________________________________________________
AliTRDrecoTask::~AliTRDrecoTask() 
{
  if(fDebugStream){ 
    delete fDebugStream;
    fDebugStream = 0x0;
  }

  if(fPlotFuncList){
    fPlotFuncList->Delete();
    delete fPlotFuncList;
    fPlotFuncList = 0x0;
  }
  
  if(fContainer){
    //fContainer->Delete();
    delete fContainer;
    fContainer = 0x0;
  }
}

//_______________________________________________________
void AliTRDrecoTask::ConnectInputData(Option_t *)
{
  //
  // Connect input data
  //

  fTracks = dynamic_cast<TObjArray *>(GetInputData(0));
}

//_______________________________________________________
void AliTRDrecoTask::Exec(Option_t *)
{
  if(!fPlotFuncList){
    AliWarning("No functor list defined for the reference plots");
    return;
  }
  if(!fTracks) return;
  if(!fTracks->GetEntriesFast()) return;
  
  AliTRDtrackInfo *trackInfo = 0x0;
  TIter plotIter(fPlotFuncList);
  TObjArrayIter trackIter(fTracks);
  while((trackInfo = dynamic_cast<AliTRDtrackInfo*>(trackIter()))){
    fTrack = trackInfo->GetTrack();
    fMC    = trackInfo->GetMCinfo();
    fESD   = trackInfo->GetESDinfo();

    TMethodCall *plot = 0x0;
    plotIter.Reset();
    while((plot=dynamic_cast<TMethodCall*>(plotIter()))){
      plot->Execute(this);
    }
  }
  PostData(0, fContainer);
}

//_______________________________________________________
void AliTRDrecoTask::GetRefFigure(Int_t /*ifig*/)
{
  AliWarning("Retrieving reference figures not implemented.");
}

//_______________________________________________________
void AliTRDrecoTask::InitFunctorList()
{
  TClass *c = this->IsA();

  TMethod *m = 0x0;
  TIter methIter(c->GetListOfMethods());
  while((m=dynamic_cast<TMethod*>(methIter()))){
    TString name(m->GetName());
    if(!name.BeginsWith("Plot")) continue;
    if(!fPlotFuncList) fPlotFuncList = new TList();
    fPlotFuncList->AddLast(new TMethodCall(c, (const char*)name, ""));
  }
}

//_______________________________________________________
Bool_t AliTRDrecoTask::Load(const Char_t *filename)
{
  if(!TFile::Open(filename)){
    AliWarning(Form("Couldn't open file %s.", filename));
    return kFALSE;
  }
  TObjArray *o = 0x0;
  if(!(o = (TObjArray*)gFile->Get(GetName()))){
    AliWarning("Missing histogram container.");
    return kFALSE;
  }
  fContainer = (TObjArray*)o->Clone(GetName());
  gFile->Close();
  return kTRUE;
}

//_______________________________________________________
Bool_t AliTRDrecoTask::PostProcess()
{
  AliWarning("Post processing of reference histograms not implemented.");
  return kFALSE;
}

//_______________________________________________________
void AliTRDrecoTask::SetDebugLevel(Int_t level)
{
  fDebugLevel = level;
  if(fDebugLevel>=1){
    TDirectory *savedir = gDirectory;
    fDebugStream = new TTreeSRedirector(Form("TRD.Debug%s.root", GetName()));
    savedir->cd();
  }
}
