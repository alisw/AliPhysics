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
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TObjArray.h"
#include "TDirectory.h"
#include "TTreeStream.h"
#include "TBox.h"
#include "TLatex.h"
#include "TVectorT.h"

#include "AliLog.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliExternalTrackParam.h"

#include "info/AliTRDeventInfo.h"
#include "AliTRDrecoTask.h"
#include "AliTRDtrackV1.h"
#include "AliTRDpidUtil.h"

ClassImp(AliTRDrecoTask)

Float_t AliTRDrecoTask::fgPt0[AliTRDrecoTask::fgNPt0] = {0.5, 0.8, 1.5, 5};
TList* AliTRDrecoTask::fgTrendPoint(NULL);
TTreeSRedirector* AliTRDrecoTask::fgDebugStream(NULL);
//_______________________________________________________
AliTRDrecoTask::AliTRDrecoTask()
  : AliAnalysisTaskSE()
  ,fNRefFigures(0)
  ,fDets(NULL)
  ,fDetsV(NULL)
  ,fContainer(NULL)
  ,fEvent(NULL)
  ,fTracks(NULL)
  ,fClusters(NULL)
  ,fkClusters(NULL)
  ,fkTrack(NULL)
  ,fkMC(NULL)
  ,fkESD(NULL)
  ,fSpecies(-6)
  ,fPt(-1.)
  ,fPhi(0.)
  ,fEta(0.)
  ,fPlotFuncList(NULL)
  ,fDetFuncList(NULL)
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
  ,fDetsV(NULL)
  ,fContainer(NULL)
  ,fEvent(NULL)
  ,fTracks(NULL)
  ,fClusters(NULL)
  ,fkClusters(NULL)
  ,fkTrack(NULL)
  ,fkMC(NULL)
  ,fkESD(NULL)
  ,fSpecies(-6)
  ,fPt(-1.)
  ,fPhi(0.)
  ,fEta(0.)
  ,fPlotFuncList(NULL)
  ,fDetFuncList(NULL)
  ,fRunTerminate(kFALSE)
{
// Constructor for all derived performance tasks

  SetTitle(title);
  snprintf(fNameId, 10, "no name");
  DefineInput (1, TObjArray::Class()); // track list
  DefineInput (2, AliTRDeventInfo::Class()); // event info object
  DefineInput (3, TObjArray::Class()); // cluster list object
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
  if(fDetFuncList){
    fDetFuncList->Delete();
    delete fDetFuncList;
    fDetFuncList = NULL;
  }
  
  if(fDets){
    if(fDets->IsOwner()) fDets->Delete();
    delete fDets;
    fDets = NULL;
  }
  if(fDetsV) delete fDetsV; fDetsV=NULL;

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

//____________________________________________________________________
Int_t AliTRDrecoTask::GetPtBinSignificant(Float_t pt)
{
// Get significant (very low, low, medium, high, very high) pt bin

  Int_t ipt(0);
  while(ipt<fgNPt0){
    if(pt<fgPt0[ipt]) break;
    ipt++;
  }
  return ipt-1;
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

  fTracks   = dynamic_cast<TObjArray *>(GetInputData(1));
  fEvent    = dynamic_cast<AliTRDeventInfo *>(GetInputData(2));
  fClusters = dynamic_cast<TObjArray*>(GetInputData(3));

  if(!fPlotFuncList){
    AliWarning("No track functor list defined for the task");
    return;
  }
  if(!fTracks) return;
  if(!fTracks->GetEntriesFast()) return;
  else AliDebug(2, Form("Tracks[%d] for %s", fTracks->GetEntriesFast(), GetName()));

  Int_t itrk(-1);
  AliTRDtrackInfo *trackInfo(NULL);
  TIter plotIter(fPlotFuncList);
  TObjArrayIter trackIter(fTracks);
  while((trackInfo = dynamic_cast<AliTRDtrackInfo*>(trackIter()))){
    itrk++; fPt=-1; fEta=0.; fPhi=0.; fSpecies=-6;
    fkMC    = trackInfo->GetMCinfo();
    fkESD   = trackInfo->GetESDinfo();
    if((fkTrack = trackInfo->GetTrack())){
      // cache properties of the track at TRD entrance
      // check input track status
      AliExternalTrackParam *tin(NULL);
      if(!(tin = fkTrack->GetTrackIn())) AliDebug(2, Form("Missing TRD track[%d] :: entry point.", itrk));
      else {
        fPt   = tin->Pt();
        fEta  = tin->Eta();
        Double_t xyz[3];
        if(!tin->GetXYZ(xyz)) AliDebug(2, Form("Failed TRD track[%d] :: global track postion", itrk));
        else fPhi  = TMath::ATan2(xyz[1], xyz[0]);
        fSpecies= fkTrack->Charge()*(AliTRDpidUtil::Mass2Pid(fkTrack->GetMass())+1);
      }
    } else AliDebug(2, Form("Missing TRD track[%d].", itrk));

    TMethodCall *plot(NULL);
    plotIter.Reset();
    while((plot=dynamic_cast<TMethodCall*>(plotIter()))) plot->Execute(this);
  }
  if(!fClusters) return;
  if(!fDetFuncList){
    AliDebug(1, "No detector functor list defined for task");
    return;
  }
  TIter detIter(fDetFuncList);
  for(Int_t idet(0); idet<AliTRDgeometry::kNdet; idet++){
    if(!(fkClusters = (TObjArray*)fClusters->At(idet))) continue;
    TMethodCall *det(NULL);
    detIter.Reset();
    while((det=dynamic_cast<TMethodCall*>(detIter()))) det->Execute(this);
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
  if(fDetFuncList) fDetFuncList->Clear();

  TMethod *m(NULL);
  TIter methIter(c->GetListOfMethods());
  while((m=dynamic_cast<TMethod*>(methIter()))){
    TString name(m->GetName());
    if(name.BeginsWith("Plot")){
      if(!fPlotFuncList) fPlotFuncList = new TList();
      fPlotFuncList->AddLast(new TMethodCall(c, (const char*)name, ""));
    } else if(name.BeginsWith("Det")){
      if(!fDetFuncList) fDetFuncList = new TList();
      fDetFuncList->AddLast(new TMethodCall(c, (const char*)name, ""));
    }
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
    TVector *vdets = (TVector*)info->At(4);
    if(!vdets){
      AliWarning("Missing detector map from TRDinfoGen results.");
      return kFALSE;
    } else fDetsV = (TVector*)vdets->Clone();
  } else fDets = (TObjArray*)dets->Clone("Chambers");
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
void AliTRDrecoTask::MakeDetectorPlot(Int_t ly, const Option_t *opt)
{
// Draw chamber boundaries in eta/phi plots with misalignments
// based on info collected by AliTRDinfoGen

  if(!fDets){
    AliWarning("OLD Detector map and status not available. Try NEW");
    MakeDetectorPlotNEW(ly, opt);
    return;
  }
  Float_t xmin(0.), xmax(0.);
  TBox *gdet = new TBox();
  gdet->SetLineColor(kBlack);gdet->SetFillColor(kBlack);
  Int_t style[] = {0, 3003};
  for(Int_t idet(0); idet<540; idet++){
    if(idet%6 != ly) continue;
    TVectorF *det((TVectorF*)fDets->At(idet));
    if(!det) continue;
    Int_t iopt = Int_t((*det)[4]);
    if(strcmp(opt, "eta")==0){
      xmin=(*det)[0]; xmax=(*det)[2];
    } else if(strcmp(opt, "pad")==0){
      Int_t stk(AliTRDgeometry::GetStack(idet));
      xmin=-0.6+16*(4-stk)-(stk<2?4:0); xmax=xmin+(stk==2?12:16)-0.2;
    } else continue;
    AliDebug(2, Form("det[%03d] 0[%+4.1f(%+4.1f) %+4.1f] 1[%+4.1f(%+4.1f) %+4.1f] opt[%d]", idet, xmin, (*det)[0], (*det)[1], xmax, (*det)[2], (*det)[3], iopt));
    if(iopt==1){
      gdet->SetFillStyle(style[1]);gdet->SetFillColor(kBlack);
      gdet->DrawBox(xmin, (*det)[1], xmax, (*det)[3]);
    } else {
      gdet->SetFillStyle(style[0]);
      gdet->DrawBox(xmin, (*det)[1], xmax, (*det)[3]);
      if(iopt==2){
        gdet->SetFillStyle(style[1]);gdet->SetFillColor(kGreen);
        gdet->DrawBox(xmin, (*det)[1], xmax, 0.5*((*det)[3]+(*det)[1]));
      } else if(iopt==3){
        gdet->SetFillStyle(style[1]);gdet->SetFillColor(kRed);
        gdet->DrawBox(xmin, 0.5*((*det)[3]+(*det)[1]), xmax, (*det)[3]);
      } else if(iopt!=0) AliError(Form("Wrong chmb. status[%d] for det[%03d]", iopt, idet));
    }
  }
  Float_t dsm = TMath::TwoPi()/AliTRDgeometry::kNsector;
  xmin=0.;
  if(strcmp(opt, "pad")==0) xmin=38.;
  TLatex *sm = new TLatex(); sm->SetTextAlign(22);sm->SetTextColor(0); sm->SetTextFont(32);sm->SetTextSize(0.03);
  for(Int_t is(0); is<AliTRDgeometry::kNsector; is++) sm->DrawLatex(xmin, -TMath::Pi()+(is+0.5)*dsm, Form("%02d", is>=9?(is-9):(is+9)));
}

//_______________________________________________________
void AliTRDrecoTask::MakeDetectorPlotNEW(Int_t ly, const Option_t *opt)
{
// Draw chamber boundaries in eta/phi plots with misalignments
// based on info collected by AliTRDinfoGen NEW data storage

  if(!fDetsV){
    AliWarning("NEW Detector map and status not available.");
    return;
  }
  Float_t xmin(0.), xmax(0.);
  TBox *gdet = new TBox();
  gdet->SetLineColor(kBlack);gdet->SetFillColor(kBlack);
  Int_t style[] = {0, 3003};
  for(Int_t idet(0), jdet(0); idet<AliTRDgeometry::kNdet; idet++, jdet+=5){
    if(idet%6 != ly) continue;
    Int_t iopt = Int_t((*fDetsV)[jdet+4]);
    if(strcmp(opt, "eta")==0){
      xmin=(*fDetsV)[jdet+0]; xmax=(*fDetsV)[jdet+2];
    } else if(strcmp(opt, "pad")==0){
      Int_t stk(AliTRDgeometry::GetStack(idet));
      xmin=-0.6+16*(4-stk)-(stk<2?4:0); xmax=xmin+(stk==2?12:16)-0.2;
    } else continue;
    AliDebug(2, Form("det[%03d] 0[%+4.1f(%+4.1f) %+4.1f] 1[%+4.1f(%+4.1f) %+4.1f] opt[%d]", idet, xmin, (*fDetsV)[jdet+0], (*fDetsV)[jdet+1], xmax, (*fDetsV)[jdet+2], (*fDetsV)[jdet+3], iopt));
    if(iopt==1){
      gdet->SetFillStyle(style[1]);gdet->SetFillColor(kBlack);
      gdet->DrawBox(xmin, (*fDetsV)[jdet+1], xmax, (*fDetsV)[jdet+3]);
    } else {
      gdet->SetFillStyle(style[0]);
      gdet->DrawBox(xmin, (*fDetsV)[jdet+1], xmax, (*fDetsV)[jdet+3]);
      if(iopt==2){
        gdet->SetFillStyle(style[1]);gdet->SetFillColor(kGreen);
        gdet->DrawBox(xmin, (*fDetsV)[jdet+1], xmax, 0.5*((*fDetsV)[jdet+3]+(*fDetsV)[jdet+1]));
      } else if(iopt==3){
        gdet->SetFillStyle(style[1]);gdet->SetFillColor(kRed);
        gdet->DrawBox(xmin, 0.5*((*fDetsV)[jdet+3]+(*fDetsV)[jdet+1]), xmax, (*fDetsV)[jdet+3]);
      } else if(iopt!=0) AliError(Form("Wrong chmb. status[%d] for det[%03d]", iopt, idet));
    }
  }
  Float_t dsm = TMath::TwoPi()/AliTRDgeometry::kNsector;
  xmin=0.;
  if(strcmp(opt, "pad")==0) xmin=38.;
  TLatex *sm = new TLatex(); sm->SetTextAlign(22);sm->SetTextColor(0); sm->SetTextFont(32);sm->SetTextSize(0.03);
  for(Int_t is(0); is<AliTRDgeometry::kNsector; is++) sm->DrawLatex(xmin, -TMath::Pi()+(is+0.5)*dsm, Form("%02d", is>=9?(is-9):(is+9)));
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


//________________________________________________________
void AliTRDrecoTask::SetNormZ(TH2 *h2, Int_t bxmin, Int_t bxmax, Int_t bymin, Int_t bymax, Float_t thr)
{
// Normalize histo content to the mean value in the range specified by bin ranges
// [bxmin, bxmax] on the x axis and [bymin, bymax] on the y axis.
// Optionally a threshold "thr" can be specified to disregard entries with no meaning

  Float_t s = 0., c=0.; Int_t is(0);
  for(Int_t ix(bxmin); ix<=(bxmax>0?bxmax:(h2->GetXaxis()->GetNbins())); ix++){
    for(Int_t iy(bymin); iy<=(bymax>0?bymax:(h2->GetYaxis()->GetNbins())); iy++){
      if((c = h2->GetBinContent(ix, iy))<thr) continue;
      s += c; is++;
    }
  }
  s/= (is?is:1);
  for(Int_t ix(1); ix<=h2->GetXaxis()->GetNbins(); ix++){
    for(Int_t iy(1); iy<=h2->GetYaxis()->GetNbins(); iy++){
      if((c = h2->GetBinContent(ix, iy))<thr) h2->SetBinContent(ix, iy, thr-1000);
      else h2->SetBinContent(ix, iy, 100.*(c/s-1.));
    }
  }
}

//________________________________________________________
void AliTRDrecoTask::SetRangeZ(TH2 *h2, Float_t min, Float_t max, Float_t thr)
{
// Set range on Z axis such to avoid outliers

  Float_t c(0.), dz(1.e-3*(max-min));
  for(Int_t ix(1); ix<=h2->GetXaxis()->GetNbins(); ix++){
    for(Int_t iy(1); iy<=h2->GetYaxis()->GetNbins(); iy++){
      if((c = h2->GetBinContent(ix, iy))<thr) continue;
      if(c<=min) h2->SetBinContent(ix, iy, min+dz);
    }
  }
  h2->GetZaxis()->SetRangeUser(min, max);
}

//________________________________________________________
Float_t AliTRDrecoTask::GetMeanStat(TH1 *h, Float_t cut, Option_t *opt)
{
// return mean number of entries/bin of histogram "h"
// if option "opt" is given the following values are accepted:
//   "<" : consider only entries less than "cut"
//   ">" : consider only entries greater than "cut"

  //Int_t dim(h->GetDimension());
  Int_t nbx(h->GetNbinsX()), nby(h->GetNbinsY()), nbz(h->GetNbinsZ());
  Double_t sum(0.); Int_t n(0);
  for(Int_t ix(1); ix<=nbx; ix++)
    for(Int_t iy(1); iy<=nby; iy++)
      for(Int_t iz(1); iz<=nbz; iz++){
        if(strcmp(opt, "")==0){sum += h->GetBinContent(ix, iy, iz); n++;}
        else{
          if(strcmp(opt, "<")==0) {
            if(h->GetBinContent(ix, iy, iz)<cut) {sum += h->GetBinContent(ix, iy, iz); n++;}
          } else if(strcmp(opt, ">")==0){
            if(h->GetBinContent(ix, iy, iz)>cut) {sum += h->GetBinContent(ix, iy, iz); n++;}
          } else {sum += h->GetBinContent(ix, iy, iz); n++;}
        }
      }
  return n>0?sum/n:0.;
}

//________________________________________________________
AliTRDrecoTask::AliTRDrecoProjection::AliTRDrecoProjection()
  :TNamed()
  ,fH(NULL)
  ,fNrebin(0)
  ,fRebinX(NULL)
  ,fRebinY(NULL)
{
  // constructor
  memset(fAx, 0, 3*sizeof(Int_t));
  memset(fRange, 0, 4*sizeof(Float_t));
}

//________________________________________________________
AliTRDrecoTask::AliTRDrecoProjection::~AliTRDrecoProjection()
{
  // destructor
  if(fH) delete fH;
}

//________________________________________________________
void AliTRDrecoTask::AliTRDrecoProjection::Build(const Char_t *n, const Char_t *t, Int_t ix, Int_t iy, Int_t iz, TAxis *aa[])
{
// check and build (if neccessary) projection determined by axis "ix", "iy" and "iz"
  if(!aa[ix] || !aa[iy] || !aa[iz]) return;
  TAxis *ax(aa[ix]), *ay(aa[iy]), *az(aa[iz]);
  // check ax definiton to protect against older versions of the data
  if(ax->GetNbins()<=0 || (ax->GetXmax()-ax->GetXmin())<=0.){
    AliWarning(Form("Wrong definition of axis[%d] \"%s\"[%d](%f %f).", ix, ax->GetTitle(), ax->GetNbins(), ax->GetXmin(), ax->GetXmax()));
    return;
  }
  if(ay->GetNbins()<=0 || (ay->GetXmax()-ay->GetXmin())<=0.){
    AliWarning(Form("Wrong definition of axis[%d] \"%s\"[%d](%f %f).", ix, ay->GetTitle(), ay->GetNbins(), ay->GetXmin(), ay->GetXmax()));
    return;
  }
  if(az->GetNbins()<=0 || (az->GetXmax()-az->GetXmin())<=0.){
    AliWarning(Form("Wrong definition of axis[%d] \"%s\"[%d](%f %f).", ix, az->GetTitle(), az->GetNbins(), az->GetXmin(), az->GetXmax()));
    return;
  }
  SetNameTitle(n,t);
  fH = new TH3I(n, Form("%s;%s;%s;%s", t, ax->GetTitle(), ay->GetTitle(), az->GetTitle()),
    ax->GetNbins(), ax->GetXmin(), ax->GetXmax(),
    ay->GetNbins(), ay->GetXmin(), ay->GetXmax(),
    az->GetNbins(), az->GetXmin(), az->GetXmax());
  fAx[0] = ix; fAx[1] = iy; fAx[2] = iz;
  fRange[0] = az->GetXmin()/3.; fRange[1] = az->GetXmax()/3.;
  AliDebug(2, Form("H3(%s, %s) :: %s[%3d %4.2f %4.2f]%s[%3d %4.2f %4.2f]%s[%3d %4.2f %4.2f]", n, t,
    ax->GetTitle(), ax->GetNbins(), ax->GetXmin(), ax->GetXmax(),
    ay->GetTitle(), ay->GetNbins(), ay->GetXmin(), ay->GetXmax(),
    az->GetTitle(), az->GetNbins(), az->GetXmin(), az->GetXmax()));
}

//________________________________________________________
AliTRDrecoTask::AliTRDrecoProjection& AliTRDrecoTask::AliTRDrecoProjection::operator=(const AliTRDrecoProjection& rhs)
{
// copy projections
  if(this == &rhs) return *this;

  TNamed::operator=(rhs);
  if(fNrebin){fNrebin=0; delete [] fRebinX; delete [] fRebinY;}
  if(rhs.fNrebin) SetRebinStrategy(rhs.fNrebin, rhs.fRebinX, rhs.fRebinY);
  memcpy(fAx, rhs.fAx, 3*sizeof(Int_t));
  memcpy(fRange, rhs.fRange, 4*sizeof(Float_t));
  if(fH) delete fH;
  if(rhs.fH) fH=(TH3I*)rhs.fH->Clone(Form("%s_CLONE", rhs.fH->GetName()));
  return *this;
}

//________________________________________________________
AliTRDrecoTask::AliTRDrecoProjection& AliTRDrecoTask::AliTRDrecoProjection::operator+=(const AliTRDrecoProjection& other)
{
// increment projections
  if(!fH || !other.fH) return *this;
  AliDebug(2, Form("%s+=%s [%s+=%s]", GetName(), other.GetName(), fH->GetName(), (other.fH)->GetName()));
  fH->Add(other.fH);
  return *this;
}

//________________________________________________________
void AliTRDrecoTask::AliTRDrecoProjection::Increment(Int_t bin[], Double_t v)
{
// increment bin with value "v" pointed by general coord in "bin"
  if(!fH) return;
  AliDebug(4, Form("  %s[%2d]", fH->GetName(), Int_t(v)));
  fH->AddBinContent(fH->GetBin(bin[fAx[0]],bin[fAx[1]],bin[fAx[2]]), Int_t(v));
}

//________________________________________________________
TH2* AliTRDrecoTask::AliTRDrecoProjection::Projection2D(const Int_t nstat, const Int_t ncol, const Int_t mid, Bool_t del)
{
// build the 2D projection and adjust binning

  const Char_t *title[] = {"Mean", "#mu", "MPV"};
  if(!fH){
    AliDebug(1, Form("Missing 3D in %s", GetName()));
    return NULL;
  }
  TAxis *ax(fH->GetXaxis()), *ay(fH->GetYaxis()), *az(fH->GetZaxis());
  TH2D *h2s(NULL), *hyx(NULL);
  if(!(h2s = (TH2D*)fH->Project3D("yx"))){
    AliDebug(1, Form("Failed Project3D(\"yx\") in %s", GetName()));
    return NULL;
  }
  // save a copy of the original distribution
  if(!del){
    hyx = (TH2D*)h2s->Clone();
    hyx->SetName(Form("%sEn", fH->GetName()));
  }
  Int_t irebin(0), dxBin(1), dyBin(1);
  while(irebin<fNrebin && (AliTRDrecoTask::GetMeanStat(h2s, .5, ">")<nstat)){
    h2s->Rebin2D(fRebinX[irebin], fRebinY[irebin]);
    dxBin*=fRebinX[irebin];dyBin*=fRebinY[irebin];
    irebin++;
  }
  Int_t nx(h2s->GetNbinsX()), ny(h2s->GetNbinsY());
  delete h2s;
  if(mid<0) return NULL;

  // start projection
  TH1 *h(NULL); Int_t n(0);
  Float_t dz=(fRange[1]-fRange[1])/ncol;
  TString titlez(az->GetTitle()); TObjArray *tokenTitle(titlez.Tokenize(" "));
  Int_t nt(tokenTitle->GetEntriesFast());
  TH2 *h2(NULL);
  if((h2 = (TH2*)gDirectory->Get(Form("%s_2D", fH->GetName())))) delete h2; // avoid ROOT warning messages
  h2 = new TH2F(Form("%s_2D", fH->GetName()),
            Form("%s;%s;%s;%s(%s) %s", fH->GetTitle(), ax->GetTitle(), ay->GetTitle(), title[mid], nt>0?(*tokenTitle)[0]->GetName():"", nt>1?(*tokenTitle)[1]->GetName():""),
            nx, ax->GetXmin(), ax->GetXmax(), ny, ay->GetXmin(), ay->GetXmax());
  h2->SetContour(ncol);
  h2->GetZaxis()->CenterTitle();
  h2->GetZaxis()->SetTitleOffset(1.4);
  h2->GetZaxis()->SetRangeUser(fRange[0], fRange[1]);
  AliDebug(2, Form("%s[%s] nx[%d] ny[%d]", h2->GetName(), h2->GetTitle(), nx, ny));
  for(Int_t iy(0); iy<ny; iy++){
    for(Int_t ix(0); ix<nx; ix++){
      h = fH->ProjectionZ(Form("%s_z", h2->GetName()), ix*dxBin+1, (ix+1)*dxBin, iy*dyBin+1, (iy+1)*dyBin);
      Int_t ne((Int_t)h->Integral());
      //printf("  x[%2d %2d] y[%2d %2d] ne[%4d]\n", ix*dxBin+1, (ix+1)*dxBin, iy*dyBin+1, (iy+1)*dyBin, ne);
      if(ne<nstat/2){
        h2->SetBinContent(ix+1, iy+1, -999);
        h2->SetBinError(ix+1, iy+1, 1.);
        n++;
      }else{
        // redo the projection by adding 1 bin @ left and 1 bin @ right for smoothing
        h = fH->ProjectionZ(Form("%s_z", h2->GetName()), ix*dxBin, (ix+1)*dxBin+1, iy*dyBin, (iy+1)*dyBin+1);
        Float_t v(h->GetMean()), ve(h->GetRMS());
        if(mid==1){
          TF1 fg("fg", "gaus", az->GetXmin(), az->GetXmax());
          fg.SetParameter(0, Float_t(ne)); fg.SetParameter(1, v); fg.SetParameter(2, ve);
          h->Fit(&fg, "WQ0");
          v = fg.GetParameter(1); ve = fg.GetParameter(2);
        } else if (mid==2) {
          TF1 fl("fl", "landau", az->GetXmin(), az->GetXmax());
          fl.SetParameter(0, Float_t(ne)); fl.SetParameter(1, v); fl.SetParameter(2, ve);
          h->Fit(&fl, "WQ0");
          v = fl.GetMaximumX(); ve = fl.GetParameter(2);
/*          TF1 fgle("gle", "[0]*TMath::Landau(x, [1], [2], 1)*TMath::Exp(-[3]*x/[1])", az->GetXmin(), az->GetXmax());
          fgle.SetParameter(0, fl.GetParameter(0));
          fgle.SetParameter(1, fl.GetParameter(1));
          fgle.SetParameter(2, fl.GetParameter(2));
          fgle.SetParameter(3, 1.);fgle.SetParLimits(3, 0., 5.);
          h->Fit(&fgle, "WQ");
          v = fgle.GetMaximumX(); ve = fgle.GetParameter(2);*/
        }
        if(v<fRange[0]) h2->SetBinContent(ix+1, iy+1, fRange[0]+0.1*dz);
        else h2->SetBinContent(ix+1, iy+1, v);
        h2->SetBinError(ix+1, iy+1, ve);
      }
    }
  }
  if(h) delete h;
  if(n==nx*ny){  // clean empty projections
    AliDebug(1, Form("Empty projection in %s", GetName()));
    delete h2; h2=NULL;
  }
  return h2;
}

//________________________________________________________
void AliTRDrecoTask::AliTRDrecoProjection::SetRebinStrategy(Int_t n, Int_t rebx[], Int_t reby[])
{
// define rebinning strategy for this projection
  fNrebin = n;
  fRebinX = new Int_t[n]; memcpy(fRebinX, rebx, n*sizeof(Int_t));
  fRebinY = new Int_t[n]; memcpy(fRebinY, reby, n*sizeof(Int_t));
}



