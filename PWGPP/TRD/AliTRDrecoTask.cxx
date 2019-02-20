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

#include <TClass.h>
#include <TMethod.h>
#include <TMethodCall.h>
#include <TMethodArg.h>
#include <TFile.h>
#include <TChain.h>
#include <TList.h>
#include <TMap.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TDirectory.h>
#include <TTreeStream.h>
#include <TBox.h>
#include <TLatex.h>
#include <TVectorT.h>

#include <AliLog.h>
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliExternalTrackParam.h>
#include <AliTRDtrackV1.h>
#include <AliTRDpidUtil.h>

#include "AliTRDchmbInfo.h"
#include "AliTRDeventInfo.h"
#include "AliTRDtrendingManager.h"
#include "AliTRDrecoTask.h"

ClassImp(AliTRDrecoTask)

Float_t AliTRDrecoTask::fgPt[AliTRDrecoTask::fgNPt+1] = {0.};
TTreeSRedirector* AliTRDrecoTask::fgDebugStream(NULL);
TH1* AliTRDrecoTask::fgProjector(NULL);
//_______________________________________________________
AliTRDrecoTask::AliTRDrecoTask()
  : AliAnalysisTaskSE()
  ,fNRefFigures(0)
  ,fDets(NULL)
  ,fDetsV(NULL)
  ,fContainer(NULL)
  ,fEvent(NULL)
  ,fTracks(NULL)
  ,fOnlTracklets(NULL)
  ,fClusters(NULL)
  ,fkClusters(NULL)
  ,fkTrack(NULL)
  ,fkMC(NULL)
  ,fkESD(NULL)
  ,fSpecies(-6)
  ,fTriggerSlot(0)
  ,fPt(-1.)
  ,fPhi(0.)
  ,fEta(0.)
  ,fNpt(-1)
  ,fTriggerList(NULL)
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
  ,fOnlTracklets(NULL)
  ,fClusters(NULL)
  ,fkClusters(NULL)
  ,fkTrack(NULL)
  ,fkMC(NULL)
  ,fkESD(NULL)
  ,fSpecies(-6)
  ,fTriggerSlot(0)
  ,fPt(-1.)
  ,fPhi(0.)
  ,fEta(0.)
  ,fNpt(-1)
  ,fTriggerList(NULL)
  ,fPlotFuncList(NULL)
  ,fDetFuncList(NULL)
  ,fRunTerminate(kFALSE)
{
// Constructor for all derived performance tasks

  SetTitle(title);
  snprintf(fNameId, 10, "no name");
  DefineInput (1, TObjArray::Class()); // track list
  DefineInput (2, AliTRDeventInfo::Class()); // event info object
  DefineInput (3, TObjArray::Class()); // online tracklets list object
  DefineInput (4, TObjArray::Class()); // cluster list object
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

  if(fgProjector){
    delete fgProjector;
    fgProjector = NULL;
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
  if(fTriggerList){fTriggerList->Delete(); delete fTriggerList;}

  if(fContainer && !(AliAnalysisManager::GetAnalysisManager() && AliAnalysisManager::GetAnalysisManager()->IsProofMode())){
    if(fContainer->IsOwner()) fContainer->Delete();
    delete fContainer;
    fContainer = NULL;
  }

/*  if(fgTrendPoint){
    TFile::Open("TRD.PerformanceTrend.root", "UPDATE");
    fgTrendPoint->Write();
    delete fgTrendPoint;
    fgTrendPoint=NULL;
    gFile->Close();
  }*/
}

//_______________________________________________________
Int_t AliTRDrecoTask::GetNRefFigures() const  
{ 
  if(!fNRefFigures) AliWarning("No reference plots available.");
  return fNRefFigures; 
} 

//____________________________________________________________________
Int_t AliTRDrecoTask::GetPtBin(Float_t pt)
{
// Get significant (very low, low, medium, high, very high) pt bin

  Int_t ipt(0);
  while(ipt<fNpt){
    if(pt<fgPt[ipt]) break;
    ipt++;
  }
  return ipt-1;
}

//_______________________________________________________
Bool_t AliTRDrecoTask::MakeMomSegmentation()
{
  switch(fNpt){
  case fgNPt:
    fgPt[0]=0.3;
    for(Int_t j(1); j<=fgNPt; j++) fgPt[j]=fgPt[j-1]+(TMath::Exp(j*j*2.e-3)-1.);
    AliDebug(1, "Using debug momentum segmentation");
    break;
  case 3:
    fgPt[0]=0.5; fgPt[1]=0.8; fgPt[2]=1.5; fgPt[3]=5.;
    AliDebug(1, "Using default momentum segmentation");
    break;
  default:
    AliError(Form("Momentum segmentation %d not supported.", fNpt));
    fNpt=0;
    return kFALSE;
  }
  return kTRUE;
}

//_______________________________________________________
void AliTRDrecoTask::UserCreateOutputObjects()
{
  if(!HasFunctorList()) InitFunctorList();
  if(DebugLevel()) fNpt = fgNPt;
  else fNpt = 3;
  MakeMomSegmentation();

  fContainer = Histos();
  PostData(1, fContainer);
}

//_______________________________________________________
void AliTRDrecoTask::UserExec(Option_t *)
{
// Loop over Plot functors published by particular tasks

  fTracks   = dynamic_cast<TObjArray *>(GetInputData(1));
  fEvent    = dynamic_cast<AliTRDeventInfo *>(GetInputData(2));
  fTriggerSlot=0;
  if(fTriggerList && fEvent){
    for(Int_t itrig(0); itrig<fTriggerList->GetEntries(); itrig++){
      if(!fEvent->GetFiredTriggerClasses().Contains(((TObjString*)(*fTriggerList)[itrig])->GetName())) continue;
      //printf("\"%s\" selected\n", ((TObjString*)(*fTriggerList)[itrig])->GetName());
      SETBIT(fTriggerSlot,itrig);
    }
    if(!fTriggerSlot){
      AliDebug(2, Form("Triggers[%s] not used for %s", fEvent->GetFiredTriggerClasses().Data(),  GetName()));
      return;
    }
  }
  fOnlTracklets = dynamic_cast<TObjArray*>(GetInputData(3)); // link online tracklets
  fClusters  = dynamic_cast<TObjArray*>(GetInputData(4)); // link offline clusters

  if(!fPlotFuncList){
    AliWarning("No track functor list defined for the task");
    return;
  }
  if(!fEvent || !fTracks) return;
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
Bool_t AliTRDrecoTask::PutTrendValue(const Char_t *name, Double_t val, Double_t err)
{
// Generic publisher for trend values

  AliTRDtrendingManager *tm = AliTRDtrendingManager::Instance();
  if(!tm){
    AliError("Wrong usage of the trending functionality. Could not instantiate AliTRDtrendingManager singleton.");
    return kFALSE;
  }
  tm->AddValue(Form("%s_%s", GetName(), name), val, err);
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
    gFile->Close();
    return kFALSE;
  }
  if(!(fContainer = (TObjArray*)gDirectory->Get(GetName()))){
    AliWarning("Missing histogram container.");
    gFile->Close();
    return kFALSE;
  }
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
    gFile->Close();
    return kFALSE;
  }
  TObjArray *info = NULL;
  if(!(info = (TObjArray*)gDirectory->Get("TRDinfoGen"))){
    AliWarning("Missing TRDinfoGen container.");
    gFile->Close();
    return kFALSE;
  }

  if(info->FindObject("Chambers Status"))
    fDets = (TObjArray*)((TObjArray*)info->FindObject("Chambers Status"))->Clone();

  if(!fDets){
    if(!info->At(4) || strcmp("TObjArray", info->At(4)->IsA()->GetName())) AliError("Looking for old style chamber status map. Failed.");
    else {
      AliWarning("Looking for old style chamber status map.");
      fDetsV = (TObjArray*)((TObjArray*)info->At(4))->Clone();
    }
  }
  gFile->Close();
  info->Delete(); delete info;
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
    AliWarning("NEW Detector map and status not available. Try OLD");
    MakeDetectorPlotOLD(ly, opt);
    return;
  }
  AliTRDchmbInfo *ci(NULL);
  for(Int_t idet(0); idet<fDets->GetEntriesFast(); idet++){
    if(!(ci = (AliTRDchmbInfo*)fDets->At(idet))) continue;
    if(AliTRDgeometry::GetLayer(ci->GetDetector()) != ly) continue;
    ci->Draw(opt);
  }
}

//_______________________________________________________
void AliTRDrecoTask::MakeDetectorPlotOLD(Int_t ly, const Option_t *opt)
{
// Draw chamber boundaries in eta/phi plots with misalignments
// based on info collected by AliTRDinfoGen OLD data storage

  if(!fDetsV){
    AliError("OLD Detector map and status not available.");
    return;
  }
  if(!fDetsV->GetEntries()){
    AliError("OLD Detector map and status not filled.");
    return;
  }

  Float_t xmin(0.), xmax(0.);
  TBox *gdet = new TBox();
  gdet->SetLineColor(kBlack);gdet->SetFillColor(kBlack);
  Int_t style[] = {0, 3003};
  for(Int_t idet(0); idet<540; idet++){
    if(idet%6 != ly) continue;
    TVectorF *det((TVectorF*)fDetsV->At(idet));
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
  TLatex *sm = new TLatex(); sm->SetTextAlign(22);sm->SetTextColor(kBlack); sm->SetTextFont(32);sm->SetTextSize(0.03);
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
  if(DebugLevel()>=1 && !fgDebugStream){
    AliInfo(Form("Debug Level for Task %s set to %d", GetName(), level));
    TDirectory *savedir = gDirectory;
    fgDebugStream = new TTreeSRedirector("TRD.DebugPerformance.root", "RECREATE");
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
Float_t AliTRDrecoTask::SetNormZ(TH2 *h2, Int_t bxmin, Int_t bxmax, Int_t bymin, Int_t bymax, Float_t thr)
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
  return s;
}

//________________________________________________________
void AliTRDrecoTask::SetRangeZ(TH2 *h2, Float_t min, Float_t max, Float_t thr, Float_t scale)
{
// Set range on Z axis such to avoid outliers. Optionally a scale factor can be applied 

  Float_t c(0.), dz(1.e-3*(max-min));
  for(Int_t ix(1); ix<=h2->GetXaxis()->GetNbins(); ix++){
    for(Int_t iy(1); iy<=h2->GetYaxis()->GetNbins(); iy++){
      c = h2->GetBinContent(ix, iy)*scale;
      if(c<thr) continue;
      if(c<=min) h2->SetBinContent(ix, iy, min+dz);
      else h2->SetBinContent(ix, iy, c);
    }
  }
  h2->GetZaxis()->SetRangeUser(min, max);
}

//________________________________________________________
Float_t AliTRDrecoTask::GetMeanStat(TH1 *h, Float_t cut, Int_t opt, Float_t *sigma)
{
// Return mean number of entries/bin of histogram "h".
// If optionally sigma is allocated than it is also filled with sigma paramter of the gauss fit 
//
// Option "opt" is given the following values are accepted:
//   -1 : consider only entries less than "cut"
//   1  : consider only entries greater than "cut"
//   0  : no "cut" [dafault]
// Error codes
//   -999. : statistics too low [20]
//   -998. : fit failed

  const Int_t kvd(100000);
  Float_t v[kvd];
  Int_t nbx(h->GetNbinsX()), nby(h->GetNbinsY()), nbz(h->GetNbinsZ());
  Int_t nv(0); Float_t xmin(1.e5), xmax(-xmin);
  for(Int_t ix(1); ix<=nbx; ix++){
    for(Int_t iy(1); iy<=nby; iy++){
      for(Int_t iz(1); iz<=nbz; iz++){
        Float_t c = h->GetBinContent(ix, iy, iz);
        if(opt*(c-cut)<0.) continue;
        v[nv++] = c;
        if(c<xmin) xmin = c;
        if(c>xmax) xmax = c;
        if(nv==kvd){
          printf("W - AliTRDrecoTask::GetMeanStat() :: Unreliable results for %s[%s]. Statical allocation exceeded.\n", h->GetName(), h->GetTitle());
          break;
        }
      }
      if(nv==kvd) break;
    }
    if(nv==kvd) break;
  }
  if(nv<10){
    //printf("W - AliTRDrecoTask::GetMeanStat() :: Failed for %s[%s]. Statical undefined [%d].\n", h->GetName(), h->GetTitle(), nv);
    return -999.;
  }
  if(fgProjector) delete fgProjector;
  fgProjector = new TH1F("hProjector", "", 20, 0.5*(3*xmin-xmax), 0.5*(3*xmax - xmin));
  for(Int_t iv(0); iv<nv; iv++) fgProjector->Fill(v[iv]);
  TF1 f("f", "gaus", xmin, xmax);
  f.SetParameter(0, fgProjector->Integral());
  f.SetParameter(1, fgProjector->GetMean()); f.SetParLimits(1, xmin, xmax);
  f.SetParameter(2, fgProjector->GetRMS());
  if(fgProjector->Fit(&f, "WQ0", "goff")) return -998.;
  if(sigma) *sigma = f.GetParameter(2);
  return f.GetParameter(1);
}

//________________________________________________________
Int_t AliTRDrecoTask::Rebin(TH2 *h, Int_t n, Int_t rebinX[], Int_t rebinY[], Int_t nstat)
{
// Rebin histo "h" according to "rebinning" strategy "n" steps such to obtain mean statistics per bin over "nstat"

  Int_t irebin(0);
  while(irebin<n && GetMeanStat(h, .5, 1)<nstat){
    h->Rebin2D(rebinX[irebin], rebinY[irebin]);
    irebin++;
  }
  return irebin;
}

//________________________________________________________
AliTRDrecoTask::AliTRDrecoProjection::AliTRDrecoProjection()
  :TNamed()
  ,fH(NULL)
  ,fNrebin(0)
{
  // constructor
  fRebin[0] = NULL;fRebin[1] = NULL;
  memset(fAx, 0, 3*sizeof(Int_t));
  memset(fRange, 0, 4*sizeof(Float_t));
}

//________________________________________________________
AliTRDrecoTask::AliTRDrecoProjection::~AliTRDrecoProjection()
{
  // destructor
  if(fH) delete fH;
  if(fRebin[0]) delete [] fRebin[0];
  if(fRebin[1]) delete [] fRebin[1];
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
  if(fNrebin){fNrebin=0; delete [] fRebin[0]; delete [] fRebin[1];}
  if(rhs.fNrebin) SetRebinStrategy(rhs.fNrebin, rhs.fRebin[0], rhs.fRebin[1]);
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
  //fH->AddBinContent(fH->GetBin(bin[fAx[0]],bin[fAx[1]],bin[fAx[2]]), Int_t(v));
  TAxis *ax(fH->GetXaxis()),  *ay(fH->GetYaxis()),  *az(fH->GetZaxis());
  fH->Fill(ax->GetBinCenter(bin[fAx[0]]), ay->GetBinCenter(bin[fAx[1]]), az->GetBinCenter(bin[fAx[2]]), v);
}

//________________________________________________________
Double_t AliTRDrecoTask::AliTRDrecoProjection::GetTrendValue(const Int_t mid, Double_t *e, Double_t *s, Double_t *se) const
{
//   Return result of fitting the main distribution (represented on the z axis) with the function selected
// "mid". Optionally return the Mean and RMS of the distribution pointing to "m" and "s"

  if(!fH){
    AliDebug(1, Form("Missing 3D in %s", GetName()));
    return -999.;
  }
  TH1 *h1s(NULL);
  if(!(h1s = (TH1D*)fH->Project3D("z"))){
    AliDebug(1, Form("Failed Project3D(\"z\") in %s", GetName()));
    return -999.;
  }
  Int_t ne((Int_t)h1s->Integral());
  if(ne<30){
    AliDebug(1, Form("Statistics too low[%2d] in %s", ne, GetName()));
    delete h1s;
    return -999.;
  }
  TAxis *az(h1s->GetXaxis());
  Float_t mn(h1s->GetMean()), rms(h1s->GetRMS()),
          v(mn),  // main trending value (mean, mu, MPV)
          ve(rms),// dispersion (RMS, sigma, landau 2nd param)
          ev(h1s->GetMeanError()), // error on v
          eve(h1s->GetRMSError());// error on ve
  if(mid==1){
    TF1 fg("fg", "gaus", az->GetXmin(), az->GetXmax());
    fg.SetParameter(0, Float_t(ne)); fg.SetParameter(1, mn); fg.SetParameter(2, rms);
    h1s->Fit(&fg, "WQ0");
    v = fg.GetParameter(1); ev = fg.GetParError(1);
    ve= fg.GetParameter(2); eve= fg.GetParError(2);
  } else if (mid==2) {
    TF1 fl("fl", "landau", az->GetXmin(), az->GetXmax());
    fl.SetParameter(0, Float_t(ne)); fl.SetParameter(1, mn); fl.SetParameter(2, rms);
    h1s->Fit(&fl, "WQ0");
    v = fl.GetMaximumX(); ev = fl.GetParError(1);
    ve= fl.GetParameter(2);eve= fl.GetParError(2);
  }
  if(e)  *e  = ev;
  if(s)  *s  = ve;
  if(se) *se = eve;
  AliDebug(2, Form("[%d] %s(%4d) = M{%f+-%f} S{%f+-%f}", mid, fH->GetName(), (Int_t)h1s->Integral(), v, ev, ve, eve));

  delete h1s;
  return v;
}

//________________________________________________________
TH2* AliTRDrecoTask::AliTRDrecoProjection::Projection2Dbin(Int_t bin, Bool_t mc)
{
// dumb 2D projection for bin including under/over flow. Default all [bin==-1]

  TAxis *ax(fH->GetXaxis()), *ay(fH->GetYaxis()), *az(fH->GetZaxis());
  Int_t nbins(az->GetNbins());
  TH2F *h2(NULL);
  if(bin<0) h2 = new TH2F(Form("%s_2D", fH->GetName()),
                Form("%s;%s;%s;Entries", fH->GetTitle(), ax->GetTitle(), ay->GetTitle()),
                ax->GetNbins(), ax->GetXmin(), ax->GetXmax(),
                ay->GetNbins(), ay->GetXmin(), ay->GetXmax());
  else h2 = new TH2F(Form("%s%d_2D", fH->GetName(), bin),
                Form("%s | #it{%4.2f<=p_{t}^{%s}[GeV/c]<%4.2f};%s;%s;Entries", fH->GetTitle(),
                bin?fgPt[bin-1]:0., mc?"MC":"", bin>nbins?99.99:fgPt[bin], ax->GetTitle(), ay->GetTitle()),
                ax->GetNbins(), ax->GetXmin(), ax->GetXmax(),
                ay->GetNbins(), ay->GetXmin(), ay->GetXmax());
  for(Int_t ix(1); ix<=ax->GetNbins(); ix++){
    for(Int_t iy(1); iy<=ay->GetNbins(); iy++){
      Int_t ibin = h2->GetBin(ix, iy);
      for(Int_t iz(0); iz<=az->GetNbins()+1; iz++){
        if(bin<0) h2->AddBinContent(ibin, fH->GetBinContent(ix, iy, iz));
        else if(bin==iz){
          h2->AddBinContent(ibin, fH->GetBinContent(ix, iy, iz));
          break;
        }
      }
    }
  }
  return h2;
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
  TH2D *h2s(NULL);
  if(!(h2s = (TH2D*)fH->Project3D("yx"))){
    AliDebug(1, Form("Failed Project3D(\"yx\") in %s", GetName()));
    return NULL;
  }
  Int_t irebin(Rebin(h2s, fNrebin, fRebin[0], fRebin[1], nstat)), dxBin(1), dyBin(1);
  for(Int_t ir(0); ir<irebin; ir++){dxBin*=fRebin[0][ir]; dyBin*=fRebin[1][ir];}
  Int_t nx(h2s->GetNbinsX()), ny(h2s->GetNbinsY());
  // save a copy of the original distribution
  if(!del) h2s->SetNameTitle(Form("%sEn", fH->GetName()), Form("%s statistics", fH->GetTitle()));
  else delete h2s;
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
  tokenTitle->Delete(); delete tokenTitle;
  h2->SetContour(ncol);
  h2->GetZaxis()->CenterTitle();
  h2->GetZaxis()->SetTitleOffset(1.4);
  h2->GetZaxis()->SetRangeUser(fRange[0], fRange[1]);
  AliDebug(2, Form("%s[%s] nx[%d] ny[%d]", h2->GetName(), h2->GetTitle(), nx, ny));
  for(Int_t iy(0); iy<ny; iy++){
    for(Int_t ix(0); ix<nx; ix++){
      h = fH->ProjectionZ(Form("%s_z", h2->GetName()), ix*dxBin+1, (ix+1)*dxBin, iy*dyBin+1, (iy+1)*dyBin);
      Int_t ne((Int_t)h->Integral());
      //printf("  %s :: x[%2d %2d] y[%2d %2d] ne[%4d] nstat[%2d %2d]\n", h2->GetName(), ix*dxBin+1, (ix+1)*dxBin, iy*dyBin+1, (iy+1)*dyBin, ne, nstat, nstat/4);
      if(ne<nstat/4){
        h2->SetBinContent(ix+1, iy+1, -999);
        h2->SetBinError(ix+1, iy+1, 1.);
        n++;
      }else{
        // redo the projection by adding 1 bin @ left and 1 bin @ right for smoothing
        h = fH->ProjectionZ(Form("%s_z", h2->GetName()), ix*dxBin, (ix+1)*dxBin+1, iy*dyBin, (iy+1)*dyBin+1);
        Float_t v(h->GetMean()), ve(h->GetRMS());
        if(ne<h->GetNbinsX()) h->Rebin(2);
        if(mid==1){
          TF1 fg("fg", "gaus", az->GetXmin(), az->GetXmax());
          fg.SetParameter(0, h->GetBinContent(h->GetMaximumBin())); 
          fg.SetParameter(1, v);fg.SetParLimits(1, v-0.5*ve, v+0.5*ve); 
          fg.SetParameter(2, ve);fg.SetParLimits(2, 0.5*ve, 1.5*ve);
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
  fRebin[0] = new Int_t[n]; memcpy(fRebin[0], rebx, n*sizeof(Int_t));
  fRebin[1] = new Int_t[n]; memcpy(fRebin[1], reby, n*sizeof(Int_t));
}

//________________________________________________________
void AliTRDrecoTask::SetTriggerList(const Char_t *tl)
{
// Store list of triggers to be monitored
  TString stl(tl);
  if(fTriggerList){ fTriggerList->Delete(); delete fTriggerList;}
  TObjArray *atl = stl.Tokenize(" ");
  fTriggerList = (TObjArray*)atl->Clone("");
  atl->Delete(); delete atl;
  AliInfo("Running only for triggers::");
  fTriggerList->Print();
}


