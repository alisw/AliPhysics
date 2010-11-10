#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSeqCollection.h"
#include "TObjArray.h"
#include "TObjArray.h"
#include "TChain.h"
#include "TMCProcess.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TNtuple.h"

#include "AliLog.h"
#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventCuts.h"
#include "AliMultiplicity.h"
#include "AliESDtrackCuts.h"
#include "AliVertex.h"
#include "AliFlowEventSimple.h"
#include "AliFlowVector.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskQAflow.h"

ClassImp(AliAnalysisTaskQAflow)

//________________________________________________________________________
AliAnalysisTaskQAflow::AliAnalysisTaskQAflow()
  : AliAnalysisTaskSE(),
    fOutput(NULL),
    fFillNtuple(kFALSE),
    fNtuple(NULL),
    fEventCuts(NULL),
    fTrackCuts(NULL)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskQAflow::AliAnalysisTaskQAflow(const char* name)
  : AliAnalysisTaskSE(name),
    fOutput(NULL),
    fFillNtuple(kFALSE),
    fNtuple(NULL),
    fEventCuts(NULL),
    fTrackCuts(NULL)
{
  // Constructor
  DefineInput(1, AliFlowEventSimple::Class());
  DefineOutput(1, TObjArray::Class());
}

//________________________________________________________________________
void AliAnalysisTaskQAflow::UserCreateOutputObjects()
{
  // Called once at the beginning
  fOutput=new TObjArray();
  fNtuple = new TNtuple("flowQAtree","flowQAtree","meanpt:qx:qy:mult:refmult:trig:tpcvtxx:tpcvtxy:tpcvtxz:ntrackletsA:ntrackletsC");
  //TDirectory* before = gDirectory->mkdir("before cuts","before cuts");
  //TDirectory* after = gDirectory->mkdir("after cuts","after cuts");
  TObjArray* before = new TObjArray();
  TObjArray* after = new TObjArray();
  fOutput->Add(before);
  fOutput->Add(after);

  //define histograms

  TH1* hist;  
  hist = new TH1I("tracklet_multiplicity","tracklet multiplicity",1000,0,20000);
  before->Add(hist); after->Add(hist->Clone()); //0
  hist = new TH1I("track_multiplicity", "standard TPC track multiplicity",1000,0,10000);
  before->Add(hist); after->Add(hist->Clone()); //1
  hist = new TH1I("refmult","refmult",1000, 0,10000);
  before->Add(hist); after->Add(hist->Clone()); //2
  hist = new TH1I("SPD_clusters","SPD clusters",1000,0,100000);
  before->Add(hist); after->Add(hist->Clone()); //3
  hist = new TH1D("primary_vertexZ","primary vertex z",100,-20,20);
  before->Add(hist); after->Add(hist->Clone()); //4
  hist = new TH1I("ITS_clusters_on_track", "ITS clusters on track", 8, 0, 8);
  before->Add(hist); after->Add(hist->Clone()); //5
  hist = new TH1I("TPC_clusters_on_track", "TPC clusters on track", 159, 1, 160);
  before->Add(hist); after->Add(hist->Clone()); //6
  hist = new TH1D("TPC_chi2_per_cluster","TPC #chi^{2}/cluster",100,0.0,5.0);
  before->Add(hist); after->Add(hist->Clone()); //7
  hist = new TH1D("DCA_xy","DCA xy", 1000, -5.0, 5.0 );
  before->Add(hist); after->Add(hist->Clone()); //8
  hist = new TH1D("DCA_z","DCA z", 1000, -5.0, 5.0 );
  before->Add(hist); after->Add(hist->Clone()); //9
  hist = new TH1D("phi_tracklets","#phi tracklets", 1000, 0.0, TMath::TwoPi() );
  before->Add(hist); after->Add(hist->Clone()); //10
  hist = new TH1D("phi_tracks","#phi tracks", 1000, 0.0, TMath::TwoPi() );
  before->Add(hist); after->Add(hist->Clone()); //11
  hist = new TH1D("eta_tracklets","#eta tracklets", 1000, -2.0, 2.0 );
  before->Add(hist); after->Add(hist->Clone()); //12
  hist = new TH1D("eta_tracks","#eta tracks", 1000, -2.0, 2.0 );
  before->Add(hist); after->Add(hist->Clone()); //13
  hist = new TH1D("TPC_vertex_z", "TPC vertex z", 100,-20.0,20.0);
  before->Add(hist); after->Add(hist->Clone()); //14
  hist = new TH1D("ptyield", "p_{t} spectrum", 10000,0.0,10.0);
  before->Add(hist); after->Add(hist->Clone()); //15
  
  //post data here as it doesn't change anyway (the pointer to list anyway)

  //restore dir add status
  PostData(0, fNtuple);
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskQAflow::UserExec(Option_t *)
{

  //get teh input data
  AliESDEvent* event = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!event)
  {
    AliFatal("no ESD event");
    return;
  }

  //TObjArray* before = ((TDirectory*)fOutput->At(0))->GetList();
  //TObjArray* after = ((TDirectory*)fOutput->At(1))->GetList();
  TObjArray* before = (TObjArray*)fOutput->At(0);
  TObjArray* after = (TObjArray*)fOutput->At(1);
  TH1* htrackletmultB = dynamic_cast<TH1*>(before->At(0));
  TH1* htrackletmultA = dynamic_cast<TH1*>(after->At(0));
  TH1* htrackmultB = dynamic_cast<TH1*>(before->At(1));
  TH1* htrackmultA = dynamic_cast<TH1*>(after->At(1));
  TH1* hrefmultB = dynamic_cast<TH1*>(before->At(2));
  TH1* hrefmultA = dynamic_cast<TH1*>(after->At(2));
  TH1* hspdclustersB = dynamic_cast<TH1*>(before->At(3));
  TH1* hspdclustersA = dynamic_cast<TH1*>(after->At(3));
  TH1* hprimvtxzB = dynamic_cast<TH1*>(before->At(4));
  TH1* hprimvtxzA = dynamic_cast<TH1*>(after->At(4));

  TH1* hITSclsB = dynamic_cast<TH1*>(before->At(5));
  TH1* hITSclsA = dynamic_cast<TH1*>(after->At(5));
  TH1* hTPCclsB = dynamic_cast<TH1*>(before->At(6));
  TH1* hTPCclsA = dynamic_cast<TH1*>(after->At(6));
  TH1* hTPCchi2B = dynamic_cast<TH1*>(before->At(7));
  TH1* hTPCchi2A = dynamic_cast<TH1*>(after->At(7));
  TH1* hdcaxyB = dynamic_cast<TH1*>(before->At(8));
  TH1* hdcaxyA = dynamic_cast<TH1*>(after->At(8));
  TH1* hdcazB = dynamic_cast<TH1*>(before->At(9));
  TH1* hdcazA = dynamic_cast<TH1*>(after->At(9));
  TH1* hphitrackletsB = dynamic_cast<TH1*>(before->At(10));
  TH1* hphitrackletsA = dynamic_cast<TH1*>(after->At(10));
  TH1* hphitracksB = dynamic_cast<TH1*>(before->At(11));
  TH1* hphitracksA = dynamic_cast<TH1*>(after->At(11));
  TH1* hetatrackletsB = dynamic_cast<TH1*>(before->At(12));
  TH1* hetatrackletsA = dynamic_cast<TH1*>(after->At(12));
  TH1* hetatracksB = dynamic_cast<TH1*>(before->At(13));
  TH1* hetatracksA = dynamic_cast<TH1*>(after->At(13));
  TH1* hprimvtxzTPCB = dynamic_cast<TH1*>(before->At(14));
  TH1* hprimvtxzTPCA = dynamic_cast<TH1*>(after->At(14));
  TH1* hptyieldB = dynamic_cast<TH1*>(before->At(15));
  TH1* hptyieldA = dynamic_cast<TH1*>(after->At(15));

  AliMultiplicity* tracklets = const_cast<AliMultiplicity*>(event->GetMultiplicity());
  Int_t ntracklets=0;
  Int_t nspdclusters=0;
  Int_t nspd1clusters=0;
  Int_t ntrackletsA=0;
  Int_t ntrackletsC=0;
  if (tracklets)
  {
    ntracklets = tracklets->GetNumberOfTracklets();
    nspdclusters = tracklets->GetNumberOfITSClusters(0,1);
    nspd1clusters = tracklets->GetNumberOfITSClusters(1);
    for (Int_t i=0; i<tracklets->GetNumberOfTracklets(); i++)
    {
      Bool_t pass=fTrackCuts->IsSelected(tracklets,i);
      Float_t phi=tracklets->GetPhi(i);
      Float_t eta=tracklets->GetEta(i);
      hphitrackletsB->Fill(phi); if (pass) hphitrackletsA->Fill(phi);
      hetatrackletsB->Fill(eta); if (pass) hetatrackletsA->Fill(eta); 
      if (eta>0) ntrackletsC++;
      else ntrackletsA++;
    }
    
  }
  //Int_t trackmult = AliESDtrackCuts::GetReferenceMultiplicity(event,kTRUE);
  AliFlowEventCuts newcuts;
  newcuts.SetRefMultCuts(fTrackCuts);
  Int_t trackmult = newcuts.RefMult(event);
  Int_t refmult = fEventCuts->RefMult(event);
  Bool_t passevent = fEventCuts->IsSelected(event);
  htrackletmultB->Fill(ntracklets); if (passevent) htrackletmultA->Fill(ntracklets); 
  htrackmultB->Fill(trackmult); if (passevent) htrackmultA->Fill( trackmult); 
  hrefmultB->Fill(refmult); if (passevent) hrefmultA->Fill( refmult); 
  hspdclustersB->Fill(nspdclusters); if (passevent) hspdclustersA->Fill( nspdclusters);
  const AliVertex* vertex = event->GetPrimaryVertex();
  Float_t vtxz=0.0;
  Float_t vtxx=0.0;
  Float_t vtxy=0.0;
  if (vertex)
  {
    vtxz = vertex->GetZ();
    vtxx = vertex->GetX();
    vtxy = vertex->GetY();
    hprimvtxzB->Fill(vtxz); if (passevent) hprimvtxzA->Fill(vtxz);
  }
  const AliVertex* vertextpc = event->GetPrimaryVertexTPC();
  Float_t vtxTPCx=0.0;
  Float_t vtxTPCy=0.0;
  Float_t vtxTPCz=0.0;
  if (vertextpc)
  {
    vtxTPCx = vertextpc->GetX();
    vtxTPCy = vertextpc->GetY();
    vtxTPCz = vertextpc->GetZ();
    hprimvtxzTPCB->Fill(vtxTPCz); if (passevent) hprimvtxzTPCA->Fill(vtxTPCz);
  }
  fTrackCuts->SetEvent(event);
  Float_t meanpt=0;
  Int_t ntracks=fTrackCuts->GetNumberOfInputObjects();
  for (Int_t i=0; i<ntracks; i++)
  {
    TObject* obj = fTrackCuts->GetInputObject(i);
    Bool_t pass = fTrackCuts->IsSelected(obj,i);
    Float_t dcaxy=0.0;
    Float_t dcaz=0.0;
    Float_t tpcchi2=0.0;
    Float_t tpcchi2percls=0.0;
    Int_t ntpccls=0;
    Int_t nitscls=0;
    Float_t eta=0.0;
    Float_t phi=0.0;
    Float_t pt=0.0;
    AliESDtrack* track = dynamic_cast<AliESDtrack*>(fTrackCuts->GetTrack());
    if (track)
    {
      track->GetImpactParameters(dcaxy,dcaz);
      tpcchi2=track->GetTPCchi2();
      ntpccls=track->GetTPCNcls();
      eta=track->Eta();
      phi=track->Phi();
      pt=track->Pt();
      meanpt+=pt;
      tpcchi2percls= (ntpccls==0)?0.0:tpcchi2/ntpccls;
      nitscls = track->GetNcls(0);
      hITSclsB->Fill(nitscls); if (pass) hITSclsA->Fill( nitscls);
      hTPCclsB->Fill(ntpccls); if (pass) hTPCclsA->Fill( ntpccls);
      hTPCchi2B->Fill(tpcchi2percls); if (pass) hTPCchi2A->Fill( tpcchi2percls);
      hdcaxyB->Fill(dcaxy); if (pass) hdcaxyA->Fill( dcaxy);
      hdcazB->Fill(dcaz); if (pass) hdcazA->Fill(dcaz);
      hetatracksB->Fill(eta); if (pass) hetatracksA->Fill(eta);
      hphitracksB->Fill(phi); if (pass) hphitracksA->Fill(phi);
      hptyieldB->Fill(pt); if (pass) hptyieldA->Fill(pt);
    }
  }
  meanpt = meanpt/ntracks;

  ///////////////////////////////////////////////////////////////////////
  //the ntuple part/////////////
  AliFlowEventSimple* flowevent = dynamic_cast<AliFlowEventSimple*>(GetInputData(1));
  if (flowevent && fFillNtuple)
  {
    AliFlowVector qvec = flowevent->GetQ(2);
    Double_t qx = qvec.X();
    Double_t qy = qvec.Y();
    TString triggers = event->GetFiredTriggerClasses();
    Int_t trig=0;
    if (triggers.Contains("CMBAC-B-NOPF-ALL")) trig+=1;
    if (triggers.Contains("CMBS2C-B-NOPF-ALL")) trig+=10;
    if (triggers.Contains("CMBS2A-B-NOPF-ALL")) trig+=100;
    fNtuple->Fill(meanpt,qx,qy,trackmult,refmult,trig,vtxTPCx,vtxTPCy,vtxTPCz,ntrackletsA,ntrackletsC);
  }//if flowevent
}

//________________________________________________________________________
void AliAnalysisTaskQAflow::Terminate(Option_t *)
{
  
}

//________________________________________________________________________
AliAnalysisTaskQAflow::~AliAnalysisTaskQAflow()
{
  //dtor
  delete fTrackCuts;
  delete fEventCuts;
}
