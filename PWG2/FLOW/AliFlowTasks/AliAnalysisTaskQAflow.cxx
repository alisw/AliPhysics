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
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliESDtrack.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventCuts.h"
#include "AliMultiplicity.h"
#include "AliESDtrackCuts.h"
#include "AliVertex.h"
#include "AliFlowEventSimple.h"
#include "AliFlowEvent.h"
#include "AliFlowVector.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

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
  DefineOutput(2,TNtuple::Class());
}

//________________________________________________________________________
void AliAnalysisTaskQAflow::UserCreateOutputObjects()
{
  // Called once at the beginning
  fOutput=new TObjArray();
  fNtuple = new TNtuple("flowQAtree","flowQAtree","mpt:qx:qy:mul:rmul:phys:vtxtpcx:vtxtpcy:vtxtpcz:ntra:ntrc:mv0a:mv0c:zdcp1:zdcn1:zdcp2:zdcn2:zdcpart1:zdcpart2:t1:t2:t3:t4:t5:vtxspdx:vtxspdy:vtxspdz:vtxx:vtxy:vtxz:rawmeanpt:maxpt");
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
  hist = new TH1D("ptyieldplus", "p_{t} spectrum +", 10000,0.0,10.0);
  before->Add(hist); after->Add(hist->Clone()); //16
  hist = new TH1D("ptyieldneg", "p_{t} spectrum -", 10000,0.0,10.0);
  before->Add(hist); after->Add(hist->Clone()); //17
  hist = new TH1D("charges", "charge distribution", 5,-2.5,2.5);
  before->Add(hist); after->Add(hist->Clone()); //17
  
  //post data here as it doesn't change anyway (the pointer to list anyway)

  //restore dir add status
  PostData(0, fNtuple);
  PostData(1, fOutput);
  PostData(2, fNtuple);
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
  TH1* hptyieldplusB = dynamic_cast<TH1*>(before->At(16));
  TH1* hptyieldplusA = dynamic_cast<TH1*>(after->At(16));
  TH1* hptyieldnegB = dynamic_cast<TH1*>(before->At(17));
  TH1* hptyieldnegA = dynamic_cast<TH1*>(after->At(17));
  TH1* hchargesB = dynamic_cast<TH1*>(before->At(18));
  TH1* hchargesA = dynamic_cast<TH1*>(after->At(18));

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
  const AliVertex* vertexspd = event->GetPrimaryVertexSPD();
  Float_t vtxSPDx=0.0;
  Float_t vtxSPDy=0.0;
  Float_t vtxSPDz=0.0;
  if (vertexspd)
  {
    vtxSPDx = vertexspd->GetX();
    vtxSPDy = vertexspd->GetY();
    vtxSPDz = vertexspd->GetZ();
  }
  AliESDVZERO* vzero=event->GetVZEROData();
  AliESDZDC* zdc=event->GetESDZDC();
  Float_t mv0a=vzero->GetMTotV0A();
  Float_t mv0c=vzero->GetMTotV0C();
  Float_t zdcp1=zdc->GetZDCP1Energy();
  Float_t zdcn1=zdc->GetZDCN1Energy();
  Float_t zdcp2=zdc->GetZDCP2Energy();
  Float_t zdcn2=zdc->GetZDCN2Energy();
  Float_t zdcpart1=zdc->GetZDCPartSideA();
  Float_t zdcpart2=zdc->GetZDCPartSideC();
  fTrackCuts->SetEvent(event);
  Float_t meanpt=0;
  Float_t rawmeanpt=0;
  Float_t maxpt=0.0;
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
    Short_t charge=0;
    AliESDtrack* track = dynamic_cast<AliESDtrack*>(fTrackCuts->GetTrack());
    if (track)
    {
      track->GetImpactParameters(dcaxy,dcaz);
      tpcchi2=track->GetTPCchi2();
      ntpccls=track->GetTPCNcls();
      eta=track->Eta();
      phi=track->Phi();
      pt=track->Pt();
      if (pt>maxpt) maxpt=pt;
      charge=track->Charge();
      rawmeanpt+=pt;
      if (pass) meanpt+=pt;
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
      if (charge>0) {hptyieldplusB->Fill(pt); if (pass) hptyieldplusA->Fill(pt);}
      if (charge<0) {hptyieldnegB->Fill(pt); if (pass) hptyieldnegA->Fill(pt);}
      if (charge<0) {hchargesB->Fill(charge); if (pass) hchargesA->Fill(charge);}
    }
  }
  meanpt = meanpt/ntracks;
  
  Bool_t isSelectedEventSelection = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  
  ///////////////////////////////////////////////////////////////////////
  //the ntuple part/////////////
  if (fFillNtuple)
  {
    Double_t qx = 0.0;
    Double_t qy = 0.0;
    AliFlowVector qvec;
    AliFlowEventSimple* flowevent = dynamic_cast<AliFlowEventSimple*>(GetInputData(1));
    if (!flowevent) flowevent = new AliFlowEvent(fTrackCuts,fTrackCuts);
    qvec = flowevent->GetQ(2);
    qx = qvec.X();
    qy = qvec.Y();
    Float_t x[32];
    x[0]=meanpt; x[1]=qx; x[2]=qy; x[3]=trackmult; x[4]=refmult; x[5]=(isSelectedEventSelection)?1:0;
    x[6]=vtxTPCx; x[7]=vtxTPCy; x[8]=vtxTPCz; x[9]=ntrackletsA; x[10]=ntrackletsC;
    x[11]=mv0a; x[12]=mv0c; x[13]=zdcp1; x[14]=zdcn1; x[15]=zdcp2; x[16]=zdcn2;
    x[17]=zdcpart1; x[18]=zdcpart2;
    x[19]=0; if (event->IsTriggerClassFired("CMBAC-B-NOPF-ALL"))  x[19]=1;
    x[20]=0; if (event->IsTriggerClassFired("CMBS2C-B-NOPF-ALL")) x[20]=1;
    x[21]=0; if (event->IsTriggerClassFired("CMBS2A-B-NOPF-ALL")) x[21]=1;
    x[22]=0; if (event->IsTriggerClassFired("CMBAC-A-NOPF-ALL"))  x[22]=1;
    x[23]=0; if (event->IsTriggerClassFired("CMBAC-C-NOPF-ALL"))  x[23]=1;
    x[24]=vtxSPDx;
    x[25]=vtxSPDy;
    x[26]=vtxSPDz;
    x[27]=vtxx;
    x[28]=vtxy;
    x[29]=vtxz;
    x[30]=rawmeanpt;
    x[31]=maxpt;
    fNtuple->Fill(x);
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
