// "$MYDATA/gsi/gsi.root"

const char* fileName = "$MYDATA/gsi02/gsi02a.root";
const char* listName = "Cont_cloizide_DhcAna";

// const char* fileName = "dhctask_output_etamax1.root";
// const char* listName = "dhclist";

// From AliDhcTask.h...keep current!
const Int_t nPairAxes = 6;
const Int_t nEvtAxes  = 3;
const Int_t nTrkAxes  = 2;
//const Int_t nExtra    = 3;

TString sCorrAxes[] = {"#Delta#eta", "assoc. p_{T}", "trig. p_{T}", 
		       "Centrality", "#Delta#phi (rad)",
		       "z-vertex (cm)", "Charge Comb."};
TString sEvtAxes[] = {"z-vertex", "V0M Centrality (%)", "CL1 Centrality"};
TString sTrkAxes[] = {"phi", "eta"};

enum ePairHistAxes  {kDeta, kPtAssc, kPtTrig, kCent, kDphi,
		     kZvtx};
enum eEventHistAxes {kZvtxEvt, kCentV0M, kCentCL1};
enum eTrackHistAxes {kPhiTrk, kEtaTrk};

TObjArray* hists = new TObjArray();

//gROOT->LoadMacro("../../../common/ProjectionUtils.C+");
void CheckOutput()
{
  gStyle->SetOptTitle(0);

  TFile* inFile = new TFile(fileName, "read");
  TList* dhcList = inFile->Get(listName);
  if (!dhcList) {
    Error("","No List %s in %s", fileName, listName);
    return;
  }

  THnSparse* hs = dhcList->FindObject("fHS");
  THnSparse* hm = dhcList->FindObject("fHM");

  if (!hs) {
    Error("","No hs in list %s", listName);
    return;
  }

  TH2F* hEvt = dhcList->FindObject("fHEvt");
  TH2F* hTrk = dhcList->FindObject("fHTrk");
  hists->Add(hEvt);
  hists->Add(hTrk);
  TH1D* heta = hTrk->ProjectionY();
  hists->Add(heta);

  TH1D* hzvtx = hEvt->ProjectionX("hzvtx", 1, 90);
  hists->Add(hzvtx);
  TH1D* hcent = hEvt->ProjectionY("hcent");
  hists->Add(hcent);
  // Correlation projections

  // hs->GetAxis(kZvtx)->SetRangeUser(-1.9, 1.9);
  // hm->GetAxis(kZvtx)->SetRangeUser(-1.9, 1.9);

  // hs->GetAxis(kCent)->SetRangeUser(0, 4.99);
  // hm->GetAxis(kCent)->SetRangeUser(0, 4.99);

  TH2D* sig = hs->Projection(kDeta,kDphi); sig->SetName("sig");
  TH2D* bkg = hm->Projection(kDeta,kDphi); sig->SetName("bkg");
  sig->Sumw2(); bkg->Sumw2();
  sig->Scale(1./sig->Integral()); 
  bkg->Scale(1./bkg->Integral()); 
  sig->Divide(bkg);
  sig->GetYaxis()->SetRangeUser(-1.999, 1.999); 
  TH2D* hcz = hs->Projection(kZvtx, kCent);
  hcz->SetName("zvtx_vs_cent");

  // 1D
  TH1D* hzr = 0;

  TH1D *hsVar[nPairAxes], *hmVar[nPairAxes], *evtVar[nEvtAxes], 
    *trkVar[nTrkAxes];

  for (int i=0; i<nPairAxes; i++) {

    hsVar[i] = hs->Projection(i);
    hsVar[i]->SetFillColor(kRed-7);
    hsVar[i]->Sumw2();
    hsVar[i]->SetTitle(Form("title;%s;same-evt. pairs", sCorrAxes[i].Data()));
    hmVar[i] = hm->Projection(i);
    hmVar[i]->SetFillColor(kAzure-9);
    hmVar[i]->Sumw2();
    hmVar[i]->SetTitle(Form("title;%s;mixed-evt. pairs", sCorrAxes[i].Data()));
    hists->Add(hsVar[i]);
    hists->Add(hmVar[i]);
    if(i==kZvtx) {
      hzr = (TH1D*) hsVar[i]->Clone();
      hzr->SetFillColor(kYellow);
      hzr->Divide(hmVar[i]);
      hzr->SetTitle("track count ratio;z-vertex (cm);same/mixed evt track multiplicity");
      hists->Add(hzr);
    }
  }
  hists->Add(sig);

  /*
  for (int i=0; i<nEvtAxes; i++) {
    evtVar[i] = hEvt->Projection(i);
    evtVar[i]->SetFillColor(kGray);
    evtVar[i]->SetTitle(Form("title;%s;events", sEvtAxes[i].Data()));
    hists->Add(evtVar[i]);
    if (i==1) cout << evtVar[i]->Integral(1, 89) << endl;
  }
  for (int i=0; i<nTrkAxes; i++) {
    trkVar[i] = hTrk->Projection(i);
    trkVar[i]->SetTitle(Form("title;%s;tracks", sTrkAxes[i].Data()));
    hists->Add(trkVar[i]);
  }

  //hists->Add(hcz);
  */

  for (int n=0; n<hists->GetEntries(); n++) {
    TCanvas* c = new TCanvas(Form("c%d",n), Form("c%d",n), 
		       10*n, 10*n, 700, 500);
    TObject* obj = hists->At(n);
    if (obj->InheritsFrom("TH2")) {
      TH2* h2 = (TH2*) obj;
      h2->Draw("surf1");
    }
    else if (obj->InheritsFrom("TH1")) {
      TH1* h = (TH1*) obj;
      h->SetLineWidth(2);
      h->SetFillStyle(1001); // solid
      if(h->GetFillColor()==kNone) {
	h->SetFillColor(kGreen-9);
      }
      h->Draw("histf");
      h->Draw("histepsame");
    }
  }

  return;
}
