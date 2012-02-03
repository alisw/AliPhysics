// Useful on the prompt:
// TH2* h = (TH2*)Cont_cloizide_DhcAna->FindObject("hS482"); h->Draw("surf1"); h->Integral()
// To quickly count events:
// TH2* h = (TH2*)Cont_cloizide_DhcAna->FindObject("fHEvt");
// h->Draw("colz"); h->Integral()
// TH1D* hc = h->ProjectionY("hc"); hc->Draw();                                             
// hc->Integral(1, 90)

bool ispp = 1;

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TMath.h"
#include "TKey.h"
#include "TROOT.h"
#include "TCutG.h"
#include "TFormula.h"
#include "TProfile2D.h"

TFile* inFile  = new TFile("$MYDATA/lhc11a_28nov2011_sum.root", "read"); 
TFile* outFile = 0;
TString listName = "Cont_cloizide_DhcAna";

// histograms
const int maxHists  = 99999;
TH2 *sHists[maxHists];
TH2 *mHists[maxHists];
TH2 *cHists[maxHists];
TH2 *hEvt  = 0;
TH2 *hTrk  = 0;
TH1 *hDefi = 0; //fHPtTrg
TH1 *hDefj = 0; //fHPtAss
TH1 *hDefk = 0; //fHCent
TH1 *hDefz = 0; //fHZvtx
TFormula *fIndex = 0;

// The angular regions defining the 1d projections. The arrays
// k{Phi,Eta}{Min,Max}[] define these regions.
const int kNRegions = 5;
const char* regionStr[] = 
  {"NSJET", "RIDGE", "ALL", "ETA_NS", "ETA_AS"};
enum eRegion {NSJET, RIDGE, ALL, ETA_NS, ETA_AS};

// Correlation types:
// 0. Same 
// 1. Mixed 
// 2. CF def "A" (proj 2Ds to dphi then divide)
// 3. CF def "B" (divide 2D's, then project to dhpi). The def. B CFs
// already exist in the input file, just pass them through to the new
// output file, renaming "c" --> "cB".
const int nCorrTypes = 4;
enum eCorrType {kS=0, kM=1, kCA=2, kCB=3};
TString sCorrType[] = {"s","m","c","c"};

// There are two sets of centrality binning: the cbh bins are the
// original ones from the task. The cb1,2 bins are a superset
// including cbh1,2 + any desired combinations. In Loop 2 (the
// cent. combination loop), cb1,2 are indexed as "k", while cbh1,2 are
// indexed as "cb". If no combined bins are needed, Loop 2 just copies
// the existing bins to the output file.

// Binning for gsi03+ generation (TH2 train output)
int nCentBins = 19; // 10 orig + 8 new combined bins
double cb1[] =  {0, 0, 0,2,2, 1,3,0,    0, 1, 2, 3, 4,  5, 10, 20, 30, 40, 60 };
double cb2[] =  {10,20,2,5,10,3,5,5,    1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 90 };

// Just for reference, this should be equivalent to cbh1,2[] in gsi03:
// double cb1[] =  { 0, 1, 2, 3, 4,  5, 10, 20, 30, 40 };
// double cb2[] =  { 1, 2, 3, 4, 5, 10, 20, 30, 40, 50 };

int nZbins = 0;        // zvtx binning from task - read from input file
double z1[99]   = { 0 };
double z2[99]   = { 0 };
int nCentBinsHist = 0; // The exact binning from task - read from input file
double cbh1[99] =  { 0 };
double cbh2[99] =  { 0 };

int colors[] = {kYellow+2, kOrange+2, kBlack, kMagenta, kRed, kOrange, 
		kGreen+2, kCyan, kBlue, kViolet, kGray};

// TGraphs to contain complete bin information. Convention: i,j,k,z
// always means pt_trig, pt_assoc, centrality, and z_vtx bin index,
// respectively.
TGraphErrors *gi, *gj, *gk, *gz;
TObjArray* hists = new TObjArray();
TObjArray* profs = new TObjArray(); // profiles storing <pt> info

// Histos used to compute weight factors for summation.
TH1D* hz[101]; // z-vertex dist. at each cent bin
TH1D* hcent = 0; // No longer used - can remove.

double kPhiMin[kNRegions];
double kPhiMax[kNRegions];
double kEtaMin[kNRegions];
double kEtaMax[kNRegions];

// Reassign as fn. arguments
Double_t dEtaOuter = 1.799;
Double_t dEtaInner = 0.8;

// Functions ---------------------------------------------------------
const char* Region(int r);
char* BinLabel(int ijk, int i);
TGraphErrors* BinGraph(int ijk);
void PrintBins();
TObjArray* GetHistList(TFile& fromfile, TString clname, TString dir);
void MakeWeightHistos();
TH1* ProjectTH2(TH2& h, TString xy, float x1, float y1, float x2, float y2,
	   TString name, TString opt);
void LoadHistos();
TH2 *GetHist(int i, int j, int z, int k, char type);

// Implementation ----------------------------------------------------
void MakeWeightHistos()
{
  // make a z-vertex distribution for each centrality bin
  for (int cb=0; cb<nCentBinsHist; cb++) {
    TAxis* ax = hEvt->GetYaxis();
    int bin1 = ax->FindBin(cbh1[cb]);
    int bin2 = ax->FindBin(cbh2[cb]-0.001);

    hz[cb] = hEvt->ProjectionX(Form("hz%d", cb), bin1, bin2);

    if (0)
      hz[cb]->Draw(cb==0? "" : "same");

    // now filled from file
    // hcent->SetBinContent(cb+1, hz[cb]->Integral());
    // hcent->GetXaxis()->SetBinLabel(cb+1,Form("%.2g-%.2g", cbh1[cb],cbh2[cb]));
  }

  //gk->Draw("aep");
  return;
}

void LoadHistos()
{
  if (!inFile) {
    Error("LoadHistos()","TFile ptr zero");
    gSystem->Exit(123);
  }

  TList *list = 0;
  if (ispp) 
    list = (TList*)inFile->Get("Cont_loizides_DhcAna");
  else
    list = (TList*)inFile->Get("Cont_cloizide_DhcAna");
  if (!list) {
    Error("LoadHistos()","TList ptr zero");
    gSystem->Exit(123);
  }

  hEvt = (TH2*)list->FindObject("fHEvt");
  hTrk = (TH2*)list->FindObject("fHTrk");
  if (!hEvt)
    Error("LoadHistos()", "!hEvt");
  if (!hTrk)
    Error("LoadHistos()", "!hTrk");

  hDefi   = (TH1*)list->FindObject("fHPtTrg");
  hDefj   = (TH1*)list->FindObject("fHPtAss");
  hDefk   = (TH1*)list->FindObject("fHCent");
  hDefz   = (TH1*)list->FindObject("fHZvtx");

  gi = new TGraphErrors();
  gj = new TGraphErrors();
  gk = new TGraphErrors();
  gz = new TGraphErrors();
  gi->SetName("TrigPtBins");
  gj->SetName("AsscPtBins");
  gk->SetName("EvCentBins");
  gz->SetName("EvZvtxBins");

  for (int i=0; i<hDefi->GetNbinsX(); i++) {
    gi->SetPoint(i, hDefi->GetBinCenter(i+1), 1.0);
    gi->SetPointError(i, 0.5*hDefi->GetBinWidth(i+1), 0);
  }
  for (int i=0; i<hDefj->GetNbinsX(); i++) {
    gj->SetPoint(i, hDefj->GetBinCenter(i+1), 1.0);
    gj->SetPointError(i, 0.5*hDefj->GetBinWidth(i+1), 0);
  }
  nCentBinsHist = hDefk->GetNbinsX();
  for (int i=0; i<hDefk->GetNbinsX(); i++) {
    cbh1[i] = hDefk->GetBinLowEdge(i+1);
    cbh2[i] = hDefk->GetBinLowEdge(i+2);
  }

  // Get TProfiles with <pt> and <pt^2> info into memory
  for (int n=0; n<list->GetEntries(); n++) {
    TObject* obj = list->At(n);
    TString clName = obj->ClassName();
    if (clName.Contains("TProfile")) {
      profs->Add(obj);
    }
  }

  //  hcent = new TH1D("hcent", "hcent", nCentBinsHist,0,nCentBinsHist);

  for (int i=0; i<nCentBins; i++) {
    gk->SetPoint(i, (cb1[i]+cb2[i])/2., 1.0);
    gk->SetPointError(i, TMath::Abs(cb2[i]-cb1[i])/2, 0);
    //cout << i << " " <<  (cb1[i]+cb2[i])/2. << endl;
  }
         
  nZbins = hDefz->GetNbinsX();
  for (int i=0; i<hDefz->GetNbinsX(); i++) {
    gz->SetPoint(i, hDefz->GetBinCenter(i+1), 1.0);
    gz->SetPointError(i, 0.5*hDefz->GetBinWidth(i+1), 0);
    z1[i] = hDefz->GetBinLowEdge(i+1);
    z2[i] = hDefz->GetBinLowEdge(i+2);
  }

  fIndex = new TFormula("GlobIndex", //this is now relative to 0 (and not to 1 as for histos)
			"(t)*[0]*[1]*[2]+(z)*[0]*[1]+(x)*[0]+(y)+0*[4]");
  fIndex->SetParameters(gi->GetN(),
			gj->GetN(),
			gz->GetN(),
			gk->GetN());
  fIndex->SetParNames("NTrigBins (i)","NAssocBins (j)", "NZvertexBins (z)", "NCentBins (k)"); 

  for (int i=0;i<maxHists;++i) {
    mHists[i]=0;
    sHists[i]=0;
  }

  Int_t ent=list->GetEntries();
  for (Int_t i=0;i<ent;++i) {
    TH2 *obj = dynamic_cast<TH2*>(list->At(i));
    if (!obj)
      continue;
    TString name(obj->GetName());
    if (!name.BeginsWith("hM") && !name.BeginsWith("hS"))
      continue;
    const char *ptr = name.Data()+2;
    Int_t num = atoi(ptr);
    if (num>=maxHists) {
      Error("LoadHistos()", "Found object with too high index: %s %d",name.Data(),maxHists);
      gSystem->Exit(123);
    }
    //obj->SetDirectory(0);
    if (name[1]=='S')
      sHists[num]=obj;
    else if (name[1]=='M')
      mHists[num]=obj;
    else {
      Error("LoadHistos()", "Found object with unknown: %s %d",name.Data(),maxHists);
      gSystem->Exit(123);
    }
  }

  for (int i=0;i<maxHists;++i) {
    if (mHists[i]==0)
      continue;
    if (sHists[i]==0) {
      Error("LoadHistos()", "Both histogram ptr should be set: %d",i);
      continue;
    }

    int nEmptyBinsS = 0;
    int nEmptyBinsM = 0;
    for (int nx=1; nx<36; nx++) {
      for (int ny=1; ny<20; ny++) {
	if(sHists[i]->GetBinContent(nx,ny)==0)
	  nEmptyBinsS++;
	if(mHists[i]->GetBinContent(nx,ny)==0)
	  nEmptyBinsM++;
      }
    }
    double emptyFracS = double(nEmptyBinsS)/720.;
    double emptyFracM = double(nEmptyBinsM)/720.;
    if (emptyFracS > 0.50)
      Warning("LoadHistos","%d %s: %.0f%% of bins are EMPTY", i, sHists[i]->GetTitle(), 100*emptyFracS);
    if (emptyFracM > 0.50)
      Warning("LoadHistos","%d %s: %.0f%% of bins are EMPTY", i, mHists[i]->GetTitle(), 100*emptyFracM);

    cHists[i] = (TH2F*)sHists[i]->Clone(Form("hC%d",i));
    cHists[i]->Divide(sHists[i], mHists[i], 1./sHists[i]->Integral(), 1./mHists[i]->Integral());
  }
}

TH2 *GetHist(int i, int j, int z, int k, char type)
{
  if (!fIndex)
    return 0;

  Int_t ind = fIndex->Eval(i,j,z,k);
  TH2* h = 0;

  if ((ind>=maxHists)||(ind<0)) {
    Error("GetHist", "Input %d %d %d %d gives to large index %d",i,j,z,k,ind);
    return 0;
  }

  if (type == 's') 
    h = sHists[ind];
  else if (type == 'm')
    h = mHists[ind];
  else if (type == 'c')
    h = cHists[ind];
  else {
    Error("GetHist", "Unknown type %c",type);
    return 0;
  }

  if (h)
    if (h->GetSumw2N() == 0)
      h->Sumw2();

  return h;
}

void SumAndProjectTH2s(Double_t dEtaInnerArg = 1.2, Double_t dEtaOuterArg = 1.799)
{
  if (ispp) {
    nCentBins = 1;
    cb1[0] = 0;
    cb2[0] = 0;
  }
  
  dEtaInner = dEtaInnerArg; 
  dEtaOuter = dEtaOuterArg;

  outFile = new TFile(Form("$MYDATA/lhc11a_etamin%02d.root", (int)(10*dEtaInner)), "recreate");

  LoadHistos();
  MakeWeightHistos();
  PrintBins();

  int nk = ispp ? 1 : 10;
  for (int k=0; k<nk; k++) {
    TH2* h = GetHist(0,0,4, k, 's');
    cout << Form("%s %.2g",h->GetTitle(), h->Integral()) << endl;
  }

  // Dude, this is such a hack. TODO: calculate safely the number of y
  // bins included in an x projection, or vice versa.
  double nProjectedBins[] = {8, 12, 18, 18, 18};

  kPhiMin[NSJET] = -TMath::PiOver2();
  kPhiMax[NSJET] = 3*TMath::PiOver2();
  kEtaMin[NSJET] = -dEtaInner;
  kEtaMax[NSJET] = +dEtaInner;
  kPhiMin[RIDGE] = -TMath::PiOver2();
  kPhiMax[RIDGE] = 3*TMath::PiOver2();
  kEtaMin[RIDGE] = -dEtaInner; // selection gets inverted
  kEtaMax[RIDGE] = +dEtaInner;
  kPhiMin[ALL]   = -TMath::PiOver2();
  kPhiMax[ALL]   = 3*TMath::PiOver2();
  kEtaMin[ALL]   = -dEtaOuter;
  kEtaMax[ALL]   = +dEtaOuter;

  kPhiMin[ETA_NS]   = -TMath::PiOver2();
  kPhiMax[ETA_NS]   = +TMath::PiOver2();
  kEtaMin[ETA_NS]   = -dEtaOuter;
  kEtaMax[ETA_NS]   = +dEtaOuter;

  kPhiMin[ETA_AS]   =   TMath::PiOver2();
  kPhiMax[ETA_AS]   = 3*TMath::PiOver2();
  kEtaMin[ETA_AS]   = -dEtaOuter;
  kEtaMax[ETA_AS]   = +dEtaOuter;


  // First loop: sum z-vertex bins - 2D
  for (int cb=0; cb<nCentBinsHist; cb++) {

    double z_int = hz[cb]->Integral();
    double centLo = cbh1[cb];
    double centHi = cbh2[cb];
    if (ispp) {
      centLo = -1;
      centHi = 101;
    }
    cout << Form("Cent bin %d/%d: %.2g to %.2g%%", cb, nCentBinsHist, centLo, centHi)
	 << endl;

    for (int i=0; i<gi->GetN(); i++) {
      for (int j=0; j<=i; j++) {

	for (int ict=0; ict<nCorrTypes-1; ict++) {

	  TH2* hz2 = 0;
	  double wtSum = 0;
	  bool printZWeights = (0 && i==0 && j==0 && ict==kCA);

	  if (printZWeights) 
	    cout << Form("i%d j%d %.0fto%.0f: ", i, j, centLo, centHi);

	  for (int iz=0; iz<nZbins; iz++) {

	    int zbin1 = hz[cb]->FindBin(z1[iz]);
	    int zbin2 = hz[cb]->FindBin(z2[iz]-0.001);

	    double weight = hz[cb]->Integral(zbin1, zbin2) / z_int;

	    if (ict==kS || ict==kM) 
	      weight = 1.0; // Only correlation functions should be weighted
	    
	    wtSum += weight;

	    // There was an off-by-1 bug here!! Getting iz offset from
	    // zero, but loop was from iz=1 up. Fixed.
	    TH2* hc2 = GetHist(i,j,iz,cb,sCorrType[ict][0]);
	    if (!hc2) {
	      if (iz==0 && i==0 && j==0)
		Info("Loop 1", "Histo not found - skipping");
	      continue;
	    }
	    // if (i==7 && j==6 && ict < kCA) {
	    //   cout << Form("%s ===>  %.3g", hc2->GetTitle(), hc2->Integral()) << endl;
	    // }
	    if (printZWeights) {
	      cout << Form("%.2gto%.2g:%.2g ",z1[iz], z2[iz], weight);
	    }
	    
	    int cLo = ispp ? 0 : centLo;
	    int cHi = ispp ? 0 : centHi;
	    const char* newName = Form("%s_%s_%d_%d_%dto%d", 
				       "ETAPHI", sCorrType[ict].Data(),
				       i, j, cLo, cHi);
	    
	    if (iz==0) {
	      hz2 = (TH2*)hc2->Clone(newName);
	      hz2->Reset();
	    }
	    hz2->Add(hc2, weight);
	  } // iz loop
	  
	  if (printZWeights)  
	    cout << endl;
	  
	  hists->Add(hz2);
	  if((ict==kCA || ict==kCB) && TMath::Abs(wtSum - 1.) > 0.01)
	    Warning("MakeWeightHistos()", "weight sum %.2g", wtSum);
	  if (ict == kCA) {
	    double normInt = hz2->Integral()/720;
	    if (TMath::Abs(1-normInt) > 0.05)
	      Warning("Loop 1", "%s (%d) integral not 1: %.3g", 
		      hz2->GetName(), (int)fIndex->Eval(i,j,3,cb), normInt);
	    
	  }
	  
	}
      }
    }
  }
  
  cout<<"\nFinished z-vertex bin sums\n"<<endl;

  /*
  for (int n=0; n<hists->GetEntries(); n++) {
    TH2* h = (TH2*)hists->At(n);
    TString name(h->GetName());
    if (name.Contains("s_0_0"))
      cout<<Form("%s %s %.2g",name.Data(), h->GetTitle(), h->Integral())<<endl;
  }
  */

  outFile->mkdir("PbPb");
  outFile->cd("PbPb");

  // Second loop: sum centrality bins - 2D histos. Purpose: make the
  // desired centrality bin k in the TGraph from smaller centrality
  // bins cb that already exist.
  for (int i=0; i<gi->GetN(); i++) {
    for (int j=0; j<=i; j++) {
      for (int k=0; k<gk->GetN(); k++) { // the target centrality bin
	for (int ict=0; ict<nCorrTypes-1; ict++) {



	  TObjArray* arr = 0;
	  std::vector<double> weights;
	  std::vector<TString> centIntervals;
	  double wtsum = 0;
	  // double centLo = cb1[k];
	  // double centHi = cb2[k]; same as below
	  double centLo = ispp? 0 : gk->GetX()[k] - gk->GetEX()[k];
	  double centHi = ispp? 0 : gk->GetX()[k] + gk->GetEX()[k];
	  const char* name = Form("%s_%s_%d_%d_%.0fto%.0f", 
				  "ETAPHI", sCorrType[ict].Data(), i, j, centLo,centHi);
	  const char* newName = Form("%s_%s_%d_%d_%d", "ETAPHI", sCorrType[ict].Data(), i, j, k);
	  const char* title   = Form("p_{T}^{t} %s, p_{T}^{a} %s, %.0f-%.0f%%", 
				     BinLabel(0, i), 
				     BinLabel(1, j), 
				     centLo, centHi);

	  //	  Printf("%g %g", centLo, centHi);

	  TH2* h2 = 0;
	  h2 = (TH2*) hists->FindObject(name);

	  if (h2) {
	    // Give correct overall normalization for CFs.
	    if (ict == kCA || ict == kCB)
	      h2->Scale(h2->GetNbinsX()*h2->GetNbinsY()/h2->Integral());
	   
	    h2->SetNameTitle(newName, title);
	    h2->Write(newName);

	    if (1) cout<<h2->GetName()<<" "<<endl;
	  }
	  else { // If cent. bin k doesn't exist, it needs to be created as a sum.
	    for (int cb=0; cb<nCentBinsHist; cb++) { // pre-existing cent bins
	      const char* name2 = Form("%s_%s_%d_%d_%.0fto%.0f", 
				       "ETAPHI", sCorrType[ict].Data(), i, j, cbh1[cb], cbh2[cb]);
	      TH2* hcb2 = (TH2*) hists->FindObject(name2);
	      //	      TH2* h = GetHist(0,0,4, k, 's');
	      if (!hcb2) {
		Error("Loop 2","!hcb2");
		continue;
	      }
	      cout << name2 <<" "<< cbh1[cb] << " " << centLo << " " << cbh2[cb] << " " << centHi << endl;

	      // Store relevant histos and their weights
	      if ( cbh1[cb] >= centLo && cbh2[cb] <= centHi ) {
		if (1) cout<<"Adding " << hcb2->GetName()<<" "<<endl;
		if (!arr) arr = new TObjArray();
		arr->Add(hcb2);

		// weight by integral of same-event 2D distribution rather than event count
		// small but visible correction from old way. 
		/* old way:
		weights.push_back(hcent->GetBinContent(cb+1));
		wtsum += hcent->GetBinContent(cb+1);
		*/
		// new way
		const char* ss = Form("ETAPHI_s_%d_%d_%.0fto%.0f", i, j, cbh1[cb], cbh2[cb]);
		TH2* hs = (TH2*) hists->FindObject(ss);
		if (!hs)
		  Error("Loop 2","Problem finding %s", ss);
		double wt = hs->Integral();
		weights.push_back(wt);
		centIntervals.push_back(TString(Form("%.0fto%.0f", cbh1[cb], cbh2[cb])));
		wtsum += wt;

	      }
	    } // cb loop

	    // Add up histos in arr. This part of code only gets used if
	    // centrality combination is requested by including new
	    // bins in cb1,2[] arrays.

	    bool printWeights = (ict == kCA); // print once/bin, not 3x
	    double tot=0;

	    if (printWeights) cout<<name<<" "<<flush;
	    
  	  cout<< arr <<endl;
	  if (!arr)
	    continue;

	    for (int n=0; n<arr->GetEntries(); n++) {

	      tot += weights.at(n)/wtsum;
	      TH2* hn = (TH2*)arr->At(n);

	      if (n==0) {
		h2 = (TH2*)hn->Clone(newName);
		h2->Scale(weights.at(n)/wtsum);
	      }
	      else {
		h2->Add(hn, weights.at(n)/wtsum);
	      }
	      if (printWeights) 
		cout << centIntervals.at(n).Data() << ":" 
		     << weights.at(n)/wtsum << " " << flush;
	    }
	    if (printWeights) cout << " = " << tot << endl;

	    // Give correct overall normalization for CFs.
	    if (ict == kCA || ict == kCB)
	      h2->Scale(h2->GetNbinsX()*h2->GetNbinsY()/h2->Integral());

	    h2->SetNameTitle(newName, title);
	    h2->Write(newName);

	    if (arr) 
	      delete arr; // arr doesn't own its histos, so ok to delete.
	  }
	} // ict
      } // k
    } // j
  } // i

  cout<<"\nFinished centrality bin sums\n"<<endl;

  outFile->cd();
  for (int ijkz=0; ijkz<4; ijkz++)
    BinGraph(ijkz)->Write();
  profs->Write();

  TObjArray* summed = GetHistList(*outFile, "TH2", "PbPb");
  TObjArray* dphis = new TObjArray();

  outFile->cd("PbPb");
  for (int n=0; n<summed->GetEntries(); n++) {
    TH2* h = (TH2*)summed->At(n);
    for (int r=0; r<5; r++) {
      TString name1D(h->GetName());
      name1D.ReplaceAll("ETAPHI", Region(r));
      TString projStr = "x";
      if (r==RIDGE)
	projStr.Prepend("-");
      else if (r==ETA_NS || r==ETA_AS)
	projStr = "y";

      TH1* hp = ProjectTH2(*h, projStr, 
      			   kPhiMin[r], kEtaMin[r], kPhiMax[r], kEtaMax[r], 
      			   name1D, "e");
      hp->SetTitle(h->GetTitle());

      // TODO check this
      if (name1D.Contains("_c_"))
	hp->Scale(1./nProjectedBins[r]);

      if (r==NSJET || r==RIDGE || r==ALL)
	dphis->Add(hp);

      hp->Write();
    }
  }

  // CF defs A and B
  // First find existing "Def. B" CF _c_. Then create "Def. A" CFs as well.
  for (int n=0; n<dphis->GetEntries(); n++) {
    TObject* obj = dphis->At(n);
    TString hname =  obj->GetName();
    TH1* hcB = (TH1*)obj;
    if (hname.Contains("_c_") ) {
      hcB = (TH1*)obj;
    }
    else
      continue;

    TString s(hname);
    TString m(hname);
    TString cA(hname);
    TString cB(hname);
    s.ReplaceAll("_c_", "_s_");
    m.ReplaceAll("_c_", "_m_");
    cA.ReplaceAll("_c_", "_cA_");
    cB.ReplaceAll("_c_", "_cB_");
    
    TH1* hsame  = (TH1*)dphis->FindObject(s.Data());
    TH1* hmixed = (TH1*)dphis->FindObject(m.Data());
    TH1* hcA = (TH1*)hcB->Clone(cA.Data());

    if (0)
      Info("", "%s %s %s", hsame->GetName(), hmixed->GetName(), hcA->GetName());

    hcA->Divide(hsame, hmixed, 1./hsame->Integral(), 1./hmixed->Integral());
    hcB->SetName(cB.Data());

    hcA->Write();
    hcB->Write();

  }
  
  outFile->cd("PbPb");
  int nObjectsInFile = gDirectory->GetListOfKeys()->GetEntries();
  cout <<"\nWrote " << nObjectsInFile 
       << " objects to PbPb/ in " << outFile->GetName() << endl;

  cout << "Closing file..." << endl; 
  outFile->Close();

  return;
}

const char* Region(int r)
{
  return regionStr[r];
}

char* BinLabel(int ijkz, int i)
{
  // Returns bin extents like "a-b" to <= 4 sig. fig. precision
  TGraphErrors* g = BinGraph(ijkz);
  if (i<0 || i>=g->GetN()) {
    Error("bin_label((TGraphErrors&, Int_t)", "Error: no bin %d", i);
    return 0;
  }

  double x1 = g->GetX()[i] - g->GetEX()[i];
  double x2 = g->GetX()[i] + g->GetEX()[i];

  return Form("%.*g-%.*g", 4, x1, 4, x2);
}

TGraphErrors* BinGraph(int ijkz)
{
  if (ijkz < 0 || ijkz > 3)
    return 0;
  else if (ijkz==0) return gi;
  else if (ijkz==1) return gj;
  else if (ijkz==2) return gk;
  else if (ijkz==3) return gz;
  else 
    return 0;
}

TObjArray* GetHistList(TFile& file, TString clname, TString dir)
{
  file.cd(dir.Data());

  TObjArray* hList = new TObjArray();
  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  
  while ((key=(TKey*)next())) {
    TString className(key->GetClassName());
    TString keyName(key->GetName());
    if (0) 
      printf("%10s %20s\n", className.Data(), keyName.Data());
    
    if (className.Contains(clname) && clname.Contains("TH2")) {
      hList->Add((TH2*)gDirectory->Get(keyName.Data()));
    }
  }

  cout << hList->GetEntries() << " objects retrieved from "
       << file.GetName()  << "/" << gDirectory->GetName() 
       << endl;

  return hList;
}

TH1* 
ProjectTH2(TH2& h, TString xy, float x1, float y1, float x2, float y2,
	   TString name, TString opt)
{
  TH1* hp = 0;
  double x[] = {x1,x2,x2,x1,x1};
  double y[] = {y1,y1,y2,y2,y1};
  TCutG *cutg = new TCutG("cutg", 5, x, y);
  double area = cutg->Area();
  TString projOpt = "[cutg]";

  // Invert selection: fill bins outside rectangle
  if (xy.Contains("-"))
    projOpt.ReplaceAll("[cutg]", "[-cutg]");
  
  // "e" errors, "d" draw in current pad, "o" original full axis
  projOpt.Append(opt);

  if (xy.Contains("x")) {
    hp = h.ProjectionX(name.Data(), 0, -1, projOpt.Data());
  }
  else if (xy.Contains("y")) {
    hp = h.ProjectionY(name.Data(), 0, -1, projOpt.Data());
  }
  if (!hp)
    Error("ProjectTH2()","!hp");


  if (0) { // For test/debug purposes
    int bx1, bx2, by1, by2;
    double wx, hy;
    TAxis *ax = h.GetXaxis(), *ay = h.GetYaxis();

    bx1 = ax->FindBin(x1);
    bx2 = ax->FindBin(x2);
    by1 = ay->FindBin(y1);
    by2 = ay->FindBin(y2);
    wx = ax->GetBinUpEdge(bx2) - ax->GetBinLowEdge(bx1);
    hy = ay->GetBinUpEdge(by2) - ay->GetBinLowEdge(by1);
    
    if (xy.Contains("-"))
      cout << Form("Selection is inverted.") << endl;
    
    cout << Form("Selection area: %f", area) << endl;
    cout << Form("Actual projected area %f", wx*hy) << endl;
  }

  return hp;
}

void PrintBins()
{
  TGraphErrors* g = 0;

  cout << endl;
  cout <<"-----------------Bin Information-----------------"<< endl;
  for (int n=0; n<4; n++) {
    g = BinGraph(n);
    cout << g->GetName() << ": ";
    for (int i=0; i<g->GetN(); i++) {
      cout << Form("%s ", BinLabel(n, i));
    }
    cout << endl;
  }
  cout <<"-------------------------------------------------"<< endl;
 return;
}
