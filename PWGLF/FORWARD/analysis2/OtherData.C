#ifndef __CINT__
# include <TFile.h>
# include <TString.h>
# include <TGraphAsymmErrors.h>
# include <TMultiGraph.h>
# include <TROOT.h>
# include <TError.h>
# include <TSystem.h>
# include <TH1F.h>
# include <TStyle.h>
# include <TMath.h>
#else 
class TFile;
class TGraphAsymmErrors;
class TMultiGraph;
class TGraph;
class TH1;
#endif 


struct RefData
{
  //____________________________________________________________________
  /**
   * Values used 
   * 
   * @ingroup pwglf_forward_otherdata 
   */
  enum { 
    UA5, 
    CMS, 
    ALICE, 
    WIP,
    PYTHIA,
    INEL, 
    INELGt0, 
    NSD
  };
  
  //____________________________________________________________________
  /** 
   * Get a pointer to our data file.  @a path is the path to the file
   * containing the data.  If it is null, then first search the
   * current directory, and if not found there, search in specific
   * AliROOT directory. If @a rw is try the file is (re-)opened in
   * UPDATE mode.
   * 
   * @param path Path to file.
   * @param rw   If true, open read/write.  
   * 
   * @return Pointer to file, or null
   */
  static TFile* GetFile(const char* path=0, Bool_t rw=false) { 
    TString base = ((path && path[0] != '\0') ? 
		    gSystem->BaseName(path) : "other.root");
      
    TObject* o = gROOT->GetListOfFiles()->FindObject(base);
    
    TFile* f  = 0;
    if (o) {
      f = static_cast<TFile*>(o);
      if (!rw) return f;
      if (!f->IsWritable() && f->ReOpen("UPDATE") < 0) return 0;
      return f;
    }

    const char* mode = (rw ? "UPDATE" : "READ");
    if (path && path[0] != '\0' && !gSystem->AccessPathName(path)) 
      f = TFile::Open(path, mode);
    if (!f && !gSystem->AccessPathName("other.root")) 
      f = TFile::Open("other.root", mode);
    if (!f) 
      f = TFile::Open("$ALICE_ROOT/PWGLF/FORWARD/analysis2/other.root",mode);
    if (!f) 
      ::Error("", "Failed to open file");
    return f;
  }
  //____________________________________________________________________
  /** 
   * Get the collision system name from the identifier (1: pp, 2:
   * PbPb, 3: pPb)
   * 
   * @param sys Collision system identifier 
   * 
   * @return String or null
   */
  static const char* SysName(UShort_t sys) { 
    switch (sys) { 
    case 1: return "pp";
    case 2: return "PbPb";
    case 3: return "pPb";
    }
    ::Error("", "Unknown system: %d", sys);
    return 0;
  }
  //____________________________________________________________________
  /** 
   * Get the zero-padded collision energy name 
   * 
   * @param sNN Collision energy (in GeV)
   * 
   * @return Zero-padded collision energy
   */
  static const char* SNNName(UShort_t sNN) 
  { 
    return Form("%05d", sNN);
  }
  //____________________________________________________________________
  /** 
   * Get the centrality method prefix.  @a trg is a bit mask, of which
   * bits 4-7 are used here.  The flags are 
   *
   * - 0x10 V0M 
   * - 0x20 V0A
   * - 0x40 ZNA
   * - 0X80 ZNC
   * 
   * @param trg Bits
   * 
   * @return Prefix string or empty string 
   */
  static const char* CntName(UShort_t trg) 
  {
    switch (trg >> 4) { 
    case 1: return "V0M_";
    case 2: return "V0A_";
    case 4: return "ZNA_";
    case 8: return "ZNC_";
    }
    return "";
  }
  //____________________________________________________________________
  /** 
   * Get the trigger name.  If @f$ c_2 > c_1@f$ then the bits of @a
   * trg are ignored, save for the bits 4-7 which are interpreted
   * according to CntName.
   *
   * The meaning of the bits are 
   *
   * - 0x01 INEL
   * - 0x02 INEL>0
   * - 0x04 NSD 
   * - 0xf0 Mask for centrality estimator. 
   * 
   * @param trg Trigger mask.  
   * @param c1  Least centrality @f$ c_1@f$ 
   * @param c2  Largest centrality @f$ c_2@f$
   * 
   * @return Trigger string or null
   */
  static const char* TrgName(UShort_t trg, UShort_t c1, UShort_t c2)
  {
    if (c2 > c1) return Form("%s%03d_%03d", CntName(trg), c1, c2);

    switch (trg) { 
    case 1: return "INEL";
    case 2: return "INELGT0";
    case 4: return "NSD";
    }
    ::Error("", "Unknown trigger: %d", trg);
    return 0;
  }
  //____________________________________________________________________
  /** 
   * Get the experiment name.  
   *
   * - 0: UA5 
   * - 1: CMS 
   * - 2: ALICE (published or pre-print)
   * - 3: ALICE Work-in-progress 
   * - 4: PYTHIA (or MC)
   * 
   * @param which Experiment identifier 
   * 
   * @return Experiment name or null
   */
  static const char* ExpName(UShort_t which) 
  {
    switch (which) { 
    case 0: return "UA5";
    case 1: return "CMS";
    case 2: return "ALICE";
    case 3: return "WIP";
    case 4: return "PYTHIA";
    }
    ::Error("", "Unknown experiment: %d", which); 
    return 0;
  }
  // _________________________________________________________________
  /** 
   * Get graphs for selected experiments under @a d.  
   * 
   * @param d        Directory to search
   * @param which    Which experiments to get data from 
   * @param verbose  Whether to be verbose or not 
   * 
   * @return Graph of data, or null
   */
  static TMultiGraph* GetExps(TDirectory* d, UShort_t which, Bool_t verbose)
  {
    TMultiGraph* ret = 0;
    for (UShort_t w = UA5; w <= PYTHIA; w++) { 
      if (!(which & (1 << w))) continue;

      const char* expName = ExpName(w);
      TDirectory* expDir  = 0;
      if (!expName || !(expDir = d->GetDirectory(expName))) continue;

      TObject* o = expDir->Get("data");
      if (o) {
	if (!ret) ret = new TMultiGraph();
	TMultiGraph* mg = static_cast<TMultiGraph*>(o);
	if (w == WIP) { 
	  TIter   next(mg->GetListOfGraphs());
	  TGraph* g = 0;
	  while ((g = static_cast<TGraph*>(next()))) 
	    if (g->GetMarkerColor() == kMagenta+3) {
	      g->SetMarkerColor(kCyan-6);
	      g->SetLineColor(kCyan-6);
	    }
	}
	ret->Add(mg);
      }
    }
    if (!ret && verbose)
      ::Error("GetExps", "Didn't get any data for exp=0x%x in dir %s", 
	      which, d->GetPath());
    return ret;
  }

  // _________________________________________________________________
  /** 
   * Get data for selected trigger and experiments under @a d.  
   * 
   * @param d        Directory to seach 
   * @param type     Which triggers 
   * @param which    Which experiments
   * @param verbose  Whether to be verbose 
   * 
   * @return Graph of data, or null
   */
  static TMultiGraph* GetTrigs(TDirectory* d, UShort_t type, 
			       UShort_t which, Bool_t verbose) 
  {
    TMultiGraph* ret = 0;
    for (UShort_t t = INEL; t <= NSD; t++) { 
      UShort_t trg = (1 << (t-INEL));
      if (!(type & trg)) {
	if (verbose) 
	  ::Info("GetTrigs", "Skipping trigger 0x%x (0x%x)", trg, type);
	continue;
      }

      const char* trgName = TrgName(trg, 0, 0);
      TDirectory* trgDir  = 0;
      if (!trgName || !(trgDir = d->GetDirectory(trgName))) {
	if (verbose) 
	  ::Warning("GetTrigs", "No directory %s for 0x%x in %s", 
		    trgName, trg, d->GetPath());
	continue;
      }

      TMultiGraph* g = GetExps(trgDir, which, verbose);
      if (g) { 
	if (!ret) ret = new TMultiGraph();
	ret->Add(g);
      }
    }
    if (!ret) 
      ::Error("GetTrigs", 
	      "Didn't get any data for trigger=0x%x and exp=0x%x in dir %s", 
	      type, which, d->GetPath());
    return ret;
  }

  // _________________________________________________________________
  /** 
   * Get data for selected centrality range 
   * 
   * @param d            Directory to search
   * @param experiment   Which experiment 
   * @param trigger      Which centrality estimator (possibly 0)
   * @param centLow      Least centrality 
   * @param centHigh     Largetst centrality
   * @param verbose      Whether to be verbose 
   * 
   * @return Graph of data or null
   */
  static TMultiGraph* GetCents(TDirectory* d, UShort_t experiment, 
			       UShort_t trigger, UShort_t centLow, 
			       UShort_t centHigh, Bool_t verbose) 
  {
    // We need to find the bins we can and check for the
    // experiments
    TMultiGraph* ret = 0;
    TIter    next(d->GetListOfKeys());
    TObject* obj = 0;
    const char* cntPre = CntName(trigger);
    // Info("", "trigger=0x%x pre=%s", trigger, cntPre);
    while ((obj = next())) { 
      TString n(obj->GetName());
      if (n.EqualTo("all")) continue;
      UShort_t off = 0;
      if (cntPre && cntPre[0] != '\0') {
	if (!n.BeginsWith(cntPre, TString::kIgnoreCase)) continue;
	off = 4;
      }
      if (n.Length() != 7+off) continue;
      
      TString l1(n(off,  3)); TString l2 = l1.Strip(TString::kLeading, '0');
      TString h1(n(off+4,3)); TString h2 = h1.Strip(TString::kLeading, '0');
      UShort_t c1 = l2.Atoi();
      UShort_t c2 = h2.Atoi();
      
      // Info("", "n=%s off=%d c1=%d c2=%d", n.Data(), off, c1, c2);
      if (c1 < centLow || c2 > centHigh) {
	if (verbose) ::Info("", "Skipping %s in %s", n.Data(),d->GetPath());
	continue;
      }
      
      TDirectory* centDir = d->GetDirectory(obj->GetName());
      if (!centDir) continue;
      
      TMultiGraph* exps = GetExps(centDir, experiment, verbose);
      if (exps) {
	if (!ret) ret = new TMultiGraph();
	ret->Add(exps);
      }
    } // experiment (key)
    if (!ret && verbose) 
      ::Error("GetCents", "No graphs for centralities %d-%d%% in %s", 
	      centLow, centHigh, d->GetPath());
    return ret;
  }
  //____________________________________________________________________
  /** 
   * Get a multi graph of data for a given energy and trigger type 
   * 
   * @param sys      Collision system (1: pp, 2: PbPb, 3: pPb)
   * @param energy   Energy in GeV (900, 2360, 2760, 7000, 8000)
   * @param triggers Bit pattern of trigger type 
   *   - 0x01 INEL 
   *   - 0x02 INEL>0
   *   - 0x04 NSD 
   *   - 0xF0 Mask for centrality estimator 
   * @param centLow     Low centrality cut (not pp)
   * @param centHigh    High centrality cut (not pp)
   * @param experiments From which experiments 
   * 
   * @return A multi graph with the selected data. 
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TMultiGraph* GetData(UShort_t sys, 
			      UShort_t sNN,
			      UShort_t triggers=0x1, 
			      UShort_t centLow=0, 
			      UShort_t centHigh=0, 
			      UShort_t experiments=0x7)
  {
    Bool_t verbose = false;
    UShort_t trg = (triggers & 0xF7);
    if (triggers & 0x2000) trg |= 0x4;

    TFile* f = GetFile(0,false);
    if (!f) return 0;

    TDirectory* sysDir = 0;
    const char* sysName = SysName(sys);
    if (!sysName || !(sysDir = f->GetDirectory(sysName))) { 
      ::Error("", "Invalid system %d (%s)", sys, sysName);
      return 0;
    }

    TDirectory* sNNDir = 0;
    const char* sNNName = SNNName(sNN);
    if (!sNNName || !(sNNDir = sysDir->GetDirectory(sNNName))) { 
      ::Error("", "Invalid CMS energy %d (%s)", sNN, sNNName);
      return 0;
    }

    TMultiGraph* ret = 0;
    // If we have a centrality request 
    if (centHigh > centLow) { 
      if (centLow == 0 && centHigh >= 100) 
	ret = GetCents(sNNDir, experiments, trg, 
		       centLow, centHigh, verbose);
      else {
	// Look for specific centrality bin 
	TDirectory* centDir = sNNDir->GetDirectory(TrgName(trg, 
							   centLow,centHigh));
	if (!centDir) {
	  Warning("", "No directory '%s' (0x%x,%d%d)", 
		  TrgName(trg, centLow,centHigh), trg, centLow, centHigh);
	  return 0;
	}

	return GetExps(centDir, experiments, verbose);
      }
    } // centHigh > centLow
    else 
      ret = GetTrigs(sNNDir, trg, experiments, verbose);

    if (ret) {
      TString title;
      FormatTitle(title, sys, sNN, trg, centLow, centHigh);
      ret->SetTitle(title);
    }
    return ret;
  }
  //__________________________________________________________________
  /** 
   * Format title of a plot
   * 
   * @param title     On return, the title 
   * @param sys      Collision system (1: pp, 2: PbPb, 3: pPb)
   * @param energy   Energy in GeV (900, 2360, 2760, 7000, 8000)
   * @param triggers Bit pattern of trigger type 
   *   - 0x01 INEL 
   *   - 0x02 INEL>0
   *   - 0x04 NSD 
   *   - 0xF0 Mask for centrality estimator 
   * @param centLow  Low centrality cut (not pp)
   * @param centHigh High centrality cut (not pp)
   * @param seenUA5  If true and sys=1, then put in p-pbar
   */
  static void FormatTitle(TString& title,
			  UShort_t sys, 
			  UShort_t sNN,
			  UShort_t triggers, 
			  UShort_t centLow, 
			  UShort_t centHigh,
			  Bool_t   seenUA5=false)
  {
    TString sn(SysName(sys));
    if (seenUA5) sn.Append("(p#bar{p})");

    TString en(Form("#sqrt{s%s}=", (sys==1 ? "" : "_{NN}")));
    if (sNN < 1000)             en.Append(Form("%dGeV", sNN));
    else if ((sNN % 1000) == 0) en.Append(Form("%dTeV",  (sNN/1000)));
    else                        en.Append(Form("%.2fTeV",  Float_t(sNN)/1000));
    TString tn;
    if (centHigh > centLow) 
      tn = Form("%d%% - %d%% central", centLow, centHigh);
    else { 
      for (UShort_t t = INEL; t <= NSD; t++) { 
	UShort_t trg = (1 << (t-INEL));
	if (!(triggers & trg)) continue;
	if (!tn.IsNull()) tn.Append("|");
	switch (t) { 
	case INEL:    tn.Append("INEL"); break;
	case INELGt0: tn.Append("INEL>0"); break;
	case NSD:     tn.Append("NSD"); break;
	}
      } // for 
    }
    if (!en.IsNull()) en.Prepend(", ");
    if (!tn.IsNull()) tn.Prepend(", ");
    title.Form("%s%s%s", sn.Data(), en.Data(), tn.Data());
  }
  //=== Importing ====================================================
  enum { 
    /** Style used for UA5 data */
    UA5Style   = 21, 
    /** Style used for CMS data */
    CMSStyle   = 29, 
    /** Style used for ALICE published data */
    ALICEStyle = 27,
    /** Color used for ALICE work-in-progress data */
    WIPStyle = 33,
    /** Style used for Pythia data */
    PYTHIAStyle = 28,
    /** Color used for UA5 data */
    UA5Color   = kBlue+1,
    /** Color used for Pytia data */
    PYTHIAColor = kGray+2,
    /** Color used for CMS data */
    CMSColor   = kGreen+1,
    /** Color used for ALICE data */
    ALICEColor = kMagenta+1,
    /** Color used for ALICE work-in-progress data */
    WIPColor = kCyan+2
  }; 
  enum { 
    /** Marker style INEL data */
    INELStyle   = 22,
    /** Marker style INEL>0 data */
    INELGt0Style= 29,
    /** Marker style NSD data */
    NSDStyle    = 23,
    /** Color used for UA5 data */
    INELColor   = kBlue+1,
    /** Color used for CMS data */
    INELGt0Color = kGreen+1,
    /** Color used for ALICE data */
    NSDColor     = kMagenta+1
  };
  enum {
    /** Style offset for mirror data */
    MirrorOff  = 4
  };
  //____________________________________________________________________
  /** 
   * Set graph attributes based on trigger type and experiment. 
   * 
   * @param g        Graph
   * @param exp      Experiment 
   * @param mirror   True if mirrored data 
   * @param name     Name of graph 
   * @param title    Title of graph 
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static void SetGraphAttributes(TGraph* g, 
				 Int_t /*trig*/, 
				 Int_t exp, 
				 bool mirror,
				 const Char_t* name, 
				 const Char_t* title)
  {
    Int_t color = 0;
    Int_t style = 0;
    switch (exp) { 
    case UA5:       color = UA5Color;    style = UA5Style;     break;
    case CMS:       color = CMSColor;    style = CMSStyle;     break;
    case ALICE:     color = ALICEColor;  style = ALICEStyle;   break;
    case WIP:       color = WIPColor;    style = WIPStyle;     break;
    case PYTHIA:    color = PYTHIAColor; style = PYTHIAStyle;  break;
    }
    Float_t size = g->GetMarkerSize();
    switch (style) {
    case 21: // fall-through
    case 25: size *= 0.8; break;
    case 27: size *= 1.4; break;
    case 33: size *= 1.4; break;
    }
    
    if (mirror) style += MirrorOff;

    if (name)  g->SetName(name);
    if (title) g->SetTitle(title);
    g->SetMarkerStyle(style);
    g->SetMarkerSize(size);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetFillColor(0);
    g->SetFillStyle(0);
    g->GetHistogram()->SetStats(kFALSE);
    g->GetHistogram()->SetXTitle("#eta");
    g->GetHistogram()->SetYTitle("#frac{1}{N} #frac{dN_{ch}}{#eta}");
  }
  //__________________________________________________________________
  /** 
   * Get the color for a centrality bin
   * 
   * @param centLow  Centrality bin 
   * @param centHigh Centrality bin 
   * 
   * @return Color 
   */
  static Int_t CentralityColor(UShort_t centLow, 
			       UShort_t centHigh,
			       UShort_t /*nBins*/=0)
  {
#if 0
    if (nBins > 0 && nBins < 6) { 
      switch (bin) { 
      case 1: return kRed+2;
      case 2: return kGreen+2;
      case 3: return kBlue+1;
      case 4: return kCyan+1;
      case 5: return kMagenta+1;
      case 6: return kYellow+2;
      }
    }
#endif
    gStyle->SetPalette(1);
    Float_t  fc       = (centLow+double(centHigh-centLow)/2) / 100;
    Int_t    nCol     = gStyle->GetNumberOfColors();
    Int_t    icol     = TMath::Min(nCol-1,int(fc * nCol + .5));
    Int_t    col      = gStyle->GetColorPalette(icol);
    //Info("GetCentralityColor","%3d: %3d-%3d -> %3d",bin,centLow,centHigh,col);
    return col;
  }
  /** 
   * Import a histogram into the data base
   * 
   * @param h            Histogram
   * @param title        Title on plot
   * @param experiment   Which experiement
   * @param sys          Collision system
   * @param sNN          Collision energy (in GeV)
   * @param trigger      Trigger type 
   * @param centLow      Lease centrality
   * @param centHigh     Largest centrality 
   * @param path         Possible path to database 
   * 
   * @return true on success
   */
  static Bool_t Import(TH1* h, 
		       const char* title, 
		       UShort_t experiment,
		       UShort_t sys, 
		       UShort_t sNN,
		       UShort_t trigger, 
		       UShort_t centLow=0, 
		       UShort_t centHigh=0,
		       const char* path=0)
  {
    TGraphAsymmErrors* g = new TGraphAsymmErrors();
    Int_t nx = h->GetNbinsX();
    Int_t j  = 0;
    for (Int_t i = 1; i <= nx; i++) { 
      Double_t x  = h->GetXaxis()->GetBinCenter(i);
      Double_t ex = h->GetXaxis()->GetBinWidth(i)/2;
      Double_t y  = h->GetBinContent(i);
      Double_t ey = h->GetBinError(i);
      
      if (TMath::Abs(y) < 1e-6 || ey < 1e-6) continue;
      
      g->SetPoint(j, x, y);
      g->SetPointError(j, ex, ex, ey, ey);
      j++;
    }
    if (j <= 0) return false;

    return Import(g, title, experiment, sys, sNN, trigger, centLow, centHigh, 
		  path);
  }
  /** 
   * Import a graph into the data base
   * 
   * @param g            Graph
   * @param title        Title on plot
   * @param experiment   Which experiement
   * @param sys          Collision system
   * @param sNN          Collision energy (in GeV)
   * @param trigger      Trigger type 
   * @param centLow      Lease centrality
   * @param centHigh     Largest centrality 
   * @param path         Possible path to database 
   * 
   * @return true on success
   */
  static Bool_t Import(TGraphAsymmErrors* g,
		       const char* title, 
		       UShort_t experiment,
		       UShort_t sys, 
		       UShort_t sNN,
		       UShort_t trigger, 
		       UShort_t centLow=0, 
		       UShort_t centHigh=0,
		       const char* path=0) 
  {
    if (!g) return false;

    TString     expName = ExpName(experiment); expName.ToLower();
    const char* sNNName = SNNName(sNN);
    const char* sysName = SysName(sys);        
    const char* trgName = TrgName(trigger,centLow,centHigh); 
    TString     name    = Form("%s%s%s%s", 
			       expName.Data(), sysName, sNNName, trgName);
    
    SetGraphAttributes(g, trigger, experiment, false, name, title);
    if (centLow < centHigh) { 
      Int_t col = CentralityColor(centLow, centHigh);
      g->SetMarkerColor(col);
      g->SetLineColor(col);
      g->SetFillColor(col);
    }
      
    TMultiGraph* mg = new TMultiGraph("data","");
    mg->Add(g);

    return Import(mg, experiment, sys, sNN, trigger, centLow, centHigh, path);
  }
  /** 
   * Import a graph into the data base
   * 
   * @param g            Graph
   * @param title        Title on plot
   * @param experiment   Which experiement
   * @param sys          Collision system
   * @param sNN          Collision energy (in GeV)
   * @param trigger      Trigger type 
   * @param centLow      Lease centrality
   * @param centHigh     Largest centrality 
   * @param path         Possible path to database 
   * 
   * @return true on success
   */
  static Bool_t Import(TMultiGraph* g,
		       UShort_t experiment,
		       UShort_t sys, 
		       UShort_t sNN,
		       UShort_t trigger, 
		       UShort_t centLow=0, 
		       UShort_t centHigh=0,
		       const char* path=0)
  {
    TFile* file = GetFile(path, true);
    
    const char* sysName = SysName(sys);
    const char* sNNName = SNNName(sNN);
    const char* trgName = TrgName(trigger, centLow, centHigh);
    const char* expName = ExpName(experiment);
    
    if (!sysName || !sNNName || !trgName || !expName) return false;
    
    TString dirName;
    dirName = Form("%s/%s/%s/%s", sysName, sNNName, trgName, expName);
    
    TDirectory* dir = file->GetDirectory(dirName);
    if (!dir) dir = file->mkdir(dirName);
    file->cd(dirName);

    if (dir->Get("data")) { 
      ::Warning("", "Already have data in %s", dirName.Data());
      // return false;
    }

    g->SetName("data");
    g->Write();
    
    file->cd();
    file->Close();

    return true;
  }
};

// EOF



    
      
