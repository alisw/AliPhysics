#include "SummaryDrawer.C"

/**
 * Class to draw a summary of the AOD production
 *
 * @par Input: 
 * - The merged <tt>forward.root</tt> file.
 *   If the file isn't merged, it should still work. 
 *
 * @par Output:
 * - A PDF file named after the input, but with <tt>.root</tt>
 *   replaced with <tt>pdf</tt>
 * 
 */
class SummaryAODDrawer : public SummaryDrawer
{
public:
  enum EFlags { 
    kEventInspector    = 0x001, 
    kSharingFilter     = 0x002, 
    kDensityCalculator = 0x004,
    kCorrector         = 0x008,
    kHistCollector     = 0x010,
    kSteps             = 0x020, 
    kResults           = 0x040, 
    kPause             = 0x080,
    kLandscape         = 0x100, 
    kCentral           = 0x200,
    kNormal            = 0x27F
  };
  SummaryAODDrawer() 
    : SummaryDrawer(),
      fSums(0),
      fResults(0)
  {}
  
  //__________________________________________________________________
  /** 
   * 
   * 
   * @param fname 
   * @param what 
   */
  void Run(const char* fname, UShort_t what=kNormal)
  {
    // --- Open the file ---------------------------------------------
    TString filename(fname);
    TFile*  file = TFile::Open(filename, "READ");
    if (!file) { 
      Error("Run", "Failed to open \"%s\"", filename.Data());
      return;
    }


    // --- Get top-level collection ----------------------------------
    fSums = GetCollection(file, "Forward");
    if (!fSums) return;

    // --- Make our canvas -------------------------------------------
    TString pdfName(filename);
    pdfName.ReplaceAll(".root", ".pdf");
    CreateCanvas(pdfName, what & 0x100);

    // --- Possibly make a chapter here ------------------------------
    if (what & kCentral && GetCollection(file, "Central")) 
      MakeChapter("Forward");
    
    // --- Set pause flag --------------------------------------------
    fPause = what & 0x80;

    // --- Do each sub-algorithm -------------------------------------
    if (what & 0x01) DrawEventInspector();
    if (what & 0x02) DrawSharingFilter();
    if (what & 0x04) DrawDensityCalculator();
    if (what & 0x08) DrawCorrector();
    if (what & 0x10) DrawHistCollector();
  
    // --- Do the results ----------------------------------------------
    fResults = GetCollection(file, "ForwardResults");
    if (!fResults) fResults = fSums; // Old-style
    
    if (what & 0x20) DrawSteps();
    if (what & 0x40) DrawResults();

    // --- SPD clusters ----------------------------------------------
    if (what & 0x200) { 
      // --- Get top-level collection --------------------------------
      fSums = GetCollection(file, "Central");
      if (fSums) {
	MakeChapter("Central");
	DrawCentral();
	if (what & 0x01) DrawEventInspector();

      }
    }

    CloseCanvas();
  }
protected:
  //____________________________________________________________________
  void DrawEventInspector()
  {
    Info("DrawEventInspector", "Drawing event inspector");
    TCollection* c = GetCollection(fSums, "fmdEventInspector");
    if (!c) return;

    Int_t sys=0, sNN=0, field=0, runNo=0, lowFlux=0, nPileUp=0;
    Int_t aliRev=0, aliBra=0;
    Bool_t fpVtx=false, v0and=false;
    Double_t dPileUp=0.;
    Double_t y = .8;

    fBody->cd();

    Double_t save = fParName->GetTextSize();
    fParName->SetTextSize(0.03);
    fParVal->SetTextSize(0.03);

    if (GetParameter(c, "sys", sys))
      DrawParameter(y, "System", (sys == 1 ? "pp" : sys == 2 ? "PbPb" : 
				  sys == 3 ? "pPb" : "unknown"));
    if (GetParameter(c, "sNN", sNN)) 
      DrawParameter(y, "#sqrt{s_{NN}}", Form("%5dGeV", sNN));

    if (GetParameter(c, "field", field))
      DrawParameter(y, "L3 B field", Form("%+2dkG", field));

    if (GetParameter(c, "runNo", runNo))
      DrawParameter(y, "Run #", Form("%6d", runNo));

    if (GetParameter(c, "lowFlux", lowFlux))
      DrawParameter(y, "Low flux cut", Form("%6d", lowFlux));

    if (GetParameter(c, "fpVtx", fpVtx))
      DrawParameter(y, "Use PWG-UD vertex", (fpVtx ? "yes" : "no"));

    if (GetParameter(c, "v0and", v0and))
      DrawParameter(y, "Use V0AND for NSD", (v0and ? "yes" : "no"));

    if (GetParameter(c, "nPileUp", nPileUp))
      DrawParameter(y, "Least # of pile-up vertex", Form("%d", nPileUp));
      
    if (GetParameter(c, "dPileup", dPileUp))
      DrawParameter(y, "Least distance of pile-up vertex",
		    Form("%fcm", dPileUp));

    if (GetParameter(c, "alirootRev", aliRev) || 
	GetParameter(c, "alirootBranch", aliBra))
      DrawParameter(y, "AliROOT", Form("%7lu/0x%08lx", ULong_t(aliRev), 
				       ULong_t(aliBra)));

    PrintCanvas("Event Inspector");
    fParName->SetTextSize(save);
    fParVal->SetTextSize(save);

    fBody->Divide(2,4);
    
    TH1*    nEventsTr    = GetH1(c, "nEventsTr");
    TH1*    nEventsTrVtx = GetH1(c, "nEventsTrVtx");
    TH1*    vertex       = GetH1(c, "vertex", false);
    Bool_t  mc           = (vertex != 0);
    if (nEventsTr)    nEventsTr->Rebin(2);
    if (nEventsTrVtx) nEventsTrVtx->Rebin(2);
    if (vertex) {
      // vertex->Rebin(2);
      vertex->SetFillColor(kMagenta+2);
    }
    DrawInPad(fBody, 1, nEventsTr, "", 0x2);
    DrawInPad(fBody, 1, vertex, "same");
    DrawInPad(fBody, 1, nEventsTrVtx, "same"); 
    DrawInPad(fBody, 1, GetH1(c, "nEventsAccepted"), "same", 0x10);


    DrawInPad(fBody, 2, GetH2(c, "nEventsAcceptedXY"), "colz", 0x4);
    DrawInPad(fBody, 3, GetH1(c, "triggers"),          "hist text");
    if (GetH1(c, "trgStatus"))
      DrawInPad(fBody, 4, GetH1(c, "trgStatus"),         "hist text");
    else  // Old one 
      DrawInPad(fBody, 4, GetH2(c, "triggerCorr"),       "colz", 0x4);
    DrawInPad(fBody, 5, GetH1(c, "status"),            "hist text");
    if (GetH1(c, "vtxStatus"))
      DrawInPad(fBody, 6, GetH1(c, "vtxStatus"),         "hist text");
    else // old 
      DrawInPad(fBody, 6, GetH1(c, "type"),              "hist text");
    DrawInPad(fBody, 7, GetH1(c, "cent"));
    DrawInPad(fBody, 8, GetH2(c, "centVsQuality"), "colz", 0x4);

    PrintCanvas("EventInspector - Histograms");  

    if (!mc) return; // not MC 
  
    TH1* phiR         = GetH1(c, "phiR");
    TH1* b            = GetH1(c, "b");
    TH2* bVsNpart     = GetH2(c, "bVsParticipants");
    TH2* bVsNbin      = GetH2(c, "bVsBinary");
    TH2* bVsCent      = GetH2(c, "bVsCentrality");
    TH2* vzComparison = GetH2(c, "vzComparison");
    TH2* centVsNpart  = GetH2(c, "centralityVsParticipans");// Spelling!
    TH2* centVsNbin   = GetH2(c, "centralityVsBinary");
  
    fBody->Divide(2,3);

    DrawInPad(fBody, 1, phiR);
    DrawInPad(fBody, 2, vzComparison, "colz", 0x4);
    DrawInPad(fBody, 3, b);

    TProfile* nPartB = bVsNpart->ProfileX("nPartB",1,-1,"s");
    TProfile* nBinB  = bVsNbin->ProfileX("nBinB",1,-1,"s");
    nPartB->SetMarkerColor(kBlue+2);
    nPartB->SetMarkerStyle(20);
    nPartB->SetLineColor(kBlue+2);
    nPartB->SetFillColor(kBlue-10);
    nPartB->SetFillStyle(1001);
    nPartB->SetMarkerSize(0.7);
    nBinB->SetMarkerColor(kRed+2);
    nBinB->SetMarkerStyle(21);
    nBinB->SetLineColor(kRed+2);
    nBinB->SetFillColor(kRed-10);
    nBinB->SetMarkerSize(0.7);
    nBinB->SetFillStyle(1001);

    DrawTwoInPad(fBody, 4, nPartB, nBinB, "e3 p", 0x10);

    DrawInPad(fBody, 5, bVsCent, "colz", 0x4);

    TProfile* nPartC = centVsNpart->ProfileY("nPartC",1,-1,"s");
    TProfile* nBinC  = centVsNbin->ProfileY("nBinC",1,-1,"s");
    nPartC->SetMarkerColor(kBlue+2);
    nPartC->SetMarkerStyle(20);
    nPartC->SetLineColor(kBlue+2);
    nPartC->SetFillColor(kBlue-10);
    nPartC->SetFillStyle(1001);
    nPartC->SetMarkerSize(0.7);
    nBinC->SetMarkerColor(kRed+2);
    nBinC->SetMarkerStyle(21);
    nBinC->SetLineColor(kRed+2);
    nBinC->SetFillColor(kRed-10);
    nBinC->SetMarkerSize(0.7);
    nBinC->SetFillStyle(1001);

    DrawTwoInPad(fBody, 6, nPartC, nBinC, "e3 p", 0x10);

    PrintCanvas("EventInspector - Monte-Carlo");  
  }
  //____________________________________________________________________
  void DrawSharingFilter()
  {
    Info("DrawEventInspector", "Drawing sharing filter");
    TCollection* c = GetCollection(fSums, "fmdSharingFilter");
    if (!c) return;

    fBody->Divide(1, 3);
    fBody->cd(1);
  
    Double_t y = .8;
    Bool_t   angle=false, lowSignal=false, simple=false;
  
    if (GetParameter(c, "angle", angle))
      DrawParameter(y, "Angle correct", (angle ? "yes" : "no")); 
    if (GetParameter(c, "lowSignal", lowSignal))
      DrawParameter(y, "Lower signal",  (lowSignal ? "yes" : "no"));

    if (GetParameter(c, "simple", simple))
      DrawParameter(y, "Simple method", (simple ? "yes" : "no"));


    DrawInPad(fBody, 2, GetH2(c, "lowCuts"), "colz");
    DrawInPad(fBody, 3, GetH2(c, "highCuts"), "colz");
  
    PrintCanvas("Sharing filter");

    const char* subs[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3O", "FMD3I", 0 };
    const char** ptr   = subs;
    while (*ptr) { 
      TCollection* sc = GetCollection(c, *ptr);
      if (!sc) { ptr++; continue; }
    
      fBody->Divide(2,3);
      DrawInPad(fBody, 1, GetH1(sc, "esdEloss"),       "",     0x2);
      DrawInPad(fBody, 1, GetH1(sc, "anaEloss"),       "same", 0x12);
      DrawInPad(fBody, 2, GetH1(sc, "singleEloss"),    "",     0x2);
      DrawInPad(fBody, 2, GetH1(sc, "doubleEloss"),    "same", 0x2);
      DrawInPad(fBody, 2, GetH1(sc, "tripleEloss"),    "same", 0x12);  
      DrawInPad(fBody, 3, GetH2(sc, "singlePerStrip"), "colz", 0x4);
      // DrawInPad(fBody, 4, GetH1(sc, "distanceBefore"), "",     0x2);
      // DrawInPad(fBody, 4, GetH1(sc, "distanceAfter"),  "same", 0x12);
      DrawInPad(fBody, 4, GetH2(sc, "summed"),         "colz", 0x0);

      TH2* nB = GetH2(sc, "neighborsBefore");
      if (nB) { 
	nB->GetXaxis()->SetRangeUser(0,8); 
	nB->GetYaxis()->SetRangeUser(0,8); 
      }
      DrawInPad(fBody, 5, nB, "colz", 0x4);
      DrawInPad(fBody, 5, GetH2(sc, "neighborsAfter"), "p same", 0x4);
      DrawInPad(fBody, 6, GetH2(sc, "beforeAfter"),    "colz",   0x4);

      PrintCanvas(Form("Sharing filter - %s", *ptr));
      ptr++;
    }

    // --- MC --------------------------------------------------------
    TCollection* cc = GetCollection(c, "esd_mc_comparion", false); // Spelling!
    if (!cc) return; // Not MC 

    DivideForRings(false, false);
    DrawInRingPad(1, 'I', GetH2(cc, "FMD1i_corr"), "colz", 0x4);
    DrawInRingPad(2, 'I', GetH2(cc, "FMD2i_corr"), "colz", 0x4);
    DrawInRingPad(2, 'O', GetH2(cc, "FMD2o_corr"), "colz", 0x4);
    DrawInRingPad(3, 'I', GetH2(cc, "FMD3i_corr"), "colz", 0x4);
    DrawInRingPad(3, 'O', GetH2(cc, "FMD3o_corr"), "colz", 0x4);

    PrintCanvas("Sharing filter - MC vs Reco");

    // --- MC --------------------------------------------------------
    TCollection* mc = GetCollection(c, "mcTrackDensity", false);
    if (!mc) return; // Not MC 

    fBody->Divide(2,3);
    DrawInPad(fBody, 1, GetH2(mc, "binFlow"),    "colz", 0x4);
    DrawInPad(fBody, 2, GetH2(mc, "binFlowEta"), "colz", 0x4);
    DrawInPad(fBody, 3, GetH2(mc, "binFlowPhi"), "colz", 0x4);
    DrawInPad(fBody, 4, GetH1(mc, "nRefs"),       "",    0x2);
    DrawInPad(fBody, 4, GetH1(mc, "clusterRefs"), "same");
    DrawInPad(fBody, 4, GetH1(mc, "clusterSize"), "same");
    DrawInPad(fBody, 4, GetH1(mc, "nClusters"),    "same", 0x10);
    DrawInPad(fBody, 5, GetH2(mc, "clusterVsRefs"),"colz", 0x4);

    PrintCanvas("Sharing filter - MC");  
  }

  //____________________________________________________________________
  void DrawDensityCalculator()
  {
    Info("DrawDensityCalculator", "Drawing density calculator");
    TCollection* c = GetCollection(fSums, "fmdDensityCalculator");
    if (!c) return;

    fBody->Divide(2, 2);
    fBody->cd(1);
  
    Double_t y = .8;
    Int_t maxParticles=0, phiAcceptance=0, etaLumping=0, phiLumping=0;
    Bool_t method=false, recalcEta=false, recalcPhi=false;
    Double_t size = fLandscape ? 0.06 : 0.04;
  
    GetParameter(c, "maxParticle", maxParticles);

    if (GetParameter(c, "phiAcceptance", phiAcceptance))
      DrawParameter(y, "#phi acceptance method", 
		    (phiAcceptance == 1 ? "N_{ch}" : 
		     phiAcceptance == 2 ? "#DeltaE" : "none"),       size);
    
    if (GetParameter(c, "etaLumping", etaLumping) &&
	GetParameter(c, "phiLumping", phiLumping))
      DrawParameter(y, "Region size (sector#timesstrip)", 
		    Form("%2d #times %2d", phiLumping, etaLumping),  size);
    
    if (GetParameter(c, "method", method))
      DrawParameter(y, "Method", (method ? "Poisson" : "#DeltaE"),   size); 

    if (GetParameter(c, "recalcEta", recalcEta))
      DrawParameter(y, "Recalculate #eta",(recalcEta ? "yes" : "no"),size); 

    if (GetParameter(c, "recalcPhi", recalcPhi))
      DrawParameter(y, "Recalculate #phi",(recalcPhi ? "yes" : "no"),size); 



    TVirtualPad* p = fBody; // fBody->cd(2);
    // p->Divide(3,1);

    TH1* accI = GetH1(c, "accI");
    TH1* accO = GetH1(c, "accO");
    if (accI) { 
      Double_t scale = 1./accI->GetMaximum();
      accI->Scale(scale); 
      accO->Scale(scale);
      accI->SetMinimum(0); 
    }
    DrawInPad(p, 2, accI); 
    DrawInPad(p, 2, accO, "same", 0x10); 
    DrawInPad(p, 3, GetH2(c, "lowCuts"), "colz");
    DrawInPad(p, 4, GetH2(c, "maxWeights"), "colz");
  
    PrintCanvas("Density calculator");

    const char* subs[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3O", "FMD3I", 0 };
    const char** ptr   = subs;
    while (*ptr) { 
      TCollection* sc = GetCollection(c, *ptr);
      if (!sc) { ptr++; continue; }
    
      fBody->Divide(2,3);
    
      DrawInPad(fBody, 1, GetH2(sc, "elossVsPoisson"), "colz",   0x4);
      DrawInPad(fBody, 2, GetH1(sc, "diffElossPoisson"), "",     0x2);
      DrawInPad(fBody, 3, GetH1(sc, "occupancy"),        "",     0x2);
      DrawInPad(fBody, 4, GetH1(sc, "eloss"),            "",     0x2);
      DrawInPad(fBody, 4, GetH1(sc, "elossUsed"),        "same", 0x12);
      TH1* phiB = GetH1(sc, "phiBefore");
      TH1* phiA = GetH1(sc, "phiAfter");
      if (phiB && phiA) { 
	phiA->Add(phiB, -1);
	phiA->Divide(phiB);
	phiA->SetTitle("#Delta#phi from Ip (x,y) correction");
	phiA->SetYTitle("(#phi_{after}-#phi_{before})/#phi_{before}");
      }
      DrawInPad(fBody, 5, phiA);
      DrawInPad(fBody, 6, GetH2(sc, "phiAcc"), "colz",   0x4);
    
      PrintCanvas(Form("Density calculator - %s", *ptr));
      ptr++;    
    }

    TCollection* cc = GetCollection(c, "esd_mc_comparison", false); 
    if (!cc) return; // Not MC 

    fBody->Divide(2,5);
    DrawInPad(fBody, 1, GetH2(cc, "FMD1I_corr_mc_esd"), "colz", 0x4);
    DrawInPad(fBody, 3, GetH2(cc, "FMD2I_corr_mc_esd"), "colz", 0x4);
    DrawInPad(fBody, 5, GetH2(cc, "FMD2O_corr_mc_esd"), "colz", 0x4);
    DrawInPad(fBody, 7, GetH2(cc, "FMD3O_corr_mc_esd"), "colz", 0x4);
    DrawInPad(fBody, 9, GetH2(cc, "FMD3I_corr_mc_esd"), "colz", 0x4);
    DrawInPad(fBody, 2,  GetH1(cc, "FMD1I_diff_mc_esd"), "", 0x2);
    DrawInPad(fBody, 4,  GetH1(cc, "FMD2I_diff_mc_esd"), "", 0x2);
    DrawInPad(fBody, 6,  GetH1(cc, "FMD2O_diff_mc_esd"), "", 0x2);
    DrawInPad(fBody, 8,  GetH1(cc, "FMD3O_diff_mc_esd"), "", 0x2);
    DrawInPad(fBody, 10, GetH1(cc, "FMD3I_diff_mc_esd"), "", 0x2);

    PrintCanvas("Density calculator - MC vs Reco");
  }

  //____________________________________________________________________
  void DrawCorrector()
  {
    Info("DrawCorrector", "Drawing corrector"); 
    TCollection* c = GetCollection(fSums, "fmdCorrector"); 
    if (!c) return;
  
    fBody->cd();
  
    Double_t y = .8;  
    Bool_t secondary=false, vertexBias=false, acceptance=false, merging=false;  
    if (GetParameter(c, "secondary", secondary))
      DrawParameter(y, "Secondary corr.", secondary ? "yes" : "no");
    if (GetParameter(c, "acceptance", acceptance))
      DrawParameter(y, "Acceptance corr.", acceptance ? "yes" : "no");
    if (GetParameter(c, "vertexBias", vertexBias))
      DrawParameter(y, "Vertex bias corr.", vertexBias ? "yes" : "no");
    if (GetParameter(c, "merging", merging))  
      DrawParameter(y, "Merging eff.", merging ? "yes" : "no");
    
    PrintCanvas("Corrector");

    TCollection* cc = GetCollection(c, "esd_mc_comparison", false); 
    if (!cc) return; // Not MC 
    
    DivideForRings(false, false);

    DrawInRingPad(1, 'I', GetH2(cc, "FMD1I_esd_vs_mc"), "colz", 0x0);
    DrawInRingPad(2, 'I', GetH2(cc, "FMD2I_esd_vs_mc"), "colz", 0x0);
    DrawInRingPad(2, 'O', GetH2(cc, "FMD2O_esd_vs_mc"), "colz", 0x0);
    DrawInRingPad(3, 'O', GetH2(cc, "FMD3O_esd_vs_mc"), "colz", 0x0);
    DrawInRingPad(3, 'I', GetH2(cc, "FMD3I_esd_vs_mc"), "colz", 0x0);

    PrintCanvas("Corrector - MC vs Reco");
  }

  //____________________________________________________________________
  void DrawHistCollector()
  {
    Info("DrawHistCollector", "Drawing histogram collector");  
    TCollection* c = GetCollection(fSums, "fmdHistCollector");
    if (!c) return;

    fBody->Divide(1, 3);
    fBody->cd(1);

    Double_t y = .8;
    Int_t nCutBins=0, fiducial=0, merge=0, skipRings=0;
    Double_t fiducialCut=0.;
    Bool_t  bgAndHits=false;

    if (GetParameter(c, "nCutBins", nCutBins))
      DrawParameter(y, "# of bins to cut",      Form("%d", nCutBins));

    if (GetParameter(c, "skipRings", skipRings)) {
      TString skipped;
      if (skipRings & 0x05) skipped.Append("FMD1i ");
      if (skipRings & 0x09) skipped.Append("FMD2i ");
      if (skipRings & 0x0a) skipped.Append("FMD2o ");
      if (skipRings & 0x11) skipped.Append("FMD3i ");
      if (skipRings & 0x12) skipped.Append("FMD3o ");
      DrawParameter(y, "Skipped rings", skipped);
    }

    if (GetParameter(c, "bgAndHits", bgAndHits))
      DrawParameter(y, "Bg & hit maps stored.", bgAndHits?"yes":"no");

    if (GetParameter(c, "merge", merge))
      DrawParameter(y, "Merge method", 
		    (merge == 0 ? "straight mean" :
		     merge == 1 ? "straight mean, no zeroes" : 
		     merge == 2 ? "weighted mean" : 
		     merge == 3 ? "least error" : 
		     merge == 4 ? "sum" : "unknown"));

    if (GetParameter(c, "fiducial", fiducial))
      DrawParameter(y, "Fiducial method.", 
		    fiducial == 0 ? "cut" : "distance");

    if (GetParameter(c, "fiducialCut", fiducialCut))
      DrawParameter(y, "Fiducial cut.", Form("%f", fiducialCut));

		 
    DrawInPad(fBody, 2, GetH2(c, "sumRings"), "colz"); 
    DrawInPad(fBody, 3, GetH2(c, "coverage"), "colz");

    fBody->cd(1)->Modified();
    fBody->cd(2)->Modified();
    fBody->cd(3)->Modified();
    fBody->cd(1)->Update();
    fBody->cd(2)->Update();
    fBody->cd(3)->Update();
    PrintCanvas("Histogram collector");
		
    
    TIter next(c);
    TObject* o = 0;
    TRegexp regexp("[pm][0-9]+_[pm][0-9]+");
    while ((o = next())) { 
      TString name(o->GetName());
      if (name.Index(regexp) == kNPOS) continue;
      
      TList* vl = static_cast<TList*>(o);

      DivideForRings(false, false);
      
      DrawInRingPad(1, 'I', GetH2(vl, "secMapFMD1I"), "colz", 0x0);
      DrawInRingPad(2, 'I', GetH2(vl, "secMapFMD2I"), "colz", 0x0);
      DrawInRingPad(2, 'O', GetH2(vl, "secMapFMD2O"), "colz", 0x0);
      DrawInRingPad(3, 'O', GetH2(vl, "secMapFMD3O"), "colz", 0x0);
      DrawInRingPad(3, 'I', GetH2(vl, "secMapFMD3I"), "colz", 0x0);
      DrawInRingPad(1, 'I', GetH2(vl, "hitMapFMD1I"), "box same", 0x0);
      DrawInRingPad(2, 'I', GetH2(vl, "hitMapFMD2I"), "box same", 0x0);
      DrawInRingPad(2, 'O', GetH2(vl, "hitMapFMD2O"), "box same", 0x0);
      DrawInRingPad(3, 'O', GetH2(vl, "hitMapFMD3O"), "box same", 0x0);
      DrawInRingPad(3, 'I', GetH2(vl, "hitMapFMD3I"), "box same", 0x0);

      PrintCanvas(Form("Histogram Collector - Vertex bin %s", vl->GetName()));
    }

    o = c->FindObject("byCentrality");
    if (!o) return;
    TList* bc = static_cast<TList*>(o);

    DrawInPad(fBody, GetH3(bc, "FMD1I"), "box", 0);
    DrawInPad(fBody, GetH3(bc, "FMD2I"), "box same", 0);
    DrawInPad(fBody, GetH3(bc, "FMD2O"), "box same", 0);
    DrawInPad(fBody, GetH3(bc, "FMD3O"), "box same", 0);
    DrawInPad(fBody, GetH3(bc, "FMD3I"), "box same", 0x10);
  }

  //____________________________________________________________________
  void DrawCentral()
  {
    Info("DrawCentral", "Drawing central (SPD)");  
    TCollection* c = fSums; 
    if (!c) return;

    fBody->Divide(1, 3);
    fBody->cd(1);

		 
    DrawInPad(fBody, 1, GetH2(c, "coverage"), "col", 0);
    DrawInPad(fBody, 2, GetH2(c, "nClusterVsnTracklet"), "colz", 0x03); 
    DrawInPad(fBody, 3, GetH2(c, "clusterPerTracklet"), "colz", 0x0); 

    fBody->cd(1)->Modified();
    fBody->cd(2)->Modified();
    fBody->cd(3)->Modified();
    fBody->cd(1)->Update();
    fBody->cd(2)->Update();
    fBody->cd(3)->Update();
    PrintCanvas("Central - overview");
		
    
    TIter next(c);
    TObject* o = 0;
    TRegexp regexp("[pm][0-9]+_[pm][0-9]+");
    while ((o = next())) { 
      TString name(o->GetName());
      if (name.Index(regexp) == kNPOS) continue;
      
      TList* vl = static_cast<TList*>(o);

      fBody->Divide(1, 3);
    
      DrawInPad(fBody, 1, GetH1(vl, "acceptance"), "", 0);

      TH1* sec = GetH1(vl, "secondary");
      sec->SetMarkerStyle(21);
      sec->SetMarkerSize(1.2);
      DrawInPad(fBody, 2, sec, "", 0);
      DrawInPad(fBody, 2, GetH1(vl, "secondaryFiducial"),    "same", 0x0);
      DrawInPad(fBody, 3, GetH2(vl, "secondaryMapFiducial"), "colz", 0);
      DrawInPad(fBody, 3, GetH2(vl, "hitMap"),               "box same", 0x0);

      fBody->cd(1)->Modified();
      fBody->cd(2)->Modified();
      fBody->cd(3)->Modified();
      fBody->cd(1)->Update();
      fBody->cd(2)->Update();
      fBody->cd(3)->Update();
      PrintCanvas(Form("Central - Vertex bin %s", vl->GetName()));
    }
  }

  
  //____________________________________________________________________
  void AddToAll(THStack* all, const THStack* stack, Int_t curr, Int_t step)
  {
    if (!stack) return;

    TIter   next(stack->GetHists());
    TH1*    h = 0;
    while ((h = static_cast<TH1*>(next()))) {
      TH1* copy = static_cast<TH1*>(h->Clone(Form("%s_copy", h->GetName())));
      copy->SetDirectory(0);
      if (curr != step) {
	copy->SetMarkerColor(kGray);
	copy->SetLineColor(kGray);
      }
      all->Add(copy);
    }
  }
  //____________________________________________________________________
  void AddToAll(THStack* all, const THStack* stack)
  {
    if (!stack) return;

    TIter   next(stack->GetHists());
    TH1*    h = 0;
    while ((h = static_cast<TH1*>(next()))) {
      TH1* copy = static_cast<TH1*>(h->Clone(Form("%s_copy", h->GetName())));
      copy->SetDirectory(0);
      copy->SetMarkerColor(kGray);
      copy->SetLineColor(kGray);
      all->Add(copy);
    }
  }

  //____________________________________________________________________
  void DrawStep(Int_t        step,
		THStack*     all,
		TObject*     cur,
		TLegend*     leg,
		const char*  title,
		TVirtualPad* can)
  {
    if (all->GetHists()->GetEntries() <= 0 || !cur) return;

    // Info("", "Drawing step # %d", step);
    Bool_t left = (step % 2) == 1; 
    TVirtualPad* p = can->cd(step);
    gStyle->SetOptTitle(0);
    p->SetTitle(Form("Step # %d", step));
    p->SetFillColor(kWhite);
    p->SetRightMargin(left ? 0 : 0.02);
    p->SetTopMargin(0); // 0.02);

    p->cd();
    all->Draw("nostack");
    all->GetHistogram()->SetXTitle("#eta");
    all->GetHistogram()->SetYTitle("signal");

    // p->cd();
    gROOT->SetSelectedPad(p);
    cur->DrawClone("same nostack");
    leg->DrawClone("");

    TLatex* ltx = new TLatex(.97, .97, title);
    ltx->SetNDC();
    ltx->SetTextSize(.06);
    ltx->SetTextAlign(33);
    ltx->Draw();

    ltx = new TLatex((left ? .12 : .02), .97, p->GetTitle());
    ltx->SetNDC();
    ltx->SetTextSize(.06);
    ltx->SetTextAlign(13);
    ltx->Draw();

    p->Modified();
    p->Update();
    p->cd();
  }

  //____________________________________________________________________
  void FixStack(THStack* stack, const TString& title, 
		const TString& extra, Int_t marker)
  {
    if (!stack) return;
    stack->SetTitle(title);
    TIter next(stack->GetHists());
    TH1*  h = 0;
    while ((h = static_cast<TH1*>(next())))  {
      h->SetMarkerStyle(marker);
      TString tit(h->GetTitle());
      tit.ReplaceAll("cache", "");
      tit.Append(extra);
      h->SetTitle(tit);
    }
  }
  void AddLegendEntry(TLegend* l, 
		      const TH1* h, 
		      const TString& title)
  {
    if (!h) return;

    TLegendEntry* e = l->AddEntry("dummy", title.Data(), "pl");
    e->SetMarkerStyle(h->GetMarkerStyle());
    e->SetMarkerColor(kGray);
    e->SetLineColor(kGray);
    e->SetTextColor(kGray);
  }
		      

  //____________________________________________________________________
  void DrawSteps()
  {
    // MakeChapter(can, "Steps");

    THStack* deltas  = GetStack(GetCollection(fResults, "fmdSharingFilter"), 
				"sums", "summed");
    THStack* nchs    = GetStack(GetCollection(fResults, 
					      "fmdDensityCalculator"), 
				"sums", "inclDensity");
    THStack* prims   = GetStack(GetCollection(fResults, "fmdCorrector"), 
				"sums", "primaryDensity");
    THStack* rings   = GetStack(GetCollection(fResults, "ringResults"), "all");
    THStack* mcRings = GetStack(GetCollection(fResults, "mcRingResults", false),
				"all","dndeta_eta", false);
    TH1*     dndeta  = GetH1(fResults, "dNdeta");
    if (dndeta) dndeta->SetMarkerColor(kBlack);

    FixStack(deltas, "#sum_{} #Delta/#Delta_{mip}",  "",     20);
    FixStack(nchs,   "#sum_{} N_{ch,incl}", 	     "",     21);
    FixStack(prims,  "#sum_{} N_{ch,primary}",       "",     22);
    FixStack(rings,  "dN/d#eta per ring",            "",     23);
    FixStack(mcRings,"dN/d#eta per ring (MC)",       "(MC)", 34);

    THStack* all = new THStack;
    AddToAll(all, mcRings);
    AddToAll(all, deltas);
    AddToAll(all, nchs);
    AddToAll(all, prims);
    AddToAll(all, rings);

    TH1* res = 0;
    if (dndeta) {
      res = static_cast<TH1*>(dndeta->Clone("dNdeta"));
      res->SetTitle("dN/d#eta");
      res->SetMarkerColor(kGray);
      res->SetLineColor(kGray);
      res->SetDirectory(0);
      all->Add(res);
    }

    TLegend* l = new TLegend(.35, .2, .55, .9);
    l->SetFillColor(kWhite);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    TLegendEntry* e = 0;

    TH1* h = 0;
    if (mcRings) {
      h = static_cast<TH1*>(mcRings->GetHists()->At(0));
      AddLegendEntry(l, h, mcRings->GetTitle());
    }

    if (deltas) {
      h = static_cast<TH1*>(deltas->GetHists()->At(0));
      AddLegendEntry(l, h, deltas->GetTitle());    
    }

    if (nchs) {
      h = static_cast<TH1*>(nchs->GetHists()->At(0));
      AddLegendEntry(l, h, nchs->GetTitle());    
    }

    if (prims) {
      h = static_cast<TH1*>(prims->GetHists()->At(0));
      AddLegendEntry(l, h, prims->GetTitle());    
    }

    if (rings) {
      h = static_cast<TH1*>(rings->GetHists()->At(0));
      AddLegendEntry(l, h, rings->GetTitle());    
    }

    if (res) {
      h = res;
      AddLegendEntry(l, h, h->GetTitle());
    }
    
    TObject* objs[] = { mcRings, 
			deltas, 
			nchs, 
			prims, 
			rings, 
			dndeta };
    const char* titles[] = { "MC", 
			     "After merging", 
			     "After particle counting", 
			     "After corrections", 
			     "After normalization", 
			     "After combining" };
    
    fBody->Divide(2, 3, 0, 0); // mcRings ? 4 : 3, 0, 0);
    Int_t step = 0;
    for (Int_t i = 0; i < 6; i++) { 
      TObject* obj = objs[i];
      if (!obj) continue;
      e = static_cast<TLegendEntry*>(l->GetListOfPrimitives()->At(step));
      e->SetMarkerColor(kBlack);
      e->SetLineColor(kBlack);
      e->SetTextColor(kBlack);
      step++;
      
      DrawStep(step, all, obj, l, titles[i], fBody);
      e->SetMarkerColor(kGray);
      e->SetLineColor(kGray);
      e->SetTextColor(kGray);
    }

    if (!mcRings && deltas) { 
      fBody->cd(6);
      TLegend* ll = new TLegend(0.01, 0.11, 0.99, 0.99);
      // ll->SetNDC();
      ll->SetFillColor(kWhite);
      ll->SetFillStyle(0);
      ll->SetBorderSize(0);

      TIter next(deltas->GetHists());
      TH1*  hh = 0;
      while ((hh = static_cast<TH1*>(next()))) {
	e = ll->AddEntry("dummy", hh->GetTitle(), "pl");
	e->SetMarkerColor(hh->GetMarkerColor());
	e->SetMarkerStyle(hh->GetMarkerStyle());
	e->SetLineColor(kBlack);
      }
      ll->Draw();
    }
    PrintCanvas("Steps");
  }


  //____________________________________________________________________
  void DrawResults()
  {
    // MakeChapter(can, "Results");

    fBody->Divide(2,2);

    TCollection* c = GetCollection(fResults, "ringResults");
    if (!c) return;
  
    THStack* mcRings = GetStack(GetCollection(fResults, "mcRingResults", false),
				"all", "dndeta_eta", false);

    TH1* dndeta_phi = GetH1(fResults, "dNdeta");
    TH1* dndeta_eta = GetH1(fResults, "dNdeta_");
    dndeta_phi->SetTitle("1/N_{ev}dN_{ch}/d#eta (#varphi norm)");
    dndeta_eta->SetTitle("1/N_{ev}dN_{ch}/d#eta (#eta norm)");
    dndeta_eta->SetMarkerSize(0.7);

    DrawInPad(fBody, 1, GetStack(c, "all"), "nostack", mcRings ? 0 : 0x10);
    DrawInPad(fBody, 2, dndeta_phi);
    DrawInPad(fBody, 2, dndeta_eta, "Same", 0x10);
    DrawInPad(fBody, 3, GetH1(fResults, "norm"));
    DrawInPad(fBody, 3, GetH1(fResults, "phi"), "same", 0x10);
    DrawInPad(fBody, 4, GetH1(fSums,    "d2Ndetadphi"), "colz");
    DrawInPad(fBody, 1, mcRings, "nostack same", 0x10);

    fBody->cd(1);
    TLatex* l = new TLatex(.5, .2, "Ring results");
    l->SetNDC();
    l->SetTextAlign(21);
    l->Draw();

    fBody->cd(2);
    l->DrawLatex(.5, .2, "1/N_{ev}dN_{ch}/d#eta");

    // fBody->cd(3);
    // l->DrawLatex(.5, .2, "1/N_{ev}dN_{ch}/d#eta (#vta norm.)");
    
    PrintCanvas("Results");
  }

  TCollection* fSums;
  TCollection* fResults;
};

// #endif
