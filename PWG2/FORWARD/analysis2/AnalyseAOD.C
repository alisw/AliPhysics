/**
 * @file 
 * 
 * @ingroup pwg2_forward_scripts
 */
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TError.h>
#include <TSystemDirectory.h>
#include <TROOT.h>
#include <TStyle.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TSystem.h>
#include "AliAODForwardMult.h"
#include "OtherData.C"

/** 
 * Draw the data stored in the AOD 
 *
 * To use this, do 
 * 
 * @code 
 *   Root> .L $ALICE_ROOT/PWG2/FORWARD/analysis2/Compile.C
 *   Root> Compile("$ALICE_ROOT/PWG2/FORWARD/analysis2/AnalyseAOD.C")
 *   Root> AnalyseAOD dr
 *   Root> dr.Run("AliAODs.root",-10,10,5,AliAODForwardMult::kInel,900)
 * @endcode 
 * 
 * The output is stored in a ROOT file 
 * 
 * See also the script Pass2.C 
 * 
 * @ingroup pwg2_forward_scripts
 */
struct AnalyseAOD 
{
public: 
  /** AOD tree */
  TTree*             fTree;
  /** AOD object */
  AliAODForwardMult* fAOD;
  /** AOD object */
  AliAODForwardMult* fMCAOD;
  /** Output file */
  TFile*             fOut;
  /** Summed histogram */
  TH2D*              fSum;
  /** Summed histogram */
  TH2D*              fMCSum;
  /** Primary information */
  TH2D*              fPrimary;
  /** Sum primary information */
  TH2D*              fSumPrimary;
  /** Central region */
  TH2D*              fCentral;
  /** Central region */
  TH2D*              fSumCentral;
  /** Vertex efficiency */
  Double_t           fVtxEff;
  /** Title to put on the plot */
  TString            fTitle;
  /** Do HHD comparison */
  Bool_t             fDoHHD;
  /** Number of events with a trigger */
  Int_t              fNTriggered;
  /** Number of events with a vertex */
  Int_t              fNWithVertex;
  /** Number of events accepted */
  Int_t              fNAccepted;
  /** Number of events accepted */
  Int_t              fNAll;
  /** Number of events accepted */
  Int_t              fNMB;
  /** Number of B triggers */
  Int_t              fNB;
  /** Number of A triggers */
  Int_t              fNA;
  /** Number of C triggers */
  Int_t              fNC;
  /** Number of E triggers */
  Int_t              fNE;
  /** Min vertex */
  Double_t           fVtxMin;
  /** Max vertex */
  Double_t           fVtxMax;
  /** Whether to cut edges in the rebinning */
  Bool_t             fCutEdges;
    
  //__________________________________________________________________
  /** 
   * Constructor 
   * 
   */
  AnalyseAOD()
    : fTree(0), 
      fAOD(0),
      fMCAOD(0),
      fOut(0), 
      fSum(0),
      fMCSum(0),
      fPrimary(0), 
      fSumPrimary(0),
      fCentral(0),
      fSumCentral(0),
      fVtxEff(0),
      fTitle(""),
      fDoHHD(kTRUE),
      fNTriggered(0),
      fNWithVertex(0),
      fNAccepted(0),
      fNAll(0),
      fNMB(0),
      fNB(0),
      fNA(0),
      fNC(0),
      fNE(0),
      fVtxMin(0),
      fVtxMax(0),
      fCutEdges(false)
  {}
  //__________________________________________________________________
  /** 
   * Reset internal structures 
   * 
   */
  void Clear(Option_t* )
  {
    if (fTree && fTree->GetCurrentFile()) { 
      fTree->GetCurrentFile()->Close();
      delete fTree;
    }
    if (fOut) {
      fOut->Close();
      delete fOut;
    }
    if (fSum)        delete fSum;
    if (fMCSum)      delete fMCSum;
    if (fSumPrimary) delete fSumPrimary;
    if (fCentral)    delete fCentral;
    if (fSumCentral) delete fSumCentral;
    fTree       = 0;
    fOut        = 0;
    fSum        = 0;
    fMCSum      = 0;
    fPrimary    = 0;
    fSumPrimary = 0;
    fCentral    = 0;
    fSumCentral = 0;
    fVtxEff     = 0;
    fNTriggered	= 0;
    fNWithVertex= 0;
    fNAccepted	= 0;
    fNAll	= 0;
    fNMB        = 0;
    fNB		= 0;
    fNA		= 0;
    fNC		= 0;
    fNE		= 0;
    // fVtxMin	= 0;
    // fVtxMax	= 0;
    
  }
  //__________________________________________________________________
  /** 
   * Run 
   * 
   * @param file   Input path for files with AOD trees
   * @param vzMin  Minimum interaction point z coordinate 
   * @param vzMax  Maximum interaction point z coordinate 
   * @param rebin  How many times to re-bin the @f$ dN_{ch}/d\eta@f$
   * @param mask   Trigger mask 
   * @param energy Collision energy @f$ \sqrt{s}@f$ or @f$\sqrt{s_{NN}}@f$ 
   * @param title  Title to put on the plot
   * @param doHHD  Whether to do HHD comparison
   * @param doComp Whether to do comparisons 
   * 
   * @return True on success, false otherwise 
   */
  Bool_t Run(const char* file=".", 
	     Double_t vzMin=-10, Double_t vzMax=10, Int_t rebin=1,
	     Int_t mask=AliAODForwardMult::kInel, Int_t energy=900,
	     const char* title="", Bool_t doHHD=false, Bool_t doComp=false)
  {
    TString trgName;
    if    (mask & AliAODForwardMult::kInel)    trgName.Append("inel-");
    if    (mask & AliAODForwardMult::kInelGt0) trgName.Append("inelgt0-");
    if    (mask & AliAODForwardMult::kNSD)     trgName.Append("nsd-");
    if (trgName.EndsWith("-")) trgName.Remove(trgName.Length()-1);
    if (trgName.IsNull()) trgName = "unknown";
    TString outName = 
      TString::Format("dndeta_%04dGeV_%c%02d-%c%02dcm_rb%02d_%s.root",
		      energy, 
		      vzMin < 0 ? 'm' : 'p', Int_t(TMath::Abs(vzMin)+.5), 
		      vzMax < 0 ? 'm' : 'p', Int_t(TMath::Abs(vzMax)+.5), 
		      rebin, trgName.Data());
    fTitle  = title;
    fVtxMin = vzMin;
    fVtxMax = vzMax;
    if (!Open(file, outName)) return kFALSE;
    if (!Process(mask)) return kFALSE;
    if (!Finish(rebin, mask, energy,doHHD,doComp)) return kFALSE;

    return kTRUE;
  }
  //__________________________________________________________________
  /** 
   * Scan directory @a dir for AODs 
   * 
   * @param dir        Directory to scan
   * @param chain      Chain to add to
   * @param recursive  Whether to scan recusively 
   */
  void ScanDirectory(TSystemDirectory* dir, TChain* chain, bool recursive)
  {
    // gROOT->IndentLevel();
    // printf("Scanning %s ...\n", dir->GetName());
    // gROOT->IncreaseDirLevel();

    // Get list of files, and go back to old working directory
    TString oldDir(gSystem->WorkingDirectory());
    TList* files = dir->GetListOfFiles();
    gSystem->ChangeDirectory(oldDir);

    // Sort list of files and check if we should add it 
    files->Sort();
    TIter next(files);
    TSystemFile* file = 0;
    while ((file = static_cast<TSystemFile*>(next()))) {
      TString name(file->GetName());
    
      // Ignore special links 
      if (name == "." || name == "..") continue;

      // Check if this is a directory 
      if (file->IsDirectory()) { 
	if (recursive) 
	  ScanDirectory(static_cast<TSystemDirectory*>(file),chain,recursive);
	continue;
      }
    
      // If this is not a root file, ignore 
      if (!name.EndsWith(".root")) continue;

      // If this file does not contain AliESDs, ignore 
      if (!name.Contains("AliAOD")) continue;
    
      // Get the path 
      TString aod(Form("%s/%s", file->GetTitle(), name.Data()));

      // Print and add 
      // gROOT->IndentLevel();
      // printf("adding %s\n", aod.Data());
      chain->Add(aod);

    }
    // gROOT->DecreaseDirLevel();
  }
  
  //__________________________________________________________________
  /** 
   * Make a chain of AOD files 
   * 
   * @param aoddir     Directory so search
   * @param recursive  Whether to scan recusively 
   * 
   * @return The new chain
   */
  TChain*
  MakeChain(const char* aoddir, bool recursive=false)
  {
    // --- Our data chain ----------------------------------------------
    TChain* chain = new TChain("aodTree");
    
    // --- Get list of AODs --------------------------------------------
    // Open source directory, and make sure we go back to were we were 
    TString oldDir(gSystem->WorkingDirectory());
    TSystemDirectory d(aoddir, aoddir);
    ScanDirectory(&d, chain, recursive);
    
    chain->GetListOfFiles()->ls();
    if (chain->GetListOfFiles()->GetEntries() <= 0) { 
      delete chain;
      chain = 0;
    }
    return chain;
  }

  //__________________________________________________________________
  /** 
   * Open data.  Make a chain of all AOD files in the given path 
   * and below.  
   * 
   * @param path    Path to search
   * @param outname Output name 
   * 
   * @return 
   */
  Bool_t Open(const char* path, const char* outname) 
  {
    Clear("");
    
    // Get the AOD tree 
    fTree = MakeChain(path, true);
    if (!fTree) {
      Error("Init", "Couldn't get the tree");
      return kFALSE;
    }

    // Set the branch pointer 
    fTree->SetBranchAddress("Forward", &fAOD);

    // Set the branch pointer 
    if (fTree->GetBranch("ForwardMC"))
      fTree->SetBranchAddress("ForwardMC", &fMCAOD);

    // Set the branch pointer 
    if (fTree->GetBranch("primary"))
      fTree->SetBranchAddress("primary", &fPrimary);
    
    // Set branch pointer 
    if (fTree->GetBranch("Central")) 
      fTree->SetBranchAddress("Central", &fCentral);

    fOut = TFile::Open(outname, "RECREATE");
    if (!fOut) { 
      Error("Open", "Couldn't open %s", outname);
      return kFALSE;
    }
    return kTRUE;
    
  }


  //__________________________________________________________________
  /** 
   * Process the events 
   * 
   * @param vzMin  Minimum interaction point z coordinate 
   * @param vzMax  Maximum interaction point z coordinate 
   * @param mask   Trigger mask 
   *
   * @return true on success, false otherwise 
   */
  Bool_t Process(Int_t mask) 
  {
    fNTriggered       = 0;                    // # of triggered ev.
    fNWithVertex      = 0;                    // # of ev. w/vertex
    fNAccepted        = 0;                    // # of ev. used
    Int_t nAvailable  = fTree->GetEntries();  // How many entries
 
    for (int i = 0; i < nAvailable; i++) { 
      fTree->GetEntry(i);
      if (((i+1) % 100) == 0) {
	fprintf(stdout,"Event # %9d of %9d, %9d accepted so far\r", 
	       i+1, nAvailable, fNAccepted);
	fflush(stdout);
      }

      // Create sum histogram on first event - to match binning to input
      if (!fSum) 
	fSum = static_cast<TH2D*>(fAOD->GetHistogram().Clone("d2ndetadphi"));
      if (!fMCSum && fTree->GetBranch("ForwardMC")) 
	fMCSum = 
	  static_cast<TH2D*>(fMCAOD->GetHistogram().Clone("d2ndetadphiMC"));
      if (!fSumPrimary && fTree->GetBranch("primary")) 
	fSumPrimary = static_cast<TH2D*>(fPrimary->Clone("primarySum"));
      if (!fSumCentral && fTree->GetBranch("Central"))
	fSumCentral = static_cast<TH2D*>(fCentral->Clone("centralSum"));
      
      // Add contribution from this event 
      if (fSumPrimary) fSumPrimary->Add(fPrimary);
     
      fNAll++;
      if (fAOD->IsTriggerBits(AliAODForwardMult::kB)) fNB++;
      if (fAOD->IsTriggerBits(AliAODForwardMult::kA)) fNA++;
      if (fAOD->IsTriggerBits(AliAODForwardMult::kC)) fNC++;
      if (fAOD->IsTriggerBits(AliAODForwardMult::kE)) fNE++;
      
      // Other trigger/event requirements could be defined 
      if (fAOD->IsTriggerBits(AliAODForwardMult::kInel)) fNMB++;

      // fAOD->Print();
      // Other trigger/event requirements could be defined 
      if (!fAOD->IsTriggerBits(mask)) continue; 
      fNTriggered++;

      // Check if there's a vertex 
      if (!fAOD->HasIpZ()) continue; 
      fNWithVertex++;

      // Select vertex range (in centimeters) 
      if (!fAOD->InRange(fVtxMin, fVtxMax)) continue; 
      fNAccepted++;
 
      // Add contribution from this event
      fSum->Add(&(fAOD->GetHistogram()));

      // Add contribution from this event
      if (fMCSum) fMCSum->Add(&(fMCAOD->GetHistogram()));

      // Add contribution from this event
      if (fSumCentral) fSumCentral->Add(fCentral);      
    }
    printf("\n");
    // fVtxEff = Double_t(fNWithVertex)/fNTriggered;
    Int_t nGood = (fNB-fNA-fNC+2*fNE);
    fVtxEff = (Double_t(fNMB) / fNTriggered) * 
      (nGood > 0 ? Double_t(fNAccepted) / nGood : 0);
    
    Info("Process", "Total of %9d events\n"
	 "                   of these %9d are minimum bias\n"
	 "                   of these %9d has a trigger\n" 
	 "                   of these %9d has a vertex\n" 
	 "                   of these %9d was used\n"
	 "                   Triggers by type:\n"
	 "                     B   = %9d\n"
	 "                     A|C = %9d (%9d+%-9d)\n"
	 "                     E   = %9d\n"
	 "                   Implies %9d good triggers\n"
	 "                   Vertex efficiency: %f",
	 nAvailable, fNMB, fNTriggered, fNWithVertex, fNAccepted,
	 fNB, fNA+fNC, fNA, fNC, fNE, nGood, fVtxEff);

    return kTRUE;
  }
  //__________________________________________________________________
  /** 
   * Finish the stuff and draw 
   * 
   * @param rebin  How many times to rebin
   * @param energy Collision energy 
   * @param mask   Trigger mask 
   * @param doHHD  Whether to do HHD comparison
   * @param doComp Whether to do comparisons 
   * 
   * @return true on success, false otherwise 
   */
  Bool_t Finish(Int_t rebin, Int_t mask, Int_t energy, 
		Bool_t doHHD, Bool_t doComp)
  {
    fOut->cd();

    // Get acceptance normalisation from underflow bins 
    TH1D* norm   = fSum->ProjectionX("norm", 0, 1, "");
    // Project onto eta axis - _ignoring_underflow_bins_!
    TH1D* dndeta = fSum->ProjectionX("dndeta", 1, -1, "e");
    dndeta->SetTitle("ALICE Forward");
    // Normalize to the acceptance 
    dndeta->Divide(norm);
    // Scale by the vertex efficiency 
    dndeta->Scale(fVtxEff, "width");
    dndeta->SetMarkerColor(kRed+1);
    dndeta->SetMarkerStyle(20);
    dndeta->SetMarkerSize(1);
    dndeta->SetFillStyle(0);
    Rebin(dndeta, rebin);

    TH1D* dndetaMC = 0;
    if (fMCSum) { 
      // Get acceptance normalisation from underflow bins 
      norm = fMCSum->ProjectionX("norm", 0, 1, "");
      // Project onto eta axis - _ignoring_underflow_bins_!
      dndetaMC = fMCSum->ProjectionX("dndetaMC", 1, -1, "e");
      dndetaMC->SetTitle("ALICE Forward (MC)");
      // Normalize to the acceptance 
      dndetaMC->Divide(norm);
      // Scale by the vertex efficiency 
      dndetaMC->Scale(fVtxEff, "width");
      dndetaMC->SetMarkerColor(kRed+3);
      dndetaMC->SetMarkerStyle(21);
      dndetaMC->SetMarkerSize(1);
      dndetaMC->SetFillStyle(0);
      Rebin(dndetaMC, rebin);
    }

    TH1D* dndetaTruth = 0;
    if (fSumPrimary) { 
      dndetaTruth = fSumPrimary->ProjectionX("dndetaTruth", -1, -1, "e");
      //dndetaTruth->Scale(1./fNTriggered, "width");
      dndetaTruth->Scale(1./fNAll, "width");
      dndetaTruth->SetMarkerColor(kGray+3);
      dndetaTruth->SetMarkerStyle(22);
      dndetaTruth->SetMarkerSize(1);
      dndetaTruth->SetFillStyle(0);
      Rebin(dndetaTruth, rebin);
    }
    TH1D* dndetaCentral = 0;
    if (fSumCentral) { 
      dndetaCentral = fSumCentral->ProjectionX("dndetaCentral", -1, -1, "e");
      dndetaCentral->SetTitle("ALICE Central - track(let)s");
      // dndetaCentral->Scale(1./fNTriggered, "width");
      dndetaCentral->Scale(1./(fNB-fNA-fNC+2*fNE), "width");
      dndetaCentral->SetMarkerColor(kGray+3);
      dndetaCentral->SetMarkerStyle(22);
      dndetaCentral->SetMarkerSize(1);
      dndetaCentral->SetFillStyle(0);
      dndetaCentral->GetXaxis()->SetRangeUser(-1,1);
      Rebin(dndetaCentral, rebin <= 1 ? 1 : 2*(rebin/2));
      // 1 -> 1
      // 2 -> 2*2/2 -> 2*1 -> 2
      // 3 -> 2*3/2 -> 2*1 -> 2
      // 4 -> 2*4/2 -> 2*2 -> 4
      // 5 -> 2*5/2 -> 2*2 -> 4
      // 6 -> 2*6/2 -> 2*3 -> 6
    }

    DrawIt(dndeta, dndetaMC, dndetaTruth, dndetaCentral, 
	   mask, energy, doHHD, doComp);

    return kTRUE;
  }
  //__________________________________________________________________
  /** 
   */
  void DrawIt(TH1* dndeta, TH1* dndetaMC, TH1* dndetaTruth, TH1* dndetaCentral,
	      Int_t mask, Int_t energy, Bool_t doHHD, 
	      Bool_t doComp)
  {
    // --- 1st part - prepare data -----------------------------------
    TH1* dndetaSym = Symmetrice(dndeta);

    Double_t max = dndeta->GetMaximum();

    // Make our histogram stack 
    THStack* stack = new THStack("results", "Results");

    TH1* dndetaTruthSym = 0;
    if (dndetaTruth) {
      dndetaTruth->SetFillColor(kGray);
      dndetaTruth->SetFillStyle(3001);
      dndetaTruthSym = Symmetrice(dndetaTruth);
      stack->Add(dndetaTruthSym, "e5 p");
      stack->Add(dndetaTruth, "e5 p");
      Info("DrawIt", "Maximum of MC dndeta (truth): %f, was %f", 
	   dndetaTruth->GetMaximum(),dndeta->GetMaximum());
      max = TMath::Max(dndetaTruth->GetMaximum(),max);
    }

    // Get the data from HHD's analysis - if available 
    TH1* dndetaHHD    = 0;
    TH1* dndetaHHDSym = 0;
    // Info("DrawIt", "Will %sdraw HHD results", (doHHD ? "" : "not "));
    if (doHHD) {
      TString hhdName(fOut->GetName());
      hhdName.ReplaceAll("dndeta", "hhd");    
      dndetaHHD    = GetHHD(hhdName.Data(), mask & AliAODForwardMult::kNSD);
      dndetaHHDSym = 0;
      if (dndetaHHD) { 
	// Symmetrice and add to stack 
	dndetaHHD->SetTitle("ALICE Forward (HHD)");
	dndetaHHDSym = Symmetrice(dndetaHHD);
	stack->Add(dndetaHHDSym);
	stack->Add(dndetaHHD);
	max = TMath::Max(dndetaHHD->GetMaximum(),max);
      }
      else 
	Warning("DrawIt", "No HHD data found");
      fOut->cd();
    }

    // If we have MC analysed data, plot it 
    if (dndetaMC) { 
      TH1* dndetaMCSym = Symmetrice(dndetaMC);
      stack->Add(dndetaMCSym);
      stack->Add(dndetaMC);
      max = TMath::Max(dndetaMC->GetMaximum(),max);
    }

    // Add the analysis results to the list 
    stack->Add(dndetaSym);
    stack->Add(dndeta);

    // If we have central region data, add that 
    if (dndetaCentral) { 
      stack->Add(dndetaCentral);
      max = TMath::Max(dndetaCentral->GetMaximum(),max);
    }


    // Get graph of 'other' data - e.g., UA5, CMS, ... - and check if
    // there's any graphs.  Set the pad division based on that.
    // Info("DrawIt", "Will %sdraw other results", (doComp ? "" : "not "));
    TMultiGraph* other    = (doComp ? GetData(energy, mask) : 0);
    THStack*     ratios   = MakeRatios(dndeta,      dndetaSym, 
				       dndetaHHD,   dndetaHHDSym, 
				       dndetaTruth, dndetaTruthSym, 
				       dndetaCentral,
				       other);

    // Check if we have ratios 

    // --- 2nd part - Draw in canvas ---------------------------------
    // Make a canvas 
    gStyle->SetOptTitle(0);
    gStyle->SetTitleFont(132, "xyz");
    gStyle->SetLabelFont(132, "xyz");
    TCanvas* c = new TCanvas("c", "C", 900, 800);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(0);
    c->cd();

    Double_t yd = (ratios ? 0.3 : 0);

    // Make a sub-pad for the result itself
    TPad* p1 = new TPad("p1", "p1", 0, yd, 1.0, 1.0, 0, 0);
    p1->SetTopMargin(0.05);
    p1->SetBottomMargin(ratios ? 0.001 : 0.1);
    p1->SetRightMargin(0.05);
    p1->SetGridx();
    p1->SetTicks(1,1);
    p1->Draw();
    p1->cd();

    // Fix the apperance of the stack and redraw. 
    if (other) {
      TGraphAsymmErrors* o      = 0;
      TIter              nextG(other->GetListOfGraphs());
      while ((o = static_cast<TGraphAsymmErrors*>(nextG()))) {
	Double_t gmax = TMath::MaxElement(o->GetN(), o->GetY());
	//Info("DrawIt", "Maximum of %s is %f, was %f", o->GetName(),gmax,max);
	max = TMath::Max(max, gmax);
      }
    }
    max *= 1.1;
    // Info("DrawIt", "Setting maximum to %f", max);
    stack->SetMinimum(ratios ? -0.1 : 0);
    stack->SetMaximum(max);
    FixAxis(stack, 1/(1-yd)/1.5, 
	    "#frac{1}{#it{N}} #frac{#it{dN}_{ch}}{#it{d#eta}}");
    p1->Clear();
    stack->DrawClone("nostack e1");
    if (dndetaCentral) 
      dndetaCentral->GetXaxis()->SetRangeUser(-1,1);
      
    // Draw other data 
    if (other) {
      TGraphAsymmErrors* o      = 0;
      TIter              nextG(other->GetListOfGraphs());
      while ((o = static_cast<TGraphAsymmErrors*>(nextG())))
	o->Draw("same p");
    }

    // Make a legend in the main result pad 
    TString trigs(AliAODForwardMult::GetTriggerString(mask));
    TLegend* l = p1->BuildLegend(.15,p1->GetBottomMargin()+.01,.90,.35);
    l->SetNColumns(2);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(132);
    p1->cd();

    // Put a title on top 
    TLatex* tit = new TLatex(0.10, 0.95, fTitle.Data());
    tit->SetNDC();
    tit->SetTextFont(132);
    tit->SetTextSize(0.05);
    tit->Draw();

    // Put a nice label in the plot 
    TString eS;
    if (energy > 1000) eS = Form("%4.2fTeV", float(energy)/1000);
    else               eS = Form("%3dGeV", energy);
    TLatex* tt = new TLatex(.93, .93, 
			    Form("#sqrt{s}=%s, %s", eS.Data(),
				 AliAODForwardMult::GetTriggerString(mask)));
    tt->SetNDC();
    tt->SetTextFont(132);
    tt->SetTextAlign(33);
    tt->Draw();

    // Put number of accepted events on the plot 
    TLatex* et = new TLatex(.93, .83, Form("%d events", fNAccepted));
    et->SetNDC();
    et->SetTextFont(132);
    et->SetTextAlign(33);
    et->Draw();

    // Put number of accepted events on the plot 
    TLatex* vt = new TLatex(.93, .88, 
			    Form("v_{z}#in[%+5.1f,%+5.1f]cm",fVtxMin,fVtxMax));
    vt->SetNDC();
    vt->SetTextFont(132);
    vt->SetTextAlign(33);
    vt->Draw();

    // Mark the plot as preliminary
    TLatex* pt = new TLatex(.12, .93, "Preliminary");
    pt->SetNDC();
    pt->SetTextFont(22);
    pt->SetTextSize(0.07);
    pt->SetTextColor(kRed+1);
    pt->SetTextAlign(13);
    pt->Draw();
    c->cd();

    // If we do not have the ratios, fix up the display 
    // p1->SetPad(0, 0, 1, 1);
    // p1->SetBottomMargin(0.1);
    // l->SetY1(0.11);
    // stack->SetMinimum(0);
    // FixAxis(stack, (1-yd)/1,  "#frac{1}{N} #frac{dN_{ch}}{#eta}",10,false);
    if (ratios) {
      // If we do have the ratios, then make a new pad and draw the 
      // ratios there 
      c->cd();
      TPad* p2 = new TPad("p2", "p2", 0, 0.0, 1.0, yd, 0, 0);
      p2->SetTopMargin(0.001);
      p2->SetRightMargin(0.05);
      p2->SetBottomMargin(1/yd * 0.07);
      p2->SetGridx();
      p2->SetTicks(1,1);
      p2->Draw();
      p2->cd();

      // Fix up axis 
      FixAxis(ratios, 1/yd/1.5, "Ratios", 5);

      // Fix up y range and redraw 
      ratios->SetMinimum(.58);
      ratios->SetMaximum(1.22);
      p2->Clear();
      ratios->DrawClone("nostack e1");
      
      // Make a legend 
      TLegend* l2 = p2->BuildLegend(.15,p2->GetBottomMargin()+.01,.9,.6);
      l2->SetNColumns(2);
      l2->SetFillColor(0);
      l2->SetFillStyle(0);
      l2->SetBorderSize(0);
      l2->SetTextFont(132);

      // Make a nice band from 0.9 to 1.1 
      TGraphErrors* band = new TGraphErrors(2);
      band->SetPoint(0, dndetaSym->GetXaxis()->GetXmin(), 1);
      band->SetPoint(1, dndeta->GetXaxis()->GetXmax(), 1);
      band->SetPointError(0, 0, .1);
      band->SetPointError(1, 0, .1);
      band->SetFillColor(kYellow+2);
      band->SetFillStyle(3002);
      band->SetLineStyle(2);
      band->SetLineWidth(1);
      band->Draw("3 same");
      band->DrawClone("X L same");

      // Replot the ratios on top 
      ratios->DrawClone("nostack e1 same");

      c->cd();
    }
    
    // Plot to disk
    TString imgName(fOut->GetName());
    imgName.ReplaceAll(".root", ".png");
    c->SaveAs(imgName.Data());
    imgName.ReplaceAll(".png", ".C");
    c->SaveAs(imgName.Data());
    
    fOut->cd();
    stack->Write();
    if (other)  other->Write();
    if (ratios) ratios->Write();

    // Close our file 
    fOut->Close();
  }
  //__________________________________________________________________
  /** 
   * Get the result from previous analysis code 
   * 
   * @param fn  File to open 
   * @param nsd Whether this is NSD
   * 
   * @return null or result of previous analysis code 
   */
  TH1* GetHHD(const char* fn="fmd_dNdeta_mult.root", bool nsd=false) 
  {
    TDirectory* savdir = gDirectory;
    if (gSystem->AccessPathName(fn)) { 
      Warning("GetHHD", "Output of HHD analysis (%s) not available", fn);
      return 0;
    }
    TFile* file = TFile::Open(fn, "READ");
    if (!file) { 
      Warning("GetHHD", "couldn't open HHD file %s", fn);
      savdir->cd();
      return 0;
    }
    TString hist(Form("dNdeta_dNdeta%s", (nsd ? "NSD" : "")));
    TH1* h = static_cast<TH1*>(file->Get(hist.Data()));
    if (!h) { 
      Warning("GetHHD", "Couldn't find HHD histogram %s in %s", 
	      hist.Data(), fn);
      file->Close();
      savdir->cd();
      return 0;
    }
    TH1* r = static_cast<TH1*>(h->Clone("dndeta_hhd"));
    r->SetTitle("1/N dN_{ch}/d#eta (HHD)");
    r->SetFillStyle(0);
    r->SetFillColor(0);
    r->SetMarkerStyle(21);
    r->SetMarkerColor(kPink+1);
    r->SetDirectory(savdir);

    file->Close();
    savdir->cd();
    return r;
  }
  //__________________________________________________________________
  /** 
   */ 
  THStack* MakeRatios(const TH1* dndeta, const TH1* dndetaSym, 
		      const TH1* hhd,    const TH1* hhdSym, 
		      const TH1* mc,     const TH1* mcSym,
		      const TH1* central, 
		      TMultiGraph* other) const 
  {
    // If we have 'other' data, then do the ratio of the results to that
    Bool_t hasOther = (other && other->GetListOfGraphs() && 
		       other->GetListOfGraphs()->GetEntries() > 0);
    Bool_t hasHhd   = (hhd && hhdSym);
    if (!hasOther && !hasHhd && !mc && !mcSym) return 0;

    THStack* ratios = new THStack("ratios", "Ratios");
    if (hasOther) {
      TGraphAsymmErrors* ua5    = 0;
      TGraphAsymmErrors* alice  = 0;
      TGraphAsymmErrors* cms    = 0;
      TGraphAsymmErrors* o      = 0;
      TIter              nextG(other->GetListOfGraphs());
      while ((o = static_cast<TGraphAsymmErrors*>(nextG()))) {
	ratios->Add(Ratio(dndeta,    o));
	ratios->Add(Ratio(dndetaSym, o));
	ratios->Add(Ratio(hhd,       o));
	ratios->Add(Ratio(hhdSym,    o));
	ratios->Add(Ratio(central,   o));
	TString oName(o->GetName());
	oName.ToLower();
	if (oName.Contains("ua5"))   ua5 = o;
	if (oName.Contains("alice")) alice = o;
	if (oName.Contains("cms"))   cms = o;
      }
      if (ua5 && alice) 
	ratios->Add(Ratio(alice, ua5));
      if (cms && alice) 
	ratios->Add(Ratio(alice, cms));
    }
    

    // If we have data from HHD's analysis, then do the ratio of 
    // our result to that data. 
    if (hasHhd) { 
      ratios->Add(Ratio(dndeta,    hhd));
      ratios->Add(Ratio(dndetaSym, hhdSym));
    }

    // Do comparison to MC 
    if (mc) { 
      ratios->Add(Ratio(dndeta,    mc));
      ratios->Add(Ratio(dndetaSym, mcSym));
    }

    // Check if we have ratios 
    Bool_t   hasRatios = (ratios->GetHists() && 
			  (ratios->GetHists()->GetEntries() > 0));
#if 0
    Info("MakeRatios", "Got a total of %d ratios", !hasRatios ? 0 :
	 ratios->GetHists()->GetEntries());
#endif

    if (!hasRatios) { delete ratios; ratios = 0; }
    return ratios;
  }

  //__________________________________________________________________
  /** 
   * Fix the apperance of the axis in a stack 
   * 
   * @param stack  stack of histogram
   * @param s      Scaling factor 
   * @param ytitle Y axis title 
   * @param force  Whether to draw the stack first or not 
   * @param ynDiv  Divisions on Y axis 
   */
  void FixAxis(THStack* stack, Double_t s, const char* ytitle, 
	       Int_t ynDiv=210, Bool_t force=true) 
  {
    if (force) stack->Draw("nostack e1");

    TH1* h = stack->GetHistogram();
    if (!h) return;

    h->SetXTitle("#eta");
    h->SetYTitle(ytitle);
    TAxis* xa = h->GetXaxis();
    TAxis* ya = h->GetYaxis();
    if (xa) { 
      xa->SetTitle("#eta");
      // xa->SetTicks("+-");
      xa->SetTitleSize(s*xa->GetTitleSize());
      xa->SetLabelSize(s*xa->GetLabelSize());
      xa->SetTickLength(s*xa->GetTickLength());
    }
    if (ya) { 
      ya->SetTitle(ytitle);
      ya->SetDecimals();
      // ya->SetTicks("+-");
      ya->SetNdivisions(ynDiv);
      ya->SetTitleSize(s*ya->GetTitleSize());
      ya->SetLabelSize(s*ya->GetLabelSize());
    }      
  }
  //__________________________________________________________________
  /** 
   * Compute the ratio of @a h to @a g.  @a g is evaluated at the bin
   * centers of @a h
   * 
   * @param h  Numerator 
   * @param g  Divisor 
   * 
   * @return h/g 
   */
  TH1* Ratio(const TH1* h, const TGraph* g) const 
  {
    if (!h || !g) return 0;

    TH1* ret = static_cast<TH1*>(h->Clone("tmp"));
    ret->SetName(Form("%s_over_%s", h->GetName(), g->GetName()));
    ret->SetTitle(Form("%s / %s", h->GetTitle(), g->GetTitle()));
    ret->Reset();
    ret->SetMarkerStyle(g->GetMarkerStyle());
    ret->SetMarkerColor(h->GetMarkerColor());
    ret->SetMarkerSize(0.9*g->GetMarkerSize());
    Double_t xlow  = g->GetX()[0];
    Double_t xhigh = g->GetX()[g->GetN()-1];
    if (xlow > xhigh) { Double_t t = xhigh; xhigh = xlow; xlow = t; }

    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t c = h->GetBinContent(i);
      if (c <= 0) continue;

      Double_t x = h->GetBinCenter(i);
      if (x < xlow || x > xhigh) continue; 

      Double_t f = g->Eval(x);
      if (f <= 0) continue; 

      ret->SetBinContent(i, c / f);
      ret->SetBinError(i, h->GetBinError(i) / f);
    }
    if (ret->GetEntries() <= 0) { delete ret; ret = 0; }
    return ret;
  }
  //__________________________________________________________________
  /** 
   * Make the ratio of h1 to h2 
   * 
   * @param h1 First histogram (numerator) 
   * @param h2 Second histogram (denominator)
   * 
   * @return h1 / h2
   */
  TH1* Ratio(const TH1* h1, const TH1* h2) const
  {
    TH1* t1 = static_cast<TH1*>(h1->Clone(Form("%s_%s", 
					       h1->GetName(), 
					       h2->GetName())));
    t1->SetTitle(Form("%s / %s", h1->GetTitle(), h2->GetTitle()));
    t1->Divide(h2);
    t1->SetMarkerColor(h1->GetMarkerColor());
    t1->SetMarkerStyle(h2->GetMarkerStyle());
    return t1;
  }
  //__________________________________________________________________
  /** 
   * Calculate the ratio of two graphs - g1 / g2
   * 
   * @param g1 Numerator 
   * @param g2 Denominator
   * 
   * @return g1 / g2 in a histogram 
   */
  TH1* Ratio(const TGraphAsymmErrors* g1, const TGraphAsymmErrors* g2) const
  {
    Int_t    nBins = g1->GetN();
    TArrayF  bins(nBins+1);
    Double_t dx = 0;
    for (Int_t i = 0; i < nBins; i++) { 
      Double_t x   = g1->GetX()[i];
      Double_t exl = g1->GetEXlow()[i];
      Double_t exh = g1->GetEXhigh()[i];
      bins.fArray[i]   = x-exl;
      bins.fArray[i+1] = x+exh;
      Double_t dxi = exh+exl;
      if (i == 0) dx  = dxi;
      else if (dxi != dx) dx = 0;
    }
    TString name(Form("%s_%s", g1->GetName(), g2->GetName()));
    TString title(Form("%s / %s", g1->GetTitle(), g2->GetTitle()));
    TH1* h = 0;
    if (dx != 0) {
      h = new TH1F(name.Data(), title.Data(), nBins, bins[0], bins[nBins]);
    }
    else {
      h = new TH1F(name.Data(), title.Data(), nBins, bins.fArray);
    }
    h->SetMarkerStyle(g2->GetMarkerStyle());
    h->SetMarkerColor(g1->GetMarkerColor());
    h->SetMarkerSize(0.9*g2->GetMarkerSize());

    for (Int_t i = 0; i < nBins; i++) { 
      Double_t x  = g1->GetX()[i];
      Double_t c1 = g1->GetY()[i];
      Double_t e1 = g1->GetErrorY(i);
      Double_t c2 = g2->Eval(x);
      
      h->SetBinContent(i+1, c1 / c2);
      h->SetBinError(i+1, e1 / c2);
    }
    return h;
  }
  //__________________________________________________________________
  /** 
   * Make an extension of @a h to make it symmetric about 0 
   * 
   * @param h Histogram to symmertrice 
   * 
   * @return Symmetric extension of @a h 
   */
  TH1* Symmetrice(const TH1* h) const
  {
    fOut->cd();

    Int_t nBins = h->GetNbinsX();
    TH1*  s     = new TH1D(Form("%s_mirror", h->GetName()),
			   Form("%s (mirrored)", h->GetTitle()), 
			   nBins, 
			   -h->GetXaxis()->GetXmax(), 
			   -h->GetXaxis()->GetXmin());
    s->SetMarkerColor(h->GetMarkerColor());
    s->SetMarkerSize(h->GetMarkerSize());
    s->SetMarkerStyle(h->GetMarkerStyle()+4);
    s->SetFillColor(h->GetFillColor());
    s->SetFillStyle(h->GetFillStyle());
    // s->SetDirectory(0);

    // Find the first and last bin with data 
    Int_t first = nBins+1;
    Int_t last  = 0;
    for (Int_t i = 1; i <= nBins; i++) { 
      if (h->GetBinContent(i) <= 0) continue; 
      first = TMath::Min(first, i);
      last  = TMath::Max(last,  i);
    }
    
    Double_t xfirst = h->GetBinCenter(first-1);
    Int_t    f1     = h->GetXaxis()->FindBin(-xfirst);
    Int_t    l2     = s->GetXaxis()->FindBin(xfirst);
    for (Int_t i = f1, j=l2; i <= last; i++,j--) { 
      s->SetBinContent(j, h->GetBinContent(i));
      s->SetBinError(j, h->GetBinError(i));
    }
    // Fill in overlap bin 
    s->SetBinContent(l2+1, h->GetBinContent(first));
    s->SetBinError(l2+1, h->GetBinError(first));
    return s;
  }
  //__________________________________________________________________
  /** 
   * Rebin a histogram 
   * 
   * @param h     Histogram to rebin
   * @param rebin Rebinning factor 
   * 
   * @return 
   */
  virtual void Rebin(TH1* h, Int_t rebin) const
  { 
    if (rebin <= 1) return;

    Int_t nBins = h->GetNbinsX();
    if(nBins % rebin != 0) {
      Warning("Rebin", "Rebin factor %d is not a devisor of current number "
	      "of bins %d in the histogram %s", rebin, nBins, h->GetName());
      return;
    }
    
    // Make a copy 
    TH1* tmp = static_cast<TH1*>(h->Clone("tmp"));
    tmp->Rebin(rebin);
    tmp->SetDirectory(0);

    // The new number of bins 
    Int_t nBinsNew = nBins / rebin;
    for(Int_t i = 1;i<= nBinsNew; i++) {
      Double_t content = 0;
      Double_t sumw    = 0;
      Double_t wsum    = 0;
      Int_t    nbins   = 0;
      for(Int_t j = 1; j<=rebin;j++) {
	Int_t    bin = (i-1)*rebin + j;
	Double_t c   =  h->GetBinContent(bin);

	if (c <= 0) continue;

	if (fCutEdges) {
	  if (h->GetBinContent(bin+1)<=0 || 
	      h->GetBinContent(bin-1)) {
	    Warning("Rebin", "removing bin %d=%f of %s (%d=%f,%d=%f)", 
		    bin, c, h->GetName(), 
		    bin+1, h->GetBinContent(bin+1), 
		    bin-1, h->GetBinContent(bin-1));
	    continue;
	  }	
	}
	Double_t e =  h->GetBinError(bin);
	Double_t w =  1 / (e*e); // 1/c/c
	content    += c;
	sumw       += w;
	wsum       += w * c;
	nbins++;
      }
      
      if(content > 0 && nbins > 0) {
	tmp->SetBinContent(i, wsum / sumw);
	tmp->SetBinError(i,1./TMath::Sqrt(sumw));
      }
    }

    // Finally, rebin the histogram, and set new content
    h->Rebin(rebin);
    h->Reset();
    for(Int_t i = 1; i<= nBinsNew; i++) {
      h->SetBinContent(i,tmp->GetBinContent(i));
      h->SetBinError(i,  tmp->GetBinError(i));
    }
    
    delete tmp;
  }
};

//____________________________________________________________________
//
// EOF
//

