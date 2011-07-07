/**
 * A task to do a comparison between tracklets and clusers in the SPD
 * 
 * Since the class SPDComparisonTask derives from a compiled class
 * (AliAnalysisTaskSE) we need to compile that code.  The script will,
 * when executed in the AliROOT prompt load it self again and byte
 * compile it with the preprocessor flag BUILD defined.  \
 *
 * The preprocessor define BUILD is not defined at first, when the
 * script is loaded using 
 * 
 * @verbatim 
 *   Root> .x SPDComparison.C 
 * @endverbatim 
 * 
 * which means that CINT will only see the function SPDComparison.
 * In that function, we define the BUILD preprocessor symbol 
 *
 * @code 
 *   gSystem->AddIncludePath("-DBUILD=1 ...");
 * @endcode 
 * 
 * and then ACLic compile ourselves 
 * 
 * @code 
 *   gROOT->LoadMacro("./SPDComparison.C++");
 * @endcode 
 * 
 * But since BUILD is now defined, it means that ACLic will only see
 * the class and not the function SPDComparison. 
 * 
 * This trick hinges on that when you initially load the script and
 * when it is done inside the script it is done using two distinct
 * paths - otherwise ROOT will try to unload the script first, and
 * that fails.   
 * 
 */

#ifdef BUILD
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliLog.h>
#include <AliMultiplicity.h>
#include "AliFMDEventInspector.h"
#include "AliAODForwardMult.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TList.h>
#include <TObjArray.h>
#include <TMath.h>


/**
 * @class SPDComparisonTask 
 *
 * Task to do forward/backward correlations using FMD and SPD
 * (cluster) data. 
 * 
 * The class contains the sub-structure Bin that represents an
 * @f$\eta@f$ range.  One can add (possibly overlapping) @f$\eta@f$
 * ranges by calling the member function AddBin 
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
class SPDComparisonTask : public AliAnalysisTaskSE
{
public:
  /**
   * A single vertex bin - contains summed histogram of clusters and
   * tracklets 
   * 
   */
  struct VtxBin : public TNamed
  {
    /** 
     * Static function to get a bin name 
     * 
     * @param low   Low limit 
     * @param high  High limit 
     * 
     * @return Pointer to static array
     */
    static const char* BinName(Double_t low, Double_t high)
    {
      static TString n;
      n = (Form("vtx%+05.1f_%+05.1f", low, high));
      n.ReplaceAll("+", "p");
      n.ReplaceAll("-", "m");
      n.ReplaceAll(".", "d");
      return n.Data();
    }
    /** 
     * Constructor
     * 
     * @param low   Low limit 
     * @param high  High limit 
     */
    VtxBin(Double_t low, Double_t high)
      : TNamed(BinName(low,high), ""),
	fClusters(0), 
	fTracklets(0)
    {
    }
    /** 
     * Constructor
     */
    VtxBin() : TNamed(), fClusters(0), fTracklets(0) {}
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    VtxBin(const VtxBin& o) : 
      TNamed(o), fClusters(o.fClusters), fTracklets(o.fTracklets) 
    {}
    /** 
     * Assignment operator 
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this.
     */
    VtxBin& operator=(const VtxBin& o) 
    {
      TNamed::operator=(o);
      fClusters = o.fClusters;
      fTracklets = o.fTracklets;
      return *this;
    }
    /** 
     * Define outputs
     * 
     * @param l    Output list
     * @param eta  Eta axis to use 
     */
    void DefineOutput(TList* l, const TAxis& eta)
    {
      TList* ll = new TList;
      ll->SetName(GetName());
      ll->SetOwner();
      l->Add(ll);

      fClusters = new TH1D("clusters", "Clusters", 
			   eta.GetNbins(), 
			   eta.GetXmin(), 
			   eta.GetXmax());
      fClusters->SetXTitle("#eta");
      fClusters->SetYTitle("dN/d#eta");
      fClusters->SetDirectory(0);
      fClusters->Sumw2();
      ll->Add(fClusters);

      fTracklets = static_cast<TH1D*>(fClusters->Clone("tracklets"));
      fTracklets->SetTitle("Tracklets");
      fTracklets->SetDirectory(0);
      ll->Add(fTracklets);
    }
    /** 
     * Process the input 
     * 
     * @param spdmult  Input 
     */
    void Process(const AliMultiplicity* spdmult)
    {
      //Filling clusters in layer 1 used for tracklets...
      for(Int_t j = 0; j< spdmult->GetNumberOfTracklets();j++) {
	Double_t eta = spdmult->GetEta(j);
	fClusters->Fill(eta);
	fTracklets->Fill(eta);
      }

      //...and then the unused ones in layer 1 
      for(Int_t j = 0; j< spdmult->GetNumberOfSingleClusters();j++) {
	Double_t eta = -TMath::Log(TMath::Tan(spdmult->GetThetaSingle(j)/2.));
	fClusters->Fill(eta);
      }
    }
    /** 
     * End of job processing
     * 
     * @param l        Output list
     * @param nEvents  Number of events to normalise to 
     */
    void Finish(TList* l, Int_t nEvents)
    {
      TList* ll  = static_cast<TList*>(l->FindObject(GetName()));
      fClusters  = static_cast<TH1D*>(ll->FindObject("clusters"));
      fTracklets = static_cast<TH1D*>(ll->FindObject("tracklets"));
      fClusters->Scale(1. / nEvents, "width");
      fTracklets->Scale(1. / nEvents, "width");
      TH1* ratio = static_cast<TH1D*>(fClusters->Clone("ratio"));
      ratio->SetDirectory(0);
      ratio->SetTitle("Clusters/Tracklets");
      ratio->Divide(fTracklets);
      ll->Add(ratio);
    }
    TH1D* fClusters;  // Cluster histogram
    TH1D* fTracklets; // Tracklet histogram 
  };
    
  //__________________________________________________________________
  /** 
   * Default constructor - do not use 
   */
  SPDComparisonTask() 
    : AliAnalysisTaskSE(), 
      fInspector(),
      fBins(0),
      fEvents(0),
      fList(0), 
      fFirstEvent(true), 
      fVertexAxis(20, -20, 20)
  {}
  /** 
   * Constructor with a single argument. 
   */
  SPDComparisonTask(const char*) 
    : AliAnalysisTaskSE("spdComparision"), 
      fInspector("eventInspector"),
      fBins(0),
      fEvents(0),
      fList(0), 
      fFirstEvent(true), 
      fVertexAxis(20, -20, 20)
  {
    // Declare our output container 
    DefineOutput(1,TList::Class());
  }
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  SPDComparisonTask(const SPDComparisonTask& o) 
    : AliAnalysisTaskSE(o),
      fInspector(o.fInspector),
      fBins(o.fBins),
      fEvents(o.fEvents),
      fList(o.fList), 
      fFirstEvent(o.fFirstEvent), 
      fVertexAxis(20, -20, 20)
  {
    SetVertexAxis(o.fVertexAxis);
  }
  /** 
   * Assigment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  SPDComparisonTask& operator=(const SPDComparisonTask& o) 
  { 
    AliAnalysisTaskSE::operator=(o);
    fInspector = o.fInspector;
    fEvents    = o.fEvents;
    fList      = o.fList;
    fFirstEvent = o.fFirstEvent;
    SetVertexAxis(o.fVertexAxis);
    return *this; 
  }
  /** 
   * Destructor 
   */
  ~SPDComparisonTask() {}
  /** 
   * Set the vertex axis to use 
   * 
   * @param a Axis to set from 
   */
  void SetVertexAxis(const TAxis& a)
  {
    SetVertexAxis(a.GetNbins(), 
		  a.GetXmin(),
		  a.GetXmax());
  }
  /** 
   * Set the vertex axis to use 
   * 
   * @param n     Number of bins
   * @param xmin  Least @f$v_z@f$
   * @param xmax  Most @f$v_z@f$
   */
  void SetVertexAxis(Int_t n, Double_t xmin, Double_t xmax)
  {
    fVertexAxis.Set(n, xmin, xmax);
  }
  /** 
   * Create output objects 
   * 
   */
  void UserCreateOutputObjects()
  {
    fList = new TList;
    fList->SetName(GetName());
    fList->SetOwner();


    fEvents = new TH1D("events", "Events", 
		       fVertexAxis.GetNbins(), 
		       fVertexAxis.GetXmin(), 
		       fVertexAxis.GetXmax());
    fEvents->SetDirectory(0);
    fEvents->SetXTitle("v_{z} [cm]");
    fList->Add(fEvents);

    TAxis& vtxAxis = fVertexAxis;
    fBins = new TObjArray(vtxAxis.GetNbins(),1);
    
    TAxis etaAxis(120, -3, 3);
    for (Int_t i = 1; i <= vtxAxis.GetNbins(); i++) {
      VtxBin* bin = new VtxBin(vtxAxis.GetBinLowEdge(i), 
			       vtxAxis.GetBinUpEdge(i));
      bin->DefineOutput(fList, etaAxis);
      fBins->AddAt(bin,i);
    }


    fInspector.DefineOutput(fList);
    fInspector.Init(*(fEvents->GetXaxis()));
    
    PostData(1, fList);
  }
  /** 
   * Process one event
   * 
   * @param option Not used 
   */    
  void UserExec(Option_t* option="")
  {
    // Get the input data - ESD event
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!esd) { 
      AliWarning("No ESD event found for input event");
      return;
    }
    
    //--- Read run information -----------------------------------------
    if (fFirstEvent && esd->GetESDRun()) {
      fInspector.ReadRunDetails(esd);
      
      AliInfo(Form("Initializing with parameters from the ESD:\n"
		   "         AliESDEvent::GetBeamEnergy()   ->%f\n"
		   "         AliESDEvent::GetBeamType()     ->%s\n"
		   "         AliESDEvent::GetCurrentL3()    ->%f\n"
		   "         AliESDEvent::GetMagneticField()->%f\n"
		   "         AliESDEvent::GetRunNumber()    ->%d\n",
		   esd->GetBeamEnergy(),
		   esd->GetBeamType(),
		   esd->GetCurrentL3(),
		   esd->GetMagneticField(),
		   esd->GetRunNumber()));
      
      Print();
      fFirstEvent = false;
    }
    
    // Some variables 
    UInt_t   triggers; // Trigger bits
    Bool_t   lowFlux;  // Low flux flag
    UShort_t iVz;      // Vertex bin from ESD
    Double_t vZ;       // Z coordinate from ESD
    Double_t cent;     // Centrality 
    UShort_t nClusters;// Number of SPD clusters 
    // Process the data 
    UInt_t retESD = fInspector.Process(esd, triggers, lowFlux, iVz, vZ, 
				       cent, nClusters);
    Bool_t isInel   = triggers & AliAODForwardMult::kB;
    Bool_t hasVtx   = retESD == AliFMDEventInspector::kOk;

    if (!isInel || !hasVtx) return;
    fEvents->Fill(vZ);

    const AliMultiplicity* spdmult = esd->GetMultiplicity();

    VtxBin* bin = static_cast<VtxBin*>(fBins->At(iVz));
    if (!bin) { 
      AliError(Form("No bin @ %d (%fcm)", iVz, vZ));
      return;
    }
    bin->Process(spdmult);
    
    // Post data to output container 
    PostData(1,fList);
  }
  /** 
   * Called at the end of the processing 
   * 
   * @param option 
   */
  void Terminate(Option_t* option="")
  {
    fList = dynamic_cast<TList*>(GetOutputData(1));
    if (!fList) {
      AliError("No output list defined");
      return;
    }

    fEvents    = static_cast<TH1D*>(fList->FindObject("events"));

    Int_t nEvents = fEvents->GetEntries();
    AliInfo(Form("Got a total of %d events", nEvents));

    TH1* clusters = 0;
    TH1* tracklets = 0;
    TIter next(fBins);
    VtxBin* bin = 0;
    Int_t i = 1;
    while ((bin = static_cast<VtxBin*>(next()))) { 
      bin->Finish(fList, fEvents->GetBinContent(i++));
      if (!clusters) clusters = static_cast<TH1D*>(bin->fClusters->Clone());
      else clusters->Add(bin->fClusters);
      if (!tracklets) tracklets = static_cast<TH1D*>(bin->fTracklets->Clone());
      else tracklets->Add(bin->fTracklets);
    }
    clusters->SetDirectory(0);
    tracklets->SetDirectory(0);
    clusters->Scale(1. / i);
    tracklets->Scale(1. / i);
    

    TH1D* ratio = static_cast<TH1D*>(clusters->Clone("ratio"));
    ratio->SetDirectory(0);
    ratio->SetTitle("Clusters/Tracklets");
    ratio->Divide(tracklets);

    fList->Add(clusters);
    fList->Add(tracklets);
    fList->Add(ratio);
  }
protected: 
  AliFMDEventInspector fInspector; // Inspector
  TObjArray*           fBins;      // Vertex bins 
  TH1D*                fEvents;    // Events 
  TList*               fList;      // Output list
  Bool_t               fFirstEvent;// First event flag 
  TAxis                fVertexAxis;

  ClassDef(SPDComparisonTask,1); // BF analysis
};
#else

//====================================================================
/** 
 * Run the analysis 
 * 
 * @param esddir  Input directory 
 * @param nEvents Number of events, negative means all 
 */
void
SPDComparison(const char* esddir, Int_t nEvents=-1)
{
  // --- Libraries to load -------------------------------------------
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  // --- Our data chain ----------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/MakeChain.C");
  TChain* chain = MakeChain("ESD", esddir, true);
  // If 0 or less events is select, choose all 
  if (nEvents <= 0) nEvents = chain->GetEntries();
  
  // --- Manager -----------------------------------------------------
  AliAnalysisManager* mgr = new AliAnalysisManager("FB", "FB train");
  AliAnalysisManager::SetCommonFileName("spd_comps.root");

  // --- AOD input handler -------------------------------------------
  AliESDInputHandler *inputHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(inputHandler);      
   
  // --- compile our code --------------------------------------------
  gSystem->AddIncludePath("-I${ALICE_ROOT}/PWG2/FORWARD/analysis2 "
                          "-I${ALICE_ROOT}/ANALYSIS "
                          "-I${ALICE_ROOT}/include -DBUILD=1");
  gROOT->LoadMacro("./SPDComparison.C++g");
  
  // --- Make our object ---------------------------------------------
  SPDComparisonTask* task = new SPDComparisonTask("SPD_COMP");
  task->SetVertexAxis(10, -10, 10);
  mgr->AddTask(task);

  // --- create containers for input/output --------------------------
  AliAnalysisDataContainer *sums = 
    mgr->CreateContainer("spdComp", TList::Class(), 
                         AliAnalysisManager::kOutputContainer, 
                         AliAnalysisManager::GetCommonFileName());
  // --- connect input/output ----------------------------------------
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, sums);

  // --- Run the analysis --------------------------------------------
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("SPDComparison", "Failed to initialize analysis train!");
    return;
  }
  // Some informative output 
  mgr->PrintStatus();
  // mgr->SetDebugLevel(3);
  if (mgr->GetDebugLevel() < 1) mgr->SetUseProgressBar(kTRUE,100);

  // Run the train 
  t.Start();
  Printf("=== RUNNING ANALYSIS on %9d events ==================", nEvents);
  mgr->StartAnalysis("local", chain, nEvents);
  t.Stop();
  t.Print();

  DrawSPDComparison();
}

//====================================================================
/** 
 * Draw results
 * 
 * @param filename 
 */
void
DrawSPDComparison(const char* filename="spd_comps.root")
{
  // --- Open the file -----------------------------------------------
  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("DrawSPDComparison", "Failed to open file %s", filename);
    return;
  }
  
  // --- Get our list ------------------------------------------------
  TList* spd = static_cast<TList*>(file->Get("spdComp"));
  if (!spd) { 
    Error("DrawSPDComparison", "Failed to get list SPD_COMP from file %s",
	  filename);
    return;
  }

  // --- Loop over list and get correlation plots --------------------
  TH1* clusters  = static_cast<TH1*>(spd->FindObject("clusters"));
  TH1* tracklets = static_cast<TH1*>(spd->FindObject("tracklets"));
  TH1* ratio     = static_cast<TH1*>(spd->FindObject("ratio"));
  TH1* events    = static_cast<TH1*>(spd->FindObject("events"));

  // --- Set style parameters ----------------------------------------
  gStyle->SetPalette(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleX(0.69);
  gStyle->SetTitleY(0.99);
  gStyle->SetTitleW(0.30);
  gStyle->SetTitleH(0.10);

  // --- Make canvas for correlation plots and plot them -------------
  TCanvas* c1 = new TCanvas("c1", "c1", 900, 600);
  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.05);
  c1->SetFillColor(0);
  c1->SetFillStyle(0);
  c1->SetBorderSize(0);
  c1->SetBorderMode(0);
  
  c1->cd();
  TPad* p1 = new TPad("p1","p1", 0.0, 0.3, 1.0, 1.0, 0, 0, 0);
  p1->SetFillColor(0);
  p1->SetFillStyle(0);
  p1->SetBorderSize(0);
  p1->SetBorderMode(0);
  p1->SetTopMargin(0.05);
  p1->SetBottomMargin(0.001);
  p1->SetRightMargin(0.05);
  p1->Draw();
  p1->cd();
  clusters->SetMarkerStyle(20);
  clusters->SetMarkerColor(kRed+1);
  clusters->Draw();

  tracklets->SetMarkerStyle(20);
  tracklets->SetMarkerColor(kBlue+1);
  tracklets->Draw("same");

  TString v(Form("%+5.1f<|v_{z}|<%+5.1f",
		 events->GetXaxis()->GetXmin(),
		 events->GetXaxis()->GetXmax()));
  TLatex* ltx = new TLatex(.2,.80, v.Data());
  ltx->Draw();

  TLegend* l = p1->BuildLegend(.6,.6,.94,.94);
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);

  c1->cd();
  TPad* p2 = new TPad("p2","p2", 0.0, 0.0, 1.0, 0.3, 0, 0, 0);
  p2->SetFillColor(0);
  p2->SetFillStyle(0);
  p2->SetBorderSize(0);
  p2->SetBorderMode(0);
  p2->SetTopMargin(0.001);
  p2->SetBottomMargin(0.2);
  p2->SetRightMargin(0.05);
  p2->Draw();
  p2->cd();
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerColor(kGray+1);
  ratio->Draw();
  ratio->GetYaxis()->SetRangeUser(0,2);

  c1->SaveAs("SPDComparison.png");
  c1->SaveAs("SPDComparison.root");
}

    
  
  
  
#endif
//
// EOF
//

