#ifdef BUILD
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliMultiplicity.h"
#include "AliFMDMCEventInspector.h"
#include "AliAODForwardMult.h"
#include "AliLog.h"
#include <TH1I.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TList.h>
#include <TObjArray.h>
#include <TParameter.h>
#include <TStopwatch.h>
#include <TROOT.h>
#include <THStack.h>
#include <TStyle.h>

//====================================================================
/**
 * Task to evaluate trigger bias in pp 
 * 
 */
class EvaluateTrigger : public AliAnalysisTaskSE
{
public:
  enum { 
    kNone = 0x0, 
    kESD  = 0x1, 
    kMC   = 0x2
  };
  /** 
   * Constructor
   */
  EvaluateTrigger() 
    : AliAnalysisTaskSE(),
      fInel(),
      fInelGt0(),
      fNSD(),
      fInspector(), 
      fFirstEvent(true),
      fData(0), 
      fTriggers(0), 
      fTrackletRequirement(kESD),
      fVertexRequirement(kESD), 
      fVertexAxis(0, 0, 0), 
      fVertexESD(0),
      fVertexMC(0), 
      fM(0)
  {}
  /** 
   * Constructor 
   */
  EvaluateTrigger(const char* /*name*/) 
    : AliAnalysisTaskSE("evaluateTrigger"),
      fInel(AliAODForwardMult::kInel),
      fInelGt0(AliAODForwardMult::kInelGt0),
      fNSD(AliAODForwardMult::kNSD),
      fInspector("eventInspector"), 
      fFirstEvent(true), 
      fData(0), 
      fTriggers(0),
      fTrackletRequirement(kESD),
      fVertexRequirement(kESD), 
      fVertexAxis(10, -10, 10), 
      fVertexESD(0),
      fVertexMC(0), 
      fM(0)
  {
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
  }
  void SetVertexRequirement(UShort_t m) { fVertexRequirement = m; }
  void SetTrackletRequirement(UShort_t m) { fTrackletRequirement = m; }
  void SetVertexAxis(Int_t nBins, Double_t low, Double_t high) 
  {
    fVertexAxis.Set(nBins, low, high);
  }
  /** 
   * Intialize 
   */
  void Init() {}
  /** 
   * Create output objects 
   */
  void UserCreateOutputObjects()
  {
    fList = new TList;
    fList->SetOwner();
    fList->SetName(GetName());

    Double_t mb[] = { 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 1000 };
    Int_t    nM   = 10;
    TAxis mAxis(nM, mb);
    TAxis eAxis(200, -4, 6);
    TAxis pAxis(40, 0, 2*TMath::Pi());

    fData = new TH2D("data", "Cache", 
		     eAxis.GetNbins(), eAxis.GetXmin(), eAxis.GetXmax(), 
		     pAxis.GetNbins(), pAxis.GetXmin(), pAxis.GetXmax());
    fData->SetDirectory(0);
    fData->SetXTitle("#eta");
    fData->SetYTitle("#varphi [radians]");
    fData->SetZTitle("N_{ch}(#eta,#varphi)");
    fData->Sumw2();
    
    fM = new TH1D("m", "Distribution of N_{ch}|_{|#eta|<1}",nM+1, -0.5, nM+.5); 
    fM->SetXTitle("N_{ch}|_{|#eta|<1}");
    fM->SetYTitle("Events");
    fM->SetFillColor(kRed+1);
    fM->SetFillStyle(3001);
    fM->SetDirectory(0);
    fList->Add(fM);

    for (Int_t i = 0; i <= nM; i++) { 
      TString lbl;
      if (i == 0)       lbl = "all";
      else if (i == nM) lbl = Form("%d+",int(mAxis.GetBinLowEdge(i)+.5));
      else              lbl = Form("<%d",int(mAxis.GetBinUpEdge(i)+.5));
      fM->GetXaxis()->SetBinLabel(i+1, lbl);
    }

    fTriggers = new TH1I("triggers", "Triggers", 6, -.5, 5.5);
    fTriggers->SetDirectory(0);
    fTriggers->GetXaxis()->SetBinLabel(1, "INEL (MC)");
    fTriggers->GetXaxis()->SetBinLabel(2, "INEL (ESD)");
    fTriggers->GetXaxis()->SetBinLabel(3, "INEL>0 (MC)");
    fTriggers->GetXaxis()->SetBinLabel(4, "INEL>0 (ESD)");
    fTriggers->GetXaxis()->SetBinLabel(5, "NSD (MC)");
    fTriggers->GetXaxis()->SetBinLabel(6, "NSD (ESD)");
    fTriggers->SetFillColor(kYellow+1);
    fTriggers->SetFillStyle(3001);
    fList->Add(fTriggers);

    fVertexESD = new TH1D("vertexESD", "ESD vertex distribution", 
			  fVertexAxis.GetNbins(), 
			  fVertexAxis.GetXmin(), 
			  fVertexAxis.GetXmax());
    fVertexESD->SetDirectory(0);
    fVertexESD->SetFillColor(kRed+1);
    fVertexESD->SetFillStyle(3001);
    fVertexESD->SetXTitle("v_{z} [cm]");
    fVertexESD->SetYTitle("P(v_{z})");
    fList->Add(fVertexESD);

    fVertexMC = new TH1D("vertexMC", "MC vertex distribution", 
			  fVertexAxis.GetNbins(), 
			  fVertexAxis.GetXmin(), 
			  fVertexAxis.GetXmax());
    fVertexMC->SetDirectory(0);
    fVertexMC->SetFillColor(kBlue+1);
    fVertexMC->SetFillStyle(3001);
    fVertexMC->SetXTitle("v_{z} [cm]");
    fVertexMC->SetYTitle("P(v_{z})");
    fList->Add(fVertexMC);

    fInel.CreateObjects(fList, fM, fData);
    fInelGt0.CreateObjects(fList, fM, fData);
    fNSD.CreateObjects(fList, fM, fData);


    fInspector.DefineOutput(fList);
    fInspector.Init(fVertexAxis);

    PostData(1, fList);
  }
  /** 
   * Event processing 
   */
  void UserExec(Option_t*) 
  {
    // Get the input data - MC event
    AliMCEvent*  mcEvent = MCEvent();
    if (!mcEvent) { 
      AliWarning("No MC event found");
      return;
    }
    
    // Get the input data - ESD event
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!esd) { 
      AliWarning("No ESD event found for input event");
      return;
    }

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
      
      fFirstEvent = false;
    }

    // Get the particle stack 
    AliStack* stack = mcEvent->Stack();

    // Some variables 
    UInt_t   triggers; // Trigger bits
    Bool_t   lowFlux;  // Low flux flag
    UShort_t iVz;      // Vertex bin from ESD
    Double_t vZ;       // Z coordinate from ESD
    Double_t cent;     // Centrality 
    UShort_t iVzMc;    // Vertex bin from MC
    Double_t vZMc;     // Z coordinate of IP vertex from MC
    Double_t b;        // Impact parameter
    Int_t    nPart;    // Number of participants 
    Int_t    nBin;     // Number of binary collisions 
    Double_t phiR;     // Reaction plane from MC
    
    // Process the data 
    Int_t retESD = fInspector.Process(esd, triggers, lowFlux, iVz, vZ, cent);
    Int_t retMC  = fInspector.ProcessMC(mcEvent, triggers, iVzMc, 
					vZMc, b, nPart, nBin, phiR);

    Bool_t hasESDVtx = retESD == AliFMDEventInspector::kOk;
    Bool_t hasMCVtx  = retMC  == AliFMDEventInspector::kOk;
    if (hasESDVtx) fVertexESD->Fill(vZ);
    if (hasMCVtx)  fVertexMC->Fill(vZMc);

    Bool_t isMcInel = true; // (triggers & AliAODForwardMult::kB);
    Bool_t isMcNSD  = (triggers & AliAODForwardMult::kMCNSD);

    Int_t mESD = 0;
    const AliMultiplicity* spdmult = esd->GetMultiplicity();
    if (!spdmult) {
      AliWarning("No SPD multiplicity");
    }
    else { 
      // Check if we have one or more tracklets 
      // in the range -1 < eta < 1 to set the INEL>0 
      // trigger flag. 
      Int_t n = spdmult->GetNumberOfTracklets();
      for (Int_t j = 0; j < n; j++) 
	if(TMath::Abs(spdmult->GetEta(j)) < 1) mESD++;
    }

    // Reset cache 
    fData->Reset();
    Int_t mMC = 0; // Number of particles in |eta|<1

    // Loop over all tracks 
    Int_t nTracks = mcEvent->GetNumberOfTracks();
    for (Int_t iTr = 0; iTr < nTracks; iTr++) { 
      AliMCParticle* particle = 
	static_cast<AliMCParticle*>(mcEvent->GetTrack(iTr));
    
      // Check the returned particle 
      if (!particle) continue;

      // Check if this charged and a primary 
      Bool_t isCharged = particle->Charge() != 0;
      Bool_t isPrimary = stack->IsPhysicalPrimary(iTr);

      if (!isCharged || !isPrimary) continue;

      
      // Fill (eta,phi) of the particle into histograsm for b
      Double_t eta = particle->Eta();
      Double_t phi = particle->Phi();
      
      fData->Fill(eta, phi);
      if (TMath::Abs(eta) <= 1) mMC++;
    }
    Int_t m = mESD;
    if (fTrackletRequirement == kMC) m = mMC;
    fM->Fill(m);

    bool isMcInelGt0 = isMcInel && (mMC > 0);
    
    bool hasVertex   = true;
    if (fVertexRequirement & kMC)  hasVertex = hasVertex && hasMCVtx;
    if (fVertexRequirement & kESD) hasVertex = hasVertex && hasESDVtx;

    if (isMcInel) {
      fTriggers->Fill(0);
      bool triggered = (triggers & AliAODForwardMult::kInel);
      if (triggered) fTriggers->Fill(1);
      fInel.AddEvent(triggered, hasVertex, m, fData);
    }
    if (isMcInelGt0) {
      fTriggers->Fill(2);
      bool triggered = (triggers & AliAODForwardMult::kInelGt0);
      if (triggered) fTriggers->Fill(3);
      fInelGt0.AddEvent(triggered, hasVertex, m, fData);
    }
    if (isMcNSD) {
      fTriggers->Fill(4);
      bool triggered = (triggers & AliAODForwardMult::kNSD);
      if (triggered) fTriggers->Fill(5);
      fNSD.AddEvent(triggered, hasVertex, m, fData);
    }
    PostData(1, fList);
  }
  /** 
   * End of job processing 
   */
  void Terminate(Option_t*)
  {
    fList = dynamic_cast<TList*>(GetOutputData(1));
    if (!fList) {
      AliError(Form("No output list defined (%p)", GetOutputData(1)));
      if (GetOutputData(1)) GetOutputData(1)->Print();
      return;
    }


    TList* output = new TList;
    output->SetName(GetName());
    output->SetOwner();

    fVertexMC = static_cast<TH1D*>(fList->FindObject("vertexMC"));
    fVertexESD = static_cast<TH1D*>(fList->FindObject("vertexESD"));
    fM         = static_cast<TH1D*>(fList->FindObject("m"));
    if (fVertexMC) { 
      TH1D* vtxMC = static_cast<TH1D*>(fVertexMC->Clone("vertexMC"));
      vtxMC->SetDirectory(0);
      vtxMC->Scale(1. / vtxMC->GetEntries());
      output->Add(vtxMC);
    }
    if (fVertexESD) { 
      TH1D* vtxESD = static_cast<TH1D*>(fVertexESD->Clone("vertexESD"));
      vtxESD->SetDirectory(0);
      vtxESD->Scale(1. / vtxESD->GetEntries());
      output->Add(vtxESD);
    }
    if (fM) { 
      TH1D* m = static_cast<TH1D*>(fM->Clone("m"));
      m->SetDirectory(0);
      m->SetYTitle("P(N_{ch}|_{|#eta|<1} < X)");
      m->Scale(1. / m->GetBinContent(1));
      output->Add(m);
    }      

    fInel.Finish(fList, output);
    fInelGt0.Finish(fList, output);
    fNSD.Finish(fList, output);

    PostData(2, output);
  }
    
protected:
  //__________________________________________________________________
  /** 
   * Structure to hold per trigger type information 
   */
  struct TriggerType : public TNamed
  {
    //________________________________________________________________
    /** 
     * Structure to hold per multiplicity bin information 
     */
    struct MBin : public TNamed
    {
      TH2D*     fTruth;
      TH2D*     fTriggered; 
      TH2D*     fAccepted;
      TH1I*     fCounts;
      UShort_t  fLow;
      UShort_t  fHigh;
      /** 
       * Constructor 
       */
      MBin() 
	: fTruth(0), fTriggered(0), fAccepted(0), 
	  fCounts(0), fLow(0), fHigh(1000) {}
      /** 
       * Constructor 
       * 
       * @param p      Parent list 
       * @param low    Low cut 
       * @param high   High cut 
       * @param eAxis  Eta axis 
       * @param pAxis  Phi axis 
       */
      MBin(TList* p, UShort_t low, UShort_t high, const TH2D* dHist) 
	: fTruth(0), 
	  fTriggered(0), 
	  fAccepted(0),
	  fCounts(0), 
	  fLow(low), 
	  fHigh(high)
      {
	if (low >= high) SetName("all");
	else             SetName(Form("m%03d_%03d", fLow, fHigh));
	TList* l = new TList;
	l->SetOwner();
	l->SetName(GetName());
	p->Add(l);

	fTruth = static_cast<TH2D*>(dHist->Clone(("truth")));
	fTruth->SetTitle("MC truth");
	fTruth->SetDirectory(0);
	fTruth->SetZTitle("#sum_i^{N_X} N_{ch}(#eta,#varphi)");
	fTruth->Sumw2();
	fTruth->Reset();
	l->Add(fTruth);

	fTriggered = static_cast<TH2D*>(fTruth->Clone(("triggered")));
	fTriggered->SetTitle("Triggered");
	fTriggered->SetDirectory(0);
	fTriggered->SetZTitle("#sum_i^{N_T} N_{ch}(#eta,#varphi)");
	fTriggered->Sumw2();
	fTriggered->Reset();
	l->Add(fTriggered);

	fAccepted = static_cast<TH2D*>(fTruth->Clone(("accepted")));
	fAccepted->SetTitle("Accepted");
	fAccepted->SetDirectory(0);
	fAccepted->SetZTitle("#sum_i^{N_T} N_{ch}(#eta,#varphi)");
	fAccepted->Sumw2();
	fAccepted->Reset();
	l->Add(fAccepted);
	
	fCounts = new TH1I("counts", "Event counts", 3, -.5, 2.5);
	fCounts->SetDirectory(0);
	fCounts->GetXaxis()->SetBinLabel(1, "Truth");
	fCounts->GetXaxis()->SetBinLabel(2, "Triggered");
	fCounts->GetXaxis()->SetBinLabel(3, "Accepted");
	fCounts->SetYTitle("# events");
	l->Add(fCounts);
      }
      /** 
       * Add event observation
       * 
       * @param triggered Whether the event was triggered
       * @param event     Data for this event 
       */
      void AddEvent(Bool_t triggered, Bool_t hasVtx, const TH2D* event) 
      {
	fCounts->Fill(0);
	fTruth->Add(event);
	if (triggered) { 
	  fCounts->Fill(1);
	  fTriggered->Add(event);
	  if (hasVtx) {
	    fCounts->Fill(2);
	    fAccepted->Add(event);
	  }
	}
      }
      /** 
       * End of job processing 
       * 
       * @param p      Parent list
       * @param o      Output parent list
       * @param stack  Stack of histograms
       * 
       * @return Trigger efficiency
       */ 
     Double_t Finish(const TList* p, TList* o, THStack* stack) 
      {
	TList* l = dynamic_cast<TList*>(p->FindObject(GetName()));
	if (!l) { 
	  Warning("Finish", "Cannot find %s in %s", GetName(), p->GetName());
	  return 0;
	}
	fTruth     = static_cast<TH2D*>(l->FindObject("truth"));
	fTriggered = static_cast<TH2D*>(l->FindObject("triggered"));
	fAccepted  = static_cast<TH2D*>(l->FindObject("accepted"));
	fCounts    = static_cast<TH1I*>(l->FindObject("counts"));
	
	Int_t    nTruth     = fCounts->GetBinContent(1);
	Int_t    nTriggered = fCounts->GetBinContent(2);
	Int_t    nAccepted  = fCounts->GetBinContent(3);
	Double_t eff        = 0;
	if (nTruth > 0) eff = double(nTriggered) / nTruth;
	else if (nTriggered == nTruth) eff = 1;

	if (nTruth > 0)     fTruth->Scale(1. / nTruth);
	if (nTriggered > 0) fTriggered->Scale(1. / nTriggered);
	if (nAccepted > 0)  fAccepted->Scale(1. / nAccepted);

	if (fLow >= fHigh) 
	  Info("Finish", "%-6s  [all]  E_X=N_T/N_X=%9d/%-9d=%f "
	       "E_V=N_A/N_T=%9d/%-9d=%f", 
	       p->GetName(), nTriggered, nTruth, eff, nAccepted, nTriggered, 
	       (nTriggered > 0 ? double(nAccepted) / nTriggered: 0));
	else
	  Info("Finish", "%-6s [%2d-%2d] E_X=N_T/N_X=%9d/%-9d=%f "
	       "E_V=N_A/N_T=%9d/%-9d=%f", 
	       p->GetName(), fLow, fHigh, nTriggered, nTruth, eff, 
	       nAccepted, nTriggered, 
	       (nTriggered > 0 ? double(nAccepted) / nTriggered: 0));
	
	TList* out = new TList;
	out->SetName(GetName());
	out->SetOwner();
	o->Add(out);
	
	out->Add(fTruth);
	out->Add(fTriggered);
	out->Add(new TParameter<double>("eff", eff));
	
	TH2D* bias = static_cast<TH2D*>(fAccepted->Clone("bias"));
	bias->Divide(fTruth);
	bias->SetDirectory(0);
	bias->SetZTitle("Trigger bias (accepted/truth)");
	out->Add(bias);

	TH1D* truth_px = static_cast<TH1D*>(fTruth->ProjectionX("truth_eta"));
	truth_px->SetTitle(Form("%d #leq N_{tracklets}|_{|#eta|<1} < %d", 
			       fLow, fHigh));
	truth_px->Scale(1. / fTruth->GetNbinsY());
	truth_px->SetDirectory(0);
	truth_px->SetLineColor(kRed+1);
	truth_px->SetMarkerColor(kRed+1);
	truth_px->SetFillColor(kRed+1);
	truth_px->SetFillStyle(0);
	out->Add(truth_px);

	TH1D* triggered_px = 
	  static_cast<TH1D*>(fTriggered->ProjectionX("triggered_eta"));
	triggered_px->SetTitle(Form("%d #leq N_{tracklets}|_{|#eta|<1} < %d", 
			       fLow, fHigh));
	triggered_px->Scale(1. / fTriggered->GetNbinsY());
	triggered_px->SetDirectory(0);
 	triggered_px->SetLineColor(kGreen+1);
	triggered_px->SetMarkerColor(kGreen+1);
	triggered_px->SetFillColor(kGreen+1);
	triggered_px->SetFillStyle(0);
	out->Add(triggered_px);

	TH1D* accepted_px = 
	  static_cast<TH1D*>(fAccepted->ProjectionX("accepted_eta"));
	accepted_px->SetTitle(Form("%d #leq N_{tracklets}|_{|#eta|<1} < %d", 
				   fLow, fHigh));
	accepted_px->Scale(1. / fAccepted->GetNbinsY());
 	accepted_px->SetLineColor(kBlue+1);
	accepted_px->SetMarkerColor(kBlue+1);
	accepted_px->SetFillColor(kBlue+1);
	accepted_px->SetDirectory(0);
	out->Add(accepted_px);

	THStack* data = new THStack("data", "Data distributions");
	data->Add(truth_px);
	data->Add(triggered_px);
	data->Add(accepted_px);
	out->Add(data);

#if 0
	TH1D* bias_px = static_cast<TH1D*>(bias->ProjectionX("bias_eta"));
	bias_px->SetTitle(Form("%d #leq N_{tracklets}|_{|#eta|<1} < %d", 
			       fLow, fHigh));
	bias_px->Scale(1. / bias->GetNbinsY());
#else
	TH1D* bias_px = static_cast<TH1D*>(accepted_px->Clone("bias_px"));
	bias_px->Divide(truth_px);
	bias_px->SetYTitle("Trigger bias (triggered/truth)");
#endif
	bias_px->SetDirectory(0);
	bias_px->SetMarkerStyle(20);
	bias_px->SetFillStyle(0);
	bias_px->SetMinimum(0);
	out->Add(bias_px);

	stack->Add(bias_px);

	return eff;
      }	
      ClassDef(MBin,1);
    };

    /** 
     * Constructor 
     */
    TriggerType() 
      : TNamed(), 
	fMask(0),
	fM(0),
	fBins(0)
    {}
    //--- Constructor ------------------------------------------------
    /** 
     * Constructor 
     * 
     * @param mask  Trigger mask
     */
    TriggerType(UShort_t mask) 
      : TNamed(AliAODForwardMult::GetTriggerString(mask), ""),
	fMask(mask), 
	fM(0), 
	fBins(0)
    {
    }
    /** 
     * Create objects 
     * 
     * @param list   PArent list
     * @param mAxis  Multiplicity axis 
     * @param eAxis  Eta axis
     * @param pAxis  Phi axis
     */
    void CreateObjects(TList* list, 
		       const TH1D* mHist, 
		       const TH2D* dHist)
    {
      TList* ours  = new TList;
      ours->SetName(GetName());
      ours->SetOwner();
      list->Add(ours);

      fM = static_cast<TH1D*>(mHist->Clone("m"));
      fM->SetDirectory(0);
      ours->Add(fM);
      
      fBins = new TObjArray;
      fBins->AddAt(new MBin(ours, 0, 0, dHist), 0);
      for (Int_t i = 1; i <= fM->GetNbinsX(); i++) { 
	Double_t low  = fM->GetXaxis()->GetBinLowEdge(i);
	Double_t high = fM->GetXaxis()->GetBinUpEdge(i);

	fBins->AddAt(new MBin(ours, low, high, dHist), i);
      }
    }
    /** 
     * Find bin corresponding to m
     * 
     * @param m Multiplicity 
     * 
     * @return Bin. 
     */
    MBin* FindBin(UShort_t m) 
    { 
      Int_t b = fM->GetXaxis()->FindBin(m);
      if (b <= 0) return 0;
      if (b >= fM->GetNbinsX()+1) b = fM->GetNbinsX();
      return static_cast<MBin*>(fBins->At(b));
    }
    /** 
     * Add event observation 
     * 
     * @param triggered  IF this is triggered
     * @param m          Multiplicity 
     * @param data       Observation 
     */
    void AddEvent(Bool_t triggered, Bool_t hasVtx, UShort_t m, const TH2D* data)
    {
      fM->AddBinContent(1);
      fM->AddBinContent(TMath::Min(fM->GetNbinsX(), m+2));

      MBin* all = static_cast<MBin*>(fBins->At(0));
      all->AddEvent(triggered, hasVtx, data);
      
      MBin* bin = FindBin(m);
      bin->AddEvent(triggered, hasVtx, data);      
    }      
    /** 
     * End of job processing 
     * 
     * @param p Parent list 
     * @param o Parent output list 
     */
    void Finish(const TList* p, TList* o)
    {
      TList* l = dynamic_cast<TList*>(p->FindObject(GetName()));
      if (!l) { 
	Warning("Finish", "Cannot find %s in %s", GetName(), p->GetName());
	return;
      }
      
      TList* ours  = new TList;
      ours->SetName(GetName());
      ours->SetOwner();
      o->Add(ours);

      fM = static_cast<TH1D*>(l->FindObject("m"));
      if (!fM) { 
	Warning("Finish", "Didn't find object 'm' in %s", l->GetName());
	return;
      }
      TH1D* m = static_cast<TH1D*>(fM->Clone("m"));
      m->SetDirectory(0);
      m->SetYTitle("P(N_{ch}|_{|#eta|<1} < X)");
      if (m->GetBinContent(1) > 0) 
	m->Scale(1. / m->GetBinContent(1));
      ours->Add(m);

      Int_t nBin = fM->GetNbinsX();
      TH1D* effs = static_cast<TH1D*>(fM->Clone("effs"));
      effs->SetYTitle("#epsilon_{X}");
      effs->SetFillColor(kRed+1);
      effs->SetDirectory(0);
      effs->SetMinimum(0);

      gStyle->SetPalette(1);
      Int_t   nCol = gStyle->GetNumberOfColors();
      THStack* stack = new THStack("biases", "Trigger biases");
      for (Int_t i = 0; i <= nBin; i++) { 
	MBin* bin = static_cast<MBin*>(fBins->At(i));
	effs->SetBinContent(i+1, bin->Finish(l, ours, stack));
	TH1* h = static_cast<TH1*>(stack->GetHists()->At(i));
	Int_t col = kBlack;
	if (i != 0) { 
	  Int_t icol = TMath::Min(nCol-1,int(double(i)/nBin * nCol + .5));
	  col        = gStyle->GetColorPalette(icol);
	}
	h->SetMarkerColor(col);
	h->SetFillColor(col);
	h->SetLineColor(col);
      }

      ours->Add(stack);
      ours->Add(effs);
    } 
    UShort_t   fMask;
    TH1D*      fM;
    TObjArray* fBins;
    ClassDef(TriggerType,1);
  };
  TriggerType            fInel;
  TriggerType            fInelGt0;
  TriggerType            fNSD;
  AliFMDMCEventInspector fInspector;
  TList*                 fList;
  Bool_t                 fFirstEvent;
  TH2D*                  fData;
  TH1I*                  fTriggers;
  UInt_t                 fTrackletRequirement;
  UInt_t                 fVertexRequirement;
  TAxis                  fVertexAxis;
  TH1D*                  fVertexESD;
  TH1D*                  fVertexMC;
  TH1D*                  fM;
  ClassDef(EvaluateTrigger,1);
};
#else 
//====================================================================
void MakeEvaluateTriggers(const char* esddir, 
			  Int_t       nEvents    = -1, 
			  UInt_t      vtx        = 0x1,
			  UInt_t      trk        = 0x1,
			  UInt_t      vz         = 10,
			  Int_t       proof      = 0)
{
  // --- Libraries to load -------------------------------------------
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  // --- Check for proof mode, and possibly upload pars --------------
  if (proof> 0) { 
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadPars.C");
    if (!LoadPars(proof)) { 
      Error("MakeAOD", "Failed to load PARs");
      return;
    }
  }
  
  // --- Our data chain ----------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/MakeChain.C");
  TChain* chain = MakeChain("ESD", esddir,true);
  // If 0 or less events is select, choose all 
  if (nEvents <= 0) nEvents = chain->GetEntries();
  
  // --- Set the macro path ------------------------------------------
  gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2:"
			   "$ALICE_ROOT/ANALYSIS/macros",
			   gROOT->GetMacroPath()));

  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Triggers", 
						    "Forward multiplicity");
  AliAnalysisManager::SetCommonFileName("triggers.root");

  // --- ESD input handler -------------------------------------------
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdHandler);      
       
  // --- Monte Carlo handler -----------------------------------------
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  mcHandler->SetReadTR(true);    

  // --- Add tasks ---------------------------------------------------
  // Physics selection 
  gROOT->Macro(Form("AddTaskPhysicsSelection.C(1,1,0)"));

#if 0
  // --- Fix up physics selection to give proper A,C, and E triggers -
  AliInputEventHandler* ih =
    static_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
  AliPhysicsSelection* ps = 
    static_cast<AliPhysicsSelection*>(ih->GetEventSelection());
  // Ignore trigger class when selecting events.  This mean that we
  // get offline+(A,C,E) events too
  ps->SetSkipTriggerClassSelection(true);
#endif

  // --- compile our code --------------------------------------------
  gSystem->AddIncludePath("-I${ALICE_ROOT}/PWG2/FORWARD/analysis2 "
			  "-I${ALICE_ROOT}/ANALYSIS "
			  "-I${ALICE_ROOT}/include -DBUILD=1");
  gROOT->LoadMacro("${ALICE_ROOT}/PWG2/FORWARD/analysis2/MakeEvaluateTriggers.C++g");
  
  // --- Make our object ---------------------------------------------
  EvaluateTrigger* task = new EvaluateTrigger("triggers");
  mgr->AddTask(task);
  task->SetVertexRequirement(vtx);
  task->SetTrackletRequirement(trk);
  task->SetVertexAxis(10, -vz, vz);

  // --- create containers for input/output --------------------------
  AliAnalysisDataContainer *sums = 
    mgr->CreateContainer("triggerSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("triggerResults", TList::Class(), 
			 AliAnalysisManager::kParamContainer, 
			 AliAnalysisManager::GetCommonFileName());
  
  // --- connect input/output ----------------------------------------
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, sums);
  mgr->ConnectOutput(task, 2, output);
  
  // --- Run the analysis --------------------------------------------
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("MakeAOD", "Failed to initialize analysis train!");
    return;
  }
  // Skip terminate if we're so requested and not in Proof or full mode
  mgr->SetSkipTerminate(false);
  // Some informative output 
  mgr->PrintStatus();
  if (proof) mgr->SetDebugLevel(3);
  if (mgr->GetDebugLevel() < 1 && !proof) 
    mgr->SetUseProgressBar(kTRUE,100);

  // Run the train 
  t.Start();
  Printf("=== RUNNING ANALYSIS on %9d events ==========================",
	 nEvents);
  mgr->StartAnalysis(proof ? "proof" : "local", chain, nEvents);
  t.Stop();
  t.Print();
}
#endif
