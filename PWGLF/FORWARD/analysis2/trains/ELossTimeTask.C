#ifndef __CINT__
# include <AliFMDEventInspector.h>
# include <AliForwardCorrectionManager.h>
# include <TH2.h>
# include <AliESDFMD.h>
# include <AliAODForwardMult.h>
# include <AliESDEvent.h>
# include <TFile.h>
# include <TSystem.h>
#else
class AliFMDEventInspector;
class TH2;
class AliESDFMD;
#endif
#include <AliBaseESDTask.h>
#include <AliForwardUtil.h>
#define NO_TASK
#define NO_SORTER
#include "EventTimeTask.C"

/**
 * @defgroup pwglf_forward_eventtime Investigation of time-between-events 
 * @ingroup pwglf_forward_trains_specific
 */
/**
 * Task to analyse the energy loss in the FMD rings as function of the
 * time to previous event.
 * 
 * @ingroup pwglf_forward_eventtime
 */
struct ELossTimeTask : public AliBaseESDTask
{
  /** 
   * Constructor - for I/O only
   */
  ELossTimeTask() 
    : AliBaseESDTask(),
      fEventInspector(),
      fFMD1i(),
      fFMD2i(),
      fFMD2o(), 
      fFMD3i(),
      fFMD3o(),
      fMap(0),
      fDt(0)
  {}
  /** 
   * Constructor 
   * 
   * @param name Name of task 
   */     
  ELossTimeTask(const char* name) 
    : AliBaseESDTask(name, "ELossTimeTask", 
		     &(AliForwardCorrectionManager::Instance())),
      fEventInspector("fmdEventInspector"),
      fFMD1i(1, 'i'),
      fFMD2i(2, 'i'),
      fFMD2o(2, 'o'), 
      fFMD3i(3, 'i'),
      fFMD3o(3, 'o'),
      fMap(0),
      fDt(0)
  {
  }
  /** 
   * Book output objects. Derived class should define this to book
   * output objects on the processing output list @c fList before the
   * actual event processing.  This is called on the master and on
   * each slave.
   * 
   * If this member function returns false, the execution is stopped
   * with a fatal signal.
   *
   * @return true on success. 
   */
  virtual Bool_t Book() 
  {
    fNeededCorrections = 0;
    fExtraCorrections  = 0;

    fDt = new TH1D("dt", "Time-to-last event (PS-triggered)", 60, 0, 15);
    fDt->SetXTitle("log_{10}(#Deltat)");
    fDt->SetFillColor(kYellow+2);
    fDt->SetLineColor(kYellow+2);
    fDt->SetFillStyle(3001);
    fDt->SetDirectory(0);
    fList->Add(fDt);

    fFMD1i.Book(fList, fDt);
    fFMD2i.Book(fList, fDt);
    fFMD2o.Book(fList, fDt);
    fFMD3i.Book(fList, fDt);
    fFMD3o.Book(fList, fDt);

    // Possibly re-read map
    ReadMap("map.root");

    return true;
  }
  /** 
   * Process a single event
   * 
   * @param esd Input event 
   * 
   * @return true on success 
   */
  virtual Bool_t Event(AliESDEvent& esd)
  {
    Bool_t   lowFlux   = kFALSE;
    UInt_t   triggers  = 0;
    UShort_t ivz       = 0;
    TVector3 ip;
    Double_t cent      = 0;
    UShort_t nClusters = 0;
    UInt_t   found     = fEventInspector.Process(&esd, triggers, lowFlux, 
						 ivz, ip, cent, nClusters);
    if (found & AliFMDEventInspector::kNoEvent)    return false;
    if (found & AliFMDEventInspector::kNoTriggers) return false;
    if (found & AliFMDEventInspector::kNoSPD)      return false;
    if (found & AliFMDEventInspector::kNoFMD)      return false;
    if (found & AliFMDEventInspector::kNoVertex)   return false;
    if (found & AliFMDEventInspector::kBadVertex)  return false;

    // do not process pile-up, A, C, and E events 
    if (triggers & AliAODForwardMult::kPileUp)     return false;
    if (triggers & AliAODForwardMult::kA)          return false;
    if (triggers & AliAODForwardMult::kC)          return false;
    if (triggers & AliAODForwardMult::kE)          return false;
  
    // We want only the events found by off-line 
    if (!(triggers & AliAODForwardMult::kOffline)) return false;
    
    // Perhaps we should also insist on MB only 
    // if (fOnlyMB && (!(triggers & AliAODForwardMult::kInel))) return false;
	  
    // Get FMD data 
    AliESDFMD* esdFMD = esd.GetFMDData();
    ULong64_t  period = esd.GetPeriodNumber();
    ULong64_t  orbit  = esd.GetOrbitNumber();
    ULong64_t  bc     = esd.GetBunchCrossNumber();
    ULong64_t  full   = EventTimeData::EncodeFull(period, orbit, bc);
    ULong64_t  dt     = fMap->Get(full);
    Double_t   logDt  = TMath::Log10(25. * dt);
    if (dt == EventTimeMap::kInvalidTime) {
      logDt = 0;
      Printf("!!! Failed to find dT for 0x%016llu", full);
    }
    // else 
    //   Printf("=== 0x%016llu -> 0x%016llu", full, dt);
    fDt->Fill(logDt);
    
    for (UShort_t d = 1; d <= 3; d++) { 
      UShort_t nQ = d == 1 ? 1 : 2;
      for (UShort_t q = 0; q < nQ; q++) { 
	RingHistos* r = 0;
	switch (d) { 
	case 1: r = &fFMD1i; break;
	case 2: r = (q == 0 ? &fFMD2i : &fFMD2o); break;
	case 3: r = (q == 0 ? &fFMD3i : &fFMD3o); break;
	}
	r->Event(*esdFMD, logDt);
      }
    }
    return true;
  }
  /** 
   * Do the final analysis on the merged output. 
   * 
   * @return true on success
   */
  virtual Bool_t Finalize() 
  { 
    fDt = static_cast<TH1*>(fList->FindObject("dt"));

    fFMD1i.Finalize(fList, fResults, fDt);
    fFMD2i.Finalize(fList, fResults, fDt);
    fFMD2o.Finalize(fList, fResults, fDt);
    fFMD3i.Finalize(fList, fResults, fDt);
    fFMD3o.Finalize(fList, fResults, fDt);
    return true; 
  }
  /** 
   * Get a reference to the event inspector. User must override this
   * to return proper object
   * 
   * @return Reference to the event inspector 
   */
  virtual AliFMDEventInspector& GetEventInspector() { return fEventInspector; }
  /** 
   * Get a reference to the event inspector. User must override this
   * to return proper object
   * 
   * @return Reference to the event inspector 
   */
  virtual const AliFMDEventInspector& GetEventInspector() const 
  {
    return fEventInspector;
  }
  /** 
   * Read the map from timestamp to time-to-previous event 
   * 
   * @param filename File to read the map from 
   * 
   * @return true on success, false otherwise 
   */
  Bool_t ReadMap(const char* filename)
  {
    if (gSystem->AccessPathName(filename, kReadPermission)) {
      // TSystem::AccessPathName returns false if file can be accessed!
      Error("ReadMap", "File \"%s\" cannot be open for reading", filename);
      return false;
    }
    Printf("Opening \"%s\" ...", filename);
    TFile* file = TFile::Open(filename, "READ");
    if (!file) { 
      Error("ReadMap", "Failed to open file \"%s\"", filename);
      return false;
    }
    Printf("Opened \"%s\" ...", filename);
    TObject* o = file->Get("eventTimeMap");
    if (!o) { 
      Error("ReadMap", "Failed to get \"eventTimeMap\" from %s", filename);
      return false;
    }
    Printf("Got object \"eventTimeMap\" from \"%s\" ...", filename);
    if (!o->IsA()->InheritsFrom(EventTimeMap::Class())) { 
      Error("ReadMap", "Object \"%s\" is not an EventTimeMap, but a %s", 
	    o->GetName(), o->ClassName());
      return false;
    }
    Printf("Set the \"eventTimeMap\" to use");
    if (fMap) { 
      delete fMap;
      fMap = 0;
    }
    fMap = static_cast<EventTimeMap*>(o);
    file->Close();
    return true;
  }
  /** 
   * Create and connect the task 
   * 
   * @param mapfile File name of file containing timestamp map
   * 
   * @return true on connect
   */
  static Bool_t Create(const char* mapfile) 
  {
    ELossTimeTask* task = new ELossTimeTask("elossTime");
    if (!task->ReadMap(mapfile)) return false;
    task->Connect();
    return true;
  }
protected:
  /** 
   * Dummy copy constructor 
   * 
   * @param o Object to copy from 
   */
  ELossTimeTask(const ELossTimeTask& o) : AliBaseESDTask(o) {}
  /** Our event inspector */
  AliFMDEventInspector fEventInspector;
  /** 
   * Structure to hold per-ring histograms 
   */
  struct RingHistos : public AliForwardUtil::RingHistos
  {
    /** 
     * Constructor - for I/O only
     */
    RingHistos() 
      : AliForwardUtil::RingHistos(),
	fDtVsELoss(0)
    {}
    /** 
     * Constructor 
     * 
     * @param d detector number 
     * @param r ring identifier 
     */
    RingHistos(UShort_t d, Char_t r) 
      : AliForwardUtil::RingHistos(d,r),
	fDtVsELoss(0)
    {}
    /** 
     * Book histograms 
     * 
     * @param dir Parent list to add our list to
     * @param dt  Histogram of time differences 
     * 
     * @return true on success 
     */
    Bool_t Book(TList* dir, TH1* dt)
    {
      TList* l = DefineOutputList(dir);
      // Double_t dtBins[] = { 0, 1, 2e3, 5e5, 1.1e6, 5e6, 1e12, 1e20 };
      fDtVsELoss = new TH2D("dtVsELoss", 
			    Form("#Deltat vs #Delta/#Delta_{mip} - %s", 
				 GetName()), 450, 0, 15, 
			    dt->GetXaxis()->GetNbins(), 
			    dt->GetXaxis()->GetXmin(), 
			    dt->GetXaxis()->GetXmax());
      fDtVsELoss->SetXTitle("#Delta/#Delta_{mip}");
      fDtVsELoss->SetYTitle("log_{10}(#Deltat) [ns]");
      fDtVsELoss->SetMarkerColor(Color());
      fDtVsELoss->SetDirectory(0);
      l->Add(fDtVsELoss);

      return true;
    }
    /** 
     * Process a single event
     * 
     * @param fmd   FMD ESD data
     * @param logDt Logarithm (base 10) of Time-to-previous event
     */
    void Event(AliESDFMD& fmd, Double_t logDt)
    {
      const UShort_t nSec = NSector();
      const UShort_t nStr = NStrip();
      for (UShort_t s = 0; s < nSec; s++) { 
	for (UShort_t t = 0; t < nStr; t++) { 
	  Float_t mult = fmd.Multiplicity(fDet,fRing,s,t);
	  if(mult == AliESDFMD::kInvalidMult || mult <= 0) 
	    continue;
	  fDtVsELoss->Fill(mult, logDt);
	}
      }
    }
    void Finalize(const TList* input, TList* output, TH1* dt)
    {
      TList* in  = GetOutputList(input);
      TList* out = DefineOutputList(output);

      fDtVsELoss = static_cast<TH2*>(in->FindObject("dtVsELoss"));

      TH2* dtVsELoss = static_cast<TH2*>(fDtVsELoss->Clone());
      dtVsELoss->SetDirectory(0);
      dtVsELoss->SetZTitle("1/N_{events}");
      dtVsELoss->Reset();
      for (UShort_t j = 1; j <= dt->GetNbinsX(); j++) { 
	Double_t norm = dt->GetBinContent(j);
	if (norm <= 0) {
	  // Warning("Finalize", "Bin %d empty in dT", j);
	  continue;
	}
	for (UShort_t i = 1; i <= fDtVsELoss->GetNbinsX(); i++) { 
	  Double_t c = fDtVsELoss->GetBinContent(i, j);
	  Double_t e = fDtVsELoss->GetBinError(i, j);
	  dtVsELoss->SetBinContent(i,j,c/norm);
	  dtVsELoss->SetBinError(i,j,e/norm);
	}
      }
      out->Add(dtVsELoss);
    }
    /** Our histogram */
    TH2* fDtVsELoss;
    ClassDef(RingHistos,1)
  };
  /** Container of FMD1i histograms */
  RingHistos    fFMD1i;
  /** Container of FMD2i histograms */
  RingHistos    fFMD2i;
  /** Container of FMD2o histograms */
  RingHistos    fFMD2o;
  /** Container of FMD3i histograms */
  RingHistos    fFMD3i;
  /** Container of FMD3o histograms */
  RingHistos    fFMD3o;
  /** Map from timestamp to time-to-previous event*/
  EventTimeMap* fMap;
  /** Distribution of log10(dt) */
  TH1* fDt;
  ClassDef(ELossTimeTask,2)
};
//
// EOF
//
