//====================================================================
#include "AliBasedNdetaTask.h"
#include <TMath.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TF1.h>
#include <THStack.h>
#include <TList.h>
#include <TProfile.h>
#include <AliAnalysisManager.h>
#include <AliAODEvent.h>
#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include "AliForwardUtil.h"
#include "AliAODForwardMult.h"
#include "AliAODMultEventClass.h"
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TParameter.h>
#include <TColor.h>
#include <TLatex.h>

//====================================================================
namespace {
  /** 
   * Get bin number corresponding to the passed centrality 
   * 
   * @param c1 Least centrality 
   * @param c2 Largest centrality 
   * 
   * @return Bin number 
   */
  Int_t PbPbBin(Double_t c1, Double_t c2)
  {
    Double_t c = (c1+c2) / 2;
    if      (c <  5) return 0;
    else if (c < 10) return 1;
    else if (c < 20) return 2;
    else if (c < 30) return 3;
    else if (c < 40) return 4;
    else if (c < 50) return 5;
    else if (c < 60) return 6;
    else if (c < 70) return 7;
    else if (c < 80) return 8;
    else if (c < 90) return 9;
    return                  10;
  }
  //__________________________________________________________________
  /** 
   * Get the centrality color for PbPb 
   * 
   * @param c1 Lower edge
   * @param c2 Upper edge
   * 
   * @return Color 
   */
  Color_t PbPbColor(Double_t c1, Double_t c2)
  {
    const Color_t cc[] = { kMagenta+2,           // 0
			   kBlue+2,              // 1
			   kAzure-1,             // 2
			   kCyan+2,              // 3
			   kGreen+1,             // 4
			   kSpring+5,            // 5
			   kYellow+1,            // 6
			   kOrange+5,            // 7
			   kRed+1,               // 8
			   kPink+5,              // 9
			   kBlack };             // 10
    Int_t bin = PbPbBin(c1,c2);
    return cc[bin];
  }
  //__________________________________________________________________
  /** 
   * Get bin number corresponding to the passed centrality 
   * 
   * @param c1 Least centrality 
   * @param c2 Largest centrality 
   * 
   * @return Bin number 
   */
  Int_t pPbBin(Double_t c1, Double_t c2)
  {
    Double_t c   = (c1+c2) / 2;
    if      (c <   5) return 0;
    else if (c <  10) return 1;
    else if (c <  20) return 2;
    else if (c <  40) return 3;
    else if (c <  60) return 4;
    else if (c <  80) return 5;
    else if (c < 100) return 6;
    return                   7;
  }
  //__________________________________________________________________
  /** 
   * Get the centrality color pPb 
   * 
   * @param c1 Lower edge
   * @param c2 Upper edge
   * 
   * @return Color 
   */
  Color_t pPbColor(Double_t c1, Double_t c2)
  {
    const Color_t cc[] = { kMagenta+2, // 0
                           kBlue+2,    // 1
                           kCyan+2,    // 2
                           kGreen+1,   // 3
                           kYellow+1,  // 4
                           kRed+1,     // 5 
                           kPink+5,    // 6 
                           kBlack };   // 7
    Int_t bin = pPbBin(c1,c2);
    return cc[bin];
  }
}

//____________________________________________________________________
AliBasedNdetaTask::AliBasedNdetaTask()
  : AliBaseAODTask(), 
    fCorrEmpty(true), 
    fUseROOTProj(false),
    fTriggerEff(1),
    fTriggerEff0(1),
    fListOfCentralities(0),
    fNormalizationScheme(kFull), 
    fFinalMCCorrFile(""),
    fSatelliteVertices(0),
    fEmpiricalCorrection(0),
    fMeanVsC(0),
    fSeenCent(0),
    fTakenCent(0),
    fCentMethod("default"),
    fAnaUtil(),
    fUseUtilPileup(false),
    fIpzReweight(0),
    fCacheCent(-10),
    fCacheQual(0xFFFF)
{
  // 
  // Constructor
  // 
  DGUARD(fDebug,3,"Default CTOR of AliBasedNdetaTask");
}

//____________________________________________________________________
AliBasedNdetaTask::AliBasedNdetaTask(const char* name)
  : AliBaseAODTask(Form("%sdNdeta", name),"AliBasedNdetaTask"), 
    fCorrEmpty(true), 
    fUseROOTProj(false),
    fTriggerEff(1),
    fTriggerEff0(1),
    fListOfCentralities(0),
    fNormalizationScheme(kFull), 
    fFinalMCCorrFile(""),
    fSatelliteVertices(0),
    fEmpiricalCorrection(0),
    fMeanVsC(0),	
    fSeenCent(0),
    fTakenCent(0),
    fCentMethod("default"),
    fAnaUtil(),
    fUseUtilPileup(false),
    fIpzReweight(0),
    fCacheCent(-10),
    fCacheQual(0xFFFF)
{
  // 
  // Constructor
  // 
  DGUARD(fDebug, 3,"Named CTOR of AliBasedNdetaTask: %s", name);

  SetIPzAxis(10,-10,+10);
  fTriggerMask        = AliAODForwardMult::kInel;
  fListOfCentralities = new TObjArray(1);
  fListOfCentralities->SetName("centralityBins");
  // fListOfCentralities->SetTitle("Centrality bins");
  
  // Set the normalisation scheme 
  SetNormalizationScheme(kFull);
}


//____________________________________________________________________
AliBasedNdetaTask::~AliBasedNdetaTask()
{
  // 
  // Destructor
  // 
  DGUARD(fDebug,3,"Destruction of AliBasedNdetaTask");
}

//________________________________________________________________________
void 
AliBasedNdetaTask::SetDebugLevel(Int_t lvl)
{
  AliAnalysisTaskSE::SetDebugLevel(lvl);
  for (Int_t i = 0; i < fListOfCentralities->GetEntries(); i++) { 
    CentralityBin* bin = 
      static_cast<CentralityBin*>(fListOfCentralities->At(i));
    bin->SetDebugLevel(lvl);
  }
}

//________________________________________________________________________
void 
AliBasedNdetaTask::AddCentralityBin(UShort_t at, Float_t low, Float_t high)
{
  // 
  // Add a centrality bin 
  // 
  // Parameters:
  //    low  Low cut
  //    high High cut
  //
  DGUARD(fDebug,3,"Add a centrality bin [%6.2f,%6.2f] @ %d", low, high, at);
  CentralityBin* bin = MakeCentralityBin(GetName(), low, high);
  if (!bin) { 
    Error("AddCentralityBin", 
	  "Failed to create centrality bin for %s [%6.2f,%6.2f] @ %d", 
	  GetName(), low, high, at);
    return;
  }
  bin->SetSatelliteVertices(fSatelliteVertices);
  bin->SetDebugLevel(fDebug);
  Color_t color = (bin->IsAllBin() ? kBlack : kGray+1);
  if (HasCentrality() && !(high <= low) && (high <= 100)) {
    if  (fCentAxis.GetNbins() > 7) color = PbPbColor(low,high);
    else                           color = pPbColor (low,high);
  }
  DMSG(fDebug, 3, "Color of bin %d", color);
  bin->SetColor(color);
  fListOfCentralities->AddAtAndExpand(bin, at);
}

//________________________________________________________________________
AliBasedNdetaTask::CentralityBin*
AliBasedNdetaTask::MakeCentralityBin(const char* name, 
				     Float_t low, Float_t high) const
{
  // 
  // Make a centrality bin 
  // 
  // Parameters:
  //    name  Name used for histograms
  //    low   Low cut in percent
  //    high  High cut in percent
  // 
  // Return:
  //    A newly created centrality bin 
  //
  DGUARD(fDebug,3,"Make a centrality bin %s [%6.2f,%6.2f]", name, low, high);
  return new CentralityBin(name, low, high);
}

#define TESTAPPEND(SCHEME,BIT,STRING) \
  do { if (!(SCHEME & BIT)) break;					\
    if (!s.IsNull()) s.Append(","); s.Append(STRING); } while(false) 
  
//________________________________________________________________________
const Char_t*
AliBasedNdetaTask::NormalizationSchemeString(UShort_t scheme)
{
  // Create a string from normalization scheme bits 
  static TString s;
  s = "";

  if (scheme == kNone) 
    return s.Data();
  if (scheme == kFull) { 
    s = "FULL";
    return s.Data();
  }
  TESTAPPEND(scheme, kEventLevel, 	 "EVENT");
  TESTAPPEND(scheme, kBackground, 	 "BACKGROUND");
  TESTAPPEND(scheme, kTriggerEfficiency, "TRIGGER");
  TESTAPPEND(scheme, kZeroBin, 		 "ZEROBIN");

  return s.Data();
}
//________________________________________________________________________
UShort_t
AliBasedNdetaTask::ParseNormalizationScheme(const char* what)
{
  UShort_t    scheme = 0;
  TString     twhat(what);
  twhat.ToUpper();
  if (twhat.EqualTo("DEFAULT"))
    return kTriggerEfficiency|kEventLevel;
  TObjString* opt;
  TObjArray* token = twhat.Tokenize(" ,|");
  TIter       next(token);
  while ((opt = static_cast<TObjString*>(next()))) { 
    TString s(opt->GetString());
    if      (s.IsNull()) continue;
    Bool_t add = true;
    switch (s[0]) { 
    case '-': add = false; // Fall through 
    case '+': s.Remove(0,1);  // Remove character 
    }
    UShort_t bit = 0;
    if (s.EqualTo("SHAPE")) { 
      AliWarningGeneral("AliBasedNdetaTask",
			Form("SHAPE correction no longer supported (%s)",
			     what));
      continue;
    }
    if      (s.CompareTo("EVENT")     == 0) bit = kEventLevel;
    else if (s.CompareTo("BACKGROUND")== 0) bit = kBackground;
    else if (s.CompareTo("TRIGGER")   == 0) bit = kTriggerEfficiency;
    else if (s.CompareTo("FULL")      == 0) bit = kFull;
    else if (s.CompareTo("NONE")      == 0) bit = kNone;
    else if (s.CompareTo("ZEROBIN")   == 0) bit = kZeroBin;
    else 
      ::Warning("SetNormalizationScheme", "Unknown option %s", s.Data());
    if (add) scheme |= bit;
    else     scheme ^= bit;
  }
  delete token;
  return scheme;
}  
//________________________________________________________________________
void 
AliBasedNdetaTask::SetNormalizationScheme(const char* what)
{
  // 
  // Set normalisation scheme 
  // 
  DGUARD(fDebug,3,"Set the normalization scheme: %s", what);
  SetNormalizationScheme(ParseNormalizationScheme(what));
}
//________________________________________________________________________
void 
AliBasedNdetaTask::SetNormalizationScheme(UShort_t scheme) 
{
  DGUARD(fDebug,3,"Set the normalization scheme: 0x%x", scheme);
  fNormalizationScheme = scheme; 
}
//____________________________________________________________________
Bool_t
AliBasedNdetaTask::SetCentralityMethod(const TString& method)
{
  Int_t id = GetCentMethodID(method);
  if (id < -1) {
    AliErrorF("Unknown centrality estimator: %s", method.Data());
    fCentMethod = "";
    return false;
  }
  if (id < 0) { 
    // Do not use any estimator 
    AliInfoF("No centrality estimator: \"%s\"", method.Data());
    fCentMethod = "";
    return false;
  }

  TString meth = GetCentMethod(id);
  if (fName.Contains("Forward", TString::kIgnoreCase) && meth.Contains("FMD")) 
    AliWarningF("Centrality estimator %s used by %s - beware of auto-corr",
		meth.Data(), fName.Data());  
  else if (fName.Contains("Central", TString::kIgnoreCase) && 
	   (meth.Contains("CL0") || meth.Contains("TKL")))
    AliWarningF("Centrality estimator %s used by %s - beware of auto-corr",
		meth.Data(), fName.Data());

  fCentMethod = meth;
  AliInfoF("Centrality estimator set to %s", fCentMethod.Data());
  return true;
}

//________________________________________________________________________
Int_t
AliBasedNdetaTask::GetCentMethodID(const TString& meth)
{
  Int_t ret = -2;
  TString m(meth);
  m.ToUpper();
  if (m.EqualTo("NONE") || m.EqualTo("NO") || m.EqualTo("FALSE"))
    ret = kCentNone;
  else if (m.IsNull())               ret = kCentDefault;
  else if (m.BeginsWith("DEFAULT"))  ret = kCentDefault;
  else if (m.BeginsWith("ZEMVSZDC")) ret = kCentZEMvsZDC; 
  else if (m.BeginsWith("TKLVSV0M")) ret = kCentTklvsV0M; 
  else if (m.BeginsWith("V0MVSFMD")) ret = kCentV0MvsFMD; 
  else if (m.BeginsWith("NPA"))      ret = kCentNPA;
  else if (m.BeginsWith("ZPC"))      ret = kCentZPC;
  else if (m.BeginsWith("ZPA"))      ret = kCentZPA;
  else if (m.BeginsWith("ZNC"))      ret = kCentZNC;
  else if (m.BeginsWith("ZNA"))      ret = kCentZNA;
  else if (m.BeginsWith("CND"))      ret = kCentCND;
  else if (m.BeginsWith("CL1"))      ret = kCentCL1;
  else if (m.BeginsWith("CL0"))      ret = kCentCL0;
  else if (m.BeginsWith("TKL"))      ret = kCentTkl;
  else if (m.BeginsWith("TRK"))      ret = kCentTrk;
  else if (m.BeginsWith("FMD"))      ret = kCentFMD;
  else if (m.BeginsWith("V0C"))      ret = kCentV0C;
  else if (m.BeginsWith("V0A123"))   ret = kCentV0A123;
  else if (m.BeginsWith("V0A"))      ret = kCentV0A;
  else if (m.BeginsWith("V0M"))      ret = kCentV0M;
  else if (m.BeginsWith("MULTV0A"))  ret = kMultV0A;
  else if (m.BeginsWith("MULTV0M"))  ret = kMultV0M;
  else if (m.BeginsWith("MULTV0C"))  ret = kMultV0C;
  else if (m.BeginsWith("MULT"))     ret = kMult;
  if (m.Contains("TRUE"))            ret |= kCentTrue;
  if (m.Contains("EQ"))              ret |= kCentEq;
  
  return ret;
}
//________________________________________________________________________
const char*
AliBasedNdetaTask::GetCentMethod(UShort_t id)
{
  static TString ret("");
  UShort_t base = (id & 0xFF);
  switch (base) {
  case kCentNone:       ret = "none";           break;
  case kCentDefault:    ret = "";               break;
  case kCentV0M:	ret = "V0M";		break;
  case kCentV0A:	ret = "V0A";		break;
  case kCentV0A123:	ret = "V0A123";		break;
  case kCentV0C:	ret = "V0C";		break;
  case kCentFMD:	ret = "FMD";		break;
  case kCentTrk:	ret = "TRK";		break;
  case kCentTkl:	ret = "TKL";		break;
  case kCentCL0:	ret = "CL0";		break;
  case kCentCL1:	ret = "CL1";		break;
  case kCentCND:	ret = "CND";		break;
  case kCentZNA:	ret = "ZNA";		break;
  case kCentZNC:	ret = "ZNC";		break;
  case kCentZPA:	ret = "ZPA";		break;
  case kCentZPC:	ret = "ZPC";		break;
  case kCentNPA:	ret = "NPA";		break;
  case kCentV0MvsFMD:	ret = "V0MvsFMD";	break;
  case kCentTklvsV0M:	ret = "TKLvsV0M";	break;
  case kCentZEMvsZDC:	ret = "ZEMvsZDC";	break;
  case kMult:           ret = "MULT";           break;
  case kMultV0A:        ret = "MULTV0A";        break;
  case kMultV0M:        ret = "MULTV0M";        break;
  case kMultV0C:        ret = "MULTV0C";        break;
  default:              ret = "";               break;
  }
  Bool_t tru = id & kCentTrue;
  Bool_t eq  = id & kCentEq;
  if (eq) { 
    if (!tru) ret.Append("eq");
    else      ret.Append("Eq");
  }
  if (tru)    ret.Append("true");
  
  return ret.Data();
}


//________________________________________________________________________
void 
AliBasedNdetaTask::InitializeCentBins()
{
  if (fListOfCentralities->GetEntries() > 0) return;

  // Automatically add 'all' centrality bin if nothing has been defined. 
  AddCentralityBin(0, 0, 0);

  // Check if the centrality method was set to none, and in that case
  // remove the centrality axis.
  if (fCentMethod.EqualTo("none", TString::kIgnoreCase)) {
    fCentAxis.Set(0,0,0);
    TH1* h = static_cast<TH1*>(fSums->FindObject(fCentAxis.GetName()));
    if (h) {
      Info("InitializeCentBins",
	   "No centrality, forcing centrality axis to null");
      h->GetXaxis()->Set(0,0,0);
      h->Rebuild();
    }
  }
  if (HasCentrality()) { 
    const TArrayD* bins = fCentAxis.GetXbins();
    Int_t          nbin = fCentAxis.GetNbins(); 
    for (Int_t i = 0; i < nbin; i++) 
      AddCentralityBin(i+1,  (*bins)[i], (*bins)[i+1]);
    AddCentralityBin(nbin+1, 0, 101);
  }
}

//________________________________________________________________________
Bool_t
AliBasedNdetaTask::Book()
{
  // 
  // Create output objects.  
  //
  // This is called once per slave process 
  //
  DGUARD(fDebug,1,"Create user ouput object");

  // Modify axis of centrality histograms 
  fCent->SetXTitle(Form("Centrality (%s) [%%]", fCentMethod.Data()));
  fAccCent->SetXTitle(Form("Centrality (%s) [%%]", fCentMethod.Data()));

  fSums->Add(AliForwardUtil::MakeParameter("empirical", 
					   fEmpiricalCorrection != 0));
  fSums->Add(AliForwardUtil::MakeParameter("scheme", fNormalizationScheme));
  fSums->Add(AliForwardUtil::MakeParameter("centEstimator", 
					   GetCentMethodID(fCentMethod)));
  if (fIpzReweight)  fSums->Add(fIpzReweight->Clone("ipZw")); 
  // fSums->Add(new TNamed("centEstimator", fCentMethod.Data()));

  // Make our centrality bins 
  InitializeCentBins();

  // Loop over centrality bins 
  TIter next(fListOfCentralities);
  CentralityBin* bin = 0;
  while ((bin = static_cast<CentralityBin*>(next()))) 
    bin->CreateOutputObjects(fSums, fTriggerMask);
  

  if (HasCentrality()) {
    if (fCentAxis.GetXbins()->GetArray()) {
      Int_t nbin = fCentAxis.GetNbins();
      TArrayD a(nbin+1, fCentAxis.GetXbins()->GetArray());
      a.Set(nbin+2);
      a[nbin+1] = fCentAxis.GetXmax()+fCentAxis.GetBinWidth(nbin+1);
      
      fSeenCent = new TH1D("centSeen", "Centralities seen",
			    a.GetSize()-1, a.GetArray());
      fMeanVsC=new TProfile("sumVsC",
			    "Integral vs centrality",
			    a.GetSize()-2, a.GetArray());	
    }
    else {
      fSeenCent = new TH1D("centSeen", "Centralities seen",
			   fCentAxis.GetNbins()+1,
			   fCentAxis.GetXmin(),
			   fCentAxis.GetXmax()+fCentAxis.GetBinWidth(1));
      fMeanVsC=new TProfile("sumVsC",
			    "Integral vs centrality",
			    fCentAxis.GetNbins()+1,
			    fCentAxis.GetXmin(),
			    fCentAxis.GetXmax()+fCentAxis.GetBinWidth(1));
    }
  }
  else {
    fSeenCent = new TH1D("centSeen", "Null",1,0,100);
    fMeanVsC = new TProfile("sumVsC", "Null",1,0,100);
  }
  fSeenCent->SetDirectory(0);
  fSeenCent->SetXTitle(Form("Centrality (%s) [%%]",fCentMethod.Data()));
  fSeenCent->SetYTitle("Events");
  fSeenCent->SetFillStyle(3004);
  fSeenCent->SetFillColor(kYellow-2);
  fMeanVsC->SetDirectory(0);
  fMeanVsC->SetXTitle(fSeenCent->GetXaxis()->GetTitle());
  fMeanVsC->SetYTitle("#LT#int signal#GT");
  fMeanVsC->SetMarkerStyle(20);
  fMeanVsC->SetMarkerColor(kRed+1);
  fMeanVsC->SetLineColor(kRed+1);
  

  fTakenCent = static_cast<TH1D*>(fSeenCent->Clone("centTaken"));
  fTakenCent->SetTitle("Centralities taken by bins");
  fTakenCent->SetFillColor(kCyan-2);
  fTakenCent->SetFillStyle(3005);
  fTakenCent->SetDirectory(0);
  
  fSums->Add(fSeenCent);
  fSums->Add(fTakenCent);
  fSums->Add(fMeanVsC);

  // fSums->ls();
  return true;
}

//____________________________________________________________________
Bool_t
AliBasedNdetaTask::CheckEvent(const AliAODForwardMult& fwd) 
{
  if (!fCentMethod.BeginsWith("MULT")) 
    AliBaseAODTask::CheckEvent(fwd);
  else
    fwd.CheckEvent(fTriggerMask,
		   fIPzAxis.GetXmin(),
		   fIPzAxis.GetXmax(),
		   0,
		   0,
		   fTriggers,
		   fEventStatus,
		   fFilterMask);

  // Here, we always return true, as the centrality bins will do 
  // their own checks on the events - this is needed for event 
  // normalization etc. 
  return true;
}

//____________________________________________________________________
Double_t
AliBasedNdetaTask::GetCentrality(AliAODEvent&       event,
				 AliAODForwardMult* forward,
				 Int_t&             qual)
{
  DGUARD(fDebug,2,"Getting centrality from event of object: %s",
	 fCentMethod.Data());
  if (fCacheCent > -1) {
    // In case we already got the centrality, don't do it again
    DMSG(fDebug,1,"Returning cached value: %5.1f%%", fCacheCent);
    qual = fCacheQual;
    return fCacheCent;
  }
  if (fCentMethod.EqualTo("none")) {
    fCacheQual = 0;
    return fCacheCent = 0;
  }
  qual            = 0;
  Float_t cent    = AliBaseAODTask::GetCentrality(event,forward,qual);
  DMSG(fDebug,1,"Centrality stored in AOD forward: %5.1f%%", cent);
  if (!fCentMethod.IsNull()) {
    // Clear bad cent if already set
    cent = AliForwardUtil::GetCentrality(event,fCentMethod,qual,(fDebug > 1));
    if      (cent < 0)                   cent = -.5; // Bad centrality 
    else if (TMath::Abs(cent-100) < 1.1) cent = 100; // Special centralities
    DMSG(fDebug,1,"Centrality from mult: %5.1f%% (%d)", cent, qual);
  }
  fCacheQual        = qual;
  return fCacheCent = cent;  
}
//____________________________________________________________________
Double_t
AliBasedNdetaTask::GetCentrality(AliAODEvent& event,
				 AliAODForwardMult* forward)
{
  Int_t    qual = 0;
  Double_t cent = GetCentrality(event, forward, qual);
  if (qual > 0)   forward->SetTriggerBits(AliAODForwardMult::kCentNoCalib);
  return cent;
}
  
//____________________________________________________________________
Bool_t
AliBasedNdetaTask::Event(AliAODEvent& aod) 
{
  // 
  // Process a single event 
  // 
  // Parameters:
  //    option Not used
  //
  // Main loop
  DGUARD(fDebug,1,"Analyse the AOD event");
  if (fUseUtilPileup && fAnaUtil.IsPileUpEvent(&aod)) return false;

  AliAODForwardMult* forward = GetForward(aod);
  if (!forward) return false;
  
  // Fill centrality histogram 
    
  Double_t vtx    = forward->GetIpZ();
  TH2D*    data   = GetHistogram(aod, false);
  TH2D*    dataMC = GetHistogram(aod, true);
  if (!data) return false;

  CheckEventData(vtx, data, dataMC);
  
  if (!ApplyEmpiricalCorrection(forward,data))
    return false;

  Double_t ipzW = 1;
  if (fIpzReweight)  {
    ipzW = fIpzReweight->Eval(vtx);
    DMSG(fDebug,5,"IPz=%f -> Weight %f", vtx, ipzW);
  }
  
  Bool_t isZero = ((fNormalizationScheme & kZeroBin) &&
		   !forward->IsTriggerBits(AliAODForwardMult::kNClusterGt0));
  Bool_t taken  = false;
  // Loop over centrality bins 
  CentralityBin* allBin = 
    static_cast<CentralityBin*>(fListOfCentralities->At(0));
  if (allBin->ProcessEvent(forward, fTriggerMask, isZero,
			   fIPzAxis.GetXmin(), fIPzAxis.GetXmax(),
			   data, dataMC, fFilterMask, ipzW)) taken = true;
  
  // Find this centrality bin
  if (HasCentrality()) {
    taken                  = false;
    // After this call, the if the event isn't covered by the
    // centrality calibrations, a flag will be set in the forward
    // object.
    Double_t       cent    = GetCentrality(aod, forward);
    fMeanVsC->Fill(cent, data->Integral());

    // fSeenCent->Fill(cent);
    DMSG(fDebug,1,"Got event centrality %f (%s)", cent,
	 forward->IsTriggerBits(AliAODForwardMult::kInclusive|
				AliAODForwardMult::kCentNoCalib) ?
	 "out-of-calib" : "in-calib");
    
    Int_t          icent   = fCentAxis.FindBin(cent);
    DMSG(fDebug,1,"Centrality %5.2f%% maps to bin %d", cent, icent);
    if (icent == (fCentAxis.GetNbins()+1) &&
	TMath::Abs(cent-fCentAxis.GetXmax()) < 1e-6)
      // Make sure that upper edge is analysed 
      icent = fCentAxis.GetNbins();
    CentralityBin* thisBin = 0;
    if (icent >= 1 && icent <= fCentAxis.GetNbins()) 
      thisBin = static_cast<CentralityBin*>(fListOfCentralities->At(icent));
    if (thisBin && thisBin->ProcessEvent(forward, fTriggerMask,
					 isZero,
					 fIPzAxis.GetXmin(),
					 fIPzAxis.GetXmax(), 
					 data, dataMC, fFilterMask, ipzW)) 
      taken = true;
    if (taken) fTakenCent->Fill(cent);
    
    Int_t nbins = fCentAxis.GetNbins();
    CentralityBin* fullBin = 
      static_cast<CentralityBin*>(fListOfCentralities->At(nbins+1));
    if (fullBin && fullBin->ProcessEvent(forward, fTriggerMask, isZero,
					 fIPzAxis.GetXmin(),
					 fIPzAxis.GetXmax(), 
					 data, dataMC,
					 fFilterMask, ipzW))
      fSeenCent->Fill(cent);
  }
  
  return taken;
}

//________________________________________________________________________
void AliBasedNdetaTask::CheckEventData(Double_t,
				       TH2*,
				       TH2*) 
{
}

//________________________________________________________________________
void 
AliBasedNdetaTask::SetHistogramAttributes(TH1D* h, Int_t colour, Int_t marker,
					  const char* title, const char* ytitle)
{
  // 
  // Set histogram graphical options, etc. 
  // 
  // Parameters:
  //    h       Histogram to modify
  //    colour  Marker color 
  //    marker  Marker style
  //    title   Title of histogram
  //    ytitle  Title on y-axis. 
  //
  h->SetTitle(title);
  h->SetMarkerColor(colour);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(marker == 29 || marker == 30 ? 1.2 : 1);
  h->SetFillStyle(0);
  TString ytit;
  if (ytitle && ytitle[0] != '\0') ytit = ytitle;
  ytit = "#frac{1}{#it{N}}#frac{d#it{N}_{ch}}{d#it{#eta}}";
  h->SetYTitle(ytit);
  h->GetXaxis()->SetTitleFont(132);
  h->GetXaxis()->SetLabelFont(132);
  h->GetXaxis()->SetNdivisions(10);
  h->GetYaxis()->SetTitleFont(132);
  h->GetYaxis()->SetLabelFont(132);
  h->GetYaxis()->SetNdivisions(10);
  h->GetYaxis()->SetDecimals();
  h->SetStats(0);
}

//________________________________________________________________________
void
AliBasedNdetaTask::ScaleToCoverage(TH2D* copy, const TH1D* norm) 
{
  // Normalize to the acceptance -
  // dndeta->Divide(accNorm);
  for (Int_t i = 1; i <= copy->GetNbinsX(); i++) { 
    Double_t a = norm->GetBinContent(i);
    for (Int_t j = 1; j <= copy->GetNbinsY(); j++) { 
      if (a <= 0) { 
	copy->SetBinContent(i,j,0);
	copy->SetBinError(i,j,0);
	continue;
      }
      Double_t c = copy->GetBinContent(i, j);
      Double_t e = copy->GetBinError(i, j);
      copy->SetBinContent(i, j, c / a);
      copy->SetBinError(i, j, e / a);
    }
  }
}
//________________________________________________________________________
void
AliBasedNdetaTask::ScaleToCoverage(TH1D* copy, const TH1D* norm) 
{
  // Normalize to the acceptance -
  // dndeta->Divide(accNorm);
  for (Int_t i = 1; i <= copy->GetNbinsX(); i++) { 
    Double_t a = norm->GetBinContent(i);
    if (a <= 0) { 
      copy->SetBinContent(i,0);
      copy->SetBinError(i,0);
      continue;
    }
    Double_t c = copy->GetBinContent(i);
    Double_t e = copy->GetBinError(i);
    copy->SetBinContent(i, c / a);
    copy->SetBinError(i, e / a);
  }
}

//________________________________________________________________________
TH1D*
AliBasedNdetaTask::ProjectX(const TH2D* h, 
			    const char* name,
			    Int_t firstbin, 
			    Int_t lastbin, 
			    bool  useRoot,
			    bool  corr,
			    bool  error)
{
  // 
  // Project onto the X axis 
  // 
  // Parameters:
  //    h         2D histogram 
  //    name      New name 
  //    firstbin  First bin to use 
  //    lastbin   Last bin to use
  //    error     Whether to calculate errors
  // 
  // Return:
  //    Newly created histogram or null
  //
  if (!h) return 0;
  if (useRoot) 
    return h->ProjectionX(name, firstbin, lastbin, (error ? "e" : ""));
  
  const TAxis* xaxis = h->GetXaxis();
  const TAxis* yaxis = h->GetYaxis();
  TH1D*  ret   = new TH1D(name, h->GetTitle(), xaxis->GetNbins(), 
			  xaxis->GetXmin(), xaxis->GetXmax());
  static_cast<const TAttLine*>(h)->Copy(*ret);
  static_cast<const TAttFill*>(h)->Copy(*ret);
  static_cast<const TAttMarker*>(h)->Copy(*ret);
  ret->GetXaxis()->ImportAttributes(xaxis);

  Int_t first = firstbin; 
  Int_t last  = lastbin;
  if      (first < 0)                    first = 1;
  else if (first >= yaxis->GetNbins()+2) first = yaxis->GetNbins()+1;
  if      (last  < 0)                    last  = yaxis->GetNbins();
  else if (last  >= yaxis->GetNbins()+2) last  = yaxis->GetNbins()+1;
  if (last-first < 0) { 
    AliWarningGeneral("AliBasedNdetaTask", 
		      Form("Nothing to project [%d,%d]", first, last));
    return 0;
    
  }

  // Loop over X bins 
  //DMSG(fDebug,3,"Projecting bins [%d,%d] of %s", first, last, h->GetName()));
  Int_t ybins = (last-first+1);
  for (Int_t xbin = 0; xbin <= xaxis->GetNbins()+1; xbin++) { 
    Double_t content = 0;
    Double_t error2  = 0;
    Int_t    nbins   = 0;
    
    
    for (Int_t ybin = first; ybin <= last; ybin++) { 
      Double_t c1 = h->GetBinContent(h->GetBin(xbin, ybin));
      Double_t e1 = h->GetBinError(h->GetBin(xbin, ybin));

      // Ignore empty bins 
      if (c1 < 1e-12) continue;
      if (e1 < 1e-12) {
	if (error) continue; 
	e1 = 1;
      }

      content    += c1;
      error2     += e1*e1;
      nbins++;
    } // for (ybin)
    if(content > 0 && nbins > 0) {
      Double_t factor = (corr ? Double_t(ybins) / nbins : 1);
#if 0
      AliWarningGeneral(ret->GetName(), 
			Form("factor @ %d is %d/%d -> %f", 
			     xbin, ybins, nbins, factor));
#endif
      if (error) { 
	// Calculate weighted average
	ret->SetBinContent(xbin, content * factor);
	ret->SetBinError(xbin, factor * TMath::Sqrt(error2));
      }
      else 
	ret->SetBinContent(xbin, factor * content);
    }
  } // for (xbin)
  
  return ret;
}
 
//________________________________________________________________________
Bool_t 
AliBasedNdetaTask::Finalize() 
{
  // 
  // Called at end of event processing.. 
  //
  // This is called once in the master 
  // 
  // Parameters:
  //    option Not used 
  //
  // Draw result to screen, or perform fitting, normalizations Called
  // once at the end of the query
  DGUARD(fDebug,1,"Process final merged results");

  UShort_t sNN;
  UShort_t sys; 
  ULong_t  trig;
  ULong_t  filter;
  UShort_t scheme;
  Int_t    centID;
  AliForwardUtil::GetParameter(fSums->FindObject("sNN"), sNN);
  AliForwardUtil::GetParameter(fSums->FindObject("sys"), sys);
  AliForwardUtil::GetParameter(fSums->FindObject("trigger"), trig);
  AliForwardUtil::GetParameter(fSums->FindObject("filter"), filter);
  AliForwardUtil::GetParameter(fSums->FindObject("scheme"), scheme);
  AliForwardUtil::GetParameter(fSums->FindObject("centEstimator"), centID);

  TH1*   cH = static_cast<TH1*>(fSums->FindObject("centAxis"));
  Info("", "centAxis: %p (%s)", cH, (cH ? cH->ClassName() : "null"));
  // TAxis* cA = static_cast<TAxis*>(fSums->FindObject("centAxis"));
  TAxis*  cA = (cH ? cH->GetXaxis() : 0);
  if (cA) cA->Copy(fCentAxis);
  fCentAxis.SetName("centAxis");
  fCentAxis.SetTitle("Centrality [%]");
  
  // (Re)Create our centrality bins 
  InitializeCentBins();

  // Print before we loop
  Print();

  // Loop over centrality bins 
  TIter next(fListOfCentralities);
  CentralityBin* bin = 0;
  gStyle->SetPalette(1);
  THStack* dndetaStack        = new THStack("dndeta", "dN/d#eta");
  THStack* dndetaMCStack      = new THStack("dndetaMC", "dN_{ch}/d#eta");
  THStack* dndetaEmpStack     = new THStack("dndetaEmp", "dN_{ch}/d#eta");
  THStack* leftRightStack     = new THStack("leftRight", "#eta>0/#eta<0");
  
  TList* mclist = 0;
  TList* truthlist = 0;
  
  if (fFinalMCCorrFile.Contains(".root")) {
    AliForwardUtil::SuppressGuard g;
    TFile* ftest = TFile::Open(fFinalMCCorrFile.Data());
    if(ftest) {
      mclist    = dynamic_cast<TList*>(ftest->Get(Form("%sResults",GetName())));
      truthlist = dynamic_cast<TList*>(ftest->Get("MCTruthResults"));
    }
    else 
      AliWarning("MC analysis file invalid - no final MC correction possible");
  }
  Int_t style = GetMarker();
  Int_t color = GetColor();
  
  DMSG(fDebug,3,"Marker style=%d, color=%d", style, color);
  while ((bin = static_cast<CentralityBin*>(next()))) {
    if (!bin->End(fSums, fResults, fNormalizationScheme,
		 fTriggerEff, fTriggerEff0, 
		 fUseROOTProj, fCorrEmpty, fTriggerMask, 
		 style, color, mclist, truthlist))
      continue;
    if (HasCentrality() && bin->IsAllBin()) {
      // If we have centrality bins, do not add the min-bias
      // distribution to the output stack.
      AliWarning("Skipping MB bin since we have centrality");
      continue;
    }
#if 0
    if (HasCentrality() && bin->IsInclusiveBin()) {
      AliWarning("Skipping 0-100 bin, as we have centrality");
      continue;
    }
#endif 
    TH1* dndeta      = bin->GetResult("");
    TH1* dndetaMC    = bin->GetResult("MC", false);
    TH1* dndetaEmp   = bin->GetResult("Emp", false);
    TH1* leftRight   = Asymmetry(dndetaEmp ? dndetaEmp : dndeta);
    // if (leftRight) bin->fOutput->Add(leftRight);
    DMSG(fDebug,2,"Results: bare=%p mcbare=%p", dndeta, dndetaMC);
    if (dndeta)      dndetaStack->Add(dndeta);
    if (dndetaMC)    dndetaMCStack->Add(dndetaMC);
    if (dndetaEmp)   dndetaEmpStack->Add(dndetaEmp);
    if (leftRight)   leftRightStack->Add(leftRight);
  }
  // Output the stack
  fResults->Add(dndetaStack);

  // If available, output track-ref stack
  if (!dndetaMCStack->GetHists() || 
      dndetaMCStack->GetHists()->GetEntries() <= 0) {
    // AliWarning("No MC histograms found");
    delete dndetaMCStack;
    dndetaMCStack = 0;
  }
  if (dndetaMCStack) fResults->Add(dndetaMCStack);

  // If available, output track-ref stack
  DMSG(0,fDebug,"Emp stack: %p (%d)", dndetaEmpStack,
       dndetaEmpStack->GetHists() && dndetaEmpStack->GetHists()->GetEntries()
       ? dndetaEmpStack->GetHists()->GetEntries() : -1);
  if (!dndetaEmpStack->GetHists() || 
      dndetaEmpStack->GetHists()->GetEntries() <= 0) {
    // AliWarning("No MC histograms found");
    delete dndetaEmpStack;
    dndetaEmpStack = 0;
  }
  if (dndetaEmpStack) fResults->Add(dndetaEmpStack);
  

  // If available, output track-ref stack
  if (!leftRightStack->GetHists() || 
      leftRightStack->GetHists()->GetEntries() <= 0) {
    // AliWarning("No MC histograms found");
    delete leftRightStack;
    leftRightStack = 0;
  }
  if (leftRightStack) fResults->Add(leftRightStack);

  // Output collision energy string 
  if (sNN > 0) {
    TNamed* sNNObj = new TNamed("sNN", 
				AliForwardUtil::CenterOfMassEnergyString(sNN));
    sNNObj->SetUniqueID(sNN);
    fResults->Add(sNNObj); 
  }

  // Output collision system string 
  if (sys > 0) { 
    TNamed* sysObj = new TNamed("sys", 
				AliForwardUtil::CollisionSystemString(sys));
    sysObj->SetUniqueID(sys);
    fResults->Add(sysObj); 
  }

  // Output centrality axis
  TNamed* centEstimator = new TNamed("centEstimator", fCentMethod.Data());
  centEstimator->SetUniqueID(centID);
  fResults->Add(&fCentAxis);
  fResults->Add(centEstimator);
  fResults->Add(AliForwardUtil::MakeParameter("absMinCent", fAbsMinCent));
  
  // Output trigger string 
  if (trig) {
    TString tstr = AliAODForwardMult::GetTriggerString(trig);
    if (fTriggerEff > 0 && fTriggerEff < 1) {
      if      (tstr.EqualTo("MBOR"))  tstr = "INEL";
      else if (tstr.EqualTo("V0AND")) tstr = "NSD";
    }
    else if (tstr.EqualTo("V0AND"))   tstr = "VISX";
    TNamed* maskObj = new TNamed("trigger",tstr.Data());
    maskObj->SetUniqueID(trig);
    fResults->Add(maskObj); 
  }
  // Output filter string 
  if (filter) { 
    TNamed* maskObj = new TNamed("filter",
				 AliAODForwardMult::GetTriggerString(filter,
								     " | "));
    maskObj->SetUniqueID(filter);
    fResults->Add(maskObj); 
  }
  
  // Normalization string 
  if (scheme > 0) {
    TNamed* schemeObj = new TNamed("scheme",
				   NormalizationSchemeString(scheme));
    schemeObj->SetUniqueID(scheme);
    fResults->Add(schemeObj);
  }

  // Output vertex axis 
  fIPzAxis.SetName("vtxAxis");
  fIPzAxis.SetTitle(Form("v_{z}#in[%+5.1f,%+5.1f]cm",
			 fIPzAxis.GetXmin(),
			 fIPzAxis.GetXmax()));
  TAxis* ipzAxis = static_cast<TAxis*>(fIPzAxis.Clone());
  fResults->Add(ipzAxis);

  // Output trigger efficiency 
  fResults->Add(AliForwardUtil::MakeParameter("triggerEff", fTriggerEff));
  fResults->Add(AliForwardUtil::MakeParameter("triggerEff0", fTriggerEff0));

  TNamed* options = new TNamed("options","");
  TString str;
  str.Append(Form("Empty bins %scorrected for, ", fCorrEmpty ? "" : "not "));
  str.Append(Form("TH2::ProjectionX %sused", fUseROOTProj ? "" : "not "));
  options->SetTitle(str);
  fResults->Add(options);

  TObject* sumVsC = fSums->FindObject("sumVsC");
  if (sumVsC) {
    fMeanVsC = static_cast<TProfile*>(sumVsC->Clone());
    fMeanVsC->SetDirectory(0);
    fResults->Add(fMeanVsC);
  }
  
  return true;
}

#define PF(N,V,...)					\
  AliForwardUtil::PrintField(N,V, ## __VA_ARGS__)
#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)
#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)
#if 0
namespace {
  void appendBit(TString& str, const char* bit)
  {
    if (!str.IsNull()) str.Append("|");
    str.Append(bit);
  }  
}
#endif

//________________________________________________________________________
void 
AliBasedNdetaTask::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  AliBaseAODTask::Print(option);
  TString schemeString(NormalizationSchemeString(fNormalizationScheme));

  gROOT->IncreaseDirLevel();
  PFB("Use AnaUtil pile-up",     fUseUtilPileup);
  PFB("Use TH2::ProjectionX",    fUseROOTProj);
  PFB("Correct for empty",       fCorrEmpty);
  PFV("Normalization scheme",	 schemeString );
  PFV("Trigger efficiency",	 fTriggerEff);
  PFV("Bin-0 Trigger efficiency", fTriggerEff0);
  PFV("Centrality estimator",    (fCentMethod.IsNull() ? 
				  "-default-" : fCentMethod.Data()));
  PFB("Re-weight from IPz",      fIpzReweight!=0);
  if (fIpzReweight) fIpzReweight->Print();
    
  
  TString opt(option);
  opt.ToUpper();
  if (opt.Contains("R") &&
      fListOfCentralities &&
      fListOfCentralities->GetEntries() > 0) {
    TIter next(fListOfCentralities);
    TObject* bin = 0;
    while ((bin = next())) bin->Print(option);
  }
  gROOT->DecreaseDirLevel();  
}

//__________________________________________________________________
Int_t
AliBasedNdetaTask::GetMarkerStyle(UShort_t bits)
{
  Int_t  base   = bits & (0xFE);
  Bool_t hollow = bits & kHollow;
  switch (base) { 
  case kCircle:       return (hollow ? 24 : 20);
  case kSquare:       return (hollow ? 25 : 21);
  case kUpTriangle:   return (hollow ? 26 : 22);
  case kDownTriangle: return (hollow ? 32 : 23);
  case kDiamond:      return (hollow ? 27 : 33); 
  case kCross:        return (hollow ? 28 : 34); 
  case kStar:         return (hollow ? 30 : 29); 
  }
  return 1;
}
//__________________________________________________________________
UShort_t
AliBasedNdetaTask::GetMarkerBits(Int_t style) 
{ 
  UShort_t bits = 0;
  switch (style) { 
  case 24: case 25: case 26: case 27: case 28: case 30: case 32: 
    bits |= kHollow; break;
  }
  switch (style) { 
  case 20: case 24: bits |= kCircle;       break;
  case 21: case 25: bits |= kSquare;       break;
  case 22: case 26: bits |= kUpTriangle;   break;
  case 23: case 32: bits |= kDownTriangle; break;
  case 27: case 33: bits |= kDiamond;      break;
  case 28: case 34: bits |= kCross;        break;
  case 29: case 30: bits |= kStar;         break;
  }
  return bits;
}
//__________________________________________________________________
Int_t
AliBasedNdetaTask::FlipHollowStyle(Int_t style) 
{
  UShort_t bits = GetMarkerBits(style);
  Int_t    ret  = GetMarkerStyle(bits ^ kHollow);
  return ret;
}

//====================================================================
void
AliBasedNdetaTask::Sum::Init(TList* list, const TH2D* data, Int_t col)
{
  DGUARD(fDebug,1,"Initializing sums with %s", data->GetName());
  TString n(GetHistName(0));
  TString n0(GetHistName(1));
  const char* postfix = GetTitle();

  fSum = static_cast<TH2D*>(data->Clone(n));
  if (postfix) fSum->SetTitle(Form("%s (%s)", data->GetTitle(), postfix));
  fSum->SetDirectory(0);
  fSum->SetMarkerColor(col);
  fSum->SetMarkerStyle(GetMarkerStyle(kCircle|kSolid));
  fSum->Reset();
  list->Add(fSum);

  fSum0 = static_cast<TH2D*>(data->Clone(n0));
  if (postfix) 
    fSum0->SetTitle(Form("%s 0-bin (%s)", data->GetTitle(), postfix));
  else   
    fSum0->SetTitle(Form("%s 0-bin", data->GetTitle()));
  fSum0->SetDirectory(0);
  fSum0->SetMarkerColor(col);
  fSum0->SetMarkerStyle(GetMarkerStyle(kCross|kHollow));
  fSum0->Reset();
  list->Add(fSum0);

  fEvents = new TH1I(GetHistName(2), "Event types", 2, -.5, 1.5);
  fEvents->SetDirectory(0);
  fEvents->SetFillColor(kRed+1);
  fEvents->SetFillStyle(3002);
  fEvents->GetXaxis()->SetBinLabel(1, "Non-zero");
  fEvents->GetXaxis()->SetBinLabel(2, "Zero");
  list->Add(fEvents);
}

//____________________________________________________________________
TString
AliBasedNdetaTask::Sum::GetHistName(const char* /*name*/, 
				    Int_t what, const char* post)
{
  TString n;
  switch (what) { 
  case 0: n = "sum"; break;
  case 1: n = "sum0"; break;
  case 2: n = "events"; break;
  }
  if (post && post[0] != '\0')  n.Append(post);
  return n;
}

//____________________________________________________________________
TString
AliBasedNdetaTask::Sum::GetHistName(Int_t what) const
{
  return GetHistName(GetName(), what, GetTitle());
}

//____________________________________________________________________
void
AliBasedNdetaTask::Sum::Add(const TH2D* data, Bool_t isZero, Double_t weight)
{
  DGUARD(fDebug,2,"Adding %s to sums w/weight %f (%f)",
	 data->GetName(),weight,data->Integral());
  if (isZero) fSum0->Add(data, weight);
  else        fSum->Add(data, weight);
  fEvents->Fill(isZero ? 1 : 0);
}

//____________________________________________________________________
TH2D*
AliBasedNdetaTask::Sum::CalcSum(TList*       output, 
				Double_t&    ntotal,
				Double_t     epsilon0, 
				Double_t     epsilon,
				Int_t        marker,
				Bool_t       rootProj, 
				Bool_t       corrEmpty) const
{
  DGUARD(fDebug,2,"Calculating final summed histogram %s", fSum->GetName());

  // The return value `ret' is not scaled in anyway
  TH2D* ret      = static_cast<TH2D*>(fSum->Clone(fSum->GetName()));
  ret->SetDirectory(0);
  Int_t n        = Int_t(fEvents->GetBinContent(1));
  Int_t n0       = Int_t(fEvents->GetBinContent(2));
  ntotal         = n;
  if (n0 > 0) { 
    ret->Reset();
    DMSG(fDebug,1, 
	 "Adding histograms %s(%d) and %s(%d) w/weights %f and %f resp.",
	 fSum0->GetName(), n, fSum->GetName(), n0, 1./epsilon,1./epsilon0);
    ret->Add(fSum0, fSum, 1. / epsilon0, 1. / epsilon); 
    ntotal = n / epsilon + n0 / epsilon0;
  }

  TList* out = new TList;
  out->SetOwner();
  const char* postfix = GetTitle();
  if (!postfix) postfix = "";
  out->SetName(Form("partial%s", postfix));
  output->Add(out);

  // Now make copies, normalize them, and store in output list 
  // Note, these are the only ones normalized here
  // These are mainly for diagnostics 
  TH2D* sumCopy  = static_cast<TH2D*>(fSum->Clone("sum"));
  TH2D* sum0Copy = static_cast<TH2D*>(fSum0->Clone("sum0"));
  TH2D* retCopy  = static_cast<TH2D*>(ret->Clone("sumAll"));
  sumCopy->SetMarkerStyle(FlipHollowStyle(marker));
  sumCopy->SetDirectory(0);
  sum0Copy->SetMarkerStyle(GetMarkerStyle(GetMarkerBits(marker)+4));
  sum0Copy->SetDirectory(0);
  retCopy->SetMarkerStyle(marker);
  retCopy->SetDirectory(0);

  Int_t nY      = fSum->GetNbinsY();
  Int_t o       = 0; // nY+1;
  TH1D* norm    = ProjectX(fSum,  "norm",    o, o, rootProj, corrEmpty, false);
  TH1D* norm0   = ProjectX(fSum0, "norm0",   o, o, rootProj, corrEmpty, false);
  TH1D* normAll = ProjectX(ret,   "normAll", o, o, rootProj, corrEmpty, false);
  norm->SetTitle("#eta coverage - >0-bin");
  norm0->SetTitle("#eta coverage - 0-bin");
  normAll->SetTitle("#eta coverage");
  norm->SetDirectory(0);
  norm0->SetDirectory(0);
  normAll->SetDirectory(0);
  
  TH1D* sumCopyPx  = ProjectX(sumCopy,  "average",    1,nY,rootProj,corrEmpty);
  TH1D* sum0CopyPx = ProjectX(sum0Copy, "average0",   1,nY,rootProj,corrEmpty);
  TH1D* retCopyPx  = ProjectX(retCopy,  "averageAll", 1,nY,rootProj,corrEmpty);
  sumCopyPx-> SetTitle(Form("#sum_{i}^{N_{#phi}}%s", sumCopy->GetTitle()));
  sum0CopyPx->SetTitle(Form("#sum_{i}^{N_{#phi}}%s", sum0Copy->GetTitle()));
  retCopyPx-> SetTitle(Form("#sum_{i}^{N_{#phi}}%s", retCopy->GetTitle()));
  sumCopyPx-> SetDirectory(0);
  sum0CopyPx->SetDirectory(0);
  retCopyPx-> SetDirectory(0);

  TH1D* phi    = ProjectX(fSum,  "phi",    nY+1,nY+1,rootProj,corrEmpty,false);
  TH1D* phi0   = ProjectX(fSum0, "phi0",   nY+1,nY+1,rootProj,corrEmpty,false);
  TH1D* phiAll = ProjectX(ret,   "phiAll", nY+1,nY+1,rootProj,corrEmpty,false);
  phi   ->SetTitle("#phi acceptance from dead strips - >0-bin");
  phi0  ->SetTitle("#phi acceptance from dead strips - 0-bin");
  phiAll->SetTitle("#phi acceptance from dead strips");
  phi   ->SetDirectory(0);
  phi0  ->SetDirectory(0);
  phiAll->SetDirectory(0);

  const TH1D* cov    = (corrEmpty ? norm    : phi);
  const TH1D* cov0   = (corrEmpty ? norm0   : phi0);
  const TH1D* covAll = (corrEmpty ? normAll : phiAll);

  // Here, we scale to the coverage (or phi acceptance)
  ScaleToCoverage(sumCopy,  cov);
  ScaleToCoverage(sum0Copy, cov0);
  ScaleToCoverage(retCopy,  covAll);

  // Scale our 1D histograms
  sumCopyPx ->Scale(1., "width");
  sum0CopyPx->Scale(1., "width");  
  retCopyPx ->Scale(1., "width");  

  DMSG(fDebug,2,"Maximum %f,%f changed to %f", sumCopyPx->GetMaximum(), 
       sum0CopyPx->GetMaximum(), retCopyPx->GetMaximum());

  // Scale the normalization - they should be 1 at the maximum
  norm   ->Scale(n > 0   ? 1. / n  : 1);
  norm0  ->Scale(n0 > 0 ? 1. / n0 : 1);
  normAll->Scale(ntotal > 0 ? 1. / ntotal : 1);

  // Scale the normalization - they should be 1 at the maximum
  phi   ->Scale(n > 0   ? 1. / n  : 1);
  phi0  ->Scale(n0 > 0 ? 1. / n0 : 1);
  phiAll->Scale(ntotal > 0 ? 1. / ntotal : 1);

  out->Add(sumCopy);
  out->Add(sum0Copy);
  out->Add(retCopy);
  out->Add(sumCopyPx);
  out->Add(sum0CopyPx);
  out->Add(retCopyPx);
  out->Add(norm);
  out->Add(norm0);
  out->Add(normAll);
  out->Add(phi);
  out->Add(phi0);
  out->Add(phiAll);

  if (fDebug >= 1) {
    if (n0 > 0) 
      DMSG(fDebug,1,"Returning  (1/%f * %s + 1/%f * %s), "
	   "1/%f * %d + 1/%f * %d = %d", 
	   epsilon0, fSum0->GetName(), epsilon, fSum->GetName(), 
	   epsilon0, n0, epsilon, n, int(ntotal));
    else 
      DMSG(fDebug,1, "Returning %s, %d", fSum->GetName(), int(ntotal));
  }

#if 0
  for (Int_t i = 1; i <= ret->GetNbinsX(); i++) { 
    Double_t nc  = sum->GetBinContent(i, 0);
    Double_t nc0 = sum0->GetBinContent(i, 0);
    ret->SetBinContent(i, 0, nc + nc0); // Just count events 
  }
#endif
 
  return ret;
}
//____________________________________________________________________
void
AliBasedNdetaTask::Sum::Print(Option_t*) const
{
  PFV("dN/deta sum bin", GetName());
  gROOT->IncreaseDirLevel();
  PF("Normal sum", "%s (%p)", GetHistName(0).Data(), fSum);
  PF("0-bin sum",  "%s (%p)", GetHistName(1).Data(), fSum0);
  PF("Event count","%s (%p)", GetHistName(2).Data(), fEvents);
  gROOT->DecreaseDirLevel();
}

//====================================================================
AliBasedNdetaTask::CentralityBin::CentralityBin()
  : TNamed("", ""), 
    fSums(0), 
    fOutput(0),
    fSum(0), 
    fSumMC(0), 
    fTriggers(0), 
    fStatus(0),
    fLow(0), 
    fHigh(0),
    fDoFinalMCCorrection(false),
    fSatelliteVertices(false),
    fDebug(0)
{
  // 
  // Constructor 
  //
  DGUARD(fDebug,3,"Default CTOR of AliBasedNdeta::CentralityBin");
}
#define TRUNC(X) (Int_t(X) + Float_t(Int_t(X*100)%100)/100)

//____________________________________________________________________
AliBasedNdetaTask::CentralityBin::CentralityBin(const char* name, 
						Float_t low, Float_t high)
  : TNamed(name, ""), 
    fSums(0), 
    fOutput(0),
    fSum(0), 
    fSumMC(0), 
    fTriggers(0),
    fStatus(0),
    fLow(TRUNC(low)), 
    fHigh(TRUNC(high)),
    fDoFinalMCCorrection(false), 
    fSatelliteVertices(false),
    fDebug(0)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name used for histograms (e.g., Forward)
  //    low  Lower centrality cut in percent 
  //    high Upper centrality cut in percent 
  //
  DGUARD(fDebug,3,"Named CTOR of AliBasedNdeta::CentralityBin: "
	 "%s [%6.2f,%6.2f]",name, fLow, fHigh);
  if (low <= 0 && high <= 0) { 
    fLow  = 0; 
    fHigh = 0;
    SetTitle("All centralities");
  }
  else {
    fLow  = TRUNC(low);
    fHigh = TRUNC(high);
    SetTitle(Form("Centrality bin from %6.2f%% to %6.2f%%", low, high));
  }
}
//____________________________________________________________________
AliBasedNdetaTask::CentralityBin::CentralityBin(const CentralityBin& o)
  : TNamed(o), 
    fSums(o.fSums), 
    fOutput(o.fOutput),
    fSum(o.fSum), 
    fSumMC(o.fSumMC), 
    fTriggers(o.fTriggers), 
    fStatus(o.fStatus),
    fLow(o.fLow), 
    fHigh(o.fHigh),
    fDoFinalMCCorrection(o.fDoFinalMCCorrection),
    fSatelliteVertices(o.fSatelliteVertices),
    fDebug(o.fDebug)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    other Object to copy from 
  //
  DGUARD(fDebug,3,"Copy CTOR of AliBasedNdeta::CentralityBin");
}
//____________________________________________________________________
AliBasedNdetaTask::CentralityBin::~CentralityBin()
{
  // 
  // Destructor 
  //
  DGUARD(fDebug,3,"DTOR of AliBasedNdeta::CentralityBin");

  // if (fSums) fSums->Delete();
  // if (fOutput) fOutput->Delete();
}

//____________________________________________________________________
AliBasedNdetaTask::CentralityBin&
AliBasedNdetaTask::CentralityBin::operator=(const CentralityBin& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    other Object to assign from 
  // 
  // Return:
  //    Reference to this 
  //
  DGUARD(fDebug,3,"Centrality bin assignment");
  if (&o == this) return *this; 
  SetName(o.GetName());
  SetTitle(o.GetTitle());
  fSums      = o.fSums;
  fOutput    = o.fOutput;
  fSum       = o.fSum;
  fSumMC     = o.fSumMC;
  fTriggers  = o.fTriggers;
  fStatus    = o.fStatus;
  fLow       = o.fLow;
  fHigh      = o.fHigh;
  fDoFinalMCCorrection = o.fDoFinalMCCorrection;
  fSatelliteVertices = o.fSatelliteVertices;

  return *this;
}
#if 0
namespace {
  Color_t Brighten(Color_t origNum, Int_t nTimes=2)
  {
    TColor* col   = gROOT->GetColor(origNum);
    if (!col) return origNum;
    Int_t origR  = Int_t(0xFF * col->GetRed());
    Int_t origG  = Int_t(0xFF * col->GetGreen());
    Int_t origB  = Int_t(0xFF * col->GetBlue());
    Int_t off    = nTimes*0x33;
    Int_t newR   = TMath::Min((origR+off),0xff);
    Int_t newG   = TMath::Min((origG+off),0xff);
    Int_t newB   = TMath::Min((origB+off),0xff);
    Int_t newNum = TColor::GetColor(newR, newG, newB);
    return newNum;
  }
}
#endif
//____________________________________________________________________
Int_t
AliBasedNdetaTask::CentralityBin::GetColor(Int_t fallback) const 
{
  return fColor;
#if 0
  if (IsAllBin()) return kBlack; // fallback;
  Float_t  fc       = /*(fLow+double(fHigh-fLow)/2)*/ fHigh / 100;
  Int_t    nCol     = gStyle->GetNumberOfColors();
  Int_t    icol     = TMath::Min(nCol-1,int(fc * nCol + .5));
  Int_t    col      = gStyle->GetColorPalette(icol);
#if 0
  Color_t orig = col;
  col = Brighten(orig);
#endif
  return col;
#endif 
}

//____________________________________________________________________
const char* 
AliBasedNdetaTask::CentralityBin::GetListName() const
{
  // 
  // Get the list name 
  // 
  // Return:
  //    List Name 
  //
  if (IsAllBin()) return "all";
  return Form("cent%03dd%02d_%03dd%02d",
	      Int_t(fLow),  Int_t(fLow*100)  % 100, 
	      Int_t(fHigh), Int_t(fHigh*100) % 100);

}
//____________________________________________________________________
void
AliBasedNdetaTask::CentralityBin::CreateOutputObjects(TList* dir, Int_t mask)
{
  // 
  // Create output objects 
  // 
  // Parameters:
  //    dir   Parent list
  //
  DGUARD(fDebug,1,"Create centrality bin output objects");
  fSums = new TList;
  fSums->SetName(GetListName());
  fSums->SetOwner();
  dir->Add(fSums);

  fTriggers = AliAODForwardMult::MakeTriggerHistogram("triggers", mask);
  fTriggers->SetDirectory(0);

  fStatus = AliAODForwardMult::MakeStatusHistogram("status");
  fStatus->SetDirectory(0);

  fSums->Add(fTriggers);
  fSums->Add(fStatus);

  fSums->Add(new TParameter<float>("low",  fLow,  'f'));
  fSums->Add(new TParameter<float>("high", fHigh, 'f'));
}
//____________________________________________________________________
void
AliBasedNdetaTask::CentralityBin::SetDebugLevel(Int_t lvl)
{
  fDebug = lvl;
  if (fSum) fSum->fDebug = lvl;
  if (fSumMC) fSumMC->fDebug = lvl;
}

//____________________________________________________________________
Bool_t
AliBasedNdetaTask::CentralityBin::ReadSum(TList* list, bool mc)
{
  const char* post   = (mc ? "MC" : "");
  TString     sn     = Sum::GetHistName(GetName(),0,post);
  TString     sn0    = Sum::GetHistName(GetName(),1,post);
  TString     ev     = Sum::GetHistName(GetName(),2,post);
  TH2D*       sum    = static_cast<TH2D*>(list->FindObject(sn));
  TH2D*       sum0   = static_cast<TH2D*>(list->FindObject(sn0));
  TH1I*       events = static_cast<TH1I*>(list->FindObject(ev));
  if (!sum || !sum0 || !events) {
    if (!mc)
      AliWarningF("Failed to find one or more histograms: "
		  "%s (%p) %s (%p) %s (%p)", 
		  sn.Data(), sum, 
		  sn0.Data(), sum0, 
		  ev.Data(), events); 
    return false;
  }
  Sum* ret     = new Sum(GetName(), post);
  ret->fSum    = sum;
  ret->fSum0   = sum0;
  ret->fEvents = events;
  ret->fDebug  = fDebug;
  if (mc) fSumMC = ret;
  else    fSum   = ret;

  return true;
}

//____________________________________________________________________
void
AliBasedNdetaTask::CentralityBin::CreateSums(const TH2D* data, const TH2D* mc)
{
  // 
  // Create sum histogram 
  // 
  // Parameters:
  //    data  Data histogram to clone 
  //    mc    (optional) MC histogram to clone 
  //
  DGUARD(fDebug,1,"Create centrality bin sums from %s", 
	 data ? data->GetName() : "(null)");
  if (data) {
    fSum = new Sum(GetName(),"");
    fSum->Init(fSums, data, GetColor());
    fSum->fDebug = fDebug;
  }

  // If no MC data is given, then do not create MC sum histogram 
  if (!mc) return;

  fSumMC = new Sum(GetName(), "MC");
  fSumMC->Init(fSums, mc, GetColor());
  fSumMC->fDebug = fDebug;
}

//____________________________________________________________________
Bool_t
AliBasedNdetaTask::CentralityBin::CheckEvent(const AliAODForwardMult* forward,
					     Int_t     triggerMask,
					     Double_t  vzMin, 
					     Double_t  vzMax,
					     Int_t     filter)
{
  // 
  // Check the trigger, vertex, and centrality
  // 
  // Parameters:
  //    aod Event input 
  // 
  // Return:
  //    true if the event is to be used 
  //
  if (!forward) return false;

  DGUARD(fDebug,2,"Check the event "
	 "IPz: %f < %f < %f  Trigger: 0x%08x (masked 0x%08x  vetoed 0x%08x)",
	 vzMin, forward->GetIpZ(), vzMax,
	 forward->GetTriggerBits(),
	 forward->GetTriggerBits() & triggerMask,
	 forward->GetTriggerBits() & filter);
  // We do not check for centrality here - it's already done 
  Bool_t ret = forward->CheckEvent(triggerMask, vzMin, vzMax, 0, 0, 
				   fTriggers, fStatus, filter);
  DMSG(fDebug, 2, "%s", (ret ? "Accepted" : "Rejected"));
  return ret;
}
  

//____________________________________________________________________
Bool_t
AliBasedNdetaTask::CentralityBin::ProcessEvent(const AliAODForwardMult* forward,
					       UInt_t      triggerMask,
					       Bool_t      isZero,
					       Double_t    vzMin, 
					       Double_t    vzMax,
					       const TH2D* data, 
					       const TH2D* mc,
					       UInt_t      filter,
					       Double_t    weight)
{
  // 
  // Process an event
  // 
  // Parameters:
  //    forward     Forward data (for trigger, vertex, & centrality)
  //    triggerMask Trigger mask 
  //    vzMin       Minimum IP z coordinate
  //    vzMax       Maximum IP z coordinate
  //    data        Data histogram 
  //    mc          MC histogram
  //    weight      Event weight 
  //
  DGUARD(fDebug,1,"Process one event for %s a given centrality bin " 
	 "[%5.1f%%,%5.1f%%) w/weight %f",
	 data ? data->GetName() : "(null)", fLow, fHigh, weight);
  if (!CheckEvent(forward, triggerMask, vzMin, vzMax, filter)) 
    return false;
  if (!data) return false;
  if (!fSum) CreateSums(data, mc);

  fSum->Add(data, isZero, weight);
  if (mc) fSumMC->Add(mc, isZero, weight);

  return true;
}

//________________________________________________________________________
Double_t 
AliBasedNdetaTask::CentralityBin::Normalization(const TH1I& t,
						UShort_t    scheme,
						Double_t    trigEff,
						Double_t&   ntotal,
						TString*    text) const
{
  // 
  // Calculate normalization 
  // 
  // Parameters: 
  //    t       Trigger histogram
  //    scheme  Normaliztion scheme 
  //    trigEff From MC
  //    ntotal  On return, contains the number of events. 
  //
  DGUARD(fDebug,1,"Normalize centrality bin %s [%6.2f-%6.2f%%] with %s", 
	 GetName(), fLow, fHigh, t.GetName());
  Double_t nAll        = t.GetBinContent(AliAODForwardMult::kBinAll);
  Double_t nB          = t.GetBinContent(AliAODForwardMult::kBinB);
  Double_t nA          = t.GetBinContent(AliAODForwardMult::kBinA);
  Double_t nC          = t.GetBinContent(AliAODForwardMult::kBinC);
  Double_t nE          = t.GetBinContent(AliAODForwardMult::kBinE);
  Double_t nOffline    = t.GetBinContent(AliAODForwardMult::kBinOffline);
  Double_t nTriggered  = t.GetBinContent(AliAODForwardMult::kWithTrigger);
  Double_t nWithVertex = t.GetBinContent(AliAODForwardMult::kWithVertex);
  Double_t nAccepted   = ntotal; 
  ntotal               = 0;
  
  if (nTriggered <= 0.1) { 
    AliError("Number of triggered events <= 0");
    return -1;
  }
  if (nWithVertex <= 0.1) { 
    AliError("Number of events with vertex <= 0");
    return -1;
  }
  ntotal          = nAccepted;
  Double_t vtxEff = nWithVertex / nTriggered;
  Double_t scaler = 1;
  Double_t beta   = nA + nC - 2*nE;


  TString rhs("N = N_acc");
  if (!(scheme & kZeroBin)) {
    if (scheme & kEventLevel) {
      ntotal = nAccepted / vtxEff;
      scaler = vtxEff;
      DMSG(fDebug,0,"Calculating event normalisation as\n"
	   " N = N_A * N_T / N_V = %d * %d / %d = %f (%f)",
	   Int_t(nAccepted), Int_t(nTriggered), Int_t(nWithVertex), 
	   ntotal, scaler);    
      if (scheme & kBackground) {
	//          1            E_V             E_V
	//   s = --------- = ------------- = ------------ 
	//        1 - beta   1 - beta E_V    1 - beta N_V 
	//       ---  ----       --------        ---- ---
	//       E_V  N_V          N_V           N_V  N_T 
	// 
	//          E_V
	//     = ------------
	//        1 - beta 
	//            ----
	//             N_T 
	// 
	ntotal -= nAccepted * beta / nWithVertex;
	// This one is direct and correct. 
	// scaler = 1. / (1. / vtxEff - beta / nWithVertex);
	// A simpler expresion
	scaler /= (1 - beta / nTriggered); // 0.831631 -> 0.780689
	DMSG(fDebug,0,"Calculating event normalisation as\n"
	     " beta = N_a + N_c + 2 N_e = %d + %d - 2 * %d = %d\n"
	     " N = N - N_A * beta / N_V = %f - %d * %d / %d = %f (%f)",
	     Int_t(nA), Int_t(nC), Int_t(nE), Int_t(beta),
	     nAccepted / vtxEff, Int_t(nAccepted), Int_t(beta), 
	     Int_t(nWithVertex), ntotal, scaler);
	rhs.Append("(1/eps_V - beta/N_vtx)");
      } // Background 
      else 
	rhs.Append("/eps_V");
    } // Event-level
    if (scheme & kTriggerEfficiency) {
      Int_t old =  Int_t(ntotal);
      ntotal    /= trigEff;
      scaler    *= trigEff;
      DMSG(fDebug,0,"Correcting for trigger efficiency:\n"
	   " N = 1 / E_X * N = 1 / %f * %d = %f (%f)", 
	   trigEff, old, ntotal, scaler);
      rhs.Append("/eps_T");
    } // Trigger efficiency
  } 
  else  {
    // Calculate as 
    // 
    //  N = N_A + 1/E_X * N_A / N_V (N_T - N_V - beta)
    //    = N_A (1 + 1/E_X (N_T/N_V - 1 - beta / N_V))
    //    = N_A (1 + 1/E_X (1/E_V - 1 - beta / N_V))
    // 
    //  s = N_A/N = 1 / (1 + 1/E_X (N_T/N_V - 1 - beta / N_V))
    //    = N_V / (N_V + 1/E_X (N_T - N_V - beta)) 
    // 
    if (!(scheme & kBackground)) beta = 0;
    ntotal = nAccepted * (1 + 1/trigEff * (nTriggered / nWithVertex - 1 
					 - beta / nWithVertex));
    scaler = nWithVertex / (nWithVertex + 
			    1/trigEff * (nTriggered-nWithVertex-beta));
    DMSG(fDebug,0,"Calculating event normalisation as\n"
	 "  beta = N_a + N_c + 2 N_e = %d + %d - 2 * %d = %d\n"
	 "  N = N_A (1 + 1/E_X (N_T/N_V - 1 - beta / N_V)) = "
	 "%d (1 + 1 / %f (%d / %d - 1 - %d / %d)) = %f (%f)",
	 Int_t(nA), Int_t(nC), Int_t(nE), Int_t(beta),
	 Int_t(nAccepted), trigEff, Int_t(nTriggered), 
	 Int_t(nWithVertex), Int_t(beta), Int_t(nWithVertex), 
	 ntotal, scaler);
    rhs.Append("(1+1/eps_T(1/eps_V-1-beta/N_vtx))");
  }

  if (text) {
    text->Append(Form("%-40s = %d\n", "N_all",	      UInt_t(nAll)));
    text->Append(Form("%-40s = %d\n", "N_acc",	      UInt_t(nAccepted)));
    text->Append(Form("%-40s = %d\n", "N_trg",        UInt_t(nTriggered)));
    text->Append(Form("%-40s = %d\n", "N_vtx",	      UInt_t(nWithVertex)));
    text->Append(Form("%-40s = %d\n", "N_B",	      UInt_t(nB)));
    text->Append(Form("%-40s = %d\n", "N_A",	      UInt_t(nA)));
    text->Append(Form("%-40s = %d\n", "N_C",  	      UInt_t(nC)));
    text->Append(Form("%-40s = %d\n", "N_E",	      UInt_t(nE)));
    text->Append(Form("%-40s = %d\n", "beta = N_A + N_C - 2N_E",UInt_t(beta)));
    text->Append(Form("%-40s = %f\n", "eps_V = N_vtx/N_trg",vtxEff));
    text->Append(Form("%-40s = %f\n", "eps_T",	      trigEff));
    text->Append(Form("%-40s = %f\n", rhs.Data(),     ntotal));
  }
  TString tN = t.GetXaxis()->GetBinLabel(AliAODForwardMult::kWithTrigger);
  tN.ReplaceAll("w/Selected trigger (","");
  tN.ReplaceAll(")", "");
  DMSG(fDebug,0,"\n"
       " Total of        %9d events for %s\n"
       "  of these       %9d have an offline trigger\n"
       "  of these N_T = %9d has the selected trigger (%s)\n"
       "  of these N_V = %9d has a vertex\n"
       "  of these N_A = %9d were in the selected range",
       Int_t(nAll), GetTitle(),
       Int_t(nOffline),
       Int_t(nTriggered), tN.Data(),
       Int_t(nWithVertex),
       Int_t(nAccepted));       
  DMSG(fDebug,0,"\n"
       "  Triggers by hardware type:\n"
       "    N_b   =  %9d\n"
       "    N_ac  =  %9d (%d+%d)\n"
       "    N_e   =  %9d\n"
       "  Vertex efficiency:          %f\n"
       "  Trigger efficiency:         %f\n"
       "  Total number of events: N = %f\n"
       "  Scaler (N_A/N):             %f\n"
       "  %-25s = %f",
       Int_t(nB),
       Int_t(nA+nC),Int_t(nA),Int_t(nC),
       Int_t(nE),
       vtxEff,
       trigEff,
       ntotal,
       scaler,
       rhs.Data(), ntotal);
  return scaler;
}

//________________________________________________________________________
const char* 
AliBasedNdetaTask::CentralityBin::GetResultName(const char* postfix) const
{
  static TString n;
  n = GetName();
  n.ReplaceAll("dNdeta", "");
  n.Prepend("dndeta");
  n.Append(postfix);
  return n.Data();
}
//________________________________________________________________________
TH1* 
AliBasedNdetaTask::CentralityBin::GetResult(const char* postfix,
					    Bool_t      verbose) const
{
  if (!fOutput) { 
    AliWarningF("No output list defined in %s [%6.2f,%6.2f]", GetName(), 
		fLow, fHigh);
    return 0;
  }
  TString  n = GetResultName(postfix);
  TObject* o = fOutput->FindObject(n.Data());
  if (!o) { 
    if (verbose)
      AliWarningF("Object %s not found in output list of %s", 
		  n.Data(), GetName());
    return 0;
  }
  return static_cast<TH1*>(o);
}

//________________________________________________________________________
void 
AliBasedNdetaTask::CentralityBin::MakeResult(const TH2D* sum,  
					     const char* postfix, 
					     bool        rootProj, 
					     bool        corrEmpty,
					     Double_t    scaler,
					     Int_t       marker,
					     Int_t       color,
					     TList*      mclist, 
					     TList*      truthlist)
{
  // 
  // Generate the dN/deta result from input 
  // 
  // Parameters: 
  //     sum        Sum of 2D hists 
  //     postfix    Post fix on names
  //     rootProj   Whether to use ROOT TH2::ProjectionX
  //     corrEmpty  Correct for empty bins 
  //     scaler     Event-level normalization scaler  
  //
  DGUARD(fDebug,1,"Make centrality bin result from %s", sum->GetName());
  TString base(GetName());
  base.ReplaceAll("dNdeta", "");
  base.Append(postfix);
  TH2D* copy    = static_cast<TH2D*>(sum->Clone(Form("d2Ndetadphi%s",
						     base.Data())));
  
  TH1D* accNorm = 0;
  Int_t nY      = sum->GetNbinsY();
  // Hack HHD Hans test code to manually remove FMD2 dead channel (but
  // it is on outer?)
  // 
  // cholm comment: The original hack has been moved to
  // AliForwarddNdetaTask::CheckEventData - this simplifies things a
  // great deal, and we could in principle use the new phi acceptance.
  // 
  // However, since we may have filtered out the dead sectors in the
  // AOD production already, we can't be sure we can recalculate the
  // phi acceptance correctly, so for now, we rely on fCorrEmpty being set. 
  Int_t o       = (corrEmpty ? 0 : nY+1);
  accNorm = ProjectX(sum, Form("norm%s",base.Data()), o, o, 
		     rootProj, corrEmpty, false);
  accNorm->SetDirectory(0);

  // --- Normalize to the coverage -----------------------------------
  if (corrEmpty) {
    ScaleToCoverage(copy, accNorm);
    // --- Event-level normalization ---------------------------------
    copy->Scale(scaler);
  }

  // --- Project on X axis -------------------------------------------
  TH1D* dndeta = ProjectX(copy, Form("dndeta%s",base.Data()),
			  1, nY, rootProj, corrEmpty);
  dndeta->SetDirectory(0);
  // Event-level normalization 
  if (!corrEmpty) {
    ScaleToCoverage(dndeta, accNorm);
    dndeta->Scale(scaler);
  }
  dndeta->Scale(1., "width");
  copy->Scale(1., "width");
  
  TH1D*  dndetaMCCorrection = 0;
  TH1D*  dndetaMCtruth      = 0;
  TList* centlist           = 0;
  TList* truthcentlist      = 0;
  
  // --- Possible final correction to <MC analysis> / <MC truth> -----
  // we get the rebinned distribution for satellite to make the correction
  TString rebinSuf(fSatelliteVertices ? "_rebin05" : "");
  if(mclist) {
    centlist = static_cast<TList*> (mclist->FindObject(GetListName()));
    if(centlist)
      dndetaMCCorrection = 
	static_cast<TH1D*>(centlist->FindObject(Form("dndeta%s%s",
						     base.Data(),
						     rebinSuf.Data())));
  }
  if (truthlist) {
    truthcentlist = static_cast<TList*>(truthlist->FindObject(GetListName()));
    if (truthcentlist)
      // TODO here new is "dndetaTruth"
      dndetaMCtruth = 
	static_cast<TH1D*>(truthcentlist->FindObject(Form("dndetaMCTruth%s",
							  rebinSuf.Data())));
  }
  
  if (dndetaMCCorrection && dndetaMCtruth) {
    AliInfo("Correcting with final MC correction");
    TString testString(dndetaMCCorrection->GetName());

    // We take the measured MC dN/deta and divide with truth 
    dndetaMCCorrection->Divide(dndetaMCtruth);
    dndetaMCCorrection->SetTitle("Final MC correction");
    dndetaMCCorrection->SetName("finalMCCorr");
    for(Int_t m = 1; m <= dndetaMCCorrection->GetNbinsX(); m++) {
      if(dndetaMCCorrection->GetBinContent(m) < 0.5 || 
	 dndetaMCCorrection->GetBinContent(m) > 1.75) {
	dndetaMCCorrection->SetBinContent(m,1.);
	dndetaMCCorrection->SetBinError(m,0.1);
      }
    }
    // Applying the correction
    if (!fSatelliteVertices)
      // For non-satellites we took the same binning, so we do a straight 
      // division 
      dndeta->Divide(dndetaMCCorrection);
    else {
      // For satelitte events, we took coarser binned histograms, so 
      // we need to do a bit more 
      for(Int_t m = 1; m <= dndeta->GetNbinsX(); m++) {
	if(dndeta->GetBinContent(m) <= 0.01 ) continue;
	
	Double_t eta     = dndeta->GetXaxis()->GetBinCenter(m);
	Int_t    bin     = dndetaMCCorrection->GetXaxis()->FindBin(eta);
	Double_t mccorr  = dndetaMCCorrection->GetBinContent(bin);
	Double_t mcerror = dndetaMCCorrection->GetBinError(bin);
	if (mccorr < 1e-6) {
	  dndeta->SetBinContent(m, 0);
	  dndeta->SetBinError(m, 0);
	}
	Double_t value   = dndeta->GetBinContent(m);
	Double_t error   = dndeta->GetBinError(m);
	Double_t sumw2   = (error   * error   * mccorr * mccorr +
			    mcerror * mcerror * value  * value);
	dndeta->SetBinContent(m,value/mccorr) ;
	dndeta->SetBinError(m,TMath::Sqrt(sumw2)/mccorr/mccorr);
      }
    }
  }
  else 
    DMSG(fDebug,1,"No final MC correction applied");
  
  // --- Set some histogram attributes -------------------------------
  TString post;
  Int_t rColor = GetColor(color);
  if (postfix && postfix[0] != '\0') post = Form(" (%s)", postfix);
  SetHistogramAttributes(dndeta,  rColor, marker, 
			 Form("ALICE %s%s", GetName(), post.Data()));
  SetHistogramAttributes(accNorm, rColor, marker, 
			 Form("ALICE %s normalisation%s", 
			      GetName(), post.Data())); 

  // --- Add the histograms to output list ---------------------------
  fOutput->Add(dndeta);
  fOutput->Add(accNorm);
  fOutput->Add(copy);
  if (dndetaMCCorrection) fOutput->Add(dndetaMCCorrection);
  
  // HHD Test of dN/deta in phi bins add flag later?
  // 
  // cholm comment: We disable this for now 
#if 0
  for (Int_t nn=1; nn <= sum->GetNbinsY(); nn++) {
    TH1D* dndeta_phi = ProjectX(copy, Form("dndeta%s_phibin%d",
					   base.Data(), nn), 
				nn, nn, rootProj, corrEmpty);
    dndeta_phi->SetDirectory(0);
    // Event-level normalization 
    dndeta_phi->Scale(TMath::Pi()/10., "width");
     
    if(centlist)
      dndetaMCCorrection = 
	static_cast<TH1D*>(centlist->FindObject(Form("dndeta%s_phibin%d",
						     base.Data(),nn)));
    if(truthcentlist)
      dndetaMCtruth 
	= static_cast<TH1D*>(truthcentlist->FindObject("dndetaMCTruth"));

    if (dndetaMCCorrection && dndetaMCtruth) {
      AliInfo("Correcting with final MC correction");
      TString testString(dndetaMCCorrection->GetName());
      dndetaMCCorrection->Divide(dndetaMCtruth);
      dndetaMCCorrection->SetTitle(Form("Final_MC_correction_phibin%d",nn));
      dndetaMCCorrection->SetName(Form("Final_MC_correction_phibin%d",nn));
      for(Int_t m = 1; m <= dndetaMCCorrection->GetNbinsX(); m++) {
	if(dndetaMCCorrection->GetBinContent(m) < 0.25 || 
	   dndetaMCCorrection->GetBinContent(m) > 1.75) {
	  dndetaMCCorrection->SetBinContent(m,1.);
	  dndetaMCCorrection->SetBinError(m,0.1);
	}
      }
      //Applying the correction
      dndeta_phi->Divide(dndetaMCCorrection);
    }
    fOutput->Add(dndeta_phi);
    if(dndetaMCCorrection) fOutput->Add(dndetaMCCorrection);
  } // End of phi
#endif
}  

//________________________________________________________________________
bool
AliBasedNdetaTask::CentralityBin::End(TList*      sums, 
				      TList*      results,
				      UShort_t    scheme,
				      Double_t    trigEff,
				      Double_t    trigEff0,
				      Bool_t      rootProj,
				      Bool_t      corrEmpty, 
				      Int_t       triggerMask,
				      Int_t       marker,
				      Int_t       color, 
				      TList*      mclist,
				      TList*      truthlist) 
{
  // 
  // End of processing 
  // 
  // Parameters:
  //    sums        List of sums
  //    results     Output list of results
  //    trigEff     Trigger efficiency 
  //    corrEmpty   Whether to correct for empty bins
  //    triggerMask Trigger mask 
  //
  DGUARD(fDebug,1,"End centrality bin procesing");

  fSums = dynamic_cast<TList*>(sums->FindObject(GetListName()));
  if(!fSums) {
    AliErrorF("Could not retrieve TList fSums: %s", GetListName());
    sums->ls();
    return false; 
  }
  
  fOutput = new TList;
  fOutput->SetName(GetListName());
  fOutput->SetOwner();
  results->Add(fOutput);

  if (!fSum) { 
    if (!ReadSum(fSums, false)) {
      AliInfo("This task did not produce any output");
      return false;
    }
  }
  if (!fSumMC) ReadSum(fSums, true);

  fTriggers = static_cast<TH1I*>(fSums->FindObject("triggers"));

  if (!fTriggers) { 
    AliError("Couldn't find histogram 'triggers' in list");
    return false;
  }

  // --- Get normalization scaler ------------------------------------
  Double_t epsilonT  = trigEff;
  Double_t epsilonT0 = trigEff0;
  DMSG(fDebug,2,"Using epsilonT=%f, epsilonT0=%f for 0x%x", 
       epsilonT, epsilonT0, triggerMask);

  // Get our histograms 
  Double_t nSum   = 0;
  TH2D*    sum    = fSum->CalcSum(fOutput, nSum, epsilonT0, epsilonT, 
				  marker, rootProj, corrEmpty);
  Double_t nSumMC = 0;
  TH2D*    sumMC  = 0;
  if (fSumMC) sumMC = fSumMC->CalcSum(fOutput, nSumMC, 
				      epsilonT0, epsilonT, marker,
				      rootProj, corrEmpty);
  if (!sum) { 
    AliError("Failed to get sum from summer - bailing out");
    return false;
  }
    
  TString  text;
  Double_t ntotal = nSum;
  Double_t scaler = Normalization(*fTriggers, scheme, epsilonT, ntotal, &text);
  if (scaler < 0) { 
    AliError("Failed to calculate normalization - bailing out");
    return false;
  }
  fOutput->Add(fTriggers->Clone());
  fOutput->Add(new TNamed("normCalc", text.Data()));
  fOutput->Add(new TParameter<float>("low", fLow,  'f'));
  fOutput->Add(new TParameter<float>("high",fHigh, 'f'));

  // --- Make result and store ---------------------------------------
  MakeResult(sum, "", rootProj, corrEmpty, scaler, marker, color, 
	     mclist, truthlist);

  // --- Process result from TrackRefs -------------------------------
  if (sumMC) 
    MakeResult(sumMC, "MC", rootProj, corrEmpty, scaler,
	       GetMarkerStyle(GetMarkerBits(marker)+4), color, 
	       mclist, truthlist);
  
  // Temporary stuff 
  // if (!IsAllBin()) return;
  return true;
}

//____________________________________________________________________
void
AliBasedNdetaTask::CentralityBin::Print(Option_t* option) const
{
  PFV("Centrality bin", GetTitle());
  gROOT->IncreaseDirLevel();
  PF("Range", "%6.2f - %6.2f%%", fLow, fHigh);
  PFB("All bin", IsAllBin());
  PFB("Final MC", fDoFinalMCCorrection);
  PFB("Satellite collisions", fSatelliteVertices);

  TString opt(option);
  opt.ToUpper();
  if (opt.Contains("R")) {
    if (fSum)   fSum->Print(option);
    if (fSumMC) fSumMC->Print(option);
  }

  gROOT->DecreaseDirLevel();
}

//====================================================================
Bool_t 
AliBasedNdetaTask::ApplyEmpiricalCorrection(const AliAODForwardMult* aod,
					    TH2D* data)
{
  if (!fEmpiricalCorrection || !data)
    return true;
  Float_t zvertex=aod->GetIpZ();
  Int_t binzvertex=fEmpiricalCorrection->GetXaxis()->FindBin(zvertex);
  if(binzvertex<1||binzvertex>fEmpiricalCorrection->GetNbinsX())
    return false;
  for (int i=1;i<=data->GetNbinsX();i++) {
    Int_t bincorrection=fEmpiricalCorrection->GetYaxis()
      ->FindBin(data->GetXaxis()->GetBinCenter(i));
    if(bincorrection<1||bincorrection>fEmpiricalCorrection->GetNbinsY())
      return false;
    Float_t correction=fEmpiricalCorrection
      ->GetBinContent(binzvertex,bincorrection);
    if(correction<0.001) {
      data->SetBinContent(i,0,0);
      data->SetBinContent(i,data->GetNbinsY()+1,0);
    }	
    for(int j=1;j<=data->GetNbinsY();j++) {
      if (data->GetBinContent(i,j)>0.0) {
	data->SetBinContent(i,j,data->GetBinContent(i,j)*correction);
	data->SetBinError(i,j,data->GetBinError(i,j)*correction);
      }	
    }
  }
  return true;
}

//____________________________________________________________________
TH1*
AliBasedNdetaTask::Asymmetry(TH1* h)
{
  if (!h) return 0;
  
  TH1* ret = static_cast<TH1*>(h->Clone(Form("%s_leftright", h->GetName())));
  // Int_t    oBins = h->GetNbinsX();
  // Double_t high  = h->GetXaxis()->GetXmax();
  // Double_t low   = h->GetXaxis()->GetXmin();
  // Double_t dBin  = (high - low) / oBins;
  // Int_t    tBins = Int_t(2*high/dBin+.5);
  // ret->SetBins(tBins, -high, high);
  ret->SetDirectory(0);
  ret->Reset();
  ret->SetTitle(Form("%s (+/-)", h->GetTitle()));
  ret->SetYTitle("Right/Left");
  Int_t nBins = h->GetNbinsX();
  for (Int_t i = 1; i <= nBins; i++) { 
    Double_t x   = h->GetBinCenter(i);
    if (x > 0) break;
    
    Double_t c1  = h->GetBinContent(i);
    Double_t e1  = h->GetBinError(i);
    if (c1 <= 0) continue; 
    
    Int_t    j   = h->FindBin(-x);
    if (j <= 0 || j > nBins) continue;
    
    Double_t c2  = h->GetBinContent(j);
    Double_t e2  = h->GetBinError(j);    
    Double_t c12 = c1*c1;
    Double_t e   = TMath::Sqrt((e2*e2*c1*c1+e1*e1*c2*c2)/(c12*c12));    
    Int_t    k   = ret->FindBin(x);
    ret->SetBinContent(k, c2/c1);
    ret->SetBinError(k, e);
  }
  TF1* fit = new TF1("fit", "pol0");
  ret->Fit(fit,"QN+");
  // TF1* fit = ret->GetFunction("pol0");
  if (fit) {
    TLatex* ltx = new TLatex(-1, fit->GetParameter(0),
			     Form("%6.4f#pm%6.4f #chi^{2}/#nu=%5.2f",
				  fit->GetParameter(0),
				  fit->GetParError(0),
				  fit->GetChisquare()/fit->GetNDF()));
    ltx->SetTextColor(ret->GetMarkerColor());
    ltx->SetTextAlign(12);
    if (!ret->GetListOfFunctions()->FindObject(fit))
      ret->GetListOfFunctions()->Add(fit);
    ret->GetListOfFunctions()->Add(ltx);
  }
  
  return ret;
}

//
// EOF
//
