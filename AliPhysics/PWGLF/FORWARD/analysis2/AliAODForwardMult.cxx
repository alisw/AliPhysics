//
// Class that contains the forward multiplicity data per event 
//
// This class contains a histogram of 
// @f[
//   \frac{d^2N_{ch}}{d\eta d\phi}\quad,
// @f]
// as well as a trigger mask for each analysed event.  
// 
// The eta acceptance of the event is stored in the underflow bins of
// the histogram.  So to build the final histogram, one needs to
// correct for this acceptance (properly weighted by the events), and
// the vertex efficiency.  This simply boils down to defining a 2D
// histogram and summing the event histograms in that histogram.  One
// should of course also do proper book-keeping of the accepted event.
//
#include "AliAODForwardMult.h"
#include <TBrowser.h>
#include <iostream>
#include <TMath.h>
#include <TObjString.h>
#include <TObjArray.h>
#include "AliLog.h"
ClassImp(AliAODForwardMult)
#if 0 // def DOXY_INPUT
; // For Emacs 
#endif

//____________________________________________________________________
const Float_t AliAODForwardMult::fgkInvalidIpZ = 1e6;

//____________________________________________________________________
AliAODForwardMult::AliAODForwardMult()
  : fIsMC(false),
    fHist(),
    fTriggers(0),
    fIpZ(fgkInvalidIpZ), 
    fCentrality(-1),				
    fNClusters(0)
{
  // 
  // Constructor 
  // 
}

//____________________________________________________________________
AliAODForwardMult::AliAODForwardMult(Bool_t isMC) 
  : fIsMC(isMC),
    fHist("forwardMult", "d^{2}N_{ch}/d#etad#varphi in the forward regions", 
	  200, -4, 6, 20, 0, 2*TMath::Pi()),
    fTriggers(0),
    fIpZ(fgkInvalidIpZ), 
    fCentrality(-1),				
    fNClusters(0)
{
  // 
  // Constructor 
  // 
  // Parameters: 
  //  isMC   If set to true this is for MC data (effects branch name)
  // 
  fHist.SetXTitle("#eta");
  fHist.SetYTitle("#varphi [radians]");
  fHist.SetZTitle("#frac{d^{2}N_{ch}}{d#etad#varphi}");
  fHist.SetDirectory(0);
  fHist.Sumw2();
}

//____________________________________________________________________
void
AliAODForwardMult::Init(const TAxis& etaAxis)
{
  // Initialize the histogram with an eta axis 
  // 
  // Parameters: 
  //   etaAxis       Eta axis to use 
  // 
  fHist.SetBins(etaAxis.GetNbins(), etaAxis.GetXmin(), etaAxis.GetXmax(), 
		20, 0, 2*TMath::Pi());
}

//____________________________________________________________________
void
AliAODForwardMult::Clear(Option_t* option)
{
  // Clear (or reset) internal values 
  // 
  // Parameters: 
  //  option   Passed to TH1::Reset 
  // 
  fHist.Reset(option);
  fTriggers  = 0;
  fIpZ       = fgkInvalidIpZ;
  fNClusters = 0;
}

//____________________________________________________________________
TH1*
AliAODForwardMult::GetEtaCoverage() const
{
  Int_t bin = GetEtaCoverageBin();
  TH1*  ret = fHist.ProjectionX(Form("%s_etacov", fHist.GetName()),
				bin, bin, "e");
  ret->SetDirectory(0);
  ret->SetYTitle("#it{#eta} coverage");
  return ret;
}
//____________________________________________________________________
TH1*
AliAODForwardMult::GetPhiAcceptance() const
{
  Int_t bin = GetPhiAcceptanceBin();
  TH1*  ret = fHist.ProjectionX(Form("%s_phiacc", fHist.GetName()),
				bin, bin, "e");
  ret->SetDirectory(0);
  ret->SetYTitle("#it{#varphi} acceptance");
  return ret;
}
//____________________________________________________________________
void
AliAODForwardMult::FillEtaCoverage(TH1& h) const
{
  TH1* hh = GetEtaCoverage();
  h.Add(hh);
}
//____________________________________________________________________
void
AliAODForwardMult::FillPhiAcceptance(TH1& h) const
{
  TH1* hh = GetPhiAcceptance();
  h.Add(hh);
}
//____________________________________________________________________
void
AliAODForwardMult::SetSNN(UShort_t snn)
{
  // set the center of mass energy per nucleon pair (GeV). 
  // This is stored in bin (0,0) of the histogram 
  // 
  // Parameters: 
  //   sNN   Center of mass energy per nuclean 
  fHist.SetBinContent(0,0,snn);
}
//____________________________________________________________________
void
AliAODForwardMult::SetSystem(UShort_t sys)
{
  // set the center of mass energy per nucleon pair (GeV). 
  // This is stored in bin (N+1,0) of the histogram 
  // 
  // Parameters: 
  //   sys   Collision system number 
  fHist.SetBinContent(fHist.GetNbinsX()+1,0,sys);
}

//____________________________________________________________________
Bool_t
AliAODForwardMult::HasIpZ() const
{
  // Check if we have valid z coordinate of the interaction point 
  // 
  // Return:
  //   true if the z coordinate of the interaction point is valid 
  // 
  return TMath::Abs(fIpZ - fgkInvalidIpZ) > 1;
}

//____________________________________________________________________
UShort_t
AliAODForwardMult::GetSNN() const
{
  // set the center of mass energy per nucleon pair (GeV). 
  // This is stored in bin (0,0) of the histogram 
  // 
  // Parameters: 
  //   sNN   Center of mass energy per nuclean 
  return UShort_t(fHist.GetBinContent(0,0));
}

//____________________________________________________________________
UShort_t
AliAODForwardMult::GetSystem() const
{
  // set the center of mass energy per nucleon pair (GeV). 
  // This is stored in bin (N+1,0) of the histogram 
  // 
  // Parameters: 
  //   sNN   Center of mass energy per nuclean 
  return UShort_t(fHist.GetBinContent(fHist.GetNbinsX()+1,0));
}

//____________________________________________________________________
void
AliAODForwardMult::Browse(TBrowser* b)
{
  // Browse this object 
  // 
  // Parameters: 
  //   b   Browser to use 
  static TObjString ipz;
  static TObjString trg;
  static TObjString cnt;
  static TObjString ncl;
  ipz = Form("ip_z=%fcm", fIpZ);
  trg = GetTriggerString(fTriggers);
  cnt = Form("%+6.1f%%", fCentrality);
  ncl = Form("%d clusters", fNClusters);
  b->Add(&fHist);
  b->Add(&ipz);
  b->Add(&trg);
  b->Add(&cnt);
  b->Add(&ncl);
}

namespace { 
  void AppendAnd(TString& trg, const char* sep, const TString& what)
  {
    if (!trg.IsNull()) trg.Append(Form(" %s ", sep));
    trg.Append(what);
  }
}
//____________________________________________________________________
const Char_t*
AliAODForwardMult::GetTriggerString(UInt_t mask, const char* sep)
{
  // Get a string that describes the triggers 
  // 
  // Parameters: 
  //   mask  Bit pattern of triggers 
  // Return: 
  //   Character string representation of mask 
  static TString trg;
  trg = "";
  if (mask == 0) return "none";

  bool   isOr = mask & kInclusive;
  UInt_t tmp  = 0x7FFFFFFF & mask;
  TString s(sep);
  if (s.IsNull()) s = isOr ? "|" : "&";
  
  if ((tmp & kInel)        != 0x0) AppendAnd(trg, s, "MBOR");
  if ((tmp & kInelGt0)     != 0x0) AppendAnd(trg, s, "INEL>0");
  if ((tmp & kNSD)         != 0x0) AppendAnd(trg, s, "MBAND");
  if ((tmp & kV0AND)       != 0x0) AppendAnd(trg, s, "V0AND");
  if ((tmp & kA)           != 0x0) AppendAnd(trg, s, "A");
  if ((tmp & kB)           != 0x0) AppendAnd(trg, s, "B");
  if ((tmp & kC)           != 0x0) AppendAnd(trg, s, "C");
  if ((tmp & kE)           != 0x0) AppendAnd(trg, s, "E");
  if ((tmp & kMCNSD)       != 0x0) AppendAnd(trg, s, "MCNSD");
  if ((tmp & kNClusterGt0) != 0x0) AppendAnd(trg, s, "NCluster>0");
  if ((tmp & kSatellite)   != 0x0) AppendAnd(trg, s, "Satellite");
  if ((tmp & kOffline)     != 0x0) AppendAnd(trg, s, "Offline");
  if ((tmp & kSPDOutlier)  != 0x0) AppendAnd(trg, s, "Outlier");
  if ((tmp & kPileUp)      != 0x0) AppendAnd(trg, s, "Pileup");
  if ((tmp & kPileupSPD)   != 0x0) AppendAnd(trg, s, "Pileup-SPD");
  if ((tmp & kPileupTrack) != 0x0) AppendAnd(trg, s, "Pileup-TRK");
  if ((tmp & kPileupBC)    != 0x0) AppendAnd(trg, s, "Pileup-BC");
  if ((tmp & kPileupBins)  != 0x0) AppendAnd(trg, s, "Pileup-BIN");
  if ((tmp & kADOR)        != 0x0) AppendAnd(trg, s, "ADOR");
  if ((tmp & kADAND)       != 0x0) AppendAnd(trg, s, "ADAND");
  return trg.Data();
}
  
//____________________________________________________________________
TH1I*
AliAODForwardMult::MakeTriggerHistogram(const char* name, UInt_t triggerMask) 
{
  // 
  // Make a histogram to record triggers in. 
  //
  // The bins defined by the trigger enumeration in this class.  One
  // can use this enumeration to retrieve the number of triggers for
  // each class.
  // 
  // Parameters:
  //    name Name of the histogram 
  // 
  // Return:
  //    Newly allocated histogram 
  //
  TString sel = GetTriggerString(triggerMask, "");
  TH1I*   ret = new TH1I(name,
			 Form("Triggers (%s)", sel.Data()),
			 kAccepted, .5, kAccepted+.5);
  ret->SetYTitle("Events");
  ret->SetFillColor(kRed-2);
  ret->SetFillStyle(3002);
  ret->GetXaxis()->SetBinLabel(kBinAll,         "All events");
  ret->GetXaxis()->SetBinLabel(kBinInel,        "Coll. & MBOR");
  ret->GetXaxis()->SetBinLabel(kBinInelGt0,     "Coll. & MBOR&&nTracklet>0");
  ret->GetXaxis()->SetBinLabel(kBinNSD,         "Coll. & V0AND||FASTOR>5");
  ret->GetXaxis()->SetBinLabel(kBinV0AND,       "Coll. & V0AND");
  ret->GetXaxis()->SetBinLabel(kBinADOR,        "Coll. & ADOR");
  ret->GetXaxis()->SetBinLabel(kBinADAND,       "Coll. & ADAND");
  ret->GetXaxis()->SetBinLabel(kBinB,           "B (Coll.) & Sel.");
  ret->GetXaxis()->SetBinLabel(kBinA,           "A & Sel.");
  ret->GetXaxis()->SetBinLabel(kBinC,           "C & Sel.");
  ret->GetXaxis()->SetBinLabel(kBinE,           "E & Sel.");
  ret->GetXaxis()->SetBinLabel(kBinSatellite,   "Satellite");
  ret->GetXaxis()->SetBinLabel(kBinMCNSD,       "NSD (MC truth)");
  ret->GetXaxis()->SetBinLabel(kBinPileUp,      "w/Pileup");
  ret->GetXaxis()->SetBinLabel(kBinOffline,     "w/Offline");
  ret->GetXaxis()->SetBinLabel(kBinNClusterGt0, "w/N_{cluster}>1");
  ret->GetXaxis()->SetBinLabel(kWithVertex,     "w/Vertex");
  ret->GetXaxis()->SetBinLabel(kWithTrigger,    "w/Selected trigger");
  ret->GetXaxis()->SetBinLabel(kAccepted,       "Accepted by cut");
  ret->GetXaxis()->SetNdivisions(kAccepted, false);
  ret->SetStats(0);

  return ret;
}

//____________________________________________________________________
TH1I*
AliAODForwardMult::MakeStatusHistogram(const char* name) 
{
  // 
  // Make a histogram to record status in. 
  //
  // The bins defined by the status enumeration in this class.  
  // 
  // Parameters:
  //    name Name of the histogram 
  // 
  // Return:
  //    Newly allocated histogram 
  //
  Int_t nBins = kOutlierEvent;
  TH1I* ret = new TH1I(name, "Event selection status", nBins+1, -.5, nBins+.5);
  ret->SetYTitle("Events");
  ret->SetFillColor(kBlue+1);
  ret->SetFillStyle(3002);
  ret->GetXaxis()->SetBinLabel(kGoodEvent+1,       "Good");
  ret->GetXaxis()->SetBinLabel(kWrongCentrality+1, "Out-of-range centrality");
  ret->GetXaxis()->SetBinLabel(kWrongTrigger+1,    "Wrong trigger");
  ret->GetXaxis()->SetBinLabel(kIsPileup+1,        "Pile-up");
  ret->GetXaxis()->SetBinLabel(kIsFilterOut+1,     "Filtered out");
  ret->GetXaxis()->SetBinLabel(kNoVertex+1,        "No IP_{z}");
  ret->GetXaxis()->SetBinLabel(kWrongVertex+1,     "Out-or-range IP_{z}");
  ret->GetXaxis()->SetBinLabel(kOutlierEvent+1,    "SPD Outlier");
  ret->GetXaxis()->SetNdivisions(nBins, false);
  ret->SetStats(0);
  return ret;
}

//____________________________________________________________________
UInt_t 
AliAODForwardMult::MakeTriggerMask(const char* what, const char* sep)
{
  UInt_t      trgMask = 0;
  TString     trgs(what);
  trgs.ToUpper();
  if (trgs.EqualTo("NONE") || trgs.EqualTo("ALL")) return trgMask;
  TObjArray*  parts = trgs.Tokenize(sep);
  TObjString* trg;
  TIter       next(parts);
  if      (sep && sep[0] == '|') trgMask |= kInclusive; 
  else if (sep && sep[0] == '&') trgMask &= ~kInclusive;  
  while ((trg = static_cast<TObjString*>(next()))) { 
    TString s(trg->GetString());
    s.Strip(TString::kBoth, ' ');
    s.ToUpper();
    // Printf("Full: %s, part: %s", what, s.Data());
    if      (s.IsNull()) continue;
    if      (s.CompareTo("INEL")       == 0) trgMask |= kInel;
    else if (s.CompareTo("MBOR")       == 0) trgMask |= kInel;
    else if (s.CompareTo("INEL>0")     == 0) trgMask |= kInelGt0;
    else if (s.CompareTo("INELGT0")    == 0) trgMask |= kInelGt0;
    else if (s.CompareTo("MBAND")      == 0) trgMask |= kNSD;
    else if (s.CompareTo("NSD")        == 0) trgMask |= kV0AND;
    else if (s.CompareTo("V0AND")      == 0) trgMask |= kV0AND;
    else if (s.CompareTo("MCNSD")      == 0) trgMask |= kMCNSD;
    else if (s.CompareTo("B")          == 0) trgMask |= kB;
    else if (s.CompareTo("A")          == 0) trgMask |= kA;
    else if (s.CompareTo("C")          == 0) trgMask |= kC;
    else if (s.CompareTo("SAT")        == 0) trgMask |= kSatellite;
    else if (s.CompareTo("E")          == 0) trgMask |= kE;
    else if (s.CompareTo("NCLUSTER>0") == 0) trgMask |= kNClusterGt0;
    else if (s.CompareTo("CENT")       == 0) trgMask |= kInel;
    else if (s.CompareTo("MULT")       == 0) trgMask |= kInel;
    else if (s.CompareTo("OFFLINE")    == 0) trgMask |= kOffline;
    else if (s.CompareTo("OUTLIER")    == 0) trgMask |= kSPDOutlier;
    else if (s.CompareTo("PILEUP-SPD") == 0) trgMask |= kPileupSPD;
    else if (s.CompareTo("PILEUP-TRK") == 0) trgMask |= kPileupTrack;
    else if (s.CompareTo("PILEUP-BIN") == 0) trgMask |= kPileupBins;
    else if (s.CompareTo("PILEUP-BC")  == 0) trgMask |= kPileupBC;
    else if (s.CompareTo("PILEUP")     == 0) trgMask |= kPileUp;
    else if (s.CompareTo("ADOR")       == 0) trgMask |= kADOR;
    else if (s.CompareTo("ADAND")      == 0) trgMask |= kADAND;
    else if (s.CompareTo("CENTCALIB")  == 0) trgMask |= kCentNoCalib;
    else if (s.BeginsWith("INCL"))           trgMask |= kInclusive;
    // trgMask &= ~(kInel|kInelGt0|kNSD|kV0AND|kMCNSD);
    else 
      AliWarningGeneral("MakeTriggerMask", 
			Form("Unknown trigger %s", s.Data()));
  }
  delete parts;
  return trgMask;
}


//____________________________________________________________________
void
AliAODForwardMult::FillTriggerHistogram(UInt_t triggerMask,
					UInt_t trg,
					TH1* hist)
{
  if (!hist) return;

  UInt_t tmp = triggerMask & ~kB;
  hist->Fill(kBinAll);
  if (IsTriggerBits(kB|kInel,    trg)) hist->Fill(kBinInel);
  if (IsTriggerBits(kB|kInelGt0, trg)) hist->Fill(kBinInelGt0);
  if (IsTriggerBits(kB|kNSD,     trg)) hist->Fill(kBinNSD);
  if (IsTriggerBits(kB|kV0AND,   trg)) hist->Fill(kBinV0AND);
  if (IsTriggerBits(kB|kADOR,    trg)) hist->Fill(kBinADOR);
  if (IsTriggerBits(kB|kADAND,   trg)) hist->Fill(kBinADAND);
  if (IsTriggerBits(kB,          trg) &&
      IsTriggerBits(tmp,         trg)) hist->Fill(kBinB);
  if (IsTriggerBits(kA,          trg) &&
      IsTriggerBits(tmp,         trg)) hist->Fill(kBinA);
  if (IsTriggerBits(kC,          trg) &&
      IsTriggerBits(tmp,         trg)) hist->Fill(kBinC);
  if (IsTriggerBits(kE,          trg) &&
      IsTriggerBits(tmp,         trg)) hist->Fill(kBinE);
  if (IsTriggerBits(kSatellite,  trg)) hist->Fill(kBinSatellite);
  if (IsTriggerBits(kMCNSD,      trg)) hist->Fill(kBinMCNSD);
  if (IsTriggerBits(kPileUp,     trg)) hist->Fill(kBinPileUp);
  if (IsTriggerBits(kOffline,    trg)) hist->Fill(kBinOffline);
  if (IsTriggerBits(kNClusterGt0,trg)) hist->Fill(kBinNClusterGt0);
  if (IsTriggerBits(triggerMask,trg) &&
      !IsTriggerBits(kB|tmp,trg))
    ::Warning("FillTriggerHistogram",
	      "event: 0x%x, mask: 0x%x, tmp: 0x%x, tmp|b: 0x%x",
	      trg, triggerMask, tmp, tmp|kB);
}

//____________________________________________________________________
Bool_t
AliAODForwardMult::FilterEvent(Double_t vzMin, Double_t vzMax,
			       TH1* hist, TH1* status,
			       UInt_t filterMask) const
{
  // Check if event is SPD outlier
  if (filterMask & kSPDOutlier && IsTriggerBits(kSPDOutlier)) {
    if (status) status->Fill(kOutlierEvent);
    return false;
  }
  // Check for pileup
  UInt_t checkPileup = (filterMask & (kPileUp      |
				      kPileupSPD   |
				      kPileupTrack |
				      kPileupBC    |
				      kPileupBins));
  if (checkPileup && IsTriggerBits(kInclusive | checkPileup)) {
    if (status) status->Fill(kIsPileup);
    return false;
  }
  UInt_t other = filterMask & ~(kPileUp      |
				kPileupSPD   |
				kPileupTrack |
				kPileupBC    |
				kPileupBins  | 
				kSPDOutlier);
  // Unspecified stuff 
  if (other && IsTriggerBits(kInclusive | other)) {
    if (status) status->Fill(kIsFilterOut);
    return false;
  }
  // Got a non-vetoed trigger 
  if (hist) hist->Fill(kWithTrigger);
  
  // Check that we have a valid vertex
  if (vzMin < vzMax && !HasIpZ()) {
    if (status) status->Fill(kNoVertex);
    return false;
  }
  // Got valid vertex 
  if (hist) hist->Fill(kWithVertex);

  // Check that vertex is within cuts 
  if (vzMin < vzMax && !InRange(vzMin, vzMax)) {
    if (status) status->Fill(kWrongVertex);
    return false;
  }
  // Trigger and vertex 
  if (hist)   hist->Fill(kAccepted);
  if (status) status->Fill(kGoodEvent);

  return true;
}

//____________________________________________________________________
Bool_t
AliAODForwardMult::CheckEvent(UInt_t   triggerMask,
			      Double_t vzMin, Double_t vzMax,
			      Double_t cMin,  Double_t cMax, 
			      TH1*     hist,  TH1*     status,
			      UInt_t   filterMask) const 
{
  // 
  // Check if event meets the passses requirements.   
  //
  // It returns true if @e all of the following is true 
  //
  // - The trigger is within the bit mask passed.
  // - The vertex is within the specified limits. 
  // - The centrality is within the specified limits, or if lower
  //   limit is equal to or larger than the upper limit.
  // 
  // If a histogram is passed in the last parameter, then that
  // histogram is filled with the trigger bits. 
  // 
  // Parameters:
  //    triggerMask  Trigger mask
  //    vzMin        Minimum @f$ v_z@f$ (in centimeters)
  //    vzMax        Maximum @f$ v_z@f$ (in centimeters) 
  //    cMin         Minimum centrality (in percent)
  //    cMax         Maximum centrality (in percent)
  //    hist         Histogram to fill 
  // 
  // Return:
  //    @c true if the event meets the requirements 
  //
  if (cMin < cMax && (cMin > fCentrality || cMax <= fCentrality)) {
    if (status) status->Fill(kWrongCentrality);
    return false;
  }
  FillTriggerHistogram(triggerMask, fTriggers, hist);
  
  // Check if we have an event of interest. 
  if (!IsTriggerBits(triggerMask)) {
    if (status) status->Fill(kWrongTrigger);
    return false;
  }
  return FilterEvent(vzMin,vzMax,hist,status,filterMask);
}

//____________________________________________________________________
void
AliAODForwardMult::Print(Option_t* option) const
{
  // Print this object 
  // 
  // Parameters: 
  //  option   Passed to TH1::Print 
  fHist.Print(option);
  UShort_t sys = GetSystem();
  TString  str = "unknown";
  switch (sys) { 
  case 1:  str = "pp"; break;
  case 2:  str = "PbPb"; break;
  case 3:  str = "pPb" ; break;
  }
  std::cout << "Ipz:         " << fIpZ << "cm " << (HasIpZ() ? "" : "in") 
	    << "valid\n"
	    << "Triggers:    " << GetTriggerString(fTriggers)  << "\n"
	    << "sNN:         " << GetSNN() << "GeV\n" 
	    << "System:      " << str << "\n"
	    << "Centrality:  " << fCentrality << "%" 
	    << std::endl;
}

//____________________________________________________________________
//
// EOF
//
