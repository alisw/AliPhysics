/**
 * @file   AliBaseMultTask
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Feb  7 01:01:24 2013
 * 
 * @brief  Base class for multiplicity distribution tasks 
 * 
 * 
 * @ingroup pwglf_forward_multdist
 */

#include <TH1D.h>
#include "AliBaseMultTask.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"


ClassImp(AliBaseMultTask)
#if 0
; // This is for Emacs - do not delete
#endif
//_____________________________________________________________________
AliBaseMultTask::AliBaseMultTask() 
  : AliBaseAODTask(),
    fBins(), 
    fIsSelected(false)
{
  //
  // Default Constructor
  //
}
 
//_____________________________________________________________________
AliBaseMultTask::AliBaseMultTask(const char* name) 
  : AliBaseAODTask(name,"AliBaseMultTask"),
    fBins(), 
    fIsSelected(false)
{
  //
  // Constructor
  //
}

//_____________________________________________________________________
void
AliBaseMultTask::DefaultBins()
{
  //Add Full eta-ranges
  AddBin(-3.4,5.1);
  
  //Add Symmetric eta bins.
  Double_t limits[] = { 3.4, 3.0, 2.5, 2.4, 2.0, 1.5, 1.4, 1.0, 0.5, 0. };
  Double_t* limit = limits;
  while ((*limit) > 0.1) {
    AddBin(-(*limit), +(*limit));
    limit++;
  }
}
//_____________________________________________________________________
Bool_t 
AliBaseMultTask::CheckEvent(const AliAODForwardMult& fwd)
{
  fIsSelected = AliBaseAODTask::CheckEvent(fwd);
  // We always return true, so that we can process MC truth;
  return true;
}
//_____________________________________________________________________
Bool_t AliBaseMultTask::Book()
{
  //
  // Create Output Objects
  //
  //Loop over all individual eta bins, and define their hisograms 
  TIter next(&fBins);
  Bin * bin = 0;
  while ((bin = static_cast<Bin*>(next()))) { 
    bin->CreateOutputObjects(fSums, 400);
  } 
  return true;  
}

//_____________________________________________________________________
Bool_t AliBaseMultTask::IsESDClass(AliAODForwardMult* m) const
{
  return m->IsTriggerBits(AliAODForwardMult::kNSD);
}
//_____________________________________________________________________
Bool_t AliBaseMultTask::IsMCClass(AliAODForwardMult* m) const
{
  return m->IsTriggerBits(AliAODForwardMult::kMCNSD);
}

//_____________________________________________________________________
Bool_t AliBaseMultTask::Event(AliAODEvent& aod)
{
  //
  // User Exec
  //
  // Here, we used to get ForwardMC!
  AliAODForwardMult* aodForward = GetForward(aod);
  if (!aodForward) return false;

  // Here, we used to get CentralClusters!
  AliAODCentralMult* aodCentral = GetCentral(aod);
  if (!aodCentral) return false;

  TH2D& forward = aodForward->GetHistogram();
  TH2D& central = aodCentral->GetHistogram();
  
  //check if event is NSD according to either MC truth or analysis (ESD)
  Bool_t isMCClass  = IsMCClass(aodForward);
  Bool_t isESDClass = IsESDClass(aodForward);
  Bool_t isPileUp   = aodForward->IsTriggerBits(AliAODForwardMult::kPileUp);
  //primary dndeta dist
  TH2D* primHist            = GetPrimary(aod);
  TH1D* dndetaForward  = forward.ProjectionX("dndetaForward",1,
					     forward.GetNbinsY(),"");
  TH1D* dndetaCentral  = central.ProjectionX("dndetaCentral",1,
					     central.GetNbinsY(),"");
  TH1D* dndetaMC       = (primHist ?
			  primHist->ProjectionX("dndetaMC",1,
						primHist->GetNbinsY(),"") :
			  0);
   
  // underflow eta bin of forward/central carry information on whether
  // there is acceptance. 1= acceptance, 0= no acceptance
  TH1D* normForward       = forward.ProjectionX("dndetaForwardNorm",0,0,"");
  TH1D* normCentral       = central.ProjectionX("dndetaCentralNorm",0,0,"");

  // Interaction point 
  Double_t ipZ = aodForward->GetIpZ();

  Process(dndetaForward,
	  dndetaCentral, 
	  normForward,
	  normCentral, 
	  dndetaMC,
	  ipZ,
	  isPileUp,
	  fIsSelected,
	  isMCClass,
	  isESDClass,
	  aod);
  // Clean up after projections 
  if (dndetaForward) delete dndetaForward;
  if (dndetaCentral) delete dndetaCentral;
  if (dndetaMC)      delete dndetaMC;
  
  return true;    
}

//_____________________________________________________________________
void AliBaseMultTask::Process(TH1D*              dndetaForward,
			      TH1D*              dndetaCentral,
			      TH1D*              normForward,
			      TH1D*              normCentral,
			      TH1D*              dndetaMC,
			      Double_t           ipZ,
			      Bool_t             pileup, 
			      Bool_t             selectedTrigger,
			      Bool_t             isMCClass,
			      Bool_t             isESDClass,
			      const AliAODEvent& aodevent)
{
  // loop over all eta bins, 
  TIter next(&fBins);
  Bin * bin = 0;
  while ((bin = static_cast<Bin*>(next()))) { 
    bin->Process(dndetaForward,
		 dndetaCentral, 
		 normForward,
		 normCentral, 
		 dndetaMC,
		 ipZ,
		 pileup,
		 fIsSelected,
		 isMCClass,
		 isESDClass,
		 aodevent,
		 fMinIpZ,
		 fMaxIpZ);
  }
}

//=====================================================================
AliBaseMultTask::Bin::Bin(Double_t etaLow, Double_t etaHigh) 
  : TNamed("", ""), 
    fEtaLow(etaLow),   // low eta limit 
    fEtaHigh(etaHigh), // high eta limit 
    fHist(),           // multiplicity histogram 
    fHistMC(),         // multiplicity histogram MC truth primaries
    fAcceptance(),     // histogram showing the 'holes' in
		       // acceptance. BinContent of 1 shows a hole,
		       // and BinContent of 10 shows data coverage
    fVtxZvsNdataBins() // VtxZ vs. number of data acceptance bins
{
  //
  // Constructor
  //
  const char* name = FormBinName(etaLow,etaHigh);
  SetName(name);
  SetTitle(Form("%+4.1f < #eta < %+4.1f", fEtaLow, fEtaHigh));
  
}
//_____________________________________________________________________
AliBaseMultTask::Bin::Bin() 
  : TNamed(), 
    fEtaLow(),          // low eta limit 
    fEtaHigh(),         // high eta limit 
    fHist(),            // multiplicity histogram 
    fHistMC(),          // multiplicity histogram MC truth primaries
    fAcceptance(),      // histogram showing the 'holes' in
			// acceptance. BinContent of 1 shows a hole,
			// and BinContent of 10 shows data coverage
    fVtxZvsNdataBins()  // VtxZ vs. number of data acceptance bins
			// (normalised to the eta range)
{
  //
  // Default Constructor
  //
}
//_____________________________________________________________________
void AliBaseMultTask::Bin::CreateOutputObjects(TList* cont,  Int_t max)
{
  //
  // Define eta bin output histos
  //
  TList* out = new TList;
  out->SetName(GetName());
  out->SetOwner();
  cont->Add(out);
  
  fHist             = new TH1D("mult", GetTitle(), max, -0.5, max-.5);
  fHistMC           = new TH1D("multTruth", GetTitle(), max, -0.5, max-.5);
  fVtxZvsNdataBins  = new TH2D("VtxZvsNdataBins",
			       "VtxZ vs dataAcceptance/etaRange;"
			       "z-vtz;dataAcceptance/etaRange",
			       20, -10,10, 130,0,1.3);
  fAcceptance       = new TH2D("Acceptance","Acceptance;#eta;z-vtx",
			       200,-4, 6 , 20,-10,10);

  out->Add(fHist);
  out->Add(fHistMC);
  out->Add(fVtxZvsNdataBins);
  out->Add(fAcceptance);  
}
 

//_____________________________________________________________________
Double_t
AliBaseMultTask::Bin::CalcMult(TH1D*              dndetaForward, 
			       TH1D*              dndetaCentral,
			       TH1D*              normForward,   
			       TH1D*              normCentral,
			       TH1D*              mc, 
			       Double_t           ipZ,
			       Double_t&          statErr,
			       Double_t&          sysErr,
			       Double_t&          mcMult, 
			       Double_t&          mcErr)
{
  // Find bin bounds 
  Int_t    first = dndetaForward->GetXaxis()->FindBin(fEtaLow);
  Int_t    last  = dndetaForward->GetXaxis()->FindBin(fEtaHigh-.0001);

  // Fill acceptance histogram 
  Double_t acceptanceBins=0;
  for(Int_t n= first;n<=last;n++){
    if(normForward->GetBinContent(n)>0||normCentral->GetBinContent(n)>0){
      acceptanceBins++;
    }
    fAcceptance->SetBinContent(n, fAcceptance->GetYaxis()->FindBin(ipZ), 1);
    if(normForward->GetBinContent(n)>0||normCentral->GetBinContent(n)>0)
      fAcceptance->SetBinContent(n, fAcceptance->GetYaxis()->FindBin(ipZ),10);
  }
  fVtxZvsNdataBins->Fill(ipZ, (Double_t)acceptanceBins/(last-first+1));

  // Zero output 
  statErr = 0;
  sysErr  = 0;
  mcMult  = -1;
  mcErr   = -1;
  // Sum up measurements
  Double_t c        = 0;
  Double_t e2       = 0;
  Double_t cPrimary = 0;
  Double_t e2Primary= 0;    
  for (Int_t i = first; i <= last; i++){ 
    Double_t cForward  = 0;
    Double_t cCentral  = 0;
    Double_t e2Forward = 0;
    Double_t e2Central = 0;
    Bool_t   aForward  = normForward->GetBinContent(i) > 0;
    Bool_t   aCentral  = normCentral->GetBinContent(i) > 0;
    if (aForward) {
      cForward  = dndetaForward->GetBinContent(i);
      e2Forward = TMath::Power(dndetaForward->GetBinError(i),2);
    }
    if (aCentral) { 
      cCentral  = dndetaCentral->GetBinContent(i);
      e2Central = TMath::Power(dndetaCentral->GetBinError(i),2);
    }
    Double_t cc  = 0;
    Double_t ee2 = 0;
    if (aCentral && aForward) { 
      cc  = 0.5 * (cForward  + cCentral);
      ee2 = 0.5 * (e2Forward + e2Central);
    }
    else if (aCentral) { 
      cc  = cCentral;
      ee2 = e2Central;
    }
    else if (aForward) { 
      cc  = cForward;
      ee2 = e2Forward;
    }
    cPrimary  += (mc ? mc->GetBinContent(i) : 0);
    e2Primary += (mc ? TMath::Power(mc->GetBinError(i),2) : 0);
    c         += cc;
    e2        += ee2;
  }

  // Systematic errors from here
  Int_t fmd=0;
  Int_t spd=0;
  Int_t overlap=0;
  // number of eta bins in fmd, spd and overlap respectively
  for(Int_t i = first;i<=last;i++){
    if(normForward->GetBinContent(i)>0&&normCentral->GetBinContent(i)<1)
      fmd++;
    if(normForward->GetBinContent(i)>0&&normCentral->GetBinContent(i)>0)
      overlap++;
    if(normCentral->GetBinContent(i)>0&&normForward->GetBinContent(i)<1)
      spd++;
  }
  
  Double_t sysErrorSquared = 0;  
  // estimate of systematic uncertainty on the event multiplicity -
  // hardcoded :(. estimates taken from Hans Hjersing Dalsgaard or
  // Casper Nygaard phd theses.
  Double_t fmdSysError= 0.08;
  Double_t spdSysError= 0.04;
  Double_t total = 0;
  total= fmd + spd + overlap; 
  if(total){  
    // Combined systematc event uncertainty, by weighting with the
    // number of eta-bins of fmd-only, spd-only and the overlap.
    sysErrorSquared= (fmd*TMath::Power(fmdSysError,2)+
		      spd*TMath::Power(spdSysError,2)+
		      0.5*overlap*TMath::Power(fmdSysError,2)+
		      0.5*overlap*TMath::Power(spdSysError,2))/total;
  }

  // Return values
  statErr = TMath::Sqrt(e2);
  sysErr  = TMath::Sqrt(sysErrorSquared);
  mcMult  = mc ? cPrimary : -1;
  mcErr   = mc ? TMath::Sqrt(e2Primary) : -1;
  
  return c;
}

//_____________________________________________________________________
void 
AliBaseMultTask::Bin::Process(TH1D*              dndetaForward, 
			      TH1D*              dndetaCentral,
			      TH1D*              normForward,   
			      TH1D*              normCentral, 
			      TH1D*              mc, 
			      Double_t           ipZ,
			      Bool_t             pileup,
			      Bool_t             selectedTrigger, 
			      Bool_t             isMCClass, 
			      Bool_t             isESDClass, 
			      const AliAODEvent& aodevent,
			      Double_t           minIPz,
			      Double_t           maxIPz) 
{
  //
  // Process a single eta bin
  //
  if (pileup) return;
  if (!selectedTrigger) return;
  if (ipZ < minIPz || ipZ > maxIPz) return;
  
  Double_t mcMult, mcErr, statErr, sysErr;
  Double_t mult = CalcMult(dndetaForward,
			   dndetaCentral,
			   normForward,
			   normCentral,
			   mc,
			   ipZ,
			   statErr,
			   sysErr,
			   mcMult,
			   mcErr);
  fHist->Fill(mult);
  if (mc) fHistMC->Fill(mcMult);
}
//______________________________________________________________________
const Char_t*
AliBaseMultTask::Bin::FormBinName(Double_t etaLow, Double_t etaHigh)
{
  //
  //  Form name of eta bin
  //
  TString sLow(TString::Format("%+03d",int(10*etaLow)));
  TString sHigh(TString::Format("%+03d",int(10*etaHigh)));
  sLow.ReplaceAll("+", "plus");
  sLow.ReplaceAll("-", "minus");
  sHigh.ReplaceAll("+", "plus");
  sHigh.ReplaceAll("-", "minus");
  static TString tmp;
  tmp = TString::Format("%s_%s", sLow.Data(), sHigh.Data());
  return tmp.Data();
}
//______________________________________________________________________
//
// EOF
//
