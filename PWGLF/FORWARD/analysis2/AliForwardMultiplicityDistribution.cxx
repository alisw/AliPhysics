/**
 * @file   AliForwardMultiplicityDistribution.cxx
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Feb  7 01:03:52 2013
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_multdist
 * 
 */
#include <TH1D.h>
#include "AliForwardMultiplicityDistribution.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"
#include "TString.h"

ClassImp(AliForwardMultiplicityDistribution)
#if 0
; // This is for Emacs - do not delete
#endif
//______________________________________________________________________
AliForwardMultiplicityDistribution::AliForwardMultiplicityDistribution() 
  : AliBasedNdetaTask(), 
    fTrigger(0),   // trigger histogram
    fBins(),       // eta bin list
    fOutput(0),    // output list
    fLowCent(-1),  // lower centrality limit
    fHighCent(-1), // upper centrality limit
    fNBins(-1),    // multiplicity axis' runs from 0 to fNbins
    fCent(0)       // centrality
{
  //
  // Default Constructor
  //
}
 
//_____________________________________________________________________
AliForwardMultiplicityDistribution::AliForwardMultiplicityDistribution(const char* name) 
  : AliBasedNdetaTask(name),
    fTrigger(0),   // trigger histogram
    fBins(),       // eta bin list
    fOutput(0),    // output list
    fLowCent(-1),  // lower centrality limit
    fHighCent(-1), // upper centrality limit
    fNBins(-1),    // multiplicity axis' runs from 0 to fNbins
    fCent(0)       // centrality
{
  //
  // Constructor
  //
  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________
void AliForwardMultiplicityDistribution::UserCreateOutputObjects()
{
  //
  // Create Output Objects
  //
  fOutput = new TList;
  fOutput->SetOwner();
  fOutput->SetName(GetName());
  
  TH1D* dndetaSumForward = new TH1D("dndetaSumForward","dndetaSumForward", 200,-4,6);
  TH1D* dndetaSumCentral = new TH1D("dndetaSumCentral","dndetaSumCentral", 200,-4,6);
  fOutput->Add(dndetaSumForward);
  fOutput->Add(dndetaSumCentral);
  
  fCent = new TH1D("cent", "Centrality", 100, 0, 100);
  fCent->SetDirectory(0);
  fCent->SetXTitle(0);
  fOutput->Add(fCent);
  
  TH1D* dndetaEventForward = new TH1D();
  TH1D* dndetaEventCentral = new TH1D();
  fOutput->Add(dndetaEventForward);
  fOutput->Add(dndetaEventCentral);
  
  //make trigger diagnostics histogram
  fTrigger= new TH1I();
  fTrigger= AliAODForwardMult::MakeTriggerHistogram("triggerHist");
  fOutput->Add(fTrigger);
  
  //Loop over all individual eta bins, and define their hisograms 
  TIter next(&fBins);
  Bin * bin = 0;
  while ((bin = static_cast<Bin*>(next()))) { 
    bin->CreateOutputObjectss(fOutput, fNBins);
  }
  //fOutput->ls();
  PostData(1, fOutput);
}


//_____________________________________________________________________
void AliForwardMultiplicityDistribution::UserExec(Option_t */*option*/)
{
  //
  // User Exec
  //
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aod) {
    AliError("Cannot get the AOD event");
    return;
  } 
  // Check AOD event
  AliAODForwardMult* AODforward = static_cast<AliAODForwardMult*>(aod->FindListObject("Forward"));
  if (!AODforward->CheckEvent(fTriggerMask, fVtxMin, fVtxMax, 0,0,fTrigger)){
    AliError("Invalid Event");
    return;
  }
  
  // Fill centrality histogram, only useful for PbPb 
  Float_t cent = AODforward->GetCentrality();
  if(fLowCent<fHighCent){
    if(cent<fLowCent||cent>fHighCent) return;
  }
  fCent->Fill(cent);

  TH2D forward;
  TH2D central;
  
  // Get 2D eta-phi histograms for each event
  GetHistograms(aod, forward, central);
  
  TH1D* dndetaSumForward   = (TH1D*)fOutput->FindObject("dndetaSumForward");
  TH1D* dndetaSumCentral   = (TH1D*)fOutput->FindObject("dndetaSumCentral");
  TH1D* dndetaEventForward = 0;//(TH1D*)fOutput->FindObject("dndetaEventForward");
  TH1D* dndetaEventCentral = 0;//(TH1D*)fOutput->FindObject("dndetaEventCentral");
  TH1D* normEventForward   = 0;
  TH1D* normEventCentral   = 0;

  dndetaEventForward = forward.ProjectionX("dndetaForward",1,forward.GetNbinsY(),"");
  dndetaEventCentral = central.ProjectionX("dndetaCentral",1,central.GetNbinsY(),"");
  dndetaSumForward->Add(dndetaEventForward);
  dndetaSumCentral->Add(dndetaEventCentral);
  
  // underflow eta bin of forward/central carry information on whether there is acceptance. 1= acceptance, 0= no acceptance
  normEventForward   = forward.ProjectionX("dndetaEventForwardNorm",0,0,"");
  normEventCentral   = central.ProjectionX("dndetaEventCentralNorm",0,0,"");
  
  Double_t VtxZ = AODforward->GetIpZ();

  
  // loop over all eta bins, and fill multiplicity histos
  TIter next(&fBins);
  Bin * bin = 0;
  while ((bin = static_cast<Bin*>(next()))) { 
    bin->Process(dndetaEventForward, dndetaEventCentral, 
		 normEventForward,   normEventCentral, VtxZ);
  }
  
  PostData(1, fOutput);
  
}

//_____________________________________________________________________
void AliForwardMultiplicityDistribution::Terminate(Option_t */*option*/)
{
  //
  // Terminate
  //
}

//_________________________________________________________________-
TH2D* AliForwardMultiplicityDistribution::GetHistogram(const AliAODEvent*, Bool_t)
{
  //
  // implementation of pure virtual function, always returning 0
  //
  return 0;
}

void AliForwardMultiplicityDistribution::GetHistograms(const AliAODEvent* aod, TH2D& forward, TH2D& central, Bool_t /*mc*/)
{
  //
  // Get single event forward and central dN²/dEta dPhi histograms 
  //
  TObject* forwardObj = 0;
  TObject* centralObj = 0;

    
  forwardObj = aod->FindListObject("Forward");
  centralObj = aod->FindListObject("CentralClusters");

  AliAODForwardMult* AODforward = static_cast<AliAODForwardMult*>(forwardObj);
  AliAODCentralMult* AODcentral = static_cast<AliAODCentralMult*>(centralObj);
  
  forward= AODforward->GetHistogram();
  central= AODcentral->GetHistogram();
  
}

//______________________________________________________________________
#if 1
const Char_t*
#else
TString
#endif
AliForwardMultiplicityDistribution::FormBinName(Double_t etaLow, Double_t etaHigh)
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
 #if 0
  TString* combined= new TString();
  combined->Append(TString::Format("%s_%s", sLow.Data(), sHigh.Data()));
  return combined->Data();
#else 
  static TString tmp;
  tmp = TString::Format("%s_%s", sLow.Data(), sHigh.Data());
  return tmp.Data();
#endif
}

//=====================================================================
AliForwardMultiplicityDistribution::Bin::Bin() 
  : TNamed(), 
    fEtaLow(-1111),     // low eta limit 
    fEtaHigh(-1111),    // high eta limit 
    fHist(0),           // multiplicity distribution hist 
    fHistPlus05(0),     // multiplicity distribution hist scaled up with 5%
    fHistPlus075(0),    // multiplicity distribution hist scaled up with 7.5%
    fHistPlus10(0),     // multiplicity distribution hist scaled up with 10%
    fHistMinus05(0),    // multiplicity distribution hist scaled down with 5%
    fHistMinus075(0),   // multiplicity distribution hist scaled down with 7.5%
    fHistMinus10(0),    // multiplicity distribution hist scaled down with 10% 
    fHistPlusSys(0),    // multiplicity distribution hist scaled up with the event uncertainty
    fHistMinusSys(0),   // multiplicity distribution hist scaled down with the event uncertainty
    fAcceptance(0),     // histogram showing the 'holes' in acceptance. 
    fVtxZvsNdataBins(0) // VtxZ vs. number of data acceptance bins (normalised to the eta range) 
{
  // 
  // Default Constructor
  //
}
//_____________________________________________________________________
AliForwardMultiplicityDistribution::Bin::Bin(Double_t etaLow, Double_t etaHigh) 
  : TNamed("", ""), 
    fEtaLow(etaLow),    // low eta limit 
    fEtaHigh(etaHigh),  // high eta limit 
    fHist(0),           // multiplicity distribution hist 
    fHistPlus05(0),     // multiplicity distribution hist scaled up with 5%
    fHistPlus075(0),    // multiplicity distribution hist scaled up with 7.5%
    fHistPlus10(0),     // multiplicity distribution hist scaled up with 10%
    fHistMinus05(0),    // multiplicity distribution hist scaled down with 5%
    fHistMinus075(0),   // multiplicity distribution hist scaled down with 7.5%
    fHistMinus10(0),    // multiplicity distribution hist scaled down with 10% 
    fHistPlusSys(0),    // multiplicity distribution hist scaled up with the event uncertainty
    fHistMinusSys(0),   // multiplicity distribution hist scaled down with the event uncertainty
    fAcceptance(0),     // histogram showing the 'holes' in acceptance. 
    fVtxZvsNdataBins(0) // VtxZ vs. number of data acceptance bins (normalised to the eta range) 
{
  // 
  // Constructor
  //
  const char* name = AliForwardMultiplicityDistribution::FormBinName(fEtaLow,fEtaHigh);

  SetName(name);
  SetTitle(TString::Format("%+4.1f < #eta < %+4.1f", fEtaLow, fEtaHigh));

}
//_____________________________________________________________________
void AliForwardMultiplicityDistribution::Bin::CreateOutputObjectss(TList* cont,  Int_t max)
{
  //
  // Define eta bin output histograms
  //
  TList* out = new TList;
  out->SetName(GetName());
  cont->Add(out);
 
  fHist            = new TH1D("mult", GetTitle(), max, -0.5, max-.5);
  fHistPlus05      = new TH1D("multPlus05", GetTitle(), max, -0.5, max-.5);
  fHistPlus075     = new TH1D("multPlus075", GetTitle(), max, -0.5, max-.5);
  fHistPlus10      = new TH1D("multPlus10", GetTitle(), max, -0.5, max-.5);
  fHistMinus05     = new TH1D("multMinus05", GetTitle(), max, -0.5, max-.5);
  fHistMinus075    = new TH1D("multMinus075", GetTitle(), max, -0.5, max-.5);
  fHistMinus10     = new TH1D("multMinus10", GetTitle(), max, -0.5, max-.5);
  fHistPlusSys     = new TH1D("multPlusSys", GetTitle(), max, -0.5, max-.5);
  fHistMinusSys    = new TH1D("multMinusSys", GetTitle(), max, -0.5, max-.5);
  fVtxZvsNdataBins = new TH2D("VtxZvsNdataBins", "VtxZ vs dataAcceptance/etaRange;z-vtz;dataAcceptance/etaRange", 20, -10,10, 130,0,1.3);
  fAcceptance      = new TH2D("Acceptance","Acceptance;#eta;z-vtx",200,-4, 6 , 20,-10,10);

  out->Add(fHist);
  out->Add(fHistPlus05);
  out->Add(fHistPlus075);
  out->Add(fHistPlus10);
  out->Add(fHistMinus05);
  out->Add(fHistMinus075);
  out->Add(fHistMinus10);
  out->Add(fHistPlusSys);
  out->Add(fHistMinusSys);
  out->Add(fVtxZvsNdataBins);
  out->Add(fAcceptance);

}
 
//_____________________________________________________________________
void AliForwardMultiplicityDistribution::Bin::Process(TH1D* dndetaForward, TH1D* dndetaCentral,
						      TH1D* normForward,   TH1D* normCentral, Double_t VtxZ) 
{
  //
  // Process single eta bin
  //
  
  // Diagnostics on event acceptance
  Int_t    first = dndetaForward->GetXaxis()->FindBin(fEtaLow);
  Int_t    last  = dndetaForward->GetXaxis()->FindBin(fEtaHigh-.0001);
  Double_t acceptanceBins=0;
  for(Int_t n= first;n<=last;n++){
    if(normForward->GetBinContent(n)>0||normCentral->GetBinContent(n)>0){
      acceptanceBins++;
    }
    fAcceptance->SetBinContent(n, fAcceptance->GetYaxis()->FindBin(VtxZ), 1);
    if(normForward->GetBinContent(n)>0||normCentral->GetBinContent(n)>0)
      fAcceptance->SetBinContent(n, fAcceptance->GetYaxis()->FindBin(VtxZ),10);
  }

  //std::cout << "% of bins covered in eta-range: " << acceptanceBins/(last-first+1) << std::endl;
  fVtxZvsNdataBins->Fill(VtxZ, (Double_t)acceptanceBins/(last-first+1));
  

  Double_t c     = 0;
  Double_t e2    = 0;
  for (Int_t i = first; i <= last; i++) { 
    Double_t cForward  = 0;
    Double_t cCentral  = 0;
    Double_t e2Forward = 0;
    Double_t e2Central = 0;
    if (normForward->GetBinContent(i) > 0) {
      cForward  = dndetaForward->GetBinContent(i);
      e2Forward = dndetaForward->GetBinError(i) * dndetaForward->GetBinError(i);
    }
    if (normCentral->GetBinContent(i) > 0) { 
      cCentral  = dndetaCentral->GetBinContent(i);
      e2Central = dndetaCentral->GetBinError(i) * dndetaCentral->GetBinError(i);
    }
    Double_t cc  = 0;
    Double_t ee2 = 0;
    // if overlap between central/forward, use average of the two
    if (cCentral > 0 && cForward > 0) { 
      cc  = 0.5 * (cForward  + cCentral);
      ee2 = 0.5 * (e2Forward + e2Central);
    }
    else if (cCentral > 0) { 
      cc  = cCentral;
      ee2 = e2Central;
    }
    else if (cForward > 0) { 
      cc  = cForward;
      ee2 = e2Forward;
    }
    c  += cc;
    e2 += ee2;
  }
  
  //Systematic errors from here
  
  Double_t fmd=0;
  Double_t spd=0;
  Double_t overlap=0;

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
  // estimate of systematic uncertainty on the event multiplicity - hardcoded :(. estimates taken from Hans Hjersing Dalsgaard or Casper Nygaard phd theses.
  Double_t fmdSysError= 0.08;
  Double_t spdSysError= 0.067;
  
  Double_t total = 0;
  total= fmd+ spd + overlap; 
  // Combined systematc event uncertainty, by weighting with the number of eta-bins of fmd-only, spd-only and the overlap.
  if(total){
    sysErrorSquared= (fmd*TMath::Power(fmdSysError,2)+ spd*TMath::Power(spdSysError,2)+
		      0.5*overlap*TMath::Power(fmdSysError,2)+ 0.5*overlap*TMath::Power(spdSysError,2))/total;
  }
    
  //Strangeness correction, assumed to be 2.0 +- 1% (syst),  taken from HHD and CN thesis
  c=0.98*c;
  

  // correction for missing material description. Correction is based on separate PbPb analysis of difference between nominel and displaced vertices
  if(fEtaHigh<1.55&&fEtaHigh>1.45)
    c=0.98*c;
  if(fEtaHigh<2.05&&fEtaHigh>1.95)
    c=0.963*c;
  if(fEtaHigh<2.45&&fEtaHigh>2.35)
    c=0.956*c;
  if(fEtaHigh>2.95)
    c=0.955*c;
  
  fHist->Fill(c);
  fHistPlusSys->Fill((1+TMath::Sqrt(sysErrorSquared))*c);
  fHistMinusSys->Fill((1-TMath::Sqrt(sysErrorSquared))*c);
  fHistPlus05->Fill(1.05*c);
  fHistPlus075->Fill(1.075*c);
  fHistPlus10->Fill(1.1*c);
  fHistMinus05->Fill(0.95*c);
  fHistMinus075->Fill(0.925*c);
  fHistMinus10->Fill(0.9*c);

  
}
 



//_____________________________________________________________________
//
//
// EOF
