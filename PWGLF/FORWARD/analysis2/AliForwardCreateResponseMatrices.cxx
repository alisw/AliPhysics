/**
 * @file   AliForwardCreateResponseMatrices.cxx
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Feb  7 01:01:24 2013
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_multdist
 */

#include <TH1D.h>
#include "AliForwardCreateResponseMatrices.h"
#include "AliForwardMultiplicityDistribution.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"
#include "AliFMDMCEventInspector.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"


ClassImp(AliForwardCreateResponseMatrices)
#if 0
; // This is for Emacs - do not delete
#endif
//_____________________________________________________________________
AliForwardCreateResponseMatrices::AliForwardCreateResponseMatrices() 
  : AliBaseAODTask(),
    fBins(), 
    fIsSelected(false)
{
  //
  // Default Constructor
  //
}
 
//_____________________________________________________________________
AliForwardCreateResponseMatrices::AliForwardCreateResponseMatrices(const char* name) 
  : AliBaseAODTask(name),
    fBins(), 
    fIsSelected(false)
{
  //
  // Constructor
  //
}

//_____________________________________________________________________
Bool_t 
AliForwardCreateResponseMatrices::CheckEvent(const AliAODForwardMult& fwd)
{
  fIsSelected = AliBaseAODTask::CheckEvent(fwd);
  // We always return true, so that we can process MC truth;
  return true;
}

//_____________________________________________________________________
Bool_t AliForwardCreateResponseMatrices::Book()
{
  //
  // Create Output Objects
  //
  TH1D* dndetaSumForward = new TH1D("dndetaSumForward","dndetaSumForward", 
				    200,-4,6);
  TH1D* dndetaSumCentral = new TH1D("dndetaSumCentral","dndetaSumCentral", 
				    200,-4,6);
  TH1D* dndetaSumMC = new TH1D("dndetaSumMC","dndetaSumMC", 200,-4,6);
 
  fSums->Add(dndetaSumForward);
  fSums->Add(dndetaSumCentral);
  fSums->Add(dndetaSumMC);

  TH1D* dndetaEventForward = new TH1D();
  TH1D* dndetaEventCentral = new TH1D();
  TH1D* dndetaEventMC = new TH1D();
  
  fSums->Add(dndetaEventForward);
  fSums->Add(dndetaEventCentral);
  fSums->Add(dndetaEventMC);

  //Loop over all individual eta bins, and define their hisograms 
  TIter next(&fBins);
  Bin * bin = 0;
  while ((bin = static_cast<Bin*>(next()))) { 
    bin->CreateOutputObjectss(fSums, 400);
  }
 
  return true;
  
}


//_____________________________________________________________________
Bool_t AliForwardCreateResponseMatrices::Event(AliAODEvent& aod)
{
  //
  // User Exec
  //
  // Here, we used to get ForwardMC!
  AliAODForwardMult* AODforward = GetForward(aod);
  if (!AODforward) return false;

  // Here, we used to get CentralClusters!
  AliAODCentralMult* AODcentral = GetCentral(aod);
  if (!AODcentral) return false;

  TH2D& forward = AODforward->GetHistogram();
  TH2D& central = AODcentral->GetHistogram();
  
  //check if event is NSD according to either MC truth or analysis (ESD)
  Bool_t isMCNSD  = kFALSE;
  Bool_t isESDNSD = kFALSE;
  if(AODforward->IsTriggerBits(AliAODForwardMult::kMCNSD))   isMCNSD=kTRUE; 
  if(AODforward->IsTriggerBits(AliAODForwardMult::kNSD))     isESDNSD=kTRUE;
  
  //primary dndeta dist
  TH2D*    primHist   = GetPrimary(aod);
  
  TH1D* dndetaSumForward    = (TH1D*)fSums->FindObject("dndetaSumForward");
  TH1D* dndetaSumCentral    = (TH1D*)fSums->FindObject("dndetaSumCentral");
  TH1D* dndetaSumMC         = (TH1D*)fSums->FindObject("dndetaSumMC");
  TH1D* dndetaEventForward  = forward.ProjectionX("dndetaForward",1,
						  forward.GetNbinsY(),"");
  TH1D* dndetaEventCentral  = central.ProjectionX("dndetaCentral",1,
						  central.GetNbinsY(),"");
  TH1D* dndetaEventMC       = primHist->ProjectionX("dndetaMC",1,
						    primHist->GetNbinsY(),"");
   
  // underflow eta bin of forward/central carry information on whether
  // there is acceptance. 1= acceptance, 0= no acceptance
  TH1D* normEventForward = 0;
  TH1D* normEventCentral = 0;
  normEventForward       = forward.ProjectionX("dndetaEventForwardNorm",0,0,"");
  normEventCentral       = central.ProjectionX("dndetaEventCentralNorm",0,0,"");

  dndetaSumForward->Add(dndetaEventForward);
  dndetaSumCentral->Add(dndetaEventCentral);
  dndetaSumMC->Add(dndetaEventMC);
     
  
  Double_t VtxZ = AODforward->GetIpZ();

  // loop over all eta bins, and create response matrices,
  // trigger-vertex bias histograms etc
  TIter next(&fBins);
  Bin * bin = 0;
  while ((bin = static_cast<Bin*>(next()))) { 
    bin->Process(dndetaEventForward, dndetaEventCentral, 
		 normEventForward,   normEventCentral, 
		 dndetaEventMC, VtxZ, fIsSelected,
		 isMCNSD, isESDNSD,  aod);
  }
    
  return true;  
}


//=====================================================================
AliForwardCreateResponseMatrices::Bin::Bin(Double_t etaLow, Double_t etaHigh) 
  : TNamed("", ""), 
    fEtaLow(etaLow),          // low eta limit 
    fEtaHigh(etaHigh),        // high eta limit 
    fHist(),                  // multiplicity histogram 
    fHistMC(),                // multiplicity histogram MC truth primaries
    fAcceptance(),            // histogram showing the 'holes' in acceptance. BinContent of 1 shows a hole, and BinContent of 10 shows data coverage
    fVtxZvsNdataBins(),       // VtxZ vs. number of data acceptance bins (normalised to the eta range)
    fResponseMatrix(),        // Response matrix (MC truth vs. analysed multiplicity)
    fResponseMatrixPlus05(),  // Response matrix with analysed multiplicity scaled up by 5%
    fResponseMatrixPlus075(), // Response matrix  with analysed multiplicity scaled up by 7.5%
    fResponseMatrixPlus10(),  // Response matrix with analysed multiplicity scaled up by 10%
    fResponseMatrixMinus05(), // Response matrix with analysed multiplicity scaled down by 5%
    fResponseMatrixMinus075(),// Response matrix with analysed multiplicity scaled down by 7.55%
    fResponseMatrixMinus10(), // Response matrix with analysed multiplicity scaled down by 10%
    fResponseMatrixMinusSys(),// Response matrix with analysed multiplicity scaled up by event mult uncertainty
    fResponseMatrixPlusSys(), // Response matrix with analysed multiplicity scaled down by event mult uncertainty
    fESDNSD(),                // number of events found as NSD by the analysis vs. multiplicity
    fMCNSD(),                 // number of events found as NSD by the MC truth vs. multiplicity
    fMCESDNSD(),              // number of events found as NSD by both analysis and MC truth vs. multiplicity   
    fTriggerBias()             // histogram for trigger vertex bias correction
{
  //
  // Constructor
  //
  const char* name = AliForwardMultiplicityDistribution::FormBinName(etaLow,etaHigh);

  SetName(name);
  SetTitle(Form("%+4.1f < #eta < %+4.1f", fEtaLow, fEtaHigh));
  
}
//_____________________________________________________________________
AliForwardCreateResponseMatrices::Bin::Bin() 
  : TNamed(), 
    fEtaLow(),          // low eta limit 
    fEtaHigh(),        // high eta limit 
    fHist(),                  // multiplicity histogram 
    fHistMC(),                // multiplicity histogram MC truth primaries
    fAcceptance(),            // histogram showing the 'holes' in acceptance. BinContent of 1 shows a hole, and BinContent of 10 shows data coverage
    fVtxZvsNdataBins(),       // VtxZ vs. number of data acceptance bins (normalised to the eta range)
    fResponseMatrix(),        // Response matrix (MC truth vs. analysed multiplicity)
    fResponseMatrixPlus05(),  // Response matrix with analysed multiplicity scaled up by 5%
    fResponseMatrixPlus075(), // Response matrix  with analysed multiplicity scaled up by 7.5%
    fResponseMatrixPlus10(),  // Response matrix with analysed multiplicity scaled up by 10%
    fResponseMatrixMinus05(), // Response matrix with analysed multiplicity scaled down by 5%
    fResponseMatrixMinus075(),// Response matrix with analysed multiplicity scaled down by 7.55%
    fResponseMatrixMinus10(), // Response matrix with analysed multiplicity scaled down by 10%
    fResponseMatrixMinusSys(),// Response matrix with analysed multiplicity scaled up by event mult uncertainty
    fResponseMatrixPlusSys(), // Response matrix with analysed multiplicity scaled down by event mult uncertainty
    fESDNSD(),                // number of events found as NSD by the analysis vs. multiplicity
    fMCNSD(),                 // number of events found as NSD by the MC truth vs. multiplicity
    fMCESDNSD(),              // number of events found as NSD by both analysis and MC truth vs. multiplicity        
    fTriggerBias()             // histogram for trigger vertex bias correction
{
  //
  // Default Constructor
  //
}
//_____________________________________________________________________
void AliForwardCreateResponseMatrices::Bin::CreateOutputObjectss(TList* cont,  Int_t max)
{
  //
  // Define eta bin output histos
  //
  TList* out = new TList;
  out->SetName(GetName());
  cont->Add(out);
  
  fHist                    = new TH1D("mult", GetTitle(), max, -0.5, max-.5);
  fHistMC                  = new TH1D("multTruth", GetTitle(), max, -0.5, max-.5);
  fVtxZvsNdataBins         = new TH2D("VtxZvsNdataBins", "VtxZ vs dataAcceptance/etaRange;z-vtz;dataAcceptance/etaRange", 20, -10,10, 130,0,1.3);
  fAcceptance              = new TH2D("Acceptance","Acceptance;#eta;z-vtx",200,-4, 6 , 20,-10,10);
  fResponseMatrix          = new TH2D("responseMatrix","Response Matrix;MC_{truth};MC_{measured}", max, -0.5, max-.5, max, -0.5, max-.5);
  fResponseMatrixPlus05    = new TH2D("responseMatrixPlus05","Response Matrix;MC_{truth};MC_{measured}", max, -0.5, max-.5, max, -0.5, max-.5);
  fResponseMatrixPlus075   = new TH2D("responseMatrixPlus075","Response Matrix;MC_{truth};MC_{measured}", max, -0.5, max-.5, max, -0.5, max-.5);
  fResponseMatrixPlus10    = new TH2D("responseMatrixPlus10","Response Matrix;MC_{truth};MC_{measured}", max, -0.5, max-.5, max, -0.5, max-.5);
  fResponseMatrixMinus05   = new TH2D("responseMatrixMinus05","Response Matrix;MC_{truth};MC_{measured}", max, -0.5, max-.5, max, -0.5, max-.5);
  fResponseMatrixMinus075  = new TH2D("responseMatrixMinus075","Response Matrix;MC_{truth};MC_{measured}", max, -0.5, max-.5, max, -0.5, max-.5);
  fResponseMatrixMinus10   = new TH2D("responseMatrixMinus10","Response Matrix;MC_{truth};MC_{measured}", max, -0.5, max-.5, max, -0.5, max-.5);
  fResponseMatrixMinusSys  = new TH2D("responseMatrixMinusSys","Response Matrix;MC_{truth};MC_{measured}", max, -0.5, max-.5, max, -0.5, max-.5);
  fResponseMatrixPlusSys   = new TH2D("responseMatrixPlusSys","Response Matrix;MC_{truth};MC_{measured}", max, -0.5, max-.5, max, -0.5, max-.5);
  fTriggerBias             = new TH1D("triggerBias","triggerBias;MC_{truth};Correction}", max, -0.5, max-.5);
  fMCNSD                   = new TH1D("fMCNSD","fMCNSD", 300,-0.5,299.5);
  fESDNSD                  = new TH1D("fESDNSD","fESDNSD", 300,-0.5,299.5);
  fMCESDNSD                = new TH1D("fMCESDNSD","fMCESDNSD", 300,-0.5,299.5);
  
  out->Add(fMCNSD);
  out->Add(fESDNSD);
  out->Add(fMCESDNSD);
  out->Add(fHist);
  out->Add(fHistMC);
  out->Add(fVtxZvsNdataBins);
  out->Add(fAcceptance);
  out->Add(fResponseMatrix);
  out->Add(fResponseMatrixPlus05);
  out->Add(fResponseMatrixPlus075);
  out->Add(fResponseMatrixPlus10);
  out->Add(fResponseMatrixMinus05);
  out->Add(fResponseMatrixMinus075);
  out->Add(fResponseMatrixMinus10);
  out->Add(fResponseMatrixPlusSys);
  out->Add(fResponseMatrixMinusSys);
  out->Add(fTriggerBias);
  
}
 

//_____________________________________________________________________
void 
AliForwardCreateResponseMatrices::Bin::Process(TH1D* dndetaForward, 
					       TH1D* dndetaCentral,
					       TH1D* normForward,   
					       TH1D* normCentral, 
					       TH1D* mc, 
					       Double_t VtxZ, 
					       Bool_t selectedTrigger, 
					       Bool_t isMCNSD, 
					       Bool_t isESDNSD, 
					       const AliAODEvent& aodevent) 
{
  //
  // Process a single eta bin
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
  fVtxZvsNdataBins->Fill(VtxZ, (Double_t)acceptanceBins/(last-first+1));



  Double_t c        = 0;
  Double_t e2       = 0;
  Double_t cPrimary = 0;
  Double_t e2Primary= 0;
    
  for (Int_t i = first; i <= last; i++){ 
    Double_t cForward  = 0;
    Double_t cCentral  = 0;
    Double_t e2Forward = 0;
    Double_t e2Central = 0;
    Double_t cMC  = 0;
    Double_t e2MC = 0;
    cMC= mc->GetBinContent(i);
    e2MC= mc->GetBinError(i) * mc->GetBinError(i);
    cPrimary+=cMC;
    e2Primary+=e2MC;
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
  
  // retreive MC particles from event
  TClonesArray* mcArray = (TClonesArray*)aodevent.FindListObject(AliAODMCParticle::StdBranchName());
  if(!mcArray){
    AliWarning("No MC array found in AOD. Try making it again.");
    return;
  }
  AliAODMCHeader* header = 
    dynamic_cast<AliAODMCHeader*>(aodevent.
				  FindListObject(AliAODMCHeader::StdBranchName()));
  if (!header) {
    AliWarning("No header file found.");
    return;
  }
  
  
  Int_t ntracks = mcArray->GetEntriesFast();
  // Track loop: find MC truth multiplicity in selected eta bin
  Double_t mult = 0;
  for (Int_t it = 0; it < ntracks; it++) {
    AliAODMCParticle* particle = (AliAODMCParticle*)mcArray->At(it);
    if (!particle) {
      AliError(Form("Could not receive track %d", it));
      continue;
    }
    if (!particle->IsPhysicalPrimary()) continue;
    if (particle->Charge() == 0) continue;
    if(particle->Eta()>fEtaLow&&particle->Eta()<fEtaHigh-0.0001)
      mult++;
  }
  //fill fMCNSD with multiplicity of MC truth NSD events with MC truth |vtxz|<4
  if(header->GetVtxZ()>-4&&header->GetVtxZ()<4){
    if(isMCNSD){
      fMCNSD->Fill(mult);
    }
  }
  //fill fESDNSD with multiplicity from events with a reconstructed NSD trigger and reconstructed |vtxz|<4
  if(VtxZ>-4&&VtxZ<4){
    if(isESDNSD){
      fESDNSD->Fill(mult);
    }
  }
  // fullfilling both requirements of fMCNSD and fESDNSD
  if(header->GetVtxZ()>-4&&header->GetVtxZ()<4&&VtxZ>-4&&VtxZ<4&&isMCNSD&&isESDNSD){
    fMCESDNSD->Fill(mult);
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
  Double_t spdSysError= 0.04;
  Double_t total = 0;
  total= fmd + spd + overlap; 
  if(total){  
    // Combined systematc event uncertainty, by weighting with the number of eta-bins of fmd-only, spd-only and the overlap.
    sysErrorSquared= (fmd*TMath::Power(fmdSysError,2)+ spd*TMath::Power(spdSysError,2)+
		      0.5*overlap*TMath::Power(fmdSysError,2)+ 0.5*overlap*TMath::Power(spdSysError,2))/total;
  }
  
  
  if(selectedTrigger){
    fHist->Fill(c);
    fHistMC->Fill(cPrimary);
    fResponseMatrix->Fill(cPrimary, c);
    fResponseMatrixPlusSys->Fill(cPrimary,(1+TMath::Sqrt(sysErrorSquared))*c);
    fResponseMatrixMinusSys->Fill(cPrimary,(1-TMath::Sqrt(sysErrorSquared))*c);
    fResponseMatrixPlus05->Fill(cPrimary, 1.05*c);
    fResponseMatrixPlus075->Fill(cPrimary, 1.075*c);
    fResponseMatrixPlus10->Fill(cPrimary, 1.1*c);
    fResponseMatrixMinus05->Fill(cPrimary, 0.95*c);
    fResponseMatrixMinus075->Fill(cPrimary, 0.925*c);
    fResponseMatrixMinus10->Fill(cPrimary, 0.9*c);
     
  }
  
}




//_____________________________________________________________________
//
//
// EOF
