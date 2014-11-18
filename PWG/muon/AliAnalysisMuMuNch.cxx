#include "AliAnalysisMuMuNch.h"

/**
 *
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuNch
 *
 * SPD tracklet to Nch analysis
 *
 * The idea is that this sub-analysis is used first (within the AliAnalysisTaskMuMu sub-framework)
 * in order to compute the number of charged particles within an event, so that information
 * can be used in subsequent sub-analysis, like the invariant mass or mean pt ones.
 *
 */

#include "AliInputEventHandler.h"
#include "AliMultiplicity.h"
#include "AliAODTracklets.h"
#include "AliAnalysisMuonUtility.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliESDUtils.h"
#include "TMath.h"
#include "AliAnalysisMuMuCutRegistry.h"
#include "AliAnalysisMuMuEventCutter.h"
#include "AliAnalysisMuMuCutElement.h"
#include "AliAnalysisMuMuBinning.h"
#include "Riostream.h"
#include "TParameter.h"
#include <set>
#include <utility>
#include "AliMergeableCollection.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TObjString.h"
#include "TRandom3.h"

#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"


namespace {

  Double_t SPDgeomR(Double_t* x,Double_t* par) // Eta position of the SPD right edge eta as "seen" from a z vertex position
  {
    // par[0] = radius of SPD layer
    
    Double_t d = x[0];
    Double_t z(0.);
    Double_t theta(0.);
    
    if( d < 14.09999 )
    {
      z = 14.1 - d;
      theta = TMath::ATan(par[0]/z);
    }
    
    else if( d > 14.09999 )
    {
      z = d - 14.1;
      theta = TMath::Pi() - TMath::ATan(par[0]/z);
    }
    
    return -TMath::Log(TMath::Tan(theta/2.));
  }
  
  Double_t SPDgeomL(Double_t* x,Double_t* par) // Eta position of the SPD left edge eta as "seen" from a z vertex position
  {
    // par[0] = radius of SPD layer
    
    Double_t d = x[0];
    Double_t z(0.);
    Double_t theta(0.);
    
    if( d > -14.09999 )
    {
      z = 14.1 + d;
      theta = TMath::Pi() - TMath::ATan(par[0]/z);
    }
    
    if( d < -14.09999 )
    {
      z = -14.1 - d;
      theta = TMath::ATan(par[0]/z);
    }
    
    return -TMath::Log(TMath::Tan(theta/2.));
  }

}

////_____________________________________________________________________________
//AliAnalysisMuMuNch::AliAnalysisMuMuNch(TH2* spdCorrection, Double_t etaMin, Double_t etaMax
//                                      , Double_t zMin, Double_t zMax,Bool_t disableHistos, Bool_t computeResolution)
//: AliAnalysisMuMuBase(),
//fSPDCorrection(0x0),
//fEtaAxis(new TAxis(TMath::Nint(10./0.1),-5.,5.)),
//fZAxis(new TAxis(TMath::Nint(80/0.25),-40.,40.)),
//fCurrentEvent(0x0),
//fEtaMin(etaMin),
//fEtaMax(etaMax),
//fZMin(zMin),
//fZMax(zMax),
//fResolution(computeResolution)
//{
//  if ( spdCorrection )
//  {
//    fSPDCorrection = static_cast<TH2F*>(spdCorrection->Clone());
//    fSPDCorrection->SetDirectory(0);
//  }
//  DefineSPDAcceptance();
//  
//  if ( disableHistos )
//  {
//    DisableHistograms("*");
//  }
//}

//FIXME: First and second constructor may be ambiguous if we do not set all the arguments

//_____________________________________________________________________________
AliAnalysisMuMuNch::AliAnalysisMuMuNch(TH2F* spdCorrection, TProfile* spdMeanCorrection, Double_t meanTrRef, Double_t etaMin, Double_t etaMax
                                       , Double_t zMin, Double_t zMax, Bool_t disableHistos, Bool_t computeResolution)
: AliAnalysisMuMuBase(),
fSPDOneOverAccxEff(0x0),
fSPDCorrectionMap(0x0),
fSPDCorrectionList(0x0),
fSPDMeanTracklets(0x0),
fSPDMeanTrackletsCorrToCompare(0x0),
fEtaAxis(new TAxis(TMath::Nint(10./0.1),-5.,5.)),
fZAxis(new TAxis(TMath::Nint(80/0.25),-40.,40.)),
fCurrentEvent(0x0),
fMeanTrRef(meanTrRef),
fEtaMin(etaMin),
fEtaMax(etaMax),
fEtaMinToCompare(0.),
fEtaMaxToCompare(0.),
fetaRange(),
fZMin(zMin),
fZMax(zMax),
fResolution(computeResolution),
frand(0x0),
fGeneratorHeaderClass(new TString("AliGenDPMjetEventHeader"))
{
  //FIXME: Add a protection to avoid an etamin or etamax non multiple of the eta bin size
  
  //FIXME: Add fluctuations to the SPD AccxEff correction method
  
  if ( spdCorrection && spdMeanCorrection )
  {
    AliWarning("Two methods used to correct tracklets: Data-driven and SPDAccxEff corrected tracklets will be available in the event list, but the histograms in this class will be filled just with the AccxEff corrected values");
//    AliFatal("Cannot use 2 different methods for tracklets correction");
  }
//  else if ( spdCorrection )
  if ( spdCorrection )
  {
    AliWarning("SPD AccxEff correction method not ready to be used, fluctuations are missing. Do not trust the results");
    
    if ( fMeanTrRef > 0. && !spdMeanCorrection ) AliWarning("Reference mean nof tracklets argument will not be used: SPD AccxEff correction method in use");
    
    fSPDOneOverAccxEff = static_cast<TH2F*>(spdCorrection->Clone());
    fSPDOneOverAccxEff->SetDirectory(0);
  }
//  else if ( spdMeanCorrection )
  if ( spdMeanCorrection )
  {
    if ( fMeanTrRef < 0. ) AliWarning("Reference mean nof tracklets argument set to -1: Maximum of the spdMeanCorrection will be used");
    else AliWarning(Form("Using %f as reference mean nof tracklets for correction",fMeanTrRef));
    
    frand = new TRandom3();
    fSPDMeanTracklets = static_cast<TProfile*>(spdMeanCorrection->Clone());
    fSPDMeanTracklets->SetDirectory(0);
  }
  
  DefineSPDAcceptance();
  
  if ( disableHistos ) // FIXME: Is this really useful? it breaks when setting to ktrue due to non existence of histos in SetEvent(). Answer: It is useful, it will speed up the task when we want to execute only SetEvent and not fill all the multiplicity histos (i.e. doing J/psi vs multiplicity analysis). The problem is that disabling the histos the method CreateHisto() does not create histosm thats why the task breaks in SetEvent, so fix this
  {
    DisableHistograms("*");
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuNch::AliAnalysisMuMuNch(TProfile* spdMeanCorrection, TProfile* spdMeanCorrectionToCompare, Double_t meanTrRef, Double_t etaMin,
                                       Double_t etaMax, Double_t zMin, Double_t zMax, Double_t etaMinToCompare, Double_t etaMaxToCompare,Bool_t disableHistos, Bool_t computeResolution)
: AliAnalysisMuMuBase(),
fSPDOneOverAccxEff(0x0),
fSPDCorrectionMap(0x0),
fSPDCorrectionList(0x0),
fSPDMeanTracklets(0x0),
fSPDMeanTrackletsCorrToCompare(0x0),
fEtaAxis(new TAxis(TMath::Nint(10./0.1),-5.,5.)),
fZAxis(new TAxis(TMath::Nint(80/0.25),-40.,40.)),
fCurrentEvent(0x0),
fMeanTrRef(meanTrRef),
fEtaMin(etaMin),
fEtaMax(etaMax),
fEtaMinToCompare(etaMinToCompare),
fEtaMaxToCompare(etaMaxToCompare),
fetaRange(),
fZMin(zMin),
fZMax(zMax),
fResolution(computeResolution),
frand(0x0),
fGeneratorHeaderClass(new TString("AliGenDPMjetEventHeader"))
{
  //FIXME: Add a protection to avoid an etamin or etamax non multiple of the eta bin size
  
  /// This construction is designed to compute everything in etaMin < eta < etaMax and correct with spdMeanCorrection (main eta tange and correction), but when filling the histos in FillHistosForEvent is able to compute the number of tracklets in etaMinToCompare < eta < etaMaxToCompare and correct them with spdMeanCorrectionToCompare (secondary eta range and correction) in order to compare N_{tr}^{|eta|< etaMax} vs N_{tr}^{|eta|< etaMaxToCompare}
  
  if ( !spdMeanCorrection || !spdMeanCorrectionToCompare )
  {
    AliFatal("Need the two corrections to compare. Maybe you are using the wrong constructor");
  }
  if ( fEtaMinToCompare < fEtaMin || fEtaMaxToCompare > fEtaMax )
  {
    AliFatal("Cannot select a eta range to compare wider than the main eta range");
  }
  
  if ( fMeanTrRef < 0. ) AliWarning("Reference mean nof tracklets argument set to -1: Maximum of the spdMeanCorrection will be used");
  else AliWarning(Form("Using %f as reference mean nof tracklets for correction",fMeanTrRef));
  
  frand = new TRandom3();
  fSPDMeanTracklets = static_cast<TProfile*>(spdMeanCorrection->Clone());
  fSPDMeanTracklets->SetDirectory(0);
  
  fSPDMeanTrackletsCorrToCompare = static_cast<TProfile*>(spdMeanCorrectionToCompare->Clone());
  fSPDMeanTrackletsCorrToCompare->SetDirectory(0);
    
  DefineSPDAcceptance();
  
  if ( disableHistos ) // FIXME: Is this really useful? it breaks when setting to ktrue due to non existence of histos in SetEvent()
  {
    DisableHistograms("*");
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuNch::AliAnalysisMuMuNch(TObjArray* spdCorrectionList, Double_t meanTrRef, Double_t etaMin, Double_t etaMax
                                       , Double_t zMin, Double_t zMax,Bool_t disableHistos, Bool_t computeResolution)
: AliAnalysisMuMuBase(),
fSPDOneOverAccxEff(0x0),
fSPDCorrectionMap(0x0),
fSPDCorrectionList(0x0),
fSPDMeanTracklets(0x0),
fSPDMeanTrackletsCorrToCompare(0x0),
fEtaAxis(new TAxis(TMath::Nint(10./0.1),-5.,5.)),
fZAxis(new TAxis(TMath::Nint(80/0.25),-40.,40.)),
fCurrentEvent(0x0),
fMeanTrRef(meanTrRef),
fEtaMin(etaMin),
fEtaMax(etaMax),
fEtaMinToCompare(0.),
fEtaMaxToCompare(0.),
fetaRange(),
fZMin(zMin),
fZMax(zMax),
fResolution(computeResolution),
frand(0x0),
fGeneratorHeaderClass(new TString("AliGenDPMjetEventHeader"))
{
  //FIXME: Add a protection to avoid an etamin or etamax non multiple of the eta bin size
  
  // Uses a different correction for each group of runs (both SPD AccxEff OR mean tracklets are supported)
  
  if ( spdCorrectionList ) DefineSPDCorrectionMap(spdCorrectionList);
  else AliFatal("No SPD correction list provided");
  
  if ( fMeanTrRef < 0. ) AliWarning("Reference mean nof tracklets argument set to -1: Maximum of the spdMeanCorrection will be used");
  else AliWarning(Form("Using %f as reference mean nof tracklets for correction",fMeanTrRef));

  frand = new TRandom3();
  
  DefineSPDAcceptance();
  
  if ( disableHistos ) // FIXME: Is this really useful? it breaks when setting to ktrue due to non existence of histos in SetEvent()
  {
    DisableHistograms("*");
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuNch::~AliAnalysisMuMuNch()
{
  delete fSPDOneOverAccxEff;
  delete fSPDCorrectionMap;
  delete fSPDCorrectionList;
  delete fSPDMeanTracklets;
  delete fSPDMeanTrackletsCorrToCompare;
  delete fEtaAxis;
  delete fZAxis;
  delete fSPD1LR;
  delete fSPD1LL;
  delete fSPD2LR;
  delete fSPD2LL;
  delete frand;
  delete fGeneratorHeaderClass;
}

//_____________________________________________________________________________
void AliAnalysisMuMuNch::DefineGeneratorName(const char* genName)
{
//  TString sgenName(genName);
//  if ( sgenNam.Contains("pythia") ) fGeneratorHeaderClass = "AliGenPythiaEventHeader";
//  else if ( sgenNam.Contains("dpmjet") ) fGeneratorHeaderClass = "AliGenDPMjetEventHeader";
//  else if ( sgenNam.Contains("dpmjet") ) fGeneratorHeaderClass = "AliGenHijingEventHeader";
  
  fGeneratorHeaderClass->Form("%s",genName);
  
  std::cout << " Will use " << fGeneratorHeaderClass->Data() << " tracks" << std::endl;
}
                                                   
//_____________________________________________________________________________
void AliAnalysisMuMuNch::DefineHistogramCollection(const char* eventSelection,
                                                   const char* triggerClassName,
                                                   const char* centrality)
{
  // Define multiplicity histos
  
  if ( Histo(eventSelection,triggerClassName,centrality,"AliAnalysisMuMuNch") )
  {
    return;
  }

  // dummy histogram to signal that we already defined all our histograms (see above)
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"AliAnalysisMuMuNch","Dummy semaphore",1,0,1);

  Double_t multMin = -0.5;  //Tracklets multiplicity range
  Double_t multMax = 500.5;
  Int_t nbinsMult = GetNbins(multMin,multMax,1.);
  
  Double_t phimin = 0.; //Phi range
  Double_t phimax = 2*TMath::Pi();
  Int_t nphibins = GetNbins(phimin,phimax,0.05);

  if ( !fSPDMeanTracklets && !fSPDOneOverAccxEff && fResolution ) // Resolution histos
  {
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"EtaRes","#eta resolution;#eta_{Reco} - #eta_{MC};Counts",(fEtaAxis->GetNbins())*2,fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"PhiRes","#phi resolution;#phi_{Reco} - #phi_{MC};Counts",nphibins*2,-phimax,phimax);
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"PhiResShifted","#phi resolution;#phi_{Reco} - #phi_{MC};Counts",nphibins*4,-phimax/4,phimax/4);
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"EtaResVsZ","#eta resolution vs MC Zvertex;ZVertex (cm);#eta_{Reco} - #eta_{MC}",(fZAxis->GetNbins())*20,fZAxis->GetXmin(),fZAxis->GetXmax(),(fEtaAxis->GetNbins())*8,fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"PhiResVsZ","#phi resolution vs MC Zvertex;ZVertex (cm);#phi_{Reco} - #phi_{MC}",(fZAxis->GetNbins())*20,fZAxis->GetXmin(),fZAxis->GetXmax(),nphibins*4,-phimax/4,phimax/4);
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"EtaResVsnC","#eta resolution vs Nof contributors to SPD vertex;NofContributors;#eta_{Reco} - #eta_{MC}",200,0,200,(fEtaAxis->GetNbins())*8,fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"PhiResVsnC","#phi resolution vs Nof contributors to SPD vertex;NofContributors;#phi_{Reco} - #phi_{MC}",200,0,200,nphibins*4,-phimax/4,phimax/4);
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"SPDZvResVsnC","SPD Zvertex resolution vs Nof contributors;NofContributors;Zvertex_{Reco} - Zvertex_{MC} (cm)",200,0,200,(fZAxis->GetNbins())*20,fZAxis->GetXmin(),fZAxis->GetXmax());
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"SPDZvResVsMCz","SPD Zvertex resolution vs MC z vertex;MC Zvertex (cm);Zvertex_{Reco} - Zvertex_{MC} (cm)",(fZAxis->GetNbins())*20,fZAxis->GetXmin(),fZAxis->GetXmax(),(fZAxis->GetNbins())*10,fZAxis->GetXmin(),fZAxis->GetXmax());
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"Phi","Reco #phi distribution;#phi;Counts",nphibins,phimin,phimax);
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"MCPhi","MC #phi distribution;#phi;Counts",nphibins,phimin,phimax);
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"Eta","Reco #eta distribution;#phi;Counts",(fEtaAxis->GetNbins())*2,fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"MCEta","MC #eta distribution;#phi;Counts",(fEtaAxis->GetNbins())*2,fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
    
    return; // When computing resolutions we don't want to create the rest of histos
  }

  TString TrackletsCorrectedName("");
  TString TrackletsCorrectedAxisName("");
  if ( fSPDMeanTracklets && !fSPDOneOverAccxEff )
  {
    TrackletsCorrectedName = "Number of corrected tracklets";
    TrackletsCorrectedAxisName = "N_{tr}^{corr}";
  }
  else if( fSPDOneOverAccxEff )
  {
    TrackletsCorrectedName = "Number of charged particles";
    TrackletsCorrectedAxisName = "N_{ch}";
  }
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsPhi","Number of tracklets vs Z vertex vs #phi;Z vertex;#phi",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),nphibins,phimin,phimax);
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsEta","Number of tracklets vs ZVertex vs #eta;ZVertex;#eta",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
  AttachSPDAcceptance(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsEta");
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Tracklets",Form("Number of tracklets in |#eta| < %1.1f distribution;N_{Tracklets};N_{events}",fEtaMax),nbinsMult,multMin,multMax);
    
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"NBkgTrackletsVsZVertexVsEta","Number of background tracklets vs ZVertex vs #eta;ZVertex;#eta",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
  AttachSPDAcceptance(kHistoForMCInput,eventSelection,triggerClassName,centrality,"NBkgTrackletsVsZVertexVsEta");
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"NchVsZVertexVsEta","Number of charged particles vs ZVertex vs #eta;ZVertex;#eta",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"NchVsRecoZVertexVsEta","Number of gen. charged particles vs reco. ZVertex vs gen. #eta;ZVertex;#eta",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"NchVsZVertexVsEta",Form("%s vs ZVertex vs #eta;ZVertex;#eta",TrackletsCorrectedName.Data()),fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
  AttachSPDAcceptance(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"NchVsZVertexVsEta");
  AttachSPDAcceptance(kHistoForMCInput,eventSelection,triggerClassName,centrality,"NchVsRecoZVertexVsEta");
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"NchVsZVertexVsPhi","Number of charged particles vs ZVertex vs #phi;ZVertex;#phi",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),nphibins,phimin,phimax);
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"EventsVsZVertexVsEta","Effective number of events vs ZVertex vs #eta;ZVertex;#eta",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax()); // Fill 1 unit in each "touched" eta bin per event (represents the eta bins in which each event contributes)
  AttachSPDAcceptance(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"EventsVsZVertexVsEta");
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"Nch",Form("Number of charged particles in |#eta| < %1.1f distribution;N_{ch};N_{events}",fEtaMax),nbinsMult,multMin,multMax);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"Nch",Form("%s in |#eta| < %1.1f distribution;%s;N_{events}",TrackletsCorrectedName.Data(),fEtaMax,TrackletsCorrectedAxisName.Data()),nbinsMult,multMin,multMax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"dNchdEta",Form("<dNchd#eta> in |#eta| < %1.1f distribution;dN_{ch}/d#eta;N_{events}",fEtaMax),nbinsMult,multMin,multMax);
   CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"dNchdEta",Form("<d%sd#eta> in |#eta| < %1.1f distribution;d%s/d#eta;N_{events}",TrackletsCorrectedAxisName.Data(),fEtaMax,TrackletsCorrectedAxisName.Data()),nbinsMult,multMin,multMax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"TrackletsVsNch",Form("Number of tracklets vs number of generated charged particles in |#eta| < %1.1f;N_{ch};N_{tracklets}",fEtaMax),nbinsMult,multMin,multMax,nbinsMult,multMin,multMax); //Response matrix
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"TrackletsVsNch",Form("Number of tracklets vs %s in |#eta| < %1.1f;%s;N_{tracklets}",TrackletsCorrectedName.Data(),fEtaMax,TrackletsCorrectedAxisName.Data()),nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
  
   CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"dNchdetaVsMCdNchdeta",Form("Corrected dN_{ch}/d#eta vs MC dN_{ch}/d#eta in |#eta| < %1.1f;(dN_{ch}/d#eta)_{MC};dN_{ch}/d#eta",fEtaMax),nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"V0AMultVsNch",Form("V0A multiplicity vs number of generated charged particles in |#eta| < %1.1f;N_{ch};V0A Mult",fEtaMax),nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"V0AMultVsNch",Form("V0A multiplicity vs %s in |#eta| < %1.1f;%s;V0A Mult",TrackletsCorrectedName.Data(),fEtaMax,TrackletsCorrectedAxisName.Data()),nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"V0CMultVsNch",Form("V0C multiplicity vs number of generated charged particles in |#eta| < %1.1f;N_{ch};V0C Mult",fEtaMax),nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"V0CMultVsNch",Form("V0C multiplicity vs %s in |#eta| < %1.1f;%s;V0A Mult",TrackletsCorrectedName.Data(),fEtaMax,TrackletsCorrectedAxisName.Data()),nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
  
//  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"NchVsTracklets","Number of generated charged particles vs reco tracklets;N_{tracklets};N_{ch}",nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
  
  if ( fSPDMeanTrackletsCorrToCompare )
  {
    CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"CorrTrackletsEtaSecVsCorrTrackletsEtaPrim",Form("N_{tr}^{corr} in |#eta|< %1.1f vs N_{tr}^{corr} in |#eta|< %1.1f;N_{tracklets}^{|#eta|<%1.1f};N_{tracklets}^{|#eta|<%1.1f}",fEtaMaxToCompare,fEtaMax,fEtaMax,fEtaMaxToCompare),nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
    CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"MeanNchEtaSecVsZVertex",Form("Mean number of corrected tracklets in |#eta| < %1.1f vs Z vertex;Z vertex;<N_{tr}^{corr}>",fEtaMaxToCompare),fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),0);
    CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"MeanTrackletsEtaSecVsZVertex",Form("Mean number of tracklets in |#eta| < %1.1f vs Z vertex;Z vertex;<N_{tr}^{corr}>",fEtaMaxToCompare),fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),0);
    CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"TrackletsSecVsZVertexVsEta","Number of tracklets vs ZVertex vs #eta;ZVertex;#eta",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
  }
  
  if ( fSPDMeanTracklets && fSPDOneOverAccxEff )
  {
    CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"dNchdetaComparison2Corrections","dN_{ch}/d#eta comparison with SPD AccxEff and data-driven corrections;dN_{ch}/d#eta;A*N^{corr}_{Tracklets}",151,multMin,150.5,151,multMin,150.5);
    CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"DispersiondNchdetaComparison2Corrections","dN_{ch}/d#eta dispersion with SPD AccxEff and data-driven corrections;(dN_{ch}/d#eta - A*N^{corr}_{Tracklets}) /dN_{ch}/d#eta;N_{events}",201,-1.005,1.005);
    CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"CheckMeanNtrCorrVsZVertex","Check for NtrCorr",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),0);
    CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"CheckNtrCorr","Check for NtrCorr distribution;N^{corr}_{tr};N_{events}",nbinsMult,multMin,multMax);
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"RelDispersiondNchdetaFromNtrCorrVsdNchdEtaMC",Form("Relative dispersion of dN_{ch}/d#eta from N_{tr}^{corr} vs MC dN_{ch}/d#eta in |#eta| < %1.1f;((dN_{ch}/d#eta)_{gen} - A*N^{corr}_{Tracklets}) /(dN_{ch}/d#eta)_{gen}",fEtaMax),201,-1.005,1.005);
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"DispersiondNchdetaFromNtrCorrVsdNchdEtaMC",Form("Dispersion of dN_{ch}/d#eta from N_{tr}^{corr} vs MC dN_{ch}/d#eta in |#eta| < %1.1f;(dN_{ch}/d#eta)_{gen} - A*N^{corr}_{Tracklets}",fEtaMax),402,-100.5,100.5);
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"dNchdetaFromNtrCorrVsdNchdEtaMC",Form("dN_{ch}/d#eta from N_{tr}^{corr} vs MC dN_{ch}/d#eta in |#eta| < %1.1f;(dN_{ch}/d#eta)_{gen};dN_{ch}/d#eta",fEtaMax),nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
    
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"RelDispersiondNchdetaFromAccEffVsdNchdEtaMC",Form("Relative dispersion of dN_{ch}/d#eta from N_{tr}^{AccxEff} vs MC dN_{ch}/d#eta in |#eta| < %1.1f;((dN_{ch}/d#eta)_{gen} - dN_{ch}/d#eta) /(dN_{ch}/d#eta)_{gen};N_{events}",fEtaMax),201,-1.005,1.005);
     CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"DispersiondNchdetaFromAccEffVsdNchdEtaMC",Form("Dispersion of dN_{ch}/d#eta from N_{tr}^{AccxEff} vs MC dN_{ch}/d#eta in |#eta| < %1.1f;(dN_{ch}/d#eta)_{gen} - dN_{ch}/d#eta;N_{events}",fEtaMax),402,-100.5,100.5);
    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"dNchdetaFromAccEffVsdNchdEtaMC",Form("dN_{ch}/d#eta from N_{tr}^{AccxEff} vs MC dN_{ch}/d#eta in |#eta| < %1.1f;(dN_{ch}/d#eta)_{gen};dN_{ch}/d#eta",fEtaMax),nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
  }
  
//   CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"CorrTrackletsVsNch","Reco tracklets corrected vs number of generated charged particles;N_{ch};N_{tracklets}^{corr}",nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"CorrTrackletsVsNch",Form("%s vs number of generated charged particles in |#eta| < %1.1f;N_{ch};%s",TrackletsCorrectedName.Data(),fEtaMax,TrackletsCorrectedAxisName.Data()),nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
  
  // profile histograms
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"MeanTrackletsVsEta","Mean number of tracklets vs #eta;#eta;<N_{Tracklets}>",fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax(),0);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"MeandNchdEtaVsEta","Mean dN_{ch}/d#eta vs #eta;#eta;<dN_{ch}/d#eta>",fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax(),0);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"MeanTrackletsVsZVertex",Form("Mean number of tracklets in |#eta| < %1.1f vs Z vertex;Z vertex;<N_{Tracklets}>",fEtaMax),fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),0);

  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"MeanNchVsEta","Mean number of generated charged particles vs #eta;#eta;<N_{ch}>",fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax(),0); // Each bin has to be divided by the binwidth to became dNch/dEta (Done in the terminate and stored in MeandNchdEta histo)
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"MeanNchVsEta",Form("Mean %s vs #eta;#eta;<%s>",TrackletsCorrectedName.Data(),TrackletsCorrectedAxisName.Data()),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax(),0); // Each bin has to be divided by the binwidth to became dNch/dEta (Done in the terminate and stored in MeandNchdEta histo)
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"MeanNchVsZVertex",Form("Mean number of generated charged particles in |#eta| < %1.1f vs Z vertex;Z vertex;<N_{ch}>",fEtaMax),fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),0);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"MeanNchVsZVertex",Form("Mean %s in |#eta| < %1.1f vs Z vertex;Z vertex;<%s>",TrackletsCorrectedName.Data(),fEtaMax,TrackletsCorrectedAxisName.Data()),fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),0);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"MeandNchdEtaVsZVertex",Form("Mean d%s/d#eta in |#eta| < %1.1f vs Z vertex;Z vertex;<%s>",TrackletsCorrectedAxisName.Data(),fEtaMax,TrackletsCorrectedAxisName.Data()),fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),0);
  
  
  TObjArray* bins = Binning()->CreateBinObjArray("psi","dnchdeta","");
  TIter nextBin(bins);
  AliAnalysisMuMuBinning::Range* r;
  
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
  {
    CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,Form("EventsIn%s",r->AsString().Data()),Form("Number of %s events in %s bin",triggerClassName,r->AsString().Data()),1,0.,2);
  }

  delete bins;
  
  TObjArray* binsntr = Binning()->CreateBinObjArray("psi","ntrcorr","");
  TIter nextBinNtr(binsntr);
  AliAnalysisMuMuBinning::Range* rNtr;
  
  while ( ( rNtr = static_cast<AliAnalysisMuMuBinning::Range*>(nextBinNtr()) ) )
  {
    CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,Form("MeanNchVsZVertex%s",rNtr->AsString().Data()),Form("Mean %s in |#eta| < %1.1f vs Z vertex in bin %s;Z vertex;<%s>",TrackletsCorrectedName.Data(),fEtaMax,rNtr->AsString().Data(),TrackletsCorrectedAxisName.Data()),fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),0);
    
    CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,Form("MeanTrackletsVsZVertex%s",rNtr->AsString().Data()),Form("Mean number of tracklets in |#eta| < %1.1f vs Z vertex in bin %s;Z vertex;<N_{Tracklets}>",fEtaMax,rNtr->AsString().Data()),fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),0);
  }
  
  delete binsntr;
  
  TObjArray* binsRel = Binning()->CreateBinObjArray("psi","relntrcorr","");
  TIter nextBinRel(binsRel);
  AliAnalysisMuMuBinning::Range* rRel;
  
  while ( ( rRel = static_cast<AliAnalysisMuMuBinning::Range*>(nextBinRel()) ) )
  {
    CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,Form("EventsIn%s",rRel->AsString().Data()),Form("Number of %s events in %s bin",triggerClassName,rRel->AsString().Data()),1,0.,2);
  }
  
  delete binsRel;
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuNch::DefineSPDAcceptance()
{
  // Defines the functions ( eta = f(z) ) of the edges (right/left) of the inner and outer SPD layers
  // R_in = 3.9 cm ; R_out = 7.6 cm
  
  fSPD1LR = new TF1("fSPD1LR",SPDgeomR,-40,40,1);
  fSPD1LR->SetParameter(0,3.9);
  fSPD1LL = new TF1("fSPD1LL",SPDgeomL,-40,40,1);
  fSPD1LL->SetParameter(0,3.9);
  fSPD2LR = new TF1("fSPD2LR",SPDgeomR,-40,40,1);
  fSPD2LR->SetParameter(0,7.6);
  fSPD2LL = new TF1("fSPD2LL",SPDgeomL,-40,40,1);
  fSPD2LL->SetParameter(0,7.6);
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuNch::DefineSPDCorrectionMap(TObjArray* spdCorrectionList)
{
  // Defines a TMap of the SPD periods and the corresponding correction
  // Each SPD correction in the list must contain in the name the 1st run for which is valid and respect the naming convention:
  // (SPDCorrection_1stValidRunNumber_lastValidRunNumber).
  // The list must be ordered from smaller to bigger run number
  // This is usable with both types of corrections
  
  if ( !fSPDCorrectionMap ) // Creates a TMap for the run->Correction correspondance
  {
    fSPDCorrectionMap = new TMap;
    fSPDCorrectionMap->SetOwnerKeyValue(kTRUE,kTRUE);
  }
  
  TIter next(spdCorrectionList);
  TH1* SPDCorrection(0x0);
  Int_t runRef(0);
  Int_t i(0); // SPDCorrection index in the list
  while ( (SPDCorrection = static_cast<TH1*>(next())) ) // Checks if the correction list format is ok and if it is ordered
  {
    if ( static_cast<TH1*>(SPDCorrection)->IsA() != TH2::Class() && static_cast<TH1*>(SPDCorrection)->IsA() != TProfile::Class() )
    {
      AliFatal("Unrecognized SPD correction Class");
    }
      
    TString name = SPDCorrection->GetName();
    if ( !name.BeginsWith("SPDCorrection_") ) AliFatal(Form("Incorrect SPD correction at %d format: Objects in list must named as 'SPDCorrection_1stValidRunNumber_lastValidRunNumber'",i));
    
    name.Remove(0,name.First("_") + 1);
    TString nameLast = name;
    name.Remove(name.First("_"),name.Length());
    nameLast.Remove(0,nameLast.First("_") + 1);
    if ( !name.IsDigit() ) AliFatal(Form("Incorrect SPD correction at %d format: Impossible to retrieve first valid run number",i));
    if ( !nameLast.IsDigit() ) AliFatal(Form("Incorrect SPD correction at %d format: Impossible to retrieve last valid run number",i));
    
    Int_t runLow = name.Atoi();
    Int_t runHigh = nameLast.Atoi();
    if ( runHigh < runLow ) AliFatal(Form("SPD correction at %d validity range not valid",i));
    if ( runLow <= runRef ) AliFatal("SPD correction list not in correct order. Should be ordered from smaller to bigger run number");
    else runRef = runHigh;
    
    for ( Int_t j = runLow ; j <= runHigh ; j++ )
    {
      fSPDCorrectionMap->Add(new TObjString(Form("%d",j)),new TObjString(Form("%d",i))); // The mapping is done between the runs and the index of the correction in the spdCorrection List
    }
    
    i++;
  }
  
  fSPDCorrectionList = static_cast<TObjArray*>(spdCorrectionList->Clone());
}

//_____________________________________________________________________________
void AliAnalysisMuMuNch::AddHisto(const char* eventSelection,
                                  const char* triggerClassName,
                                  const char* centrality,
                                  const char* histoname,
                                  Double_t z,
                                  TH1* h,Bool_t isMC)
{
  // Adds the content of a 1D histo to a 2D histo a the z position
  
  Int_t zbin = fZAxis->FindBin(z);
  TH2F* h2;
  if (isMC) h2 = static_cast<TH2F*>(MCHisto(eventSelection,triggerClassName,centrality,histoname));
  else h2 = static_cast<TH2F*>(Histo(eventSelection,triggerClassName,centrality,histoname));
  
  for ( Int_t i = 1; i <= h->GetXaxis()->GetNbins(); ++i )
  {
    Double_t content = h2->GetCellContent(zbin,i);
    
    if ( h->GetBinContent(i) >  0 )
    {
      h2->SetCellContent(zbin,i,content + h->GetBinContent(i));
    }
  }
  
  h2->SetEntries(h2->GetSumOfWeights());
}


//_____________________________________________________________________________
void AliAnalysisMuMuNch::AttachSPDAcceptance(UInt_t dataType,
                                             const char* eventSelection,
                                             const char* triggerClassName,
                                             const char* centrality,const char* histoname)
{
  // Attachs AccxEff curves to the histogram
  
  if ( dataType & kHistoForData )
  {
    if( !Histo(eventSelection,triggerClassName,centrality,histoname) )
    {
      AliError(Form("ERROR: SPD Acceptance curves attach failed. Histo /%s/%s/%s/%s not found",eventSelection,triggerClassName,centrality,histoname));
      return;
    }
    
    Histo(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD1LR->Clone());
    Histo(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD1LL->Clone());
    Histo(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD2LR->Clone());
    Histo(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD2LL->Clone());
  }
  if ( (dataType & kHistoForMCInput) && HasMC() )
  {
    if( !MCHisto(eventSelection,triggerClassName,centrality,histoname) )
    {
      AliError(Form("ERROR: SPD Acceptance curves attach failed. MC Histo /%s/%s/%s/%s not found",eventSelection,triggerClassName,centrality,histoname));
      return;
    }
    
    MCHisto(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD1LR->Clone());
    MCHisto(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD1LL->Clone());
    MCHisto(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD2LR->Clone());
    MCHisto(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD2LL->Clone());
    
  }
  
}


//_____________________________________________________________________________
void AliAnalysisMuMuNch::FillHistosForEvent(const char* eventSelection,
                                            const char* triggerClassName,
                                            const char* centrality)
{
  // Fills the data (or reco if simu) multiplicity histos
  
  if ( IsHistogrammingDisabled() ) return;
  
  if ( fResolution ) return; //When computing resolutions we skip this method

//  if ( !AliAnalysisMuonUtility::IsAODEvent(Event()) )
//  {
//    AliError("Don't know how to deal with ESDs...");
//    return;
//  }
  
//  if ( HasMC() && !fSPDOneOverAccxEff && !fSPDMeanTracklets ) // We have MC but no correction (SPD correction computation mode)so we skip the method
//  {
//    return;
//  }

//  AliAODEvent* aod = static_cast<AliAODEvent*>(Event());
//
//  AliVVertex* vertex = aod->GetPrimaryVertexSPD();
  
  AliVVertex* vertex = AliAnalysisMuonUtility::GetVertexSPD(Event());
  
  TList* nchList = static_cast<TList*>(Event()->FindListObject("NCH"));
  
  if (!nchList || nchList->IsEmpty() ) // Empty NCH means that there is no SPD vertex ( see SetEvent() ) when runing on data.
  {
    AliError("Empty Nch list in event");
    return;
  }
//  nchList->Print();
  Double_t SPDZv;
  
  if ( !vertex || vertex->GetZ() == 0.0 ) // Running in Simu the spdZ == 0 means no SPD info. In data, avoid breaks in events w/o SPD vertex
  {
    SPDZv = -40.;
  } 
  else SPDZv = vertex->GetZ();
  
  TH1* hSPDcorrectionVsEta = static_cast<TH1*>(nchList->FindObject("SPDcorrectionVsEta"));
  TH1* hNTrackletVsEta = static_cast<TH1*>(nchList->FindObject("NTrackletVsEta"));
  TH1* hNTrackletVsPhi = static_cast<TH1*>(nchList->FindObject("NTrackletVsPhi"));
  
  TH1* hNTrackletSecVsEta(0x0);
  if ( fSPDMeanTrackletsCorrToCompare ) // In case we have a secondary eta range we clone the tracklets vs eta histo
  {
    hNTrackletSecVsEta = static_cast<TH1*>(hNTrackletVsEta->Clone());
    hNTrackletSecVsEta->Reset();
  }
  
  TH1* hNchVsEta = static_cast<TH1*>(hNTrackletVsEta->Clone("NchVsEta"));
  
  TProfile* hMeanTrackletsVsEta = static_cast<TProfile*>(Histo(eventSelection,triggerClassName,centrality,"MeanTrackletsVsEta"));
  TProfile* hMeanNchVsEta = static_cast<TProfile*>(Histo(eventSelection,triggerClassName,centrality,"MeanNchVsEta"));
  TProfile* hMeandNchdEtaVsEta = static_cast<TProfile*>(Histo(eventSelection,triggerClassName,centrality,"MeandNchdEtaVsEta"));

  TH2* hEventsVsZVertexVsEta = static_cast<TH2*>(Histo(eventSelection,triggerClassName,centrality,"EventsVsZVertexVsEta"));
    
  Int_t nBins(0);
  
  Int_t nTracklets[2] = {0,0}; // {nTracklets_Eta_fetaRange,nTracklets_Eta_fetaRangeComp} The first component is the number of tracklets in the primary eta range and the second is  the number of tracklets in the secondary eta range
  Double_t nch[2] = {0.,0.}; // {nCorrTracklets_Eta_fetaRange,nCorrTracklets_Eta_fetaRangeComp} The first component is the corrected(with the primary correction) number of tracklets in the primary eta range and the second is  the number of corrected(with the secondary correction) tracklets in the secondary eta range
  
  Int_t etaBinMin[2] = {fEtaAxis->FindBin(fetaRange[0]),fEtaAxis->FindBin(fEtaMinToCompare)}; // Primary and secondary (if any) eta ranges 
  Int_t etaBinMax[2] = {fEtaAxis->FindBin(fetaRange[1]),fEtaAxis->FindBin(fEtaMaxToCompare)-1}; // Note that the binMax is the bin-1 because the upper extreme of a bin belongs to the next bin (i.e. the bin containing -0.5 has as center -0.45, the bin is [-0.5,-0.4) but the bin containing 0.5 has as center 1.05, the bin is [0.5,1.1) so we have to take the previous bin [0.4,0.5)) This is already taken into account in the GetEtaRangeSPD method
  
  for (Int_t j = etaBinMin[0] ; j <= etaBinMax[0] ; j++) // Loop over eta bins.
  {
    Double_t correction = hSPDcorrectionVsEta->GetBinContent(j);
    
    Double_t eta = fEtaAxis->GetBinCenter(j);
    
    if ( correction < 0 ) continue;
    else if ( correction == 0.0 ) // No tracklets found in this eta bin.
    {
      correction = GetSPDCorrection(SPDZv,eta); //FIXME: Here eta is the bincenter not the exact one, will be a problem if the SPD AccxEff binning is not the same as the fEtaAxis one

      if ( correction == 0. || correction > 2.5) continue; // We need to know if the eta bin is in a region we want to take into account to count or not the zero
    }
  
    Int_t ntr = hNTrackletVsEta->GetBinContent(j); // Tracklets in eta bin
    nTracklets[0] += ntr;
    if ( j >= etaBinMin[1] && j <= etaBinMax[1] )
    {
      nTracklets[1] += ntr; // number of tracklets in the secondary eta range
      hNTrackletSecVsEta->SetBinContent(j,ntr);
    }
    
    hMeanTrackletsVsEta->Fill(eta,ntr); // Fill the number of tracklets of each eta bin in the profile

    if ( fSPDOneOverAccxEff )
    {
      nch[0] += ntr * correction; // Number of charged particles (SPD AccxEff corrected tracklets)
      hMeanNchVsEta->Fill(eta,ntr*correction);
      hMeandNchdEtaVsEta->Fill(eta,ntr*correction/fEtaAxis->GetBinWidth(5));
      hNchVsEta->SetBinContent(j,hNchVsEta->GetBinContent(j) * correction ); // Fill the number of charged particles of each eta bin in the profile
    }

    ++nBins; // We sum up the number of bins entering in the computation

    hEventsVsZVertexVsEta->Fill(SPDZv,eta,1.0); // Fill 1 count each eta bin where the events contributes
  }
  
  AddHisto(eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsPhi",SPDZv,hNTrackletVsPhi);
  AddHisto(eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsEta",SPDZv,hNTrackletVsEta);

  AddHisto(eventSelection,triggerClassName,centrality,"NchVsZVertexVsEta",SPDZv,hNchVsEta);
  
  Histo(eventSelection,triggerClassName,centrality,"Tracklets")->Fill(nTracklets[0]);
  Histo(eventSelection,triggerClassName,centrality,"MeanTrackletsVsZVertex")->Fill(SPDZv,nTracklets[0]);
  
    
  if ( fSPDMeanTracklets && !fSPDOneOverAccxEff ) // Just correct by SPDMeanTracklets method if there is no SPD AccxEff correction
  {
    Double_t SPDr = GetTrackletsMeanCorrection(SPDZv,nTracklets[0]); // Get 'mean correction' for the zvtx
    
    if ( SPDr < -999.) nch[0] = -1;
    else nch[0] = nTracklets[0] + SPDr; // In case of 'mean correction' nch has not be filled in the eta bins loop
    
    if ( fSPDMeanTrackletsCorrToCompare ) // Comparison of corrected tracklets in the primary and secondary eta ranges
    {
      AddHisto(eventSelection,triggerClassName,centrality,"TrackletsSecVsZVertexVsEta",SPDZv,hNTrackletSecVsEta);
      
      Double_t SPDrEtaComp = GetTrackletsMeanCorrection(SPDZv,nTracklets[1],kTRUE); // Get secondary 'mean correction' for the zvtx
      
      if ( SPDrEtaComp < -999.) nch[1] = -1;
      else nch[1] = nTracklets[1] + SPDrEtaComp; // In case of 'mean correction' nch has not be filled in the eta bins loop
      
      Histo(eventSelection,triggerClassName,centrality,"CorrTrackletsEtaSecVsCorrTrackletsEtaPrim")->Fill(nch[0],nch[1]);
      Histo(eventSelection,triggerClassName,centrality,"MeanNchEtaSecVsZVertex")->Fill(SPDZv,nch[1]); // Control plot to check if the secondary correction is applied correctly
      Histo(eventSelection,triggerClassName,centrality,"MeanTrackletsEtaSecVsZVertex")->Fill(SPDZv,nTracklets[1]);
    }
    
  }
  else if ( fSPDMeanTracklets && fSPDOneOverAccxEff ) // Here we compare the Corrected tracklets with the two correction methods
  {
    Int_t i(-1);
    Bool_t parFound1(kFALSE),parFound2(kFALSE);
    Double_t NtrCorr(0.),dNchdeta(0.);
    while ( i < nchList->GetEntries() - 1 && !(parFound1 && parFound2) )
    {
      i++;
      while ( nchList->At(i)->IsA() != TParameter<Double_t>::Class() && i < nchList->GetEntries() - 1 ) // In case there is a diferent object, just to skip it
      {
        i++;
      }
      
      TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(nchList->At(i));
      
      if ( TString(p->GetName()).Contains("NtrCorr") )
      {
        parFound1 = kTRUE;
        NtrCorr = p->GetVal();
      }
      else if ( TString(p->GetName()).Contains("MeandNchdEta") )
      {
        parFound2 = kTRUE;
        dNchdeta = p->GetVal();
      }
    }
    
    if ( parFound1 && parFound2 )
    {
//      Double_t dNchdetaPubli = 17.35; //FIXME: hardcoded (pPb value)
      Double_t ctToNch = 1.11; //FIXME: hardcoded (value for Nch vs NtrCorr(eta<0.5) in pPb)
      
      Histo(eventSelection,triggerClassName,centrality,"dNchdetaComparison2Corrections")->Fill(dNchdeta,ctToNch*NtrCorr/(2*fEtaMax));
      Histo(eventSelection,triggerClassName,centrality,"CheckMeanNtrCorrVsZVertex")->Fill(SPDZv,NtrCorr);
      if ( dNchdeta !=0 )
      {
        Histo(eventSelection,triggerClassName,centrality,"DispersiondNchdetaComparison2Corrections")->Fill((dNchdeta - ctToNch*NtrCorr/(2*fEtaMax)) / dNchdeta);
      }
      
      Histo(eventSelection,triggerClassName,centrality,"CheckNtrCorr")->Fill(NtrCorr);
      
    }
  }
  
  
  Histo(eventSelection,triggerClassName,centrality,"TrackletsVsNch")->Fill(nch[0],nTracklets[0]);
  Histo(eventSelection,triggerClassName,centrality,"Nch")->Fill(nch[0]);
  Histo(eventSelection,triggerClassName,centrality,"MeanNchVsZVertex")->Fill(SPDZv,nch[0]);
  
  
  Double_t V0AMult = 0.;
  Double_t V0CMult = 0.;
  
//  AliVVZERO* vzero = aod->GetVZEROData();
  AliVVZERO* vzero = Event()->GetVZEROData();
  if (vzero)
  {
    Double_t multV0A = vzero->GetMTotV0A();
    V0AMult = AliESDUtils::GetCorrV0A(multV0A,SPDZv);
    Double_t multV0C = vzero->GetMTotV0C();
    V0CMult = AliESDUtils::GetCorrV0C(multV0C,SPDZv);
    
    Histo(eventSelection,triggerClassName,centrality,"V0AMultVsNch")->Fill(V0AMult,nch[0]);
    Histo(eventSelection,triggerClassName,centrality,"V0CMultVsNch")->Fill(V0CMult,nch[0]);

  }

  
  delete hNchVsEta;
  
  // Mean dNch/dEta computation. In case the correction is the 'mean' one, this is not really dNch/deta, we have to multiply it by a factor <Nch>/<Ntrkls_corr> (we get <Nch> from a MC)
  Double_t meandNchdEta(0.);
  
  if ( nBins >  0 )
  {
    if ( fSPDOneOverAccxEff ) meandNchdEta = nch[0] / (nBins*fEtaAxis->GetBinWidth(5)); // Divide by nBins to get the mean and by the binWidht to get the d/dEta
    else if ( fSPDMeanTracklets ) meandNchdEta = nch[0] / (2.*fEtaMax); //fEtaAxis->GetBinWidth(5);
  }
  
  Histo(eventSelection,triggerClassName,centrality,"dNchdEta")->Fill(meandNchdEta);
  Histo(eventSelection,triggerClassName,centrality,"MeandNchdEtaVsZVertex")->Fill(SPDZv,meandNchdEta);
  
  
  //_____________These were tests //FIXME: Check if this tests are still neccesary_____________
  TObjArray* bins = Binning()->CreateBinObjArray("psi","dnchdeta","");
  TIter nextBin(bins);
  AliAnalysisMuMuBinning::Range* r;
  
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
  {
    Bool_t ok(kFALSE);
    if ( r->Quantity() == "DNCHDETA" )
    {
      ok = r->IsInRange(meandNchdEta);
    }
    
    if ( ok )
    {
      Histo(eventSelection,triggerClassName,centrality,Form("EventsIn%s",r->AsString().Data()))->Fill(1.);
    }
  }
  
  delete bins;
  
  
  TObjArray* binsNtrRaw = Binning()->CreateBinObjArray("psi","ntrraw","");
  TIter nextBinNtrRaw(binsNtrRaw);
  AliAnalysisMuMuBinning::Range* rNtrRaw;
  
  while ( ( rNtrRaw = static_cast<AliAnalysisMuMuBinning::Range*>(nextBinNtrRaw()) ) )
  {
    Bool_t ok(kFALSE);
    if ( rNtrRaw->Quantity() == "NTRRAW" )
    {
      ok = rNtrRaw->IsInRange(nTracklets[0]);
    }
    
    if ( ok )
    {
      Histo(eventSelection,triggerClassName,centrality,Form("MeanTrackletsVsZVertex%s",rNtrRaw->AsString().Data()))->Fill(SPDZv,nTracklets[0]);
    }
  }

  delete binsNtrRaw;
  
  
  TObjArray* binsNtr = Binning()->CreateBinObjArray("psi","ntrcorr","");
  TIter nextBinNtr(binsNtr);
  AliAnalysisMuMuBinning::Range* rNtr;

  while ( ( rNtr = static_cast<AliAnalysisMuMuBinning::Range*>(nextBinNtr()) ) )
  {
    Bool_t ok(kFALSE);
    if ( rNtr->Quantity() == "NTRCORR" )
    {
      ok = rNtr->IsInRange(nch[0]);
    }
    
    if ( ok )
    {
      Histo(eventSelection,triggerClassName,centrality,Form("MeanNchVsZVertex%s",rNtr->AsString().Data()))->Fill(SPDZv,nch[0]);
    }
  }
  
  delete binsNtr;
  
  
  TObjArray* binsRel = Binning()->CreateBinObjArray("psi","relntrcorr","");
  TIter nextBinRel(binsRel);
  AliAnalysisMuMuBinning::Range* rRel;
  while ( ( rRel = static_cast<AliAnalysisMuMuBinning::Range*>(nextBinRel()) ) )
  {
    Bool_t ok(kFALSE);
    if ( rRel->Quantity() == "RELNTRCORR" )
    {
      ok = rRel->IsInRange(nch[0]/fMeanTrRef);
    }
    
    if ( ok )
    {
      Histo(eventSelection,triggerClassName,centrality,Form("EventsIn%s",rRel->AsString().Data()))->Fill(1.);
    }
  }
  
  delete binsRel;
  //_________________________________________________________________________________________
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuNch::FillHistosForMCEvent(const char* eventSelection,const char* triggerClassName,const char* centrality)
{
  /// Fill input MC multiplicity histos
  
  if ( IsHistogrammingDisabled() ) return;
  
  TList* nchList = static_cast<TList*>(Event()->FindListObject("NCH"));
  
  if (!nchList || nchList->IsEmpty())
  {
    AliError("Empty Nch list in event");
    return;
  }
  
  Double_t MCZv = AliAnalysisMuonUtility::GetMCVertexZ(Event(),MCEvent()); // Definition of MC generated z vertex
//  AliVVertex* vertex = static_cast<AliAODEvent*>(Event())->GetPrimaryVertexSPD();
  AliVVertex* vertex = AliAnalysisMuonUtility::GetVertexSPD(Event());
  
  Double_t SPDZv(0.);
  TParameter<Double_t>* p(0x0);
  
  //____Resolution Histos___
  if ( !fSPDOneOverAccxEff && !fSPDMeanTracklets && fResolution )
  {
    Int_t nContributors(0);
    if ( vertex )
    {
      nContributors  = vertex->GetNContributors();
      SPDZv = vertex->GetZ();
    }
    MCHisto(eventSelection,triggerClassName,centrality,"SPDZvResVsnC")->Fill(nContributors,SPDZv - MCZv);
    MCHisto(eventSelection,triggerClassName,centrality,"SPDZvResVsMCz")->Fill(MCZv,SPDZv - MCZv);
    
    Double_t EtaReco(0.),EtaMC(0.),PhiReco(0.),PhiMC(0.);
    Int_t i(-1),labelEtaReco(-1),labelEtaMC(-1),labelPhiReco(-1),labelPhiMC(-1);
    
    while ( i < nchList->GetEntries() - 1 )
    {
      i++;
      while ( nchList->At(i)->IsA() != TParameter<Double_t>::Class() ) // In case there is a diferent object, just to skip it
      {
        i++;
      }
      
      p = static_cast<TParameter<Double_t>*>(nchList->At(i));

      
         //Eta Resolution
      if ( TString(p->GetName()).Contains("EtaReco") ) // We take the reco eta
      {
        sscanf(p->GetName(),"EtaReco%d",&labelEtaReco);
        EtaReco = p->GetVal();
        MCHisto(eventSelection,triggerClassName,centrality,"Eta")->Fill(EtaReco);
      }
      else if ( TString(p->GetName()).Contains("EtaMC") ) // We take the generated eta
      {
        sscanf(p->GetName(),"EtaMC%d",&labelEtaMC);
        EtaMC = p->GetVal();
        MCHisto(eventSelection,triggerClassName,centrality,"MCEta")->Fill(EtaMC);
      }
      if ( labelEtaReco > 0 && labelEtaReco == labelEtaMC ) // To be sure we compute the difference for the same particle
      {
        labelEtaReco = -1; // Restart of the label value to avoid double count the eta difference when computing the phi one 
        Double_t EtaDif = EtaReco - EtaMC;
        MCHisto(eventSelection,triggerClassName,centrality,"EtaRes")->Fill(EtaDif);
        MCHisto(eventSelection,triggerClassName,centrality,"EtaResVsZ")->Fill(MCZv,EtaDif);
        MCHisto(eventSelection,triggerClassName,centrality,"EtaResVsnC")->Fill(nContributors,EtaDif);
      }
      
         //Phi Resolution
      else if ( TString(p->GetName()).Contains("PhiReco") ) // We take the reco phi
      {
        sscanf(p->GetName(),"PhiReco%d",&labelPhiReco);
        PhiReco = p->GetVal();
        MCHisto(eventSelection,triggerClassName,centrality,"Phi")->Fill(PhiReco);
      }
      else if ( TString(p->GetName()).Contains("PhiMC") ) // We take the generated phi
      {
        sscanf(p->GetName(),"PhiMC%d",&labelPhiMC);
        PhiMC = p->GetVal();
        MCHisto(eventSelection,triggerClassName,centrality,"MCPhi")->Fill(PhiMC);
      }
      
      if ( labelPhiReco > 0 && labelPhiReco == labelPhiMC ) // To be sure we compute the difference for the same particle
      {
        labelPhiReco = -1; // Restart of the label value to avoid double count the phi difference when computing the eta one
        Double_t PhiDif = PhiReco - PhiMC;
        MCHisto(eventSelection,triggerClassName,centrality,"PhiRes")->Fill(PhiDif);
        
        //___With the following algorithm we refer the differences to the interval [-Pi/2,Pi/2]
        if ( PhiDif < -TMath::PiOver2() && PhiDif > -TMath::Pi() )
        {
          PhiDif = -TMath::Pi() - PhiDif;
        }
        else if ( PhiDif < -TMath::Pi() )
        {
          PhiDif = -2.*TMath::Pi() - PhiDif;
          if ( PhiDif < -TMath::PiOver2() )
          {
            PhiDif = -TMath::Pi() - PhiDif;
          }
        }
        
        else if ( PhiDif > TMath::PiOver2() && PhiDif < TMath::Pi() )
        {
          PhiDif = TMath::Pi() - PhiDif;
        }  
        else if ( PhiDif > TMath::Pi() )
        {
          PhiDif = 2.*TMath::Pi() - PhiDif;
          if ( PhiDif > TMath::PiOver2() )
          {
            PhiDif = TMath::Pi() - PhiDif;
          }
        }
        //___
        
        MCHisto(eventSelection,triggerClassName,centrality,"PhiResVsZ")->Fill(MCZv,PhiDif);
        MCHisto(eventSelection,triggerClassName,centrality,"PhiResShifted")->Fill(PhiDif);
        MCHisto(eventSelection,triggerClassName,centrality,"PhiResVsnC")->Fill(nContributors,PhiDif);
      }
    }
    
    return; // When computing resolutions we skip the rest of the method
  }
  //_____
  
  //___Multiplicity histos (just when not computing resolutions)___
  TH1* hSPDcorrectionVsEta = static_cast<TH1*>(nchList->FindObject("MCSPDcorrectionVsEta"));
  TH1* hNBkgTrackletsVSEta = static_cast<TH1*>(nchList->FindObject("NBkgTrackletsVSEta"));
  TH1* hNchVsPhi = static_cast<TH1*>(nchList->FindObject("MCNchVsPhi"));
  TH1* hNchVsEta = static_cast<TH1*>(nchList->FindObject("MCNchVsEta"));
  TH1* hNTrackletVsEta = static_cast<TH1*>(nchList->FindObject("MCNTrackletVsEta"));
  TH1* hNTrackletVsPhi = static_cast<TH1*>(nchList->FindObject("MCNTrackletVsPhi"));
  
  TProfile* hMeanNchVsEta = static_cast<TProfile*>(MCHisto(eventSelection,triggerClassName,centrality,"MeanNchVsEta"));
  
  TH2* hEventsVsZVertexVsEta = static_cast<TH2*>(MCHisto(eventSelection,triggerClassName,centrality,"EventsVsZVertexVsEta"));

  Int_t nBins(0);
  
  Double_t nchSum(0.),nTracklets(0.);
  for (Int_t j = 1 ; j <= fEtaAxis->GetNbins() ; j++) //Loop over all eta bins 
  {
    Double_t correction = hSPDcorrectionVsEta->GetBinContent(j);
    
    if ( correction < 0 ) continue; // We count just the particles in the SPD acceptance. (MCSPDcorrectionVsEta is filled to -1 in SetEvent() for those bins out of range)
    
    nTracklets += hNTrackletVsEta->GetBinContent(j); // Reco tracklets
    
    ++nBins; // We sum up the number of bins entering in the computation
    
    Double_t eta = fEtaAxis->GetBinCenter(j);
    Double_t nch = hNchVsEta->GetBinContent(j); // Generated particles
    
    nchSum += nch;
    
    hMeanNchVsEta->Fill(eta,nch); // Fill the number of charged particles of each eta bin in the profile
    
    hEventsVsZVertexVsEta->Fill(MCZv,eta,1.0); // Fill 1 count each eta bin where the events contributes
  }
  
  MCHisto(eventSelection,triggerClassName,centrality,"Tracklets")->Fill(nTracklets); // Note that these are NOT the same tracklets as in the FillHistosForEvent() since here the SPD "dead" zones (the ones where correction > threshold) are not rejected
  
  if ( vertex )
  {
    SPDZv = vertex->GetZ();
  }

  Bool_t isMChisto = kTRUE; // Used to get the MC histos in the Add method
  
  AddHisto(eventSelection,triggerClassName,centrality,"NBkgTrackletsVsZVertexVsEta",SPDZv,hNBkgTrackletsVSEta,isMChisto);
  AddHisto(eventSelection,triggerClassName,centrality,"NchVsZVertexVsPhi",MCZv,hNchVsPhi,isMChisto);
  AddHisto(eventSelection,triggerClassName,centrality,"NchVsZVertexVsEta",MCZv,hNchVsEta,isMChisto);
  AddHisto(eventSelection,triggerClassName,centrality,"NchVsRecoZVertexVsEta",SPDZv,hNchVsEta,isMChisto);
  AddHisto(eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsEta",SPDZv,hNTrackletVsEta,isMChisto);
  AddHisto(eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsPhi",SPDZv,hNTrackletVsPhi,isMChisto);
  
  MCHisto(eventSelection,triggerClassName,centrality,"MeanNchVsZVertex")->Fill(MCZv,nchSum);
  MCHisto(eventSelection,triggerClassName,centrality,"Nch")->Fill(nchSum);
  
  
  // Mean dNch/dEta computation
  Double_t meandNchdEta(0.);
  if ( nBins >  0 )
  {
    meandNchdEta = nchSum / (nBins*fEtaAxis->GetBinWidth(5)); // Divide by nBins to get the mean and by the binWidht to get the d/dEta
  }
  
  MCHisto(eventSelection,triggerClassName,centrality,"dNchdEta")->Fill(meandNchdEta);
  
  
  Int_t i(-1);
  p = 0x0;
  Double_t nTrCorr(0.);
  Double_t dNchdetaReco(0.);
  Double_t ctToNch = 1.11; //FIXME: hardcoded (value for Nch vs NtrCorr(eta<0.5) in pPb)
  while ( i < nchList->GetEntries() - 1 )
  {
    i++;
    while ( nchList->At(i)->IsA() != TParameter<Double_t>::Class() && i < nchList->GetEntries() - 1 ) // In case there is a diferent object, just to skip it
    {
      i++;
    }
    
    p = static_cast<TParameter<Double_t>*>(nchList->At(i));
    
    if ( ( TString(p->GetName()).Contains("NtrCorr") || TString(p->GetName()).BeginsWith("Nch") ) && ( (fSPDMeanTracklets && !fSPDOneOverAccxEff) || (!fSPDMeanTracklets && fSPDOneOverAccxEff) ))
    {
      MCHisto(eventSelection,triggerClassName,centrality,"CorrTrackletsVsNch")->Fill(nchSum,p->GetVal());
    }
    else if ( TString(p->GetName()).Contains("NtrCorr") )
    {
      nTrCorr = p->GetVal();
      MCHisto(eventSelection,triggerClassName,centrality,"dNchdetaFromNtrCorrVsdNchdEtaMC")->Fill(meandNchdEta,ctToNch*nTrCorr/(2*fEtaMax));
      MCHisto(eventSelection,triggerClassName,centrality,"RelDispersiondNchdetaFromNtrCorrVsdNchdEtaMC")->Fill((meandNchdEta - ctToNch*nTrCorr/(2*fEtaMax)) / meandNchdEta);
      MCHisto(eventSelection,triggerClassName,centrality,"DispersiondNchdetaFromNtrCorrVsdNchdEtaMC")->Fill((meandNchdEta - ctToNch*nTrCorr/(2*fEtaMax)));
    }
    else if ( TString(p->GetName()).Contains("Ntr") && !TString(p->GetName()).Contains("SPDOk") && !TString(p->GetName()).Contains("Corr"))
    {
      MCHisto(eventSelection,triggerClassName,centrality,"TrackletsVsNch")->Fill(nchSum,p->GetVal());
    }
    else if ( TString(p->GetName()).Contains("MeandNchdEta") )
    {
      dNchdetaReco = p->GetVal();
      MCHisto(eventSelection,triggerClassName,centrality,"dNchdetaVsMCdNchdeta")->Fill(meandNchdEta,dNchdetaReco);
      MCHisto(eventSelection,triggerClassName,centrality,"dNchdetaFromAccEffVsdNchdEtaMC")->Fill(meandNchdEta,dNchdetaReco);
      MCHisto(eventSelection,triggerClassName,centrality,"RelDispersiondNchdetaFromAccEffVsdNchdEtaMC")->Fill((meandNchdEta - dNchdetaReco) / meandNchdEta);
      MCHisto(eventSelection,triggerClassName,centrality,"DispersiondNchdetaFromAccEffVsdNchdEtaMC")->Fill((meandNchdEta - dNchdetaReco));
    }
  }
  
  
  Double_t V0AMult = 0.;
  Double_t V0CMult = 0.;
  
  AliVVZERO* vzero = Event()->GetVZEROData();
  if (vzero)
  {
    Double_t multV0A = vzero->GetMTotV0A();
    V0AMult = AliESDUtils::GetCorrV0A(multV0A,MCZv);
    Double_t multV0C = vzero->GetMTotV0C();
    V0CMult = AliESDUtils::GetCorrV0C(multV0C,MCZv);
    
    MCHisto(eventSelection,triggerClassName,centrality,"V0AMultVsNch")->Fill(V0AMult,nchSum);
    MCHisto(eventSelection,triggerClassName,centrality,"V0CMultVsNch")->Fill(V0CMult,nchSum);
    
  }
  //____
}


//_____________________________________________________________________________
Bool_t AliAnalysisMuMuNch::GetEtaRangeSPD(Double_t spdZVertex, Double_t etaRange[])
{
  // Fill the SPD eta range for a given z vertex and returns whether the range is valid or not.
  
  if ( fSPDMeanTracklets || (!fSPDMeanTracklets && !fSPDOneOverAccxEff) )
  {
    // If we are using the mean correction or we want the raw distributions we should not apply the algorithm to restric the eta range, in order to avoid the effect of the eta-Zv bins on the distributions ( the algorithm makes that when an eta bin is partially out of the SPD acceptance we remove it. This makes that the mean raw tracklets distribution is not smooth for Zv in which the eta range is partially out of the SPD )
    
    //FIXME: Implement something to warning or reject events for which the selected eta range is completely out of the SPD acceptance
    etaRange[0] = fEtaMin;
    etaRange[1] = fEtaMax - 1E-6;
    
    return kTRUE;
  }
  
  Double_t etaMax(fEtaMax),etaMin(fEtaMin);
  
  Double_t vf1LR = fSPD1LR->Eval(spdZVertex); //Eta values for the z vertex over the SPD acceptance curves
  Double_t vf2LR = fSPD2LR->Eval(spdZVertex);
  Double_t vf1LL = fSPD1LL->Eval(spdZVertex);
  Double_t vf2LL = fSPD2LL->Eval(spdZVertex);
  
  //____We start by asigning the eta range as the eta values in the SPD acceptance curve
  if ( spdZVertex < 0. )
  {
    etaMax = vf2LR;
    etaMin = TMath::Max(vf1LL,vf2LL);
  }
  else
  {
    etaMax = TMath::Min(vf1LR,vf2LR);
    etaMin = vf2LL;
  }
  //____
  
  
  //____Algorithm to avoid taking bins which are crossed by the SPD acceptance curves
  Int_t binYMin = fEtaAxis->FindBin(etaMin); // Find the corresponding bins for eta max & min and z vertex
  Int_t binYMax = fEtaAxis->FindBin(etaMax);
  Int_t binX = fZAxis->FindBin(spdZVertex);
  
  // Define the values for the relevant edges of the eta and z bins
  Double_t upEdge = fEtaAxis->GetBinUpEdge(binYMax);  // up edge of the top eta bin
  Double_t lowEdge = fEtaAxis->GetBinLowEdge(binYMin); // low edge of the bottom eta bin
  Double_t leftEdge = fZAxis->GetBinLowEdge(binX); // left edge of the z bin
  Double_t rightEdge = fZAxis->GetBinUpEdge(binX); // right edge of the z bin
  
  Double_t etaMaxTemp(0.),etaMinTemp(0.);
  if ( spdZVertex < 0. )
  {
    etaMaxTemp = fSPD2LR->Eval(rightEdge); // Define the temporary eta max as the value of the curve int the righ edge of the bin
    etaMinTemp = TMath::Max(fSPD1LL->Eval(leftEdge),fSPD2LL->Eval(leftEdge));
  }
  else
  {
    etaMaxTemp = TMath::Min(fSPD1LR->Eval(rightEdge),fSPD2LR->Eval(rightEdge));
    etaMinTemp = fSPD2LL->Eval(leftEdge);
  }
  
  while ( upEdge > etaMaxTemp ) //We take eta max as the up edge of the 1st bin which is inside the SPD acceptance but not crossed by the curve
  {
    binYMax = binYMax - 1; // Since the up edge of the current bin is bigger than the SPD acc curve we move 1 bin below
    upEdge = fEtaAxis->GetBinUpEdge(binYMax); // Take the up edge of the new bin
    etaMax = upEdge - 1E-6; // We substract 1E-6 cause the up edge of the a bin belongs to the bin+1
    
  }
  
  //We take eta min as the low edge of the 1st bin which is inside the SPD acceptance but not crossed by the curve
  while ( lowEdge < etaMinTemp )
  {
    binYMin = binYMin + 1; // Since the low edge of the current bin is smaller than the SPD acc curve we move 1 bin above
    lowEdge = fEtaAxis->GetBinLowEdge(binYMin); // Take the low edge of the new bin
    etaMin = lowEdge + 1E-6;
  }
  //____
  
  // In case the eta range we want (given in the constructor) is smaller than the one found by the algorithm (max range) we redefine the values
  if ( etaMin < fEtaMin ) etaMin = fEtaMin + 1E-6;
  if ( etaMax > fEtaMax ) etaMax = fEtaMax - 1E-6;
  
  if ( etaMin > etaMax ) // If this happens, the event is not able to cover the desired eta range, so we set an invalid range
  {
    etaRange[0] = 999.;
    etaRange[1] = -999.;
    
    return kFALSE; //FIXME: Shall we add something to invalidate this event in SetEvent to dont count it in AliAnalysisTaskMuMu::FillCounters?
  }
  else
  {
    etaRange[0] = etaMin;
    etaRange[1] = etaMax;
    
    return kTRUE;
  }
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuNch::GetSPDCorrection(Double_t zvert, Double_t eta) const
{
  if ( !fSPDOneOverAccxEff )
  {
    return 1.; // Like this we will simply compute the raw tracklets or if we are correctiong by the mean tracklets the eta bin will be counted as valid
  }
  else
  {
    Int_t bin =fSPDOneOverAccxEff->FindBin(zvert,eta);
    return fSPDOneOverAccxEff->GetBinContent(bin);// Like this we will compute Nch from the SPD AccxEFf
  }
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuNch::GetTrackletsMeanCorrection(Double_t zvert, Int_t nTracklets, Bool_t corrToCompare) const
{
  if ( !fSPDMeanTracklets )
  {
    AliFatal("ERROR: No tracklets mean correction");
    return 0;
  }
  else
  {
    TProfile* h(0x0);
    if ( corrToCompare && fSPDMeanTrackletsCorrToCompare ) h = static_cast<TProfile*>(fSPDMeanTrackletsCorrToCompare);
    else h = static_cast<TProfile*>(fSPDMeanTracklets);

    Int_t bin = h->FindBin(zvert);
    Double_t mNTtklsZ = h->GetBinContent(bin);
    
    Double_t mNTtklsZref(0.);
    if ( fSPDMeanTrackletsCorrToCompare )
    {
      if ( corrToCompare ) mNTtklsZref = h->GetMaximum();
      else
      {
        if ( fMeanTrRef > 0. ) mNTtklsZref = fMeanTrRef;
        else mNTtklsZref = h->GetMaximum();
      }
    }
    else
    {
      if ( fMeanTrRef > 0. ) mNTtklsZref = fMeanTrRef;
      else mNTtklsZref = h->GetMaximum();
    }

    Double_t deltaN;
    if ( mNTtklsZ == 0. ) deltaN = -1000.; // If the zvertex is out of the correction range the correction is < -999 (non valid)
    else
    {
      Double_t coef(1.);
      if ( mNTtklsZref < mNTtklsZ ) coef = -1.;
      
      Double_t meanPoiss = nTracklets*(mNTtklsZref - mNTtklsZ)/mNTtklsZ;
      deltaN = coef*frand->PoissonD(TMath::Abs(meanPoiss));
    }
    return deltaN;
  }
}


//_____________________________________________________________________________
AliAODTracklets* AliAnalysisMuMuNch::GetTracklets(AliVEvent* event)
{
  /// Return AliAODTracklets if the event is an AOD event. If event is ESD event creates AliAODTracklets object from ESD AliMultiplicity and returns it
  
  AliAODTracklets* tracklets(0x0);
  
  if ( event->IsA() == AliAODEvent::Class() )
  {
    tracklets = static_cast<const AliAODEvent*>(Event())->GetTracklets(); // Tracklets object
  }
  else if ( event->IsA() == AliESDEvent::Class() )
  {
    const AliMultiplicity* mult = static_cast<const AliESDEvent*>(Event())->GetMultiplicity();
    
    if (mult)
    {
      tracklets = new AliAODTracklets();
      
      if ( mult->GetNumberOfTracklets() > 0 )
      {
        tracklets->CreateContainer(mult->GetNumberOfTracklets());
        for (Int_t n = 0 ; n < mult->GetNumberOfTracklets() ; n++)
        {
          tracklets->SetTracklet(n, mult->GetTheta(n), mult->GetPhi(n), mult->GetDeltaPhi(n), mult->GetLabel(n, 0),mult->GetLabel(n, 1));
        }
      }
    }
  }
  else AliFatal("Unrecognized Event Type");
  
  return tracklets;

}

//==============================These methods are useless here, anyway they should go in the EventCutterClass=================
//_____________________________________________________________________________
Bool_t AliAnalysisMuMuNch::HasAtLeastNTrackletsInEtaRange(const AliVEvent& event, Int_t n, Double_t& etaMin, Double_t& etaMax) const
{
  if ( event.IsA() != AliAODEvent::Class() )
  {
    return kFALSE;
  }
  
  AliAODTracklets* tracklets = static_cast<const AliAODEvent&>(event).GetTracklets();
  
  if (!tracklets)
  {
    return kFALSE;
  }
  
  Int_t nTrackletsInRange(0);
  
  Int_t nTracklets = tracklets->GetNumberOfTracklets();
  
  for (Int_t i = 0 ; i < nTracklets && nTrackletsInRange < n; i++)
  {
    Double_t eta = -TMath::Log(TMath::Tan(tracklets->GetTheta(i)/2.0));
    
    if ( eta > etaMin && eta < etaMax )
    {
      ++nTrackletsInRange;
    }
  }
  
  return (nTrackletsInRange>=n);
}

//_____________________________________________________________________________
void AliAnalysisMuMuNch::NameOfHasAtLeastNTrackletsInEtaRange(TString& name, Int_t n, Double_t& etaMin, Double_t& etaMax) const
{
  if ( TMath::AreEqualAbs(TMath::Abs(etaMin),TMath::Abs(etaMax),1E-9 ) )
  {
    name.Form("ATLEAST%dTRKLINABSETALT%3.1f",n,TMath::Abs(etaMin));
  }
  else
  {
    name.Form("ATLEAST%dTRKLINETA%3.1f-%3.1f",n,etaMin,etaMax);
  }
}
//=======================================================================================

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuNch::IsMCtrackFromGenerator(Int_t indexMC) const
{
  ///Checks if the MC particle corresponding to a track has been generated with a given generator.
  
  if ( !HasMC() )
  {
    AliWarning("There is no MC information in the event");
    return kFALSE;
  }
  
  if ( indexMC >= MCEvent()->GetNumberOfTracks() )
  {
    AliWarning("This MC track does not exist. Index out of range");
//    return kFALSE;
  }

//______________"MANUAL" METHOD___________//
  //It does not work since NProduced gives just the stable dpmjet produced particles. See AliMCEvent::GetGenerator
//  Int_t labelTemp(labelMC);
//  while ( labelTemp != -1 )
//  {
//    AODpartMCTemp = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(labelMC));
//    labelTemp = AODpartMCTemp->GetMother();
//    if ( labelTemp != -1 ) labelMC = labelTemp;
//    std::cout << labelTemp << " ===>" << labelMC <<std::endl;
//  }
//  
//  AliAODMCHeader* mcHeader = static_cast<AliAODMCHeader*>(Event()->FindListObject(AliAODMCHeader::StdBranchName()));
//  if(mcHeader)
//  {
//    TList* lheaders = mcHeader->GetCocktailHeaders();
//    lheaders->Print();
//
//    AliGenEventHeader* mcGenH(0x0);
//    TIter next(lheaders); // Get the iterator on the list of cocktail headers
//    Int_t nProduced(0);
//    while ( labelMC >= nProduced ) // Loop over the cocktail headers
//    {
//      mcGenH = static_cast<AliGenEventHeader*>(next());
//      if ( !mcGenH )
//      {
//        std::cout << "Generator header not found" << std::endl;
//        return kFALSE;
//      }
//      else nProduced += mcGenH->NProduced();
//    }
//  }
//  TString genName(mcGenH->GetName());
//  return genName.Contains("dpmjet_0");
 //_________________________________________
  
  TString genName;
  if ( !MCEvent()->GetCocktailGenerator(indexMC,genName) )
  {
    AliWarning("No cocktail generator found for this event");
    return kFALSE;
  }
//  std::cout << genName.Data() << std::endl;
  return genName.Contains(fGeneratorHeaderClass->Data());
}

//==============================This method is not used=====================
//_____________________________________________________________________________
Double_t AliAnalysisMuMuNch::NumberOfTrackletsInEtaRange(const AliVEvent& event, Double_t& etaMin, Double_t& etaMax, Bool_t /*corrected*/) const
{
  if ( event.IsA() != AliAODEvent::Class() )
  {
    AliError("Not working for ESDs...");
    return 0;
  }
  
  AliAODTracklets* tracklets = static_cast<const AliAODEvent&>(event).GetTracklets();
  
  if (!tracklets)
  {
    return 0;
  }
  
  Double_t nTrackletsInRange(0);

  Int_t nTracklets = tracklets->GetNumberOfTracklets();
  
  for (Int_t i = 0 ; i < nTracklets ; i++)
  {
    Double_t thetaTracklet = tracklets->GetTheta(i);
    Double_t etaTracklet = -TMath::Log(TMath::Tan(thetaTracklet/2.));
    
    if ( etaTracklet > etaMin && etaTracklet < etaMax )
    {
      nTrackletsInRange += 1.0;
    }
  }
  
  return nTrackletsInRange;
}
//==============================


//_____________________________________________________________________________
void AliAnalysisMuMuNch::SetEvent(AliVEvent* event, AliMCEvent* mcEvent)
{
  /// Set the event, compute multiplicities and add them as a list to the event
  
//  if ( event->IsA() != AliAODEvent::Class() )
//  {
//    AliError("Only working for AODs for the moment.");
//    return;
//  }
  
  AliAnalysisMuMuBase::SetEvent(event,mcEvent); // To have Event() and MCEvent() method working
  
  TList* nchList = static_cast<TList*>(event->FindListObject("NCH")); // Define the list with the NCH info for each event
  if (!nchList)
  {
    nchList = new TList;
    nchList->SetOwner(kTRUE);
    nchList->SetName("NCH");
    event->AddObject(nchList);
  }
  
  nchList->Clear(); // We clear the NCH list for this new event
  
  Double_t SPDZv(0.);
  Bool_t vertexSPDFound(kFALSE);
//  AliAODTracklets* tracklets(0x0);
  
  AliVVertex* vertexSPD = AliAnalysisMuonUtility::GetVertexSPD(event);
  if ( vertexSPD )
  {
    vertexSPDFound = kTRUE;
    SPDZv = vertexSPD->GetZ();
    if ( SPDZv == 0. ) SPDZv = -40.; // If SPDZv = 0. means that there is no reconstructed vertex. Setting the vertex to -40 allow us to get an invalid correction
  }
  
  AliAODTracklets* tracklets = GetTracklets(event);

//  if ( event->IsA() == AliAODEvent::Class() )
//  {
//    const AliAODVertex* vertexSPD = static_cast<const AliAODEvent*>(Event())->GetPrimaryVertexSPD(); // SPD vertex object
//    if ( vertexSPD )
//    {
//      vertexSPDFound = kTRUE;
//      SPDZv = vertexSPD->GetZ();
//      if ( SPDZv == 0. ) SPDZv = -40.; // If SPDZv = 0. means that there is no reconstructed vertex. Setting the vertex to -40 allow us to get an invalid correction
//    }
////    else return;
//    
//    tracklets = static_cast<const AliAODEvent*>(Event())->GetTracklets(); // Tracklets object
//  }
//
//  else if ( event->IsA() == AliESDEvent::Class() )
//  {
//    const AliESDVertex* vertexSPD = static_cast<const AliESDEvent*>(Event())->GetPrimaryVertexSPD(); // SPD vertex object
//    if ( vertexSPD )
//    {
//      vertexSPDFound = kTRUE;
//      SPDZv = vertexSPD->GetZ();
//    }
////    else return;
//    
//    tracklets = new AliAODTracklets();
//    const AliMultiplicity* mult = static_cast<const AliESDEvent*>(Event())->GetMultiplicity();
//    
//    if (mult)
//    {
//      if (mult->GetNumberOfTracklets()>0)
//      {
//        tracklets->CreateContainer(mult->GetNumberOfTracklets());
//        for (Int_t n=0; n<mult->GetNumberOfTracklets(); n++)
//        {
////          if(fMChandler)
////          {
////            fMChandler->SelectParticle(mult->GetLabel(n, 0));
////            fMChandler->SelectParticle(mult->GetLabel(n, 1));
////          }
//          tracklets->SetTracklet(n, mult->GetTheta(n), mult->GetPhi(n), mult->GetDeltaPhi(n), mult->GetLabel(n, 0),mult->GetLabel(n, 1));
//        }
//      }
//    }
//  }
//  else AliFatal("Unrecognized Event Type");
  
  TH1* hSPDcorrectionVsEta(0x0); // Pointers for the individual event histos
  TH1* hNchVsEta(0x0);
  TH1* hNTrackletVsEta(0x0);
  TH1* hNTrackletVsPhi(0x0);
  
  // Variables we will use in the multiplicity computation
  Int_t binMin,binMax;
  Int_t nTracklets(0);
  Double_t thetaTracklet(0.),phiTracklet(0.),etaTracklet(0.);
  Double_t nch(0.0);
  Int_t nTrackletsInRange(0);
  Int_t nTrackletsInRangeSPDOk(0);
  Int_t nBins(0);

  
  Bool_t isEtaRangeValid = GetEtaRangeSPD(SPDZv,fetaRange);
  
  //____Data (or reco) multiplicity computation ___
  if ( !fResolution ) // When computing resolutions we dont do anything else
  {
    //_______Create once the histos with the individual event "properties" (cleared at the beginning of each new event)_______
    if ( !Histo("AliAnalysisMuMuNch","SPDcorrectionVsEta") )
    {
      CreateEventHistos(kHistoForData | kHistoForMCInput,"AliAnalysisMuMuNch","SPDcorrectionVsEta","SPD correction-like vs #eta;#eta",fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
      
      CreateEventHistos(kHistoForMCInput,"AliAnalysisMuMuNch","NchVsEta","Nch vs #eta;#eta",fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());

      CreateEventHistos(kHistoForData | kHistoForMCInput,"AliAnalysisMuMuNch","NTrackletVsEta","Ntracklet vs #eta;#eta",fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
      
      CreateEventHistos(kHistoForData,"AliAnalysisMuMuNch","NBkgTrackletsVSEta","NBkg Tracklets vs #eta;#eta",fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
      
      Double_t phimin = 0.; //Phi range
      Double_t phimax = 2*TMath::Pi();
      Int_t nphibins = GetNbins(phimin,phimax,0.05);
      
      CreateEventHistos(kHistoForData | kHistoForMCInput,"AliAnalysisMuMuNch","NTrackletVsPhi","Ntracklet vs #phi;#phi",nphibins,phimin,phimax);
      
      CreateEventHistos(kHistoForMCInput,"AliAnalysisMuMuNch","NchVsPhi","Nch vs #phi;#phi",nphibins,phimin,phimax);
      
      CreateEventHistos(kHistoForData | kHistoForMCInput,"AliAnalysisMuMuNch","test","test",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
      
      Histo("AliAnalysisMuMuNch","test")->GetListOfFunctions()->Add(fSPD1LR->Clone());
      Histo("AliAnalysisMuMuNch","test")->GetListOfFunctions()->Add(fSPD1LL->Clone());
      Histo("AliAnalysisMuMuNch","test")->GetListOfFunctions()->Add(fSPD2LR->Clone());
      Histo("AliAnalysisMuMuNch","test")->GetListOfFunctions()->Add(fSPD2LL->Clone());
    }
    //_________
        
    hSPDcorrectionVsEta = Histo("AliAnalysisMuMuNch","SPDcorrectionVsEta"); // Set the individual event histos
    hNTrackletVsEta = Histo("AliAnalysisMuMuNch","NTrackletVsEta");
    hNTrackletVsPhi = Histo("AliAnalysisMuMuNch","NTrackletVsPhi");
        
    hSPDcorrectionVsEta->Reset(); // Reset of the individual event histos
    hNTrackletVsEta->Reset();
    hNTrackletVsPhi->Reset();

    //____Data (or Reco if simu) tracklets computation___
    if ( tracklets && vertexSPDFound && isEtaRangeValid )
    {
      
      Histo("AliAnalysisMuMuNch","test")->Fill(SPDZv,fetaRange[1]);
      Histo("AliAnalysisMuMuNch","test")->Fill(SPDZv,fetaRange[0]);
      
      nTracklets = tracklets->GetNumberOfTracklets();
      
      Double_t SPDr; // SPD efficiency (If we use the 'mean correction' or no correction it will be always 1)
      for (Int_t i = 0 ; i < nTracklets ; i++)
      {
        thetaTracklet = tracklets->GetTheta(i);
        etaTracklet = -TMath::Log(TMath::Tan(thetaTracklet/2.));
        if ( etaTracklet < fetaRange[0] || etaTracklet > fetaRange[1] ) continue; // Avoid tracklets out of the eta SPD acceptance or out of the eta cut
        
        SPDr = GetSPDCorrection(SPDZv,etaTracklet); // Get SPD AccxEff for the corresponding [zvtx,eta] bin (If we use the 'mean correction' or no correction it will be always 1)
        
        Int_t bin = fEtaAxis->FindBin(etaTracklet);
        
        nTrackletsInRange++;
        if ( SPDr!=0. && SPDr <= 2.5) // Threshold to reduce border effects
        {
          hSPDcorrectionVsEta->SetBinContent(bin,SPDr);
          hNTrackletVsEta->Fill(etaTracklet);
          hNTrackletVsPhi->Fill(tracklets->GetPhi(i));
          
          nTrackletsInRangeSPDOk++;
        }
        else // If the correction is above the threshold we store a -1 in the correction to skip this eta bin int the fill method (can just happen using SPD AccxEff correction)
        {
          hSPDcorrectionVsEta->SetBinContent(bin,-1.0);
        }
      }
      
      nchList->Add(new TParameter<Double_t>("Ntr",nTrackletsInRange));
      nchList->Add(new TParameter<Double_t>("NtrSPDOk",nTrackletsInRangeSPDOk)); // if there is not SPD AccxEff correction this will be equal to Ntr
      
      //___ Fill the out-of-eta-range bins with -1.0
      
      binMin = fEtaAxis->FindBin(fetaRange[0]);
      binMax = fEtaAxis->FindBin(fetaRange[1]);
      
      for ( Int_t i = 1; i < binMin; ++i )
      {
        hSPDcorrectionVsEta->SetBinContent(i,-1.0);
      }
      for ( Int_t i = binMax + 1 ; i <= fEtaAxis->GetNbins(); ++i )
      {
        hSPDcorrectionVsEta->SetBinContent(i,-1.0);
      }

      if ( fSPDOneOverAccxEff || fSPDMeanTracklets ) // Here we correct the raw tracklets
      {
        if ( fSPDOneOverAccxEff ) // In this case the correction is the SPD accxEff
        {
          
          //----- Mean dNchdEta computation to set it into the event
          for (Int_t j = binMin/*1*/ ; j <= binMax/*fEtaAxis->GetNbins()*/ ; j++) // Loop over eta bins
          {
            Double_t correction = hSPDcorrectionVsEta->GetBinContent(j);
            
            Double_t eta = fEtaAxis->GetBinCenter(j);
            
            if ( correction < 0 ) continue; // If the correction is < 0 we skip that eta bin
            else if ( correction == 0.0 ) // If the correction is 0 we have no tracklets in that eta bin
            {
              Double_t spdCorr = GetSPDCorrection(SPDZv,eta);
              if ( spdCorr == 0. || spdCorr > 2.5) continue; // If the correction in the eta bin is not within the threshold we do not count the "0"(that eta bin will not count for nBins)
            }
            
            nch += hNTrackletVsEta->GetBinContent(j) * correction; // Number of charged particles (tracklets*(1/SPDAccxEff))
            
            ++nBins; // We sum up the number of bins entering in the computation
          }
          
          nchList->Add(new TParameter<Double_t>("Nch",nch));
          
          Double_t meandNchdEta(0.);// We compute the mean dNch/dEta in the event
          if ( nBins >  0 )
          {
            meandNchdEta = nch / (nBins*fEtaAxis->GetBinWidth(5));
          }
          
          nchList->Add(new TParameter<Double_t>("MeandNchdEta",meandNchdEta)); // We add the mean dNch/dEta to the event. It will serve us as a multiplicity estimator.
        }
//        else if( fSPDMeanTracklets ) // In this case the correction is the 'mean correction'
        if( fSPDMeanTracklets ) // In this case the correction is the 'mean correction'
        {
          SPDr = GetTrackletsMeanCorrection(SPDZv,nTrackletsInRange); // Get 'mean correction' for the zvtx
          
          if ( SPDr < -999.) nch = -1;
          else nch = nTrackletsInRange + SPDr;
          
          nchList->Add(new TParameter<Double_t>("NtrCorr",nch));// We add the corrected number of tracklets to the event. It will serve us as a multiplicity estimator.
        }
//        else AliFatal("Unknown tracklets correction class");
      }
      
    }
    else // If there is no tracklets or vertex object or the eta range is not valid everything is invalidated
    {
      for ( Int_t i = 1 ; i <= hSPDcorrectionVsEta->GetNbinsX() ; ++i )
      {
        hSPDcorrectionVsEta->SetBinContent(i,-1.0);
      }
      
      nchList->Add(new TParameter<Double_t>("Ntr",-1)); //We set an invalid Ntracklets
      
      if ( fSPDOneOverAccxEff ) // In this case the correction is the SPD accxEff
      {
        nchList->Add(new TParameter<Double_t>("MeandNchdEta",-1)); //We set an invalid dNch/deta
      }
      else if( fSPDMeanTracklets ) // In this case the correction is the 'mean correction'
      {
        nchList->Add(new TParameter<Double_t>("NtrCorr",-1)); //We set an invalid Ntracklets corrected
      }       
    }
    //_______

    //____ We add the tracklets histos to the event
    nchList->Add(hSPDcorrectionVsEta->Clone());
    nchList->Add(hNTrackletVsEta->Clone());
    nchList->Add(hNTrackletVsPhi->Clone());
    //_______
  }
  //_______
  
  //____Input MC multiplicity computation ___
  if ( HasMC() )
  {
    Double_t etaRange[2];
    if ( !fResolution ) //When computing resolutions we dont do anything else
    {
      Double_t MCZv = AliAnalysisMuonUtility::GetMCVertexZ(Event(),MCEvent());
      GetEtaRangeSPD(MCZv,etaRange);
      
      hNchVsEta = MCHisto("AliAnalysisMuMuNch","NchVsEta");
      hSPDcorrectionVsEta = MCHisto("AliAnalysisMuMuNch","SPDcorrectionVsEta");
      TH1* hNchVsPhi = MCHisto("AliAnalysisMuMuNch","NchVsPhi");
      hNTrackletVsEta = MCHisto("AliAnalysisMuMuNch","NTrackletVsEta");
      hNTrackletVsPhi = MCHisto("AliAnalysisMuMuNch","NTrackletVsPhi");
      
      hNTrackletVsEta->Reset();
      hNchVsEta->Reset();
      hNchVsPhi->Reset();
      hSPDcorrectionVsEta->Reset();
      hNTrackletVsPhi->Reset();
      
      //___ Fill the out-of-eta-range bins with -1.0
      binMin = fEtaAxis->FindBin(etaRange[0]);
      binMax = fEtaAxis->FindBin(etaRange[1]);
      
      for ( Int_t i = 1; i < binMin; ++i )
      {
        hSPDcorrectionVsEta->SetBinContent(i,-1.0);
      }
      for ( Int_t i = binMax + 1 ; i <= fEtaAxis->GetNbins(); ++i )
      {
        hSPDcorrectionVsEta->SetBinContent(i,-1.0);
      }
      for ( Int_t i = binMin; i <= binMax; ++i ) // Fill the bins inside the eta range with +1
      {
        hSPDcorrectionVsEta->SetBinContent(i,1.0); // FIXME: Do we really want all the tracklets to be counted for MC purposes? or we want to count only those in active regions of the SPD as for Data computations?
      }
      //___
      
      Int_t nMCTracks = MCEvent()->GetNumberOfTracks(); // number of MC tracks
      
      for ( Int_t i = 0; i < nMCTracks ; ++i ) //Loop over generated tracks
      {
//        AliAODMCParticle* AODpart = static_cast<AliAODMCParticle*>(mcEvent->GetTrack(i));
        AliMCParticle* MCpart = static_cast<AliMCParticle*>(mcEvent->GetTrack(i));
        
        Bool_t isPP(kFALSE);
        Int_t nGentorDex(-2); //Generator Index: -1 as default;
        //                  0 is for DPMJET;
        //                  1, 2 ... are for other generators
        
        if ( static_cast<const AliVEvent*>(event)->IsA() == AliAODEvent::Class() )
        {
          AliAODMCParticle* AODpart = static_cast<AliAODMCParticle*>(mcEvent->GetTrack(i));
          isPP = AODpart->IsPhysicalPrimary();
          nGentorDex = AODpart->GetGeneratorIndex();
        }
        else if ( static_cast<const AliVEvent*>(event)->IsA() == AliESDEvent::Class() )
        {
          isPP = MCEvent()->IsPhysicalPrimary(i);
          nGentorDex = MCpart->GetGeneratorIndex();
        }
        
        
//        if ( AODpart->IsPhysicalPrimary() ) // We take only particles produced in the collision (Particles produced in the collision including products of strong and electromagnetic decay and excluding feed-down from weak decays of strange particles)
        if ( isPP && nGentorDex==0 )
        {
//          if ( !IsMCtrackFromGenerator(i) ) continue; // Select only the particles generated by the desired generator
          
          if ( MCpart->Charge()!=0 ) // We take only charged particles
          {
            hNchVsEta->Fill(MCpart->Eta());
            hNchVsPhi->Fill(MCpart->Phi());
            
          }
        }
      }
      
      nchList->Add(hNchVsEta->Clone("MCNchVsEta"));
      nchList->Add(hNchVsPhi->Clone("MCNchVsPhi"));
      nchList->Add(hSPDcorrectionVsEta->Clone("MCSPDcorrectionVsEta"));
    }
              //__Bkg tracklets and Resolution estimation __
    if ( tracklets ) // We can compute the Bkg tracklets and resolution only if we have the tracklets object in the event
    {
      TH1* hNBkgTrackletsVSEta(0x0); // Pointer for the Bkg histo
      if ( !fResolution )
      {
        hNBkgTrackletsVSEta = Histo("AliAnalysisMuMuNch","NBkgTrackletsVSEta");
        
        hNBkgTrackletsVSEta->Reset();
      }
      
      nTracklets = tracklets->GetNumberOfTracklets();
      
      for (Int_t i = 0 ; i < nTracklets ; i++) // Loop over tracklets to check if they come or not from the same MC particle
      {
        Int_t label1 = tracklets->GetLabel(i,0);
        Int_t label2 = tracklets->GetLabel(i,1);

//        if ( !IsMCtrackFromGenerator(label1) || !IsMCtrackFromGenerator(label2) ) continue; // Select only the particles generated by the desired generator
        
        thetaTracklet = tracklets->GetTheta(i);
        etaTracklet = -TMath::Log(TMath::Tan(thetaTracklet/2.));
        phiTracklet = tracklets->GetPhi(i);
        
        if ( !fResolution )
        {
          hNTrackletVsEta->Fill(etaTracklet);
          hNTrackletVsPhi->Fill(phiTracklet);
        }
        
        if ( label1 != label2 && !fResolution ) hNBkgTrackletsVSEta->Fill(etaTracklet); // Tracklets not comming from the same MC particle are Bkg
        else if ( !fSPDOneOverAccxEff && !fSPDMeanTracklets && fResolution ) // Compute the resolutions with the tracklets comming from the same MC particle
        {
//          AliAODMCParticle* AODpartMC = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(label1));
          AliMCParticle* partMC = static_cast<AliMCParticle*>(MCEvent()->GetTrack(label1));
          Double_t etaTrackletMC = partMC->Eta();
          Double_t phiTrackletMC = partMC->Phi();
          
          // Resolution variables
          nchList->Add(new TParameter<Double_t>(Form("EtaReco%d",label1),etaTracklet));
          nchList->Add(new TParameter<Double_t>(Form("EtaMC%d",label1),etaTrackletMC));
          
          nchList->Add(new TParameter<Double_t>(Form("PhiReco%d",label1),phiTracklet));
          nchList->Add(new TParameter<Double_t>(Form("PhiMC%d",label1),phiTrackletMC));
                 
        }
      }
      
      if (!fResolution )
      {
        nchList->Add(hNTrackletVsEta->Clone("MCNTrackletVsEta"));
        nchList->Add(hNTrackletVsPhi->Clone("MCNTrackletVsPhi"));
        nchList->Add(hNBkgTrackletsVSEta->Clone());
      }
      
    }
    
  }
  //_______
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuNch::SetRun(const AliInputEventHandler* eventHandler)
{
  // For each new run this method sets the corresponding SPD correction from SPDCorrection list using the SPDCorrection map

  if ( !fSPDCorrectionList ) return;
  
//  AliAODEvent* event = static_cast<AliAODEvent*>(eventHandler->GetEvent());
   AliVEvent* event = static_cast<AliVEvent*>(eventHandler->GetEvent());

  TString run(Form("%d",event->GetRunNumber())); // Get the run number
  TObjString*  SPDCorrectionKey = static_cast<TObjString*>(fSPDCorrectionMap->GetValue(run)); // Get the corresponding SPDCorrection map key for the run
  
  Int_t SPDCorrectionIndex(0);
  if ( SPDCorrectionKey ) SPDCorrectionIndex = SPDCorrectionKey->String().Atoi(); // Converts the key into Int, if found
  else AliFatal(Form("No SPD correction found for run %d",run.Atoi()));
  
  TH1* hSPDCorrection = static_cast<TH1*>(fSPDCorrectionList->At(SPDCorrectionIndex)); // Gets the SPD correction at key position from the list
  
  if ( static_cast<TH1*>(hSPDCorrection)->IsA() == TH2::Class() ) fSPDOneOverAccxEff = static_cast<TH2F*>(hSPDCorrection); // Sets the corresponding correction data member
  else if ( static_cast<TH1*>(hSPDCorrection)->IsA() == TProfile::Class() ) fSPDMeanTracklets = static_cast<TProfile*>(hSPDCorrection);
  else AliFatal("Unrecognized correction class");
  
  AliInfo(Form("Using correction %s for run %s",hSPDCorrection->GetName(),run.Data()));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuNch::Terminate(Option_t *)
{
  /// Called once at the end of the query
  if ( !HistogramCollection() ) return;

  if ( HistogramCollection()->FindObject(Form("/%s/AliAnalysisMuMuNch/NTrackletVsEta",MCInputPrefix())) )
  {
    HistogramCollection()->Remove(Form("/%s/AliAnalysisMuMuNch/NTrackletVsEta",MCInputPrefix()));
    HistogramCollection()->Remove(Form("/%s/AliAnalysisMuMuNch/NTrackletVsPhi",MCInputPrefix()));
    HistogramCollection()->Remove(Form("/%s/AliAnalysisMuMuNch/NchVsEta",MCInputPrefix()));
    HistogramCollection()->Remove(Form("/%s/AliAnalysisMuMuNch/NchVsPhi",MCInputPrefix()));
    HistogramCollection()->Remove(Form("/%s/AliAnalysisMuMuNch/SPDcorrectionVsEta",MCInputPrefix()));
    HistogramCollection()->Remove("/AliAnalysisMuMuNch/NBkgTrackletsVSEta");
  }
  
  if ( HistogramCollection()->FindObject("/AliAnalysisMuMuNch/NTrackletVsEta") )
  {
    HistogramCollection()->Remove("/AliAnalysisMuMuNch/NTrackletVsEta");
    HistogramCollection()->Remove("/AliAnalysisMuMuNch/test");
    HistogramCollection()->Remove("/AliAnalysisMuMuNch/NTrackletVsPhi");
    HistogramCollection()->Remove("/AliAnalysisMuMuNch/SPDcorrectionVsEta");
  }
  //____ Compute dNchdEta histo
  TObjArray* idArr =  HistogramCollection()->SortAllIdentifiers();

  TIter next(idArr);
  TObjString* id;
  
  while ( (id = static_cast<TObjString*>(next())) )
  {
    TProfile* p = static_cast<TProfile*>(HistogramCollection()->FindObject(Form("%s%s",id->GetName(),"MeanNchVsEta")));

    if ( !p ) continue;
    
    TH1* h = new TH1F("MeandNchdEta","Event averaged dN_{ch}/d#eta ;#eta;<dN_{ch}/d#eta>",fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());

    if ( p->GetNbinsX() != h->GetNbinsX() || p->GetXaxis()->GetXmin() != h->GetXaxis()->GetXmin() || p->GetXaxis()->GetXmax() != h->GetXaxis()->GetXmax() )
    {
      AliError("ERROR: Cannot compute MeandNchdEta since the binning doesn't match with MeanNchVsEta histo");
      continue;
    }
    
    for ( Int_t i = 1 ; i < h->GetNbinsX() ; i++ )
    {
      h->SetBinContent(i,p->GetBinContent(i)/p->GetBinWidth(i));
      h->SetBinError(i,p->GetBinError(i)/p->GetBinWidth(i));
      h->SetEntries(p->GetEntries());
    }
  
    HistogramCollection()->Adopt(Form("%s",id->GetName()),h);
  }
    
  delete idArr;
  //___
}
