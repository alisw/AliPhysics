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

#include "AliAODTracklets.h"
#include "AliAnalysisMuonUtility.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliESDUtils.h"
#include "TMath.h"
#include "AliAnalysisMuMuCutRegistry.h"
#include "AliAnalysisMuMuCutElement.h"
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

//_____________________________________________________________________________
AliAnalysisMuMuNch::AliAnalysisMuMuNch(TH2* spdCorrection, Double_t etaMin, Double_t etaMax
                                      , Double_t zMin, Double_t zMax,Bool_t disableHistos, Bool_t computeResolution)
: AliAnalysisMuMuBase(),
fSPDCorrection(0x0),
fEtaAxis(new TAxis(TMath::Nint(10./0.1),-5.,5.)),
fZAxis(new TAxis(TMath::Nint(80/0.25),-40.,40.)),
fCurrentEvent(0x0),
fEtaMin(etaMin),
fEtaMax(etaMax),
fZMin(zMin),
fZMax(zMax),
fResolution(computeResolution)
{
  if ( spdCorrection )
  {
    fSPDCorrection = static_cast<TH2F*>(spdCorrection->Clone());
    fSPDCorrection->SetDirectory(0);
  }
  DefineSPDAcceptance();
  
  if ( disableHistos )
  {
    DisableHistograms("*");
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuNch::~AliAnalysisMuMuNch()
{
  delete fSPDCorrection;
  delete fEtaAxis;
  delete fZAxis;
  delete fSPD1LR;
  delete fSPD1LL;
  delete fSPD2LR;
  delete fSPD2LL;
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

  Double_t multMin = 0.;  //Tracklets multiplicity range
  Double_t multMax = 500.;
  Int_t nbinsMult = GetNbins(multMin,multMax,1.);
  
  Double_t phimin = 0.; //Phi range
  Double_t phimax = 2*TMath::Pi();
  Int_t nphibins = GetNbins(phimin,phimax,0.05);
  
  if ( !fSPDCorrection && fResolution ) // Resolution histos
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

  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsPhi","Number of tracklets vs Z vertex vs #phi;Z vertex;#phi",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),nphibins,phimin,phimax);
  AttachSPDAcceptance(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsPhi");// Attach the SPD acc curves
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsEta","Number of tracklets vs ZVertex vs #eta;ZVertex;#eta",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
  AttachSPDAcceptance(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsEta");
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Tracklets","Number of tracklets distribution;N_{Tracklets};N_{events}",nbinsMult,multMin,multMax);
    
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"NBkgTrackletsVsZVertexVsEta","Number of background tracklets vs ZVertex vs #eta;ZVertex;#eta",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
  AttachSPDAcceptance(kHistoForMCInput,eventSelection,triggerClassName,centrality,"NBkgTrackletsVsZVertexVsEta");
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"NchVsZVertexVsEta","Number of charged particles vs ZVertex vs #eta;ZVertex;#eta",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax());
  AttachSPDAcceptance(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"NchVsZVertexVsEta");
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"NchVsZVertexVsPhi","Number of charged particles vs ZVertex vs #phi;ZVertex;#eta",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),nphibins,phimin,phimax);
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"EventsVsZVertexVsEta","Effective number of events vs ZVertex vs #eta;ZVertex;#eta",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax()); // Fill 1 unit in each "touched" eta bin per event (represents the eta bins in which each event contributes)
  AttachSPDAcceptance(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"EventsVsZVertexVsEta");
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Nch","Number of charged particles distribution;N_{ch};N_{events}",nbinsMult,multMin,multMax);
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"dNchdEta","<dNchdEta> distribution;dN_{ch}/d#eta;N_{events}",nbinsMult,multMin,multMax);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"TrackletsVsNch","Number of tracklets vs Number of charged particles;N_{ch};N_{tracklets}",nbinsMult,multMin,multMax,nbinsMult,multMin,multMax); //Response matrix
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"V0AMultVsNch","V0A multiplicity vs Number of charged particles;N_{ch};V0A Mult",nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"V0CMultVsNch","V0C multiplicity vs Number of charged particles;N_{ch};V0C Mult",nbinsMult,multMin,multMax,nbinsMult,multMin,multMax);
  
  // profile histograms
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"MeanTrackletsVsEta","Mean number of tracklets vs #eta;#eta;<N_{Tracklets}>",fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax(),0);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"MeanTrackletsVsZVertex","Mean number of tracklets vs Z vertex;Z vertex;<N_{Tracklets}>",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),0);

  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"MeanNchVsEta","Mean number of charged particles vs #eta;#eta;<N_{ch}>",fEtaAxis->GetNbins(),fEtaAxis->GetXmin(),fEtaAxis->GetXmax(),0); // Each bin has to be divided by the binwidth to became dNch/dEta (Done in the terminate and stored in MeandNchdEta histo)
  
   CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"MeanNchVsZVertex","Mean number of charged particles vs Z vertex;Z vertex;<N_{ch}>",fZAxis->GetNbins(),fZAxis->GetXmin(),fZAxis->GetXmax(),0);
    
}

//_____________________________________________________________________________
void AliAnalysisMuMuNch::DefineSPDAcceptance()
{
  // Defines the functions ( eta = f(z) )of the edges (right/left) of the inner and outer SPD layers
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
  if ( dataType & kHistoForData )
  {
    if( !Histo(eventSelection,triggerClassName,centrality,histoname) )
    {
      AliError(Form("ERROR: SPD Acceptance curves attach failed. Histo /%s/%s/%s/%s not found",eventSelection,triggerClassName,centrality,histoname));
      return;
    }
    
    Histo(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD1LR);
    Histo(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD1LL);
    Histo(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD2LR);
    Histo(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD2LL);
  }
  if ( dataType & kHistoForMCInput )
  {
    if( !MCHisto(eventSelection,triggerClassName,centrality,histoname) )
    {
      AliError(Form("ERROR: SPD Acceptance curves attach failed. MC Histo /%s/%s/%s/%s not found",eventSelection,triggerClassName,centrality,histoname));
      return;
    }
    
    MCHisto(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD1LR);
    MCHisto(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD1LL);
    MCHisto(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD2LR);
    MCHisto(eventSelection,triggerClassName,centrality,histoname)->GetListOfFunctions()->Add(fSPD2LL);
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

  if ( !AliAnalysisMuonUtility::IsAODEvent(Event()) )
  {
    AliError("Don't know how to deal with ESDs...");
    return;
  }
  
  if ( HasMC() && !fSPDCorrection ) // We have MC but no correction (SPD correction computation mode)so we skip the method
  {
    return;
  }

  AliAODEvent* aod = static_cast<AliAODEvent*>(Event());

  AliVVertex* vertex = aod->GetPrimaryVertexSPD();
  
  TList* nchList = static_cast<TList*>(Event()->FindListObject("NCH"));
  
  if (!nchList || nchList->IsEmpty() ) // Empty NCH means that there is no SPD vertex ( see SetEvent() ) when runing on data.
  {
    return;
  }
  
  Double_t SPDZv;
  
  if ( !vertex || vertex->GetZ() == 0.0 ) // Running in Simu the spdZ == 0 means no SPD info. In data, avoid breaks in events w/o SPD vertex
  {
    SPDZv = -40.;
  } 
  
  else SPDZv = vertex->GetZ();
  
  TH1* hSPDcorrectionVsEta = static_cast<TH1*>(nchList->FindObject("SPDcorrectionVsEta"));
  TH1* hNTrackletVsEta = static_cast<TH1*>(nchList->FindObject("NTrackletVsEta"));
  TH1* hNTrackletVsPhi = static_cast<TH1*>(nchList->FindObject("NTrackletVsPhi"));
  
  TH1* hNchVsEta =static_cast<TH1*>(hNTrackletVsEta->Clone("NchVsEta"));
  
  TProfile* hMeanTrackletsVsEta = static_cast<TProfile*>(Histo(eventSelection,triggerClassName,centrality,"MeanTrackletsVsEta"));
  TProfile* hMeanNchVsEta = static_cast<TProfile*>(Histo(eventSelection,triggerClassName,centrality,"MeanNchVsEta"));

  TH2* hEventsVsZVertexVsEta = static_cast<TH2*>(Histo(eventSelection,triggerClassName,centrality,"EventsVsZVertexVsEta"));
    
  Int_t nBins(0);

  Double_t nch(0.0);
  Double_t nTracklets(0.0);
  
  for (Int_t j = 1 ; j <= fEtaAxis->GetNbins() ; j++) // Loop over eta bins
  {
    Double_t correction = hSPDcorrectionVsEta->GetBinContent(j);
    
    Double_t eta = fEtaAxis->GetBinCenter(j);
    
    if ( correction < 0 ) continue;
    else if ( correction == 0.0 ) // No tracklets found in this eta bin.
    {
      correction = GetSPDCorrection(SPDZv,eta);

      if ( correction == 0. || correction > 2.5) continue; // We need to know if the eta bin is in a region we want to take into account to count or not the zero
    }
  
    Double_t ntr = hNTrackletVsEta->GetBinContent(j); // Tracklets in eta bin
    
    nch += ntr * correction; // Number of charged particles (corrected tracklets)
    
    nTracklets += ntr;
    
    hMeanTrackletsVsEta->Fill(eta,ntr); // Fill the number of tracklets of each eta bin in the profile

    hMeanNchVsEta->Fill(eta,ntr*correction);

    ++nBins; // We sum up the number of bins entering in the computation
    
    // Fill the number of charged particles of each eta bin in the profile
    hNchVsEta->SetBinContent( j, hNchVsEta->GetBinContent(j) * correction );
    
    hEventsVsZVertexVsEta->Fill(SPDZv,eta,1.0); // Fill 1 count each eta bin where the events contributes
  }
  
  AddHisto(eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsPhi",SPDZv,hNTrackletVsPhi);
  AddHisto(eventSelection,triggerClassName,centrality,"TrackletsVsZVertexVsEta",SPDZv,hNTrackletVsEta);

  AddHisto(eventSelection,triggerClassName,centrality,"NchVsZVertexVsEta",SPDZv,hNchVsEta);
  
  Histo(eventSelection,triggerClassName,centrality,"Tracklets")->Fill(nTracklets);
  Histo(eventSelection,triggerClassName,centrality,"MeanTrackletsVsZVertex")->Fill(SPDZv,nTracklets);
  Histo(eventSelection,triggerClassName,centrality,"MeanNchVsZVertex")->Fill(SPDZv,nch);
  
  Histo(eventSelection,triggerClassName,centrality,"TrackletsVsNch")->Fill(nTracklets,nch);
  Histo(eventSelection,triggerClassName,centrality,"Nch")->Fill(nch);
  
  Double_t V0AMult = 0.;
  Double_t V0CMult = 0.;
  
  AliVVZERO* vzero = aod->GetVZEROData();
  if (vzero)
  {
    Double_t multV0A = vzero->GetMTotV0A();
    V0AMult = AliESDUtils::GetCorrV0A(multV0A,SPDZv);
    Double_t multV0C = vzero->GetMTotV0C();
    V0CMult = AliESDUtils::GetCorrV0C(multV0C,SPDZv);
    
    Histo(eventSelection,triggerClassName,centrality,"V0AMultVsNch")->Fill(V0AMult,nch);
    Histo(eventSelection,triggerClassName,centrality,"V0CMultVsNch")->Fill(V0CMult,nch);

  }

  
  delete hNchVsEta;
  
  // Mean dNch/dEta computation
  Double_t meandNchdEta(0.);
  
  if ( nBins >  0 )
  {
    meandNchdEta = nch / (nBins*fEtaAxis->GetBinWidth(5)); // Divide by nBins to get the mean and by the binWidht to get the d/dEta
  }
  
  Histo(eventSelection,triggerClassName,centrality,"dNchdEta")->Fill(meandNchdEta);
  
  
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuNch::FillHistosForMCEvent(const char* eventSelection,const char* triggerClassName,const char* centrality)
{
  /// Fill input MC multiplicity histos
  
  if ( IsHistogrammingDisabled() ) return;
  
  TList* nchList = static_cast<TList*>(Event()->FindListObject("NCH"));
  
  if (!nchList || nchList->IsEmpty())
  {
    return;
  }

  Double_t MCZv = AliAnalysisMuonUtility::GetMCVertexZ(Event(),MCEvent()); // Definition of MC generated z vertex
  AliVVertex* vertex = static_cast<AliAODEvent*>(Event())->GetPrimaryVertexSPD();
  
  Double_t SPDZv(0.);
  
  //____Resolution Histos___
  if ( !fSPDCorrection && fResolution )
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
      while ( nchList->At(i)->IsA() != TObjString::Class() ) // In case there is a diferent object, just to skip it
      {
        i++;
      }
      
      TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(nchList->At(i));

      
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
  for (Int_t j = 1 ; j <= fEtaAxis->GetNbins() ; j++) //Loop over eta bins
  {
    Double_t correction = hSPDcorrectionVsEta->GetBinContent(j);
    
    if ( correction < 0 ) continue; // We count just the particles in the SPD acceptance.
    
    nTracklets += hNTrackletVsEta->GetBinContent(j); // Reco tracklets
    
    ++nBins; // We sum up the number of bins entering in the computation
    
    Double_t eta = fEtaAxis->GetBinCenter(j);
    Double_t nch = hNchVsEta->GetBinContent(j); // Generated particles
    
    nchSum += nch;
    
    hMeanNchVsEta->Fill(eta,nch); // Fill the number of charged particles of each eta bin in the profile
    
    hEventsVsZVertexVsEta->Fill(MCZv,eta,1.0); // Fill 1 count each eta bin where the events contributes
  }
  
  MCHisto(eventSelection,triggerClassName,centrality,"Tracklets")->Fill(nTracklets);
  
  if ( vertex )
  {
    SPDZv = vertex->GetZ();
  }

  Bool_t isMChisto = kTRUE; // Used to get the MC histos in the Add method
  
  AddHisto(eventSelection,triggerClassName,centrality,"NBkgTrackletsVsZVertexVsEta",SPDZv,hNBkgTrackletsVSEta,isMChisto);
  AddHisto(eventSelection,triggerClassName,centrality,"NchVsZVertexVsPhi",MCZv,hNchVsPhi,isMChisto);
  AddHisto(eventSelection,triggerClassName,centrality,"NchVsZVertexVsEta",MCZv,hNchVsEta,isMChisto);
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
void AliAnalysisMuMuNch::GetEtaRangeSPD(Double_t spdZVertex, Double_t etaRange[])
{
  // Returns the SPD eta range for a given z vertex.
  
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
  
  etaRange[0] = etaMin;
  etaRange[1] = etaMax;
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuNch::GetSPDCorrection(Double_t zvert, Double_t eta) const
{
  if (!fSPDCorrection)
  {
    AliFatal("ERROR: No SPD Correction");
    return 0;
  }
  Int_t bin = fSPDCorrection->FindBin(zvert,eta);
  return fSPDCorrection->GetBinContent(bin);
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

//_____________________________________________________________________________
void AliAnalysisMuMuNch::SetEvent(AliVEvent* event, AliMCEvent* mcEvent)
{
  /// Set the event, compute multiplicities and add them to the event

  if ( event->IsA() != AliAODEvent::Class() )
  {
    AliError("Only working for AODs for the moment.");
    return;
  }
  
  AliAnalysisMuMuBase::SetEvent(event,mcEvent); // To have Event() and MCEvent() method working
  
  TList* nchList = static_cast<TList*>(event->FindListObject("NCH")); // Define the list with the NCH info for each event
  if (!nchList)
  {
    nchList = new TList;
    nchList->SetOwner(kTRUE);
    nchList->SetName("NCH");
    event->AddObject(nchList);
  }
  
  const AliAODVertex* vertexSPD = static_cast<const AliAODEvent*>(Event())->GetPrimaryVertexSPD(); // SPD vertex object
  AliAODTracklets* tracklets = static_cast<const AliAODEvent*>(Event())->GetTracklets(); // Tracklets object
  
  TH1* hSPDcorrectionVsEta(0x0); // Pointers for the individual event histos
  TH1* hNchVsEta(0x0);
  TH1* hNTrackletVsEta(0x0);
  TH1* hNTrackletVsPhi(0x0);
  
  //_______Create once the histos with the individual event "properties" (cleared at the beginning of each new event)_______
  if ( !fResolution ) // When computing resolutions we dont do anything else
  {
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
      
      Histo("AliAnalysisMuMuNch","test")->GetListOfFunctions()->Add(fSPD1LR);
      Histo("AliAnalysisMuMuNch","test")->GetListOfFunctions()->Add(fSPD1LL);
      Histo("AliAnalysisMuMuNch","test")->GetListOfFunctions()->Add(fSPD2LR);
      Histo("AliAnalysisMuMuNch","test")->GetListOfFunctions()->Add(fSPD2LL);
    }
    //_________
    
    
    hSPDcorrectionVsEta = Histo("AliAnalysisMuMuNch","SPDcorrectionVsEta"); // Set the individual event histos
    hNTrackletVsEta = Histo("AliAnalysisMuMuNch","NTrackletVsEta");
    hNTrackletVsPhi = Histo("AliAnalysisMuMuNch","NTrackletVsPhi");
    
    
    hSPDcorrectionVsEta->Reset(); // Reset of the individual event histos
    hNTrackletVsEta->Reset();
    hNTrackletVsPhi->Reset();
  }
  
  Double_t etaRange[2]; // Variables we will use in the multiplicity computation
  Int_t binMin,binMax;
  Int_t nTracklets(0);
  Double_t thetaTracklet(0.),phiTracklet(0.),etaTracklet(0.);
  Double_t nch(0.0);
  Int_t nBins(0);
  
  
  //____Data (or Reco if simu) multiplicity computation___
  if ( fSPDCorrection && !fResolution ) // When computing resolutions we dont do anything else
  {
    if ( tracklets && vertexSPD )
    {
      Double_t SPDZv = vertexSPD->GetZ();

      GetEtaRangeSPD(SPDZv,etaRange);
      
      Histo("AliAnalysisMuMuNch","test")->Fill(SPDZv,etaRange[1]);
      Histo("AliAnalysisMuMuNch","test")->Fill(SPDZv,etaRange[0]);
      
      nTracklets = tracklets->GetNumberOfTracklets();
      
      Double_t SPDr;
      for (Int_t i = 0 ; i < nTracklets ; i++)
      {
        thetaTracklet = tracklets->GetTheta(i);
        etaTracklet = -TMath::Log(TMath::Tan(thetaTracklet/2.));
        if ( etaTracklet < etaRange[0] || etaTracklet > etaRange[1] ) continue; // Avoid tracklets out of the eta SPD acceptance or out of the eta cut
        
        SPDr = GetSPDCorrection(SPDZv,etaTracklet);
        
        Int_t bin = fEtaAxis->FindBin(etaTracklet);
        
        if ( SPDr!=0. && SPDr <= 2.5) // Threshold to reduce border effects
        {
          hSPDcorrectionVsEta->SetBinContent(bin,SPDr);
          hNTrackletVsEta->Fill(etaTracklet);
          hNTrackletVsPhi->Fill(tracklets->GetPhi(i));
        }
        else // If the correction is above the threshold we store a -1 in the correction to skip this eta bin int the fill method
        {
          hSPDcorrectionVsEta->SetBinContent(bin,-1.0);
        }
      }
      
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
      //____
      nchList->Clear(); // We clear the NCH list for this new event
      
      if ( !IsHistogrammingDisabled() )
      {
        nchList->Add(hSPDcorrectionVsEta->Clone());
        nchList->Add(hNTrackletVsEta->Clone());
        nchList->Add(hNTrackletVsPhi->Clone());
      }
      
          //----- Mean dNchdEta computation to set it into the event
      for (Int_t j = 1 ; j <= fEtaAxis->GetNbins() ; j++) // Loop over eta bins
      {
        Double_t correction = hSPDcorrectionVsEta->GetBinContent(j);
        
        Double_t eta = fEtaAxis->GetBinCenter(j);
        
        if ( correction < 0 ) continue; // If the correction is < 0 we skip that eta bin 
        else if ( correction == 0.0 ) // If the correction is 0 we have no tracklets in that eta bin
        {
          Double_t spdCorrection = GetSPDCorrection(SPDZv,eta);
          if ( spdCorrection == 0. || spdCorrection > 2.5) continue; // If the correction in the eta bin is not within the threshold we do not count the "0"(that eta bin will not count for nBins)
        }
        
        nch += hNTrackletVsEta->GetBinContent(j) * correction; // Number of charged particles (tracklets*SPDcorrection)
        
        ++nBins; // We sum up the number of bins entering in the computation
      }
      
      Double_t meandNchdEta(0.); // We compute the mean dNch/dEta in the event
      if ( nBins >  0 ) 
      {
        meandNchdEta = nch / (nBins*fEtaAxis->GetBinWidth(5));
      }
      
      nchList->Add(new TParameter<Double_t>("MeandNchdEta",meandNchdEta)); // We add the mean dNch/dEta to the event. It will serve us as a multiplicity estimator.
      //------
      
    }
    
    else  nchList->Clear(); //To clear the NCH list for MC in case the event has no reconstructed SPD vertex
  }
  //_______
  
  else nchList->Clear(); // To clear the NCH list for MC in case we have MC and no correction (SPD correction computation mode)

  
  //____Input MC multiplicity computation ___
  if ( HasMC() )
  {
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
        hSPDcorrectionVsEta->SetBinContent(i,1.0);
      }
      //___
      
      Int_t nMCTracks = MCEvent()->GetNumberOfTracks(); // MC number of MC tracks
      
      for ( Int_t i = 0; i < nMCTracks ; ++i ) //Loop over generated tracks
      {
        AliAODMCParticle* AODpart = static_cast<AliAODMCParticle*>(mcEvent->GetTrack(i));
        
        if ( AODpart->IsPhysicalPrimary() ) // We take only particles produced in the collision (Particles produced in the collision including products of strong and electromagnetic decay and excluding feed-down from weak decays of strange particles)
        {
          if ( AODpart->Charge()!=0 ) // We take only charged particles
          {
            hNchVsEta->Fill(AODpart->Eta());
            hNchVsPhi->Fill(AODpart->Phi());
            
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
        thetaTracklet = tracklets->GetTheta(i);
        etaTracklet = -TMath::Log(TMath::Tan(thetaTracklet/2.));
        phiTracklet = tracklets->GetPhi(i);
        
        
        if ( !fResolution )
        {
          hNTrackletVsEta->Fill(etaTracklet);
          hNTrackletVsPhi->Fill(phiTracklet);
        }
        
        Int_t label1 = tracklets->GetLabel(i,0);
        Int_t label2 = tracklets->GetLabel(i,1);  
        
        if (label1 != label2 ) // Tracklets not comming from the same MC particle are Bkg
        {
         if (!fResolution ) hNBkgTrackletsVSEta->Fill(etaTracklet);
        }
        else if ( !fSPDCorrection && fResolution ) // Compute the resolutions with the tracklets comming from the same MC particle
        {
          AliAODMCParticle* AODpartMC = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(label1));
          Double_t etaTrackletMC = AODpartMC->Eta();
          Double_t phiTrackletMC = AODpartMC->Phi();
          
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
void AliAnalysisMuMuNch::Terminate(Option_t *)
{
  /// Called once at the end of the query
  if ( !HistogramCollection() ) return;

  if ( HistogramCollection()->FindObject("/MCINPUT/AliAnalysisMuMuNch/NTrackletVsEta") )
  {
    HistogramCollection()->Remove("/MCINPUT/AliAnalysisMuMuNch/NTrackletVsEta");
    HistogramCollection()->Remove("/MCINPUT/AliAnalysisMuMuNch/NTrackletVsPhi");
    HistogramCollection()->Remove("/MCINPUT/AliAnalysisMuMuNch/NchVsEta");
    HistogramCollection()->Remove("/MCINPUT/AliAnalysisMuMuNch/NchVsPhi");
    HistogramCollection()->Remove("/MCINPUT/AliAnalysisMuMuNch/SPDcorrectionVsEta");
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
