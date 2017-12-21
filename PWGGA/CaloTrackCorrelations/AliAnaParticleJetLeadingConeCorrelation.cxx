/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// --- ROOT system ---
#include "TH2F.h"
#include "TClonesArray.h"
#include "TClass.h"
//#include "Riostream.h"

// --- Analysis system ---
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliCaloTrackReader.h"
#include "AliNeutralMesonSelection.h"
#include "AliAnaParticleJetLeadingConeCorrelation.h"
#include "AliCaloPID.h"
#include "AliCaloTrackParticleCorrelation.h"
#include "AliFiducialCut.h"

/// \cond CLASSIMP
ClassImp(AliAnaParticleJetLeadingConeCorrelation) ;
/// \endcond

//_________________________________________________________________________________
/// Default constructor. Initialize parameters.
//_________________________________________________________________________________
AliAnaParticleJetLeadingConeCorrelation::AliAnaParticleJetLeadingConeCorrelation() :
AliAnaCaloTrackCorrBaseClass(), fJetsOnlyInCTS(kFALSE), fPbPb(kFALSE),
fSeveralConeAndPtCuts(0),  fReMakeJet(0),
fDeltaPhiMaxCut(0.), fDeltaPhiMinCut(0.),
fLeadingRatioMaxCut(0.),  fLeadingRatioMinCut(0.),
fJetCTSRatioMaxCut(0.), fJetCTSRatioMinCut(0.),
fJetRatioMaxCut(0.),  fJetRatioMinCut(0.),
fJetNCone(0),fJetNPt(0), fJetCone(0),
fJetPtThreshold(0),fJetPtThresPbPb(0),
fPtTriggerSelectionCut(0.0), fSelect(0),fSelectIsolated(0),
fTrackVector(),fBkgMom(),fJetMom(),fJetConstMom(),
fLeadingMom(),fLeadingPi0Mom(),fLeadingPhoMom1(),fLeadingPhoMom2(),fLeadingChargeMom(),
// Histograms
fOutCont(0x0),
// Leading
fhChargedLeadingPt(0),fhChargedLeadingPhi(0),fhChargedLeadingEta(0),
fhChargedLeadingDeltaPt(0),fhChargedLeadingDeltaPhi(0),fhChargedLeadingDeltaEta(0),
fhChargedLeadingRatioPt(0),
fhNeutralLeadingPt(0),fhNeutralLeadingPhi(0),fhNeutralLeadingEta(0),
fhNeutralLeadingDeltaPt(0),fhNeutralLeadingDeltaPhi(0),fhNeutralLeadingDeltaEta(0),
fhNeutralLeadingRatioPt(0),fhChargedLeadingXi(0), fhNeutralLeadingXi(0),
fhChargedLeadingDeltaPhiRatioPt30(0), fhNeutralLeadingDeltaPhiRatioPt30(0),
fhChargedLeadingDeltaPhiRatioPt50(0), fhNeutralLeadingDeltaPhiRatioPt50(0),
//Jet
fhJetPt(0),fhJetRatioPt(0),fhJetDeltaPhi(0), fhJetDeltaEta(0),
fhJetLeadingRatioPt(0),fhJetLeadingDeltaPhi(0),fhJetLeadingDeltaEta(0),
fhJetFFz(0),fhJetFFxi(0),fhJetFFpt(0),fhJetNTracksInCone(0),
fhBkgPt(0),fhBkgRatioPt(0),fhBkgDeltaPhi(0), fhBkgDeltaEta(0),
fhBkgLeadingRatioPt(0),fhBkgLeadingDeltaPhi(0),fhBkgLeadingDeltaEta(0),
fhBkgFFz(0),fhBkgFFxi(0),fhBkgFFpt(0),fhBkgNTracksInCone(0),
// Several cones and thres histograms
fhJetPts(),fhJetRatioPts(),fhJetDeltaPhis(), fhJetDeltaEtas(),
fhJetLeadingRatioPts(),fhJetLeadingDeltaPhis(),fhJetLeadingDeltaEtas(),
fhJetFFzs(),fhJetFFxis(),fhJetFFpts(),fhJetNTracksInCones(),
fhBkgPts(),fhBkgRatioPts(),fhBkgDeltaPhis(), fhBkgDeltaEtas(),
fhBkgLeadingRatioPts(),fhBkgLeadingDeltaPhis(),fhBkgLeadingDeltaEtas(),
fhBkgFFzs(),fhBkgFFxis(),fhBkgFFpts(),fhBkgNTracksInCones()
{
  for(Int_t i = 0; i<6; i++)
  {
    fJetXMin1[i]     = 0.0 ;
    fJetXMin2[i]     = 0.0 ;
    fJetXMax1[i]     = 0.0 ;
    fJetXMax2[i]     = 0.0 ;
    fBkgMean[i]      = 0.0 ;
    fBkgRMS[i]       = 0.0 ;
    if( i < 2 ){
      fJetE1[i]        = 0.0 ;
      fJetE2[i]        = 0.0 ;
      fJetSigma1[i]    = 0.0 ;
      fJetSigma2[i]    = 0.0 ;
    }
  }
  
  // Several cones and thres histograms
  for(Int_t i = 0; i<5; i++)
  {
    fJetCones[i]         = 0.0 ;
    fJetNameCones[i]     = ""  ;
    fJetPtThres[i]      = 0.0 ;
    fJetNamePtThres[i]  = ""  ;
    for(Int_t j = 0; j<5; j++)
    {
      fhJetPts[i][j]=0 ;
      fhJetRatioPts[i][j]=0 ;
      fhJetDeltaPhis[i][j]=0 ;
      fhJetDeltaEtas[i][j]=0 ;
      fhJetLeadingRatioPts[i][j]=0 ;
      fhJetLeadingDeltaPhis[i][j]=0 ;
      fhJetLeadingDeltaEtas[i][j]=0 ;
      fhJetFFzs[i][j]=0 ;
      fhJetFFxis[i][j]=0 ;
      fhJetFFpts[i][j]=0 ;
      fhJetNTracksInCones[i][j]=0 ;
      fhBkgPts[i][j]=0 ;
      fhBkgRatioPts[i][j]=0 ;
      fhBkgDeltaPhis[i][j]=0 ;
      fhBkgDeltaEtas[i][j]=0 ;
      fhBkgLeadingRatioPts[i][j]=0 ;
      fhBkgLeadingDeltaPhis[i][j]=0 ;
      fhBkgLeadingDeltaEtas[i][j]=0 ;
      fhBkgFFzs[i][j]=0 ;
      fhBkgFFxis[i][j]=0 ;
      fhBkgFFpts[i][j]=0 ;
      fhBkgNTracksInCones[i][j]=0 ;
    }
  }
  
  InitParameters();
}

//____________________________________________________________________________
/// Calculate the ratio of the jet and trigger particle limit for the selection
/// WARNING: need to check what it does.
//____________________________________________________________________________
Double_t AliAnaParticleJetLeadingConeCorrelation::CalculateJetRatioLimit(Double_t ptg,
                                                                         const Double_t *par, const Double_t *x) const
{
  //printf("CalculateLimit: x1 %2.3f, x2%2.3f\n",x[0],x[1]);
  Double_t ePP = par[0] + par[1] * ptg ;
  Double_t sPP = par[2] + par[3] * ptg ;
  Double_t f   = x[0]   + x[1]   * ptg ;
  Double_t ePbPb = ePP + par[4] ;
  Double_t sPbPb = TMath::Sqrt(sPP*sPP+ par[5]*par[5]) ;
  Double_t rat = (ePbPb - sPbPb * f) / ptg ;
  //printf("CalculateLimit: ePP %2.3f, sPP %2.3f, f %2.3f\n", ePP, sPP, f);
  //printf("CalculateLimit: ePbPb %2.3f, sPbPb %2.3f, rat %2.3f\n", ePbPb, sPbPb, rat);
  return rat ;
}

//___________________________________________________________________________________________________
/// Fill jet and background histograms.
//___________________________________________________________________________________________________
void AliAnaParticleJetLeadingConeCorrelation::FillJetHistos(AliCaloTrackParticleCorrelation * particle,
                                                            const TLorentzVector jet,
                                                            const TString & type, const TString & lastname)
{
  Double_t ptTrig  = particle->Pt();
  Double_t ptJet   = jet.Pt();
  Double_t ptLead  = fLeadingMom.Pt();
  Double_t phiTrig = particle->Phi();
  Double_t phiJet  = jet.Phi();
  if(phiJet < 0) phiJet+=TMath::TwoPi();
  Double_t phiLead = fLeadingMom.Phi();
  if(phiLead < 0) phiLead+=TMath::TwoPi();
  Double_t etaTrig = particle->Eta();
  Double_t etaJet  = jet.Eta();
  Double_t etaLead = fLeadingMom.Eta();
  
  TH2F *h1 = 0x0;
  h1 = dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sPt%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())));
  if(h1) h1->Fill(ptTrig, ptJet, GetEventWeight());
  
  TH2F *h2 = 0x0;
  h2 = dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sRatioPt%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())));
  if(h2) h2->Fill(ptTrig, ptJet/ptTrig, GetEventWeight());
  
  TH2F *h3 = 0x0;
  h3 = dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sLeadingRatioPt%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())));
  if(h3) h3->Fill(ptTrig, ptLead/ptJet, GetEventWeight());
  
  //   dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sPhi%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())))->
  //     Fill(ptTrig,phiJet);
  TH2F *h4 = 0x0;
  h4 = dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sDeltaPhi%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())));
  if(h4) h4->Fill(ptTrig, phiJet-phiTrig, GetEventWeight());
    
  TH2F *h5 = 0x0;
  h5 = dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sLeadingDeltaPhi%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())));
  if(h5) h5->Fill(ptTrig, phiJet-phiLead, GetEventWeight());
  
  //   dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sEta%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())))->
  //     Fill(ptTrig,etaJet);
  TH2F *h6 = 0x0;
  h6 = dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sDeltaEta%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())));
  if(h6) h6->Fill(ptTrig, etaJet-etaTrig, GetEventWeight());
    
  TH2F *h7 = 0x0;
  h7 = dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sLeadingDeltaEta%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())));
  if(h7) h7->Fill(ptTrig, etaJet-etaLead, GetEventWeight());
  
  // Construct fragmentation function
  TObjArray * pl = new TObjArray;
  
  if     (type == "Jet") pl = particle->GetObjArray(Form("%sTracks",GetAODObjArrayName().Data()));
  else if(type == "Bkg") particle->GetObjArray(Form("%sTracksBkg",GetAODObjArrayName().Data()));
  
  if(!pl) return ;
  
  // Different pt cut for jet particles in different collisions systems
  // Only needed when jet is recalculated from AODs
  // Float_t ptcut = fJetPtThreshold;
  // if(fPbPb && !fSeveralConeAndPtCuts && ptTrig > fPtTriggerSelectionCut)  ptcut = fJetPtThresPbPb ;
  
  Int_t nTracksInCone = 0;
  
  for(Int_t ipr = 0;ipr < pl->GetEntriesFast() ; ipr ++ )
  {
    AliVTrack* track = dynamic_cast<AliVTrack *>(pl->At(ipr)) ;
    if(track)fTrackVector.SetXYZ(track->Px(),track->Py(),track->Pz());
    else AliDebug(1,"Track not available");
    
    //Recheck if particle is in jet cone
    if(fReMakeJet || fSeveralConeAndPtCuts)
      if(!IsParticleInJetCone(fTrackVector.Eta(), fTrackVector.Phi(), fLeadingMom.Eta(), fLeadingMom.Phi()) ) continue ;
    
    nTracksInCone++;
    
    TH2F *ha =dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sFFz%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())));
    if(ha) ha->Fill(ptTrig, fTrackVector.Pt()/ptTrig, GetEventWeight());
      
    TH2F *hb  =dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sFFxi%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())));
    if(hb) hb->Fill(ptTrig, TMath::Log(ptTrig/fTrackVector.Pt()), GetEventWeight());
      
    TH2F *hc =dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sFFpt%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())));
    if(hc) hc->Fill(ptTrig, fTrackVector.Pt(), GetEventWeight());
    
  }//track loop
  
  if(nTracksInCone > 0)
  {
    TH2F *hd = dynamic_cast<TH2F*>(GetOutputContainer()->FindObject(Form("%s%sNTracksInCone%s",GetAddedHistogramsStringToName().Data(),type.Data(),lastname.Data())));
    if(hd) hd->Fill(ptTrig, nTracksInCone, GetEventWeight());
  }
}

//________________________________________________________________________
/// Create histograms to be saved in output file and
/// store them in fOutCont.
//________________________________________________________________________
TList *  AliAnaParticleJetLeadingConeCorrelation::GetCreateOutputObjects()
{
  fOutCont = new TList() ;
  fOutCont->SetName("ParticleJetLeadingInConeCorrelationHistograms") ;
  
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins();
  Int_t nphibins = GetHistogramRanges()->GetHistoPhiBins();
  Int_t netabins = GetHistogramRanges()->GetHistoEtaBins();
  Float_t ptmax  = GetHistogramRanges()->GetHistoPtMax();
  Float_t phimax = GetHistogramRanges()->GetHistoPhiMax();
  Float_t etamax = GetHistogramRanges()->GetHistoEtaMax();
  Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin();
  Float_t phimin = GetHistogramRanges()->GetHistoPhiMin();
  Float_t etamin = GetHistogramRanges()->GetHistoEtaMin();
  
  fhChargedLeadingPt  = new TH2F("ChargedLeadingPt","p_{T leading charge} vs p_{T trigger}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
  fhChargedLeadingPt->SetYTitle("p_{T leading charge}");
  fhChargedLeadingPt->SetXTitle("p_{T trigger} (GeV/c)");
  
  fhChargedLeadingPhi  = new TH2F("ChargedLeadingPhi","#phi_{h^{#pm}}  vs p_{T trigger}", nptbins,ptmin,ptmax,nphibins,phimin,phimax);
  fhChargedLeadingPhi->SetYTitle("#phi_{h^{#pm}} (rad)");
  fhChargedLeadingPhi->SetXTitle("p_{T trigger} (GeV/c)");
  
  fhChargedLeadingEta  = new TH2F("ChargedLeadingEta","#eta_{h^{#pm}}  vs p_{T trigger}",nptbins,ptmin,ptmax,netabins,etamin,etamax);
  fhChargedLeadingEta->SetYTitle("#eta_{h^{#pm}} ");
  fhChargedLeadingEta->SetXTitle("p_{T trigger} (GeV/c)");
  
  fhChargedLeadingDeltaPt  = new TH2F("ChargedLeadingDeltaPt","p_{T trigger} - p_{T h^{#pm}} vs p_{T trigger}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
  fhChargedLeadingDeltaPt->SetYTitle("#Delta p_{T} (GeV/c)");
  fhChargedLeadingDeltaPt->SetXTitle("p_{T trigger} (GeV/c)");
  
  fhChargedLeadingDeltaPhi  = new TH2F("ChargedLeadingDeltaPhi","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,TMath::TwoPi());
  fhChargedLeadingDeltaPhi->SetYTitle("#Delta #phi (rad)");
  fhChargedLeadingDeltaPhi->SetXTitle("p_{T trigger} (GeV/c)");
  
  fhChargedLeadingDeltaEta  = new TH2F("ChargedLeadingDeltaEta","#eta_{trigger} - #eta_{h^{#pm}} vs p_{T trigger}",nptbins,ptmin,ptmax,120,-2,2);
  fhChargedLeadingDeltaEta->SetYTitle("#Delta #eta");
  fhChargedLeadingDeltaEta->SetXTitle("p_{T trigger} (GeV/c)");
  
  fhChargedLeadingRatioPt  = new TH2F("ChargedLeadingRatioPt","p_{T leading charge} /p_{T trigger} vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,2);
  fhChargedLeadingRatioPt->SetYTitle("p_{T lead charge} /p_{T trigger}");
  fhChargedLeadingRatioPt->SetXTitle("p_{T trigger} (GeV/c)");
  
  fhChargedLeadingXi  = new TH2F("ChargedLeadingXi","ln(p_{T trigger} / p_{T leading charge} ) vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,10);
  fhChargedLeadingXi->SetYTitle("#xi");
  fhChargedLeadingXi->SetXTitle("p_{T trigger} (GeV/c)");
  
  fOutCont->Add(fhChargedLeadingPt) ;
  fOutCont->Add(fhChargedLeadingPhi) ;
  fOutCont->Add(fhChargedLeadingEta) ;
  fOutCont->Add(fhChargedLeadingDeltaPt) ;
  fOutCont->Add(fhChargedLeadingDeltaPhi) ;
  fOutCont->Add(fhChargedLeadingDeltaEta) ;
  fOutCont->Add(fhChargedLeadingRatioPt) ;
  fOutCont->Add(fhChargedLeadingXi) ;
  
  fhChargedLeadingDeltaPhiRatioPt30  = new TH2F("ChargedLeadingDeltaPhiRatioPt30","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T leading}/p_{T trigger}, charged leading, p_{T trigger} > 30 GeV/c",120,0,TMath::TwoPi(),nptbins,0,1);
  fhChargedLeadingDeltaPhiRatioPt30->SetXTitle("#Delta #phi (rad)");
  fhChargedLeadingDeltaPhiRatioPt30->SetYTitle("p_{T leading} / p_{T trigger}");
  
  fhChargedLeadingDeltaPhiRatioPt50  = new TH2F("ChargedLeadingDeltaPhiRatioPt50","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T leading}/p_{T trigger}, charged leading, p_{T trigger} > 50 GeV/c",120,0,TMath::TwoPi(),nptbins,0,1);
  fhChargedLeadingDeltaPhiRatioPt50->SetXTitle("#Delta #phi (rad)");
  fhChargedLeadingDeltaPhiRatioPt50->SetYTitle("p_{T leading} / p_{T trigger}");
  
  fOutCont->Add(fhChargedLeadingDeltaPhiRatioPt30) ;
  fOutCont->Add(fhChargedLeadingDeltaPhiRatioPt50) ;
  
  if(!fJetsOnlyInCTS){
    
    fhNeutralLeadingPt  = new TH2F("NeutralLeadingPt","p_{T leading #pi^{0}} vs p_{T trigger}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhNeutralLeadingPt->SetYTitle("p_{T leading #pi^{0}}");
    fhNeutralLeadingPt->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhNeutralLeadingPhi  = new TH2F("NeutralLeadingPhi","#phi_{#pi^{0}}  vs p_{T trigger}",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhNeutralLeadingPhi->SetYTitle("#phi_{#pi^{0}} (rad)");
    fhNeutralLeadingPhi->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhNeutralLeadingEta  = new TH2F("NeutralLeadingEta","#eta_{#pi^{0}}  vs p_{T trigger}",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhNeutralLeadingEta->SetYTitle("#eta_{#pi^{0}} ");
    fhNeutralLeadingEta->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhNeutralLeadingDeltaPt  = new TH2F("NeutralLeadingDeltaPt","p_{T trigger} - p_{T #pi^{0}} vs p_{T trigger}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhNeutralLeadingDeltaPt->SetYTitle("#Delta p_{T} (GeV/c)");
    fhNeutralLeadingDeltaPt->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhNeutralLeadingDeltaPhi  = new TH2F("NeutralLeadingDeltaPhi","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,TMath::TwoPi());
    fhNeutralLeadingDeltaPhi->SetYTitle("#Delta #phi (rad)");
    fhNeutralLeadingDeltaPhi->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhNeutralLeadingDeltaEta  = new TH2F("NeutralLeadingDeltaEta","#eta_{trigger} - #eta_{#pi^{0}} vs p_{T trigger}",nptbins,ptmin,ptmax,120,-2,2);
    fhNeutralLeadingDeltaEta->SetYTitle("#Delta #eta");
    fhNeutralLeadingDeltaEta->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhNeutralLeadingRatioPt  = new TH2F("NeutralLeadingRatioPt","p_{T leading #pi^{0}} /p_{T trigger} vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,2);
    fhNeutralLeadingRatioPt->SetYTitle("p_{T lead #pi^{0}} /p_{T trigger}");
    fhNeutralLeadingRatioPt->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhNeutralLeadingXi  = new TH2F("NeutralLeadingXi","ln(p_{T trigger} / p_{T leading #pi^{0}} ) vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,10);
    fhNeutralLeadingXi->SetYTitle("#xi");
    fhNeutralLeadingXi->SetXTitle("p_{T trigger} (GeV/c)");
    
    fOutCont->Add(fhNeutralLeadingPt) ;
    fOutCont->Add(fhNeutralLeadingPhi) ;
    fOutCont->Add(fhNeutralLeadingEta) ;
    fOutCont->Add(fhNeutralLeadingDeltaPt) ;
    fOutCont->Add(fhNeutralLeadingDeltaPhi) ;
    fOutCont->Add(fhNeutralLeadingDeltaEta) ;
    fOutCont->Add(fhNeutralLeadingRatioPt) ;
    fOutCont->Add(fhNeutralLeadingXi) ;
    
    fhNeutralLeadingDeltaPhiRatioPt30  = new TH2F("NeutralLeadingDeltaPhiRatioPt30","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T leading}/p_{T trigger}, neutral leading, p_{T trigger} > 30 GeV/c",120,0,TMath::TwoPi(),nptbins,0,1);
    fhNeutralLeadingDeltaPhiRatioPt30->SetXTitle("#Delta #phi (rad)");
    fhNeutralLeadingDeltaPhiRatioPt30->SetYTitle("p_{T leading} / p_{T trigger}");
    
    fhNeutralLeadingDeltaPhiRatioPt50  = new TH2F("NeutralLeadingDeltaPhiRatioPt50","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T leading}/p_{T trigger}, neutral leading, p_{T trigger} > 50 GeV/c",120,0,TMath::TwoPi(),nptbins,0,1);
    fhNeutralLeadingDeltaPhiRatioPt50->SetXTitle("#Delta #phi (rad)");
    fhNeutralLeadingDeltaPhiRatioPt50->SetYTitle("p_{T leading} / p_{T trigger}");
    fOutCont->Add(fhNeutralLeadingDeltaPhiRatioPt30) ;
    fOutCont->Add(fhNeutralLeadingDeltaPhiRatioPt50) ;
    
  }
  
  if(!fSeveralConeAndPtCuts){// not several cones
    
    //Jet Distributions
    fhJetPt  = new TH2F("JetPt","p_{T  jet} vs p_{T trigger}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhJetPt->SetYTitle("p_{T  jet}");
    fhJetPt->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhJetRatioPt  = new TH2F("JetRatioPt","p_{T  jet}/p_{T trigger} vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,2);
    fhJetRatioPt->SetYTitle("p_{T  jet}/p_{T trigger}");
    fhJetRatioPt->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhJetDeltaPhi  = new TH2F("JetDeltaPhi","#phi_{jet} - #phi_{trigger} vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,TMath::TwoPi());
    fhJetDeltaPhi->SetYTitle("#Delta #phi (rad)");
    fhJetDeltaPhi->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhJetDeltaEta  = new TH2F("JetDeltaEta","#eta_{jet} - #eta_{trigger} vs p_{T trigger}",nptbins,ptmin,ptmax,120,-2,2);
    fhJetDeltaEta->SetYTitle("#Delta #eta");
    fhJetDeltaEta->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhJetLeadingRatioPt  = new TH2F("JetLeadingRatioPt","p_{T  jet} vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,2);
    fhJetLeadingRatioPt->SetYTitle("p_{T  leading}/p_{T jet}");
    fhJetLeadingRatioPt->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhJetLeadingDeltaPhi  = new TH2F("JetLeadingDeltaPhi","#phi_{jet} - #phi_{leading} vs p_{T trigger}",nptbins,ptmin,ptmax,120,-TMath::Pi(),TMath::Pi());
    fhJetLeadingDeltaPhi->SetYTitle("#Delta #phi (rad)");
    fhJetLeadingDeltaPhi->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhJetLeadingDeltaEta  = new TH2F("JetLeadingDeltaEta","#eta_{jet} - #eta_{leading} vs p_{T trigger}",nptbins,ptmin,ptmax,120,-2,2);
    fhJetLeadingDeltaEta->SetYTitle("#Delta #eta");
    fhJetLeadingDeltaEta->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhJetFFz  = new TH2F("JetFFz","z = p_{T i charged}/p_{T trigger} vs p_{T trigger}",nptbins,ptmin,ptmax,200,0.,2);
    fhJetFFz->SetYTitle("z");
    fhJetFFz->SetXTitle("p_{T trigger}");
    
    fhJetFFxi  = new TH2F("JetFFxi","#xi = ln(p_{T trigger}/p_{T i charged}) vs p_{T trigger}",nptbins,ptmin,ptmax,100,0.,10.);
    fhJetFFxi->SetYTitle("#xi");
    fhJetFFxi->SetXTitle("p_{T trigger}");
    
    fhJetFFpt  = new TH2F("JetFFpt","#xi = p_{T i charged}) vs p_{T trigger}",nptbins,ptmin,ptmax,200,0.,50.);
    fhJetFFpt->SetYTitle("p_{T charged hadron}");
    fhJetFFpt->SetXTitle("p_{T trigger}");
    
    fhJetNTracksInCone  = new TH2F("JetNTracksInCone","N particles in cone vs p_{T trigger}",nptbins,ptmin,ptmax,5000,0, 5000);
    fhJetNTracksInCone->SetYTitle("N tracks in jet cone");
    fhJetNTracksInCone->SetXTitle("p_{T trigger} (GeV/c)");
    
    fOutCont->Add(fhJetPt) ;
    fOutCont->Add(fhJetRatioPt) ;
    fOutCont->Add(fhJetDeltaPhi) ;
    fOutCont->Add(fhJetDeltaEta) ;
    fOutCont->Add(fhJetLeadingRatioPt) ;
    fOutCont->Add(fhJetLeadingDeltaPhi) ;
    fOutCont->Add(fhJetLeadingDeltaEta) ;
    fOutCont->Add(fhJetFFz) ;
    fOutCont->Add(fhJetFFxi) ;
    fOutCont->Add(fhJetFFpt) ;
    fOutCont->Add(fhJetNTracksInCone) ;
    
    //Bkg Distributions
    fhBkgPt  = new TH2F("BkgPt","p_{T  bkg} vs p_{T trigger}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhBkgPt->SetYTitle("p_{T  bkg}");
    fhBkgPt->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhBkgRatioPt  = new TH2F("BkgRatioPt","p_{T  bkg}/p_{T trigger} vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,2);
    fhBkgRatioPt->SetYTitle("p_{T  bkg}/p_{T trigger}");
    fhBkgRatioPt->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhBkgDeltaPhi  = new TH2F("BkgDeltaPhi","#phi_{bkg} - #phi_{trigger} vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,TMath::TwoPi());
    fhBkgDeltaPhi->SetYTitle("#Delta #phi (rad)");
    fhBkgDeltaPhi->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhBkgDeltaEta  = new TH2F("BkgDeltaEta","#eta_{bkg} - #eta_{trigger} vs p_{T trigger}",nptbins,ptmin,ptmax,120,-2,2);
    fhBkgDeltaEta->SetYTitle("#Delta #eta");
    fhBkgDeltaEta->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhBkgLeadingRatioPt  = new TH2F("BkgLeadingRatioPt","p_{T  bkg} vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,2);
    fhBkgLeadingRatioPt->SetYTitle("p_{T  leading}/p_{T bkg}");
    fhBkgLeadingRatioPt->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhBkgLeadingDeltaPhi  = new TH2F("BkgLeadingDeltaPhi","#phi_{bkg} - #phi_{leading} vs p_{T trigger}",nptbins,ptmin,ptmax,120,0,TMath::TwoPi());
    fhBkgLeadingDeltaPhi->SetYTitle("#Delta #phi (rad)");
    fhBkgLeadingDeltaPhi->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhBkgLeadingDeltaEta  = new TH2F("BkgLeadingDeltaEta","#eta_{bkg} - #eta_{leading} vs p_{T trigger}",nptbins,ptmin,ptmax,120,-2,2);
    fhBkgLeadingDeltaEta->SetYTitle("#Delta #eta");
    fhBkgLeadingDeltaEta->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhBkgFFz  = new TH2F("BkgFFz","z = p_{T i charged}/p_{T trigger} vs p_{T trigger}", nptbins,ptmin,ptmax,200,0.,2);
    fhBkgFFz->SetYTitle("z");
    fhBkgFFz->SetXTitle("p_{T trigger}");
    
    fhBkgFFxi  = new TH2F("BkgFFxi","#xi = ln(p_{T trigger}/p_{T i charged}) vs p_{T trigger}", nptbins,ptmin,ptmax,100,0.,10.);
    fhBkgFFxi->SetYTitle("#xi");
    fhBkgFFxi->SetXTitle("p_{T trigger}");
    
    fhBkgFFpt  = new TH2F("BkgFFpt","p_{T charged hadron } vs p_{T trigger}", nptbins,ptmin,ptmax,200,0.,50.);
    fhBkgFFpt->SetYTitle("p_{T charged} hadron");
    fhBkgFFpt->SetXTitle("p_{T trigger}");
    
    fhBkgNTracksInCone  = new TH2F("BkgNTracksInCone","N particles in cone vs p_{T trigger}",nptbins,ptmin,ptmax,5000,0, 5000);
    fhBkgNTracksInCone->SetYTitle("N tracks in bkg cone");
    fhBkgNTracksInCone->SetXTitle("p_{T trigger} (GeV/c)");
    
    fOutCont->Add(fhBkgPt) ;
    fOutCont->Add(fhBkgRatioPt) ;
    fOutCont->Add(fhBkgDeltaPhi) ;
    fOutCont->Add(fhBkgDeltaEta) ;
    fOutCont->Add(fhBkgLeadingRatioPt) ;
    fOutCont->Add(fhBkgLeadingDeltaPhi) ;
    fOutCont->Add(fhBkgLeadingDeltaEta) ;
    fOutCont->Add(fhBkgFFz) ;
    fOutCont->Add(fhBkgFFxi) ;
    fOutCont->Add(fhBkgFFpt) ;
    fOutCont->Add(fhBkgNTracksInCone) ;
    
  }//not several cones
  else{ //If we want to study the jet for different cones and pt
    for(Int_t icone = 0; icone<fJetNCone; icone++){//icone
      for(Int_t ipt = 0; ipt<fJetNPt;ipt++){ //ipt
        
        TString lastnamehist ="Cone"+ fJetNameCones[icone]+"Pt"+ fJetNamePtThres[ipt];
        TString lastnametitle =", cone ="+fJetNameCones[icone]+", pt > " +fJetNamePtThres[ipt]+" GeV/c";
        
        //Jet Distributions
        fhJetPts[icone][ipt] = new TH2F(Form("JetPt%s",lastnamehist.Data()),Form("p_{T  jet} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhJetPts[icone][ipt]->SetYTitle("p_{T  jet}");
        fhJetPts[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhJetRatioPts[icone][ipt] = new TH2F(Form("JetRatioPt%s",lastnamehist.Data()),Form("p_{T  jet}/p_{T trigger} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,120,0,2);
        fhJetRatioPts[icone][ipt]->SetYTitle("p_{T  jet}/p_{T trigger}");
        fhJetRatioPts[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhJetDeltaPhis[icone][ipt] = new TH2F(Form("JetDeltaPhi%s",lastnamehist.Data()),Form("#phi_{jet} - #phi_{trigger} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,120,0,TMath::TwoPi());
        fhJetDeltaPhis[icone][ipt]->SetYTitle("#Delta #phi (rad)");
        fhJetDeltaPhis[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhJetDeltaEtas[icone][ipt] = new TH2F(Form("JetDeltaEta%s",lastnamehist.Data()),Form("#eta_{jet} - #eta_{trigger} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,120,-2,2);
        fhJetDeltaEtas[icone][ipt]->SetYTitle("#Delta #eta");
        fhJetDeltaEtas[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhJetLeadingRatioPts[icone][ipt] = new TH2F(Form("JetLeadingRatioPt%s",lastnamehist.Data()),Form("p_{T  jet} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,120,0,2);
        fhJetLeadingRatioPts[icone][ipt]->SetYTitle("p_{T  leading}/p_{T jet}");
        fhJetLeadingRatioPts[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhJetLeadingDeltaPhis[icone][ipt] = new TH2F(Form("JetLeadingDeltaPhi%s",lastnamehist.Data()),Form("#phi_{jet} - #phi_{leading} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,120,0,TMath::TwoPi());
        fhJetLeadingDeltaPhis[icone][ipt]->SetYTitle("#Delta #phi (rad)");
        fhJetLeadingDeltaPhis[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhJetLeadingDeltaEtas[icone][ipt] = new TH2F(Form("JetLeadingDeltaEta%s",lastnamehist.Data()),Form("#eta_{jet} - #eta_{leading} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,120,-2,2);
        fhJetLeadingDeltaEtas[icone][ipt]->SetYTitle("#Delta #eta");
        fhJetLeadingDeltaEtas[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhJetFFzs[icone][ipt] = new TH2F(Form("JetFFz%s",lastnamehist.Data()),"z = p_{T i charged}/p_{T trigger} vs p_{T trigger}", 120,0.,120.,200,0.,2);
        fhJetFFzs[icone][ipt]->SetYTitle("z");
        fhJetFFzs[icone][ipt]->SetXTitle("p_{T trigger}");
        
        fhJetFFxis[icone][ipt] = new TH2F(Form("JetFFxi%s",lastnamehist.Data()),"#xi = ln(p_{T trigger}/p_{T i charged}) vs p_{T trigger}", 120,0.,120.,100,0.,10.);
        fhJetFFxis[icone][ipt]->SetYTitle("#xi");
        fhJetFFxis[icone][ipt]->SetXTitle("p_{T trigger}");
        
        fhJetFFpts[icone][ipt] = new TH2F(Form("JetFFpt%s",lastnamehist.Data()),"p_{T charged hadron } in jet vs p_{T trigger}", 120,0.,120.,200,0.,50.);
        fhJetFFpts[icone][ipt]->SetYTitle("p_{T charged hadron}");
        fhJetFFpts[icone][ipt]->SetXTitle("p_{T trigger}");
        
        fhJetNTracksInCones[icone][ipt] = new TH2F(Form("JetNTracksInCone%s",lastnamehist.Data()),Form("N particles in cone vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,5000,0, 5000);
        fhJetNTracksInCones[icone][ipt]->SetYTitle("N tracks in jet cone");
        fhJetNTracksInCones[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fOutCont->Add(fhJetPts[icone][ipt]) ;
        fOutCont->Add(fhJetRatioPts[icone][ipt]) ;
        fOutCont->Add(fhJetDeltaPhis[icone][ipt]) ;
        fOutCont->Add(fhJetDeltaEtas[icone][ipt]) ;
        fOutCont->Add(fhJetLeadingRatioPts[icone][ipt]) ;
        fOutCont->Add(fhJetLeadingDeltaPhis[icone][ipt]) ;
        fOutCont->Add(fhJetLeadingDeltaEtas[icone][ipt]) ;
        fOutCont->Add(fhJetFFzs[icone][ipt]) ;
        fOutCont->Add(fhJetFFxis[icone][ipt]) ;
        fOutCont->Add(fhJetFFpts[icone][ipt]) ;
        fOutCont->Add(fhJetNTracksInCones[icone][ipt]) ;
        
        //Bkg Distributions
        fhBkgPts[icone][ipt] = new TH2F(Form("BkgPt%s",lastnamehist.Data()),Form("p_{T  bkg} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhBkgPts[icone][ipt]->SetYTitle("p_{T  bkg}");
        fhBkgPts[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhBkgRatioPts[icone][ipt] = new TH2F(Form("BkgRatioPt%s",lastnamehist.Data()),Form("p_{T  bkg}/p_{T trigger} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,120,0,2);
        fhBkgRatioPts[icone][ipt]->SetYTitle("p_{T  bkg}/p_{T trigger}");
        fhBkgRatioPts[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhBkgDeltaPhis[icone][ipt] = new TH2F(Form("BkgDeltaPhi%s",lastnamehist.Data()),Form("#phi_{bkg} - #phi_{trigger} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,120,0,TMath::TwoPi());
        fhBkgDeltaPhis[icone][ipt]->SetYTitle("#Delta #phi (rad)");
        fhBkgDeltaPhis[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhBkgDeltaEtas[icone][ipt] = new TH2F(Form("BkgDeltaEta%s",lastnamehist.Data()),Form("#eta_{bkg} - #eta_{trigger} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,120,-2,2);
        fhBkgDeltaEtas[icone][ipt]->SetYTitle("#Delta #eta");
        fhBkgDeltaEtas[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhBkgLeadingRatioPts[icone][ipt] = new TH2F(Form("BkgLeadingRatioPt%s",lastnamehist.Data()),Form("p_{T  bkg} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,120,0,2);
        fhBkgLeadingRatioPts[icone][ipt]->SetYTitle("p_{T  leading}/p_{T bkg}");
        fhBkgLeadingRatioPts[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhBkgLeadingDeltaPhis[icone][ipt] = new TH2F(Form("BkgLeadingDeltaPhi%s",lastnamehist.Data()),Form("#phi_{bkg} - #phi_{leading} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,120,0,TMath::TwoPi());
        fhBkgLeadingDeltaPhis[icone][ipt]->SetYTitle("#Delta #phi (rad)");
        fhBkgLeadingDeltaPhis[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhBkgLeadingDeltaEtas[icone][ipt] = new TH2F(Form("BkgLeadingDeltaEta%s",lastnamehist.Data()),Form("#eta_{bkg} - #eta_{leading} vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,120,-2,2);
        fhBkgLeadingDeltaEtas[icone][ipt]->SetYTitle("#Delta #eta");
        fhBkgLeadingDeltaEtas[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhBkgFFzs[icone][ipt] = new TH2F(Form("BkgFFz%s",lastnamehist.Data()),"z = p_{T i charged}/p_{T trigger} vs p_{T trigger}", 120,0.,120.,200,0.,2);
        fhBkgFFzs[icone][ipt]->SetYTitle("z");
        fhBkgFFzs[icone][ipt]->SetXTitle("p_{T trigger}");
        
        fhBkgFFxis[icone][ipt] = new TH2F(Form("BkgFFxi%s",lastnamehist.Data()),"#xi = ln(p_{T trigger}/p_{T i charged}) vs p_{T trigger}", 120,0.,120.,100,0.,10.);
        fhBkgFFxis[icone][ipt]->SetYTitle("#xi");
        fhBkgFFxis[icone][ipt]->SetXTitle("p_{T trigger}");
        
        fhBkgFFpts[icone][ipt] = new TH2F(Form("BkgFFpt%s",lastnamehist.Data()),"p_{T charged hadron} in jet vs p_{T trigger}", 120,0.,120.,200,0.,50.);
        fhBkgFFpts[icone][ipt]->SetYTitle("p_{T charged hadron}");
        fhBkgFFpts[icone][ipt]->SetXTitle("p_{T trigger}");
        
        fhBkgNTracksInCones[icone][ipt] = new TH2F(Form("BkgNTracksInCone%s",lastnamehist.Data()),Form("N particles in cone vs p_{T trigger} %s",lastnametitle.Data()),nptbins,ptmin,ptmax,5000,0, 5000);
        fhBkgNTracksInCones[icone][ipt]->SetYTitle("N tracks in bkg cone");
        fhBkgNTracksInCones[icone][ipt]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fOutCont->Add(fhBkgPts[icone][ipt]) ;
        fOutCont->Add(fhBkgRatioPts[icone][ipt]) ;
        fOutCont->Add(fhBkgDeltaPhis[icone][ipt]) ;
        fOutCont->Add(fhBkgDeltaEtas[icone][ipt]) ;
        fOutCont->Add(fhBkgLeadingRatioPts[icone][ipt]) ;
        fOutCont->Add(fhBkgLeadingDeltaPhis[icone][ipt]) ;
        fOutCont->Add(fhBkgLeadingDeltaEtas[icone][ipt]) ;
        fOutCont->Add(fhBkgFFzs[icone][ipt]) ;
        fOutCont->Add(fhBkgFFxis[icone][ipt]) ;
        fOutCont->Add(fhBkgFFpts[icone][ipt]) ;
        fOutCont->Add(fhBkgNTracksInCones[icone][ipt]) ;
        
      }//ipt
    } //icone
  }//If we want to study any cone or pt threshold
  
  //Keep neutral meson selection histograms if requiered
  //Setting done in AliNeutralMesonSelection
  if(GetNeutralMesonSelection()){
    TList * nmsHistos = GetNeutralMesonSelection()->GetCreateOutputObjects() ;
    if(GetNeutralMesonSelection()->AreNeutralMesonSelectionHistosKept())
      for(Int_t i = 0; i < nmsHistos->GetEntries(); i++) fOutCont->Add(nmsHistos->At(i)) ;
    delete nmsHistos;
  }
  
  //  if(GetDebug() > 2)
  //  {
  //    printf("AliAnaParticleJetLeadingConeCorrelation::GetCreateOutputObjects() - All histograms names : \n");
  //    for(Int_t i  = 0 ;  i<  fOutCont->GetEntries(); i++)
  //      printf("Histo i %d name %s\n",i,((fOutCont->At(i))->GetName()));
  //    //cout<< (fOutCont->At(i))->GetName()<<endl;
  //  }
  
  return fOutCont;
}

//__________________________________________________________________________________________________________
/// Search Charged or Neutral leading particle, select the highest one and fill AOD.
//__________________________________________________________________________________________________________
Bool_t  AliAnaParticleJetLeadingConeCorrelation::GetLeadingParticle(AliCaloTrackParticleCorrelation * particle)
{
  GetLeadingCharge(particle) ;
  if(!fJetsOnlyInCTS) GetLeadingPi0(particle) ;
  
  Double_t ptch = fLeadingChargeMom.Pt();
  Double_t ptpi = fLeadingPi0Mom   .Pt();
  
  if (ptch > 0 || ptpi > 0)
  {
    if((ptch >= ptpi))
    {
      AliDebug(1,"Leading found in CTS");
      
      fLeadingMom = fLeadingChargeMom;
      
      AliDebug(1,Form("Found Leading: pt %2.3f, phi %2.3f deg, eta %2.3f",fLeadingMom.Pt(),fLeadingMom.Phi()*TMath::RadToDeg(),fLeadingMom.Eta())) ;
      
      //Put leading in AOD
      particle->SetLeading(fLeadingChargeMom);
      particle->SetLeadingDetector(kCTS);
      return kTRUE;
    }
    else
    {
      AliDebug(1,"Leading found in EMCAL");
      
      fLeadingMom = fLeadingPi0Mom;
      
      AliDebug(1,Form("Found Leading: pt %2.3f, phi %2.3f, eta %2.3f",fLeadingMom.Pt(),fLeadingMom.Phi()*TMath::RadToDeg(),fLeadingMom.Eta())) ;
      //Put leading in AOD
      particle->SetLeading(fLeadingPi0Mom);
      particle->SetLeadingDetector(kEMCAL);
      return kTRUE;
    }
  }
  
  AliDebug(1,"NO LEADING PARTICLE FOUND");
  
  return kFALSE;
}

//_______________________________________________________________________________________________________
/// Search for the charged particle with highest pt and with
/// Phi=Phi_trigger-Pi and pT=0.1E_gamma.
//_______________________________________________________________________________________________________
void  AliAnaParticleJetLeadingConeCorrelation::GetLeadingCharge(AliCaloTrackParticleCorrelation * particle)
{
  if(!GetCTSTracks()) return;
  
  Double_t ptTrig  = particle->Pt();
  Double_t phiTrig = particle->Phi();
  Double_t rat     = -100 ;
  Double_t ptl     = -100 ;
  Double_t phil    = -100 ;
  Double_t pt      = -100.;
  Double_t phi     = -100.;
  
  for(Int_t ipr = 0;ipr < GetCTSTracks()->GetEntriesFast() ; ipr ++ )
  {
    AliVTrack* track = (AliVTrack *)(GetCTSTracks()->At(ipr)) ;
    fTrackVector.SetXYZ(track->Px(),track->Py(),track->Pz());
    pt   = fTrackVector.Pt();
    phi  = fTrackVector.Phi() ;
    if(phi < 0) phi+=TMath::TwoPi();
    rat  = pt/ptTrig ;
    
    //printf("AliAnaParticleJetLeadingConeCorrelation::GetLeadingCharge() - Tracks: pt %2.3f eta %2.3f phi %2.3f pt/ptTrig %2.3f \n",
    //	   pt, fTrackVector.Eta(), phi,pt/ptTrig) ;
    
    Float_t deltaphi = TMath::Abs(phiTrig-phi);
    if((deltaphi > fDeltaPhiMinCut)     && (deltaphi < fDeltaPhiMaxCut) &&
       (rat      > fLeadingRatioMinCut) && (rat      < fLeadingRatioMaxCut)  && (pt  > ptl))
    {
      phil = phi ;
      ptl  = pt ;
      fLeadingChargeMom.SetVect(fTrackVector);
    }
  }// track loop
  
  if( ptl > 0 )AliDebug(1,Form("Leading in CTS: pt %2.3f eta %2.3f phi %2.3f pt/ptTrig %2.3f, |phiTrig-phi| %2.3f",
                               ptl, fLeadingChargeMom.Eta(), phil,ptl/ptTrig, TMath::Abs(phiTrig-phil))) ;
}

//____________________________________________________________________________________________________
/// Search for the neutral pion with highest pt and with
/// Phi=Phi_trigger-Pi and pT=0.1E_gamma.
//____________________________________________________________________________________________________
void  AliAnaParticleJetLeadingConeCorrelation::GetLeadingPi0(AliCaloTrackParticleCorrelation * particle)
{
  if(!GetEMCALClusters()) return ;
  
  Double_t ptTrig  = particle->Pt();
  Double_t phiTrig = particle->Phi();
  Double_t rat     = -100 ;
  Double_t ptl     = -100 ;
  Double_t phil    = -100 ;
  Double_t pt      = -100.;
  Double_t phi     = -100.;
  
  //Get vertex for photon momentum calculation
  Double_t vertex [] = {0,0,0} ; //vertex
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
  {
    GetVertex(vertex);
  }
  
  //Cluster loop, select pairs with good pt, phi and fill AODs or histograms
  for(Int_t iclus = 0;iclus < GetEMCALClusters()->GetEntriesFast() ; iclus ++ )
  {
    AliVCluster * calo = (AliVCluster *)(GetEMCALClusters()->At(iclus)) ;
    
    //Cluster selection, not charged, with photon or pi0 id and in fiducial cut
    Int_t pdgi=0;
    if(!SelectCluster(calo, vertex,  fLeadingPhoMom1, pdgi))  continue ;
    
    AliDebug(1,Form("Neutral cluster: pt %2.3f, phi %2.3f",fLeadingPhoMom1.Pt(),fLeadingPhoMom1.Phi()));
    
    //2 gamma overlapped, found with PID
    if(pdgi == AliCaloPID::kPi0)
    {
      AliDebug(1,"Neutral cluster ID as Pi0");
      
      pt  = fLeadingPhoMom1.Pt();
      rat = pt/ptTrig;
      phi = fLeadingPhoMom1.Phi();
      if(phi < 0) phi+=TMath::TwoPi();
      
      //Selection within angular and energy limits
      Float_t deltaphi = TMath::Abs(phiTrig-phi);
      if(pt > ptl  && rat > fLeadingRatioMinCut  && rat < fLeadingRatioMaxCut  &&
         deltaphi > fDeltaPhiMinCut && deltaphi < fDeltaPhiMaxCut )
      {
        phil = phi ;
        ptl  = pt ;
        fLeadingPi0Mom.SetPxPyPzE(fLeadingPhoMom1.Px(),fLeadingPhoMom1.Py(),fLeadingPhoMom1.Pz(),fLeadingPhoMom1.E());
      }// cuts
    }// pdg = AliCaloPID::kPi0
    //Make invariant mass analysis
    else if(pdgi == AliCaloPID::kPhoton)
    {
      //Search the photon companion in case it comes from  a Pi0 decay
      //Apply several cuts to select the good pair
      for(Int_t jclus = iclus+1; jclus < GetEMCALClusters()->GetEntriesFast() ; jclus++ )
      {
        AliVCluster * calo2 = (AliVCluster *) (GetEMCALClusters()->At(jclus)) ;
        
        //Cluster selection, not charged with photon or pi0 id and in fiducial cut
        Int_t pdgj=0;
        
        if     (!SelectCluster(calo2, vertex,  fLeadingPhoMom2, pdgj))  continue ;
        
        if(pdgj == AliCaloPID::kPhoton )
        {
          pt  = (fLeadingPhoMom1+fLeadingPhoMom2).Pt();
          phi = (fLeadingPhoMom1+fLeadingPhoMom2).Phi();
          if(phi < 0) phi+=TMath::TwoPi();
          rat = pt/ptTrig;
          
          //Selection within angular and energy limits
          Float_t deltaphi = TMath::Abs(phiTrig-phi);
          
          AliDebug(1,Form("Neutral Hadron Correlation: gamma pair: pt %2.2f, phi %2.2f, eta %2.2f, |phiTrig-phi| %2.3f, pt/ptTrig %2.3f, M %2.3f",
                          pt,phi,(fLeadingPhoMom1+fLeadingPhoMom2).Eta(), deltaphi, rat, (fLeadingPhoMom1+fLeadingPhoMom2).M()));
          
          if(pt > ptl  && rat > fLeadingRatioMinCut  && rat < fLeadingRatioMaxCut  &&
             deltaphi > fDeltaPhiMinCut && deltaphi < fDeltaPhiMaxCut ){
            //Select good pair (aperture and invariant mass)
            if(GetNeutralMesonSelection()->SelectPair(fLeadingPhoMom1, fLeadingPhoMom2,kEMCAL)){
              phil = phi ;
              ptl  = pt ;
              fLeadingPi0Mom=(fLeadingPhoMom1+fLeadingPhoMom2);
              
              AliDebug(1,Form("Neutral Hadron Correlation: Selected gamma pair: pt %2.2f, phi %2.2f, eta %2.2f, M %2.3f",
                              ptl,phil,(fLeadingPhoMom1+fLeadingPhoMom2).Eta(), (fLeadingPhoMom1+fLeadingPhoMom2).M()));
            }//pi0 selection
          }//Pair selected as leading
        }//if pair of gammas
      }//2nd loop
    }// if pdg = 22
  }// 1st Loop
  
  if(fLeadingPi0Mom.Pt() > 0 )
    AliDebug(1,Form("Leading EMCAL: pt %2.3f eta %2.3f phi %2.3f pt/Eg %2.3f",
                    fLeadingPi0Mom.Pt(), fLeadingPi0Mom.Eta(), phil,  fLeadingPi0Mom.Pt()/ptTrig)) ;
}

//____________________________________________________________________________
/// Initialize the parameters of the analysis.
//____________________________________________________________________________
void AliAnaParticleJetLeadingConeCorrelation::InitParameters()
{
  SetInputAODName("CaloTrackParticle");
  SetAODObjArrayName("JetLeadingCone");
  AddToHistogramsName("AnaJetLCCorr_");
  
  fJetsOnlyInCTS      = kFALSE ;
  fPbPb               = kFALSE ;
  fReMakeJet          = kFALSE ;
  
  //Leading selection parameters
  fDeltaPhiMinCut     = 2.9 ;
  fDeltaPhiMaxCut     = 3.4 ;
  fLeadingRatioMinCut = 0.1;
  fLeadingRatioMaxCut = 1.5;
  
  //Jet selection parameters
  //Fixed cut
  fJetRatioMaxCut     = 1.2 ;
  fJetRatioMinCut     = 0.3 ;
  fJetCTSRatioMaxCut  = 1.2 ;
  fJetCTSRatioMinCut  = 0.3 ;
  fSelect               = 0  ; //0, Accept all jets, 1, selection depends on energy, 2 fixed selection
  
  fSelectIsolated = kFALSE;
  
  //Cut depending on gamma energy
  fPtTriggerSelectionCut = 10.; //For Low pt jets+BKG, another limits applied
  //Reconstructed jet energy dependence parameters
  //e_jet = a1+e_gamma b2.
  //Index 0-> Pt>2 GeV r = 0.3; Index 1-> Pt>0.5 GeV r = 0.3
  fJetE1[0] = -5.75; fJetE1[1] = -4.1;
  fJetE2[0] = 1.005; fJetE2[1] = 1.05;
  
  //Reconstructed sigma of jet energy dependence parameters
  //s_jet = a1+e_gamma b2.
  //Index 0-> Pt>2 GeV r = 0.3; Index 1-> Pt>0.5 GeV r = 0.3
  fJetSigma1[0] = 2.65;   fJetSigma1[1] = 2.75;
  fJetSigma2[0] = 0.0018; fJetSigma2[1] = 0.033;
  
  //Background mean energy and RMS
  //Index 0-> No BKG; Index 1-> BKG > 2 GeV;
  //Index 2-> (low pt jets)BKG > 0.5 GeV;
  //Index > 2, same for CTS conf
  fBkgMean[0] = 0.; fBkgMean[1] = 8.8 ; fBkgMean[2] = 69.5;
  fBkgMean[3] = 0.; fBkgMean[4] = 6.4;  fBkgMean[5] = 48.6;
  fBkgRMS[0]  = 0.; fBkgRMS[1]  = 7.5;  fBkgRMS[2]  = 22.0;
  fBkgRMS[3]  = 0.; fBkgRMS[4]  = 5.4;  fBkgRMS[5]  = 13.2;
  
  //Factor x of min/max = E -+ x * sigma. Obtained after selecting the
  //limits for monoenergetic jets.
  //Index 0-> No BKG; Index 1-> BKG > 2 GeV;
  //Index 2-> (low pt jets) BKG > 0.5 GeV;
  //Index > 2, same for CTS conf
  
  fJetXMin1[0] =-0.69 ; fJetXMin1[1] = 0.39 ; fJetXMin1[2] =-0.88 ;
  fJetXMin1[3] =-2.0  ; fJetXMin1[4] =-0.442 ; fJetXMin1[5] =-1.1  ;
  fJetXMin2[0] = 0.066; fJetXMin2[1] = 0.038; fJetXMin2[2] = 0.034;
  fJetXMin2[3] = 0.25 ; fJetXMin2[4] = 0.113; fJetXMin2[5] = 0.077 ;
  fJetXMax1[0] =-3.8  ; fJetXMax1[1] =-0.76 ; fJetXMax1[2] =-3.6  ;
  fJetXMax1[3] =-2.7  ; fJetXMax1[4] =-1.21 ; fJetXMax1[5] =-3.7  ;
  fJetXMax2[0] =-0.012; fJetXMax2[1] =-0.022; fJetXMax2[2] = 0.016;
  fJetXMax2[3] =-0.024; fJetXMax2[4] =-0.008; fJetXMax2[5] = 0.027;
  
  
  //Different cones and pt thresholds to construct the jet
  
  fJetCone        = 0.3  ;
  fJetPtThreshold = 0.5   ;
  fJetPtThresPbPb = 2.   ;
  fJetNCone       = 4    ;
  fJetNPt         = 4    ;
  fJetCones[0]    = 0.2  ; fJetNameCones[0]   = "02" ;
  fJetCones[1]    = 0.3  ; fJetNameCones[1]   = "03" ;
  fJetCones[2]    = 0.4  ; fJetNameCones[2]   = "04" ;
  fJetCones[2]    = 0.5  ; fJetNameCones[2]   = "05" ;
  
  fJetPtThres[0]  = 0.0  ; fJetNamePtThres[0] = "00" ;
  fJetPtThres[1]  = 0.5  ; fJetNamePtThres[1] = "05" ;
  fJetPtThres[2]  = 1.0  ; fJetNamePtThres[2] = "10" ;
  fJetPtThres[3]  = 2.0  ; fJetNamePtThres[3] = "20" ;
}

//__________________________________________________________________________________________________
/// Given the pt of the jet and the trigger particle, select the jet or not
/// 3 options, fSelect=0 accepts all, fSelect=1 selects jets depending on a
/// function energy dependent and fSelect=2 selects on simple fixed cuts.
//__________________________________________________________________________________________________
Bool_t AliAnaParticleJetLeadingConeCorrelation::IsJetSelected(Double_t ptTrig, Double_t ptjet) const
{
  if(ptjet == 0) return kFALSE;
  
  Double_t rat = ptTrig / ptjet ;
  
  //###############################################################
  if(fSelect == 0)
    return kTRUE; //Accept all jets, no restriction
  //###############################################################
  else if(fSelect == 1)
  {
    //Check if the energy of the reconstructed jet is within an energy window
    //WARNING: to be rechecked, don't remember what all the steps mean
    Double_t par[6];
    Double_t xmax[2];
    Double_t xmin[2];
    
    Int_t iCTS = 0;
    if(fJetsOnlyInCTS)
      iCTS = 3 ;
    
    if(!fPbPb){
      //Phythia alone, jets with pt_th > 0.2, r = 0.3
      par[0] = fJetE1[0]; par[1] = fJetE2[0];
      //Energy of the jet peak
      //e_jet = fJetE1[0]+fJetE2[0]*e_gamma, simulation fit
      par[2] = fJetSigma1[0]; par[3] = fJetSigma2[0];
      //Sigma  of the jet peak
      //sigma_jet = fJetSigma1[0]+fJetSigma2[0]*e_gamma, simulation fit
      par[4] = fBkgMean[0 + iCTS]; par[5] = fBkgRMS[0 + iCTS];
      //Parameters reserved for PbPb bkg.
      xmax[0] = fJetXMax1[0 + iCTS]; xmax[1] = fJetXMax2[0 + iCTS];
      xmin[0] = fJetXMin1[0 + iCTS]; xmin[1] = fJetXMin2[0 + iCTS];
      //Factor that multiplies sigma to obtain the best limits,
      //by observation, of mono jet ratios (ptjet/ptTrig)
      //X_jet = fJetX1[0]+fJetX2[0]*e_gamma
      
    }
    else
    {
      if(ptTrig > fPtTriggerSelectionCut)
      {
        //Phythia +PbPb with  pt_th > 2 GeV/c, r = 0.3
        par[0] = fJetE1[0]; par[1] = fJetE2[0];
        //Energy of the jet peak, same as in pp
        //e_jet = fJetE1[0]+fJetE2[0]*e_gamma, simulation fit
        par[2] = fJetSigma1[0]; par[3] = fJetSigma2[0];
        //Sigma  of the jet peak, same as in pp
        //sigma_jet = fJetSigma1[0]+fJetSigma2[0]*e_gamma, simulation fit
        par[4] = fBkgMean[1 + iCTS]; par[5] = fBkgRMS[1 + iCTS];
        //Mean value and RMS of PbPb Bkg
        xmax[0] = fJetXMax1[1 + iCTS]; xmax[1] = fJetXMax2[1 + iCTS];
        xmin[0] = fJetXMin1[1 + iCTS]; xmin[1] = fJetXMin2[1 + iCTS];
        //Factor that multiplies sigma to obtain the best limits,
        //by observation, of mono jet ratios (ptjet/ptTrig) mixed with PbPb Bkg,
        //pt_th > 2 GeV, r = 0.3
        //X_jet = fJetX1[0]+fJetX2[0]*e_gamma
      }
      else
      {
        //Phythia + PbPb with  pt_th > 0.5 GeV/c, r = 0.3
        par[0] = fJetE1[1]; par[1] = fJetE2[1];
        //Energy of the jet peak, pt_th > 2 GeV/c, r = 0.3
        //e_jet = fJetE1[0]+fJetE2[0]*e_gamma, simulation fit
        par[2] = fJetSigma1[1]; par[3] = fJetSigma2[1];
        //Sigma  of the jet peak, pt_th > 2 GeV/c, r = 0.3
        //sigma_jet = fJetSigma1[0]+fJetSigma2[0]*e_gamma, simulation fit
        par[4] = fBkgMean[2 + iCTS]; par[5] = fBkgRMS[2 + iCTS];
        //Mean value and RMS of PbPb Bkg in a 0.3 cone, pt > 2 GeV.
        xmax[0] = fJetXMax1[2 + iCTS]; xmax[1] = fJetXMax2[2 + iCTS];
        xmin[0] = fJetXMin1[2 + iCTS]; xmin[1] = fJetXMin2[2 + iCTS];
        //Factor that multiplies sigma to obtain the best limits,
        //by observation, of mono jet ratios (ptjet/ptTrig) mixed with PbPb Bkg,
        //pt_th > 2 GeV, r = 0.3
        //X_jet = fJetX1[0]+fJetX2[0]*e_gamma
      }//If low pt jet in bkg
    }//if Bkg
    
    //Calculate minimum and maximum limits of the jet ratio.
    Double_t min = CalculateJetRatioLimit(ptTrig, par, xmin);
    Double_t max = CalculateJetRatioLimit(ptTrig, par, xmax);
    
    AliDebug(3,Form("Jet selection? : Limits min %2.3f, max %2.3f,  pt_jet %2.3f,  pt_gamma %2.3f, pt_jet / pt_gamma %2.3f",min,max,ptjet,ptTrig,rat));
    
    if(( min < rat ) && ( max > ptjet/rat))
      return kTRUE;
    else
      return kFALSE;
  }//fSelect = 1
  //###############################################################
  else if(fSelect == 2)
  {
    //Simple selection
    if(!fJetsOnlyInCTS)
    {
      if((rat <  fJetRatioMaxCut) && (rat > fJetRatioMinCut )) return kTRUE;
    }
    else
    {
      if((rat <  fJetCTSRatioMaxCut) && (rat > fJetCTSRatioMinCut )) return kTRUE;
    }
  }// fSelect = 2
  //###############################################################
  else
  {
    AliWarning(")Jet selection option larger than 2, DON'T SELECT JETS");
    return kFALSE ;
  }
  
  return kFALSE;
}

//___________________________________________________________________
/// Check if the particle is inside the cone defined by the leading particle.
/// WARNING: To be rechecked.
//___________________________________________________________________
Bool_t AliAnaParticleJetLeadingConeCorrelation::IsParticleInJetCone(Double_t eta , Double_t phi,
                                                                    Double_t etal, Double_t phil) const
{
  if(phi < 0) phi+=TMath::TwoPi();
  if(phil < 0) phil+=TMath::TwoPi();
  Double_t  rad = 10000 + fJetCone;
  
  if(TMath::Abs(phi-phil) <= (TMath::TwoPi() - fJetCone))
    rad = TMath::Sqrt(TMath::Power(eta-etal,2)+TMath::Power(phi-phil,2));
  else{
    if(phi-phil > TMath::TwoPi() - fJetCone)
      rad = TMath::Sqrt(TMath::Power(eta-etal,2)+TMath::Power((phi-TMath::TwoPi())-phil,2));
    if(phi-phil < -(TMath::TwoPi() - fJetCone))
      rad = TMath::Sqrt(TMath::Power(eta-etal,2)+TMath::Power((phi+TMath::TwoPi())-phil,2));
  }
  
  if(rad < fJetCone) return kTRUE ;
  else               return kFALSE ;
}

//__________________________________________________________________
// Particle-Jet/Hadron Correlation Analysis, fill AODs.
//__________________________________________________________________
void  AliAnaParticleJetLeadingConeCorrelation::MakeAnalysisFillAOD()
{
  if(!GetInputAODBranch())
  {
    AliFatal(Form("No input particles in AOD with name branch < %s >",GetInputAODName().Data()));
    return;// coverity
  }
  
  if(strcmp(GetInputAODBranch()->GetClass()->GetName(), "AliCaloTrackParticleCorrelation"))
    AliFatal(Form("Wrong type of AOD object, change AOD class name in input AOD: It should be <AliCaloTrackParticleCorrelation> and not <%s>",
                  GetInputAODBranch()->GetClass()->GetName()));
  
  
  AliDebug(1,"Begin jet leading cone  correlation analysis, fill AODs");
  AliDebug(1,Form("In particle branch aod entries %d", GetInputAODBranch()->GetEntriesFast()));
  AliDebug(1,Form("In CTS aod entries %d",             GetCTSTracks()     ->GetEntriesFast()));
  AliDebug(1,Form("In EMCAL aod entries %d",           GetEMCALClusters() ->GetEntriesFast()));
  
  //Loop on stored AOD particles, trigger
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliCaloTrackParticleCorrelation* particle =  (AliCaloTrackParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    //  printf("AliAnaParticleJetLeadingConeCorrelation::MakeAnalysisFillAOD() - Trigger : pt %3.2f, phi %2.2f, eta %2.2f\n",particle->Pt(), particle->Phi(), particle->Eta());
    
    //Search leading particles in CTS and EMCAL
    if(GetLeadingParticle(particle))
    {
      //Construct the jet around the leading, Fill AOD jet particle list, select jet
      //and fill AOD with jet and background
      MakeAODJet(particle);
      
    }//Leading found
  }//AOD trigger particle loop
  
  AliDebug(1,"End of jet leading cone analysis, fill AODs");
}

//_________________________________________________________________________
/// Particle-Jet/Hadron Correlation Analysis, fill histograms.
//_________________________________________________________________________
void  AliAnaParticleJetLeadingConeCorrelation::MakeAnalysisFillHistograms()
{
  if(!GetInputAODBranch())
  {
    printf("AliAnaParticleJetLeadingConeCorrelation::MakeAnalysisFillHistograms() - No input particles in AOD with name branch < %s >",
           GetInputAODName().Data());
    return;
  }
  
  AliDebug(1,"Begin jet leading cone  correlation analysis, fill histograms");
  AliDebug(1,Form("In particle branch aod entries %d", GetInputAODBranch()->GetEntriesFast()));
  AliDebug(1,Form("In CTS aod entries %d",             GetCTSTracks()     ->GetEntriesFast()));
  AliDebug(1,Form("In EMCAL aod entries %d",           GetEMCALClusters() ->GetEntriesFast()));
  
  //Loop on stored AOD particles, trigger
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliCaloTrackParticleCorrelation* particle =  (AliCaloTrackParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    if(OnlyIsolated() && !particle->IsIsolated()) continue;
    
    Double_t pt  = particle->Pt();
    Double_t phi = particle->Phi();
    Double_t eta = particle->Eta();
    
    //Get leading particle, fill histograms
    fLeadingMom = particle->GetLeading();
    Int_t det   = particle->GetLeadingDetector();
    
    if(det > 0 && fLeadingMom.Pt() > 0)
    {
      Double_t ptL = fLeadingMom.Pt();
      Double_t phiL = fLeadingMom.Phi();
      if(phiL < 0 ) phiL+=TMath::TwoPi();
      Double_t etaL = fLeadingMom.Eta();
      
      AliDebug(1,Form("Trigger with pt %3.2f, phi %2.2f, eta %2.2f", pt, phi, eta));
      
      if(det == kCTS)
      {
        fhChargedLeadingPt      ->Fill(pt, ptL     , GetEventWeight());
        fhChargedLeadingPhi     ->Fill(pt, phiL    , GetEventWeight());
        fhChargedLeadingEta     ->Fill(pt, etaL    , GetEventWeight());
        fhChargedLeadingDeltaPt ->Fill(pt, pt-ptL  , GetEventWeight());
        fhChargedLeadingDeltaPhi->Fill(pt, TMath::Abs(phi-phiL), GetEventWeight());
        fhChargedLeadingDeltaEta->Fill(pt, eta-etaL, GetEventWeight());
        fhChargedLeadingRatioPt ->Fill(pt, ptL/pt  , GetEventWeight());
        fhChargedLeadingXi      ->Fill(pt, TMath::Log(pt/ptL), GetEventWeight());
        if(pt > 30) fhChargedLeadingDeltaPhiRatioPt30->Fill(TMath::Abs(phi-phiL), ptL/pt, GetEventWeight());
        if(pt > 50) fhChargedLeadingDeltaPhiRatioPt50->Fill(TMath::Abs(phi-phiL), ptL/pt, GetEventWeight());
      }
      else if(det == kEMCAL)
      {
        fhNeutralLeadingPt      ->Fill(pt, ptL   , GetEventWeight());
        fhNeutralLeadingPhi     ->Fill(pt, phiL  , GetEventWeight());
        fhNeutralLeadingEta     ->Fill(pt, etaL  , GetEventWeight());
        fhNeutralLeadingDeltaPt ->Fill(pt, pt-ptL, GetEventWeight());
        fhNeutralLeadingDeltaPhi->Fill(pt, TMath::Abs(phi-phiL), GetEventWeight());
        fhNeutralLeadingDeltaEta->Fill(pt, eta-etaL, GetEventWeight());
        fhNeutralLeadingRatioPt ->Fill(pt, ptL/pt  , GetEventWeight());
        fhNeutralLeadingXi      ->Fill(pt, TMath::Log(pt/ptL), GetEventWeight());
        if(pt > 30) fhNeutralLeadingDeltaPhiRatioPt30->Fill(TMath::Abs(phi-phiL), ptL/pt, GetEventWeight());
        if(pt > 50) fhNeutralLeadingDeltaPhiRatioPt50->Fill(TMath::Abs(phi-phiL), ptL/pt, GetEventWeight());
      }
      
      // Fill Jet histograms
      
      if(!fSeveralConeAndPtCuts)
      {//just fill histograms
        if(!fReMakeJet)
        {
          fJetMom=particle->GetCorrelatedJet();
          fBkgMom=particle->GetCorrelatedBackground();
        }
        else  MakeJetFromAOD(particle);
        
        if(fJetMom.Pt() > 0)
        {//Jet was found
          FillJetHistos(particle, fJetMom,"Jet","");
          FillJetHistos(particle, fBkgMom,"Bkg","");
        }
      }
      else if(fSeveralConeAndPtCuts)
      {
        for(Int_t icone = 0; icone<fJetNCone; icone++)
        {
          fJetCone=fJetCones[icone];
          for(Int_t ipt = 0; ipt<fJetNPt;ipt++)
          {
            TString lastname ="Cone"+ fJetNameCones[icone]+"Pt"+ fJetNamePtThres[ipt];
            fJetPtThreshold=fJetPtThres[ipt];
            
            MakeJetFromAOD(particle);
            
            if(fJetMom.Pt() > 0)
            {//Jet was found
              FillJetHistos(particle, fJetMom,"Jet",lastname);
              FillJetHistos(particle, fBkgMom,"Bkg",lastname);
            }
          }//icone
        }//ipt
      }//fSeveralConeAndPtCuts
    }//Leading
  }//AOD trigger particle loop
  
  AliDebug(1,"End of jet leading cone analysis, fill histograms");
}

//_______________________________________________________________________________________________
/// Fill the jet with the particles around the leading particle with
/// R=fJetCone and pt_th = fJetPtThres. Calculate the energy of the jet and
/// fill aod with found information.
//_______________________________________________________________________________________________
void AliAnaParticleJetLeadingConeCorrelation::MakeAODJet(AliCaloTrackParticleCorrelation *particle)
{
  Double_t ptTrig   = particle->Pt();
  Double_t phiTrig  = particle->Phi();
  Double_t phil     = fLeadingMom.Phi();
  if(phil<0) phil+=TMath::TwoPi();
  Double_t etal     = fLeadingMom.Eta();
  
  // Different pt cut for jet particles in different collisions systems
  Float_t ptcut = fJetPtThreshold;
  if(fPbPb && !fSeveralConeAndPtCuts && ptTrig > fPtTriggerSelectionCut)  ptcut = fJetPtThresPbPb ;
  
  // Add charged particles to jet if they are in cone around the leading particle
  if(!GetCTSTracks())
  {
    AliFatal("Cannot construct jets without tracks, STOP analysis");
    return;
  }
  
  // Fill jet with tracks
  // Initialize reference arrays that will contain jet and background tracks
  TObjArray * reftracks     = new TObjArray;
  TObjArray * reftracksbkg  = new TObjArray;
  
  for(Int_t ipr = 0;ipr < (GetCTSTracks())->GetEntriesFast() ; ipr ++ ){
    AliVTrack* track = (AliVTrack *)((GetCTSTracks())->At(ipr)) ;
    fTrackVector.SetXYZ(track->Px(),track->Py(),track->Pz());
    
    //Particles in jet
    if(IsParticleInJetCone(fTrackVector.Eta(), fTrackVector.Phi(), etal, phil)){
      
      reftracks->Add(track);
      
      if(fTrackVector.Pt() > ptcut )
      {
        fJetConstMom.SetVect(fTrackVector);
        fJetMom+=fJetConstMom;
      }
    }
    
    // Background around (phi_gamma-pi, eta_leading)
    else if(IsParticleInJetCone(fTrackVector.Eta(),fTrackVector.Phi(),etal, phiTrig)) {
      
      reftracksbkg->Add(track);
      
      if(fTrackVector.Pt() > ptcut ){
        fJetConstMom.SetVect(fTrackVector);
        fBkgMom+=fJetConstMom;
      }
    }
  }// Track loop
  
  // Add referenced tracks to AOD
  if(reftracks->GetEntriesFast() > 0 )
  {
    reftracks->SetName(Form("%sTracks",GetAODObjArrayName().Data()));
    particle->AddObjArray(reftracks);
  }
  else  AliDebug(2,"No tracks in jet cone");
  
  if(reftracksbkg->GetEntriesFast() > 0 )
  {
    reftracksbkg->SetName(Form("%sTracksBkg",GetAODObjArrayName().Data()));
    particle->AddObjArray(reftracksbkg);
  }
  else AliDebug(1,"No background tracks in jet cone");
  
  // Add neutral particles to jet
  // Initialize reference arrays that will contain jet and background tracks
  TObjArray * refclusters     = new TObjArray;
  TObjArray * refclustersbkg  = new TObjArray;
  if(!fJetsOnlyInCTS && GetEMCALClusters())
  {
    // Get vertex for photon momentum calculation
    Double_t vertex[]  = {0,0,0} ; //vertex
    if(GetReader()->GetDataType()!= AliCaloTrackReader::kMC)
    {
      GetReader()->GetVertex(vertex);
      //if(GetReader()->GetSecondInputAODTree()) GetReader()->GetSecondInputAODVertex(vertex2);
    }
    
    for(Int_t iclus = 0;iclus < (GetEMCALClusters())->GetEntriesFast() ; iclus ++ )
    {
      AliVCluster * calo = (AliVCluster *) (GetEMCALClusters()->At(iclus)) ;
      
      // Cluster selection, not charged
      if(IsTrackMatched(calo,GetReader()->GetInputEvent())) continue ;
      
      // Get Momentum vector,
      calo->GetMomentum(fJetConstMom,vertex) ;//Assume that come from vertex in straight line
      
      // Particles in jet
      if(IsParticleInJetCone(fJetConstMom.Eta(),fJetConstMom.Phi(), etal, phil))
      {
        refclusters->Add(calo);
        
        if(fJetConstMom.Pt() > ptcut )  fJetMom+=fJetConstMom;
      }
      //Background around (phi_gamma-pi, eta_leading)
      else if(IsParticleInJetCone(fJetConstMom.Eta(),fJetConstMom.Phi(),etal, phiTrig))
      {
        refclustersbkg->Add(calo);
        
        if(fJetConstMom.Pt() > ptcut )  fBkgMom+=fJetConstMom;
      }
    }//cluster loop
  }//jets with neutral particles
  
  //Add referenced clusters to AOD
  if(refclusters->GetEntriesFast() > 0 )
  {
    refclusters->SetName(Form("%sClusters",GetAODObjArrayName().Data()));
    particle->AddObjArray(refclusters);
  }
  else  AliDebug(2,"No clusters in jet cone");
  
  if(refclustersbkg->GetEntriesFast() > 0 )
  {
    refclustersbkg->SetName(Form("%sClustersBkg",GetAODObjArrayName().Data()));
    particle->AddObjArray(refclustersbkg);
  }
  else AliDebug(1,"No background clusters in jet cone");
  
  // If there is any jet found, select after some criteria and
  // and fill AOD with corresponding TLorentzVector kinematics
  if(IsJetSelected(particle->Pt(), fJetMom.Pt()))
  {
    particle->SetCorrelatedJet(fJetMom);
    particle->SetCorrelatedBackground(fBkgMom);
    AliDebug(1,Form("Found jet: Trigger  pt %2.3f, Jet pt %2.3f, Bkg pt %2.3f",ptTrig,fJetMom.Pt(),fBkgMom.Pt()));
  }
}

//______________________________________________________________________________________________________
/// Fill the jet with the particles around the leading particle with
/// R=fJetCone and pt_th = fJetPtThres. Calculate the energy of the jet and
/// fill aod tlorentzvectors with jet and bakcground found.
//______________________________________________________________________________________________________
void AliAnaParticleJetLeadingConeCorrelation::MakeJetFromAOD(AliCaloTrackParticleCorrelation *particle)
{
  Double_t ptTrig  = particle->Pt();
  Double_t phiTrig = particle->Phi();
  Double_t etal    = fLeadingMom.Eta();
  Double_t phil    = fLeadingMom.Phi();
  if(phil < 0) phil+=TMath::TwoPi();
  
  TObjArray * refclusters    = particle->GetObjArray(Form("Clusters%s"   ,GetAODObjArrayName().Data()));
  TObjArray * reftracks      = particle->GetObjArray(Form("Tracks%s"     ,GetAODObjArrayName().Data()));
  TObjArray * refclustersbkg = particle->GetObjArray(Form("ClustersBkg%s",GetAODObjArrayName().Data()));
  TObjArray * reftracksbkg   = particle->GetObjArray(Form("TracksBkg%s"  ,GetAODObjArrayName().Data()));
  
  //Different pt cut for jet particles in different collisions systems
  Float_t ptcut = fJetPtThreshold;
  if(fPbPb && !fSeveralConeAndPtCuts && ptTrig > fPtTriggerSelectionCut)  ptcut = fJetPtThresPbPb ;
  
  //Fill jet with tracks
  //Particles in jet
  if(reftracks)
  {
    for(Int_t ipr = 0;ipr < reftracks->GetEntriesFast() ; ipr ++ )
    {
      AliVTrack* track = (AliVTrack *) reftracks->At(ipr) ;
      fTrackVector.SetXYZ(track->Px(),track->Py(),track->Pz());
      Float_t phi = fTrackVector.Phi();
      if(phi < 0) phi+=TMath::TwoPi();
      if(fTrackVector.Pt() > ptcut && IsParticleInJetCone(fTrackVector.Eta(), phi, etal, phil) )
      {
        fJetConstMom.SetVect(fTrackVector);
        fJetMom+=fJetConstMom;
      }
    }//jet Track loop
  }
    
  // Particles in background
  if(reftracksbkg){
    for(Int_t ipr = 0;ipr < reftracksbkg->GetEntriesFast() ; ipr ++ )
    {
      AliVTrack* track = (AliVTrack *) reftracksbkg->At(ipr) ;
      fTrackVector.SetXYZ(track->Px(),track->Py(),track->Pz());
      if(fTrackVector.Pt() > ptcut && IsParticleInJetCone(fTrackVector.Eta(),fTrackVector.Phi(),etal, phiTrig) )
      {
        fJetConstMom.SetVect(fTrackVector);
        fBkgMom+=fJetConstMom;
      }
    }//background Track loop
  }
  
  // Add neutral particles to jet
  if(!fJetsOnlyInCTS && refclusters)
  {
    // Get vertex for photon momentum calculation
    Double_t vertex[]  = {0,0,0} ; //vertex
    if(GetReader()->GetDataType()!= AliCaloTrackReader::kMC)
    {
      GetReader()->GetVertex(vertex);
    }
    
    // Loop on jet particles
    if(refclusters){
      for(Int_t iclus = 0;iclus < refclusters->GetEntriesFast() ; iclus ++ )
      {
        AliVCluster * calo = (AliVCluster *) refclusters->At(iclus) ;
        
        calo->GetMomentum(fJetConstMom,vertex) ;//Assume that come from vertex in straight line
        
        if(fJetConstMom.Pt() > ptcut && IsParticleInJetCone(fJetConstMom.Eta(),fJetConstMom.Phi(), etal, phil)) fJetMom+=fJetConstMom;
      }//jet cluster loop
    }
    
    // Loop on background particles
    if(refclustersbkg)
    {
      for(Int_t iclus = 0;iclus < refclustersbkg->GetEntriesFast() ; iclus ++ )
      {
        AliVCluster * calo = (AliVCluster *) refclustersbkg->At(iclus) ;
        
        calo->GetMomentum(fJetConstMom,vertex) ;//Assume that come from vertex in straight line
        
        if( fJetConstMom.Pt() > ptcut && IsParticleInJetCone(fJetConstMom.Eta(),fJetConstMom.Phi(),etal, phiTrig)) fBkgMom+=fJetConstMom;
      }// background cluster loop
    }
  }// clusters in jet
  
  // If there is any jet found, leave jet and bkg as they are,
  // if not set them to 0.
  if(!IsJetSelected(particle->Pt(), fJetMom.Pt()))
  {
    fJetMom.SetPxPyPzE(0.,0.,0.,0.);
    fBkgMom.SetPxPyPzE(0.,0.,0.,0.);
  }
  else AliDebug(1,Form("Found jet: Trigger  pt %2.3f, Jet pt %2.3f, Bkg pt %2.3f",ptTrig,fJetMom.Pt(),fBkgMom.Pt()));
}

//__________________________________________________________________
// Print some relevant parameters set for the analysis.
//__________________________________________________________________
void AliAnaParticleJetLeadingConeCorrelation::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print  %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  if(fJetsOnlyInCTS)printf("Jets reconstructed in CTS \n");
  else printf("Jets reconstructed in CTS+EMCAL \n");
  
  if(fPbPb) printf("PbPb events, pT cut in jet cone energy reconstruction %2.1f \n", fJetPtThreshold);
  else printf("pp events, pT cut in jet cone energy reconstruction %2.1f \n", fJetPtThresPbPb);
  
  printf("If pT of trigger < %2.3f, select jets as in pp? \n", fPtTriggerSelectionCut);
  
  printf("Phi gamma-Leading        <     %3.2f\n", fDeltaPhiMaxCut) ;
  printf("Phi gamma-Leading        >     %3.2f\n", fDeltaPhiMinCut) ;
  printf("pT Leading / pT Trigger  <     %3.2f\n", fLeadingRatioMaxCut) ;
  printf("pT Leading / pT Trigger  >     %3.2f\n", fLeadingRatioMinCut) ;
  
  if(fSelect == 2){
    printf("pT Jet / pT Gamma                     <    %3.2f\n", fJetRatioMaxCut) ;
    printf("pT Jet / pT Gamma                     >    %3.2f\n", fJetRatioMinCut) ;
    printf("pT Jet (Only CTS)/ pT Trigger   <    %3.2f\n", fJetCTSRatioMaxCut) ; 
    printf("pT Jet (Only CTS)/ pT Trigger   >    %3.2f\n", fJetCTSRatioMinCut) ;
  }
  else if(fSelect == 0)
    printf("Accept all reconstructed jets \n") ;
  else   if(fSelect == 1)
    printf("Accept jets depending on trigger energy \n") ;
  else 
    printf("Wrong jet selection option:   %d \n", fSelect) ;
  
  printf("Isolated Trigger?  %d\n", fSelectIsolated) ;
} 

