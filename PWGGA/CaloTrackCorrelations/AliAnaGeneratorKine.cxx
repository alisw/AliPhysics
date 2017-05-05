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
#include "TParticle.h"
#include "TDatabasePDG.h"

//---- ANALYSIS system ----
#include "AliAnaGeneratorKine.h" 
#include "AliMCEvent.h"  
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCParticle.h"

/// \cond CLASSIMP
ClassImp(AliAnaGeneratorKine) ;
/// \endcond

//__________________________________________
/// Default Constructor. Initialize parameters with default values.
//__________________________________________
AliAnaGeneratorKine::AliAnaGeneratorKine() :
AliAnaCaloTrackCorrBaseClass(), 
fTriggerDetector(),    fTriggerDetectorString(),
fFidCutTrigger(0),
fMinChargedPt(0),      fMinNeutralPt(0),
//fParton2(0),           fParton3(0),
fParton6(),            fParton7(),
fParton6PDG(0),        fParton7PDG(0),
fJet6(),               fJet7(),
fTrigger(),            fLVTmp(),
fNPrimaries(0),        fPtHard(0),
fhPtHard(0),           fhPtParton(0),         fhPtJet(0),
fhPtPartonPtHard(0),   fhPtJetPtHard(0),      fhPtJetPtParton(0),
fhPtOtherDecayMesonId(0),
fhPtPi0Status(0),      fhPtEtaStatus(0),
fhPtPi0DecayStatus(0), fhPtEtaDecayStatus(0), fhPtOtherDecayStatus(0),
fhPtPi0Not2Gamma(0),   fhPtEtaNot2Gamma(0)
{
  InitParameters();

  for(Int_t p = 0; p < fgkNmcPrimTypes; p++)
  {
    fhPt      [p] = 0;
    fhPtOrigin[p] = 0;
    fhPtOriginNotFinal[p] = 0;
    fhPhi     [p] = 0;
    fhEta     [p] = 0;
    fhEtaPhi  [p] = 0;

    fhPhiStatus[p] = 0;
    fhEtaStatus[p] = 0;
    
    for(Int_t i = 0; i < fgkNIso; i++)
    {
      fhPtLeading[p][i]          = 0;
      fhPtLeadingSumPt[p][i]     = 0;
      fhPtLeadingIsolated[p][i]  = 0;
      for(Int_t j = 0; j < 2; j++)
      {
        fhZHard[p][j][i]           = 0;
        fhZHardIsolated[p][j][i]   = 0;
        fhZParton[p][j][i]         = 0;
        fhZPartonIsolated[p][j][i] = 0;
        fhZJet[p][j][i]            = 0;
        fhZJetIsolated[p][j][i]    = 0;
        fhXE[p][j][i]              = 0;
        fhXEIsolated[p][j][i]      = 0;
        fhXEUE[p][j][i]            = 0;
        fhXEUEIsolated[p][j][i]    = 0;
        
        fhPtPartonTypeNear[p][j][i]         = 0;
        fhPtPartonTypeNearIsolated[p][j][i] = 0;
        
        fhPtPartonTypeAway[p][j][i]         = 0;
        fhPtPartonTypeAwayIsolated[p][j][i] = 0;
        
        if( p == 0 )  fhPtAcceptedGammaJet[j][i] = 0;
      }
    }
  }
}

//___________________________________________________________________________
/// Correlate trigger with partons or jets, get z.
//___________________________________________________________________________
Bool_t  AliAnaGeneratorKine::CorrelateWithPartonOrJet(Int_t   indexTrig,
                                                      Int_t   partType,
                                                      Bool_t  leading [fgkNIso],
                                                      Bool_t  isolated[fgkNIso],
                                                      Int_t & iparton )  
{
  AliDebug(1,"Start");
  
  if( fNPrimaries < 7 )
  {
    AliDebug(1,Form("End, not enough partons, n primaries %d",fNPrimaries));
    return kFALSE ;
  }
  
  //
  // Get the index of the mother
  //
  if     ( GetReader()->ReadStack() )
  {
    iparton =  (GetMC()->Particle(indexTrig))->GetFirstMother();
    TParticle * mother = GetMC()->Particle(iparton);
    while (iparton > 7)
    {
      iparton = mother->GetFirstMother();
      if(iparton < 0)
      {
        AliDebug(1,"Negative index, skip ESD event");
        return kFALSE;
      }
      mother = GetMC()->Particle(iparton);
    }
  }
  else if( GetReader()->ReadAODMCParticles() )
  {
    iparton =  ((AliAODMCParticle*) GetMC()->GetTrack(indexTrig))->GetMother();
    AliAODMCParticle * mother = (AliAODMCParticle*) GetMC()->GetTrack(iparton);
    while (iparton > 7)
    {
      iparton   = mother->GetMother();
      if(iparton < 0)
      {
        AliDebug(1,"Negative index, skip AOD event");
        return kFALSE;
      }
      mother = (AliAODMCParticle*) GetMC()->GetTrack(iparton);
    }
  }
  
  //printf("Mother is parton %d with pdg %d\n",imom,mother->GetPdgCode());
  
  if(iparton < 6)
  {
    AliDebug(1,Form("This particle is not from hard process - pdg %d, parton index %d\n",partType, iparton));
    return kFALSE; 
  }
  
  //
  // Get the kinematics
  //
  Float_t  ptTrig  = fTrigger.Pt();
  Float_t phiTrig  = fTrigger.Phi();
  Float_t etaTrig  = fTrigger.Eta();
  
  AliDebug(1,Form("Trigger pdg %d pT %2.2f, phi %2.2f, eta %2.2f", partType,ptTrig,phiTrig*TMath::RadToDeg(),etaTrig));
  
  Float_t jetPt    = fJet6.Pt();
  Float_t jetPhi   = fJet6.Phi();
  Float_t jetEta   = fJet6.Eta();
  
  AliDebug(1,Form("Jet 6 pT %2.2f, phi %2.2f, eta %2.2f",jetPt,jetPhi*TMath::RadToDeg(),jetEta));
  
  Float_t awayJetPt  = fJet7.Pt();
  Float_t awayJetEta = fJet7.Eta();
  Float_t awayJetPhi = fJet7.Phi();
  
  AliDebug(1,Form("Jet 7 pT %2.2f, phi %2.2f, eta %2.2f",awayJetPt,awayJetPhi*TMath::RadToDeg(),awayJetEta));
  
  Float_t partonPt = fParton6.Pt();

  Int_t nearPDG = fParton6PDG;
  Int_t awayPDG = fParton7PDG;
  
  AliDebug(1,Form("Parton6 pT pT %2.2f, pdg %d",fParton6.Pt(),fParton6PDG));
  AliDebug(1,Form("Parton7 pT pT %2.2f, pdg %d",fParton7.Pt(),fParton7PDG));
  
  if ( iparton == 7 )
  {
    partonPt = fParton7.Pt();
    
    jetPt  = fJet7.Pt();
    jetPhi = fJet7.Phi();
    jetEta = fJet7.Eta();
    
    awayJetPt  = fJet6.Pt();
    awayJetEta = fJet6.Eta();
    awayJetPhi = fJet6.Phi();

    nearPDG = fParton7PDG;
    awayPDG = fParton6PDG;
  }

  Float_t deltaPhi = TMath::Abs(phiTrig-awayJetPhi) *TMath::RadToDeg();
  AliDebug(1,Form("Trigger Away jet phi %2.2f\n",deltaPhi));
  
  //
  // Get id of parton in near and away side
  //
  Int_t away = -1;
  Int_t near = -1;
  if     (nearPDG == 22) near = 0;
  else if(nearPDG == 21) near = 1;
  else                   near = 2;
  
  if     (awayPDG == 22) away = 0;
  else if(awayPDG == 21) away = 1;
  else                   away = 2;
  
  // RATIOS
  
  Float_t zHard = -1;
  Float_t zPart = -1;
  Float_t zJet  = -1;
  
  if( fPtHard  > 0 ) zHard = ptTrig / fPtHard ;
  if( partonPt > 0 ) zPart = ptTrig / partonPt;
  if( jetPt    > 0 ) zJet  = ptTrig / jetPt   ;
  
  //if(zHard > 1 ) printf("*** Particle energy larger than pT hard z=%f\n",zHard);
  
  //printf("Z: hard %2.2f, parton %2.2f, jet %2.2f\n",zHard,zPart,zJet);

  // conditions loop
  for( Int_t i = 0; i < fgkNIso ; i++ )
  {
    fhPtPartonTypeNear[partType][leading[i]][i]->Fill(ptTrig, near, GetEventWeight());
    fhPtPartonTypeAway[partType][leading[i]][i]->Fill(ptTrig, away, GetEventWeight());
    
    fhZHard  [partType][leading[i]][i]->Fill(ptTrig, zHard, GetEventWeight());
    fhZParton[partType][leading[i]][i]->Fill(ptTrig, zPart, GetEventWeight());
    fhZJet   [partType][leading[i]][i]->Fill(ptTrig, zJet , GetEventWeight());
    
    if(isolated[i])
    {
      fhPtPartonTypeNearIsolated[partType][leading[i]][i]->Fill(ptTrig, near, GetEventWeight());
      fhPtPartonTypeAwayIsolated[partType][leading[i]][i]->Fill(ptTrig, away, GetEventWeight());
      
      fhZHardIsolated  [partType][leading[i]][i]->Fill(ptTrig, zHard, GetEventWeight());
      fhZPartonIsolated[partType][leading[i]][i]->Fill(ptTrig, zPart, GetEventWeight());
      fhZJetIsolated   [partType][leading[i]][i]->Fill(ptTrig, zJet , GetEventWeight());
    }
    
    if(partType == kmcPrimPhoton && deltaPhi < 220 && deltaPhi > 140 && TMath::Abs(awayJetEta) < 0.6)
    {
      //printf("Accept jet\n");
      fhPtAcceptedGammaJet[leading[i]][i]->Fill(ptTrig, away, GetEventWeight());
    }
    //else printf("Reject jet!!!\n");
    
  } // conditions loop
  
  AliDebug(1,"End TRUE");
  
  return kTRUE;
}

//____________________________________________________
/// Create histograms to be saved in output file.
//____________________________________________________
TList *  AliAnaGeneratorKine::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("GenKineHistos") ; 
  
  Int_t   nptbins    = GetHistogramRanges()->GetHistoPtBins();
  Float_t ptmax      = GetHistogramRanges()->GetHistoPtMax();
  Float_t ptmin      = GetHistogramRanges()->GetHistoPtMin();

  Int_t   nptsumbins = GetHistogramRanges()->GetHistoNPtSumBins();
  Float_t ptsummax   = GetHistogramRanges()->GetHistoPtSumMax();
  Float_t ptsummin   = GetHistogramRanges()->GetHistoPtSumMin();
  
  fhPtHard  = new TH1F("hPtHard"," pt hard for selected triggers",nptbins,ptmin,ptmax); 
  fhPtHard->SetXTitle("#it{p}_{T}^{hard} (GeV/#it{c})");
  outputContainer->Add(fhPtHard);
  
  fhPtParton  = new TH1F("hPtParton"," pt parton for selected triggers",nptbins,ptmin,ptmax); 
  fhPtParton->SetXTitle("#it{p}_{T}^{parton} (GeV/#it{c})");
  outputContainer->Add(fhPtParton);
  
  fhPtJet  = new TH1F("hPtJet"," pt jet for selected triggers",nptbins,ptmin,ptmax); 
  fhPtJet->SetXTitle("#it{p}_{T}^{jet} (GeV/#it{c})");
  outputContainer->Add(fhPtJet);
  
  fhPtPartonPtHard  = new TH2F("hPtPartonPtHard","parton pt / pt hard for selected triggers",nptbins,ptmin,ptmax,200,0,2); 
  fhPtPartonPtHard->SetXTitle("#it{p}_{T}^{hard} (GeV/#it{c})");
  fhPtPartonPtHard->SetYTitle("#it{p}_{T}^{parton}/#it{p}_{T}^{hard}");
  outputContainer->Add(fhPtPartonPtHard);
  
  fhPtJetPtHard  = new TH2F("hPtJetPtHard","jet pt / pt hard for selected triggers",nptbins,ptmin,ptmax,200,0,2); 
  fhPtJetPtHard->SetXTitle("#it{p}_{T}^{hard} (GeV/#it{c})");
  fhPtJetPtHard->SetYTitle("#it{p}_{T}^{jet}/#it{p}_{T}^{hard}");
  outputContainer->Add(fhPtJetPtHard);
  
  fhPtJetPtParton  = new TH2F("hPtJetPtParton","parton pt / pt hard for selected triggers",nptbins,ptmin,ptmax,200,0,2); 
  fhPtJetPtParton->SetXTitle("#it{p}_{T}^{hard} (GeV/#it{c})");
  fhPtJetPtParton->SetYTitle("#it{p}_{T}^{jet}/#it{p}_{T}^{parton}");
  outputContainer->Add(fhPtJetPtParton);

  fhPtPi0Not2Gamma  = new TH1F("hPtPi0Not2Gamma","#pi^{0} decay other than 2 #gamma",nptbins,ptmin,ptmax);
  fhPtPi0Not2Gamma->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtPi0Not2Gamma);

  fhPtEtaNot2Gamma  = new TH1F("hPtEtaNot2Gamma","#eta decay other than 2 #gamma",nptbins,ptmin,ptmax);
  fhPtEtaNot2Gamma->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtEtaNot2Gamma);

  fhPtGammaFromPi0Not2Gamma  = new TH1F("hPtGammaFromPi0Not2Gamma","#gamma from #pi^{0} decay other than 2 #gamma",nptbins,ptmin,ptmax);
  fhPtGammaFromPi0Not2Gamma->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtGammaFromPi0Not2Gamma);
  
  fhPtGammaFromEtaNot2Gamma  = new TH1F("hPtGammaFromEtaNot2Gamma","#gamma from #eta decay other than 2 #gamma",nptbins,ptmin,ptmax);
  fhPtGammaFromEtaNot2Gamma->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtGammaFromEtaNot2Gamma);
  
  fhPtPi0Status  = new TH2F("hPtPi0Status","#pi^{0} status",nptbins,ptmin,ptmax,101,-50,50);
  fhPtPi0Status->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fhPtPi0Status->SetYTitle("status");
  outputContainer->Add(fhPtPi0Status);

  fhPtEtaStatus  = new TH2F("hPtEtaStatus","#eta status",nptbins,ptmin,ptmax,101,-50,50);
  fhPtEtaStatus->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fhPtEtaStatus->SetYTitle("status");
  outputContainer->Add(fhPtEtaStatus);

  fhPtPi0DecayStatus  = new TH2F("hPtPi0DecayStatus","#gamma from #pi^{0}, mother status",nptbins,ptmin,ptmax,101,-50,50);
  fhPtPi0DecayStatus->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fhPtPi0DecayStatus->SetYTitle("status");
  outputContainer->Add(fhPtPi0DecayStatus);
  
  fhPtEtaDecayStatus  = new TH2F("hPtEtaStatus","#gamma from #eta, mother status",nptbins,ptmin,ptmax,101,-50,50);
  fhPtEtaDecayStatus->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fhPtEtaDecayStatus->SetYTitle("status");
  outputContainer->Add(fhPtEtaDecayStatus);

  fhPtOtherDecayStatus  = new TH2F("hPtOtherDecayStatus","#gamma from other decay particle, mother status",nptbins,ptmin,ptmax,101,-50,50);
  fhPtOtherDecayStatus->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fhPtOtherDecayStatus->SetYTitle("status");
  outputContainer->Add(fhPtOtherDecayStatus);
  
  fhPtOtherDecayMesonId     = new TH2F("hPtOtherDecayMesonId","Particle decaying into #gamma, Id",nptbins,ptmin,ptmax,8,0,8) ;
  fhPtOtherDecayMesonId->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fhPtOtherDecayMesonId->SetYTitle("Origin");
  fhPtOtherDecayMesonId->GetYaxis()->SetBinLabel(1 ,"Resonances");
  fhPtOtherDecayMesonId->GetYaxis()->SetBinLabel(2 ,"#rho");
  fhPtOtherDecayMesonId->GetYaxis()->SetBinLabel(3 ,"#omega");
  fhPtOtherDecayMesonId->GetYaxis()->SetBinLabel(4 ,"K");
  fhPtOtherDecayMesonId->GetYaxis()->SetBinLabel(5 ,"Other");
  fhPtOtherDecayMesonId->GetYaxis()->SetBinLabel(6 ,"#eta");
  fhPtOtherDecayMesonId->GetYaxis()->SetBinLabel(7 ,"#eta prime");
  outputContainer->Add(fhPtOtherDecayMesonId) ;
  
  TString name   [] = {"","_EMC","_Photon","_EMC_Photon"};
  TString title  [] = {"",", neutral in EMCal",", neutral only #gamma-like",", neutral in EMCal and only #gamma-like"};
  TString leading[] = {"NotLeading","Leading"};
  
  TString partTitl[] = {"#gamma_{direct}","#gamma_{decay}^{#pi}","#gamma_{decay}^{#eta}","#gamma_{decay}^{other}","#pi^{0}","#eta"};
  TString particle[] = {"DirectPhoton"   ,"Pi0DecayPhoton"      ,"EtaDecayPhoton"       ,"OtherDecayPhoton"      ,"Pi0"    ,"Eta"};

  for(Int_t p = 0; p < fgkNmcPrimTypes; p++)
  {
    fhPt[p]  = new TH1F(Form("h%sPt",particle[p].Data()),Form("Input %s p_{T}",partTitl[p].Data()),nptbins,ptmin,ptmax);
    fhPt[p]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPt[p]);

    fhPhi[p]  = new TH1F(Form("h%sPhi",particle[p].Data()),Form("Input %s #phi with #it{p}_{T} > %f GeV/c",partTitl[p].Data(),GetMinPt()),180,0,TMath::TwoPi());
    fhPhi[p]->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhi[p]);

    fhEta[p]  = new TH1F(Form("h%sEta",particle[p].Data()),Form("Input %s #eta with #it{p}_{T} > %f GeV/c",partTitl[p].Data(),GetMinPt()),200,-2,2);
    fhEta[p]->SetXTitle("#eta");
    outputContainer->Add(fhEta[p]);

    fhEtaPhi[p]  = new TH2F(Form("h%sEtaPhi",particle[p].Data()),
                            Form("Input %s #eta vs #phi with #it{p}_{T} > %f GeV/c",partTitl[p].Data(),GetMinPt()),
                            200,-2,2,180,0,TMath::TwoPi());
    fhEtaPhi[p]->SetXTitle("#eta");
    fhEtaPhi[p]->SetYTitle("#phi (rad)");
    outputContainer->Add(fhEtaPhi[p]);

    fhPhiStatus[p]  = new TH2F(Form("h%sPhiStatus",particle[p].Data()),
                               Form("Input %s #phi vs status code with #it{p}_{T} > %f GeV/c",partTitl[p].Data(),GetMinPt()),
                               180,0,TMath::TwoPi(),101,-50,50);
    fhPhiStatus[p]->SetXTitle("#phi (rad)");
    fhPhiStatus[p]->SetYTitle("status code");
    outputContainer->Add(fhPhiStatus[p]);
    
    fhEtaStatus[p]  = new TH2F(Form("h%sEtaStatus",particle[p].Data()),
                               Form("Input %s #eta vs status code with #it{p}_{T} > %f GeV/c",partTitl[p].Data(),GetMinPt()),
                               200,-2,2,101,-50,50);
    fhEtaStatus[p]->SetXTitle("#eta");
    fhEtaStatus[p]->SetYTitle("status code");
    outputContainer->Add(fhEtaStatus[p]);

    
    fhPtOrigin[p]     = new TH2F(Form("h%sPtOrigin",particle[p].Data()),Form("Input %s p_{T} vs origin",partTitl[p].Data()),nptbins,ptmin,ptmax,11,0,11) ;
    fhPtOrigin[p]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhPtOrigin[p]->SetYTitle("Origin");
    fhPtOrigin[p]->GetYaxis()->SetBinLabel(1 ,"Status 21");
    fhPtOrigin[p]->GetYaxis()->SetBinLabel(2 ,"Quark");
    fhPtOrigin[p]->GetYaxis()->SetBinLabel(3 ,"qq Resonances ");
    fhPtOrigin[p]->GetYaxis()->SetBinLabel(4 ,"Resonances");
    fhPtOrigin[p]->GetYaxis()->SetBinLabel(5 ,"#rho");
    fhPtOrigin[p]->GetYaxis()->SetBinLabel(6 ,"#omega");
    fhPtOrigin[p]->GetYaxis()->SetBinLabel(7 ,"K");
    fhPtOrigin[p]->GetYaxis()->SetBinLabel(8 ,"Other");
    fhPtOrigin[p]->GetYaxis()->SetBinLabel(9 ,"#eta");
    fhPtOrigin[p]->GetYaxis()->SetBinLabel(10 ,"#eta prime");
    outputContainer->Add(fhPtOrigin[p]) ;

    fhPtOriginNotFinal[p]     = new TH2F(Form("h%sPtOriginNotFinal",particle[p].Data()),Form("Input %s p_{T} vs origin, status 0",partTitl[p].Data()),nptbins,ptmin,ptmax,11,0,11) ;
    fhPtOriginNotFinal[p]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhPtOriginNotFinal[p]->SetYTitle("Origin");
    fhPtOriginNotFinal[p]->GetYaxis()->SetBinLabel(1 ,"Status 21");
    fhPtOriginNotFinal[p]->GetYaxis()->SetBinLabel(2 ,"Quark");
    fhPtOriginNotFinal[p]->GetYaxis()->SetBinLabel(3 ,"qq Resonances ");
    fhPtOriginNotFinal[p]->GetYaxis()->SetBinLabel(4 ,"Resonances");
    fhPtOriginNotFinal[p]->GetYaxis()->SetBinLabel(5 ,"#rho");
    fhPtOriginNotFinal[p]->GetYaxis()->SetBinLabel(6 ,"#omega");
    fhPtOriginNotFinal[p]->GetYaxis()->SetBinLabel(7 ,"K");
    fhPtOriginNotFinal[p]->GetYaxis()->SetBinLabel(8 ,"Other");
    fhPtOriginNotFinal[p]->GetYaxis()->SetBinLabel(9 ,"#eta");
    fhPtOriginNotFinal[p]->GetYaxis()->SetBinLabel(10 ,"#eta prime");
    outputContainer->Add(fhPtOriginNotFinal[p]) ;
    
    for(Int_t i = 0; i < fgkNIso; i++)
    {
      // Pt
      
      fhPtLeading[p][i]  = new TH1F(Form("h%sPtLeading%s", particle[p].Data(), name[i].Data()),
                                    Form("%s: Leading of all particles%s",partTitl[p].Data(),title[i].Data()),
                                    nptbins,ptmin,ptmax);
      fhPtLeading[p][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtLeading[p][i]);
      
      fhPtLeadingIsolated[p][i]  = new TH1F(Form("h%sPtLeadingIsolated%s", particle[p].Data(), name[i].Data()),
                                            Form("%s: Leading of all particles%s, isolated",partTitl[p].Data(),title[i].Data()),
                                            nptbins,ptmin,ptmax);
      fhPtLeadingIsolated[p][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtLeadingIsolated[p][i]);
      
      fhPtLeadingSumPt[p][i]  = new TH2F(Form("h%sPtLeadingSumPt%s", particle[p].Data(), name[i].Data()),
                                         Form("%s: Leading of all particles%s",partTitl[p].Data(),title[i].Data()),
                                         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhPtLeadingSumPt[p][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtLeadingSumPt[p][i]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtLeadingSumPt[p][i]);
      
      
      // Leading or not loop
      for(Int_t j = 0; j < fgkNLead; j++)
      {
        if(p==0)
        {
          fhPtAcceptedGammaJet[j][i]  = new TH2F(Form("hPtAcceptedGammaJet%s%s",           leading[j].Data(), name[i].Data()),
                                                 Form("#gamma-jet: %s of all particles%s", leading[j].Data(), title[i].Data()),
                                                 nptbins,ptmin,ptmax,3,0,3);
          fhPtAcceptedGammaJet[j][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhPtAcceptedGammaJet[j][i]->SetYTitle("Parton type");
          fhPtAcceptedGammaJet[j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
          fhPtAcceptedGammaJet[j][i]->GetYaxis()->SetBinLabel(2,"g");
          fhPtAcceptedGammaJet[j][i]->GetYaxis()->SetBinLabel(3,"q");
          outputContainer->Add(fhPtAcceptedGammaJet[j][i]);
        }
        // Near side parton
        
        fhPtPartonTypeNear[p][j][i]  = new TH2F(Form("h%sPtPartonTypeNear%s%s",   particle[p].Data(), leading[j].Data(), name[i].Data()),
                                                Form("%s: %s of all particles%s", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                                nptbins,ptmin,ptmax,3,0,3);
        fhPtPartonTypeNear[p][j][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPtPartonTypeNear[p][j][i]->SetYTitle("Parton type");
        fhPtPartonTypeNear[p][j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
        fhPtPartonTypeNear[p][j][i]->GetYaxis()->SetBinLabel(2,"g");
        fhPtPartonTypeNear[p][j][i]->GetYaxis()->SetBinLabel(3,"q");
        outputContainer->Add(fhPtPartonTypeNear[p][j][i]);
        
        fhPtPartonTypeNearIsolated[p][j][i]  = new TH2F(Form("h%sPtPartonTypeNear%sIsolated%s",     particle[p].Data(), leading[j].Data(), name[i].Data()),
                                                        Form("%s: %s of all particles%s, isolated", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                                        nptbins,ptmin,ptmax,3,0,3);
        fhPtPartonTypeNearIsolated[p][j][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPtPartonTypeNearIsolated[p][j][i]->SetYTitle("Parton type");
        fhPtPartonTypeNearIsolated[p][j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
        fhPtPartonTypeNearIsolated[p][j][i]->GetYaxis()->SetBinLabel(2,"g");
        fhPtPartonTypeNearIsolated[p][j][i]->GetYaxis()->SetBinLabel(3,"q");
        outputContainer->Add(fhPtPartonTypeNearIsolated[p][j][i]);
        
        
        // Away side parton
        
        fhPtPartonTypeAway[p][j][i]  = new TH2F(Form("h%sPtPartonTypeAway%s%s",   particle[p].Data(), leading[j].Data(), name[i].Data()),
                                                Form("%s: %s of all particles%s", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                                nptbins,ptmin,ptmax,3,0,3);
        fhPtPartonTypeAway[p][j][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPtPartonTypeAway[p][j][i]->SetYTitle("Parton type");
        fhPtPartonTypeAway[p][j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
        fhPtPartonTypeAway[p][j][i]->GetYaxis()->SetBinLabel(2,"g");
        fhPtPartonTypeAway[p][j][i]->GetYaxis()->SetBinLabel(3,"q");
        outputContainer->Add(fhPtPartonTypeAway[p][j][i]);
        
        fhPtPartonTypeAwayIsolated[p][j][i]  = new TH2F(Form("h%sPtPartonTypeAway%sIsolated%s",     particle[p].Data(), leading[j].Data(), name[i].Data()),
                                                        Form("%s: %s of all particles%s, isolated", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                                        nptbins,ptmin,ptmax,3,0,3);
        fhPtPartonTypeAwayIsolated[p][j][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPtPartonTypeAwayIsolated[p][j][i]->SetYTitle("Parton type");
        fhPtPartonTypeAwayIsolated[p][j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
        fhPtPartonTypeAwayIsolated[p][j][i]->GetYaxis()->SetBinLabel(2,"g");
        fhPtPartonTypeAwayIsolated[p][j][i]->GetYaxis()->SetBinLabel(3,"q");
        outputContainer->Add(fhPtPartonTypeAwayIsolated[p][j][i]);
        
        // zHard
        
        fhZHard[p][j][i]  = new TH2F(Form("h%sZHard%s%s",                               particle[p].Data(), leading[j].Data(), name[i].Data()),
                                     Form("#it{z}_{Hard} of %s: %s of all particles%s", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                     nptbins,ptmin,ptmax,200,0,2);
        fhZHard[p][j][i]->SetYTitle("#it{p}_{T}^{particle}/#it{p}_{T}^{hard}");
        fhZHard[p][j][i]->SetXTitle("#it{p}_{T}^{particle} (GeV/#it{c})");
        outputContainer->Add(fhZHard[p][j][i]);
        
        fhZHardIsolated[p][j][i]  = new TH2F(Form("h%sZHard%sIsolated%s",                                 particle[p].Data(), leading[j].Data(), name[i].Data()),
                                             Form("#it{z}_{Hard} of %s: %s of all particles%s, isolated", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                             nptbins,ptmin,ptmax,200,0,2);
        fhZHardIsolated[p][j][i]->SetYTitle("#it{p}_{T}^{particle}/#it{p}_{T}^{hard}");
        fhZHardIsolated[p][j][i]->SetXTitle("#it{p}_{T}^{particle} (GeV/#it{c})");
        outputContainer->Add(fhZHardIsolated[p][j][i]);
        
        // zHard
        
        fhZParton[p][j][i]  = new TH2F(Form("h%sZParton%s%s",                               particle[p].Data(), leading[j].Data(), name[i].Data()),
                                       Form("#it{z}_{Parton} of %s: %s of all particles%s", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                       nptbins,ptmin,ptmax,200,0,2);
        fhZParton[p][j][i]->SetYTitle("#it{p}_{T}^{particle}/#it{p}_{T}^{hard}");
        fhZParton[p][j][i]->SetXTitle("#it{p}_{T}^{particle} (GeV/#it{c})");
        outputContainer->Add(fhZParton[p][j][i]);
        
        fhZPartonIsolated[p][j][i]  = new TH2F(Form("h%sZParton%sIsolated%s",                                 particle[p].Data(), leading[j].Data(), name[i].Data()),
                                               Form("#it{z}_{Parton} of %s: %s of all particles%s, isolated", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                               nptbins,ptmin,ptmax,200,0,2);
        fhZPartonIsolated[p][j][i]->SetYTitle("#it{p}_{T}^{particle}/#it{p}_{T}^{hard}");
        fhZPartonIsolated[p][j][i]->SetXTitle("#it{p}_{T}^{particle} (GeV/#it{c})");
        outputContainer->Add(fhZPartonIsolated[p][j][i]);
        
        
        // zJet
        
        fhZJet[p][j][i]  = new TH2F(Form("h%sZJet%s%s",                               particle[p].Data(), leading[j].Data(), name[i].Data()),
                                    Form("#it{z}_{Jet} of %s: %s of all particles%s", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                    nptbins,ptmin,ptmax,200,0,2);
        fhZJet[p][j][i]->SetYTitle("#it{p}_{T}^{particle}/#it{p}_{T}^{hard}");
        fhZJet[p][j][i]->SetXTitle("#it{p}_{T}^{particle} (GeV/#it{c})");
        outputContainer->Add(fhZJet[p][j][i]);
        
        fhZJetIsolated[p][j][i]  = new TH2F(Form("h%sZJet%sIsolated%s",                                 particle[p].Data(), leading[j].Data(), name[i].Data()),
                                            Form("#it{z}_{Jet} of %s: %s of all particles%s, isolated", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                            nptbins,ptmin,ptmax,200,0,2);
        fhZJetIsolated[p][j][i]->SetYTitle("#it{p}_{T}^{particle}/#it{p}_{T}^{hard}");
        fhZJetIsolated[p][j][i]->SetXTitle("#it{p}_{T}^{particle} (GeV/#it{c})");
        outputContainer->Add(fhZJetIsolated[p][j][i]);
        
        
        // XE
        
        fhXE[p][j][i]  = new TH2F(Form("h%sXE%s%s",                                 particle[p].Data(), leading[j].Data(), name[i].Data()),
                                  Form("#it{z}_{Jet} of %s: %s of all particles%s", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                  nptbins,ptmin,ptmax,200,0,2);
        fhXE[p][j][i]->SetYTitle("#it{p}_{T}^{particle}/#it{p}_{T}^{hard}");
        fhXE[p][j][i]->SetXTitle("#it{p}_{T}^{particle} (GeV/#it{c})");
        outputContainer->Add(fhXE[p][j][i]);
        
        fhXEIsolated[p][j][i]  = new TH2F(Form("h%sXE%sIsolated%s",                                   particle[p].Data(), leading[j].Data(), name[i].Data()),
                                          Form("#it{z}_{Jet} of %s: %s of all particles%s, isolated", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                          nptbins,ptmin,ptmax,200,0,2);
        fhXEIsolated[p][j][i]->SetYTitle("#it{p}_{T}^{particle}/#it{p}_{T}^{hard}");
        fhXEIsolated[p][j][i]->SetXTitle("#it{p}_{T}^{particle} (GeV/#it{c})");
        outputContainer->Add(fhXEIsolated[p][j][i]);
        
        
        // XE from UE
        
        fhXEUE[p][j][i]  = new TH2F(Form("h%sXEUE%s%s",                               particle[p].Data(), leading[j].Data(), name[i].Data()),
                                    Form("#it{z}_{Jet} of %s: %s of all particles%s", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                    nptbins,ptmin,ptmax,200,0,2);
        fhXEUE[p][j][i]->SetYTitle("#it{p}_{T}^{particle}/#it{p}_{T}^{hard}");
        fhXEUE[p][j][i]->SetXTitle("#it{p}_{T}^{particle} (GeV/#it{c})");
        outputContainer->Add(fhXEUE[p][j][i]);
        
        fhXEUEIsolated[p][j][i]  = new TH2F(Form("h%sXEUE%sIsolated%s",                                 particle[p].Data(), leading[j].Data(), name[i].Data()),
                                            Form("#it{z}_{Jet} of %s: %s of all particles%s, isolated", partTitl[p].Data(), leading[j].Data(), title[i].Data()),
                                            nptbins,ptmin,ptmax,200,0,2); 
        fhXEUEIsolated[p][j][i]->SetYTitle("#it{p}_{T}^{particle}/#it{p}_{T}^{hard}");
        fhXEUEIsolated[p][j][i]->SetXTitle("#it{p}_{T}^{particle} (GeV/#it{c})");
        outputContainer->Add(fhXEUEIsolated[p][j][i]);
      }
    }
  }
  
  return outputContainer;
}

//____________________________________________
/// Fill data members with partons,jets and generated pt hard.
//____________________________________________
void  AliAnaGeneratorKine::GetPartonsAndJets()
{
  AliDebug(1,"Start");
  
//  if( nPrimary > 2 ) fParton2 = GetMC()->Particle(2) ;
//  if( nPrimary > 3 ) fParton3 = GetMC()->Particle(3) ;
  
  Float_t p6phi = -1 ;
  Float_t p6eta = -10;
  Float_t p6pt  =  0 ;
    
  if( fNPrimaries > 6 )
  {
    p6pt  = fParton6.Pt();
    p6eta = fParton6.Eta();
    p6phi = fParton6.Phi();
    if(p6phi < 0) p6phi +=TMath::TwoPi();
  }

  Float_t p7phi = -1 ;
  Float_t p7eta = -10;
  Float_t p7pt  =  0 ;
    
  if( fNPrimaries > 7 )
  {
    p7pt  = fParton7.Pt();
    p7phi = fParton7.Eta();
    p7phi = fParton7.Phi();
    if(p7phi < 0) p7phi +=TMath::TwoPi();
  }
  
  //printf("parton6: pt %2.2f, eta %2.2f, phi %2.2f with pdg %d\n",p6pt,p6eta,p6phi, fParton6PDG);
  //printf("parton7: pt %2.2f, eta %2.2f, phi %2.2f with pdg %d\n",p7pt,p7eta,p7phi, fParton7PDG);
  
  // Get the jets, only for pythia
  if(!strcmp(GetReader()->GetGenEventHeader()->ClassName(), "AliGenPythiaEventHeader"))
  {
    AliGenPythiaEventHeader* pygeh= (AliGenPythiaEventHeader*) GetReader()->GetGenEventHeader();
    
    fPtHard = pygeh->GetPtHard();
    
    //printf("pt Hard %2.2f\n",fPtHard);
    
    const Int_t nTriggerJets =  pygeh->NTriggerJets();
        
    Float_t tmpjet[]={0,0,0,0};
    
    // select the closest jet to parton
    Float_t jet7R = 100;
    Float_t jet6R = 100;
    
    for(Int_t ijet = 0; ijet< nTriggerJets; ijet++)
    {
      pygeh->TriggerJet(ijet, tmpjet);
      
      fLVTmp.SetPxPyPzE(tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3]);
      Float_t jphi = fLVTmp.Phi();
      if(jphi < 0) jphi +=TMath::TwoPi();
      
      Double_t radius6 = GetIsolationCut()->Radius(p6eta, p6phi, fLVTmp.Eta() , jphi) ;
      Double_t radius7 = GetIsolationCut()->Radius(p7eta, p7phi, fLVTmp.Eta() , jphi) ;
      
      //printf("jet %d: pt %2.2f, eta %2.2f, phi %2.2f, r6 %2.2f, r7 %2.2f\n",ijet,jet.Pt(),jet.Eta(),jphi,radius6, radius7);
      
      if (radius6 < jet6R)
      {
        jet6R = radius6;
        fJet6 = fLVTmp;
        
      }
      
      if (radius7 < jet7R) 
      {
        jet7R = radius7;
        fJet7 = fLVTmp;
      }
    } // jet loop
    
    //printf("jet6: pt %2.2f, eta %2.2f, phi %2.2f\n",fJet6.Pt(),fJet6.Eta(),fJet6.Phi());
    //printf("jet7: pt %2.2f, eta %2.2f, phi %2.2f\n",fJet7.Pt(),fJet7.Eta(),fJet6.Phi());
  } // pythia header
  
  fhPtHard   ->Fill(fPtHard   , GetEventWeight());
  fhPtJet    ->Fill(fJet6.Pt(), GetEventWeight());
  fhPtJet    ->Fill(fJet7.Pt(), GetEventWeight());
  fhPtParton ->Fill(p6pt      , GetEventWeight());
  fhPtParton ->Fill(p7pt      , GetEventWeight());

  if( fPtHard > 0 )
  {
    fhPtPartonPtHard->Fill(fPtHard, p6pt/fPtHard      , GetEventWeight());
    fhPtPartonPtHard->Fill(fPtHard, p7pt/fPtHard      , GetEventWeight());
    fhPtJetPtHard   ->Fill(fPtHard, fJet6.Pt()/fPtHard, GetEventWeight());
    fhPtJetPtHard   ->Fill(fPtHard, fJet7.Pt()/fPtHard, GetEventWeight());
  }
  
  if( p6pt > 0 ) fhPtJetPtParton ->Fill(fPtHard, fJet6.Pt()/p6pt, GetEventWeight());
  if( p7pt > 0 ) fhPtJetPtParton ->Fill(fPtHard, fJet7.Pt()/p7pt, GetEventWeight());
  
  AliDebug(1,"End");
}

//_____________________________________________________
/// Calculate the real XE and the UE XE.
//_____________________________________________________
void AliAnaGeneratorKine::GetXE(Int_t   indexTrig,
                                Int_t   partType,
                                Bool_t  leading [fgkNIso],
                                Bool_t  isolated[fgkNIso],
                                Int_t   iparton)
{
  AliDebug(1,"Start");
  
  Float_t ptTrig  = fTrigger.Pt();
  Float_t phiTrig = fTrigger.Phi();
  if(phiTrig < 0 ) phiTrig += TMath::TwoPi();
  
  Int_t  pdg    = 0;
  Int_t  status = 0;
  Int_t  ipartonAway = 0;
  Int_t  charge = 0;
  //Loop on primaries, start from position 8, no partons
  for(Int_t ipr = 8; ipr < fNPrimaries; ipr ++ )
  {
    if(ipr==indexTrig) continue;

    // Get ESD particle kinematics
    if     ( GetReader()->ReadStack() )
    {
      TParticle * particle = GetMC()->Particle(ipr) ;
      
      pdg    = particle->GetPdgCode();
      status = particle->GetStatusCode();
      
      // Compare trigger with final state particles
      if( status != 1) continue ; // do it here to avoid crashes
      
      particle->Momentum(fLVTmp);
      
      charge = (Int_t) TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
    }
    else if( GetReader()->ReadAODMCParticles() )
    {
      AliAODMCParticle * particle = (AliAODMCParticle*) GetMC()->GetTrack(ipr) ;
      
      pdg    = particle->GetPdgCode();
      status = particle->GetStatus();
      
      // Compare trigger with final state particles
      if( status != 1) continue ; // do it here to avoid crashes
      
      fLVTmp.SetPxPyPzE(particle->Px(),particle->Py(),particle->Pz(),particle->E());
      
      charge = particle->Charge();
    }
    
    // construct xe only with charged
    if( charge == 0 ) continue;
    
    Float_t pt  = fLVTmp.Pt();
    Float_t phi = fLVTmp.Phi();
    if(phi < 0 ) phi += TMath::TwoPi();
    
    if( pt < fMinChargedPt)    continue ;
    
    Bool_t inTPC = GetFiducialCut()->IsInFiducialCut(fLVTmp.Eta(),fLVTmp.Phi(),kCTS) ;
    
    if(!inTPC) continue;
    
    // ---------------------------------------------------
    // Get the index of the mother, get from what parton
    // ESD
    if     ( GetReader()->ReadStack() )
    {
      ipartonAway =  GetMC()->Particle(ipr)->GetFirstMother();
      if(ipartonAway < 0)
      {
        AliDebug(1,"End, no mother index");
        return;
      }
      
      TParticle * mother = GetMC()->Particle(ipartonAway);
      while (ipartonAway > 7)
      {
        ipartonAway   = mother->GetFirstMother();
        if(ipartonAway < 0) break;
        mother = GetMC()->Particle(ipartonAway);
      }
    }
    else if( GetReader()->ReadAODMCParticles() )
    {
      ipartonAway =  ((AliAODMCParticle*) GetMC()->GetTrack(ipr))->GetMother();
      if(ipartonAway < 0)
      {
        AliDebug(1,"End, no mother index");
        return;
      }
      
      AliAODMCParticle * mother = (AliAODMCParticle*) GetMC()->GetTrack(ipartonAway);
      while (ipartonAway > 7)
      {
        ipartonAway   = mother->GetMother();
        if(ipartonAway < 0) break;
        mother = (AliAODMCParticle*) GetMC()->GetTrack(ipartonAway);
      }
    }
    
    //-----------------------------------------
    // Get XE of particles belonging to the jet
    // on the opposite side of the trigger
    
    Float_t xe = -pt/ptTrig*TMath::Cos(phi-phiTrig);

    if((ipartonAway==6 || ipartonAway==7) && iparton!=ipartonAway)
    {
      for( Int_t i = 0; i < fgkNIso; i++ )
      {
        fhXE[partType][leading[i]][i]  ->Fill(ptTrig, xe, GetEventWeight());
        
        if(isolated[i])
        {
          fhXEIsolated[partType][leading[i]][i]  ->Fill(ptTrig, xe, GetEventWeight());
        }
      } // conditions loop
    } // Away side
    
    //----------------------------------------------------------
    // Get the XE from particles not attached to any of the jets
    if(ipartonAway!=6 && ipartonAway!=7)
    {
      for( Int_t i = 0; i < fgkNIso; i++ )
      {
        fhXEUE[partType][leading[i]][i]  ->Fill(ptTrig, xe, GetEventWeight());
        
        if(isolated[i])
        {
          fhXEUEIsolated[partType][leading[i]][i]  ->Fill(ptTrig, xe, GetEventWeight());
        }
      } // conditions loop
    } // Away side
    
  } // primary loop
  
  AliDebug(1,"End");
}

//________________________________________
/// Initialize the parameters of the analysis.
//________________________________________
void AliAnaGeneratorKine::InitParameters()
{
  AddToHistogramsName("AnaGenKine_");
  
  fTriggerDetector = kEMCAL;
  
  fMinChargedPt    = 0.2;
  fMinNeutralPt    = 0.3;
}

//_____________________________________________________________________
/// Check if the trigger is the leading particle and if it is isolated.
/// In case of neutral particles check all neutral or neutral in EMCAL acceptance.
//_____________________________________________________________________
void  AliAnaGeneratorKine::IsLeadingAndIsolated(Int_t indexTrig,
                                                Int_t partType,
                                                Bool_t leading[fgkNIso],
                                                Bool_t isolated[fgkNIso]) 
{
  AliDebug(1,"Start");

  Float_t ptMaxCharged       = 0; // all charged
  Float_t ptMaxNeutral       = 0; // all neutral
  Float_t ptMaxNeutEMCAL     = 0; // for neutral, select them in EMCAL acceptance
  Float_t ptMaxNeutPhot      = 0; // for neutral, take only photons
  Float_t ptMaxNeutEMCALPhot = 0; // for neutral, take only photons in EMCAL acceptance 
  
  leading[0] = 0;
  leading[1] = 0;
  leading[2] = 0;
  leading[3] = 0;
  
  isolated[0] = 0;
  isolated[1] = 0;
  isolated[2] = 0;
  isolated[3] = 0;
  
  Float_t ptTrig  = fTrigger.Pt();
  Float_t etaTrig = fTrigger.Eta();
  Float_t phiTrig = fTrigger.Phi();
  if(phiTrig < 0 ) phiTrig += TMath::TwoPi();

  // Minimum track or cluster energy

  //Isolation cuts
  Float_t ptThresIC    = GetIsolationCut()->GetPtThreshold();
  Float_t sumThresIC   = GetIsolationCut()->GetPtThreshold();
  Float_t rThresIC     = GetIsolationCut()->GetConeSize();
  Float_t isoMethod    = GetIsolationCut()->GetICMethod();
  
  // Counters
  Int_t   nICTrack     = 0;
  Int_t   nICNeutral   = 0;
  Int_t   nICNeutEMCAL = 0;
  Int_t   nICNeutPhot  = 0;
  Int_t   nICNeutEMCALPhot = 0;
  
  // Sum of pT
  Float_t sumNePt         = 0;
  Float_t sumNePtPhot     = 0;
  Float_t sumNePtEMC      = 0;
  Float_t sumNePtEMCPhot  = 0;
  Float_t sumChPt         = 0;
  
  // Loop on primaries, start from position 8, no partons
  
  Int_t imother = -1;
  Int_t pdg     =  0;
  Int_t status  =  0;
  Int_t charge  =  0;
  for(Int_t ipr = 8; ipr < fNPrimaries; ipr ++ )
  {
    if(ipr == indexTrig) continue;
    
    if     ( GetReader()->ReadStack() )
    {
      TParticle * particle = GetMC()->Particle(ipr) ;
      
      imother = particle->GetFirstMother();
      pdg     = particle->GetPdgCode();
      status  = particle->GetStatusCode();
      
      // Compare trigger with final state particles
      if( status != 1) continue ; // do it here to avoid crashes
      
      charge  = (Int_t) TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
      particle->Momentum(fLVTmp);
    }
    else if( GetReader()->ReadAODMCParticles() )
    {
      AliAODMCParticle * particle = (AliAODMCParticle*) GetMC()->GetTrack(ipr) ;
      
      imother = particle->GetMother();
      pdg     = particle->GetPdgCode();
      status  = particle->GetStatus();
      
      // Compare trigger with final state particles
      if( status != 1) continue ; // do it here to avoid crashes
      
      charge  = particle->Charge();
      fLVTmp.SetPxPyPzE(particle->Px(),particle->Py(),particle->Pz(),particle->E());
    }
    
    // Do not consider the photon decays from pi0 and eta
    //printf("Leading ipr %d - mother %d - iTrig\n",ipr, imother,indexTrig);
    if( imother == indexTrig)  continue ;
    
    Float_t pt  = fLVTmp.Pt();
    Float_t eta = fLVTmp.Eta();
    Float_t phi = fLVTmp.Phi();
    if(phi < 0 ) phi += TMath::TwoPi();
    
    // Select all particles in at least the TPC acceptance
    Bool_t inTPC = GetFiducialCut()->IsInFiducialCut(eta,phi,kCTS) ;
    if(!inTPC) continue;
    
    //Isolation
    Double_t radius = GetIsolationCut()->Radius(etaTrig, phiTrig, eta , phi) ;
    
    if(charge==0)
    {
      if(pt < fMinNeutralPt)  continue ;
      
      if( ptMaxNeutral < pt ) ptMaxNeutral = pt;
      
      if( radius < rThresIC )
      {
        if( pt > ptThresIC ) nICNeutral++ ;
        sumNePt+= pt;
      }
      
      Bool_t phPDG = kFALSE;
      if(pdg==22 || pdg==111) phPDG = kTRUE;
    
//      if(pt > ptTrig) printf(" --- pdg %d, phPDG %d pT %2.2f, pTtrig %2.2f, eta %2.2f, phi %2.2f\n",pdg,phPDG,pt,ptTrig,eta, phi*TMath::RadToDeg());
      if(phPDG)
      {
        if( ptMaxNeutPhot < pt) ptMaxNeutPhot = pt;
        
        if( radius < rThresIC )
        {
          if(pt > ptThresIC) nICNeutPhot++ ;
          sumNePtPhot += pt;
        }
      }
      
      // Calorimeter acceptance
      Bool_t inCalo = GetFiducialCut()->IsInFiducialCut(eta,phi,GetCalorimeter()) ;
      if(!inCalo) continue;
      
      if( ptMaxNeutEMCAL < pt ) ptMaxNeutEMCAL = pt;
      if( radius < rThresIC )
      {
        if( pt > ptThresIC ) nICNeutEMCAL++ ;
        sumNePtEMC += pt;
      }
      
      if(phPDG)
      {
        if( ptMaxNeutEMCALPhot < pt ) ptMaxNeutEMCALPhot = pt;
        if(  radius < rThresIC )
        {
          if (pt > ptThresIC) nICNeutEMCALPhot++ ;
          sumNePtEMCPhot += pt;
        }
      }
    }
    else
    {
      if( pt < fMinChargedPt)  continue ;
      
      if( ptMaxCharged < pt )   ptMaxCharged   = pt;
      
      if( radius < rThresIC )
      {
//        printf("UE track? pTtrig %2.2f, pt %2.2f, etaTrig %2.2f,  eta %2.2f, phiTrig %2.2f,  phi %2.2f, radius %2.2f\n",
//               ptTrig, pt,etaTrig, eta, phiTrig*TMath::RadToDeg(), phi*TMath::RadToDeg(),radius);
        if( pt > ptThresIC ) nICTrack++ ;
        sumChPt += pt;
      }
    }
  } // particle loop
  
  // Leding decision
  if(ptTrig > ptMaxCharged)
  {
//    printf("pt charged %2.2f, pt neutral %2.2f, pt neutral emcal %2.2f, pt photon %2.2f, pt photon emcal %2.2f\n",
//           ptMaxCharged, ptMaxNeutral, ptMaxNeutEMCAL,ptMaxNeutPhot, ptMaxNeutEMCALPhot);
    if(ptTrig > ptMaxNeutral      ) leading[0] = kTRUE ;
    if(ptTrig > ptMaxNeutEMCAL    ) leading[1] = kTRUE ;
    if(ptTrig > ptMaxNeutPhot     ) leading[2] = kTRUE ;
    if(ptTrig > ptMaxNeutEMCALPhot) leading[3] = kTRUE ;
  }
  
//  printf("N in cone over threshold: tracks  %d, neutral %d, neutral emcal %d, photon %d, photon emcal %d\n",
//         nICTrack, nICNeutral ,nICNeutEMCAL,nICNeutPhot, nICNeutEMCALPhot);
  
  //------------------
  // Isolation decision
  if( isoMethod == AliIsolationCut::kPtThresIC )
  {
    if( nICTrack == 0 )
    {
      if(nICNeutral       == 0 ) isolated[0] = kTRUE ;
      if(nICNeutEMCAL     == 0 ) isolated[1] = kTRUE ;
      if(nICNeutPhot      == 0 ) isolated[2] = kTRUE ;
      if(nICNeutEMCALPhot == 0 ) isolated[3] = kTRUE ;
    }
  }
  else if( isoMethod == AliIsolationCut::kSumPtIC )
  {
    if(sumChPt + sumNePt        < sumThresIC ) isolated[0] = kTRUE ;
    if(sumChPt + sumNePtEMC     < sumThresIC ) isolated[1] = kTRUE ;
    if(sumChPt + sumNePtPhot    < sumThresIC ) isolated[2] = kTRUE ;
    if(sumChPt + sumNePtEMCPhot < sumThresIC ) isolated[3] = kTRUE ;
  }
  
  //----------------------------------------------------
  // Fill histograms if conditions apply for all 4 cases
  for( Int_t i = 0; i < fgkNIso; i++ )
  {
    if(leading[i])
    {
      fhPtLeading[partType][i]->Fill(ptTrig, GetEventWeight());
      
      if     (i == 0) fhPtLeadingSumPt[partType][i]->Fill(ptTrig, sumChPt + sumNePt       , GetEventWeight());
      else if(i == 1) fhPtLeadingSumPt[partType][i]->Fill(ptTrig, sumChPt + sumNePtEMC    , GetEventWeight());
      else if(i == 2) fhPtLeadingSumPt[partType][i]->Fill(ptTrig, sumChPt + sumNePtPhot   , GetEventWeight());
      else if(i == 3) fhPtLeadingSumPt[partType][i]->Fill(ptTrig, sumChPt + sumNePtEMCPhot, GetEventWeight());
      
      if(isolated[i]) fhPtLeadingIsolated[partType][i]->Fill(ptTrig, GetEventWeight());
    }
  } // conditions loop
 
  AliDebug(1,"End");
}
  
//_____________________________________________________
/// Particle-Parton/Jet/Hadron Correlation Analysis, main method.
//_____________________________________________________
void  AliAnaGeneratorKine::MakeAnalysisFillHistograms()
{
  if( !GetMC() )
  {
    AliFatal("MCEvent not available, is the MC handler called? STOP");
    return;
  }
  
  AliDebug(1,"Start");
  
  fParton6.SetPxPyPzE(0,0,0,0);
  fParton7.SetPxPyPzE(0,0,0,0);
  fParton6PDG = 0;
  fParton7PDG = 0;

  AliAODMCParticle * primAOD = 0;
  
  fNPrimaries = GetMC()->GetNumberOfPrimaries(); // GetNtrack();

  //
  // Get the ESD MC particles container
  if( GetReader()->ReadStack() )
  {
    if(fNPrimaries > 6)
    {
      (GetMC()->Particle(6))->Momentum(fParton6) ;
      fParton6PDG =  (GetMC()->Particle(6))->GetPdgCode();
    }
    
    if(fNPrimaries > 7)
    {
      (GetMC()->Particle(7))->Momentum(fParton7) ;
      fParton7PDG =  (GetMC()->Particle(7))->GetPdgCode();
    }
  }
  else if( GetReader()->ReadAODMCParticles() )
  {
    if(fNPrimaries > 6)
    {
      primAOD = (AliAODMCParticle *) GetMC()->GetTrack(6);
      fParton6.SetPxPyPzE(primAOD->Px(),primAOD->Py(),primAOD->Pz(),primAOD->E());
      
      fParton6PDG =  primAOD->GetPdgCode();
    }
    
    if(fNPrimaries > 7)
    {
      primAOD = (AliAODMCParticle *) GetMC()->GetTrack(7);
      fParton7.SetPxPyPzE(primAOD->Px(),primAOD->Py(),primAOD->Pz(),primAOD->E());
      
      fParton7PDG =  primAOD->GetPdgCode();
    }
  }

  GetPartonsAndJets();
  
  // Main particle loop
  Int_t   pdgTrig    = 0;
  Int_t   statusTrig = 0;
  Int_t   imother    = 0;
  Float_t ptTrig     = 0;
  Int_t   momStatus  = 0;
  Int_t   momPdg     = 0;
  Int_t   momNDaugh  = 0;
  Int_t   momImom    = 0;
  Int_t   pdg0       = 0;
  Int_t   pdg1       = 0;
  Int_t   id0        = 0;
  Int_t   id1        = 0;
  Int_t   nDaughters = 0;
  
  for(Int_t ipr = 0; ipr < fNPrimaries; ipr ++ )
  {
    if( GetReader()->ReadStack() )
    {
      TParticle * particle = GetMC()->Particle(ipr) ;
      
      pdgTrig    = particle->GetPdgCode();
      statusTrig = particle->GetStatusCode();
      imother    = particle->GetFirstMother();
      nDaughters = particle->GetNDaughters();
      id0        = particle->GetDaughter(0);
      id1        = particle->GetDaughter(1);
      // Recover the kinematics:
      particle->Momentum(fTrigger);
    }
    else if( GetReader()->ReadAODMCParticles() )
    {
      AliAODMCParticle* particle = (AliAODMCParticle*) GetMC()->GetTrack(ipr) ;
      
      pdgTrig    = particle->GetPdgCode();
      statusTrig = particle->GetStatus();
      imother    = particle->GetMother();
      nDaughters = particle->GetNDaughters();
      id0        = particle->GetDaughter(0);
      id1        = particle->GetDaughter(1);
      // Recover the kinematics:
      fTrigger.SetPxPyPzE(particle->Px(),particle->Py(),particle->Pz(),particle->E());
    }
    
    //
    // Select final state photons, pi0s or eta's (mesons not final state)
    //
    // Final state particles have status code 1
    // pi0 and etas can have status code 0, (not decayed by pythia but geant?)
    // In the current code photons are always final state with status code 1,
    // avoid complications with conversions although pi0 with status code 0
    // generate photons with status code 0
    
    if( pdgTrig == 22 && statusTrig != 1 ) continue ;
    
    if( pdgTrig != 111 && pdgTrig != 22 && pdgTrig !=221 ) continue ;
    
    //
    // Pt cut
    //
    ptTrig  = fTrigger.Pt();
    
    if( ptTrig < GetMinPt() ) continue ;
    
    //
    // Select and tag the particles being, direct photon (prompt, fragmentation or isr)
    // decay photon from pi0, eta or other, and pi0 or eta
    //
    Int_t partType = -1;
    
    if     (pdgTrig==22 )
    {
      if(imother > 0 )
      {
        if( GetReader()->ReadStack() )
        {
          momStatus = (GetMC()->Particle(imother))->GetStatusCode();
          momPdg    = (GetMC()->Particle(imother))->GetPdgCode();
          momNDaugh = (GetMC()->Particle(imother))->GetNDaughters();
          momImom   = (GetMC()->Particle(imother))->GetFirstMother();
        }
        else if( GetReader()->ReadAODMCParticles() )
        {
          momStatus = ((AliAODMCParticle*) GetMC()->GetTrack(imother))->GetStatus();
          momPdg    = ((AliAODMCParticle*) GetMC()->GetTrack(imother))->GetPdgCode();
          momNDaugh = ((AliAODMCParticle*) GetMC()->GetTrack(imother))->GetNDaughters();
          momImom   = ((AliAODMCParticle*) GetMC()->GetTrack(imother))->GetMother();
        }
        
        if     (imother < 8 && statusTrig == 1)
        {
          partType = kmcPrimPhoton ;
        }
        else if(momPdg == 111 )
        {
          partType = kmcPrimPi0Decay;
        }
        else if(momPdg == 221 )
        {
          partType = kmcPrimEtaDecay;
        }
        else if(TMath::Abs(momStatus) > 0 )
        {
          partType = kmcPrimOtherDecay ;
        }
      }
    }
    else if( (pdgTrig==111 || pdgTrig==221) && nDaughters == 2 )
    {
      if( GetReader()->ReadStack() )
      {
        pdg0 = GetMC()->Particle(id0)->GetPdgCode();
        pdg1 = GetMC()->Particle(id1)->GetPdgCode();
      }
      else if( GetReader()->ReadAODMCParticles() )
      {
        pdg0 = ((AliAODMCParticle*) GetMC()->GetTrack(id0))->GetPdgCode();
        pdg1 = ((AliAODMCParticle*) GetMC()->GetTrack(id1))->GetPdgCode();
      }
      
      if( pdg0 == 22 && pdg1== 22 )
      {
        if     ( pdgTrig==111 ) partType = kmcPrimPi0;
        else if( pdgTrig==221 ) partType = kmcPrimEta;
      }
    }
    else if( (pdgTrig==111 || pdgTrig==221) )
    {
      // Fill histogram to see how many pi0/eta decay other than 2 photons in trigger detector acceptance
      Bool_t in = GetFiducialCutForTrigger()->IsInFiducialCut(fTrigger.Eta(),fTrigger.Phi(),fTriggerDetector) ;
      if(! in )  continue ;
      
      if(pdgTrig==111) fhPtPi0Not2Gamma->Fill(ptTrig, GetEventWeight());
      if(pdgTrig==221) fhPtEtaNot2Gamma->Fill(ptTrig, GetEventWeight());
    }
    
    if(partType < 0 ) continue ;
    
    //
    // Fill particle acceptance histograms
    //
    Float_t eta = fTrigger.Eta();
    Float_t phi = fTrigger.Phi();
    if(phi < 0) phi+=TMath::TwoPi();
    
    fhPhi   [partType]->Fill(phi,      GetEventWeight());
    fhEta   [partType]->Fill(eta,      GetEventWeight());
    fhEtaPhi[partType]->Fill(eta, phi, GetEventWeight());
    
    if(partType < 4 &&  partType!=0)
    {
      fhPhiStatus[partType]->Fill(phi, momStatus, GetEventWeight());
      fhEtaStatus[partType]->Fill(eta, momStatus, GetEventWeight());
    }
    else
    {
      fhPhiStatus[partType]->Fill(phi, statusTrig, GetEventWeight());
      fhEtaStatus[partType]->Fill(eta, statusTrig, GetEventWeight());
    }
    
    //
    // Select particles in trigger detector acceptance
    //
    Bool_t in = GetFiducialCutForTrigger()->IsInFiducialCut(eta,phi,fTriggerDetector) ;
    if(! in )  continue ;
    
    
    AliDebug(1,Form("Select trigger particle %d: pdg %d, type %d, status %d, mother index %d, pT %2.2f, eta %2.2f, phi %2.2f",
                    ipr, pdgTrig, partType, statusTrig, imother, ptTrig, fTrigger.Eta(), fTrigger.Phi()*TMath::RadToDeg()));
    
    //
    // Fill particle pT histograms, check also status
    //
    
    fhPt[partType]->Fill(ptTrig, GetEventWeight());
    
    if     (partType==kmcPrimPi0)
    {
      fhPtPi0Status->Fill(ptTrig, statusTrig, GetEventWeight());
    }
    else if(partType==kmcPrimEta)
    {
      fhPtEtaStatus->Fill(ptTrig, statusTrig, GetEventWeight());
    }
    else if(partType == kmcPrimPi0Decay )
    {
      fhPtPi0DecayStatus  ->Fill(ptTrig, momStatus, GetEventWeight());
      
      if(momNDaugh!=2) fhPtGammaFromPi0Not2Gamma->Fill(ptTrig, GetEventWeight());
    }
    else if(partType == kmcPrimEtaDecay )
    {
      fhPtEtaDecayStatus  ->Fill(ptTrig, momStatus, GetEventWeight());
      
      if(momNDaugh!=2) fhPtGammaFromEtaNot2Gamma->Fill(ptTrig, GetEventWeight());
    }
    else if(partType == kmcPrimOtherDecay)
    {
      fhPtOtherDecayStatus->Fill(ptTrig, momStatus, GetEventWeight());
      
      //Fill histogram with origin of this kind of photon
      if     (momPdg     > 2100  && momPdg   < 2210)
                                fhPtOtherDecayMesonId->Fill(ptTrig, 0.5, GetEventWeight());// resonances
      else if(momPdg    == 221) fhPtOtherDecayMesonId->Fill(ptTrig, 5.5, GetEventWeight());//eta
      else if(momPdg    == 331) fhPtOtherDecayMesonId->Fill(ptTrig, 6.5, GetEventWeight());//eta prime
      else if(momPdg    == 213) fhPtOtherDecayMesonId->Fill(ptTrig, 1.5, GetEventWeight());//rho
      else if(momPdg    == 223) fhPtOtherDecayMesonId->Fill(ptTrig, 2.5, GetEventWeight());//omega
      else if(momPdg    >= 310   && momPdg    <= 323)
                                fhPtOtherDecayMesonId->Fill(ptTrig, 3.5, GetEventWeight());//k0S, k+-,k*
      else if(momPdg    == 130) fhPtOtherDecayMesonId->Fill(ptTrig, 3.5, GetEventWeight());//k0L
      else                      fhPtOtherDecayMesonId->Fill(ptTrig, 4.5, GetEventWeight());//other?
    }
    
    //
    // Origin of particle, which mother of the pi0 or eta or
    // the eta or pi0 or other meson at the origin of the decay photon
    //
    Int_t momOrgPdg    = -1;
    Int_t momOrgStatus = -1;
    if( GetReader()->ReadStack() )
    {
      if(momImom >=0 || imother >=0)
      {
        TParticle* mother = 0;
        if(partType < 4 && partType!=0 )
          mother = GetMC()->Particle(momImom);
        else
          mother = GetMC()->Particle(imother);
        
        momOrgPdg    = TMath::Abs(mother->GetPdgCode());
        momOrgStatus = mother->GetStatusCode();
      }
    }
    else if( GetReader()->ReadAODMCParticles() )
    {
      if(momImom >=0 || imother >=0)
      {
        AliAODMCParticle* mother = 0;
        if(partType < 4 && partType!=0 )
          mother = (AliAODMCParticle*) GetMC()->GetTrack(momImom);
        else
          mother = (AliAODMCParticle*) GetMC()->GetTrack(imother);
        
        momOrgPdg    = TMath::Abs(mother->GetPdgCode());
        momOrgStatus = mother->GetStatus();
      }
    }
    
    if     (momOrgStatus  == 21) fhPtOrigin[partType]->Fill(ptTrig, 0.5, GetEventWeight());//parton
    else if(momOrgPdg     < 22 ) fhPtOrigin[partType]->Fill(ptTrig, 1.5, GetEventWeight());//quark
    else if(momOrgPdg     > 2100  && momOrgPdg   < 2210)
                                 fhPtOrigin[partType]->Fill(ptTrig, 2.5, GetEventWeight());// resonances
    else if(momOrgPdg    == 221) fhPtOrigin[partType]->Fill(ptTrig, 8.5, GetEventWeight());//eta
    else if(momOrgPdg    == 331) fhPtOrigin[partType]->Fill(ptTrig, 9.5, GetEventWeight());//eta prime
    else if(momOrgPdg    == 213) fhPtOrigin[partType]->Fill(ptTrig, 4.5, GetEventWeight());//rho
    else if(momOrgPdg    == 223) fhPtOrigin[partType]->Fill(ptTrig, 5.5, GetEventWeight());//omega
    else if(momOrgPdg    >= 310   && momOrgPdg    <= 323)
                                 fhPtOrigin[partType]->Fill(ptTrig, 6.5, GetEventWeight());//k0S, k+-,k*
    else if(momOrgPdg    == 130) fhPtOrigin[partType]->Fill(ptTrig, 6.5, GetEventWeight());//k0L
    else if(momOrgStatus == 11 || momOrgStatus  == 12 )
                                 fhPtOrigin[partType]->Fill(ptTrig, 3.5, GetEventWeight());//resonances
    else                         fhPtOrigin[partType]->Fill(ptTrig, 7.5, GetEventWeight());//other?
    
    if(statusTrig == 0)
    {
      // Histogram will not be filled for photons, leave it like this for now
      // in case we leave not final photons in the future
      if     (momOrgStatus  == 21) fhPtOriginNotFinal[partType]->Fill(ptTrig, 0.5, GetEventWeight());//parton
      else if(momOrgPdg     < 22 ) fhPtOriginNotFinal[partType]->Fill(ptTrig, 1.5, GetEventWeight());//quark
      else if(momOrgPdg     > 2100  && momOrgPdg   < 2210)
                                   fhPtOriginNotFinal[partType]->Fill(ptTrig, 2.5, GetEventWeight());// resonances
      else if(momOrgPdg    == 221) fhPtOriginNotFinal[partType]->Fill(ptTrig, 8.5, GetEventWeight());//eta
      else if(momOrgPdg    == 331) fhPtOriginNotFinal[partType]->Fill(ptTrig, 9.5, GetEventWeight());//eta prime
      else if(momOrgPdg    == 213) fhPtOriginNotFinal[partType]->Fill(ptTrig, 4.5, GetEventWeight());//rho
      else if(momOrgPdg    == 223) fhPtOriginNotFinal[partType]->Fill(ptTrig, 5.5, GetEventWeight());//omega
      else if(momOrgPdg    >= 310   && momOrgPdg    <= 323)
                                   fhPtOriginNotFinal[partType]->Fill(ptTrig, 6.5, GetEventWeight());//k0S, k+-,k*
      else if(momOrgPdg    == 130) fhPtOriginNotFinal[partType]->Fill(ptTrig, 6.5, GetEventWeight());//k0L
      else if(momOrgStatus == 11 || momOrgStatus  == 12 )
                                   fhPtOriginNotFinal[partType]->Fill(ptTrig, 3.5, GetEventWeight());//resonances
      else                         fhPtOriginNotFinal[partType]->Fill(ptTrig, 7.5, GetEventWeight());//other?
    }
    
    //
    // Check if it is leading or isolated
    //
    Bool_t leading [fgkNIso] ;
    Bool_t isolated[fgkNIso] ;
    
    IsLeadingAndIsolated(ipr, partType, leading, isolated);
    
    //
    // Correlate trigger particle with partons or jets
    //
    Int_t iparton = -1;
    Int_t ok = CorrelateWithPartonOrJet(ipr, partType, leading, isolated, iparton);
    if(!ok) continue;
    
    //
    // Correlate trigger particle with hadrons
    //
    GetXE(ipr,partType,leading,isolated,iparton) ;
    
  }
  
  AliDebug(1,"End fill histograms");
}

//_________________________________________________________
/// Set the calorimeter for the analysis.
//_________________________________________________________
void AliAnaGeneratorKine::SetTriggerDetector(TString & det)
{
  fTriggerDetectorString = det;
  
  if     (det=="EMCAL") fTriggerDetector = kEMCAL;
  else if(det=="PHOS" ) fTriggerDetector = kPHOS;
  else if(det=="CTS")   fTriggerDetector = kCTS;
  else if(det=="DCAL")  fTriggerDetector = kDCAL;
  else if(det.Contains("DCAL") && det.Contains("PHOS")) fTriggerDetector = kDCALPHOS;
  else AliFatal(Form("Detector < %s > not known!", det.Data()));
}

//_____________________________________________________
/// Set the detrimeter for the analysis.
//_____________________________________________________
void AliAnaGeneratorKine::SetTriggerDetector(Int_t det)
{
  fTriggerDetector = det;
  
  if     (det==kEMCAL)    fTriggerDetectorString = "EMCAL";
  else if(det==kPHOS )    fTriggerDetectorString = "PHOS";
  else if(det==kCTS)      fTriggerDetectorString = "CTS";
  else if(det==kDCAL)     fTriggerDetectorString = "DCAL";
  else if(det==kDCALPHOS) fTriggerDetectorString = "DCAL_PHOS";
  else AliFatal(Form("Detector < %d > not known!", det));
}

