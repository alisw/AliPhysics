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

//_________________________________________________________________________
// Do photon/pi0 analysis for isolation and correlation
// at the generator level. Only for kine stack (ESDs)
//
//
// -- Author: Gustavo Conesa (LPSC-CNRS-Grenoble) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TH2F.h"
#include "TParticle.h"
#include "TDatabasePDG.h"

//---- ANALYSIS system ----
#include "AliAnaGeneratorKine.h" 
#include "AliStack.h"  
#include "AliGenPythiaEventHeader.h"

ClassImp(AliAnaGeneratorKine)


//__________________________________________
AliAnaGeneratorKine::AliAnaGeneratorKine() : 
AliAnaCaloTrackCorrBaseClass(), 
fTriggerDetector(),  fTriggerDetectorString(),
fFidCutTrigger(0),
fMinChargedPt(0),    fMinNeutralPt(0),
fStack(0),
//fParton2(0),         fParton3(0),
fParton6(0),         fParton7(0),
fJet6(),             fJet7(),
fTrigger(),          fLVTmp(),
fPtHard(0),
fhPtHard(0),         fhPtParton(0),    fhPtJet(0),
fhPtPartonPtHard(0), fhPtJetPtHard(0), fhPtJetPtParton(0)
{
  //Default Ctor
  
  //Initialize parameters
  InitParameters();
  
 

  for(Int_t p = 0; p < fgkNmcPrimTypes; p++)
  {
    fhPt[p] = 0;
    
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
Bool_t  AliAnaGeneratorKine::CorrelateWithPartonOrJet(Int_t   indexTrig,
                                                      Int_t   partType,
                                                      Bool_t  leading [fgkNIso],
                                                      Bool_t  isolated[fgkNIso],
                                                      Int_t & iparton )  
{
  // Correlate trigger with partons or jets, get z
  
  AliDebug(1,"Start");
  
  if( fStack->GetNprimary() < 7 ) return kFALSE ;
  
  //Get the index of the mother
  iparton =  (fStack->Particle(indexTrig))->GetFirstMother();
  TParticle * mother = fStack->Particle(iparton);
  while (iparton > 7) 
  {
    iparton   = mother->GetFirstMother();
    if(iparton < 0)
    {
      AliWarning("Negative index, skip event");
      return kFALSE;
    }
    mother = fStack->Particle(iparton);
  }
  
  //printf("Mother is parton %d with pdg %d\n",imom,mother->GetPdgCode());
  
  if(iparton < 6)
  {
    //printf("This particle is not from hard process - pdg %d, parton index %d\n",partType, iparton);
    return kFALSE; 
  }
  
  Float_t  ptTrig  = fTrigger.Pt();
  Float_t phiTrig  = fTrigger.Phi();
  Float_t etaTrig  = fTrigger.Eta();
  
  Float_t partonPt = fParton6->Pt();
  
  Float_t jetPt    = fJet6.Pt();
  Float_t jetPhi   = fJet6.Phi();
  Float_t jetEta   = fJet6.Eta();
  
  Float_t awayJetPt  = fJet7.Pt();
  Float_t awayJetEta = fJet7.Eta();
  Float_t awayJetPhi = fJet7.Phi();
  
  Int_t nearPDG = fParton6->GetPdgCode();
  Int_t awayPDG = fParton7->GetPdgCode();
  
  if ( iparton == 7 )
  {
    partonPt = fParton7->Pt();
    
    jetPt  = fJet7.Pt();
    jetPhi = fJet7.Phi();
    jetEta = fJet7.Eta();
    
    awayJetPt  = fJet6.Pt();
    awayJetEta = fJet6.Eta();
    awayJetPhi = fJet6.Phi();

    nearPDG = fParton7->GetPdgCode();
    awayPDG = fParton6->GetPdgCode();
    
  }

  Float_t deltaPhi = TMath::Abs(phiTrig-awayJetPhi) *TMath::RadToDeg();
  //printf("Trigger Away jet phi %2.2f\n",deltaPhi);
  
  //Get id of parton in near and away side
  Int_t away = -1;
  Int_t near = -1;
  if     (nearPDG == 22) near = 0;
  else if(nearPDG == 21) near = 1;
  else                   near = 2;
  
  if     (awayPDG == 22) away = 0;
  else if(awayPDG == 21) away = 1;
  else                   away = 2;
  //printf("parton 6 pdg = %d, parton 7 pdg = %d\n",fParton6->GetPdgCode(),fParton7->GetPdgCode());

  AliDebug(1,Form("Trigger pdg %d, (Trigger,JetNear, JetAway): pT (%2.2f,%2.2f,%2.2f), phi (%2.2f,%2.2f,%2.2f), eta (%2.2f,%2.2f,%2.2f); Delta phi trigger-away jet %f; parton type: away side %d; near side %d",
                  partType,
                  ptTrig,jetPt,awayJetPt,
                  phiTrig*TMath::RadToDeg(),jetPhi*TMath::RadToDeg(),awayJetPhi*TMath::RadToDeg(),
                  etaTrig,jetEta,awayJetEta,
                  deltaPhi,away,near));
  
  // RATIOS
  
  Float_t zHard = ptTrig / fPtHard;
  Float_t zPart = ptTrig / partonPt;
  Float_t zJet  = ptTrig / jetPt;
  
  //if(zHard > 1 ) printf("*** Particle energy larger than pT hard z=%f\n",zHard);
  
  //printf("Z: hard %2.2f, parton %2.2f, jet %2.2f\n",zHard,zPart,zJet);

  // conditions loop
  for( Int_t i = 0; i < fgkNIso ; i++ )
  {
    fhPtPartonTypeNear[partType][leading[i]][i]->Fill(ptTrig,near);
    fhPtPartonTypeAway[partType][leading[i]][i]->Fill(ptTrig,away);
    
    fhZHard  [partType][leading[i]][i]->Fill(ptTrig,zHard);
    fhZParton[partType][leading[i]][i]->Fill(ptTrig,zPart);
    fhZJet   [partType][leading[i]][i]->Fill(ptTrig,zJet );
    
    if(isolated[i])
    {
      fhPtPartonTypeNearIsolated[partType][leading[i]][i]->Fill(ptTrig,near);
      fhPtPartonTypeAwayIsolated[partType][leading[i]][i]->Fill(ptTrig,away);
      
      fhZHardIsolated  [partType][leading[i]][i]->Fill(ptTrig,zHard);
      fhZPartonIsolated[partType][leading[i]][i]->Fill(ptTrig,zPart);
      fhZJetIsolated   [partType][leading[i]][i]->Fill(ptTrig,zJet);
    }
    
    if(partType == kmcPrimPhoton && deltaPhi < 220 && deltaPhi > 140 && TMath::Abs(awayJetEta) < 0.6)
    {
      //printf("Accept jet\n");
      fhPtAcceptedGammaJet[leading[i]][i]->Fill(ptTrig,away);
    }
    //else printf("Reject jet!!!\n");
    
  } // conditions loop
  
  AliDebug(1,"End TRUE");
  
  return kTRUE;
}


//____________________________________________________
TList *  AliAnaGeneratorKine::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file 
  
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
void  AliAnaGeneratorKine::GetPartonsAndJets() 
{
  // Fill data members with partons,jets and generated pt hard 
  
  AliDebug(1,"Start");
  
  Int_t nPrimary = fStack->GetNprimary();
  
//  if( nPrimary > 2 ) fParton2 = fStack->Particle(2) ;
//  if( nPrimary > 3 ) fParton3 = fStack->Particle(3) ;
  
  Float_t p6phi = -1 ;
  Float_t p6eta = -10;
  Float_t p6pt  =  0 ;
  fParton6 = 0x0;
  if( nPrimary > 6 )
  {
    fParton6 = fStack->Particle(6) ;
    p6pt  = fParton6->Pt();
    p6eta = fParton6->Eta();
    p6phi = fParton6->Phi();
    if(p6phi < 0) p6phi +=TMath::TwoPi();
  }
  
  Float_t p7phi = -1 ;
  Float_t p7eta = -10;
  Float_t p7pt  =  0 ;
  fParton7 = 0x0;
  if( nPrimary > 7 )
  {
    fParton7 = fStack->Particle(7) ;
    p7pt  = fParton7->Pt();
    p7phi = fParton7->Eta();
    p7phi = fParton7->Phi();
    if(p7phi < 0) p7phi +=TMath::TwoPi();
  }


  //printf("parton6: pt %2.2f, eta %2.2f, phi %2.2f with pdg %d\n",fParton6->Pt(),fParton6->Eta(),p6phi, fParton6->GetPdgCode());
  //printf("parton7: pt %2.2f, eta %2.2f, phi %2.2f with pdg %d\n",fParton7->Pt(),fParton7->Eta(),p7phi, fParton7->GetPdgCode());
  
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
  
  fhPtHard   ->Fill(fPtHard);
  fhPtJet    ->Fill(fJet6.Pt());
  fhPtJet    ->Fill(fJet7.Pt());
  fhPtParton ->Fill(p6pt);
  fhPtParton ->Fill(p7pt);

  fhPtPartonPtHard->Fill(fPtHard, p6pt/fPtHard);
  fhPtPartonPtHard->Fill(fPtHard, p7pt/fPtHard);
  fhPtJetPtHard   ->Fill(fPtHard, fJet6.Pt()/fPtHard);
  fhPtJetPtHard   ->Fill(fPtHard, fJet7.Pt()/fPtHard);
  if( p6pt > 0 ) fhPtJetPtParton ->Fill(fPtHard, fJet6.Pt()/p6pt);
  if( p7pt > 0 ) fhPtJetPtParton ->Fill(fPtHard, fJet7.Pt()/p7pt);
  
  AliDebug(1,"End");

}

//_____________________________________________________
void AliAnaGeneratorKine::GetXE(Int_t   indexTrig,
                                Int_t   partType,
                                Bool_t  leading [fgkNIso],
                                Bool_t  isolated[fgkNIso],
                                Int_t   iparton)
{

  // Calculate the real XE and the UE XE

  AliDebug(1,"Start");
  
  Float_t ptTrig  = fTrigger.Pt();
  Float_t phiTrig = fTrigger.Phi();
  if(phiTrig < 0 ) phiTrig += TMath::TwoPi();
  
  //Loop on primaries, start from position 8, no partons
  for(Int_t ipr = 8; ipr < fStack->GetNprimary(); ipr ++ )
  {
    TParticle * particle = fStack->Particle(ipr) ;
    
    if(ipr==indexTrig) continue;
    
    Int_t   pdg    = particle->GetPdgCode();
    Int_t   status = particle->GetStatusCode();
        
    // Compare trigger with final state particles
    if( status != 1 ) continue ;
    
    Double_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
    
    if(charge == 0 ) continue; // construct xe only with charged
    
    Float_t pt     = particle->Pt();
    Float_t phi    = particle->Phi();
    if(phi < 0 ) phi += TMath::TwoPi();
    
    if( pt < fMinChargedPt)    continue ;
    
    particle->Momentum(fLVTmp);
    Bool_t inTPC = GetFiducialCut()->IsInFiducialCut(fLVTmp.Eta(),fLVTmp.Phi(),kCTS) ;
    
    if(!inTPC) continue;
    
    Float_t xe = -pt/ptTrig*TMath::Cos(phi-phiTrig);
    
    // ---------------------------------------------------
    // Get the index of the mother, get from what parton
    Int_t ipartonAway =  particle->GetFirstMother();
    if(ipartonAway < 0) return;
    
    TParticle * mother = fStack->Particle(ipartonAway);
    while (ipartonAway > 7) 
    {
      ipartonAway   = mother->GetFirstMother();
      if(ipartonAway < 0) break;
      mother = fStack->Particle(ipartonAway);
    }
    
    //-----------------------------------------
    // Get XE of particles belonging to the jet
    // on the opposite side of the trigger
    if((ipartonAway==6 || ipartonAway==7) && iparton!=ipartonAway) 
    {
      for( Int_t i = 0; i < fgkNIso; i++ )
      {
        fhXE[partType][leading[i]][i]  ->Fill(ptTrig,xe);
        
        if(isolated[i])
        {
          fhXEIsolated[partType][leading[i]][i]  ->Fill(ptTrig,xe);
        }
      } // conditions loop
    } // Away side

    //----------------------------------------------------------
    // Get the XE from particles not attached to any of the jets
    if(ipartonAway!=6 && ipartonAway!=7)
    {
      for( Int_t i = 0; i < fgkNIso; i++ )
      {
        fhXEUE[partType][leading[i]][i]  ->Fill(ptTrig,xe);
        
        if(isolated[i])
        {
          fhXEUEIsolated[partType][leading[i]][i]  ->Fill(ptTrig,xe);
        }
      } // conditions loop
    } // Away side    
    
  } // primary loop

  AliDebug(1,"End");

}

//________________________________________
void AliAnaGeneratorKine::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaGenKine_");
  
  fTriggerDetector = kEMCAL;
  
  fMinChargedPt    = 0.2;
  fMinNeutralPt    = 0.3;
  
}

//_____________________________________________________________________
void  AliAnaGeneratorKine::IsLeadingAndIsolated(Int_t indexTrig,
                                                Int_t partType,
                                                Bool_t leading[fgkNIso],
                                                Bool_t isolated[fgkNIso]) 
{
  // Check if the trigger is the leading particle and if it is isolated
  // In case of neutral particles check all neutral or neutral in EMCAL acceptance
  
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
  
  //Counters
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
  
  //Loop on primaries, start from position 8, no partons
  for(Int_t ipr = 8; ipr < fStack->GetNprimary(); ipr ++ )
  {
    if(ipr == indexTrig) continue;
    TParticle * particle = fStack->Particle(ipr) ;
    
    Int_t   imother = particle->GetFirstMother();
    //printf("Leading ipr %d - mother %d\n",ipr, imother);
    
    // Do not consider the photon decays from pi0 and eta
    if( imother == indexTrig)  continue ;
    
    Int_t   pdg    = particle->GetPdgCode();
    Int_t   status = particle->GetStatusCode();
     
    // Compare trigger with final state particles
    if( status != 1) continue ;

    // Select all particles in at least the TPC acceptance
    Bool_t inTPC = GetFiducialCut()->IsInFiducialCut(fTrigger.Eta(),fTrigger.Phi(),kCTS) ;
    if(!inTPC) continue;
    
    Float_t pt     = particle->Pt();
    Float_t eta    = particle->Eta();
    Float_t phi    = particle->Phi();
    if(phi < 0 ) phi += TMath::TwoPi();
    
    Double_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
    
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
    
      //if(pt > ptTrig) printf(" --- pdg %d, phPDG %d pT %2.2f, pTtrig %2.2f, eta %2.2f, phi %2.2f\n",pdg,phPDG,pt,ptTrig,particle->Eta(), particle->Phi()*TMath::RadToDeg());
      if(phPDG)
      {
        if( ptMaxNeutPhot < pt) ptMaxNeutPhot = pt;
        
        if( radius < rThresIC )
        {
          if(pt > ptThresIC) nICNeutPhot++ ;
          sumNePtPhot += pt;
        }
      }
      
      //Calorimeter acceptance
      Bool_t inCalo = GetFiducialCut()->IsInFiducialCut(fTrigger.Eta(),fTrigger.Phi(),GetCalorimeter()) ;
      if(!inCalo) continue;
      
      if( ptMaxNeutEMCAL < pt ) ptMaxNeutEMCAL = pt;
      if( radius < rThresIC )
      {
        if( pt > ptThresIC ) nICNeutEMCAL++ ;
        sumNePtEMC += pt;
      }
      
      if(phPDG)
      {
        if( ptMaxNeutEMCALPhot < pt )             ptMaxNeutEMCALPhot = pt;
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
        //printf("UE track? pTtrig %2.2f, pt %2.2f, etaTrig %2.2f,  eta %2.2f, phiTrig %2.2f,  phi %2.2f, radius %2.2f\n",
        //       ptTrig, pt,etaTrig, eta, phiTrig*TMath::RadToDeg(), phi*TMath::RadToDeg(),radius);
        if( pt > ptThresIC ) nICTrack++ ;
        sumChPt += pt;
      }
    }

  } // particle loop
  
  // Leding decision
  if(ptTrig > ptMaxCharged)
  {
    //printf("pt charged %2.2f, pt neutral %2.2f, pt neutral emcal %2.2f, pt photon %2.2f, pt photon emcal %2.2f\n", 
    //       ptMaxCharged, ptMaxNeutral, ptMaxNeutEMCAL,ptMaxNeutPhot, ptMaxNeutEMCALPhot);
    if(ptTrig > ptMaxNeutral      ) leading[0] = kTRUE ;
    if(ptTrig > ptMaxNeutEMCAL    ) leading[1] = kTRUE ;
    if(ptTrig > ptMaxNeutPhot     ) leading[2] = kTRUE ;
    if(ptTrig > ptMaxNeutEMCALPhot) leading[3] = kTRUE ;
  }
  
  //printf("N in cone over threshold: tracks  %d, neutral %d, neutral emcal %d, photon %d, photon emcal %d\n", 
  //       nICTrack, nICNeutral ,nICNeutEMCAL,nICNeutPhot, nICNeutEMCALPhot);
  
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
      fhPtLeading[partType][i]->Fill(ptTrig);
      
      if     (i == 0) fhPtLeadingSumPt[partType][i]->Fill(ptTrig, sumChPt + sumNePt);
      else if(i == 1) fhPtLeadingSumPt[partType][i]->Fill(ptTrig, sumChPt + sumNePtEMC);
      else if(i == 2) fhPtLeadingSumPt[partType][i]->Fill(ptTrig, sumChPt + sumNePtPhot);
      else if(i == 3) fhPtLeadingSumPt[partType][i]->Fill(ptTrig, sumChPt + sumNePtEMCPhot);
      
      if(isolated[i]) fhPtLeadingIsolated[partType][i]->Fill(ptTrig);
    }
  } // conditions loop
 
  AliDebug(1,"End");
  
}
  
//_____________________________________________________
void  AliAnaGeneratorKine::MakeAnalysisFillHistograms() 
{
  //Particle-Parton Correlation Analysis, fill histograms
  
  AliDebug(1,"Start");
  
  fStack =  GetMCStack() ;
  
  if( !fStack )
    AliFatal("No Stack available, STOP");

  GetPartonsAndJets();
  
  for(Int_t ipr = 0; ipr < fStack->GetNprimary(); ipr ++ )
  {
    TParticle * particle = fStack->Particle(ipr) ;
    
    Int_t   pdgTrig    = particle->GetPdgCode();
    Int_t   statusTrig = particle->GetStatusCode();
    Int_t   imother    = particle->GetFirstMother();
    Float_t ptTrig     = particle->Pt(); 

    // Acceptance and kinematical cuts
    if( ptTrig < GetMinPt() ) continue ;

    // Select final state photons or pi0s or eta's
    
    if( pdgTrig == 22 && statusTrig != 1 ) continue ;
    
    if( pdgTrig != 111 && pdgTrig != 22 && pdgTrig !=221 ) continue ;
    
    // Acceptance and kinematical cuts
    // Recover the kinematics:
    particle->Momentum(fTrigger);
    
    Bool_t in = GetFiducialCutForTrigger()->IsInFiducialCut(fTrigger.Eta(),fTrigger.Phi(),fTriggerDetector) ;
    if(! in )  continue ;

    // Identify the particle to fill appropriate histogram
    Int_t partType = -1;
    
    if     (pdgTrig==22 )
    {
      if(imother > 0 )
      {
        Int_t momStatus = (fStack->Particle(imother))->GetStatusCode();
        Int_t momPdg    = (fStack->Particle(imother))->GetPdgCode();
        
        if     (imother < 8 && statusTrig == 1)  partType = kmcPrimPhoton ;
        else if(momPdg == 111 ) partType = kmcPrimPi0Decay   ;
        else if(momPdg == 221 ) partType = kmcPrimEtaDecay   ;
        else if(momStatus > 0 ) partType = kmcPrimOtherDecay ;
      }
    }
    else if(pdgTrig==111 && particle->GetNDaughters() == 2)
    {
      Int_t pdg0 = fStack->Particle(particle->GetDaughter(0))->GetPdgCode();
      Int_t pdg1 = fStack->Particle(particle->GetDaughter(1))->GetPdgCode();
      
      if(pdg0 == 22 && pdg1== 22 ) partType = kmcPrimPi0;
    }
    else if(pdgTrig==221 && particle->GetNDaughters() == 2)
    {
      Int_t pdg0 = fStack->Particle(particle->GetDaughter(0))->GetPdgCode();
      Int_t pdg1 = fStack->Particle(particle->GetDaughter(1))->GetPdgCode();
      
      if(pdg0 == 22 && pdg1== 22 ) partType = kmcPrimEta;
    }
    
    if(partType < 0 ) continue ;

    AliDebug(1,Form("Select trigger particle %d: pdg %d, type %d, status %d, mother index %d, pT %2.2f, eta %2.2f, phi %2.2f",
                    ipr, pdgTrig, partType, statusTrig, imother, ptTrig, particle->Eta(), particle->Phi()*TMath::RadToDeg()));
    
    //    if(pdgTrig==111)
    //    {
    //      printf("\t pi0 daughters %d, %d\n", particle->GetDaughter(0), particle->GetDaughter(1));
    //    }

    // Fill histograms do analysis
    
    fhPt[partType]->Fill(ptTrig);
    
    // Check if it is leading
    Bool_t leading [fgkNIso] ;
    Bool_t isolated[fgkNIso] ;

    IsLeadingAndIsolated(ipr, partType, leading, isolated);
    
    Int_t iparton = -1;
    Int_t ok = CorrelateWithPartonOrJet(ipr, partType, leading, isolated, iparton);
    if(!ok) continue;
    
    GetXE(ipr,partType,leading,isolated,iparton) ;
    
  }
  
  AliDebug(1,"End fill histograms");
  
}

//_________________________________________________________
void AliAnaGeneratorKine::SetTriggerDetector(TString & det)
{
  // Set the detrimeter for the analysis
  
  fTriggerDetectorString = det;
  
  if     (det=="EMCAL") fTriggerDetector = kEMCAL;
  else if(det=="PHOS" ) fTriggerDetector = kPHOS;
  else if(det=="CTS")   fTriggerDetector = kCTS;
  else if(det=="DCAL")  fTriggerDetector = kDCAL;
  else if(det.Contains("DCAL") && det.Contains("PHOS")) fTriggerDetector = kDCALPHOS;
  else AliFatal(Form("Detector < %s > not known!", det.Data()));
  
}

//_____________________________________________________
void AliAnaGeneratorKine::SetTriggerDetector(Int_t det)
{
  // Set the detrimeter for the analysis
  
  fTriggerDetector = det;
  
  if     (det==kEMCAL)    fTriggerDetectorString = "EMCAL";
  else if(det==kPHOS )    fTriggerDetectorString = "PHOS";
  else if(det==kCTS)      fTriggerDetectorString = "CTS";
  else if(det==kDCAL)     fTriggerDetectorString = "DCAL";
  else if(det==kDCALPHOS) fTriggerDetectorString = "DCAL_PHOS";
  else AliFatal(Form("Detector < %d > not known!", det));
  
}

