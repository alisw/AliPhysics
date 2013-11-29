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
fStack(0),
fParton2(0),         fParton3(0), 
fParton6(0),         fParton7(0),   
fJet6(),             fJet7(),
fPtHard(0),
fhPtHard(0),         fhPtParton(0),    fhPtJet(0),
fhPtPartonPtHard(0), fhPtJetPtHard(0), fhPtJetPtParton(0),
fhPtPhoton(0),       fhPtPi0(0)
{
  //Default Ctor
  
  //Initialize parameters
  InitParameters();
  
  for(Int_t i = 0; i < 4; i++)
  {
    fhPtPhotonLeading[i]          = fhPtPi0Leading[i]          = 0;
    fhPtPhotonLeadingIsolated[i]  = fhPtPi0LeadingIsolated[i]  = 0;
    for(Int_t j = 0; j < 2; j++)
    {
    fhZHardPhoton[j][i]           = fhZHardPi0[j][i]           = 0;            
    fhZHardPhotonIsolated[j][i]   = fhZHardPi0Isolated[j][i]   = 0; 
    fhZPartonPhoton[j][i]         = fhZPartonPi0[j][i]         = 0;            
    fhZPartonPhotonIsolated[j][i] = fhZPartonPi0Isolated[j][i] = 0; 
    fhZJetPhoton[j][i]            = fhZJetPi0[j][i]            = 0;            
    fhZJetPhotonIsolated[j][i]    = fhZJetPi0Isolated[j][i]    = 0; 
    fhXEPhoton[j][i]              = fhXEPi0[j][i]              = 0;            
    fhXEPhotonIsolated[j][i]      = fhXEPi0Isolated[j][i]      = 0; 
    fhXEUEPhoton[j][i]            = fhXEUEPi0[j][i]            = 0;            
    fhXEUEPhotonIsolated[j][i]    = fhXEUEPi0Isolated[j][i]    = 0; 

    fhPtPartonTypeNearPhoton[j][i]         = fhPtPartonTypeNearPi0[j][i]         = 0;            
    fhPtPartonTypeNearPhotonIsolated[j][i] = fhPtPartonTypeNearPi0Isolated[j][i] = 0; 

    fhPtPartonTypeAwayPhoton[j][i]         = fhPtPartonTypeAwayPi0[j][i]         = 0;            
    fhPtPartonTypeAwayPhotonIsolated[j][i] = fhPtPartonTypeAwayPi0Isolated[j][i] = 0; 
    }
  }
  
}

//___________________________________________________________________________
Bool_t  AliAnaGeneratorKine::CorrelateWithPartonOrJet(TLorentzVector trigger,
                                                      Int_t   indexTrig,
                                                      Int_t   pdgTrig,
                                                      Bool_t  leading[4],
                                                      Bool_t  isolated[4],
                                                      Int_t & iparton )  
{
  //Correlate trigger with partons or jets, get z
  
  //Get the index of the mother
  iparton =  (fStack->Particle(indexTrig))->GetFirstMother();
  TParticle * mother = fStack->Particle(iparton);
  while (iparton > 7) 
  {
    iparton   = mother->GetFirstMother();
    if(iparton < 0) { printf("AliAnaGeneratorKine::CorrelateWithPartonOrJet() - Negative index, skip event\n"); return kFALSE; }
    mother = fStack->Particle(iparton);
  }
  
  //printf("Mother is parton %d with pdg %d\n",imom,mother->GetPdgCode());
  
  if(iparton < 6)
  {
    //printf("This particle is not from hard process - pdg %d, parton index %d\n",pdgTrig, iparton);
    return kFALSE; 
  }
  
  Float_t ptTrig   = trigger.Pt(); 
  Float_t partonPt = fParton6->Pt();
  Float_t jetPt    = fJet6.Pt();
  if(iparton==7)
  {
    partonPt = fParton6->Pt();
    jetPt    = fJet6.Pt();
  }
  
  //Get id of parton in near and away side
  
  Int_t away = -1;
  Int_t near = -1;
  Int_t nearPDG = -1;
  Int_t awayPDG = -1;
  
  //printf("parton 6 pdg = %d, parton 7 pdg = %d\n",fParton6->GetPdgCode(),fParton7->GetPdgCode());
  
  if(iparton==6)
  {
    nearPDG = fParton6->GetPdgCode();
    awayPDG = fParton7->GetPdgCode();
  }
  else 
  {
    nearPDG = fParton7->GetPdgCode();
    awayPDG = fParton6->GetPdgCode();
  }

  if     (nearPDG == 22) near = 0;
  else if(nearPDG == 21) near = 1;
  else                   near = 2;
  
  if     (awayPDG == 22) away = 0;
  else if(awayPDG == 21) away = 1;
  else                   away = 2;
  
  for( Int_t i = 0; i < 4; i++ )
  {
    if(pdgTrig==111)
    {
      
      fhPtPartonTypeNearPi0[leading[i]][i]->Fill(ptTrig,near);
      fhPtPartonTypeAwayPi0[leading[i]][i]->Fill(ptTrig,away);
      if(isolated[i])
      {
        fhPtPartonTypeNearPi0Isolated[leading[i]][i]->Fill(ptTrig,near);
        fhPtPartonTypeAwayPi0Isolated[leading[i]][i]->Fill(ptTrig,away);
      }
      
    }// pi0
    else if(pdgTrig==22)
    {
      fhPtPartonTypeNearPhoton[leading[i]][i]->Fill(ptTrig,near);
      fhPtPartonTypeAwayPhoton[leading[i]][i]->Fill(ptTrig,away);
      if(isolated[i])
      {
        fhPtPartonTypeNearPhotonIsolated[leading[i]][i]->Fill(ptTrig,near);
        fhPtPartonTypeAwayPhotonIsolated[leading[i]][i]->Fill(ptTrig,away);
      }
      
    } // photon
  } // conditions loop
  
  
  // RATIOS
  
  fhPtPartonPtHard->Fill(fPtHard, partonPt/fPtHard);
  fhPtJetPtHard   ->Fill(fPtHard, jetPt/fPtHard);
  fhPtJetPtParton ->Fill(fPtHard, jetPt/partonPt);

  Float_t zHard = ptTrig / fPtHard;
  Float_t zPart = ptTrig / partonPt;
  Float_t zJet  = ptTrig / jetPt;

  //if(zHard > 1 ) printf("*** Particle energy larger than pT hard z=%f\n",zHard); 
  
  //printf("Z : hard %2.2f, parton %2.2f, jet %2.2f\n",zHard,zPart,zJet);
  
  for( Int_t i = 0; i < 4; i++ )
  {
    if(pdgTrig==111)
    {
      fhZHardPi0[leading[i]][i]  ->Fill(ptTrig,zHard);
      fhZPartonPi0[leading[i]][i]->Fill(ptTrig,zPart);
      fhZJetPi0[leading[i]][i]   ->Fill(ptTrig,zJet );
      
      if(isolated[i])
      {
        fhZHardPi0Isolated[leading[i]][i]  ->Fill(ptTrig,zHard);
        fhZPartonPi0Isolated[leading[i]][i]->Fill(ptTrig,zPart);
        fhZJetPi0Isolated[leading[i]][i]   ->Fill(ptTrig,zJet);
      }
      
    }// pi0
    else if(pdgTrig==22)
    {
      
      fhZHardPhoton[leading[i]][i]  ->Fill(ptTrig,zHard);
      fhZPartonPhoton[leading[i]][i]->Fill(ptTrig,zPart);
      fhZJetPhoton[leading[i]][i]   ->Fill(ptTrig,zJet );
      
      if(isolated[i])
      {
        fhZHardPhotonIsolated[leading[i]][i]  ->Fill(ptTrig,zHard);
        fhZPartonPhotonIsolated[leading[i]][i]->Fill(ptTrig,zPart);
        fhZJetPhotonIsolated[leading[i]][i]   ->Fill(ptTrig,zJet);
      }
      
    } // photon
  } // conditions loop
  
  return kTRUE;
}


//____________________________________________________
TList *  AliAnaGeneratorKine::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file 
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("GenKineHistos") ; 
  
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins();  
  Float_t ptmax  = GetHistogramRanges()->GetHistoPtMax();  
  Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin(); 
  
  fhPtHard  = new TH1F("hPtHard"," pt hard for selected triggers",nptbins,ptmin,ptmax); 
  fhPtHard->SetXTitle("p_{T}^{hard} (GeV/c)");
  outputContainer->Add(fhPtHard);
  
  fhPtParton  = new TH1F("hPtParton"," pt parton for selected triggers",nptbins,ptmin,ptmax); 
  fhPtParton->SetXTitle("p_{T}^{parton} (GeV/c)");
  outputContainer->Add(fhPtParton);
  
  fhPtJet  = new TH1F("hPtJet"," pt jet for selected triggers",nptbins,ptmin,ptmax); 
  fhPtJet->SetXTitle("p_{T}^{jet} (GeV/c)");
  outputContainer->Add(fhPtJet);
  
  fhPtPartonPtHard  = new TH2F("hPtPartonPtHard","parton pt / pt hard for selected triggers",nptbins,ptmin,ptmax,200,0,2); 
  fhPtPartonPtHard->SetXTitle("p_{T}^{hard} (GeV/c)");
  fhPtPartonPtHard->SetYTitle("p_{T}^{parton}/p_{T}^{hard}");
  outputContainer->Add(fhPtPartonPtHard);
  
  fhPtJetPtHard  = new TH2F("hPtJetPtHard","jet pt / pt hard for selected triggers",nptbins,ptmin,ptmax,200,0,2); 
  fhPtJetPtHard->SetXTitle("p_{T}^{hard} (GeV/c)");
  fhPtJetPtHard->SetYTitle("p_{T}^{jet}/p_{T}^{hard}");
  outputContainer->Add(fhPtJetPtHard);
  
  fhPtJetPtParton  = new TH2F("hPtJetPtParton","parton pt / pt hard for selected triggers",nptbins,ptmin,ptmax,200,0,2); 
  fhPtJetPtParton->SetXTitle("p_{T}^{hard} (GeV/c)");
  fhPtJetPtParton->SetYTitle("p_{T}^{jet}/p_{T}^{parton}");
  outputContainer->Add(fhPtJetPtParton);
  
  
  fhPtPhoton  = new TH1F("hPtPhoton","Input Photon",nptbins,ptmin,ptmax); 
  fhPtPhoton->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtPhoton);

  fhPtPi0  = new TH1F("hPtPi0","Input Pi0",nptbins,ptmin,ptmax); 
  fhPtPi0->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtPi0);
  
  TString name   [] = {"","_EMC","_Photon","_EMC_Photon"};
  TString title  [] = {"",", neutral in EMCal",", neutral only photon like",", neutral in EMCal and only photon like"};
  TString leading[] = {"NotLeading","Leading"};
  
  for(Int_t i = 0; i < 4; i++)
  {
    
    // Pt
    
    fhPtPhotonLeading[i]  = new TH1F(Form("hPtPhotonLeading%s",name[i].Data()),
                                     Form("Photon : Leading of all particles%s",title[i].Data()),
                                     nptbins,ptmin,ptmax);
    fhPtPhotonLeading[i]->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtPhotonLeading[i]);
    
    fhPtPi0Leading[i]  = new TH1F(Form("hPtPi0Leading%s",name[i].Data()),
                                  Form("Pi0 : Leading of all particles%s",title[i].Data()),
                                  nptbins,ptmin,ptmax);
    fhPtPi0Leading[i]->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtPi0Leading[i]);
    
    fhPtPhotonLeadingIsolated[i]  = new TH1F(Form("hPtPhotonLeadingIsolated%s",name[i].Data()),
                                             Form("Photon : Leading of all particles%s, isolated",title[i].Data()),
                                             nptbins,ptmin,ptmax);
    fhPtPhotonLeadingIsolated[i]->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtPhotonLeadingIsolated[i]);
    
    fhPtPi0LeadingIsolated[i]  = new TH1F(Form("hPtPi0LeadingIsolated%s",name[i].Data()),
                                          Form("Pi0 : Leading of all particles%s, isolated",title[i].Data()),
                                          nptbins,ptmin,ptmax);
    fhPtPi0LeadingIsolated[i]->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtPi0LeadingIsolated[i]);
    
    // Leading or not loop
    for(Int_t j = 0; j < 2; j++)
    {
      // Near side parton
      
      fhPtPartonTypeNearPhoton[j][i]  = new TH2F(Form("hPtPartonTypeNearPhoton%s%s",leading[j].Data(),name[i].Data()),
                                                 Form("Photon : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                                 nptbins,ptmin,ptmax,3,0,3);
      fhPtPartonTypeNearPhoton[j][i]->SetXTitle("p_{T} (GeV/c)");
      fhPtPartonTypeNearPhoton[j][i]->SetYTitle("Parton type");
      fhPtPartonTypeNearPhoton[j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
      fhPtPartonTypeNearPhoton[j][i]->GetYaxis()->SetBinLabel(2,"g");
      fhPtPartonTypeNearPhoton[j][i]->GetYaxis()->SetBinLabel(3,"q");
      outputContainer->Add(fhPtPartonTypeNearPhoton[j][i]);
      
      fhPtPartonTypeNearPi0[j][i]  = new TH2F(Form("hPtPartonTypeNearPi0%s%s",leading[j].Data(),name[i].Data()),
                                              Form("Pi0 : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                              nptbins,ptmin,ptmax,3,0,3);
      fhPtPartonTypeNearPi0[j][i]->SetXTitle("p_{T} (GeV/c)");
      fhPtPartonTypeNearPi0[j][i]->SetYTitle("Parton type");
      fhPtPartonTypeNearPi0[j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
      fhPtPartonTypeNearPi0[j][i]->GetYaxis()->SetBinLabel(2,"g");
      fhPtPartonTypeNearPi0[j][i]->GetYaxis()->SetBinLabel(3,"q");
      outputContainer->Add(fhPtPartonTypeNearPi0[j][i]);
      
      fhPtPartonTypeNearPhotonIsolated[j][i]  = new TH2F(Form("hPtPartonTypeNearPhoton%sIsolated%s",leading[j].Data(),name[i].Data()),
                                                         Form("Photon : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                                         nptbins,ptmin,ptmax,3,0,3);
      fhPtPartonTypeNearPhotonIsolated[j][i]->SetXTitle("p_{T} (GeV/c)");
      fhPtPartonTypeNearPhotonIsolated[j][i]->SetYTitle("Parton type");
      fhPtPartonTypeNearPhotonIsolated[j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
      fhPtPartonTypeNearPhotonIsolated[j][i]->GetYaxis()->SetBinLabel(2,"g");
      fhPtPartonTypeNearPhotonIsolated[j][i]->GetYaxis()->SetBinLabel(3,"q");
      outputContainer->Add(fhPtPartonTypeNearPhotonIsolated[j][i]);
      
      fhPtPartonTypeNearPi0Isolated[j][i]  = new TH2F(Form("hPtPartonTypeNearPi0%sIsolated%s",leading[j].Data(),name[i].Data()),
                                                      Form("Pi0 : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                                      nptbins,ptmin,ptmax,3,0,3);
      fhPtPartonTypeNearPi0Isolated[j][i]->SetXTitle("p_{T} (GeV/c)");
      fhPtPartonTypeNearPi0Isolated[j][i]->SetYTitle("Parton type");
      fhPtPartonTypeNearPi0Isolated[j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
      fhPtPartonTypeNearPi0Isolated[j][i]->GetYaxis()->SetBinLabel(2,"g");
      fhPtPartonTypeNearPi0Isolated[j][i]->GetYaxis()->SetBinLabel(3,"q");
      outputContainer->Add(fhPtPartonTypeNearPi0Isolated[j][i]);
      
      
      // Away side parton
      
      fhPtPartonTypeAwayPhoton[j][i]  = new TH2F(Form("hPtPartonTypeAwayPhoton%s%s",leading[j].Data(),name[i].Data()),
                                                 Form("Photon : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                                 nptbins,ptmin,ptmax,3,0,3);
      fhPtPartonTypeAwayPhoton[j][i]->SetXTitle("p_{T} (GeV/c)");
      fhPtPartonTypeAwayPhoton[j][i]->SetYTitle("Parton type");
      fhPtPartonTypeAwayPhoton[j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
      fhPtPartonTypeAwayPhoton[j][i]->GetYaxis()->SetBinLabel(2,"g");
      fhPtPartonTypeAwayPhoton[j][i]->GetYaxis()->SetBinLabel(3,"q");
      outputContainer->Add(fhPtPartonTypeAwayPhoton[j][i]);
      
      fhPtPartonTypeAwayPi0[j][i]  = new TH2F(Form("hPtPartonTypeAwayPi0%s%s",leading[j].Data(),name[i].Data()),
                                              Form("Pi0 : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                              nptbins,ptmin,ptmax,3,0,3);
      fhPtPartonTypeAwayPi0[j][i]->SetXTitle("p_{T} (GeV/c)");
      fhPtPartonTypeAwayPi0[j][i]->SetYTitle("Parton type");
      fhPtPartonTypeAwayPi0[j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
      fhPtPartonTypeAwayPi0[j][i]->GetYaxis()->SetBinLabel(2,"g");
      fhPtPartonTypeAwayPi0[j][i]->GetYaxis()->SetBinLabel(3,"q");
      outputContainer->Add(fhPtPartonTypeAwayPi0[j][i]);
      
      fhPtPartonTypeAwayPhotonIsolated[j][i]  = new TH2F(Form("hPtPartonTypeAwayPhoton%sIsolated%s",leading[j].Data(),name[i].Data()),
                                                         Form("Photon : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                                         nptbins,ptmin,ptmax,3,0,3);
      fhPtPartonTypeAwayPhotonIsolated[j][i]->SetXTitle("p_{T} (GeV/c)");
      fhPtPartonTypeAwayPhotonIsolated[j][i]->SetYTitle("Parton type");
      fhPtPartonTypeAwayPhotonIsolated[j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
      fhPtPartonTypeAwayPhotonIsolated[j][i]->GetYaxis()->SetBinLabel(2,"g");
      fhPtPartonTypeAwayPhotonIsolated[j][i]->GetYaxis()->SetBinLabel(3,"q");
      outputContainer->Add(fhPtPartonTypeAwayPhotonIsolated[j][i]);
      
      fhPtPartonTypeAwayPi0Isolated[j][i]  = new TH2F(Form("hPtPartonTypeAwayPi0%sIsolated%s",leading[j].Data(),name[i].Data()),
                                                      Form("Pi0 : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                                      nptbins,ptmin,ptmax,3,0,3);
      fhPtPartonTypeAwayPi0Isolated[j][i]->SetXTitle("p_{T} (GeV/c)");
      fhPtPartonTypeAwayPi0Isolated[j][i]->SetYTitle("Parton type");
      fhPtPartonTypeAwayPi0Isolated[j][i]->GetYaxis()->SetBinLabel(1,"#gamma");
      fhPtPartonTypeAwayPi0Isolated[j][i]->GetYaxis()->SetBinLabel(2,"g");
      fhPtPartonTypeAwayPi0Isolated[j][i]->GetYaxis()->SetBinLabel(3,"q");
      outputContainer->Add(fhPtPartonTypeAwayPi0Isolated[j][i]);
      
      // zHard
      
      fhZHardPhoton[j][i]  = new TH2F(Form("hZHardPhoton%s%s",leading[j].Data(),name[i].Data()),
                                      Form("Z-Hard of Photon : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                      nptbins,ptmin,ptmax,200,0,2);
      fhZHardPhoton[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhZHardPhoton[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhZHardPhoton[j][i]);
      
      fhZHardPi0[j][i]  = new TH2F(Form("hZHardPi0%s%s",leading[j].Data(),name[i].Data()),
                                   Form("Z-Hard of Pi0 : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                   nptbins,ptmin,ptmax,200,0,2);
      fhZHardPi0[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhZHardPi0[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhZHardPi0[j][i]);
      
      fhZHardPhotonIsolated[j][i]  = new TH2F(Form("hZHardPhoton%sIsolated%s",leading[j].Data(),name[i].Data()),
                                              Form("Z-Hard of Photon : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                              nptbins,ptmin,ptmax,200,0,2);
      fhZHardPhotonIsolated[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhZHardPhotonIsolated[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhZHardPhotonIsolated[j][i]);
      
      fhZHardPi0Isolated[j][i]  = new TH2F(Form("hZHardPi0%sIsolated%s",leading[j].Data(),name[i].Data()),
                                           Form("Z-Hard of Pi0 : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                           nptbins,ptmin,ptmax,200,0,2);
      fhZHardPi0Isolated[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhZHardPi0Isolated[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhZHardPi0Isolated[j][i]);
      
      // zHard
      
      fhZPartonPhoton[j][i]  = new TH2F(Form("hZPartonPhoton%s%s",leading[j].Data(),name[i].Data()),
                                        Form("Z-Parton of Photon : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                        nptbins,ptmin,ptmax,200,0,2);
      fhZPartonPhoton[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhZPartonPhoton[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhZPartonPhoton[j][i]);
      
      fhZPartonPi0[j][i]  = new TH2F(Form("hZPartonPi0%s%s",leading[j].Data(),name[i].Data()),
                                     Form("Z-Parton of Pi0 : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                     nptbins,ptmin,ptmax,200,0,2);
      fhZPartonPi0[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhZPartonPi0[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhZPartonPi0[j][i]);
      
      fhZPartonPhotonIsolated[j][i]  = new TH2F(Form("hZPartonPhoton%sIsolated%s",leading[j].Data(),name[i].Data()),
                                                Form("Z-Parton of Photon : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                                nptbins,ptmin,ptmax,200,0,2);
      fhZPartonPhotonIsolated[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhZPartonPhotonIsolated[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhZPartonPhotonIsolated[j][i]);
      
      fhZPartonPi0Isolated[j][i]  = new TH2F(Form("hZPartonPi0%sIsolated%s",leading[j].Data(),name[i].Data()),
                                             Form("Z-Parton of Pi0 : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                             nptbins,ptmin,ptmax,200,0,2);
      fhZPartonPi0Isolated[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhZPartonPi0Isolated[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhZPartonPi0Isolated[j][i]);
      
      
      // zJet
      
      fhZJetPhoton[j][i]  = new TH2F(Form("hZJetPhoton%s%s",leading[j].Data(),name[i].Data()),
                                     Form("Z-Jet of Photon : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                     nptbins,ptmin,ptmax,200,0,2);
      fhZJetPhoton[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhZJetPhoton[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhZJetPhoton[j][i]);
      
      fhZJetPi0[j][i]  = new TH2F(Form("hZJetPi0%s%s",leading[j].Data(),name[i].Data()),
                                  Form("Z-Jet of Pi0 : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                  nptbins,ptmin,ptmax,200,0,2);
      fhZJetPi0[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhZJetPi0[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhZJetPi0[j][i]);
      
      fhZJetPhotonIsolated[j][i]  = new TH2F(Form("hZJetPhoton%sIsolated%s",leading[j].Data(),name[i].Data()),
                                             Form("Z-Jet of Photon : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                             nptbins,ptmin,ptmax,200,0,2);
      fhZJetPhotonIsolated[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhZJetPhotonIsolated[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhZJetPhotonIsolated[j][i]);
      
      fhZJetPi0Isolated[j][i]  = new TH2F(Form("hZJetPi0%sIsolated%s",leading[j].Data(),name[i].Data()),
                                          Form("Z-Jet of Pi0 : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                          nptbins,ptmin,ptmax,200,0,2);
      fhZJetPi0Isolated[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhZJetPi0Isolated[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhZJetPi0Isolated[j][i]);
      
      
      // XE
      
      fhXEPhoton[j][i]  = new TH2F(Form("hXEPhoton%s%s",leading[j].Data(),name[i].Data()),
                                   Form("Z-Jet of Photon : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                   nptbins,ptmin,ptmax,200,0,2);
      fhXEPhoton[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhXEPhoton[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhXEPhoton[j][i]);
      
      fhXEPi0[j][i]  = new TH2F(Form("hXEPi0%s%s",leading[j].Data(),name[i].Data()),
                                Form("Z-Jet of Pi0 : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                nptbins,ptmin,ptmax,200,0,2);
      fhXEPi0[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhXEPi0[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhXEPi0[j][i]);
      
      fhXEPhotonIsolated[j][i]  = new TH2F(Form("hXEPhoton%sIsolated%s",leading[j].Data(),name[i].Data()),
                                           Form("Z-Jet of Photon : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                           nptbins,ptmin,ptmax,200,0,2);
      fhXEPhotonIsolated[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhXEPhotonIsolated[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhXEPhotonIsolated[j][i]);
      
      fhXEPi0Isolated[j][i]  = new TH2F(Form("hXEPi0%sIsolated%s",leading[j].Data(),name[i].Data()),
                                        Form("Z-Jet of Pi0 : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                        nptbins,ptmin,ptmax,200,0,2);
      fhXEPi0Isolated[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhXEPi0Isolated[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhXEPi0Isolated[j][i]);
      
      
      // XE from UE
      
      fhXEUEPhoton[j][i]  = new TH2F(Form("hXEUEPhoton%s%s",leading[j].Data(),name[i].Data()),
                                     Form("Z-Jet of Photon : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                     nptbins,ptmin,ptmax,200,0,2);
      fhXEUEPhoton[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhXEUEPhoton[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhXEUEPhoton[j][i]);
      
      fhXEUEPi0[j][i]  = new TH2F(Form("hXEUEPi0%s%s",leading[j].Data(),name[i].Data()),
                                  Form("Z-Jet of Pi0 : %s of all particles%s",leading[j].Data(),title[i].Data()),
                                  nptbins,ptmin,ptmax,200,0,2);
      fhXEUEPi0[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhXEUEPi0[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhXEUEPi0[j][i]);
      
      fhXEUEPhotonIsolated[j][i]  = new TH2F(Form("hXEUEPhoton%sIsolated%s",leading[j].Data(),name[i].Data()),
                                             Form("Z-Jet of Photon : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                             nptbins,ptmin,ptmax,200,0,2);
      fhXEUEPhotonIsolated[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhXEUEPhotonIsolated[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhXEUEPhotonIsolated[j][i]);
      
      fhXEUEPi0Isolated[j][i]  = new TH2F(Form("hXEUEPi0%sIsolated%s",leading[j].Data(),name[i].Data()),
                                          Form("Z-Jet of Pi0 : %s of all particles%s, isolated",leading[j].Data(),title[i].Data()),
                                          nptbins,ptmin,ptmax,200,0,2); 
      fhXEUEPi0Isolated[j][i]->SetYTitle("p_{T}^{particle}/p_{T}^{hard}");
      fhXEUEPi0Isolated[j][i]->SetXTitle("p_{T}^{particle} (GeV/c)");
      outputContainer->Add(fhXEUEPi0Isolated[j][i]);          
    }
  }
  
  return outputContainer;
  
}

//____________________________________________
void  AliAnaGeneratorKine::GetPartonsAndJets() 
{
  // Fill data members with partons,jets and generated pt hard 
  
  fStack =  GetMCStack() ;
  
  if(!fStack) 
  {
    printf("AliAnaGeneratorKine::MakeAnalysisFillHistograms() - No Stack available, STOP\n");
    abort();
  }
  
  fParton2 = fStack->Particle(2) ;
  fParton3 = fStack->Particle(3) ;
  fParton6 = fStack->Particle(6) ;
  fParton7 = fStack->Particle(7) ;
  
  Float_t p6phi = fParton6->Phi();
  if(p6phi < 0) p6phi +=TMath::TwoPi();
  Float_t p7phi = fParton7->Phi();
  if(p7phi < 0) p7phi +=TMath::TwoPi();  
  
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
      
      TLorentzVector jet(tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3]);
      Float_t jphi = jet.Phi();
      if(jphi < 0) jphi +=TMath::TwoPi();
      
      Double_t radius6 = GetIsolationCut()->Radius(fParton6->Eta(), p6phi, jet.Eta() , jphi) ;
      Double_t radius7 = GetIsolationCut()->Radius(fParton7->Eta(), p7phi, jet.Eta() , jphi) ;
      
      //printf("jet %d: pt %2.2f, eta %2.2f, phi %2.2f, r6 %2.2f, r7 %2.2f\n",ijet,jet.Pt(),jet.Eta(),jphi,radius6, radius7);
      
      if (radius6 < jet6R)
      {
        jet6R = radius6;
        fJet6 = jet;
        
      }
      if (radius7 < jet7R) 
      {
        jet7R = radius7;
        fJet7 = jet;
      }
            
    } // jet loop
    
    //printf("jet6: pt %2.2f, eta %2.2f, phi %2.2f\n",fJet6.Pt(),fJet6.Eta(),fJet6.Phi());
    //printf("jet7: pt %2.2f, eta %2.2f, phi %2.2f\n",fJet7.Pt(),fJet7.Eta(),fJet6.Phi());
    
  } // pythia header
  
  fhPtHard   ->Fill(fPtHard);
  fhPtJet    ->Fill(fJet6.Pt());
  fhPtJet    ->Fill(fJet7.Pt());
  fhPtParton ->Fill(fParton6->Pt());
  fhPtParton ->Fill(fParton7->Pt());

}

//_____________________________________________________
void AliAnaGeneratorKine::GetXE(TLorentzVector trigger,
                                Int_t   indexTrig,
                                Int_t   pdgTrig,
                                Bool_t  leading[4],
                                Bool_t  isolated[4],
                                Int_t   iparton)
{

  // Calculate the real XE and the UE XE

  Float_t ptThresTrack = 0.2;

  Float_t ptTrig  = trigger.Pt();
  Float_t etaTrig = trigger.Eta();
  Float_t phiTrig = trigger.Phi();
  if(phiTrig < 0 ) phiTrig += TMath::TwoPi();
  
  //Loop on primaries, start from position 8, no partons
  for(Int_t ipr = 8; ipr < fStack->GetNprimary(); ipr ++ )
  {
    TParticle * particle = fStack->Particle(ipr) ;
    
    if(ipr==indexTrig) continue;
    
    
    Int_t   pdg    = particle->GetPdgCode();
    Int_t   status = particle->GetStatusCode();
        
    // Compare trigger with final state particles
    if( status != 1) continue ;
    
    Double_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
    
    if(charge==0) continue; // construct xe only with charged        
    
    Float_t pt     = particle->Pt();
    Float_t eta    = particle->Eta();
    Float_t phi    = particle->Phi();
    if(phi < 0 ) phi += TMath::TwoPi();
    
    if( pt < ptThresTrack)    continue ;
    
    if(TMath::Abs(eta) > 0.8) continue ; // TPC acceptance cut

    //Isolation
    Double_t radius = GetIsolationCut()->Radius(etaTrig, phiTrig, eta , phi) ;
    
    Float_t xe = -pt/ptTrig*TMath::Cos(phi-phiTrig);
    
    //Get the index of the mother
    Int_t ipartonAway =  particle->GetFirstMother();
    if(ipartonAway < 0) return;
    TParticle * mother = fStack->Particle(ipartonAway);
    while (ipartonAway > 7) 
    {
      ipartonAway   = mother->GetFirstMother();
      if(ipartonAway < 0) break;
      mother = fStack->Particle(ipartonAway);
    }
    
    if((ipartonAway==6 || ipartonAway==7) && iparton!=ipartonAway) 
    {
      //printf("xE : iparton %d, ipartonAway %d\n",iparton,ipartonAway);
      if(radius > 1 ) continue; // avoid particles too far from trigger
      
      for( Int_t i = 0; i < 4; i++ )
      {
        if(pdgTrig==111)
        {
          
          fhXEPi0[leading[i]][i]  ->Fill(ptTrig,xe);
          
          if(isolated[i])
          {
            fhXEPi0Isolated[leading[i]][i]  ->Fill(ptTrig,xe);
          }
          
        }// pi0
        else if(pdgTrig==22)
        {
          
          fhXEPhoton[leading[i]][i]  ->Fill(ptTrig,xe);
          
          if(isolated[i])
          {
            fhXEPhotonIsolated[leading[i]][i]  ->Fill(ptTrig,xe);
          }
          
        } // photon
      } // conditions loop
    } // Away side

    if(ipartonAway!=6 && ipartonAway!=7)
    {
      
      //printf("xE UE : iparton %d, ipartonAway %d\n",iparton,ipartonAway);
      
      for( Int_t i = 0; i < 4; i++ )
      {
        if(pdgTrig==111)
        {
          
          fhXEUEPi0[leading[i]][i]  ->Fill(ptTrig,xe);
          
          if(isolated[i])
          {
            fhXEUEPi0Isolated[leading[i]][i]  ->Fill(ptTrig,xe);
          }
          
        }// pi0
        else if(pdgTrig==22)
        {
          
          fhXEUEPhoton[leading[i]][i]  ->Fill(ptTrig,xe);
          
          if(isolated[i])
          {
            fhXEUEPhotonIsolated[leading[i]][i]  ->Fill(ptTrig,xe);
          }
          
        } // photon
      } // conditions loop  
    } // Away side    
    
  } // primary loop


}

//________________________________________
void AliAnaGeneratorKine::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaGenKine_");
  
}

//_____________________________________________________________________
void  AliAnaGeneratorKine::IsLeadingAndIsolated(TLorentzVector trigger,
                                                Int_t indexTrig,
                                                Int_t pdgTrig,
                                                Bool_t leading[4],
                                                Bool_t isolated[4]) 
{
  // Check if the trigger is the leading particle and if it is isolated
  // In case of neutral particles check all neutral or neutral in EMCAL acceptance
  
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
  
  Float_t ptTrig  = trigger.Pt();
  Float_t etaTrig = trigger.Eta();
  Float_t phiTrig = trigger.Phi();
  if(phiTrig < 0 ) phiTrig += TMath::TwoPi();

  // Minimum track or cluster energy
  Float_t ptThresTrack = 0.2;
  Float_t ptThresCalo  = 0.3;

  //Isolation cuts
  Float_t ptThresIC    = 0.5;
  Float_t rThresIC     = 0.4;
  Int_t   nICTrack     = 0;
  Int_t   nICNeutral   = 0;
  Int_t   nICNeutEMCAL = 0;
  Int_t   nICNeutPhot  = 0;
  Int_t   nICNeutEMCALPhot = 0;
  
  //Loop on primaries, start from position 8, no partons
  for(Int_t ipr = 8; ipr < fStack->GetNprimary(); ipr ++ )
  {
    if(ipr == indexTrig) continue;
    TParticle * particle = fStack->Particle(ipr) ;
    
    Int_t   imother= particle->GetFirstMother();
    //printf("Leading ipr %d - mother %d\n",ipr, imother);
    
    if(imother==indexTrig)  continue ;
    
    Int_t   pdg    = particle->GetPdgCode();
    Int_t   status = particle->GetStatusCode();
     
    // Compare trigger with final state particles
    if( status != 1) continue ;
    
    Float_t pt     = particle->Pt();
    Float_t eta    = particle->Eta();
    Float_t phi    = particle->Phi();
    if(phi < 0 ) phi += TMath::TwoPi();
    
    if(TMath::Abs(eta) > 0.8) continue ; // TPC acceptance cut

    Double_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
    
    //Isolation
    Double_t radius = GetIsolationCut()->Radius(etaTrig, phiTrig, eta , phi) ;
    
    if(charge==0)
    {
      if(pt < ptThresCalo)  continue ;
      
      if( ptMaxNeutral < pt )                   ptMaxNeutral = pt;
      if( pt > ptThresIC && radius < rThresIC ) nICNeutral++ ;

      Bool_t phPDG = kFALSE;
      if(pdg==22 || pdg==111) phPDG = kTRUE;
    
      //if(pt > ptTrig) printf(" --- pdg %d, phPDG %d pT %2.2f, pTtrig %2.2f, eta %2.2f, phi %2.2f\n",pdg,phPDG,pt,ptTrig,particle->Eta(), particle->Phi()*TMath::RadToDeg());
      if(phPDG)
      {
        if( ptMaxNeutPhot < pt)                   ptMaxNeutPhot = pt;
        if( pt > ptThresIC && radius < rThresIC ) nICNeutPhot++ ;
      }
      
      //EMCAL acceptance
      Bool_t inEMCAL = GetFiducialCut()->IsInFiducialCut(trigger,"EMCAL") ;
      
      if(inEMCAL)
      {
        if( ptMaxNeutEMCAL < pt )                 ptMaxNeutEMCAL = pt;
        if( pt > ptThresIC && radius < rThresIC ) nICNeutEMCAL++ ;

        if(phPDG)
        {
          if( ptMaxNeutEMCALPhot < pt )             ptMaxNeutEMCALPhot = pt;
          if( pt > ptThresIC && radius < rThresIC ) nICNeutEMCALPhot++ ;
        }
      }
    }
    else  
    {
      if( pt < ptThresTrack)  continue ;

      if( ptMaxCharged < pt )   ptMaxCharged   = pt;
      
      if( pt > ptThresIC && radius < rThresIC )  
      {
        //printf("UE track? pTtrig %2.2f, pt %2.2f, etaTrig %2.2f,  eta %2.2f, phiTrig %2.2f,  phi %2.2f, radius %2.2f\n",
        //       ptTrig, pt,etaTrig, eta, phiTrig*TMath::RadToDeg(), phi*TMath::RadToDeg(),radius);
        nICTrack++ ;
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
  
  //printf("N in cone over threshold : tracks  %d, neutral %d, neutral emcal %d, photon %d, photon emcal %d\n", 
  //       nICTrack, nICNeutral ,nICNeutEMCAL,nICNeutPhot, nICNeutEMCALPhot);
  
  // Isolation decision
  if(nICTrack == 0)
  {
    if(nICNeutral       == 0 ) isolated[0] = kTRUE ;
    if(nICNeutEMCAL     == 0 ) isolated[1] = kTRUE ;
    if(nICNeutPhot      == 0 ) isolated[2] = kTRUE ;
    if(nICNeutEMCALPhot == 0 ) isolated[3] = kTRUE ;
  }
  
  // Fill histograms if conditions apply for all 4 cases
  for( Int_t i = 0; i < 4; i++ )
  {
    if(pdgTrig==111)
    {
      if(leading[i])
      { 
        fhPtPi0Leading[i]->Fill(ptTrig);
        if(isolated[i]) fhPtPi0LeadingIsolated[i]->Fill(ptTrig);
      }
    }// pi0
    else if(pdgTrig==22)
    {
      if(leading[i])
      { 
        fhPtPhotonLeading[i]->Fill(ptTrig);
        if(isolated[i]) fhPtPhotonLeadingIsolated[i]->Fill(ptTrig);
      }      
    } // photon
  } // conditions loop
    
}
  
//_____________________________________________________
void  AliAnaGeneratorKine::MakeAnalysisFillHistograms() 
{
  //Particle-Parton Correlation Analysis, fill histograms
  
  TLorentzVector trigger;
  
  GetPartonsAndJets();
  
  for(Int_t ipr = 0; ipr < fStack->GetNprimary(); ipr ++ )
  {
    TParticle * particle = fStack->Particle(ipr) ;
    
    Int_t   pdgTrig    = particle->GetPdgCode();
    Int_t   statusTrig = particle->GetStatusCode();
    Int_t   imother    = particle->GetFirstMother();
    Float_t ptTrig     = particle->Pt(); 

    // Select final state photons (prompt, fragmentation) or pi0s
    
    //Check the origin of the photon, accept if prompt or initial/final state radiation
    if(pdgTrig == 22 && statusTrig == 1)
    {
      if(imother > 8) continue;
    }
    // If not photon, trigger on pi0
    else if(pdgTrig != 111) continue;
    
    // Acceptance and kinematical cuts
    if( ptTrig < GetMinPt() )    continue ;
    
    //EMCAL acceptance, a bit less
    if(TMath::Abs(particle->Eta()) > 0.6) continue ;
    if(particle->Phi() > TMath::DegToRad()*176) continue ;
    if(particle->Phi() < TMath::DegToRad()*74 ) continue ;
    
    particle->Momentum(trigger);

//    Bool_t in = GetFiducialCut()->IsInFiducialCu(trigger,"EMCAL") ;
//    if(! in ) continue ;
    
//    printf("Particle %d : pdg %d status %d, mother index %d, pT %2.2f, eta %2.2f, phi %2.2f \n",
//           ipr, pdgTrig, statusTrig, imother, ptTrig, particle->Eta(), particle->Phi()*TMath::RadToDeg());
    
//    if(pdgTrig==111)
//    {
//      printf("\t pi0 daughters %d, %d\n", particle->GetDaughter(0), particle->GetDaughter(1));
//    }

    if     (pdgTrig==22 ) fhPtPhoton->Fill(ptTrig);
    else if(pdgTrig==111) fhPtPi0   ->Fill(ptTrig);
    
    // Check if it is leading
    Bool_t leading[4] ;
    Bool_t isolated[4] ;

    IsLeadingAndIsolated(trigger, ipr, pdgTrig, leading, isolated);
    
    Int_t iparton = -1;
    Int_t ok = CorrelateWithPartonOrJet(trigger, ipr, pdgTrig, leading, isolated, iparton); 
    if(!ok) continue;
    
    GetXE(trigger,ipr,pdgTrig,leading,isolated,iparton) ;    
    
  }
  
  if(GetDebug() > 1) printf("AliAnaGeneratorKine::MakeAnalysisFillHistograms() - End fill histograms \n");
  
} 
