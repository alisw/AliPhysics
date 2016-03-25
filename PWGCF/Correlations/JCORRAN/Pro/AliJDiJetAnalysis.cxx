/*************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

// Comment describing what this class does needed!

//===========================================================
// Dummy comment, should be replaced by a real one
// comment
// comment
//===========================================================

#include <TRandom.h>
#include <TMath.h>
#include "AliJJet.h"
#include "AliJDiJetAnalysis.h"
#include "AliJHistManager.h"

AliJDiJetAnalysis::AliJDiJetAnalysis():
  fInputList(NULL),
  fJetList(NULL),
  fJetListOfList(),
  fTrackOrMCParticle(),
  fChargedOrFull(),
  fIsFullAcceptance(0),
  fJetEtaRange(0),
  fIsMartasDiJet(0),
  fTrackJetMap(NULL),
  fJetPtBins(NULL),
  fJetPtPairBins(NULL),
  fJetPtMinCut(0),
  fMeanDiJetPt(0),
  fParton23Pt(0),
  fCard(NULL),
  fHMG(NULL),
  fJetPtBin(),
  fJetPtPairBin(),
  fInvMBin(),
  fJetTypeBin(),
  fJetSelectionBin(),
  fDiJetSelectionBin(),
  fJetRBin(),
  fJetDRBin(),
  fDiJetBinTypeBin(),
  fPYJetTypeBin(),
  fBin2(),
  fhJetPt (),
  fhJetPhi(),
  fhJetEta(),
  fhJetNConst(),
  fhJetEnergy(),
  fhJetEnergyComp(),
  fhJetPtComp(),
  fhJetDRToRef(),
  fhJetDPhiToRef(),
  fhDiJetPtPair(),
  fhDiJetInvM(),
  fhDiJetKtA(),
  fhDiJetDeltaR(),
  fhDiJetDeltaPhi(),
  fhDiJetDeltaEta(),
  fhDiJetPtAsymm(),
  fhDiJetEAsymm(),
  fhDiJetMultiplicity(),
  fhPythiaJetPtPair(),
  fhPythiaJetSum()
{
  // Constaractor
}

AliJDiJetAnalysis::AliJDiJetAnalysis( AliJBaseCard * card ):
  fInputList(NULL),
  fJetList(NULL),
  fJetListOfList(),
  fTrackOrMCParticle(),
  fChargedOrFull(),
  fIsFullAcceptance(0),
  fJetEtaRange(0),
  fIsMartasDiJet(0),
  fTrackJetMap(NULL),
  fJetPtBins(NULL),
  fJetPtPairBins(NULL),
  fJetPtMinCut(0),
  fMeanDiJetPt(0),
  fParton23Pt(0),
  fCard(card),
  fHMG(NULL),
  fJetPtBin(),
  fJetPtPairBin(),
  fInvMBin(),
  fJetTypeBin(),
  fJetSelectionBin(),
  fDiJetSelectionBin(),
  fJetRBin(),
  fJetDRBin(),
  fDiJetBinTypeBin(),
  fPYJetTypeBin(),
  fBin2(),
  fhJetPt (),
  fhJetPhi(),
  fhJetEta(),
  fhJetNConst(),
  fhJetEnergy(),
  fhJetEnergyComp(),
  fhJetPtComp(),
  fhJetDRToRef(),
  fhJetDPhiToRef(),
  fhDiJetPtPair(),
  fhDiJetInvM(),
  fhDiJetKtA(),
  fhDiJetDeltaR(),
  fhDiJetDeltaPhi(),
  fhDiJetDeltaEta(),
  fhDiJetPtAsymm(),
  fhDiJetEAsymm(),
  fhDiJetMultiplicity(),
  fhPythiaJetPtPair(),
  fhPythiaJetSum()
{
  // Constaractor
}

AliJDiJetAnalysis::AliJDiJetAnalysis( const AliJDiJetAnalysis & obj ):
  fInputList(obj.fInputList),
  fJetList(obj.fJetList),
  fJetListOfList(obj.fJetListOfList),
  fTrackOrMCParticle(obj.fTrackOrMCParticle),
  fChargedOrFull(obj.fChargedOrFull),
  fIsFullAcceptance(obj.fIsFullAcceptance),
  fJetEtaRange(obj.fJetEtaRange),
  fIsMartasDiJet(obj.fIsMartasDiJet),
  fTrackJetMap(obj.fTrackJetMap),
  fJetPtBins(obj.fJetPtBins),
  fJetPtPairBins(obj.fJetPtPairBins),
  fJetPtMinCut(obj.fJetPtMinCut),
  fMeanDiJetPt(obj.fMeanDiJetPt),
  fParton23Pt(obj.fParton23Pt),
  fCard(obj.fCard),
  fHMG(obj.fHMG),
  fJetPtBin(obj.fJetPtBin),
  fJetPtPairBin(obj.fJetPtPairBin),
  fInvMBin(obj.fInvMBin),
  fJetTypeBin(obj.fJetTypeBin),
  fJetSelectionBin(obj.fJetSelectionBin),
  fDiJetSelectionBin(obj.fDiJetSelectionBin),
  fJetRBin(obj.fJetRBin),
  fJetDRBin(obj.fJetDRBin),
  fDiJetBinTypeBin(obj.fDiJetBinTypeBin),
  fPYJetTypeBin(obj.fPYJetTypeBin),
  fBin2(obj.fBin2),
  fhJetPt(obj.fhJetPt),
  fhJetPhi(obj.fhJetPhi),
  fhJetEta(obj.fhJetEta),
  fhJetNConst(obj.fhJetNConst),
  fhJetEnergy(obj.fhJetEnergy),
  fhJetEnergyComp(obj.fhJetEnergyComp),
  fhJetPtComp(obj.fhJetPtComp),
  fhJetDRToRef(obj.fhJetDRToRef),
  fhJetDPhiToRef(obj.fhJetDPhiToRef),
  fhDiJetPtPair(obj.fhDiJetPtPair),
  fhDiJetInvM(obj.fhDiJetInvM),
  fhDiJetKtA(obj.fhDiJetKtA),
  fhDiJetDeltaR(obj.fhDiJetDeltaR),
  fhDiJetDeltaPhi(obj.fhDiJetDeltaPhi),
  fhDiJetDeltaEta(obj.fhDiJetDeltaEta),
  fhDiJetPtAsymm(obj.fhDiJetPtAsymm),
  fhDiJetEAsymm(obj.fhDiJetEAsymm),
  fhDiJetMultiplicity(obj.fhDiJetMultiplicity),
  fhPythiaJetPtPair(obj.fhPythiaJetPtPair),
  fhPythiaJetSum(obj.fhPythiaJetSum)
{
  // Constaractor
}

AliJDiJetAnalysis& AliJDiJetAnalysis::operator=(const AliJDiJetAnalysis & obj){
  if( this != &obj ){
    this->~AliJDiJetAnalysis();
    new(this) AliJDiJetAnalysis(obj);
  }
  return *this;
}

bool AliJDiJetAnalysis::CompareTrackByPt( AliJBaseTrack * a, AliJBaseTrack* b ){
  return (a->Pt()>b->Pt());
}

vector<AliJBaseTrack*> AliJDiJetAnalysis::SortTrackByPt( TObjArray * os ){
  vector<AliJBaseTrack*> a(os->GetEntriesFast());
  for( int i=0;i<os->GetEntriesFast();i++ ){
    a[i] = (AliJBaseTrack*)os->At(i);
  }
  sort( a.begin(), a.end(), CompareTrackByPt );
  return a;
}

void AliJDiJetAnalysis::UserCreateOutputObjects(){
  TH1D::StatOverflows();
  fJetListOfList.Clear();
  fTrackOrMCParticle.clear(),
  fChargedOrFull.clear(),
  // comment needed
  fJetPtBins = fCard->GetVector("Jet:PtBins");
  fJetPtPairBins = fCard->GetVector("Jet:PtPairBins");
  fJetPtMinCut = fCard->Get("Jet:PtMinCut");
  CreateHistos();

  for( int i=0;i<kJNJetTypeMax;i++ ){
    for( int j=0;j<kJNJetSelection;j++ ){
      fJets[i][j] = NULL;
      for( int k=0;k<kJNDiJetSelection;k++ ){
        fDiJets[i][j][k]=NULL;
        fDiJetBin[i][j][k][0] = -1; 
        fDiJetBin[i][j][k][1] = -1; 
      }
    }
  }
}

void AliJDiJetAnalysis::ClearBeforeEvent(){
  //fJetListOfList.Clear();
  for( int i=0;i<fJetListOfList.GetEntriesFast();i++ ){
    for( int j=0;j<kJNJetSelection;j++ ){
      if( fJets[i][j] ){
        delete fJets[i][j];
        fJets[i][j] = NULL;
      }
      for( int k=0;k<kJNDiJetSelection;k++ ){
        if( fDiJets[i][j][k] ){
          delete fDiJets[i][j][k];
          fDiJets[i][j][k] = NULL;
        }
      }
    }
  }
}

void AliJDiJetAnalysis::UserExec(){
  //cout<<"DEBUG_B1 "<<fJetListOfList.GetEntries()<<endl;
  ClearBeforeEvent();
  for( int i=0;i<fJetListOfList.GetEntriesFast();i++ ){
    TObjArray * ar = (TObjArray*) fJetListOfList[i];
    if(!ar) {
      //cout<<"DEBUG_B2 no array at "<<i<<endl;
      continue;
    }
    for( int j=0;j<ar->GetEntriesFast();j++ ){
      AliJJet * jet = (AliJJet*) ar->At(j);
      //if( i > 1 && i!=kJMultiPartonAll && jet->LeadingParticleId() < 0 ) jet->ReSum();
      jet->ReSum2(); // TODO temp
    }
  }
  FillHistosJets();
  FillHistosDiJet();
  FillPythiaDiJet(0);
}
void AliJDiJetAnalysis::Terminate() const{
  // comment needed
  fHMG->Write();
  //fHMG->WriteConfig();
}


AliJDiJet * AliJDiJetAnalysis::GetDiJets(int type, int js, int djs){
  if( fDiJets[type][js][djs] ) return fDiJets[type][js][djs];
  TObjArray * jets = GetJets( type, js );
  AliJDiJet * dijet = fDiJets[type][js][djs] = new AliJDiJet();
  if( !jets ) return dijet;

  fDiJets[type][js][djs] = dijet;
  //==== DIJET DiJet - Leading - sub leading
  if( djs == kJDiJetLeadingSubLeading5 || djs==kJDiJetLeadingSubLeading10 || djs==kJDiJetLeadingSubLeading20 || djs==kJDiJetLeadingSubLeading50 ){
    for( int ij=0;ij<jets->GetEntriesFast();ij++ ){
      AliJJet * jet = (AliJJet*) jets->At(ij);
      dijet->Add( jet );
      if( dijet->GetEntries() >= 2 ) break;
    }
    double jetptcut= 0;
    if( djs == kJDiJetLeadingSubLeading5 ) jetptcut = 5;
    if( djs == kJDiJetLeadingSubLeading10 ) jetptcut = 10;
    if( djs == kJDiJetLeadingSubLeading20 ) jetptcut = 20;
    if( djs == kJDiJetLeadingSubLeading50 ) jetptcut = 50;
    if(dijet->GetEntries() > 1 && dijet->jet(0).Pt()>jetptcut && dijet->jet(1).Pt()>jetptcut){
    }else{
      dijet->Clear();
    }
  }
  //==== DIJET DiJet - leading - subleading - opposite 
  if( djs == kLeadingSubLeadingOpposite ){
    for( int ij=0;ij<jets->GetEntriesFast();ij++ ){
      AliJJet * jet = (AliJJet*) jets->At(ij);
      if( jet->Pt() < fJetPtMinCut ) break;
      dijet->Add(jet);
      if( dijet->GetEntries() >= 2 ) break;
    }
    if( dijet->GetEntries() > 1){
      double dphi = dijet->DPhi();
      if( fabs(dphi) < TMath::Pi()/2 ) dijet->Clear();
    }

  }
  //==== DIJET DiJet - leading - subleading - opposite 
  if( djs == kLeadingSubLeadingEtaWindow ){
    for( int ij=0;ij<jets->GetEntriesFast();ij++ ){
      AliJJet * jet = (AliJJet*) jets->At(ij);
      if( jet->Pt() < fJetPtMinCut ) break;
      dijet->Add(jet);
      if( dijet->GetEntries() >= 2 ) break;
    }
    if( dijet->GetEntries() < 2){

    }else{
      AliJJet & j0 = dijet->jet(0);
      AliJJet & j1 = dijet->jet(1);
      TLorentzVector lcms = (j0+j1)*0.5;
      AliJJet jj0 = j0 - lcms;
      AliJJet jj1 = j1 - lcms;
      if( fabs(jj0.Eta()) > 0.4 || fabs(jj1.Eta()) > 0.4 ){
        dijet->Clear();
      }
    }
  }
  //==== DIJET 2 : Marta - Find DiJet in Eta
  if( djs == kJDiJetMarta ){
    if(  !jets ) return dijet;
    AliJJet * trigJet = NULL;
    for( int ij=0;ij<jets->GetEntriesFast();ij++ ){
      AliJJet * jet = (AliJJet*) jets->At(ij);
      if( jet->Area() < 0.3 ) continue;
      if( jet->Pt() < 20 ) break;
      if( jet->LeadingParticlePt() < 5 || jet->LeadingParticlePt() > 100 ) continue;
      trigJet = jet;break;
    }
    if( trigJet == NULL ) return dijet;
    double ptt = trigJet->Pt();
    AliJJet * asocJet = NULL;
    for( int ij=0;ij<jets->GetEntriesFast();ij++ ){
      AliJJet * jet = (AliJJet*) jets->At(ij);
      if( !jet ){
        cout<<"JWARN_GETDIJET13: "<<dijet->GetEntries()<<endl;
        continue;
      }
      if( jet->Area() < 0.3 ) continue;
      double pta = jet->Pt();
      if( pta < 20 ) break;
      if( pta > ptt ) break;
      if( jet->LeadingParticlePt() < 5 || jet->LeadingParticlePt() > 100 ) continue;
      if( fabs(jet->DeltaPhi( *trigJet ))<TMath::Pi()/2 ) continue;
      asocJet = jet;break;
    }
    if( asocJet ){
      dijet->SetJets( trigJet,asocJet );
    }
  }
  //===== DIJET CMS 
  if( djs == kJDiJetCMS ){
    const double kRapidityCut = 0.5;
    AliJJet * trigJet = NULL;
    double trigJetPtMax = 0; 
    for( int ij=0;ij<jets->GetEntriesFast();ij++ ){
      AliJJet * jet = (AliJJet*) jets->At(ij);
      double jetPt = jet->Pt();
      if( jetPt < 60 ) continue;
      if( fabs(jet->Rapidity()) < kRapidityCut ) continue;
      if( trigJetPtMax < jetPt ){
        trigJet = jet;
        trigJetPtMax = jetPt;
      }
    }
    if( trigJet == NULL ) return dijet;
    double ptt = trigJet->Pt();
    AliJJet * asocJet = NULL;
    double asocJetPtMax = 0;
    for( int ij=0;ij<jets->GetEntriesFast();ij++ ){
      AliJJet * jet = (AliJJet*) jets->At(ij);
      if( !jet ){
        cout<<"JWARN_GETDIJET13: "<<dijet->GetEntries()<<endl;
        continue;
      }
      double pta = jet->Pt();
      if( fabs(jet->Rapidity()) < kRapidityCut ) continue;
      if( pta < 30 ) continue;
      if( pta > ptt ) break;
      if( asocJetPtMax < pta ){
        asocJet = jet;
        asocJetPtMax = pta;
      }
    }
    if( asocJet ){
      dijet->SetJets( trigJet,asocJet );
    }
  }

  return fDiJets[type][js][djs];
}

TObjArray * AliJDiJetAnalysis::GetJets(int type, int js){
  if( fJets[type][js] ) return fJets[type][js];
  TObjArray * na = new TObjArray;
  //=========================================
  //== Jet Selection
  //=========================================
  TObjArray * jetsOrg =  (TObjArray*)fJetListOfList[type];
  if( !jetsOrg ) return na;
  vector<AliJBaseTrack*> jets = SortTrackByPt(jetsOrg);
  int nEntries = jets.size();
  for( int i=0;i<nEntries;i++ ){
    AliJJet * jet = (AliJJet*)  jets[i];
    //if(  jet->LeadingParticleId() < 0 ) jet->ReSum();
     if( js == kJEtaAll ){
      //=========================================
      //== kJEtaAll
      //=========================================
      if( fabs(jet->Eta()) < 5 && jet->Pt()<5 ){
        na->Add(jet);
      }
    }else if( js == kJEtaAlice ){
      //=========================================
      //== kJEtaAlice
      //=========================================
      if( fabs(jet->Eta()) < 0.4 ){
        if( jet->Pt()<5 ) continue;
        na->Add(jet);
      }
    }
  }// jets

  // DEBUG
  if( nEntries > 2 ){
    AliJJet * jet0 = (AliJJet*)  jets[0];
    AliJJet * jet1 = (AliJJet*)  jets[1];
    AliJJet * jet2 = (AliJJet*)  jets[2];
    if( jet0->Pt() < jet1->Pt() ){
      cout<<"JERROR : jet0 is weaker than jet1 "<< jet0->Pt()<<"\t"<<jet1->Pt()<<"\t"<<jet2->Pt()<<endl;
      cout<<"JERROR : jet0 is weaker than jet1 "<< jet0->E()<<"\t"<<jet1->E()<<"\t"<<jet2->E()<<endl;
      gSystem->Exit(2);
    }
  }

  for( int i=0;i<na->GetEntriesFast();i++ ){


  }
  //=========================================
  //== NOJET
  //=========================================
  //if( na->GetEntriesFast() < 1 ){
  //    delete na;
  //    na=NULL;
  //}
  return fJets[type][js]=na;
}

void AliJDiJetAnalysis::FillHistosJets(){
  //==== TYPE 
  for( int type=0;type<fJetListOfList.GetEntriesFast();type++ ){
    if( ! fJetListOfList[type] ) continue;
    //==== Jet Selection
    for( int js=0;js<kJNJetSelection;js++ ){
      TObjArray *jets = GetJets( type, js );
      if( !jets ) continue;
      //==== jets
      for( int i=0;i<jets->GetEntriesFast();i++ ){
        AliJJet * jet = dynamic_cast<AliJJet*>( jets->At(i) );
        if( !jet ) continue;
        double jpt = jet->Pt();
        int ipt = fJetPtBin.GetBin( jpt );
        int iE = fJetPtBin.GetBin( jet->E() );
        if( ipt >= 0 ){
          fhJetInvM[type][js][ipt]->Fill( jet->M() );
        }
        fhJetPt[type][js]->Fill( jet->Pt() );
        fhJetPhi[type][js]->Fill( jet->Phi() );
        fhJetEta[type][js]->Fill( jet->Eta() );
        fhJetNConst[type][js]->Fill( jet->GetNConstituents() );
        fhJetEnergy[type][js]->Fill( jet->E() );
        if( iE >= 0 ){
          fhJetMjjEbin[type][js][iE]->Fill( jet->M() );
          fhJetMtEbin[type][js][iE]->Fill( jet->Mt() );
        }
        //== DiJet Multiplicity
        for( int j=i+1;j<jets->GetEntriesFast();j++ ){
          AliJJet * jet1 = dynamic_cast<AliJJet*>( jets->At(j) );
          AliJDiJet dijet(jet,jet1);
          fhDiJetMultiplicity[type][js]->Fill( dijet.InvM() );
          //fhDiJetInvMMultiplicityPt[type][js][ipt]->Fill( dijet.InvM() );
        }
      }
    }
  }
}

void AliJDiJetAnalysis::FillHistosDiJet(){
  //==== TYPE 
  //for( int type=0;type<fJetListOfList.GetEntriesFast();type++ ){
  for( int type=0;type<fJetListOfList.GetEntriesFast();type++ ){
    //==== Jet Selection
    for( int js=0;js<kJNJetSelection;js++ ){
      //==== DiJet Selection
      for( int djs=0;djs<kJNDiJetSelection;djs++ ){
        AliJDiJet & dijet = *GetDiJets( type, js, djs );
        fDiJetBin[type][js][djs][0] = -1;
        fDiJetBin[type][js][djs][1] = -1;
        fDiJetBin[type][js][djs][2] = -1;
        fDiJetBin[type][js][djs][3] = -1;
        fDiJetBin[type][js][djs][4] = -1;
        fDiJetMass[type][js][djs][0] = -1;
        fDiJetMass[type][js][djs][1] = -1;
        fDiJetMass[type][js][djs][2] = -1;
        fDiJetMass[type][js][djs][3] = -1;
        fDiJetMass[type][js][djs][4] = -1;
        if( dijet.GetEntries() < 2 ) continue;
        //TODO error check
        double mjj = dijet.InvM();
        double lpt = dijet.LeadingPt();
        double mT0  = dijet.Mt0(); // Mt^2 with Massless Jet = (pt1+pt2)^2 - PtPair^2
        double mT1  = dijet.Mt1(); // Mt^2 with Massive Jet  = ((m1^2+pt1^2)^.5+(m2^2+pt2^2)^.5)^2 - PtPair^2
        double mT2  = dijet.Mt2(); // Mt^2 with Massive Jet  = ((m1^2+pt1^2)^.5+(m2^2+pt2^2)^.5)^2 - PtPair^2
        fDiJetMass[type][js][djs][0] = mjj;
        fDiJetMass[type][js][djs][1] = mT0;
        fDiJetMass[type][js][djs][2] = mT1;
        fDiJetMass[type][js][djs][3] = mT2;
        fDiJetMass[type][js][djs][4] = lpt;
        fDiJetBin[type][js][djs][0] = fInvMBin.GetBin( mjj );
        fDiJetBin[type][js][djs][1] = fInvMBin.GetBin( mT0 );
        fDiJetBin[type][js][djs][2] = fInvMBin.GetBin( mT1 );
        fDiJetBin[type][js][djs][3] = fInvMBin.GetBin( mT2 );
        fDiJetBin[type][js][djs][4] = fJetPtBin.GetBin( lpt );
        fhInvMPttCor[type][js][djs]->Fill( mjj, lpt );
        fhMtMjjCor[type][js][djs][0]->Fill( mjj, mT1 );
        fhMtMjjCor[type][js][djs][1]->Fill( mT2, mT0 );
        fhMtMjjCor[type][js][djs][2]->Fill( mT2, mT1 );
        fhDiJetInvMInclu[type][js][djs]->Fill( mjj );
        fhDiJetMtInclu[type][js][djs]->Fill( mT1 );
        //if( type==0 ) continue;
        //cout<<Form("DEBUG_M1 : (%d,%d,%d) mT0=%10.4f mT1=%10.4f mT2=%10.4f",type,js,djs,mT0,mT1,mT2)<<endl;
      }
    }
  }
  //==== TYPE 
  for( int type=0;type<fJetListOfList.GetEntriesFast();type++ ){
    //==== Jet Selection
    for( int js=0;js<kJNJetSelection;js++ ){
      //==== DiJet Selection
      for( int djs=0;djs<kJNDiJetSelection;djs++ ){
        //==== Select Outgoing Jet without cut
        int ojs = js;int odjs = djs;
        //==== Get DiJet
        AliJDiJet & dijet = *GetDiJets( type, js, djs );
        if( dijet.GetEntries() < 2 ) continue;
        double pt1    = dijet(0).Pt();
        double pt2    = dijet(1).Pt();
        double ptpair = dijet.PtPair();
        double invm   = dijet.InvM();
        double dR     = dijet.DR();
        double dPhi   = dijet.DPhi();
        double dPhi2  = dPhi<0?dPhi+TMath::TwoPi():dPhi;
        double dEta   = dijet.DEta();
        double eAsymm = dijet.EAsymm();
        double ptAsymm= dijet.PtAsymm();
        double kTA    = dijet.kTA();
        double mT     = dijet.Mt();
        double mT0     = dijet.Mt0();
        double mT1     = dijet.Mt1();

        for( int i=0; i<12 ; i++ ){
          int bin = -1;
          int t,j,d,b=-1;
          if( i==4 ) { t=type;j=js;d=djs;b=0; }
          if( i==5 ) { t=type;j=js;d=djs;b=1; }
          if( i==6 ) { t=type;j=js;d=djs;b=2; }
          if( i==7 ) { t=type;j=js;d=djs;b=3; }
          //if( i==5 ) bin = fDiJetBin[type][ojs][odjs][1];
          if( js == 1 && i==8){
            AliJDiJet & dj = *GetDiJets( type, 0, djs );
            if( dj.GetEntries() < 2 ){ continue; }
            if( ( dj(0).GetID() == dijet(0).GetID() && dj(1).GetID()==dijet(1).GetID() ) 
                || ( dj(0).GetID() == dijet(1).GetID() && dj(1).GetID()==dijet(0).GetID() )
              ) { t=type;j=js;d=djs;b=0; }
          }
          if( js == 1 && i==9){
            AliJDiJet & dj = *GetDiJets( type, 0, djs );
            if( dj.GetEntries() < 2 ){ continue; }
            if( ( dj(0).GetID() == dijet(0).GetID() && dj(1).GetID()==dijet(1).GetID() ) 
                || ( dj(0).GetID() == dijet(1).GetID() && dj(1).GetID()==dijet(0).GetID() )
              ) { t=type;j=js;d=djs;b=1; }
          }
          if( js == 1 && i==10){
            AliJDiJet & dj = *GetDiJets( type, 0, djs );
            if( dj.GetEntries() < 2 ){ continue; }
            if( ( dj(0).GetID() == dijet(0).GetID() && dj(1).GetID()==dijet(1).GetID() ) 
                || ( dj(0).GetID() == dijet(1).GetID() && dj(1).GetID()==dijet(0).GetID() )
              ) { t=type;j=js;d=djs;b=2; }
          }
          if( js == 1 && i==11){
            AliJDiJet & dj = *GetDiJets( type, 0, djs );
            if( dj.GetEntries() < 2 ){ continue; }
            if( ( dj(0).GetID() == dijet(0).GetID() && dj(1).GetID()==dijet(1).GetID() ) 
                || ( dj(0).GetID() == dijet(1).GetID() && dj(1).GetID()==dijet(0).GetID() )
              ) { t=type;j=js;d=djs;b=3; }
          }
          if( b < 0 ) continue;
          bin = fDiJetBin[t][j][d][b];
          if( bin < 0 ) continue;
          double mass = fDiJetMass[type][js][djs][b];

          fhDiJetPtPair[type][js][djs][i][bin]->Fill( ptpair, ptpair<1e-8?1./1e-8:1./ptpair );
          fhDiJetPtPairRaw[type][js][djs][i][bin]->Fill( ptpair );
          fhDiJetPt1[type][js][djs][i][bin]->Fill(pt1);
          fhDiJetPt2[type][js][djs][i][bin]->Fill(pt2);
          fhDiJetKtA[type][js][djs][i][bin]->Fill( kTA );
          fhDiJetInvM[type][js][djs][i][bin]->Fill( mass );
          fhDiJetDeltaR[type][js][djs][i][bin]->Fill( dR );
          fhDiJetDeltaPhi[type][js][djs][i][bin]->Fill( dPhi2 );
          fhDiJetDeltaEta[type][js][djs][i][bin]->Fill( dEta );
          //cout<<"DEBUG_J1 : "<<ptAsymm<<endl;
          fhDiJetPtAsymm[type][js][djs][i][bin]->Fill( ptAsymm );
          fhDiJetEAsymm[type][js][djs][i][bin]->Fill( eAsymm );

          fhDiJetSingleJetMass[type][js][djs][i]->Fill( dijet(0).M() );
          fhDiJetSingleJetMass[type][js][djs][i]->Fill( dijet(1).M() );
          fhDiJetSingleJetArea[type][js][djs][i]->Fill( dijet(0).Area() );
          fhDiJetSingleJetArea[type][js][djs][i]->Fill( dijet(1).Area() );
          fhDiJetSingleJetNCont[type][js][djs][i]->Fill( dijet(0).GetNConstituentsRef() );
          fhDiJetSingleJetNCont[type][js][djs][i]->Fill( dijet(1).GetNConstituentsRef() );
          fhDiJetSingleJetActivity[type][js][djs][i]->Fill( dijet(0).LeadingParticlePt() / dijet(0).Pt() );
          fhDiJetSingleJetActivity[type][js][djs][i]->Fill( dijet(1).LeadingParticlePt() / dijet(1).Pt() );
          // TODO rapidity

        }
          AliJDiJet & dijet0 = * GetDiJets(1, js, 0);
          AliJDiJet & dijet1 = * GetDiJets(2, js, djs);
          for( int j0=0;j0<2;j0++){
            //int bin=-1;
            for( int j1=0;j1<dijet0.GetEntries();j1++ ){
              int rbin = fJetDRBin.GetBin( dijet(j0).DeltaR(dijet0(j1) ));
              if( rbin < 0 ) continue;
              fhJetEnergyComp[1][js][djs][rbin]->Fill( dijet(j0).E(), dijet0(j1).E() );
              fhJetPtComp[1][js][djs][rbin]->Fill( dijet(j0).Pt(), dijet0(j1).Pt() );
            }
            for( int j1=0;j1<dijet1.GetEntries();j1++ ){
              int rbin = fJetDRBin.GetBin( dijet(j0).DeltaR(dijet1(j1) ));
              if( rbin < 0 ) continue;
              fhJetEnergyComp[2][js][djs][rbin]->Fill( dijet(j0).E(), dijet1(j1).E() );
              fhJetPtComp[2][js][djs][rbin]->Fill( dijet(j0).Pt(), dijet1(j1).Pt() );
            }
          }
      }
    }
  }

  if(1)
    for( int js=0;js<kJNJetSelection;js++ ){
      for( int djs=0;djs<kJNDiJetSelection;djs++ ){
        for( int type0=0;type0<fJetListOfList.GetEntriesFast()-1;type0++ ){
          AliJDiJet & dijet0 = * GetDiJets(1, type0<2?0:js, type0<2?0:djs);
          if( dijet0.GetEntries() < 2 ) continue;
          for( int type1=type0+1;type1<fJetListOfList.GetEntriesFast();type1++ ){
            if( type0<2 && type1<2 && (js > 0 || djs >0) ) continue;
            AliJDiJet & dijet1 = * GetDiJets(1, type1<2?0:js, type1<2?0:djs);
            if( dijet1.GetEntries() < 2 ) continue;
            fhDiJetTypeCor[type0][type1][js][djs]->Fill( dijet0.InvM(), dijet1.InvM() );
          }
        }
      }
    }
}

/*
   void AliJDiJetAnalysis::FillFullChargedComparision( int js, int djs ){
   AliJDiJet & cjets = *GetDiJets( 4, js, djs );
   AliJDiJet & fjets = *GetDiJets( 3, js, djs );

   if( cjets.GetEntries() <1 && fjets.GetEntries() < 1 ) return;

   }
   */

void AliJDiJetAnalysis::FillPythiaDiJet(int js ){

  AliJDiJet & jetsIn  = *GetDiJets( 0, 0, 0 );
  AliJDiJet & jetsOut = *GetDiJets( 1, js, 0 );

  if( jetsIn.GetEntries() < 2 || jetsOut.GetEntries() < 2 ) return;

  AliJJet&  jet00 =  jetsIn(0);
  AliJJet&  jet01 =  jetsIn(1);
  AliJJet&  jet10 =  jetsOut(0);
  AliJJet&  jet11 =  jetsOut(1);

  AliJJet  jet0 = jet00+jet01;
  AliJJet  jet1 = jet10+jet11;

  AliJJet jetsum = jet0+jet1;
  int iM0 = fInvMBin.GetBin( jetsIn.InvM() ); 
  int iM1 = fInvMBin.GetBin( jetsOut.InvM() ); 

  if( iM0 > -1 && iM1 > -1 ){
    double ptpair0 = jetsIn.PtPair();
    double ptpair1 = jetsOut.PtPair();
    fhPythiaJetPtPair[0][iM0]->Fill( ptpair0 );
    fhPythiaJetPtPair[1][iM1]->Fill( ptpair1 );
    fhPythiaJetSum[0]->Fill( jetsum.P() );
    fhPythiaJetSum[1]->Fill( jetsum.E() );
  }
}

void AliJDiJetAnalysis::CreateHistos(){
  // Create Histograms

  //==== Hist Manager
  fHMG = new AliJHistManager( "AliJDiJetAnalysisHistManager", "AliJDiJetAnalysisHistManager");

  //==== BIN
  fJetPtPairBin   .Set("JetPtPair",   "P", "p_{Tpair}:%2.0f-%2.0f").SetBin( fCard->GetVector("Jet:PtPairBins") );
  fJetPtBin       .Set("JetPt",       "P", "p_{T}:%2.0f-%2.0f").SetBin( fCard->GetVector("Jet:PtBins") );
  fInvMBin        .Set("InvMB",       "M", "M:%2.0f-%2.0f").SetBin( fCard->GetVector("InvMBin"));
  fJetTypeBin     .Set("JetType",     "T", "T%1.0f", AliJBin::kSingle).SetBin( kJNJetTypeMax );
  fJetSelectionBin.Set("JetSelection","J", "J%1.0f", AliJBin::kSingle).SetBin( kJNJetSelection );
  fDiJetSelectionBin.Set("DiJetSelection","J", "J%1.0f",AliJBin::kSingle).SetBin( kJNDiJetSelection);

  fJetRBin        .Set("JetR",        "R", "R:%2.1f",AliJBin::kSingle ).SetBin( "0.4, 0.5");
  fDiJetBinTypeBin.Set("BinType",     "B", "B:%1.f", AliJBin::kSingle ).SetBin( 12 );
  fJetDRBin	    .Set("JetDR",   	"R", "R:%1.f" ).SetBin( "0.5, 1, 2, 4");
  fPYJetTypeBin.Set("PYJetType",   "T", "T%1.0f", AliJBin::kSingle).SetBin( 2 );
  fBin2		.Set("Bin2",	"B",  "B%1.0f", AliJBin::kSingle).SetBin( 2 );

  //==== Log Bins
  int nBINS2=100;
  double logBinsX2[nBINS2+2], limL2=0.1, LimH2=2000;
  double logBW2 = (log(LimH2)-log(limL2))/nBINS2;
  for(int ij=0;ij<=nBINS2;ij++) logBinsX2[ij+1]=limL2*exp(ij*logBW2);
  logBinsX2[0]=0;

  int nBINS3=100;
  double logBinsX3[nBINS3+2], limL3=0.1, LimH3=2000;
  double logBW3 = (log(LimH3)-log(limL3))/nBINS3;
  for(int ij=0;ij<=nBINS3;ij++) logBinsX3[ij+1]=limL3*exp(ij*logBW3);
  logBinsX3[0]=0;

  //==== Histogram
  int nDR = 1000;double xDR0= -10; double xDR1 = 10;
  //int nDPhi=1000;double xDPhi0=-TMath::Pi(); double xDPhi1=-xDPhi0;
  int nDPhi=1000;double xDPhi0=-1; double xDPhi1=TMath::TwoPi()+1;
  //== Jets QA
  fhJetPt 
    << TH1D("hJetPt","",nBINS2, logBinsX2 )
    << fJetTypeBin << fJetSelectionBin  <<"END";
  fhJetPhi 
    << TH1D("hJetPhi","",nDPhi, xDPhi0, xDPhi1 )
    << fJetTypeBin << fJetSelectionBin  <<"END";
  //fJetTypeBin.Print();
  //fJetSelectionBin.Print();
  fhJetEta 
    << TH1D("hJetEta","",nDR, xDR0, xDR1 )
    << fJetTypeBin << fJetSelectionBin  <<"END";
  fhJetNConst 
    << TH1D("hJetNConst","",100, 0-.5, 100-.5 )
    << fJetTypeBin << fJetSelectionBin  <<"END";
  fhJetNConst 
    << TH1D("hJetNConst","",nBINS2, logBinsX2 )
    << fJetTypeBin << fJetSelectionBin  <<"END";
  fhJetEnergy 
    << TH1D("hJetEnergy","",nBINS2, logBinsX2 )
    << fJetTypeBin << fJetSelectionBin  <<"END";
  fhJetInvM
    << TH1D("hJetInvM","",nBINS2, logBinsX2 )
    << fJetTypeBin << fJetSelectionBin  << fJetPtBin << "END";
  fhJetMjjEbin
    << TH1D("fhJetMjjEbin","",nBINS2, logBinsX2 )
    << fJetTypeBin << fJetSelectionBin  << fJetPtBin << "END";
  fhJetMtEbin
    << TH1D("fhJetMtEbin","",nBINS2, logBinsX2 )
    << fJetTypeBin << fJetSelectionBin  << fJetPtBin << "END";

  //== Jet Comparision
  fhJetEnergyComp
    << TProfile("hJetEnergyComp","",nBINS2, logBinsX2 )
    << fJetTypeBin << fJetSelectionBin << fDiJetSelectionBin << fJetDRBin <<"END";
  fhJetPtComp
    << TProfile("hJetPtComp","",nBINS2, logBinsX2 )
    << fJetTypeBin << fJetSelectionBin << fDiJetSelectionBin << fJetDRBin <<"END";
  fhJetDRToRef
    << TH1D("hJetDRToRef","", nDR, xDR0, xDR1 )
    << fJetTypeBin << fJetSelectionBin << fJetPtBin <<"END";
  fhJetDPhiToRef
    << TH1D("hJetDPhiToRef","", nDPhi, xDPhi0, xDPhi1 )
    << fJetTypeBin << fJetSelectionBin << fJetPtBin <<"END";



  //== DiJets
  fhDiJetPtPair // Inclusive PtPair
    << TH1D("hDiJetPtPair","",nBINS2, logBinsX2 )
    << fJetTypeBin << fJetSelectionBin << fDiJetSelectionBin << fDiJetBinTypeBin << fInvMBin<<"END";
  fhDiJetPtPairRaw  .SetWith( fhDiJetPtPair, "hDiJetPtPaiRaw" );
  fhDiJetPt1        .SetWith( fhDiJetPtPair, "hDiJetPt1" ); // Inclusive PtPair
  fhDiJetPt2        .SetWith( fhDiJetPtPair, "hDiJetPt2" ); // Inclusive PtPair
  fhDiJetInvM       .SetWith( fhDiJetPtPair, "hDiJetInvM" );
  fhDiJetKtA        .SetWith( fhDiJetPtPair, "hDiJetKtA");

    fhDiJetSingleJetMass      .SetWith( fhDiJetPtPair, "hDiJetSingleJetMass");
    fhDiJetSingleJetArea      .SetWith( fhDiJetPtPair, TH1D("fhDiJetSingleJetArea","",100, 0,1 ) );
  fhDiJetSingleJetActivity  .SetWith( fhDiJetSingleJetArea , "hDiJetSingleJetActivity");
    fhDiJetSingleJetNCont     .SetWith( fhDiJetPtPair, TH1D("fhDiJetSingleJetNCont","",100, -.5,100-.5 ) );

  fhDiJetDeltaR .SetWith( fhDiJetPtPair, TH1D("hDiJetDeltaR","",nDR, xDR0, xDR1 ) );// TODO change bins
  fhDiJetDeltaPhi
    << TH1D("hDiJetDeltaPhi","",nDPhi, xDPhi0, xDPhi1) // TODO change bins
    << fJetTypeBin << fJetSelectionBin << fDiJetSelectionBin << fDiJetBinTypeBin << fInvMBin<<"END";
  fhDiJetDeltaEta
    << TH1D("hDiJetDeltaEta","",nDR, xDR0, xDR1 ) // TODO change bins
    << fJetTypeBin << fJetSelectionBin << fDiJetSelectionBin << fDiJetBinTypeBin << fInvMBin<<"END";
  fhDiJetPtAsymm
    << TH1D("hDiJetPtAsymm","",100,0,1 ) // TODO change bins
    << fJetTypeBin << fJetSelectionBin << fDiJetSelectionBin << fDiJetBinTypeBin << fInvMBin<<"END";
  fhDiJetEAsymm
    << TH1D("hDiJetEAsymm","",100,0,1 ) // TODO change bins
    << fJetTypeBin << fJetSelectionBin << fDiJetSelectionBin << fDiJetBinTypeBin << fInvMBin<<"END";


  fhDiJetInvMInclu
    << TH1D("hDiJetInvM","",nBINS2, logBinsX2 )
    << fJetTypeBin << fJetSelectionBin << fDiJetSelectionBin <<"END";
  fhDiJetMtInclu    .SetWith( fhDiJetInvMInclu, "hDiJetMtInclu" );
  fhDiJetMultiplicity
    << TH1D("hDiJetMultiplicity","",nBINS2, logBinsX2)// TODO change bins
    << fJetTypeBin << fJetSelectionBin << "END";


  //== DiJet Correlation
  fhDiJetTypeCor
    << TH2D("hDiJetTypeCor","", nBINS3, logBinsX3, nBINS3, logBinsX3 )
    << fJetTypeBin << fJetTypeBin << fJetSelectionBin << fDiJetSelectionBin << "END";

  fhMtMjjCor
    << TH2D("hDiJetMtMjjCor","", nBINS3, logBinsX3, nBINS3, logBinsX3 )
    << fJetTypeBin << fJetSelectionBin << fDiJetSelectionBin << 3 <<"END";
  fhInvMPttCor
    << TH2D("hDiJetInvMPttCor","", nBINS3, logBinsX3, nBINS3, logBinsX3 )
    << fJetTypeBin << fJetSelectionBin << fDiJetSelectionBin <<"END";

  //== PYTHIA
  fhPythiaJetPtPair
    << TH1D("hPythiaJetPtPair", "", nBINS2, logBinsX2 )
    << fPYJetTypeBin << fInvMBin << "END";
  fhPythiaJetSum
    << TH1D("hPythiaJetSum", "", 1000,-10,10 )
    << fBin2 << fInvMBin << "END";

  //fHMG->Print();
  fHMG->WriteConfig();
}



