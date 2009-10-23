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
//
// QA class of Heavy Flavor quark and fragmeted/decayed particles
// -Check kinematics of Heavy Quarks/hadrons, and decayed leptons
//    pT, rapidity
//    decay lepton kinematics w/wo acceptance
//    heavy hadron decay length, electron pT fraction carried from decay
// -Check yield of Heavy Quarks/hadrons
//    Number of produced heavy quark
//    Number of produced hadron of given pdg code
//
//
// Authors:
//   MinJung Kweon <minjung@physi.uni-heidelberg.de>
//

#include <TH2F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TParticle.h>

#include <AliLog.h>
#include <AliStack.h>
#include <AliAODMCParticle.h>

#include "AliHFEmcQA.h"

const Int_t AliHFEmcQA::fgkGluon=21;
const Int_t AliHFEmcQA::fgkMaxGener=10;
const Int_t AliHFEmcQA::fgkMaxIter=100;
const Int_t AliHFEmcQA::fgkqType = 7;    // number of species waiting for QA done

ClassImp(AliHFEmcQA)

//_______________________________________________________________________________________________
AliHFEmcQA::AliHFEmcQA() : 
        fStack(0x0) 
        ,fMCArray(0x0)
        ,fNparents(0) 
{
        // Default constructor
        
}

//_______________________________________________________________________________________________
AliHFEmcQA::AliHFEmcQA(const AliHFEmcQA&p):
        TObject(p)
        ,fStack(0x0) 
        ,fMCArray(0x0)
        ,fNparents(p.fNparents) 
{
        // Copy constructor
}

//_______________________________________________________________________________________________
AliHFEmcQA&
AliHFEmcQA::operator=(const AliHFEmcQA &)
{
  // Assignment operator

  AliInfo("Not yet implemented.");
  return *this;
}

//_______________________________________________________________________________________________
AliHFEmcQA::~AliHFEmcQA()
{
        // Destructor

        AliInfo("Analysis Done.");
}

//_______________________________________________________________________________________________
void AliHFEmcQA::PostAnalyze() const
{
}

//__________________________________________
void AliHFEmcQA::CreateHistograms(const Int_t kquark, Int_t icut, TString hnopt) 
{
  // create histograms

  if (kquark != kCharm && kquark != kBeauty) {
    AliDebug(1, "This task is only for heavy quark QA, return\n");
    return; 
  }
  Int_t iq = kquark - kCharm; 

  TString kqTypeLabel[fgkqType];
  if (kquark == kCharm){
    kqTypeLabel[kQuark]="c";
    kqTypeLabel[kantiQuark]="cbar";
    kqTypeLabel[kHadron]="cHadron";
    kqTypeLabel[keHadron]="ceHadron";
    kqTypeLabel[kDeHadron]="nullHadron";
    kqTypeLabel[kElectron]="ce";
    kqTypeLabel[kElectron2nd]="nulle";
  } else if (kquark == kBeauty){
    kqTypeLabel[kQuark]="b";
    kqTypeLabel[kantiQuark]="bbar";
    kqTypeLabel[kHadron]="bHadron";
    kqTypeLabel[keHadron]="beHadron";
    kqTypeLabel[kDeHadron]="bDeHadron";
    kqTypeLabel[kElectron]="be";
    kqTypeLabel[kElectron2nd]="bce";
  }


  TString hname; 
  for (Int_t iqType = 0; iqType < fgkqType; iqType++ ){
     if (iqType < keHadron && icut > 0) continue; // don't duplicate histogram for quark and hadron
     hname = hnopt+"PdgCode_"+kqTypeLabel[iqType];
     fHist[iq][iqType][icut].fPdgCode = new TH1F(hname,hname,20001,-10000.5,10000.5);
     hname = hnopt+"Pt_"+kqTypeLabel[iqType];
     fHist[iq][iqType][icut].fPt = new TH1F(hname,hname+";p_{T} (GeV/c)",250,0,50);
     hname = hnopt+"Y_"+kqTypeLabel[iqType];
     fHist[iq][iqType][icut].fY = new TH1F(hname,hname,150,-7.5,7.5);
     hname = hnopt+"Eta_"+kqTypeLabel[iqType];
     fHist[iq][iqType][icut].fEta = new TH1F(hname,hname,150,-7.5,7.5);
  }

  if (icut == 0){ 
    hname = hnopt+"Nq_"+kqTypeLabel[kQuark];
    fHistComm[iq][icut].fNq = new TH1F(hname,hname,10,-0.5,9.5);
    hname = hnopt+"ProcessID_"+kqTypeLabel[kQuark];
    fHistComm[iq][icut].fProcessID = new TH1F(hname,hname,21,-10.5,10.5);
  }
  hname = hnopt+"ePtRatio_"+kqTypeLabel[kQuark];
  fHistComm[iq][icut].fePtRatio = new TH2F(hname,hname+";p_{T} (GeV/c);momentum fraction",100,0,30,100,0,1);
  hname = hnopt+"DePtRatio_"+kqTypeLabel[kQuark];
  fHistComm[iq][icut].fDePtRatio = new TH2F(hname,hname+";p_{T} (GeV/c);momentum fraction",100,0,30,100,0,1);
  hname = hnopt+"eDistance_"+kqTypeLabel[kQuark];
  fHistComm[iq][icut].feDistance= new TH2F(hname,hname+";p_{T} (GeV/c);distance (cm)",100,0,30,200,0,2);
  hname = hnopt+"DeDistance_"+kqTypeLabel[kQuark];
  fHistComm[iq][icut].fDeDistance= new TH2F(hname,hname+";p_{T} (GeV/c);distance (cm)",100,0,30,200,0,2);

}

//__________________________________________
void AliHFEmcQA::Init()
{
  // called at begining every event
  
  for (Int_t i=0; i<2; i++){
     fIsHeavy[i] = 0;
  } 

  fNparents = 7;

  fParentSelect[0][0] =  411;  
  fParentSelect[0][1] =  421;
  fParentSelect[0][2] =  431;
  fParentSelect[0][3] = 4122;
  fParentSelect[0][4] = 4132;
  fParentSelect[0][5] = 4232;
  fParentSelect[0][6] = 4332;

  fParentSelect[1][0] =  511;
  fParentSelect[1][1] =  521;
  fParentSelect[1][2] =  531;
  fParentSelect[1][3] = 5122;
  fParentSelect[1][4] = 5132;
  fParentSelect[1][5] = 5232;
  fParentSelect[1][6] = 5332;

}

//__________________________________________
void AliHFEmcQA::GetQuarkKine(TParticle *part, Int_t iTrack, const Int_t kquark) 
{
  // get heavy quark kinematics

    if (kquark != kCharm && kquark != kBeauty) {
      AliDebug(1, "This task is only for heavy quark QA, return\n");
      return; 
    }
    Int_t iq = kquark - kCharm; 

    if (iTrack < 0) { 
      AliDebug(1, "Stack label is negative, return\n");
      return; 
    }

    //TParticle *part = fStack->Particle(iTrack); 
    Int_t partPdgcode = TMath::Abs(part->GetPdgCode());

    // select heavy hadron or not fragmented heavy quark 
    if ( int(partPdgcode/100.)==kquark || int(partPdgcode/1000.)==kquark || (partPdgcode==kquark && (part->GetNDaughters()==0 && iTrack>5)) ){ 

      TParticle *partMother;
      Int_t iLabel;

      if (partPdgcode == kquark){ // in case of not fragmented heavy quark  
        partMother = part; 
        iLabel = iTrack;
      } else{ // in case of heavy hadron, start to search for mother heavy parton 
        iLabel = part->GetFirstMother(); 
        if (iLabel>-1) { partMother = fStack->Particle(iLabel); }
        else {
          AliDebug(1, "Stack label is negative, return\n");
          return; 
        }
      }

      // heavy parton selection as a mother of heavy hadron 
      // if the heavy particle comes from string which is denoted as particle status 12|12|12...12|11,[PYTHIA p.60]
      // in this case, the mother of heavy particle can be one of the fragmented parton of the string
      // should I make a condition that partMother should be quark or diquark? -> not necessary
      if ( abs(partMother->GetPdgCode()) == kquark || (partMother->GetStatusCode() == 11) ){
      //if ( abs(partMother->GetPdgCode()) == kquark || (partMother->GetStatusCode() == 11 || partMother->GetStatusCode() == 12) ){

        if ( abs(partMother->GetPdgCode()) != kquark ){
          // search fragmented partons in the same string
          Bool_t isSameString = kTRUE; 
          for (Int_t i=1; i<fgkMaxIter; i++){
             iLabel = iLabel - 1;
             if (iLabel>-1) { partMother = fStack->Particle(iLabel); }
             else {
               AliDebug(1, "Stack label is negative, return\n");
               return; 
             }
             if ( abs(partMother->GetPdgCode()) == kquark ) break;
             if ( partMother->GetStatusCode() != 12 ) isSameString = kFALSE;
             if (!isSameString) return; 
          }
        }
        AliDebug(1, "Can not find heavy parton of this heavy hadron in the string, return\n");
        if (abs(partMother->GetPdgCode()) != kquark) return; 

        if (fIsHeavy[iq] >= 50) return;  
        fHeavyQuark[fIsHeavy[iq]] = partMother;
        fIsHeavy[iq]++;

        // fill kinematics for heavy parton
        if (partMother->GetPdgCode() > 0) { // quark
          fHist[iq][kQuark][0].fPdgCode->Fill(partMother->GetPdgCode());
          fHist[iq][kQuark][0].fPt->Fill(partMother->Pt());
          fHist[iq][kQuark][0].fY->Fill(GetRapidity(partMother));
          fHist[iq][kQuark][0].fEta->Fill(partMother->Eta());
        } else{ // antiquark
          fHist[iq][kantiQuark][0].fPdgCode->Fill(partMother->GetPdgCode());
          fHist[iq][kantiQuark][0].fPt->Fill(partMother->Pt());
          fHist[iq][kantiQuark][0].fY->Fill(GetRapidity(partMother));
          fHist[iq][kantiQuark][0].fEta->Fill(partMother->Eta());
        }

      } // end of heavy parton slection loop 

    } // end of heavy hadron or quark selection

}

//__________________________________________
void AliHFEmcQA::EndOfEventAna(const Int_t kquark)
{
  // end of event analysis

  if (kquark != kCharm && kquark != kBeauty) {
    AliDebug(1, "This task is only for heavy quark QA, return\n");
    return; 
  }
  Int_t iq = kquark - kCharm; 


  // # of heavy quark per event
  AliDebug(1,Form("Number of heavy quark in this event = %d \n",fIsHeavy[iq]));
  fHistComm[iq][0].fNq->Fill(fIsHeavy[iq]);

  Int_t motherID[fgkMaxGener];
  Int_t motherType[fgkMaxGener];
  Int_t motherLabel[fgkMaxGener];
  Int_t ancestorPdg[fgkMaxGener];
  Int_t ancestorLabel[fgkMaxGener];

  for (Int_t i = 0; i < fgkMaxGener; i++){ // initialization
     motherID[i] = 0;
     motherType[i] = 0;
     motherLabel[i] = 0;
     ancestorPdg[i] = 0;
     ancestorLabel[i] = 0;
  }

  // check history of found heavy quarks
  for (Int_t i = 0; i < fIsHeavy[iq]; i++){

     ancestorLabel[0] = i;
     ancestorPdg[0] = fHeavyQuark[i]->GetPdgCode(); 
     ancestorLabel[1] = fHeavyQuark[i]->GetFirstMother(); 

     AliDebug(1,Form("pdg code= %d\n",ancestorPdg[0]));
     AliDebug(1,Form("ancestor label= %d\n",ancestorLabel[1]));

     Int_t ig = 1;
     while (ancestorLabel[ig] != -1){
          // in case there is mother, get mother's pdg code and grandmother's label
          IdentifyMother(ancestorLabel[ig], ancestorPdg[ig], ancestorLabel[ig+1]); 
          // if mother is still heavy, find again mother's ancestor
          if (ancestorPdg[ig-1] == ancestorPdg[ig]) {
            ig++;
            continue; // if it is from same heavy
          }
          // if the heavy's mother is not heavy, check the mother's label to know if it comes from inital or final parton shower
          if (IsFromInitialShower(ancestorLabel[ig],motherID[i],motherType[i],motherLabel[i])) break;
          if (IsFromFinalParton(ancestorLabel[ig],motherID[i],motherType[i],motherLabel[i])) break;
          // if it is not the above case, something is strange
          ReportStrangeness(motherID[i],motherType[i],motherLabel[i]);
          break;
     } 
     if (ancestorLabel[ig] == -1){ // from hard scattering
       HardScattering(kquark, motherID[i],motherType[i], motherLabel[i]);
     }

  } // end of found heavy quark loop


  // check process type
  Int_t processID = 0;
  for (Int_t i = 0; i < fIsHeavy[iq]; i++){
     AliDebug(1,Form("Mother ID= %d type= %d label= %d\n",motherID[i],motherType[i],motherLabel[i]));
  }


  Int_t nheavypair = Int_t(fIsHeavy[iq]/2.); 
  for (Int_t ipair = 0; ipair < nheavypair; ipair++){

     Int_t id1 = ipair*2;
     Int_t id2 = ipair*2 + 1;

     if (motherType[id1] == 2 && motherType[id2] == 2){
       if (motherLabel[id1] == motherLabel[id2]) processID = kGluonSplitting; // gluon spliting
       else processID = -9;
     }
     else if (motherType[id1] == -1 && motherType[id2] == -1) {
       if (motherLabel[id1] == -1 && motherLabel[id2] == -1) {
         if (motherID[id1] == fgkGluon) processID = kPairCreationFromg; // gluon fusion
         else processID = kPairCreationFromq; // q-qbar pair creation
       }
       else processID = -8;
     }
     else if (motherType[id1] == -1 || motherType[id2] == -1) {
       if ((motherLabel[id1] == -1 || motherLabel[id2] == -1) && (motherLabel[id1]*motherLabel[id2] == -2 || motherLabel[id1]*motherLabel[id2] == -3)) {
         if(motherID[id1]*motherID[id2] == kquark*fgkGluon) processID = kFlavourExitation; // flavour exitation 
         else processID = kLightQuarkShower;
       }
       else processID = -7;
     }
     else if (motherType[id1] == -2 || motherType[id2] == -2) {
       if (motherLabel[id1] == motherLabel[id2]) processID = kInitialPartonShower; // initial parton shower
       else processID = -6;
       
     }
     else processID = -5;

     if (nheavypair >1) AliDebug(1,Form("Multi pair found : process ID = %d\n",processID));
     else fHistComm[iq][0].fProcessID->Fill(processID);
     AliDebug(1,Form("Process ID = %d\n",processID));
  } // end of # heavy quark pair loop

}

//__________________________________________
void AliHFEmcQA::GetHadronKine(TParticle* mcpart, const Int_t kquark)
{
    // decay electron kinematics

    if (kquark != kCharm && kquark != kBeauty) {
      AliDebug(1, "This task is only for heavy quark QA, return\n");
      return;
    }
    Int_t iq = kquark - kCharm;

    //TParticle* mcpart = fStack->Particle(iTrack);

    Int_t iLabel = mcpart->GetFirstMother();
    if (iLabel<0){
      AliDebug(1, "Stack label is negative, return\n");
      return;
    }

    TParticle *partCopy = mcpart;
    Int_t pdgcode = mcpart->GetPdgCode();
    Int_t pdgcodeCopy = pdgcode;

    // if the mother is charmed hadron  
    Bool_t isDirectCharm = kFALSE;
    if ( int(abs(pdgcode)/100.) == kCharm || int(abs(pdgcode)/1000.) == kCharm ) {

          // iterate until you find B hadron as a mother or become top ancester 
          for (Int_t i=1; i<fgkMaxIter; i++){

             Int_t jLabel = mcpart->GetFirstMother();
             if (jLabel == -1){
               isDirectCharm = kTRUE;
               break; // if there is no ancester
             }
             if (jLabel < 0){ // safety protection
               AliDebug(1, "Stack label is negative, return\n");
               return;
             }
             // if there is an ancester
             TParticle* mother = fStack->Particle(jLabel);
             Int_t motherPDG = mother->GetPdgCode();
    
             for (Int_t j=0; j<fNparents; j++){
                if (abs(motherPDG)==fParentSelect[1][j]) return; // return if this hadron is originated from b
             }

             mcpart = mother;
          } // end of iteration 
    } // end of if
    if((isDirectCharm == kTRUE && kquark == kCharm) || kquark == kBeauty) {
         for (Int_t i=0; i<fNparents; i++){
            if (abs(pdgcodeCopy)==fParentSelect[iq][i]){

              // fill hadron kinematics
              fHist[iq][kHadron][0].fPdgCode->Fill(pdgcodeCopy);
              fHist[iq][kHadron][0].fPt->Fill(partCopy->Pt());
              fHist[iq][kHadron][0].fY->Fill(GetRapidity(partCopy));
              fHist[iq][kHadron][0].fEta->Fill(partCopy->Eta());

            }
         }
    } // end of if
}

//__________________________________________
void AliHFEmcQA::GetDecayedKine(TParticle* mcpart, const Int_t kquark, Int_t kdecayed, Int_t icut) 
{
    // decay electron kinematics
    
    if (kquark != kCharm && kquark != kBeauty) {
      AliDebug(1, "This task is only for heavy quark QA, return\n");
      return; 
    }
    Int_t iq = kquark - kCharm; 
    Bool_t isFinalOpenCharm = kFALSE;

/*
    if (iTrack < 0) { 
      AliDebug(1, "Stack label is negative, return\n");
      return; 
    }
    */

    //TParticle* mcpart = fStack->Particle(iTrack);

    if ( abs(mcpart->GetPdgCode()) != kdecayed ) return;

    Int_t iLabel = mcpart->GetFirstMother(); 
    if (iLabel<0){
      AliDebug(1, "Stack label is negative, return\n");
      return; 
    }

    TParticle *partMother = fStack->Particle(iLabel);
    TParticle *partMotherCopy = partMother;
    Int_t maPdgcode = partMother->GetPdgCode();
    Int_t maPdgcodeCopy = maPdgcode;

    // get electron production vertex   
    TLorentzVector ePoint;
    mcpart->ProductionVertex(ePoint);

    // if the mother is charmed hadron  
    Bool_t isMotherDirectCharm = kFALSE;
    if ( int(abs(maPdgcode)/100.) == kCharm || int(abs(maPdgcode)/1000.) == kCharm ) { 

         for (Int_t i=0; i<fNparents; i++){
            if (abs(maPdgcode)==fParentSelect[0][i]){
              isFinalOpenCharm = kTRUE;
            } 
         }  
         if (!isFinalOpenCharm) return ;

          // iterate until you find B hadron as a mother or become top ancester 
          for (Int_t i=1; i<fgkMaxIter; i++){

             Int_t jLabel = partMother->GetFirstMother(); 
             if (jLabel == -1){
               isMotherDirectCharm = kTRUE;
               break; // if there is no ancester
             }
             if (jLabel < 0){ // safety protection
               AliDebug(1, "Stack label is negative, return\n");
               return; 
             }

             // if there is an ancester
             TParticle* grandMa = fStack->Particle(jLabel);
             Int_t grandMaPDG = grandMa->GetPdgCode();

             for (Int_t j=0; j<fNparents; j++){
                if (abs(grandMaPDG)==fParentSelect[1][j]){

                  if (kquark == kCharm) return;
                  // fill electron kinematics
                  fHist[iq][kElectron2nd][icut].fPdgCode->Fill(mcpart->GetPdgCode());
                  fHist[iq][kElectron2nd][icut].fPt->Fill(mcpart->Pt());
                  fHist[iq][kElectron2nd][icut].fY->Fill(GetRapidity(mcpart));
                  fHist[iq][kElectron2nd][icut].fEta->Fill(mcpart->Eta());

                  // fill mother hadron kinematics
                  fHist[iq][kDeHadron][icut].fPdgCode->Fill(grandMaPDG); 
                  fHist[iq][kDeHadron][icut].fPt->Fill(grandMa->Pt());
                  fHist[iq][kDeHadron][icut].fY->Fill(GetRapidity(grandMa));
                  fHist[iq][kDeHadron][icut].fEta->Fill(grandMa->Eta());
                 
                  // ratio between pT of electron and pT of mother B hadron 
                  if(grandMa->Pt()) fHistComm[iq][icut].fDePtRatio->Fill(grandMa->Pt(),mcpart->Pt()/grandMa->Pt());

                  // distance between electron production point and mother hadron production point
                  TLorentzVector debPoint;
                  grandMa->ProductionVertex(debPoint);
                  TLorentzVector dedistance = ePoint - debPoint;
                  fHistComm[iq][icut].fDeDistance->Fill(grandMa->Pt(),dedistance.P());
                  return;
                }
             } 

             partMother = grandMa;
          } // end of iteration 
    } // end of if
    if((isMotherDirectCharm == kTRUE && kquark == kCharm) || kquark == kBeauty) {
         for (Int_t i=0; i<fNparents; i++){
            if (abs(maPdgcodeCopy)==fParentSelect[iq][i]){

              // fill electron kinematics
              fHist[iq][kElectron][icut].fPdgCode->Fill(mcpart->GetPdgCode());
              fHist[iq][kElectron][icut].fPt->Fill(mcpart->Pt());
              fHist[iq][kElectron][icut].fY->Fill(GetRapidity(mcpart));
              fHist[iq][kElectron][icut].fEta->Fill(mcpart->Eta());  

              // fill mother hadron kinematics
              fHist[iq][keHadron][icut].fPdgCode->Fill(maPdgcodeCopy); 
              fHist[iq][keHadron][icut].fPt->Fill(partMotherCopy->Pt());
              fHist[iq][keHadron][icut].fY->Fill(GetRapidity(partMotherCopy));
              fHist[iq][keHadron][icut].fEta->Fill(partMotherCopy->Eta());

              // ratio between pT of electron and pT of mother B hadron 
              if(partMotherCopy->Pt()) fHistComm[iq][icut].fePtRatio->Fill(partMotherCopy->Pt(),mcpart->Pt()/partMotherCopy->Pt());

              // distance between electron production point and mother hadron production point
              TLorentzVector ebPoint;
              partMotherCopy->ProductionVertex(ebPoint);
              TLorentzVector edistance = ePoint - ebPoint;
              fHistComm[iq][icut].feDistance->Fill(partMotherCopy->Pt(),edistance.P());
            }
         }
    } // end of if
}

//____________________________________________________________________
void  AliHFEmcQA::GetDecayedKine(AliAODMCParticle *mcpart, const Int_t kquark, Int_t kdecayed, Int_t icut)
{
  // decay electron kinematics

  if (kquark != kCharm && kquark != kBeauty) {
    AliDebug(1, "This task is only for heavy quark QA, return\n");
    return;
  }

  Int_t iq = kquark - kCharm;
  Bool_t isFinalOpenCharm = kFALSE;

  if ( abs(mcpart->GetPdgCode()) != kdecayed ) return;

  // mother
  Int_t iLabel = mcpart->GetMother();
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return;
  }

  AliAODMCParticle *partMother = (AliAODMCParticle*)fMCArray->At(iLabel);
  AliAODMCParticle *partMotherCopy = partMother;
  Int_t maPdgcode = partMother->GetPdgCode();
  Int_t maPdgcodeCopy = maPdgcode;

  Bool_t isMotherDirectCharm = kFALSE;
  if ( int(abs(maPdgcode)/100.) == kCharm || int(abs(maPdgcode)/1000.) == kCharm ) {

    for (Int_t i=0; i<fNparents; i++){
       if (abs(maPdgcode)==fParentSelect[0][i]){
         isFinalOpenCharm = kTRUE;
       }
    } 
    if (!isFinalOpenCharm) return;

    for (Int_t i=1; i<fgkMaxIter; i++){

       Int_t jLabel = partMother->GetMother();
       if (jLabel == -1){
         isMotherDirectCharm = kTRUE;
         break; // if there is no ancester
       }
       if (jLabel < 0){ // safety protection
         AliDebug(1, "Stack label is negative, return\n");
         return;
       }

       // if there is an ancester
       AliAODMCParticle* grandMa = (AliAODMCParticle*)fMCArray->At(jLabel);
       Int_t grandMaPDG = grandMa->GetPdgCode();

       for (Int_t j=0; j<fNparents; j++){
          if (abs(grandMaPDG)==fParentSelect[1][j]){

            if (kquark == kCharm) return;
            // fill electron kinematics
            fHist[iq][kElectron2nd][icut].fPdgCode->Fill(mcpart->GetPdgCode());
            fHist[iq][kElectron2nd][icut].fPt->Fill(mcpart->Pt());
            fHist[iq][kElectron2nd][icut].fY->Fill(GetRapidity(mcpart));
            fHist[iq][kElectron2nd][icut].fEta->Fill(mcpart->Eta());

            // fill mother hadron kinematics
            fHist[iq][kDeHadron][icut].fPdgCode->Fill(grandMaPDG);
            fHist[iq][kDeHadron][icut].fPt->Fill(grandMa->Pt());
            fHist[iq][kDeHadron][icut].fY->Fill(GetRapidity(grandMa));
            fHist[iq][kDeHadron][icut].fEta->Fill(grandMa->Eta());

            return;
          }
       }

       partMother = grandMa;
    } // end of iteration 
  } // end of if
  if ((isMotherDirectCharm == kTRUE && kquark == kCharm) || kquark == kBeauty) {
    for (Int_t i=0; i<fNparents; i++){
       if (abs(maPdgcodeCopy)==fParentSelect[iq][i]){

         // fill electron kinematics
         fHist[iq][kElectron][icut].fPdgCode->Fill(mcpart->GetPdgCode());
         fHist[iq][kElectron][icut].fPt->Fill(mcpart->Pt());
         fHist[iq][kElectron][icut].fY->Fill(GetRapidity(mcpart));
         fHist[iq][kElectron][icut].fEta->Fill(mcpart->Eta());

         // fill mother hadron kinematics
         fHist[iq][keHadron][icut].fPdgCode->Fill(maPdgcodeCopy);
         fHist[iq][keHadron][icut].fPt->Fill(partMotherCopy->Pt());
         fHist[iq][keHadron][icut].fY->Fill(GetRapidity(partMotherCopy));
         fHist[iq][keHadron][icut].fEta->Fill(partMotherCopy->Eta());

       }
    }
  } // end of if

}

//__________________________________________
void AliHFEmcQA::IdentifyMother(Int_t motherlabel, Int_t &motherpdg, Int_t &grandmotherlabel)
{
       // find mother pdg code and label 

       if (motherlabel < 0) { 
         AliDebug(1, "Stack label is negative, return\n");
         return; 
       }
       TParticle *heavysMother = fStack->Particle(motherlabel);
       motherpdg = heavysMother->GetPdgCode();
       grandmotherlabel = heavysMother->GetFirstMother();
       AliDebug(1,Form("ancestor pdg code= %d\n",motherpdg));
}

//__________________________________________
void AliHFEmcQA::HardScattering(const Int_t kquark, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
       // mothertype -1 means this heavy quark coming from hard vertex

       TParticle *afterinitialrad1  = fStack->Particle(4);
       TParticle *afterinitialrad2  = fStack->Particle(5);
           
       motherlabel = -1;

       if (abs(afterinitialrad1->GetPdgCode()) == fgkGluon && abs(afterinitialrad2->GetPdgCode()) == fgkGluon){
         AliDebug(1,"heavy from gluon gluon pair creation!\n");
         mothertype = -1;
         motherID = fgkGluon;
       }
       else if (abs(afterinitialrad1->GetPdgCode()) == kquark || abs(afterinitialrad2->GetPdgCode()) == kquark){ // one from Q and the other from g
         AliDebug(1,"heavy from flavor exitation!\n");
         mothertype = -1;
         motherID = kquark;
       }
       else if  (abs(afterinitialrad1->GetPdgCode()) == abs(afterinitialrad2->GetPdgCode())){
         AliDebug(1,"heavy from q-qbar pair creation!\n");
         mothertype = -1;
         motherID = 1;
       }
       else {
         AliDebug(1,"something strange!\n");
         mothertype = -999;
         motherlabel = -999;
         motherID = -999;
       }
}

//__________________________________________
Bool_t AliHFEmcQA::IsFromInitialShower(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
       // mothertype -2 means this heavy quark coming from initial state 

       if (inputmotherlabel==2 || inputmotherlabel==3){ // mother exist before initial state radiation
         TParticle *heavysMother = fStack->Particle(inputmotherlabel);
         motherID = heavysMother->GetPdgCode(); 
         mothertype = -2; // there is mother before initial state radiation
         motherlabel = inputmotherlabel;
         AliDebug(1,"initial parton shower! \n");

         return kTRUE;
       }

       return kFALSE;
}

//__________________________________________
Bool_t AliHFEmcQA::IsFromFinalParton(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
       // mothertype 2 means this heavy quark coming from final state 

       if (inputmotherlabel > 5){ // mother exist after hard scattering
         TParticle *heavysMother = fStack->Particle(inputmotherlabel);
         motherID = heavysMother->GetPdgCode(); 
         mothertype = 2; // 
         motherlabel = inputmotherlabel;
         AliDebug(1,Form("heavy quark from %d after hard scattering! \n",motherID));

         return kTRUE;
       }
       return kFALSE;
}

//__________________________________________
void AliHFEmcQA::ReportStrangeness(Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
      // mark strange behavior  

       mothertype = -888;
       motherlabel = -888;
       motherID = -888;
       AliDebug(1,"something strange!\n");
}

//__________________________________________
Float_t AliHFEmcQA::GetRapidity(TParticle *part) const
{
      // return rapidity

       Float_t rapidity;        
       if(!((part->Energy() - part->Pz())*(part->Energy() + part->Pz())>0)) rapidity=-999; 
       else rapidity = 0.5*(TMath::Log((part->Energy()+part->Pz()) / (part->Energy()-part->Pz()))); 
       return rapidity;
}

//__________________________________________
Float_t AliHFEmcQA::GetRapidity(AliAODMCParticle *part) const
{
      // return rapidity

       Float_t rapidity;        
       if(!((part->E() - part->Pz())*(part->E() + part->Pz())>0)) rapidity=-999; 
       else rapidity = 0.5*(TMath::Log((part->E()+part->Pz()) / (part->E()-part->Pz()))); 
       return rapidity;
}

//__________________________________________
Int_t AliHFEmcQA::GetElectronSource(AliAODMCParticle *mcpart)
{        
  // decay electron's origin 

  if ( abs(mcpart->GetPdgCode()) != AliHFEmcQA::kElectronPDG ) return -1;
       
  Int_t origin = -1;
  Bool_t isFinalOpenCharm = kFALSE;

  // mother
  Int_t iLabel = mcpart->GetMother();
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return -1;
  } 
       
  AliAODMCParticle *partMother = (AliAODMCParticle*)fMCArray->At(iLabel);
  Int_t maPdgcode = partMother->GetPdgCode();
  
  // if the mother is charmed hadron  
  if ( int(abs(maPdgcode)/100.) == kCharm || int(abs(maPdgcode)/1000.) == kCharm ) {
    
    for (Int_t i=0; i<fNparents; i++){
       if (abs(maPdgcode)==fParentSelect[0][i]){
         isFinalOpenCharm = kTRUE;
       }
    }
    if (!isFinalOpenCharm) return -1;

    // iterate until you find B hadron as a mother or become top ancester 
    for (Int_t i=1; i<fgkMaxIter; i++){

       Int_t jLabel = partMother->GetMother();
       if (jLabel == -1){
         origin = kDirectCharm;
         return origin;
       }
       if (jLabel < 0){ // safety protection
         AliDebug(1, "Stack label is negative, return\n");
         return -1;
       }

       // if there is an ancester
       AliAODMCParticle* grandMa = (AliAODMCParticle*)fMCArray->At(jLabel);
       Int_t grandMaPDG = grandMa->GetPdgCode();

       for (Int_t i=0; i<fNparents; i++){
          if (abs(grandMaPDG)==fParentSelect[1][i]){
            origin = kBeautyCharm;
            return origin;
          }
       }

       partMother = grandMa;
    } // end of iteration 
  } // end of if
  else if ( int(abs(maPdgcode)/100.) == kBeauty || int(abs(maPdgcode)/1000.) == kBeauty ) {
    for (Int_t i=0; i<fNparents; i++){
       if (abs(maPdgcode)==fParentSelect[1][i]){
         origin = kDirectBeauty;
         return origin;
       }
    }
  } // end of if
  else if ( abs(maPdgcode) == 22 ) {
    origin = kGamma;
    return origin;
  } // end of if
  else if ( abs(maPdgcode) == 111 ) {
    origin = kPi0;
    return origin;
  } // end of if

  return origin;
}

//__________________________________________
Int_t AliHFEmcQA::GetElectronSource(TParticle* mcpart)
{
  // decay electron's origin 

  if ( abs(mcpart->GetPdgCode()) != AliHFEmcQA::kElectronPDG ) return -1;

  Int_t origin = -1;
  Bool_t isFinalOpenCharm = kFALSE;

  Int_t iLabel = mcpart->GetFirstMother();
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return -1;
  }

  TParticle *partMother = fStack->Particle(iLabel);
  Int_t maPdgcode = partMother->GetPdgCode();

   // if the mother is charmed hadron  
   if ( int(abs(maPdgcode)/100.) == kCharm || int(abs(maPdgcode)/1000.) == kCharm ) {

     for (Int_t i=0; i<fNparents; i++){
        if (abs(maPdgcode)==fParentSelect[0][i]){
          isFinalOpenCharm = kTRUE;
        }
     }
     if (!isFinalOpenCharm) return -1;

     // iterate until you find B hadron as a mother or become top ancester 
     for (Int_t i=1; i<fgkMaxIter; i++){

        Int_t jLabel = partMother->GetFirstMother();
        if (jLabel == -1){
          origin = kDirectCharm;
          return origin;
        }
        if (jLabel < 0){ // safety protection
          AliDebug(1, "Stack label is negative, return\n");
          return -1;
        }

        // if there is an ancester
        TParticle* grandMa = fStack->Particle(jLabel);
        Int_t grandMaPDG = grandMa->GetPdgCode();

        for (Int_t i=0; i<fNparents; i++){
           if (abs(grandMaPDG)==fParentSelect[1][i]){
             origin = kBeautyCharm;
             return origin;
           }
        }

        partMother = grandMa;
     } // end of iteration 
   } // end of if
   else if ( int(abs(maPdgcode)/100.) == kBeauty || int(abs(maPdgcode)/1000.) == kBeauty ) {
     for (Int_t i=0; i<fNparents; i++){
        if (abs(maPdgcode)==fParentSelect[1][i]){
          origin = kDirectBeauty;
          return origin;
        }
     }
   } // end of if
   else if ( abs(maPdgcode) == 22 ) {
     origin = kGamma;
     return origin;
   } // end of if
   else if ( abs(maPdgcode) == 111 ) {
     origin = kPi0;
     return origin;
   } // end of if
   else origin = kElse;

   return origin;
}
