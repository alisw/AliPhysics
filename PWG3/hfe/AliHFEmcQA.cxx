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
/**************************************************************************
 *                                                                        *
 * QA class of Heavy Flavor quark and fragmeted/decayed particles         *
 *                                                                        *
 * Authors:                                                               *
 *   MinJung Kweon <minjung@physi.uni-heidelberg.de>                      *
 *                                                                        *
 **************************************************************************/


#include <iostream>
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCut.h>
#include <TRandom.h>
#include <TDatabasePDG.h>
#include <TROOT.h>
#include <TParticle.h>

#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliESDInputHandler.h>
#include <AliESDVertex.h>
#include <AliStack.h>

#include "AliHFEmcQA.h"

const Int_t AliHFEmcQA::fgkGluon=21;
const Int_t AliHFEmcQA::fgkMaxGener=10;
const Int_t AliHFEmcQA::fgkMaxIter=1000;
const Int_t AliHFEmcQA::fgkqType = 6;    // number of species waiting for QA done

ClassImp(AliHFEmcQA)

//_______________________________________________________________________________________________
AliHFEmcQA::AliHFEmcQA() : 
        fVerbos(kFALSE) 
        ,fStack(0x0) 
        ,fNparents(0) 
{
        // Default constructor
        
        if (fVerbos) AliInfo("***** Warning! fVerbos is set to TRUE! ******");
}

//_______________________________________________________________________________________________
AliHFEmcQA::AliHFEmcQA(const AliHFEmcQA&p):
        TObject(p)
        ,fVerbos(p.fVerbos)
        ,fStack(0x0) 
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
void AliHFEmcQA::PostAnalyze()
{
}

//__________________________________________
void AliHFEmcQA::CreateHistograms(const Int_t kquark, TString hnopt) 
{
  // create histograms

  if (kquark != fkCharm && kquark != fkBeauty) {
    AliDebug(1, "This task is only for heavy quark QA, return\n");
    return; 
  }
  Int_t iq = kquark - fkCharm; 

  TString kqTypeLabel[fgkqType];
  if (kquark == fkCharm){
    kqTypeLabel[fkQuark]="c";
    kqTypeLabel[fkantiQuark]="cbar";
    kqTypeLabel[fkElectron]="ce";
    kqTypeLabel[fkElectron2nd]="nulle";
    kqTypeLabel[fkeHadron]="ceHadron";
    kqTypeLabel[fkDeHadron]="nullHadron";
  } else if (kquark == fkBeauty){
    kqTypeLabel[fkQuark]="b";
    kqTypeLabel[fkantiQuark]="bbar";
    kqTypeLabel[fkElectron]="be";
    kqTypeLabel[fkElectron2nd]="bce";
    kqTypeLabel[fkeHadron]="beHadron";
    kqTypeLabel[fkDeHadron]="bDeHadron";
  }


  TString hname; 
  for (Int_t iqType = 0; iqType < fgkqType; iqType++ ){
     hname = hnopt+"PdgCode_"+kqTypeLabel[iqType];
     fHist[iq][iqType].fPdgCode = new TH1F(hname,hname,20001,-10000.5,10000.5);
     hname = hnopt+"Pt_"+kqTypeLabel[iqType];
     fHist[iq][iqType].fPt = new TH1F(hname,hname+";p_{T} (GeV/c)",150,0,30);
     hname = hnopt+"Y_"+kqTypeLabel[iqType];
     fHist[iq][iqType].fY = new TH1F(hname,hname,150,-7.5,7.5);
     hname = hnopt+"Eta_"+kqTypeLabel[iqType];
     fHist[iq][iqType].fEta = new TH1F(hname,hname,150,-7.5,7.5);
  }

  hname = hnopt+"Nq_"+kqTypeLabel[fkQuark];
  fHistComm[iq].fNq = new TH1F(hname,hname,10,-0.5,9.5);
  hname = hnopt+"ProcessID_"+kqTypeLabel[fkQuark];
  fHistComm[iq].fProcessID = new TH1F(hname,hname,21,-10.5,10.5);
  hname = hnopt+"ePtRatio_"+kqTypeLabel[fkQuark];
  fHistComm[iq].fePtRatio = new TH2F(hname,hname+";p_{T} (GeV/c);momentum fraction",100,0,30,100,0,1);
  hname = hnopt+"DePtRatio_"+kqTypeLabel[fkQuark];
  fHistComm[iq].fDePtRatio = new TH2F(hname,hname+";p_{T} (GeV/c);momentum fraction",100,0,30,100,0,1);
  hname = hnopt+"eDistance_"+kqTypeLabel[fkQuark];
  fHistComm[iq].feDistance= new TH2F(hname,hname+";p_{T} (GeV/c);distance (cm)",100,0,30,200,0,2);
  hname = hnopt+"DeDistance_"+kqTypeLabel[fkQuark];
  fHistComm[iq].fDeDistance= new TH2F(hname,hname+";p_{T} (GeV/c);distance (cm)",100,0,30,200,0,2);

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
void AliHFEmcQA::GetQuarkKine(Int_t iTrack, const Int_t kquark) 
{
  // get heavy quark kinematics

    if (kquark != fkCharm && kquark != fkBeauty) {
      AliDebug(1, "This task is only for heavy quark QA, return\n");
      return; 
    }
    Int_t iq = kquark - fkCharm; 


    if (iTrack < 0) { 
      AliDebug(1, "Stack label is negative, return\n");
      return; 
    }

    TParticle *part = fStack->Particle(iTrack); 
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
          Bool_t IsSameString = kTRUE; 
          for (Int_t i=1; i<fgkMaxIter; i++){
             iLabel = iLabel - 1;
             if (iLabel>-1) { partMother = fStack->Particle(iLabel); }
             else {
               AliDebug(1, "Stack label is negative, return\n");
               return; 
             }
             if ( abs(partMother->GetPdgCode()) == kquark ) break;
             if ( partMother->GetStatusCode() != 12 ) IsSameString = kFALSE;
             if (!IsSameString) return; 
          }
        }
        AliDebug(1, "Can not find heavy parton of this heavy hadron in the string, return\n");
        if (abs(partMother->GetPdgCode()) != kquark) return; 

        if (fIsHeavy[iq] >= 50) return;  
        fHeavyQuark[fIsHeavy[iq]] = partMother;
        fIsHeavy[iq]++;

        // fill kinematics for heavy parton
        if (partMother->GetPdgCode() > 0) { // quark
          fHist[iq][fkQuark].fPdgCode->Fill(partMother->GetPdgCode());
          fHist[iq][fkQuark].fPt->Fill(partMother->Pt());
          fHist[iq][fkQuark].fY->Fill(GetRapidity(partMother));
          fHist[iq][fkQuark].fEta->Fill(partMother->Eta());
        } else{ // antiquark
          fHist[iq][fkantiQuark].fPdgCode->Fill(partMother->GetPdgCode());
          fHist[iq][fkantiQuark].fPt->Fill(partMother->Pt());
          fHist[iq][fkantiQuark].fY->Fill(GetRapidity(partMother));
          fHist[iq][fkantiQuark].fEta->Fill(partMother->Eta());
        }

      } // end of heavy parton slection loop 

    } // end of heavy hadron or quark selection

}

//__________________________________________
void AliHFEmcQA::EndOfEventAna(const Int_t kquark)
{
  // end of event analysis

  if (kquark != fkCharm && kquark != fkBeauty) {
    AliDebug(1, "This task is only for heavy quark QA, return\n");
    return; 
  }
  Int_t iq = kquark - fkCharm; 


  // # of heavy quark per event
  AliDebug(1,Form("Number of heavy quark in this event = %d \n",fIsHeavy[iq]));
  fHistComm[iq].fNq->Fill(fIsHeavy[iq]);

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
          reportStrangeness(motherID[i],motherType[i],motherLabel[i]);
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
       if (motherLabel[id1] == motherLabel[id2]) processID = fkGluonSplitting; // gluon spliting
       else processID = -9;
     }
     else if (motherType[id1] == -1 && motherType[id2] == -1) {
       if (motherLabel[id1] == -1 && motherLabel[id2] == -1) {
         if (motherID[id1] == fgkGluon) processID = fkPairCreationFromg; // gluon fusion
         else processID = fkPairCreationFromq; // q-qbar pair creation
       }
       else processID = -8;
     }
     else if (motherType[id1] == -1 || motherType[id2] == -1) {
       if ((motherLabel[id1] == -1 || motherLabel[id2] == -1) && (motherLabel[id1]*motherLabel[id2] == -2 || motherLabel[id1]*motherLabel[id2] == -3)) {
         if(motherID[id1]*motherID[id2] == kquark*fgkGluon) processID = fkFlavourExitation; // flavour exitation 
         else processID = fkLightQuarkShower;
       }
       else processID = -7;
     }
     else if (motherType[id1] == -2 || motherType[id2] == -2) {
       if (motherLabel[id1] == motherLabel[id2]) processID = fkInitialPartonShower; // initial parton shower
       else processID = -6;
       
     }
     else processID = -5;

     if (nheavypair >1) AliDebug(1,Form("Multi pair found : process ID = %d\n",processID));
     else fHistComm[iq].fProcessID->Fill(processID);
     AliDebug(1,Form("Process ID = %d\n",processID));
  } // end of # heavy quark pair loop

}

//__________________________________________
void AliHFEmcQA::GetDecayedKine(Int_t iTrack, const Int_t kquark, Int_t kdecayed, Bool_t isbarrel) 
{
    // decay electron kinematics
    
    if (kquark != fkCharm && kquark != fkBeauty) {
      AliDebug(1, "This task is only for heavy quark QA, return\n");
      return; 
    }
    Int_t iq = kquark - fkCharm; 

    if (iTrack < 0) { 
      AliDebug(1, "Stack label is negative, return\n");
      return; 
    }

    TParticle* mcpart = fStack->Particle(iTrack);

    if ( abs(mcpart->GetPdgCode()) != kdecayed ) return;
    if ( isbarrel && TMath::Abs(mcpart->Eta()) > 0.9 ) return;

    Int_t iLabel = mcpart->GetFirstMother(); 
    if (iLabel<0){
      AliDebug(1, "Stack label is negative, return\n");
      return; 
    }

    TParticle *partMother = fStack->Particle(iLabel);
    Int_t maPdgcode = partMother->GetPdgCode();
    TParticle *partMotherCopy = fStack->Particle(iLabel);
    Int_t maPdgcodeCopy = partMotherCopy->GetPdgCode();

    // get electron production vertex   
    TLorentzVector ePoint;
    mcpart->ProductionVertex(ePoint);


    // if the mother is charmed hadron  
    Bool_t IsMotherDirectCharm = kFALSE;
    if ( int(abs(maPdgcode)/100.) == fkCharm || int(abs(maPdgcode)/1000.) == fkCharm ) { 

          // iterate until you find B hadron as a mother or become top ancester 
          for (Int_t i=1; i<fgkMaxIter; i++){

             Int_t jLabel = partMother->GetFirstMother(); 
             if (jLabel == -1){
               IsMotherDirectCharm = kTRUE;
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

                  if (kquark == fkCharm) return;
                  // fill electron kinematics
                  fHist[iq][fkElectron2nd].fPdgCode->Fill(mcpart->GetPdgCode());
                  fHist[iq][fkElectron2nd].fPt->Fill(mcpart->Pt());
                  fHist[iq][fkElectron2nd].fY->Fill(GetRapidity(mcpart));
                  fHist[iq][fkElectron2nd].fEta->Fill(mcpart->Eta());

                  // fill mother hadron kinematics
                  fHist[iq][fkDeHadron].fPdgCode->Fill(grandMaPDG); 
                  fHist[iq][fkDeHadron].fPt->Fill(grandMa->Pt());
                  fHist[iq][fkDeHadron].fY->Fill(GetRapidity(grandMa));
                  fHist[iq][fkDeHadron].fEta->Fill(grandMa->Eta());
                 
                  // ratio between pT of electron and pT of mother B hadron 
                  if(grandMa->Pt()) fHistComm[iq].fDePtRatio->Fill(grandMa->Pt(),mcpart->Pt()/grandMa->Pt());

                  // distance between electron production point and mother hadron production point
                  TLorentzVector debPoint;
                  grandMa->ProductionVertex(debPoint);
                  TLorentzVector dedistance = ePoint - debPoint;
                  fHistComm[iq].fDeDistance->Fill(grandMa->Pt(),dedistance.P());
                  return;
                }
             }

             partMother = grandMa;
          } // end of iteration 
    } // end of if
    if((IsMotherDirectCharm == kTRUE && kquark == fkCharm) || kquark == fkBeauty) {
         for (Int_t i=0; i<fNparents; i++){
            if (abs(maPdgcodeCopy)==fParentSelect[iq][i]){

              // fill electron kinematics
              fHist[iq][fkElectron].fPdgCode->Fill(mcpart->GetPdgCode());
              fHist[iq][fkElectron].fPt->Fill(mcpart->Pt());
              fHist[iq][fkElectron].fY->Fill(GetRapidity(mcpart));
              fHist[iq][fkElectron].fEta->Fill(mcpart->Eta());  

              // fill mother hadron kinematics
              fHist[iq][fkeHadron].fPdgCode->Fill(maPdgcodeCopy); 
              fHist[iq][fkeHadron].fPt->Fill(partMotherCopy->Pt());
              fHist[iq][fkeHadron].fY->Fill(GetRapidity(partMotherCopy));
              fHist[iq][fkeHadron].fEta->Fill(partMotherCopy->Eta());

              // ratio between pT of electron and pT of mother B hadron 
              if(partMotherCopy->Pt()) fHistComm[iq].fePtRatio->Fill(partMotherCopy->Pt(),mcpart->Pt()/partMotherCopy->Pt());

              // distance between electron production point and mother hadron production point
              TLorentzVector ebPoint;
              partMotherCopy->ProductionVertex(ebPoint);
              TLorentzVector edistance = ePoint - ebPoint;
              fHistComm[iq].feDistance->Fill(partMotherCopy->Pt(),edistance.P());
            }
         }
    } // end of if
}


//__________________________________________
void AliHFEmcQA::IdentifyMother(Int_t mother_label, Int_t &mother_pdg, Int_t &grandmother_label)
{
       if (mother_label < 0) { 
         AliDebug(1, "Stack label is negative, return\n");
         return; 
       }
       TParticle *heavysMother = fStack->Particle(mother_label);
       mother_pdg = heavysMother->GetPdgCode();
       grandmother_label = heavysMother->GetFirstMother();
       AliDebug(1,Form("ancestor pdg code= %d\n",mother_pdg));
}

//__________________________________________
// mothertype -1 means this heavy quark coming from hard vertex
void AliHFEmcQA::HardScattering(const Int_t kquark, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
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
// mothertype -2 means this heavy quark coming from initial state 
Bool_t AliHFEmcQA::IsFromInitialShower(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
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
// mothertype 2 means this heavy quark coming from final state 
Bool_t AliHFEmcQA::IsFromFinalParton(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
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
void AliHFEmcQA::reportStrangeness(Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
       mothertype = -888;
       motherlabel = -888;
       motherID = -888;
       AliDebug(1,"something strange!\n");
}

//__________________________________________
Float_t AliHFEmcQA::GetRapidity(TParticle *part)
{
       Float_t rapidity;        
       if(part->Energy() - part->Pz() == 0 || part->Energy() + part->Pz() == 0) rapidity=-999; 
       else rapidity = 0.5*(TMath::Log((part->Energy()+part->Pz()) / (part->Energy()-part->Pz()))); 
       return rapidity;
}
