/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
//*************************************************************************
//                   Class AliAnalysisTaskSEBtoJPSItoEle
//                AliAnalysisTaskSE for the reconstruction 
//                   of heavy-flavour decay candidates
//            Author: C.Di Giglio, carmelo.digiglio@ba.infn.it
//*************************************************************************
class AliAnalysisTaskSE;
class AliESDEvent;
class AliVEvent;
class iostream;
class TSystem;
class TROOT;
#include <TFile.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODRecoDecayHF2Prong.h"

#include "AliAnalysisBtoJPSItoEle.h"
#include "AliAnalysisTaskSEBtoJPSItoEle.h"

ClassImp(AliAnalysisTaskSEBtoJPSItoEle)

//______________________________________________________________________________
AliAnalysisTaskSEBtoJPSItoEle::AliAnalysisTaskSEBtoJPSItoEle() :
AliAnalysisTaskSE(),
fOutput(0),
fCdfFit(0),     
fNtupleJPSI(0),
fhDecayTimeMCjpsifromB(0),
fhDecayTime(0),                         
fhInvMass(0),                           
fhD0(0),                                
fhD0D0(0),                              
fhCosThetaStar(0),                      
fhCosThetaPointing(0),                  
fhCtsVsD0D0(0),
fOkAODMC(kFALSE),
fNameMCfile(""),
fOkMinimize(kFALSE)
{
  // default constructor
}                 
//_________________________________________________________________________________________________
AliAnalysisTaskSEBtoJPSItoEle::AliAnalysisTaskSEBtoJPSItoEle(const char* name) :
AliAnalysisTaskSE(name),
fOutput(0),
fCdfFit(0),
fNtupleJPSI(0),
fhDecayTimeMCjpsifromB(0),
fhDecayTime(0),                         
fhInvMass(0),                           
fhD0(0),                                
fhD0D0(0),                              
fhCosThetaStar(0),                      
fhCosThetaPointing(0),                  
fhCtsVsD0D0(0),
fOkAODMC(kFALSE),
fNameMCfile(""),
fOkMinimize(kFALSE)
{
  // standard constructor

  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
}
//_________________________________________________________________________________________________
AliAnalysisTaskSEBtoJPSItoEle::~AliAnalysisTaskSEBtoJPSItoEle()
{
  // destructor

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEBtoJPSItoEle::Init()
{
  // Initialization
  //Double_t ptCuts[2] = {1.,1.5}; // 1 GeV < pt < 1.5 GeV
  Double_t ptCuts[2] = {1.,100}; // 
  SetPtCuts(ptCuts);
  SetMinimize(kTRUE);
  SetAODMCInfo(kTRUE);

  if(fDebug > 1) printf("AliAnalysisTaskSEBtoJPSItoEle::Init() \n");

  return;

}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEBtoJPSItoEle::UserCreateOutputObjects() 
{
  // Create the output container

  if(fDebug > 1) printf("AliAnalysisTaskSEBtoJPSItoEle::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();

  if(fOkAODMC){
  fhDecayTimeMCjpsifromB = new TH1F("fhDecayTimeMCjpsifromB", "Secondary J/#Psi Monte Carlo pseudo proper decay time; X [#mu m]; Entries",100,-4000.,4000.);
  fhDecayTimeMCjpsifromB->Sumw2();
  fhDecayTimeMCjpsifromB->SetMinimum(0);
  fOutput->Add(fhDecayTimeMCjpsifromB);
  }

  fhDecayTime = new TH1F("fhDecayTime", "J/#Psi pseudo proper decay time; X [#mu m]; Entries",100,-4000.,4000.);
  fhDecayTime->Sumw2();
  fhDecayTime->SetMinimum(0);
  fOutput->Add(fhDecayTime);

  fhInvMass = new TH1F("fhInvMass", "J/#Psi invariant mass; M [GeV]; Entries",100,0.,3.2);
  fhInvMass->Sumw2();
  fhInvMass->SetMinimum(0);
  fOutput->Add(fhInvMass);

  fhD0 = new TH1F("fhD0", "Impact parameter; D0 [#mu m]; Entries",100,-5000.,5000.);
  fhD0->Sumw2();
  fhD0->SetMinimum(0);
  fOutput->Add(fhD0);

  fhD0D0 = new TH1F("fhD0D0", "Product of impact parameters; D0D0 [#mu m^2]; Entries",100,-100000.,100000.);
  fhD0D0->Sumw2();
  fhD0D0->SetMinimum(0);
  fOutput->Add(fhD0D0);

  fhCosThetaStar = new TH1F("fhCosThetaStar", "Cosine of decay angle; Cosine Theta Star; Entries",50,-1.2,1.2);
  fhCosThetaStar->Sumw2();
  fhCosThetaStar->SetMinimum(0);
  fOutput->Add(fhCosThetaStar);

  fhCosThetaPointing = new TH1F("fhCosThetaPointing", "Cosine of pointing angle; Cosine Pointing Angle; Entries",100,-1,1);
  fhCosThetaPointing->Sumw2();
  fhCosThetaPointing->SetMinimum(0);
  fOutput->Add(fhCosThetaPointing);

  fhCtsVsD0D0 = new TH2F("fhCtsVsD0D0", "Cosine theta star Vs. D0; Cosine theta star; D0D0",50,-1,1,100,-100000,100000);
  fhCtsVsD0D0->Sumw2();
  fhCtsVsD0D0->SetMinimum(0);
  fOutput->Add(fhCtsVsD0D0);

  fNtupleJPSI = new TNtuple("fNtupleJPSI","J/#Psi pseudo-proper decay time & invariant mass","Xdecaytime:Mass");
  fOutput->Add(fNtupleJPSI);

  return;
}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEBtoJPSItoEle::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  
  // load JPSI candidates                                                   
  TClonesArray *inputArrayJPSItoEle =
    (TClonesArray*)aod->GetList()->FindObject("JPSItoEle");
  if(!inputArrayJPSItoEle) {
    printf("AliAnalysisTaskSECompareHF::UserExec: JPSItoEle branch not found!\n");
    return;
  } 

  // Read AOD Monte Carlo information
  if(fOkAODMC) ReadAODMCInfo(aod,inputArrayJPSItoEle);

  // loop over J/Psi->ee candidates
  Int_t nInJPSItoEle = inputArrayJPSItoEle->GetEntriesFast();
  printf("Number of B->JPSI->e+e-: %d\n",nInJPSItoEle);

  for (Int_t iJPSItoEle = 0; iJPSItoEle < nInJPSItoEle; iJPSItoEle++) {
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)inputArrayJPSItoEle->UncheckedAt(iJPSItoEle);
    //Secondary vertex
    Double_t vtxSec[3] = {0.,0.,0.};
    Double_t vtxPrim[3] = {0.,0.,0.};
    d->GetSecondaryVtx(vtxSec); 
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
    //Int_t okJPSI=0;
    // here analyze only if cuts are passed
    //if(d->SelectBtoJPSI(fCuts,okJPSI)) {
    if(d->Pt() < fPtCuts[1] && d->Pt() > fPtCuts[0]){ // apply pt cut only
      fhInvMass->Fill(d->InvMassJPSIee()); 
      fhD0->Fill(10000*d->ImpParXY());
      fhD0D0->Fill(1e8*d->Prodd0d0());
      fhCosThetaStar->Fill(d->CosThetaStarJPSI());      
      fhCtsVsD0D0->Fill(d->CosThetaStarJPSI(),1e8*d->Prodd0d0());
      fhCosThetaPointing->Fill(d->CosPointingAngle());

      // compute X variable
      AliAODVertex* primVertex = d->GetOwnPrimaryVtx();
      vtxPrim[0] = primVertex->GetX();
      vtxPrim[1] = primVertex->GetY();
      vtxPrim[2] = primVertex->GetZ();
      Double_t lxy = ((vtxSec[0]-vtxPrim[0])*(d->Px()) + (vtxSec[1]-vtxPrim[1])*(d->Py()))/d->Pt();
      Double_t pseudoProperDecayTime = lxy*(TDatabasePDG::Instance()->GetParticle(443)->Mass())/d->Pt();

      fhDecayTime->Fill(10000*pseudoProperDecayTime);
      // Post the data already here
      PostData(1,fOutput);
    
      fNtupleJPSI->Fill(pseudoProperDecayTime,d->InvMassJPSIee()); // fill N-tuple with X,M values

     } // end of JPSItoEle candidates selection according to cuts

    if(unsetvtx) d->UnsetOwnPrimaryVtx();

   }// end loop on JPSI to ele candidates

}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEBtoJPSItoEle::Terminate(Option_t */*option*/)
{
  //
  // Terminate analysis
  //

  if(fDebug > 1) printf("AliAnalysisTaskSEBtoJPSItoEle: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }

  if(fOkAODMC) fhDecayTimeMCjpsifromB = dynamic_cast<TH1F*>(fOutput->FindObject("fhDecayTimeMCjpsifromB"));
  fhDecayTime = dynamic_cast<TH1F*>(fOutput->FindObject("fhDecayTime"));
  fhInvMass = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMass"));
  fhD0 = dynamic_cast<TH1F*>(fOutput->FindObject("fhD0"));
  fhD0D0 = dynamic_cast<TH1F*>(fOutput->FindObject("fhD0D0"));
  fhCosThetaStar = dynamic_cast<TH1F*>(fOutput->FindObject("fhCosThetaStar"));
  fhCosThetaPointing = dynamic_cast<TH1F*>(fOutput->FindObject("fhCosThetaPointing"));
  fhCtsVsD0D0 = dynamic_cast<TH2F*>(fOutput->FindObject("fhCtsVsD0D0"));
  fNtupleJPSI = dynamic_cast<TNtuple*>(fOutput->FindObject("fNtupleJPSI"));

  if(fOkMinimize){
     Double_t* pseudoProperValues=0x0;
     Double_t* invMassValues=0x0;
     Int_t ncand=0;
     fCdfFit = new AliAnalysisBtoJPSItoEle();

     printf("+++\n+++  AliAnalysisBtoJPSItoEle object instantiated ---> OK\n+++\n");
     fCdfFit->ReadCandidates(fNtupleJPSI,pseudoProperValues,invMassValues,ncand);
     printf("+++\n+++  X and M vectors filled starting from N-tuple ---> OK\n+++\n");
     printf("+++\n+++  Number of candidates ---> %d J/psi \n+++\n",ncand);
     fCdfFit->SetPtBin(0);
     printf("+++\n+++  Pt bin setted n. ---> %d  \n+++\n",fCdfFit->GetPtBin());

     if(fOkAODMC) {fCdfFit->CloneMCtemplate(fhDecayTimeMCjpsifromB);} 
        else { 
          //SetPathMCfile(".");
          TH1F* hLocal = OpenLocalMCFile(fNameMCfile,fCdfFit->GetPtBin()); 
          fCdfFit->CloneMCtemplate(hLocal); 
        }

     printf("+++\n+++  MC template histo copied ---> OK\n+++\n");
     fCdfFit->DoMinimization(pseudoProperValues,invMassValues,ncand);
   } 

  return;
}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEBtoJPSItoEle::SetPtCuts(const Double_t ptCuts[2])
{
  //
  // apply Pt cuts (lower cut, upper cut)
  //
  for(Int_t i = 0; i < 2; i++) fPtCuts[i] = ptCuts[i];
}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEBtoJPSItoEle::SetCutsD0(const Double_t cuts[9])
{
  //
  // apply D0 and JPSI cuts 
  //
  for(Int_t i = 0; i < 9; i++) fCuts[i] = cuts[i];
}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEBtoJPSItoEle::ReadAODMCInfo(const AliAODEvent* aodEvent, const TClonesArray* inputArray) 
{
  //
  // Reads MC information if needed (for example in case of parametrized PID)
  //

  // load MC particles
  TClonesArray *mcArray =
    (TClonesArray*)aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(!mcArray) {
    printf("AliAnalysisTaskSEBtoJPSItoEle::UserExec: MC particles branch not found!\n");
    return;
  }

  // load MC header 
  AliAODMCHeader* mcHeader =
    (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader){
    printf("AliAnalysisTaskSEBtoJPSItoEle::UserExec: MC AOD header branch not found!\n");
  }

  Double_t rmax = 0.000005;
  Double_t mcPrimVtx[3];

  // get MC primary Vtx
  mcHeader->GetVertex(mcPrimVtx);

  Int_t nInCandidates = inputArray->GetEntriesFast();
  printf("Number of Candidates for MC analysis: %d\n",nInCandidates);

  Int_t lab0, lab1, pdgMother, labMother, pdgJPSI, pdg0, pdg1, labJPSIMother, pdgJPSIMother;
  Int_t labJPSIdaugh0=0;
  Int_t labJPSIdaugh1=0;

  for (Int_t iCandidate = 0; iCandidate < nInCandidates; iCandidate++) {
    AliAODRecoDecayHF2Prong *dd = (AliAODRecoDecayHF2Prong*)inputArray->UncheckedAt(iCandidate);
    if(dd->Pt() < fPtCuts[1] && dd->Pt() > fPtCuts[0]){ // apply pt cut only

      AliAODTrack *trk0 = (AliAODTrack*)dd->GetDaughter(0);
      AliAODTrack *trk1 = (AliAODTrack*)dd->GetDaughter(1);
      lab0 = trk0->GetLabel();
      lab1 = trk1->GetLabel();

      AliAODMCParticle* part0 = (AliAODMCParticle*)mcArray->At(lab0);
      if(!part0) {
        printf("no MC particle\n");
        continue;
      }

      while(part0->GetMother()>=0) {
       labMother=part0->GetMother();
       part0 = (AliAODMCParticle*)mcArray->At(labMother);
       if(!part0) {
         printf("no MC mother particle\n");
         break;
       }
       pdgMother = TMath::Abs(part0->GetPdgCode());
       if(pdgMother==443) {//this for JPSI
         labJPSIdaugh0=labMother;
         break;
       }
      }

      AliAODMCParticle* part1 = (AliAODMCParticle*)mcArray->At(lab1);
      if(!part1) {
        printf("no MC particle\n");
        continue;
      }

      while(part1->GetMother()>=0) {
       labMother=part1->GetMother();
       part1 = (AliAODMCParticle*)mcArray->At(labMother);
       if(!part1) {
         printf("no MC mother particle\n");
         break;
       }
       pdgMother = TMath::Abs(part1->GetPdgCode());
       if(pdgMother==443) {//this for JPSI
         labJPSIdaugh1=labMother;
         break;
       }
      }

      if(labJPSIdaugh0>=0 && labJPSIdaugh1>=0 && labJPSIdaugh0==labJPSIdaugh1) { // check if JPSI is signal
        AliAODMCParticle* partJPSI = (AliAODMCParticle*)mcArray->At(labJPSIdaugh0);
        pdgJPSI = partJPSI->GetPdgCode();
        if(pdgJPSI == 443){//this is for MC JPSI
           if(TMath::Sqrt((partJPSI->Xv()-mcPrimVtx[0])*(partJPSI->Xv()-mcPrimVtx[0])+
                          (partJPSI->Yv()-mcPrimVtx[1])*(partJPSI->Yv()-mcPrimVtx[1])+
                          (partJPSI->Zv()-mcPrimVtx[2])*(partJPSI->Zv()-mcPrimVtx[2])>rmax)){ //this is for MC displaced JPSI
              Int_t pdaugh0 = partJPSI->GetDaughter(0);
              Int_t pdaugh1 = partJPSI->GetDaughter(1);
              AliAODMCParticle* partDaugh0 = (AliAODMCParticle*)mcArray->At(pdaugh0);
              AliAODMCParticle* partDaugh1 = (AliAODMCParticle*)mcArray->At(pdaugh1);
              pdg0 = partDaugh0->GetPdgCode();
              pdg1 = partDaugh1->GetPdgCode();
              if(TMath::Abs(pdg0) == 11 && TMath::Abs(pdg1) == 11){ // this is for MC JPSI -> ee
                labJPSIMother = partJPSI->GetMother();
                AliAODMCParticle* partJPSIMother = (AliAODMCParticle*)mcArray->At(labJPSIMother);
                pdgJPSIMother = partJPSIMother->GetPdgCode();
               if(pdgJPSIMother==511   || pdgJPSIMother==521   ||
                  pdgJPSIMother==10511 || pdgJPSIMother==10521 ||
                  pdgJPSIMother==513   || pdgJPSIMother==523   ||
                  pdgJPSIMother==10513 || pdgJPSIMother==10523 ||
                  pdgJPSIMother==20513 || pdgJPSIMother==20523 ||
                  pdgJPSIMother==515   || pdgJPSIMother==525   ||
                  pdgJPSIMother==531   || pdgJPSIMother==10531 ||
                  pdgJPSIMother==533   || pdgJPSIMother==10533 ||
                  pdgJPSIMother==20533 || pdgJPSIMother==535   ||
                  pdgJPSIMother==541   || pdgJPSIMother==10541 ||
                  pdgJPSIMother==543   || pdgJPSIMother==10543 ||
                  pdgJPSIMother==20543 || pdgJPSIMother==545)
                  { //this is for MC JPSI -> ee from B-hadrons

                  Double_t mcLxy = ((partJPSI->Xv()-mcPrimVtx[0])*(partJPSI->Px()) + (partJPSI->Yv()-mcPrimVtx[1])*(partJPSI->Py()))/partJPSI->Pt();
                  Double_t mcPseudoProperDecayTime = mcLxy*(TDatabasePDG::Instance()->GetParticle(443)->Mass())/partJPSI->Pt();
                  fhDecayTimeMCjpsifromB->Fill(10000*mcPseudoProperDecayTime,1);

                  // Post the data already here
                  PostData(1,fOutput);

                } //this is for MC JPSI -> ee from B-hadrons
              } //this is for MC JPSI -> ee
            } //this is for MC displaced JPSI
          } //this is for MC JPSI 
        } // end of check if JPSI is signal
      }
    }

}
//_________________________________________________________________________________________________
TH1F* AliAnalysisTaskSEBtoJPSItoEle::OpenLocalMCFile(TString datadir, Int_t binNum) 
{
  //
  // Open a local file with MC x distribution for JPSI from B-hadron
  //

  TH1F* hCsiMC = new TH1F();
  TString useFile = datadir.Data();
  useFile.Append("/CsiMCfromKine_10PtBins.root");
  TFile *f = new TFile(useFile);
  Double_t scale = 0.;
  char processCsiMCXHisto[1024];
  if(binNum == 0) sprintf(processCsiMCXHisto,"hPseudoProper");
  else sprintf(processCsiMCXHisto,"hPseudoProper%d",binNum);
  hCsiMC = (TH1F*)f->Get(processCsiMCXHisto);
  scale=1/hCsiMC->Integral();
  hCsiMC->Scale(scale);
  printf ("Opening Histo with Csi_MC(x) template with n. %f entries", hCsiMC->GetEntries());

  return hCsiMC;

}

