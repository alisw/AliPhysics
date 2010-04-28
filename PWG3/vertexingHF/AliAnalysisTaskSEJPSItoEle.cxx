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
//                    Class AliAnalysisTaskSEJPSItoEle
// AliAnalysisTaskSE class to read both J/psi -> e+e- candidates and 
// like-sign pair candidates and to store them into a specific
// stand-alone AOD file for J/psi into dieletrons analysis. 
// The current task:
//
//  1) produces several histograms (including invariant mass distributions)
//     for both unlike sign (US) and like sign (LS) pairs.
//
//  2) selects only J/Psi to dielectrons candidates and copies it to a 
//     specific AOD file, namely AliAODjpsi.root. The references
//     to AliAODTrack objects are copied as well. 
//     The specific AliAODjpsi.root is thought  as the input file 
//     to the AliAnalysisBtoJPSItoEle.h class, which performs (locally)
//     the analysis of J/Psi from beauty hadrons decays. 
//     
//     Author: C.Di Giglio, carmelo.digiglio@ba.infn.it
//*************************************************************************
class AliAnalysisTaskSE;
class AliESDEvent;
class AliVEvent;
class iostream;
class TSystem;
class TROOT;
#include "TFile.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TROOT.h"
#include "TCanvas.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSEJPSItoEle.h"
#include "AliAODTrack.h"
#include "AliExternalTrackParam.h"

ClassImp(AliAnalysisTaskSEJPSItoEle)

//______________________________________________________________________________
AliAnalysisTaskSEJPSItoEle::AliAnalysisTaskSEJPSItoEle() :
AliAnalysisTaskSE(),
fOutput(0),
fhDecayTimeMCjpsifromB(0),
fhDecayTime(0),                         
fhDecayTimeOut(0),
fhInvMass(0),                           
fhD0(0),                                
fhD0D0(0),                              
fhCosThetaStar(0),                      
fhCosThetaPointing(0),                  
fhDCA(0),
fhCtsVsD0D0(0),
fHistMassLS(0),
fHistCtsLS(0),
fHistCtsLSpos(0),
fHistCtsLSneg(0),
fHistCPtaLS(0),
fHistd0d0LS(0),
fHistDCALS(0),
fVHF(0),
fTotPosPairs(0),
fTotNegPairs(0),
fLsNormalization(1.),
fOkAODMC(kFALSE),
fOkLikeSign(kFALSE),
fVerticesHFTClArr(0),
fJpsiToEleTClArr(0),
fLikeSignTClArr(0),
fTracksTClArr(0),
fChain(0),
fOrigAOD(0),
fNewAOD(0)
{
  // default constructor
}                 
//_________________________________________________________________________________________________
AliAnalysisTaskSEJPSItoEle::AliAnalysisTaskSEJPSItoEle(const char* name) :
AliAnalysisTaskSE(name),
fOutput(0),
fhDecayTimeMCjpsifromB(0),
fhDecayTime(0),
fhDecayTimeOut(0),
fhInvMass(0),                           
fhD0(0),                                
fhD0D0(0),                              
fhCosThetaStar(0),                      
fhCosThetaPointing(0),                  
fhDCA(0),
fhCtsVsD0D0(0),
fHistMassLS(0),
fHistCtsLS(0),
fHistCtsLSpos(0),
fHistCtsLSneg(0),
fHistCPtaLS(0),
fHistd0d0LS(0),
fHistDCALS(0),
fVHF(0),
fTotPosPairs(0),
fTotNegPairs(0),
fLsNormalization(1.),
fOkAODMC(kFALSE),
fOkLikeSign(kFALSE),
fVerticesHFTClArr(0),
fJpsiToEleTClArr(0),
fLikeSignTClArr(0),
fTracksTClArr(0),
fChain(0),
fOrigAOD(0),
fNewAOD(0)
{
  // standard constructor

  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());

  // Output slot #0 writes into a TTree container
  // Output slot #1 writes into a TList container
  DefineOutput(0, TTree::Class());
  DefineOutput(1,TList::Class());  //My private output
}
//_________________________________________________________________________________________________
AliAnalysisTaskSEJPSItoEle::~AliAnalysisTaskSEJPSItoEle()
{
  // destructor

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  } 

  if (fVHF) {
    delete fVHF;
    fVHF = 0;
  }

}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEJPSItoEle::Init()
{
  // Initialization

  if(fDebug > 1) printf("AliAnalysisTaskSEJPSItoEle::Init() \n");

  gROOT->LoadMacro("ConfigVertexingHF.C");

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");
  fVHF->PrintStatus();

  return;
}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEJPSItoEle::UserCreateOutputObjects() 
{
  // Create the output container

  if(fDebug > 1) printf("AliAnalysisTaskSEJPSItoEle::UserCreateOutputObjects() \n");

  if (!AODEvent()) {
     Fatal("UserCreateOutputObjects", "This task needs an AOD handler");
     return;
  }   

  if(!fNewAOD) fNewAOD = new AliAODEvent();
  fNewAOD->CreateStdContent();

  TString filename = "AliAOD.Jpsitoele.root";
  if (!IsStandardAOD()) filename = "";

  // Array of HF vertices  
  fVerticesHFTClArr = new TClonesArray("AliAODVertex", 0);
  fVerticesHFTClArr->SetName("VerticesHF");
  AddAODBranch("TClonesArray", &fVerticesHFTClArr);

  // Array of J/psi -> e+e- candidates 
  fJpsiToEleTClArr = new TClonesArray("AliAODRecoDecayHF2Prong",0);
  fJpsiToEleTClArr->SetName("JpsiToEleEvents");
  AddAODBranch("TClonesArray", &fJpsiToEleTClArr, filename);

  // Array of like sign candidates 
  fLikeSignTClArr = new TClonesArray("AliAODRecoDecayHF2Prong",0);
  fLikeSignTClArr->SetName("LikeSignEvents");
  AddAODBranch("TClonesArray", &fLikeSignTClArr, filename);

  // Array of tracks from J/psi -> e+e- candidates
  fTracksTClArr = new TClonesArray("AliAODTrack", 0);
  fTracksTClArr->SetName("Tracks");
  AddAODBranch("TClonesArray", &fTracksTClArr, filename);

  fOutput = new TList();
  fOutput->SetOwner();

  if(fOkAODMC){
  fhDecayTimeMCjpsifromB = new TH1F("fhDecayTimeMCjpsifromB", "Secondary J/#Psi Monte Carlo pseudo proper decay time; X [#mu m]; Entries",100,-4000.,4000.);
  fhDecayTimeMCjpsifromB->Sumw2();
  fhDecayTimeMCjpsifromB->SetMinimum(0);
  fOutput->Add(fhDecayTimeMCjpsifromB);
  }

  // Invariant mass
  fhInvMass = new TH1F("fhInvMass", "J/#Psi invariant mass; M [GeV]; Entries",100,0.,3.2);
  fhInvMass->Sumw2();
  fhInvMass->SetMinimum(0);
  fOutput->Add(fhInvMass);

  fHistMassLS = new TH1F("fHistMassLS", "Like sign pairs invariant mass; M [GeV]; Entries",200,2.8,3.25);
  fHistMassLS->Sumw2();
  fHistMassLS->SetMinimum(0);
  fOutput->Add(fHistMassLS);

  // Pseudo proper decay time
  fhDecayTime = new TH1F("fhDecayTime", "J/#Psi pseudo proper decay time; X [#mu m]; Entries",200,-2000.,2000.);
  fhDecayTime->Sumw2();
  fhDecayTime->SetMinimum(0);
  fOutput->Add(fhDecayTime);
 
  // Pseudo proper decay time
  fhDecayTimeOut = new TH1F("fhDecayTimeOut", "J/#Psi pseudo proper decay time (output standalone AOD); X [#mu m]; Entries",200,-2000.,2000.);
  fhDecayTimeOut->Sumw2();
  fhDecayTimeOut->SetMinimum(0);
  fOutput->Add(fhDecayTimeOut);

  // Dictance of closest approach
  fhDCA = new TH1F("fhDCA", "J/#Psi distance of closest approach; dca [10^{2}#mu m]; Entries",100,0.,5.);
  fhDCA->Sumw2();
  fhDCA->SetMinimum(0);
  fOutput->Add(fhDCA);

  fHistDCALS = new TH1F("fHistDCALS", "Like sign pairs distance of closest approach; dca [10^{2}#mu m]; Entries",100,0.,5.);
  fHistDCALS->Sumw2();
  fHistDCALS->SetMinimum(0);
  fOutput->Add(fHistDCALS);

  // Impact parameter
  fhD0 = new TH1F("fhD0", "Impact parameter; d0 [#mu m]; Entries",100,-5000.,5000.);
  fhD0->Sumw2();
  fhD0->SetMinimum(0);
  fOutput->Add(fhD0);

  // Product of impact parameters
  fhD0D0 = new TH1F("fhD0D0", "Product of impact parameters; D0D0 [#mu m^2]; Entries",100,-100000.,100000.);
  fhD0D0->Sumw2();
  fhD0D0->SetMinimum(0);
  fOutput->Add(fhD0D0);

  fHistd0d0LS = new TH1F("fHistd0d0LS", "Like sign pairs product of impact parameters; d0xd0 [#mu m^{2}]; Entries",200,-100000.,100000.);
  fHistd0d0LS->Sumw2();
  fHistd0d0LS->SetMinimum(0);
  fOutput->Add(fHistd0d0LS);

  // Cosine of the decay angle
  fhCosThetaStar = new TH1F("fhCosThetaStar", "Cosine of decay angle; Cosine Theta Star; Entries",50,-1.2,1.2);
  fhCosThetaStar->Sumw2();
  fhCosThetaStar->SetMinimum(0);
  fOutput->Add(fhCosThetaStar);

  fHistCtsLS = new TH1F("fHistCtsLS", "Like sign pairs cosine of decay angle; Cos#Theta*; Entries",200,-1.,1.);
  fHistCtsLS->Sumw2();
  fHistCtsLS->SetMinimum(0);
  fOutput->Add(fHistCtsLS);

  fHistCtsLSpos = new TH1F("fHistCtsLSpos", "Like sign ++ pairs cosine of decay angle; Cos#Theta*; Entries",200,-1.,1.);
  fHistCtsLSpos->Sumw2();
  fHistCtsLSpos->SetMinimum(0);
  fOutput->Add(fHistCtsLSpos);

  fHistCtsLSneg = new TH1F("fHistCtsLSneg", "Like sign -- pairs cosine of decay angle; Cos#Theta*; Entries",200,-1.,1.);
  fHistCtsLSneg->Sumw2();
  fHistCtsLSneg->SetMinimum(0);
  fOutput->Add(fHistCtsLSneg);

  // Cosine of pointing angle
  fhCosThetaPointing = new TH1F("fhCosThetaPointing", "Cosine of pointing angle; Cosine Pointing Angle; Entries",100,-1,1);
  fhCosThetaPointing->Sumw2();
  fhCosThetaPointing->SetMinimum(0);
  fOutput->Add(fhCosThetaPointing);

  fHistCPtaLS = new TH1F("fHistCPtaLS", "Like sign pairs cosine of pointing angle; Cos#Theta_{point}; Entries",200,-1.,1.);
  fHistCPtaLS->Sumw2();
  fHistCPtaLS->SetMinimum(0);
  fOutput->Add(fHistCPtaLS);

  // CosThetaStar vs. d0xd0
  fhCtsVsD0D0 = new TH2F("fhCtsVsD0D0", "Cosine theta star Vs. D0; Cosine theta star; D0D0",50,-1,1,100,-100000,100000);
  fhCtsVsD0D0->Sumw2();
  fhCtsVsD0D0->SetMinimum(0);
  fOutput->Add(fhCtsVsD0D0);

  return;
}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEJPSItoEle::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  TClonesArray *arrayJPSItoEle = 0;
  TClonesArray *arrayLikeSign = 0;
  TClonesArray *arrayTracks = 0;

  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();

      // load tracks from AOD event   
      arrayTracks=(TClonesArray*)aod->GetList()->FindObject("tracks");
      //Int_t totTracks = arrayTracks->GetEntriesFast();

      // load Jpsi candidates from AOD friend  
      arrayJPSItoEle=(TClonesArray*)aodFromExt->GetList()->FindObject("JPSItoEle");
      //Int_t totJPSItoEleCand = arrayJPSItoEle->GetEntriesFast();

      // load like sign candidates from AOD friend
      arrayLikeSign=(TClonesArray*)aodFromExt->GetList()->FindObject("LikeSign2Prong");
      //Int_t totLikeSignCand = arrayLikeSign->GetEntriesFast();

    }

  } else {
    // load Jpsi candidates                                                   
    arrayJPSItoEle=(TClonesArray*)aod->GetList()->FindObject("JPSItoEle");
    //Int_t totJPSItoEleCand = arrayJPSItoEle->GetEntriesFast();
    // load like sign candidates
    arrayLikeSign=(TClonesArray*)aod->GetList()->FindObject("LikeSign2Prong");
    //Int_t totLikeSignCand = arrayLikeSign->GetEntriesFast();
  }

  fOrigAOD = aod; // copy pointer to the current AliAODEvent in the data member fOrigAOD
  if (!aod) return;
  Int_t nTracks = fOrigAOD->GetNumberOfTracks();
  printf("+++\n+++ Number of tracks in Event---> %d\n+++\n",nTracks);

  if(!arrayJPSItoEle) {
    printf("AliAnalysisTaskSEJPSItoEle::UserExec: JPSItoEle branch not found!\n");
    return;
  }
  if(!arrayLikeSign) {
    printf("AliAnalysisTaskSEJPSItoEle::UserExec: LikeSign2Prong branch not found!\n");
    return;
  }

  // load MC particles and read MC info (for sim only)
  TClonesArray* mcArray =  dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcArray) AliError("Could not find Monte-Carlo in AOD");
     //AliInfo(Form("+++\n+++ MC particles found in mcArray ---> %d \n+++\n",mcArray->GetEntriesFast()));
  if(fOkAODMC) ReadAODMCInfo(aod,arrayJPSItoEle);

  // retrieve AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();

  Int_t iOutVerticesHF = 0, iOutJPSItoEle = 0, iOutLikeSign = 0, iOutTracks = 0; 

  fVerticesHFTClArr->Delete();
  iOutVerticesHF = fVerticesHFTClArr->GetEntriesFast();
  TClonesArray &verticesHFRef = *fVerticesHFTClArr;

  fJpsiToEleTClArr->Delete();
  iOutJPSItoEle = fJpsiToEleTClArr->GetEntriesFast();
  TClonesArray &arrayJpsiToEleRef = *fJpsiToEleTClArr;

  fLikeSignTClArr->Delete();
  iOutLikeSign = fLikeSignTClArr->GetEntriesFast();
  TClonesArray &arrayLikeSignRef = *fLikeSignTClArr;

  fTracksTClArr->Delete();
  iOutTracks = fTracksTClArr->GetEntriesFast();
  //TClonesArray &arrayTrackRef = *fTracksTClArr;

  //
  // LOOP OVER LIKE SIGN CANDIDATES
  //

  // Access to the new AOD container of tracks
  Int_t trkIDtoEntry[100000];
  for(Int_t it=0;it<aod->GetNumberOfTracks();it++) {
    AliAODTrack *track = aod->GetTrack(it);
    trkIDtoEntry[track->GetID()]=it;
  }

  Int_t nPosPairs=0,nNegPairs=0;
  Int_t nLikeSign = arrayLikeSign->GetEntriesFast();
  if(fDebug>1) printf("+++\n+++ Number of total like sign pairs (before cuts)---> %d \n+++\n", nLikeSign);

  for(Int_t iLikeSign = 0; iLikeSign < nLikeSign; iLikeSign++) {
    AliAODRecoDecayHF2Prong *dlsin = (AliAODRecoDecayHF2Prong*)arrayLikeSign->UncheckedAt(iLikeSign);
    Bool_t unsetvtx=kFALSE;
    if(!dlsin->GetOwnPrimaryVtx()) {
        dlsin->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
        unsetvtx=kTRUE;
    }
    Int_t okBtoJPSIls=0;
    if(dlsin->SelectBtoJPSI(fVHF->GetBtoJPSICuts(),okBtoJPSIls)) {
       AliAODTrack *trk0 = (AliAODTrack*)dlsin->GetDaughter(0);
       fHistMassLS->Fill(dlsin->InvMassJPSIee());
       fHistCPtaLS->Fill(dlsin->CosPointingAngle());
       fHistd0d0LS->Fill(1e8*dlsin->Prodd0d0());
       fHistDCALS->Fill(100*dlsin->GetDCA());
       if(!trk0) {
          trk0=aod->GetTrack(trkIDtoEntry[dlsin->GetProngID(0)]);
          printf("references to standard AOD not available \n");
       }
       if((trk0->Charge())==1) {
          nPosPairs++;
          fHistCtsLS->Fill(dlsin->CosThetaStar(0,443,11,11));
          fHistCtsLSpos->Fill(dlsin->CosThetaStar(0,443,11,11));
        } else {
          nNegPairs++;
          fHistCtsLS->Fill(dlsin->CosThetaStarJPSI());
          fHistCtsLSneg->Fill(dlsin->CosThetaStarJPSI());
        }
       PostData(1,fOutput);
     
       // Clone like sign candidate for output AOD
       new(arrayLikeSignRef[iOutLikeSign++]) AliAODRecoDecayHF2Prong(*dlsin);

    }
    if(unsetvtx) dlsin->UnsetOwnPrimaryVtx();
  }

  if(fDebug>1) printf("+++\n+++ N. of positive pairs passing cuts in Event ----- %d \n+++\n", nPosPairs);
  if(fDebug>1) printf("+++\n+++ N. of negative pairs passing cuts in Event ----- %d \n+++\n", nNegPairs);

  fTotPosPairs += nPosPairs;
  fTotNegPairs += nNegPairs;

  //
  // LOOP OVER J/psi->e+e- CANDIDATES
  //

  Int_t nInJPSItoEle = arrayJPSItoEle->GetEntriesFast();
  if(fDebug>1) printf("+++\n+++ Number of total like JPSI -> ee candidates (before cuts)---> %d \n+++\n", nInJPSItoEle);

  //totJPSIin +=  nInJPSItoEle;

  for (Int_t iJPSItoEle = 0; iJPSItoEle < nInJPSItoEle; iJPSItoEle++) {

    AliAODRecoDecayHF2Prong *dIn = (AliAODRecoDecayHF2Prong*)arrayJPSItoEle->UncheckedAt(iJPSItoEle);
    Int_t mcLabel = 0;
    if(fOkAODMC) mcLabel = dIn->MatchToMC(443,mcArray) ; // select only reco JPSI that are true JPSI

    //Secondary vertex
    Double_t vtxSec[3] = {0.,0.,0.};
    Double_t vtxPrim[3] = {0.,0.,0.};
    Double_t vtxSecOut[3] = {0.,0.,0.};
    Double_t vtxPrimOut[3] = {0.,0.,0.};
    dIn->GetSecondaryVtx(vtxSec); 
    Bool_t unsetvtx=kFALSE;
    if(!dIn->GetOwnPrimaryVtx()) {
      dIn->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }

    Int_t okBtoJPSI=0;
    if(dIn->SelectBtoJPSI(fVHF->GetBtoJPSICuts(),okBtoJPSI)) {
      if ( fOkAODMC && mcLabel == -1){AliDebug(2,"No MC particle found");} else {

         fhInvMass->Fill(dIn->InvMassJPSIee()); 
         fhD0->Fill(10000*dIn->ImpParXY());
         fhD0D0->Fill(1e8*dIn->Prodd0d0());
         fhCosThetaStar->Fill(dIn->CosThetaStarJPSI());      
         fhCtsVsD0D0->Fill(dIn->CosThetaStarJPSI(),1e8*dIn->Prodd0d0());
         fhCosThetaPointing->Fill(dIn->CosPointingAngle());
         fhDCA->Fill(100*dIn->GetDCA());

         // compute X variable
         AliAODVertex* primVertex = dIn->GetOwnPrimaryVtx();
         vtxPrim[0] = primVertex->GetX();
         vtxPrim[1] = primVertex->GetY();
         vtxPrim[2] = primVertex->GetZ();
         Double_t lxy = ((vtxSec[0]-vtxPrim[0])*(dIn->Px()) + (vtxSec[1]-vtxPrim[1])*(dIn->Py()))/dIn->Pt();
         Double_t pseudoProperDecayTime = lxy*(TDatabasePDG::Instance()->GetParticle(443)->Mass())/dIn->Pt();
         fhDecayTime->Fill(10000*pseudoProperDecayTime);
  
         // clone candidate for output AOD
         AliAODVertex *v = new(verticesHFRef[iOutVerticesHF++]) 
           AliAODVertex(*(dIn->GetSecondaryVtx()));
         AliAODRecoDecayHF2Prong* dOut = new(arrayJpsiToEleRef[iOutJPSItoEle++]) 
           AliAODRecoDecayHF2Prong(*dIn); // copy constructor used
         dOut->SetSecondaryVtx(v);
         dOut->SetOwnPrimaryVtx((AliAODVertex*)((dIn->GetOwnPrimaryVtx())->Clone()));
         v->SetParent(dOut);

         // compute X variable from the cloned candidate in the stand-alone AOD file
         dOut->GetSecondaryVtx(vtxSecOut);
         AliAODVertex* primVertexOut = dOut->GetOwnPrimaryVtx();
         vtxPrimOut[0] = primVertexOut->GetX();
         vtxPrimOut[1] = primVertexOut->GetY();
         vtxPrimOut[2] = primVertexOut->GetZ();
         Double_t lxyOut = ((vtxSecOut[0]-vtxPrimOut[0])*(dOut->Px()) + (vtxSecOut[1]-vtxPrimOut[1])*(dOut->Py()))/dOut->Pt();
         Double_t pseudoProperDecayTimeOut = lxyOut*(TDatabasePDG::Instance()->GetParticle(443)->Mass())/dOut->Pt();
         fhDecayTimeOut->Fill(10000*pseudoProperDecayTimeOut);

         // retrieve tracks using the clone J/psi-->e+e- candidate stored in the output AOD
         //AliAODTrack* daugh0Out = (AliAODTrack*)dOut->GetSecondaryVtx()->GetDaughter(0);
         //AliAODTrack* daugh1Out= (AliAODTrack*)dOut->GetSecondaryVtx()->GetDaughter(1);
         //printf("pt of positive track: %f\n",daugh0Out->Pt());
         //printf("pt of negative track: %f\n",daugh1Out->Pt());

        }

     } // end of JPSItoEle candidates selection according to cuts

    if(unsetvtx) dIn->UnsetOwnPrimaryVtx();

    PostData(0,fOutput);

 }// end loop on JPSI to ele candidates

 printf("+++\n+++ Number of selected J/psi->e+e-: %d\n+++\n",iOutJPSItoEle);
 
 //totJPSIout += iOutJPSItoEle;

 return;

}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEJPSItoEle::Terminate(Option_t */*option*/)
{
  //
  // Terminate analysis
  //

  if(fDebug > 1) printf("AliAnalysisTaskSEJPSItoEle: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }

  fLsNormalization = 2.*TMath::Sqrt(fTotPosPairs*fTotNegPairs);

  if(fOkAODMC) fhDecayTimeMCjpsifromB = dynamic_cast<TH1F*>(fOutput->FindObject("fhDecayTimeMCjpsifromB"));
  fhDecayTime = dynamic_cast<TH1F*>(fOutput->FindObject("fhDecayTime"));
  fhDecayTimeOut = dynamic_cast<TH1F*>(fOutput->FindObject("fhDecayTimeOut"));
  fhInvMass = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMass"));
  fhD0 = dynamic_cast<TH1F*>(fOutput->FindObject("fhD0"));
  fhD0D0 = dynamic_cast<TH1F*>(fOutput->FindObject("fhD0D0"));
  fhCosThetaStar = dynamic_cast<TH1F*>(fOutput->FindObject("fhCosThetaStar"));
  fhCosThetaPointing = dynamic_cast<TH1F*>(fOutput->FindObject("fhCosThetaPointing"));
  fhCtsVsD0D0 = dynamic_cast<TH2F*>(fOutput->FindObject("fhCtsVsD0D0"));
  fHistMassLS = dynamic_cast<TH1F*>(fOutput->FindObject("fHistMassLS"));
  fHistCtsLS = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCtsLS"));
  fHistCtsLSpos = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCtsLSpos"));
  fHistCtsLSneg = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCtsLSneg"));
  fHistCPtaLS = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCPtaLS"));
  fHistd0d0LS = dynamic_cast<TH1F*>(fOutput->FindObject("fHistd0d0LS"));
  fhDCA = dynamic_cast<TH1F*>(fOutput->FindObject("fhDCA"));
  fHistDCALS = dynamic_cast<TH1F*>(fOutput->FindObject("fHistDCALS"));

  if(fLsNormalization>0.) {
    fHistMassLS->Scale((1/fLsNormalization)*fHistMassLS->GetEntries());
    fHistCtsLS->Scale((1/fLsNormalization)*fHistCtsLS->GetEntries());
    fHistCtsLSpos->Scale((1/fLsNormalization)*fHistCtsLSpos->GetEntries());
    fHistCtsLSneg->Scale((1/fLsNormalization)*fHistCtsLSneg->GetEntries());
    fHistCPtaLS->Scale((1/fLsNormalization)*fHistCPtaLS->GetEntries());
    fHistd0d0LS->Scale((1/fLsNormalization)*fHistd0d0LS->GetEntries());
    fHistDCALS->Scale((1/fLsNormalization)*fHistDCALS->GetEntries());
  }

  //printf("Tot JPSI in %d\n", totJPSIin);
  //printf("Tot JPSI out %d\n", totJPSIout);

  return;
}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEJPSItoEle::SetPtCuts(const Double_t ptCuts[2])
{
  //
  // apply Pt cuts (lower cut, upper cut)
  //
  for(Int_t i = 0; i < 2; i++) fPtCuts[i] = ptCuts[i];
}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEJPSItoEle::SetCutsJPSI(const Double_t cuts[9])
{
  //
  // apply JPSI cuts 
  //
  for(Int_t i = 0; i < 9; i++) fCuts[i] = cuts[i];
}
//_________________________________________________________________________________________________
void AliAnalysisTaskSEJPSItoEle::ReadAODMCInfo(AliAODEvent* aodEvent, const TClonesArray* inputArray) 
{
  //
  // Reads MC information if needed (for example in case of parametrized PID)
  //

  // load MC particles
  TClonesArray* mcArray =
     dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcArray) AliError("Could not find Monte-Carlo in AOD");
     AliInfo(Form("MC particles found in mcArray ---> %d ",mcArray->GetEntriesFast()));

  // load MC header 
  AliAODMCHeader* mcHeader =
    (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader){
    printf("AliAnalysisTaskSEJPSItoEle::UserExec: MC AOD header branch not found!\n");
  }

  AliAODVertex *vtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();

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
    Int_t mcLabelSec = dd->MatchToMC(443,mcArray) ; // select only reco JPSI that are true JPSI

    Double_t vtxPrim[3] = {0.,0.,0.};
    Double_t vtxSec[3] = {0.,0.,0.};
    dd->GetSecondaryVtx(vtxSec);
    Bool_t unsetvtx=kFALSE;
    if(!dd->GetOwnPrimaryVtx()) {
      dd->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }

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

      if (mcLabelSec == -1)
         {
          AliDebug(2,"No MC particle found");
         }
      else {

      if(labJPSIdaugh0>=0 && labJPSIdaugh1>=0 && labJPSIdaugh0==labJPSIdaugh1) { // check if JPSI is signal
          AliAODMCParticle* partJPSI = (AliAODMCParticle*)mcArray->At(labJPSIdaugh0);
          pdgJPSI = partJPSI->GetPdgCode();
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

        // chech if JPSI comes from B 
           if( pdgJPSIMother==511   || pdgJPSIMother==521   ||
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

                  fhInvMass->Fill(dd->InvMassJPSIee()); 
                  fhD0->Fill(10000*dd->ImpParXY());
                  fhD0D0->Fill(1e8*dd->Prodd0d0());
                  fhCosThetaStar->Fill(dd->CosThetaStarJPSI());      
                  fhCtsVsD0D0->Fill(dd->CosThetaStarJPSI(),1e8*dd->Prodd0d0());
                  fhCosThetaPointing->Fill(dd->CosPointingAngle());
  
                  // compute X variable
                  AliAODVertex* primVertex = dd->GetOwnPrimaryVtx();
                  vtxPrim[0] = primVertex->GetX();
                  vtxPrim[1] = primVertex->GetY();
                  vtxPrim[2] = primVertex->GetZ();
                  Double_t lxy = ((vtxSec[0]-vtxPrim[0])*(dd->Px()) + (vtxSec[1]-vtxPrim[1])*(dd->Py()))/dd->Pt();
                  Double_t pseudoProperDecayTime = lxy*(TDatabasePDG::Instance()->GetParticle(443)->Mass())/dd->Pt();
  
                  fhDecayTime->Fill(10000*pseudoProperDecayTime);
                  // Post the data already here
                  PostData(1,fOutput);

                  Double_t mcLxy = ((partJPSI->Xv()-mcPrimVtx[0])*(partJPSI->Px()) + (partJPSI->Yv()-mcPrimVtx[1])*(partJPSI->Py()))/partJPSI->Pt();
                  Double_t mcPseudoProperDecayTime = mcLxy*(TDatabasePDG::Instance()->GetParticle(443)->Mass())/partJPSI->Pt();
                  fhDecayTimeMCjpsifromB->Fill(10000*mcPseudoProperDecayTime,1);

                  // Post the data already here
                  PostData(1,fOutput);

            } //this is for MC JPSI -> ee from B-hadrons
          } //this is for MC JPSI -> ee
        }
      } // end of check if JPSI is signal
    }
  }
}
