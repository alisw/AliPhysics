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

/////////////////////////////////////////////////////////////
//
// Class AliAnalysisTaskSECharmFraction
// AliAnalysisTaskSE for the extraction of the fraction of prompt charm
// using the charm hadron impact parameter to the primary vertex
//
// Author: Andrea Rossi, andrea.rossi@ts.infn.it
/////////////////////////////////////////////////////////////


#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TROOT.h>

#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecay.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSECharmFraction.h"


class TCanvas;
class TTree;
class TChain;
class AliAODInputHandler;
class AliAnalysisManager;
class AliAnalysisTaskSE;


ClassImp(AliAnalysisTaskSECharmFraction)
 
//________________________________________________________________________
  AliAnalysisTaskSECharmFraction::AliAnalysisTaskSECharmFraction() 
    : AliAnalysisTaskSE(),
      fVHFloose(0),
      fVHFtight(0),
      fmD0PDG(),
      fnbins(),
      fptbins(0),
      fAcceptanceCuts(),
      fsignalInvMassCut(),
      flargeInvMassCut(),
      fsidebandInvMassCut(),
      fsidebandInvMassWindow(),
      fUseMC(kTRUE),
      fNentries(0),
      fSignalType(0),
      fSignalTypeLsCuts(0),
      fSignalTypeTghCuts(0),
      flistNoCutsSignal(0),
      flistNoCutsBack(0),
      flistNoCutsFromB(0),
      flistNoCutsFromDstar(0),
      flistNoCutsOther(0),
      flistLsCutsSignal(0),
      flistLsCutsBack(0),
      flistLsCutsFromB(0),
      flistLsCutsFromDstar(0),
      flistLsCutsOther(0),
      flistTghCutsSignal(0),
      flistTghCutsBack(0),
      flistTghCutsFromB(0),
      flistTghCutsFromDstar(0),
      flistTghCutsOther(0)
   
{
  //Default constructor
}
//________________________________________________________________________
  AliAnalysisTaskSECharmFraction::AliAnalysisTaskSECharmFraction(const char *name) 
    : AliAnalysisTaskSE(name),
      fVHFloose(0),
      fVHFtight(0),
      fmD0PDG(),
      fnbins(),
      fptbins(0),
      fAcceptanceCuts(),
      fsignalInvMassCut(),
      flargeInvMassCut(),
      fsidebandInvMassCut(),
      fsidebandInvMassWindow(),
      fUseMC(kTRUE),
      fNentries(0),
      fSignalType(0),
      fSignalTypeLsCuts(0),
      fSignalTypeTghCuts(0),
      flistNoCutsSignal(0),
      flistNoCutsBack(0),
      flistNoCutsFromB(0),
      flistNoCutsFromDstar(0),
      flistNoCutsOther(0),
      flistLsCutsSignal(0),
      flistLsCutsBack(0),
      flistLsCutsFromB(0),
      flistLsCutsFromDstar(0),
      flistLsCutsOther(0),
      flistTghCutsSignal(0),
      flistTghCutsBack(0),
      flistTghCutsFromB(0),
      flistTghCutsFromDstar(0),
      flistTghCutsOther(0)
   
{
  // Constructor
 
  // Define input and output slots here
  // Input slot #0 works with a TChain
  // Output slot #0 writes into a TH1 container

  //Standard pt bin
  fnbins=4;
  fptbins=new Double_t[fnbins+1];
  fptbins[0]=0.;
  fptbins[1]=1.;
  fptbins[2]=3.;
  fptbins[3]=5.;
  fptbins[4]=1000.;
  //fAcceptanceCuts=new Double_t[3];
  SetAcceptanceCut();
  SetStandardMassSelection();
  DefineOutput(1, TH1F::Class());
  DefineOutput(2, TH1F::Class());
  DefineOutput(3, TH1F::Class());
  DefineOutput(4, TH1F::Class());
  for(Int_t j=5;j<20;j++){
    DefineOutput(j, TList::Class());
  }


}


AliAnalysisTaskSECharmFraction::AliAnalysisTaskSECharmFraction(const char *name,Int_t nptbins,Double_t *ptbins) 
  : AliAnalysisTaskSE(name),
    fVHFloose(0),
    fVHFtight(0),
    fmD0PDG(),
    fnbins(),
    fptbins(0),
    fAcceptanceCuts(),
    fsignalInvMassCut(),
    flargeInvMassCut(),
    fsidebandInvMassCut(),
    fsidebandInvMassWindow(),
    fUseMC(kTRUE),
    fNentries(0),
    fSignalType(0),
    fSignalTypeLsCuts(0),
    fSignalTypeTghCuts(0),
    flistNoCutsSignal(0),
    flistNoCutsBack(0),
    flistNoCutsFromB(0),
    flistNoCutsFromDstar(0),
    flistNoCutsOther(0),
    flistLsCutsSignal(0),
    flistLsCutsBack(0),
    flistLsCutsFromB(0),
    flistLsCutsFromDstar(0),
    flistLsCutsOther(0),
    flistTghCutsSignal(0),
    flistTghCutsBack(0),
    flistTghCutsFromB(0),
    flistTghCutsFromDstar(0),
    flistTghCutsOther(0)
{
  // Constructor
  // ptbins must be of dimension nptbins +1
  
  SetNPtBins(nptbins,ptbins);
  SetStandardMassSelection();
  //  fAcceptanceCuts=new Double_t[3];
  SetAcceptanceCut();
  // Define input and output slots here
 
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TH1F::Class());
  DefineOutput(2, TH1F::Class());
  DefineOutput(3, TH1F::Class());
  DefineOutput(4, TH1F::Class());
  for(Int_t j=5;j<20;j++){

    DefineOutput(j, TList::Class());
  }

 
}

//________________________________________________________________________
AliAnalysisTaskSECharmFraction::~AliAnalysisTaskSECharmFraction()
{ //Destructor 
  
  if (fVHFtight) {
    delete fVHFtight;
    fVHFtight = 0;
  }
  if (fVHFloose) {
    delete fVHFloose;
    fVHFloose = 0;
  }
  if(fptbins){
    delete fptbins;
    fptbins =0;
  }
  /*  if(fAcceptanceCuts){
    delete fAcceptanceCuts;
    fAcceptanceCuts=0;
    }*/
  if (fNentries) {
    delete fNentries;
    fNentries = 0;
  }   
  if (fSignalType) {
    delete fSignalType;
    fSignalType = 0;
  } 
  if (fSignalTypeLsCuts) {
    delete fSignalTypeLsCuts;
    fSignalTypeLsCuts = 0;
  } 
  if (fSignalTypeTghCuts) {
    delete fSignalTypeTghCuts;
    fSignalTypeTghCuts = 0;
  } 
  if(flistNoCutsSignal){
    delete flistNoCutsSignal;
    flistNoCutsSignal=0;
  }
  if(flistNoCutsBack){
    delete flistNoCutsBack;
    flistNoCutsBack=0;
  }
  if(flistNoCutsFromB){
    delete flistNoCutsFromB;
    flistNoCutsFromB=0;
  }
  if(flistNoCutsFromDstar){
    delete flistNoCutsFromDstar;
    flistNoCutsFromDstar=0;
  }
  if(flistNoCutsOther){
    delete flistNoCutsOther;
    flistNoCutsOther=0;
  }
  
 if(flistLsCutsSignal){
    delete flistLsCutsSignal;
    flistLsCutsSignal=0;
  }
  if(flistLsCutsBack){
    delete flistLsCutsBack;
    flistLsCutsBack=0;
  }
  if(flistLsCutsFromB){
    delete flistLsCutsFromB;
    flistLsCutsFromB=0;
  }
  if(flistLsCutsFromDstar){
    delete flistLsCutsFromDstar;
    flistLsCutsFromDstar=0;
  }
  if(flistLsCutsOther){
    delete flistLsCutsOther;
    flistLsCutsOther=0;
  }
  
 if(flistTghCutsSignal){
    delete flistTghCutsSignal;
    flistTghCutsSignal=0;
  }
  if(flistTghCutsBack){
    delete flistTghCutsBack;
    flistTghCutsBack=0;
  }
  if(flistTghCutsFromB){
    delete flistTghCutsFromB;
    flistTghCutsFromB=0;
  }
  if(flistTghCutsFromDstar){
    delete flistTghCutsFromDstar;
    flistTghCutsFromDstar=0;
  }
  if(flistTghCutsOther){
    delete flistTghCutsOther;
    flistTghCutsOther=0;
  }
  
  
}  


//________________________________________________________________________
void AliAnalysisTaskSECharmFraction::Init()
{
  // Initialization
  
  if(fDebug > 1) printf("AnalysisTaskSED0Mass::Init() \n");
  fmD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF.C");
  
  // 2 sets of dedidcated cuts -- defined in UserExec
  //  the config file and the way the cuts are set is for further development
  //   (to be interfaced with AliAnalsysTaskSETuneCuts)

  fVHFtight = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");
  fVHFloose = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");  
  if(!fptbins){
    //SET STANDARD PT BINNING
    fnbins=4;
    fptbins=new Double_t[fnbins+1];
    fptbins[0]=0.;
    fptbins[1]=1.;
    fptbins[2]=3.;
    fptbins[3]=5.;
    fptbins[4]=1000.;
  }
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmFraction::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  // ################ NAMING SCHEME ###################################
  //            LISTS NAMING SCHEME
  // "list" + cut selection string + MC selection string
  //      cut strings:  "NC" =nocuts, "LSC"= loose cuts, "TGHC"= tight cuts
  //      MC sel. strings: "sign"= D0 from c quark
  //                       "fromDstar" = D0 from Dstar from c quark
  //                       "fromB"= D0from B decay (->from b quark) + D0from Dstar from B
  //                       "back"= backgroun, generic except the cas "other"
  //                       "other"= background case for candidates made of a pion and a kaon coming from the same D0 (in 4 prong) or from D+
  //
  //           HISTS NAMING SCHEME 
  // 
  //  "h" + specific name + cut selection string + MC selection string + (InvMass region string) + (pt string)
  //
  //        cut selection strings    = those for lists
  //        MC selection strings     = those for lists
  //        inv mass region strings  : "PM" or "SB" for global properties and pt integrated histos
  //                                   "_PkMss" or "_SBMss" for impact par. pt dependent histos
  //                   pt string     : "_pt" + integer number of ptbin
  //
  //###################################################################

  TString namehist;
  TString titlehist;
  TString strnamept,strtitlept;
 
  fNentries=new TH1F("nentriesChFr", "Look at the number of entries! = number of AODs", 2,1.,2.);
  fSignalType=new TH1F("hsignaltype", "Histo for type of MC signal", 21,-1.,20.);
  fSignalTypeLsCuts=new TH1F("hsignaltypeLsCuts", "Histo for type of MC signal with loose cuts", 21,-1.,20.);
  fSignalTypeTghCuts=new TH1F("hsignaltypeTghCuts", "Histo for type of MC signal with tight cuts", 21,-1.,20.);

  //##########  DEFINE THE TLISTS ##################
  
  flistNoCutsSignal = new TList();
  flistNoCutsSignal->SetOwner();
  flistNoCutsSignal->SetName("listNCsign");

  flistNoCutsBack = new TList();
  flistNoCutsBack->SetOwner();
  flistNoCutsBack->SetName("listNCback");

  flistNoCutsFromB = new TList();
  flistNoCutsFromB->SetOwner();
  flistNoCutsFromB->SetName("listNCfromB");

  flistNoCutsFromDstar = new TList();
  flistNoCutsFromDstar->SetOwner();
  flistNoCutsFromDstar->SetName("listNCfromDstar");

  flistNoCutsOther = new TList();
  flistNoCutsOther->SetOwner();
  flistNoCutsOther->SetName("listNCother");


  flistLsCutsSignal = new TList();
  flistLsCutsSignal->SetOwner();
  flistLsCutsSignal->SetName("listLSCsign");

  flistLsCutsBack = new TList();
  flistLsCutsBack->SetOwner();
  flistLsCutsBack->SetName("listLSCback");

  flistLsCutsFromB = new TList();
  flistLsCutsFromB->SetOwner();
  flistLsCutsFromB->SetName("listLSCfromB");

  flistLsCutsFromDstar = new TList();
  flistLsCutsFromDstar->SetOwner();
  flistLsCutsFromDstar->SetName("listLSCfromDstar");

  flistLsCutsOther = new TList();
  flistLsCutsOther->SetOwner();
  flistLsCutsOther->SetName("listLSCother");


  flistTghCutsSignal = new TList();
  flistTghCutsSignal->SetOwner();
  flistTghCutsSignal->SetName("listTGHCsign");

  flistTghCutsBack = new TList();
  flistTghCutsBack->SetOwner();
  flistTghCutsBack->SetName("listTGHCback");

  flistTghCutsFromB = new TList();
  flistTghCutsFromB->SetOwner();
  flistTghCutsFromB->SetName("listTGHCfromB");

  flistTghCutsFromDstar = new TList();
  flistTghCutsFromDstar->SetOwner();
  flistTghCutsFromDstar->SetName("listTGHCfromDstar");

  flistTghCutsOther = new TList();
  flistTghCutsOther->SetOwner();
  flistTghCutsOther->SetName("listTGHCother");




  //################################################################################################
  //                                                                                               #
  //                         HISTOS FOR NO CUTS CASE                                               #
  //                                                                                               #
  //################################################################################################


  //############ NO CUTS SIGNAL HISTOGRAMS ###############
  //
  // ####### global properties histo ############

  TH2F *hCPtaVSd0d0NCsign=new TH2F("hCPtaVSd0d0NCsign","hCPtaVSd0d0_NoCuts_Signal",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZNCsign=new TH1F("hSecVtxZNCsign","hSecVtxZ_NoCuts_Signal",1000,-8.,8.);
  TH1F *hSecVtxXNCsign=new TH1F("hSecVtxXNCsign","hSecVtxX_NoCuts_Signal",1000,-3000.,3000.);
  TH1F *hSecVtxYNCsign=new TH1F("hSecVtxYNCsign","hSecVtxY_NoCuts_Signal",1000,-3000.,3000.);
  TH2F *hSecVtxXYNCsign=new TH2F("hSecVtxXYNCsign","hSecVtxXY_NoCuts_Signal",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiNCsign=new TH1F("hSecVtxPhiNCsign","hSecVtxPhi_NoCuts_Signal",180,-180.1,180.1);
  TH1F *hCPtaNCsign=new TH1F("hCPtaNCsign","hCPta_NoCuts_Signal",100,0.,1.);
  TH1F *hd0xd0NCsign=new TH1F("hd0xd0NCsign","hd0xd0_NoCuts_Signal",1000,-100000.,100000.);
  TH1F *hMassTrueNCsign=new TH1F("hMassTrueNCsign","D^{0} MC inv. Mass No Cuts Signal(All momenta)",600,1.600,2.200);
  TH1F *hMassNCsign=new TH1F("hMassNCsign","D^{0} inv. Mass No Cuts Signal (All momenta)",600,1.600,2.200);
  hMassNCsign->Sumw2();
  TH1F *hMassTrueNCsignPM=new TH1F("hMassTrueNCsignPM","D^{0} MC inv. Mass No Cuts Signal, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassNCsignPM=new TH1F("hMassNCsignPM","D^{0} inv. Mass No Cuts Signal (All momenta), MassPeak",600,1.600,2.200);
  hMassNCsignPM->Sumw2();

  TH1F *hMassTrueNCsignSB=new TH1F("hMassTrueNCsignSB","D^{0} MC inv. Mass in Side Bands No Cuts Signal(All momenta)",600,1.600,2.200);
  TH1F *hMassNCsignSB=new TH1F("hMassNCsignSB","D^{0} inv. Mass in Side Bands No Cuts Signal (All momenta)",600,1.600,2.200);
  hMassNCsignSB->Sumw2();

  flistNoCutsSignal->Add(hCPtaVSd0d0NCsign);
  flistNoCutsSignal->Add(hSecVtxZNCsign);
  flistNoCutsSignal->Add(hSecVtxYNCsign);
  flistNoCutsSignal->Add(hSecVtxXNCsign);
  flistNoCutsSignal->Add(hSecVtxXYNCsign);
  flistNoCutsSignal->Add(hSecVtxPhiNCsign);
  flistNoCutsSignal->Add(hCPtaNCsign);
  flistNoCutsSignal->Add(hd0xd0NCsign);
  flistNoCutsSignal->Add(hMassTrueNCsign);
  flistNoCutsSignal->Add(hMassNCsign);
  flistNoCutsSignal->Add(hMassTrueNCsignPM);
  flistNoCutsSignal->Add(hMassNCsignPM);
  flistNoCutsSignal->Add(hMassTrueNCsignSB);
  flistNoCutsSignal->Add(hMassNCsignSB);

  // ####### d0 D0 histos ############
  TH1F *hd0D0NCsignPM = new TH1F("hd0D0NCsignPM","D^{0} impact par. plot , No Cuts ,Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0NCsignPM->SetXTitle("Impact parameter [#mum]");
  hd0D0NCsignPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueNCsignPM = new TH1F("hd0D0VtxTrueNCsignPM","D^{0} impact par. w.r.t. True Vtx, No Cuts, Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueNCsignPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueNCsignPM->SetYTitle("Entries");

  TH1F *hMCd0D0NCsignPM = new TH1F("hMCd0D0NCsignPM","D^{0} impact par. plot, No Cuts, Signal,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0NCsignPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0NCsignPM->SetYTitle("Entries");

  TH1F *hd0D0NCsignSB = new TH1F("hd0D0NCsignSB","D^{0} impact par. plot , No Cuts ,Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0NCsignSB->SetXTitle("Impact parameter [#mum]");
  hd0D0NCsignSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueNCsignSB = new TH1F("hd0D0VtxTrueNCsignSB","D^{0} impact par. w.r.t. True Vtx, No Cuts, Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueNCsignSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueNCsignSB->SetYTitle("Entries");

  TH1F *hMCd0D0NCsignSB = new TH1F("hMCd0D0NCsignSB","D^{0} impact par. plot, No Cuts, Signal,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0NCsignSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0NCsignSB->SetYTitle("Entries");

  flistNoCutsSignal->Add(hd0D0NCsignPM);
  flistNoCutsSignal->Add(hd0D0VtxTrueNCsignPM);
  flistNoCutsSignal->Add(hMCd0D0NCsignPM);
  flistNoCutsSignal->Add(hd0D0NCsignSB);
  flistNoCutsSignal->Add(hd0D0VtxTrueNCsignSB);
  flistNoCutsSignal->Add(hMCd0D0NCsignSB);
  
  TH1F **hd0D0ptNCsignPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptNCsignPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptNCsignPM=new TH1F*[fnbins];
  TH1F **hd0D0ptNCsignSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptNCsignSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptNCsignSB=new TH1F*[fnbins];
  namehist="hd0D0ptNCsign_";
  titlehist="D^{0} impact par. plot, No Cuts, Signal, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCsignPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCsignPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCsignPM[i]->SetYTitle("Entries");
    flistNoCutsSignal->Add(hd0D0ptNCsignPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCsignPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCsignPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCsignPM[i]->SetYTitle("Entries");
    flistNoCutsSignal->Add(hMCd0D0ptNCsignPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCsignPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCsignPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCsignPM[i]->SetYTitle("Entries");
    flistNoCutsSignal->Add(hd0D0VtxTrueptNCsignPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCsignSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCsignSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCsignSB[i]->SetYTitle("Entries");
    flistNoCutsSignal->Add(hd0D0ptNCsignSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCsignSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCsignSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCsignSB[i]->SetYTitle("Entries");
    flistNoCutsSignal->Add(hMCd0D0ptNCsignSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCsignSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCsignSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCsignSB[i]->SetYTitle("Entries");
    flistNoCutsSignal->Add(hd0D0VtxTrueptNCsignSB[i]);
  }


  //############ NO CUTS BACKGROUND HISTOGRAMS ###########
  //
  //   ######## global properties histos #######
  TH2F *hCPtaVSd0d0NCback=new TH2F("hCPtaVSd0d0NCback","hCPtaVSd0d0_NoCuts_Background",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZNCback=new TH1F("hSecVtxZNCback","hSecVtxZ_NoCuts_Background",1000,-8.,8.);
  TH1F *hSecVtxXNCback=new TH1F("hSecVtxXNCback","hSecVtxX_NoCuts_Background",1000,-3000.,3000.);
  TH1F *hSecVtxYNCback=new TH1F("hSecVtxYNCback","hSecVtxY_NoCuts_Background",1000,-3000.,3000.);
  TH2F *hSecVtxXYNCback=new TH2F("hSecVtxXYNCback","hSecVtxXY_NoCuts_Background",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiNCback=new TH1F("hSecVtxPhiNCback","hSecVtxPhi_NoCuts_Background",180,-180.1,180.1);
  TH1F *hCPtaNCback=new TH1F("hCPtaNCback","hCPta_NoCuts_Background",100,0.,1.);
  TH1F *hd0xd0NCback=new TH1F("hd0xd0NCback","hd0xd0_NoCuts_Background",1000,-100000.,100000.);
  TH1F *hMassTrueNCback=new TH1F("hMassTrueNCback","D^{0} MC inv. Mass No Cuts Background(All momenta)",600,1.600,2.200);
  TH1F *hMassNCback=new TH1F("hMassNCback","D^{0} inv. Mass No Cuts Background (All momenta)",600,1.600,2.200);
  hMassNCback->Sumw2();
  TH1F *hMassTrueNCbackPM=new TH1F("hMassTrueNCbackPM","D^{0} MC inv. Mass No Cuts Background, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassNCbackPM=new TH1F("hMassNCbackPM","D^{0} inv. Mass No Cuts Background (All momenta), MassPeak",600,1.600,2.200);
  hMassNCbackPM->Sumw2();
  TH1F *hMassTrueNCbackSB=new TH1F("hMassTrueNCbackSB","D^{0} MC inv. Mass in Side Bands No Cuts Background(All momenta)",600,1.600,2.200);
  TH1F *hMassNCbackSB=new TH1F("hMassNCbackSB","D^{0} inv. Mass in Side Bands No Cuts Background (All momenta)",600,1.600,2.200);
  hMassNCbackSB->Sumw2();

  flistNoCutsBack->Add(hCPtaVSd0d0NCback);
  flistNoCutsBack->Add(hSecVtxZNCback);
  flistNoCutsBack->Add(hSecVtxYNCback);
  flistNoCutsBack->Add(hSecVtxXNCback);
  flistNoCutsBack->Add(hSecVtxXYNCback);
  flistNoCutsBack->Add(hSecVtxPhiNCback);
  flistNoCutsBack->Add(hCPtaNCback);
  flistNoCutsBack->Add(hd0xd0NCback);
  flistNoCutsBack->Add(hMassTrueNCback);
  flistNoCutsBack->Add(hMassNCback);
  flistNoCutsBack->Add(hMassTrueNCbackPM);
  flistNoCutsBack->Add(hMassNCbackPM);
  flistNoCutsBack->Add(hMassTrueNCbackSB);
  flistNoCutsBack->Add(hMassNCbackSB);


  // ####### d0 D0 histos ############
  
 TH1F *hd0D0NCbackPM = new TH1F("hd0D0NCbackPM","D^{0} impact par. plot , No Cuts ,Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0NCbackPM->SetXTitle("Impact parameter [#mum]");
  hd0D0NCbackPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueNCbackPM = new TH1F("hd0D0VtxTrueNCbackPM","D^{0} impact par. w.r.t. True Vtx, No Cuts, Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueNCbackPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueNCbackPM->SetYTitle("Entries");

  TH1F *hMCd0D0NCbackPM = new TH1F("hMCd0D0NCbackPM","D^{0} impact par. plot, No Cuts, Background,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0NCbackPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0NCbackPM->SetYTitle("Entries");

  TH1F *hd0D0NCbackSB = new TH1F("hd0D0NCbackSB","D^{0} impact par. plot , No Cuts ,Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0NCbackSB->SetXTitle("Impact parameter [#mum]");
  hd0D0NCbackSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueNCbackSB = new TH1F("hd0D0VtxTrueNCbackSB","D^{0} impact par. w.r.t. True Vtx, No Cuts, Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueNCbackSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueNCbackSB->SetYTitle("Entries");

  TH1F *hMCd0D0NCbackSB = new TH1F("hMCd0D0NCbackSB","D^{0} impact par. plot, No Cuts, Background,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0NCbackSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0NCbackSB->SetYTitle("Entries");

  flistNoCutsBack->Add(hd0D0NCbackPM);
  flistNoCutsBack->Add(hd0D0VtxTrueNCbackPM);
  flistNoCutsBack->Add(hMCd0D0NCbackPM);
  flistNoCutsBack->Add(hd0D0NCbackSB);
  flistNoCutsBack->Add(hd0D0VtxTrueNCbackSB);
  flistNoCutsBack->Add(hMCd0D0NCbackSB);
  
  TH1F **hd0D0ptNCbackPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptNCbackPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptNCbackPM=new TH1F*[fnbins];
  TH1F **hd0D0ptNCbackSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptNCbackSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptNCbackSB=new TH1F*[fnbins];
  namehist="hd0D0ptNCback_";
  titlehist="D^{0} impact par. plot, No Cuts, Background, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCbackPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCbackPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCbackPM[i]->SetYTitle("Entries");
    flistNoCutsBack->Add(hd0D0ptNCbackPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCbackPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCbackPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCbackPM[i]->SetYTitle("Entries");
    flistNoCutsBack->Add(hMCd0D0ptNCbackPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCbackPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCbackPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCbackPM[i]->SetYTitle("Entries");
    flistNoCutsBack->Add(hd0D0VtxTrueptNCbackPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCbackSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCbackSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCbackSB[i]->SetYTitle("Entries");
    flistNoCutsBack->Add(hd0D0ptNCbackSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCbackSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCbackSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCbackSB[i]->SetYTitle("Entries");
    flistNoCutsBack->Add(hMCd0D0ptNCbackSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCbackSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCbackSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCbackSB[i]->SetYTitle("Entries");
    flistNoCutsBack->Add(hd0D0VtxTrueptNCbackSB[i]);
  }



 //############ NO CUTS FROMB HISTOGRAMS ###########
  //
  //#######  global properties histos

  TH2F *hCPtaVSd0d0NCfromB=new TH2F("hCPtaVSd0d0NCfromB","hCPtaVSd0d0_NoCuts_FromB",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZNCfromB=new TH1F("hSecVtxZNCfromB","hSecVtxZ_NoCuts_FromB",1000,-8.,8.);
  TH1F *hSecVtxXNCfromB=new TH1F("hSecVtxXNCfromB","hSecVtxX_NoCuts_FromB",1000,-3000.,3000.);
  TH1F *hSecVtxYNCfromB=new TH1F("hSecVtxYNCfromB","hSecVtxY_NoCuts_FromB",1000,-3000.,3000.);
  TH2F *hSecVtxXYNCfromB=new TH2F("hSecVtxXYNCfromB","hSecVtxXY_NoCuts_FromB",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiNCfromB=new TH1F("hSecVtxPhiNCfromB","hSecVtxPhi_NoCuts_FromB",180,-180.1,180.1);
  TH1F *hCPtaNCfromB=new TH1F("hCPtaNCfromB","hCPta_NoCuts_FromB",100,0.,1.);
  TH1F *hd0xd0NCfromB=new TH1F("hd0xd0NCfromB","hd0xd0_NoCuts_FromB",1000,-100000.,100000.);
  TH1F *hMassTrueNCfromB=new TH1F("hMassTrueNCfromB","D^{0} MC inv. Mass No Cuts FromB(All momenta)",600,1.600,2.200);
  TH1F *hMassNCfromB=new TH1F("hMassNCfromB","D^{0} inv. Mass No Cuts FromB (All momenta)",600,1.600,2.200);
  hMassNCfromB->Sumw2();
  TH1F *hMassTrueNCfromBPM=new TH1F("hMassTrueNCfromBPM","D^{0} MC inv. Mass No Cuts FromB, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassNCfromBPM=new TH1F("hMassNCfromBPM","D^{0} inv. Mass No Cuts FromB (All momenta), MassPeak",600,1.600,2.200);
  hMassNCfromB->Sumw2();
  TH1F *hMassTrueNCfromBSB=new TH1F("hMassTrueNCfromBSB","D^{0} MC inv. Mass in Side Bands No Cuts FromB(All momenta)",600,1.600,2.200);
  TH1F *hMassNCfromBSB=new TH1F("hMassNCfromBSB","D^{0} inv. Mass in Side Bands No Cuts FromB (All momenta)",600,1.600,2.200);
  hMassNCfromBSB->Sumw2();

  flistNoCutsFromB->Add(hCPtaVSd0d0NCfromB);
  flistNoCutsFromB->Add(hSecVtxZNCfromB);
  flistNoCutsFromB->Add(hSecVtxYNCfromB);
  flistNoCutsFromB->Add(hSecVtxXNCfromB);
  flistNoCutsFromB->Add(hSecVtxXYNCfromB);
  flistNoCutsFromB->Add(hSecVtxPhiNCfromB);
  flistNoCutsFromB->Add(hCPtaNCfromB);
  flistNoCutsFromB->Add(hd0xd0NCfromB);
  flistNoCutsFromB->Add(hMassTrueNCfromB);
  flistNoCutsFromB->Add(hMassNCfromB);
  flistNoCutsFromB->Add(hMassTrueNCfromBPM);
  flistNoCutsFromB->Add(hMassNCfromBPM);
  flistNoCutsFromB->Add(hMassTrueNCfromBSB);
  flistNoCutsFromB->Add(hMassNCfromBSB);

  // ######### d0 D0 histos ##############
  TH1F *hd0D0NCfromBPM = new TH1F("hd0D0NCfromBPM","D^{0} impact par. plot , No Cuts ,FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0NCfromBPM->SetXTitle("Impact parameter [#mum]");
  hd0D0NCfromBPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueNCfromBPM = new TH1F("hd0D0VtxTrueNCfromBPM","D^{0} impact par. w.r.t. True Vtx, No Cuts, FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueNCfromBPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueNCfromBPM->SetYTitle("Entries");

  TH1F *hMCd0D0NCfromBPM = new TH1F("hMCd0D0NCfromBPM","D^{0} impact par. plot, No Cuts, FromB,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0NCfromBPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0NCfromBPM->SetYTitle("Entries");

  TH1F *hd0D0NCfromBSB = new TH1F("hd0D0NCfromBSB","D^{0} impact par. plot , No Cuts ,FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0NCfromBSB->SetXTitle("Impact parameter [#mum]");
  hd0D0NCfromBSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueNCfromBSB = new TH1F("hd0D0VtxTrueNCfromBSB","D^{0} impact par. w.r.t. True Vtx, No Cuts, FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueNCfromBSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueNCfromBSB->SetYTitle("Entries");

  TH1F *hMCd0D0NCfromBSB = new TH1F("hMCd0D0NCfromBSB","D^{0} impact par. plot, No Cuts, FromB,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0NCfromBSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0NCfromBSB->SetYTitle("Entries");

  flistNoCutsFromB->Add(hd0D0NCfromBPM);
  flistNoCutsFromB->Add(hd0D0VtxTrueNCfromBPM);
  flistNoCutsFromB->Add(hMCd0D0NCfromBPM);
  flistNoCutsFromB->Add(hd0D0NCfromBSB);
  flistNoCutsFromB->Add(hd0D0VtxTrueNCfromBSB);
  flistNoCutsFromB->Add(hMCd0D0NCfromBSB);
  
  TH1F **hd0D0ptNCfromBPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptNCfromBPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptNCfromBPM=new TH1F*[fnbins];
  TH1F **hd0D0ptNCfromBSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptNCfromBSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptNCfromBSB=new TH1F*[fnbins];
  namehist="hd0D0ptNCfromB_";
  titlehist="D^{0} impact par. plot, No Cuts, FromB, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCfromBPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCfromBPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCfromBPM[i]->SetYTitle("Entries");
    flistNoCutsFromB->Add(hd0D0ptNCfromBPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCfromBPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCfromBPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCfromBPM[i]->SetYTitle("Entries");
    flistNoCutsFromB->Add(hMCd0D0ptNCfromBPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCfromBPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCfromBPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCfromBPM[i]->SetYTitle("Entries");
    flistNoCutsFromB->Add(hd0D0VtxTrueptNCfromBPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCfromBSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCfromBSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCfromBSB[i]->SetYTitle("Entries");
    flistNoCutsFromB->Add(hd0D0ptNCfromBSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCfromBSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCfromBSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCfromBSB[i]->SetYTitle("Entries");
    flistNoCutsFromB->Add(hMCd0D0ptNCfromBSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCfromBSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCfromBSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCfromBSB[i]->SetYTitle("Entries");
    flistNoCutsFromB->Add(hd0D0VtxTrueptNCfromBSB[i]);
  }



  //############ NO CUTS FROM DSTAR HISTOGRAMS ###########
  //
  //#############  global properties histos #######

  TH2F *hCPtaVSd0d0NCfromDstar=new TH2F("hCPtaVSd0d0NCfromDstar","hCPtaVSd0d0_NoCuts_FromDStar",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZNCfromDstar=new TH1F("hSecVtxZNCfromDstar","hSecVtxZ_NoCuts_FromDStar",1000,-8.,8.);
  TH1F *hSecVtxXNCfromDstar=new TH1F("hSecVtxXNCfromDstar","hSecVtxX_NoCuts_FromDStar",1000,-3000.,3000.);
  TH1F *hSecVtxYNCfromDstar=new TH1F("hSecVtxYNCfromDstar","hSecVtxY_NoCuts_FromDStar",1000,-3000.,3000.);
  TH2F *hSecVtxXYNCfromDstar=new TH2F("hSecVtxXYNCfromDstar","hSecVtxXY_NoCuts_FromDStar",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiNCfromDstar=new TH1F("hSecVtxPhiNCfromDstar","hSecVtxPhi_NoCuts_FromDStar",180,-180.1,180.1);
  TH1F *hCPtaNCfromDstar=new TH1F("hCPtaNCfromDstar","hCPta_NoCuts_FromDStar",100,0.,1.);
  TH1F *hd0xd0NCfromDstar=new TH1F("hd0xd0NCfromDstar","hd0xd0_NoCuts_FromDStar",1000,-100000.,100000.);
  TH1F *hMassTrueNCfromDstar=new TH1F("hMassTrueNCfromDstar","D^{0} MC inv. Mass No Cuts FromDStar(All momenta)",600,1.600,2.200);
  TH1F *hMassNCfromDstar=new TH1F("hMassNCfromDstar","D^{0} inv. Mass No Cuts FromDStar (All momenta)",600,1.600,2.200);
  hMassNCfromDstar->Sumw2();
  TH1F *hMassTrueNCfromDstarPM=new TH1F("hMassTrueNCfromDstarPM","D^{0} MC inv. Mass No Cuts FromDStar, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassNCfromDstarPM=new TH1F("hMassNCfromDstarPM","D^{0} inv. Mass No Cuts FromDStar (All momenta), MassPeak",600,1.600,2.200);
  hMassNCfromDstarPM->Sumw2();
  TH1F *hMassTrueNCfromDstarSB=new TH1F("hMassTrueNCfromDstarSB","D^{0} MC inv. Mass in Side Bands No Cuts FromDStar(All momenta)",600,1.600,2.200);
  TH1F *hMassNCfromDstarSB=new TH1F("hMassNCfromDstarSB","D^{0} inv. Mass in Side Bands No Cuts FromDStar (All momenta)",600,1.600,2.200);
  hMassNCfromDstarSB->Sumw2();

  flistNoCutsFromDstar->Add(hCPtaVSd0d0NCfromDstar);
  flistNoCutsFromDstar->Add(hSecVtxZNCfromDstar);
  flistNoCutsFromDstar->Add(hSecVtxYNCfromDstar);
  flistNoCutsFromDstar->Add(hSecVtxXNCfromDstar);
  flistNoCutsFromDstar->Add(hSecVtxXYNCfromDstar);
  flistNoCutsFromDstar->Add(hSecVtxPhiNCfromDstar);
  flistNoCutsFromDstar->Add(hCPtaNCfromDstar);
  flistNoCutsFromDstar->Add(hd0xd0NCfromDstar);
  flistNoCutsFromDstar->Add(hMassTrueNCfromDstar);
  flistNoCutsFromDstar->Add(hMassNCfromDstar);
  flistNoCutsFromDstar->Add(hMassTrueNCfromDstarPM);
  flistNoCutsFromDstar->Add(hMassNCfromDstarPM);
  flistNoCutsFromDstar->Add(hMassTrueNCfromDstarSB);
  flistNoCutsFromDstar->Add(hMassNCfromDstarSB);

  //########## d0 D0 histos #############  
  TH1F *hd0D0NCfromDstPM = new TH1F("hd0D0NCfromDstarPM","D^{0} impact par. plot , No Cuts ,FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0NCfromDstPM->SetXTitle("Impact parameter [#mum]");
  hd0D0NCfromDstPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueNCfromDstPM = new TH1F("hd0D0VtxTrueNCfromDstarPM","D^{0} impact par. w.r.t. True Vtx, No Cuts, FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueNCfromDstPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueNCfromDstPM->SetYTitle("Entries");

  TH1F *hMCd0D0NCfromDstPM = new TH1F("hMCd0D0NCfromDstarPM","D^{0} impact par. plot, No Cuts, FromDStar,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0NCfromDstPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0NCfromDstPM->SetYTitle("Entries");

  TH1F *hd0D0NCfromDstSB = new TH1F("hd0D0NCfromDstarSB","D^{0} impact par. plot , No Cuts ,FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0NCfromDstSB->SetXTitle("Impact parameter [#mum]");
  hd0D0NCfromDstSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueNCfromDstSB = new TH1F("hd0D0VtxTrueNCfromDstarSB","D^{0} impact par. w.r.t. True Vtx, No Cuts, FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueNCfromDstSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueNCfromDstSB->SetYTitle("Entries");

  TH1F *hMCd0D0NCfromDstSB = new TH1F("hMCd0D0NCfromDstarSB","D^{0} impact par. plot, No Cuts, FromDStar,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0NCfromDstSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0NCfromDstSB->SetYTitle("Entries");

  flistNoCutsFromDstar->Add(hd0D0NCfromDstPM);
  flistNoCutsFromDstar->Add(hd0D0VtxTrueNCfromDstPM);
  flistNoCutsFromDstar->Add(hMCd0D0NCfromDstPM);
  flistNoCutsFromDstar->Add(hd0D0NCfromDstSB);
  flistNoCutsFromDstar->Add(hd0D0VtxTrueNCfromDstSB);
  flistNoCutsFromDstar->Add(hMCd0D0NCfromDstSB);
  
  TH1F **hd0D0ptNCfromDstPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptNCfromDstPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptNCfromDstPM=new TH1F*[fnbins];
  TH1F **hd0D0ptNCfromDstSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptNCfromDstSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptNCfromDstSB=new TH1F*[fnbins];
  namehist="hd0D0ptNCfromDstar_";
  titlehist="D^{0} impact par. plot, No Cuts, FromDStar, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCfromDstPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCfromDstPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCfromDstPM[i]->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hd0D0ptNCfromDstPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCfromDstPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCfromDstPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCfromDstPM[i]->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hMCd0D0ptNCfromDstPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCfromDstPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCfromDstPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCfromDstPM[i]->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hd0D0VtxTrueptNCfromDstPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCfromDstSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCfromDstSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCfromDstSB[i]->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hd0D0ptNCfromDstSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCfromDstSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCfromDstSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCfromDstSB[i]->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hMCd0D0ptNCfromDstSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCfromDstSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCfromDstSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCfromDstSB[i]->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hd0D0VtxTrueptNCfromDstSB[i]);
  }


  //############ NO CUTS OTHER HISTOGRAMS ###########
  //
  //########### global properties histos ###########

  TH2F *hCPtaVSd0d0NCother=new TH2F("hCPtaVSd0d0NCother","hCPtaVSd0d0_NoCuts_other",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZNCother=new TH1F("hSecVtxZNCother","hSecVtxZ_NoCuts_other",1000,-8.,8.);
  TH1F *hSecVtxXNCother=new TH1F("hSecVtxXNCother","hSecVtxX_NoCuts_other",1000,-3000.,3000.);
  TH1F *hSecVtxYNCother=new TH1F("hSecVtxYNCother","hSecVtxY_NoCuts_other",1000,-3000.,3000.);
  TH2F *hSecVtxXYNCother=new TH2F("hSecVtxXYNCother","hSecVtxXY_NoCuts_other",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiNCother=new TH1F("hSecVtxPhiNCother","hSecVtxPhi_NoCuts_other",180,-180.1,180.1);
  TH1F *hCPtaNCother=new TH1F("hCPtaNCother","hCPta_NoCuts_other",100,0.,1.);
  TH1F *hd0xd0NCother=new TH1F("hd0xd0NCother","hd0xd0_NoCuts_other",1000,-100000.,100000.);
  TH1F *hMassTrueNCother=new TH1F("hMassTrueNCother","D^{0} MC inv. Mass No Cuts other(All momenta)",600,1.600,2.200);
  TH1F *hMassNCother=new TH1F("hMassNCother","D^{0} inv. Mass No Cuts other (All momenta)",600,1.600,2.200);
  hMassNCother->Sumw2();
  TH1F *hMassTrueNCotherPM=new TH1F("hMassTrueNCotherPM","D^{0} MC inv. Mass No Cuts Other, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassNCotherPM=new TH1F("hMassNCotherPM","D^{0} inv. Mass No Cuts Other (All momenta), MassPeak",600,1.600,2.200);
  hMassNCotherPM->Sumw2();
  TH1F *hMassTrueNCotherSB=new TH1F("hMassTrueNCotherSB","D^{0} MC inv. Mass in Side Bands No Cuts other(All momenta)",600,1.600,2.200);
  TH1F *hMassNCotherSB=new TH1F("hMassNCotherSB","D^{0} inv. Mass in Side Bands No Cuts other (All momenta)",600,1.600,2.200);
  hMassNCotherSB->Sumw2();

  flistNoCutsOther->Add(hCPtaVSd0d0NCother);
  flistNoCutsOther->Add(hSecVtxZNCother);
  flistNoCutsOther->Add(hSecVtxYNCother);
  flistNoCutsOther->Add(hSecVtxXNCother);
  flistNoCutsOther->Add(hSecVtxXYNCother);
  flistNoCutsOther->Add(hSecVtxPhiNCother);
  flistNoCutsOther->Add(hCPtaNCother);
  flistNoCutsOther->Add(hd0xd0NCother);
  flistNoCutsOther->Add(hMassTrueNCother);
  flistNoCutsOther->Add(hMassNCother);
  flistNoCutsOther->Add(hMassTrueNCotherPM);
  flistNoCutsOther->Add(hMassNCotherPM);
  flistNoCutsOther->Add(hMassTrueNCotherSB);
  flistNoCutsOther->Add(hMassNCotherSB);

  //############# d0 D0 histos ###############à
  TH1F *hd0D0NCotherPM = new TH1F("hd0D0NCotherPM","D^{0} impact par. plot , No Cuts ,Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0NCotherPM->SetXTitle("Impact parameter [#mum]");
  hd0D0NCotherPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueNCotherPM = new TH1F("hd0D0VtxTrueNCotherPM","D^{0} impact par. w.r.t. True Vtx, No Cuts, Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueNCotherPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueNCotherPM->SetYTitle("Entries");

  TH1F *hMCd0D0NCotherPM = new TH1F("hMCd0D0NCotherPM","D^{0} impact par. plot, No Cuts, Other,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0NCotherPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0NCotherPM->SetYTitle("Entries");

  TH1F *hd0D0NCotherSB = new TH1F("hd0D0NCotherSB","D^{0} impact par. plot , No Cuts ,Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0NCotherSB->SetXTitle("Impact parameter [#mum]");
  hd0D0NCotherSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueNCotherSB = new TH1F("hd0D0VtxTrueNCotherSB","D^{0} impact par. w.r.t. True Vtx, No Cuts, Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueNCotherSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueNCotherSB->SetYTitle("Entries");

  TH1F *hMCd0D0NCotherSB = new TH1F("hMCd0D0NCotherSB","D^{0} impact par. plot, No Cuts, Other,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0NCotherSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0NCotherSB->SetYTitle("Entries");

  flistNoCutsOther->Add(hd0D0NCotherPM);
  flistNoCutsOther->Add(hd0D0VtxTrueNCotherPM);
  flistNoCutsOther->Add(hMCd0D0NCotherPM);
  flistNoCutsOther->Add(hd0D0NCotherSB);
  flistNoCutsOther->Add(hd0D0VtxTrueNCotherSB);
  flistNoCutsOther->Add(hMCd0D0NCotherSB);
  
  TH1F **hd0D0ptNCotherPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptNCotherPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptNCotherPM=new TH1F*[fnbins];
  TH1F **hd0D0ptNCotherSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptNCotherSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptNCotherSB=new TH1F*[fnbins];
  namehist="hd0D0ptNCother_";
  titlehist="D^{0} impact par. plot, No Cuts, Other, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCotherPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCotherPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCotherPM[i]->SetYTitle("Entries");
    flistNoCutsOther->Add(hd0D0ptNCotherPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCotherPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCotherPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCotherPM[i]->SetYTitle("Entries");
    flistNoCutsOther->Add(hMCd0D0ptNCotherPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCotherPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCotherPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCotherPM[i]->SetYTitle("Entries");
    flistNoCutsOther->Add(hd0D0VtxTrueptNCotherPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCotherSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCotherSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCotherSB[i]->SetYTitle("Entries");
    flistNoCutsOther->Add(hd0D0ptNCotherSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCotherSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCotherSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCotherSB[i]->SetYTitle("Entries");
    flistNoCutsOther->Add(hMCd0D0ptNCotherSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCotherSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCotherSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCotherSB[i]->SetYTitle("Entries");
    flistNoCutsOther->Add(hd0D0VtxTrueptNCotherSB[i]);
  }


  //################################################################################################
  //                                                                                               #
  //                         HISTOS FOR LOOSE CUTS                                                 #
  //                                                                                               #
  //################################################################################################

  //############ LOOSE CUTS SIGNAL HISTOGRAMS ###############
  //
  // ####### global properties histo ############

  TH2F *hCPtaVSd0d0LSCsign=new TH2F("hCPtaVSd0d0LSCsign","hCPtaVSd0d0_LooseCuts_Signal",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZLSCsign=new TH1F("hSecVtxZLSCsign","hSecVtxZ_LooseCuts_Signal",1000,-8.,8.);
  TH1F *hSecVtxXLSCsign=new TH1F("hSecVtxXLSCsign","hSecVtxX_LooseCuts_Signal",1000,-3000.,3000.);
  TH1F *hSecVtxYLSCsign=new TH1F("hSecVtxYLSCsign","hSecVtxY_LooseCuts_Signal",1000,-3000.,3000.);
  TH2F *hSecVtxXYLSCsign=new TH2F("hSecVtxXYLSCsign","hSecVtxXY_LooseCuts_Signal",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiLSCsign=new TH1F("hSecVtxPhiLSCsign","hSecVtxPhi_LooseCuts_Signal",180,-180.1,180.1);
  TH1F *hCPtaLSCsign=new TH1F("hCPtaLSCsign","hCPta_LooseCuts_Signal",100,0.,1.);
  TH1F *hd0xd0LSCsign=new TH1F("hd0xd0LSCsign","hd0xd0_LooseCuts_Signal",1000,-100000.,100000.);
  TH1F *hMassTrueLSCsign=new TH1F("hMassTrueLSCsign","D^{0} MC inv. Mass Loose Cuts Signal(All momenta)",600,1.600,2.200);
  TH1F *hMassLSCsign=new TH1F("hMassLSCsign","D^{0} inv. Mass Loose Cuts Signal (All momenta)",600,1.600,2.200);
  hMassLSCsign->Sumw2();
  TH1F *hMassTrueLSCsignPM=new TH1F("hMassTrueLSCsignPM","D^{0} MC inv. Mass Loose Cuts Signal, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassLSCsignPM=new TH1F("hMassLSCsignPM","D^{0} inv. Mass Loose Cuts Signal (All momenta), MassPeak",600,1.600,2.200);
  hMassLSCsignPM->Sumw2();
  TH1F *hMassTrueLSCsignSB=new TH1F("hMassTrueLSCsignSB","D^{0} MC inv. Mass in Side Bands Loose Cuts Signal(All momenta)",600,1.600,2.200);
  TH1F *hMassLSCsignSB=new TH1F("hMassLSCsignSB","D^{0} inv. Mass in Side Bands Loose Cuts Signal (All momenta)",600,1.600,2.200);
  hMassLSCsignSB->Sumw2();

  flistLsCutsSignal->Add(hCPtaVSd0d0LSCsign);
  flistLsCutsSignal->Add(hSecVtxZLSCsign);
  flistLsCutsSignal->Add(hSecVtxYLSCsign);
  flistLsCutsSignal->Add(hSecVtxXLSCsign);
  flistLsCutsSignal->Add(hSecVtxXYLSCsign);
  flistLsCutsSignal->Add(hSecVtxPhiLSCsign);
  flistLsCutsSignal->Add(hCPtaLSCsign);
  flistLsCutsSignal->Add(hd0xd0LSCsign);
  flistLsCutsSignal->Add(hMassTrueLSCsign);
  flistLsCutsSignal->Add(hMassLSCsign);
  flistLsCutsSignal->Add(hMassTrueLSCsignPM);
  flistLsCutsSignal->Add(hMassLSCsignPM);
  flistLsCutsSignal->Add(hMassTrueLSCsignSB);
  flistLsCutsSignal->Add(hMassLSCsignSB);

  // ####### d0 D0 histos ############
  TH1F *hd0D0LSCsignPM = new TH1F("hd0D0LSCsignPM","D^{0} impact par. plot , Loose Cuts ,Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0LSCsignPM->SetXTitle("Impact parameter [#mum]");
  hd0D0LSCsignPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueLSCsignPM = new TH1F("hd0D0VtxTrueLSCsignPM","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueLSCsignPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueLSCsignPM->SetYTitle("Entries");

  TH1F *hMCd0D0LSCsignPM = new TH1F("hMCd0D0LSCsignPM","D^{0} impact par. plot, Loose Cuts, Signal,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0LSCsignPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0LSCsignPM->SetYTitle("Entries");

  TH1F *hd0D0LSCsignSB = new TH1F("hd0D0LSCsignSB","D^{0} impact par. plot , Loose Cuts ,Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0LSCsignSB->SetXTitle("Impact parameter [#mum]");
  hd0D0LSCsignSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueLSCsignSB = new TH1F("hd0D0VtxTrueLSCsignSB","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueLSCsignSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueLSCsignSB->SetYTitle("Entries");

  TH1F *hMCd0D0LSCsignSB = new TH1F("hMCd0D0LSCsignSB","D^{0} impact par. plot, Loose Cuts, Signal,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0LSCsignSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0LSCsignSB->SetYTitle("Entries");

  flistLsCutsSignal->Add(hd0D0LSCsignPM);
  flistLsCutsSignal->Add(hd0D0VtxTrueLSCsignPM);
  flistLsCutsSignal->Add(hMCd0D0LSCsignPM);
  flistLsCutsSignal->Add(hd0D0LSCsignSB);
  flistLsCutsSignal->Add(hd0D0VtxTrueLSCsignSB);
  flistLsCutsSignal->Add(hMCd0D0LSCsignSB);
  
  TH1F **hd0D0ptLSCsignPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptLSCsignPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptLSCsignPM=new TH1F*[fnbins];
  TH1F **hd0D0ptLSCsignSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptLSCsignSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptLSCsignSB=new TH1F*[fnbins];
  namehist="hd0D0ptLSCsign_";
  titlehist="D^{0} impact par. plot, Loose Cuts, Signal, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCsignPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCsignPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCsignPM[i]->SetYTitle("Entries");
    flistLsCutsSignal->Add(hd0D0ptLSCsignPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCsignPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCsignPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCsignPM[i]->SetYTitle("Entries");
    flistLsCutsSignal->Add(hMCd0D0ptLSCsignPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCsignPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCsignPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCsignPM[i]->SetYTitle("Entries");
    flistLsCutsSignal->Add(hd0D0VtxTrueptLSCsignPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCsignSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCsignSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCsignSB[i]->SetYTitle("Entries");
    flistLsCutsSignal->Add(hd0D0ptLSCsignSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCsignSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCsignSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCsignSB[i]->SetYTitle("Entries");
    flistLsCutsSignal->Add(hMCd0D0ptLSCsignSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCsignSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCsignSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCsignSB[i]->SetYTitle("Entries");
    flistLsCutsSignal->Add(hd0D0VtxTrueptLSCsignSB[i]);
  }


  //############ LOOSE CUTS BACKGROUND HISTOGRAMS ###########
  //
  //   ######## global properties histos #######
  TH2F *hCPtaVSd0d0LSCback=new TH2F("hCPtaVSd0d0LSCback","hCPtaVSd0d0_LooseCuts_Background",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZLSCback=new TH1F("hSecVtxZLSCback","hSecVtxZ_LooseCuts_Background",1000,-8.,8.);
  TH1F *hSecVtxXLSCback=new TH1F("hSecVtxXLSCback","hSecVtxX_LooseCuts_Background",1000,-3000.,3000.);
  TH1F *hSecVtxYLSCback=new TH1F("hSecVtxYLSCback","hSecVtxY_LooseCuts_Background",1000,-3000.,3000.);
  TH2F *hSecVtxXYLSCback=new TH2F("hSecVtxXYLSCback","hSecVtxXY_LooseCuts_Background",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiLSCback=new TH1F("hSecVtxPhiLSCback","hSecVtxPhi_LooseCuts_Background",180,-180.1,180.1);
  TH1F *hCPtaLSCback=new TH1F("hCPtaLSCback","hCPta_LooseCuts_Background",100,0.,1.);
  TH1F *hd0xd0LSCback=new TH1F("hd0xd0LSCback","hd0xd0_LooseCuts_Background",1000,-100000.,100000.);
  TH1F *hMassTrueLSCback=new TH1F("hMassTrueLSCback","D^{0} MC inv. Mass Loose Cuts Background(All momenta)",600,1.600,2.200);
  TH1F *hMassLSCback=new TH1F("hMassLSCback","D^{0} inv. Mass Loose Cuts Background (All momenta)",600,1.600,2.200);
  hMassLSCback->Sumw2();
  TH1F *hMassTrueLSCbackPM=new TH1F("hMassTrueLSCbackPM","D^{0} MC inv. Mass Loose Cuts Background, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassLSCbackPM=new TH1F("hMassLSCbackPM","D^{0} inv. Mass Loose Cuts Background (All momenta), MassPeak",600,1.600,2.200);
  hMassLSCbackPM->Sumw2();
  TH1F *hMassTrueLSCbackSB=new TH1F("hMassTrueLSCbackSB","D^{0} MC inv. Mass in Side Bands Loose Cuts Background(All momenta)",600,1.600,2.200);
  TH1F *hMassLSCbackSB=new TH1F("hMassLSCbackSB","D^{0} inv. Mass in Side Bands Loose Cuts Background (All momenta)",600,1.600,2.200);
  hMassLSCbackSB->Sumw2();

  flistLsCutsBack->Add(hCPtaVSd0d0LSCback);
  flistLsCutsBack->Add(hSecVtxZLSCback);
  flistLsCutsBack->Add(hSecVtxYLSCback);
  flistLsCutsBack->Add(hSecVtxXLSCback);
  flistLsCutsBack->Add(hSecVtxXYLSCback);
  flistLsCutsBack->Add(hSecVtxPhiLSCback);
  flistLsCutsBack->Add(hCPtaLSCback);
  flistLsCutsBack->Add(hd0xd0LSCback);
  flistLsCutsBack->Add(hMassTrueLSCback);
  flistLsCutsBack->Add(hMassLSCback);
  flistLsCutsBack->Add(hMassTrueLSCbackPM);
  flistLsCutsBack->Add(hMassLSCbackPM);
  flistLsCutsBack->Add(hMassTrueLSCbackSB);
  flistLsCutsBack->Add(hMassLSCbackSB);


  // ####### d0 D0 histos ############
  
 TH1F *hd0D0LSCbackPM = new TH1F("hd0D0LSCbackPM","D^{0} impact par. plot , Loose Cuts ,Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0LSCbackPM->SetXTitle("Impact parameter [#mum]");
  hd0D0LSCbackPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueLSCbackPM = new TH1F("hd0D0VtxTrueLSCbackPM","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueLSCbackPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueLSCbackPM->SetYTitle("Entries");

  TH1F *hMCd0D0LSCbackPM = new TH1F("hMCd0D0LSCbackPM","D^{0} impact par. plot, Loose Cuts, Background,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0LSCbackPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0LSCbackPM->SetYTitle("Entries");

  TH1F *hd0D0LSCbackSB = new TH1F("hd0D0LSCbackSB","D^{0} impact par. plot , Loose Cuts ,Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0LSCbackSB->SetXTitle("Impact parameter [#mum]");
  hd0D0LSCbackSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueLSCbackSB = new TH1F("hd0D0VtxTrueLSCbackSB","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueLSCbackSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueLSCbackSB->SetYTitle("Entries");

  TH1F *hMCd0D0LSCbackSB = new TH1F("hMCd0D0LSCbackSB","D^{0} impact par. plot, Loose Cuts, Background,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0LSCbackSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0LSCbackSB->SetYTitle("Entries");

  flistLsCutsBack->Add(hd0D0LSCbackPM);
  flistLsCutsBack->Add(hd0D0VtxTrueLSCbackPM);
  flistLsCutsBack->Add(hMCd0D0LSCbackPM);
  flistLsCutsBack->Add(hd0D0LSCbackSB);
  flistLsCutsBack->Add(hd0D0VtxTrueLSCbackSB);
  flistLsCutsBack->Add(hMCd0D0LSCbackSB);
  
  TH1F **hd0D0ptLSCbackPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptLSCbackPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptLSCbackPM=new TH1F*[fnbins];
  TH1F **hd0D0ptLSCbackSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptLSCbackSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptLSCbackSB=new TH1F*[fnbins];
  namehist="hd0D0ptLSCback_";
  titlehist="D^{0} impact par. plot, Loose Cuts, Background, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCbackPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCbackPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCbackPM[i]->SetYTitle("Entries");
    flistLsCutsBack->Add(hd0D0ptLSCbackPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCbackPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCbackPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCbackPM[i]->SetYTitle("Entries");
    flistLsCutsBack->Add(hMCd0D0ptLSCbackPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCbackPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCbackPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCbackPM[i]->SetYTitle("Entries");
    flistLsCutsBack->Add(hd0D0VtxTrueptLSCbackPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCbackSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCbackSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCbackSB[i]->SetYTitle("Entries");
    flistLsCutsBack->Add(hd0D0ptLSCbackSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCbackSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCbackSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCbackSB[i]->SetYTitle("Entries");
    flistLsCutsBack->Add(hMCd0D0ptLSCbackSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCbackSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCbackSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCbackSB[i]->SetYTitle("Entries");
    flistLsCutsBack->Add(hd0D0VtxTrueptLSCbackSB[i]);
  }



 //############ LOOSE CUTS FROMB HISTOGRAMS ###########
  //
  //#######  global properties histos

  TH2F *hCPtaVSd0d0LSCfromB=new TH2F("hCPtaVSd0d0LSCfromB","hCPtaVSd0d0_LooseCuts_FromB",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZLSCfromB=new TH1F("hSecVtxZLSCfromB","hSecVtxZ_LooseCuts_FromB",1000,-8.,8.);
  TH1F *hSecVtxXLSCfromB=new TH1F("hSecVtxXLSCfromB","hSecVtxX_LooseCuts_FromB",1000,-3000.,3000.);
  TH1F *hSecVtxYLSCfromB=new TH1F("hSecVtxYLSCfromB","hSecVtxY_LooseCuts_FromB",1000,-3000.,3000.);
  TH2F *hSecVtxXYLSCfromB=new TH2F("hSecVtxXYLSCfromB","hSecVtxXY_LooseCuts_FromB",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiLSCfromB=new TH1F("hSecVtxPhiLSCfromB","hSecVtxPhi_LooseCuts_FromB",180,-180.1,180.1);
  TH1F *hCPtaLSCfromB=new TH1F("hCPtaLSCfromB","hCPta_LooseCuts_FromB",100,0.,1.);
  TH1F *hd0xd0LSCfromB=new TH1F("hd0xd0LSCfromB","hd0xd0_LooseCuts_FromB",1000,-100000.,100000.);
  TH1F *hMassTrueLSCfromB=new TH1F("hMassTrueLSCfromB","D^{0} MC inv. Mass Loose Cuts FromB(All momenta)",600,1.600,2.200);
  TH1F *hMassLSCfromB=new TH1F("hMassLSCfromB","D^{0} inv. Mass Loose Cuts FromB (All momenta)",600,1.600,2.200);
  hMassLSCfromB->Sumw2();
  TH1F *hMassTrueLSCfromBPM=new TH1F("hMassTrueLSCfromBPM","D^{0} MC inv. Mass Loose Cuts FromB, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassLSCfromBPM=new TH1F("hMassLSCfromBPM","D^{0} inv. Mass Loose Cuts FromB (All momenta), MassPeak",600,1.600,2.200);
  hMassLSCfromBPM->Sumw2();
  TH1F *hMassTrueLSCfromBSB=new TH1F("hMassTrueLSCfromBSB","D^{0} MC inv. Mass in Side Bands Loose Cuts FromB(All momenta)",600,1.600,2.200);
  TH1F *hMassLSCfromBSB=new TH1F("hMassLSCfromBSB","D^{0} inv. Mass in Side Bands Loose Cuts FromB (All momenta)",600,1.600,2.200);
  hMassLSCfromBSB->Sumw2();

  flistLsCutsFromB->Add(hCPtaVSd0d0LSCfromB);
  flistLsCutsFromB->Add(hSecVtxZLSCfromB);
  flistLsCutsFromB->Add(hSecVtxYLSCfromB);
  flistLsCutsFromB->Add(hSecVtxXLSCfromB);
  flistLsCutsFromB->Add(hSecVtxXYLSCfromB);
  flistLsCutsFromB->Add(hSecVtxPhiLSCfromB);
  flistLsCutsFromB->Add(hCPtaLSCfromB);
  flistLsCutsFromB->Add(hd0xd0LSCfromB);
  flistLsCutsFromB->Add(hMassTrueLSCfromB);
  flistLsCutsFromB->Add(hMassLSCfromB);
  flistLsCutsFromB->Add(hMassTrueLSCfromBPM);
  flistLsCutsFromB->Add(hMassLSCfromBPM);
  flistLsCutsFromB->Add(hMassTrueLSCfromBSB);
  flistLsCutsFromB->Add(hMassLSCfromBSB);

  // ######### d0 D0 histos ##############
  TH1F *hd0D0LSCfromBPM = new TH1F("hd0D0LSCfromBPM","D^{0} impact par. plot , Loose Cuts ,FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0LSCfromBPM->SetXTitle("Impact parameter [#mum]");
  hd0D0LSCfromBPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueLSCfromBPM = new TH1F("hd0D0VtxTrueLSCfromBPM","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueLSCfromBPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueLSCfromBPM->SetYTitle("Entries");

  TH1F *hMCd0D0LSCfromBPM = new TH1F("hMCd0D0LSCfromBPM","D^{0} impact par. plot, Loose Cuts, FromB,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0LSCfromBPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0LSCfromBPM->SetYTitle("Entries");

  TH1F *hd0D0LSCfromBSB = new TH1F("hd0D0LSCfromBSB","D^{0} impact par. plot , Loose Cuts ,FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0LSCfromBSB->SetXTitle("Impact parameter [#mum]");
  hd0D0LSCfromBSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueLSCfromBSB = new TH1F("hd0D0VtxTrueLSCfromBSB","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueLSCfromBSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueLSCfromBSB->SetYTitle("Entries");

  TH1F *hMCd0D0LSCfromBSB = new TH1F("hMCd0D0LSCfromBSB","D^{0} impact par. plot, Loose Cuts, FromB,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0LSCfromBSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0LSCfromBSB->SetYTitle("Entries");

  flistLsCutsFromB->Add(hd0D0LSCfromBPM);
  flistLsCutsFromB->Add(hd0D0VtxTrueLSCfromBPM);
  flistLsCutsFromB->Add(hMCd0D0LSCfromBPM);
  flistLsCutsFromB->Add(hd0D0LSCfromBSB);
  flistLsCutsFromB->Add(hd0D0VtxTrueLSCfromBSB);
  flistLsCutsFromB->Add(hMCd0D0LSCfromBSB);
  
  TH1F **hd0D0ptLSCfromBPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptLSCfromBPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptLSCfromBPM=new TH1F*[fnbins];
  TH1F **hd0D0ptLSCfromBSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptLSCfromBSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptLSCfromBSB=new TH1F*[fnbins];
  namehist="hd0D0ptLSCfromB_";
  titlehist="D^{0} impact par. plot, Loose Cuts, FromB, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCfromBPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCfromBPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCfromBPM[i]->SetYTitle("Entries");
    flistLsCutsFromB->Add(hd0D0ptLSCfromBPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCfromBPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCfromBPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCfromBPM[i]->SetYTitle("Entries");
    flistLsCutsFromB->Add(hMCd0D0ptLSCfromBPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCfromBPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCfromBPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCfromBPM[i]->SetYTitle("Entries");
    flistLsCutsFromB->Add(hd0D0VtxTrueptLSCfromBPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCfromBSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCfromBSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCfromBSB[i]->SetYTitle("Entries");
    flistLsCutsFromB->Add(hd0D0ptLSCfromBSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCfromBSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCfromBSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCfromBSB[i]->SetYTitle("Entries");
    flistLsCutsFromB->Add(hMCd0D0ptLSCfromBSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCfromBSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCfromBSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCfromBSB[i]->SetYTitle("Entries");
    flistLsCutsFromB->Add(hd0D0VtxTrueptLSCfromBSB[i]);
  }



 //############ LOOSE CUTS FROM DSTAR HISTOGRAMS ###########
 //
  //############## global properties histos
  TH2F *hCPtaVSd0d0LSCfromDstar=new TH2F("hCPtaVSd0d0LSCfromDstar","hCPtaVSd0d0_LooseCuts_FromDStar",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZLSCfromDstar=new TH1F("hSecVtxZLSCfromDstar","hSecVtxZ_LooseCuts_FromDStar",1000,-8.,8.);
  TH1F *hSecVtxXLSCfromDstar=new TH1F("hSecVtxXLSCfromDstar","hSecVtxX_LooseCuts_FromDStar",1000,-3000.,3000.);
  TH1F *hSecVtxYLSCfromDstar=new TH1F("hSecVtxYLSCfromDstar","hSecVtxY_LooseCuts_FromDStar",1000,-3000.,3000.);
  TH2F *hSecVtxXYLSCfromDstar=new TH2F("hSecVtxXYLSCfromDstar","hSecVtxXY_LooseCuts_FromDStar",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiLSCfromDstar=new TH1F("hSecVtxPhiLSCfromDstar","hSecVtxPhi_LooseCuts_FromDStar",180,-180.1,180.1);
  TH1F *hCPtaLSCfromDstar=new TH1F("hCPtaLSCfromDstar","hCPta_LooseCuts_FromDStar",100,0.,1.);
  TH1F *hd0xd0LSCfromDstar=new TH1F("hd0xd0LSCfromDstar","hd0xd0_LooseCuts_FromDStar",1000,-100000.,100000.);
  TH1F *hMassTrueLSCfromDstar=new TH1F("hMassTrueLSCfromDstar","D^{0} MC inv. Mass Loose Cuts FromDStar(All momenta)",600,1.600,2.200);
  TH1F *hMassLSCfromDstar=new TH1F("hMassLSCfromDstar","D^{0} inv. Mass Loose Cuts FromDStar (All momenta)",600,1.600,2.200);
  hMassLSCfromDstar->Sumw2();
  TH1F *hMassTrueLSCfromDstarPM=new TH1F("hMassTrueLSCfromDstarPM","D^{0} MC inv. Mass Loose Cuts FromDStar, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassLSCfromDstarPM=new TH1F("hMassLSCfromDstarPM","D^{0} inv. Mass Loose Cuts FromDStar (All momenta), MassPeak",600,1.600,2.200);
  hMassLSCfromDstarPM->Sumw2();
  TH1F *hMassTrueLSCfromDstarSB=new TH1F("hMassTrueLSCfromDstarSB","D^{0} MC inv. Mass in Side Bands Loose Cuts FromDStar(All momenta)",600,1.600,2.200);
  TH1F *hMassLSCfromDstarSB=new TH1F("hMassLSCfromDstarSB","D^{0} inv. Mass in Side Bands Loose Cuts FromDStar (All momenta)",600,1.600,2.200);
  hMassLSCfromDstarSB->Sumw2();

  flistLsCutsFromDstar->Add(hCPtaVSd0d0LSCfromDstar);
  flistLsCutsFromDstar->Add(hSecVtxZLSCfromDstar);
  flistLsCutsFromDstar->Add(hSecVtxYLSCfromDstar);
  flistLsCutsFromDstar->Add(hSecVtxXLSCfromDstar);
  flistLsCutsFromDstar->Add(hSecVtxXYLSCfromDstar);
  flistLsCutsFromDstar->Add(hSecVtxPhiLSCfromDstar);
  flistLsCutsFromDstar->Add(hCPtaLSCfromDstar);
  flistLsCutsFromDstar->Add(hd0xd0LSCfromDstar);
  flistLsCutsFromDstar->Add(hMassTrueLSCfromDstar);
  flistLsCutsFromDstar->Add(hMassLSCfromDstar);
 flistLsCutsFromDstar->Add(hMassTrueLSCfromDstarPM);
  flistLsCutsFromDstar->Add(hMassLSCfromDstarPM);
  flistLsCutsFromDstar->Add(hMassTrueLSCfromDstarSB);
  flistLsCutsFromDstar->Add(hMassLSCfromDstarSB);

  //########## d0 D0 histos #############  
  TH1F *hd0D0LSCfromDstPM = new TH1F("hd0D0LSCfromDstarPM","D^{0} impact par. plot , Loose Cuts ,FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0LSCfromDstPM->SetXTitle("Impact parameter [#mum]");
  hd0D0LSCfromDstPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueLSCfromDstPM = new TH1F("hd0D0VtxTrueLSCfromDstarPM","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueLSCfromDstPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueLSCfromDstPM->SetYTitle("Entries");

  TH1F *hMCd0D0LSCfromDstPM = new TH1F("hMCd0D0LSCfromDstarPM","D^{0} impact par. plot, Loose Cuts, FromDStar,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0LSCfromDstPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0LSCfromDstPM->SetYTitle("Entries");

  TH1F *hd0D0LSCfromDstSB = new TH1F("hd0D0LSCfromDstarSB","D^{0} impact par. plot , Loose Cuts ,FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0LSCfromDstSB->SetXTitle("Impact parameter [#mum]");
  hd0D0LSCfromDstSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueLSCfromDstSB = new TH1F("hd0D0VtxTrueLSCfromDstarSB","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueLSCfromDstSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueLSCfromDstSB->SetYTitle("Entries");

  TH1F *hMCd0D0LSCfromDstSB = new TH1F("hMCd0D0LSCfromDstarSB","D^{0} impact par. plot, Loose Cuts, FromDStar,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0LSCfromDstSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0LSCfromDstSB->SetYTitle("Entries");

  flistLsCutsFromDstar->Add(hd0D0LSCfromDstPM);
  flistLsCutsFromDstar->Add(hd0D0VtxTrueLSCfromDstPM);
  flistLsCutsFromDstar->Add(hMCd0D0LSCfromDstPM);
  flistLsCutsFromDstar->Add(hd0D0LSCfromDstSB);
  flistLsCutsFromDstar->Add(hd0D0VtxTrueLSCfromDstSB);
  flistLsCutsFromDstar->Add(hMCd0D0LSCfromDstSB);
  
  TH1F **hd0D0ptLSCfromDstPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptLSCfromDstPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptLSCfromDstPM=new TH1F*[fnbins];
  TH1F **hd0D0ptLSCfromDstSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptLSCfromDstSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptLSCfromDstSB=new TH1F*[fnbins];
  namehist="hd0D0ptLSCfromDstar_";
  titlehist="D^{0} impact par. plot, Loose Cuts, FromDStar, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCfromDstPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCfromDstPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCfromDstPM[i]->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hd0D0ptLSCfromDstPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCfromDstPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCfromDstPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCfromDstPM[i]->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hMCd0D0ptLSCfromDstPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCfromDstPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCfromDstPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCfromDstPM[i]->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hd0D0VtxTrueptLSCfromDstPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCfromDstSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCfromDstSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCfromDstSB[i]->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hd0D0ptLSCfromDstSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCfromDstSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCfromDstSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCfromDstSB[i]->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hMCd0D0ptLSCfromDstSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCfromDstSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCfromDstSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCfromDstSB[i]->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hd0D0VtxTrueptLSCfromDstSB[i]);
  }


  //############ LOOSE CUTS OTHER HISTOGRAMS ###########
  //
  //########### global properties histos ###########

  TH2F *hCPtaVSd0d0LSCother=new TH2F("hCPtaVSd0d0LSCother","hCPtaVSd0d0_LooseCuts_other",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZLSCother=new TH1F("hSecVtxZLSCother","hSecVtxZ_LooseCuts_other",1000,-8.,8.);
  TH1F *hSecVtxXLSCother=new TH1F("hSecVtxXLSCother","hSecVtxX_LooseCuts_other",1000,-3000.,3000.);
  TH1F *hSecVtxYLSCother=new TH1F("hSecVtxYLSCother","hSecVtxY_LooseCuts_other",1000,-3000.,3000.);
  TH2F *hSecVtxXYLSCother=new TH2F("hSecVtxXYLSCother","hSecVtxXY_LooseCuts_other",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiLSCother=new TH1F("hSecVtxPhiLSCother","hSecVtxPhi_LooseCuts_other",180,-180.1,180.1);
  TH1F *hCPtaLSCother=new TH1F("hCPtaLSCother","hCPta_LooseCuts_other",100,0.,1.);
  TH1F *hd0xd0LSCother=new TH1F("hd0xd0LSCother","hd0xd0_LooseCuts_other",1000,-100000.,100000.);
  TH1F *hMassTrueLSCother=new TH1F("hMassTrueLSCother","D^{0} MC inv. Mass Loose Cuts other(All momenta)",600,1.600,2.200);
  TH1F *hMassLSCother=new TH1F("hMassLSCother","D^{0} inv. Mass Loose Cuts other (All momenta)",600,1.600,2.200);
  hMassLSCother->Sumw2();
  TH1F *hMassTrueLSCotherPM=new TH1F("hMassTrueLSCotherPM","D^{0} MC inv. Mass Loose Cuts other, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassLSCotherPM=new TH1F("hMassLSCotherPM","D^{0} inv. Mass Loose Cuts other (All momenta), MassPeak",600,1.600,2.200);
  hMassLSCotherPM->Sumw2();
  TH1F *hMassTrueLSCotherSB=new TH1F("hMassTrueLSCotherSB","D^{0} MC inv. Mass in Side Bands Loose Cuts other(All momenta)",600,1.600,2.200);
  TH1F *hMassLSCotherSB=new TH1F("hMassLSCotherSB","D^{0} inv. Mass in Side Bands Loose Cuts other (All momenta)",600,1.600,2.200);
  hMassLSCotherSB->Sumw2();

  flistLsCutsOther->Add(hCPtaVSd0d0LSCother);
  flistLsCutsOther->Add(hSecVtxZLSCother);
  flistLsCutsOther->Add(hSecVtxYLSCother);
  flistLsCutsOther->Add(hSecVtxXLSCother);
  flistLsCutsOther->Add(hSecVtxXYLSCother);
  flistLsCutsOther->Add(hSecVtxPhiLSCother);
  flistLsCutsOther->Add(hCPtaLSCother);
  flistLsCutsOther->Add(hd0xd0LSCother);
  flistLsCutsOther->Add(hMassTrueLSCother);
  flistLsCutsOther->Add(hMassLSCother);
  flistLsCutsOther->Add(hMassTrueLSCotherPM);
  flistLsCutsOther->Add(hMassLSCotherPM);
  flistLsCutsOther->Add(hMassTrueLSCotherSB);
  flistLsCutsOther->Add(hMassLSCotherSB);

  //############# d0 D0 histos ###############à
  TH1F *hd0D0LSCotherPM = new TH1F("hd0D0LSCotherPM","D^{0} impact par. plot , Loose Cuts ,Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0LSCotherPM->SetXTitle("Impact parameter [#mum]");
  hd0D0LSCotherPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueLSCotherPM = new TH1F("hd0D0VtxTrueLSCotherPM","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueLSCotherPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueLSCotherPM->SetYTitle("Entries");

  TH1F *hMCd0D0LSCotherPM = new TH1F("hMCd0D0LSCotherPM","D^{0} impact par. plot, Loose Cuts, Other,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0LSCotherPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0LSCotherPM->SetYTitle("Entries");

  TH1F *hd0D0LSCotherSB = new TH1F("hd0D0LSCotherSB","D^{0} impact par. plot , Loose Cuts ,Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0LSCotherSB->SetXTitle("Impact parameter [#mum]");
  hd0D0LSCotherSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueLSCotherSB = new TH1F("hd0D0VtxTrueLSCotherSB","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueLSCotherSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueLSCotherSB->SetYTitle("Entries");

  TH1F *hMCd0D0LSCotherSB = new TH1F("hMCd0D0LSCotherSB","D^{0} impact par. plot, Loose Cuts, Other,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0LSCotherSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0LSCotherSB->SetYTitle("Entries");

  flistLsCutsOther->Add(hd0D0LSCotherPM);
  flistLsCutsOther->Add(hd0D0VtxTrueLSCotherPM);
  flistLsCutsOther->Add(hMCd0D0LSCotherPM);
  flistLsCutsOther->Add(hd0D0LSCotherSB);
  flistLsCutsOther->Add(hd0D0VtxTrueLSCotherSB);
  flistLsCutsOther->Add(hMCd0D0LSCotherSB);
  
  TH1F **hd0D0ptLSCotherPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptLSCotherPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptLSCotherPM=new TH1F*[fnbins];
  TH1F **hd0D0ptLSCotherSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptLSCotherSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptLSCotherSB=new TH1F*[fnbins];
  namehist="hd0D0ptLSCother_";
  titlehist="D^{0} impact par. plot, Loose Cuts, Other, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCotherPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCotherPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCotherPM[i]->SetYTitle("Entries");
    flistLsCutsOther->Add(hd0D0ptLSCotherPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCotherPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCotherPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCotherPM[i]->SetYTitle("Entries");
    flistLsCutsOther->Add(hMCd0D0ptLSCotherPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCotherPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCotherPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCotherPM[i]->SetYTitle("Entries");
    flistLsCutsOther->Add(hd0D0VtxTrueptLSCotherPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCotherSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCotherSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCotherSB[i]->SetYTitle("Entries");
    flistLsCutsOther->Add(hd0D0ptLSCotherSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCotherSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCotherSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCotherSB[i]->SetYTitle("Entries");
    flistLsCutsOther->Add(hMCd0D0ptLSCotherSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCotherSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCotherSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCotherSB[i]->SetYTitle("Entries");
    flistLsCutsOther->Add(hd0D0VtxTrueptLSCotherSB[i]);
  }




  //################################################################################################
  //                                                                                               #
  //                         HISTOS FOR TIGHT CUTS                                                 #
  //                                                                                               #
  //################################################################################################

  //############ TIGHT CUTS SIGNAL HISTOGRAMS ###############
  //
  // ####### global properties histo ############

  TH2F *hCPtaVSd0d0TGHCsign=new TH2F("hCPtaVSd0d0TGHCsign","hCPtaVSd0d0_TightCuts_Signal",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZTGHCsign=new TH1F("hSecVtxZTGHCsign","hSecVtxZ_TightCuts_Signal",1000,-8.,8.);
  TH1F *hSecVtxXTGHCsign=new TH1F("hSecVtxXTGHCsign","hSecVtxX_TightCuts_Signal",1000,-3000.,3000.);
  TH1F *hSecVtxYTGHCsign=new TH1F("hSecVtxYTGHCsign","hSecVtxY_TightCuts_Signal",1000,-3000.,3000.);
  TH2F *hSecVtxXYTGHCsign=new TH2F("hSecVtxXYTGHCsign","hSecVtxXY_TightCuts_Signal",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiTGHCsign=new TH1F("hSecVtxPhiTGHCsign","hSecVtxPhi_TightCuts_Signal",180,-180.1,180.1);
  TH1F *hCPtaTGHCsign=new TH1F("hCPtaTGHCsign","hCPta_TightCuts_Signal",100,0.,1.);
  TH1F *hd0xd0TGHCsign=new TH1F("hd0xd0TGHCsign","hd0xd0_TightCuts_Signal",1000,-100000.,100000.);
  TH1F *hMassTrueTGHCsign=new TH1F("hMassTrueTGHCsign","D^{0} MC inv. Mass Tight Cuts Signal(All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCsign=new TH1F("hMassTGHCsign","D^{0} inv. Mass Tight Cuts Signal (All momenta)",600,1.600,2.200);
  hMassTGHCsign->Sumw2();
  TH1F *hMassTrueTGHCsignPM=new TH1F("hMassTrueTGHCsignPM","D^{0} MC inv. Mass Tight Cuts Signal, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCsignPM=new TH1F("hMassTGHCsignPM","D^{0} inv. Mass Tight Cuts Signal (All momenta), MassPeak",600,1.600,2.200);
  hMassTGHCsignPM->Sumw2();
  TH1F *hMassTrueTGHCsignSB=new TH1F("hMassTrueTGHCsignSB","D^{0} MC inv. Mass in Side Bands Tight Cuts Signal(All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCsignSB=new TH1F("hMassTGHCsignSB","D^{0} inv. Mass in Side Bands Tight Cuts Signal (All momenta)",600,1.600,2.200);
  hMassTGHCsignSB->Sumw2();

  flistTghCutsSignal->Add(hCPtaVSd0d0TGHCsign);
  flistTghCutsSignal->Add(hSecVtxZTGHCsign);
  flistTghCutsSignal->Add(hSecVtxYTGHCsign);
  flistTghCutsSignal->Add(hSecVtxXTGHCsign);
  flistTghCutsSignal->Add(hSecVtxXYTGHCsign);
  flistTghCutsSignal->Add(hSecVtxPhiTGHCsign);
  flistTghCutsSignal->Add(hCPtaTGHCsign);
  flistTghCutsSignal->Add(hd0xd0TGHCsign);
  flistTghCutsSignal->Add(hMassTrueTGHCsign);
  flistTghCutsSignal->Add(hMassTGHCsign);
  flistTghCutsSignal->Add(hMassTrueTGHCsignPM);
  flistTghCutsSignal->Add(hMassTGHCsignPM);
  flistTghCutsSignal->Add(hMassTrueTGHCsignSB);
  flistTghCutsSignal->Add(hMassTGHCsignSB);

  // ####### d0 D0 histos ############
  TH1F *hd0D0TGHCsignPM = new TH1F("hd0D0TGHCsignPM","D^{0} impact par. plot , Tight Cuts ,Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0TGHCsignPM->SetXTitle("Impact parameter [#mum]");
  hd0D0TGHCsignPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueTGHCsignPM = new TH1F("hd0D0VtxTrueTGHCsignPM","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueTGHCsignPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueTGHCsignPM->SetYTitle("Entries");

  TH1F *hMCd0D0TGHCsignPM = new TH1F("hMCd0D0TGHCsignPM","D^{0} impact par. plot, Tight Cuts, Signal,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0TGHCsignPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0TGHCsignPM->SetYTitle("Entries");

  TH1F *hd0D0TGHCsignSB = new TH1F("hd0D0TGHCsignSB","D^{0} impact par. plot , Tight Cuts ,Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0TGHCsignSB->SetXTitle("Impact parameter [#mum]");
  hd0D0TGHCsignSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueTGHCsignSB = new TH1F("hd0D0VtxTrueTGHCsignSB","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueTGHCsignSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueTGHCsignSB->SetYTitle("Entries");

  TH1F *hMCd0D0TGHCsignSB = new TH1F("hMCd0D0TGHCsignSB","D^{0} impact par. plot, Tight Cuts, Signal,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0TGHCsignSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0TGHCsignSB->SetYTitle("Entries");

  flistTghCutsSignal->Add(hd0D0TGHCsignPM);
  flistTghCutsSignal->Add(hd0D0VtxTrueTGHCsignPM);
  flistTghCutsSignal->Add(hMCd0D0TGHCsignPM);
  flistTghCutsSignal->Add(hd0D0TGHCsignSB);
  flistTghCutsSignal->Add(hd0D0VtxTrueTGHCsignSB);
  flistTghCutsSignal->Add(hMCd0D0TGHCsignSB);
  
  TH1F **hd0D0ptTGHCsignPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptTGHCsignPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptTGHCsignPM=new TH1F*[fnbins];
  TH1F **hd0D0ptTGHCsignSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptTGHCsignSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptTGHCsignSB=new TH1F*[fnbins];
  namehist="hd0D0ptTGHCsign_";
  titlehist="D^{0} impact par. plot, Tight Cuts, Signal, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCsignPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCsignPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCsignPM[i]->SetYTitle("Entries");
    flistTghCutsSignal->Add(hd0D0ptTGHCsignPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCsignPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCsignPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCsignPM[i]->SetYTitle("Entries");
    flistTghCutsSignal->Add(hMCd0D0ptTGHCsignPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCsignPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCsignPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCsignPM[i]->SetYTitle("Entries");
    flistTghCutsSignal->Add(hd0D0VtxTrueptTGHCsignPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCsignSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCsignSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCsignSB[i]->SetYTitle("Entries");
    flistTghCutsSignal->Add(hd0D0ptTGHCsignSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCsignSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCsignSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCsignSB[i]->SetYTitle("Entries");
    flistTghCutsSignal->Add(hMCd0D0ptTGHCsignSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCsignSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCsignSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCsignSB[i]->SetYTitle("Entries");
    flistTghCutsSignal->Add(hd0D0VtxTrueptTGHCsignSB[i]);
  }


  //############ TIGHT CUTS BACKGROUND HISTOGRAMS ###########
  //
  //   ######## global properties histos #######
  TH2F *hCPtaVSd0d0TGHCback=new TH2F("hCPtaVSd0d0TGHCback","hCPtaVSd0d0_TightCuts_Background",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZTGHCback=new TH1F("hSecVtxZTGHCback","hSecVtxZ_TightCuts_Background",1000,-8.,8.);
  TH1F *hSecVtxXTGHCback=new TH1F("hSecVtxXTGHCback","hSecVtxX_TightCuts_Background",1000,-3000.,3000.);
  TH1F *hSecVtxYTGHCback=new TH1F("hSecVtxYTGHCback","hSecVtxY_TightCuts_Background",1000,-3000.,3000.);
  TH2F *hSecVtxXYTGHCback=new TH2F("hSecVtxXYTGHCback","hSecVtxXY_TightCuts_Background",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiTGHCback=new TH1F("hSecVtxPhiTGHCback","hSecVtxPhi_TightCuts_Background",180,-180.1,180.1);
  TH1F *hCPtaTGHCback=new TH1F("hCPtaTGHCback","hCPta_TightCuts_Background",100,0.,1.);
  TH1F *hd0xd0TGHCback=new TH1F("hd0xd0TGHCback","hd0xd0_TightCuts_Background",1000,-100000.,100000.);
  TH1F *hMassTrueTGHCback=new TH1F("hMassTrueTGHCback","D^{0} MC inv. Mass Tight Cuts Background(All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCback=new TH1F("hMassTGHCback","D^{0} inv. Mass Tight Cuts Background (All momenta)",600,1.600,2.200);
  hMassTGHCback->Sumw2();
  TH1F *hMassTrueTGHCbackPM=new TH1F("hMassTrueTGHCbackPM","D^{0} MC inv. Mass Tight Cuts Background, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCbackPM=new TH1F("hMassTGHCbackPM","D^{0} inv. Mass Tight Cuts Background (All momenta), MassPeak",600,1.600,2.200);
  hMassTGHCbackPM->Sumw2();
  TH1F *hMassTrueTGHCbackSB=new TH1F("hMassTrueTGHCbackSB","D^{0} MC inv. Mass in Side Bands Tight Cuts Backgrround(All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCbackSB=new TH1F("hMassTGHCbackSB","D^{0} inv. Mass in Side Bands Tight Cuts Background (All momenta)",600,1.600,2.200);
  hMassTGHCbackSB->Sumw2();

  flistTghCutsBack->Add(hCPtaVSd0d0TGHCback);
  flistTghCutsBack->Add(hSecVtxZTGHCback);
  flistTghCutsBack->Add(hSecVtxYTGHCback);
  flistTghCutsBack->Add(hSecVtxXTGHCback);
  flistTghCutsBack->Add(hSecVtxXYTGHCback);
  flistTghCutsBack->Add(hSecVtxPhiTGHCback);
  flistTghCutsBack->Add(hCPtaTGHCback);
  flistTghCutsBack->Add(hd0xd0TGHCback);
  flistTghCutsBack->Add(hMassTrueTGHCback);
  flistTghCutsBack->Add(hMassTGHCback);
  flistTghCutsBack->Add(hMassTrueTGHCbackPM);
  flistTghCutsBack->Add(hMassTGHCbackPM);
  flistTghCutsBack->Add(hMassTrueTGHCbackSB);
  flistTghCutsBack->Add(hMassTGHCbackSB);


  // ####### d0 D0 histos ############
  
 TH1F *hd0D0TGHCbackPM = new TH1F("hd0D0TGHCbackPM","D^{0} impact par. plot , Tight Cuts ,Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0TGHCbackPM->SetXTitle("Impact parameter [#mum]");
  hd0D0TGHCbackPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueTGHCbackPM = new TH1F("hd0D0VtxTrueTGHCbackPM","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueTGHCbackPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueTGHCbackPM->SetYTitle("Entries");

  TH1F *hMCd0D0TGHCbackPM = new TH1F("hMCd0D0TGHCbackPM","D^{0} impact par. plot, Tight Cuts, Background,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0TGHCbackPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0TGHCbackPM->SetYTitle("Entries");

  TH1F *hd0D0TGHCbackSB = new TH1F("hd0D0TGHCbackSB","D^{0} impact par. plot , Tight Cuts ,Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0TGHCbackSB->SetXTitle("Impact parameter [#mum]");
  hd0D0TGHCbackSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueTGHCbackSB = new TH1F("hd0D0VtxTrueTGHCbackSB","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueTGHCbackSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueTGHCbackSB->SetYTitle("Entries");

  TH1F *hMCd0D0TGHCbackSB = new TH1F("hMCd0D0TGHCbackSB","D^{0} impact par. plot, Tight Cuts, Background,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0TGHCbackSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0TGHCbackSB->SetYTitle("Entries");

  flistTghCutsBack->Add(hd0D0TGHCbackPM);
  flistTghCutsBack->Add(hd0D0VtxTrueTGHCbackPM);
  flistTghCutsBack->Add(hMCd0D0TGHCbackPM);
  flistTghCutsBack->Add(hd0D0TGHCbackSB);
  flistTghCutsBack->Add(hd0D0VtxTrueTGHCbackSB);
  flistTghCutsBack->Add(hMCd0D0TGHCbackSB);
  
  TH1F **hd0D0ptTGHCbackPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptTGHCbackPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptTGHCbackPM=new TH1F*[fnbins];
  TH1F **hd0D0ptTGHCbackSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptTGHCbackSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptTGHCbackSB=new TH1F*[fnbins];
  namehist="hd0D0ptTGHCback_";
  titlehist="D^{0} impact par. plot, Tight Cuts, Background, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCbackPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCbackPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCbackPM[i]->SetYTitle("Entries");
    flistTghCutsBack->Add(hd0D0ptTGHCbackPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCbackPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCbackPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCbackPM[i]->SetYTitle("Entries");
    flistTghCutsBack->Add(hMCd0D0ptTGHCbackPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCbackPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCbackPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCbackPM[i]->SetYTitle("Entries");
    flistTghCutsBack->Add(hd0D0VtxTrueptTGHCbackPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCbackSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCbackSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCbackSB[i]->SetYTitle("Entries");
    flistTghCutsBack->Add(hd0D0ptTGHCbackSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCbackSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCbackSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCbackSB[i]->SetYTitle("Entries");
    flistTghCutsBack->Add(hMCd0D0ptTGHCbackSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCbackSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCbackSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCbackSB[i]->SetYTitle("Entries");
    flistTghCutsBack->Add(hd0D0VtxTrueptTGHCbackSB[i]);
  }



 //############ TIGHT CUTS FROMB HISTOGRAMS ###########
  //
  //#######  global properties histos

  TH2F *hCPtaVSd0d0TGHCfromB=new TH2F("hCPtaVSd0d0TGHCfromB","hCPtaVSd0d0_TightCuts_FromB",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZTGHCfromB=new TH1F("hSecVtxZTGHCfromB","hSecVtxZ_TightCuts_FromB",1000,-8.,8.);
  TH1F *hSecVtxXTGHCfromB=new TH1F("hSecVtxXTGHCfromB","hSecVtxX_TightCuts_FromB",1000,-3000.,3000.);
  TH1F *hSecVtxYTGHCfromB=new TH1F("hSecVtxYTGHCfromB","hSecVtxY_TightCuts_FromB",1000,-3000.,3000.);
  TH2F *hSecVtxXYTGHCfromB=new TH2F("hSecVtxXYTGHCfromB","hSecVtxXY_TightCuts_FromB",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiTGHCfromB=new TH1F("hSecVtxPhiTGHCfromB","hSecVtxPhi_TightCuts_FromB",180,-180.1,180.1);
  TH1F *hCPtaTGHCfromB=new TH1F("hCPtaTGHCfromB","hCPta_TightCuts_FromB",100,0.,1.);
  TH1F *hd0xd0TGHCfromB=new TH1F("hd0xd0TGHCfromB","hd0xd0_TightCuts_FromB",1000,-100000.,100000.);
  TH1F *hMassTrueTGHCfromB=new TH1F("hMassTrueTGHCfromB","D^{0} MC inv. Mass Tight Cuts FromB(All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCfromB=new TH1F("hMassTGHCfromB","D^{0} inv. Mass Tight Cuts FromB (All momenta)",600,1.600,2.200);
  hMassTGHCfromB->Sumw2();
  TH1F *hMassTrueTGHCfromBPM=new TH1F("hMassTrueTGHCfromBPM","D^{0} MC inv. Mass Tight Cuts FromB, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCfromBPM=new TH1F("hMassTGHCfromBPM","D^{0} inv. Mass Tight Cuts FromB (All momenta), MassPeak",600,1.600,2.200);
  hMassTGHCfromBPM->Sumw2();
  TH1F *hMassTrueTGHCfromBSB=new TH1F("hMassTrueTGHCfromBSB","D^{0} MC inv. Mass in Side Bands Tight Cuts FromB(All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCfromBSB=new TH1F("hMassTGHCfromBSB","D^{0} inv. Mass in Side Bands Tight Cuts FromB (All momenta)",600,1.600,2.200);
  hMassTGHCfromBSB->Sumw2();

  flistTghCutsFromB->Add(hCPtaVSd0d0TGHCfromB);
  flistTghCutsFromB->Add(hSecVtxZTGHCfromB);
  flistTghCutsFromB->Add(hSecVtxYTGHCfromB);
  flistTghCutsFromB->Add(hSecVtxXTGHCfromB);
  flistTghCutsFromB->Add(hSecVtxXYTGHCfromB);
  flistTghCutsFromB->Add(hSecVtxPhiTGHCfromB);
  flistTghCutsFromB->Add(hCPtaTGHCfromB);
  flistTghCutsFromB->Add(hd0xd0TGHCfromB);
  flistTghCutsFromB->Add(hMassTrueTGHCfromB);
  flistTghCutsFromB->Add(hMassTGHCfromB);
  flistTghCutsFromB->Add(hMassTrueTGHCfromBPM);
  flistTghCutsFromB->Add(hMassTGHCfromBPM);
  flistTghCutsFromB->Add(hMassTrueTGHCfromBSB);
  flistTghCutsFromB->Add(hMassTGHCfromBSB);

  // ######### d0 D0 histos ##############
  TH1F *hd0D0TGHCfromBPM = new TH1F("hd0D0TGHCfromBPM","D^{0} impact par. plot , Tight Cuts ,FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0TGHCfromBPM->SetXTitle("Impact parameter [#mum]");
  hd0D0TGHCfromBPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueTGHCfromBPM = new TH1F("hd0D0VtxTrueTGHCfromBPM","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueTGHCfromBPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueTGHCfromBPM->SetYTitle("Entries");

  TH1F *hMCd0D0TGHCfromBPM = new TH1F("hMCd0D0TGHCfromBPM","D^{0} impact par. plot, Tight Cuts, FromB,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0TGHCfromBPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0TGHCfromBPM->SetYTitle("Entries");

  TH1F *hd0D0TGHCfromBSB = new TH1F("hd0D0TGHCfromBSB","D^{0} impact par. plot , Tight Cuts ,FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0TGHCfromBSB->SetXTitle("Impact parameter [#mum]");
  hd0D0TGHCfromBSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueTGHCfromBSB = new TH1F("hd0D0VtxTrueTGHCfromBSB","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueTGHCfromBSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueTGHCfromBSB->SetYTitle("Entries");

  TH1F *hMCd0D0TGHCfromBSB = new TH1F("hMCd0D0TGHCfromBSB","D^{0} impact par. plot, Tight Cuts, FromB,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0TGHCfromBSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0TGHCfromBSB->SetYTitle("Entries");

  flistTghCutsFromB->Add(hd0D0TGHCfromBPM);
  flistTghCutsFromB->Add(hd0D0VtxTrueTGHCfromBPM);
  flistTghCutsFromB->Add(hMCd0D0TGHCfromBPM);
  flistTghCutsFromB->Add(hd0D0TGHCfromBSB);
  flistTghCutsFromB->Add(hd0D0VtxTrueTGHCfromBSB);
  flistTghCutsFromB->Add(hMCd0D0TGHCfromBSB);
  
  TH1F **hd0D0ptTGHCfromBPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptTGHCfromBPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptTGHCfromBPM=new TH1F*[fnbins];
  TH1F **hd0D0ptTGHCfromBSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptTGHCfromBSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptTGHCfromBSB=new TH1F*[fnbins];
  namehist="hd0D0ptTGHCfromB_";
  titlehist="D^{0} impact par. plot, Tight Cuts, FromB, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCfromBPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCfromBPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCfromBPM[i]->SetYTitle("Entries");
    flistTghCutsFromB->Add(hd0D0ptTGHCfromBPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCfromBPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCfromBPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCfromBPM[i]->SetYTitle("Entries");
    flistTghCutsFromB->Add(hMCd0D0ptTGHCfromBPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCfromBPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCfromBPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCfromBPM[i]->SetYTitle("Entries");
    flistTghCutsFromB->Add(hd0D0VtxTrueptTGHCfromBPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCfromBSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCfromBSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCfromBSB[i]->SetYTitle("Entries");
    flistTghCutsFromB->Add(hd0D0ptTGHCfromBSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCfromBSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCfromBSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCfromBSB[i]->SetYTitle("Entries");
    flistTghCutsFromB->Add(hMCd0D0ptTGHCfromBSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCfromBSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCfromBSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCfromBSB[i]->SetYTitle("Entries");
    flistTghCutsFromB->Add(hd0D0VtxTrueptTGHCfromBSB[i]);
  }



 //############ TIGHT CUTS FROM DSTAR HISTOGRAMS ###########
 //
  //############## global properties histos
  TH2F *hCPtaVSd0d0TGHCfromDstar=new TH2F("hCPtaVSd0d0TGHCfromDstar","hCPtaVSd0d0_TightCuts_FromDStar",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZTGHCfromDstar=new TH1F("hSecVtxZTGHCfromDstar","hSecVtxZ_TightCuts_FromDStar",1000,-8.,8.);
  TH1F *hSecVtxXTGHCfromDstar=new TH1F("hSecVtxXTGHCfromDstar","hSecVtxX_TightCuts_FromDStar",1000,-3000.,3000.);
  TH1F *hSecVtxYTGHCfromDstar=new TH1F("hSecVtxYTGHCfromDstar","hSecVtxY_TightCuts_FromDStar",1000,-3000.,3000.);
  TH2F *hSecVtxXYTGHCfromDstar=new TH2F("hSecVtxXYTGHCfromDstar","hSecVtxXY_TightCuts_FromDStar",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiTGHCfromDstar=new TH1F("hSecVtxPhiTGHCfromDstar","hSecVtxPhi_TightCuts_FromDStar",180,-180.1,180.1);
  TH1F *hCPtaTGHCfromDstar=new TH1F("hCPtaTGHCfromDstar","hCPta_TightCuts_FromDStar",100,0.,1.);
  TH1F *hd0xd0TGHCfromDstar=new TH1F("hd0xd0TGHCfromDstar","hd0xd0_TightCuts_FromDStar",1000,-100000.,100000.);
  TH1F *hMassTrueTGHCfromDstar=new TH1F("hMassTrueTGHCfromDstar","D^{0} MC inv. Mass Tight Cuts FromDStar(All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCfromDstar=new TH1F("hMassTGHCfromDstar","D^{0} inv. Mass Tight Cuts FromDStar (All momenta)",600,1.600,2.200);
  hMassTGHCfromDstar->Sumw2();
  TH1F *hMassTrueTGHCfromDstarPM=new TH1F("hMassTrueTGHCfromDstarPM","D^{0} MC inv. Mass Tight Cuts FromDStar, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCfromDstarPM=new TH1F("hMassTGHCfromDstarPM","D^{0} inv. Mass Tight Cuts FromDStar (All momenta), MassPeak",600,1.600,2.200);
  hMassTGHCfromDstarPM->Sumw2();
  TH1F *hMassTrueTGHCfromDstarSB=new TH1F("hMassTrueTGHCfromDstarSB","D^{0} MC inv. Mass in Side Bands Tight Cuts FromDStar(All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCfromDstarSB=new TH1F("hMassTGHCfromDstarSB","D^{0} inv. Mass in Side Bands Tight Cuts FromDStar (All momenta)",600,1.600,2.200);
  hMassTGHCfromDstarSB->Sumw2();

  flistTghCutsFromDstar->Add(hCPtaVSd0d0TGHCfromDstar);
  flistTghCutsFromDstar->Add(hSecVtxZTGHCfromDstar);
  flistTghCutsFromDstar->Add(hSecVtxYTGHCfromDstar);
  flistTghCutsFromDstar->Add(hSecVtxXTGHCfromDstar);
  flistTghCutsFromDstar->Add(hSecVtxXYTGHCfromDstar);
  flistTghCutsFromDstar->Add(hSecVtxPhiTGHCfromDstar);
  flistTghCutsFromDstar->Add(hCPtaTGHCfromDstar);
  flistTghCutsFromDstar->Add(hd0xd0TGHCfromDstar);
  flistTghCutsFromDstar->Add(hMassTrueTGHCfromDstar);
  flistTghCutsFromDstar->Add(hMassTGHCfromDstar);
  flistTghCutsFromDstar->Add(hMassTrueTGHCfromDstarPM);
  flistTghCutsFromDstar->Add(hMassTGHCfromDstarPM);
  flistTghCutsFromDstar->Add(hMassTrueTGHCfromDstarSB);
  flistTghCutsFromDstar->Add(hMassTGHCfromDstarSB);

  //########## d0 D0 histos #############  
  TH1F *hd0D0TGHCfromDstPM = new TH1F("hd0D0TGHCfromDstarPM","D^{0} impact par. plot , Tight Cuts ,FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0TGHCfromDstPM->SetXTitle("Impact parameter [#mum]");
  hd0D0TGHCfromDstPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueTGHCfromDstPM = new TH1F("hd0D0VtxTrueTGHCfromDstarPM","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueTGHCfromDstPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueTGHCfromDstPM->SetYTitle("Entries");

  TH1F *hMCd0D0TGHCfromDstPM = new TH1F("hMCd0D0TGHCfromDstarPM","D^{0} impact par. plot, Tight Cuts, FromDStar,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0TGHCfromDstPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0TGHCfromDstPM->SetYTitle("Entries");

  TH1F *hd0D0TGHCfromDstSB = new TH1F("hd0D0TGHCfromDstarSB","D^{0} impact par. plot , Tight Cuts ,FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0TGHCfromDstSB->SetXTitle("Impact parameter [#mum]");
  hd0D0TGHCfromDstSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueTGHCfromDstSB = new TH1F("hd0D0VtxTrueTGHCfromDstarSB","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueTGHCfromDstSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueTGHCfromDstSB->SetYTitle("Entries");

  TH1F *hMCd0D0TGHCfromDstSB = new TH1F("hMCd0D0TGHCfromDstarSB","D^{0} impact par. plot, Tight Cuts, FromDStar,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0TGHCfromDstSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0TGHCfromDstSB->SetYTitle("Entries");

  flistTghCutsFromDstar->Add(hd0D0TGHCfromDstPM);
  flistTghCutsFromDstar->Add(hd0D0VtxTrueTGHCfromDstPM);
  flistTghCutsFromDstar->Add(hMCd0D0TGHCfromDstPM);
  flistTghCutsFromDstar->Add(hd0D0TGHCfromDstSB);
  flistTghCutsFromDstar->Add(hd0D0VtxTrueTGHCfromDstSB);
  flistTghCutsFromDstar->Add(hMCd0D0TGHCfromDstSB);
  
  TH1F **hd0D0ptTGHCfromDstPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptTGHCfromDstPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptTGHCfromDstPM=new TH1F*[fnbins];
  TH1F **hd0D0ptTGHCfromDstSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptTGHCfromDstSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptTGHCfromDstSB=new TH1F*[fnbins];
  namehist="hd0D0ptTGHCfromDstar_";
  titlehist="D^{0} impact par. plot, Tight Cuts, FromDStar, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCfromDstPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCfromDstPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCfromDstPM[i]->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hd0D0ptTGHCfromDstPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCfromDstPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCfromDstPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCfromDstPM[i]->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hMCd0D0ptTGHCfromDstPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCfromDstPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCfromDstPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCfromDstPM[i]->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hd0D0VtxTrueptTGHCfromDstPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCfromDstSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCfromDstSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCfromDstSB[i]->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hd0D0ptTGHCfromDstSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCfromDstSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCfromDstSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCfromDstSB[i]->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hMCd0D0ptTGHCfromDstSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCfromDstSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCfromDstSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCfromDstSB[i]->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hd0D0VtxTrueptTGHCfromDstSB[i]);
  }


  //############ TIGHT CUTS OTHER HISTOGRAMS ###########
  //
  //########### global properties histos ###########

  TH2F *hCPtaVSd0d0TGHCother=new TH2F("hCPtaVSd0d0TGHCother","hCPtaVSd0d0_TightCuts_other",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZTGHCother=new TH1F("hSecVtxZTGHCother","hSecVtxZ_TightCuts_other",1000,-8.,8.);
  TH1F *hSecVtxXTGHCother=new TH1F("hSecVtxXTGHCother","hSecVtxX_TightCuts_other",1000,-3000.,3000.);
  TH1F *hSecVtxYTGHCother=new TH1F("hSecVtxYTGHCother","hSecVtxY_TightCuts_other",1000,-3000.,3000.);
  TH2F *hSecVtxXYTGHCother=new TH2F("hSecVtxXYTGHCother","hSecVtxXY_TightCuts_other",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiTGHCother=new TH1F("hSecVtxPhiTGHCother","hSecVtxPhi_TightCuts_other",180,-180.1,180.1);
  TH1F *hCPtaTGHCother=new TH1F("hCPtaTGHCother","hCPta_TightCuts_other",100,0.,1.);
  TH1F *hd0xd0TGHCother=new TH1F("hd0xd0TGHCother","hd0xd0_TightCuts_other",1000,-100000.,100000.);
  TH1F *hMassTrueTGHCother=new TH1F("hMassTrueTGHCother","D^{0} MC inv. Mass Tight Cuts other(All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCother=new TH1F("hMassTGHCother","D^{0} inv. Mass Tight Cuts other (All momenta)",600,1.600,2.200);
  hMassTGHCother->Sumw2();
  TH1F *hMassTrueTGHCotherPM=new TH1F("hMassTrueTGHCotherPM","D^{0} MC inv. Mass Tight Cuts other, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCotherPM=new TH1F("hMassTGHCotherPM","D^{0} inv. Mass Tight Cuts other (All momenta), MassPeak",600,1.600,2.200);
  hMassTGHCotherPM->Sumw2();
  TH1F *hMassTrueTGHCotherSB=new TH1F("hMassTrueTGHCotherSB","D^{0} MC inv. Mass in Side Bands Tight Cuts other(All momenta)",600,1.600,2.200);
  TH1F *hMassTGHCotherSB=new TH1F("hMassTGHCotherSB","D^{0} inv. Mass in Side Bands Tight Cuts other (All momenta)",600,1.600,2.200);
  hMassTGHCotherSB->Sumw2();

  flistTghCutsOther->Add(hCPtaVSd0d0TGHCother);
  flistTghCutsOther->Add(hSecVtxZTGHCother);
  flistTghCutsOther->Add(hSecVtxYTGHCother);
  flistTghCutsOther->Add(hSecVtxXTGHCother);
  flistTghCutsOther->Add(hSecVtxXYTGHCother);
  flistTghCutsOther->Add(hSecVtxPhiTGHCother);
  flistTghCutsOther->Add(hCPtaTGHCother);
  flistTghCutsOther->Add(hd0xd0TGHCother);
  flistTghCutsOther->Add(hMassTrueTGHCother);
  flistTghCutsOther->Add(hMassTGHCother);
  flistTghCutsOther->Add(hMassTrueTGHCotherPM);
  flistTghCutsOther->Add(hMassTGHCotherPM);
  flistTghCutsOther->Add(hMassTrueTGHCotherSB);
  flistTghCutsOther->Add(hMassTGHCotherSB);

  //############# d0 D0 histos ###############à
  TH1F *hd0D0TGHCotherPM = new TH1F("hd0D0TGHCotherPM","D^{0} impact par. plot , Tight Cuts ,Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0TGHCotherPM->SetXTitle("Impact parameter [#mum]");
  hd0D0TGHCotherPM->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueTGHCotherPM = new TH1F("hd0D0VtxTrueTGHCotherPM","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueTGHCotherPM->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueTGHCotherPM->SetYTitle("Entries");

  TH1F *hMCd0D0TGHCotherPM = new TH1F("hMCd0D0TGHCotherPM","D^{0} impact par. plot, Tight Cuts, Other,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0TGHCotherPM->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0TGHCotherPM->SetYTitle("Entries");

  TH1F *hd0D0TGHCotherSB = new TH1F("hd0D0TGHCotherSB","D^{0} impact par. plot , Tight Cuts ,Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0TGHCotherSB->SetXTitle("Impact parameter [#mum]");
  hd0D0TGHCotherSB->SetYTitle("Entries");

  TH1F *hd0D0VtxTrueTGHCotherSB = new TH1F("hd0D0VtxTrueTGHCotherSB","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrueTGHCotherSB->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrueTGHCotherSB->SetYTitle("Entries");

  TH1F *hMCd0D0TGHCotherSB = new TH1F("hMCd0D0TGHCotherSB","D^{0} impact par. plot, Tight Cuts, Other,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0TGHCotherSB->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0TGHCotherSB->SetYTitle("Entries");

  flistTghCutsOther->Add(hd0D0TGHCotherPM);
  flistTghCutsOther->Add(hd0D0VtxTrueTGHCotherPM);
  flistTghCutsOther->Add(hMCd0D0TGHCotherPM);
  flistTghCutsOther->Add(hd0D0TGHCotherSB);
  flistTghCutsOther->Add(hd0D0VtxTrueTGHCotherSB);
  flistTghCutsOther->Add(hMCd0D0TGHCotherSB);
  
  TH1F **hd0D0ptTGHCotherPM=new TH1F*[fnbins];
  TH1F **hMCd0D0ptTGHCotherPM=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptTGHCotherPM=new TH1F*[fnbins];
  TH1F **hd0D0ptTGHCotherSB=new TH1F*[fnbins];
  TH1F **hMCd0D0ptTGHCotherSB=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTrueptTGHCotherSB=new TH1F*[fnbins];
  namehist="hd0D0ptTGHCother_";
  titlehist="D^{0} impact par. plot, Tight Cuts, Other, ";
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCotherPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCotherPM[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCotherPM[i]->SetYTitle("Entries");
    flistTghCutsOther->Add(hd0D0ptTGHCotherPM[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCotherPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCotherPM[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCotherPM[i]->SetYTitle("Entries");
    flistTghCutsOther->Add(hMCd0D0ptTGHCotherPM[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCotherPM[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCotherPM[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCotherPM[i]->SetYTitle("Entries");
    flistTghCutsOther->Add(hd0D0VtxTrueptTGHCotherPM[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCotherSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCotherSB[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCotherSB[i]->SetYTitle("Entries");
    flistTghCutsOther->Add(hd0D0ptTGHCotherSB[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCotherSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCotherSB[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCotherSB[i]->SetYTitle("Entries");
    flistTghCutsOther->Add(hMCd0D0ptTGHCotherSB[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCotherSB[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCotherSB[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCotherSB[i]->SetYTitle("Entries");
    flistTghCutsOther->Add(hd0D0VtxTrueptTGHCotherSB[i]);
  }

  
}



//________________________________________________________________________
void AliAnalysisTaskSECharmFraction::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  TClonesArray *arrayD0toKpi=0;

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
      arrayD0toKpi=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
    }
  } else {
    arrayD0toKpi=(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
  }

  if(!arrayD0toKpi) {
    Printf("AliAnalysisTaskSECharmFraction::UserExec: D0toKpi branch not found!\n");
    return;
  }
  
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  
  
  // load MC particles
  TClonesArray *arrayMC = 
    (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(!arrayMC) {
    Printf("AliAnalysisTaskSECharmFraction::UserExec: MC particles branch not found!\n");
    return;
  }
  
  // load MC header
  AliAODMCHeader *aodmcHeader = 
    (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!aodmcHeader) {
    Printf("AliAnalysisTaskSECharmFraction::UserExec: MC header branch not found!\n");
    return;
  }
  
  // MC primary vertex
  Double_t vtxTrue[3];
  aodmcHeader->GetVertex(vtxTrue);
  
  if (!aod) {
    Printf("ERROR: aod not available");
    return;
  }
  //histogram filled with 1 for every AOD
  fNentries->Fill(1);
  PostData(1,fNentries);


  //Printf("There are %d tracks in this event", aod->GetNumberOfTracks());
  //  Int_t nTotHF=0,nTotDstar=0,nTot3Prong=0;
  Int_t nTotD0toKpi=0;
  Int_t okd0tight,okd0bartight,okd0loose,okd0barloose;
  Bool_t isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,isSideBand;
  Bool_t isinacceptance;
  Int_t signallevel=-1;
  Int_t ptbin; 
  //  const  Int_t nptbins=10;
  Double_t invMassD0,invMassD0bar,ptD0,massmumtrue;
  
 
  AliAODRecoDecayHF *aodDMC=0x0;// to be used to create a fake true sec vertex
  // make trkIDtoEntry register (temporary)
  Int_t trkIDtoEntry[100000];
  for(Int_t it=0;it<aod->GetNumberOfTracks();it++) {
    AliAODTrack *track = aod->GetTrack(it);
    if(track->GetID()<0) {
      //printf("Track ID <0, id= %d\n",track->GetID());
      return;
    }
    trkIDtoEntry[track->GetID()]=it;
  }
  
  
  // loop over D0->Kpi candidates
  Int_t nD0toKpi = arrayD0toKpi->GetEntriesFast();
  nTotD0toKpi += nD0toKpi;
  //	cout<<"Number of D0->Kpi: "<<nD0toKpi<<endl;
  
  for (Int_t iD0toKpi = 0; iD0toKpi < nD0toKpi; iD0toKpi++) {
    if(aodDMC!=0x0)delete aodDMC;
      
    isPeakD0=kFALSE;
    isPeakD0bar=kFALSE;
    isSideBandD0=kFALSE;
    isSideBandD0bar=kFALSE;
    isSideBand=kFALSE;
    isinacceptance=kFALSE;
    okd0tight=0;
    okd0bartight=0;
    okd0loose=0;
    okd0barloose=0;
  
    signallevel=-1;
    

    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->UncheckedAt(iD0toKpi);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
      
    
    //####### DATA SELECTION ####################################
    //
    // ######## CHECK FOR ACCEPTANCE ##########
    ptD0=d->Pt();
    isinacceptance = (TMath::Abs(d->EtaProng(0))<fAcceptanceCuts[0]&&TMath::Abs(d->EtaProng(1))<fAcceptanceCuts[0]); //eta acceptance
    
    //######## INVARIANT MASS SELECTION ###############
    CheckInvMassD0(d,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar);
    if((isSideBandD0||isSideBandD0bar)&&!(isPeakD0||isPeakD0bar))isSideBand=kTRUE;// TEMPORARY, NOT DONE IN THE METHOD CALLED ABOVE ONLY FOR FURTHER SIDE BAND STUDY
  
    // INVESTIGATE SIGNAL TYPE : ACCESS TO MC INFORMATION
    aodDMC=GetD0toKPiSignalType(d,arrayMC,signallevel,massmumtrue,vtxTrue);
    if(!isinacceptance)signallevel=9;
    fSignalType->Fill(signallevel);
  
    // END OF BACKGROUND TYPE SELECTION

    // NOW APPLY CUTS
    //NO CUTS CASE IS FOR FREE
    
    // CHECK TIGHTER CUTS 
    ptbin=SetStandardCuts(ptD0,flargeInvMassCut);
    d->SelectD0(fVHFtight->GetD0toKpiCuts(),okd0tight,okd0bartight);
    if(((isPeakD0&&okd0tight)||(isPeakD0bar&&okd0bartight))&&isinacceptance)fSignalTypeTghCuts->Fill(signallevel);
    
    // CHECK LOOSER CUTS 
    ptbin=SetStandardCuts(ptD0,flargeInvMassCut);
    d->SelectD0(fVHFloose->GetD0toKpiCuts(),okd0loose,okd0barloose);
    if(((isPeakD0&&okd0loose)||(isPeakD0bar&&okd0barloose))&&isinacceptance)fSignalTypeLsCuts->Fill(signallevel);
   
    
    //###################    FILL HISTOS      ########################
    //################################################################
    //
    //######## improvement: SPEED HERE CAN BE IMPROVED: CALCULATE ONCE AND FOR ALL 
    //            CANDIDATE VARIABLES   

    //NO CUTS Case: force okD0 and okD0bar = kTRUE
    if(signallevel==1)FillHistos(d,flistNoCutsSignal,ptbin,kTRUE,kTRUE,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==2)FillHistos(d,flistNoCutsFromDstar,ptbin,kTRUE,kTRUE,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==3||signallevel==4)FillHistos(d,flistNoCutsFromB,ptbin,kTRUE,kTRUE,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==-1||signallevel==7||signallevel==8||signallevel==10||signallevel==9)FillHistos(d,flistNoCutsBack,ptbin,kTRUE,kTRUE,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==5||signallevel==6)FillHistos(d,flistNoCutsOther,ptbin,kTRUE,kTRUE,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);

    //LOOSE CUTS Case
    if(signallevel==1)FillHistos(d,flistLsCutsSignal,ptbin,okd0loose,okd0barloose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==2)FillHistos(d,flistLsCutsFromDstar,ptbin,okd0loose,okd0barloose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==3||signallevel==4)FillHistos(d,flistLsCutsFromB,ptbin,okd0loose,okd0barloose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==-1||signallevel==7||signallevel==8||signallevel==10)FillHistos(d,flistLsCutsBack,ptbin,okd0loose,okd0barloose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==5||signallevel==6)FillHistos(d,flistLsCutsOther,ptbin,okd0loose,okd0barloose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);

    //TIGHT CUTS Case
    if(signallevel==1)FillHistos(d,flistTghCutsSignal,ptbin,okd0tight,okd0bartight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==2)FillHistos(d,flistTghCutsFromDstar,ptbin,okd0tight,okd0bartight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==3||signallevel==4)FillHistos(d,flistTghCutsFromB,ptbin,okd0tight,okd0bartight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==-1||signallevel==7||signallevel==8||signallevel==10)FillHistos(d,flistTghCutsBack,ptbin,okd0tight,okd0bartight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==5||signallevel==6)FillHistos(d,flistTghCutsOther,ptbin,okd0tight,okd0bartight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    
    
    if(aodDMC!=0x0){
      delete aodDMC;
      aodDMC=0x0;
    }
    
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
    
  }
  
  // ####################### POST OUTPUT TLIST DATA #########################
  // ####### histo for #AOD entries already posted
  
  PostData(2,fSignalType);
  PostData(3,fSignalTypeLsCuts);
  PostData(4,fSignalTypeTghCuts);
  PostData(5,flistNoCutsSignal);
  PostData(6,flistNoCutsBack);
  PostData(7,flistNoCutsFromB);
  PostData(8,flistNoCutsFromDstar);
  PostData(9,flistNoCutsOther);
  PostData(10,flistLsCutsSignal);
  PostData(11,flistLsCutsBack);
  PostData(12,flistLsCutsFromB);
  PostData(13,flistLsCutsFromDstar);
  PostData(14,flistLsCutsOther);
  PostData(15,flistTghCutsSignal);
  PostData(16,flistTghCutsBack);
  PostData(17,flistTghCutsFromB);
  PostData(18,flistTghCutsFromDstar);
  PostData(19,flistTghCutsOther);

  return;
}

//_________________________________________
Int_t AliAnalysisTaskSECharmFraction::SetStandardCuts(Double_t pt,Double_t invMassCut){
  //#############
  // TEMPORARY: to be change in :
  //             for(j<nptbins)
  //                       if pt < standardptbin[j+1]
  //                            SetCuts, bin=j
  //                            break 
  //                            
  // the way the cuts are set is for further development
  //   (to be interfaced with AliAnalsysTaskSETuneCuts)
  //
  // Cuts: 
  // 0 = inv. mass half width [GeV]
  // 1 = dca [cm]
  // 2 = cosThetaStar
  // 3 = pTK [GeV/c]
  // 4 = pTPi [GeV/c]
  // 5 = d0K [cm]   upper limit!
  // 6 = d0Pi [cm]  upper limit!
  // 7 = d0d0 [cm^2]
  // 8 = cosThetaPoint

  Int_t ptbin=-1;
  if (pt>0. && pt<=1.) {
    ptbin=0;
    fVHFtight->SetD0toKpiCuts(invMassCut,0.04,0.8,0.5,0.5,0.05,0.05,-0.0002,0.5);
    fVHFloose->SetD0toKpiCuts(invMassCut,0.04,0.8,0.5,0.5,0.05,0.05,-0.00025,0.7);
  }
  
  if(pt>1. && pt<=3.) {
    ptbin=1;  
    fVHFtight->SetD0toKpiCuts(invMassCut,0.03,0.8,0.6,0.6,0.05,0.05,-0.0002,0.6);
    fVHFloose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,1,1,-0.00025,0.8);
    //printf("I'm in the bin %d\n",ptbin);
  }
  
  if(pt>3. && pt<=5.){
    ptbin=2;  
    fVHFtight->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.0001,0.8);
    fVHFloose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.8);
    //printf("I'm in the bin %d\n",ptbin);
  }
  if(pt>5.){
    ptbin=3;
    fVHFtight->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00005,0.8);
    fVHFloose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.9);
  }//if(pt>5)
  return ptbin;
}

//__________________________________________________________
void AliAnalysisTaskSECharmFraction::CheckInvMassD0(AliAODRecoDecayHF2Prong *d,Double_t &invMassD0,Double_t &invMassD0bar,Bool_t &isPeakD0,Bool_t &isPeakD0bar,Bool_t &isSideBandD0,Bool_t &isSideBandD0bar){
  //Check wheter the candidate inv. mass is compatible with signal or sideband inv. mass selection
      
  d->InvMassD0(invMassD0,invMassD0bar);
  //CHECK if ISPEAK 
  if(TMath::Abs(invMassD0-fmD0PDG)<fsignalInvMassCut)isPeakD0=kTRUE;
  if(TMath::Abs(invMassD0bar-fmD0PDG)<fsignalInvMassCut)isPeakD0bar=kTRUE;
  //CHECK if ISSIDEBAND: no constraint is present between signal region definition and side band definition
  // ######## TO BE CHANGED the distinction between sidebandD0 and sidebandD0bar is meaningless 
  //               and it is present only for side band region study (see which inv mass has the D0(D0bar) 
  //               in case the D0bar(D0) is in the sideband) #######
  if(TMath::Abs(invMassD0-fmD0PDG)>fsidebandInvMassCut&&TMath::Abs(invMassD0-fmD0PDG)<fsidebandInvMassCut+fsidebandInvMassWindow){
    isSideBandD0=kTRUE;
  }
  if(TMath::Abs(invMassD0bar-fmD0PDG)>fsidebandInvMassCut&&TMath::Abs(invMassD0bar-fmD0PDG)<fsidebandInvMassCut+fsidebandInvMassWindow){
    isSideBandD0bar=kTRUE;
  }
  
}
	


//_______________________
AliAODRecoDecayHF* AliAnalysisTaskSECharmFraction::GetD0toKPiSignalType(const AliAODRecoDecayHF2Prong *d,TClonesArray *arrayMC,Int_t &signaltype,Double_t &massMumTrue,Double_t *primaryVtx){
  //THIS METHOD CHECK THE TYPE OF SIGNAL/BACKGROUND THE CANDIDATE IS. 
  //  IF (!AND ONLY IF) THE TWO DAUGHTERS COME FROM A COMMONE MOTHER A FAKE TRUE SECONDARY VERTEX IS CONSTRUCTED (aodDMC)  
  //
  // THE FOLLOWING SCHEME IS ADOPTED: signaltype is set to
                        //  1:signal (D0 prompt); 2: signal D0 from Dstar; 3: D0 fromB 4: D0 from Dstar fromB
                        // then background categories: -1: one or both daughters is a fake track
                        //                             5: both daughters come from a D meson != D0
                        //                             6: both daughters come from a D0->4prongs  
                        //                             7: both daughetrs are primaries
                        //                             8: generic background (can include one of the previous if desired)
                        //                             9: daughters out of acceptance
                        //                            10: pathologic cases (not clear)
                        //                            11: end of the method without output
                        //                            12: different result than MatchToMC method

  AliAODMCParticle *mum1=0x0;
  AliAODMCParticle *b1=0x0,*b2=0x0;
  AliAODMCParticle *grandmoth1=0x0;
  massMumTrue=-1;
  
  Int_t pdgmum,dglabels[2],matchtoMC;
  Int_t pdgdaughters[2]={211,321};
  // get daughter AOD tracks
  AliAODTrack *trk0 = (AliAODTrack*)d->GetDaughter(0);
  AliAODTrack *trk1 = (AliAODTrack*)d->GetDaughter(1);
  AliAODRecoDecayHF *aodDMC=0x0;
  if(trk0==0x0||trk1==0x0){
    AliDebug(2,"Delete tracks I AM \n");
  
    signaltype=-1;
    return aodDMC;
   
  }
  dglabels[0]=trk0->GetLabel();
  dglabels[1]=trk1->GetLabel();
  if(dglabels[0]==-1||dglabels[1]==-1){
    AliDebug(2,"HERE I AM \n");

    //fake tracks
    
    signaltype=-1;
    return aodDMC;

  }
  //      printf("Before entering the MC checks \n");
  
  b1=(AliAODMCParticle*)arrayMC->At(trk0->GetLabel());
  b2=(AliAODMCParticle*)arrayMC->At(trk1->GetLabel());
  
  if(b1->GetMother()==-1||b2->GetMother()==-1){
    //Tracks with no mother  ????? FAKE DECAY VERTEX
 
    
    signaltype=10;
    return aodDMC;
  }
  
  mum1=(AliAODMCParticle*)arrayMC->At(b1->GetMother());
  //  mum2=(AliAODMCParticle*)arrayMC->At(b2->GetMother());//FOR FURTHER USE
  
  if(b1->GetMother()!=b2->GetMother()){
    //Check the label of the mother is the same
    // NOT SAME MOTHER
   

    signaltype=8;
    return aodDMC;
  }
  massMumTrue=mum1->GetCalcMass();
  
  matchtoMC=d->MatchToMC(421,arrayMC,2,pdgdaughters);
  aodDMC=ConstructFakeTrueSecVtx(b1,b2,mum1,primaryVtx);
 
  if(aodDMC==0x0){
    signaltype=10;
    return aodDMC;
  }

  // if((mum1->GetPdgCode()!=mum2->GetPdgCode()))continue; //Check the mother is the same particle
  // printf("Particle codes: tr1: %d, tr2: %d, mum1: %d, mum 2: %d \n",b1->GetPdgCode(),b2->GetPdgCode(),mum1->GetPdgCode(),mum2->GetPdgCode());
  if(!((TMath::Abs(b1->GetPdgCode())==321&&TMath::Abs(b1->GetPdgCode())!=211)||(TMath::Abs(b1->GetPdgCode())==211&&TMath::Abs(b1->GetPdgCode()!=321)))){
    // Not a Kaon and a Pion
   
    signaltype=8;
    return aodDMC;
  }
  
  pdgmum=mum1->GetPdgCode();
  if(TMath::Abs(pdgmum)!=421){
    if(TMath::Abs(pdgmum)==411||TMath::Abs(pdgmum)==431||TMath::Abs(pdgmum)==443){
      // IT IS A SECONDARY VERTEX FROM CHARM BUT NOT A D0
      
      signaltype=5;
      return aodDMC;
    }
    else {
       signaltype=8;
       return aodDMC;
    }
  }

  if(mum1->GetDaughter(1)-mum1->GetDaughter(0)+1!=2){
    // from D0 but NOT A 2 PRONG DECAY
    signaltype=6;
    return aodDMC;
   
  }
  
  if(mum1->GetMother()==-1){
    // A particle coming from nothing
    signaltype=10;
    return aodDMC;
    
  }
  Bool_t isfromDstar=kFALSE;
  //  matchtoMC=d->MatchToMC(421,arrayMC,2,pdgdaughters);
  grandmoth1=(AliAODMCParticle*)arrayMC->At(mum1->GetMother());
  if(TMath::Abs(grandmoth1->GetPdgCode())==413||TMath::Abs(grandmoth1->GetPdgCode())==423)isfromDstar=kTRUE;// D0 COMING FROM A D0*
  
  /*
  //CHECK FOR CABIBBO SUPPRESSED DECAY
  Int_t isCabibSup=0,pdgKaon;
 
  pdgKaon=b1->GetPdgCode();
  if(TMath::Abs(pdgKaon)!=321)pdgKaon=b2->GetPdgCode();
  if(pdgmum>0&&pdgKaon>0)isCabibSup=1;
  if(pdgmum<0&&pdgKaon<0)isCabibSup=1;
  if(isCabibSup){
    signaltype=0;
    return aodDMC;
  }
  */
  //	  if(fcheckMCD0){//THIS CHECK IS NEEDED TO AVOID POSSIBLE FAILURE IN THE SECOND WHILE, FOR DEBUGGING
  while(TMath::Abs(grandmoth1->GetPdgCode())!=4&&TMath::Abs(grandmoth1->GetPdgCode())!=5){
    if(grandmoth1->GetMother()==-1){
      //### THE FOLLOWING IN CASE OF DEBUGGING ##########à
      /*printf("mother=-1, pdgcode: %d \n",grandmoth1->GetPdgCode());
	Int_t son=grandmoth1->GetDaughter(0);
	sonpart=(AliAODMCParticle*)arrayMC->At(son);
	while(TMath::Abs(sonpart->GetPdgCode())!=421){
	printf("mother=-1, pdgcode: %d \n",sonpart->GetPdgCode());
	son++;
	sonpart=(AliAODMCParticle*)arrayMC->At(son);
	}*/
   
      signaltype=10;
      return aodDMC;
    }
    grandmoth1=(AliAODMCParticle*)arrayMC->At(grandmoth1->GetMother());
  }
  
  if(TMath::Abs(grandmoth1->GetPdgCode())==4){
    if(matchtoMC!=-1){
      
      if(isfromDstar)signaltype=2;
      else signaltype=1;
      return aodDMC;
    }
    else {
      signaltype=12;
      return aodDMC;
      
    }
  }
  else if(TMath::Abs(grandmoth1->GetPdgCode())==5){
    if(matchtoMC!=-1){
      if(isfromDstar)signaltype=4;
      else signaltype=3;
      return aodDMC;
      
    }
    else {
     
      signaltype=12;
      return aodDMC;
    }
  }
  signaltype=11;// JUST FOR SAFETY: SHOULD NEVER REACH THIS POINT
  return aodDMC;
  //  return 11;
}

//___________________________________
AliAODRecoDecayHF* AliAnalysisTaskSECharmFraction::ConstructFakeTrueSecVtx(const AliAODMCParticle *b1, const AliAODMCParticle *b2, const AliAODMCParticle *mum,Double_t *primaryVtxTrue){
  // CONSTRUCT A FAKE TRUE SECONDARY VERTEX (aodDMC)  
  //!!!NOTE THAT ONLY ONE MOTHER IS CONSIDERED: THE METHOD REQUIRES THE DAUGHTERS COME FROM THE SAME MOTHER !!
  if(b1==0x0||b2==0x0)return 0x0;
  if(mum==0x0)return 0x0;
  Double_t pD[3],xD[3],pXtrTrue[2],pYtrTrue[2],pZtrTrue[2],xtr1[3],xtr2[3];
  Int_t charge[2]={0,0};
  if(b1->Charge()==-1)charge[0]=1;
  else {
    if(b2->Charge()==-1){
      //printf("Same charges for prongs \n");
      return 0x0;
    }
    charge[1]=1;
  }
  
  pXtrTrue[charge[0]]=b1->Px();
  pYtrTrue[charge[0]]=b1->Py();
  pZtrTrue[charge[0]]=b1->Pz();
  if(!b1->XvYvZv(xtr1)){
    return 0x0;
  }
  
  pXtrTrue[charge[1]]=b2->Px();
  pYtrTrue[charge[1]]=b2->Py();
  pZtrTrue[charge[1]]=b2->Pz();
  
  if(!mum->PxPyPz(pD)){
    //printf("!D from B:Get momentum failed \n");
    return 0x0;
  }
  if(!mum->XvYvZv(xD)){
    //printf("!D from B:Get position failed \n");
    return 0x0;
  }
  /* ############ THIS HAPPENS FROM TIME TO TIME: NUMERIC PROBLEM KNOWN #################
     if(pXtrTrue[0]+pXtrTrue[1]!=pD[0]){
     }*/
  
  
  if(!b2->XvYvZv(xtr2)){
    return 0x0;
  }
  Double_t d0dummy[2]={0.,0.};//TEMPORARY : dummy d0 for AliAODRecoDeay constructor
  AliAODRecoDecayHF* aodDMC=new AliAODRecoDecayHF(primaryVtxTrue,xD,2,0,pXtrTrue,pYtrTrue,pZtrTrue,d0dummy);
  
  /*   ######## THE FOLLOWINF FOR DEBUGGING ############
       Printf("testing the Fake vertex: SecVtxX: %f, Px: %f, Py: %f, Pz:%f \n ",aodDMC->GetSecVtxX(),aodDMC->Px(),aodDMC->Py(),aodDMC->Pz());
       Printf("pD: x=%f, y=%f,z=%f\n",pD[0],pD[1],pD[2]);
       Printf("Daughters :px1:%f, px2:%f \n",pXtrTrue[0],pXtrTrue[1]);
       Printf("Daughters :py1:%f, py2:%f \n",pYtrTrue[0],pYtrTrue[1]);
       Printf("Daughters :pz1:%f, pz2:%f \n",pZtrTrue[0],pZtrTrue[1]);
       Printf("Mother pdg: %d",mum->GetPdgCode());
       Printf("Impact Par Prod: %f\n",aodDMC->ImpParXY());
  */

  return aodDMC;
}

//________________________________________________________
Bool_t AliAnalysisTaskSECharmFraction::FillHistos(AliAODRecoDecayHF2Prong *d,TList *&list,Int_t ptbin,Int_t okD0,Int_t okD0bar,Double_t invMassD0,Double_t invMassD0bar,Bool_t isPeakD0,Bool_t isPeakD0bar,Bool_t isSideBand,Double_t massmumtrue,AliAODRecoDecayHF *aodDMC,Double_t *vtxTrue){//FILL THE HISTOGRAMS: TAKE THE HISTOS FROM THE list NAME

  
  if((!okD0)&&(!okD0bar))return kTRUE;
  
  // ######### Get Standard label for hist in tlist ###############
  TString namehist=list->GetName(),str;
  namehist.ReplaceAll("list","");

  //  ######### Global properties histos #################
  // ####### take care: only for candidates which pass the cuts !! not for side band ########
  if((isPeakD0&&okD0)||(isPeakD0bar&&okD0bar)){
    str="hCPtaVSd0d0";
    str.Append(namehist.Data());
    ((TH1F*)list->FindObject(str.Data()))->Fill(1e8*d->Prodd0d0(),d->CosPointingAngle());
    
    str="hSecVtxZ";
    str.Append(namehist.Data());
    ((TH1F*)list->FindObject(str.Data()))->Fill(d->GetSecVtxZ());
    
    str="hSecVtxX";
    str.Append(namehist.Data());
    ((TH1F*)list->FindObject(str.Data()))->Fill(d->GetSecVtxX()*10000.);
    
    str="hSecVtxY";
    str.Append(namehist.Data());
    ((TH1F*)list->FindObject(str.Data()))->Fill(d->GetSecVtxY()*10000.);
    
    str="hSecVtxXY";
    str.Append(namehist.Data());
    ((TH2F*)list->FindObject(str.Data()))->Fill(d->GetSecVtxX()*10000.,d->GetSecVtxY()*10000.);
    
    str="hSecVtxPhi";
    str.Append(namehist.Data());
    ((TH1F*)list->FindObject(str.Data()))->Fill(TMath::ATan2(d->GetSecVtxY()*10000.,d->GetSecVtxX()*10000.)*TMath::RadToDeg());
    
    str="hCPta";
    str.Append(namehist.Data());
    ((TH1F*)list->FindObject(str.Data()))->Fill(d->CosPointingAngle());
    
    str="hd0xd0";
    str.Append(namehist.Data());
    ((TH1F*)list->FindObject(str.Data()))->Fill(1e8*d->Prodd0d0());
  }
  
  
  //  ######### Invariant mass histos #################
  str="hMass";
  str.Append(namehist.Data());
  ((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0);
  ((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0bar);
   
      
  if(isPeakD0||isPeakD0bar){
    str="hMass";
    str.Append(namehist.Data());
    str.Append("PM");
    if(isPeakD0&&okD0)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0);
    if(isPeakD0bar&&okD0bar)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0bar);
  }
  if(isSideBand){
    str="hMass";
    str.Append(namehist.Data());
    str.Append("SB");
    if(okD0)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0);
    if(okD0bar)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0bar);
  }
  
  if(massmumtrue>0.){
    str="hMassTrue";
    str.Append(namehist.Data());
    ((TH1F*)list->FindObject(str.Data()))->Fill(massmumtrue);
    
    if(isPeakD0||isPeakD0bar){
      str="hMassTrue";
      str.Append(namehist.Data());
      str.Append("PM");
      ((TH1F*)list->FindObject(str.Data()))->Fill(massmumtrue);
    }
    if(isSideBand){
      str="hMassTrue";
      str.Append(namehist.Data());
      str.Append("SB");
      ((TH1F*)list->FindObject(str.Data()))->Fill(massmumtrue);
    }
  }
  
  // ################ D0 Impact Parameter Histos #####################
  if((isPeakD0&&okD0)||(isPeakD0bar&&okD0bar)){
    str="hd0D0";
    str.Append(namehist.Data());
    str.Append("PM");
    ((TH1F*)list->FindObject(str.Data()))->Fill(d->ImpParXY()*10000.);
    
    str="hd0D0pt";
    str.Append(namehist.Data());
    str.Append("_PkMss_pt");
    str+=ptbin;
    ((TH1F*)list->FindObject(str.Data()))->Fill(d->ImpParXY()*10000.);
     
    
    if(vtxTrue){
      str="hd0D0VtxTrue";
      str.Append(namehist.Data());
      str.Append("PM");
      ((TH1F*)list->FindObject(str.Data()))->Fill(d->AliAODRecoDecay::ImpParXY(vtxTrue)*10000.);
      
      str="hd0D0VtxTruept";
      str.Append(namehist.Data());
      str.Append("_PkMss_pt");
      str+=ptbin;
      ((TH1F*)list->FindObject(str.Data()))->Fill(d->AliAODRecoDecay::ImpParXY(vtxTrue)*10000.);
    }
    
    if(aodDMC!=0x0){
      aodDMC->Print("");
      aodDMC->ImpParXY();
      aodDMC->Print("");
      str="hMCd0D0";
      str.Append(namehist.Data());
      str.Append("PM");
      ((TH1F*)list->FindObject(str.Data()))->Fill(aodDMC->ImpParXY()*10000.);
      
      str="hMCd0D0pt";
      str.Append(namehist.Data());
      str.Append("_PkMss_pt");
      str+=ptbin;
      ((TH1F*)list->FindObject(str.Data()))->Fill(aodDMC->ImpParXY()*10000.);
    }
    
  }
  else if(isSideBand){
    str="hd0D0";
    str.Append(namehist.Data());
    str.Append("SB");
    ((TH1F*)list->FindObject(str.Data()))->Fill(d->ImpParXY()*10000.);
    
    str="hd0D0pt";
    str.Append(namehist.Data());
    str.Append("_SBMss_pt");
    str+=ptbin;
    ((TH1F*)list->FindObject(str.Data()))->Fill(d->ImpParXY()*10000.);
    
    
    if(vtxTrue){
      str="hd0D0VtxTrue";
      str.Append(namehist.Data());
      str.Append("SB");
      ((TH1F*)list->FindObject(str.Data()))->Fill(d->AliAODRecoDecay::ImpParXY(vtxTrue)*10000.);
      
      str="hd0D0VtxTruept";
      str.Append(namehist.Data());
      str.Append("_SBMss_pt");
      str+=ptbin;
      ((TH1F*)list->FindObject(str.Data()))->Fill(d->AliAODRecoDecay::ImpParXY(vtxTrue)*10000.);
      
    }
    
    if(aodDMC!=0x0){
      str="hMCd0D0";
      str.Append(namehist.Data());
      str.Append("SB");
      ((TH1F*)list->FindObject(str.Data()))->Fill(aodDMC->ImpParXY()*10000.);
      
      str="hMCd0D0pt";
      str.Append(namehist.Data());
      str.Append("_SBMss_pt");
      str+=ptbin;
      ((TH1F*)list->FindObject(str.Data()))->Fill(aodDMC->ImpParXY()*10000.);
    }
    
  }
     
   return kTRUE;
}


void AliAnalysisTaskSECharmFraction::SetNPtBins(Int_t nbins,const Double_t *ptbins){
  if((fptbins)!=0x0)delete fptbins;
  fnbins=nbins;fptbins=new Double_t[fnbins];
  memcpy(fptbins,ptbins,fnbins*sizeof(Double_t));
  return;
}

void AliAnalysisTaskSECharmFraction::SetStandardMassSelection(){
  //SET THE DEFAULT VALUES FOR INVARIANT MASS SELECTION
  SetSignalInvMassCut();
  SetLargeInvMassCut();
  SetSideBandInvMassCut();
  SetSideBandInvMassWindow();
  return;
  }


void AliAnalysisTaskSECharmFraction::Terminate(const Option_t*){
  //TERMINATE METHOD: NOTHING TO DO


}
