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

#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TROOT.h>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecay.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisTaskSECharmFraction.h"

ClassImp(AliAnalysisTaskSECharmFraction)
 
//________________________________________________________________________
  AliAnalysisTaskSECharmFraction::AliAnalysisTaskSECharmFraction() 
    : AliAnalysisTaskSE(),
      fVHFloose(0),
      fVHFtight(0),
      fmD0PDG(),
      fnbins(),
      fptbins(0),
      fsignalInvMassCut(),
      flargeInvMassCut(),
      fsidebandInvMassCut(),
      fsidebandInvMassWindow(),
      fUseMC(kTRUE),
      fNentries(0),
      fSignalType(0),
      fSignalTypeLsCuts(0),
      fSignalTypeTghCuts(0),
      flist_NoCuts_Signal(0),
      flist_NoCuts_Back(0),
      flist_NoCuts_FromB(0),
      flist_NoCuts_FromDstar(0),
      flist_NoCuts_Other(0),
      flist_LsCuts_Signal(0),
      flist_LsCuts_Back(0),
      flist_LsCuts_FromB(0),
      flist_LsCuts_FromDstar(0),
      flist_LsCuts_Other(0),
      flist_TghCuts_Signal(0),
      flist_TghCuts_Back(0),
      flist_TghCuts_FromB(0),
      flist_TghCuts_FromDstar(0),
      flist_TghCuts_Other(0)
   
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
      fsignalInvMassCut(),
      flargeInvMassCut(),
      fsidebandInvMassCut(),
      fsidebandInvMassWindow(),
      fUseMC(kTRUE),
      fNentries(0),
      fSignalType(0),
      fSignalTypeLsCuts(0),
      fSignalTypeTghCuts(0),
      flist_NoCuts_Signal(0),
      flist_NoCuts_Back(0),
      flist_NoCuts_FromB(0),
      flist_NoCuts_FromDstar(0),
      flist_NoCuts_Other(0),
      flist_LsCuts_Signal(0),
      flist_LsCuts_Back(0),
      flist_LsCuts_FromB(0),
      flist_LsCuts_FromDstar(0),
      flist_LsCuts_Other(0),
      flist_TghCuts_Signal(0),
      flist_TghCuts_Back(0),
      flist_TghCuts_FromB(0),
      flist_TghCuts_FromDstar(0),
      flist_TghCuts_Other(0)
   
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
    fsignalInvMassCut(),
    flargeInvMassCut(),
    fsidebandInvMassCut(),
    fsidebandInvMassWindow(),
    fUseMC(kTRUE),
    fNentries(0),
    fSignalType(0),
    fSignalTypeLsCuts(0),
    fSignalTypeTghCuts(0),
    flist_NoCuts_Signal(0),
    flist_NoCuts_Back(0),
    flist_NoCuts_FromB(0),
    flist_NoCuts_FromDstar(0),
    flist_NoCuts_Other(0),
    flist_LsCuts_Signal(0),
    flist_LsCuts_Back(0),
    flist_LsCuts_FromB(0),
    flist_LsCuts_FromDstar(0),
    flist_LsCuts_Other(0),
    flist_TghCuts_Signal(0),
    flist_TghCuts_Back(0),
    flist_TghCuts_FromB(0),
    flist_TghCuts_FromDstar(0),
    flist_TghCuts_Other(0)
{
  // Constructor
  // ptbins must be of dimension nptbins +1
  
  SetNPtBins(nptbins,ptbins);
  SetStandardMassSelection();
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
  if(flist_NoCuts_Signal){
    delete flist_NoCuts_Signal;
    flist_NoCuts_Signal=0;
  }
  if(flist_NoCuts_Back){
    delete flist_NoCuts_Back;
    flist_NoCuts_Back=0;
  }
  if(flist_NoCuts_FromB){
    delete flist_NoCuts_FromB;
    flist_NoCuts_FromB=0;
  }
  if(flist_NoCuts_FromDstar){
    delete flist_NoCuts_FromDstar;
    flist_NoCuts_FromDstar=0;
  }
  if(flist_NoCuts_Other){
    delete flist_NoCuts_Other;
    flist_NoCuts_Other=0;
  }
  
 if(flist_LsCuts_Signal){
    delete flist_LsCuts_Signal;
    flist_LsCuts_Signal=0;
  }
  if(flist_LsCuts_Back){
    delete flist_LsCuts_Back;
    flist_LsCuts_Back=0;
  }
  if(flist_LsCuts_FromB){
    delete flist_LsCuts_FromB;
    flist_LsCuts_FromB=0;
  }
  if(flist_LsCuts_FromDstar){
    delete flist_LsCuts_FromDstar;
    flist_LsCuts_FromDstar=0;
  }
  if(flist_LsCuts_Other){
    delete flist_LsCuts_Other;
    flist_LsCuts_Other=0;
  }
  
 if(flist_TghCuts_Signal){
    delete flist_TghCuts_Signal;
    flist_TghCuts_Signal=0;
  }
  if(flist_TghCuts_Back){
    delete flist_TghCuts_Back;
    flist_TghCuts_Back=0;
  }
  if(flist_TghCuts_FromB){
    delete flist_TghCuts_FromB;
    flist_TghCuts_FromB=0;
  }
  if(flist_TghCuts_FromDstar){
    delete flist_TghCuts_FromDstar;
    flist_TghCuts_FromDstar=0;
  }
  if(flist_TghCuts_Other){
    delete flist_TghCuts_Other;
    flist_TghCuts_Other=0;
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


  TString namehist;
  TString titlehist;
  TString strnamept,strtitlept;
 
  fNentries=new TH1F("nentriesChFr", "Look at the number of entries! = number of AODs", 2,1.,2.);
  fSignalType=new TH1F("hsignaltype", "Histo for type of MC signal", 21,-1.,20.);
  fSignalTypeLsCuts=new TH1F("hsignaltypeLsCuts", "Histo for type of MC signal with loose cuts", 21,-1.,20.);
  fSignalTypeTghCuts=new TH1F("hsignaltypeTghCuts", "Histo for type of MC signal with tight cuts", 21,-1.,20.);

  //##########  DEFINE THE TLISTS ##################
  
  flist_NoCuts_Signal = new TList();
  flist_NoCuts_Signal->SetOwner();
  flist_NoCuts_Signal->SetName("list_nc_sign");

  flist_NoCuts_Back = new TList();
  flist_NoCuts_Back->SetOwner();
  flist_NoCuts_Back->SetName("list_nc_back");

  flist_NoCuts_FromB = new TList();
  flist_NoCuts_FromB->SetOwner();
  flist_NoCuts_FromB->SetName("list_nc_fromB");

  flist_NoCuts_FromDstar = new TList();
  flist_NoCuts_FromDstar->SetOwner();
  flist_NoCuts_FromDstar->SetName("list_nc_fromDstar");

  flist_NoCuts_Other = new TList();
  flist_NoCuts_Other->SetOwner();
  flist_NoCuts_Other->SetName("list_nc_other");


  flist_LsCuts_Signal = new TList();
  flist_LsCuts_Signal->SetOwner();
  flist_LsCuts_Signal->SetName("list_ls_sign");

  flist_LsCuts_Back = new TList();
  flist_LsCuts_Back->SetOwner();
  flist_LsCuts_Back->SetName("list_ls_back");

  flist_LsCuts_FromB = new TList();
  flist_LsCuts_FromB->SetOwner();
  flist_LsCuts_FromB->SetName("list_ls_fromB");

  flist_LsCuts_FromDstar = new TList();
  flist_LsCuts_FromDstar->SetOwner();
  flist_LsCuts_FromDstar->SetName("list_ls_fromDstar");

  flist_LsCuts_Other = new TList();
  flist_LsCuts_Other->SetOwner();
  flist_LsCuts_Other->SetName("list_ls_other");


  flist_TghCuts_Signal = new TList();
  flist_TghCuts_Signal->SetOwner();
  flist_TghCuts_Signal->SetName("list_tgh_sign");

  flist_TghCuts_Back = new TList();
  flist_TghCuts_Back->SetOwner();
  flist_TghCuts_Back->SetName("list_tgh_back");

  flist_TghCuts_FromB = new TList();
  flist_TghCuts_FromB->SetOwner();
  flist_TghCuts_FromB->SetName("list_tgh_fromB");

  flist_TghCuts_FromDstar = new TList();
  flist_TghCuts_FromDstar->SetOwner();
  flist_TghCuts_FromDstar->SetName("list_tgh_fromDstar");

  flist_TghCuts_Other = new TList();
  flist_TghCuts_Other->SetOwner();
  flist_TghCuts_Other->SetName("list_tgh_other");




  //################################################################################################
  //                                                                                               #
  //                         HISTOS FOR NO CUTS CASE                                               #
  //                                                                                               #
  //################################################################################################


  //############ NO CUTS SIGNAL HISTOGRAMS ###############
  //
  // ####### global properties histo ############

  TH2F *hCPtaVSd0d0_nc_sign=new TH2F("hCPtaVSd0d0_nc_sign","hCPtaVSd0d0_NoCuts_Signal",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_nc_sign=new TH1F("hSecVtxZ_nc_sign","hSecVtxZ_NoCuts_Signal",1000,-8.,8.);
  TH1F *hSecVtxX_nc_sign=new TH1F("hSecVtxX_nc_sign","hSecVtxX_NoCuts_Signal",1000,-3000.,3000.);
  TH1F *hSecVtxY_nc_sign=new TH1F("hSecVtxY_nc_sign","hSecVtxY_NoCuts_Signal",1000,-3000.,3000.);
  TH2F *hSecVtxXY_nc_sign=new TH2F("hSecVtxXY_nc_sign","hSecVtxXY_NoCuts_Signal",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_nc_sign=new TH1F("hSecVtxPhi_nc_sign","hSecVtxPhi_NoCuts_Signal",180,-180.1,180.1);
  TH1F *hCPta_nc_sign=new TH1F("hCPta_nc_sign","hCPta_NoCuts_Signal",100,0.,1.);
  TH1F *hd0xd0_nc_sign=new TH1F("hd0xd0_nc_sign","hd0xd0_NoCuts_Signal",1000,-100000.,100000.);
  TH1F *hMassTrue_nc_sign=new TH1F("hMassTrue_nc_sign","D^{0} MC inv. Mass No Cuts Signal(All momenta)",600,1.600,2.200);
  TH1F *hMass_nc_sign=new TH1F("hMass_nc_sign","D^{0} inv. Mass No Cuts Signal (All momenta)",600,1.600,2.200);
  hMass_nc_sign->Sumw2();
  TH1F *hMassTrue_nc_sign_pm=new TH1F("hMassTrue_nc_sign_pm","D^{0} MC inv. Mass No Cuts Signal, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_nc_sign_pm=new TH1F("hMass_nc_sign_pm","D^{0} inv. Mass No Cuts Signal (All momenta), MassPeak",600,1.600,2.200);
  hMass_nc_sign_pm->Sumw2();

  TH1F *hMassTrue_SB_nc_sign=new TH1F("hMassTrue_nc_sign_sb","D^{0} MC inv. Mass in Side Bands No Cuts Signal(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_nc_sign=new TH1F("hMass_nc_sign_sb","D^{0} inv. Mass in Side Bands No Cuts Signal (All momenta)",600,1.600,2.200);
  hMass_SB_nc_sign->Sumw2();

  flist_NoCuts_Signal->Add(hCPtaVSd0d0_nc_sign);
  flist_NoCuts_Signal->Add(hSecVtxZ_nc_sign);
  flist_NoCuts_Signal->Add(hSecVtxY_nc_sign);
  flist_NoCuts_Signal->Add(hSecVtxX_nc_sign);
  flist_NoCuts_Signal->Add(hSecVtxXY_nc_sign);
  flist_NoCuts_Signal->Add(hSecVtxPhi_nc_sign);
  flist_NoCuts_Signal->Add(hCPta_nc_sign);
  flist_NoCuts_Signal->Add(hd0xd0_nc_sign);
  flist_NoCuts_Signal->Add(hMassTrue_nc_sign);
  flist_NoCuts_Signal->Add(hMass_nc_sign);
  flist_NoCuts_Signal->Add(hMassTrue_nc_sign_pm);
  flist_NoCuts_Signal->Add(hMass_nc_sign_pm);
  flist_NoCuts_Signal->Add(hMassTrue_SB_nc_sign);
  flist_NoCuts_Signal->Add(hMass_SB_nc_sign);

  // ####### d0 D0 histos ############
  TH1F *hd0D0_nc_sign_pm = new TH1F("hd0D0_nc_sign_pm","D^{0} impact par. plot , No Cuts ,Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_nc_sign_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_nc_sign_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_nc_sign_pm = new TH1F("hd0D0VtxTrue_nc_sign_pm","D^{0} impact par. w.r.t. True Vtx, No Cuts, Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_nc_sign_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_nc_sign_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_nc_sign_pm = new TH1F("hMCd0D0_nc_sign_pm","D^{0} impact par. plot, No Cuts, Signal,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_nc_sign_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_nc_sign_pm->SetYTitle("Entries");

  TH1F *hd0D0_nc_sign_sb = new TH1F("hd0D0_nc_sign_sb","D^{0} impact par. plot , No Cuts ,Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_nc_sign_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_nc_sign_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_nc_sign_sb = new TH1F("hd0D0VtxTrue_nc_sign_sb","D^{0} impact par. w.r.t. True Vtx, No Cuts, Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_nc_sign_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_nc_sign_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_nc_sign_sb = new TH1F("hMCd0D0_nc_sign_sb","D^{0} impact par. plot, No Cuts, Signal,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_nc_sign_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_nc_sign_sb->SetYTitle("Entries");

  flist_NoCuts_Signal->Add(hd0D0_nc_sign_pm);
  flist_NoCuts_Signal->Add(hd0D0VtxTrue_nc_sign_pm);
  flist_NoCuts_Signal->Add(hMCd0D0_nc_sign_pm);
  flist_NoCuts_Signal->Add(hd0D0_nc_sign_sb);
  flist_NoCuts_Signal->Add(hd0D0VtxTrue_nc_sign_sb);
  flist_NoCuts_Signal->Add(hMCd0D0_nc_sign_sb);
  
  TH1F **hd0D0pt_nc_sign_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_nc_sign_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_nc_sign_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_nc_sign_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_nc_sign_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_nc_sign_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_nc_sign_";
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
    
    hd0D0pt_nc_sign_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_nc_sign_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_nc_sign_pm[i]->SetYTitle("Entries");
    flist_NoCuts_Signal->Add(hd0D0pt_nc_sign_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_nc_sign_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_nc_sign_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_nc_sign_pm[i]->SetYTitle("Entries");
    flist_NoCuts_Signal->Add(hMCd0D0pt_nc_sign_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_nc_sign_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_nc_sign_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_nc_sign_pm[i]->SetYTitle("Entries");
    flist_NoCuts_Signal->Add(hd0D0VtxTruept_nc_sign_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_nc_sign_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_nc_sign_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_nc_sign_sb[i]->SetYTitle("Entries");
    flist_NoCuts_Signal->Add(hd0D0pt_nc_sign_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_nc_sign_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_nc_sign_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_nc_sign_sb[i]->SetYTitle("Entries");
    flist_NoCuts_Signal->Add(hMCd0D0pt_nc_sign_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_nc_sign_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_nc_sign_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_nc_sign_sb[i]->SetYTitle("Entries");
    flist_NoCuts_Signal->Add(hd0D0VtxTruept_nc_sign_sb[i]);
  }


  //############ NO CUTS BACKGROUND HISTOGRAMS ###########
  //
  //   ######## global properties histos #######
  TH2F *hCPtaVSd0d0_nc_back=new TH2F("hCPtaVSd0d0_nc_back","hCPtaVSd0d0_NoCuts_Background",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_nc_back=new TH1F("hSecVtxZ_nc_back","hSecVtxZ_NoCuts_Background",1000,-8.,8.);
  TH1F *hSecVtxX_nc_back=new TH1F("hSecVtxX_nc_back","hSecVtxX_NoCuts_Background",1000,-3000.,3000.);
  TH1F *hSecVtxY_nc_back=new TH1F("hSecVtxY_nc_back","hSecVtxY_NoCuts_Background",1000,-3000.,3000.);
  TH2F *hSecVtxXY_nc_back=new TH2F("hSecVtxXY_nc_back","hSecVtxXY_NoCuts_Background",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_nc_back=new TH1F("hSecVtxPhi_nc_back","hSecVtxPhi_NoCuts_Background",180,-180.1,180.1);
  TH1F *hCPta_nc_back=new TH1F("hCPta_nc_back","hCPta_NoCuts_Background",100,0.,1.);
  TH1F *hd0xd0_nc_back=new TH1F("hd0xd0_nc_back","hd0xd0_NoCuts_Background",1000,-100000.,100000.);
  TH1F *hMassTrue_nc_back=new TH1F("hMassTrue_nc_back","D^{0} MC inv. Mass No Cuts Background(All momenta)",600,1.600,2.200);
  TH1F *hMass_nc_back=new TH1F("hMass_nc_back","D^{0} inv. Mass No Cuts Background (All momenta)",600,1.600,2.200);
  hMass_nc_back->Sumw2();
  TH1F *hMassTrue_nc_back_pm=new TH1F("hMassTrue_nc_back_pm","D^{0} MC inv. Mass No Cuts Background, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_nc_back_pm=new TH1F("hMass_nc_back_pm","D^{0} inv. Mass No Cuts Background (All momenta), MassPeak",600,1.600,2.200);
  hMass_nc_back_pm->Sumw2();
  TH1F *hMassTrue_SB_nc_back=new TH1F("hMassTrue_nc_back_sb","D^{0} MC inv. Mass in Side Bands No Cuts Background(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_nc_back=new TH1F("hMass_nc_back_sb","D^{0} inv. Mass in Side Bands No Cuts Background (All momenta)",600,1.600,2.200);
  hMass_SB_nc_back->Sumw2();

  flist_NoCuts_Back->Add(hCPtaVSd0d0_nc_back);
  flist_NoCuts_Back->Add(hSecVtxZ_nc_back);
  flist_NoCuts_Back->Add(hSecVtxY_nc_back);
  flist_NoCuts_Back->Add(hSecVtxX_nc_back);
  flist_NoCuts_Back->Add(hSecVtxXY_nc_back);
  flist_NoCuts_Back->Add(hSecVtxPhi_nc_back);
  flist_NoCuts_Back->Add(hCPta_nc_back);
  flist_NoCuts_Back->Add(hd0xd0_nc_back);
  flist_NoCuts_Back->Add(hMassTrue_nc_back);
  flist_NoCuts_Back->Add(hMass_nc_back);
  flist_NoCuts_Back->Add(hMassTrue_nc_back_pm);
  flist_NoCuts_Back->Add(hMass_nc_back_pm);
  flist_NoCuts_Back->Add(hMassTrue_SB_nc_back);
  flist_NoCuts_Back->Add(hMass_SB_nc_back);


  // ####### d0 D0 histos ############
  
 TH1F *hd0D0_nc_back_pm = new TH1F("hd0D0_nc_back_pm","D^{0} impact par. plot , No Cuts ,Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_nc_back_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_nc_back_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_nc_back_pm = new TH1F("hd0D0VtxTrue_nc_back_pm","D^{0} impact par. w.r.t. True Vtx, No Cuts, Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_nc_back_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_nc_back_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_nc_back_pm = new TH1F("hMCd0D0_nc_back_pm","D^{0} impact par. plot, No Cuts, Background,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_nc_back_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_nc_back_pm->SetYTitle("Entries");

  TH1F *hd0D0_nc_back_sb = new TH1F("hd0D0_nc_back_sb","D^{0} impact par. plot , No Cuts ,Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_nc_back_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_nc_back_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_nc_back_sb = new TH1F("hd0D0VtxTrue_nc_back_sb","D^{0} impact par. w.r.t. True Vtx, No Cuts, Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_nc_back_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_nc_back_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_nc_back_sb = new TH1F("hMCd0D0_nc_back_sb","D^{0} impact par. plot, No Cuts, Background,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_nc_back_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_nc_back_sb->SetYTitle("Entries");

  flist_NoCuts_Back->Add(hd0D0_nc_back_pm);
  flist_NoCuts_Back->Add(hd0D0VtxTrue_nc_back_pm);
  flist_NoCuts_Back->Add(hMCd0D0_nc_back_pm);
  flist_NoCuts_Back->Add(hd0D0_nc_back_sb);
  flist_NoCuts_Back->Add(hd0D0VtxTrue_nc_back_sb);
  flist_NoCuts_Back->Add(hMCd0D0_nc_back_sb);
  
  TH1F **hd0D0pt_nc_back_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_nc_back_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_nc_back_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_nc_back_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_nc_back_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_nc_back_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_nc_back_";
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
    
    hd0D0pt_nc_back_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_nc_back_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_nc_back_pm[i]->SetYTitle("Entries");
    flist_NoCuts_Back->Add(hd0D0pt_nc_back_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_nc_back_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_nc_back_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_nc_back_pm[i]->SetYTitle("Entries");
    flist_NoCuts_Back->Add(hMCd0D0pt_nc_back_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_nc_back_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_nc_back_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_nc_back_pm[i]->SetYTitle("Entries");
    flist_NoCuts_Back->Add(hd0D0VtxTruept_nc_back_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_nc_back_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_nc_back_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_nc_back_sb[i]->SetYTitle("Entries");
    flist_NoCuts_Back->Add(hd0D0pt_nc_back_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_nc_back_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_nc_back_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_nc_back_sb[i]->SetYTitle("Entries");
    flist_NoCuts_Back->Add(hMCd0D0pt_nc_back_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_nc_back_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_nc_back_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_nc_back_sb[i]->SetYTitle("Entries");
    flist_NoCuts_Back->Add(hd0D0VtxTruept_nc_back_sb[i]);
  }



 //############ NO CUTS FROMB HISTOGRAMS ###########
  //
  //#######  global properties histos

  TH2F *hCPtaVSd0d0_nc_fromB=new TH2F("hCPtaVSd0d0_nc_fromB","hCPtaVSd0d0_NoCuts_FromB",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_nc_fromB=new TH1F("hSecVtxZ_nc_fromB","hSecVtxZ_NoCuts_FromB",1000,-8.,8.);
  TH1F *hSecVtxX_nc_fromB=new TH1F("hSecVtxX_nc_fromB","hSecVtxX_NoCuts_FromB",1000,-3000.,3000.);
  TH1F *hSecVtxY_nc_fromB=new TH1F("hSecVtxY_nc_fromB","hSecVtxY_NoCuts_FromB",1000,-3000.,3000.);
  TH2F *hSecVtxXY_nc_fromB=new TH2F("hSecVtxXY_nc_fromB","hSecVtxXY_NoCuts_FromB",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_nc_fromB=new TH1F("hSecVtxPhi_nc_fromB","hSecVtxPhi_NoCuts_FromB",180,-180.1,180.1);
  TH1F *hCPta_nc_fromB=new TH1F("hCPta_nc_fromB","hCPta_NoCuts_FromB",100,0.,1.);
  TH1F *hd0xd0_nc_fromB=new TH1F("hd0xd0_nc_fromB","hd0xd0_NoCuts_FromB",1000,-100000.,100000.);
  TH1F *hMassTrue_nc_fromB=new TH1F("hMassTrue_nc_fromB","D^{0} MC inv. Mass No Cuts FromB(All momenta)",600,1.600,2.200);
  TH1F *hMass_nc_fromB=new TH1F("hMass_nc_fromB","D^{0} inv. Mass No Cuts FromB (All momenta)",600,1.600,2.200);
  hMass_nc_fromB->Sumw2();
  TH1F *hMassTrue_nc_fromB_pm=new TH1F("hMassTrue_nc_fromB_pm","D^{0} MC inv. Mass No Cuts FromB, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_nc_fromB_pm=new TH1F("hMass_nc_fromB_pm","D^{0} inv. Mass No Cuts FromB (All momenta), MassPeak",600,1.600,2.200);
  hMass_nc_fromB->Sumw2();
  TH1F *hMassTrue_SB_nc_fromB=new TH1F("hMassTrue_nc_fromB_sb","D^{0} MC inv. Mass in Side Bands No Cuts FromB(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_nc_fromB=new TH1F("hMass_nc_fromB_sb","D^{0} inv. Mass in Side Bands No Cuts FromB (All momenta)",600,1.600,2.200);
  hMass_SB_nc_fromB->Sumw2();

  flist_NoCuts_FromB->Add(hCPtaVSd0d0_nc_fromB);
  flist_NoCuts_FromB->Add(hSecVtxZ_nc_fromB);
  flist_NoCuts_FromB->Add(hSecVtxY_nc_fromB);
  flist_NoCuts_FromB->Add(hSecVtxX_nc_fromB);
  flist_NoCuts_FromB->Add(hSecVtxXY_nc_fromB);
  flist_NoCuts_FromB->Add(hSecVtxPhi_nc_fromB);
  flist_NoCuts_FromB->Add(hCPta_nc_fromB);
  flist_NoCuts_FromB->Add(hd0xd0_nc_fromB);
  flist_NoCuts_FromB->Add(hMassTrue_nc_fromB);
  flist_NoCuts_FromB->Add(hMass_nc_fromB);
  flist_NoCuts_FromB->Add(hMassTrue_nc_fromB_pm);
  flist_NoCuts_FromB->Add(hMass_nc_fromB_pm);
  flist_NoCuts_FromB->Add(hMassTrue_SB_nc_fromB);
  flist_NoCuts_FromB->Add(hMass_SB_nc_fromB);

  // ######### d0 D0 histos ##############
  TH1F *hd0D0_nc_fromB_pm = new TH1F("hd0D0_nc_fromB_pm","D^{0} impact par. plot , No Cuts ,FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_nc_fromB_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_nc_fromB_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_nc_fromB_pm = new TH1F("hd0D0VtxTrue_nc_fromB_pm","D^{0} impact par. w.r.t. True Vtx, No Cuts, FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_nc_fromB_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_nc_fromB_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_nc_fromB_pm = new TH1F("hMCd0D0_nc_fromB_pm","D^{0} impact par. plot, No Cuts, FromB,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_nc_fromB_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_nc_fromB_pm->SetYTitle("Entries");

  TH1F *hd0D0_nc_fromB_sb = new TH1F("hd0D0_nc_fromB_sb","D^{0} impact par. plot , No Cuts ,FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_nc_fromB_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_nc_fromB_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_nc_fromB_sb = new TH1F("hd0D0VtxTrue_nc_fromB_sb","D^{0} impact par. w.r.t. True Vtx, No Cuts, FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_nc_fromB_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_nc_fromB_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_nc_fromB_sb = new TH1F("hMCd0D0_nc_fromB_sb","D^{0} impact par. plot, No Cuts, FromB,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_nc_fromB_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_nc_fromB_sb->SetYTitle("Entries");

  flist_NoCuts_FromB->Add(hd0D0_nc_fromB_pm);
  flist_NoCuts_FromB->Add(hd0D0VtxTrue_nc_fromB_pm);
  flist_NoCuts_FromB->Add(hMCd0D0_nc_fromB_pm);
  flist_NoCuts_FromB->Add(hd0D0_nc_fromB_sb);
  flist_NoCuts_FromB->Add(hd0D0VtxTrue_nc_fromB_sb);
  flist_NoCuts_FromB->Add(hMCd0D0_nc_fromB_sb);
  
  TH1F **hd0D0pt_nc_fromB_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_nc_fromB_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_nc_fromB_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_nc_fromB_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_nc_fromB_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_nc_fromB_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_nc_fromB_";
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
    
    hd0D0pt_nc_fromB_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_nc_fromB_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_nc_fromB_pm[i]->SetYTitle("Entries");
    flist_NoCuts_FromB->Add(hd0D0pt_nc_fromB_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_nc_fromB_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_nc_fromB_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_nc_fromB_pm[i]->SetYTitle("Entries");
    flist_NoCuts_FromB->Add(hMCd0D0pt_nc_fromB_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_nc_fromB_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_nc_fromB_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_nc_fromB_pm[i]->SetYTitle("Entries");
    flist_NoCuts_FromB->Add(hd0D0VtxTruept_nc_fromB_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_nc_fromB_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_nc_fromB_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_nc_fromB_sb[i]->SetYTitle("Entries");
    flist_NoCuts_FromB->Add(hd0D0pt_nc_fromB_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_nc_fromB_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_nc_fromB_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_nc_fromB_sb[i]->SetYTitle("Entries");
    flist_NoCuts_FromB->Add(hMCd0D0pt_nc_fromB_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_nc_fromB_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_nc_fromB_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_nc_fromB_sb[i]->SetYTitle("Entries");
    flist_NoCuts_FromB->Add(hd0D0VtxTruept_nc_fromB_sb[i]);
  }



  //############ NO CUTS FROM DSTAR HISTOGRAMS ###########
  //
  //#############  global properties histos #######

  TH2F *hCPtaVSd0d0_nc_fromDstar=new TH2F("hCPtaVSd0d0_nc_fromDstar","hCPtaVSd0d0_NoCuts_FromDStar",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_nc_fromDstar=new TH1F("hSecVtxZ_nc_fromDstar","hSecVtxZ_NoCuts_FromDStar",1000,-8.,8.);
  TH1F *hSecVtxX_nc_fromDstar=new TH1F("hSecVtxX_nc_fromDstar","hSecVtxX_NoCuts_FromDStar",1000,-3000.,3000.);
  TH1F *hSecVtxY_nc_fromDstar=new TH1F("hSecVtxY_nc_fromDstar","hSecVtxY_NoCuts_FromDStar",1000,-3000.,3000.);
  TH2F *hSecVtxXY_nc_fromDstar=new TH2F("hSecVtxXY_nc_fromDstar","hSecVtxXY_NoCuts_FromDStar",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_nc_fromDstar=new TH1F("hSecVtxPhi_nc_fromDstar","hSecVtxPhi_NoCuts_FromDStar",180,-180.1,180.1);
  TH1F *hCPta_nc_fromDstar=new TH1F("hCPta_nc_fromDstar","hCPta_NoCuts_FromDStar",100,0.,1.);
  TH1F *hd0xd0_nc_fromDstar=new TH1F("hd0xd0_nc_fromDstar","hd0xd0_NoCuts_FromDStar",1000,-100000.,100000.);
  TH1F *hMassTrue_nc_fromDstar=new TH1F("hMassTrue_nc_fromDstar","D^{0} MC inv. Mass No Cuts FromDStar(All momenta)",600,1.600,2.200);
  TH1F *hMass_nc_fromDstar=new TH1F("hMass_nc_fromDstar","D^{0} inv. Mass No Cuts FromDStar (All momenta)",600,1.600,2.200);
  hMass_nc_fromDstar->Sumw2();
  TH1F *hMassTrue_nc_fromDstar_pm=new TH1F("hMassTrue_nc_fromDstar_pm","D^{0} MC inv. Mass No Cuts FromDStar, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_nc_fromDstar_pm=new TH1F("hMass_nc_fromDstar_pm","D^{0} inv. Mass No Cuts FromDStar (All momenta), MassPeak",600,1.600,2.200);
  hMass_nc_fromDstar_pm->Sumw2();
  TH1F *hMassTrue_SB_nc_fromDstar=new TH1F("hMassTrue_nc_fromDstar_sb","D^{0} MC inv. Mass in Side Bands No Cuts FromDStar(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_nc_fromDstar=new TH1F("hMass_nc_fromDstar_sb","D^{0} inv. Mass in Side Bands No Cuts FromDStar (All momenta)",600,1.600,2.200);
  hMass_SB_nc_fromDstar->Sumw2();

  flist_NoCuts_FromDstar->Add(hCPtaVSd0d0_nc_fromDstar);
  flist_NoCuts_FromDstar->Add(hSecVtxZ_nc_fromDstar);
  flist_NoCuts_FromDstar->Add(hSecVtxY_nc_fromDstar);
  flist_NoCuts_FromDstar->Add(hSecVtxX_nc_fromDstar);
  flist_NoCuts_FromDstar->Add(hSecVtxXY_nc_fromDstar);
  flist_NoCuts_FromDstar->Add(hSecVtxPhi_nc_fromDstar);
  flist_NoCuts_FromDstar->Add(hCPta_nc_fromDstar);
  flist_NoCuts_FromDstar->Add(hd0xd0_nc_fromDstar);
  flist_NoCuts_FromDstar->Add(hMassTrue_nc_fromDstar);
  flist_NoCuts_FromDstar->Add(hMass_nc_fromDstar);
  flist_NoCuts_FromDstar->Add(hMassTrue_nc_fromDstar_pm);
  flist_NoCuts_FromDstar->Add(hMass_nc_fromDstar_pm);
  flist_NoCuts_FromDstar->Add(hMassTrue_SB_nc_fromDstar);
  flist_NoCuts_FromDstar->Add(hMass_SB_nc_fromDstar);

  //########## d0 D0 histos #############  
  TH1F *hd0D0_nc_fromDst_pm = new TH1F("hd0D0_nc_fromDstar_pm","D^{0} impact par. plot , No Cuts ,FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_nc_fromDst_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_nc_fromDst_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_nc_fromDst_pm = new TH1F("hd0D0VtxTrue_nc_fromDstar_pm","D^{0} impact par. w.r.t. True Vtx, No Cuts, FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_nc_fromDst_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_nc_fromDst_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_nc_fromDst_pm = new TH1F("hMCd0D0_nc_fromDstar_pm","D^{0} impact par. plot, No Cuts, FromDStar,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_nc_fromDst_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_nc_fromDst_pm->SetYTitle("Entries");

  TH1F *hd0D0_nc_fromDst_sb = new TH1F("hd0D0_nc_fromDstar_sb","D^{0} impact par. plot , No Cuts ,FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_nc_fromDst_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_nc_fromDst_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_nc_fromDst_sb = new TH1F("hd0D0VtxTrue_nc_fromDstar_sb","D^{0} impact par. w.r.t. True Vtx, No Cuts, FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_nc_fromDst_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_nc_fromDst_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_nc_fromDst_sb = new TH1F("hMCd0D0_nc_fromDstar_sb","D^{0} impact par. plot, No Cuts, FromDStar,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_nc_fromDst_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_nc_fromDst_sb->SetYTitle("Entries");

  flist_NoCuts_FromDstar->Add(hd0D0_nc_fromDst_pm);
  flist_NoCuts_FromDstar->Add(hd0D0VtxTrue_nc_fromDst_pm);
  flist_NoCuts_FromDstar->Add(hMCd0D0_nc_fromDst_pm);
  flist_NoCuts_FromDstar->Add(hd0D0_nc_fromDst_sb);
  flist_NoCuts_FromDstar->Add(hd0D0VtxTrue_nc_fromDst_sb);
  flist_NoCuts_FromDstar->Add(hMCd0D0_nc_fromDst_sb);
  
  TH1F **hd0D0pt_nc_fromDst_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_nc_fromDst_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_nc_fromDst_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_nc_fromDst_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_nc_fromDst_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_nc_fromDst_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_nc_fromDstar_";
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
    
    hd0D0pt_nc_fromDst_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_nc_fromDst_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_nc_fromDst_pm[i]->SetYTitle("Entries");
    flist_NoCuts_FromDstar->Add(hd0D0pt_nc_fromDst_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_nc_fromDst_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_nc_fromDst_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_nc_fromDst_pm[i]->SetYTitle("Entries");
    flist_NoCuts_FromDstar->Add(hMCd0D0pt_nc_fromDst_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_nc_fromDst_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_nc_fromDst_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_nc_fromDst_pm[i]->SetYTitle("Entries");
    flist_NoCuts_FromDstar->Add(hd0D0VtxTruept_nc_fromDst_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_nc_fromDst_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_nc_fromDst_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_nc_fromDst_sb[i]->SetYTitle("Entries");
    flist_NoCuts_FromDstar->Add(hd0D0pt_nc_fromDst_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_nc_fromDst_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_nc_fromDst_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_nc_fromDst_sb[i]->SetYTitle("Entries");
    flist_NoCuts_FromDstar->Add(hMCd0D0pt_nc_fromDst_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_nc_fromDst_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_nc_fromDst_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_nc_fromDst_sb[i]->SetYTitle("Entries");
    flist_NoCuts_FromDstar->Add(hd0D0VtxTruept_nc_fromDst_sb[i]);
  }


  //############ NO CUTS OTHER HISTOGRAMS ###########
  //
  //########### global properties histos ###########

  TH2F *hCPtaVSd0d0_nc_other=new TH2F("hCPtaVSd0d0_nc_other","hCPtaVSd0d0_NoCuts_other",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_nc_other=new TH1F("hSecVtxZ_nc_other","hSecVtxZ_NoCuts_other",1000,-8.,8.);
  TH1F *hSecVtxX_nc_other=new TH1F("hSecVtxX_nc_other","hSecVtxX_NoCuts_other",1000,-3000.,3000.);
  TH1F *hSecVtxY_nc_other=new TH1F("hSecVtxY_nc_other","hSecVtxY_NoCuts_other",1000,-3000.,3000.);
  TH2F *hSecVtxXY_nc_other=new TH2F("hSecVtxXY_nc_other","hSecVtxXY_NoCuts_other",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_nc_other=new TH1F("hSecVtxPhi_nc_other","hSecVtxPhi_NoCuts_other",180,-180.1,180.1);
  TH1F *hCPta_nc_other=new TH1F("hCPta_nc_other","hCPta_NoCuts_other",100,0.,1.);
  TH1F *hd0xd0_nc_other=new TH1F("hd0xd0_nc_other","hd0xd0_NoCuts_other",1000,-100000.,100000.);
  TH1F *hMassTrue_nc_other=new TH1F("hMassTrue_nc_other","D^{0} MC inv. Mass No Cuts other(All momenta)",600,1.600,2.200);
  TH1F *hMass_nc_other=new TH1F("hMass_nc_other","D^{0} inv. Mass No Cuts other (All momenta)",600,1.600,2.200);
  hMass_nc_other->Sumw2();
  TH1F *hMassTrue_nc_other_pm=new TH1F("hMassTrue_nc_other_pm","D^{0} MC inv. Mass No Cuts Other, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_nc_other_pm=new TH1F("hMass_nc_other_pm","D^{0} inv. Mass No Cuts Other (All momenta), MassPeak",600,1.600,2.200);
  hMass_nc_other_pm->Sumw2();
  TH1F *hMassTrue_SB_nc_other=new TH1F("hMassTrue_nc_other_sb","D^{0} MC inv. Mass in Side Bands No Cuts other(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_nc_other=new TH1F("hMass_nc_other_sb","D^{0} inv. Mass in Side Bands No Cuts other (All momenta)",600,1.600,2.200);
  hMass_SB_nc_other->Sumw2();

  flist_NoCuts_Other->Add(hCPtaVSd0d0_nc_other);
  flist_NoCuts_Other->Add(hSecVtxZ_nc_other);
  flist_NoCuts_Other->Add(hSecVtxY_nc_other);
  flist_NoCuts_Other->Add(hSecVtxX_nc_other);
  flist_NoCuts_Other->Add(hSecVtxXY_nc_other);
  flist_NoCuts_Other->Add(hSecVtxPhi_nc_other);
  flist_NoCuts_Other->Add(hCPta_nc_other);
  flist_NoCuts_Other->Add(hd0xd0_nc_other);
  flist_NoCuts_Other->Add(hMassTrue_nc_other);
  flist_NoCuts_Other->Add(hMass_nc_other);
  flist_NoCuts_Other->Add(hMassTrue_nc_other_pm);
  flist_NoCuts_Other->Add(hMass_nc_other_pm);
  flist_NoCuts_Other->Add(hMassTrue_SB_nc_other);
  flist_NoCuts_Other->Add(hMass_SB_nc_other);

  //############# d0 D0 histos ###############à
  TH1F *hd0D0_nc_other_pm = new TH1F("hd0D0_nc_other_pm","D^{0} impact par. plot , No Cuts ,Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_nc_other_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_nc_other_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_nc_other_pm = new TH1F("hd0D0VtxTrue_nc_other_pm","D^{0} impact par. w.r.t. True Vtx, No Cuts, Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_nc_other_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_nc_other_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_nc_other_pm = new TH1F("hMCd0D0_nc_other_pm","D^{0} impact par. plot, No Cuts, Other,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_nc_other_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_nc_other_pm->SetYTitle("Entries");

  TH1F *hd0D0_nc_other_sb = new TH1F("hd0D0_nc_other_sb","D^{0} impact par. plot , No Cuts ,Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_nc_other_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_nc_other_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_nc_other_sb = new TH1F("hd0D0VtxTrue_nc_other_sb","D^{0} impact par. w.r.t. True Vtx, No Cuts, Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_nc_other_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_nc_other_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_nc_other_sb = new TH1F("hMCd0D0_nc_other_sb","D^{0} impact par. plot, No Cuts, Other,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_nc_other_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_nc_other_sb->SetYTitle("Entries");

  flist_NoCuts_Other->Add(hd0D0_nc_other_pm);
  flist_NoCuts_Other->Add(hd0D0VtxTrue_nc_other_pm);
  flist_NoCuts_Other->Add(hMCd0D0_nc_other_pm);
  flist_NoCuts_Other->Add(hd0D0_nc_other_sb);
  flist_NoCuts_Other->Add(hd0D0VtxTrue_nc_other_sb);
  flist_NoCuts_Other->Add(hMCd0D0_nc_other_sb);
  
  TH1F **hd0D0pt_nc_other_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_nc_other_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_nc_other_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_nc_other_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_nc_other_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_nc_other_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_nc_other_";
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
    
    hd0D0pt_nc_other_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_nc_other_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_nc_other_pm[i]->SetYTitle("Entries");
    flist_NoCuts_Other->Add(hd0D0pt_nc_other_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_nc_other_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_nc_other_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_nc_other_pm[i]->SetYTitle("Entries");
    flist_NoCuts_Other->Add(hMCd0D0pt_nc_other_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_nc_other_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_nc_other_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_nc_other_pm[i]->SetYTitle("Entries");
    flist_NoCuts_Other->Add(hd0D0VtxTruept_nc_other_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_nc_other_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_nc_other_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_nc_other_sb[i]->SetYTitle("Entries");
    flist_NoCuts_Other->Add(hd0D0pt_nc_other_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_nc_other_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_nc_other_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_nc_other_sb[i]->SetYTitle("Entries");
    flist_NoCuts_Other->Add(hMCd0D0pt_nc_other_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_nc_other_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_nc_other_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_nc_other_sb[i]->SetYTitle("Entries");
    flist_NoCuts_Other->Add(hd0D0VtxTruept_nc_other_sb[i]);
  }


  //################################################################################################
  //                                                                                               #
  //                         HISTOS FOR LOOSE CUTS                                                 #
  //                                                                                               #
  //################################################################################################

  //############ LOOSE CUTS SIGNAL HISTOGRAMS ###############
  //
  // ####### global properties histo ############

  TH2F *hCPtaVSd0d0_ls_sign=new TH2F("hCPtaVSd0d0_ls_sign","hCPtaVSd0d0_LooseCuts_Signal",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_ls_sign=new TH1F("hSecVtxZ_ls_sign","hSecVtxZ_LooseCuts_Signal",1000,-8.,8.);
  TH1F *hSecVtxX_ls_sign=new TH1F("hSecVtxX_ls_sign","hSecVtxX_LooseCuts_Signal",1000,-3000.,3000.);
  TH1F *hSecVtxY_ls_sign=new TH1F("hSecVtxY_ls_sign","hSecVtxY_LooseCuts_Signal",1000,-3000.,3000.);
  TH2F *hSecVtxXY_ls_sign=new TH2F("hSecVtxXY_ls_sign","hSecVtxXY_LooseCuts_Signal",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_ls_sign=new TH1F("hSecVtxPhi_ls_sign","hSecVtxPhi_LooseCuts_Signal",180,-180.1,180.1);
  TH1F *hCPta_ls_sign=new TH1F("hCPta_ls_sign","hCPta_LooseCuts_Signal",100,0.,1.);
  TH1F *hd0xd0_ls_sign=new TH1F("hd0xd0_ls_sign","hd0xd0_LooseCuts_Signal",1000,-100000.,100000.);
  TH1F *hMassTrue_ls_sign=new TH1F("hMassTrue_ls_sign","D^{0} MC inv. Mass Loose Cuts Signal(All momenta)",600,1.600,2.200);
  TH1F *hMass_ls_sign=new TH1F("hMass_ls_sign","D^{0} inv. Mass Loose Cuts Signal (All momenta)",600,1.600,2.200);
  hMass_ls_sign->Sumw2();
  TH1F *hMassTrue_ls_sign_pm=new TH1F("hMassTrue_ls_sign_pm","D^{0} MC inv. Mass Loose Cuts Signal, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_ls_sign_pm=new TH1F("hMass_ls_sign_pm","D^{0} inv. Mass Loose Cuts Signal (All momenta), MassPeak",600,1.600,2.200);
  hMass_ls_sign_pm->Sumw2();
  TH1F *hMassTrue_SB_ls_sign=new TH1F("hMassTrue_ls_sign_sb","D^{0} MC inv. Mass in Side Bands Loose Cuts Signal(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_ls_sign=new TH1F("hMass_ls_sign_sb","D^{0} inv. Mass in Side Bands Loose Cuts Signal (All momenta)",600,1.600,2.200);
  hMass_SB_ls_sign->Sumw2();

  flist_LsCuts_Signal->Add(hCPtaVSd0d0_ls_sign);
  flist_LsCuts_Signal->Add(hSecVtxZ_ls_sign);
  flist_LsCuts_Signal->Add(hSecVtxY_ls_sign);
  flist_LsCuts_Signal->Add(hSecVtxX_ls_sign);
  flist_LsCuts_Signal->Add(hSecVtxXY_ls_sign);
  flist_LsCuts_Signal->Add(hSecVtxPhi_ls_sign);
  flist_LsCuts_Signal->Add(hCPta_ls_sign);
  flist_LsCuts_Signal->Add(hd0xd0_ls_sign);
  flist_LsCuts_Signal->Add(hMassTrue_ls_sign);
  flist_LsCuts_Signal->Add(hMass_ls_sign);
  flist_LsCuts_Signal->Add(hMassTrue_ls_sign_pm);
  flist_LsCuts_Signal->Add(hMass_ls_sign_pm);
  flist_LsCuts_Signal->Add(hMassTrue_SB_ls_sign);
  flist_LsCuts_Signal->Add(hMass_SB_ls_sign);

  // ####### d0 D0 histos ############
  TH1F *hd0D0_ls_sign_pm = new TH1F("hd0D0_ls_sign_pm","D^{0} impact par. plot , Loose Cuts ,Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_ls_sign_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_ls_sign_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_ls_sign_pm = new TH1F("hd0D0VtxTrue_ls_sign_pm","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_ls_sign_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_ls_sign_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_ls_sign_pm = new TH1F("hMCd0D0_ls_sign_pm","D^{0} impact par. plot, Loose Cuts, Signal,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_ls_sign_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_ls_sign_pm->SetYTitle("Entries");

  TH1F *hd0D0_ls_sign_sb = new TH1F("hd0D0_ls_sign_sb","D^{0} impact par. plot , Loose Cuts ,Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_ls_sign_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_ls_sign_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_ls_sign_sb = new TH1F("hd0D0VtxTrue_ls_sign_sb","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_ls_sign_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_ls_sign_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_ls_sign_sb = new TH1F("hMCd0D0_ls_sign_sb","D^{0} impact par. plot, Loose Cuts, Signal,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_ls_sign_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_ls_sign_sb->SetYTitle("Entries");

  flist_LsCuts_Signal->Add(hd0D0_ls_sign_pm);
  flist_LsCuts_Signal->Add(hd0D0VtxTrue_ls_sign_pm);
  flist_LsCuts_Signal->Add(hMCd0D0_ls_sign_pm);
  flist_LsCuts_Signal->Add(hd0D0_ls_sign_sb);
  flist_LsCuts_Signal->Add(hd0D0VtxTrue_ls_sign_sb);
  flist_LsCuts_Signal->Add(hMCd0D0_ls_sign_sb);
  
  TH1F **hd0D0pt_ls_sign_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_ls_sign_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_ls_sign_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_ls_sign_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_ls_sign_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_ls_sign_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_ls_sign_";
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
    
    hd0D0pt_ls_sign_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_ls_sign_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_ls_sign_pm[i]->SetYTitle("Entries");
    flist_LsCuts_Signal->Add(hd0D0pt_ls_sign_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_ls_sign_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_ls_sign_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_ls_sign_pm[i]->SetYTitle("Entries");
    flist_LsCuts_Signal->Add(hMCd0D0pt_ls_sign_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_ls_sign_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_ls_sign_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_ls_sign_pm[i]->SetYTitle("Entries");
    flist_LsCuts_Signal->Add(hd0D0VtxTruept_ls_sign_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_ls_sign_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_ls_sign_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_ls_sign_sb[i]->SetYTitle("Entries");
    flist_LsCuts_Signal->Add(hd0D0pt_ls_sign_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_ls_sign_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_ls_sign_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_ls_sign_sb[i]->SetYTitle("Entries");
    flist_LsCuts_Signal->Add(hMCd0D0pt_ls_sign_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_ls_sign_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_ls_sign_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_ls_sign_sb[i]->SetYTitle("Entries");
    flist_LsCuts_Signal->Add(hd0D0VtxTruept_ls_sign_sb[i]);
  }


  //############ LOOSE CUTS BACKGROUND HISTOGRAMS ###########
  //
  //   ######## global properties histos #######
  TH2F *hCPtaVSd0d0_ls_back=new TH2F("hCPtaVSd0d0_ls_back","hCPtaVSd0d0_LooseCuts_Background",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_ls_back=new TH1F("hSecVtxZ_ls_back","hSecVtxZ_LooseCuts_Background",1000,-8.,8.);
  TH1F *hSecVtxX_ls_back=new TH1F("hSecVtxX_ls_back","hSecVtxX_LooseCuts_Background",1000,-3000.,3000.);
  TH1F *hSecVtxY_ls_back=new TH1F("hSecVtxY_ls_back","hSecVtxY_LooseCuts_Background",1000,-3000.,3000.);
  TH2F *hSecVtxXY_ls_back=new TH2F("hSecVtxXY_ls_back","hSecVtxXY_LooseCuts_Background",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_ls_back=new TH1F("hSecVtxPhi_ls_back","hSecVtxPhi_LooseCuts_Background",180,-180.1,180.1);
  TH1F *hCPta_ls_back=new TH1F("hCPta_ls_back","hCPta_LooseCuts_Background",100,0.,1.);
  TH1F *hd0xd0_ls_back=new TH1F("hd0xd0_ls_back","hd0xd0_LooseCuts_Background",1000,-100000.,100000.);
  TH1F *hMassTrue_ls_back=new TH1F("hMassTrue_ls_back","D^{0} MC inv. Mass Loose Cuts Background(All momenta)",600,1.600,2.200);
  TH1F *hMass_ls_back=new TH1F("hMass_ls_back","D^{0} inv. Mass Loose Cuts Background (All momenta)",600,1.600,2.200);
  hMass_ls_back->Sumw2();
  TH1F *hMassTrue_ls_back_pm=new TH1F("hMassTrue_ls_back_pm","D^{0} MC inv. Mass Loose Cuts Background, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_ls_back_pm=new TH1F("hMass_ls_back_pm","D^{0} inv. Mass Loose Cuts Background (All momenta), MassPeak",600,1.600,2.200);
  hMass_ls_back_pm->Sumw2();
  TH1F *hMassTrue_SB_ls_back=new TH1F("hMassTrue_ls_back_sb","D^{0} MC inv. Mass in Side Bands Loose Cuts Background(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_ls_back=new TH1F("hMass_ls_back_sb","D^{0} inv. Mass in Side Bands Loose Cuts Background (All momenta)",600,1.600,2.200);
  hMass_SB_ls_back->Sumw2();

  flist_LsCuts_Back->Add(hCPtaVSd0d0_ls_back);
  flist_LsCuts_Back->Add(hSecVtxZ_ls_back);
  flist_LsCuts_Back->Add(hSecVtxY_ls_back);
  flist_LsCuts_Back->Add(hSecVtxX_ls_back);
  flist_LsCuts_Back->Add(hSecVtxXY_ls_back);
  flist_LsCuts_Back->Add(hSecVtxPhi_ls_back);
  flist_LsCuts_Back->Add(hCPta_ls_back);
  flist_LsCuts_Back->Add(hd0xd0_ls_back);
  flist_LsCuts_Back->Add(hMassTrue_ls_back);
  flist_LsCuts_Back->Add(hMass_ls_back);
  flist_LsCuts_Back->Add(hMassTrue_ls_back_pm);
  flist_LsCuts_Back->Add(hMass_ls_back_pm);
  flist_LsCuts_Back->Add(hMassTrue_SB_ls_back);
  flist_LsCuts_Back->Add(hMass_SB_ls_back);


  // ####### d0 D0 histos ############
  
 TH1F *hd0D0_ls_back_pm = new TH1F("hd0D0_ls_back_pm","D^{0} impact par. plot , Loose Cuts ,Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_ls_back_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_ls_back_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_ls_back_pm = new TH1F("hd0D0VtxTrue_ls_back_pm","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_ls_back_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_ls_back_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_ls_back_pm = new TH1F("hMCd0D0_ls_back_pm","D^{0} impact par. plot, Loose Cuts, Background,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_ls_back_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_ls_back_pm->SetYTitle("Entries");

  TH1F *hd0D0_ls_back_sb = new TH1F("hd0D0_ls_back_sb","D^{0} impact par. plot , Loose Cuts ,Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_ls_back_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_ls_back_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_ls_back_sb = new TH1F("hd0D0VtxTrue_ls_back_sb","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_ls_back_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_ls_back_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_ls_back_sb = new TH1F("hMCd0D0_ls_back_sb","D^{0} impact par. plot, Loose Cuts, Background,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_ls_back_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_ls_back_sb->SetYTitle("Entries");

  flist_LsCuts_Back->Add(hd0D0_ls_back_pm);
  flist_LsCuts_Back->Add(hd0D0VtxTrue_ls_back_pm);
  flist_LsCuts_Back->Add(hMCd0D0_ls_back_pm);
  flist_LsCuts_Back->Add(hd0D0_ls_back_sb);
  flist_LsCuts_Back->Add(hd0D0VtxTrue_ls_back_sb);
  flist_LsCuts_Back->Add(hMCd0D0_ls_back_sb);
  
  TH1F **hd0D0pt_ls_back_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_ls_back_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_ls_back_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_ls_back_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_ls_back_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_ls_back_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_ls_back_";
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
    
    hd0D0pt_ls_back_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_ls_back_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_ls_back_pm[i]->SetYTitle("Entries");
    flist_LsCuts_Back->Add(hd0D0pt_ls_back_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_ls_back_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_ls_back_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_ls_back_pm[i]->SetYTitle("Entries");
    flist_LsCuts_Back->Add(hMCd0D0pt_ls_back_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_ls_back_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_ls_back_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_ls_back_pm[i]->SetYTitle("Entries");
    flist_LsCuts_Back->Add(hd0D0VtxTruept_ls_back_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_ls_back_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_ls_back_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_ls_back_sb[i]->SetYTitle("Entries");
    flist_LsCuts_Back->Add(hd0D0pt_ls_back_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_ls_back_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_ls_back_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_ls_back_sb[i]->SetYTitle("Entries");
    flist_LsCuts_Back->Add(hMCd0D0pt_ls_back_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_ls_back_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_ls_back_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_ls_back_sb[i]->SetYTitle("Entries");
    flist_LsCuts_Back->Add(hd0D0VtxTruept_ls_back_sb[i]);
  }



 //############ LOOSE CUTS FROMB HISTOGRAMS ###########
  //
  //#######  global properties histos

  TH2F *hCPtaVSd0d0_ls_fromB=new TH2F("hCPtaVSd0d0_ls_fromB","hCPtaVSd0d0_LooseCuts_FromB",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_ls_fromB=new TH1F("hSecVtxZ_ls_fromB","hSecVtxZ_LooseCuts_FromB",1000,-8.,8.);
  TH1F *hSecVtxX_ls_fromB=new TH1F("hSecVtxX_ls_fromB","hSecVtxX_LooseCuts_FromB",1000,-3000.,3000.);
  TH1F *hSecVtxY_ls_fromB=new TH1F("hSecVtxY_ls_fromB","hSecVtxY_LooseCuts_FromB",1000,-3000.,3000.);
  TH2F *hSecVtxXY_ls_fromB=new TH2F("hSecVtxXY_ls_fromB","hSecVtxXY_LooseCuts_FromB",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_ls_fromB=new TH1F("hSecVtxPhi_ls_fromB","hSecVtxPhi_LooseCuts_FromB",180,-180.1,180.1);
  TH1F *hCPta_ls_fromB=new TH1F("hCPta_ls_fromB","hCPta_LooseCuts_FromB",100,0.,1.);
  TH1F *hd0xd0_ls_fromB=new TH1F("hd0xd0_ls_fromB","hd0xd0_LooseCuts_FromB",1000,-100000.,100000.);
  TH1F *hMassTrue_ls_fromB=new TH1F("hMassTrue_ls_fromB","D^{0} MC inv. Mass Loose Cuts FromB(All momenta)",600,1.600,2.200);
  TH1F *hMass_ls_fromB=new TH1F("hMass_ls_fromB","D^{0} inv. Mass Loose Cuts FromB (All momenta)",600,1.600,2.200);
  hMass_ls_fromB->Sumw2();
  TH1F *hMassTrue_ls_fromB_pm=new TH1F("hMassTrue_ls_fromB_pm","D^{0} MC inv. Mass Loose Cuts FromB, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_ls_fromB_pm=new TH1F("hMass_ls_fromB_pm","D^{0} inv. Mass Loose Cuts FromB (All momenta), MassPeak",600,1.600,2.200);
  hMass_ls_fromB_pm->Sumw2();
  TH1F *hMassTrue_SB_ls_fromB=new TH1F("hMassTrue_ls_fromB_sb","D^{0} MC inv. Mass in Side Bands Loose Cuts FromB(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_ls_fromB=new TH1F("hMass_ls_fromB_sb","D^{0} inv. Mass in Side Bands Loose Cuts FromB (All momenta)",600,1.600,2.200);
  hMass_SB_ls_fromB->Sumw2();

  flist_LsCuts_FromB->Add(hCPtaVSd0d0_ls_fromB);
  flist_LsCuts_FromB->Add(hSecVtxZ_ls_fromB);
  flist_LsCuts_FromB->Add(hSecVtxY_ls_fromB);
  flist_LsCuts_FromB->Add(hSecVtxX_ls_fromB);
  flist_LsCuts_FromB->Add(hSecVtxXY_ls_fromB);
  flist_LsCuts_FromB->Add(hSecVtxPhi_ls_fromB);
  flist_LsCuts_FromB->Add(hCPta_ls_fromB);
  flist_LsCuts_FromB->Add(hd0xd0_ls_fromB);
  flist_LsCuts_FromB->Add(hMassTrue_ls_fromB);
  flist_LsCuts_FromB->Add(hMass_ls_fromB);
  flist_LsCuts_FromB->Add(hMassTrue_ls_fromB_pm);
  flist_LsCuts_FromB->Add(hMass_ls_fromB_pm);
  flist_LsCuts_FromB->Add(hMassTrue_SB_ls_fromB);
  flist_LsCuts_FromB->Add(hMass_SB_ls_fromB);

  // ######### d0 D0 histos ##############
  TH1F *hd0D0_ls_fromB_pm = new TH1F("hd0D0_ls_fromB_pm","D^{0} impact par. plot , Loose Cuts ,FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_ls_fromB_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_ls_fromB_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_ls_fromB_pm = new TH1F("hd0D0VtxTrue_ls_fromB_pm","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_ls_fromB_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_ls_fromB_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_ls_fromB_pm = new TH1F("hMCd0D0_ls_fromB_pm","D^{0} impact par. plot, Loose Cuts, FromB,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_ls_fromB_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_ls_fromB_pm->SetYTitle("Entries");

  TH1F *hd0D0_ls_fromB_sb = new TH1F("hd0D0_ls_fromB_sb","D^{0} impact par. plot , Loose Cuts ,FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_ls_fromB_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_ls_fromB_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_ls_fromB_sb = new TH1F("hd0D0VtxTrue_ls_fromB_sb","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_ls_fromB_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_ls_fromB_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_ls_fromB_sb = new TH1F("hMCd0D0_ls_fromB_sb","D^{0} impact par. plot, Loose Cuts, FromB,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_ls_fromB_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_ls_fromB_sb->SetYTitle("Entries");

  flist_LsCuts_FromB->Add(hd0D0_ls_fromB_pm);
  flist_LsCuts_FromB->Add(hd0D0VtxTrue_ls_fromB_pm);
  flist_LsCuts_FromB->Add(hMCd0D0_ls_fromB_pm);
  flist_LsCuts_FromB->Add(hd0D0_ls_fromB_sb);
  flist_LsCuts_FromB->Add(hd0D0VtxTrue_ls_fromB_sb);
  flist_LsCuts_FromB->Add(hMCd0D0_ls_fromB_sb);
  
  TH1F **hd0D0pt_ls_fromB_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_ls_fromB_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_ls_fromB_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_ls_fromB_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_ls_fromB_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_ls_fromB_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_ls_fromB_";
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
    
    hd0D0pt_ls_fromB_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_ls_fromB_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_ls_fromB_pm[i]->SetYTitle("Entries");
    flist_LsCuts_FromB->Add(hd0D0pt_ls_fromB_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_ls_fromB_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_ls_fromB_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_ls_fromB_pm[i]->SetYTitle("Entries");
    flist_LsCuts_FromB->Add(hMCd0D0pt_ls_fromB_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_ls_fromB_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_ls_fromB_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_ls_fromB_pm[i]->SetYTitle("Entries");
    flist_LsCuts_FromB->Add(hd0D0VtxTruept_ls_fromB_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_ls_fromB_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_ls_fromB_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_ls_fromB_sb[i]->SetYTitle("Entries");
    flist_LsCuts_FromB->Add(hd0D0pt_ls_fromB_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_ls_fromB_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_ls_fromB_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_ls_fromB_sb[i]->SetYTitle("Entries");
    flist_LsCuts_FromB->Add(hMCd0D0pt_ls_fromB_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_ls_fromB_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_ls_fromB_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_ls_fromB_sb[i]->SetYTitle("Entries");
    flist_LsCuts_FromB->Add(hd0D0VtxTruept_ls_fromB_sb[i]);
  }



 //############ LOOSE CUTS FROM DSTAR HISTOGRAMS ###########
 //
  //############## global properties histos
  TH2F *hCPtaVSd0d0_ls_fromDstar=new TH2F("hCPtaVSd0d0_ls_fromDstar","hCPtaVSd0d0_LooseCuts_FromDStar",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_ls_fromDstar=new TH1F("hSecVtxZ_ls_fromDstar","hSecVtxZ_LooseCuts_FromDStar",1000,-8.,8.);
  TH1F *hSecVtxX_ls_fromDstar=new TH1F("hSecVtxX_ls_fromDstar","hSecVtxX_LooseCuts_FromDStar",1000,-3000.,3000.);
  TH1F *hSecVtxY_ls_fromDstar=new TH1F("hSecVtxY_ls_fromDstar","hSecVtxY_LooseCuts_FromDStar",1000,-3000.,3000.);
  TH2F *hSecVtxXY_ls_fromDstar=new TH2F("hSecVtxXY_ls_fromDstar","hSecVtxXY_LooseCuts_FromDStar",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_ls_fromDstar=new TH1F("hSecVtxPhi_ls_fromDstar","hSecVtxPhi_LooseCuts_FromDStar",180,-180.1,180.1);
  TH1F *hCPta_ls_fromDstar=new TH1F("hCPta_ls_fromDstar","hCPta_LooseCuts_FromDStar",100,0.,1.);
  TH1F *hd0xd0_ls_fromDstar=new TH1F("hd0xd0_ls_fromDstar","hd0xd0_LooseCuts_FromDStar",1000,-100000.,100000.);
  TH1F *hMassTrue_ls_fromDstar=new TH1F("hMassTrue_ls_fromDstar","D^{0} MC inv. Mass Loose Cuts FromDStar(All momenta)",600,1.600,2.200);
  TH1F *hMass_ls_fromDstar=new TH1F("hMass_ls_fromDstar","D^{0} inv. Mass Loose Cuts FromDStar (All momenta)",600,1.600,2.200);
  hMass_ls_fromDstar->Sumw2();
  TH1F *hMassTrue_ls_fromDstar_pm=new TH1F("hMassTrue_ls_fromDstar_pm","D^{0} MC inv. Mass Loose Cuts FromDStar, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_ls_fromDstar_pm=new TH1F("hMass_ls_fromDstar_pm","D^{0} inv. Mass Loose Cuts FromDStar (All momenta), MassPeak",600,1.600,2.200);
  hMass_ls_fromDstar_pm->Sumw2();
  TH1F *hMassTrue_SB_ls_fromDstar=new TH1F("hMassTrue_ls_fromDstar_sb","D^{0} MC inv. Mass in Side Bands Loose Cuts FromDStar(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_ls_fromDstar=new TH1F("hMass_ls_fromDstar_sb","D^{0} inv. Mass in Side Bands Loose Cuts FromDStar (All momenta)",600,1.600,2.200);
  hMass_SB_ls_fromDstar->Sumw2();

  flist_LsCuts_FromDstar->Add(hCPtaVSd0d0_ls_fromDstar);
  flist_LsCuts_FromDstar->Add(hSecVtxZ_ls_fromDstar);
  flist_LsCuts_FromDstar->Add(hSecVtxY_ls_fromDstar);
  flist_LsCuts_FromDstar->Add(hSecVtxX_ls_fromDstar);
  flist_LsCuts_FromDstar->Add(hSecVtxXY_ls_fromDstar);
  flist_LsCuts_FromDstar->Add(hSecVtxPhi_ls_fromDstar);
  flist_LsCuts_FromDstar->Add(hCPta_ls_fromDstar);
  flist_LsCuts_FromDstar->Add(hd0xd0_ls_fromDstar);
  flist_LsCuts_FromDstar->Add(hMassTrue_ls_fromDstar);
  flist_LsCuts_FromDstar->Add(hMass_ls_fromDstar);
 flist_LsCuts_FromDstar->Add(hMassTrue_ls_fromDstar_pm);
  flist_LsCuts_FromDstar->Add(hMass_ls_fromDstar_pm);
  flist_LsCuts_FromDstar->Add(hMassTrue_SB_ls_fromDstar);
  flist_LsCuts_FromDstar->Add(hMass_SB_ls_fromDstar);

  //########## d0 D0 histos #############  
  TH1F *hd0D0_ls_fromDst_pm = new TH1F("hd0D0_ls_fromDstar_pm","D^{0} impact par. plot , Loose Cuts ,FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_ls_fromDst_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_ls_fromDst_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_ls_fromDst_pm = new TH1F("hd0D0VtxTrue_ls_fromDstar_pm","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_ls_fromDst_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_ls_fromDst_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_ls_fromDst_pm = new TH1F("hMCd0D0_ls_fromDstar_pm","D^{0} impact par. plot, Loose Cuts, FromDStar,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_ls_fromDst_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_ls_fromDst_pm->SetYTitle("Entries");

  TH1F *hd0D0_ls_fromDst_sb = new TH1F("hd0D0_ls_fromDstar_sb","D^{0} impact par. plot , Loose Cuts ,FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_ls_fromDst_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_ls_fromDst_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_ls_fromDst_sb = new TH1F("hd0D0VtxTrue_ls_fromDstar_sb","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_ls_fromDst_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_ls_fromDst_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_ls_fromDst_sb = new TH1F("hMCd0D0_ls_fromDstar_sb","D^{0} impact par. plot, Loose Cuts, FromDStar,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_ls_fromDst_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_ls_fromDst_sb->SetYTitle("Entries");

  flist_LsCuts_FromDstar->Add(hd0D0_ls_fromDst_pm);
  flist_LsCuts_FromDstar->Add(hd0D0VtxTrue_ls_fromDst_pm);
  flist_LsCuts_FromDstar->Add(hMCd0D0_ls_fromDst_pm);
  flist_LsCuts_FromDstar->Add(hd0D0_ls_fromDst_sb);
  flist_LsCuts_FromDstar->Add(hd0D0VtxTrue_ls_fromDst_sb);
  flist_LsCuts_FromDstar->Add(hMCd0D0_ls_fromDst_sb);
  
  TH1F **hd0D0pt_ls_fromDst_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_ls_fromDst_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_ls_fromDst_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_ls_fromDst_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_ls_fromDst_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_ls_fromDst_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_ls_fromDstar_";
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
    
    hd0D0pt_ls_fromDst_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_ls_fromDst_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_ls_fromDst_pm[i]->SetYTitle("Entries");
    flist_LsCuts_FromDstar->Add(hd0D0pt_ls_fromDst_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_ls_fromDst_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_ls_fromDst_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_ls_fromDst_pm[i]->SetYTitle("Entries");
    flist_LsCuts_FromDstar->Add(hMCd0D0pt_ls_fromDst_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_ls_fromDst_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_ls_fromDst_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_ls_fromDst_pm[i]->SetYTitle("Entries");
    flist_LsCuts_FromDstar->Add(hd0D0VtxTruept_ls_fromDst_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_ls_fromDst_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_ls_fromDst_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_ls_fromDst_sb[i]->SetYTitle("Entries");
    flist_LsCuts_FromDstar->Add(hd0D0pt_ls_fromDst_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_ls_fromDst_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_ls_fromDst_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_ls_fromDst_sb[i]->SetYTitle("Entries");
    flist_LsCuts_FromDstar->Add(hMCd0D0pt_ls_fromDst_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_ls_fromDst_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_ls_fromDst_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_ls_fromDst_sb[i]->SetYTitle("Entries");
    flist_LsCuts_FromDstar->Add(hd0D0VtxTruept_ls_fromDst_sb[i]);
  }


  //############ LOOSE CUTS OTHER HISTOGRAMS ###########
  //
  //########### global properties histos ###########

  TH2F *hCPtaVSd0d0_ls_other=new TH2F("hCPtaVSd0d0_ls_other","hCPtaVSd0d0_LooseCuts_other",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_ls_other=new TH1F("hSecVtxZ_ls_other","hSecVtxZ_LooseCuts_other",1000,-8.,8.);
  TH1F *hSecVtxX_ls_other=new TH1F("hSecVtxX_ls_other","hSecVtxX_LooseCuts_other",1000,-3000.,3000.);
  TH1F *hSecVtxY_ls_other=new TH1F("hSecVtxY_ls_other","hSecVtxY_LooseCuts_other",1000,-3000.,3000.);
  TH2F *hSecVtxXY_ls_other=new TH2F("hSecVtxXY_ls_other","hSecVtxXY_LooseCuts_other",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_ls_other=new TH1F("hSecVtxPhi_ls_other","hSecVtxPhi_LooseCuts_other",180,-180.1,180.1);
  TH1F *hCPta_ls_other=new TH1F("hCPta_ls_other","hCPta_LooseCuts_other",100,0.,1.);
  TH1F *hd0xd0_ls_other=new TH1F("hd0xd0_ls_other","hd0xd0_LooseCuts_other",1000,-100000.,100000.);
  TH1F *hMassTrue_ls_other=new TH1F("hMassTrue_ls_other","D^{0} MC inv. Mass Loose Cuts other(All momenta)",600,1.600,2.200);
  TH1F *hMass_ls_other=new TH1F("hMass_ls_other","D^{0} inv. Mass Loose Cuts other (All momenta)",600,1.600,2.200);
  hMass_ls_other->Sumw2();
  TH1F *hMassTrue_ls_other_pm=new TH1F("hMassTrue_ls_other_pm","D^{0} MC inv. Mass Loose Cuts other, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_ls_other_pm=new TH1F("hMass_ls_other_pm","D^{0} inv. Mass Loose Cuts other (All momenta), MassPeak",600,1.600,2.200);
  hMass_ls_other_pm->Sumw2();
  TH1F *hMassTrue_SB_ls_other=new TH1F("hMassTrue_ls_other_sb","D^{0} MC inv. Mass in Side Bands Loose Cuts other(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_ls_other=new TH1F("hMass_ls_other_sb","D^{0} inv. Mass in Side Bands Loose Cuts other (All momenta)",600,1.600,2.200);
  hMass_SB_ls_other->Sumw2();

  flist_LsCuts_Other->Add(hCPtaVSd0d0_ls_other);
  flist_LsCuts_Other->Add(hSecVtxZ_ls_other);
  flist_LsCuts_Other->Add(hSecVtxY_ls_other);
  flist_LsCuts_Other->Add(hSecVtxX_ls_other);
  flist_LsCuts_Other->Add(hSecVtxXY_ls_other);
  flist_LsCuts_Other->Add(hSecVtxPhi_ls_other);
  flist_LsCuts_Other->Add(hCPta_ls_other);
  flist_LsCuts_Other->Add(hd0xd0_ls_other);
  flist_LsCuts_Other->Add(hMassTrue_ls_other);
  flist_LsCuts_Other->Add(hMass_ls_other);
  flist_LsCuts_Other->Add(hMassTrue_ls_other_pm);
  flist_LsCuts_Other->Add(hMass_ls_other_pm);
  flist_LsCuts_Other->Add(hMassTrue_SB_ls_other);
  flist_LsCuts_Other->Add(hMass_SB_ls_other);

  //############# d0 D0 histos ###############à
  TH1F *hd0D0_ls_other_pm = new TH1F("hd0D0_ls_other_pm","D^{0} impact par. plot , Loose Cuts ,Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_ls_other_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_ls_other_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_ls_other_pm = new TH1F("hd0D0VtxTrue_ls_other_pm","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_ls_other_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_ls_other_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_ls_other_pm = new TH1F("hMCd0D0_ls_other_pm","D^{0} impact par. plot, Loose Cuts, Other,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_ls_other_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_ls_other_pm->SetYTitle("Entries");

  TH1F *hd0D0_ls_other_sb = new TH1F("hd0D0_ls_other_sb","D^{0} impact par. plot , Loose Cuts ,Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_ls_other_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_ls_other_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_ls_other_sb = new TH1F("hd0D0VtxTrue_ls_other_sb","D^{0} impact par. w.r.t. True Vtx, Loose Cuts, Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_ls_other_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_ls_other_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_ls_other_sb = new TH1F("hMCd0D0_ls_other_sb","D^{0} impact par. plot, Loose Cuts, Other,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_ls_other_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_ls_other_sb->SetYTitle("Entries");

  flist_LsCuts_Other->Add(hd0D0_ls_other_pm);
  flist_LsCuts_Other->Add(hd0D0VtxTrue_ls_other_pm);
  flist_LsCuts_Other->Add(hMCd0D0_ls_other_pm);
  flist_LsCuts_Other->Add(hd0D0_ls_other_sb);
  flist_LsCuts_Other->Add(hd0D0VtxTrue_ls_other_sb);
  flist_LsCuts_Other->Add(hMCd0D0_ls_other_sb);
  
  TH1F **hd0D0pt_ls_other_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_ls_other_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_ls_other_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_ls_other_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_ls_other_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_ls_other_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_ls_other_";
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
    
    hd0D0pt_ls_other_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_ls_other_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_ls_other_pm[i]->SetYTitle("Entries");
    flist_LsCuts_Other->Add(hd0D0pt_ls_other_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_ls_other_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_ls_other_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_ls_other_pm[i]->SetYTitle("Entries");
    flist_LsCuts_Other->Add(hMCd0D0pt_ls_other_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_ls_other_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_ls_other_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_ls_other_pm[i]->SetYTitle("Entries");
    flist_LsCuts_Other->Add(hd0D0VtxTruept_ls_other_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_ls_other_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_ls_other_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_ls_other_sb[i]->SetYTitle("Entries");
    flist_LsCuts_Other->Add(hd0D0pt_ls_other_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_ls_other_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_ls_other_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_ls_other_sb[i]->SetYTitle("Entries");
    flist_LsCuts_Other->Add(hMCd0D0pt_ls_other_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_ls_other_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_ls_other_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_ls_other_sb[i]->SetYTitle("Entries");
    flist_LsCuts_Other->Add(hd0D0VtxTruept_ls_other_sb[i]);
  }




  //################################################################################################
  //                                                                                               #
  //                         HISTOS FOR TIGHT CUTS                                                 #
  //                                                                                               #
  //################################################################################################

  //############ TIGHT CUTS SIGNAL HISTOGRAMS ###############
  //
  // ####### global properties histo ############

  TH2F *hCPtaVSd0d0_tgh_sign=new TH2F("hCPtaVSd0d0_tgh_sign","hCPtaVSd0d0_TightCuts_Signal",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_tgh_sign=new TH1F("hSecVtxZ_tgh_sign","hSecVtxZ_TightCuts_Signal",1000,-8.,8.);
  TH1F *hSecVtxX_tgh_sign=new TH1F("hSecVtxX_tgh_sign","hSecVtxX_TightCuts_Signal",1000,-3000.,3000.);
  TH1F *hSecVtxY_tgh_sign=new TH1F("hSecVtxY_tgh_sign","hSecVtxY_TightCuts_Signal",1000,-3000.,3000.);
  TH2F *hSecVtxXY_tgh_sign=new TH2F("hSecVtxXY_tgh_sign","hSecVtxXY_TightCuts_Signal",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_tgh_sign=new TH1F("hSecVtxPhi_tgh_sign","hSecVtxPhi_TightCuts_Signal",180,-180.1,180.1);
  TH1F *hCPta_tgh_sign=new TH1F("hCPta_tgh_sign","hCPta_TightCuts_Signal",100,0.,1.);
  TH1F *hd0xd0_tgh_sign=new TH1F("hd0xd0_tgh_sign","hd0xd0_TightCuts_Signal",1000,-100000.,100000.);
  TH1F *hMassTrue_tgh_sign=new TH1F("hMassTrue_tgh_sign","D^{0} MC inv. Mass Tight Cuts Signal(All momenta)",600,1.600,2.200);
  TH1F *hMass_tgh_sign=new TH1F("hMass_tgh_sign","D^{0} inv. Mass Tight Cuts Signal (All momenta)",600,1.600,2.200);
  hMass_tgh_sign->Sumw2();
  TH1F *hMassTrue_tgh_sign_pm=new TH1F("hMassTrue_tgh_sign_pm","D^{0} MC inv. Mass Tight Cuts Signal, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_tgh_sign_pm=new TH1F("hMass_tgh_sign_pm","D^{0} inv. Mass Tight Cuts Signal (All momenta), MassPeak",600,1.600,2.200);
  hMass_tgh_sign_pm->Sumw2();
  TH1F *hMassTrue_SB_tgh_sign=new TH1F("hMassTrue_tgh_sign_sb","D^{0} MC inv. Mass in Side Bands Tight Cuts Signal(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_tgh_sign=new TH1F("hMass_tgh_sign_sb","D^{0} inv. Mass in Side Bands Tight Cuts Signal (All momenta)",600,1.600,2.200);
  hMass_SB_tgh_sign->Sumw2();

  flist_TghCuts_Signal->Add(hCPtaVSd0d0_tgh_sign);
  flist_TghCuts_Signal->Add(hSecVtxZ_tgh_sign);
  flist_TghCuts_Signal->Add(hSecVtxY_tgh_sign);
  flist_TghCuts_Signal->Add(hSecVtxX_tgh_sign);
  flist_TghCuts_Signal->Add(hSecVtxXY_tgh_sign);
  flist_TghCuts_Signal->Add(hSecVtxPhi_tgh_sign);
  flist_TghCuts_Signal->Add(hCPta_tgh_sign);
  flist_TghCuts_Signal->Add(hd0xd0_tgh_sign);
  flist_TghCuts_Signal->Add(hMassTrue_tgh_sign);
  flist_TghCuts_Signal->Add(hMass_tgh_sign);
  flist_TghCuts_Signal->Add(hMassTrue_tgh_sign_pm);
  flist_TghCuts_Signal->Add(hMass_tgh_sign_pm);
  flist_TghCuts_Signal->Add(hMassTrue_SB_tgh_sign);
  flist_TghCuts_Signal->Add(hMass_SB_tgh_sign);

  // ####### d0 D0 histos ############
  TH1F *hd0D0_tgh_sign_pm = new TH1F("hd0D0_tgh_sign_pm","D^{0} impact par. plot , Tight Cuts ,Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_tgh_sign_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_tgh_sign_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_tgh_sign_pm = new TH1F("hd0D0VtxTrue_tgh_sign_pm","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_tgh_sign_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_tgh_sign_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_tgh_sign_pm = new TH1F("hMCd0D0_tgh_sign_pm","D^{0} impact par. plot, Tight Cuts, Signal,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_tgh_sign_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_tgh_sign_pm->SetYTitle("Entries");

  TH1F *hd0D0_tgh_sign_sb = new TH1F("hd0D0_tgh_sign_sb","D^{0} impact par. plot , Tight Cuts ,Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_tgh_sign_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_tgh_sign_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_tgh_sign_sb = new TH1F("hd0D0VtxTrue_tgh_sign_sb","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, Signal,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_tgh_sign_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_tgh_sign_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_tgh_sign_sb = new TH1F("hMCd0D0_tgh_sign_sb","D^{0} impact par. plot, Tight Cuts, Signal,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_tgh_sign_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_tgh_sign_sb->SetYTitle("Entries");

  flist_TghCuts_Signal->Add(hd0D0_tgh_sign_pm);
  flist_TghCuts_Signal->Add(hd0D0VtxTrue_tgh_sign_pm);
  flist_TghCuts_Signal->Add(hMCd0D0_tgh_sign_pm);
  flist_TghCuts_Signal->Add(hd0D0_tgh_sign_sb);
  flist_TghCuts_Signal->Add(hd0D0VtxTrue_tgh_sign_sb);
  flist_TghCuts_Signal->Add(hMCd0D0_tgh_sign_sb);
  
  TH1F **hd0D0pt_tgh_sign_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_tgh_sign_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_tgh_sign_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_tgh_sign_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_tgh_sign_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_tgh_sign_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_tgh_sign_";
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
    
    hd0D0pt_tgh_sign_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_tgh_sign_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_tgh_sign_pm[i]->SetYTitle("Entries");
    flist_TghCuts_Signal->Add(hd0D0pt_tgh_sign_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_tgh_sign_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_tgh_sign_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_tgh_sign_pm[i]->SetYTitle("Entries");
    flist_TghCuts_Signal->Add(hMCd0D0pt_tgh_sign_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_tgh_sign_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_tgh_sign_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_tgh_sign_pm[i]->SetYTitle("Entries");
    flist_TghCuts_Signal->Add(hd0D0VtxTruept_tgh_sign_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_tgh_sign_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_tgh_sign_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_tgh_sign_sb[i]->SetYTitle("Entries");
    flist_TghCuts_Signal->Add(hd0D0pt_tgh_sign_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_tgh_sign_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_tgh_sign_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_tgh_sign_sb[i]->SetYTitle("Entries");
    flist_TghCuts_Signal->Add(hMCd0D0pt_tgh_sign_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_tgh_sign_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_tgh_sign_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_tgh_sign_sb[i]->SetYTitle("Entries");
    flist_TghCuts_Signal->Add(hd0D0VtxTruept_tgh_sign_sb[i]);
  }


  //############ TIGHT CUTS BACKGROUND HISTOGRAMS ###########
  //
  //   ######## global properties histos #######
  TH2F *hCPtaVSd0d0_tgh_back=new TH2F("hCPtaVSd0d0_tgh_back","hCPtaVSd0d0_TightCuts_Background",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_tgh_back=new TH1F("hSecVtxZ_tgh_back","hSecVtxZ_TightCuts_Background",1000,-8.,8.);
  TH1F *hSecVtxX_tgh_back=new TH1F("hSecVtxX_tgh_back","hSecVtxX_TightCuts_Background",1000,-3000.,3000.);
  TH1F *hSecVtxY_tgh_back=new TH1F("hSecVtxY_tgh_back","hSecVtxY_TightCuts_Background",1000,-3000.,3000.);
  TH2F *hSecVtxXY_tgh_back=new TH2F("hSecVtxXY_tgh_back","hSecVtxXY_TightCuts_Background",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_tgh_back=new TH1F("hSecVtxPhi_tgh_back","hSecVtxPhi_TightCuts_Background",180,-180.1,180.1);
  TH1F *hCPta_tgh_back=new TH1F("hCPta_tgh_back","hCPta_TightCuts_Background",100,0.,1.);
  TH1F *hd0xd0_tgh_back=new TH1F("hd0xd0_tgh_back","hd0xd0_TightCuts_Background",1000,-100000.,100000.);
  TH1F *hMassTrue_tgh_back=new TH1F("hMassTrue_tgh_back","D^{0} MC inv. Mass Tight Cuts Background(All momenta)",600,1.600,2.200);
  TH1F *hMass_tgh_back=new TH1F("hMass_tgh_back","D^{0} inv. Mass Tight Cuts Background (All momenta)",600,1.600,2.200);
  hMass_tgh_back->Sumw2();
  TH1F *hMassTrue_tgh_back_pm=new TH1F("hMassTrue_tgh_back_pm","D^{0} MC inv. Mass Tight Cuts Background, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_tgh_back_pm=new TH1F("hMass_tgh_back_pm","D^{0} inv. Mass Tight Cuts Background (All momenta), MassPeak",600,1.600,2.200);
  hMass_tgh_back_pm->Sumw2();
  TH1F *hMassTrue_SB_tgh_back=new TH1F("hMassTrue_tgh_back_sb","D^{0} MC inv. Mass in Side Bands Tight Cuts Backgrround(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_tgh_back=new TH1F("hMass_tgh_back_sb","D^{0} inv. Mass in Side Bands Tight Cuts Background (All momenta)",600,1.600,2.200);
  hMass_SB_tgh_back->Sumw2();

  flist_TghCuts_Back->Add(hCPtaVSd0d0_tgh_back);
  flist_TghCuts_Back->Add(hSecVtxZ_tgh_back);
  flist_TghCuts_Back->Add(hSecVtxY_tgh_back);
  flist_TghCuts_Back->Add(hSecVtxX_tgh_back);
  flist_TghCuts_Back->Add(hSecVtxXY_tgh_back);
  flist_TghCuts_Back->Add(hSecVtxPhi_tgh_back);
  flist_TghCuts_Back->Add(hCPta_tgh_back);
  flist_TghCuts_Back->Add(hd0xd0_tgh_back);
  flist_TghCuts_Back->Add(hMassTrue_tgh_back);
  flist_TghCuts_Back->Add(hMass_tgh_back);
  flist_TghCuts_Back->Add(hMassTrue_tgh_back_pm);
  flist_TghCuts_Back->Add(hMass_tgh_back_pm);
  flist_TghCuts_Back->Add(hMassTrue_SB_tgh_back);
  flist_TghCuts_Back->Add(hMass_SB_tgh_back);


  // ####### d0 D0 histos ############
  
 TH1F *hd0D0_tgh_back_pm = new TH1F("hd0D0_tgh_back_pm","D^{0} impact par. plot , Tight Cuts ,Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_tgh_back_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_tgh_back_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_tgh_back_pm = new TH1F("hd0D0VtxTrue_tgh_back_pm","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_tgh_back_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_tgh_back_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_tgh_back_pm = new TH1F("hMCd0D0_tgh_back_pm","D^{0} impact par. plot, Tight Cuts, Background,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_tgh_back_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_tgh_back_pm->SetYTitle("Entries");

  TH1F *hd0D0_tgh_back_sb = new TH1F("hd0D0_tgh_back_sb","D^{0} impact par. plot , Tight Cuts ,Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_tgh_back_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_tgh_back_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_tgh_back_sb = new TH1F("hd0D0VtxTrue_tgh_back_sb","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, Background,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_tgh_back_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_tgh_back_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_tgh_back_sb = new TH1F("hMCd0D0_tgh_back_sb","D^{0} impact par. plot, Tight Cuts, Background,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_tgh_back_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_tgh_back_sb->SetYTitle("Entries");

  flist_TghCuts_Back->Add(hd0D0_tgh_back_pm);
  flist_TghCuts_Back->Add(hd0D0VtxTrue_tgh_back_pm);
  flist_TghCuts_Back->Add(hMCd0D0_tgh_back_pm);
  flist_TghCuts_Back->Add(hd0D0_tgh_back_sb);
  flist_TghCuts_Back->Add(hd0D0VtxTrue_tgh_back_sb);
  flist_TghCuts_Back->Add(hMCd0D0_tgh_back_sb);
  
  TH1F **hd0D0pt_tgh_back_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_tgh_back_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_tgh_back_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_tgh_back_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_tgh_back_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_tgh_back_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_tgh_back_";
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
    
    hd0D0pt_tgh_back_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_tgh_back_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_tgh_back_pm[i]->SetYTitle("Entries");
    flist_TghCuts_Back->Add(hd0D0pt_tgh_back_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_tgh_back_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_tgh_back_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_tgh_back_pm[i]->SetYTitle("Entries");
    flist_TghCuts_Back->Add(hMCd0D0pt_tgh_back_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_tgh_back_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_tgh_back_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_tgh_back_pm[i]->SetYTitle("Entries");
    flist_TghCuts_Back->Add(hd0D0VtxTruept_tgh_back_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_tgh_back_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_tgh_back_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_tgh_back_sb[i]->SetYTitle("Entries");
    flist_TghCuts_Back->Add(hd0D0pt_tgh_back_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_tgh_back_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_tgh_back_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_tgh_back_sb[i]->SetYTitle("Entries");
    flist_TghCuts_Back->Add(hMCd0D0pt_tgh_back_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_tgh_back_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_tgh_back_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_tgh_back_sb[i]->SetYTitle("Entries");
    flist_TghCuts_Back->Add(hd0D0VtxTruept_tgh_back_sb[i]);
  }



 //############ TIGHT CUTS FROMB HISTOGRAMS ###########
  //
  //#######  global properties histos

  TH2F *hCPtaVSd0d0_tgh_fromB=new TH2F("hCPtaVSd0d0_tgh_fromB","hCPtaVSd0d0_TightCuts_FromB",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_tgh_fromB=new TH1F("hSecVtxZ_tgh_fromB","hSecVtxZ_TightCuts_FromB",1000,-8.,8.);
  TH1F *hSecVtxX_tgh_fromB=new TH1F("hSecVtxX_tgh_fromB","hSecVtxX_TightCuts_FromB",1000,-3000.,3000.);
  TH1F *hSecVtxY_tgh_fromB=new TH1F("hSecVtxY_tgh_fromB","hSecVtxY_TightCuts_FromB",1000,-3000.,3000.);
  TH2F *hSecVtxXY_tgh_fromB=new TH2F("hSecVtxXY_tgh_fromB","hSecVtxXY_TightCuts_FromB",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_tgh_fromB=new TH1F("hSecVtxPhi_tgh_fromB","hSecVtxPhi_TightCuts_FromB",180,-180.1,180.1);
  TH1F *hCPta_tgh_fromB=new TH1F("hCPta_tgh_fromB","hCPta_TightCuts_FromB",100,0.,1.);
  TH1F *hd0xd0_tgh_fromB=new TH1F("hd0xd0_tgh_fromB","hd0xd0_TightCuts_FromB",1000,-100000.,100000.);
  TH1F *hMassTrue_tgh_fromB=new TH1F("hMassTrue_tgh_fromB","D^{0} MC inv. Mass Tight Cuts FromB(All momenta)",600,1.600,2.200);
  TH1F *hMass_tgh_fromB=new TH1F("hMass_tgh_fromB","D^{0} inv. Mass Tight Cuts FromB (All momenta)",600,1.600,2.200);
  hMass_tgh_fromB->Sumw2();
  TH1F *hMassTrue_tgh_fromB_pm=new TH1F("hMassTrue_tgh_fromB_pm","D^{0} MC inv. Mass Tight Cuts FromB, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_tgh_fromB_pm=new TH1F("hMass_tgh_fromB_pm","D^{0} inv. Mass Tight Cuts FromB (All momenta), MassPeak",600,1.600,2.200);
  hMass_tgh_fromB_pm->Sumw2();
  TH1F *hMassTrue_SB_tgh_fromB=new TH1F("hMassTrue_tgh_fromB_sb","D^{0} MC inv. Mass in Side Bands Tight Cuts FromB(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_tgh_fromB=new TH1F("hMass_tgh_fromB_sb","D^{0} inv. Mass in Side Bands Tight Cuts FromB (All momenta)",600,1.600,2.200);
  hMass_SB_tgh_fromB->Sumw2();

  flist_TghCuts_FromB->Add(hCPtaVSd0d0_tgh_fromB);
  flist_TghCuts_FromB->Add(hSecVtxZ_tgh_fromB);
  flist_TghCuts_FromB->Add(hSecVtxY_tgh_fromB);
  flist_TghCuts_FromB->Add(hSecVtxX_tgh_fromB);
  flist_TghCuts_FromB->Add(hSecVtxXY_tgh_fromB);
  flist_TghCuts_FromB->Add(hSecVtxPhi_tgh_fromB);
  flist_TghCuts_FromB->Add(hCPta_tgh_fromB);
  flist_TghCuts_FromB->Add(hd0xd0_tgh_fromB);
  flist_TghCuts_FromB->Add(hMassTrue_tgh_fromB);
  flist_TghCuts_FromB->Add(hMass_tgh_fromB);
  flist_TghCuts_FromB->Add(hMassTrue_tgh_fromB_pm);
  flist_TghCuts_FromB->Add(hMass_tgh_fromB_pm);
  flist_TghCuts_FromB->Add(hMassTrue_SB_tgh_fromB);
  flist_TghCuts_FromB->Add(hMass_SB_tgh_fromB);

  // ######### d0 D0 histos ##############
  TH1F *hd0D0_tgh_fromB_pm = new TH1F("hd0D0_tgh_fromB_pm","D^{0} impact par. plot , Tight Cuts ,FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_tgh_fromB_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_tgh_fromB_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_tgh_fromB_pm = new TH1F("hd0D0VtxTrue_tgh_fromB_pm","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_tgh_fromB_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_tgh_fromB_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_tgh_fromB_pm = new TH1F("hMCd0D0_tgh_fromB_pm","D^{0} impact par. plot, Tight Cuts, FromB,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_tgh_fromB_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_tgh_fromB_pm->SetYTitle("Entries");

  TH1F *hd0D0_tgh_fromB_sb = new TH1F("hd0D0_tgh_fromB_sb","D^{0} impact par. plot , Tight Cuts ,FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_tgh_fromB_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_tgh_fromB_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_tgh_fromB_sb = new TH1F("hd0D0VtxTrue_tgh_fromB_sb","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, FromB,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_tgh_fromB_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_tgh_fromB_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_tgh_fromB_sb = new TH1F("hMCd0D0_tgh_fromB_sb","D^{0} impact par. plot, Tight Cuts, FromB,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_tgh_fromB_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_tgh_fromB_sb->SetYTitle("Entries");

  flist_TghCuts_FromB->Add(hd0D0_tgh_fromB_pm);
  flist_TghCuts_FromB->Add(hd0D0VtxTrue_tgh_fromB_pm);
  flist_TghCuts_FromB->Add(hMCd0D0_tgh_fromB_pm);
  flist_TghCuts_FromB->Add(hd0D0_tgh_fromB_sb);
  flist_TghCuts_FromB->Add(hd0D0VtxTrue_tgh_fromB_sb);
  flist_TghCuts_FromB->Add(hMCd0D0_tgh_fromB_sb);
  
  TH1F **hd0D0pt_tgh_fromB_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_tgh_fromB_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_tgh_fromB_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_tgh_fromB_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_tgh_fromB_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_tgh_fromB_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_tgh_fromB_";
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
    
    hd0D0pt_tgh_fromB_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_tgh_fromB_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_tgh_fromB_pm[i]->SetYTitle("Entries");
    flist_TghCuts_FromB->Add(hd0D0pt_tgh_fromB_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_tgh_fromB_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_tgh_fromB_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_tgh_fromB_pm[i]->SetYTitle("Entries");
    flist_TghCuts_FromB->Add(hMCd0D0pt_tgh_fromB_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_tgh_fromB_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_tgh_fromB_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_tgh_fromB_pm[i]->SetYTitle("Entries");
    flist_TghCuts_FromB->Add(hd0D0VtxTruept_tgh_fromB_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_tgh_fromB_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_tgh_fromB_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_tgh_fromB_sb[i]->SetYTitle("Entries");
    flist_TghCuts_FromB->Add(hd0D0pt_tgh_fromB_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_tgh_fromB_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_tgh_fromB_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_tgh_fromB_sb[i]->SetYTitle("Entries");
    flist_TghCuts_FromB->Add(hMCd0D0pt_tgh_fromB_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_tgh_fromB_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_tgh_fromB_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_tgh_fromB_sb[i]->SetYTitle("Entries");
    flist_TghCuts_FromB->Add(hd0D0VtxTruept_tgh_fromB_sb[i]);
  }



 //############ TIGHT CUTS FROM DSTAR HISTOGRAMS ###########
 //
  //############## global properties histos
  TH2F *hCPtaVSd0d0_tgh_fromDstar=new TH2F("hCPtaVSd0d0_tgh_fromDstar","hCPtaVSd0d0_TightCuts_FromDStar",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_tgh_fromDstar=new TH1F("hSecVtxZ_tgh_fromDstar","hSecVtxZ_TightCuts_FromDStar",1000,-8.,8.);
  TH1F *hSecVtxX_tgh_fromDstar=new TH1F("hSecVtxX_tgh_fromDstar","hSecVtxX_TightCuts_FromDStar",1000,-3000.,3000.);
  TH1F *hSecVtxY_tgh_fromDstar=new TH1F("hSecVtxY_tgh_fromDstar","hSecVtxY_TightCuts_FromDStar",1000,-3000.,3000.);
  TH2F *hSecVtxXY_tgh_fromDstar=new TH2F("hSecVtxXY_tgh_fromDstar","hSecVtxXY_TightCuts_FromDStar",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_tgh_fromDstar=new TH1F("hSecVtxPhi_tgh_fromDstar","hSecVtxPhi_TightCuts_FromDStar",180,-180.1,180.1);
  TH1F *hCPta_tgh_fromDstar=new TH1F("hCPta_tgh_fromDstar","hCPta_TightCuts_FromDStar",100,0.,1.);
  TH1F *hd0xd0_tgh_fromDstar=new TH1F("hd0xd0_tgh_fromDstar","hd0xd0_TightCuts_FromDStar",1000,-100000.,100000.);
  TH1F *hMassTrue_tgh_fromDstar=new TH1F("hMassTrue_tgh_fromDstar","D^{0} MC inv. Mass Tight Cuts FromDStar(All momenta)",600,1.600,2.200);
  TH1F *hMass_tgh_fromDstar=new TH1F("hMass_tgh_fromDstar","D^{0} inv. Mass Tight Cuts FromDStar (All momenta)",600,1.600,2.200);
  hMass_tgh_fromDstar->Sumw2();
  TH1F *hMassTrue_tgh_fromDstar_pm=new TH1F("hMassTrue_tgh_fromDstar_pm","D^{0} MC inv. Mass Tight Cuts FromDStar, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_tgh_fromDstar_pm=new TH1F("hMass_tgh_fromDstar_pm","D^{0} inv. Mass Tight Cuts FromDStar (All momenta), MassPeak",600,1.600,2.200);
  hMass_tgh_fromDstar_pm->Sumw2();
  TH1F *hMassTrue_SB_tgh_fromDstar=new TH1F("hMassTrue_tgh_fromDstar_sb","D^{0} MC inv. Mass in Side Bands Tight Cuts FromDStar(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_tgh_fromDstar=new TH1F("hMass_tgh_fromDstar_sb","D^{0} inv. Mass in Side Bands Tight Cuts FromDStar (All momenta)",600,1.600,2.200);
  hMass_SB_tgh_fromDstar->Sumw2();

  flist_TghCuts_FromDstar->Add(hCPtaVSd0d0_tgh_fromDstar);
  flist_TghCuts_FromDstar->Add(hSecVtxZ_tgh_fromDstar);
  flist_TghCuts_FromDstar->Add(hSecVtxY_tgh_fromDstar);
  flist_TghCuts_FromDstar->Add(hSecVtxX_tgh_fromDstar);
  flist_TghCuts_FromDstar->Add(hSecVtxXY_tgh_fromDstar);
  flist_TghCuts_FromDstar->Add(hSecVtxPhi_tgh_fromDstar);
  flist_TghCuts_FromDstar->Add(hCPta_tgh_fromDstar);
  flist_TghCuts_FromDstar->Add(hd0xd0_tgh_fromDstar);
  flist_TghCuts_FromDstar->Add(hMassTrue_tgh_fromDstar);
  flist_TghCuts_FromDstar->Add(hMass_tgh_fromDstar);
  flist_TghCuts_FromDstar->Add(hMassTrue_tgh_fromDstar_pm);
  flist_TghCuts_FromDstar->Add(hMass_tgh_fromDstar_pm);
  flist_TghCuts_FromDstar->Add(hMassTrue_SB_tgh_fromDstar);
  flist_TghCuts_FromDstar->Add(hMass_SB_tgh_fromDstar);

  //########## d0 D0 histos #############  
  TH1F *hd0D0_tgh_fromDst_pm = new TH1F("hd0D0_tgh_fromDstar_pm","D^{0} impact par. plot , Tight Cuts ,FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_tgh_fromDst_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_tgh_fromDst_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_tgh_fromDst_pm = new TH1F("hd0D0VtxTrue_tgh_fromDstar_pm","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_tgh_fromDst_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_tgh_fromDst_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_tgh_fromDst_pm = new TH1F("hMCd0D0_tgh_fromDstar_pm","D^{0} impact par. plot, Tight Cuts, FromDStar,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_tgh_fromDst_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_tgh_fromDst_pm->SetYTitle("Entries");

  TH1F *hd0D0_tgh_fromDst_sb = new TH1F("hd0D0_tgh_fromDstar_sb","D^{0} impact par. plot , Tight Cuts ,FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_tgh_fromDst_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_tgh_fromDst_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_tgh_fromDst_sb = new TH1F("hd0D0VtxTrue_tgh_fromDstar_sb","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, FromDStar,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_tgh_fromDst_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_tgh_fromDst_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_tgh_fromDst_sb = new TH1F("hMCd0D0_tgh_fromDstar_sb","D^{0} impact par. plot, Tight Cuts, FromDStar,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_tgh_fromDst_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_tgh_fromDst_sb->SetYTitle("Entries");

  flist_TghCuts_FromDstar->Add(hd0D0_tgh_fromDst_pm);
  flist_TghCuts_FromDstar->Add(hd0D0VtxTrue_tgh_fromDst_pm);
  flist_TghCuts_FromDstar->Add(hMCd0D0_tgh_fromDst_pm);
  flist_TghCuts_FromDstar->Add(hd0D0_tgh_fromDst_sb);
  flist_TghCuts_FromDstar->Add(hd0D0VtxTrue_tgh_fromDst_sb);
  flist_TghCuts_FromDstar->Add(hMCd0D0_tgh_fromDst_sb);
  
  TH1F **hd0D0pt_tgh_fromDst_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_tgh_fromDst_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_tgh_fromDst_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_tgh_fromDst_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_tgh_fromDst_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_tgh_fromDst_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_tgh_fromDstar_";
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
    
    hd0D0pt_tgh_fromDst_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_tgh_fromDst_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_tgh_fromDst_pm[i]->SetYTitle("Entries");
    flist_TghCuts_FromDstar->Add(hd0D0pt_tgh_fromDst_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_tgh_fromDst_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_tgh_fromDst_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_tgh_fromDst_pm[i]->SetYTitle("Entries");
    flist_TghCuts_FromDstar->Add(hMCd0D0pt_tgh_fromDst_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_tgh_fromDst_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_tgh_fromDst_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_tgh_fromDst_pm[i]->SetYTitle("Entries");
    flist_TghCuts_FromDstar->Add(hd0D0VtxTruept_tgh_fromDst_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_tgh_fromDst_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_tgh_fromDst_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_tgh_fromDst_sb[i]->SetYTitle("Entries");
    flist_TghCuts_FromDstar->Add(hd0D0pt_tgh_fromDst_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_tgh_fromDst_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_tgh_fromDst_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_tgh_fromDst_sb[i]->SetYTitle("Entries");
    flist_TghCuts_FromDstar->Add(hMCd0D0pt_tgh_fromDst_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_tgh_fromDst_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_tgh_fromDst_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_tgh_fromDst_sb[i]->SetYTitle("Entries");
    flist_TghCuts_FromDstar->Add(hd0D0VtxTruept_tgh_fromDst_sb[i]);
  }


  //############ TIGHT CUTS OTHER HISTOGRAMS ###########
  //
  //########### global properties histos ###########

  TH2F *hCPtaVSd0d0_tgh_other=new TH2F("hCPtaVSd0d0_tgh_other","hCPtaVSd0d0_TightCuts_other",1000,-100000.,100000.,100,0.,1.);
  TH1F *hSecVtxZ_tgh_other=new TH1F("hSecVtxZ_tgh_other","hSecVtxZ_TightCuts_other",1000,-8.,8.);
  TH1F *hSecVtxX_tgh_other=new TH1F("hSecVtxX_tgh_other","hSecVtxX_TightCuts_other",1000,-3000.,3000.);
  TH1F *hSecVtxY_tgh_other=new TH1F("hSecVtxY_tgh_other","hSecVtxY_TightCuts_other",1000,-3000.,3000.);
  TH2F *hSecVtxXY_tgh_other=new TH2F("hSecVtxXY_tgh_other","hSecVtxXY_TightCuts_other",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhi_tgh_other=new TH1F("hSecVtxPhi_tgh_other","hSecVtxPhi_TightCuts_other",180,-180.1,180.1);
  TH1F *hCPta_tgh_other=new TH1F("hCPta_tgh_other","hCPta_TightCuts_other",100,0.,1.);
  TH1F *hd0xd0_tgh_other=new TH1F("hd0xd0_tgh_other","hd0xd0_TightCuts_other",1000,-100000.,100000.);
  TH1F *hMassTrue_tgh_other=new TH1F("hMassTrue_tgh_other","D^{0} MC inv. Mass Tight Cuts other(All momenta)",600,1.600,2.200);
  TH1F *hMass_tgh_other=new TH1F("hMass_tgh_other","D^{0} inv. Mass Tight Cuts other (All momenta)",600,1.600,2.200);
  hMass_tgh_other->Sumw2();
  TH1F *hMassTrue_tgh_other_pm=new TH1F("hMassTrue_tgh_other_pm","D^{0} MC inv. Mass Tight Cuts other, Mass Peak. (All momenta)",600,1.600,2.200);
  TH1F *hMass_tgh_other_pm=new TH1F("hMass_tgh_other_pm","D^{0} inv. Mass Tight Cuts other (All momenta), MassPeak",600,1.600,2.200);
  hMass_tgh_other_pm->Sumw2();
  TH1F *hMassTrue_SB_tgh_other=new TH1F("hMassTrue_tgh_other_sb","D^{0} MC inv. Mass in Side Bands Tight Cuts other(All momenta)",600,1.600,2.200);
  TH1F *hMass_SB_tgh_other=new TH1F("hMass_tgh_other_sb","D^{0} inv. Mass in Side Bands Tight Cuts other (All momenta)",600,1.600,2.200);
  hMass_SB_tgh_other->Sumw2();

  flist_TghCuts_Other->Add(hCPtaVSd0d0_tgh_other);
  flist_TghCuts_Other->Add(hSecVtxZ_tgh_other);
  flist_TghCuts_Other->Add(hSecVtxY_tgh_other);
  flist_TghCuts_Other->Add(hSecVtxX_tgh_other);
  flist_TghCuts_Other->Add(hSecVtxXY_tgh_other);
  flist_TghCuts_Other->Add(hSecVtxPhi_tgh_other);
  flist_TghCuts_Other->Add(hCPta_tgh_other);
  flist_TghCuts_Other->Add(hd0xd0_tgh_other);
  flist_TghCuts_Other->Add(hMassTrue_tgh_other);
  flist_TghCuts_Other->Add(hMass_tgh_other);
  flist_TghCuts_Other->Add(hMassTrue_tgh_other_pm);
  flist_TghCuts_Other->Add(hMass_tgh_other_pm);
  flist_TghCuts_Other->Add(hMassTrue_SB_tgh_other);
  flist_TghCuts_Other->Add(hMass_SB_tgh_other);

  //############# d0 D0 histos ###############à
  TH1F *hd0D0_tgh_other_pm = new TH1F("hd0D0_tgh_other_pm","D^{0} impact par. plot , Tight Cuts ,Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_tgh_other_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0_tgh_other_pm->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_tgh_other_pm = new TH1F("hd0D0VtxTrue_tgh_other_pm","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_tgh_other_pm->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_tgh_other_pm->SetYTitle("Entries");

  TH1F *hMCd0D0_tgh_other_pm = new TH1F("hMCd0D0_tgh_other_pm","D^{0} impact par. plot, Tight Cuts, Other,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_tgh_other_pm->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_tgh_other_pm->SetYTitle("Entries");

  TH1F *hd0D0_tgh_other_sb = new TH1F("hd0D0_tgh_other_sb","D^{0} impact par. plot , Tight Cuts ,Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0_tgh_other_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0_tgh_other_sb->SetYTitle("Entries");

  TH1F *hd0D0VtxTrue_tgh_other_sb = new TH1F("hd0D0VtxTrue_tgh_other_sb","D^{0} impact par. w.r.t. True Vtx, Tight Cuts, Other,Mass Peak (All momenta)",1000,-1000.,1000.);
  hd0D0VtxTrue_tgh_other_sb->SetXTitle("Impact parameter [#mum]");
  hd0D0VtxTrue_tgh_other_sb->SetYTitle("Entries");

  TH1F *hMCd0D0_tgh_other_sb = new TH1F("hMCd0D0_tgh_other_sb","D^{0} impact par. plot, Tight Cuts, Other,Mass Peak  (All momenta)",1000,-1000.,1000.);
  hMCd0D0_tgh_other_sb->SetXTitle("MC Impact parameter [#mum]");
  hMCd0D0_tgh_other_sb->SetYTitle("Entries");

  flist_TghCuts_Other->Add(hd0D0_tgh_other_pm);
  flist_TghCuts_Other->Add(hd0D0VtxTrue_tgh_other_pm);
  flist_TghCuts_Other->Add(hMCd0D0_tgh_other_pm);
  flist_TghCuts_Other->Add(hd0D0_tgh_other_sb);
  flist_TghCuts_Other->Add(hd0D0VtxTrue_tgh_other_sb);
  flist_TghCuts_Other->Add(hMCd0D0_tgh_other_sb);
  
  TH1F **hd0D0pt_tgh_other_pm=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_tgh_other_pm=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_tgh_other_pm=new TH1F*[fnbins];
  TH1F **hd0D0pt_tgh_other_sb=new TH1F*[fnbins];
  TH1F **hMCd0D0pt_tgh_other_sb=new TH1F*[fnbins];
  TH1F ** hd0D0VtxTruept_tgh_other_sb=new TH1F*[fnbins];
  namehist="hd0D0pt_tgh_other_";
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
    
    hd0D0pt_tgh_other_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_tgh_other_pm[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_tgh_other_pm[i]->SetYTitle("Entries");
    flist_TghCuts_Other->Add(hd0D0pt_tgh_other_pm[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_tgh_other_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_tgh_other_pm[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_tgh_other_pm[i]->SetYTitle("Entries");
    flist_TghCuts_Other->Add(hMCd0D0pt_tgh_other_pm[i]);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_tgh_other_pm[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_tgh_other_pm[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_tgh_other_pm[i]->SetYTitle("Entries");
    flist_TghCuts_Other->Add(hd0D0VtxTruept_tgh_other_pm[i]);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0pt_tgh_other_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0pt_tgh_other_sb[i]->SetXTitle("Impact parameter [#mum] ");
    hd0D0pt_tgh_other_sb[i]->SetYTitle("Entries");
    flist_TghCuts_Other->Add(hd0D0pt_tgh_other_sb[i]);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0pt_tgh_other_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0pt_tgh_other_sb[i]->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0pt_tgh_other_sb[i]->SetYTitle("Entries");
    flist_TghCuts_Other->Add(hMCd0D0pt_tgh_other_sb[i]);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTruept_tgh_other_sb[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTruept_tgh_other_sb[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTruept_tgh_other_sb[i]->SetYTitle("Entries");
    flist_TghCuts_Other->Add(hd0D0VtxTruept_tgh_other_sb[i]);
  }

  
}



//________________________________________________________________________
void AliAnalysisTaskSECharmFraction::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

 
  // load D0->Kpi candidates                                                   
  TClonesArray *arrayD0toKpi =
    (TClonesArray*)aod->GetList()->FindObject("D0toKpi");
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
  Int_t okD0_tight,okD0bar_tight,okD0_loose,okD0bar_loose;
  Bool_t isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,isSideBand;
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
    okD0_tight=0;
    okD0bar_tight=0;
    okD0_loose=0;
    okD0bar_loose=0;
  
    signallevel=-1;
    

    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->UncheckedAt(iD0toKpi);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
    
    ptD0=d->Pt();
      
    //####### DATA SELECTION ###########
    
    
    //######## INVARIANT MASS SELECTION ###############
    CheckInvMassD0(d,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar);
    if(isSideBandD0&&isSideBandD0bar)isSideBand=kTRUE;// TEMPORARY, NOT DONE IN THE METHOD CALLED ABOVE ONLY FOR FURTHER SIDE BAND STUDY
  
    // INVESTIGATE SIGNAL TYPE : ACCESS TO MC INFORMATION
    aodDMC=GetD0toKPiSignalType(d,arrayMC,signallevel,massmumtrue,vtxTrue);
    fSignalType->Fill(signallevel);
    // END OF BACKGROUND TYPE SELECTION

    // NOW APPLY CUTS
    //NO CUTS CASE IS FOR FREE
    
    // CHECK TIGHTER CUTS 
    ptbin=SetStandardCuts(ptD0,flargeInvMassCut);
    d->SelectD0(fVHFtight->GetD0toKpiCuts(),okD0_tight,okD0bar_tight);
    if((isPeakD0&&okD0_tight)||(isPeakD0bar&&okD0bar_tight))fSignalTypeTghCuts->Fill(signallevel);
    
    // CHECK LOOSER CUTS 
    ptbin=SetStandardCuts(ptD0,flargeInvMassCut);
    d->SelectD0(fVHFloose->GetD0toKpiCuts(),okD0_loose,okD0bar_loose);
    if((isPeakD0&&okD0_loose)||(isPeakD0bar&&okD0bar_loose))fSignalTypeLsCuts->Fill(signallevel);
   
    
    //###################    FILL HISTOS      ########################
    //################################################################
    //
    //######## improvement: SPEED HERE CAN BE IMPROVED: CALCULATE ONCE AND FOR ALL 
    //            CANDIDATE VARIABLES   

    //NO CUTS Case: force okD0 and okD0bar = kTRUE
    if(signallevel==1)FillHistos(d,flist_NoCuts_Signal,ptbin,kTRUE,kTRUE,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==2)FillHistos(d,flist_NoCuts_FromDstar,ptbin,kTRUE,kTRUE,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==3||signallevel==4)FillHistos(d,flist_NoCuts_FromB,ptbin,kTRUE,kTRUE,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==-1||signallevel==7||signallevel==8)FillHistos(d,flist_NoCuts_Back,ptbin,kTRUE,kTRUE,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==5||signallevel==6)FillHistos(d,flist_NoCuts_Other,ptbin,kTRUE,kTRUE,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);

    //LOOSE CUTS Case
    if(signallevel==1)FillHistos(d,flist_LsCuts_Signal,ptbin,okD0_loose,okD0bar_loose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==2)FillHistos(d,flist_LsCuts_FromDstar,ptbin,okD0_loose,okD0bar_loose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==3||signallevel==4)FillHistos(d,flist_LsCuts_FromB,ptbin,okD0_loose,okD0bar_loose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==-1||signallevel==7||signallevel==8)FillHistos(d,flist_LsCuts_Back,ptbin,okD0_loose,okD0bar_loose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==5||signallevel==6)FillHistos(d,flist_LsCuts_Other,ptbin,okD0_loose,okD0bar_loose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);

    //TIGHT CUTS Case
    if(signallevel==1)FillHistos(d,flist_TghCuts_Signal,ptbin,okD0_tight,okD0bar_tight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==2)FillHistos(d,flist_TghCuts_FromDstar,ptbin,okD0_tight,okD0bar_tight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==3||signallevel==4)FillHistos(d,flist_TghCuts_FromB,ptbin,okD0_tight,okD0bar_tight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==-1||signallevel==7||signallevel==8)FillHistos(d,flist_TghCuts_Back,ptbin,okD0_tight,okD0bar_tight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==5||signallevel==6)FillHistos(d,flist_TghCuts_Other,ptbin,okD0_tight,okD0bar_tight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBand,massmumtrue,aodDMC,vtxTrue);
   
    
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
  PostData(5,flist_NoCuts_Signal);
  PostData(6,flist_NoCuts_Back);
  PostData(7,flist_NoCuts_FromB);
  PostData(8,flist_NoCuts_FromDstar);
  PostData(9,flist_NoCuts_Other);
  PostData(10,flist_LsCuts_Signal);
  PostData(11,flist_LsCuts_Back);
  PostData(12,flist_LsCuts_FromB);
  PostData(13,flist_LsCuts_FromDstar);
  PostData(14,flist_LsCuts_Other);
  PostData(15,flist_TghCuts_Signal);
  PostData(16,flist_TghCuts_Back);
  PostData(17,flist_TghCuts_FromB);
  PostData(18,flist_TghCuts_FromDstar);
  PostData(19,flist_TghCuts_Other);

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
  Int_t ptbin=-1;
  if (pt>0. && pt<=1.) {
    ptbin=0;
    fVHFtight->SetD0toKpiCuts(invMassCut,0.04,0.8,0.5,0.5,0.05,0.05,-0.0003,0.7);
    fVHFloose->SetD0toKpiCuts(invMassCut,0.04,0.8,0.5,0.5,0.05,0.05,-0.00025,0.7);
  }
  
  if(pt>1. && pt<=3.) {
    ptbin=1;  
    fVHFtight->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.0003,0.9);
    fVHFloose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,1,1,-0.00025,0.8);
    //printf("I'm in the bin %d\n",ptbin);
  }
  
  if(pt>3. && pt<=5.){
    ptbin=2;  
    fVHFtight->SetD0toKpiCuts(invMassCut,0.015,0.8,0.7,0.7,0.05,0.05,-0.0002,0.9);
    fVHFloose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.8);
    //printf("I'm in the bin %d\n",ptbin);
  }
  if(pt>5.){
    ptbin=3;
    fVHFtight->SetD0toKpiCuts(invMassCut,0.015,0.8,0.7,0.7,0.05,0.05,-0.0002,0.95);
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
AliAODRecoDecayHF* AliAnalysisTaskSECharmFraction::GetD0toKPiSignalType(AliAODRecoDecayHF2Prong *d,TClonesArray *arrayMC,Int_t &signaltype,Double_t &massMumTrue,Double_t *primaryVtx){
  //THIS METHOD CHECK THE TYPE OF SIGNAL/BACKGROUND THE CANDIDATE IS. 
  //  IF (!AND ONLY IF) THE TWO DAUGHTERS COME FROM A COMMONE MOTHER A FAKE TRUE SECONDARY VERTEX IS CONSTRUCTED (aodDMC)  
  //
  // THE FOLLOWING SCHEME IS ADOPTED: signaltype is set to
                        //  1:signal (D0 prompt); 2: signal D0 from Dstar; 3: D0 fromB 4: D0 from Dstar fromB
                        // then background categories: -1: one or both daughters is a fake track
                        //                             5: one or both daughters come from a D meson != D0
                        //                             6: both daughters come from a D0->4prongs  
                        //                             7: both daughetrs are primaries
                        //                             8: generic background (can include one of the previous if desired)
                        //
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
AliAODRecoDecayHF* AliAnalysisTaskSECharmFraction::ConstructFakeTrueSecVtx(AliAODMCParticle *b1, AliAODMCParticle *b2, AliAODMCParticle *mum,Double_t *primaryVtxTrue){
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
Bool_t AliAnalysisTaskSECharmFraction::FillHistos(AliAODRecoDecayHF2Prong *d,TList *&list,Int_t ptbin,Int_t okD0,Int_t okD0bar,Double_t invMassD0,Double_t invMassD0bar,Bool_t isPeakD0,Bool_t isPeakD0bar,Bool_t isSideBand,Double_t massmumtrue,AliAODRecoDecayHF *aodDMC,Double_t *vtxTrue){
  
  if((!okD0)&&(!okD0bar))return kTRUE;
  
  // ######### Get Standard label for hist in tlist ###############
  TString namehist=list->GetName(),str;
  namehist.ReplaceAll("list","");

  //  ######### Global properties histos #################
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


   //  ######### Invariant mass histos #################
   str="hMass";
   str.Append(namehist.Data());
   ((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0);
   ((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0bar);
   
      
   if(isPeakD0||isPeakD0bar){
     str="hMass";
     str.Append(namehist.Data());
     str.Append("_pm");
     if(isPeakD0&&okD0)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0);
     if(isPeakD0bar&&okD0bar)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0bar);
   }
   if(isSideBand){
     str="hMass";
     str.Append(namehist.Data());
     str.Append("_sb");
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
       str.Append("_pm");
       ((TH1F*)list->FindObject(str.Data()))->Fill(massmumtrue);
     }
     if(isSideBand){
       str="hMassTrue";
       str.Append(namehist.Data());
       str.Append("_sb");
       ((TH1F*)list->FindObject(str.Data()))->Fill(massmumtrue);
     }
   }

   // ################ D0 Impact Parameter Histos #####################
   if(isPeakD0||isPeakD0bar){
     str="hd0D0";
     str.Append(namehist.Data());
     str.Append("_pm");
     ((TH1F*)list->FindObject(str.Data()))->Fill(d->ImpParXY()*10000.);

     str="hd0D0pt";
     str.Append(namehist.Data());
     str.Append("_PkMss_pt");
     str+=ptbin;
     ((TH1F*)list->FindObject(str.Data()))->Fill(d->ImpParXY()*10000.);
     
   
     if(vtxTrue){
       str="hd0D0VtxTrue";
       str.Append(namehist.Data());
       str.Append("_pm");
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
       str.Append("_pm");
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
     str.Append("_sb");
     ((TH1F*)list->FindObject(str.Data()))->Fill(d->ImpParXY()*10000.);
 
     str="hd0D0pt";
     str.Append(namehist.Data());
     str.Append("_SBMss_pt");
     str+=ptbin;
     ((TH1F*)list->FindObject(str.Data()))->Fill(d->ImpParXY()*10000.);
     
  
     if(vtxTrue){
       str="hd0D0VtxTrue";
       str.Append(namehist.Data());
       str.Append("_sb");
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
       str.Append("_sb");
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

void AliAnalysisTaskSECharmFraction::Terminate(const Option_t*){
  //TERMINATE METHOD: NOTHING TO DO


}
