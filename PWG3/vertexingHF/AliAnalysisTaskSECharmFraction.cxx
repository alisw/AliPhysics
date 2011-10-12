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
// Author: Andrea Rossi, andrea.rossi@pd.infn.it
/////////////////////////////////////////////////////////////


#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TROOT.h>
#include "AliAODEvent.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecay.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliESDtrack.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliAODMCParticle.h"
#include "AliAODPid.h"
#include "AliTPCPIDResponse.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSECharmFraction.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliNormalizationCounter.h"

class TCanvas;
class TTree;
class TChain;
class AliAnalysisTaskSE;


ClassImp(AliAnalysisTaskSECharmFraction)
 
//________________________________________________________________________
  AliAnalysisTaskSECharmFraction::AliAnalysisTaskSECharmFraction() 
    : AliAnalysisTaskSE(),
      fCutsLoose(0),
      fCutsTight(0),
      fFastAnalysis(1),
      fReadMC(kFALSE),
      fsplitMassD0D0bar(kTRUE),
      fLikeSign(kFALSE),
      fusePID(kTRUE),
      fmD0PDG(),
      fnbins(1),
      fptbins(0),
      fNtrMaxforVtx(-1),
      fptAll(),                          
      fptAllSq(),                        
      fptMax(),
      fAcceptanceCuts(),
      fsignalInvMassCut(),
      flargeInvMassCut(),
      fsidebandInvMassCut(),
      fsidebandInvMassWindow(),
      fUseMC(kTRUE),
      fCleanCandOwnVtx(kFALSE),
      fNentries(0),
      fSignalType(0),
      fSignalTypeLsCuts(0),
      fSignalTypeTghCuts(0),
      fCounter(0),
      flistMCproperties(0),
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
      fCutsLoose(0x0),
      fCutsTight(0x0),
      fFastAnalysis(1),
      fReadMC(kFALSE),
      fsplitMassD0D0bar(kTRUE),
      fLikeSign(kFALSE),
      fusePID(kTRUE),
      fmD0PDG(),
      fnbins(1),
      fptbins(0),
      fNtrMaxforVtx(-1),
      fptAll(),                          
      fptAllSq(),                        
      fptMax(),
      fAcceptanceCuts(),
      fsignalInvMassCut(-1.),
      flargeInvMassCut(-1.),
      fsidebandInvMassCut(-1.),
      fsidebandInvMassWindow(-1.),
      fUseMC(kFALSE),
      fCleanCandOwnVtx(kFALSE),
      fNentries(0),
      fSignalType(0),
      fSignalTypeLsCuts(0),
      fSignalTypeTghCuts(0),
      fCounter(0),
      flistMCproperties(0),
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
  fnbins=SetStandardCuts(fptbins);// THIS TO SET NBINS AND BINNING
 
  DefineOutput(1, TH1F::Class());
  DefineOutput(2, TH1F::Class());
  DefineOutput(3, TH1F::Class());
  DefineOutput(4, TH1F::Class());
  DefineOutput(5, AliNormalizationCounter::Class());

  for(Int_t j=6;j<22;j++){
    DefineOutput(j, TList::Class());
  }

  // Output slot for the Cut Objects 
  DefineOutput(22,AliRDHFCutsD0toKpi::Class());  //My private output
  DefineOutput(23,AliRDHFCutsD0toKpi::Class());  //My private output

}


AliAnalysisTaskSECharmFraction::AliAnalysisTaskSECharmFraction(const char *name,AliRDHFCutsD0toKpi *cutsA,AliRDHFCutsD0toKpi *cutsB) 
  : AliAnalysisTaskSE(name),
    fCutsLoose(0),
    fCutsTight(0),
    fFastAnalysis(1),
    fReadMC(kFALSE),
    fsplitMassD0D0bar(kTRUE),
    fLikeSign(kFALSE),
    fusePID(kTRUE),
    fmD0PDG(),
    fnbins(1),
    fptbins(0),
    fNtrMaxforVtx(-1),
    fptAll(),                          
    fptAllSq(),                        
    fptMax(),
    fAcceptanceCuts(),
    fsignalInvMassCut(-1.),
    flargeInvMassCut(-1.),
    fsidebandInvMassCut(-1.),
    fsidebandInvMassWindow(-1.),
    fUseMC(kFALSE),
    fCleanCandOwnVtx(kFALSE),
    fNentries(0),
    fSignalType(0),
    fSignalTypeLsCuts(0),
    fSignalTypeTghCuts(0),
    fCounter(0),
    flistMCproperties(0),
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
  if(fCutsTight){
    delete fCutsTight;fCutsTight=NULL;
  }
  if(fCutsLoose){
    delete fCutsLoose;fCutsLoose=NULL;
  }
  
  //Check consistency between sets of cuts:
  if(cutsA->GetNPtBins()!=cutsB->GetNPtBins()){
    printf("Different number of pt bins between the two sets of cuts: SWITCH TO STANDARD CUTS \n");
    fnbins=SetStandardCuts(fptbins);
  }
  else{
    fCutsTight=new AliRDHFCutsD0toKpi(*cutsA);
    fCutsLoose=new AliRDHFCutsD0toKpi(*cutsB);
    for(Int_t j=0;j<cutsA->GetNPtBins();j++){
      if(TMath::Abs(cutsA->GetPtBinLimits()[j]-cutsB->GetPtBinLimits()[j])>1.e-7){
	printf("Different pt bin limits in the two set of cuts: use the first as reference \n");
	fCutsLoose->SetPtBins(cutsA->GetNPtBins(),cutsA->GetPtBinLimits());
	break;
      }     
    }
    SetPtBins(fCutsTight->GetNPtBins(),fCutsTight->GetPtBinLimits());   
  }
  
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TH1F::Class());
  DefineOutput(2, TH1F::Class());
  DefineOutput(3, TH1F::Class());
  DefineOutput(4, TH1F::Class());
  DefineOutput(5, AliNormalizationCounter::Class());

  for(Int_t j=6;j<22;j++){

    DefineOutput(j, TList::Class());
  }
 // Output slot for the Cut Objects 
  DefineOutput(22,AliRDHFCutsD0toKpi::Class());  //My private output
  DefineOutput(23,AliRDHFCutsD0toKpi::Class());  //My private output
 
}

//________________________________________________________________________
AliAnalysisTaskSECharmFraction::~AliAnalysisTaskSECharmFraction()
{ //Destructor 
  
  if (fCutsTight) {   
    delete fCutsTight;
    fCutsTight = 0;
  }
  if (fCutsLoose) {  
    delete fCutsLoose;
    fCutsLoose = 0;
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
  
  if (fCounter) {
    delete fCounter;
    fCounter = 0;
  } 
  
  if(flistMCproperties){
    flistMCproperties->Delete();
    delete flistMCproperties;
    flistMCproperties=0;
  }
  
  if(flistNoCutsSignal){
    flistNoCutsSignal->Delete();
    delete flistNoCutsSignal;
    flistNoCutsSignal=0;
  }
  if(flistNoCutsBack){
    flistNoCutsBack->Delete();
    delete flistNoCutsBack;
    flistNoCutsBack=0;
  }
  if(flistNoCutsFromB){
    flistNoCutsFromB->Delete();
    delete flistNoCutsFromB;
    flistNoCutsFromB=0;
  }
  if(flistNoCutsFromDstar){
    flistNoCutsFromDstar->Delete();
    delete flistNoCutsFromDstar;
    flistNoCutsFromDstar=0;
  }
  if(flistNoCutsOther){
    flistNoCutsOther->Delete();
    delete flistNoCutsOther;
    flistNoCutsOther=0;
  }
  
 if(flistLsCutsSignal){
   flistLsCutsSignal->Delete();
    delete flistLsCutsSignal;
    flistLsCutsSignal=0;
  }
  if(flistLsCutsBack){
    flistLsCutsBack->Delete();
    delete flistLsCutsBack;
    flistLsCutsBack=0;
  }
  if(flistLsCutsFromB){
    flistLsCutsFromB->Delete();
    delete flistLsCutsFromB;
    flistLsCutsFromB=0;
  }
  if(flistLsCutsFromDstar){
    flistLsCutsFromDstar->Delete();
    delete flistLsCutsFromDstar;
    flistLsCutsFromDstar=0;
  }
  if(flistLsCutsOther){
    flistLsCutsOther->Delete();
    delete flistLsCutsOther;
    flistLsCutsOther=0;
  }
  
 if(flistTghCutsSignal){
   flistTghCutsSignal->Delete();
    delete flistTghCutsSignal;
    flistTghCutsSignal=0;
  }
  if(flistTghCutsBack){
    flistTghCutsBack->Delete();
    delete flistTghCutsBack;
    flistTghCutsBack=0;
  }
  if(flistTghCutsFromB){
    flistTghCutsFromB->Delete();
    delete flistTghCutsFromB;
    flistTghCutsFromB=0;
  }
  if(flistTghCutsFromDstar){
    flistTghCutsFromDstar->Delete();
    delete flistTghCutsFromDstar;
    flistTghCutsFromDstar=0;
  }
  if(flistTghCutsOther){
    flistTghCutsOther->Delete();
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
  
  //  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF.C");
  // gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/D0fromBSetCuts.C");
  // 2 sets of dedicated cuts: fCutsTight is assumed as the standard cut object
  
  // SetAcceptanceCut();
  if(fNtrMaxforVtx<0)SetNMaxTrForVtx(3); //DEFAULT : NO SELECTION
  if(fsignalInvMassCut<0.||flargeInvMassCut<0.||fsidebandInvMassCut<0.||fsidebandInvMassWindow<0.){
    printf("AliAnalysisTaskSECharmFraction: Not All info for mass selection provided: switch to default values \n");
    SetStandardMassSelection();
  }
  
  AliRDHFCutsD0toKpi* copyfCutsTight=new AliRDHFCutsD0toKpi(*fCutsTight);
  const char* nameoutputTight=GetOutputSlot(22)->GetContainer()->GetName();
  copyfCutsTight->SetName(nameoutputTight);
  AliRDHFCutsD0toKpi* copyfCutsLoose=new AliRDHFCutsD0toKpi(*fCutsLoose);
  const char* nameoutputLoose=GetOutputSlot(23)->GetContainer()->GetName();
  copyfCutsLoose->SetName(nameoutputLoose);

  // Post the data
  PostData(22,copyfCutsTight);  
  PostData(23,copyfCutsLoose);
  
  
  fCleanCandOwnVtx=kFALSE;
  if(fCutsTight->GetIsPrimaryWithoutDaughters()^fCutsLoose->GetIsPrimaryWithoutDaughters()) {
    printf("Two cut objects have different selection for primary vertex recalculation w/o daughters:\n Dangerous for variable drawing!! \n");   
  }
  else{
    if(fCutsTight->GetIsPrimaryWithoutDaughters()){
      fCleanCandOwnVtx=kTRUE;
      fCutsTight->SetRemoveDaughtersFromPrim(kFALSE);
      fCutsLoose->SetRemoveDaughtersFromPrim(kFALSE);
    }
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
  Printf("INSIDE USER CREATE \n");
  
  // fNentries=new TH1F("nentriesChFr", "Look at the number of entries! = number of AODs", 2,1.,2.);
  
  fNentries=new TH1F("nentriesChFr", "Analyzed sample properties", 21,-0.5,20.5);

  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");

  fNentries->GetXaxis()->SetBinLabel(2,"nEvTGHTsel");
  fNentries->GetXaxis()->SetBinLabel(3,"nEvTGHTPile-up Rej");
  fNentries->GetXaxis()->SetBinLabel(4,"nEvTGHTGoodVtxS");
  fNentries->GetXaxis()->SetBinLabel(5,"nEvTGHTRejVtxZ");
  fNentries->GetXaxis()->SetBinLabel(6,"nTracksTGHTEv");
  fNentries->GetXaxis()->SetBinLabel(7,"nCandTGHTEv");
  fNentries->GetXaxis()->SetBinLabel(8,"nCandSelTGHTEv");
  fNentries->GetXaxis()->SetBinLabel(20,"nUnexpErrorTGHT");

  fNentries->GetXaxis()->SetBinLabel(9,"nEvLSsel");
  fNentries->GetXaxis()->SetBinLabel(10,"nEvLSPile-up Rej");
  fNentries->GetXaxis()->SetBinLabel(11,"nEvLSGoodVtxS");
  fNentries->GetXaxis()->SetBinLabel(12,"nEvLSRejVtxZ");
  fNentries->GetXaxis()->SetBinLabel(13,"nTracksLSEv");
  fNentries->GetXaxis()->SetBinLabel(14,"nCandLSEv");
  fNentries->GetXaxis()->SetBinLabel(15,"nCandSelLSEv");
  fNentries->GetXaxis()->SetBinLabel(21,"nUnexpErrorTGHT");

  /*   -----------------  NOT ACTIVATED YET ------------------
    fNentries->GetXaxis()->SetBinLabel(5,"ptbin = -1");
    fNentries->GetXaxis()->SetBinLabel(6,"no daughter");
    fNentries->GetXaxis()->SetBinLabel(7,"nCandSel(Tr)");
    fNentries->GetXaxis()->SetBinLabel(8,"PID=0");
    fNentries->GetXaxis()->SetBinLabel(9,"PID=1");
    fNentries->GetXaxis()->SetBinLabel(10,"PID=2");
    fNentries->GetXaxis()->SetBinLabel(11,"PID=3");
    fNentries->GetXaxis()->SetBinLabel(12,"K");
    fNentries->GetXaxis()->SetBinLabel(13,"Lambda");
    fNentries->GetXaxis()->SetBinLabel(14,"Pile-up Rej");
    fNentries->GetXaxis()->SetBinLabel(15,"N. of 0SMH");
  */

  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  fSignalType=new TH1F("hsignaltype", "Histo for type of MC signal", 61,-1.,60.);
  fSignalTypeLsCuts=new TH1F("hsignaltypeLsCuts", "Histo for type of MC signal with loose cuts", 61,-1.,60.);
  fSignalTypeTghCuts=new TH1F("hsignaltypeTghCuts", "Histo for type of MC signal with tight cuts", 61,-1.,60.);



  fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(5)->GetContainer()->GetName()));
 
  //##########  DEFINE THE TLISTS ##################
  flistMCproperties=new TList();
  flistMCproperties->SetOwner();
  flistMCproperties->SetName("listMCproperties");

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



  Float_t ptbinsD0arr[35]={0.,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.,1.25,1.5,1.75,2.,2.3,2.6,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,9.,10.,12.,14.,16.,20.,25.,30.,40.,50.,100.};
  Float_t dumbinning[201];
  for(Int_t j=0;j<201;j++){
    dumbinning[j]=(Float_t)j*0.5;
  }

  // DEFINE EDGES FOR SPARSE HISTOS
  const Int_t nPtbinsForSparse=91;//nuber of edges, -1 to get number of bins
  Double_t ptbinsForNsparse[nPtbinsForSparse];//Binning in pt: step: 0.2 GeV/c up to 8 GeV/c, 0.5 Upto 20, 1 GeV/c up to 40, 5 upt to 70-> 8/0.2+12/0.5+20/1.+30./5
  Double_t pT=0;
  Double_t massbins[186],impparbins[401];  
  Double_t massHypoBins[4]={1.,2.,3.,4.};
  Int_t nbinsSparse[5]={185,185,nPtbinsForSparse-1,400,3};
  for(Int_t nBins=0;nBins<nPtbinsForSparse;nBins++){    
    ptbinsForNsparse[nBins]=pT;
    if(pT<8.)pT+=0.2;
    else if(pT<20)pT+=0.5;
    else if(pT<40)pT+=1;
    else if(pT<70)pT+=5.;  
  }  
  for(Int_t nBins=0;nBins<186;nBins++){
    massbins[nBins]=1.680+nBins*(2.050-1.680)/185.;
  }
  for(Int_t nBins=0;nBins<401;nBins++){
    impparbins[nBins]=-1000+nBins*(2000.)/400.;
  }

  
  // Lxy and CosPointXY study
  Int_t nbinsSparsCxyLxy[4]={84,fCutsTight->GetNPtBins(),10,25}; 
  Double_t binLowLimitSparseCxyLxy[4]={1.680,fCutsTight->GetPtBinLimits()[0],0.99,0.};// Use OverFlow/UnderFlow to get other cases
  Double_t binUpLimitSparseCxyLxy[4]={2.100,fCutsTight->GetPtBinLimits()[fCutsTight->GetNPtBins()],1.,50.};
  Double_t *ptbinlimitsCxyLxy=new Double_t[fCutsTight->GetNPtBins()+1];
  for(Int_t nBins=0;nBins<=fCutsTight->GetNPtBins();nBins++){
    ptbinlimitsCxyLxy[nBins]=fCutsTight->GetPtBinLimits()[nBins];
  }
  

  //################################################################################################
  //                                                                                               #
  //                HISTO FOR MC PROPERTIES OF D0, c quarks and B mesons                           #
  //                                                                                               # 
  //################################################################################################
  TH1F *hMCcquarkAllPt=new TH1F("hMCcquarkAllPt","c quark Pt (all cquarks produced)",34,ptbinsD0arr);
  TH1F *hMCcquarkAllEta=new TH1F("hMCcquarkAllEta","c quark Eta (all cquarks produced)",50,-3.,3.);
  TH1F *hMCcquarkAllEnergy=new TH1F("hMCcquarkAllEnergy","c quark Pt (all cquarks produced)",200,0.,100.);
  TH1F *hMCcquarkNdaught=new TH1F("hMCcquarkNdaught","N cquark daughters (all cquarks produced)",100,0.,100.);
  TH1F *hMCD0fromcPt=new TH1F("hMCD0fromcPt","D0 from c Pt",34,ptbinsD0arr);
  TH1F *hMCD0fromcEta=new TH1F("hMCD0fromcEta","D0 from c Eta",50,-3.,3.);
  TH1F *hMCD0fromcEnergy=new TH1F("hMCD0fromcEnergy","D0 from c Energy",200,0.,100.);
  
  TH2F *hMCD0VscquarkPt=new TH2F("hMCD0VscquarkPt","D0 pt Vs cquark pt",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hMCD0VscquarkEnergy=new TH2F("hMCD0VscquarkEnergy","D0 Energy Vs cquark Energy",200,0.,50.,200,0.,50.);
  TH1F *hMCD0deltacquarkEnergy=new TH1F("hMCD0deltacquarkEnergy","Fractional D0 Energy w.r.t. cquark Energy",20,0.,1.);
  TH1F *hMCD0EnergyVsAvcquarkDaughtEn=new TH1F("hMCD0EnergyVsAvcquarkDaughtEn","#Delta(E^{D^0}-E_{avg})/E_{cquark}",40,-1.,1.);
  TH1F *hMCD0cquarkAngle=new TH1F("hMCD0cquarkAngle","cosine of the angle between D0 and c quark particle",40,-1.,1.);
  TH2F *hMCD0cquarkAngleEnergy=new TH2F("hMCD0cquarkAngleEnergy","cosine of the angle between D0 and c quark particle as a function of Energy",25,0.,50.,40,-1.,1.);

  TH1I *hMCfromBpdgB=new TH1I("hMCfromBpdgB","hMCfromBpdgB",10000,0.,10000);
  TH1F *hMCBhadrPt=new TH1F("hMCBhadrPt","B hadr Pt",34,ptbinsD0arr);
  TH1F *hMCBhadrEta=new TH1F("hMCBhadrEta","B hadr Eta",50,-3.,3.);
  TH1F *hMCBhadrEnergy=new TH1F("hMCBhadrEnergy","B hadr Pt",200,0.,100.);
  TH1F *hMCBhadrNdaught=new TH1F("hMCBhadrNdaught","N Bhadr daughters",100,0.,100.);
  TH1F *hMCD0fromBPt=new TH1F("hMCD0fromBPt","D0 from B Pt",34,ptbinsD0arr);
  TH1F *hMCD0fromBEta=new TH1F("hMCD0fromBEta","D0 from B Eta",50,-3.,3.);
  TH1F *hMCD0fromBEnergy=new TH1F("hMCD0fromBEnergy","D0 from B Energy",200,0.,100.);

  TH2F *hMCD0VsBhadrPt=new TH2F("hMCD0VsBhadrPt","D0 pt Vs Bhadr pt",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hMCD0VsBhadrEnergy=new TH2F("hMCD0VsBhadrEnergy","D0 Energy Vs Bhadr Energy",200,0.,50.,200,0.,50.);
  TH1F *hMCD0deltaBhadrEnergy=new TH1F("hMCD0deltaBhadrEnergy","Fractional D0 Energy w.r.t. Bhadr Energy",20,0.,1.);
  TH1F *hMCD0EnergyVsAvBDaughtEn=new TH1F("hMCD0EnergyVsAvBDaughtEn","#Delta(E^{D^0}-E_{avg})/E_{Bahdr}",40,-1.,1.);
  TH1F *hMCD0BhadrAngle=new TH1F("hMCD0BhadrAngle","cosine of the angle between D0 and Bhadr particle",40,-1.,1.);
  TH2F *hMCD0BhadrAngleEnergy=new TH2F("hMCD0BhadrAngleEnergy","cosine of the angle between D0 and Bhadr particle as a function of Energy",25,0.,50.,40,-1.,1.);

  TH1I *hMCPartFound=new TH1I("hMCPartFound","1=c,2=D0,3=fromBall,4=fromBmeson,5=fromBbaryon",6,0,6); 
 

  flistMCproperties->Add(hMCcquarkAllPt);
  flistMCproperties->Add(hMCcquarkAllEta);
  flistMCproperties->Add(hMCcquarkAllEnergy);
  flistMCproperties->Add(hMCcquarkNdaught);
  flistMCproperties->Add(hMCD0fromcPt);
  flistMCproperties->Add(hMCD0fromcEta);
  flistMCproperties->Add(hMCD0fromcEnergy);
  flistMCproperties->Add(hMCD0VscquarkPt);
  flistMCproperties->Add(hMCD0VscquarkEnergy);
  flistMCproperties->Add(hMCD0deltacquarkEnergy);
  flistMCproperties->Add(hMCD0EnergyVsAvcquarkDaughtEn);
  flistMCproperties->Add(hMCD0cquarkAngle);
  flistMCproperties->Add(hMCD0cquarkAngleEnergy);
  
  flistMCproperties->Add(hMCfromBpdgB);
  flistMCproperties->Add(hMCBhadrPt);
  flistMCproperties->Add(hMCBhadrEta);
  flistMCproperties->Add(hMCBhadrEnergy);
  flistMCproperties->Add(hMCBhadrNdaught);
  flistMCproperties->Add(hMCD0fromBPt);
  flistMCproperties->Add(hMCD0fromBEta);
  flistMCproperties->Add(hMCD0fromBEnergy);
  flistMCproperties->Add(hMCD0VsBhadrPt);
  flistMCproperties->Add(hMCD0VsBhadrEnergy);
  flistMCproperties->Add(hMCD0deltaBhadrEnergy);
  flistMCproperties->Add(hMCD0EnergyVsAvBDaughtEn);
  flistMCproperties->Add(hMCD0BhadrAngle);
  flistMCproperties->Add(hMCD0BhadrAngleEnergy);
  flistMCproperties->Add(hMCPartFound);

  //################################################################################################
  //                                                                                               #
  //                         HISTOS FOR NO CUTS CASE                                               #
  //                                                                                               #
  //################################################################################################
  Printf("AFTER MC HISTOS \n");

  //############ NO CUTS SIGNAL HISTOGRAMS ###############
  //
  // ####### global properties histo ############

  TH2F *hCPtaVSd0d0NCsign=new TH2F("hCPtaVSd0d0NCsign","hCPtaVSd0d0_NoCuts_Signal",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZNCsign=new TH1F("hSecVtxZNCsign","hSecVtxZ_NoCuts_Signal",1000,-8.,8.);
  TH1F *hSecVtxXNCsign=new TH1F("hSecVtxXNCsign","hSecVtxX_NoCuts_Signal",1000,-3000.,3000.);
  TH1F *hSecVtxYNCsign=new TH1F("hSecVtxYNCsign","hSecVtxY_NoCuts_Signal",1000,-3000.,3000.);
  TH2F *hSecVtxXYNCsign=new TH2F("hSecVtxXYNCsign","hSecVtxXY_NoCuts_Signal",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiNCsign=new TH1F("hSecVtxPhiNCsign","hSecVtxPhi_NoCuts_Signal",180,-180.1,180.1);
  TH1F *hd0singlTrackNCsign=new TH1F("hd0singlTrackNCsign","hd0singlTrackNoCuts_Signal",1000,-5000.,5000.);
  TH1F *hCPtaNCsign=new TH1F("hCPtaNCsign","hCPta_NoCuts_Signal",100,-1.,1.);
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
  flistNoCutsSignal->Add(hd0singlTrackNCsign);
  flistNoCutsSignal->Add(hCPtaNCsign);
  flistNoCutsSignal->Add(hd0xd0NCsign);
  flistNoCutsSignal->Add(hMassTrueNCsign);
  flistNoCutsSignal->Add(hMassNCsign);
  flistNoCutsSignal->Add(hMassTrueNCsignPM);
  flistNoCutsSignal->Add(hMassNCsignPM);
  flistNoCutsSignal->Add(hMassTrueNCsignSB);
  flistNoCutsSignal->Add(hMassNCsignSB);
  
  //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaNCsign=new TH1F("hdcaNCsign","hdca_NoCuts_Signal",100,0.,1000.);
  hdcaNCsign->SetXTitle("dca   [#mum]");
  hdcaNCsign->SetYTitle("Entries");
  TH1F *hcosthetastarNCsign=new TH1F("hcosthetastarNCsign","hCosThetaStar_NoCuts_Signal",50,-1.,1.);
  hcosthetastarNCsign->SetXTitle("cos #theta^{*}");
  hcosthetastarNCsign->SetYTitle("Entries");
  TH1F *hptD0NCsign=new TH1F("hptD0NCsign","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0NCsign->SetXTitle("p_{t}  [GeV/c]");
  hptD0NCsign->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtNCsign=new TH1F("hptD0VsMaxPtNCsign","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtNCsign=new TH2F("hptD0PTallsqrtNCsign","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallNCsign=new TH2F("hptD0PTallNCsign","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBNCsign=new TH2F("hptD0vsptBNCsign","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBNCsign=new TH2F("hpD0vspBNCsign","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkNCsign=new TH2F("hptD0vsptcquarkNCsign","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkNCsign=new TH2F("hpD0vspcquarkNCsign","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistNoCutsSignal->Add(hdcaNCsign);
  flistNoCutsSignal->Add(hcosthetastarNCsign);
  flistNoCutsSignal->Add(hptD0NCsign);
  flistNoCutsSignal->Add(hptD0VsMaxPtNCsign);
  flistNoCutsSignal->Add(hptD0PTallsqrtNCsign);
  flistNoCutsSignal->Add(hptD0PTallNCsign);
  flistNoCutsSignal->Add(hptD0vsptBNCsign);
  flistNoCutsSignal->Add(hpD0vspBNCsign);
  flistNoCutsSignal->Add(hptD0vsptcquarkNCsign);
  flistNoCutsSignal->Add(hpD0vspcquarkNCsign);

  TH1F *hd0zD0ptNCsign;
  TH1F *hInvMassD0NCsign,*hInvMassD0barNCsign;
  TH2F *hInvMassPtNCsign=new TH2F("hInvMassPtNCsign","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
  flistNoCutsSignal->Add(hInvMassPtNCsign);
  THnSparseF *hSparseNCsign=new THnSparseF("hSparseNCsign","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseNCsign->SetBinEdges(0,massbins);
  hSparseNCsign->SetBinEdges(1,massbins);
  hSparseNCsign->SetBinEdges(2,ptbinsForNsparse);
  hSparseNCsign->SetBinEdges(3,impparbins);
  hSparseNCsign->SetBinEdges(4,massHypoBins); 
  flistNoCutsSignal->Add(hSparseNCsign);




  THnSparseF *hSparseCxyLxyNCsign=new THnSparseF("hSparseCxyLxyNCsign","Candidate Mass;massD0;Pt;CosXY;Lxy",4,nbinsSparsCxyLxy,binLowLimitSparseCxyLxy,binUpLimitSparseCxyLxy);
  hSparseCxyLxyNCsign->SetBinEdges(1,ptbinlimitsCxyLxy);
  hSparseCxyLxyNCsign->GetAxis(0)->SetName("mass");
  hSparseCxyLxyNCsign->GetAxis(0)->SetTitle("Invariant Mass (K#pi) [GeV/c^{2}]");
  hSparseCxyLxyNCsign->GetAxis(1)->SetName("pt");
  hSparseCxyLxyNCsign->GetAxis(1)->SetTitle("p_{t} [GeV/c]");
  hSparseCxyLxyNCsign->GetAxis(2)->SetName("CosPointXY");
  hSparseCxyLxyNCsign->GetAxis(2)->SetTitle("Cos#theta_{point}^{XY}");
  hSparseCxyLxyNCsign->GetAxis(3)->SetName("NormDecLengthXY");
  hSparseCxyLxyNCsign->GetAxis(3)->SetTitle("Normalized XY decay length");
  
  flistNoCutsSignal->Add(hSparseCxyLxyNCsign);

  

  TH1F *hetaNCsign;
  TH1F *hCosPDPBNCsign;
  TH1F *hCosPcPDNCsign;
  // ADDITIONAL HISTOS
  TH2F *hd0D0VSd0xd0NCsignpt;
  TH2F *hangletracksVSd0xd0NCsignpt;
  TH2F *hangletracksVSd0D0NCsignpt;
  TH1F *hd0xd0NCsignpt;
  // AZIMUHAL HISTOS
  TH1F *hPhiHistPMNCsignpt,*hPhiHistSBNCsignpt;

  
  TH2F *hTOFpidNCsign=new TH2F("hTOFpidNCsign","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistNoCutsSignal->Add(hTOFpidNCsign);
 
  

  //##################
  for(Int_t i=0;i<fnbins;i++){
    //Printf("INSIDE HISTOS CREATION LOOP: %d \n",fnbins);
    
    namehist="hPhiHistPMNCsign_pt";
    namehist+=i;
    titlehist="Azimuthal correlation No Cuts Sign PM ptbin=";
    titlehist+=i;
    hPhiHistPMNCsignpt=new TH1F(namehist.Data(),titlehist.Data(),100,-3.15,3.15);
    hPhiHistPMNCsignpt->Sumw2();
    flistNoCutsSignal->Add(hPhiHistPMNCsignpt);

    namehist="hPhiHistSBNCsign_pt";
    namehist+=i;
    titlehist="Azimuthal correlation No Cuts Sign SB ptbin=";
    titlehist+=i;
    hPhiHistSBNCsignpt=new TH1F(namehist.Data(),titlehist.Data(),100,-3.15,3.15);
    hPhiHistSBNCsignpt->Sumw2();
    flistNoCutsSignal->Add(hPhiHistSBNCsignpt);
    

    namehist="hd0zD0ptNCsign_pt";
    namehist+=i;
    titlehist="d0(z) No Cuts Signalm ptbin=";
    titlehist+=i;
    hd0zD0ptNCsign=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptNCsign->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptNCsign->SetYTitle("Entries");
    flistNoCutsSignal->Add(hd0zD0ptNCsign);

    namehist="hInvMassD0NCsign_pt";
    namehist+=i;
    titlehist="Invariant Mass D0 No Cuts Signal ptbin=";
    titlehist+=i;
    hInvMassD0NCsign=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0NCsign->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0NCsign->SetYTitle("Entries");
    flistNoCutsSignal->Add(hInvMassD0NCsign);


    namehist="hInvMassD0barNCsign_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar No Cuts Signal ptbin=";
    titlehist+=i;
    hInvMassD0barNCsign=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barNCsign->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barNCsign->SetYTitle("Entries");
    flistNoCutsSignal->Add(hInvMassD0barNCsign);


    namehist="hetaNCsign_pt";
    namehist+=i;
    titlehist="eta No Cuts Signal ptbin=";
    titlehist+=i;
    hetaNCsign=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaNCsign->SetXTitle("Pseudorapidity");
    hetaNCsign->SetYTitle("Entries");
    flistNoCutsSignal->Add(hetaNCsign);

    namehist="hCosPDPBNCsign_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBNCsign=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBNCsign->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBNCsign->SetYTitle("Entries");
    flistNoCutsSignal->Add(hCosPDPBNCsign);

    namehist="hCosPcPDNCsign_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDNCsign=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDNCsign->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDNCsign->SetYTitle("Entries");
    flistNoCutsSignal->Add(hCosPcPDNCsign);
    

    // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0NCsign_pt";
    namehist+=i;
    titlehist="d0xd0 No Cuts Signal ptbin=";
    titlehist+=i;
    hd0xd0NCsignpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0NCsignpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0NCsignpt->SetYTitle("Entries");
    flistNoCutsSignal->Add(hd0xd0NCsignpt);


    namehist="hd0D0VSd0xd0NCsign_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} No Cuts Signal ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0NCsignpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0NCsignpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0NCsignpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistNoCutsSignal->Add(hd0D0VSd0xd0NCsignpt);
    
    
    namehist="hangletracksVSd0xd0NCsign_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} No Cuts Signal ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0NCsignpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0NCsignpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0NCsignpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistNoCutsSignal->Add(hangletracksVSd0xd0NCsignpt);
    

    namehist="hangletracksVSd0D0NCsign_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} No Cuts Signal ptbin=";
    titlehist+=i;
    hangletracksVSd0D0NCsignpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0NCsignpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0NCsignpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistNoCutsSignal->Add(hangletracksVSd0D0NCsignpt);

  }
  Printf("AFTER LOOP HISTOS CREATION \n");
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  
  TH1F *hd0D0ptNCsignPM;
  TH1F *hMCd0D0ptNCsignPM;
  TH1F *hd0D0VtxTrueptNCsignPM;
  TH1F *hd0D0ptNCsignSB;
  TH1F *hMCd0D0ptNCsignSB;
  TH1F *hd0D0VtxTrueptNCsignSB;
  namehist="hd0D0ptNCsign_";
  titlehist="D^{0} impact par. plot, No Cuts, Signal, ";
  for(Int_t i=0;i<fnbins;i++){
    //Printf("IN HISTOS CREATION USING PTBINS VALUES for NAMES \n");
    strnamept=namehist;
    strnamept.Append("PkMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Mass Peak, ");

    strtitlept+=fptbins[i];
    //Printf("IN HISTOS CREATION USING PTBINS VALUES for NAMES %d: %f\n",i,fptbins[i]);
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCsignPM= new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCsignPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCsignPM->SetYTitle("Entries");
    flistNoCutsSignal->Add(hd0D0ptNCsignPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCsignPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCsignPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCsignPM->SetYTitle("Entries");
    flistNoCutsSignal->Add(hMCd0D0ptNCsignPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCsignPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCsignPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCsignPM->SetYTitle("Entries");
    flistNoCutsSignal->Add(hd0D0VtxTrueptNCsignPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCsignSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCsignSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCsignSB->SetYTitle("Entries");
    flistNoCutsSignal->Add(hd0D0ptNCsignSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCsignSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCsignSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCsignSB->SetYTitle("Entries");
    flistNoCutsSignal->Add(hMCd0D0ptNCsignSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCsignSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCsignSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCsignSB->SetYTitle("Entries");
    flistNoCutsSignal->Add(hd0D0VtxTrueptNCsignSB);
  }

  //Printf("AFTER SIGNAL HISTOS CREATION for NOCUTS\n");


  //############ NO CUTS BACKGROUND HISTOGRAMS ###########
  //
  //   ######## global properties histos #######
  TH2F *hCPtaVSd0d0NCback=new TH2F("hCPtaVSd0d0NCback","hCPtaVSd0d0_NoCuts_Background",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZNCback=new TH1F("hSecVtxZNCback","hSecVtxZ_NoCuts_Background",1000,-8.,8.);
  TH1F *hSecVtxXNCback=new TH1F("hSecVtxXNCback","hSecVtxX_NoCuts_Background",1000,-3000.,3000.);
  TH1F *hSecVtxYNCback=new TH1F("hSecVtxYNCback","hSecVtxY_NoCuts_Background",1000,-3000.,3000.);
  TH2F *hSecVtxXYNCback=new TH2F("hSecVtxXYNCback","hSecVtxXY_NoCuts_Background",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiNCback=new TH1F("hSecVtxPhiNCback","hSecVtxPhi_NoCuts_Background",180,-180.1,180.1);
  TH1F *hd0singlTrackNCback=new TH1F("hd0singlTrackNCback","hd0singlTrackNoCuts_Back",1000,-5000.,5000.);
  TH1F *hCPtaNCback=new TH1F("hCPtaNCback","hCPta_NoCuts_Background",100,-1.,1.);
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
  flistNoCutsBack->Add(hd0singlTrackNCback);
  flistNoCutsBack->Add(hCPtaNCback);
  flistNoCutsBack->Add(hd0xd0NCback);
  flistNoCutsBack->Add(hMassTrueNCback);
  flistNoCutsBack->Add(hMassNCback);
  flistNoCutsBack->Add(hMassTrueNCbackPM);
  flistNoCutsBack->Add(hMassNCbackPM);
  flistNoCutsBack->Add(hMassTrueNCbackSB);
  flistNoCutsBack->Add(hMassNCbackSB);


 //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaNCback=new TH1F("hdcaNCback","hdca_NoCuts_Backgr",100,0.,1000.);
  hdcaNCback->SetXTitle("dca   [#mum]");
  hdcaNCback->SetYTitle("Entries");
  TH1F *hcosthetastarNCback=new TH1F("hcosthetastarNCback","hCosThetaStar_NoCuts_Backgr",50,-1.,1.);
  hcosthetastarNCback->SetXTitle("cos #theta^{*}");
  hcosthetastarNCback->SetYTitle("Entries");
  TH1F *hptD0NCback=new TH1F("hptD0NCback","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0NCback->SetXTitle("p_{t}  [GeV/c]");
  hptD0NCback->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtNCback=new TH1F("hptD0VsMaxPtNCback","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtNCback=new TH2F("hptD0PTallsqrtNCback","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallNCback=new TH2F("hptD0PTallNCback","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBNCback=new TH2F("hptD0vsptBNCback","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBNCback=new TH2F("hpD0vspBNCback","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkNCback=new TH2F("hptD0vsptcquarkNCback","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkNCback=new TH2F("hpD0vspcquarkNCback","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistNoCutsBack->Add(hdcaNCback);
  flistNoCutsBack->Add(hcosthetastarNCback);
  flistNoCutsBack->Add(hptD0NCback);
  flistNoCutsBack->Add(hptD0VsMaxPtNCback);
  flistNoCutsBack->Add(hptD0PTallsqrtNCback);
  flistNoCutsBack->Add(hptD0PTallNCback);
  flistNoCutsBack->Add(hptD0vsptBNCback);
  flistNoCutsBack->Add(hpD0vspBNCback);
  flistNoCutsBack->Add(hptD0vsptcquarkNCback);
  flistNoCutsBack->Add(hpD0vspcquarkNCback);
 
  TH1F *hd0zD0ptNCback;
  TH1F *hInvMassD0NCback,*hInvMassD0barNCback;
  TH2F *hInvMassPtNCback=new TH2F("hInvMassPtNCback","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
  THnSparseF *hSparseNCback=new THnSparseF("hSparseNCback","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseNCback->SetBinEdges(0,massbins);
  hSparseNCback->SetBinEdges(1,massbins);
  hSparseNCback->SetBinEdges(2,ptbinsForNsparse);
  hSparseNCback->SetBinEdges(3,impparbins);
  hSparseNCback->SetBinEdges(4,massHypoBins); 
  flistNoCutsBack->Add(hSparseNCback);

  TH1F *hetaNCback;
  TH1F *hCosPDPBNCback;
  TH1F *hCosPcPDNCback;
   // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0NCbackpt;
  TH2F *hangletracksVSd0xd0NCbackpt;
  TH2F *hangletracksVSd0D0NCbackpt;
  TH1F *hd0xd0NCbackpt;
  flistNoCutsBack->Add(hInvMassPtNCback);

  TH2F *hTOFpidNCback=new TH2F("hTOFpidNCback","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistNoCutsBack->Add(hTOFpidNCback);

  for(Int_t i=0;i<fnbins;i++){
    namehist="hd0zD0ptNCback_pt";
    namehist+=i;
    titlehist="d0(z) No Cuts Backgrm ptbin=";
    titlehist+=i;
    hd0zD0ptNCback=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptNCback->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptNCback->SetYTitle("Entries");
    flistNoCutsBack->Add(hd0zD0ptNCback);

    namehist="hInvMassD0NCback_pt";
    namehist+=i;
    titlehist="Invariant Mass No Cuts Backgr ptbin=";
    titlehist+=i;
    hInvMassD0NCback=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0NCback->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0NCback->SetYTitle("Entries");
    flistNoCutsBack->Add(hInvMassD0NCback);

    
    namehist="hInvMassD0barNCback_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar No Cuts Back ptbin=";
    titlehist+=i;
    hInvMassD0barNCback=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barNCback->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barNCback->SetYTitle("Entries");
    flistNoCutsBack->Add(hInvMassD0barNCback);


    namehist="hetaNCback_pt";
    namehist+=i;
    titlehist="eta No Cuts Backgr ptbin=";
    titlehist+=i;
    hetaNCback=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaNCback->SetXTitle("Pseudorapidity");
    hetaNCback->SetYTitle("Entries");
    flistNoCutsBack->Add(hetaNCback);

    namehist="hCosPDPBNCback_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBNCback=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBNCback->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBNCback->SetYTitle("Entries");
    flistNoCutsBack->Add(hCosPDPBNCback);

    namehist="hCosPcPDNCback_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDNCback=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDNCback->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDNCback->SetYTitle("Entries");
    flistNoCutsBack->Add(hCosPcPDNCback);


    // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0NCback_pt";
    namehist+=i;
    titlehist="d0xd0 No Cuts Background ptbin=";
    titlehist+=i;
    hd0xd0NCbackpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0NCbackpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0NCbackpt->SetYTitle("Entries");
    flistNoCutsBack->Add(hd0xd0NCbackpt);


    namehist="hd0D0VSd0xd0NCback_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} No Cuts Back ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0NCbackpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0NCbackpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0NCbackpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistNoCutsBack->Add(hd0D0VSd0xd0NCbackpt);
    
    
    namehist="hangletracksVSd0xd0NCback_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} No Cuts Back ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0NCbackpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0NCbackpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0NCbackpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistNoCutsBack->Add(hangletracksVSd0xd0NCbackpt);
    

    namehist="hangletracksVSd0D0NCback_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} No Cuts Back ptbin=";
    titlehist+=i;
    hangletracksVSd0D0NCbackpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0NCbackpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0NCbackpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistNoCutsBack->Add(hangletracksVSd0D0NCbackpt);


    
  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
  
  TH1F *hd0D0ptNCbackPM;
  TH1F *hMCd0D0ptNCbackPM;
  TH1F *hd0D0VtxTrueptNCbackPM;
  TH1F *hd0D0ptNCbackSB;
  TH1F *hMCd0D0ptNCbackSB;
  TH1F *hd0D0VtxTrueptNCbackSB;
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
    
    hd0D0ptNCbackPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCbackPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCbackPM->SetYTitle("Entries");
    flistNoCutsBack->Add(hd0D0ptNCbackPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCbackPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCbackPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCbackPM->SetYTitle("Entries");
    flistNoCutsBack->Add(hMCd0D0ptNCbackPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCbackPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCbackPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCbackPM->SetYTitle("Entries");
    flistNoCutsBack->Add(hd0D0VtxTrueptNCbackPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCbackSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCbackSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCbackSB->SetYTitle("Entries");
    flistNoCutsBack->Add(hd0D0ptNCbackSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCbackSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCbackSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCbackSB->SetYTitle("Entries");
    flistNoCutsBack->Add(hMCd0D0ptNCbackSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCbackSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCbackSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCbackSB->SetYTitle("Entries");
    flistNoCutsBack->Add(hd0D0VtxTrueptNCbackSB);
  }



 //############ NO CUTS FROMB HISTOGRAMS ###########
  //
  //#######  global properties histos

  TH2F *hCPtaVSd0d0NCfromB=new TH2F("hCPtaVSd0d0NCfromB","hCPtaVSd0d0_NoCuts_FromB",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZNCfromB=new TH1F("hSecVtxZNCfromB","hSecVtxZ_NoCuts_FromB",1000,-8.,8.);
  TH1F *hSecVtxXNCfromB=new TH1F("hSecVtxXNCfromB","hSecVtxX_NoCuts_FromB",1000,-3000.,3000.);
  TH1F *hSecVtxYNCfromB=new TH1F("hSecVtxYNCfromB","hSecVtxY_NoCuts_FromB",1000,-3000.,3000.);
  TH2F *hSecVtxXYNCfromB=new TH2F("hSecVtxXYNCfromB","hSecVtxXY_NoCuts_FromB",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiNCfromB=new TH1F("hSecVtxPhiNCfromB","hSecVtxPhi_NoCuts_FromB",180,-180.1,180.1);
  TH1F *hd0singlTrackNCfromB=new TH1F("hd0singlTrackNCfromB","hd0singlTrackNoCuts_FromB",1000,-5000.,5000.);
  TH1F *hCPtaNCfromB=new TH1F("hCPtaNCfromB","hCPta_NoCuts_FromB",100,-1.,1.);
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
  flistNoCutsFromB->Add(hd0singlTrackNCfromB);
  flistNoCutsFromB->Add(hCPtaNCfromB);
  flistNoCutsFromB->Add(hd0xd0NCfromB);
  flistNoCutsFromB->Add(hMassTrueNCfromB);
  flistNoCutsFromB->Add(hMassNCfromB);
  flistNoCutsFromB->Add(hMassTrueNCfromBPM);
  flistNoCutsFromB->Add(hMassNCfromBPM);
  flistNoCutsFromB->Add(hMassTrueNCfromBSB);
  flistNoCutsFromB->Add(hMassNCfromBSB);





 //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaNCfromB=new TH1F("hdcaNCfromB","hdca_NoCuts_FromB",100,0.,1000.);
  hdcaNCfromB->SetXTitle("dca   [#mum]");
  hdcaNCfromB->SetYTitle("Entries");
  TH1F *hcosthetastarNCfromB=new TH1F("hcosthetastarNCfromB","hCosThetaStar_NoCuts_FromB",50,-1.,1.);
  hcosthetastarNCfromB->SetXTitle("cos #theta^{*}");
  hcosthetastarNCfromB->SetYTitle("Entries");
  TH1F *hptD0NCfromB=new TH1F("hptD0NCfromB","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0NCfromB->SetXTitle("p_{t}  [GeV/c]");
  hptD0NCfromB->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtNCfromB=new TH1F("hptD0VsMaxPtNCfromB","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtNCfromB=new TH2F("hptD0PTallsqrtNCfromB","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallNCfromB=new TH2F("hptD0PTallNCfromB","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBNCfromB=new TH2F("hptD0vsptBNCfromB","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBNCfromB=new TH2F("hpD0vspBNCfromB","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkNCfromB=new TH2F("hptD0vsptcquarkNCfromB","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkNCfromB=new TH2F("hpD0vspcquarkNCfromB","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistNoCutsFromB->Add(hdcaNCfromB);
  flistNoCutsFromB->Add(hcosthetastarNCfromB);
  flistNoCutsFromB->Add(hptD0NCfromB);
  flistNoCutsFromB->Add(hptD0VsMaxPtNCfromB);
  flistNoCutsFromB->Add(hptD0PTallsqrtNCfromB);
  flistNoCutsFromB->Add(hptD0PTallNCfromB);
  flistNoCutsFromB->Add(hptD0vsptBNCfromB);
  flistNoCutsFromB->Add(hpD0vspBNCfromB);
  flistNoCutsFromB->Add(hptD0vsptcquarkNCfromB);
  flistNoCutsFromB->Add(hpD0vspcquarkNCfromB);
 
  TH1F *hd0zD0ptNCfromB;
  TH1F *hInvMassD0NCfromB,*hInvMassD0barNCfromB;
  TH2F *hInvMassPtNCfromB=new TH2F("hInvMassPtNCfromB","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
  THnSparseF *hSparseNCfromB=new THnSparseF("hSparseNCfromB","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseNCfromB->SetBinEdges(0,massbins);
  hSparseNCfromB->SetBinEdges(1,massbins);
  hSparseNCfromB->SetBinEdges(2,ptbinsForNsparse);
  hSparseNCfromB->SetBinEdges(3,impparbins);
  hSparseNCfromB->SetBinEdges(4,massHypoBins); 
  flistNoCutsFromB->Add(hSparseNCfromB);



  THnSparseF *hSparseRecoNCfromB=new THnSparseF("hSparseRecoNCfromB","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseRecoNCfromB->SetBinEdges(0,massbins);
  hSparseRecoNCfromB->SetBinEdges(1,massbins);
  hSparseRecoNCfromB->SetBinEdges(2,ptbinsForNsparse);
  hSparseRecoNCfromB->SetBinEdges(3,impparbins);
  hSparseRecoNCfromB->SetBinEdges(4,massHypoBins); 
  flistNoCutsFromB->Add(hSparseRecoNCfromB);


  TH1F *hetaNCfromB;
  TH1F *hCosPDPBNCfromB;
  TH1F *hCosPcPDNCfromB;

   // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0NCfromBpt;
  TH2F *hangletracksVSd0xd0NCfromBpt;
  TH2F *hangletracksVSd0D0NCfromBpt;
  TH1F *hd0xd0NCfromBpt;
  flistNoCutsFromB->Add(hInvMassPtNCfromB);

  TH2F *hTOFpidNCfromB=new TH2F("hTOFpidNCfromB","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistNoCutsFromB->Add(hTOFpidNCfromB);

  for(Int_t i=0;i<fnbins;i++){
    namehist="hd0zD0ptNCfromB_pt";
    namehist+=i;
    titlehist="d0(z) No Cuts FromB ptbin=";
    titlehist+=i;
    hd0zD0ptNCfromB=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptNCfromB->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptNCfromB->SetYTitle("Entries");
    flistNoCutsFromB->Add(hd0zD0ptNCfromB);

    namehist="hInvMassD0NCfromB_pt";
    namehist+=i;
    titlehist="Invariant Mass No Cuts FromB ptbin=";
    titlehist+=i;
    hInvMassD0NCfromB=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0NCfromB->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0NCfromB->SetYTitle("Entries");
    flistNoCutsFromB->Add(hInvMassD0NCfromB);


    namehist="hInvMassD0barNCfromB_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar No Cuts FromB ptbin=";
    titlehist+=i;
    hInvMassD0barNCfromB=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barNCfromB->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barNCfromB->SetYTitle("Entries");
    flistNoCutsFromB->Add(hInvMassD0barNCfromB);



    namehist="hetaNCfromB_pt";
    namehist+=i;
    titlehist="eta No Cuts FromB ptbin=";
    titlehist+=i;
    hetaNCfromB=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaNCfromB->SetXTitle("Pseudorapidity");
    hetaNCfromB->SetYTitle("Entries");
    flistNoCutsFromB->Add(hetaNCfromB);

    namehist="hCosPDPBNCfromB_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBNCfromB=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBNCfromB->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBNCfromB->SetYTitle("Entries");
    flistNoCutsFromB->Add(hCosPDPBNCfromB);

    namehist="hCosPcPDNCfromB_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDNCfromB=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDNCfromB->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDNCfromB->SetYTitle("Entries");
    flistNoCutsFromB->Add(hCosPcPDNCfromB);

// %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0NCfromB_pt";
    namehist+=i;
    titlehist="d0xd0 No Cuts FromB ptbin=";
    titlehist+=i;
    hd0xd0NCfromBpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0NCfromBpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0NCfromBpt->SetYTitle("Entries");
    flistNoCutsFromB->Add(hd0xd0NCfromBpt);


    namehist="hd0D0VSd0xd0NCfromB_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} No Cuts FromB ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0NCfromBpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0NCfromBpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0NCfromBpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistNoCutsFromB->Add(hd0D0VSd0xd0NCfromBpt);
    
    
    namehist="hangletracksVSd0xd0NCfromB_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} No Cuts FromB ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0NCfromBpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0NCfromBpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0NCfromBpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistNoCutsFromB->Add(hangletracksVSd0xd0NCfromBpt);
    

    namehist="hangletracksVSd0D0NCfromB_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} No Cuts FromB ptbin=";
    titlehist+=i;
    hangletracksVSd0D0NCfromBpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0NCfromBpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0NCfromBpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistNoCutsFromB->Add(hangletracksVSd0D0NCfromBpt);

    
  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
  
  TH1F *hd0D0ptNCfromBPM;
  TH1F *hMCd0D0ptNCfromBPM;
  TH1F *hd0D0VtxTrueptNCfromBPM;
  TH1F *hd0D0ptNCfromBSB;
  TH1F *hMCd0D0ptNCfromBSB;
  TH1F *hd0D0VtxTrueptNCfromBSB;
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
    
    hd0D0ptNCfromBPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCfromBPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCfromBPM->SetYTitle("Entries");
    flistNoCutsFromB->Add(hd0D0ptNCfromBPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCfromBPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCfromBPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCfromBPM->SetYTitle("Entries");
    flistNoCutsFromB->Add(hMCd0D0ptNCfromBPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCfromBPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCfromBPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCfromBPM->SetYTitle("Entries");
    flistNoCutsFromB->Add(hd0D0VtxTrueptNCfromBPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCfromBSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCfromBSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCfromBSB->SetYTitle("Entries");
    flistNoCutsFromB->Add(hd0D0ptNCfromBSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCfromBSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCfromBSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCfromBSB->SetYTitle("Entries");
    flistNoCutsFromB->Add(hMCd0D0ptNCfromBSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCfromBSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCfromBSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCfromBSB->SetYTitle("Entries");
    flistNoCutsFromB->Add(hd0D0VtxTrueptNCfromBSB);
  }



  //############ NO CUTS FROM DSTAR HISTOGRAMS ###########
  //
  //#############  global properties histos #######

  TH2F *hCPtaVSd0d0NCfromDstar=new TH2F("hCPtaVSd0d0NCfromDstar","hCPtaVSd0d0_NoCuts_FromDStar",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZNCfromDstar=new TH1F("hSecVtxZNCfromDstar","hSecVtxZ_NoCuts_FromDStar",1000,-8.,8.);
  TH1F *hSecVtxXNCfromDstar=new TH1F("hSecVtxXNCfromDstar","hSecVtxX_NoCuts_FromDStar",1000,-3000.,3000.);
  TH1F *hSecVtxYNCfromDstar=new TH1F("hSecVtxYNCfromDstar","hSecVtxY_NoCuts_FromDStar",1000,-3000.,3000.);
  TH2F *hSecVtxXYNCfromDstar=new TH2F("hSecVtxXYNCfromDstar","hSecVtxXY_NoCuts_FromDStar",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiNCfromDstar=new TH1F("hSecVtxPhiNCfromDstar","hSecVtxPhi_NoCuts_FromDStar",180,-180.1,180.1);
  TH1F *hd0singlTrackNCfromDstar=new TH1F("hd0singlTrackNCfromDstar","hd0singlTrackNoCuts_fromDstar",1000,-5000.,5000.);
  TH1F *hCPtaNCfromDstar=new TH1F("hCPtaNCfromDstar","hCPta_NoCuts_FromDStar",100,-1.,1.);
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
  flistNoCutsFromDstar->Add(hd0singlTrackNCfromDstar);
  flistNoCutsFromDstar->Add(hCPtaNCfromDstar);
  flistNoCutsFromDstar->Add(hd0xd0NCfromDstar);
  flistNoCutsFromDstar->Add(hMassTrueNCfromDstar);
  flistNoCutsFromDstar->Add(hMassNCfromDstar);
  flistNoCutsFromDstar->Add(hMassTrueNCfromDstarPM);
  flistNoCutsFromDstar->Add(hMassNCfromDstarPM);
  flistNoCutsFromDstar->Add(hMassTrueNCfromDstarSB);
  flistNoCutsFromDstar->Add(hMassNCfromDstarSB);




//%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaNCfromDstar=new TH1F("hdcaNCfromDstar","hdca_NoCuts_FromDstar",100,0.,1000.);
  hdcaNCfromDstar->SetXTitle("dca   [#mum]");
  hdcaNCfromDstar->SetYTitle("Entries");
  TH1F *hcosthetastarNCfromDstar=new TH1F("hcosthetastarNCfromDstar","hCosThetaStar_NoCuts_FromDstar",50,-1.,1.);
  hcosthetastarNCfromDstar->SetXTitle("cos #theta^{*}");
  hcosthetastarNCfromDstar->SetYTitle("Entries");
  TH1F *hptD0NCfromDstar=new TH1F("hptD0NCfromDstar","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0NCfromDstar->SetXTitle("p_{t}  [GeV/c]");
  hptD0NCfromDstar->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtNCfromDstar=new TH1F("hptD0VsMaxPtNCfromDstar","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtNCfromDstar=new TH2F("hptD0PTallsqrtNCfromDstar","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallNCfromDstar=new TH2F("hptD0PTallNCfromDstar","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBNCfromDstar=new TH2F("hptD0vsptBNCfromDstar","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBNCfromDstar=new TH2F("hpD0vspBNCfromDstar","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkNCfromDstar=new TH2F("hptD0vsptcquarkNCfromDstar","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkNCfromDstar=new TH2F("hpD0vspcquarkNCfromDstar","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistNoCutsFromDstar->Add(hdcaNCfromDstar);
  flistNoCutsFromDstar->Add(hcosthetastarNCfromDstar);
  flistNoCutsFromDstar->Add(hptD0NCfromDstar);
  flistNoCutsFromDstar->Add(hptD0VsMaxPtNCfromDstar);
  flistNoCutsFromDstar->Add(hptD0PTallsqrtNCfromDstar);
  flistNoCutsFromDstar->Add(hptD0PTallNCfromDstar);
  flistNoCutsFromDstar->Add(hptD0vsptBNCfromDstar);
  flistNoCutsFromDstar->Add(hpD0vspBNCfromDstar);
  flistNoCutsFromDstar->Add(hptD0vsptcquarkNCfromDstar);
  flistNoCutsFromDstar->Add(hpD0vspcquarkNCfromDstar);
 
  TH1F *hd0zD0ptNCfromDstar;
  TH1F *hInvMassD0NCfromDstar,*hInvMassD0barNCfromDstar;
  TH2F *hInvMassPtNCfromDstar=new TH2F("hInvMassPtNCfromDstar","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
 THnSparseF *hSparseNCfromDstar=new THnSparseF("hSparseNCfromDstar","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseNCfromDstar->SetBinEdges(0,massbins);
  hSparseNCfromDstar->SetBinEdges(1,massbins);
  hSparseNCfromDstar->SetBinEdges(2,ptbinsForNsparse);
  hSparseNCfromDstar->SetBinEdges(3,impparbins);
  hSparseNCfromDstar->SetBinEdges(4,massHypoBins); 
  flistNoCutsFromDstar->Add(hSparseNCfromDstar);
  TH1F *hetaNCfromDstar;
  TH1F *hCosPDPBNCfromDstar;
  TH1F *hCosPcPDNCfromDstar;
  flistNoCutsFromDstar->Add(hInvMassPtNCfromDstar);
  // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0NCfromDstarpt;
  TH2F *hangletracksVSd0xd0NCfromDstarpt;
  TH2F *hangletracksVSd0D0NCfromDstarpt;
  TH1F *hd0xd0NCfromDstarpt;

  TH2F *hTOFpidNCfromDstar=new TH2F("hTOFpidNCfromDstar","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistNoCutsFromDstar->Add(hTOFpidNCfromDstar);

  for(Int_t i=0;i<fnbins;i++){
    namehist="hd0zD0ptNCfromDstar_pt";
    namehist+=i;
    titlehist="d0(z) No Cuts FromDstarm ptbin=";
    titlehist+=i;
    hd0zD0ptNCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptNCfromDstar->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptNCfromDstar->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hd0zD0ptNCfromDstar);

    namehist="hInvMassD0NCfromDstar_pt";
    namehist+=i;
    titlehist="Invariant Mass No Cuts FromDstar ptbin=";
    titlehist+=i;
    hInvMassD0NCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0NCfromDstar->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0NCfromDstar->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hInvMassD0NCfromDstar);


   namehist="hInvMassD0barNCfromDstar_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar No Cuts FromDstar ptbin=";
    titlehist+=i;
    hInvMassD0barNCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barNCfromDstar->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barNCfromDstar->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hInvMassD0barNCfromDstar);



    namehist="hetaNCfromDstar_pt";
    namehist+=i;
    titlehist="eta No Cuts FromDstar ptbin=";
    titlehist+=i;
    hetaNCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaNCfromDstar->SetXTitle("Pseudorapidity");
    hetaNCfromDstar->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hetaNCfromDstar);

    namehist="hCosPDPBNCfromDstar_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBNCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBNCfromDstar->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBNCfromDstar->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hCosPDPBNCfromDstar);

    namehist="hCosPcPDNCfromDstar_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDNCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDNCfromDstar->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDNCfromDstar->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hCosPcPDNCfromDstar);
  
 // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0NCfromDstar_pt";
    namehist+=i;
    titlehist="d0xd0 No Cuts FromDstar ptbin=";
    titlehist+=i;
    hd0xd0NCfromDstarpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0NCfromDstarpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0NCfromDstarpt->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hd0xd0NCfromDstarpt);


    namehist="hd0D0VSd0xd0NCfromDstar_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} No Cuts FromDstar ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0NCfromDstarpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0NCfromDstarpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0NCfromDstarpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistNoCutsFromDstar->Add(hd0D0VSd0xd0NCfromDstarpt);
    
    
    namehist="hangletracksVSd0xd0NCfromDstar_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} No Cuts FromDstar ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0NCfromDstarpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0NCfromDstarpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0NCfromDstarpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistNoCutsFromDstar->Add(hangletracksVSd0xd0NCfromDstarpt);
    

    namehist="hangletracksVSd0D0NCfromDstar_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} No Cuts FromDstar ptbin=";
    titlehist+=i;
    hangletracksVSd0D0NCfromDstarpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0NCfromDstarpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0NCfromDstarpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistNoCutsFromDstar->Add(hangletracksVSd0D0NCfromDstarpt);
  
  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
  
  TH1F *hd0D0ptNCfromDstPM;
  TH1F *hMCd0D0ptNCfromDstPM;
  TH1F *hd0D0VtxTrueptNCfromDstPM;
  TH1F *hd0D0ptNCfromDstSB;
  TH1F *hMCd0D0ptNCfromDstSB;
  TH1F *hd0D0VtxTrueptNCfromDstSB;
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
    
    hd0D0ptNCfromDstPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCfromDstPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCfromDstPM->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hd0D0ptNCfromDstPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCfromDstPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCfromDstPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCfromDstPM->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hMCd0D0ptNCfromDstPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCfromDstPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCfromDstPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCfromDstPM->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hd0D0VtxTrueptNCfromDstPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCfromDstSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCfromDstSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCfromDstSB->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hd0D0ptNCfromDstSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCfromDstSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCfromDstSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCfromDstSB->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hMCd0D0ptNCfromDstSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCfromDstSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCfromDstSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCfromDstSB->SetYTitle("Entries");
    flistNoCutsFromDstar->Add(hd0D0VtxTrueptNCfromDstSB);
  }


  //############ NO CUTS OTHER HISTOGRAMS ###########
  //
  //########### global properties histos ###########

  TH2F *hCPtaVSd0d0NCother=new TH2F("hCPtaVSd0d0NCother","hCPtaVSd0d0_NoCuts_other",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZNCother=new TH1F("hSecVtxZNCother","hSecVtxZ_NoCuts_other",1000,-8.,8.);
  TH1F *hSecVtxXNCother=new TH1F("hSecVtxXNCother","hSecVtxX_NoCuts_other",1000,-3000.,3000.);
  TH1F *hSecVtxYNCother=new TH1F("hSecVtxYNCother","hSecVtxY_NoCuts_other",1000,-3000.,3000.);
  TH2F *hSecVtxXYNCother=new TH2F("hSecVtxXYNCother","hSecVtxXY_NoCuts_other",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiNCother=new TH1F("hSecVtxPhiNCother","hSecVtxPhi_NoCuts_other",180,-180.1,180.1);
  TH1F *hd0singlTrackNCother=new TH1F("hd0singlTrackNCother","hd0singlTrackNoCuts_Other",1000,-5000.,5000.);
  TH1F *hCPtaNCother=new TH1F("hCPtaNCother","hCPta_NoCuts_other",100,-1.,1.);
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
  flistNoCutsOther->Add(hd0singlTrackNCother);
  flistNoCutsOther->Add(hCPtaNCother);
  flistNoCutsOther->Add(hd0xd0NCother);
  flistNoCutsOther->Add(hMassTrueNCother);
  flistNoCutsOther->Add(hMassNCother);
  flistNoCutsOther->Add(hMassTrueNCotherPM);
  flistNoCutsOther->Add(hMassNCotherPM);
  flistNoCutsOther->Add(hMassTrueNCotherSB);
  flistNoCutsOther->Add(hMassNCotherSB);



 //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaNCother=new TH1F("hdcaNCother","hdca_NoCuts_Other",100,0.,1000.);
  hdcaNCother->SetXTitle("dca   [#mum]");
  hdcaNCother->SetYTitle("Entries");
  TH1F *hcosthetastarNCother=new TH1F("hcosthetastarNCother","hCosThetaStar_NoCuts_Other",50,-1.,1.);
  hcosthetastarNCother->SetXTitle("cos #theta^{*}");
  hcosthetastarNCother->SetYTitle("Entries");
  TH1F *hptD0NCother=new TH1F("hptD0NCother","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0NCother->SetXTitle("p_{t}  [GeV/c]");
  hptD0NCother->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtNCother=new TH1F("hptD0VsMaxPtNCother","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtNCother=new TH2F("hptD0PTallsqrtNCother","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallNCother=new TH2F("hptD0PTallNCother","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBNCother=new TH2F("hptD0vsptBNCother","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBNCother=new TH2F("hpD0vspBNCother","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkNCother=new TH2F("hptD0vsptcquarkNCother","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkNCother=new TH2F("hpD0vspcquarkNCother","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistNoCutsOther->Add(hdcaNCother);
  flistNoCutsOther->Add(hcosthetastarNCother);
  flistNoCutsOther->Add(hptD0NCother);
  flistNoCutsOther->Add(hptD0VsMaxPtNCother);
  flistNoCutsOther->Add(hptD0PTallsqrtNCother);
  flistNoCutsOther->Add(hptD0PTallNCother);
  flistNoCutsOther->Add(hptD0vsptBNCother);
  flistNoCutsOther->Add(hpD0vspBNCother);
  flistNoCutsOther->Add(hptD0vsptcquarkNCother);
  flistNoCutsOther->Add(hpD0vspcquarkNCother);

  TH1F *hd0zD0ptNCother;
  TH1F *hInvMassD0NCother,*hInvMassD0barNCother;
  TH2F *hInvMassPtNCother=new TH2F("hInvMassPtNCother","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
 THnSparseF *hSparseNCother=new THnSparseF("hSparseNCother","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseNCother->SetBinEdges(0,massbins);
  hSparseNCother->SetBinEdges(1,massbins);
  hSparseNCother->SetBinEdges(2,ptbinsForNsparse);
  hSparseNCother->SetBinEdges(3,impparbins);
  hSparseNCother->SetBinEdges(4,massHypoBins); 
  flistNoCutsOther->Add(hSparseNCother);
  TH1F *hetaNCother;
  TH1F *hCosPDPBNCother;
  TH1F *hCosPcPDNCother;
  flistNoCutsOther->Add(hInvMassPtNCother);
  // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0NCotherpt;
  TH2F *hangletracksVSd0xd0NCotherpt;
  TH2F *hangletracksVSd0D0NCotherpt;
  TH1F *hd0xd0NCotherpt;

  TH2F *hTOFpidNCother=new TH2F("hTOFpidNCother","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistNoCutsOther->Add(hTOFpidNCother);

  for(Int_t i=0;i<fnbins;i++){
    namehist="hd0zD0ptNCother_pt";
    namehist+=i;
    titlehist="d0(z) No Cuts Otherm ptbin=";
    titlehist+=i;
    hd0zD0ptNCother=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptNCother->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptNCother->SetYTitle("Entries");
    flistNoCutsOther->Add(hd0zD0ptNCother);

    namehist="hInvMassD0NCother_pt";
    namehist+=i;
    titlehist="Invariant Mass No Cuts Other ptbin=";
    titlehist+=i;
    hInvMassD0NCother=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0NCother->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0NCother->SetYTitle("Entries");
    flistNoCutsOther->Add(hInvMassD0NCother);


   namehist="hInvMassD0barNCother_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar No Cuts Other ptbin=";
    titlehist+=i;
    hInvMassD0barNCother=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barNCother->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barNCother->SetYTitle("Entries");
    flistNoCutsOther->Add(hInvMassD0barNCother);


    namehist="hetaNCother_pt";
    namehist+=i;
    titlehist="eta No Cuts Other ptbin=";
    titlehist+=i;
    hetaNCother=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaNCother->SetXTitle("Pseudorapidity");
    hetaNCother->SetYTitle("Entries");
    flistNoCutsOther->Add(hetaNCother);

    namehist="hCosPDPBNCother_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBNCother=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBNCother->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBNCother->SetYTitle("Entries");
    flistNoCutsOther->Add(hCosPDPBNCother);

    namehist="hCosPcPDNCother_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDNCother=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDNCother->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDNCother->SetYTitle("Entries");
    flistNoCutsOther->Add(hCosPcPDNCother);
    

 // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0NCother_pt";
    namehist+=i;
    titlehist="d0xd0 No Cuts Other ptbin=";
    titlehist+=i;
    hd0xd0NCotherpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0NCotherpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0NCotherpt->SetYTitle("Entries");
    flistNoCutsOther->Add(hd0xd0NCotherpt);


    namehist="hd0D0VSd0xd0NCother_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} No Cuts Other ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0NCotherpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0NCotherpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0NCotherpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistNoCutsOther->Add(hd0D0VSd0xd0NCotherpt);
    
    
    namehist="hangletracksVSd0xd0NCother_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} No Cuts Other ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0NCotherpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0NCotherpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0NCotherpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistNoCutsOther->Add(hangletracksVSd0xd0NCotherpt);
    

    namehist="hangletracksVSd0D0NCother_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} No Cuts Other ptbin=";
    titlehist+=i;
    hangletracksVSd0D0NCotherpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0NCotherpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0NCotherpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistNoCutsOther->Add(hangletracksVSd0D0NCotherpt);

  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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
  
  TH1F *hd0D0ptNCotherPM;
  TH1F *hMCd0D0ptNCotherPM;
  TH1F *hd0D0VtxTrueptNCotherPM;
  TH1F *hd0D0ptNCotherSB;
  TH1F *hMCd0D0ptNCotherSB;
  TH1F *hd0D0VtxTrueptNCotherSB;
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
    
    hd0D0ptNCotherPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCotherPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCotherPM->SetYTitle("Entries");
    flistNoCutsOther->Add(hd0D0ptNCotherPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCotherPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCotherPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCotherPM->SetYTitle("Entries");
    flistNoCutsOther->Add(hMCd0D0ptNCotherPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCotherPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCotherPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCotherPM->SetYTitle("Entries");
    flistNoCutsOther->Add(hd0D0VtxTrueptNCotherPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptNCotherSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptNCotherSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptNCotherSB->SetYTitle("Entries");
    flistNoCutsOther->Add(hd0D0ptNCotherSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptNCotherSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptNCotherSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptNCotherSB->SetYTitle("Entries");
    flistNoCutsOther->Add(hMCd0D0ptNCotherSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptNCotherSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptNCotherSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptNCotherSB->SetYTitle("Entries");
    flistNoCutsOther->Add(hd0D0VtxTrueptNCotherSB);
  }


  //################################################################################################
  //                                                                                               #
  //                         HISTOS FOR LOOSE CUTS                                                 #
  //                                                                                               #
  //################################################################################################

  //############ LOOSE CUTS SIGNAL HISTOGRAMS ###############
  //
  // ####### global properties histo ############

  TH2F *hCPtaVSd0d0LSCsign=new TH2F("hCPtaVSd0d0LSCsign","hCPtaVSd0d0_LooseCuts_Signal",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZLSCsign=new TH1F("hSecVtxZLSCsign","hSecVtxZ_LooseCuts_Signal",1000,-8.,8.);
  TH1F *hSecVtxXLSCsign=new TH1F("hSecVtxXLSCsign","hSecVtxX_LooseCuts_Signal",1000,-3000.,3000.);
  TH1F *hSecVtxYLSCsign=new TH1F("hSecVtxYLSCsign","hSecVtxY_LooseCuts_Signal",1000,-3000.,3000.);
  TH2F *hSecVtxXYLSCsign=new TH2F("hSecVtxXYLSCsign","hSecVtxXY_LooseCuts_Signal",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiLSCsign=new TH1F("hSecVtxPhiLSCsign","hSecVtxPhi_LooseCuts_Signal",180,-180.1,180.1);
  TH1F *hd0singlTrackLSCsign=new TH1F("hd0singlTrackLSCsign","hd0singlTrackLooseCuts_Signal",1000,-5000.,5000.);
  TH1F *hCPtaLSCsign=new TH1F("hCPtaLSCsign","hCPta_LooseCuts_Signal",100,-1.,1.);
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
  flistLsCutsSignal->Add(hd0singlTrackLSCsign);
  flistLsCutsSignal->Add(hCPtaLSCsign);
  flistLsCutsSignal->Add(hd0xd0LSCsign);
  flistLsCutsSignal->Add(hMassTrueLSCsign);
  flistLsCutsSignal->Add(hMassLSCsign);
  flistLsCutsSignal->Add(hMassTrueLSCsignPM);
  flistLsCutsSignal->Add(hMassLSCsignPM);
  flistLsCutsSignal->Add(hMassTrueLSCsignSB);
  flistLsCutsSignal->Add(hMassLSCsignSB);


  //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaLSCsign=new TH1F("hdcaLSCsign","hdca_LooseCuts_Sign",100,0.,1000.);
  hdcaLSCsign->SetXTitle("dca   [#mum]");
  hdcaLSCsign->SetYTitle("Entries");
  TH1F *hcosthetastarLSCsign=new TH1F("hcosthetastarLSCsign","hCosThetaStar_LooseCuts_Sign",50,-1.,1.);
  hcosthetastarLSCsign->SetXTitle("cos #theta^{*}");
  hcosthetastarLSCsign->SetYTitle("Entries");
  TH1F *hptD0LSCsign=new TH1F("hptD0LSCsign","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0LSCsign->SetXTitle("p_{t}  [GeV/c]");
  hptD0LSCsign->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtLSCsign=new TH1F("hptD0VsMaxPtLSCsign","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtLSCsign=new TH2F("hptD0PTallsqrtLSCsign","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallLSCsign=new TH2F("hptD0PTallLSCsign","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBLSCsign=new TH2F("hptD0vsptBLSCsign","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBLSCsign=new TH2F("hpD0vspBLSCsign","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkLSCsign=new TH2F("hptD0vsptcquarkLSCsign","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkLSCsign=new TH2F("hpD0vspcquarkLSCsign","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistLsCutsSignal->Add(hdcaLSCsign);
  flistLsCutsSignal->Add(hcosthetastarLSCsign);
  flistLsCutsSignal->Add(hptD0LSCsign);
  flistLsCutsSignal->Add(hptD0VsMaxPtLSCsign);
  flistLsCutsSignal->Add(hptD0PTallsqrtLSCsign);
  flistLsCutsSignal->Add(hptD0PTallLSCsign);
  flistLsCutsSignal->Add(hptD0vsptBLSCsign);
  flistLsCutsSignal->Add(hpD0vspBLSCsign);
  flistLsCutsSignal->Add(hptD0vsptcquarkLSCsign);
  flistLsCutsSignal->Add(hpD0vspcquarkLSCsign);
 
  TH1F *hd0zD0ptLSCsign;
  TH1F *hInvMassD0LSCsign,*hInvMassD0barLSCsign;
  TH2F *hInvMassPtLSCsign=new TH2F("hInvMassPtLSCsign","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
 THnSparseF *hSparseLSCsign=new THnSparseF("hSparseLSCsign","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseLSCsign->SetBinEdges(0,massbins);
  hSparseLSCsign->SetBinEdges(1,massbins);
  hSparseLSCsign->SetBinEdges(2,ptbinsForNsparse);
  hSparseLSCsign->SetBinEdges(3,impparbins);
  hSparseLSCsign->SetBinEdges(4,massHypoBins); 
  flistLsCutsSignal->Add(hSparseLSCsign);
  TH1F *hetaLSCsign;
  TH1F *hCosPDPBLSCsign;
  TH1F *hCosPcPDLSCsign;
  flistLsCutsSignal->Add(hInvMassPtLSCsign);


  
 THnSparseF *hSparseCxyLxyLSCsign=new THnSparseF("hSparseCxyLxyLSCsign","Candidate Mass;massD0;Pt;CosXY;Lxy",4,nbinsSparsCxyLxy,binLowLimitSparseCxyLxy,binUpLimitSparseCxyLxy);
  hSparseCxyLxyLSCsign->SetBinEdges(1,ptbinlimitsCxyLxy); 
  hSparseCxyLxyLSCsign->GetAxis(0)->SetName("mass");
  hSparseCxyLxyLSCsign->GetAxis(0)->SetTitle("Invariant Mass (K#pi) [GeV/c^{2}]");
  hSparseCxyLxyLSCsign->GetAxis(1)->SetName("pt");
  hSparseCxyLxyLSCsign->GetAxis(1)->SetTitle("p_{t} [GeV/c]");
  hSparseCxyLxyLSCsign->GetAxis(2)->SetName("CosPointXY");
  hSparseCxyLxyLSCsign->GetAxis(2)->SetTitle("Cos#theta_{point}^{XY}");
  hSparseCxyLxyLSCsign->GetAxis(3)->SetName("NormDecLengthXY");
  hSparseCxyLxyLSCsign->GetAxis(3)->SetTitle("Normalized XY decay length");

  flistLsCutsSignal->Add(hSparseCxyLxyLSCsign);
  // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0LSCsignpt;
  TH2F *hangletracksVSd0xd0LSCsignpt;
  TH2F *hangletracksVSd0D0LSCsignpt;
  TH1F *hd0xd0LSCsignpt;
  TH1F *hPhiHistPMLSCsignpt,*hPhiHistSBLSCsignpt;

  TH2F *hTOFpidLSCsign=new TH2F("hTOFpidLSCsign","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistLsCutsSignal->Add(hTOFpidLSCsign);

  for(Int_t i=0;i<fnbins;i++){
 
    namehist="hPhiHistPMLSCsign_pt";
    namehist+=i;
    titlehist="Azimuthal correlation LS Cuts Sign PM ptbin=";
    titlehist+=i;
    hPhiHistPMLSCsignpt=new TH1F(namehist.Data(),titlehist.Data(),100,-3.15,3.15);
    hPhiHistPMLSCsignpt->Sumw2();
    flistLsCutsSignal->Add(hPhiHistPMLSCsignpt);

    namehist="hPhiHistSBLSCsign_pt";
    namehist+=i;
    titlehist="Azimuthal correlation LS Cuts Sign SB ptbin=";
    titlehist+=i;
    hPhiHistSBLSCsignpt=new TH1F(namehist.Data(),titlehist.Data(),100,-3.15,3.15);
    hPhiHistSBLSCsignpt->Sumw2();
    flistLsCutsSignal->Add(hPhiHistSBLSCsignpt);

  

    namehist="hd0zD0ptLSCsign_pt";
    namehist+=i;
    titlehist="d0(z) Loose Cuts Signm ptbin=";
    titlehist+=i;
    hd0zD0ptLSCsign=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptLSCsign->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptLSCsign->SetYTitle("Entries");
    flistLsCutsSignal->Add(hd0zD0ptLSCsign);

    namehist="hInvMassD0LSCsign_pt";
    namehist+=i;
    titlehist="Invariant Mass Loose Cuts Sign ptbin=";
    titlehist+=i;
    hInvMassD0LSCsign=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0LSCsign->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0LSCsign->SetYTitle("Entries");
    flistLsCutsSignal->Add(hInvMassD0LSCsign);


    namehist="hInvMassD0barLSCsign_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar Loose Cuts Signal ptbin=";
    titlehist+=i;
    hInvMassD0barLSCsign=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barLSCsign->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barLSCsign->SetYTitle("Entries");
    flistLsCutsSignal->Add(hInvMassD0barLSCsign);

    namehist="hetaLSCsign_pt";
    namehist+=i;
    titlehist="eta Loose Cuts Sign ptbin=";
    titlehist+=i;
    hetaLSCsign=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaLSCsign->SetXTitle("Pseudorapidity");
    hetaLSCsign->SetYTitle("Entries");
    flistLsCutsSignal->Add(hetaLSCsign);

    namehist="hCosPDPBLSCsign_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBLSCsign=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBLSCsign->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBLSCsign->SetYTitle("Entries");
    flistLsCutsSignal->Add(hCosPDPBLSCsign);

    namehist="hCosPcPDLSCsign_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDLSCsign=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDLSCsign->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDLSCsign->SetYTitle("Entries");
    flistLsCutsSignal->Add(hCosPcPDLSCsign);


 // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0LSCsign_pt";
    namehist+=i;
    titlehist="d0xd0 Loose Cuts Sign ptbin=";
    titlehist+=i;
    hd0xd0LSCsignpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0LSCsignpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0LSCsignpt->SetYTitle("Entries");
    flistLsCutsSignal->Add(hd0xd0LSCsignpt);


    namehist="hd0D0VSd0xd0LSCsign_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} Loose Cuts Sign ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0LSCsignpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0LSCsignpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0LSCsignpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistLsCutsSignal->Add(hd0D0VSd0xd0LSCsignpt);
    
    
    namehist="hangletracksVSd0xd0LSCsign_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} Loose Cuts Sign ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0LSCsignpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0LSCsignpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0LSCsignpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistLsCutsSignal->Add(hangletracksVSd0xd0LSCsignpt);
    

    namehist="hangletracksVSd0D0LSCsign_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} Loose Cuts Sign ptbin=";
    titlehist+=i;
    hangletracksVSd0D0LSCsignpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0LSCsignpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0LSCsignpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistLsCutsSignal->Add(hangletracksVSd0D0LSCsignpt);

    
  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
  
  TH1F *hd0D0ptLSCsignPM;
  TH1F *hMCd0D0ptLSCsignPM;
  TH1F *hd0D0VtxTrueptLSCsignPM;
  TH1F *hd0D0ptLSCsignSB;
  TH1F *hMCd0D0ptLSCsignSB;
  TH1F *hd0D0VtxTrueptLSCsignSB;
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
    
    hd0D0ptLSCsignPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCsignPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCsignPM->SetYTitle("Entries");
    flistLsCutsSignal->Add(hd0D0ptLSCsignPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCsignPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCsignPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCsignPM->SetYTitle("Entries");
    flistLsCutsSignal->Add(hMCd0D0ptLSCsignPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCsignPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCsignPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCsignPM->SetYTitle("Entries");
    flistLsCutsSignal->Add(hd0D0VtxTrueptLSCsignPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCsignSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCsignSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCsignSB->SetYTitle("Entries");
    flistLsCutsSignal->Add(hd0D0ptLSCsignSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCsignSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCsignSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCsignSB->SetYTitle("Entries");
    flistLsCutsSignal->Add(hMCd0D0ptLSCsignSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCsignSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCsignSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCsignSB->SetYTitle("Entries");
    flistLsCutsSignal->Add(hd0D0VtxTrueptLSCsignSB);
  }


  //############ LOOSE CUTS BACKGROUND HISTOGRAMS ###########
  //
  //   ######## global properties histos #######
  TH2F *hCPtaVSd0d0LSCback=new TH2F("hCPtaVSd0d0LSCback","hCPtaVSd0d0_LooseCuts_Background",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZLSCback=new TH1F("hSecVtxZLSCback","hSecVtxZ_LooseCuts_Background",1000,-8.,8.);
  TH1F *hSecVtxXLSCback=new TH1F("hSecVtxXLSCback","hSecVtxX_LooseCuts_Background",1000,-3000.,3000.);
  TH1F *hSecVtxYLSCback=new TH1F("hSecVtxYLSCback","hSecVtxY_LooseCuts_Background",1000,-3000.,3000.);
  TH2F *hSecVtxXYLSCback=new TH2F("hSecVtxXYLSCback","hSecVtxXY_LooseCuts_Background",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiLSCback=new TH1F("hSecVtxPhiLSCback","hSecVtxPhi_LooseCuts_Background",180,-180.1,180.1);
  TH1F *hd0singlTrackLSCback=new TH1F("hd0singlTrackLSCback","hd0singlTrackLooseCuts_Back",1000,-5000.,5000.);
  TH1F *hCPtaLSCback=new TH1F("hCPtaLSCback","hCPta_LooseCuts_Background",100,-1.,1.);
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
  flistLsCutsBack->Add(hd0singlTrackLSCback);
  flistLsCutsBack->Add(hCPtaLSCback);
  flistLsCutsBack->Add(hd0xd0LSCback);
  flistLsCutsBack->Add(hMassTrueLSCback);
  flistLsCutsBack->Add(hMassLSCback);
  flistLsCutsBack->Add(hMassTrueLSCbackPM);
  flistLsCutsBack->Add(hMassLSCbackPM);
  flistLsCutsBack->Add(hMassTrueLSCbackSB);
  flistLsCutsBack->Add(hMassLSCbackSB);








 //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaLSCback=new TH1F("hdcaLSCback","hdca_LooseCuts_Backgr",100,0.,1000.);
  hdcaLSCback->SetXTitle("dca   [#mum]");
  hdcaLSCback->SetYTitle("Entries");
  TH1F *hcosthetastarLSCback=new TH1F("hcosthetastarLSCback","hCosThetaStar_LooseCuts_Backgr",50,-1.,1.);
  hcosthetastarLSCback->SetXTitle("cos #theta^{*}");
  hcosthetastarLSCback->SetYTitle("Entries");
  TH1F *hptD0LSCback=new TH1F("hptD0LSCback","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0LSCback->SetXTitle("p_{t}  [GeV/c]");
  hptD0LSCback->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtLSCback=new TH1F("hptD0VsMaxPtLSCback","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtLSCback=new TH2F("hptD0PTallsqrtLSCback","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallLSCback=new TH2F("hptD0PTallLSCback","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBLSCback=new TH2F("hptD0vsptBLSCback","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBLSCback=new TH2F("hpD0vspBLSCback","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkLSCback=new TH2F("hptD0vsptcquarkLSCback","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkLSCback=new TH2F("hpD0vspcquarkLSCback","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistLsCutsBack->Add(hdcaLSCback);
  flistLsCutsBack->Add(hcosthetastarLSCback);
  flistLsCutsBack->Add(hptD0LSCback);
  flistLsCutsBack->Add(hptD0VsMaxPtLSCback);
  flistLsCutsBack->Add(hptD0PTallsqrtLSCback);
  flistLsCutsBack->Add(hptD0PTallLSCback);
  flistLsCutsBack->Add(hptD0vsptBLSCback);
  flistLsCutsBack->Add(hpD0vspBLSCback);
  flistLsCutsBack->Add(hptD0vsptcquarkLSCback);
  flistLsCutsBack->Add(hpD0vspcquarkLSCback);
 
  TH1F *hd0zD0ptLSCback;
  TH1F *hInvMassD0LSCback,*hInvMassD0barLSCback;
  TH2F *hInvMassPtLSCback=new TH2F("hInvMassPtLSCback","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
 THnSparseF *hSparseLSCback=new THnSparseF("hSparseLSCback","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseLSCback->SetBinEdges(0,massbins);
  hSparseLSCback->SetBinEdges(1,massbins);
  hSparseLSCback->SetBinEdges(2,ptbinsForNsparse);
  hSparseLSCback->SetBinEdges(3,impparbins);
  hSparseLSCback->SetBinEdges(4,massHypoBins); 
  flistLsCutsBack->Add(hSparseLSCback);
  TH1F *hetaLSCback;
  TH1F *hCosPDPBLSCback;
  TH1F *hCosPcPDLSCback;
  flistLsCutsBack->Add(hInvMassPtLSCback);
 // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0LSCbackpt;
  TH2F *hangletracksVSd0xd0LSCbackpt;
  TH2F *hangletracksVSd0D0LSCbackpt;
  TH1F *hd0xd0LSCbackpt;

  TH2F *hTOFpidLSCback=new TH2F("hTOFpidLSCback","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistLsCutsBack->Add(hTOFpidLSCback);

  for(Int_t i=0;i<fnbins;i++){
    namehist="hd0zD0ptLSCback_pt";
    namehist+=i;
    titlehist="d0(z) Loose Cuts Backgr ptbin=";
    titlehist+=i;
    hd0zD0ptLSCback=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptLSCback->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptLSCback->SetYTitle("Entries");
    flistLsCutsBack->Add(hd0zD0ptLSCback);

    namehist="hInvMassD0LSCback_pt";
    namehist+=i;
    titlehist="Invariant Mass Loose Cuts Backgr ptbin=";
    titlehist+=i;
    hInvMassD0LSCback=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0LSCback->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0LSCback->SetYTitle("Entries");
    flistLsCutsBack->Add(hInvMassD0LSCback);
    
    namehist="hInvMassD0barLSCback_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar Loose Cuts Back ptbin=";
    titlehist+=i;
    hInvMassD0barLSCback=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barLSCback->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barLSCback->SetYTitle("Entries");
    flistLsCutsBack->Add(hInvMassD0barLSCback);


    namehist="hetaLSCback_pt";
    namehist+=i;
    titlehist="eta Loose Cuts Backgr ptbin=";
    titlehist+=i;
    hetaLSCback=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaLSCback->SetXTitle("Pseudorapidity");
    hetaLSCback->SetYTitle("Entries");
    flistLsCutsBack->Add(hetaLSCback);

    namehist="hCosPDPBLSCback_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBLSCback=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBLSCback->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBLSCback->SetYTitle("Entries");
    flistLsCutsBack->Add(hCosPDPBLSCback);

    namehist="hCosPcPDLSCback_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDLSCback=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDLSCback->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDLSCback->SetYTitle("Entries");
    flistLsCutsBack->Add(hCosPcPDLSCback);

 // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0LSCback_pt";
    namehist+=i;
    titlehist="d0xd0 Loose Cuts Back ptbin=";
    titlehist+=i;
    hd0xd0LSCbackpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0LSCbackpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0LSCbackpt->SetYTitle("Entries");
    flistLsCutsBack->Add(hd0xd0LSCbackpt);


    namehist="hd0D0VSd0xd0LSCback_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} Loose Cuts Back ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0LSCbackpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0LSCbackpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0LSCbackpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistLsCutsBack->Add(hd0D0VSd0xd0LSCbackpt);
    
    
    namehist="hangletracksVSd0xd0LSCback_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} Loose Cuts Back ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0LSCbackpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0LSCbackpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0LSCbackpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistLsCutsBack->Add(hangletracksVSd0xd0LSCbackpt);
    

    namehist="hangletracksVSd0D0LSCback_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} Loose Cuts Back ptbin=";
    titlehist+=i;
    hangletracksVSd0D0LSCbackpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0LSCbackpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0LSCbackpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistLsCutsBack->Add(hangletracksVSd0D0LSCbackpt);
    
  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







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
  
  TH1F *hd0D0ptLSCbackPM;
  TH1F *hMCd0D0ptLSCbackPM;
  TH1F *hd0D0VtxTrueptLSCbackPM;
  TH1F *hd0D0ptLSCbackSB;
  TH1F *hMCd0D0ptLSCbackSB;
  TH1F *hd0D0VtxTrueptLSCbackSB;
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
    
    hd0D0ptLSCbackPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCbackPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCbackPM->SetYTitle("Entries");
    flistLsCutsBack->Add(hd0D0ptLSCbackPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCbackPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCbackPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCbackPM->SetYTitle("Entries");
    flistLsCutsBack->Add(hMCd0D0ptLSCbackPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCbackPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCbackPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCbackPM->SetYTitle("Entries");
    flistLsCutsBack->Add(hd0D0VtxTrueptLSCbackPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCbackSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCbackSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCbackSB->SetYTitle("Entries");
    flistLsCutsBack->Add(hd0D0ptLSCbackSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCbackSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCbackSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCbackSB->SetYTitle("Entries");
    flistLsCutsBack->Add(hMCd0D0ptLSCbackSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCbackSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCbackSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCbackSB->SetYTitle("Entries");
    flistLsCutsBack->Add(hd0D0VtxTrueptLSCbackSB);
  }



 //############ LOOSE CUTS FROMB HISTOGRAMS ###########
  //
  //#######  global properties histos

  TH2F *hCPtaVSd0d0LSCfromB=new TH2F("hCPtaVSd0d0LSCfromB","hCPtaVSd0d0_LooseCuts_FromB",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZLSCfromB=new TH1F("hSecVtxZLSCfromB","hSecVtxZ_LooseCuts_FromB",1000,-8.,8.);
  TH1F *hSecVtxXLSCfromB=new TH1F("hSecVtxXLSCfromB","hSecVtxX_LooseCuts_FromB",1000,-3000.,3000.);
  TH1F *hSecVtxYLSCfromB=new TH1F("hSecVtxYLSCfromB","hSecVtxY_LooseCuts_FromB",1000,-3000.,3000.);
  TH2F *hSecVtxXYLSCfromB=new TH2F("hSecVtxXYLSCfromB","hSecVtxXY_LooseCuts_FromB",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiLSCfromB=new TH1F("hSecVtxPhiLSCfromB","hSecVtxPhi_LooseCuts_FromB",180,-180.1,180.1);
  TH1F *hd0singlTrackLSCfromB=new TH1F("hd0singlTrackLSCfromB","hd0singlTrackLooseCuts_FromB",1000,-5000.,5000.);
  TH1F *hCPtaLSCfromB=new TH1F("hCPtaLSCfromB","hCPta_LooseCuts_FromB",100,-1.,1.);
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
  flistLsCutsFromB->Add(hd0singlTrackLSCfromB);
  flistLsCutsFromB->Add(hCPtaLSCfromB);
  flistLsCutsFromB->Add(hd0xd0LSCfromB);
  flistLsCutsFromB->Add(hMassTrueLSCfromB);
  flistLsCutsFromB->Add(hMassLSCfromB);
  flistLsCutsFromB->Add(hMassTrueLSCfromBPM);
  flistLsCutsFromB->Add(hMassLSCfromBPM);
  flistLsCutsFromB->Add(hMassTrueLSCfromBSB);
  flistLsCutsFromB->Add(hMassLSCfromBSB);




 //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaLSCfromB=new TH1F("hdcaLSCfromB","hdca_LooseCuts_FromB",100,0.,1000.);
  hdcaLSCfromB->SetXTitle("dca   [#mum]");
  hdcaLSCfromB->SetYTitle("Entries");
  TH1F *hcosthetastarLSCfromB=new TH1F("hcosthetastarLSCfromB","hCosThetaStar_LooseCuts_FromB",50,-1.,1.);
  hcosthetastarLSCfromB->SetXTitle("cos #theta^{*}");
  hcosthetastarLSCfromB->SetYTitle("Entries");
  TH1F *hptD0LSCfromB=new TH1F("hptD0LSCfromB","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0LSCfromB->SetXTitle("p_{t}  [GeV/c]");
  hptD0LSCfromB->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtLSCfromB=new TH1F("hptD0VsMaxPtLSCfromB","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtLSCfromB=new TH2F("hptD0PTallsqrtLSCfromB","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallLSCfromB=new TH2F("hptD0PTallLSCfromB","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBLSCfromB=new TH2F("hptD0vsptBLSCfromB","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBLSCfromB=new TH2F("hpD0vspBLSCfromB","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkLSCfromB=new TH2F("hptD0vsptcquarkLSCfromB","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkLSCfromB=new TH2F("hpD0vspcquarkLSCfromB","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistLsCutsFromB->Add(hdcaLSCfromB);
  flistLsCutsFromB->Add(hcosthetastarLSCfromB);
  flistLsCutsFromB->Add(hptD0LSCfromB);
  flistLsCutsFromB->Add(hptD0VsMaxPtLSCfromB);
  flistLsCutsFromB->Add(hptD0PTallsqrtLSCfromB);
  flistLsCutsFromB->Add(hptD0PTallLSCfromB);
  flistLsCutsFromB->Add(hptD0vsptBLSCfromB);
  flistLsCutsFromB->Add(hpD0vspBLSCfromB);
  flistLsCutsFromB->Add(hptD0vsptcquarkLSCfromB);
  flistLsCutsFromB->Add(hpD0vspcquarkLSCfromB);
 
  TH1F *hd0zD0ptLSCfromB;
  TH1F *hInvMassD0LSCfromB,*hInvMassD0barLSCfromB;
  TH2F *hInvMassPtLSCfromB=new TH2F("hInvMassPtLSCfromB","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
 THnSparseF *hSparseLSCfromB=new THnSparseF("hSparseLSCfromB","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseLSCfromB->SetBinEdges(0,massbins);
  hSparseLSCfromB->SetBinEdges(1,massbins);
  hSparseLSCfromB->SetBinEdges(2,ptbinsForNsparse);
  hSparseLSCfromB->SetBinEdges(3,impparbins);
  hSparseLSCfromB->SetBinEdges(4,massHypoBins); 
  flistLsCutsFromB->Add(hSparseLSCfromB);


  THnSparseF *hSparseRecoLSCfromB=new THnSparseF("hSparseRecoLSCfromB","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseRecoLSCfromB->SetBinEdges(0,massbins);
  hSparseRecoLSCfromB->SetBinEdges(1,massbins);
  hSparseRecoLSCfromB->SetBinEdges(2,ptbinsForNsparse);
  hSparseRecoLSCfromB->SetBinEdges(3,impparbins);
  hSparseRecoLSCfromB->SetBinEdges(4,massHypoBins); 
  flistLsCutsFromB->Add(hSparseRecoLSCfromB);


  TH1F *hetaLSCfromB;
  TH1F *hCosPDPBLSCfromB;
  TH1F *hCosPcPDLSCfromB;
  flistLsCutsFromB->Add(hInvMassPtLSCfromB);
   // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0LSCfromBpt;
  TH2F *hangletracksVSd0xd0LSCfromBpt;
  TH2F *hangletracksVSd0D0LSCfromBpt;
  TH1F *hd0xd0LSCfromBpt;


  TH2F *hTOFpidLSCfromB=new TH2F("hTOFpidLSCfromB","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistLsCutsFromB->Add(hTOFpidLSCfromB);

  for(Int_t i=0;i<fnbins;i++){
    namehist="hd0zD0ptLSCfromB_pt";
    namehist+=i;
    titlehist="d0(z) Loose Cuts FromBm ptbin=";
    titlehist+=i;
    hd0zD0ptLSCfromB=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptLSCfromB->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptLSCfromB->SetYTitle("Entries");
    flistLsCutsFromB->Add(hd0zD0ptLSCfromB);

    namehist="hInvMassD0LSCfromB_pt";
    namehist+=i;
    titlehist="Invariant Mass Loose Cuts FromB ptbin=";
    titlehist+=i;
    hInvMassD0LSCfromB=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0LSCfromB->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0LSCfromB->SetYTitle("Entries");
    flistLsCutsFromB->Add(hInvMassD0LSCfromB);

    namehist="hInvMassD0barLSCfromB_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar Loose Cuts FromB ptbin=";
    titlehist+=i;
    hInvMassD0barLSCfromB=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barLSCfromB->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barLSCfromB->SetYTitle("Entries");
    flistLsCutsFromB->Add(hInvMassD0barLSCfromB);

    namehist="hetaLSCfromB_pt";
    namehist+=i;
    titlehist="eta Loose Cuts FromB ptbin=";
    titlehist+=i;
    hetaLSCfromB=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaLSCfromB->SetXTitle("Pseudorapidity");
    hetaLSCfromB->SetYTitle("Entries");
    flistLsCutsFromB->Add(hetaLSCfromB);

    namehist="hCosPDPBLSCfromB_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBLSCfromB=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBLSCfromB->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBLSCfromB->SetYTitle("Entries");
    flistLsCutsFromB->Add(hCosPDPBLSCfromB);

    namehist="hCosPcPDLSCfromB_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDLSCfromB=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDLSCfromB->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDLSCfromB->SetYTitle("Entries");
    flistLsCutsFromB->Add(hCosPcPDLSCfromB);

 // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0LSCfromB_pt";
    namehist+=i;
    titlehist="d0xd0 Loose Cuts FromB ptbin=";
    titlehist+=i;
    hd0xd0LSCfromBpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0LSCfromBpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0LSCfromBpt->SetYTitle("Entries");
    flistLsCutsFromB->Add(hd0xd0LSCfromBpt);


    namehist="hd0D0VSd0xd0LSCfromB_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} Loose Cuts FromB ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0LSCfromBpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0LSCfromBpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0LSCfromBpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistLsCutsFromB->Add(hd0D0VSd0xd0LSCfromBpt);
    
    
    namehist="hangletracksVSd0xd0LSCfromB_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} Loose Cuts FromB ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0LSCfromBpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0LSCfromBpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0LSCfromBpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistLsCutsFromB->Add(hangletracksVSd0xd0LSCfromBpt);
    

    namehist="hangletracksVSd0D0LSCfromB_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} Loose Cuts FromB ptbin=";
    titlehist+=i;
    hangletracksVSd0D0LSCfromBpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0LSCfromBpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0LSCfromBpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistLsCutsFromB->Add(hangletracksVSd0D0LSCfromBpt);
    
  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





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
  
  TH1F *hd0D0ptLSCfromBPM;
  TH1F *hMCd0D0ptLSCfromBPM;
  TH1F *hd0D0VtxTrueptLSCfromBPM;
  TH1F *hd0D0ptLSCfromBSB;
  TH1F *hMCd0D0ptLSCfromBSB;
  TH1F *hd0D0VtxTrueptLSCfromBSB;
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
    
    hd0D0ptLSCfromBPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCfromBPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCfromBPM->SetYTitle("Entries");
    flistLsCutsFromB->Add(hd0D0ptLSCfromBPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCfromBPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCfromBPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCfromBPM->SetYTitle("Entries");
    flistLsCutsFromB->Add(hMCd0D0ptLSCfromBPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCfromBPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCfromBPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCfromBPM->SetYTitle("Entries");
    flistLsCutsFromB->Add(hd0D0VtxTrueptLSCfromBPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCfromBSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCfromBSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCfromBSB->SetYTitle("Entries");
    flistLsCutsFromB->Add(hd0D0ptLSCfromBSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCfromBSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCfromBSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCfromBSB->SetYTitle("Entries");
    flistLsCutsFromB->Add(hMCd0D0ptLSCfromBSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCfromBSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCfromBSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCfromBSB->SetYTitle("Entries");
    flistLsCutsFromB->Add(hd0D0VtxTrueptLSCfromBSB);
  }



 //############ LOOSE CUTS FROM DSTAR HISTOGRAMS ###########
 //
  //############## global properties histos
  TH2F *hCPtaVSd0d0LSCfromDstar=new TH2F("hCPtaVSd0d0LSCfromDstar","hCPtaVSd0d0_LooseCuts_FromDStar",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZLSCfromDstar=new TH1F("hSecVtxZLSCfromDstar","hSecVtxZ_LooseCuts_FromDStar",1000,-8.,8.);
  TH1F *hSecVtxXLSCfromDstar=new TH1F("hSecVtxXLSCfromDstar","hSecVtxX_LooseCuts_FromDStar",1000,-3000.,3000.);
  TH1F *hSecVtxYLSCfromDstar=new TH1F("hSecVtxYLSCfromDstar","hSecVtxY_LooseCuts_FromDStar",1000,-3000.,3000.);
  TH2F *hSecVtxXYLSCfromDstar=new TH2F("hSecVtxXYLSCfromDstar","hSecVtxXY_LooseCuts_FromDStar",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiLSCfromDstar=new TH1F("hSecVtxPhiLSCfromDstar","hSecVtxPhi_LooseCuts_FromDStar",180,-180.1,180.1);
  TH1F *hd0singlTrackLSCfromDstar=new TH1F("hd0singlTrackLSCfromDstar","hd0singlTrackLooseCuts_FromDstar",1000,-5000.,5000.);
  TH1F *hCPtaLSCfromDstar=new TH1F("hCPtaLSCfromDstar","hCPta_LooseCuts_FromDStar",100,-1.,1.);
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
  flistLsCutsFromDstar->Add(hd0singlTrackLSCfromDstar);
  flistLsCutsFromDstar->Add(hCPtaLSCfromDstar);
  flistLsCutsFromDstar->Add(hd0xd0LSCfromDstar);
  flistLsCutsFromDstar->Add(hMassTrueLSCfromDstar);
  flistLsCutsFromDstar->Add(hMassLSCfromDstar);
  flistLsCutsFromDstar->Add(hMassTrueLSCfromDstarPM);
  flistLsCutsFromDstar->Add(hMassLSCfromDstarPM);
  flistLsCutsFromDstar->Add(hMassTrueLSCfromDstarSB);
  flistLsCutsFromDstar->Add(hMassLSCfromDstarSB);







 //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaLSCfromDstar=new TH1F("hdcaLSCfromDstar","hdca_LooseCuts_FromDstar",100,0.,1000.);
  hdcaLSCfromDstar->SetXTitle("dca   [#mum]");
  hdcaLSCfromDstar->SetYTitle("Entries");
  TH1F *hcosthetastarLSCfromDstar=new TH1F("hcosthetastarLSCfromDstar","hCosThetaStar_LooseCuts_FromDstar",50,-1.,1.);
  hcosthetastarLSCfromDstar->SetXTitle("cos #theta^{*}");
  hcosthetastarLSCfromDstar->SetYTitle("Entries");
  TH1F *hptD0LSCfromDstar=new TH1F("hptD0LSCfromDstar","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0LSCfromDstar->SetXTitle("p_{t}  [GeV/c]");
  hptD0LSCfromDstar->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtLSCfromDstar=new TH1F("hptD0VsMaxPtLSCfromDstar","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtLSCfromDstar=new TH2F("hptD0PTallsqrtLSCfromDstar","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallLSCfromDstar=new TH2F("hptD0PTallLSCfromDstar","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBLSCfromDstar=new TH2F("hptD0vsptBLSCfromDstar","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBLSCfromDstar=new TH2F("hpD0vspBLSCfromDstar","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkLSCfromDstar=new TH2F("hptD0vsptcquarkLSCfromDstar","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkLSCfromDstar=new TH2F("hpD0vspcquarkLSCfromDstar","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistLsCutsFromDstar->Add(hdcaLSCfromDstar);
  flistLsCutsFromDstar->Add(hcosthetastarLSCfromDstar);
  flistLsCutsFromDstar->Add(hptD0LSCfromDstar);
  flistLsCutsFromDstar->Add(hptD0VsMaxPtLSCfromDstar);
  flistLsCutsFromDstar->Add(hptD0PTallsqrtLSCfromDstar);
  flistLsCutsFromDstar->Add(hptD0PTallLSCfromDstar);
  flistLsCutsFromDstar->Add(hptD0vsptBLSCfromDstar);
  flistLsCutsFromDstar->Add(hpD0vspBLSCfromDstar);
  flistLsCutsFromDstar->Add(hptD0vsptcquarkLSCfromDstar);
  flistLsCutsFromDstar->Add(hpD0vspcquarkLSCfromDstar);
 
  TH1F *hd0zD0ptLSCfromDstar;
  TH1F *hInvMassD0LSCfromDstar,*hInvMassD0barLSCfromDstar;
  TH2F *hInvMassPtLSCfromDstar=new TH2F("hInvMassPtLSCfromDstar","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
 THnSparseF *hSparseLSCfromDstar=new THnSparseF("hSparseLSCfromDstar","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseLSCfromDstar->SetBinEdges(0,massbins);
  hSparseLSCfromDstar->SetBinEdges(1,massbins);
  hSparseLSCfromDstar->SetBinEdges(2,ptbinsForNsparse);
  hSparseLSCfromDstar->SetBinEdges(3,impparbins);
  hSparseLSCfromDstar->SetBinEdges(4,massHypoBins); 
  flistLsCutsFromDstar->Add(hSparseLSCfromDstar);
  TH1F *hetaLSCfromDstar;
  TH1F *hCosPDPBLSCfromDstar;
  TH1F *hCosPcPDLSCfromDstar;
  flistLsCutsFromDstar->Add(hInvMassPtLSCfromDstar);
  // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0LSCfromDstarpt;
  TH2F *hangletracksVSd0xd0LSCfromDstarpt;
  TH2F *hangletracksVSd0D0LSCfromDstarpt;
  TH1F *hd0xd0LSCfromDstarpt;

  TH2F *hTOFpidLSCfromDstar=new TH2F("hTOFpidLSCfromDstar","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistLsCutsFromDstar->Add(hTOFpidLSCfromDstar);

  for(Int_t i=0;i<fnbins;i++){
    namehist="hd0zD0ptLSCfromDstar_pt";
    namehist+=i;
    titlehist="d0(z) Loose Cuts FromDstarm ptbin=";
    titlehist+=i;
    hd0zD0ptLSCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptLSCfromDstar->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptLSCfromDstar->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hd0zD0ptLSCfromDstar);

    namehist="hInvMassD0LSCfromDstar_pt";
    namehist+=i;
    titlehist="Invariant Mass Loose Cuts FromDstar ptbin=";
    titlehist+=i;
    hInvMassD0LSCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0LSCfromDstar->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0LSCfromDstar->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hInvMassD0LSCfromDstar);

    namehist="hInvMassD0barLSCfromDstar_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar Loose Cuts FromDstar ptbin=";
    titlehist+=i;
    hInvMassD0barLSCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barLSCfromDstar->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barLSCfromDstar->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hInvMassD0barLSCfromDstar);

    namehist="hetaLSCfromDstar_pt";
    namehist+=i;
    titlehist="eta Loose Cuts FromDstar ptbin=";
    titlehist+=i;
    hetaLSCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaLSCfromDstar->SetXTitle("Pseudorapidity");
    hetaLSCfromDstar->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hetaLSCfromDstar);

    namehist="hCosPDPBLSCfromDstar_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBLSCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBLSCfromDstar->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBLSCfromDstar->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hCosPDPBLSCfromDstar);

    namehist="hCosPcPDLSCfromDstar_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDLSCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDLSCfromDstar->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDLSCfromDstar->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hCosPcPDLSCfromDstar);
    
 // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0LSCfromDstar_pt";
    namehist+=i;
    titlehist="d0xd0 Loose Cuts FromDstar ptbin=";
    titlehist+=i;
    hd0xd0LSCfromDstarpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0LSCfromDstarpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0LSCfromDstarpt->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hd0xd0LSCfromDstarpt);


    namehist="hd0D0VSd0xd0LSCfromDstar_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} Loose Cuts FromDstar ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0LSCfromDstarpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0LSCfromDstarpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0LSCfromDstarpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistLsCutsFromDstar->Add(hd0D0VSd0xd0LSCfromDstarpt);
    
    
    namehist="hangletracksVSd0xd0LSCfromDstar_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} Loose Cuts FromDstar ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0LSCfromDstarpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0LSCfromDstarpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0LSCfromDstarpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistLsCutsFromDstar->Add(hangletracksVSd0xd0LSCfromDstarpt);
    

    namehist="hangletracksVSd0D0LSCfromDstar_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} Loose Cuts FromDstar ptbin=";
    titlehist+=i;
    hangletracksVSd0D0LSCfromDstarpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0LSCfromDstarpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0LSCfromDstarpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistLsCutsFromDstar->Add(hangletracksVSd0D0LSCfromDstarpt);


  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






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
  
  TH1F *hd0D0ptLSCfromDstPM;
  TH1F *hMCd0D0ptLSCfromDstPM;
  TH1F *hd0D0VtxTrueptLSCfromDstPM;
  TH1F *hd0D0ptLSCfromDstSB;
  TH1F *hMCd0D0ptLSCfromDstSB;
  TH1F *hd0D0VtxTrueptLSCfromDstSB;
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
    
    hd0D0ptLSCfromDstPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCfromDstPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCfromDstPM->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hd0D0ptLSCfromDstPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCfromDstPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCfromDstPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCfromDstPM->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hMCd0D0ptLSCfromDstPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCfromDstPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCfromDstPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCfromDstPM->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hd0D0VtxTrueptLSCfromDstPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCfromDstSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCfromDstSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCfromDstSB->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hd0D0ptLSCfromDstSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCfromDstSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCfromDstSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCfromDstSB->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hMCd0D0ptLSCfromDstSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCfromDstSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCfromDstSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCfromDstSB->SetYTitle("Entries");
    flistLsCutsFromDstar->Add(hd0D0VtxTrueptLSCfromDstSB);
  }


  //############ LOOSE CUTS OTHER HISTOGRAMS ###########
  //
  //########### global properties histos ###########

  TH2F *hCPtaVSd0d0LSCother=new TH2F("hCPtaVSd0d0LSCother","hCPtaVSd0d0_LooseCuts_other",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZLSCother=new TH1F("hSecVtxZLSCother","hSecVtxZ_LooseCuts_other",1000,-8.,8.);
  TH1F *hSecVtxXLSCother=new TH1F("hSecVtxXLSCother","hSecVtxX_LooseCuts_other",1000,-3000.,3000.);
  TH1F *hSecVtxYLSCother=new TH1F("hSecVtxYLSCother","hSecVtxY_LooseCuts_other",1000,-3000.,3000.);
  TH2F *hSecVtxXYLSCother=new TH2F("hSecVtxXYLSCother","hSecVtxXY_LooseCuts_other",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiLSCother=new TH1F("hSecVtxPhiLSCother","hSecVtxPhi_LooseCuts_other",180,-180.1,180.1);
  TH1F *hd0singlTrackLSCother=new TH1F("hd0singlTrackLSCother","hd0singlTrackLooseCuts_Other",1000,-5000.,5000.);
  TH1F *hCPtaLSCother=new TH1F("hCPtaLSCother","hCPta_LooseCuts_other",100,-1.,1.);
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
  flistLsCutsOther->Add(hd0singlTrackLSCother);
  flistLsCutsOther->Add(hCPtaLSCother);
  flistLsCutsOther->Add(hd0xd0LSCother);
  flistLsCutsOther->Add(hMassTrueLSCother);
  flistLsCutsOther->Add(hMassLSCother);
  flistLsCutsOther->Add(hMassTrueLSCotherPM);
  flistLsCutsOther->Add(hMassLSCotherPM);
  flistLsCutsOther->Add(hMassTrueLSCotherSB);
  flistLsCutsOther->Add(hMassLSCotherSB);




 //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaLSCother=new TH1F("hdcaLSCother","hdca_LooseCuts_Other",100,0.,1000.);
  hdcaLSCother->SetXTitle("dca   [#mum]");
  hdcaLSCother->SetYTitle("Entries");
  TH1F *hcosthetastarLSCother=new TH1F("hcosthetastarLSCother","hCosThetaStar_LooseCuts_Other",50,-1.,1.);
  hcosthetastarLSCother->SetXTitle("cos #theta^{*}");
  hcosthetastarLSCother->SetYTitle("Entries");
  TH1F *hptD0LSCother=new TH1F("hptD0LSCother","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0LSCother->SetXTitle("p_{t}  [GeV/c]");
  hptD0LSCother->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtLSCother=new TH1F("hptD0VsMaxPtLSCother","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtLSCother=new TH2F("hptD0PTallsqrtLSCother","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallLSCother=new TH2F("hptD0PTallLSCother","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBLSCother=new TH2F("hptD0vsptBLSCother","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBLSCother=new TH2F("hpD0vspBLSCother","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkLSCother=new TH2F("hptD0vsptcquarkLSCother","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkLSCother=new TH2F("hpD0vspcquarkLSCother","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistLsCutsOther->Add(hdcaLSCother);
  flistLsCutsOther->Add(hcosthetastarLSCother);
  flistLsCutsOther->Add(hptD0LSCother);
  flistLsCutsOther->Add(hptD0VsMaxPtLSCother);
  flistLsCutsOther->Add(hptD0PTallsqrtLSCother);
  flistLsCutsOther->Add(hptD0PTallLSCother);
  flistLsCutsOther->Add(hptD0vsptBLSCother);
  flistLsCutsOther->Add(hpD0vspBLSCother);
  flistLsCutsOther->Add(hptD0vsptcquarkLSCother);
  flistLsCutsOther->Add(hpD0vspcquarkLSCother);

  TH1F *hd0zD0ptLSCother;
  TH1F *hInvMassD0LSCother,*hInvMassD0barLSCother;
  TH2F *hInvMassPtLSCother=new TH2F("hInvMassPtLSCother","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
 THnSparseF *hSparseLSCother=new THnSparseF("hSparseLSCother","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseLSCother->SetBinEdges(0,massbins);
  hSparseLSCother->SetBinEdges(1,massbins);
  hSparseLSCother->SetBinEdges(2,ptbinsForNsparse);
  hSparseLSCother->SetBinEdges(3,impparbins);
  hSparseLSCother->SetBinEdges(4,massHypoBins); 
  flistLsCutsOther->Add(hSparseLSCother);
  TH1F *hetaLSCother;
  TH1F *hCosPDPBLSCother;
  TH1F *hCosPcPDLSCother;
  flistLsCutsOther->Add(hInvMassPtLSCother);
 // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0LSCotherpt;
  TH2F *hangletracksVSd0xd0LSCotherpt;
  TH2F *hangletracksVSd0D0LSCotherpt;
  TH1F *hd0xd0LSCotherpt;

  TH2F *hTOFpidLSCother=new TH2F("hTOFpidLSCother","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistLsCutsOther->Add(hTOFpidLSCother);

  for(Int_t i=0;i<fnbins;i++){
    namehist="hd0zD0ptLSCother_pt";
    namehist+=i;
    titlehist="d0(z) Loose Cuts Otherm ptbin=";
    titlehist+=i;
    hd0zD0ptLSCother=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptLSCother->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptLSCother->SetYTitle("Entries");
    flistLsCutsOther->Add(hd0zD0ptLSCother);

    namehist="hInvMassD0LSCother_pt";
    namehist+=i;
    titlehist="Invariant Mass Loose Cuts Other ptbin=";
    titlehist+=i;
    hInvMassD0LSCother=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0LSCother->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0LSCother->SetYTitle("Entries");
    flistLsCutsOther->Add(hInvMassD0LSCother);

    namehist="hInvMassD0barLSCother_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar Loose Cuts Other ptbin=";
    titlehist+=i;
    hInvMassD0barLSCother=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barLSCother->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barLSCother->SetYTitle("Entries");
    flistLsCutsOther->Add(hInvMassD0barLSCother);

    namehist="hetaLSCother_pt";
    namehist+=i;
    titlehist="eta Loose Cuts Other ptbin=";
    titlehist+=i;
    hetaLSCother=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaLSCother->SetXTitle("Pseudorapidity");
    hetaLSCother->SetYTitle("Entries");
    flistLsCutsOther->Add(hetaLSCother);

    namehist="hCosPDPBLSCother_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBLSCother=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBLSCother->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBLSCother->SetYTitle("Entries");
    flistLsCutsOther->Add(hCosPDPBLSCother);

    namehist="hCosPcPDLSCother_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDLSCother=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDLSCother->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDLSCother->SetYTitle("Entries");
    flistLsCutsOther->Add(hCosPcPDLSCother);

 // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0LSCother_pt";
    namehist+=i;
    titlehist="d0xd0 Loose Cuts Other ptbin=";
    titlehist+=i;
    hd0xd0LSCotherpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0LSCotherpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0LSCotherpt->SetYTitle("Entries");
    flistLsCutsOther->Add(hd0xd0LSCotherpt);


    namehist="hd0D0VSd0xd0LSCother_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} Loose Cuts Other ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0LSCotherpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0LSCotherpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0LSCotherpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistLsCutsOther->Add(hd0D0VSd0xd0LSCotherpt);
    
    
    namehist="hangletracksVSd0xd0LSCother_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} Loose Cuts Other ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0LSCotherpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0LSCotherpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0LSCotherpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistLsCutsOther->Add(hangletracksVSd0xd0LSCotherpt);
    

    namehist="hangletracksVSd0D0LSCother_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} Loose Cuts Other ptbin=";
    titlehist+=i;
    hangletracksVSd0D0LSCotherpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0LSCotherpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0LSCotherpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistLsCutsOther->Add(hangletracksVSd0D0LSCotherpt);

    
  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
  
  TH1F *hd0D0ptLSCotherPM;
  TH1F *hMCd0D0ptLSCotherPM;
  TH1F *hd0D0VtxTrueptLSCotherPM;
  TH1F *hd0D0ptLSCotherSB;
  TH1F *hMCd0D0ptLSCotherSB;
  TH1F *hd0D0VtxTrueptLSCotherSB;
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
    
    hd0D0ptLSCotherPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCotherPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCotherPM->SetYTitle("Entries");
    flistLsCutsOther->Add(hd0D0ptLSCotherPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCotherPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCotherPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCotherPM->SetYTitle("Entries");
    flistLsCutsOther->Add(hMCd0D0ptLSCotherPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCotherPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCotherPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCotherPM->SetYTitle("Entries");
    flistLsCutsOther->Add(hd0D0VtxTrueptLSCotherPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptLSCotherSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptLSCotherSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptLSCotherSB->SetYTitle("Entries");
    flistLsCutsOther->Add(hd0D0ptLSCotherSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptLSCotherSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptLSCotherSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptLSCotherSB->SetYTitle("Entries");
    flistLsCutsOther->Add(hMCd0D0ptLSCotherSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptLSCotherSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptLSCotherSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptLSCotherSB->SetYTitle("Entries");
    flistLsCutsOther->Add(hd0D0VtxTrueptLSCotherSB);
  }

  //Printf("END OF LSCUTS HISTOS CREATION \n");


  //################################################################################################
  //                                                                                               #
  //                         HISTOS FOR TIGHT CUTS                                                 #
  //                                                                                               #
  //################################################################################################

  //############ TIGHT CUTS SIGNAL HISTOGRAMS ###############
  //
  // ####### global properties histo ############

  TH2F *hCPtaVSd0d0TGHCsign=new TH2F("hCPtaVSd0d0TGHCsign","hCPtaVSd0d0_TightCuts_Signal",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZTGHCsign=new TH1F("hSecVtxZTGHCsign","hSecVtxZ_TightCuts_Signal",1000,-8.,8.);
  TH1F *hSecVtxXTGHCsign=new TH1F("hSecVtxXTGHCsign","hSecVtxX_TightCuts_Signal",1000,-3000.,3000.);
  TH1F *hSecVtxYTGHCsign=new TH1F("hSecVtxYTGHCsign","hSecVtxY_TightCuts_Signal",1000,-3000.,3000.);
  TH2F *hSecVtxXYTGHCsign=new TH2F("hSecVtxXYTGHCsign","hSecVtxXY_TightCuts_Signal",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiTGHCsign=new TH1F("hSecVtxPhiTGHCsign","hSecVtxPhi_TightCuts_Signal",180,-180.1,180.1);
  TH1F *hd0singlTrackTGHCsign=new TH1F("hd0singlTrackTGHCsign","hd0singlTrackTightCuts_Signal",1000,-5000.,5000.);
  TH1F *hCPtaTGHCsign=new TH1F("hCPtaTGHCsign","hCPta_TightCuts_Signal",100,-1.,1.);
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
  flistTghCutsSignal->Add(hd0singlTrackTGHCsign);
  flistTghCutsSignal->Add(hCPtaTGHCsign);
  flistTghCutsSignal->Add(hd0xd0TGHCsign);
  flistTghCutsSignal->Add(hMassTrueTGHCsign);
  flistTghCutsSignal->Add(hMassTGHCsign);
  flistTghCutsSignal->Add(hMassTrueTGHCsignPM);
  flistTghCutsSignal->Add(hMassTGHCsignPM);
  flistTghCutsSignal->Add(hMassTrueTGHCsignSB);
  flistTghCutsSignal->Add(hMassTGHCsignSB);






 //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaTGHCsign=new TH1F("hdcaTGHCsign","hdca_TightCuts_Signal",100,0.,1000.);
  hdcaTGHCsign->SetXTitle("dca   [#mum]");
  hdcaTGHCsign->SetYTitle("Entries");
  TH1F *hcosthetastarTGHCsign=new TH1F("hcosthetastarTGHCsign","hCosThetaStar_TightCuts_Signal",50,-1.,1.);
  hcosthetastarTGHCsign->SetXTitle("cos #theta^{*}");
  hcosthetastarTGHCsign->SetYTitle("Entries");
  TH1F *hptD0TGHCsign=new TH1F("hptD0TGHCsign","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0TGHCsign->SetXTitle("p_{t}  [GeV/c]");
  hptD0TGHCsign->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtTGHCsign=new TH1F("hptD0VsMaxPtTGHCsign","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtTGHCsign=new TH2F("hptD0PTallsqrtTGHCsign","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallTGHCsign=new TH2F("hptD0PTallTGHCsign","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBTGHCsign=new TH2F("hptD0vsptBTGHCsign","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBTGHCsign=new TH2F("hpD0vspBTGHCsign","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkTGHCsign=new TH2F("hptD0vsptcquarkTGHCsign","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkTGHCsign=new TH2F("hpD0vspcquarkTGHCsign","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistTghCutsSignal->Add(hdcaTGHCsign);
  flistTghCutsSignal->Add(hcosthetastarTGHCsign);
  flistTghCutsSignal->Add(hptD0TGHCsign);
  flistTghCutsSignal->Add(hptD0VsMaxPtTGHCsign);
  flistTghCutsSignal->Add(hptD0PTallsqrtTGHCsign);
  flistTghCutsSignal->Add(hptD0PTallTGHCsign);
  flistTghCutsSignal->Add(hptD0vsptBTGHCsign);
  flistTghCutsSignal->Add(hpD0vspBTGHCsign);
  flistTghCutsSignal->Add(hptD0vsptcquarkTGHCsign);
  flistTghCutsSignal->Add(hpD0vspcquarkTGHCsign);
 
  TH1F *hd0zD0ptTGHCsign;
  TH1F *hInvMassD0TGHCsign,*hInvMassD0barTGHCsign;
  TH2F *hInvMassPtTGHCsign=new TH2F("hInvMassPtTGHCsign","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
 THnSparseF *hSparseTGHCsign=new THnSparseF("hSparseTGHCsign","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseTGHCsign->SetBinEdges(0,massbins);
  hSparseTGHCsign->SetBinEdges(1,massbins);
  hSparseTGHCsign->SetBinEdges(2,ptbinsForNsparse);
  hSparseTGHCsign->SetBinEdges(3,impparbins);
  hSparseTGHCsign->SetBinEdges(4,massHypoBins); 
  flistTghCutsSignal->Add(hSparseTGHCsign);

  THnSparseF *hSparseRecoTGHCfromB=new THnSparseF("hSparseRecoTGHCfromB","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseRecoTGHCfromB->SetBinEdges(0,massbins);
  hSparseRecoTGHCfromB->SetBinEdges(1,massbins);
  hSparseRecoTGHCfromB->SetBinEdges(2,ptbinsForNsparse);
  hSparseRecoTGHCfromB->SetBinEdges(3,impparbins);
  hSparseRecoTGHCfromB->SetBinEdges(4,massHypoBins); 
  flistTghCutsFromB->Add(hSparseRecoTGHCfromB);



  THnSparseF *hSparseCxyLxyTGHCsign=new THnSparseF("hSparseCxyLxyTGHCsign","Candidate Mass;massD0;Pt;CosXY;Lxy",4,nbinsSparsCxyLxy,binLowLimitSparseCxyLxy,binUpLimitSparseCxyLxy);
  hSparseCxyLxyTGHCsign->SetBinEdges(1,ptbinlimitsCxyLxy);
  hSparseCxyLxyTGHCsign->GetAxis(0)->SetName("mass");
  hSparseCxyLxyTGHCsign->GetAxis(0)->SetTitle("Invariant Mass (K#pi) [GeV/c^{2}]");
  hSparseCxyLxyTGHCsign->GetAxis(1)->SetName("pt");
  hSparseCxyLxyTGHCsign->GetAxis(1)->SetTitle("p_{t} [GeV/c]");
  hSparseCxyLxyTGHCsign->GetAxis(2)->SetName("CosPointXY");
  hSparseCxyLxyTGHCsign->GetAxis(2)->SetTitle("Cos#theta_{point}^{XY}");
  hSparseCxyLxyTGHCsign->GetAxis(3)->SetName("NormDecLengthXY");
  hSparseCxyLxyTGHCsign->GetAxis(3)->SetTitle("Normalized XY decay length");


  flistTghCutsSignal->Add(hSparseCxyLxyTGHCsign);
  
  
  TH1F *hetaTGHCsign;
  TH1F *hCosPDPBTGHCsign;
  TH1F *hCosPcPDTGHCsign;
  flistTghCutsSignal->Add(hInvMassPtTGHCsign);
// %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0TGHCsignpt;
  TH2F *hangletracksVSd0xd0TGHCsignpt;
  TH2F *hangletracksVSd0D0TGHCsignpt;
  TH1F *hd0xd0TGHCsignpt;
  TH1F *hPhiHistPMTGHCsignpt,*hPhiHistSBTGHCsignpt;

  TH2F *hTOFpidTGHCsign=new TH2F("hTOFpidTGHCsign","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistTghCutsSignal->Add(hTOFpidTGHCsign);

  for(Int_t i=0;i<fnbins;i++){
    // Printf("INSIDE FIRST LOOP FOR  TIGHT CUTS HISTO CREATION %d\n", fnbins);
    namehist="hPhiHistPMTGHCsign_pt";
    namehist+=i;
    titlehist="Azimuthal correlation TGH Cuts Sign PM ptbin=";
    titlehist+=i;
    hPhiHistPMTGHCsignpt=new TH1F(namehist.Data(),titlehist.Data(),100,-3.15,3.15);
    hPhiHistPMTGHCsignpt->Sumw2();
    flistTghCutsSignal->Add(hPhiHistPMTGHCsignpt);

    namehist="hPhiHistSBTGHCsign_pt";
    namehist+=i;
    titlehist="Azimuthal correlation TGH Cuts Sign SB ptbin=";
    titlehist+=i;
    hPhiHistSBTGHCsignpt=new TH1F(namehist.Data(),titlehist.Data(),100,-3.15,3.15);
    hPhiHistSBTGHCsignpt->Sumw2();
    flistTghCutsSignal->Add(hPhiHistSBTGHCsignpt);

    namehist="hd0zD0ptTGHCsign_pt";
    namehist+=i;
    titlehist="d0(z) Tight Cuts Signal ptbin=";
    titlehist+=i;
    hd0zD0ptTGHCsign=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptTGHCsign->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptTGHCsign->SetYTitle("Entries");
    flistTghCutsSignal->Add(hd0zD0ptTGHCsign);

    namehist="hInvMassD0TGHCsign_pt";
    namehist+=i;
    titlehist="Invariant Mass Tight Cuts Signal ptbin=";
    titlehist+=i;
    hInvMassD0TGHCsign=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0TGHCsign->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0TGHCsign->SetYTitle("Entries");
    flistTghCutsSignal->Add(hInvMassD0TGHCsign);

    namehist="hInvMassD0barTGHCsign_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar Tight Cuts Signal ptbin=";
    titlehist+=i;
    hInvMassD0barTGHCsign=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barTGHCsign->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barTGHCsign->SetYTitle("Entries");
    flistTghCutsSignal->Add(hInvMassD0barTGHCsign);


    namehist="hetaTGHCsign_pt";
    namehist+=i;
    titlehist="eta Tight Cuts Signal ptbin=";
    titlehist+=i;
    hetaTGHCsign=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaTGHCsign->SetXTitle("Pseudorapidity");
    hetaTGHCsign->SetYTitle("Entries");
    flistTghCutsSignal->Add(hetaTGHCsign);

    namehist="hCosPDPBTGHCsign_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBTGHCsign=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBTGHCsign->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBTGHCsign->SetYTitle("Entries");
    flistTghCutsSignal->Add(hCosPDPBTGHCsign);

    namehist="hCosPcPDTGHCsign_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDTGHCsign=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDTGHCsign->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDTGHCsign->SetYTitle("Entries");
    flistTghCutsSignal->Add(hCosPcPDTGHCsign);
   
 // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0TGHCsign_pt";
    namehist+=i;
    titlehist="d0xd0 Tight Cuts Signal ptbin=";
    titlehist+=i;
    hd0xd0TGHCsignpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0TGHCsignpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0TGHCsignpt->SetYTitle("Entries");
    flistTghCutsSignal->Add(hd0xd0TGHCsignpt);


    namehist="hd0D0VSd0xd0TGHCsign_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} Tight Cuts Signal ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0TGHCsignpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0TGHCsignpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0TGHCsignpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistTghCutsSignal->Add(hd0D0VSd0xd0TGHCsignpt);
    
    
    namehist="hangletracksVSd0xd0TGHCsign_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} Tight Cuts Signal ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0TGHCsignpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0TGHCsignpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0TGHCsignpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistTghCutsSignal->Add(hangletracksVSd0xd0TGHCsignpt);
    

    namehist="hangletracksVSd0D0TGHCsign_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} Tight Cuts Signal ptbin=";
    titlehist+=i;
    hangletracksVSd0D0TGHCsignpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0TGHCsignpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0TGHCsignpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistTghCutsSignal->Add(hangletracksVSd0D0TGHCsignpt);

  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








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
  
  TH1F *hd0D0ptTGHCsignPM;
  TH1F *hMCd0D0ptTGHCsignPM;
  TH1F *hd0D0VtxTrueptTGHCsignPM;
  TH1F *hd0D0ptTGHCsignSB;
  TH1F *hMCd0D0ptTGHCsignSB;
  TH1F *hd0D0VtxTrueptTGHCsignSB;
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
    
    hd0D0ptTGHCsignPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCsignPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCsignPM->SetYTitle("Entries");
    flistTghCutsSignal->Add(hd0D0ptTGHCsignPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCsignPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCsignPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCsignPM->SetYTitle("Entries");
    flistTghCutsSignal->Add(hMCd0D0ptTGHCsignPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCsignPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCsignPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCsignPM->SetYTitle("Entries");
    flistTghCutsSignal->Add(hd0D0VtxTrueptTGHCsignPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCsignSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCsignSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCsignSB->SetYTitle("Entries");
    flistTghCutsSignal->Add(hd0D0ptTGHCsignSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCsignSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCsignSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCsignSB->SetYTitle("Entries");
    flistTghCutsSignal->Add(hMCd0D0ptTGHCsignSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCsignSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCsignSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCsignSB->SetYTitle("Entries");
    flistTghCutsSignal->Add(hd0D0VtxTrueptTGHCsignSB);
  }


  //############ TIGHT CUTS BACKGROUND HISTOGRAMS ###########
  //
  //   ######## global properties histos #######
  TH2F *hCPtaVSd0d0TGHCback=new TH2F("hCPtaVSd0d0TGHCback","hCPtaVSd0d0_TightCuts_Background",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZTGHCback=new TH1F("hSecVtxZTGHCback","hSecVtxZ_TightCuts_Background",1000,-8.,8.);
  TH1F *hSecVtxXTGHCback=new TH1F("hSecVtxXTGHCback","hSecVtxX_TightCuts_Background",1000,-3000.,3000.);
  TH1F *hSecVtxYTGHCback=new TH1F("hSecVtxYTGHCback","hSecVtxY_TightCuts_Background",1000,-3000.,3000.);
  TH2F *hSecVtxXYTGHCback=new TH2F("hSecVtxXYTGHCback","hSecVtxXY_TightCuts_Background",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiTGHCback=new TH1F("hSecVtxPhiTGHCback","hSecVtxPhi_TightCuts_Background",180,-180.1,180.1);
  TH1F *hd0singlTrackTGHCback=new TH1F("hd0singlTrackTGHCback","hd0singlTrackTightCuts_Back",1000,-5000.,5000.);
  TH1F *hCPtaTGHCback=new TH1F("hCPtaTGHCback","hCPta_TightCuts_Background",100,-1.,1.);
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
  flistTghCutsBack->Add(hd0singlTrackTGHCback);
  flistTghCutsBack->Add(hCPtaTGHCback);
  flistTghCutsBack->Add(hd0xd0TGHCback);
  flistTghCutsBack->Add(hMassTrueTGHCback);
  flistTghCutsBack->Add(hMassTGHCback);
  flistTghCutsBack->Add(hMassTrueTGHCbackPM);
  flistTghCutsBack->Add(hMassTGHCbackPM);
  flistTghCutsBack->Add(hMassTrueTGHCbackSB);
  flistTghCutsBack->Add(hMassTGHCbackSB);







 //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaTGHCback=new TH1F("hdcaTGHCback","hdca_TightCuts_Backgr",100,0.,1000.);
  hdcaTGHCback->SetXTitle("dca   [#mum]");
  hdcaTGHCback->SetYTitle("Entries");
  TH1F *hcosthetastarTGHCback=new TH1F("hcosthetastarTGHCback","hCosThetaStar_TightCuts_Backgr",50,-1.,1.);
  hcosthetastarTGHCback->SetXTitle("cos #theta^{*}");
  hcosthetastarTGHCback->SetYTitle("Entries");
  TH1F *hptD0TGHCback=new TH1F("hptD0TGHCback","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0TGHCback->SetXTitle("p_{t}  [GeV/c]");
  hptD0TGHCback->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtTGHCback=new TH1F("hptD0VsMaxPtTGHCback","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtTGHCback=new TH2F("hptD0PTallsqrtTGHCback","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallTGHCback=new TH2F("hptD0PTallTGHCback","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBTGHCback=new TH2F("hptD0vsptBTGHCback","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBTGHCback=new TH2F("hpD0vspBTGHCback","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkTGHCback=new TH2F("hptD0vsptcquarkTGHCback","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkTGHCback=new TH2F("hpD0vspcquarkTGHCback","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistTghCutsBack->Add(hdcaTGHCback);
  flistTghCutsBack->Add(hcosthetastarTGHCback);
  flistTghCutsBack->Add(hptD0TGHCback);
  flistTghCutsBack->Add(hptD0VsMaxPtTGHCback);
  flistTghCutsBack->Add(hptD0PTallsqrtTGHCback);
  flistTghCutsBack->Add(hptD0PTallTGHCback);
  flistTghCutsBack->Add(hptD0vsptBTGHCback);
  flistTghCutsBack->Add(hpD0vspBTGHCback);
  flistTghCutsBack->Add(hptD0vsptcquarkTGHCback);
  flistTghCutsBack->Add(hpD0vspcquarkTGHCback);
 
  TH1F *hd0zD0ptTGHCback;
  TH1F *hInvMassD0TGHCback,*hInvMassD0barTGHCback;
  TH2F *hInvMassPtTGHCback=new TH2F("hInvMassPtTGHCback","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
 THnSparseF *hSparseTGHCback=new THnSparseF("hSparseTGHCback","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseTGHCback->SetBinEdges(0,massbins);
  hSparseTGHCback->SetBinEdges(1,massbins);
  hSparseTGHCback->SetBinEdges(2,ptbinsForNsparse);
  hSparseTGHCback->SetBinEdges(3,impparbins);
  hSparseTGHCback->SetBinEdges(4,massHypoBins); 
  flistTghCutsBack->Add(hSparseTGHCback);
  TH1F *hetaTGHCback;
  TH1F *hCosPDPBTGHCback;
  TH1F *hCosPcPDTGHCback;
  flistTghCutsBack->Add(hInvMassPtTGHCback);
// %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0TGHCbackpt;
  TH2F *hangletracksVSd0xd0TGHCbackpt;
  TH2F *hangletracksVSd0D0TGHCbackpt;
  TH1F *hd0xd0TGHCbackpt;

  TH2F *hTOFpidTGHCback=new TH2F("hTOFpidTGHCback","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistTghCutsBack->Add(hTOFpidTGHCback);


  for(Int_t i=0;i<fnbins;i++){
    namehist="hd0zD0ptTGHCback_pt";
    namehist+=i;
    titlehist="d0(z) Tight Cuts Backgrm ptbin=";
    titlehist+=i;
    hd0zD0ptTGHCback=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptTGHCback->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptTGHCback->SetYTitle("Entries");
    flistTghCutsBack->Add(hd0zD0ptTGHCback);

    namehist="hInvMassD0TGHCback_pt";
    namehist+=i;
    titlehist="Invariant Mass Tight Cuts Backgr ptbin=";
    titlehist+=i;
    hInvMassD0TGHCback=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0TGHCback->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0TGHCback->SetYTitle("Entries");
    flistTghCutsBack->Add(hInvMassD0TGHCback);

    namehist="hInvMassD0barTGHCback_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar Tight Cuts Back ptbin=";
    titlehist+=i;
    hInvMassD0barTGHCback=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barTGHCback->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barTGHCback->SetYTitle("Entries");
    flistTghCutsBack->Add(hInvMassD0barTGHCback);

    namehist="hetaTGHCback_pt";
    namehist+=i;
    titlehist="eta Tight Cuts Backgr ptbin=";
    titlehist+=i;
    hetaTGHCback=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaTGHCback->SetXTitle("Pseudorapidity");
    hetaTGHCback->SetYTitle("Entries");
    flistTghCutsBack->Add(hetaTGHCback);

    namehist="hCosPDPBTGHCback_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBTGHCback=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBTGHCback->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBTGHCback->SetYTitle("Entries");
    flistTghCutsBack->Add(hCosPDPBTGHCback);

    namehist="hCosPcPDTGHCback_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDTGHCback=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDTGHCback->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDTGHCback->SetYTitle("Entries");
    flistTghCutsBack->Add(hCosPcPDTGHCback);

 // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0TGHCback_pt";
    namehist+=i;
    titlehist="d0xd0 Tight Cuts Back ptbin=";
    titlehist+=i;
    hd0xd0TGHCbackpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0TGHCbackpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0TGHCbackpt->SetYTitle("Entries");
    flistTghCutsBack->Add(hd0xd0TGHCbackpt);


    namehist="hd0D0VSd0xd0TGHCback_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} Tight Cuts Back ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0TGHCbackpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0TGHCbackpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0TGHCbackpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistTghCutsBack->Add(hd0D0VSd0xd0TGHCbackpt);
    
    
    namehist="hangletracksVSd0xd0TGHCback_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} Tight Cuts Back ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0TGHCbackpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0TGHCbackpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0TGHCbackpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistTghCutsBack->Add(hangletracksVSd0xd0TGHCbackpt);
    

    namehist="hangletracksVSd0D0TGHCback_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} Tight Cuts Back ptbin=";
    titlehist+=i;
    hangletracksVSd0D0TGHCbackpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0TGHCbackpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0TGHCbackpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistTghCutsBack->Add(hangletracksVSd0D0TGHCbackpt);

    
  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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
  
  TH1F *hd0D0ptTGHCbackPM;
  TH1F *hMCd0D0ptTGHCbackPM;
  TH1F *hd0D0VtxTrueptTGHCbackPM;
  TH1F *hd0D0ptTGHCbackSB;
  TH1F *hMCd0D0ptTGHCbackSB;
  TH1F *hd0D0VtxTrueptTGHCbackSB;
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
    
    hd0D0ptTGHCbackPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCbackPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCbackPM->SetYTitle("Entries");
    flistTghCutsBack->Add(hd0D0ptTGHCbackPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCbackPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCbackPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCbackPM->SetYTitle("Entries");
    flistTghCutsBack->Add(hMCd0D0ptTGHCbackPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCbackPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCbackPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCbackPM->SetYTitle("Entries");
    flistTghCutsBack->Add(hd0D0VtxTrueptTGHCbackPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCbackSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCbackSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCbackSB->SetYTitle("Entries");
    flistTghCutsBack->Add(hd0D0ptTGHCbackSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCbackSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCbackSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCbackSB->SetYTitle("Entries");
    flistTghCutsBack->Add(hMCd0D0ptTGHCbackSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCbackSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCbackSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCbackSB->SetYTitle("Entries");
    flistTghCutsBack->Add(hd0D0VtxTrueptTGHCbackSB);
  }



 //############ TIGHT CUTS FROMB HISTOGRAMS ###########
  //
  //#######  global properties histos

  TH2F *hCPtaVSd0d0TGHCfromB=new TH2F("hCPtaVSd0d0TGHCfromB","hCPtaVSd0d0_TightCuts_FromB",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZTGHCfromB=new TH1F("hSecVtxZTGHCfromB","hSecVtxZ_TightCuts_FromB",1000,-8.,8.);
  TH1F *hSecVtxXTGHCfromB=new TH1F("hSecVtxXTGHCfromB","hSecVtxX_TightCuts_FromB",1000,-3000.,3000.);
  TH1F *hSecVtxYTGHCfromB=new TH1F("hSecVtxYTGHCfromB","hSecVtxY_TightCuts_FromB",1000,-3000.,3000.);
  TH2F *hSecVtxXYTGHCfromB=new TH2F("hSecVtxXYTGHCfromB","hSecVtxXY_TightCuts_FromB",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiTGHCfromB=new TH1F("hSecVtxPhiTGHCfromB","hSecVtxPhi_TightCuts_FromB",180,-180.1,180.1);
  TH1F *hd0singlTrackTGHCfromB=new TH1F("hd0singlTrackTGHCfromB","hd0singlTrackTightCuts_FromB",1000,-5000.,5000.);
  TH1F *hCPtaTGHCfromB=new TH1F("hCPtaTGHCfromB","hCPta_TightCuts_FromB",100,-1.,1.);
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
  flistTghCutsFromB->Add(hd0singlTrackTGHCfromB);
  flistTghCutsFromB->Add(hCPtaTGHCfromB);
  flistTghCutsFromB->Add(hd0xd0TGHCfromB);
  flistTghCutsFromB->Add(hMassTrueTGHCfromB);
  flistTghCutsFromB->Add(hMassTGHCfromB);
  flistTghCutsFromB->Add(hMassTrueTGHCfromBPM);
  flistTghCutsFromB->Add(hMassTGHCfromBPM);
  flistTghCutsFromB->Add(hMassTrueTGHCfromBSB);
  flistTghCutsFromB->Add(hMassTGHCfromBSB);



 //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaTGHCfromB=new TH1F("hdcaTGHCfromB","hdca_TightCuts_FromB",100,0.,1000.);
  hdcaTGHCfromB->SetXTitle("dca   [#mum]");
  hdcaTGHCfromB->SetYTitle("Entries");
  TH1F *hcosthetastarTGHCfromB=new TH1F("hcosthetastarTGHCfromB","hCosThetaStar_TightCuts_FromB",50,-1.,1.);
  hcosthetastarTGHCfromB->SetXTitle("cos #theta^{*}");
  hcosthetastarTGHCfromB->SetYTitle("Entries");
  TH1F *hptD0TGHCfromB=new TH1F("hptD0TGHCfromB","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0TGHCfromB->SetXTitle("p_{t}  [GeV/c]");
  hptD0TGHCfromB->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtTGHCfromB=new TH1F("hptD0VsMaxPtTGHCfromB","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtTGHCfromB=new TH2F("hptD0PTallsqrtTGHCfromB","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallTGHCfromB=new TH2F("hptD0PTallTGHCfromB","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBTGHCfromB=new TH2F("hptD0vsptBTGHCfromB","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBTGHCfromB=new TH2F("hpD0vspBTGHCfromB","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkTGHCfromB=new TH2F("hptD0vsptcquarkTGHCfromB","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkTGHCfromB=new TH2F("hpD0vspcquarkTGHCfromB","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistTghCutsFromB->Add(hdcaTGHCfromB);
  flistTghCutsFromB->Add(hcosthetastarTGHCfromB);
  flistTghCutsFromB->Add(hptD0TGHCfromB);
  flistTghCutsFromB->Add(hptD0VsMaxPtTGHCfromB);
  flistTghCutsFromB->Add(hptD0PTallsqrtTGHCfromB);
  flistTghCutsFromB->Add(hptD0PTallTGHCfromB);
  flistTghCutsFromB->Add(hptD0vsptBTGHCfromB);
  flistTghCutsFromB->Add(hpD0vspBTGHCfromB);
  flistTghCutsFromB->Add(hptD0vsptcquarkTGHCfromB);
  flistTghCutsFromB->Add(hpD0vspcquarkTGHCfromB);
 
  TH1F *hd0zD0ptTGHCfromB;
  TH1F *hInvMassD0TGHCfromB,*hInvMassD0barTGHCfromB;
  TH2F *hInvMassPtTGHCfromB=new TH2F("hInvMassPtTGHCfromB","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
 THnSparseF *hSparseTGHCfromB=new THnSparseF("hSparseTGHCfromB","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseTGHCfromB->SetBinEdges(0,massbins);
  hSparseTGHCfromB->SetBinEdges(1,massbins);
  hSparseTGHCfromB->SetBinEdges(2,ptbinsForNsparse);
  hSparseTGHCfromB->SetBinEdges(3,impparbins);
  hSparseTGHCfromB->SetBinEdges(4,massHypoBins); 
  flistTghCutsFromB->Add(hSparseTGHCfromB);
  TH1F *hetaTGHCfromB;
  TH1F *hCosPDPBTGHCfromB;
  TH1F *hCosPcPDTGHCfromB;
  flistTghCutsFromB->Add(hInvMassPtTGHCfromB);
// %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0TGHCfromBpt;
  TH2F *hangletracksVSd0xd0TGHCfromBpt;
  TH2F *hangletracksVSd0D0TGHCfromBpt;
  TH1F *hd0xd0TGHCfromBpt;

  TH2F *hTOFpidTGHCfromB=new TH2F("hTOFpidTGHCfromB","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistTghCutsFromB->Add(hTOFpidTGHCfromB);


  for(Int_t i=0;i<fnbins;i++){
    namehist="hd0zD0ptTGHCfromB_pt";
    namehist+=i;
    titlehist="d0(z) Tight Cuts FromBm ptbin=";
    titlehist+=i;
    hd0zD0ptTGHCfromB=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptTGHCfromB->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptTGHCfromB->SetYTitle("Entries");
    flistTghCutsFromB->Add(hd0zD0ptTGHCfromB);

    namehist="hInvMassD0TGHCfromB_pt";
    namehist+=i;
    titlehist="Invariant Mass Tight Cuts FromB ptbin=";
    titlehist+=i;
    hInvMassD0TGHCfromB=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0TGHCfromB->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0TGHCfromB->SetYTitle("Entries");
    flistTghCutsFromB->Add(hInvMassD0TGHCfromB);

    namehist="hInvMassD0barTGHCfromB_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar Tight Cuts FromB ptbin=";
    titlehist+=i;
    hInvMassD0barTGHCfromB=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barTGHCfromB->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barTGHCfromB->SetYTitle("Entries");
    flistTghCutsFromB->Add(hInvMassD0barTGHCfromB);

    namehist="hetaTGHCfromB_pt";
    namehist+=i;
    titlehist="eta Tight Cuts FromB ptbin=";
    titlehist+=i;
    hetaTGHCfromB=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaTGHCfromB->SetXTitle("Pseudorapidity");
    hetaTGHCfromB->SetYTitle("Entries");
    flistTghCutsFromB->Add(hetaTGHCfromB);

    namehist="hCosPDPBTGHCfromB_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBTGHCfromB=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBTGHCfromB->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBTGHCfromB->SetYTitle("Entries");
    flistTghCutsFromB->Add(hCosPDPBTGHCfromB);

    namehist="hCosPcPDTGHCfromB_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDTGHCfromB=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDTGHCfromB->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDTGHCfromB->SetYTitle("Entries");
    flistTghCutsFromB->Add(hCosPcPDTGHCfromB);

// %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0TGHCfromB_pt";
    namehist+=i;
    titlehist="d0xd0 Tight Cuts FromB ptbin=";
    titlehist+=i;
    hd0xd0TGHCfromBpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0TGHCfromBpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0TGHCfromBpt->SetYTitle("Entries");
    flistTghCutsFromB->Add(hd0xd0TGHCfromBpt);


    namehist="hd0D0VSd0xd0TGHCfromB_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} Tight Cuts FromB ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0TGHCfromBpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0TGHCfromBpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0TGHCfromBpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistTghCutsFromB->Add(hd0D0VSd0xd0TGHCfromBpt);
    
    
    namehist="hangletracksVSd0xd0TGHCfromB_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} Tight Cuts FromB ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0TGHCfromBpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0TGHCfromBpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0TGHCfromBpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistTghCutsFromB->Add(hangletracksVSd0xd0TGHCfromBpt);
    

    namehist="hangletracksVSd0D0TGHCfromB_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} Tight Cuts FromB ptbin=";
    titlehist+=i;
    hangletracksVSd0D0TGHCfromBpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0TGHCfromBpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0TGHCfromBpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistTghCutsFromB->Add(hangletracksVSd0D0TGHCfromBpt);
    
  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





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
  
  TH1F *hd0D0ptTGHCfromBPM;
  TH1F *hMCd0D0ptTGHCfromBPM;
  TH1F *hd0D0VtxTrueptTGHCfromBPM;
  TH1F *hd0D0ptTGHCfromBSB;
  TH1F *hMCd0D0ptTGHCfromBSB;
  TH1F *hd0D0VtxTrueptTGHCfromBSB;
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
    
    hd0D0ptTGHCfromBPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCfromBPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCfromBPM->SetYTitle("Entries");
    flistTghCutsFromB->Add(hd0D0ptTGHCfromBPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCfromBPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCfromBPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCfromBPM->SetYTitle("Entries");
    flistTghCutsFromB->Add(hMCd0D0ptTGHCfromBPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCfromBPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCfromBPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCfromBPM->SetYTitle("Entries");
    flistTghCutsFromB->Add(hd0D0VtxTrueptTGHCfromBPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCfromBSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCfromBSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCfromBSB->SetYTitle("Entries");
    flistTghCutsFromB->Add(hd0D0ptTGHCfromBSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCfromBSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCfromBSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCfromBSB->SetYTitle("Entries");
    flistTghCutsFromB->Add(hMCd0D0ptTGHCfromBSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCfromBSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCfromBSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCfromBSB->SetYTitle("Entries");
    flistTghCutsFromB->Add(hd0D0VtxTrueptTGHCfromBSB);
  }



 //############ TIGHT CUTS FROM DSTAR HISTOGRAMS ###########
 //
  //############## global properties histos
  TH2F *hCPtaVSd0d0TGHCfromDstar=new TH2F("hCPtaVSd0d0TGHCfromDstar","hCPtaVSd0d0_TightCuts_FromDStar",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZTGHCfromDstar=new TH1F("hSecVtxZTGHCfromDstar","hSecVtxZ_TightCuts_FromDStar",1000,-8.,8.);
  TH1F *hSecVtxXTGHCfromDstar=new TH1F("hSecVtxXTGHCfromDstar","hSecVtxX_TightCuts_FromDStar",1000,-3000.,3000.);
  TH1F *hSecVtxYTGHCfromDstar=new TH1F("hSecVtxYTGHCfromDstar","hSecVtxY_TightCuts_FromDStar",1000,-3000.,3000.);
  TH2F *hSecVtxXYTGHCfromDstar=new TH2F("hSecVtxXYTGHCfromDstar","hSecVtxXY_TightCuts_FromDStar",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiTGHCfromDstar=new TH1F("hSecVtxPhiTGHCfromDstar","hSecVtxPhi_TightCuts_FromDStar",180,-180.1,180.1);
  TH1F *hd0singlTrackTGHCfromDstar=new TH1F("hd0singlTrackTGHCfromDstar","hd0singlTrackTightCuts_FromDstar",1000,-5000.,5000.);
  TH1F *hCPtaTGHCfromDstar=new TH1F("hCPtaTGHCfromDstar","hCPta_TightCuts_FromDStar",100,-1.,1.);
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
  flistTghCutsFromDstar->Add(hd0singlTrackTGHCfromDstar);
  flistTghCutsFromDstar->Add(hCPtaTGHCfromDstar);
  flistTghCutsFromDstar->Add(hd0xd0TGHCfromDstar);
  flistTghCutsFromDstar->Add(hMassTrueTGHCfromDstar);
  flistTghCutsFromDstar->Add(hMassTGHCfromDstar);
  flistTghCutsFromDstar->Add(hMassTrueTGHCfromDstarPM);
  flistTghCutsFromDstar->Add(hMassTGHCfromDstarPM);
  flistTghCutsFromDstar->Add(hMassTrueTGHCfromDstarSB);
  flistTghCutsFromDstar->Add(hMassTGHCfromDstarSB);





 //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaTGHCfromDstar=new TH1F("hdcaTGHCfromDstar","hdca_TightCuts_FromDstar",100,0.,1000.);
  hdcaTGHCfromDstar->SetXTitle("dca   [#mum]");
  hdcaTGHCfromDstar->SetYTitle("Entries");
  TH1F *hcosthetastarTGHCfromDstar=new TH1F("hcosthetastarTGHCfromDstar","hCosThetaStar_TightCuts_FromDstar",50,-1.,1.);
  hcosthetastarTGHCfromDstar->SetXTitle("cos #theta^{*}");
  hcosthetastarTGHCfromDstar->SetYTitle("Entries");
  TH1F *hptD0TGHCfromDstar=new TH1F("hptD0TGHCfromDstar","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0TGHCfromDstar->SetXTitle("p_{t}  [GeV/c]");
  hptD0TGHCfromDstar->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtTGHCfromDstar=new TH1F("hptD0VsMaxPtTGHCfromDstar","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtTGHCfromDstar=new TH2F("hptD0PTallsqrtTGHCfromDstar","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallTGHCfromDstar=new TH2F("hptD0PTallTGHCfromDstar","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBTGHCfromDstar=new TH2F("hptD0vsptBTGHCfromDstar","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBTGHCfromDstar=new TH2F("hpD0vspBTGHCfromDstar","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkTGHCfromDstar=new TH2F("hptD0vsptcquarkTGHCfromDstar","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkTGHCfromDstar=new TH2F("hpD0vspcquarkTGHCfromDstar","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistTghCutsFromDstar->Add(hdcaTGHCfromDstar);
  flistTghCutsFromDstar->Add(hcosthetastarTGHCfromDstar);
  flistTghCutsFromDstar->Add(hptD0TGHCfromDstar);
  flistTghCutsFromDstar->Add(hptD0VsMaxPtTGHCfromDstar);
  flistTghCutsFromDstar->Add(hptD0PTallsqrtTGHCfromDstar);
  flistTghCutsFromDstar->Add(hptD0PTallTGHCfromDstar);
  flistTghCutsFromDstar->Add(hptD0vsptBTGHCfromDstar);
  flistTghCutsFromDstar->Add(hpD0vspBTGHCfromDstar);
  flistTghCutsFromDstar->Add(hptD0vsptcquarkTGHCfromDstar);
  flistTghCutsFromDstar->Add(hpD0vspcquarkTGHCfromDstar);
 
  TH1F *hd0zD0ptTGHCfromDstar;
  TH1F *hInvMassD0TGHCfromDstar,*hInvMassD0barTGHCfromDstar;
  TH1F *hetaTGHCfromDstar;
  TH2F *hInvMassPtTGHCfromDstar=new TH2F("hInvMassPtTGHCfromDstar","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
 THnSparseF *hSparseTGHCfromDstar=new THnSparseF("hSparseTGHCfromDstar","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseTGHCfromDstar->SetBinEdges(0,massbins);
  hSparseTGHCfromDstar->SetBinEdges(1,massbins);
  hSparseTGHCfromDstar->SetBinEdges(2,ptbinsForNsparse);
  hSparseTGHCfromDstar->SetBinEdges(3,impparbins);
  hSparseTGHCfromDstar->SetBinEdges(4,massHypoBins); 
  flistTghCutsFromDstar->Add(hSparseTGHCfromDstar);
  TH1F *hCosPDPBTGHCfromDstar;
  TH1F *hCosPcPDTGHCfromDstar;
  flistTghCutsFromDstar->Add(hInvMassPtTGHCfromDstar);
// %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0TGHCfromDstarpt;
  TH2F *hangletracksVSd0xd0TGHCfromDstarpt;
  TH2F *hangletracksVSd0D0TGHCfromDstarpt;
  TH1F *hd0xd0TGHCfromDstarpt;

  TH2F *hTOFpidTGHCfromDstar=new TH2F("hTOFpidTGHCfromDstar","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistTghCutsFromDstar->Add(hTOFpidTGHCfromDstar);

  for(Int_t i=0;i<fnbins;i++){
    namehist="hd0zD0ptTGHCfromDstar_pt";
    namehist+=i;
    titlehist="d0(z) Tight Cuts FromDstarm ptbin=";
    titlehist+=i;
    hd0zD0ptTGHCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptTGHCfromDstar->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptTGHCfromDstar->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hd0zD0ptTGHCfromDstar);

    namehist="hInvMassD0TGHCfromDstar_pt";
    namehist+=i;
    titlehist="Invariant Mass Tight Cuts FromDstar ptbin=";
    titlehist+=i;
    hInvMassD0TGHCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0TGHCfromDstar->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0TGHCfromDstar->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hInvMassD0TGHCfromDstar);

    namehist="hInvMassD0barTGHCfromDstar_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar Tight Cuts FromDstar ptbin=";
    titlehist+=i;
    hInvMassD0barTGHCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barTGHCfromDstar->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barTGHCfromDstar->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hInvMassD0barTGHCfromDstar);

    namehist="hetaTGHCfromDstar_pt";
    namehist+=i;
    titlehist="eta Tight Cuts FromDstar ptbin=";
    titlehist+=i;
    hetaTGHCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaTGHCfromDstar->SetXTitle("Pseudorapidity");
    hetaTGHCfromDstar->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hetaTGHCfromDstar);

    namehist="hCosPDPBTGHCfromDstar_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBTGHCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBTGHCfromDstar->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBTGHCfromDstar->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hCosPDPBTGHCfromDstar);

    namehist="hCosPcPDTGHCfromDstar_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDTGHCfromDstar=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDTGHCfromDstar->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDTGHCfromDstar->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hCosPcPDTGHCfromDstar);

 // %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0TGHCfromDstar_pt";
    namehist+=i;
    titlehist="d0xd0 Tight Cuts FromDstar ptbin=";
    titlehist+=i;
    hd0xd0TGHCfromDstarpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0TGHCfromDstarpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0TGHCfromDstarpt->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hd0xd0TGHCfromDstarpt);


    namehist="hd0D0VSd0xd0TGHCfromDstar_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} Tight Cuts FromDstar ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0TGHCfromDstarpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0TGHCfromDstarpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0TGHCfromDstarpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistTghCutsFromDstar->Add(hd0D0VSd0xd0TGHCfromDstarpt);
    
    
    namehist="hangletracksVSd0xd0TGHCfromDstar_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} Tight Cuts FromDstar ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0TGHCfromDstarpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0TGHCfromDstarpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0TGHCfromDstarpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistTghCutsFromDstar->Add(hangletracksVSd0xd0TGHCfromDstarpt);
    

    namehist="hangletracksVSd0D0TGHCfromDstar_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} Tight Cuts FromDstar ptbin=";
    titlehist+=i;
    hangletracksVSd0D0TGHCfromDstarpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0TGHCfromDstarpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0TGHCfromDstarpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistTghCutsFromDstar->Add(hangletracksVSd0D0TGHCfromDstarpt);

    
  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
  
  TH1F *hd0D0ptTGHCfromDstPM;
  TH1F *hMCd0D0ptTGHCfromDstPM;
  TH1F *hd0D0VtxTrueptTGHCfromDstPM;
  TH1F *hd0D0ptTGHCfromDstSB;
  TH1F *hMCd0D0ptTGHCfromDstSB;
  TH1F *hd0D0VtxTrueptTGHCfromDstSB;
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
    
    hd0D0ptTGHCfromDstPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCfromDstPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCfromDstPM->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hd0D0ptTGHCfromDstPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCfromDstPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCfromDstPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCfromDstPM->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hMCd0D0ptTGHCfromDstPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCfromDstPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCfromDstPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCfromDstPM->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hd0D0VtxTrueptTGHCfromDstPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCfromDstSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCfromDstSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCfromDstSB->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hd0D0ptTGHCfromDstSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCfromDstSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCfromDstSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCfromDstSB->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hMCd0D0ptTGHCfromDstSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCfromDstSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCfromDstSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCfromDstSB->SetYTitle("Entries");
    flistTghCutsFromDstar->Add(hd0D0VtxTrueptTGHCfromDstSB);
  }


  //############ TIGHT CUTS OTHER HISTOGRAMS ###########
  //
  //########### global properties histos ###########

  TH2F *hCPtaVSd0d0TGHCother=new TH2F("hCPtaVSd0d0TGHCother","hCPtaVSd0d0_TightCuts_other",1000,-100000.,100000.,100,-1.,1.);
  TH1F *hSecVtxZTGHCother=new TH1F("hSecVtxZTGHCother","hSecVtxZ_TightCuts_other",1000,-8.,8.);
  TH1F *hSecVtxXTGHCother=new TH1F("hSecVtxXTGHCother","hSecVtxX_TightCuts_other",1000,-3000.,3000.);
  TH1F *hSecVtxYTGHCother=new TH1F("hSecVtxYTGHCother","hSecVtxY_TightCuts_other",1000,-3000.,3000.);
  TH2F *hSecVtxXYTGHCother=new TH2F("hSecVtxXYTGHCother","hSecVtxXY_TightCuts_other",1000,-3000.,3000.,1000,-3000.,3000.);
  TH1F *hSecVtxPhiTGHCother=new TH1F("hSecVtxPhiTGHCother","hSecVtxPhi_TightCuts_other",180,-180.1,180.1);
  TH1F *hd0singlTrackTGHCother=new TH1F("hd0singlTrackTGHCother","hd0singlTrackTightCuts_Other",1000,-5000.,5000.);
  TH1F *hCPtaTGHCother=new TH1F("hCPtaTGHCother","hCPta_TightCuts_other",100,-1.,1.);
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
  flistTghCutsOther->Add(hd0singlTrackTGHCother);
  flistTghCutsOther->Add(hCPtaTGHCother);
  flistTghCutsOther->Add(hd0xd0TGHCother);
  flistTghCutsOther->Add(hMassTrueTGHCother);
  flistTghCutsOther->Add(hMassTGHCother);
  flistTghCutsOther->Add(hMassTrueTGHCotherPM);
  flistTghCutsOther->Add(hMassTGHCotherPM);
  flistTghCutsOther->Add(hMassTrueTGHCotherSB);
  flistTghCutsOther->Add(hMassTGHCotherSB);




 //%%% NEW HISTOS %%%%%%%%%%%%%%%%
  TH1F *hdcaTGHCother=new TH1F("hdcaTGHCother","hdca_TightCuts_Other",100,0.,1000.);
  hdcaTGHCother->SetXTitle("dca   [#mum]");
  hdcaTGHCother->SetYTitle("Entries");
  TH1F *hcosthetastarTGHCother=new TH1F("hcosthetastarTGHCother","hCosThetaStar_TightCuts_Other",50,-1.,1.);
  hcosthetastarTGHCother->SetXTitle("cos #theta^{*}");
  hcosthetastarTGHCother->SetYTitle("Entries");
  TH1F *hptD0TGHCother=new TH1F("hptD0TGHCother","D^{0} transverse momentum distribution",34,ptbinsD0arr);
  hptD0TGHCother->SetXTitle("p_{t}  [GeV/c]");
  hptD0TGHCother->SetYTitle("Entries");
  TH1F *hptD0VsMaxPtTGHCother=new TH1F("hptD0VsMaxPtTGHCother","Difference between D^{0} pt and highest (or second) pt",400,-50.,50.);
  TH2F *hptD0PTallsqrtTGHCother=new TH2F("hptD0PTallsqrtTGHCother","D^{0} pt Vs Sqrt(Sum pt square)",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0PTallTGHCother=new TH2F("hptD0PTallTGHCother","D^{0} pt Vs Sum pt ",34,ptbinsD0arr,200,dumbinning);
  TH2F *hptD0vsptBTGHCother=new TH2F("hptD0vsptBTGHCother","D^{0} pt Vs B pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspBTGHCother=new TH2F("hpD0vspBTGHCother","D^{0} tot momentum Vs B tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hptD0vsptcquarkTGHCother=new TH2F("hptD0vsptcquarkTGHCother","D^{0} pt Vs cquark pt distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  TH2F *hpD0vspcquarkTGHCother=new TH2F("hpD0vspcquarkTGHCother","D^{0} tot momentum Vs cquark tot momentum distribution",34,ptbinsD0arr,34,ptbinsD0arr);
  flistTghCutsOther->Add(hdcaTGHCother);
  flistTghCutsOther->Add(hcosthetastarTGHCother);
  flistTghCutsOther->Add(hptD0TGHCother);
  flistTghCutsOther->Add(hptD0VsMaxPtTGHCother);
  flistTghCutsOther->Add(hptD0PTallsqrtTGHCother);
  flistTghCutsOther->Add(hptD0PTallTGHCother);
  flistTghCutsOther->Add(hptD0vsptBTGHCother);
  flistTghCutsOther->Add(hpD0vspBTGHCother);
  flistTghCutsOther->Add(hptD0vsptcquarkTGHCother);
  flistTghCutsOther->Add(hpD0vspcquarkTGHCother);

  TH1F *hd0zD0ptTGHCother;
  TH1F *hInvMassD0TGHCother,*hInvMassD0barTGHCother;
  TH2F *hInvMassPtTGHCother=new TH2F("hInvMassPtTGHCother","Candidate p_{t} Vs invariant mass",330,1.700,2.030,200,0.,20.);
  THnSparseF *hSparseTGHCother=new THnSparseF("hSparseTGHCother","Candidate Masses, pt, Imp Par;massD0;massD0bar;pt;impactpar;selcase",5,nbinsSparse);
  hSparseTGHCother->SetBinEdges(0,massbins);
  hSparseTGHCother->SetBinEdges(1,massbins);
  hSparseTGHCother->SetBinEdges(2,ptbinsForNsparse);
  hSparseTGHCother->SetBinEdges(3,impparbins);
  hSparseTGHCother->SetBinEdges(4,massHypoBins); 
  flistTghCutsOther->Add(hSparseTGHCother);
  TH1F *hetaTGHCother;
  TH1F *hCosPDPBTGHCother;
  TH1F *hCosPcPDTGHCother;
  flistTghCutsOther->Add(hInvMassPtTGHCother);
// %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
  TH2F *hd0D0VSd0xd0TGHCotherpt;
  TH2F *hangletracksVSd0xd0TGHCotherpt;
  TH2F *hangletracksVSd0D0TGHCotherpt;
  TH1F *hd0xd0TGHCotherpt;

  TH2F *hTOFpidTGHCother=new TH2F("hTOFpidTGHCother","TOF time VS momentum",10,0.,4.,50,-50000.,50000.);
  flistTghCutsOther->Add(hTOFpidTGHCother);

  for(Int_t i=0;i<fnbins;i++){
    namehist="hd0zD0ptTGHCother_pt";
    namehist+=i;
    titlehist="d0(z) Tight Cuts Otherm ptbin=";
    titlehist+=i;
    hd0zD0ptTGHCother=new TH1F(namehist.Data(),titlehist.Data(),1000,-3000,3000.);
    hd0zD0ptTGHCother->SetXTitle("d_{0}(z)    [#mum]");
    hd0zD0ptTGHCother->SetYTitle("Entries");
    flistTghCutsOther->Add(hd0zD0ptTGHCother);

    namehist="hInvMassD0TGHCother_pt";
    namehist+=i;
    titlehist="Invariant Mass Tight Cuts Other ptbin=";
    titlehist+=i;
    hInvMassD0TGHCother=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0TGHCother->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0TGHCother->SetYTitle("Entries");
    flistTghCutsOther->Add(hInvMassD0TGHCother);

    namehist="hInvMassD0barTGHCother_pt";
    namehist+=i;
    titlehist="Invariant Mass D0bar Tight Cuts Other ptbin=";
    titlehist+=i;
    hInvMassD0barTGHCother=new TH1F(namehist.Data(),titlehist.Data(),600,1.600,2.200);
    hInvMassD0barTGHCother->SetXTitle("Invariant Mass    [GeV]");
    hInvMassD0barTGHCother->SetYTitle("Entries");
    flistTghCutsOther->Add(hInvMassD0barTGHCother);

    namehist="hetaTGHCother_pt";
    namehist+=i;
    titlehist="eta Tight Cuts Other ptbin=";
    titlehist+=i;
    hetaTGHCother=new TH1F(namehist.Data(),titlehist.Data(),100,-3.,3.);
    hetaTGHCother->SetXTitle("Pseudorapidity");
    hetaTGHCother->SetYTitle("Entries");
    flistTghCutsOther->Add(hetaTGHCother);

    namehist="hCosPDPBTGHCother_pt";
    namehist+=i;
    titlehist="Cosine between D0 momentum and B momentum, ptbin=";
    titlehist+=i;
    hCosPDPBTGHCother=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPDPBTGHCother->SetXTitle("Cosine between D0 momentum and B momentum");
    hCosPDPBTGHCother->SetYTitle("Entries");
    flistTghCutsOther->Add(hCosPDPBTGHCother);

    namehist="hCosPcPDTGHCother_pt";
    namehist+=i;
    titlehist="Cosine between cquark momentum and D0 momentum, ptbin=";
    titlehist+=i;
    hCosPcPDTGHCother=new TH1F(namehist.Data(),titlehist.Data(),50,-1.,1.);
    hCosPcPDTGHCother->SetXTitle("Cosine between c quark momentum and D0 momentum");
    hCosPcPDTGHCother->SetYTitle("Entries");
    flistTghCutsOther->Add(hCosPcPDTGHCother);

// %%%%%%%%%%% HERE THE ADDITIONAL HISTS %%%%%%%%%%%%%%%%%
    namehist="hd0xd0TGHCother_pt";
    namehist+=i;
    titlehist="d0xd0 Tight Cuts Other ptbin=";
    titlehist+=i;
    hd0xd0TGHCotherpt=new TH1F(namehist.Data(),titlehist.Data(),1000,-50000.,10000.);
    hd0xd0TGHCotherpt->SetXTitle("d_{0}^{K}xd_{0}^{#pi}   [#mum^2]");
    hd0xd0TGHCotherpt->SetYTitle("Entries");
    flistTghCutsOther->Add(hd0xd0TGHCotherpt);


    namehist="hd0D0VSd0xd0TGHCother_pt";
    namehist+=i;
    titlehist="d_{0}^{D^{0}} Vs d_{0}^{K}xd_{0}^{#pi} Tight Cuts Other ptbin=";
    titlehist+=i;
    hd0D0VSd0xd0TGHCotherpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,100,-300,300);
    hd0D0VSd0xd0TGHCotherpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hd0D0VSd0xd0TGHCotherpt->SetYTitle(" d_{0}^{D^{0}}    [#mum]");
    flistTghCutsOther->Add(hd0D0VSd0xd0TGHCotherpt);
    
    
    namehist="hangletracksVSd0xd0TGHCother_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{K}xd_{0}^{#pi} Tight Cuts Other ptbin=";
    titlehist+=i;
    hangletracksVSd0xd0TGHCotherpt=new TH2F(namehist.Data(),titlehist.Data(),200,-50000.,30000.,40,-0.1,3.24);
    hangletracksVSd0xd0TGHCotherpt->SetXTitle(" d_{0}^{K}xd_{0}^{#pi} [#mum]");
    hangletracksVSd0xd0TGHCotherpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistTghCutsOther->Add(hangletracksVSd0xd0TGHCotherpt);
    

    namehist="hangletracksVSd0D0TGHCother_pt";
    namehist+=i;
    titlehist="Angle between K and #pi tracks Vs d_{0}^{D^{0}} Tight Cuts Other ptbin=";
    titlehist+=i;
    hangletracksVSd0D0TGHCotherpt=new TH2F(namehist.Data(),titlehist.Data(),200,-400.,400.,40,-0.12,3.24);
    hangletracksVSd0D0TGHCotherpt->SetXTitle(" d_{0}^{D^{0}} [#mum]");
    hangletracksVSd0D0TGHCotherpt->SetYTitle(" angle between K and #p tracks  [rad]");
    flistTghCutsOther->Add(hangletracksVSd0D0TGHCotherpt);
    
  }
  // %%%%%%%% END OF NEW HISTOS %%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





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
  
  TH1F *hd0D0ptTGHCotherPM;
  TH1F *hMCd0D0ptTGHCotherPM;
  TH1F *hd0D0VtxTrueptTGHCotherPM;
  TH1F *hd0D0ptTGHCotherSB;
  TH1F *hMCd0D0ptTGHCotherSB;
  TH1F *hd0D0VtxTrueptTGHCotherSB;
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
    
    hd0D0ptTGHCotherPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCotherPM->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCotherPM->SetYTitle("Entries");
    flistTghCutsOther->Add(hd0D0ptTGHCotherPM);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCotherPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCotherPM->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCotherPM->SetYTitle("Entries");
    flistTghCutsOther->Add(hMCd0D0ptTGHCotherPM);
 

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCotherPM = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCotherPM->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCotherPM->SetYTitle("Entries");
    flistTghCutsOther->Add(hd0D0VtxTrueptTGHCotherPM);
    
    strnamept=namehist;
    strnamept.Append("SBMss_pt");
    strnamept+=i;

    strtitlept=titlehist;
    strtitlept.Append(" Side Bands, ");
    strtitlept+=fptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=fptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    hd0D0ptTGHCotherSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0ptTGHCotherSB->SetXTitle("Impact parameter [#mum] ");
    hd0D0ptTGHCotherSB->SetYTitle("Entries");
    flistTghCutsOther->Add(hd0D0ptTGHCotherSB);

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    hMCd0D0ptTGHCotherSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hMCd0D0ptTGHCotherSB->SetXTitle("MC Impact parameter [#mum] ");
    hMCd0D0ptTGHCotherSB->SetYTitle("Entries");
    flistTghCutsOther->Add(hMCd0D0ptTGHCotherSB);

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    hd0D0VtxTrueptTGHCotherSB = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    hd0D0VtxTrueptTGHCotherSB->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    hd0D0VtxTrueptTGHCotherSB->SetYTitle("Entries");
    flistTghCutsOther->Add(hd0D0VtxTrueptTGHCotherSB);
  }
  Printf("AFTER DATA HISTOS CREATION \n");

  delete ptbinlimitsCxyLxy;


  PostData(1,fNentries);
  PostData(2,fSignalType);
  PostData(3,fSignalTypeLsCuts);
  PostData(4,fSignalTypeTghCuts);
  PostData(5,fCounter);
  PostData(6,flistMCproperties);
  PostData(7,flistNoCutsSignal);
  PostData(8,flistNoCutsBack);
  PostData(9,flistNoCutsFromB);
  PostData(10,flistNoCutsFromDstar);
  PostData(11,flistNoCutsOther);
  PostData(12,flistLsCutsSignal);
  PostData(13,flistLsCutsBack);
  PostData(14,flistLsCutsFromB);
  PostData(15,flistLsCutsFromDstar);
  PostData(16,flistLsCutsOther);
  PostData(17,flistTghCutsSignal);
  PostData(18,flistTghCutsBack);
  PostData(19,flistTghCutsFromB);
  PostData(20,flistTghCutsFromDstar);
  PostData(21,flistTghCutsOther);

  return;

}





//________________________________________________________________________
void AliAnalysisTaskSECharmFraction::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if (!aod) {
    Printf("ERROR: aod not available");
    return;
  }
  TClonesArray *arrayD0toKpi;
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
      AliAODEvent* aodFromExt = ext->GetAOD();
      if(fLikeSign){
	// load 2Prong Like Sign                                                   
	arrayD0toKpi =(TClonesArray*)aodFromExt->GetList()->FindObject("LikeSign2Prong");
	if(!arrayD0toKpi) {
	  Printf("AliAnalysisTaskSECharmFraction::UserExec: LikeSign branch not found!\n");
	  return;
	}
      }
      else {
	// load D0->Kpi candidates                                                   
	arrayD0toKpi = (TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
	if(!arrayD0toKpi) {
	  Printf("AliAnalysisTaskSECharmFraction::UserExec: D0toKpi branch not found!\n");
	  return;
	}
      }  
    }
  } else {
    if(fLikeSign){
      // load 2Prong Like Sign                                                   
      arrayD0toKpi =(TClonesArray*)aod->GetList()->FindObject("LikeSign2Prong");
      if(!arrayD0toKpi) {
	Printf("AliAnalysisTaskSECharmFraction::UserExec: LikeSign branch not found!\n");
	return;
      }
    }
    else {
      // load D0->Kpi candidates                                                   
      arrayD0toKpi = (TClonesArray*)aod->GetList()->FindObject("D0toKpi");
      if(!arrayD0toKpi) {
	Printf("AliAnalysisTaskSECharmFraction::UserExec: D0toKpi branch not found!\n");
	return;
      }
    }     
  }


  if(!arrayD0toKpi) {
    printf("AliAnalysisTaskSECharmFraction::UserExec: input branch not found!\n");
    return;
  }

  //histogram filled with 1 for every AOD
  fNentries->Fill(0);
  fCounter->StoreEvent(aod,fCutsLoose,fReadMC); 
  
  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  //  TString trigclass=aod->GetFiredTriggerClasses();
  // if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNentries->Fill(14);
 
  Int_t nSelectedloose=0, nSelectedtight=0; 

  Bool_t isEventSelTGHT=kTRUE,isEventSelLOOSE=kTRUE;
  if(!fCutsTight->IsEventSelected(aod)){
    isEventSelTGHT=kFALSE;
    if(fCutsTight->GetWhyRejection()==1){ 
      // rejected for pileup
      fNentries->Fill(2);    
    }
    if(fCutsTight->GetWhyRejection()==6){
      // |prim Vtx Zspd| > acceptable
      fNentries->Fill(4);      
    }
  }
  else {
    fNentries->Fill(1);    
  }
  if(!fCutsLoose->IsEventSelected(aod)){
    isEventSelLOOSE=kFALSE;
    if(fCutsLoose->GetWhyRejection()==1){
      // rejected for pileup    
      fNentries->Fill(9);
      
    }
    if(fCutsLoose->GetWhyRejection()==6){
      // |prim Vtx Z| > acceptable
      fNentries->Fill(11);      
    }
  }
  else {
    fNentries->Fill(8);    
  }
  
  if(!(isEventSelTGHT||isEventSelLOOSE)){
      PostData(1,fNentries);
      return;
  }
  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001){
    if(isEventSelTGHT)fNentries->Fill(19);
    if(isEventSelLOOSE)fNentries->Fill(20);
    PostData(1,fNentries);
    return;
  }
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) { 
    
    if(isEventSelTGHT)fNentries->Fill(3);
    if(isEventSelLOOSE)fNentries->Fill(10);
  }
  else {
    PostData(1,fNentries);
    return;
    
  }

  // FILL n-tracks counter
  if(isEventSelTGHT)fNentries->Fill(5,aod->GetNumberOfTracks());
  if(isEventSelLOOSE)fNentries->Fill(12,aod->GetNumberOfTracks());
  

  Bool_t aziListIsFilled=kFALSE;
  Double_t azilist[30000];
  Int_t trkIDlist[30000],nprim=0;
  
 
  for(Int_t ephi=0;ephi<30000;ephi++){
    azilist[ephi]=-999.;
    trkIDlist[ephi]=-999;
  }
  //aziListIsFilled=kFALSE;
    



  
  TClonesArray *arrayMC=0x0;
  AliAODMCHeader *aodmcHeader=0x0;
  Double_t vtxTrue[3];
 
  
  if(fReadMC){
    // load MC particles
    arrayMC = 
      (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      Printf("AliAnalysisTaskSECharmFraction::UserExec: MC particles branch not found!\n");
      return;
    }
    // load MC header
    aodmcHeader = 
      (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!aodmcHeader) {
      Printf("AliAnalysisTaskSECharmFraction::UserExec: MC header branch not found!\n");
      return;
    }
    // MC primary vertex
    aodmcHeader->GetVertex(vtxTrue);
    // FILL HISTOS FOR D0 mesons, c quarks and D0 from B properties
    FillHistoMCproperties(arrayMC);
  }
 


  //Printf("There are %d tracks in this event", aod->GetNumberOfTracks());
  //  Int_t nTotHF=0,nTotDstar=0,nTot3Prong=0;
  Int_t nTotD0toKpi=0;
  Int_t okd0tight,okd0bartight,okd0loose,okd0barloose,okd0tightnopid,okd0bartightnopid;
  Bool_t defaultNC=kTRUE;
  Bool_t isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,isSideBand;
  Bool_t isinacceptance;
  Int_t signallevel=-1;
  Int_t ptbin,nVtx; 
  //  const  Int_t nptbins=10;
  Double_t invMassD0,invMassD0bar,ptD0,massmumtrue;
  
 
  AliAODRecoDecayHF *aodDMC=0x0;// to be used to create a fake true sec vertex
  // make trkIDtoEntry register (temporary)  
  
  if(fFastAnalysis<1){
    Int_t trkIDtoEntry[100000];
    fptAll=0.;
    fptAllSq=0.;
    fptMax[0]=0.;
    fptMax[1]=0.;
    fptMax[2]=0.;
    for(Int_t it=0;it<aod->GetNumberOfTracks();it++) {
      AliAODTrack *track = aod->GetTrack(it);
      fptAll+=track->Pt();
      fptAllSq+=track->Pt()*track->Pt();
      if(track->Pt()>fptMax[0]){
	fptMax[2]=fptMax[1];
	fptMax[1]=fptMax[0];
	fptMax[0]=track->Pt();
      }
      else if(track->Pt()>fptMax[1]){
	fptMax[2]=fptMax[1];
	fptMax[1]=track->Pt();
      }
      else if(track->Pt()>fptMax[2])fptMax[2]=track->Pt();
      if(track->GetID()<0) {
	if(isEventSelTGHT)fNentries->Fill(19);
	if(isEventSelLOOSE)fNentries->Fill(20);
	//printf("Track ID <0, id= %d\n",track->GetID());
	PostData(1,fNentries);
	return;
      }
      trkIDtoEntry[track->GetID()]=it;
    }
  }

  // loop over D0->Kpi candidates
  Int_t nD0toKpi = arrayD0toKpi->GetEntriesFast();
  nTotD0toKpi += nD0toKpi;
  // fille D0 candidate counter 
  if(isEventSelTGHT)fNentries->Fill(6,nD0toKpi);
  if(isEventSelLOOSE)fNentries->Fill(13,nD0toKpi);
  
  //	cout<<"Number of D0->Kpi: "<<nD0toKpi<<endl;
  
  for (Int_t iD0toKpi = 0; iD0toKpi < nD0toKpi; iD0toKpi++) {
    if(aodDMC!=0x0)delete aodDMC;
      
    defaultNC=kTRUE;
    isPeakD0=kFALSE;
    isPeakD0bar=kFALSE;
    isSideBandD0=kFALSE;
    isSideBandD0bar=kFALSE;
    isSideBand=kFALSE;
    isinacceptance=kFALSE;
    okd0tight=0;
    okd0bartight=0;
    okd0tightnopid=0;
    okd0bartightnopid=0;
    okd0loose=0;
    okd0barloose=0;
  
    signallevel=-1;
    

    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->UncheckedAt(iD0toKpi);
   //  Bool_t unsetvtx=kFALSE;
//     if(!d->GetOwnPrimaryVtx()) {
//       d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
//       unsetvtx=kTRUE;
//     }

    //recalculate vertex w/o daughters
    AliAODVertex *origownvtx=0x0;
    if(fCleanCandOwnVtx){      
      if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
      if(!fCutsTight->RecalcOwnPrimaryVtx(d,aod))  defaultNC=kFALSE;   
      
    }

  
    // ############ MISALIGN HERE: TEMPORARY SOLUTION ##########
    //    d->Misalign("resC");

    
    //############# SIGNALLEVEL DESCRIPTION #####################
    // TO BE CONSIDERED WITH GREAT CARE, only =0 and =1 (and MC selection when possible) are reliable
    //                                    For the other signallevel numbers the order in which cut are applied is relevant 
    // signallevel =0,1: is selected as signal,is signal (MC)
    //               from 2 to 20: MC information
    //              from 21 to 29: "detector" selection (acceptance, pt, refits) 
    //                                                  21: acceptance, eta (N.B. before 24 May was signallevel=9)     
    //                                                  22: isinfiducialacceptance
    //                                                  23: single track p
    //                                                  25: ITS cluster selection
    //                                                  26: TPC refit
    //                                                  27: ITS refit
    //                                                  28: no (TOF||TPC) pid information (no kTOFpid,kTOFout,kTIME,kTPCpid,...)
    //
    //              from 30 to 39: PID selection
    //                                                  31: no Kaon compatible tracks found between daughters
    //                                                  32: no Kaon identified tracks found (strong sel. at low momenta)
    //                                                  33: both mass hypotheses are rejected 
    //              from 40 to 45: standard cut selection
    //              from 45 to 49: special cut signal kinematic selection
    //                                                  46: pstar cut
    //              from 50 to 60: special cut selection
    //                                                  51: Nvtx contributors
    //                                                  52: angle between tracks
    //                                                  53: vtx not reconstructed when excludind daughters 
    //                                                  54: track not propagated to dca when the vtx is recalculated
    //                                                  55: single track normalized impact par.
    //                                                  56: normalized d0xd0 
    //                                                  57: d0xd0 cut with vtx on the fly
    //                                                  58,59: cut normalized decay lenght and decay lenght
    //####### DATA SELECTION ####################################
    //
    // ######## CHECK FOR ACCEPTANCE ##########
    ptD0=d->Pt();
    ptbin=fCutsTight->PtBin(ptD0);
    //    Double_t relangle=d->ProngsRelAngle(0,1);
    // UPV: HERE TO CHANGE WITH:
     //  isinacceptance = (TMath::Abs(d->EtaProng(0))<fAcceptanceCuts[0]&&TMath::Abs(d->EtaProng(1))<fAcceptanceCuts[0]); //eta acceptance 

    // INVESTIGATE SIGNAL TYPE : ACCESS TO MC INFORMATION
    if(fReadMC){
      aodDMC=GetD0toKPiSignalType(d,arrayMC,signallevel,massmumtrue,vtxTrue);
    }
    else signallevel=0;
    //   ACCOUNT FOR ETA & ITS CLUSTER SELECTION
    if(fFastAnalysis<1){  // ALREADY DONE BY AliRDHFCutsD0ToKPi selection
      isinacceptance=fCutsTight->AreDaughtersSelected(d); 
      if(!isinacceptance)signallevel=21;
    }
    if(!fCutsTight->IsInFiducialAcceptance(ptD0,d->Y(421))){
      isinacceptance=kFALSE;
      signallevel=22; 
    }
    else{
      nSelectedloose++;
    }
   
    //###################################################################################
    //
    // ################ SPECIAL ADDITIONAL CUTS FOR SIGNAL SEARCH #######################
    //  UPV: ITS CLUSTER SELECTION ALREADY DONE , COUNT THEM THE SAME
    //  
    Int_t nlayers=0,nSPD=0,nSSD=0;
    Bool_t spd1=kFALSE;
    if(fFastAnalysis<1){
      for(Int_t idg=0;idg<d->GetNDaughters();idg++){
	
	AliAODTrack *dgTrack = (AliAODTrack*)d->GetDaughter(idg);
	if(TESTBIT(dgTrack->GetITSClusterMap(),5)){
	  nlayers++;
	  nSSD++;
	}
      if(TESTBIT(dgTrack->GetITSClusterMap(),4)){
	nlayers++;
	nSSD++;
      }
      if(TESTBIT(dgTrack->GetITSClusterMap(),3))nlayers++;
      if(TESTBIT(dgTrack->GetITSClusterMap(),2))nlayers++;
      if(TESTBIT(dgTrack->GetITSClusterMap(),1)){
	nlayers++;
	nSPD++;
      }
      if(TESTBIT(dgTrack->GetITSClusterMap(),0)){
	nlayers++;
	nSPD++;
	spd1=kTRUE;
      }
      }
    }
    /*
      // ######## NOW SELECTION ##########
      if(dgTrack->Pt()<0.5){
	// ########## k-Both selection ##############
	if(nlayers<5)signallevel=25;
	if(nSPD<2)signallevel=25;
      }
      else if(dgTrack->Pt()<1.){
	// ########## 4 its clust (1 SSD,1 SPD) && k-Any selection ##############
	if(nlayers<4)signallevel=25;
	if(nSSD<1)signallevel=25;
	if(nSPD<1)signallevel=25;
      }
      else{ 
	// ########## 3 its clust (1 SPD, 1 SSD) && k-Any selection ##########
	if(nlayers<3)signallevel=25;
	if(nSSD<1)signallevel=25;
	if(nSPD<1)signallevel=25;	    
      }
    }
  */


   

    //###########   END OF SPECIAL CUTS        ######################
    //
    //###############################################################

    // NOW APPLY CUTS
    // Check tighter cuts w/o PID:
    // 
    //    Int_t ncont=vtx1->GetNContributors();   
    if(fFastAnalysis<1)if(vtx1->GetNContributors()<1)signallevel=51;
  
    if(defaultNC&&fFastAnalysis<1&&fNtrMaxforVtx<2)defaultNC=SpecialSelD0(d,nVtx);// No More in USE!
    if(isEventSelTGHT&&defaultNC){
      Bool_t iscutusingpid=fCutsTight->GetIsUsePID();
      fCutsTight->SetUsePID(kFALSE);
      Int_t isSelectedTightNoPid=fCutsTight->IsSelected(d,AliRDHFCuts::kAll,aod);
      switch(isSelectedTightNoPid){
      case 0:
	okd0tightnopid=kFALSE;
	okd0bartightnopid=kFALSE;
	break;
      case 1:
	okd0tightnopid=kTRUE;
	okd0bartightnopid=kFALSE;
	break;
      case 2:
	okd0tightnopid=kFALSE;
	okd0bartightnopid=kTRUE;
	break;
      case 3:
	okd0tightnopid=kTRUE;
	okd0bartightnopid=kTRUE;
	break;
      default:
	okd0tightnopid=kTRUE;
	okd0bartightnopid=kTRUE;
	break;
      }
      
    
      //    signallevel=fCutsTight->GetSelectionStep();
      fSignalType->Fill(signallevel);
      
    
      
      
      // ######### SPECIAL SELECTION PID ##############
      fCutsTight->SetUsePID(iscutusingpid);
      if(okd0tightnopid||okd0bartightnopid){
	Int_t isSelectedPid=fCutsTight->IsSelected(d,AliRDHFCuts::kPID,aod);
	Int_t isSelectedTight=fCutsTight->CombineSelectionLevels(3,isSelectedTightNoPid,isSelectedPid);
	switch(isSelectedTight){
	case 0:
	  okd0tight=kFALSE;
	  okd0bartight=kFALSE;
	  break;
	case 1:
	  okd0tight=kTRUE;
	  okd0bartight=kFALSE;
	  break;
	case 2:
	  okd0tight=kFALSE;
	  okd0bartight=kTRUE;
	  break;
	case 3:
	  okd0tight=kTRUE;
	  okd0bartight=kTRUE;
	  break;
	default:
	  okd0tight=kTRUE;
	  okd0bartight=kTRUE;
	  break;
	}
      }
    }
    else{
      fSignalType->Fill(signallevel);
    }

    
  
    
    if(isEventSelLOOSE&&defaultNC){
      // CHECK LOOSER CUTS 
      //ptbin=SetStandardCuts(ptD0,flargeInvMassCut);
      
      //    d->SelectD0(fCutsLoose->GetD0toKpiCuts(),okd0loose,okd0barloose);
      Int_t isSelectedLoose=fCutsLoose->IsSelected(d,AliRDHFCuts::kAll,aod);
      switch(isSelectedLoose){
      case 0:
	okd0loose=kFALSE;
	okd0barloose=kFALSE;
	break;
      case 1:
	okd0loose=kTRUE;
	okd0barloose=kFALSE;
	break;
      case 2:
	okd0loose=kFALSE;
	okd0barloose=kTRUE;
      break;
      case 3:
	okd0loose=kTRUE;
	okd0barloose=kTRUE;
	break;
      default:
	okd0loose=kTRUE;
	okd0barloose=kTRUE;
	break;
      }
     
    }
    
  
    //NO CUTS Case: force okD0 and okD0bar = kTRUE
    //     special cuts are applied also in the "NO Cuts" case
    //
    //
    // SPECIAL modification:
    // IMPORTANT!!!!!!  ONLY FOR TEMPORARY CONVENIENCE
    // IF fusePID is kTRUE, NO CUTS BECOMES NO PID case with tight cuts (fill signal histos when 30<=signallevel<40 )!!!!!!!!!!  
    if((!fusePID)&&isEventSelTGHT){
      okd0tightnopid=defaultNC;
      okd0bartightnopid=defaultNC;
    }        

    if(okd0loose||okd0barloose||okd0tight||okd0bartight||okd0tightnopid||okd0bartightnopid){
      //######## INVARIANT MASS SELECTION ###############
      CheckInvMassD0(d,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar);
      if((isSideBandD0||isSideBandD0bar)){
	if(!(isPeakD0||isPeakD0bar))isSideBand=kTRUE; //(isSideBand no more used in the following, can remove it)
	else {// AVOID TO CONSIDER IN THE SIDEBAND  THOSE  CANDIDATES FOR WHICH 1 MASS HYPOTHESIS IS IN THE PEAK REGION
	  isSideBand=kFALSE; 
	  isSideBandD0=kFALSE;
	  isSideBandD0bar=kFALSE;	   
	}
      }
      if(fFastAnalysis<2){
	if(!aziListIsFilled){
	  FillAziList(aod,azilist,trkIDlist,nprim);
	  aziListIsFilled=kTRUE;
	}
	
	if(signallevel==1||signallevel==0){
	  if(nprim!=0){
	    FillAziHistos(d,flistNoCutsSignal,ptbin,azilist,trkIDlist,nprim,okd0tightnopid,okd0bartightnopid,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar); 
	    FillAziHistos(d,flistTghCutsSignal,ptbin,azilist,trkIDlist,nprim,okd0tight,okd0bartight,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar); 
	    FillAziHistos(d,flistLsCutsSignal,ptbin,azilist,trkIDlist,nprim,okd0loose,okd0barloose,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar); 	  
	  }
	}
	
      }
    }
    if(((isPeakD0&&okd0tight)||(isPeakD0bar&&okd0bartight))&&isinacceptance)fSignalTypeTghCuts->Fill(signallevel);
    if(((isPeakD0&&okd0loose)||(isPeakD0bar&&okd0barloose))&&isinacceptance)fSignalTypeLsCuts->Fill(signallevel);
  
   
    
  
    //###################    FILL HISTOS      ########################
    //################################################################
    //
    //######## improvement: SPEED HERE CAN BE IMPROVED: CALCULATE ONCE AND FOR ALL 
    //            CANDIDATE VARIABLES   

   
    if(signallevel==1||signallevel==0)FillHistos(d,flistNoCutsSignal,ptbin,okd0tightnopid,okd0bartightnopid,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);    // else if(fusePID&&signallevel>=30&&signallevel<40)FillHistos(d,flistNoCutsSignal,ptbin,okd0tightnopid,okd0bartightnopid,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);// OLD LINE, COULD BE REMOVED 
    else if(signallevel==2)FillHistos(d,flistNoCutsFromDstar,ptbin,okd0tightnopid,okd0bartightnopid,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==3||signallevel==4)FillHistos(d,flistNoCutsFromB,ptbin,okd0tightnopid,okd0bartightnopid,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==-1||signallevel==7||signallevel==8||signallevel==10||signallevel==9||signallevel==11)FillHistos(d,flistNoCutsBack,ptbin,okd0tightnopid,okd0bartightnopid,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==5||signallevel==6)FillHistos(d,flistNoCutsOther,ptbin,okd0tightnopid,okd0bartightnopid,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    

    
    //LOOSE CUTS Case
    if(okd0loose||okd0barloose)fNentries->Fill(14);

    if(signallevel==1||signallevel==0)FillHistos(d,flistLsCutsSignal,ptbin,okd0loose,okd0barloose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==2)FillHistos(d,flistLsCutsFromDstar,ptbin,okd0loose,okd0barloose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==3||signallevel==4)FillHistos(d,flistLsCutsFromB,ptbin,okd0loose,okd0barloose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==-1||signallevel==7||signallevel==8||signallevel==10||signallevel==11)FillHistos(d,flistLsCutsBack,ptbin,okd0loose,okd0barloose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==5||signallevel==6)FillHistos(d,flistLsCutsOther,ptbin,okd0loose,okd0barloose,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);

    //TIGHT CUTS Case
    if(okd0tight||okd0bartight){
      fNentries->Fill(7);
      nSelectedtight++; 
    }
    
    if(signallevel==1||signallevel==0)FillHistos(d,flistTghCutsSignal,ptbin,okd0tight,okd0bartight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==2)FillHistos(d,flistTghCutsFromDstar,ptbin,okd0tight,okd0bartight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==3||signallevel==4)FillHistos(d,flistTghCutsFromB,ptbin,okd0tight,okd0bartight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==-1||signallevel==7||signallevel==8||signallevel==10||signallevel==11)FillHistos(d,flistTghCutsBack,ptbin,okd0tight,okd0bartight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    else if(signallevel==5||signallevel==6)FillHistos(d,flistTghCutsOther,ptbin,okd0tight,okd0bartight,invMassD0,invMassD0bar,isPeakD0,isPeakD0bar,isSideBandD0,isSideBandD0bar,massmumtrue,aodDMC,vtxTrue);
    

    // ######## PRINTING INFO FOR D0-like candidate

    if(nSPD==2&&ptD0>2.){
      if((okd0tight&&TMath::Abs(invMassD0-1.864)<0.01)||(okd0bartight&&TMath::Abs(invMassD0bar-1.864)<0.01)){
	//printf("INFO FOR DRAWING: \n pt: %f \n Rapidity: %f \n Period Number: %d \n Run Number: %d \n BunchCrossNumb: %d \n OrbitNumber: %d \n",ptD0,d->Y(421),aod->GetPeriodNumber(),aod->GetRunNumber(),aod->GetBunchCrossNumber(),aod->GetOrbitNumber());
	//printf("PrimVtx NContributors: %d \n Prongs Rel Angle: %f \n \n",ncont,relangle);
      }
    }
    if(aodDMC!=0x0){
      delete aodDMC;
      aodDMC=0x0;
    }
    
    if(fCleanCandOwnVtx)fCutsTight->CleanOwnPrimaryVtx(d,aod,origownvtx);
    

    //   if(unsetvtx) d->UnsetOwnPrimaryVtx();
    
  }


  fCounter->StoreCandidates(aod,nSelectedloose,kTRUE);  
  fCounter->StoreCandidates(aod,nSelectedtight,kFALSE); 
  
  
  // ####################### POST OUTPUT TLIST DATA #########################
  // ####### histo for #AOD entries already posted
  PostData(1,fNentries);
  PostData(2,fSignalType);
  PostData(3,fSignalTypeLsCuts);
  PostData(4,fSignalTypeTghCuts);
  PostData(5,fCounter);
  PostData(6,flistMCproperties);
  PostData(7,flistNoCutsSignal);
  PostData(8,flistNoCutsBack);
  PostData(9,flistNoCutsFromB);
  PostData(10,flistNoCutsFromDstar);
  PostData(11,flistNoCutsOther);
  PostData(12,flistLsCutsSignal);
  PostData(13,flistLsCutsBack);
  PostData(14,flistLsCutsFromB);
  PostData(15,flistLsCutsFromDstar);
  PostData(16,flistLsCutsOther);
  PostData(17,flistTghCutsSignal);
  PostData(18,flistTghCutsBack);
  PostData(19,flistTghCutsFromB);
  PostData(20,flistTghCutsFromDstar);
  PostData(21,flistTghCutsOther);

  return;
}



//_________________________________________
Int_t  AliAnalysisTaskSECharmFraction::SetStandardCuts(Float_t *&ptbinlimits){
  //
  // creating cuts for D0 -> Kpi
  //

  Printf("Using Default Cuts as set in AliAnalysisTaskSECharmFraction \n");
  //  const Double_t ptmin = 0.1;
  const Double_t ptmax = 9999.;
  const Int_t nptbins =13;
  const Int_t nvars=9;
  Int_t varycuts=-1;
  
  if(fCutsTight){
    delete fCutsTight;fCutsTight=NULL;
  }
  if(fCutsLoose){
    delete fCutsLoose;fCutsLoose=NULL;
  }


  fCutsTight = new AliRDHFCutsD0toKpi();
  fCutsTight->SetName("D0toKpiCutsStandard");
  fCutsTight->SetTitle("Standard Cuts for D0 analysis");
  
  fCutsLoose = new AliRDHFCutsD0toKpi();
  fCutsLoose->SetName("D0toKpiCutsLoose");
  fCutsLoose->SetTitle("Loose Cuts for D0 analysis");
  
  // EVENT CUTS
  fCutsTight->SetMinVtxContr(1);
  fCutsLoose->SetMinVtxContr(1);
  
  // TRACKS ON SINGLE TRACKS
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //  esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  

  fCutsTight->AddTrackCuts(esdTrackCuts);
  fCutsLoose->AddTrackCuts(esdTrackCuts);
  
  

  Float_t ptbins[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=0.5;	
  ptbins[2]=1.;
  ptbins[3]=2.;
  ptbins[4]=3.;
  ptbins[5]=4.;
  ptbins[6]=5.;
  ptbins[7]=6.;
  ptbins[8]=8.;
  ptbins[9]=12.;
  ptbins[10]=16.;
  ptbins[11]=20.;
  ptbins[12]=24.;
  ptbins[13]=ptmax;

  fCutsTight->SetGlobalIndex(nvars,nptbins);
  fCutsLoose->SetGlobalIndex(nvars,nptbins);
  fCutsTight->SetPtBins(nptbins+1,ptbins);
  fCutsLoose->SetPtBins(nptbins+1,ptbins);
  
  /*	Float_t cutsArrayD0toKpiStand_1[9]={0.200,300.*1E-4,0.8,0.3,0.3,1000.*1E-4,1000.*1E-4,-35000.*1E-8,0.7};   // pt<1 
		Float_t cutsArrayD0toKpiStand_2[9]={0.200,200.*1E-4,0.8,0.4,0.4,1000.*1E-4,1000.*1E-4,-30000.*1E-8,0.8}; // 1<=pt<2 
		Float_t cutsArrayD0toKpiStand_3[9]={0.200,200.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-26000.*1E-8,0.94};   // 2<=pt<3 
		Float_t cutsArrayD0toKpiStand_4[9]={0.200,200.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-15000.*1E-8,0.88};   // 3<=pt<5 
		Float_t cutsArrayD0toKpiStand_5[9]={0.200,150.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-10000.*1E-8,0.9};   // 5<=pt<8
		Float_t cutsArrayD0toKpiStand_6[9]={0.200,150.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-10000.*1E-8,0.9};   // 8<pt<12
		Float_t cutsArrayD0toKpiStand_7[9]={0.200,150.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-10000.*1E-8,0.9}; // pt>12
	*/

	const Int_t nvary=3;
	Float_t varyd0xd0[nptbins][nvary]={{-35000.*1E-8,-40000.*1E-8,-50000.*1E-8},/* pt<0.5*/
					   {-35000.*1E-8,-40000.*1E-8,-50000.*1E-8},/* 0.5<pt<1*/
					   {-25000.*1E-8,-32000.*1E-8,-38000.*1E-8},/* 1<pt<2 */
					   {-22000.*1E-8,-26000.*1E-8,-30000.*1E-8},/* 2<pt<3 */
					   {-12000.*1E-8,-15000.*1E-8,-20000.*1E-8},/* 3<pt<4 */
					   {-12000.*1E-8,-15000.*1E-8,-20000.*1E-8},/* 4<pt<5 */
					   {-5000.*1E-8,-10000.*1E-8,-15000.*1E-8},/* 5<pt<6 */
					   {-5000.*1E-8,-10000.*1E-8,-15000.*1E-8},/* 6<pt<8 */
					   {-0.*1E-8,-10000.*1E-8,-12000.*1E-8},/* 8<pt<12 */
					   {5000.*1E-8,-5000.*1E-8,-10000.*1E-8},/* 12<pt<16 */
					   {5000.*1E-8,-5000.*1E-8,-10000.*1E-8},/* 16<pt<20 */
					   {5000.*1E-8,-5000.*1E-8,-10000.*1E-8},/* 20<pt<24 */
					   {5000.*1E-8,-5000.*1E-8,-10000.*1E-8}};/* pt>24 */
       

	Float_t varyCosPoint[nptbins][nvary]={{0.75,0.80,0.85},/* 0<pt<0.5 */
					      {0.75,0.80,0.85},/* 0.5<pt<1*/
					      {0.75,0.80,0.85},/* 1<pt<2 */
					      {0.92,0.94,0.95},/* 2<pt<3 */
					      {0.85,0.88,0.91},/* 3<pt<4 */
					      {0.85,0.88,0.91},/* 4<pt<5 */
					      {0.88,0.90,0.92},/* 5<pt<6 */
					      {0.88,0.90,0.92},/* 6<pt<8 */
					      {0.85,0.90,0.92},/* 8<pt<12 */
					      {0.85,0.90,0.92},/* 12<pt<16 */
					      {0.8,0.85,0.9},/* 16<pt<20 */
					      {0.8,0.85,0.9},/* 20<pt<24 */
					      {0.75,0.82,0.9}};/* pt>24 */
	

	
	if(varycuts==-1){//DEFAULT CUTS
	  varycuts=11;
	  varyd0xd0[9][1]=-10000.*1E-8;
	  varyd0xd0[10][1]=-10000.*1E-8;
	  varyd0xd0[11][1]=-10000.*1E-8;	  
	  varyd0xd0[12][1]=-10000.*1E-8;
	}
	Int_t vcd0xd0=varycuts/10;
	Int_t vccospoint=varycuts%10;
	// ######################## STAND VARY CUTS  ###########################################	
	Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.400,300.*1E-4,0.8,0.3,0.3,1000.*1E-4,1000.*1E-4,varyd0xd0[0][vcd0xd0],varyCosPoint[0][vccospoint]},/* 0<pt<0.5*/
							{0.400,300.*1E-4,0.8,0.3,0.3,1000.*1E-4,1000.*1E-4,varyd0xd0[1][vcd0xd0],varyCosPoint[1][vccospoint]},/* 0.5<pt<1*/
							{0.400,300.*1E-4,0.8,0.4,0.4,1000.*1E-4,1000.*1E-4,varyd0xd0[2][vcd0xd0],varyCosPoint[2][vccospoint]},/* 1<pt<2 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[3][vcd0xd0],varyCosPoint[3][vccospoint]},/* 2<pt<3 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[4][vcd0xd0],varyCosPoint[4][vccospoint]},/* 3<pt<4 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[5][vcd0xd0],varyCosPoint[5][vccospoint]},/* 4<pt<5*/     
							{0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[6][vcd0xd0],varyCosPoint[6][vccospoint]},/* 5<pt<6 */
							{0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[7][vcd0xd0],varyCosPoint[7][vccospoint]},/* 6<pt<8 */
							{0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[8][vcd0xd0],varyCosPoint[8][vccospoint]},/* 8<pt<12 */
							{0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[9][vcd0xd0],varyCosPoint[9][vccospoint]},/*12< pt <16*/
							{0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[10][vcd0xd0],varyCosPoint[10][vccospoint]}, /*16< pt <20*/
							{0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[11][vcd0xd0],varyCosPoint[11][vccospoint]}, /*20< pt <24*/
							{0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[12][vcd0xd0],varyCosPoint[12][vccospoint]}
	};/* pt > 24*/
	
	Float_t cutsMatrixD0toKpiLoose[nptbins][nvars]={{0.400,300.*1E-4,0.8,0.3,0.3,1000.*1E-4,1000.*1E-4,-35000.*1E-8,0.73},/* pt<0.5*/
							{0.400,300.*1E-4,0.8,0.3,0.3,1000.*1E-4,1000.*1E-4,-35000.*1E-8,0.73},/* 0.5<pt<1*/
							{0.400,300.*1E-4,0.8,0.4,0.4,1000.*1E-4,1000.*1E-4,-25000.*1E-8,0.75},/* 1<pt<2 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-15000.*1E-8,0.8},/* 2<pt<3 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.85},/* 3<pt<4 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.85},/* 4<pt<5 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.85},/* 5<pt<6 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.85},/* 6<pt<8 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-5000.*1E-8,0.85},/* 8<pt<12 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-0.*1E-8,0.85},/* 12<pt<16 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,0.*1E-8,0.85},/* 16<pt<20 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,0.*1E-8,0.85},/* 20<pt<24 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,0.*1E-8,0.85}};/* pt>24 */
	
	
	//CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
	Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
	for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
	Float_t **cutsMatrixTransposeLoose=new Float_t*[nvars];
	for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeLoose[iv]=new Float_t[nptbins];

	for (Int_t ibin=0;ibin<nptbins;ibin++){
	  for (Int_t ivar = 0; ivar<nvars; ivar++){
	    cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];
	    cutsMatrixTransposeLoose[ivar][ibin]=cutsMatrixD0toKpiLoose[ibin][ivar];
	    //printf("cutsMatrixD0toKpi[%d][%d] = %f\n",ibin, ivar,cutsMatrixD0toKpiStand[ibin][ivar]);
	  }
	}



	fCutsTight->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);
	fCutsLoose->SetCuts(nvars,nptbins,cutsMatrixTransposeLoose);


	for (Int_t ivar = 0; ivar<nvars; ivar++){
	  delete [] cutsMatrixTransposeStand[ivar];
	  delete [] cutsMatrixTransposeLoose[ivar];
	}
	delete [] cutsMatrixTransposeStand;
	cutsMatrixTransposeStand=NULL;
	delete [] cutsMatrixTransposeLoose;
	cutsMatrixTransposeLoose=NULL;



	fCutsTight->SetUseSpecialCuts(kTRUE);
	fCutsLoose->SetUseSpecialCuts(kTRUE);
	fCutsTight->SetRemoveDaughtersFromPrim(kTRUE);
	fCutsLoose->SetRemoveDaughtersFromPrim(kTRUE);
	// PID SETTINGS
	AliAODPidHF* pidObj=new AliAODPidHF();
	//pidObj->SetName("pid4D0");
	Int_t mode=1;
	const Int_t nlims=2;
	Double_t plims[nlims]={0.6,0.8}; //TPC limits in momentum [GeV/c]
	Bool_t compat=kTRUE; //effective only for this mode
	Bool_t asym=kTRUE;
	Double_t sigmas[5]={2.,1.,0.,3.,0.}; //to be checked and to be modified with new implementation of setters by Rossella
	pidObj->SetAsym(asym);// if you want to use the asymmetric bands in TPC
	pidObj->SetMatch(mode);
	pidObj->SetPLimit(plims,nlims);
	pidObj->SetSigma(sigmas);
	pidObj->SetCompat(compat);
	pidObj->SetTPC(kTRUE);
	pidObj->SetTOF(kTRUE);

	fCutsTight->SetPidHF(pidObj);
	fCutsLoose->SetPidHF(pidObj);
	delete pidObj; pidObj=NULL;
	fCutsTight->SetUsePID(kTRUE);
	fCutsLoose->SetUsePID(kTRUE);

	fCutsTight->SetUseDefaultPID(kFALSE);
	fCutsLoose->SetUseDefaultPID(kFALSE);

	// PILE UP REJECTION
	fCutsTight->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);
	fCutsLoose->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

	ptbinlimits=ptbins;
	fCutsTight->PrintAll();


	return nptbins;
 
}


//_________________________________________
Int_t AliAnalysisTaskSECharmFraction::SetStandardCuts(Double_t pt,Double_t invMassCut){
  // UPV: this should set the cut object

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

 

  /*//#######################################################################
  //###########################################################################
  //                    STANDARD SETS OF CUTS ("tight"~PPR like;  commented loose are more stringent than "tight")
  // #########################################################################
  Int_t ptbin=-1;
  if(pt>0. && pt<=1.) {
  ptbin=0;
    fCutsTight->SetD0toKpiCuts(invMassCut,0.04,0.8,0.5,0.5,0.05,0.05,-0.0002,0.5);
    // fCutsLoose->SetD0toKpiCuts(invMassCut,0.04,0.8,0.5,0.5,0.05,0.05,-0.00025,0.7);
    fCutsLoose->SetD0toKpiCuts(invMassCut,0.04,0.8,0.3,0.3,1.,1.,-0.0002,0.7);
  }   
  if(pt>1. && pt<=2.) {
    ptbin=1;  
    fCutsTight->SetD0toKpiCuts(invMassCut,0.03,0.8,0.6,0.6,0.05,0.05,-0.0002,0.6);
    //fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,1,1,-0.00025,0.8);
    fCutsLoose->SetD0toKpiCuts(invMassCut,0.03,0.8,0.3,0.3,1.,1.,-0.0001,0.7);
    //printf("I'm in the bin %d\n",ptbin);
  }
  if(pt>2. && pt<=3.) {
    ptbin=2;  
    fCutsTight->SetD0toKpiCuts(invMassCut,0.03,0.8,0.6,0.6,0.05,0.05,-0.0002,0.6);
    //fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,1,1,-0.00025,0.8);
    fCutsLoose->SetD0toKpiCuts(invMassCut,0.03,0.8,0.3,0.3,1.,1.,-0.0001,0.7);
    //printf("I'm in the bin %d\n",ptbin);
  } 
  if(pt>3. && pt<=5.){
    ptbin=3;  
    fCutsTight->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.0001,0.8);
    //    fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.8);
    fCutsLoose->SetD0toKpiCuts(invMassCut,0.03,0.8,0.4,0.4,1.,1.,-0.0001,0.75);
    //printf("I'm in the bin %d\n",ptbin);
  }
  if(pt>5.){
    ptbin=4;
    fCutsTight->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00005,0.8);
    //    fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.9);
    fCutsLoose->SetD0toKpiCuts(invMassCut,0.03,0.8,0.5,0.5,1.,1.,-0.00005,0.75);
  }//if(pt>5)
  return ptbin;
  //############################################################################
  */



  /* //#######################################################################
     //################# VARY CUTS for d0xd0 STUDY  ##########################

if(pt>0. && pt<=1.) {
     ptbin=0;
     fCutsTight->SetD0toKpiCuts(invMassCut,0.04,0.8,0.5,0.5,0.05,0.05,-0.0002,0.5);
     // fCutsLoose->SetD0toKpiCuts(invMassCut,0.04,0.8,0.5,0.5,0.05,0.05,-0.00025,0.7);
     fCutsLoose->SetD0toKpiCuts(invMassCut,0.04,0.8,0.3,0.3,1.,1.,-0.0002,0.7);
     }  
     if(pt>1. && pt<=2.) {
     ptbin=1;  
     fCutsTight->SetD0toKpiCuts(invMassCut,0.03,0.8,0.6,0.6,0.05,0.05,0.2,0.6);
     //fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,1,1,-0.00025,0.8);
     fCutsLoose->SetD0toKpiCuts(invMassCut,0.03,0.8,0.3,0.3,1.,1.,-0.0001,0.1);
     //printf("I'm in the bin %d\n",ptbin);
     }
     if(pt>2. && pt<=3.) {
     ptbin=2;  
     fCutsTight->SetD0toKpiCuts(invMassCut,0.03,0.8,0.6,0.6,0.05,0.05,0.2,0.6);
     //fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,1,1,-0.00025,0.8);
     fCutsLoose->SetD0toKpiCuts(invMassCut,0.03,0.8,0.3,0.3,1.,1.,-0.0001,0.1);
     //printf("I'm in the bin %d\n",ptbin);
     }  
     if(pt>3. && pt<=5.){
     ptbin=3;  
     fCutsTight->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,0.2,0.8);
     //    fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.8);
     fCutsLoose->SetD0toKpiCuts(invMassCut,0.03,0.3,0.4,0.4,1.,1.,-0.0001,0.1);
     //printf("I'm in the bin %d\n",ptbin);
     }
     if(pt>5.){
     ptbin=4;
     fCutsTight->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,0.2,0.8);
     //    fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.9);
     fCutsLoose->SetD0toKpiCuts(invMassCut,0.03,0.8,0.5,0.5,1.,1.,-0.00005,0.1);
     }//if(pt>5)
     return ptbin;
     //     #################################################################
  */    

  //##########################################################################
  //################## CUTS with d0xd0 cut released  #########################
  //###                    and TGHC cuts d0K and d0Pi to 0.1 instead of 0.05
  //### USED FOR PHDthesis
  //##########################################################################
 
  /* if(pt>0. && pt<=1.) {
     ptbin=0;
     fCutsTight->SetD0toKpiCuts(invMassCut,0.04,0.8,0.5,0.5,0.1,0.1,-0.000,0.5);
     // fCutsLoose->SetD0toKpiCuts(invMassCut,0.04,0.8,0.5,0.5,0.05,0.05,-0.00025,0.7);
     fCutsLoose->SetD0toKpiCuts(invMassCut,0.04,0.8,0.3,0.3,1.,1.,-0.000,0.7);
     }   
     if(pt>1. && pt<=2.) {
     ptbin=1;  
     fCutsTight->SetD0toKpiCuts(invMassCut,0.03,0.8,0.6,0.6,0.1,0.1,-0.000,0.6);
     //fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,1,1,-0.00025,0.8);
     fCutsLoose->SetD0toKpiCuts(invMassCut,0.03,0.8,0.4,0.4,1.,1.,-0.0000,0.7);
     //printf("I'm in the bin %d\n",ptbin);
     }
     if(pt>2. && pt<=3.) {
     ptbin=2;  
     fCutsTight->SetD0toKpiCuts(invMassCut,0.03,0.8,0.6,0.6,0.1,0.1,-0.000,0.6);
     //fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,1,1,-0.00025,0.8);
     fCutsLoose->SetD0toKpiCuts(invMassCut,0.03,0.8,0.4,0.4,1.,1.,-0.000,0.7);
     //printf("I'm in the bin %d\n",ptbin);
     } 
     if(pt>3. && pt<=5.){
     ptbin=3;  
     fCutsTight->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.1,0.1,-0.000,0.8);
     //    fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.8);
     fCutsLoose->SetD0toKpiCuts(invMassCut,0.03,0.8,0.5,0.5,1.,1.,-0.000,0.75);
     //printf("I'm in the bin %d\n",ptbin);
     }
     if(pt>5.){
     ptbin=4;
     fCutsTight->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.1,0.1,-0.0000,0.8);
     //    fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.9);
     fCutsLoose->SetD0toKpiCuts(invMassCut,0.03,0.8,0.5,0.5,1.,1.,-0.0000,0.75);
     }//if(pt>5)
     return ptbin;
  */




  //########## LOOKING FOR SIGNAL #####################
  /*  
  if(pt>0. && pt<=1.) {
    ptbin=0;
    fCutsTight->SetD0toKpiCuts(5*invMassCut,0.03,0.8,0.3,0.3,0.1,0.1,-0.00035,0.7);
    // fCutsLoose->SetD0toKpiCuts(invMassCut,0.04,0.8,0.5,0.5,0.05,0.05,-0.00025,0.7);
    fCutsLoose->SetD0toKpiCuts(5*invMassCut,0.04,0.8,0.3,0.3,0.1,0.1,-0.00025,0.7);
  }   
  if(pt>1. && pt<=2.) {
    ptbin=1;  
    fCutsTight->SetD0toKpiCuts(5*invMassCut,0.02,0.8,0.4,0.4,0.1,0.1,-0.00035,0.8);
    //fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,1,1,-0.00025,0.8);
    fCutsLoose->SetD0toKpiCuts(5*invMassCut,0.03,0.8,0.3,0.3,0.1,0.1,-0.0025,0.75);
    //printf("I'm in the bin %d\n",ptbin);
  }
  if(pt>2. && pt<=3.) {
    ptbin=2;  
    fCutsTight->SetD0toKpiCuts(5*invMassCut,0.02,0.8,0.7,0.7,0.1,0.1,-0.00026,0.94);
    //fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,1,1,-0.00025,0.8);
    fCutsLoose->SetD0toKpiCuts(5*invMassCut,0.02,0.8,0.7,0.7,0.1,0.1,-0.0002,0.92);
    //printf("I'm in the bin %d\n",ptbin);
  } 
  if(pt>3. && pt<=5.){
    ptbin=3;  
    fCutsTight->SetD0toKpiCuts(5*invMassCut,0.02,0.8,0.7,0.7,0.1,0.1,-0.00015,0.88);
    //    fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.8);
    fCutsLoose->SetD0toKpiCuts(5*invMassCut,0.02,0.8,0.7,0.7,0.1,0.1,-0.00015,0.9);
    //printf("I'm in the bin %d\n",ptbin);
  }
  if(pt>5.&& pt<=8.){
    ptbin=4;
    fCutsTight->SetD0toKpiCuts(5*invMassCut,0.015,0.8,0.7,0.7,0.1,0.1,-0.0001,0.9);
    //    fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.9);
    fCutsLoose->SetD0toKpiCuts(5*invMassCut,0.02,0.8,0.7,0.7,0.1,0.1,-0.0000,0.88);
  }//if(pt>5)
  if(pt>8.&&pt<=12.){
    ptbin=5;
    fCutsTight->SetD0toKpiCuts(5*invMassCut,0.015,0.8,0.7,0.7,0.1,0.1,-0.0001,0.9);
    //    fCutsLoose->SetD0toKpiCuts(invMassCut,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.9);
    fCutsLoose->SetD0toKpiCuts(5*invMassCut,0.015,0.8,0.7,0.7,0.1,0.1,-0.0005,0.88);
  }//if(pt>5)
  
  return ptbin;
  */
  printf("AliAnalysisTaskSECharmFraction::Obsolete method! Parameters pt=%f,invmasscut=%f not used \n",pt,invMassCut);
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
	

//__________________________________________________________________
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
  AliAODRecoDecayHF *aodDMC=0x0;
  Int_t pdgdaughters[2]={211,321};
  Int_t labMum=d->MatchToMC(421,arrayMC,2,pdgdaughters);
  if(labMum==-1){
    signaltype=-1;
    return aodDMC;    
  }

  // get daughter AOD tracks
  AliAODTrack *trk0 = (AliAODTrack*)d->GetDaughter(0);
  AliAODTrack *trk1 = (AliAODTrack*)d->GetDaughter(1);
 
  if(trk0==0x0||trk1==0x0){
    AliDebug(2,"Delete tracks I AM \n");
  
    signaltype=-1;
    return aodDMC;
   
  }

  b1=(AliAODMCParticle*)arrayMC->At(trk0->GetLabel());
  b2=(AliAODMCParticle*)arrayMC->At(trk1->GetLabel());
  mum1=(AliAODMCParticle*)arrayMC->At(labMum);  
  massMumTrue=mum1->GetCalcMass();
  aodDMC=ConstructFakeTrueSecVtx(b1,b2,mum1,primaryVtx);
 
  if(aodDMC==0x0){
    signaltype=10;
    return aodDMC;
  }

  Bool_t isfromDstar=kFALSE;
  grandmoth1=(AliAODMCParticle*)arrayMC->At(mum1->GetMother());
  if(TMath::Abs(grandmoth1->GetPdgCode())==413||TMath::Abs(grandmoth1->GetPdgCode())==423)isfromDstar=kTRUE;// D0 COMING FROM A D0*
  
  Int_t origin=CheckOrigin(arrayMC,mum1);
  if(origin==4){
      if(isfromDstar)signaltype=2;
      else signaltype=1;
      return aodDMC;
  }
  else if(origin==5){
    if(isfromDstar)signaltype=4;
    else signaltype=3;
    return aodDMC;
  }
  else if(origin==-1){
    signaltype=11;
    return aodDMC;
  }
  else if(origin==-2){
    signaltype=-1;
    return aodDMC;
  }
  
  signaltype=11;// JUST FOR SAFETY: SHOULD NEVER REACH THIS POINT
  return aodDMC;

}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSECharmFraction::CheckOrigin(const TClonesArray* arrayMC, const AliAODMCParticle *mcPartCandidate) const {		
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();

  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;
  while (mother >0 ){
    
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcGranma->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      return -2;
    }
  }
  
  if(!isQuarkFound)return -1;
  if(isFromB) return 5;  
  else return 4;
}

//__________________________________________________
AliAODRecoDecayHF* AliAnalysisTaskSECharmFraction::GetD0toKPiSignalTypeObsolete(const AliAODRecoDecayHF2Prong *d,TClonesArray *arrayMC,Int_t &signaltype,Double_t &massMumTrue,Double_t *primaryVtx){// OBSOLETE METHOD!!!!!
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
  if(dglabels[0]<0||dglabels[1]<0){
    AliDebug(2,"HERE I AM \n");

    //fake tracks
    
    signaltype=-1;
    return aodDMC;

  }
  //      printf("Before entering the MC checks \n");
  
  b1=(AliAODMCParticle*)arrayMC->At(trk0->GetLabel());
  b2=(AliAODMCParticle*)arrayMC->At(trk1->GetLabel());
  if(!b1||!b2){
    //Tracks with no mother  ??? FAKE DECAY VERTEX
    signaltype=10;
    return aodDMC;
  }
  if(b1->GetMother()<0||b2->GetMother()<0){
    //Tracks with no mother  ??? FAKE DECAY VERTEX
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
  
  if(mum1->GetMother()<0){
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
    if(grandmoth1->GetMother()<0){
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
      if(!fLikeSign)return 0x0;
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
Bool_t AliAnalysisTaskSECharmFraction::FillHistos(AliAODRecoDecayHF2Prong *d,TList *&list,Int_t ptbin,Int_t okD0,Int_t okD0bar,Double_t invMassD0,Double_t invMassD0bar,Bool_t isPeakD0,Bool_t isPeakD0bar,Bool_t isSideBandD0,Bool_t isSideBandD0bar,Double_t massmumtrue,AliAODRecoDecayHF *aodDMC,Double_t *vtxTrue){//FILL THE HISTOGRAMS: TAKE THE HISTOS FROM THE list NAME

  
  if((!okD0)&&(!okD0bar))return kTRUE;
  if(ptbin==-1)return kTRUE;
  //  flistNoCutsSignal->Add(hptD0NCsign);
  // flistNoCutsSignal->Add(hptD0VsMaxPtNCsign);
  // flistNoCutsSignal->Add(hptD0PTallsqrtNCsign);
  //  flistNoCutsSignal->Add(hptD0PTallNCsign);
  
  // %%%%%% TO BE DONE 
  //    flistNoCutsSignal->Add(hptD0vsptBNCsign);
  // flistNoCutsSignal->Add(hpD0vspBNCsign);
  //flistNoCutsSignal->Add(hptD0vsptcquarkNCsign);
  //flistNoCutsSignal->Add(hpD0vspcquarkNCsign);
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  // DONE
  //hd0zD0ptLSCsign_pt
  //hInvMassD0LSCsign_pt
  //hetaLSCsign_pt
  //
  // %%% TO BE DONE %% 
  //hCosPDPBLSCsign_pt
  //hCosPcPDLSCsign_pt
  
  Double_t pt=d->Pt();
  Double_t impparxy=d->ImpParXY()*10000.;



  // ######### Get Standard label for hist in tlist ###############
  TString namehist=list->GetName(),str;
  namehist.ReplaceAll("list","");

  //  ######### Global properties histos #################
  // ####### take care: only for candidates which pass the cuts !! not for side band ########
  if(fFastAnalysis<2){ 
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
      
      
      str="hd0singlTrack";
      str.Append(namehist.Data());
      ((TH1F*)list->FindObject(str.Data()))->Fill(d->Getd0Prong(0)*10000.);
      ((TH1F*)list->FindObject(str.Data()))->Fill(d->Getd0Prong(1)*10000.);
      
      str="hCPta";
      str.Append(namehist.Data());
      ((TH1F*)list->FindObject(str.Data()))->Fill(d->CosPointingAngle());
      
      str="hd0xd0";
      str.Append(namehist.Data());
      ((TH1F*)list->FindObject(str.Data()))->Fill(1e8*d->Prodd0d0());
      
      //%%%%%%%% NEW HISTO %%%%%%%%%%
      str="hdca";
      str.Append(namehist.Data());
      ((TH1F*)list->FindObject(str.Data()))->Fill(1e4*d->GetDCA());
      
      str="hcosthetastar";
      str.Append(namehist.Data());
      if(okD0)((TH1F*)list->FindObject(str.Data()))->Fill(d->CosThetaStarD0());
      if(okD0bar)((TH1F*)list->FindObject(str.Data()))->Fill(d->CosThetaStarD0bar());
      
      str="hptD0";
      str.Append(namehist.Data());
      ((TH1F*)list->FindObject(str.Data()))->Fill(pt);
      
      str="hptD0VsMaxPt";
      str.Append(namehist.Data());
      Int_t pr=0;
      if(d->PtProng(1)>d->PtProng(0))pr=1;
      if(d->PtProng(pr)<fptMax[0]) ((TH1F*)list->FindObject(str.Data()))->Fill(pt-fptMax[0]);
      else if(d->PtProng(TMath::Abs(pr-1))<fptMax[1])((TH1F*)list->FindObject(str.Data()))->Fill(pt-fptMax[1]);
      else ((TH1F*)list->FindObject(str.Data()))->Fill(pt-fptMax[2]);
      
      
      str="hptD0PTallsqrt";
      str.Append(namehist.Data());
      Double_t sumsqrpt=fptAllSq-d->PtProng(1)*d->PtProng(1)-d->PtProng(0)*d->PtProng(0);
      if(sumsqrpt>0.)((TH1F*)list->FindObject(str.Data()))->Fill(pt,TMath::Sqrt(sumsqrpt));
      
      str="hptD0PTall";
      str.Append(namehist.Data());
      ((TH1F*)list->FindObject(str.Data()))->Fill(pt,fptAll-d->PtProng(1)-d->PtProng(0));
      
      
      str="hd0zD0pt";
      str.Append(namehist.Data());
      str.Append("_pt");
      str+=ptbin;
      if(d->GetPrimaryVtx()!=0x0)((TH1F*)list->FindObject(str.Data()))->Fill(1e4*(d->Zv()-d->GetPrimaryVtx()->GetZ()));
      
      str="heta";
      str.Append(namehist.Data());
      str.Append("_pt");
      str+=ptbin;
      ((TH1F*)list->FindObject(str.Data()))->Fill(d->Eta());
      
      // OTHER NEW ADDITIONAL HISTOS
      
      str="hd0xd0";
      str.Append(namehist.Data());
      str.Append("_pt");
      str+=ptbin;
      //printf("Hist name: %s \n",str.Data());
      ((TH1F*)list->FindObject(str.Data()))->Fill(1e8*d->Prodd0d0());
      
      
      str="hd0D0VSd0xd0";
      str.Append(namehist.Data());
      str.Append("_pt");
      str+=ptbin;
      //printf("Hist name: %s \n",str.Data());
      ((TH2F*)list->FindObject(str.Data()))->Fill(1e8*d->Prodd0d0(),impparxy);
      
      
      str="hangletracksVSd0xd0";
      str.Append(namehist.Data());
      str.Append("_pt");
      str+=ptbin;
      //printf("Hist name: %s \n",str.Data());
      ((TH2F*)list->FindObject(str.Data()))->Fill(1e8*d->Prodd0d0(),d->ProngsRelAngle(0,1));
      
      str="hangletracksVSd0D0";
      str.Append(namehist.Data());
      str.Append("_pt");
      str+=ptbin;
      //  printf("Hist name: %s \n",str.Data());
      ((TH2F*)list->FindObject(str.Data()))->Fill(impparxy,d->ProngsRelAngle(0,1));
    // ####################################################
    }
  }  
  
  //  ######### Invariant mass histos #################
  if(fFastAnalysis<1){ 
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
  }

  // The Following is a NEW HISTO  
  str="hInvMassD0";
  str.Append(namehist.Data());
  str.Append("_pt");
  str+=ptbin;
  if(okD0)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0);
  if((!fsplitMassD0D0bar)&&okD0bar)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0bar);
  str="hInvMassD0bar";
  str.Append(namehist.Data());
  str.Append("_pt");
  str+=ptbin;
  if(fsplitMassD0D0bar&&okD0bar)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0bar);
  

  // FILLING OF THE SPARSE HISTO
  if(fFastAnalysis<=2){ // ONLY IF NOT VERY FAST ANALYSIS
    str="hSparse";
    str.Append(namehist.Data());

    Double_t point[5]={invMassD0,invMassD0bar,pt,impparxy,0.};
    if(okD0&&okD0bar)point[4]=3.5;
    else if(okD0)point[4]=1.5;
    else if(okD0bar)point[4]=2.5;
    if(fReadMC&&aodDMC!=0x0&&namehist.Contains("fromB")){     
      point[3]=aodDMC->ImpParXY()*10000.;
    }
    ((THnSparseF*)list->FindObject(str.Data()))->Fill(point);
    if(fReadMC&&aodDMC!=0x0&&namehist.Contains("fromB")){     
      point[3]=impparxy;
      str="hSparseReco";
      str.Append(namehist.Data());
      ((THnSparseF*)list->FindObject(str.Data()))->Fill(point);
    }

    
    str="hInvMassPt";
    str.Append(namehist.Data());
    if(okD0)((TH2F*)list->FindObject(str.Data()))->Fill(invMassD0,pt);
    if(okD0bar)((TH2F*)list->FindObject(str.Data()))->Fill(invMassD0bar,pt);

  }




  if(fFastAnalysis<=3&&namehist.Contains("sign")){
    str="hSparseCxyLxy";
    str.Append(namehist.Data()); 
    Double_t nLxy=d->NormalizedDecayLengthXY()*d->P()/pt;
    Double_t cosPxy=TMath::Abs(d->CosPointingAngleXY());
    Double_t point[4]={invMassD0,pt,cosPxy,nLxy};
    if(okD0){
      //      printf("Listname: %s, Here the histo : %p \n",namehist.Data(),((THnSparseF*)list->FindObject(str.Data())));
      ((THnSparseF*)list->FindObject(str.Data()))->Fill(point);
    }
    point[0]=invMassD0bar;
    if(okD0bar){
      ((THnSparseF*)list->FindObject(str.Data()))->Fill(point);
    }
  }



  /* if(isPeakD0||isPeakD0bar){
    str="hMass";
    str.Append(namehist.Data());
    str.Append("PM");
    if(isPeakD0&&okD0)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0);
    if(isPeakD0bar&&okD0bar)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0bar);
    // The Following is a NEW HISTO
    str="hInvMassD0";
    str.Append(namehist.Data());
    str.Append("_pt");
    str+=ptbin;
    if(isPeakD0&&okD0)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0);
    if(isPeakD0bar&&okD0bar)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0bar);
    }*/
  if(fFastAnalysis<2){ 
    if(isSideBandD0||isSideBandD0bar){
      str="hMass";
      str.Append(namehist.Data());
      str.Append("SB");
      if(okD0)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0);
      if(okD0bar)((TH1F*)list->FindObject(str.Data()))->Fill(invMassD0bar);
    }
  }
  
  if(fReadMC){
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
      if(isSideBandD0||isSideBandD0bar){
	str="hMassTrue";
	str.Append(namehist.Data());
	str.Append("SB");
	((TH1F*)list->FindObject(str.Data()))->Fill(massmumtrue);
      }
    }
  }

  // ################ D0 Impact Parameter Histos #####################
  if(isPeakD0||isPeakD0bar){    
   
    str="hd0D0";
    str.Append(namehist.Data());
    str.Append("PM");
    if(!fReadMC){
      // WE COUNT TWICE A CANDIDATE UNDER THE INV MASS PEAK BOTH AS a D0 and a D0bar (if selected) for DATA ONLY
      // THIS BECAUSE WE SUBTRACT a "BACKGROUND" AMOUNT ESTIMATED USING THE INV MASS FIT INFORMATION
      // WHICH SHOULD ACCOUNT FOR REFLECTIONS
      if(isPeakD0&&okD0){
	((TH1F*)list->FindObject(str.Data()))->Fill(impparxy);
      }
      if(isPeakD0bar&&okD0bar){
	((TH1F*)list->FindObject(str.Data()))->Fill(impparxy);
      }
    }
    else {
      if((isPeakD0&&okD0)||(isPeakD0bar&&okD0bar))((TH1F*)list->FindObject(str.Data()))->Fill(impparxy);
    }
    
    str="hd0D0pt";
    str.Append(namehist.Data());
    str.Append("_PkMss_pt");
    str+=ptbin;     
    if(!fReadMC){
      // WE COUNT TWICE A CANDIDATE UNDER THE INV MASS PEAK BOTH AS a D0 and a D0bar (if selected) for DATA ONLY
      // THIS BECAUSE WE SUBTRACT a "BACKGROUND" AMOUNT ESTIMATED USING THE INV MASS FIT INFORMATION
      // WHICH SHOULD ACCOUNT FOR REFLECTIONS
      if(isPeakD0&&okD0){
	((TH1F*)list->FindObject(str.Data()))->Fill(impparxy);
      }
      if(isPeakD0bar&&okD0bar){
	((TH1F*)list->FindObject(str.Data()))->Fill(impparxy);
      }
    }
    else {
      if((isPeakD0&&okD0)||(isPeakD0bar&&okD0bar))((TH1F*)list->FindObject(str.Data()))->Fill(impparxy);
    }
    
    
    if(fReadMC&&vtxTrue){
      // ONLY AN HISTO FOR QA: WE DO NOT CONSIDER THE IMPACT PARAMETER FOR EACH MASS HYPOTHESIS
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
    
    if(fReadMC&&aodDMC!=0x0){
      // WE NEED JUST THE SHAPE: AVOID TAKING THE IMPACT PAR FOR EACH MASS HYPO PASSING THE CUTS
      // aodDMC->Print("");
      //aodDMC->ImpParXY();
      // aodDMC->Print("");
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
  else if(isSideBandD0||isSideBandD0bar){
    // WE ASSUME THE IMPACT PARAMETER DISTRIBUION FOR BACKGROUND(SIDEBANDS) CANDIDATES
    // IS NOT CORRELATED TO THE INVARIANT MASSES. THEREFORE WE JUST TAKE ONE TIME
    // THE IMPACT PARAMETER AND NOT ONE FOR EACH MASS HYPOTHESIS PASSING THE CUTS

    str="hd0D0";
    str.Append(namehist.Data());
    str.Append("SB");
    ((TH1F*)list->FindObject(str.Data()))->Fill(impparxy);
    
    str="hd0D0pt";
    str.Append(namehist.Data());
    str.Append("_SBMss_pt");
    str+=ptbin;
    ((TH1F*)list->FindObject(str.Data()))->Fill(impparxy);
    
    
    if(fReadMC&&vtxTrue){
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
    
    if(fReadMC&&aodDMC!=0x0){
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


 void  AliAnalysisTaskSECharmFraction::FillHistoMCproperties(TClonesArray *arrayMC){ 
    //#############################################################
    //            HERE LOOK AT global properties of D0 mesons, c quarks and B
    // 
    //#############################################################
   Double_t pxyzMum[3],pxyzDaught[3],cosOpenAngle=-1.1,ptmum,ptdaught;
   Int_t ncdaught=0,cquarksMC=0,nD0all=0,nD0FromB=0,nBdaught=0,nD0bquark=0,nD0bMeson=0,nD0bBaryon=0;
   for (Int_t iPart=0; iPart<arrayMC->GetEntriesFast(); iPart++) { 
     AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iPart));
     if (!mcPart) {
       AliWarning("Particle not found in tree, skipping"); 
       continue;
     } 
     if (TMath::Abs(mcPart->GetPdgCode()) == 4){
       cquarksMC++;  
       mcPart->PxPyPz(pxyzMum);
       ptmum=mcPart->Pt();
       
       ((TH1F*)flistMCproperties->FindObject("hMCcquarkAllPt"))->Fill(ptmum);
       ((TH1F*)flistMCproperties->FindObject("hMCcquarkAllEta"))->Fill(mcPart->Eta());
       ((TH1F*)flistMCproperties->FindObject("hMCcquarkAllEnergy"))->Fill(mcPart->E());
       //NOW LOOK FOR A D0 among cquark daughters
       ncdaught=mcPart->GetDaughter(1)-mcPart->GetDaughter(0)+1;
       ((TH1F*)flistMCproperties->FindObject("hMCcquarkNdaught"))->Fill(ncdaught);
       if(ncdaught>1){
	 for(Int_t iDaught=mcPart->GetDaughter(0);iDaught<mcPart->GetDaughter(1);iDaught++){
	   AliAODMCParticle* mcPartD0 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iDaught));
	   if(mcPartD0==0x0)continue;
	   if(TMath::Abs(mcPartD0->GetPdgCode()) == 421){
	     // a D0 coming from a c quark
	     mcPartD0->PxPyPz(pxyzDaught);
	     ptdaught=mcPartD0->Pt();
	     ((TH1F*)flistMCproperties->FindObject("hMCD0fromcPt"))->Fill(ptdaught);
	     ((TH1F*)flistMCproperties->FindObject("hMCD0fromcEta"))->Fill(mcPartD0->Eta());
	     ((TH1F*)flistMCproperties->FindObject("hMCD0fromcEnergy"))->Fill(mcPartD0->E());
	     // ##############################################################################################
	     //                            Compare D0 momentum and c quarks: 
	     //              NB: here ALL D0 are considered, also those not decaying in KPi !!!
	     // ##############################################################################################
	     ((TH2F*)flistMCproperties->FindObject("hMCD0VscquarkPt"))->Fill(mcPart->Pt(),mcPartD0->Pt());
	     ((TH2F*)flistMCproperties->FindObject("hMCD0VscquarkEnergy"))->Fill(mcPart->E(),mcPartD0->E());
	     ((TH1F*)flistMCproperties->FindObject("hMCD0deltacquarkEnergy"))->Fill(mcPartD0->E()/mcPart->E());
	     ((TH1F*)flistMCproperties->FindObject("hMCD0EnergyVsAvcquarkDaughtEn"))->Fill((mcPartD0->E()-(mcPart->E()/ncdaught))/mcPart->E());
	     //calculate open angle
	     if((pxyzMum[0]!=0.||pxyzMum[1]!=0.||pxyzMum[2]!=0.)&&(pxyzDaught[0]!=0.||pxyzDaught[1]!=0.||pxyzDaught[2]!=0.))cosOpenAngle=(pxyzDaught[0]*pxyzMum[0]+pxyzDaught[1]*pxyzMum[1]+pxyzDaught[2]*pxyzMum[2])/(TMath::Sqrt(pxyzDaught[0]*pxyzDaught[0]+pxyzDaught[1]*pxyzDaught[1]+pxyzDaught[2]*pxyzDaught[2])*TMath::Sqrt(pxyzDaught[0]*pxyzDaught[0]+pxyzDaught[1]*pxyzDaught[1]+pxyzDaught[2]*pxyzDaught[2]));
	     ((TH1F*)flistMCproperties->FindObject("hMCD0cquarkAngle"))->Fill(cosOpenAngle);
	     ((TH2F*)flistMCproperties->FindObject("hMCD0cquarkAngleEnergy"))->Fill(mcPart->E(),cosOpenAngle);
	   }
	 }
       }
     }
     
     // NOW LOOK FOR D0 not coming from cquarks
      if (TMath::Abs(mcPart->GetPdgCode()) == 421){
	nD0all++;  
	if(mcPart->GetMother()<0)continue;
	AliAODMCParticle* mcD0Parent = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mcPart->GetMother()));
	if(mcD0Parent==0x0)continue;
	Bool_t notfound=kFALSE,bMeson=kFALSE,bBaryon=kFALSE;
	//CheckOrigin
	while(TMath::Abs(mcD0Parent->GetPdgCode())!=4&&TMath::Abs(mcD0Parent->GetPdgCode())!=5){
	  if(500<TMath::Abs(mcD0Parent->GetPdgCode())%10000&&TMath::Abs(mcD0Parent->GetPdgCode())<600){
	    bMeson=kTRUE;
	    break;
	  }
	  else if (5000<TMath::Abs(mcD0Parent->GetPdgCode())&&TMath::Abs(mcD0Parent->GetPdgCode())<6000){
	    bBaryon=kTRUE;
	    break;
	  }
	  if(mcD0Parent->GetMother()==0x0){
	    notfound=kTRUE;
	    break;
	  };
	  if(mcD0Parent->GetMother()<0){
	    notfound=kTRUE;
	    break;
	  }
	  mcD0Parent=dynamic_cast<AliAODMCParticle*>(arrayMC->At(mcD0Parent->GetMother()));
	  if(mcD0Parent==0x0) break;
	}

	if(mcD0Parent==0x0)continue;
	if(notfound)continue;
	if(TMath::Abs(mcD0Parent->GetPdgCode())==4)continue;//D0 from c quarks already counted
	((TH1F*)flistMCproperties->FindObject("hMCfromBpdgB"))->Fill(TMath::Abs(mcD0Parent->GetPdgCode()));
	if(bBaryon)nD0bBaryon++;
	else if(bMeson)nD0bMeson++;
	else nD0bquark++;
	nD0FromB++;
	mcD0Parent->PxPyPz(pxyzMum);
	ptmum=mcD0Parent->Pt();
	((TH1F*)flistMCproperties->FindObject("hMCBhadrPt"))->Fill(ptmum);
	((TH1F*)flistMCproperties->FindObject("hMCBhadrEta"))->Fill(mcD0Parent->Eta());
	((TH1F*)flistMCproperties->FindObject("hMCBhadrEnergy"))->Fill(mcD0Parent->E());
	
	nBdaught=mcD0Parent->GetDaughter(1)-mcD0Parent->GetDaughter(0)+1;
	((TH1F*)flistMCproperties->FindObject("hMCBhadrNdaught"))->Fill(nBdaught);

	
	// Now take properties of this D0 coming from a B
	mcPart->PxPyPz(pxyzDaught);
	ptdaught=mcPart->Pt();
	((TH1F*)flistMCproperties->FindObject("hMCD0fromBPt"))->Fill(ptdaught);
	((TH1F*)flistMCproperties->FindObject("hMCD0fromBEta"))->Fill(mcPart->Eta());
	((TH1F*)flistMCproperties->FindObject("hMCD0fromBEnergy"))->Fill(mcPart->E());
	// ##############################################################################################
	//                            Compare D0 momentum and b hadron: 
	//              NB: here ALL D0 are considered, also those not decaying in KPi !!!
	// ##############################################################################################
	((TH2F*)flistMCproperties->FindObject("hMCD0VsBhadrPt"))->Fill(mcD0Parent->Pt(),mcPart->Pt());
	((TH2F*)flistMCproperties->FindObject("hMCD0VsBhadrEnergy"))->Fill(mcD0Parent->E(),mcPart->E());
	((TH1F*)flistMCproperties->FindObject("hMCD0deltaBhadrEnergy"))->Fill(mcPart->E()/mcD0Parent->E());
	((TH1F*)flistMCproperties->FindObject("hMCD0EnergyVsAvBDaughtEn"))->Fill((mcPart->E()-(mcD0Parent->E()/nBdaught))/mcD0Parent->E());
	//calculate open angle
	if((pxyzMum[0]!=0.||pxyzMum[1]!=0.||pxyzMum[2]!=0.)&&(pxyzDaught[0]!=0.||pxyzDaught[1]!=0.||pxyzDaught[2]!=0.))cosOpenAngle=(pxyzDaught[0]*pxyzMum[0]+pxyzDaught[1]*pxyzMum[1]+pxyzDaught[2]*pxyzMum[2])/(TMath::Sqrt(pxyzDaught[0]*pxyzDaught[0]+pxyzDaught[1]*pxyzDaught[1]+pxyzDaught[2]*pxyzDaught[2])*TMath::Sqrt(pxyzDaught[0]*pxyzDaught[0]+pxyzDaught[1]*pxyzDaught[1]+pxyzDaught[2]*pxyzDaught[2]));
	((TH1F*)flistMCproperties->FindObject("hMCD0BhadrAngle"))->Fill(cosOpenAngle);
	((TH2F*)flistMCproperties->FindObject("hMCD0BhadrAngleEnergy"))->Fill(mcPart->E(),cosOpenAngle);
      }
   }
   ((TH1F*)flistMCproperties->FindObject("hMCPartFound"))->Fill(1,cquarksMC);
   ((TH1F*)flistMCproperties->FindObject("hMCPartFound"))->Fill(2,nD0all);
   ((TH1F*)flistMCproperties->FindObject("hMCPartFound"))->Fill(3,nD0FromB);
   ((TH1F*)flistMCproperties->FindObject("hMCPartFound"))->Fill(4,nD0bMeson);
   ((TH1F*)flistMCproperties->FindObject("hMCPartFound"))->Fill(5,nD0bBaryon);
   
 }


void AliAnalysisTaskSECharmFraction::SetPtBins(Int_t nbins,const Float_t *ptbins){
  if((fptbins)!=0x0)delete fptbins;
  fnbins=nbins;fptbins=new Float_t[fnbins];
  memcpy(fptbins,ptbins,fnbins*sizeof(Float_t));
  return;
}

void AliAnalysisTaskSECharmFraction::SetStandardMassSelection(){
  //SET THE DEFAULT VALUES FOR INVARIANT MASS SELECTION

  /*HERE DEFAULT
    SetSignalInvMassCut();
    SetLargeInvMassCut();
    SetSideBandInvMassCut();
    SetSideBandInvMassWindow();
  */

  // HERE FOR SEARCH FOR SIGNAL
  SetSignalInvMassCut();
  SetLargeInvMassCut();
  SetSideBandInvMassCut();
  SetSideBandInvMassWindow();
  return;
}

Bool_t AliAnalysisTaskSECharmFraction::SpecialSelD0(AliAODRecoDecayHF2Prong *d,Int_t &nusedforVtx){
  
  AliAODTrack *trk0 = (AliAODTrack*)d->GetDaughter(0);
  AliAODTrack *trk1 = (AliAODTrack*)d->GetDaughter(1);
  nusedforVtx=0;
  if(trk0->GetUsedForPrimVtxFit())nusedforVtx++;
  if(trk1->GetUsedForPrimVtxFit())nusedforVtx++;
  if(nusedforVtx>fNtrMaxforVtx)return kFALSE;
  
  if(TMath::Abs(d->Getd0Prong(1)) < -99999.  || 
     TMath::Abs(d->Getd0Prong(0)) < -99999.) return kFALSE;
  
  return kTRUE;
}



AliAODVertex* AliAnalysisTaskSECharmFraction::GetPrimaryVtxSkipped(AliAODEvent *aodev,AliAODRecoDecayHF2Prong *d){
  //Calculate the primary vertex w/o the daughter tracks of the candidate
  
  AliESDVertex *vertexESD=0x0;
  AliAODVertex *vertexAOD=0x0;
  AliVertexerTracks *vertexer = new AliVertexerTracks(aodev->GetMagneticField());
  
  Int_t skipped[2];
  Int_t nTrksToSkip=2;
  AliAODTrack *dgTrack = (AliAODTrack*)d->GetDaughter(0);
  skipped[0]=dgTrack->GetID();
  dgTrack = (AliAODTrack*)d->GetDaughter(1);
  skipped[1]=dgTrack->GetID();

 
  //
  vertexer->SetSkipTracks(nTrksToSkip,skipped);
  vertexESD = (AliESDVertex*)vertexer->FindPrimaryVertex(aodev); 
  vertexer->SetMinClusters(4);  
  if(!vertexESD) return vertexAOD;
  if(vertexESD->GetNContributors()<=0) { 
    AliDebug(2,"vertexing failed"); 
    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }
  
  delete vertexer; vertexer=NULL;
  
  
  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  delete vertexESD; vertexESD=NULL;
  
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF);
  return vertexAOD;
  
}



 Bool_t AliAnalysisTaskSECharmFraction::FillAziList(AliAODEvent *aod,Double_t azilist[30000],Int_t trkIDlist[30000],Int_t &nprim)const{
   Int_t ntracks=aod->GetNumberOfTracks();
   Double_t ptmin=1.;
   if(ntracks>30000){
     nprim=1;
     return kFALSE;	  
   }
   nprim=0;
   for(Int_t it=0;it<ntracks;it++) {
     AliAODTrack *track = aod->GetTrack(it);
     
     if(track->IsPrimaryCandidate()){
       if(track->Pt()>ptmin){
	 
	 azilist[nprim]=track->Phi();
	 trkIDlist[nprim]=track->GetID();
	 nprim++;
       }
     }
   }
   return kTRUE;
 }
 
 
 void AliAnalysisTaskSECharmFraction::FillAziHistos(AliAODRecoDecayHF2Prong *d,TList *&list,Int_t ptbin,Double_t azilist[30000],Int_t trkIDlist[30000],Int_t nprim,Int_t okD0,Int_t okD0bar,Bool_t isPeakD0,Bool_t isPeakD0bar,Bool_t isSideBandD0,Bool_t isSideBandD0bar)const{
   
   if((!okD0)&&(!okD0bar))return;
   if(ptbin==-1)return;
   TString namehist=list->GetName(),str;
   namehist.ReplaceAll("list","");
   //   Double_t ptD=d->Pt();
 
   str="hPhiHist";
   if(isPeakD0||isPeakD0bar)str.Append("PM");
   else if(isSideBandD0||isSideBandD0bar)str.Append("SB");
   else return;
   str.Append(namehist.Data());
   str.Append("_pt");
   str+=ptbin;
   
   AliAODTrack *dtr;
   dtr=(AliAODTrack*)d->GetDaughter(0);
   Int_t id1=dtr->GetID();
   dtr=(AliAODTrack*)d->GetDaughter(1);
   Int_t id2=dtr->GetID();
   
   Double_t phi=d->Phi();	
   Double_t weight=1./nprim;
   Double_t azi;
   for(Int_t j=0;j<nprim;j++){
     if(trkIDlist[j]!=id1&&trkIDlist[j]!=id2){
       azi=azilist[j]-phi;
       if(azi>TMath::Pi())azi-=2.*TMath::Pi();
       else if(azi<-TMath::Pi())azi+=2.*TMath::Pi();
       
       ((TH1F*)list->FindObject(str.Data()))->Fill(azi,weight);
     }
   }
   
   
 }









void AliAnalysisTaskSECharmFraction::Terminate(const Option_t*){
  //TERMINATE METHOD: NOTHING TO DO



}
