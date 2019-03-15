// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Hege Erdal       <hege.erdal@gmail.com>               *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliDxHFEParticleSelectionMCEl.cxx
/// @author Hege Erdal, Matthias Richter
/// @date   2012-07-19
/// @brief  MC El selection for D0-HFE correlation
///

#include "AliDxHFEParticleSelectionMCEl.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "THnSparse.h"
#include "AliAODMCParticle.h"
#include "TH1F.h"
#include "TAxis.h"
#include "AliAODTrack.h"
#include "AliReducedParticle.h"
#include "TList.h"
#include <iostream>
#include <cerrno>
#include <memory>

#include "AliGenHijingEventHeader.h"
#include "AliVEvent.h"
#include "AliAODMCHeader.h"
#include "TClonesArray.h"
using namespace std;
using std::vector;


/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEParticleSelectionMCEl)

AliDxHFEParticleSelectionMCEl::AliDxHFEParticleSelectionMCEl(const char* opt)
  : AliDxHFEParticleSelectionEl(opt)
  , fMCTools()
  , fHistoList(NULL)
  , fOriginMother(0)
  , fResultMC(0)
  , fUseKine(kFALSE)
  , fMotherPDGs(0)
  , fUseMCReco(kFALSE)
  , fSelectionStep(AliDxHFEParticleSelectionEl::kNoCuts)
  , fStoreCutStepInfo(kFALSE)
  , fElSelection(kAllPassingSelection)
  , fStoreOnlyMCElectrons(kFALSE)
  , fMCInfo(kMCLast)
  , fRemoveEfromD0(kFALSE)
  , fRemoveSecondary(kFALSE)
  , fUseGenerator(kFALSE)
  , fGenerator(-1)
  , fUseHIJINGOnly(kFALSE)
  , fSystem(0)
{
  // constructor
  // 
  // 
  // 
  // 

  fMCTools.~AliDxHFEToolsMC();

  ParseArguments(opt);

  // This is also checked in AddTasks, but just for safety! 
  if(fUseMCReco && fUseKine) AliFatal("CAN'T SET BOTH usekine AND elmcreco AT THE SAME TIME");

  bool kineTools=kFALSE;
  TString strOption(opt);
  cout << strOption << endl;
  if (strOption.Contains("OnlyKineTools")) kineTools=kTRUE;

  // TODO: argument scan, build tool options accordingly
  // e.g. set mc mode first/last, skip control histograms
  TString toolopt("pdg=11");
  if(fMCInfo==kMCLast) toolopt+=" mc-last";
  if(fMCInfo==kMCOnly) toolopt+=" mc-first";
  if(fMCInfo==kMCFirst) toolopt+=" mc-first";
  if(fRemoveSecondary )toolopt+=" removesecondary";
  else toolopt+=" keepsecondary";
  if(fUseKine || kineTools){cout << "sending usekine to tools " << endl; toolopt+=" usekine";}
  new (&fMCTools) AliDxHFEToolsMC(toolopt);
}

//at the moment not used, keep for now (to be used when plotting)
const char* AliDxHFEParticleSelectionMCEl::fgkPDGMotherBinLabels[]={
  "d",
  "u",
  "s",
  "c",
  "b",
  "gluon",
  "gamma",
  "#pi^{0}",
  "#eta",
  "proton",
  "others"
};

                             
const char*  AliDxHFEParticleSelectionMCEl::fgkPDGBinLabels[]={
  "positron",
  "electron",
  "#mu+",
  "#mu-",
  "#pi+",
  "#pi-",
  "K+",
  "K-",
  "proton",
  "antiproton",
  "others"
};

AliDxHFEParticleSelectionMCEl::~AliDxHFEParticleSelectionMCEl()
{
  // destructor
  if(fHistoList){
    delete fHistoList;
    fHistoList=NULL;
  }
}

int AliDxHFEParticleSelectionMCEl::Init()
{
  //
  // init function
  // 

  // Particles considered HFE background. can be expanded
  fMotherPDGs.push_back(AliDxHFEToolsMC::kPDGpi0); 
  fMotherPDGs.push_back(AliDxHFEToolsMC::kPDGeta);
  fMotherPDGs.push_back(AliDxHFEToolsMC::kPDGgamma);
  fMotherPDGs.push_back(AliDxHFEToolsMC::kPDGJpsi);

  int iResult=0;
  iResult=AliDxHFEParticleSelectionEl::Init();
  if (iResult<0) return iResult;

  // Histo containing PDG of track which was not MC truth electron
  fHistoList=new TList;
  fHistoList->SetName("HFE MC Histograms");
  fHistoList->SetOwner();

  fHistoList->Add(CreateControlHistogram("fPDGnotMCElectron","PDG of track not MC truth electron",AliDxHFEToolsMC::kNofPDGLabels,fgkPDGBinLabels));
  fHistoList->Add(CreateControlHistogram("fPDGNotHFMother","PDG of mother not HF",5000));  
  fHistoList->Add(CreateControlHistogram("fPDGHadronMother","PDG of mother of hadrons",5000));  

  AddControlObject(fHistoList);

  return 0;
}

THnSparse* AliDxHFEParticleSelectionMCEl::DefineTHnSparse()
{
  //
  // Defines the THnSparse. 

  const double Pi=TMath::Pi();
  TString name;
  THnSparse* thn=NULL;
  name.Form("%s info", GetName());

  int nrMotherEl=AliDxHFEToolsMC::kNrOrginMother+kNrBackground;

  if(fStoreCutStepInfo){
    const int thnSizeExt =6;

    InitTHnSparseArray(thnSizeExt);

    // TODO: Redo binning of distributions more?
    //     		              0    1      2     3              4             5
    // 	 	                      Pt   Phi   Eta   mother       generator    CutstepInfo
    int    thnBinsExt[thnSizeExt] = { 50,  100, 100,nrMotherEl+1,      4,        kNCutLabels-1 };
    double thnMinExt [thnSizeExt] = {   0,    0, -1.,  -1.5,        -1.5,    kRecKineITSTPC-0.5};
    double thnMaxExt [thnSizeExt] = {  10, 2*Pi,  1.,nrMotherEl-0.5, 2.5,       kSelected-0.5};
    const char* thnNamesExt[thnSizeExt]={
      "Pt",
      "Phi",
      "Eta", 
      "Mother", //bin==-1: Not MC truth electron
      "Generator",
      "Last survived cut step"
    };
    thn=(THnSparse*)CreateControlTHnSparse(name,thnSizeExt,thnBinsExt,thnMinExt,thnMaxExt,thnNamesExt);
  }
  else{
    if(!fUseGenerator){
      const int thnSize =4;
      InitTHnSparseArray(thnSize);
      
      // TODO: Redo binning of distributions more?
      //     		       0       1      2     3            
      // 	 	               Pt     Phi   Eta   mother
      int    thnBins[thnSize] = { 50,  100, 100, nrMotherEl+1 };
      double thnMin [thnSize] = {   0,    0, -1.,  -1.5 };
      double thnMax [thnSize] = {  10, 2*Pi,  1.,nrMotherEl-0.5 };
      const char* thnNames[thnSize]={
	"Pt",
	"Phi",
	"Eta", 
	"Mother" //bin==-1: Not MC truth electron

      };
      thn=(THnSparse*)CreateControlTHnSparse(name,thnSize,thnBins,thnMin,thnMax,thnNames);
    }else{
      const int thnSizeGen =5;
      InitTHnSparseArray(thnSizeGen);
      
      // TODO: Redo binning of distributions more?
      //     		       0       1      2     3           4 
      // 	 	               Pt     Phi   Eta   mother     generator 
      int    thnBins[thnSizeGen] = { 50,  100, 100, nrMotherEl+1,      4  };
      double thnMin [thnSizeGen] = {   0,    0, -1.,  -1.5    ,     -1.5 };
      double thnMax [thnSizeGen] = {  10, 2*Pi,  1.,nrMotherEl-0.5,  2.5};
      const char* thnNames[thnSizeGen]={
	"Pt",
	"Phi",
	"Eta", 
	"Mother", //bin==-1: Not MC truth electron
	"Generator"
      };
      thn=(THnSparse*)CreateControlTHnSparse(name,thnSizeGen,thnBins,thnMin,thnMax,thnNames);
    }
  }
  return thn;
}

int AliDxHFEParticleSelectionMCEl::FillParticleProperties(AliVParticle* p, Double_t* data, int dimension) const
{
  // fill the data array from the particle data
  if (!data) return -EINVAL;
  if (!p) return -ENODATA;
  // handle different types of tracks, can be extended
  AliReducedParticle *trRP=dynamic_cast<AliReducedParticle*>(p);
  int i=0;
  if (dimension!=GetDimTHnSparse()) {
    // TODO: think about filling only the available data and throwing a warning
    return -ENOSPC;
  }
  memset(data, 0, dimension*sizeof(data[0]));
  data[i++]=p->Pt();
  data[i++]=p->Phi();
  data[i++]=p->Eta();
  if (trRP) data[i++]=trRP->GetOriginMother();
  else data[i++]=fOriginMother;
  if(fUseGenerator)data[i++]=fGenerator;//i++ or only i here?
  i++; // take out of conditionals to be save
  if (i<dimension) {
    if(fStoreCutStepInfo) data[i]=GetLastSurvivedCutsStep();
    i++; // take out of conditionals to be save
  }
  i++; // take out of conditionals to be save
  return i;
}

int AliDxHFEParticleSelectionMCEl::IsSelected(AliVParticle* p, const AliVEvent* pEvent)
{
  /// overloaded from AliDxHFEParticleSelection: check particle
  /// H: Have changed function. Now doing particle selection first, then run MC over 
  /// selected tracks. Could configure it to be configurable, but not sure if it
  /// is needed.  
  /// Result from normal track selection is returned, result from MC is stored in
  /// THnSparse. 

  int iResult=0;
  fResultMC=0;
  if (!p || !pEvent){
    return -EINVAL;
  }
  fOriginMother=-1;

  if(!fUseKine && !fUseMCReco){
  // step 1:
  // optional MC selection before the particle selection
    if (fMCTools.MCFirst()){
      fResultMC=CheckMC(p, pEvent);
      // histograming?
      if(fMCInfo==kMCOnly) return fResultMC;
      if(fResultMC==0) return fResultMC;
    }

  // step 2 or 1, depending on sequence:
  // normal particle selection
  iResult=AliDxHFEParticleSelectionEl::IsSelected(p, pEvent);
  if (fMCTools.MCFirst() || iResult==0) return iResult;


  }

  // step 2, only executed if MC check is last
  // optional MC selection after the particle selection
  // result stored to be filled into THnSparse
  // TODO: strictly speaken the particles should be rejected
  // if not mc selected, however skip this for the moment, because of
  // the logic outside
  // This line will run always when running directly over the stack or over MC reconstructed tracks
  // For MC reconstructed tracks, should consider adding more constrictions on the tracks,
  // maybe do the track selection first (which means could call AliDxHFEParticleSelectionEl::IsSelected()
  // and don't do PID
  
  if(fUseMCReco) {// && ( fSelectionStep > AliDxHFEParticleSelectionEl::kNoCuts )){
    // always check the base class method in mode fUseMCReco
    iResult=AliDxHFEParticleSelectionEl::IsSelected(p, pEvent);
    if(iResult == 0) return iResult;
  }
  fResultMC=CheckMC(p, pEvent);
  
  // Return only result of MC checks (rather than result from reconstruction) for kinematical studies, 
  // when looking at the various selection steps (want only electrons) or when we only want to store specific
  // sources from the reconstruction (hadrons/nonHFE/HFE/onlyc/onlyb)
  if(fUseKine || fUseMCReco || fElSelection>kAllPassingSelection || fStoreOnlyMCElectrons)
    return fResultMC;

  return iResult;
   
}

int AliDxHFEParticleSelectionMCEl::CheckMC(AliVParticle* p, const AliVEvent* pEvent)
{
  /// check if MC criteria are fulfilled
  if (!p || !pEvent){
    return -EINVAL;
  }
  fOriginMother=-1;
  int iResult=0;

  if (!fMCTools.IsInitialized() && (iResult=fMCTools.InitMCParticles(pEvent))<0) {
    // TODO: message? but has to be filtered in order to avoid message flood
    return 0; // no meaningful filtering on mc possible
  }

  if(fUseKine){
    // If run on kinematical level, check if particle is electron (only ones who are relevant) 
    // and if pass IsPhysicalPrimary()
    bool test = fMCTools.CheckMCParticle(p);
    if(!test) return 0;
    // Add additional contraints on the electrons? 
    // TODO: remove delta e? also remove E<300MeV?
    // Should also mark dalitz decay and gamma conversion..
    AliAODMCParticle *mcp=dynamic_cast<AliAODMCParticle*>(p);
    if(!mcp ) {return 0;}
    //if(fRemoveSecondary)
    //  if(!mcp->IsPhysicalPrimary()){cout <<"not physprim" << endl; return 0; }
    if(TMath::Abs(mcp->Eta())>0.8){ return 0;}
  }

  int pdgFirstMother=0;  
  int pdgOriginMother=0;
  int pdgParticle=-1;
  if (!fUseKine && !fMCTools.CheckMCParticle(p, &pdgParticle)) {
    // rejected by pdg (and maybe IsPhysicalPrimary
    // TODO: Move this to fMCTools???? Can this be part of the statistics in the MC class?
    ((TH1D*)fHistoList->FindObject("fPDGnotMCElectron"))->Fill(fMCTools.MapPDGLabel(pdgParticle));
    ((TH1D*)fHistoList->FindObject("fPDGHadronMother"))->Fill(fMCTools.FindMotherPDG(p,AliDxHFEToolsMC::kGetFirstMother));

    // If the hadrons (and only hadrons) are going to be stored
    if(fElSelection==kHadron) return 1;
    else return 0;
  }

  // Find PDG of first mother
  pdgFirstMother=fMCTools.FindMotherPDG(p,AliDxHFEToolsMC::kGetFirstMother);

  if(fRemoveEfromD0 && ( TMath::Abs(pdgFirstMother)==AliDxHFEToolsMC::kPDGD0)){
    AliDebug(2,"Electron comes from D0");
    return 0;
  }

  // Check if first mother is counted as background
  Bool_t isNotBackground=fMCTools.RejectByPDG(pdgFirstMother,fMotherPDGs);
  Int_t selection=-1;

  // loops back to find original quark + if there was a gluon. 
  pdgOriginMother=fMCTools.FindMotherPDG(p,AliDxHFEToolsMC::kGetOriginMother);
  int motherCode=fMCTools.GetOriginMother();
  bool partHF=kFALSE;

  // Test if the origin mother is HF
  if(motherCode>=AliDxHFEToolsMC::kOriginCharm && motherCode<=AliDxHFEToolsMC::kOriginBeauty)
    partHF=kTRUE;

  if(motherCode>=AliDxHFEToolsMC::kOriginGluonCharm && motherCode<=AliDxHFEToolsMC::kOriginGluonBeauty)
    partHF=kTRUE;


  if(!isNotBackground){
    // Set fOriginMother if mother is counted as background
    // TODO: Could this be done in a more elegant way?
    switch(pdgFirstMother){
    case(AliDxHFEToolsMC::kPDGpi0): if(partHF) fOriginMother=AliDxHFEToolsMC::kNrOrginMother+kPi0HF; else fOriginMother=AliDxHFEToolsMC::kNrOrginMother+kPi0; break;
    case(AliDxHFEToolsMC::kPDGeta): if(partHF) fOriginMother=AliDxHFEToolsMC::kNrOrginMother+kEtaHF; else fOriginMother=AliDxHFEToolsMC::kNrOrginMother+kEta; break;
    case(AliDxHFEToolsMC::kPDGgamma): if(partHF) fOriginMother=AliDxHFEToolsMC::kNrOrginMother+kGammaHF; else fOriginMother=AliDxHFEToolsMC::kNrOrginMother+kGamma;break;
    case(AliDxHFEToolsMC::kPDGJpsi): fOriginMother=AliDxHFEToolsMC::kNrOrginMother+kJPsi; break;
    }
    selection=kNonHFE;
  }
  else{
    // If loop over Stack, also checks if First mother is a HF meson
    Bool_t isHFmeson =fMCTools.TestMotherHFMeson(TMath::Abs(pdgFirstMother));

    if(isHFmeson){
      // If first mother is a HF meson,
      //Result is  stored in fOriginMother
      //pdgMother=fMCTools.FindMotherPDG(p,AliDxHFEToolsMC::kGetOriginMother);
      fOriginMother=motherCode; //fMCTools.GetOriginMother();
    }
    else{
      //NotHFmother
      ((TH1D*)fHistoList->FindObject("fPDGNotHFMother"))->Fill(pdgFirstMother);
      if(partHF) fOriginMother=AliDxHFEToolsMC::kNrOrginMother+kOtherHF;
      else fOriginMother=AliDxHFEToolsMC::kNrOrginMother+kOther;
      selection=kNonHFE;
    }

  }


  // Checks on the electrons, to return only specific sources
  if(fOriginMother >= AliDxHFEToolsMC::kOriginNone && fOriginMother <= AliDxHFEToolsMC::kOriginStrange)
    selection=kNonHFE;

  if(fOriginMother == AliDxHFEToolsMC::kOriginCharm || fOriginMother == AliDxHFEToolsMC::kOriginGluonCharm)
    selection=kOnlyc;

  if(fOriginMother == AliDxHFEToolsMC::kOriginBeauty || fOriginMother == AliDxHFEToolsMC::kOriginGluonBeauty)
    selection=kOnlyb;

  if(fElSelection==kNonHFE && ( selection!= kNonHFE)){
    AliDebug(2,Form("Particle selected as: %d, want only to select %d",selection, fElSelection));
    return 0;
  }
  if(fElSelection==kHFE && (selection != kOnlyc && selection != kOnlyb)) {
    AliDebug(2,Form("Particle selected as: %d, want only to select %d",selection, fElSelection));
    return 0;
  }
  if(fElSelection==kOnlyc && selection !=kOnlyc){
    AliDebug(2,Form("Particle selected as: %d, want only to select %d",selection, fElSelection));
    return 0;
  }
  if(fElSelection==kOnlyb && selection !=kOnlyb){
    AliDebug(2,Form("Particle selected as: %d, want only to select %d",selection, fElSelection));
    return 0;
  }

  TString nameGen;	
  fGenerator=0;// adding per track: -1 hijing; 0-pythiaHF; 2-the rest  --> -2= pure hijing -> back ok (but check identity); 0= phytia,pythia-> ok for checking signal; reject the rest: -1 (Hij,pyt), 1 (hij, rest), 2 (pythia,rest) ,4 (rest,rest) 
      
  if(fUseGenerator){
    //    originvsGen=-1; 
    if(fSystem!=2){   
      AliAODTrack *aodtrack = static_cast<AliAODTrack *>(p);
      if(!aodtrack) {cout << "no AODtrack" << endl; return -1;}
      
      if(fUseKine){
	fMCTools.GetTrackPrimaryGenerator(aodtrack,nameGen,kTRUE);
      }
      else{
	fMCTools.GetTrackPrimaryGenerator(aodtrack,nameGen);
      }
      
      
      if(nameGen.Contains("ijing")){
	//cout << "hijing generator" << endl;
	fGenerator=1;
      }
      else if(nameGen.Contains("ythia")){
	//cout << "pythia generator" << endl;
	fGenerator=2;
      }
    }else{ //pPb case
      //Remove tracks that are not from the HIJING generator, e.g., in LHC13d3 which has an enhaced sample from PYTHIA. 
      //Selecting only HIJING lets the sample be used as a minimum bias sample
      //Be vary of how this is set for D0 as well
      
	AliAODMCHeader* mcHeader = dynamic_cast<AliAODMCHeader*>(pEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
	AliAODMCParticle *mcpart=dynamic_cast<AliAODMCParticle*>(p);
	AliAODTrack *aodtrack = static_cast<AliAODTrack *>(p);
	TClonesArray *mcArray = dynamic_cast<TClonesArray*>(pEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	
	if(fUseKine){ //Kine
	  Bool_t fromhijing = kFALSE;
	  Int_t Prim = GetPrimary((mcpart->GetLabel()), mcArray);
	  AliAODMCParticle *AODMCtrackPrim = (AliAODMCParticle*)mcArray->At(Prim);
	  Int_t trkIndexPrim = AODMCtrackPrim->GetLabel();//gives index of the particle in original MCparticle
	  if(IsFromBGEventAOD(trkIndexPrim, pEvent)) fromhijing = kTRUE;//check if the particle comes from hijing or from enhanced event
	  if(!fromhijing){
	    fGenerator=2;
	    if(fUseHIJINGOnly){return 0;}
	  }else{fGenerator=1;}
	}else{ //Reco
	  Bool_t MCElectron=kFALSE;
	  Int_t trkLabel = TMath::Abs(p->GetLabel());
	  AliAODMCParticle *MCtrack = (AliAODMCParticle*)mcArray->At(trkLabel);
	  Int_t PrimElectron = GetPrimary(trkLabel, mcArray);
	  AliAODMCParticle *AODMCElectron = (AliAODMCParticle*)mcArray->At(PrimElectron);
	  Int_t trkIndexElectron = AODMCElectron->GetLabel();//gives index of the particle in original MCparticle array
	  MCElectron = IsFromBGEventAOD(trkIndexElectron, pEvent);//check if the particle comes from hijing or from enhanced event
	  if(!MCElectron){
	    fGenerator=2;
	    if(fUseHIJINGOnly){return 0;}
	  }else{fGenerator=1;}
	}
    }
  }
  return 1;
}

void AliDxHFEParticleSelectionMCEl::Clear(const char* option)
{
  /// clear internal memory
  fMCTools.Clear(option);
}

AliVParticle *AliDxHFEParticleSelectionMCEl::CreateParticle(AliVParticle* track)
{
  //AliReducedParticle *part = new AliReducedParticle(track->Eta(), track->Phi(), track->Pt(),track->Charge(),fOriginMother,fGenerator);
  AliReducedParticle *part = new AliReducedParticle(track->Eta(), track->Phi(), track->Pt(),track->Charge(),fOriginMother);
  return part;

}

int AliDxHFEParticleSelectionMCEl::ParseArguments(const char* arguments)
{
  // parse arguments and set internal flags
  TString strArguments(arguments);
  unique_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
  if (!tokens.get()) return 0;

  AliInfo(strArguments);
  TIter next(tokens.get());
  TObject* token;
  while ((token=next())) {
    TString argument=token->GetName();
    if(argument.CompareTo("usekine")==0){
      fUseKine=kTRUE;
      fRemoveSecondary=kTRUE;
      continue;
    }
    if(argument.CompareTo("usekine-keepsecondary")==0){
      fUseKine=kTRUE;
      fRemoveSecondary=kFALSE;
      continue;
    }
    if (argument.BeginsWith("elmcreco")){
      fUseMCReco=kTRUE;
      if(argument.BeginsWith("elmcreco=")){
	argument.ReplaceAll("elmcreco=", "");
	if(argument.CompareTo("alltracks")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kNoCuts;
	else if(argument.CompareTo("afterreckineitstpc")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kRecKineITSTPC;
	else if(argument.CompareTo("afterrecprim")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kRecPrim;
	else if(argument.CompareTo("afterhfeits")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kHFEcutsITS;
	else if(argument.CompareTo("afterhfetof")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kHFEcutsTOF;
	else if(argument.CompareTo("afterhfetpc")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kHFEcutsTPC;
	else if(argument.CompareTo("aftertrackcuts")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kHFEcutsTPC;
	else if(argument.CompareTo("aftertofpid")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kPIDTOF;
	else if(argument.CompareTo("aftertpcpid")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kPIDTPC;
	else if(argument.CompareTo("afterfullpid")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kPIDTOFTPC;
	else if(argument.CompareTo("afterpirej")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kRejPi;
	else if(argument.CompareTo("afterprotonrej")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kRejProton;
	else if(argument.CompareTo("afterinvmass")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kINVMASS;
	else if(argument.CompareTo("afterselection")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kSelected;
	else AliFatal(Form("unknown argument '%s'", argument.Data()));

	AliDxHFEParticleSelectionEl::SetFinalCutStep(fSelectionStep);
      }
      continue;
    }
    if(argument.BeginsWith("storeonlyMCelectrons")){
      AliInfo("Store only MC truth electrons");
      fStoreOnlyMCElectrons=kTRUE;
      continue;
    }
    if(argument.BeginsWith("storelastcutstep")){
      AliInfo("Stores the last cut step");
      fUseMCReco=kTRUE;
      fStoreCutStepInfo=kTRUE;
      AliDxHFEParticleSelectionEl::SetStoreLastCutStep(kTRUE);
      continue;
    }
    if(argument.BeginsWith("removeEfromD0")){
      AliInfo("Removing electrons decaying from D0");
      fRemoveEfromD0=kTRUE;
      continue;
    }
    if(argument.BeginsWith("ElSelection=")){
      argument.ReplaceAll("ElSelection=","");
      if(argument.CompareTo("hadron")==0){ fElSelection=kHadron;}
      else if(argument.CompareTo("nonHFE")==0) {fElSelection=kNonHFE;}
      else if(argument.CompareTo("HFE")==0){ fElSelection=kHFE;}
      else if(argument.CompareTo("Onlyc")==0){ fElSelection=kOnlyc;}
      else if(argument.CompareTo("Onlyb")==0){ fElSelection=kOnlyb;}
      else AliFatal(Form("unknown argument '%s'", argument.Data()));
      AliInfo(Form("Selecting only source %d",fElSelection));
      continue;
    }
    if (argument.BeginsWith("notusePhysPrim")){
      AliInfo("Not Use IsPhysicalPrimary()");
      fRemoveSecondary=kFALSE;
      continue;
    }
    if (argument.BeginsWith("onlyPhysPrim")){
      AliInfo("Use IsPhysicalPrimary()");
      fRemoveSecondary=kTRUE;
      continue;
    }
    if(argument.BeginsWith("mc-only")){
      AliInfo("Do only test on MC info");
      fMCInfo=kMCOnly;
      continue;
    }
    if(argument.BeginsWith("mc-first")){
      AliInfo("Do test on MC info first");
      fMCInfo=kMCFirst;
      continue;
    }
    if(argument.BeginsWith("mc-last")){
      AliInfo("Do test on MC info last");
      fMCInfo=kMCLast;
      continue;
    }
    if(argument.BeginsWith("usegenerator")){
      AliInfo("Select source also based on generator");
      fUseGenerator=kTRUE;	    
      continue;
    }
    if(argument.BeginsWith("useHIJINGOnly")){
      AliInfo("Select source also based on generator");
      fUseHIJINGOnly=kTRUE;
      fUseGenerator=kTRUE;
      continue;
    }
    if (argument.BeginsWith("system=")) {
      argument.ReplaceAll("system=", "");
      if (argument.CompareTo("pp")==0) {fSystem=0;}
      else if (argument.CompareTo("Pb-Pb")==0) {fSystem=1;}
      else if (argument.CompareTo("p-Pb")==0) {fSystem=2;}
      else {
	AliWarning(Form("can not set collision system, unknown parameter '%s'", argument.Data()));
	// TODO: check what makes sense
	fSystem=0;
      }
      continue;
    }
    // forwarding of single argument works, unless key-option pairs separated
    // by blanks are introduced
    AliDxHFEParticleSelection::ParseArguments(argument);
  }
  
  return 0;
}

//___________________________________________________________________________
Bool_t AliDxHFEParticleSelectionMCEl::IsFromBGEventAOD(Int_t Index,  const AliVEvent* aodEvent)
{
  //    AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(event);
    //Check if the particle is from Hijing or Enhanced event
    AliAODMCHeader *mcHeader;
    Int_t fNBG =-1;
    
    mcHeader = dynamic_cast<AliAODMCHeader*>(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!mcHeader) {
        AliError("Could not find MC Header in AOD");
        return (0);
    }
    
    TList *List = mcHeader->GetCocktailHeaders();
    //        List->Print();
    AliGenHijingEventHeader* hijingH = dynamic_cast<AliGenHijingEventHeader*>(List->FindObject("Hijing UE_1"));
    if (!hijingH){
        AliError("no GenHijing header");
        return (0);
    }
    fNBG = hijingH->NProduced();
    return (Index < fNBG);
}
Int_t AliDxHFEParticleSelectionMCEl::GetPrimary(Int_t id, TClonesArray *mcArray){
    
    // Return number of primary that has generated track
    int current, parent;
    parent=id;
    while (1) {
        current=parent;
        AliAODMCParticle *Part = (AliAODMCParticle*)mcArray->At(current);
        parent=Part->GetMother();
        //  cout << "GetPartArr momid :"  << parent << endl;
        if(parent<0) return current;
    }
}
