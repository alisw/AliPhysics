
//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Sedat Altinpinar <Sedat.Altinpinar@cern.ch>           *
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

// $Id$

/// @file   AliDxHFECorrelationMC.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-11-02
/// @brief  Worker class for DxHFE correlation
///

#include "AliDxHFECorrelationMC.h"
#include "AliVParticle.h"
#include "AliLog.h"
//#include "AliAnalysisCuts.h"         // required dependency libANALYSISalice.so
//#include "AliFlowTrackSimple.h"      // required dependency libPWGflowBase.so
//#include "AliFlowCandidateTrack.h"   // required dependency libPWGflowTasks.so
//#include "AliCFContainer.h"          // required dependency libCORRFW.so
//#include "AliAODRecoDecayHF2Prong.h"   // libPWGHFvertexingHF
#include "TObjArray.h"
#include "AliDxHFECorrelation.h"
#include "AliHFAssociatedTrackCuts.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "AliReducedParticle.h"
#include "AliDxHFEParticleSelection.h"
#include "AliDxHFEParticleSelectionMCD0.h"
#include "AliDxHFEParticleSelectionMCEl.h"
#include "AliDxHFEToolsMC.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

ClassImp(AliDxHFECorrelationMC)

AliDxHFECorrelationMC::AliDxHFECorrelationMC(const char* name)
  : AliDxHFECorrelation(name?name:"AliDxHFECorrelationMC")
  , fMCEventType(0)
  , fStoreOriginEl(kAll)
  , fStoreOriginD(kAll)
  , fRunMode(kReducedModeFullMCInfo)
  , fSystem(0)
  , fUseReducedOrigin(kFALSE)
  , fReducedOriginEl(0)
  , fReducedOriginD0(0)
  , fStorePoolbin(kFALSE)
{
  // default constructor
  // 
  //

}


AliDxHFECorrelationMC::~AliDxHFECorrelationMC()
{
  // destructor
  //
  //
}

THnSparse* AliDxHFECorrelationMC::DefineTHnSparse()
{
  //
  // Defines the THnSparse. 
  if(fRunMode==kFullMode)
    AliDebug(1, "Creating Full Corr THnSparse");
  else
    AliDebug(1, "Creating Reduced Corr THnSparse");
  // here is the only place to change the dimension
  static const int sizeEventdphi = 10;  
  static const int sizeEventdphiReduced = 5;  
  static const int sizeEventdphiRedMC = 8;  
  static const int sizeEventdphiMCpPb = 8;  
  static const int sizeEventdphiRedpPb = 6;  
  static const int sizeEventdphiRedpPbnoPb = 5;  
  const double pi=TMath::Pi();
  Double_t minPhi=GetMinPhi();
  Double_t maxPhi=GetMaxPhi();
  THnSparse* thn=NULL;
  int nrMotherEl=AliDxHFEToolsMC::kNrOrginMother+AliDxHFEParticleSelectionMCEl::kNrBackground;

  cout << "nrMotherEl:  " << nrMotherEl << "     " <<AliDxHFEToolsMC::kNrOrginMother << "  " <<AliDxHFEParticleSelectionMCEl::kNrBackground<< endl;
  // TODO: Everything here needed for eventmixing? 
  TString name;
  name.Form("%s info", GetName());

  if(fRunMode==kFullMode){
    InitTHnSparseArray(sizeEventdphi);

    // 			                        0        1     2      3      4      5       6      7    8            9
    // 			                      D0invmass  PtD0 PhiD0 PtbinD0 Pte   dphi    dEta OrigD0 origEl      process
    int         binsEventdphi[sizeEventdphi] = {   200,    100,  100,  21,   100,  64,     100,   10,  nrMotherEl+1 ,  100 };
    double      minEventdphi [sizeEventdphi] = { 1.5648,     0,    0,   0,     0, minPhi, -2.0, -1.5,   -1.5,       -0.5 };
    double      maxEventdphi [sizeEventdphi] = { 2.1648,    50, 2*pi,  20,    10, maxPhi,  2.0,  8.5,  nrMotherEl-0.5, 99.5 };
    const char* nameEventdphi[sizeEventdphi] = {
      "D0InvMass",
      "PtD0",
      "PhiD0",
      "PtBinD0",
      "PtEl",
      "#Delta#Phi",
      "#Delta#eta", 
      "Origin D0", 
      "Origin Electron",
      "Original Process"
    };
    thn=(THnSparse*)CreateControlTHnSparse(name,sizeEventdphi,binsEventdphi,minEventdphi,maxEventdphi,nameEventdphi);
  }
  else if(fRunMode==kReducedModeFullMCInfo){
    if(2==fSystem){//Reduced bins for pPb due to memory consumption
      InitTHnSparseArray(sizeEventdphiMCpPb);
      // 			                                 0        1    2     3        4     5           6       7              8
      // 			                             D0invmass  PtD0   Pte   dphi    dEta poolbin   OrigD0   origEl        generatorEl
      int         binsEventdphiRedMC[sizeEventdphiMCpPb] = {   150,      28,    50,  16,     20,    6,      10,  nrMotherEl+1};//,       4 };
      double      minEventdphiRedMC [sizeEventdphiMCpPb] = { 1.5848,      2,     0, minPhi, -2.0,  -0.5,  -1.5,    -1.5};//,          -1.5 };
      double      maxEventdphiRedMC [sizeEventdphiMCpPb] = { 2.1848,     16,    10, maxPhi,  2.0,   5.5,    8.5,  nrMotherEl-0.5};//,   2.5 };
      const char* nameEventdphiRedMC[sizeEventdphiMCpPb] = {
	"D0InvMass",
	"PtD0",
	"PtEl",
	"#Delta#Phi",
	"#Delta#eta", 
	"PoolBin",
	"Origin D0", 
	"Origin Electron"
      };
      thn=(THnSparse*)CreateControlTHnSparse(name,sizeEventdphiMCpPb,binsEventdphiRedMC,minEventdphiRedMC,maxEventdphiRedMC,nameEventdphiRedMC);
    }
    else{
      InitTHnSparseArray(sizeEventdphiRedMC);
      // 			                                       0        1    2     3      4      5       6          7
      // 			                                  D0invmass  PtD0   Pte   dphi    dEta OrigD0 origEl    generatorEl
      int         binsEventdphiRedMC[sizeEventdphiRedMC] = {   200,    100,   100,  64,     100,   10, nrMotherEl+1,    4 };
      double      minEventdphiRedMC [sizeEventdphiRedMC] = { 1.5648,     0,     0, minPhi, -2.0, -1.5,   -1.5,       -1.5 };
      double      maxEventdphiRedMC [sizeEventdphiRedMC] = { 2.1648,    50,    10, maxPhi,  2.0,  8.5, nrMotherEl-0.5,2.5 };
      const char* nameEventdphiRedMC[sizeEventdphiRedMC] = {
	"D0InvMass",
	"PtD0",
	"PtEl",
	"#Delta#Phi",
	"#Delta#eta", 
	"Origin D0", 
	"Origin Electron",
      };
      thn=(THnSparse*)CreateControlTHnSparse(name,sizeEventdphiRedMC,binsEventdphiRedMC,minEventdphiRedMC,maxEventdphiRedMC,nameEventdphiRedMC);
    }
  }
  else{//Reduced mode
      if(2==fSystem){//Tailored binning for pPb

	if(fStorePoolbin){
	  InitTHnSparseArray(sizeEventdphiRedpPb);
	  // 			                                 0        1    2     3        4     5           6       7              8
	  // 			                             D0invmass  PtD0   Pte   dphi    dEta poolbin   OrigD0   origEl        generatorEl
	  int         binsEventdphiRedMC[sizeEventdphiRedpPb] = {   150,      28,    50,  16,     20,    6};//,      10,  nrMotherEl+1};//,       4 };
	  double      minEventdphiRedMC [sizeEventdphiRedpPb] = { 1.5848,      2,     0, minPhi, -2.0,  -0.5};//,  -1.5,    -1.5};//,          -1.5 };
	  double      maxEventdphiRedMC [sizeEventdphiRedpPb] = { 2.1848,     16,    10, maxPhi,  2.0,   5.5};//,    8.5,  nrMotherEl-0.5};//,   2.5 };
	  const char* nameEventdphiRedMC[sizeEventdphiRedpPb] = {
	    "D0InvMass",
	    "PtD0",
	    "PtEl",
	    "#Delta#Phi",
	    "#Delta#eta", 
	    "PoolBin"
	  };
	  thn=(THnSparse*)CreateControlTHnSparse(name,sizeEventdphiRedpPb,binsEventdphiRedMC,minEventdphiRedMC,maxEventdphiRedMC,nameEventdphiRedMC);
	}else{
	  InitTHnSparseArray(sizeEventdphiRedpPbnoPb);
	  // 			                                 0        1    2     3        4   
	  // 			                             D0invmass  PtD0   Pte   dphi    dEta 
	  int         binsEventdphiRedMC[sizeEventdphiRedpPbnoPb] = {   150,      28,    50,  16,     20};
	  double      minEventdphiRedMC [sizeEventdphiRedpPbnoPb] = { 1.5848,      2,     0, minPhi, -2.0};
	  double      maxEventdphiRedMC [sizeEventdphiRedpPbnoPb] = { 2.1848,     16,    10, maxPhi,  2.0};
	  const char* nameEventdphiRedMC[sizeEventdphiRedpPbnoPb] = {
	    "D0InvMass",
	    "PtD0",
	    "PtEl",
	    "#Delta#Phi",
	    "#Delta#eta" 
	  };
	  thn=(THnSparse*)CreateControlTHnSparse(name,sizeEventdphiRedpPbnoPb,binsEventdphiRedMC,minEventdphiRedMC,maxEventdphiRedMC,nameEventdphiRedMC);
	}
      }else{
	InitTHnSparseArray(sizeEventdphiReduced);
	//           			                    0        1      2    3      4  
	// 		                	                 D0invmass  PtD0  Pte   dphi    dEta
	int         binsEventdphiRed[sizeEventdphiReduced] = {   200,    100,  100,   64,    100 };
	double      minEventdphiRed [sizeEventdphiReduced] = { 1.5648,     0,   0, minPhi, -2.0 };
	double      maxEventdphiRed [sizeEventdphiReduced] = { 2.1648,    50,  10, maxPhi,  2.0 };
	const char* nameEventdphiRed[sizeEventdphiReduced] = {
	  "D0InvMass",
	  "PtD0",
	  "PtEl",
	  "#Delta#Phi",
	  "#Delta#eta" 
	};
	thn=(THnSparse*)CreateControlTHnSparse(name,sizeEventdphiReduced,binsEventdphiRed,minEventdphiRed,maxEventdphiRed,nameEventdphiRed);
      }	
  }
  return thn;
      

}

Bool_t AliDxHFECorrelationMC::TestParticle(AliVParticle* p, Int_t id){

  AliReducedParticle *part=(AliReducedParticle*)p;
  if (!part) return -ENODATA;

  Bool_t selected = kTRUE;
  Bool_t isCharm=(part->GetOriginMother()==AliDxHFEToolsMC::kOriginCharm || 
  		  part->GetOriginMother()==AliDxHFEToolsMC::kOriginGluonCharm);
  Bool_t isBeauty=(part->GetOriginMother()==AliDxHFEToolsMC::kOriginBeauty || 
		   part->GetOriginMother()==AliDxHFEToolsMC::kOriginGluonBeauty);

  if(id==kD){
    if(isCharm) fReducedOriginD0=kOriginC;
    if(isBeauty) fReducedOriginD0=kOriginB;
    if(!(isCharm || isBeauty)) fReducedOriginD0=kOriginNonHF;
    if(!(fStoreOriginD==kAll)){
      if(fStoreOriginD==kC ){
	// Test if particle really is from C
	if(!isCharm)
	  selected =kFALSE;
      }
      else if(fStoreOriginD==kB ){
	// Test if particle really is from B
	if(!isBeauty)
	  selected =kFALSE;
      }
      else if(fStoreOriginD==kHF ){
	// Test if particle really is HF
	if(!(isCharm || isBeauty))
	  selected=kFALSE;
      }
    }
  }
  //TODO add for background from HF
  // Test to see if test for D/el and whether to do further selection
  if(id==kElectron){
    if(isCharm) fReducedOriginEl=kOriginC; //c and gluon-to-c
    if(isBeauty) fReducedOriginEl=kOriginB; //b and gluon-to-b
    if((part->GetOriginMother() >= AliDxHFEToolsMC::kOriginNone && part->GetOriginMother()<=AliDxHFEToolsMC::kOriginStrange)
       || (part->GetOriginMother() > AliDxHFEToolsMC::kOriginGluonBeauty )) fReducedOriginEl=kOriginNonHF; //NonHF
    if(part->GetOriginMother()<AliDxHFEToolsMC::kOriginNone) fReducedOriginEl=kOriginHadron; //Hadrons [Warning] This might overwrite NonHF if some cases.
    
    if(!(fStoreOriginEl==kAll)){
      if(fStoreOriginEl==kC){ // Test if particle really is from C
	if(!isCharm)
	  selected =kFALSE;
      }
      else if(fStoreOriginEl==kB){ // Test if particle really is from B

	if(!isBeauty)
	  selected =kFALSE;
      }
      else if(fStoreOriginEl==kHF){ // Test if particle really is HF	
	if((!isCharm) || (!isBeauty))
	  selected=kFALSE;
      }
  
      // Test extra for source of el
      if(fStoreOriginEl==kNonHF){
	// Test if particle really is from NonHF
	if(!((part->GetOriginMother() >= AliDxHFEToolsMC::kOriginNone && part->GetOriginMother()<=AliDxHFEToolsMC::kOriginStrange)
	     || (part->GetOriginMother() > AliDxHFEToolsMC::kOriginGluonBeauty )))
	  selected =kFALSE;
      }
      else if(fStoreOriginEl==kHadrons){
	// Test if particle really is from hadrons
	if(! (part->GetOriginMother()<AliDxHFEToolsMC::kOriginNone))
	  selected =kFALSE;
      }
    }
  }

  return selected;

}

int AliDxHFECorrelationMC::FillParticleProperties(AliVParticle* tr, AliVParticle *as, Double_t* data, int dimension) const
{
  // fill the data array from the particle data
  if (!data) return -EINVAL;
  AliReducedParticle *ptrigger=(AliReducedParticle*)tr;
  AliReducedParticle *assoc=(AliReducedParticle*)as;
  if (!ptrigger || !assoc) return -ENODATA;
  int i=0;
  if (dimension!=GetDimTHnSparse()) {
    // TODO: think about filling only the available data and throwing a warning
    return -ENOSPC;
  }
  if(AliDxHFECorrelation::GetTriggerParticleType()==kD){
    data[i++]=ptrigger->GetInvMass();
    data[i++]=ptrigger->Pt();
    if(fRunMode==kFullMode)data[i++]=ptrigger->Phi();
    if(fRunMode==kFullMode) {data[i++]=ptrigger->GetPtBin(); }
    data[i++]=assoc->Pt();
  } 
  else{
    data[i++]=assoc->GetInvMass();
    data[i++]=assoc->Pt();
    data[i++]=assoc->Phi();
    if(fRunMode==kFullMode) data[i++]=assoc->GetPtBin(); 
    data[i++]=ptrigger->Pt();
  }
  data[i++]=AliDxHFECorrelation::GetDeltaPhi();
  data[i++]=AliDxHFECorrelation::GetDeltaEta();
  if(fSystem==2 && fStorePoolbin)data[i++]=AliDxHFECorrelation::GetPoolBin(); //Remove "fSystem==2" if to be used with other systems. Only pPb has poolbin storage enabled for now, so this works as a double safety
  if(fRunMode==kFullMode || fRunMode==kReducedModeFullMCInfo){
    if(AliDxHFECorrelation::GetTriggerParticleType()==kD){
      if(!fUseReducedOrigin){data[i++]=ptrigger->GetOriginMother();}else{data[i++]=fReducedOriginD0;}
      if(!fUseReducedOrigin){data[i++]=assoc->GetOriginMother();}else{data[i++]=fReducedOriginEl;}
    }
    else {
      if(!fUseReducedOrigin){data[i++]=assoc->GetOriginMother();}else{data[i++]=fReducedOriginD0;}
      if(!fUseReducedOrigin){data[i++]=ptrigger->GetOriginMother();}else{data[i++]=fReducedOriginEl;}
    }
    if(AliDxHFECorrelation::GetTriggerParticleType()==kD){
      if(fSystem!=2) data[i++]=assoc->GetGeneratorIndex();//[FIXME] This should be "GetGenerator()" once the changes in AliReducedParticle are in place
    }
    else {
      if(fSystem!=2) data[i++]=ptrigger->GetGeneratorIndex();//[FIXME] This should be "GetGenerator()" once the changes in AliReducedParticle are in place
    }
  }
  if(fRunMode==kFullMode ) data[i++]=fMCEventType;
  return i;
}

int AliDxHFECorrelationMC::Fill(const TObjArray* triggerCandidates, TObjArray* associatedTracks, const AliVEvent* pEvent)
{
  // TODO: Implement more on MC?? (Needed?)
  return AliDxHFECorrelation::Fill(triggerCandidates,associatedTracks,pEvent);
}

double AliDxHFECorrelationMC::GetD0Eff(AliVParticle* tr, Double_t evMult){

  Double_t D0eff=1;
  AliReducedParticle *track=(AliReducedParticle*)tr;
  if (!track) return -ENODATA;
  Double_t pt=track->Pt();
  Double_t origin=track->GetOriginMother();

  //  Bool_t isCharm=(origin==AliDxHFEToolsMC::kOriginCharm || 
  //		  origin==AliDxHFEToolsMC::kOriginGluonCharm);
  Bool_t isBeauty=(origin==AliDxHFEToolsMC::kOriginBeauty || 
		   origin==AliDxHFEToolsMC::kOriginGluonBeauty);

  AliHFAssociatedTrackCuts* cuts=dynamic_cast<AliHFAssociatedTrackCuts*>(GetCuts());
  if (!cuts) {
    if (GetCuts())
      AliError(Form("cuts object of wrong type %s, required AliHFAssociatedTrackCuts", cuts->ClassName()));
    else
      AliError("mandatory cuts object missing");
    return -EINVAL;
  }

  if(isBeauty)
    D0eff=cuts->GetTrigWeightB(pt,evMult);
  else
    D0eff=cuts->GetTrigWeight(pt,evMult);

  return D0eff;
}


int AliDxHFECorrelationMC::ParseArguments(const char* arguments)
{
  // parse arguments and set internal flags
  TString strArguments(arguments);
  unique_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
  if (!tokens.get()) return -ENOMEM;
  TIter next(tokens.get());
  TObject* token;
  while ((token=next())) {
    TString argument=token->GetName();
    if (argument.BeginsWith("storeoriginD=")){
      argument.ReplaceAll("storeoriginD=", "");
      if (argument.CompareTo("all")==0) { fStoreOriginD=kAll; AliInfo("Store all Ds"); }
      else if (argument.CompareTo("charm")==0) { fStoreOriginD=kC; AliInfo("Store only D from charm");}
      else if (argument.CompareTo("beauty")==0){ fStoreOriginD=kB; AliInfo("Store only D from beauty");}
      else if (argument.CompareTo("HF")==0){ fStoreOriginD=kHF; AliInfo("Store D from HF");}
      continue;
    }  
    if (argument.BeginsWith("storeoriginEl=")){
      argument.ReplaceAll("storeoriginEl=", "");
      if (argument.CompareTo("all")==0) { fStoreOriginEl=kAll; AliInfo("Store all electrons"); }
      else if (argument.CompareTo("charm")==0) { fStoreOriginEl=kC; AliInfo("Store only electrons from charm");}
      else if (argument.CompareTo("beauty")==0){ fStoreOriginEl=kB; AliInfo("Store only electrons from beauty");}
      else if (argument.CompareTo("HF")==0){ fStoreOriginEl=kHF; AliInfo("Store electrons from HF");}
      else if (argument.CompareTo("nonHF")==0){ fStoreOriginEl=kNonHF; AliInfo("Store electrons from nonHF");}
      else if (argument.CompareTo("hadrons")==0){ fStoreOriginEl=kHadrons; AliInfo("Store electrons candidates from hadrons");}
      continue;
    }  
    if ((argument.CompareTo("reducedMode")==0) || (argument.CompareTo("reducedmode")==0)){
      fRunMode=kReducedMode;
      AliInfo("Running in Reduced mode");
      continue;
    }
    if (argument.BeginsWith("FullMode") || argument.BeginsWith("fullmode")){
      fRunMode=kFullMode;
      AliInfo("Running in Full mode");
      continue;
    }
    if (argument.BeginsWith("reducedModewithMC") || argument.BeginsWith("reducedmodewithmc")){
      fRunMode=kReducedModeFullMCInfo;
      AliInfo("Running in Reduced mode with MC info stored as well");
      continue;
    }
    if (argument.BeginsWith("storeReducedOrigin")){
      fUseReducedOrigin=kTRUE;
      AliInfo("Storing reduced origin information for electrons and D0 in correlation sparse");
      continue;
    }
    if (argument.BeginsWith("StorePoolbin")){
      fStorePoolbin=kTRUE;
      AliInfo("Storing poolbin information");
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
    //    AliWarning(Form("unknown argument '%s'", argument.Data()));
    AliDxHFECorrelation::ParseArguments(argument);      
  }

  return 0;
}
