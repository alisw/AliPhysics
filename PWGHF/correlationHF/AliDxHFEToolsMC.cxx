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

/// @file   AliDxHFEToolsMC.cxx
/// @author Hege Erdal, Matthias Richter
/// @date   2012-07-19
/// @brief  Common Tools for MC particle selection
///

#include "AliDxHFEToolsMC.h"
#include "AliAODMCParticle.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "TObjArray.h"
#include "TString.h"
#include "TH1D.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEToolsMC)

AliDxHFEToolsMC::AliDxHFEToolsMC(const char* option)
  : fSequence(kMCLast)
  , fMCParticles(NULL)
  , fPDGs()
  , fMotherPDGs()
  , fHistPDG(NULL)
  , fHistPDGMother(NULL)
  , fOriginMother(kOriginNone)
  , fMClabel(-1)
  , fNrMCParticles(-1)
  , fUseKine(kFALSE)
{
  // constructor
  // 
  // 
  // 
  // 
  Init(option);
}

const char*  AliDxHFEToolsMC::fgkPDGBinLabels[]={
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

const char*  AliDxHFEToolsMC::fgkPDGMotherBinLabels[]={
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

const char*  AliDxHFEToolsMC::fgkStatisticsBinLabels[]={
  "all",   // all MC particles
  "found", // selected particles with correct MC
  "fake"   // selected particles without corresponding MC
};

AliDxHFEToolsMC::~AliDxHFEToolsMC()
{
  // destructor
  if (fHistPDG) delete fHistPDG;
  fHistPDG=NULL;
  if (fHistPDGMother) delete fHistPDGMother;
  fHistPDGMother=NULL;
}

int AliDxHFEToolsMC::Init(const char* option)
{
  // initialize according to options
  TString strOption(option);
  bool bControlHist=true;
  std::auto_ptr<TObjArray> tokens(strOption.Tokenize(" "));
  if (tokens.get() && tokens->GetEntriesFast()>0) {
    for (int itoken=0; itoken<tokens->GetEntriesFast(); itoken++) {
      if (tokens->At(itoken)==NULL) continue;
      TString arg=tokens->At(itoken)->GetName();
      const char* key="";
      key="pdg=";
      if (arg.BeginsWith(key)) {
	arg.ReplaceAll(key, "");
	fPDGs.push_back(arg.Atoi());
      }
      key="mother-pdg=";
      if (arg.BeginsWith(key)) {
	arg.ReplaceAll(key, "");
	fMotherPDGs.push_back(arg.Atoi());
      }
      key="control-hist=";
      if (arg.BeginsWith(key)) {
	arg.ReplaceAll(key, "");
	bControlHist=arg.CompareTo("off")!=0;
      }
      key="mc-first";
      if (arg.BeginsWith(key)) {
	fSequence=kMCFirst;	
      }
      key="mc-last";
      if (arg.BeginsWith(key)) {
	fSequence=kMCLast;	
      }
      key="usekine";
      if(arg.BeginsWith(key)) {
	printf("AliDxHFEToolsMC::Init()  Using Kinematical\n");
	fUseKine=kTRUE;
      }
    }
  }
  
  if (bControlHist) {
    fHistPDG=CreateControlHistogram("histPDG",
				    "pdg code of selected particle",
				    sizeof(fgkPDGBinLabels)/sizeof(const char*),
				    fgkPDGBinLabels);
    fHistPDGMother=CreateControlHistogram("histPDGMother",
					  "pdg code of first mother of selected particle",
					  sizeof(fgkPDGMotherBinLabels)/sizeof(fgkPDGMotherBinLabels[0]),
					  fgkPDGMotherBinLabels);
  }
  return 0;
}

int AliDxHFEToolsMC::InitMCParticles(const AliVEvent* pEvent)
{

  // init MC info from event object
  if (!pEvent) return -EINVAL;

  // TODO: choose branch name depending on VEvent type; configurable?
  TString branchname(AliAODMCParticle::StdBranchName());
  TObject* o=pEvent->FindListObject(branchname);
  if (!o) {
    AliWarningClass(Form("can not find MC info '%s' in event of type '%s'", branchname.Data(), pEvent->ClassName()));
    return -ENOENT;
  }
  fMCParticles = dynamic_cast<TObjArray*>(o);
  if (!fMCParticles) {
    AliWarningClass(Form("ignoring MC info '%s' of wrong type '%s', expecting TObjArray", branchname.Data(), o->ClassName()));
    return -ENODATA;
  }
  fNrMCParticles=fMCParticles->GetEntriesFast();
  return 0;
}

bool AliDxHFEToolsMC::RejectByPDG(int pdg, const vector<int> &list) const
{
  // check if pdg should be rejected, particle is not rejected
  // if it is in the list, returns always false if list is empty
  if (list.size()==0) return false;
  for (vector<int>::const_iterator i=list.begin();
       i!=list.end(); i++) {
    if (*i==pdg) return false;
  }
  return true;
}

bool AliDxHFEToolsMC::RejectByPDG(AliVParticle* p, bool doStatistics, int* pdgParticleResult)
{
  // check if pdg should be rejected
  // always false if not pdg list is initialized
  if (!p) return false;
  Int_t MClabel = p->GetLabel();
  if(MClabel<0) return 0;

  int pdgPart=-1;
  // TODO: there might be more types of particles to be checked
  AliAODMCParticle* aodmcp=0;
  aodmcp=dynamic_cast<AliAODMCParticle*>(fMCParticles->At(MClabel));
  if (aodmcp)
    pdgPart=TMath::Abs(aodmcp->GetPdgCode());
  if (pdgPart<0) return 0;

  if (pdgParticleResult)
    *pdgParticleResult=pdgPart;

  bool bReject=RejectByPDG(pdgPart, fPDGs);
  if (doStatistics && fHistPDG) {
    // TODO: think about histogramming mode, e.g. histogramming of rejected particles?
    if (!bReject) {
      fHistPDG->Fill(MapPDGLabel(p->PdgCode()));
    }
  }
  return bReject;
}

bool AliDxHFEToolsMC::RejectByMotherPDG(AliVParticle* p, bool doStatistics)
{
  // check if pdg should be rejected by mother
  // always false if not mother pdg list is initialized
  // H: think maybe this is particle specific, and should be moved to PartSelMCEl
  if (!p) return false;
  int motherpdg=FindPdgOriginMother(p);
  // TODO: This should be tweaked. Want to 
  bool bReject=RejectByPDG(motherpdg, fMotherPDGs);
  if (doStatistics && fHistPDGMother) {
    // TODO: think about histogramming mode, e.g. histogramming of rejected particles?
    if (!bReject) {
      fHistPDGMother->Fill(MapPDGMotherLabel(p->PdgCode()));
    }
  }
  return bReject;
}

int AliDxHFEToolsMC::FindMotherPDG(AliVParticle* p, bool bReturnFirstMother)
{
  // Will find and return pdg of the mother, either first or loop down to the
  // initial quark

  // To reset fOriginMother. 
  fOriginMother=kOriginNone;

  if (!p) return false;
  int motherpdg=FindPdgOriginMother(p, bReturnFirstMother);

  return motherpdg;
}

TH1* AliDxHFEToolsMC::CreateControlHistogram(const char* name,
					     const char* title,
					     int nBins,
					     const char** binLabels) const
{
  /// create control histogram
  std::auto_ptr<TH1> h(new TH1D(name, title, nBins, -0.5, nBins-0.5));
  if (!h.get()) return NULL;
  for (int iLabel=0; iLabel<nBins; iLabel++) {
    h->GetXaxis()->SetBinLabel(iLabel+1, binLabels[iLabel]);    
  }
  
  return h.release();
}

int AliDxHFEToolsMC::MapPDGLabel(int pdg) const
{
  /// mapping of pdg code to enum
  switch (pdg) {
  case kPDGelectron  : return kPDGLabelElectron;
  case -kPDGelectron : return kPDGLabelPositron;
  case kPDGmuon	     : return kPDGLabelMuPlus;
  case -kPDGmuon     : return kPDGLabelMuMinus;
  case kPDGpion	     : return kPDGLabelPiPlus;
  case -kPDGpion     : return kPDGLabelPiMinus;
  case kPDGkaon	     : return kPDGLabelKPlus;
  case -kPDGkaon     : return kPDGLabelKMinus;
  case kPDGproton    : return kPDGLabelProton;
  case -kPDGproton   : return kPDGLabelAntiproton;
  default:
    return kPDGLabelOthers;
  }
}

int AliDxHFEToolsMC::MapPDGMotherLabel(int pdg) const
{
  /// mapping of pdg code to enum
  switch (pdg) {
  case kPDGd     : return kPDGMotherLabelD;
  case kPDGu     : return kPDGMotherLabelU;
  case kPDGs     : return kPDGMotherLabelS;
  case kPDGc     : return kPDGMotherLabelC;
  case kPDGb     : return kPDGMotherLabelB;
  case kPDGgluon : return kPDGMotherLabelGluon;
  case kPDGgamma : return kPDGMotherLabelGamma;
  case kPDGpi0   : return kPDGMotherLabelPi0;
  case kPDGeta   : return kPDGMotherLabelEta;
  case kPDGproton: return kPDGMotherLabelProton;
  default:
    return kPDGLabelOthers;
  }
}

int AliDxHFEToolsMC::FindPdgOriginMother(AliVParticle* p, bool bReturnFirstMother) 
{
  // Return the pgd of original mother particle
  // TODO: need also to have specific for D0, electron etc
  // for instance to mark when you have gluon, charm or beauty
  // among the mothers. Or maybe this will be the same anyway?
  // TODO: implement tests on origin, if charm/beauty quark and if
  // they came from gluon. use booleans to set this which can be accessed from 
  // outside? Something like fSequence. 

  if (!p) return kPDGnone;

  Int_t imother=-1;
  Int_t MClabel=0;

  // Either using MClabel set from outside (needed for Dmesons), or find it
  if(fMClabel<0){
    MClabel = p->GetLabel();
    if(MClabel<0){
      return kPDGnone;
    }
  }
  else MClabel=fMClabel;

  // try different classes, unfortunately there is no common base class
  AliAODMCParticle* aodmcp=0;
  if(fUseKine)
    aodmcp=dynamic_cast<AliAODMCParticle*>(p);
  else
    aodmcp=dynamic_cast<AliAODMCParticle*>(fMCParticles->At(MClabel));

  if (!aodmcp) {
    return kPDGnone;
  }
  imother = aodmcp->GetMother();

  // Absolute or +/- on pgd ???
  Int_t pdg=TMath::Abs(aodmcp->GetPdgCode());

  if (imother<0){
    // also check this particle
    CheckOriginMother(pdg);
    return pdg;
  }

  if (!fMCParticles->At(imother)) {
    AliErrorClass(Form("no mc particle with label %d", imother));
    return pdg;
  }

  AliAODMCParticle * mother=dynamic_cast<AliAODMCParticle*>(fMCParticles->At(imother));
  if (!mother) {
    AliErrorClass(Form("mc mother particle of wrong class type %s", fMCParticles->At(imother)->ClassName()));
    return pdg;
  }

  if(bReturnFirstMother){
    return mother->GetPdgCode();
  }

  CheckOriginMother(pdg);
  
  if (mother->GetPdgCode()==kPDGproton){
    // why? - H: This is the proton, can't get further back
    // To be discussed whether to do this a different way
    return pdg;
  }

  //Reset fMClabel to find pdg of mother if looping on D
  fMClabel=-1;
  pdg=FindPdgOriginMother(mother);

  return pdg;  
}

void AliDxHFEToolsMC::CheckOriginMother(int pdg)
{

  // Checking if the particle is a quark or gluon and setting fOriginMother accordingly.
  // Right now only check on quark. Need also check on whether it comes from D or B meson?

  switch(pdg){
  case(kPDGc):
    fOriginMother = kOriginCharm; break;
  case(kPDGb): 
    fOriginMother = kOriginBeauty; break;
  case(kPDGgluon): 
    if(fOriginMother==kOriginCharm) fOriginMother=kOriginGluonCharm;
    else if(fOriginMother==kOriginBeauty) fOriginMother=kOriginGluonBeauty;
    else fOriginMother=kOriginGluon;
    break;
  case(kPDGd):
    if(!TestIfHFquark(fOriginMother))
      fOriginMother=kOriginDown; 
    break;
  case(kPDGu):
    if(!TestIfHFquark(fOriginMother))
      fOriginMother=kOriginUp; 
    break;
  case(kPDGs):
    if(!TestIfHFquark(fOriginMother))
      fOriginMother=kOriginStrange; 
    break;
  }
}

Bool_t AliDxHFEToolsMC::TestIfHFquark(int origin)
{

  // Checking if particle has been marked as charm/beauty quark
  
  Bool_t test=kFALSE;
  switch(origin){
  case(kOriginCharm):
    test=kTRUE; break;
  case(kOriginBeauty): 
    test=kTRUE; break;
  case(kOriginGluonCharm): 
    test=kTRUE; break;
  case(kOriginGluonBeauty): 
    test=kTRUE; break;
  }
  return test; 
}

Bool_t AliDxHFEToolsMC::TestMotherHFMeson(int pdg)
{

  // Checking if the pdg corresponds to a HF meson (D or B meson)

  Bool_t isD = kFALSE;
  Bool_t isB = kFALSE;
  if((pdg>=400 && pdg <500) || (pdg>=4000 && pdg<5000 )) 
    isD = kTRUE;
  if((pdg>=500 && pdg <600) || (pdg>=5000 && pdg<6000 )) 
    {isD = kFALSE; isB = kTRUE;}

  if(isD || isB) return kTRUE;
  else return kFALSE;

}




void AliDxHFEToolsMC::Clear(const char* /*option*/)
{
  // clear internal memory
  fMCParticles=NULL;
  fNrMCParticles=-1;
}


int AliDxHFEToolsMC::CheckMCParticle(AliVParticle* p){

  // Checks if MC particle is desired particle

  AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(p);
  if (!mcPart) {
    AliInfoClass("MC Particle not found in tree, skipping"); 
    return -1;
  }
			
  Int_t PDG =TMath::Abs(mcPart->PdgCode()); 
  int bReject=RejectByPDG(PDG, fPDGs);

  return bReject;

}
