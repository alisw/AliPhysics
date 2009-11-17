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

#include <TROOT.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TList.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TArrayI.h>
#include <TF1.h>

#include <AliHeader.h>
#include <AliStack.h>
#include <AliLog.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliESDEvent.h>
#include <AliMCEvent.h>
#include <AliESDVertex.h>
#include <AliVertexerTracks.h>

#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenCocktailEventHeader.h>
#include <AliGenDPMjetEventHeader.h>

#include <AliMathBase.h>
#include <AliESDtrackCuts.h>
#include "dNdPt/AlidNdPtEventCuts.h"
#include "dNdPt/AlidNdPtAcceptanceCuts.h"
#include "dNdPt/AlidNdPtHelper.h"

//____________________________________________________________________
ClassImp(AlidNdPtHelper)

Int_t AlidNdPtHelper::fgLastProcessType = -1;

//____________________________________________________________________
Bool_t AlidNdPtHelper::IsEventTriggered(const AliESD* aEsd, Trigger trigger)
{
  // see function with ULong64_t argument

  ULong64_t triggerMask = aEsd->GetTriggerMask();
  return IsEventTriggered(triggerMask, trigger);
}

//____________________________________________________________________
Bool_t AlidNdPtHelper::IsEventTriggered(ULong64_t triggerMask, Trigger trigger)
{
  // check if the event was triggered
  //
  // this function needs the branch fTriggerMask
  
  // definitions from p-p.cfg
  ULong64_t spdFO = (1 << 14);
  ULong64_t v0left = (1 << 11);
  ULong64_t v0right = (1 << 12);

  switch (trigger)
  {
    case kMB1:
    {
      if (triggerMask & spdFO || ((triggerMask & v0left) || (triggerMask & v0right)))
        return kTRUE;
      break;
    }
    case kMB2:
    {
      if (triggerMask & spdFO && ((triggerMask & v0left) || (triggerMask & v0right)))
        return kTRUE;
      break;
    }
    case kSPDFASTOR:
    {
      if (triggerMask & spdFO)
        return kTRUE;
      break;
    }
  }

  return kFALSE;
}

//____________________________________________________________________
Bool_t AlidNdPtHelper::TestVertex(const AliESDVertex* vertex, AnalysisMode analysisMode, Bool_t debug)
{
  // Checks if a vertex meets the needed quality criteria
  if(!vertex) return kFALSE;

  Float_t requiredZResolution = -1;
  if (analysisMode == kSPD || analysisMode == kTPCITS || analysisMode == kTPCSPDvtx)
  {
    requiredZResolution = 0.1;
  }
  else if (analysisMode == kTPC || analysisMode == kMCRec || 
           analysisMode == kMCPion || analysisMode == kMCKaon || 
	   analysisMode == kMCProton || analysisMode ==kPlus || analysisMode ==kMinus)
    requiredZResolution = 10.;

  // check Ncontributors
  if (vertex->GetNContributors() <= 0) {
    if (debug){
      Printf("AlidNdPtHelper::GetVertex: NContributors() <= 0: %d",vertex->GetNContributors());
      Printf("AlidNdPtHelper::GetVertex: NIndices(): %d",vertex->GetNIndices());
      vertex->Print();
    }
    return kFALSE;
  }

  // check resolution
  Double_t zRes = vertex->GetZRes();
  if (zRes == 0) {
    Printf("AlidNdPtHelper::GetVertex: UNEXPECTED: resolution is 0.");
    return kFALSE;
  }

  if (zRes > requiredZResolution) {
    if (debug)
      Printf("AlidNdPtHelper::TestVertex: Resolution too poor %f (required: %f", zRes, requiredZResolution);
    return kFALSE;
  }

  return kTRUE;
}

//____________________________________________________________________
const AliESDVertex* AlidNdPtHelper::GetVertex(AliESDEvent* aEsd, AlidNdPtEventCuts *evtCuts, AlidNdPtAcceptanceCuts *accCuts, AliESDtrackCuts *trackCuts, AnalysisMode analysisMode, Bool_t debug, Bool_t bRedoTPC, Bool_t bUseMeanVertex)
{
  // Get the vertex from the ESD and returns it if the vertex is valid
  //
  // Second argument decides which vertex is used (this selects
  // also the quality criteria that are applied)

  if(!aEsd) 
  { 
    ::Error("AlidNdPtHelper::GetVertex()","esd event is NULL");
    return NULL;  
  }
 
  if(!evtCuts || !accCuts || !trackCuts) 
  { 
    ::Error("AlidNdPtHelper::GetVertex()","cuts not available");
    return NULL;  
  }

  const AliESDVertex* vertex = 0;
  AliESDVertex *initVertex = 0;
  if (analysisMode == kSPD || analysisMode == kTPCITS || analysisMode == kTPCSPDvtx ||  
      analysisMode == kMCRec || analysisMode == kMCPion || analysisMode == kMCKaon || 
      analysisMode == kMCProton || analysisMode ==kPlus || analysisMode ==kMinus )
  {
    vertex = aEsd->GetPrimaryVertexSPD();
    if (debug)
      Printf("AlidNdPtHelper::GetVertex: Returning SPD vertex");
  }
  else if (analysisMode == kTPC || analysisMode == kMCRec || 
           analysisMode == kMCPion || analysisMode == kMCKaon || analysisMode == kMCProton || 
	   analysisMode == kPlus || analysisMode == kMinus)
  {
    if(bRedoTPC) {

      Double_t kBz = aEsd->GetMagneticField();
      AliVertexerTracks vertexer(kBz);

      if(bUseMeanVertex) {
	 Double_t pos[3]={evtCuts->GetMeanXv(),evtCuts->GetMeanYv(),evtCuts->GetMeanZv()};
	 Double_t err[3]={evtCuts->GetSigmaMeanXv(),evtCuts->GetSigmaMeanYv(),evtCuts->GetSigmaMeanZv()};
	 //printf("pos[0] %f, pos[1] %f, pos[2] %f \n", pos[0], pos[1], pos[2]);
	 initVertex = new AliESDVertex(pos,err);
	 vertexer.SetVtxStart(initVertex);
	 vertexer.SetConstraintOn();
      }

      //vertexer.SetTPCMode(Double_t dcacut=0.1, Double_t dcacutIter0=1.0, Double_t maxd0z0=5.0, Int_t minCls=10, Int_t mintrks=1, Double_t nsigma=3., Double_t mindetfitter=0.1, Double_t maxtgl=1.5, Double_t fidR=3., Double_t fidZ=30., Int_t finderAlgo=1, Int_t finderAlgoIter0=4);

      Double_t maxDCAr = accCuts->GetMaxDCAr();
      Double_t maxDCAz = accCuts->GetMaxDCAz();
      Int_t minTPCClust = trackCuts->GetMinNClusterTPC();

      vertexer.SetTPCMode(0.1,1.0,5.0,minTPCClust,1,3.,0.1,2.0,maxDCAr,maxDCAz,1,4);

      // TPC track preselection
      Int_t ntracks = aEsd->GetNumberOfTracks();
      TObjArray array(ntracks);
      UShort_t *id = new UShort_t[ntracks];

      Int_t count=0;
      for (Int_t i=0;i <ntracks; i++) {
        AliESDtrack *t = aEsd->GetTrack(i);
        if (!t) continue;
        if (t->Charge() == 0) continue;
        if (!t->GetTPCInnerParam()) continue;
        if (t->GetTPCNcls()<vertexer.GetMinClusters()) continue;
        AliExternalTrackParam  *tpcTrack  = new AliExternalTrackParam(*(t->GetTPCInnerParam()));
	if(tpcTrack) { 
	  array.AddLast(tpcTrack);
	  id[count] = (UShort_t)t->GetID();
	  count++;
	}
      } 
      AliESDVertex *vTPC = vertexer.VertexForSelectedTracks(&array,id, kTRUE, kTRUE, bUseMeanVertex);
      aEsd->SetPrimaryVertexTPC(vTPC);

      for (Int_t i=0; i<aEsd->GetNumberOfTracks(); i++) {
	AliESDtrack *t = aEsd->GetTrack(i);
	t->RelateToVertexTPC(vTPC, kBz, kVeryBig);
      }
      
      delete vTPC;
      array.Delete();
      delete [] id; id=NULL;

    }
    vertex = aEsd->GetPrimaryVertexTPC();
    if (debug)
      Printf("AlidNdPtHelper::GetVertex: Returning vertex from tracks");
    }
      else
       Printf("AlidNdPtHelper::GetVertex: ERROR: Invalid second argument %d", analysisMode);

    if (!vertex) {
     if (debug)
      Printf("AlidNdPtHelper::GetVertex: No vertex found in ESD");
      return 0;
    }

  if (debug)
  {
    Printf("AlidNdPtHelper::GetVertex: Returning valid vertex: %s", vertex->GetTitle());
    vertex->Print();
  }
  
  if(initVertex) delete initVertex; initVertex=NULL;
  return vertex;
}

//____________________________________________________________________
Bool_t AlidNdPtHelper::IsPrimaryParticle(AliStack* stack, Int_t idx, AnalysisMode analysisMode)
{
// check primary particles 
// depending on the analysis mode
//
  if(!stack) return kFALSE;

  TParticle* particle = stack->Particle(idx);
  if (!particle) return  kFALSE;

  // only charged particles
  Double_t charge = particle->GetPDG()->Charge()/3.;
  if (charge == 0.0) return kFALSE;

  Int_t pdg = TMath::Abs(particle->GetPdgCode());

  // physical primary
  Bool_t prim = stack->IsPhysicalPrimary(idx);

  if(analysisMode==kMCPion) {
    if(prim && pdg==kPiPlus) return kTRUE;
    else return kFALSE;
  } 

  if (analysisMode==kMCKaon) {
    if(prim && pdg==kKPlus) return kTRUE;
    else return kFALSE;
  }
    
  if (analysisMode==kMCProton) {
    if(prim && pdg==kProton) return kTRUE;
    else return kFALSE;
  }

return prim;
}

//____________________________________________________________________
Bool_t AlidNdPtHelper::IsPrimaryCharged(TParticle* aParticle, Int_t aTotalPrimaries, Bool_t adebug)
{
  //
  // this function checks if a particle from the event generator (i.e. among the nPrim particles in the stack)
  // shall be counted as a primary particle
  //
  // This function or a equivalent should be available in some common place of AliRoot
  //
  // WARNING: Call this function only for particles that are among the particles from the event generator!
  // --> stack->Particle(id) with id < stack->GetNprimary()

  // if the particle has a daughter primary, we do not want to count it
  if (aParticle->GetFirstDaughter() != -1 && aParticle->GetFirstDaughter() < aTotalPrimaries)
  {
    if (adebug)
      printf("Dropping particle because it has a daughter among the primaries.\n");
    return kFALSE;
  }

  Int_t pdgCode = TMath::Abs(aParticle->GetPdgCode());
  

  // skip quarks and gluon
  if (pdgCode <= 10 || pdgCode == 21)
  {
    if (adebug)
      printf("Dropping particle because it is a quark or gluon.\n");
    return kFALSE;
  }

  Int_t status = aParticle->GetStatusCode();
  // skip non final state particles..
  if(status!=1){
    if (adebug)
      printf("Dropping particle because it is not a final state particle.\n");
    return kFALSE;
  }

  if (strcmp(aParticle->GetName(),"XXX") == 0)
  {
    Printf("WARNING: There is a particle named XXX (pdg code %d).", pdgCode);
    return kFALSE;
  }

  TParticlePDG* pdgPart = aParticle->GetPDG();

  if (strcmp(pdgPart->ParticleClass(),"Unknown") == 0)
  {
    Printf("WARNING: There is a particle with an unknown particle class (pdg code %d).", pdgCode);
    return kFALSE;
  }

  if (pdgPart->Charge() == 0)
  {
    if (adebug)
      printf("Dropping particle because it is not charged.\n");
    return kFALSE;
  }

  return kTRUE;
}

//____________________________________________________________________
void AlidNdPtHelper::CreateProjections(TH3* hist, Bool_t save)
{
  // create projections of 3d hists to all 2d combinations
  // the histograms are not returned, just use them from memory or use this to create them in a file

  TH1* proj = hist->Project3D("yx");
  proj->SetXTitle(hist->GetXaxis()->GetTitle());
  proj->SetYTitle(hist->GetYaxis()->GetTitle());
  if (save)
    proj->Write();

  proj = hist->Project3D("zx");
  proj->SetXTitle(hist->GetXaxis()->GetTitle());
  proj->SetYTitle(hist->GetZaxis()->GetTitle());
  if (save)
    proj->Write();

  proj = hist->Project3D("zy");
  proj->SetXTitle(hist->GetYaxis()->GetTitle());
  proj->SetYTitle(hist->GetZaxis()->GetTitle());
  if (save)
    proj->Write();
}

//____________________________________________________________________
void AlidNdPtHelper::CreateDividedProjections(TH3* hist, TH3* hist2, const char* axis, Bool_t putErrors, Bool_t save)
{
  // create projections of the 3d hists divides them
  // axis decides to which plane, if axis is 0 to all planes
  // the histograms are not returned, just use them from memory or use this to create them in a file

  if (axis == 0)
  {
    CreateDividedProjections(hist, hist2, "yx", putErrors, save);
    CreateDividedProjections(hist, hist2, "zx", putErrors, save);
    CreateDividedProjections(hist, hist2, "zy", putErrors, save);

    return;
  }

  TH1* proj = hist->Project3D(axis);

  if (strlen(axis) == 2)
  {
    proj->SetYTitle(GetAxisTitle(hist, axis[0]));
    proj->SetXTitle(GetAxisTitle(hist, axis[1]));
  }
  else if (strlen(axis) == 1)
    proj->SetXTitle(GetAxisTitle(hist, axis[0]));

  TH1* proj2 = hist2->Project3D(axis);
  if (strlen(axis) == 2)
  {
    proj2->SetYTitle(GetAxisTitle(hist2, axis[0]));
    proj2->SetXTitle(GetAxisTitle(hist2, axis[1]));
  }
  else if (strlen(axis) == 1)
    proj2->SetXTitle(GetAxisTitle(hist2, axis[0]));

  TH1* division = dynamic_cast<TH1*> (proj->Clone(Form("%s_div_%s", proj->GetName(), proj2->GetName())));
  //printf("doing axis: %s, x axis has %d %d bins, min %f %f max %f %f\n", axis, division->GetNbinsX(), proj2->GetNbinsX(), division->GetXaxis()->GetBinLowEdge(1), proj2->GetXaxis()->GetBinLowEdge(1), division->GetXaxis()->GetBinUpEdge(division->GetNbinsX()), proj2->GetXaxis()->GetBinUpEdge(proj2->GetNbinsX()));
  //printf("doing axis: %s, y axis has %d %d bins, min %f %f max %f %f\n", axis, division->GetNbinsY(), proj2->GetNbinsY(), division->GetYaxis()->GetBinLowEdge(1), proj2->GetYaxis()->GetBinLowEdge(1), division->GetYaxis()->GetBinUpEdge(division->GetNbinsY()), proj2->GetYaxis()->GetBinUpEdge(proj2->GetNbinsY()));
  division->Divide(proj, proj2, 1, 1, "B");
  division->SetTitle(Form("%s divided %s", proj->GetTitle(), proj2->GetTitle()));

  if (putErrors)
  {
    division->Sumw2();
    if (division->GetDimension() == 1)
    {
      Int_t nBins = division->GetNbinsX();
      for (Int_t i = 1; i <= nBins; ++i)
        if (proj2->GetBinContent(i) != 0)
          division->SetBinError(i, TMath::Sqrt(proj->GetBinContent(i)) / proj2->GetBinContent(i));
    }
    else if (division->GetDimension() == 2)
    {
      Int_t nBinsX = division->GetNbinsX();
      Int_t nBinsY = division->GetNbinsY();
      for (Int_t i = 1; i <= nBinsX; ++i)
        for (Int_t j = 1; j <= nBinsY; ++j)
          if (proj2->GetBinContent(i, j) != 0)
            division->SetBinError(i, j, TMath::Sqrt(proj->GetBinContent(i, j)) / proj2->GetBinContent(i, j));
    }
  }

  if (save)
  {
    proj->Write();
    proj2->Write();
    division->Write();
  }
}

//____________________________________________________________________
const char* AlidNdPtHelper::GetAxisTitle(TH3* hist, const char axis)
{
  // returns the title of the axis given in axis (x, y, z)

  if (axis == 'x')
    return hist->GetXaxis()->GetTitle();
  else if (axis == 'y')
    return hist->GetYaxis()->GetTitle();
  else if (axis == 'z')
    return hist->GetZaxis()->GetTitle();

  return 0;
}


AlidNdPtHelper::MCProcessType AlidNdPtHelper::GetPythiaEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug) {

  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader);

  if (!pythiaGenHeader) {
    printf("AlidNdPtHelper::GetProcessType : Unknown gen Header type). \n");
    return kInvalidProcess;
  }


  Int_t pythiaType = pythiaGenHeader->ProcessType();
  fgLastProcessType = pythiaType;
  MCProcessType globalType = kInvalidProcess;  


  if (adebug) {
    printf("AlidNdPtHelper::GetProcessType : Pythia process type found: %d \n",pythiaType);
  }


  if(pythiaType==92||pythiaType==93){
    globalType = kSD;
  }
  else if(pythiaType==94){
    globalType = kDD;
  }
  //else if(pythiaType != 91){ // also exclude elastic to be sure... CKB??
  else {
    globalType = kND;
  }
  return globalType;
}


AlidNdPtHelper::MCProcessType AlidNdPtHelper::GetDPMjetEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug) {
  //
  // get the process type of the event.
  //

  // can only read pythia headers, either directly or from cocktalil header
  AliGenDPMjetEventHeader* dpmJetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(aHeader);

  if (!dpmJetGenHeader) {
    printf("AlidNdPtHelper::GetDPMjetProcessType : Unknown header type (not DPMjet or). \n");
    return kInvalidProcess;
  }

  Int_t dpmJetType = dpmJetGenHeader->ProcessType();
  fgLastProcessType = dpmJetType;
  MCProcessType globalType = kInvalidProcess;  


  if (adebug) {
    printf("AlidNdPtHelper::GetDPMJetProcessType : DPMJet process type found: %d \n",dpmJetType);
  }


  if(dpmJetType == 1){ // this is explicitly inelastic
    globalType = kND;
  }  
  else if(dpmJetType==5||dpmJetType==6){
    globalType = kSD;
  }
  else if(dpmJetType==7||dpmJetType==4){// DD and double pomeron
    globalType = kDD;
  }
  return globalType;
}


AlidNdPtHelper::MCProcessType AlidNdPtHelper::GetEventProcessType(AliHeader* aHeader, Bool_t adebug) {
  //
  // get the process type of the event.
  //


  // Check for simple headers first

  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader->GenEventHeader());
  if (pythiaGenHeader) {
    return GetPythiaEventProcessType(pythiaGenHeader,adebug);
  }

  AliGenDPMjetEventHeader* dpmJetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(aHeader->GenEventHeader());
  if (dpmJetGenHeader) {
    return GetDPMjetEventProcessType(dpmJetGenHeader,adebug);
  }
  

  // check for cocktail

  AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(aHeader->GenEventHeader());
  if (!genCocktailHeader) {
    printf("AlidNdPtHelper::GetProcessType : Unknown header type (not Pythia or Cocktail). \n");
    return kInvalidProcess;
  }

  TList* headerList = genCocktailHeader->GetHeaders();
  if (!headerList) {
    return kInvalidProcess;
  }

  for (Int_t i=0; i<headerList->GetEntries(); i++) {

    pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
    if (pythiaGenHeader) {
      return GetPythiaEventProcessType(pythiaGenHeader,adebug);
    }

    dpmJetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(headerList->At(i));
    if (dpmJetGenHeader) {
      return GetDPMjetEventProcessType(dpmJetGenHeader,adebug);
    }
  }
  return kInvalidProcess;
}



//____________________________________________________________________
TParticle* AlidNdPtHelper::FindPrimaryMother(AliStack* stack, Int_t label)
{
  //
  // Finds the first mother among the primary particles of the particle identified by <label>,
  // i.e. the primary that "caused" this particle
  //

  Int_t motherLabel = FindPrimaryMotherLabel(stack, label);
  if (motherLabel < 0)
    return 0;

  return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AlidNdPtHelper::FindPrimaryMotherLabel(AliStack* stack, Int_t label)
{
  //
  // Finds the first mother among the primary particles of the particle identified by <label>,
  // i.e. the primary that "caused" this particle
  //
  // returns its label
  //

  Int_t nPrim  = stack->GetNprimary();

  while (label >= nPrim)
  {
    //printf("Particle %d (pdg %d) is not a primary. Let's check its mother %d\n", label, mother->GetPdgCode(), mother->GetMother(0));

    TParticle* particle = stack->Particle(label);
    if (!particle)
    {
      AliDebugGeneral("FindPrimaryMother", AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack.", label));
      return -1;
    }
 
    // find mother
    if (particle->GetMother(0) < 0)
    {
      AliDebugGeneral("FindPrimaryMother", AliLog::kError, Form("UNEXPECTED: Could not find mother of secondary particle %d.", label));
      return -1;
    }

    label = particle->GetMother(0);
  }

  return label;
}

//____________________________________________________________________
void AlidNdPtHelper::NormalizeToBinWidth(TH1* hist)
{
  //
  // normalizes a 1-d histogram to its bin width
  //

  for (Int_t i=1; i<=hist->GetNbinsX(); ++i)
  {
    hist->SetBinContent(i, hist->GetBinContent(i) / hist->GetBinWidth(i));
    hist->SetBinError(i, hist->GetBinError(i) / hist->GetBinWidth(i));
  }
}

//____________________________________________________________________
void AlidNdPtHelper::NormalizeToBinWidth(TH2* hist)
{
  //
  // normalizes a 2-d histogram to its bin width (x width * y width)
  //

  for (Int_t i=1; i<=hist->GetNbinsX(); ++i)
    for (Int_t j=1; j<=hist->GetNbinsY(); ++j)
    {
      Double_t factor = hist->GetXaxis()->GetBinWidth(i) * hist->GetYaxis()->GetBinWidth(j);
      hist->SetBinContent(i, j, hist->GetBinContent(i, j) / factor);
      hist->SetBinError(i, j, hist->GetBinError(i, j) / factor);
    }
}

//____________________________________________________________________
void AlidNdPtHelper::PrintConf(AnalysisMode analysisMode, Trigger trigger)
{
  //
  // Prints the given configuration
  //

  TString str(">>>> Running with ");

  switch (analysisMode)
  {
    case kInvalid: str += "invalid setting"; break;
    case kSPD : str += "SPD-only"; break;
    case kTPC : str += "TPC-only"; break;
    case kTPCITS : str += "Global tracking"; break;
    case kTPCSPDvtx : str += "TPC tracking + SPD event vertex"; break;
    case kMCRec : str += "TPC tracking + Replace rec. with MC values"; break;
    case kMCPion : str += "TPC tracking + only pion MC tracks"; break;
    case kMCKaon : str += "TPC tracking + only kaon MC tracks"; break;
    case kMCProton : str += "TPC tracking + only proton MC tracks"; break;
    case kPlus: str += "TPC tracking + only positive charged tracks"; break;
    case kMinus : str += "TPC tracking + only negative charge tracks"; break;
  }
  str += " and trigger ";

  switch (trigger)
  {
    case kMB1 : str += "MB1"; break;
    case kMB2 : str += "MB2"; break;
    case kSPDFASTOR : str += "SPD FASTOR"; break;
  }

  str += " <<<<";

  Printf("%s", str.Data());
}

//____________________________________________________________________
Int_t AlidNdPtHelper::ConvertPdgToPid(TParticle *particle) {
//
// Convert Abs(pdg) to pid 
// (0 - e, 1 - muons, 2 - pions, 3 - kaons, 4 - protons, 5 -all rest)
//
Int_t pid=-1;

  if (TMath::Abs(particle->GetPdgCode()) == kElectron)         { pid = 0; }
  else if (TMath::Abs(particle->GetPdgCode()) == kMuonMinus) { pid = 1; }
  else if (TMath::Abs(particle->GetPdgCode()) == kPiPlus)    { pid = 2; }
  else if (TMath::Abs(particle->GetPdgCode()) == kKPlus)     { pid = 3; }
  else if (TMath::Abs(particle->GetPdgCode()) == kProton)    { pid = 4; }
  else                                                       { pid = 5; }

return pid;
}

//_____________________________________________________________________________
TH1F* AlidNdPtHelper::CreateResHisto(TH2F* hRes2, TH1F **phMean, Int_t integ,  Bool_t drawBinFits, Int_t minHistEntries)
{
  TVirtualPad* currentPad = gPad;
  TAxis* axis = hRes2->GetXaxis();
  Int_t nBins = axis->GetNbins();
  Bool_t overflowBinFits = kFALSE;
  TH1F* hRes, *hMean;
  if (axis->GetXbins()->GetSize()){
    hRes = new TH1F("hRes", "", nBins, axis->GetXbins()->GetArray());
    hMean = new TH1F("hMean", "", nBins, axis->GetXbins()->GetArray());
  }
  else{
    hRes = new TH1F("hRes", "", nBins, axis->GetXmin(), axis->GetXmax());
    hMean = new TH1F("hMean", "", nBins, axis->GetXmin(), axis->GetXmax());

  }
  hRes->SetStats(false);
  hRes->SetOption("E");
  hRes->SetMinimum(0.);
  //
  hMean->SetStats(false);
  hMean->SetOption("E");
 
  // create the fit function
  TF1 * fitFunc = new TF1("G","[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",-3,3);
  
  fitFunc->SetLineWidth(2);
  fitFunc->SetFillStyle(0);
  // create canvas for fits
  TCanvas* canBinFits = NULL;
  Int_t nPads = (overflowBinFits) ? nBins+2 : nBins;
  Int_t nx = Int_t(sqrt(nPads-1.));// + 1;
  Int_t ny = (nPads-1) / nx + 1;
  if (drawBinFits) {
    canBinFits = (TCanvas*)gROOT->FindObject("canBinFits");
    if (canBinFits) delete canBinFits;
    canBinFits = new TCanvas("canBinFits", "fits of bins", 200, 100, 500, 700);
    canBinFits->Divide(nx, ny);
  }

  // loop over x bins and fit projection
  Int_t dBin = ((overflowBinFits) ? 1 : 0);
  for (Int_t bin = 1-dBin; bin <= nBins+dBin; bin++) {
    if (drawBinFits) canBinFits->cd(bin + dBin);
    Int_t bin0=TMath::Max(bin-integ,0);
    Int_t bin1=TMath::Min(bin+integ,nBins);
    TH1D* hBin = hRes2->ProjectionY("hBin", bin0, bin1);
    //    
    if (hBin->GetEntries() > minHistEntries) {
      fitFunc->SetParameters(hBin->GetMaximum(),hBin->GetMean(),hBin->GetRMS());
      hBin->Fit(fitFunc,"s");
      Double_t sigma = TMath::Abs(fitFunc->GetParameter(2));

      if (sigma > 0.){
	hRes->SetBinContent(bin, TMath::Abs(fitFunc->GetParameter(2)));
	hMean->SetBinContent(bin, fitFunc->GetParameter(1));	
      }
      else{
	hRes->SetBinContent(bin, 0.);
	hMean->SetBinContent(bin,0);
      }
      hRes->SetBinError(bin, fitFunc->GetParError(2));
      hMean->SetBinError(bin, fitFunc->GetParError(1));
      
      //
      //

    } else {
      hRes->SetBinContent(bin, 0.);
      hRes->SetBinError(bin, 0.);
      hMean->SetBinContent(bin, 0.);
      hMean->SetBinError(bin, 0.);
    }
    

    if (drawBinFits) {
      char name[256];
      if (bin == 0) {
	sprintf(name, "%s < %.4g", axis->GetTitle(), axis->GetBinUpEdge(bin));
      } else if (bin == nBins+1) {
	sprintf(name, "%.4g < %s", axis->GetBinLowEdge(bin), axis->GetTitle());
      } else {
	sprintf(name, "%.4g < %s < %.4g", axis->GetBinLowEdge(bin),
		axis->GetTitle(), axis->GetBinUpEdge(bin));
      }
      canBinFits->cd(bin + dBin);
      hBin->SetTitle(name);
      hBin->SetStats(kTRUE);
      hBin->DrawCopy("E");
      canBinFits->Update();
      canBinFits->Modified();
      canBinFits->Update();
    }
    
    delete hBin;
  }

  delete fitFunc;
  currentPad->cd();
  *phMean = hMean;
  return hRes;
}

//_____________________________________________________________________________
TH1F* AlidNdPtHelper::MakeResol(TH2F * his, Int_t integ, Bool_t type, Bool_t drawBins, Int_t minHistEntries){
// Create resolution histograms
  
     TH1F *hisr=0, *hism=0;
     if (!gPad) new TCanvas;
         //hisr = AliTreeDraw::CreateResHistoI(his,&hism,integ);
         hisr = CreateResHisto(his,&hism,integ,drawBins,minHistEntries);
         if (type) return hism;
         else return hisr;

return hisr;	 
}

//_____________________________________________________________________________
AliESDtrack* AlidNdPtHelper::GetTPCOnlyTrack(AliESDEvent* esd, const AliESDVertex *vtx, Int_t iTrack)
{
  // creates a TPC only track from the given esd track
  // the track has to be deleted by the user
  //
  // NB. most of the functionality to get a TPC only track from an ESD track is in AliESDtrack, where it should be
  // there are only missing propagations here that are needed for old data
  // this function will therefore become obsolete
  //
  // adapted from code provided by CKB

  // no vertex
  if (!vtx) return 0;  
  if(!vtx->GetStatus()) return 0; 

  AliESDtrack* track = esd->GetTrack(iTrack);
  if (!track)
    return 0;

  AliESDtrack *tpcTrack = new AliESDtrack();

  // This should have been done during the reconstruction
  // fixed by Juri in r26675
  // but recalculate for older data CKB
  Float_t p[2],cov[3];
  track->GetImpactParametersTPC(p,cov);
  if(p[0]==0&&p[1]==0)
    //track->RelateToVertexTPC(esd->GetPrimaryVertexTPC(),esd->GetMagneticField(),kVeryBig);
    track->RelateToVertexTPC(vtx,esd->GetMagneticField(),kVeryBig);
  // BKC

  // only true if we have a tpc track
  if (!track->FillTPCOnlyTrack(*tpcTrack))
  {
    delete tpcTrack;
    return 0;
  }

  // propagate to Vertex
  // not needed for normal reconstructed ESDs...
  // Double_t pTPC[2],covTPC[3];
  //tpcTrack->PropagateToDCA(esd->GetPrimaryVertexTPC(), esd->GetMagneticField(), 10000,  pTPC, covTPC);
  //tpcTrack->PropagateToDCA(vtx, esd->GetMagneticField(), 10000,  pTPC, covTPC);

  return tpcTrack;
}

//_____________________________________________________________________________
TObjArray* AlidNdPtHelper::GetAllChargedTracks(AliESDEvent *esdEvent, const AliESDVertex *vtx, AnalysisMode analysisMode)
{
  //
  // all charged TPC particles 
  //
  TObjArray *allTracks = new TObjArray();
  if(!allTracks) return allTracks;

  AliESDtrack *track=0;
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) 
  { 
    if(analysisMode == AlidNdPtHelper::kTPC || analysisMode == AlidNdPtHelper::kTPCSPDvtx || analysisMode == AlidNdPtHelper::kMCRec || analysisMode == kMCPion || analysisMode == kMCKaon || analysisMode == kMCProton || analysisMode ==kPlus || analysisMode ==kMinus) { 
      // track must be deleted by the user 
      track = GetTPCOnlyTrack(esdEvent,vtx,iTrack);
      //track = AliESDtrackCuts::GetTPCOnlyTrack(esdEvent,iTrack);
    } else {
      track=esdEvent->GetTrack(iTrack);
    }
    if(!track) continue;

    if(track->Charge()==0) { 
      if(analysisMode == AlidNdPtHelper::kTPC || analysisMode == AlidNdPtHelper::kTPCSPDvtx ||  analysisMode == AlidNdPtHelper::kMCRec || analysisMode == kMCPion || analysisMode == kMCKaon || analysisMode == kMCProton || analysisMode ==kPlus || analysisMode ==kMinus) {
        delete track; continue; 
      } else {
        continue;
      } 
    }

    allTracks->Add(track);
  }
  if(analysisMode == AlidNdPtHelper::kTPC || analysisMode == AlidNdPtHelper::kTPCSPDvtx || analysisMode == AlidNdPtHelper::kMCRec || analysisMode == kMCPion || analysisMode == kMCKaon || analysisMode == kMCProton || analysisMode ==kPlus || analysisMode ==kMinus) allTracks->SetOwner(kTRUE);

return allTracks;
}

//_____________________________________________________________________________
Int_t AlidNdPtHelper::GetTPCMBTrackMult(AliESDEvent *esdEvent, AlidNdPtEventCuts *evtCuts, AlidNdPtAcceptanceCuts *accCuts, AliESDtrackCuts *trackCuts)
{
  //
  // get MB event track multiplicity
  //
  if(!esdEvent) 
  { 
    ::Error("AlidNdPtHelper::GetTPCMBTrackMult()","esd event is NULL");
    return 0;  
  }
 
  if(!evtCuts || !accCuts || !trackCuts) 
  { 
    ::Error("AlidNdPtHelper::GetTPCMBTrackMult()","cuts not available");
    return 0;  
  }

  //
  Double_t pos[3]={evtCuts->GetMeanXv(),evtCuts->GetMeanYv(),evtCuts->GetMeanZv()};
  Double_t err[3]={evtCuts->GetSigmaMeanXv(),evtCuts->GetSigmaMeanYv(),evtCuts->GetSigmaMeanZv()};
  AliESDVertex vtx0(pos,err);

  //
  Float_t maxDCAr = accCuts->GetMaxDCAr();
  Float_t maxDCAz = accCuts->GetMaxDCAz();
  Float_t minTPCClust = trackCuts->GetMinNClusterTPC();
  //
  Int_t ntracks = esdEvent->GetNumberOfTracks();
  Double_t dca[2],cov[3];
  Int_t mult=0;
  for (Int_t i=0;i <ntracks; i++){
    AliESDtrack *t = esdEvent->GetTrack(i);
    if (!t) continue;
    if (t->Charge() == 0) continue;
    if (!t->GetTPCInnerParam()) continue;
    if (t->GetTPCNcls()<minTPCClust) continue;
    //
    AliExternalTrackParam  *tpcTrack  = new AliExternalTrackParam(*(t->GetTPCInnerParam()));
    if (!tpcTrack->PropagateToDCA(&vtx0,esdEvent->GetMagneticField(),100.,dca,cov)) 
    {
      if(tpcTrack) delete tpcTrack; 
      continue;
    }
    //
    if (TMath::Abs(dca[0])>maxDCAr || TMath::Abs(dca[1])>maxDCAz) {
      if(tpcTrack) delete tpcTrack; 
      continue;
    }

    mult++;    

    if(tpcTrack) delete tpcTrack; 
  }

return mult;  
}

//_____________________________________________________________________________
Int_t AlidNdPtHelper::GetTPCMBPrimTrackMult(AliESDEvent *esdEvent, AliStack * stack, AlidNdPtEventCuts *evtCuts, AlidNdPtAcceptanceCuts *accCuts, AliESDtrackCuts *trackCuts)
{
  //
  // get MB primary event track multiplicity
  //
  if(!esdEvent) 
  { 
    ::Error("AlidNdPtHelper::GetTPCMBPrimTrackMult()","esd event is NULL");
    return 0;  
  }

  if(!stack) 
  { 
    ::Error("AlidNdPtHelper::GetTPCMBPrimTrackMult()","esd event is NULL");
    return 0;  
  }
 
  if(!evtCuts || !accCuts || !trackCuts) 
  { 
    ::Error("AlidNdPtHelper::GetTPCMBPrimTrackMult()","cuts not available");
    return 0;  
  }

  //
  Double_t pos[3]={evtCuts->GetMeanXv(),evtCuts->GetMeanYv(),evtCuts->GetMeanZv()};
  Double_t err[3]={evtCuts->GetSigmaMeanXv(),evtCuts->GetSigmaMeanYv(),evtCuts->GetSigmaMeanZv()};
  AliESDVertex vtx0(pos,err);

  //
  Float_t maxDCAr = accCuts->GetMaxDCAr();
  Float_t maxDCAz = accCuts->GetMaxDCAz();
  Float_t minTPCClust = trackCuts->GetMinNClusterTPC();

  //
  Int_t ntracks = esdEvent->GetNumberOfTracks();
  Double_t dca[2],cov[3];
  Int_t mult=0;
  for (Int_t i=0;i <ntracks; i++){
    AliESDtrack *t = esdEvent->GetTrack(i);
    if (!t) continue;
    if (t->Charge() == 0) continue;
    if (!t->GetTPCInnerParam()) continue;
    if (t->GetTPCNcls()<minTPCClust) continue;
    //
    AliExternalTrackParam  *tpcTrack  = new AliExternalTrackParam(*(t->GetTPCInnerParam()));
    if (!tpcTrack->PropagateToDCA(&vtx0,esdEvent->GetMagneticField(),100.,dca,cov)) 
    {
      if(tpcTrack) delete tpcTrack; 
      continue;
    }
    //
    if (TMath::Abs(dca[0])>maxDCAr || TMath::Abs(dca[1])>maxDCAz) {
      if(tpcTrack) delete tpcTrack; 
      continue;
    }

    Int_t label = TMath::Abs(t->GetLabel());
    TParticle *part = stack->Particle(label);
    if(!part) { 
      if(tpcTrack) delete tpcTrack; 
      continue;
    }
    if(!stack->IsPhysicalPrimary(label)) 
    { 
      if(tpcTrack) delete tpcTrack; 
      continue;
    }

    mult++;    

    if(tpcTrack) delete tpcTrack; 
  }

return mult;  
}





//_____________________________________________________________________________
Int_t AlidNdPtHelper::GetMCTrueTrackMult(AliMCEvent *mcEvent, AlidNdPtEventCuts *evtCuts, AlidNdPtAcceptanceCuts *accCuts)
{
  //
  // calculate mc event true track multiplicity
  //
  if(!mcEvent) return 0;

  AliStack* stack = 0;
  Int_t mult = 0;

  // MC particle stack
  stack = mcEvent->Stack();
  if (!stack) return 0;

  Bool_t isEventOK = evtCuts->AcceptMCEvent(mcEvent);
  if(!isEventOK) return 0; 

  Int_t nPart  = stack->GetNtrack();
  for (Int_t iMc = 0; iMc < nPart; ++iMc) 
  {
     TParticle* particle = stack->Particle(iMc);
     if (!particle)
     continue;

     // only charged particles
     Double_t charge = particle->GetPDG()->Charge()/3.;
     if (charge == 0.0)
     continue;
      
     // physical primary
     Bool_t prim = stack->IsPhysicalPrimary(iMc);
     if(!prim) continue;

     // checked accepted
     if(accCuts->AcceptTrack(particle)) 
     {
       mult++;
     }
  }

return mult;  
}

//_______________________________________________________________________
void  AlidNdPtHelper::PrintMCInfo(AliStack *pStack,Int_t label)
{
// print information about particles in the stack

  if(!pStack)return;
  label = TMath::Abs(label);
  TParticle *part = pStack->Particle(label);
  Printf("########################");
  Printf("%s:%d %d UniqueID %d PDG %d P %3.3f",(char*)__FILE__,__LINE__,label,part->GetUniqueID(),part->GetPdgCode(),part->P());
  part->Print();
  TParticle* mother = part;
  Int_t imo = part->GetFirstMother();
  Int_t nprim = pStack->GetNprimary();

  while((imo >= nprim)) {
      mother =  pStack->Particle(imo);
      Printf("Mother %s:%d Label %d UniqueID %d PDG %d P %3.3f",(char*)__FILE__,__LINE__,imo,mother->GetUniqueID(),mother->GetPdgCode(),mother->P());
      mother->Print();
      imo =  mother->GetFirstMother();
 }

 Printf("########################");
}


//_____________________________________________________________________________
TH1* AlidNdPtHelper::GetContCorrHisto(TH1 *hist) 
{
//
// get contamination histogram
//
 if(!hist) return 0;

 Int_t nbins = hist->GetNbinsX();
 TH1 *h_cont = (TH1D *)hist->Clone();

 for(Int_t i=0; i<=nbins+1; i++) {
   Double_t binContent = hist->GetBinContent(i);
   Double_t binError = hist->GetBinError(i);

   h_cont->SetBinContent(i,1.-binContent);
   h_cont->SetBinError(i,binError);
 }

return h_cont;
}


//_____________________________________________________________________________
TH1* AlidNdPtHelper::ScaleByBinWidth(TH1 *hist) 
{
//
// scale by bin width
//
 if(!hist) return 0;

 TH1 *h_scale = (TH1D *)hist->Clone();
 h_scale->Scale(1.,"width");

return h_scale;
}

//_____________________________________________________________________________
TH1* AlidNdPtHelper::CalcRelativeDifference(TH1 *hist1, TH1 *hist2) 
{
//
// calculate rel. difference 
//

 if(!hist1) return 0;
 if(!hist2) return 0;

 TH1 *h1_clone = (TH1D *)hist1->Clone();
 h1_clone->Sumw2();

 // (rec-mc)/mc
 h1_clone->Add(hist2,-1);
 h1_clone->Divide(hist2);

return h1_clone;
}

//_____________________________________________________________________________
TH1* AlidNdPtHelper::CalcRelativeDifferenceFun(TH1 *hist1, TF1 *fun) 
{
//
// calculate rel. difference
// between histogram and function
//
 if(!hist1) return 0;
 if(!fun) return 0;

 TH1 *h1_clone = (TH1D *)hist1->Clone();
 h1_clone->Sumw2();

 // 
 h1_clone->Add(fun,-1);
 h1_clone->Divide(hist1);

return h1_clone;
}

//_____________________________________________________________________________
TH1* AlidNdPtHelper::NormalizeToEvent(TH2 *hist1, TH1 *hist2) 
{
// normalise to event for a given multiplicity bin
// return pt histogram 

 if(!hist1) return 0;
 if(!hist2) return 0;
 char name[256];

 Int_t nbinsX = hist1->GetNbinsX();
 //Int_t nbinsY = hist1->GetNbinsY();

 TH1D *hist_norm = 0;
 for(Int_t i=0; i<=nbinsX+1; i++) {
   sprintf(name,"mom_%d",i);
   TH1D *hist = (TH1D*)hist1->ProjectionY(name,i+1,i+1);

   sprintf(name,"mom_norm");
   if(i==0) { 
     hist_norm = (TH1D *)hist->Clone(name);
     hist_norm->Reset();
   }

   Double_t nbEvents = hist2->GetBinContent(i);
   if(!nbEvents) { nbEvents = 1.; };

   hist->Scale(1./nbEvents);
   hist_norm->Add(hist);
 }

return hist_norm;
}

//_____________________________________________________________________________
THnSparse* AlidNdPtHelper::GenerateCorrMatrix(THnSparse *hist1, THnSparse *hist2, char *name) {
// generate correction matrix
if(!hist1 || !hist2) return 0; 

THnSparse *h =(THnSparse*)hist1->Clone(name);;
h->Divide(hist1,hist2,1,1,"B");

return h;
}

//_____________________________________________________________________________
TH2* AlidNdPtHelper::GenerateCorrMatrix(TH2 *hist1, TH2 *hist2, char *name) {
// generate correction matrix
if(!hist1 || !hist2) return 0; 

TH2D *h =(TH2D*)hist1->Clone(name);;
h->Divide(hist1,hist2,1,1,"B");

return h;
}

//_____________________________________________________________________________
TH1* AlidNdPtHelper::GenerateCorrMatrix(TH1 *hist1, TH1 *hist2, char *name) {
// generate correction matrix
if(!hist1 || !hist2) return 0; 

TH1D *h =(TH1D*)hist1->Clone(name);;
h->Divide(hist1,hist2,1,1,"B");

return h;
}

//_____________________________________________________________________________
THnSparse* AlidNdPtHelper::GenerateContCorrMatrix(THnSparse *hist1, THnSparse *hist2, char *name) {
// generate contamination correction matrix
if(!hist1 || !hist2) return 0; 

THnSparse *hist =  GenerateCorrMatrix(hist1, hist2, name);
if(!hist) return 0;

// only for non ZERO bins!!!!

Int_t* coord = new Int_t[hist->GetNdimensions()];
memset(coord, 0, sizeof(Int_t) * hist->GetNdimensions());

  for (Long64_t i = 0; i < hist->GetNbins(); ++i) {
    Double_t v = hist->GetBinContent(i, coord);
    hist->SetBinContent(coord, 1.0-v);
    //printf("v %f, hist->GetBinContent(i, coord) %f \n",v,hist->GetBinContent(i, coord));
    Double_t err = hist->GetBinError(coord);
    hist->SetBinError(coord, err);
  }

delete [] coord;

return hist;
}

//_____________________________________________________________________________
TH2* AlidNdPtHelper::GenerateContCorrMatrix(TH2 *hist1, TH2 *hist2, char *name) {
// generate contamination correction matrix
if(!hist1 || !hist2) return 0; 

TH2 *hist = GenerateCorrMatrix(hist1, hist2, name);
if(!hist) return 0;

Int_t nBinsX = hist->GetNbinsX();
Int_t nBinsY = hist->GetNbinsY();

  for (Int_t i = 0; i < nBinsX+1; i++) {
  for (Int_t j = 0; j < nBinsY+1; j++) {
     Double_t cont = hist->GetBinContent(i,j);
     hist->SetBinContent(i,j,1.-cont);
     Double_t err = hist->GetBinError(i,j);
     hist->SetBinError(i,j,err);
  }
  }

return hist;
}

//_____________________________________________________________________________
TH1* AlidNdPtHelper::GenerateContCorrMatrix(TH1 *hist1, TH1 *hist2, char *name) {
// generate contamination correction matrix
if(!hist1 || !hist2) return 0; 

TH1 *hist = GenerateCorrMatrix(hist1, hist2, name);
if(!hist) return 0;

Int_t nBinsX = hist->GetNbinsX();

  for (Int_t i = 0; i < nBinsX+1; i++) {
     Double_t cont = hist->GetBinContent(i);
     hist->SetBinContent(i,1.-cont);
     Double_t err = hist->GetBinError(i);
     hist->SetBinError(i,err);
  }

return hist;
}

//_____________________________________________________________________________
const AliESDVertex* AlidNdPtHelper::GetTPCVertexZ(AliESDEvent* esdEvent, Float_t sigmaXYcut, Float_t distXYcut, Float_t distZcut, Int_t nclCut, Float_t fraction, Int_t ntracksMin){
  //
  // TPC Z vertexer
  //
  Double_t vtxpos[3]={0.,0.,0.};
  Double_t vtxsigma[3]={.2,.2,100.};
  AliESDVertex vtx0(vtxpos,vtxsigma);
  //
  Int_t ntracks = esdEvent->GetNumberOfTracks();
  TVectorD ztrack(ntracks);
  //Float_t dcar, dcaz;
  //Float_t point[2],cov[3];
  Double_t dca[2],cov[3];
  Int_t counter=0;
  for (Int_t i=0;i <ntracks; i++){
    AliESDtrack *t = esdEvent->GetTrack(i);
    if (!t) continue;
    if (!t->GetTPCInnerParam()) continue;
    if (t->GetTPCNcls()<nclCut) continue;
    //
    AliExternalTrackParam  *tpcTrack  = new AliExternalTrackParam(*(t->GetTPCInnerParam()));
    if (!tpcTrack->PropagateToDCA(&vtx0,esdEvent->GetMagneticField(),100.,dca,cov)) continue;

    //
    if (TMath::Abs(dca[0])>distXYcut) continue;
    if (TMath::Sqrt(cov[0])>sigmaXYcut) continue;    
    if (TMath::Abs(tpcTrack->GetZ())>distZcut) continue;

    /*
    t->GetImpactParametersTPC(dcar,dcaz);
    if (TMath::Abs(dcar)>distXYcut) continue;
    //
    t->GetImpactParametersTPC(point,cov);
    if (TMath::Sqrt(cov[0])>sigmaXYcut) continue;    
    //
    AliExternalTrackParam  tpcTrack(*(t->GetTPCInnerParam()));
    if (!tpcTrack.PropagateToDCA(&vtx0,esdEvent->GetMagneticField(), 100)) continue;
    if (TMath::Abs(tpcTrack.GetZ())>distZcut) continue;
    */
    ztrack[counter]=tpcTrack->GetZ();
    counter++;    

    if(tpcTrack) delete tpcTrack;
  }
  //
  // Find LTM z position
  //
  Double_t mean=0, sigma=0;
  if (counter<ntracksMin) return 0;
  //
  Int_t nused = TMath::Nint(counter*fraction);
  if (nused==counter) nused=counter-1;  
  if (nused>1){
    AliMathBase::EvaluateUni(counter, ztrack.GetMatrixArray(), mean,sigma, TMath::Nint(counter*fraction));
    sigma/=TMath::Sqrt(nused);
  }else{
    mean  = TMath::Mean(counter, ztrack.GetMatrixArray());
    sigma = TMath::RMS(counter, ztrack.GetMatrixArray());
    sigma/=TMath::Sqrt(counter-1);
  }
  vtxpos[2]=mean;
  vtxsigma[2]=sigma;
  const AliESDVertex* vertex = new AliESDVertex(vtxpos, vtxsigma);
  return vertex;
}

//_____________________________________________________________________________
Int_t  AlidNdPtHelper::GetSPDMBTrackMult(AliESDEvent* esdEvent, Float_t deltaThetaCut, Float_t deltaPhiCut) 
{
  //
  // SPD track multiplicity
  //

  // get tracklets
  const AliMultiplicity* mult = esdEvent->GetMultiplicity();
  if (!mult)
     return 0;

  // get multiplicity from SPD tracklets
  Int_t inputCount = 0; 
  for (Int_t i=0; i<mult->GetNumberOfTracklets(); ++i)
  {
    //printf("%d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), mult->GetDeltaPhi(i));

     Float_t phi = mult->GetPhi(i);
     if (phi < 0)
       phi += TMath::Pi() * 2;
     Float_t deltaPhi = mult->GetDeltaPhi(i);
     Float_t deltaTheta = mult->GetDeltaTheta(i);

     if (TMath::Abs(deltaPhi) > 1)
       printf("WARNING: Very high Delta Phi: %d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), deltaPhi);

     if (deltaThetaCut > 0. && TMath::Abs(deltaTheta) > deltaThetaCut)
        continue;

     if (deltaPhiCut > 0. && TMath::Abs(deltaPhi) > deltaPhiCut)
        continue;
      
     ++inputCount;
  }

return inputCount;
}

//_____________________________________________________________________________
Int_t  AlidNdPtHelper::GetSPDMBPrimTrackMult(AliESDEvent* esdEvent, AliStack* stack, Float_t deltaThetaCut, Float_t deltaPhiCut) 
{
  //
  // SPD track multiplicity
  //

  // get tracklets
  const AliMultiplicity* mult = esdEvent->GetMultiplicity();
  if (!mult)
     return 0;

  // get multiplicity from SPD tracklets
  Int_t inputCount = 0; 
  for (Int_t i=0; i<mult->GetNumberOfTracklets(); ++i)
  {
    //printf("%d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), mult->GetDeltaPhi(i));

     Float_t phi = mult->GetPhi(i);
     if (phi < 0)
       phi += TMath::Pi() * 2;
     Float_t deltaPhi = mult->GetDeltaPhi(i);
     Float_t deltaTheta = mult->GetDeltaTheta(i);

     if (TMath::Abs(deltaPhi) > 1)
       printf("WARNING: Very high Delta Phi: %d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), deltaPhi);

     if (deltaThetaCut > 0. && TMath::Abs(deltaTheta) > deltaThetaCut)
        continue;

     if (deltaPhiCut > 0. && TMath::Abs(deltaPhi) > deltaPhiCut)
        continue;


     if (mult->GetLabel(i, 0) < 0 || mult->GetLabel(i, 0) != mult->GetLabel(i, 1) || 
         !stack->IsPhysicalPrimary(mult->GetLabel(i, 0)))
        continue;

      
     ++inputCount;
  }

return inputCount;
}


