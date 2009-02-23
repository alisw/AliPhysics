/* $Id$ */

#include <AliPWG0Helper.h>

#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TList.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>

#include <AliHeader.h>
#include <AliStack.h>
#include <AliLog.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliVertexerTracks.h>
#include <AliMultiplicity.h>

#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenCocktailEventHeader.h>
#include <AliGenDPMjetEventHeader.h>
#include <AliESDVZERO.h>

//____________________________________________________________________
ClassImp(AliPWG0Helper)

Int_t AliPWG0Helper::fgLastProcessType = -1;

//____________________________________________________________________
Bool_t AliPWG0Helper::IsEventTriggered(const AliESDEvent* aEsd, Trigger trigger)
{
  // checks if an event has been triggered
  // this function implements the "offline" methods that use the ESD, other trigger requests are passed to the function prototype with ULong_t

  Int_t firedChips = 0;
  Bool_t v0A = kFALSE;
  Bool_t v0C = kFALSE;

  // offline triggers have to be dealt with here, because we need the esd pointer
  if (trigger == kOfflineFASTOR || trigger == kOfflineMB1 || trigger == kOfflineMB2 || trigger == kOfflineMB3)
  {
      const AliMultiplicity* mult = aEsd->GetMultiplicity();
      if (!mult)
      {
        Printf("AliPWG0Helper::IsEventTriggered: ERROR: AliMultiplicity not available");
        return kFALSE;
      }
      firedChips = mult->GetNumberOfFiredChips(0) + mult->GetNumberOfFiredChips(1);
  }
  if (trigger == kOfflineMB1 || trigger == kOfflineMB2 || trigger == kOfflineMB3)
  {
    AliESDVZERO* v0Data = aEsd->GetVZEROData();
    if (!v0Data)
    {
      Printf("AliPWG0Helper::IsEventTriggered: ERROR: AliESDVZERO not available");
      return kFALSE;
    }
    for (Int_t i=0; i<32; i++)
    {
      if (v0Data->BBTriggerV0A(i))
        v0A = kTRUE;
      if (v0Data->BBTriggerV0C(i))
        v0C = kTRUE;
    }
  }
      
  switch (trigger)
  {
    case kOfflineFASTOR:
    {
      if (firedChips > 0)
        return kTRUE;
      break;
    }
    case kOfflineMB1:
    {
      if ((firedChips > 0) || v0A || v0C)
        return kTRUE;
      break;
    }
    case kOfflineMB2:
    {
      if ((firedChips > 0) && (v0A || v0C))
        return kTRUE;
      break;
    }
    case kOfflineMB3:
    {
      if ((firedChips > 0) && v0A && v0C)
        return kTRUE;
      break;
    }
    default:
    {
      return IsEventTriggered(aEsd->GetTriggerMask(), trigger);
      break;
    }
  }
  
  return kFALSE;
}

//____________________________________________________________________
Bool_t AliPWG0Helper::IsEventTriggered(ULong64_t triggerMask, Trigger trigger)
{
  // check if the event was triggered
  //
  // this function needs the branch fTriggerMask
  
  // definitions from p-p.cfg
  ULong64_t spdFO = (1 << 14);
  ULong64_t v0left = (1 << 10);
  ULong64_t v0right = (1 << 11);

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
    case kMB3:
    {
      if (triggerMask & spdFO && (triggerMask & v0left) && (triggerMask & v0right))
        return kTRUE;
      break;
    }
    case kSPDFASTOR:
    {
      if (triggerMask & spdFO)
        return kTRUE;
      break;
    }
    default:
      Printf("IsEventTriggered: ERROR: Trigger type %d not implemented in this method", (Int_t) trigger);
      break;
  }

  return kFALSE;
}

//____________________________________________________________________
Bool_t AliPWG0Helper::TestVertex(const AliESDVertex* vertex, AnalysisMode analysisMode, Bool_t debug)
{
    // Checks if a vertex meets the needed quality criteria

  Float_t requiredZResolution = -1;
  if (analysisMode == kSPD || analysisMode == kTPCITS)
  {
    requiredZResolution = 0.1;
  }
  else if (analysisMode == kTPC)
    requiredZResolution = 10.;

  // check resolution
  Double_t zRes = vertex->GetZRes();

  if (zRes > requiredZResolution) {
    if (debug)
      Printf("AliPWG0Helper::TestVertex: Resolution too poor %f (required: %f", zRes, requiredZResolution);
    return kFALSE;
  }

  return kTRUE;
}

//____________________________________________________________________
const AliESDVertex* AliPWG0Helper::GetVertex(AliESDEvent* aEsd, AnalysisMode analysisMode, Bool_t debug, Bool_t bRedoTPC)
{
  // Get the vertex from the ESD and returns it if the vertex is valid
  //
  // Second argument decides which vertex is used (this selects
  // also the quality criteria that are applied)

  const AliESDVertex* vertex = 0;
  if (analysisMode == kSPD || analysisMode == kTPCITS)
  {
    vertex = aEsd->GetPrimaryVertexSPD();
    if (debug)
      Printf("AliPWG0Helper::GetVertex: Returning SPD vertex");
  }
  else if (analysisMode == kTPC)
  {
    if(bRedoTPC){
      if (debug)
        Printf("AliPWG0Helper::GetVertex: Redoing vertex");
      Double_t kBz = aEsd->GetMagneticField();
      AliVertexerTracks vertexer(kBz);
      vertexer.SetTPCMode();
      AliESDVertex *vTPC = vertexer.FindPrimaryVertex(aEsd);
      aEsd->SetPrimaryVertexTPC(vTPC);
      for (Int_t i=0; i<aEsd->GetNumberOfTracks(); i++) {
	AliESDtrack *t = aEsd->GetTrack(i);
	t->RelateToVertexTPC(vTPC, kBz, kVeryBig);
      }
      delete vTPC;
    }

    vertex = aEsd->GetPrimaryVertexTPC();
    if (debug)
      Printf("AliPWG0Helper::GetVertex: Returning vertex from tracks");
  }
  else
    Printf("AliPWG0Helper::GetVertex: ERROR: Invalid second argument %d", analysisMode);

  if (!vertex) {
    if (debug)
      Printf("AliPWG0Helper::GetVertex: No vertex found in ESD");
    return 0;
  }

  // check Ncontributors
  if (vertex->GetNContributors() <= 0) {
    if (debug){
      Printf("AliPWG0Helper::GetVertex: NContributors() <= 0: %d",vertex->GetNContributors());
      Printf("AliPWG0Helper::GetVertex: NIndices(): %d",vertex->GetNIndices());
      vertex->Print();
    }
    return 0;
  }

  // check resolution
  Double_t zRes = vertex->GetZRes();
  if (zRes == 0) {
    Printf("AliPWG0Helper::GetVertex: UNEXPECTED: resolution is 0.");
    return 0;
  }

  if (debug)
  {
    Printf("AliPWG0Helper::GetVertex: Returning valid vertex: %s", vertex->GetTitle());
    vertex->Print();
  }

  return vertex;
}

//____________________________________________________________________
Bool_t AliPWG0Helper::IsPrimaryCharged(TParticle* aParticle, Int_t aTotalPrimaries, Bool_t adebug)
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
void AliPWG0Helper::CreateProjections(TH3* hist, Bool_t save)
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
void AliPWG0Helper::CreateDividedProjections(TH3* hist, TH3* hist2, const char* axis, Bool_t putErrors, Bool_t save)
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
const char* AliPWG0Helper::GetAxisTitle(TH3* hist, const char axis)
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


AliPWG0Helper::MCProcessType AliPWG0Helper::GetPythiaEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug) {

  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader);

  if (!pythiaGenHeader) {
    printf("AliPWG0Helper::GetProcessType : Unknown gen Header type). \n");
    return kInvalidProcess;
  }


  Int_t pythiaType = pythiaGenHeader->ProcessType();
  fgLastProcessType = pythiaType;
  MCProcessType globalType = kInvalidProcess;  


  if (adebug) {
    printf("AliPWG0Helper::GetProcessType : Pythia process type found: %d \n",pythiaType);
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


AliPWG0Helper::MCProcessType AliPWG0Helper::GetDPMjetEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug) {
  //
  // get the process type of the event.
  //

  // can only read pythia headers, either directly or from cocktalil header
  AliGenDPMjetEventHeader* dpmJetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(aHeader);

  if (!dpmJetGenHeader) {
    printf("AliPWG0Helper::GetDPMjetProcessType : Unknown header type (not DPMjet or). \n");
    return kInvalidProcess;
  }

  Int_t dpmJetType = dpmJetGenHeader->ProcessType();
  fgLastProcessType = dpmJetType;
  MCProcessType globalType = kInvalidProcess;  


  if (adebug) {
    printf("AliPWG0Helper::GetDPMJetProcessType : DPMJet process type found: %d \n",dpmJetType);
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


AliPWG0Helper::MCProcessType AliPWG0Helper::GetEventProcessType(AliHeader* aHeader, Bool_t adebug) {
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
    printf("AliPWG0Helper::GetProcessType : Unknown header type (not Pythia or Cocktail). \n");
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
TParticle* AliPWG0Helper::FindPrimaryMother(AliStack* stack, Int_t label)
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
Int_t AliPWG0Helper::FindPrimaryMotherLabel(AliStack* stack, Int_t label)
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
void AliPWG0Helper::NormalizeToBinWidth(TH1* hist)
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
void AliPWG0Helper::NormalizeToBinWidth(TH2* hist)
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
void AliPWG0Helper::PrintConf(AnalysisMode analysisMode, Trigger trigger)
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
  }

  str += " and trigger ";

  switch (trigger)
  {
    case kMB1 : str += "MB1"; break;
    case kMB2 : str += "MB2"; break;
    case kMB3 : str += "MB3"; break;
    case kSPDFASTOR : str += "SPD FASTOR"; break;
    case kOfflineMB1 : str += "Offline MB1"; break;
    case kOfflineMB2 : str += "Offline MB2"; break;
    case kOfflineMB3 : str += "Offline MB3"; break;
    case kOfflineFASTOR : str += "Offline SPD FASTOR"; break;
  }

  str += " <<<<";

  Printf("%s", str.Data());
}

