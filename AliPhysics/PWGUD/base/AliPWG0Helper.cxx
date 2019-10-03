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
#include <AliESD.h>
#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliESDRun.h>
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
Bool_t AliPWG0Helper::TestVertex(const AliESDVertex* vertex, AnalysisMode analysisMode, Bool_t debug)
{
    // Checks if a vertex meets the needed quality criteria

  Float_t requiredZResolution = -1;
  if (analysisMode & kSPD || analysisMode & kTPCITS || analysisMode & kTPCSPD)
  {
    // disable cut on resolution
    requiredZResolution = 1000;
  }
  else if (analysisMode & kTPC)
    requiredZResolution = 10.;

  // check resolution
  Double_t zRes = vertex->GetZRes();

  if (zRes > requiredZResolution) {
    if (debug)
      Printf("AliPWG0Helper::TestVertex: Resolution too poor %f (required: %f", zRes, requiredZResolution);
    return kFALSE;
  }
  
  if (vertex->IsFromVertexerZ())
  {
    if (vertex->GetDispersion() > 0.02) 
    {
      if (debug)
        Printf("AliPWG0Helper::TestVertex: Delta Phi too large in Vertexer Z: %f (required: %f", vertex->GetDispersion(), 0.02);
      return kFALSE;
    }
  }

  return kTRUE;
}

//____________________________________________________________________
const AliESDVertex* AliPWG0Helper::GetVertex(const AliESDEvent* aEsd, AnalysisMode analysisMode, Bool_t debug)
{
  // Get the vertex from the ESD and returns it if the vertex is valid
  //
  // Second argument decides which vertex is used (this selects
  // also the quality criteria that are applied)

  const AliESDVertex* vertex = 0;
  if (analysisMode & kSPD)
  {
    vertex = aEsd->GetPrimaryVertexSPD();
    if (debug)
      Printf("AliPWG0Helper::GetVertex: Returning SPD vertex");
  }
  else if (analysisMode & kTPCITS || analysisMode & kTPCSPD)
  {
    vertex = aEsd->GetPrimaryVertexTracks();
    if (debug)
      Printf("AliPWG0Helper::GetVertex: Returning vertex from tracks");
    if (!vertex || vertex->GetNContributors() <= 0)
    {
      if (debug)
        Printf("AliPWG0Helper::GetVertex: Vertex from tracks has no contributors. Falling back to SPD vertex.");
      vertex = aEsd->GetPrimaryVertexSPD();
    }
  }
  else if (analysisMode & kTPC)
  {
    vertex = aEsd->GetPrimaryVertexTPC();
    if (debug)
      Printf("AliPWG0Helper::GetVertex: Returning vertex from TPC-only tracks");
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

  TH1* division = static_cast<TH1*> (proj->Clone(Form("%s_div_%s", proj->GetName(), proj2->GetName())));
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


  if (dpmJetType == 1 || dpmJetType == 4) { // explicitly inelastic plus central diffraction
    globalType = kND;
  }  
  else if (dpmJetType==5 || dpmJetType==6) {
    globalType = kSD;
  }
  else if (dpmJetType==7) {
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
void AliPWG0Helper::PrintConf(AnalysisMode analysisMode, AliTriggerAnalysis::Trigger trigger, AliPWG0Helper::DiffTreatment diffTreatment)
{
  //
  // Prints the given configuration
  //

  TString str(">>>> Running with >");

  if (analysisMode & kSPD)
    str += "SPD-only";
    
  if (analysisMode & kSPDOnlyL0)
    str += " (only L0 clusters)";
  
  if (analysisMode & kTPC)
     str += "TPC-only";
    
  if (analysisMode & kTPCITS)
     str += "Global tracking";
  
  if (analysisMode & kTPCSPD) 
    str += "Tracks and tracklets";

  if (analysisMode & kFieldOn)
  {
     str += " (with field)";
  }
  else
     str += " (WITHOUT field)";
  
  str += "< and trigger >";
  str += AliTriggerAnalysis::GetTriggerName(trigger);
  str += "< and diffractive treatment >";
  
  switch (diffTreatment)
  {
    case kMCFlags:
      str += "MC flags";
      break;
    
    case kUA5Cuts:
      str += "UA5 cuts";
      break;
       
    case kE710Cuts:
      str += "E710 cuts";
      break;
    
    case kALICEHadronLevel:
      str += "ALICE hadron level";
      break;
  }
  
  str += "< <<<<";

  Printf("%s", str.Data());
}

//____________________________________________________________________
Double_t AliPWG0Helper::Rapidity(Double_t pt, Double_t pz, Double_t m)
{
  //
  // calculates rapidity keeping the sign in case E == pz
  //

  Double_t energy = TMath::Sqrt(pt*pt+pz*pz+m*m);
  if (energy != TMath::Abs(pz))
    return 0.5*TMath::Log((energy+pz)/(energy-pz));

  Printf("W- mt=0");
  return TMath::Sign(1.e30,pz);
}

//____________________________________________________________________
Bool_t AliPWG0Helper::IsHadronLevelSingleDiffractive(AliStack* stack, Float_t cms, Float_t xiMin, Float_t xiMax)
{
  //
  // return if process is single diffractive on hadron level
  // 
  // xiMax and xiMin cut on M^2/s
  //
  // Based on code from Martin Poghoysan
  //
  
  TParticle* part1 = 0;
  TParticle* part2 = 0;
  
  Double_t smallestY = 1e10;
  Double_t largestY  = -1e10;
  
  for (Int_t iParticle = 0; iParticle < stack->GetNprimary(); iParticle++)
  {
    TParticle* part = stack->Particle(iParticle);
    if (!part)
      continue;

    Int_t pdg = TMath::Abs(part->GetPdgCode());

    Int_t child1 = part->GetFirstDaughter();
    Int_t ist    = part->GetStatusCode();

    Int_t mfl  = Int_t (pdg / TMath::Power(10, Int_t(TMath::Log10(pdg))));
    if (child1 > -1 || ist != 1)
      mfl = 0; // select final state charm and beauty
    if (!(stack->IsPhysicalPrimary(iParticle) || pdg == 111 || pdg == 3212 || pdg==3124 || mfl >= 4)) 
      continue;
    Int_t imother = part->GetFirstMother();
    if (imother>0)
    {
      TParticle *partM = stack->Particle(imother);
      Int_t pdgM=TMath::Abs(partM->GetPdgCode());
      if (pdgM==111 || pdgM==3124 || pdgM==3212) 
        continue;
    }
    
    Double_t y = 0;

    // fix for problem with getting mass of particle 3124
    if (pdg != 3124)
      y = Rapidity(part->Pt(), part->Pz(), part->GetMass());
    else 
      y = Rapidity(part->Pt(), part->Pz(), 1.5195);
      
    if (y < smallestY)
    {
      smallestY = y;
      part1 = part;
    }
    
    if (y > largestY)
    {
      largestY = y;
      part2 = part;
    }
  }
  
  if (part1 == 0 || part2 == 0)
    return kFALSE;

  Int_t pdg1 = part1->GetPdgCode();
  Int_t pdg2 = part2->GetPdgCode();

  Double_t pt1 = part1->Pt();
  Double_t pt2 = part2->Pt();
  Double_t pz1 = part1->Pz();
  Double_t pz2 = part2->Pz();
  
  Double_t y1 = TMath::Abs(Rapidity(pt1, pz1, 0.938));
  Double_t y2 = TMath::Abs(Rapidity(pt2, pz2, 0.938));
  
  Int_t arm = -99999;
  if (pdg1 == 2212 && pdg2 == 2212)
  {
    if (y1 > y2) 
      arm = 0;
    else
      arm = 1;
  }
  else if (pdg1 == 2212) 
    arm = 0;
  else if (pdg2 == 2212) 
    arm = 1;

  Double_t M02s = 1. - 2 * part1->Energy() / cms;
  Double_t M12s = 1. - 2 * part2->Energy() / cms;

  if (arm == 0 && M02s > xiMin && M02s < xiMax)
    return kTRUE;
  else if (arm == 1 && M12s > xiMin && M12s < xiMax)
    return kTRUE;

  return kFALSE;
}

//____________________________________________________________________
AliPWG0Helper::MCProcessType AliPWG0Helper::GetEventProcessType(AliESDEvent* esd, AliHeader* header, AliStack* stack, AliPWG0Helper::DiffTreatment diffTreatment)
{
  // 
  // get process type
  //
  // diffTreatment:
  //   kMCFlags: use MC flags
  //   kUA5Cuts: SD events are those that have the MC SD flag and fulfill M^2/s < 0.05; DD from MC flag; Remainder is ND
  //   kE710Cuts: SD events are those that have the MC SD flag and fulfill 2 < M^2 < 0.05s; DD from MC flag; Remainder is ND
  //   kALICEHadronLevel: SD events are those that fulfill M^2/s < 0.05; DD from MC flag; Remainder is ND
  //

  MCProcessType mcProcessType = GetEventProcessType(header);
  
  if (diffTreatment == kMCFlags)
    return mcProcessType;
    
  if (!esd)
  {
    Printf("ERROR: AliPWG0Helper::GetEventProcessType: diffTreatment != kMCFlags and esd == 0");
    return kInvalidProcess;
  }
    
  Float_t cms = esd->GetESDRun()->GetBeamEnergy();
  if (esd->GetESDRun()->IsBeamEnergyIsSqrtSHalfGeV())
    cms *= 2;
  //Printf("cms = %f", cms);

  if (diffTreatment == kUA5Cuts && mcProcessType == kSD)
  {
    if (IsHadronLevelSingleDiffractive(stack, cms, 0, 0.05))
      return kSD;
  }
  else if (diffTreatment == kE710Cuts && mcProcessType == kSD)
  {
    if (IsHadronLevelSingleDiffractive(stack, cms, 2. / cms / cms, 0.05))
      return kSD;
  }
  else if (diffTreatment == kALICEHadronLevel)
  {
    if (IsHadronLevelSingleDiffractive(stack, cms, 0, 0.05))
      return kSD;
  }
  
  if (mcProcessType == kSD)
    return kND;
    
  return mcProcessType;
}
