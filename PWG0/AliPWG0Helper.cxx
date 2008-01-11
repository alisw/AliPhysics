/* $Id$ */

#include <AliPWG0Helper.h>

#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1.h>
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
#include <AliESDVertex.h>

#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenCocktailEventHeader.h>

//____________________________________________________________________
ClassImp(AliPWG0Helper)

//____________________________________________________________________
Bool_t AliPWG0Helper::IsEventTriggered(const AliESD* aEsd, Trigger trigger)
{
  // see function with ULong64_t argument

  ULong64_t triggerMask = aEsd->GetTriggerMask();
  return IsEventTriggered(triggerMask, trigger);
}

//____________________________________________________________________
Bool_t AliPWG0Helper::IsEventTriggered(ULong64_t triggerMask, Trigger trigger)
{
  // check if the event was triggered
  //
  // this function needs the branch fTriggerMask
  //
  // MB should be
  // ITS_SPD_GFO_L0  : 32
  // VZERO_OR_LEFT   : 1
  // VZERO_OR_RIGHT  : 2

  switch (trigger)
  {
    case kMB1:
    {
      if (triggerMask&32 || ((triggerMask&1) || (triggerMask&2)))
        return kTRUE;
      break;
    }
    case kMB2:
    {
      if (triggerMask&32 && ((triggerMask&1) || (triggerMask&2)))
        return kTRUE;
      break;
    }
  }

  return kFALSE;
}

//____________________________________________________________________
Bool_t AliPWG0Helper::IsVertexReconstructed(const AliESD* aEsd)
{
  // see function with AliESDVertex argument

  const AliESDVertex* vtxESD = aEsd->GetVertex();
  return IsVertexReconstructed(vtxESD);
}

//____________________________________________________________________
Bool_t AliPWG0Helper::IsVertexReconstructed(const AliESDVertex* vtxESD)
{
  // checks if the vertex is reasonable
  //
  // this function needs the branches fSPDVertex*

  if (!vtxESD)
    return kFALSE;

  // the vertex should be reconstructed
  if (strcmp(vtxESD->GetName(), "default")==0)
    return kFALSE;

  Double_t vtx_res[3];
  vtx_res[0] = vtxESD->GetXRes();
  vtx_res[1] = vtxESD->GetYRes();
  vtx_res[2] = vtxESD->GetZRes();

  if (vtx_res[2]==0 || vtx_res[2]>0.1)
    return kFALSE;

  // check Ncontributors, if <0 it means error *gna*

  return kTRUE;
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

  if (strcmp(aParticle->GetName(),"XXX") == 0)
  {
    if (adebug)
      printf("WARNING: There is a particle named XXX.\n");
    return kFALSE;
  }

  TParticlePDG* pdgPart = aParticle->GetPDG();

  if (strcmp(pdgPart->ParticleClass(),"Unknown") == 0)
  {
    if (adebug)
      printf("WARNING: There is a particle with an unknown particle class (pdg code %d).\n", pdgCode);
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

//____________________________________________________________________
Int_t AliPWG0Helper::GetPythiaEventProcessType(AliHeader* aHeader, Bool_t adebug) {
  //
  // get the process type of the event.
  //

  // can only read pythia headers, either directly or from cocktalil header
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader->GenEventHeader());

  if (!pythiaGenHeader) {

    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(aHeader->GenEventHeader());
    if (!genCocktailHeader) {
      printf("AliPWG0Helper::GetProcessType : Unknown header type (not Pythia or Cocktail). \n");
      return -1;
    }

    TList* headerList = genCocktailHeader->GetHeaders();
    if (!headerList) {
      return -1;
    }

    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
      if (pythiaGenHeader)
        break;
    }

    if (!pythiaGenHeader) {
      printf("AliPWG0Helper::GetProcessType : Could not find Pythia header. \n");
      return -1;
    }
  }

  if (adebug) {
    printf("AliPWG0Helper::GetProcessType : Pythia process type found: %d \n",pythiaGenHeader->ProcessType());
  }

  return pythiaGenHeader->ProcessType();
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

void AliPWG0Helper::SetBranchStatusRecursive(TTree* tree, char *bname, Bool_t status, Bool_t debug)
{
  // Function  to switch on/off all data members of a top level branch
  // this is needed for branches without a trailing dot ".", for those
  // the root functionality with regular expressions does not work.
  // Usage e.g.
  // chain->SetBranchStatus("*", 0); 
  // SetBranchStatusRecursive(chain,"SPDVertex",1);
  // You need to give the full name of the top level branch zou want to access
  //==========================================================================
  // Author Christian.Klein-Boesing@cern.ch

  if (!tree)
    return;

  TBranch *br = tree->GetBranch(bname);
  if(!br) {
    Printf("AliPWG0Helper::SetBranchStatusRecursive: Branch %s not found", bname);
  }

  TObjArray *leaves = tree->GetListOfLeaves();
  Int_t nleaves = leaves->GetEntries();
  TLeaf *leaf = 0;
  TBranch *branch = 0;
  TBranch *mother = 0;
  for (Int_t i=0;i<nleaves;i++)  {
    // the matched entry e.g. SPDVertex is its own Mother
    leaf = (TLeaf*)leaves->UncheckedAt(i);
    branch = (TBranch*)leaf->GetBranch();
    mother = branch->GetMother();
    if (mother==br){
      if (debug)
        Printf(">>>> AliPWG0Helper::SetBranchStatusRecursive: Setting status %s %s to %d", mother->GetName(), leaf->GetName(), status);
      if (status) branch->ResetBit(kDoNotProcess);
      else        branch->SetBit(kDoNotProcess);
    }
  }
}
