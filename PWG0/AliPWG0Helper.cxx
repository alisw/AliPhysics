/* $Id$ */

#include <AliPWG0Helper.h>

#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1.h>
#include <TH3.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliESDVertex.h>

//____________________________________________________________________
ClassImp(AliPWG0Helper)

//____________________________________________________________________
Bool_t AliPWG0Helper::IsEventTriggered(AliESD* aEsd)
{
  // check if the event was triggered
  //
  // this function needs the branch fTriggerMask
  //
  // MB should be
  // ITS_SPD_GFO_L0  : 32
  // VZERO_OR_LEFT   : 1
  // VZERO_OR_RIGHT  : 2

  ULong64_t triggerMask = aEsd->GetTriggerMask();

  if (triggerMask&32 && ((triggerMask&1) || (triggerMask&2)))
    return kTRUE;

  return kFALSE;
}

//____________________________________________________________________
Bool_t AliPWG0Helper::IsVertexReconstructed(AliESD* aEsd)
{
  // checks if the vertex is reasonable
  //
  // this function needs the branches fSPDVertex*


  const AliESDVertex* vtxESD = aEsd->GetVertex();

  // the vertex should be reconstructed
  if (strcmp(vtxESD->GetName(), "default")==0)
    return kFALSE;

  Double_t vtx_res[3];
  vtx_res[0] = vtxESD->GetXRes();
  vtx_res[1] = vtxESD->GetYRes();
  vtx_res[2] = vtxESD->GetZRes();

  if (vtx_res[2]==0 || vtx_res[2]>0.1)
    return kFALSE;

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
void AliPWG0Helper::CreateProjections(TH3* hist)
{
  // create projections of 3d hists to all 2d combinations
  // the histograms are not returned, just use them from memory or use this to create them in a file

  TH1* proj = hist->Project3D("yx");
  proj->SetXTitle(hist->GetXaxis()->GetTitle());
  proj->SetYTitle(hist->GetYaxis()->GetTitle());

  proj = hist->Project3D("zx");
  proj->SetXTitle(hist->GetXaxis()->GetTitle());
  proj->SetYTitle(hist->GetZaxis()->GetTitle());

  proj = hist->Project3D("zy");
  proj->SetXTitle(hist->GetYaxis()->GetTitle());
  proj->SetYTitle(hist->GetZaxis()->GetTitle());
}

//____________________________________________________________________
void AliPWG0Helper::CreateDividedProjections(TH3* hist, TH3* hist2, const char* axis, Bool_t putErrors)
{
  // create projections of the 3d hists divides them
  // axis decides to which plane, if axis is 0 to all planes
  // the histograms are not returned, just use them from memory or use this to create them in a file

  if (axis == 0)
  {
    CreateDividedProjections(hist, hist2, "yx", putErrors);
    CreateDividedProjections(hist, hist2, "zx", putErrors);
    CreateDividedProjections(hist, hist2, "zy", putErrors);

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
  division->Divide(proj2);

  if (putErrors)
  {
    division->Sumw2();
    if (division->GetDimension() == 1)
    {
      Int_t nBins = division->GetNbinsX();
      for (Int_t i = 0; i <= nBins; ++i)
        if (proj2->GetBinContent(i) != 0)
          division->SetBinError(i, TMath::Sqrt(proj->GetBinContent(i)) / proj2->GetBinContent(i));
    }
    else if (division->GetDimension() == 2)
    {
      Int_t nBinsX = division->GetNbinsX();
      Int_t nBinsY = division->GetNbinsY();
      for (Int_t i = 0; i <= nBinsX; ++i)
        for (Int_t j = 0; j <= nBinsY; ++j)
          if (proj2->GetBinContent(i, j) != 0)
            division->SetBinError(i, j, TMath::Sqrt(proj->GetBinContent(i, j)) / proj2->GetBinContent(i, j));
    }
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
