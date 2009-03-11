//
// Class AliRsnFunctionDef
//
// This class defines a base classe to implement a function
// which uses the internal RSN package event format (AliRsnEvent).
// It contains some default flags which turn out to be useful:
//  - a flag to select only the "true" pairs (tracks from same resonance)
//  - a flag to know if the computation is done over two events (mixing)
//
// Any kind of analysis object should be implemented as inheriting from this
// because the AliRsnAnalyzer which executes the analysis will accept a collection
// of such objects, in order to have a unique format of processing method
//
// The user who implements a kind of computation type should inherit from
// this class and override the virtual functions defined in it, which
// initialize the final output histogram and define how to process data.
//
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//

#include <Riostream.h>

#include <TH1.h>
#include <TList.h>
#include <TString.h>

#include "AliLog.h"

#include "AliRsnPairParticle.h"
#include "AliRsnPairDef.h"
//#include "AliRsnCut.h"
//#include "AliRsnMCInfo.h"
//#include "AliRsnDaughter.h"
//#include "AliRsnEvent.h"


#include "AliRsnFunctionDef.h"

ClassImp(AliRsnFunctionDef)

//________________________________________________________________________________________
AliRsnFunctionDef::AliRsnFunctionDef() :
  fFcnType(kFcnTypes),
  fNBinsX(0),
  fXmin(0.0),
  fXmax(0.0),
  fYUsed(0),
  fRotAngle(0.0)
{
  //
  // Constructor.
  // The histogram data member cannot be passed externally,
  // its initialization MUST be defined inside the Init() method,
  // which must be overridden in any derivate implementation.
  //

  Int_t j;
  for (j = 0 ; j < kBinMax; j++)
  {
    fNBinsY[j] = 0;
    fYbins[j].Set(0);
  }
}

//________________________________________________________________________________________
AliRsnFunctionDef::AliRsnFunctionDef
(EFcnType type, Int_t nbins, Double_t min, Double_t max) :
  fFcnType(type),
  fNBinsX(nbins),
  fXmin(min),
  fXmax(max),
  fYUsed(0),
  fRotAngle(0.0)
{
  //
  // Constructor.
  // The histogram data member cannot be passed externally,
  // its initialization MUST be defined inside the Init() method,
  // which must be overridden in any derivate implementation.
  //

  Int_t j;
  for (j = 0 ; j < kBinMax; j++)
  {
    fNBinsY[j] = 0;
    fYbins[j].Set(0);
  }
}

//________________________________________________________________________________________
AliRsnFunctionDef::AliRsnFunctionDef(const AliRsnFunctionDef &copy) :
  TObject(copy),
  fFcnType(copy.fFcnType),
  fNBinsX(copy.fNBinsX),
  fXmin(copy.fXmin),
  fXmax(copy.fXmax),
  fYUsed(copy.fYUsed),
  fRotAngle(copy.fRotAngle)
{
  //
  // Copy constructor.
  // Calls the functions to define binning.
  //

  Int_t i, j, nbins;

  for (i = 0; i < fYUsed; i++)
  {
    nbins = copy.fNBinsY[i];
    if (!nbins) continue;
    fYbins[i].Set(nbins);
    for (j = 0; j < nbins; j++) fYbins[i][j] = copy.fYbins[i][j];
  }
}
//________________________________________________________________________________________
const AliRsnFunctionDef& AliRsnFunctionDef::operator=(const AliRsnFunctionDef& copy)
{
  //
  // Assignment operator.
  // Behaves like copy constructor.
  // Also in this case, the histogram is not copied, and,
  // if it was present, it is destroyed and will need to be recreated.
  //

  fFcnType = copy.fFcnType;
  fNBinsX = copy.fNBinsX;
  fXmin = copy.fXmin;
  fXmax = copy.fXmax;
  fYUsed = copy.fYUsed;
  fRotAngle = copy.fRotAngle;

  Int_t i, j, nbins;

  for (i = 0; i < fYUsed; i++)
  {
    nbins = copy.fNBinsY[i];
    if (!nbins) continue;
    fYbins[i].Set(nbins);
    for (j = 0; j < nbins; j++) fYbins[i][j] = copy.fYbins[i][j];
  }

  return (*this);
}

//________________________________________________________________________________________
void AliRsnFunctionDef::SetBinningY
(Int_t i, EFcnBinType type, Double_t xmin, Double_t xmax, Double_t step)
{
//
// Define a secondary binning with fixed bins of given step
//

  if (i >= fYUsed) return;

  Int_t ib;

  fBinType[i] = type;
  fNBinsY[i] = (Int_t)((xmax - xmin) / step) + 2;
  fYbins[i].Set(fNBinsY[i]);

  for (ib = 0; ib < fNBinsY[i]; ib++)
  {
    fYbins[i][ib] = xmin + (Double_t)ib * step;
  }
}

//________________________________________________________________________________________
void AliRsnFunctionDef::SetBinningY
(Int_t i, EFcnBinType type, Int_t nbins, Double_t *bins)
{
//
// Define a secondary binning with variable bins
//

  if (i >= fYUsed) return;

  fBinType[i] = type;

  Int_t ib;

  fYbins[i].Set(nbins);
  fNBinsY[i] = nbins;
  for (ib = 0; ib < nbins; ib++)
  {
    fYbins[i][ib] = bins[ib];
  }
}

//________________________________________________________________________________________
Bool_t AliRsnFunctionDef::AddBinningY(EFcnBinType type, Double_t xmin, Double_t xmax, Double_t step)
{
//
// Add a new secondary binning with fixed bins of given step
//

  if (fYUsed >= kBinMax)
  {
    AliError(Form("Cannot create a new binning: maximum allowed of %d is reached", kBinMax));
    return kFALSE;
  }

  fYUsed++;
  SetBinningY(fYUsed - 1, type, xmin, xmax, step);

  return kTRUE;
}

//________________________________________________________________________________________
Bool_t AliRsnFunctionDef::AddBinningY(EFcnBinType type, Int_t nbins, Double_t *bins)
{
//
// Define a secondary binning with variable bins
//

  if (fYUsed >= kBinMax)
  {
    AliError(Form("Cannot create a new binning: maximum allowed of %d is reached", kBinMax));
    return kFALSE;
  }

  fYUsed++;
  SetBinningY(fYUsed - 1, type, nbins, bins);

  return kTRUE;
}

//________________________________________________________________________________________
TString AliRsnFunctionDef::GetFcnName()
{
  //
  // Return a string which names the function type
  //

  TString text("Undef");

  switch (fFcnType)
  {
    case kInvMass:
      text = "IM";
      break;
    case kInvMassMC:
      text = "IM_MC";
      break;
    case kInvMassRotated:
      text = Form("IMR%.2f", fRotAngle);
      break;
    case kResolution:
      text = "RES";
      break;
    case kPtSpectrum:
      text = "PT";
      break;
    case kEtaSpectrum:
      text = "ETA";
      break;
    default:
      AliError("Type not defined");
  }

  return text;
}

//________________________________________________________________________________________
TString AliRsnFunctionDef::GetFcnTitle()
{
  //
  // Return a string which names the function type
  //

  TString text("Undef");

  switch (fFcnType)
  {
    case kInvMass:
      text = "Invariant mass";
      break;
    case kInvMassMC:
      text = "Invariant mass (MC)";
      break;
    case kResolution:
      text = "Resolution";
      break;
    case kPtSpectrum:
      text = "p_{#perp} distribution";
      break;
    case kEtaSpectrum:
      text = "#eta distribution";
      break;
    default:
      AliError("Type not defined");
  }

  return text;
}

//________________________________________________________________________________________
TString AliRsnFunctionDef::GetBinningName(Int_t i)
{
  //
  // Return a string which names the function type
  //

  TString text("Undef");

  if (i >= kBinMax) return text;

  switch (fBinType[i])
  {
    case kPt:
      text = "PT";
      break;
    case kEta:
      text = "ETA";
      break;
    default:
      AliError("Type not defined");
  }

  return text;
}

//________________________________________________________________________________________
Double_t AliRsnFunctionDef::EvalX(AliRsnPairParticle *pair, AliRsnPairDef *ref)
{
//
// Evaluate X value according to selected type
//

  Double_t recInvMass = pair->GetInvMass(ref->GetMass(0), ref->GetMass(1));
  Double_t simInvMass = pair->GetInvMassMC(ref->GetMass(0), ref->GetMass(1));

  switch (fFcnType)
  {
    case kInvMass:
      return recInvMass;
    case kInvMassMC:
      return simInvMass;
    case kInvMassRotated:
      pair->GetDaughter(1)->RotateP(fRotAngle * TMath::DegToRad());
      pair->ResetPair();
      return pair->GetInvMass(ref->GetMass(0), ref->GetMass(1));
    case kResolution:
      return (simInvMass - recInvMass) / simInvMass;
    case kPtSpectrum:
      return pair->GetPt();
    case kEtaSpectrum:
      return pair->GetEta();
    default:
      AliError("Type not defined");
      return 1000.0;
  }
}

//________________________________________________________________________________________
Double_t AliRsnFunctionDef::EvalY(Int_t i, AliRsnPairParticle *pair)
{
//
// Evaluate Y value according to selected type
//

  if (i >= fYUsed) return 1000.0;

  switch (fBinType[i])
  {
    case kPt:
      return pair->GetPt();
    case kEta:
      return pair->GetEta();
    default:
      AliError("Binning type not defined");
      return 1000.0;
  }
}

//________________________________________________________________________________________
void AliRsnFunctionDef::Print(Option_t* /*opt*/)
{
//
// Print details on this function
//

  Int_t i, j;

  cout << "DEFINITION X TYPE: " << GetFcnName().Data() << endl;
  cout << "DEFINITION Y TYPES = " << fYUsed << endl;

  for (i = 0; i < fYUsed; i++) {
    cout << "Y TYPE #" << i << ": " << GetBinningName(i) << " [bins = ";
    for (j = 0; j < fYbins[i].GetSize(); j++) {
      cout << fYbins[i][j] << ' ';
    }
    cout << endl;
  }
}
