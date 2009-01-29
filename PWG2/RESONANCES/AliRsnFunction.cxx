//
// Class AliRsnFunction
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

#include "AliRsnMCInfo.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnPairDef.h"
#include "AliRsnCut.h"

#include "AliRsnFunction.h"

ClassImp(AliRsnFunction)

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction() :
    fFcnType(kFcnTypes),
    fRotAngle(0.0),
    fUseBins(kFALSE),
    fSkipOutsideInterval(kFALSE),
    fNumberOfBinTypes(0),
//     fBins(0),
//     fBinningCut(),
//     fBinningCutType(AliRsnCut::kLastCutType),
    fHistoDef(0x0)
{
  //
  // Constructor.
  // The histogram data member cannot be passed externally,
  // its initialization MUST be defined inside the Init() method,
  // which must be overridden in any derivate implementation.
  //

  Int_t i, j;
  for (j = 0 ; j < kFcnBinTypes; j++)
    for (i = 0; i < 100; i++)
      fHisto[j][i] = 0x0;

}

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction
(EFcnType type, AliRsnHistoDef *hd, Bool_t skipOut) :
    fFcnType(type),
    fRotAngle(0.0),
    fUseBins(kFALSE),
    fSkipOutsideInterval(skipOut),
    fNumberOfBinTypes(0),
//     fBins(0),
//     fBinningCut(),
//     fBinningCutType(AliRsnCut::kLastCutType),
    fHistoDef(hd)
{
  //
  // Constructor.
  // The histogram data member cannot be passed externally,
  // its initialization MUST be defined inside the Init() method,
  // which must be overridden in any derivate implementation.
  //

  Int_t i, j;
  for (j = 0 ; j < kFcnBinTypes; j++)
    for (i = 0; i < 100; i++)
      fHisto[j][i] = 0x0;
}

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction(const AliRsnFunction &copy) :
    TObject(copy),
    fFcnType(copy.fFcnType),
    fRotAngle(copy.fRotAngle),
    fUseBins(copy.fUseBins),
    fSkipOutsideInterval(copy.fSkipOutsideInterval),
    fNumberOfBinTypes(copy.fNumberOfBinTypes),
//     fBins(0),
//     fBinningCut(),
//     fBinningCutType(AliRsnCut::kLastCutType),
    fHistoDef(copy.fHistoDef)
{
  //
  // Copy constructor.
  // Calls the function to define binning.
  //

  Int_t i, j, n;
  for (j = 0 ; j < kFcnBinTypes; j++)
    for (i = 0; i < 100; i++)
      fHisto[j][i] = 0x0;

  if (fUseBins)
  {
    for (i = 0 ; i < kFcnBinTypes; i++){
      if (fNumberOfBinTypes<=i) continue;
      n = copy.fBins[i].GetSize();
      Double_t *array = new Double_t[n];
      for (j = 0; j < n; j++) array[j] = copy.fBins[i][j];
      SetBinningCut(copy.fBinningCutType[i], copy.fBins[i].GetSize(), array,i,kTRUE);
      delete [] array;
    }
  }
}
//________________________________________________________________________________________
const AliRsnFunction& AliRsnFunction::operator=(const AliRsnFunction& /*copy*/)
{
  //
  // Assignment operator.
  // Behaves like copy constructor.
  // Also in this case, the histogram is not copied, and,
  // if it was present, it is destroyed and will need to be recreated.
  //

  return (*this);
}
//________________________________________________________________________________________
void AliRsnFunction::Clear(Option_t* /*option*/)
{
  //
  // Clear arrays and histogram.
  // For the sake of security, all pointers are also set explicitly to NULL.
  //

  Int_t i, j;
  for (j = 0 ; j < kFcnBinTypes; j++)
    for (i = 0; i < 100; i++)
    {
      delete fHisto[j][i];
      fHisto[j][i] = 0x0;
    }
}

//________________________________________________________________________________________
TList* AliRsnFunction::Init(const char *histoName, const char *histoTitle)
{
  //
  // Initialization function.
  // By default, it initializes the owned histogram using the method
  // from AliRsnHistoDef class, giving the same name and title of this.
  // A user can override this behaviour, if necessary.
  // Before creating, the HistoDef is checked for proper initialization.
  //

  Clear();

  Int_t i, ibin, nbins = fHistoDef->GetNBins();
  Double_t min = fHistoDef->GetMin(), max = fHistoDef->GetMax();

  // list is created and named after the general
  // settings used for the contained histograms
  TList *histos = new TList;
  histos->SetName(Form("%s", GetFcnName().Data()));

  // a general histogram is always added,
  // which overrides the binning and collects everything

  fHisto[0][0] = new TH1D(histoName, histoTitle, nbins, min, max);
  fHisto[0][0]->Sumw2();
  histos->AddLast(fHisto[0][0]);

  // if requested a binning w.r. to some cut variable, histograms are added
  // for that in this part of the method (one per each bin)
  Char_t hName[255];
  Char_t hTitle[255];
  Int_t j;
  if (fUseBins)
  {
    for (j = 0 ; j < kFcnBinTypes; j++){
      if (fNumberOfBinTypes<=j) continue;
      for (ibin = 0, i = 1; ibin < fBins[j].GetSize() - 1; ibin++, i++)
      {
	sprintf(hName, "%s_%d%02d_[%.2f-%.2f]", histoName, j,i,fBins[j][ibin], fBins[j][ibin+1]);
	sprintf(hTitle, "%s [%.2f-%.2f]", histoTitle, fBins[j][ibin], fBins[j][ibin+1]);
// 	AliInfo(Form("Adding %s",hName));
	fHisto[j][i] = new TH1D(hName, hTitle, nbins, min, max);
	fHisto[j][i]->Sumw2();
	histos->AddLast(fHisto[j][i]);
      }
    }
  }

  // returns the full list at the end
  return histos;
}

//________________________________________________________________________________________
void AliRsnFunction::Init(const char *histoName, const char *histoTitle, TList *histos)
{
  //
  // Initialization function.
  // By default, it initializes the owned histogram using the method
  // from AliRsnHistoDef class, giving the same name and title of this.
  // A user can override this behaviour, if necessary.
  // Before creating, the HistoDef is checked for proper initialization.
  //

  Clear();

  Int_t i, ibin, nbins = fHistoDef->GetNBins();
  Double_t min = fHistoDef->GetMin(), max = fHistoDef->GetMax();

  // list is created and named after the general
  // settings used for the contained histograms
  if (!histos)
  {
    AliError("NULL target list!");
    return;
  }

  // a general histogram is always added,
  // which overrides the binning and collects everything
  fHisto[0][0] = new TH1D(histoName, histoTitle, nbins, min, max);
  histos->AddLast(fHisto[0][0]);

  // if requested a binning w.r. to some cut variable, histograms are added
  // for that in this part of the method (one per each bin)
  Char_t hName[255];
  Char_t hTitle[255];
  Int_t j;
  if (fUseBins)
  {
    for (j = 0 ; j < kFcnBinTypes; j++){
      if (fNumberOfBinTypes<=j) continue;
      for (ibin = 0, i = 1; ibin < fBins[j].GetSize() - 1; ibin++, i++)
      {

	sprintf(hName, "%s_%d_%02d[%.2f-%.2f]", histoName,j,i, fBins[j][ibin], fBins[j][ibin+1]);
	sprintf(hTitle, "%s [%.2f-%.2f]", histoTitle, fBins[j][ibin], fBins[j][ibin+1]);
// 	AliInfo(Form("Adding %s",hName));
	fHisto[j][i] = new TH1D(hName, hTitle, nbins, min, max);
	histos->AddLast(fHisto[j][i]);
      }
    }
  }
}

//________________________________________________________________________________________
void AliRsnFunction::SetBinningCut
(AliRsnCut::EType type, Double_t min, Double_t max, Double_t step,Int_t index,Bool_t IsCopyConstructor)
{
  //
  // Set fixed bins
  //

  if (index >= kFcnBinTypes) {
    AliError(Form("We support only %d Binning cuts(0-%d). Skipping...",kFcnBinTypes,kFcnBinTypes-1));
    return;
  }
  
  if (!IsCopyConstructor){
    // TODO if some one sets indexes 0,2,3 it is a bug here(i'll solve it)
    if (index == fNumberOfBinTypes)
      fNumberOfBinTypes++;
    else {
      AliError(Form("Wrong index %d. fUseBins is set to kFALSE",index));
//       fUseBins = kFALSE;
      return;
    }
  }

  fUseBins = kTRUE;
  Int_t i, nBins = (Int_t)((max - min) / step) + 1;
  fBinningCutType[index] = type;
  fBins[index].Set(nBins);
  for (i = 0; i < nBins; i++)
  {
    fBins[index][i] = min + (Double_t)i * step;
  }
}

//________________________________________________________________________________________
void AliRsnFunction::SetBinningCut
(AliRsnCut::EType type, Int_t nbins, Double_t *bins,Int_t index,Bool_t IsCopyConstructor)
{
  //
  // Set variable bins
  //

  if (index >= kFcnBinTypes) {
    AliError(Form("We support only %d Binning cuts(0-%d). Skipping...",kFcnBinTypes,kFcnBinTypes-1));
    return;
  }
   if (!IsCopyConstructor){
     // TODO if some one sets indexes 0,2,3 it is a bug here(i'll solve it)
     if (index >= fNumberOfBinTypes)
      fNumberOfBinTypes++;
     else {
        AliError(Form("Wrong index %d (%d). fUseBins is set to kFALSE",index,fNumberOfBinTypes));
 //        fUseBins = kFALSE;
        return;
      }
   }

  fUseBins = kTRUE;
  Int_t i;
  fBinningCutType[index] = type;
  fBins[index].Set(nbins);
  for (i = 0; i < nbins; i++)
  {
    fBins[index][i] = bins[i];
  }
}

//________________________________________________________________________________________
TString AliRsnFunction::GetFcnName()
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
TString AliRsnFunction::GetFcnTitle()
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
Bool_t AliRsnFunction::Fill(AliRsnPairParticle *pair, AliRsnPairDef *ref)
{
  //
  // Fillse the histogram with data contained in a defined pair.
  // This method must be overidden by an appropriate definition in each inheriting class.
  //

  Double_t value = FcnValue(pair, ref);
  if (fSkipOutsideInterval)
  {
    if (value < fHistoDef->GetMin()) return kFALSE;
    if (value > fHistoDef->GetMax()) return kFALSE;
  }

  // fill global histogram
  fHisto[0][0]->Fill(value);

  // if bins are allocated, find right one and fill it
  if (fUseBins)
  {
    Int_t i, j, ibin;
    for (j = 0 ; j < kFcnBinTypes; j++){
      if (fNumberOfBinTypes<=j) continue;
      for (ibin = 0, i = 1; ibin < fBins[j].GetSize() - 1; ibin++, i++)
      {
	if (!fHisto[j][i]) continue;
	fBinningCut[j].SetCutValues(fBinningCutType[j], fBins[j][ibin], fBins[j][ibin+1]);
	if (fBinningCut[j].IsSelected(AliRsnCut::kPair, pair))
	{
	  fHisto[j][i]->Fill(value);
	  break;
	}
      }
    }
  }

  return kTRUE;
}

//________________________________________________________________________________________
Double_t AliRsnFunction::FcnValue(AliRsnPairParticle *pair, AliRsnPairDef *ref)
{
  //
  // This method must be overridden in all inheritin functions.
  // It computes the value which must be used to fill the histogram.
  //

  switch (fFcnType)
  {
    case kInvMass:
      return pair->GetInvMass(ref->GetMass(0), ref->GetMass(1));
    case kInvMassMC:
      return pair->GetInvMassMC(ref->GetMass(0), ref->GetMass(1));
    case kInvMassRotated:
      //AliInfo(Form("*** ROTATION ANGLE = %f ***", fRotAngle));
      //AliInfo(Form("UNROTATED INV MASS = %f", pair->GetInvMass(ref->GetMass(0), ref->GetMass(1))));
      //pair->GetDaughter(1)->Print("P");
      pair->GetDaughter(1)->RotateP(fRotAngle * TMath::DegToRad());
      pair->ResetPair();
      //AliInfo(Form("  ROTATED INV MASS = %f", pair->GetInvMass(ref->GetMass(0), ref->GetMass(1))));
      //pair->GetDaughter(1)->Print("P");
      return pair->GetInvMass(ref->GetMass(0), ref->GetMass(1));
    case kResolution:
      return FcnResolution(pair, ref);
    case kPtSpectrum:
      return pair->GetPt();
    case kEtaSpectrum:
      return pair->GetEta();
    default:
      AliError("Type not defined");
  }

  return 0.0;
}

//________________________________________________________________________________________
inline Double_t AliRsnFunction::FcnResolution(AliRsnPairParticle *pair, AliRsnPairDef *ref)
{
  //
  // Invariant mass resolution (compared between reconstructed and montecarlo)
  //

  Double_t recInvMass = pair->GetInvMass(ref->GetMass(0), ref->GetMass(1));
  Double_t simInvMass = pair->GetInvMassMC(ref->GetMass(0), ref->GetMass(1));

  return (simInvMass - recInvMass) / simInvMass;
}
