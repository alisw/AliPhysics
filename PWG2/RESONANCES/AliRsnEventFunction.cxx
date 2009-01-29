//
// Class AliRsnEventFunction
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
#include "AliRsnCutSet.h"

#include "AliRsnEventFunction.h"

ClassImp(AliRsnEventFunction)

//________________________________________________________________________________________
AliRsnEventFunction::AliRsnEventFunction() :
  fType(kTypes),
  fPIDMethod(AliRsnDaughter::kNoPID),
  fPIDType(AliRsnPID::kUnknown),
  fCharge('0'),
  fLeadPtMin(0.0),
  fAccept(kFALSE),
  fUseBins(kFALSE),
  fBins(0),
  fBinningCut(),
  fBinningCutType(AliRsnCut::kLastCutType),
  fEventCuts(0x0),
  fTrackCuts(0x0),
  fHistoDef(0x0)
{
  //
  // Constructor.
  // The histogram data member cannot be passed externally,
  // its initialization MUST be defined inside the Init() method,
  // which must be overridden in any derivate implementation.
  //

  Int_t i;
  for (i = 0; i < 100; i++)
  {
    fHisto[i] = 0x0;
  }
}

//________________________________________________________________________________________
AliRsnEventFunction::AliRsnEventFunction
(EType type, AliRsnHistoDef *hd, AliRsnDaughter::EPIDMethod pidMethod,
 AliRsnPID::EType pidType, Char_t sign) :
  fType(type),
  fPIDMethod(pidMethod),
  fPIDType(pidType),
  fCharge(sign),
  fLeadPtMin(0.0),
  fAccept(kFALSE),
  fUseBins(kFALSE),
  fBins(0),
  fBinningCut(),
  fBinningCutType(AliRsnCut::kLastCutType),
  fEventCuts(0x0),
  fTrackCuts(0x0),
  fHistoDef(hd)
{
  //
  // Constructor.
  // The histogram data member cannot be passed externally,
  // its initialization MUST be defined inside the Init() method,
  // which must be overridden in any derivate implementation.
  //

  Int_t i;
  for (i = 0; i < 100; i++)
  {
    fHisto[i] = 0x0;
  }
}

//________________________________________________________________________________________
AliRsnEventFunction::AliRsnEventFunction(const AliRsnEventFunction &copy) :
  TObject(copy),
  fType(copy.fType),
  fPIDMethod(copy.fPIDMethod),
  fPIDType(copy.fPIDType),
  fCharge(copy.fCharge),
  fLeadPtMin(copy.fLeadPtMin),
  fAccept(kFALSE),
  fUseBins(copy.fUseBins),
  fBins(0),
  fBinningCut(),
  fBinningCutType(AliRsnCut::kLastCutType),
  fEventCuts(copy.fEventCuts),
  fTrackCuts(copy.fTrackCuts),
  fHistoDef(copy.fHistoDef)
{
  //
  // Constructor.
  // The histogram data member cannot be passed externally,
  // its initialization MUST be defined inside the Init() method,
  // which must be overridden in any derivate implementation.
  //

  Int_t i, n = 100;
  for (i = 0; i < n; i++)
  {
    fHisto[i] = 0x0;
  }

  if (fUseBins)
  {
    n = copy.fBins.GetSize();
    Double_t *array = new Double_t[n];
    for (i = 0; i < n; i++) array[i] = copy.fBins[i];
    SetBinningCut(copy.fBinningCutType, copy.fBins.GetSize(), array);
    delete [] array;
  }
}

//________________________________________________________________________________________
void AliRsnEventFunction::Clear(Option_t* /*option*/)
{
  //
  // Clear arrays and histogram.
  // For the sake of security, all pointers are also set explicitly to NULL.
  //

  Int_t i;
  for (i = 0; i < 100; i++)
  {
    delete fHisto[i];
    fHisto[i] = 0x0;
  }
}

//________________________________________________________________________________________
void AliRsnEventFunction::Init(TList *histos)
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
  fHisto[0] = new TH1D(GetFcnName(), "", nbins, min, max);
  histos->AddLast(fHisto[0]);

  // if requested a binning w.r. to some cut variable, histograms are added
  // for that in this part of the method (one per each bin)
  Char_t hName[255];
  if (fUseBins)
  {
    for (ibin = 0, i = 1; ibin < fBins.GetSize() - 1; ibin++, i++)
    {
      sprintf(hName, "%s[%.2f-%.2f]", GetFcnName().Data(), fBins[ibin], fBins[ibin+1]);
      fHisto[i] = new TH1D(hName, "", nbins, min, max);
      histos->AddLast(fHisto[i]);
    }
  }
}

//________________________________________________________________________________________
void AliRsnEventFunction::SetBinningCut
(AliRsnCut::EType type, Double_t min, Double_t max, Double_t step)
{
  //
  // Set fixed bins
  //

  fUseBins = kTRUE;

  Int_t i, nBins = (Int_t)((max - min) / step) + 1;
  fBinningCutType = type;
  fBins.Set(nBins);
  for (i = 0; i < nBins; i++)
  {
    fBins[i] = min + (Double_t)i * step;
  }
}

//________________________________________________________________________________________
void AliRsnEventFunction::SetBinningCut
(AliRsnCut::EType type, Int_t nbins, Double_t *bins)
{
  //
  // Set variable bins
  //

  fUseBins = kTRUE;

  Int_t i;
  fBinningCutType = type;
  fBins.Set(nbins);
  for (i = 0; i < nbins; i++)
  {
    fBins[i] = bins[i];
  }
}

//________________________________________________________________________________________
TString AliRsnEventFunction::GetFcnName()
{
  //
  // Return a string which names the function type
  //

  TString text("Undef");

  switch (fType)
  {
    case kMultiplicity:
      text = "MULT";
      break;
    case kLeadingMomentum:
      text = "PLEAD";
      break;
    case kLeadingTheta:
      text = "LEADTHETA";
      break;
    case kAverageMomentum:
      text = "PAVG";
      break;
    case kAngleLeadingMean:
      text = "ANGLEADMEAN";
      break;
    case kAngleLeadingRMS:
      text = "ANGLEADRMS";
      break;
    case kVtResolution:
      text = "VTRES";
      break;
    case kVzResolution:
      text = "VZRES";
      break;
    default:
      AliError("Type not defined");
  }

  switch (fPIDMethod)
  {
    case AliRsnDaughter::kNoPID:
      text += "_NOPID_";
      break;
    case AliRsnDaughter::kPerfect:
      text += "_PERFECT_";
      break;
    case AliRsnDaughter::kRealistic:
      text += "_REALISTIC_";
      break;
    default:
      AliError("PID method not defined");
  }

  text += AliRsnPID::ParticleName(fPIDType);
  text += fCharge;

  if (fEventCuts) {
    text += '_';
    text += fEventCuts->GetName();
  }
  if (fTrackCuts) {
    text += '_';
    text += fTrackCuts->GetName();
  }

  return text;
}

//________________________________________________________________________________________
Bool_t AliRsnEventFunction::Fill(AliRsnEvent *event)
{
  //
  // Fillse the histogram with data contained in a defined pair.
  // This method must be overidden by an appropriate definition in each inheriting class.
  //

  if (fEventCuts) if (!fEventCuts->IsSelected(AliRsnCut::kEvent, event)) return kFALSE;

  // first of all, set all selection definitions, using the ones in this object
  event->SetSelectionPIDType(fPIDType);
  event->SetSelectionCharge(fCharge);
  event->SetSelectionPIDMethod(fPIDMethod);
  event->SetSelectionTrackCuts(fTrackCuts);

  Double_t value = FcnValue(event);
  if (!fAccept) return kFALSE;

  // fill global histogram
  fHisto[0]->Fill(value);

  // if bins are allocated, find right one and fill it
  if (fUseBins)
  {
    Int_t i, ibin;
    for (ibin = 0, i = 1; ibin < fBins.GetSize() - 1; ibin++, i++)
    {
      if (!fHisto[i]) continue;
      fBinningCut.SetCutValues(fBinningCutType, (Double_t)fBins[ibin], (Double_t)fBins[ibin+1]);
      fBinningCut.SetCutValues(fBinningCutType, (Int_t)fBins[ibin], (Int_t)fBins[ibin+1]);
      if (fBinningCut.IsSelected(AliRsnCut::kEvent, event))
      {
        fHisto[i]->Fill(value);
        break;
      }
    }
  }

  return kTRUE;
}

//________________________________________________________________________________________
Double_t AliRsnEventFunction::FcnValue(AliRsnEvent *event)
{
  //
  // This method must be overridden in all inheritin functions.
  // It computes the value which must be used to fill the histogram.
  //

  Int_t count;
  Double_t output, v[3], vMC[3], vt, vtMC, mean = 0.0, rms = 0.0;
  AliRsnDaughter *trk = 0x0;

  switch (fType)
  {
    case kMultiplicity:
      fAccept = kTRUE;
      output = (Double_t)event->GetMultiplicity();
      break;
    case kLeadingMomentum:
      trk = event->GetLeadingParticle(fLeadPtMin);
      if (trk) {
        fAccept = kTRUE;
        output = trk->P();
      }
      else {
        fAccept = kFALSE;
        output = 0.0;
      }
      break;
    case kLeadingTheta:
      trk = event->GetLeadingParticle(fLeadPtMin);
      if (trk) {
        fAccept = kTRUE;
        output = trk->Theta();
      }
      else {
        fAccept = kFALSE;
        output = 0.0;
      }
      break;
    case kAverageMomentum:
      output = event->GetAverageMomentum(count);
      fAccept = (count > 0);
      break;
    case kAngleLeadingMean:
    case kAngleLeadingRMS:
      fAccept = event->GetAngleDistrWRLeading(mean, rms, fLeadPtMin);
      if (fType == kAngleLeadingMean) output = mean; else output = rms;
      break;
    case kVzResolution:
    case kVtResolution:
      fAccept = kTRUE;
      v[0] = event->GetPrimaryVertexX();
      v[1] = event->GetPrimaryVertexY();
      v[2] = event->GetPrimaryVertexZ();
      vMC[0] = event->GetPrimaryVertexXMC();
      vMC[1] = event->GetPrimaryVertexYMC();
      vMC[2] = event->GetPrimaryVertexZMC();
      vt = TMath::Sqrt(v[0]*v[0] + v[1]*v[1]);
      vtMC = TMath::Sqrt(vMC[0]*vMC[0] + vMC[1]*vMC[1]);
      if (fType == kVtResolution) return (vt - vtMC); else return (v[2] - vMC[2]);
      break;
    default:
      AliError(Form("Type '%d' not supported for EVENT functions", fType));
      fAccept = kFALSE;
      output = 0.0;
  }

  return output;
}

