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
#include <TString.h>

#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnPairDef.h"
#include "AliRsnCut.h"

#include "AliRsnFunction.h"

ClassImp(AliRsnFunction)

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction() :
  TNamed(),
  fFcnType(kFcnTypes),
  fPairDef(0x0),
  fTrack(0x0),
  fPair(0x0),
  fEvent(0x0),
  fHistogram(0x0)
{
//
// Constructor for 1D functions.
// Requires only the binning of the output function,
// which is stored as 'main' histoDef in fHistoDef[0]
//

  fBinType[0] = kNoBins;
  fBinType[1] = kNoBins;

  fHistoDef[0] = 0x0;
  fHistoDef[1] = 0x0;
  fHistoDef[2] = 0x0;
}

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction
(EFcnType fcnType, AliRsnHistoDef *hd) :
  TNamed(),
  fFcnType(fcnType),
  fPairDef(0x0),
  fTrack(0x0),
  fPair(0x0),
  fEvent(0x0),
  fHistogram(0x0)
{
//
// Constructor for 1D functions.
// Requires only the binning of the output function,
// which is stored as 'main' histoDef in fHistoDef[0]
//

  fBinType[0] = kNoBins;
  fBinType[1] = kNoBins;

  fHistoDef[0] = hd;
  fHistoDef[1] = 0x0;
  fHistoDef[2] = 0x0;

  DefineName();
}

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction
(EFcnType fcnType, EBinType binType, AliRsnHistoDef *hdMain, AliRsnHistoDef *hdBin) :
  TNamed(),
  fFcnType(fcnType),
  fPairDef(0x0),
  fTrack(0x0),
  fPair(0x0),
  fEvent(0x0),
  fHistogram(0x0)
{
//
// Constructor for 2D functions.
// Requires the binning of the output function,
// which is stored as 'main' histoDef in fHistoDef[0],
// and a definition for a secondary binning, stored in fHistoDef[1]
//

  fBinType[0] = binType;
  fBinType[1] = kNoBins;

  fHistoDef[0] = hdMain;
  fHistoDef[1] = hdBin;
  fHistoDef[2] = 0x0;

  DefineName();
}

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction
(EFcnType fcnType, EBinType binType1, EBinType binType2,
 AliRsnHistoDef *hdMain, AliRsnHistoDef *hdBin1, AliRsnHistoDef *hdBin2) :
  fFcnType(fcnType),
  fPairDef(0x0),
  fTrack(0x0),
  fPair(0x0),
  fEvent(0x0),
  fHistogram(0x0)
{
//
// Constructor for 3D functions.
// Requires the binning of the output function,
// which is stored as 'main' histoDef in fHistoDef[0],
// and a definition for two secondary binnings, stored in fHistoDef[1,2]
//

  fBinType[0] = binType1;
  fBinType[1] = binType2;

  fHistoDef[0] = hdMain;
  fHistoDef[1] = hdBin1;
  fHistoDef[2] = hdBin2;

  DefineName();
}    

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction(const AliRsnFunction &copy) :
  TNamed(copy),
  fFcnType(copy.fFcnType),
  fPairDef(copy.fPairDef),
  fTrack(copy.fTrack),
  fPair(copy.fPair),
  fEvent(copy.fEvent),
  fHistogram(0x0)
{
//
// Copy constructor.
//

  Int_t i;
  for (i = 0; i < 3; i++) {
    fHistoDef[i] = copy.fHistoDef[i];
    if (i < 2) fBinType[i] = copy.fBinType[i];
  }

  DefineName();
}

//________________________________________________________________________________________
const AliRsnFunction& AliRsnFunction::operator=(const AliRsnFunction& copy)
{
//
// Assignment operator.
//

  SetName(copy.GetName());
  SetTitle(copy.GetTitle());

  fFcnType = copy.fFcnType;

  fPairDef = copy.fPairDef;

  Int_t i;
  for (i = 0; i < 3; i++) {
    fHistoDef[i] = copy.fHistoDef[i];
    if (i < 2) fBinType[i] = copy.fBinType[i];
  }

  fTrack = copy.fTrack;
  fPair = copy.fPair;
  fEvent = copy.fEvent;

  if (fHistogram) delete fHistogram;
  fHistogram = 0x0;

  DefineName();

  return (*this);
}

//________________________________________________________________________________________
const char* AliRsnFunction::FcnName()
{
//
// Defines the name of this object according to
// the function type and binning
//

  switch (fFcnType)
  {
    case kTrackPt:
      return  "TRKPT";
      break;
    case kTrackEta:
      return  "TRKETA";
      break;
    case kInvMass:
      return  "IM";
      break;
    case kInvMassMC:
      return  "IMMC";
      break;
    case kResolution:
      return  "RES";
      break;
    case kPairPt:
      return  "PT";
      break;
    case kPairEta:
      return  "ETA";
      break;
    case kEventMult:
      return  "MULT";
      break;
    default:
      return  "UNDEF";
  }
}

//________________________________________________________________________________________
void AliRsnFunction::DefineName()
{
//
// Defines the name of this object according to
// the function type and binning
//

  Int_t  dim = CheckDim();

  switch (dim)
  {
    case 1:
      SetName(FcnName());
      break;
    case 2:
      SetName(Form("%s_%s", FcnName(), BinName(fBinType[0])));
      break;
    case 3:
      SetName(Form("%s_%s_%s", FcnName(), BinName(fBinType[0]), BinName(fBinType[1])));
      break;
    default:
      SetName("UNDEF");
  }
}

//________________________________________________________________________________________
Double_t AliRsnFunction::Eval()
{
//
// Compute value for functions with 'event' argument type
//

  Double_t value;

  switch (fFcnType)
  {
    case kTrackPt:
      return fTrack->Pt();
    case kTrackEta:
      return fTrack->Eta();
    case kInvMass:
      return fPair->GetInvMass(fPairDef->GetMass(0), fPairDef->GetMass(1));
    case kInvMassMC:
      return fPair->GetInvMassMC(fPairDef->GetMass(0), fPairDef->GetMass(1));
    case kResolution:
      value  = fPair->GetInvMass(fPairDef->GetMass(0), fPairDef->GetMass(1));
      value -= fPair->GetInvMassMC(fPairDef->GetMass(0), fPairDef->GetMass(1));
      value /= fPair->GetInvMassMC(fPairDef->GetMass(0), fPairDef->GetMass(1));
      return value;
    case kPairPt:
      return fPair->GetPt();
    case kPairEta:
      return fPair->GetEta();
    case kEventMult:
      return fEvent->GetMultiplicity();
    default:
      AliWarning("Function type not supported");
      return -999.0;
  }
}

//_________________________________________________________________________________________________
Bool_t AliRsnFunction::CheckInput(Option_t *option)
{
//
// Checks if the argument type is coherent with
// the function type required
//

  TString opt(option);
  opt.ToUpper();

  if (opt.Contains("TRACK")) {
    if (!fTrack) {
      AliError("Input track object is NULL");
      return kFALSE;
    }
  }

  if (opt.Contains("PAIR")) {
    if (!fPair) {
      AliError("Input pair object is NULL");
      return kFALSE;
    }
  }

  if (opt.Contains("EVENT")) {
    if (!fEvent) {
      AliError("Input event object is NULL");
      return kFALSE;
    }
  }

  return kTRUE;
}

//________________________________________________________________________________________
TH1* AliRsnFunction::CreateHistogram(const char *histoName, const char *histoTitle)
{
//
// Creates and returns the histogram defined using
// arguments fo name and title, and the first histoDef for binning.
// Variable-sized histogram binning is always called, due to use of histoDef,
// even if the bins are equal, since they are defined in this class.
// Eventually present histoDef's in other slots of array (1, 2) are ignored.
//

  // first binning is required
  if (!fHistoDef[0]) return 0;

  // retrieve binnings for main and secondary axes
  Int_t    i, nbins[3] = {0, 0, 0};
  Double_t min[3] = {0., 0., 0.}, max[3] = {0., 0., 0.};
  for (i = 0; i < 3; i++)
  {
    if (fHistoDef[i])
    {
      nbins[i] = fHistoDef[i]->GetNBins();
      min[i] = fHistoDef[i]->GetMin();
      max[i] = fHistoDef[i]->GetMax();
    }
  }
    
  // define the kind of output according to the number of histoDefs
  if (fHistogram) delete fHistogram;
  if (!nbins[1] && !nbins[2]) {
    fHistogram = new TH1D(histoName, histoTitle, nbins[0], min[0], max[0]);
    fHistogram->SetXTitle(FcnName());
  }
  else if (nbins[1] > 0 && !nbins[2]) {
    fHistogram = new TH2D(histoName, histoTitle, nbins[0], min[0], max[0], nbins[1], min[1], max[1]);
    fHistogram->SetXTitle(FcnName());
    fHistogram->SetYTitle(BinName(fBinType[0]));
  }
  else {
    fHistogram = new TH3D(histoName, histoTitle, nbins[0], min[0], max[0], nbins[1], min[1], max[1], nbins[2], min[2], max[2]);
    fHistogram->SetXTitle(FcnName());
    fHistogram->SetYTitle(BinName(fBinType[0]));
    fHistogram->SetZTitle(BinName(fBinType[1]));
  }
  
  fHistogram->Sumw2();

  return fHistogram;
}

//________________________________________________________________________________________
Bool_t AliRsnFunction::Fill()
{
//
// Fill function histogram with values computed from given input object.
//
  AliDebug(AliLog::kDebug +2,"->");

  // checks coherence between fcn type and passed argument
  switch (fFcnType) {
    case kTrackPt:
    case kTrackEta:
      if (!CheckInput("TRACK")) return kFALSE;
      break;
    case kInvMass:
    case kInvMassMC:
    case kResolution:
    case kPairPt:
    case kPairEta:
      if (!CheckInput("PAIR")) return kFALSE;
      break;
    case kEventMult:
      if (!CheckInput("EVENT")) return kFALSE;
      break;
    default:
      AliError(Form("Input type %d not defined", (Int_t)fFcnType));
      return kFALSE;
  }

  // check presence of output histogram
  if (!fHistogram) {
    AliError("Histogram is not yet initialized");
    return kFALSE;
  }

  // compute value and stores into histogram
  Int_t    dim = CheckDim();
  Double_t mainValue, binValue[2];

  TH1D *h1 = dynamic_cast<TH1D*>(fHistogram);
  TH2D *h2 = dynamic_cast<TH2D*>(fHistogram);
  TH3D *h3 = dynamic_cast<TH3D*>(fHistogram);
  
  mainValue = Eval();

  switch (dim)
  {
    case 1:
      if (h1) h1->Fill(mainValue);
      break;
    case 2:
      binValue[0] = BinValue(fBinType[0]);
      if (h2) h2->Fill(mainValue, binValue[0]);
      break;
    case 3:
      binValue[0] = BinValue(fBinType[0]);
      binValue[1] = BinValue(fBinType[1]);
      if (h3) h3->Fill(mainValue, binValue[0], binValue[1]);
      break;
    default:
      AliError("Wrong number of dimensions in the histogram. Check HD initialization");
      return kFALSE;
  }

  AliDebug(AliLog::kDebug +2,"->");
  return kTRUE;
}

//________________________________________________________________________________________
Double_t AliRsnFunction::BinValue(EBinType binType)
{
//
// Computes the value for binning from the argument.
// For each kind of binning type, the object is expected
// to be of a given type, otherwise an error is raised.
//

  // checks coherence between bin type and passed argument
  switch (binType) {
    case kBinPairPt:
      if (!CheckInput("PAIR")) return 0.0;
      return fPair->GetPt();
    case kBinPairEta:
      if (!CheckInput("PAIR")) return 0.0;
      return fPair->GetEta();
    case kBinEventMult:
      if (!CheckInput("EVENT")) return 0.0;
      return fEvent->GetMultiplicity();
    default:
      AliError(Form("%s: Binning type not defined", GetName()));
      return 0.0;
  }
}

//________________________________________________________________________________________
Int_t AliRsnFunction::CheckDim()
{
//
// Checks number of dimensions.
// Makes sure that eventual binnings are coherent and well defined
//

  if (!fHistoDef[0]) return 0;
  if (fHistoDef[0] && !fHistoDef[1] && !fHistoDef[2]) return 1;
  if (fHistoDef[0] && fHistoDef[1] && !fHistoDef[2]) return 2;

  return 3;
}

//________________________________________________________________________________________
const char* AliRsnFunction::BinName(EBinType binType)
{
//
// Defines the name of binning
//

  switch (binType)
  {
    case kBinPairPt:
      return "PT";
      break;
    case kBinPairEta:
      return "ETA";
      break;
    case kBinEventMult:
      return "MULT";
      break;
    default:
      return "UNDEF";
  }
}