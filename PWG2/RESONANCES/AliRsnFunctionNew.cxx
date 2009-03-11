//
// Class AliRsnFunctionNew
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
#include <TH2.h>
#include <TList.h>
#include <TString.h>

#include "AliLog.h"

#include "AliRsnMCInfo.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnPairDef.h"
#include "AliRsnCut.h"

#include "AliRsnFunctionNew.h"

ClassImp(AliRsnFunctionNew)

//________________________________________________________________________________________
AliRsnFunctionNew::AliRsnFunctionNew() :
  fFcnDef(0x0),
  fHistoTot(0x0)
{
  //
  // Constructor.
  // The histogram data member cannot be passed externally,
  // its initialization MUST be defined inside the Init() method,
  // which must be overridden in any derivate implementation.
  //

  Int_t i;
  for (i = 0 ; i < AliRsnFunctionDef::kBinMax; i++) fHistoBin[i] = 0x0;
}

//________________________________________________________________________________________
AliRsnFunctionNew::AliRsnFunctionNew(AliRsnFunctionDef *fd) :
  fFcnDef(fd),
  fHistoTot(0x0)
{
  //
  // Constructor.
  // The histogram data member cannot be passed externally,
  // its initialization MUST be defined inside the Init() method,
  // which must be overridden in any derivate implementation.
  //

  Int_t i;
  for (i = 0 ; i < AliRsnFunctionDef::kBinMax; i++) fHistoBin[i] = 0x0;
}

//________________________________________________________________________________________
AliRsnFunctionNew::AliRsnFunctionNew(const AliRsnFunctionNew &copy) :
  TObject(copy),
  fFcnDef(copy.fFcnDef),
  fHistoTot(0x0)
{
  //
  // Copy constructor.
  // Calls the function to define binning.
  //

  Int_t i;
  for (i = 0 ; i < AliRsnFunctionDef::kBinMax; i++) fHistoBin[i] = 0x0;
}
//________________________________________________________________________________________
const AliRsnFunctionNew& AliRsnFunctionNew::operator=(const AliRsnFunctionNew& /*copy*/)
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
void AliRsnFunctionNew::Clear(Option_t* /*option*/)
{
  //
  // Clear arrays and histogram.
  // For the sake of security, all pointers are also set explicitly to NULL.
  //

  Int_t i;
  for (i = 0 ; i < AliRsnFunctionDef::kBinMax; i++)
  {
    delete fHistoBin[i];
    fHistoBin[i] = 0x0;
  }

  delete fHistoTot;
  fHistoTot = 0x0;
}

//________________________________________________________________________________________
TList* AliRsnFunctionNew::Init(const char *histoName, const char *histoTitle)
{
  //
  // Initialization function.
  // By default, it initializes the owned histogram using the method
  // from AliRsnHistoDef class, giving the same name and title of this.
  // A user can override this behaviour, if necessary.
  // Before creating, the HistoDef is checked for proper initialization.
  //

  Clear();

  Int_t nbins = fFcnDef->GetNBinsX();
  Double_t min = fFcnDef->GetXmin(), max = fFcnDef->GetXmax();

  // list is created and named after the general
  // settings used for the contained histograms
  TList *histos = new TList;
  histos->SetName(Form("%s", fFcnDef->GetFcnName().Data()));

  // a general histogram is always added,
  // which overrides the binning and collects everything
  fHistoTot = new TH1D(histoName, histoTitle, nbins, min, max);
  fHistoTot->Sumw2();
  histos->AddLast(fHistoTot);

  // if requested a binning w.r. to some cut variable, histograms are added
  // for that in this part of the method (one per each bin)
  Char_t hName[255];
  Char_t hTitle[255];
  Int_t j, binMax = fFcnDef->GetMaxYBinUsed();
  if (binMax > 0)
  {
    for (j = 0 ; j < binMax; j++)
    {
      Int_t     nY    = fFcnDef->GetNBinsY(j);
      TArrayD  *binsY = fFcnDef->GetYbins(j);
      Double_t *arrY  = binsY->GetArray();
      sprintf(hName, "%s_%02d_%s", histoName, j, fFcnDef->GetBinningName(j).Data());
      sprintf(hTitle, "%s [%s]", histoTitle, fFcnDef->GetBinningName(j).Data());
      fHistoBin[j] = new TH2D(hName, hTitle, nbins, min, max, nY - 1, arrY);
      fHistoBin[j]->Sumw2();
      histos->AddLast(fHistoBin[j]);
    }
  }

  // returns the full list at the end
  return histos;
}

//________________________________________________________________________________________
void AliRsnFunctionNew::Init(const char *histoName, const char *histoTitle, TList *histos)
{
  //
  // Initialization function.
  // By default, it initializes the owned histogram using the method
  // from AliRsnHistoDef class, giving the same name and title of this.
  // A user can override this behaviour, if necessary.
  // Before creating, the HistoDef is checked for proper initialization.
  //

  Clear();

  Int_t nbins = fFcnDef->GetNBinsX();
  Double_t min = fFcnDef->GetXmin(), max = fFcnDef->GetXmax();

  // list is not created but its existence is checked
  if (!histos) return;

  // a general histogram is always added,
  // which overrides the binning and collects everything
  fHistoTot = new TH1D(histoName, histoTitle, nbins, min, max);
  fHistoTot->Sumw2();
  histos->AddLast(fHistoTot);

  // if requested a binning w.r. to some cut variable, histograms are added
  // for that in this part of the method (one per each bin)
  Char_t hName[255];
  Char_t hTitle[255];
  Int_t j, binMax = fFcnDef->GetMaxYBinUsed();
  if (binMax > 0)
  {
    for (j = 0 ; j < binMax; j++)
    {
      Int_t    nY    = fFcnDef->GetNBinsY(j);
      TArrayD *binsY = fFcnDef->GetYbins(j);
      sprintf(hName, "%s_%02d_%s", histoName, j, fFcnDef->GetBinningName(j).Data());
      sprintf(hTitle, "%s [%s]", histoTitle, fFcnDef->GetBinningName(j).Data());
      fHistoBin[j] = new TH2D(hName, hTitle, nbins, min, max, nY, binsY->GetArray());
      fHistoBin[j]->Sumw2();
      histos->AddLast(fHistoBin[j]);
    }
  }
}

//________________________________________________________________________________________
Bool_t AliRsnFunctionNew::Fill(AliRsnPairParticle *pair, AliRsnPairDef *ref)
{
  //
  // Fillse the histogram with data contained in a defined pair.
  // This method must be overidden by an appropriate definition in each inheriting class.
  //

  Double_t valueY, valueX = fFcnDef->EvalX(pair, ref);

  // fill global histogram
  fHistoTot->Fill(valueX);

  // if bins are allocated, fill them
  Int_t j, binMax = fFcnDef->GetMaxYBinUsed();
  if (binMax > 0)
  {
    for (j = 0 ; j < binMax; j++)
    {
      valueY = fFcnDef->EvalY(j, pair);
      if (fHistoBin[j]) fHistoBin[j]->Fill(valueX, valueY);
    }
  }

  return kTRUE;
}
