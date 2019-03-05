///
/// \file AliFemtoCorrFctn3DLCMSSym.cxx
///


#include "AliFemtoCorrFctn3DLCMSSym.h"
#include "AliFemtoPairCut.h"

#include <TH3F.h>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoCorrFctn3DLCMSSym);
  /// \endcond
#endif


AliFemtoCorrFctn3DLCMSSym
  ::AliFemtoCorrFctn3DLCMSSym(const char* title,
                              const int nbins,
                              const float QHi)
  : AliFemtoCorrFctn3DLCMSSym(title, nbins, QHi, true)
{}


AliFemtoCorrFctn3DLCMSSym::AliFemtoCorrFctn3DLCMSSym(const char* title,
                                                     const int nbins,
                                                     const float QHi,
                                                     bool enable_errs):
  AliFemtoCorrFctn()
  , fNumerator(nullptr)
  , fDenominator(nullptr)
  , fNumeratorW(nullptr)
  , fDenominatorW(nullptr)
  , fUseLCMS(1)
{
  // Basic constructor
  TString hist_title = TString::Format("%s; q_{out} (GeV); q_{side} (GeV); q_{long} (GeV)", title);

  // set up numerator
  fNumerator = new TH3F(TString("Num") + title,
                        hist_title,
                        nbins, -QHi, QHi,
                        nbins, -QHi, QHi,
                        nbins, -QHi, QHi);
  // set up denominator
  fDenominator = new TH3F(TString("Den") + title,
                          hist_title,
                          nbins, -QHi, QHi,
                          nbins, -QHi, QHi,
                          nbins, -QHi, QHi);

  // Weighted by qinv histos
  // set up numerator
  fNumeratorW = new TH3F(TString("NumWqinv") + title,
                         hist_title,
                         nbins, -QHi, QHi,
                         nbins, -QHi, QHi,
                         nbins, -QHi, QHi);

  // set up denominator
  fDenominatorW = new TH3F(TString("DenWqinv") + title,
                           hist_title,
                           nbins, -QHi, QHi,
                           nbins, -QHi, QHi,
                           nbins, -QHi, QHi);

  // to enable error bar calculation...
  if (enable_errs) {
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  fNumeratorW->Sumw2();
  fDenominatorW->Sumw2();
  }
}

AliFemtoCorrFctn3DLCMSSym::AliFemtoCorrFctn3DLCMSSym(const AliFemtoCorrFctn3DLCMSSym& aCorrFctn):
  AliFemtoCorrFctn(aCorrFctn)
  , fNumerator(new TH3F(*aCorrFctn.fNumerator))
  , fDenominator(new TH3F(*aCorrFctn.fDenominator))
  , fNumeratorW(new TH3F(*aCorrFctn.fNumeratorW))
  , fDenominatorW(new TH3F(*aCorrFctn.fDenominatorW))
  , fUseLCMS(aCorrFctn.fUseLCMS)
{
  // Copy constructor
}

//____________________________
AliFemtoCorrFctn3DLCMSSym::~AliFemtoCorrFctn3DLCMSSym()
{
  // Destructor
  delete fNumerator;
  delete fDenominator;
  delete fNumeratorW;
  delete fDenominatorW;
}
//_________________________
AliFemtoCorrFctn3DLCMSSym& AliFemtoCorrFctn3DLCMSSym::operator=(const AliFemtoCorrFctn3DLCMSSym& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  AliFemtoCorrFctn::operator=(aCorrFctn);

  *fNumerator = *aCorrFctn.fNumerator;
  *fDenominator = *aCorrFctn.fDenominator;
  *fNumeratorW = *aCorrFctn.fNumeratorW;
  *fDenominatorW = *aCorrFctn.fDenominatorW;

  fUseLCMS = aCorrFctn.fUseLCMS;

  return *this;
}

//_________________________
void AliFemtoCorrFctn3DLCMSSym::WriteOutHistos()
{
  // Write out all histograms to file
  fNumerator->Write();
  fDenominator->Write();
  fNumeratorW->Write();
  fDenominatorW->Write();
}
//______________________________
TList* AliFemtoCorrFctn3DLCMSSym::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumerator);
  tOutputList->Add(fDenominator);
  tOutputList->Add(fNumeratorW);
  tOutputList->Add(fDenominatorW);

  return tOutputList;
}

//_________________________
void AliFemtoCorrFctn3DLCMSSym::Finish()
{
}

//____________________________
AliFemtoString AliFemtoCorrFctn3DLCMSSym::Report()
{
  AliFemtoString report = "LCMS Frame Bertsch-Pratt 3D Correlation Function Report:\n";
  report += Form("Number of entries in numerator:\t%E\n", fNumerator->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n", fDenominator->GetEntries());

  if (fPairCut) {
    report += "Here is the PairCut specific to this CorrFctn\n";
    report += fPairCut->Report();
  } else {
    report += "No PairCut specific to this CorrFctn\n";
  }

  return report;
}
//____________________________
void AliFemtoCorrFctn3DLCMSSym::AddRealPair(AliFemtoPair* pair)
{
  // perform operations on real pairs
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  const Double_t qout = (fUseLCMS) ? pair->QOutCMS() : pair->QOutPf(),
                 qside = (fUseLCMS) ? pair->QSideCMS() : pair->QSidePf(),
                 qlong = (fUseLCMS) ? pair->QLongCMS() : pair->QLongPf();

  Int_t bin = fNumerator->FindBin(qout, qside, qlong);

  // avoid overflow bins
  if (!(fNumerator->IsBinOverflow(bin) or fNumerator->IsBinUnderflow(bin))) {
    fNumerator->Fill(qout, qside, qlong, 1.0);
    fNumeratorW->Fill(qout, qside, qlong, pair->QInv());
  }
}
//____________________________
void AliFemtoCorrFctn3DLCMSSym::AddMixedPair(AliFemtoPair* pair)
{
  // perform operations on mixed pairs
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  const Double_t qout = (fUseLCMS) ? pair->QOutCMS() : pair->QOutPf(),
                 qside = (fUseLCMS) ? pair->QSideCMS() : pair->QSidePf(),
                 qlong = (fUseLCMS) ? pair->QLongCMS() : pair->QLongPf();

  Int_t bin = fDenominator->FindBin(qout, qside, qlong);

  // avoid overflow bins
  if (!(fDenominator->IsBinOverflow(bin) or fDenominator->IsBinUnderflow(bin))) {
    fDenominator->Fill(qout, qside, qlong, 1.0);
    fDenominatorW->Fill(qout, qside, qlong, pair->QInv());
  }
}

void AliFemtoCorrFctn3DLCMSSym::SetUseLCMS(int aUseLCMS)
{
  fUseLCMS = aUseLCMS;
}

int  AliFemtoCorrFctn3DLCMSSym::GetUseLCMS()
{
  return fUseLCMS;
}
