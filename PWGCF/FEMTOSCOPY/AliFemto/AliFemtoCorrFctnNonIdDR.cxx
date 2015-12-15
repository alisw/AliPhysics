///
/// \file AliFemtoCorrFctnNonIdDR.cxx
///


#include "AliFemtoCorrFctnNonIdDR.h"
#include <TH1D.h>
//#include "AliFemtoHisto.h"
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoCorrFctnNonIdDR);
  /// \endcond
#endif

//____________________________
AliFemtoCorrFctnNonIdDR::AliFemtoCorrFctnNonIdDR(const char* title,
                                                 const int nbins,
                                                 const float QinvLo,
                                                 const float QinvHi):
  AliFemtoCorrFctn()
  , fNumOutP(NULL)
  , fNumOutN(NULL)
  , fNumSideP(NULL)
  , fNumSideN(NULL)
  , fNumLongP(NULL)
  , fNumLongN(NULL)
  , fDenOutP(NULL)
  , fDenOutN(NULL)
  , fDenSideP(NULL)
  , fDenSideN(NULL)
  , fDenLongP(NULL)
  , fDenLongN(NULL)
  , fkTMonitor(NULL)
{
  // Default constructor
  // set up numerators
  fNumOutP = new TH1D(TString("NumOutP") + title, title, nbins, QinvLo, QinvHi);
  fNumOutN = new TH1D(TString("NumOutN") + title, title, nbins, QinvLo, QinvHi);
  fNumSideP = new TH1D(TString("NumSideP") + title, title, nbins, QinvLo, QinvHi);
  fNumSideN = new TH1D(TString("NumSideN") + title, title, nbins, QinvLo, QinvHi);
  fNumLongP = new TH1D(TString("NumLongP") + title, title, nbins, QinvLo, QinvHi);
  fNumLongN = new TH1D(TString("NumLongN") + title, title, nbins, QinvLo, QinvHi);

  // set up denominators
  fDenOutP = new TH1D(TString("DenOutP") + title, title, nbins, QinvLo, QinvHi);
  fDenOutN = new TH1D(TString("DenOutN") + title, title, nbins, QinvLo, QinvHi);
  fDenSideP = new TH1D(TString("DenSideP") + title, title, nbins, QinvLo, QinvHi);
  fDenSideN = new TH1D(TString("DenSideN") + title, title, nbins, QinvLo, QinvHi);
  fDenLongP = new TH1D(TString("DenLongP") + title, title, nbins, QinvLo, QinvHi);
  fDenLongN = new TH1D(TString("DenLongN") + title, title, nbins, QinvLo, QinvHi);

  fkTMonitor = new TH1D(TString("kTDep") + title, title, 250, 0.0, 5.0);

  // to enable error bar calculation...
  fNumOutP->Sumw2();
  fNumOutN->Sumw2();
  fNumSideP->Sumw2();
  fNumSideN->Sumw2();
  fNumLongP->Sumw2();
  fNumLongN->Sumw2();
  fDenOutP->Sumw2();
  fDenOutN->Sumw2();
  fDenSideP->Sumw2();
  fDenSideN->Sumw2();
  fDenLongP->Sumw2();
  fDenLongN->Sumw2();

  fkTMonitor->Sumw2();
}

//____________________________
AliFemtoCorrFctnNonIdDR::AliFemtoCorrFctnNonIdDR(const AliFemtoCorrFctnNonIdDR& aCorrFctn):
  AliFemtoCorrFctn(aCorrFctn)
  , fNumOutP(new TH1D(*aCorrFctn.fNumOutP))
  , fNumOutN(new TH1D(*aCorrFctn.fNumOutN))
  , fNumSideP(new TH1D(*aCorrFctn.fNumSideP))
  , fNumSideN(new TH1D(*aCorrFctn.fNumSideN))
  , fNumLongP(new TH1D(*aCorrFctn.fNumLongP))
  , fNumLongN(new TH1D(*aCorrFctn.fNumLongN))
  , fDenOutP(new TH1D(*aCorrFctn.fDenOutP))
  , fDenOutN(new TH1D(*aCorrFctn.fDenOutN))
  , fDenSideP(new TH1D(*aCorrFctn.fDenSideP))
  , fDenSideN(new TH1D(*aCorrFctn.fDenSideN))
  , fDenLongP(new TH1D(*aCorrFctn.fDenLongP))
  , fDenLongN(new TH1D(*aCorrFctn.fDenLongN))
  , fkTMonitor(new TH1D(*aCorrFctn.fkTMonitor))
{
  // copy constructor
  fNumOutP->Sumw2();
  fNumOutN->Sumw2();
  fNumSideP->Sumw2();
  fNumSideN->Sumw2();
  fNumLongP->Sumw2();
  fNumLongN->Sumw2();
  fDenOutP->Sumw2();
  fDenOutN->Sumw2();
  fDenSideP->Sumw2();
  fDenSideN->Sumw2();
  fDenLongP->Sumw2();
  fDenLongN->Sumw2();

  fkTMonitor->Sumw2();
}

//____________________________
AliFemtoCorrFctnNonIdDR::~AliFemtoCorrFctnNonIdDR()
{
  delete fNumOutP;
  delete fNumOutN;
  delete fNumSideP;
  delete fNumSideN;
  delete fNumLongP;
  delete fNumLongN;
  delete fDenOutP;
  delete fDenOutN;
  delete fDenSideP;
  delete fDenSideN;
  delete fDenLongP;
  delete fDenLongN;
  delete fkTMonitor;
}

//_________________________
AliFemtoCorrFctnNonIdDR& AliFemtoCorrFctnNonIdDR::operator=(const AliFemtoCorrFctnNonIdDR& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn) {
    return *this;
  }

  AliFemtoCorrFctn::operator=(aCorrFctn);

  // delete all histograms
  delete fNumOutP;
  delete fNumOutN;
  delete fNumSideP;
  delete fNumSideN;
  delete fNumLongP;
  delete fNumLongN;
  delete fDenOutP;
  delete fDenOutN;
  delete fDenSideP;
  delete fDenSideN;
  delete fDenLongP;
  delete fDenLongN;
  delete fkTMonitor;

  fNumOutP = new TH1D(*aCorrFctn.fNumOutP);
  fNumOutN = new TH1D(*aCorrFctn.fNumOutN);
  fNumSideP = new TH1D(*aCorrFctn.fNumSideP);
  fNumSideN = new TH1D(*aCorrFctn.fNumSideN);
  fNumLongP = new TH1D(*aCorrFctn.fNumLongP);
  fNumLongN = new TH1D(*aCorrFctn.fNumLongN);

  fDenOutP = new TH1D(*aCorrFctn.fDenOutP);
  fDenOutN = new TH1D(*aCorrFctn.fDenOutN);
  fDenSideP = new TH1D(*aCorrFctn.fDenSideP);
  fDenSideN = new TH1D(*aCorrFctn.fDenSideN);
  fDenLongP = new TH1D(*aCorrFctn.fDenLongP);
  fDenLongN = new TH1D(*aCorrFctn.fDenLongN);

  fkTMonitor = new TH1D(*aCorrFctn.fkTMonitor);

  fNumOutP->Sumw2();
  fNumOutN->Sumw2();
  fNumSideP->Sumw2();
  fNumSideN->Sumw2();
  fNumLongP->Sumw2();
  fNumLongN->Sumw2();
  fDenOutP->Sumw2();
  fDenOutN->Sumw2();
  fDenSideP->Sumw2();
  fDenSideN->Sumw2();
  fDenLongP->Sumw2();
  fDenLongN->Sumw2();

  fkTMonitor->Sumw2();

  return *this;
}

//_________________________
void AliFemtoCorrFctnNonIdDR::Finish()
{
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  fNumerator->Draw();
  //fDenominator->Draw();
  //fRatio->Draw();
  //  fRatio->Divide(fNumerator,fDenominator,1.0,1.0);

}

//____________________________
AliFemtoString AliFemtoCorrFctnNonIdDR::Report()
{
  // construct report
  TString report = "Non-identical particles Correlation Function Report:\n";
  report += TString::Format("Number of entries in numerators:\t%E\n", fNumOutP->GetEntries() + fNumOutN->GetEntries());
  report += TString::Format("Number of entries in denominators:\t%E\n", fDenOutP->GetEntries() + fDenOutN->GetEntries());

  //  report += mCoulombWeight->Report();
  return AliFemtoString(report);
}
//____________________________
void AliFemtoCorrFctnNonIdDR::AddRealPair(AliFemtoPair* pair)
{ // add true pair

  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  double tKStar = pair->KStar();
  if (pair->KOut()>0.0)
    fNumOutP->Fill(tKStar);
  else
    fNumOutN->Fill(tKStar);

  if (pair->KSide()>0.0)
    fNumSideP->Fill(tKStar);
  else
    fNumSideN->Fill(tKStar);

  if (pair->KLong()>0.0)
    fNumLongP->Fill(tKStar);
  else
    fNumLongN->Fill(tKStar);

  fkTMonitor->Fill(pair->KT());
}
//____________________________
void AliFemtoCorrFctnNonIdDR::AddMixedPair(AliFemtoPair* pair)
{ // add mixed (background) pair

  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  double tKStar = pair->KStar();
  if (pair->KOut() > 0.0)
    fDenOutP->Fill(tKStar);
  else
    fDenOutN->Fill(tKStar);

  if (pair->KSide() > 0.0)
    fDenSideP->Fill(tKStar);
  else
    fDenSideN->Fill(tKStar);

  if (pair->KLong() > 0.0)
    fDenLongP->Fill(tKStar);
  else
    fDenLongN->Fill(tKStar);
}
//____________________________
void AliFemtoCorrFctnNonIdDR::Write()
{
  fNumOutP->Write();
  fNumOutN->Write();
  fNumSideP->Write();
  fNumSideN->Write();
  fNumLongP->Write();
  fNumLongN->Write();
  fDenOutP->Write();
  fDenOutN->Write();
  fDenSideP->Write();
  fDenSideN->Write();
  fDenLongP->Write();
  fDenLongN->Write();
  fkTMonitor->Write();
}

TList* AliFemtoCorrFctnNonIdDR::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumOutP);
  tOutputList->Add(fNumOutN);
  tOutputList->Add(fNumSideP);
  tOutputList->Add(fNumSideN);
  tOutputList->Add(fNumLongP);
  tOutputList->Add(fNumLongN);
  tOutputList->Add(fDenOutP);
  tOutputList->Add(fDenOutN);
  tOutputList->Add(fDenSideP);
  tOutputList->Add(fDenSideN);
  tOutputList->Add(fDenLongP);
  tOutputList->Add(fDenLongN);
  tOutputList->Add(fkTMonitor);

  return tOutputList;
}
