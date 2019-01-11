///
/// \file AliFemto/AliFemtoCorrFctn3DSpherical.cxx
///

#include "AliFemtoCorrFctn3DSpherical.h"
#include <TMath.h>
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoCorrFctn3DSpherical);
  /// \endcond
#endif

//____________________________
AliFemtoCorrFctn3DSpherical::AliFemtoCorrFctn3DSpherical(const char* title,
                                                         const int nqbins,
                                                         const float QLo,
                                                         const float QHi,
                                                         const int nphibins,
                                                         const int ncthetabins)

  : AliFemtoCorrFctn()
  , fNumerator(nullptr)
  , fDenominator(nullptr)
{
  // set up numerator
  TString tTitNum = TString("Num") + title;
  fNumerator = new TH3D(tTitNum,title,nqbins,QLo,QHi,nphibins,-TMath::Pi(),TMath::Pi(),ncthetabins,-1.0,1.0);
  // set up denominator
  TString tTitDen = TString("Den") + title;
  fDenominator = new TH3D(tTitDen,title,nqbins,QLo,QHi,nphibins,-TMath::Pi(),TMath::Pi(),ncthetabins,-1.0,1.0);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();
}

AliFemtoCorrFctn3DSpherical::AliFemtoCorrFctn3DSpherical(const AliFemtoCorrFctn3DSpherical& aCorrFctn)
  : AliFemtoCorrFctn(aCorrFctn)
  , fNumerator(nullptr)
  , fDenominator(nullptr)
{
  // Copy constructor
  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  fDenominator = new TH3D(*aCorrFctn.fDenominator);
}
//____________________________
AliFemtoCorrFctn3DSpherical::~AliFemtoCorrFctn3DSpherical()
{
  // Destructor
  delete fNumerator;
  delete fDenominator;
}
//_________________________
AliFemtoCorrFctn3DSpherical& AliFemtoCorrFctn3DSpherical::operator=(const AliFemtoCorrFctn3DSpherical& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  AliFemtoCorrFctn::operator=(aCorrFctn);

  if (fNumerator) delete fNumerator;
  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  if (fDenominator) delete fDenominator;
  fDenominator = new TH3D(*aCorrFctn.fDenominator);

  return *this;
}

//_________________________
void AliFemtoCorrFctn3DSpherical::WriteOutHistos()
{
  // Write out all histograms to file
  fNumerator->Write();
  fDenominator->Write();
}
//______________________________
TList* AliFemtoCorrFctn3DSpherical::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumerator);
  tOutputList->Add(fDenominator);

  return tOutputList;
}

//_________________________
void AliFemtoCorrFctn3DSpherical::Finish()
{
  // here is where we should normalize, fit, etc...
}

//____________________________
AliFemtoString AliFemtoCorrFctn3DSpherical::Report()
{
  // Construct the report
  AliFemtoString report = "PRF Frame Spherical 3D Correlation Function Report:\n";
  report += Form("Number of entries in numerator:\t%E\n", fNumerator->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n", fDenominator->GetEntries());

  if (fPairCut) {
    report += "Here is the PairCut specific to this CorrFctn\n" + fPairCut->Report();
  } else {
    report += "No PairCut specific to this CorrFctn\n";
  }

  return report;
}
//____________________________
void AliFemtoCorrFctn3DSpherical::AddRealPair(AliFemtoPair *pair)
{
  // perform operations on real pairs
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  double tKO = pair->KOut();
  double tKS = pair->KSide();
  double tKL = pair->KLong();

  double tKR = sqrt(tKO*tKO + tKS*tKS + tKL*tKL);
  double tKC;
  if ( fabs(tKR) < 1e-10 ) tKC = 0.0;
  else tKC=tKL/tKR;
  double tKP=atan2(tKS,tKO);

  fNumerator->Fill(tKR,tKP,tKC);
}
//____________________________
void AliFemtoCorrFctn3DSpherical::AddMixedPair(AliFemtoPair *pair)
{
  // perform operations on mixed pairs
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  double tKO = pair->KOut();
  double tKS = pair->KSide();
  double tKL = pair->KLong();

  double tKR = sqrt(tKO*tKO + tKS*tKS + tKL*tKL);
  double tKC;
  if ( fabs(tKR) < 1e-10 ) tKC = 0.0;
  else tKC=tKL/tKR;
  double tKP=atan2(tKS,tKO);

  fDenominator->Fill(tKR,tKP,tKC);
}
