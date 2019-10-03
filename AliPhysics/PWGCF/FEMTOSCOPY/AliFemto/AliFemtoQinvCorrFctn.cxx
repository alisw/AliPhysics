///
/// \file AliFemtoQinvCorrFctn.cxx
///

#include "AliFemtoQinvCorrFctn.h"
// #include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoQinvCorrFctn);
  /// \endcond
#endif

//____________________________
AliFemtoQinvCorrFctn::AliFemtoQinvCorrFctn(const char* title,
                                           const int& nbins,
                                           const float& QinvLo,
                                           const float& QinvHi):
  AliFemtoCorrFctn()
  , fNumerator(nullptr)
  , fDenominator(nullptr)
  , fRatio(nullptr)
  , fkTMonitor(nullptr)
  , fDetaDphiscal(kFALSE)
  , fPairKinematics(kFALSE)
  , fRaddedps(1.2)
  , fNumDEtaDPhiS(nullptr)
  , fDenDEtaDPhiS(nullptr)
  , PairReader(nullptr)
{
  // set up numerator
  //  title = "Num Qinv (MeV/c)";
  char tTitNum[101] = "Num";
  strncat(tTitNum,title, 100);
  fNumerator = new TH1D(tTitNum,title,nbins,QinvLo,QinvHi);
  // set up denominator
  //title = "Den Qinv (MeV/c)";
  char tTitDen[101] = "Den";
  strncat(tTitDen,title, 100);
  fDenominator = new TH1D(tTitDen,title,nbins,QinvLo,QinvHi);
  // set up ratio
  //title = "Ratio Qinv (MeV/c)";
  char tTitRat[101] = "Rat";
  strncat(tTitRat,title, 100);
  fRatio = new TH1D(tTitRat,title,nbins,QinvLo,QinvHi);

  char tTitkT[101] = "kTDep";
  strncat(tTitkT,title, 100);
  fkTMonitor = new TH1D(tTitkT,title,250,0.0,5.0);

  char tTitNumDeDp[101] = "NumDEtaDPhiS";
  strncat(tTitNumDeDp,title, 100);
  fNumDEtaDPhiS = new TH2D(tTitNumDeDp,title,500,-0.2*TMath::Pi(),0.2*TMath::Pi(),500,-0.5,0.5);

  char tTitDenDeDp[101] = "DenDEtaDPhiS";
  strncat(tTitDenDeDp,title, 100);
  fDenDEtaDPhiS = new TH2D(tTitDenDeDp,title,500,-0.2*TMath::Pi(),0.2*TMath::Pi(),500,-0.5,0.5);

  char tTitPair[101] = "Pair";
  strncat(tTitPair,title, 100);
  PairReader = new TNtuple(tTitPair,title,  "px1:py1:pz1:e1:px2:py2:pz2:e2");


  // this next bit is unfortunately needed so that we can have many histos of same "title"
  // it is neccessary if we typedef TH1D to TH1d (which we do)
  //fNumerator->SetDirectory(0);
  //fDenominator->SetDirectory(0);
  //fRatio->SetDirectory(0);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  fRatio->Sumw2();
  fkTMonitor->Sumw2();

  fNumDEtaDPhiS->Sumw2();
  fDenDEtaDPhiS->Sumw2();
}

//____________________________
AliFemtoQinvCorrFctn::AliFemtoQinvCorrFctn(const AliFemtoQinvCorrFctn& aCorrFctn):
  AliFemtoCorrFctn(aCorrFctn)
  , fNumerator(nullptr)
  , fDenominator(nullptr)
  , fRatio(nullptr)
  , fkTMonitor(nullptr)
  , fDetaDphiscal(aCorrFctn.fDetaDphiscal)
  , fPairKinematics(aCorrFctn.fPairKinematics)
  , fRaddedps(aCorrFctn.fRaddedps)
  , fNumDEtaDPhiS(nullptr)
  , fDenDEtaDPhiS(nullptr)
  , PairReader(nullptr)
{
  // copy constructor
  fNumerator = new TH1D(*aCorrFctn.fNumerator);
  fDenominator = new TH1D(*aCorrFctn.fDenominator);
  fRatio = new TH1D(*aCorrFctn.fRatio);
  fkTMonitor = new TH1D(*aCorrFctn.fkTMonitor);

  fNumDEtaDPhiS = new TH2D(*aCorrFctn.fNumDEtaDPhiS);
  fDenDEtaDPhiS = new TH2D(*aCorrFctn.fDenDEtaDPhiS);

  if (aCorrFctn.PairReader) {
    PairReader = (TNtuple*)aCorrFctn.PairReader->Clone();
  }
}
//____________________________
AliFemtoQinvCorrFctn::~AliFemtoQinvCorrFctn()
{
  // destructor
  delete fNumerator;
  delete fDenominator;
  delete fRatio;
  delete fkTMonitor;
  delete fNumDEtaDPhiS;
  delete fDenDEtaDPhiS;
  delete PairReader;
}
//_________________________
AliFemtoQinvCorrFctn& AliFemtoQinvCorrFctn::operator=(const AliFemtoQinvCorrFctn& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn) {
    return *this;
  }

  AliFemtoCorrFctn::operator=(aCorrFctn);

  *fNumerator = *aCorrFctn.fNumerator;
  *fDenominator = *aCorrFctn.fDenominator;
  *fRatio = *aCorrFctn.fRatio;
  *fkTMonitor = *aCorrFctn.fkTMonitor;

  *fNumDEtaDPhiS = *aCorrFctn.fNumDEtaDPhiS;
  *fDenDEtaDPhiS = *aCorrFctn.fDenDEtaDPhiS;

  fDetaDphiscal = aCorrFctn.fDetaDphiscal;
  fRaddedps = aCorrFctn.fRaddedps;

  fPairKinematics = aCorrFctn.fPairKinematics;

  delete PairReader;
  PairReader = (TNtuple*)aCorrFctn.PairReader->Clone();

  return *this;
}

void AliFemtoQinvCorrFctn::Finish()
{
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  fNumerator->Draw();
  //fDenominator->Draw();
  //fRatio->Draw();
  fRatio->Divide(fNumerator,fDenominator,1.0,1.0);
}

AliFemtoString AliFemtoQinvCorrFctn::Report()
{
  // construct report
  AliFemtoString report = "Qinv Correlation Function Report:\n";
  report += Form("Number of entries in numerator:\t%E\n", fNumerator->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n", fDenominator->GetEntries());
  report += Form("Number of entries in ratio:\t%E\n", fRatio->GetEntries());

  return report;
}


// Function used by both AddRealPair & AddMixedPair to do the appropriate calculations
// If the deta_dphi_hist argument is null, the (Δη, Δϕ*) calculation will be skipped,
// otherwise the results are stored in the histogram at which it points.
//
static
void AddPair(const AliFemtoPair &pair, TH1 &qinv_hist, TH2 *deta_dphi_hist, double rad)
{
  double qinv = fabs(pair.QInv());
  qinv_hist.Fill(qinv);

  if (deta_dphi_hist) {

    const AliAODInputHandler *aodH = static_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!aodH) {
      std::cerr << "E-AliFemtoQCorrFctn: Could not get AODInputHandler\n";
    }

    double magfield = aodH->GetEvent()->GetMagneticField();
    Int_t magsign = magfield == 0.0 ? 0 : std::copysign(magfield, 1);

    const AliFemtoTrack
      &track1 = *pair.Track1()->Track(),
      &track2 = *pair.Track2()->Track();

    const AliFemtoThreeVector
      &p1 = track1.P(),
      &p2 = track2.P();

    const double
      delta_eta = p2.PseudoRapidity() - p1.PseudoRapidity(),
      delta_phi = p2.Phi() - p1.Phi(),

      afsi0b = 0.07510020733 * magsign * rad * track1.Charge() / track1.Pt(),
      afsi1b = 0.07510020733 * magsign * rad * track2.Charge() / track2.Pt(),

      delta_phistar = TVector2::Phi_mpi_pi(delta_phi + TMath::ASin(afsi1b) - TMath::ASin(afsi0b));

    deta_dphi_hist->Fill(delta_phistar, delta_eta);
  }
}

void AliFemtoQinvCorrFctn::AddRealPair(AliFemtoPair* pair)
{
  // add true pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  // if null, the (Δη, Δϕ*) calculation will be skipped
  AddPair(*pair, *fNumerator, fDetaDphiscal ? fNumDEtaDPhiS : nullptr, fRaddedps);

  fkTMonitor->Fill(pair->KT());
}

//____________________________
void AliFemtoQinvCorrFctn::AddMixedPair(AliFemtoPair* pair)
{
  // add mixed (background) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  AddPair(*pair, *fDenominator, fDetaDphiscal ? fDenDEtaDPhiS : nullptr, fRaddedps);

  if (fPairKinematics) {
    const AliFemtoParticle &fTrack1 = *pair->Track1(),
                           &fTrack2 = *pair->Track2();

    const auto &p1 = fTrack1.FourMomentum(),
               &p2 = fTrack2.FourMomentum();

    PairReader->Fill(p1.x(), p1.y(), p1.z(), p1.e(),
		     p2.x(), p2.y(), p2.z(), p2.e());
  }
}

void AliFemtoQinvCorrFctn::Write()
{
  // Write out neccessary objects
  fNumerator->Write();
  fDenominator->Write();
  fkTMonitor->Write();
  if (fDetaDphiscal) {
    fNumDEtaDPhiS->Write();
    fDenDEtaDPhiS->Write();
  }
  if (fPairKinematics) {
    PairReader->Write();
  }
}

TList* AliFemtoQinvCorrFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumerator);
  tOutputList->Add(fDenominator);
  tOutputList->Add(fkTMonitor);
  if (fDetaDphiscal) {
    tOutputList->Add(fNumDEtaDPhiS);
    tOutputList->Add(fDenDEtaDPhiS);
  }
  if (fPairKinematics) {
    tOutputList->Add(PairReader);
  }
  return tOutputList;
}
