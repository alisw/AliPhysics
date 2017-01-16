///
/// \file AliFemtoCorrFctnKStar.cxx
///

#include "AliFemtoCorrFctnKStar.h"


const float DFLT_KStarMin = 0.0, DFLT_KStarMax = 1.0,
            DFLT_kTMin    = 0.0, DFLT_kTMax    = 5.0,
            DFLT_mTMin    = 0.0, DFLT_mTMax    = 5.0;

const int DFLT_NbinsKStar = 200,
          DFLT_NbinskT    = 250,
          DFLT_NbinsmT    = 250;

const float DFLT_KStarOutMin  = -1.0, DFLT_KStarOutMax  = 1.0,
            DFLT_KStarSideMin = -1.0, DFLT_KStarSideMax = 1.0,
            DFLT_KStarLongMin = -1.0, DFLT_KStarLongMax = 1.0;

const int DFLT_NbinsKStarOut  = 400,
          DFLT_NbinsKStarSide = 400,
          DFLT_NbinsKStarLong = 400;


const float DFLT_DEtaMin  = -1.0,             DFLT_DEtaMax  = 1.0,
            DFLT_DPhiSMin = -0.5*TMath::Pi(), DFLT_DPhiSMax = 0.5*TMath::Pi();

const int DFLT_NbinsDEta = 1000,
          DFLT_NbinsDPhiS = 1000;



//____________________________
AliFemtoCorrFctnKStar::AliFemtoCorrFctnKStar():
  AliFemtoCorrFctnKStar("CorrFctnKStar", DFLT_NbinsKStar, DFLT_KStarMin, DFLT_KStarMax)
{
// no-op
}


//____________________________
AliFemtoCorrFctnKStar::AliFemtoCorrFctnKStar(const char* title,
                                             const int nbins,
                                             const float KStarLo,
                                             const float KStarHi):
  AliFemtoCorrFctn(),
  fTitle(title),
  fNbinsKStar(nbins),
  fKStarLow(KStarLo),
  fKStarHigh(KStarHi),
  fNumerator(nullptr),
  fDenominator(nullptr),
  fRatio(nullptr),
  fkTMonitor(nullptr),
  fDetaDphiscal(kFALSE),
  fPairKinematics(kFALSE),
  fRaddedps(1.2),
  fNumDEtaDPhiS(nullptr),
  fDenDEtaDPhiS(nullptr),
  fPairKStar(nullptr),

  fBuildkTBinned(kFALSE),
  fNumerator_kT(nullptr),
  fDenominator_kT(nullptr),

  fBuildmTBinned(kFALSE),
  fNumerator_mT(nullptr),
  fDenominator_mT(nullptr),
  fNumeratorv2_mT(nullptr),
  fDenominatorv2_mT(nullptr),

  fBuild3d(kFALSE),
  fNumerator3d(nullptr),
  fDenominator3d(nullptr)
{
  fNumerator = new TH1D(TString::Format("Num%s", fTitle.Data()),
                        "KStar - Numerator; k*(GeV/c);",
                        nbins, KStarLo, KStarHi);
  fDenominator = new TH1D(TString::Format("Den%s", fTitle.Data()),
                          "KStar - Denominator; k*(GeV/c);",
                          nbins, KStarLo, KStarHi);
  fRatio = new TH1D(TString::Format("Rat%s", fTitle.Data()),
                    "KStar - Ratio; k*(GeV/c);",
                    nbins, KStarLo, KStarHi);
  fkTMonitor = new TH1D(TString::Format("kTDep%s", fTitle.Data()),
                        "kT Dependence; kT(GeV/c)",
                        DFLT_NbinskT, DFLT_kTMin, DFLT_kTMax);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  fRatio->Sumw2();
  fkTMonitor->Sumw2();
}

//____________________________
AliFemtoCorrFctnKStar::AliFemtoCorrFctnKStar(const AliFemtoCorrFctnKStar& aCorrFctn):
  AliFemtoCorrFctn(aCorrFctn),
  fTitle(aCorrFctn.fTitle),
  fNbinsKStar(aCorrFctn.fNbinsKStar),

  fNumerator(aCorrFctn.fNumerator ? new TH1D(*aCorrFctn.fNumerator) : nullptr),
  fDenominator(aCorrFctn.fDenominator ? new TH1D(*aCorrFctn.fDenominator) : nullptr),

  fKStarLow(aCorrFctn.fKStarLow),
  fKStarHigh(aCorrFctn.fKStarHigh),
  fDetaDphiscal(aCorrFctn.fDetaDphiscal),
  fPairKinematics(aCorrFctn.fPairKinematics),
  fRaddedps(aCorrFctn.fRaddedps),
  fBuildkTBinned(aCorrFctn.fBuildkTBinned),
  fBuildmTBinned(aCorrFctn.fBuildmTBinned),
  fBuild3d(aCorrFctn.fBuild3d)
{
  // copy constructor
  fRatio = (aCorrFctn.fRatio) ? new TH1D(*aCorrFctn.fRatio) : nullptr;
  fkTMonitor =(aCorrFctn.fkTMonitor) ? static_cast<TH1D*>(aCorrFctn.fkTMonitor->Clone()) : nullptr;

  if(aCorrFctn.fNumDEtaDPhiS) fNumDEtaDPhiS = new TH2D(*aCorrFctn.fNumDEtaDPhiS);
    else fNumDEtaDPhiS = 0;

  if(aCorrFctn.fDenDEtaDPhiS) fDenDEtaDPhiS = new TH2D(*aCorrFctn.fDenDEtaDPhiS);
    else fDenDEtaDPhiS = 0;

  if(aCorrFctn.fPairKStar) fPairKStar = (TNtuple*)aCorrFctn.fPairKStar->CloneTree();
    else fPairKStar = 0;

  if(aCorrFctn.fNumerator_kT) fNumerator_kT = new TH2F(*aCorrFctn.fNumerator_kT);
    else fNumerator_kT = 0;
  if(aCorrFctn.fDenominator_kT) fDenominator_kT = new TH2F(*aCorrFctn.fDenominator_kT);
    else fDenominator_kT = 0;

  if(aCorrFctn.fNumerator_mT) fNumerator_mT = new TH2F(*aCorrFctn.fNumerator_mT);
    else fNumerator_mT = 0;
  if(aCorrFctn.fDenominator_mT) fDenominator_mT = new TH2F(*aCorrFctn.fDenominator_mT);
    else fDenominator_mT = 0;

  if(aCorrFctn.fNumeratorv2_mT) fNumeratorv2_mT = new TH2F(*aCorrFctn.fNumeratorv2_mT);
    else fNumeratorv2_mT = 0;
  if(aCorrFctn.fDenominatorv2_mT) fDenominatorv2_mT = new TH2F(*aCorrFctn.fDenominatorv2_mT);
    else fDenominatorv2_mT = 0;

  if(aCorrFctn.fNumerator3d) fNumerator3d = new TH3F(*aCorrFctn.fNumerator3d);
    else fNumerator3d = 0;
  if(aCorrFctn.fDenominator3d) fDenominator3d = new TH3F(*aCorrFctn.fDenominator3d);
    else fDenominator3d = 0;

}
//____________________________
AliFemtoCorrFctnKStar::~AliFemtoCorrFctnKStar()
{
  // destructor
  delete fNumerator;
  delete fDenominator;
  delete fRatio;
  delete fkTMonitor;
  delete fNumDEtaDPhiS;
  delete fDenDEtaDPhiS;
  delete fPairKStar;
  delete fNumerator_kT;
  delete fDenominator_kT;
  delete fNumerator_mT;
  delete fDenominator_mT;
  delete fNumeratorv2_mT;
  delete fDenominatorv2_mT;
  delete fNumerator3d;
  delete fDenominator3d;
}
//_________________________
AliFemtoCorrFctnKStar& AliFemtoCorrFctnKStar::operator=(const AliFemtoCorrFctnKStar& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn) {
    return *this;
  }

  AliFemtoCorrFctn::operator=(aCorrFctn);

  fNbinsKStar = aCorrFctn.fNbinsKStar;
  fKStarLow = aCorrFctn.fKStarLow;
  fKStarHigh = aCorrFctn.fKStarHigh;
  fDetaDphiscal = aCorrFctn.fDetaDphiscal;
  fRaddedps = aCorrFctn.fRaddedps;
  fPairKinematics = aCorrFctn.fPairKinematics;
  fBuildkTBinned = aCorrFctn.fBuildkTBinned;
  fBuildmTBinned = aCorrFctn.fBuildmTBinned;
  fBuild3d = aCorrFctn.fBuild3d;

  fTitle = aCorrFctn.fTitle;

  if(fNumerator) delete fNumerator;
    fNumerator = new TH1D(*aCorrFctn.fNumerator);
  if(fDenominator) delete fDenominator;
    fDenominator = new TH1D(*aCorrFctn.fDenominator);
  if(fRatio) delete fRatio;
    fRatio = new TH1D(*aCorrFctn.fRatio);
  if(fkTMonitor) delete fkTMonitor;
    fkTMonitor = new TH1D(*aCorrFctn.fkTMonitor);

  if(fNumDEtaDPhiS) delete fNumDEtaDPhiS;
    fNumDEtaDPhiS = new TH2D(*aCorrFctn.fNumDEtaDPhiS);
  if(fDenDEtaDPhiS) delete fDenDEtaDPhiS;
    fDenDEtaDPhiS = new TH2D(*aCorrFctn.fDenDEtaDPhiS);

  if(aCorrFctn.fPairKStar) delete fPairKStar;
    fPairKStar = (TNtuple*)aCorrFctn.fPairKStar->CloneTree();

  if(fNumerator_kT) delete fNumerator_kT;
    fNumerator_kT = new TH2F(*aCorrFctn.fNumerator_kT);
  if(fDenominator_kT) delete fDenominator_kT;
    fDenominator_kT = new TH2F(*aCorrFctn.fDenominator_kT);
  if(fNumerator_mT) delete fNumerator_mT;
    fNumerator_mT = new TH2F(*aCorrFctn.fNumerator_mT);
  if(fDenominator_mT) delete fDenominator_mT;
    fDenominator_mT = new TH2F(*aCorrFctn.fDenominator_mT);
  if(fNumeratorv2_mT) delete fNumeratorv2_mT;
    fNumeratorv2_mT = new TH2F(*aCorrFctn.fNumeratorv2_mT);
  if(fDenominatorv2_mT) delete fDenominatorv2_mT;
    fDenominatorv2_mT = new TH2F(*aCorrFctn.fDenominatorv2_mT);

  if(fNumerator3d) delete fNumerator3d;
    fNumerator3d = new TH3F(*aCorrFctn.fNumerator3d);
  if(fDenominator3d) delete fDenominator3d;
    fDenominator3d = new TH3F(*aCorrFctn.fDenominator3d);

  return *this;
}


//____________________________
AliFemtoString AliFemtoCorrFctnKStar::Report()
{
  // construct report
  TString report = "KStar Correlation Function Report:\n";
  report += TString::Format("Number of entries in numerator:\t%E\n", fNumerator->GetEntries());
  report += TString::Format("Number of entries in denominator:\t%E\n", fDenominator->GetEntries());
  report += TString::Format("Number of entries in ratio:\t%E\n", fRatio->GetEntries());
  return AliFemtoString(report);
}

//______________________________
TList* AliFemtoCorrFctnKStar::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumerator);
  tOutputList->Add(fDenominator);
  tOutputList->Add(fkTMonitor);
  tOutputList->Add(fRatio);
  if(fDetaDphiscal) {
    tOutputList->Add(fNumDEtaDPhiS);
    tOutputList->Add(fDenDEtaDPhiS);
  }
  if(fPairKinematics) {
    tOutputList->Add(fPairKStar);
  }
  if(fBuildkTBinned) {
    tOutputList->Add(fNumerator_kT);
    tOutputList->Add(fDenominator_kT);
  }
  if(fBuildmTBinned) {
    tOutputList->Add(fNumerator_mT);
    tOutputList->Add(fDenominator_mT);
    tOutputList->Add(fNumeratorv2_mT);
    tOutputList->Add(fDenominatorv2_mT);
  }
  if(fBuild3d) {
    tOutputList->Add(fNumerator3d);
    tOutputList->Add(fDenominator3d);
  }

  return tOutputList;
}

//_________________________
void AliFemtoCorrFctnKStar::Finish()
{
  fRatio->Divide(fNumerator, fDenominator, 1.0, 1.0);
}

//____________________________
void AliFemtoCorrFctnKStar::Write()
{
  // Write out neccessary objects
  fNumerator->Write();
  fDenominator->Write();
  fkTMonitor->Write();
  if(fDetaDphiscal) {
    fNumDEtaDPhiS->Write();
    fDenDEtaDPhiS->Write();
  }
  if(fPairKinematics) {
    fPairKStar->Write();
  }
  if(fBuildkTBinned) {
    fNumerator_kT->Write();
    fDenominator_kT->Write();
  }
  if(fBuildmTBinned) {
    fNumerator_mT->Write();
    fDenominator_mT->Write();
    fNumeratorv2_mT->Write();
    fDenominatorv2_mT->Write();
  }
  if(fBuild3d) {
    fNumerator3d->Write();
    fDenominator3d->Write();
  }
}


//____________________________
void AliFemtoCorrFctnKStar::AddRealPair(AliFemtoPair* aPair)
{
  // add true pair
  if (fPairCut && !fPairCut->Pass(aPair)) {
    return;
  }

  double tKStar = fabs(aPair->KStar());

  fNumerator->Fill(tKStar);
  fkTMonitor->Fill(aPair->KT());

  if(fDetaDphiscal) FillDEtaDPhiS(fNumDEtaDPhiS,aPair);
  if(fBuildkTBinned) fNumerator_kT->Fill(aPair->KStar(), aPair->KT());
  if(fBuildmTBinned) fNumerator_mT->Fill(aPair->KStar(), CalcMt(aPair));
  if(fBuildmTBinned) fNumeratorv2_mT->Fill(aPair->KStar(), CalcMtv2(aPair));
  if(fBuild3d) fNumerator3d->Fill(aPair->KStarOut(),aPair->KStarSide(),aPair->KStarLong());
}

//____________________________
void AliFemtoCorrFctnKStar::AddMixedPair(AliFemtoPair* aPair)
{
  // add mixed (background) pair
  if (fPairCut && !fPairCut->Pass(aPair)) {
    return;
  }

  double weight = 1.0;
  double tKStar = fabs(aPair->KStar());
  fDenominator->Fill(tKStar,weight);

  if(fPairKinematics) fPairKStar->Fill(aPair->KStar(),aPair->KOut(),aPair->KSide(),aPair->KLong());
  if(fDetaDphiscal) FillDEtaDPhiS(fDenDEtaDPhiS,aPair);
  if(fBuildkTBinned) fDenominator_kT->Fill(aPair->KStar(), aPair->KT());
  if(fBuildmTBinned) fDenominator_mT->Fill(aPair->KStar(), CalcMt(aPair));
  if(fBuildmTBinned) fDenominatorv2_mT->Fill(aPair->KStar(), CalcMtv2(aPair));
  if(fBuild3d) fDenominator3d->Fill(aPair->KStarOut(),aPair->KStarSide(),aPair->KStarLong());
}


//____________________________
void AliFemtoCorrFctnKStar::FillDEtaDPhiS(TH2D* aHist, AliFemtoPair* aPair)
{
  double phi1 = aPair->Track1()->Track()->P().Phi();
  double phi2 = aPair->Track2()->Track()->P().Phi();
  double chg1 = aPair->Track1()->Track()->Charge();
  double chg2 = aPair->Track2()->Track()->Charge();
  double ptv1 = aPair->Track1()->Track()->Pt();
  double ptv2 = aPair->Track2()->Track()->Pt();
  double eta1 = aPair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = aPair->Track2()->Track()->P().PseudoRapidity();

  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  Int_t magsign = 0;

  if (!aodH) {
    //AliWarning("Could not get AODInputHandler");
  }
  else {
    AliAODEvent *fAOD;
    fAOD = aodH->GetEvent();
    magsign = fAOD->GetMagneticField();
  }

  Int_t fMagSign;
  if (magsign > 1) fMagSign = 1;
  else if ( magsign < 1) fMagSign = -1;
  else fMagSign = magsign;

  Double_t rad = fRaddedps;

  double afsi0b = 0.07510020733*chg1*fMagSign*rad/ptv1;
  double afsi1b = 0.07510020733*chg2*fMagSign*rad/ptv2;
  Double_t dps6 =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
  dps6 = TVector2::Phi_mpi_pi(dps6);
  double etad = eta2 - eta1;

  aHist->Fill(dps6,etad);
}

//____________________________
void AliFemtoCorrFctnKStar::SetkTMonitorBins(int aNbinskT, double akTMin, double akTMax)
{
  fkTMonitor->SetBins(aNbinskT,akTMin,akTMax);
}

//____________________________
void AliFemtoCorrFctnKStar::SetKStarVskTBins(int aNbinsKStar, double aKStarMin, double aKStarMax,
                                             int aNbinskT,    double akTMin,    double akTMax)
{
  if(!fNumerator_kT || !fDenominator_kT) SetBuildkTBinned(true);
  fNumerator_kT->SetBins(aNbinsKStar,aKStarMin,aKStarMax, aNbinskT,akTMin,akTMax);
  fDenominator_kT->SetBins(aNbinsKStar,aKStarMin,aKStarMax, aNbinskT,akTMin,akTMax);
}

//____________________________
void AliFemtoCorrFctnKStar::SetKStarVsmTBins(int aNbinsKStar, double aKStarMin, double aKStarMax,
                                             int aNbinsmT,    double amTMin,    double amTMax)
{
  if(!fNumerator_mT || !fDenominator_mT) SetBuildmTBinned(true);
  fNumerator_mT->SetBins(aNbinsKStar,aKStarMin,aKStarMax, aNbinsmT,amTMin,amTMax);
  fDenominator_mT->SetBins(aNbinsKStar,aKStarMin,aKStarMax, aNbinsmT,amTMin,amTMax);
  fNumeratorv2_mT->SetBins(aNbinsKStar,aKStarMin,aKStarMax, aNbinsmT,amTMin,amTMax);
  fDenominatorv2_mT->SetBins(aNbinsKStar,aKStarMin,aKStarMax, aNbinsmT,amTMin,amTMax);
}

//____________________________
void AliFemtoCorrFctnKStar::SetDEtaDPhiSBins(int aNbinsDEta,  double aDEtaMin,  double aDEtaMax,
                                             int aNbinsDPhis, double aDPhiSMin, double aDPhiSMax)
{
  if(!fNumDEtaDPhiS || !fDenDEtaDPhiS) SetCalculateDetaDphis(true,fRaddedps);
  fNumDEtaDPhiS->SetBins(aNbinsDPhis,aDPhiSMin,aDPhiSMax, aNbinsDEta,aDEtaMin,aDEtaMax);
  fDenDEtaDPhiS->SetBins(aNbinsDPhis,aDPhiSMin,aDPhiSMax, aNbinsDEta,aDEtaMin,aDEtaMax);
}

//____________________________
void AliFemtoCorrFctnKStar::Set3dBins(int aNbinsKStarOut,  double aKStarOutMin,  double aKStarOutMax,
                                      int aNbinsKStarSide, double aKStarSideMin, double aKStarSideMax,
                                      int aNbinsKStarLong, double aKStarLongMin, double aKStarLongMax)
{
  if(!fNumerator3d || !fDenominator3d) SetBuild3d(true);
  fNumerator3d->SetBins(aNbinsKStarOut, aKStarOutMin, aKStarOutMax,
                        aNbinsKStarSide,aKStarSideMin,aKStarSideMax,
                        aNbinsKStarLong,aKStarLongMin,aKStarLongMax);

  fDenominator3d->SetBins(aNbinsKStarOut, aKStarOutMin, aKStarOutMax,
                        aNbinsKStarSide,aKStarSideMin,aKStarSideMax,
                        aNbinsKStarLong,aKStarLongMin,aKStarLongMax);
}

//____________________________
float AliFemtoCorrFctnKStar::CalcMt(const AliFemtoPair* aPair)
{
  const double mass_1 = aPair->Track1()->FourMomentum().m(),
               mass_2 = aPair->Track2()->FourMomentum().m();
  const double avg_mass = (mass_1 + mass_2) / 2.0;
  return TMath::Sqrt(avg_mass * avg_mass + ::pow(aPair->KT(), 2));
}

//____________________________
float AliFemtoCorrFctnKStar::CalcMtv2(const AliFemtoPair* aPair)
{
//  return aPair->FourMomentumSum().mt();  //This is not correct for this application
  const double mass_1 = aPair->Track1()->FourMomentum().m(),
               mass_2 = aPair->Track2()->FourMomentum().m();
  const double avg_mass = (mass_1 * mass_2)/(mass_1 + mass_2);
  return TMath::Sqrt(avg_mass * avg_mass + ::pow(aPair->KT(), 2));
}

//____________________________
void AliFemtoCorrFctnKStar::SetCalculateDetaDphis(Bool_t dedpsc, Double_t rad)
{
  fDetaDphiscal = dedpsc;
  fRaddedps = rad;

  if(fDetaDphiscal && (!fNumDEtaDPhiS || !fDenDEtaDPhiS))
  {
    fNumDEtaDPhiS = new TH2D(TString::Format("NumDEtaDPhiS%s", fTitle.Data()),
                             "dPhiS vs dEta - Numerator",
                             DFLT_NbinsDPhiS,DFLT_DPhiSMin,DFLT_DPhiSMax,
                             DFLT_NbinsDEta,DFLT_DEtaMin,DFLT_DEtaMax);
    fDenDEtaDPhiS = new TH2D(TString::Format("DenDEtaDPhiS%s", fTitle.Data()),
                             "dPhiS vs dEta - Denominator",
                             DFLT_NbinsDPhiS,DFLT_DPhiSMin,DFLT_DPhiSMax,
                             DFLT_NbinsDEta,DFLT_DEtaMin,DFLT_DEtaMax);
    fNumDEtaDPhiS->Sumw2();
    fDenDEtaDPhiS->Sumw2();
  }
}

//____________________________
void AliFemtoCorrFctnKStar::SetCalculatePairKinematics(Bool_t pk)
{
  fPairKinematics = pk;

  if (fPairKinematics && !fPairKStar) {
    fPairKStar = new TNtuple(TString::Format("PairKStar%s", fTitle.Data()),
                             "Pair Kinematics",
                             "KStarMag:KStarOut:KStarSide:KStarLong");
  }
}

//____________________________
void AliFemtoCorrFctnKStar::SetBuildkTBinned(Bool_t aBuild)
{
  fBuildkTBinned = aBuild;

  if(fBuildkTBinned && (!fNumerator_kT || !fDenominator_kT))
  {
    fNumerator_kT = new TH2F(TString::Format("KStarVskTNum%s", fTitle.Data()),
                             "KStar vs kT Numerator; k*(GeV); k_T (GeV);",
                             fNbinsKStar, fKStarLow, fKStarHigh,
                             DFLT_NbinskT, DFLT_kTMin, DFLT_kTMax);

    fDenominator_kT = new TH2F(TString::Format("KStarVskTDen%s", fTitle.Data()),
                           "KStar vs kT Denominator; k* (GeV); k_T (GeV)",
                           fNbinsKStar, fKStarLow, fKStarHigh,
                           DFLT_NbinskT, DFLT_kTMin, DFLT_kTMax);
    fNumerator_kT->Sumw2();
    fDenominator_kT->Sumw2();
  }
}

//____________________________
void AliFemtoCorrFctnKStar::SetBuildmTBinned(Bool_t aBuild)
{
  fBuildmTBinned = aBuild;

  if(fBuildmTBinned && (!fNumerator_mT || !fDenominator_mT))
  {
    fNumerator_mT = new TH2F(TString::Format("KStarVsmTNum%s", fTitle.Data()),
                             "KStar vs m_{T} Numerator; k*(GeV); m_{T} (GeV);",
                             fNbinsKStar, fKStarLow, fKStarHigh,
                             DFLT_NbinsmT, DFLT_mTMin, DFLT_mTMax);
    fDenominator_mT = new TH2F(TString::Format("KStarVsmTDen%s", fTitle.Data()),
                               "KStar vs m_{T} Denominator; k* (GeV); m_{T} (GeV)",
                                fNbinsKStar, fKStarLow, fKStarHigh,
                                DFLT_NbinsmT, DFLT_mTMin, DFLT_mTMax);
    fNumerator_mT->Sumw2();
    fDenominator_mT->Sumw2();

    fNumeratorv2_mT = new TH2F(TString::Format("KStarVsmTNumv2%s", fTitle.Data()),
                             "KStar vs m_{T} Numerator; k*(GeV); m_{T} (GeV);",
                             fNbinsKStar, fKStarLow, fKStarHigh,
                             DFLT_NbinsmT, DFLT_mTMin, DFLT_mTMax);
    fDenominatorv2_mT = new TH2F(TString::Format("KStarVsmTDenv2%s", fTitle.Data()),
                               "KStar vs m_{T} Denominator; k* (GeV); m_{T} (GeV)",
                                fNbinsKStar, fKStarLow, fKStarHigh,
                                DFLT_NbinsmT, DFLT_mTMin, DFLT_mTMax);
    fNumeratorv2_mT->Sumw2();
    fDenominatorv2_mT->Sumw2();
  }
}

//____________________________
void AliFemtoCorrFctnKStar::SetBuild3d(Bool_t aBuild)
{
  fBuild3d = aBuild;

  if(fBuild3d && (!fNumerator3d || !fDenominator3d))
  {
    fNumerator3d = new TH3F(TString::Format("Num3d%s", fTitle.Data()),
                            "KStar3d - Numerator; k*(GeV/c);",
                            DFLT_NbinsKStarOut, DFLT_KStarOutMin, DFLT_KStarOutMax,
                            DFLT_NbinsKStarSide, DFLT_KStarSideMin, DFLT_KStarSideMax,
                            DFLT_NbinsKStarLong, DFLT_KStarLongMin, DFLT_KStarLongMax);
    fDenominator3d = new TH3F(TString::Format("Den3d%s", fTitle.Data()),
                              "KStar3d - Denominator; k*(GeV/c);",
                              DFLT_NbinsKStarOut, DFLT_KStarOutMin, DFLT_KStarOutMax,
                              DFLT_NbinsKStarSide, DFLT_KStarSideMin, DFLT_KStarSideMax,
                              DFLT_NbinsKStarLong, DFLT_KStarLongMin, DFLT_KStarLongMax);
    fNumerator3d->Sumw2();
    fDenominator3d->Sumw2();
  }
}
