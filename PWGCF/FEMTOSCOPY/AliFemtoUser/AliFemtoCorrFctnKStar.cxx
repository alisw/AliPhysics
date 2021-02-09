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

  fNumerator_RotatePar2(nullptr),

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

  fBuildIndmTBinned(kFALSE),
  fNumerator_IndmT1(nullptr),
  fDenominator_IndmT1(nullptr),
  fNumerator_IndmT2(nullptr),
  fDenominator_IndmT2(nullptr),

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

  fNumerator_RotatePar2 = new TH1D(TString::Format("Num_RotatePar2%s", fTitle.Data()),
                                   "KStar - Numerator - Rotate Particle 2; k*(GeV/c);",
                                    nbins, KStarLo, KStarHi);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  fRatio->Sumw2();
  fkTMonitor->Sumw2();

  fNumerator_RotatePar2->Sumw2();
}

//____________________________
AliFemtoCorrFctnKStar::AliFemtoCorrFctnKStar(const AliFemtoCorrFctnKStar& aCorrFctn):
  AliFemtoCorrFctn(aCorrFctn),
  fTitle(aCorrFctn.fTitle),
  fNbinsKStar(aCorrFctn.fNbinsKStar),
  fKStarLow(aCorrFctn.fKStarLow),
  fKStarHigh(aCorrFctn.fKStarHigh),

  fNumerator(aCorrFctn.fNumerator ? new TH1D(*aCorrFctn.fNumerator) : nullptr),
  fDenominator(aCorrFctn.fDenominator ? new TH1D(*aCorrFctn.fDenominator) : nullptr),

  fDetaDphiscal(aCorrFctn.fDetaDphiscal),
  fPairKinematics(aCorrFctn.fPairKinematics),
  fRaddedps(aCorrFctn.fRaddedps),
  fBuildkTBinned(aCorrFctn.fBuildkTBinned),
  fBuildmTBinned(aCorrFctn.fBuildmTBinned),
  fBuildIndmTBinned(aCorrFctn.fBuildIndmTBinned),
  fBuild3d(aCorrFctn.fBuild3d)
{
  // copy constructor
  fRatio = (aCorrFctn.fRatio) ? new TH1D(*aCorrFctn.fRatio) : nullptr;
  fkTMonitor =(aCorrFctn.fkTMonitor) ? static_cast<TH1D*>(aCorrFctn.fkTMonitor->Clone()) : nullptr;

  fNumerator_RotatePar2 =(aCorrFctn.fNumerator_RotatePar2) ? static_cast<TH1D*>(aCorrFctn.fNumerator_RotatePar2->Clone()) : nullptr;

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

  if(aCorrFctn.fNumerator_IndmT1) fNumerator_IndmT1 = new TH2F(*aCorrFctn.fNumerator_IndmT1);
    else fNumerator_IndmT1 = 0;
  if(aCorrFctn.fDenominator_IndmT1) fDenominator_IndmT1 = new TH2F(*aCorrFctn.fDenominator_IndmT1);
    else fDenominator_IndmT1 = 0;

  if(aCorrFctn.fNumerator_IndmT2) fNumerator_IndmT2 = new TH2F(*aCorrFctn.fNumerator_IndmT2);
    else fNumerator_IndmT2 = 0;
  if(aCorrFctn.fDenominator_IndmT2) fDenominator_IndmT2 = new TH2F(*aCorrFctn.fDenominator_IndmT2);
    else fDenominator_IndmT2 = 0;

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
  delete fNumerator_RotatePar2;
  delete fNumDEtaDPhiS;
  delete fDenDEtaDPhiS;
  delete fPairKStar;
  delete fNumerator_kT;
  delete fDenominator_kT;
  delete fNumerator_mT;
  delete fDenominator_mT;
  delete fNumerator_IndmT1;
  delete fDenominator_IndmT1;
  delete fNumerator_IndmT2;
  delete fDenominator_IndmT2;
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
  fBuildIndmTBinned = aCorrFctn.fBuildIndmTBinned;
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

  if(fNumerator_RotatePar2) delete fNumerator_RotatePar2;
    fNumerator_RotatePar2 = new TH1D(*aCorrFctn.fNumerator_RotatePar2);

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

  if(fNumerator_IndmT1) delete fNumerator_IndmT1;
    fNumerator_IndmT1 = new TH2F(*aCorrFctn.fNumerator_IndmT1);
  if(fDenominator_IndmT1) delete fDenominator_IndmT1;
    fDenominator_IndmT1 = new TH2F(*aCorrFctn.fDenominator_IndmT1);
  if(fNumerator_IndmT2) delete fNumerator_IndmT2;
    fNumerator_IndmT2 = new TH2F(*aCorrFctn.fNumerator_IndmT2);
  if(fDenominator_IndmT2) delete fDenominator_IndmT2;
    fDenominator_IndmT2 = new TH2F(*aCorrFctn.fDenominator_IndmT2);

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
  return AliFemtoString((const char *)report);
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

  tOutputList->Add(fNumerator_RotatePar2);

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
  }
  if(fBuildIndmTBinned) {
    tOutputList->Add(fNumerator_IndmT1);
    tOutputList->Add(fDenominator_IndmT1);
    tOutputList->Add(fNumerator_IndmT2);
    tOutputList->Add(fDenominator_IndmT2);
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

  fNumerator_RotatePar2->Write();

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
  }
  if(fBuildIndmTBinned) {
    fNumerator_IndmT1->Write();
    fDenominator_IndmT1->Write();
    fNumerator_IndmT2->Write();
    fDenominator_IndmT2->Write();
  }
  if(fBuild3d) {
    fNumerator3d->Write();
    fDenominator3d->Write();
  }
}


//____________________________
void AliFemtoCorrFctnKStar::AddRealPair(AliFemtoPair* aPair)
{
  if(PassPairCut_RotatePar2(aPair)) fNumerator_RotatePar2->Fill(fabs(CalcKStar_RotatePar2(aPair)));

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
  if(fBuildIndmTBinned)
  {
    fNumerator_IndmT1->Fill(aPair->KStar(), aPair->Track1()->FourMomentum().mt());
    fNumerator_IndmT2->Fill(aPair->KStar(), aPair->Track2()->FourMomentum().mt());
  }
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
  if(fBuildIndmTBinned)
  {
    fDenominator_IndmT1->Fill(aPair->KStar(), aPair->Track1()->FourMomentum().mt());
    fDenominator_IndmT2->Fill(aPair->KStar(), aPair->Track2()->FourMomentum().mt());
  }
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
}

//____________________________
void AliFemtoCorrFctnKStar::SetKStarVsIndmTBins(int aNbinsKStar, double aKStarMin, double aKStarMax,
                                                int aNbinsmT1,    double amTMin1,    double amTMax1,
                                                int aNbinsmT2,    double amTMin2,    double amTMax2)
{
  if(!fNumerator_IndmT1 || !fDenominator_IndmT1 || !fNumerator_IndmT2 || !fDenominator_IndmT2) SetBuildIndmTBinned(true);
  fNumerator_IndmT1->SetBins(aNbinsKStar,aKStarMin,aKStarMax, aNbinsmT1,amTMin1,amTMax1);
  fDenominator_IndmT1->SetBins(aNbinsKStar,aKStarMin,aKStarMax, aNbinsmT1,amTMin1,amTMax1);
  fNumerator_IndmT2->SetBins(aNbinsKStar,aKStarMin,aKStarMax, aNbinsmT2,amTMin2,amTMax2);
  fDenominator_IndmT2->SetBins(aNbinsKStar,aKStarMin,aKStarMax, aNbinsmT2,amTMin2,amTMax2);
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
void AliFemtoCorrFctnKStar::RotateThreeVecBy180InTransversePlane(AliFemtoThreeVector &a3Vec)
{
  a3Vec.SetX(-1.*a3Vec.x());
  a3Vec.SetY(-1.*a3Vec.y());
}

//____________________________
bool AliFemtoCorrFctnKStar::PassPairCut_RotatePar2(const AliFemtoPair* aPair)
{
  if(!fPairCut) return true;

  AliFemtoParticle *tPart1 = new AliFemtoParticle(*aPair->Track1());
  AliFemtoParticle *tPart2 = new AliFemtoParticle(*aPair->Track2());

  AliFemtoPair* tPair = new AliFemtoPair;
  tPair->SetTrack1(tPart1);

  //------------------------------
  AliFemtoLorentzVector tFourMom2 = AliFemtoLorentzVector(tPart2->FourMomentum());
  tFourMom2.SetPx(-1.*tFourMom2.px());
  tFourMom2.SetPy(-1.*tFourMom2.py());
  tPart2->ResetFourMomentum(tFourMom2);
  //------------------------------

  AliFemtoThreeVector tTpcThreeVec;

  if(tPart2->Track()) //tPart2 is a track
  {
    double **tTpcPositions;
    tTpcPositions = new double*[9];
    for (int i = 0; i < 9; i++) tTpcPositions[i] = new double[3];

    tTpcThreeVec = tPart2->Track()->NominalTpcEntrancePoint();
    RotateThreeVecBy180InTransversePlane(tTpcThreeVec);
    tPart2->Track()->SetNominalTPCEntrancePoint(tTpcThreeVec);

    tTpcThreeVec = tPart2->Track()->NominalTpcExitPoint();
    RotateThreeVecBy180InTransversePlane(tTpcThreeVec);
    tPart2->Track()->SetNominalTPCExitPoint(tTpcThreeVec);

    for (int i = 0; i < 9; i++)
    {
      tTpcThreeVec = tPart2->Track()->NominalTpcPoint(i);
      tTpcPositions[i][0] = -1.*tTpcThreeVec.x();
      tTpcPositions[i][1] = -1.*tTpcThreeVec.y();
      tTpcPositions[i][2] = tTpcThreeVec.z();
    }
    tPart2->Track()->SetNominalTPCPoints(tTpcPositions);

    for (int i = 0; i < 9; i++) delete [] tTpcPositions[i];
    delete [] tTpcPositions;
  }
  else if(tPart2->V0()) //tPart2 is a V0
  {
    //-----Positive daughter
    tTpcThreeVec = tPart2->V0()->NominalTpcEntrancePointPos();
    RotateThreeVecBy180InTransversePlane(tTpcThreeVec);
    tPart2->V0()->SetNominalTpcEntrancePointPos(tTpcThreeVec);

    tTpcThreeVec = tPart2->V0()->NominalTpcExitPointPos();
    RotateThreeVecBy180InTransversePlane(tTpcThreeVec);
    tPart2->V0()->SetNominalTpcExitPointPos(tTpcThreeVec);

    AliFemtoThreeVector tTpcThreeVecArr[9];
    for (int i = 0; i < 9; i++)
    {
      tTpcThreeVecArr[i].SetX(-1.*tPart2->V0()->NominalTpcPointPos(i).x());
      tTpcThreeVecArr[i].SetY(-1.*tPart2->V0()->NominalTpcPointPos(i).y());
      tTpcThreeVecArr[i].SetZ(tPart2->V0()->NominalTpcPointPos(i).z());
    }
    tPart2->V0()->SetNominalTpcPointPos(tTpcThreeVecArr);

    //-----Negative daughter
    tTpcThreeVec = tPart2->V0()->NominalTpcEntrancePointNeg();
    RotateThreeVecBy180InTransversePlane(tTpcThreeVec);
    tPart2->V0()->SetNominalTpcEntrancePointNeg(tTpcThreeVec);

    tTpcThreeVec = tPart2->V0()->NominalTpcExitPointNeg();
    RotateThreeVecBy180InTransversePlane(tTpcThreeVec);
    tPart2->V0()->SetNominalTpcExitPointNeg(tTpcThreeVec);

    for (int i = 0; i < 9; i++)
    {
      tTpcThreeVecArr[i].SetX(-1.*tPart2->V0()->NominalTpcPointNeg(i).x());
      tTpcThreeVecArr[i].SetY(-1.*tPart2->V0()->NominalTpcPointNeg(i).y());
      tTpcThreeVecArr[i].SetZ(tPart2->V0()->NominalTpcPointNeg(i).z());
    }
    tPart2->V0()->SetNominalTpcPointNeg(tTpcThreeVecArr);
  }
  else if(tPart2->Xi()) //tPart2 is a Xi
  {
    //-----Positive V0-daughter
    tTpcThreeVec = tPart2->Xi()->NominalTpcEntrancePointPos();
    RotateThreeVecBy180InTransversePlane(tTpcThreeVec);
    tPart2->Xi()->SetNominalTpcEntrancePointPos(tTpcThreeVec);

    tTpcThreeVec = tPart2->Xi()->NominalTpcExitPointPos();
    RotateThreeVecBy180InTransversePlane(tTpcThreeVec);
    tPart2->Xi()->SetNominalTpcExitPointPos(tTpcThreeVec);

    AliFemtoThreeVector tTpcThreeVecArr[9];
    for (int i = 0; i < 9; i++)
    {
      tTpcThreeVecArr[i].SetX(-1.*tPart2->Xi()->NominalTpcPointPos(i).x());
      tTpcThreeVecArr[i].SetY(-1.*tPart2->Xi()->NominalTpcPointPos(i).y());
      tTpcThreeVecArr[i].SetZ(tPart2->Xi()->NominalTpcPointPos(i).z());
    }
    tPart2->Xi()->SetNominalTpcPointPos(tTpcThreeVecArr);

    //-----Negative V0-daughter
    tTpcThreeVec = tPart2->Xi()->NominalTpcEntrancePointNeg();
    RotateThreeVecBy180InTransversePlane(tTpcThreeVec);
    tPart2->Xi()->SetNominalTpcEntrancePointNeg(tTpcThreeVec);

    tTpcThreeVec = tPart2->Xi()->NominalTpcExitPointNeg();
    RotateThreeVecBy180InTransversePlane(tTpcThreeVec);
    tPart2->Xi()->SetNominalTpcExitPointNeg(tTpcThreeVec);

    for (int i = 0; i < 9; i++)
    {
      tTpcThreeVecArr[i].SetX(-1.*tPart2->Xi()->NominalTpcPointNeg(i).x());
      tTpcThreeVecArr[i].SetY(-1.*tPart2->Xi()->NominalTpcPointNeg(i).y());
      tTpcThreeVecArr[i].SetZ(tPart2->Xi()->NominalTpcPointNeg(i).z());
    }
    tPart2->Xi()->SetNominalTpcPointNeg(tTpcThreeVecArr);

    //-----Bachelor daughter
    tTpcThreeVec = tPart2->Xi()->NominalTpcEntrancePointBac();
    RotateThreeVecBy180InTransversePlane(tTpcThreeVec);
    tPart2->Xi()->SetNominalTpcEntrancePointBac(tTpcThreeVec);

    tTpcThreeVec = tPart2->Xi()->NominalTpcExitPointBac();
    RotateThreeVecBy180InTransversePlane(tTpcThreeVec);
    tPart2->Xi()->SetNominalTpcExitPointBac(tTpcThreeVec);

    for (int i = 0; i < 9; i++)
    {
      tTpcThreeVecArr[i].SetX(-1.*tPart2->Xi()->NominalTpcPointBac(i).x());
      tTpcThreeVecArr[i].SetY(-1.*tPart2->Xi()->NominalTpcPointBac(i).y());
      tTpcThreeVecArr[i].SetZ(tPart2->Xi()->NominalTpcPointBac(i).z());
    }
    tPart2->Xi()->SetNominalTpcPointBac(tTpcThreeVecArr);
  }
  else return false;

  //------------------------------
  tPair->SetTrack2(tPart2);
  bool tPass = fPairCut->Pass(tPair);

  delete tPart1;
  delete tPart2;
  delete tPair;

  return tPass;
}

//____________________________
double AliFemtoCorrFctnKStar::CalcKStar_RotatePar2(const AliFemtoPair* aPair)
{
  double tKStarCalc = 0.;

  const AliFemtoLorentzVector p1 = aPair->Track1()->FourMomentum();
  const AliFemtoLorentzVector p2 = aPair->Track2()->FourMomentum();

  AliFemtoLorentzVector p2_Rot = AliFemtoLorentzVector(p2);
  p2_Rot.SetPx(-1.*p2_Rot.px());
  p2_Rot.SetPy(-1.*p2_Rot.py());

  const double p_inv = (p1 + p2_Rot).m2(),
               q_inv = (p1 - p2_Rot).m2(),
           mass_diff = p1.m2() - p2_Rot.m2();

  const double tQ = ::pow(mass_diff, 2) / p_inv - q_inv;
  tKStarCalc = ::sqrt(tQ) / 2.0;

  return tKStarCalc;
}


//____________________________
float AliFemtoCorrFctnKStar::CalcMt(const AliFemtoPair* aPair)
{
  return 0.5*aPair->FourMomentumSum().mt();
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
  }
}

//____________________________
void AliFemtoCorrFctnKStar::SetBuildIndmTBinned(Bool_t aBuild)
{
  fBuildIndmTBinned = aBuild;

  if(fBuildIndmTBinned && (!fNumerator_IndmT1 || !fDenominator_IndmT1 || !fNumerator_IndmT2 || !fDenominator_IndmT2))
  {
    fNumerator_IndmT1 = new TH2F(TString::Format("KStarVsIndmT1Num%s", fTitle.Data()),
                                 "KStar vs m_{T} of particle1 Numerator; k*(GeV); m_{T} (GeV);",
                                 fNbinsKStar, fKStarLow, fKStarHigh,
                                 DFLT_NbinsmT, DFLT_mTMin, DFLT_mTMax);
    fDenominator_IndmT1 = new TH2F(TString::Format("KStarVsIndmT1Den%s", fTitle.Data()),
                                   "KStar vs m_{T} of particle1 Denominator; k* (GeV); m_{T} (GeV)",
                                   fNbinsKStar, fKStarLow, fKStarHigh,
                                   DFLT_NbinsmT, DFLT_mTMin, DFLT_mTMax);
    fNumerator_IndmT1->Sumw2();
    fDenominator_IndmT1->Sumw2();

    fNumerator_IndmT2 = new TH2F(TString::Format("KStarVsIndmT2Num%s", fTitle.Data()),
                                 "KStar vs m_{T} of particle2 Numerator; k*(GeV); m_{T} (GeV);",
                                 fNbinsKStar, fKStarLow, fKStarHigh,
                                 DFLT_NbinsmT, DFLT_mTMin, DFLT_mTMax);
    fDenominator_IndmT2 = new TH2F(TString::Format("KStarVsIndmT2Den%s", fTitle.Data()),
                                   "KStar vs m_{T} of particle2 Denominator; k* (GeV); m_{T} (GeV)",
                                   fNbinsKStar, fKStarLow, fKStarHigh,
                                   DFLT_NbinsmT, DFLT_mTMin, DFLT_mTMax);
    fNumerator_IndmT2->Sumw2();
    fDenominator_IndmT2->Sumw2();
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
