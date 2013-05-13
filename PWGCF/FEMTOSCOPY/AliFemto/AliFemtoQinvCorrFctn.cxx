///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoQinvCorrFctn:                                                 //
// a simple Q-invariant correlation function                             //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliFemtoQinvCorrFctn.h"
//#include "AliFemtoHisto.h"
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoQinvCorrFctn)
#endif

//____________________________
AliFemtoQinvCorrFctn::AliFemtoQinvCorrFctn(char* title, const int& nbins, const float& QinvLo, const float& QinvHi):
  fNumerator(0),
  fDenominator(0),
  fRatio(0),
  fkTMonitor(0),
  fDetaDphiscal(kFALSE),
  fRaddedps(1.2),
  fNumDEtaDPhiS(0),
  fDenDEtaDPhiS(0)
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
AliFemtoQinvCorrFctn::AliFemtoQinvCorrFctn(const AliFemtoQinvCorrFctn& aCorrFctn) :
  AliFemtoCorrFctn(),
  fNumerator(0),
  fDenominator(0),
  fRatio(0),
  fkTMonitor(0),
  fDetaDphiscal(kFALSE),
  fRaddedps(1.2),
  fNumDEtaDPhiS(0),
  fDenDEtaDPhiS(0)
{
  // copy constructor
  fNumerator = new TH1D(*aCorrFctn.fNumerator);
  fDenominator = new TH1D(*aCorrFctn.fDenominator);
  fRatio = new TH1D(*aCorrFctn.fRatio);
  fkTMonitor = new TH1D(*aCorrFctn.fkTMonitor);

  fNumDEtaDPhiS = new TH2D(*aCorrFctn.fNumDEtaDPhiS);
  fDenDEtaDPhiS = new TH2D(*aCorrFctn.fDenDEtaDPhiS);

  fDetaDphiscal = aCorrFctn.fDetaDphiscal;
  fRaddedps = aCorrFctn.fRaddedps;

}
//____________________________
AliFemtoQinvCorrFctn::~AliFemtoQinvCorrFctn(){
  // destructor
  delete fNumerator;
  delete fDenominator;
  delete fRatio;
  delete fkTMonitor;
  delete fNumDEtaDPhiS;
  delete fDenDEtaDPhiS;

}
//_________________________
AliFemtoQinvCorrFctn& AliFemtoQinvCorrFctn::operator=(const AliFemtoQinvCorrFctn& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (fNumerator) delete fNumerator;
  fNumerator = new TH1D(*aCorrFctn.fNumerator);
  if (fDenominator) delete fDenominator;
  fDenominator = new TH1D(*aCorrFctn.fDenominator);
  if (fRatio) delete fRatio;
  fRatio = new TH1D(*aCorrFctn.fRatio);
  if (fkTMonitor) delete fkTMonitor;
  fkTMonitor = new TH1D(*aCorrFctn.fkTMonitor);

  if (fNumDEtaDPhiS) delete fNumDEtaDPhiS;
  fNumDEtaDPhiS = new TH2D(*aCorrFctn.fNumDEtaDPhiS);
  if (fDenDEtaDPhiS) delete fDenDEtaDPhiS;
  fDenDEtaDPhiS = new TH2D(*aCorrFctn.fDenDEtaDPhiS);

  fDetaDphiscal = aCorrFctn.fDetaDphiscal;
  fRaddedps = aCorrFctn.fRaddedps;

  return *this;
}

//_________________________
void AliFemtoQinvCorrFctn::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  fNumerator->Draw();
  //fDenominator->Draw();
  //fRatio->Draw();
  fRatio->Divide(fNumerator,fDenominator,1.0,1.0);

}

//____________________________
AliFemtoString AliFemtoQinvCorrFctn::Report(){
  // construct report
  string stemp = "Qinv Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fNumerator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fDenominator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in ratio:\t%E\n",fRatio->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoQinvCorrFctn::AddRealPair(AliFemtoPair* pair){
  // add true pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...

  fNumerator->Fill(tQinv);
  fkTMonitor->Fill(pair->KT());


//_______________________________________
  if (fDetaDphiscal) {

    double phi1 = pair->Track1()->Track()->P().Phi();
    double phi2 = pair->Track2()->Track()->P().Phi();
    double chg1 = pair->Track1()->Track()->Charge();
    double chg2 = pair->Track2()->Track()->Charge();
    double ptv1 = pair->Track1()->Track()->Pt();
    double ptv2 = pair->Track2()->Track()->Pt();
    double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
    double eta2 = pair->Track2()->Track()->P().PseudoRapidity();

    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    AliAODEvent *fAOD;

    if (!aodH) {
      //AliWarning("Could not get AODInputHandler");
    }
    else {
      fAOD = aodH->GetEvent();
    }

    Int_t magsign = fAOD->GetMagneticField();
    Int_t fMagSign;

    if (magsign > 1)
      fMagSign = 1;
    else if ( magsign < 1)
      fMagSign = -1;
    else
      fMagSign = magsign;

    Double_t rad = fRaddedps;

    double afsi0b = 0.07510020733*chg1*fMagSign*rad/ptv1;
    double afsi1b = 0.07510020733*chg2*fMagSign*rad/ptv2;
    Double_t dps6 =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
    dps6 = TVector2::Phi_mpi_pi(dps6);

    // Double_t dps = (phi1-phi2+(TMath::ASin(-0.075*chg1*fMagSign*rad/ptv1))-(TMath::ASin(-0.075*chg2*fMagSign*rad/ptv2)));
    // dps = TVector2::Phi_mpi_pi(dps);
    double etad = eta2 - eta1;

    fNumDEtaDPhiS->Fill(dps6,etad);
  }
//_______________________________________________________________

  //  cout << "AliFemtoQinvCorrFctn::AddRealPair : " << pair->qInv() << " " << tQinv <<
  //" " << pair->track1().FourMomentum() << " " << pair->track2().FourMomentum() << endl;
}

//____________________________
void AliFemtoQinvCorrFctn::AddMixedPair(AliFemtoPair* pair){
  // add mixed (background) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  double weight = 1.0;
  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  fDenominator->Fill(tQinv,weight);

//_______________________________________
  if (fDetaDphiscal) {

    double phi1 = pair->Track1()->Track()->P().Phi();
    double phi2 = pair->Track2()->Track()->P().Phi();
    double chg1 = pair->Track1()->Track()->Charge();
    double chg2 = pair->Track2()->Track()->Charge();
    double ptv1 = pair->Track1()->Track()->Pt();
    double ptv2 = pair->Track2()->Track()->Pt();
    double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
    double eta2 = pair->Track2()->Track()->P().PseudoRapidity();

    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    AliAODEvent *fAOD;

    if (!aodH) {
      //AliWarning("Could not get AODInputHandler");
      cout << "Could not get AODInputHandler" << endl;
    }
    else {
      fAOD = aodH->GetEvent();
    }

    Int_t magsign = fAOD->GetMagneticField();
    Int_t fMagSign;
    if (magsign > 1)
      fMagSign = 1;
    else if ( magsign < 1)
      fMagSign = -1;
    else
      fMagSign = magsign;

    Double_t rad = fRaddedps;

    double afsi0b = 0.07510020733*chg1*fMagSign*rad/ptv1;
    double afsi1b = 0.07510020733*chg2*fMagSign*rad/ptv2;
    Double_t dps6 =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
    dps6 = TVector2::Phi_mpi_pi(dps6);
    double etad = eta2 - eta1;

    fDenDEtaDPhiS->Fill(dps6,etad);
  }
//_______________________________________________________________

}
//____________________________
void AliFemtoQinvCorrFctn::Write(){
  // Write out neccessary objects
  fNumerator->Write();
  fDenominator->Write();
  fkTMonitor->Write();
  if (fDetaDphiscal) {
    fNumDEtaDPhiS->Write();
    fDenDEtaDPhiS->Write();
  }
}
//______________________________
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
  return tOutputList;
}

void AliFemtoQinvCorrFctn::CalculateDetaDphis(Bool_t dedpsc, Double_t rad) {
  fDetaDphiscal = dedpsc;
  fRaddedps = rad;
}
