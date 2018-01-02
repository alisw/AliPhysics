////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelBPLCMSCorrFctnKK - the class for correlation function which   ///
/// uses the model framework and weight generation and calculated the 3D     ///
/// correlation function in the Bertsh-Pratt LCMS system                     ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoModelBPLCMSCorrFctnKK.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoAnalysisReactionPlane.h"
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoModelBPLCMSCorrFctnKK)
#endif

//____________________________
AliFemtoModelBPLCMSCorrFctnKK::AliFemtoModelBPLCMSCorrFctnKK(const char* title, const int& nbins, const float& QLo, const float& QHi)
  :
  AliFemtoModelCorrFctn(title, nbins, QLo, QHi),
  fNumerator3DTrue(0),
  fNumerator3DFake(0),
  fDenominator3D(0),
  fQinvHisto(0),
  fPairCut(0),
  fUseRPSelection(0),
  fNumerator3DTrueIdeal(0),
  fNumerator3DFakeIdeal(0),
  fDenominator3DIdeal(0)
{
  // set up true numerator
  char tTitNumT[101] = "Num3DTrue";
  strncat(tTitNumT,title, 100);
  fNumerator3DTrue = new TH3D(tTitNumT,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up fake numerator
  char tTitNumF[101] = "Num3DFake";
  strncat(tTitNumF,title, 100);
  fNumerator3DFake = new TH3D(tTitNumF,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up denominator
  char tTitDen[101] = "Den3D";
  strncat(tTitDen,title, 100);
  fDenominator3D = new TH3D(tTitDen,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up ave qInv
  char tTitQinv[101] = "Qinv";
  strncat(tTitQinv,title, 100);
  fQinvHisto = new TH3D(tTitQinv,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);

 // set up true Ideal numerator
  char tTitNumTI[101] = "Num3DTrueIdeal";
  strncat(tTitNumTI,title, 100);
  fNumerator3DTrueIdeal = new TH3D(tTitNumTI,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);

 // set up fake Ideal numerator
  char tTitNumFI[101] = "Num3DFakeIdeal";
  strncat(tTitNumFI,title, 100);
  fNumerator3DFakeIdeal = new TH3D(tTitNumFI,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);

// set up denominator
  char tTitDenI[101] = "Den3DIdeal";
  strncat(tTitDenI,title, 100);
  fDenominator3DIdeal = new TH3D(tTitDenI,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);


  // to enable error bar calculation...
  fNumerator3DTrue->Sumw2();
  fNumerator3DFake->Sumw2();
  fDenominator3D->Sumw2();
  fNumerator3DTrueIdeal->Sumw2();
  fNumerator3DFakeIdeal->Sumw2();
  fDenominator3DIdeal->Sumw2();
}

AliFemtoModelBPLCMSCorrFctnKK::AliFemtoModelBPLCMSCorrFctnKK(const AliFemtoModelBPLCMSCorrFctnKK& aCorrFctn) :
  AliFemtoModelCorrFctn(aCorrFctn),
  fNumerator3DTrue(0),
  fNumerator3DFake(0),
  fDenominator3D(0),
  fQinvHisto(0),
  fPairCut(0),
  fUseRPSelection(0),
  fNumerator3DTrueIdeal(0),
  fNumerator3DFakeIdeal(0),
  fDenominator3DIdeal(0)
{
  // Copy constructor
  fNumerator3DTrue = new TH3D(*aCorrFctn.fNumerator3DTrue);
  fNumerator3DFake = new TH3D(*aCorrFctn.fNumerator3DFake);
  fDenominator3D   = new TH3D(*aCorrFctn.fDenominator3D);
  fNumerator3DTrueIdeal = new TH3D(*aCorrFctn.fNumerator3DTrueIdeal);
  fNumerator3DFakeIdeal = new TH3D(*aCorrFctn.fNumerator3DFakeIdeal);
  fDenominator3DIdeal   = new TH3D(*aCorrFctn.fDenominator3DIdeal);
  fQinvHisto       = new TH3D(*aCorrFctn.fQinvHisto);
  fPairCut         = aCorrFctn.fPairCut->Clone();
}
//____________________________
AliFemtoModelBPLCMSCorrFctnKK::~AliFemtoModelBPLCMSCorrFctnKK()
{
  // destructor
  if (fNumeratorTrue) delete fNumeratorTrue;
  if (fNumeratorFake) delete fNumeratorFake;
  if (fDenominator) delete fDenominator;
  delete fNumerator3DTrue;
  delete fNumerator3DFake;
  delete fDenominator3D;
  delete fNumerator3DTrueIdeal;
  delete fNumerator3DFakeIdeal;
  delete fDenominator3DIdeal;
  delete fQinvHisto;
  if (fPairCut) delete fPairCut;
}
//_________________________
AliFemtoModelBPLCMSCorrFctnKK& AliFemtoModelBPLCMSCorrFctnKK::operator=(const AliFemtoModelBPLCMSCorrFctnKK& aCorrFctn)
{
  // Assignment operator
  if (this == &aCorrFctn)
    return *this;
  if (fNumerator3DTrue) delete fNumerator3DTrue;
  fNumerator3DTrue = new TH3D(*aCorrFctn.fNumerator3DTrue);
  if (fNumerator3DFake) delete fNumerator3DFake;
  fNumerator3DFake = new TH3D(*aCorrFctn.fNumerator3DFake);
  if (fDenominator3D) delete fDenominator3D;
  fDenominator3D = new TH3D(*aCorrFctn.fDenominator3D);

 if (fNumerator3DTrueIdeal) delete fNumerator3DTrueIdeal;
  fNumerator3DTrueIdeal = new TH3D(*aCorrFctn.fNumerator3DTrueIdeal);
 if (fNumerator3DFakeIdeal) delete fNumerator3DFakeIdeal;
  fNumerator3DFakeIdeal = new TH3D(*aCorrFctn.fNumerator3DFakeIdeal);
if (fDenominator3DIdeal) delete fDenominator3DIdeal;
  fDenominator3DIdeal = new TH3D(*aCorrFctn.fDenominator3DIdeal);
  if (fQinvHisto) delete fQinvHisto;
  fQinvHisto = new TH3D(*aCorrFctn.fQinvHisto);
  fPairCut = aCorrFctn.fPairCut->Clone();

  return *this;
}

//_________________________
void AliFemtoModelBPLCMSCorrFctnKK::Write(){
  // Write out data histograms
  AliFemtoModelCorrFctn::Write();
  fNumerator3DTrue->Write();
  fNumerator3DFake->Write();
  fDenominator3D->Write();
  fNumerator3DTrueIdeal->Write();
  fNumerator3DFakeIdeal->Write();
  fDenominator3DIdeal->Write();
  fQinvHisto->Write();
}
//________________________
TList* AliFemtoModelBPLCMSCorrFctnKK::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = AliFemtoModelCorrFctn::GetOutputList();

  tOutputList->Add(fNumerator3DTrue);
  tOutputList->Add(fNumerator3DFake);
  tOutputList->Add(fDenominator3D);
  tOutputList->Add(fNumerator3DTrueIdeal);
  tOutputList->Add(fNumerator3DFakeIdeal);
  tOutputList->Add(fDenominator3DIdeal);
  tOutputList->Add(fQinvHisto);

  return tOutputList;
}

//_________________________
void AliFemtoModelBPLCMSCorrFctnKK::Finish(){
  fQinvHisto->Divide(fDenominator);
}

//____________________________
AliFemtoString AliFemtoModelBPLCMSCorrFctnKK::Report(){
  // Prepare a report from the execution
  string stemp = "LCMS Frame Bertsch-Pratt 3D Model Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fNumeratorTrue->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fDenominator->GetEntries());
  stemp += ctemp;

  /*  if (fCorrection)
      {
      float radius = fCorrection->GetRadius();
      snprintf(ctemp , 100, "Coulomb correction used radius of\t%E\n",radius);
      }
      else
      {
      snprintf(ctemp , 100, "No Coulomb Correction applied to this CorrFctn\n");
      }
      stemp += ctemp;
  */

  //
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoModelBPLCMSCorrFctnKK::AddRealPair( AliFemtoPair* pair)
{
  // Store a real pair in numerator
  if (fPairCut){
    if (fUseRPSelection) {
      AliFemtoKTPairCut *ktc = dynamic_cast<AliFemtoKTPairCut *>(fPairCut);
      if (!ktc) {
	cout << "RP aware cut requested, but not connected to the CF" << endl;
	if (!(fPairCut->Pass(pair))) return;
      }
      else {
	AliFemtoAnalysisReactionPlane *arp = dynamic_cast<AliFemtoAnalysisReactionPlane *> (HbtAnalysis());
	if (!arp) {
	  cout << "RP aware cut requested, but not connected to the CF" << endl;
	  if (!(fPairCut->Pass(pair))) return;
	}
	else if (!(ktc->Pass(pair, arp->GetCurrentReactionPlane()))) return;
      }
    }
    else
      if (!(fPairCut->Pass(pair))) return;
  }
//   if (fPairCut){
//     if (!(fPairCut->Pass(pair))) return;
//   }

  Double_t weight = fManager->GetWeight(pair);

  double qOut = (pair->QOutCMS());
  double qSide = (pair->QSideCMS());
  double qLong = (pair->QLongCMS());

  double qOutTrue = GetQoutTrue(pair);
  double qSideTrue = GetQsideTrue(pair);
  double qLongTrue =  GetQlongTrue(pair);


  fNumerator3DTrue->Fill(qOut, qSide, qLong, weight);
  fNumeratorTrue->Fill(pair->QInv(), weight);

  fNumerator3DTrueIdeal->Fill(qOutTrue, qSideTrue, qLongTrue, weight);

}
//____________________________
void AliFemtoModelBPLCMSCorrFctnKK::AddMixedPair( AliFemtoPair* pair){
  // store mixed pair in denominator
  if (fPairCut){
    if (fUseRPSelection) {
      AliFemtoKTPairCut *ktc = dynamic_cast<AliFemtoKTPairCut *>(fPairCut);
      if (!ktc) {
	cout << "RP aware cut requested, but not connected to the CF" << endl;
	if (!(fPairCut->Pass(pair))) return;
      }
      else {
	AliFemtoAnalysisReactionPlane *arp = dynamic_cast<AliFemtoAnalysisReactionPlane *> (HbtAnalysis());
	if (!arp) {
	  cout << "RP aware cut requested, but not connected to the CF" << endl;
	  if (!(fPairCut->Pass(pair))) return;
	}
	else if (!(ktc->Pass(pair, arp->GetCurrentReactionPlane()))) return;
      }
    }
    else
      if (!(fPairCut->Pass(pair))) return;
  }
//   if (fPairCut){
//     if (!(fPairCut->Pass(pair))) return;
//   }

  Double_t weight = fManager->GetWeight(pair);

  double qOut = (pair->QOutCMS());
  double qSide = (pair->QSideCMS());
  double qLong = (pair->QLongCMS());

  double qOutTrue = GetQoutTrue(pair);
  double qSideTrue = GetQsideTrue(pair);
  double qLongTrue =  GetQlongTrue(pair);

  fNumerator3DFake->Fill(qOut, qSide, qLong, weight);
  fDenominator3D->Fill(qOut, qSide, qLong, 1.0);
  fNumeratorFake->Fill(pair->QInv(), weight);
  fDenominator->Fill(pair->QInv(), 1.0);

  fNumerator3DFakeIdeal->Fill(qOutTrue, qSideTrue, qLongTrue, weight);
  fDenominator3DIdeal->Fill(qOutTrue, qSideTrue, qLongTrue, 1.0);

}
//_______________________
AliFemtoModelCorrFctn* AliFemtoModelBPLCMSCorrFctnKK::Clone()
{
  // Clone the correlation function
  AliFemtoModelBPLCMSCorrFctnKK *tCopy = new AliFemtoModelBPLCMSCorrFctnKK(*this);

  return tCopy;
}

void AliFemtoModelBPLCMSCorrFctnKK::SetSpecificPairCut(AliFemtoPairCut* aCut)
{
  fPairCut = aCut;
}

void AliFemtoModelBPLCMSCorrFctnKK::SetUseRPSelection(unsigned short aRPSel)
{
  fUseRPSelection = aRPSel;
}


//_______________________
Double_t AliFemtoModelBPLCMSCorrFctnKK::GetQinvTrue(AliFemtoPair* aPair)
{
  AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
  AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();

  AliFemtoLorentzVector fm1;
  AliFemtoThreeVector* temp = inf1->GetTrueMomentum();
  fm1.SetVect(*temp);
  double ener = TMath::Sqrt(temp->Mag2()+(aPair->Track1()->Track()->GetMass())*(aPair->Track1()->Track()->GetMass()));
  fm1.SetE(ener);

  AliFemtoLorentzVector fm2;
  AliFemtoThreeVector* temp2 = inf2->GetTrueMomentum();
  fm2.SetVect(*temp2);
  ener = TMath::Sqrt(temp2->Mag2()+(aPair->Track2()->Track()->GetMass())*(aPair->Track2()->Track()->GetMass()));
  fm2.SetE(ener);

  AliFemtoLorentzVector tQinvTrueVec = (fm1-fm2);
  Double_t tQinvTrue = -1.* tQinvTrueVec.m();

  return tQinvTrue;
}

//_______________________
Double_t AliFemtoModelBPLCMSCorrFctnKK::GetQoutTrue(AliFemtoPair* aPair)
{

 // relative momentum out component in lab frame
  AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
  AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();

  // AliFemtoLorentzVector fm1;
  AliFemtoThreeVector* tmp1 = inf1->GetTrueMomentum();
  /*
 fm1.SetVect(*tmp1);
  double ener = TMath::Sqrt(tmp1->Mag2()+(aPair->Track1()->Track()->GetMass())*(aPair->Track1()->Track()->GetMass()));
  fm1.SetE(ener);
  */
  //  AliFemtoLorentzVector fm2;
  AliFemtoThreeVector* tmp2 = inf2->GetTrueMomentum();
  /*
  fm2.SetVect(*tmp2);
  ener = TMath::Sqrt(tmp2->Mag2()+(aPair->Track2()->Track()->GetMass())*(aPair->Track2()->Track()->GetMass()));
  fm2.SetE(ener);
  */
    double dx = tmp1->x() - tmp2->x();
    double xt = tmp1->x() + tmp2->x();

    double dy = tmp1->y() - tmp2->y();
    double yt = tmp1->y() + tmp2->y();

    double k1 = (::sqrt(xt*xt+yt*yt));
    double k2 = (dx*xt+dy*yt);
    double tmp;

    if(k1!=0) tmp= k2/k1;
    else tmp=0;

    return (tmp);

}



//_________________
Double_t AliFemtoModelBPLCMSCorrFctnKK::GetQsideTrue(AliFemtoPair* aPair)
{

  AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
  AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();
  AliFemtoThreeVector* tmp1 = inf1->GetTrueMomentum();
  AliFemtoThreeVector* tmp2 = inf2->GetTrueMomentum();

  // relative momentum side component in lab frame

    double x1 = tmp1->x();  double y1 = tmp1->y();
    double x2 = tmp2->x();  double y2 = tmp2->y();

    double xt = x1+x2;  double yt = y1+y2;
    double k1 = ::sqrt(xt*xt+yt*yt);

    double tmp;
    if(k1!=0) tmp= 2.0*(x2*y1-x1*y2)/k1;
    else tmp=0;

    return (tmp);
}

//_________________________
double AliFemtoModelBPLCMSCorrFctnKK::GetQlongTrue(AliFemtoPair* aPair)
{
  // relative momentum component in lab frame

  AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
  AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();
  AliFemtoThreeVector* temp1 = inf1->GetTrueMomentum();
  AliFemtoThreeVector* temp2 = inf2->GetTrueMomentum();

  AliFemtoLorentzVector tmp1;
  tmp1.SetVect(*temp1);
 double ener = TMath::Sqrt(temp1->Mag2()+(aPair->Track1()->Track()->GetMass())*(aPair->Track1()->Track()->GetMass()));
  tmp1.SetE(ener);

AliFemtoLorentzVector tmp2;
  tmp2.SetVect(*temp2);
 double ener2 = TMath::Sqrt(temp2->Mag2()+(aPair->Track2()->Track()->GetMass())*(aPair->Track2()->Track()->GetMass()));
  tmp2.SetE(ener2);

    double dz = tmp1.z() - tmp2.z();
    double zz = tmp1.z() + tmp2.z();

    double dt = tmp1.t() - tmp2.t();
    double tt = tmp1.t() + tmp2.t();

    double beta = zz/tt;
    double gamma = 1.0/TMath::Sqrt((1.-beta)*(1.+beta));

    double temp = gamma*(dz - beta*dt);
    return (temp);
}
