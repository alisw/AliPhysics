///
/// \file AliFemtoBPLCMS3DCorrFctnKK.cxx
///

#include "AliFemtoBPLCMS3DCorrFctnKK.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoAnalysisReactionPlane.h"
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoBPLCMS3DCorrFctnKK);
  /// \endcond
#endif

//____________________________
AliFemtoBPLCMS3DCorrFctnKK::AliFemtoBPLCMS3DCorrFctnKK(
  const char* title,
  const int& nbins,
  const float& QLo,
  const float& QHi
):
  AliFemtoCorrFctn(),
//   fIDNumHisto(0),
//   fIDDenHisto(0),
//   fIDRatHisto(0),
//   fSMNumHisto(0),
//   fSMDenHisto(0),
//   fSMRatHisto(0),
//   fCorrectionHisto(0),
//   fCorrCFHisto(0),
  fNumerator(NULL),
  fDenominator(NULL),
  fRatio(NULL),
  fQinvHisto(NULL),
  fLambda(0.6),
  fRout2(6.0 * 6.0),
  fRside2(6.0 * 6.0),
  fRlong2(7.0 * 7.0),
  fQinvNormLo((QHi - QLo) * 0.8),
  fQinvNormHi((QHi - QLo) * 0.8),
  fNumRealsNorm(0),
  fNumMixedNorm(0),
  fUseRPSelection(0)
{
  /// Basic constructor
  /// set some stuff...

  //  fCorrection = 0;  // pointer to Coulomb Correction object

  //  fSmearPair = 0; // no resolution correction unless user sets SmearPair

  // set up numerator
  char tTitNum[101] = "Num";
  strncat(tTitNum, title, 100);
  fNumerator = new TH3D(tTitNum, title, nbins, -QHi, QHi, nbins, -QHi, QHi, nbins, -QHi, QHi);
  // set up denominator
  char tTitDen[101] = "Den";
  strncat(tTitDen, title, 100);
  fDenominator = new TH3D(tTitDen, title, nbins, -QHi, QHi, nbins, -QHi, QHi, nbins, -QHi, QHi);
  // set up uncorrected denominator
  char tTitDenUncoul[101] = "DenNoCoul";
  strncat(tTitDenUncoul, title, 100);
  //  fUncorrectedDenominator = new TH3D(tTitDenUncoul,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
  // set up ratio
  char tTitRat[101] = "Rat";
  strncat(tTitRat, title, 100);
  fRatio = new TH3D(tTitRat, title, nbins, -QHi, QHi, nbins, -QHi, QHi, nbins, -QHi, QHi);
  // set up ave qInv
  char tTitQinv[101] = "Qinv";
  strncat(tTitQinv, title, 100);
  fQinvHisto = new TH3D(tTitQinv, title, nbins, -QHi, QHi, nbins, -QHi, QHi, nbins, -QHi, QHi);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  //  fUncorrectedDenominator->Sumw2();
  fRatio->Sumw2();

//   // Following histos are for the momentum resolution correction
//   // they are filled only if a AliFemtoSmear object is plugged in
//   // here comes the "idea" numerator and denominator and ratio...
//   char tTitNumID[101] = "IDNum";
//   strncat(tTitNumID,title, 100);
//   fIDNumHisto = new TH3D(tTitNumID,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   char tTitDenID[101] = "IDDen";
//   strncat(tTitDenID,title, 100);
//   fIDDenHisto = new TH3D(tTitDenID,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   char tTitRatID[101] = "IDRat";
//   strncat(tTitRatID,title, 100);
//   fIDRatHisto = new TH3D(tTitRatID,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);

//   fIDNumHisto->Sumw2();
//   fIDDenHisto->Sumw2();
//   fIDRatHisto->Sumw2();

//   //
//   // here comes the "smeared" numerator and denominator...
//   char tTitNumSM[101] = "SMNum";
//   strncat(tTitNumSM,title, 100);
//   fSMNumHisto = new TH3D(tTitNumSM,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   char tTitDenSM[101] = "SMDen";
//   strncat(tTitDenSM,title, 100);
//   fSMDenHisto = new TH3D(tTitDenSM,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   char tTitRatSM[101] = "SMRat";
//   strncat(tTitRatSM,title, 100);
//   fSMRatHisto = new TH3D(tTitRatSM,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   //
//   fSMNumHisto->Sumw2();
//   fSMDenHisto->Sumw2();
//   fSMRatHisto->Sumw2();
//   //
//   // here comes the correction factor (which is just ratio of ideal ratio to smeared ratio)
//   char tTitCorrection[101] = "CorrectionFactor";
//   strncat(tTitCorrection,title, 100);
//   fCorrectionHisto = new TH3D(tTitCorrection,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   fCorrectionHisto->Sumw2();
//   // here comes the fully corrected correlation function
//   char tTitCorrCF[101] = "CorrectedCF";
//   strncat(tTitCorrCF,title, 100);
//   fCorrCFHisto = new TH3D(tTitCorrCF,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   fCorrCFHisto->Sumw2();

}

AliFemtoBPLCMS3DCorrFctnKK::AliFemtoBPLCMS3DCorrFctnKK(const AliFemtoBPLCMS3DCorrFctnKK& aCorrFctn):
  AliFemtoCorrFctn(aCorrFctn),
//   fIDNumHisto(0),
//   fIDDenHisto(0),
//   fIDRatHisto(0),
//   fSMNumHisto(0),
//   fSMDenHisto(0),
//   fSMRatHisto(0),
//   fCorrectionHisto(0),
//   fCorrCFHisto(0),
  fNumerator(NULL),
  fDenominator(NULL),
  fRatio(NULL),
  fQinvHisto(NULL),
  fLambda(aCorrFctn.fLambda),
  fRout2(aCorrFctn.fRout2),
  fRside2(aCorrFctn.fRside2),
  fRlong2(aCorrFctn.fRlong2),
  fQinvNormLo(aCorrFctn.fQinvNormLo),
  fQinvNormHi(aCorrFctn.fQinvNormHi),
  fNumRealsNorm(aCorrFctn.fNumRealsNorm),
  fNumMixedNorm(aCorrFctn.fNumMixedNorm),
  fUseRPSelection(aCorrFctn.fUseRPSelection)
{
/// Copy constructor
///   fIDNumHisto = new TH3D(*aCorrFctn.fIDNumHisto);
///   fIDDenHisto = new TH3D(*aCorrFctn.fIDDenHisto);
///   fIDRatHisto = new TH3D(*aCorrFctn.fIDRatHisto);
///   fSMNumHisto = new TH3D(*aCorrFctn.fSMNumHisto);
///   fSMDenHisto = new TH3D(*aCorrFctn.fSMDenHisto);
///   fSMRatHisto = new TH3D(*aCorrFctn.fSMRatHisto);
///   fCorrectionHisto = new TH3D(*aCorrFctn.fCorrectionHisto);
///   fCorrCFHisto = new TH3D(*aCorrFctn.fCorrCFHisto);

  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  fDenominator = new TH3D(*aCorrFctn.fDenominator);
  fRatio = new TH3D(*aCorrFctn.fRatio);
  fQinvHisto = new TH3D(*aCorrFctn.fQinvHisto);
}
//____________________________
AliFemtoBPLCMS3DCorrFctnKK::~AliFemtoBPLCMS3DCorrFctnKK()
{
  /// Destructor

  delete fNumerator;
  delete fDenominator;
  delete fRatio;
  delete fQinvHisto;
//   delete fIDNumHisto;
//   delete fIDDenHisto;
//   delete fIDRatHisto;
//   delete fSMNumHisto;
//   delete fSMDenHisto;
//   delete fSMRatHisto;
//   delete fCorrectionHisto;
//   delete fCorrCFHisto;
}
//_________________________
AliFemtoBPLCMS3DCorrFctnKK& AliFemtoBPLCMS3DCorrFctnKK::operator=(const AliFemtoBPLCMS3DCorrFctnKK& aCorrFctn)
{
  /// assignment operator

  if (this == &aCorrFctn)
    return *this;
//   if (fIDNumHisto) delete fIDNumHisto;
//   fIDNumHisto = new TH3D(*aCorrFctn.fIDNumHisto);
//   if (fIDDenHisto) delete fIDDenHisto;
//   fIDDenHisto = new TH3D(*aCorrFctn.fIDDenHisto);
//   if (fIDRatHisto) delete fIDRatHisto;
//   fIDRatHisto = new TH3D(*aCorrFctn.fIDRatHisto);
//   if (fSMNumHisto) delete fSMNumHisto;
//   fSMNumHisto = new TH3D(*aCorrFctn.fSMNumHisto);
//   if (fSMDenHisto) delete fSMDenHisto;
//   fSMDenHisto = new TH3D(*aCorrFctn.fSMDenHisto);
//   if (fSMRatHisto) delete fSMRatHisto;
//   fSMRatHisto = new TH3D(*aCorrFctn.fSMRatHisto);

//   if (fCorrectionHisto) delete fCorrectionHisto;
//   fCorrectionHisto = new TH3D(*aCorrFctn.fCorrectionHisto);
//   if (fCorrCFHisto) delete fCorrCFHisto;
//   fCorrCFHisto = new TH3D(*aCorrFctn.fCorrCFHisto);
  if (fNumerator) delete fNumerator;
  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  if (fDenominator) delete fDenominator;
  fDenominator = new TH3D(*aCorrFctn.fDenominator);
  if (fRatio) delete fRatio;
  fRatio = new TH3D(*aCorrFctn.fRatio);
  if (fQinvHisto) delete fQinvHisto;
  fQinvHisto = new TH3D(*aCorrFctn.fQinvHisto);

  fLambda = aCorrFctn.fLambda;
  fRout2 = aCorrFctn.fRout2;
  fRside2 = aCorrFctn.fRside2;
  fRlong2 = aCorrFctn.fRlong2;
  fQinvNormLo = aCorrFctn.fQinvNormLo;
  fQinvNormHi = aCorrFctn.fQinvNormHi;
  fNumRealsNorm = aCorrFctn.fNumRealsNorm;
  fNumMixedNorm = aCorrFctn.fNumMixedNorm;
  fUseRPSelection = aCorrFctn.fUseRPSelection;

  return *this;
}

//_________________________
void AliFemtoBPLCMS3DCorrFctnKK::WriteOutHistos()
{
  /// Write out all histograms to file

  fNumerator->Write();
  fDenominator->Write();
  //  fUncorrectedDenominator->Write();
  //ml fRatio->Write();
  fQinvHisto->Write();

  /*
    if (fSmearPair){
    fIDNumHisto->Write();
    fIDDenHisto->Write();
    fIDRatHisto->Write();
    //
    fSMNumHisto->Write();
    fSMDenHisto->Write();
    fSMRatHisto->Write();
    //
    fCorrectionHisto->Write();
    fCorrCFHisto->Write();
    }
  */
}
//______________________________
TList* AliFemtoBPLCMS3DCorrFctnKK::GetOutputList()
{
  /// Prepare the list of objects to be written to the output

  TList *tOutputList = new TList();

  tOutputList->Add(fNumerator);
  tOutputList->Add(fDenominator);
  tOutputList->Add(fQinvHisto);

  return tOutputList;
}

//_________________________
void AliFemtoBPLCMS3DCorrFctnKK::Finish()
{
  /// here is where we should normalize, fit, etc...

  double tNumFact, tDenFact;
  if ((fNumRealsNorm != 0) && (fNumMixedNorm != 0)) {
    tNumFact = double(fNumRealsNorm);
    tDenFact = double(fNumMixedNorm);
  }
  // can happen that the fNumRealsNorm and fNumMixedNorm = 0 if you do non-standard
  //   things like making a new CorrFctn and just setting the Numerator and Denominator
  //   from OTHER CorrFctns which you read in (like when doing parallel processing)
  //ml else{
  //ml    cout << "Warning! - no normalization constants defined - I do the best I can..." << endl;
  //ml    int nbins = fNumerator->GetNbinsX();
  //ml   int half_way = nbins/2;
  //ml   tNumFact = fNumerator->Integral(half_way,nbins,half_way,nbins,half_way,nbins);
  //ml  tDenFact = fDenominator->Integral(half_way,nbins,half_way,nbins,half_way,nbins);
  //ml    }

  ////ml  fRatio->Divide(fNumerator,fDenominator,tDenFact,tNumFact);
  //  fQinvHisto->Divide(fUncorrectedDenominator);
  ///ml fQinvHisto->Divide(fDenominator);

  /*
  // now do all the resolution correction stuff..
  if (fSmearPair){  // but only do it if we have been working with a SmearPair
  fIDRatHisto->Divide(fIDNumHisto,fIDDenHisto);
  fSMRatHisto->Divide(fSMNumHisto,fSMDenHisto);
  fCorrectionHisto->Divide(fIDRatHisto,fSMRatHisto);
  fCorrCFHisto->Multiply(fRatio,fCorrectionHisto);
  }
  */

}

//____________________________
AliFemtoString AliFemtoBPLCMS3DCorrFctnKK::Report()
{
  /// Construct the report
  AliFemtoString report("LCMS Frame Bertsch-Pratt 3D Correlation Function Report:\n");

  report += Form("Number of entries in numerator:\t%E\n", fNumerator->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n", fDenominator->GetEntries());
  report += Form("Number of entries in ratio:\t%E\n", fRatio->GetEntries());
  report += Form("Normalization region in Qinv was:\t%E\t%E\n", fQinvNormLo, fQinvNormHi);
  report += Form("Number of pairs in Normalization region was:\n");
  report += Form("In numerator:\t%lu\t In denominator:\t%lu\n", fNumRealsNorm, fNumMixedNorm);

  /*  if (fCorrection)
      {
      report += Form("Coulomb correction used radius of\t%E\n",fCorrection->GetRadius());
      }
      else
      {
      report += "No Coulomb Correction applied to this CorrFctn\n";
      }
  */

  if (fPairCut) {
    report += "Here is the PairCut specific to this CorrFctn\n";
    report += fPairCut->Report();
  } else {
    report += "No PairCut specific to this CorrFctn\n";
  }

  return report;
}
//____________________________
void AliFemtoBPLCMS3DCorrFctnKK::AddRealPair( AliFemtoPair* pair)
{
  /// perform operations on real pairs

  if (fPairCut) {
    if (fUseRPSelection) {
      AliFemtoKTPairCut *ktc = dynamic_cast<AliFemtoKTPairCut *>(fPairCut);
      if (!ktc) {
        cout << "RP aware cut requested, but not connected to the CF" << endl;
        if (!(fPairCut->Pass(pair))) return;
      } else {
        AliFemtoAnalysisReactionPlane *arp = dynamic_cast<AliFemtoAnalysisReactionPlane *> (HbtAnalysis());
        if (!arp) {
          cout << "RP aware cut requested, but not connected to the CF" << endl;
          if (!(fPairCut->Pass(pair))) return;
        } else if (!(ktc->Pass(pair, arp->GetCurrentReactionPlane()))) return;
      }
    } else if (!(fPairCut->Pass(pair))) return;
  }

  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  if ((tQinv < fQinvNormHi) && (tQinv > fQinvNormLo)) {
    fNumRealsNorm++;
  }
  const double qOut = (pair->QOutCMS()),
               qSide = (pair->QSideCMS()),
               qLong = (pair->QLongCMS());

  fNumerator->Fill(qOut, qSide, qLong);
}
//____________________________
void AliFemtoBPLCMS3DCorrFctnKK::AddMixedPair( AliFemtoPair* pair)
{
/// perform operations on mixed pairs
///   if (fPairCut){
///     if (!(fPairCut->Pass(pair))) return;
///   }

  if (fPairCut) {
    if (fUseRPSelection) {
      AliFemtoKTPairCut *ktc = dynamic_cast<AliFemtoKTPairCut *>(fPairCut);
      if (!ktc) {
        cout << "RP aware cut requested, but not connected to the CF" << endl;
        if (!(fPairCut->Pass(pair))) return;
      } else {
        AliFemtoAnalysisReactionPlane *arp = dynamic_cast<AliFemtoAnalysisReactionPlane *> (HbtAnalysis());
        if (!arp) {
          cout << "RP aware cut requested, but not connected to the CF" << endl;
          if (!(fPairCut->Pass(pair))) return;
        } else if (!(ktc->Pass(pair, arp->GetCurrentReactionPlane()))) return;
      }
    } else if (!(fPairCut->Pass(pair))) return;
  }

  //  double CoulombWeight = (fCorrection ? fCorrection->CoulombCorrect(pair) : 1.0);
  double tCoulombWeight = 1.0;

  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  if ((tQinv < fQinvNormHi) && (tQinv > fQinvNormLo)) fNumMixedNorm++;
  double qOut = (pair->QOutCMS());
  double qSide = (pair->QSideCMS());
  double qLong = (pair->QLongCMS());

  fDenominator->Fill(qOut, qSide, qLong, tCoulombWeight);
  //  fUncorrectedDenominator->Fill(qOut,qSide,qLong,1.0);
  fQinvHisto->Fill(qOut, qSide, qLong, tQinv);

  /*
  // now for the momentum resolution stuff...
  if (fSmearPair){
      double CorrWeight =  1.0 +
      fLambda*exp((-qOut*qOut*fRout2 -qSide*qSide*fRside2 -qLong*qLong*fRlong2)/0.038936366329);
    CorrWeight *= CoulombWeight;  // impt.

    fIDNumHisto->Fill(qOut,qSide,qLong,CorrWeight);
    fIDDenHisto->Fill(qOut,qSide,qLong,CoulombWeight);

    fSmearPair->SetUnsmearedPair(pair);
    double qOut_prime = fabs(fSmearPair->SmearedPair().qOutCMS());
    double qSide_prime = fabs(fSmearPair->SmearedPair().qSideCMS());
    double qLong_prime = fabs(fSmearPair->SmearedPair().qLongCMS());

    fSMNumHisto->Fill(qOut_prime,qSide_prime,qLong_prime,CorrWeight);

    double SmearedCoulombWeight = ( fCorrection ?
            fCorrection->CoulombCorrect(&(fSmearPair->SmearedPair())) :
            1.0);

    fSMDenHisto->Fill(qOut_prime,qSide_prime,qLong_prime,SmearedCoulombWeight);
  }
  */
}


void AliFemtoBPLCMS3DCorrFctnKK::SetUseRPSelection(unsigned short aRPSel)
{
  fUseRPSelection = aRPSel;
}
