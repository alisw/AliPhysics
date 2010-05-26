///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoBPLCMS3DCorrFctn: a class to calculate 3D correlation         //
// for pairs of identical particles.                                     //
// It also stored the weighted qinv per bin histogram for the coulomb    //
// correction.                                                           //
// In analysis the function should be first created in a macro, then     //
// added to the analysis, and at the end of the macro the procedure to   //
// write out histograms should be called.                                //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliFemtoBPLCMS3DCorrFctn.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoAnalysisReactionPlane.h"
//#include "AliFemtoHisto.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoBPLCMS3DCorrFctn)
#endif

//____________________________
AliFemtoBPLCMS3DCorrFctn::AliFemtoBPLCMS3DCorrFctn(char* title, const int& nbins, const float& QLo, const float& QHi)
  :
  AliFemtoCorrFctn(),
//   fIDNumHisto(0),
//   fIDDenHisto(0),
//   fIDRatHisto(0),
//   fSMNumHisto(0),
//   fSMDenHisto(0),
//   fSMRatHisto(0),
//   fCorrectionHisto(0),
//   fCorrCFHisto(0),
  fNumerator(0),
  fDenominator(0),
  fRatio(0),
  fQinvHisto(0),
  fLambda(0),
  fRout2(0),
  fRside2(0),
  fRlong2(0),
  fQinvNormLo(0),
  fQinvNormHi(0),
  fNumRealsNorm(0),
  fNumMixedNorm(0),
  fUseRPSelection(0)
{
  // Basic constructor
  // set some stuff...
  fQinvNormLo = (QHi-QLo)*0.8;
  fQinvNormHi = (QHi-QLo)*0.8;
  fNumRealsNorm = 0;
  fNumMixedNorm = 0;
  //  fCorrection = 0;  // pointer to Coulomb Correction object

  //  fSmearPair = 0; // no resolution correction unless user sets SmearPair

  // set up numerator
  char tTitNum[100] = "Num";
  strcat(tTitNum,title);
  fNumerator = new TH3D(tTitNum,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
  // set up denominator
  char tTitDen[100] = "Den";
  strcat(tTitDen,title);
  fDenominator = new TH3D(tTitDen,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
  // set up uncorrected denominator
  char tTitDenUncoul[100] = "DenNoCoul";
  strcat(tTitDenUncoul,title);
  //  fUncorrectedDenominator = new TH3D(tTitDenUncoul,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
  // set up ratio
  char tTitRat[100] = "Rat";
  strcat(tTitRat,title);
  fRatio = new TH3D(tTitRat,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
  // set up ave qInv
  char tTitQinv[100] = "Qinv";
  strcat(tTitQinv,title);
  fQinvHisto = new TH3D(tTitQinv,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  //  fUncorrectedDenominator->Sumw2();
  fRatio->Sumw2();

//   // Following histos are for the momentum resolution correction
//   // they are filled only if a AliFemtoSmear object is plugged in
//   // here comes the "idea" numerator and denominator and ratio...
//   char tTitNumID[100] = "IDNum";
//   strcat(tTitNumID,title);
//   fIDNumHisto = new TH3D(tTitNumID,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   char tTitDenID[100] = "IDDen";
//   strcat(tTitDenID,title);
//   fIDDenHisto = new TH3D(tTitDenID,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   char tTitRatID[100] = "IDRat";
//   strcat(tTitRatID,title);
//   fIDRatHisto = new TH3D(tTitRatID,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);

//   fIDNumHisto->Sumw2();
//   fIDDenHisto->Sumw2();
//   fIDRatHisto->Sumw2();

//   //
//   // here comes the "smeared" numerator and denominator...
//   char tTitNumSM[100] = "SMNum";
//   strcat(tTitNumSM,title);
//   fSMNumHisto = new TH3D(tTitNumSM,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   char tTitDenSM[100] = "SMDen";
//   strcat(tTitDenSM,title);
//   fSMDenHisto = new TH3D(tTitDenSM,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   char tTitRatSM[100] = "SMRat";
//   strcat(tTitRatSM,title);
//   fSMRatHisto = new TH3D(tTitRatSM,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   //
//   fSMNumHisto->Sumw2();
//   fSMDenHisto->Sumw2();
//   fSMRatHisto->Sumw2();
//   //
//   // here comes the correction factor (which is just ratio of ideal ratio to smeared ratio)
//   char tTitCorrection[100] = "CorrectionFactor";
//   strcat(tTitCorrection,title);
//   fCorrectionHisto = new TH3D(tTitCorrection,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);  
//   fCorrectionHisto->Sumw2();
//   // here comes the fully corrected correlation function
//   char tTitCorrCF[100] = "CorrectedCF";
//   strcat(tTitCorrCF,title);
//   fCorrCFHisto = new TH3D(tTitCorrCF,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
//   fCorrCFHisto->Sumw2();

  // user can (and should) override these defaults...
  fLambda = 0.6;
  fRout2 = 6.0*6.0;
  fRside2 = 6.0*6.0;
  fRlong2 = 7.0*7.0;

}

AliFemtoBPLCMS3DCorrFctn::AliFemtoBPLCMS3DCorrFctn(const AliFemtoBPLCMS3DCorrFctn& aCorrFctn) :
  AliFemtoCorrFctn(aCorrFctn),
//   fIDNumHisto(0),
//   fIDDenHisto(0),
//   fIDRatHisto(0),
//   fSMNumHisto(0),
//   fSMDenHisto(0),
//   fSMRatHisto(0),
//   fCorrectionHisto(0),
//   fCorrCFHisto(0),
  fNumerator(0),
  fDenominator(0),
  fRatio(0),
  fQinvHisto(0),
  fLambda(0),
  fRout2(0),
  fRside2(0),
  fRlong2(0),
  fQinvNormLo(0),
  fQinvNormHi(0),
  fNumRealsNorm(0),
  fNumMixedNorm(0),
  fUseRPSelection(0)
{
  // Copy constructor
//   fIDNumHisto = new TH3D(*aCorrFctn.fIDNumHisto);
//   fIDDenHisto = new TH3D(*aCorrFctn.fIDDenHisto);
//   fIDRatHisto = new TH3D(*aCorrFctn.fIDRatHisto);
//   fSMNumHisto = new TH3D(*aCorrFctn.fSMNumHisto);
//   fSMDenHisto = new TH3D(*aCorrFctn.fSMDenHisto);
//   fSMRatHisto = new TH3D(*aCorrFctn.fSMRatHisto);
//   fCorrectionHisto = new TH3D(*aCorrFctn.fCorrectionHisto);
//   fCorrCFHisto = new TH3D(*aCorrFctn.fCorrCFHisto);
  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  fDenominator = new TH3D(*aCorrFctn.fDenominator);
  fRatio = new TH3D(*aCorrFctn.fRatio);
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
}
//____________________________
AliFemtoBPLCMS3DCorrFctn::~AliFemtoBPLCMS3DCorrFctn(){
  // Destructor
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
AliFemtoBPLCMS3DCorrFctn& AliFemtoBPLCMS3DCorrFctn::operator=(const AliFemtoBPLCMS3DCorrFctn& aCorrFctn)
{
  // assignment operator
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
void AliFemtoBPLCMS3DCorrFctn::WriteOutHistos(){
  // Write out all histograms to file
  fNumerator->Write();
  fDenominator->Write();
  //  fUncorrectedDenominator->Write();
  fRatio->Write();
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
TList* AliFemtoBPLCMS3DCorrFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumerator); 
  tOutputList->Add(fDenominator);  
  tOutputList->Add(fQinvHisto);  

  return tOutputList;
}

//_________________________
void AliFemtoBPLCMS3DCorrFctn::Finish(){
  // here is where we should normalize, fit, etc...
  double tNumFact,tDenFact;
  if ((fNumRealsNorm !=0) && (fNumMixedNorm !=0)){
    tNumFact = double(fNumRealsNorm);
    tDenFact = double(fNumMixedNorm);
  }
  // can happen that the fNumRealsNorm and fNumMixedNorm = 0 if you do non-standard
  //   things like making a new CorrFctn and just setting the Numerator and Denominator
  //   from OTHER CorrFctns which you read in (like when doing parallel processing) 
  else{
    cout << "Warning! - no normalization constants defined - I do the best I can..." << endl;
    int nbins = fNumerator->GetNbinsX();
    int half_way = nbins/2;
    tNumFact = fNumerator->Integral(half_way,nbins,half_way,nbins,half_way,nbins);
    tDenFact = fDenominator->Integral(half_way,nbins,half_way,nbins,half_way,nbins);
  }

  fRatio->Divide(fNumerator,fDenominator,tDenFact,tNumFact);
  //  fQinvHisto->Divide(fUncorrectedDenominator);
  fQinvHisto->Divide(fDenominator);

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
AliFemtoString AliFemtoBPLCMS3DCorrFctn::Report(){
  // Construct the report
  string stemp = "LCMS Frame Bertsch-Pratt 3D Correlation Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerator:\t%E\n",fNumerator->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominator:\t%E\n",fDenominator->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in ratio:\t%E\n",fRatio->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Normalization region in Qinv was:\t%E\t%E\n",fQinvNormLo,fQinvNormHi);
  stemp += ctemp;
  sprintf(ctemp,"Number of pairs in Normalization region was:\n");
  stemp += ctemp;
  sprintf(ctemp,"In numerator:\t%lu\t In denominator:\t%lu\n",fNumRealsNorm,fNumMixedNorm);
  stemp += ctemp;
  /*  if (fCorrection)
      {
      float radius = fCorrection->GetRadius();
      sprintf(ctemp,"Coulomb correction used radius of\t%E\n",radius);
      }
      else
      {
      sprintf(ctemp,"No Coulomb Correction applied to this CorrFctn\n");
      }
      stemp += ctemp;
  */

  if (fPairCut){
    sprintf(ctemp,"Here is the PairCut specific to this CorrFctn\n");
    stemp += ctemp;
    stemp += fPairCut->Report();
  }
  else{
    sprintf(ctemp,"No PairCut specific to this CorrFctn\n");
    stemp += ctemp;
  }

  //  
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoBPLCMS3DCorrFctn::AddRealPair( AliFemtoPair* pair){
  // perform operations on real pairs
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
	if (!(ktc->Pass(pair, arp->GetCurrentReactionPlane()))) return;
      }
    }
    else
      if (!(fPairCut->Pass(pair))) return;
  }

  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  if ((tQinv < fQinvNormHi) && (tQinv > fQinvNormLo)) fNumRealsNorm++;
  double qOut = (pair->QOutCMS());
  double qSide = (pair->QSideCMS());
  double qLong = (pair->QLongCMS());

  fNumerator->Fill(qOut,qSide,qLong);
}
//____________________________
void AliFemtoBPLCMS3DCorrFctn::AddMixedPair( AliFemtoPair* pair){
  // perform operations on mixed pairs
//   if (fPairCut){
//     if (!(fPairCut->Pass(pair))) return;
//   }
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
	if (!(ktc->Pass(pair, arp->GetCurrentReactionPlane()))) return;
      }
    }
    else
      if (!(fPairCut->Pass(pair))) return;
  }

  //  double CoulombWeight = (fCorrection ? fCorrection->CoulombCorrect(pair) : 1.0);
  double tCoulombWeight = 1.0;

  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  if ((tQinv < fQinvNormHi) && (tQinv > fQinvNormLo)) fNumMixedNorm++;
  double qOut = (pair->QOutCMS());
  double qSide = (pair->QSideCMS());
  double qLong = (pair->QLongCMS());

  fDenominator->Fill(qOut,qSide,qLong,tCoulombWeight);
  //  fUncorrectedDenominator->Fill(qOut,qSide,qLong,1.0);
  fQinvHisto->Fill(qOut,qSide,qLong,tQinv);

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


void AliFemtoBPLCMS3DCorrFctn::SetUseRPSelection(unsigned short aRPSel)
{
  fUseRPSelection = aRPSel;
}
