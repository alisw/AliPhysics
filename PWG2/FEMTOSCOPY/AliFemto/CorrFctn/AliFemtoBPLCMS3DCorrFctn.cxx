/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   3D Bertsch-Pratt decomposition in the LCMS frame
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.7  2003/01/31 19:20:54  magestro
 * Cleared up simple compiler warnings on i386_linux24
 *
 * Revision 1.6  2002/06/07 22:51:38  lisa
 * Widely used AliFemtoBPLCMS3DCorrFctn class now accumulates UNcorrected denominator and has a WriteOutHistos method
 *
 * Revision 1.5  2001/05/23 00:19:04  lisa
 * Add in Smearing classes and methods needed for momentum resolution studies and correction
 *
 * Revision 1.4  2000/10/26 19:48:49  rcwells
 * Added functionality for Coulomb correction of <qInv> in 3D correltions
 *
 * Revision 1.3  2000/09/14 18:36:53  lisa
 * Added Qinv and ExitSep pair cuts and AliFemtoBPLCMS3DCorrFctn_SIM CorrFctn
 *
 * Revision 1.2  2000/08/23 19:43:43  lisa
 * added alternate normalization algorithm to 3d CorrFctns in case normal one fails
 *
 * Revision 1.1  2000/08/17 20:48:39  lisa
 * Adding correlationfunction in LCMS frame
 *
 *
 **************************************************************************/

#include "CorrFctn/AliFemtoBPLCMS3DCorrFctn.h"
//#include "Infrastructure/AliFemtoHisto.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoBPLCMS3DCorrFctn)
#endif

//____________________________
AliFemtoBPLCMS3DCorrFctn::AliFemtoBPLCMS3DCorrFctn(char* title, const int& nbins, const float& QLo, const float& QHi)
  :
  fIDNumHisto(0),
  fIDDenHisto(0),
  fIDRatHisto(0),
  fSMNumHisto(0),
  fSMDenHisto(0),
  fSMRatHisto(0),
  fCorrectionHisto(0),
  fCorrCFHisto(0),
  fNumerator(0),
  fDenominator(0),
  fRatio(0),
  fQinvHisto(0),
  fLambda(0),
  fRout2(0),
  fRside2(0),
  fRlong2(0),
  fPairCut(0), 
  fQinvNormLo(0),
  fQinvNormHi(0),
  fNumRealsNorm(0),
  fNumMixedNorm(0)
{

  // set some stuff...
  fQinvNormLo = 0.15;
  fQinvNormHi = 0.18;
  fNumRealsNorm = 0;
  fNumMixedNorm = 0;
  //  fCorrection = 0;  // pointer to Coulomb Correction object

  fPairCut = 0; // added Sept2000 - CorrFctn-specific PairCut

  //  fSmearPair = 0; // no resolution correction unless user sets SmearPair

  // set up numerator
  char TitNum[100] = "Num";
  strcat(TitNum,title);
  fNumerator = new TH3D(TitNum,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up denominator
  char TitDen[100] = "Den";
  strcat(TitDen,title);
  fDenominator = new TH3D(TitDen,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up uncorrected denominator
  char TitDenUncoul[100] = "DenNoCoul";
  strcat(TitDenUncoul,title);
  //  fUncorrectedDenominator = new TH3D(TitDenUncoul,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up ratio
  char TitRat[100] = "Rat";
  strcat(TitRat,title);
  fRatio = new TH3D(TitRat,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up ave qInv
  char TitQinv[100] = "Qinv";
  strcat(TitQinv,title);
  fQinvHisto = new TH3D(TitQinv,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  //  fUncorrectedDenominator->Sumw2();
  fRatio->Sumw2();

  // Following histos are for the momentum resolution correction
  // they are filled only if a AliFemtoSmear object is plugged in
  // here comes the "idea" numerator and denominator and ratio...
  char TitNumID[100] = "IDNum";
  strcat(TitNumID,title);
  fIDNumHisto = new TH3D(TitNumID,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  char TitDenID[100] = "IDDen";
  strcat(TitDenID,title);
  fIDDenHisto = new TH3D(TitDenID,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  char TitRatID[100] = "IDRat";
  strcat(TitRatID,title);
  fIDRatHisto = new TH3D(TitRatID,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);

  fIDNumHisto->Sumw2();
  fIDDenHisto->Sumw2();
  fIDRatHisto->Sumw2();

  //
  // here comes the "smeared" numerator and denominator...
  char TitNumSM[100] = "SMNum";
  strcat(TitNumSM,title);
  fSMNumHisto = new TH3D(TitNumSM,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  char TitDenSM[100] = "SMDen";
  strcat(TitDenSM,title);
  fSMDenHisto = new TH3D(TitDenSM,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  char TitRatSM[100] = "SMRat";
  strcat(TitRatSM,title);
  fSMRatHisto = new TH3D(TitRatSM,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  //
  fSMNumHisto->Sumw2();
  fSMDenHisto->Sumw2();
  fSMRatHisto->Sumw2();
  //
  // here comes the correction factor (which is just ratio of ideal ratio to smeared ratio)
  char TitCorrection[100] = "CorrectionFactor";
  strcat(TitCorrection,title);
  fCorrectionHisto = new TH3D(TitCorrection,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);  
  fCorrectionHisto->Sumw2();
  // here comes the fully corrected correlation function
  char TitCorrCF[100] = "CorrectedCF";
  strcat(TitCorrCF,title);
  fCorrCFHisto = new TH3D(TitCorrCF,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  fCorrCFHisto->Sumw2();

  // user can (and should) override these defaults...
  fLambda = 0.6;
  fRout2 = 6.0*6.0;
  fRside2 = 6.0*6.0;
  fRlong2 = 7.0*7.0;

}

AliFemtoBPLCMS3DCorrFctn::AliFemtoBPLCMS3DCorrFctn(const AliFemtoBPLCMS3DCorrFctn& aCorrFctn) :
  fIDNumHisto(0),
  fIDDenHisto(0),
  fIDRatHisto(0),
  fSMNumHisto(0),
  fSMDenHisto(0),
  fSMRatHisto(0),
  fCorrectionHisto(0),
  fCorrCFHisto(0),
  fNumerator(0),
  fDenominator(0),
  fRatio(0),
  fQinvHisto(0),
  fLambda(0),
  fRout2(0),
  fRside2(0),
  fRlong2(0),
  fPairCut(0), 
  fQinvNormLo(0),
  fQinvNormHi(0),
  fNumRealsNorm(0),
  fNumMixedNorm(0)
{
  fIDNumHisto = aCorrFctn.fIDNumHisto;
  fIDDenHisto = aCorrFctn.fIDDenHisto;
  fIDRatHisto = aCorrFctn.fIDRatHisto;
  fSMNumHisto = aCorrFctn.fSMNumHisto;
  fSMDenHisto = aCorrFctn.fSMDenHisto;
  fSMRatHisto = aCorrFctn.fSMRatHisto;
  fCorrectionHisto = aCorrFctn.fCorrectionHisto;
  fCorrCFHisto = aCorrFctn.fCorrCFHisto;
  fNumerator = aCorrFctn.fNumerator;
  fDenominator = aCorrFctn.fDenominator;
  fRatio = aCorrFctn.fRatio;
  fQinvHisto = aCorrFctn.fQinvHisto;
  fLambda = aCorrFctn.fLambda;
  fRout2 = aCorrFctn.fRout2;
  fRside2 = aCorrFctn.fRside2;
  fRlong2 = aCorrFctn.fRlong2;
  fPairCut = aCorrFctn.fPairCut; 
  fQinvNormLo = aCorrFctn.fQinvNormLo;
  fQinvNormHi = aCorrFctn.fQinvNormHi;
  fNumRealsNorm = aCorrFctn.fNumRealsNorm;
  fNumMixedNorm = aCorrFctn.fNumMixedNorm;
}
//____________________________
AliFemtoBPLCMS3DCorrFctn::~AliFemtoBPLCMS3DCorrFctn(){
  delete fNumerator;
  delete fDenominator;
  delete fRatio;
  delete fQinvHisto;
  delete fIDNumHisto;
  delete fIDDenHisto;
  delete fIDRatHisto;
  delete fSMNumHisto;
  delete fSMDenHisto;
  delete fSMRatHisto;
  delete fCorrectionHisto;
  delete fCorrCFHisto;
}
//_________________________
AliFemtoBPLCMS3DCorrFctn& AliFemtoBPLCMS3DCorrFctn::operator=(const AliFemtoBPLCMS3DCorrFctn& aCorrFctn)
{
  if (this == &aCorrFctn)
    return *this;
  if (fIDNumHisto) delete fIDNumHisto;
  fIDNumHisto = new TH3D(*aCorrFctn.fIDNumHisto);
  if (fIDDenHisto) delete fIDDenHisto;
  fIDDenHisto = new TH3D(*aCorrFctn.fIDDenHisto);
  if (fIDRatHisto) delete fIDRatHisto;
  fIDRatHisto = new TH3D(*aCorrFctn.fIDRatHisto);
  if (fSMNumHisto) delete fSMNumHisto;
  fSMNumHisto = new TH3D(*aCorrFctn.fSMNumHisto);
  if (fSMDenHisto) delete fSMDenHisto;
  fSMDenHisto = new TH3D(*aCorrFctn.fSMDenHisto);
  if (fSMRatHisto) delete fSMRatHisto;
  fSMRatHisto = new TH3D(*aCorrFctn.fSMRatHisto);

  if (fCorrectionHisto) delete fCorrectionHisto;
  fCorrectionHisto = new TH3D(*aCorrFctn.fCorrectionHisto);
  if (fCorrCFHisto) delete fCorrCFHisto;
  fCorrCFHisto = new TH3D(*aCorrFctn.fCorrCFHisto);
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
  fPairCut = aCorrFctn.fPairCut; 
  fQinvNormLo = aCorrFctn.fQinvNormLo;
  fQinvNormHi = aCorrFctn.fQinvNormHi;
  fNumRealsNorm = aCorrFctn.fNumRealsNorm;
  fNumMixedNorm = aCorrFctn.fNumMixedNorm;
  
  return *this;
}

//_________________________
void AliFemtoBPLCMS3DCorrFctn::WriteOutHistos(){

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

//_________________________
void AliFemtoBPLCMS3DCorrFctn::Finish(){
  // here is where we should normalize, fit, etc...
  double NumFact,DenFact;
  if ((fNumRealsNorm !=0) && (fNumMixedNorm !=0)){
    NumFact = double(fNumRealsNorm);
    DenFact = double(fNumMixedNorm);
  }
  // can happen that the fNumRealsNorm and fNumMixedNorm = 0 if you do non-standard
  //   things like making a new CorrFctn and just setting the Numerator and Denominator
  //   from OTHER CorrFctns which you read in (like when doing parallel processing) 
  else{
    cout << "Warning! - no normalization constants defined - I do the best I can..." << endl;
    int nbins = fNumerator->GetNbinsX();
    int half_way = nbins/2;
    NumFact = fNumerator->Integral(half_way,nbins,half_way,nbins,half_way,nbins);
    DenFact = fDenominator->Integral(half_way,nbins,half_way,nbins,half_way,nbins);
  }

  fRatio->Divide(fNumerator,fDenominator,DenFact,NumFact);
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
void AliFemtoBPLCMS3DCorrFctn::AddRealPair(const AliFemtoPair* pair){

  if (fPairCut){
    if (!(fPairCut->Pass(pair))) return;
  }

  double Qinv = fabs(pair->qInv());   // note - qInv() will be negative for identical pairs...
  if ((Qinv < fQinvNormHi) && (Qinv > fQinvNormLo)) fNumRealsNorm++;
  double qOut = fabs(pair->qOutCMS());
  double qSide = fabs(pair->qSideCMS());
  double qLong = fabs(pair->qLongCMS());

  fNumerator->Fill(qOut,qSide,qLong);
}
//____________________________
void AliFemtoBPLCMS3DCorrFctn::AddMixedPair(const AliFemtoPair* pair){

  if (fPairCut){
    if (!(fPairCut->Pass(pair))) return;
  }

  //  double CoulombWeight = (fCorrection ? fCorrection->CoulombCorrect(pair) : 1.0);
  double CoulombWeight = 1.0;

  double Qinv = fabs(pair->qInv());   // note - qInv() will be negative for identical pairs...
  if ((Qinv < fQinvNormHi) && (Qinv > fQinvNormLo)) fNumMixedNorm++;
  double qOut = fabs(pair->qOutCMS());
  double qSide = fabs(pair->qSideCMS());
  double qLong = fabs(pair->qLongCMS());

  fDenominator->Fill(qOut,qSide,qLong,CoulombWeight);
  //  fUncorrectedDenominator->Fill(qOut,qSide,qLong,1.0);
  fQinvHisto->Fill(qOut,qSide,qLong,Qinv);

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


