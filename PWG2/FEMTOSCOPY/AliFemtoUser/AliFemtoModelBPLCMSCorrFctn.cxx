////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelBPLCMSCorrFctn - the class for correlation function which   ///
/// uses the model framework and weight generation and calculated the 3D     ///
/// correlation function in the Bertsh-Pratt LCMS system                     ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoModelBPLCMSCorrFctn.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoModelBPLCMSCorrFctn)
#endif

//____________________________
AliFemtoModelBPLCMSCorrFctn::AliFemtoModelBPLCMSCorrFctn(char* title, const int& nbins, const float& QLo, const float& QHi)
  :
  fNumerator(0),
  fDenominator(0),
  fQinvHisto(0)
{

  // set up numerator
  char TitNum[100] = "Num";
  strcat(TitNum,title);
  fNumerator = new TH3D(TitNum,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up denominator
  char TitDen[100] = "Den";
  strcat(TitDen,title);
  fDenominator = new TH3D(TitDen,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up ave qInv
  char TitQinv[100] = "Qinv";
  strcat(TitQinv,title);
  fQinvHisto = new TH3D(TitQinv,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();
}

AliFemtoModelBPLCMSCorrFctn::AliFemtoModelBPLCMSCorrFctn(const AliFemtoModelBPLCMSCorrFctn& aCorrFctn) :
  fNumerator(0),
  fDenominator(0),
  fQinvHisto(0)
{
  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  fDenominator = new TH3D (*aCorrFctn.fDenominator);
  fQinvHisto = new TH3D(*aCorrFctn.fQinvHisto);
}
//____________________________
AliFemtoModelBPLCMSCorrFctn::~AliFemtoModelBPLCMSCorrFctn(){
  delete fNumerator;
  delete fDenominator;
  delete fQinvHisto;
}
//_________________________
AliFemtoModelBPLCMSCorrFctn& AliFemtoModelBPLCMSCorrFctn::operator=(const AliFemtoModelBPLCMSCorrFctn& aCorrFctn)
{
  if (this == &aCorrFctn)
    return *this;
  if (fNumerator) delete fNumerator;
  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  if (fDenominator) delete fDenominator;
  fDenominator = new TH3D(*aCorrFctn.fDenominator);
  if (fQinvHisto) delete fQinvHisto;
  fQinvHisto = new TH3D(*aCorrFctn.fQinvHisto);

  return *this;
}

//_________________________
void AliFemtoModelBPLCMSCorrFctn::WriteOutHistos(){

  fNumerator->Write();
  fDenominator->Write();
  fQinvHisto->Write();
}

//_________________________
void AliFemtoModelBPLCMSCorrFctn::Finish(){
  fQinvHisto->Divide(fDenominator);
}

//____________________________
AliFemtoString AliFemtoModelBPLCMSCorrFctn::Report(){
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
void AliFemtoModelBPLCMSCorrFctn::AddRealPair( AliFemtoPair* pair){

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
void AliFemtoModelBPLCMSCorrFctn::AddMixedPair( AliFemtoPair* pair){

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


