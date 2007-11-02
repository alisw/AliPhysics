////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelBPLCMSCorrFctn - the class for correlation function which   ///
/// uses the model framework and weight generation and calculated the 3D     ///
/// correlation function in the Bertsh-Pratt LCMS system                     ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoModelBPLCMSCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoModelBPLCMSCorrFctn)
#endif

//____________________________
AliFemtoModelBPLCMSCorrFctn::AliFemtoModelBPLCMSCorrFctn(char* title, const int& nbins, const float& QLo, const float& QHi)
  :
  AliFemtoModelCorrFctn(title, nbins, QLo, QHi),
  fNumerator3DTrue(0),
  fNumerator3DFake(0),
  fDenominator3D(0),
  fQinvHisto(0)
{

  // set up true numerator
  char TitNumT[100] = "Num3DTrue";
  strcat(TitNumT,title);
  fNumerator3DTrue = new TH3D(TitNumT,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up fake numerator
  char TitNumF[100] = "Num3DFake";
  strcat(TitNumF,title);
  fNumerator3DFake = new TH3D(TitNumF,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up denominator
  char TitDen[100] = "Den3D";
  strcat(TitDen,title);
  fDenominator3D = new TH3D(TitDen,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up ave qInv
  char TitQinv[100] = "Qinv";
  strcat(TitQinv,title);
  fQinvHisto = new TH3D(TitQinv,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);

  // to enable error bar calculation...
  fNumerator3DTrue->Sumw2();
  fNumerator3DFake->Sumw2();
  fDenominator3D->Sumw2();
}

AliFemtoModelBPLCMSCorrFctn::AliFemtoModelBPLCMSCorrFctn(const AliFemtoModelBPLCMSCorrFctn& aCorrFctn) :
  AliFemtoModelCorrFctn(aCorrFctn),
  fNumerator3DTrue(0),
  fNumerator3DFake(0),
  fDenominator3D(0),
  fQinvHisto(0)
{
  fNumerator3DTrue = new TH3D(*aCorrFctn.fNumerator3DTrue);
  fNumerator3DFake = new TH3D(*aCorrFctn.fNumerator3DFake);
  fDenominator3D   = new TH3D(*aCorrFctn.fDenominator3D);
  fQinvHisto       = new TH3D(*aCorrFctn.fQinvHisto);
}
//____________________________
AliFemtoModelBPLCMSCorrFctn::~AliFemtoModelBPLCMSCorrFctn()
{
  if (fNumeratorTrue) delete fNumeratorTrue;
  if (fNumeratorFake) delete fNumeratorFake;
  if (fDenominator) delete fDenominator;
  delete fNumerator3DTrue;
  delete fNumerator3DFake;
  delete fDenominator3D;
  delete fQinvHisto;
}
//_________________________
AliFemtoModelBPLCMSCorrFctn& AliFemtoModelBPLCMSCorrFctn::operator=(const AliFemtoModelBPLCMSCorrFctn& aCorrFctn)
{
  if (this == &aCorrFctn)
    return *this;
  if (fNumerator3DTrue) delete fNumerator3DTrue;
  fNumerator3DTrue = new TH3D(*aCorrFctn.fNumerator3DTrue);
  if (fNumerator3DFake) delete fNumerator3DFake;
  fNumerator3DFake = new TH3D(*aCorrFctn.fNumerator3DFake);
  if (fDenominator3D) delete fDenominator3D;
  fDenominator3D = new TH3D(*aCorrFctn.fDenominator3D);
  if (fQinvHisto) delete fQinvHisto;
  fQinvHisto = new TH3D(*aCorrFctn.fQinvHisto);

  return *this;
}

//_________________________
void AliFemtoModelBPLCMSCorrFctn::Write(){
  AliFemtoModelCorrFctn::Write();
  fNumerator3DTrue->Write();
  fNumerator3DFake->Write();
  fDenominator3D->Write();
  fQinvHisto->Write();
}

//_________________________
void AliFemtoModelBPLCMSCorrFctn::Finish(){
  fQinvHisto->Divide(fDenominator);
}

//____________________________
AliFemtoString AliFemtoModelBPLCMSCorrFctn::Report(){
  string stemp = "LCMS Frame Bertsch-Pratt 3D Model Correlation Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerator:\t%E\n",fNumeratorTrue->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominator:\t%E\n",fDenominator->GetEntries());
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

  //  
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoModelBPLCMSCorrFctn::AddRealPair( AliFemtoPair* pair){
  Double_t weight = fManager->GetWeight(pair);

  double qOut = fabs(pair->QOutCMS());
  double qSide = fabs(pair->QSideCMS());
  double qLong = fabs(pair->QLongCMS());

  fNumerator3DTrue->Fill(qOut, qSide, qLong, weight);
  fNumeratorTrue->Fill(pair->QInv(), weight);
}
//____________________________
void AliFemtoModelBPLCMSCorrFctn::AddMixedPair( AliFemtoPair* pair){
  Double_t weight = fManager->GetWeight(pair);

  double qOut = fabs(pair->QOutCMS());
  double qSide = fabs(pair->QSideCMS());
  double qLong = fabs(pair->QLongCMS());

  fNumerator3DFake->Fill(qOut, qSide, qLong, weight);
  fDenominator3D->Fill(qOut, qSide, qLong, 1.0);
  fNumeratorFake->Fill(pair->QInv(), weight);
  fDenominator->Fill(pair->QInv(), 1.0);

}
//_______________________
AliFemtoModelCorrFctn* AliFemtoModelBPLCMSCorrFctn::Clone()
{
  // Clone the correlation function
  AliFemtoModelBPLCMSCorrFctn *tCopy = new AliFemtoModelBPLCMSCorrFctn(*this);
  
  return tCopy;
}
