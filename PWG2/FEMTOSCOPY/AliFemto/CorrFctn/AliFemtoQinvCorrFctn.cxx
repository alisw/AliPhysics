/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   a simple Q-invariant correlation function           
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.4  2000/01/25 17:34:45  laue
 * I. In order to run the stand alone version of the AliFemtoMaker the following
 * changes have been done:
 * a) all ClassDefs and ClassImps have been put into #ifdef __ROOT__ statements
 * b) unnecessary includes of StMaker.h have been removed
 * c) the subdirectory AliFemtoMaker/doc/Make has been created including everything
 * needed for the stand alone version
 *
 * II. To reduce the amount of compiler warning
 * a) some variables have been type casted
 * b) some destructors have been declared as virtual
 *
 * Revision 1.3  1999/07/29 02:47:09  lisa
 * 1) add OpeningAngle correlation function 2) add AliFemtoMcEventReader 3) make histos in CorrFctns do errors correctly
 *
 * Revision 1.2  1999/07/06 22:33:20  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#include "CorrFctn/AliFemtoQinvCorrFctn.h"
//#include "Infrastructure/AliFemtoHisto.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoQinvCorrFctn)
#endif

//____________________________
AliFemtoQinvCorrFctn::AliFemtoQinvCorrFctn(char* title, const int& nbins, const float& QinvLo, const float& QinvHi){
  // set up numerator
  //  title = "Num Qinv (MeV/c)";
  char TitNum[100] = "Num";
  strcat(TitNum,title);
  fNumerator = new TH1D(TitNum,title,nbins,QinvLo,QinvHi);
  // set up denominator
  //title = "Den Qinv (MeV/c)";
  char TitDen[100] = "Den";
  strcat(TitDen,title);
  fDenominator = new TH1D(TitDen,title,nbins,QinvLo,QinvHi);
  // set up ratio
  //title = "Ratio Qinv (MeV/c)";
  char TitRat[100] = "Rat";
  strcat(TitRat,title);
  fRatio = new TH1D(TitRat,title,nbins,QinvLo,QinvHi);
  // this next bit is unfortunately needed so that we can have many histos of same "title"
  // it is neccessary if we typedef TH1D to TH1d (which we do)
  //fNumerator->SetDirectory(0);
  //fDenominator->SetDirectory(0);
  //fRatio->SetDirectory(0);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  fRatio->Sumw2();

}

//____________________________
AliFemtoQinvCorrFctn::~AliFemtoQinvCorrFctn(){
  delete fNumerator;
  delete fDenominator;
  delete fRatio;
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
  string stemp = "Qinv Correlation Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerator:\t%E\n",fNumerator->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominator:\t%E\n",fDenominator->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in ratio:\t%E\n",fRatio->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoQinvCorrFctn::AddRealPair(const AliFemtoPair* pair){
  double Qinv = fabs(pair->qInv());   // note - qInv() will be negative for identical pairs...
  fNumerator->Fill(Qinv);
  //  cout << "AliFemtoQinvCorrFctn::AddRealPair : " << pair->qInv() << " " << Qinv <<
  //" " << pair->track1().FourMomentum() << " " << pair->track2().FourMomentum() << endl;
}
//____________________________
void AliFemtoQinvCorrFctn::AddMixedPair(const AliFemtoPair* pair){
  double weight = 1.0;
  double Qinv = fabs(pair->qInv());   // note - qInv() will be negative for identical pairs...
  fDenominator->Fill(Qinv,weight);
}


