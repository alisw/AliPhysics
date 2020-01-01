////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoGammaMonitorAlpha - A correlation function that analyzes           //
// two particle correlations with respect to the angle between two particles  //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Elena Rogochaya Adam.Kisiel@cern.ch                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnGammaMonitorAlpha.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoCorrFctnGammaMonitorAlpha)
#endif

//____________________________
AliFemtoCorrFctnGammaMonitorAlpha::AliFemtoCorrFctnGammaMonitorAlpha(const char* title, const int& aMinvBins=20, const int& aDAlphaBins=20):
  AliFemtoCorrFctn(),
  fNumPMinvDAlpha(0),
  fDenPMinvDAlpha(0),
  fNumNMinvDAlpha(0),
  fDenNMinvDAlpha(0)
{
  // set up numerator
  char tTitNumD[101] = "NumPMinvAlpha";
  strncat(tTitNumD,title, 100);
  fNumPMinvDAlpha = new TH2D(tTitNumD,title,aMinvBins,0.0,1.0,aDAlphaBins,0.996,1.0);
  // set up denominator
  char tTitDenD[101] = "DenPMinvAlpha";
  strncat(tTitDenD,title, 100);
  fDenPMinvDAlpha = new TH2D(tTitDenD,title,aMinvBins,0.0,1.0,aDAlphaBins,0.996,1.0);

  // set up numerator
  char tTitNumR[101] = "NumNMinvAlpha";
  strncat(tTitNumR,title, 100);
  fNumNMinvDAlpha = new TH2D(tTitNumR,title,aMinvBins,0.0,1.0,aDAlphaBins,0.996,1.0);
  // set up denominator
  char tTitDenR[101] = "DenNMinvAlpha";
  strncat(tTitDenR,title, 100);
  fDenNMinvDAlpha = new TH2D(tTitDenR,title,aMinvBins,0.0,1.0,aDAlphaBins,0.996,1.0);

  // to enable error bar calculation...
  fNumPMinvDAlpha->Sumw2();
  fDenPMinvDAlpha->Sumw2();
  fNumNMinvDAlpha->Sumw2();
  fDenNMinvDAlpha->Sumw2();
}

//____________________________
AliFemtoCorrFctnGammaMonitorAlpha::AliFemtoCorrFctnGammaMonitorAlpha(const AliFemtoCorrFctnGammaMonitorAlpha& aCorrFctn) :
  AliFemtoCorrFctn(),
  fNumPMinvDAlpha(nullptr),
  fDenPMinvDAlpha(nullptr),
  fNumNMinvDAlpha(nullptr),
  fDenNMinvDAlpha(nullptr)
{
  // copy constructor
  fNumPMinvDAlpha = new TH2D(*aCorrFctn.fNumPMinvDAlpha);
  fDenPMinvDAlpha = new TH2D(*aCorrFctn.fDenPMinvDAlpha);
  fNumNMinvDAlpha = new TH2D(*aCorrFctn.fNumNMinvDAlpha);
  fDenNMinvDAlpha = new TH2D(*aCorrFctn.fDenNMinvDAlpha);
}
//____________________________
AliFemtoCorrFctnGammaMonitorAlpha::~AliFemtoCorrFctnGammaMonitorAlpha(){
  // destructor
  delete fNumPMinvDAlpha;
  delete fDenPMinvDAlpha;
  delete fNumNMinvDAlpha;
  delete fDenNMinvDAlpha;
}
//_________________________
AliFemtoCorrFctnGammaMonitorAlpha& AliFemtoCorrFctnGammaMonitorAlpha::operator=(const AliFemtoCorrFctnGammaMonitorAlpha& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  *fNumPMinvDAlpha = *aCorrFctn.fNumPMinvDAlpha;
  *fDenPMinvDAlpha = *aCorrFctn.fDenPMinvDAlpha;
  *fNumNMinvDAlpha = *aCorrFctn.fNumNMinvDAlpha;
  *fDenNMinvDAlpha = *aCorrFctn.fDenNMinvDAlpha;

  return *this;
}
//_________________________
void AliFemtoCorrFctnGammaMonitorAlpha::Finish()
{
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();
}

//____________________________
AliFemtoString AliFemtoCorrFctnGammaMonitorAlpha::Report()
{
  // create report
  AliFemtoString report = "Gamma MonitorAlpha Function Report:\n";
  report += Form("Number of entries in numerator:\t%E\n",fNumPMinvDAlpha->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n",fDenPMinvDAlpha->GetEntries());
  //  report += mCoulombWeight->Report();

  return report;
}
//____________________________
void AliFemtoCorrFctnGammaMonitorAlpha::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  double me = 0.000511;

  double dalpha = TMath::Abs((  pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x()
		   + pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y()
		   + pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z()
			       )/pair->Track1()->Track()->P().Mag()/pair->Track2()->Track()->P().Mag());

  double e1 = TMath::Sqrt(me*me + pair->Track1()->Track()->P().Mag2());
  double e2 = TMath::Sqrt(me*me + pair->Track2()->Track()->P().Mag2());

  double minv = (2*me*me + 2*(e1*e2 -
			     pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			     pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
					 pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z()));

  double sminv = TMath::Sqrt(minv);

  if (pair->KSide()>0.0) {
    fNumPMinvDAlpha->Fill(sminv, dalpha);
  }
  else {
    fNumNMinvDAlpha->Fill(sminv, dalpha);
  }
}
//____________________________
void AliFemtoCorrFctnGammaMonitorAlpha::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  double me = 0.000511;

  double dalpha = TMath::Abs((  pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x()
		   + pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y()
		   + pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z()
				)/pair->Track1()->Track()->P().Mag()/pair->Track2()->Track()->P().Mag());

  double e1 = TMath::Sqrt(me*me + pair->Track1()->Track()->P().Mag2());
  double e2 = TMath::Sqrt(me*me + pair->Track2()->Track()->P().Mag2());

  double minv = (2*me*me + 2*(e1*e2 -
			     pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			     pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
					 pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z()));

  double sminv = TMath::Sqrt(minv);

  if (pair->KSide()>0.0) {
    fDenPMinvDAlpha->Fill(sminv, dalpha);
  }
  else {
    fDenNMinvDAlpha->Fill(sminv, dalpha);
  }
}


void AliFemtoCorrFctnGammaMonitorAlpha::WriteHistos()
{
  // Write out result histograms
  fNumPMinvDAlpha->Write();
  fDenPMinvDAlpha->Write();
  fNumNMinvDAlpha->Write();
  fDenNMinvDAlpha->Write();
}

TList* AliFemtoCorrFctnGammaMonitorAlpha::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumPMinvDAlpha);
  tOutputList->Add(fDenPMinvDAlpha);
  tOutputList->Add(fNumNMinvDAlpha);
  tOutputList->Add(fDenNMinvDAlpha);

  return tOutputList;
}
