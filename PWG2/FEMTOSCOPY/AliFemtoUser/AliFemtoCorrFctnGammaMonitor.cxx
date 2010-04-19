////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoGammaMonitor - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnGammaMonitor.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__ 
ClassImp(AliFemtoCorrFctnGammaMonitor)
#endif

//____________________________
AliFemtoCorrFctnGammaMonitor::AliFemtoCorrFctnGammaMonitor(char* title, const int& aMinvBins=20, const int& aDThetaBins=20):
  AliFemtoCorrFctn(),
  fNumPMinvDTheta(0),
  fDenPMinvDTheta(0),
  fNumNMinvDTheta(0),
  fDenNMinvDTheta(0)
{
  // set up numerator
  char tTitNumD[100] = "NumPMinvTheta";
  strcat(tTitNumD,title);
  fNumPMinvDTheta = new TH2D(tTitNumD,title,aMinvBins,0.0,0.2,aDThetaBins,0.0,0.2);
  // set up denominator
  char tTitDenD[100] = "DenPMinvTheta";
  strcat(tTitDenD,title);
  fDenPMinvDTheta = new TH2D(tTitDenD,title,aMinvBins,0.0,0.2,aDThetaBins,0.0,0.2);

  // set up numerator
  char tTitNumR[100] = "NumNMinvTheta";
  strcat(tTitNumR,title);
  fNumNMinvDTheta = new TH2D(tTitNumR,title,aMinvBins,0.0,0.2,aDThetaBins,0.0,0.2);
  // set up denominator
  char tTitDenR[100] = "DenNMinvTheta";
  strcat(tTitDenR,title);
  fDenNMinvDTheta = new TH2D(tTitDenR,title,aMinvBins,0.0,0.2,aDThetaBins,0.0,0.2);

  // to enable error bar calculation...
  fNumPMinvDTheta->Sumw2();
  fDenPMinvDTheta->Sumw2();
  fNumNMinvDTheta->Sumw2();
  fDenNMinvDTheta->Sumw2();
}

//____________________________
AliFemtoCorrFctnGammaMonitor::AliFemtoCorrFctnGammaMonitor(const AliFemtoCorrFctnGammaMonitor& aCorrFctn) :
  AliFemtoCorrFctn(),
  fNumPMinvDTheta(0),
  fDenPMinvDTheta(0),
  fNumNMinvDTheta(0),
  fDenNMinvDTheta(0)
{
  // copy constructor
  if (aCorrFctn.fNumPMinvDTheta)
    fNumPMinvDTheta = new TH2D(*aCorrFctn.fNumPMinvDTheta);
  else
    fNumPMinvDTheta = 0;
  if (aCorrFctn.fDenPMinvDTheta)
    fDenPMinvDTheta = new TH2D(*aCorrFctn.fDenPMinvDTheta);
  else
    fDenPMinvDTheta = 0;

  if (aCorrFctn.fNumNMinvDTheta)
    fNumNMinvDTheta = new TH2D(*aCorrFctn.fNumNMinvDTheta);
  else
    fNumNMinvDTheta = 0;
  if (aCorrFctn.fDenNMinvDTheta)
    fDenNMinvDTheta = new TH2D(*aCorrFctn.fDenNMinvDTheta);
  else
    fDenNMinvDTheta = 0;

}
//____________________________
AliFemtoCorrFctnGammaMonitor::~AliFemtoCorrFctnGammaMonitor(){
  // destructor
  delete fNumPMinvDTheta;
  delete fDenPMinvDTheta;
  delete fNumNMinvDTheta;
  delete fDenNMinvDTheta;
}
//_________________________
AliFemtoCorrFctnGammaMonitor& AliFemtoCorrFctnGammaMonitor::operator=(const AliFemtoCorrFctnGammaMonitor& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fNumPMinvDTheta)
    fNumPMinvDTheta = new TH2D(*aCorrFctn.fNumPMinvDTheta);
  else
    fNumPMinvDTheta = 0;
  if (aCorrFctn.fDenPMinvDTheta)
    fDenPMinvDTheta = new TH2D(*aCorrFctn.fDenPMinvDTheta);
  else
    fDenPMinvDTheta = 0;

  if (aCorrFctn.fNumNMinvDTheta)
    fNumNMinvDTheta = new TH2D(*aCorrFctn.fNumNMinvDTheta);
  else
    fNumNMinvDTheta = 0;
  if (aCorrFctn.fDenNMinvDTheta)
    fDenNMinvDTheta = new TH2D(*aCorrFctn.fDenNMinvDTheta);
  else
    fDenNMinvDTheta = 0;


  return *this;
}
//_________________________
void AliFemtoCorrFctnGammaMonitor::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnGammaMonitor::Report(){
  // create report
  string stemp = "Gamma Monitor Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerator:\t%E\n",fNumPMinvDTheta->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominator:\t%E\n",fDenPMinvDTheta->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctnGammaMonitor::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  double me = 0.000511;

  double theta1 = pair->Track1()->Track()->P().Theta();
  double theta2 = pair->Track2()->Track()->P().Theta();
  double dtheta = TMath::Abs(theta1 - theta2);

  double e1 = TMath::Sqrt(me*me + pair->Track1()->Track()->P().Mag2());
  double e2 = TMath::Sqrt(me*me + pair->Track2()->Track()->P().Mag2());

  double minv = 2*me*me + 2*(e1*e2 - 
			     pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			     pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
			     pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z());

  if (pair->KSide()>0.0) 
    fNumPMinvDTheta->Fill(minv, dtheta);
  else
    fNumNMinvDTheta->Fill(minv, dtheta);
}
//____________________________
void AliFemtoCorrFctnGammaMonitor::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  double me = 0.000511;

  double theta1 = pair->Track1()->Track()->P().Theta();
  double theta2 = pair->Track2()->Track()->P().Theta();
  double dtheta = TMath::Abs(theta1 - theta2);

  double e1 = TMath::Sqrt(me*me + pair->Track1()->Track()->P().Mag2());
  double e2 = TMath::Sqrt(me*me + pair->Track2()->Track()->P().Mag2());

  double minv = 2*me*me + 2*(e1*e2 - 
			     pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			     pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
			     pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z());

  if (pair->KSide()>0.0) 
    fDenPMinvDTheta->Fill(minv, dtheta);
  else
    fDenNMinvDTheta->Fill(minv, dtheta);
}


void AliFemtoCorrFctnGammaMonitor::WriteHistos()
{
  // Write out result histograms
  fNumPMinvDTheta->Write();
  fDenPMinvDTheta->Write();
  fNumNMinvDTheta->Write();
  fDenNMinvDTheta->Write();
}

TList* AliFemtoCorrFctnGammaMonitor::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumPMinvDTheta);
  tOutputList->Add(fDenPMinvDTheta);
  tOutputList->Add(fNumNMinvDTheta);
  tOutputList->Add(fDenNMinvDTheta);

  return tOutputList;
}
