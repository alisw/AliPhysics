////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnGammaMonitor - A correlation function that analyzes        //
// two particle mass minvariant with various mass assumptions                 //
//                                                                            //
// Authors: Małgorzata Janik majanik@cern.ch
//          Anna Zaborowska azaborow@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnMinvMonitor.h"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoCorrFctnMinvMonitor)
#endif

//____________________________
AliFemtoCorrFctnMinvMonitor::AliFemtoCorrFctnMinvMonitor(const char* title):
   AliFemtoCorrFctn(),
   fMinveeFail(0),
   fMinvee(0),
   fMinv2piFail(0),
   fMinv2pi(0),
   fMinvppiFail(0),
   fMinvppi(0)
{
  fMinveeFail = new TH1D(Form("MinveeGamma%s",title), "ee mass assumption GAMMA, minv",1000, 0.0, 10.0);
   fMinvee = new TH1D(Form("Minvee%s",title), "ee mass assumption, minv",1000, 0.0, 10.0);
   fMinv2piFail = new TH1D(Form("Minv2piResonances%s",title), "pipi mass assumption RESONANCES, minv",1000, 0.0, 10.0);
   fMinv2pi = new TH1D(Form("Minv2pi%s",title), "pipi mass assumption, minv",1000, 0.0, 10.0);
   fMinvppiFail = new TH1D(Form("MinvppiResonances%s",title), "ppi mass assumption RESONANCES, minv",1000, 0.0, 10.0);
   fMinvppi = new TH1D(Form("Minvppi%s",title), "ppi mass assumption, minv",1000, 0.0, 10.0);
}

//____________________________
AliFemtoCorrFctnMinvMonitor::AliFemtoCorrFctnMinvMonitor(const AliFemtoCorrFctnMinvMonitor& aCorrFctn) :
  AliFemtoCorrFctn(),
   fMinveeFail(0),
   fMinvee(0),
   fMinv2piFail(0),
   fMinv2pi(0),
   fMinvppiFail(0),
   fMinvppi(0)
{
  // copy constructor
  if (fMinveeFail) delete fMinveeFail;
  fMinveeFail = new TH1D(*aCorrFctn.fMinveeFail);
  if (fMinvee) delete fMinvee;
  fMinvee = new TH1D(*aCorrFctn.fMinvee);
  if (fMinv2piFail) delete fMinv2piFail;
  fMinv2piFail = new TH1D(*aCorrFctn.fMinv2piFail);
  if (fMinv2pi) delete fMinv2pi;
  fMinv2pi = new TH1D(*aCorrFctn.fMinv2pi);
  if (fMinvppiFail) delete fMinvppiFail;
  fMinvppiFail = new TH1D(*aCorrFctn.fMinvppiFail);
  if (fMinvppi) delete fMinvppi;
  fMinvppi = new TH1D(*aCorrFctn.fMinvppi);
}
//____________________________
AliFemtoCorrFctnMinvMonitor::~AliFemtoCorrFctnMinvMonitor(){
  // destructor
    delete fMinveeFail;
  delete fMinvee;
  delete fMinv2piFail;
  delete fMinv2pi;
  delete fMinvppiFail;
  delete fMinvppi;
}
//_________________________
AliFemtoCorrFctnMinvMonitor& AliFemtoCorrFctnMinvMonitor::operator=(const AliFemtoCorrFctnMinvMonitor& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (fMinveeFail) delete fMinveeFail;
  fMinveeFail = new TH1D(*aCorrFctn.fMinveeFail);
  if (fMinvee) delete fMinvee;
  fMinvee = new TH1D(*aCorrFctn.fMinvee);
  if (fMinv2piFail) delete fMinv2piFail;
  fMinv2piFail = new TH1D(*aCorrFctn.fMinv2piFail);
  if (fMinv2pi) delete fMinv2pi;
  fMinv2pi = new TH1D(*aCorrFctn.fMinv2pi);
  if (fMinvppiFail) delete fMinvppiFail;
  fMinvppiFail = new TH1D(*aCorrFctn.fMinvppiFail);
  if (fMinvppi) delete fMinvppi;
  fMinvppi = new TH1D(*aCorrFctn.fMinvppi);

  return *this;
}
//_________________________
void AliFemtoCorrFctnMinvMonitor::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnMinvMonitor::Report(){
  // create report
  string stemp = "Mass invariant Monitor Function Report\n";
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctnMinvMonitor::AddRealPair( AliFemtoPair* pair){
   double me = 0.000511;
  double mPi = 0.13957018;
  double mp = 0.938272046;

  double mgammamax = 0.04;
  double mK0min = 0.00049;
  double mK0max = 0.00051;
  //double mK0 = 0.000497614;
  double mRhomin = 0.000765;
  double mRhomax = 0.000785;
  //double mRho = 0.00077526;
  double mLmin = 1.095;
  double mLmax = 1.135;
  //double mL = 1.115683;

  if ((pair->Track1()->Track()->Charge() * pair->Track2()->Track()->Charge()) < 0.0) {

    // check on ee pairs (gamma)
    double e1 = TMath::Sqrt(me*me + pair->Track1()->Track()->P().Mag2());
    double e2 = TMath::Sqrt(me*me + pair->Track2()->Track()->P().Mag2());
    double minvGamma = 2*me*me + 2*(e1*e2 -
			       pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			       pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
			       pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z());
    if ( minvGamma < mgammamax )
    {
       fMinveeFail->Fill(minvGamma);
    }
    else fMinvee->Fill(minvGamma);
    //check on resonances
    double pi1 =  TMath::Sqrt(mPi*mPi + pair->Track1()->Track()->P().Mag2());
    double pi2 =  TMath::Sqrt(mPi*mPi + pair->Track2()->Track()->P().Mag2());
    double p1 =  TMath::Sqrt(mp*mp + pair->Track1()->Track()->P().Mag2());
    double p2 =  TMath::Sqrt(mp*mp + pair->Track2()->Track()->P().Mag2());
    //check on K0 and Rho
    double minv2pi = 2*mPi*mPi + 2*(pi1*pi2 -
			       pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			       pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
			       pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z());
    if ( ((minv2pi>mK0min && minv2pi<mK0max) || (minv2pi>mRhomin && minv2pi<mRhomax)) )
    {
       fMinv2piFail->Fill(minv2pi);
    }
    else fMinv2pi->Fill(minv2pi);
    //check on L0
    double minvpPi = 2*mp*mPi + 2*(p1*pi2 -
			       pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			       pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
			       pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z());
    double minvPip = 2*mPi*mp + 2*(pi1*p2 -
			       pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			       pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
			       pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z());
    if( (minvpPi>mLmin) && (minvpPi<mLmax) )
    {
       fMinvppiFail->Fill(minvpPi);
    }
    else fMinvppi->Fill(minvpPi);
    if( (minvPip>mLmin) && (minvPip<mLmax) )
    {
       fMinvppiFail->Fill(minvPip);
    }
    else fMinvppi->Fill(minvPip);
  }
}
//____________________________

void AliFemtoCorrFctnMinvMonitor::WriteHistos()
{
  // Write out result histograms
  fMinveeFail->Write();
  fMinvee->Write();
  fMinv2piFail->Write();
  fMinv2pi->Write();
  fMinvppiFail->Write();
  fMinvppi->Write();
}

TList* AliFemtoCorrFctnMinvMonitor::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fMinveeFail);
  tOutputList->Add(fMinvee);
  tOutputList->Add(fMinv2piFail);
  tOutputList->Add(fMinv2pi);
  tOutputList->Add(fMinvppiFail);
  tOutputList->Add(fMinvppi);

  return tOutputList;
}
