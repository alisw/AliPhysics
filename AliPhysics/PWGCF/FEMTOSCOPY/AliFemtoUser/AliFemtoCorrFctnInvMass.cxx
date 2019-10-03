////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhi - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnInvMass.h"
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoCorrFctnInvMass)
#endif

//____________________________
AliFemtoCorrFctnInvMass::AliFemtoCorrFctnInvMass(const char* title, const int& aBins=3000, const double& aMin = 1000, const double& aMax = 5000, const double& aMass1 = 0, const double& aMass2 = 0):
  AliFemtoCorrFctn(),
  fNumInvMass(nullptr),
  fDenInvMass(nullptr),
  fBins(aBins),
  fMin(aMin),
  fMax(aMax),
  fMass1(aMass1),
  fMass2(aMass2)
{
  // set up numerator
  char *num_name = Form("NumInvMass%s", title);
  char *num_title = Form("Minv Numerator - (masses = %g, %g); M_{inv} GeV", fMass1, fMass2);
  fNumInvMass = new TH1D(num_name, num_title, aBins,aMin,aMax);

  // set up denominator
  char *den_name = Form("DenInvMass%s", title);
  char *den_title = Form("Minv Denominator - (masses = %g, %g); M_{inv} GeV", fMass1, fMass2);
  fDenInvMass = new TH1D(den_name, den_title, aBins,aMin,aMax);

  // to enable error bar calculation...
  fNumInvMass->Sumw2();
  fDenInvMass->Sumw2();
}

//____________________________
AliFemtoCorrFctnInvMass::AliFemtoCorrFctnInvMass(const AliFemtoCorrFctnInvMass& aCorrFctn) :
  AliFemtoCorrFctn(aCorrFctn),
  fNumInvMass(nullptr),
  fDenInvMass(nullptr),
  fBins(aCorrFctn.fBins),
  fMin(aCorrFctn.fMin),
  fMax(aCorrFctn.fMax),
  fMass1(aCorrFctn.fMass1),
  fMass2(aCorrFctn.fMass2)
{
  // copy constructor
  fNumInvMass = new TH1D(*aCorrFctn.fNumInvMass);
  fDenInvMass = new TH1D(*aCorrFctn.fDenInvMass);
}

//____________________________
AliFemtoCorrFctnInvMass::~AliFemtoCorrFctnInvMass()
{
  // destructor

  delete fNumInvMass;
  delete fDenInvMass;
}

//_________________________
AliFemtoCorrFctnInvMass& AliFemtoCorrFctnInvMass::operator=(const AliFemtoCorrFctnInvMass& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  fBins = aCorrFctn.fBins;
  fMin = aCorrFctn.fMax;
  fMax = aCorrFctn.fMax;

  fMass1 = aCorrFctn.fMass1;
  fMass2 = aCorrFctn.fMass2;

  *fNumInvMass = *aCorrFctn.fNumInvMass;
  *fDenInvMass = *aCorrFctn.fDenInvMass;

  return *this;
}
//_________________________
void AliFemtoCorrFctnInvMass::Finish()
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
AliFemtoString AliFemtoCorrFctnInvMass::Report()
{
  // create report
  AliFemtoString report = "Invariant Mass Correlation Function Report:\n";
  report += Form("Number of entries in numerator:\t%E\n",fNumInvMass->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n",fDenInvMass->GetEntries());
  report += Form("Assumed masses:\t%g\t%g\n", fMass1, fMass2);

  return report;
}
//____________________________
void AliFemtoCorrFctnInvMass::AddRealPair(AliFemtoPair* pair)
{
  // add real (effect) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  double px1 = 999999;
  double py1 = 999999;
  double pz1 = 999999;
  if(pair->Track1()->Track())
  {
    px1 = pair->Track1()->Track()->P().x();
    py1 = pair->Track1()->Track()->P().y();
    pz1 = pair->Track1()->Track()->P().z();
  }
  else if(pair->Track1()->V0())
  {
    px1 = pair->Track1()->V0()->MomV0().x();
    py1 = pair->Track1()->V0()->MomV0().y();
    pz1 = pair->Track1()->V0()->MomV0().z();
  }
  double mass1 = fMass1;

  double E1 = TMath::Sqrt(px1*px1+py1*py1+pz1*pz1+mass1*mass1);

  double px2 = 999999;
  double py2 = 999999;
  double pz2 = 999999;
  if(pair->Track2()->Track())
  {
    px2 = pair->Track2()->Track()->P().x();
    py2 = pair->Track2()->Track()->P().y();
    pz2 = pair->Track2()->Track()->P().z();
  }
  else if(pair->Track2()->V0())
  {
    px2 = pair->Track2()->V0()->MomV0().x();
    py2 = pair->Track2()->V0()->MomV0().y();
    pz2 = pair->Track2()->V0()->MomV0().z();
  }
  double mass2 = fMass2;

  double E2 = TMath::Sqrt(px2*px2+py2*py2+pz2*pz2+mass2*mass2);

  double minv = TMath::Sqrt((E1+E2)*(E1+E2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));
  fNumInvMass->Fill(minv);
}
//____________________________
void AliFemtoCorrFctnInvMass::AddMixedPair(AliFemtoPair* pair)
{
  // add mixed (background) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  double px1 = 999999;
  double py1 = 999999;
  double pz1 = 999999;
  if(pair->Track1()->Track())
  {
    px1 = pair->Track1()->Track()->P().x();
    py1 = pair->Track1()->Track()->P().y();
    pz1 = pair->Track1()->Track()->P().z();
  }
  else if(pair->Track1()->V0())
  {
    px1 = pair->Track1()->V0()->MomV0().x();
    py1 = pair->Track1()->V0()->MomV0().y();
    pz1 = pair->Track1()->V0()->MomV0().z();
  }
  double mass1 = fMass1;

  double E1 = TMath::Sqrt(px1*px1+py1*py1+pz1*pz1+mass1*mass1);


  double px2 = 999999;
  double py2 = 999999;
  double pz2 = 999999;
  if(pair->Track2()->Track())
  {
    px2 = pair->Track2()->Track()->P().x();
    py2 = pair->Track2()->Track()->P().y();
    pz2 = pair->Track2()->Track()->P().z();
  }
  else if(pair->Track2()->V0())
  {
    px2 = pair->Track2()->V0()->MomV0().x();
    py2 = pair->Track2()->V0()->MomV0().y();
    pz2 = pair->Track2()->V0()->MomV0().z();
  }
  double mass2 = fMass2;

  double E2 = TMath::Sqrt(px2*px2+py2*py2+pz2*pz2+mass2*mass2);

  double minv = TMath::Sqrt((E1+E2)*(E1+E2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));
  fDenInvMass->Fill(minv);
}


void AliFemtoCorrFctnInvMass::WriteHistos()
{
  // Write out result histograms
  fNumInvMass->Write();
  fDenInvMass->Write();
}

TList* AliFemtoCorrFctnInvMass::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumInvMass);
  tOutputList->Add(fDenInvMass);

  return tOutputList;
}

void AliFemtoCorrFctnInvMass::SetParticleMasses(double mass1, double mass2)
{
  fMass1=mass1;
  fMass2=mass2;
  fNumInvMass->SetTitle(Form("Minv Numerator - masses = %g, %g", fMass1, fMass2));
  fDenInvMass->SetTitle(Form("Minv Denominator - masses = %g, %g", fMass1, fMass2));
}

void AliFemtoCorrFctnInvMass::SetParticle1Mass(double mass)
{
  SetParticleMasses(mass, fMass2);
}

void AliFemtoCorrFctnInvMass::SetParticle2Mass(double mass)
{
  SetParticleMasses(fMass1, mass);
}
