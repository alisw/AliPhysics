////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhi - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
// Edits: Daniel Rodak                                                                           //
//																			  //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnMezonPhi.h"
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoCorrFctnMezonPhi)
#endif

//____________________________
AliFemtoCorrFctnMezonPhi::AliFemtoCorrFctnMezonPhi(const char* title, const int& aBins=3000, const double& aMin = 1000, const double& aMax = 5000, const double& aMass1 = 0, const double& aMass2 = 0, const int& aSections = 1, const double& pMin = 0, const double& pMax = 0):
  AliFemtoCorrFctn(),
  fNumInvMass(nullptr),
  fDenInvMass(nullptr),
  fBins(aBins),
  fMin(aMin),
  fMax(aMax),
  fMass1(aMass1),
  fMass2(aMass2),
  fSections(aSections),
  fpMin(pMin),
  fpMax(pMax)
 {
  fNumInvMass = new TH1D*[fSections];
  fDenInvMass = new TH1D*[fSections]; 
  fNumTransvMom = new TH1D*[fSections]; 
  fDenTransvMom = new TH1D*[fSections];
  fIntersection = (fpMax - fpMin)/fSections;
  
  char *num_name_2D = Form("NumInvMass2D%s", title);
  char *num_title_2D = Form("Minv Numerator 2D - (masses = %g, %g); P_{t} GeV; M_{inv} GeV", fMass1, fMass2);
  
  fNumInvMass2D = new TH2D(num_name_2D, num_title_2D, aBins, fpMin, fpMax, aBins, fMin, fMax) ;
  
  char *den_name_2D = Form("DenInvMass2D%s", title);
  char *den_title_2D = Form("Minv Denominator 2D - (masses = %g, %g); P_{t} GeV; M_{inv}_GeV", fMass1, fMass2);
  
  fDenInvMass2D = new TH2D(den_name_2D, den_title_2D, aBins, fpMin, fpMax, aBins, fMin, fMax) ;

	for(int i = 0; i<fSections;i++){
		// set up numerator
		char *num_name = Form("NumInvMass%s%d", title, i);
		char *num_title = Form("Minv Numerator %d - (masses = %g, %g); M_{inv} GeV", i, fMass1, fMass2);
		fNumInvMass[i] = new TH1D(num_name, num_title, aBins,aMin,aMax);
		
		// set up denominator
		char *den_name = Form("DenInvMass%s%d", title, i);
		char *den_title = Form("Minv Denominator %d - (masses = %g, %g); M_{inv} GeV", i, fMass1, fMass2);
		fDenInvMass[i] = new TH1D(den_name, den_title, aBins,aMin,aMax);
		
		// set up transverse momentum
		char *num_mom_name = Form("TransMomNum%s%d", title, i);
		char *num_mom_title = Form("Numerator Transverse Momentum of Parent particle (masses = %g, %g) GeV %d", fMass1, fMass2, i);
		fNumTransvMom[i] = new TH1D(num_mom_name, num_mom_title, aBins, fpMin, fpMax);

		// set up transverse momentum
		char *den_mom_name = Form("TransMomDen%s%d", title, i);
		char *den_mom_title = Form("Denominator Transverse Momentum of Parent particle (masses = %g, %g) GeV %d", fMass1, fMass2, i);
		fDenTransvMom[i] = new TH1D(den_mom_name, den_mom_title, aBins, fpMin, fpMax);
		// to enable error bar calculation...
		fNumInvMass[i]->Sumw2();
		fDenInvMass[i]->Sumw2();
		fNumTransvMom[i]->Sumw2();
		fDenTransvMom[i]->Sumw2();
		}
/*		
  // set up numerator
  char *num_name = Form("NumInvMass%s%d", title);
  char *num_title = Form("Minv Numerator - (masses = %g, %g); M_{inv} GeV", fMass1, fMass2);
  fNumInvMass = new TH1D(num_name, num_title, aBins,aMin,aMax);

  // set up denominator
  char *den_name = Form("DenInvMass%s", title);
  char *den_title = Form("Minv Denominator - (masses = %g, %g); M_{inv} GeV", fMass1, fMass2);
  fDenInvMass = new TH1D(den_name, den_title, aBins,aMin,aMax);
  
  // set up transverse momentum
  char *num_mom_name = Form("TransMomNum%s", title);
  char *num_mom_title = Form("Numerator Transverse Momentum of Parent particle (masses = %g, %g) GeV", fMass1, fMass2);
  fNumTransvMom = new TH1D(num_mom_name, num_mom_title, aBins, 0, 5);

	// set up transverse momentum
  char *den_mom_name = Form("TransMomDen%s", title);
  char *den_mom_title = Form("Denominator Transverse Momentum of Parent particle (masses = %g, %g) GeV", fMass1, fMass2);
  fDenTransvMom = new TH1D(den_mom_name, den_mom_title, aBins, 0, 5);
  // to enable error bar calculation...
  fNumInvMass->Sumw2();
  fDenInvMass->Sumw2();
  fNumTransvMom->Sumw2();
  fDenTransvMom->Sumw2();
  */
}

//____________________________
AliFemtoCorrFctnMezonPhi::AliFemtoCorrFctnMezonPhi(const AliFemtoCorrFctnMezonPhi& aCorrFctn) :
  AliFemtoCorrFctn(aCorrFctn),
  fNumInvMass(nullptr),
  fDenInvMass(nullptr),
  fNumTransvMom(nullptr),
  fDenTransvMom(nullptr),
  fBins(aCorrFctn.fBins),
  fMin(aCorrFctn.fMin),
  fMax(aCorrFctn.fMax),
  fMass1(aCorrFctn.fMass1),
  fMass2(aCorrFctn.fMass2)
{
  // copy constructor
  for(int i = 0; i < fSections; i++){
	fNumInvMass[i] = new TH1D(*aCorrFctn.fNumInvMass[i]);
	fDenInvMass[i] = new TH1D(*aCorrFctn.fDenInvMass[i]);
	fNumTransvMom[i] = new TH1D(*aCorrFctn.fNumTransvMom[i]);
	fDenTransvMom[i] = new TH1D(*aCorrFctn.fNumTransvMom[i]);
  }
  
}

//____________________________
AliFemtoCorrFctnMezonPhi::~AliFemtoCorrFctnMezonPhi()
{
  // destructor

  delete fNumInvMass;
  delete fDenInvMass;
  delete fNumTransvMom;
  delete fDenTransvMom;
  delete fNumInvMass2D;
  delete fDenInvMass2D;
}

//_________________________
AliFemtoCorrFctnMezonPhi& AliFemtoCorrFctnMezonPhi::operator=(const AliFemtoCorrFctnMezonPhi& aCorrFctn)
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
  *fNumTransvMom = *aCorrFctn.fNumTransvMom;
  *fDenTransvMom = *aCorrFctn.fDenTransvMom;
  
  *fNumInvMass2D = *aCorrFctn.fNumInvMass2D;
  *fDenInvMass2D = *aCorrFctn.fDenInvMass2D;

  return *this;
}
//_________________________
void AliFemtoCorrFctnMezonPhi::Finish()
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
AliFemtoString AliFemtoCorrFctnMezonPhi::Report()
{
  // create report
  AliFemtoString report = "Invariant Mass Correlation Function Report:\n";
  for(int i = 0; i<fSections; i++){
	report += Form("Number of entries in numerator %d:\t%E\n", i, fNumInvMass[i]->GetEntries());
	report += Form("Number of entries in denominator %d:\t%E\n", i, fDenInvMass[i]->GetEntries());
	report += Form("Number of entries in momentum numerator %d:\t%E\n", i, fNumTransvMom[i]->GetEntries());
	report += Form("Number of entries in momentum denominator %d:\t%E\n", i, fDenTransvMom[i]->GetEntries());
	report += Form("Assumed masses %d :\t%g\t%g\n", i, fMass1, fMass2);
  }
  

  return report;
}
//____________________________
void AliFemtoCorrFctnMezonPhi::AddRealPair(AliFemtoPair* pair)
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

  double pt = TMath::Sqrt((E1+E2)*(E1+E2) - (minv*minv) - (pz1+pz2)*(pz1+pz2));
	//fNumInvMass->Fill(minv);
	//fNumTransvMom->Fill(pt);
	fNumInvMass2D->Fill(pt,minv);
  for(int i = 0; i < fSections; i++){
		if(pt > (fMin + i*fIntersection) && pt < (fMin + (i+1)*fIntersection)){
			 fNumInvMass[i]->Fill(minv);
			 fNumTransvMom[i]->Fill(pt);
		}
	}
}
//____________________________
void AliFemtoCorrFctnMezonPhi::AddMixedPair(AliFemtoPair* pair)
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
  double pt = TMath::Sqrt((E1+E2)*(E1+E2) - (minv*minv) - (pz1+pz2)*(pz1+pz2));
  //fDenInvMass->Fill(minv);
	//fDenTransvMom->Fill(minv);
	fDenInvMass2D->Fill(pt, minv);
	for(int i = 0; i < fSections; i++){
		if(pt > (fMin + i*fIntersection) && pt < (fMin + (i+1)*fIntersection)){
			fDenInvMass[i]->Fill(minv);
			fDenTransvMom[i]->Fill(minv);
		}
	}
 
}



void AliFemtoCorrFctnMezonPhi::WriteHistos()
{
  // Write out result histograms
  fNumInvMass2D->Write();
  fDenInvMass2D->Write();
  for(int i = 0; i<fSections; i++){
	fNumInvMass[i]->Write();
	fDenInvMass[i]->Write();
	fNumTransvMom[i]->Write();
	fDenTransvMom[i]->Write();
  }
  
}

TList* AliFemtoCorrFctnMezonPhi::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();
  tOutputList->Add(fNumInvMass2D);
  tOutputList->Add(fDenInvMass2D);
	for(int i = 0; i<fSections; i++){
		tOutputList->Add(fNumInvMass[i]);
		tOutputList->Add(fDenInvMass[i]);
		tOutputList->Add(fNumTransvMom[i]);
		tOutputList->Add(fDenTransvMom[i]);
	}
  

  return tOutputList;
}

void AliFemtoCorrFctnMezonPhi::SetParticleMasses(double mass1, double mass2)
{
  fMass1=mass1;
  fMass2=mass2;
  fNumInvMass2D->SetTitle(Form("Minv Numerator 2D - masses = %g, %g", fMass1, fMass2));
  fDenInvMass2D->SetTitle(Form("Minv Denominator 2D - masses = %g, %g", fMass1, fMass2));
  for(int i =0; i<fSections; i++){
	fNumInvMass[i]->SetTitle(Form("Minv Numerator %d - masses = %g, %g",i, fMass1, fMass2));
	fDenInvMass[i]->SetTitle(Form("Minv Denominator %d - masses = %g, %g",i, fMass1, fMass2));
	fNumTransvMom[i]->SetTitle(Form("Transverse Momentum Numerator %d = %g, %g",i,  fMass1, fMass2));
	fDenTransvMom[i]->SetTitle(Form("Transverse Momentum Denominator %d = %g, %g",i,  fMass1, fMass2));
  }
  
}

void AliFemtoCorrFctnMezonPhi::SetParticle1Mass(double mass)
{
  SetParticleMasses(mass, fMass2);
}

void AliFemtoCorrFctnMezonPhi::SetParticle2Mass(double mass)
{
  SetParticleMasses(fMass1, mass);
}
