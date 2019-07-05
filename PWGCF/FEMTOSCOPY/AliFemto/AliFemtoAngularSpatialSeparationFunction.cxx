////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoAngularSpatialSeparationFunction - A calss which calculates        //
// average theta and phi of baryons and antibaryons in given event and        //
// builds a histogram of angles difference.                                   //
//                                                                            //
// Authors: Jeremi Niedziela jeremi.niedziela@cern.ch                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoAngularSpatialSeparationFunction.h"
#include "AliFemtoModelHiddenInfo.h"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoAngularSpatialSeparationFunction)
#endif

AliFemtoAngularSpatialSeparationFunction::AliFemtoAngularSpatialSeparationFunction(const char* title, const int numberOfBins) : AliFemtoCorrFctn()
{
  fAlphaNum = new TH1D(Form("NumSpatialSeparation_%s",title),Form("NumSpatialSeparation_%s",title),numberOfBins,0,TMath::Pi());

  fAlphaDen = new TH1D(Form("DenSpatialSeparation_%s",title),Form("DenSpatialSeparation_%s",title),numberOfBins,0,TMath::Pi());

  fAlphaNum->Sumw2();
  fAlphaDen->Sumw2();
}

AliFemtoAngularSpatialSeparationFunction::AliFemtoAngularSpatialSeparationFunction(const AliFemtoAngularSpatialSeparationFunction& aFunction) : AliFemtoCorrFctn()
{
  if (aFunction.fAlphaNum)    fAlphaNum = new TH1D(*aFunction.fAlphaNum);
  else                        fAlphaNum = nullptr;

  if (aFunction.fAlphaDen)    fAlphaDen = new TH1D(*aFunction.fAlphaDen);
  else                        fAlphaDen = nullptr;
}

AliFemtoAngularSpatialSeparationFunction::~AliFemtoAngularSpatialSeparationFunction()
{
  if(fAlphaNum) delete fAlphaNum;
  if(fAlphaDen) delete fAlphaDen;
}

AliFemtoAngularSpatialSeparationFunction& AliFemtoAngularSpatialSeparationFunction::operator=(const AliFemtoAngularSpatialSeparationFunction& aFunction)
{
  if (aFunction.fAlphaNum)  fAlphaNum = new TH1D(*aFunction.fAlphaNum);
  else                      fAlphaNum = nullptr;

  if (aFunction.fAlphaDen)  fAlphaNum = new TH1D(*aFunction.fAlphaDen);
  else                      fAlphaDen = nullptr;

  return *this;
}

void AliFemtoAngularSpatialSeparationFunction::Finish()
{
}

AliFemtoString AliFemtoAngularSpatialSeparationFunction::Report()
{
  AliFemtoString report = "BBbar spatial separation report:\n";
  report += Form("Number of entries in numerator:\t%E\n",fAlphaNum->GetEntries());

  return report;
}

void AliFemtoAngularSpatialSeparationFunction::AddFirstParticle(AliFemtoParticle *particle, bool mixing)
{
  AliFemtoLorentzVector fourMomentum = particle->FourMomentum();

  if(mixing){
    phi1mixed.push_back(fourMomentum.Phi());
    theta1mixed.push_back(fourMomentum.Theta());
  }
  else{
    phi1real.push_back(fourMomentum.Phi());
    theta1real.push_back(fourMomentum.Theta());
  }
}

void AliFemtoAngularSpatialSeparationFunction::AddSecondParticle(AliFemtoParticle *particle)
{
  AliFemtoLorentzVector fourMomentum = particle->FourMomentum();

  phi2real.push_back(fourMomentum.Phi());
  theta2real.push_back(fourMomentum.Theta());
}

void AliFemtoAngularSpatialSeparationFunction::CalculateAnglesForEvent()
{
  double avgPhi1real=0;
  double avgPhi2real=0;
  double avgPhi1mixed=0;

  double avgTheta1real=0;
  double avgTheta2real=0;
  double avgTheta1mixed=0;

  for(UInt_t i=0;i<phi1real.size();i++){
    avgPhi1real+=phi1real[i];
  }
  avgPhi1real /= phi1real.size();

  for(UInt_t i=0;i<phi2real.size();i++){
    avgPhi2real+=phi2real[i];
  }
  avgPhi2real /= phi2real.size();

  for(UInt_t i=0;i<phi1mixed.size();i++){
    avgPhi1mixed+=phi1mixed[i];
  }
  avgPhi1mixed /= phi1mixed.size();

  for(UInt_t i=0;i<theta1real.size();i++){
    avgTheta1real+=theta1real[i];
  }
  avgTheta1real /= theta1real.size();

  for(UInt_t i=0;i<theta2real.size();i++){
    avgTheta2real+=theta2real[i];
  }
  avgTheta2real /= theta2real.size();

  for(UInt_t i=0;i<theta1mixed.size();i++){
    avgTheta1mixed+=theta1mixed[i];
  }
  avgTheta1mixed /= theta1mixed.size();


  double d1real[3];
  double d2real[3];
  double d1mixed[3];

  d1real[0] = sin(avgTheta1real) * cos(avgPhi1real);
  d1real[1] = sin(avgTheta1real) * sin(avgPhi1real);
  d1real[2] = cos(avgTheta1real);

  d2real[0] = sin(avgTheta2real) * cos(avgPhi2real);
  d2real[1] = sin(avgTheta2real) * sin(avgPhi2real);
  d2real[2] = cos(avgTheta2real);

  d1mixed[0] = sin(avgTheta1mixed) * cos(avgPhi1mixed);
  d1mixed[1] = sin(avgTheta1mixed) * sin(avgPhi1mixed);
  d1mixed[2] = cos(avgTheta1mixed);

  double alphaReal  = acos(d1real[0]*d2real[0] + d1real[1]*d2real[1] + d1real[2]*d2real[2]);
  double alphaMixed = acos(d1mixed[0]*d2real[0] + d1mixed[1]*d2real[1] + d1mixed[2]*d2real[2]);

  fAlphaNum->Fill(fabs(alphaReal));
  fAlphaDen->Fill(fabs(alphaMixed));

  phi1real.clear();
  phi2real.clear();
  phi1mixed.clear();

  theta1real.clear();
  theta2real.clear();
  theta1mixed.clear();

}

void AliFemtoAngularSpatialSeparationFunction::WriteHistos()
{
  fAlphaNum->Write();
  fAlphaDen->Write();
}

TList* AliFemtoAngularSpatialSeparationFunction::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fAlphaNum);
  tOutputList->Add(fAlphaDen);
  return tOutputList;
}
