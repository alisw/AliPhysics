////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoSpatialSeparationFunction - A calss which calculates          //
// total momentum of baryons and antibaryons in given event and builds        //
// a histogram of anges between those vectors.                                //
//                                                                            //
// Authors: Jeremi Niedziela jeremi.niedziela@cern.ch                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoSpatialSeparationFunction.h"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__ 
ClassImp(AliFemtoSpatialSeparationFunction)
#endif

AliFemtoSpatialSeparationFunction::AliFemtoSpatialSeparationFunction(const char* title, const int numberOfBins) : AliFemtoCorrFctn()
{
  fAlpha = new TH1D(Form("SpatialSeparation_%s",title),Form("SpatialSeparation_%s",title),numberOfBins,0,TMath::Pi());

  fAlpha->Sumw2();
  
  for(int i=0;i<3;i++){
    p1[i] = 0.0;
    p2[i] = 0.0;
  }
}

AliFemtoSpatialSeparationFunction::AliFemtoSpatialSeparationFunction(const AliFemtoSpatialSeparationFunction& aFunction) : AliFemtoCorrFctn()
{
  if (aFunction.fAlpha)   fAlpha = new TH1D(*aFunction.fAlpha);
  else                    fAlpha = nullptr;
}

AliFemtoSpatialSeparationFunction::~AliFemtoSpatialSeparationFunction()
{
  if(fAlpha)    delete fAlpha;
}

AliFemtoSpatialSeparationFunction& AliFemtoSpatialSeparationFunction::operator=(const AliFemtoSpatialSeparationFunction& aFunction)
{
  if (aFunction.fAlpha)   fAlpha = new TH1D(*aFunction.fAlpha);
  else                        fAlpha = nullptr;
  
  return *this;
}

void AliFemtoSpatialSeparationFunction::Finish()
{
}

AliFemtoString AliFemtoSpatialSeparationFunction::Report()
{
  string stemp = "BBbar spatial separation report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fAlpha->GetEntries());
  stemp += ctemp;

  AliFemtoString returnThis = stemp;
  return returnThis;
}

void AliFemtoSpatialSeparationFunction::AddFirstParticle(AliFemtoParticle *particle)
{
  AliFemtoLorentzVector momentum = particle->FourMomentum();
  
  p1[0] += momentum.x();
  p1[1] += momentum.y();
  p1[2] += momentum.z();
}

void AliFemtoSpatialSeparationFunction::AddSecondParticle(AliFemtoParticle *particle)
{
  AliFemtoLorentzVector momentum = particle->FourMomentum();
  
  p2[0] += momentum.x();
  p2[1] += momentum.y();
  p2[2] += momentum.z();
}

void AliFemtoSpatialSeparationFunction::CalculateAnglesForEvent()
{
  double mod1 = sqrt( p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2] );
  double mod2 = sqrt( p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2] );
  
  if(fabs(mod1) < 0.0000001 || fabs(mod2) < 0.0000001) return;
  
  double alpha = acos((p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2])/(mod1*mod2));
  
  fAlpha->Fill(fabs(alpha));
  
  for(int i=0;i<3;i++){
    p1[i] = 0.0;
    p2[i] = 0.0;
  }
}

void AliFemtoSpatialSeparationFunction::WriteHistos()
{
  fAlpha->Write();
}

TList* AliFemtoSpatialSeparationFunction::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fAlpha);
  return tOutputList;
}

