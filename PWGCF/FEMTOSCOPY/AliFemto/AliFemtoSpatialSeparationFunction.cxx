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
#include "AliFemtoModelHiddenInfo.h"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__ 
ClassImp(AliFemtoSpatialSeparationFunction)
#endif

AliFemtoSpatialSeparationFunction::AliFemtoSpatialSeparationFunction(const char* title, const int numberOfBins) : AliFemtoCorrFctn()
{
  fAlphaNum = new TH1D(Form("NumSpatialSeparation_%s",title),Form("NumSpatialSeparation_%s",title),numberOfBins,0,TMath::Pi());

  fAlphaDen = new TH1D(Form("DenSpatialSeparation_%s",title),Form("DenSpatialSeparation_%s",title),numberOfBins,0,TMath::Pi());
  
  fAlphaNum->Sumw2();
  fAlphaDen->Sumw2();
  
  fFirstParticleId      = new TH1D(Form("IdFirst_%s",title), Form("IdFirst_%s",title), 6000, 0.0, 6000.0);
  fSecondParticleId     = new TH1D(Form("IdSecond_%s",title), Form("IdSecond_%s",title), 6000, 0.0, 6000.0);
  fFirstParticleOrigin  = new TH1D(Form("OriginFirst_%s",title), Form("OriginFirst_%s",title), 6000, 0.0, 6000.0);
  fSecondParticleOrigin = new TH1D(Form("OriginSecond_%s",title), Form("OriginSecond_%s",title), 6000, 0.0, 6000.0);
  
  
  
  for(int i=0;i<3;i++){
    p1real[i] = 0.0;
    p2real[i] = 0.0;
    p1mixed[i] = 0.0;
  }
}

AliFemtoSpatialSeparationFunction::AliFemtoSpatialSeparationFunction(const AliFemtoSpatialSeparationFunction& aFunction) : AliFemtoCorrFctn()
{
  if (aFunction.fAlphaNum)    fAlphaNum = new TH1D(*aFunction.fAlphaNum);
  else                        fAlphaNum = nullptr;
  
  if (aFunction.fAlphaDen)    fAlphaDen = new TH1D(*aFunction.fAlphaDen);
  else                        fAlphaDen = nullptr;
  
  if (aFunction.fFirstParticleId)    fFirstParticleId = new TH1D(*aFunction.fFirstParticleId);
  else                        fFirstParticleId = nullptr;
  
  if (aFunction.fSecondParticleId)    fSecondParticleId = new TH1D(*aFunction.fSecondParticleId);
  else                        fSecondParticleId = nullptr;
  
  if (aFunction.fFirstParticleOrigin)    fFirstParticleOrigin = new TH1D(*aFunction.fFirstParticleOrigin);
  else                        fFirstParticleOrigin = nullptr;
  
  if (aFunction.fSecondParticleOrigin)    fSecondParticleOrigin = new TH1D(*aFunction.fSecondParticleOrigin);
  else                        fSecondParticleOrigin = nullptr;
  
}

AliFemtoSpatialSeparationFunction::~AliFemtoSpatialSeparationFunction()
{
  if(fAlphaNum) delete fAlphaNum;
  if(fAlphaDen) delete fAlphaDen;
  
  if(fFirstParticleId) delete fFirstParticleId;
  if(fSecondParticleId) delete fSecondParticleId;
  if(fFirstParticleOrigin) delete fFirstParticleOrigin;
  if(fSecondParticleOrigin) delete fSecondParticleOrigin;
}

AliFemtoSpatialSeparationFunction& AliFemtoSpatialSeparationFunction::operator=(const AliFemtoSpatialSeparationFunction& aFunction)
{
  if (aFunction.fAlphaNum)  fAlphaNum = new TH1D(*aFunction.fAlphaNum);
  else                      fAlphaNum = nullptr;
  
  if (aFunction.fAlphaDen)  fAlphaNum = new TH1D(*aFunction.fAlphaDen);
  else                      fAlphaDen = nullptr;
  
  
  if (aFunction.fFirstParticleId)  fFirstParticleId = new TH1D(*aFunction.fFirstParticleId);
  else                      fFirstParticleId = nullptr;
  
  if (aFunction.fSecondParticleId)  fSecondParticleId = new TH1D(*aFunction.fSecondParticleId);
  else                      fSecondParticleId = nullptr;
  
  if (aFunction.fFirstParticleOrigin)  fFirstParticleOrigin = new TH1D(*aFunction.fFirstParticleOrigin);
  else                      fFirstParticleOrigin = nullptr;
  
  if (aFunction.fSecondParticleOrigin)  fSecondParticleOrigin = new TH1D(*aFunction.fSecondParticleOrigin);
  else                      fSecondParticleOrigin = nullptr;

  return *this;
}

void AliFemtoSpatialSeparationFunction::Finish()
{
}

AliFemtoString AliFemtoSpatialSeparationFunction::Report()
{
  string stemp = "BBbar spatial separation report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fAlphaNum->GetEntries());
  stemp += ctemp;

  AliFemtoString returnThis = stemp;
  return returnThis;
}

void AliFemtoSpatialSeparationFunction::AddFirstParticle(AliFemtoParticle *particle, bool mixing)
{
  AliFemtoLorentzVector momentum = particle->FourMomentum();
  
  if(mixing)
  {
    p1mixed[0] += momentum.x();
    p1mixed[1] += momentum.y();
    p1mixed[2] += momentum.z();
  }
  else
  {
    p1real[0] += momentum.x();
    p1real[1] += momentum.y();
    p1real[2] += momentum.z();
  }
  
  AliFemtoModelHiddenInfo *hiddenInfo = (AliFemtoModelHiddenInfo*)particle->GetHiddenInfo();
  
  if(hiddenInfo)
  {
    int partID = TMath::Abs(hiddenInfo->GetPDGPid());
    int motherID = TMath::Abs(hiddenInfo->GetMotherPdgCode());
    
    fFirstParticleId->Fill(partID);
    fFirstParticleOrigin->Fill(motherID);
  }
}

void AliFemtoSpatialSeparationFunction::AddSecondParticle(AliFemtoParticle *particle)
{
  AliFemtoLorentzVector momentum = particle->FourMomentum();
  
  p2real[0] += momentum.x();
  p2real[1] += momentum.y();
  p2real[2] += momentum.z();
  
  AliFemtoModelHiddenInfo *hiddenInfo = (AliFemtoModelHiddenInfo*)particle->GetHiddenInfo();
  
  if(hiddenInfo)
  {
    int partID = TMath::Abs(hiddenInfo->GetPDGPid());
    int motherID = TMath::Abs(hiddenInfo->GetMotherPdgCode());
    
    fSecondParticleId->Fill(partID);
    fSecondParticleOrigin->Fill(motherID);
  }
}

void AliFemtoSpatialSeparationFunction::CalculateAnglesForEvent()
{
  double mod1real = sqrt( p1real[0]*p1real[0] + p1real[1]*p1real[1] + p1real[2]*p1real[2] );
  double mod2real = sqrt( p2real[0]*p2real[0] + p2real[1]*p2real[1] + p2real[2]*p2real[2] );
  double mod1mixed = sqrt( p1mixed[0]*p1mixed[0] + p1mixed[1]*p1mixed[1] + p1mixed[2]*p1mixed[2] );
  
  
  if(fabs(mod1real) < 0.0000001 || fabs(mod2real) < 0.0000001 || fabs(mod1mixed) < 0.0000001) return;
  
  double alphaReal = acos((p1real[0]*p2real[0] + p1real[1]*p2real[1] + p1real[2]*p2real[2])/(mod1real*mod2real));
  double alphaMixed = acos((p1mixed[0]*p2real[0] + p1mixed[1]*p2real[1] + p1mixed[2]*p2real[2])/(mod1mixed*mod2real));
  
  fAlphaNum->Fill(fabs(alphaReal));
  fAlphaDen->Fill(fabs(alphaMixed));
  
  for(int i=0;i<3;i++){
    p1real[i] = 0.0;
    p2real[i] = 0.0;
    p1mixed[i] = 0.0;
  }
}

void AliFemtoSpatialSeparationFunction::WriteHistos()
{
  fAlphaNum->Write();
  fAlphaDen->Write();
  
  fFirstParticleId->Write();
  fSecondParticleId->Write();
  fFirstParticleOrigin->Write();
  fSecondParticleOrigin->Write();
}

TList* AliFemtoSpatialSeparationFunction::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fAlphaNum);
  tOutputList->Add(fAlphaDen);
  tOutputList->Add(fFirstParticleId);
  tOutputList->Add(fSecondParticleId);
  tOutputList->Add(fFirstParticleOrigin);
  tOutputList->Add(fSecondParticleOrigin);
  return tOutputList;
}

