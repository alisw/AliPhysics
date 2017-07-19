////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoSpatialSeparationFunction - A calss which calculates          //
// total momentum of baryons and antibaryons in given event and builds        //
// a histogram of anges between those vectors.                                //
//                                                                            //
// Authors: Jeremi Niedziela jeremi.niedziela@cern.ch                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoSpatialSeparationFunction_H
#define AliFemtoSpatialSeparationFunction_H

#include "TH1D.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoSpatialSeparationFunction : public AliFemtoCorrFctn
{
public:
  AliFemtoSpatialSeparationFunction(const char* title, const int numberOfBins=50);
  AliFemtoSpatialSeparationFunction(const AliFemtoSpatialSeparationFunction& aFunction);
  AliFemtoSpatialSeparationFunction& operator=(const AliFemtoSpatialSeparationFunction& aFunction);
  virtual ~AliFemtoSpatialSeparationFunction();

  virtual AliFemtoString Report();
  
  void AddFirstParticle(AliFemtoParticle* particle, bool mixing=false);
  void AddSecondParticle(AliFemtoParticle* particle);

  void CalculateAnglesForEvent();
  
  virtual void Finish();
  void WriteHistos();
  virtual TList* GetOutputList();
  
private:
  TH1D *fAlphaNum; // distibution of angles between momenta of particle collections from the same event
  TH1D *fAlphaDen; // distibution of angles between momenta of particle collections from different events
  
  double p1real[3]; // sum of momenta of particles of first type
  double p2real[3]; // sum of momenta of particles of second type
  
  double p1mixed[3]; // sum of momenta of particles of first type (other event)
  
#ifdef __ROOT__
  ClassDef(AliFemtoSpatialSeparationFunction, 1)
#endif
};


#endif
