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
  
  void AddFirstParticle(AliFemtoParticle* particle);
  void AddSecondParticle(AliFemtoParticle* particle);

  void CalculateAnglesForEvent();
  
  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
private:
  
  TH1D *fAlpha;          // angles between momenta of particle collections
  
  double p1[3];
  double p2[3];
  
#ifdef __ROOT__
  ClassDef(AliFemtoSpatialSeparationFunction, 1)
#endif
};


#endif
