////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoAngularSpatialSeparationFunction - A calss which calculates        //
// average theta and phi of baryons and antibaryons in given event and        //
// builds a histogram of angles difference.                                   //
//                                                                            //
// Authors: Jeremi Niedziela jeremi.niedziela@cern.ch                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoAngularSpatialSeparationFunction_H
#define AliFemtoAngularSpatialSeparationFunction_H

#include "TH1D.h"
#include "AliFemtoCorrFctn.h"

#include <vector>

class AliFemtoAngularSpatialSeparationFunction : public AliFemtoCorrFctn
{
public:
  AliFemtoAngularSpatialSeparationFunction(const char* title, const int numberOfBins=50);
  AliFemtoAngularSpatialSeparationFunction(const AliFemtoAngularSpatialSeparationFunction& aFunction);
  AliFemtoAngularSpatialSeparationFunction& operator=(const AliFemtoAngularSpatialSeparationFunction& aFunction);
  virtual ~AliFemtoAngularSpatialSeparationFunction();

  virtual AliFemtoString Report();

  void AddFirstParticle(AliFemtoParticle* particle, bool mixing=false);
  void AddSecondParticle(AliFemtoParticle* particle);

  void CalculateAnglesForEvent();

  virtual void Finish();
  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoAngularSpatialSeparationFunction(*this); }

private:
  TH1D *fAlphaNum; // distibution of angles between direction of particle collections from the same event
  TH1D *fAlphaDen; // distibution of angles between direction of particle collections from different events

  std::vector<double> phi1real;   // vector of all phi angles for particles of the first type
  std::vector<double> phi2real;   // vector of all phi angles for particles of the second type
  std::vector<double> phi1mixed;  // vector of all phi angles for particles of the first type (other event)

  std::vector<double> theta1real; // vector of all theta angles for particles of the first type
  std::vector<double> theta2real; // vector of all theta angles for particles of the second type
  std::vector<double> theta1mixed;// vector of all theta angles for particles of the first type (other event)

#ifdef __ROOT__
  ClassDef(AliFemtoAngularSpatialSeparationFunction, 1)
#endif
};


#endif
