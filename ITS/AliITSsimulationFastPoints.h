#ifndef ALIITSSIMULATIONFASTPOINTS_H
#define ALIITSSIMULATIONFASTPOINTS_H

#include "AliITSsimulation.h"

class AliITSmodule;
class AliITSstatistics;

class AliITSsimulationFastPoints : public AliITSsimulation
{

public:
  AliITSsimulationFastPoints(); // default constructor
  virtual ~AliITSsimulationFastPoints(); 
  void CreateFastRecPoints(AliITSmodule *mod);
private:
  void AddSPD(Float_t &e,AliITSmodule *mod,Int_t trackNumber);
  void AddSDD(Float_t &e,AliITSmodule *mod,Int_t trackNumber);
  void AddSSD(Float_t &e,AliITSmodule *mod,Int_t trackNumber);

private:

  AliITSstatistics *fSx;   // pointer to AliITSstatistics class
  AliITSstatistics *fSz;   // pointer to AliITSstatistics class

  ClassDef(AliITSsimulationFastPoints,1) // Fast point simulator.

};

#endif
