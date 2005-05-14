#ifndef ALIITSSIMULATIONFASTPOINTSV0_H
#define ALIITSSIMULATIONFASTPOINTSV0_H

#include "AliITSsimulation.h"
/////////////////////////////////////////////////////////
//  fast simulation V0
/////////////////////////////////////////////////////////
class AliITSmodule;
class AliITSstatistics;
class TRandom;

class AliITSsimulationFastPointsV0 : public AliITSsimulation
{

public:
  AliITSsimulationFastPointsV0(); // default constructor
  AliITSsimulationFastPointsV0(const char *dataType); // standard constructor
  virtual ~AliITSsimulationFastPointsV0(); 
  void CreateFastRecPoints(AliITSmodule *mod,Int_t module,TRandom *rndm);
private:
  virtual AliITSsimulation& operator=(const AliITSsimulation &)
    {return *this;};
  AliITSsimulationFastPointsV0(const AliITSsimulationFastPointsV0 &);
  AliITSsimulationFastPointsV0 & operator=(const AliITSsimulationFastPointsV0 &)
    { return *this;}
  void AddSPD(Float_t &e,AliITSmodule *mod,Int_t trackNumber);
  void AddSDD(Float_t &e,AliITSmodule *mod,Int_t trackNumber);
  void AddSSD(Float_t &e,AliITSmodule *mod,Int_t trackNumber);

private:

  AliITSstatistics *fSx;   // pointer to AliITSstatistics class
  AliITSstatistics *fSz;   // pointer to AliITSstatistics class

  ClassDef(AliITSsimulationFastPointsV0,1) // Fast point simulator.

};

#endif
