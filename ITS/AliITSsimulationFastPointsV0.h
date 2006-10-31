#ifndef ALIITSSIMULATIONFASTPOINTSV0_H
#define ALIITSSIMULATIONFASTPOINTSV0_H

#include "AliITSsimulation.h"
/////////////////////////////////////////////////////////
//  fast simulation V0
/////////////////////////////////////////////////////////
class AliITSmodule;
class AliITSstatistics;
class TRandom;
class TClonesArray;

class AliITSsimulationFastPointsV0 : public AliITSsimulation
{

public:
  AliITSsimulationFastPointsV0(); // default constructor
  AliITSsimulationFastPointsV0(const char *dataType); // standard constructor
  AliITSsimulationFastPointsV0(const AliITSsimulationFastPointsV0 &rec);
  AliITSsimulationFastPointsV0 & operator=(const AliITSsimulationFastPointsV0 &);  
  virtual AliITSsimulation& operator=(const AliITSsimulation &)
    {return *this;};

  virtual ~AliITSsimulationFastPointsV0(); 
  void CreateFastRecPoints(AliITSmodule *mod,Int_t module,TRandom *rndm,
			   TClonesArray* recp);
private:

  void AddSPD(Float_t &e,AliITSmodule *mod,Int_t trackNumber,TClonesArray* recp);
  void AddSDD(Float_t &e,AliITSmodule *mod,Int_t trackNumber,TClonesArray* recp);
  void AddSSD(Float_t &e,AliITSmodule *mod,Int_t trackNumber,TClonesArray* recp);

private:

  Int_t fNrecp;            //current number of  fast point
  AliITSstatistics *fSx;   // pointer to AliITSstatistics class
  AliITSstatistics *fSz;   // pointer to AliITSstatistics class

  ClassDef(AliITSsimulationFastPointsV0,2) // Fast point simulator.

};

#endif
