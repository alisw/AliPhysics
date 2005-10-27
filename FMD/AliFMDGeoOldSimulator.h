#ifndef ALIFMDGEOOLDSIMULATOR_H
#define ALIFMDGEOOLDSIMULATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ALIFMDGEOSIMULATOR
# include <AliFMDGeoSimulator.h>
#endif
class AliFMD;
class AliFMDRing;
class AliFMDDetector;
class AliFMD1;
class AliFMD2;
class AliFMD3;

//____________________________________________________________________
class AliFMDGeoOldSimulator : public AliFMDGeoSimulator
{
public:
  AliFMDGeoOldSimulator();
  /** CTOR */
  AliFMDGeoOldSimulator(AliFMD* fmd, Bool_t detailed=kTRUE);
  virtual ~AliFMDGeoOldSimulator() {}
  virtual void UseDivided(Bool_t) { fUseDivided = kTRUE; }
protected:
  /** Make a ring volume 
      @param r Ring geometry 
      @return  Ring volume */
  TGeoVolume* RingGeometry(AliFMDRing* r);
  ClassDef(AliFMDGeoOldSimulator,1);
};


#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//

