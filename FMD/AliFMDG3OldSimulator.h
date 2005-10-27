#ifndef ALIFMDG3OLDSIMULATOR_H
#define ALIFMDG3OLDSIMULATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ALIFMDG3SIMULATOR
# include <AliFMDG3Simulator.h>
#endif
class AliFMD;
class AliFMDRing;
class AliFMDDetector;
class AliFMD1;
class AliFMD2;
class AliFMD3;

//____________________________________________________________________
class AliFMDG3OldSimulator : public AliFMDG3Simulator
{
public:
  AliFMDG3OldSimulator();
  /** CTOR */
  AliFMDG3OldSimulator(AliFMD* fmd, Bool_t detailed=kTRUE);
  virtual ~AliFMDG3OldSimulator() {}
  virtual void UseDivided(Bool_t) { fUseDivided = kTRUE; }
protected:
  /** Make a ring volume 
      @param r Ring geometry 
      @return  Ring volume */
  Bool_t RingGeometry(AliFMDRing* r);
  ClassDef(AliFMDG3OldSimulator,1);
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

