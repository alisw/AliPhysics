#ifndef ALIFMDG3SIMULATOR_H
#define ALIFMDG3SIMULATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ALIFMDSIMULATOR
# include <AliFMDSimulator.h>
#endif
class AliFMD;
class AliFMDRing;
class AliFMDDetector;
class AliFMD1;
class AliFMD2;
class AliFMD3;

//____________________________________________________________________
class AliFMDG3Simulator : public AliFMDSimulator
{
public:
  AliFMDG3Simulator();
  /** CTOR */
  AliFMDG3Simulator(AliFMD* fmd, Bool_t detailed=kTRUE);
  virtual ~AliFMDG3Simulator() {}
  /** Register */
  virtual void DefineGeometry();
protected:
  /** Make a ring volume 
      @param r Ring geometry 
      @return  Ring volume */
  Bool_t RingGeometry(AliFMDRing* r);
  /** Make a detector volume 
      @param d Detector geometry 
      @param mother Mother volume (detector volume)
      @param zmother Z position of mother 
      @param inner Inner ring volume 
      @param outer Outer ring volume 
      @return  Detector volume */
  Bool_t DetectorGeometry(AliFMDDetector* d, Double_t zmother);
  /** Make FMD1 volume 
      @param d Detector geometry 
      @param inner Inner ring volume 
      @return FMD1 volume  */
  Bool_t FMD1Geometry(AliFMD1* d);
  /** Make FMD2 volume 
      @param d Detector geometry 
      @param inner Inner ring volume 
      @param outer Outer ring volume 
      @return FMD2 volume  */
  Bool_t FMD2Geometry(AliFMD2* d);
  /** Make FMD3 volume 
      @param d Detector geometry 
      @param inner Inner ring volume 
      @param outer Outer ring volume 
      @return FMD3 volume  */
  Bool_t FMD3Geometry(AliFMD3* d);

  ClassDef(AliFMDG3Simulator,1);
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

