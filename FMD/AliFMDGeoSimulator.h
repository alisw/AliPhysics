#ifndef ALIFMDGEOSIMULATOR_H
#define ALIFMDGEOSIMULATOR_H
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
class TGeoVolume;
class TGeoMedium;
class AliFMD;
class AliFMDRing;
class AliFMDDetector;
class AliFMD1;
class AliFMD2;
class AliFMD3;

//____________________________________________________________________
class AliFMDGeoSimulator : public AliFMDSimulator
{
public:
  AliFMDGeoSimulator();
  /** CTOR */
  AliFMDGeoSimulator(AliFMD* fmd, Bool_t detailed=kTRUE);
  virtual ~AliFMDGeoSimulator() {}
  /** Initialize */
  virtual void DefineMaterials();
  /** Register */
  virtual void DefineGeometry();
protected:
  /** Make a ring volume 
      @param r Ring geometry 
      @return  Ring volume */
  virtual TGeoVolume* RingGeometry(AliFMDRing* r);
  /** Make a detector volume 
      @param d Detector geometry 
      @param mother Mother volume (detector volume)
      @param zmother Z position of mother 
      @param inner Inner ring volume 
      @param outer Outer ring volume 
      @return  Detector volume */
  virtual TGeoVolume* DetectorGeometry(AliFMDDetector* d, 
				       TGeoVolume* mother, 
				       Double_t zmother, 
				       TGeoVolume* inner, 
				       TGeoVolume* outer=0);
  /** Make FMD1 volume 
      @param d Detector geometry 
      @param inner Inner ring volume 
      @return FMD1 volume  */
  virtual TGeoVolume* FMD1Geometry(AliFMD1* d, TGeoVolume* inner);
  /** Make FMD2 volume 
      @param d Detector geometry 
      @param inner Inner ring volume 
      @param outer Outer ring volume 
      @return FMD2 volume  */
  virtual TGeoVolume* FMD2Geometry(AliFMD2* d, TGeoVolume* inner, 
			   TGeoVolume* outer);
  /** Make FMD3 volume 
      @param d Detector geometry 
      @param inner Inner ring volume 
      @param outer Outer ring volume 
      @return FMD3 volume  */
  virtual TGeoVolume* FMD3Geometry(AliFMD3* d, TGeoVolume* inner, 
			   TGeoVolume* outer);
  TGeoMedium* fSi;	 //! Si Medium
  TGeoMedium* fC;	 //! C Medium
  TGeoMedium* fAl;	 //! Al Medium
  TGeoMedium* fPCB;	 //! PCB Medium
  TGeoMedium* fChip;	 //! Chip Medium
  TGeoMedium* fAir;	 //! Air Medium
  TGeoMedium* fPlastic;	 //! Plastic Medium
  TGeoMedium* fCopper;	 //! Copper Medium

  ClassDef(AliFMDGeoSimulator,1)
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

