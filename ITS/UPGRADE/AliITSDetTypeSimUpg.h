#ifndef ALIITSDETTYPESIMUPG_H
#define ALIITSDETTYPESIMUPG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id: AliITSDetTypeSimUpg.h 53025 2011-11-19 22:50:51Z masera $ 
*/

/////////////////////////////////////////////////////////////////////////
// * This class contains all of the "external" information needed to do//
// * detector specific simulations for the ITS.                        //
/////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include "AliITSLoader.h"
#include "AliITSSimuParam.h"
#include "AliITSDetTypeSim.h"

class TObjArray;
class TClonesArray;
class TTree;
class AliCDBMetaData;
class AliITSdigit;
class AliITSdigitPixUpg;
class AliITSmodule;
class AliITSpListItem;
class AliITSsimulation;
class AliITSsegmentation;
class AliITSresponse;
class AliITSCalibration;
class AliITSgeom;
class AliITSTriggerConditions;


class AliITSDetTypeSimUpg : public AliITSDetTypeSim {
 public:
  //
  enum {kDetPixUpg, kNDetTypes};
  //  
  AliITSDetTypeSimUpg();
  virtual ~AliITSDetTypeSimUpg(); 
  
  virtual AliITSCalibration* GetCalibrationModel(Int_t iMod) const;
  virtual AliITSTriggerConditions* GetTriggerConditions() {return 0;} // tmp
  
  virtual void SetDefaults();
  virtual void SetDefaultSimulation();
  //
  virtual void AddSimDigit(Int_t branch, const AliITSdigit *d);
  virtual void AddSimDigit(Int_t branch,Float_t phys,Int_t* digits,
			   Int_t* tracks,Int_t *hits,Float_t* trkcharges,
			   Int_t sigexpanded=-1000);

 protected:
  virtual Bool_t GetCalibration();
    
 private:
  AliITSDetTypeSimUpg(const AliITSDetTypeSimUpg &source);
  AliITSDetTypeSimUpg& operator=(const AliITSDetTypeSimUpg &source);
  void SetDefaultSegmentation(Int_t idet);  // creates def segm.
  //
  Int_t              fDetNModules[kNDetTypes]; // Total numbers of modules of each type, RS
  //       
  static const char* fgkDetTypeName[kNDetTypes]; // detector names
  //
  ClassDef(AliITSDetTypeSimUpg,1) // ITS Upg Simulation structure
 
};

#endif
