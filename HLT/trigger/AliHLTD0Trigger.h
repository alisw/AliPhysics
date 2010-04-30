//-*- Mode: C++ -*-
// $Id: AliHLTD0Trigger.h 
#ifndef ALIHLTD0TRIGGER_H
#define ALIHLTD0TRIGGER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTD0Trigger.h
/// @author Gaute Ovrebekk
/// @date   2009-10-28
/// @brief  HLT trigger component for D0->Kpi

#include "AliHLTTrigger.h"
#include <vector>
#include "AliHLTD0toKpi.h"

class TH1F;
class TObjArray;
class AliESDVertex;
class AliExternalTrackParam;
class AliHLTMCEvent;

/**
 * @class  AliHLTD0Trigger
 *
 * HLT trigger component for D0->Kpi
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b D0Trigger                             <br>
 * Library: \b libAliHLTTrigger.so                                        <br>
 * Input Data Types:  kAliHLTDataTypeESDObject, kAliHLTDataTypeESDTree
 *                    kAliHLTDataTypeTrack                                <br>
 * Output Data Types: ::kAliHLTAnyDataType                                <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -pt    <i> pt cut for decay products </i> <br>
 * \li -dca    <i> dca cut for distance between decay tracks  </i> <br>
 * \li -invmass    <i> inv. mass half width of D0  </i> <br>
 * \li -costhetastar    <i> cos of decay angle  </i> <br>
 * \li -d0    <i> Impact parameter for decay products  </i> <br>
 * \li -d0d0    <i> Product of impact parameter for decay products  </i> <br>
 * \li -cospoint    <i> pointing angle  </i> <br>
 * \li -plothistogram    <i> ploting the inv. mass and pt of D0  </i> <br>
 * \li -useKF    <i> will use partilce KF for vertex finding  </i> <br>
 *
 * By default, configuration is loaded from OCDB, can be overridden by
 * component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * HLT/ConfigHLT/D0Trigger: TObjString storing the arguments
 *
 * <h2>Performance:</h2>
 * 
 *
 * <h2>Memory consumption:</h2>
 * 
 *
 * <h2>Output size:</h2>
 * 
 *
 * \ingroup alihlt_trigger_components
 */
class AliHLTD0Trigger : public AliHLTTrigger
{
 public:
  AliHLTD0Trigger();
  ~AliHLTD0Trigger();

  /// inherited from AliHLTTrigger: name of this trigger
  virtual const char* GetTriggerName() const;
  /// inherited from AliHLTComponent: create an instance
  virtual AliHLTComponent* Spawn();

  /// inherited from AliHLTComponent: return OCDB requirements
  void GetOCDBObjectDescription( TMap* const targetMap);

 protected:
  /// inherited from AliHLTComponent: handle the initialization
  int DoInit(int argc, const char** argv);

  /// inherited from AliHLTComponent: handle cleanup
  int DoDeinit();

  /// inherited from AliHLTComponent: handle re-configuration event
  int Reconfigure(const char* cdbEntry, const char* chainId);

  /// inherited from AliHLTComponent, scan one argument and
  /// its parameters
  int ScanConfigurationArgument(int argc, const char** argv);

 private:
  /// Not implemented. Do not allow copying of this object.
  AliHLTD0Trigger(const AliHLTD0Trigger& );
  /// Not implemented. Do not allow copying of this object.
  AliHLTD0Trigger& operator=(const  AliHLTD0Trigger& );

  /// inherited from AliHLTTrigger: calculate the trigger
  virtual int DoTrigger();
  
  void SingleTrackSelect(AliExternalTrackParam*);
  Int_t RecV0(const TObject* iter);
  void RecD0(Int_t&,Int_t&);
  bool CheckTrackMC(AliExternalTrackParam* pt, AliExternalTrackParam* pn);

  /// pt cut for decay, minimum [GeV/c]
  float fPtMin;                                            //! transient
  /// Distance between decay tracks [cm] ??
  float fdca;                                              //! transient
  /// Inv. mass half width [GeV]
  float finvMass;                                          //! transient
  /// Decay angle
  float fcosThetaStar;                                     //! transient  
  /// Distance from primary vertex for decay tracks [cm] 
  float fd0;                                               //! transient
  /// Product of d0 for the two decay tracks [cm^2]
  float fd0d0;                                             //! transient
  /// Pionting angle
  float fcosPoint;                                         //! transient 

  bool fplothisto;                                         //! transient 
  bool fUseV0;                                             //! transient 

  Double_t mD0PDG;                                         //! transient

  /// D0 inv. mass plot
  TH1F *fD0mass;                                           //! transient  
  TH1F *fD0pt;                                             //! transient  

  vector<AliExternalTrackParam*> fPos;                       //! transient
  vector<AliExternalTrackParam*> fNeg;                       //! transient

  AliHLTD0toKpi *fd0calc;                                   //! transient
  TObjArray *ftwoTrackArray;                                //! transient

  Int_t fTotalD0;                                           //! transient
  Int_t fTotalD0true;                                       //! transient
  AliESDVertex *fVertex;                                    //! transient
  Double_t fField;                                          //!transient
  
  AliHLTMCEvent* fEvent;                                    //!transient

  bool fuseKF;                                               //!transient

  /// the default configuration entry for this component
  static const char* fgkOCDBEntry; //!transient

  ClassDef(AliHLTD0Trigger, 0)
};
#endif //ALIHLTD0TRIGGER_H
