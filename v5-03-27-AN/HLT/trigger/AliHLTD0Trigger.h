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
class TTree;
class TClonesArray;

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
 * TODO: code audit 2010-07-23 describe component arguments
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
 * \li -useV0    <i> will use the V0's found by the vertexer and stored in the ESD</i> <br>
 * \li -useKF    <i> will use partilce KF for vertex finding  </i> <br>
 * \li -send-candidates    <i> will send out an array of candidates for each event</i> <br>
 * \li -write-file    <i> will store a local file. Only use for small local tests.</i> <br>
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
  
  /// Adding single track cut on input tracks, and split in pos. and neg.
  void SingleTrackSelect(AliExternalTrackParam*,const AliESDVertex*,Double_t field);
  /// Useing the V0's in the ESD found by the V0 finder
  Int_t RecV0(const TObject* iter);
  /// Reconstructing the D0 from K and pi
  void RecD0(Int_t&,const AliESDVertex*,Double_t field);

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

  /// Option for ploting InvMass and Pt of D0's
  bool fplothisto;                                         //! transient
  /// Option for useing the V0' from V0 finder
  bool fUseV0;                                             //! transient 

  /// D0 mass
  Double_t fD0PDG;                                         //! transient

  /// D0 inv. mass plot
  TH1F *fD0mass;                                           //! transient  
  /// Pt plot of D0's
  TH1F *fD0pt;                                             //! transient  

  /// Vector for positive tracks
  vector<AliExternalTrackParam*> fPos;                     //! transient
  /// Vector for negative tracks
  vector<AliExternalTrackParam*> fNeg;                     //! transient

  /// Object for calculations
  AliHLTD0toKpi *fd0calc;                                  //! transient
  /// Array of the two decay products
  TObjArray *ftwoTrackArray;                               //! transient

  /// Counters for D0
  Int_t fTotalD0;                                          //! transient

  /// Option for useing KF particle for vertexing
  bool fuseKF;                                             //!transient

  /// Option for storing MC information
  bool fSendCandidate;                                     //!transient
  /// Tree for storing the MC information
  TTree *fCandidateTree;                                   //!transient
  /// Array with D0 candidates
  TClonesArray *fCandidateArray;                           //!transient
  /// Option for writing MC information to file
  bool fWriteFile;                                         //!transient

  /// the default configuration entry for this component
  static const char* fgkOCDBEntry; //!transient

  ClassDef(AliHLTD0Trigger, 0)
};
#endif //ALIHLTD0TRIGGER_H
