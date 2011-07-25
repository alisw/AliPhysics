//**************************************************************************                        
//* This file is property of and copyright by the ALICE HLT Project        *                        
//* ALICE Experiment at CERN, All rights reserved.                         *                        
//*                                                                        *                        
//* Primary Authors: marcelfigueredo@gmail.com                             *                        
//*                  for The ALICE HLT Project.                            *                        
//*                                                                        *                        
//* Permission to use, copy, modify and distribute this software and its   *                        
//* documentation strictly for non-commercial purposes is hereby granted   *                        
//* without fee, provided that the above copyright notice appears in all   *                        
//* copies and that both the copyright notice and this permission notice   *                        
//* appear in the supporting documentation. The authors make no claims     *                        
//* about the suitability of this software for any purpose. It is          *                        
//* provided "as is" without express or implied warranty.                  *                        
//**************************************************************************   


#ifndef ALIHLTTRIGGEREMCALELECTRON_H
#define ALIHLTTRIGGEREMCALELECTRON_H

#include "AliHLTTriggerEmcalClusterEnergy.h"
#include "AliHLTTrigger.h"

class AliHLTCaloClusterReader;
class TRefArray;
class AliESDEvent;
class TMap;

class AliHLTTriggerEmcalElectron : public AliHLTTrigger 
{

public:
  AliHLTTriggerEmcalElectron();
  ~AliHLTTriggerEmcalElectron();

  /// inherited from AliHLTTrigger: name of this trigger
  const char* GetTriggerName() const;
  /// inherited from AliHLTComponent: create an instance
  AliHLTComponent* Spawn();

  ///Inherited from AliHLTComponent: Get list of OCDB objects
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

  /// inherited from AliHLTComponent
  //  Get a ratio by how much the data volume is shrunken or enhanced.
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

protected :

  ///Get the clusters from the esd
  Int_t GetClustersFromEsd( const AliESDEvent * esd, TRefArray * clustersRefs );

  // FR: Set the appropriate readout list for each calo
  void SetCaloReadoutList();
  
  /// inherited from AliHLTTrigger: calculate the trigger
  Int_t DoTrigger();

  
  /// Check if cluster fullfills criteria and if so trigger
  template <class T> 
  Bool_t TriggerOnEoverP(T* cluster,AliESDEvent *esd);

  /// Threshold to trigger on EoverP
  Float_t fEThreshold;
  Float_t fEoverPThreshold;
  Float_t fEoverPLimit;

  
  ///array to hold esd clusters
  TRefArray * fClustersRefs;  //!transient

    ///Cluster data struct reader
  AliHLTCaloClusterReader * fClusterReader; //!transient

  /// the default configuration entry for this component
  const char* fOCDBEntry; //!transient
  const TString fDetector;
  
  AliHLTComponentDataType fInputDataType;   ///Input data type for calo struct input, must be set in child class
  
  
  
 private:
  
  /// Copy constructor prohibited
  AliHLTTriggerEmcalElectron(const AliHLTTriggerEmcalElectron & );
  
  /// Assignment operator prohibited
  AliHLTTriggerEmcalElectron& operator=(const AliHLTTriggerEmcalElectron &);
  
  ClassDef(AliHLTTriggerEmcalElectron, 0);

};


#endif //ALIHLTTRIGGERCALOCLUSTERENERGY_H
