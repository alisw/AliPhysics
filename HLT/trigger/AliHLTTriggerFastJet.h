//**************************************************************************                        
//* This file is property of and copyright by the ALICE HLT Project        *                        
//* ALICE Experiment at CERN, All rights reserved.                         *                        
//*                                                                        *                        
//* Primary Authors: leonidas.xaplanteris@gmail.com                        *                        
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

#ifndef ALIHLTTRIGGERFASTJET_H
#define ALIHLTTRIGGERFASTJET_H

#include "AliHLTTrigger.h"
#include "AliFJWrapper.h"

class TRefArray;
class AliESDEvent;
class TMap;

class AliHLTTriggerFastJet : public AliHLTTrigger
{
  
 public:
  AliHLTTriggerFastJet(TString detector);
  ~AliHLTTriggerFastJet();
  
  ///
  AliHLTComponent* Spawn();
  
  // inherited from AliHLTTrigger : name of this trigger
  const char* GetTriggerName() const;
  
  // inherited from AliHLTComponent : get list of OCDB objects
  void GetOCDBObjectDescription( TMap* const targetMap);

 protected:
  // inherited from AliHLTComponent: initialization
  int DoInit(int argc, const char** argv);
  
  // inherited from AliHLTComponent: de-initialization
  int DoDeInit();
  
  // inherited from AliHLTComponent: re-configuration
  int Reconfigure(const char* cdbEntry, const char* chainId);

  // inherited from AliHLTComponent: scan one argument and its parameters
  int ScanConfigurationArgument(int argc, const char** argv);

  // inherited from AliHLTComponent: Get a ratio by how much the data volume is shrunken or enhanced.
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

protected :
  // get the clusters from the esd
  Int_t GetClustersFromEsd( const AliESDEvent * esd, TRefArray * clustersRefs );

  // inherited from AliHLTTrigger: calculate the trigger
  Int_t DoTrigger();

  // default constructor prohibited
  AliHLTTriggerFastJet();


 
  // check if cluster fullfills criteria and if so trigger
  template <class T> 
  Bool_t TriggerOnJet(T Jet);

  // threshold jet energy to trigger on
  Float_t fEThreshold;

  // the detector string (PHOS or EMCAL)
  const TString fDetector;

  // fast jet wrapper
  AliFJWrapper *fFastJetWrapper;

  // array to hold esd clusters
  TRefArray * fClustersRefs;  //!transient

  // the default configuration entry for this component
  const char* fOCDBEntry; //!transient

  
  // offset for calo tracks
  Int_t fOffset;

  // input data type for calo struct input, must be set in child class
  AliHLTComponentDataType fInputDataType;   

 private:

  // copy constructor prohibited
  AliHLTTriggerFastJet(const AliHLTTriggerFastJet & );
  
  // assignment operator prohibited
  AliHLTTriggerFastJet& operator=(const AliHLTTriggerFastJet &);
  
  ClassDef(AliHLTTriggerFastJet, 0)
};

#endif //ALIHLTTRIGGERFASTJET_H
