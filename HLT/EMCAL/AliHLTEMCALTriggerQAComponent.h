/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Salvatore Aiola                                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTEMCALTRIGGERQACOMPONENT_H
#define ALIHLTEMCALTRIGGERQACOMPONENT_H

/**
 * @file   AliHLTEMCALTriggerQAComponent.h
 * @author Salvatore Aiola
 * @date   Nov. 2, 2015
 * @brief  A trigger QA component for EMCAL HLT
 */

#include "AliHLTCaloProcessor.h"

class AliEMCALTriggerQA;
class AliEMCALTriggerPatchInfo;
class AliEMCALTriggerBitConfig;
class AliEMCALGeometry;
class AliEMCALTriggerFastOR;
class AliHLTCaloTriggerPatchDataStruct;
class AliHLTEMCALGeometry;
class AliHLTEMCALCaloCells;

/**
 * @class AliHLTEMCALTriggerQAComponent
 * @brief HLT component for EMCAL/DCAL trigger QA
 */
class AliHLTEMCALTriggerQAComponent : public AliHLTCaloProcessor {
public:

  enum EBeamType_t {
    kPP,
    kPbPb
  };

  AliHLTEMCALTriggerQAComponent();
  virtual ~AliHLTEMCALTriggerQAComponent();

  void SetTriggerBitConfig(const AliEMCALTriggerBitConfig* const config) { fTriggerBitConfig = config; }
	
	/**
	 *
	 * @param evtData
	 * @param blocks
	 * @param trigData
	 * @param outputPtr
	 * @param size
	 * @param outputBlocks
	 * @return
	 */
  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
	      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
	      std::vector<AliHLTComponentBlockData>& outputBlocks);

  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();

  /** interface function, see @ref AliHLTComponent for description */
  void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();

  /** interface function, see @ref AliHLTComponent for description */
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

  bool CheckInputDataType(const AliHLTComponentDataType &datatype);

  void SetPbPb2015TriggerClasses();
  void SetPP2016TriggerClasses();

protected:
  /** interface function, see @ref AliHLTComponent for description */
   int DoInit(int argc, const char** argv);

   /** interface function, see @ref AliHLTComponent for description */
   virtual int Deinit();

   /** retrieve fired trigger classes for the current event **/
   int RetrieveFiredTriggerClasses();

   /** process trigger patches contained in block **/
   //void ProcessCells(const AliHLTComponentBlockData* block, AliHLTEMCALCaloCells& cells);

   /** process trigger patches contained in block **/
   void ProcessTriggerPatches(const AliHLTComponentBlockData* block);

   /** process trigger FastORs contained in block **/
   void ProcessTriggerFastors(const AliHLTComponentBlockData* block, const AliHLTEMCALCaloCells* cells);
   
   /** converts the HLT trigger patch flat structure into an AliEMCALTriggerPatchInfo object */
   void HLTPatch2Patch(const AliHLTCaloTriggerPatchDataStruct& htlpatch, AliEMCALTriggerPatchInfo& patch) const;

   /** converts the HLT trigger FastOR flat structure into an AliEMCALTriggerFastOR object */
   void HLTFastor2Fastor(const AliHLTCaloTriggerDataStruct& htlfastor, AliEMCALTriggerFastOR& fastor) const;

   /** push histograms contained in the list */
   void PushHistograms(TCollection* list);

   /** initialise the geometry */
   void InitialiseGeometry();

   const AliEMCALTriggerBitConfig       *fTriggerBitConfig   ;  ///<  Trigger bit configuration, aliroot-dependent
   Bool_t                                fHistoResetOnPush   ;  ///<  Reset histograms when data is pushed
   TString                               fFilterTrgClass     ;  ///<  Space-separated trigger classes to be taken into consideration
   EBeamType_t                           fBeamType           ;  ///<  Beam type
   int                                   fLocalEventCount    ;  //!<! Event counter
   AliHLTEMCALGeometry                  *fGeometry           ;  //!<! EMCal geometry
   std::vector<TString>                  fFiredTriggerClasses;  //!<! Trigger classes fired in the current event

private:
   /** Pointer to the trigger QA class */
   AliEMCALTriggerQA                          *fTriggerQAPtr;            //! Transient

   /** Copy constructor,  not implemented */
   AliHLTEMCALTriggerQAComponent(const AliHLTEMCALTriggerQAComponent&);

   /** Assignment operator, not implemented */
   AliHLTEMCALTriggerQAComponent& operator=(const AliHLTEMCALTriggerQAComponent);

   ClassDef(AliHLTEMCALTriggerQAComponent, 1);
};

#endif
