#ifndef ALIHLTTRDESDWRITERCOMPONENT_H
#define ALIHLTTRDESDWRITERCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTRDEsdWriterComponent.h
    @author Mateusz Ploskon
    @date   
    @brief  Writer component to store tracks of the HLT TRD 

                                                                          */
#include "AliHLTRootFileWriterComponent.h"

// forward declarations
class TTree;
class AliESDEvent;

/**
 * @class AliHLTTRDEsdWriterComponent
 * @see AliHLTFileWriter and AliHLTRootFileWriterComponent for more parameters
 */
class AliHLTTRDEsdWriterComponent : public AliHLTRootFileWriterComponent
{
 public:
  /** standard constructor */
  AliHLTTRDEsdWriterComponent();
  /** destructor */
  ~AliHLTTRDEsdWriterComponent();

  /**
   * The id of the component.
   * @return component id (string)
   */
  const char* GetComponentID() {return "TRDEsdWriter";};

  void GetInputDataTypes( vector<AliHLTComponent_DataType>& list);

  /**
   * Spawn function.
   * @return new class instance
   */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn() {return new AliHLTTRDEsdWriterComponent;}

 protected:
  /**
   * Data processing method for the component.
   * The function can be overloaded by specific ROOT file writer
   * components.
   * @param evtData       event data structure
   * @param blocks        input data block descriptors
   * @param trigData	  trigger data structure
   */
  virtual int DumpEvent( const AliHLTComponentEventData& evtData,
			 const AliHLTComponentBlockData* blocks, 
			 AliHLTComponentTriggerData& trigData );
			 
  using AliHLTRootFileWriterComponent::DumpEvent;

  /**
   * Scan one argument and adjacent parameters.
   * @param argc           size of the argument array
   * @param argv           agument array for component initialization
   * @return number of processed members of the argv <br>
   *         -EINVAL unknown argument <br>
   *         -EPROTO parameter for argument missing <br>
   */
  int ScanArgument(int argc, const char** argv);
  
  

 private:
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTTRDEsdWriterComponent(const AliHLTTRDEsdWriterComponent&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTTRDEsdWriterComponent& operator=(const AliHLTTRDEsdWriterComponent&);

  /**
   * Init the writer.
   * The DoInit function is not available for this child class. InitWriter is the
   * corresponding function for classes derived from AliHLTFileWriter.
   */
  int InitWriter();

  /**
   * Init the writer.
   * The DoDeinit function is not available for this child class. CloseWriter is the
   * corresponding function for classes derived from AliHLTFileWriter.
   */
  int CloseWriter();

  /** the ESD tree */
  TTree* fTree; //! transient value
  
  /** output percentage descared*/
  Int_t fOutputPercentage; 

  /** the ESD */
  AliESDEvent* fESD; //! transient value

  int DoEvent(	const AliHLTComponent_EventData& evtData,
						const AliHLTComponent_BlockData* blocks,
						AliHLTComponent_TriggerData& /*trigData*/,
						AliHLTUInt8_t* /*outputPtr*/,
						AliHLTUInt32_t& /*size*/,
						vector<AliHLTComponent_BlockData>& /*outputBlocks*/);
  						
  ClassDef(AliHLTTRDEsdWriterComponent, 1)

};
#endif
