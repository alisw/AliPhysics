#ifndef ALIHLTPGLOBALVERTEXERHISTOCOMPONENT_H
#define ALIHLTPGLOBALVERTEXERHISTOCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTGlobalVertexerComponent.h
    @author Timur Pocheptsov
    @brief  Component for monitor primary vertex
*/

#include <TH2F.h>

#include "AliHLTProcessor.h"


class AliHLTGlobalVertexerHistoComponent : public AliHLTProcessor 
{
public:
  /** default constructor */
  AliHLTGlobalVertexerHistoComponent();
  
  //No need to declare and define dtor, copy-ctor, copy-assignment operator.
  //Compiler generated versions are ok (as soon as class has no explicit 
  //resource management). And if TH2F is exception-safe, my class is
  //exception-safe too.

  // Overriders for AliHLTComponent's interface.
  // These functions are required for the registration process.
  /** interface function, see AliHLTComponent for description. why not const? */
  const char* GetComponentID();
  /** interface function, see AliHLTComponent for description. why not const? */
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  /** interface function, see AliHLTComponent for description. why not const? */
  AliHLTComponentDataType GetOutputDataType();
  /** interface function, see AliHLTComponent for description */
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);

   /** interface function, see AliHLTComponent for description */
 virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();

private:  
  /** interface function, see AliHLTComponent for description */
  int DoInit(int argc, const char** argv);
  /** interface function, see AliHLTComponent for description */
  int DoDeinit();

  //AliHLTProcess has overriden DoEvent from AliHLTComponent and 
  //added new virtual function with the same name and different
  //signature. If I override only one of them, the second will be hidden
  //and good compiler can issue a warning, for example, Comeau can 
  //(my gcc does not issue any warning). So, there are 2 ways:
  //override both or declare the second with using-declaration.

  /** interface function, see AliHLTComponent for description */
  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
              AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
              AliHLTUInt32_t& size, AliHLTComponentBlockDataList& outputBlocks);

  using AliHLTProcessor::DoEvent;

  /** interface function, see AliHLTComponent for description */
  int Reconfigure(const char* cdbEntry, const char* chainId);

  TH2F fPrimaryXY;  //X and Y coords.
  TH2F fPrimaryZX;  //Z and X coords.
  TH2F fPrimaryZY;  //Z and Y coords.

  TH2F fSPDVertexXY;// ITS SPD vertex
  TH1F fSPDVertexX;// ITS SPD vertex
  TH1F fSPDVertexY;// ITS SPD vertex
  TH1F fSPDVertexZ;// ITS SPD vertex

  //Aux. function.
  unsigned long GetOutputDataSize()const;
};

#endif
