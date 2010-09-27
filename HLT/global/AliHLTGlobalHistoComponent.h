// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTGLOBALHISTOCOMPONENT_H
#define ALIHLTGLOBALHISTOCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTGlobalHistoComponent.h
/// @author Matthias Richter
/// @date   2010-09-16
/// @brief  A histogramming component for global ESD properties based
///         on the AliHLTTTreeProcessor

#include "AliHLTTTreeProcessor.h"
#include <string>
#include <map>
#include <vector>

/**
 * @class AliHLTGlobalHistoComponent
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b GlobalHisto                                         <br>
 * Library: \b libAliHLTGlobal.so	      				<br>
 * Input Data Types: ::kAliHLTAnyDataType				<br>
 * Output Data Types: none						<br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *      
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Configuration by component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * The component loads no CDB entries.
 *
 * <h2>Performance:</h2>
 * The component does not process any event data.
 *
 * <h2>Memory consumption:</h2>
 * The component does not process any event data.
 *
 * <h2>Output size:</h2>
 * Depending on the mode.
 *
 * @ingroup alihlt_util_components
 */
class AliHLTGlobalHistoComponent : public AliHLTTTreeProcessor
{
 public:
  /// standard constructor
  AliHLTGlobalHistoComponent();
  /// destructor
  virtual ~AliHLTGlobalHistoComponent();

  /// inherited from AliHLTComponent: return id of the component.
  virtual const char* GetComponentID() {return "GlobalHisto";};
  /// inherited from AliHLTComponent: input data types
  virtual void GetInputDataTypes(AliHLTComponentDataTypeList&);

  /// inherited from AliHLTComponent: spawn function, create an instance.
  virtual AliHLTComponent* Spawn() {return new AliHLTGlobalHistoComponent;}

  void FillHistogramDefinitions();

  /// @class AliHLTGlobalHistoVariables
  /// container for the tree branch variables
  class AliHLTGlobalHistoVariables {
  public:
    AliHLTGlobalHistoVariables();
    AliHLTGlobalHistoVariables(int capacity, const char* names);
    ~AliHLTGlobalHistoVariables();

    /// init the arrays
    int Init(int capacity, const char* names);

    /// capacity for every key
    int Capacity() const {return fCapacity;}
    /// number of variables
    int Variables() const {return fArrays.size();}

    /// fill variable at index
    int Fill(unsigned index, float value);
    /// fill variable at key
    int Fill(const char* key, float value);
    /// get array at key
    float* GetArray(const char* key);
    /// get the key of an array
    const char* GetKey(int index) const;

    /// reset and cleanup arrays
    int Reset();

    /// reset the fill counts
    int ResetCount();

  private:
    int FindKey(const char* key) const;

    /// capacity of all arrays
    int fCapacity; //!
    /// pointers of arrays
    vector<float*> fArrays; //!
    /// fill count for arrays
    vector<int> fCount; //!
    /// map of keys
    map<string,int> fKeys; //!
  };

 protected:
  /// inherited from AliHLTTTreeProcessor: create the tree instance and all branches
  TTree* CreateTree(int argc, const char** argv);
  /// inherited from AliHLTTTreeProcessor: process input blocks and fill tree
  int FillTree(TTree* pTree, const AliHLTComponentEventData& evtData, 
                       AliHLTComponentTriggerData& trigData );
  /// dtOrigin for PushBack.
  AliHLTComponentDataType GetOriginDataType() const;
  /// spec for PushBack
  AliHLTUInt32_t GetDataSpec() const {return 0;}

  int ResetVariables();

private:
  /// copy constructor prohibited
  AliHLTGlobalHistoComponent(const AliHLTGlobalHistoComponent&);
  /// assignment operator prohibited
  AliHLTGlobalHistoComponent& operator=(const AliHLTGlobalHistoComponent&);

  /// the event number, tree filling variable
  int fEvent; //!
  /// track count, tree filling variable
  int fNofTracks; //!
  /// filling arrays for track parameters
  AliHLTGlobalHistoComponent::AliHLTGlobalHistoVariables fTrackVariables; //!

  ClassDef(AliHLTGlobalHistoComponent, 0) // HLT Global Histogram component
};
#endif
