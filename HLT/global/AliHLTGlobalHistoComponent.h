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
#include "TObjString.h"
#include "TObjArray.h"
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
 
  
  /// interface function, see AliHLTComponent for description
  AliHLTComponentDataType GetOutputDataType();
  
  
  /// inherited from AliHLTComponent: spawn function, create an instance.
  virtual AliHLTComponent* Spawn() {return new AliHLTGlobalHistoComponent;}

  void FillHistogramDefinitions();

  /// @class AliHLTGlobalHistoVariables
  /// container for the tree branch variables
  template <typename T>
  class AliHLTGlobalHistoVariables {
  public:
    AliHLTGlobalHistoVariables()  : fCapacity(0), fArrays(), fCount(), fKeys() {}
    ~AliHLTGlobalHistoVariables() {Reset();}

    /// init the arrays
    int Init(int capacity, const char* names)
    {
      /// init the arrays
      int iResult=0;
      TString initializer(names);
      TObjArray* pTokens=initializer.Tokenize(" ");
      fCapacity=capacity;
      if (pTokens) {
	int entries=pTokens->GetEntriesFast();
	fArrays.resize(entries);
	fCount.resize(entries);
	for (int i=0; i<entries; i++) {
	  fKeys[pTokens->At(i)->GetName()]=i;
	  fArrays[i]=new T[fCapacity];
	}
	delete pTokens;
      }
      if (fArrays.size()!=fCount.size() ||
	  fArrays.size()!=fKeys.size()) {
	return -EFAULT;
      }

      ResetCount();
      return iResult;
    }

    /// capacity for every key
    int Capacity() const {return fCapacity;}
    /// number of variables
    int Variables() const {return fArrays.size();}

    /// fill variable at index
    int Fill(unsigned index, T value)
    {
      if (index>=fArrays.size() || index>=fCount.size()) return -ENOENT;
      if (fCount[index]>=fCapacity) return -ENOSPC;

      (fArrays[index])[fCount[index]++]=value;
      return fCount[index];
    }

    /// fill variable at key
    int Fill(const char* key, T value)
    {
      int index=FindKey(key);
      if (index<0) return -ENOENT;
      return Fill(index, value);
    }

    /// get array at key
    T* GetArray(const char* key)
    {
      int index=FindKey(key); if (index<0) return NULL;
      if ((unsigned)index>=fArrays.size()) return NULL;
      return fArrays[index];
    }

    /// get the key of an array
    const char* GetKey(int index) const
    {
      for (map<string, int>::const_iterator element=fKeys.begin(); element!=fKeys.end(); element++) {
	if (element->second==index) return element->first.c_str();
      }
      return NULL;
    }

    /// reset and cleanup arrays
    int Reset()
    {
      for (unsigned i=0; i<fArrays.size(); i++) {delete fArrays[i];}
      fArrays.clear(); fCount.clear(); fKeys.clear();
      return 0;
    }

    /// reset the fill counts
    int ResetCount()
    {
      for (vector<int>::iterator element=fCount.begin(); element!=fCount.end(); element++) *element=0;
      return 0;
    }

    char GetType() const {AliHLTGlobalHistoVariablesType type(T&); return type.GetType();}

  private:
    AliHLTGlobalHistoVariables(const AliHLTGlobalHistoVariables& src);
    AliHLTGlobalHistoVariables& operator=(const AliHLTGlobalHistoVariables& src);

    int FindKey(const char* key) const
    {
      map<string, int>::const_iterator element=fKeys.find(key);
      if (element==fKeys.end()) return -ENOENT;
      return element->second;
    }

    /// internal helper class to get the type of the template
    class AliHLTGlobalHistoVariablesType {
    public:
      AliHLTGlobalHistoVariablesType(float&) : fType('f') {}
      AliHLTGlobalHistoVariablesType(int&)   : fType('i') {}
      char GetType() const {return fType;}
    private:
      char fType; //!
    };

    /// capacity of all arrays
    int fCapacity; //!
    /// pointers of arrays
    vector<T*> fArrays; //!
    /// fill count for arrays
    vector<int> fCount; //!
    /// map of keys
    map<string,int> fKeys; //!
  };

 protected:
  /// inherited from AliHLTTTreeProcessor: create the tree instance and all branches
  TTree* CreateTree(int argc, const char** argv);
  /// inherited from AliHLTTTreeProcessor: process input blocks and fill tree
  int FillTree(TTree* pTree, const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );
  /// dtOrigin for PushBack.
  AliHLTComponentDataType GetOriginDataType() const;
  /// clean up variables
  int ResetVariables();
  /// inherited from AliHLTComponent, scan argument
  int ScanConfigurationArgument(int argc, const char** argv);
  /// function for online reconfiguration
  int Reconfigure(const char* cdbEntry, const char* chainId);

private:
  /// copy constructor prohibited
  AliHLTGlobalHistoComponent(const AliHLTGlobalHistoComponent&);
  /// assignment operator prohibited
  AliHLTGlobalHistoComponent& operator=(const AliHLTGlobalHistoComponent&);
    
  /// the event number, tree filling variable
  int fEvent; //!
  /// track count, tree filling variable
  int fNofTracks; //!
  /// V0 count, tree filling variable
  int fNofV0s; //!
  /// contributors count, tree filling variable
  int fNofContributors; //!
  /// x coordinate of vertex
  float fVertexX; //!
  /// y coordinate of vertex
  float fVertexY; //!
  /// z coordinate of vertex
  float fVertexZ; //!
  /// vertex status, found or not
  bool fVertexStatus; //!
  /// maximum track multiplicity
  int fMaxTrackCount; //!
  /// maximum number of V0 entries
  int fMaxV0Count; //!
  /// activate event properties branch
  bool fFillV0; //!
 
  /// filling arrays for track parameters
  AliHLTGlobalHistoVariables<float> fTrackVariables; //!
  /// filling array for the track status
  AliHLTGlobalHistoVariables<int> fTrackVariablesInt; //!
  /// filling arrays for V0 parameters
  AliHLTGlobalHistoVariables<float> fV0Variables; //!
   
  ClassDef(AliHLTGlobalHistoComponent, 0) // HLT Global Histogram component
};
#endif
