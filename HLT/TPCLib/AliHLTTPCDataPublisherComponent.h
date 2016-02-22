//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCDATAPUBLISHERCOMPONENT_H
#define ALIHLTTPCDATAPUBLISHERCOMPONENT_H
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCDataPublisherComponent.h
/// @author Matthias Richter
/// @date   2011-11-18
/// @brief  Specific publisher for TPC raw data from the AliRawReader
///         

#include "AliHLTRawReaderPublisherComponent.h"
#include "AliHLTTPCRawCluster.h"
#include <map>

class AliHLTTPCClusterMCLabel;
class AliHLTTPCClusterMCData;
class AliHLTTPCDataCompressionDecoder;

/**
 * @class AliHLTTPCDataPublisherComponent
 * This component uses the functionality of AliHLTRawReaderPublisherComponent
 * and overloads IsSelected and GetSpecificationFromEquipmentId. Blocks are
 * only generated if the corresponding partition is missing in HLTOUT.
 *
 * It is used in an emulation chain which produces all compressed cluster
 * blocks which are missing in HLTOUT. If TPC reconstruction requires HLT
 * clusters, the emulator is automatically executed and the compressed
 * data produced if raw data is available.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCRawReaderPublisher      <br>
 * Library: \b libAliHLTTPC.so     <br>
 * Input Data Types:  <br>
 * Output Data Types: <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->

 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 *
 * <h2>Performance:</h2>
 *
 * <h2>Memory consumption:</h2>
 *
 * <h2>Output size:</h2>
 *
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCDataPublisherComponent : public AliHLTRawReaderPublisherComponent {
public:
  /// standard constructor
  AliHLTTPCDataPublisherComponent();
  /// destructor
  ~AliHLTTPCDataPublisherComponent();

  enum {
    kPublisherModeDefault = 0,
    kRegisterClusterBlocks= 0x1, // only register data blocks
    kPublishClustersAll   = 0x2, // unpack data blocks
    kPublishRawAll        = 0x4, // publish all raw data
    kPublishRawFiltered   = 0x8, // publish raw data filtered by existence of clusters
    kLastContainerMode
  };

   
  /// get the output data types of the component.
   
  AliHLTComponentDataType GetOutputDataType();  
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);			   

 /// set mode
  void SetMode(int mode) {fMode=mode;}

  bool CheckMode(int flag) const {return (fMode&flag)==flag;}

  /// inherited from AliHLTComponent: id of the component
  virtual const char* GetComponentID();

  /// inherited from AliHLTComponent: spawn function.
  virtual AliHLTComponent* Spawn();

  /**
   * @class AliRawClusterContainer
   * Container of AliHLTTPCRawCluster, The class implements the interface to be
   * used in the decoding of compressed TPC data.
   * Data is decoded into an external buffer.
   */
  class AliRawClusterContainer : public AliHLTLogging {
  public:
    AliRawClusterContainer();
    virtual ~AliRawClusterContainer();

    /// set/reset the external target buffer
    int SetTargetBuffer(AliHLTUInt8_t* pBuffer, int size);

    /// merge track model clusters into partition cluster blocks 
    int Sort();

    /// fill block descriptors of extracted partition cluster blocks to target list
    int CopyBlockDescriptors(AliHLTComponentBlockDataList& target) const;
    /// get reference to block descriptor list
    const AliHLTComponentBlockDataList& GetBlockDescriptors() const {
      return fDescriptors;
    }

    struct AliClusterIdBlock {
      AliClusterIdBlock() : fIds(NULL), fSize(0) {}
      AliHLTUInt32_t* fIds; //!
      AliHLTUInt32_t  fSize; //!
    };

    class iterator {
    public:
      iterator() : fClusterNo(-1), fCluster(NULL), fClusterId(kAliHLTVoidDataSpec), fContainer(NULL) {}
      iterator(AliRawClusterContainer* pContainer) : fClusterNo(-1), fCluster(NULL), fClusterId(kAliHLTVoidDataSpec), fContainer(pContainer) {}
      iterator(const iterator& other) : fClusterNo(other.fClusterNo), fCluster(other.fCluster), fClusterId(other.fClusterId), fContainer(other.fContainer) {}
      iterator& operator=(const iterator& other) {
	if (this==&other) return *this;
	fClusterNo=other.fClusterNo; fCluster=other.fCluster, fClusterId=other.fClusterId; fContainer=other.fContainer; return *this;
      }
      ~iterator() {fCluster=NULL; fContainer=NULL;}

      void SetPadRow(int row)          {if (fCluster) fCluster->SetPadRow(row);}
      void SetPad(float pad) 	       {if (fCluster) fCluster->SetPad(pad);}
      void SetTime(float time) 	       {if (fCluster) fCluster->SetTime(time);}
      void SetSigmaY2(float sigmaY2)   {if (fCluster) fCluster->SetSigmaPad2(sigmaY2);}
      void SetSigmaZ2(float sigmaZ2)   {if (fCluster) fCluster->SetSigmaTime2(sigmaZ2);}
      void SetCharge(unsigned charge)  {if (fCluster) fCluster->SetCharge(charge);}
      void SetQMax(unsigned qmax)      {if (fCluster) fCluster->SetQMax(qmax);}
      void SetMC(const AliHLTTPCClusterMCLabel* pMC) {
	if (!fCluster || !pMC) return;
      }

      // switch to next cluster
      iterator& Next(int slice, int partition);

    private:
      int fClusterNo; //! cluster no in the current block
      AliHLTTPCRawCluster* fCluster; //! pointer to current cluster
      AliHLTUInt32_t fClusterId; //! id of the cluster, from optional cluster id blocks
      AliRawClusterContainer* fContainer; // instance of container
    };

    /// legacy, to be removed later
    iterator& BeginRemainingClusterBlock(int count, AliHLTUInt32_t specification) {
      return BeginPartitionClusterBlock(count, specification);
    }
    /// iterator of partition clusters block of specification
    iterator& BeginPartitionClusterBlock(int count, AliHLTUInt32_t specification);
    /// iterator of track model clusters
    iterator& BeginTrackModelClusterBlock(int count);
    /// base method to start cluster iterator
    iterator& ClusterIterator(int count, AliHLTComponentDataType dt, AliHLTUInt32_t specification, AliHLTTPCRawClusterData* &pData);

    /// get block count, i.e. number of calls to create an iterator
    int GetBlockCount() const {return fBlockCount;}
    /// get number of decoded clusters
    /// Note: only if there is enough space in the target buffer the clusters
    //  will be properly written
    int GetClusterCount() const {return fTotalClusterCount;}
    /// get the state of the cluster decoding
    int GetState() const {return fState;}

    /// internal cleanup
    virtual void  Clear(Option_t * option="");
    /// print info
    virtual void Print(Option_t *option=NULL) const;

  protected:
    AliHLTTPCRawCluster* NextCluster(int slice, int partition);

  private:
    AliRawClusterContainer(const AliRawClusterContainer&);
    AliRawClusterContainer& operator=(const AliRawClusterContainer&);

    int fBlockCount; //! number of data blocks with clusters
    int fTotalClusterCount; //! total number of decoded clusters
    int fBlockClusterCount; //! number of decoded clusters in current block
    AliHLTUInt8_t* fpBuffer; //! target buffer for decoded data
    int fBufferSize; //! size of target buffer
    AliHLTComponentBlockDataList fDescriptors; //! list of block descriptors
    AliHLTTPCRawClusterData* fCurrentBlock; // current cluster block
    AliHLTTPCRawClusterData* fTrackModelClusters; //! track model cluster block
    vector<AliHLTUInt32_t>   fTrackModelClusterMap; //! slice-partition map for track model clusters
    iterator fIterator; //! iterator for filling of data
    int fState; //! state
  };

protected:
  /// inherited from AliHLTDataSource: get one event
  int GetEvent( const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& outputBlocks );

  /// inherited from AliHLTComponent: initialize
  int DoInit( int argc, const char** argv );

  /// inherited from AliHLTComponent: cleanup
  int DoDeinit();

  /// inherited from AliHLTComponent: argument scan
  int ScanConfigurationArgument(int argc, const char** argv);

  /// read cluster from HLTOUT
  int ReadClusterFromHLTOUT(AliRawClusterContainer* pContainer);

  /// inherited from AliHLTRawReaderPublisherComponent: get specification
  virtual int GetSpecificationFromEquipmentId(int id, AliHLTUInt32_t &specification) const;

  /// inherited from AliHLTRawReaderPublisherComponent: check if a block is selected or not
  virtual bool IsSelected(int equipmentId) const;

private:
  AliHLTTPCDataPublisherComponent(const AliHLTTPCDataPublisherComponent&);
  AliHLTTPCDataPublisherComponent& operator=(const AliHLTTPCDataPublisherComponent&);

  int fMode; //! operation mode
  bool* fArraySelected; //! transient
  AliRawClusterContainer* fClusters; // target for decoded clusters
  AliHLTTPCDataCompressionDecoder* fpDecoder; // decoder for compressed cluster blocks

  ClassDef(AliHLTTPCDataPublisherComponent, 0)
};

#endif //ALIHLTTPCDATAPUBLISHERCOMPONENT_H
