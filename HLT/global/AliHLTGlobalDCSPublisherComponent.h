//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTGLOBALDCSPUBLISHERCOMPONENT_H
#define ALIHLTGLOBALDCSPUBLISHERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTGlobalDCSPublisherComponent.h
    @author Matthias Richter
    @date   20010-03-10
    @brief  DIM publisher component for global HLT data
*/
#include "AliHLTDataSink.h"
#include "AliHLTDimServer.h"
#include <functional>

/**
 * @class AliHLTGlobalDCSPublisherComponent
 * DIM Publisher component for global HLT data.
 * It implements a DIM server which publishes global HLT data through the
 * following services:
 * - Vertex_X
 * - Vertex_Y
 * - Vertex_Z
 * - RmsVertex_X
 * - RmsVertex_Y
 * - RmsVertex_Z
 * 
 */
class AliHLTGlobalDCSPublisherComponent : public AliHLTDataSink {
public:
  AliHLTGlobalDCSPublisherComponent();
  ~AliHLTGlobalDCSPublisherComponent();

  /// service enumerators
  enum {
    kVertexX,
    kVertexY,
    kVertexZ,
    kDimensions,
    kVertexRmsX=kDimensions,
    kVertexRmsY,
    kVertexRmsZ,
    kLastService
  };

  virtual const char* GetComponentID();
  virtual void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  virtual AliHLTComponent* Spawn();

 protected:
  /// component initialization
  int DoInit( int argc, const char** argv );

  /// component cleanup
  int DoDeinit();

  /// Data processing method for the component.
  virtual int DumpEvent( const AliHLTComponentEventData& evtData,
			 AliHLTComponentTriggerData& trigData );

  using AliHLTDataSink::DumpEvent;

  int ScanConfigurationArgument(int argc, const char** argv);
  int Publish(bool reset=false);

private:
  /// copy constructor not permitted
  AliHLTGlobalDCSPublisherComponent(const AliHLTGlobalDCSPublisherComponent&);
  /// assignment operator not permitted
  AliHLTGlobalDCSPublisherComponent& operator=(const AliHLTGlobalDCSPublisherComponent&);

  /// binary function for the calculation of rms using std::accumulate
  template <class T> class MeanSqtOp {
  public:
    MeanSqtOp() : fMean(0) {}
    MeanSqtOp(T mean) : fMean(mean) {}
    T operator() (const T& x, const T& y) const {
      T val=x+(y-fMean)*(y-fMean);
      return val;
    }
  private:
    T fMean;
  };

  void Reset(vector<float>& sample, int size) const;
  template <class T> T Mean(const vector<T>& sample, int count) const;
  template <class T> T Rms(const vector<T>& sample, T mean, int count) const;

  /// dim server name
  string fServerName; //!
  /// dim DNS server name
  string fDimdns; //!
  /// dim server instance
  AliHLTDimServer* fpServer; //!

  /// position in the event buffer
  int fPosition; //!
  /// size of the event buffer
  int fEventBufferSize; //!
  /// last update time
  unsigned int fLastUpdate; //!
  /// update period in seconds
  unsigned int fUpdatePeriod; //!

  /// event buffers
  vector<float> fEventBuffers[kDimensions]; //!
  /// array of dim services
  AliHLTDimServer::AliHLTDimService* fDimServices[kLastService]; //!
  /// dim service for the number of samples for
  AliHLTDimServer::AliHLTDimService* fDimServiceEventCount; //!

  /// names of services
  static const char* fgkServiceNames[kLastService]; //!

  ClassDef(AliHLTGlobalDCSPublisherComponent, 0)
};
#endif
