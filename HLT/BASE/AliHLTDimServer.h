//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTDIMSERVER_H
#define ALIHLTDIMSERVER_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTDimServer.h
//  @author Matthias Richter
//  @date   20010-03-10
//  @brief  HLT DIM server implementation and dynamic access
//          to DIM library

#include "AliHLTLogging.h"
#include "TNamed.h"
#include "TObjArray.h"

class TThread;

/**
 * @class AliHLTDimServer
 * Implementation of a DIM server for HLT and the dynamic access to the
 * DIM  library.
 */
class AliHLTDimServer : public TNamed {
public:
  AliHLTDimServer();
  AliHLTDimServer(const char* servername);
  ~AliHLTDimServer();

  /// Data type identifiers for services.
  enum AliHLTDimServiceDataType{
    kDataTypeUnknown = 0, /// initializer
    kDataTypeInt,         /// Integer type
    kDataTypeFloat,       /// Float type
    kDataTypeString,      /// String type
  };

  /// The service data field.
  struct AliHLTDimServicePoint_t {
    union {
      int iVal;     /// integer value
      float fVal;   /// float value
      void* strVal; /// string value, casted to string* before use
    };
  };

  /** @class AliHLTDimService
   * Base class for DIM services
   */
  class AliHLTDimService : public TNamed {
  public:
    AliHLTDimService();
    AliHLTDimService(AliHLTDimServiceDataType type, const char* servicename);
    
    void Update(const AliHLTDimServicePoint_t& sp);
    AliHLTDimServiceDataType GetType() const {return fType;}
    void* GetLocation() {return &fData.iVal;}
    int GetId() const {return fId;}
    int SetId(int id) {fId=id;return id;}

  private:
    AliHLTDimServicePoint_t fData; /// the data point
    AliHLTDimServiceDataType fType; /// type of this service
    int fId; /// id of the service
  };

  /** @class AliHLTDimServiceFloat
   * DIM service for a float value
   */
  class AliHLTDimServiceFloat : public AliHLTDimService {
  public:
    AliHLTDimServiceFloat();
    ~AliHLTDimServiceFloat();

    void Update(float f) {
      AliHLTDimServicePoint_t sp; sp.fVal=f; AliHLTDimService::Update(sp);
    }
  };

  /** @class AliHLTDimServiceInt
   * DIM service for a int value
   */
  class AliHLTDimServiceInt : public AliHLTDimService {
  public:
    AliHLTDimServiceInt();
    ~AliHLTDimServiceInt();

    void Update(int i) {
      AliHLTDimServicePoint_t sp; sp.iVal=i; AliHLTDimService::Update(sp);
    }
  };

  /**
   * Register a service.
   * @param pService    the service to be registered
   */
  int RegisterService(AliHLTDimService* pService);

  /**
   * Create a service.
   * @param type        type of the channel, see @ref ceServiceDataType
   * @param name        unique name of the service
   * @return dim service object, needs to be cleaned by the caller
   */
  AliHLTDimService* CreateService(AliHLTDimServiceDataType type, const char* name);

  /**
   * Create a group of services.
   * The names are built from the basename and the number of services.
   * @param type        type of the channel
   * @param basename    base name of the services, the name might contain a '%d' sequence which is then
   *                    replaced by the number, number is appended if no '%d' provided
   * @param count       number of services in this group, passed to the <i>update</i> and <i>set</i> function as parameter major
   * @return            TObjArray of AliHLTDimService objects, the array needs to be cleaned by the caller
   */
  TObjArray* CreateServiceGroup(AliHLTDimServiceDataType type, const char* basename, int count);

  /// Update all services via the Dim channel
  int UpdateServices();

  /// init the server
  /// load the dim library and function pointers
  /// init dim (DNS and server name)
  int Init(const char* dimNameServer);

  /// Reset
  int Reset();

  /// start the server
  int Start();

  /// stop the server
  int Stop();

protected:
  enum AliHLTDimServerState_t {
    // server is not started
    kStateOff = 0,
    // starting, will be changed by the server thread to kStateRunning
    kStateStarting,
    // server running
    kStateRunning,
    // set by the main thread and changed by the server thread before it terminates
    kStateStopping,
    // error
    kStateError
  };

  int SetState(int state) {fState=state; return fState;}

  int GetState() const {return fState;}

  typedef void (*fctVoid)();
  typedef int (*fctDisServiceCallback)( const char*);
  typedef int (*fctDisAddService)     ( const char* service, 
					const char* type, 
					void* buffer, 
					int size, 
					fctDisServiceCallback cb, 
					long int tag);
  typedef int (*fctDisRemoveService)  ( unsigned int id);
  typedef int (*fctDisUpdateService)  ( unsigned int id);
  typedef int (*fctDisCharArg)        ( const char*);
  typedef int (*fctDisNoArg)          ( );

  /** 
   * @class AliHLTDimInterface
   * Interface to the dim library
   */
  class AliHLTDimInterface : public AliHLTLogging {
  public:
    AliHLTDimInterface();
    ~AliHLTDimInterface();

    /// load the dim library and function pointers
    int Init();

    int DisAddService(const char* service, const char* type, void* buffer, 
		      int size, fctDisServiceCallback cb, long int tag) const {
      if (fpDisAddService) return (*fpDisAddService)(service, type, buffer, size, cb, tag);
      return -ENODEV;
    }

    int DisAddService(const char* service, const char* type, void* buffer, int size) const {
      if (fpDisAddService) return (*fpDisAddService)(service, type, buffer, size, NULL, 0);
      return -ENODEV;
    }

    int DisRemoveService(unsigned int id) const {
      if (fpDisRemoveService) return (*fpDisRemoveService)(id);
      return -ENODEV;
    }

    int DisUpdateService(unsigned int id) const {
      if (fpDisUpdateService) return (*fpDisUpdateService)(id);
      return -ENODEV;
    }

    int DisStartServing(const char *server) const {
      if (fpDisStartServing) return (*fpDisStartServing)(server);
      return -ENODEV;
    }

    int DisStopServing() const {
      if (fpDisStopServing) return (*fpDisStopServing)();
      return -ENODEV;
    }

    int DisSetDnsNode(const char *server) const {
      if (fpDisSetDnsNode) return (*fpDisSetDnsNode)(server);
      return -ENODEV;
    }

  private:
    fctVoid FindSymbol(const char* library, const char* symbol) const;

    fctDisAddService      fpDisAddService;      //! transient
    fctDisRemoveService   fpDisRemoveService;   //! transient
    fctDisUpdateService   fpDisUpdateService;   //! transient
    fctDisCharArg         fpDisStartServing;    //! transient
    fctDisNoArg           fpDisStopServing;     //! transient
    fctDisCharArg         fpDisSetDnsNode;      //! transient
    static const char*    fgkDimLibraryName  ;       //!
    static const char*    fgkDisAddServiceSymbol;    //!
    static const char*    fgkDisRemoveServiceSymbol;    //!
    static const char*    fgkDisUpdateServiceSymbol; //!
    static const char*    fgkDisStartServingSymbol;  //!
    static const char*    fgkDisStopServingSymbol;  //!
    static const char*    fgkDisSetDnsNodeSymbol;  //!
  };

  static AliHLTDimInterface* Interface();

private:
  /// copy constructor not permitted
  AliHLTDimServer(const AliHLTDimServer&);
  /// assignment operator not permitted
  AliHLTDimServer& operator=(const AliHLTDimServer&);

  /// entry point for the thread, param is pointer to object
  static void* ServerLoop(void* param);

  /// the server loop
  void* ServerLoop();

  TObjArray fServices; //! list of services
  int fState; //! state of the server
  TThread* fpServerThread; //! thread
  int fUpdatePeriod; //! update period for the DIM loop in ms

  static AliHLTDimInterface* fgpInterface; //! the dim interface

  ClassDef(AliHLTDimServer, 0)
};
#endif
