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
    kDataTypeCustom,      /// Custom format maintained by the user.
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
    
    /**
     * Create a new service with a particular predefined type and name.
     * \param type  The type of the service
     * \param servicename  The name of the service.
     */
    AliHLTDimService(AliHLTDimServiceDataType type, const char* servicename);
    
    /**
     * Create a new service with a particular custom type.
     * \param type  The type of the service as a string.
     *      The format parameter specifies the contents of the structure in the
     *      form T:N[;T:N]*[;T] where T is the item type: (I)nteger, (C)arachter,
     *      (L)ong, (S)hort, (F)loat, (D)ouble, X(tra long) and N is the number
     *      of such items. The type alone at the end means all following items
     *      are of the same type. Example: "I:3;F:2;C" means 3 Integers, 2 Floats
     *      and Characters until the end. The format parameter is used for
     *      communicating between different platforms.
     * \param data  Points to a buffer maintained by the user which stores the
     *      data to publish. This buffer must exist as long as the DIM service
     *      is registered and active.
     * \param size  The size of the data structure pointed to by data.
     * \param servicename  The name of the service.
     */
    AliHLTDimService(const char* type, void* data, int size, const char* servicename);
    
    /**
     * Updates the DIM data point for custom data structures.
     * i.e. This method should be used if the service was created with:
     * AliHLTDimService(const char* type, void* data, int size, const char* servicename)
     */
    void Update();
    
    /**
     * Updates the DIM data point.
     * This method should be used if the service was created with:
     * AliHLTDimService(AliHLTDimServiceDataType type, const char* servicename)
     * \param sp  The new data point to publish via DIM.
     */
    void Update(const AliHLTDimServicePoint_t& sp);
    
    AliHLTDimServiceDataType GetType() const {return fType;}
    const char* GetTypeString() const { return fTypeString.Data(); }
    void* GetLocation() {return fDataBuffer;}
    int GetId() const {return fId;}
    int SetId(int id) {fId=id;return id;}
    void* GetDataBuffer() const { return fDataBuffer; }
    int GetDataSize() const { return fDataSize; }

  private:
          
    // Do not allow copying of this class
    AliHLTDimService(const AliHLTDimService&);
    AliHLTDimService& operator = (const AliHLTDimService&);
          
    AliHLTDimServicePoint_t fData; /// the data point
    AliHLTDimServiceDataType fType; /// type of this service
    TString fTypeString;  /// The string representing the service type.
    void* fDataBuffer;  /// Pointer to the data buffer.
    int fDataSize;  /// The size of the data buffer.
    int fId; /// id of the service
  };

  /** @class AliHLTDimServiceFloat
   * DIM service for a float value
   */
  class AliHLTDimServiceFloat : public AliHLTDimService {
  public:
    AliHLTDimServiceFloat(){}
    ~AliHLTDimServiceFloat(){}

    void Update(float f) {
      AliHLTDimServicePoint_t sp; sp.fVal=f; AliHLTDimService::Update(sp);
    }
  };

  /** @class AliHLTDimServiceInt
   * DIM service for a int value
   */
  class AliHLTDimServiceInt : public AliHLTDimService {
  public:
    AliHLTDimServiceInt(){}
    ~AliHLTDimServiceInt(){}

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
