#ifndef ALIRAWDATAERRORLOG_H
#define ALIRAWDATAERRORLOG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
// class AliRawDataErrorLog                                        //
// This is a class for logging raw-data related errors.            //
// It is used to record and retrieve of the errors                 //
// during the reading and reconstruction of raw-data and ESD       //
// analysis.                                                       //
// Further description of the methods and functionality are given  //
// inline.                                                         //
//                                                                 //
// cvetan.cheshkov@cern.ch                                         //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include <TNamed.h>

class AliRawDataErrorLog: public TNamed {

 public:

  enum ERawDataErrorType {
    kNone = 0, 
    kMinor = 1, 
    kMajor = 2, 
    kFatal = 3, 
  };

  AliRawDataErrorLog();
  AliRawDataErrorLog(Int_t eventNumber, Int_t ddlId,
		     ERawDataErrorType errorType,
		     const char *message = NULL);
  AliRawDataErrorLog(const AliRawDataErrorLog & source);
  AliRawDataErrorLog & operator=(const AliRawDataErrorLog & source);
  virtual ~AliRawDataErrorLog() {};

  Int_t             GetEventNumber() const { return fEventNumber; }
  Int_t             GetDdlID()       const { return fDdlID; }
  ERawDataErrorType GetErrorType()   const { return fErrorType; }
  const char *      GetMessage()     const { return fName.Data(); }

  Bool_t            IsSortable() const {return kTRUE;}
  Int_t             Compare(const TObject* obj) const;

 private:

  Int_t             fEventNumber; // Event number as it appears in the input raw-data file
  Int_t             fDdlID;       // ID of the DLL in which the error occured
  ERawDataErrorType fErrorType;   // Type of the raw data error

  ClassDef(AliRawDataErrorLog, 1)
};

#endif
