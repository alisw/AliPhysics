#ifndef ALISYSINFO_H
#define ALISYSINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
// This is the class which is to be used during the writing of
// simulated raw data (DDL files format).
// It is using the root functionality in order to deal correctly
// with little/big endian issue. By convention the detector raw
// data payload is stored always with little endian (this corresponds
// to the real life situation when the detector data is coming from
// the hardware).
//-------------------------------------------------------------------------

#include <TObject.h>
class TStopwatch;
class TTree;

class AliSysInfo : public TObject {
public:
  AliSysInfo();
  static AliSysInfo * Instance();
  static void AddStamp(const char *stamp);
  void Print(Option_t* option = "") const;
  static TTree * MakeTree(const char *lname);
private:
  AliSysInfo(const AliSysInfo& source);
  AliSysInfo& operator= (const AliSysInfo& rec);

  fstream        *fSysWatch;       // system watch - Memory and CPU usage 
  TStopwatch     *fTimer;          // timer
  static AliSysInfo *   fInstance; //instance pointer
  ClassDef(AliSysInfo,0)
};

#endif
