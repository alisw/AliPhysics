#ifndef ALITRDSAXHANDLER_H
#define ALITRDSAXHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * * See cxx source for full Copyright notice */

/* $Id: AliTRDSaxHandler.h 26327 2008-06-02 15:36:18Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The SAX XML file handler used in the preprocessor                     //
//                                                                        //
//  Author:                                                               //
//    Frederick Kramer (kramer@ikf.uni-frankfurt.de)                      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <TObject.h>
#include "Cal/AliTRDCalDCSGTUCtpOpc.h"
#include "Cal/AliTRDCalDCSGTUBoardInfo.h"
#include "Cal/AliTRDCalDCSGTUSegment.h"
#include "Cal/AliTRDCalDCSGTUTmu.h" 

class TObjArray;
class AliTRDCalDCSv2;
class AliTRDCalDCSFEEv2;
class AliTRDCalDCSPTR;
class AliTRDCalDCSGTU;


class AliTRDSaxHandler : public TObject {

public:
  enum { 
    kInsideFEE = 1, 
    kInsidePTR = 2,
    kInsideGTU = 3 
  }; // System level
  enum { 
    kInsideTgu = -1,
    kInsideNone = -2,
    kInsideSegment = -3,
    kInsideGainTable = -4 
  }; // The level under system (1)
  enum { 
    kInsideTmu = 10,
    kInsideSmu = 11 
  }; // The level under that   (2)

  AliTRDSaxHandler();
  AliTRDSaxHandler(const AliTRDSaxHandler &sh);
  virtual ~AliTRDSaxHandler();
  AliTRDSaxHandler &operator=(const AliTRDSaxHandler &sh);

  TObjArray*         GetDCSFEEDataArray() const { return fFEEArr;        }
  TObjArray*         GetDCSPTRDataArray() const { return fPTRArr;        }
  AliTRDCalDCSv2*    GetCalDCSObj(); // to be called by the preprocessor
  void               ParseConfigName(TString cfgname) const;


  Int_t              GetHandlerStatus() const { return fHandlerStatus; }

  // functions for all possible events
  void               OnStartDocument() const;
  void               OnEndDocument() const;
  void               OnStartElement(const char *name, const TList *attributes);
  void               OnEndElement(const char *name);
  void               OnCharacters(const char *name);
  void               OnComment(const char *name) const;
  void               OnWarning(const char *name);
  void               OnError(const char *name);
  void               OnFatalError(const char *name);
  void               OnCdataBlock(const char *name, Int_t len) const;

 private:

  bool               CompareString(TString str, const char *str2); 

  Int_t              fHandlerStatus;      // 0: everything OK, >0: error
  Int_t              fNDCSPTR;            // number of current PTR unit (to be abandonned soon)
  Int_t              fNDCSGTU;            // number of current GTU unit (to be abandonned soon)
  TObjArray*         fFEEArr;             // array of AliTRDCalDCSFEEv2 objects
  TObjArray*         fPTRArr;             // array of AliTRDCalDCSPTR objects
  //   TObjArray*       fGTUArr;        // array of AliTRDCalDCSGTU objects
  Int_t              fSystem;             // current system (FEE/PTR/GTU) (while parsing)
  Int_t              fInsideRstate;       // if we are inside rstate (while parsing)
  Int_t              fCurrentSM;          // current supermodule (while parsing)
  Int_t              fCurrentStack;       // current stack (while parsing)
  Int_t              fCurrentROB;         // current ROB (while parsing)
  Int_t              fCurrentMCM;         // current MCM (while parsing)
  Int_t              fCurrentADC;         // current ADC (while parsing)
  TString            fContent;            // content of the xml element (text)
  AliTRDCalDCSFEEv2* fDCSFEEObj;          // the calib object for one FEE DCS board
  AliTRDCalDCSPTR*   fDCSPTRObj;          // the calib object for one PTR DCS board
  AliTRDCalDCSGTU*   fDCSGTUObj;          // the calib object for one GTU DCS board
  AliTRDCalDCSv2*    fCalDCSObj;          // the complete calib obj containing all info
  Int_t              fLevel1Tag;          // 1st level in XML (while parsing)
  Int_t              fLevel2Tag;          // 2nd level in XML (while parsing)
  Bool_t             fInsideBoardInfo;    // if we are inside BoardInfo (while parsing)

  AliTRDCalDCSGTUTmu*       fTmu;       // GTU calibration data: pattern generator
  AliTRDCalDCSGTUCtpOpc*    fCtpOpc;    // GTU calibration data: OPC
  AliTRDCalDCSGTUSegment*   fSegment;   // GTU calibration data: SMU tracklets/tracks/triggers
  AliTRDCalDCSGTUBoardInfo* fBoardInfo; // GTU calibration data: hard-/software and type

  ClassDef(AliTRDSaxHandler,3);         // The XML file handler for the preprocessor
};
#endif

