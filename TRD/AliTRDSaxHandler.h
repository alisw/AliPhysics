#ifndef AliTRDSaxHandler_H
#define AliTRDSaxHandler_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * * See cxx source for full Copyright notice */
/* $Id$ */


#include "TObject.h"
#include "Cal/AliTRDCalDCS.h"
#include "Cal/AliTRDCalDCSFEE.h"
#include "Cal/AliTRDCalDCSPTR.h"
#include "Cal/AliTRDCalDCSGTU.h"
#include <TObjArray.h>
#include <TList.h>


class AliTRDSaxHandler : public TObject {
public:
  enum { kInsideFEE = 1, kInsidePTR = 2, kInsideGTU = 3 };

  AliTRDSaxHandler();
  AliTRDSaxHandler(const AliTRDSaxHandler &sh);
  virtual ~AliTRDSaxHandler();
  AliTRDSaxHandler &operator=(const AliTRDSaxHandler &sh);

  TObjArray*    GetDCSFEEDataArray()  { return fFEEArr;        }
  TObjArray*    GetDCSPTRDataArray()  { return fPTRArr;        }
  TObjArray*    GetDCSGTUDataArray()  { return fGTUArr;        }
  AliTRDCalDCS* GetCalDCSObj(); // to be called by the preprocessor

  Int_t         GetHandlerStatus() const { return fHandlerStatus; }

  // functions for all possible events
  void          OnStartDocument();
  void          OnEndDocument();
  void          OnStartElement(const char *name, const TList *attributes);
  void          OnEndElement(const char *name);
  void          OnCharacters(const char *name);
  void          OnComment(const char *name);
  void          OnWarning(const char *name);
  void          OnError(const char *name);
  void          OnFatalError(const char *name);
  void          OnCdataBlock(const char *name, Int_t len);

private:

  Int_t            fHandlerStatus; // 0: everything OK, >0: error
  Int_t            fNDCSPTR;       // number of current PTR unit (to be abandonned soon)
  Int_t            fNDCSGTU;       // number of current GTU unit (to be abandonned soon)
  TObjArray*       fFEEArr;        // array of AliTRDCalDCSFEE objects
  TObjArray*       fPTRArr;        // array of AliTRDCalDCSPTR objects
  TObjArray*       fGTUArr;        // array of AliTRDCalDCSGTU objects
  Int_t            fSystem;        // current system (FEE/PTR/GTU)
  Int_t            fCurrentSM;     // current supermodule
  Int_t            fCurrentStack;  // current stack
  TString          fContent;       // content of the xml element (text) 
  AliTRDCalDCSFEE* fDCSFEEObj;     // the calib object for one FEE DCS board
  AliTRDCalDCSPTR* fDCSPTRObj;     // the calib object for one PTR DCS board
  AliTRDCalDCSGTU* fDCSGTUObj;     // the calib object for one GTU DCS board
  AliTRDCalDCS*    fCalDCSObj;     // the complete calib obj containing all inform.

  ClassDef(AliTRDSaxHandler,1);
};
#endif
