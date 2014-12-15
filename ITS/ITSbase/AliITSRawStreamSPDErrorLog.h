#ifndef ALIITSRAWSTREAMSPDERRORLOG_H
#define ALIITSRAWSTREAMSPDERRORLOG_H

/* $Id$ */

///////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                            //
// For easier handling of error messages from AliITSRawStreamSPD.    //
// The purpose of this class is to make possible the switch to the   //
// AliRoot raw data parsing routines in the onlinte monitoring       //
// programs, like SPD-MOOD, and still keep all the old functionality.//
///////////////////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGText.h>

class AliITSRawStreamSPDErrorLog : public TObject {

 public:
  AliITSRawStreamSPDErrorLog();
  AliITSRawStreamSPDErrorLog(const AliITSRawStreamSPDErrorLog& logger);
  AliITSRawStreamSPDErrorLog& operator=(const AliITSRawStreamSPDErrorLog& logger);
  virtual ~AliITSRawStreamSPDErrorLog();

  enum    {kNrErrorCodes = 21};
  enum    {kTotal = 0};

  void    Reset();
  void    ProcessError(UInt_t errorCode, UInt_t eq, Int_t bytesRead, Int_t headersRead, const Char_t *errMess);
  void    AddMessage(const Char_t *errMess);

  void    ResetEvent();
  void    ProcessEvent(UInt_t eventNum);
  void    AddErrorMessagesFromCurrentEvent(UInt_t eventNum);
  void    SummarizeEvent(UInt_t eventNum);

  UInt_t  GetNrErrors(UInt_t errorCode, UInt_t eq);
  UInt_t  GetNrErrorsAllEq(UInt_t errorCode);
  UInt_t  GetNrErrorsTotal(UInt_t errorCode, UInt_t eq);
  UInt_t  GetNrErrorsTotalAllEq(UInt_t errorCode);

  void    SetByteOffset(UInt_t eq, Int_t size);
  void    SuppressErrorMessages(UInt_t errorCode, Bool_t suppr = kTRUE);
  void    SuppressErrorEq(UInt_t eq, Bool_t suppr = kTRUE);

  static  UInt_t GetNrErrorCodes(){return kNrErrorCodes;}

  TGraph* GetConsErrEvent(UInt_t errorCode, UInt_t eq);
  TGraph* GetConsErrPos(UInt_t errorCode, UInt_t eq);
  TH1F*   GetConsErrType(UInt_t eq);
  TH1F*   GetConsErrFraction(UInt_t eq);        // NB!!! Take care of deleting returned object later
  TH1F*   GetConsErrFractionUnScaled(UInt_t eq);
  TGText* GetText() {return fText;}
  TGText* GetTextThisEvent(UInt_t eq) {if (eq<20) return fTextTmp[eq]; else return NULL;}
  TGText* GetTextGeneralThisEvent() {return fTextTmpGeneral;}

  UInt_t  GetEventErrPosCounter(UInt_t errorCode, UInt_t eq);
  UInt_t  GetEventErrPos(UInt_t index, UInt_t errorCode, UInt_t eq);

 private:
  Int_t   fNErrors[kNrErrorCodes][20];          // number of errors for this event, for each code and eq
  Int_t   fNErrorsTotal[kNrErrorCodes][20];     // number of errors for all events, for each code and eq
  UInt_t  fNEvents[20];                         // number of events used, for each eq
  UInt_t  fErrEventCounter[kNrErrorCodes][20];  // event counter used when filling graph
  UInt_t  fErrPosCounter[kNrErrorCodes][20];    // event counter used when filling graph
  UInt_t  fErrPosTMPCounter[kNrErrorCodes][20]; // event counter used when filling graph
  Int_t   fByteOffset[20];                      // offset: how many bytes in the equipment header, for each eq
  Bool_t  fSuppressMess[kNrErrorCodes];         // do we suppress error messages for a specific error code?
  Bool_t  fSuppressEq[20];                      // do we suppress error messages for a specific eq?

  TGraph  *fConsErrEvent[kNrErrorCodes][20];    // graphs to display number of errors found in each event
  TGraph  *fConsErrPos[kNrErrorCodes][20];      // graphs to display number of bytes read for each error and event
  TGraph  *fConsErrPosTMP[kNrErrorCodes][20];   // temporary, to fill tgraph above for event
  TH1F    *fConsErrType[20];                    // histogram of how many errors for each error code
  TH1F    *fConsErrFraction[20];                // histogram of rate of events with errors for each error code

  TGText  *fText;                               // text buffer for all events analyzed so far
  TGText  *fTextTmp[20];                        // text buffer for this event (defined error codes)
  TGText  *fTextTmpGeneral;                     // text buffer for this event (general errors)

  void    InitHistograms();
  void    DeleteHistograms() ;

  ClassDef(AliITSRawStreamSPDErrorLog, 2);
};

#endif

