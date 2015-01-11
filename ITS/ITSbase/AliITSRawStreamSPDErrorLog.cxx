///////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                            //
// For easier handling of error messages from AliITSRawStreamSPD.    //
// The purpose of this class is to make possible the switch to the   //
// AliRoot raw data parsing routines in the onlinte monitoring       //
// programs, like SPD-MOOD, and still keep all the old functionality.//
///////////////////////////////////////////////////////////////////////

#include "AliITSRawStreamSPDErrorLog.h"
#include "AliLog.h"

/* $Id$ */

ClassImp(AliITSRawStreamSPDErrorLog)
//________________________________________________________________________________________________
AliITSRawStreamSPDErrorLog::AliITSRawStreamSPDErrorLog() :
  fText(NULL), fTextTmpGeneral(NULL)
{
  // default constructor
  fText=new TGText();
  fTextTmpGeneral=new TGText();
  for (UInt_t eq=0; eq<20; eq++) {
    fTextTmp[eq] = new TGText();
    fNEvents[eq] = 0;
    fByteOffset[eq] = 0;
    fSuppressEq[eq] = kFALSE;
    for (UInt_t err=0; err<kNrErrorCodes; err++) {
      fNErrors[err][eq] = 0;
      fNErrorsTotal[err][eq] = 0;
    }
  }
  for (UInt_t err=0; err<kNrErrorCodes; err++) {
    fSuppressMess[err] = kFALSE;
  }
  InitHistograms();
}
//________________________________________________________________________________________________
AliITSRawStreamSPDErrorLog::AliITSRawStreamSPDErrorLog(const AliITSRawStreamSPDErrorLog& logger) :
  TObject(logger), fText(NULL), fTextTmpGeneral(NULL)
{
  // copy constructor
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t err=0; err<kNrErrorCodes; err++) {
      fNErrors[err][eq] = logger.fNErrors[err][eq];
      fNErrorsTotal[err][eq] = logger.fNErrorsTotal[err][eq];
      fErrEventCounter[err][eq] = logger.fErrEventCounter[err][eq];
    }
    fNEvents[eq] = logger.fNEvents[eq];
    fByteOffset[eq] = logger.fByteOffset[eq];
    fSuppressEq[eq] = logger.fSuppressEq[eq];
  }
  for (UInt_t err=0; err<kNrErrorCodes; err++) {
    fSuppressMess[err] = logger.fSuppressMess[err];
  }
  
  fText = new TGText(logger.fText);
  fTextTmpGeneral = new TGText(logger.fTextTmpGeneral);
  for (UInt_t eq=0; eq<20; eq++) {
    fTextTmp[eq] = new TGText(logger.fTextTmp[eq]);
  }
  
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t err=0; err<kNrErrorCodes; err++) {
      fConsErrEvent[err][eq] = new TGraph(*logger.fConsErrEvent[err][eq]);
      fConsErrPos[err][eq] = new TGraph(*logger.fConsErrPos[err][eq]);
    }
    fConsErrType[eq] = new TH1F(*logger.fConsErrType[eq]);
    fConsErrType[eq]->SetDirectory(0);
    fConsErrFraction[eq] = new TH1F(*logger.fConsErrFraction[eq]);
    fConsErrFraction[eq]->SetDirectory(0);
  }
  
}
//________________________________________________________________________________________________
AliITSRawStreamSPDErrorLog& AliITSRawStreamSPDErrorLog::operator=(const AliITSRawStreamSPDErrorLog& logger) {
  // assignment operator
  if (this!=&logger) {
    delete fText;
    delete fTextTmpGeneral;
    for(Int_t eq=0; eq<20;eq++){
    delete fTextTmp[eq];
    }
    this->DeleteHistograms();

    for (UInt_t eq=0; eq<20; eq++) {
      for (UInt_t err=0; err<kNrErrorCodes; err++) {
	fNErrors[err][eq] = logger.fNErrors[err][eq];
	fNErrorsTotal[err][eq] = logger.fNErrorsTotal[err][eq];
	fErrEventCounter[err][eq] = logger.fErrEventCounter[err][eq];
      }
      fNEvents[eq] = logger.fNEvents[eq];
      fByteOffset[eq] = logger.fByteOffset[eq];
      fSuppressEq[eq] = logger.fSuppressEq[eq];
    }
    for (UInt_t err=0; err<kNrErrorCodes; err++) {
      fSuppressMess[err] = logger.fSuppressMess[err];
    }

    fText = new TGText(logger.fText);
    fTextTmpGeneral = new TGText(logger.fTextTmpGeneral);
    for (UInt_t eq=0; eq<20; eq++) {
      fTextTmp[eq] = new TGText(logger.fTextTmp[eq]);
    }
    
    for (UInt_t eq=0; eq<20; eq++) {
      for (UInt_t err=0; err<kNrErrorCodes; err++) {
	fConsErrEvent[err][eq] = new TGraph(*logger.fConsErrEvent[err][eq]);
	fConsErrPos[err][eq] = new TGraph(*logger.fConsErrPos[err][eq]);
      }
      fConsErrType[eq] = new TH1F(*logger.fConsErrType[eq]);
      fConsErrType[eq]->SetDirectory(0);
      fConsErrFraction[eq] = new TH1F(*logger.fConsErrFraction[eq]);
      fConsErrFraction[eq]->SetDirectory(0);
    }

  }
  return *this;
  
}
//________________________________________________________________________________________________
AliITSRawStreamSPDErrorLog::~AliITSRawStreamSPDErrorLog() {
  // destructor
  DeleteHistograms();
  delete fText;
  delete fTextTmpGeneral;
  for(Int_t eq=0; eq<20; eq++){
  delete fTextTmp[eq];
  }
}
//________________________________________________________________________________________________
void AliITSRawStreamSPDErrorLog::InitHistograms() {
  // initialize histograms
  for (UInt_t eq=0; eq<20; eq++) {
    TString histName, histTitle;
    histName = Form("ConsErrType %d",eq);
    histTitle = Form("Distribution of errors, eq %d",eq);
    fConsErrType[eq]=new TH1F(histName.Data(),histTitle.Data(),kNrErrorCodes,-0.5,kNrErrorCodes-0.5);
    fConsErrType[eq]->SetDirectory(0);
    fConsErrType[eq]->SetXTitle("Error Code");
    fConsErrType[eq]->SetYTitle("Nr Errors");
    
    histName = Form("ConsErrFraction %d",eq);
    histTitle = Form("Fraction of events with errors, eq %d",eq);
    fConsErrFraction[eq]=new TH1F(histName.Data(),histTitle.Data(),kNrErrorCodes,-0.5,kNrErrorCodes-0.5);
    fConsErrFraction[eq]->SetDirectory(0);
    fConsErrFraction[eq]->SetXTitle("Error Code");
    fConsErrFraction[eq]->SetYTitle("Fraction");
    
    for (UInt_t err=0; err<kNrErrorCodes; err++) {
      fConsErrEvent[err][eq]=new TGraph();
      fConsErrPos[err][eq]=new TGraph();
      fConsErrPosTMP[err][eq]=new TGraph();
      fErrEventCounter[err][eq] = 0;
      fErrPosCounter[err][eq] = 0;
      fErrPosTMPCounter[err][eq] = 0;
    }
  }
}
//________________________________________________________________________________________________
void AliITSRawStreamSPDErrorLog::DeleteHistograms() {
  // delete histograms
  for (UInt_t eq=0; eq<20; eq++) {
    delete fConsErrType[eq];
    delete fConsErrFraction[eq];
    for (UInt_t err=0; err<kNrErrorCodes; err++) {
      delete fConsErrEvent[err][eq];
      delete fConsErrPos[err][eq];
      delete fConsErrPosTMP[err][eq];
    }
  }
}
//________________________________________________________________________________________________
void AliITSRawStreamSPDErrorLog::Reset() {
  // Reset
  fText->Clear();
  fTextTmpGeneral->Clear();
  for (UInt_t eq=0; eq<20; eq++) {
    fTextTmp[eq]->Clear();
    fConsErrType[eq]->Reset();
    fConsErrFraction[eq]->Reset();
    fNEvents[eq] = 0;
    for (UInt_t err=0; err<kNrErrorCodes; err++) {
      fNErrors[err][eq] = 0;
      fNErrorsTotal[err][eq] = 0;
      delete fConsErrEvent[err][eq];
      delete fConsErrPos[err][eq];
      delete fConsErrPosTMP[err][eq];
      fErrEventCounter[err][eq] = 0;
      fErrPosCounter[err][eq] = 0;
      fErrPosTMPCounter[err][eq] = 0;
      fConsErrEvent[err][eq] = new TGraph();
      fConsErrPos[err][eq] = new TGraph();
      fConsErrPosTMP[err][eq] = new TGraph();
    }
  }
}
//________________________________________________________________________________________________
void AliITSRawStreamSPDErrorLog::ProcessError(UInt_t errorCode, UInt_t eq, Int_t bytesRead, Int_t headersRead, const Char_t *errMess) {
  // Process an error
  if (eq>=20) {
    AliWarning(Form("Equipment number (%d) out of bounds",eq));
    return;
  }
  // check if we want to exclude the error...
  if (!(fSuppressMess[errorCode] || fSuppressEq[eq])) {
    fNErrors[errorCode][eq]++;
    fNErrorsTotal[errorCode][eq]++;
    if (errorCode!=kTotal) {
      fNErrors[kTotal][eq]++;
      fNErrorsTotal[kTotal][eq]++;
    }

    if (bytesRead>=0) {
      fConsErrPosTMP[errorCode][eq]->SetPoint(fErrPosTMPCounter[errorCode][eq],0,bytesRead+fByteOffset[eq]);
      fErrPosTMPCounter[errorCode][eq]++;
      if (errorCode!=kTotal) {
	fConsErrPosTMP[kTotal][eq]->SetPoint(fErrPosTMPCounter[kTotal][eq],0,bytesRead+fByteOffset[eq]);
	fErrPosTMPCounter[kTotal][eq]++;
      }
    }

    TString msg;
    if (bytesRead<0) {
      msg = Form("%s",errMess);
    }
    else {
      msg = Form("%s (%d bytes read, %d chip headers read)",errMess,bytesRead+fByteOffset[eq],headersRead);
    }
    fTextTmp[eq]->InsLine(fTextTmp[eq]->RowCount(),msg.Data());
  }
}
//________________________________________________________________________________________________
void AliITSRawStreamSPDErrorLog::AddMessage(const Char_t *errMess) {
  // add a general message to the buffer
  fTextTmpGeneral->InsLine(fTextTmpGeneral->RowCount(),errMess);
}
//________________________________________________________________________________________________
void AliITSRawStreamSPDErrorLog::ResetEvent() {
  // reset error counters for current event
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t err=0; err<kNrErrorCodes; err++) {
      fNErrors[err][eq]=0;
    }
    fByteOffset[eq]=0;
    fTextTmpGeneral->Clear();
    fTextTmp[eq]->Clear();
  }

  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t err=0; err<kNrErrorCodes; err++) {
      delete fConsErrPosTMP[err][eq];
      fErrPosTMPCounter[err][eq] = 0;
      fConsErrPosTMP[err][eq] = new TGraph();
    }
  }
}
//________________________________________________________________________________________________
void AliITSRawStreamSPDErrorLog::AddErrorMessagesFromCurrentEvent(UInt_t eventNum) {
  TString msg;
  if (fText->RowCount()>5000) {
    fText->Clear();
    msg = "*** previous errors cleared ***";
    fText->InsLine(fText->RowCount(),msg.Data());
  }
  if (fTextTmpGeneral->RowCount()>1) {
    msg = Form("*** Event %d , General Errors: ***",eventNum);
    fText->InsLine(fText->RowCount(),msg.Data());
    fText->AddText(fTextTmpGeneral);
    //    fTextTmpGeneral->Clear();
  }

  for (UInt_t eq=0; eq<20; eq++) {
    if (fNErrors[kTotal][eq]>0) {
      msg = Form("*** Event %d , Eq %d: ***",eventNum,eq);
      fText->InsLine(fText->RowCount(),msg.Data());
      fText->AddText(fTextTmp[eq]);
      fText->InsLine(fText->RowCount()," ");
    }
    //    fTextTmp[eq]->Clear();
  }
}
//________________________________________________________________________________________________
void AliITSRawStreamSPDErrorLog::ProcessEvent(UInt_t eventNum) {
  // summarize the information for the current event
//  TString msg;
//  if (fText->RowCount()>5000) {
//    fText->Clear();
//    msg = "*** previous errors cleared ***";
//    fText->InsLine(fText->RowCount(),msg.Data());
//  }
//  if (fTextTmpGeneral->RowCount()>1) {
//    msg = Form("*** Event %d , General Errors: ***",eventNum);
//    fText->InsLine(fText->RowCount(),msg.Data());
//    fText->AddText(fTextTmpGeneral);
//    fTextTmpGeneral->Clear();
//  }

  AddErrorMessagesFromCurrentEvent(eventNum);

  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t err=0; err<kNrErrorCodes; err++) {
      fConsErrType[eq]->Fill(err,fNErrors[err][eq]);
      if (fNErrors[err][eq]>0) {
	fConsErrEvent[err][eq]->SetPoint(fErrEventCounter[err][eq],eventNum,fNErrors[err][eq]);
	fErrEventCounter[err][eq]++;
	fConsErrFraction[eq]->Fill(err);
      }
      for (UInt_t pind=0; pind<fErrPosTMPCounter[err][eq]; pind++) {
	Double_t x,y;
	fConsErrPosTMP[err][eq]->GetPoint(pind,x,y);
	fConsErrPos[err][eq]->SetPoint(fErrPosCounter[err][eq],eventNum,y);
	fErrPosCounter[err][eq]++;
      }
    }

    fNEvents[eq]++;
//    if (fNErrors[kTotal][eq]>0) {
//      msg = Form("*** Event %d , Eq %d: ***",eventNum,eq);
//      fText->InsLine(fText->RowCount(),msg.Data());
//      fText->AddText(fTextTmp[eq]);
//      fText->InsLine(fText->RowCount()," ");
//    }
//    fByteOffset[eq]=0;
//    fTextTmp[eq]->Clear();
  }

}
//________________________________________________________________________________________________
void AliITSRawStreamSPDErrorLog::SummarizeEvent(UInt_t eventNum) {
  // summarize the information for the current event 
  ProcessEvent(eventNum);
  ResetEvent();
}
//________________________________________________________________________________________________
void AliITSRawStreamSPDErrorLog::SetByteOffset(UInt_t eq, Int_t size) {
  // set byte offset (equipment header size)
  if (eq<20) {
    fByteOffset[eq]=size;
  }
  else {
    AliWarning(Form("Equipment number (%d) out of range",eq));
  }
}
//________________________________________________________________________________________________
void AliITSRawStreamSPDErrorLog::SuppressErrorMessages(UInt_t errorCode, Bool_t suppr) {
  // suppress error messages for specific error code
  if (errorCode<kNrErrorCodes) {
    fSuppressMess[errorCode] = suppr;
  }
  else {
    AliWarning(Form("Error code (%d) out of range",errorCode));
  }
}
//________________________________________________________________________________________________
void AliITSRawStreamSPDErrorLog::SuppressErrorEq(UInt_t eq, Bool_t suppr) {
  // suppress error messages for specific equipment
  if (eq<20) {
    fSuppressEq[eq] = suppr;
  }
  else {
    AliWarning(Form("Eq id (%d) out of range",eq));
  }
}
//________________________________________________________________________________________________
UInt_t AliITSRawStreamSPDErrorLog::GetNrErrors(UInt_t errorCode, UInt_t eq) {
  // returns the number of errors for this event for a specific error code and equipment
  if (errorCode<kNrErrorCodes && eq<20) return fNErrors[errorCode][eq];
  else {
    AliWarning(Form("Error code (%d) or equipment (%d) out of range, returning 0",errorCode,eq));
    return 0;
  }
}
//________________________________________________________________________________________________
UInt_t AliITSRawStreamSPDErrorLog::GetNrErrorsAllEq(UInt_t errorCode) {
  // returns the number of errors for this event for a specific error code and all equipments
  if (errorCode<kNrErrorCodes) {
    UInt_t returnval=0;
    for (UInt_t eq=0; eq<20; eq++) {
      returnval += fNErrors[errorCode][eq];
    }
    return returnval;
  }
  else {
    AliWarning(Form("Error code (%d) out of range, returning 0",errorCode));
    return 0;
  }
}
//________________________________________________________________________________________________
UInt_t AliITSRawStreamSPDErrorLog::GetNrErrorsTotal(UInt_t errorCode, UInt_t eq) {
  // returns the total number of errors for a specific error code and equipment
  if (errorCode<kNrErrorCodes && eq<20) {
    return fNErrorsTotal[errorCode][eq];
    //    return (UInt_t) fConsErrType[eq]->GetBinContent(errorCode+1);
  }
  else {
    AliWarning(Form("Error code (%d) or equipment (%d) out of range, returning 0",errorCode,eq));
    return 0;
  }
}
//________________________________________________________________________________________________
UInt_t AliITSRawStreamSPDErrorLog::GetNrErrorsTotalAllEq(UInt_t errorCode) {
  // returns the total number of errors for a specific error code and for all equipments
  if (errorCode<kNrErrorCodes) {
    UInt_t returnval=0;
    for (UInt_t eq=0; eq<20; eq++) {
      returnval += fNErrorsTotal[errorCode][eq];
      //      returnval += (UInt_t) fConsErrType[eq]->GetBinContent(errorCode+1);
    }
    return returnval;
  }
  else {
    AliWarning(Form("Error code (%d) out of range, returning 0",errorCode));
    return 0;
  }
}
//________________________________________________________________________________________________
TGraph* AliITSRawStreamSPDErrorLog::GetConsErrEvent(UInt_t errorCode, UInt_t eq) {
  // returns a pointer to the graph
  if (errorCode<kNrErrorCodes && eq<20) return fConsErrEvent[errorCode][eq];
  else {
    AliWarning(Form("Error code (%d) or equipment (%d) out of range, returning NULL",errorCode,eq));
    return NULL;
  }
}
//________________________________________________________________________________________________
TH1F* AliITSRawStreamSPDErrorLog::GetConsErrType(UInt_t eq) {
  // returns a pointer to the histogram
  if (eq<20) return fConsErrType[eq];
  else {
    AliWarning(Form("Eq nr (%d) out of bounds",eq));
    return NULL;
  }
}
//________________________________________________________________________________________________
TH1F* AliITSRawStreamSPDErrorLog::GetConsErrFraction(UInt_t eq) {
  // creates a new histogram and returns a pointer to it. 
  // NB!!! Take care of deleting this object later
  if (eq<20) {
    TH1F* returnhisto = (TH1F*)(fConsErrFraction[eq]->Clone());
    if (fNEvents[eq]!=0) returnhisto->Scale(1./fNEvents[eq]);
    //    returnhisto->SetMaximum(1.);
    return returnhisto;
  }
  else {
    AliWarning(Form("Eq nr (%d) out of bounds",eq));
    return NULL;
  }
}
//________________________________________________________________________________________________
TH1F* AliITSRawStreamSPDErrorLog::GetConsErrFractionUnScaled(UInt_t eq) {
  // returns a pointer to the histogram
  if (eq<20) return fConsErrFraction[eq];
  else {
    AliWarning(Form("Eq nr (%d) out of bounds",eq));
    return NULL;
  }
}
//________________________________________________________________________________________________
TGraph* AliITSRawStreamSPDErrorLog::GetConsErrPos(UInt_t errorCode, UInt_t eq) {
  // returns a pointer to the histogram
  if (errorCode<kNrErrorCodes && eq<20) return fConsErrPos[errorCode][eq];
  else {
    AliWarning(Form("Error code (%d) or equipment (%d) out of range, returning NULL",errorCode,eq));
    return NULL;
  }
}
//________________________________________________________________________________________________
UInt_t AliITSRawStreamSPDErrorLog::GetEventErrPosCounter(UInt_t errorCode, UInt_t eq) {
  // returns the number of errors with bytes information for current event
  if (errorCode<kNrErrorCodes && eq<20) return fErrPosTMPCounter[errorCode][eq];
  else {
    AliError(Form("Error code (%d) or equipment (%d) out of range, returning 0",errorCode,eq));
    return 0;
  }
}
//________________________________________________________________________________________________
UInt_t AliITSRawStreamSPDErrorLog::GetEventErrPos(UInt_t index, UInt_t errorCode, UInt_t eq) {
  // returns the bytes read for the index'th error of current event
  if (errorCode<kNrErrorCodes && eq<20) {
    if (index<fErrPosTMPCounter[errorCode][eq]) {
      Double_t x,y;
      fConsErrPosTMP[errorCode][eq]->GetPoint(index,x,y);
      return (UInt_t) y;
    }
    else {
      AliError(Form("Index %d out of range, returning 0",index));
    return 0;
    }
  }
  else {
    AliError(Form("Error code (%d) or equipment (%d) out of range, returning 0",errorCode,eq));
    return 0;
  }

}

