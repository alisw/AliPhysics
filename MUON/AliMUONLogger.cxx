/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

#include "AliMUONLogger.h"

#include "AliMUONStringIntMap.h"
#include "AliLog.h"
#include "Riostream.h"

//-----------------------------------------------------------------------------
/// \class AliMUONLogger
///
/// A logger that keeps track of the number of times a message appeared.
///
/// Typically used to print all messages to screen at once, e.g. in the
/// dtor of a worker class.
///
/// For instance, it is used in AliMUONDigitizerV3, to note which channels
/// are disabled, and this information is printed in a condensed form
/// only once when DigitizerV3 is destroyed.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONLogger)
/// \endcond

//_____________________________________________________________________________
AliMUONLogger::AliMUONLogger(Int_t maxNumberOfEntries) 
: TObject(), 
  fMaxNumberOfEntries(maxNumberOfEntries),
  fLog(new AliMUONStringIntMap)
{
    /// ctor. After maxNumberOfEntries, the log is printed and reset
}

//_____________________________________________________________________________
AliMUONLogger::~AliMUONLogger()
{
  /// dtor
  delete fLog;
}

//_____________________________________________________________________________
Int_t 
AliMUONLogger::Log(const char* message)
{
  /// Log a message

  if ( fMaxNumberOfEntries >0 && fLog->GetNofItems() >= fMaxNumberOfEntries ) 
  {
    AliWarning(Form("Reached max number of entries (%d over %d). Printing and resetting.",
                    fLog->GetNofItems(),fMaxNumberOfEntries));
    Print();
    fLog->Clear();
  }
  
  Int_t i = fLog->Get(message);
  
  fLog->Set(message,i+1);
  
  return i+1;
}

//_____________________________________________________________________________
void   
AliMUONLogger::Clear(Option_t* /*option*/) 
{  
  /// reset logger spool

  fLog->Clear();
}

//_____________________________________________________________________________
void 
AliMUONLogger::Print(Option_t* opt) const
{
  /// Print the entire log
  if ( fLog->GetNofItems() )
  {
    fLog->Print(opt);
  }
  else
  {
    cout << "No message" << endl;
  }
}
  
//_____________________________________________________________________________
void
AliMUONLogger::Print(TString& key, ofstream& out) const
{
  /// print out into a given streamer with a key word in front of the message
  fLog->Print(key, out); 

      
}
   
//_____________________________________________________________________________
void 
AliMUONLogger::ResetItr()
{
  /// call reset iterator method
  fLog->ResetItr();
     
}
 
//_____________________________________________________________________________
Bool_t 
AliMUONLogger::Next(TString& msg, Int_t& occurance)
{
  /// call next iterator method
  return fLog->Next(msg, occurance);
     
}
    
//_____________________________________________________________________________
Int_t
AliMUONLogger::NumberOfEntries() const
{
  /// Get the number of logs we have so far
  return fLog->GetNofItems();
}

