#ifndef ALIFLATESDTRIGGER_H
#define ALIFLATESDTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/*
 * See implementation file for documentation
 */

/*
*/

#include "Rtypes.h"
#include "AliVMisc.h"
#include <string>

class AliFlatESDTrigger{

 public:

  // --------------------------------------------------------------------------------
  // -- Constructor / Destructors
  AliFlatESDTrigger();   
  ~AliFlatESDTrigger();

  // constructor and method for reinitialisation of virtual table
  AliFlatESDTrigger( AliVConstructorReinitialisationFlag );
  void Reinitialize() const {} // no virtual table - do nothing

  // --------------------------------------------------------------------------------
  // -- Fill / Set methods

  Int_t SetTriggerClass(  const char *TriggerClassName, Int_t TriggerIndex, ULong64_t MaxSize );
  
  // --------------------------------------------------------------------------------
  // -- Getter methods

  Int_t GetTriggerIndex() const { return fTriggerIndex; } 
  
  const Char_t *GetTriggerClassName() const { return reinterpret_cast<const Char_t*>( fContent ); } 
  
  // --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  
  
  const AliFlatESDTrigger *GetNextTrigger() const { return reinterpret_cast<const AliFlatESDTrigger*>(fContent+fContentSize); }
 
  AliFlatESDTrigger *GetNextTriggerNonConst() { return reinterpret_cast<AliFlatESDTrigger*>(fContent+fContentSize); }
 
  // --------------------------------------------------------------------------------
  // -- Size methods

  ULong64_t GetSize()  {return fContent + fContentSize -  reinterpret_cast<Byte_t*>(this) ;}
    
 private:

  AliFlatESDTrigger(const AliFlatESDTrigger&);
  AliFlatESDTrigger& operator=(const AliFlatESDTrigger&);

  // --------------------------------------------------------------------------------
  // -- Fixed size member objects
  
  UInt_t fContentSize;                      // Size of this object
  Int_t  fTriggerIndex; // trigger index  
   
  // --------------------------------------------------------------------------------
  // -- Variable Size Object

  Byte_t fContent[1];                  // Variale size object, which contains all data

};


// _______________________________________________________________________________________________________
// inline implementation of some methods 

inline AliFlatESDTrigger::AliFlatESDTrigger():
  fContentSize(1),
  fTriggerIndex(0)
{   
  // Default constructor
  fContent[0] = '\0';
}

inline AliFlatESDTrigger::~AliFlatESDTrigger() 
{
  // Destructor  
}

#pragma GCC diagnostic ignored "-Weffc++" 
inline AliFlatESDTrigger::AliFlatESDTrigger( AliVConstructorReinitialisationFlag ) {} // do nothing
#pragma GCC diagnostic warning "-Weffc++" 

inline Int_t AliFlatESDTrigger::SetTriggerClass(  const char *TriggerClassName, Int_t TriggerIndex, ULong64_t MaxSize )
{
  // Set trigger class, returns non-zero when the memory needed exeeeds MaxSize
	
	
  size_t len = strlen( TriggerClassName ) ;
	
	// strlen does not count the terminating \0 character, but this has to be safed too
	len ++;
    
  if( ( fContent + len ) > reinterpret_cast<Byte_t*>(this) + MaxSize ) return -1;
  
  fTriggerIndex = TriggerIndex;
  fContentSize =len;
  strcpy( reinterpret_cast<char*>(fContent), TriggerClassName );
  return 0;
}

#endif
