#ifndef ALICENTRALTRIGGER_H
#define ALICENTRALTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This class for running the Central Trigger Processor                      //
//                                                                           //
//                                                                           //
//    Load Configuration                                                     //
//    Make a list the trigger detectors involved ( from the configuration)   //
//    For the each event                                                     //
//           Run the Trigger for the each detector                           //
//           Get the inputs                                                  //
//           Check the trigger classes                                       //
//           Create the class mask                                           //
//           Save result                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TObjArray.h>

class TTree;
class AliRunLoader;
class AliTriggerConfiguration;

class AliCentralTrigger : public TObject {

public:
                          AliCentralTrigger();
                          AliCentralTrigger( TString & config );
               virtual   ~AliCentralTrigger();

                Bool_t    LoadConfiguration( TString & config );
                Bool_t    RunTrigger( AliRunLoader * runloader , const char* detectors);
             ULong64_t    TriggerClasses();
                  void    Reset();
		  void    DeleteConfiguration();
                  void    MakeBranch( TString name, TTree * tree );
  //  Getters
               TString    GetDetectors();
             ULong64_t    GetClassMask() const { return fClassMask; }
	       UChar_t    GetClusterMask() const { return fClusterMask; }
 AliTriggerConfiguration* GetConfiguration() { return fConfiguration; }
             TObjArray*   GetFiredClasses() const;
                  void    Print( const Option_t* opt ="" ) const;
	       TString    GetTriggeredDetectors() const;
protected:
             ULong64_t    fClassMask;          // UID ( bitwise OR of conditions mask )
               UChar_t    fClusterMask;        // UID ( bitwise OR of clusters mask )
 AliTriggerConfiguration* fConfiguration;      // Trigger Configuration used

private:
                void      SetOwner(Bool_t x=kTRUE){SetBit(22,x);} // Bit 22 indicates that the object owns fConfiguration
                Bool_t    IsOwner() const {return TestBit(22);} // Test bit 22 to check that the object owns fConfiguration
                Bool_t    IsSelected( TString detName, TString& detectors ) const;
		AliCentralTrigger( const AliCentralTrigger& ctp ); // Implemented
		AliCentralTrigger& operator=( const AliCentralTrigger& ctp ); // Not implemented

   ClassDef( AliCentralTrigger, 3 )  // class for running the Central Trigger Processor
};


#endif
