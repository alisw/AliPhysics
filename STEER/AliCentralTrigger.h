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
//    Load Descriptors                                                       //
//    Make a list the trigger detectors involve from the descriptors         //
//    For the each event                                                     //
//           Run the Trigger for the each detector                           //
//           Get the inputs                                                  //
//           Check the condition classes                                     //
//           Create the class mask                                           //
//           Save result                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TObjArray.h>

class TTree;
class AliRunLoader;

class AliCentralTrigger : public TObject {

public:
                          AliCentralTrigger();
                          AliCentralTrigger( TString & descriptor );
                          AliCentralTrigger( const AliCentralTrigger& ctp );
               virtual   ~AliCentralTrigger();

                Bool_t    LoadDescriptor( TString & descriptor );
                Bool_t    RunTrigger( AliRunLoader * runloader );
                ULong64_t CheckConditions();
                  void    Reset();
                  void    DeleteDescriptors();
                  void    MakeBranch( TString name, TTree * tree );
  //  Getters
               TString    GetDetectors();
             ULong64_t    GetClassMask() const { return fClassMask; }
               UChar_t    GetClusterMask();
             TObjArray*   GetLoadedDescriptors() { return &fDescriptors; }
             TObjArray*   GetResultConditions();
                  void    Print( const Option_t* opt ="" ) const;
protected:
       //        TString    fRunCondition;     // Running modes Ej. Pb-Pb, p-p, p-A
             ULong64_t    fClassMask;          // UID ( bitwise OR of conditions mask )
             TObjArray    fDescriptors;        // Array of Trigger Descriptors (AliTriggerDescriptor)
             TObjArray    fInputs;             //! Array of Trigger Inputs

private:
                Bool_t    IsSelected( TString detName, TString& detectors ) const;

   ClassDef( AliCentralTrigger, 1 )  // class for running the Central Trigger Processor
};


#endif
