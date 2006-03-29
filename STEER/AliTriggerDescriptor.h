#ifndef ALITRIGGERDESCRIPTOR_H
#define ALITRIGGERDESCRIPTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This class for running a Trigger Descriptor                               //
//                                                                           //
//                                                                           //
// A Trigger Descriptor define a trigger setup for specific runnign
// condition (Pb-Pb, p-p, p-A, Calibration, etc).
// It keep:
//    - cluster detector (List of detectors involved)
//    - List of conditions                              
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class TNamed;

class TString;
class TObjArray;
class AliRunLoader;

class AliTriggerDescriptor : public TNamed {

public:
                          AliTriggerDescriptor();
                          AliTriggerDescriptor( TString & name, TString & description );
                          AliTriggerDescriptor( const AliTriggerDescriptor& des );
               virtual   ~AliTriggerDescriptor() { fConditions.SetOwner(); fConditions.Delete(); }
  AliTriggerDescriptor&   operator=(const AliTriggerDescriptor& des);

   //  Setters
                Bool_t    AddDetectorCluster( TString & cluster );
                  void    AddCondition( TString & cond,  TString & name,
                                        TString & description, Long_t mask  );
                  void    AddCondition( AliTriggerCondition* cond ) { fConditions.AddLast( cond ); }
  //  Getters
               TString    GetDetectorCluster() const { return fDetectorCluster; }
             TObjArray*   GetTriggerConditions() { return &fConditions; }
                Bool_t    CheckInputsConditions( TString & configfile );
                  void    Print( const Option_t* opt ="" ) const;
  //  Descriptors Database (root file)
                  void    WriteDescriptor( const char* filename="" );
      static TObjArray*   GetAvailableDescriptors( const char* filename="" );
      static
  AliTriggerDescriptor*   LoadDescriptor( TString & des, const char* filename="" );
  //TODO       static Bool_t    RemoveDescriptor( AliTriggerDescriptor* descriptor, const char* filename="" );
  //TODO       static Bool_t    RemoveDescriptor( TString* descriptor, const char* filename="" );

protected:
      //         TString    fRunCondition;       // Running modes Ej. Pb-Pb, p-p, p-A
               TString    fDetectorCluster;    // Array of Detector Trigger
             TObjArray    fConditions;         // Array of Trigger Condition (AliTriggerCondition)
               
    static const Int_t    fgkNDetectors = 9;              //! number possible trigger detectors
     static const char*   fgkDetectorName[fgkNDetectors]; //! names of detectors

private:
                Bool_t    IsSelected( TString detName, TString & detectors ) const;
  static const TString    fgkDescriptorFileName;        //! name of default descriptors file

   ClassDef( AliTriggerDescriptor, 1 )  // Define a trigger descriptor
};

#endif
