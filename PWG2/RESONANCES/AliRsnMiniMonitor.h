#ifndef AliRsnMiniMonitor_H
#define AliRsnMiniMonitor_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Monitors
//
////////////////////////////////////////////////////////////////////////////////

#include "TObjArray.h"

class AliRsnDaughter;
class AliRsnEvent;

class AliRsnMiniMonitor : public TNamed {
public:

   enum EType {
      kdEdxTPCvsP,        // TPC signal vs. momentum
      ktimeTOFvsPKaon,    // TOF time vs. momentum
      ktimeTOFvsPPion,    // TOF time vs. momentum
      ktimeTOFvsPProton,  // TOF time vs. momentum
      kTypes              // total number of cuts
   };

   AliRsnMiniMonitor();
   AliRsnMiniMonitor(const char *name, EType type, Int_t cutID);
   AliRsnMiniMonitor(const AliRsnMiniMonitor& copy);
   AliRsnMiniMonitor& operator=(const AliRsnMiniMonitor& copy);
   virtual ~AliRsnMiniMonitor() { }

   EType              GetType()   {return fType;}
   Int_t              GetCutID()  {return fCutID;}
   Int_t              GetListID() {return fListID;}
   
   void               SetType(EType type)  {fType = type;}
   void               SetCutID(Int_t id)   {fCutID = id;}

   static const char* Label(EType type); 
   Bool_t             Init(const char *name, TList *list);
   Bool_t             Fill(AliRsnDaughter *track, AliRsnEvent *event);

protected:

   EType      fType;     //  monitor type
   Int_t      fCutID;    //  ID for cut to be used
   Int_t      fListID;   //  histogram ID in the list
   TList     *fList;     //! global output list
                                       
   ClassDef(AliRsnMiniMonitor, 1)  //  AliRsnMiniMonitor class
};

inline const char* AliRsnMiniMonitor::Label(EType type)
{
   switch (type) {
      case kdEdxTPCvsP      : return "TPCsignal";
      case ktimeTOFvsPPion  : return "TOFsignalPi";
      case ktimeTOFvsPKaon  : return "TOFsignalK";
      case ktimeTOFvsPProton: return "TOFsignalP";
      default               : return "X";
   }
}

#endif
