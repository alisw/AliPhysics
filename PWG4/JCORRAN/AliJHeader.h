// $Id: AliJHeader.h,v 1.1 2008/05/02 11:56:23 djkim Exp $

#ifndef ALIJHEADER_H
#define ALIJHEADER_H


////////////////////////////////////////////////////
/*!
  \file AliJHeader.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.1 $
  \date $Date: 2008/05/02 11:56:39 $

*/
////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "AliPhJBaseHeader.h"

//class TObject;

class AliJHeader : public AliPhJBaseHeader {
  
public:

  AliJHeader();              // default constructor
  AliJHeader(int eventid,
             short cent,
             float vrtz,
             ULong64_t triggmaskAli,
             UInt_t triggmaskJC,
             Int_t  refmult,
             UInt_t eventType
            );

  AliJHeader(const AliJHeader& a);                           

  virtual ~AliJHeader(){;}     // destructor
  
  ULong64_t  GetTriggerMaskAlice()   const {return fTriggerMaskAlice;}  
  UInt_t     GetTriggerMaskJCorran() const {return fTriggerMaskJCorran;}  
  Int_t      GetSPDTrackletMult()    const {return fSPDTrackletMult;}
  UInt_t     GetEventType()          const {return fEventType;}
  
  void SetTriggerMaskAlice(ULong64_t mask) {fTriggerMaskAlice = mask;}
  void SetTriggerMaskJCorran(UInt_t mask) {fTriggerMaskJCorran = mask;}
  void SetSPDTrackletMult(Int_t ref) { fSPDTrackletMult = ref;}
  void SetEventType(UInt_t eventype) {fEventType = eventype;}

  AliJHeader&  operator=(const AliJHeader& header);

 private:
 
  ULong64_t   fTriggerMaskAlice;         //Alice Trigger MASK
  UInt_t      fTriggerMaskJCorran;         // JCorran Trigger MASK
  Int_t       fSPDTrackletMult;             //SPD tracklet multiplicity
  UInt_t      fEventType;           // Type of Event
 
 
  ClassDef(AliJHeader,1)

};

#endif
