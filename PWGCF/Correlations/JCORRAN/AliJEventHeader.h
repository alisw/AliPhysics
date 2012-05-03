// $Id: AliJEventHeader.h,v 1.1 2008/05/02 11:56:23 djkim Exp $
////////////////////////////////////////////////////
/*!
  \file AliJEventHeader.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.1 $
  \date $Date: 2008/05/02 11:56:39 $

*/
////////////////////////////////////////////////////

#ifndef ALIJEVENTHEADER_H
#define ALIJEVENTHEADER_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "AliJBaseEventHeader.h"

class AliJEventHeader : public AliJBaseEventHeader {
 public:

  enum { kcV0M, kcFMD, kcTRK, kcTKL, kcCL0, kcCL1, kcV0MvsFMD, kcTKLvsV0, kcZEMvsZDC, kcNTYPE };

  /*
   * V0M = V0 multiplicity
   * FMD = FMD raw multiplicity
   * TRK = N. of tracks
   * TKL = N. of tracklets
   * CL0 = N. of clusters in layer 0
   * CL1 = N. of clusters in layer 1
   * V0MvsFMD = correlation between V0 and FMD
   * TKLvsV0 = correlation between tracklets and V0
   * ZEMvsZDC = correlation between ZEM and ZDC 
   */
  AliJEventHeader();              // default constructor
  AliJEventHeader(int eventid,
                  float cent,
                  float vrtz,
                  ULong64_t triggmaskAli,
                  UInt_t triggmaskJC,
                  Int_t  refmult,
                  Float_t  v0mult,
                  UInt_t eventType
                 );

  AliJEventHeader(const AliJEventHeader& a);                           

  virtual ~AliJEventHeader(){;}     // destructor

  ULong64_t  GetTriggerMaskAlice()   const {return fTriggerMaskAlice;}  
  UInt_t     GetTriggerMaskJCorran() const {return fTriggerMaskJCorran;}  
  Int_t      GetSPDTrackletMult()    const {return fSPDTrackletMult;}
  UInt_t     GetEventType()          const {return fEventType;}
  Float_t    GetV0Mult()             const {return fV0Mult;}
  Int_t      GetVtxMult()            const { return fVtxMult; };//FK// EFF

  Float_t GetCentralityArray( UInt_t it ) const { return it<kcNTYPE ? fCentralityArray[it] : -1; }

  void SetTriggerMaskAlice(ULong64_t mask) {fTriggerMaskAlice = mask;}
  void SetTriggerMaskJCorran(UInt_t mask) {fTriggerMaskJCorran = mask;}
  void SetSPDTrackletMult(Int_t ref) { fSPDTrackletMult = ref;}
  void SetEventType(UInt_t eventype) {fEventType = eventype;}
  void SetV0Mult(Float_t multV0) {fV0Mult = multV0;}
  void SetVtxMult(Int_t m){ fVtxMult = m; };//FK// EFF
  void SetCentralityArray(UInt_t it, Float_t cen ){ if( it < kcNTYPE ) fCentralityArray[it]=cen; }

  AliJEventHeader&  operator=(const AliJEventHeader& header);

 private:

  ULong64_t   fTriggerMaskAlice;           //Alice Trigger MASK
  UInt_t      fTriggerMaskJCorran;         // JCorran Trigger MASK
  Int_t       fSPDTrackletMult;             //SPD tracklet multiplicity
  Double32_t   fV0Mult;                   // VZERO multiplicity
  UInt_t      fEventType;                 // Type of Event
  Int_t       fVtxMult;                   //FK// EFF number of vertex contributors 
  Double32_t  fCentralityArray[kcNTYPE];  //?//

  ClassDef(AliJEventHeader,1)

};

#endif
