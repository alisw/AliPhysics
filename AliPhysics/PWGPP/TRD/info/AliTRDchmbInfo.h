#ifndef ALITRDCHMBINFO_H
#define ALITRDCHMBINFO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Chamber Info Incapsulation                                            //
//                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TBox;
class TLatex;
class TCollection;
class AliTRDchmbInfo : public TNamed
{
public:  
  AliTRDchmbInfo();
  AliTRDchmbInfo(Int_t det, Int_t stat, Double_t p[4]);
  virtual ~AliTRDchmbInfo();

  void        Draw(Option_t* option = "");  // *MENU*
  Short_t     GetDetector() const          { return fDet;}
  const Double_t*   GetPosition() const    { return &fPosition[0];}
  Char_t      GetStatus() const            { return fStatus;}
  Long64_t    Merge(TCollection* /*list*/) { return 0;}
  void        Print(Option_t *o="") const;  // *MENU*
  void        SetDetector(Int_t det);
  void        SetPosition(Double_t p[4])   { memcpy(fPosition, p, 4*sizeof(Double_t));}
  void        SetStatus(Int_t stat)        { fStatus = stat;}

private:
  AliTRDchmbInfo(const AliTRDchmbInfo &ref);
  const AliTRDchmbInfo& operator=(const AliTRDchmbInfo &ref);

  Short_t     fDet;           // chamber no.
  Char_t      fStatus;        // status of chamber
  Double_t    fPosition[4];   // eta-phi position
  TBox        *fBox;          //! graph representation of chmb info
  TBox        *fShade;        //! graph representation of not OK status
  TLatex      *fLabel;        //! detector name
  ClassDef(AliTRDchmbInfo, 1) // TRD chamber position/status representation
};

#endif

