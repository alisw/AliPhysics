#ifndef ALITRDSIMPLE_H
#define ALITRDSIMPLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id$ */
 
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Simplified TRD slow simulator                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include <TObject.h>

class AliTRDsimpleGen;
 
class AliTRDsimple : public TObject {
 
 public:     

  AliTRDsimple();
  AliTRDsimple(const AliTRDsimple &s); 
                                                                                
  virtual ~AliTRDsimple();
  AliTRDsimple &operator=(const AliTRDsimple &s);    

  virtual void             Init();
  virtual void             Copy(TObject &s) const;
  virtual void             ProcessEvent(Int_t ievent);

  virtual AliTRDsimpleGen *GetGenerator()      const { return fGenerator; };

 protected:

  AliTRDsimpleGen      *fGenerator;   //  The generator class for the simple simulator

  ClassDef(AliTRDsimple,1)            //  Simplified TRD slow simulator
 
};
#endif                                                                          
