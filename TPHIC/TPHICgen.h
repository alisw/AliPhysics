#ifndef ROOT_TPHICGEN
#define ROOT_TPHICGEN
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
//------------------------------------------------------------------------
// TPHICgen is an interface class to fortran event generator of
// two-photon processes in ultraperipheral ion collisions
//%
// Yuri.Kharlov@cern.ch
// 15 April 2003
//------------------------------------------------------------------------

#include "TGenerator.h"

class TPHICgen : public TGenerator {

public:
  TPHICgen();
  virtual ~TPHICgen();

  void Initialize    ();
  void GenerateEvent ();
  void Finish        ();

  // Setters for COMMON /GGINI/
  void SetIPROC    (Int_t   iproc   );
  void SetNEVENT   (Int_t   nevent  );
  void SetILUMF    (Int_t   ilumf   );
  void SetLUMFIL   (TString lumfil  );
  void SetEBMN     (Float_t ebmn    );
  void SetIZ       (Int_t   iz      );
  void SetIA       (Int_t   ia      );
  void SetAMAS     (Float_t amas    );
  void SetAMIN     (Float_t amin    );
  void SetAMAX     (Float_t amax    );
  void SetYMIN     (Float_t ymin    );
  void SetYMAX     (Float_t ymax    );
  void SetNMAS     (Int_t   nmas    );
  void SetNY       (Int_t   ny      );
  void SetKFERM    (Int_t   kferm   );
  void SetKFONIUM  (Int_t   kfonium );
  void SetXMRES    (Float_t xmres   );
  void SetXGTRES   (Float_t xgtres  );
  void SetXGGRES   (Float_t xggres  );
  void SetMODDCY   (Int_t   moddcy  );
  void SetTHETAMIN (Float_t thetamin);
  void SetKV1      (Int_t   kv1     );
  void SetKV2      (Int_t   kv2     );

  // Getters for COMMON /GGEVNT/
  Float_t GetWSQ  () const;
  Float_t GetYGG  () const;
  Float_t GetXMG1 () const;
  Float_t GetXMG2 () const;
  Float_t GetP2G  (Int_t i) const;
  Float_t GetPTAG1(Int_t i) const;
  Float_t GetPTAG2(Int_t i) const;
  Int_t   GetNGG  () const;
  Int_t   GetKGG  (Int_t i) const;
  Float_t GetPGG  (Int_t i, Int_t j) const;

  // Getters for COMMON /GGXS/
  Float_t GetXSMAX0() const;
  Float_t GetXSCUR0() const;
  Float_t GetXSCUR () const;
  Float_t GetXSTOT () const;
  Float_t GetXSTOTE() const;

  ClassDef(TPHICgen,1)  //Interface to TPHIC Event Generator
};

#endif
