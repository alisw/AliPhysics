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
  void SetIPROC    (const Int_t   iproc   );
  void SetNEVENT   (const Int_t   nevent  );
  void SetILUMF    (const Int_t   ilumf   );
  void SetLUMFIL   (const TString lumfil  );
  void SetEBMN     (const Float_t ebmn    );
  void SetIZ       (const Int_t   iz      );
  void SetIA       (const Int_t   ia      );
  void SetAMAS     (const Float_t amas    );
  void SetAMIN     (const Float_t amin    );
  void SetAMAX     (const Float_t amax    );
  void SetYMIN     (const Float_t ymin    );
  void SetYMAX     (const Float_t ymax    );
  void SetNMAS     (const Int_t   nmas    );
  void SetNY       (const Int_t   ny      );
  void SetKFERM    (const Int_t   kferm   );
  void SetKFONIUM  (const Int_t   kfonium );
  void SetXMRES    (const Float_t xmres   );
  void SetXGTRES   (const Float_t xgtres  );
  void SetXGGRES   (const Float_t xggres  );
  void SetMODDCY   (const Int_t   moddcy  );
  void SetTHETAMIN (const Float_t thetamin);
  void SetKV1      (const Int_t   kv1     );
  void SetKV2      (const Int_t   kv2     );

  // Getters for COMMON /GGEVNT/
  Float_t GetWSQ  ()             ;
  Float_t GetYGG  ()             ;
  Float_t GetXMG1 ()             ;
  Float_t GetXMG2 ()             ;
  Float_t GetP2G  (const Int_t i);
  Float_t GetPTAG1(const Int_t i);
  Float_t GetPTAG2(const Int_t i);
  Int_t   GetNGG  ()             ;
  Int_t   GetKGG  (const Int_t i);
  Float_t GetPGG  (const Int_t i, const Int_t j);

  // Getters for COMMON /GGXS/
  Float_t GetXSMAX0()            ;
  Float_t GetXSCUR0()            ;
  Float_t GetXSCUR ()            ;
  Float_t GetXSTOT ()            ;
  Float_t GetXSTOTE()            ;

  ClassDef(TPHICgen,1)  //Interface to TPHIC Event Generator
};

#endif
