#ifndef ALIITSGEOINFO_H
#define ALIITSGEOINFO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

  #include <TROOT.h>

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class AliITSgeoinfo : public TObject {
 
 //Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.
    
	 public:
	 
    AliITSgeoinfo();
	 
    Int_t Nlad[6];
	 Int_t Ndet[6];
	 Float_t Avrad[6];
	 Float_t Detx[6];
	 Float_t Detz[6];
	 

    ClassDef(AliITSgeoinfo,1)
};

#endif
