#ifndef ALIITSTRACKING_H
#define ALIITSTRACKING_H

#include <TObject.h>
#include <TList.h>

class TObjArray;
class TVector;
class TMatrix;
class AliITStrack;
class AliITS;
class AliITSRad;

class AliITStracking : public TObject {

//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
   
  Double_t  Rlayer[6];
       
public:
  
  AliITStracking() {;}

  AliITStracking(TList *trackITSlist,AliITStrack *reference,AliITS *obj,TObjArray *fpoints,
                 Double_t Ptref, Int_t **vettid, Bool_t flagvert, AliITSRad *rl );

  Int_t NewIntersection(AliITStrack &track, Double_t rk,Int_t layer, Int_t &ladder, Int_t &detector );
  Double_t PhiDef(Double_t x, Double_t y);

  void KalmanFilter(AliITStrack *newtrack, TVector &cluster, Double_t sigma[2]);  
  void KalmanFilterVert(AliITStrack *newtrack, TVector &cluster, Double_t sigma[2]);  
  		       
  ClassDef(AliITStracking,1)
};

#endif
