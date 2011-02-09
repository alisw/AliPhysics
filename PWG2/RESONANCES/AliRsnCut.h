//
// *** Class AliRsnCut ***
//
// Cut base class: all other cuts inherit from it.
// The 'core' of the class is the method "IsSelected()" which
// must be overloaded by any specific cut implementation.
//
// This class provides some default instruments to check values
// agains a reference or an allowed range, in order to permit
// a unique way to execute such kind of checks.
//
// authors: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          Martin Vala (martin.vala@cern.ch)
//

#ifndef ALIRSNCUT_H
#define ALIRSNCUT_H

#include "AliRsnTarget.h"

class AliRsnEvent;

class AliRsnCut : public AliRsnTarget {
public:

   AliRsnCut(const char *name = "dummy", RSNTARGET target = AliRsnTarget::kTargetTypes);
   AliRsnCut(const char *name, RSNTARGET target, Int_t    imin, Int_t    imax = 0 , Double_t dmin = 0., Double_t dmax = 0.);
   AliRsnCut(const char *name, RSNTARGET target, Double_t dmin, Double_t dmax = 0., Int_t    imin = 0 , Int_t    imax = 0);
   AliRsnCut(const AliRsnCut& copy);
   AliRsnCut& operator=(const AliRsnCut& copy);
   virtual ~AliRsnCut() { };

   Int_t            GetMinI()                  {return fMinI;}
   Int_t            GetMaxI()                  {return fMaxI;}
   Double_t         GetMinD()                  {return fMinD;}
   Double_t         GetMaxD()                  {return fMaxD;}
   Int_t            GetCutValueI()             {return fCutValueI;}
   Double_t         GetCutValueD()             {return fCutValueD;}
   Bool_t           GetCutResult()             {return fCutResult;}

   void             SetRangeI(Int_t    min, Int_t    max) {fMinI = min; fMaxI = max;}
   void             SetRangeD(Double_t min, Double_t max) {fMinD = min; fMaxD = max;}

   void             SetValueI(Int_t value)    {fMinI = value;}
   void             SetValueD(Double_t value) {fMinD = value;}

   Bool_t           OkValueI();
   Bool_t           OkRangeI();
   Bool_t           OkValueD();
   Bool_t           OkRangeD();

   virtual Bool_t   IsSelected(TObject *object);
   virtual void     Print(Option_t *opt = "") const;

protected:

   Int_t        fMinI;       //  lower edge of INT range or ref. value for INT CUT
   Int_t        fMaxI;       //  upper edge of INT range (not used for value cuts)
   Double_t     fMinD;       //  lower edge of DOUBLE range or ref. value for DOUBLE CUT
   Double_t     fMaxD;       //  upper edge of DOUBLE range (not used for value cuts)

   Int_t        fCutValueI;  //  cut value INT
   Double_t     fCutValueD;  //  cut value DOUBLE

   Bool_t       fCutResult;  //  tells if the cut is passed or not

   ClassDef(AliRsnCut, 1)
};

#endif
