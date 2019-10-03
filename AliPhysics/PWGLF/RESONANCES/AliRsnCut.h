//
// *** Class AliRsnCut ***
//
// Cut base class: all other cuts inherit from it.
// The 'core' of the class is the method "IsSelected()" which
// must be overloaded by any specific cut implementation.
//
// This class provides some default instruments to check values
// against a reference or an allowed range, in order to permit
// a unique way to execute such kind of checks.
// Moreover, if one checks values and ranges using default methods
// a debug message can be printed on request.
//
// authors: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          Martin Vala (martin.vala@cern.ch)
//

#ifndef ALIRSNCUT_H
#define ALIRSNCUT_H

#include "AliRsnTarget.h"

class AliRsnCut : public AliRsnTarget {
public:

   AliRsnCut(const char *name = "dummy", RSNTARGET target = AliRsnTarget::kTargetTypes);
   AliRsnCut(const char *name, RSNTARGET target, Long_t    imin, Long_t    imax = 0 , Double_t dmin = 0., Double_t dmax = 0.);
   AliRsnCut(const char *name, RSNTARGET target, Int_t    imin, Int_t    imax = 0 , Double_t dmin = 0., Double_t dmax = 0.);
   AliRsnCut(const char *name, RSNTARGET target, Double_t dmin, Double_t dmax = 0., Int_t    imin = 0 , Int_t    imax = 0);
   AliRsnCut(const AliRsnCut &copy);
   AliRsnCut &operator=(const AliRsnCut &copy);
   virtual ~AliRsnCut() { };

   Int_t            GetMinI()      const {return fMinI;}
   Int_t            GetMaxI()      const {return fMaxI;}
   Double_t         GetMinD()      const {return fMinD;}
   Double_t         GetMaxD()      const {return fMaxD;}
   Int_t            GetCutValueI() const {return fCutValueI;}
   Double_t         GetCutValueD() const {return fCutValueD;}
   Bool_t           GetCutResult() const {return fCutResult;}

   void             SetRangeI(Int_t    min, Int_t    max) {fMinI = min; fMaxI = max;}
   void             SetRangeD(Double_t min, Double_t max) {fMinD = min; fMaxD = max;}

   void             SetValueI(Int_t value)    {fMinI = value;}
   void             SetValueD(Double_t value) {fMinD = value;}

   void             SetMaxPt(Double_t value) {fMaxPt = value;}
   void             SetMinPt(Double_t value) {fMinPt = value;}
   void		    SetPtDepCut(Bool_t flag) {fPtDepCut = flag;}
   void		    SetPtDepCutMaxFormula(const char *formula) {fPtDepCutMaxFormula = formula;}
   void		    SetPtDepCutMinFormula(const char *formula) {fPtDepCutMinFormula = formula;}

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
   
   Int_t        fMinIptdep;  //  lower edge of INT range or ref. value for INT CUT -- pt dependent
   Int_t        fMaxIptdep;  //  upper edge of INT range (not used for value cuts) -- pt dependent
   Double_t     fMinDptdep;  //  lower edge of DOUBLE range or ref. value for DOUBLE CUT -- pt dependent
   Double_t     fMaxDptdep;  //  upper edge of DOUBLE range (not used for value cuts) -- pt dependent

   Int_t        fCutValueI;  //  cut value INT
   Double_t     fCutValueD;  //  cut value DOUBLE
   Bool_t       fPtDepCut;     // flag to enable a pt dependent pair cut
   Double_t     fRefPtValueD;  //  pt value for a pt dependent pair cut
   Double_t     fMaxPt;        // maximum pt at which applying a pt dependent cut
   Double_t     fMinPt;        // minimum pt at which applying a pt dependent cut
   TString      fPtDepCutMaxFormula; //pt dependent cut formula
   TString      fPtDepCutMinFormula; //pt dependent cut formula

   Bool_t       fCutResult;  //  tells if the cut is passed or not

   ClassDef(AliRsnCut, 3)
};

#endif
