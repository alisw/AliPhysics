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

#include <TFormula.h>
#include "AliRsnCut.h"

ClassImp(AliRsnCut)

//______________________________________________________________________________
AliRsnCut::AliRsnCut(const char *name, RSNTARGET target) :
   AliRsnTarget(name, target),
   fMinI(0),
   fMaxI(0),
   fMinD(0.),
   fMaxD(0.),
   fMinIptdep(0),
   fMaxIptdep(0),
   fMinDptdep(0),
   fMaxDptdep(0),
   fCutValueI(0),
   fCutValueD(0.0),
   fPtDepCut(kFALSE),
   fRefPtValueD(0.0),
   fMaxPt(1E20),
   fMinPt(0.0),
   fPtDepCutMaxFormula(""),
   fPtDepCutMinFormula(""),
   fCutResult(kTRUE)
{
//
// Default constructor.
//
}

//______________________________________________________________________________
AliRsnCut::AliRsnCut
(const char *name, RSNTARGET target, Long_t imin, Long_t imax, Double_t dmin, Double_t dmax) :
   AliRsnTarget(name, target),
   fMinI(imin),
   fMaxI(imax),
   fMinD(dmin),
   fMaxD(dmax),
   fMinIptdep(0),
   fMaxIptdep(0),
   fMinDptdep(0),
   fMaxDptdep(0),
   fCutValueI(0),
   fCutValueD(0.0),
   fPtDepCut(kFALSE),
   fRefPtValueD(0.0),
   fMaxPt(1E20),
   fMinPt(0.0),
   fPtDepCutMaxFormula(""),
   fPtDepCutMinFormula(""),
   fCutResult(kTRUE)
{
//
// Constructor with arguments.
// This is provided to allow a quick setting of all data members.
//
}

//______________________________________________________________________________
AliRsnCut::AliRsnCut
(const char *name, RSNTARGET target, Int_t imin, Int_t imax, Double_t dmin, Double_t dmax) :
   AliRsnTarget(name, target),
   fMinI(imin),
   fMaxI(imax),
   fMinD(dmin),
   fMaxD(dmax),
   fMinIptdep(0),
   fMaxIptdep(0),
   fMinDptdep(0),
   fMaxDptdep(0),
   fCutValueI(0),
   fCutValueD(0.0),
   fPtDepCut(kFALSE),
   fRefPtValueD(0.0),
   fMaxPt(1E20),
   fMinPt(0.0),
   fPtDepCutMaxFormula(""),
   fPtDepCutMinFormula(""),
   fCutResult(kTRUE)
{
//
// Constructor with arguments.
// This is provided to allow a quick setting of all data members.
//
}

//______________________________________________________________________________
AliRsnCut::AliRsnCut
(const char *name, RSNTARGET target, Double_t dmin, Double_t dmax, Int_t imin, Int_t imax) :
   AliRsnTarget(name, target),
   fMinI(imin),
   fMaxI(imax),
   fMinD(dmin),
   fMaxD(dmax),
   fMinIptdep(0),
   fMaxIptdep(0),
   fMinDptdep(0),
   fMaxDptdep(0),
   fCutValueI(0),
   fCutValueD(0.0),
   fPtDepCut(kFALSE),
   fRefPtValueD(0.0),
   fMaxPt(1E20),
   fMinPt(0.0),
   fPtDepCutMaxFormula(""),
   fPtDepCutMinFormula(""),
   fCutResult(kTRUE)
{
//
// Constructor with arguments.
// This is provided to allow a quick setting of all data members.
//
}

//______________________________________________________________________________
AliRsnCut::AliRsnCut(const AliRsnCut &copy) :
   AliRsnTarget(copy),
   fMinI(copy.fMinI),
   fMaxI(copy.fMaxI),
   fMinD(copy.fMinD),
   fMaxD(copy.fMaxD),
   fMinIptdep(copy.fMinIptdep),
   fMaxIptdep(copy.fMaxIptdep),
   fMinDptdep(copy.fMinDptdep),
   fMaxDptdep(copy.fMaxDptdep),
   fCutValueI(copy.fCutValueI),
   fCutValueD(copy.fCutValueD),
   fPtDepCut(copy.fPtDepCut),
   fRefPtValueD(copy.fRefPtValueD),
   fMaxPt(copy.fMaxPt),
   fMinPt(copy.fMinPt),
   fPtDepCutMaxFormula(copy.fPtDepCutMaxFormula),
   fPtDepCutMinFormula(copy.fPtDepCutMinFormula),
   fCutResult(copy.fCutResult)
{
//
// Copy constructor.
// Don't duplicate memory occupancy for pointer
//
}

//______________________________________________________________________________
AliRsnCut &AliRsnCut::operator=(const AliRsnCut &copy)
{
//
// Assignment operator.
// Don't duplicate memory occupancy for pointer
//

   AliRsnTarget::operator=(copy);
   if (this == &copy)
      return *this;

   fMinI      = copy.fMinI;
   fMaxI      = copy.fMaxI;
   fMinD      = copy.fMinD;
   fMaxD      = copy.fMaxD;
   fMinIptdep = copy.fMinIptdep;
   fMaxIptdep = copy.fMaxIptdep;
   fMinDptdep = copy.fMinDptdep;
   fMaxDptdep = copy.fMaxDptdep;
   fCutValueI = copy.fCutValueI;
   fCutValueD = copy.fCutValueD;
   fPtDepCut = copy.fPtDepCut;
   fRefPtValueD = copy.fRefPtValueD;
   fMaxPt = copy.fMaxPt;
   fMinPt = copy.fMinPt;
   fPtDepCutMaxFormula = copy.fPtDepCutMaxFormula;
   fPtDepCutMinFormula = copy.fPtDepCutMinFormula;
   fCutResult = copy.fCutResult;

   return (*this);
}

//______________________________________________________________________________
Bool_t AliRsnCut::IsSelected(TObject * /*object*/)
{
//
// Virtual cut-checking method.
// In this implementation, it does nothing, and all classes
// inheriting from this, should provide a proper implementation
// which must return kTRUE if the cut is passed, and kFALSE otherwise.
//

   AliWarning("This virtual function must be implemented properly");
   return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRsnCut::OkValueI()
{
//
// This method is used to compare a value with a reference.
// In the case of integers, the equality must be exact.
//

   // eval result
   
   if (fPtDepCut){
   	if(fRefPtValueD > fMaxPt) {
    	AliDebug(2,Form("pt = %f (> %f), cutting at %d\n",fRefPtValueD, fMaxPt, fMinI)); 
    	fCutResult = (fCutValueI == fMinI);
  	} else if (fRefPtValueD < fMinPt){
	AliDebug(2,Form("pt = %f (< %f), cutting at %d\n",fRefPtValueD, fMinPt, fMinI));
	fCutResult = (fCutValueI == fMinI);
	}else{ 
	TString str(fPtDepCutMinFormula);
        str.ReplaceAll("pt", "x");
        TFormula ptdepcut(Form("%s_ptdepcut", GetName()), str.Data());
        fMinIptdep = static_cast<int> (ptdepcut.Eval(fRefPtValueD));	
    	AliDebug(2,Form("pt = %f (> %f and < %f), cutting  at %d\n",fRefPtValueD, fMinPt, fMaxPt, fMinIptdep)); 
    	fCutResult = (fCutValueI == fMinIptdep);
  	}
  }
  else fCutResult = (fCutValueI == fMinI);

   // print debug message
   AliDebug(AliLog::kDebug + 2, "=== CUT DEBUG ========================================================");
   AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
   AliDebug(AliLog::kDebug + 2, Form("Checked value: %d", fCutValueI));
   AliDebug(AliLog::kDebug + 2, Form("Cut value    : %d", fMinI));
   AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
   AliDebug(AliLog::kDebug + 2, "=== END CUT DEBUG ====================================================");

   return fCutResult;
}

//______________________________________________________________________________
Bool_t AliRsnCut::OkValueD()
{
//
// This method is used to compare a value with a reference.
// In the case of doubles, the equality consists in being very close.
//

   // eval result
   
    if (fPtDepCut){
   	if(fRefPtValueD > fMaxPt) {
    	AliDebug(2,Form("pt = %f (> %f), cutting at %f\n",fRefPtValueD, fMaxPt, fMinD)); 
    	fCutResult = (TMath::Abs(fCutValueD - fMinD) < 1E-6);
  	} else if (fRefPtValueD < fMinPt){
	AliDebug(2,Form("pt = %f (< %f), cutting at %f\n",fRefPtValueD, fMinPt, fMinD));
	fCutResult = (TMath::Abs(fCutValueD - fMinD) < 1E-6);
	}else{
	TString str(fPtDepCutMinFormula);
        str.ReplaceAll("pt", "x");
        TFormula ptdepcut(Form("%s_ptdepcut", GetName()), str.Data());
        fMinDptdep = ptdepcut.Eval(fRefPtValueD);	
    	AliDebug(2,Form("pt = %f (> %f and < %f), cutting  at %f\n",fRefPtValueD, fMinPt, fMaxPt, fMinDptdep)); 
    	fCutResult = (TMath::Abs(fCutValueD - fMinDptdep) < 1E-6);
  	}
  }
  else fCutResult = (TMath::Abs(fCutValueD - fMinD) < 1E-6);

   // print debug message
   AliDebug(AliLog::kDebug + 2, "=== CUT DEBUG =======================================================");
   AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
   AliDebug(AliLog::kDebug + 2, Form("Checked value: %f", fCutValueD));
   AliDebug(AliLog::kDebug + 2, Form("Cut value    : %f", fMinD));
   AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
   AliDebug(AliLog::kDebug + 2, "=== END CUT DEBUG ===================================================");

   return fCutResult;
}

//______________________________________________________________________________
Bool_t AliRsnCut::OkRangeI()
{
//
// This method is used to compare a value with an integer range.
//

  // eval result
  if (fPtDepCut){
    if(fRefPtValueD > fMaxPt) {
      AliDebug(2,Form("pt = %f (> %f), cutting between [%d, %d]\n",fRefPtValueD, fMaxPt, fMinI, fMaxI)); 
      fCutResult = ((fCutValueI >= fMinI) && (fCutValueD <= fMaxI));
    } else if (fRefPtValueD < fMinPt){
      AliDebug(2,Form("pt = %f (< %f), cutting between [%d, %d]\n",fRefPtValueD, fMinPt, fMinI, fMaxI));
      fCutResult = ((fCutValueI >= fMinI) && (fCutValueD <= fMaxI));
    } else {
      TString str(fPtDepCutMinFormula);
      str.ReplaceAll("pt", "x");
      TFormula ptdepcut(Form("%s_ptdepcut", GetName()), str.Data());
      fMinIptdep = static_cast<int> (ptdepcut.Eval(fRefPtValueD));
	
      TString str2(fPtDepCutMaxFormula);
      str2.ReplaceAll("pt", "x");
      TFormula ptdepcut2(Form("%s_ptdepcut", GetName()), str2.Data());
      fMaxIptdep = static_cast<int> (ptdepcut2.Eval(fRefPtValueD));
    		    
      AliDebug(2,Form("pt = %f (> %f and < %f), cutting  according to the fiducial zone [%d, %d]\n",fRefPtValueD, fMinPt, fMaxPt, fMinIptdep, fMaxIptdep)); 
      fCutResult = ((fCutValueI >= fMinIptdep) && (fCutValueI <= fMaxIptdep));
    }
  }
  else fCutResult = ((fCutValueI >= fMinI) && (fCutValueI <= fMaxI));

   // print debug message
   AliDebug(AliLog::kDebug + 2, "=== CUT DEBUG ========================================================");
   AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
   AliDebug(AliLog::kDebug + 2, Form("Checked value: %d", fCutValueI));
   AliDebug(AliLog::kDebug + 2, Form("Cut range    : %d , %d", fMinI, fMaxI));
   AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
   AliDebug(AliLog::kDebug + 2, "=== END CUT DEBUG ====================================================");

   return fCutResult;
}

//______________________________________________________________________________
Bool_t AliRsnCut::OkRangeD()
{
//
// This method is used to compare a value with a double-float range.
//

  // eval result
   
  if (fPtDepCut){
    if(fRefPtValueD > fMaxPt) {
      AliDebug(2,Form("pt = %f (> %f), cutting between [%f, %f]\n",fRefPtValueD, fMaxPt, fMinD, fMaxD)); 
      fCutResult = ((fCutValueD >= fMinD) && (fCutValueD <= fMaxD));
    } else if (fRefPtValueD < fMinPt) {
      AliDebug(2,Form("pt = %f (< %f), cutting between [%f, %f]\n",fRefPtValueD, fMinPt, fMinD, fMaxD));
      fCutResult = ((fCutValueD >= fMinD) && (fCutValueD <= fMaxD));
    } else {
      TString str(fPtDepCutMinFormula);
      str.ReplaceAll("pt", "x");
      TFormula ptdepcut(Form("%s_ptdepcut", GetName()), str.Data());
      fMinDptdep = ptdepcut.Eval(fRefPtValueD);
      
      TString str2(fPtDepCutMaxFormula);
      str2.ReplaceAll("pt", "x");
      TFormula ptdepcut2(Form("%s_ptdepcut", GetName()), str2.Data());
      fMaxDptdep = ptdepcut2.Eval(fRefPtValueD);   
      
      AliDebug(2,Form("pt = %f (> %f and < %f), cutting  according to the fiducial zone [%f, %f]\n",fRefPtValueD, fMinPt, fMaxPt, fMinDptdep, fMaxDptdep)); 
      fCutResult = ((fCutValueD >= fMinDptdep) && (fCutValueD <= fMaxDptdep));
    }
  }
  else fCutResult = ((fCutValueD >= fMinD) && (fCutValueD <= fMaxD));

   // print debug message
   AliDebug(AliLog::kDebug + 2, "=== CUT DEBUG ========================================================");
   AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
   AliDebug(AliLog::kDebug + 2, Form("Checked value: %f", fCutValueD));
   AliDebug(AliLog::kDebug + 2, Form("Cut range    : %f , %f", fMinD, fMaxD));
   AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
   AliDebug(AliLog::kDebug + 2, "=== END CUT DEBUG ====================================================");

   return fCutResult;
}


//______________________________________________________________________________
void AliRsnCut::Print(Option_t *) const
{
//
// Override TObject::Print() method,
// and print some useful info about the cut general parameters.
//

   AliInfo("=== CUT DETAILS ====================================");
   AliInfo(Form("Cut name     : [%s]", GetName()));
   AliInfo(Form("Cut target   : [%s]", GetTargetTypeName()));
   AliInfo(Form("Cut edges [D]: [%f - %f]", fMinD, fMaxD));
   AliInfo(Form("Cut edges [I]: [%d - %d]", fMinI, fMaxI));
   AliInfo(Form("Cut pt dependent: %s", (fPtDepCut ? "YES" : "NO")));
   AliInfo("====================================================");
}
