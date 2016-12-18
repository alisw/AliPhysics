#ifndef ALITRACKLETDNDETAUTILS_H
#define ALITRACKLETDNDETAUTILS_H
#ifndef __CINT__
# include <TH2.h>
# include <TH1.h>
# include <TList.h>
# include <TParameter.h>
# include <TError.h>
# include <TMath.h>
# include <TDirectory.h>
# include <THashList.h>
#else
class TList;
class TH1;
class TH2;
class TAxis;
class TDirectory;
#endif

/**
 * Class with utlity functions 
 * 
 */
class AliTrackletdNdetaUtils
{
public:
  typedef TList Container;
  /** 
   * Dummy constructor 
   */
  AliTrackletdNdetaUtils() {}
  /**
   * Dummy destructor 
   */
  virtual ~AliTrackletdNdetaUtils() {}
  //_________________________________________________________________________
  /** 
   * @{
   * @name Service functions to check histogram consistencies 
   */
  /** 
   * Check if both axis have the same number of bins 
   * 
   * @param which Which axis is being checked 
   * @param a1    First axis
   * @param a2    Second axis 
   * 
   * @return true of both axis have equal number of bins 
   */
  static Bool_t CheckAxisNBins(const char*  which,
			       const TAxis* a1,
			       const TAxis* a2);
  /** 
   * Check axis limits (min,max) are the same 
   * 
   * @param which Which axis we're checking 
   * @param a1    First axis 
   * @param a2    Second axis 
   * 
   * @return true if axis min/max are indetical 
   */
  static Bool_t CheckAxisLimits(const char*  which,
				const TAxis* a1,
				const TAxis* a2);
  /** 
   * In case of non-uniform bins, check each bin edge to see if they
   * are the same
   *  
   * @param which Which axis we're checking  
   * @param a1    First axis 
   * @param a2    Second axis 
   * 
   * @return true if all bin edges are the same 
   */
  static Bool_t CheckAxisBins(const char*  which,
			      const TAxis* a1,
			      const TAxis* a2);
  /** 
   * Check if all bin labels (if specified) are the same 
   * 
   * @param which Which axis we're checking 
   * @param a1    First axis 
   * @param a2    Second axis 
   * 
   * @return True if all bin labels are the same 
   */
  static Bool_t CheckAxisLabels(const char*  which,
				const TAxis* a1,
				const TAxis* a2);
  /** 
   * Check that axis definitions are compatible 
   * 
   * @param which    Which axis we're checking 
   * @param a1       First axis 
   * @param a2       Second axis 
   * @param alsoLbls If true, also check labels 
   *  
   * @return true if the axes are compatible 
   */
  static Bool_t CheckAxis(const char*  which, 
			  const TAxis* a1,
			  const TAxis* a2,
			  Bool_t       alsoLbls);
  /** 
   * Check if two histograms are compatible by checking each defined
   * axis.
   * 
   * @param h1 First histogram 
   * @param h2 Second histogram 
   * 
   * @return true of they are compatible 
   */
  static Bool_t CheckConsistency(const TH1* h1, const TH1* h2);
  /* @} */
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Utilities 
   */
  /** 
   * Get an object from a container. 
   * 
   * @param parent Container 
   * @param name   Name of object 
   * @param cls    If not null, check if object has this type 
   * 
   * @return Pointer to object or null
   */
  static TObject* GetO(Container* parent,const char* name,TClass* cls);
  /** 
   * Get an object from a directory. 
   * 
   * @param parent directory 
   * @param name   Name of object 
   * @param cls    If not null, check if object has this type 
   * 
   * @return Pointer to object or null
   */
  static TObject* GetO(TDirectory* parent,const char* name,TClass* cls);
  /** 
   * Get a 1D histogram from a container 
   * 
   * @param parent Container 
   * @param name   Name of histogram 
   * 
   * @return Pointer to histogram or null 
   */
  static TH1* GetH1(Container* parent, const char* name);
  /** 
   * Get a 2D histogram from a container 
   * 
   * @param parent Container 
   * @param name   Name of histogram 
   * 
   * @return Pointer to histogram or null 
   */
  static TH2* GetH2(Container* parent, const char* name);
  /** 
   * Get a container from a container 
   * 
   * @param parent Container 
   * @param name   Name of container 
   * 
   * @return Pointer to container or null 
   */
  static Container* GetC(Container* parent, const char* name);
  /** 
   * Get a container from a directory 
   * 
   * @param parent Container 
   * @param name   Name of container 
   * 
   * @return Pointer to container or null 
   */
  static Container* GetC(TDirectory* parent, const char* name);
  /** 
   * Get a double-precision value from a container 
   * 
   * @param parent Parent container 
   * @param name   Name of parameter 
   * @param def    Default value if not found 
   * 
   * @return The value (or default if not found)
   */
  static Double_t GetD(Container* parent, const char* name, Double_t def=-1);
  /** 
   * Get a integer value from a container 
   * 
   * @param parent Parent container 
   * @param name   Name of parameter 
   * @param def    Default value if not found 
   * 
   * @return The value (or default if not found)
   */  
  static Int_t GetI(Container* parent, const char* name, Int_t def=-1);
  /** 
   * Get a boolean value from a container 
   * 
   * @param parent Parent container 
   * @param name   Name of parameter 
   * @param def    Default value if not found 
   * 
   * @return The value (or default if not found)
   */  
  static Int_t GetB(Container* parent, const char* name, Bool_t def=false);
  /** 
   * Get a copy of a 1D histogram from a container 
   * 
   * @param parent  Container 
   * @param name    Name of histogram 
   * @param newName Optional new name of copy 
   * 
   * @return Pointer to histogram or null 
   */
  static TH1* CopyH1(Container* parent, const char* name,const char* newName=0);
  /** 
   * Get a copy of a 2D histogram from a container 
   * 
   * @param parent  Container 
   * @param name    Name of histogram 
   * @param newName Optional new name of copy 
   * 
   * @return Pointer to histogram or null 
   */
  static TH2* CopyH2(Container* parent, const char* name,const char* newName=0);
  /** 
   * Service function to make a 1D histogram from an axis definition 
   * 
   * @param name  Name of histogram 
   * @param title Title of histogram 
   * @param xAxis X axis to use 
   * @param c     Container 
   * @param color Color 
   * @param style Marker style 
   * 
   * @return Newly created histogram 
   */
  static TH1* Make1D(Container&     c,
		     const TString& name,
		     const TString& title,
		     Color_t        color,
		     Style_t        style,
		     const TAxis&   xAxis);
  /** 
   * Service function to make a 2D histogram from axis definitions 
   * 
   * @param name   Name of histogram 
   * @param title  Title of histogram 
   * @param xAxis  X axis definition 
   * @param yAxis  Y axis definition 
   * @param c     Container 
   * @param color Color 
   * @param style Marker style 
   * 
   * @return Newly created histogram 
   */
  static TH2* Make2D(Container&     c,
		     const TString& name,
		     const TString& title,
		     Color_t        color,
		     Style_t        style,
		     const TAxis&   xAxis,
		     const TAxis&   yAxis);

  /** 
   * Scale a histogram by factor with error propagation 
   * 
   * @param h   Histogram to scale 
   * @param x   Scalar 
   * @param xe  Error on scalar 
   * 
   * @return Pointer to histogram 
   */
  static TH1* Scale(TH1* h, Double_t x, Double_t xe);
  /** 
   * Scale a histogram by factor with error propagation 
   * 
   * @param h   Histogram to scale 
   * @param x   Scalar 
   * @param xe  Error on scalar 
   * 
   * @return Pointer to histogram 
   */
  static TH2* Scale(TH2* h, Double_t x, Double_t xe);
  /** 
   * Scale the two dimensional histogram with content of the
   * 1-dimensional histogram (over X)
   * 
   * @param h Histogram to scale 
   * @param s The scalar histogram
   * 
   * @return The scaled histogram (pointer to input)
   */
  static TH2* Scale(TH2* h, TH1* s);
  /** 
   * Fix axis attributes
   * 
   * @param axis  Axis to fix 
   * @param title Possible title for axis 
   */
  static void FixAxis(TAxis& axis, const char* title=0);
  /** 
   * Scale bins of an axis by constant factor.  The number of bins
   * remains the same.
   * 
   * @param axis Axis to modify
   * @param fact Factor to scale by
   */
  static void ScaleAxis(TAxis& axis, Double_t fact=1);
  /** 
   * Set an axis based on bin borders 
   * 
   * @param axis    Axis to set 
   * @param n       Number of bins 
   * @param borders Bin borders (n+1 entries)
   */
  static void SetAxis(TAxis& axis, Int_t n, Double_t* borders);
  /** 
   * Set an axis based on test string of specs.  The token separator
   * is given in @a sep.
   * 
   * @param axis Axis to set 
   * @param spec Specification
   * @param sep  Token separate 
   */
  static void SetAxis(TAxis& axis, const TString& spec, const char* sep=":,");
  /** 
   * Set axis with least and largest values
   * 
   * @param axis Axis to set
   * @param n    Number of bins 
   * @param l    Least value 
   * @param h    Largest value 
   */
  static void SetAxis(TAxis& axis, Int_t n, Double_t l, Double_t h);
  /** 
   * Set a symmetric axis 
   * 
   * @param axis Axis to set
   * @param n    Number of bins 
   * @param m    Maximum absolute value 
   */
  static void SetAxis(TAxis& axis, Int_t n, Double_t m);
  /** 
   * Print axis 
   * 
   * @param axis Axis to print 
   * @param nSig Number of significant digits 
   * @param alt  Atlernate name 
   */
  static void PrintAxis(const TAxis& axis, Int_t nSig=2, const char* alt=0);
  /** 
   * Scale each Y-row of input histogram to the number of events in
   * second histogram.  Also scale to bin-width in X direction. 
   * 
   * @param h    Input histogram
   * @param ipZ  Histogram to scale by 
   * @param full If true, do error propegation of IPz errors
   * 
   * @return New copy of input, scaled by second histogram
   */
  static TH2* ScaleToIPz(TH2* h, TH1* ipZ, Bool_t full=false);
  /** 
   * Average primary particle
   * @f$\mathrm{d}N_{\mathrm{ch}}/\mathrm{d}\\eta@f$ over 
   * @f$\mathrm{IP}_z@f$ 
   * 
   * Mode can be one of 
   *
   * - 0: Scale average by content of @a ipz and do error propagation 
   * - 1: Only propagate errors from @a ipz (no scale) 
   * - 2: Do not propagate errors or scale by @a ipz 
   *
   * @param h     Histogram 
   * @param name  Name of projection 
   * @param mode  Mode of operation. 
   * @param ipz   Vertex distribution
   * @param mask  Optional mask - if a bin is zero here, do not count
   *              it in average.
   * 
   * @return Newly allocated histogram or null
   */
  static TH1* AverageOverIPz(TH2*        h,
			     const char* name,
			     UShort_t    mode,
			     TH1*        ipz,
			     TH2*        mask=0);
  /** 
   * Clone an object and add to container 
   * 
   * @param c Container to add to 
   * @param o Object to clone 
   * 
   * @return The new copy of the object 
   */
  static TObject* CloneAndAdd(Container* c, TObject* o);
  /** 
   * Helper function to integrate a histogram.  Note, for symmetric
   * histograms, both sides are integrate, and the errors added in
   * quadrature.
   * 
   * @param h   Histogram to integrate 
   * @param min Least limint 
   * @param max Largest limit 
   * @param err On return, the error on the integral (95% CL)
   * 
   * @return The integral value. 
   */
  static Double_t Integrate(TH1* h, Double_t min, Double_t max, Double_t& err);
  /** 
   * Calculate ratio and error on ratio
   * 
   * @param n   Numerator value 
   * @param en  Numerator error 
   * @param d   Denominator value 
   * @param ed  Denominator error 
   * @param er  On return, ratio error 
   * 
   * @return Ratio 
   */
  static Double_t RatioE(Double_t n, Double_t en,
			 Double_t d, Double_t ed,
			 Double_t& er);
  /** 
   * Calculate ratio and error on ratio
   * 
   * @param n    Numerator value 
   * @param e2n  Squared numerator error 
   * @param d    Denominator value 
   * @param e2d  Squared denominator error 
   * @param e2r  On return, squared ratio error 
   * 
   * @return Ratio 
   */
  static Double_t RatioE2(Double_t n, Double_t e2n,
			  Double_t d, Double_t e2d,
			  Double_t& e2r);
  /* @} */
protected:
  ClassDef(AliTrackletdNdetaUtils,1); // Utilities
};

//====================================================================
// Utilities
//____________________________________________________________________
Bool_t AliTrackletdNdetaUtils::CheckAxisNBins(const char*  which,
					      const TAxis* a1,
					      const TAxis* a2)
{
  if (a1->GetNbins() != a2->GetNbins()) {
    ::Warning("CheckAxisNBins", "Incompatible number %s bins: %d vs %d",
	      which, a1->GetNbins(), a2->GetNbins());
    return false;
  }
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaUtils::CheckAxisLimits(const char*  which,
					       const TAxis* a1,
					       const TAxis* a2)
{
  if (!TMath::AreEqualRel(a1->GetXmin(), a2->GetXmin(),1.E-12) ||
      !TMath::AreEqualRel(a1->GetXmax(), a2->GetXmax(),1.E-12)) {
    Warning("CheckAxisLimits",
	    "Limits of %s axis incompatible [%f,%f] vs [%f,%f]", which,
	    a1->GetXmin(), a1->GetXmax(), a2->GetXmin(), a2->GetXmax());
    return false;
  }
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaUtils::CheckAxisBins(const char*  which,
					     const TAxis* a1,
					     const TAxis* a2)
{
  const TArrayD * h1Array = a1->GetXbins();
  const TArrayD * h2Array = a2->GetXbins();
  Int_t fN = h1Array->fN;
  if ( fN == 0 ) return true;
  if (h2Array->fN != fN) {
    // Redundant?
    Warning("CheckAxisBins", "Not equal number of %s bin limits: %d vs %d",
	    which, fN, h2Array->fN);
    return false;
  }
  else {
    for (int i = 0; i < fN; ++i) {
      if (!TMath::AreEqualRel(h1Array->GetAt(i),h2Array->GetAt(i),1E-10)) {
	Warning("CheckAxisBins",
		"%s limit # %3d incompatible: %f vs %f",
		which, i, h1Array->GetAt(i),h2Array->GetAt(i));
	return false;
      }
    }
  }
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaUtils::CheckAxisLabels(const char*  which,
					       const TAxis* a1,
					       const TAxis* a2)
{
  // check that axis have same labels
  THashList *l1 = (const_cast<TAxis*>(a1))->GetLabels();
  THashList *l2 = (const_cast<TAxis*>(a2))->GetLabels();
  
  if (!l1 && !l2) return true;
  if (!l1 ||  !l2) {
    Warning("CheckAxisLabels", "Difference in %s labels: %p vs %p",
	    which, l1, l2);
    return false;
  }
  // check now labels sizes  are the same
  if (l1->GetSize() != l2->GetSize()) {
    Warning("CheckAxisLabels", "Different number of %s labels: %d vs %d",
	    which, l1->GetSize(), l2->GetSize());
    return false;
  }
  for (int i = 1; i <= a1->GetNbins(); ++i) {
    TString label1 = a1->GetBinLabel(i);
    TString label2 = a2->GetBinLabel(i);
    if (label1 != label2) {
      Warning("CheckAxisLabels", "%s label # %d not the same: '%s' vs '%s'",
	      which, i, label1.Data(), label2.Data());
      return false;
    }
  }
    
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaUtils::CheckAxis(const char*  which, 
					 const TAxis* a1,
					 const TAxis* a2,
					 Bool_t       alsoLbls)
{
  if (!CheckAxisNBins (which, a1, a2)) return false;
  if (!CheckAxisLimits(which, a1, a2)) return false;
  if (!CheckAxisBins  (which, a1, a2)) return false;
  if (alsoLbls && !CheckAxisLabels(which, a1, a2)) return false;
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaUtils::CheckConsistency(const TH1* h1, const TH1* h2)
{
  // Check histogram compatibility
  if (h1 == h2) return true;
    
  if (h1->GetDimension() != h2->GetDimension() ) {
    Warning("CheckConsistency",
	    "%s and %s have different dimensions %d vs %d",
	    h1->GetName(), h2->GetName(),
	    h1->GetDimension(), h2->GetDimension());
    return false;
  }
  Int_t dim = h1->GetDimension(); 
    
  Bool_t alsoLbls = (h1->GetEntries() != 0 && h2->GetEntries() != 0);
  if (!CheckAxis("X", h1->GetXaxis(), h2->GetXaxis(), alsoLbls)) return false;
  if (dim > 1 &&
      !CheckAxis("Y", h1->GetYaxis(), h2->GetYaxis(), alsoLbls)) return false;
  if (dim > 2 &&
      !CheckAxis("Z", h1->GetZaxis(), h2->GetZaxis(), alsoLbls)) return false;
    
  return true;
}

//____________________________________________________________________
TH2* AliTrackletdNdetaUtils::Scale(TH2* h, TH1* s)
{
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    Double_t dc = s->GetBinContent(i);
    Double_t de = s->GetBinError  (i);
    Double_t dr = (dc > 1e-6 ? 1/dc : 0);
    Double_t dq = dr * de;
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) {
      Double_t nc  = h->GetBinContent(i,j);
      Double_t ne  = h->GetBinError  (i,j);
      Double_t ns  = (nc > 0 ? ne/nc : 0);
      Double_t sc  = dr * nc;
      Double_t se  = sc*TMath::Sqrt(ns*ns+dq*dq);
      // Printf("Setting bin %2d,%2d to %f +/- %f", sc, se);
      h->SetBinContent(i,j,sc);
      h->SetBinError  (i,j,se);
    }
  }
  return h;
}

//____________________________________________________________________
TH2* AliTrackletdNdetaUtils::Scale(TH2* h, Double_t x, Double_t xe)
{
  Double_t rr    = xe/x;
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) {
      Double_t c  = h->GetBinContent(i,j);
      Double_t e  = h->GetBinError  (i,j);
      Double_t s  = (c > 0 ? e/c : 0);
      Double_t sc = x * c;
      Double_t se = sc*TMath::Sqrt(s*s+rr*rr);
      // Printf("Setting bin %2d,%2d to %f +/- %f", sc, se);
      h->SetBinContent(i,j,sc);
      h->SetBinError  (i,j,se);
    }
  }
  return h;
}

//____________________________________________________________________
TH1* AliTrackletdNdetaUtils::Scale(TH1* h, Double_t x, Double_t xe)
{
  Double_t rr    = xe/x;
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    Double_t c  = h->GetBinContent(i);
    Double_t e  = h->GetBinError  (i);
    Double_t s  = (c > 0 ? e/c : 0);
    Double_t sc = x * c;
    Double_t se = sc*TMath::Sqrt(s*s+rr*rr);
    // Printf("Setting bin %2d,%2d to %f +/- %f", sc, se);
    h->SetBinContent(i,sc);
    h->SetBinError  (i,se);
  }
  return h;
}

//____________________________________________________________________
TObject* AliTrackletdNdetaUtils::GetO(Container*  parent,
				     const char* name,
				     TClass*     cls)
{
  if (!parent) {
    ::Warning("GetO", "No parent container passed");
    return 0;
  }
  TObject* o = parent->FindObject(name);
  if (!o) {
    ::Warning("GetO", "Object \"%s\" not found in \"%s\"",
	      name, parent->GetName());
    // parent->ls();
    return 0;
  }
  if (!cls) return o;
  if (!o->IsA()->InheritsFrom(cls)) {
    ::Warning("GetO", "\"%s\" is an object of class %s, not %s",
	      name, o->ClassName(), cls->GetName());
    return 0;
  }
  return o;
}
//____________________________________________________________________
TObject* AliTrackletdNdetaUtils::GetO(TDirectory* parent,
				      const char* name,
				      TClass*     cls)
{
  if (!parent) {
    ::Warning("GetO", "No parent directory passed");
    return 0;
  }
  TObject* o = parent->Get(name);
  if (!o) {
    ::Warning("GetO", "Object \"%s\" not found in \"%s\"",
	      name, parent->GetName());
    parent->ls();
    return 0;
  }
  if (!cls) return o;
  if (!o->IsA()->InheritsFrom(cls)) {
    ::Warning("GetO", "\"%s\" is an object of class %s, not %s",
	      name, o->ClassName(), cls->GetName());
    return 0;
  }
  return o;
}
//____________________________________________________________________
TH1* AliTrackletdNdetaUtils::GetH1(Container* parent, const char* name)
{
  return static_cast<TH1*>(GetO(parent, name, TH1::Class()));
}
//____________________________________________________________________
TH2* AliTrackletdNdetaUtils::GetH2(Container* parent, const char* name)
{
  return static_cast<TH2*>(GetO(parent, name, TH2::Class()));
}
//____________________________________________________________________
AliTrackletdNdetaUtils::Container*
AliTrackletdNdetaUtils::GetC(Container* parent, const char* name)
{
  return static_cast<Container*>(GetO(parent, name, Container::Class()));
}
//____________________________________________________________________
AliTrackletdNdetaUtils::Container*
AliTrackletdNdetaUtils::GetC(TDirectory* parent, const char* name)
{
  return static_cast<Container*>(GetO(parent, name, Container::Class()));
}
//____________________________________________________________________
Double_t AliTrackletdNdetaUtils::GetD(Container*  parent,
				     const char* name,
				     Double_t    def)
{
  TParameter<double>* p =
    static_cast<TParameter<double>*>(GetO(parent, name,
					    TParameter<double>::Class()));
  if (!p) return def;
  return p->GetVal();
}
//____________________________________________________________________  
Int_t AliTrackletdNdetaUtils::GetI(Container*  parent,
				  const char* name,
				  Int_t       def)
{
  TParameter<int>* p =
    static_cast<TParameter<int>*>(GetO(parent, name, TParameter<int>::Class()));
  if (!p) return def;
  return p->GetVal();
}
//____________________________________________________________________  
Int_t AliTrackletdNdetaUtils::GetB(Container*  parent,
				  const char* name,
				  Bool_t      def)
{
  TParameter<bool>* p =
    static_cast<TParameter<bool>*>(GetO(parent, name,
					  TParameter<bool>::Class()));
  if (!p) return def;
  return p->GetVal();
}
//____________________________________________________________________
TH1* AliTrackletdNdetaUtils::CopyH1(Container*  parent,
				    const char* name,
				    const char* newName)
{
  TH1* orig = GetH1(parent, name);
  if (!orig) return 0;
  TH1* ret  = static_cast<TH1*>(orig->Clone(newName ? newName : name));
  ret->SetDirectory(0); // Release from file container
  return ret;
}
//____________________________________________________________________
TH2* AliTrackletdNdetaUtils::CopyH2(Container*  parent,
				    const char* name,
				    const char* newName)
{
  TH2* orig = GetH2(parent, name);
  if (!orig) return 0;
  TH2* ret  = static_cast<TH2*>(orig->Clone(newName ? newName : name));
  ret->SetDirectory(0); // Release from file container
  return ret;
}
//____________________________________________________________________
TH1* AliTrackletdNdetaUtils::Make1D(Container&     c,
				   const TString& name,
				   const TString& title,
				   Color_t        color,
				   Style_t        style,
				   const TAxis&   xAxis)
{
  TString   n   = name;
  TString   t   = title;
  TH1*      ret = 0;
  Int_t     nx  = xAxis.GetNbins();
  if (t.IsNull()) t = xAxis.GetTitle();
  if (xAxis.GetXbins() && xAxis.GetXbins()->GetArray())
    ret = new TH1D(n,t,nx,xAxis.GetXbins()->GetArray());
  else
    ret = new TH1D(n,t,nx,xAxis.GetXmin(),xAxis.GetXmax());
  ret->Sumw2();
  ret->SetXTitle(xAxis.GetTitle());
  static_cast<const TAttAxis&>(xAxis).Copy(*(ret->GetXaxis()));
  ret->SetDirectory(0);
  ret->SetLineColor(color);
  ret->SetMarkerColor(color);
  ret->SetFillColor(kWhite);// color);
  ret->SetFillStyle(0);
  ret->SetMarkerStyle(style);
  if (const_cast<TAxis&>(xAxis).GetLabels()) {
    for (Int_t i = 1; i <= xAxis.GetNbins(); i++)
      ret->GetXaxis()->SetBinLabel(i, xAxis.GetBinLabel(i));
  }
  switch (style) {
  case 27:
  case 33:
    ret->SetMarkerSize(1.4);
    break;
  case 29:    
  case 30:
    ret->SetMarkerSize(1.2);
    break;
  }
  c.Add(ret);
  return ret;
}
//____________________________________________________________________
TH2* AliTrackletdNdetaUtils::Make2D(Container&     c,
				   const TString& name,
				   const TString& title,
				   Color_t        color,
				   Style_t        style,
				   const TAxis&   xAxis,
				   const TAxis&   yAxis)
{
  TString   n         = name;
  TString   t         = title;
  TH2*      ret       = 0;
  Int_t     nx        = xAxis.GetNbins();
  Int_t     ny        = yAxis.GetNbins();
  const Double_t* xb  = (xAxis.GetXbins() && xAxis.GetXbins()->GetArray() ?
			 xAxis.GetXbins()->GetArray() : 0);
  const Double_t* yb  = (yAxis.GetXbins() && yAxis.GetXbins()->GetArray() ?
			 yAxis.GetXbins()->GetArray() : 0);
  if (t.IsNull())
    t.Form("%s vs %s", yAxis.GetTitle(), xAxis.GetTitle());
  if (xb) {	  
    if   (yb) ret = new TH2D(n,t,nx,xb,ny,yb);
    else      ret = new TH2D(n,t,
			     nx,xb,
			     ny,yAxis.GetXmin(),yAxis.GetXmax());
  }
  else { 
    if   (yb) ret = new TH2D(n,t,
			     nx,xAxis.GetXmin(),xAxis.GetXmax(),
			     ny,yb);
    else      ret = new TH2D(n,t,
			     nx,xAxis.GetXmin(),xAxis.GetXmax(),
			     ny,yAxis.GetXmin(),yAxis.GetXmax());
  }
  ret->Sumw2();
  ret->SetXTitle(xAxis.GetTitle());
  ret->SetYTitle(yAxis.GetTitle());
  ret->SetLineColor(color);
  ret->SetMarkerColor(color);
  ret->SetFillColor(color);
  ret->SetMarkerStyle(style);
  static_cast<const TAttAxis&>(xAxis).Copy(*(ret->GetXaxis()));
  static_cast<const TAttAxis&>(yAxis).Copy(*(ret->GetYaxis()));
  ret->SetDirectory(0);
  if (const_cast<TAxis&>(xAxis).GetLabels()) {
    for (Int_t i = 1; i <= xAxis.GetNbins(); i++)
      ret->GetXaxis()->SetBinLabel(i, xAxis.GetBinLabel(i));
  }
  if (const_cast<TAxis&>(yAxis).GetLabels()) {
    for (Int_t i = 1; i <= yAxis.GetNbins(); i++)
      ret->GetYaxis()->SetBinLabel(i, yAxis.GetBinLabel(i));
  }
  c.Add(ret);
  return ret;
}
//____________________________________________________________________
void AliTrackletdNdetaUtils::FixAxis(TAxis& axis, const char* title)
{
  if (title && title[0] != '\0') axis.SetTitle(title);
  axis. SetNdivisions(210);
  axis. SetLabelFont(42);
  axis. SetLabelSize(0.03);
  axis. SetLabelOffset(0.005);
  axis. SetLabelColor(kBlack);
  axis. SetTitleOffset(1);
  axis. SetTitleFont(42);
  axis. SetTitleSize(0.04);
  axis. SetTitleColor(kBlack);
  axis. SetTickLength(0.03);
  axis. SetAxisColor(kBlack);
}
//____________________________________________________________________
void AliTrackletdNdetaUtils::ScaleAxis(TAxis&       ret,
				      Double_t     fact)
{
  if (ret.GetXbins()->GetArray()) {
    TArrayD bins(*ret.GetXbins());
    for (Int_t i = 0; i < bins.GetSize(); i++) bins[i] *= fact;
    ret.Set(ret.GetNbins(), bins.GetArray());
  }
  else {
    ret.Set(ret.GetNbins(), fact*ret.GetXmin(), fact*ret.GetXmax());
  }
  FixAxis(ret);
}
//____________________________________________________________________
void AliTrackletdNdetaUtils::SetAxis(TAxis& axis, Int_t n, Double_t* borders)
{
  axis.Set(n, borders);
  FixAxis(axis);
}
//____________________________________________________________________
void AliTrackletdNdetaUtils::SetAxis(TAxis&         axis,
				    const TString& spec,
				    const char*    sep)
{
  TString s(spec);
  Bool_t isRange = false, isUnit = false;
  if (s[0] == 'r' || s[0] == 'R') {
    isRange = true;
    s.Remove(0,1);
  }
  if (s[0] == 'u' || s[0] == 'U') {
    isUnit = true;
    s.Remove(0,1);
  }
  TObjArray*  tokens = s.Tokenize(sep);
  TArrayD     bins(tokens->GetEntries());
  TObjString* token = 0;
  TIter       next(tokens);
  Int_t       i = 0;
  while ((token = static_cast<TObjString*>(next()))) {
    Double_t v = token->String().Atof();
    bins[i] = v;
    i++;
  }
  delete tokens;
  if (isUnit) {
    if (bins.GetSize() > 1)
      SetAxis(axis, Int_t(bins[1]-bins[0]), bins[0], bins[1]);
    else
      SetAxis(axis, 2*Int_t(bins[0]), bins[0]);
  }      
  else if (isRange) {
    Int_t    nBins = Int_t(bins[0]);
    if (bins.GetSize() > 2) 
      SetAxis(axis, nBins, bins[1], bins[2]);
    else
      SetAxis(axis, nBins, bins[1]);
  }
  else 
    SetAxis(axis, bins.GetSize()-1,bins.GetArray());
}
//____________________________________________________________________
void AliTrackletdNdetaUtils::SetAxis(TAxis&   axis,
				    Int_t    n,
				    Double_t l,
				    Double_t h)
{
  axis.Set(n, l, h);
  FixAxis(axis);
}
//____________________________________________________________________
void AliTrackletdNdetaUtils::SetAxis(TAxis& axis, Int_t n, Double_t m)
{
  SetAxis(axis, n, -TMath::Abs(m), +TMath::Abs(m));
}
//____________________________________________________________________
void AliTrackletdNdetaUtils::PrintAxis(const TAxis& axis,
				      Int_t nSig,
				      const char* alt)
{
  printf(" %17s axis: ", alt ? alt : axis.GetTitle());
  if (axis.GetXbins() && axis.GetXbins()->GetArray()) {
    printf("%.*f", nSig, axis.GetBinLowEdge(1));
    for (Int_t i = 1; i <= axis.GetNbins(); i++)
      printf(":%.*f", nSig, axis.GetBinUpEdge(i));
  }
  else
    printf("%d bins between %.*f and %.*f",
	   axis.GetNbins(), nSig, axis.GetXmin(),nSig,axis.GetXmax());
  printf("\n");
}
//____________________________________________________________________
TH2* AliTrackletdNdetaUtils::ScaleToIPz(TH2* h, TH1* ipZ, Bool_t full)
{
  if (!h) {
    ::Warning("ScaleToIPz","Nothing to scale");
    return 0;
  }
  if (!ipZ) {
    ::Warning("ScaleToIPz","Nothing to scale by");
    return 0;
  }
  TH2* ret = static_cast<TH2*>(h->Clone());
  ret->SetDirectory(0);
  if (!ipZ) return ret;
  for (Int_t iy = 1; iy <= ret->GetNbinsY(); iy++) {
    Double_t z   = ret->GetYaxis()->GetBinCenter(iy);
    Int_t    bin = ipZ->GetXaxis()->FindBin(z);
    Double_t nEv = ipZ->GetBinContent(bin);
    Double_t eEv = ipZ->GetBinError  (bin);
    Double_t esc = (nEv > 0 ? 1./nEv : 0);
    Double_t rE2 = esc*esc*eEv*eEv;
    for (Int_t ix = 1; ix <= ret->GetNbinsX(); ix++) {
      Double_t  c   = ret->GetBinContent(ix,iy);
      Double_t  e   = ret->GetBinError  (ix,iy);
      Double_t  r   = (c > 0 ? e/c : 0);
      // Scale by number of events, and error propagate 
      Double_t  sc  = c * esc;
      Double_t  se  = 0;
      if (full) se  = sc * TMath::Sqrt(r*r+rE2);
      else      se  = e * esc;
      Double_t scl = 1 / ret->GetXaxis()->GetBinWidth(ix);
      ret->SetBinContent(ix, iy, scl*sc);
      ret->SetBinError  (ix, iy, scl*se);
    }
  }
  return ret;
}
//____________________________________________________________________
TH1* AliTrackletdNdetaUtils::AverageOverIPz(TH2*        h,
					    const char* name,
					    UShort_t    mode,
					    TH1*        ipz,
					    TH2*        other)
{
  if (!h) return 0;
  
  Int_t nIPz = h->GetNbinsY();
  Int_t nEta = h->GetNbinsX();
  TH1*  p    = h->ProjectionX(name,1,nIPz,"e");
  TH2*  mask = (other ? other : h);
  p->SetDirectory(0);
  p->SetFillColor(0);
  p->SetFillStyle(0);
  p->SetYTitle(Form("#LT%s#GT", h->GetYaxis()->GetTitle()));
  p->Reset();
  for (Int_t etaBin = 1; etaBin <= nEta; etaBin++) {
    TArrayD hv(nIPz);
    TArrayD he(nIPz);
    TArrayD hr(nIPz);
    TArrayI hb(nIPz);
    Int_t   j = 0;
    for (Int_t ipzBin = 1; ipzBin <= nIPz; ipzBin++) {
      Double_t bc = mask->GetBinContent(etaBin, ipzBin);
      if (bc < 1e-9) continue; // Low value
      Double_t be = mask->GetBinError  (etaBin, ipzBin);
      if (TMath::IsNaN(be)) continue; // Bad error value 
      hv[j] = bc;
      he[j] = be;
      hr[j] = be/bc;
      hb[j] = ipzBin;
      j++;		
    }
    // Now we have all interesting values.  Sort the relative error
    // array to get the most significant first.  Note, we sort on the
    // mask errors which may not be the same as the histogram errors.
    TArrayI idx(nIPz);
    TMath::Sort(j, hr.fArray, idx.fArray, false);
    Double_t nsum  = 0; // Running weighted sum
    Double_t nsume = 0; // Running weighted sum error
    Double_t dsum  = 0;
    Double_t dsume = 0;
    Int_t    n     = 0;
    Double_t rat   = 1e99;
    
    Int_t k = 0;
    for (k = 0; k < j; k++) {
      Int_t    l      = idx[k]; // Sorted index - ipzBin
      Int_t    ipzBin = hb[l];
      Double_t hvv    = hv[l];      
      Double_t hee    = he[l];
      Double_t x      = TMath::Sqrt(nsume+hee*hee)/(nsum+hvv);
      if (x > rat) {
	continue; // Ignore - does not help
      }

      Double_t by = mask->GetYaxis()->GetBinCenter(ipzBin);
      Int_t    ib = ipz ? ipz->FindBin(by) : 0;
      rat   =  x;
      nsum  += h->GetBinContent(etaBin, ipzBin);
      nsume += TMath::Power(h->GetBinError(etaBin, ipzBin), 2);
      // If we do not have the vertex distribution, then just count
      // number of observations. 
      dsum  += !ipz ? 1 : ipz->GetBinContent(ib);
      dsume += !ipz ? 0 : TMath::Power(ipz->GetBinError(ib),2);
      n++;
    }
    if (k == 0 || n == 0) {
      ::Warning("Average", "Eta bin # %3d has no data",etaBin);
      continue; // This eta bin empty!
    }
    Double_t norm = (mode > 0 ? n : dsum);
    Double_t rne  = nsume/nsum/nsum;
    Double_t rde  = dsume/dsum/dsum;
    Double_t av   = nsum/norm;
    Double_t ave  = 0;
    if (mode==2) ave = TMath::Sqrt(nsume)/n;
    else         ave = av*TMath::Sqrt(rne+rde);
    p->SetBinContent(etaBin, av);
    p->SetBinError  (etaBin, ave);
  }
  if (mode == 0) p->Scale(1, "width");
  return p;		   
}

//____________________________________________________________________
TObject* AliTrackletdNdetaUtils::CloneAndAdd(Container* c, TObject* o)
{
  if (!o) {
    ::Warning("CloneAndAdd", "Nothing to clone");
    return 0;
  }
  TObject* copy = o->Clone();
  if (copy->IsA()->InheritsFrom(TH1::Class()))
    // Release from underlying directory 
    static_cast<TH1*>(copy)->SetDirectory(0);
  if (c)
    // Add to output container
    c->Add(copy);
  return copy;
}
//____________________________________________________________________
Double_t AliTrackletdNdetaUtils::Integrate(TH1*     h,
					  Double_t min,
					  Double_t max,
					  Double_t& err)
{
  const Double_t epsilon = 1e-6;
  Int_t bMin = h->FindBin(min+epsilon);
  Int_t bMax = h->FindBin(max-epsilon);
  if (bMin < 1) bMin = 0; // Include underflow bin
  Double_t val = h->IntegralAndError(bMin, bMax, err);
  // For a-symmetric histograms, return 
  if (TMath::Abs(h->GetXaxis()->GetXmin()+h->GetXaxis()->GetXmax())>=epsilon)
    return val;

  // Otherwise, also integrate negative side
  Double_t err2;
  bMin =  h->FindBin(-min+epsilon);
  bMax =  h->FindBin(-max-epsilon);
  val  += h->IntegralAndError(bMin, bMax, err2);
  err  =  TMath::Sqrt(err*err+err2*err2);
  return val;
}
//____________________________________________________________________
Double_t AliTrackletdNdetaUtils::RatioE(Double_t n, Double_t en,
				       Double_t d, Double_t ed,
				       Double_t& er)
{
  Double_t r = 0;
  er = 0;
  if (TMath::Abs(n) < 1e-16 || TMath::Abs(d) < 1e-16) return 0;
  r  = n/d;
  er = TMath::Sqrt(en*en/n/n + ed*ed/d/d);
  return r;
}
//____________________________________________________________________
Double_t AliTrackletdNdetaUtils::RatioE2(Double_t n, Double_t e2n,
					Double_t d, Double_t e2d,
					Double_t& e2r)
{
  Double_t r = 0;
  e2r = 0;
  if (TMath::Abs(n) < 1e-16 || TMath::Abs(d) < 1e-16) return 0;
  r   = n/d;
  e2r = (e2n/n/n + e2d/d/d);
  return r;
}

#endif
// Local Variables:
//  mode: C++
// End:

