/**
 * @file   AliTrackletAODUtils.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 27 16:46:28 2016
 * 
 * @brief  Utilities for midrapidity analysis 
 * 
 * This contains code shared at different stages in the analysis 
 *
 * @ingroup pwglf_forward_tracklets
 */
#ifndef ALITRACKLETAODUTILS_H
#define ALITRACKLETAODUTILS_H
#ifndef __CINT__
# include <TH1.h>
# include <TH2.h>
# include <TH3.h>
# include <TList.h>
# include <TParameter.h>
# include <TError.h>
# include <TMath.h>
# include <TDirectory.h>
# include <THashList.h>
# include <TProfile.h>
# include <TProfile2D.h>
# include <TDatabasePDG.h>
#else
class TList;
class TH1;
class TH2;
class TH3;
class TProfile;
class TProfile2D;
class TAxis;
class TDirectory;
class TDatabasePDG; // Auto-load
#endif

/**
 * @defgroup pwglf_forward_tracklets  Mid-rapidity tracklet code for dN/deta
 * 
 */
/**
 * Class with utlity functions 
 * 
 * @ingroup pwglf_forward_tracklets
 */
class AliTrackletAODUtils
{
public:
  typedef TList Container;
  /** 
   * Dummy constructor 
   */
  AliTrackletAODUtils() {}
  /**
   * Dummy destructor 
   */
  virtual ~AliTrackletAODUtils() {}
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
   * @param verb   Whether to be verboe 
   * 
   * @return Pointer to object or null
   */
  static TObject* GetO(Container*  parent,
		       const char* name,
		       TClass*     cls=0,
		       Bool_t      verb=true);
  /** 
   * Get an object from a directory. 
   * 
   * @param parent directory 
   * @param name   Name of object 
   * @param cls    If not null, check if object has this type 
   * @param verb   Whether to be verboe 
   * 
   * @return Pointer to object or null
   */
  static TObject* GetO(TDirectory* parent,
		       const char* name,
		       TClass*     cls=0,
		       Bool_t      verb=true);
  /** 
   * Get a 1D histogram from a container 
   * 
   * @param parent Container 
   * @param name   Name of histogram 
   * @param verb   Whether to be verbose
   * 
   * @return Pointer to histogram or null 
   */
  static TH1* GetH1(Container* parent, const char* name, Bool_t verb=true);
  /** 
   * Get a 1D profile from a container 
   * 
   * @param parent Container 
   * @param name   Name of profile 
   * @param verb   Whether to be verbose
   * 
   * @return Pointer to profile or null 
   */
  static TProfile* GetP1(Container* parent, const char* name, Bool_t verb=true);
  /** 
   * Get a 2D histogram from a container 
   * 
   * @param parent Container 
   * @param name   Name of histogram 
   * @param verb   Whether to be verbose
   * 
   * @return Pointer to histogram or null 
   */
  static TH2* GetH2(Container* parent, const char* name, Bool_t verb=true);
  /** 
   * Get a 3D histogram from a container 
   * 
   * @param parent Container 
   * @param name   Name of histogram 
   * @param verb   Whether to be verbose
   * 
   * @return Pointer to histogram or null 
   */
  static TH3* GetH3(Container* parent, const char* name, Bool_t verb=true);
  /** 
   * Get a 2D profile from a container 
   * 
   * @param parent Container 
   * @param name   Name of profile 
   * @param verb   Whether to be verbose
   * 
   * @return Pointer to profile or null 
   */
  static TProfile2D* GetP2(Container* parent, const char* name,
			   Bool_t verb=true);
  /** 
   * Get a profile from a container 
   * 
   * @param parent Container 
   * @param name   Name of profile 
   * @param verb   Whether to be verbose
   * 
   * @return Pointer to profile or null 
   */
  static TProfile* GetP(Container* parent, const char* name, Bool_t verb=true)
  {
    return GetP1(parent, name, verb);
  }
  /** 
   * Get a container from a container 
   * 
   * @param parent Container 
   * @param name   Name of container 
   * @param verb   Whether to be verbose
   * 
   * @return Pointer to container or null 
   */
  static Container* GetC(Container* parent, const char* name, Bool_t verb=true);
  /** 
   * Get a container from a directory 
   * 
   * @param parent Container 
   * @param name   Name of container 
   * @param verb   Whether to be verbose
   * 
   * @return Pointer to container or null 
   */
  static Container* GetC(TDirectory* parent, const char* name,Bool_t verb=true);
  /** 
   * Get a double-precision value from a container 
   * 
   * @param parent Parent container 
   * @param name   Name of parameter 
   * @param def    Default value if not found 
   * @param verb   Whether to be verbose
   * 
   * @return The value (or default if not found)
   */
  static Double_t GetD(Container* parent, const char* name,
		       Double_t def=-1, Bool_t verb=true);
  /** 
   * Get a integer value from a container 
   * 
   * @param parent Parent container 
   * @param name   Name of parameter 
   * @param def    Default value if not found 
   * @param verb   Whether to be verbose
   * 
   * @return The value (or default if not found)
   */  
  static Int_t GetI(Container* parent, const char* name,
		    Int_t def=-1, Bool_t verb=true);
  /** 
   * Get a boolean value from a container 
   * 
   * @param parent Parent container 
   * @param name   Name of parameter 
   * @param def    Default value if not found 
   * @param verb   Whether to be verbose
   * 
   * @return The value (or default if not found)
   */  
  static Int_t GetB(Container* parent, const char* name,
		    Bool_t def=false, Bool_t verb=true);
  /** 
   * Get a copy of a 1D histogram from a container 
   * 
   * @param parent  Container 
   * @param name    Name of histogram 
   * @param newName Optional new name of copy 
   * @param verb    Whether to be verbose
   * 
   * @return Pointer to histogram or null 
   */
  static TH1* CopyH1(Container*  parent,
		     const char* name,
		     const char* newName=0,
		     Bool_t      verb=true);
  /** 
   * Get a copy of a 2D histogram from a container 
   * 
   * @param parent  Container 
   * @param name    Name of histogram 
   * @param newName Optional new name of copy 
   * @param verb    Whether to be verbose
   * 
   * @return Pointer to histogram or null 
   */
  static TH2* CopyH2(Container*  parent,
		     const char* name,
		     const char* newName=0,
		     Bool_t      verb=true);
  /** 
   * Get a copy of a 3D histogram from a container 
   * 
   * @param parent  Container 
   * @param name    Name of histogram 
   * @param newName Optional new name of copy 
   * @param verb    Whether to be verbose
   * 
   * @return Pointer to histogram or null 
   */
  static TH3* CopyH3(Container*  parent,
		     const char* name,
		     const char* newName=0,
		     Bool_t      verb=true);
  /** 
   * Copy attributes from one histogram to another 
   * 
   * @param src  Source histogram 
   * @param dest Destination histogram
   */
  static void CopyAttr(const TH1* src, TH1* dest);
  /** 
   * Service function to make a 1D histogram from an axis definition 
   * 
   * @param c     Container to add to 
   * @param name  Name of histogram 
   * @param title Title of histogram 
   * @param color Color 
   * @param style Marker style 
   * @param xAxis X axis to use 
   * 
   * @return Newly created histogram 
   */
  static TH1* Make1D(Container*     c,
		     const TString& name,
		     const TString& title,
		     Color_t        color,
		     Style_t        style,
		     const TAxis&   xAxis);
  /** 
   * Service function to make a 1D profile from an axis definition 
   * 
   * @param c     Container to add to 
   * @param name  Name of profile 
   * @param title Title of profile 
   * @param color Color 
   * @param style Marker style 
   * @param xAxis X axis to use 
   * 
   * @return Newly created profile
   */
  static TProfile* Make1P(Container*     c,
			  const TString& name,
			  const TString& title,
			  Color_t        color,
			  Style_t        style,
			  const TAxis&   xAxis);
  /** 
   * Service function to make a 2D histogram from axis definitions 
   * 
   * @param c     Container to add to 
   * @param name  Name of histogram 
   * @param title Title of histogram 
   * @param color Color 
   * @param style Marker style 
   * @param xAxis X axis definition 
   * @param yAxis Y axis definition 
   * 
   * @return Newly created histogram 
   */
  static TH2* Make2D(Container*     c,
		     const TString& name,
		     const TString& title,
		     Color_t        color,
		     Style_t        style,
		     const TAxis&   xAxis,
		     const TAxis&   yAxis);
  /** 
   * Service function to make a 3D histogram from axis definitions 
   * 
   * @param c     Container to add to 
   * @param name  Name of histogram 
   * @param title Title of histogram 
   * @param color Color 
   * @param style Marker style 
   * @param xAxis X axis definition 
   * @param yAxis Y axis definition 
   * @param zAxis Z axis definition 
   * 
   * @return Newly created histogram 
   */
  static TH3* Make3D(Container*     c,
		     const TString& name,
		     const TString& title,
		     Color_t        color,
		     Style_t        style,
		     const TAxis&   xAxis,
		     const TAxis&   yAxis,
		     const TAxis&   zAxis);
  /** 
   * Service function to make a 2D profile from axis definitions 
   * 
   * @param c     Container to add to 
   * @param name  Name of profile 
   * @param title Title of profile 
   * @param color Color 
   * @param style Marker style 
   * @param xAxis X axis definition 
   * @param yAxis Y axis definition 
   * 
   * @return Newly created profile 
   */
  static TProfile2D* Make2P(Container*     c,
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
   * @param ret  Axis to modify
   * @param fact Factor to scale by
   */
  static void ScaleAxis(TAxis& ret, Double_t fact=1);
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
   * @param alt  Alternative 
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
  static TH3* ScaleToIPz(TH3* h, TH1* ipZ, Bool_t full=false);
  static TH3* ScaleDelta(TH3* h, TH2* scale);
  /** 
   * Project (integrate) 
   @f[
   \frac{\mathrm{d}^3N}{\mathrm{d}\eta\mathrm{d}\Delta\mathrm{dIP}_z}
   @f] 
   * over @f$\mathrm{IP}_z@f$ to return 
   *
   @f[ 
   \frac{\mathrm{d}^2N}{\mathrm{d}\eta\mathrm{d}\Delta}
   @f] 
   * 
   * The integral is scaled to number of bins, that is, we calculate
   * the weighted average over @f$\mathrm{IP}_z@f$ 
   * 
   @f{eqnarray}
   N(\eta,\Delta) &=& \sum_{\mathrm{IP}_z} w N(\eta,\Delta,\mathrm{IP}_z)\\
   \delta N(\eta,\Delta) &=& \sqrt{1/\sum_{\mathrm{IP}_z} w}\\
   w &=&  \delta N(\eta,\Delta,\mathrm{IP}_z)
   @f{eqnarray}
   * 
   * @param h 3D histogram to integrate 
   * 
   * @return Newly allocated projection 
   */
  static TH2* ProjectEtaDelta(TH3* h);
  /** 
   * Project (integrate) 
   @f[ 
   \frac{\mathrm{d}^2N}{\mathrm{d}\eta\mathrm{d}\Delta}
   @f] 
   * over @f$\eta@f$ 
   *
   @f[ 
   \frac{\mathrm{d}N}{\mathrm{d}\Delta}
   @f] 
   * 
   * @param h 
   * 
   * @return 
   */
  static TH1* ProjectDelta(TH2* h);
  /** 
   * Project (integrate) 
   @f[ 
   \frac{\mathrm{d}^3N}{\mathrm{d}\eta\mathrm{d}\Delta\mathrm{dIP}_z}
   @f] 
   * over @f$\eta,\mathrm{IP}_z@f$ 
   *
   @f[ 
   \frac{\mathrm{d}N}{\mathrm{d}\Delta}
   @f] 
   * 
   * @param h 
   * 
   * @return 
   */
  static TH1* ProjectDeltaFull(TH3* h);
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
   * @param verb  Whether to be verbose 
   * 
   * @return Newly allocated histogram or null
   */
  static TH1* AverageOverIPz(TH2*        h,
			     const char* name,
			     UShort_t    mode,
			     TH1*        ipz,
			     TH2*        mask=0,
			     Bool_t      verb=true);
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
  /** 
   * Get attributes corresponding to a PDG code 
   * 
   * @param pdg PDG code 
   * @param nme On return, the name 
   * @param c   On return, the color 
   * @param s   On return, the marker style 
   */
  static void PdgAttr(Int_t pdg, TString& nme, Color_t& c, Style_t& s);
  /** 
   * @{ 
   * @name Particle type utilities 
   */
  static Int_t* PdgArray(Int_t& size);    
  /** 
   * Get the PDG axis 
   * 
   * @return Reference to static TAxis object
   */
  static const TAxis& PdgAxis();
  /** 
   * Get the PDG bin 
   * 
   * @param pdg Particle type 
   * 
   * @return Bin number 
   */
  static Int_t  PdgBin(Int_t pdg);
  /* @} */
protected:
  ClassDef(AliTrackletAODUtils,1); // Utilities
};

//====================================================================
// Utilities
//____________________________________________________________________
Bool_t AliTrackletAODUtils::CheckAxisNBins(const char*  which,
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
Bool_t AliTrackletAODUtils::CheckAxisLimits(const char*  which,
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
Bool_t AliTrackletAODUtils::CheckAxisBins(const char*  which,
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
Bool_t AliTrackletAODUtils::CheckAxisLabels(const char*  which,
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
Bool_t AliTrackletAODUtils::CheckAxis(const char*  which, 
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
Bool_t AliTrackletAODUtils::CheckConsistency(const TH1* h1, const TH1* h2)
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
TH2* AliTrackletAODUtils::Scale(TH2* h, TH1* s)
{
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    Double_t dc = s->GetBinContent(i);
    Double_t de = s->GetBinError  (i);
    Double_t dr = (dc > 1e-6 ? de/dc : 0);
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) {
      Double_t nc  = h->GetBinContent(i,j);
      Double_t ne  = h->GetBinError  (i,j);
      Double_t ns  = (nc > 0 ? ne/nc : 0);
      Double_t sc  = dc * nc;
      Double_t se  = sc*TMath::Sqrt(ns*ns+dr*dr);
      // Printf("Setting bin %2d,%2d to %f +/- %f", sc, se);
      h->SetBinContent(i,j,sc);
      h->SetBinError  (i,j,se);
    }
  }
  return h;
}

//____________________________________________________________________
TH2* AliTrackletAODUtils::Scale(TH2* h, Double_t x, Double_t xe)
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
TH1* AliTrackletAODUtils::Scale(TH1* h, Double_t x, Double_t xe)
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
TObject* AliTrackletAODUtils::GetO(Container*  parent,
				   const char* name,
				   TClass*     cls,
				   Bool_t      verb)
{
  if (!parent) {
    if (verb) ::Warning("GetO", "No parent container passed");
    return 0;
  }
  TObject* o = parent->FindObject(name);
  if (!o) {
    if (verb) ::Warning("GetO", "Object \"%s\" not found in \"%s\"",
			name, parent->GetName());
    // parent->ls();
    return 0;
  }
  if (!cls) return o;
  if (!o->IsA()->InheritsFrom(cls)) {
    if (verb) ::Warning("GetO", "\"%s\" is an object of class %s, not %s",
			name, o->ClassName(), cls->GetName());
    return 0;
  }
  return o;
}
//____________________________________________________________________
TObject* AliTrackletAODUtils::GetO(TDirectory* parent,
				   const char* name,
				   TClass*     cls,
				   Bool_t      verb)
{
  if (!parent) {
    if (!verb) ::Warning("GetO", "No parent directory passed");
    return 0;
  }
  TObject* o = parent->Get(name);
  if (!o) {
    if (verb) ::Warning("GetO", "Object \"%s\" not found in \"%s\"",
			name, parent->GetName());
    parent->ls();
    return 0;
  }
  if (!cls) return o;
  if (!o->IsA()->InheritsFrom(cls)) {
    if (verb) ::Warning("GetO", "\"%s\" is an object of class %s, not %s",
			name, o->ClassName(), cls->GetName());
    return 0;
  }
  return o;
}
//____________________________________________________________________
TH1* AliTrackletAODUtils::GetH1(Container* parent, const char* name, Bool_t v)
{
  return static_cast<TH1*>(GetO(parent, name, TH1::Class(), v));
}
//____________________________________________________________________
TH2* AliTrackletAODUtils::GetH2(Container* parent, const char* name, Bool_t v)
{
  return static_cast<TH2*>(GetO(parent, name, TH2::Class(), v));
}
//____________________________________________________________________
TH3* AliTrackletAODUtils::GetH3(Container* parent, const char* name, Bool_t v)
{
  return static_cast<TH3*>(GetO(parent, name, TH3::Class(), v));
}
//____________________________________________________________________
TProfile* AliTrackletAODUtils::GetP1(Container* parent,const char* name,
				     Bool_t v)
{
  return static_cast<TProfile*>(GetO(parent, name, TProfile::Class(), v));
}
//____________________________________________________________________
TProfile2D* AliTrackletAODUtils::GetP2(Container* parent, const char* name,
				       Bool_t v)
{
  return static_cast<TProfile2D*>(GetO(parent, name, TProfile2D::Class(), v));
}
//____________________________________________________________________
AliTrackletAODUtils::Container*
AliTrackletAODUtils::GetC(Container* parent, const char* name, Bool_t v)
{
  return static_cast<Container*>(GetO(parent, name, Container::Class(), v));
}
//____________________________________________________________________
AliTrackletAODUtils::Container*
AliTrackletAODUtils::GetC(TDirectory* parent, const char* name, Bool_t v)
{
  return static_cast<Container*>(GetO(parent, name, Container::Class(), v));
}
//____________________________________________________________________
Double_t AliTrackletAODUtils::GetD(Container*  parent,
				   const char* name,
				   Double_t    def,
				   Bool_t      v)
{
  TParameter<double>* p =
    static_cast<TParameter<double>*>(GetO(parent, name,
					  TParameter<double>::Class(),v));
  if (!p) return def;
  return p->GetVal();
}
//____________________________________________________________________  
Int_t AliTrackletAODUtils::GetI(Container*  parent,
				const char* name,
				Int_t       def,
				Bool_t      v)
{
  TParameter<int>* p =
    static_cast<TParameter<int>*>(GetO(parent, name,
				       TParameter<int>::Class(),v));
  if (!p) return def;
  return p->GetVal();
}
//____________________________________________________________________  
Int_t AliTrackletAODUtils::GetB(Container*  parent,
				const char* name,
				Bool_t      def,
				Bool_t      v)
{
  TParameter<bool>* p =
    static_cast<TParameter<bool>*>(GetO(parent, name,
					TParameter<bool>::Class(), v));
  if (!p) return def;
  return p->GetVal();
}
//____________________________________________________________________
TH1* AliTrackletAODUtils::CopyH1(Container*  parent,
				 const char* name,
				 const char* newName,
				 Bool_t      v)
{
  TH1* orig = GetH1(parent, name, v);
  if (!orig) return 0;
  TH1* ret  = static_cast<TH1*>(orig->Clone(newName ? newName : name));
  ret->SetDirectory(0); // Release from file container
  return ret;
}
//____________________________________________________________________
TH2* AliTrackletAODUtils::CopyH2(Container*  parent,
				 const char* name,
				 const char* newName,
				 Bool_t      v)
{
  TH2* orig = GetH2(parent, name, v);
  if (!orig) return 0;
  TH2* ret  = static_cast<TH2*>(orig->Clone(newName ? newName : name));
  ret->SetDirectory(0); // Release from file container
  return ret;
}
//____________________________________________________________________
TH3* AliTrackletAODUtils::CopyH3(Container*  parent,
				 const char* name,
				 const char* newName,
				 Bool_t      v)
{
  TH3* orig = GetH3(parent, name, v);
  if (!orig) return 0;
  TH3* ret  = static_cast<TH3*>(orig->Clone(newName ? newName : name));
  ret->SetDirectory(0); // Release from file container
  return ret;
}

//____________________________________________________________________
void AliTrackletAODUtils::CopyAttr(const TH1* src, TH1* dest)
{
  if (!src || !dest) return;
  dest->SetMarkerStyle(src->GetMarkerStyle());
  dest->SetMarkerColor(src->GetMarkerColor());
  dest->SetMarkerSize (src->GetMarkerSize());
  dest->SetLineStyle  (src->GetLineStyle());
  dest->SetLineColor  (src->GetLineColor());
  dest->SetLineWidth  (src->GetLineWidth());
  dest->SetFillStyle  (src->GetFillStyle());
  dest->SetFillColor  (src->GetFillColor());
}


//____________________________________________________________________
TH1* AliTrackletAODUtils::Make1D(Container*     c,
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
  if (c) c->Add(ret);
  return ret;
}
//____________________________________________________________________
TProfile* AliTrackletAODUtils::Make1P(Container*     c,
				      const TString& name,
				      const TString& title,
				      Color_t        color,
				      Style_t        style,
				      const TAxis&   xAxis)
{
  TString   n   = name;
  TString   t   = title;
  TProfile* ret = 0;
  Int_t     nx  = xAxis.GetNbins();
  if (t.IsNull()) t = xAxis.GetTitle();
  if (xAxis.GetXbins() && xAxis.GetXbins()->GetArray())
    ret = new TProfile(n,t,nx,xAxis.GetXbins()->GetArray());
  else
    ret = new TProfile(n,t,nx,xAxis.GetXmin(),xAxis.GetXmax());
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
  if (c) c->Add(ret);
  return ret;
}
//____________________________________________________________________
TH2* AliTrackletAODUtils::Make2D(Container*     c,
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
  if (c) c->Add(ret);
  return ret;
}
//____________________________________________________________________
TH3* AliTrackletAODUtils::Make3D(Container*     c,
				 const TString& name,
				 const TString& title,
				 Color_t        color,
				 Style_t        style,
				 const TAxis&   xAxis,
				 const TAxis&   yAxis,
				 const TAxis&   zAxis)
{
  TString   n   = name;
  TString   t   = title;
  TH3*      ret = 0;
  Int_t     nx  = xAxis.GetNbins();
  Int_t     ny  = yAxis.GetNbins();
  Int_t     nz  = zAxis.GetNbins();
  Double_t* xb  = (xAxis.GetXbins() && xAxis.GetXbins()->GetArray() ?
		   const_cast<Double_t*>(xAxis.GetXbins()->GetArray()) : 0);
  Double_t* yb  = (yAxis.GetXbins() && yAxis.GetXbins()->GetArray() ?
		   const_cast<Double_t*>(yAxis.GetXbins()->GetArray()) : 0);
  Double_t* zb  = (zAxis.GetXbins() && zAxis.GetXbins()->GetArray() ?
		   const_cast<Double_t*>(zAxis.GetXbins()->GetArray()) : 0);
  if (t.IsNull())
    t.Form("%s vs %s vs %s",
	   zAxis.GetTitle(), yAxis.GetTitle(), xAxis.GetTitle());
  if (xb || yb || zb) {
    // One or more axis are defined as arrays.  Make sure the rest are
    // also arrays.
    if (xb) {
      xb = new Double_t[nx+1];
      Double_t dx = (xAxis.GetXmax()-xAxis.GetXmin())/nx;
      xb[0] = xAxis.GetXmin();
      for (Int_t i = 1; i <= nx; i++) xb[i] = xb[i-1]+dx;
    }
    if (yb) {
      yb = new Double_t[ny+1];
      Double_t dy = (yAxis.GetXmax()-yAxis.GetXmin())/ny;
      yb[0] = yAxis.GetXmin();
      for (Int_t i = 1; i <= ny; i++) yb[i] = yb[i-1]+dy;
    }
    if (zb) {
      zb = new Double_t[nz+1];
      Double_t dz = (zAxis.GetXmax()-zAxis.GetXmin())/nz;
      zb[0] = zAxis.GetXmin();
      for (Int_t i = 1; i <= nz; i++) zb[i] = zb[i-1]+dz;
    }
    ret = new TH3D(n,t,nx,xb,ny,yb,nz,zb);
  }
  else  {
    ret = new TH3D(n,t,
		   nx, xAxis.GetXmin(), xAxis.GetXmax(),
		   ny, yAxis.GetXmin(), yAxis.GetXmax(),
		   nz, zAxis.GetXmin(), zAxis.GetXmax());
  }
  ret->Sumw2();
  ret->SetXTitle(xAxis.GetTitle());
  ret->SetYTitle(yAxis.GetTitle());
  ret->SetZTitle(zAxis.GetTitle());
  ret->SetLineColor(color);
  ret->SetMarkerColor(color);
  ret->SetFillColor(color);
  ret->SetMarkerStyle(style);
  static_cast<const TAttAxis&>(xAxis).Copy(*(ret->GetXaxis()));
  static_cast<const TAttAxis&>(yAxis).Copy(*(ret->GetYaxis()));
  static_cast<const TAttAxis&>(zAxis).Copy(*(ret->GetZaxis()));
  ret->SetDirectory(0);
  if (const_cast<TAxis&>(xAxis).GetLabels()) {
    for (Int_t i = 1; i <= xAxis.GetNbins(); i++)
      ret->GetXaxis()->SetBinLabel(i, xAxis.GetBinLabel(i));
  }
  if (const_cast<TAxis&>(yAxis).GetLabels()) {
    for (Int_t i = 1; i <= yAxis.GetNbins(); i++)
      ret->GetYaxis()->SetBinLabel(i, yAxis.GetBinLabel(i));
  }
  if (const_cast<TAxis&>(zAxis).GetLabels()) {
    for (Int_t i = 1; i <= zAxis.GetNbins(); i++)
      ret->GetZaxis()->SetBinLabel(i, zAxis.GetBinLabel(i));
  }
  if (c) c->Add(ret);
  return ret;
}  
//____________________________________________________________________
TProfile2D* AliTrackletAODUtils::Make2P(Container*     c,
					const TString& name,
					const TString& title,
					Color_t        color,
					Style_t        style,
					const TAxis&   xAxis,
					const TAxis&   yAxis)
{
  TString         n   = name;
  TString         t   = title;
  TProfile2D*     ret = 0;
  Int_t           nx  = xAxis.GetNbins();
  Int_t           ny  = yAxis.GetNbins();
  const Double_t* xb  = (xAxis.GetXbins() && xAxis.GetXbins()->GetArray() ?
			 xAxis.GetXbins()->GetArray() : 0);
  const Double_t* yb  = (yAxis.GetXbins() && yAxis.GetXbins()->GetArray() ?
			 yAxis.GetXbins()->GetArray() : 0);
  if (t.IsNull())
    t.Form("%s vs %s", yAxis.GetTitle(), xAxis.GetTitle());
  if (xb) {	  
    if   (yb) ret = new TProfile2D(n,t,nx,xb,ny,yb);
    else      ret = new TProfile2D(n,t,
				   nx,xb,
				   ny,yAxis.GetXmin(),yAxis.GetXmax());
  }
  else { 
    if   (yb) ret = new TProfile2D(n,t,
				   nx,xAxis.GetXmin(),xAxis.GetXmax(),
				   ny,yb);
    else      ret = new TProfile2D(n,t,
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
  if (c) c->Add(ret);
  return ret;
}
//____________________________________________________________________
void AliTrackletAODUtils::FixAxis(TAxis& axis, const char* title)
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
void AliTrackletAODUtils::ScaleAxis(TAxis&       ret,
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
void AliTrackletAODUtils::SetAxis(TAxis& axis, Int_t n, Double_t* borders)
{
  axis.Set(n, borders);
  FixAxis(axis);
}
//____________________________________________________________________
void AliTrackletAODUtils::SetAxis(TAxis&         axis,
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
void AliTrackletAODUtils::SetAxis(TAxis&   axis,
				    Int_t    n,
				    Double_t l,
				    Double_t h)
{
  axis.Set(n, l, h);
  FixAxis(axis);
}
//____________________________________________________________________
void AliTrackletAODUtils::SetAxis(TAxis& axis, Int_t n, Double_t m)
{
  SetAxis(axis, n, -TMath::Abs(m), +TMath::Abs(m));
}
//____________________________________________________________________
void AliTrackletAODUtils::PrintAxis(const TAxis& axis,
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
TH2* AliTrackletAODUtils::ScaleToIPz(TH2* h, TH1* ipZ, Bool_t full)
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
TH3* AliTrackletAODUtils::ScaleToIPz(TH3* h, TH1* ipZ, Bool_t full)
{
  if (!h) {
    ::Warning("ScaleToIPz","Nothing to scale");
    return 0;
  }
  if (!ipZ) {
    ::Warning("ScaleToIPz","Nothing to scale by");
    return 0;
  }
  TH3* ret = static_cast<TH3*>(h->Clone());
  ret->SetDirectory(0);
  if (!ipZ) return ret;
  for (Int_t iz = 1; iz <= ret->GetNbinsZ(); iz++) {
    Double_t z   = ret->GetZaxis()->GetBinCenter(iz);
    Int_t    bin = ipZ->GetXaxis()->FindBin(z);
    Double_t nEv = ipZ->GetBinContent(bin);
    Double_t eEv = ipZ->GetBinError  (bin);
    Double_t esc = (nEv > 0 ? 1./nEv : 0);
    Double_t rE2 = esc*esc*eEv*eEv;
    for (Int_t iy = 1; iy <= ret->GetNbinsY(); iy++) {
      for (Int_t ix = 1; ix <= ret->GetNbinsX(); ix++) {
	Double_t  c   = ret->GetBinContent(ix,iy,iz);
	Double_t  e   = ret->GetBinError  (ix,iy,iz);
	Double_t  r   = (c > 0 ? e/c : 0);
	// Scale by number of events, and error propagate 
	Double_t  sc  = c * esc;
	Double_t  se  = 0;
	if (full) se  = sc * TMath::Sqrt(r*r+rE2);
	else      se  = e * esc;
	Double_t scl = 1 / ret->GetXaxis()->GetBinWidth(ix);
	ret->SetBinContent(ix, iy, iz, scl*sc);
	ret->SetBinError  (ix, iy, iz, scl*se);
      }
    }
  }
  return ret;
}
//____________________________________________________________________
TH3* AliTrackletAODUtils::ScaleDelta(TH3* h, TH2* etaIPzScale)
{
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) { // eta
    for (Int_t j = 1; j <= h->GetNbinsZ(); j++) { // IPz
      Double_t scale  = etaIPzScale->GetBinContent(i,j);
      Double_t scaleE = etaIPzScale->GetBinError  (i,j);
      Double_t q      = (scale > 0 ? scaleE / scale : 0);
      for (Int_t k = 1; k <= h->GetNbinsY(); k++) { // delta
	Double_t c = h->GetBinContent(i,k,j);
	Double_t e = h->GetBinError  (i,k,j);
	Double_t r = (c > 0 ? e / c : 0);
	Double_t w = c * scale;
	Double_t v = w * TMath::Sqrt(r*r + q*q);
	h->SetBinContent(i,k,j,w);
	h->SetBinError  (i,k,j,v);
#if 0
	Printf("%2d,%3d,%2d=(%9g+/-%9g)*(%9g+/-%9g)=(%9g+/-%9g)",
	       i,k,j,c,e,scale,scaleE,w,v);
#endif 
      }
    }
  }
  return h;
}
//____________________________________________________________________
TH2* AliTrackletAODUtils::ProjectEtaDelta(TH3* h)
{
  TH2* etaDelta = static_cast<TH2*>(h->Project3D("yx e"));
  etaDelta->SetName("etaDelta");
  etaDelta->SetTitle(h->GetTitle());
  etaDelta->SetDirectory(0);
  etaDelta->SetZTitle("d^{2}#it{N}/(d#etad#Delta)");
#if 1
  // Reset content of projected histogram and calculate averages
  etaDelta->Reset();
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) { // Loop over eta
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) { // Loop over Delta
      Double_t sum  = 0;
      Double_t sumw = 0;
      Int_t    cnt  = 0;
      for (Int_t k = 1; k <= h->GetNbinsZ(); k++) { // Loop over IPz
	Double_t c  = h->GetBinContent(i,j,k);
	Double_t e  = h->GetBinError  (i,j,k);
	if (c < 1e-6 || e/c > 1) continue;
#if 0
	Double_t w  =  1/e/e;
	sum         += w * c;
	sumw        += w;
#else
	sum         += c;
	sumw        += e*e;
#endif 
	cnt         += 1;
      }
      if (sumw < 1e-6) continue;
#if 0
      etaDelta->SetBinContent(i,j,sum/sumw);
      etaDelta->SetBinError  (i,j,TMath::Sqrt(1/sumw));
#else
      etaDelta->SetBinContent(i,j,sum/cnt);
      etaDelta->SetBinError  (i,j,TMath::Sqrt(sumw)/cnt);
#endif 
    }
  }
#endif 
  return etaDelta;
}
//____________________________________________________________________
TH1* AliTrackletAODUtils::ProjectDelta(TH2* h)
{
  TH1* delta = h->ProjectionY("delta");
  delta->SetDirectory(0);
  delta->SetTitle(h->GetTitle());
  delta->SetYTitle("d#it{N}/d#Delta");
  delta->Scale(1./h->GetNbinsX());
  return delta;
}
//____________________________________________________________________
TH1* AliTrackletAODUtils::ProjectDeltaFull(TH3* h)
{
  TH2* tmp   = ProjectEtaDelta(h);  
  TH1* delta = ProjectDelta(tmp);
  delta->SetDirectory(0);
  tmp->SetDirectory(0);
  delete tmp;
  return delta;
}

//____________________________________________________________________
TH1* AliTrackletAODUtils::AverageOverIPz(TH2*        h,
					 const char* name,
					 UShort_t    mode,
					 TH1*        ipz,
					 TH2*        other,
					 Bool_t      verb)
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
      if (verb) 
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
TObject* AliTrackletAODUtils::CloneAndAdd(Container* c, TObject* o)
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
Double_t AliTrackletAODUtils::Integrate(TH1*     h,
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
Double_t AliTrackletAODUtils::RatioE(Double_t n, Double_t en,
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
Double_t AliTrackletAODUtils::RatioE2(Double_t n, Double_t e2n,
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

//____________________________________________________________________
Int_t* AliTrackletAODUtils::PdgArray(Int_t& size)
{
  static Int_t codes[] = {
    211,        // #pi^{+}			 
    2212, 	  // p				 
    321, 	  // K^{+}			 
    323, 	  // K^{*+}			 
    11, 	  // e^{-}			 
    13, 	  // #mu^{-}			 
    213, 	  // #rho^{+}			 
    411, 	  // D^{+}			 
    413, 	  // D^{*+}			 
    431, 	  // D_{s}^{+}			 
    433, 	  // D_{s}^{*+}			 
    1114, 	  // #Delta^{-}			 
    2214, 	  // #Delta^{+}			 
    2224, 	  // #Delta^{++}			 
    3112, 	  // #Sigma^{-}			 
    3222, 	  // #Sigma^{+}			 
    3114, 	  // #Sigma^{*-}			 
    3224, 	  // #Sigma^{*+}			 
    4214, 	  // #Sigma^{*+}_{c}		 
    4224, 	  // #Sigma^{*++}_{c}		 
    3312, 	  // #Xi^{-}			 
    3314, 	  // #Xi^{*-}			 
    4122, 	  // #Lambda^{+}_{c}		 
    2112, 	  // n				 
    2114, 	  // #Delta^{0}			 
    22, 	  // #gamma			 
    310, 	  // K^{0}_{S}			 
    130, 	  // K^{0}_{L}			 
    311, 	  // K^{0}			 
    313, 	  // K^{*}			 
    221, 	  // #eta				 
    111, 	  // #pi^{0}			 
    113, 	  // #rho^{0}			 
    333, 	  // #varphi			 
    331, 	  // #eta'			 
    223, 	  // #omega			 
    3122, 	  // #Lambda			 
    3212, 	  // #Sigma^{0}			 
    4114, 	  // #Sigma^{*0}_{c}		 
    3214, 	  // #Sigma^{*0}			 
    421, 	  // D^{0}			 
    423, 	  // D^{*0}			 
    3322, 	  // #Xi_{0}			 
    3324, 	  // #Xi^{*0}			 
    4132, 	  // #Xi^{0}_{c}			 
    4314,	  // #Xi^{*0}_{c}
    1000000000  // Nuclei			 
  };
  static Int_t  nCodes = sizeof(codes) / sizeof(Int_t);
  static Int_t* sorted = 0;
  size = nCodes;
  if (sorted) return sorted;
  
  sorted     = new Int_t[nCodes];
  Int_t* idx = new Int_t[nCodes];
  TMath::Sort(nCodes, codes, idx, false);
  for (Int_t i = 0; i < nCodes; i++) sorted[i] = codes[idx[i]];
  delete [] idx;
  return sorted;
}

//____________________________________________________________________
void AliTrackletAODUtils::PdgAttr(Int_t pdg, TString& nme,
				  Color_t& c, Style_t& s)
{
  // Leptons are black
  // Non-strange, non-charm meson are red
  // Non-strange, non-charm baryons are magenta
  // Strange-mesons are blue 
  // Strange-baryons are cyan
  // Charmed mesons are green
  // Charmed baryons are yellow
  switch (pdg) {
  case -1:        c=kPink   +7; s=20; nme="?";               break;
  case 11: 	  c=kBlack  +0; s=20; nme="e^{-}"; 	     break;
  case 13: 	  c=kBlack  +0; s=21; nme="#mu^{-}"; 	     break;
  case 22: 	  c=kBlack  +0; s=22; nme="#gamma"; 	     break;
  case 111: 	  c=kRed    +2; s=20; nme="#pi^{0}"; 	     break;
  case 211:       c=kRed    +2; s=21; nme="#pi^{+}";         break;
  case 213: 	  c=kRed    +2; s=22; nme="#rho^{+}"; 	     break;
  case 113: 	  c=kRed    +2; s=23; nme="#rho^{0}"; 	     break;
  case 223: 	  c=kRed    +2; s=24; nme="#omega"; 	     break;
  case 321: 	  c=kBlue   +2; s=20; nme="K^{+}"; 	     break;
  case 323: 	  c=kBlue   +2; s=21; nme="K^{*+}"; 	     break;
  case 310: 	  c=kBlue   +2; s=22; nme="K^{0}_{S}"; 	     break;
  case 130: 	  c=kBlue   +2; s=23; nme="K^{0}_{L}"; 	     break;
  case 311: 	  c=kBlue   +2; s=24; nme="K^{0}"; 	     break;
  case 313: 	  c=kBlue   +2; s=25; nme="K^{*}"; 	     break;
  case 221: 	  c=kBlue   +2; s=26; nme="#eta"; 	     break;
  case 333: 	  c=kBlue   +2; s=27; nme="#varphi"; 	     break;
  case 331: 	  c=kBlue   +2; s=28; nme="#eta'"; 	     break;
  case 411: 	  c=kGreen  +2; s=20; nme="D^{+}"; 	     break;
  case 413: 	  c=kGreen  +2; s=21; nme="D^{*+}"; 	     break;
  case 421: 	  c=kGreen  +2; s=22; nme="D^{0}"; 	     break;
  case 423: 	  c=kGreen  +2; s=23; nme="D^{*0}"; 	     break;
  case 431: 	  c=kGreen  +2; s=24; nme="D_{s}^{+}"; 	     break;
  case 433: 	  c=kGreen  +2; s=25; nme="D_{s}^{*+}";      break;
  case 2212: 	  c=kMagenta+2; s=20; nme="p"; 	             break;
  case 2112: 	  c=kMagenta+2; s=21; nme="n"; 	             break;
  case 2114: 	  c=kMagenta+2; s=22; nme="#Delta^{0}";      break;
  case 1114: 	  c=kMagenta+2; s=23; nme="#Delta^{-}";      break;
  case 2214: 	  c=kMagenta+2; s=24; nme="#Delta^{+}";	     break;
  case 2224: 	  c=kMagenta+2; s=25; nme="#Delta^{++}";     break;
  case 3112: 	  c=kCyan   +2; s=20; nme="#Sigma^{-}";      break;
  case 3222: 	  c=kCyan   +2; s=21; nme="#Sigma^{+}";      break;
  case 3114: 	  c=kCyan   +2; s=22; nme="#Sigma^{*-}";     break;
  case 3224: 	  c=kCyan   +2; s=23; nme="#Sigma^{*+}";     break;
  case 3312: 	  c=kCyan   +2; s=24; nme="#Xi^{-}"; 	     break;
  case 3314: 	  c=kCyan   +2; s=25; nme="#Xi^{*-}"; 	     break;
  case 3122: 	  c=kCyan   +2; s=26; nme="#Lambda"; 	     break;
  case 3212: 	  c=kCyan   +2; s=27; nme="#Sigma^{0}";      break;
  case 3214: 	  c=kCyan   +2; s=28; nme="#Sigma^{*0}";     break;
  case 3322: 	  c=kCyan   +2; s=29; nme="#Xi_{0}"; 	     break;
  case 3324: 	  c=kCyan   +2; s=30; nme="#Xi^{*0}"; 	     break;
  case 4214: 	  c=kYellow +2; s=20; nme="#Sigma^{*+}_{c}"; break;
  case 4224: 	  c=kYellow +2; s=21; nme="#Sigma^{*++}_{c}";break;
  case 4122: 	  c=kYellow +2; s=22; nme="#Lambda^{+}_{c}"; break;
  case 4114: 	  c=kYellow +2; s=23; nme="#Sigma^{*0}_{c}"; break;
  case 4132: 	  c=kYellow +2; s=24; nme="#Xi^{0}_{c}";     break;
  case 4314:	  c=kYellow +2; s=25; nme="#Xi^{*0}_{c}";    break;
  case 1000000000:c=kPink   +2; s=20; nme="Nuclei";          break;
  default:        c=kGray;      s=1;  nme.Form("%d",pdg);    break;
  };
}

//____________________________________________________________________
Int_t AliTrackletAODUtils::PdgBin(Int_t pdg)
{
  Int_t  size;
  Int_t* array = PdgArray(size);
  Int_t  apdg  = TMath::Abs(pdg);
  Int_t  idx   = TMath::BinarySearch(size, array, apdg);
  if (idx        == size - 1) return idx+1;
  if (array[idx] != apdg)     return size+1;
  return idx+1;
}

//____________________________________________________________________
const TAxis& AliTrackletAODUtils::PdgAxis()
{
  static TAxis* axis = 0;
  if (axis) return *axis;

  Int_t         size;
  Int_t*        array = PdgArray(size);
  axis                = new TAxis(size+1, +.5, size+1.5);
  // TDatabasePDG* pdgDb = TDatabasePDG::Instance();
  for (Int_t i = 1; i <= size; i++) {
    Int_t         pdg  = array[i-1];
    axis->SetBinLabel(i, Form("%d", pdg));
    // TString       name = "?";
    // TParticlePDG* type = pdgDb->GetParticle(pdg);
    // if (type)     name = type->GetName();
    // axis->SetBinLabel(i, name.Data());
  }
  axis->SetBinLabel(size+1, "-1");
  return *axis;
}
    
#endif
// Local Variables:
//  mode: C++
// End:

