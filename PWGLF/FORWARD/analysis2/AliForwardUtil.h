// 
// Utilities used in the forward multiplcity analysis 
// 
//
#ifndef ALIFORWARDUTIL_H
#define ALIFORWARDUTIL_H
/**
 * @file   AliForwardUtil.h
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Mar 23 14:06:54 2011
 * 
 * @brief  Various utilities used in PWGLF/FORWARD
 * 
 * @ingroup pwglf_forward 
 */
#include <TObject.h>
#include <TString.h>
#include <TObjArray.h>
class TH2D;
class TH1I;
class TH1;
class TF1;
class TAxis;
class TArrayD;
class TVector3;
class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliAnalysisTaskSE;

/** 
 * Utilities used in the forward multiplcity analysis 
 * 
 * @ingroup pwglf_forward 
 */
class AliForwardUtil : public TObject
{
public:
  enum {
    /** Value if things cannot be calculated */
    kInvalidValue = -9999
  };
  enum { 
    kSkipRing = (1 << 19) // Bit for skipping a histogram
  };
  /** 
   * Get the standard color for a ring  
   *
   * @param d Detector
   * @param r Ring 
   * 
   * @return 
   */
  static Color_t RingColor(UShort_t d, Char_t r)
  {
#if 0
    return ((d == 1 ? kRed : (d == 2 ? kGreen : kBlue))
	    + ((r == 'I' || r == 'i') ? 2 : -3));
#else
    // Brigher colours that look better in print 
    return ((d == 1 ? kRed : (d == 2 ? kGreen : kBlue))
	    + ((r == 'I' || r == 'i') ? -3 : -9));
#endif
  }
  //==================================================================
  /** 
   * @{ 
   * @name AliROOT version
   */
  /** 
   * Get the revision number of AliROOT
   * 
   * @return Subversion revision number of AliROOT used
   */
  static ULong_t AliROOTRevision();
  /**
   * Get the branch identifier of AliROOT.  In case of trunk, return
   * 0xFFFFFFFF, while for @b vM-N-R{-S}, we get
   *
   * @code 
   *    ((M & 0xFF) << 12 | (N & 0xFF) << 8 | (R & 0xFF) << 3 | (X))
   * @endcode 
   * where @c X is 0xAA if @b S is specified (e.g., analysis tag). 
   *
   * @return branch identifer encoded in bits 
   */
  static ULong_t AliROOTBranch();
  /* @} */
  //==================================================================
  /** 
   * @{ 
   * @name Collision/run parameters 
   */
  /**						
   * Defined collision types 
   */
  enum ECollisionSystem {
    kUnknown, 
    kPP, 
    kPbPb,
    kPPb,
    kPbp,
    kXeXe
  };
  //__________________________________________________________________
  /** 
   * Calculate the beam rapidity. 
   *
   * @b Note: The beam energy is given in GeV/charge
   * 
   * @param beam Beam energy in GeV/charge 
   * @param z    Charge number of projectile
   * @param a    Mass number of projectile 
   * 
   * @return The rapidity of the beam
   */
  static Float_t BeamRapidity(Float_t beam, UShort_t z, UShort_t a);
  /** 
   * Calculate the center of mass energy from the beam energy per
   * charge and the nucleus numbers.
   * 
   * @param beam Beam energy in GeV/charge 
   * @param z1   Charge number of projectile
   * @param a1   Mass number of projectile
   * @param z2   Charge number of projectile
   * @param a2 	 Mass number of projectile  
   * 
   * @return The center of mass energy in GeV/nucleon
   */
  static Float_t CenterOfMassEnergy(Float_t beam, UShort_t z1, UShort_t a1, 
				    Short_t z2=-1, Short_t a2=-1);
  /** 
   * Calculate the center of mass rapidity (shift)
   * 
   * @param z1 Charge number of projectile
   * @param a1 Mass number of projectile  
   * @param z2 Charge number of projectile
   * @param a2 Mass number of projectile  
   * x
   * @return Rapidity of the center of mass 
   */
  static Float_t CenterOfMassRapidity(UShort_t z1, UShort_t a1, 
				      Short_t z2=-1, Short_t a2=-1);
  /** 
   * Parse a collision system spec given in a string.   Known values are 
   * 
   *  - "pp", "p-p" which returns kPP 
   *  - "PbPb", "Pb-Pb", "A-A", which returns kPbPb
   *  - "pPb", "p-Pb", "pA", p-A" which returns kPPb
   *  - Everything else gives kUnknown 
   * 
   * @param sys Collision system spec 
   * 
   * @return Collision system id 
   */
  static UShort_t ParseCollisionSystem(const char* sys);
  /** 
   * Parse a collision system spec given in a string.   Known values are 
   * 
   *  - "pp", "p-p" which returns kPP 
   *  - "PbPb", "Pb-Pb", "A-A", which returns kPbPb
   *  - "pPb", "p-Pb", "pA", p-A" which returns kPPb
   *  - Everything else gives kUnknown 
   * 
   * @param sys Collision system spec 
   * 
   * @return Collision system id 
   */
  static UShort_t ParseCollisionSystem(Int_t b1a, Int_t b1z,
				       Int_t b2a, Int_t b2z,
				       const char* sys);
  /** 
   * Get a string representation of the collision system 
   * 
   * @param sys  Collision system 
   * - kPP -> "pp"
   * - kPbPb -> "PbPb" 
   * - kPPb -> "pPb"
   * - anything else gives "unknown"
   * 
   * @return String representation of the collision system 
   */
  static const char* CollisionSystemString(UShort_t sys);
  //__________________________________________________________________
  /** 
   * Parse the center of mass energy given as a float and return known 
   * values as a unsigned integer
   * 
   * @param sys  Collision system (needed for AA)
   * @param cms  Center of mass energy * total charge 
   * 
   * @return Center of mass energy per nucleon
   */
  static UShort_t ParseCenterOfMassEnergy(UShort_t sys, Float_t cms);
  /** 
   * Get a string representation of the center of mass energy per nuclean
   * 
   * @param cms  Center of mass energy per nucleon
   * 
   * @return String representation of the center of mass energy per nuclean
   */
  static const char* CenterOfMassEnergyString(UShort_t cms);
  //__________________________________________________________________
  /** 
   * Parse the magnetic field (in kG) as given by a floating point number
   * 
   * @param field  Magnetic field in kG 
   * 
   * @return Short integer value of magnetic field in kG 
   */
  static Short_t ParseMagneticField(Float_t field);
  /** 
   * Get a string representation of the magnetic field
   * 
   * @param field Magnetic field in kG
   * 
   * @return String representation of the magnetic field
   */
  static const char* MagneticFieldString(Short_t field);
  /* @} */

  //==================================================================
  /** 
   * @{ 
   * @name Recalculate @f$\eta@f$, @f$\phi@f$, etc. 
   */
  /** 
   * Get the radius of a strip. 
   * 
   * @param ring  Ring identifier 'I' or 'O'
   * @param strip Strip number 
   * 
   * @return Radial distance from beam of the strip 
   */
  static Double_t GetStripR(Char_t ring, UShort_t strip);
  /** 
   * Get the location of a sector in Z (in centimeters) 
   * 
   * @param det   Detector 
   * @param ring  Ring 
   * @param sec   Sector 
   * 
   * @return Z position 
   */
  static Double_t GetSectorZ(UShort_t det, Char_t ring, UShort_t sec);
  /** 
   * Get the azimythal angle of sector (in radians)
   * 
   * @param det   Detector
   * @param ring  Ring
   * @param sec   Sector 
   * 
   * @return Azimuthal angle (in radians)
   */
  static Double_t GetSectorPhi(UShort_t det, Char_t ring, UShort_t sec);
  /** 
   * Get the coordinates of a strip relative to the interaction point
   * 
   * @param det   Detector 
   * @param ring  Ring
   * @param sec   Sector 
   * @param str   Strip 
   * @param ip   The interaction point
   * @param pos  On return, the relative position 
   *
   * @return true on success 
   */
  static Bool_t GetXYZ(UShort_t det, Char_t ring, UShort_t sec,
		       UShort_t str, const TVector3& ip,
		       TVector3& pos);
  /** 
   * Get the eta and phi of a strip given an interaction point 
   * 
   * @param det   Detector 
   * @param ring  Ring
   * @param sec   Sector 
   * @param str   Strip 
   * @param ip    Interaction point 
   * @param eta   On return, the eta
   * @param phi   On return, the phi (in radians)
   *
   * @return true on success 
   */
  static Bool_t GetEtaPhi(UShort_t det, Char_t ring, UShort_t sec,
			  UShort_t str, const TVector3& ip,
			  Double_t& eta, Double_t& phi);
  /** 
   * Get eta from strip
   * 
   * @param det, ring, sec, strip, zvtx
   * 
   * @return eta
   */
  static Double_t GetEtaFromStrip(UShort_t det, Char_t ring, 
				  UShort_t sec, UShort_t strip, Double_t zvtx);
  /** 
   * Recalculate the eta of a strip.  This will correct for an off-set IP
   * 
   * @param r      Ring identifier 
   * @param strip  Strip number 
   * @param eta    Nominal eta, on return corrected
   * @param phi    Nominal phi, on return corrected
   * @param ipX    IP X coordinate 
   * @param ipY    IP Y coordinate 
   * 
   * @return true on success 
   */
  static Bool_t GetEtaPhiFromStrip(Char_t r, UShort_t strip,
				   Double_t& eta, Double_t& phi , 
				   Double_t ipX,  Double_t ipY);
  /** 
   * Get the azimuthal angle of a strip
   * 
   * @param ring  Ring identifier 'I' or 'O'
   * @param strip Strip number 
   * @param phi   Straight forward strip phi
   * @param xvtx  Ip X coordinate 
   * @param yvtx  Ip Y coordinate 
   * 
   * @return The phi angle correctef for (X,Y) off set. 
   */  
  static Double_t GetPhiFromStrip(Char_t ring, UShort_t strip, 
				  Double_t phi, Double_t xvtx, Double_t yvtx);
  /* @} */

  //==================================================================
  /** 
   * @{ 
   * @name Manager related tasks 
   */
  /** 
   * Get the AOD event - either from the input (AOD analysis) or the
   * output (ESD analysis)
   * 
   * @param task Task to do the investigation for
   * 
   * @return Found AOD event or null
   */
  static AliAODEvent* GetAODEvent(AliAnalysisTaskSE* task);
  /** 
   * Check if we have something that will provide and AOD event 
   * 
   * @return 0 if there's nothing that provide an AOD event, 1 if it
   * is provided on the input (AOD analysis) or 2 if it is provided on
   * the output (ESD analysis)
   */
  static UShort_t CheckForAOD();
  /** 
   * Check if we have a particular (kind) of task in our train
   * 
   * @param clsOrName  Class name or name of task 
   * @param cls If true, look for a task of a particular class -
   * otherwise search for a speficially name task
   * 
   * @return true if the needed task was found 
   */
  static Bool_t CheckForTask(const char* clsOrName, Bool_t cls=true);
  /* @} */

  //==================================================================
  /** 
   * @{ 
   * @name Member functions to store and retrieve analysis parameters 
   */
  /** 
   * Make a parameter 
   * 
   * @param name   Name of the parameter 
   * @param value  Value of the parameter 
   * 
   * @return Newly allocated parameter object 
   */
  static TObject* MakeParameter(const char* name, UShort_t value);
  /** 
   * Make a parameter 
   * 
   * @param name   Name of the parameter 
   * @param value  Value of the parameter 
   * 
   * @return Newly allocated parameter object 
   */
  static TObject* MakeParameter(const char* name, Int_t value);
  /** 
   * Make a parameter 
   * 
   * @param name   Name of the parameter 
   * @param value  Value of the parameter 
   * 
   * @return Newly allocated parameter object 
   */
  static TObject* MakeParameter(const char* name, Double_t value);
  /** 
   * Make a parameter 
   * 
   * @param name   Name of the parameter 
   * @param value  Value of the parameter 
   * 
   * @return Newly allocated parameter object 
   */
  static TObject* MakeParameter(const char* name, Bool_t value);
  /** 
   * Make a parameter 
   * 
   * @param name   Name of the parameter 
   * @param value  Value of the parameter 
   * 
   * @return Newly allocated parameter object 
   */
  static TObject* MakeParameter(const char* name, ULong_t value);
  /** 
   * Get a parameter value from an object
   * 
   * @param o      Object
   * @param value  On return, the parameter value 
   */
  static void GetParameter(TObject* o, UShort_t& value);
  /** 
   * Get a parameter value from an object
   * 
   * @param o      Object
   * @param value  On return, the parameter value 
   */
  static void GetParameter(TObject* o, Int_t& value);
  /** 
   * Get a parameter value from an object
   * 
   * @param o      Object
   * @param value  On return, the parameter value 
   */
  static void GetParameter(TObject* o, Double_t& value);
  /** 
   * Get a parameter value from an object
   * 
   * @param o      Object
   * @param value  On return, the parameter value 
   */
  static void GetParameter(TObject* o, Bool_t& value);
  /** 
   * Get a parameter value from an object
   * 
   * @param o      Object
   * @param value  On return, the parameter value 
   */
  static void GetParameter(TObject* o, ULong_t& value);
  /* @} */

  //==================================================================
  /** 
   * @{ 
   * @name Axis functions 
   */
  /** 
   * Make a full @f$ \mbox{IP}_z@f$ axis 
   * 
   * @param nCenter Number of bins in the center 
   * 
   * @return Newly allocated axis object. 
   */
  static TAxis* MakeFullIpZAxis(Int_t nCenter=20);
  /** 
   * Make a full @f$ \mbox{IP}_z@f$ axis 
   * 
   * @param nCenter Number of bins in the center 
   * @param bins    On return, the bin edges 
   */
  static void MakeFullIpZAxis(Int_t nCenter, TArrayD& bins);
  /** 
   * Make a log-scale axis.  That is, each bin has constant with when
   * drawn on a log scale.
   * 
   * @param nBins      Number of bins
   * @param minOrder   Least power of 10 
   * @param maxOrder   Largest power of 10
   * @param bins       On return, the bin edges
   */
  static void MakeLogScale(Int_t nBins, Int_t minOrder, Int_t maxOrder,
			   TArrayD& bins);
  /* @} */

  //==================================================================
  /** 
   * @{ 
   * @name Centrality functions 
   */
  /** 
   * Get the centrality of the event.  One can pass a both an ESD and
   * AOD event.  Also, there's fall-back to the old centrality objects
   * in case the newer AliMultSelection object isn't found.
   * 
   * Possible values of qual 
   *
   * - 0   All is good and event within calibrated sample 
   * - 198 Centrality beyond anchor point 
   * - 199 No calibration for estimator (return <0)
   * - 200 Not the same trigger as in calibration 
   * - 201 Not an INEL>0 (if used in calibration)
   * - 202 Event outside IPz cut used in calibration 
   * - 203 Event is a pile-up kind not used in calibration 
   * - 204 Not consistent SPD/TPC IP (if used in calibration)
   * - 205 SPD outlier event (if used in calibration)
   * - 206 IP has less than 1 contributor (if used in calibration)
   * - 0xFFFF No estimate of centrality could be done (return <0)
   *
   * @param event    Event 
   * @param method   Centrality method to use 
   * @param qual     On return, the quality - 0 is good 
   * @param verbose  If true, be verbose 
   * 
   * @return The centrality or -1 in case of problems. 
   */
  static Float_t GetCentrality(const AliVEvent& event, 
			       const TString&   method, 
			       Int_t&           qual, 
			       Bool_t           verbose=false);
  //____________________________________________________________________
  /** 
   * Get the centrality of the event.  One can pass a both an ESD and
   * AOD event.  
   * 
   * @param event    Event 
   * @param method   Centrality method to use 
   * @param qual     On return, the quality - 0 is good 
   * @param verbose  If true, be verbose 
   * 
   * @return The centrality or -1 in case of problems. 
   */
  static Float_t GetCentralityMult(const AliVEvent& event, 
				   const TString&   method, 
				   Int_t&           qual, 
				   Bool_t           verbose=false);
  /** 
   * Compatbility function to retrieve the centrality from the event.
   * 
   * @param event   Event 
   * @param method  Centrality method 
   * @param qual    Quality flat - 0 is good 
   * @param verbose If true, be verbose 
   * 
   * @return The centrality or -1 in case of problems 
   */
  static Float_t GetCentralityCompat(const AliVEvent& event,
				     const TString&   method,
				     Int_t&           qual,
				     Bool_t           verbose);
  //____________________________________________________________________
  /** 
   * Compatbility function to retrieve the centrality from the ESD event. 
   * 
   * @param event   Event 
   * @param method  Centrality method 
   * @param qual    Quality flat - 0 is good 
   * @param verbose If true, be verbose 
   * 
   * @return The centrality or -1 in case of problems 
   */
  static Float_t GetCentralityCompat(const AliESDEvent& event,
				     const TString&     method,
				     Int_t&             qual,
				     Bool_t             verbose);
  //____________________________________________________________________
  /** 
   * Compatbility function to retrieve the centrality from the AOD header. 
   * 
   * @param event   Event 
   * @param method  Centrality method 
   * @param qual    Quality flat - 0 is good 
   * @param verbose If true, be verbose 
   * 
   * @return The centrality or -1 in case of problems 
   */
  static Float_t GetCentralityCompat(const AliAODEvent& event,
				     const TString&     method,
				     Int_t&             qual,
				     Bool_t             verbose);
  
  //==================================================================
  /** 
   * @{ 
   * @name Printing utilities 
   */
  /** 
   * Print a task 
   * 
   * @param o The task to print 
   */
  static void PrintTask(const TObject& o);
  /** 
   * Print a name 
   * 
   * @param name The name to pring 
   */
  static void PrintName(const char* name);
  /** 
   * Print a field 
   * 
   * @param name  The name of the field
   * @param value Format for the value 
   */
  static void PrintField(const char* name, const char* value, ...);
  /* @} */

  //==================================================================
  /** 
   * @{
   * @name Convenience containers 
   */
  /** 
   * Structure to hold histograms 
   *
   * @ingroup pwglf_forward 
   */
  struct Histos : public TObject
  {	
    /** 
     * Constructor 
     */
    Histos() : fFMD1i(0), fFMD2i(0), fFMD2o(0), fFMD3i(0), fFMD3o(0) {}
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    Histos(const Histos& o) 
      : TObject(o), 
	fFMD1i(o.fFMD1i), 
	fFMD2i(o.fFMD2i), 
	fFMD2o(o.fFMD2o), 
	fFMD3i(o.fFMD3i), 
	fFMD3o(o.fFMD3o)
    {}
    /** 
     * Assignement operator 
     * 
     * @return Reference to this 
     */
    Histos& operator=(const Histos&) { return *this;}
    /** 
     * Destructor.  This does not delete the interally allocated
     * memory.  Use the member function Delete for that.
     */
    ~Histos();
    /** 
     * Clear internal memory.  Note, if the internal histograms are
     * added to an output container, then we must not free this
     * memory.
     */
    void Delete(Option_t* opt="");
    /** 
     * Initialize the object 
     * 
     * @param etaAxis Eta axis to use 
     */
    void Init(const TAxis& etaAxis);
    /** 
     * Re-initialize the object with new @f$\eta@f$ axis 
     * 
     * @param etaAxis Eta axis to use 
     */
    void ReInit(const TAxis& etaAxis);
    /** 
     * Make a histogram 
     * 
     * @param d        Detector
     * @param r        Ring 
     * @param etaAxis  Eta axis to use
     * 
     * @return Newly allocated histogram 
     */
    static TH2D* Make(UShort_t d, Char_t r, const TAxis& etaAxis);
    /** 
     * Set the @f$\eta@f$ axis 
     * 
     * @param hist    Histogram
     * @param etaAxis @f$\eta@f$ axis to use
     */
    static void RebinEta(TH2D* hist, const TAxis& etaAxis);
    /** 
     * Clear data 
     * 
     * @param option Not used 
     */
    void  Clear(Option_t* option="");
    /** 
     * Get the histogram for a particular detector,ring
     * 
     * @param d Detector 
     * @param r Ring 
     * 
     * @return Histogram for detector,ring or nul 
     */
    TH2D* Get(UShort_t d, Char_t r) const;
    TH2D* fFMD1i; // Histogram for FMD1i
    TH2D* fFMD2i; // Histogram for FMD2i
    TH2D* fFMD2o; // Histogram for FMD2o
    TH2D* fFMD3i; // Histogram for FMD3i
    TH2D* fFMD3o; // Histogram for FMD3o

    ClassDef(Histos,2);
  };

  //__________________________________________________________________
  /**
   * Base class for structure holding ring specific histograms
   * 
   * @ingroup pwglf_forward 
   */
  struct RingHistos : public TObject
  {
    /** 
     * Constructor
     */
    RingHistos() : fDet(0), fRing('\0'), fName(""), fkNSector(0), fkNStrip(0) {}
    /** 
     * Constructor
     * 
     * @param d Detector
     * @param r Ring 
     */
    RingHistos(UShort_t d, Char_t r) 
      : fDet(d), fRing(r), fName(TString::Format("FMD%d%c", d, r)),
	fkNSector(r == 'i' || r == 'I' ? 20 : 40), 
	fkNStrip(r == 'i' || r == 'I' ? 512 : 256)
    {}
    /** 
     * Copy constructor
     * 
     * @param o Object to copy from 
     */
    RingHistos(const RingHistos& o) 
      : TObject(o), fDet(o.fDet), fRing(o.fRing), fName(o.fName),
	fkNSector(o.fkNSector), fkNStrip(o.fkNStrip)
    {}
    /** 
     * 
     */
    virtual ~RingHistos() {}
    /** 
     * Assignement operator
     * 
     * @param o Object to assign from
     * 
     * @return Reference to this
     */
    RingHistos& operator=(const RingHistos& o) 
    {
      if (&o == this) return *this;
      TObject::operator=(o);
      fDet  = o.fDet;
      fRing = o.fRing;
      fName = o.fName;
      fkNSector = o.fkNSector;
      fkNStrip  = o.fkNStrip;
      return *this;
    }
    /** 
     * Define the outout list in @a d
     * 
     * @param d Where to put the output list
     * 
     * @return Newly allocated TList object or null
     */
    TList* DefineOutputList(TList* d) const;
    /** 
     * Get our output list from the container @a d
     * 
     * @param d where to get the output list from 
     * 
     * @return The found TList or null
     */
    TList* GetOutputList(const TList* d) const;
    /** 
     * Find a specific histogram in the source list @a d
     * 
     * @param d     (top)-container 
     * @param name  Name of histogram
     * 
     * @return Found histogram or null
     */
    TH1* GetOutputHist(const TList* d, const char* name) const;
    /** 
     * Get the colour of this ring 
     * 
     * @return 
     */
    Color_t Color() const 
    { 
      return AliForwardUtil::RingColor(fDet, fRing);
    }
    /** 
     * The name of this ring 
     * 
     * @return Name of this ring 
     */
    const char* GetName() const { return fName.Data(); } 
    /** 
     * Get number of sectors 
     */
    const UShort_t& NSector() const { return fkNSector; }
    /** 
     * Get number of strips 
     */
    const UShort_t& NStrip() const { return fkNStrip; }
    UShort_t fDet;   // Detector
    Char_t   fRing;  // Ring
    TString  fName;  // Name
    UShort_t fkNSector; // Number of sectors 
    UShort_t fkNStrip;  // Number of strips 

    ClassDef(RingHistos,1);
  };
  /* @} */

  //__________________________________________________________________
  /**
   * A guard idom for producing debug output 
   * 
   */
  struct DebugGuard 
  {
    /** 
     * Constructor 
     * 
     * @param lvl       Current level
     * @param msgLvl    Target level 
     * @param format    @c printf -like format
     * 
     * @return 
     */
    DebugGuard(Int_t lvl, Int_t msgLvl, const char* format, ...);
    /** 
     * Destructor
     */
    ~DebugGuard();
    /** 
     * Make a message 
     * 
     * @param lvl    Current level
     * @param msgLvl Target level 
     * @param format @c printf -like format
     */
    static void Message(Int_t lvl, Int_t msgLvl, const char* format, ...);
  private:
    /** 
     * Output the message 
     * 
     * @param in    Direction
     * @param msg   Message 
     */
    static void Output(int in, TString& msg);
    /** 
     * Format a message 
     * 
     * @param out     Output is stored here
     * @param format  @c printf -like format
     * @param ap      List of arguments
     */
    static void Format(TString& out, const char* format, va_list ap);
    TString fMsg;
  };
  //__________________________________________________________________
   /** 
   * A guard to suppress messages 
   */
  struct SuppressGuard
  {
    /** The previous message level */
    Int_t save;
    /** 
     * Constructor 
     * 
     * @param lvl Level to suppress to 
     */
    SuppressGuard(Int_t lvl=2000);
    /** 
     * Destructor 
     */
    ~SuppressGuard();
  };
private:
  /** 
   * Constructor 
   */
  AliForwardUtil() {}
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardUtil(const AliForwardUtil& o) : TObject(o) {}
  /** 
   * Assingment operator 
   * 
   * 
   * @return Reference to this object 
   */
  AliForwardUtil& operator=(const AliForwardUtil&) { return *this; }
  /** 
   * Destructor
   */
  ~AliForwardUtil() {}
  

  ClassDef(AliForwardUtil,1);// Utilities - do not make object
};

// #ifdef LOG_NO_DEBUG
// # define DGUARD(L,N,F,...) do {} while(false) 
// #else
/** 
 * Macro to declare a DebugGuard
 * 
 * @param L Current debug level
 * @param N Target debug level
 * @param F @c printf -like Format 
 */
# define DGUARD(L,N,F,...)					\
  AliForwardUtil::DebugGuard _GUARD(L,N,F, ## __VA_ARGS__)
/** 
 * Macro to make a debug message, using DebugGuard::Message
 * 
 * @param L Current debug level
 * @param N Target debug level
 * @param F @c printf -like Format 
 */
# define DMSG(L,N,F,...)					\
  AliForwardUtil::DebugGuard::Message(L,N,F, ## __VA_ARGS__)
// #endif
#endif
// Local Variables:
//  mode: C++
// End:

