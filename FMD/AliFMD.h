#ifndef ALIFMD_H
#define ALIFMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMD.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 17:59:37 2006
    @brief   Declaration of AliFMD detector driver 
*/
/** @mainpage ALICE FMD Off-line code 
    
    @b Contents 
    - @ref intro 
    - @ref structure 
      - @ref base  (see also @ref FMD_base)
      - @ref sim   (see also @ref FMD_sim)
      - @ref rec   (see also @ref FMD_rec)
      - @ref util  (see also @ref FMD_util)
    - @ref script  (see also @ref FMD_script)
    - @ref quick
    - @ref authors
    
    @section intro Introduction:

    This file contains a short overview of the FMD code.   It is by no 
    means authoritative  - only the code is good for that.  However,
    I'll try to explain things a bit here. 

    @section structure Structure:

    There are 4 libraries build for the FMD.  These are 

    - libFMDbase: This contains the basic stuff, like data classes and
      managers.
    - libFMDsim: This contains code used by the simulation only.  That
      includes the detector class AliFMD and it's derivatives.  It
      also contains the code for building the geometry, as well as the
      digitisers and raw data writer.
    - libFMDrec: Code needed for the reconstruction.  This include the
      reconstruction code itself, as well as a raw data reader.
    - libFMDutil: This is a special library that contains various
      utility classes for the FMD expert/developer/tester.  It
      includes code to read all data produced by the FMD, a simple
      event display, and code to make fake calibration and alignment
      data.  This library is @e not loaded by aliroot
      automatically. The user has to load it herself:
      @code 
      gSystem->Load("libFMDutil.so");
      @endcode
    The content of these libraries is detailed more below. 
    
    @subsection base libFMDbase:

    This currently (18th or March 2006) contains the classes 

    - AliFMDBaseDigit, AliFMDDigit, AliFMDSDigit: Base class for
      digits, real digits, and summable digits.  The are all data
      classes that holds the ADC value(s) for a single strip.
    - AliFMDBoolMap: A map of booleans, one for each strip.
    - AliFMDUShortMap: A map of unsigned short integers, one for each
      strip.
    - AliFMDDetector, AliFMD1, AliFMD2, AliFMD3: 3 classes that holds
      the parameters for each of the 3 FMD sub-detectors.  These
      derive from AliFMDDetector, and are managed by AliFMDGeometry.
      Each of these also contains the translation from sensor
      reference frame to global reference frame.
    - AliFMDRing: Parameters for the FMD rings (inner and outer type).
      These are managed by AliFMDGeometry.
    - AliFMDGeometry: Manager of FMD geometry data. All code queries
      this manager for geometry parameters, so that the data is always
      consistent.  
    - AliFMDParameters: Manager of FMD parameters, like calibration
      parameters.  This class fetches data from CDB if necessary.
      All code queries this manager for parameters, so that the data
      is always consistent. 
    - AliFMDCalibPedestal, AliFMDCalibGain, AliFMDCalibSampleRate,
      AliFMDAltroMapping: Containers of calibration parameters.  These
      correspond to the pedestal and its width, the gain of each
      strip, the oversampling rate used in read-out, and the map from
      hardware address to detector.
    - AliFMDAltroIO, AliFMDAltroReader, AliFMDAltroWriter: Low-level
      classes to do input/output on ALTRO formated buffers. 

  

    @subsection sim libFMDsim:

    This currently (18th or March 2006) contains the classes 

    - AliFMDEdepMap: Cache of energy deposited and total number of
      hits for each strip.  The digitiser AliFMDDigitizer uses this to
      store simulation data before making digits.
    - AliFMDHit: A hit in the FMD active elements, as described by the
      simulation back-end during transport.
    - AliFMD, AliFMDv0, AliFMDv1: Simulation drivers for the FMD.
      AliFMD is the base class. AliFMDv0 corresponds to a simulation
      where no hits are created, but the material distribution is
      right.  AliFMDv1 is like AliFMDv0, except that hits are
      produced.
    - AliFMDGeometryBuilder: Build the FMD geometry in terms of TGeo
      objects.  The information for building the geometry is retrieved
      from AliFMDGeometry.
    - AliFMDBaseDigitizer, AliFMDDigitizer, AliFMDSDigitizer: Base
      class for the digitisers.  AliFMDDigitizer makes `real' digits
      (AliFMDDigit) from hits, and AliFMDSDigitizer makes summable
      digits from hits.
    - AliFMDRawWriter: Writes a pseudo raw data file from the digits
      created by the digitisers.  It uses the AliFMDAltroMapping from
      AliFMDParameters to make the mapping from detector coordinates
      to hardware addresses.

    @subsection rec libFMDrec:

    This currently (18th or March 2006) contains the classes 

    - AliFMDReconstructor: Reconstruct (in a naiive way) the charged
      particle multiplicity in the FMD strips.  This also writes an
      AliESDFMD object to the ESD files (that class is in libESD).

    - AliFMDRecPoint: Reconstructed point in the FMD.  These objects
      are made AliFMDReconstructor. 

    - AliFMDRawReader: Classes to read raw data files. 

    @subsection util libFMDutil:

    This currently (18th or March 2006) contains the classes 

    - AliFMDInput, AliFMDInputHits, AliFMDInputDigits,
      AliFMDInputSDigits, AliFMDInputRecPoints: Base class, and
      concrete classes to read in FMD generated data.  These provide a
      simple and unified way of getting the data.  Hooks are defined
      to process hits, tracks, digits, and reconstructed points, as
      well as geometry and ESD data.   See for example the scripts
      @c DrawHits.C, @c DrawHitsDigits.C, @c DrawHitsRecs.C, @c
      DrawDigitsRecs.C in the @c FMD/scripts sub-directory.

    - AliFMDDisplay: Simple event display for FMD data, including
      hits, digits, reconstructed points and ESD data. 

    - AliFMDCalibFaker, AliFMDAlignFaker: Classes to write fake (or
      dummy) calibration and alignment 	data.  These derive from
      TTask.  

    @section script Scripts 
    
    Most scripts live in @c FMD/scripts.  The notiable exceptions are 
    @ref Simulate.C, @ref Reconstruct.C, and @ref Config.C 

    @section quick Quick start 

    First, install ROOT.  Then Install TGeant3: 
    @verbatim 
    > cd ~/
    > mkdir alice
    > cd alice
    > cvs -d :pserver:cvs@root.cern.ch:/user/cvs login 
    Password: cvs
    > cvs -d :pserver:cvs@root.cern.ch:/user/cvs co geant3
    > cd geant3
    > make 
    @endverbatim 

    Now you can install AliRoot 
    @verbatim 
    > cd ../
    > cvs -d :pserver:cvs@alisoft.cern.ch:/soft/cvsroot login
    Password: <empty>
    > cvs -d :pserver:cvs@alisoft.cern.ch:/soft/cvsroot co AliRoot
    > cd AliRoot
    > export ALICE_TARGET=`root-config --arch`
    > export ALICE=${HOME}/alice
    > export ALICE_ROOT=${ALICE}/AliRoot
    > export ALICE_LEVEL=new
    > export LD_LIBRARY_PATH=${ALICE_ROOT}/lib/tgt_${ALICE_TERGET}:${LD_LIBRARY_PATH}
    > export PATH=${ALICE_ROOT}/bin/tgt_${ALICE_TERGET}:${PATH}
    > export G3SYS=${ALICE}/geant3
    > make 
    @endverbatim 
    
    To simulate one event, do something like 

    @verbatim 
    > aliroot ${ALICE_ROOT}/FMD/Simulate.C
    @endverbatim 

    To reconstruct the generated event, do 
    @verbatim 
    > aliroot ${ALICE_ROOT}/FMD/Reconstruct.C
    @endverbatim 

    Now, open the file `AliESDs.root' in AliRoot, and browse through  that. 

    @section authors Authors:

    - Alla Maevskaya		<Alla.Maevskaia@cern.ch>	
    - Christian Holm Christensen 	<cholm@nbi.dk>
*/
/** @defgroup FMD_sim Simulation */

//____________________________________________________________________
//
//  Manager class for the FMD - Base class.
//  AliFMDv1, AliFMDv0, and AliFMDAlla 
//  provides concrete implementations. 
//  This class is sooooo crowded
//
#ifndef ALIDETECTOR_H  
# include <AliDetector.h>
#endif
class TBranch;
class TClonesArray;
class TBrowser;
class TMarker3DBox;
class AliDigitizer;
class AliFMDHit;

//____________________________________________________________________
/** @class AliFMD AliFMD.h <FMD/AliFMD.h>
    @brief Forward Multiplicity Detector based on Silicon wafers. 
    This class is the driver for especially simulation. 

    The Forward Multiplicity Detector consists of 3 sub-detectors FMD1,
    FMD2, and FMD3, each of which has 1 or 2 rings of silicon sensors. 
                                                          
    This is the base class for all FMD manager classes. 
                       
    The actual code is done by various separate classes.   Below is
    diagram showing the relationship between the various FMD classes
    that handles the simulation

    @verbatim
          +----------+   +----------+   
          | AliFMDv1 |   | AliFMDv0 |   
          +----------+   +----------+   
               |              |                    +-----------------+
          +----+--------------+                 +--| AliFMDDigitizer |
          |                                     |  +-----------------+
          |           +---------------------+   |
          |        +--| AliFMDBaseDigitizer |<--+
          V     1  |  +---------------------+   |
     +--------+<>--+                            |  +------------------+
     | AliFMD |                                 +--| AliFMDSDigitizer |    
     +--------+<>--+                               +------------------+       
	       1  |  +---------------------+
	          +--| AliFMDReconstructor |
	  	     +---------------------+
    @endverbatim
		     
    -  AliFMD 
       This defines the interface for the various parts of AliROOT that
       uses the FMD, like AliFMDSimulator, AliFMDDigitizer, 
       AliFMDReconstructor, and so on. 
    -  AliFMDv0
       This is a concrete implementation of the AliFMD interface. 
       It is the responsibility of this class to create the FMD
       geometry.
    -  AliFMDv1 
       This is a concrete implementation of the AliFMD interface. 
       It is the responsibility of this class to create the FMD
       geometry, process hits in the FMD, and serve hits and digits to
       the various clients. 
    -  AliFMDSimulator
       This is the base class for the FMD simulation tasks.   The
       simulator tasks are responsible to implment the geoemtry, and
       process hits. 
    -  AliFMDReconstructor
       This is a concrete implementation of the AliReconstructor that
       reconstructs pseudo-inclusive-multiplicities from digits (raw or
       from simulation)

    Calibration and geometry parameters are managed by separate
    singleton managers.  These are AliFMDGeometry and
    AliFMDParameters.  Please refer to these classes for more
    information on these.    
 */
class AliFMD : public AliDetector 
{
public:
  /** Default constructor.  Do not use. */
  AliFMD();
  /** Normal constructor 
      @param name  Name of object.
      @param title Title of object. */
  AliFMD(const char *name, const char *title);
  /** Destructor */
  virtual ~AliFMD(); 
  /** Wheter to make a detailed geometry
      @param use If true, make detailed geometry  */
  void UseDetailed(Bool_t use=kTRUE) { fDetailed = use; }
  
  /** @{*/
  /** @name GEometry ANd Tracking (GEANT :-) */
  /** Define the geometry.  This is done by asking the manager
      AliFMDGeometry to construct the geometry.  This in turn calls
      AliFMDGeometryBuilder.   */
  virtual void   CreateGeometry();
  /** Create entries for alignable volumes associating the symbolic volume
      name with the corresponding volume path. Needs to be syncronized with
      eventual changes in the geometry.   */
  virtual void  AddAlignableVolumes() const;
  /** Create the tracking mediums used by the FMD.  This associates
      the tracking mediums defined with the FMD in the
      TVirtualMCApplication (AliMC). 
      
      The defined mediums are 
      -	@c FMD @c Si$	Silicon (active medium in sensors)
      -	@c FMD @c C$	Carbon fibre (support cone for FMD3 and vacuum pipe)
      -	@c FMD @c Al$	Aluminium (honeycomb support plates)
      -	@c FMD @c PCB$	Printed Circuit Board (FEE board with VA1_3)
      -	@c FMD @c Chip$	Electronics chips (currently not used)
      -	@c FMD @c Air$	Air (Air in the FMD)
      -	@c FMD @c Plastic$ Plastic (Support legs for the hybrid cards)
  */
  virtual void   CreateMaterials(); 
  /** Initialize this detector */
  virtual void   Init();
  /** This member function is called when ever a track deposites
      energy (or similar) in an FMD tracking medium.  In this base
      class this member function is pure abstract.   In concrete
      sub-classes, the member function may make hits or other
      stuff. */ 
  virtual void   StepManager() = 0;
  /** Called at the end of each simulation event.  If the debug level
      is high enough a list of @e bad hits is printed. */
  virtual void   FinishEvent();
  /** @}*/
  
  /** @{*/
  /** @name Graphics and event display */
  /** Build simple ROOT TNode geometry for event display. With the new
      geometry modeller, TGeoManager, this seems rather redundant. */
  virtual        void   BuildGeometry();
  /** Draw a shaded view of the Forward multiplicity detector.  This 
      isn't really useful anymore. */
  virtual        void   DrawDetector();
  /** Calculate the distance from the mouse to the FMD on the screen
      Dummy routine */
  virtual        Int_t  DistancetoPrimitive(Int_t px, Int_t py);
  /** Store x, y, z of all hits in memory for display. 
      Normally, the hits are drawn using TPolyMarker3D - however, that
      is not very useful for the FMD.  Therefor, this member function
      is overloaded to make TMarker3D, via the class AliFMDPoints.
      AliFMDPoints is a local class. 
      @param track the track number to load the hits for */
  virtual        void   LoadPoints(Int_t track);
  /** @}*/
  
  /** @{ */
  /** @name Hit and digit management */
  /* Create Tree branches for the FMD.
     @param opt  One of 
     - @c H Make a branch of TClonesArray of AliFMDHit's
     - @c D Make a branch of TClonesArray of AliFMDDigit's
     - @c S Make a branch of TClonesArray of AliFMDSDigit's */
  virtual void          MakeBranch(Option_t *opt=" ");
  /** Set the TClonesArray to read hits into.
      @param b The branch to containn the hits */
  virtual void          SetHitsAddressBranch(TBranch *b);
  /** Set branch address for the Hits, Digits, and SDigits Tree. */
  virtual void          SetTreeAddress();
  /** Get the array of summable digits
      @return summable digits */
  virtual TClonesArray* SDigits() { return fSDigits; }        
  /** Reset the array of summable digits */
  virtual void          ResetSDigits();
  /** Add a hit to the hits tree 
      @param  track  Track #
      @param  vol Volume parameters, interpreted as 
      - ivol[0]  [UShort_t ] Detector # 
      - ivol[1]	 [Char_t   ] Ring ID 
      - ivol[2]	 [UShort_t ] Sector #
      - ivol[3]	 [UShort_t ] Strip # 
      @param hits Hit information 
      - hits[0]	 [Float_t  ] Track's X-coordinate at hit 
      - hits[1]	 [Float_t  ] Track's Y-coordinate at hit
      - hits[3]  [Float_t  ] Track's Z-coordinate at hit
      - hits[4]  [Float_t  ] X-component of track's momentum 	       	 
      - hits[5]	 [Float_t  ] Y-component of track's momentum	       	 
      - hits[6]	 [Float_t  ] Z-component of track's momentum	       	
      - hits[7]	 [Float_t  ] Energy deposited by track		       	
      - hits[8]	 [Int_t    ] Track's particle Id # 
      - hits[9]	 [Float_t  ] Time when the track hit */
  virtual void          AddHit(Int_t track, Int_t *vol, Float_t *hits);
  /** Add a hit to the list
      @param track     Track #
      @param detector  Detector # (1, 2, or 3)                      
      @param ring      Ring ID ('I' or 'O')
      @param sector    Sector # (For inner/outer rings: 0-19/0-39)
      @param strip     Strip # (For inner/outer rings: 0-511/0-255)
      @param x	       Track's X-coordinate at hit
      @param y	       Track's Y-coordinate at hit
      @param z	       Track's Z-coordinate at hit
      @param px	       X-component of track's momentum 
      @param py	       Y-component of track's momentum
      @param pz	       Z-component of track's momentum
      @param edep      Energy deposited by track
      @param pdg       Track's particle Id #
      @param t	       Time when the track hit 
      @param len       Track length through the material. 
      @param stopped   Whether track was stopped or disappeared */
  virtual AliFMDHit*    AddHitByFields(Int_t    track, 
				       UShort_t detector, 
				       Char_t   ring, 
				       UShort_t sector, 
				       UShort_t strip, 
				       Float_t  x=0,
				       Float_t  y=0, 
				       Float_t  z=0,
				       Float_t  px=0, 
				       Float_t  py=0, 
				       Float_t  pz=0,
				       Float_t  edep=0,
				       Int_t    pdg=0,
				       Float_t  t=0, 
				       Float_t  len=0, 
				       Bool_t   stopped=kFALSE);
  /** Add a digit to the Digit tree 
      @param digits
      - digits[0]  [UShort_t] Detector #
      - digits[1]  [Char_t]   Ring ID
      - digits[2]  [UShort_t] Sector #
      - digits[3]  [UShort_t] Strip #
      - digits[4]  [UShort_t] ADC Count 
      - digits[5]  [Short_t]  ADC Count, -1 if not used
      - digits[6]  [Short_t]  ADC Count, -1 if not used 
      @param notused Not used */
  virtual        void   AddDigit(Int_t *digits, Int_t* notused=0);
  /** add a real digit
      @param detector  Detector # (1, 2, or 3)                      
      @param ring	  Ring ID ('I' or 'O')
      @param sector	  Sector # (For inner/outer rings: 0-19/0-39)
      @param strip	  Strip # (For inner/outer rings: 0-511/0-255)
      @param count1    ADC count (a 10-bit word)
      @param count2    ADC count (a 10-bit word), or -1 if not used
      @param count3    ADC count (a 10-bit word), or -1 if not used */
  virtual        void   AddDigitByFields(UShort_t detector=0, 
					 Char_t   ring='\0', 
					 UShort_t sector=0, 
					 UShort_t strip=0, 
					 UShort_t count1=0, 
					 Short_t  count2=-1, 
					 Short_t  count3=-1);
  /** Add a digit to the Digit tree 
      @param digits
      - digits[0]  [UShort_t] Detector #
      - digits[1]  [Char_t]   Ring ID
      - digits[2]  [UShort_t] Sector #
      - digits[3]  [UShort_t] Strip #
      - digits[4]  [UShort_t] ADC Count 
      - digits[5]  [Short_t]  ADC Count, -1 if not used
      - digits[6]  [Short_t]  ADC Count, -1 if not used  */
  virtual        void   AddSDigit(Int_t *digits);
  /** add a summable digit - as coming from data
      @param detector  Detector # (1, 2, or 3)                      
      @param ring      Ring ID ('I' or 'O')
      @param sector    Sector # (For inner/outer rings: 0-19/0-39)
      @param strip     Strip # (For inner/outer rings: 0-511/0-255)
      @param edep      Energy deposited   
      @param count1    ADC count (a 10-bit word)
      @param count2    ADC count (a 10-bit word), or -1 if not used 
      @param count3    ADC count (a 10-bit word), or -1 if not used */
  virtual        void   AddSDigitByFields(UShort_t detector=0, 
					  Char_t   ring='\0', 
					  UShort_t sector=0, 
					  UShort_t strip=0, 
					  Float_t  edep=0,
					  UShort_t count1=0, 
					  Short_t  count2=-1, 
					  Short_t  count3=-1);
  /** @}*/

  /** @{ */
  /** @name Digitisation */
  /** Create a digitizer object
      @param manager Digitization manager
      @return a newly allocated AliFMDDigitizer */
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
  /** Create AliFMDDigit's from AliFMDHit's.  This is done by creating
      an AliFMDDigitizer object, and executing it.  */
  virtual        void   Hits2Digits();
  /** Create AliFMDSDigit's from AliFMDHit's.  This is done by creating
      an AliFMDSDigitizer object, and executing it.  */
  virtual        void   Hits2SDigits();
  /** @}*/

  /** @{ */
  /** @name Raw data */
  /** Turn digits into raw data. This uses the class AliFMDRawWriter
      to do the job.   Please refer to that class for more
      information. */
  virtual        void   Digits2Raw();
  /** @}*/

  /** @{ */
  /** @name Utility */
  /** Browse this object 
      @param b Browser to show this object in */
  void   Browse(TBrowser* b);
  /** @}*/
protected:
  /** Initialize hit array if not already done, and return pointert. 
      @return Hit array */
  TClonesArray*      HitsArray();
  /** Initialize digit array if not already done, and return pointert. 
      @return Digit array */
  TClonesArray*      DigitsArray();
  /** Initialize summable digit array if not already done, and return
      pointert.  
      @return Summable digit array */
  TClonesArray*      SDigitsArray();

  TClonesArray*      fSDigits;              // Summable digits
  Int_t              fNsdigits;             // Number of digits  
  Bool_t             fDetailed;             // Use detailed geometry
  Bool_t             fUseOld;               // Use old approx geometry
  Bool_t             fUseAssembly;          // Use divided volumes
  
  enum {
    kSiId,                 // ID index of Si medium
    kAirId,                // ID index of Air medium
    kPlasticId,            // ID index of Plastic medium
    kPcbId,                // ID index of PCB medium
    kSiChipId,             // ID index of Si Chip medium
    kAlId,                 // ID index of Al medium
    kCarbonId,             // ID index of Carbon medium
    kCopperId,             // ID index of Copper Medium
    kKaptonId              // ID index of Kapton Medium
  };  

  TObjArray*         fBad;                  //! debugging - bad hits 

private:  
  /** Copy constructor 
      @param other Object to copy from */
  AliFMD(const AliFMD& other);
  /** Assignment operator 
      @param other Object to assign from
      @return Reference to this object  */
  AliFMD& operator=(const AliFMD& other);

  ClassDef(AliFMD,11)     // Base class FMD entry point
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//
