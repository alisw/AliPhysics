#ifndef ALIFMDSURVEYTOALIGNOBJS_H
#define ALIFMDSURVEYTOALIGNOBJS_H
#include <AliSurveyToAlignObjs.h>
#include <TGeoMatrix.h>

// Forward decl
class TVector3;
class TGeoMatrix;


/**
 * Class to take survey data and transform that to alignment objects. 
 * 
 */
class AliFMDSurveyToAlignObjs : public AliSurveyToAlignObjs
{
public:
  /** 
   * Constructor
   * 
   */
  AliFMDSurveyToAlignObjs() : AliSurveyToAlignObjs(),
			      fFMD1Delta(0),
			      fFMD2Delta(0) {}
  /** 
   * Run the task.
   * 
   */  
  void Run();
  /** 
   * 
   * Method to create the alignment objects
   * 
   * @return @c true on success, @c false otherwise
   */  
  Bool_t CreateAlignObjs();

  
  TClonesArray* GetAlignObjArray() const { return fAlignObjArray; }
  
protected:
  /** 
   * Do the FMD1 analysis.  We have 4 survey targets on V0-A on the
   * C-side.  These are 
   *
   *  - V0A_ICT  In-side, C-side, top.
   *  - V0A_ICB  In-side, C-side, bottom.  
   *  - V0A_OCT  Out-side, C-side, top.	 
   *  - V0A_OCB	 Out-side, C-side, bottom.
   * 
   * These 4 survey targets sit 3.3mm over the V0-A C-side surface, or
   * 3.3mm over the back surface of FMD1.  
   *
   * Since these are really sitting on a plane, we can use the method
   * proposed by the CORE offline. 
   * 
   * @return @c true on success, @c false otherwise.
   */
  Bool_t DoFMD1();
  /** 
   * Get the FMD1 plane from the survey points
   * 
   * @param rot    Rotation matrix (direction cosines)
   * @param trans  Translation
   * 
   * @return @c true on success, @c false otherwise.
   */
  Bool_t GetFMD1Plane(Double_t* rot, Double_t* trans) const;
  /** 
   * Do the FMD2 calculations.  We have 6 survey points of which only
   * 5 are normally surveyed.  These are all sittings 
   *
   *  - FMD2_ITOP   - In-side, top
   *  - FMD2_IBOTM  - In-side, middle bottom
   *  - FMD2_IBOT   - In-side, bottom
   *  - FMD2_OTOP   - Out-side, top
   *  - FMD2_OBOTM  - Out-side, middle bottom
   *  - FMD2_OBOT   - Out-side, bottom
   *
   * The nominal coordinates of these retro-fitted survey stickers
   * isn't known.  Also, these stickers are put on a thin (0.3mm
   * thick) carbon cover which flexes quite easily.  This means, that
   * to rotations and xy-translation obtained from the survey data
   * cannot be used, and left is only the z-translation.
   *
   * Further more, since FMD2 to is attached to the ITS SPD thermal
   * screen, it is questionable if the FMD2 survey will ever be used. 
   * 
   * @return @c true on success, @c false otherwise.
   */
  Bool_t DoFMD2();
  /** 
   * Get the surveyed plane corresponding to the backside of FMD2.
   * The plane is done as a best fit of the plane equation to at least
   * 4 of the available survey points.
   * 
   * @param rot    Rotation matrix (direction cosines)
   * @param trans  Translation vector.
   * 
   * @return @c true on success, @c false otherwise
   */
  Bool_t GetFMD2Plane(Double_t* rot, Double_t* trans) const;
  /** 
   * Get the factor to translate current coordinates to the canonical
   * unit (centi-meters). 
   * 
   * @return Conversion factor
   */
  Double_t GetUnitFactor() const;
  /** 
   * Get the coordinates of a survey point (if available).
   * 
   * @param name Name of the survey point.
   * @param p    Coordinates.
   * @param e    Error on the measurement.
   * 
   * @return @c true if the survey data is available, @c false otherwise.
   */
  Bool_t   GetPoint(const char* name, TVector3& p, TVector3& e) const;
  /** 
   * Calculate the plane translation and rotation from 3 survey points
   * 
   * @param a     1st Survey point 
   * @param b     2nd Survey point
   * @param c     3rd Survey point
   * @param trans Translation vector
   * @param rot   Rotation matrix (direction cosines)
   * 
   * @return 
   */
  Bool_t   CalculatePlane(const     TVector3& a, 
			  const     TVector3& b,
			  const     TVector3& c,
			  Double_t  depth,
			  Double_t* trans,
			  Double_t* rot) const;
  /** 
   * Calculate the plane rotation and translation by doing a fit of
   * the plane equation to the surveyed points.  At least 4 points
   * must be passed in the @a points array with corresponding errors
   * in the array @a errors.  The arrays are assumed to contain
   * TVector3 objects.
   * 
   * @param points Array surveyed positions
   * @param errors Array of errors corresponding to @a points
   * @param depth  Survey targets depth (perpendicular to the plane)
   * @param trans  On return, translation of the plane
   * @param rot    On return, rotation (direction cosines) of the plane
   * 
   * @return @c true on success, @c false otherwise
   */
  Bool_t FitPlane(const TObjArray& points, 
		  const TObjArray& errors,
		  Double_t         depth,
		  Double_t*        trans,
		  Double_t*        rot) const;
  /** 
   * Create a delta transform from a global rotation matrix and
   * translation. 
   * 
   * @param global Global matrix of element to transform.
   * @param rot    Rotation matrix (direction cosines)
   * @param trans  Translation 
   * @param delta  On return, the delta transform
   * 
   * @return Newly 
   */
  Bool_t MakeDelta(const TGeoMatrix*  global,
		   const Double_t*    rot, 
		   const Double_t*    trans,
		   TGeoHMatrix& delta) const;
  /** 
   * Create a delta transform from a global rotation matrix and
   * translation. 
   * 
   * @param path   Path of element to transform.
   * @param rot    Rotation matrix (direction cosines)
   * @param trans  Translation 
   * @param delta  On return, the delta transform
   * 
   * @return Newly 
   */
  Bool_t MakeDelta(const char*  path, 
		   const Double_t*    rot, 
		   const Double_t*    trans,
		   TGeoHMatrix& delta) const;
  /** 
   * Service member function to print a vector
   * 
   * @param text Prefix text
   * @param v    Vector (array of 3 doubles)
   */
  static void PrintVector(const char* text, const Double_t* v);
  /** 
   * Service member function to print a vector
   * 
   * @param text Prefix text
   * @param v    Vector
   */
  static void PrintVector(const char* text, const TVector3& v);
  /** 
   * Service member function to print a rotation matrix
   * 
   * @param text Prefix text
   * @param v    Matrix (array of 9 doubles)
   */
  static void PrintRotation(const char* text, const Double_t* rot);

  TGeoHMatrix fFMD1Delta; // FMD1 delta transform
  TGeoHMatrix fFMD2Delta; // FMD2 delta transform 
  
  ClassDef(AliFMDSurveyToAlignObjs,0) // Convert FMD survey to alignments
};


#endif
//____________________________________________________________________
//
// Local Variables:
//  mode: C++
// End:
//

