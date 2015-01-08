#ifndef ALIFMDALIGNFAKER_H
#define ALIFMDALIGNFAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/** 
 * @file    AliFMDAlignFaker.h
 * @author  Christian Holm Christensen <cholm@nbi.dk>
 * @date    Sun Mar 26 18:26:39 2006
 * @brief   Make fake alignment 
 */
//____________________________________________________________________
//
//  Class 
//  to 
//  make 
//  fake 
//  alignment
//  parameters 
//
#ifndef ROOT_TTask
# include <TTask.h>
#endif
#ifndef ROOT_TVector3
# include <TVector3.h>
#endif
class TClonesArray;
class TString;

/** 
 * @class AliFMDAlignFaker 
 * @brief This task creates fake alignment. 
 * 
 * Which alignment, depends on the bit mask passed to the
 * constructor, or added by `AddAlign'.
 * 
 * The default is to write all alignment parameters to a local
 * storage @c local://cdb which is a directory in the current
 * directory.
 *
 * @ingroup FMD_util
 */
class AliFMDAlignFaker : public TTask
{
public:
  /** 
   * What to make alignments for 
   */
  enum EWhat {
    /** MAke alignment for sensors */
    kSensors =  1, 
    /** Make alignment for half-rings */
    kHalves
  };
  enum {
    /** All types of alignment */
    kAll             = (1<<kHalves|1<<kSensors)
  };
  /**
   * Constructor 
      @param mask Bit mask of what to make alignments for 
      @param geo  File to read geometry from, if needed
      @param loc  Where to store the result 
   */
  AliFMDAlignFaker(Int_t mask=kAll, 
		   const char* geo="geometry.root",
		   const char* loc="");
  /**
   * Destructor 
   */
  virtual ~AliFMDAlignFaker() {}
  /**
   * Add something to make alignment for 
   * @param w Bit of alignment mask 
   */
  void AddAlign(EWhat w) { SETBIT(fMask, w); }
  /**
   * Remove something to make alignment for 
   * @param w Bit of alignment mask 
   */
  void RemoveAlign(EWhat w) { CLRBIT(fMask, w); }
  /**
   * Set alignment select mask 
   * @param mask Bit mask 
   */
  void SetAlign(Int_t mask) { fMask = mask; }
  /**
   * Set the displacement (translation) of sensors.  The displacement
   * is selected random uniformly between the corresponding minimum
   * and maximum. 
   * @param x1 Minimum X displacement (in centimeters)
   * @param y1 Minimum Y displacement (in centimeters)
   * @param z1 Minimum Z displacement (in centimeters)
   * @param x2 Maximum X displacement (in centimeters)
   * @param y2 Maximum Y displacement (in centimeters)
   * @param z2 Maximum Z displacement (in centimeters) 
   */
  void SetSensorDisplacement(Double_t x1=0,   Double_t y1=0,   Double_t z1=0,
			     Double_t x2=.01, Double_t y2=.01, Double_t z2=0);
  /**
   * Set the rotation of sensors.  The displacement is selected
   * random uniformly between the corresponding minimum and maximum.
   * @param x1 Minimum X rotation (in degrees)
   * @param y1 Minimum Y rotation (in degrees)
   * @param z1 Minimum Z rotation (in degrees)
   * @param x2 Maximum X rotation (in degrees)
   * @param y2 Maximum Y rotation (in degrees)
   * @param z2 Maximum Z rotation (in degrees) 
   */
  void SetSensorRotation(Double_t x1=0,  Double_t y1=0,  Double_t z1=0,
			 Double_t x2=.5, Double_t y2=.5, Double_t z2=.5);
  /**
   * Set the displacement (translation) of half-rings.  The
   * displacement is selected random uniformly between the
   * corresponding minimum and maximum.
   * @param x1 Minimum X displacement
   * @param y1 Minimum Y displacement
   * @param z1 Minimum Z displacement
   * @param x2 Maximum X displacement
   * @param y2 Maximum Y displacement
   * @param z2 Maximum Z displacement 
   */
  void SetHalfDisplacement(Double_t x1=0,   Double_t y1=0,   Double_t z1=0,
			   Double_t x2=.05, Double_t y2=.05, Double_t z2=.05);
  /**
   * Set the rotation of half-rings.  The displacement is selected
   * random uniformly between the corresponding minimum and maximum.
   * @param x1 Minimum X rotation (in degrees)
   * @param y1 Minimum Y rotation (in degrees)
   * @param z1 Minimum Z rotation (in degrees)
   * @param x2 Maximum X rotation (in degrees)
   * @param y2 Maximum Y rotation (in degrees)
   * @param z2 Maximum Z rotation (in degrees) 
   */
  void SetHalfRotation(Double_t x1=0, Double_t y1=0, Double_t z1=0,
		       Double_t x2=0, Double_t y2=0, Double_t z2=0);
  /**
   * Set the output file name.  Should be a valid CDB URL.
   * @param file CDB URL 
   */
  void SetOutput(const char* file) { SetTitle(file); }
  /**
   * Set the file to read the geometry from. 
   * @param file File name 
   */
  void SetGeometryFile(const char* file) { SetName(file); }
  /**
   * Set the comment  
   */ 
  void SetComment(const Char_t* comment="dummy data") { fComment = comment; }
  /**
   * Make the alignment objects. 
   * @param option Not used. 
   */
  void Exec(Option_t* option="");
  /** 
   * Get the geometry
   * 
   * @param toCdb   Whether to store in CDB
   * @param storage Storage element to use 
   * 
   * @return true on success 
   */
  static Bool_t GetGeometry(Bool_t toCdb=kFALSE, 
			    const TString& storage=TString());
protected:
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDAlignFaker(const AliFMDAlignFaker& o) 
    : TTask(o), 
      fMask(0), 
      fSensorTransMin(0,0,0),
      fSensorTransMax(0,0,0),
      fSensorRotMin(0,0,0),
      fSensorRotMax(0,0,0),
      fHalfTransMin(0,0,0),
      fHalfTransMax(0,0,0),
      fHalfRotMin(0,0,0),
      fHalfRotMax(0,0,0),
      fRunMin(0), 
      fRunMax(0),
      fArray(0),
      fComment("")
  {}
  /** 
   * Assignment operator 
   * 
   * 
   * @return Reference to this
   */  
  AliFMDAlignFaker& operator=(const AliFMDAlignFaker&) { return *this; }
  
  /**
   * Make the alignment object for a path
   * @param path   Node path.
   * @param volID  Volume identifier 
   * @param transX Translation in X
   * @param transY Translation in Y
   * @param transZ Translation in Z
   * @param rotX   Rotation around X axis
   * @param rotY   Rotation around Y axis
   * @param rotZ   Rotation around Z axis
   * @return @c true on success 
   */
  Bool_t MakeAlign(const TString& path, Int_t volID, 
		   Double_t transX, Double_t transY, Double_t transZ,
		   Double_t rotX, Double_t rotY, Double_t rotZ);
  /**
   * Align a sensor 
   * @param path of a sensor 
   * @param id Volume id 
   */
  Bool_t MakeAlignSensor(const TString& path, Int_t id);
  /**
   * Align a half-ring
   * @param path of a sensor 
   * @param id Volume id 
   */
  Bool_t MakeAlignHalf(const TString& path, Int_t id);
  /**
   * Write to CDB 
   */
  void   WriteToCDB();
  /**
   * Write to file 
   */
  void   WriteToFile();
  Long_t        fMask;            // What to write 
  TVector3      fSensorTransMin;  // Minimum translations of a sensor
  TVector3      fSensorTransMax;  // Maximum translations of a sensor
  TVector3      fSensorRotMin;    // Minimum rotation of a sensor
  TVector3      fSensorRotMax;    // Maximum rotation of a sensor
  TVector3      fHalfTransMin;    // Minimum translations of a half-ring
  TVector3      fHalfTransMax;	  // Maximum translations of a half-ring
  TVector3      fHalfRotMin;	  // Minimum rotation of a half-ring    
  TVector3      fHalfRotMax;	  // Maximum rotation of a half-ring 
  Int_t         fRunMin;          // Run validity start 
  Int_t         fRunMax;          // Run validity end
  TClonesArray* fArray;           // Cache
  TString       fComment;         // Comment on data
  
  ClassDef(AliFMDAlignFaker,0)
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

