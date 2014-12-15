//
// $Id$
//
#ifndef ALIFMD3_H
#define ALIFMD3_H
/** @file    AliFMD3.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:23:48 2006
    @brief   Geometry parameters of the FMD3 detector. 
*/
// Geometry parameters of the FMD3 detector. FMD3 has a fairly
// complicated support structure.  The cone also supports the
// beam-pipe. 
// 
#ifndef ALIFMDDETECTOR_H
# include "AliFMDDetector.h"
#endif
#include <TObjArray.h>

/** @class AliFMD3 AliFMD3.h <FMD/AliFMD3.h> 
    @brief Geometry parameters of the FMD3 detector. 
    FMD3 has a fairly complicated support structure.  The cone also
    supports the beam-pipe.
    @image html FMD3.png 
    @ingroup FMD_base
*/
class AliFMD3 : public AliFMDDetector 
{
public: 
  /** Constructor 
      @param inner Pointer to inner ring description 
      @param outer Pointer to outer ring description */
  AliFMD3(AliFMDRing* inner, AliFMDRing* outer);
  /** Destructor */
  virtual ~AliFMD3(){}

  /** Initialize the geometry */
  virtual void Init();
  /** Get the Z offset (to inner ring) */
  Double_t GetNoseZ() const { return fNoseZ; }

  /** @return Z position of front of nose */
  Double_t GetFlangeDepth() const { return fFlangeDepth; }
  /** @return Nose inner radius */
  Double_t GetFlangeLength() const { return fFlangeLength; }
  /** @return Nose outer radius */
  Double_t GetFlangeWidth() const { return fFlangeWidth; }

  /** @return Length of nose in Z */
  Double_t GetFiducialRadius() const { return fFiducialRadius; }

  /** @return The angle of the cone on out-side */
  Double_t GetConeOuterAngle() const { return fConeOuterAngle; }
  /** @return The angle of the cone on out-side */
  Double_t GetConeInnerAngle() const { return fConeInnerAngle; }

  /** @return Hole off-set from nose */
  Double_t GetHoleOffset() const { return fHoleOffset; }
  /** @return Depth of holes  */
  Double_t GetHoleDepth() const { return fHoleDepth; }
  /** @return Length of holes  */
  Double_t GetHoleLength() const { return fHoleLength; }
  /** @return Lowest with of holes */
  Double_t GetHoleLowWidth() const { return fHoleLowWidth; }
  /** @return Highest width of holes */
  Double_t GetHoleHighWidth() const { return fHoleHighWidth; }

  /** @return Length of a bolt  */
  Double_t GetBoltLength() const { return fBoltLength; }
  /** @return Bolt radius  */
  Double_t GetBoltRadius() const { return fBoltRadius; }

  
  /** @return array of 3-vectors (z, r_low, r_high) of the cone
      radii. */
  const TObjArray& ConeRadii() const { return fConeRadii; }
  /** @return array of 2-vectors (x,y) of the fiducial holes in the
              flanges.  coordinates are in the global coordinate
              system. */
  const TObjArray& FiducialHoles() const { return fFiducialHoles; }
    
  /** Get the cone radii at @a z. 
      @param z Point to evaulate at 
      @param opt If @c "O" get the outer radii, if @c "I" get the
      inner radii. 
      @return the radius of the cone */
  Double_t ConeR(Double_t z, Option_t* opt="O") const;

protected: 
  Double_t	fNoseZ;			// Z position of front of nose
  Double_t      fFlangeDepth;           // Depth of flanges 
  Double_t      fFlangeHighR;           // Outer radius of flanges. 
  Double_t      fFlangeLength;          // Length of flanges 
  Double_t      fFlangeWidth;           // Width of flanges 

  Double_t      fFiducialRadius;        // Radius of fiducial holes.

  Double_t      fConeInnerAngle;        // Angle of cone on inside  
  Double_t      fConeOuterAngle;        // Angle of cone on outside  

  Double_t      fHoleOffset;            // Offset (from nose-tip) of
					// holes 
  Double_t      fHoleDepth;             // Depth of holes 
  Double_t      fHoleLength;            // Length of holes 
  Double_t      fHoleLowWidth;          // Lowest with of holes
  Double_t      fHoleHighWidth;         // Highest width of holes
  
  Double_t      fBoltLength;            // Length of a bolt 
  Double_t      fBoltRadius;            // Bolt radius 

  TObjArray     fConeRadii;             // Array of cone radii.
  TObjArray     fFiducialHoles;         // Array of fiducial hole (x,y)
  
  ClassDef(AliFMD3, 2);
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
//
// EOF
//
