#ifndef ALIFMDPOINTS_H
#define ALIFMDPOINTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDPoints.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Tue Apr 11 12:35:31 2006
    @brief   Specialised class for drawing hits in the FMD. 
    @ingroup FMD_sim
*/
#ifndef ALIPOINTS_H
# include <AliPoints.h>
#endif
class AliFMHit;
class TMarker3DBox;

/** @class AliFMDPoints 
    @brief Class to draw hits for the FMD 
    @ingroup FMD_sim
 */
class AliFMDPoints : public AliPoints 
{
public:
  /** Constructor 
      @param hit   Hit to draw
      @param color Color of hit */
  AliFMDPoints(AliFMDHit* hit, UInt_t color);
  /** Copy constructor 
      @param other Objct to copy from  */
  AliFMDPoints(const AliFMDPoints& other);
  /** Assignment operator
      @param other Objct to copy from  
      @return reference to this */
  AliFMDPoints& operator=(const AliFMDPoints& other);
  /** Destructor  */
  virtual ~AliFMDPoints();
  /** Set position 
      @param x @f$ x@f$ coordinate 
      @param y @f$ y@f$ coordinate
      @param z @f$ z@f$ coordinate */
  void SetXYZ(Double_t x, Double_t y, Double_t z);
  /** Calculate the distance to this object from cursor
      @param px Curser X-position
      @param py Cursor Y-position
      @return Distance to @f$(p_x,p_y)@f$ */
  Int_t DistancetoPrimitive(Int_t px, Int_t py);
  /** Attach to canvas 
      @param option See TMarker3DBox::Draw */
  void Draw(Option_t* option);
  /** Paint in canvas 
      @param option See TMarker3DBox::Paint */
  void Paint(Option_t* option);
  /** @param colour Colour of marker */
  void SetMarkerColor(Color_t colour);
private:
  /** Marker */
  TMarker3DBox* fMarker; // Marker 
};

#endif
// Local Variables:
//   mode: C++ 
// End:
//
// EOF
//
