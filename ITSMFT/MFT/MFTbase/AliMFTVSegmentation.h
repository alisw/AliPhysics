#ifndef AliMFTVSegmentation_H
#define AliMFTVSegmentation_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTVSegmentation
/// \brief Abstract base class for MFT Segmentation description
///
/// Units are cm and deg
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date June 9th, 2015

#include "TNamed.h"
#include "TVector3.h"
#include "TGeoMatrix.h"


class AliMFTVSegmentation : public TNamed {
  
public:
  
  AliMFTVSegmentation();
  AliMFTVSegmentation(const AliMFTVSegmentation& input);

  virtual ~AliMFTVSegmentation(){};
  
  /// Set Position of the Element. Unit is [cm]
  void SetPosition(const Double_t *pos){
    fTransformation->SetTranslation(pos[0],pos[1],pos[2]);
  };
  
  /// \brief Set The rotation angles. Unit is [deg].
  void SetRotationAngles(const Double_t *ang);
  
  /// \brief Rotate around X axis, ang in deg
  void RotateX(const Double_t ang) {fTransformation->RotateX(ang);};
  /// \brief Rotate around Y axis, ang in deg
  void RotateY(const Double_t ang) {fTransformation->RotateY(ang);};
  /// \brief Rotate around Z axis, ang in deg
  void RotateZ(const Double_t ang) {fTransformation->RotateZ(ang);};
  
  /// \brief Returns the Transformation Combining a Rotation followed by a Translation
  ///
  /// The rotation is a composition of : first a rotation about Z axis with
  /// angle phi, then a rotation with theta about the rotated X axis, and
  /// finally a rotation with psi about the new Z axis.
  /// [For more details see the ROOT TGeoCombiTrans documentation](https://root.cern.ch/root/htmldoc/TGeoCombiTrans.html).
  TGeoCombiTrans * GetTransformation() const {return fTransformation;};
  
private:

  TGeoCombiTrans * fTransformation; ///< \brief Represent a rotation folowed by a translation.
                                    /// The rotation is a composition of : first a rotation about Z axis with
                                    /// angle phi, then a rotation with theta about the rotated X axis, and
                                    /// finally a rotation with psi about the new Z axis.
  
  /// \cond CLASSIMP
  ClassDef(AliMFTVSegmentation, 1);
  /// \endcond

  
};


#endif

