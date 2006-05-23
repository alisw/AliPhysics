/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// Revision of includes 07/05/2004
//
/// \ingroup geometry
/// \class AliMUONGeometryConstituent
/// \brief Helper class for definititon of an assembly of volumes.
///
/// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_CONSTITUENT_H
#define ALI_MUON_GEOMETRY_CONSTITUENT_H

#include <TNamed.h>

class TGeoTranslation;
class TGeoRotation;
class TGeoCombiTrans;
class TObjArray;

class AliMUONGeometryConstituent : public TNamed
{
  public:
    AliMUONGeometryConstituent(const TString& name, Int_t copyNo, 
                         Int_t npar, Double_t* param); 
    AliMUONGeometryConstituent(const TString& name, Int_t copyNo, 
                         const TGeoTranslation& translation,
			 Int_t npar, Double_t* param); 
    AliMUONGeometryConstituent(const TString& name, Int_t copyNo, 
                         const TGeoTranslation& translation, 
	  	         const TGeoRotation& rotation,
			 Int_t npar, Double_t* param);
     AliMUONGeometryConstituent(const TString& name, Int_t copyNo, 
                         const TGeoCombiTrans& transform,
			 Int_t npar, Double_t* param);
   AliMUONGeometryConstituent();
    virtual ~AliMUONGeometryConstituent();

    // get methods
    Int_t                  GetCopyNo() const;  
    Int_t                  GetNpar() const;
    Double_t*              GetParam() const;
    const TGeoCombiTrans*  GetTransformation() const;

  protected:
    AliMUONGeometryConstituent(const AliMUONGeometryConstituent& rhs);

    // operators  
    AliMUONGeometryConstituent& operator = (const AliMUONGeometryConstituent& rhs);

  private:
    Int_t            fCopyNo;        ///< copy number
    Int_t            fNpar;          ///< number of shape parameters
    
    /// shape parameters
    Double_t*        fParam;         //[fNpar] shape parameters

    TGeoCombiTrans*  fTransformation;///< \brief the constituent transformation
                                     ///  wrt to the envelope
 
  ClassDef(AliMUONGeometryConstituent,1) // MUON chamber geometry base class
};

// inline functions

inline Int_t AliMUONGeometryConstituent::GetCopyNo() const
{ return fCopyNo; }  

inline Int_t AliMUONGeometryConstituent::GetNpar() const
{ return fNpar; }

inline Double_t* AliMUONGeometryConstituent::GetParam() const
{ return fParam; }

inline const TGeoCombiTrans* AliMUONGeometryConstituent::GetTransformation() const 
{ return fTransformation; }

#endif //ALI_MUON_GEOMETRY_CONSTITUENT_H
