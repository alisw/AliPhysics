/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpGraphContext.h,v 1.11 2006/05/24 13:58:13 ivana Exp $

/// \ingroup mpgraphics
/// \class AliMpGraphContext
/// \brief Class describing the correspondance between a given area
/// in pad, and a zone of real (cm) position
///
/// \author David GUEZ, IPN Orsay

#ifndef ALI_MP_GRAPH_CONTEXT_H
#define ALI_MP_GRAPH_CONTEXT_H

#include <TObject.h>

#include "AliMpExMap.h"

#include <TVector2.h>

class MPainter;

class AliMpGraphContext : public TObject
{
 public:
  void Push() const;
  void Pop();
  static AliMpGraphContext *Instance();

  //
  // set methods
  //
           /// Set position of the pad area where to draw
  void SetPadPosition(const TVector2 &position){fPadPosition=position;}
           /// Set dimensions of the pad area where to draw
  void SetPadDimensions(const TVector2 &dimensions){fPadDimensions=dimensions;}
           /// Set position of the real area where to draw
  void SetRealPosition(const TVector2 &position){fRealPosition=position;}
           /// Set dimensions of the real area where to draw
  void SetRealDimensions(const TVector2 &dimensions){fRealDimensions=dimensions;}
           /// Set color to use
  void SetColor(Int_t color){fColor=color;}

  //
  // get methods
  //

           /// Return position of the pad area where to draw
  TVector2 GetPadPosition() const {return fPadPosition;}
           /// Return dimensions of the pad area where to draw
  TVector2 GetPadDimensions() const {return fPadDimensions;}
           /// Return position of the real area where to draw
  TVector2 GetRealPosition() const{return fRealPosition;}
           /// Return dimensions of the real area where to draw
  TVector2 GetRealDimensions() const{return fRealDimensions;}
           /// Return color to use
  Int_t GetColor() const {return fColor;}

  //methods
  TVector2 RealToPad(const TVector2 &position) const;
  void RealToPad(const TVector2 &position
		 ,const TVector2 &dimensions,
		 TVector2 &padPosition,
		 TVector2 &padDimensions) const;
  void SetPadPosForReal(const TVector2 &position,const TVector2 &dimensions);

 protected:
  AliMpGraphContext(const AliMpGraphContext& right);
  AliMpGraphContext&  operator = (const AliMpGraphContext& right);

 private:
  //private constructor (not instanciable from outside)
  AliMpGraphContext();

  ///< static data members
  static AliMpGraphContext* fgInstance; ///< the global instance
  static TObjArray          fgStack;    ///< the object stack
  static Int_t              fgStackSize;///< the object stack size

  //data members
  Int_t    fColor;          ///< color to use
  TVector2 fPadPosition;    ///< Position of the pad area where to draw
  TVector2 fPadDimensions;  ///< Dimensions of the pad area where to draw

  TVector2 fRealPosition;   ///< Position of the real area to draw
  TVector2 fRealDimensions; ///< Dimensions of the real area to draw

  ClassDef(AliMpGraphContext,1) // Correspondance pad area/real world
};

#endif //ALI_MP_GRAPH_CONTEXT_H


