/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

///
/// \class AliMUONContourHandler
/// 
/// Class used to create contours of the muon tracker parts :
/// manu, bus patches, detection elements, chambers.
///
/// \author Laurent Aphecetche, Subatech
///

#include "AliMUONContourHandler.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONContour.h"
#include "AliMUONContourMaker.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONManuContourMaker.h"
#include "AliMUONPolygon.h"
#include "AliMUONSegment.h"
#include "AliMpArea.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpExMap.h"
#include <float.h>
#include "Riostream.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGeoMatrix.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TPolyLine.h"
#include "TSystem.h"

///\cond CLASSIMP
ClassImp(AliMUONContourHandler)
///\endcond

//_____________________________________________________________________________
AliMUONContourHandler::AliMUONContourHandler(Bool_t explodedView) 
: TObject(),
fTransformations(0x0),
fAllContourMap(0x0),
fAllContourArray(0x0)
{
  /// ctor. If explodedView=kTRUE, we generate views that will look 
  /// fine in 2D (i.e. won't show overlaps of DE as in reality)
  /// Use kFALSE if you want realistic geometry, though.
  ///
  /// IMPORTANT : as we many MUON classes, this one cannot work
  /// before you've loaded the mapping in memory (see e.g. AliMpCDB::LoadDDLStore)
  ///
  
  if ( explodedView ) 
  {
    fTransformations = GenerateTransformations(explodedView);
  }
  
  AliMUONManuContourMaker manuMaker(fTransformations);
  
  TObjArray* manus = manuMaker.GenerateManuContours(kTRUE);
  
  if (manus)
  {
    manus->SetOwner(kFALSE);
    GenerateAllContours(*manus);
  }
  delete manus;
}

//_____________________________________________________________________________
AliMUONContourHandler::~AliMUONContourHandler()
{
  /// Dtor
  delete fTransformations;
  delete fAllContourMap;
  delete fAllContourArray;
}

//______________________________________________________________________________
Bool_t 
AliMUONContourHandler::Adopt(AliMUONContour* contour)
{
  /// Adopt the given contour
  if ( GetContour(contour->GetName()) ) 
  {
    // contour already exists
    return kFALSE;
  }
  fAllContourMap->Add(new TObjString(contour->GetName()),contour);
  return kTRUE;
}

//______________________________________________________________________________
TObjArray*
AliMUONContourHandler::CreateContourList(const TObjArray& manuContours)
{    
  /// Create an array of maps of contour names
  ///
  /// Assyming that key is something like station#/chamber#/de#/buspatch#/manu#
  /// the idea here is to put one TMap for each level in mapArray :
  ///
  /// mapArray[0].key = station0
  /// mapArray[0].value = map of strings { station0/chamber0, station0/chamber1 }
  ///
  /// Then each entry in mapArray will be converted into a contour by
  /// merging its children (e.g. station0 contour will be made from the merging
  /// of station0/chamber0 and station0/chamber1 in the example above).
  ///
  
  AliCodeTimerAuto("",0);
  
  Int_t start(0);
  
  if ( !fTransformations )
  {
    start = 2; // skip chamber and station contours, as we seem to be in 3D
  }
    
  TIter next(&manuContours);
  AliMUONContour* contour;
  TObjArray* mapArray = new TObjArray;
  
  while ( ( contour = static_cast<AliMUONContour*>(next()) ) )
  {
    // Key is something like station#/chamber#/de#/buspatch#/manu#
    
    TString key(contour->GetName());
    TObjArray* s = key.Tokenize("/");
    for ( Int_t i = start; i < s->GetLast(); ++i ) 
    {
      TMap* m = static_cast<TMap*>(mapArray->At(i));
      if (!m)
      {
        m = new TMap;
        if ( i > mapArray->GetSize() ) mapArray->Expand(i);
        mapArray->AddAt(m,i);
      }
      TString parent;
      for ( Int_t k = 0; k <= i; ++k )
      {
        TObjString* str = static_cast<TObjString*>(s->At(k));
        parent += str->String();
        if ( k < i ) parent += "/";
      }
      TString child(parent);
      child += "/";
      child += static_cast<TObjString*>(s->At(i+1))->String();
      
      TObjArray* ma = static_cast<TObjArray*>(m->GetValue(parent.Data()));
      if (!ma)
      {
        ma = new TObjArray;
        m->Add(new TObjString(parent.Data()),ma);
      }
      TPair* p = static_cast<TPair*>(ma->FindObject(child.Data()));
      if ( !p ) 
      {
        ma->Add(new TObjString(child.Data()));
      }
    }
    delete s;
  }
  
  return mapArray;
}

//______________________________________________________________________________
void
AliMUONContourHandler::GenerateAllContours(const TObjArray& manuContours)
{
  /// From a map of manu contours, generate the compound contours (bp, de, etc...)
  /// by merging them.
  /// Note that manuContours should NOT be the owner of its contours,
  /// as they are adopted by the array returned by this method.
  
  AliCodeTimerAuto("",0);
  
  // Get the list of contours to create
  TObjArray* mapArray = CreateContourList(manuContours);
  
  fAllContourMap = new TMap(20000,1);
  fAllContourMap->SetOwnerKeyValue(kTRUE,kTRUE);
  
  fAllContourArray = new TObjArray;
  fAllContourArray->SetOwner(kFALSE); // the map is the real owner of the contours
  
  TIter nextContour(&manuContours);  
  AliMUONContour* contour(0x0);
  
  // start by adding the manu contours we begin with
  while ( ( contour = static_cast<AliMUONContour*>(nextContour()) ) )
  {
    fAllContourMap->Add(new TObjString(contour->GetName()),contour);
  }
  
  AliMUONContourMaker maker;
  
  for ( Int_t i = mapArray->GetLast(); i >= 1; --i ) 
    // end at 1 to avoid merging different cathodes together, which
    // would not work...
  {
    TMap* a = static_cast<TMap*>(mapArray->At(i));
    TIter next3(a);
    TObjString* str;
    while ( ( str = static_cast<TObjString*>(next3()) ) )
    {
      TObjArray* m = static_cast<TObjArray*>(a->GetValue(str->String().Data()));
      TIter next4(m);
      TObjString* k;
      TObjArray subcontours;
      subcontours.SetOwner(kFALSE);
      while ( ( k = static_cast<TObjString*>(next4()) ) )
      {
        contour = static_cast<AliMUONContour*>(fAllContourMap->GetValue(k->String().Data()));
        if ( contour ) 
        {
          subcontours.Add(contour);
        }
        else
        {
          AliError(Form("Did not find contour %s",k->String().Data()));
          continue;
        }
      }
      
      contour = maker.MergeContour(subcontours,str->String().Data());
      
      bool error(kFALSE);
      
      if (!contour)
      {
        error=kTRUE;
        AliError(Form("ERROR : could not merge into %s",str->String().Data()));
      }
      else
      {
        if ( contour->Area().IsValid() == kFALSE ) 
        {
          error=kTRUE;
          AliError(Form("ERROR : area of contour %s is invalid",str->String().Data()));
        }
      }
      
      if (!error)
      {
        fAllContourMap->Add(new TObjString(str->String().Data()),contour);
        fAllContourArray->Add(contour);
      }
    }
  }
}

//_____________________________________________________________________________
AliMpExMap* 
AliMUONContourHandler::GenerateTransformations(Bool_t exploded)
{
  /// Generate geometric transformations to be used to compute the contours
  /// If exploded=kFALSE then we generate real transformations, otherwise
  /// we generate tweaked ones that look fine on screen.
  
  AliCodeTimerAuto("",0);
  
  AliMUONGeometryTransformer transformer;
  Bool_t ok = transformer.LoadGeometryData("transform.dat");
  //  transformer.LoadGeometryData("geometry.root"); //FIXME: add a protection if geometry.root file does not exist
  if (!ok)
  {
    cout << "ERROR : cannot get geometry !" << endl;
    return 0x0;
  }
  AliMpExMap* transformations = new AliMpExMap;
  AliMpDEIterator deIt;
  deIt.First();
  while ( !deIt.IsDone() )
  {
    Int_t detElemId = deIt.CurrentDEId();
    const AliMUONGeometryDetElement* de = transformer.GetDetElement(detElemId);
    
    TGeoHMatrix* matrix = static_cast<TGeoHMatrix*>(de->GetGlobalTransformation()->Clone());

    if (exploded)
    {
      Double_t* translation = matrix->GetTranslation();
      Double_t xscale = 1.0;
      Double_t yscale = 1.5;
      Double_t shift = 5.0; // cm
          
      if ( AliMpDEManager::GetStationType(detElemId) == AliMp::kStation345 ) 
      {
        translation[0] *= xscale;
        translation[1] *= yscale; 
      }
      else
      {
        Double_t xshift[] = { shift, -shift, -shift, shift };
        Double_t yshift[] = { shift, shift, -shift, -shift };
        Int_t ishift = detElemId % 100;
        
        translation[0] += xshift[ishift];
        translation[1] += yshift[ishift];
      }
      matrix->SetTranslation(translation);
    }
    transformations->Add(detElemId,matrix);
    deIt.Next();
  }
  return transformations;
}

//_____________________________________________________________________________
AliMUONContour* 
AliMUONContourHandler::GetContour(const char* contourname) const
{
  /// Get a given contour
  return static_cast<AliMUONContour*>(fAllContourMap->GetValue(contourname));
}

//_____________________________________________________________________________
void
AliMUONContourHandler::Print(Option_t* opt) const
{
  /// printout
  
  if ( ! fAllContourMap )  return;
  
  cout << Form("Contour map : collisions = %5.3f size = %d capacity = %d", 
               fAllContourMap->AverageCollisions(),
               fAllContourMap->GetSize(),
               fAllContourMap->Capacity()) << endl;

  TString sopt(opt);
  sopt.ToUpper();
  
  if ( sopt.Contains("ALL") || sopt.Contains("FULL") )
  {
    fAllContourMap->Print();
  }
}
