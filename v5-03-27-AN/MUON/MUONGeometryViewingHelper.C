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
//
/// \ingroup macros
/// \file MUONGeometryViewingHelper.C
/// \brief Macro providing methods helping in viewing geometry
///
/// To be run from aliroot:
/// <pre>
/// root[0] .x MUONGeometryViewingHelper.C
/// root[0] buildGeometry();
/// </pre>
///
/// \author: I. Hrivnacova, IPN Orsay

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TObjArray.h>
#include <TBrowser.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>

#endif


void visibilityOff() 
{
/// Set all volumes invisible

  TObjArray* volumes = gGeoManager->GetListOfVolumes();
  for (Int_t i=0; i<volumes->GetEntriesFast(); i++) {
    if ( !((TGeoVolume*)volumes->At(i))->IsAssembly() )
      ((TGeoVolume*)volumes->At(i))->SetVisibility(kFALSE);
  }    
    
}  


void setVisibility(const TString& volumeName, Bool_t visibility= kTRUE) 
{
/// Set visibility to the volume specified by name

  TGeoVolume* volume = gGeoManager->FindVolumeFast(volumeName.Data());

  if ( ! volume ) {
    cerr << "Volume " <<  volumeName.Data() << " not found." << endl;
    return;
  }  
  
  volume->SetVisibility(visibility);
}  

void setDaughtersVisibility(const TString& volumeName, Bool_t visibility= kTRUE)  
{
/// Set visibility to daughter of the volume specified by name.
/// If the daughter volume is an assembly the visibility setting
/// is propagated to its real volumes daughters.

   TGeoVolume* volume = gGeoManager->FindVolumeFast(volumeName.Data());

   if ( ! volume ) { 
     cerr << "Volume " <<  volumeName.Data() << " not found." << endl;
     return;
   }  
     
   //for ( Int_t i=0; i<10; i++ ) {
   Int_t colourNo = 1;
   for ( Int_t i=0; i<volume->GetNdaughters(); i++ ) {
   
     TGeoVolume* dvolume = volume->GetNode(i)->GetVolume();
     if ( dvolume->IsAssembly() ) {
       // Do no set visibility to assembly but to its daughters
       setDaughtersVisibility(dvolume->GetName(), visibility);
     }
     else {  
       cout << "Setting visibility to " <<  dvolume->GetName() << endl;
       dvolume->SetVisibility(visibility);

       // change colors
       ++colourNo;
       if ( colourNo > 9 ) colourNo = 1;
       dvolume->SetLineColor(colourNo); 
     }     
   }
}

void buildGeometry(Bool_t allVisible = kFALSE ) 
{  
/// Load geometry from the file, make all volumes invisible

  TGeoManager::Import("geometry.root");
  
  new TBrowser();
  
  if ( ! allVisible ) visibilityOff();

  gGeoManager->SetVisLevel(10);
  gGeoManager->GetTopVolume()->SetVisContainers(kTRUE);
  gGeoManager->GetTopVolume()->Draw("ogl");

  cout << endl;
  cout << "You can now add volumes in the scene: " << endl;
  cout << "    setVisibility(\"MyVolume\") "  <<  endl;
  cout << "    setDaughtersVisibility(\"MyVolume\") " <<  endl;;
  cout << "or remove them from the scene: " << endl;
  cout << "    setVisibility(\"MyVolume\", kFALSE);" << endl;
  cout << "    setDaughtersVisibility(\"MyVolume\", kFALSE);" << endl;  
  cout << endl;
}  
