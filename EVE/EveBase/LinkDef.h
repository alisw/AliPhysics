// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#pragma link off all functions;
#pragma link off all globals;
#pragma link off all classes;

// AliEveEventManager
#pragma link C++ class AliEveEventManager+;
#pragma link C++ class AliEveEventManagerEditor+;
#pragma link C++ class AliEveEventManagerWindow+;

#pragma link C++ global gAliEveEvent;

#pragma link C++ class AliEveKineTools+;

#pragma link C++ class AliEveVSDCreator+;

// AliEveTrackcounter
#pragma link C++ class AliEveTrackCounter+;
#pragma link C++ class AliEveTrackCounterEditor+;

// AliEveTrackFitter
#pragma link C++ class AliEveTrackFitter+;
#pragma link C++ class AliEveTrackFitterEditor+;

#pragma link C++ class AliEveCosmicRayFitter+;
#pragma link C++ class AliEveCosmicRayFitterEditor+;

// AliEveJetPlane
#pragma link C++ class AliEveJetPlane+;
#pragma link C++ class AliEveJetPlaneGL+;
#pragma link C++ class AliEveJetPlaneEditor+;
#pragma link C++ class AliEveJetPlaneEditor::StaticDataWindow+;


// Removed. Messy code, tons of violations and incompatible with TEve
// classes. Author Ludovic Gaudichet left ALICE.
// Should be thoroughly revised.
//
// AliEveCascade
// #pragma link C++ class AliEveCascade+;
// #pragma link C++ class CascadeList+;
// #pragma link C++ class CascadeListEditor+;
//
// AliEveV0
#pragma link C++ class AliEveV0+;
#pragma link C++ class AliEveV0List+;
#pragma link C++ class AliEveV0Editor+;
#pragma link C++ class AliEveV0ListEditor+;
