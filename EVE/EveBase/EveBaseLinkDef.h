// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#pragma link off all functions;
#pragma link off all globals;
#pragma link off all classes;

// AliEveConfigManager
#pragma link C++ class AliEveConfigManager+;

// AliEveEventManager
#pragma link C++ class AliEveInit+;
#pragma link C++ class AliEveEventManager+;
#pragma link C++ class AliEveDataSourceOffline+;
#pragma link C++ class AliEveEventManagerEditor+;
#pragma link C++ class AliEveEventManagerWindow+;
#pragma link C++ class AliEveDataSource+;
#ifdef ZMQ
#pragma link C++ class AliEveDataSourceHLTZMQ+;
#pragma link C++ class AliEveDataSourceOnline+;
#endif

// AliEveSaveViews
#pragma link C++ class AliEveSaveViews+;

// AliEveEventSelector
#pragma link C++ class AliEveEventSelector+;
#pragma link C++ class AliEveEventSelectorWindow+;

// AliEveMacro and AliEveMacroExecutor
#pragma link C++ class AliEveMacro+;
#pragma link C++ class AliEveMacroExecutor+;

// AliEvePreferencesWindow
#pragma link C++ class AliEvePreferencesWindow+;

// Various
#pragma link C++ class AliEveMagField+;
#pragma link C++ class AliEveMultiView+;

// AliEveTrack
#pragma link C++ class AliEveTrack+;
#pragma link C++ class AliEveTracklet+;

// AliEveTrackcounter
#pragma link C++ class AliEveTrackCounter+;
#pragma link C++ class AliEveTrackCounterEditor+;

// AliEveCascade
#pragma link C++ class AliEveCascade+;
#pragma link C++ class AliEveCascadeEditor+;
#pragma link C++ class AliEveCascadeList+;
#pragma link C++ class AliEveCascadeListEditor+;

// AliEveV0
#pragma link C++ class AliEveV0+;
#pragma link C++ class AliEveV0List+;
#pragma link C++ class AliEveV0Editor+;
#pragma link C++ class AliEveV0ListEditor+;

// AliEveKink
#pragma link C++ class AliEveKink+;
#pragma link C++ class AliEveKinkList+;
#pragma link C++ class AliEveKinkEditor+;
#pragma link C++ class AliEveKinkListEditor+;

// Common visualisation classes
#pragma link C++ class AliEveCutsWindow+;
#pragma link C++ class AliEveESDCascades+;
#pragma link C++ class AliEveESDKinks+;
#pragma link C++ class AliEveESDMuonTracks+;
#pragma link C++ class AliEveESDTracks+;
#pragma link C++ class AliEveESDSPDTracklets+;
#pragma link C++ class AliEveESDV0s+;
#pragma link C++ class AliEveMomentumHistograms+;
#pragma link C++ class AliEveMomentumVectors+;
#pragma link C++ class AliEvePrimaryVertex+;
#pragma link C++ class AliEveKineTools+;
#pragma link C++ class AliEveKineTracks+;
