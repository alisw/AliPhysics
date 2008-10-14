#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AlienGUILinkDef.h 17141 2007-02-28 12:57:06Z hristov $ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliAnalysisGUI+;
#ifdef WITHXML
#pragma link C++ class AliAlienBrowser+;
#pragma link C++ class AliLoginFrame+;
#pragma link C++ class AliFileListFrame+;
#pragma link C++ class AliPackageFrame+;
#pragma link C++ class AliSelectorFrame+;
#pragma link C++ class AliTagFrame+;
#pragma link C++ class AliTagAnalysisFrame+;
#endif

#endif
