// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveConfigManager_H
#define AliEveConfigManager_H

#include "TObject.h"

class TGPopupMenu;

//______________________________________________________________________________
// Short description of AliEveConfigManager
//

class AliEveConfigManager : public TObject
{
public:
  static AliEveConfigManager* InitializeMaster();
  static AliEveConfigManager* GetMaster();

  virtual ~AliEveConfigManager() {}

  void AliEvePopupHandler(Int_t id);

protected:
  static AliEveConfigManager* fgMaster;  // Main instance.

  TGPopupMenu      *fAliEvePopup; // AliEve menu.

private:
  AliEveConfigManager();

  AliEveConfigManager(const AliEveConfigManager&);            // Not implemented
  AliEveConfigManager& operator=(const AliEveConfigManager&); // Not implemented

  ClassDef(AliEveConfigManager, 0); // Short description.
};

#endif
