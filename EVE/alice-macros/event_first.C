// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void event_first()
{
  if (AliEveEventManager::GetMaster() == 0) {
    printf("AliEveEventManager is not initialized!\n");
    return;
  }
  AliEveEventManager::GetMaster()->GotoEvent(0);
}
