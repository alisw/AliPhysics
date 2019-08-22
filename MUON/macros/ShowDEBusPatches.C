// macro to show all the buspatches of one detection element
//

#include "AliCDBManager.h"
#include "AliMpBusPatch.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include <algorithm>
#include <iostream>
#include <vector>

void ShowDEBusPatches(std::vector<int> deids,
                      const char *ocdb = "local://$ALIROOT_OCDB_ROOT/OCDB",
                      int runNumber = 297624) {
  AliCDBManager::Instance()->SetDefaultStorage(ocdb);
  AliCDBManager::Instance()->SetRun(runNumber);
  AliMpCDB::LoadAll(true);

  for (auto deid : deids) {
    AliMpDetElement *de = AliMpDDLStore::Instance()->GetDetElement(deid);
    std::cout << "DE " << deid << " " << de->NofChannels() << " channels\n";
    for (int i = 0; i < de->GetNofBusPatches(); i++) {
      int bpid = de->GetBusPatchId(i);
      AliMpBusPatch *bp = AliMpDDLStore::Instance()->GetBusPatch(bpid);
      bp->Print("FULL MANU");
    }
  }
}
