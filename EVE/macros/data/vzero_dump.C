#if !defined(__CINT__) || defined(__MAKECINT__)
#include <AliESDEvent.h>
#include <AliESDVZERO.h>
#include <AliEveEventManager.h>
#endif

void vzero_dump()
{
  AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();

  AliESDVZERO *vz = esd->GetVZEROData();
  if (vz)
  {
    printf("====================================================================\n");
    printf("================   VZERO DUMP    ===================================\n");
    printf("====================================================================\n");
    printf("|  id   multip     time    BB BG ||  id   multip     time    BB BG |\n");
    printf("+------------------------------++----------------------------------+\n");
    for (Int_t i=0; i<32; ++i)
    {
      Int_t j = i + 32;
      printf("|  %2d %8.1f %8.1f   %2d %2d  ||  %2d %8.1f %8.1f   %2d %2d  |\n",
             i, vz->GetMultiplicity(i), vz->GetTime(i), vz->GetBBFlag(i), vz->GetBGFlag(i),
             j, vz->GetMultiplicity(j), vz->GetTime(j), vz->GetBBFlag(j), vz->GetBGFlag(j));
    }
    printf("====================================================================\n");
  }

}
