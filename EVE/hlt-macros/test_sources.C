AliEveHOMERManager* homerM = 0;

void test_sources()
{
  homerM = new AliEveHOMERManager("./sampleConfig.xml");

  gEve->AddToListTree(homerM, kTRUE);

  homerM->CreateHOMERSourcesList();

  AliEveHOMERSourceList* srcL = new AliEveHOMERSourceList("Sources");
  srcL->SetManager(homerM);
  homerM->AddElement(srcL);

  srcL->CreateByType();

  /*
  TList* srcList = homerM->GetSourceList();

  AliEveHOMERSourceMap* smd = AliEveHOMERSourceMap::Create(AliEveHOMERSourceMap::kSG_ByDet);
  smd->FillMap(srcList, 1);
  printf(" **** ByDet  XXX ****\n");
  smd->PrintXXX();

  AliEveHOMERSourceMap* smt = AliEveHOMERSourceMap::Create(AliEveHOMERSourceMap::kSG_ByType);
  smt->FillMap(srcList, 1);
  printf(" **** ByType XXX ****\n");
  smt->PrintXXX();

  */
}
