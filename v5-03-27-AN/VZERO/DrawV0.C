void DrawV0()
{
   TGeoVolume *top = gGeoManager->GetMasterVolume();
   gGeoManager->SetNsegments(80);
   Int_t nd = top->GetNdaughters();
   for (Int_t i=0; i<nd; i++) top->GetNode(i)->GetVolume()->InvisibleAll();
   TGeoVolume *v0ri = gGeoManager->GetVolume("V0RI");  
   TGeoVolume *v0le = gGeoManager->GetVolume("V0LE");
   v0ri->SetVisibility(kTRUE);
   v0ri->VisibleDaughters(kTRUE);
   v0le->SetVisibility(kTRUE);
   v0le->VisibleDaughters(kTRUE);
   top->SetVisibility(kTRUE);
   top->Draw();
}
