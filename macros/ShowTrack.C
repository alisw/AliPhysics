//***This macro shows a track pointed by its index in the kinematics data base.
//***Use this macro only after starting a display.C macro !
void ShowTrack(int idx) {
   AliTPC *TPC=(AliTPC*)gAlice->GetDetector("TPC");
   TClonesArray *points=TPC->Points();
   int ntracks=points->GetEntriesFast();
   AliPoints *trk=0;
   for (int track=0;track<ntracks;track++) {
      AliPoints *pm = (AliPoints*)points->UncheckedAt(track);
      if (!pm) continue;
      if (idx == pm->GetIndex()) {
         pm->SetMarkerColor(2);
         pm->SetMarkerStyle(22);
         pm->Draw("same");
         trk=pm;
         break;
      }
   }
   if (!trk) return;

   TPad *pp=(TPad*)gROOT->FindObject("viewpad");
   pp->Update();
   pp->Modified();

   TClonesArray *particles=gAlice->Particles();
   GParticle *p = (GParticle*)particles->UncheckedAt(idx);
   cout<<"Paritcle ID "<<p->GetKF()<<endl;
   cout<<"Parent "<<p->GetParent()<<endl;
   cout<<"First child "<<p->GetFirstChild()<<endl;
   cout<<"Px,Py,Pz "<<p->GetPx()<<' '<<p->GetPy()<<' '<<p->GetPz()<<endl;

   cerr<<"OK ? "; char c[100]; cin.getline(c,100);
   if (tolower(c[0])!='y') return;

   trk->SetMarkerColor(5);
   trk->SetMarkerStyle(1);
   trk->Draw("same");
   pp->Update();
   pp->Modified();
}
