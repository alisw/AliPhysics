void ecut() {

//Begin_Html
/*
<img src="picts/ecut.gif">
*/
//End_Html
   AliDisplay *edisplay;
   TSlider *etacut = edisplay->EtaSlider();
   for (Int_t i=0;i<50;i++) {
      Float_t ymin = 0.02*Float_t(i);
      Float_t ymax = ymin +0.02;
      etacut->SetRange(ymin,ymax);
      //etacut->SetMinimum(ymin);
      //etacut->SetMaximum(ymax);
      edisplay->Draw();
      Canvas->Update();
   }
}
