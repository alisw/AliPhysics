void
polar()
{
   TCanvas *c1 = new TCanvas("c1","c1",500,500);
   TGraphPolar * grP1 = new TGraphPolar();
   grP1->SetTitle("TGraphPolar example");

   for (Int_t i = 0; i < 2; i++) {
     grP1->SetPoint(i, (i+1) * TMath::Pi() / 4, (i+1) * 0.05);
     grP1->SetPointError(i, 9*TMath::Pi()/180, 0.0);
     Double_t rr = grP1->GetY()[i];
     Double_t tt = grP1->GetX()[i];
     Double_t x  = rr * TMath::Cos(tt);
     Double_t y  = rr * TMath::Sin(tt);
     Printf("(x,y)=(%f,%f)", x, y);
   }
   Double_t r = 1;
   TH2* frame = new TH2F("frame", "Frame", 100, -r,r, 100, -r, r);
   frame->Draw();

   grP1->SetMarkerStyle(1);
   grP1->SetMarkerSize(1.);
   grP1->SetMarkerColor(4);
   grP1->SetLineColor(4);
   grP1->SetLineWidth(3);
   grP1->SetFillColor(kRed+1);
   grP1->SetFillStyle(3001);
   grP1->Draw("PNEF same");
   // grP1->Draw("APNEF");

   // Update, otherwise GetPolargram returns 0
   c1->Update();
   TGraphPolargram* gram = grP1->GetPolargram();
   gram->SetLineWidth(0);
   // gram->SetLineColor(kWhite);
   gram->SetNdivPolar(20);
   gram->SetNdivRadial(10);   
   gram->SetTextSize(0);
   gram->SetRadialLabelSize(0);
   gram->SetPolarLabelSize(0);
   gram->SetAxisAngle(9*TMath::Pi()/180);
   gram->SetTwoPi();
   gram->SetToRadian();
   c1->SetGridx();
   c1->SetGridy();
   c1->Update();
}
