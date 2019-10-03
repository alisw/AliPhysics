{

TFile f("output/pp-trunk010709-pdc1-00XXX.root");
//TFile f("results_V0QA_v416rev08_090709.root");

THnSparseF* sparseV0;

sparseV0=(THnSparseF*)(outputTRD->FindObject("sparseV0"));

sparseV0->GetAxis(5)->SetRangeUser(0.,3.9);
sparseV0->GetAxis(9)->SetRangeUser(0.,3.9);


TH1* ptgeantPPhys=sparseV0->Projection(5);
TH1* etageantPPhys=sparseV0->Projection(6);
TH1* phigeantPPhys=sparseV0->Projection(7);

TH1* ptResSinglePPhys=sparseV0->Projection(13);
TH1* etasinglePPhys=sparseV0->Projection(14);
TH1* phisinglePPhys=sparseV0->Projection(15);


TH1* ptgeantNPhys=sparseV0->Projection(9);
TH1* etageantNPhys=sparseV0->Projection(10);
TH1* phigeantNPhys=sparseV0->Projection(11);

TH1* ptResSingleNPhys=sparseV0->Projection(21);
TH1* etasingleNPhys=sparseV0->Projection(22);
TH1* phisingleNPhys=sparseV0->Projection(23);

TH1* efiEtaSingPPhys=(TH1*) etageantPPhys->Clone();
efiEtaSingPPhys->Sumw2();
efiEtaSingPPhys->Divide(etasinglePPhys,etageantPPhys);
efiEtaSingPPhys->SetMinimum(0.);
efiEtaSingPPhys->SetMaximum(1.2);
efiEtaSingPPhys->SetAxisRange(-1.5,1.5);
efiEtaSingPPhys->SetMarkerStyle(20);
efiEtaSingPPhys->SetMarkerColor(2);



TH1* efiEtaSingNPhys=(TH1*) etageantNPhys->Clone();
efiEtaSingNPhys->Sumw2();
efiEtaSingNPhys->Divide(etasingleNPhys,etageantNPhys);
efiEtaSingNPhys->SetMinimum(0.);
efiEtaSingNPhys->SetMaximum(1.2);
efiEtaSingNPhys->SetAxisRange(-1.5,1.5);
efiEtaSingNPhys->SetMarkerStyle(20);
efiEtaSingNPhys->SetMarkerColor(2);



TH1* efiPhiSingPPhys=(TH1*) phigeantPPhys->Clone();
efiPhiSingPPhys->Sumw2();
efiPhiSingPPhys->Divide(phisinglePPhys,phigeantPPhys);
efiPhiSingPPhys->SetMinimum(0.);
efiPhiSingPPhys->SetMaximum(1.2);
efiPhiSingPPhys->SetAxisRange(-TMath::Pi(),TMath::Pi());
efiPhiSingPPhys->SetMarkerStyle(20);
efiPhiSingPPhys->SetMarkerColor(2);


TH1* efiPhiSingNPhys=(TH1*) phigeantNPhys->Clone();
efiPhiSingNPhys->Sumw2();
efiPhiSingNPhys->Divide(phisingleNPhys,phigeantNPhys);
efiPhiSingNPhys->SetMinimum(0.);
efiPhiSingNPhys->SetMaximum(1.2);
efiPhiSingNPhys->SetAxisRange(-TMath::Pi(),TMath::Pi());
efiPhiSingNPhys->SetMarkerStyle(20);
efiPhiSingNPhys->SetMarkerColor(2);


// Reconstructable tracks

// track Lenght
sparseV0->GetAxis(8)->SetRangeUser(50,200);
sparseV0->GetAxis(12)->SetRangeUser(50,200);
// momentum
sparseV0->GetAxis(5)->SetRangeUser(0.,3.9);
sparseV0->GetAxis(9)->SetRangeUser(0.,3.9);


TH1* rgeantTrk=sparseV0->Projection(3);
TH1* ptgeantPTrk=sparseV0->Projection(5);
TH1* etageantPTrk=sparseV0->Projection(6);
TH1* phigeantPTrk=sparseV0->Projection(7);


TH1* ptgeantNTrk=sparseV0->Projection(9);
TH1* etageantNTrk=sparseV0->Projection(10);
TH1* phigeantNTrk=sparseV0->Projection(11);

// Status of the single track reconstruction. "1" means reconstructed,"-1" not reconstructed

sparseV0->GetAxis(20)->SetRangeUser(0,2.);
sparseV0->GetAxis(28)->SetRangeUser(0,2.);

TH1* rsingleTrk=sparseV0->Projection(3);
TH1* ptsinglePTrk=sparseV0->Projection(5);
TH1* etasinglePTrk=sparseV0->Projection(6);
TH1* phisinglePTrk=sparseV0->Projection(7);

TH1* ptsingleNTrk=sparseV0->Projection(9);
TH1* etasingleNTrk=sparseV0->Projection(10);
TH1* phisingleNTrk=sparseV0->Projection(11);

//Status of the V0 track reconstruction.  "1" means reconstructed,"-1" not reconstructed

sparseV0->GetAxis(32)->SetRangeUser(0,2.);
sparseV0->GetAxis(36)->SetRangeUser(0,2.);

TH1* rV0Trk=sparseV0->Projection(3);
TH1* ptV0PTrk=sparseV0->Projection(5);
TH1* etaV0PTrk=sparseV0->Projection(6);
TH1* phiV0PTrk=sparseV0->Projection(7);



TH1* ptV0NTrk=sparseV0->Projection(9);
TH1* etaV0NTrk=sparseV0->Projection(10);
TH1* phiV0NTrk=sparseV0->Projection(11);


TH1* efiSingPPhys=(TH1*) ptgeantPPhys->Clone();
efiSingPPhys->Sumw2();
efiSingPPhys->Divide(ptgeantPTrk,ptgeantPPhys);
efiSingPPhys->SetAxisRange(0.,5);
efiSingPPhys->SetMinimum(0.);
efiSingPPhys->SetMaximum(1.2);
efiSingPPhys->SetMarkerStyle(20);
efiSingPPhys->SetMarkerColor(2);


TH1* efiSingNPhys=(TH1*) ptgeantNPhys->Clone();
efiSingNPhys->Sumw2();
efiSingNPhys->Divide(ptgeantNTrk,ptgeantNPhys);
efiSingNPhys->SetAxisRange(0.,5);
efiSingNPhys->SetMinimum(0.);
efiSingNPhys->SetMaximum(1.2);
efiSingNPhys->SetMarkerStyle(20);
efiSingNPhys->SetMarkerColor(2);




TH1* efiSingPTrk=(TH1*) ptgeantPTrk->Clone();
efiSingPTrk->Sumw2();
efiSingPTrk->Divide(ptsinglePTrk,ptgeantPTrk);
efiSingPTrk->SetAxisRange(0.,5);
efiSingPTrk->SetMinimum(0.);
efiSingPTrk->SetMaximum(1.2);
efiSingPTrk->SetMarkerStyle(21);
efiSingPTrk->SetMarkerColor(4);


TH1* efiSingNTrk=(TH1*) ptgeantNTrk->Clone();
efiSingNTrk->Sumw2();
efiSingNTrk->Divide(ptsingleNTrk,ptgeantNTrk);
efiSingNTrk->SetAxisRange(0.,5);
efiSingNTrk->SetMinimum(0.);
efiSingNTrk->SetMaximum(1.2);
efiSingNTrk->SetMarkerStyle(21);
efiSingNTrk->SetMarkerColor(4);


TH1* efiV0PTrk=(TH1*) ptsinglePTrk->Clone();
efiV0PTrk->Sumw2();
efiV0PTrk->Divide(ptV0PTrk,ptsinglePTrk);
efiV0PTrk->SetAxisRange(0.,5);
efiV0PTrk->SetMinimum(0.);
efiV0PTrk->SetMaximum(1.2);
efiV0PTrk->SetMarkerStyle(21);
efiV0PTrk->SetMarkerColor(3);



TH1* efiV0NTrk=(TH1*) ptsingleNTrk->Clone();
efiV0NTrk->Sumw2();
efiV0NTrk->Divide(ptV0NTrk,ptsingleNTrk);
efiV0NTrk->SetAxisRange(0.,5);
efiV0NTrk->SetMinimum(0.);
efiV0NTrk->SetMaximum(1.2);
efiV0NTrk->SetMarkerStyle(21);
efiV0NTrk->SetMarkerColor(3);



TH1* efiEtaSingPTrk=(TH1*) etageantPTrk->Clone();
efiEtaSingPTrk->Sumw2();
efiEtaSingPTrk->Divide(etasinglePTrk,etageantPTrk);
efiEtaSingPTrk->SetMinimum(0.);
efiEtaSingPTrk->SetMaximum(1.2);
efiEtaSingPTrk->SetAxisRange(-1.5,1.5);
efiEtaSingPTrk->SetMarkerStyle(21);
efiEtaSingPTrk->SetMarkerColor(4);


TH1* efiEtaSingNTrk=(TH1*) etageantNTrk->Clone();
efiEtaSingNTrk->Sumw2();
efiEtaSingNTrk->Divide(etasingleNTrk,etageantNTrk);
efiEtaSingNTrk->SetMinimum(0.);
efiEtaSingNTrk->SetMaximum(1.2);
efiEtaSingNTrk->SetAxisRange(-1.5,1.5);
efiEtaSingNTrk->SetMarkerStyle(21);
efiEtaSingNTrk->SetMarkerColor(4);


TH1* efiEtaV0PTrk=(TH1*) etasinglePTrk->Clone();
efiEtaV0PTrk->Sumw2();
efiEtaV0PTrk->Divide(etaV0PTrk,etasinglePTrk);
efiEtaV0PTrk->SetMinimum(0.);
efiEtaV0PTrk->SetMaximum(1.2);
efiEtaV0PTrk->SetMarkerStyle(21);
efiEtaV0PTrk->SetMarkerColor(3);


TH1* efiEtaV0NTrk=(TH1*) etasingleNTrk->Clone();
efiEtaV0NTrk->Sumw2();
efiEtaV0NTrk->Divide(etaV0NTrk,etasingleNTrk);
efiEtaV0NTrk->SetMinimum(0.);
efiEtaV0NTrk->SetMaximum(1.2);
efiEtaV0NTrk->SetMarkerStyle(21);
efiEtaV0NTrk->SetMarkerColor(3);



TH1* efiPhiSingPTrk=(TH1*) phigeantPTrk->Clone();
efiPhiSingPTrk->Sumw2();
efiPhiSingPTrk->Divide(phisinglePTrk,phigeantPTrk);
efiPhiSingPTrk->SetMinimum(0.);
efiPhiSingPTrk->SetMaximum(1.2);
efiPhiSingPTrk->SetAxisRange(-TMath::Pi(),TMath::Pi());
efiPhiSingPTrk->SetMarkerStyle(20);
efiPhiSingPTrk->SetMarkerColor(4);


TH1* efiPhiSingNTrk=(TH1*) phigeantNTrk->Clone();
efiPhiSingNTrk->Sumw2();
efiPhiSingNTrk->Divide(phisingleNTrk,phigeantNTrk);
efiPhiSingNTrk->SetMinimum(0.);
efiPhiSingNTrk->SetMaximum(1.2);
efiPhiSingNTrk->SetAxisRange(-TMath::Pi(),TMath::Pi());
efiPhiSingNTrk->SetMarkerStyle(20);
efiPhiSingNTrk->SetMarkerColor(4);

TH1* efiPhiV0PTrk=(TH1*) phisinglePTrk->Clone();
efiPhiV0PTrk->Sumw2();
efiPhiV0PTrk->Divide(phiV0PTrk,phisinglePTrk);
efiPhiV0PTrk->SetMinimum(0.);
efiPhiV0PTrk->SetMaximum(1.2);
efiPhiV0PTrk->SetMarkerStyle(21);
efiPhiV0PTrk->SetMarkerColor(3);


TH1* efiPhiV0NTrk=(TH1*) phisingleNTrk->Clone();
efiPhiV0NTrk->Sumw2();
efiPhiV0NTrk->Divide(phiV0NTrk,phisingleNTrk);
efiPhiV0NTrk->SetMinimum(0.);
efiPhiV0NTrk->SetMaximum(1.2);
efiPhiV0NTrk->SetMarkerStyle(21);
efiPhiV0NTrk->SetMarkerColor(3);



sparseV0->GetAxis(32)->SetRangeUser(-2,0.);
sparseV0->GetAxis(36)->SetRangeUser(-2,0.);

TH1* sigbXYZEPlusRecSingle=sparseV0->Projection(17);
TH1* sigbXYZEMinusRecSingle=sparseV0->Projection(25);

TH1* ptsinglePTrkNo=sparseV0->Projection(5);
TH1* etasinglePTrkNo=sparseV0->Projection(6);
TH1* phisinglePTrkNo=sparseV0->Projection(7);

TH1* ptsingleNTrkNo=sparseV0->Projection(9);
TH1* etasingleNTrkNo=sparseV0->Projection(10);
TH1* phisingleNTrkNo=sparseV0->Projection(11);





TH1* efiV0PTrkNo=(TH1*) ptsinglePTrkNo->Clone();
efiV0PTrkNo->Sumw2();
efiV0PTrkNo->Divide(ptsinglePTrkNo,ptsinglePTrk);
efiV0PTrkNo->SetAxisRange(0.,5);
efiV0PTrkNo->SetMinimum(0.);
efiV0PTrkNo->SetMaximum(1.2);
efiV0PTrkNo->SetMarkerStyle(21);
efiV0PTrkNo->SetMarkerColor(6);



TH1* efiV0NTrkNo=(TH1*) ptsingleNTrkNo->Clone();
efiV0NTrkNo->Sumw2();
efiV0NTrkNo->Divide(ptsingleNTrkNo,ptsingleNTrk);
efiV0NTrkNo->SetAxisRange(0.,5);
efiV0NTrkNo->SetMinimum(0.);
efiV0NTrkNo->SetMaximum(1.2);
efiV0NTrkNo->SetMarkerStyle(21);
efiV0NTrkNo->SetMarkerColor(6);


TH1* efiEtaV0PTrkNo=(TH1*) etasinglePTrkNo->Clone();
efiEtaV0PTrkNo->Sumw2();
efiEtaV0PTrkNo->Divide(etasinglePTrkNo,etasinglePTrk);
efiEtaV0PTrkNo->SetAxisRange(-1.5,1.5);
efiEtaV0PTrkNo->SetMinimum(0.);
efiEtaV0PTrkNo->SetMaximum(1.2);
efiEtaV0PTrkNo->SetMarkerStyle(21);
efiEtaV0PTrkNo->SetMarkerColor(6);



TH1* efiEtaV0NTrkNo=(TH1*) etasingleNTrkNo->Clone();
efiEtaV0NTrkNo->Sumw2();
efiEtaV0NTrkNo->Divide(etasingleNTrkNo,etasingleNTrk);
efiEtaV0NTrkNo->SetAxisRange(-1.5,1.5);
efiEtaV0NTrkNo->SetMinimum(0.);
efiEtaV0NTrkNo->SetMaximum(1.2);
efiEtaV0NTrkNo->SetMarkerStyle(21);
efiEtaV0NTrkNo->SetMarkerColor(6);


TH1* efiPhiV0PTrkNo=(TH1*) phisinglePTrkNo->Clone();
efiPhiV0PTrkNo->Sumw2();
efiPhiV0PTrkNo->Divide(phisinglePTrkNo,phisinglePTrk);
efiPhiV0PTrkNo->SetAxisRange(-3.5,3.5);
efiPhiV0PTrkNo->SetMinimum(0.);
efiPhiV0PTrkNo->SetMaximum(1.2);
efiPhiV0PTrkNo->SetMarkerStyle(21);
efiPhiV0PTrkNo->SetMarkerColor(6);



TH1* efiPhiV0NTrkNo=(TH1*) phisingleNTrkNo->Clone();
efiPhiV0NTrkNo->Sumw2();
efiPhiV0NTrkNo->Divide(phisingleNTrkNo,phisingleNTrk);
efiPhiV0NTrkNo->SetAxisRange(-3.5,3.5);
efiPhiV0NTrkNo->SetMinimum(0.);
efiPhiV0NTrkNo->SetMaximum(1.2);
efiPhiV0NTrkNo->SetMarkerStyle(21);
efiPhiV0NTrkNo->SetMarkerColor(6);


TCanvas* c1 = new TCanvas("c1","MC Info",200,10,700,900);
c1->Divide(2,4);
c1->SetFillColor(0);
c1->GetFrame()->SetFillColor(0);
c1->SetBorderMode(0);
TPad *c1_1 = (TPad*)(c1->GetPrimitive("c1_1"));
TPad *c1_2 = (TPad*)(c1->GetPrimitive("c1_2"));
TPad *c1_3 = (TPad*)(c1->GetPrimitive("c1_3"));
TPad *c1_4 = (TPad*)(c1->GetPrimitive("c1_4"));
TPad *c1_5 = (TPad*)(c1->GetPrimitive("c1_5"));
TPad *c1_6 = (TPad*)(c1->GetPrimitive("c1_6"));


 c1_1->SetGridx(1);
 c1_3->SetGridx(1);
 c1_2->SetGridx(1);
 c1_4->SetGridx(1);
 c1_5->SetGridx(1);
 c1_6->SetGridx(1);

 c1_1->SetGridy(1);
 c1_3->SetGridy(1);
 c1_2->SetGridy(1);
 c1_4->SetGridy(1);
 c1_5->SetGridy(1);
 c1_6->SetGridy(1);


c1->cd(1);
efiSingPTrk->Draw();
efiSingPTrk->SetXTitle("1/sqrt(pt) pos (sqrt(GeV))");
efiSingPPhys->Draw("same");
efiV0PTrk->Draw("same");
efiV0PTrkNo->Draw("same");

TLegend *leg1 = new TLegend(0.6,0.75,0.9,0.9);
leg1->AddEntry(efiSingPTrk,"Single tracking eff");
leg1->AddEntry(efiV0PTrk,"V0 tracking eff");
leg1->AddEntry(efiV0PTrkNo,"V0 tracking ineff");
leg1->Draw();
c1->cd(2);
efiSingNTrk->Draw();
efiSingNTrk->SetXTitle("1/sqrt(pt) neg (sqrt(GeV))");
efiSingNPhys->Draw("same");
efiV0NTrk->Draw("same");
efiV0NTrkNo->Draw("same");

c1->cd(3);
efiEtaSingPTrk->Draw();
//efiEtaSingPPhys->Draw("same");
efiEtaV0PTrk->Draw("same");
efiEtaV0PTrkNo->Draw("same");

c1->cd(4);
efiEtaSingNTrk->Draw();
//efiEtaSingNPhys->Draw("same");
efiEtaV0NTrk->Draw("same");
efiEtaV0NTrkNo->Draw("same");

c1->cd(5);
efiPhiSingPTrk->Draw();
//efiPhiSingPPhys->Draw("same");
efiPhiV0PTrk->Draw("same");
efiPhiV0PTrkNo->Draw("same");

c1->cd(6);
efiPhiSingNTrk->Draw();
//efiPhiSingNPhys->Draw("same");
efiPhiV0NTrk->Draw("same");
efiPhiV0NTrkNo->Draw("same");

c1->cd(7);
TH1* efiRSingTrk=(TH1*) rgeantTrk->Clone();
efiRSingTrk->Sumw2();
efiRSingTrk->Divide(rsingleTrk,rgeantTrk);

//efiRSingTrk->SetAxisRange(0.,5);
efiRSingTrk->SetMinimum(0.);
efiRSingTrk->SetMaximum(1.2);
efiRSingTrk->SetMarkerStyle(21);
efiRSingTrk->SetMarkerColor(4);

efiRSingTrk->Draw();

TH1* efiRV0Trk=(TH1*) rgeantTrk->Clone();
efiRV0Trk->Sumw2();
efiRV0Trk->Divide(rV0Trk,rsingleTrk);


efiRV0Trk->SetMinimum(0.);
efiRV0Trk->SetMaximum(1.2);
efiRV0Trk->SetMarkerStyle(21);
efiRV0Trk->SetMarkerColor(3);
efiRV0Trk->Draw("same");


c1->Print("Effi-norefit-trunk010709.ps"); 
//c1->Print("Effi-norefit-v416rev08-090709.ps"); 

}
