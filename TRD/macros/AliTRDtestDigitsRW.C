//
// Test digits
// Author: Mateusz Ploskon
//
void showLegend(TH1F *h, TH1F *h2)
{
  TLegend * leg = new TLegend(0, 0, 0.2, 0.2);
  leg->AddEntry(h, h->GetTitle(), "LBF");
  leg->AddEntry(h2, h2->GetTitle(), "LBF");
  leg->Draw();
}

void AliTRDtestDigitsRW(Int_t thresh = 0, Int_t maxDet = 540)
{
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  cout << endl;

  TFile *finD = TFile::Open("TRD.Digits.root");
  TTree *treeD = (TTree*)finD->Get("Event0/TreeD");

  AliTRDdigitsManager manD;
  manD.ReadDigits(treeD);

  AliTRDdigitsManager manR;
  manR.CreateArrays();
  
  AliRawReaderRoot reader("raw.root", 0);
  reader.Select("TRD");
  AliTRDRawStreamV2 stream(&reader);
  stream.SetRawVersion(3);

  Int_t ichambers = 0;
  while (stream.NextChamber(&manR) >= 0)
    ichambers++;

  cout << "Chambers loaded - stream V2 " << ichambers << endl;

  TH1F *hsignalD = new TH1F("hsignalD", "hsignalD - stream D", 1000, 0, 1000);
  hsignalD->SetFillColor(32);
  hsignalD->SetFillStyle(1001);
  TH1F *hsignalR = new TH1F("hsignalR", "hsignalR - stream V2", 1000, 0, 1000);
  hsignalR->SetLineWidth(2);
  hsignalR->SetLineColor(kRed);

  TH1F *htbinD = new TH1F("htbinD", "htbinD - stream D", 30, 0, 30);
  htbinD->SetFillColor(32);
  htbinD->SetFillStyle(1001);
  TH1F *htbinR = new TH1F("htbinR", "htbinR - stream V2", 30, 0, 30);
  htbinR->SetLineWidth(2);
  htbinR->SetLineColor(kRed);

  TH1F *hdetD = new TH1F("hdetD", "hdetD - stream D", 540, 0, 540);
  hdetD->SetFillColor(32);
  hdetD->SetFillStyle(1001);
  TH1F *hdetR = new TH1F("hdetR", "hdetR - stream V2", 540, 0, 540);
  hdetR->SetLineWidth(2);
  hdetR->SetLineColor(kRed);
  
  TH2F *h2D = new TH2F("h2D", "h2D - stream D", 16, 0, 16, 144, 0, 144 );
  TH2F *h2R = new TH2F("h2R", "h2R - stream V2", 16, 0, 16, 144, 0, 144 );

  TH2F *h2Drct[16]; // 1 for each row
  TH2F *h2Rrct[16]; // 1 for each row
  for (Int_t ir = 0; ir < 16; ir++)
    {
      char cname[255];
      char ctitle[255];
      sprintf(cname, "h2Drct_%d", ir);
      sprintf(ctitle, "D Row %d;col;tbin;value", ir);
      h2Drct[ir] = new TH2F(cname, ctitle, 144, 0, 144, 30, 0, 30);

      sprintf(cname, "h2Rrct_%d", ir);
      sprintf(ctitle, "R Row %d;col;tbin;value", ir);
      h2Rrct[ir] = new TH2F(cname, ctitle, 144, 0, 144, 30, 0, 30);
    }

  for (Int_t idet = 0; idet < maxDet; idet++)
    {
      digitsD = manD.GetDigits(idet);
      digitsD->Expand();
      digitsR = manR.GetDigits(idet);
      
      Int_t rowMax = digitsD->GetNrow();
      Int_t colMax = digitsD->GetNcol();
      Int_t timeMax = digitsD->GetNtime();	
      
      //cout << "\r Detector " << idet << endl;
      cout << "\r Detector " << idet; cout.flush();
      
      for (Int_t irow = 0; irow < rowMax; irow++)
	{
	  for (Int_t icol = 0; icol < colMax; icol++)
	    {
	      for (Int_t itime = 0; itime < timeMax; itime++)
		{
		  Int_t ivalD = digitsD->GetDataUnchecked(irow, icol, itime);
		  Int_t ivalR = digitsR->GetDataUnchecked(irow, icol, itime);

		  //if (ivalD > thresh && ivalR > thresh)
		  if (ivalD > thresh)
		    {
		      hsignalD->Fill(ivalD);		      
		      htbinD->Fill(itime, ivalD);		      
		      hdetD->Fill(idet, ivalD);
		      h2D->Fill(irow, icol, ivalD);
		      h2Drct[irow]->Fill(icol, itime, ivalD);
		    }

		  hsignalR->Fill(ivalR);
		  htbinR->Fill(itime, ivalR);
		  hdetR->Fill(idet, ivalR);
		  h2R->Fill(irow, icol, ivalR);
		  h2Rrct[irow]->Fill(icol, itime, ivalR);		  
		} //time
	    } //col
	} //row

      digitsD->Compress(1,0);
    }

  cout << endl;

  gStyle->SetPalette(1);

  TCanvas *tc = new TCanvas("tc_d", "tc_d");
  tc->Divide(2,3);
  tc->cd(1);
  hsignalR->Draw();
  hsignalD->Draw("same");
  hsignalR->Draw("same");
  gPad->SetLogy();
  showLegend(hsignalR, hsignalD);

  tc->cd(2);
  htbinR->Draw();
  htbinD->Draw("same");
  htbinR->Draw("same");
  gPad->SetLogy();

  tc->cd(3);
  hdetR->Draw();
  hdetD->Draw("same");
  hdetR->Draw("same");
  gPad->SetLogy();

  tc->cd(5);
  h2D->Draw("colz");
  tc->cd(6);
  h2R->Draw("colz");

  TCanvas *tc1 = new TCanvas("tc_d_rct_1", "tc_d_rct_1");
  tc1->Divide(4, 4);

  TCanvas *tc2 = new TCanvas("tc_d_rct_2", "tc_d_rct_2");
  tc2->Divide(4, 4);
  for (Int_t ir = 0; ir < 16; ir++)
    {
      if (ir < 8)
	{
	  tc1->cd(ir*2 + 1);
	  h2Drct[ir]->Draw("colz");
	  tc1->cd(ir*2 + 2);
	  h2Rrct[ir]->Draw("colz");
	}
      else
	{
	  tc2->cd((ir-8)*2 + 1);
	  h2Drct[ir]->Draw("colz");
	  tc2->cd((ir-8)*2 + 2);
	  h2Rrct[ir]->Draw("colz");
	}
    }

  //TCanvas *tc = new TCanvas("tc_d_rct", "tc_d_rct");
  //tc->Divide(4, 8);

//     for (Int_t ir = 0; ir < 16; ir++)
//     {
//           tc1->cd(ir + 1);
//           h2Drct[ir]->Draw("colz");

//           tc2->cd(ir + 1);
//           h2Rrct[ir]->Draw("colz");
//     }

}

///---------------------------------------------------------------
void testDigits2streams()
{
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  //AliCDBStorage* localStorage = cdb->GetStorage("local://$ALICE_ROOT/OCDB");
  //  cout << "[I] Is storage set : " << cdb->IsDefaultStorageSet() << endl;
  cout << endl;

//   TFile *finD = TFile::Open("TRD.Digits.root");
//   TTree *treeD = (TTree*)finD->Get("Event0/TreeD");

//   AliTRDdigitsManager manD;
//   manD.ReadDigits(treeD);

  //OLD stream
  AliTRDdigitsManager manD;
  manD.CreateArrays();

  AliRawReaderRoot reader1("raw.root", 0);
  reader1.Select("TRD");
  AliTRDRawStream stream1(&reader1);
  stream1.SetRawVersion(2);

  Int_t ichambers = 0;
  while (stream1.NextChamber(&manD) >= 0)
    ichambers++;

  cout << "Chambers loaded - stream V1" << ichambers << endl;

  // NEW STREAM
  AliTRDdigitsManager manR;
  manR.CreateArrays();
  
  AliRawReaderRoot reader("raw.root", 0);
  reader.Select("TRD");
  AliTRDRawStreamV2 stream(&reader);
  stream.SetRawVersion(2);

  Int_t ichambers = 0;
  while (stream.NextChamber(&manR) >= 0)
    ichambers++;

  cout << "Chambers loaded - stream V2" << ichambers << endl;

  AliTRDdataArrayI *digitsD = 0;
  AliTRDdataArrayI *digitsR = 0;

  TH1F *hsignalD = new TH1F("hsignalD", "hsignalD - stream V1", 1000, 0, 1000);
  hsignalD->SetFillColor(32);
  hsignalD->SetFillStyle(1001);
  TH1F *hsignalR = new TH1F("hsignalR", "hsignalR - stream V2", 1000, 0, 1000);
  hsignalR->SetLineWidth(2);
  hsignalR->SetLineColor(kRed);

  TH1F *htbinD = new TH1F("htbinD", "htbinD - stream V1", 30, 0, 30);
  htbinD->SetFillColor(32);
  htbinD->SetFillStyle(1001);
  TH1F *htbinR = new TH1F("htbinR", "htbinR - stream V2", 30, 0, 30);
  htbinR->SetLineWidth(2);
  htbinR->SetLineColor(kRed);
  
  TH2F *h2D = new TH2F("h2D", "h2D - stream V1", 16, 0, 16, 144, 0, 144 );
  TH2F *h2R = new TH2F("h2R", "h2R - stream V2", 16, 0, 16, 144, 0, 144 );
  
  //for (Int_t idet = 0; idet < 540; idet++)
  for (Int_t idet = 0; idet < 10; idet++)
  //for (Int_t idet = 0; idet < 10; idet++)
    {
      digitsD = manD.GetDigits(idet);
      digitsD->Expand();
      digitsR = manR.GetDigits(idet);
      
      Int_t rowMax = digitsD->GetNrow();
      Int_t colMax = digitsD->GetNcol();
      Int_t timeMax = digitsD->GetNtime();	
      
      //cout << "\r Detector " << idet << endl;
      cout << "\r Detector " << idet; cout.flush();
      
      for (Int_t irow = 0; irow < rowMax; irow++)
	{
	  for (Int_t icol = 0; icol < colMax; icol++)
	    {
	      for (Int_t itime = 0; itime < timeMax; itime++)
		{
		  Int_t ivalD = digitsD->GetDataUnchecked(irow, icol, itime);
		  Int_t ivalR = digitsR->GetDataUnchecked(irow, icol, itime);
		  
		  hsignalD->Fill(ivalD);
		  hsignalR->Fill(ivalR);

		  htbinD->Fill(itime, ivalD);
		  htbinR->Fill(itime, ivalR);

		  h2D->Fill(irow, icol, ivalD);
		  h2R->Fill(irow, icol, ivalR);
		  
		} //time
	    } //col
	} //row

      digitsD->Compress(1,0);
    }

  cout << endl;

  gStyle->SetPalette(1);

  TCanvas *tc = new TCanvas("tc", "tc");
  tc->Divide(2,2);
  tc->cd(1);
  hsignalR->Draw();
  hsignalD->Draw("same");
  hsignalR->Draw("same");
  gPad->SetLogy();
  showLegend(hsignalR, hsignalD);

  tc->cd(2);
  htbinR->Draw();
  htbinD->Draw("same");
  htbinR->Draw("same");
  gPad->SetLogy();

  tc->cd(3);
  h2D->Draw("colz");
  tc->cd(4);
  h2R->Draw("colz");
}

///---------------------------------------------------------------
void printDigits()
{
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  //AliCDBStorage* localStorage = cdb->GetStorage("local://$ALICE_ROOT/OCDB");
  //  cout << "[I] Is storage set : " << cdb->IsDefaultStorageSet() << endl;
  cout << endl;

  TFile *finD = TFile::Open("TRD.Digits.root");
  TTree *treeD = (TTree*)finD->Get("Event0/TreeD");

  AliTRDdigitsManager manD;
  manD.SetRaw();
  manD.ReadDigits(treeD);

  AliTRDdataArrayI *digitsD = 0;

  for (Int_t idet = 0; idet < 540; idet++)
    {
      digitsD = manD.GetDigits(idet);
      digitsD.Expand();
      
      Int_t rowMax = digitsD->GetNrow();
      Int_t colMax = digitsD->GetNcol();
      Int_t timeMax = digitsD->GetNtime();	
      
      cout << "Detector " << idet << " nrows " << rowMax << endl;

      for (Int_t irow = 0; irow < rowMax; irow++)
	{
	  for (Int_t icol = 0; icol < colMax; icol++)
	    {
	      for (Int_t itime = 0; itime < timeMax; itime++)
		{
		  Int_t ivalD = digitsD->GetDataUnchecked(irow, icol, itime);
		  cout << Form("det %d rct %d %d %d valD=%d", idet, irow, icol, itime, ivalD) << endl;
		}
	    }
	}

      digitsD->Compress(1,0);
    }
}
