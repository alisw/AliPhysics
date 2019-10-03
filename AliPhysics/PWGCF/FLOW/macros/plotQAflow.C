
void plotQAflow(const char* filename="")
{
 TFile f(filename,"read");
 TKey* key = f.GetListOfKeys()->At(0);
 if (!key) return;
 TObjArray* flowQA = dynamic_cast<TObjArray*>(f.Get(key->GetName()));
 if (!flowQA) return;
 TObjArray* before = dynamic_cast<TObjArray*>(flowQA->At(0));
 TObjArray* after = dynamic_cast<TObjArray*>(flowQA->At(1));
 for (Int_t i=0; i<before->GetEntries(); i++)
 {
   TH1* hbefore = dynamic_cast<TH1*>(before->At(i));
   TH1* hafter = dynamic_cast<TH1*>(after->At(i));
   TCanvas* canvas = new TCanvas(hbefore->GetName(), hbefore->GetTitle());
   canvas->SetLogy();
   TLegend* legend = new TLegend(0.8,0.8,1.0,1.0);
   hbefore->SetAxisRange(0.1,hbefore->GetBinContent(hbefore->GetMaximumBin()),"Y");
   hbefore->Draw();
   legend->AddEntry(hbefore,"before cuts","l");
   hafter->SetLineColor(2);
   legend->AddEntry(hafter,"after cuts","l");
   hafter->Draw("same");
   legend->Draw();
 }
}
