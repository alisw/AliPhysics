// $Id$
//
// Test macro for reading motif type data and iterate over them.

void testMotifTypeIterators(AliMpStationType station = kStation1,
                            AliMpPlaneType plane = kBendingPlane)
{
  TString names="ABCDEFGHI";
  //TString names="FEG";
    
  TH2C* histos[] = new TH2C* [names.Length()];
  TCanvas* canv[] = new TCanvas* [1+(names.Length()-1)/4];
  Int_t i;
  for (i=0;i<1+(names.Length()-1)/4;++i){
    // canv[i] = new TCanvas(Form("canv%d",i),"Iterator viewing...");
               // CINT limitation on DEC
	       
    TString cname("canv"); cname += i;
    canv[i] = new TCanvas(cname.Data(),"Iterator viewing...");
    
    canv[i]->Divide(2,2); 
  }
    
  AliMpReader r(station, plane);
  //r.SetVerboseLevel(2);

  for (i=0;i<names.Length();++i){
    AliMpMotifType *mt = r.BuildMotifType(names[i]);
    canv[i/4]->cd(1+ (i%4));
    //histos[i] = new TH2C(Form("h%d",i),Form("Motif type %c",names[i]),
    //                     mt->GetNofPadsX(),-0.5,mt->GetNofPadsX()-0.5,
    //                     mt->GetNofPadsY(),-0.5,mt->GetNofPadsY()-0.5);
               // CINT limitation on DEC

    TString hname("h"); hname += i;
    TString mname = names(i,1);	       

    histos[i] = new TH2C(hname.Data(), mname.Data(),
                         mt->GetNofPadsX(),-0.5,mt->GetNofPadsX()-0.5,
                         mt->GetNofPadsY(),-0.5,mt->GetNofPadsY()-0.5);

    cout<<"Motif Type "<<mt->GetID()<<endl;
    cout<<"--------------------------------"<<endl;
    Int_t num=0;

    AliMpMotifTypePadIterator it = AliMpMotifTypePadIterator(mt);

    for (it.First(); ! it.IsDone(); it.Next()) {
      cout << "Iterator " << num << ' '<< it.CurrentItem().GetIndices() << endl;
      ++num;
      histos[i]->Fill(it.CurrentItem().GetIndices().GetFirst(),
                      it.CurrentItem().GetIndices().GetSecond(),num);
    }

    delete mt;

    histos[i]->Draw("text");
    canv[i/4]->Update();
    
  }
}
