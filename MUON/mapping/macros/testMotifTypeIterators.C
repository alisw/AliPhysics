// $Id$
// $MpId: testMotifTypeIterators.C,v 1.11 2005/09/26 16:05:25 ivana Exp $
//
// Test macro for reading motif type data and iterate over them.

void testMotifTypeIterators(AliMp::StationType station = AliMp::kStation1,
                            AliMp::PlaneType plane = AliMp::kBendingPlane,
	        	    Bool_t rootInput = false)
{
  TString names;
  TString names2;
  Int_t nv =0;
  if ( station == AliMp::kStation1 )
    if ( plane == AliMp::kBendingPlane ) 
      names ="ABCDEFGHI";
    else
      names = "ABCDEFGHIJKLMN";
  else if ( station == AliMp::kStation2 ) 
    if ( plane == AliMp::kBendingPlane ) {
      names ="ABCDEFGHIJKLMNOPQRSTUVWXY";
      names2 ="abcdefghimnptuvvvvv";
      nv = 5;
    }  
    else {
      names = "ABCEFGHIJKLMN";
      names2 ="abcdefgijklmnopqrstuwvvvvv";
      nv = 5;
    }  
  Int_t nofMotifs = names.Length() + names2.Length(); 
  // cout << " nofMotifs: " << nofMotifs << endl;   
    
  TH2C* histos[] = new TH2C* [nofMotifs];
  TCanvas* canv[] = new TCanvas* [1+(nofMotifs-1)/4];
  Int_t i;
  for (i=0;i<1+(nofMotifs-1)/4;++i){
    // canv[i] = new TCanvas(Form("canv%d",i),"Iterator viewing...");
               // CINT limitation on DEC
	       
    TString cname("canv"); cname += i;
    canv[i] = new TCanvas(cname.Data(),"Iterator viewing...");
    
    canv[i]->Divide(2,2); 
  }
    
  AliMpMotifReader r(station, plane);
  //r.SetVerboseLevel(2);

  for (i=0;i<nofMotifs;++i){

    // Get motif name
    TString mname;
    if (i<names.Length())
      mname = names(i, 1);
    else {
      mname = names2(i-names.Length(), 1); 
      if (mname == "v")
        mname += i - names.Length() - (names2.Length()-nv-1);
      else 	   
        mname += "1";
    }	
    //if (i==36) continue;  
        // break for these motifs (St2, BP) - to be investigated
   
    AliMpMotifType *mt = r.BuildMotifType(mname);

    canv[i/4]->cd(1+ (i%4));
    //histos[i] = new TH2C(Form("h%d",i),Form("Motif type %c",names[i]),
    //                     mt->GetNofPadsX(),-0.5,mt->GetNofPadsX()-0.5,
    //                     mt->GetNofPadsY(),-0.5,mt->GetNofPadsY()-0.5);
               // CINT limitation on DEC

    TString hname("h"); hname += i;

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

    //delete mt;
    histos[i]->Draw("text");
    canv[i/4]->Update();
  }
}
