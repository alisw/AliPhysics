// $Id$
// $MpId: testMotifTypeIterators.C,v 1.11 2005/09/26 16:05:25 ivana Exp $
//
// Test macro for reading motif type data and iterate over them.

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpMotifReader.h"
#include "AliMpMotifType.h"
#include "AliMpMotifTypePadIterator.h"
#include "AliMpEncodePair.h"

#include <Riostream.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TH2.h>

#include <vector>

#endif

TCanvas* CreateTCanvas(const TString& name, const TString& title,
                       AliMq::Station12Type station, AliMp::PlaneType plane)
{
  TString newName(name);
  TString newTitle(title);
  TString unique = AliMq::Station12TypeName(station) + AliMp::PlaneTypeName(plane);
  newName += unique;
  newTitle += unique;
  return new TCanvas(newName.Data(), newTitle.Data());
}                     

void testMotifTypeIterators(AliMq::Station12Type station, AliMp::PlaneType plane)
{
  TString names;
  TString names2;
  Int_t nv =0;
  if ( station == AliMq::kStation1 ) {
    if ( plane == AliMp::kBendingPlane ) 
      names ="ABCDEFGHI";
    else
      names = "ABCDEFGHIJKLMN";
  }    
  else if ( station == AliMq::kStation2 ) {
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
  }  
  Int_t nofMotifs = names.Length() + names2.Length(); 
  // cout << " nofMotifs: " << nofMotifs << endl;   
    
  std::vector<TH2C*> histos;
  std::vector<TCanvas*> cvs;
  Int_t i;
  for (i=0;i<1+(nofMotifs-1)/4;++i){
    TString cname("canv"); cname += i;
    TCanvas* canv = CreateTCanvas(cname.Data(),"Iterator viewing...", station, plane);
    canv->Divide(2,2); 
    cvs.push_back(canv);
  }
    
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpMotifReader r(dataStreams, AliMp::kStation12, station, plane);
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

    cvs[i/4]->cd(1+ (i%4));
    //histos[i] = new TH2C(Form("h%d",i),Form("Motif type %c",names[i]),
    //                     mt->GetNofPadsX(),-0.5,mt->GetNofPadsX()-0.5,
    //                     mt->GetNofPadsY(),-0.5,mt->GetNofPadsY()-0.5);
               // CINT limitation on DEC

    TString hname("h"); hname += i;

    TH2C* histo = new TH2C(hname.Data(), mname.Data(),
                         mt->GetNofPadsX(),-0.5,mt->GetNofPadsX()-0.5,
                         mt->GetNofPadsY(),-0.5,mt->GetNofPadsY()-0.5);
    histos.push_back(histo);                         

    cout<<"Motif Type "<<mt->GetID()<<endl;
    cout<<"--------------------------------"<<endl;
    Int_t num=0;

    AliMpMotifTypePadIterator it = AliMpMotifTypePadIterator(mt);

    for (it.First(); ! it.IsDone(); it.Next()) {
      cout << "Iterator " << num << ' ';
      AliMp::PairPut(cout, it.CurrentItem().GetIndices()) << endl;
      ++num;
      histo->Fill(it.CurrentItem().GetIx(), it.CurrentItem().GetIy(),num);
    }

    //delete mt;
    histo->Draw("text");
    cvs[i/4]->Update();
  }
}
void testSt12MotifTypeIterators()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testMotifTypeIterators for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testMotifTypeIterators(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
