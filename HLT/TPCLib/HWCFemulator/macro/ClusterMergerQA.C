/****************************************************************************
 *           Very important, delicate and rather obscure macro.             *
 *                                                                          *
 *               Creates list of "trackable" tracks,                        *
 *             calculates efficiency, resolutions etc.                      *
 *     There is a possibility to run this macro over several events.         *
 *                                                                          *
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 * with several nice improvements by: M.Ivanov, GSI, m.ivanov@gsi.de        *
 ****************************************************************************/

#if defined(__CINT__) && !defined(__MAKECINT__)
{
  
  gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/TPC -I$ALICE/geant3/TGeant3 -I$ALICE_ROOT/HLT/TPCLib/HWCFemulator -I$ALICE_ROOT/HLT/TPCLib -I$ALICE_ROOT/HLT/BASE ");    
  AliHLTPluginBase::GetInstance()->LoadComponentLibraries("libAliHLTTPC.so");    
  gROOT->LoadMacro("ClusterMergerQA.C++");

  runTest("NoMerge/TPC.RecPoints.root");

}
#else

  #include <TMath.h>
  #include <TError.h>
  #include <Riostream.h>
  #include <TH1F.h>
  #include <TH2F.h>
  #include <TTree.h>
  #include <TParticle.h>
  #include <TCanvas.h>
  #include <TLine.h>
  #include <TText.h>
  #include <TStyle.h>
  #include <TFile.h>
  #include <TROOT.h>
  #include "TString.h"
  #include "AliTPCclusterMI.h"

  #include "AliStack.h"
  #include "AliHeader.h"
  #include "AliTrackReference.h"
  #include "AliRunLoader.h"
  #include "AliRun.h"
  #include "AliTPCtrack.h"
  #include "AliTracker.h"
  
  #include "AliTPC.h"
  #include "AliTPCClustersArray.h"
  #include "AliTPCClustersRow.h"
  #include "AliTPCcluster.h"
  #include "AliTPCLoader.h"
  #include "TParticlePDG.h"
  #include "TDatabasePDG.h"
  #include "AliGeomManager.h"
  #include <iostream>
  #include <string>
  #include "AliCDBManager.h"
  #include "AliHLTTPCHWClusterMerger.h"
  #include "AliHLTTPCTransform.h"

extern AliRun *gAlice;
extern TROOT *gROOT;


void myDraw(TH1F *h)
{
  Int_t minc=33; 
  if (h->GetEntries()<minc) h->Draw(); else h->Fit("gaus","Q");
}

Int_t runTest(const char *recPointsName ) 
{

  // access clusters
  
  //const Char_t* recPointsName = "NoMerge/TPC.RecPoints.root";  

 
  AliRunLoader *rl = AliRunLoader::Open("galice.root","COMPARISON");
  if (!rl) return 1;
 
  rl->LoadgAlice();
  rl->LoadHeader();

  AliCDBManager *man=AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  //man->SetDefaultStorage("local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Residual/");
  man->SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));
  man->SetRun(0);
  //man->SetRun(137366);
  
  if(AliGeomManager::GetGeometry()==NULL){
    AliGeomManager::LoadGeometry();
  }

  AliGeomManager::ApplyAlignObjsFromCDB("TPC");
  
  //TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG, AliMagF::kBeamTypeAA, 1300.));
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));
 
  // Check field
  if (!TGeoGlobalMagField::Instance()) {
    ::Error("","magnetic field not initialized, please set up TGeoGlobalMagField and AliMagF");
    return 1;
  }
  

  Int_t nev=rl->GetNumberOfEvents();
 
  // ******* Loop over events *********


  AliMagF *field = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();

  if(! field) {
    ::Error("","No mag. field found!!!");
    return 1;
  }

  AliHLTTPCHWClusterMerger merger;
  merger.Init();

  int statNGood=0, statNMerged=0;

  //nev = 1;
  // ********  Loop over generated events 
  for (Int_t nEv=0; nEv<nev; nEv++) {
    cout<<"Event " <<nEv<<endl;
    
    merger.Clear();
    static int nClusters =0;
    static int nAccepted = 0;

    vector<AliTPCclusterMI> vClusters;

    // Read clusters
    { 
      rl->GetEvent(nEv);
    
      TFile* file = new TFile(recPointsName, "READ");
      {
	TString str = "Event";
	str+= nEv; 
	file->cd(str.Data());
      }
      
      TTree* tRecPoints = dynamic_cast<TTree*> (gDirectory->Get("TreeR"));
      
      AliTPCClustersRow *row = new AliTPCClustersRow();
      
      tRecPoints->GetBranch("Segment")->SetAddress(&row);

      for (Int_t iRow = 0; iRow < tRecPoints->GetEntriesFast(); iRow++) {
	tRecPoints->GetEntry(iRow);             
	for (Int_t nCl = 0; nCl < row->GetArray()->GetEntriesFast(); nCl++){
	  	  
	  AliTPCclusterMI * cl = static_cast<AliTPCclusterMI*> (row->GetArray()->UncheckedAt(nCl));
	  
	  int slice=-1, padRow=-1;
	  	  
	  if( !AliHLTTPCTransform::Sector2Slice(slice, padRow, cl->GetDetector(), cl->GetRow() ) ) continue;	  
	  int patch = AliHLTTPCTransform::GetPatch(padRow);
	  padRow -=AliHLTTPCTransform::GetFirstRow(patch);

	  double yW = AliHLTTPCTransform::GetPadPitchWidth(patch);
	  double zW = AliHLTTPCTransform::GetZWidth();

	  AliHLTTPCRawCluster rawCluster;
	  rawCluster.SetPadRow(padRow);
	  rawCluster.SetPad(cl->GetPad());
	  rawCluster.SetTime(cl->GetTimeBin());

	  rawCluster.SetSigmaY2(cl->GetSigmaY2()/yW/yW);
	  rawCluster.SetSigmaZ2(cl->GetSigmaZ2()/zW/zW);
	  rawCluster.SetCharge((UShort_t)cl->GetQ()); 
	  rawCluster.SetQMax((UShort_t)cl->GetMax());
	  
	  AliHLTTPCClusterMCLabel mc;
	  int l0 = ( cl->GetLabel(0)>=0 );
	  int l1 = ( cl->GetLabel(1)>=0 );
	  int l2 = ( cl->GetLabel(2)>=0 );

	  int nmc = l0+l1+l2;
	  float w = (nmc>0) ? cl->GetQ()/nmc :0;

	  mc.fClusterID[0].fMCID = cl->GetLabel(0);  mc.fClusterID[0].fWeight=l0*w;
	  mc.fClusterID[1].fMCID = cl->GetLabel(1);  mc.fClusterID[0].fWeight=l1*w;
	  mc.fClusterID[2].fMCID = cl->GetLabel(2);  mc.fClusterID[0].fWeight=l2*w;

	  nClusters++;

	  if( !merger.CheckCandidate(slice, patch, rawCluster) ) continue;
	  
	  
	  int id = merger.AddCandidate(slice, patch, ~AliHLTUInt32_t(0), rawCluster, &mc);
	  if( id>=0 ){
	    nAccepted++;
	    vClusters.push_back(*cl);
	  }

	  //if( (id>=0) &&(merger.GetRecords()[id].GetBorder()>=0) ) nAccepted++;
	  
	}
      } // loop over clusters                    
  
      delete row;
      file->Close();  
    }


    cout<<"Merge..."<<endl;

    static int nMerged = 0;
    nMerged+=merger.Merge();

    cout<<"Read "<<nClusters<<" clusters, accepted "<<nAccepted<<" clusters "<<endl;
    cout<<"Merged "<<2*nMerged<<" clusters "
	<<" ("<<100.*2*nMerged/nAccepted<<" % of "<<nAccepted<< " accepted border clusters )"<<endl;

    for( int iB=0; iB<merger.GetNSlices()*merger.GetNBorders(); iB+=2 ){
       const AliHLTTPCHWClusterMerger::AliBorderRecord 
	 *b1 = merger.GetBorderClusters() + merger.GetBorderFirstCluster(iB),
	 *b2 = merger.GetBorderClusters() + merger.GetBorderFirstCluster(iB+1);
       int 
	 n1 = merger.GetBorderNClusters(iB),
	 n2 = merger.GetBorderNClusters(iB+1);
       if( n1<=0 || n2<=0 ) continue;

       //cout<<iB<<" "<<n1<<" "<<n2<<endl;
       //cout<<iB<<" "<<merger.GetBorderFirstCluster(iB)<<" "<<n1<<endl;
       //cout<<iB+1<<" "<<merger.GetBorderFirstCluster(iB+1)<<" "<<n2<<endl;
       /*
       cout<<"Border "<<iB<<":"<<endl;
       cout<<"  left: "<<n1<<" : ";
       for( int i=0; i<n1; i++ ){
	 const AliTPCclusterMI &c = vClusters[b1[i].fClusterRecordID];
	 cout<<" ("<<c.GetLabel(0)<<" "<<(int)c.GetTimeBin()<<") ";
       }
       cout<<endl;
       cout<<"  right: "<<n2<<" : ";
       for( int i=0; i<n2; i++ ){
	 const AliTPCclusterMI &c = vClusters[b2[i].fClusterRecordID];
	 cout<<" ("<<c.GetLabel(0)<<" "<<(int)c.GetTimeBin()<<") ";
       }
       cout<<endl;
       */

       for( int i1=0; i1<n1; i1++ ){
	 const AliTPCclusterMI &c1 = vClusters[b1[i1].fClusterRecordID];
	 int lab = c1.GetLabel(0);
	 if( lab<0 ) continue;

	 // check if the label has been already checked
	 bool b=0;
	 for( int i=0; i<i1; i++ ){
	   const AliTPCclusterMI &c = vClusters[b1[i].fClusterRecordID];
	   if( c.GetLabel(0) == lab ){
	     b = 1;
	     break;
	   }
	 }
	 if( b ) continue;

	 // check if there is a partner with the same label
	 bool isPartner = 0;
	 bool isMerged = 0;
	 for( int i2=0; i2<n2; i2++ ){
	   const AliTPCclusterMI &c2 = vClusters[b2[i2].fClusterRecordID];
	   if( c2.GetLabel(0) != lab ) continue;
	   if( fabs(c2.GetTimeBin()-c1.GetTimeBin() ) > 5 ) continue;
 
	   isPartner = 1 ;
	   // check if it has been a merged with a cluster with the same label	     
	   int id = merger.GetRecords()[b2[i2].fClusterRecordID].IsMergedTo();
	   if( id<0 ) continue;
	   const AliTPCclusterMI &c21 = vClusters[id];
	   if( c21.GetLabel(0) == lab ) isMerged = 1;	 
	 }
	 if( isPartner ) statNGood++;
	 if( isMerged ) statNMerged++;
       }
       /*
       cout<<" Merging efficiency: merged "<<statNMerged<<" of "<<statNGood<<" tracks : "
	   <<100.*statNMerged/( statNGood>0 ?statNGood :1)<<" %"<<endl;
       */ 
    }

       cout<<" Merging efficiency: merged "<<statNMerged<<" of "<<statNGood<<" tracks : "
	   <<100.*statNMerged/( statNGood>0 ?statNGood :1)<<" %"<<endl;
 
  }// ***** End of the loop over events
 
  delete rl;
  

  /*
  gStyle->SetOptStat(111110);
  gStyle->SetOptFit(1);   
  //hphidist->Draw();
  Int_t minc=33; 
  {
    TCanvas *c1=new TCanvas("c0","Offline vs HLT clusters",0,0,1500,800);  
    c1->Divide(1,3);
    c1->cd(1);
    hQ->Draw();
    c1->cd(2);
    hQ10->Draw();
    c1->cd(3);
    hQ30->Draw();
    c1->cd(0);
    c1->Update();
  }
  */  
  return 0;
}

#endif
