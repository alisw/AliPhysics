/**
 * >> Testing Macro to compare FlatESDEvents from output files <<
 **
 * Primary Authors : Steffen Weber
 *
 * Usage:
 *  aliroot -b -l -q LoadLibs.C CompareFlatESDs.C++
 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "./AliFlatESDEvent.h"
#include "./AliFlatESDTrack.h"
#include "./AliFlatTPCCluster.h"
#include "./AliFlatExternalTrackParam.h"
#include "Riostream.h"
#endif   

void CompareFlatESDs(const char* filename1="outFlatESD1.dat",const char* filename2="outFlatESD2.dat", Bool_t verbose=kFALSE) {
  
  
  // Create output histograms
  
  
  TString outputFilename = "$PWD/compare.root";
  
	cout<< "creating histograms"<<endl;
	TH2F* hNTracks = new TH2F("nTracks","number of tracks", 100,0,100, 100,0,100);
	TH2F* hNV0s = new TH2F("nV0s","number of V0s", 10,0,10, 10,0,10);
	TH2F* hVtxTr = new TH2F("vtxTr","vtx Tracks", 2,0,2, 2,0,2);
	TH2F* hVtxSPD = new TH2F("vtxSPD","vtx SPD", 2,0,2, 2,0,2);
	TH1F * hStat = new TH1F("stat","statistics Differences",20,0,20);
	hStat->GetXaxis()->SetBinLabel(1,"All events");
	hStat->GetXaxis()->SetBinLabel(2,"no diffs");
	hStat->GetXaxis()->SetBinLabel(3,"nTracks");
	hStat->GetXaxis()->SetBinLabel(4,"nV0s");
	hStat->GetXaxis()->SetBinLabel(5,"vtxTracks");
	hStat->GetXaxis()->SetBinLabel(6,"vtxSPD");
	hStat->GetXaxis()->SetBinLabel(7,"tracks->extParams");

  
  

  ifstream is1(filename1, std::ifstream::binary | std::ifstream::in);
  ifstream is2(filename2, std::ifstream::binary | std::ifstream::in);
  if (is1 && is2 ){
    is1.seekg (0, is1.end);
    int length1 = is1.tellg();
    is1.seekg (0, is1.beg);
    char * buffer1 = new char [length1];
    
    std::cout << "Reading " << length1 << " characters... ";
    
    is1.read (buffer1,length1);
    if (is1)
      std::cout << "all characters read successfully." << endl;
    else
      std::cout << "error: only " << is1.gcount() << " could be read";
    is1.close();
	
	
    is2.seekg (0, is2.end);
    int length2 = is2.tellg();
    is2.seekg (0, is2.beg);
    char * buffer2 = new char [length2];
    
    std::cout << "Reading " << length2 << " characters... ";
    
    is2.read (buffer2,length2);
    if (is2)
      std::cout << "all characters read successfully." << endl;
    else
      std::cout << "error: only " << is2.gcount() << " could be read";
    is2.close();
    
	
	
    // ...buffer contains the entire file...
    
    char *curr1 = buffer1;
    char *endBuff1 = buffer1+length1;
	
    char *curr2 = buffer2;
    char *endBuff2 = buffer2+length2;
	
    int iEvent = 0;
	static const int nExt = 4;
    
	while( curr1 < endBuff1  && curr2 < endBuff2 ){
//cout<<" curr1 endBuff1 curr2 endBuff2 "<< static_cast<void*> (curr1)<<" "<<static_cast<void*> (endBuff1)<<" "<<static_cast<void*> (curr2)<<" "<<static_cast<void*> (endBuff2)<<endl;

	Int_t diff[6]={0};
      AliFlatESDEvent *flatEsd1 = reinterpret_cast<AliFlatESDEvent *>(curr1);
      AliFlatESDEvent *flatEsd2 = reinterpret_cast<AliFlatESDEvent *>(curr2);
	  /*
	  if(flatEsd1->GetNumberOfTracks()==0  ||  flatEsd2->GetNumberOfTracks() ==0){
		if(flatEsd1->GetNumberOfTracks()==0){
		  curr1=curr1+ flatEsd1->GetSize();
		}
		if(flatEsd2->GetNumberOfTracks()==0){
		  curr2=curr2+ flatEsd2->GetSize();
		}
		continue;
	  }
	  */
	  
      cout<<endl<<"Reading event "<<iEvent<<":"<<endl;
	  /*
	if(  flatEsd1->GetNumberOfTracks() != flatEsd2->GetNumberOfTracks() ) {
		cout<<"\t\tDIFFERENCE!: ";
		diff[0] =1;
	}
	cout<<"\t\tntracks: "<<flatEsd1->GetNumberOfTracks()<< " | " <<flatEsd2->GetNumberOfTracks();
	hNTracks->Fill(flatEsd1->GetNumberOfTracks(),flatEsd2->GetNumberOfTracks());
	  
	  if(  flatEsd1->GetNumberOfV0s() != flatEsd2->GetNumberOfV0s() ){
		cout<<"\t\tDIFFERENCE!: ";
		diff[1] =1;
	}
	  cout<<"\t\tnV0's: "<<flatEsd1->GetNumberOfV0s()<< " | " <<flatEsd2->GetNumberOfV0s()<<endl;
	  hNV0s->Fill(flatEsd1->GetNumberOfV0s(),flatEsd2->GetNumberOfV0s());
	 
	  if( (Bool_t) flatEsd1->GetPrimaryVertexTracks() != (Bool_t) flatEsd2->GetPrimaryVertexTracks()  ){
		cout<<"\t\tDIFFERENCE!: ";
		diff[2] =1;
	}


	  cout<<"\t\tvtx tracks: "<<(Bool_t) flatEsd1->GetPrimaryVertexTracks()<< " | " <<	(Bool_t) flatEsd2->GetPrimaryVertexTracks();	  
	  hVtxTr->Fill( (Bool_t) flatEsd1->GetPrimaryVertexTracks(), (Bool_t) flatEsd2->GetPrimaryVertexTracks());

	 
	  
	  if( (Bool_t) flatEsd1->GetPrimaryVertexSPD() != (Bool_t) flatEsd2->GetPrimaryVertexSPD()  ){
		cout<<"\t\tDIFFERENCE!: ";
		diff[3] =1;
	}
      cout<<"\t\tvtx SPD: "<<(Bool_t) flatEsd1->GetPrimaryVertexSPD() << " | " << (Bool_t) flatEsd2->GetPrimaryVertexSPD()<<endl;
	  hVtxSPD->Fill( (Bool_t) flatEsd1->GetPrimaryVertexSPD(), (Bool_t) flatEsd2->GetPrimaryVertexSPD());

  if(true|| (Bool_t)flatEsd1->GetPrimaryVertexSPD() && (Bool_t)flatEsd2->GetPrimaryVertexSPD()  ){
 		cout<<endl<<"\t\tvtx tracksX: "<< flatEsd1->GetPrimaryVertexSPD()<<" | " <<	flatEsd2->GetPrimaryVertexSPD();	  
		}
	  */
	  
	  // compare tracks
	  #if 1
	  AliFlatESDTrack *track1 = flatEsd1->GetTracks();
	  AliFlatESDTrack *track2 = flatEsd2->GetTracks();
    for (Int_t idxTrack = 0; idxTrack < flatEsd1->GetNumberOfTracks() && track1 && track2; ++idxTrack) { 

		AliFlatExternalTrackParam* ext[2][nExt] ={
			{
				track1->GetTrackParamRefitted(),
				track1->GetTrackParamIp(),
				track1->GetTrackParamTPCInner(),
				track1->GetTrackParamOp(),
		//		track1->GetTrackParamCp(),
		//		track1->GetTrackParamITSOut()
			},
			{
				track2->GetTrackParamRefitted(),
				track2->GetTrackParamIp(),
				track2->GetTrackParamTPCInner(),
				track2->GetTrackParamOp(),
			//	track2->GetTrackParamCp(),
			//	track2->GetTrackParamITSOut()
			}
		};
	
     	//Printf("  TEST: FlatTrack1 %d > FlatExternalTrackParam1 > %p %p %p %p", idxTrack, exp11, exp21, exp31, exp41);
     	//Printf("  TEST: FlatTrack2 %d > FlatExternalTrackParam2 > %p %p %p %p", idxTrack, exp12, exp22, exp32, exp42);


	for(int iExt=0; iExt<nExt; ++iExt){
cout<<endl<<iExt<<endl;		
if(!ext[0][iExt] && !ext[1][iExt]) continue;	
		if(!ext[0][iExt] && ext[1][iExt]){
		//	cout<<"DIFFERENCE!: ";
	 		cout<<" ext"<<iExt<<" not set in "<<filename1<<endl;
		}	
		if(ext[0][iExt] && !ext[1][iExt]){
		//	cout<<"DIFFERENCE!: ";
	 		cout<<" ext"<<iExt<<" not set in "<<filename2<<endl;
		}


		if( (!ext[0][iExt] || !ext[1][iExt])|| ext[0][iExt]->GetAlpha() != ext[1][iExt]->GetAlpha() ) {
	//		cout<<"DIFFERENCE!: ";
	 		//cout<<" alpha"<<iExt<<" :"  << (ext[0][iExt] ? ext[0][iExt]->GetAlpha() : -99.)  << "\t\t" << (ext[1][iExt] ?  ext[1][iExt]->GetAlpha(): -99.)<<endl;
			diff[4]=1;
		}	cout<<" alpha"<<iExt<<" :"  << (ext[0][iExt] ? ext[0][iExt]->GetAlpha() : -99.)  << "\t\t" << (ext[1][iExt] ?  ext[1][iExt]->GetAlpha(): -99.)<<endl;
			

		if( (!ext[0][iExt] || !ext[1][iExt])||ext[0][iExt]->GetX() != ext[1][iExt]->GetX() ) {
			//cout<<"DIFFERENCE!: ";
	 		//cout<<" GetX"<<iExt<<" :"  << (ext[0][iExt] ? ext[0][iExt]->GetX(): -99.)  << " | " << (ext[1][iExt] ?  ext[1][iExt]->GetX(): -99.)<<endl;
			diff[4]=1;
		}	
cout<<" GetX"<<iExt<<" :"  << (ext[0][iExt] ? ext[0][iExt]->GetX(): -99.)  << " | " << (ext[1][iExt] ?  ext[1][iExt]->GetX(): -99.)<<endl;


		if( (!ext[0][iExt] || !ext[1][iExt])||ext[0][iExt]->GetSigned1Pt() !=  ext[0][iExt]->GetSigned1Pt() ) {
			//cout<<"DIFFERENCE!: ";
	 		//cout<<" 1/pt"<<iExt<<" :"  <<  (ext[0][iExt] ? ext[0][iExt]->GetSigned1Pt(): -99.)  << " | " << (ext[1][iExt] ?  ext[1][iExt]->GetSigned1Pt(): -99.)<<endl;
			diff[4]=1;
		}	
	cout<<" 1/pt"<<iExt<<" :"  <<  (ext[0][iExt] ? ext[0][iExt]->GetSigned1Pt(): -99.)  << " | " << (ext[1][iExt] ?  ext[1][iExt]->GetSigned1Pt(): -99.)<<endl;
			

}
	
	  
	  /*
	if( track1->GetNumberOfTPCClusters() != track2->GetNumberOfTPCClusters() ){
		cout<<"DIFFERENCE!: ";
		cout<<" nTPCclusters: "<<track1->GetNumberOfTPCClusters()<< " | " <<track2->GetNumberOfTPCClusters()<< endl;
		diff[4]=1;
	}  
	if( track1->GetNumberOfITSClusters() != track2->GetNumberOfITSClusters() ){
		cout<<"DIFFERENCE!: ";
	 	cout<<" nITSclusters: "<<track1->GetNumberOfITSClusters()<< " | " <<track2->GetNumberOfITSClusters()<< endl;
		diff[4]=1;
	}
*/

// compare clusters
	if( verbose &&  track1->GetNumberOfTPCClusters() == track2->GetNumberOfTPCClusters()){
		for (Int_t idxCluster = 0; idxCluster < track1->GetNumberOfTPCClusters(); ++idxCluster){
			AliFlatTPCCluster * cl1 = track1->GetTPCCluster(idxCluster);
			AliFlatTPCCluster * cl2 = track2->GetTPCCluster(idxCluster);
/*
			if( cl1->GetX()&& cl2->GetX() && cl1->GetX() != cl2->GetX() ){
				cout<<"DIFFERENCE!: ";
			 	cout<<" cluster: "<<idxCluster<<" GetX :"<<cl1->GetX()<< " | " <<cl2->GetX()<< endl;
				diff=kTRUE;
			}
			 	cout<<" cluster: "<<idxCluster<<" GetX :"<<cl1->GetX()<< " | " <<cl2->GetX()<< endl;
			 	cout<<" cluster: "<<idxCluster<<" GetY :"<<cl1->GetY()<< " | " <<cl2->GetY()<< endl;

			if( cl1 && cl2 && cl1->GetY() != cl2->GetY() ){
				cout<<"DIFFERENCE!: ";
			 	cout<<" cluster: "<<idxCluster<<" GetY :"<<cl1->GetY()<< " | " <<cl2->GetY()<< endl;
				diff=kTRUE;
			}
			if( cl1->GetZ()&& cl2->GetZ() && cl1->GetZ() != cl2->GetZ() ){
				cout<<"DIFFERENCE!: ";
			 	cout<<" cluster: "<<idxCluster<<" GetZ :"<<cl1->GetZ()<< " | " <<cl2->GetZ()<< endl;
				diff=kTRUE;
			}
*/
			if( cl1->GetPadRow()&& cl2->GetPadRow() && cl1->GetPadRow() != cl2->GetPadRow() ){
				cout<<"DIFFERENCE!: ";
			 	cout<<" cluster: "<<idxCluster<<" GetPadRow :"<<cl1->GetPadRow()<< " | " <<cl2->GetPadRow()<< endl;
				diff[5]=1;
			}

		}
	  }
     
      track1 = track1->GetNextTrack();
      track2 = track2->GetNextTrack();
	  
	  
	  }
#endif
	   hStat->Fill(0);	  
	  Bool_t diffs=kFALSE;
	  for(int iDiff=0; iDiff<5;++iDiff){
		if(diff[iDiff]){
			hStat->Fill(iDiff+2);
			diffs = kTRUE;
		}
	}
	if(!diffs) hStat->Fill(1);	  


      curr1=curr1+ flatEsd1->GetSize();
      curr2=curr2+ flatEsd2->GetSize();
      iEvent++;
    }

    delete[] buffer1;
    delete[] buffer2;
  }
  else {
    cout << "File could not be read" << endl;
  }



  
	TList histosList;
	histosList.Add(hStat);
	histosList.Add(hNTracks);
	histosList.Add(hNV0s);
	histosList.Add(hVtxTr);
	histosList.Add(hVtxSPD);
  histosList.SaveAs(outputFilename);
  
  return;
}
