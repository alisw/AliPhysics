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
#include "THnSparse.h"
#endif   
Double_t printDiff(string name, double val1, double val2);
Double_t printDiff(string name, TString val1, TString val2);
void CompareFlatESDs(const char* filename1="outFlatESD1.dat",const char* filename2="outFlatESD2.dat", Bool_t verbose=kFALSE) {
  // Create output histograms
  
  
  TString outputFilename = "$PWD/compare.root";
  
	cout<< "creating histograms"<<endl;
	THnSparse * hDiff;
	const	Int_t nDim = 12;
	
		Int_t bins[nDim] = {2};
		Double_t mins[nDim] = {0};
		Double_t maxs[nDim] = {0};
		hDiff = new THnSparseD("Differences","Differences",nDim,bins,mins,maxs);
	
  

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

      AliFlatESDEvent *flatEsd1 = reinterpret_cast<AliFlatESDEvent *>(curr1);
      AliFlatESDEvent *flatEsd2 = reinterpret_cast<AliFlatESDEvent *>(curr2);
	  
	  flatEsd1->Reinitialize();
	  flatEsd2->Reinitialize();
	  
      cout<<endl<<"________________________________________________________________________"<<endl;
      cout<<endl<<"Reading event "<<iEvent<<":\t file1 | file 2\t|\t abs. diff\t|\t rel. diff"<<endl;
		
		Double_t diffs[nDim] ={
			printDiff("GetMagneticField",		flatEsd1->GetMagneticField(),		flatEsd2->GetMagneticField()),
			printDiff("GetPeriodNumber",		flatEsd1->GetPeriodNumber(),		flatEsd2->GetPeriodNumber()),
			printDiff("GetRunNumber",				flatEsd1->GetRunNumber(),				flatEsd2->GetRunNumber()),
			printDiff("GetOrbitNumber",			flatEsd1->GetOrbitNumber(),			flatEsd2->GetOrbitNumber()),
			printDiff("GetBunchCrossNumber",flatEsd1->GetBunchCrossNumber(),flatEsd2->GetBunchCrossNumber()),
			printDiff("GetTriggerMask",			flatEsd1->GetTriggerMask(),			flatEsd2->GetTriggerMask()),
			printDiff("GetTriggerMaskNext50",flatEsd1->GetTriggerMaskNext50(),			flatEsd2->GetTriggerMaskNext50()),
			printDiff("GetFiredTriggerClasses",			flatEsd1->GetFiredTriggerClasses(),			flatEsd2->GetFiredTriggerClasses()),
			printDiff("GetNumberOfTracks",	flatEsd1->GetNumberOfTracks(),	flatEsd2->GetNumberOfTracks()),
			printDiff("GetNumberOfV0s",			flatEsd1->GetNumberOfV0s(),			flatEsd2->GetNumberOfV0s()),
			printDiff("GetTimeStamp",			flatEsd1->GetTimeStamp(),			flatEsd2->GetTimeStamp()),
			printDiff("GetEventSpecie",			flatEsd1->GetEventSpecie(),			flatEsd2->GetEventSpecie())
	
		};

		hDiff->Fill(diffs);
	  
	  
	  /*
	  
	  if( (Bool_t) flatEsd1->GetFlatPrimaryVertexTracks() != (Bool_t) flatEsd2->GetFlatPrimaryVertexTracks()  ){
		cout<<"\t\tDIFFERENCE!: "<<endl;
		diff[2] =1;
	}


	  cout<<"vtx tracks:\t"<<(Bool_t) flatEsd1->GetFlatPrimaryVertexTracks()<< " | " <<	(Bool_t) flatEsd2->GetFlatPrimaryVertexTracks()<<endl;	  
	  //hVtxTr->Fill( (Bool_t) flatEsd1->GetFlatPrimaryVertexTracks(), (Bool_t) flatEsd2->GetFlatPrimaryVertexTracks());

	 
	  
	  if( (Bool_t) flatEsd1->GetFlatPrimaryVertexSPD() != (Bool_t) flatEsd2->GetFlatPrimaryVertexSPD()  ){
		cout<<"\t\tDIFFERENCE!: "<<endl;
		diff[3] =1;
	}
      cout<<"vtx SPD:\t"<<(Bool_t) flatEsd1->GetFlatPrimaryVertexSPD() << " | " << (Bool_t) flatEsd2->GetFlatPrimaryVertexSPD()<<endl;
	  //hVtxSPD->Fill( (Bool_t) flatEsd1->GetFlatPrimaryVertexSPD(), (Bool_t) flatEsd2->GetFlatPrimaryVertexSPD());

	  
	  
  if((Bool_t)flatEsd1->GetFlatPrimaryVertexTracks() && (Bool_t)flatEsd2->GetFlatPrimaryVertexTracks()  ){
 		cout<<endl<<"vtx tracks -> X,Y,Z:\t"
			<< flatEsd1->GetFlatPrimaryVertexTracks()->GetX()
			<<","<< flatEsd1->GetFlatPrimaryVertexTracks()->GetY()
			<<","<< flatEsd1->GetFlatPrimaryVertexTracks()->GetZ()
			<<" | " <<flatEsd2->GetFlatPrimaryVertexTracks()->GetX()
			<<","<< flatEsd2->GetFlatPrimaryVertexTracks()->GetY()
			<<","<< flatEsd2->GetFlatPrimaryVertexTracks()->GetZ()<<endl;	  
		}
	  */
	  
	  // compare tracks
if(verbose){
	  AliFlatESDTrack *track1 = const_cast<AliFlatESDTrack*> (flatEsd1->GetTracks());
	  AliFlatESDTrack *track2 = const_cast<AliFlatESDTrack*> (flatEsd2->GetTracks());
    for (Int_t idxTrack = 0; idxTrack < flatEsd1->GetNumberOfTracks() && track1 && track2; ++idxTrack) { 

		//track2->Reinitialize();
		const AliFlatExternalTrackParam* ext[2][nExt] ={
			{
				track1->GetFlatTrackParamRefitted(),
				track1->GetFlatTrackParamIp(),
				track1->GetFlatTrackParamTPCInner(),
				track1->GetFlatTrackParamOp(),
		//		track1->GetFlatTrackParamCp(),
		//		track1->GetFlatTrackParamITSOut()
			},
			{
				track2->GetFlatTrackParamRefitted(),
				track2->GetFlatTrackParamIp(),
				track2->GetFlatTrackParamTPCInner(),
				track2->GetFlatTrackParamOp(),
			//	track2->GetFlatTrackParamCp(),
			//	track2->GetFlatTrackParamITSOut()
			}
		};
	
     	//Printf("  TEST: FlatTrack1 %d > FlatExternalTrackParam1 > %p %p %p %p", idxTrack, exp11, exp21, exp31, exp41);
     	//Printf("  TEST: FlatTrack2 %d > FlatExternalTrackParam2 > %p %p %p %p", idxTrack, exp12, exp22, exp32, exp42);


	for(int iExt=0; iExt<nExt; ++iExt){
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
			cout<<"\t\tDIFFERENCE!: "<<endl;
	 		//cout<<" alpha"<<iExt<<" :"  << (ext[0][iExt] ? ext[0][iExt]->GetAlpha() : -99.)  << "\t\t" << (ext[1][iExt] ?  ext[1][iExt]->GetAlpha(): -99.)<<endl;
		}	cout<<" alpha"<<iExt<<" :\t"  << (ext[0][iExt] ? ext[0][iExt]->GetAlpha() : -99.)  << " | " << (ext[1][iExt] ?  ext[1][iExt]->GetAlpha(): -99.)<<endl;
			

		if( (!ext[0][iExt] || !ext[1][iExt])||ext[0][iExt]->GetX() != ext[1][iExt]->GetX() ) {
			cout<<"\t\tDIFFERENCE!: "<<endl;
	 		//cout<<" GetX"<<iExt<<" :"  << (ext[0][iExt] ? ext[0][iExt]->GetX(): -99.)  << " | " << (ext[1][iExt] ?  ext[1][iExt]->GetX(): -99.)<<endl;
		}	
cout<<" GetX"<<iExt<<" :\t"  << (ext[0][iExt] ? ext[0][iExt]->GetX(): -99.)  << " | " << (ext[1][iExt] ?  ext[1][iExt]->GetX(): -99.)<<endl;


		if( (!ext[0][iExt] || !ext[1][iExt])||ext[0][iExt]->GetSigned1Pt() !=  ext[0][iExt]->GetSigned1Pt() ) {
			cout<<"\t\tDIFFERENCE!: "<<endl;
	 		//cout<<" 1/pt"<<iExt<<" :"  <<  (ext[0][iExt] ? ext[0][iExt]->GetSigned1Pt(): -99.)  << " | " << (ext[1][iExt] ?  ext[1][iExt]->GetSigned1Pt(): -99.)<<endl;
		}	
	cout<<" 1/pt"<<iExt<<" :\t"  <<  (ext[0][iExt] ? ext[0][iExt]->GetSigned1Pt(): -99.)  << " | " << (ext[1][iExt] ?  ext[1][iExt]->GetSigned1Pt(): -99.)<<endl;
			

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
	  
#if 0

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
	  
#endif
      track1 = const_cast<AliFlatESDTrack*> (track1->GetNextTrack());
      track2 = const_cast<AliFlatESDTrack*> (track2->GetNextTrack());
	  
	  
	  }
}
	  
	  /*
	//   hStat->Fill(0);	  
	  Bool_t diffs=kFALSE;
	  for(int iDiff=0; iDiff<5;++iDiff){
		if(diff[iDiff]){
	//		hStat->Fill(iDiff+2);
			diffs = kTRUE;
		}
	}
	if(!diffs) hStat->Fill(1);	  
*/

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
	histosList.Add(hDiff);
  histosList.SaveAs(outputFilename);

  return;
}

Double_t printDiff(string name, double val1, double val2){
	double relDiff = ( val1 != 0 || val2!=0 ) ? fabs(val1-val2)/(fabs(val1) + fabs(val2)): 0;
	cout<<name<<":\t"<<val1<<" | " << val2 <<"\t|\t"<<(val1-val2)<<"\t|\t"<<relDiff<<endl;
	return relDiff > 1e-6 ? 1:0;
}
Double_t printDiff(string name, TString val1, TString val2){
	cout<<name<<":"<<endl<<"\t"<<val1<<endl<<"\t"<< val2 <<endl;
	return val1.EqualTo(val2) ?0:1;
}
