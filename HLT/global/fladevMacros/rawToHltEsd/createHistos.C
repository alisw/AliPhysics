void createHistos(){
  
	TString outputFilename = "$PWD/histosBenchmark.root";
	
	static const Int_t n=8;
	
	TString axes[n] = {"nEvents", 	"inputSize(kB)", 	"outputSize(kB)", 	"realTime(ms)", 	"cpuTime(ms)", 	"nTracks",	"timestamp (s)", 	"nV0s" };
	Int_t bins[n] 	= {100, 		200, 				200, 				2000,				1000,			200,		1600,				10};
	Double_t mins[n]= {0., 			0., 				0., 				0.,					0.,				0.,			0.,					0.};
	Double_t maxs[n]= {100., 		200., 				200., 				20.,				10.,			200.,		16000,				10.};
	
	THnSparseD * s = new THnSparseD("benchmarkInformation", "Benchmark information", n, bins, mins, maxs);
	
	for(int i=0; i<n;++i){
		s->GetAxis(i)->SetName(axes[i]);
		s->GetAxis(i)->SetTitle(axes[i]);
	}
	
	
	TTimeStamp ts;
	Int_t t = (Int_t) ts.GetSec() ;
	TNamed *time = new TNamed("time",Form("%d",t));
	
	
	TList histosList;
	histosList.Add(s);
	histosList.Add(time);
	
 // TH2F* hCpuTimeVsSize;
 // TH2F* hRealTimeVsSize;
 // TH2F* hCpuTimeVsTracks;
 // TH2F* hRealTimeVsTracks;
  
  
//	hCpuTimeVsSize = new TH2F("cpuTimeVsSize","cpu time vs. size", 1000,0,200, 1100,0,11);
//	hRealTimeVsSize = new TH2F("realTimeVsSize","real time vs. size", 1000,0,200, 1000,0,10);
//	hCpuTimeVsTracks = new TH2F("cpuTimeVsTracks","cpu time vs. number of tracks", 200,0,200, 1100,0,11);
//	hRealTimeVsTracks = new TH2F("realTimeVsTracks","real time vs. number of tracks", 200,0,200, 1000,0,10);

//	histosList.Add(hCpuTimeVsSize);
//	histosList.Add(hRealTimeVsSize);
//	histosList.Add(hCpuTimeVsTracks);
//	histosList.Add(hRealTimeVsTracks);
	histosList.SaveAs(outputFilename);
	
	
	
	
	
}
