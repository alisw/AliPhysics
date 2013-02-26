#include <fstream>
#include "TString.h"
#include "TDatime.h"
#include "TPRegexp.h"
#include "TArrayD.h"
#include "TMath.h"
#include <TCanvas.h>
#include <TH1F.h>

const int runCapacity = 500;
const int nMods =3;
const int sensorCapacity = 3;

const int nStat = 3;
enum StatTypes {kEntries=0, kMean=1, kRMS=2};

const char* csvTemperatureFileNames[3] = {"m4_matrix_temp_2013-02-19.csv", "m3_matrix_temp_2013-02-19.csv", "m2_matrix_temp_2013-02-19.csv"};
const char* runListFileName = "lhc13b.csv";
const char* periodFromDatime = "2013-01-20 12:00:00";
const char* periodToDatime = "2013-02-10 12:00:00";

TH1* hists[nMods][sensorCapacity] = {{0}};

UInt_t nRuns =0;
UInt_t runIndexes[runCapacity] = {0};
UInt_t fromTimes[runCapacity] = {0};
UInt_t toTimes[runCapacity] = {0};

std::ifstream* csvStreams[nMods] = {0};
TString* sensorNames[nMods][sensorCapacity] = {{0}};


TArrayD* temps[nMods][sensorCapacity][runCapacity] = {{{0}}};


Double_t stats[nMods][sensorCapacity][nStat][runCapacity] = {{{{0.}}}};
Double_t stat_errs[nMods][sensorCapacity][nStat][runCapacity] = {{{{0.}}}};




void ReadRuns(){
  Printf("reading runs");

  std::ifstream ifs(runListFileName);
  char cline[1024] = "";
  ifs.getline(cline, 1024);
  while(ifs.good()) {
    ifs.getline(cline, 1024);
    TString tline(cline);
    int index = tline.Index(",");
    if( index == -1) continue; // last line

    // date format fix:
    //Printf(tline.Data());
    tline[7] = cline[13];
    tline[8] = cline[14];
    tline[9] = cline[15];
    tline[10] = cline[16];
    tline[11] = '-';
    tline[12] = cline[10];
    tline[13] = cline[11];
    tline[14] = '-';
    tline[15] = cline[7];
    tline[16] = cline[8];
    
    tline[27] = cline[33];
    tline[28] = cline[34];
    tline[29] = cline[35];
    tline[30] = cline[36];
    tline[31] = '-';
    tline[32] = cline[30];
    tline[33] = cline[31];
    tline[34] = '-';
    tline[35] = cline[27];
    tline[36] = cline[28];
    
    
    //Printf(tline.Data());

    TStringToken fields(tline.Data(), ",");
    fields.NextToken();
    runIndexes[nRuns] = fields.Atoi();
    fields.NextToken();
    printf("%s ", fields.Data());
    TDatime fromTime(fields);
    fromTimes[nRuns] = fromTime.Convert();
    fields.NextToken();
    Printf(fields.Data());
    TDatime toTime(fields);
    toTimes[nRuns] = toTime.Convert();
    
    //Printf("run #%i, from=%s, to=%s", runIndexes[nRuns], fromTime.AsString(), toTime.AsString());
    ++nRuns;
  }
  
  Printf("done reading runs");
}

void ReadCSVHeaders() {
  // Should be called prior to ReadCSVTemp, opens ifstreams and read the first lines.
  Printf("\nOpening files and reading headers");
  for(int module = 0; module < nMods; ++module) {
    printf("module %i: ", module);
    csvStreams[module] = new std::ifstream(csvTemperatureFileNames[module]);
    std::ifstream* ifs = csvStreams[module];
    
    char cline[1024] = "";
    ifs->getline(cline, 1024); // Fist line is DCS details
    
    ifs->getline(cline, 1024); // Second line is sensor names
    TStringToken fnames(cline, ",", true);
    fnames.NextToken(); // first token will be emtpy
    
    int sensor = -1;
    while(++sensor < sensorCapacity && fnames.NextToken() ) {
      if(1==module && sensor>1 ) continue;
      sensorNames[module][sensor] = new TString(fnames);
      printf("%s, ", sensorNames[module][sensor]->Data());
    }
    printf("\n");
  }
  Printf("done reading CSV headers");
}

void MakePeriodHistograms()
{
  for(int mod = 0; mod < nMods; ++mod)
    for(int sens = 0; sens < sensorCapacity; ++sens)
      if( sensorNames[mod][sens] ) {
	hists[mod][sens] = new TH1I(sensorNames[mod][sens]->Data(), sensorNames[mod][sens]->Data(), 1000, -50, 0 );
      }
}

void ReadCSVTemp() {
  Printf(" reading data ");

  UInt_t fromTime = TDatime(periodFromDatime).Convert();
  UInt_t toTime = TDatime(periodToDatime).Convert();

  // Make Arrays
  for(int mod = 0; mod < nMods; ++mod)
    for(int sens = 0; sens < sensorCapacity; ++sens)
      for(int run = 0; run < nRuns; ++ run)
	temps[mod][sens][run] = new TArrayD(10);
    
  // go through lines
  char cline[1024] = "";
  for(int mod=0; mod<nMods; ++mod){
    Printf("module %i", mod);
    while( csvStreams[mod]->good() ) {
      csvStreams[mod]->getline(cline, 1024);
      TString tline(cline);
      int index = tline.Index(",");
      if( index == -1) continue; // last line
	
      // Date/Time
      TString dtline = tline(0,19);
      dtline.ReplaceAll("/", "-");
      TDatime datime = TDatime( dtline.Data() );
      UInt_t dateTime = datime.Convert();


      for(UInt_t run=0; run<nRuns; ++run){
      	// data
	TStringToken fields(tline, ",",true);
	fields.NextToken(); // first token is date
	int sensor = -1;
	while(++sensor < sensorCapacity && fields.NextToken() && fields.IsFloat() && sensorNames[mod][sensor]) {
	  if(fromTime < dateTime && dateTime < toTime ) {
	    hists[mod][sensor]->Fill(fields.Atof());
	    
	    if(sensorNames[mod][sensor]->EqualTo("Matrix_temperature_module_2_18"))
	      Printf("%s; %f", sensorNames[mod][sensor]->Data(), fields.Atof());
	  }


	  if( fromTimes[run] < dateTime && dateTime < toTimes[run] ) {
	    TArrayD* array = temps[mod][sensor][run];
	    Double_t& entries = stats[mod][sensor][kEntries][run];
	    ++entries;
	    if( entries > array->GetSize() )
	      array->Set(entries *2);
	    array->operator[](entries-1) = fields.Atof();
	  }
	} // sensor
      } // run
    } // line
  } // module
  
  // Extract Stats
  for(int mod = 0; mod < nMods; ++mod)
    for(int sens = 0; sens < sensorCapacity; ++sens)
      for(int run = 0; run < nRuns; ++ run){
	Double_t entries = stats[mod][sens][kEntries][run];
	//if ( entries < 3) continue;
	Double_t* array = temps[mod][sens][run]->GetArray();
	Double_t mean = TMath::Mean(entries, array);
	Double_t rms = TMath::RMS(entries, array);
	stats[mod][sens][kMean][run] = mean;
	stats[mod][sens][kRMS][run] = rms;
	//stat_errs[mod][sens][kMean][run] = rms;
      }
}

void PlotTemps() {

  bool first = true;
  
  for(int mod = 0; mod < nMods; ++mod)
    for(int sens = 0; sens < sensorCapacity; ++sens){
      if( ! sensorNames[mod][sens] ) continue;
      TString name = Form("hEntries_%s", sensorNames[mod][sens]->Data());
      TString title = Form("hEntries_%s", sensorNames[mod][sens]->Data());
      TH1F* hEntries = new TH1F(name.Data(), title.Data(), nRuns, 0, nRuns);
      for(int run = 0; run < nRuns; ++ run){
	hEntries->SetBinContent(run+1, stats[mod][sens][kEntries][run]);
	//hEntries->SetBinError(run+1, stat_errs[mod][sens][kEntries][run]);
	Printf("%f", stats[mod][sens][kEntries][run] );
      }
      //hEntries->GetYaxis()->SetRangeUser(-30, -20);
      if(first)
	hEntries->Draw();
      else
	hEntries->Draw("same");
      first = false;
    }

 return;  
}

void PlotPeriodHistograms()
{
  for(int mod = 0; mod < nMods; ++mod)
    for(int sens = 0; sens < sensorCapacity; ++sens)
      if( sensorNames[mod][sens] ) {
	TCanvas* canv = new TCanvas(sensorNames[mod][sens]->Data(), sensorNames[mod][sens]->Data());
	hists[mod][sens]->Draw();
	canv->Print(Form("imgs/%s.pdf", sensorNames[mod][sens]->Data()));
	canv->Print(Form("imgs/%s.png", sensorNames[mod][sens]->Data()));
      }
}

void temp(){
  ReadRuns();
  ReadCSVHeaders();
  MakePeriodHistograms();
  ReadCSVTemp();
  //PlotTemps();
  PlotPeriodHistograms();
}
