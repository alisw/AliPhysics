int total=0; int rejected=0;
Float_t fromPt=1.; Float_t toPt=5.;

namespace RawProduction {
  class Output;
}

void CountReject(RawProduction::Output& output, const char* trigger, int cent, const char* pid, const char* methode)
{
  TH1* hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, pid, methode));
  int fromBin = hist->FindBin(fromPt);
  int toBin = hist->FindBin( toPt);
  for(int bin = fromBin; bin<=toBin; ++bin) {
    ++total;
    if( 0. == hist->GetBinContent(bin))
      ++rejected;
  }
}


void CountRejected()
{
  gROOT->LoadMacro("MakeRawProduction.C+g");
  RawProduction::Output output("RawProduction.root");
  gStyle->SetOptStat(0);
  
  TStringToken pids("All Allcore Allwou Disp Disp2 Dispcore Dispwou CPV CPVcore CPV2 Both Bothcore", " ");
  while( pids.NextToken() ) {
    TStringToken fits("mr1r mr1 mr2r mr2", " ");
    while( fits.NextToken() ) {
      CountReject(output, "kMB", -1, pids.Data(), fits.Data());
      //CountReject(output, "kMB", -11, pids.Data(), fits.Data());
      CountReject(output, "kMB", -10, pids.Data(), fits.Data());
      //CountReject(output, "kMB", -6, pids.Data(), fits.Data());
      
      CountReject(output, "kCentral", -1, pids.Data(), fits.Data());

      //CountReject(output, "kSemiCentral", -11, pids.Data(), fits.Data());
    }
  }
    

  Printf("fits:%i,  rejected:%i", total, rejected);
}
