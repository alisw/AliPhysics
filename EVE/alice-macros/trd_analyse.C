//
// How to access the basic TRD data structures when a 
// pointer to an Eve opject is available. One can get a 
// valid pointer from the event display by accessing the 
// function ExportToCINT 
// 
// Usage:
// .L trd_analysis.C
// analyseXXX(ptr)
// 
// Author:
// Alex Bercuci (A.Bercuci@gsi.de)
// 
void analyseDigits(AliEveTRDDigits *digits = 0x0)
{
// Simple print digits from a detector

  if(!digits) {
    Info("analyseDigits", "Invalid digits set.");
    return;
  }

  Int_t adc ;
  AliTRDdataArrayI *data = digits->GetData();
  Int_t nrows = data->GetNrow(),
        ncols = data->GetNcol(),
        ntbs  = data->GetNtime();
  data->Expand();
  printf("nrows[%d] ncols[%d] ntbs[%d]\n", nrows, ncols, ntbs);
  for (Int_t  row = 0;  row <  nrows;  row++)
  for (Int_t  col = 0;  col <  ncols;  col++)
    for (Int_t time = 0; time < ntbs; time++) {
      if((adc = data->GetDataUnchecked(row, col, time)) <= 1) continue;
      printf("r[%d] c[%d] t[%d] ADC[%d]\n", row, col, time, adc);
    }
  data->Compress(1);
}


void analyseClusters(TEvePointSet *points = 0x0)
{
// print some info about clusters in one detector or 
// attached to tracks.
  if(!points) {
    Info("analyseClusters", "Invalid points set.");
    return;
  }

  AliTRDcluster *c = 0x0;
  for(Int_t ic=0; ic<points->Size(); ic++){
    if(!(c = (AliTRDcluster*)points->GetPointId(ic))) continue;
  
    printf("%2d[%p] Det[%d] LabelMC[%d] TB[%d]\n", ic, c, c->GetDetector(), c->GetLabel(0), c->GetLocalTimeBin());
  }
}

void analyseTracklet(TEveLine *line)
{
// print tracklet information
  if(!line) {
    Info("analyseTracklet", "Invalid line.");
    return;
  }
  
  AliTRDseedV1 *tracklet = 0x0;
  tracklet = (AliTRDseedV1*)line->GetUserData(0);
  tracklet->Print();
}

