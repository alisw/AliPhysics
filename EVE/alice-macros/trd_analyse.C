void analyseClusters(TEvePointSet *points = 0x0)
{
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
  if(!line) {
    Info("analyseTracklet", "Invalid line.");
    return;
  }
  
  AliTRDseedV1 *tracklet = 0x0;
  tracklet = (AliTRDseedV1*)line->GetPointId(0);
  tracklet->Print();
}