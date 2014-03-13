void makeEventList(const char* file, Double_t ptMinHighPt = 8., Double_t ptMinV0s = 3.)
{
  ///////////////////
  // make event lists based on the filtered ESD trees (Filter_Events.root)
  // WARNING: output of this needs to be parsed by the makeEventList.sh script!

  TFile f(file);

  TTree* c = NULL;
  Long_t nEvents=0;

  c=(TTree*)f.Get("highPt");
  if (c)
  {
    if (c->GetEntries()>0)
    { 
      printf("offlineTrigger: highPt\n");
      c->SetScanField(nEvents);
      c->Scan("esdTrack.Pt():runNumber:evtNumberInFile:fileName.GetString():gid:evtTimeStamp",Form("esdTrack.Pt()>%lf",ptMinHighPt),"col=.2f:8.d:8.d:130.s:15.lu:12.d");
    }
  }

  c=(TTree*)f.Get("V0s");
  if (c)
  {
    if (c->GetEntries()>0)
    { 
      printf("offlineTrigger: V0s\n");
      c->SetScanField(nEvents);
      c->Scan("v0.Pt():runNumber:evtNumberInFile:fileName.GetString():gid:evtTimeStamp",Form("v0.Pt()>%lf",ptMinV0s),"col=.2f:8.d:8.d:130.s:15.lu:12.d");
    }
  }

  c=(TTree*)f.Get("Laser");
  if (c)
  {
    if (c->GetEntries()>0)
    { 
      printf("offlineTrigger: Laser\n");
      c->SetScanField(nEvents);
      c->Scan("runNumber:runNumber:evtNumberInFile:fileName.GetString():gid:evtTimeStamp","","col=8.d:8.d:8.d:130.s:15.lu:12.d");
    }
  }

  c=(TTree*)f.Get("CosmicPairs");
  if (c)
  {
    if (c->GetEntries()>0)
    { 
      printf("offlineTrigger: CosmicPairs\n");
      TCut ptCut="abs(t0.fP[4])<0.33"; //cut on 1/pt < 0.33
      TCut cutDCA="abs(0.5*(t0.fD-t1.fD))>5&&abs(0.5*(t0.fD-t1.fD))<80"; //tracks crossing the inner field cage (80cm)
      TCut cutCross="t0.fOp.fP[1]*t1.fOp.fP[1]<0"; //tracks crossing central electrode
      c->SetScanField(nEvents);
      c->Scan("runNumber:runNumber:evtNumberInFile:fileName.GetString():gid:evtTimeStamp", ptCut && cutDCA && cutCross,"col=8.d:8.d:8.d:130.s:15.lu:12.d");
    }
  }
}//dumpList

