/// \file AliTPCMonitorExec.C

void AliTPCMonitorExec(Int_t val)
{
  AliTPCMonitor* moni  =0;
  moni  = (AliTPCMonitor*)gROOT->FindObjectAny("monitor");
  
  if(moni=(AliTPCMonitor*)gROOT->FindObjectAny("monitor"))
    { 
      if      (val==1){   moni->ExecPad()      ;}
      else if (val==2){   moni->ExecRow()      ;}
      else if (val==3){   moni->ExecProcess()       ;}
      else     cout << " no value specified " << endl;
    }
  else
    {
      cout << " Object not found " << endl;
    }
       
}
