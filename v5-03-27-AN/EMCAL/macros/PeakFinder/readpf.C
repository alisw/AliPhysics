
void readpf()
{
  cout << "HELLO WORLD" << endl;
  TFile *f = new TFile("peakfindervectors2.root", "read");
  
  //   f->ls();
  //   f->Print();
  //   AliCaloRawAnalyzerPeakFinder *p = (AliCaloRawAnalyzerPeakFinder*)f->GetKey("AliCaloRawAnalyzerPeakFinder");
  
  AliCaloPeakFinderVectors *p = f->Get("AliCaloPeakFinderVectors");

  if ( p == 0 )
    {
      cout << "ERROR, P == 0" << endl;
    }
  else
    {
      cout << "INFO, p = "<< p  << endl;
      p->PrintVectors();
    }

  cout << endl <<  "********** !!!!!!!!!!!!!! " << endl;
  
  
  // AliCaloRawAnalyzerPeakFinder *p2 =  new AliCaloRawAnalyzerPeakFinder();
  //  p2->PrintVectors();
}


