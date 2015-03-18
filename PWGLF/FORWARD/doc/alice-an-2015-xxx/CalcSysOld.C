TObjString* Token(TObjArray* a, Int_t i)
{
  return static_cast<TObjString*>(a->At(i));
}
void GetTwo(const TString& cell, Float_t& low, Float_t& high)
{
  TObjArray* parts = cell.Tokenize(" ");
  low  = Token(parts, 0)->String().Atof();
  high = Token(parts, 2)->String().Atof();
  parts->Delete();
  delete parts;
}

TNtuple* ReadOld()
{
  TFile*   out = TFile::Open("old.root","RECREATE");
  TNtuple* ret = new TNtuple("old", "old", 
			     "etaMin:etaMax:v0005:e0005:"
			     "v0510:e0510:v1020:e1020:v2030:e2030");
  std::ifstream in("old.dat");
  do { 
    TString line;
    line.ReadLine(in);
    TObjArray*  columns = line.Tokenize(";");
    TIter       next(columns);
    TObjString* cell = 0;
    Float_t     x[]  = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    Float_t*    px   = x;
    while ((cell = static_cast<TObjString*>(next()))) {
      GetTwo(cell->String(), *px, *(px+1));
      px += 2;
      
    }
    columns->Delete();
    delete columns;

    // std::cout << line << "\t"; 
    // for (Int_t i = 0; i < 10; i++) std::cout << x[i] << "\t";
    // std::cout << std::endl;

    if (x[9] == 0) break;

    ret->Fill(x);

  } while (!in.eof());
  
  in.close();
  out->Write(); 
  out->Close();
  
  out = TFile::Open("old.root","READ");
  ret = static_cast<TNtuple*>(out->Get("old"));

  return ret;
}

  

void CalcSysOld()
{
  TNtuple* nt  = ReadOld();

  nt->Draw("v0005:(etaMax+etaMin)/2");


  for (Int_t i = 0; i < nt->GetEntries(); i++) { 
    Float_t*  cells = nt->GetArgs();
    Double_t  eta   = (cells[1]+cells[0])/2;
    Double_t  dEta  = (cells[1]-eta);
    
    printf("%6f +/- %6f", eta, dEta);
    Double_t sumSq = 0;
    for (Int_t j = 0; j < 4; j++) {
      Int_t    k    = 2*(j+1);
      Double_t v    = cells[k+0];
      Double_t e    = cells[k+1];
      Double_t r    = e/v * 100;
      printf(" %7f", r);
      sumSq         += r*r;
    }
    Printf(" -> %f", TMath::Sqrt(sumSq));
  }
  out->Write();
}
