// AliTRDcompareDigits.C
//
// Author
//   Ken Oyama
//   Mateusz Ploskon
//
// Example and diagnostics macro to compare digits
// before and after the raw data simulation.
// Please simulate raw data in root format with keeping
// all intermediate data unerased to run the macro.

class Histograms {
public:
  Histograms() {
    adc_all = new TH1F("adc_all", "adc_all", 1025, -1.5, 1023.5 );
  }

  ~Histograms() {
    delete adc_all;
  }

  setattr( Color_t c, Color_t fc, Float_t lw ) {
    if( c  != -1 ) adc_all->SetLineColor( c );
    if( fc != -1 ) adc_all->SetFillColor( fc );
    if( lw != -1 ) adc_all->SetLineWidth( lw );
  }

  TH1F *adc_all; 
};

AliTRDgeometry *gGeo;

void AliTRDcompareDigits( Int_t id_start = 0 , Int_t id_end = 539 )
{
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  //AliCDBStorage* localStorage = cdb->GetStorage("local://$ALICE_ROOT/OCDB");
  //  cout << "[I] Is storage set : " << cdb->IsDefaultStorageSet() << endl;
  cout << endl;

  // Read raw data after zero suppression
  AliLog::SetClassDebugLevel("AliTRDRawStream", 10);
  AliLog::SetFileOutput("decodeZS.log");
  AliRawReaderRoot reader("raw.root", 0);
  reader.Select("TRD");
  AliTRDrawData data;
  AliTRDdigitsManager *dmz = data.Raw2Digits(&reader);

  gGeo = new AliTRDgeometry(); 
  Ho   = new Histograms(); // Original
  Hz   = new Histograms(); // Zero suppressed
  Ho->setattr( 46 , 46, -1);
  Hz->setattr( -1 , -1, 2 );

  // Read original data from digits tree
  AliRunLoader *runLoader  = AliRunLoader::Open("galice.root" ,AliConfig::GetDefaultEventFolderName(),"READ");
  AliLoader *loader = runLoader->GetLoader("TRDLoader");
  loader->LoadDigits();
  AliTRDdigitsManager *dmo = new AliTRDdigitsManager();  // Original digits
  dmo->ReadDigits(loader->TreeD());

  for( int id = id_start ; id <= id_end ; id++ ) {
    AliTRDdataArrayI *dao = dmo->GetDigits( id );
    AliTRDdataArrayI *daz = dmz->GetDigits( id );
    dao->Expand();
    daz->Expand();
    analyze_det( id, dao, Ho );
    analyze_det( id, daz, Hz );
    compare_signal_position( id, dao, daz );
  }

  TCanvas *c = new TCanvas("decodeZS", "decodeZS", 800, 500);
  c->Draw();
  c->SetLogy();
  Ho->adc_all->Draw();
  Hz->adc_all->Draw("SAME");
}

//
// Analyze detector digits, printout data and fill histogram
//
void analyze_det( Int_t det, AliTRDdataArrayI *d, Histograms *h )
{
  Int_t        plan = gGeo->GetPlane( det );   // Plane
  Int_t        cham = gGeo->GetChamber( det ); // Chamber
  Int_t        sect = gGeo->GetSector( det );  // Sector (=iDDL)
  Int_t        nRow = gGeo->GetRowMax( plan, cham, sect );
  Int_t        nCol = gGeo->GetColMax( plan );
  Int_t       nTBin = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
  Int_t       ctype = 0;                       // Chamber type (0:C0, 1:C1)
  Int_t        adc;

  printf("Det=%03d Sector=%02d Stack=%d Layer=%d\n", det,sect, cham, plan );
  // printf("Size=%d\n", d->GetSize());

  for (Int_t irow = 0; irow < nRow; irow++ ) {
    for (Int_t icol = 0; icol < nCol; icol++ ) {
      // printf(" Row=%02d Col=%03d :", irow, icol);
      Int_t signal = 0;
      for (Int_t it = 0; it < nTBin; it++ ) {
 	adc = d->GetData(irow, icol, it);
	// printf(" % 3d", adc );
	h->adc_all->Fill( adc );
	if( adc > 5 ) signal++;
      }
      // if( signal != 0 ) printf(" ... signal!");
      // printf("\n");
    }
  } 
}

//
// Compare two digits array
//
void compare_signal_position( Int_t det, AliTRDdataArrayI *dao, AliTRDdataArrayI *daz )
{
  Int_t        plan = gGeo->GetPlane( det );   // Plane
  Int_t        cham = gGeo->GetChamber( det ); // Chamber
  Int_t        sect = gGeo->GetSector( det );  // Sector (=iDDL)
  Int_t        nRow = gGeo->GetRowMax( plan, cham, sect );
  Int_t        nCol = gGeo->GetColMax( plan );
  Int_t       nTBin = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
  Int_t       wrong = 0;
  Int_t       adc_o;
  Int_t       adc_z;

  for (Int_t irow = 0; irow < nRow; irow++ ) {
    for (Int_t icol = 0; icol < nCol; icol++ ) {
      for (Int_t it = 0; it < nTBin; it++ ) {
	adc_o = dao->GetData(irow, icol, it);
	adc_z = daz->GetData(irow, icol, it);
	if( adc_o > 5 && adc_o != adc_z ) {
	  printf("              can not find adc=%d at Row=%02d Col=%03d TimeBin=%02d in ZS digits.\n", adc_o, irow, icol, it ); 
	  wrong++;
	}
      }
    }
  } 
  if( wrong != 0 ) printf("               %d wrong ADC hit(s) found !!!\n", wrong);
}
