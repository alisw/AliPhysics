#include <iomanip>
void
MakeFakeHits()
{
  if (!gAlice) { 
    std::cerr << "This script must be run in AliROOT" <<std::edl;
    return;
  }
  
  // Initialize 
  gAlice->InitMC("$ALICE_ROOT/FMD/Config.C");

  // Get our detector 
  AliFMD* fmd = gAlice->GetDetector("FMD");
  if (!fmd) { 
    std::cerr << "FMD object not defined" << std::endl;
    return;
  }
  // Get the runloader
  AliRunLoader* runLoader = gAlice->GetRunLoader();
  if (!runLoader) { 
    std::cerr << "No run loader defined" << std::end;
    // return;
  }

  // OCDB manager
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB")
  cdb->SetRunNumber(0)

  // Geometry database 
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();

  // Parameter database 
  // AliFMDParameters* param = AliFMDParameters::Instance();
  // param->Init(kFALSE, AliFMDParameters::kAltroMap);
  AliFMDAltroMapping map;

  // Monte-carlo application
  AliMC* mc = gAlice->GetMCApp();
  if (!mc) { 
    std::cerr << "No MC application defined" << std::endl;
    return;
  }

  // Make primaries - one for each strip
  mc->BeginEvent();
  Int_t ntr = 0;
  for (size_t i = 0; i < 51200; i++) {
    mc->PushTrack(1, -1, 211, 
		  0, 0, 0, 0, 
		  0, 0, 0, 0, 
		  0, 0, 0, 
		  kPPrimary, ntr);
    // std::cout << "Made track # " << ntr << std::endl;
  }

  // Make hits
  ntr = 0;
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nrng = (d == 1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nrng; ir++) { 
      Char_t   r    = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ?  20 :  40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      for (UShort_t s = 0; s < nsec; s++) { 
	for (UShort_t t = 0; t < nstr; t++) { 
	  mc->BeginPrimary();
	  
	  Double_t x, y, z;
	  geom->Detector2XYZ(d, r, s, t, x, y, z);
	  UInt_t ddl, board, altro, channel;
	  UShort_t timebin;
	  map.Detector2Hardware(d, r, s, t, 0, 0, 1, 
				ddl, board, altro, channel, timebin);
	  Float_t e = Float_t(s) / nsec + Float_t(t)/(100*nstr);
	  e = Float_t(timebin % 64) / 64 * 1.5;

	  std::cout << "FMD" << d << r << "[" << std::setfill('0') 
		    << std::setw(2) << s << "," << std::setw(3) << t 
		    << "] " << std::setfill(' ')
		    << std::setw(4) << timebin << " ("
		    << std::setw(10) << x << "," 
		    << std::setw(10) << y << "," 
		    << std::setw(10) << z << ") -> " 
		    << std::setw(10) << e << "\t\r"
		    << std::flush;
	  fmd->AddHitByFields(ntr,     // Int_t    track, 		
			      d,       // UShort_t detector, 	
			      r,       // Char_t   ring, 		
			      s,       // UShort_t sector, 	
			      t,       // UShort_t strip, 		
			      x,       // Float_t  x=0,		
			      y,       // Float_t  y=0, 		
			      z,       // Float_t  z=0,		
			      0,       // Float_t  px=0, 		
			      0,       // Float_t  py=0, 		
			      0,       // Float_t  pz=0,		
			      e,       // Float_t  edep=0,		
			      211,     // Int_t    pdg=0,		
			      0,       // Float_t  t=0, 			
			      0.03);   // Float_t  len=0, 		
	  mc->FinishPrimary();	       
	  ntr++;		      
	} // Strip
      } // Sector
      std::cout << std::endl;
    } // Ring 
  } // Detector

  // End of event.
  mc->FinishEvent();
  
  // End of run
  gAlice->FinishRun();
  gGeoManager->Export("geometry.root");
  
}

  
    
      
	  
