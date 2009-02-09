#include <iomanip>

void
FullMapping()
{
  AliGeomManager::LoadGeometry("geometry.root");
  AliCDBManager*    cdb   = AliCDBManager::Instance();
  AliFMDParameters* param = AliFMDParameters::Instance();
  AliFMDGeometry*   geom  = AliFMDGeometry::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  param->Init();
  geom->Init();
  geom->InitTransformations();
  AliFMDAltroMapping* map = param->GetAltroMap();

  std::ostream& out = cerr;
  
  out << "  DDL | Board | Altro | Channel |"
      << " Detector | Ring | Sector | Base-Strip |  Phi |"
      << "    X    |    Y    |    Z    " << std::endl;
  for (UShort_t ddl = 0; ddl < 3; ddl++) { 
    out << "------+-------+-------+---------+"
	<< "----------+------+--------+------------+------+"
	<< "---------+---------+---------" << std::endl;
    for (UShort_t ib = 0; ib < 4; ib++) { 
      if (ddl == 0 && (ib == 1 || ib == 3)) continue;
      UShort_t board = ib + (ib < 0x2 ? 0 : 0x10-2);
      for (UShort_t altro = 0; altro < 3; altro++) { 
	UShort_t nCh = (altro == 1 ? 8 : 16);
	for (UShort_t chan = 0; chan < nCh; chan++) { 
	  UShort_t det, sec;
	  Short_t baseStr;
	  Char_t  ring;
	  det = map->DDL2Detector(ddl);
	  if (!map->Channel2StripBase(board, altro, chan, ring, sec, baseStr))
	    continue;
	  
	  Double_t x, y, z;
	  geom->Detector2XYZ(det, ring, sec, baseStr, x, y, z);

	  Double_t phi = TMath::ATan2(y,x) * 180 / TMath::Pi();
	  if (phi < 0) phi += 360;
	  
	  out << " " << std::setprecision(4)
	      << std::setw(4)  << (3072) + ddl << " | " 
	      << std::setw(5)  << board   << " | " 
	      << std::setw(5)  << altro   << " | " 
	      << std::setw(7)  << chan    << " | " 
	      << std::setw(8)  << det     << " | " 
	      << std::setw(4)  << ring    << " | " 
	      << std::setw(6)  << sec     << " | " 
	      << std::setw(10) << baseStr << " | " 
	      << std::setw(4)  << phi     << " | " 
	      << std::setw(7)  << x       << " | " 
	      << std::setw(7)  << y       << " | " 
	      << std::setw(7)  << z       << std::endl;
	} // Chan
      } // Altro
    } // Board 
  } // DDL
}

  
    
      
	
	  
	  
	  
	  
      
