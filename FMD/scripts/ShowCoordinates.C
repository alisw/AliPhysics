/**
 * A script to dump the physical coordinates as given by the
 * geometry. 
 * 
 */
#include <iomanip>

/** 
 * Get the physical coordinates of a strip 
 * 
 * @param det Detector
 * @param rng Ring
 * @param sec Sector
 * @param str Strip
 */
void
PhysicalCoordinates(UShort_t det, Char_t rng, UShort_t sec, UShort_t str)
{
  Double_t x, y, z;
  AliFMDGeometry::Instance()->Detector2XYZ(det, rng, sec, str, x, y, z);
  Double_t phi   = TMath::ATan2(y, x);
  Double_t r     = TMath::Sqrt(x * x + y * y);
  Double_t theta = TMath::ATan2(r, z);
  if (theta < 0) theta += TMath::Pi();
  Double_t eta   = -TMath::Log(TMath::Tan(theta / 2));
  Double_t deg   = 180. / TMath::Pi();
  
  std::cout << det << rng << "[" 
	    << std::setw(2) << sec << "," 
	    << std::setw(3) << str << "] | "
	    << std::setw(9) << x << "," 
	    << std::setw(9) << y << "," 
	    << std::setw(9) << z << " | " 
	    << std::setw(9) << phi * deg << "," 
	    << std::setw(9) << theta * deg << "," 
	    << std::setw(9) << eta << std::endl;
}

/** 
 * Show coordinates of all strips 
 * 
 */
void
ShowCoordinates()
{
  AliFMDGeometry::Instance()->Init();
  AliFMDGeometry::Instance()->InitTransformations();
  std::cout << std::setw(1+1+1+2+1+3+1) << "Detector" << " | " 
	    << std::setw(9+1+9+1+9) << "Cartisian Coords" << " | " 
	    << std::setw(9+1+9+1+9) << "phi,theta,eta" << "\n"
	    << std::setfill('-') 
	    << std::setw(1+1+1+2+1+3+1+1+1)  << "+"
	    << std::setw(1+9+1+9+1+9+1+1)  << "+"
	    << std::setw(1+9+1+9+1+9+1+1)  << "+" 
	    << std::setfill(' ') << std::endl;
  for (UShort_t d = 1; d <= 3; d++) {
    UShort_t nrng = (d == 1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nrng; ir++) { 
      Char_t r = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (r == 'I' ?  20 :  40);
      UShort_t nstr = 1; // (r == 'I' ? 512 : 256);
      for (UShort_t s = 0; s < nsec; s++) { 
	for (UShort_t t = 0; t < nstr; t++) {
	  PhysicalCoordinates(d, r, s, t);
	}
      }
    }
  }
}
//
// EOF
//
