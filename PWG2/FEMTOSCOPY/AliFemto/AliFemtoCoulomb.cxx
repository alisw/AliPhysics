///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoCoulomb: This is a Coulomb correction class which             //
//  1. Reads in the dat from a file                                      //  
//  2. Performs a linear interpolation in R and creates any array of     //
//     interpolations                                                    //
//  3. Interpolates in eta and returns the Coulomb correction to user    //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliFemtoCoulomb.h"
//#include "Stiostream.h"
#include <stdio.h>
#include <cassert>
//#include "PhysicalConstants.h"
#define fine_structure_const 0.00729735

#ifdef __ROOT__
ClassImp(AliFemtoCoulomb)
#endif

AliFemtoCoulomb::AliFemtoCoulomb() :
  fFile(""),
  fRadius(-1.0),
  fZ1Z2(1.0),
  fNLines(0)
{
  // Default constructor
  fFile = "/afs/rhic/star/hbt/coul/AliFemtoCorrectionFiles/correctionpp.dat";
  if (!fFile) {
    cout << " No file, dummy!" << endl;
    assert(0);
  }
  cout << "You have 1 default Coulomb correction!" << endl;
}

AliFemtoCoulomb::AliFemtoCoulomb(const AliFemtoCoulomb& aCoul) :
  fFile(aCoul.fFile),
  fRadius(aCoul.fRadius),
  fZ1Z2(aCoul.fZ1Z2),
  fNLines(0)
{
  // copy constructor
  CreateLookupTable(fRadius);
}

AliFemtoCoulomb::AliFemtoCoulomb(const char* readFile, const double& radius, const double& charge) :
  fFile(readFile),
  fRadius(radius),
  fZ1Z2(0),
  fNLines(0)
{
  // constructor with explicit filename
  fFile = readFile;
  fRadius = radius;
  CreateLookupTable(fRadius);
  fZ1Z2 = charge;
  cout << "You have 1 Coulomb correction!" << endl;
}

AliFemtoCoulomb::~AliFemtoCoulomb() {
  // destructor
}

AliFemtoCoulomb& AliFemtoCoulomb::operator=(const AliFemtoCoulomb& aCoul)
{
  // assignment operator
  if (this == &aCoul)
    return *this;

  fFile = aCoul.fFile;
  fRadius = aCoul.fRadius;
  fZ1Z2 = aCoul.fZ1Z2;

  CreateLookupTable(fRadius);
  
  return *this;
}


void AliFemtoCoulomb::SetRadius(const double& radius) {
  // set the coulomb radius
  cout << " AliFemtoCoulomb::setRadius() " << endl;
  fRadius = radius;
  CreateLookupTable(fRadius);
}

double AliFemtoCoulomb::GetRadius() const {
  // return coulomb radius
  return (fRadius);
}

void AliFemtoCoulomb::SetFile(const char* readFile) {
  // set the filename with coulomb calculations
  cout << " AliFemtoCoulomb::SetFile() " << endl;
  fFile = readFile;
  // Create new lookup table since file has changed
  if (fRadius>0.0) {
    CreateLookupTable(fRadius);
  }
}

void AliFemtoCoulomb::SetChargeProduct(const double& charge) {
  // set pair charge
  cout << " AliFemtoCoulomb::SetChargeProduct() " << endl;
  if ( fZ1Z2!=charge ) { 
    fZ1Z2 = charge;
    if ( fZ1Z2>0 ) {
      fFile = "/afs/rhic/star/hbt/coul/AliFemtoCorrectionFiles/correctionpp.dat";
    }
    else {
      fFile = "/afs/rhic/star/hbt/coul/AliFemtoCorrectionFiles/correctionpm.dat";
    }
    CreateLookupTable(fRadius);
  }
}

void AliFemtoCoulomb::CreateLookupTable(const double& radius) {
  // Read radii from fFile
  // Create array(pair) of linear interpolation between radii
  cout << " AliFemtoCoulomb::CreateLookupTable() " << endl;

  if (radius<0.0) {
    cout << " AliFemtoCoulomb::CreateLookupTable -> NEGATIVE RADIUS " << endl;
    cout << "  call AliFemtoCoulomb::SetRadius(r) with positive r " << endl;
    cerr << " AliFemtoCoulomb::CreateLookupTable -> NEGATIVE RADIUS " << endl;
    cerr << "  call AliFemtoCoulomb::SetRadius(r) with positive r " << endl;
    assert(0);
  }
  ifstream mystream(fFile);
  if (!mystream) {
    cout << "Could not open file" << endl;
    assert(0);
  }
  else {
    cout << "Input correction file opened" << endl;
  }

  static char tempstring[2001];
  static float radii[2000];
  static int tNRadii = 0;
  tNRadii = 0;
  if (!mystream.getline(tempstring,2000)) {
    cout << "Could not read radii from file" << endl;
    assert(0);
  }
  for (unsigned int ii=0; ii<strlen(tempstring); ii++) {
    while (tempstring[ii]==' ') ii++;
    sscanf(&tempstring[ii++],"%f",&radii[++tNRadii]);
    while ( tempstring[ii]!=' ' && (ii)<strlen(tempstring) )ii++;
  }
  cout << " Read " << tNRadii << " radii from file" << endl;

  static double tLowRadius = -1.0;
  static double tHighRadius = -1.0;
  static int tLowIndex = 0;
  tLowRadius = -1.0;
  tHighRadius = -1.0;
  tLowIndex = 0;
  for(int iii=1; iii<=tNRadii-1; iii++) { // Loop to one less than #radii
    if ( radius >= radii[iii] && radius <= radii[iii+1] ) {
      tLowRadius = radii[iii];
      tHighRadius = radii[iii+1];
      tLowIndex = iii;
    }
  }
  if ( (tLowRadius < 0.0) || (tHighRadius < 0.0) ) {
    cout << "AliFemtoCoulomb::CreateLookupTable --> Problem interpolating radius" << endl;
    cout << "  Check range of radii in lookup file...." << endl;
    cerr << "AliFemtoCoulomb::CreateLookupTable --> Problem interpolating radius" << endl;
    cerr << "  Check range of radii in lookup file...." << endl;
    assert(0);
  }

  static double corr[100];           // array of corrections ... must be > tNRadii
  fNLines = 0;
  static double tempEta = 0;
  tempEta = 0;
  while (mystream >> tempEta) {
    for (int i=1; i<=tNRadii; i++) {
      mystream >> corr[i];
    }
    static double tLowCoulomb = 0;
    static double tHighCoulomb = 0;
    static double nCorr = 0;
    tLowCoulomb = corr[tLowIndex];
    tHighCoulomb = corr[tLowIndex+1];
    nCorr = ( (radius-tLowRadius)*tHighCoulomb+(tHighRadius-radius)*tLowCoulomb )/(tHighRadius-tLowRadius);
      fEta[fNLines] = tempEta;     // Eta
      fCoulomb[fNLines] = nCorr;   // Interpolated Coulomb correction for radius
      fNLines++;
  }
  mystream.close();
  cout << "Lookup Table is created with " << fNLines << " points" << endl;
}

double AliFemtoCoulomb::CoulombCorrect(const double& eta) {
  // Interpolates in eta
  if (fRadius < 0.0) {
    cout << "AliFemtoCoulomb::CoulombCorrect(eta) --> Trying to correct for negative radius!" << endl;
    cerr << "AliFemtoCoulomb::CoulombCorrect(eta) --> Trying to correct for negative radius!" << endl;
    assert(0);
  }
  static int middle=0;
  middle=int( (fNLines-1)/2 );
  if (eta*fEta[middle]<0.0) {
    cout << "AliFemtoCoulomb::CoulombCorrect(eta) --> eta: " << eta << " has wrong sign for data file! " << endl;
    cerr << "AliFemtoCoulomb::CoulombCorrect(eta) --> eta: " << eta << " has wrong sign for data file! " << endl;
    assert(0);
  }

  static double tCorr = 0;
  tCorr = -1.0;
  
  if ( (eta>fEta[0]) && (fEta[0]>0.0) ) {
    tCorr = fCoulomb[0];
    return (tCorr);
  }
  if ( (eta<fEta[fNLines-1]) && (fEta[fNLines-1]<0.0) ) {
    tCorr = fCoulomb[fNLines-1];
    return (tCorr);
  }
  // This is a binary search for the bracketing pair of data points
  static int high = 0;
  static int low = 0;
  static int width = 0;
  high = fNLines-1;
  low = 0;
  width = high-low;
  middle = int(width/2.0); // Was instantiated above
  while (middle > 0) {
    if (fEta[low+middle] < eta) {
      // eta is in the 1st half
      high-=middle;
      width = high-low;
      middle = int(width/2.0);
    }
    else {
      // eta is in the 2nd half
      low+=middle;
      width = high-low;
      middle = int(width/2.0);
    }
  }
  // Make sure we found the right one
  if ( (fEta[low] >= eta) && (eta >= fEta[low+1]) ) {
    static double tLowEta = 0;
    static double tHighEta = 0;    
    static double tLowCoulomb = 0;
    static double tHighCoulomb = 0;
    tLowEta = fEta[low];
    tHighEta = fEta[low+1];    
    tLowCoulomb = fCoulomb[low];
    tHighCoulomb = fCoulomb[low+1];
    //      cout << tLowEta << " *** Eta *** " << tHighEta << endl;
    //      cout << tLowCoulomb << " *** Coulomb *** " << tHighCoulomb << endl;
    tCorr = ( (eta-tLowEta)*tHighCoulomb+(tHighEta-eta)*tLowCoulomb )/(tHighEta-tLowEta);
  }
  if (tCorr<0.0) {
    cout << "AliFemtoCoulomb::CoulombCorrect(eta) --> No correction" << endl;
    cout << "  Check range of eta in file: Input eta  " << eta << endl;
    cerr << "AliFemtoCoulomb::CoulombCorrect(eta) --> No correction" << endl;
    cerr << "  Check range of eta in file: Input eta  " << eta << endl;
    assert(0);
  } 
  return (tCorr);

}

double AliFemtoCoulomb::CoulombCorrect(const double& eta,
				const double& radius) {
  // Checks radii ... input radius and fRadius
  // Calls createLookupTable if neccessary
  // Interpolate(linear) between etas in the created lookup table

  if (radius < 0.0) {
    if (fRadius < 0.0) {
      // Both radii are negative
      cout << "AliFemtoCoulomb::CoulombCorrect(eta,r) --> input and member radii are negative!" << endl;
      cerr << "AliFemtoCoulomb::CoulombCorrect(eta,r) --> input and member radii are negative!" << endl;
      assert(0);
    }
  }
  else {
    // radius > 0.0
    if (radius == fRadius) {
      // Both radii are positive and equal
      //      cout << "Radii are the same!!!" << endl;
    }
    else {
      // Both radii are positive but not equal
      fRadius = radius;
      CreateLookupTable(fRadius);
    }
  }

  // Interpolate in eta
  return ( CoulombCorrect(eta) );
}

double AliFemtoCoulomb::CoulombCorrect(const AliFemtoPair* pair) {
  return ( CoulombCorrect( Eta(pair) ) );;
}

double AliFemtoCoulomb::CoulombCorrect(const AliFemtoPair* pair, const double& radius) {
  return ( CoulombCorrect( Eta(pair),radius ) );
}

double AliFemtoCoulomb::Eta(const AliFemtoPair* pair) {
  // calculate eta
  static double px1,py1,pz1,px2,py2,pz2;
  static double px1new,py1new,pz1new;
  static double px2new,py2new,pz2new;
  static double vx1cms,vy1cms,vz1cms;
  static double vx2cms,vy2cms,vz2cms;
  static double tVcmsX,tVcmsY,tVcmsZ;
  static double dv = 0.0;
  static double e1,e2,e1new,e2new;
  static double psi,theta;
  static double beta,gamma;
  static double tVcmsXnew;

  px1 = pair->Track1()->FourMomentum().px();
  py1 = pair->Track1()->FourMomentum().py();
  pz1 = pair->Track1()->FourMomentum().pz();
  e1 = pair->Track1()->FourMomentum().e();
  px2 = pair->Track2()->FourMomentum().px();
  py2 = pair->Track2()->FourMomentum().py();
  pz2 = pair->Track2()->FourMomentum().pz();
  e2 = pair->Track2()->FourMomentum().e();
  
  tVcmsX = ( px1+px2 )/( e1+e2 );
  tVcmsY = ( py1+py2 )/( e1+e2 );
  tVcmsZ = ( pz1+pz2 )/( e1+e2 );
  // Rotate tVcms to x-direction
  psi = atan(tVcmsY/tVcmsX);
  tVcmsXnew = tVcmsX*cos(psi)+tVcmsY*sin(psi);
  tVcmsX = tVcmsXnew;
  theta = atan(tVcmsZ/tVcmsX);
  tVcmsXnew = tVcmsX*cos(theta)+tVcmsZ*sin(theta);
  tVcmsX = tVcmsXnew;
  // Gamma and Beta
  beta = tVcmsX;
  gamma = 1.0/::sqrt( 1.0-beta*beta );

  // Rotate p1 and p2 to new frame
  px1new = px1*cos(psi)+py1*sin(psi);
  py1new = -px1*sin(psi)+py1*cos(psi);
  px1 = px1new;
  px1new = px1*cos(theta)+pz1*sin(theta);
  pz1new = -px1*sin(theta)+pz1*cos(theta);
  px1 = px1new;
  py1 = py1new;
  pz1 = pz1new;

  px2new = px2*cos(psi)+py2*sin(psi);
  py2new = -px2*sin(psi)+py2*cos(psi);
  px2 = px2new;
  px2new = px2*cos(theta)+pz2*sin(theta);
  pz2new = -px2*sin(theta)+pz2*cos(theta);
  px2 = px2new;
  py2 = py2new;
  pz2 = pz2new;

  // Lorentz transform the x component and energy
  e1new = gamma*e1 - gamma*beta*px1;
  px1new = -gamma*beta*e1 + gamma*px1;
  e2new = gamma*e2 - gamma*beta*px2;
  px2new = -gamma*beta*e2 + gamma*px2;
  px1 = px1new;
  px2 = px2new;

  // New velocities
  vx1cms = px1/e1new;
  vy1cms = py1/e1new;
  vz1cms = pz1/e1new;
  vx2cms = px2/e2new;
  vy2cms = py2/e2new;
  vz2cms = pz2/e2new;

  // Velocity difference in CMS frame
  dv = ::sqrt( (vx1cms-vx2cms)*(vx1cms-vx2cms) +
	     (vy1cms-vy2cms)*(vy1cms-vy2cms) +
	     (vz1cms-vz2cms)*(vz1cms-vz2cms) );

  return ( fZ1Z2*fine_structure_const/(dv) );
}

TH1D* AliFemtoCoulomb::CorrectionHistogram(const double& mass1, const double& mass2, const int& nBins, 
						const double& low, const double& high) {
  // return correction histogram

  if ( mass1!=mass2 ) {
    cout << "Masses not equal ... try again.  No histogram created." << endl;
    assert(0);
  }
  TH1D* correction = new TH1D("correction","Coulomb correction",nBins,low,high);
  const double kReducedMass = mass1*mass2/(mass1+mass2);
  double qInv = low;
  //double dQinv = (high-low)/( (double)nBins );
  double eta;
  for (int ii=1; ii<=nBins; ii++) 
    {
      qInv = correction->GetBinCenter(ii);
      eta = 2.0*fZ1Z2*kReducedMass*fine_structure_const/( qInv );
      CoulombCorrect( eta );
      correction->Fill( qInv, CoulombCorrect(eta,fRadius) );
    }

  return (correction);
}

#ifdef __ROOT__
TH1D* AliFemtoCoulomb::CorrectionHistogram(const TH1D* histo, const double mass) {
  // return correction histogram - 1D case
  TH1D* correction = (TH1D*) ((TH1D*)histo)->Clone();
  correction->Reset();
  correction->SetDirectory(0);
  int    nBins = correction->GetXaxis()->GetNbins();
  const double kReducedMass = 0.5*mass;
  double qInv;
  double eta;
  for (int ii=1; ii<=nBins; ii++) 
    {
      qInv = correction->GetBinCenter(ii);
      eta = 2.0*fZ1Z2*kReducedMass*fine_structure_const/( qInv );
      correction->Fill( qInv, CoulombCorrect(eta,fRadius) );
    }

  return (correction);
}

TH3D* AliFemtoCoulomb::CorrectionHistogram(const TH3D* histo, const double mass) {
  // return correction histogram - 3D case
  TH3D* correction = (TH3D*) ((TH3D*)histo)->Clone();
  correction->Reset();
  correction->SetDirectory(0);
  int    nBinsX = correction->GetXaxis()->GetNbins();
  int    nBinsY = correction->GetYaxis()->GetNbins();
  int    nBinsZ = correction->GetZaxis()->GetNbins();
  const double kReducedMass = 0.5*mass;
  double eta;
  double qInv;
  int binNumber;
  for (int ii=1; ii<=nBinsX; ii++) { 
    for (int iii=1; iii<=nBinsY; iii++) {
      for (int iv=1; iv<=nBinsZ; iv++) {
	binNumber = histo->GetBin(ii,iii,iv);
	qInv = histo->GetBinContent(binNumber);
	eta = 2.0*fZ1Z2*kReducedMass*fine_structure_const/( qInv );
	correction->SetBinContent(binNumber, CoulombCorrect(eta,fRadius) );
      }
    }
  }
  return (correction);
}
#endif

double AliFemtoCoulomb::CoulombCorrect(const double& mass, const double& charge,
				    const double& radius, const double& qInv) {
  // return correction factor
  fRadius = radius;
  fZ1Z2 = charge;
  const double kReducedMass = 0.5*mass; // must be same mass particles
  double eta = 2.0*fZ1Z2*kReducedMass*fine_structure_const/( qInv );
  return ( CoulombCorrect(eta,fRadius) );
}
