#include <cstdlib>
#include <cstring>
#include "Hypersurface.h"

Hypersurface::Hypersurface(const char *dirname) {
// Read FOHSI.txt file and init parameters
  char *dirfull;
  dirfull = (char *) malloc(sizeof(char)*strlen(dirname) + 5);
  sprintf(dirfull, "%s", dirname);
  if (strlen(dirfull) > 0) 
    if (dirfull[strlen(dirfull)-1] != '/')
      strcat(dirfull, "/");
  FName = (char *) malloc(sizeof(char)* strlen(dirfull)+ 50);
  sprintf(FName,"%sFOHSI.txt", dirfull);
  HSFile = new ifstream;
  HSFile->open(FName);
  if(HSFile->is_open()) {
    HSFile->seekg (0, ios::beg);
// phi
    HSFile->getline(buff,100);
    Np=atoi(buff);
    HSFile->getline(buff,100);
    ip=atof(buff);
    HSFile->getline(buff,100);
    fp=atof(buff);
// zeta
    HSFile->getline(buff,100);
    Nz=atoi(buff);
    HSFile->getline(buff,100);
    iz=atof(buff);
    HSFile->getline(buff,100);
    fz=atof(buff);
// freeze-out temperature
    HSFile->getline(buff,100);
    TFO=atof(buff);
// hydro initial time
    HSFile->getline(buff,100);
    tau0=atof(buff);

    HSFile->close();
    dp = (fp-ip)/(Np-1);
    dz = (fz-iz)/(Nz-1);

// create 2D array
    aArr   = new double* [Np];
    vArr   = new double* [Np];
    dArr   = new double* [Np];
    DpdArr = new double* [Np];
    DzdArr = new double* [Np];
    for(i=0;i<Np;i++) {
      aArr[i]   = new double [Nz];
      vArr[i]   = new double [Nz];
      dArr[i]   = new double [Nz];
      DpdArr[i] = new double [Nz];
      DzdArr[i] = new double [Nz];
    }
  } else {
    Np = 0; ip = 0.0; fp = 0.0; dp = 1.0;
    Nz = 0; iz = 0.0; fz = 0.0; dz = 1.0;
    TFO = 0.0;
    aArr   = NULL;
    vArr   = NULL;
    dArr   = NULL;
    DpdArr = NULL;
    DzdArr = NULL;
  }
  delete HSFile;

  // Read FOHSa.txt file and fill array
  sprintf(FName,"%sFOHSa.txt", dirfull);
  //  FName = "FOHSa.txt";
  HSFile = new ifstream;
  HSFile->open(FName);
  if(HSFile->is_open()) {
    HSFile->seekg (0, ios::beg);
    for(i=0;i<Np;i++) {
      for(j=0;j<Nz;j++) {
        HSFile->getline(buff,100);
        aArr[i][j]=atof(buff);
      }
    }
    HSFile->close();
  }
  delete HSFile;

  // Read FOHSv.txt file and fill array
  sprintf(FName,"%sFOHSv.txt", dirfull);
  //  FName = "FOHSv.txt";
  HSFile = new ifstream;
  HSFile->open(FName);
  if(HSFile->is_open()) {
    HSFile->seekg (0, ios::beg);
    for(i=0;i<Np;i++) {
      for(j=0;j<Nz;j++) {
        HSFile->getline(buff,100);
        vArr[i][j]=atof(buff);
      }
    }
    HSFile->close();
  }
  delete HSFile;

  // Read FOHSd.txt file and fill array
  sprintf(FName,"%sFOHSd.txt", dirfull);
  //  FName = "FOHSd.txt";
  HSFile = new ifstream;
  HSFile->open(FName);
  if(HSFile->is_open()) {
    HSFile->seekg (0, ios::beg);
    for(i=0;i<Np;i++) {
      for(j=0;j<Nz;j++) {
        HSFile->getline(buff,100);
        dArr[i][j]=atof(buff);
      }
    }
    HSFile->close();
  }
  delete HSFile;

  // Read FOHSDpd.txt file and fill array
  sprintf(FName,"%sFOHSDpd.txt", dirfull);
  //  FName = "FOHSDpd.txt";
  HSFile = new ifstream;
  HSFile->open(FName);
  if(HSFile->is_open()) {
    HSFile->seekg (0, ios::beg);
    for(i=0;i<Np;i++) {
      for(j=0;j<Nz;j++) {
        HSFile->getline(buff,100);
        DpdArr[i][j]=atof(buff);
      }
    }
    HSFile->close();
  }
  delete HSFile;

  // Read FOHSDzd.txt file and fill array
  sprintf(FName,"%sFOHSDzd.txt", dirfull);
  //  FName = "FOHSDzd.txt";
  HSFile = new ifstream;
  HSFile->open(FName);
  if(HSFile->is_open()) {
    HSFile->seekg (0, ios::beg);
    for(i=0;i<Np;i++) {
      for(j=0;j<Nz;j++) {
        HSFile->getline(buff,100);
        DzdArr[i][j]=atof(buff);
      }
    }
    HSFile->close();
  }
  delete HSFile;
  free (FName);
  free (dirfull);
}

Hypersurface::Hypersurface(void) {
  Hypersurface("./");
}

Hypersurface::Hypersurface(const Hypersurface &aSurf)
{
  i = aSurf.i;
  Hypersurface("./");
  
}

Hypersurface& Hypersurface::operator=(const Hypersurface &aSurf)
{
  if (this != &aSurf) {
    i = aSurf.i;
    Hypersurface("./");
  }
  
  return *this;
}

Hypersurface::~Hypersurface(void) {
  for(i=0;i<Np;i++) {
    delete[] aArr[i];
    delete[] vArr[i];
    delete[] dArr[i];
    delete[] DpdArr[i];
    delete[] DzdArr[i];
  }
  delete[] aArr;
  delete[] vArr;
  delete[] dArr;
  delete[] DpdArr;
  delete[] DzdArr;
}

double Hypersurface::fahs(double p, double z) {
  i=int(p/dp);
  j=int(z/dz);
  if(i>=Np-1) i=Np-2;
  if(j>=Nz-1) j=Nz-2;
  return (aArr[i][j]   * (i+1-p/dp) + aArr[i+1][j]   * (p/dp-i)) * (j+1-z/dz) +
         (aArr[i][j+1] * (i+1-p/dp) + aArr[i+1][j+1] * (p/dp-i)) * (z/dz-j);
}

double Hypersurface::fvhs(double p, double z) {
  i=int(p/dp);
  j=int(z/dz);
  if(i>=Np-1) i=Np-2;
  if(j>=Nz-1) j=Nz-2;
  return (vArr[i][j]   * (i+1-p/dp) + vArr[i+1][j]   * (p/dp-i)) * (j+1-z/dz) +
         (vArr[i][j+1] * (i+1-p/dp) + vArr[i+1][j+1] * (p/dp-i)) * (z/dz-j);
}

double Hypersurface::fdhs(double p, double z) {
  i=int(p/dp);
  j=int(z/dz);
  if(i>=Np-1) i=Np-2;
  if(j>=Nz-1) j=Nz-2;
  return (dArr[i][j]   * (i+1-p/dp) + dArr[i+1][j]   * (p/dp-i)) * (j+1-z/dz) +
         (dArr[i][j+1] * (i+1-p/dp) + dArr[i+1][j+1] * (p/dp-i)) * (z/dz-j);
}

double Hypersurface::fDpdhs(double p, double z) {
  i=int(p/dp);
  j=int(z/dz);
  if(i>=Np-1) i=Np-2;
  if(j>=Nz-1) j=Nz-2;
  return (DpdArr[i][j]   * (i+1-p/dp) + DpdArr[i+1][j]   * (p/dp-i)) * (j+1-z/dz) +
         (DpdArr[i][j+1] * (i+1-p/dp) + DpdArr[i+1][j+1] * (p/dp-i)) * (z/dz-j);
}

double Hypersurface::fDzdhs(double p, double z) {
  i=int(p/dp);
  j=int(z/dz);
  if(i>=Np-1) i=Np-2;
  if(j>=Nz-1) j=Nz-2;
  return (DzdArr[i][j]   * (i+1-p/dp) + DzdArr[i+1][j]   * (p/dp-i)) * (j+1-z/dz) +
         (DzdArr[i][j+1] * (i+1-p/dp) + DzdArr[i+1][j+1] * (p/dp-i)) * (z/dz-j);
}
