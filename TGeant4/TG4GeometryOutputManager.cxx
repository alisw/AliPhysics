// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "TG4GeometryOutputManager.h"
#include <iostream.h>
#include <iomanip.h>

TG4GeometryOutputManager::TG4GeometryOutputManager() {
//
}

TG4GeometryOutputManager::~TG4GeometryOutputManager() {
//
}

// public methods

void TG4GeometryOutputManager::OpenFile(G4String filePath)
{ 
// Opens output files.
// ---

  G4cout << "TG4GeometryOutputManager::OpenFile: " << filePath << endl;
  
  fOutFile.open(filePath, ios::out, filebuf::openprot); 
  
  if (!fOutFile) {
    G4String text = "Cannot open ";
    text = text + filePath;
    TG4Globals::Warning(text);  
  }
  
  // use FORTRAN compatibility output
  fOutFile << setiosflags(ios::showpoint | ios::uppercase);
}


void TG4GeometryOutputManager::CloseFile()
{ 
// Closess output files.
// ---

  fOutFile.close(); 
}


void TG4GeometryOutputManager::WriteGsvolu( 
              G4String vname, G4String shape, G4int nmed, G4double* Rpar,
              G4int npar)
{
// from fortran (g3routines.F)
// write(fmt,'(A,I2,A)')'(a4,1x,a6,1x,a4,1x,a4,2i5,',max(npar,1),
//>      '(1x,e16.8))'
// write(lunlist,fmt) context, rname, name, shape, nmed, npar,
//+      (par(k),k=1,npar)
// ---

  G4String context("----");
  G4String rname("GSVOLU");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << vname    << space
	  << shape
	  << setw(5) << nmed  
	  << setw(5) << npar;
  for (G4int i=0; i<npar; i++)	      
    fOutFile << space << setw(16) << setprecision(8) << Rpar[i];
  fOutFile << G4endl;   
}


void TG4GeometryOutputManager::WriteGspos(
             G4String vname, G4int num, G4String vmoth, G4double x,
             G4double y, G4double z, G4int irot, G4String vonly)
{	     
// from fortran (g3routines.F)
// write(lunlist,
//+      '(a4,1x,a6,1x,a4,i5,1x,a4,3(1x,e16.8),i5,1x,a4)')
//+      context, rname, name, num, moth, x, y, z, irot, only
// ---

  G4String context("----");
  G4String rname("GSPOS");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << vname   << space
	  << setw(5) << num << space 
	  << vmoth   << space
          << setw(16) << setprecision(8) << x << space
          << setw(16) << setprecision(8) << y << space
          << setw(16) << setprecision(8) << z
	  << setw(5) << irot << space
	  << vonly 
	  << G4endl;
}	     
  

void TG4GeometryOutputManager::WriteGsposp(
              G4String vname, G4int num, G4String vmoth, G4double x,
              G4double y, G4double z, G4int irot, G4String vonly,
              G4double pars[], G4int npar)
{
// from fortran (g3routines.F)
// write(fmt,'(A,A,I2,A)')
//+      '(a4,1x,a6,1x,a4,i5,1x,a4,3(1x,e16.8),',       
//+      'i5,1x,a4,i5,',max(npar,1),'(1x,e16.8))'
// write(lunlist,fmt)
//+      context, rname, name, num, moth, x, y, z, irot, only,
//+      npar,
//+      (par(k),k=1,npar)
// ---

  G4String context("----");
  G4String rname("GSPOSP");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << vname   << space
	  << setw(5) << num << space 
	  << vmoth   << space
          << setw(16) << setprecision(8) << x << space
          << setw(16) << setprecision(8) << y << space
          << setw(16) << setprecision(8) << z
	  << setw(5) << irot << space
	  << vonly 
	  << setw(5) << npar;
  for (G4int i=0; i<npar; i++)	      
    fOutFile << space << setw(16) << setprecision(8) << pars[i];
  fOutFile << G4endl;
}	      


void TG4GeometryOutputManager::WriteGsrotm(
              G4int irot, G4double theta1, G4double phi1,
              G4double theta2, G4double phi2, G4double theta3, G4double phi3)
{
// from fortran (g3routines.F)
// write(lunlist,
//+      '(a4,1x,a6,i5,6f11.5)')
//+      context, rname, irot, theta1, phi1, theta2, phi2,
//+      theta3, phi3
// ---
  
  G4String context("----");
  G4String rname("GSROTM");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << setw(5) << irot
          << setw(11) << setprecision(5) << theta1
          << setw(11) << setprecision(5) << phi1
          << setw(11) << setprecision(5) << theta2
          << setw(11) << setprecision(5) << phi2
          << setw(11) << setprecision(5) << theta3
          << setw(11) << setprecision(5) << phi3
	  << G4endl;
}	  


void TG4GeometryOutputManager::WriteGsdvn(
              G4String vname, G4String vmoth, G4int ndiv, G4int iaxis)
{
// from fortran (g3routines.F)
// write(lunlist,
//+      '(a4,1x,a6,1x,a4,1x,a4,i5,i3)')
//+      context, rname, name, moth, ndiv, iaxis
// ---

  G4String context("----");
  G4String rname("GSDVN");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << vname    << space
	  << vmoth   << space
	  << setw(5) << ndiv  
	  << setw(5) << iaxis
          << G4endl;
}	     


void TG4GeometryOutputManager::WriteGsdvn2(
               G4String vname, G4String vmoth, G4int ndiv, G4int iaxis,
               G4double c0, G4int numed)
{
// from fortran (g3routines.F)
// write(lunlist,
//+      '(a4,1x,a6,1x,a4,1x,a4,i5,i3,(1x,e16.8),i5)')
//+      context, rname, name, moth, ndiv, iaxis, c0, numed
// ---

  G4String context("----");
  G4String rname("GSDVN2");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << vname    << space
	  << vmoth   << space
	  << setw(5) << ndiv  
	  << setw(5) << iaxis << " "
          << setw(16) << setprecision(8) << c0
	  << setw(5) << numed
	  << G4endl;
}	     


void TG4GeometryOutputManager::WriteGsdvt(
             G4String vname, G4String vmoth, G4double step, G4int iaxis,
             G4int numed, G4int ndvmx) 
{
// from fortran (g3routines.F)
// write(lunlist,
// +    '(a4,1x,a6,1x,a4,1x,a4,(1x,e16.8),3i5)')
// +    context, rname, name, moth, step, iaxis, numed, ndvmx
// ---

  G4String context("----");
  G4String rname("GSDVT");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << vname    << space
	  << vmoth   << space
          << setw(16) << setprecision(8) << step
	  << setw(5) << iaxis
	  << setw(5) << numed
	  << setw(5) << ndvmx
	  << G4endl;
}	     


void TG4GeometryOutputManager::WriteGsdvt2(
               G4String vname, G4String vmoth, G4double step, G4int iaxis,
               G4double c0, G4int numed, G4int ndvmx)
{	       
// from fortran (g3routines.F)
// write(lunlist,
//+      '(a4,1x,a6,1x,a4,1x,a4,(1x,e16.8),i3,(1x,e16.8),2i5)')
//+      context, rname, name, moth, step, iaxis, c0, numed, ndvmx
// ---

  G4String context("----");
  G4String rname("GSDVT2");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << vname    << space
	  << vmoth   << space
          << setw(16) << setprecision(8) << step
	  << setw(3) << iaxis << space
          << setw(16) << setprecision(8) << c0
	  << setw(5) << numed
	  << setw(5) << ndvmx
	  << endl;
}	     


void TG4GeometryOutputManager::WriteGsdvx(
             G4String name, G4String moth, G4int ndiv, G4int iaxis,
             G4double step, G4double c0, G4int numed, G4int ndvmx)
{
// from fortran (g3routines.F)
// write(lunlist,
// +     '(a4,1x,a6,1x,a4,1x,a4,i5,i3,2(1x,e16.8),2i5)')
// +     context, rname, name, moth, ndiv, iaxis,step, c0,
// +     numed, ndvmx
// ---

  G4String context("----");
  G4String rname("GSDVX");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << name    << space
	  << moth    << space
	  << setw(5) << ndiv
	  << setw(3) << iaxis << space
          << setw(16) << setprecision(8) << step << space
          << setw(16) << setprecision(8) << c0
	  << setw(5) << numed
	  << setw(5) << ndvmx
	  << endl;
}	     


void TG4GeometryOutputManager::WriteGsmate(
              G4int imate, G4String name, G4double ain, G4double zin,
              G4double densin, G4double radl, G4int nwbf, G4double* ubuf)
{
// from fortran (g3routines.F)
// write(fmt,'(A,I3,A)')
//+      '(a4,1x,a6,i5,1x,''"'',a,''"'',4(1x,e16.8),i3,',
//+      max(nwbf,1),'(1x,e16.8))'
// write(lunlist,fmt)
//+      context, rname, imate, name, a, z, dens, radl,
//+      nwbf, (ubf(k), k=1,nwbf)
// ---

  G4String context("----");
  G4String rname("GSMATE");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << setw(5) << imate << space 
	  << '"' << name << '"' << space
          << setw(16) << setprecision(8) << ain << space
          << setw(16) << setprecision(8) << zin << space
          << setw(16) << setprecision(8) << densin << space
          << setw(16) << setprecision(8) << radl
	  << setw(3) << nwbf;
  for (G4int i=0; i<nwbf; i++)	      
    fOutFile << space << setw(16) << setprecision(8) << ubuf[i];
  fOutFile << endl;
}	      


void TG4GeometryOutputManager::WriteGsmixt(
              G4int imate, G4String name, G4double* a, G4double* z,
              G4double dens, G4int nlmat, G4double* wmat)
{
// from fortran (g3routines.F)
// write(fmt,'(A,I3,A,I3,A,I3,A)')
//+      '(a4,1x,a6,i5,1x,''"'',a,''"'',1x,e16.8,1x,i3,',
//+      max(nlmata,1),
//+      '(1x,e16.8),',max(nlmata,1),'(1x,e16.8),',
//+      max(nlmata,1),'(1x,e16.8))'
// write(lunlist,fmt)
//+      context, rname, imate, name, dens,
//+      nlmat, 
//+      (a(k), k=1,abs(nlmat)),
//+      (z(k), k=1,abs(nlmat)),
//+      (wmat(k), k=1,abs(nlmat))
// ---

  G4String context("----");
  G4String rname("GSMIXT");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << setw(5) << imate << space 
	  << '"' << name << '"' << space
          << setw(16) << setprecision(8) << dens << space
	  << setw(3) << nlmat;
  G4int i;	  
  for (i=0; i<abs(nlmat); i++)	      
    fOutFile << space << setw(16) << setprecision(8) << a[i];
  for (i=0; i<abs(nlmat); i++)	      
    fOutFile << space << setw(16) << setprecision(8) << z[i];
  for (i=0; i<abs(nlmat); i++)	      
    fOutFile << space << setw(16) << setprecision(8) << wmat[i];
  fOutFile << endl;
}	      


void TG4GeometryOutputManager::WriteGstmed(
              G4int itmed, G4String name, G4int nmat, G4int isvol,
              G4int ifield, G4double fieldm, G4double tmaxfd,
              G4double stemax, G4double deemax, G4double epsil,
              G4double stmin, G4double* ubuf, G4int nwbuf)
{
// from fortran (g3routines.F)
// write(fmt,'(A,I3,A)') 
//>      '(a4,1x,a6,i5,1x,''"'',a,''"'',3i3,6(1x,e16.8),i3,',
//>      max(nwbuf,1),'(1x,e16.8))'
// write(lunlist,fmt)
//+      context, rname, itmed, name, nmat, isvol, ifield, fieldm,
//+      tmaxfd, stemax, deemax, epsil, stmin,
//+      nwbuf, (ubuf(k),k=1,nwbuf)
// ---

  G4String context("----");
  G4String rname("GSTMED");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << setw(5) << itmed << space
	  << '"' << name << '"' 
	  << setw(3) << nmat
	  << setw(3) << isvol
	  << setw(3) << ifield << space
          << setw(16) << setprecision(8) << fieldm << space
          << setw(16) << setprecision(8) << tmaxfd << space
          << setw(16) << setprecision(8) << stemax << space
          << setw(16) << setprecision(8) << deemax << space
          << setw(16) << setprecision(8) << epsil << space
          << setw(16) << setprecision(8) << stmin << space
	  << setw(3) << nwbuf;
  for (G4int i=0; i<nwbuf; i++)	      
    fOutFile << space << setw(16) << setprecision(8) << ubuf[i];
  fOutFile << endl;
}	  


void TG4GeometryOutputManager::WriteGstpar(
               G4int itmed, G4String param, G4double parval)
{	       
// from fortran (g3routines.F)
// write(lunlist,
//+     '(a4,1x,a6,i5,1x,a4,(1x,e16.8))')
//+     context, rname, itmed, chpar, parval
// ---

  G4String context("----");
  G4String rname("GSTPAR");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   << space
	  << setw(5) << itmed << space 
	  << param   << space
          << setw(16) << setprecision(8) << parval
          << endl;
}	      


void TG4GeometryOutputManager::WriteGgclos()
{
// Writes GGCLOS token
// ---

  G4String context("----");
  G4String rname("GGCLOS");
  G4String space(" "); 
  fOutFile << context << space 
          << rname   
	  << endl;
}	  
