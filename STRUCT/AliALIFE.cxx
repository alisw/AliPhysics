/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.1  2000/07/26 15:10:57  morsch
Helper class to write geometry in ALIFE format in parallel with Geant geometry definition.

*/

#include <AliALIFE.h>

ClassImp(AliALIFE)

    AliALIFE::AliALIFE(const char* file1, const char* file2)
{  
// Constructor
    fNBodies = 0; 
    fNVolumes= 0; 
    fBodyFile   = file1;   // File for Fluka bodies
    fVolumeFile = file2;
    fFile1=fopen(fBodyFile,"w");
    fFile2=fopen(fVolumeFile,"w");
    BodyHeader();
    VolumeHeader();
    fDefaultVolume1 = "DEFAU1";
    fDefaultVolume2 = "DEFAU2";    
}

    AliALIFE::AliALIFE()
{
// Default constructor
    fBodyFile   = "FlukaBody.inp";  
    fVolumeFile = "FlukaVolume.inp";
    fFile1=fopen(fBodyFile,"w");
    fFile2=fopen(fVolumeFile,"w");
    BodyHeader();
    VolumeHeader();
}

void AliALIFE::BodyHeader()
{
// Write header for body definitions
    fprintf(fFile1,"*\n");        
    fprintf(fFile1,"GEOBEGIN                                                              COMBINAT\n");
    fprintf(fFile1,"    0    0         AliRoot Generated\n");
    fprintf(fFile1,"*\n");
    fprintf(fFile1,"*== Body Definitions ================================================\n");
    fprintf(fFile1,"*\n");
                            
}


void AliALIFE::VolumeHeader()
{
// Write header for region (volume)  definitions
    fprintf(fFile2,"REGIONS\n");
    fprintf(fFile2,"*\n");
    fprintf(fFile2,"*== Region Definitions =================================================\n");
    fprintf(fFile2,"*\n");
}


void AliALIFE:: Cylinder(Float_t rmin, Float_t rmax, 
			 Float_t zmin, Float_t zmax, 
			 Float_t pos[3], 
			 char* Material, char* Field, char* Cuts) 
{
// Simple cylinder
//
// Bodies
// ^^^^^^
    char name1[5], name2[5], name3[5], name4[5];

//  outer radius
    sprintf(name1, "R%4.4d", fNBodies++);
    fprintf(fFile1,"%5s ZCC%10.3f%10.3f%10.3f\n",name1, pos[0], pos[1], rmax); // inner radius
    sprintf(name2, "R%4.4d", fNBodies++);
    fprintf(fFile1,"%5s ZCC%10.3f%10.3f%10.3f\n",name2, pos[0], pos[1], rmin); 
// z-max
    sprintf(name3, "Z%4.4d", fNBodies++);
    fprintf(fFile1,"%5s XYP%10.3f\n",name3, zmax); 
// z-min
    sprintf(name4, "Z%4.4d", fNBodies++);
    fprintf(fFile1,"%5s XYP%10.3f\n",name4, zmin); 
//
// Volumes
// ^^^^^^^

    fprintf(fFile2,">%s:%s\n", Material, Field);
    fprintf(fFile2,"EMFCUT=%s\n", Cuts);
    fprintf(fFile2,"WW-FACTOR=%s\n", Cuts);
    if (rmin >0) {
	fprintf(fFile2,"+%s-%s+%s-%s\n", name1, name2, name3, name4);
    } else {
	fprintf(fFile2,"+%s+%s-%s\n", name1, name3, name4);
    }
    
    fprintf(fFile2,"\n");    
}

void AliALIFE::OnionCylinder(Float_t* r, Int_t nr, Float_t zmin, Float_t zmax,
			     Float_t pos[3], char** Materials, 
			     char** Fields, char** Cuts) 
{
//
// Concentric cylinders
//

// Bodies
// ^^^^^^
    char nameRin[6], nameRou[6], nameZin[6], nameZou[6];
// z-limits
// z-max
    sprintf(nameZou, "Z%4.4d", fNBodies++);
    nameZou[5]='\0';
    
    fprintf(fFile1,"%5s XYP%10.3f\n",nameZou, zmax); 
// z-min
    sprintf(nameZin, "Z%4.4d", fNBodies++);
    nameZin[5]='\0';
    fprintf(fFile1,"%5s XYP%10.3f\n",nameZin, zmin); 

// inner radius
    sprintf(nameRin, "R%4.4d", fNBodies++);
    nameRin[5]='\0';
    fprintf(fFile1,"%5s ZCC%10.3f%10.3f%10.3f\n",nameRin, pos[0], pos[1], r[0]); 

    Int_t i;
    for (i=1; i<nr; i++) {
//  outer radius
	sprintf(nameRou, "R%4.4d", fNBodies++);
	nameRou[5]='\0';
	fprintf(fFile1,"%5s ZCC%10.3f%10.3f%10.3f\n",
		nameRou, pos[0], pos[1], r[i]); 	
//
// Volumes
// ^^^^^^^
	if (Fields && Cuts) {
	    fprintf(fFile2,">%s:%s\n", Materials[i-1], Fields[i-1]);
	    fprintf(fFile2,"EMFCUT=%s\n", Cuts[i-1]);
	    fprintf(fFile2,"WW-FACTOR=%s\n", Cuts[i-1]);
	} else {
	    fprintf(fFile2,">%s:%s\n", Materials[i-1], "MF");
	    fprintf(fFile2,"EMFCUT=%s\n", "$UNSHIELDED");
	    fprintf(fFile2,"WW-FACTOR=%s\n", "$UNSHIELDED");
	}
	if (r[i-1] != 0.) {
	    fprintf(fFile2,"+%5s-%5s+%5s-%5s\n", nameZou, nameZin, nameRou, nameRin);
	} else {
	    fprintf(fFile2,"+%5s-%5s+%5s\n", nameZou, nameZin, nameRou);
	}
	fprintf(fFile2,"\n");
	strcpy(nameRin,nameRou);
    }
}


void AliALIFE::Cone(Float_t rmin1, Float_t rmin2, 
		     Float_t rmax1, Float_t rmax2,
		     Float_t zmin, Float_t zmax, 
		     Float_t pos[3], 
		     char* Material, char* Field, char* Cuts) 
{
// Simple cone 

//
// Bodies
// ^^^^^^
    char nameCou[6], nameCin[6];
    Float_t d, r1, r2;
    char nameZin[6], nameZou[6];
// z-max
    sprintf(nameZou, "Z%4.4d", fNBodies++);
    fprintf(fFile1,"%5s XYP%10.3f\n",nameZou, zmax); 
// z-min
    sprintf(nameZin, "Z%4.4d", fNBodies++);
    fprintf(fFile1,"%5s XYP%10.3f\n",nameZin, zmin); 
    
//  outer radius
    d=zmax-zmin;
    if (rmax1 >= 0. && rmax2 >= 0.) {
	if (rmax1!=rmax2) {
	    if (rmax1 > rmax2) {
		pos[2]=zmin;
		r1=rmax1;
		r2=rmax2;
	    } else {
		d=-d;
		pos[2]=zmax;
		r1=rmax2;
		r2=rmax1;
	    }
	    sprintf(nameCou, "C%4.4d", fNBodies++);
	    fprintf(fFile1,"%5s TRC%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n",
		    nameCou, pos[0], pos[1], pos[2], 0., 0., d); 
	    fprintf(fFile1,"         %10.3f%10.3f\n",r1,r2);
	} else {
	    sprintf(nameCou, "C%4.4d", fNBodies++);
	    fprintf(fFile1,"%5s ZCC%10.3f%10.3f%10.3f\n",
		    nameCou, pos[0], pos[1], rmax1); 
	} 
    }else {
	strcpy(nameCou,fDefaultVolume1);
    }

    
    
// inner radius
    if (rmin1 >= 0. && rmin2 >= 0.) {
	if (rmin1!=rmin2) {
	    if (rmin1 != 0 && rmin2 !=0) {
		d=zmax-zmin;
		if (rmin1 > rmin2) {
		    pos[2]=zmin;
		    r1=rmin1;
		    r2=rmin2;
		} else {
		    pos[2]=zmax;
		    r1=rmin2;
		    r2=rmin1;
		    d=-d;
		}
		sprintf(nameCin, "C%4.4d", fNBodies++);
		fprintf(fFile1,"%5s TRC%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n",
			nameCin, pos[0], pos[1], pos[2], 0., 0., d); 
		fprintf(fFile1,"         %10.3f%10.3f\n",r1,r2);
	    } 
	} else {
	    sprintf(nameCin, "C%4.4d", fNBodies++);
	    fprintf(fFile1,"%5s ZCC%10.3f%10.3f%10.3f\n",
		    nameCin, pos[0], pos[1], rmin1); 
	} 
    }else {
	strcpy(nameCin,fDefaultVolume2);
    }



//
// Volumes
// ^^^^^^^
    fprintf(fFile2,">%s:%s\n", Material, Field);
    fprintf(fFile2,"EMFCUT=%s\n", Cuts);
    fprintf(fFile2,"WW-FACTOR=%s\n", Cuts);
    if (rmin1 != 0 && rmin2 !=0) {
	fprintf(fFile2,"+%s-%s+%s-%s\n", nameCou, nameCin, nameZou, nameZin);
    } else {
	fprintf(fFile2,"+%s+%s-%s\n", nameCou, nameZou, nameZin);
    }
    fprintf(fFile2,"\n");        
}

void AliALIFE::OnionCone (Float_t* r1, Float_t* r2, Int_t nr, 
			  Float_t zmin, Float_t zmax,
			  Float_t pos[3], char** Materials, char** Fields,
			  char** Cuts) 
{
// Concentric cones
//
    char nameCin[6], nameCou[6], nameZin[6], nameZou[6];
//
// Bodies
// ^^^^^^
    Float_t ri1, ri2, ro1, ro2, d;

// z-max
    sprintf(nameZou, "Z%4.4d", fNBodies++);
    fprintf(fFile1,"%5s XYP%10.3f\n",nameZou, zmax); 
// z-min
    sprintf(nameZin, "Z%4.4d", fNBodies++);
    fprintf(fFile1,"%5s XYP%10.3f\n",nameZin, zmin); 
     
//  inner radius
    d=zmax-zmin;
    Bool_t hasInner=kFALSE;
    Bool_t isCylinder=kFALSE;
    if (r1[0]>=0 && r2[0]>=0) {
	if (r1[0] != 0. &&  r2[0] !=0.) {
	    if (r1[0]!=r2[0]) {
		hasInner=kTRUE;
		if (r1[0] > r2[0]) {
		    pos[2]=zmin;
		    ri1=r1[0];
		    ri2=r2[0];
		} else {
		    d=-d;
		    pos[2]=zmax;
		    ri1=r2[0];
		    ri2=r1[0];
		}
		sprintf(nameCin, "C%4.4d", fNBodies++);
		nameCin[5]='\0';
		fprintf(fFile1,"%5s TRC%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n",
			nameCin, pos[0], pos[1], pos[2], 0., 0., d); 
		fprintf(fFile1,"         %10.3f%10.3f\n",ri1,ri2);
	    } else {
		isCylinder=kTRUE;
		sprintf(nameCin, "C%4.4d", fNBodies++);
		nameCin[5]='\0';
		fprintf(fFile1,"%5s ZCC%10.3f%10.3f%10.3f\n",
			nameCin, pos[0], pos[1], r1[0]); 
	    }
	}
    } else {
	strcpy(nameCin,fDefaultVolume1);
    }
    
    
    
// outer radius
    Int_t i;
    for (i=1; i<nr; i++) {
	if (r1[i] >= 0. && r2[i] >=0) {
	    if (r1[i]!=r2[i]) {
		d=zmax-zmin;
		if (r1[i] > r2[i]) {
		    pos[2]=zmin;
		    ro1=r1[i];
		    ro2=r2[i];
		} else {
		    pos[2]=zmax;
		    ro1=r2[i];
		    ro2=r1[i];
		    d=-d;
		}
		sprintf(nameCou, "C%4.4d", fNBodies++);
		nameCou[5]='\0';
		fprintf(fFile1,"%5s TRC%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n",
			nameCou, pos[0], pos[1], pos[2], 0., 0., d); 
		fprintf(fFile1,"         %10.3f%10.3f\n",ro1,ro2);
	    } else {
		isCylinder=kTRUE;
		sprintf(nameCou, "C%4.4d", fNBodies++);
		nameCou[5]='\0';
		fprintf(fFile1,"%5s ZCC%10.3f%10.3f%10.3f\n",
			nameCou, pos[0], pos[1], r1[i]); 
	    }
	} else {
	    strcpy(nameCou,fDefaultVolume1);
	}
	
// Volumes
// ^^^^^^^
	if (Fields && Cuts) {
	    fprintf(fFile2,">%s:%s\n", Materials[i-1], Fields[i-1]);
	    fprintf(fFile2,"EMFCUT=%s\n", Cuts[i-1]);
	    fprintf(fFile2,"WW-FACTOR=%s\n", Cuts[i-1]);
	} else {
	    fprintf(fFile2,">%s:%s\n", Materials[i-1], "MF");
	    fprintf(fFile2,"EMFCUT=%s\n", "$UNSHIELDED");
	    fprintf(fFile2,"WW-FACTOR=%s\n", "$UNSHIELDED");
	}
	if (hasInner) {
	    fprintf(fFile2,"+%5s-%5s+%5s-%5s\n", 
		    nameCou, nameCin, nameZou, nameZin);
	} else {
	    fprintf(fFile2,"+%5s\n", nameCou);
	    hasInner=kTRUE;
	}
	fprintf(fFile2,"\n");        
	strcpy(nameCin,nameCou);
    }
}


void AliALIFE::PolyCone(Float_t* rmin, Float_t* rmax, Float_t* z, 
			Int_t nz,
			Float_t pos[3], 
			char* Material, char* Field, char* Cuts) 
{
//
// Equivalent to the Geant3 PCON
    Int_t i;
    for (i=0; i<nz-1; i++) {
// skip (z_i = z_i+1)
	if (z[i] == z[i+1]) continue;
	Cone(rmin[i], rmin[i+1], rmax[i], rmax[i+1], z[i], z[i+1], pos,
	     Material, Field, Cuts);
    }
}

void AliALIFE::OnionPolyCone(Float_t** r, Float_t* z,
			     Int_t nr, Int_t nz,
			     Float_t pos[3], 
			     char** Materials, char** Fields, char** Cuts)
{
//
// Concentric PCONS
    Int_t i, j;
    for (i=0; i<nz-1; i++) {
// skip (z_i = z_i+1)
	if (z[i] == z[i+1]) continue;
	for (j=0; j<nr-1; j++) {
	    if (Fields && Cuts) {
		Cone(r[i][j], r[i+1][j], r[i][j+1], r[i+1][j+1], 
		     z[i], z[i+1], pos,
		     Materials[i], Fields[i], Cuts[i]);
	    } else {
		Cone(r[i][j], r[i+1][j], r[i][j+1], r[i+1][j+1], 
		     z[i], z[i+1], pos,
		     Materials[i]);
	    }
	}
    }
}

void AliALIFE::Comment(char* Comment)
{
// Insert comment line
    fprintf(fFile1,"*%s\n", Comment);        
    fprintf(fFile2,"*%s\n", Comment);        
}


void AliALIFE::Finish()
{
// Finish geometry definition
    char s[BUFSIZ];
    fprintf(fFile1,"      END\n");
    fclose(fFile2);
    fFile2=fopen(fVolumeFile,"r");
    while (fgets(s, BUFSIZ, fFile2)) {
	fputs(s,fFile1);
    }
    
    fclose(fFile1);
    fclose(fFile2);    
}











