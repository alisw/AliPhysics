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

/* $Id$*/

#include "TFlukaScoringOption.h"
#include "TFlukaMCGeometry.h"
#include <TGeoManager.h>
#include <TGeoVolume.h>

ClassImp(TFlukaScoringOption)

FILE*              TFlukaScoringOption::fgFile(0x0);
TFlukaMCGeometry*  TFlukaScoringOption::fgGeom(0x0);

TFlukaScoringOption::TFlukaScoringOption()
   : fNopfl(0),
     fOutFile(""),
     fLun(0)
{

    // Default constructor
}

TFlukaScoringOption::TFlukaScoringOption(const char* name, const char* sdum, Int_t nopfl, char* outfile, Float_t* what)
    : TNamed(name, sdum),
      fNopfl(nopfl),
      fOutFile(outfile),
      fLun(0)
{
    // Constructor
    fNopfl   = nopfl;
    fOutFile = outfile;
    Int_t npar = (nopfl == 0)? 6 : 12;
    for (Int_t i = 0; i < npar; i++)  fWhat[i] = what[i];
}

TFlukaScoringOption::TFlukaScoringOption(const char* name, const char* sdum, Int_t nopfl, char* outfile, Float_t* what, 
                                         const char* det1, const char* det2, const char* det3)
    : TNamed(name, sdum),
      fNopfl(nopfl + 2),
      fOutFile(outfile),
      fLun(0)
{
    // Constructor
//    fNopfl   = nopfl + 2;
//    fOutFile = outfile;
    Int_t npar = (nopfl == 0)? 6 : 12;
    for (Int_t i = 0; i < npar; i++)  fWhat[i] = what[i];

    fName[0] = det1;
    fName[1] = det2;
    fName[2] = det3;
}


//--- GET METHODS
//______________________________________________
const char* TFlukaScoringOption::GetRegName(Int_t ndet) 
{
    // Return ndet'th region name
    return fName[ndet - 1];
}

//
// Write Fluka Input Cards
//
void TFlukaScoringOption::WriteOpenFlukaFile()
{
    //
    // Write Fluka input card for output file opening
    // 
    fprintf(fgFile, "OPEN      %10.1f                                                    %s\n", 
            GetLun(), "NEW");
    fprintf(fgFile, "%s\n", GetFileName());
}


void TFlukaScoringOption::WriteFlukaInputCards()
{
    //
    // Write the Fluka Input Cards for this scoring option
    //

    const char* cont_line = "&";

//***************************************************************************
//*    Par()==0  only 6 parameter && regions by float                       *
//*    Par()==1  only 6 parameter && regions by Name                        *
//*    Par()==2      12 parameter && regions by float                       *
//*    Par()==3      12 parameter && regions by Name                        *
//*    Par()==4      12 parameter && regions by Name (for USRBIN only)      *
//***************************************************************************
//
// USRBIN
//
    if(strncmp(GetName(), "USRBIN", 6) == 0){
        if (Par() == 0) {
            fprintf(fgFile, "USRBIN    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%s\n",
                    What(1), What(2), GetLun(), What(4), What(5), What(6), GetTitle());
        } else if (Par() == 1) {
            fprintf(fgFile, "USRBIN    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%s\n",
                    What(1), What(2), GetLun(), What(4), What(5),  What(6), GetTitle());
            fprintf(fgFile, "USRBIN    %10.1f%10.4g%10.1f%10.1f%10.1f%10.1f  %s\n",
                    What(7), What(8), What(9), What(10), What(11), What(12), cont_line);
        } else if (Par() == 2) {
            if(What(1) == 2.0 || What(1) == 12){
                fprintf(fgFile, "USRBIN    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%s\n",
                        What(1), What(2), GetLun(), Float_t(GetRegionByName(GetRegName(1))),
                        Float_t(GetRegionByName(GetRegName(2))), Float_t(GetRegionByName(GetRegName(3))), GetTitle());
                fprintf(fgFile, "USRBIN    %10.1f%10.4g%10.1f%10.1f%10.1f%10.1f  %s\n",
                        Float_t(GetRegionByName(GetRegName(1))), Float_t(GetRegionByName(GetRegName(2))),
                        Float_t(GetRegionByName(GetRegName(3))), 1., 1., 1., cont_line);
            } else {
                printf("Check consistency of SetUserScoring values \n");
            }
        }
    }
   
//
// USRBDX
//
    if(strncmp(GetName(), "USRBDX", 6) == 0){
        if (Par() == 2) {
            fprintf(fgFile, "USRBDX    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n",
                    What(1), What(2), GetLun(), Float_t(GetRegionByName(GetRegName(1))),
                    Float_t(GetRegionByName(GetRegName(2))), What(6));
        } else if (Par() == 3) {
            fprintf(fgFile, "USRBDX    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n",
                    What(1), What(2), GetLun(), Float_t(GetRegionByName(GetRegName(1))),
                    Float_t(GetRegionByName(GetRegName(2))), What(6));
            fprintf(fgFile, "USRBDX    %10.1f%10.4g%10.1f%10.1f%10.1f%10.1f  %s\n",
                    What(7), What(8), What(9), What(10), What(11), What(12), cont_line);
        }
    }
    
//
// USTRACK
//
    if(strncmp(GetName(), "USRTRACK", 6) == 0){
        fprintf(fgFile, "USRTRACK  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n",
                What(1), What(2), GetLun(), Float_t(GetRegionByName(GetRegName(1))), What(4), What(5));
        fprintf(fgFile, "USRTRACK  %10.1f%10.4g                                          %s\n",
                What(7), What(8), cont_line);
    }
    
//
// USRCOLL
// 
    if(strncmp(GetName(), "USRCOLL", 6) == 0){
        fprintf(fgFile, "USRCOLL   %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n",
                What(1), What(2), GetLun(), Float_t(GetRegionByName(GetRegName(1))), What(4), What(5));
        fprintf(fgFile, "USRCOLL   %10.1f%10.4g                                          %s\n",
                What(7), What(8), cont_line);
    }
}


Int_t TFlukaScoringOption::GetRegionByName(const char* detname)
{
// Get number of region for a given detector name
    TGeoVolume* vol = dynamic_cast<TGeoVolume*>((gGeoManager->GetListOfVolumes())->FindObject(detname));
    Int_t ireg = (vol)? vol->GetNumber() : -999;
    return ireg;
}         
