/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

// Comment describing what this class does needed!

#include "AliJTrackCounter.h"

// track counter
// blah
// blah
// blah
// blah

ClassImp(AliJTrackCounter)

        AliJTrackCounter::AliJTrackCounter() :
        findex(-1),
        fptBin(-1),
        fpt(-1),
        fdphi(-1)
        {
          Reset();
        }


        AliJTrackCounter::AliJTrackCounter(int inind0, int inind1, int inlist0, int inlist1, double indphi) :
        findex(-1),
        fptBin(-1),
        fpt(-1),
        fdphi(indphi)
        {
          // constructor
            Reset();
            fpairTrackID[0]=inind0; 
            fpairTrackID[1]=inind1;
            fpairListID[0]=inlist0;
            fpairListID[1]=inlist1;
            fdphi = indphi;
        }    //constructor



        void AliJTrackCounter::Reset(){
          // reset
            findex=-1; 
            fpairTrackID[0]=-1; fpairTrackID[1]=-1;  
            fpairListID[0]=-1;  fpairListID[1]=-1;  
            fpt=-1; fptBin=-1; fdphi=-1;
        }



       AliJTrackCounter& AliJTrackCounter::operator=(const AliJTrackCounter& counter){
         // copy constructor
           fpairTrackID[0] = counter.fpairTrackID[0];
           fpairTrackID[1] = counter.fpairTrackID[1];
           fpairListID[0] = counter.fpairListID[0];
           fpairListID[1] = counter.fpairListID[1];
           fdphi           = counter.fdphi;
           return *this; 
       }















