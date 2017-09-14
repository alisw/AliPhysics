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

#include "AliJRunTable.h"

#include <iostream>
#include <cstdlib>
#include <TSystem.h>
#include <TPRegexp.h>
using namespace std;

AliJRunTable::AliJRunTable() :
  fCPeriod(0),
  fCRunNumber(0),
  fCPeriodMCName("undefined")
{
  Init();
}

void AliJRunTable::Init(){
  // comment needed
    SetPeriodInformation( kUnknownPeriod, "UnknownPeriod", kPP, kRE, 0, -1, -1, "UnKnownPeriod" );
    SetPeriodInformation( kLHC10b, "LHC10b", kPP, kRE, 7000, 114737, 117223, "LHC10d1" );
    SetPeriodInformation( kLHC10c, "LHC10c", kPP, kRE, 7000, 118503, 121040, "LHC10d4" );
    SetPeriodInformation( kLHC10d, "LHC10d", kPP, kRE, 7000, 122195, 126437, "LHC10f6a" );
    SetPeriodInformation( kLHC10e, "LHC10e", kPP, kRE, 7000, 127712, 130850, "LHC10e20" );
    SetPeriodInformation( kLHC10h, "LHC10h", kPbPb, kRE, 2760, 136851, 139517, "LHC11a10a_bis" );
    SetPeriodInformation( kLHC11h, "LHC11h", kPbPb, kRE, 2760, 167813, 170595, "LHC12a17" );
    SetPeriodInformation( kLHC15o, "LHC15o", kPbPb, kRE, 5020, 244640, 247173, "LHC16g" ); // HIJING MB 16g1

    //LHC11a
    SetPeriodInformation(kLHC11a, "LHC11a", kPP, kRE, 2760, 144871, 146860, "LHC11b10a" );

    // pp 7TeV LHC11bcde
    SetPeriodInformation( kLHC11b, "LHC11b", kPP, kRE, 7000, 0, 0, "LHC12d2_plus" );
    SetPeriodInformation( kLHC11c, "LHC11c", kPP, kRE, 7000, 153533, 154789, "LHC12d2_plus" );
    SetPeriodInformation( kLHC11d, "LHC11d", kPP, kRE, 7000, 156620, 159580, "LHC12d2_plus" );
    SetPeriodInformation( kLHC11e, "LHC11e", kPP, kRE, 7000, 0, 0, "LHC12d2_plus" );

    // LHC12g - TODO
    SetPeriodInformation( kLHC12g, "LHC12g",kPA, kRE, 5020, 188356,188503, "LHC13b2" );
    SetPeriodInformation( kLHC12h, "LHC12h",kPA, kRE, 5020, 189122,192732, "LHC13b2" );
    SetPeriodInformation( kLHC13b, "LHC13b",kPA, kRE, 5020, 195344,195483, "LHC13b2-efix_p1" );
    SetPeriodInformation( kLHC13c, "LHC13c",kPA, kRE, 5020, 195529,195677, "LHC13b2-efix_p1" );
    SetPeriodInformation( kLHC13d, "LHC13d",kPA, kRE, 5020, 195724,195872, "LHC13b2-efix_p1" );
    SetPeriodInformation( kLHC13e, "LHC13e",kPA, kRE, 5020, 195955,196310, "LHC13b2-efix_p1" );
    SetPeriodInformation( kLHC13g, "LHC13g",kPA, kRE, 5020, 197669,200000, "LHC13b2-efix_p1" );
    // p-Pb 5TeV
    SetPeriodInformation( kLHC16q, "LHC16q",kPA, kRE, 5020, 264896,265533, "LHC17f2" );
}


TString AliJRunTable::GetPeriodName( int period ) const {
    // TODO 
    if( period < 0 ) period = fCPeriod;
    return fPeriodName[period];
}

int AliJRunTable::GetRunNumberToPeriod( int runnumber ){
  // comment needed
    int period = -1;
    for( int ip=0;ip<kJNPeriod;ip++ ){
        if(fDataType[ip] == kMC ) continue;
        if( runnumber >= fRunRange[ip][0] && runnumber <= fRunRange[ip][1] ){
            //cout<< fPeriodName[ip] <<"\t"<<fRunRange[ip][0]<<"~"<<fRunRange[ip][1]<<endl;
            period = ip;
            break;
        }
    }
    if( period < 0 ){ 
        cout<<"J_ERROR : no period for runnumber "<<runnumber<<endl;
        gSystem->Exit(1);
    }
    return period;
}

void AliJRunTable::SetPeriodInformation(int period, TString name, int beamtype, int datatype, int energy, int run0, int run1, TString MCPeriod){
  // comment needed
    fPeriodName[period] = name;
    fBeamType[period] = beamtype;
    fDataType[period] = datatype;
    fEnergy[period] = energy;
    fRunRange[period][0] = run0;
    fRunRange[period][1] = run1;
    fMCPeriod[period] = MCPeriod;
}


int AliJRunTable::GetPeriodCode( TString perstr ) const{
  // comment needed
    int period = kUnknownPeriod;
    for( int ip=0;ip<kJNPeriod;ip++ ){
        if( perstr == GetPeriodName(ip) ){
            period = ip;
            break;
        }
    }
    if( period <0 ){
        cout<<"J_ERROR : no period for "<<perstr<<endl;
    }
    return period;

}

int AliJRunTable::GetRunNumberFromString(const char * tstr ){
  // comment needed
    TPMERegexp rexRunNumber( "1\\d{5}" );
    rexRunNumber.Match(tstr);
    return TString(rexRunNumber[0]).Atoi();
}

TString AliJRunTable::GetPeriodFromString(const char * tstr ) const{
  // comment needed
    TPMERegexp rexPeriod( "LHC1\\d[a-zA-Z]" );
    rexPeriod.Match(tstr);
    return rexPeriod[0];
}

TString AliJRunTable::GetMCPeriodFromString(const char * tstr ) const{
  // comment needed
    TPMERegexp rexPeriod( "LHC1\\d[a-zA-Z0-9]{2,}(\\w+)?" );
    rexPeriod.Match(tstr);
    return rexPeriod[0];
}

int AliJRunTable::GetPeriodFromMCPeriod( const char * tstr ){
  // comment needed
    int period = -1;
    for( int ip=0;ip<kJNPeriod;ip++ ){
        if(fDataType[ip] == kMC ) continue;
        if( fMCPeriod[ip] == TString(tstr) ){
            period = ip;
            break;
        }
    }
    if( period < 0 ){ 
        cout<<"J_ERROR : no period for MCPeriod "<<tstr<<endl;
        exit(1);
    }
    return period;
}

bool AliJRunTable::ParseString( const char * tstr ){
  // comment needed
    fCRunNumber = GetRunNumberFromString(tstr);
    fCPeriodMCName = GetMCPeriodFromString(tstr);
    if( fCRunNumber > 0 ){
        fCPeriod = GetRunNumberToPeriod( fCRunNumber );
    }else
    if( fCPeriodMCName.Length() > 0 ){
        fCPeriod = GetPeriodFromMCPeriod( fCPeriodMCName );
    }else{
        fCPeriod = kUnknownPeriod;
    }
    return true;
}

// GetBeamStr is never used anywhere. Should it be removed altogether?
const char * AliJRunTable::GetBeamStr( int ib ) const {
  // comment needed
    if( ib < 0 ) ib = fBeamType[fCPeriod]; 
    switch (ib){
        case kPP : return "pp";
        case kPbPb: return "PbPb";
        case kPA:   return "PA";
    }
    return NULL;
}

const char * AliJRunTable::GetDataStr( int ib ) const {
  // comment needed
    if( ib < 0 ) ib = fDataType[fCPeriod]; 
    switch (ib){
        case kRE: return "REAL";
        case kMC: return "MC";
    }
    return NULL;
}
AliJRunTable& AliJRunTable::GetSpecialInstance(){
  // comment needed
    static AliJRunTable instance;
    return instance;
}

// GetEnergyStr is never used anywhere. Should it be removed altogether?
const char * AliJRunTable::GetEnergyStr( int ib ) const {
  // comment needed
    if( ib < 0 ) ib = fEnergy[fCPeriod]; 
    if( ib < 1000 ) return Form("%dGeV",ib);
    if( ib == 2760 ) return "2.76TeV";
    return Form("%dTeV",ib/1000);
}

const AliJRunTable& AliJRunTable::GetInstance(){
  // comment needed
    return GetSpecialInstance();
}
