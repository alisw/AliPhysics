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
$Id$
*/
#include <Riostream.h>
#include <TRandom.h>
#include <TH1.h>
#include <TMath.h>
#include <TString.h>
#include <TParticle.h>

#include "AliRun.h"
#include "AliITS.h"
#include "AliITShit.h"
#include "AliITSdigitSPD.h"
#include "AliITSmodule.h"
#include "AliITSMapA2.h" 
#include "AliITSpList.h"
#include "AliITSsimulationSPDdubna.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSresponseSPD.h"

//#define DEBUG

ClassImp(AliITSsimulationSPDdubna)
////////////////////////////////////////////////////////////////////////
// Version: 1
// Modified by Bjorn S. Nilsen
// Version: 0
// Written by Boris Batyunya
// December 20 1999
//
// AliITSsimulationSPDdubna is to do the simulation of SPDs.
//______________________________________________________________________
AliITSsimulationSPDdubna::AliITSsimulationSPDdubna():
AliITSsimulation(),
fHis(0),
fSPDname(),
fCoupling(0){
    // Default constructor.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A default constructed AliITSsimulationSPDdubna class.

    if(GetDebug(1)) Info("AliITSsimulationSPDdubda()",
                         "Calling degault constructor");
}
//______________________________________________________________________
AliITSsimulationSPDdubna::AliITSsimulationSPDdubna(AliITSsegmentation *seg,
						                 AliITSresponse *resp,Int_t cup):
AliITSsimulation(seg,resp),
fHis(0),
fSPDname(),
fCoupling(cup){
    // standard constructor
    // Inputs:
    //    AliITSsegmentation *seg  A pointer to the segmentation class
    //                             to be used for this simulation
    //    AliITSresponse     *resp A pointer to the responce class to
    //                             be used for this simulation
    //    Int_t              cup   The type of coupling to be used
    //                             =1 uses SetCoupling, =2 uses SetCouplingOld
    //                             With diffusion tured off
    //                             =3 uses SetCoupling, =4 uses SetCouplingOld
    //                             with diffusion on other, no coupling.
    // Outputs:
    //    none.
    // Return:
    //    A default constructed AliITSsimulationSPDdubna class.

    if(GetDebug(1)) Info("AliITSsimulationSPDdubda",
                         "Calling degault constructor seg=%p resp=%p cup=%d",
                         seg,resp,cup);
    if(cup==1||cup==2){ // For the moment, remove defusion if Coupling is
        // set.
        resp->SetTemperature(0.0);
        resp->SetDistanceOverVoltage(0.0);
    } // end if
    Init();
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::Init(){
    // Initilization
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    const Double_t kmictocm = 1.0e-4; // convert microns to cm.

    SetModuleNumber(0);
    SetEventNumber(0);
    SetMap(new AliITSpList(GetNPixelsZ(),GetNPixelsX()));
    GetResp(0,0)->SetDistanceOverVoltage(kmictocm*GetSeg()->Dy(),50.0);
}
//______________________________________________________________________
AliITSsimulationSPDdubna::~AliITSsimulationSPDdubna(){
    // destructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //     none.

    if (fHis) {
        fHis->Delete(); 
        delete fHis;     
    } // end if fHis
}
//______________________________________________________________________
AliITSsimulationSPDdubna::AliITSsimulationSPDdubna(const 
						   AliITSsimulationSPDdubna 
						   &s) : AliITSsimulation(s){
    //     Copy Constructor
    // Inputs:
    //    AliITSsimulationSPDdubna &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:

    *this = s;
    return;
}
//______________________________________________________________________
AliITSsimulationSPDdubna&  AliITSsimulationSPDdubna::operator=(const 
                                           AliITSsimulationSPDdubna &s){
    //    Assignment operator
    // Inputs:
    //    AliITSsimulationSPDdubna &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:

    if(&s == this) return *this;
    this->fHis = s.fHis;
    fCoupling  = s.fCoupling;
    fSPDname   = s.fSPDname;
    return *this;
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::InitSimulationModule(Int_t module, Int_t event){
    //  This function creates maps to build the list of tracks for each
    //  summable digit. Inputs defined by base class.
    //  Inputs:
    //    Int_t module   // Module number to be simulated
    //    Int_t event    // Event number to be simulated
    //  Outputs:
    //    none
    //  Returns:
    //    none

    if(GetDebug(1)) Info("InitSimulationModule","(module=%d,event=%d)",
                         module,event);
    SetModuleNumber(module);
    SetEventNumber(event);
    ClearMap();
}
//_____________________________________________________________________
void AliITSsimulationSPDdubna::SDigitiseModule(AliITSmodule *mod,Int_t,
                                               Int_t event){
    //  This function begins the work of creating S-Digits.  Inputs defined
    //  by base class.
    //  Inputs:
    //    AliITSmodule *mod  //  module
    //    Int_t              //  not used
    //    Int_t event        //  Event number
    //  Outputs:
    //    none
    //  Return:
    //    test              //  test returns kTRUE if the module contained hits
    //                      //  test returns kFALSE if it did not contain hits

    if(GetDebug(1)) Info("SDigitiseModule","(mod=%p, ,event=%d)",mod,event);
    if(!(mod->GetNhits())){
        if(GetDebug(1)) Info("SDigitiseModule","In event %d module %d there "
                             "are %d hits returning.",event,
                             mod->GetIndex(),mod->GetNhits());
        return;// if module has no hits don't create Sdigits
    } // end if
    SetModuleNumber(mod->GetIndex());
    SetEventNumber(event);
    HitToSDigit(mod);
    WriteSDigits();
    ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::WriteSDigits(){
    //  This function adds each S-Digit to pList
    //  Inputs:
    //    none.
    //  Outputs:
    //    none.
    //  Return:
    //    none
    Int_t ix, nix, iz, niz;
    static AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");

    if(GetDebug(1))Info("WriteSDigits","Writing SDigits for module %d",
                        GetModuleNumber());
    GetMap()->GetMaxMapIndex(niz, nix);
    for(iz=0; iz<niz; iz++)for(ix=0; ix<nix; ix++){
        if(GetMap()->GetSignalOnly(iz,ix)>0.0){
            aliITS->AddSumDigit(*(GetMap()->GetpListItem(iz,ix)));
            if(GetDebug(1)){
                cout <<"AliITSsimulationSPDdubna:WriteSDigits " << iz << "," 
                     << ix << "," << *(GetMap()->GetpListItem(iz,ix)) << endl;
            } // end if GetDebug
        } // end if GetMap()->GetSignalOnly(iz,ix)>0.0
    } // end for iz,ix
    return; 
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::FinishSDigitiseModule(){
    //  This function calls SDigitsToDigits which creates Digits from SDigits
    //  Inputs:
    //    none
    //  Outputs:
    //    none
    //  Return
    //    none

    if(GetDebug(1)) Info("SDigitiseModule","()");
    pListToDigits(); // Charge To Signal both adds noise and
    ClearMap();
    return;
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::DigitiseModule(AliITSmodule *mod,Int_t,
                                              Int_t){
    //  This function creates Digits straight from the hits and then adds
    //  electronic noise to the digits before adding them to pList
    //  Each of the input variables is passed along to HitToSDigit
    //  Inputs:
    //    AliITSmodule *mod     module
    //    Int_t                 Dummy.
    //    Int_t                 Dummy
    //  Outputs:
    //     none.
    //  Return:
    //    none.

    if(GetDebug(1)) Info("DigitiseModule","(mod=%p,,)",mod);
    HitToSDigit(mod);
    pListToDigits();
    ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::HitToSDigit(AliITSmodule *mod){
    // Does the charge distributions using Gaussian diffusion charge charing.
    // Inputs:
    //    AliITSmodule *mod  Pointer to this module
    // Output:
    //    none.
    // Return:
    //    none.
    const Double_t kmictocm = 1.0e-4; // convert microns to cm.
    TObjArray *hits = mod->GetHits();
    Int_t nhits = hits->GetEntriesFast();
    Int_t h,ix,iz,i;
    Int_t idtrack;
    Double_t x0=0.0,x1=0.0,y0=0.0,y1=0.0,z0=0.0,z1=0.0,de=0.0;
    Double_t x,y,z,t,tp,st,dt=0.2,el,sig;
    Double_t thick = kmictocm*GetSeg()->Dy();

    if(GetDebug(1)) Info("HitsToSDigits","(mod=%p) fCoupling=%d",
                         mod,fCoupling);
    if(nhits<=0) return;
    for(h=0;h<nhits;h++){
        if(GetDebug(1)){
            cout << "Hits=" << h << "," << *(mod->GetHit(h)) << endl;
        } // end if GetDebug
        if(!mod->LineSegmentL(h,x0,x1,y0,y1,z0,z1,de,idtrack)) continue;
        st = TMath::Sqrt(x1*x1+y1*y1+z1*z1);
        if(st>0.0){
            st = (Double_t)((Int_t)(1.0E+04*st)); // number of microns
            if(st<=0.0) st = 1.0;
            dt = 1.0/st;
            for(t=0.0;t<1.0;t+=dt){ // Integrate over t
                tp  = t+0.5*dt;
                x   = x0+x1*tp;
                y   = y0+y1*tp;
                z   = z0+z1*tp;
                if(!(GetSeg()->LocalToDet(x,z,ix,iz))) continue; // outside
                el  = GetResp(ix,iz)->GeVToCharge((Double_t)(dt*de));
                if(GetDebug(1)){
                    if(el<=0.0) cout<<"el="<<el<<" dt="<<dt
                                    <<" de="<<de<<endl;
                } // end if GetDebug
                sig = GetResp(ix,iz)->SigmaDiffusion1D(thick + y);
                SpreadCharge(x,z,ix,iz,el,sig,idtrack,h);
            } // end for t
        } else { // st == 0.0 deposit it at this point
            x   = x0;
            y   = y0;
            z   = z0;
            if(!(GetSeg()->LocalToDet(x,z,ix,iz))) continue; // outside
            el  = GetResp(ix,iz)->GeVToCharge((Double_t)de);
            sig = GetResp(ix,iz)->SigmaDiffusion1D(thick + y);
            SpreadCharge(x,z,ix,iz,el,sig,idtrack,h);
        } // end if st>0.0
        // Coupling
        switch (fCoupling) {
        default:
            break;
        case 1: case 3:
            // x is column and z is row (see AliITSsegmentationSPD::GetPadIxz)
            for(i=0;i<GetMap()->GetEntries();i++) 
                if(GetMap()->GetpListItem(i)==0) continue;
                else{
                    GetMap()->GetMapIndex(
                              GetMap()->GetpListItem(i)->GetIndex(),iz,ix);
                    SetCoupling(iz,ix,idtrack,h);
                } // end for i
            break;
        case 2: case 4:
            // x is column and z is row (see AliITSsegmentationSPD::GetPadIxz)
            for(i=0;i<GetMap()->GetEntries();i++) 
                if(GetMap()->GetpListItem(i)==0) continue;
                else{
                    GetMap()->GetMapIndex(
                                GetMap()->GetpListItem(i)->GetIndex(),iz,ix);
                    SetCouplingOld(iz,ix,idtrack,h);
                } // end for i
            break;
        } // end switch
    } // Loop over all hits h
    if(GetDebug(2))Info("HitToSDigit","Finished fCoupling=%d",fCoupling);
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::SpreadCharge(Double_t x0,Double_t z0,
                                            Int_t ix0,Int_t iz0,
					    Double_t el,Double_t sig,Int_t t,
					    Int_t hi){
    // Spreads the charge over neighboring cells. Assume charge is distributed
    // as charge(x,z) = (el/2*pi*sig*sig)*exp(-arg)
    // arg=((x-x0)*(x-x0)/2*sig*sig)+((z-z0*z-z0)/2*sig*sig)
    // Defined this way, the integral over all x and z is el.
    // Inputs:
    //    Double_t x0   x position of point where charge is liberated
    //    Double_t y0   y position of point where charge is liberated
    //    Double_t z0   z position of point where charge is liberated
    //    Int_t    ix0  row of cell corresponding to point x0
    //    Int_t    iz0  columb of cell corresponding to point z0
    //    Double_t el   number of electrons liberated in this step
    //    Double_t sig  Sigma difusion for this step (y0 dependent)
    //    Int_t    t    track number
    //    Int_t    ti   hit track index number
    //    Int_t    hi   hit "hit" index number
    // Outputs:
    //     none.
    // Return:
    //     none.
    const Int_t knx = 3,knz = 2;
    const Double_t kRoot2 = 1.414213562; // Sqrt(2).
    const Double_t kmictocm = 1.0e-4; // convert microns to cm.
    Int_t ix,iz,ixs,ixe,izs,ize;
    Float_t x,z;
    Double_t x1,x2,z1,z2,s,sp;

    if(GetDebug(4)) Info("SpreadCharge","(x0=%e,z0=%e,ix0=%d,iz0=%d,el=%e,"
                         "sig=%e,t=%d,i=%d)",x0,z0,ix0,iz0,el,sig,t,hi);
    if(sig<=0.0) { // if sig<=0 No diffusion to simulate.
        GetMap()->AddSignal(iz0,ix0,t,hi,GetModuleNumber(),el);
        if(GetDebug(2)){
            cout << "sig<=0.0=" << sig << endl;
        } // end if GetDebug
        return;
    } // end if
    sp = 1.0/(sig*kRoot2);
    if(GetDebug(2)){
        cout << "sig=" << sig << " sp=" << sp << endl;
    } // end if GetDebug
    ixs = TMath::Max(-knx+ix0,0);
    ixe = TMath::Min(knx+ix0,GetSeg()->Npx()-1);
    izs = TMath::Max(-knz+iz0,0);
    ize = TMath::Min(knz+iz0,GetSeg()->Npz()-1);
    for(ix=ixs;ix<=ixe;ix++) for(iz=izs;iz<=ize;iz++){
        GetSeg()->DetToLocal(ix,iz,x,z); // pixel center
        x1  = x;
        z1  = z;
        x2  = x1 + 0.5*kmictocm*GetSeg()->Dpx(ix); // Upper
        x1 -= 0.5*kmictocm*GetSeg()->Dpx(ix);  // Lower
        z2  = z1 + 0.5*kmictocm*GetSeg()->Dpz(iz); // Upper
        z1 -= 0.5*kmictocm*GetSeg()->Dpz(iz);  // Lower
        x1 -= x0; // Distance from where track traveled
        x2 -= x0; // Distance from where track traveled
        z1 -= z0; // Distance from where track traveled
        z2 -= z0; // Distance from where track traveled
        s   = 0.25; // Correction based on definision of Erfc
        s  *= TMath::Erfc(sp*x1) - TMath::Erfc(sp*x2);
        if(GetDebug(3)){
            cout <<"el="<<el<<" ix0="<<ix0<<" ix="<<ix<<" x0="<<x<<
                " iz0="<<iz0<<" iz="<<iz<<" z0="<<z<< 
                " sp*x1="<<sp*x1<<" sp*x2="<<sp*x2<<" s="<<s;
        } // end if GetDebug
        s  *= TMath::Erfc(sp*z1) - TMath::Erfc(sp*z2);
        if(GetDebug(3)){
            cout<<" sp*z1="<<sp*z1<<" sp*z2="<<sp*z2<<" s="<<s<< endl;
        } // end if GetDebug
        GetMap()->AddSignal(iz,ix,t,hi,GetModuleNumber(),s*el);
    } // end for ix, iz
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::pListToDigits(){
    // add noise and electronics, perform the zero suppression and add the
    // digit to the list
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    static AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");
    Int_t j,ix,iz;
    Double_t  electronics;
    Double_t sig;
    const Int_t    nmaxtrk=AliITSdigitSPD::GetNTracks();
    static AliITSdigitSPD dig;

    if(GetDebug(1)) Info("pListToDigits","()");
    for(iz=0; iz<GetNPixelsZ(); iz++) for(ix=0; ix<GetNPixelsX(); ix++){
        // Apply Noise/Dead channals and the like
        if(GetResp(ix,iz)->IsPixelDead(GetModuleNumber(),ix,iz)) continue;
        electronics = GetResp(ix,iz)->ApplyBaselineAndNoise();
        UpdateMapNoise(ix,iz,electronics);
        //
        // Apply Threshold and write Digits.
        sig = GetMap()->GetSignalOnly(iz,ix);
        FillHistograms(ix,iz,sig+electronics);
        if(GetDebug(3)){
            cout<<sig<<"+"<<electronics<<">threshold("<<ix<<","<<iz
                <<")="<<GetThreshold(ix,iz) <<endl;
        } // end if GetDebug
        if (sig+electronics <= GetThreshold(ix,iz)) continue;
        dig.SetCoord1(iz);
        dig.SetCoord2(ix);
        dig.SetSignal(1);
        dig.SetSignalSPD((Int_t) GetMap()->GetSignal(iz,ix));
        for(j=0;j<nmaxtrk;j++){
            if (j<GetMap()->GetNEnteries()) {
                dig.SetTrack(j,GetMap()->GetTrack(iz,ix,j));
                dig.SetHit(j,GetMap()->GetHit(iz,ix,j));
            }else { // Default values
                dig.SetTrack(j,-3);
                dig.SetHit(j,-1);
            } // end if GetMap()
        } // end for j
        if(GetDebug(3)){
            cout<<iz<<","<<ix<<","<<*(GetMap()->GetpListItem(iz,ix))<<endl;
        } // end if GetDebug
        aliITS->AddSimDigit(0,&dig);
    } //  for ix/iz
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::CreateHistograms(){
    // create 1D histograms for tests
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //     none.

    if(GetDebug(1)) Info("CreateHistograms","create histograms");

    fHis = new TObjArray(GetNPixelsZ());
    TString fSPDname("spd_");
    for(Int_t i=0;i<GetNPixelsZ();i++) {
        Char_t pixelz[4];
        sprintf(pixelz,"%d",i);
        fSPDname.Append(pixelz);
        fHis->AddAt(new TH1F(fSPDname.Data(),"SPD maps",
                             GetNPixelsX(),0.,(Double_t)GetNPixelsX()),i);
    } // end for i
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::FillHistograms(Int_t ix,Int_t iz,Double_t v){
    // Fill the histogram
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //     none.

    if(!GetHistArray()) return; // Only fill if setup.
    if(GetDebug(2)) Info("FillHistograms","fill histograms");
    GetHistogram(iz)->Fill(ix,v);
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::ResetHistograms(){
    // Reset histograms for this detector
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //     none.

    if(!GetHistArray()) return; // Only fill if setup.
    if(GetDebug(2)) Info("FillHistograms","fill histograms");
    for ( int i=0;i<GetNPixelsZ();i++ ) {
        if (fHis->At(i))    ((TH1F*)fHis->At(i))->Reset();
    } // end for i
}

//______________________________________________________________________
void AliITSsimulationSPDdubna::SetCoupling(Int_t row, Int_t col, Int_t ntrack,
				      Int_t idhit) {
    //  Take into account the coupling between adiacent pixels.
    //  The parameters probcol and probrow are the probability of the
    //  signal in one pixel shared in the two adjacent pixels along
    //  the column and row direction, respectively.
    //  Note pList is goten via GetMap() and module is not need any more.
    //  Otherwise it is identical to that coded by Tiziano Virgili (BSN).
    //Begin_Html
    /*
      <img src="picts/ITS/barimodel_3.gif">
      </pre>
      <br clear=left>
      <font size=+2 color=red>
      <a href="mailto:tiziano.virgili@cern.ch"></a>.
      </font>
      <pre>
    */
    //End_Html
    // Inputs:
    //    Int_t row            z cell index
    //    Int_t col            x cell index
    //    Int_t ntrack         track incex number
    //    Int_t idhit          hit index number
    // Outputs:
    //    none.
    // Return:
    //     none.
    Int_t j1,j2,flag=0;
    Double_t pulse1,pulse2;
    Double_t couplR=0.0,couplC=0.0;
    Double_t xr=0.;

    GetCouplings(couplR,couplC);
    if(GetDebug(3)) Info("SetCoupling","(row=%d,col=%d,ntrack=%d,idhit=%d) "
                         "Calling SetCoupling couplR=%e couplC=%e",
                         row,col,ntrack,idhit,couplR,couplC);
    j1 = row;
    j2 = col;
    pulse1 = GetMap()->GetSignalOnly(row,col);
    pulse2 = pulse1;
    for (Int_t isign=-1;isign<=1;isign+=2){// loop in row direction
        do{
            j1 += isign;
            //   pulse1 *= couplR; 
            xr = gRandom->Rndm();
            //if ((j1<0)||(j1>GetNPixelsZ()-1)||(pulse1<GetThreshold(j1,col))){
            if ((j1<0) || (j1>GetNPixelsZ()-1) || (xr>couplR)){
                j1 = row;
                flag = 1;
            }else{
                UpdateMapSignal(col,j1,ntrack,idhit,pulse1);
                //  flag = 0;
                flag = 1; // only first next!!
            } // end if
        } while(flag == 0);
        // loop in column direction
        do{
            j2 += isign;
            // pulse2 *= couplC; 
            xr = gRandom->Rndm();
            //if((j2<0)||j2>(GetNPixelsX()-1)||pulse2<GetThreshold(row,j2)){
            if ((j2<0) || (j2>GetNPixelsX()-1) || (xr>couplC)){
                j2 = col;
                flag = 1;
            }else{
                UpdateMapSignal(j2,row,ntrack,idhit,pulse2);
                //  flag = 0;
                flag = 1; // only first next!!
            } // end if
        } while(flag == 0);
    } // for isign
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::SetCouplingOld(Int_t row, Int_t col,
                Int_t ntrack,Int_t idhit) {
    //  Take into account the coupling between adiacent pixels.
    //  The parameters probcol and probrow are the fractions of the
    //  signal in one pixel shared in the two adjacent pixels along
    //  the column and row direction, respectively.
    //Begin_Html
    /*
      <img src="picts/ITS/barimodel_3.gif">
      </pre>
      <br clear=left>
      <font size=+2 color=red>
      <a href="mailto:Rocco.Caliandro@ba.infn.it"></a>.
      </font>
      <pre>
    */
    //End_Html
    // Inputs:
    //    Int_t row            z cell index
    //    Int_t col            x cell index
    //    Int_t ntrack         track incex number
    //    Int_t idhit          hit index number
    //    Int_t module         module number
    // Outputs:
    //    none.
    // Return:
    //     none.
    Int_t j1,j2,flag=0;
    Double_t pulse1,pulse2;
    Double_t couplR=0.0,couplC=0.0;

    GetCouplings(couplR,couplC);
    if(GetDebug(3)) Info("SetCouplingOld","(row=%d,col=%d,ntrack=%d,idhit=%d) "
                         "Calling SetCoupling couplR=%e couplC=%e",
                         row,col,ntrack,idhit,couplR,couplC);
    j1 = row;
    j2 = col;
    pulse1 = GetMap()->GetSignalOnly(row,col);
    pulse2 = pulse1;
    for (Int_t isign=-1;isign<=1;isign+=2){// loop in row direction
        do{
            j1 += isign;
            pulse1 *= couplR;
            if ((j1<0)||(j1>GetNPixelsZ()-1)||(pulse1<GetThreshold(j1,col))){
                pulse1 = GetMap()->GetSignalOnly(row,col);
                j1 = row;
                flag = 1;
            }else{
                UpdateMapSignal(col,j1,ntrack,idhit,pulse1);
                flag = 0;
            } // end if
        } while(flag == 0);
        // loop in column direction
        do{
            j2 += isign;
            pulse2 *= couplC;
            if((j2<0)||(j2>(GetNPixelsX()-1))||(pulse2<GetThreshold(row,j2))){
                pulse2 = GetMap()->GetSignalOnly(row,col);
                j2 = col;
                flag = 1;
            }else{
                UpdateMapSignal(j2,row,ntrack,idhit,pulse2);
                flag = 0;
            } // end if
        } while(flag == 0);
    } // for isign
}
