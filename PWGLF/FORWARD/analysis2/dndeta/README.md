# Forward dN/deta for all systems

Here, you will find the Work-in-Progress results for the Forward
dN/deta for all available systems and triggers.

## DISCLAIMER

These results are _Work-in-Progress_ and meant for internal ALICE use
_only_. Points and error bars are subject to change. 

## Website

https://twiki.cern.ch/twiki/bin/view/ALICE/FMDdNdetaAllSystems

## Structure

The gzipped tar-ball contains the results in the form of ROOT scripts in sub-directories like 


> FWD_dNdeta/CORRECTION/COLLISION SYSTEM/COLLISION ENERGY/TRIGGER/RESOLUTION.C


where 

* CORRECTION
  * *normal* MC based secondary correction
  * *nosec* Empirical secondary correction 
* COLLISION SYSTEM
  * *pp*
  * *pPb* 
  * *Pbp* 
  * *PbPb* 
* COLLISION ENERGY
  * *0900* pp only 
  * *2760* pp and PbPb
  * *5023* pPb Pbp
  * *7000* pp only
  * *8000* pp only
* TRIGGER
  * *INEL*  pp only
  * *INELGt0* pp only
  * *V0AND* pp only (NSD)
  * *CENTMB* pPb and Pbp
  * *CENTV0M* pPb and Pbp 
  * *CENTV0A* pPb only 
  * *CENTZNA* pPb only
  * *CENTV0C* Pbp only 
  * *CENTZNC* Pbp only
  * *CENT* PbPb only
* RESOLUTION
  * *full* Full &eta; resolution
  * *rebin* Rebinned &eta; resolution

Each script accepts pointers to a THStack and a TLegend (optionally
null) and a marker style.

## Scripts

* `Draw.C(SYSTEM,SNN,TRIGGER,EMPIRICAL,REBINNED)` will draw a given
  result.
* `Scalebypp.C(SYSTEM,SNN,PPSNN,TRIGGER,PPTRIGGER,ETASHIFT,EMPIRICAL,REBINNED)`
  will scale a given result (SYSTEM,SNN,TRIGGER) by the given pp
  result (PPSNN,PPTRIGGER), optionally shifting in eta by ETASHIFT
* Both of these accepts and optional last arguement SYMMETRICE which
  for pPb will draw the combined pPb/Pbp result
* `DrawATLAS.C` draws the ALICE and ATLAS results for pPb with the V0A
  centrality selector.
* `DrawATLAS.C` also defines the function `DrawAvgATLAS()` which
  averages the ATLAS centrality bins to agree with the ALICE
  centrality bins.
* `OtherData.C` contains the class `RefData.C` to access data from
  sources than the ALICE forward detectors.
* `Drawer.C` contains the class `Drawer` to draw the results.

## Usage

To use from an interactive session do 


    Root> THStack* stack = new THStack("stack", "");
    Root> TLegend* leg   = new TLegend(0.35,0.2,0.55,0.8);
    Root> .L FWD_dNdeta/nosec/pPb/5023/CENTMB/full.C
    Root> full(stack, leg, 0);
    Root> stack->Draw("nostack");
    Root> if (leg) leg->Draw();
    

or in code 


    #ifndef __CINT__
    #include <TLegend.h>
    #include <THStack.h>
    #include <TROOT.h>
    #include <TString.h>
    #endif
    
    void DrawFWDdNdetapPbCENTMB() {
      Int_t    marker = 0;
      TCanvas* canvas = new TCanvas("c","C");
      THStack* stack  = new THStack("stack", "");
      TLegend* leg    = new TLegend(0.35,0.2,0.55,0.8);
      TString  args   = Form("(THStack*)%p,(TLegend*)%p,%d",stack,leg,marker);
      gROOT->Macro(Form("FWD_dNdeta/nosec/pPb/5023/CENTMB/full.C(%s)", args.Data())));
      stack->Draw("nostack");
      if (leg) leg->Draw();
    }


## Contact

Christian Holm Christensen <cholm@nbi.dk>

