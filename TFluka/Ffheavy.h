extern "C" {
//*$ create fheavy.add
//*copy fheavy
//*
//*=== fheavy ===========================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     include file: fheavy                                             *
//*                                                                      *
//*     created  on  5 april 1990     by   alfredo ferrari, infn milan   *
//*                                                                      *
//*     last change on   26-jul-97    by   alfredo ferrari, infn milan   *
//*                                                                      *
//*     included in the following subroutines or functions: not updated  *
//*                                                                      *
//*     description of the common block(s) and variable(s)               *
//*                                                                      *
//*     /fheavy/ is the storage for heavy secondaries created in the     *
//*              nuclear evaporation                                     *
//*        npheav     = number of secondaries                            *
//*        kheavy(ip) = type of the secondary ip                         *
//*                   ( 3 = deuteron, 4 = 3-h, 5 = 3-he, 6 = 4-he,       *
//*                     7-12 = "heavy" fragment specified by ibheav and  *
//*                     icheav )                                         *
//*        cxheav(ip) = direction cosine of the secondary ip             *
//*                     with respect to x-axis                           *
//*        cyheav(ip) = direction cosine of the secondary ip             *
//*                     with respect to y-axis                           *
//*        czheav(ip) = direction cosine of the secondary ip             *
//*                     with respect to z-axis                           *
//*        tkheav(ip) = kinetic energy of secondary ip                   *
//*        pheavy(ip) = momentum of the secondary ip                     *
//*        wheavy(ip) = weight of the secondary ip                       *
//*        agheav(ip) = "age" of the secondary ip with respect to the    *
//*                     interaction time                                 *
//*        amheav(kp) = atomic masses of the twelve types of evaporated  *
//*                     or fragmented or fissioned particles             *
//*        amnhea(kp) = nuclear masses of the twelve types of evaporated *
//*                     or fragmented or fissioned particles             *
//*     bhheav(jp,kp) = (nuclear) binding energy of the jp_th hyperon of *
//*                     the kp-type heavy particle                       *
//*        anheav(kp) = name of the kp-type heavy particle               *
//*        icheav(kp) = charge of the kp-type heavy particle             *
//*        ibheav(kp) = mass number of the kp-type heavy particle        *
//*        imheav(kp) = isomeric state of the kp-type heavy particle     *
//*        ihheav(kp) = number of hyperons of the kp-type heavy particle *
//*     khheav(jp,kp) = id of the jp_th hyperon of the kp-type heavy     *
//*                     particle                                         *
//*        infhea(ip) = possible extra infos for the ip_th secondary     * 2006.3
//*   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   *
//*   !!! there is now the possibility to produce up to 6 "heavy" !!!!   *
//*   !!! fragments besides the residual nucleus recorded in      !!!!   *
//*   !!! resnuc: they are identified by indeces 7-12, of course  !!!!   *
//*   !!! the corresponding physical properties (z,a,m..) must be !!!!   *
//*   !!! updated every time they are produced                    !!!!   *
//*   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   *
//*----------------------------------------------------------------------*
//*
const Int_t mxheav = 100;
const Int_t kxheav = 30;

typedef struct {
    Double_t cxheav[mxheav];
    Double_t cyheav[mxheav];
    Double_t czheav[mxheav];
    Double_t tkheav[mxheav];
    Double_t pheavy[mxheav];
    Double_t wheavy[mxheav];
    Double_t agheav[mxheav];
    Double_t bhheav[kxheav][ihypmx];
    Double_t amheav[kxheav];
    Double_t amnhea[kxheav];
    Int_t    kheavy[mxheav];
    Int_t    infhea[mxheav]; // 2006.3
    Int_t    icheav[kxheav];
    Int_t    ibheav[kxheav];
    Int_t    imheav[kxheav];
    Int_t    ihheav[kxheav];
    Int_t    khheav[kxheav][ihypmx];
    Int_t    npheav;
} fheavyCommon;
#define FHEAVY COMMON_BLOCK(FHEAVY,fheavy)
COMMON_BLOCK_DEF(fheavyCommon,FHEAVY);

typedef struct {
   Char_t   anheav[kxheav][8];
} fheavcCommon;
#define FHEAVC COMMON_BLOCK(FHEAVC,fheavc)
COMMON_BLOCK_DEF(fheavcCommon,FHEAVC);
}
