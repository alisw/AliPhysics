extern "C" {
//*$ create iounit.add
//*copy iounit
//*                                                                      *
//*=== iounit ===========================================================*
//*                                                                      *
//*----------------------------------------------------------------------*
//*                                                                      *
//*      iounit: included in any routine                                 *
//*                                                                      *
//*     created on   01 june 1990    by    alfredo ferrari, infn -milan  *
//*                                                                      *
//*     last change on 20-jan-99     by    alfredo ferrari               *
//*                                                                      *
//*         lunin  = standard input unit                                 *
//*         lunout = standard output unit                                *
//*         lunerr = standard error unit                                 *
//*         lunber = input file for bertini nuclear data                 *
//*         lunech = echo file for pegs dat                              *
//*         lunflu = input file for photoelectric edges and x-ray fluo-  *
//*                  rescence data                                       *
//*         lungeo = scratch file for combinatorial geometry             *
//*         lunpmf = input file for pegs material data                   *
//*         lunran = output file for the final random number seed        *
//*         lunxsc = input file for low energy neutron cross sections    *
//*         lundet = output file for the detect option                   *
//*         lunray = output file for ray-tracing options                 *
//*         lunrdb = unit number for reading (extra) auxiliary external  *
//*                  files to be closed just after reading               *
//*         lunrd2 = unit number for reading (extra) auxiliary external  *
//*                  files to be closed just after reading               *
//*         lunscr = unit number to be used for temporary scratch files  *
//*         lunpgo = output file for plotgeom                            *
//*         lunpgs = store (formatted/unformatted) file for plotgeom     *
//*                                                                      *
//*----------------------------------------------------------------------*
//*                                                                      *
const Int_t lunin = 5;
//* start_vax_seq
//*     parameter ( lunout =  6 )
//* end_vax_seq
//* start_ibm_seq
//*     parameter ( lunout =  6 )
//* end_ibm_seq
//* start_unix_seq
const Int_t lunout = 11;
//* end_unix_seq
const Int_t lunerr = 15;
const Int_t lunber = 14;
const Int_t lunech = 8;
const Int_t lunflu = 13;
const Int_t lungeo = 16;
const Int_t lunpmf = 12;
const Int_t lunran = 2;
const Int_t lunxsc = 9;
const Int_t lundet = 17;
const Int_t lunray = 10;
const Int_t lunrdb = 1;
const Int_t lunrd2 = 18;
const Int_t lundpm = 19;
const Int_t lunpgo = 7;
const Int_t lunpgs = 4;
const Int_t lunscr = 3;
}
