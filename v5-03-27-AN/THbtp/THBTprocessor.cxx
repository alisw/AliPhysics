#include "THBTprocessor.h"
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_____________________________________________________________________________
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class THBTprocessor                                                       //
//                                                                           //
// Wrapper class to HBT processor fortran code.                              //
// For more information see AliGenHBTprocessor class                         //
// HBT processor is written by Lanny Ray                                     //
//                                                                           //
// Comunication is done via COMMON BLOCKS declared in HBTprocCOMMON.h        //
// using cfortran.h routines                                                 //
// User should use class AliGenHBTprocessor and all their interface          //
// see there for more description                                            //
//                                                                           //
// Wrapper class written by                                                  //
// Piotr Krzysztof Skowronski (Piotr.Skowronski@cern.ch)                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <Riostream.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TParticle.h>

#include "HBTprocCOMMON.h"

#ifndef WIN32
# define hbtprocessor hbtprocessor_
# define ctest ctest_
# define type_of_call
#else
# define hbtprocessor HBTPROCESSOR
# define ctest CTEST
# define type_of_call _stdcall
#endif    


ClassImp(THBTprocessor)
 
extern "C" void  type_of_call hbtprocessor();   
extern "C" void  type_of_call ctest();


/*************************************************/
THBTprocessor::THBTprocessor()// it is better not to intialize it:TGenerator("THBTprocessor","THBTprocessor")
                              //because it allocates memmory for TObjArray::fParticles which is not used in our case
                              //and we are using TClonesArray in import paerticles
 {
 //constructor
  PARAMETERS.ALICE            = 1; //flag that we are working in AliRoot (0==STAR stand alone mode)
  PARAMETERS.pi               = TMath::Pi();//3.141592654;
  PARAMETERS.rad              = 180.0/TMath::Pi();
  PARAMETERS.hbc              = 0.19732891;
  
  fParticles = 0;
  Initialize(); //Enforce initialization (cleaning all commons)
 }
/*****************************************************************************************/

#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
void THBTprocessor::GenerateEvent() const
#else
void THBTprocessor::GenerateEvent()
#endif
{
//Starts processing

  Info("GenerateEvent","Entering Fortran");
  ctest();
  hbtprocessor();
  Info("GenerateEvent","Exited Fortran");
  if(PARAMETERS.errorcode)
    {
      TString message("HBT Processor (fortran part) finished with errors\n");
      message+="Error code is ";
      message+=PARAMETERS.errorcode;
      message+="\n";
      message+="See hbt_simulation.out file for more detailed information.";
      Fatal("GenerateEvent","%s",message.Data());
    }
  else
    Info("GenerateEvent","GOOD ! HBT Processor finished without errors");
}
/*****************************************************************************************/

void THBTprocessor::Initialize() const
{ 
  //IT RESETS ALL THE PREVIOUS SETTINGS
  //Call this method to set default values in PARAMETERS & MESH
  //and zero other common block
  
  if(gDebug) 
   Info("Initialize","Setting Default valuses in all COMMON BLOCKS");
  
  PARAMETERS.ref_control      = 2;
  PARAMETERS.switch_1d        = 3;
  PARAMETERS.switch_3d        = 1;
  PARAMETERS.switch_type      = 3;
  PARAMETERS.switch_coherence = 0;
  PARAMETERS.switch_coulomb   = 3;
  PARAMETERS.switch_fermi_bose= 1;
  PARAMETERS.trk_accep        = 1.0;
  PARAMETERS.print_full       = 1;
  PARAMETERS.print_sector_data= 1;

  PARAMETERS.n_pid_types      = 2;
  PARAMETERS.pid[0]           = 8;
  PARAMETERS.pid[1]           = 9;
  PARAMETERS.deltap           = 0.1;
  PARAMETERS.maxit            = 50;
  PARAMETERS.delchi           = 1.0;
  PARAMETERS.irand            = 76564;
  PARAMETERS.lambda           = 0.6;
  PARAMETERS.R_1d             = 7.0;
  PARAMETERS.Rside            = 6.0;

  PARAMETERS.Rout             = 7.0;
  PARAMETERS.Rlong            = 4.0;
  PARAMETERS.Rperp            = 6.0;
  PARAMETERS.Rparallel        = 4.0;
  PARAMETERS.R0               = 4.0;
  PARAMETERS.Q0               = 9.0;

  PARAMETERS.n_part_1_trk         = 0;
  PARAMETERS.n_part_2_trk         = 0;
  PARAMETERS.n_part_tot_trk       = 0;
  PARAMETERS.n_part_used_1_trk    = 0;

  PARAMETERS.n_part_used_2_trk    = 0;
  PARAMETERS.n_part_1_trk2        = 0;
  PARAMETERS.n_part_2_trk2        = 0;
  PARAMETERS.n_part_tot_trk2      = 0;
  PARAMETERS.n_part_used_1_trk2   = 0;
  PARAMETERS.n_part_used_2_trk2   = 0;
  PARAMETERS.n_part_used_1_ref    = 0;
  PARAMETERS.n_part_used_2_ref    = 0;
  PARAMETERS.n_part_used_1_inc    = 0;
  PARAMETERS.n_part_used_2_inc    = 0;
 
  PARAMETERS.num_pairs_like       = 0;
  PARAMETERS.num_pairs_unlike     = 0;
  PARAMETERS.num_pairs_like_ref   = 0;
  PARAMETERS.num_pairs_like_inc   = 0;
  PARAMETERS.num_pairs_unlike_inc = 0;
  PARAMETERS.event_line_counter   = 0;
  PARAMETERS.file10_line_counter  = 0;

  PARAMETERS.chisq_wt_like_1d         = 1.0;
  PARAMETERS.chisq_wt_unlike_1d       = 1.0;
  PARAMETERS.chisq_wt_like_3d_fine    = 1.0;

  PARAMETERS.chisq_wt_unlike_3d_fine  = 1.0;
  PARAMETERS.chisq_wt_like_3d_coarse  = 1.0;
  PARAMETERS.chisq_wt_unlike_3d_coarse= 1.0;
  PARAMETERS.chisq_wt_hist1_1         = 1.0;
  PARAMETERS.chisq_wt_hist1_2         = 1.0; // /////internal comment 25 fields

/*********************************************/

 
  MESH.n_pt_bins        = 50;                            //OK
  MESH.pt_min           = 0.1;  //Pt in GeV/c            //OK
  MESH.pt_max           = 0.98; //Pt in GeV/c
  MESH.n_phi_bins       = 50;                          //OK
  MESH.phi_min          = 0.0;                              //OK
  MESH.phi_max          = 360.0;                          //OK
  MESH.n_eta_bins       = 50;                          //OK
  MESH.eta_min          =-1.5;                          //OK
  MESH.eta_max          = 1.5;                          //OK
  MESH.n_1d_fine        = 10;                          //OK
  MESH.binsize_1d_fine  = 0.01;                         //ok 
  MESH.n_1d_coarse      = 2;                          //O
  MESH.binsize_1d_coarse= 0.05;                       //ok
  MESH.n_3d_fine        = 8;                          //OK
  MESH.binsize_3d_fine  = 0.01;                       //ok
  MESH.n_3d_coarse      = 2;                          //OK
  MESH.binsize_3d_coarse= 0.08;                       //ok
  MESH.n_3d_fine_project= 3;                          //OK
  MESH.n_px_bins        = 20;                          //OK
  MESH.px_min           =-1.0;                          //OK
  MESH.px_max           = 1.0;                          //OK
  MESH.n_py_bins        = 20;                          //OK
  MESH.py_min           =-1.0;                          //OK
  MESH.py_max           = 1.0;                          //OK
  MESH.n_pz_bins        = 70;                          //OK
  MESH.pz_min           =-3.6;                          //OK
  MESH.pz_max           = 3.6;                          //OK
  
/*********************************************/

  Int_t i; //loop variable
  
  for (i =0; i<TRK_MAXLEN;i++)
   {
    TRACK1.trk_id[i] = 0;
    TRACK1.trk_px_sec[i] = 0; 
    TRACK1.trk_py_sec[i] = 0;
    TRACK1.trk_pz_sec[i] = 0;
    TRACK1.trk_sector[i] = 0;
    TRACK1.trk_flag[i] = 0;
    TRACK1.trk_out_flag[i] = 0;
    TRACK1.trk_merge_flag[i] = 0;
    TRACK1.trk_ge_pid[i] = 0;
    TRACK1.trk_start_vertex[i] = 0;
    TRACK1.trk_stop_vertex[i] = 0;
    TRACK1.trk_event_line[i] = 0;
    TRACK1.trk_px[i] = 0;
    TRACK1.trk_py[i] = 0;
    TRACK1.trk_pz[i] = 0;
    TRACK1.trk_E[i] = 0;
    TRACK1.trk_pt[i] = 0;
    TRACK1.trk_phi[i] = 0;
    TRACK1.trk_eta[i] = 0; 
   }
  
/*********************************************/

  for (i =0; i<TRK2_MAXLEN;i++)
   {
    TRACK2.trk2_id[i] = 0;
    TRACK2.trk2_px_sec[i] = 0; 
    TRACK2.trk2_py_sec[i] = 0;
    TRACK2.trk2_pz_sec[i] = 0;
    TRACK2.trk2_sector[i] = 0;
    TRACK2.trk2_flag[i] = 0;
    TRACK2.trk2_out_flag[i] = 0;
    TRACK2.trk2_merge_flag[i] = 0;
    TRACK2.trk2_ge_pid[i] = 0;
    TRACK2.trk2_start_vertex[i] = 0;
    TRACK2.trk2_stop_vertex[i] = 0;
    TRACK2.trk2_event_line[i] = 0;
    TRACK2.trk2_px[i] = 0;
    TRACK2.trk2_py[i] = 0;
    TRACK2.trk2_pz[i] = 0;
    TRACK2.trk2_E[i] = 0;
    TRACK2.trk2_pt[i] = 0;
    TRACK2.trk2_phi[i] = 0;
    TRACK2.trk2_eta[i] = 0; 
   }

/*********************************************/

  for (i =0; i<SEC_MAXLEN;i++)
   {
    SEC_TRK_MAP.stm_sec_id [i] = 0;
    SEC_TRK_MAP.stm_n_trk_sec[i] = 0; 
    SEC_TRK_MAP.stm_flag[i] = 0;
    
    for (Int_t j=0; j<MAX_TRK_SEC;j++)
       SEC_TRK_MAP.stm_track_id[i][j] = 0;
   }

/*********************************************/

  for (i =0; i<SEC_MAXLEN2;i++)
   {
    SEC_TRK_MAP2.stm_sec_id2[i] = 0;
    SEC_TRK_MAP2.stm_n_trk_sec2[i] = 0; 
    SEC_TRK_MAP2.stm_flag2[i] = 0;
    
    for (Int_t j=0; j<MAX_TRK_SEC;j++)
       SEC_TRK_MAP2.stm_track_id2[i][j] = 0;
   }

/*********************************************/

  for (i =0; i<PART_MAXLEN;i++)
   {
     PARTICLE.part_id[i] = 0;
     PARTICLE.part_charge[i] = 0;
     PARTICLE.part_mass[i] = 0;
     PARTICLE.part_lifetime[i] = 0;
   }


/*********************************************/
  for (i =0; i<MAX_C2_1D;i++)
   {
     CORRELATIONS.c2mod_like_1d[i] = 0; 
     CORRELATIONS.c2mod_unlike_1d[i] = 0;
     CORRELATIONS.c2fit_like_1d[i] = 0;
     CORRELATIONS.c2fit_unlike_1d[i] = 0;
     CORRELATIONS.c2err_like_1d[i] = 0;
     CORRELATIONS.c2err_unlike_1d[i] = 0;
   }
/*********************************************/
  for (i =0; i<MAX_C2_3D;i++)
   for (Int_t j =0; j<MAX_C2_3D;j++)
    for (Int_t k =0; k<MAX_C2_3D;k++)
     {
      CORRELATIONS.c2mod_like_3d_fine[i][j][k] = 0.0;
      CORRELATIONS.c2mod_unlike_3d_fine[i][j][k] = 0.0;
      CORRELATIONS.c2mod_like_3d_coarse[i][j][k] = 0.0;
      CORRELATIONS.c2mod_unlike_3d_coarse[i][j][k] = 0.0;
      CORRELATIONS.c2fit_like_3d_fine[i][j][k] = 0.0;
      CORRELATIONS.c2fit_unlike_3d_fine[i][j][k] = 0.0;
      CORRELATIONS.c2fit_like_3d_coarse[i][j][k] = 0.0;
      CORRELATIONS.c2fit_unlike_3d_coarse[i][j][k] = 0.0;
      CORRELATIONS.c2err_like_3d_fine[i][j][k] = 0.0;
      CORRELATIONS.c2err_unlike_3d_fine[i][j][k] = 0.0;
      CORRELATIONS.c2err_like_3d_coarse[i][j][k] = 0.0;
      CORRELATIONS.c2err_unlike_3d_coarse[i][j][k] = 0.0;
     }
/*********************************************/

   EVENT_SUMMARY.niter_mean = 0.0; 
   EVENT_SUMMARY.niter_rms = 0.0;
 
   EVENT_SUMMARY.npart1_mean = 0.0; 
   EVENT_SUMMARY.npart1_rms = 0.0;
 
   EVENT_SUMMARY.npart2_mean = 0.0; 
   EVENT_SUMMARY.npart2_rms = 0.0;
 
   EVENT_SUMMARY.npart_tot_mean = 0.0; 
   EVENT_SUMMARY.npart_tot_rms = 0.0;
 
   EVENT_SUMMARY.nsec_flag_mean = 0.0; 
   EVENT_SUMMARY.nsec_flag_rms = 0.0;
 
   EVENT_SUMMARY.frac_trks_out_mean = 0.0;
   EVENT_SUMMARY.frac_trks_out_rms = 0.0;
 
   EVENT_SUMMARY.frac_trks_flag_mean = 0.0;
   EVENT_SUMMARY.frac_trks_flag_rms = 0.0;
 
   EVENT_SUMMARY.chi_l1d_mean = 0.0;
   EVENT_SUMMARY.chi_l1d_rms = 0.0;

   EVENT_SUMMARY.chi_u1d_mean = 0.0;
   EVENT_SUMMARY.chi_u1d_rms = 0.0;

 for (i =0; i<MAX_EVENTS;i++) 
  {
    EVENT_SUMMARY.n_part_used_1_store[i] = 0.0; 
    EVENT_SUMMARY.n_part_used_2_store[i] = 0.0; 
    EVENT_SUMMARY.n_part_tot_store[i] = 0.0; 
    EVENT_SUMMARY.num_sec_flagged_store[i] = 0.0; 
    EVENT_SUMMARY.frac_trks_out[i] = 0.0; 
    EVENT_SUMMARY.frac_trks_flag[i] = 0.0;
    EVENT_SUMMARY.chisq_like_1d_store[i] = 0.0;         
    EVENT_SUMMARY.num_iter[i] = 0.0; 
    EVENT_SUMMARY.chisq_unlike_1d_store[i] = 0.0;       
  } 
/*********************************************/
 for (i =0; i<MAX_C2_COUL;i++) 
  {
      COULMB.c2_coul_like[i] = 0.0;
      COULMB.c2_coul_unlike[i] = 0.0;
      COULMB.q_coul[i] = 0.0; 
 
  }
/*********************************************/
 for (i =0; i<MAX_H_1D;i++) 
   {  
      HISTOGRAMS.hist1_pt_1[i] = 0;
      HISTOGRAMS.hist1_phi_1[i] = 0;
      HISTOGRAMS.hist1_eta_1[i] = 0;
      HISTOGRAMS.hist1_pt_2[i] = 0;
      HISTOGRAMS.hist1_phi_2[i] = 0;
      HISTOGRAMS.hist1_eta_2[i] = 0;
      HISTOGRAMS.htmp1_pt_1[i] = 0;
      HISTOGRAMS.htmp1_phi_1[i] = 0;
      HISTOGRAMS.htmp1_eta_1[i] = 0;
      HISTOGRAMS.htmp1_pt_2[i] = 0;
      HISTOGRAMS.htmp1_phi_2[i] = 0;
      HISTOGRAMS.htmp1_eta_2[i] = 0;
      HISTOGRAMS.href1_pt_1[i] = 0;
      HISTOGRAMS.href1_phi_1[i] = 0;
      HISTOGRAMS.href1_eta_1[i] = 0;
      HISTOGRAMS.href1_pt_2[i] = 0;
      HISTOGRAMS.href1_phi_2[i] = 0;
      HISTOGRAMS.href1_eta_2[i] = 0;
      HISTOGRAMS.hinc1_pt_1[i] = 0;
      HISTOGRAMS.hinc1_phi_1[i] = 0;
      HISTOGRAMS.hinc1_eta_1[i] = 0;
      HISTOGRAMS.hinc1_pt_2[i] = 0;
      HISTOGRAMS.hinc1_phi_2[i] = 0;
      HISTOGRAMS.hinc1_eta_2[i] = 0;
      HISTOGRAMS.hist_like_1d[i] = 0;
      HISTOGRAMS.hist_unlike_1d[i] = 0;
      HISTOGRAMS.htmp_like_1d[i] = 0;
      HISTOGRAMS.htmp_unlike_1d[i] = 0;
      HISTOGRAMS.href_like_1d[i] = 0;
      HISTOGRAMS.href_unlike_1d[i] = 0;
      HISTOGRAMS.hinc_like_1d[i] = 0;
      HISTOGRAMS.hinc_unlike_1d[i] = 0;
  }

 for (i =0; i<MAX_H_3D;i++) 
   for (Int_t j =0; j<MAX_H_3D;j++) 
     for (Int_t k =0; k<MAX_H_3D;k++) 
       {
          HISTOGRAMS.hist_like_3d_fine[i][j][k] = 0;
          HISTOGRAMS.hist_unlike_3d_fine[i][j][k] = 0;
          HISTOGRAMS.hist_like_3d_coarse[i][j][k] = 0;
          HISTOGRAMS.hist_unlike_3d_coarse[i][j][k] = 0;
          HISTOGRAMS.htmp_like_3d_fine[i][j][k] = 0;
          HISTOGRAMS.htmp_unlike_3d_fine[i][j][k] = 0;
          HISTOGRAMS.htmp_like_3d_coarse[i][j][k] = 0;
          HISTOGRAMS.htmp_unlike_3d_coarse[i][j][k] = 0;
          HISTOGRAMS.href_like_3d_fine[i][j][k] = 0;
          HISTOGRAMS.href_unlike_3d_fine[i][j][k] = 0;
          HISTOGRAMS.href_like_3d_coarse[i][j][k] = 0;
          HISTOGRAMS.href_unlike_3d_coarse[i][j][k] = 0;
          HISTOGRAMS.hinc_like_3d_fine[i][j][k] = 0;
          HISTOGRAMS.hinc_unlike_3d_fine[i][j][k] = 0;
          HISTOGRAMS.hinc_like_3d_coarse[i][j][k] = 0;
          HISTOGRAMS.hinc_unlike_3d_coarse[i][j][k] = 0;
       }
/*********************************************/


/*********************************************/

//  cout<<" FINISHED"<<endl;
  
}


/*****************************************************************************************/


Int_t THBTprocessor::ImportParticles(TClonesArray *particles, Option_t */*option*/)
 {
  //Copy particle data into TClonesArray
  if (particles == 0) return 0;
  TClonesArray &rparticles = *particles;
  rparticles.Clear();
 
  Int_t nrpart = 0;
  for (Int_t i=0; i < TRK_MAXLEN; i++) 
   {
   
    
      if (TRACK1.trk_E[i] == 0.) continue;
    
      Float_t px   = TRACK1.trk_px[i];
      Float_t py   = TRACK1.trk_py[i];
      Float_t pz   = TRACK1.trk_pz[i];
//    Float_t pE   = TRACK.trk_E[i];
      Float_t mass = PARTICLE.part_mass[TRACK1.trk_ge_pid[i]];
    
      new(rparticles[nrpart++]) TParticle(0,0,0,0,0,0,px,py,pz,
                                         TMath::Sqrt(mass*mass+px*px+py*py+pz*pz),
                                         0,0,0,0);
   }
  return nrpart;        
 }

/*****************************************************************************************/
TObjArray * THBTprocessor::ImportParticles(Option_t */*option*/)
 {
  //Copy particle data into TObjArray
  fParticles->Clear();
 
  for (Int_t i=0; i < TRK_MAXLEN; i++) 
   {
   
    
      if (TRACK1.trk_E[i] == 0.) continue;
    
      Float_t px   = TRACK1.trk_px[i];
      Float_t py   = TRACK1.trk_py[i];
      Float_t pz   = TRACK1.trk_pz[i];
//    Float_t pE   = TRACK.trk_E[i];
      Float_t mass = PARTICLE.part_mass[TRACK1.trk_ge_pid[i]];
    
      TParticle * p =new TParticle(0,0,0,0,0,0,px,py,pz,
				   TMath::Sqrt(mass*mass+px*px+py*py+pz*pz),
				   0,0,0,0);
      fParticles->Add(p);
   }
  return fParticles;        
 }

/*****************************************************************************************/

void THBTprocessor::PrintEvent() const
 {
 //Prints all particles (common block data)  
   cout<<"Print Event"<<endl;
   for (Int_t i=0; i<TRK_MAXLEN;i++)
    {
      if(TRACK1.trk_E[i]==0.) continue;

      cout<<"trk_id: "<<TRACK1.trk_id[i]<<"  trk_px :"<<TRACK1.trk_px[i]<<"  trk_py :"<<TRACK1.trk_py[i]<<"  trk_pz :"<<TRACK1.trk_pz[i]<<endl;
      cout<<"                trk_E: "<<TRACK1.trk_E[i]<<"  trk_pt: "<<TRACK1.trk_pt[i]<<"  trk_phi: "<<TRACK1.trk_phi[i]<<"  trk_eta: "<<TRACK1.trk_eta[i]<<endl;
    }
 }


/*****************************************************************************************/
void THBTprocessor::DumpSettings() const
{
 //prints values set in common blocks
  ctest();
}
void THBTprocessor::SetTrackRejectionFactor(Float_t trf) const 
{
 //setter
 PARAMETERS.trk_accep = trf;
}
void THBTprocessor::SetRefControl(Int_t rc) const 
{
 //setter
 PARAMETERS.ref_control = rc;
}
void THBTprocessor::SetPIDs(Int_t pid1,Int_t pid2) const 
{
 //setter
 PARAMETERS.pid[0]=pid1; PARAMETERS.pid[1]=pid2;
}
void THBTprocessor::SetNPIDtypes(Int_t npidt)const 
{
 //setter
 PARAMETERS.n_pid_types = npidt;
}
void THBTprocessor::SetDeltap(Float_t deltp) const 
{
 //setter
 PARAMETERS.deltap=deltp; 
}
void THBTprocessor::SetMaxIterations(Int_t maxiter) const 
{
 //setter
 PARAMETERS.maxit = maxiter;
}
void THBTprocessor::SetDelChi(Float_t dc)const 
{
 //setter
 PARAMETERS.delchi = dc;
}
void THBTprocessor::SetIRand(Int_t irnd) const 
{
 //setter
 PARAMETERS.irand = irnd;
}
void THBTprocessor::SetLambda(Float_t lam) const 
{ 
 //setter
 PARAMETERS.lambda = lam;
}
void THBTprocessor::SetR1d(Float_t r) const 
{
 //setter
 PARAMETERS.R_1d=r;
}
void THBTprocessor::SetRSide(Float_t rs) const 
{
 //setter
 PARAMETERS.Rside=rs;
}
void THBTprocessor::SetROut(Float_t ro) const 
{
 //setter
 PARAMETERS.Rout=ro;
}
void THBTprocessor::SetRLong(Float_t rl) const 
{
 //setter
 PARAMETERS.Rlong=rl;
}
void THBTprocessor::SetRPerp(Float_t rp) const 
{
 //setter
 PARAMETERS.Rperp=rp;
}
void THBTprocessor::SetRParallel(Float_t rprl) const 
{
 //setter
 PARAMETERS.Rparallel=rprl;
}
void THBTprocessor::SetR0(Float_t r0) const 
{
 //setter
 PARAMETERS.R0=r0;
}
void THBTprocessor::SetQ0(Float_t q0) const 
{
 //setter
PARAMETERS.Q0=q0;
}
void THBTprocessor::SetSwitch1D(Int_t s1d) const 
{
 //setter
 PARAMETERS.switch_1d = s1d;
}
void THBTprocessor::SetSwitch3D(Int_t s3d) const 
{
 //setter
 PARAMETERS.switch_3d = s3d;
}
void THBTprocessor::SetSwitchType(Int_t st) const 
{
 //setter
 PARAMETERS.switch_type = st;
}

void THBTprocessor::SetSwitchCoherence(Int_t sc) const 
{
 //setter
 PARAMETERS. switch_coherence = sc;
}
void THBTprocessor::SetSwitchCoulomb(Int_t scol) const 
{
 //setter
 PARAMETERS. switch_coulomb = scol;
}
void THBTprocessor::SetSwitchFermiBose(Int_t sfb) const 
{
 //setter
 PARAMETERS.switch_fermi_bose = sfb;
}
      
void THBTprocessor::SetPtRange(Float_t ptmin, Float_t ptmax)const  //Pt in GeV/c
{ 
 //setter
  MESH.pt_min=ptmin;MESH.pt_max=ptmax;
}

void THBTprocessor::SetPxRange(Float_t pxmin, Float_t pxmax) const
{ 
 //setter
 MESH.px_min=pxmin;MESH.px_max=pxmax;
}

void THBTprocessor::SetPyRange(Float_t pymin, Float_t pymax)const 
{ 
 //setter
 MESH.py_min=pymin;MESH.py_max=pymax;
}

void THBTprocessor::SetPzRange(Float_t pzmin, Float_t pzmax) const
{ 
 //setter
 MESH.pz_min=pzmin;MESH.pz_max=pzmax;
}

void THBTprocessor::SetPhiRange(Float_t phimin, Float_t phimax)const //Angle in degrees
{ 
 //setter
 MESH.phi_min=phimin;MESH.phi_max=phimax;
}

void THBTprocessor::SetEtaRange(Float_t etamin, Float_t etamax)const //Pseudorapidity  !!!!!!!!!
{ 
 //setter
 MESH.eta_min=etamin;MESH.eta_max=etamax;
}

void THBTprocessor::SetNPtBins(Int_t nptbin)const
{
//setter
MESH.n_pt_bins=nptbin;
}
void THBTprocessor::SetNPhiBins(Int_t nphibin)const
{
//setter
MESH.n_phi_bins=nphibin;
}
void THBTprocessor::SetNEtaBins(Int_t netabin)const
{
//setter
MESH.n_eta_bins=netabin;
}
      
void THBTprocessor::SetNPxBins(Int_t npxbin)const
{
//setter
MESH.n_px_bins=npxbin;
}
void THBTprocessor::SetNPyBins(Int_t npybin)const
{
//setter
MESH.n_py_bins=npybin;
}
void THBTprocessor::SetNPzBins(Int_t npzbin)const
{
//setter
MESH.n_pz_bins=npzbin;
}

void THBTprocessor::SetNBins1DFineMesh(Int_t n)const
{
//setter
MESH.n_1d_fine=n;
}
void THBTprocessor::SetBinSize1DFineMesh(Float_t x)const
{
//setter
MESH.binsize_1d_fine=x;
}
      
void THBTprocessor::SetNBins1DCoarseMesh(Int_t n)const
{
//setter
MESH.n_1d_coarse =n;
}
void THBTprocessor::SetBinSize1DCoarseMesh(Float_t x)const
{
//setter
MESH.binsize_1d_coarse=x;
}
      
void THBTprocessor::SetNBins3DFineMesh(Int_t n)const
{
//setter
MESH.n_3d_fine =n;
}
void THBTprocessor::SetBinSize3DFineMesh(Float_t x)const
{
//setter
MESH.binsize_3d_fine=x;
}
      
void THBTprocessor::SetNBins3DCoarseMesh(Int_t n)const
{
//setter
MESH.n_3d_coarse=n;
}
void THBTprocessor::SetBinSize3DCoarseMesh(Float_t x)const
{
//setter
MESH.binsize_3d_coarse=x;
}      
      
void THBTprocessor::SetNBins3DFineProjectMesh(Int_t n )const
{
//setter
MESH.n_3d_fine_project=n;
}

void THBTprocessor::SetPrintFull(Int_t flag) const 
{
 //setter
 PARAMETERS.print_full=flag;
} 
