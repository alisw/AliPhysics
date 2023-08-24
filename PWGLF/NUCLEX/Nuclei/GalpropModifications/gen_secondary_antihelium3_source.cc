// serksnyte


using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galprop_internal.h"

//#include <fort_interface.h>

#include <Processes_Interface.h>

#include <cstring>
#include <sstream>
#include <fstream>

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galprop::gen_secondary_antihelium3_source(Particle &particle)
{
   INFO("ENTRY");

   ostringstream buf;
   buf<<"Generating "<<particle.name<<" source function for n_spatial_dimensions="
       <<gcr[0].n_spatial_dimensions;
   INFO(buf.str());

   if ("secondary_antihelium3" != particle.name) 
   {  
      buf.str("");
      buf<<"invalid particle "<<particle.name; 
      WARNING(buf.str());
      return 2; 
   }
   if (galdef.secondary_antihelium3 != 1  ) 
   {
      ERROR("option for secondary_antihelium3 is not known");
      ERROR("Available options: 1");
      return 3;
   }

   int stat=0, iprotons=-1, iHelium =-1,  Z1, A1, Z2, A2;

   //Store the cross sections and indices in vectors
   //An index of -1 means use protons distribution
   //The first value in all vectors should always exist and the corresponding index should be to protons.
   std::vector<double> cs_HI, cs_He;
   std::vector<int> ind_HI, ind_He;
   cs_HI.reserve(n_species);
   cs_He.reserve(n_species);
   ind_HI.reserve(n_species);
   ind_He.reserve(n_species);

   Distribution protons;                 // IMOS20000606.6

// identify CR protons                   // IMOS20000606.7
   if(galdef.n_spatial_dimensions==2) protons.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
   if(galdef.n_spatial_dimensions==3) protons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
   protons=0.;
   for(int i=0; i<n_species; i++)  
      if(101==100*gcr[i].Z+gcr[i].A)
      {
         iprotons=i;
       protons+=gcr[iprotons].cr_density;

         buf.str("");
         buf<<"  CR protons found as species #"<<iprotons;
         INFO(buf.str());
      }
   if(iprotons==-1) { ERROR("CR protons not found!"); return 1; }

// identify CR Helium
   for(int i=0; i<n_species; i++) if(204 == 100*gcr[i].Z+gcr[i].A) iHelium =i;
   if(iHelium ==-1) { ERROR("CR Helium  not found!"); return 1; }
   else if(galdef.verbose>=1) {
      buf.str("");
      buf<<"  CR Helium  found as species #"<<iHelium;
      INFO(buf.str());
   }


// The production cross section of antihelium 3 in p-p collisions from Shukla et al  https://doi.org/10.1103/PhysRevD.102.063004

ifstream fileAnirvan;
fileAnirvan.open("antihelium3CSAnirvan.txt");
double AnirvanCSMatrix[27][27];

for(int i = 0; i < 27; ++i)
{
   for(int j = 0; j < 27; ++j)
   {
      fileAnirvan >> AnirvanCSMatrix[i][j];
   }
}

fileAnirvan.close();

//Gulli20070821 
   for(int ip_sec=0; ip_sec<particle.n_pgrid; ip_sec++)
   {
      for(int ip=0; ip<gcr[iprotons].n_pgrid; ip++)
      {
         cs_HI.clear();
         cs_He.clear();
         ind_HI.clear();
         ind_He.clear();
         if (galdef.secondary_antihelium3 == 1) {
            double AnirvanCs = AnirvanCSMatrix[ip_sec][ip]; // given as dsigma/dekin in cm² c/GeV 
            double antihelium3_cs = AnirvanCs*1.e24*particle.p[ip_sec]/particle.Etot[ip_sec]; // transform to dsigma/dp, and from cm² to barn       

            cs_HI.push_back(antihelium3_cs);
            ind_HI.push_back(-1); //Use the protons

            cs_He.push_back(antihelium3_cs*pow(4., 2.2/3.)); //beam+target: p+He
            ind_He.push_back(-1);

            cs_HI.push_back(antihelium3_cs*pow(4., 2.2/3.)); // beam+target: He+HI
            ind_HI.push_back(iHelium);

            cs_He.push_back(antihelium3_cs*pow(16., 2.2/3.)); // beam+target: He+He
            ind_He.push_back(iHelium);
         } 

         if(galaxy.n_spatial_dimensions==2)
         {
            for(int ir=0; ir<gcr[iprotons].n_rgrid; ir++)
            {
               for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
               {
                  const double gas = galaxy.n_HI.d2[ir][iz].s[0]+2.0*galaxy.n_H2.d2[ir][iz].s[0]+galaxy.n_HII.d2[ir][iz].s[0];
                  particle.secondary_source_function.d2[ir][iz].s[ip_sec ] +=  
                     gas * (cs_HI[0] + galdef.He_H_ratio * cs_He[0]) * 
                     protons.d2[ir][iz].s[ip] * gcr[iprotons].Ekin[ip];
                  for (size_t ics=1; ics < cs_HI.size(); ++ics) {
                     particle.secondary_source_function.d2[ir][iz].s[ip_sec ]+=  
                        gas * cs_HI[ics] * gcr[ind_HI[ics]].cr_density.d2[ir][iz].s[ip] * 
                        gcr[ind_HI[ics]].Ekin[ip] * gcr[ind_HI[ics]].A;
                  }
                  for (size_t ics=1; ics < cs_He.size(); ++ics) {
                     particle.secondary_source_function.d2[ir][iz].s[ip_sec ]+=  
                        gas * galdef.He_H_ratio * cs_He[ics] * gcr[ind_He[ics]].cr_density.d2[ir][iz].s[ip] * 
                        gcr[ind_He[ics]].Ekin[ip] * gcr[ind_He[ics]].A;
                  }
                  //std::cout<<"END He gas "<<std::endl;
                  
                        
               }  //  iz
            }  //  ir
         }  //  particle.n_spatial_dimensions==2
         //std::cout<<max<<std::endl;
         if(galaxy.n_spatial_dimensions==3)
         {
            for(int ix=0; ix<gcr[iprotons].n_xgrid; ix++)
            {
               for(int iy=0; iy<gcr[iprotons].n_ygrid; iy++)
               {
                  for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
                  {
                     const double gas = galaxy.n_HI.d3[ix][iy][iz].s[0]+2.0*galaxy.n_H2.d3[ix][iy][iz].s[0]+galaxy.n_HII.d3[ix][iy][iz].s[0];

                     //First is always protons and always exists
                     particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ] +=  
                        gas * (cs_HI[0] + galdef.He_H_ratio * cs_He[0]) * 
                        protons.d3[ix][iy][iz].s[ip] * gcr[iprotons].Ekin[ip];

                     //Loop over the rest
                     for (size_t ics=1; ics < cs_HI.size(); ++ics) {
                        particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ]+=  
                           gas * cs_HI[ics] * gcr[ind_HI[ics]].cr_density.d3[ix][iy][iz].s[ip] * 
                           gcr[ind_HI[ics]].Ekin[ip] * gcr[ind_HI[ics]].A;
                     }
                     for (size_t ics=1; ics < cs_He.size(); ++ics) {
                        particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ]+=  
                           gas * galdef.He_H_ratio * cs_He[ics] * gcr[ind_He[ics]].cr_density.d3[ix][iy][iz].s[ip] * 
                           gcr[ind_He[ics]].Ekin[ip] * gcr[ind_HI[ics]].A;
                     }
                  }  //  iz
               }  //  iy
            }  //  ix
         }  //  particle.n_spatial_dimensions==3
      }  //  ip
   }  //  ip_sec
   
   const double factor=1.e-24 *1.e-3 *C *log(galdef.Ekin_factor); // transformation to cm2/MeV and constant factors
   particle.secondary_source_function *= factor;

   protons.delete_array();                  // IMOS20000606.10

   if(galdef.verbose>=2)
   {
      cout<<"   particle.secondary_source_function for "<<particle.name<<endl;
      particle.secondary_source_function.print();
    }
   INFO("Exit");
   return stat;
}
