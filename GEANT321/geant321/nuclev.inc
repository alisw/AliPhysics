*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:41  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_NUCLEV_INC
#define CERNLIB_GEANT321_NUCLEV_INC
*
*
* nuclev.inc
*
*
*=== nuclev ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     NUClear LEVel common:                                            *
*                                                                      *
*     Created on  10 january 1992  by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 26-mar-93     by    Alfredo Ferrari               *
*                                                                      *
*     Included in the following routines:                              *
*                                                                      *
*                        BERTTP                                        *
*                        BIMSEL                                        *
*                        NCLVIN                                        *
*                        NCLVGT                                        *
*                        NCLVST                                        *
*                        NUCNUC                                        *
*                        NWISEL                                        *
*                        PFNCLV                                        *
*                        PIOABS                                        *
*                                                                      *
*     Description of the variables (nuclev):                           *
*                                                                      *
*                   paenuc(k,i) = pairing energies (GeV), for i=1 pro- *
*                                 ton pairing energy for Z=k, for i=2  *
*                                 neutron pairing energy for N=k       *
*                   shenuc(k,i) = shell energies (GeV), for i=1 pro-   *
*                                 ton shell energy for Z=k, for i=2    *
*                                 neutron shell energy for N=k         *
*                     defrmi(i) = level spacing at the Fermi level ac- *
*                                 cording to a Fermi gas               *
*                     defmag(i) = extra energy needed to overcome the  *
*                                 gap between the closed shell and the *
*                                 next level for nuclei with magic     *
*                                 numbers of protons or neutrons       *
*                   cumrad(k,i) = sum for j=1,k of 1 /Rj where Ri is   *
*                                 the maximum radius at which can be   *
*                                 found the jth nucleon of type i      *
*                                 (i=1 proton, i=2 neutron)            *
*                        ruspro = maximum radius at which protons can  *
*                                 still be found                       *
*                        rusneu = maximum radius at which neutrons can *
*                                 still be found                       *
*                        juspro = index for flagging that the jth      *
*                                 proton has been already used if      *
*                                 = inuclv                             *
*                        jusneu = index for flagging that the jth      *
*                                 neutron has been already used if     *
*                                 = inuclv                             *
*                        ntapro = number of protons in the target      *
*                        ntaneu = number of neutrons in the target     *
*                        navpro = number of available proton levels in *
*                                 the target                           *
*                        navneu = number of available neutron levels   *
*                                 in the target                        *
*                        nlspro = number of protons in the last level  *
*                        nlsneu = number of neutrons in the last level *
*                        ncopro = number of protons in the core        *
*                        nconeu = number of neutrons in the core       *
*                        nskpro = number of protons in the skin        *
*                        nskneu = number of neutrons in the skin       *
*                        nhapro = number of protons in the halo        *
*                        nhaneu = number of neutrons in the halo       *
*                        nuspro = number of "used" protons             *
*                        nusneu = number of "used" neutrons            *
*                        jmxpro = index  of the highest energy avai-   *
*                                 lable proton                         *
*                        jmxneu = index  of the highest energy avai-   *
*                                 lable neutron                        *
*                        nprnuc = number of presently used nucleons    *
*                                 (max 2)                              *
*                        iprnuc = type of the presently considered     *
*                                 nucleon                              *
*                        jprnuc = index of the presently considered    *
*                                 nucleon                              *
*                        magnum = magic numbers                        *
*                        magpro = last (possibly not filled) proton    *
*                                 shell                                *
*                        magneu = last (possibly not filled) neutron   *
*                                 shell                                *
*                     mgspro(k) = status of the kth proton shell:      *
*                                 -2 not filled but paired             *
*                                 -1 not filled and unpaired           *
*                                  0 filled (paired of course)         *
*                                 +1 originally  filled but now unpai- *
*                                    red                               *
*                     mgsneu(k) = status of the kth neutron shell      *
*                        inuclv = index of the present event           *
*                                                                      *
*----------------------------------------------------------------------*
*
*  NUClear LEVel
      COMMON / FKNULE / PAENUC (200,2), SHENUC (200,2), DEFRMI (2),
     &                  DEFMAG (2), ENNCLV (160,2), RANCLV (160,2),
     &                  CUMRAD (0:160,2), RUSNUC (2),
     &                  ENPLVL (114), ENNLVL(164), JUSNUC (160,2),
     &                  NTANUC (2), NAVNUC (2), NLSNUC (2), NCONUC (2),
     &                  NSKNUC (2), NHANUC (2), NUSNUC (2), JMXNUC (2),
     &                  IPRNUC (2), JPRNUC (2), MAGNUM (8), MAGNUC (2),
     &                  MGSNUC (8,2), MGSSNC (25,2), NSBSHL (2),
     &                  NPRNUC, INUCLV
      DIMENSION JUSPRO (160), JUSNEU (160), MGSPRO (8), MGSNEU (8),
     &          MGSSPR (19) , MGSSNE (25)
      EQUIVALENCE ( RUSNUC (1), RUSPRO )
      EQUIVALENCE ( RUSNUC (2), RUSNEU )
      EQUIVALENCE ( JUSNUC (1,1), JUSPRO (1) )
      EQUIVALENCE ( JUSNUC (1,2), JUSNEU (1) )
      EQUIVALENCE ( MGSNUC (1,1), MGSPRO (1) )
      EQUIVALENCE ( MGSNUC (1,2), MGSNEU (1) )
      EQUIVALENCE ( MGSSNC (1,1), MGSSPR (1) )
      EQUIVALENCE ( MGSSNC (1,2), MGSSNE (1) )
      EQUIVALENCE ( NTANUC (1), NTAPRO )
      EQUIVALENCE ( NTANUC (2), NTANEU )
      EQUIVALENCE ( NAVNUC (1), NAVPRO )
      EQUIVALENCE ( NAVNUC (2), NAVNEU )
      EQUIVALENCE ( NLSNUC (1), NLSPRO )
      EQUIVALENCE ( NLSNUC (2), NLSNEU )
      EQUIVALENCE ( NCONUC (1), NCOPRO )
      EQUIVALENCE ( NCONUC (2), NCONEU )
      EQUIVALENCE ( NSKNUC (1), NSKPRO )
      EQUIVALENCE ( NSKNUC (2), NSKNEU )
      EQUIVALENCE ( NHANUC (1), NHAPRO )
      EQUIVALENCE ( NHANUC (2), NHANEU )
      EQUIVALENCE ( NUSNUC (1), NUSPRO )
      EQUIVALENCE ( NUSNUC (2), NUSNEU )
      EQUIVALENCE ( JMXNUC (1), JMXPRO )
      EQUIVALENCE ( JMXNUC (2), JMXNEU )
      EQUIVALENCE ( MAGNUC (1), MAGPRO )
      EQUIVALENCE ( MAGNUC (2), MAGNEU )
 

#endif
