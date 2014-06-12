#ifndef _f_Init_included_
#define _f_Init_included_

namespace Photospp
{

const static int NMXHEP = 10000;
const static double PI    = 3.14159265358979324;
const static double TWOPI = 6.28318530717958648;


// extern "C"
//{

	/** Definition of the PHOEVT common block */
	extern struct HEPEVT
	{
		int    nevhep;
		int    nhep;
		int    isthep[NMXHEP];
		int    idhep[NMXHEP];
		int    jmohep[NMXHEP][2];
		int    jdahep[NMXHEP][2];
		double phep[NMXHEP][5];
		double vhep[NMXHEP][4];
	        int    qedrad[NMXHEP];  // to be bool once compatibility with F77 removed       
	        const static int nmxhep=NMXHEP;
		//      NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
		//  JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
	        //   int qedrad[NMXHEP]  was an add up 
                //   for  HEPEVT in F77 times. Separate common PH_PHOQED
	        //   also phoif_.chkif[NMXPHO] was add up for PHOEVT
                //   now it is pho.qedrad
	} hep,pho;
        //ph_hepevt_,phoevt_;



	extern struct PHOCOP
	{
		double alpha;
		double xphcut;
	} phocop_;

	extern struct PHNUM
	{
		double iev;
	} phnum_;



	extern struct PHOKEY
	{
		double fsec;
		double fint;
		double expeps;
		int interf;
		int isec;
		int itre;
		int iexp;
		int iftop;
		int ifw;
	} phokey_;

	extern struct PHOSTA
	{
		int status[10];
	        int ifstop;
	} phosta_;

	extern struct PHOLUN
	{
		int phlun;
	} pholun_;

	extern struct PHOPHS
	{
		double xphmax;
		double xphoto;
		double costhg;
		double sinthg;

	} phophs_;
	extern struct TOFROM
	{
		double QQ[4];
		double XM;
		double th1;
		double fi1;

	} tofrom_;

	extern struct PHOPRO
	{
		double probh;
		double corwt;
		double xf;
		int irep;
	} phopro_;

	extern struct PHOREST
	{
		double fi3;
		double fi1;
		double th1;
        	int irep;     //    provably line to be removed

	} phorest_;

	extern struct PHWT
	{
		double beta;
		double wt1;
		double wt2;
		double wt3;

	} phwt_;
	extern struct PHOCORWT
	{
		double phocorwt3;
		double phocorwt2;
		double phocorwt1;

	} phocorwt_;

	extern struct PHOMOM
	{
		double mchsqr;
		double mnesqr;
		double pneutr[5];
	} phomom_;
	extern struct PHOCMS
	{
		double bet[3];
		double gam;
	} phocms_;

	extern struct PHOEXP
	{
	        const static int NX = 10;
		double pro[NX];
		int nchan;
	        int expini;    // bool
	} phoexp_;

	//debug mode on if ipoin <  1 and ipoinm > 1
	extern struct PHLUPY
	{
		int ipoin;
		int ipoinm;
	} phlupy_;

	/** Initialize kinematic corrections */
	void PHCORK(int modcor);

	/** Single branch processing */
	void PHOTOS_MAKE_C(int id);

	/* Central management routine. Defines what action
	   will be performed at point ID. */
	void PHTYPE(int ID);
// }

} // namespace Photospp
#endif
