class PMD_dig_reco {

public:
     Int_t ncl[5001],icl[331] ;
     Int_t iorder [3] [5185],infocl [3] [73] [73],inford [4] [5185] ;
     Int_t incr,clno ,nmx1;
     Int_t module1 [31] [73] [73] ;
     Int_t module2 [31] [73] [73] ;
     Float_t u[98], c__, cd, cm ;
     Float_t det1 [73] [73],det2 [73] [73] ;
     Float_t coord [3] [73] [73],pts [3] [26],wts [26] ;
     Float_t r0[501], x0[501], y0[501] , rr[501] ;
     Float_t h0[501];
     Float_t x[501], y[501], z__[501];
     Float_t xx[501], yy[501], zz[501],hh[501];
};

