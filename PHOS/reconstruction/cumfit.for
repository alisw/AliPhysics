	function cumfit(x)
c fit for ALICE 94 RUN DATA
c	data p1,p2 /2.526, 1.104/
c94	data p1,p2 /2.300, 2.044/
c crystals for my geant
c	data p1,p2 /1.480, 1.600/
c 0rad
	data p1,p2 /1.700, 1.150/
c 0.5rad +
c	data p1,p2 /1.700, 0.800/
c 0.5rad -
c	data p1,p2 /5.500, 0.700/
c 0.3rad +
c	data p1,p2 /1.700, 0.900/
c 0.3rad -
c	data p1,p2 /4.00, 0.700/
c 0.2rad +
c	data p1,p2 /1.700, 1.000/
c 0.2rad -
c	data p1,p2 /3.50, 0.800/
c 0.1rad +
c	data p1,p2 /1.700, 1.100/
c 0.1rad -
c	data p1,p2 /2.50, 0.900/
c        data angl /0./
        data angl /-0.2/

        aangl=abs(angl)
        if(angl.ge.0.) then
                pp1=p1
                pp2=p2-0.7*aangl
        else
                pp1=p1+7.6*aangl
                pp2=p2-0.9*aangl
        endif
	x=abs(x)
c lead glass for my geant
c	x=x*0.55

	a=0.5*exp((pp1-sqrt(pp1**2+4.*pp2*x))/2./pp2)
	if(x.ge.0.) then
		cumfit=a
	else
		cumfit=1.-a
	endif
	end
