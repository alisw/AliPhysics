
	subroutine extrgw(nt,w,ng,m,e,eth)
	real w(10),e(10)
	integer m(10)
	ng=0
	do i=1,nt
	if(w(i).gt.eth) then
		ng=ng+1
		m(ng)=i
		e(ng)=w(i)
	else
		w(i)=0.
	endif
	enddo
	return
	end

