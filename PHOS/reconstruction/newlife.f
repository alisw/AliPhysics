
	function mmconvert(mmm)
	mmconvert=mmm	
	return
	end	

	function newtipcell(mmnew)
	include 'comarray.for'
	if(mmnew.lt.0) then
		write (*,*) ' newtipcell: STRANGE mmnew ',mmnew
		return
	endif
	idc=mmcells(mmnew)
	newtipcell=idcelmat(idc)
	return
	end
