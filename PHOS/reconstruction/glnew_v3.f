
	subroutine clrgl
	include 'comgl.for'
	include 'comalwl.for'
	include 'comwlgen.for'

	ipgl=0
	ipcl=0
	ntclst=0
	do nw=1,nwall
	nglw(nw)=0
	ncls(nw)=0
	enddo
	return
	end

	subroutine wlfn(iwl)
	include 'comgl.for'
	include 'comalwl.for'
	include 'comwlgen.for'

	miwl=1
	mawl=nwall
	if(iwl.gt.0.and.iwl.le.nwall) then
		miwl=iwl
		mawl=iwl
	endif
	do iw=miwl,mawl
	ip=ipwgl(iw)+1
	ng=nglw(iw)
c	ipwal=madwl(iw)+1
	call wlf(wl,ng,mgl(ip),egl(ip))
	enddo
	return
	end

	subroutine cfgln(iwl)
	include 'comgl.for'
	include 'comalwl.for'
	include 'comwlgen.for'

	miwl=1
	mawl=nwall
	if(iwl.gt.0.and.iwl.le.nwall) then
		miwl=iwl
		mawl=iwl
	endif
	do iw=miwl,mawl
	ip=ipwgl(iw)+1
	ng=nglw(iw)
c	ipwal=madwl(iw)+1
	call glcfml(cf,ng,mgl(ip),egl(ip))
	enddo
	return
	end


	subroutine clrwln(iwl)
	include 'comgl.for'
	include 'comalwl.for'
	include 'comwlgen.for'

	miwl=1
	mawl=nwall
	if(iwl.gt.0.and.iwl.le.nwall) then
		miwl=iwl
		mawl=iwl
	endif
	do iw=miwl,mawl
	ip=ipwgl(iw)+1
	ng=nglw(iw)
c	ipwal=madwl(iw)+1
	call clrwal(wl,ng,mgl(ip))
	enddo
	return
	end


        SUBROUTINE CLRWAL(W,NG,MG)
        DIMENSION W(*),MG(*)

        IF(NG.EQ.0) RETURN
        DO L=1,NG
        M=MG(L)
        W(M)=0.
        ENDDO
        RETURN
        END

      SUBROUTINE WLF(W,NG,MG,EG)
      DIMENSION W(*),MG(*),EG(*)

        IF(NG.EQ.0) RETURN
        DO  L=1,NG
        M=MG(L)
        W(M)=EG(L)
        ENDDO
        RETURN
        END

       SUBROUTINE GLCFML(cf,NG,MG,EG)
       DIMENSION MG(*),EG(*),cf(*)

        IF(NG.EQ.0) RETURN
        DO  L=1,NG
        M=MG(L)
        EG(L)=EG(L)*CF(M)
	enddo
        RETURN
        END

