      subroutine hbcdf(y,n,ia)
      character *132 string
      dimension ia(n)
      character*12 form
      write(form,'(''(F'',i3.3,''.3)'')') n
      write(string(1:n),form) y
      do k=n,n-2,-1
         if(string(k:k).eq.'0') then
            string(k:k)=' '
            if(k.eq.n-2) string(k-1:k-1)=' '
         endif
      enddo
      call uctoh1(string,ia,n)
      end
