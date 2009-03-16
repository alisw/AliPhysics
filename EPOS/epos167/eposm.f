      subroutine NumberModel(cmodel,model)
      character cmodel*21
      stop'   ***** This program can only run EPOS *****      '
      end
      
      subroutine IniModel(model)      
      end      
      
      subroutine IniEvtModel 
      end  
      
      subroutine emsaaaModel(model,iret)
      end    

      function crseModel(model,ekin,maproj,matarg,idtarg)
      crseModel=0d0
      end	      
      
      subroutine m2XXFZ( a,b)
      double precision a, b(2)
      end	

      subroutine m3SIGMA(ek,idpro,idtar,latar,matar,sigi,sige)
      end

      subroutine m6SIGMA(icl,engy,stot,sela,sine,sdifr,slela,Rho)
      end	

      subroutine m7SIGMA(stot,scut,sine,slela)
      end	
