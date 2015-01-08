      program main
c     intrinsic functions
      double precision dabs, dacos, dble
c     variables
      double precision pi, referenz, x
c     functions
      f1, f2, integrate
      
      pi = dacos(dble(-1.))
      referenz = dble(0.)
      
      x = dble(1.)      
      result = f1(x)
      print*, result
      
      
      print*, pi, referenz, dabs(pi - referenz)
      
            
      end
      
      double precision function integrate(a,b,n)
        double precision a, b
        integer n
        return
      end
      
      double precision function f1(x)
        double precision dsin, x
        f1 = x*dsin(x)        
        return    
      end
      
      double precision function f2(x)
        double precision x
        f2 = 4/(1+x*x)  
        return      
      end
      
      
