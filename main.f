      program main
c     intrinsic functions
      double precision dabs, dacos, dble
      integer timenowhigh, timediff
c     variables
      double precision pi, reference, x
      integer starttime, endtime, duration, rank, p, n
c     functions
      double precision f1, f2, integrate
      
      starttime = timenowhigh()
      
      rank = myprocid()
      p = nprocs()
      n = 2
      
      pi = dacos(dble(-1.))
      reference = dble(0.)
      
      x = dble(1.)      
      reference = integrate(dble(0.), dble(1.), n, p, rank, 2)
      
      endtime = timenowhigh()
      duration = timediff(endtime, starttime)

      if (rank.eq.0) then
        print*, 'refernce pi=', pi
        print*, 'calculated pi=', reference
        print*, 'difference=', dabs(pi - reference)
        print*, 'execution time=', duration / 10e6, ' s'
        print*, 'n=', n
        print*, 'p=', p
      endif
            
      end
      
      double precision function integrate(a, b, n, p, rank, f)
        double precision a, b, h, local_a, local_b, integrate, param
        double precision f1, f2
        integer n, p, rank, f, local_n, i, mod
        
c       width of slices
        h = (b - a) / n
c       number of slices per processor
        local_n = n / p
        local_a = a + rank * local_n * h
        local_b = local_a + local_n * h
        
        if (f.eq.1) then
            integrate = dble(0.)
            do i=0, n, 1
                param = dble(a + dble(i) * h)
                if ((i.eq.0).or.(i.eq.n)) then
                    integrate = integrate + f1(param)
                else if (mod(i,2)) then
                    integrate = integrate + 4 * f1(param)
                else
                    integrate = integrate + 2 * f1(param)
                endif
            enddo
            integrate = integrate * h / dble(3.)
        endif
        if (f.eq.2) then
            integrate = dble(0.)
            do i=0, n, 1
                param = dble(a + dble(i) * h)
                if ((i.eq.0).or.(i.eq.n)) then
                    integrate = integrate + f2(param)
                else if (mod(i,2)) then
                    integrate = integrate + 4 * f2(param)
                else
                    integrate = integrate + 2 * f2(param)
                endif
            enddo
            integrate = integrate * h / dble(3.)
        endif
        return
      end
      
      double precision function f1(x)
        double precision dsin, x, f1
        f1 = x * dsin(x)        
        return    
      end
      
      double precision function f2(x)
        double precision x, f2
        f2 = (4 / (1 + x * x))
c        print*, 'f2=', 4 / (1 + x * x), ' x=', x
        return      
      end
      
      
