      program main
c     intrinsic functions
      double precision dabs, dacos, dble
      integer timenowhigh, timediff, addnewlink, newtop
c     variables
      double precision pi, referencef1, referencef2, rresult
      integer starttime, endtime, duration, rank, p, n, topid,
     / link(7), link1
c     functions
      double precision integrate

      rank = myprocid()
      p = nprocs()      
      
     
      topid = newtop(p-1)
      
      if(p.gt.1) then 
        if(rank.eq.0) then
           do i = 1,p-1,1    
               link(i) = addnewlink(topid,i,1)
           enddo
        else
           link1 = addnewlink(topid,0,1)        
        endif
      endif
      

      n = 8000
      referencef1 = dble(0.)
      referencef2 = dble(0.)
      
      pi = dacos(dble(-1.))
      
      if (rank.eq.0) then
        print*, 'n=', n
        print*, 'p=', p
      endif
      
c     x*sin(x)      
      starttime = timenowhigh()
            
      referencef1 = integrate(dble(0.), pi, n, p, rank, 1)
      
      if(p.gt.1) then 
        if(rank.eq.0) then
c           print*, 'ref', referencef1
           do i = 1,p-1,1                   
               call recv(topid, link(i), rresult, 8)
c               print*, 'rec', rresult
               referencef1 = referencef1 + rresult
           enddo        
        else
           call send(topid, link1, referencef1, 8)       
        endif
      endif
  
      endtime = timenowhigh()
      duration = timediff(endtime, starttime)

      if (rank.eq.0) then
        print*, 'refernce pi=', pi
        print*, 'calculated pi (f1)=', referencef1
        print*, 'difference=', dabs(pi - referencef1)
        print*, 'execution time=', duration / 10e6, ' s'
      endif
      
c     4/(1+x^2)
      starttime = timenowhigh()
         
      referencef2 = integrate(dble(0.), dble(1.), n, p, rank, 2)
      
      if(p.gt.1) then
        if(rank.eq.0) then
c           print*, 'ref', referencef2
           do i = 1,p-1,1    
               call recv(topid, link(i), rresult, 8)
c               print*, 'rec', rresult
               referencef2 = referencef2 + rresult
           enddo        
        else
           call send(topid, link1, referencef2, 8)       
        endif
      endif
      
      endtime = timenowhigh()
      duration = timediff(endtime, starttime)

      if (rank.eq.0) then
        print*, 'refernce pi=', pi
        print*, 'calculated pi (f2)=', referencef2
        print*, 'difference=', dabs(pi - referencef2)
        print*, 'execution time=', duration / 10e6, ' s'
      endif
      
c      print*, integrate(dble(0.), dble(1.), n, 3, 0, 2) + integrate(
c     \ dble(0.), dble(1.), n, 3, 1, 2) + integrate(dble(0.),
c     \ dble(1.), n, 3, 2, 2)
      
      call freetop(topid)
      end
      
      
      
      double precision function integrate(a, b, n, p, rank, f)
        double precision a, b, h, local_a, local_b, integrate, param
        double precision f1, f2
        integer n, p, rank, f, local_n, i, mod
        
c       width of slices
        h = (b - a) / dble(n)
c       number of slices per processor (rounded up)
        local_n = (n + p - 1)/ p
            
        local_a = a + dble(rank) * dble(local_n) * h
c       last processor has less work to do
        if ((rank.eq.(p - 1)).and.(mod(n, local_n).ne.0)) then
            local_n = mod(n, local_n)
        endif
        local_b = local_a + dble(local_n) * h
       
        if (f.eq.1) then
            integrate = dble(0.)
            do i=0, local_n, 1
                param = dble(local_a + dble(i) * h)                
                if ((i.eq.0).or.(i.eq.local_n)) then
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
            do i=0, local_n, 1
                param = dble(local_a + dble(i) * h)
                if ((i.eq.0).or.(i.eq.local_n)) then
                    integrate = integrate + f2(param)
                else if (mod(i,2)) then
                    integrate = integrate + dble(4.) * f2(param)
                else
                    integrate = integrate + dble(2.) * f2(param)
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
        return      
      end
      
      
