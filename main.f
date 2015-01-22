      program main

c     int(x * sin(x), 0, pi)
      call execute(1)

c     int(4 / (1 + x^2), 0, 1)
      call execute(2)

      end

      subroutine execute(f)
c     intrinsic functions
      double precision dabs, dacos, dble
      integer timenowhigh, timediff, addnewlink, newtop
c     variables
      double precision pi, calculated, rresult
      integer starttime, endtime, duration, rank, p, n, topid,
     / link(7), link1, f
c     functions
      double precision integrate

c     machine pi
      pi = dacos(dble(-1.))

c     number of intervals
      n = 8000
      
c     get information about processes
      rank = myprocid()
      p = nprocs()

c     create star topology
      topid = newtop(p-1)
      if (p.gt.1) then
        if (rank.eq.0) then
          do i = 1, p-1
            link(i) = addnewlink(topid, i, 1)
          enddo
        else
          link1 = addnewlink(topid, 0, 1)
        endif
      endif

      starttime = timenowhigh()

c     calculate local portion
      calculated = 0
      if (f.eq.1) then
        calculated = integrate(dble(0.), pi, n, p, rank, f)
      endif
      if (f.eq.2) then
        calculated = integrate(dble(0.), dble(1.), n, p, rank, f)
      endif

c     recieve and sum up all portions on first process
      if (p.gt.1) then
        if (rank.eq.0) then
          do i = 1, p-1
            call recv(topid, link(i), rresult, 8)
            calculated = calculated + rresult
          enddo
        else
          call send(topid, link1, calculated, 8)
        endif
      endif

      endtime = timenowhigh()
      duration = timediff(endtime, starttime)

      if (rank.eq.0) then
        print*, 'n=', n
        print*, 'p=', p
        print*, 'refernce pi=', pi
c        print*, 'calculated pi (f', f, ')=', calculated
        write(*, 900) 'calculated pi (f', f, ')=', calculated
900     format (X1 A I1 A x2 F17.15)
        print*, 'difference=', dabs(pi - calculated)
        print*, 'execution time=', duration / 10e6, ' s'
      endif

      call freetop(topid)

      end

      double precision function integrate(a, b, n, p, rank, f)
      double precision a, b, h, local_a, local_b, integrate, param
      double precision f1, f2
      integer n, p, rank, f, local_n, i, mod

c     width of slices
      h = (b - a) / dble(n)
c     number of slices per processor (rounded up)
      local_n = (n + p - 1) / p

      if (mod(local_n, 2).ne.0) then
        print*, 'number of slices must be even, exiting!'
        stop 'end of run'
        return
      endif

      local_a = a + dble(rank) * dble(local_n) * h
c     last processor has less work to do
      if ((rank.eq.(p - 1)).and.(mod(n, local_n).ne.0)) then
        local_n = mod(n, local_n)
      endif
      local_b = local_a + dble(local_n) * h

      if (f.eq.1) then
        integrate = dble(0.)
        do i=0, local_n
          param = dble(local_a + dble(i) * h)
          if ((i.eq.0).or.(i.eq.local_n)) then
            integrate = integrate + f1(param)
          else if (mod(i, 2)) then
            integrate = integrate + 4 * f1(param)
          else
            integrate = integrate + 2 * f1(param)
          endif
        enddo
        integrate = integrate * h / dble(3.)
      endif

      if (f.eq.2) then
        integrate = dble(0.)
        do i=0, local_n
          param = dble(local_a + dble(i) * h)
          if ((i.eq.0).or.(i.eq.local_n)) then
            integrate = integrate + f2(param)
          else if (mod(i, 2)) then
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
