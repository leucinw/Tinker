c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kpolar  --  assign polarizability parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kpolar" assigns atomic dipole polarizabilities to the atoms
c     within the structure and processes any new or changed values
c
c     literature reference:
c
c     A. C. Simmonett, F. C. Pickard IV, J. W. Ponder and B. R. Brooks,
c     "An Empirical Extrapolation Scheme for Efficient Treatment of
c     Induced Dipoles", Journal of Chemical Physics, 145, 164101 (2016)
c     [OPT coefficients]
c
c
      subroutine kpolar
      use sizes
      use atoms
      use cflux
      use inform
      use iounit
      use keys
      use kpolr
      use mpole
      use neigh
      use polar
      use polpot
      use potent
      use usolve
      implicit none
      integer i,j,k,next
      integer nlist,npg
      integer pg(maxval)
      integer, allocatable :: list(:)
      real*8 pol,thl,dird
      real*8 sixth
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(polarity))  deallocate (polarity)
      allocate (polarity(n))
c
c     process keywords containing polarizability parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            k = 0
            pol = 0.0d0
            thl = -1.0d0
            dird = -1.0d0
            do j = 1, maxval
               pg(j) = 0
            end do
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=30,end=30)pol,thl,dird,(pg(j),j=1,maxval)
   30       continue
            if (k .gt. 0) then
               if (k .le. maxtyp) then
                  polr(k) = pol
                  athl(k) = thl
                  adird(k) = dird
               end if
            end if
         end if
      end do
c
c     find and store the atomic dipole polarizability parameters
c
      do i = 1, n
         polarity(i) = polr(type(i))
         thole(i) = athl(type(i))
         dirdamp(i) = adird(type(i))
      end do
c
c     process keywords containing atom specific polarizabilities
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            k = 0
            pol = 0.0d0
            thl = 0.0d0
            dird = 0.0d0
            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:240)
               read (string,*,err=80,end=80)  pol,thl,dird
   80          continue
               polarity(k) = pol
               thole(k) = thl
               dirdamp(k) = dird 
            end if
         end if
      end do
c
c     remove zero and undefined polarizable sites from the list
c
      npolar = 0
      if (use_polar) then
         npole = 0
         do i = 1, n
            if (polsiz(i).ne.0 .or. polarity(i).ne.0.0d0) then
               npole = npole + 1
               ipole(npole) = i
               pollist(i) = npole
               zaxis(npole) = zaxis(i)
               xaxis(npole) = xaxis(i)
               yaxis(npole) = yaxis(i)
               polaxe(npole) = polaxe(i)
               do k = 1, maxpole
                  pole(k,npole) = pole(k,i)
               end do
               if (polarity(i) .ne. 0.0d0)  npolar = npolar + 1
               polarity(npole) = polarity(i)
               thole(npole) = thole(i)
               dirdamp(npole) = dirdamp(i)
            end if
         end do
      end if
      return
      end
