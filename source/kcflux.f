c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kbond  --  bond stretch parameter assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kbond" assigns a force constant and ideal bond length
c     to each bond in the structure and processes any new or
c     changed parameter values
c
c
      subroutine kcflux
      use sizes
      use atoms
      use angbnd
      use atomid
      use bndstr
      use cflux 
      use inform
      use iounit
      use kcfluxes
      use keys
      use potent
      use usage
      implicit none
      integer i,j
      integer ia,ib,ic,ita,itb,itc
      integer na,nb
      integer size,next
      real*8 fc,bd
      real*8 tta,jtt1,jtt2,bd10,bd20,jb1,jb2
      character*4 pa,pb,pc
      character*8 blank2,pt2
      character*12 blank3,pt3
      character*20 keyword
      character*240 record
      character*240 string
c
c     process keywords containing bond stretch parameters
c
      blank2 = '        '
      blank3 = '            '
      size = 4
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'CFLUX-B ') then
            ia = 0
            ib = 0
            bd = 0.0d0
            fc = 0.0d0
            dobond = .true.
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,bd,fc
   10       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt2 = pa//pb
            else
               pt2 = pb//pa
            end if
            do j = 1, maxnbcf
               if (kcfb(j).eq.blank2 .or. kcfb(j).eq.pt2) then
                  kcfb(j) = pt2
                  jbnd(j) = fc
                  beq(j) = bd
                  goto 15
               end if
            end do
   15       continue
         else if (keyword(1:8) .eq. 'CFLUX-A ') then
            ia = 0
            ib = 0
            ic = 0
            tta = 0.0d0
            jtt1 = 0.0d0
            jtt2 = 0.0d0
            bd10 = 0.0d0
            bd20 = 0.0d0
            jb1 = 0.0d0
            jb2 = 0.0d0
            doangle = .true.
            string = record(next:240)
            read (string,*,err=20,end=20)  ia,ib,ic,tta,jtt1,jtt2,
     &          bd10,jb1,bd20,jb2
  20        continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt3 = pa//pb//pc
            else
               pt3 = pc//pb//pa
            end if
            do j = 1, maxnacf
               if (kcfa(j).eq.blank3 .or. kcfa(j).eq.pt3) then
                  kcfa(j) = pt3
                  theta0l(j) = tta
                  jthetal(1,j) = jtt1 
                  jthetal(2,j) = jtt2 
                  jbpl(1,j) = jb1 
                  jbpl(2,j) = jb2 
                  bp0l(1,j) = bd10 
                  bp0l(2,j) = bd20
                  goto 25 
               end if
            end do
   25       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nb = maxnbcf
      do i = maxnbcf, 1, -1
         if (kcfb(i) .eq. blank2)  nb = i - 1
      end do
c
c     determine the total number of forcefield parameters
c
      na = maxnacf
      do i = maxnacf, 1, -1
         if (kcfa(i) .eq. blank3)  na = i - 1
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(b0))  deallocate (b0)
      if (allocated(jb))  deallocate (jb)
      allocate (b0(nbond))
      allocate (jb(nbond))
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(theta0))  deallocate (theta0)
      if (allocated(bp0))  deallocate (bp0)
      if (allocated(jbp))  deallocate (jbp)
      if (allocated(jtheta))  deallocate (jtheta)
      allocate (theta0(nangle))
      allocate (bp0(2,nangle))
      allocate (jbp(2,nangle))
      allocate (jtheta(2,nangle))
c
c     assign ideal bond length and force constant for each bond
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ita = class(ia)
         itb = class(ib)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            pt2 = pa//pb
         else
            pt2 = pb//pa
         end if
         b0(i) = 0.0d0
         jb(i) = 0.0d0
c
c     assign bond stretching parameters for each bond
c
         do j = 1, nb
           if (kcfb(j) .eq. pt2) then
              b0(i) = beq(j)
              jb(i) = jbnd(j)
           end if
         end do
      end do
c
c    assign ideal bond angle and force constant for each angle
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pt3 = pa//pb//pc
         else
            pt3 = pc//pb//pa
         end if
         theta0(i) = 0.0d0
         bp0(1,i) = 0.0d0
         bp0(2,i) = 0.0d0
         jbp(1,i) = 0.0d0
         jbp(2,i) = 0.0d0
         jtheta(1,i) = 0.0d0
         jtheta(2,i) = 0.0d0
         do j = 1, na
           if (kcfa(j) .eq. pt3) then
             theta0(i) = theta0l(j) 
             bp0(1,i) = bp0l(1,j) 
             bp0(2,i) = bp0l(2,j)
             jbp(1,i) = jbpl(1,j) 
             jbp(2,i) = jbpl(2,j)
             jtheta(1,i) =jthetal(1,j)  
             jtheta(2,i) =jthetal(2,j)  
           end if
         end do
      end do
c
c     turn off the charge flux if bond and angle are not used 
c
      if (nbond .eq. 0 .and. nangle .eq. 0)  use_cflux = .false.
      return
      end
