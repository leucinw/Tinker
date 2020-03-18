c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ect0  --  charge transfer energy                 ##
c     ##                                                              ##
c     ##################################################################
c
c
c
      subroutine ect0
      use sizes
      use atoms
      use energi
      use limits
      use warp
      implicit none
c
c     choose the method for summing over pairwise interactions
c
      if (use_ctlist) then
         call ect0b
      else
         call ect0a
      end if
      return
      end 
      
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ect0  --  charge transfer energy via double loop ##
c     ##                                                              ##
c     ##################################################################

      subroutine ect0a
      use sizes
      use ctran
      use atomid
      use atoms
      use bound
      use cell
      use couple
      use energi
      use group
      use mutant
      use polgrp
      use polpot
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,it
      integer kk,kt
      real*8 e
      real*8 apre,bexp
      real*8 aprei,bexpi
      real*8 aprek,bexpk
      real*8 fgrp
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8, allocatable :: ctscale(:)
      logical proceed,usei
      logical muti,mutk
      character*6 mode
c
c
c     zero out the charge transfer energy contribution
c
c
      ect = 0.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (ctscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         ctscale(i) = 1.0d0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'CT'
      call switch (mode)
c
c     find the charge transfer energy via double loop search
c
      do ii = 1, nct-1
         i = ict(ii)
         it = jct(i)
         usei = use(i)
         muti = mut(i) 
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            ctscale(i12(j,i)) = p12scale
            do k = 1, np11(i)
               if (i12(j,i) .eq. ip11(k,i))
     &            ctscale(i12(j,i)) = p21scale
            end do
         end do
         do j = 1, n13(i)
            ctscale(i13(j,i)) = p13scale
            do k = 1, np11(i)
               if (i13(j,i) .eq. ip11(k,i))
     &            ctscale(i13(j,i)) = p31scale
            end do
         end do
         do j = 1, n14(i)
            ctscale(i14(j,i)) = p14scale
            do k = 1, np11(i)
               if (i14(j,i) .eq. ip11(k,i))
     &            ctscale(i14(j,i)) = p41scale
            end do
         end do
         do j = 1, n15(i)
            ctscale(i15(j,i)) = p15scale
            do k = 1, np11(i)
               if (i15(j,i) .eq. ip11(k,i))
     &            ctscale(i15(j,i)) = p51scale
            end do
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nct
            k = ict(kk)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jct(k)
               xr = x(i) - x(k)
               yr = y(i) - y(k)
               zr = z(i) - z(k)
               call image (xr,yr,zr)

               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               aprei = abs(aprmct(it))
               aprek = abs(aprmct(kt))
               bexpi = abs(bprmct(it))
               bexpk = abs(bprmct(kt))
               apre = sqrt(aprei*aprek)
               bexp = 0.5d0*(bexpi + bexpk)
               apre = apre*ctscale(k)
               if ((muti .and. .not.mutk) .or.
     &            (mutk .and. .not.muti)) then
                  apre = apre * elambda  
               endif

               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  e = -apre*1000.0d0*exp(-bexp*rik)

c
c     use energy switching if near the cutoff distance
c
               if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
               end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overallcharge transfer energy component
c
                     ect = ect + e
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            ctscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            ctscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            ctscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            ctscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nct
         i = ict(ii)
         it = jct(i)
         usei = use(i)
         muti = mut(i) 
c
c     set interaction scaling coefficients for connected atoms
c

         do j = 1, n12(i)
            ctscale(i12(j,i)) = p12scale
            do k = 1, np11(i)
               if (i12(j,i) .eq. ip11(k,i))
     &            ctscale(i12(j,i)) = p21scale
            end do
         end do
         do j = 1, n13(i)
            ctscale(i13(j,i)) = p13scale
            do k = 1, np11(i)
               if (i13(j,i) .eq. ip11(k,i))
     &            ctscale(i13(j,i)) = p31scale
            end do
         end do
         do j = 1, n14(i)
            ctscale(i14(j,i)) = p14scale
            do k = 1, np11(i)
               if (i14(j,i) .eq. ip11(k,i))
     &            ctscale(i14(j,i)) = p41scale
            end do
         end do
         do j = 1, n15(i)
            ctscale(i15(j,i)) = p15scale
            do k = 1, np11(i)
               if (i15(j,i) .eq. ip11(k,i))
     &            ctscale(i15(j,i)) = p51scale
            end do
         end do
c
c    decide whether to compute the current interaction
c
         do kk = ii, nct
            k = ict(kk)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jct(k)
               do j = 1, ncell
                  xr = x(i) - x(k)
                  yr = y(i) - y(k)
                  zr = z(i) - z(k)
                  call imager (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
              
                  aprei = abs(aprmct(it))
                  aprek = abs(aprmct(kt))
                  bexpi = abs(bprmct(it))
                  bexpk = abs(bprmct(kt))
              
                  
                  apre = sqrt(aprei*aprek)
                  bexp = 0.5d0*(bexpi + bexpk)

                  apre = apre*ctscale(k)
                  if ((muti .and. .not.mutk) .or.
     &               (mutk .and. .not.muti)) then
                     apre = apre * elambda  
                  endif

                  if (rik2 .le. off2) then
                     rik = sqrt(rik2)
                     e = -apre*1000.0d0*exp(-bexp*rik)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                        e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge transfer energy component;
c     interaction of an atom with its own image counts half
c
                     if (i .eq. k) e = 0.5d0 * e
                     ect = ect + e
                  end if
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            ctscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            ctscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            ctscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            ctscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (ctscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ect0b  --  Charge transfer energy via list    ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ect0b" calculates the charge transfer energy
c     using a pairwise neighbor list
c
c
      subroutine ect0b
      use sizes
      use atomid
      use atoms
      use bound
      use couple
      use ctran
      use energi
      use group
      use mutant
      use neigh
      use polgrp
      use polpot
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,it
      integer kk,kt
      real*8 e
      real*8 apre,bexp
      real*8 aprei,bexpi
      real*8 aprek,bexpk
      real*8 fgrp
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8, allocatable :: ctscale(:)
      logical proceed,usei
      logical muti, mutk
      character*6 mode
c
c
c     zero out the CT energy contribution
c
      ect = 0.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (ctscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         ctscale(i) = 1.0d0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'CT'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nct,ict,
!$OMP& jct,use,x,y,z,nctlst,ctlst,n12,n13,n14,n15,
!$OMP& i12,i13,i14,i15,p12scale,p13scale,p14scale,
!$OMP& p15scale,p21scale,p31scale,p41scale,p51scale,
!$OMP& use_group,off2,aprmct,bprmct,mut,elambda,
!$OMP& cut2,c0,c1,c2,c3,c4,c5) firstprivate(ctscale)
!$OMP& shared(ect)
!$OMP DO reduction(+:ect) schedule(guided)
c
c     find the van der Waals energy via neighbor list search
c
      do ii = 1, nct
         i = ict(ii)
         it = jct(i)
         usei = use(i)
         muti = mut(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            ctscale(i12(j,i)) = p12scale
            do k = 1, np11(i)
               if (i12(j,i) .eq. ip11(k,i))
     &            ctscale(i12(j,i)) = p21scale
            end do
         end do
         do j = 1, n13(i)
            ctscale(i13(j,i)) = p13scale
            do k = 1, np11(i)
               if (i13(j,i) .eq. ip11(k,i))
     &            ctscale(i13(j,i)) = p31scale
            end do
         end do
         do j = 1, n14(i)
            ctscale(i14(j,i)) = p14scale
            do k = 1, np11(i)
               if (i14(j,i) .eq. ip11(k,i))
     &            ctscale(i14(j,i)) = p41scale
            end do
         end do
         do j = 1, n15(i)
            ctscale(i15(j,i)) = p15scale
            do k = 1, np11(i)
               if (i15(j,i) .eq. ip11(k,i))
     &            ctscale(i15(j,i)) = p51scale
            end do
         end do
c
c     decide whether to compute the current interaction
c
         do kk = 1, nctlst(ii)
            k = ict(ctlst(kk,ii))
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jct(k)
               xr = x(i) - x(k)
               yr = y(i) - y(k)
               zr = z(i) - z(k)
               call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               aprei = abs(aprmct(it))
               aprek = abs(aprmct(kt))
               bexpi = abs(bprmct(it))
               bexpk = abs(bprmct(kt))

               !alternative is to use mixing rules 
               apre = sqrt(aprei*aprek)
               bexp = 0.5d0*(bexpi + bexpk)

               apre = apre*ctscale(k)
               if ((muti .and. .not.mutk) .or.
     &            (mutk .and. .not.muti)) then
                  apre = apre * elambda  
               endif

               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  e = -apre*1000.0d0*exp(-bexp*rik)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall CT energy components
c
                  ect = ect + e
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            ctscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            ctscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            ctscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            ctscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (ctscale)
      return
      end
