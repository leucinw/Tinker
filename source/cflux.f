c
c     ##################################################################
c     ##                                                              ##
c     ##  module cflux -- charge flux term forcefield parameters      ##
c     ##                                                              ##
c     ##################################################################
c
c
c     jb        charge flux over unit bond length
c     b0        equilibrium bond length in charge flux
c     ja        charge flux over unit angle 
c     a0        equilibrium angle in charge flux
c
      module cflux
      use sizes
      implicit none
      real*8, allocatable :: theta0(:)
      real*8, allocatable :: bp0(:,:) 
      real*8, allocatable :: jbp(:,:)
      real*8, allocatable :: jtheta(:,:)
      real*8, allocatable :: jb(:)
      real*8, allocatable :: b0(:)
      real*8, allocatable :: pchrgflux(:)
      logical dobond,doangle 
      save
      end
