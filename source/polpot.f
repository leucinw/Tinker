c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module polpot  --  polarization functional form details  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     poleps    induced dipole convergence criterion (rms Debyes/atom)
c     p2scale   scale factor for 1-2 polarization energy interactions
c     p3scale   scale factor for 1-3 polarization energy interactions
c     p4scale   scale factor for 1-4 polarization energy interactions
c     p5scale   scale factor for 1-5 polarization energy interactions
c     p41scale  additional factor for 1-4 intragroup polarization
c     d1scale   scale factor for intra-group direct induction
c     d2scale   scale factor for 1-2 group direct induction
c     d3scale   scale factor for 1-3 group direct induction
c     d4scale   scale factor for 1-4 group direct induction
c     u1scale   scale factor for intra-group mutual induction
c     u2scale   scale factor for 1-2 group mutual induction
c     u3scale   scale factor for 1-3 group mutual induction
c     u4scale   scale factor for 1-4 group mutual induction
c========for AMOEBA+=============
c     p12scale  scale factor for inter-group 1-2 polarization energy interactions
c     p13scale  scale factor for inter-group 1-3 polarization energy interactions
c     p14scale  scale factor for inter-group 1-4 polarization energy interactions
c     p15scale  scale factor for inter-group 1-5 polarization energy interactions
c     p21scale  scale factor for intra-group 1-2 polarization energy interactions
c     p31scale  scale factor for intra-group 1-3 polarization energy interactions
c     p41scale  scale factor for intra-group 1-4 polarization energy interactions
c     p51scale  scale factor for intra-group 1-5 polarization energy interactions
c==================================
c     udiag     acceleration factor for induced dipole SCF iterations
c     politer   maximum number of induced dipole SCF iterations
c     poltyp    type of polarization potential (direct or mutual)
c
c
      module polpot
      implicit none
      integer politer
      real*8 poleps
      real*8 p2scale,p3scale
      real*8 p4scale,p5scale
      real*8 d1scale,d2scale
      real*8 d3scale,d4scale
      real*8 u1scale,u2scale
      real*8 u3scale,u4scale
      real*8 p12scale,p13scale
      real*8 p14scale,p15scale
      real*8 p21scale,p31scale
      real*8 p41scale,p51scale
      real*8 udiag
      character*6 poltyp
      save
      end
