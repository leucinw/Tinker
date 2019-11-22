c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kcp  --  charge pene   parameter assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kcp" assigns parameters to any additional user defined
c
c
      subroutine kcp
      use atoms
      use atomid
      use keys
      use chgpen
      implicit none
      integer i,k,next
      real*8 pen 
      character*20 keyword
      character*120 record
      character*120 string
c
c     check for keyword for charge penetration form
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:3) .eq. 'CP ') then
            call getnumb (record,k,next)
            if (k.ge.1 .and. k.le.maxtyp) then
               pen = penalpha(k)
               string = record(next:120)
               read (string,*,err=10,end=10) pen
  10           continue
               penalpha(k) = pen 
            end if
         end if
      end do
      return
      end
