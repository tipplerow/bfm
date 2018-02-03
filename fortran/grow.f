c ----------------------------------------------------------------------

      program GROW2

c GROW2: This program grows a system of monodisperse three-dimensional
c lattice polymer chains, with excluded volume, for use in the ONE-SITE
c bond-fluctuation model.  Periodic boundary conditions are applied in
c all coordinate directions.  

c One bead of these lattice polymers occupies a single lattice site.
c The set of bond vectors along the chain are randomly generated from
c the permutations of (1,0,0) and (1,1,0).  [Bond vectors from the set
c (1,1,1) are avoided in the initial configurations for reasons
c discussed in JCP 103, We do not necessarily avoid
c generating configurations in which two bonds may cross.  The dynamics
c to which these initial chains are subject must naturally relax so that
c these bond crossings and not allow further bond crossings to occur.
c This is easily done as explained in the dynamic code.
c
c See Carmesin and Kremer, Macromolecules 21, 2819 (1988), for some
c details about this model in two-dimensions.
c
c We grow all chains simultaneously...

      implicit none

      integer          N, npoly, nbead, nx, ny, nz
      integer          maxtry

      parameter        (N     =  40)
      parameter        (npoly = 100)
      parameter        (nbead = N*npoly)

      parameter        (nx = 20)
      parameter        (ny = 20)
      parameter        (nz = 20)

      parameter        (maxtry = 1000)

      logical          success
      integer          lpoly(3,nbead)
      integer          locc(nx,ny,nz)
      integer          iseed, itry, iu, J, jbead, jpoly

c ----------------------------------------------------------------------

c ==================
c RANDOM NUMBER SEED
c ==================

      iseed = -1

c ==============
c INITIALIZATION
c ==============

      itry = 1
      success = .false.

      do while ((.not. success) .and. (itry .le. maxtry))

         write (*, *) ' *** ATTEMPT ', itry, ' ***'

         itry = itry + 1

         call grow_all ( N, npoly, nbead, nx, ny, nz, lpoly, locc, 
     :                   iseed, success )

      enddo

c =================
c WRITE COORDINATES
c =================

      iu = 20
      open (unit = iu, file = 'init.coord', status = 'unknown')

      write (iu, 90) N,     '  :  Chain length    '
      write (iu, 90) npoly, '  :  Number of chains'

      do jpoly = 1, npoly
         do J = 1, N
            jbead = J + N*(jpoly-1)

            write (iu,100) lpoly(1,jbead),lpoly(2,jbead),lpoly(3,jbead)
         enddo
      enddo

      close (unit = iu)

      stop
   90 format (1x, i6, a)
  100 format (1x, i4, 1x, i4, 1x, i4)
      end
c ----------------------------------------------------------------------

      subroutine grow_all ( N, npoly, nbead, nx, ny, nz, lpoly, locc, 
     :                      iseed, success )

c GROW_ALL: This subroutine grows a system of monodisperse,
c three-dimensional lattice polymer chains, with excluded volume, for
c use in the ONE-SITE bond-fluctuation model.  Periodic boundary
c conditions are applied in all coordinate directions.

c One bead of these lattice polymers occupies a single lattice site.
c The set of bond vectors along the chain are randomly generated from
c the permutations of (1,0,0) and (1,1,0).  We do not necessarily avoid
c generating configurations in which two bonds may cross.  The dynamics
c to which these initial chains are subject must naturally relax so that
c these bond crossings and not allow further bond crossings to occur.
c This is easily done as explained in the dynamic code.
c
c See Carmesin and Kremer, Macromolecules 21, 2819 (1988), for some
c details about this model in two-dimensions.
c
c In this subroutine, we grow all chains simultaneously...
c
c INPUT
c     N................Chain length: number of beads per polymer
c     npoly............Number of polymers
c     nbead............Total number of polymer beads
c     nx,ny,nz.........Size of the periodic lattice box
c     iseed............Random number seed
c
c OUTPUT
c     lpoly............Lattice coordinates of all polymer beads
c     locc.............Occupation numbers for all lattice sites
c     success..........Outcome of the growth
c
c ---------------------

      implicit none

c -- ARGUMENTS --------
      logical          success
      integer          N, npoly, nbead, nx, ny, nz, iseed
      integer          lpoly(3,nbead)
      integer          locc(nx,ny,nz)

c -- LOCAL ------------
      integer          maxplace
      parameter        (maxplace = 1000)

      logical          goodplace
      integer          J, jbead, jpoly, jx, jy, jz, jxxx, jyyy, jzzz,
     :                 kplace, lper
      integer          lbond(3)
      double precision ran3
      external         ran3, lper

c ---------------------

c ==========
c INITIALIZE
c ==========

      success = .true.

      do jz = 1, nz
         do jy = 1, ny
            do jx = 1, nx
               locc(jx,jy,jz) = 0
            enddo
         enddo
      enddo

c ==================
c PLACE CHAIN STARTS
c ==================

C ****      write (*, *) '    Placing chain starts...'

      J = 1
      jpoly = 1

      do while (success .and. (jpoly .le. npoly))
         jbead = J + N*(jpoly-1)

         kplace = 1
         goodplace = .false.

         do while ((.not. goodplace) .and. (kplace .le. maxplace))

c Put the position of the chain start at a random location in the box...

            jx = 1 + int(dble(nx)*ran3(iseed))
            jy = 1 + int(dble(ny)*ran3(iseed))
            jz = 1 + int(dble(nz)*ran3(iseed))

c Check excluded volume...

            if (locc(jx,jy,jz) .eq. 0) then
               goodplace = .true.

               lpoly(1,jbead) = jx
               lpoly(2,jbead) = jy
               lpoly(3,jbead) = jz

               locc(jx,jy,jz) = jbead
            endif

            kplace = kplace + 1
         enddo

         if (.not. goodplace) then
            write (*, *) ' COULD NOT PLACE CHAIN START ', jpoly
            success = .false.
         endif
         
         jpoly = jpoly + 1
      enddo

c =====================
c LOOP OVER OTHER BEADS
c =====================

      J = 2
      
      do while (success .and. (J .le. N))

C ****         write (*, *) '    Placing segment ', J, '...'

         jpoly = 1

         do while (success .and. (jpoly .le. npoly))
            jbead = J + N*(jpoly-1)

            kplace = 1
            goodplace = .false.

            do while ((.not. goodplace) .and. (kplace .le. maxplace))

c Add a bond vector chosen at random...

               call genbond2 ( lbond, iseed )
               
               jx = lpoly(1,jbead-1) + lbond(1)
               jy = lpoly(2,jbead-1) + lbond(2)
               jz = lpoly(3,jbead-1) + lbond(3)

c Check excluded volume...

               jxxx = lper(jx,nx)
               jyyy = lper(jy,ny)
               jzzz = lper(jz,nz)

               if (locc(jxxx,jyyy,jzzz) .eq. 0) then
                  goodplace = .true.
                     
                  lpoly(1,jbead) = jx
                  lpoly(2,jbead) = jy
                  lpoly(3,jbead) = jz
                  
                  locc(jxxx,jyyy,jzzz) = jbead
               endif

               kplace = kplace + 1
            enddo

            if (.not. goodplace) then
               write (*, *) 'COULD NOT PLACE SEGMENT ', J, 
     :                      ' OF CHAIN ', jpoly
               success = .false.
            endif

            jpoly = jpoly + 1
         enddo

         J = J + 1
      enddo

      return
      end
c ----------------------------------------------------------------------

      subroutine genbond2 ( lbond, iseed )

c GENBOND: This subroutine randomly generates a bond vector for a
c three-dimensional lattice polymer in the ONE-SITE bond-fluctuation
c model.  This routine chooses bonds from the
c following set of bond vectors:
c
c    (1,0,0),  (1,1,0).
c
c Note that the statistical weight of the vector set (1,0,0) is 6/18 and
c the weight of set (1,1,0) is 12/18.
c
c INPUT
c     iseed............Random number seed
c
c OUTPUT
c     lbond............Randomly generated bond vector
c     iseed............Random number seed
c
c ---------------------

      implicit none

c -- ARGUMENTS --------
      integer          lbond(3)
      integer*4        iseed

c -- LOCAL ------------
      double precision P1
      parameter        (P1 =  6.D0/18.D0)

      integer          idx, idy, idz, n1, n2, n3, nx, ny, nz
      double precision rannum1, rannum2, ran3
      external         ran3

c ----------------------------------------------------------------------

c =================
c Pick the bond set
c =================

      rannum1 = ran3(iseed)

      if (rannum1 .le. P1) then
         n1 = 1
         n2 = 0
         n3 = 0
      else 
         n1 = 1
         n2 = 1
         n3 = 0
      endif

c =========================
c Permute the displacements
c =========================

      rannum1 = ran3(iseed)
      rannum2 = ran3(iseed)

      if (rannum1 .lt. 0.333333333333333D0) then
         nx = n1
         if (rannum2 .lt. 0.5D0) then
            ny = n2
            nz = n3
         else
            nz = n2
            ny = n3
         endif
      else if (rannum1 .lt. 0.666666666666667D0) then
         ny = n1
         if (rannum2 .lt. 0.5D0) then
            nx = n2
            nz = n3
         else
            nz = n2
            nx = n3
         endif
      else
         nz = n1
         if (rannum2 .lt. 0.5D0) then
            ny = n2
            nx = n3
         else
            nx = n2
            ny = n3
         endif
      endif

c ======================
c Sign the displacements
c ======================

      if (ran3(iseed) .lt. 0.5D0) then
         idx =  nx
      else
         idx = -nx
      endif

      if (ran3(iseed) .lt. 0.5D0) then
         idy =  ny
      else
         idy = -ny
      endif

      if (ran3(iseed) .lt. 0.5D0) then
         idz =  nz
      else
         idz = -nz
      endif

c ========================
c Assemble the bond vector
c ========================

      lbond(1) = idx
      lbond(2) = idy
      lbond(3) = idz

      return
      end
c ----------------------------------------------------------------------

      integer function lper ( lcoord, n )

c LPER: Given an uncorrected lattice coordinate LCOORD, and the size of
c the lattice (in the corresponding coordinate direction) N, this
c function returns the corrected lattice coordinates by applying
c periodic boundary conditions...

      implicit none

c -- ARGUMENTS --------
      integer          lcoord, n

c ---------------------

      if (lcoord .lt. 1) then
         lper = n + mod(lcoord,n)
      else
         lper = 1 + mod(lcoord-1,n)
      endif

      return
      end
c ----------------------------------------------------------------------

      double precision function ran3 (idum)

c RAN3: Returns a uniform random deviate between 0.D0 and 1.D0.  Set
c IDUM to any negative value to initialize or reinitialize the sequence.
c
c See "Numerical Recipes", Chapter 7.

      implicit none

      integer          idum

      double precision mbig, mseed, mz, fac
      parameter (mbig =  4000000.D0)
      parameter (mseed = 1618033.D0)
      parameter (mz = 0.D0)
      parameter (fac = 2.5D-07)

      integer          i, ii, iff, inext, inextp, k
      double precision mj, mk
      double precision ma(55)

      save inext,inextp,ma

      data iff /0/

      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      return
      end

