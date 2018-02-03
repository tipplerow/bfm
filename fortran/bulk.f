c ----------------------------------------------------------------------

      program BULK

c BULK: Stochastic lattice dynamics simulation of a system of linear,
c monodisperse polymer chains with excluded volume in a periodic box.

c In this implementation, the polymer chains are modeled with the
c ONE-SITE bond fluctuation model.  The stochastic dynamics is a simple
c dynamic Monte Carlo process.

c In our implementation of the ONE-SITE bond fluctuation model, each
c polymer bead occupies one site on the lattice.  The bonds between
c connected beads are restricted to the set (1,0,0), (1,1,0), and
c (1,1,1).  Excluded volume conditions are enforced for all polymer
c beads by restricting the occupancy of the lattice sites to at most one
c polymer bead.  An array of occupation numbers for the primary lattice,
c on which the beads reside, is maintained to check for bead overlaps.

c To forbid the crossing of polymer chains during the dynamics, we also
c track the MIDPOINTS of the polymer bonds.  These coordinates of bond
c midpoints are always integer multiples of 1/2, so we maintain an
c additional lattice of occupation numbers for bond midpoints.  This
c lattice has twice the number of entries in each coordinate direction
c as the lattice of occupation numbers for polymer beads.  Any bond
c crossing must involve the simultaneous overlap of the coordinates of
c two bond midpoints, so this bond crossing condition is easy to detect.

c Note that the occupation numbers for the beads must be either zero or
c one because the excluded volume condition is always enforced.  In
c contrast, the occupation numbers for the bond midpoints may be greater
c than one, because several bond midpoints can overlap when the chains
c are allowed to cross.

c The connectivity information is provided by the array KCONN, described
c below, so that polydisperse systems can be simulated with relative
c ease, if we choose do so later.

c PARAMETERS
c ----------
c
c     N................Number of beads per polymer chain
c     npoly............Number of polymer chains
c     nbead............Total number of polymer beads
c     nx,ny,nz.........Size of the lattice box for polymer segments
c     nxmid,
c     nymid,
c     nzmid............Size of the lattice box for bond midpoints
c
c     nrelax...........Estimate of the relaxation time
c     ninit............Number of MCS to be taken during initial
c                      equilibration period
c     kdata............Number of MCS between accumulation of data 
c     ndata............Number of data points that will be calculated.
c
c                      (The total number of MCS will be KDATA*NDATA.)
c
c     noncross.........Logical flag: 
c                        TRUE,  if the non-crossability dynamical
c                               constraint is enforced,
c                        FALSE, if it is ignored.
c 
c ----------

      implicit none

      integer          N, npoly, nbead, nx,ny,nz, nxmid,nymid,nzmid
      integer          nrelax, ninit, kdata, ndata
      logical          noncross

      parameter        (N     =  40)
      parameter        (npoly = 100)

      parameter        (nbead = N*npoly)

      parameter        (nx = 20)
      parameter        (ny = 20)
      parameter        (nz = 20)

      parameter        (nxmid = 2*nx)
      parameter        (nymid = 2*ny)
      parameter        (nzmid = 2*nz)

      parameter        (nrelax = 20000)
      parameter        (ninit  = 20*nrelax)
      parameter        (kdata  = 200)
      parameter        (ndata  = 1000)

      parameter        (noncross = .true.)


c VARIABLES
c ---------
c
c     lpoly............Lattice coordinates of all polymer beads
c     locc.............Lattice occupation numbers for polymer beads:
c
c                        LOCC(jx,jy,jz) = 0 if (jx,jy,jz) is unoccupied;
c                        LOCC(jx,jy,jz) = 1 if (jx,jy,jz) is occupied.
c
c     loccmid..........Lattice occupation numbers for bond midpoints:
c
c                        LOCCMID(jxm,jym,jzm) = 0 if no bond midpoint
c                        is present at position (jxm/2,jym/2,jzm/2).
c
c                        LOCCMID(jxm,jym,jzm) is a positive integer, M,
c                        if M bond midpoints are present at position
c                        (jxm/2,jym/2,jzm/2) on the lattice; 
c
c
c     kconn............Array containing the connectivity information: 
c
c                        For polymer beads not at chain ends, 
c                        KCONN(1,I) = I-1 and KCONN(2,I) = I+1, 
c                        signifying that the bead I is bonded to 
c                        its two neighbors,
c
c                        for the first bead of a polymer chain,
c                        KCONN(1,I) = 0 and KCONN(2,I) = I+1, 
c                        because there is no bond to the "left",
c
c                        for the last bead of a polymer chain,
c                        KCONN(1,I) = I-1 and KCONN(2,I) = 0, 
c                        because there is no bond to the "right".
c
c
c ---------------------

      integer          imc1, imc2, jpoly, iseed, itrash, iu, 
     :                 J, jbead, jx, jy, jz, jxmid, jymid, jzmid,
     :                 jxxx, jyyy, jzzz, kcoord, jcon2
      integer          lper
      integer          kconn(2,nbead)
      integer          lpoly(3,nbead)
      integer          locc(nx,ny,nz), loccmid(nxmid,nymid,nzmid)
      external         lper

      double precision clock

c ----------------------------------------------------------------------

c ==================
c RANDOM NUMBER SEED
c ==================

      iseed = -1

c ==========================
c READ INITIAL CONFIGURATION
c ==========================

      iu = 10
      open (unit = iu, file = 'init.coord', status = 'old')

      read (iu, *) itrash
      if (itrash .ne. N) then
         write (*, *) 'Number of beads does not match input'
         stop
      endif

      read (iu, *) itrash
      if (itrash .ne. npoly) then
         write (*, *) 'Number of polymers does not match input'
         stop
      endif

      do jbead = 1, nbead
         read (iu, *) jx, jy, jz

         lpoly(1,jbead) = jx
         lpoly(2,jbead) = jy
         lpoly(3,jbead) = jz
      enddo

      close (unit = iu)

c ============
c CONNECTIVITY
c ============

      do jpoly = 1, npoly
         J = 1
         jbead = J + N*(jpoly-1)

         kconn(1,jbead) = 0
         kconn(2,jbead) = jbead + 1

         do J = 2, (N-1)
            jbead = J + N*(jpoly-1)

            kconn(1,jbead) = jbead - 1
            kconn(2,jbead) = jbead + 1
         enddo

         J = N
         jbead = J + N*(jpoly-1)

         kconn(1,jbead) = jbead - 1
         kconn(2,jbead) = 0
      enddo

c ==============
c RUN PARAMETERS
c ==============

      iu = 9
      open (unit = iu, file = 'vital.stats', status = 'unknown')
      write (iu,  *) ' '
      write (iu, 55) '               Chain length: ', N
      write (iu, 55) '           Number of chains: ', npoly
      write (iu,  *) ' '
      write (iu, 60) '        Simulation box size: ', nx, ny, nz
      write (iu,  *) ' '
      write (iu, 50) '            Average density: ', 
     :               dble(nbead)/dble(nx*ny*nz)
      write (iu,  *) ' '
      write (iu, 55) ' MCS between data recording: ', kdata
      write (iu, 55) '     Number of data records: ', ndata
      write (iu, 55) '        Total number of MCS: ', kdata*ndata
      write (iu,  *) ' '

      close (unit = iu)

c =============================
c Initialize occupation numbers
c =============================

      write (*, *) ' '
      write (*, *) ' Initializing occupation numbers...'
      write (*, *) ' '

      do jz = 1, nz
         do jy = 1, ny
            do jx = 1, nx
               locc(jx,jy,jz) = 0
            enddo
         enddo
      enddo
      
      do jzmid = 1, nzmid
         do jymid = 1, nymid
            do jxmid = 1, nxmid
               loccmid(jxmid,jymid,jzmid) = 0
            enddo
         enddo
      enddo
      
      do jbead = 1, nbead
         jx = lpoly(1,jbead)
         jy = lpoly(2,jbead)
         jz = lpoly(3,jbead)

         jxxx = lper(jx,nx)
         jyyy = lper(jy,ny)
         jzzz = lper(jz,nz)

         if (locc(jxxx,jyyy,jzzz) .ne. 0) then
            write (*, *) 'ERROR: Bead overlap detected'
            stop
         endif

         locc(jxxx,jyyy,jzzz) = 1

         jcon2 = kconn(2,jbead)

         if (jcon2 .ne. 0) then
            jxmid = lpoly(1,jbead) + lpoly(1,jcon2)
            jymid = lpoly(2,jbead) + lpoly(2,jcon2)
            jzmid = lpoly(3,jbead) + lpoly(3,jcon2)

            jxxx = lper(jxmid,nxmid)
            jyyy = lper(jymid,nymid)
            jzzz = lper(jzmid,nzmid)

            if (loccmid(jxxx,jyyy,jzzz) .ne. 0) then
               write (*, *) 'WARNING: Bond crossing detected: ',
     :                      jxxx,jyyy,jzzz, jbead
            endif

c NOTE: There may be more than one bond midpoint at a given location... 

            loccmid(jxxx,jyyy,jzzz) = loccmid(jxxx,jyyy,jzzz) + 1
         endif
      enddo

c ====================
c EQUILIBRATION PERIOD
c ====================

      write (*, *) ' '
      write (*, *) ' Equilibrating...'
      write (*, *) ' '

      do imc1 = 1, ninit
         call advance ( nbead, nx, ny, nz, nxmid, nymid, nzmid,
     :                  noncross, kconn, lpoly, locc, loccmid, iseed )
      enddo

c Write equilibrated coordinates...

      open (unit = 77, file = 'equil.coord', status = 'unknown')
      do jbead = 1, nbead
         write (77, 100) lpoly(1,jbead),lpoly(2,jbead),lpoly(3,jbead)
      enddo
      close (unit = 77)

c -----------------------
c Check for bond crossing
c -----------------------

      if (noncross) then

         do jbead = 1, nbead
            jcon2 = kconn(2,jbead)

            if (jcon2 .ne. 0) then
               jxmid = lpoly(1,jbead) + lpoly(1,jcon2)
               jymid = lpoly(2,jbead) + lpoly(2,jcon2)
               jzmid = lpoly(3,jbead) + lpoly(3,jcon2)

               jxxx = lper(jxmid,nxmid)
               jyyy = lper(jymid,nymid)
               jzzz = lper(jzmid,nzmid)

               if (loccmid(jxxx,jyyy,jzzz) .ne. 1) then
                  write (*, *) 'ERROR: Bond crossing detected: ',
     :                         jxxx,jyyy,jzzz, loccmid(jxxx,jyyy,jzzz)
                  stop
               endif
            endif
         enddo

         write (*, *) ' '
         write (*, *) ' No bond crossings remain...'
         write (*, *) ' '

      endif

c =================
c DATA ACCUMULATION
c =================

      write (*, *) ' '
      write (*, *) ' The real thing...'
      write (*, *) ' '

      clock = 0.D0

      kcoord = 20
      open (unit = kcoord, file = 'coord.stream', status = 'unknown')

      do imc1 = 1, ndata
         do imc2 = 1, kdata
            call advance ( nbead, nx, ny, nz, nxmid, nymid, nzmid,
     :                     noncross, kconn, lpoly, locc, loccmid,
     :                     iseed )
            clock = clock + 1.D0
         enddo

         write (*, 10) clock

c Write coordinates to output stream:

         write (kcoord, 10) clock
         do jbead = 1, nbead
            write (kcoord, 100) lpoly(1,jbead), 
     :                          lpoly(2,jbead), 
     :                          lpoly(3,jbead)
         enddo

c Write re-start information:

         open (unit = 30, file = 'restart.coord', status = 'unknown')

         write (30, 90) N,     '  :  Chain length'
         write (30, 90) npoly, '  :  Number of chains'

         do jbead = 1, nbead
            write (30, 100) lpoly(1,jbead),
     :                      lpoly(2,jbead),
     :                      lpoly(3,jbead)
         enddo

         close (unit = 30)

         open (unit = 31, file = 'restart.clock', status = 'unknown')
         write (31, 10) clock
         close (unit = 31)

      enddo

      close (unit = kcoord)

      stop
   10 format (1x, f12.2)
   50 format (a, f8.4)
   55 format (a, i10)
   60 format (a, i4, 1x, i4, 1x, i4)
   90 format (1x, i6, a)
  100 format (1x, i4, 1x, i4, 1x, i4)
      end
c ----------------------------------------------------------------------

      subroutine advance ( nbead, nx, ny, nz, nxmid, nymid, nzmid,
     :                     noncross, kconn, lpoly, locc, loccmid, 
     :                     iseed )

c ADVANCE: This subroutine advances the dynamic Monte Carlo time
c evolution of a lattice polymer system by one MC step. A total of NBEAD
c time steps are attempted.  The polymer chains are modeled by the
c ONE-SITE bond fluctuation model.

c Periodic boundary conditions are applied in all coordinate directions.

c INPUT
c     nbead............Total number of polymer beads
c     nx,ny,nz.........Size of the lattice box for polymer beads
c     nxmid,
c     nymid,
c     nzmid............Size of the lattice box for bond midpoints
c     noncross.........Flag for enforcement of non-crossability
c     lpoly............Initial lattice coordinates of all polymer beads
c     locc.............Lattice occupation numbers for polymer beads:
c
c                        LOCC(jx,jy,jz) = 0 if (jx,jy,jz) is unoccupied;
c                        LOCC(jx,jy,jz) = 1 if (jx,jy,jz) is occupied.
c
c     loccmid..........Lattice occupation numbers for bond midpoints:
c
c                        LOCCMID(jxm,jym,jzm) = 0 if no bond midpoint
c                        is present at position (jxm/2,jym/2,jzm/2).
c
c                        LOCCMID(jxm,jym,jzm) is a positive integer, M,
c                        if M bond midpoints are present at position
c                        (jxm/2,jym/2,jzm/2) on the lattice; 
c
c     kconn............Array containing the connectivity information: 
c
c                        For polymer beads not at chain ends, 
c                        KCONN(1,I) = I-1 and KCONN(2,I) = I+1, 
c                        signifying that the bead I is bonded to 
c                        its two neighbors,
c
c                        for the first bead of a polymer chain,
c                        KCONN(1,I) = 0 and KCONN(2,I) = I+1, 
c                        because there is no bond to the "left",
c
c                        for the last bead of a polymer chain,
c                        KCONN(1,I) = I-1 and KCONN(2,I) = 0, 
c                        because there is no bond to the "right".
c
c
c     iseed............Initial random number seed
c
c OUTPUT
c     lpoly............Updated lattice coordinates of all beads
c     locc.............Updated occupation numbers for beads
c     loccmid..........Updated occupation numbers for bond midpoints
c     iseed............New random number seed
c
c ---------------------

      implicit none

c -- ARGUMENTS --------
      logical          noncross
      integer          nbead, nx, ny, nz, nxmid, nymid, nzmid, iseed
      integer          kconn(2,nbead)
      integer          lpoly(3,nbead)
      integer          locc(nx,ny,nz)
      integer          loccmid(nxmid,nymid,nzmid)

c -- LOCAL ------------
      integer          iatt, jbead, jmove, jcon1, jcon2, jx, jy, jz,
     :                 jxmid1, jymid1, jzmid1, jxmid2, jymid2, jzmid2, 
     :                 jxxx, jyyy, jzzz, LSQR

      integer          lper
      double precision ran3
      external         lper, ran3

c ----------------------------------------------------------------------

c =====================
c   -----------------
c     MOVE ATTEMPTS
c   -----------------
c =====================

      do iatt = 1, nbead

c ================
c Select the event
c ================

         jbead = 1 + int(dble(nbead)*ran3(iseed))
         jmove = 1 + int(6.D0*ran3(iseed))

c ===================
c Generate trial move
c ===================

         jx = lpoly(1,jbead)
         jy = lpoly(2,jbead)
         jz = lpoly(3,jbead)

         if (jmove .eq. 1) then
            jx = jx + 1
         else if (jmove .eq. 2) then
            jx = jx - 1
         else if (jmove .eq. 3) then
            jy = jy + 1
         else if (jmove .eq. 4) then
            jy = jy - 1
         else if (jmove .eq. 5) then
            jz = jz + 1
         else
            jz = jz - 1
         endif


         jcon1 = kconn(1,jbead)
         jcon2 = kconn(2,jbead)

         if (jcon1 .ne. 0) then
            jxmid1 = jx + lpoly(1,jcon1)
            jymid1 = jy + lpoly(2,jcon1)
            jzmid1 = jz + lpoly(3,jcon1)
         endif

         if (jcon2 .ne. 0) then
            jxmid2 = jx + lpoly(1,jcon2)
            jymid2 = jy + lpoly(2,jcon2)
            jzmid2 = jz + lpoly(3,jcon2)
         endif


c ===============
c Test trial move
c ===============

c As soon as we find a violation of a necessary condition for the
c acceptance of the move, we jump to the end of the loop over trials... 

c ---------------------------
c Excluded volume constraints
c ---------------------------

         jxxx = lper(jx,nx)
         jyyy = lper(jy,ny)
         jzzz = lper(jz,nz)

         if (locc(jxxx,jyyy,jzzz) .ne. 0) goto 1000

c --------------------
c Allowed bond vectors
c --------------------

c We have already checked for bond lengths of zero 
c by enforcing the excluded volume condition...

         if (jcon1 .ne. 0) then
            LSQR = (jx - lpoly(1,jcon1))**2
     :           + (jy - lpoly(2,jcon1))**2
     :           + (jz - lpoly(3,jcon1))**2

            if (LSQR .gt. 3) goto 1000
         endif

         if (jcon2 .ne. 0) then
            LSQR = (jx - lpoly(1,jcon2))**2
     :           + (jy - lpoly(2,jcon2))**2
     :           + (jz - lpoly(3,jcon2))**2

            if (LSQR .gt. 3) goto 1000
         endif

c -------------
c Bond crossing
c -------------

         if (noncross) then

            if (jcon1 .ne. 0) then
               jxxx = lper(jxmid1,nxmid)
               jyyy = lper(jymid1,nymid)
               jzzz = lper(jzmid1,nzmid)

               if (loccmid(jxxx,jyyy,jzzz) .ne. 0) goto 1000
            endif

            if (jcon2 .ne. 0) then
               jxxx = lper(jxmid2,nxmid)
               jyyy = lper(jymid2,nymid)
               jzzz = lper(jzmid2,nzmid)

               if (loccmid(jxxx,jyyy,jzzz) .ne. 0) goto 1000
            endif

         endif


c =============
c UPDATE SYSTEM
c =============

c If we get to this point, then the move has been accepted...

c Update occupation numbers:

         jxxx = lper(lpoly(1,jbead),nx)
         jyyy = lper(lpoly(2,jbead),ny)
         jzzz = lper(lpoly(3,jbead),nz)
         locc(jxxx,jyyy,jzzz) = 0

         jxxx = lper(jx,nx)
         jyyy = lper(jy,ny)
         jzzz = lper(jz,nz)
         locc(jxxx,jyyy,jzzz) = 1


c Update bond midpoints:

         if (jcon1 .ne. 0) then
            jxxx = lper(lpoly(1,jbead) + lpoly(1,jcon1), nxmid)
            jyyy = lper(lpoly(2,jbead) + lpoly(2,jcon1), nymid)
            jzzz = lper(lpoly(3,jbead) + lpoly(3,jcon1), nzmid)
            loccmid(jxxx,jyyy,jzzz) = loccmid(jxxx,jyyy,jzzz) - 1

            jxxx = lper(jxmid1,nxmid)
            jyyy = lper(jymid1,nymid)
            jzzz = lper(jzmid1,nzmid)
            loccmid(jxxx,jyyy,jzzz) = loccmid(jxxx,jyyy,jzzz) + 1
         endif

         if (jcon2 .ne. 0) then
            jxxx = lper(lpoly(1,jbead) + lpoly(1,jcon2), nxmid)
            jyyy = lper(lpoly(2,jbead) + lpoly(2,jcon2), nymid)
            jzzz = lper(lpoly(3,jbead) + lpoly(3,jcon2), nzmid)
            loccmid(jxxx,jyyy,jzzz) = loccmid(jxxx,jyyy,jzzz) - 1

            jxxx = lper(jxmid2,nxmid)
            jyyy = lper(jymid2,nymid)
            jzzz = lper(jzmid2,nzmid)
            loccmid(jxxx,jyyy,jzzz) = loccmid(jxxx,jyyy,jzzz) + 1
         endif


c Update coordinates:

         lpoly(1,jbead) = jx
         lpoly(2,jbead) = jy
         lpoly(3,jbead) = jz

 1000 enddo

      return
      end
c ----------------------------------------------------------------------

      integer function lper ( lcoord, n )

c LPER: Given an uncorrected lattice coordinate LCOORD, and the size of
c the lattice (in the corresponding coordinate direction) N, this
c function returns the corrected lattice coordinates by applying
c periodic boundary conditions...

c Here we assume that the primary lattice box runs from 1 through N.

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

c RAN3: Returns a uniform random deviate on [0,1).  Set IDUM to any
c negative value to initialize or reinitialize the sequence.  
c
c From "Numerical Recipes", Chapter 7.

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

      save             iff, inext, inextp
      save             ma

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
   11    continue
         do 13 k=1,4
            do 12 i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.mz)ma(i)=ma(i)+mbig
   12       continue
   13    continue
         inext=0
         inextp=31
         idum=1
      endif

      inext = inext + 1
      if (inext .eq. 56) inext = 1
      inextp = inextp + 1
      if (inextp .eq. 56) inextp = 1

      mj = ma(inext) - ma(inextp)
      if (mj .lt. mz) mj = mj + mbig
      ma(inext) = mj

      ran3 = mj*fac

      return
      end

