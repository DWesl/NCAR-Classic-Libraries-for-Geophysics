c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                  copyright (c) 1998 by UCAR                   *
c     *                                                               *
c     *       University Corporation for Atmospheric Research         *
c     *                                                               *
c     *                      all rights reserved                      *
c     *                                                               *
c     *                      SPHEREPACK version 3.2                   *
c     *                                                               *
c     *       A Package of Fortran77 Subroutines and Programs         *
c     *                                                               *
c     *              for Modeling Geophysical Processes               *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *                  John Adams and Paul Swarztrauber             *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the National Center for Atmospheric Research          *
c     *                                                               *
c     *                Boulder, Colorado  (80307)  U.S.A.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the National Science Foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c
c     file alf.f contains subroutines alfk,lfim,lfim1,lfin,lfin1,lfpt
c     for computing normalized associated legendre polynomials
c
c subroutine alfk (n,m,cp)
c
c dimension of           real cp(n/2 + 1)
c arguments
c
c purpose                routine alfk computes single precision fourier
c                        coefficients in the trigonometric series
c                        representation of the normalized associated
c                        legendre function pbar(n,m,theta) for use by
c                        routines lfp and lfpt in calculating single
c                        precision pbar(n,m,theta).
c
c                        first define the normalized associated
c                        legendre functions
c
c                        pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)
c                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
c                        factorial(n)) times the (n+m)th derivative of
c                        (x**2-1)**n with respect to x=cos(theta)
c
c                        where theta is colatitude.
c
c                        then subroutine alfk computes the coefficients
c                        cp(k) in the following trigonometric
c                        expansion of pbar(m,n,theta).
c
c                        1) for n even and m even, pbar(m,n,theta) =
c                           .5*cp(1) plus the sum from k=1 to k=n/2
c                           of cp(k+1)*cos(2*k*th)
c
c                        2) for n even and m odd, pbar(m,n,theta) =
c                           the sum from k=1 to k=n/2 of
c                           cp(k)*sin(2*k*th)
c
c                        3) for n odd and m even, pbar(m,n,theta) =
c                           the sum from k=1 to k=(n+1)/2 of
c                           cp(k)*cos((2*k-1)*th)
c
c                        4) for n odd and m odd,  pbar(m,n,theta) =
c                           the sum from k=1 to k=(n+1)/2 of
c                           cp(k)*sin((2*k-1)*th)
c
c
c usage                  call alfk(n,m,cp)
c
c arguments
c
c on input               n
c                          nonnegative integer specifying the degree of
c                          pbar(n,m,theta)
c
c                        m
c                          is the order of pbar(n,m,theta). m can be
c                          any integer however cp is computed such that
c                          pbar(n,m,theta) = 0 if abs(m) is greater
c                          than n and pbar(n,m,theta) = (-1)**m*
c                          pbar(n,-m,theta) for negative m.
c
c on output              cp
c                          single precision array of length (n/2)+1
c                          which contains the fourier coefficients in
c                          the trigonometric series representation of
c                          pbar(n,m,theta)
c
c
c special conditions     none
c
c precision              single
c
c algorithm              the highest order coefficient is determined in
c                        closed form and the remainig coefficients are
c                        determined as the solution of a backward
c                        recurrence relation.
c
c accuracy               comparison between routines alfk and double
c                        precision dalfk on the cray1 indicates
c                        greater accuracy for smaller values
c                        of input parameter n.  agreement to 14
c                        places was obtained for n=10 and to 13
c                        places for n=100.
c
      subroutine alfk (n,m,cp)
      dimension       cp(n/2+1)
      parameter (sc10=1024.)
      parameter (sc20=sc10*sc10)
      parameter (sc40=sc20*sc20)
c
      cp(1) = 0.
      ma = iabs(m)
      if(ma .gt. n) return
      if(n-1) 2,3,5
    2 cp(1) = sqrt(2.)
      return
    3 if(ma .ne. 0) go to 4
      cp(1) = sqrt(1.5)
      return
    4 cp(1) = sqrt(.75)
      if(m .eq. -1) cp(1) = -cp(1)
      return
    5 if(mod(n+ma,2) .ne. 0) go to 10
      nmms2 = (n-ma)/2
      fnum = n+ma+1
      fnmh = n-ma+1
      pm1 = 1.
      go to 15
   10 nmms2 = (n-ma-1)/2
      fnum = n+ma+2
      fnmh = n-ma+2
      pm1 = -1.
 15   t1 = 1./sc20
      nex = 20
      fden = 2.
      if(nmms2 .lt. 1) go to 20
      do 18 i=1,nmms2
      t1 = fnum*t1/fden
      if(t1 .gt. sc20) then
      t1 = t1/sc40
      nex = nex+40
      end if
      fnum = fnum+2.
      fden = fden+2.
   18 continue
   20 t1 = t1/2.**(n-1-nex)
      if(mod(ma/2,2) .ne. 0) t1 = -t1
      t2 = 1. 
      if(ma .eq. 0) go to 26
      do 25 i=1,ma
      t2 = fnmh*t2/(fnmh+pm1)
      fnmh = fnmh+2.
   25 continue
   26 cp2 = t1*sqrt((n+.5)*t2)
      fnnp1 = n*(n+1)
      fnmsq = fnnp1-2.*ma*ma
      l = (n+1)/2
      if(mod(n,2) .eq. 0 .and. mod(ma,2) .eq. 0) l = l+1
      cp(l) = cp2
      if(m .ge. 0) go to 29
      if(mod(ma,2) .ne. 0) cp(l) = -cp(l)
   29 if(l .le. 1) return
      fk = n
      a1 = (fk-2.)*(fk-1.)-fnnp1
      b1 = 2.*(fk*fk-fnmsq)
      cp(l-1) = b1*cp(l)/a1
   30 l = l-1
      if(l .le. 1) return
      fk = fk-2.
      a1 = (fk-2.)*(fk-1.)-fnnp1
      b1 = -2.*(fk*fk-fnmsq)
      c1 = (fk+1.)*(fk+2.)-fnnp1
      cp(l-1) = -(b1*cp(l)+c1*cp(l+1))/a1
      go to 30
      end
