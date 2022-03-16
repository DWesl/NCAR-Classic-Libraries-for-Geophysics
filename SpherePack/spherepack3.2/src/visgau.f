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
c ... file visgau.f
c
c     contains documentation and code for subroutine visgau
c
      SUBROUTINE VISGAU (NLAT,NLON,H,LEN,EYER,EYELAT,EYELON,
     1                   THETA,WK,LWK,IWK,LIWK,IERROR)
c
c     subroutine visgau produces a three dimensional visible rendering
c     of the function h(i,j) which is tabulated on a gauss distributed
c     colatitude grid.
c
c     requires  setgau alfk lfpt gaqd drst dintql dpytha 
c               visgau embed intrpg sptcg diag stride
c               trigau vsurf vsurf1 prjct box icvmg projct 
c
c     tvisgau uses the ncar graphics package.
c     compile with: ncargf77 (all programs above)
c
c     execute with:  a.out
c
c     on screen display with:  ctrans -d x11 gmeta
c                          
c     print with:  ctrans -d ps.color gmeta > gmeta.ps
c                  lpr -p(your printer) gmeta.ps 
c
c     input parameters
c
c     nlat   the number of gauss colatitudes.
c            if nlat is odd the equator is located at
c            grid point i=(nlat+1)/2. if nlat is even the equator is
c            located half way between points i=nlat/2 and i=nlat/2+1.
c            nlat must be at least 3. note: on the half sphere, the
c            number of grid points in the colatitudinal direction is
c            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than or equal to 4. the efficiency of the computation is
c            improved when nlon is a product of small prime numbers.
c
c
c     h      a two dimensional array that contains the discrete
c            function to be displayed. h(i,j) is the distance from the
c            center of the sphere to the surface at gauss colatitude 
c            point theta(i) and longitude point 
c            phi(j) = (j-1)*2*pi/nlon.
c
c     len    the first dimension of the array h as it appears in the
c            program that calls sphere.
c
c     eyer   the distance from the center of the sphere to the eye.
c
c     eyelat the colatitudinal coordinate of the eye (in degrees).
c
c     eyelon the longitudinal  coordinate of the eye (in degrees).
c
c     theta  a double precision  array with nlat gauss colatitudes
c            computed by subroutine gaqd
c
c     wk     a real work array 
c
c     lwk    the dimension of the array wk as it appears in the
c            program that calls visgau. lwk must be at least 
c                       46*(nlat+2)*(nlon+1).
c
c     iwk    an integer work array
c
c     liwk   the dimension of the array iwk as it appears in the
c            program that calls visgau. liwk must be at least 
c                       14*(nlat+2)*(nlon+1).
c
c     ierror = 0    no error
c            = 1    the eye is positioned inside the sphere
c            = 2    lwk  is less than 46*(nlat+2)*(nlon+1)
c            = 3    liwk is less than 14*(nlat+2)*(nlon+1)
c
c  
      dimension h(len,nlon),wk(*)
      INTEGER IWK(*)
      double precision theta(nlat)
      n = nlat+2
      m = nlon+1
      mn = m*n
      ierror = 2
      if(lwk .lt. 46*mn) return
      ierror = 3
      if(liwk .lt. 14*mn) return
      ierror = 1
      do 10 j=1,nlon 
      do 10 i=1,nlat
      if(eyer .le. h(i,j)) return
   10 continue
      ierror = 0
c     ****     set up pointers to sub work arrays in wk
      ntri = mn+mn
      nw1 = 1
      nw2 = nw1+mn
      nclat = 1
      nslat = nclat+n
      nxp = 1
      nyp = nxp+mn
      nx1 = 1
      ny1 = nx1+ntri
      nz1 = ny1+ntri
      nx2 = nz1+ntri
      ny2 = nx2+ntri
      nz2 = ny2+ntri
      nx3 = nz2+ntri
      ny3 = nx3+ntri
      nz3 = ny3+ntri
      nx  = nz3+ntri
      ny = nx+mn
      nz = ny+mn
      nwrk = nx
      nitype = 1
      niflag = ntri+1
      nmst = niflag+mn
      nmfac = nmst+n
c     **** embed h in a larger array
      call embed(nlat,nlon,h,len,wk(nz1))
c     ****     mid-cell interpolation
      call intrpg(wk(nz1),m,n,wk(nw1),wk(nw2),iwk(niflag))
c     ****     transform grid points to cartesian coordinates
      call sptcg(wk(nz1),m,n,theta,wk(nclat),wk(nslat),wk(nx),wk(ny),
     1                                                      wk(nz))
c     ****     transform eye position to cartesian coordinates
      pi = 4.*atan(1.)
      dtr = pi/180.
      xeye=eyer*sin(dtr*eyelat)
      yeye=xeye*sin(dtr*eyelon)
      xeye=xeye*cos(dtr*eyelon)
      zeye=eyer*cos(dtr*eyelat)
c     ****     project grid points
      call projct(m,n,xeye,yeye,zeye,wk(nx),wk(ny),wk(nz),wk(nxp),
     1            wk(nyp)) 
c     ****     check for visibility of cell boundaries
      call diag(m,n,wk(nxp),wk(nyp),iwk(niflag))
c     ****     compute longitude stride as a function of latitude
      call stride(m,n,iwk(nmst),iwk(nmfac))
c     ****     perform triangulation
      call trigau(m,n,wk(nx),wk(ny),wk(nz),itri,wk(nx1),wk(ny1),
     1wk(nz1),wk(nx2),wk(ny2),wk(nz2),wk(nx3),wk(ny3),wk(nz3),
     2iwk(nitype),iwk(niflag),iwk(nmst))
c     ****     call surface plotting routine
      call vsurf(xeye,yeye,zeye,itri,wk(nx1),wk(ny1),wk(nz1),wk(nx2),
     1wk(ny2),wk(nz2),wk(nx3),wk(ny3),wk(nz3),iwk(nitype),wk(nwrk),
     1iwk(niflag))
      return
      end
      subroutine embed(nlat,nlon,h,len,hg)
      dimension h(len,nlon),hg(nlat+2,nlon+1)
      do 10 i=1,nlat
      do 10 j=1,nlon
      hg(i+1,j) = h(i,j)
   10 continue  
      sumn = 0.
      sums = 0.
      do 15 j=1,nlon
      sumn = sumn+h(1,j)
      sums = sums+h(nlat,j)
   15 continue  
      sumn = sumn/nlon
      sums = sums/nlon
      do 20 j=1,nlon
      hg(1,j) = sumn
      hg(nlat+2,j) = sums
   20 continue  
      do 25 i=1,nlat+2
      hg(i,nlon+1) = hg(i,1)
   25 continue  
      return
      end
      subroutine intrpg(h,m,n,w1,w2,iflag)
c     ****     interpolates to mid points of grid cells using second
c     ****     order formula
      dimension h(n,m),w1(n,m),w2(n,m+2),iflag(n,m),sten(4,4)
      data sten/.015625,2*-.078125,.015625,-.078125,2*.390625,
     12*-.078125,2*.390625,-.078125,.015625,2*-.078125,.015625/
c     ****     copy h to w2
      mm1 = m-1
      do 1 i=1,mm1
      do 1 j=1,n
      w2(j,i+1)=h(j,i)
    1 continue
c     ****     add periodic points
      do 2 j=1,n
      w2(j,1)=w2(j,m)
      w2(j,m+1)=w2(j,2)
      w2(j,m+2)=w2(j,3)
    2 continue
c     ****     perform interpolation
c     ****     set w1 to zero
      do 7 i=1,m
      do 7 j=1,n
      w1(j,i)=0.
    7 continue
c     ****     interpolate
      do 8 k=1,4
      do 8 l=1,4
      do 8 i=1,m-1
      do 8 j=2,n-2
      w1(j,i)=w1(j,i)+w2(j+l-2,i+k-1)*sten(k,l)
    8 continue
c     ****     set up iflag array
c     ****     iflag(j,i)=0  if diagonal is (j,i) to (j+1,i+1)
c     ****     iflag(j,i)=16 if diagonal is (j+1,i), (j,i+1)
      do 9 i=1,m-1
      do 9 j=2,n-2
      iflag(j,i)=icvmg(16,0,abs(.5*(w2(j,i+1)+w2(j+1,i+2))-w1(j,i))-
     1abs(.5*(w2(j,i+2)+w2(j+1,i+1))-w1(j,i)))
    9 continue
      return
      end
      subroutine sptcg(r,m,n,theta,clat,slat,x,y,z)
c     ****     transforms from spherical to cartesian coordinates
      dimension r(n,m),clat(n),slat(n),x(n,m),y(n,m),z(n,m)
      double precision theta(*)
      pi = 4.*atan(1.)
      dp = (pi+pi)/(m-1)
      clat(1) = 1.
      slat(1) = 0.
      do 10 j=2,n-1
      thet = theta(j-1)
      clat(j) = cos(thet)
      slat(j) = sin(thet)
   10 continue
      clat(n) = -1.
      slat(n) = 0.
      do 20 i=1,m-1
      clon = cos((i-1)*dp)
      slon = sin((i-1)*dp)
      do 20 j=1,n
      x(j,i)=r(j,i)*slat(j)
      y(j,i)=x(j,i)*slon
      x(j,i)=x(j,i)*clon
      z(j,i)=r(j,i)*clat(j)
   20 continue
      do 30 j=1,n
      x(j,m)=x(j,1)
      y(j,m)=y(j,1)
      z(j,m)=z(j,1)
   30 continue
      return
      end
      subroutine diag(m,n,xp,yp,iflag)
c
c     ****     label visibility of cell sides
c
c     north side corresponds to j
c     south side corresponds to j+1
c     west  side corresponds to i
c     east  side corresponds to i+1
c
c     let iflag = b4 b3 b2 b1 b0 (in binary) then b0 through b3 are 
c     either o or 1 depeending on whether the east, south, north
c     or west side is either invisible or visible, respectively.
c
c     b4 is o if the diagonal is from (i,j) to (i+1,j+1) and 1 if
c     the diagonal is from (i,j+1) to (i+1,j).
c
      dimension xp(n,m),yp(n,m),iflag(n,m)
c     ****     arithmetic statement function
      cp(j1,i1,j2,i2,j3,i3)=((xp(j1,i1)-xp(j2,i2))*(yp(j3,i3)-yp(j2,i2))
     1-(xp(j3,i3)-xp(j2,i2))*(yp(j1,i1)-yp(j2,i2)))
      do 100 j=2,n-2
      do 100 i=1,m-1
      if(iflag(j,i) .ge. 16) go to 20
      if(cp(j+1,i+1,j+1,i,j,i) .le. 0) go to 10 
c     west and south are visible
      iflag(j,i) = iflag(j,i)+10  
   10 if(cp(j,i,j,i+1,j+1,i+1) .le. 0) go to 100 
c     east and north are visible
      iflag(j,i) = iflag(j,i)+5  
      go to 100
   20 if(cp(j+1,i,j,i,j,i+1) .le. 0) go to 30 
c     west and north are visible
      iflag(j,i) = iflag(j,i)+12  
   30 if(cp(j,i+1,j+1,i+1,j+1,i) .le. 0) go to 100 
c     east and south are visible
      iflag(j,i) = iflag(j,i)+3  
  100 continue
c 
c     classify the poles
c
      do 200 i=1,m-1
      iflag(1,i) = 0
      if(cp(2,i+1,2,i,1,i) .gt. 0) iflag(1,i) = 15
      iflag(n-1,i) = 0
      if(cp(n,i,n-1,i,n-1,i+1) .gt. 0) iflag(n-1,i) = 31
  200 continue
      do 250 j=1,n-1
      iflag(j,m) = iflag(j,1)
  250 continue 
      return
      end
      subroutine stride(m,n,mst,mfac)
      dimension mfac(*),mtryh(3),mst(n),icl(8)
      data mtryh(1),mtryh(2),mtryh(3)/2,3,5/                        
      data icl(1),icl(2),icl(3),icl(4),icl(5),icl(6),icl(7),icl(8)
     1     /0,1,2,12,3,13,23,123/
c
c     find prime factors of m-1
c
      ml = m-1
      nf = 0           
      j = 0            
  101 j = j+1          
      if (j-3) 102,102,103                
  102 mtry = mtryh(j)                     
      go to 104        
  103 mtry = mtry+2    
  104 mq = ml/mtry     
      mr = ml-mtry*mq                     
      if (mr) 101,105,101                 
  105 nf = nf+1        
      mfac(nf) = mtry                    
      ml = mq          
      if (ml .ne. 1) go to 104            
      if(mfac(nf) .gt. 2) go to 106
      nf = nf-1
      mfac(nf) = 4
  106 tphi = .707/float(m-1)
      ns2 = n/2
      mf1 = mfac(nf)
      mst(1) = (m-1)/mf1
      pi = 4.*atan(1.)
      dt = pi/float(n-1)
      jf = nf-1
      do 110 jdo=2,ns2 
      j = jdo
      theta = (j-1)*dt
      st = sin(theta)
      mf2 = mf1*mfac(jf)
      if(abs(st/mf1-tphi) .gt. abs(st/mf2-tphi)) go to 115
      mst(j) = mst(j-1)
      go to 110
  115 mst(j) = (m-1)/mf2
      mf1 = mf2
      jf = jf-1
      if(jf .eq. 0) go to 120
  110 continue
  120 do 125 jdo=j,ns2
      mst(jdo) = 1
  125 continue
      do 130 jdo=1,ns2
      mst(n-jdo) = mst(jdo)
  130 continue
c      write (6,135) (mst(j),j=1,n)
  135 format(' colatitude strides'/(15i5))
c
      return
      end
      subroutine trigau(m,n,x,y,z,itri,x1,y1,z1,x2,y2,z2,x3,y3,z3,
     1                  ityp,iflag,mst)
c     ****     performs triangulation
      dimension x(n,m),y(n,m),z(n,m),x1(1),y1(1),z1(1),
     1x2(1),y2(1),z2(1),x3(1),y3(1),z3(1),ityp(1),iflag(n,m),
     2mst(n),icl(8) 
      data icl(1),icl(2),icl(3),icl(4),icl(5),icl(6),icl(7),icl(8)
     1     /0,1,2,12,3,13,23,123/
      itri = 0
      n1=2
      n2=n-2
      do 100 j=n1,n2
      do 100 i=1,m-1
      if(iflag(j,i) .ge. 16) go to 50
      if(mod(iflag(j,i),16) .lt. 8) go to 70
      itri = itri+1
      x1(itri) = x(j,i)
      y1(itri) = y(j,i)
      z1(itri) = z(j,i)
      x2(itri) = x(j+1,i)
      y2(itri) = y(j+1,i)
      z2(itri) = z(j+1,i)
      x3(itri) = x(j+1,i+1)
      y3(itri) = y(j+1,i+1)
      z3(itri) = z(j+1,i+1)
      ityph = 3
      if(mod(i-1,mst(j)) .eq. 0) go to 60
      if(mod(iflag(j,i-1),2) .eq. 0) go to 60
      ityph = ityph-1
   60 if(mod(iflag(j,i),2) .eq. 0) ityph = ityph+4
      ityp(itri) = icl(ityph+1)
   70 if(mod(iflag(j,i),2) .eq. 0) go to 100
      itri = itri+1
      x1(itri) = x(j,i)
      y1(itri) = y(j,i)
      z1(itri) = z(j,i)
      x2(itri) = x(j+1,i+1)
      y2(itri) = y(j+1,i+1)
      z2(itri) = z(j+1,i+1)
      x3(itri) = x(j,i+1)
      y3(itri) = y(j,i+1)
      z3(itri) = z(j,i+1)
      ityph = 0
      if(mod(iflag(j,i),16) .lt. 8) ityph = ityph+1
      if(mod(iflag(j,i+1),16) .lt. 8) ityph = ityph+2
      if(mod(iflag(j-1,i),4) .lt. 2) ityph = ityph+4
      ityp(itri) = icl(ityph+1)
      go to 100
   50 if(mod(iflag(j,i),16) .lt. 8) go to 20
      itri = itri+1
      x1(itri) = x(j,i)
      y1(itri) = y(j,i)
      z1(itri) = z(j,i)
      x2(itri) = x(j+1,i)
      y2(itri) = y(j+1,i)
      z2(itri) = z(j+1,i)
      x3(itri) = x(j,i+1)
      y3(itri) = y(j,i+1)
      z3(itri) = z(j,i+1)
      ityph = 1
      if(mod(i-1,mst(j)) .eq. 0) go to 10
      if(mod(iflag(j,i-1),2) .eq. 0) go to 10
      ityph = 0
   10 if(mod(iflag(j,i),2) .eq. 0) ityph = ityph+2
      if(mod(iflag(j-1,i),4) .lt. 2) ityph = ityph+4
      ityp(itri) = icl(ityph+1)
   20 if(mod(iflag(j,i),2) .eq. 0) go to 100
      itri = itri+1
      x1(itri) = x(j+1,i)
      y1(itri) = y(j+1,i)
      z1(itri) = z(j+1,i)
      x2(itri) = x(j+1,i+1)
      y2(itri) = y(j+1,i+1)
      z2(itri) = z(j+1,i+1)
      x3(itri) = x(j,i+1)
      y3(itri) = y(j,i+1)
      z3(itri) = z(j,i+1)
      ityph = 1
      if(mod(iflag(j,i+1),16) .lt. 8) ityph = ityph+2
      if(mod(iflag(j,i),16) .lt. 8) ityph = ityph+4
      ityp(itri) = icl(ityph+1)
  100 continue
c
c     ****     triangles around north and south poles
c
      do 200 i=1,m-1
      if(mod(iflag(1,i),16) .lt. 8) go to 250
      itri = itri+1
      x1(itri) = x(1,i)
      y1(itri) = y(1,i)
      z1(itri) = z(1,i)
      x2(itri) = x(2,i)
      y2(itri) = y(2,i)
      z2(itri) = z(2,i)
      x3(itri) = x(2,i+1)
      y3(itri) = y(2,i+1)
      z3(itri) = z(2,i+1)
      ityp(itri) = icl(3)
  250 if(mod(iflag(n-1,i),16) .lt. 8) go to 200 
      itri = itri+1
      x1(itri)=x(n-1,i)
      y1(itri)=y(n-1,i)
      z1(itri)=z(n-1,i)
      x2(itri)=x(n,i)
      y2(itri)=y(n,i)
      z2(itri)=z(n,i)
      x3(itri)=x(n-1,i+1)
      y3(itri)=y(n-1,i+1)
      z3(itri)=z(n-1,i+1)
      ityp(itri) = icl(1)
  200 continue
      return
      end
      call vsurf1(xeye,yeye,zeye,ntri,x1,y1,z1,x2,y2,z2,x3,y3,z3,
     1 itype,work,work(ntri+1),work(2*ntri+1),work(3*ntri+1),
     2 work(4*ntri+1),work(5*ntri+1),work(6*ntri+1),work(7*ntri+1),
     3 work(8*ntri+1),work(9*ntri+1),work(10*ntri+1),work(11*ntri+1),
     3 work(12*ntri+1),work(13*ntri+1),IWORK,IWORK(NTRI+1),
     4 IWORK(2*NTRI+1),IWORK(4*NTRI+1))
      return
      end
