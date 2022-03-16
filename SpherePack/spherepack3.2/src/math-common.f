C> Cartesian to Spherical coordinate conversion
      subroutine ctos(x,y,z,r,theta,phi)
C> Inputs: x, y, and z coordinates
      real x, y, z
Cf2py intent(in) x, y, z
C> Outputs: r, theta and phi coordinates
      real r
C> Theta is angle down from z axis (colatitude)
      real theta
C> Phi is angle around z axis (longitude)
      real phi
Cf2py intent(out) r, theta, phi
      real r1
      r1 = x*x+y*y
c$$$      if (r1.ne.0.) go to 10
      if (r1 .eq. 0.) then
         phi = 0.
         theta = 0.
         if(z .lt. 0.) theta = 4.*atan(1.)
         r = abs(z)
         return
      else
 10      r = sqrt(r1+z*z)
         r1 = sqrt(r1) 
         phi = atan2(y,x)
         theta = atan2(r1,z)
         return
      endif
      end

C> Spherical to Cartesian coordinate conversion
      subroutine stoc(r,theta,phi,x,y,z)
C> Input, spherical radius
      real r
C> Input, angle down from z axis (colatitude)
      real theta
C> Input, angle around z axis (longitude)
      real phi
Cf2py intent(in) r, theta, phi
C> Output, x, y, and z coordinates
      real x, y, z
Cf2py intent(out) x, y, z
      real st
      st = sin(theta)
      x = r*st*cos(phi)
      y = r*st*sin(phi)
      z = r*cos(theta)
      return
      end

C> atan2 with arguments reversed
C> Returns 0 if x = y = 0
      function atanxy(x,y)
      atanxy = 0.
      if(x.eq.0. .and. y.eq.0.) return
      atanxy = atan2(y,x)
      return
      end

      subroutine coe(moe,n,x,dmax)
      double precision x(n),dmax
      integer moe, n
Cf2py intent(in) moe, n
Cf2py intent(out) dmax
Cf2py intent(in,out) x
      nh = (n+1)/2
      dmax = 0.
c$$$      if(moe.ne.0) go to 1
      if (moe .eq. 0) then
         do i=1,nh
            dmax = max(dmax,dabs(x(i)-x(n-i+1)))
            x(i) = .5*(x(i)+x(n-i+1))
            x(n-i+1) = x(i)
         end do
         return
      else
 1       do i=1,nh
            dmax = max(dmax,dabs(x(i)+x(n-i+1)))
            x(i) = .5*(x(i)-x(n-i+1))
            x(n-i+1) = -x(i)
         end do
         if(mod(n,2).ne.0) x(nh) = 0.
         return
      endif
      end

C> Double-precision version of mxmx
      subroutine dmxmx(lr,lc,ld,a,mc,md,b,x,y)
      double precision a(ld,*),b(md,*),x(ld,2),y(ld,2),
     1                 sum1,sum2
Cf2py intent(in) lr, lc, ld, a, mc, md, b, x
Cf2py intent(out) y
      do k=1,lr
         y(k,1) = 0.
         y(k,2) = 0.
      end do
c
      if(lc.le.0) return
      do i=1,lc
         sum1 = 0.
         sum2 = 0.
         do j=1,mc
            sum1 = sum1 + b(i,j)*x(j,1)
            sum2 = sum2 + b(i,j)*x(j,2)
         end do
         do k=1,lr
            y(k,1) = y(k,1)+sum1*a(k,i)
            y(k,2) = y(k,2)+sum2*a(k,i)
         end do
      end do
      return
      end

      subroutine mxmx(lr,lc,ld,a,mc,md,b,x,y)
      dimension a(ld,*),b(md,*),x(ld,2),y(ld,2)
Cf2py intent(in) lr, lc, ld, a, mc, md, b, x
Cf2py intent(out) y
      do k=1,lr
         y(k,1) = 0.
         y(k,2) = 0.
      end do
c
      if(lc.le.0) return
      do i=1,lc
         sum1 = 0.
         sum2 = 0.
         do j=1,mc
            sum1 = sum1 + b(i,j)*x(j,1)
            sum2 = sum2 + b(i,j)*x(j,2)
         end do
         do k=1,lr
            y(k,1) = y(k,1)+sum1*a(k,i)
            y(k,2) = y(k,2)+sum2*a(k,i)
         end do
      end do
      return
      end

      subroutine mxm(lr,lc,ld,a,mc,md,b,nd,c)
      double precision a(ld,*),b(md,*),c(nd,*)
Cf2py intent(in) lr, lc, ld, a, mc, md, b, nd
Cf2py intent(in,out) c
      do i=1,lr
         do j=1,mc
            c(i,j) = 0.
            do k=1,lc 
               c(i,j) = c(i,j)+a(i,k)*b(k,j)
            end do
         end do
      end do
      return
      end

C> Single-precision version of mxm
      subroutine smxm(lr,lc,ld,a,mc,md,b,nd,c)
      dimension a(ld,*),b(md,*),c(nd,*)
Cf2py 
      do i=1,lr
         do j=1,mc
            c(i,j) = 0.
            do k=1,lc 
               c(i,j) = c(i,j)+a(i,k)*b(k,j)
            end do
         end do
      end do
      return
      end

      subroutine trunc(irc,n,idp,a,nrc,ijs)
      double precision a,eps
      parameter (eps=5.d-8)
      dimension a(idp,*),ijs(n)
Cf2py intent(in) irc, n, idp, a, nrc
Cf2py intent(out) ijs
c
c     irc = 0 for columns , or irc = 1 for rows
c
c$$$  if(irc.ne.0) go to 30
      if (irc .eq. 0) then
         do 20 j=1,nrc
            do i=1,n
               ijs(j) = i
c$$$               if(dabs(a(i,j)) .gt. eps) go to 20
               if (dabs(a(i, j)) .gt. eps) then
                  exit
               end if
            end do
 20      continue
         return
      else
 30      do 50 i=1,nrc
            do j=1,n
               ijs(i) = j
c$$$               if(abs(a(i,j)) .gt. eps) go to 50
               if (abs(a(i, j)) .gt. eps) then
                  exit
               end if
            end do
 50      continue
         return
      end if
      end

C> Normalize x
      subroutine normal(n,x,id,q)
      dimension x(n),q(n)
      double precision x,q,sqs
Cf2py intent(in) n, id, q
Cf2py intent(in,out) x
c
c     normalize x
c
      sqs = 0.
      do i=1,n
c$$$      The below lines were commented out in one version of normal
c$$$      but not the other: I took what looked like the later (and more
c$$$      correct) version.
c      sum = 0.
c      do j=1,n
c      sum = sum+q(i,j)*x(j)
c      end do
c      sqs = sqs+sum*x(i)
         sqs = sqs+q(i)*x(i)*x(i)
      end do
c
c$$$  if(sqs .ne. 0) go to 4
      if (sqs .eq. 0.) then
         write(*,3)
 3       format(' norm of z is zero in subroutine normal')
         return
      else
 4       sqs = dsqrt(sqs)
         do i=1,n
            x(i) = x(i)/sqs
         end do
         return
      end if
      end

C> Accumulate inner products of x with respect to y
      subroutine gs(n,x,y,z)
      dimension x(n),y(n),z(n)
      double precision x,y,z,sum
c
c     accumulate innerproducts of x with respect to y.
c
Cf2py intent(in) n, x, y
Cf2py intent(in,out) z
      sum = 0.
      do i=1,n
         sum = sum+x(i)*y(i)
      end do
      do i=1,n
         z(i) = z(i)+sum*y(i)
      end do
      return
      end
