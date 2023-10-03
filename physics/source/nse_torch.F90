module nse_torch
   use array_sizes
   use nuc_data
   use utils, only: r8
   use constants
   use errors
   implicit none
   ! FIXME: alignment broken because of numbered do loops
   private
   public testnse, indexx

   integer :: icount_nse = 0

contains

      subroutine testnse(t9,rho,ye,nvar,yps,ierr)

      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

!..declare
      integer :: i,indx(nsp),newguess,iprint,ipartition,iniab,imion,icoulomb,ierr
#if pIDX_RCLB == 3
      integer :: niso1(i325dim+1,iCfdim+2,3)
#else
      integer :: niso1(i282dim+1,iAtdim+2,3)
#endif
      real(r8) :: xtemp, xden, xye

! *** to merge in ppn....
      integer j,nvar
      real(r8) :: yps(1,nsp),t9,rho,ye,zn,an
      common/cnetw/an(nsp),zn(nsp)
! ^_^ utilities for trying to speed up
      integer cnt


!..formats
03    format(1x,a5,'=',1pe10.3,'  ',a5,'=',1pe10.3,'  ', a5,'=',1pe10.3,'  ',a5,'=',1pe10.3)

! *** rename variables called from ppn

        xtemp = t9 * 1.e9_r8
        xden  = rho
        xye   = ye

        ! unphysical Ye values can occur in the runge kutta integration, so just return an error
        ! so it can reduce its step size
        ierr = 0
        if ( ye < ZERO .or. ye > ONE ) then
           ierr = BAD_YE
           return
        end if


! *** check consistency of nsp versus abignet
        ! ^_^ abignet is defined to be == nsp, so we should not need to
        ! check this ... ?
        if (nsp /= abignet) then
                print*,'nsp must be equal to abignet:'
                print*,'check abignet in network.dek'
                stop
        end if 

! *** define parameters of nse calculations

! 1 *** ipartition = 1 partition functions not temperature dependent
!   *** ipartition = 2 partition functions temperature dependent
      ipartition = 2

! 2 *** iniab = 1 initial network as standard public_torch
!   *** iniab = 2 initial network according to ppn
      iniab = 2

! 3 *** imion = 1 mion = aion * amu
!   *** imion = 2 mion = (zion*m_p + (aion - zion)*m_n + exion)* amu
!   *** exion is the mass excess coming from REACLIB-FRDM formula
      imion = 1

! 4 *** icoulomb = 1 screening effect not included
!   *** icoulomb = 2 screening effect includedaccording to Calder et al. 2007
      icoulomb = 2

!..initialize the network
      call init_network(t9,xden,xye,niso1,ipartition,iniab,imion,icoulomb)

!..print some statistics of the root find
!..and always use new initial guesses

      iprint   = 0
      newguess = 1
      icount_nse=icount_nse+newguess
      if (icount_nse.gt.1)newguess=0

      call nse(xtemp,xden,xye,newguess,niso1,xmass_nse,iprint)
! *** here I put as a minimum value of xmass_nse 1.d-99,
! *** to avoid crazy low values like 1.d-255.
       xmass_nse(1:ionmax)=max(xmass_nse(1:ionmax),1.0d-99)

      if ( iprint  ==  1 ) then
         write(*,*) ' '
         write(*,*) 'top 20 abundances:'
      end if

      call indexx(ionmax,xmass_nse,indx)

      if ( iprint  ==  1 ) then
         write(*,03) (ionam(indx(i)),xmass_nse(indx(i)), i=ionmax,ionmax-19,-1)
         write(*,*) ' '
      end if

!***   set abundances according to nse_calculations
       yps(:,:) = 0._r8

       do i = 1, nsp
          if (considerisotope(i))then
             do j = 1, ionmax
                if (zis(i) == ionam(j)) then
                   yps(1,i) = xmass_nse(j)
                   exit
                end if
             end do
          end if
       end do

      end subroutine testnse



      subroutine nse(tt,dd,yye,newguess,niso1,xmass_out,iprint)

      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

!..given the temperature tt, density dd, and electron mole number ye
!..this routine puts a chosen reaction network into its nse distribution.

!..input:
!..tt = temperature
!..dd = density
!..ye = electron mol number
!..newguess = 1 = a new initial guess is made
!..         = 0 = values from the previous call are used as the initial guess
!..iprint = print flag

!..output:
!..xmass_out  = nse mass fraction


!..declare the pass
      integer          newguess,iprint
      double precision tt,dd,yye,xmass_out(abignet)

!..communicate
      double precision temp,den,ye_want,beta
      common /nsec1/   temp,den,ye_want,beta


!..locals
!      external         nsefunc
      logical          check
      integer          ntrial,nfev,ntaken,n,i !,j
      parameter        (ntrial = 200, n = 2)
      double precision x(n),amass,fac1,fac2,tolf,tolx,twopi,dum,resid(n),xmassini56,zni56,ani56
      parameter        (tolf = 1.0d-8, tolx = 1.0d-14, twopi=2.0d0*pi)

#if pIDX_RCLB == 3
      integer          niso1(i325dim+1,iCfdim+2,3)
#else
      integer          niso1(i282dim+1,iAtdim+2,3)
#endif

!..fill the common block
      temp    = tt
      den     = dd
      ye_want = yye
      beta    = 1.0d0/(kerg * temp)


! *** partition functions are calculated in init_network subroutine
! *** together with binding energies, not anymore here!

! ***   here x1 and x2 are defined. Some assumption is made.
!..here is an initial guess for the neutron and proton chemical potentials,
!..(x1) and (x2) respectively. obtained by setting xmass(ini56) = 1,
!..setting mup = mun, and inverting the saha equation.
!..all nse networks should at least have ni56 and this appears to be a
!..decent guess for all temp, rho, ye combinations.

! *** old fashion public_torch
!       amass  = aion(i) * amu
!       fac1   = aion(i)/(avo * den) * wpart(i)

      if (newguess  ==  1) then
       xmassini56 = 1.d0
       zni56      = 28.d0
       ani56      = 56.d0
       newguess = 0
       i      = niso1(int(ani56)+1,int(zni56)+2,1)
       if (iprint == 1 )print*,niso1(int(ani56)+1,int(zni56)+2,1)
!       do j=1,ionmax
!        print*,j,bion(j),mion(j),aion(j),wpart(j)
!       end do
!       stop'sono qui'
!       bion(i) = 484.003

!       print*,'wpart ni56',wpart(i)

! *** this is like I think should be done
       amass  = mion(i)
       fac1   = amass/den * wpart(i)
! *** this is how Frank is doing in test_nse_eos.f
!       amass = aion(i) * amu
!       fac1 = amass/den

       fac2   = (twopi/(beta*h) * amass/h )**1.5d0
       !print *,amass, fac1, fac2
!       x(1)   = -(log(fac1*fac2)/beta + bion(i)*ev2erg*1.0d6)/aion(i)
       x(1)   = (log(xmassini56) - log(fac1*fac2)/beta - bion(i)*ev2erg*1.0d6)/aion(i)
       x(2)   = x(1)
      end if
      if (iprint == 1) then
         print *,'newguess=',newguess
         print *,'x(1) guess:',x(1),xmassini56,amass,wpart(i)
         print *,'fac=',fac1,fac2,bion(i),mion(i)
      end if

!..root find on mass and charge conservation for
!..the chemical potentials of protons and neutrons
      call xnewt_nse(ntrial,x,n,tolx,tolf,ntaken,check,nfev,nsefunc)

!..be sure we converged
      if (check .or. ntaken  ==  ntrial) then
       write(*,*)
       write(*,*) 'check convergence of root finder'
       write(*,*)
      end if

!..some optional diagnostics
      if (iprint  ==  1) then
       write(*,*)
       write(*,110) 'iterations taken             =',ntaken
       write(*,110) 'function evals               =',nfev
       write(*,111) 'roots                        =',x(1),x(2)
       call nsefunc(dum,x,resid)
       write(*,111) 'mass conservation   residual =',resid(1)
       write(*,111) 'charge conservation residual =',resid(2)

 110   format(1x,a,i4)
 111   format(1x,a,1p2e14.6)
      end if


!..fill the output array using the converged values
! *** in nsefunc sub. xmass is calculated.
      call nsefunc(dum,x,resid)
      do i=1,ionmax
         xmass_out(i) = xmass(i)
      enddo

      return

      end subroutine nse

      subroutine nsefunc(x,y,f)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

!..this routine returns the root find functions.
!..input is the point x and y a vector of the unknowns.
!..output is the vector of root find functions f, which should be the
!..zero vector upon convergence.

!..y(1) is input as the neutron chemical potential
!..y(2) is input as the proton chemical potential


!..declare the pass
      double precision x,y(*),f(*)


!..locals
      integer          i,indx(nsp),j
      double precision ye,mu,amass,fac1,fac2,fac3,sum,twopi,sum2
      parameter        (twopi = 2.0d0 * pi)


!..communicate
      double precision temp,den,ye_want,beta
      common /nsec1/   temp,den,ye_want,beta
      common /mu1/   mu_c_p,mu_c

      double precision mu_c(nsp),mu_c_p


!..chemical potential and mass fraction of each isotope
!..hartmann et al, 1985 apj 297 837, eq 2
!..take the mass of each isotope to be amu * aion, otherwise a mass formula
! ***   Maxwell-Boltzmann equation in NSE

! *** old fashion public_torch
!       amass  = aion(i) * amu
!       fac1   = aion(i)/(avo * den) * wpart(i)

      do i=1,ionmax
       mu       = (aion(i) - zion(i))*y(1) + zion(i)*y(2)
       amass  = mion(i)
       fac1   = amass/den * wpart(i)
       fac2     = ( twopi/(beta*h) * amass/h )**1.5d0
       fac3     = exp( beta * (mu + bion(i)*ev2erg*1.0d6) - mu_c(i) + zion(i)*mu_c_p)
       xmass(i) = fac1 * fac2 * fac3
      enddo


!..sum the mass fractions in ascending order to minimize roundoff
      call indexx(ionmax,xmass,indx)
      sum   = 0.0d0
      do i=1,ionmax
       sum   = sum + xmass(indx(i))
      enddo

!..sum the mass fractions to form ye
! *** here below I calculate sum2 like in research nse program of Frank.
! *** In principle, it should be applicable also if mion(i) = aion * amu,
! *** and used as a general formula. Read Ivo's article
! *** seitenzahl et al 2007. Then discuss with Frank.

! *** another question for Frank:
! *** why when ye is calculated in this point instead of j = indx(i)
! *** integer i is used ? Does it make any difference ? Why ?
      sum2 = 0.0d0
      ye = 0.0d0
      do i=1,ionmax
       j = indx(i)
!       ye = ye + zion(j)/aion(j) * xmass(j)
! *** to debug imion sum2!
       sum2 = sum2 + amu/mion(j) * ((ye_want - 1.0d0)*zion(j) + ye_want * (aion(j) - zion(j))) * xmass(j)
      enddo

!..mass and charge conservation are the requirements
      f(1) = sum - 1.0d0
!      f(2) = ye - ye_want
! *** to debug imion!
      f(2) = sum2
      return
      end



      subroutine nsejac(y,f,dfdy,np)

! *** I think that in this subroutine the value of f an dfdy for protons
! *** and neutrons are calculated.

      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

!..this routine returns the functions and the jacobian to do the root find on
!..input is x, and y(n) a vector of the unknowns. output is f(n)
!..and its jacobian dfdy(np,np).

!..y(1) is the neutron chemical potential
!..y(2) is the proton chemical potential


!..declare the pass
      integer          np
      double precision y(2),f(2),dfdy(np,np)


!..locals
      integer          indx(ionmax),i,j
      double precision mu,mubn,mubp,amass,fac1,fac2,fac3,fac4,fac5, &
                       xmbn(ionmax),xmbp(ionmax),sum,sumbn,sumbp,   &
                       ye,yebn,yebp,twopi,sum2,sum2bn,sum2bp
      parameter        (twopi = 2.0d0 * pi)


!..communicate
      double precision temp,den,ye_want,beta
      common /nsec1/   temp,den,ye_want,beta
      common /mu1/   mu_c_p,mu_c

      double precision mu_c(nsp),mu_c_p

!..chemical potential and mass fraction of each isotope
!..hartmann et al, 1985 apj 297 837, eq 2
!..take the mass of each isotope to be amu * aion, otherwise a mass formula
! *** old fashion public_torch
!       amass  = aion(i) * amu
!       fac1   = aion(i)/(avo * den) * wpart(i)
!        print *,'chem neut',y(1),'chem prot',y(2)

      do i=1,ionmax
       mu       = (aion(i) - zion(i)) * y(1) + zion(i) * y(2)
       mubn     = aion(i) - zion(i)
       mubp     = zion(i)

       amass  = mion(i)
       fac1   = amass/den * wpart(i)
       fac2     = ( twopi/(beta*h) * amass/h )**1.5d0
       fac3     = exp( beta * (mu + bion(i) * ev2erg * 1.0d6) - mu_c(i) + zion(i)*mu_c_p)
       fac4     = fac1 * fac2 * fac3

       xmass(i) = fac4
       xmbn(i)  = fac4 * beta * mubn
       xmbp(i)  = fac4 * beta * mubp

      enddo

!..sum the mass fractions in ascending order to minimize roundoff
      call indexx(ionmax,xmass,indx)
      sum   = 0.0d0
      sumbn = 0.0d0
      sumbp = 0.0d0
      do i=1,ionmax
       j     = indx(i)
       sum   = sum   + xmass(j)
       sumbn = sumbn + xmbn(j)
       sumbp = sumbp + xmbp(j)
      enddo


!..sum the mass fractions to form ye
! *** fac 5 calculated considering the possibility that real mass of the
! *** isotope is considered in the calculations.
! *** In principle, it should be applicable also if mion(i) = aion * amu,
! *** and used as a general formula. Read Ivo's article
! *** seitenzahl et al 2007. Then discuss with Frank.

      ye   = 0.0d0
      yebn = 0.0d0
      yebp = 0.0d0
      sum2   = 0.0d0
      sum2bn = 0.0d0
      sum2bp = 0.0d0
      do i=1,ionmax
! *** originally was used index instead of i, but does not change anything.
!       j    = indx(i)
!       fac5 = zion(j)/aion(j)
! *** to debug imion for fac5!
       fac5 = amu/mion(i) * ((ye_want - 1.0d0)*zion(i) + ye_want * (aion(i) - zion(i)))
!   print*,'fac5',fac5,zion(i),aion(i)
!       ye   = ye   + fac5 * xmass(j)
!       yebn = yebn + fac5 * xmbn(j)
!       yebp = yebp + fac5 * xmbp(j)
! *** to debug imion sum2!
       sum2   = sum2 + fac5 * xmass(i)
       sum2bn = sum2bn + fac5 * xmbn(i)
       sum2bp = sum2bp + fac5 * xmbp(i)
      enddo

!..mass and charge conservation are the requirements
      f(1) = sum - 1.0d0
!      f(2) = ye - ye_want
! *** to debug imion f(2)!
      f(2) = sum2


!..jacobian
      dfdy(1,1) = sumbn
      dfdy(1,2) = sumbp
!      dfdy(2,1) = yebn
!      dfdy(2,2) = yebp
! *** to debug imion dfdy!
      dfdy(2,1) = sum2bn
      dfdy(2,2) = sum2bp

      return
      end subroutine nsejac



      subroutine xnewt_nse(ntrial,x,n,tolx,tolf,ntaken,check,nfev,func)
      include 'implno.dek'

!..given an initial guess x(1:n) for the root of n equations, this routine
!..finds the root by a globally convergent newtons method. the vector of
!..functions to be zeroed, called fvec(1:n) in the routine below, is
!..returned by the user supplied routine func. the output quantity check
!..is false on nomal return and true if the routine has converged to a
!..local minimum of the function xfminx_nse. if so, try restarting from a
!..different initial guess.

!..np is the maximum number of equations n
!..ntrial is the maximum number of iterations to try
!..ntaken is the number of iterations done
!..tolf sets the convergence on function values
!..tolmin sets the criterion for deciding wheather spurious convergence to
!..       a false minimum of xfminx_nse has occured
!..tolx is the convergence criteria on deltax
!..stpmx is the scaled maximum step length allowed in the line searches
!..nfev is the number of function evaluations


!..declare the pass
      external         func
      logical          check
      integer          ntrial,n,ntaken,nfev
      double precision x(n),tolx,tolf


!..common block communicates values from routine xfminx_nse
      integer          nn,np
      parameter        (np = 4)
      double precision fvec(np)
      common /newtnse/ fvec,nn

!..locals
      integer          i,its,j,indx(np)
      double precision tolmin,stpmx,d,den,f,fold,stpmax,sum,temp,test, &
                       !fjac(np,np),g(np),p(np),xold(np),xfminx_nse
                       fjac(np,np),g(np),p(np),xold(np)

      parameter        (tolmin = 1.0d-12, &
                        stpmx = 2.0d0)

!..initialize
      if (n .gt. np) stop 'n > np in routine xnewt'
      nn     = n
      f      = xfminx_nse(x,func)
      nfev   = 1
      ntaken = 0


!.. test for the initial guess being a root, using a more stringent tolf
      test = 0.0d0
      do i=1,n
       if (abs(fvec(i)) .gt. test) test = abs(fvec(i))
      enddo

      if (test .lt. 0.01_r8*tolf) then
       check = .false.
       return
      end if


!..get stpmax for the line search
      sum = 0.0d0
      do i=1,n
       sum = sum + x(i)*x(i)
      enddo
      stpmax = stpmx * max(sqrt(sum),real(n,r8))


!..start of iteration loop; get the jacobian
      do its = 1, ntrial
       ntaken = its


!..second order accurate numerical jacobian
!       call jac_nse(dum,x,fjac,n,n,np,np,func)
!       nfev = nfev + 2*n + 1


!..analytic jacobian
       call nsejac(x,fvec,fjac,np)
       nfev = nfev + 1


!..compute grad f for the line searches
       do i=1,n
        sum = 0.0d0
        do j=1,n
         sum = sum + fjac(j,i)*fvec(j)
        enddo
        g(i) = sum
       enddo


!..store x, and f and form right hand sides
       do i=1,n
        xold(i) = x(i)
       enddo
       fold = f
       do i=1,n
        p(i) = -fvec(i)
       enddo

!..solve the linear systems
       call ludcmp_nse(fjac,n,np,indx,d)
       call lubksb_nse(fjac,n,np,indx,p)

!..line search returns new x and f
!..it also gets fvec at the new x when it calls xfminx_nse
!   print *,'xold',xold
!   print *,'x',x
       call lnsrch_nse(n,xold,fold,g,p,x,f,stpmax,check,nfev,func)
!   print *,'xold',xold
!   print *,'x',x

!..test for convergence on function value
       test = 0.0d0
       do i=1,n
        if (abs(fvec(i)) .gt. test) test = abs(fvec(i))
       enddo
       if (test .lt. tolf) then
        check = .false.
        return
       end if

!..check for zero gradiant, i.e. spurious convergence
       if (check) then
        test = 0.0d0
        den  = max(f, 0.5d0 * n)
        do i=1,n
         temp = abs(g(i)) * max(abs(x(i)),1.0d0)/den
         if (temp .gt. test) test = temp
        enddo
        if (test .lt. tolmin) then
         check = .true.
        else
         check = .false.
        end if
        return
       end if

!..test for convergence on deltax
       test = 0.0d0
       do i=1,n
        temp = (abs(x(i)-xold(i)))/max(abs(x(i)),1.0d0)
        if (temp .gt. test) test = temp
       enddo
       if (test .lt. tolx) return

!..back for another iteration
      enddo
      check = .true.
      return
      end subroutine xnewt_nse



      subroutine lnsrch_nse(n,xold,fold,g,p,x,f,stpmax,check,nfev,func)
      include 'implno.dek'

!..given an n dimensional point xold(1:n), the value of the function fold
!..and the gradient g(1:n) at the point, and a direction p(1:n), this routine
!..finds a new point x(1:n) along the direction of p from xold where the
!..function xfminx_nse has decreased "sufficiently". the new function value is
!..returned in f. stpmax is an input quanity that limits the length of the
!..steps so that the function is not evaluated in regions where it is
!..undefined or subject to overflow. p is usually the newton direction. the
!..output quantity check is false on normal exit, and true when x is too
!..close to xold. in a minimization routine, this usually signals
!..convergence and can be ignored. however, in a root finding routine, the
!..calling routine should check wheather the convergence is spurious.


!..declare the pass
      external         func
      logical          check
      integer          n,nfev
      double precision f,fold,stpmax,g(n),p(n),x(n),xold(n)


!..locals
      integer          i
      !double precision xfminx_nse,a,alam,alam2,alamin,b,disc,f2,rhs1, &
      double precision a,alam,alam2,alamin,b,disc,f2,rhs1, &
                       rhs2,slope,sum,temp,test,tmplam, &
                       alf,tolx
      parameter        (alf  = 1.0d-4, &
                        tolx = 3.0d-15)


!..alf ensures sufficient decrease in the function value, tolx is the
!..convergence criterion on deltax


!..initialize and scale if the attempted step is too big
      check = .false.
      sum   = 0.0d0
      do i=1,n
       sum = sum + p(i)*p(i)
      enddo
      sum = sqrt(sum)
      if (sum .gt. stpmax) then
       do i=1,n
        p(i) = p(i) * stpmax/sum
       enddo
      end if
      slope = 0.0d0
      do i=1,n
       slope = slope + g(i)*p(i)
      enddo
      if (slope .ge. 0.0_r8) stop 'roundoff problem in lnsrch_nse'


!..compute lambda_min
      test = 0.0d0
      do i=1,n
       temp = abs(p(i))/max(abs(xold(i)),1.0d0)
       if (temp .gt. test) test = temp
      enddo
      alamin = tolx/test


!..always try a full newton step, start of iteration loop
      alam = 1.0d0
1     continue
      do i=1,n
       x(i) = xold(i) + alam*p(i)
      enddo

      f    = xfminx_nse(x,func)
      nfev = nfev + 1



!..convergence on deltax, for root finding, the calling routine
!..should verify the convergence
      if (alam .lt. alamin) then
       do i=1,n
        x(i) = xold(i)
       enddo
       check = .true.
       return

!..sufficient function decrease
      else if (f .le. fold + alf*alam*slope) then
       return

!..backtrack
      else
       if (alam  ==  1.0_r8) then
        tmplam = -slope / (2.0d0 * (f-fold-slope))
       else
        rhs1 = f  - fold - alam*slope
        rhs2 = f2 - fold - alam2*slope
        a    = (rhs1/alam**2 - rhs2/alam2**2)/(alam-alam2)
        b    = (-alam2*rhs1/alam**2 + alam*rhs2/alam2**2) / (alam-alam2)
        if (a  ==  0.0_r8) then
         tmplam = -slope/(2.0d0 * b)
        else
         disc = b*b - 3.0d0 * a * slope
         if (disc .lt. 0.0_r8) then
          tmplam = 0.5d0 * alam
         else if (b .le. 0.0_r8) then
          tmplam = (-b + sqrt(disc)) / (3.0d0 * a)
         else
          tmplam = -slope/(b + sqrt(disc))
         end if
        end if
        if (tmplam .gt. 0.5d0*alam) tmplam = 0.5d0*alam
       end if
      end if

!..store for the next trip through
      alam2 = alam
      f2    = f
      alam  = max(tmplam, 0.1d0*alam)
      goto 1
      end






      double precision function xfminx_nse(x,func)
      include 'implno.dek'


!..returns f = 0.5 f dot f at x. func is a user supplied routine of the
!..functions to be root found.

!..declare the pass
      external         func
      double precision x(*)


!..locals
      integer          i
      double precision sum,dum


!..common block communicates values back to routine xnewt
      integer          nn,np
      parameter        (np = 4)
      double precision fvec(np)
      common /newtnse/ fvec,nn


!      print *,'dum',dum
!      print *,'x',x(1)
!      print *,'fvec',fvec
      call func(dum,x,fvec)
!      print *,'dum',dum
!      print *,'x',x(1)
!      print *,'fvec',fvec

      sum = 0.0d0
      do i=1,nn
       sum = sum + fvec(i)*fvec(i)
      enddo

      xfminx_nse = 0.5d0 * sum
      return
      end






      subroutine jac_nse(x,y,dfdy,mcol,nrow,mmax,nmax,derivs)
      include 'implno.dek'

!..this routine computes a second order accurate jacobian matrix
!..of the function contained in the routine derivs.
!..
!..input is the point x and the the vector y(nrow) at which to compute the
!..jacobian dfdy(mcol,nrow).
!..
!..uses 2*nrow + 1 function evaluations


!..declare the pass
      external         derivs
      integer          mcol,nrow,mmax,nmax
      double precision x,y(nmax),dfdy(mmax,nmax)


!..locals
      integer          i,j,imax
      parameter        (imax = 4)
      double precision fminus(imax),fplus(imax),rel,ax,temp,h,hinv
      parameter        (rel = 3.162278d-8, &
                        ax  = 1.0d-16)


!..check
       if (nrow .gt. imax) stop 'nrow > imax in jacobian2'


!..for each row, get the right stepsize
      do j=1,nrow
       temp = y(j)
       h    = rel * max(abs(y(j)),ax)
       y(j) = temp + h
       h    = y(j) - temp
       call derivs(x,y,fplus)
       y(j) = temp

       temp = y(j)
       y(j) = temp - h
       h    = temp - y(j)
       call derivs(x,y,fminus)
       y(j) = temp

!..compute the jth row of the jacobian
        hinv = 1.0d0/(2.0d0 * h)
        do i=1,mcol
         dfdy(i,j) = (fplus(i) - fminus(i)) * hinv
        enddo
       enddo

!..restore the original state
      call derivs(x,y,fplus)
      return
      end




!..lu decomposition:
!..routine ludcmp does a pivoting lower-upper decomposition
!..routine lubksb does the backsubstitution from ludcmp



      subroutine ludcmp_nse(a,n,np,indx,d)
      implicit none
      save

!..given the matrix a(n,n), with physical dimensions a(np,ap) this routine
!..replaces a by the lu decompostion of a row-wise permutation of itself.
!..input are a,n,np. output is a, indx which records the row
!..permutations effected by the partial pivoting, and d which is 1 if
!..the number of interchanges is even, -1 if odd.
!..use routine lubksb to solve a system of linear equations.
!..
!..nmax is the largest expected value of n

!..declare
      integer          n,np,indx(np),nmax,i,j,k,imax
      parameter        (nmax=500)
      double precision a(np,np),d,tiny,vv(nmax),aamax,sum,dum
      parameter        (tiny=1.0d-20)


!..vv stores the implicit scaling of each row
!..loop over the rows to get the scaling information
      d = 1.0d0
      do i=1,n
       aamax = 0.0d0
       do j=1,n
        if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
       enddo
       if (aamax  ==  0.0_r8) stop 'singular matrix in ludcmp'
       vv(i) = 1.0d0/aamax
      enddo

!..for each column apply crouts method; see equation 2.3.12
      do j=1,n
       do i=1,j-1
        sum = a(i,j)
        do k=1,i-1
         sum = sum - a(i,k)*a(k,j)
        enddo
        a(i,j) = sum
       enddo

!..find the largest pivot element
       aamax = 0.0d0
       do i=j,n
        sum=a(i,j)
        do k=1,j-1
         sum = sum - a(i,k)*a(k,j)
        enddo
        a(i,j) = sum
        dum = vv(i)*abs(sum)
        if (dum .ge. aamax) then
         imax  = i
         aamax = dum
        end if
       enddo

!..if we need to interchange rows
       if (j .ne. imax) then
        do k=1,n
         dum       = a(imax,k)
         a(imax,k) = a(j,k)
         a(j,k)    = dum
        enddo
        d          = -d
        vv(imax)   = vv(j)
       end if

!..divide by the pivot element
       indx(j) = imax
       if (a(j,j)  ==  0.0_r8) a(j,j) = tiny
       if (j .ne. n) then
        dum = 1.0d0/a(j,j)
        do i=j+1,n
         a(i,j) = a(i,j)*dum
        enddo
       end if

!..and go back for another column of crouts method
      enddo
      return
      end





      subroutine lubksb_nse(a,n,np,indx,b)
      implicit none
      save

!..solves a set of n linear equations ax=b. a is input in its lu decomposition
!..form, determined by the routine above ludcmp. indx is input as the
!..permutation vector also returned by ludcmp. b is input as the right hand
!..side vector and returns with the solution vector x.
!..a,n ans np are not modified by this routine and thus can be left in place
!..for successive calls (i.e matrix inversion)

!..declare
      integer           n,np,indx(np),i,ii,j,ll
      double precision  a(np,np),b(np),sum

!..when ii is > 0, ii becomes the index of the first nonzero element of b
!..this is forward substitution of equation 2.3.6, and unscamble in place
      ii = 0
      do i=1,n
       ll = indx(i)
       sum = b(ll)
       b(ll) = b(i)
       if (ii .ne. 0) then
        do j=ii,i-1
         sum = sum - a(i,j) * b(j)
        enddo

!..nonzero element was found, so dos the sums in the loop above
       else if (sum .ne. 0.0_r8) then
        ii  = i
       end if
       b(i) = sum
      enddo

!..back substitution equation 2.3.7
      do i = n,1,-1
       sum = b(i)
       if (i .lt. n) then
        do j=i+1,n
         sum = sum - a(i,j) * b(j)
        enddo
       end if
       b(i) = sum/a(i,i)
      enddo
      return
      end







      subroutine indexx(n,arr,indx)
      include 'implno.dek'
!..
!..indexes an array arr(1:n). that is it outputs the array indx(1:n) such
!..that arr(indx(j)) is in ascending order for j=1...n. the input quantities
!..are not changed.
!..
!..declare
      integer          n,indx(n),m,nstack
      parameter        (m=7, nstack = 50)
      integer          i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
      double precision arr(n),a
!..
!..initialize
      do 11 j=1,n
       indx(j) = j
11    continue
      jstack = 0
      l      = 1
      ir     = n
!..
!..insertion sort when subbarray small enough
1     if (ir - l .lt. m) then
       do 13 j=l+1,ir
        indxt = indx(j)
        a     = arr(indxt)
        do 12 i=j-1,l,-1
         if (arr(indx(i)) .le. a) go to 2
         indx(i+1) = indx(i)
12      continue
        i = l - 1
2       indx(i+1) = indxt
13     continue
!..
!..pop stack and begin a new round of partitioning
       if (jstack  ==  0) return
       ir     = istack(jstack)
       l      = istack(jstack-1)
       jstack = jstack - 2
!..
!..choose median of left, center and right elements as partitioning element
!..also rearrange so that a(l+1) < a(l) < a(ir)
      else
       k         = (l + ir)/2
       itemp     = indx(k)
       indx(k)   = indx(l+1)
       indx(l+1) = itemp

       if (arr(indx(l)) .gt. arr(indx(ir))) then
        itemp    = indx(l)
        indx(l)  = indx(ir)
        indx(ir) = itemp
       end if


       if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
       endif
       if(arr(indx(l)).gt.arr(indx(l+1)))then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
       endif

!..
!..initialize pointers for partitioning
       i     = l + 1
       j     = ir
       indxt = indx(l+1)
       a     = arr(indxt)
3      continue
       i = i + 1
       if (arr(indx(i)) .lt. a) go to 3
4      continue
       j = j - 1
       if (arr(indx(j)) .gt. a) go to 4
       if (j .lt. i) go to 5
       itemp   = indx(i)
       indx(i) = indx(j)
       indx(j) = itemp
       go to 3
!..
5      indx(l+1) = indx(j)
       indx(j)   = indxt
       jstack    = jstack + 2
!..
!..push pointers to larger subarray on stack
       if (jstack .gt. nstack) stop 'jstack > nstack in routine indexx'
       if (ir - i + 1  .ge.  j - l) then
        istack(jstack)   = ir
        istack(jstack-1) = i
        ir               = j - 1
       else
        istack(jstack)   = j-1
        istack(jstack-1) = l
        l                = i
       end if
      end if
      go to 1
      end




      subroutine init_network(t9,den,ye_want,niso1,ipartition,iniab,imion,icoulomb)
      use reaclib, only: anumpart, znumpart, istate, spinpart, partnum, &
                         exmass, partinter, &
                         reaclib_interpolate_partition_functions
      use nuc_data, only: isomeric_state
      include 'implno.dek'
      include 'network.dek'
      include 'const.dek'

!..this routine initializes the nse network

!..declare
      integer          i, k,j
! *** note: I put 3 for isomer state instead of 2 because reaclib is using
! *** notation different for therm/ground/isomer. ppn uses 2.
#if pIDX_RCLB == 3
      integer          niso1(i325dim+1,iCfdim+2,3)
#else
      integer          niso1(i282dim+1,iAtdim+2,3)
#endif
      double precision zn,an,t9
!..for easy zeroing of the isotope pointers
      integer          isotp(nisotp)
      equivalence      (isotp(1),ih1)   ! equivalence ???
      common/cnetw/an(nsp),zn(nsp)
      integer isomer_state(nsp)
      integer ipartition,iniab,imion,icoulomb
! *** fac1 = factor of conversion amu in MeV
      double precision m_p,m_n,fact1
      DATA fact1/931.494d0/
!      DATA fact1/931.4d0/
      double precision mu_c(nsp),mu_c_p,xne,ge,tolf, &
        tolx,twopi,esqu,forthpi,third,fivth,a1,a2,a2inv,a3, &
        rt3,half_rt3,gi,sqrtgi,den,beta,ye_want,deltap,deltan, &
        mev2erg,mev2gr
      common /mu1/   mu_c_p,mu_c
      parameter        (tolf     = 1.0d-10, &
                        tolx     = 1.0d-12, &
                        twopi    = 2.0d0*pi, &
                        esqu     = qe*qe, &
                        forthpi  = 4.0d0 * pi/3.0d0, &
                        third    = 1.0d0/3.0d0, &
                        fivth    = 5.0d0/3.0d0, &
                        a1       = -0.9052d0, &
                        a2       = 0.6322d0, &
                        a2inv    = 1.0d0/a2, &
                        rt3      = 1.7320508075688772d0, &
                        half_rt3 = 0.5d0 * rt3, &
                        mev2erg = ev2erg*1.0d6, &
                        mev2gr  = mev2erg/clight**2)

      m_p = mp/amu   ! proton mass in amu
      m_n = mn/amu   ! neutron mass in amu
! *** note : mp,mn,amu given in const.dek
! *** atomic mass excess for neutron and proton
! *** Audi & Wapstra 1995 (by Frank July 2008)
      deltap = 7.288969d0  ! Mev
      deltan = 8.071323d0    !Mev


!..zero all the isotope pointers
      do i=1,nisotp
       isotp(i)   = 0
      enddo


      if(iniab == 1)then

!..set the size of the network and the number of rates
       ionmax  = 47

!..set the id numbers of the elements
! ***  the part below is obsolete

      ihe3  = 1
      ihe4  = 2
      ic12  = 3
      ic13  = 4
      in13  = 5
      in14  = 6
      io14  = 7
      io15  = 8
      io16  = 9
      io17  = 10
      io18  = 11
      if17  = 12
      if18  = 13
      if19  = 14
      ine18 = 15
      ine19 = 16
      ine20 = 17
      img22 = 18
      img24 = 19
      ial27 = 20
      isi28 = 21
      ip31  = 22
      is30  = 23
      is32  = 24
      icl35 = 25
      iar36 = 26
      ik39  = 27
      ica40 = 28
      iti44 = 29
      icr48 = 30
      icr49 = 31
      icr50 = 32
      icr51 = 33
      icr52 = 34
      icr53 = 35
      icr54 = 36
      ife52 = 37
      ife54 = 38
      ife55 = 39
      ife56 = 40
      ife57 = 41
      ife58 = 42
      ico55 = 43
      ini56 = 44
      ini58 = 45
      ineut = 46
      iprot = 47


!..set the names of the elements
      ionam(ihe3)  = 'HE  3'
      ionam(ihe4)  = 'HE  4'
      ionam(ic12)  = 'C  12'
      ionam(ic13)  = 'C  13'
      ionam(in13)  = 'N  13'
      ionam(in14)  = 'N  14'
      ionam(io14)  = 'O  14'
      ionam(io15)  = 'O  15'
      ionam(io16)  = 'O  16'
      ionam(io17)  = 'O  17'
      ionam(io18)  = 'O  18'
      ionam(if17)  = 'F  17'
      ionam(if18)  = 'F  18'
      ionam(if19)  = 'F  19'
      ionam(ine18) = 'NE 18'
      ionam(ine19) = 'NE 19'
      ionam(ine20) = 'NE 20'
      ionam(img22) = 'MG 22'
      ionam(img24) = 'MG 24'
      ionam(ial27) = 'AL 27'
      ionam(isi28) = 'SI 28'
      ionam(ip31)  = 'P  31'
      ionam(is30)  = 'S  30'
      ionam(is32)  = 'S  32'
      ionam(icl35) = 'CL 35'
      ionam(iar36) = 'AR 36'
      ionam(ik39)  = 'K  39'
      ionam(ica40) = 'CA 40'
      ionam(iti44) = 'TI 44'
      ionam(icr48) = 'CR 48'
      ionam(icr49) = 'CR 49'
      ionam(icr50) = 'CR 50'
      ionam(icr51) = 'CR 51'
      ionam(icr52) = 'CR 52'
      ionam(icr53) = 'CR 53'
      ionam(icr54) = 'CR 54'
      ionam(ife52) = 'FE 52'
      ionam(ife54) = 'FE 54'
      ionam(ife55) = 'FE 55'
      ionam(ife56) = 'FE 56'
      ionam(ife57) = 'FE 57'
      ionam(ife58) = 'FE 58'
      ionam(ico55) = 'CO 55'
      ionam(ini56) = 'NI 56'
      ionam(ini58) = 'NI 58'
      ionam(ineut) = 'NEUT '
      ionam(iprot) = 'PROT '



!..set the number of nucleons in the element
      aion(ihe3)  = 3.0d0
      aion(ihe4)  = 4.0d0
      aion(ic12)  = 12.0d0
      aion(ic13)  = 13.0d0
      aion(in13)  = 13.0d0
      aion(in14)  = 14.0d0
      aion(io14)  = 14.0d0
      aion(io15)  = 15.0d0
      aion(io16)  = 16.0d0
      aion(io17)  = 17.0d0
      aion(io18)  = 18.0d0
      aion(if17)  = 17.0d0
      aion(if18)  = 18.0d0
      aion(if19)  = 19.0d0
      aion(ine18) = 18.0d0
      aion(ine19) = 19.0d0
      aion(ine20) = 20.0d0
      aion(img22) = 22.0d0
      aion(img24) = 24.0d0
      aion(ial27) = 27.0d0
      aion(isi28) = 28.0d0
      aion(ip31)  = 31.0d0
      aion(is30)  = 30.0d0
      aion(is32)  = 32.0d0
      aion(icl35) = 35.0d0
      aion(iar36) = 36.0d0
      aion(ik39)  = 39.0d0
      aion(ica40) = 40.0d0
      aion(iti44) = 44.0d0
      aion(icr48) = 48.0d0
      aion(icr49) = 49.0d0
      aion(icr50) = 50.0d0
      aion(icr51) = 51.0d0
      aion(icr52) = 52.0d0
      aion(icr53) = 53.0d0
      aion(icr54) = 54.0d0
      aion(ife52) = 52.0d0
      aion(ife54) = 54.0d0
      aion(ife55) = 55.0d0
      aion(ife56) = 56.0d0
      aion(ife57) = 57.0d0
      aion(ife58) = 58.0d0
      aion(ico55) = 55.0d0
      aion(ini56) = 56.0d0
      aion(ini58) = 58.0d0
      aion(ineut) = 1.0d0
      aion(iprot) = 1.0d0




!..set the number of protons in the element
      zion(ihe3)  = 2.0d0
      zion(ihe4)  = 2.0d0
      zion(ic12)  = 6.0d0
      zion(ic13)  = 6.0d0
      zion(in13)  = 7.0d0
      zion(in14)  = 7.0d0
      zion(io14)  = 8.0d0
      zion(io15)  = 8.0d0
      zion(io16)  = 8.0d0
      zion(io17)  = 8.0d0
      zion(io18)  = 8.0d0
      zion(if17)  = 9.0d0
      zion(if18)  = 9.0d0
      zion(if19)  = 9.0d0
      zion(ine18) = 10.0d0
      zion(ine19) = 10.0d0
      zion(ine20) = 10.0d0
      zion(img22) = 12.0d0
      zion(img24) = 12.0d0
      zion(ial27) = 13.0d0
      zion(isi28) = 14.0d0
      zion(ip31)  = 15.0d0
      zion(is30)  = 16.0d0
      zion(is32)  = 16.0d0
      zion(icl35) = 17.0d0
      zion(iar36) = 18.0d0
      zion(ik39)  = 19.0d0
      zion(ica40) = 20.0d0
      zion(iti44) = 22.0d0
      zion(icr48) = 24.0d0
      zion(icr49) = 24.0d0
      zion(icr50) = 24.0d0
      zion(icr51) = 24.0d0
      zion(icr52) = 24.0d0
      zion(icr53) = 24.0d0
      zion(icr54) = 24.0d0
      zion(ife52) = 26.0d0
      zion(ife54) = 26.0d0
      zion(ife55) = 26.0d0
      zion(ife56) = 26.0d0
      zion(ife57) = 26.0d0
      zion(ife58) = 26.0d0
      zion(ico55) = 27.0d0
      zion(ini56) = 28.0d0
      zion(ini58) = 28.0d0
      zion(ineut) = 0.0d0
      zion(iprot) = 1.0d0


!..set the number of neutrons and niso1 index
       do i=1,ionmax
        nion(i) = aion(i) - zion(i)
        niso1(int(aion(i))+1,int(zion(i))+2,1) = i
       enddo


!..set the binding energy of the element
! *** I have now reaclib bion.
! *** I leave here if for v&v one want to use them
! *** for comparison.

      bion(ihe3)  = 7.71819d0
      bion(ihe4)  = 28.29603d0
      bion(ic12)  = 92.16294d0
      bion(ic13)  = 97.1060d0
      bion(in13)  = 94.1030d0
      bion(in14)  = 104.65998d0
      bion(io14)  = 98.7310d0
      bion(io15)  = 111.9530d0
      bion(io16)  = 127.62093d0
      bion(io17)  = 131.7600d0
      bion(io18)  = 139.8040d0
      bion(if17)  = 128.2170d0
      bion(if18)  = 137.3670d0
      bion(if19)  = 147.7980d0
      bion(ine18) = 132.1390d0
      bion(ine19) = 143.7780d0
      bion(ine20) = 160.64788d0
      bion(img22) = 168.5750d0
      bion(img24) = 198.2579d0
      bion(ial27) = 224.9480d0
      bion(isi28) = 236.5379d0
      bion(ip31)  = 262.9120d0
      bion(is30)  = 243.6810d0
      bion(is32)  = 271.7825d0
      bion(icl35) = 298.2050d0
      bion(iar36) = 306.7202d0
      bion(ik39)  = 333.7180d0
      bion(ica40) = 342.0568d0
      bion(iti44) = 375.4772d0
      bion(icr48) = 411.469d0
      bion(icr49) = 422.0370d0
      bion(icr50) = 435.0370d0
      bion(icr51) = 444.2980d0
      bion(icr52) = 456.3370d0
      bion(icr53) = 464.2760d0
      bion(icr54) = 473.9950d0
      bion(ife52) = 447.708d0
      bion(ife54) = 471.7696d0
      bion(ife55) = 481.0480d0
      bion(ife56) = 492.2450d0
      bion(ife57) = 499.8910d0
      bion(ife58) = 509.9350d0
      bion(ico55) = 476.8150d0
!      bion(ini56) = 465.98788d0 !value from reaclib JINA
      bion(ini56) = 484.003d0   ! value original from torch
      bion(ini58) = 506.4450d0
      bion(ineut) = 0.0d0
      bion(iprot) = 0.0d0


   else if(iniab == 2)then


      !..set the size of the network and niso1 according to ppn
      j = 0
      do i=1,nsp
         if(considerisotope(i))then
            j = j+1
            if(j.gt.0)then
               ionam(j) = zis(i)
               aion(j)  = an(i)
               zion(j)  = zn(i)
               nion(j)  = an(i) - zn(i)
               isomer_state(j) = isomeric_state(i)
               niso1(int(aion(j))+1,int(zion(j))+2,isomer_state(j)) = j
            end if
         end if
      end do
      ionmax = j
   else
      stop 'unknown iniab function value'
   end if

! *** here I calculate partition functions and binding energies from
! *** reaclib_basel for the correct temperature.
! *** please consider that at this point I already have many info for iniab = 1.
! *** Here below I mainly define many things for iniab = 2.

   if(ipartition == 1)then

! *** in the case below bion is not defined yet.
!    if(iniab == 2)stop'new network with old partition function'

!          do i=1,ionmax
!           wpart(i) = 1.0d0
!          enddo

         do i=1,nnpartdim
          if (istate(i).ge.1)then
! *** I assume that therm/ground/isomer has the same partition function. 
! *** Is this true? Probably not, since the spin is different....
! *** However, I think it would be better to leave isomers, see e.g. 
! *** to follow production of Ta180.
! *** Marco 30 Sep 2010 - updated 7 May 2012.
          j = niso1(anumpart(i)+1,znumpart(i)+2,istate(i))
           if(j.gt.0)then
            wpart(j) = 1.0d0 * (2.d0*spinpart(i)+1.0d0)
            if(j.gt.ionmax)exit
           end if
          end if
         end do



!..set the partition functions
!..these are generally temperature dependent
!..here i'll just take the ground state
         if(iniab == 1)then
          wpart(ineut) = 2.0d0
          wpart(iprot) = 2.0d0
         else if (iniab == 2) then
          wpart(ispe('NEUT ')) = 2.0d0
          wpart(ispe('PROT ')) = 2.0d0
         end if

        else if (ipartition == 2)then
! *** I include wpart = 1. as default, as it is done
! *** for ipartition = 1
         wpart = 1.0d0

! ^_^ INTERPOLATE THE PARTITION FUNCTIONS AT THE GIVEN TEMPERATURE
        call reaclib_interpolate_partition_functions(t9)

! *** partinter is the partition function at the temperature t9 of the isotope
! *** (Z_i,A_i). Now using niso1 index I find out the partition function for
! *** each isotope in the network.

         do i=1,nnpartdim
          if (istate(i).ge.1)then
! *** I assume that therm/ground/isomer has the same partition function. 
! *** Is this true? Probably not, since the spin is different....
! *** However, I think it would be better to leave isomers, see e.g. 
! *** to follow production of Ta180.
! *** Marco 30 Sep 2010 - updated 7 May 2012.
          j = niso1(anumpart(i)+1,znumpart(i)+2,istate(i))
           if(j.gt.0)then
            wpart(j) = partinter(i) * (2.d0*spinpart(i)+1.0d0)
!           exion(j) = exmass(i)
            if(j.gt.ionmax)exit
           end if
          end if
         end do

         if(iniab == 1)then
          wpart(ineut) = 2.0d0
          wpart(iprot) = 2.0d0
         else if (iniab == 2) then
          wpart(ispe('NEUT ')) = 2.0d0
          wpart(ispe('PROT ')) = 2.0d0
         end if


       else
         stop 'unknown ipartition function value'
       end if


! *** note: binding energy = DM(amu) = Z*m_p + (A-Z)*m_n - M_tot
! *** DM(amu) * 931.494 (Mev/amu) = DM(MeV)
! *** DM(MeV)/A --> binding energy per nucleon

! *** definition of binding energies
! *** if iniab = 1 I already have them!

       if (iniab == 2)then
!   print *, 'passa di qui'
        j = 0
        do i=1,nnpartdim
         !print*,'istate', istate(i)
    if (istate(i).ge.1)then
! *** I assume that thermlaized/gr/isomeric states have the same binding energies.
! *** in case there are isomers in these conditions. Marco Sep 2010
          !print *,niso1(anumpart(i)+1,znumpart(i)+2,istate(i))
          j = niso1(anumpart(i)+1,znumpart(i)+2,istate(i))
     if(j.gt.0)then
           !print *,'qui:',anumpart(i),exmass(i)
      exion(j) = exmass(i)
!           print *,'j',j,ionmax,
!     1 niso1(anumpart(i)+1,znumpart(i)+2,istate(i)),
!     2 anumpart(i),znumpart(i),istate(i)
      if(j.gt.ionmax)exit
     end if
    end if
   end do


        do i=1,ionmax
! *** this provides really similar with what Frank has
    bion(i) = zion(i)*deltap + (aion(i) - zion(i))*deltan - exion(i)
!         print*,'bion: ',bion(i),ionam(i) ,'exion',exion(i)
! *** I calculated this, but is a little bit different than
! *** given from Frank in original network
!         bion(i) = (- (aion(i) + exion(i)/fact1) + zion(i)*m_p +
!     1    (aion(i) - zion(i))*m_n) * fact1
!         print*,'bion: ',bion(i),ionam(i) ,'exion',exion(i)
        end do
!        stop 'eccomi'
! ***  I set bion for neut and prot = 0.
         bion(ispe('NEUT ')) = 0.0d0
         bion(ispe('PROT ')) = 0.0d0
       end if

! *** here below I define imion. Assumtion mion(i) = aion(i)*amu
! *** or I consider mass formula ? exion (mass excess) below
! *** is taken from FRDM formula.
! *** note: if mion(i) = aion(i)*amu, I am assuming also that m_p = m_n!
! *** note: in the mion formula of mion it must be +exion and not -exion...
! *** (check the file winvn in REACLIB)
      if(imion == 1)then
       do i=1,ionmax
        mion(i) = aion(i) * amu
       end do
      else if(imion == 2)then
!   print *,'passo di qui'
       do i=1,ionmax

! *** mion here is calculated
   !print *,aion(i),exion(i),fact1,amu
!    mion(i) = (aion(i) + exion(i)/fact1) * amu
       mion(i) = nion(i)*mn + zion(i)*mp - bion(i)*mev2gr
       ! ^_^:
!       mion(i) =  nion(i)*mn + zion(i)*mp - ( aion(i) * amu + exion(i) * mev2gr - zion(i) * 0.511d0*mev2gr )
       end do
      else
         stop 'unknown imion function value'
      end if


! *** here below I define screening.
!..set the coulomb corrections
      mu_c_p = 0.0d0
      do i=i,ionmax
       mu_c(i) = 0.0d0
      enddo
!..calder et al, apj 656 313 2007, eq a1
!..number density of free electrons from matter
       beta    = 1.0d0/(kerg * t9 * 1.d+9)
       xne    = ye_want * avo * den
       ge     = esqu * beta * (forthpi * xne)**third
       a3     = -0.5d0*sqrt(3.0d0) - a1 / sqrt(a2)

      if(icoulomb == 1)then
   continue
      else if(icoulomb == 2)then
       do i=1,ionmax
        gi      = zion(i)**(fivth)  * ge
        sqrtgi  = sqrt(gi)
        mu_c(i) = a1*(sqrt(gi*(a2+gi)) - a2*log(sqrt(gi*a2inv) + sqrt(1.0d0 + gi*a2inv)) + &
              2.0d0*a3*(sqrtgi - atan(sqrtgi)))
!   print *,'mu_c',mu_c(i),ionam(i)
       enddo
       mu_c_p = mu_c(niso1(1+1,1+2,1))
      else
         stop 'unknown icoulomb function value'
      end if
!   print *,'mu_c_p =',mu_c_p,niso1(1+1,1+2)
! *** printing for debugging
!        do i=1,ionmax
!         print *,i,'ionam:',ionam(i),zion(i),aion(i)
!    print *,i,'bion',bion(i),'exion',exion(i)
!    print *,i,'wpart',wpart(i),'mion',mion(i)
!    print *,i,'mu_c',mu_c(i)
!        enddo
!   stop

      return
      end subroutine init_network

end module nse_torch
