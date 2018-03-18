!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     -----------------
!c     Moleular dynamics
!c     -----------------
!c         Part IV
!c
!c     This program integrates the Newton's eq. in one time step,
!c     by means of the Velocity Verlet Algoritim (VVA).
!c
!c     This fourth part calculates the velocity after hopping
!c     in the interpolation process
!c
!c     After hopping, it updates the gradient file to continue the 
!c     timestep in the right surface. Besides that it projects the 
!c     velocities to the first and last sub-timestep on the new 
!c     surface.
!c
!c     Mario Barbatti, May-Dec 2006.
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
   program moldyn04
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      real*8, dimension(:,:), allocatable :: x,v,g,a,h,gi,gf,vi,vf
      real*8, dimension(:)  , allocatable :: aM,Z,Epot
      character(2), dimension(:), allocatable :: S
      common /first/ Nat,istep,nstat,ndamp,nintc 
      common /second/ dt,t,tmax,thres,adjmom
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv
      common /fifth/ Ms,j,irk,mom,newsurf,noldsurf

      open(unit= 7,file='control.d',status='unknown')
      open(unit=14,file='md04.inp', status='unknown')
      open(unit=20,file='moldyn04.log',status='unknown',&
      position='append')
      read(7,*) Nat,istep,nstat,nstatdyn,ndamp,kt,dt,t,tmax,nintc,mem,&
            nxrestart,thres,killstat,timekill,prog,lvprt,Etot_jump,Etot_drift
      dt=dt/timeunit

      read(14,*) Ms,j,irk,mom,adjmom,newsurf,noldsurf

      close(7)
      close(14)

      allocate (x(Nat,3),v(Nat,3),g(Nat,3),a(Nat,3),STAT=istat)
      if (istat /= 0)  STOP "*** Not enough memory ***"
      allocate (vi(Nat,3),vf(Nat,3),gi(Nat,3),gf(Nat,3),STAT=istat1)
      if (istat1 /= 0) STOP "*** Not enough memory ***"
      allocate (aM(Nat),Z(Nat),Epot(nstat),S(Nat),h(Nat,3),STAT=istat2)
      if (istat2 /= 0) STOP "*** Not enough memory ***"

! ..  Read files
      call readfiles(x,v,g,a,aM,S,Z,vi,vf,Epot,h)

! ..  Get velocities
      call velverlet(x,v,g,a,aM,vi,vf)

! ..  Write velocity
      call writefiles(x,v,aM,S,Z,vi,vf)

      deallocate (x,v,g,a)
      deallocate (aM,Z,Epot,S)

  100 format(A2,F8.1,4F14.8)
    end program moldyn04

! .............................. READ

      subroutine readfiles(x,v,g,a,aM,S,Z,vi,vf,Epot,h)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      dimension x(nat,3),v(nat,3),g(nat,3),a(nat,3)
      dimension gi(nat,3),gf(nat,3),h(nat,3)
      dimension aM(nat),Z(nat),Epot(nstat),wf(nstat,2),vi(nat,3),vf(nat,3)
      character(2) S(nat)
      common /first/ Nat,istep,nstat,ndamp,nintc 
      common /second/ dt,t,tmax,thres,adjmom
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv
      common /fifth/ Ms,j,irk,mom,newsurf,noldsurf

      if ((irk == 2) .and. (mom == -1)) then
! ... read veloc (a.u.)
        open(unit=12, file='oldveloc', status='unknown')
        open(unit=13, file='newveloc', status='unknown')
        do n=1,nat
           read(12,*) (vi(n,i), i=1,3)
           read(13,*) (vf(n,i), i=1,3)
        enddo
        close(12)
        close(13)
      endif
      
      if (irk == 1) then
! ... read geom (a.u.)
        open(unit= 8, file='geom' ,    status='unknown')
        do n=1,nat
          read(8,*) S(n),Z(n),(x(n,i),i=1,3),aM(n)
          aM(n)=aM(n)*proton
        enddo
        close(8)

! ... read veloc (a.u.)
        if (j == Ms) then
          write(20,*) "istep j Ms",istep,j,Ms
          open(unit= 9, file='veloc.bk',status='unknown')
        else
          open(unit= 9, file='veloc',   status='unknown')
        endif
        do n=1,nat
          read(9,*) (v(n,i), i=1,3)
          write(20,*) (v(n,i), i=1,3)
        enddo
        close(9)

! ... kinetic energy (a.u.)
        Ekin = 0.0_dpr
        do n=1,Nat
          do i=1,3
            Ekin = Ekin + aM(n)*v(n,i)**2
          enddo
        enddo
        Ekin = 0.5_dpr*Ekin
 
! ... read potential energy (a.u.)
        open(unit=15, file='epot.bk', status='unknown')
        do k=1,nstat
          read(15,*) Epot(k)
        enddo
        close(15)

! ... potential energy difference
        delu = Epot(newsurf) - Epot(noldsurf) ! current - previous
        discri = 1.0_dpr - delu/Ekin
        write(20,*) "istep discri, delu, Ekin:",istep,discri,delu,Ekin

! ... read old grad (a.u.)
        open(unit=11, file='grad.all' ,status='unknown')
        if (noldsurf > 1) then
           do n=1,nat*(noldsurf-1)
             read(11,*) 
           enddo
        endif

        do n=1,nat
          read(11,*) (gi(n,i), i=1,3)
        enddo
        close(11)

! ... read new grad (a.u.)
        open(unit=11, file='grad.all' ,status='unknown')
        rewind 11
        if (newsurf > 1) then
           do n=1,nat*(newsurf-1)
             read(11,*) 
           enddo
        endif

        do n=1,nat
          read(11,*) (gf(n,i), i=1,3)
        enddo
        close(11)

! ... update gradient
        open(unit=11, file='grad' ,status='unknown')
        do n=1,nat
          write(11,*) (gf(n,i), i=1,3)
        enddo
        close(11)

! ... new acceleration
        if (j /= Ms) then
          do i=1,3
            a(:,i) = -gf(:,i)/aM(:)
          enddo
        endif

! ... gradient difference
        g = (gf-gi)/2.0_dpr

! ... Nonadiabatic couplings
!     l -> k (noldsurf -> newsurf)
        l = noldsurf
        k = newsurf
        if(k.gt.l) then
          kl=((k-2)*(k-1))/2+l
        else if (k.lt.l) then
          kl=((l-2)*(l-1))/2+k
        endif

        open(unit=16, file='nad_vectors.bk', status='unknown')
        if (kl > 1) then
           do n=1,nat*(kl-1)
             read(16,*)
           enddo
        endif

        do n=1,nat
          read(16,*) (h(n,i), i=1,3)
        enddo
        close(16)

! ... Norms
        gn = 0.0_dpr
        hn = 0.0_dpr
        do n = 1,nat
          do i = 1,3
            gn = gn + g(n,i)**2
            hn = hn + h(n,i)**2
          enddo
        enddo      
        gn = dsqrt(gn)
        hn = dsqrt(hn)

! ... adjustment direction
        if (j == Ms) then
          write(20,*) "j = Ms: adjusting momentum."
          if (adjmom < 0) then
            do i = 1,3
              a(:,i) = aM(:)*v(:,i)
            enddo
          else
            a = dsin(adjmom*deg2rad)*g/gn + dcos(adjmom*deg2rad)*h/hn
          endif

! ... change momentum
!     In some very rare cases, discri can be > 0 at the time
!     of hopping, but < 0 at the end of the interpolation 
!     period. When this happens, the hopping is allowed, but
!     is not possible to rescale the energy here.

          if (discri >= 0.0_dpr) then
            call varv(discri,delu,1.0_dpr,1,a,v,nat,aM)
          else
            if (j == Ms) then
              open(unit= 9, file='veloc.bk', status='unknown')
            else
              open(unit= 9, file='veloc', status='unknown')
            endif
            do n=1,nat
              read(9,*) (v(n,i), i=1,3)
              ! >>>>>>
              write(20,*) "discri < 0 ..."
              v(n,:)=v(n,:)*discri
              ! >>>>>>
            enddo
            close(9)
          endif

! ... new kinetic energy (a.u.)
          Ekin1 = 0.0_dpr
          do n=1,Nat
            do i=1,3
              Ekin1 = Ekin1 + aM(n)*v(n,i)**2
            enddo
          enddo
          Ekin1 = 0.5_dpr*Ekin1

! ... new total energy:
          Etot  = Ekin  + Epot(noldsurf)
          Etot1 = Ekin1 + Epot(newsurf)

        write(20,'(/,A)') &
     "      Ekin old          Ekin new          Etot old          Etot new"
        write(20,'(4F18.10)') Ekin, Ekin1, Etot, Etot1
 
        endif
      endif

      endsubroutine

! ............................... INTEGRATION

      subroutine velverlet(x,v,g,a,aM,vi,vf)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      dimension x(nat,3),v(nat,3),g(nat,3),a(nat,3),aM(nat),vi(nat,3),vf(nat,3)
      common /first/ Nat,istep,nstat,ndamp,nintc 
      common /second/ dt,t,tmax,thres,adjmom
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv
      common /fifth/ Ms,j,irk,mom,newsurf,noldsurf

      if ((irk == 2) .and. (mom == -1)) then
         vi = -vi
         vf = -vf
      endif

! ... projected velocities
      if (j /= Ms) then
        if (irk == 1) then
           vi = v - a*j*dt/Ms
           vf = v + a*(Ms-j)*dt/Ms
        endif
      endif

      endsubroutine

! ............................... WRITEFILES

      subroutine writefiles(x,v,aM,S,Z,vi,vf)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      dimension x(nat,3),v(nat,3),g(nat,3),a(nat,3),vi(nat,3),vf(nat,3)
      dimension aM(nat),Z(nat),Epot(nstat),wf(nstat,2)
      character(2) S(nat)
      common /first/ Nat,istep,nstat,ndamp,nintc 
      common /second/ dt,t,tmax,thres,adjmom
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv
      common /fifth/ Ms,j,irk,mom,newsurf,noldsurf

! ... write veloc (a.u.)
      if (irk == 1) then
        open(unit=17, file='veloc', status='unknown')
        open(unit=18, file='veloc.adj', status='unknown')
        do n=1,nat
          write(17,*) (v(n,i), i=1,3)
          write(18,*) (v(n,i), i=1,3)
        enddo
        close(17)
        close(18)
      endif

      if (((irk == 2).and.(mom == -1)).or.(irk==1)) then
        open(unit=12, file='oldveloc', status='unknown')
        open(unit=13, file='newveloc', status='unknown')
        do n=1,nat
          write(12,*) (vi(n,i), i=1,3)
          write(13,*) (vf(n,i), i=1,3)  ! write new veloc to interpolate
        enddo
        close(12)
        close(13)
      endif

      endsubroutine

! ............................... CHANGE MOMENTUM

      subroutine varv(discri,delu,bb,aa,a,v,nat,aM)
! Change velocity after a hopping or after a frustated hopping
! Usage varv(discri,delu,bb,aa,vec,v,numat,ams) where
! discri = 1+DE/Ekin, DE = Epot(before) - Epot(after)
! delu   = -DE
! bb     = 1 (allowed hopping); -1 (forbiden)
! aa     = 1 (keep momentum); -1 (invert momentum) (only for bb=-1)
! a      = vector along which the momentum will be adjusted
! v      = velocity vector
! nat    = number of atoms
! aM     = atomic masses
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      dimension a(nat,3),v(nat,3),aM(nat)
      integer :: aa

      if (bb > 0.0_dpr) then
      ! Hopping is allowed: get new velocities
         dep = -delu
         c = 0.0_dpr
         b = 0.0_dpr
         do n = 1,nat
            do i = 1,3
               c = c + a(n,i)*a(n,i)/aM(n)
               b = b + v(n,i)*a(n,i)
            enddo
         enddo
         Delta = b**2 + 2.0_dpr*c*dep
         if (Delta >= 0.0_dpr) then
       ! Delta >= 0: adjust new veloc along vec
          sgm1  = 1.0_dpr/c*(-b+dsqrt(Delta))
          sgm2  = 1.0_dpr/c*(-b-dsqrt(Delta))
          alpha = min(dabs(sgm1),dabs(sgm2))
          eps = 1E-6_dpr
          if (dabs(alpha-dabs(sgm1)) < eps) alpha = sgm1
          if (dabs(alpha-dabs(sgm2)) < eps) alpha = sgm2
          !write(20,*) "dep,c,b,D,s1,s2,ap= ",dep,c,b,Delta,sgm1,sgm2,alpha
          do i = 1,3
             v(:,i) = v(:,i) + alpha*a(:,i)/aM(:)
          enddo
       else
       ! Delta < 0: only rescale velocities
          discri = dsqrt(discri)
          alpha = (discri-bb)
          !write(20,*) "dep,c,b,D,ap= ",dep,c,b,Delta,alpha
             v = v + alpha*v
       endif
      end if
      end subroutine varv
