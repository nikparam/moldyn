!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     -----------------
!c     Moleular dynamics
!c     -----------------
!c         Part II
!c
!c     This program integrates the Newton's eq. in one time step,
!c     by means of the Velocity Verlet Algoritim (VVA).
!c
!c     This second part calculates the velocity.
!c
!c     Mario Barbatti, March 2005.
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
   program moldyn02
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)

      real(kind=dpr), dimension(:,:), allocatable :: x,v,g,a,wf
      real(kind=dpr), dimension(:)  , allocatable :: aM,Z,Epot
      character(2), dimension(:), allocatable :: S
      common /first/ Nat,istep,nstat,ndamp,nintc 
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv

      open(unit= 7, file='control.d',status='unknown')
      open(unit= 8, file='geom' ,    status='unknown')
      open(unit= 9, file='veloc',    status='unknown')
      open(unit=10, file='grad' ,    status='unknown')
      open(unit=11, file='epot' ,    status='unknown')
!      open(unit=12, file='auxp' ,    status='old')
!      open(unit=13, file='etot0',    status='unknown')
!      open(unit=14, file='last' ,    status='unknown')
!      open(unit=15, file='../RESULTS/dyn.out',  status='unknown',position='append')
!      open(unit=16, file='../RESULTS/dyn.mld',  status='unknown',position='append')
!      open(unit=17, file='worse',    status='unknown')
!      open(unit=18, file='../RESULTS/en.dat' ,  status='unknown',position='append')
!      open(unit=21, file='typeofdyn.out',status='unknown')
!      open(unit=22, file='nhopsold', status='unknown')
      open(unit=23, file='wfrun',    status='unknown')

      read(7,*) Nat,istep,nstat,nstatdyn,ndamp,kt,dt,t,tmax,nintc,mem,&
                nxrestart,thres,killstat,timekill,prog,lvprt,Etot_jump,Etot_drift
      dt=dt/timeunit

      allocate (x(Nat,3),v(Nat,3),g(Nat,3),a(Nat,3),STAT=istat)
      if (istat /= 0)  STOP "*** Not enough memory ***"
      allocate (aM(Nat),Z(Nat),Epot(nstat),S(Nat),wf(nstat,2),STAT=istat2)
      if (istat2 /= 0) STOP "*** Not enough memory ***"

! ..  Read files
      call readfiles(x,v,g,a,aM,S,Z,Epot,wf)

! ..  Integrate in one time step
      call velverlet(x,v,g,a,aM)

! ..  Write velocity
      call writefiles(x,v,aM,S,Z,Etot,Ekin,Epot)

      deallocate (x,v,g,a)
      deallocate (aM,Z,Epot,S,wf)

  100 format(A2,F8.1,4F14.8)
    end

! .............................. READ

      subroutine readfiles(x,v,g,a,aM,S,Z,Epot,wf)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      dimension x(nat,3),v(nat,3),g(nat,3),a(nat,3)
      dimension aM(nat),Z(nat),Epot(nstat),wf(nstat,2)
      character(2) S(nat)
      common /first/ Nat,istep,nstat,ndamp,nintc 
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv

! ... read geom (a.u.)
      do n=1,nat
      read(8,*) S(n),Z(n),(x(n,i),i=1,3),aM(n)
      aM(n)=aM(n)*proton
      enddo

! ... read veloc (a.u.)
      open(24,file='oldveloc',status='unknown')
      do n=1,nat
      read(9,*) (v(n,i), i=1,3)
      !
      write(24,*) (v(n,i), i=1,3)  ! write old veloc to interpolate
      enddo
      close(24)

! ... read grad (a.u.)
      do n=1,nat
      read(10,*) (g(n,i), i=1,3)
      enddo

! ... read potential energy (a.u.)
      do k=1,nstat
      read(11,*) Epot(k)
      enddo

! ... aceleration
      do n=1,nat
      do i=1,3
      a(n,i)=-(1.0_dpr/aM(n))*g(n,i)
      enddo
      enddo

      do k=1,nstat                                    ! real and imaginary parts of the wave function
      read(23,*) wf(k,1),wf(k,2)
      enddo

      endsubroutine

! ............................... INTEGRATION

      subroutine velverlet(x,v,g,a,aM)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      dimension x(nat,3),v(nat,3),g(nat,3),a(nat,3),aM(nat)
      common /first/ Nat,istep,nstat,ndamp,nintc 
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv

      do n=1,nat
      do i=1,3

      if (ndamp .eq. 0) then                   ! ndamp different from 0: damped dynamics
      v(n,i) = v(n,i) + 0.5_dpr*a(n,i)*dt        ! v(t+dt)
      endif

      enddo
      enddo

      endsubroutine

! ............................... WRITEFILES

      subroutine writefiles(x,v,aM,S,Z,Etot,Ekin,Epot)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      dimension x(nat,3),v(nat,3),g(nat,3),a(nat,3)
      dimension aM(nat),Z(nat),Epot(nstat),wf(nstat,2)
      character(2) S(nat)
      common /first/ Nat,istep,nstat,ndamp,nintc 
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv

! ... write veloc (a.u.)
      open(24,file='newveloc',status='unknown')
      rewind 9
      do n=1,nat
      write(9,*) (v(n,i), i=1,3)
      !
      write(24,*) (v(n,i), i=1,3)  ! write new veloc to interpolate
      enddo
      close(24)

      endsubroutine

