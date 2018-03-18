!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     -----------------
!c     Moleular dynamics
!c     -----------------
!c         Part I
!c
!c     This program integrates the Newton's eq. in one time step,
!c     by means of the Velocity Verlet Algoritim (VVA).
!c
!c     Mario Barbatti, March 2005.
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Part I
!     Compute r(t+dt)
!     Compute v(t+dt/2)
!
!     Part II
!     Compute v(t+dt)
!     Compute Etot, Ekin
!
    program moldyn01
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)

      real(kind=dpr) :: iap
      real(kind=dpr), dimension(:,:), allocatable :: x,v,g,a,wf

      DOUBLE COMPLEX, ALLOCATABLE :: a_coeff(:,:)  ! weight coeffs

      real(kind=dpr), dimension(:)  , allocatable :: aM,Z,Epot
      character(2), dimension(:), allocatable :: S
      common /first/ Nat,istep,nstat,ndamp,nintc
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv
      common /fifth/ t_old,timekill,prog,Etot_jump,Etot_drift
      common /sixth/ killstat,lvprt

      open(unit= 7, file='control.d',status='unknown')
      open(unit= 8, file='geom' ,    status='unknown')
      open(unit= 9, file='veloc',    status='unknown')

      OPEN(UNIT=10, FILE='weights',  STATUS='unknown')

      open(unit=11, file='grad' ,    status='unknown')
      open(unit=12, file='epot' ,    status='unknown')
      open(unit=13, file='auxp' ,    status='old')
      open(unit=14, file='etot0',    status='unknown')
      open(unit=15, file='last' ,    status='unknown')
      open(unit=16, file='../RESULTS/dyn.out',  status='unknown',position='append')
!       open(unit=16, file='../RESULTS/dyn.mld',  status='unknown',position='append')
      open(unit=17, file='worse',    status='unknown')
      open(unit=18, file='../RESULTS/en.dat' ,  status='unknown',position='append')
      open(unit=21, file='typeofdyn.out',status='unknown')
      open(unit=22, file='nhopsold', status='unknown')
      open(unit=23, file='wfrun',    status='unknown')

      read(7,*) Nat,istep,nstat,nstatdyn,ndamp,kt,dt,t,tmax,nintc,mem,&
                nxrestart,thres,killstat,timekill,prog,lvprt,Etot_jump,Etot_drift
      dt=dt/timeunit

      allocate (x(Nat,3),v(Nat,3),g(Nat,3),a(Nat,3),STAT=istat)
      if (istat /= 0)  STOP "*** Not enough memory ***"
      allocate (aM(Nat),Z(Nat),Epot(nstat),S(Nat),wf(nstat,2),STAT=istat2)
      if (istat2 /= 0) STOP "*** Not enough memory ***"

      ALLOCATE( a_coeff(nstat, nat, 3) ) ! allocate weight coeff

! ..  Read files
!     call readfiles(x,v,g,a,aM,S,Z,Epot,wf)
      call readfiles(x,v,a_coeff,g,a,aM,S,Z,Epot,wf)

      if (istep .eq. 0) then

! ..    Printing information
! mbr        iap=1
! mbr        write(12,*) iap

! ..    Get total initial energy
        call totalenergy(v,a_coeff,aM,Epot,Ekin,Etot)

! ..    Write general initial output
        call writedyn(x,v,a_coeff,aM,S,Z,Etot,Ekin,Epot,wf)

      endif

! ..  Integrate in one time step
      call velverlet(x,v,g,a,aM)
      call coeff_int(a_coeff,v,Epot)

! ..  Write input for the next step calculation
      call writefiles(x,v,a_coeff,aM,S,Z,Etot,Ekin,Epot)

! ..  Write general output
      backspace(12)
      read(12,*) iap
! mbr      if (iap .eq. kt) then
      if (iap .lt. 1.d-9) then
        call writedyn(x,v,a_coeff,aM,S,Z,Etot,Ekin,Epot,wf)
      endif

      deallocate (x,v,g,a,a_coeff)
      deallocate (aM,Z,Epot,S,wf)

  100 format(A2,F8.1,4F14.8)

    end program moldyn01

! .............................. READ

      subroutine readfiles(x,v,a_coeff,g,a,aM,S,Z,Epot,wf)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      dimension x(nat,3),v(nat,3),g(nat,3),a(nat,3)

      DOUBLE COMPLEX :: a_coeff(nstat,nat,3)

      dimension aM(nat),Z(nat),Epot(nstat),wf(nstat,2)
      character(2) S(nat)
      common /first/ Nat,istep,nstat,ndamp,nintc
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv
      common /fifth/ t_old,timekill,prog,Etot_jump,Etot_drift
      common /sixth/ killstat,lvprt

! ... read type of dynamics
      read(21,*) ntype

! ... read geom (a.u.)
      do n=1,nat
!      read(8,'(A2,F8.1,4F14.8)') S(n),Z(n),(x(n,i),i=1,3),aM(n)
!      read(8,'(1x,a2,2x,f5.1,4f14.8)') S(n),Z(n),(x(n,i),i=1,3),aM(n)
      read(8,*) S(n),Z(n),(x(n,i),i=1,3),aM(n)
      aM(n)=aM(n)*proton
      enddo

! ... read veloc (a.u.)
      do n=1,nat
      read(9,*) (v(n,i), i=1,3)
      enddo

! ... read veloc (a.u.)
      do ns = 1, nstat
      do n=1,nat
      read(10,*) (a_coeff(ns,n,i), i=1,3)
      enddo
      enddo

! ... read grad (a.u.)
      do n=1,nat
      read(11,*) (g(n,i), i=1,3)
      enddo

! ... read potential energy (a.u.)
      do k=1,nstat
      read(12,*) Epot(k)
      enddo

! ... aceleration
      do n=1,nat
      do i=1,3
      a(n,i)=-(1.0_dpr/aM(n))*force
      enddo
      enddo

      do k=1,nstat                                    ! real and imaginary parts of the wave function
      read(23,*) wf(k,1),wf(k,2)
      enddo

!     Read old number of hoppings
      read(22,*) nhopsold,nrejhopsold

      endsubroutine

! ............................... INTEGRATION

      subroutine velverlet(x,v,g,a,aM)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      dimension x(nat,3),v(nat,3), g(nat,3),a(nat,3),aM(nat)
      common /first/ Nat,istep,nstat,ndamp,nintc
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fifth/ t_old,timekill,prog,Etot_jump,Etot_drift
      common /sixth/ killstat,lvprt

      do n=1,nat
      do i=1,3

      x(n,i) = x(n,i)+v(n,i)*dt+0.5_dpr*a(n,i)*dt**2   ! r(t+dt)

      if (ndamp .eq. 0) then                         ! if ndamp different from 0, damped dynamics
      v(n,i) = v(n,i)+0.5_dpr*a(n,i)*dt                ! v(t+dt/2)
      endif

      enddo
      enddo

      endsubroutine

! .............................. Weights integration

	SUBROUTINE coeff_int(a_coeff,V,Epot)
	! ... Nonadiabatic couplings
	!     l -> k (noldsurf -> newsurf)
        l = noldsurf
        k = newsurf
        if(k.gt.l) then
          kl=((k-2)*(k-1))/2+l
        else if (k.lt.l) then
          kl=((l-2)*(l-1))/2+k
        endif

        open(unit=17, file='nad_vectors.bk', status='unknown')
        if (kl > 1) then
           do n=1,nat*(kl-1)
             read(17,*)
           enddo
        endif

        do n=1,nat
          read(17,*) (h(n,i), i=1,3)
        enddo
        close(17)

	END SUBROUTINE coeff_int

! .............................. TOTAL ENERGY

      subroutine totalenergy(v,a_coeff,aM,Epot,Ekin,Etot)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      dimension x(nat,3),v(nat,3),g(nat,3),a(nstat,nat,3),aM(nat),Epot(nstat), vec()
      DOUBLE COMPLEX :: a_coeff(nat,3)
      common /first/ Nat,istep,nstat,ndamp,nintc
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fifth/ t_old,timekill,prog,Etot_jump,Etot_drift
      common /sixth/ killstat,lvprt

! ...... Kinetic energy

      Ekin=0.0_dpr
      do n=1,nat
        v2=( dsqrt( v(n,1)**2+v(n,2)**2+v(n,3)**2 ) )**2
	do ns = 1,nstat
        	TK=0.5_dpr*aM(n)*v2*abs( a_coeff(ns,n,1) )**2*abs( a_coeff(ns,n,2) )**2*abs( a_coeff(ns,n,3) )**2
        	Ekin=Ekin+TK
	enddo
      enddo

! ...... Potential energy

      DO n = 1, nat
      DO ns = 1, nstat
      Etot = abs( a_coeff(ns,n,1) )**2*abs( a_coeff(ns,n,2) )**2*abs( a_coeff(ns,n,3) )**2 * Epot(ns)
      END DO
      END DO
! ...... Total energy

      Etot=Etot+Ekin

      endsubroutine

! ............................... WRITEFILES

      subroutine writefiles(x,v,a_coeff,aM,S,Z,Etot,Ekin,Epot)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      dimension x(nat,3),v(nat,3),g(nat,3),a(nat,3)

	DOUBLE COMPLEX :: a_coeff(nstat,nat,3)

      dimension aM(nat),Z(nat),Epot(nstat)
      character(2) S(nat)
      common /first/ Nat,istep,nstat,ndamp,nintc
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fifth/ t_old,timekill,prog,Etot_jump,Etot_drift
      common /sixth/ killstat,lvprt

! ... write geom (a.u.)
      rewind 8
      do n=1,nat
        write(8,'(1x,a2,2x,f5.1,3f25.14,f14.8)') S(n),Z(n),(x(n,i),i=1,3),&
        aM(n)/proton
      enddo

! ... write veloc (a.u.)
      rewind 9
      do n=1,nat
        write(9,*) (v(n,i), i=1,3)
      enddo

! ... write weights
      rewind 10
      do n=1,nat
      do ns=1,nstat
        write(10,*) (a_coeff(ns,n,i), i=1,3)
      enddo
      enddo

! ... update control.d
      istep=istep+1
      t=t+(dt*timeunit)
      rewind 7
      write(7,205) Nat,istep,nstat,nstatdyn,ndamp,kt,dt*timeunit,t,tmax,&
        nintc,mem,nxrestart,thres,killstat,timekill,prog,lvprt,Etot_jump,Etot_drift
  205 format(6(I8,","),3(F12.4,","),2(I8,","),I8,",",F8.2, &
              ",",I8,",",F12.4,",",F9.1,",",I8,2(",",F8.3))

      endsubroutine

! ............................... WRITEDYN

      subroutine writedyn(x,v,a_coeff,aM,S,Z,Etot,Ekin,Epot,wf)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
!
      logical :: ext, write_vel
      integer :: nindices, istat, istat2, i
      integer, dimension(:), allocatable :: indices_towrite
      logical, dimension(:), allocatable :: write_pos

	DOUBLE COMPLEX :: a_coeff(nstat,nat,3)

      dimension x(nat,3),v(nat,3),g(nat,3),a(nat,3)
      dimension aM(nat),Z(nat),Epot(nstat),wf(nstat,2)
      character(2) S(nat)
      common /first/ Nat,istep,nstat,ndamp,nintc
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv
      common /fifth/ t_old,timekill,prog,Etot_jump,Etot_drift
      common /sixth/ killstat,lvprt
!
!     read input for output reduction, if present
!     write_pos for positions
!      this files contains on the first line these
!      number of entries and on the second line the blank
!      separated list of indices to write to dyn.out
!     existence of file "write_vel" indicates that the
!     velocities are to be printed out
!
      nindices=-1
      allocate (write_pos(1:Nat),STAT=istat)
      if (istat /= 0)  STOP "*** Not enough memory ***"
      inquire(file='cutback_pos', exist=ext)
      IF(ext)THEN
        do i=1,Nat
          write_pos(i)=.false.
        end do
        open(unit=24, file='cutback_pos',status='unknown')
        read(24,*) nindices
        IF(nindices.GT.0)THEN
          allocate (indices_towrite(1:nindices),STAT=istat)
          if (istat /= 0)  STOP "*** Not enough memory ***"
          read(24,*)(indices_towrite(i),i=1,nindices)
          do i=1,nindices
            write_pos(indices_towrite(i))=.true.
          end do
          deallocate(indices_towrite)
        END IF
        close(24)
      ELSE
        do i=1,Nat
          write_pos(i)=.true.
        end do
      END IF
!
      write_vel=.false.
      inquire(file='write_vel', exist=ext)
      IF(ext)THEN
        write_vel=.true.
      END IF
!
      if (istep .eq. 0) then
        if ((nxrestart == 0) .and. (lvprt < 2)) then
          write(16,105)
          write(16,102)
          write(16,105)
          write(16,106) kt,kt*dt*timeunit
          write(16,105)

          write(16,101) istep,nstatdyn,t
          write(16,113) ntype,nhopsold,nrejhopsold

          write(16,105)

!           write(16,'(I6)') nat
!           write(16,105)
          write(16,*) ' Initial geometry: '
          do n=1,nat
            if (write_pos(n)) then
            write(16,'(1x,a2,2x,f5.1,4f14.8)') S(n),Z(n),(x(n,i),i=1,3),aM(n)/proton
!             write(16,'(A2,3F14.8)') S(n),(x(n,i)*au2ang,i=1,3)
            endif
          enddo

          write(16,105)
          if (ndamp .eq. 1) then
            write(16,*) " Damped dynamics: velocity set to zero."
          endif
          if(write_vel)then
            write(16,*) ' Initial velocity: '
            do n=1,nat
              if (write_pos(n)) then
                write(16,'(3F14.8)') (v(n,i), i=1,3)
              endif
            enddo
          else
            write(16,*) ' Initial velocity not printed '
          endif
          write(16,105)

          write(16,*) "    Time  ","  Etot    ","     Ekin    ",&
          "   Epot E0,      E1, ..."
          write(16,110) t,Etot,Ekin,(Epot(k),k=1,nstat)
          write(18,112) t,(Epot(k),k=1,nstat),Epot(nstatdyn),Etot
          write(16,111) 0.0_dpr,0.0_dpr

          write(16,105)
          do k=1,nstat
            write(16,'(A21,I2,A,2F25.14)') " Wave function state ",k,":",&
            wf(k,1),wf(k,2)
          enddo

          write(16,105)
          write(16,*) "------------"
          write(16,105)
        endif
        write(13,*) Etot
        write(14,*) Etot
        write(17,*) 0.0_dpr

      else      ! istep /= 0

        write(16,101) istep,nstatdyn,t
        write(16,113) ntype,nhopsold,nrejhopsold
        write(16,105)

!         write(16,'(I6)') nat
!         write(16,105)
        write(16,*) ' New geometry: '
        do n=1,nat
          if (write_pos(n)) then
            write(16,'(1x,a2,2x,f5.1,4f14.8)') S(n),Z(n),(x(n,i),i=1,3),aM(n)/proton
	    do ns=1,nstat
		    write(16,'100("( ",F14.6,SP,F14.6," ) ")') ( a_coeff(ns,n,i),i=1,3 )
	    end do
!           write(16,'(A2,3F14.8)') S(n),(x(n,i)*au2ang,i=1,3)
          endif
        enddo

      endif

  101 format("STEP ",I8,4X,"Molecular dynamics on state ",I2,4X,"TIME = ",&
      F10.2," fs")
  102 format(6X,"*** Molecular Dynamics ***")
  105 format(" ")                                     !write blank line
  106 format("Output printed at each kt= ",I3," steps (",F5.2," fs).")
  110 format("% ",F7.1,2F14.6,100F14.6)                 !This format allow to write only Epot to 100 states.
  111 format("Etot variation = ",F10.6,                                 &
      " au    Worse conservation = ",F10.6," au")
  112 format(F10.2,2F14.6,100F14.6)                      !This format allow to write only Epot to 100 states.
  113 format("Type of dyn. = ",I4,4X,"N. of hoppings = ",I4,4X,&
             "N. of rejec. hoppings = ",I4)
  114 format("Surface hopping occurred: ", I3,"-->",I3)

      endsubroutine
