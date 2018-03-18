!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     -----------------
!c     Moleular dynamics
!c     -----------------
!c         Part III
!c
!c     This program integrates the Newton's eq. in one time step,
!c     by means of the Velocity Verlet Algoritim (VVA).
!c
!c     This part writes outputs.
!c
!c     Mario Barbatti, 2005.
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
    program moldyn03
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      real(kind=dpr) :: iap
      real(kind=dpr), dimension(:,:), allocatable :: x,v,g,a,wf
      real(kind=dpr), dimension(:)  , allocatable :: aM,Z,Epot
      character(2), dimension(:), allocatable :: S
      common /first/ Nat,istep,nstat,ndamp,nintc
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv
      common /fifth/ Etot_jump,Etot_drift

      open(unit= 7, file='control.d',status='unknown')
      open(unit= 8, file='geom' ,    status='unknown')
      open(unit= 9, file='veloc',    status='unknown')
      open(unit=10, file='grad' ,    status='unknown')
      open(unit=11, file='epot' ,    status='unknown')
      open(unit=12, file='auxp' ,    status='old')
      open(unit=13, file='etot0',    status='unknown')
      open(unit=14, file='last' ,    status='unknown')
      open(unit=15, file='../RESULTS/dyn.out',  status='unknown',position='append')
!       open(unit=16, file='../RESULTS/dyn.mld',  status='unknown',position='append')
      open(unit=17, file='worse',    status='unknown')
      open(unit=18, file='../RESULTS/en.dat' ,  status='unknown',position='append')
      open(unit=21, file='typeofdyn.out',status='unknown')
      open(unit=22, file='nhopsold', status='unknown')
      open(unit=23, file='wfrun',    status='unknown')
      open(unit=24, file='inderr',   status='unknown')

      read(7,*) Nat,istep,nstat,nstatdyn,ndamp,kt,dt,t,tmax,nintc,mem,&
                nxrestart,thres,killstat,timekill,prog,lvprt,Etot_jump,Etot_drift

      !Convert
      dt        =dt/timeunit
      Etot_jump =Etot_jump/au2ev
      Etot_drift=Etot_drift/au2ev

      allocate (x(Nat,3),v(Nat,3),g(Nat,3),a(Nat,3),STAT=istat)
      if (istat /= 0)  STOP "*** Not enough memory ***"
      allocate (aM(Nat),Z(Nat),Epot(nstat),S(Nat),wf(nstat,2),STAT=istat2)
      if (istat2 /= 0) STOP "*** Not enough memory ***"

! ..  Read files
      call readfiles(x,v,g,a,aM,S,Z,Epot,wf)

! ..  Get total energy
      call totalenergy(v,aM,Epot,Ekin,Etot)

! ..  Write general output
! mbr the file auxp is from now on written by moldyn.pl
! mbr 0 for printing >0 for non printing
      backspace(12)
      read(12,*) iap
      if (iap .lt. 1.d-9) then
        call writedyn(x,v,aM,S,Z,Etot,Ekin,Epot,wf)
! mbr     iap=1
! mbr    backspace 12
! mbr   write(12,*) iap
! mbr  else
! mbr    iap=iap+1
! mbr    backspace 12
! mbr    write(12,*) iap
      endif

! ..  Write input for the next step calculation
      call writefiles(x,v,aM,S,Z,Etot,Ekin,Epot)

      deallocate (x,v,g,a)
      deallocate (aM,Z,Epot,S,wf)

  100 format(A2,F8.1,4F14.8)
      end program moldyn03

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
      do n=1,nat
      read(9,*) (v(n,i), i=1,3)
      enddo

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

      end subroutine readfiles

! .............................. TOTAL ENERGY

      subroutine totalenergy(v,aM,Epot,Ekin,Etot)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
      dimension x(nat,3),v(nat,3),g(nat,3),a(nat,3),aM(nat),Epot(nstat)
      common /first/ Nat,istep,nstat,ndamp,nintc
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv

! ...... Kinetic energy

      Ekin=0.0_dpr
      do n=1,nat
        v2=( dsqrt( v(n,1)**2+v(n,2)**2+v(n,3)**2 ) )**2
        TK=0.5_dpr*aM(n)*v2
        Ekin=Ekin+TK
      enddo

! ...... Total energy

      Etot=Epot(nstatdyn)+Ekin

      end subroutine totalenergy

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
      rewind 9
      do n=1,nat
        write(9,*) (v(n,i), i=1,3)
      enddo

      end subroutine writefiles

! ............................... WRITEDYN

      subroutine writedyn(x,v,aM,S,Z,Etot,Ekin,Epot,wf)
      use prec_mod
      use units_mod
      implicit real(kind=dpr) (a-h,o-z)
!
      logical :: ext, write_vel
      integer :: nindices, istat, istat2, i
      integer, dimension(:), allocatable :: indices_towrite
      logical, dimension(:), allocatable :: write_pos
!
      dimension x(nat,3),v(nat,3),g(nat,3),a(nat,3)
      dimension aM(nat),Z(nat),Epot(nstat),wf(nstat,2)
      character(2) S(nat)
      common /first/ Nat,istep,nstat,ndamp,nintc
      common /second/ dt,t,tmax,thres
      common /third/ kt,nstatdyn,mem,nxrestart,nhopsold,nrejhopsold
      common /fourth/ ntype,nstatuse,istpj,kpes,nhops,nrejhops,nstatdynprv
      common /fifth/ Etot_jump,Etot_drift
!
!     read input for output reduction, if present
!     write_pos for positions
!     write_vel for velocites
!      each of these files contains on the first line these
!      number of entries and on the second line the blank
!      separated list of indices to write to dyn.out
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
      write(15,105)
      if (ndamp .eq. 1) then
        write(15,*) " Damped dynamics: velocity set to zero."
      endif
      if(write_vel)then
        write(15,*) ' New velocity: '
        do n=1,nat
          if (write_pos(n)) then
            write(15,'(3F14.8)') (v(n,i), i=1,3)
          endif
        enddo
      else
        write(15,*) ' New velocity not printed '
      endif

      write(15,105)

      write(15,*) "    Time  ","  Etot    ","     Ekin    ",&
      "   Epot E0,      E1, ..."
      write(15,110) t,Etot,Ekin,(Epot(k),k=1,nstat)
      write(18,112) t,(Epot(k),k=1,nstat),Epot(nstatdyn),Etot

      read(13,*) Etot0
      read(14,*) Elast
      read(17,*) Eworse
      Evar =Etot-Etot0
      backspace 14
      write(14,*) Etot
      if (dabs(Evar) .gt. dabs(Eworse) ) then
       Eworse=Evar
       backspace 17
       write(17,*) Eworse
      endif
      write(15,111) Etot-Elast,Eworse

      write(15,105)
      do k=1,nstat
      write(15,'(A21,I2,A,2F25.14)') " Wave function state ",k,":",&
      wf(k,1),wf(k,2)
      enddo
      write(15,105)

      write(15,*) "------------"
      write(15,105)

      ! Total energy conservation check ..
      icode = 0
      if (dabs(Etot-Elast) > Etot_jump) then
         write(15,113) Etot_jump*au2ev,dabs(Etot-Elast)*au2ev
         icode = 666
         write(24,*) icode
         STOP 666
      endif
      if (dabs(Eworse) > Etot_drift) then
         write(15,114) Etot_drift*au2ev,dabs(Eworse)*au2ev
         icode = 667
         write(24,*) icode
         STOP 667
      endif

  105 format(" ")                                     !write blank line
  110 format("% ",F10.2,2F14.6,100F14.6)               !This format allow too write only Epot to 100 states.
  111 format("Etot variation = ",F10.6,                                 &
      " au    Worse conservation = ",F10.6," au")
  112 format(F10.2,2F14.6,100F14.6)                    !This format allow too write only Epot to 100 states.
  113 format(/,&
             "************************************************************",&
           /,"ERROR TERMINATION OF MOLDYN03",/,&
             "Total energy changed more than ",F7.3," eV in 1 timestep.",/,&
             "The variation was ",F7.3," eV.",/,&
             "************************************************************",&
           /)
  114 format(/,&
             "************************************************************",&
           /,"ERROR TERMINATION OF MOLDYN03",/,&
             "Total energy changed more than ",F7.3," eV since the 1st step.",/,&
             "The variation was ",F7.3," eV.",/,&
             "************************************************************",&
           /)
      end subroutine writedyn
