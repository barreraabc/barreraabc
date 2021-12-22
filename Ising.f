	program ising_exp_crit

c	el objetivo de este programa es simular espines

	implicit none
	real*8 r,Ti,Tf,Term,temp,Tc,nu
	integer*8 L,i,j,steps,lattice,frames,energia,nT,nskip,cero
	integer*8 nescr,Jota,nskap,m,E,iarr,iaba,jder,jizq,h
	character(2) mide
	character(4) anumer

	parameter (L=100)

	dimension lattice(1:L,1:L), r(1:L,1:L)

c---------------------------- Introduciendo datos iniciales
	call entrada(Tc,nu,Term,steps,frames,nT,nescr,Jota)
	Ti=Tc*(1.d0-4.d0*L**(-1/nu))
	Tf=Tc*(1.d0+6.d0*L**(-1/nu))
	Term=Ti
	nskip=steps/nescr
	nskap=steps/frames

	


c---------------------------- Creamos el cuadro de espines
	CALL init_random_seed()

	call random_number(r)

	do i=1,L
	  do j=1,L
!	    lattice(i,j)=nint(r(i,j))*2-1
	    lattice(i,j)=1
	  end do
	end do

  	m=0
	do i=1,L
	  do j=1,L
		m=m+lattice(i,j)
	  end do
	end do
	m=m
	print*, 'Magnetizacion inical=', m

	E=0
	do i=1,L
	  do j=1,L

	  if (i.EQ.1) then
	    iarr=L
	  else 
	    iarr=i-1
	  end if

	  if (i.EQ.L) then
	    iaba=1
	  else 
	    iaba=i+1
	  end if

	  if (j.EQ.1) then
	    jizq=L
	  else
	    jizq=j-1
	  end if

	  if (j.EQ.L) then
	    jder=1
	  else
	    jder=j+1
	  end if

	  energia=-jota*(lattice(iarr,j)+lattice(i,jder)
     &+lattice(i,jizq)+lattice(iaba,j))*lattice(i,j)
	 
	  E=E+energia
	  end do
	end do
	E=E/2 !Dividimos la energia entre dos porque al calcularla
!se tiene en cuenta 2 veces la interaccion de cada espin con su vecino
	print*, 'Energia inicial=', E

c---------------------------- Termalizacion
	cero=0
!	write(anumer,'(i4.4)')cero
!	open(15,file='lattice'//anumer)

	mide='no'
	call move(steps/10,nskip,nskap,lattice,L,Ti,mide,
     &m,E,nescr,h)

	close(15)

c---------------------------- Barrido de temperaturas
	mide='si'
	temp=0.d0

	open(70,file='medias')
	do i=1,nT
!	  write(anumer,'(i4.4)')i
!	  open(69,file='EyM'//anumer)
!	  open(15,file='lattice'//anumer)
	  temp=Ti+(i-1)*(Tf-Ti)/(nT-1)
	  h=i
	  call move(steps,nskip,nskap,lattice,L,temp,mide,
     &m,E,nescr,h)
!	  close(15)
!	  close(69)

	end do
	close(70)

	stop
	end


c---------------------------- Leyendo datos iniciales

	subroutine entrada(Tc,nu,Term,steps,frames,nT,nescr,Jota)
	implicit none
	integer*8 steps,frames,nT,nescr,Jota
	real*8 Tc,nu,Term

	open(22,file='entrada.txt')
	read(22,*)Tc
	write(*,'(a24,f8.3)')'T inicial = ',Tc
	read(22,*)nu
	write(*,'(a24,f8.3)')'T final = ',nu
	read(22,*)Term
	write(*,'(a24,f8.3)')'T de termalizacion = ',Term
	read(22,*)steps
	write(*,'(a24,i12)')'Iteraciones = ',steps
	read(22,*)frames
	write(*,'(a24,i6)')'frames = ',frames
	read(22,*)nescr
	write(*,'(a24,i12)')'medidas = ',nescr
	read(22,*)nT
	write(*,'(a32,i6)')'nº de temperaturas barridas = ',nT
	read(22,*)Jota
	write(*,'(a32,i6)')'J = ',Jota
	close(22)

	return
	end




c-----------Subrutina

	subroutine move(steps,nskip,nskap,lattice,L,temp,mide,
     &m,E,nescr,h)
	implicit none
	real*8 temp,c,r,prob,Enor,Emed,mmed,mnor,mcuad,mcuadmed,Ecuad
	real*8 Ecuadmed,d, mabs, mabsmed,cesp,ji
	integer*8 i,j,k,L,Jota,steps,lattice,z,n,nescr,h
	integer*8 frames,nskip,energia,E,m,nskap,iarr,iaba,jder,jizq
	character(2) mide
!	character(len = 100) fmt
	dimension prob(-8:8)
	dimension lattice(1:L,1:L),c(1:L,1:L),r(1:L,1:L),d(1:L,1:L)
	parameter (Jota=1)

!	write(fmt,'(a, I4, 2x, a )') '(', L, ')'

	print*, 'T(',h,')=', temp

	mnor=0.d0
	mcuad=0.d0
	Enor=0.d0
	Ecuad=0.d0
	mabs=0.d0

	do i=-2,2
	    prob(4*i)=exp(-(4*i)/temp)
	end do

	do z=1,steps
	call random_number(r)
	call random_number(c)
	call random_number(d)
	do k=1,L
	  do n=1,L
	  i=int(c(k,n)*L)+1
	  j=int(d(k,n)*L)+1

	  if (i.EQ.1) then
	    iarr=L
	  else 
	    iarr=i-1
	  end if

	  if (i.EQ.L) then
	    iaba=1
	  else 
	    iaba=i+1
	  end if

	  if (j.EQ.1) then
	    jizq=L
	  else
	    jizq=j-1
	  end if

	  if (j.EQ.L) then
	    jder=1
	  else
	    jder=j+1
	  end if

	   energia=2*jota*(lattice(iarr,j)+lattice(i,jder)
     &+lattice(i,jizq)+lattice(iaba,j))*lattice(i,j)

	  if (prob(energia).GE.r(k,n)) then
	    lattice(i,j)=lattice(i,j)*(-1)
	    m=m+2*(lattice(i,j))
	    E=E+energia
	  end if

	  end do
	end do

!Esto es por si se quiere ver un gif de los espines
!	if (mide.eq.'si') then
!	  if (mod(z,nskap).eq.0)then
!	    write(15,300) lattice
!         end if
!	end if

	if (mide.eq.'si') then
	  if (mod(z,nskip).eq.0)then
	    mnor=mnor+real(m)
	    mabs=mabs+abs(real(m))
	    mcuad=mcuad+(real(m))**2
	    Enor=Enor+real(E)
 	    Ecuad=Ecuad+(real(E))**2
!	    write(69,*) j,m,mnor,mabs,E,Enor
          end if
	end if

	end do

	mabsmed=mabs/((L*L)*nescr) !valor absoluto de m por espin
	mmed=mnor/((L*L)*nescr) !m media por espin
	mcuadmed=mcuad/(nescr*(L*L)**2) !m al cuadrado media por espin
	Emed=Enor/(nescr) !Energia media total
	Ecuadmed=Ecuad/(nescr) !Energia al cuadrado media total
	cesp=(1/temp**2)*(Ecuadmed-Emed**2)/(L*L) !calor especifico
	ji=(L*L/temp)*(mcuadmed-mabsmed**2) !susceptibilidad

	if (mide.eq.'si') then
	write(70,*) temp,mabsmed,Emed/(L*L), cesp,ji,mcuadmed
	end if

!300	format(25(I2, 2x))

!time ./a.out

	return
	end


c-----subrutina para inicializar el generador de numeros aleatorios

      subroutine init_random_seed()
      implicit none
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid, t(2), s
      integer(8) :: count, tms

      call random_seed(size = n)
      allocate(seed(n))
! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream",
     +  form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
        read(un) seed
        close(un)
      else
! Fallback to XOR:ing the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.
        call system_clock(count)
        if (count /= 0) then
          t = transfer(count, t)
        else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000
     -         + dt(2) * 31_8 * 24 * 60 * 60 * 1000
     -         + dt(3) * 24 * 60 * 60 * 60 * 1000
     -         + dt(5) * 60 * 60 * 1000
     -         + dt(6) * 60 * 1000 + dt(7) * 1000
     -         + dt(8)
          t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))
        pid = getpid() + 1099279 ! Add a prime
        s = ieor(s, pid)
        if (n.ge.3) then
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          if (n > 3) then
            seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          end if
        else
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
      end if
      call random_seed(put=seed)
      end subroutine init_random_seed
