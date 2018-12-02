        program ferramentasmatematicas
        implicit none
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                            !
        !     FELIPE FERNANDES TEPEDINO              !
        !     PEDRO LINS DE MOURA MARTINS DA COSTA   !
        !                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        !!!!!!!!!!!!
        !PARÂMETROS!
        !!!!!!!!!!!!
        
        integer, parameter :: JJ = 100000  !distance steps  (Valor a definir)
        integer, parameter :: NN = 20  !time steps (Valor a definir)
        double precision, parameter :: dx = 0.1 !m
        double precision, parameter :: Ttot = 200! Valor a definir
        double precision, parameter :: PI = 3.14159265359
        double precision, parameter :: grav = 0.6  !m/s²
        double precision, parameter :: D = 0.6  !m A definir
        double precision, parameter :: f = 0.003  ! pq?
        double precision, parameter :: ro = 0.73 !kg/m³  (temp. 25º)
        double precision, parameter :: c = 340.29 !m/s  (temp. 25º)
        double precision, parameter :: teta = 0 !tubulação horizontal
        double precision, parameter :: Q = 10 !m3/s
        double precision, parameter :: p0 = 7000000.0 

        
        !!!!!!!!!!!
        !VÁRIAVEIS!
        !!!!!!!!!!!
        real(8), allocatable :: x(:),g(:,:),asup(:),aprin(:),ainf(:)
        real(8), allocatable :: bcol(:),g_novo(:)
        double precision :: zeta, S, a, b, dt
        double precision :: beta, alfa 
        integer :: i,j, N
        allocate (x(0:JJ), g(0:NN+1,0:JJ),asup(0:JJ),aprin(0:JJ))
        allocate (ainf(0:JJ),bcol(0:JJ),g_novo(0:JJ))

         
        dt = Ttot/NN
        S = PI*((d/2.0)**2)
        alfa = (16.0*f*Q)/((D**3.0)*(c**2.0)*PI)
        beta = (2.0*grav*sin(teta))/(c**2.0)

        a = dt/alfa*(dx)**2.0
        b = a*beta*dx/2.0
        zeta = (f/D)*((2.0*ro*c*Q/S)**2.0)
        
        print *, a, b , zeta
        !!!!!!!!!!!!!!!!!!
        !CONDIÇÃO INICIAL!
        !!!!!!!!!!!!!!!!!!

        do i=0, JJ
           x(i)= 0 + i*dx
           do n = 0,NN
           g(n,i) = 1
           end do     
        end do
        

        !!!!!!!!!!!!!!!!!!!!!!!
        !CONDIÇÕES DE CONTORNO!
        !!!!!!!!!!!!!!!!!!!!!!!

        do n=0, NN
           g(n,0)=7 !pressão constante na saída do compressor
           enddo

           do n=0, NN  
           g(n,JJ)= SQRT((p0)**2.0-(zeta*JJ)) !assumimos que o último elemento discreto é regime permanente.
        end do

       
        !!!!!!!!!!!!!!!!!!!!!!!!
        !RESOLVENDO AS MATRIZES!
        !!!!!!!!!!!!!!!!!!!!!!!!
        
        do n = 0 , NN !loop
        do j = 1, JJ-1  ! a,b,c are the coefficients of C-N scheme and d is the right part
           asup(j) = (-a-b)
           aprin(j) = -2.0+2.0*a
           ainf(j) = (-a+b)
           bcol(j) = (a-b)*g(n,j+1)+(-2.0-2.0*a)*g(n,j)+(a+b)*g(n,j-1)  
        end do  
   
        call  solve_tridiag(ainf,aprin,asup,bcol,g_novo, JJ)     
        g(n+1,j)= g_novo(j)
        
        print *, asup(j), aprin(j), ainf(j), bcol(j)
        end do        
        
        !GRAVANDO OS RESULTADOS
        
        open(10,file='crank_nicolsan.txt',status='REPLACE')
        write(10,*)  'ApproximateSolution =[',x(0),g(0,0)
        do j =0, JJ
           do n= 0, NN
                write(10,*) j, n, x(j),g(n,j)
           end do
        end do
         
        close(10)
        end program ferramentasmatematicas 

        
        
        !SUBROTINA METODO DE THOMAS BOLADO

            
         subroutine solve_tridiag(ainf,aprin,asup,bcol,g_novo,JJ)
      implicit none
!	 a - sub-diagonal (means it is the diagonal below the main diagonal)
!	 b - the main diagonal
!	 c - sup-diagonal (means it is the diagonal above the main diagonal)
!	 d - right part
!	 x - the answer
!	 n - number of equations

        integer,intent(in) :: JJ
        real(8),dimension(JJ),intent(in) :: ainf,aprin,asup,bcol
        real(8),dimension(JJ),intent(out) :: g_novo
        real(8),dimension(JJ) :: cp,dp
        real(8) :: m
        integer i

! initialize c-prime and d-prime
        cp(1) = asup(1)/aprin(1)
        dp(1) = bcol(1)/aprin(1)

        print *, cp, dp
! solve for vectors c-prime and d-prime
         do i = 2,JJ
           m = aprin(i)-cp(i-1)*ainf(i)
           cp(i) = asup(i)/m
           dp(i) = (bcol(i)-dp(i-1)*ainf(i))/m
         enddo
! initialize x
         g_novo(JJ) = dp(JJ)
! solve for x from the vectors c-prime and d-prime
        do i = JJ-1, 1, -1
          g_novo(i) = dp(i)-cp(i)*g_novo(i+1)
        end do

        end subroutine solve_tridiag


               



