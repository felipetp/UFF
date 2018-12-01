        program ferramentasmatematicas
        implicit none
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                            !
        !    PEDRO LINS DE MOURA MARTINS DA COSTA    !
        !                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        !!!!!!!!!!!!
        !PARÂMETROS!
        !!!!!!!!!!!!
        
        integer, parameter :: JJ = 20  !distance steps  (Valor a definir)
        integer, parameter :: NN = 20  !time steps (Valor a definir)
        integer, parameter :: dx =  ! Valor a definir
        real, parameter :: Ttot = ! Valor a definir
        real, parameter :: PI = 3.14159265359
        real, parameter :: g = 9.80665  !m/s²
        real, parameter :: D = 0.07  !m A definir
        real, parameter :: f = 0.002  ! pq?
        real, parameter :: ro = 1.184 !kg/m³  (temp. 25º)
        real, parameter :: c = 340.29 !m/s  (temp. 25º)
        real, parameter :: teta = 0 !tubulação horizontal
        real, parameter :: Q = 16 !m3/s

        
        !!!!!!!!!!!
        !VÁRIAVEIS!
        !!!!!!!!!!!
        
        real :: p0, pL, L, zeta, S, a, b, dt
        real :: beta, alfa, dx, 
        integer :: i,j, N, M
        real, dimension :: x(0:JJ),g(0:NN+1,0:JJ),asup(0:JJ),aprin(0:JJ),ainf(0:JJ), bcol(0:JJ),g_novo(0:JJ)
        
        
        dt = Ttot/N
        S = PI*((d/2.0)**2)
        alfa = (16.0*f*Q)/((D**3.0)*(c**2.0)*PI)
        beta = (2.0*g*sin(teta))/(c**2.0)

        a = dt/alfa*(dx)**2.0
        b = a*beta*dx/2.0
        zeta = (f/D)*((2.0*ro*c*Q/S)**2.0)
        
        !!!!!!!!!!!!!!!!!!
        !CONDIÇÃO INICIAL!
        !!!!!!!!!!!!!!!!!!

        do i=0, JJ
           g(0,i) = 1
        end do
        
        !!!!!!!!!!!!!!!!!!!!!!!
        !CONDIÇÕES DE CONTORNO!
        !!!!!!!!!!!!!!!!!!!!!!!

        do n=0, NN
           g(n,0)=7 !pressão constante na saída do compressor
           g(n,JJ)= SQRT((p0)**2.0-(zeta*L)) !assumimos que o último elemento discreto é regime permanente.
        end do

       
        !!!!!!!!!!!!!!!!!!!!!!!!
        !RESOLVENDO AS MATRIZES!
        !!!!!!!!!!!!!!!!!!!!!!!!
        
        do n = 0 , NN !loop
        do j = 1, JJ-1  ! a,b,c are the coefficients of C-N scheme and d is the right part
                a(j) = (a-b)
                b(j) = -2.0-2.0*a
                c(j) = (a+b)
                d(j) = (-a-b)*g(n,j+1)+(2.0-2.0*a)*g(n,j)+(-a+b)*g(n,j-1)  
        end do  
   
        call thomas_algorithm(a,b,c,d,JJ,g_novo)     
        g(n+1,j)= g_novo(j)
        end do        
        
        open(10,file='crank_nicolsan.txt')
        write(10,*)  'ApproximateSolution =[',x(0),u(0,0)
        do j =0, JI
          do n= 0, NI
        write(10,*) x(j),u(n,j)
         end do
         end do
        write(10,*)x(JI),u(NI,JI),']'
         close(10)
        end program ferramentasmatematicas 

        subroutine thomas_algorithm(a,b,c,d,JJ,g_novo)
 
!	 a - sub-diagonal (means it is the diagonal below the main diagonal)
!	 b - the main diagonal
!	 c - sup-diagonal (means it is the diagonal above the main diagonal)
!	 d - right part
!	 u_new - the answer

        implicit none
        integer,intent(in) :: JJ
        real,intent(in) :: a(JJ),c(JJ)
        real,intent(inout),dimension(JJ) :: b,d
        real,intent(out) :: g_novo(JJ)
        integer j
        real ::q
 
        do j = 2,JJ                 !combined decomposition and forward substitution
            q = a(j)/b(j-1)
              b(j) = b(j)-q*c(j-1)
              d(j) = d(j)-q*d(j-1)
           end do 
 
        !back substitution
           q = d(JJ)/b(JJ)
           g_novo(JJ)=q
           do j = JJ-1,1,-1
              q = (d(j)-c(j)*d(j+1))/b(j)
             g_novo(j)=q
           end do
          RETURN
        end subroutine thomas_algorithm



        
