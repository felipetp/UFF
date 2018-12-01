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
        
        
        integer, parameter :: N = 20  !time steps Valor a definir
        integer, parameter :: dx =  ! Valor a definir
        real, parameter :: Ttot = ! Valor a definir
        real, parameter :: PI = 3.14159265359
        real, parameter :: g = 9.80665  !m/s²
        real, parameter :: D = 0.07  !m A definir
        real, parameter :: f = 0.002  ! pq?
        real, parameter :: ro = 1.184 !kg/m³  (temp. 25º)
        real, parameter :: c = 340.29 !m/s  (temp. 25º)
        real, parameter :: teta = 0 !tubulação horizontal
        
        
        !!!!!!!!!!!
        !VÁRIAVEIS!
        !!!!!!!!!!!
        
        real :: p0, pL, L, zeta, Q, S, a, b, dt
        real :: beta, alfa, dx, 
        integer :: i,j, N, M
        real, dimension (n,n) :: A, B
        real, dimension :: x(:),g(:,:),asup(:),aprin(:),ainf(:), bcol(:),g_novo(0:)
        
        
        
        dt = Ttot/N
        S = PI*((d/2.0)**2)
        alfa = (16.0*f*Q)/((D**3.0)*(c**2.0)*PI)
        beta = (2.0*g*sin(teta))/(c**2.0)

        a = dt/alfa*(dx)**2.0
        b = a*beta*dx/2.0
        zeta = (f/D)*((2.0*ro*c*Q/S)**2.0)
        
        !!!!!!!!!!!!!!!!!!!!!!!
        !CONDIÇÕES DE CONTORNO!
        !!!!!!!!!!!!!!!!!!!!!!!

        !pressão constante na saída do compressor
        p(0) = 7.0
        
        !assumimos que o último elemento discreto é regime permanente.
        p(L) = SQRT((p0)**2.0-(zeta*L))


        !!!!!!!!!!!!!!!!!!
        !CONDIÇÃO INICIAL!
        !!!!!!!!!!!!!!!!!!

        do i=2, n-1
           p(i) = 1
        end do
       
        !!!!!!!!!!!!!!!!!!!!!!!!
        !RESOLVENDO AS MATRIZES!
        !!!!!!!!!!!!!!!!!!!!!!!!
        
        do n = 0 , N !loop
        do j = 1, JI-1  ! a,b,c are the coefficients of C-N scheme and d is the right part
                a(j) = (a-b)
                b(j) = -2.0-2.0*a
                c(j) = (a+b)
                d(j) = (-a-b)*g(n,j+1)+(2.0-2.0*a)*g(n,j)+(-a+b)*g(n,j-1)  
        end do  
   
        call thomas_algorithm(a,b,c,d,JI,u_new)     
        g(n+1,j)= g_novo(j)
        end do        

        end program ferramentasmatematicas 
