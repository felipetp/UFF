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
        
        integer, parameter :: n = 19 ! Valor a definir
        real, parameter :: PI = 3.14159265359
        real, parameter :: g = 9.80665  !m/s²
        real, parameter :: D = 10  !m A definir
        
        !!!!!!!!!!!
        !VÁRIAVEIS!
        !!!!!!!!!!!
        
        real :: p0, pL, L, zeta, f, D, ro, c, Q, S, a, b
        real :: beta, alfa, dt, dx, 
        integer :: i,j
        real, dimension (n,n) :: A, B

        alfa = (16*f*Q)/((D**3)*(c**2)*PI)
        beta = (2*g*sin(teta))/(c**2)

        a = dt/alfa*(dx)**2.0
        b = a*beta*dx/2.0
        zeta = (f/D)*((2.0*ro*c*Q/S)**2.0)
        

        !Construção de malha uniforme
        !dx = (b - a) / (N - 1)

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
       
        do i=2, n-1
        x(i) = 0 + (i)*dx
        p(i) = 
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !DIAGONAL DE COEFICIENTES DA MATRIX A!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !preencher matrix com 0
        do i = 1,n
           do j = 1,n
           A(i,j) = 0
           endo do
        end do

        !diagonal principal
        do i = 1,n
          A(i,i) = (-2.0*a - 2.0)
        end do

        !diagonal inferior
        do i = 2,n
          A(i,i-1) = (a - b)
        end do

        !diagonal superior
        do i = 2, n
           A(i-1,i) = (a + b)
        end do
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !DIAGONAL DE COEFICIENTES DA MATRIX B!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !preencher matrix com 0
        do i = 1,n
           do j = 1,n
           B(i,j) = 0
           endo do
        end do

        !diagonal principal
        do i = 1,n
          B(i,i) = (2.0*a - 2.0)
        end do

        !diagonal inferior
        do i = 2,n
          B(i,i-1) = (-a + b)
        end do

        !diagonal superior
        do i = 2, n
           B(i-1,i) = (-a - b)
        end do

        !contruindo os coeficiente de conhecidos u, matrix k

        !computando o valor da equação de pressão para o tempo desejado

        !algoritmo de Thomas

        !contruindo matrix d

        !construindo c'

        !construindo d'

        !calculando

        !condição de contorno de Neuman

        !saída de data

        end program ferramentasmatematicas 
