
! Date: 2020/02/17
! Version: 2.2
! Description:  
!               Phi 4 model implementation. It includes a simple measurement of the order parameter Q and the energy,
!               and it calculates the propability distribution of each Q value (histogram). Energy units have been 
!               chosen so that E0=1.0.

! Restrictions: cubic lattice (z=6).


 program Phi4_3D

  ! Parameters
    integer, parameter :: L=10, sps=2E7, Ndat=1E4     ! Nº of sites in one dimension, Steps per site, Nº of measured Q
    integer, parameter :: therm=5E5                   ! Thermalization steps
    real, parameter :: d=0.55, CE0=1.0                ! Width param, Coupling const/Energy barrier for local wells
    real, parameter :: beta=1.0/3.0                   ! Inverse temperature
  ! Variables
    integer :: N                                      ! Last site's number: N=L*L*L - 1
    real :: E, dE, Qavrg                              ! Inverse temp, Energy, Energy diff, Average Q value
    real, dimension(0:L*L*L-1) :: X                   ! L*L*L array containing spins (HELICAL B.C.)
    real, dimension(sps) :: Qarr                      ! Array containning Q values 
  ! Dummy variables
    integer :: k, kdat

    N = L*L*l - 1
    
    open(unit=11, file="histogram_T3_p1_d055.txt", status="replace", action="write")
    Histogram_loop: do kdat = 0, Ndat
     ! Initialize 'X'
       X = 1.0
       E = 0.0

       Time: do k = 1, (sps + therm)
          call sweep(L, N, beta, d, CE0, X, E)
          Qarr(k) = sum(X) ! Qi*(N+1)
       enddo Time

     ! Average over the last measures
       call average(sps, sps+therm, Qarr, Qavrg)
       write(unit=11, fmt=*) Qavrg/(N+1)
    enddo Histogram_loop
    close(unit=11)


! -------------------------------------------------------------------------------------------------------------- !


 contains

!  subroutine initial_state

  subroutine random(random_result, low, high)
! Returns a real random value in the range [low, high]. From "Fortran 95 Using F".

     real, intent(in) :: low, high
     real, intent(out) :: random_result
     real :: random_real

     call random_number(random_real)
     random_result = (high-low)*random_real + low
  endsubroutine random


  subroutine nn_values(i, X, L, N, nn_val)
! Returns an array containing the x values of the nearest neighbours. It assumes that the 
! lattice is cubic, and thus the number of nearest neighbours (z) is 6. Helical boundary 
! conditions are implemented here. 

     integer, intent(in) :: i, L, N
     real, dimension(0:N), intent(in) :: X 
     real, dimension(6), intent(out) :: nn_val
     integer :: L2
     
     L2 = L*L
     nn_val = (/ X(modulo(i+1, N+1)), X(modulo(i-1, N+1)),&
                &X(modulo(i+L, N+1)), X(modulo(i-L, N+1)),& 
                &X(modulo(i+L2, N+1)), X(modulo(i-L2, N+1)) /)
  endsubroutine nn_values


  subroutine hamiltonian(xi, dx, nn_val, CE0, H)
! Returns the value of the model's Hamiltonian (divided by E0).

     real, intent(in) :: xi, dx, CE0
     real, dimension(6), intent(in) :: nn_val
     real, intent(out) :: H
     integer :: j
     real :: summation

     summation = 0.0
     Nearest_neigh: do j = 1, 6
        summation = summation + (nn_val(j)-(xi+dx))**2
     enddo Nearest_neigh
     
     H = ( (xi+dx)**2 - 1)**2 + 0.5*CE0*summation
  endsubroutine
    

  subroutine energy_diff(xi, dx, CE0, rand_int, X, L, N, dE)
! Returns the energy difference 'dE' between xi and xi + dx, which is calculated 
! from Phi4 model's Hamiltonian.  

     real, intent(in) :: xi, dx, CE0
     integer, intent(in) :: rand_int, L, N
     real, dimension(0:N), intent(in out) :: X
     real, intent(out) :: dE
     real, dimension(6) :: nn_val
     real :: Ei, Ef

     call nn_values(rand_int, X, L, N, nn_val)
     call hamiltonian(xi, 0.0, nn_val, CE0, Ei)
     call hamiltonian(xi, dx, nn_val, CE0, Ef)

     dE = Ef - Ei   
  endsubroutine energy_diff
        

  subroutine sweep(L, N, beta, d, CE0, X, dE)
! "Performs one 'sweep' of the lattice, i.e., one Monte Carlo step/spin, using Metropolis algorithm."

     integer, intent(in) :: L, N
     real, intent(in) :: beta, d, CE0
     real, dimension(0:N), intent(in out) :: X
     real, intent(out) :: dE
     integer :: i, rand_int
     real :: rand, xi, dx, r
     real, dimension(6) :: nn_val

     Lattice_sweep: do i = 0, N
      ! Generate random site i
        call random(rand, 0.0, real(N))
        rand_int = nint(rand)
        xi = X(rand_int)

      ! Random dx 
        call random(dx, -d, d)

      ! Decide wether to change value from xi to xi + dx: dE positive or negative?
        call energy_diff(xi, dx, CE0, rand_int, X, L, N, dE)
        call random_number(r)

        if (dE <= 0.0) then  
           X(rand_int) = xi + dx

        elseif (r < exp(-dE*beta) ) then
           X(rand_int) = xi + dx

        else ! There is no change in xi, returned dE = 0.0
           dE = 0.0 

        endif
     enddo Lattice_sweep
  endsubroutine sweep


  subroutine average(last, sps, arr, avrg)
! It calculates the average value among the 'last' elements of an array.
  
     integer, intent(in) :: last, sps
     real, dimension(sps), intent(in) :: arr
     real, intent(out) :: avrg
     integer :: h

     avrg = 0 
     Summation: do h = 0, last-1
        avrg = avrg + arr(sps-h)
     enddo Summation
     avrg = avrg/last 

  endsubroutine average


  subroutine histogram(arr, length, round_dec)
! Given any array and its length, it returns its histogram in a .txt file. The array contains real
! values and there's the possibility to round them to any decimal place (using 'round_dec').

     integer, intent(in) :: length, round_dec
     real, dimension(0:length) :: arr
     integer :: i, j, c
     real :: r, numb, rep=9999999.0 ! 'rep' detects the values that are already written in the histogram, we use
                                    ! this value because we know that Q values will be little

     r = 10**round_dec
     arr = nint(arr*r)/r

     open(unit=11, file="histogram.txt", status="replace", action="write")
     Numbers: do i = 0, length
        numb = arr(i)
        c = 1

        if ( numb==rep ) then ! If this happens the value is already in the histogram
           cycle
        endif

        Counter: do j = i+1, length
        if ( arr(j)==numb ) then
           c = c+1
           arr(j) = rep
        endif
        enddo Counter

        write(unit=11, fmt=*) numb, c
     enddo Numbers
     close(unit=11)

  endsubroutine histogram

 
endprogram Phi4_3D 
