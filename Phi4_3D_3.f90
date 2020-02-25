
! Date: 2020/01/15
! Version: 2.3
! Description:  
!               Phi 4 model implementation. It includes a simple measurement of the order parameter Q,
!               the energy and the partition function with varying temperature. 

! Restrictions: cubic lattice (z=6).


 program Phi4_3D

  ! Parameters
    integer, parameter :: L=50, sps=2E7, therm=5E5      ! Number of sites in one dimension, Steps per site, Therm. steps
    real, parameter :: d=0.55, CE0=1.0                  ! Width param, Coupling consty/Energy barrier for local wells
    real, parameter :: Tf=5.0, dT=0.05                  ! Final temp, Temp. step
  ! Variables
    integer :: N, Nt                                    ! Last site's number: N=L*L*L - 1, Number of temp loops
    real :: beta, E                                     ! Inverse temp, Energy
    real :: Qavrg                                       ! Average Q value and energy
    real, dimension(0:L*L*L-1) :: X                     ! L*L*L array containing spins (HELICAL B.C.)
    real, dimension(sps) :: Qarr                        ! Array containing Q and energy values
  ! Dummy variables 
    integer :: k, tk

    N = L*L*l - 1 
    Nt = Tf/dT

    open(unit=11, file="temp_L50.txt", status="replace", action="write")
    Temp_sweep: do tk = 1, Nt
     ! Initialize 'X' and energy
       call initial_state(X, L, N, 0, E)
       beta = 1.0/(dT*tk)

       Time: do k = 1, sps
          call sweep(L, N, beta, d, CE0, X, E)
          Qarr(k) = sum(X) ! Qi*(N+1)
       enddo Time

     ! Average over the last 17000 measures
       call average(sps-500, sps, Qarr, Qavrg)
       write(unit=11, fmt=*) dT*tk, Qavrg/(N+1)
    enddo Temp_sweep
    close(unit=11)
    

! -------------------------------------------------------------------------------------------------------------- !

 contains

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


  subroutine initial_state(X, L, N, state, E)
! Returns the array X corresponding to the initial state. The variable 'state' corresponds to its
! temperature: if state=0 it's the X=1 state, and if state=1 it's an initial random state (X=rand).
  
     integer, intent(in) :: L, N, state
     real, dimension(0:N), intent(in out) :: X
     real, intent(out) :: E
     integer :: ix
     real :: r, dE
     real, dimension(6) :: nn_val

     select case(state)
        case(0)
           X = 1.0
           E = 0.0
        case(1)
           Lattice: do ix = 0, N
              call random(r, -1.0, 1.0)
              X(ix) = r
           enddo Lattice
           E = 0.0
           Energy: do ix = 0, N
              call nn_values(ix, X, L, N, nn_val)
              call hamiltonian(X(ix), 0.0, nn_val, CE0, dE)
              E = E + dE
           enddo Energy
     endselect
   endsubroutine initial_state


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


