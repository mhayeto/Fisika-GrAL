
! Date: 2020/02/27
! Version: 3.2
! Description: 
!               Ising  3D model (B=0), from "Monte Carlo Methods in Statistical Physics" book (pages
!               47 and 433), with helical boundary conditions (page 334). It includes simple 
!               measurements of magnetization and internal energy, the calculation of the 
!               autocorrelation functions (page 59) and the integrated autocorrelation time (page 62), 
!               and measurements with varying temperature (pages 20 and 76).
!
! Restrictions: cubic lattice (z=6), B=0, single-spin-flip. 


 program Ising2D

  ! Parameters
    integer, parameter :: L=100, sps=10000            ! Number of sites in one dimension, Steps per site, Initial temp
    integer, parameter :: tmax=600                    ! Number of time loops to calculate the autocorr. function
    real, parameter :: J=1.0, B=0.0, Tf=5.0, dT=0.2   ! Exchange energy, Inverse temperature, Magnetic field, Final temp, Temp step     
  ! Variables
    integer :: N, Nt                                  ! Last site's number: N=L*L - 1, Number of temp. loops
    integer, dimension(0:L*L*L-1) :: S                ! L*L*L array containing spins (1D because of helical b.c.)
    real :: beta, E, deltaE                           ! Inverse temp, Energy of the system, Energy diff. between two states
    real :: Eavrg, Mavrg, E2avrg, M2avrg              ! Average energy and magnetization values, and their squares
    real :: tau                                       ! Correlation time
    real, dimension(3) :: prob                        ! Array containing possible acceptance ratios
    real, dimension(0:sps) :: Earr, Marr              ! Arrays containing energy and magnetization values
  ! Dummy variables
    integer :: k, tk

    N = L*L*L - 1 
    Nt = Tf/dT
    open(unit=11, file="simulation.txt", status="replace", action="write")
    open(unit=22, file="sim_autocorr.txt", status="replace", action="write")

  ! Initialize 'S', 'prob', 'E' and 'M'
    call initial_state(S, N, 0)
    call energy_from_state(S, L, N, J, B, E)
    write(unit=11, fmt=*) 0, E/(N+1), real(sum(S))/(N+1)

    Temp_sweep: do tk = 1, Nt
       call initial_state(S, N, 0)
       call energy_from_state(S, L, N, J, B, E)
       beta = 1.0/(dT*tk) 
       call acceptance_ratios(beta, J, prob)

       Earr(0) = E
       Marr(0) = real(sum(S))

       Time: do k = 1, sps
          call sweep(L, N, beta, J, prob, S, deltaE)
          E = E + deltaE
          Earr(k) = E
          Marr(k) = real(sum(S))
       enddo Time
       
     ! Average over the last 7000 measures  
       call average(7000, sps, Earr, Eavrg)
       call average(7000, sps, Marr, Mavrg)
       call average(7000, sps, Earr**2, E2avrg)
       call average(7000, sps, Marr**2, M2avrg)

       c = beta*beta*(E2avrg - Eavrg**2)/(N+1)
       susc = beta*(M2avrg - Mavrg**2)/(N+1)
       write(unit=11, fmt=*) dT*tk, Eavrg/(N+1), Mavrg/(N+1), c, susc

       ! Magnetization's autoccorrelation function and correlation time (tau)
       call autocorrelation(tmax, sps, Marr, tau)
       write(unit=22, fmt=*) dT*tk, tau

    enddo Temp_sweep
    close(unit=11)
    close(unit=22)


! -------------------------------------------------------------------------------------------------------------- !

 contains

  subroutine initial_state(S, N, state)
! Returns the array S corresponding to the initial state. The variable 'state' corresponds to its 
! temperature: if state=0 it's the T=0 state, and if state=1 its the T=inf state (random-spin).

     integer, intent(in) :: N, state
     integer, dimension(0:N), intent(in out) :: S
     integer, dimension(0:1) :: spins = (/-1, 1/)
     integer :: is
     real :: rand_ind

     select case (state)
        case(0)
           S = 1
        case(1)
           Lattice: do is = 0, N
              call random_number(rand_ind)
              S(is) = spins(nint(rand_ind))
           enddo Lattice
     endselect
  endsubroutine initial_state


  subroutine acceptance_ratios(beta, J, prob)
! Returns an array containing possible values for acceptance ratios exp(-\beta\deltaE). It assumes 
! that the lattice is cubic and thus z=6. Then it can be seen that there are just three possible values.

     real, intent(in) :: beta, J
     real, dimension(3), intent(out) :: prob

     prob(1) = exp(-12.0*beta*J)
     prob(2) = exp(-8.0*beta*J)
     prob(3) = exp(-4.0*beta*J)
  endsubroutine acceptance_ratios


  subroutine energy_from_state(S, L, N, J, B, E) 
! Given any S state, it returns the energy of the system, which is calculated by using Ising model's Hamiltonian.
! It can calculate any state's energy, but it will be used only to calculate the energy of the initial state, as 
! "the clever thing is to calculate the energy of the system from the Hamiltonian at the very first of the simulation"
! (page 58). Boundar conditions: helical.

     integer, intent(in) :: L, N
     real, intent(in) :: J, B
     integer, dimension(0:N), intent(in) :: S
     real, intent(out) :: E
     integer :: i, spin_sum

     L2 = L*L
     spin_sum = 0
     Lattice_sweep: do i = 0, N
        spin_sum = spin_sum + S(i)*(S(modulo(i+1, N+1)) + S(modulo(i-1, N+1)) +& 
                                   &S(modulo(i+L, N+1)) + S(modulo(i-L, N+1)) +&
                                   &S(modulo(i+L2, N+1)) + S(modulo(i-L2, N+1)))
     enddo Lattice_sweep
     E = -J*spin_sum/2.0
  endsubroutine energy_from_state


  subroutine random_int(random_result, low, high)
! Returns an integer random value in the range [low, high]. From "Fortran 95 Using F".

     integer, intent(in) :: low, high
     integer, intent(out) :: random_result
     real :: random_real

     call random_number(random_real)
     random_result = int( (high-low+1)*random_real + low )
  endsubroutine random_int


  subroutine nn_up_spins(i, S, L, N, spin_sum)
! Returns the sum of the nearest neighbour spins. It assumes that the lattice is square, and thus
! the number of nearest neighbours (z) is 4. Helical boundary conditions are implemented here. 

     integer, intent(in) :: i, L, N
     integer, dimension(0:N), intent(in) :: S 
     integer, intent(out) :: spin_sum

     L2 = L*L
     spin_sum = S(modulo(i+1, N+1)) + S(modulo(i-1, N+1)) +&
               &S(modulo(i+L, N+1)) + S(modulo(i-L, N+1)) +&
               &S(modulo(i+L2, N+1)) + S(modulo(i-L2, NN+1))
  endsubroutine nn_up_spins
     

  subroutine sweep(L, N, beta, J, prob, S, deltaE)
! "Performs one 'sweep' of the lattice, i.e., one Monte Carlo step/spin, using Metropolis algorithm."
! It returns the energy difference (deltaE) between the previous and the next states.  

     integer, intent(in) :: L, N
     real, intent(in) :: beta, J
     real, intent(out) :: deltaE
     integer, dimension(0:N), intent(in out) :: S
     real, dimension(3), intent(in) :: prob
     integer :: rand_int, sk, spin_sum
     real :: r

      deltaE = 0.0
      Lattice_sweep: do i = 0, N
      ! Generate random site k
        call random_int(rand_int, 0, N)
        sk = S(rand_int)

      ! Decide wether to flip spin
        call nn_up_spins(rand_int, S, L, N, spin_sum)
        select case (sk*spin_sum)

           case(6) ! \deltaE = 12J > 0, then we will "flip it with probability A = exp(-\beta\deltaE)"
              call random_number(r)
              if (r < prob(1)) then
                 S(rand_int) = -sk ! "if r < A, then we flip the spin"
                 deltaE = deltaE + 2.0*J*6.0 ! 6=sk*spin_sum
              endif

           case(4) ! \deltaE = 8J > 0, then we will "flip it with probability A = exp(-\beta\deltaE)"
              call random_number(r) 
              if (r < prob(2)) then
                 S(rand_int) = -sk ! "if r < A, then we flip the spin"
                 deltaE = deltaE + 2.0*J*4.0 ! 4=sk*spin_sum
              endif

           case(2) ! \deltaE = 4J > 0, then we will "flip it with probability A = exp(-\beta\deltaE)"
              call random_number(r)
              if (r < prob(1)) then
                 S(rand_int) = -sk ! "if r < A, then we flip the spin"
                 deltaE = deltaE + 2.0*J*2.0 ! 2=sk*spin_sum
              endif

           case default ! \deltaE <= 0, then A=1 and "we accept the move and flip the spin"
              S(rand_int) = -sk
              deltaE = deltaE + 2.0*J*sk*spin_sum

        endselect
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


! ----------------------------- AUTOCORRELATION ZUZENDU ------------------------------- !
!  Tc baino tenperatura altuagoetarako aceptable funzionatzen du, bajuagoetarako fatal  !

  subroutine autocorrelation(tmax, sps, arr, tau)
! It calculates the autocorrelation function and the integrated autocorrelation function 
! (or correlation time, 'tau') of the data given by array 'arr' (energy or magnetization). 
! It must be used once after subroutine 'sweep' is called). 

     integer, intent(in) :: tmax, sps
     real, dimension(0:sps), intent(in) :: arr
     real, intent(out) :: tau
     integer :: t_prime, t
     real :: summation, integral, integral0
     real, dimension(0:tmax) :: acorr_f

     summation = 0.0
     Integrate_t0: do t_prime = 0, sps
        summation = summation + arr(t_prime)**2
     enddo Integrate_t0
     integral0 = (summation - (summation/sps))/sps
     acorr_f(0) = integral0

     Time: do t = 1, tmax
        summation = 0.0
        Integrate: do t_prime = 0, sps-t
           summation = summation + (arr(t_prime)*arr(t_prime+t))
        enddo Integrate
        integral = (summation - (summation/(sps-t)))/(sps-t)
        acorr_f(t) = integral
     enddo Time
     acorr_f = acorr_f/integral0

     tau = 0.0
     Integration: do t = 0, tmax-1
        tau = tau + (acorr_f(t) + acorr_f(t+1))/2.0
     enddo Integration

  endsubroutine autocorrelation
 
 endprogram Ising2D  

