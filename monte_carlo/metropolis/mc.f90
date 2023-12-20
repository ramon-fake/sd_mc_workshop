program SpinMonteCarlo
    implicit none
    ! Precision
    integer, parameter :: dp = kind(1.0d0)
    ! Parameters
    integer, parameter :: num_atoms = 10
    integer, parameter :: num_steps = 100000 
    real(dp), parameter :: k_bolt_mry = 0.00633
    real(dp), parameter :: pi = 3.1415926535897931_dp
    ! Variables
    real(dp), allocatable :: energies(:), deltae(:), lattice(:,:), spins(:,:), Jij(:,:), Dij(:,:,:)
    real(dp) :: energy, T, dE
    integer :: step, i, j


    allocate(spins(3, num_atoms))
    allocate(lattice(3, num_atoms))
    allocate(Jij(num_atoms, num_atoms)) 
    allocate(Dij(3, num_atoms, num_atoms))
    allocate(energies(num_steps))
    allocate(deltae(num_steps))


    Jij = 0.0d0
    Dij = 0.0d0
    ! Read inputs and initialize
    call read_lattice(lattice)
    call read_jij(Jij)
    call read_dij(Dij)

    ! Initialize spins with random values
    call init_config(spins)

    T = 0.000
    energies = 0.0D0

    ! Monte Carlo loop
    do step = 1, num_steps
        call monte_carlo_step(spins, lattice, Jij, Dij, T, dE)
        call calculate_energy(spins, lattice, Jij, Dij, energy)

        energies(step) = energy
    end do

    ! Write final spin configuration for visualization
    call write_jmol(spins, lattice)

    ! Print energy as a function of MC steps
    call print_energy(energies, num_steps) 

    call calculate_angles(spins, num_atoms)

    ! End of the main program
    deallocate(spins,lattice,Jij,Dij,energies)
contains

    subroutine read_lattice(lattice)
        real(dp), intent(out) :: lattice(3, num_atoms)
        integer :: unit_num, i, atom_index
        unit_num = 10
        open(unit_num, file="lattice.in", status="old")
        do i = 1, num_atoms
            read(unit_num, *) atom_index, lattice(:, i)
        end do
        close(unit_num)
    end subroutine read_lattice

    subroutine read_jij(Jij)
        real(dp), intent(out) :: Jij(num_atoms, num_atoms)
        integer :: unit_num, i, j
        unit_num = 20
        open(unit_num, file="jij.in", status="old")
        do
            read(unit_num, *, end=100) i, j, Jij(i, j)
        end do
    100 close(unit_num)
    end subroutine read_jij

    subroutine read_dij(Dij)
        real(dp), intent(out) :: Dij(3, num_atoms, num_atoms)
        integer :: unit_num, i, j
        unit_num = 30
        open(unit_num, file="dij.in", status="old")
        do
            read(unit_num, *, end=100) i, j, Dij(:, i, j)
        end do
    100 close(unit_num)
    end subroutine read_dij

    subroutine random_spin(spin)
        real(dp), intent(out) :: spin(3)
        real(dp) :: theta, phi
        theta = 2 * pi * rand()
        phi = acos(2 * rand() - 1.0D0)
        spin = [sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi)]
    end subroutine random_spin

    subroutine init_config(spins)
        real(dp), intent(out) :: spins(3, num_atoms)
        integer :: i
        do i = 1, num_atoms
            call random_spin(spins(:, i))
        end do
    end subroutine init_config

    function energy_difference(i, spins, original_spin, lattice, Jij, Dij) result(dE)
        integer, intent(in) :: i
        real(dp), intent(in) :: spins(3, num_atoms), original_spin(3, num_atoms), lattice(3, num_atoms)
        real(dp), intent(in) :: Jij(num_atoms, num_atoms), Dij(3, num_atoms, num_atoms)
        real(dp) :: dE, E_new, E_old
        integer :: j
        real(dp) :: r_ij(3), norm_r_ij, dot_prod_new, dot_prod_old
    
        E_new = 0.0D0
        E_old = 0.0D0
    
        call calculate_energy(spins, lattice, Jij, Dij, E_new)
        call calculate_energy(original_spin, lattice, Jij, Dij, E_old)

        dE = E_new - E_old
    end function energy_difference

    subroutine calculate_energy(spins, lattice, Jij, Dij, energy)
        real(dp), intent(in) :: spins(3, num_atoms), lattice(3, num_atoms)
        real(dp), intent(in) :: Jij(num_atoms, num_atoms), Dij(3, num_atoms, num_atoms)
        real(dp), intent(out) :: energy
        integer :: i, j
        real(dp) :: r_ij(3), norm_r_ij, dot_prod
    
        energy = 0.0D0
    
        do i = 1, num_atoms  
            do j = 1, num_atoms
              if(i/=j)then 
                dot_prod = dot_product(spins(:, i), spins(:, j))
                energy = energy - Jij(i, j) * dot_prod - dot_product(Dij(:, i, j), cross(spins(:, i), spins(:, j)))
              end if
            end do
        end do
    end subroutine calculate_energy


    function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross

    subroutine monte_carlo_step(spins, lattice, Jij, Dij, T, dE)
        real(dp), intent(inout) :: spins(3, num_atoms)
        real(dp), intent(in) :: lattice(3, num_atoms)
        real(dp), intent(in) :: Jij(num_atoms, num_atoms), Dij(3, num_atoms, num_atoms)
        real(dp), intent(in) :: T
        real(dp), intent(out) :: dE
        integer :: i
        real(dp) :: original_spin(3,num_atoms)
        call random_seed()
        i = int(rand() * num_atoms) + 1
        original_spin(:,:) = spins(:,:)
        call random_spin(spins(:, i))
        dE = energy_difference(i, spins, original_spin, lattice, Jij, Dij)
        write(10,*) exp(-dE / (k_bolt_mry*T)), dE, k_bolt_mry*T,  -(dE / k_bolt_mry*T)
        if (dE > 0.0D0 .and. rand() >= exp(-dE / (k_bolt_mry*T))) then
            spins(:, :) = original_spin(:, :)
        endif
    end subroutine monte_carlo_step

    subroutine write_jmol(spins_in, lattice_in)
        real(dp), intent(in) :: spins_in(3, num_atoms), lattice_in(3, num_atoms)
        integer :: unit_num, i
        unit_num = 40
        open(unit_num, file="jmol-mc.xyz", status="unknown")

        write(unit_num, '(I10)') num_atoms
        write(unit_num, '(A)') "Final Spin Configuration"
        do i = 1, num_atoms
            write(unit_num, '(A, 3F12.6, A, 3F12.6)') "Fe", lattice_in(:, i), " ", spins_in(:, i)
        end do
        close(unit_num)
    end subroutine write_jmol

    subroutine print_energy(energies, num_steps)
        implicit none
        integer, intent(in) :: num_steps
        real(dp), intent(in) :: energies(num_steps)
        integer :: i
    
        ! Open a file for writing
        open(unit=15, file='energy_vs_step-mc.dat', status='replace')
    
        ! Print headers
        write(15, '(A20, A20)') "# MC_Step", "Energy"
    
        ! Print energies for each MC step
        do i = 1, num_steps
            write(15, '(I10, 2X, E15.6)') i, energies(i)
        end do
    
        ! Close the file
        close(15)
    
    end subroutine print_energy

    subroutine calculate_angles(spins, num_atoms)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        implicit none
        real(dp), intent(in) :: spins(3, num_atoms)
        integer, intent(in) :: num_atoms
        integer :: i, j, unit_num
        real(dp) :: angle, angle_degrees, dot_prod
        real(dp), parameter :: pi = 3.14159265358979323846D0
    
        open(unit_num, file="angles-mc.out", status="unknown")
    
        do i = 1, num_atoms - 1
            do j = i + 1, num_atoms
                dot_prod = dot_product(spins(:, i), spins(:, j))
                dot_prod = max(-1.0D0, min(1.0D0, dot_prod))  ! Clamp values to [-1, 1] to avoid numerical issues
                angle = acos(dot_prod)
                angle_degrees = angle * 180.0D0 / pi
                write(unit_num, '(I4, I4, F12.6)') i, j, angle_degrees
            end do
        end do

        close(unit_num)

    end subroutine calculate_angles

end program SpinMonteCarlo
