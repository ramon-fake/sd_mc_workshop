program poormans_spin_dynamics

    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: gamma = 1.76e11_dp                  ! gyromagnetic ratio in rad/(s·T)
    real(dp), parameter :: dt = 1.0e-15
    real(dp), parameter :: pi = 3.1415926535897931_dp
    real(dp) :: alpha = 1.00_dp
    integer, parameter :: num_atoms = 10
    integer, parameter :: num_steps = 1000000 
    real(dp), allocatable :: spins(:,:), lattice(:,:), H_eff(:,:), Jij(:,:), temp_spin(:), energies(:)
    real(dp), allocatable, dimension(:, :, :) :: Dij
    real(dp) :: energy
    integer :: step, i, j

    allocate(spins(3, num_atoms))
    allocate(lattice(3, num_atoms))
    allocate(H_eff(3, num_atoms))
    allocate(temp_spin(3))
    allocate(Jij(num_atoms, num_atoms))
    allocate(Dij(3, num_atoms, num_atoms))
    allocate(energies(num_steps))

    Jij = 0.0d0
    Dij = 0.0d0
    ! Initialize lattice, Jij, Dij
    call read_lattice(lattice)
    call read_jij(Jij)
    call read_dij(Dij)

    ! Initialize spins with random values
    call init_spins(spins, num_atoms)

    ! Dynamics loop using LLG equation
    do step = 1, num_steps
        call compute_Heff(spins, lattice, Jij, Dij, H_eff, num_atoms)
        call llg_update(spins, H_eff, dt, gamma, alpha, num_atoms)
        call write_spin_moments(spins, step)
        call calculate_energy(spins, lattice, Jij, Dij, energy)

        energies(step) = energy   
    end do

    call write_jmol(spins, lattice)

    ! Print energy as a function of MC steps
    call print_energy(energies, num_steps)

    call calculate_angles(spins, num_atoms)


    ! End of main program
    deallocate(spins, lattice, H_eff, temp_spin, Jij, Dij)

contains

    subroutine init_spins(spins, num_atoms)
        integer, intent(in) :: num_atoms
        real(dp), intent(out) :: spins(3, num_atoms)
        integer :: i

        do i = 1, num_atoms
            spins(:, i) = random_spin()
        end do
    end subroutine init_spins

    function random_spin() result(spin)
        real(dp) :: spin(3), theta, phi
        theta = 2 * dp * pi * rand()
        phi = acos(2 * rand() - 1.0_dp)
        spin = [sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi)]
    end function random_spin

    subroutine compute_heff(spins, lattice, Jij, Dij, H_eff, num_atoms)
        integer, intent(in) :: num_atoms
        real(dp), intent(in) :: spins(3, num_atoms), lattice(3, num_atoms), Jij(num_atoms, num_atoms), Dij(3, num_atoms, num_atoms)
        real(dp), intent(inout) :: H_eff(3, num_atoms)
        integer :: i, j
        real(dp) :: rij(3), rij_norm
    
        H_eff = 0.0_dp
        do i = 1, num_atoms
            do j = 1, num_atoms
                if(i/=j)then
                    H_eff(:, i) = H_eff(:, i) + Jij(i, j) * spins(:, j)  +  cross_product(spins(:,j),Dij(:, i, j))
                end if
            end do
        end do
    end subroutine compute_heff

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
                energy = energy - Jij(i, j) * dot_prod - dot_product(Dij(:, i, j), cross_product(spins(:, i), spins(:, j)))
              end if
            end do
        end do
    end subroutine calculate_energy

    function cross_product(a, b) result(cross)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: cross(3)
        cross = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross_product

    subroutine llg_update(spins, H_eff, dt, gamma, alpha, num_atoms)
        integer :: i
        real(dp), intent(in) :: H_eff(3, num_atoms), dt, gamma, alpha
        integer, intent(in) :: num_atoms
        real(dp), intent(inout) :: spins(3, num_atoms)
        real(dp) :: dm(3), torque_term1(3), torque_term2(3), predicted_spins(3, num_atoms), avg_torque(3)
        real(dp) :: H_eff_updated(3, num_atoms) 

        ! Predictor Step
        do i = 1, num_atoms
            torque_term1 = -gamma * cross_product(spins(:, i), H_eff(:, i))
            torque_term2 = -alpha * gamma * cross_product(spins(:, i), cross_product(spins(:, i), H_eff(:, i)))
            dm = torque_term1 + torque_term2
            predicted_spins(:, i) = spins(:, i) + dt * dm
            predicted_spins(:, i) = predicted_spins(:, i) / sqrt(dot_product(predicted_spins(:, i), predicted_spins(:, i)))  ! Normalize
        end do
    
        call compute_heff(predicted_spins, lattice, Jij, Dij, H_eff_updated, num_atoms)
    
        ! Corrector Step
        do i = 1, num_atoms
            torque_term1 = -gamma * cross_product(predicted_spins(:, i), H_eff_updated(:, i))
            torque_term2 = -alpha * gamma * cross_product(predicted_spins(:, i), cross_product(predicted_spins(:, i)&
                                                                                               &, H_eff_updated(:, i)))
            avg_torque = 0.5_dp * (dm + (torque_term1 + torque_term2))
            spins(:, i) = spins(:, i) + dt * avg_torque
            spins(:, i) = spins(:, i) / sqrt(dot_product(spins(:, i), spins(:, i)))  ! Normalize
        end do
    
    end subroutine llg_update

    subroutine read_lattice(lattice)
        integer :: unit_num, i
        real(dp), intent(out) :: lattice(3, num_atoms)
        integer :: atom_index
        unit_num = 10
        open(unit_num, file="lattice.in", status="old")
        do i = 1, num_atoms
            read(unit_num, *) atom_index, lattice(:, i)
        end do
        close(unit_num)
    end subroutine read_lattice

    subroutine read_jij(Jij)
        integer :: unit_num, i, j, num
        real(dp), intent(out) :: Jij(num_atoms, num_atoms)
        unit_num = 20
        open(unit_num, file="jij.in", status="old")
        do
            read(unit_num, *, end=100) i, j, Jij(i, j)
        end do
    100 close(unit_num)
    end subroutine read_jij

    subroutine read_dij(Dij)
        integer :: unit_num, i, j
        real(dp), intent(out) :: Dij(3, num_atoms, num_atoms)
        unit_num = 30
        open(unit_num, file="dij.in", status="old")
        do
            read(unit_num, *, end=100) i, j, Dij(:, i, j)
        end do
    100 close(unit_num)
    end subroutine read_dij

    subroutine write_spin_moments(spins, timestep)
        implicit none
        integer, intent(in) :: timestep
        real(dp), intent(in) :: spins(3, num_atoms)
        integer :: i
        integer, parameter :: unit_num = 40
    
        ! Open the file for the first timestep, otherwise append
        if (timestep == 1) then
            open(unit_num, file="spin_moments.out", status="unknown")
            ! Write header
            write(unit_num, '(A)') "# Timestep Spin1_x Spin1_y Spin1_z ... SpinN_x SpinN_y SpinN_z"
        else
            open(unit_num, file="spin_moments.out", status="old", position="append")
        end if
    
        write(unit_num, '(I6)', advance='no') timestep
        do i = 1, num_atoms
            write(unit_num, '(3F12.6)', advance='no') spins(:, i)
        end do
        write(unit_num, *)  ! Add a newline after writing all spins for this timestep
    
        close(unit_num)
    end subroutine write_spin_moments

    subroutine write_jmol(spins, lattice)
        implicit none
        real(dp), intent(in) :: spins(3, num_atoms), lattice(3, num_atoms)
        integer :: i
        integer, parameter :: unit_num = 50
    
        open(unit_num, file="final_config.jmol", status="unknown")
        
        write(unit_num,*) num_atoms
        write(unit_num,*) 
        do i = 1, num_atoms
            write(unit_num, '(A4, 3F12.6, 3F12.6)') 'Fe', lattice(:, i), spins(:, i)
        end do
    
        close(unit_num)
    end subroutine write_jmol

    subroutine print_energy(energies, num_steps)
        implicit none
        integer, intent(in) :: num_steps
        real(dp), intent(in) :: energies(num_steps)
        integer :: i
    
        ! Open a file for writing
        open(unit=15, file='energy_vs_step-sd.dat', status='replace')
    
        ! Print headers
        write(15, '(A20, A20)') "# SD_Step", "Energy"
    
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
    
        open(unit_num, file="angles.out", status="unknown")
    
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

end program poormans_spin_dynamics
