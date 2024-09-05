// Trapezoid Kuralý Ýntegral Hesabý Gather ve Scatter fonksiyonlarýna gerek yok çünkü 
// her iþlemci ayný miktarda hesaplama yapýyor
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <chrono>
#include <iomanip>

double f(double x) {
    return x * x;
}

double trapezoid_rule(double start, double end, int num_steps) {
    double sum = 0.0;
    double step_size = (end - start) / num_steps;
    double x = start;

    for (int i = 0; i < num_steps; ++i) {
        sum += (f(x) + f(x + step_size)) * step_size / 2.0;
        x += step_size;
    }

    return sum;
}

int main() {
    int myid, numprocs;
    int num_steps = 575000000;
    double start = 0.0;
    double end = 1.0;
    double local_sum = 0.0;
    double total_sum = 0.0;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    int steps_per_process = num_steps / numprocs;
    double process_start = myid * steps_per_process * (end - start) / num_steps;
    double process_end = (myid + 1) * steps_per_process * (end - start) / num_steps;

    auto start_time = std::chrono::steady_clock::now();

    // Daðýtýlan aralýklarý tüm iþlemcilerle paylaþmak için MPI_Bcast kullanýlýyor.
    MPI_Bcast(&process_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&process_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Tüm iþlemcilerin ayný anda integralleri hesaplamasý için MPI_Allreduce kullanýlýyor.
    local_sum = trapezoid_rule(process_start, process_end, steps_per_process);
    MPI_Allreduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    auto end_time = std::chrono::steady_clock::now();

    if (myid == 0) {
        std::cout << "Toplam alan: " << std::fixed << std::setprecision(10) << total_sum << std::endl;

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Toplam sure: " << duration.count() / 1000.0 << " saniye" << std::endl;
    }

    MPI_Finalize();
    return 0;
}



////////////////////////////////////////////////////////////////////////////////////////////



 //MONTE CARLO YÖNTEMÝ PÝ SAYISI Burada da ayný miktarda hesaplama yapýldýðý için sadece 
 // Bcast ve Reduce Fonksiyonlarý kullanýlmýþtýr

#include <iostream>
#include <mpi.h>
#include <random>
#include <chrono>

double calculate_pi(int num_points) {
    int num_inside_circle = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    for (int i = 0; i < num_points; ++i) {
        double x = dis(gen);
        double y = dis(gen);
        if (x * x + y * y <= 1.0) {
            num_inside_circle++;
        }
    }

    return 4.0 * static_cast<double>(num_inside_circle) / static_cast<double>(num_points);
}

int main() {
    int myid, numprocs;
    int num_points = 100000000;
    double pi = 0.0;
    double total_pi = 0.0;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    int points_per_process = num_points / numprocs;

    auto start_time = std::chrono::steady_clock::now();

    // Her iþlemciye ayrý nokta sayýsýný göndermek için MPI_Bcast kullanýlýyor.
    MPI_Bcast(&points_per_process, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Her iþlemci kendi nokta sayýsýný kullanarak Pi'yi hesaplýyor.
    pi = calculate_pi(points_per_process);

    // Tüm iþlemcilerin hesapladýðý Pi deðerlerini toplamak için MPI_Reduce kullanýlýyor.
    MPI_Reduce(&pi, &total_pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    auto end_time = std::chrono::steady_clock::now();

    if (myid == 0) {
        std::cout << "Pi sayisi: " << total_pi / numprocs << std::endl;

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Toplam sure: " << duration.count() / 1000.0 << " saniye" << std::endl;
    }

    MPI_Finalize();
    return 0;
}



/////////////////////////////////////////////////////////
//Bariyer Ekleme
#include <iostream>
#include <mpi.h>

int main() {
    MPI_Init(NULL, NULL);

    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    std::cout << "Hello world from process " << my_rank << "/" << world_size << std::endl;

     //Bariyer fonksiyonu ekleniyor
    MPI_Barrier(MPI_COMM_WORLD);

     //Tüm iþlemciler bariyere ulaþtýðýnda bu kýsýma geçilir
    if (my_rank == 0) {
        std::cout << "All processes reached the barrier." << std::endl;
    }

    MPI_Finalize();
    return 0;
}

