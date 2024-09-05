
// Trapezoid Kuralý Ýntegral Hesabý
#include <iostream>
#include <mpi.h>
#include <cmath>
//Süre Hesaplamak Ýçin Chrono Kütüphanesi
#include <chrono>
#include <iomanip>

double f(double x) {
    //Polinom Fonksiyonu
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
    int num_steps = 34747222111; // Toplam iþlem sayýsý
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
    // Slave Corelar integralleri hesaplar
    local_sum = trapezoid_rule(process_start, process_end, steps_per_process);
    // Slave Corelar integrallerin sonuçlarýný toplar 
    MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    auto end_time = std::chrono::steady_clock::now();
    // Master Core integralýn sonucunu yazdýrýr
    if (myid == 0) {
        std::cout << "Toplam alan: " << std::fixed << std::setprecision(10) << total_sum << std::endl;

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Toplam sure: " << duration.count() / 1000.0 << " saniye" << std::endl;
    }

    MPI_Finalize();
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////



 //MONTE CARLO YÖNTEMÝ PÝ SAYISI
#include <iostream>
#include <mpi.h>
#include <random>
#include <chrono> // Çalýþma Süresini Hesaplamak için
// MONTE CARLO YÖNTEMÝ !!!
double calculate_pi(int num_points) {
    int num_inside_circle = 0;

   //Rastgele Sayýlar Üretmek Ýçin
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

   //Oluþturulan noktalar Birim Çemberin içinde mi diye kontrol etme.
    for (int i = 0; i < num_points; ++i) {
        double x = dis(gen);
        double y = dis(gen);
        if (x * x + y * y <= 1.0) {
            num_inside_circle++;
        }
    }

    //Pi Sayýsýný Hesapla
    return 4.0 * static_cast<double>(num_inside_circle) / static_cast<double>(num_points);
}

int main() {
    int myid, numprocs;
    int num_points = 1000000; // Hesaplamanýn hassaslýðý için toplam nokta sayýsý
    double pi = 0.0;
    double total_pi = 0.0;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    
    int points_per_process = num_points / numprocs;

    
    auto start_time = std::chrono::steady_clock::now();

    
    pi = calculate_pi(points_per_process);

    
    MPI_Reduce(&pi, &total_pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   
    auto end_time = std::chrono::steady_clock::now();

    // Master Core Pi Sayýsýný Yazdýrýr
    if (myid == 0) {
        std::cout << "Pi sayisi: " << total_pi / numprocs << std::endl;

        // Geçen süreyi hesapla ve yazdýr
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Toplam sure: " << duration.count() / 1000.0 << " saniye" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
