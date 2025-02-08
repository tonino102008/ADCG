#include <iostream>
#include <chrono>
#include <array>

#include "adcg_ascalar.hpp"
#include "adcg_binary_ops.hpp"
#include "adcg_math_ops.hpp"

typedef ascalar<ascalar<real>> scalar;

#define N_STEPS     501
#define N_STATES    7
#define N_CONTROLS  2

int main() {

    double g = 9.81;
    double m = 1360.0 + 2 * 45.85506 + 2 * 51.3019091;
    double Jx = 265;
    double Jy = 1040;
    double Jz = 1690;
    double a = 1.52;
    double b = 1.22;
    double wf = 1.678;
    double wr = 1.691;
    double h = 0.486;
    double R = 0.326;
    double cDrag = 0.460118343195266;
    double Fzf = 0.5 * m * g * b / (a + b);
    double Fzr = 0.5 * m * g * a / (a + b);
    double KRollF = 44915.62446386859;
    double DRollF = 9529.744091087500;
    double KRollR = 89993.18337893847;
    double DRollR = 15159.55639325625;
    double ratio = 0.5;
    double Camax = 585*3.49;
    double scaleBrake = 1.7633;
    double deltamax = 1.047197551196598;

    double Fz0fPAC = 4000;
    double Fz0rPAC = 4000;
    double Kf = 19.7969586572;
    double Kf2 = 1.79989978802;
    double Cf = 1.32229820594; 
    double Df = 1.1;
    double Ef = -0.637732885602;
    double Kr = 19.7969586572;
    double Kr2 = 1.79989978802;
    double Cr = 1.22229820594; 
    double Dr = 1.3;
    double Er = -0.637732885602;

    double nmax = 5;
    double smax = 1.0/2.627990665203855e3;
    double k = 0;
    double ds = 0.002;

    std::array<scalar, N_STEPS * (N_STATES + N_CONTROLS)> vars;
    for (size_t i = 0; i < N_STEPS; i++) {
        vars[i * (N_STATES + N_CONTROLS)] = scalar(94 / 3.6, i * (N_STATES + N_CONTROLS));
        vars[i * (N_STATES + N_CONTROLS) + 1] = scalar(0.0, i * (N_STATES + N_CONTROLS) + 1);
        vars[i * (N_STATES + N_CONTROLS) + 2] = scalar(0.0, i * (N_STATES + N_CONTROLS) + 2);
        vars[i * (N_STATES + N_CONTROLS) + 3] = scalar(0.0, i * (N_STATES + N_CONTROLS) + 3);
        vars[i * (N_STATES + N_CONTROLS) + 4] = scalar(0.0, i * (N_STATES + N_CONTROLS) + 4);
        vars[i * (N_STATES + N_CONTROLS) + 5] = scalar(0.0, i * (N_STATES + N_CONTROLS) + 5);
        vars[i * (N_STATES + N_CONTROLS) + 6] = scalar(0.0, i * (N_STATES + N_CONTROLS) + 6);

        vars[i * (N_STATES + N_CONTROLS) + 7] = scalar(0.0, i * (N_STATES + N_CONTROLS) + 7);
        vars[i * (N_STATES + N_CONTROLS) + 8] = scalar(0.0, i * (N_STATES + N_CONTROLS) + 8);

        vars[i * (N_STATES + N_CONTROLS)].getValue().getDerivative()[i * (N_STATES + N_CONTROLS)] = 1;
        vars[i * (N_STATES + N_CONTROLS) + 1].getValue().getDerivative()[i * (N_STATES + N_CONTROLS) + 1] = 1;
        vars[i * (N_STATES + N_CONTROLS) + 2].getValue().getDerivative()[i * (N_STATES + N_CONTROLS) + 2] = 1;
        vars[i * (N_STATES + N_CONTROLS) + 3].getValue().getDerivative()[i * (N_STATES + N_CONTROLS) + 3] = 1;
        vars[i * (N_STATES + N_CONTROLS) + 4].getValue().getDerivative()[i * (N_STATES + N_CONTROLS) + 4] = 1;
        vars[i * (N_STATES + N_CONTROLS) + 5].getValue().getDerivative()[i * (N_STATES + N_CONTROLS) + 5] = 1;
        vars[i * (N_STATES + N_CONTROLS) + 6].getValue().getDerivative()[i * (N_STATES + N_CONTROLS) + 6] = 1;

        vars[i * (N_STATES + N_CONTROLS) + 7].getValue().getDerivative()[i * (N_STATES + N_CONTROLS) + 7] = 1;
        vars[i * (N_STATES + N_CONTROLS) + 8].getValue().getDerivative()[i * (N_STATES + N_CONTROLS) + 8] = 1;
    }

    std::array<std::map<std::size_t, ascalar<real>>, N_STEPS * N_STATES> J;

    #ifdef TIMING
        auto start = std::chrono::steady_clock::now();
    #endif

    #pragma omp parallel for shared(J)

    for (size_t i = 0; i < N_STEPS - 1; i++) {

        scalar* Vx = &vars[i * (N_STATES + N_CONTROLS)];
        scalar* Vy = &vars[i * (N_STATES + N_CONTROLS) + 1];
        scalar* phip = &vars[i * (N_STATES + N_CONTROLS) + 2];
        scalar* rollp = &vars[i * (N_STATES + N_CONTROLS) + 3];
        scalar* roll = &vars[i * (N_STATES + N_CONTROLS) + 4];
        scalar* n = &vars[i * (N_STATES + N_CONTROLS) + 5];
        scalar* alpha = &vars[i * (N_STATES + N_CONTROLS) + 6];

        scalar* Cm = &vars[i * (N_STATES + N_CONTROLS) + 7];
        scalar* delta = &vars[i * (N_STATES + N_CONTROLS) + 8];

        scalar dFzas = -(KRollF * *roll + DRollF * *rollp) / wf;
        scalar dFzad = -dFzas;
        scalar dFzps = -(KRollR * *roll + DRollR * *rollp) / wr;
        scalar dFzpd = -dFzps;

        scalar Fzas = Fzf + dFzas;
        scalar Fzad = Fzf + dFzad;
        scalar Fzps = Fzr + dFzps;
        scalar Fzpd = Fzr + dFzpd;

        scalar Bfd = Kf * Fz0fPAC * sin(2 * atan(Fzad / (Kf2 * Fz0fPAC))) / (Cf * Df * Fzad);
        scalar Bfs = Kf * Fz0fPAC * sin(2 * atan(Fzas / (Kf2 * Fz0fPAC))) / (Cf * Df * Fzas);
        scalar Brd = Kr * Fz0rPAC * sin(2 * atan(Fzpd / (Kr2 * Fz0rPAC))) / (Cr * Dr * Fzpd);
        scalar Brs = Kr * Fz0rPAC * sin(2 * atan(Fzps / (Kr2 * Fz0rPAC))) / (Cr * Dr * Fzps);

        scalar Balphapd = Brd * atan((-*Vy + *phip * b) / (*Vx + *phip * wr));
        scalar Balphaad = Bfd * (*delta * deltamax + atan((-*Vy - *phip * a) / (*Vx + *phip * wf)));
        scalar Balphaps = Brs * atan((-*Vy + *phip * b) / (*Vx - *phip * wr));
        scalar Balphaas = Bfs * (*delta * deltamax + atan((-*Vy - *phip * a) / (*Vx - *phip * wf)));

        scalar ellFyad = 1 - pow((0.5 * *Cm * Camax * (1 - ratio) / (R * Fzad)), 2);
        scalar ellFyas = 1 - pow((0.5 * *Cm * Camax * (1 - ratio) / (R * Fzas)), 2);
        scalar ellFypd = 1 - pow((0.5 * *Cm * Camax * ratio / (R * Fzpd)), 2);
        scalar ellFyps = 1 - pow((0.5 * *Cm * Camax * ratio / (R * Fzps)), 2);

        scalar Fyad = Df * Fzad * sin(Cf * atan(Balphaad - Ef * (Balphaad - atan(Balphaad)))) * pow(ellFyad, 0.5);
        scalar Fyas = Df * Fzas * sin(Cf * atan(Balphaas - Ef * (Balphaas - atan(Balphaas)))) * pow(ellFyas, 0.5);
        scalar Fypd = Dr * Fzpd * sin(Cr * atan(Balphapd - Er * (Balphapd - atan(Balphapd)))) * pow(ellFypd, 0.5);
        scalar Fyps = Dr * Fzps * sin(Cr * atan(Balphaps - Er * (Balphaps - atan(Balphaps)))) * pow(ellFyps, 0.5);

        scalar cdelta = cos(*delta * deltamax);
        scalar sdelta = sin(*delta * deltamax);

        scalar croll = cos(*roll);
        scalar sroll = sin(*roll);

        scalar calpha = cos(*alpha);
        scalar salpha = sin(*alpha);

        scalar Fx = *Cm * Camax * ((1 - ratio) * cdelta + ratio) / R - (Fyad + Fyas) * sdelta;
        scalar Fy = Fypd + Fyps + (Fyad + Fyas) * cdelta + *Cm * Camax * (1 - ratio) * sdelta / R;
        scalar Mz = -(Fypd + Fyps) * b + (Fyad + Fyas) * a * cdelta + *Cm * Camax * (1 - ratio) * a * sdelta / R + (-Fyad + Fyas) * sdelta * wf;

        scalar nu = smax * (*Vx * calpha - *Vy * salpha) / (1 - *n * k);

        scalar f3 = (Mz - Fx * h * sroll) / (Jz * pow(croll, 2) + Jy * pow(sroll, 2));
        scalar f4 = Fy * h * croll / Jx + m * g * h * sroll / Jx + (Jy - Jz) * sroll * croll * pow(*phip, 2) / Jx - (KRollF + KRollR) * *roll / Jx - (DRollF + DRollR) * *rollp / Jx;
        scalar f1 = *phip * *Vy + Fx / m - cDrag * pow(*Vx, 2) / m - h * sroll * f3 - 2.0 * h * croll * *rollp * *phip;
        scalar f2 = -*phip * *Vx + Fy / m - h * sroll * pow(*phip, 2) + h * f4 * croll - h * sroll * pow(*rollp, 2);
        scalar f5 = *rollp;
        scalar f6 = *Vx * salpha + *Vy * calpha;
        scalar f7 = *phip - k * (*Vx * calpha - *Vy * salpha) / (1 - *n * k);

        scalar g1 = -vars[i * (N_STATES + N_CONTROLS) + (N_STATES + N_CONTROLS)] + *Vx + ds * f1 / nu;
        scalar g2 = -vars[i * (N_STATES + N_CONTROLS) + (N_STATES + N_CONTROLS) + 1] + *Vy + ds * f2 / nu;
        scalar g3 = -vars[i * (N_STATES + N_CONTROLS) + (N_STATES + N_CONTROLS) + 2] + *phip + ds * f3 / nu;
        scalar g4 = -vars[i * (N_STATES + N_CONTROLS) + (N_STATES + N_CONTROLS) + 3] + *rollp + ds * f4 / nu;
        scalar g5 = -vars[i * (N_STATES + N_CONTROLS) + (N_STATES + N_CONTROLS) + 4] + *roll + ds * f5 / nu;
        scalar g6 = -vars[i * (N_STATES + N_CONTROLS) + (N_STATES + N_CONTROLS) + 5] + *n + ds * f6 / nu;
        scalar g7 = -vars[i * (N_STATES + N_CONTROLS) + (N_STATES + N_CONTROLS) + 6] + *alpha + ds * f7 / nu;

        J[i * N_STATES] = g1.getDerivative();
        J[i * N_STATES + 1] = g2.getDerivative();
        J[i * N_STATES + 2] = g3.getDerivative();
        J[i * N_STATES + 3] = g4.getDerivative();
        J[i * N_STATES + 4] = g5.getDerivative();
        J[i * N_STATES + 5] = g6.getDerivative();
        J[i * N_STATES + 6] = g7.getDerivative();

    }

    #ifdef TIMING
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Total Time: " << elapsed.count() << " [ms]" << std::endl << std::endl;
    #endif
    
    for (auto j : J) {
        for (auto i = j.begin(); i != j.end(); i++)
            std::cout << i->first << " : " << i->second.getValue() << std::endl;
        std::cout << std::endl;
    }

    std::map<std::size_t, ascalar<real>> H0 = J[0];
    for (auto i = H0.begin(); i != H0.end(); i++) {
        for (auto j = i->second.getDerivative().begin(); j != i->second.getDerivative().end(); j++)
            std::cout << j->first << " : " << j->second << '\n';
        std::cout << std::endl;
    }
    
    return 0;

}