#include <iostream>
#include <chrono>
#include <array>
#include <iomanip>

#include "adcg_ascalar.hpp"
#include "adcg_binary_ops.hpp"
#include "adcg_math_ops.hpp"

typedef ascalar<real> scalar;

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
    double k = 0.1;
    double ds = 0.002;
    
    scalar Vx(94 / 3.6, 0);
    scalar Vy(0.0, 1);
    scalar phip(0.0, 2);
    scalar rollp(0.0, 3);
    scalar roll(0.0, 4);
    scalar n(0.0, 5);
    scalar alpha(0.0, 6);

    scalar Cm(0.0, 7);
    scalar delta(0.0, 8);

    std::cout << "Vx:       " << Vx.getValue() << std::endl;
    std::cout << "Vy:       " << Vy.getValue() << std::endl;
    std::cout << "Phip:     " << phip.getValue() << std::endl;
    std::cout << "Rollp:    " << rollp.getValue() << std::endl;
    std::cout << "Roll:     " << roll.getValue() << std::endl;
    std::cout << "N:        " << n.getValue() << std::endl;
    std::cout << "Alpha:    " << alpha.getValue() << std::endl;
    std::cout << "Cm:       " << Cm.getValue() << std::endl;
    std::cout << "Delta:    " << delta.getValue() << std::endl;
    std::cout << std::endl;

    #ifdef TIMING
        auto start = std::chrono::steady_clock::now();
    #endif

    scalar dFzas = -(KRollF * roll + DRollF * rollp) / wf;
    scalar dFzad = -dFzas;
    scalar dFzps = -(KRollR * roll + DRollR * rollp) / wr;
    scalar dFzpd = -dFzps;

    scalar Fzas = Fzf + dFzas;
    scalar Fzad = Fzf + dFzad;
    scalar Fzps = Fzr + dFzps;
    scalar Fzpd = Fzr + dFzpd;

    scalar Bfd = Kf * Fz0fPAC * sin(2 * atan(Fzad / (Kf2 * Fz0fPAC))) / (Cf * Df * Fzad);
    scalar Bfs = Kf * Fz0fPAC * sin(2 * atan(Fzas / (Kf2 * Fz0fPAC))) / (Cf * Df * Fzas);
    scalar Brd = Kr * Fz0rPAC * sin(2 * atan(Fzpd / (Kr2 * Fz0rPAC))) / (Cr * Dr * Fzpd);
    scalar Brs = Kr * Fz0rPAC * sin(2 * atan(Fzps / (Kr2 * Fz0rPAC))) / (Cr * Dr * Fzps);

    scalar Balphapd = Brd * atan((-Vy + phip * b) / (Vx + phip * wr));
    scalar Balphaad = Bfd * (delta * deltamax + atan((-Vy - phip * a) / (Vx + phip * wf)));
    scalar Balphaps = Brs * atan((-Vy + phip * b) / (Vx - phip * wr));
    scalar Balphaas = Bfs * (delta * deltamax + atan((-Vy - phip * a) / (Vx - phip * wf)));

    scalar ellFyad = 1 - pow((0.5 * Cm * Camax * (1 - ratio) / (R * Fzad)), 2);
    scalar ellFyas = 1 - pow((0.5 * Cm * Camax * (1 - ratio) / (R * Fzas)), 2);
    scalar ellFypd = 1 - pow((0.5 * Cm * Camax * ratio / (R * Fzpd)), 2);
    scalar ellFyps = 1 - pow((0.5 * Cm * Camax * ratio / (R * Fzps)), 2);

    scalar Fyad = Df * Fzad * sin(Cf * atan(Balphaad - Ef * (Balphaad - atan(Balphaad)))) * pow(ellFyad, 0.5);
    scalar Fyas = Df * Fzas * sin(Cf * atan(Balphaas - Ef * (Balphaas - atan(Balphaas)))) * pow(ellFyas, 0.5);
    scalar Fypd = Dr * Fzpd * sin(Cr * atan(Balphapd - Er * (Balphapd - atan(Balphapd)))) * pow(ellFypd, 0.5);
    scalar Fyps = Dr * Fzps * sin(Cr * atan(Balphaps - Er * (Balphaps - atan(Balphaps)))) * pow(ellFyps, 0.5);

    scalar cdelta = cos(delta * deltamax);
    scalar sdelta = sin(delta * deltamax);

    scalar croll = cos(roll);
    scalar sroll = sin(roll);

    scalar calpha = cos(alpha);
    scalar salpha = sin(alpha);

    scalar Fx = Cm * Camax * ((1 - ratio) * cdelta + ratio) / R - (Fyad + Fyas) * sdelta;
    scalar Fy = Fypd + Fyps + (Fyad + Fyas) * cdelta + Cm * Camax * (1 - ratio) * sdelta / R;
    scalar Mz = -(Fypd + Fyps) * b + (Fyad + Fyas) * a * cdelta + Cm * Camax * (1 - ratio) * a * sdelta / R + (-Fyad + Fyas) * sdelta * wf;

    scalar nu = smax * (Vx * calpha - Vy * salpha) / (1 - n * k);

    scalar f3 = (Mz - Fx * h * sroll) / (Jz * pow(croll, 2) + Jy * pow(sroll, 2));
    scalar f4 = Fy * h * croll / Jx + m * g * h * sroll / Jx + (Jy - Jz) * sroll * croll * pow(phip, 2) / Jx - (KRollF + KRollR) * roll / Jx - (DRollF + DRollR) * rollp / Jx;
    scalar f1 = phip * Vy + Fx / m - cDrag * pow(Vx, 2) / m - h * sroll * f3 - 2.0 * h * croll * rollp * phip;
    scalar f2 = -phip * Vx + Fy / m - h * sroll * pow(phip, 2) + h * f4 * croll - h * sroll * pow(rollp, 2);
    scalar f5 = rollp;
    scalar f6 = Vx * salpha + Vy * calpha;
    scalar f7 = phip - k * (Vx * calpha - Vy * salpha) / (1 - n * k);

    scalar g1 = f1 / nu;
    scalar g2 = f2 / nu;
    scalar g3 = f3 / nu;
    scalar g4 = f4 / nu;
    scalar g5 = f5 / nu;
    scalar g6 = f6 / nu;
    scalar g7 = f7 / nu;

    #ifdef TIMING
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Total Time: " << elapsed.count() << " [μs]" << std::endl << std::endl;
    #endif

    std::cout << g1.getValue() << std::endl;
    std::cout << g2.getValue() << std::endl;
    std::cout << g3.getValue() << std::endl;
    std::cout << g4.getValue() << std::endl;
    std::cout << g5.getValue() << std::endl;
    std::cout << g6.getValue() << std::endl;
    std::cout << g7.getValue() << std::endl;
    std::cout << std::endl;

    std::map<std::size_t, real> J1 = g1.getDerivative();
    std::map<std::size_t, real> J2 = g2.getDerivative();
    std::map<std::size_t, real> J3 = g3.getDerivative();
    std::map<std::size_t, real> J4 = g4.getDerivative();
    std::map<std::size_t, real> J5 = g5.getDerivative();
    std::map<std::size_t, real> J6 = g6.getDerivative();
    std::map<std::size_t, real> J7 = g7.getDerivative();
    
    for (auto i = J1.begin(); i != J1.end(); i++)
        std::cout << i->first << " : " << i->second << std::endl;
    std::cout << std::endl;

    for (auto i = J2.begin(); i != J2.end(); i++)
        std::cout << i->first << " : " << i->second << std::endl;
    std::cout << std::endl;

    for (auto i = J3.begin(); i != J3.end(); i++)
        std::cout << i->first << " : " << i->second << std::endl;
    std::cout << std::endl;

    for (auto i = J4.begin(); i != J4.end(); i++)
        std::cout << i->first << " : " << i->second << std::endl;
    std::cout << std::endl;

    for (auto i = J5.begin(); i != J5.end(); i++)
        std::cout << i->first << " : " << i->second << std::endl;
    std::cout << std::endl;

    for (auto i = J6.begin(); i != J6.end(); i++)
        std::cout << i->first << " : " << i->second << std::endl;
    std::cout << std::endl;

    for (auto i = J7.begin(); i != J7.end(); i++)
        std::cout << i->first << " : " << i->second << std::endl;
    std::cout << std::endl;

    std::array<std::map<std::size_t, real>, 7> spJ = {J1, J2, J3, J4, J5, J6, J7};
    std::array<real, 7 * 9> J{0};

    densifyDerivative(spJ, J);

    for (auto i = 0; i < 7; i++) {
        for (auto j = 0 + i * 9; j < 9 + i * 9; j++)
            std::cout << std::fixed << std::setprecision(6) << std::setw(9) << std::setfill(' ') << J[j] << '\t';
        std::cout << std::endl;
    }
    
    return 0;

}