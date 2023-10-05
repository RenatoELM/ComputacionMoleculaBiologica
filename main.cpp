#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>
#include <stdlib.h>
#include <fstream>

#define MAX_N 1001

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef long long lld;
typedef unsigned long long llu;

int n, m;
int match_score, mismatch_score, gap_score;
string s, t, u, v, w, x;
string sr, tr, ur, vr, wr, xr;
int dp[MAX_N][MAX_N];

void lec_txt() {
    std::ifstream archivoEntrada("D:\\UCSP\\2023-02\\Computación Molecular Biológica\\Tarea 1\\Sars-Cov(1080).txt");
    std::ofstream archivoSalida("D:\\UCSP\\2023-02\\Computación Molecular Biológica\\Tarea 1\\Sars-CovFin(1080).txt");
    if (archivoEntrada.is_open() && archivoSalida.is_open()) {
        std::string linea;
        while (std::getline(archivoEntrada, linea)) {
            // Elimina los espacios en blanco de la línea
            linea.erase(std::remove_if(linea.begin(), linea.end(), ::isspace), linea.end());

            // Escribe la línea procesada en el archivo de salida
            archivoSalida << linea << std::endl;
        }

        // Cierra los archivos
        archivoEntrada.close();
        archivoSalida.close();
    }
    else {
        std::cerr << "No se pudo abrir uno o ambos archivos." << std::endl;
    }
}

inline int needleman_wunsch(string s,string t)
{
    for (int i = 0; i <= n; i++) dp[i][0] = dp[0][i] = -i * gap_score;
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= m; j++)
        {
            int S = (s[i - 1] == t[j - 1]) ? match_score : -mismatch_score;
            dp[i][j] = max(dp[i - 1][j - 1] + S, max(dp[i - 1][j] - gap_score, dp[i][j - 1] - gap_score));
        }
    }
    return dp[n][m];
}

inline pair<string, string> get_optimal_alignment(string s, string t)
{
    string retA, retB;
    stack<char> SA, SB;
    int ii = n, jj = m;
    while (ii != 0 || jj != 0)
    {
        if (ii == 0)
        {
            SA.push('-');
            SB.push(t[jj - 1]);
            jj--;
        }
        else if (jj == 0)
        {
            SA.push(s[ii - 1]);
            SB.push('-');
            ii--;
        }
        else
        {
            int S = (s[ii - 1] == t[jj - 1]) ? match_score : -mismatch_score;
            if (dp[ii][jj] == dp[ii - 1][jj - 1] + S)
            {
                SA.push(s[ii - 1]);
                SB.push(t[jj - 1]);
                ii--; jj--;
            }
            else if (dp[ii - 1][jj] > dp[ii][jj - 1])
            {
                SA.push(s[ii - 1]);
                SB.push('-');
                ii--;
            }
            else
            {
                SA.push('-');
                SB.push(t[jj - 1]);
                jj--;
            }
        }
    }
    while (!SA.empty())
    {
        retA += SA.top();
        retB += SB.top();
        SA.pop();
        SB.pop();
    }
    return make_pair(retA, retB);
}

int main()
{
    //lec_txt();

    //n = 8, m = 7;

    match_score = 2, mismatch_score = 1, gap_score = 1;

    //s = "ATCTTCTT";
    //t = "ACTGACC";

    //ifstream archivo1("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 1/Pruebas/SecuenciaA.txt");
    //ifstream archivo2("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 1/Pruebas/SecuenciaB.txt");

    ifstream archivo1("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 3/Cadenas/A BRAC1 F.txt");
    ifstream archivo2("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 3/Cadenas/B BRAC1 F.txt");
    ifstream archivo3("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 3/Cadenas/C BRAC1 F.txt");
    ifstream archivo4("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 3/Cadenas/D BRAC1 F.txt");
    ifstream archivo5("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 3/Cadenas/E BRAC1 F.txt");
    ifstream archivo6("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 3/Cadenas/F BRAC1 F.txt");

    ifstream archivo7("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 3/Cadenas/A BRAC1 R.txt");
    ifstream archivo8("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 3/Cadenas/B BRAC1 R.txt");
    ifstream archivo9("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 3/Cadenas/C BRAC1 R.txt");
    ifstream archivo10("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 3/Cadenas/D BRAC1 R.txt");
    ifstream archivo11("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 3/Cadenas/E BRAC1 R.txt");
    ifstream archivo12("D:/UCSP/2023-02/Computación Molecular Biológica/Tarea 3/Cadenas/F BRAC1 R.txt");

    while (!archivo1.eof()) {
        getline(archivo1, s);
    }

    while (!archivo2.eof()) {
        getline(archivo2, t);
    }

    while (!archivo3.eof()) {
        getline(archivo3, u);
    }

    while (!archivo4.eof()) {
        getline(archivo4, v);
    }

    while (!archivo5.eof()) {
        getline(archivo5, w);
    }

    while (!archivo6.eof()) {
        getline(archivo6, x);
    }

    //

    while (!archivo7.eof()) {
        getline(archivo7, sr);
    }

    while (!archivo8.eof()) {
        getline(archivo8, tr);
    }

    while (!archivo9.eof()) {
        getline(archivo9, ur);
    }

    while (!archivo10.eof()) {
        getline(archivo10, vr);
    }

    while (!archivo11.eof()) {
        getline(archivo11, wr);
    }

    while (!archivo12.eof()) {
        getline(archivo12, xr);
    }

    archivo1.close();
    archivo2.close();
    archivo3.close();
    archivo4.close();
    archivo5.close();
    archivo6.close();
    archivo7.close();
    archivo8.close();
    archivo9.close();
    archivo10.close();
    archivo11.close();
    archivo12.close();

    n = size(t), m = size(s);

    //printf("%d\n", needleman_wunsch(t, s));
    needleman_wunsch(t, s);
    pair<string, string> alignment = get_optimal_alignment(t, s);
    printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());

    string a = alignment.first.c_str();
    string b = alignment.second.c_str();

    m = size(u);

    //printf("%d\n", needleman_wunsch(t, u));
    needleman_wunsch(t, u);
    alignment = get_optimal_alignment(t, u);
    printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());

    string c = alignment.second.c_str();

    m = size(v);

    //printf("%d\n", needleman_wunsch(t, v));
    needleman_wunsch(t, v);
    alignment = get_optimal_alignment(t, v);
    printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());

    string d = alignment.second.c_str();

    m = size(w);

    //printf("%d\n", needleman_wunsch(t, w));
    needleman_wunsch(t, w);
    alignment = get_optimal_alignment(t, w);
    printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());

    string e = alignment.second.c_str();

    m = size(x);

    //printf("%d\n", needleman_wunsch(t, x));
    needleman_wunsch(t, x);
    alignment = get_optimal_alignment(t, x);
    printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());

    string f = alignment.second.c_str();

    //

    m = size(sr);

    //printf("%d\n", needleman_wunsch(t, s));
    needleman_wunsch(t, sr);
    alignment = get_optimal_alignment(t, sr);
    printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());

    string g = alignment.second.c_str();

    m = size(tr);

    //printf("%d\n", needleman_wunsch(t, s));
    needleman_wunsch(t, tr);
    alignment = get_optimal_alignment(t, tr);
    printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());

    string h = alignment.second.c_str();

    m = size(ur);

    //printf("%d\n", needleman_wunsch(t, u));
    needleman_wunsch(t, ur);
    alignment = get_optimal_alignment(t, ur);
    printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());

    string i = alignment.second.c_str();

    m = size(vr);

    //printf("%d\n", needleman_wunsch(t, v));
    needleman_wunsch(t, vr);
    alignment = get_optimal_alignment(t, vr);
    printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());

    string j = alignment.second.c_str();

    m = size(wr);

    //printf("%d\n", needleman_wunsch(t, w));
    needleman_wunsch(t, wr);
    alignment = get_optimal_alignment(t, wr);
    printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());

    string k = alignment.second.c_str();

    m = size(xr);

    //printf("%d\n", needleman_wunsch(t, x));
    needleman_wunsch(t, xr);
    alignment = get_optimal_alignment(t, xr);
    printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());

    string l = alignment.second.c_str();

    vector<string> vec1 = { a, b, c, d, e, f, g, h, i, j, k, l };
    vector<string> vec2 = { s, t, u, v, w, x, sr, tr, ur, vr, wr, xr };

    cout << endl;
    cout << "Cadenas:" << endl;

    for (int z = 0; z < 12; z++) {
        cout << vec2[z] << endl;
    }

    cout << "\nMSA:" << endl;

    int tam;
    for (int z = 0; z < 12; z++) {
        if (z == 0) {
            tam = size(vec1[z]);
        }
        else {
            if (tam < size(vec1[z])) {
                tam = size(vec1[z]);
            }
        }
    }

    for (int z = 0; z < 12; z++) {
        if (size(vec1[z]) < tam) {
            while (size(vec1[z]) != tam) {
                vec1[z] = vec1[z] + "-";
            }
        }
    }

    for (int z = 0; z < 12; z++) {
        cout << vec1[z] << "\n";
    }
    return 0;
}