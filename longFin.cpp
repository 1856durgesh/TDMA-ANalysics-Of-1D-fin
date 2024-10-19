#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

// Class for solving heat conduction problem using TDMA
class HeatConductionFin
{
private:
  double rho, Cp, k, h, r, P, A, T_inf, L, dx, dt, del_p;
  int N, Nt;
  double a_w, a_e, a_p;
  vector<double> T;

public:
  // Constructor to initialize the parameters
  HeatConductionFin(double rho, double Cp, double k, double h, double r, double T_inf, double L, double dx, double dt, int Nt)
      : rho(rho), Cp(Cp), k(k), h(h), r(r), T_inf(T_inf), L(L), dx(dx), dt(dt), Nt(Nt)
  {
    P = 2 * 3.18 * r; // Perimeter of the circular fin
    A = 3.18 * r * r; // Area of the circular fin
    del_p = h * P * T_inf / (rho * Cp * A);
    N = L / dx; // Number of nodes
    a_w = -(k * dt) / (rho * Cp * dx * dx);
    a_e = a_w;
    a_p = 1.0 - 2 * a_w + h * P / (rho * Cp * A);

    // Initialize temperature distribution
    T.resize(N + 1, T_inf);
    T[0] = 100; // Base temperature of the fin
  }

  // Function to perform TDMA algorithm for solving tridiagonal matrix
  vector<double> TDMA(const vector<double> &a, const vector<double> &b, const vector<double> &c, const vector<double> &d)
  {
    int n = b.size();
    vector<double> Tn1(n);    // Solution vector
    vector<double> cp(n - 1); // Modified upper diagonal
    vector<double> dp(n);     // Modified RHS vector

    // Forward sweep
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    for (int i = 1; i < n; i++)
    {
      double m = 1.0 / (b[i] - a[i - 1] * cp[i - 1]); // Making the diagonal element to 1

      if (i < n - 1)
      {
        cp[i] = c[i] * m;
      }
      dp[i] = (d[i] - a[i - 1] * dp[i - 1]) * m;
    }

    // Backward substitution
    Tn1[n - 1] = dp[n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
      Tn1[i] = dp[i] - cp[i] * Tn1[i + 1];
    }

    return Tn1;
  }

  // Function to solve the heat conduction problem
  void solve()
  {
    vector<double> a(N - 1, 0); // Lower row of the TDMA matrix
    vector<double> b(N, 0);     // Middle row of the TDMA matrix
    vector<double> c(N - 1, 0); // Upper row of the TDMA matrix
    vector<double> d(N, 0);     // RHS matrix

    // Time loop
    for (int n = 0; n < Nt; n++)
    {
      // Interior nodes
      for (int i = 0; i < N - 1; i++)
      {
        a[i] = a_w;
        c[i] = a_e;
        b[i] = a_p;
        d[i] = T[i] + del_p;
      }

      b[0] -= a_w;
      b[N - 1] = a_p - a_e;
      d[0] -= 2 * a_w * 100;
      int T2=T_inf;
      d[N - 1] = T[N-1]+del_p-2*T2*a_e; // set temperature at the last node to ambient temperature 

      // Solving the TDMA system
      vector<double> Tn1 = TDMA(a, b, c, d);

      // Update temperature distribution
      for (int i = 0; i < N; i++)
      {
        T[i] = Tn1[i];
      }
    }
  }



  // Function to display results
void displayResults() 
  {
    // Create and open a text file
    ofstream outFile("LongFinResults.txt");

    if (outFile.is_open())
    {
      for (int i = 0; i < N; i++)
      {
        // Write the results to the file
        outFile << "x = " << i * dx +dx/2 << ", T = " << T[i] << endl;
      }

      outFile.close(); // Close the file when done
      cout << "Results have been saved to LongFinResults.txt" << endl;
    }
    else
    {
      cout << "Unable to open the file." << endl;
    }
  }
};

int main()
{
  // Initializing the parameters
  double rho = 7850;  // Density of aluminum
  double Cp = 473;    // Specific heat capacity
  double k = 40;      // Thermal conductivity
  double h = 10;      // Convective heat transfer coefficient
  double r = 0.00318; // Radius of the fin
  double T_inf = 25;  // Ambient temperature
  double L = 0.1;     // Fin length
  double dx = 0.01;   // Spatial step size
  double dt = 0.01;   // Time step size
  int Nt = 10000;   // Number of time steps

  // Creating an object of the HeatConductionFin class
  HeatConductionFin fin(rho, Cp, k, h, r, T_inf, L, dx, dt, Nt);

  // Solving the heat conduction problem
  fin.solve();

  // Displaying the results
  fin.displayResults();

  return 0;
}
