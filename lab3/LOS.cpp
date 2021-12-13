#include "LOS.h"

void LOS::r0() {
   matMul(x);
   real* Ax = new real[N];
   copyVec(Ax, res, N);
   for (int i = 0; i < N; i++)
   {
      r[i] = pr[i] - Ax[i];
   }
}

real LOS::scalar(real* a, real* b, int size) {
   real res = 0;
   for (int i = 0; i < size; i++)
      res += a[i] * b[i];
   return res;
}

void LOS::xk() {
   for (int i = 0; i < N; i++) {
      x[i] += alpha * z[i];
   }
}

void LOS::rk() {
   for (int i = 0; i < N; i++) {
      r[i] -= alpha * p[i];
   }
}

void LOS::zk() {
   for (int i = 0; i < N; i++) {
      z[i] = r[i] + beta * z[i];
   }
}

void LOS::pk() {
   matMul(r);
   real* Ak = new real[N];
   copyVec(Ak, res, N);
   for (int i = 0; i < N; i++){
      p[i] = Ak[i] + beta * p[i];
   }
}

real LOS::norm(real* vec) {
   real res = 0;
   for (int i = 0; i < N; i++) {
      res += vec[i] * vec[i];
   }
   return sqrt(res);
}

void LOS::calcAlpha() {
   alpha = scalar(p, r, N) / scalar(p, p, N);
}

void LOS::calcBeta() {
   matMul(r);
   beta = scalar(p, res, N) / scalar(p, p, N);
}

void LOS::matMul(real* vec) {
   for (int i = 0; i < N; i++) {
      res[i] = 0;
   }
   for (int i = 0; i < N; i++) {
      res[i] = di[i] * vec[i];
      for (int j = ig[i]; j < ig[i + 1]; j++) {
         res[i] += ggl[j] * vec[jg[j]];
         res[jg[j]] += ggu[j] * vec[i];
      }
   }
}

void LOS::copyVec(real* a, real* b, int size) {
   for (int i = 0; i < size; i++) {
      a[i] = b[i];
   }
}

void LOS::los() {
   res = new real[N];
   x = new real[N];
   r = new real[N];
   p = new real[N];
   z = new real[N];
   for (int i = 0; i < N; i++) {
      x[i] = 0;
   }
   r0();
   copyVec(z, r, N);
   matMul(z);
   copyVec(p, res, N);
   nev = scalar(r, r, N);
   for (int i = 0; i < maxiter && nev > eps; i++) {
      calcAlpha();
      xk();
      cout << "Iteration number: " << i << " | " << "Squared norm residuals: " << nev << endl;
      nev = scalar(r, r, N) - alpha * alpha * scalar(p, p, N);
      rk();
      calcBeta();
      zk();
      pk();
   }
}

void LOS::input(ifstream& fkuslau, ifstream& fig, ifstream& fjg, ifstream& fggl, ifstream& fggu, ifstream& fdi, ifstream& fpr) {
   fkuslau >> N;
   fkuslau >> maxiter;
   fkuslau >> eps;

   ig = new int[N + 1];
   for (int i = 0; i < N + 1; i++) {
      fig >> ig[i];
      ig[i]--;
   }

   di = new real[N];
   for (int i = 0; i < N; i++) {
      fdi >> di[i];
   }

   jg = new int[ig[N]];
   for (int i = 0; i < ig[N]; i++) {
      fjg >> jg[i];
      jg[i]--;
   }

   ggl = new real[ig[N]];
   ggu = new real[ig[N]];
   for (int j = 0; j < ig[N]; j++) {
      fggl >> ggl[j];
      fggu >> ggu[j];
   }
 

   pr = new real[N];
   for (int i = 0; i < N; i++) {
      fpr >> pr[i];
   }
}

void LOS::output(ofstream& output) {
   for (int i = 0; i < N; i++) {
      output << x[i] << endl;
   }
}


