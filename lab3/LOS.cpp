#include "LOS.h"

void LOS::r0() {
   real* Ax = new real[N];
   matMul(x, Ax);
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

void LOS::zk(real* buf) {
   for (int i = 0; i < N; i++) {
      z[i] = buf[i] + beta * z[i];
   }
}

void LOS::pk(real* Ak) {
   
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

void LOS::calcBeta(real* Ar) {
   beta = - scalar(p, Ar, N) / scalar(p, p, N);
}

void LOS::matMul(real* vec, real* res) {
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


void LOS::vecDivD(real* res) {
   for (int i = 0; i < N; i++) {
      res[i] = res[i] / sqrt(di[i]);
   }
}

void LOS::Gilbert(int n) {
   N = n;
   matrixG = new real * [N];
   for (int i = 0; i < N; i++) {
      matrixG[i] = new real[N];
   }
   pr = new real[N];
   for (int i = 1; i <= N; i++) {
      pr[i - 1] = 0;
      for (int j = 1; j <= N; j++) {
         matrixG[i - 1][j - 1] = 1 / double(i + j - 1);
         pr[i - 1] += matrixG[i - 1][j - 1] * (j);
      }
   }
   eps = 1e-15;
   maxiter = 10000;
   prof();
}

void LOS::prof() {
   ig = new int[100];
   jg = new int[100];
   ggl = new real[100];
   ggu = new real[100];
   di = new real[100];
   ig[0] = 0;
   for (int i = 0; i < N; i++) {
      ig[i + 1] = ig[i] + i;
      for (int j = 0; j < i; j++) {
         ggl[ig[i] + j] = matrixG[i][j];
         ggu[ig[i] + j] = ggl[ig[i] + j];
         jg[ig[i] + j] = j;
      }
      di[i] = matrixG[i][i];
   }
}

void LOS::LU(real* L, real* U, real* D) {
   copyVec(L, ggl, ig[N]);
   copyVec(U, ggu, ig[N]);
   copyVec(D, di, N);

   for (int i = 0; i < N; i++) {
      real sumdi = 0.0;	

      int i0 = ig[i];
      int i1 = ig[i + 1];

    
      for (int k = i0; k < i1; k++) {
         int j = jg[k];
         int j0 = ig[j];
                                
         int j1 = ig[j + 1];	
                                

         int ik = i0;			
         int kj = j0;			

         real suml = 0.0;		
         real sumu = 0.0;		

         while (ik < k && kj < j1) {
            
            if (jg[ik] == jg[kj]) {
               
               suml += L[ik] * U[kj];
               sumu += U[ik] * L[kj];
               ik++;
               kj++;
            }
            
            else
               jg[ik] > jg[kj] ? kj++ : ik++;
         }

         
         L[k] -= suml;
         U[k] = (U[k] - sumu) / D[j];
         sumdi += L[k] * U[k];
      }

      		
      D[i] -= sumdi;
   }
}

void LOS::Forward(real* vec, real* res) {
   copyVec(res, vec, N);
   for (int i = 0; i < N; i++) {
      real sum = 0.0;
      for (int j = ig[i]; j < ig[i + 1]; j++)
         sum += L[j] * res[jg[j]];

      res[i] -= sum;
      res[i] /= D[i];
   }
}

void LOS::Backward(real* vec, real* res) {
   copyVec(res, vec, N);
   for (int i = N - 1; i >= 0; i--) {
      for (int j = ig[i]; j < ig[i + 1]; j++)
         res[jg[j]] -= U[j] * res[i];
   }
}

void LOS::los(ofstream& iteration) {
   x = new real[N];
   r = new real[N];
   p = new real[N];
   z = new real[N];
   real diff = 1;
   real* Ak = new real[N];
   real* Ar = new real[N];
   for (int i = 0; i < N; i++) {
      x[i] = 0;
   }
   int i = 0;
   r0();
   copyVec(z, r, N);
   matMul(z, p);
   nev = scalar(r, r, N);
   for (; i < maxiter && nev > eps; i++) {
      calcAlpha();
      xk();
      diff = abs(nev - nev_next);
      iteration << "Iteration number: " << i << " | " << "Squared norm residuals: " << nev << endl;
      nev_next = nev;
      nev = scalar(r, r, N) - alpha * alpha * scalar(p, p, N);
      rk();
      matMul(r, Ar);
      calcBeta(Ar);
      zk(r);
      matMul(r, Ak);
      pk(Ak);
   }
   iteration << "Iteration number: " << i << " | " << "Squared norm residuals: " << nev << endl;
}

void LOS::los_d(ofstream& iteration) {
   x = new real[N];
   r = new real[N];
   p = new real[N];
   z = new real[N];
   real diff = 1;
   real* Ar = new real[N];
   for (int i = 0; i < N; i++) {
      x[i] = 0;
   }
   int i = 0;
   r0();
   vecDivD(r);
   copyVec(z, r, N);
   matMul(z, p);
   vecDivD(p);
   nev = scalar(r, r, N);
   for (; i < maxiter && nev > eps; i++) {
      calcAlpha();
      xk();
      diff = abs(nev - nev_next);
      iteration <<  "Iteration number: " << i << " | " << "Squared norm residuals: " << nev << endl;
      nev_next = nev;
      nev = scalar(r, r, N) - alpha * alpha * scalar(p, p, N);
      rk();
      matMul(r, Ar);
      vecDivD(Ar);
      calcBeta(Ar);
      zk(r);
      pk(Ar);
   }
   iteration << "Iteration number: " << i << " | " << "Squared norm residuals: " << nev << endl;
}

void LOS::los_LU(ofstream& iteration) {
   L = new real[ig[N]];
   U = new real[ig[N]];
   D = new real[N];
   x = new real[N];
   r = new real[N];
   p = new real[N];
   z = new real[N];
   real* bufL = new real[N];
   real* bufU = new real[N];
   real* buf = new real[N];
   real diff = 1;
   for (int i = 0; i < N; i++) {
      x[i] = 0;
   }
   int i = 0;
   LU(L, U, D);
   r0();
   Forward(r, r);
   Backward(r, z);
   matMul(z, p);
   Forward(p, p);
   nev = scalar(r, r, N);
   for (; i < maxiter && nev > eps; i++) {
      calcAlpha();
      xk();
      diff = abs(nev - nev_next);
      iteration << "Iteration number: " << i << " | " << "Squared norm residuals: " << nev << endl;
      nev_next = nev;
      nev = scalar(r, r, N) - alpha * alpha * scalar(p, p, N);
      rk();
      Backward(r, bufU);
      matMul(bufU, buf);
      Forward(buf, bufU);
      calcBeta(bufU);
      Backward(r, bufL);
      zk(bufL);
      pk(bufU);
   }
   iteration << "Iteration number: " << i << " | " << "Squared norm residuals: " << nev << endl;
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
      
      output << setprecision(14) << x[i] << endl;
   }
}


