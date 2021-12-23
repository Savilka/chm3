#include <iostream>
#include "LOS.h"
#include "Timer.h"
#include <iomanip>

using namespace std;

int main() {
	
   ifstream* kuslau = new ifstream("BigTests/4545/kuslau.txt");
   ifstream* ig = new ifstream("BigTests/4545/ig.txt");
   ifstream* jg = new ifstream("BigTests/4545/jg.txt");
   ifstream* ggl = new ifstream("BigTests/4545/ggl.txt");
   ifstream* ggu = new ifstream("BigTests/4545/ggu.txt");
   ifstream* di = new ifstream("BigTests/4545/di.txt");
   ifstream* pr = new ifstream("BigTests/4545/pr.txt");
   ofstream* output = new ofstream("output.txt");
   ofstream* iteration = new ofstream("iteration.txt");

   cout << "1 - LOS" << endl << "2 - LOS-D" << endl << "3 - LOS-LU" << endl;
   LOS Los;
   int a,n;
   cin >> a;
   Timer time;
   switch (a) {
   case 1:
      //Los.input(*kuslau, *ig, *jg, *ggl, *ggu, *di, *pr);
      cin >> n;
      Los.Gilbert(n);
      Los.los(*iteration);
      Los.output(*output);
      break;
   case 2:
      //Los.input(*kuslau, *ig, *jg, *ggl, *ggu, *di, *pr);
      cin >> n;
      Los.Gilbert(n);
      Los.los_d(*iteration);
      Los.output(*output);
      break;
   case 3:
      //Los.input(*kuslau, *ig, *jg, *ggl, *ggu, *di, *pr);
      cin >> n;
      Los.Gilbert(n);
      Los.los_LU(*iteration);
      Los.output(*output);
      break;
   default:
      break;
   }
   return 0;
}