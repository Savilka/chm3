#include <iostream>
#include "LOS.h"
#include "Timer.h"
#include <iomanip>

using namespace std;

int main() {
	
   ifstream* kuslau = new ifstream("kuslau.txt");
   ifstream* ig = new ifstream("ig.txt");
   ifstream* jg = new ifstream("jg.txt");
   ifstream* ggl = new ifstream("ggl.txt");
   ifstream* ggu = new ifstream("ggu.txt");
   ifstream* di = new ifstream("di.txt");
   ifstream* pr = new ifstream("pr.txt");
   ofstream* output = new ofstream("output.txt");
   ofstream* iteration = new ofstream("iteration.txt");
   Timer time;
   LOS Los;
   Los.input(*kuslau, *ig, *jg, *ggl, *ggu, *di, *pr);
   Los.los_d(*iteration);
   Los.output(*output);
   return 0;
}