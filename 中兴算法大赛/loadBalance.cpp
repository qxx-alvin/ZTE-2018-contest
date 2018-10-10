#include "readData.h"
#include "ga.h"
#include <iostream>
#include <time.h>
using namespace std;

void main()
{
	clock_t start_t = clock();
	readTopoAndRequest("F:\\中兴算法大赛\\迪杰斯特拉\\\case1\\grditopoAndRequest.txt");
	ga();
	clock_t end_t = clock();
	cout << "time elapsed: " << double(end_t - start_t) / CLOCKS_PER_SEC << "s" << endl;


	cin.get();
}