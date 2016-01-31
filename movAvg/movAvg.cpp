#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <queue>

using namespace std;

int main (int argc, char *argv[]) {

	if (argc != 3) {
		cerr << "Usage Error: <window size> <file name>" << endl;
		exit(EXIT_FAILURE);
	}
	
	ifstream in(argv[2]);
	
	if (!in.is_open()) {
		cerr << "I/O Error: Unable to open " << argv[2] << "!" << endl;
		exit(EXIT_FAILURE);
	}
	
	int win_size = atoi(argv[1]);
	queue<double> avg;
	double val, sum = 0;
	for (int i = 0; i < win_size; i++) {
		in >> val;
		//cout << "NAN" << endl;

		if (in.eof() && (i < (win_size-1))) {
			cerr << "Logic Error: File shorter than the window size!" << endl;
			exit(EXIT_FAILURE);
		}
		sum += val;
		avg.push(val);
		cout << sum/(double)(i+1) << endl;
	}
	
	while (!in.eof()) {
		in >> val;
		if (in.eof()) break;
		avg.push(val);
		sum -= avg.front();
		sum += val;
		avg.pop();
		cout << sum/(double)win_size << endl;
	}
	in.close();
}
