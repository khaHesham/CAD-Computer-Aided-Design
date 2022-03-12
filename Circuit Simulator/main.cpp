#include <Eigen/Dense>
#include<string.h>
#include <vector>
#include <complex>
#include<iostream>
#include <Windows.h>
#include<fstream>
#include <iomanip>

#define _USE_MATH_DEFINES
#include <math.h>


using namespace Eigen;
using namespace std;


struct Component
{
	string type;
	string name = "";
	complex<double> voltage = 0;
	complex<double> current = 0;
	complex<double> Z = 0;
	double coeff = 0;
	int node1 = 0, node2 = 0, node3 = 0, node4 = 0;
};

vector<Component*> Components;

int WhichPassive(int n1, int n2)
{
	for (int i = 0; i < Components.size(); i++)
	{
		if ((Components[i]->node1 == n1 && Components[i]->node2 == n2) || (Components[i]->node1 == n2 && Components[i]->node2 == n1))
			return i;
	}
	return -1;
}



int main()
{
	//Opening the file
	string FileName;
	cout << "Please enter the INPUT file name: ";
	cin >> FileName;
	ifstream InputFile;
	string path = "../Circuits/" + FileName + ".txt";
	InputFile.open(path);
	if (!InputFile.is_open())
	{
		cout << "File not found";
		return 1;
	}

	//Reading data from file
	char ignore;
	double W;
	InputFile >> ignore >> W;

	int NodeCount = 0;
	while (!InputFile.eof())
	{
		Component* cmp = new Component;
		Components.push_back(cmp);


		InputFile >> cmp->type >> cmp->name >> cmp->node1 >> cmp->node2;

		if (cmp->node1 > NodeCount) NodeCount++;
		if (cmp->node2 > NodeCount) NodeCount++;

		if (cmp->type == "vsrc")
		{
			double volt, phase;
			InputFile >> volt >> phase;
			cmp->voltage = polar(volt, phase * M_PI / 180);
		}

		else if (cmp->type == "isrc")
		{
			double current, phase;
			InputFile >> current >> phase;
			cmp->current = polar(current, phase * M_PI / 180);
		}
		else if (cmp->type == "cccs")
		{
			string ignore;
			InputFile >> cmp->node3 >> cmp->node4 >> ignore >> cmp->coeff;
		}
		else if (cmp->type == "ccvs")
		{
			string ignore;
			InputFile >> cmp->node3 >> cmp->node4 >> ignore >> cmp->coeff;
		}
		else if (cmp->type == "vcvs")
		{
			InputFile >> cmp->node3 >> cmp->node4 >> cmp->coeff;
		}
		else if (cmp->type == "vccs")
		{
			InputFile >> cmp->node3 >> cmp->node4 >> cmp->coeff;
		}
		else if (cmp->type == "res")
		{
			double resistance;
			InputFile >> resistance;
			cmp->Z = complex<double>(resistance, 0);
		}
		else if (cmp->type == "cap")
		{
			double conductance;
			InputFile >> conductance;
			cmp->Z = complex<double>(0, -1 / (W * conductance));

		}
		else if (cmp->type == "ind")
		{
			double inductance;
			InputFile >> inductance;
			cmp->Z = complex<double>(0, W * inductance);
		}
	}
	InputFile.close();


	//Applying STA Algorithm for circuit analysis and computing the matrix 
	int Compcount = Components.size();
	int size = 2 * Compcount + NodeCount;
	MatrixXcd Matrix = Matrix.Zero(size, size);
	VectorXcd Coeff = Coeff.Zero(size);

	for (int i = 0; i < Compcount; i++)
	{
		int n1 = Components[i]->node1;
		int n2 = Components[i]->node2;
		string type = Components[i]->type;
		int index = NodeCount + Compcount + i;
		int PassiveIndex = WhichPassive(Components[i]->node3, Components[i]->node4);
		
	
		//KCL Loop & Branch Voltage
		if (n1 != 0)
		{
			Matrix(n1 - 1, i) = 1;
			Matrix(i + NodeCount, n1 + 2 * Compcount - 1) = -1;
		}
		if (n2 != 0)
		{
			Matrix(n2 - 1, i) = -1;
			Matrix(i + NodeCount, n2 + 2 * Compcount - 1) = 1;
		}

		Matrix(i + NodeCount, i + Compcount) = 1; //1's in main diagonal in the middle section

		Matrix(index, i + Compcount) = 1; //1's in main diagonal in the lower section

		if (type == "isrc" || type == "cccs" || type == "vccs")
		{
			Matrix(index, i) = 1;
			Matrix(index, i + Compcount) = 0; //excluding current sources from lower 1's
		}

		//General cases
		if (type == "res" || type == "cap" || type == "ind")
			Matrix(index, i) = -(Components[i]->Z);
		else if (type == "vsrc")
			Coeff(index, 1) = Components[i]->voltage;
		else if (type == "isrc")
			Coeff(index, 1) = -(Components[i]->current);
		else if (type == "cccs")
			Matrix(index, PassiveIndex) = Components[i]->coeff;
		else if (type == "vcvs")
			Matrix(index, PassiveIndex + Compcount) = -(Components[i]->coeff);
		else if (type == "ccvs")
			Matrix(index, PassiveIndex) = -(Components[i]->coeff);
		else if (type == "vccs")
			Matrix(index, PassiveIndex + Compcount) = Components[i]->coeff;
	}

	//Calculating the output
	VectorXcd Result(size);
	Result = Matrix.colPivHouseholderQr().solve(Coeff);

	//Saving the output in a file
	cout << "Please enter the OUTPUT file name: ";
	cin >> FileName;
	ofstream OutputFile;
	path = "../Circuits/" + FileName + ".txt";
	OutputFile.open(path);
	if (!OutputFile.is_open())
	{
		cout << "File couldn't be created";
		return 2;
	}

	OutputFile << "\nComplex Representation		" << "Phasor representation		" <<"Time Domain Representation\n\n\n";
	OutputFile << "Branch Currents: \n";

	OutputFile << fixed << setprecision(4) << endl;
	cout << W;
	for (int i = 0; i < Compcount; i++)
	{
		int index = i;
		double theta = arg(Result(index)) * 180 / M_PI;
		double magnitude = abs(Result(index));
		OutputFile << "I" << i <<" = " << Result(index) << "	|	" << magnitude << " < " << theta << "	|	" << magnitude << "Cos(" << W << "t + " << theta << ")\n";
	}

	OutputFile << "\nBranch Voltages: \n";
	for (int i = 0; i < Compcount; i++)
	{
		int index = i + Compcount;
		double theta = arg(Result(index)) * 180 / M_PI;
		double magnitude = abs(Result(index));
		OutputFile << "v" << i << " = " << Result(index) << "	|	" << magnitude << " < " << theta << "	|	" << magnitude << "Cos(" << W << "t + " << theta << ")\n";
	}

	OutputFile << "\nNodes Voltages: \n";
	for (int i = 0; i < NodeCount; i++)
	{
		int index = i + 2 * Compcount;
		double theta = arg(Result(index)) * 180 / M_PI;
		double magnitude = abs(Result(index));
		OutputFile << "V" << i << " = " << Result(index) << "	|	" << magnitude << " < " << theta << "	|	" << magnitude << "Cos(" << W << "t + " << theta << ")\n";
	}
	OutputFile.close();
}

 