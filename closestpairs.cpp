// Description:  Finding the closest pair of points in O(n^2) and O(nlogn) time

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <iterator>
#include <chrono>

using namespace std;

// Structure for holding pair of x and y values
struct Point
{
	int x, y, index;
	Point(int xval, int yval, int line_no)
	{
		x = xval;
		y = yval;
		index = line_no;
	}
};

// Global Declarations
int numPoints, DEBUG, point1, point2, point3, point4;
float time1, time2, dist1, dist2, finish1, finish2;
vector<Point> table;
vector<int> points1, points2;
ofstream output;

// Forward declarations
void inputData();
bool compareX(const Point &a, const Point &b);
bool compareY(const Point &a, const Point &b);
float compareDistance(float x, float y);
float dist(Point a, Point b);
float brute(vector<Point> v);
float optimized(vector<Point> x, vector<Point> y);
vector<int> returnPoints(int p1, int p2);
void outputData(vector<int> p1, vector<int> p2, float d1, float d2, float t1, float t2);

// Reading in points from "data"
void inputData()
{
	int x, y;
	int line_no = 1;
	std::ifstream data;
	data.open("data");
	// File can't open
	if (!data.is_open())
	{
		std::exit(EXIT_FAILURE);
	}
	data >> numPoints >> DEBUG;
	while (data >> x >> y)
	{
		table.push_back(Point(x, y, line_no));
		line_no++;
	}
	data.close();
}

// Functions for sorting and calculating distance
bool compareX(const Point &a, const Point &b)
{
	return a.x < b.x;
}

bool compareY(const Point &a, const Point &b)
{
	return a.y < b.y;
}

float compareDistance(float x, float y)
{
	return (x < y) ? x : y;
}

float dist(Point a, Point b)
{
	return sqrt(((b.x - a.x) * (b.x - a.x)) + ((b.y - a.y) * (b.y - a.y)));
}

// Function for calculating closest pairs in O(n^2)
float brute(vector<Point> v)
{
	float min = numeric_limits<float>::max();

	for (int i = 0; i < v.size() - 1; i++)
	{
		for (int j = i + 1; j < v.size(); j++)
		{
			if (dist(v[i], v[j]) < min)
			{
				min = dist(v[i], v[j]);
				point1 = v[i].index - 1;
				point2 = v[j].index - 1;
			}
		}
	}
	return min;
}

// Functions for calculating closest pairs in O(nlogn)
float optimized(vector<Point> x, vector<Point> y)
{
	if (DEBUG == 1)
	{
		output << "pX(Lx) is :" << endl;
		for (int i = 0; i < x.size(); i++)
			output << x[i].index << ": " << x[i].x << " " << x[i].y << endl;
		output << endl;
		output << "pY(Ly) is :" << endl;
		for (int i = 0; i < y.size(); i++)
			output << y[i].index << ": " << y[i].x << " " << y[i].y << endl;
		output << endl;
	}

	if (x.size() <= 3)
	{
		return brute(x);
	}

	int N = x.size();
	vector<Point> xL, xR;
	xL.reserve(N / 2 + 1);
	xR.reserve(N / 2 + 1);

	copy(x.begin(), x.begin() + (N / 2), back_inserter(xL));
	copy(x.begin() + (N / 2), x.end(), back_inserter(xR));

	int xM = x[N / 2].x;

	vector<Point> yL, yR;
	yL.reserve(N / 2 + 1);
	yR.reserve(N / 2 + 1);

	copy_if(y.begin(), y.end(), back_inserter(yL), [&xM](const Point& a)
	{return a.x <= xM; });
	copy_if(y.begin(), y.end(), back_inserter(yR), [&xM](const Point& a)
	{return a.x > xM; });

	float deltaLeft = optimized(xL, yL);
	float deltaRight = optimized(xR, yR);
	float delta = compareDistance(deltaLeft, deltaRight);

	if (DEBUG == 1)
	{
		output << "delta: " << setprecision(6) << fixed << delta << endl;
		output << endl;
	}
	vector<Point> strip;
	strip.reserve(y.size());

	for (int i = 0; i < y.size(); i++)
	{
		if (abs(y[i].x - xM) < delta)
			strip.push_back(y[i]);
	}

	if (DEBUG == 1)
	{
		output << "strip is: " << endl;
		for (int i = 0; i < strip.size(); i++)
			output << strip[i].index << ": " << strip[i].x << " " << strip[i].y << endl;

		output << endl;
	}

	for (int i = 0; i < strip.size(); i++)
	{
		for (int j = i + 1; j < strip.size() && (strip[j].y - strip[i].y) < delta; j++)
		{
			float temp = dist(strip[i], strip[j]);

			if (temp < delta)
			{
				delta = temp;
				if (DEBUG == 1)
				{
					output << "min_dst in strip: " << setprecision(6) << fixed << delta << endl;
					output << endl;
				}
			}
		}
	}
	return delta;
}

// Calculating closest pair of points
vector<int> returnPoints(int p1, int p2)
{
	vector<int> pairofpoints;
	pairofpoints.reserve(4);
	int x1, x2, y1, y2;
	for (int i = 0; i < numPoints; i++)
	{
		if (table[i].index - 1 == p1)
		{
			x1 = table[i].x;
			y1 = table[i].y;
		}
		if (table[i].index - 1 == p2)
		{
			x2 = table[i].x;
			y2 = table[i].y;
		}
	}
	if (x1 < x2)
	{
		pairofpoints.push_back(x1);
		pairofpoints.push_back(y1);
		pairofpoints.push_back(x2);
		pairofpoints.push_back(y2);
	}
	else
	{
		pairofpoints.push_back(x2);
		pairofpoints.push_back(y2);
		pairofpoints.push_back(x1);
		pairofpoints.push_back(y1);
	}
	return pairofpoints;
}

void nsquare()
{
	auto start = chrono::high_resolution_clock::now();
	dist1 = (brute(table));
	points1 = returnPoints(point1, point2);
	auto stop = chrono::high_resolution_clock::now();
	time1 = chrono::duration_cast<chrono::microseconds> (stop - start).count();
}

void nlogn()
{
	vector<Point> x, y;
	x = table;
	x.reserve(table.size());
	y = table;
	y.reserve(table.size());
	auto start1 = chrono::high_resolution_clock::now();
	sort(x.begin(), x.end(), compareX);
	sort(y.begin(), y.end(), compareY);
	dist2 = (optimized(x, y));
	points2 = returnPoints(point1, point2);
	auto stop1 = chrono::high_resolution_clock::now();
	time2 = chrono::duration_cast<chrono::microseconds> (stop1 - start1).count();
}

// Output the data
void outputData(vector<int> p1, vector<int> p2, float d1, float d2, float t1, float t2)
{
	float temp1, temp2;
	temp1 = t1;
	temp2 = t2;
	output << "Closest pair by O(nlog(n)) algorithm: ";
	output << "(" << p1[0] << ", " << p1[1] << ") (" << p1[2] << ", " << p1[3] << ")" << endl;
	output << "Distance: ";
	output << setprecision(6) << fixed << d2 << endl;
	output << "Time\t";
	output << setprecision(0) << fixed << temp2 << " ms. " << endl;
	output << endl;
	output << "n: " << numPoints << endl;
	output << "time2/(n*log(n)): ";
	output << setprecision(6) << fixed << t2 / (numPoints * log2(numPoints)) << endl;
	output << "time2/(n*log(n)^2): ";
	output << setprecision(6) << fixed << t2 / (numPoints * pow(log2(numPoints), 2));
}

int main()
{
	output.open("output");
	// Read from "data"
	inputData();

	// Closest pairs in O(n^2)
	nsquare();
	output << "Closest pair by O(n^2) algorithm: ";
	output << "(" << points1[0] << ", " << points1[1] << ") (" << points1[2] << ", " << points1[3] << ")" << endl;
	output << "Distance: ";
	output << setprecision(6) << fixed << dist1 << endl;
	output << "Time\t";
	output << setprecision(0) << fixed << time1 << " ms. " << endl;
	output << endl;

	// Closest pairs in O(nlogn)
	nlogn();

	// Output to "output"
	outputData(points1, points2, dist1, dist2, time1, time2);
	output.close();

	return 0;
}
