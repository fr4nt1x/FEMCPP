#include "FiniteElementSolver.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;
using namespace arma;

int main()
{
    //get the user input of how many points the user 
    //wants to use in each direction     
    int n;
    string input = "";
    
    while (true)
    {
        cout << "How many points do you want to use in x and y direction?";
        getline(cin, input);

        // This code converts from string to number safely.
        stringstream myStream(input);
        if (myStream >> n)
            break;
        cout << "Invalid number, please try again" << endl;
    }

    //vetor with the point index,pairs
    vector< pair<Point,unsigned> > points;

    //holds Vector of dirichlet boundary condition  
    vector<pair<double,unsigned> > prescribed;
    
    //exact solution, used for calculating the boundary values
    function<double (double,double)> sol = [](double x, double y)
        {
            return 1+x*x +2*y*y;
        };

    //right Hand side given to the class instance for calculating the Right
    //Hand Side
    function<double (double,double)> rightFunction = [](double x, double y)
        {
            return -6;
        };
    double boundXMin,boundXMax,boundYMin,boundYMax;
    boundXMin = 0;
    boundYMin = 0;
    boundXMax = 1;
    boundYMax = 1;

    //Use an equal spacing of x and y of the box
    //[boundXMin,boundXMin]x[boundYMin,boundYMax]
    vec x_Axis = linspace<vec> (boundXMin,boundXMax,n);
    vec y_Axis = linspace<vec> (boundYMin,boundYMax,n);
    unsigned int counter = 0;
    
    for(const auto& valx : x_Axis)
    {
        for(const auto& valy :y_Axis)
        {
            points.push_back(make_pair(Point(valx,valy),counter));
            //add the index to the boundary condition, if it is at left, right,
            //top or bottom of the domain
            if(valx == boundXMin || valy == boundYMin || valx== boundXMax || valy== boundYMax )
            {
                prescribed.push_back(make_pair(sol(valx,valy),counter));
            }
            counter++; 
        }
    }   

    
    FiniteElementSolver solve(points,prescribed,rightFunction);

    solve.calculateGlobalStiffnessMatrix();
    solve.calculateRightHandSide();
    solve.solveSystem();
    
    //put solution and points in file for plotting in python
    vec solh = solve.getSolution();
    
    
    ofstream solutionFile;
    solutionFile.open("solution",ios::trunc);

    int i= 0;
    vector< pair<Point,unsigned> >::iterator pointIt;
    
    for (pointIt= points.begin() ; pointIt < points.end();pointIt++,i++)
    {
        pair<Point,unsigned> PointPair= *pointIt ;
        solutionFile<< setprecision(64)<<PointPair.first<<" "<<solh(i)<<endl;
    }
    solutionFile.close();
    
    return 0;  

}

