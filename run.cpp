#include "FiniteElementSolver.h"
#include <iostream>
#include <string>
#include <sstream>

using namespace std;
using namespace arma;
int main()
{
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
    cout << "You entered: " << n << endl << endl;
    
    mat PDEMatrix(2,2,fill::eye);
    vector< pair<Point,unsigned> > points;  
    vector<pair<double,unsigned> > prescribed;
    vec x_Axis = linspace<vec> (0,1,n);
    vec y_Axis = linspace<vec> (0,1,n);
    unsigned int counter = 0;
    
    function<double (double,double)> sol = [](double x, double y)
        {
            return exp(x+0.2*y);
        };
    function<double (double,double)> rightFunction = [](double x, double y)
        {
            return -1.04 * exp(x+0.2*y);
        };

    for(const auto& valx : x_Axis)
    {
        for(const auto& valy :y_Axis)
        {
            points.push_back(make_pair(Point(valx,valy),counter));

            if(valx == 0 || valy == 0 || valx==1 || valy==1 )
            {
                prescribed.push_back(make_pair(sol(valx,valy),counter));
                cout<<sol(valx,valy)<<endl;
            }
            counter++; 
        }
    }   

    FiniteElementSolver solve(points,prescribed,rightFunction);

    solve.calculateGlobalStiffnessMatrix();
    solve.calculateRightHandSide();
    solve.solveSystem();
    return 0;  

}

