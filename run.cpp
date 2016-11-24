#include "FiniteElementSolver.h"
using namespace std;
using namespace arma;

int main()
{
    int n = 3;
    mat PDEMatrix(2,2,fill::eye);
    vector< pair<Point,unsigned> > points;  
    vector<pair<double,unsigned> > prescribed;
    vec x_Axis = linspace<vec> (0,1,n);
    vec y_Axis = linspace<vec> (0,1,n);

    unsigned int counter = 0;
    
    for(const auto& valx : x_Axis)
    {
        for(const auto& valy :y_Axis)
        {
            points.push_back(make_pair(Point(valx,valy),counter));
            if(valx == 0 || valy == 0 || valx==1 || valy==1 )
            {
                prescribed.push_back(make_pair(1,counter));
            }
            counter++; 
        }
    }   
    FiniteElementSolver solve(points,prescribed); solve.calculateGlobalStiffnessMatrix();
    solve.calculateRightHandSide();
    solve.solveSystem();
    return 0;  

}
