#include "FiniteElementSolver.h"

using namespace arma;
using namespace std;

FiniteElementSolver::FiniteElementSolver(
        vector<pair<Point,unsigned> >& points,
        vector<pair <double,unsigned> > prescribed
        )
{
    prescribedValuePair = prescribed;
    
    triangulation.insert(points.begin(),points.end());
    PDEMatrix= mat(2,2,fill::eye);
    
    numberDof=points.size();
    rightHandSide = vec(numberDof,fill::zeros);
    solution= vec(numberDof,fill::zeros);
    globalStiffnessMatrix = mat(numberDof,numberDof,fill::zeros);
    gradBasis<< -1<< 1<< 0<<endr
             << -1<< 0<< 1<<endr; 

};  

void FiniteElementSolver::calculateTransform(   Delaunay::Face_handle face ,
                                                mat& transformMatrix)
{
    mat referenceCoord(3,3); 

    referenceCoord  <<0<<1<<0<<endr
                    <<0<<0<<1<<endr
                    <<1<<1<<1<<endr;
    //get transformmatrix 
    for (int i =0 ; i<3 ; ++i)
    {
        vec column(3,fill::zeros);
        Point point = face->vertex(i)->point();
        column(0) = point.x(); 
        column(1) = point.y();
        column(2) = 1;
        transformMatrix.col(i) = column;
    }
    transformMatrix = transformMatrix * referenceCoord.i();
};

mat FiniteElementSolver::calculateElementStiffnessMatrix(Delaunay::Face_handle element)
{
    mat elementStiffnessMatrix(3,3);
    mat transformMatrix(3,3);
    calculateTransform(element,transformMatrix);
    double determinant = abs(det(transformMatrix(span(0,1),span(0,1))));
    for (int row =0 ; row <3 ; ++row)
    {
        for(int column=0 ; column<3 ; ++column)

        { 
            vec left(transformMatrix(span(0,1),span(0,1)).i().t()*gradBasis.col(column));
            vec right(PDEMatrix * transformMatrix(span(0,1),span(0,1)).i().t()* gradBasis.col(row));
            elementStiffnessMatrix(row,column) =0.5*determinant* dot(left,right);
        }
    }
    return elementStiffnessMatrix;
};

void FiniteElementSolver::calculateGlobalStiffnessMatrix()
{
    vector<unsigned long long> RowIndices;
    vector<unsigned long long> ColumnIndices;    
    vector<double>   ValueVector;
      
     
    for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); 
        fit != triangulation.finite_faces_end(); ++fit) 
    {
        Delaunay::Face_handle face = fit;
        mat elementStiffness(3,3);
        elementStiffness = calculateElementStiffnessMatrix(face);
        
        for (int row=0; row<3;++row)
        {
            unsigned int globalRowIndex = face->vertex(row)->info();

            for (int column=0; column<3; ++column)
            {
                unsigned int globalColumnIndex = face->vertex(column)->info();
                RowIndices.push_back(globalRowIndex);
                ColumnIndices.push_back(globalColumnIndex);
                ValueVector.push_back(elementStiffness(row,column));
            }
        }
    };
    unsigned int LengthOfVector = RowIndices.size();  
    umat locations(2,LengthOfVector);
    locations.row(0) = urowvec(RowIndices);
    locations.row(1) = urowvec(ColumnIndices);

    globalStiffnessMatrix = sp_mat(true,locations,vec(ValueVector),numberDof,numberDof);
    
    for(const pair<double,unsigned>& pa: prescribedValuePair)
    {
      rowvec insertRow(numberDof,fill::zeros);
      insertRow(pa.second) = 1.0;
      globalStiffnessMatrix.row(pa.second)=insertRow;
    }
    //mat(globalStiffnessMatrix).print();
};       

void FiniteElementSolver::calculateRightHandSide()
{
    for(const pair<double,unsigned>& pa:prescribedValuePair)
    {
        rightHandSide(pa.second) = pa.first;
    }
    //rightHandSide.print();
};
void FiniteElementSolver::solveSystem()
    {
        solution = spsolve(globalStiffnessMatrix,rightHandSide); 
        solution.print();
    };

