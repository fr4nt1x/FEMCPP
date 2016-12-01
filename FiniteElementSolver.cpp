#include "FiniteElementSolver.h"

using namespace arma;
using namespace std;

FiniteElementSolver::FiniteElementSolver(
        vector<pair<Point,unsigned> >& points,
        vector<pair <double,unsigned> >& prescribed,
        std::function<double (double,double)>& rightHandF
        )
    /*
     *  Class for solving the 2D elliptic PDE with Dirichlet boundary condition
     *  Uses linear lagrangian finite elements(basis functions are called /phi_i)
     *  on an triangular grid.
     *  
     *  The grid is calculated via the Delaunaytriangulation from the CGAL Package
     *  
     *
     *  PDE :   div( PDE_Matrix * grad u) = f
     *
     * INPUT :  points      hold a pair with a 2D coordinates and the global
     *                      index of the node
     *
     *          prescribed  hold Dirichlet boundary conditions, as a pair
     *                      of values and the global index
     *
     *          rightHandF  the std function container of the Function f used for
     *                      the RHS
     */

{
    //Store the Boundary conditions
    prescribedValuePair = prescribed;

    triangulation.insert(points.begin(),points.end());
    
    //Store the Matrix of the PDE 
    PDEMatrix= mat(2,2,fill::eye);

    //Store the Right hand funtion
    rightHandFunction = rightHandF;

    //hold integral over the standard element of /phi_i * /phi_j
    //used for calculating the RightHandside without Numerical Integration rule
    elementaryBasisMatrix ={{1,0.5,0.5},
                            {0.5,1,0.5},
                            {0.5,0.5,1}};
    elementaryBasisMatrix = elementaryBasisMatrix/12.0; 
    
    //store the number of degrees of freedom, for initiating the stiffness
    //matrix and the righthandside
    numberDof=points.size();

    //initiate the righthandside, the solution and the global stiffness matrix
    //with zeros
    rightHandSide  = vec(numberDof,fill::zeros);
    solution= vec(numberDof,fill::zeros);
    
    //Gradients of the linear basis functions, used to calculate the Stiffness
    //matrices
    gradBasis<< -1<< 1<< 0<<endr
             << -1<< 0<< 1<<endr; 

};  

void FiniteElementSolver::calculateTransform(   Delaunay::Face_handle face ,
                                                mat& transformMatrix)
    /*Calculate the affine transformation Ax +b from the reference element to the
     * given element. A is a 2x2 matrix and b vector.
     *
     * face                 holds the handle of the actual element
     *
     * transformMatrix      a reference where the result is stored
     *                      the upper left submatrix is A, the the last column
     *                      wihtout the last entry holds b
     */
{
    //initiate the the Points of the reference element in matrix form, with
    //ones appended to the last row
    mat referenceCoord(3,3); 
    referenceCoord  <<0<<1<<0<<endr
                    <<0<<0<<1<<endr
                    <<1<<1<<1<<endr;

    //Do the same as for the reference element with the actual element points.
    for (int i =0 ; i<3 ; ++i)
    {
        vec column(3,fill::zeros);
        Point point = face->vertex(i)->point();
        column(0) = point.x(); 
        column(1) = point.y();
        column(2) = 1;
        transformMatrix.col(i) = column;
    }

    //Calculated the  transform matrix
    transformMatrix = transformMatrix * referenceCoord.i();
};

mat FiniteElementSolver::calculateElementStiffnessMatrix(Delaunay::Face_handle element)
/*Calculated the element stiffness matrix for the given element.
 * First calculate the transform matrix and vector.
 *
 */
{

    mat elementStiffnessMatrix(3,3);
    mat transformMatrix(3,3);
    calculateTransform(element,transformMatrix);
    
    //the integral is calculated via the transformation formula so the
    //determinant of the Jacobian is needed
    double determinant = abs(det(transformMatrix(span(0,1),span(0,1))));
    
    //the i,j -th entry of the element stiffness matrix comes from the weak
    //form of the PDE:: 
    //a(/phi_i,/phi_j) = integral_ele PDEMatrix * grad /phi_i dot /grad phi_j
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
    /*  Calculate the global stiffness matrix.
     *  Calls calculateElementStiffnessMatrix for each element and assembles
     *  the entries into the global frame.
     *  The stiffness matrix is a sparse matrix which is initiated via
     *  coordinate form.
     */
{
    //Vectors for the initiating in coordinate form
    vector<unsigned long long> RowIndices;
    vector<unsigned long long> ColumnIndices;    
    vector<double>   ValueVector;
      
    //Loop over all elements and call the element function
    for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); 
        fit != triangulation.finite_faces_end(); ++fit) 
    {
        Delaunay::Face_handle face = fit;
        mat elementStiffness(3,3);
        elementStiffness = calculateElementStiffnessMatrix(face);
        
        //assemble the element stiffness matrix into a global index frame
        for (int row=0; row<3;++row)
        {
            unsigned int globalRowIndex = face->vertex(row)->info();
            for (int column=0; column<3; ++column)
            {
                //info holds the global index of the given vertex
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

    //Change the rows for which dirichlet boundary conditions are prescribed into
    //a row with only a one at the global index. The Right hand side is later
    //put to only have the given value at the according global index.    
    for(const pair<double,unsigned>& pa: prescribedValuePair)
    {
      rowvec insertRow(numberDof,fill::zeros);
      insertRow(pa.second) = 1.0;
      globalStiffnessMatrix.row(pa.second)=insertRow;
    }
    //mat(globalStiffnessMatrix).print();
};       

mat FiniteElementSolver::calculateElementRightHandSide( Delaunay::Face_handle element)
    /*Calculate the element RightHandSide for the given element.
     * Uses the approximates the function on the right hand side with the same
     * linear lagrangian elements.
     *
     * The value of the integral can then be written via the elematry basis
     * matrix, as product of this matrix and the function evaluated at the vertices of
     * the element (assembled in a vector).
     *
     * Again the transformation formula is used. So only the reference basis
     * function are needed.
     */
{
    vec elementRightHandSide(3,fill::zeros);
    mat transformMatrix(3,3);
    calculateTransform(element,transformMatrix);
    double determinant = abs(det(transformMatrix(span(0,1),span(0,1))));
    
    
    vec functionEval(3,fill::zeros);
    for (int i = 0; i<3 ; ++i)
    {
        Point evalpoint = element->vertex(i)->point();
        functionEval(i) = rightHandFunction(evalpoint.x(),evalpoint.y()); 
        //cout<<evalpoint<<"\t"<<element->vertex(i)->info()<<endl;
    }

    elementRightHandSide = determinant * (elementaryBasisMatrix * functionEval);
    //cout<<"ERHS";
    //elementRightHandSide.print();
    
    return elementRightHandSide;
};


void FiniteElementSolver::calculateRightHandSide()
    /*  Loop over all elements, call calculateElementRightHandSide for each
     *  element, and assemble it at the global place
     */
{
    for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); 
        fit != triangulation.finite_faces_end(); ++fit) 
    {
        Delaunay::Face_handle face = fit;
        vec eleRightHandSide = calculateElementRightHandSide(face);
        for (int i =0; i<3; ++i)
        {
            rightHandSide(face->vertex(i)->info())+=eleRightHandSide(i);
        }
    }

    
    //The vertices at which values are prescribed get replaced with the
    //prescribed values
    for(const pair<double,unsigned>& pa:prescribedValuePair)
    {
        rightHandSide(pa.second) = pa.first;
    }

    //cout<<"RightHandSide:"<<endl;
    //rightHandSide.print();
};

void FiniteElementSolver::solveSystem()

    /* Solve the System K u = RHS.
     * Only working after global Stiffness and RightHandside have been
     * assembled.
     */

{
    solution = spsolve(globalStiffnessMatrix,rightHandSide); 
    //cout<<"solution:"<<endl;
    //solution.print();
};

vec FiniteElementSolver::getSolution()
{
    return solution;
}
