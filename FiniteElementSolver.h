#ifndef FINITEELEMENT_H
#define FINITEELEMENT_H

#include <math.h>
#include <armadillo>
#include <iostream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
    


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                       Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>                    Delaunay;
typedef Kernel::Point_2                                                Point;


class FiniteElementSolver
/*Class for solving the 2D Poisson equation Lu = f.
 * Constructor takes a pair of points and indices, which specifie the domain
 * and the global index of the vertex associated with this point.
 *
 * prescribed values is a pair of indices and values, which give the dirichlet
 * boundary conditions.
 *
 * rightHandF holds the function of the right hand side of the PDE.
 *
 * Initating the class does not call any function.
 *
 * Genral process of solving the equation looks like this:
 *  1. call calculateGlobalStiffnessMatrix()
 *  2. call calculateRightHandSide()
 *  3. solveSystem()
 *
 */
{
    private:
        int numberDof;
        Delaunay triangulation;
        arma::mat elementaryBasisMatrix;
        arma::mat PDEMatrix;
        arma::vec rightHandSide;
        arma::sp_mat globalStiffnessMatrix;
        arma::mat gradBasis;
        std::vector<std::pair<double,unsigned> > prescribedValuePair;
        arma::vec solution;
        std::function<double (double,double)> rightHandFunction;
        
        void calculateTransform( Delaunay::Face_handle face,arma::mat& transformMatrix);
        arma::mat calculateElementStiffnessMatrix(Delaunay::Face_handle element);
        arma::mat calculateElementRightHandSide( Delaunay::Face_handle element );

    public:
        
        FiniteElementSolver(std::vector <std::pair<Point,unsigned> >& points,std::vector<std::pair<double,unsigned> >& prescribed,std::function<double (double,double)>& rightHandF);   
        void calculateGlobalStiffnessMatrix();       
        void calculateRightHandSide();
        void solveSystem();
        arma::vec getSolution();

};

#endif
