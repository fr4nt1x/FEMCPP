#ifndef FINITEELEMENT_H
#define FINITEELEMENT_H
#include <math.h>
#include <armadillo>
#include <iostream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
    

//you Need to link agains -lgmp

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                       Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>                    Delaunay;
typedef Kernel::Point_2                                                Point;


class FiniteElementSolver
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

    public:
        FiniteElementSolver(std::vector <std::pair<Point,unsigned> >& points,std::vector<std::pair<double,unsigned> > prescribed);   
        void calculateTransform( Delaunay::Face_handle face,arma::mat& transformMatrix);
        arma::mat calculateElementStiffnessMatrix(Delaunay::Face_handle element);
        void calculateGlobalStiffnessMatrix();       
        void calculateRightHandSide();
        void solveSystem();

};

#endif
