#include <nanoflann.hpp>
using namespace nanoflann;

#include "KDTreeVectorOfVectorsAdaptor.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

typedef std::vector<std::vector<double> > my_vector_of_vectors_t;

class kdtree{
    private:
        std::vector<size_t> ind_mic;
        std::vector<std::vector<double> > positions;
        std::vector<double> cell;

    public:
        std::vector<std::vector<size_t> > indices;
        std::vector<std::vector<double> > distances;
        size_t num_neighbors;

        kdtree(int num_neighbors_in) : num_neighbors(num_neighbors_in){}

        void set_cell(double *cell_in)
        {
            cell.resize(3);
            for(size_t dim=0; dim<3; dim++)
                cell[dim] = cell_in[dim];
        }

        double get_distance(size_t ind, size_t n_neigh){
            if(ind>distances.size())
                throw std::invalid_argument("Positions not defined or index invalid");
            if(n_neigh>num_neighbors)
                throw std::invalid_argument("Neighbor index invalid");
            if(distances[ind].size()==0)
                run();
            return sqrt(distances[ind][n_neigh+1]);
        }

        int get_index(size_t ind, size_t n_neigh){
            if(indices.size()==0)
                throw std::invalid_argument("Positions not defined");
            if(n_neigh>num_neighbors)
                throw std::invalid_argument("Neighbor index invalid");
            if(indices[ind].size()==0)
                run();
            return ind_mic[indices[ind][n_neigh+1]];
        }

        void run(){
            nanoflann::KNNResultSet<double> resultSet(num_neighbors+1);
	        KDTreeVectorOfVectorsAdaptor< my_vector_of_vectors_t, double > mat_index(3 /*dim*/, positions, 10 /* max leaf */ );
            mat_index.index->buildIndex();
            for(size_t ind=0; ind<indices.size(); ind++)
            {
                indices[ind].resize(num_neighbors+1);
                distances[ind].resize(num_neighbors+1);
                resultSet.init(&(indices[ind][0]), &(distances[ind][0]));
                mat_index.index->findNeighbors(resultSet, &(positions[ind][0]), nanoflann::SearchParams(10));
            }
        }

        void set_positions(double * positions_in, int number_of_atoms, double skin=10){
            if(cell.size()!=3)
                throw std::invalid_argument("Cell not defined");
            indices.resize(number_of_atoms);
            distances.resize(number_of_atoms);
            ind_mic.resize(number_of_atoms);
            positions.resize(number_of_atoms);
            std::vector<double> pos_tmp(3);
            int ix[3];
            for(size_t ind=0; ind<number_of_atoms; ind++)
            {
                positions[ind].resize(3);
                for (int dim=0; dim<3; dim++)
                    positions[ind][dim] = positions_in[ind*3+dim];
                ind_mic[ind] = ind;
                for(ix[0]=-1; ix[0]<=1; ix[0]++)
                    for(ix[1]=-1; ix[1]<=1; ix[1]++)
                        for(ix[2]=-1; ix[2]<=1; ix[2]++)
                        {
                            if(abs(ix[0])+abs(ix[1])+abs(ix[2])==0)
                                continue;
                            for(int dim=0; dim<3; dim++)
                                pos_tmp[dim] = positions_in[ind*3+dim]+cell[dim]*ix[dim];
                            if(get_dist_from_box(pos_tmp.data())<skin)
                            {
                                // std::cout<<"hit"<<std::endl;
                                positions.push_back(pos_tmp);
                                ind_mic.push_back(ind);
                            }
                        }
            }
            std::cout<<"Tree initialization finished; number of points: "<<positions.size()<<std::endl;
        }

        double get_dist_from_box(double *point){
            double length[3];
            for(int ii=0; ii<3; ii++)
                length[ii] = -point[ii]*(point[ii]<0)+(point[ii]-cell[ii])*((point[ii]-cell[ii])>0);
            return *std::max_element(length, length+3); 
        }
};

