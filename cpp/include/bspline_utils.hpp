/*
* bspline_utils.hpp
*
* ---------------------------------------------------------------------
* Copyright (C) 2022 Matthew (matthewoots at gmail.com)
*
*  This program is free software; you can redistribute it and/or
*  modify it under the terms of the GNU General Public License
*  as published by the Free Software Foundation; either version 2
*  of the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
* ---------------------------------------------------------------------
*
* 
* 
*/

#ifndef BSPLINE_UTILS_H
#define BSPLINE_UTILS_H

#include <string>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

namespace trajectory
{
    class bspline_trajectory
    {
        
        /*
        * General matrix representations for B-splines See Theorem 1 of page 182 for getting M
        * https://link.springer.com/article/10.1007/s003710050206
        */

        public:

        struct bs_utils
        {
            int k; // Order + 1
            int k_d; // Number of knot division
            double dt; // dt is the length/span of 1 knot
            Eigen::MatrixXd M; // M matrix that is general
            vector<int> r; // Vector of control points that we have to consider
        };

        struct bs_pva_state
        {
            vector<double> rts; // Relative time span 
            vector<double> pos; // Position vector
            vector<double> vel; // Velocity vector
            vector<double> acc; // Acceleration vector
        };
        
        inline bs_pva_state get_uni_bspline_1d(
            int order, vector<double> timespan, 
            vector<double> ctrlpt, int knotdiv)
        {
            bs_utils b;
            bs_pva_state s;
            b.k = order + 1;
            b.M = create_m(order);
            b.k_d = knotdiv;
            int n = (int)ctrlpt.size() - 1;
            int m = n + order + 1; 
            // int n_k = m + 1; // Number of knots

            // Range of index to evaluate accordingly (order to length of control points)
            b.r = int_range_to_vector(order, m - order - 1); 

            b.dt = (timespan[1] - timespan[0]) / (b.r.size());
            vector<double> t = linspace(timespan[0], timespan[1], b.r.size() + 1);

            // since we know the size of the vector we should reserve first
            int vector_final_size = b.r.size() * (b.k_d-1);
            s.rts.reserve(vector_final_size);
            s.pos.reserve(vector_final_size);
            s.vel.reserve(vector_final_size);
            s.acc.reserve(vector_final_size);

            for (int i = 0; i < (int)b.r.size(); i++)
            {
                // Evaluate from the current idx to next idx
                int idx = b.r[i] - order;
                int nxt_idx = idx + 1;

                // Relative to the start time as 0 regardless of the time 
                vector<double> span = linspace((double)idx, (double)nxt_idx, b.k_d); 
                // Time in abs time (simulation time / given time)
                vector<double> actualspan = linspace(t[i], t[i+1], b.k_d); 

                // Control Points in a Span Column vector
                Eigen::VectorXd p = Eigen::VectorXd::Zero(b.k);
                // Position Row Vector, Velocity Row Vector, Acceleration Row Vector, Snap Row Vector
                Eigen::RowVectorXd u, du, ddu;
                u = du = ddu = Eigen::RowVectorXd::Zero(b.k); 

                // We are only considering [0,1) hence not including 1
                for (int j = 0; j < (int)span.size()-1; j++)
                {
                    // current time in index form, of course we dont like to play with conversion
                    double time = span[j]; 
                    // using index is the same as using time since u_t is a factor
                    double u_t = (time - idx)/((idx+1) - idx); 

                    // Make the u, du, ddu and p matrix
                    for (int l = 0; l < b.k; l++)
                    {
                        u(l) = pow(u_t, l);
                        p(l) = ctrlpt[idx + l];
                        if (l >= 1)
                            du(l) = (l) * pow(u_t, l-1);
                        if (l >= 2)
                            ddu(l) = (l) * (l-1) * pow(u_t, l-2);
                    }

                    s.pos.push_back(
                        position_at_time_segment(b.dt, b.M, u, p));
                    s.vel.push_back(
                        velocity_at_time_segment(b.dt, b.M, du, p));
                    s.acc.push_back(
                        acceleration_at_time_segment(b.dt, b.M, ddu, p));
                    s.rts.push_back(actualspan[j]);
                }
            }

            return s;

        }

        /** @brief Create the range of values within the array of index **/
        inline vector<int> int_range_to_vector(int min, int max)
        {
            vector<int> v;
            for (int i = min; i <= max; i++)
                v.push_back(i);
            
            return v;
        }
        
        /** @brief Creating Bspline recursive M matrix **/
        inline Eigen::MatrixXd create_m(int order)
        {
            int k = order + 1;
            Eigen::MatrixXd M = Eigen::MatrixXd::Zero(k,k);
            double f = 1 / factorial(k - 1);

            for (int i = 0; i < M.rows(); ++i)
            {
                for (int j = 0; j < M.cols(); ++j)
                {
                    double fac = 0;
                    for (int s = j; s <= (k-1) ; s++)
                    {
                        double p21 = (double)(k-s-1);
                        double p22 = (double)(k-i-1);
                        double p11 = (double)(s-j);
                        fac += (pow(-1.0, p11) * get_c(s-j, k) * pow(p21, p22)); 
                    }  
                    double m = get_c(k-i-1, k-1) * fac;
                    m = f * m;
                    M(i, j) = m;            
                }
            }
            return M;
        }

        /** @brief Creating Bspline single element of the M matrix **/
        inline double get_c(int i, int n)
        {
            return factorial(n)/(factorial(i) * factorial(n-i));
        }

        /** @brief Simple factorial function **/
        inline double factorial(int n)
        {
            double factorial = 1.0;
            for (int i = n; i > 1 ; i--)
            {
                factorial = factorial * (double)i;
            }
            return factorial;
        }

        /** @brief linspace function with given min and max **/
        inline vector<double> linspace(
            double min, double max, double n)
        {
            vector<double> linspaced;
            double delta = (max - min) / (n - 1.0);
            linspaced.push_back(min);
            
            for (int i = 1; i < (int)n; i++)
            {
                linspaced.push_back(linspaced[i-1] + delta);
            }

            return linspaced;
        }

        /** @brief 
         * Calculate position value 
         * @param
         * M : order+1 x order+1 matrix
         * u : Row vector of position association vector
         * p : Control points in that segment
         * u * M * p  
        **/
        inline double position_at_time_segment(
            double dt, Eigen::MatrixXd M, Eigen::RowVectorXd u, Eigen::VectorXd p)
        {
            return (u * M * p)(0,0);
        }

        /** @brief 
         * Calculate velocity value
         * @param
         * M : order+1 x order+1 matrix
         * du : Row vector of velocity association vector
         * p : Control points in that segment
         * inv_dt * du * M * p  
        **/
        inline double velocity_at_time_segment(
            double dt, Eigen::MatrixXd M, Eigen::RowVectorXd du, Eigen::VectorXd p)
        {
            return (1/dt) * (du * M * p)(0,0);
        }

        /** @brief 
         * Calculate acceleration value
         * @param
         * M : order+1 x order+1 matrix
         * ddu : Row vector of acceleration association vector
         * p : Control points in that segment
         * pow(inv_dt,2) * ddu * M * p  
        **/
        inline double acceleration_at_time_segment(
            double dt, Eigen::MatrixXd M, Eigen::RowVectorXd ddu, Eigen::VectorXd p)
        {
            return pow((1/dt),2) * (ddu * M * p)(0,0);
        }

    };
}

#endif