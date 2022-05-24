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
    class trajectory_math_fn
    {
        public:

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
    };

    class common_trajectory_tool
    {
        private:

        trajectory_math_fn tmf;

        public:

        /** @brief Uniform Distribution */
        inline vector<Eigen::Vector3d> uniform_distribution(
            Eigen::Vector3d start, vector<Eigen::Vector3d> wp, 
            double max_vel, double knot_span)
        {
            vector<Eigen::Vector3d> keypoints;
            keypoints.push_back(start);
            for (int j = 0; j < (int)wp.size(); j++)
                keypoints.push_back(wp[j]);
            
            vector<double> diff;
            vector<double> segment;
            double total_dist = 0;
            // Until cp_tmp.cols() - order - 1 since that is the last change in a clamped spline
            // But will be different for non-clamped splines 
            for (int j = 0; j < (int)wp.size(); j++)
            {
                Eigen::Vector3d diff_tmp = keypoints[j+1] - keypoints[j];
                diff.push_back(
                    sqrt(pow(diff_tmp[0],2) + 
                    pow(diff_tmp[1],2) + 
                    pow(diff_tmp[2],2)));
                total_dist =+ diff[j];
            }

            double est_dist_knot = max_vel * knot_span;

            vector<double> time_waypoint;
            double total_segment = 0;
            for (int j = 0; j < (int)diff.size(); j++)
            {
                // ceil helps to push values above 0 to 1 or more
                // or else segment count is 0 and causes an error
                double knot_in_seg_count = ceil(diff[j]/est_dist_knot);
                total_segment += knot_in_seg_count;
                
                segment.push_back(knot_in_seg_count);
                time_waypoint.push_back(total_segment * knot_span);
            }
            
            // Raw control points
            // total_segment + 1 is the total count since + 1 refers to the last missing point after linspace
            vector<Vector3d> cp;
            cp.reserve(total_segment + 1);
            for (int j = 0; j < total_segment + 1; j++)
                cp.push_back(Eigen::Vector3d::Zero());

            for (int i = 0; i < 3; i++)
            {
                double segment_counter = 0;
                for (int j = 0; j < (int)wp.size(); j++)
                {
                    double div = segment[j] + 1.0;
                    
                    vector<double> sub_cp_tmp = 
                        tmf.linspace(keypoints[j](i), keypoints[j+1](i), div);

                    for (int k = 0; k < segment[j]; k++)
                    {
                        cp[k + segment_counter](i) = sub_cp_tmp[k];
                    }
                    segment_counter += segment[j];
                }
                cp[segment_counter](i) = 
                    keypoints[keypoints.size()-1](i);
            }

            return cp;
        }

    };

    class bspline_trajectory
    {
        private:

        trajectory_math_fn tmf;
        
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

        struct bs_pva_state_1d
        {
            vector<double> rts; // Relative time span 
            vector<double> pos; // Position vector 1d
            vector<double> vel; // Velocity vector 1d
            vector<double> acc; // Acceleration vector 1d
        };

        struct bs_pva_state_3d
        {
            vector<double> rts; // Relative time span 
            vector<Eigen::Vector3d> pos; // Position vector 3d
            vector<Eigen::Vector3d> vel; // Velocity vector 3d
            vector<Eigen::Vector3d> acc; // Acceleration vector 3d
        };

        inline bs_pva_state_3d get_uni_bspline_3d(
            int order, vector<double> timespan, 
            vector<Vector3d> ctrlpt3, int knotdiv)
        {
            // Assign the control point data into row vectors
            vector<vector<double>> ctrlpt;
            for (int i = 0; i < 3; i++)
            {
                vector<double> axis_seperation;
                for (int j = 0; j < (int)ctrlpt3.size(); j++)
                {
                    axis_seperation.push_back(ctrlpt3[j](i));
                }
                ctrlpt.push_back(axis_seperation);
            }

            vector<bs_pva_state_1d> s1;
            // Pass the row vectors into the 1d bspline creation
            for (int i = 0; i < 3; i++)
                s1.push_back(get_uni_bspline_1d(order, timespan, ctrlpt[i], knotdiv));

            // Reassemble from row vectors into Vector3d columns
            bs_pva_state_3d s3;
            s3.rts = s1[0].rts;
            int total_array_size = s1[0].rts.size();
            for (int i = 0; i < total_array_size; i++)
            {
                Eigen::Vector3d pos_vect = Eigen::Vector3d::Zero();
                Eigen::Vector3d vel_vect = Eigen::Vector3d::Zero();
                Eigen::Vector3d acc_vect = Eigen::Vector3d::Zero();
                for (int j = 0; j < 3; j++)
                {
                    pos_vect(j) = s1[j].pos[i];
                    vel_vect(j) = s1[j].vel[i];
                    acc_vect(j) = s1[j].acc[i];
                }
                s3.pos.push_back(pos_vect);
                s3.vel.push_back(vel_vect);
                s3.acc.push_back(acc_vect);
            }

            return s3;
        };

        
        inline bs_pva_state_1d get_uni_bspline_1d(
            int order, vector<double> timespan, 
            vector<double> ctrlpt, int knotdiv)
        {
            bs_utils b;
            bs_pva_state_1d s;
            b.k = order + 1;
            b.M = create_m(order);
            b.k_d = knotdiv;
            int n = (int)ctrlpt.size() - 1;
            int m = n + order + 1; 
            // int n_k = m + 1; // Number of knots

            // Range of index to evaluate accordingly (order to length of control points)
            b.r = int_range_to_vector(order, m - order - 1); 

            b.dt = (timespan[1] - timespan[0]) / (b.r.size());
            vector<double> t = tmf.linspace(timespan[0], timespan[1], b.r.size() + 1);

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
                vector<double> span = tmf.linspace((double)idx, (double)nxt_idx, b.k_d); 
                // Time in abs time (simulation time / given time)
                vector<double> actualspan = tmf.linspace(t[i], t[i+1], b.k_d); 

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