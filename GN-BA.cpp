#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "sophus/se3.h"


using namespace std;


typedef vector<Vector3d, Eigen::aligned_allocator<Vector3d>> VecVector3d;
typedef vector<Vector2d, Eigen::aligned_allocator<Vector3d>> VecVector2d;
typedef Matrix<double, 6, 1> Vector6d;

string p3d_file = "../p3d.txt";
string p2d_file = "../p2d.txt";


int main(int argc, char **argv) {

    VecVector2d p2d;
    VecVector3d p3d;
    Matrix3d K;
    double fx = 520.9, fy = 521.0, cx = 325.1, cy = 249.7;
    K << fx, 0, cx, 0, fy, cy, 0, 0, 1;

    ifstream fin1(p3d_file);
    while(fin1)
    {
        double data1[3] = {0};  //其中data1[0]表示时间，跳过（略去）
        for (auto &d1:data1)
            fin1 >> d1;
        Eigen::Vector3d temp1(data1[0], data1[1], data1[2]);    //获取平移数据
        //cout<<data1[0]<<data1[1]<<data1[2]<<endl;
        p3d.push_back(temp1);
    }

    ifstream fin2(p2d_file);
    while(fin2)
    {
        double data2[2] = {0};  //其中data1[0]表示时间，跳过（略去）
        for (auto &d2:data2)
            fin2 >> d2;
        Eigen::Vector2d temp2(data2[0], data2[1]);    //获取平移数据
        //cout<<data2[0]<<data2[1]<<endl;
        p2d.push_back(temp2);
    }
        p3d.pop_back();  //删除迭代器最后一个元素（最后一个值为0）
        p2d.pop_back();


    assert(p3d.size() == p2d.size());  //作用是如果它的条件返回错误，则终止程序执行

    int iterations = 100;
    double cost = 0, lastCost = 0;
    int nPoints = p3d.size();
    cout << "points: " << nPoints << endl;

    Sophus::SE3 T_esti; // estimated pose   默认初始化是所有值均为0

    for (int iter = 0; iter < iterations; iter++) {

        Matrix<double, 6, 6> H = Matrix<double, 6, 6>::Zero();
        Vector6d b = Vector6d::Zero();
        cost = 0;
        // compute cost
        for (int i = 0; i < nPoints; i++) {
            Vector4d p4d_world_i (p3d[i](0,0),p3d[i](1,0),p3d[i](2,0),1);
            Vector2d p2d_i = p2d[i];
            Vector4d bianhuan_4d=T_esti.matrix()*p4d_world_i;
            Vector3d p3d_carmera_i(bianhuan_4d(0,0),bianhuan_4d(1,0),bianhuan_4d(2,0));
            Vector3d u_v_1 = (1/bianhuan_4d(2,0))*K*p3d_carmera_i;
            Vector2d u_v(u_v_1(0,0),u_v_1(1,0));

            Vector2d e(0,0);  //（0,0）初始化,使二个维度为0，而不是取坐标号，注意误差e应该是个二维的量
            e = p2d_i - u_v;


/*****************观测SE中的.matrix(),.rotation_matrix()和.translation()的关系************************/
            /*Matrix4d temp_matrix4d=Matrix4d::Identity(4,4);  //临时单位阵
            cout<<"temp_matrix4d=\n"<<temp_matrix4d<<endl;  //输出是否单位化
            Matrix4d matrix_guancha_4d = T_esti.matrix()* temp_matrix4d;  //不能直接cout<<T_esti.matrix()
            cout<<"matrix_guancha_4d=\n"<<matrix_guancha_4d<<endl;  //输出变换阵


            Matrix3d temp_matrix3d=Matrix3d::Identity(3,3); //临时单位阵
            Matrix3d matrix_ratation_3d = T_esti.rotation_matrix()*temp_matrix3d;
            cout<<"matrix_ratation_3d=\n"<<matrix_ratation_3d<<endl;  //输出旋转阵

            Vector3d temp_vector3d(0,0,0); //临时零阵
            Vector3d vector_translation_3d=T_esti.translation()+temp_vector3d;
            cout<<"vector_translation_3d=\n"<<vector_translation_3d<<endl; //输出平移阵 */

/****************************************************************************************************/

	    // compute jacobian
            Matrix<double, 2, 6> J;
            J(0,0)=-(fx/p3d_carmera_i(2,0));
            J(0,1)=0;
            J(0,2)=fx*p3d_carmera_i(0,0)/(p3d_carmera_i(2,0)*p3d_carmera_i(2,0));
            J(0,3)=fx*p3d_carmera_i(0,0)*p3d_carmera_i(1,0)/(p3d_carmera_i(2,0)*p3d_carmera_i(2,0));
            J(0,4)=-(fx+fx*p3d_carmera_i(0,0)*p3d_carmera_i(0,0)/(p3d_carmera_i(2,0)*p3d_carmera_i(2,0)));
            J(0,5)=fx*p3d_carmera_i(1,0)/p3d_carmera_i(2,0);

            J(1,0)=0;
            J(1,1)=-fy/p3d_carmera_i(2,0);
            J(1,2)=fy*p3d_carmera_i(1,0)/(p3d_carmera_i(2,0)*p3d_carmera_i(2,0));
            J(1,3)=fy+fy*p3d_carmera_i(1,0)*p3d_carmera_i(1,0)/(p3d_carmera_i(2,0)*p3d_carmera_i(2,0));
            J(1,4)=-fy*p3d_carmera_i(0,0)*p3d_carmera_i(1,0)/(p3d_carmera_i(2,0)*p3d_carmera_i(2,0));
            J(1.5)=-fy*p3d_carmera_i(0,0)/p3d_carmera_i(2,0);
            
            H += J.transpose() * J;
            b += -J.transpose() * e;
            cost += e.norm();
        }

        Vector6d dx;

        dx=H.ldlt().solve(b);

        if (iter > 0 && cost > lastCost) {
            // cost increase, update is not good
            cout << "cost: " << cost << ", last cost: " << lastCost << endl;
            break;
        }

        T_esti = Sophus::SE3::exp(dx)*T_esti;  //这是SE3的乘法
        
        lastCost = cost;

        cout << "iteration " << iter << " cost=" << cout.precision(12) << cost << endl;
    }

    cout << "estimated pose: \n" << T_esti.matrix() << endl;
    return 0;
}
